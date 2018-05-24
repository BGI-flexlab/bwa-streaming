//
// Created by huangzhibo on 2017/4/16.
//

#include <strings.h>
#include <math.h>
#include "filter.h"
#include "ksort.h"
#include "utils.h"
#include "bwa.h"
extern int rg_number;

//// filter step (SOAPnuke)
static void worker(void *data, int i, int tid)
{
	filter_worker_t *w = (filter_worker_t*)data;
	if (!(w->filter_opt->is_pe)) {
        if (bwa_verbose >= 4) printf("=====> Processing read '%s' <=====\n", w->seqs[i].name);
//        w->is_clean = statistics_se(w->seqs[i<<1|0],w->seqs[i<<1|1], w->filter_opt, w->aux[tid]);
    } else{
        if (bwa_verbose >= 4) printf("=====> Processing read '%s'/1/2 <=====\n", w->seqs[i<<1|0].name);
        int is_clean = statistics_pe(&w->seqs[i<<1|0], &w->seqs[i<<1|1], w->filter_opt, w->read_info[tid]);
        w->seqs[i<<1|0].filter = !is_clean;
        w->seqs[i<<1|1].filter = !is_clean;
        if (bwa_verbose >= 4) printf("=====> Processing read filter stat: '%d' <=====\n", is_clean);
    }
}

void soapnuke_filter(const filter_opt_t *opt, int64_t n_processed, int n, bseq1_t *seqs, read_info_t *read_info[])
{
    extern void kt_for(int n_threads, void (*func)(void*,int,int), void *data, int n);
    filter_worker_t *w;
    double ctime, rtime;
    int i, j;

    ctime = cputime(); rtime = realtime();
    w = calloc(1, sizeof(read_info_t));
    w->filter_opt = opt;
    w->seqs = seqs;
    w->n_processed = n_processed;
//    w->read_info = read_info;  //fixme

    w->read_info = malloc(opt->n_threads * sizeof(read_info_t**));
    for (i = 0; i < opt->n_threads; ++i) {
        if(rg_number < 0){
            w->read_info[i] = (read_info_t**)malloc(sizeof(read_info_t *));
            w->read_info[i][0] = read_info_init(opt->is_pe);
        }else{
            w->read_info[i] = (read_info_t**)malloc(rg_number*sizeof(read_info_t *));
            for(j=0; j<rg_number; j++){
                w->read_info[i][j] = read_info_init(opt->is_pe);
            }
        }
    }

    kt_for(opt->n_threads, worker, w, (opt->is_pe)? n>>1 : n); // generate alignment

    for (i = 0; i < opt->n_threads; ++i) {
        if(rg_number < 0){
            merge_report(read_info[0], w->read_info[i][0]);
            read_info_destroy(w->read_info[i][0]);
            free(w->read_info[i]);
        }else{
            for(j=0; j<rg_number; j++){
                merge_report(read_info[j], w->read_info[i][j]);
                read_info_destroy(w->read_info[i][j]);
            }
            free(w->read_info[i]);
        }
    }
    free(w->read_info);

    if (bwa_verbose >= 3)
        fprintf(stderr, "[M::%s] Processed %d reads in %.3f CPU sec, %.3f real sec\n", __func__, n, cputime() - ctime, realtime() - rtime);
    free(w);
}

// todo:  resize w.seqs
void remove_bad_reads(ktp_data_t *data)
{
    int i, index;
    bseq1_t* filter_seqs=malloc(data->n_seqs*sizeof(bseq1_t));
    for(i=0, index=0; i < data->n_seqs; i++){
        if(!data->seqs[i].filter){
            filter_seqs[index] = data->seqs[i];
            index++;
        } else{
            free(data->seqs[i].name); free(data->seqs[i].comment);
            free(data->seqs[i].seq); free(data->seqs[i].qual);
        }
    }
    free(data->seqs);
    data->n_seqs = index;
    data->seqs = filter_seqs;
}

int statistics_pe(bseq1_t *read1, bseq1_t *read2, const filter_opt_t *opt, read_info_t *read_info[])
{
    StatisInfo si1, si2;
    read_info_t *report = NULL;
    memset(&si1, 0, sizeof(StatisInfo));
    memset(&si2, 0, sizeof(StatisInfo));

    if(rg_number>0){
        if(read1->comment){
            int sample_id = atoi(read1->comment);
            if(sample_id > rg_number){
                fprintf(stderr,"sample id is large than read group number\n");
                exit(1);
            }else{
                report = read_info[sample_id];
            }
        } else{
            fprintf(stderr,"comment info is empty\n");
            exit(1);
        }
    }else{
        report = read_info[0];
    }

    int trim_tail1 = opt->trim[1];
    int trim_tail2 = opt->trim[3];
    if(opt->adp1 != NULL){
        int index1 = adapter_align(read1, opt->adp1, opt);
        if(index1 == -2){
            si1.hasAdpt = 1;
        } else if(index1 >= 0) {
            int cut_adapter_len = read1->l_seq - index1;
            trim_tail1 = cut_adapter_len > trim_tail1 ? cut_adapter_len : trim_tail1;
            report->total_cut_adapter_num++;
        }
    }

    if(opt->adp2 != NULL) {
        int index2 = adapter_align(read2, opt->adp2, opt);

        if (index2 == -2) {
            si2.hasAdpt = 1;
        } else if (index2 >= 0) {
            int cut_adapter_len = read2->l_seq - index2;
            trim_tail2 = cut_adapter_len > trim_tail2 ? cut_adapter_len : trim_tail2;
            report->total_cut_adapter_num++;
        }
    }

    report->total_raw_read_num++;

    FqInfo *info1 = report->read1_info;
    FqInfo *info2 = report->read2_info;

    //fq1
    seq_stat(read1, opt, opt->trim[0], trim_tail1, info1, &si1);
    //fq2
    seq_stat(read2, opt, opt->trim[2], trim_tail2, info2, &si2);

    // filter adapter read
    if(si1.hasAdpt || si2.hasAdpt){
        if (si1.hasAdpt){
            info1->adapter_num++;
        }
        if (si2.hasAdpt){
            info2->adapter_num++;
        }
        report->total_adapter_num++;
        return 0;
    }

    // filter short length read
    if (read1->l_seq < opt->min_read_len || read2->l_seq < opt->min_read_len){
        if(read1->l_seq < opt->min_read_len){
            info1->short_length_n++;
        }
        if(read2->l_seq < opt->min_read_len){
            info2->short_length_n++;
        }
        report->total_short_length_num++;
        return 0;
    }

    if(si1.nExceed || si2.nExceed){
        if (si1.nExceed){
            info1->n_exceed_num++;
        }
        if (si2.nExceed) {
            info2->n_exceed_num++;
        }
        report->total_n_exceed_num++;
        return 0;
    }

    if(si1.isLowQual || si2.isLowQual){
        if (si1.isLowQual){
            info1->low_qual_num++;
        }
        if (si2.isLowQual) {
            info2->low_qual_num++;
        }
        report->total_low_qual_num++;
        return 0;
    }

    if(si1.sumQuality < opt->mean * read1->l_seq || si2.sumQuality < opt->mean * read2->l_seq){
        if(si1.sumQuality < opt->mean * read1->l_seq){
            info1->low_mean_num++;
        }
        if(si2.sumQuality < opt->mean * read2->l_seq){
            info2->low_mean_num++;
        }
        report->total_low_mean_num++;
        return 0;
    }

//    // todo polyA
//    if ((polyAType_==0 && (!si1.isPolyA || !si2.isPolyA)) ||(polyAType_==1 && !si1.isPolyA && !si2.isPolyA)) //not polyA
//    {
//    }
//    else //has polyA
//    {
//        if(polyAType_==0){
//            info->polyANum++;
//            info2->polyANum++;
//        }else{
//            if(si1.isPolyA)
//                info->polyANum++;
//            if(si2.isPolyA)
//                info2->polyANum++;
//        }
//        info->totalPolyANum++;
//    }

    info1->clean_base_A += si1.a;
    info1->clean_base_C += si1.c;
    info1->clean_base_G += si1.g;
    info1->clean_base_T += si1.t;
    info1->clean_base_N += si1.n;
    info1->clean_q20 += si1.q20;
    info1->clean_q30 += si1.q30;
    info1->clean_base_num += read1->l_seq;
    calculate_base_distribution(read1, info1);

    info2->clean_base_A += si2.a;
    info2->clean_base_C += si2.c;
    info2->clean_base_G += si2.g;
    info2->clean_base_T += si2.t;
    info2->clean_base_N += si2.n;
    info2->clean_q20 += si2.q20;
    info2->clean_q30 += si2.q30;
    info2->clean_base_num += read2->l_seq;
    calculate_base_distribution(read2, info2);

    report->total_clean_read_num++;

    return 1;
}

void calculate_base_distribution(bseq1_t *read, FqInfo *info)
{
    int qual, i;

    for (i=0; i<read->l_seq; ++i) {
        switch (read->seq[i])
        {
            case 'A':
                ++info->clean_base[i][0];
                break;
            case 'C':
                ++info->clean_base[i][1];
                break;
            case 'G':
                ++info->clean_base[i][2];
                break;
            case 'T':
                ++info->clean_base[i][3];
                break;
            case 'N':
                ++info->clean_base[i][4];
                break;
            default:break;
        }

        qual = read->qual[i] - 33;
        ++info->clean_qual[i][qual];

        if (qual >= 20)
        {
            ++info->clean_q20q30[i][0];
            if (qual >= 30)
            {
                ++info->clean_q20q30[i][1];
            }
        }
    }
}

void seq_stat(bseq1_t *read, const filter_opt_t *opt, int head_trim_n, int tail_trim_n, FqInfo *info, StatisInfo *si)
{
    int qual, i;

    info->max_raw_read_len = read->l_seq > info->max_raw_read_len ? read->l_seq : info->max_raw_read_len;
    info->raw_base_num += read->l_seq;

    //todo 切除低质量末端

    int a = 0, g = 0, c = 0, t = 0, n = 0;
    int q20 = 0, q30 = 0;

    int right = read->l_seq - tail_trim_n;
    int sumQual = 0;

    for (i=0; i<read->l_seq; ++i)
    {
        switch (read->seq[i])
        {
            case 'A':
                a++;
                info->base[i][0]++;
                if (i>=head_trim_n && i<right)
                {
                    si->a++;
                }
                break;
            case 'C':
                c++;
                info->base[i][1]++;
                if (i>=head_trim_n && i<right)
                {
                    si->c++;
                }
                break;
            case 'G':
                g++;
                info->base[i][2]++;
                if (i>=head_trim_n && i<right)
                {
                    si->g++;
                }
                break;
            case 'T':
                t++;
                info->base[i][3]++;
                if (i>=head_trim_n && i<right)
                {
                    si->t++;
                }
                break;
            case 'N':
                n++;
                info->base[i][4]++;
                if (i>=head_trim_n && i<right)
                {
                    si->n++;
                }
                break;
            default:break;
        }

        qual = read->qual[i] - 33;

        if (qual > MAX_QUALITY)
        {
            fprintf(stderr, "[W::%s] some bases' quality larger than %d, they have been set to %d.\n", __func__, MAX_QUALITY, MAX_QUALITY);
            qual = MAX_QUALITY;
        }

        if (qual < 0)
        {
            fprintf(stderr, "[W::%s] some bases' quality smaller than 0, they have been set to 0.\n", __func__);
            qual = 0;
        }

        if (qual > info->max_quality_value)
        {
            info->max_quality_value = qual;
        }


        if (i>=head_trim_n && i<right)
        {
            sumQual += qual;
        }

        info->qual[i][qual]++;

        if (qual >= 20)
        {
            q20++;
            info->q20q30[i][0]++;
            if (qual >= 30)
            {
                q30++;
                info->q20q30[i][1]++;
            }
        }

        if (i>=head_trim_n && i<right)
        {
            if (qual <= opt->lowQual)
            {
                si->lowQual++;
            }
            else if (qual >=20)
            {
                si->q20++;
                if (qual >= 30)
                {
                    si->q30++;
                }
            }
        }
    }

    info->raw_base_A += a;
    info->raw_base_C += c;
    info->raw_base_G += g;
    info->raw_base_T += t;
    info->raw_base_N += n;
    info->raw_q20 += q20;
    info->raw_q30 += q30;

    //clean data read length
    read->l_seq = read->l_seq - head_trim_n - tail_trim_n;
    info->max_clean_read_len = read->l_seq > info->max_clean_read_len ? read->l_seq : info->max_clean_read_len;

    if(read->l_seq < 0){
        read->l_seq = 0;
        read->seq[0] = '\0';
        read->seq[0] = '\0';
    }else{
        //截断read的两端
        char* seq_tmp = strdup(read->seq+head_trim_n);
        seq_tmp[read->l_seq] = '\0';
        free(read->seq);
        read->seq = seq_tmp;

        char* qual_tmp = strdup(read->qual+head_trim_n);
        seq_tmp[read->l_seq] = '\0';
        free(read->qual);
        read->qual = qual_tmp;

//        strncpy(read->seq,read->seq+head_trim_n , (size_t) read->l_seq);
//        read->seq[read->l_seq] = '\0';

//        strncpy(read->qual,read->qual+head_trim_n , (size_t) read->l_seq);
//        read->qual[read->l_seq] = '\0';
    }

    si->nExceed = (si->n >= read->l_seq * opt->nRate);
    si->isLowQual = (si->lowQual >= read->l_seq * opt->qualRate);
    si->sumQuality = sumQual;

    if (opt->polyA >= 1E-6)
    {
        si->isPolyA = (1.0 * si->a / read->l_seq) >= (opt->polyA - 1E-6);
    }
}

//int has_adapter();  todo adapterList 情况
//int find_adappter(bseq1_t *read, char* adapter, int adpLen,StatisInfo si)
//{
//    int adp_index = adapter_align(read->seq, readLen, adapter.c_str(), adpLen);
//    if(adp_index != -1)
//    {
//        if(!cutAdaptor || (cutAdaptor && index < minReadLength))
//        {
//            si.hasAdpt = 1;
//        }
//    }
//    return adp_index;
//}

// -1: no adapter  -2: filter adapter due to the adapter is too long  >0: adapter index to trim
int adapter_align(bseq1_t *read, const char *adapter, const filter_opt_t *opt)
{
    int find = -1, c;
    int adptLen = (int) strlen(adapter);
    int minMatchLen = (int) ceilf(adptLen * opt->matchRatio);
    int a1 = adptLen - minMatchLen;
    int r1 = 0;
    int len, mis;
//    return -1; //fixme

    int right = read->l_seq - minMatchLen;

    for (r1 = 0; r1 <= right;)
    {
        int len1 = adptLen - a1;
        int len2 = read->l_seq - r1;
        len = (len1 < len2) ? len1 : len2;
        mis = 0;
        int map[MAX_LENGTH];
        map[0] = 0;
        for (c = 0; c < len; ++c)
        {
            if (adapter[a1 + c] == read->seq[r1 + c])
            {
                map[mis]++;
            }
            else
            {
                mis++;
                map[mis] = 0;
            }
        }
        int max_map = 0;
        for (c = 0; c <= mis; ++c)
        {
            if (map[c] > max_map)
            {
                max_map = map[c];
            }
        }
        if ((mis <= opt->misMatch) || (max_map >= minMatchLen))
        {
            find = r1;
            break;
        }
        if (a1 > 0)
        {
            a1--;
        }
        else
        {
            r1++;
        }
    }

    if(find != -1)
    {
        if(!opt->cutAdapter || (opt->cutAdapter && find < opt->min_read_len))
        {
            find = -2;
        }
    }

    return find;
}

filter_opt_t *filter_opt_init()
{
    int i;
    filter_opt_t *o;
    o = calloc(1, sizeof(filter_opt_t));
    o->skip_filter = 1;
    o->hold_reads = 0;
    o->is_phred64 = 0;
    o->is_pe = 1;
    o->n_threads = 1;
    o->adp1 = NULL;
    o->adp2 = NULL;
    o->tile = NULL;
    o->misMatch = 1;
    o->cutAdapter = 0;
    o->matchRatio = 0.5;
    o->lowQual = 5;
    o->qualRate = 0.5;
    o->nRate = 0.05;
    o->mean = 0;
    o->small = -1;
    o->min_read_len = 20;
    o->polyA = 0;
    o->polyAType = 0;
    for(i=0; i < 4; i++){
        o->trim[i] = 0;
    }
    return o;
}

void filter_opt_destroy(filter_opt_t *o)
{
    if(NULL != o->adp1)
            free(o->adp1);
    if(NULL != o->adp2)
        free(o->adp2);
    if(NULL != o->tile)
        free(o->tile);
    free(o);
}

read_info_t *read_info_init(int is_pe)
{
    read_info_t *o = calloc(1, sizeof(read_info_t));
    o->read1_info = calloc(1, sizeof(FqInfo));
    o->read1_info->max_quality_value = 41;
    if(is_pe) {
        o->read2_info = calloc(1, sizeof(FqInfo));
        o->read2_info->max_quality_value = 41;
    }
    else
        o->read2_info = NULL;

    return o;
}

void merge_fq_info(FqInfo *fqinfo, FqInfo *ptr) {
    int i, j;

    fqinfo->max_quality_value = fqinfo->max_quality_value > ptr->max_quality_value ? fqinfo->max_quality_value : ptr->max_quality_value;
    fqinfo->max_raw_read_len = fqinfo->max_raw_read_len > ptr->max_raw_read_len ? fqinfo->max_raw_read_len : ptr->max_raw_read_len;
    fqinfo->max_clean_read_len = fqinfo->max_clean_read_len > ptr->max_clean_read_len ? fqinfo->max_clean_read_len : ptr->max_clean_read_len;

    fqinfo->raw_base_num += ptr->raw_base_num;
    fqinfo->clean_base_num += ptr->clean_base_num;
    fqinfo->raw_base_A += ptr->raw_base_A;
    fqinfo->clean_base_A += ptr->clean_base_A;
    fqinfo->raw_base_C += ptr->raw_base_C;
    fqinfo->clean_base_C += ptr->clean_base_C;
    fqinfo->raw_base_G += ptr->raw_base_G;
    fqinfo->clean_base_G += ptr->clean_base_G;
    fqinfo->raw_base_T += ptr->raw_base_T;
    fqinfo->clean_base_T += ptr->clean_base_T;
    fqinfo->raw_base_N += ptr->raw_base_N;
    fqinfo->clean_base_N += ptr->clean_base_N;
    fqinfo->raw_q20 += ptr->raw_q20;
    fqinfo->clean_q20 += ptr->clean_q20;
    fqinfo->raw_q30 += ptr->raw_q30;
    fqinfo->clean_q30 += ptr->clean_q30;

    fqinfo->adapter_num += ptr->adapter_num;
    fqinfo->n_exceed_num += ptr->n_exceed_num;
    fqinfo->low_qual_num += ptr->low_qual_num;
    fqinfo->low_mean_num += ptr->low_mean_num;
    fqinfo->small_insert_num += ptr->small_insert_num;
    fqinfo->polyA_num += ptr->polyA_num;
    fqinfo->short_length_n += ptr->short_length_n;

    for(i=0; i < MAX_LENGTH; i++){
        fqinfo->clean_read_len_distribution[i] += fqinfo->clean_read_len_distribution[i];
        for(j=0; j < 5; j++){
            fqinfo->base[i][j] += ptr->base[i][j];
            fqinfo->clean_base[i][j] += ptr->clean_base[i][j];
        }
        for(j=0; j < 2; j++){
            fqinfo->q20q30[i][j] += ptr->q20q30[i][j];
            fqinfo->clean_q20q30[i][j] += ptr->clean_q20q30[i][j];
        }
        for(j=0; j <= MAX_QUALITY; j++){
            fqinfo->qual[i][j] += ptr->qual[i][j];
            fqinfo->clean_qual[i][j] += ptr->clean_qual[i][j];
        }
    }
}

void merge_report(read_info_t *info, read_info_t *ptr) {
    info->total_raw_read_num += ptr->total_raw_read_num;
    info->total_clean_read_num += ptr->total_clean_read_num;
    info->total_short_length_num += ptr->total_short_length_num;
    info->total_n_exceed_num += ptr->total_n_exceed_num;
    info->total_low_qual_num += ptr->total_low_qual_num;
    info->total_low_mean_num += ptr->total_low_mean_num;
    info->total_adapter_num += ptr->total_adapter_num;
    info->total_cut_adapter_num += ptr->total_cut_adapter_num;
    info->total_small_insert_num += ptr->total_small_insert_num;
    info->total_polyA_num += ptr->total_polyA_num;

    merge_fq_info(info->read1_info, ptr->read1_info);
    if(NULL != info->read2_info)
        merge_fq_info(info->read2_info, ptr->read2_info);

}


void read_info_destroy(read_info_t *read_info) {
    free(read_info->read1_info);
    if(NULL != read_info->read2_info)
        free(read_info->read2_info);
    free(read_info);
}

void print_fq_info(FqInfo *fqInfo)
{
    err_printf("%d\t%d\t",fqInfo->max_raw_read_len, fqInfo->max_clean_read_len);
    err_printf("%lu\t%lu\t%lu\t%lu\t%lu\t",fqInfo->raw_base_num, fqInfo->clean_base_num, fqInfo->raw_base_A, fqInfo->clean_base_A, fqInfo->raw_base_C);
    err_printf("%lu\t%lu\t%lu\t%lu\t%lu\t",fqInfo->clean_base_C, fqInfo->raw_base_G, fqInfo->clean_base_G, fqInfo->raw_base_T, fqInfo->clean_base_T);
    err_printf("%lu\t%lu\t%lu\t%lu\t%lu\t",fqInfo->raw_base_N, fqInfo->clean_base_N, fqInfo->raw_q20, fqInfo->clean_q20, fqInfo->raw_q30);
    err_printf("%lu\t%lu\t%lu\t%lu\t",fqInfo->clean_q30, fqInfo->adapter_num, fqInfo->n_exceed_num, fqInfo->low_qual_num);
    err_printf("%lu\t%lu\t%lu\t",fqInfo->low_mean_num, fqInfo->small_insert_num, fqInfo->polyA_num);
    err_printf("%d#S\n",fqInfo->max_quality_value);

    //base distributions by read position
    err_printf("%s\t#S\n", "#Base_distributions_by_read_position");
    int i, j;
    for (i=0; i<fqInfo->max_raw_read_len; ++i)
    {
        err_printf("%lu", fqInfo->base[i][0]);
        for (j=1; j<5; j++)
        {
            err_printf("\t%lu", fqInfo->base[i][j]);
        }

        err_printf("\t%lu", fqInfo->clean_base[i][0]);
        for (j=1; j<5; j++)
        {
            err_printf("\t%lu", fqInfo->clean_base[i][j]);
        }
        err_fputs("#S\n", stdout);
    }

    //Raw Base_quality_value_distribution_by_read_position && Distribution_of_Q20_Q30_bases_by_read_position
    err_printf("%s\t#S\n", "#Raw_Base_quality_value_distribution_by_read_position");
    for (i=0; i<fqInfo->max_raw_read_len; ++i)
    {
        err_printf("%lu\t%lu", fqInfo->q20q30[i][0],fqInfo->q20q30[i][1]);
        for (j=0; j<=fqInfo->max_quality_value; j++)
        {
            err_printf("\t%lu", fqInfo->qual[i][j]);
        }
        err_fputs("#S\n", stdout);
    }

    //Clean Base_quality_value_distribution_by_read_position && Distribution_of_Q20_Q30_bases_by_read_position
    err_printf("%s\t#S\n", "#Clean_Base_quality_value_distribution_by_read_position");
    for (i=0; i<fqInfo->max_clean_read_len; ++i)
    {
        err_printf("%lu\t%lu", fqInfo->clean_q20q30[i][0],fqInfo->clean_q20q30[i][1]);
        for (j=0; j<=fqInfo->max_quality_value; j++)
        {
            err_printf("\t%lu", fqInfo->clean_qual[i][j]);
        }
        err_fputs("#S\n", stdout);
    }
}

void report_print(read_info_t *read_info)
{
    err_printf(">%d\t%s\t#S\n", read_info->id, read_info->rg_id);
    err_fputs("#Total_statistical_information\t#S\n", stdout);
    err_printf("%lu\t%lu\t",read_info->total_raw_read_num, read_info->total_clean_read_num);
    err_printf("%lu\t%lu\t%lu\t",read_info->total_short_length_num, read_info->total_n_exceed_num, read_info->total_low_qual_num);
    err_printf("%lu\t%lu\t%lu\t",read_info->total_low_mean_num, read_info->total_adapter_num, read_info->total_cut_adapter_num);
    err_printf("%lu\t%lu#S\n",read_info->total_small_insert_num, read_info->total_polyA_num);
    err_fputs("#Fq1_statistical_information\t#S\n", stdout);
    print_fq_info(read_info->read1_info);
    if (NULL != read_info->read2_info){
        err_fputs("#Fq2_statistical_information\t#S\n", stdout);
        print_fq_info(read_info->read2_info);
    }
}
