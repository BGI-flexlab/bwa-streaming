//
// Created by huangzhibo on 2017/4/16.
//

#include <strings.h>
#include <math.h>
#include "filter.h"
#include "ksort.h"
#include "utils.h"
#include "bwa.h"

//// filter step (SOAPnuke)
static void worker(void *data, int i, int tid)
{
	filter_worker_t *w = (filter_worker_t*)data;
	if (!(w->filter_opt->is_pe)) {
        if (bwa_verbose >= 4) printf("=====> Processing read '%s' <=====\n", w->seqs[i].name);
//        w->is_clean = statistics_se(w->seqs[i<<1|0],w->seqs[i<<1|1], w->filter_opt, w->aux[tid]);
    } else{
        if (bwa_verbose >= 4) printf("=====> Processing read '%s'/1/2 <=====\n", w->seqs[i<<1|0].name);
        int is_clean = statistics_pe(&w->seqs[i<<1|0], &w->seqs[i<<1|1], w->filter_opt, w->fq_info);
        w->seqs[i<<1|0].filter = !is_clean;
        w->seqs[i<<1|1].filter = !is_clean;
        if (bwa_verbose >= 4) printf("=====> Processing read filter stat: '%d' <=====\n", is_clean);
    }
}

void soapnuke_filter(const filter_opt_t *opt, int64_t n_processed, int n, bseq1_t *seqs, FqInfo *fq_info)
{
    extern void kt_for(int n_threads, void (*func)(void*,int,int), void *data, int n);
    filter_worker_t w;
    double ctime, rtime;

    ctime = cputime(); rtime = realtime();
    w.filter_opt = opt;
    w.seqs = seqs;
    w.n_processed = n_processed;
    w.fq_info = fq_info;  //fixme

    kt_for(opt->n_threads, worker, &w, (opt->is_pe)? n>>1 : n); // generate alignment
    if (bwa_verbose >= 3)
        fprintf(stderr, "[M::%s] Processed %d reads in %.3f CPU sec, %.3f real sec\n", __func__, n, cputime() - ctime, realtime() - rtime);
}

// todo:  resize w.seqs
void remove_bad_reads(ktp_data_t *data)
{
    int i, index;
    for(i=0, index=-1; i < data->n_seqs;i++){
        if(!data->seqs[i].filter){
            index ++;
            if(i != index){
                memcpy(data->seqs + index, data->seqs + i, sizeof(bseq1_t));
            }
        }
    }
    data->n_seqs = index + 1;
    realloc(data->seqs, data->n_seqs* sizeof(bseq1_t));
}

int statistics_pe(bseq1_t *read1, bseq1_t *read2, const filter_opt_t *opt, FqInfo *info)
{
    FqInfo *info2 = info + 1;
    StatisInfo si1, si2;
    memset(&si1, 0, sizeof(StatisInfo));
    memset(&si2, 0, sizeof(StatisInfo));

    int trim_tail1 = opt->trim[1];
    int trim_tail2 = opt->trim[3];
    if(opt->adp1 != NULL){
        int index1 = adapter_align(read1, opt->adp1, opt);
        if(index1 == -2){
            si1.hasAdpt = 1;
        } else if(index1 >= 0) {
            int cut_adapter_len = read1->l_seq - index1;
            trim_tail1 = cut_adapter_len > trim_tail1 ? cut_adapter_len : trim_tail1;
            info->totalCutAdaptorNum++;
        }
    }

    if(opt->adp2 != NULL) {
        int index2 = adapter_align(read2, opt->adp2, opt);

        if (index2 == -2) {
            si2.hasAdpt = 1;
        } else if (index2 >= 0) {
            int cut_adapter_len = read2->l_seq - index2;
            trim_tail2 = cut_adapter_len > trim_tail2 ? cut_adapter_len : trim_tail2;
            info2->totalCutAdaptorNum++;
        }
    }

    //fq1
    seq_stat(read1, opt, opt->trim[0], trim_tail1, info, &si1);
    //fq2
    seq_stat(read2, opt, opt->trim[2], trim_tail2, info2, &si2);

    // filter adapter read
    if(si1.hasAdpt || si2.hasAdpt){
        if (si1.hasAdpt){
            info->adapterNum++;
        }
        if (si2.hasAdpt){
            info2->adapterNum++;
        }
        info->totalAdapterNum++;
        return 0;
    }

    // filter short length read
    if (read1->l_seq < opt->min_read_len || read2->l_seq < opt->min_read_len){
        if(read1->l_seq < opt->min_read_len){
            info->short_length_n++;
        }
        if(read2->l_seq < opt->min_read_len){
            info2->short_length_n++;
        }
        info->total_short_length_n++;
        return 0;
    }

    if(si1.nExceed || si2.nExceed){
        if (si1.nExceed){
            info->nExceedNum++;
        }
        if (si2.nExceed) {
            info2->nExceedNum++;
        }
        info->totalNExceedNum++;
        return 0;
    }

    if(si1.isLowQual || si2.isLowQual){
        if (si1.isLowQual){
            info->lowQualNum++;
        }
        if (si2.isLowQual) {
            info2->lowQualNum++;
        }
        info->totalLowQualNum++;
        return 0;
    }

    if(si1.sumQuality < opt->mean * read1->l_seq || si2.sumQuality < opt->mean * read2->l_seq){
        if(si1.sumQuality < opt->mean * read1->l_seq){
            info->lowMeanNum++;
        }
        if(si2.sumQuality < opt->mean * read2->l_seq){
            info2->lowMeanNum++;
        }
        info->totalLowMeanNum++;
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

    info->cleanBaseA += si1.a;
    info->cleanBaseC += si1.c;
    info->cleanBaseG += si1.g;
    info->cleanBaseT += si1.t;
    info->cleanBaseN += si1.n;
    info->cleanQ20 += si1.q20;
    info->cleanQ30 += si1.q30;
    info->cleanTotalReadNum++;
    info->cleanTotalBaseNum += read1->l_seq;
    calculate_base_distribution(read1, info);

    info2->cleanBaseA += si2.a;
    info2->cleanBaseC += si2.c;
    info2->cleanBaseG += si2.g;
    info2->cleanBaseT += si2.t;
    info2->cleanBaseN += si2.n;
    info2->cleanQ20 += si2.q20;
    info2->cleanQ30 += si2.q30;
    info2->cleanTotalReadNum++;
    info2->cleanTotalBaseNum += read2->l_seq;
    calculate_base_distribution(read2, info2);

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

        qual = read->qual[i];
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

    info->rawReadLength = read->l_seq > info->rawReadLength ? read->l_seq : info->rawReadLength;
    info->rawTotalReadNum++;
    info->rawTotalBaseNum += read->l_seq;

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

        if (qual > info->maxQualityValue)
        {
            info->maxQualityValue = qual;
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

    info->rawBaseA += a;
    info->rawBaseC += c;
    info->rawBaseG += g;
    info->rawBaseT += t;
    info->rawBaseN += n;
    info->rawQ20 += q20;
    info->rawQ30 += q30;

    //clean data read length
    read->l_seq = read->l_seq - head_trim_n - tail_trim_n;
    info->cleanReadLength = read->l_seq > info->cleanReadLength ? read->l_seq : info->cleanReadLength;

    if(read->l_seq < 0){
        read->l_seq = 0;
        read->seq[0] = '\0';
        read->seq[0] = '\0';
    }else{
        //截断read的两端
        read->seq[right] = '\0';
        read->seq = read->seq + head_trim_n;
        read->qual[right] = '\0';
        read->qual = read->qual + head_trim_n;
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
    o->is_phred64 = 0;
    o->is_pe = 1;
    o->n_threads = 1;
    o->adp1 = NULL;
    o->adp2 = NULL;
    o->tile = NULL;
    o->misMatch = 1;
    o->cutAdapter = 1;
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

FqInfo *fq_info_init()
{
    int i, j;
    FqInfo *o = calloc(2, sizeof(FqInfo));
    FqInfo *o2 = o + 1;

    o->rawReadLength = 0;
    o->rawReadLength = 0;     //raw data 读长
    o->cleanReadLength = 0;   //clean data 读长
    o->rawTotalReadNum = 0;  //raw data read个数
    o->cleanTotalReadNum = 0;//clean data read 个数
    o->rawTotalBaseNum = 0;  //raw data total base number
    o->cleanTotalBaseNum = 0; //clean data total base number
    o->rawBaseA = 0;     //raw data base A number
    o->cleanBaseA = 0;   //clean data base A number
    o->rawBaseC = 0;     //raw data base C number
    o->cleanBaseC = 0; //clean data base C number
    o->rawBaseG = 0; //raw data base G number
    o->cleanBaseG = 0; //clean data base G number
    o->rawBaseT = 0; //raw data base T number
    o->cleanBaseT = 0; //clean data base T number
    o->rawBaseN = 0; //raw data base N number
    o->cleanBaseN = 0; //clean data base N number
    o->rawQ20 = 0; //rawfq文件中碱基质量>=20的碱基总数
    o->cleanQ20 = 0; //cleanfq文件中碱基质量>=20的碱基总数
    o->rawQ30 = 0; //rawfq文件中碱基质量>=30的碱基总数
    o->cleanQ30 = 0; //cleanfq文件中碱基质量>=30的碱基总数

    o->adapterNum = 0;  //the number of read which contain adapter in raw data
    o->nExceedNum = 0;  //the number of read which n rate was exceed in raw data
    o->lowQualNum = 0;  //low qualtiy read number in raw data
    o->lowMeanNum = 0;  //low mean quality read number in raw data
    o->smallInsertNum = 0;  //samll inert number in raw data
    o->polyANum = 0;    //polyA number in raw data
    o->short_length_n = 0;    //the read is too short

    o->total_short_length_n = 0;
    o->totalAdapterNum = 0;
    o->totalNExceedNum = 0;
    o->totalLowQualNum = 0;
    o->totalLowMeanNum = 0;
    o->totalSmallInsertNum = 0;
    o->totalPolyANum = 0;

    o->totalCutAdaptorNum = 0;

    o->maxQualityValue = 41;

    for(i=0; i < MAX_LENGTH; i++){
        o->cleanReadLengthDistribution[i] = 0;
        for(j=0; j < 5; j++){
            o->base[i][j] = 0;
            o->clean_base[i][j] = 0;
        }
        for(j=0; j < 2; j++){
            o->q20q30[i][j] = 0;
            o->clean_q20q30[i][j] = 0;
        }
        for(j=0; j <= MAX_QUALITY; j++){
            o->qual[i][j] = 0;
            o->clean_qual[i][j] = 0;
        }
    }

    o2->rawReadLength = 0;
    o2->rawReadLength = 0;     //raw data 读长
    o2->cleanReadLength = 0;   //clean data 读长
    o2->rawTotalReadNum = 0;  //raw data read个数
    o2->cleanTotalReadNum = 0;//clean data read 个数
    o2->rawTotalBaseNum = 0;  //raw data total base number
    o2->cleanTotalBaseNum = 0; //clean data total base number
    o2->rawBaseA = 0;     //raw data base A number
    o2->cleanBaseA = 0;   //clean data base A number
    o2->rawBaseC = 0;     //raw data base C number
    o2->cleanBaseC = 0; //clean data base C number
    o2->rawBaseG = 0; //raw data base G number
    o2->cleanBaseG = 0; //clean data base G number
    o2->rawBaseT = 0; //raw data base T number
    o2->cleanBaseT = 0; //clean data base T number
    o2->rawBaseN = 0; //raw data base N number
    o2->cleanBaseN = 0; //clean data base N number
    o2->rawQ20 = 0; //rawfq文件中碱基质量>=20的碱基总数
    o2->cleanQ20 = 0; //cleanfq文件中碱基质量>=20的碱基总数
    o2->rawQ30 = 0; //rawfq文件中碱基质量>=30的碱基总数
    o2->cleanQ30 = 0; //cleanfq文件中碱基质量>=30的碱基总数

    o2->adapterNum = 0;  //the number of read which contain adapter in raw data
    o2->nExceedNum = 0;  //the number of read which n rate was exceed in raw data
    o2->lowQualNum = 0;  //low qualtiy read number in raw data
    o2->lowMeanNum = 0;  //low mean quality read number in raw data
    o2->smallInsertNum = 0;  //samll inert number in raw data
    o2->polyANum = 0;    //polyA number in raw data
    o2->short_length_n = 0;    //the read is too short

    o2->total_short_length_n = 0;
    o2->totalAdapterNum = 0;
    o2->totalNExceedNum = 0;
    o2->totalLowQualNum = 0;
    o2->totalLowMeanNum = 0;
    o2->totalSmallInsertNum = 0;
    o2->totalPolyANum = 0;

    o2->totalCutAdaptorNum = 0;

    o2->maxQualityValue = 41;

    for(i=0; i < MAX_LENGTH; i++){
        o2->cleanReadLengthDistribution[i] = 0;
        for(j=0; j < 5; j++){
            o2->base[i][j] = 0;
            o2->clean_base[i][j] = 0;
        }
        for(j=0; j < 2; j++){
            o2->q20q30[i][j] = 0;
            o2->clean_q20q30[i][j] = 0;
        }
        for(j=0; j <= MAX_QUALITY; j++){
            o2->qual[i][j] = 0;
            o2->clean_qual[i][j] = 0;
        }
    }
    return o;
}

void print_fq_info(const FqInfo *fqInfo)
{
    err_printf("%d\t%d\t%lu\t%lu\t%d\t",fqInfo->rawReadLength, fqInfo->cleanReadLength, fqInfo->rawTotalReadNum, fqInfo->cleanTotalReadNum, 0);
    err_printf("%lu\t%lu\t%lu\t%lu\t%lu\t",fqInfo->rawTotalBaseNum, fqInfo->cleanTotalBaseNum, fqInfo->rawBaseA, fqInfo->cleanBaseA, fqInfo->rawBaseC);
    err_printf("%lu\t%lu\t%lu\t%lu\t%lu\t",fqInfo->cleanBaseC, fqInfo->rawBaseG, fqInfo->cleanBaseG, fqInfo->rawBaseT, fqInfo->cleanBaseT);
    err_printf("%lu\t%lu\t%lu\t%lu\t%lu\t",fqInfo->rawBaseN, fqInfo->cleanBaseN, fqInfo->rawQ20, fqInfo->cleanQ20, fqInfo->rawQ30);
    err_printf("%lu\t%d\t%lu\t%lu\t%lu\t",fqInfo->cleanQ30, 0, fqInfo->adapterNum, fqInfo->nExceedNum, fqInfo->lowQualNum);
    err_printf("%lu\t%lu\t%lu\t%d\t%lu\t",fqInfo->lowMeanNum, fqInfo->smallInsertNum, fqInfo->polyANum, 0, fqInfo->totalAdapterNum);
    err_printf("%lu\t%lu\t%lu\t%lu\t%lu\t",fqInfo->totalNExceedNum, fqInfo->totalLowQualNum, fqInfo->totalSmallInsertNum, fqInfo->totalPolyANum, fqInfo->totalCutAdaptorNum);
    err_printf("%d#S\n",fqInfo->maxQualityValue);

    //base distributions by read position
    err_printf("%s\t#S\n", "#Base_distributions_by_read_position");
    int i, j;
    for (i=0; i<fqInfo->rawReadLength; ++i)
    {
        err_printf("%lu", fqInfo->base[i][0]);
        for (j=1; j<5; j++)
        {
            err_printf("\t%lu", fqInfo->base[i][j]);
        }
        err_fputs("#S\n", stdout);

        err_printf("%lu", fqInfo->clean_base[i][0]);
        for (j=1; j<5; j++)
        {
            err_printf("\t%lu", fqInfo->clean_base[i][j]);
        }
        err_fputs("#S\n", stdout);
    }

    //Raw Base_quality_value_distribution_by_read_position && Distribution_of_Q20_Q30_bases_by_read_position
    err_printf("%s\t#S\n", "#Raw_Base_quality_value_distribution_by_read_position");
    for (i=0; i<fqInfo->rawReadLength; ++i)
    {
        err_printf("%lu\t%lu", fqInfo->q20q30[i][0],fqInfo->q20q30[i][1]);
        for (j=0; j<=fqInfo->maxQualityValue; j++)
        {
            err_printf("\t%lu", fqInfo->qual[i][j]);
        }
        err_fputs("#S\n", stdout);
    }

    //Clean Base_quality_value_distribution_by_read_position && Distribution_of_Q20_Q30_bases_by_read_position
    err_printf("%s\t#S\n", "#Clean_Base_quality_value_distribution_by_read_position");
    for (i=0; i<fqInfo->cleanReadLength; ++i)
    {
        err_printf("%lu\t%lu", fqInfo->clean_q20q30[i][0],fqInfo->clean_q20q30[i][1]);
        for (j=0; j<=fqInfo->maxQualityValue; j++)
        {
            err_printf("\t%lu", fqInfo->clean_qual[i][j]);
        }
        err_fputs("#S\n", stdout);
    }
}

void report_print(const FqInfo *fqInfo1, const FqInfo *fqInfo2)
{
    err_fputs("#Fq1_statistical_information\t#S\n", stdout);
    print_fq_info(fqInfo1);
    if (fqInfo2 != NULL){
        err_fputs("#Fq2_statistical_information\t#S\n", stdout);
        print_fq_info(fqInfo2);
    }
}
