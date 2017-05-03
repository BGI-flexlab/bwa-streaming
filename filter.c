//
// Created by huangzhibo on 2017/4/16.
//

#include <strings.h>
#include <math.h>
#include "filter.h"
#include "ksort.h"
#include "utils.h"

//// filter step (SOAPnuke)
static void worker(void *data, int i, int tid)
{
	filter_worker_t *w = (filter_worker_t*)data;
	if (!(w->filter_opt->is_pe)) {
        if (bwa_verbose >= 4) printf("=====> Processing read '%s' <=====\n", w->seqs[i].name);
//        w->is_clean = statistics_se(w->seqs[i<<1|0],w->seqs[i<<1|1], w->filter_opt, w->aux[tid]);
    } else{
        if (bwa_verbose >= 4) printf("=====> Processing read '%s'/1/2 <=====\n", w->seqs[i<<1|0].name);
        w->is_clean = statistics_pe(&w->seqs[i<<1|0], &w->seqs[i<<1|1], w->filter_opt, w->fqInfo);
    }
}

void soapnuke_filter(const filter_opt_t *opt,int64_t n_processed, int n, bseq1_t *seqs);

void soapnuke_filter(const filter_opt_t *opt, int64_t n_processed, int n, bseq1_t *seqs, FqInfo fqInfo[2]) {
    extern void kt_for(int n_threads, void (*func)(void*,int,int), void *data, int n);
    filter_worker_t w;
    mem_pestat_t pes[4];
    double ctime, rtime;
    int i;

    ctime = cputime(); rtime = realtime();
    w.filter_opt = opt;
    w.seqs = seqs;
    w.n_processed = n_processed;
    w.fqInfo = fqInfo; //fixme

    kt_for(opt->n_threads, worker, &w, (opt->is_pe)? n>>1 : n); // generate alignment
    if (bwa_verbose >= 3)
        fprintf(stderr, "[M::%s] Processed %d reads in %.3f CPU sec, %.3f real sec\n", __func__, n, cputime() - ctime, realtime() - rtime);
}

int statistics_pe(bseq1_t *read1, bseq1_t *read2, filter_opt_t *opt, FqInfo *info[2])
{
    StatisInfo si1, si2;
    int trim_tail1 = 0;
    int trim_tail2 = 0;

    int index1 = adapter_align(read1, opt->adp1, opt);
    int index2 = adapter_align(read2, opt->adp2, opt);

    if(index1 == -2){
        si1.hasAdpt = 1;
    } else if(index1 >= 0) {
        int cut_adapter_len = read1->l_seq - index1;
        trim_tail1 = cut_adapter_len > opt->trim[1] ? cut_adapter_len : opt->trim[1];
    }

    if(index2 == -2){
        si2.hasAdpt = 1;
    } else if(index2 >= 0){
        int cut_adapter_len = read2->l_seq - index2;
        trim_tail2 = cut_adapter_len > opt->trim[3] ? cut_adapter_len : opt->trim[3];
        info[1]->totalCutAdaptorNum++;
    }

    //fq1
    seq_stat(read1, opt, opt->trim[0], trim_tail1, info[0], &si1);
    //fq2
    seq_stat(read2, opt, opt->trim[2], trim_tail2, info[1], &si2);

    // filter short length read
    if (read1->l_seq < opt->min_read_len || read2->l_seq < opt->min_read_len){
        if(read1->l_seq < opt->min_read_len){
            info[0]->short_length_n++;
        }
        if(read2->l_seq < opt->min_read_len){
            info[1]->short_length_n++;
        }
        info[0]->total_short_length_n++;
        return 0;
    }

    // filter adapter read
    if(si1.hasAdpt || si2.hasAdpt){
        if (si1.hasAdpt){
            info[0]->adapterNum++;
        }
        if (si2.hasAdpt){
            info[1]->adapterNum++;
        }
        info[0]->totalAdapterNum++;
        return 0;
    }

    if(si1.nExceed || si2.nExceed){
        if (si1.nExceed){
            info[0]->nExceedNum++;
        }
        if (si2.nExceed) {
            info[1]->nExceedNum++;
        }
        info[0]->totalNExceedNum++;
        return 0;
    }

    if(si1.isLowQual || si2.isLowQual){
        if (si1.isLowQual){
            info[0]->lowQualNum++;
        }
        if (si2.isLowQual) {
            info[1]->lowQualNum++;
        }
        info[0]->totalLowQualNum++;
        return 0;
    }

    if(si1.sumQuality < opt->mean * read1->l_seq || si2.sumQuality < opt->mean * read2->l_seq){
        if(si1.sumQuality < opt->mean * read1->l_seq){
            info[0]->lowMeanNum++;
        }
        if(si2.sumQuality < opt->mean * read2->l_seq){
            info[1]->lowMeanNum++;
        }
        info[0]->totalLowMeanNum++;
        return 0;
    }

//    // todo polyA
//    if ((polyAType_==0 && (!si1.isPolyA || !si2.isPolyA)) ||(polyAType_==1 && !si1.isPolyA && !si2.isPolyA)) //not polyA
//    {
//    }
//    else //has polyA
//    {
//        if(polyAType_==0){
//            info[0]->polyANum++;
//            info[1]->polyANum++;
//        }else{
//            if(si1.isPolyA)
//                info[0]->polyANum++;
//            if(si2.isPolyA)
//                info[1]->polyANum++;
//        }
//        info[0]->totalPolyANum++;
//    }

    info[0]->cleanBaseA += si1.a;
    info[0]->cleanBaseC += si1.c;
    info[0]->cleanBaseG += si1.g;
    info[0]->cleanBaseT += si1.t;
    info[0]->cleanBaseN += si1.n;
    info[0]->cleanQ20 += si1.q20;
    info[0]->cleanQ30 += si1.q30;
    info[0]->cleanTotalReadNum++;
    info[0]->cleanTotalBaseNum += read1->l_seq;
    calculate_base_distribution(read1, info[0]);

    info[1]->cleanBaseA += si2.a;
    info[1]->cleanBaseC += si2.c;
    info[1]->cleanBaseG += si2.g;
    info[1]->cleanBaseT += si2.t;
    info[1]->cleanBaseN += si2.n;
    info[1]->cleanQ20 += si2.q20;
    info[1]->cleanQ30 += si2.q30;
    info[1]->cleanTotalReadNum++;
    info[1]->cleanTotalBaseNum += read2->l_seq;
    calculate_base_distribution(read2, info[1]);

    return 1;
}

void calculate_base_distribution(bseq1_t *read, FqInfo *info) {
    int qual;

    for (int i=0; i<read->l_seq; ++i) {
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

void seq_stat(bseq1_t *read, filter_opt_t *opt, int head_trim_n, int tail_trim_n, FqInfo *info, StatisInfo *si) {
    int qual;

    info->rawTotalReadNum ++;
    info->rawTotalBaseNum += read->l_seq;

    //todo 切除低质量末端

    int a = 0, g = 0, c = 0, t = 0, n = 0;
    int q20 = 0, q30 = 0;

    int right = read->l_seq - tail_trim_n;
    int sumQual = 0;

    for (int i=0; i<read->l_seq; ++i)
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

        qual = read->qual[i];

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
int adapter_align(bseq1_t *read, const char *adapter, filter_opt_t *opt) {
    int find = -1;
    int adptLen = (int) strlen(adapter);
    int minMatchLen = (int) ceilf(adptLen * opt->matchRatio);
    int a1 = adptLen - minMatchLen;
    int r1 = 0;
    int len, mis;

    int right = read->l_seq - minMatchLen;

    for (r1 = 0; r1 <= right;)
    {
        int len1 = adptLen - a1;
        int len2 = read->l_seq - r1;
        len = (len1 < len2) ? len1 : len2;
        mis = 0;
        int map[MAX_LENGTH];
        map[0] = 0;
        for (int c = 0; c < len; ++c)
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
        for (int c = 0; c <= mis; ++c)
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

filter_opt_t *filter_opt_init() {
    filter_opt_t *o;
    o = calloc(1, sizeof(filter_opt_t));
    o->is_phred64 = 1;
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
    for(int i=0; i < 4; i++){
        o->trim[i] = 0;
    }
    return o;
}
