//
// Created by huangzhibo on 2017/4/16.
//

#ifndef BWA_FILTER_H
#define BWA_FILTER_H

#include "bwa.h"
#include "bwtaln.h"
#include "bwamem.h"

#define MAX_LENGTH 1000
#define MAX_QUALITY 50

//统计read截断了两端后的信息
typedef struct
{
    int q20, q30;
    int a, c, g, t, n;
    int lowQual;

    int hasAdpt;  // non-zero if has adaptor seq
    int isLowQual; // non-zero if qual is low
    int nExceed; // non-zero if has too many N base
    int isPolyA; // non-zero if is polyA read
    int sumQuality; // for mean quality filter
}StatisInfo;

typedef struct {
    int is_pe;
    int n_threads;
    char *adp1;
    char *adp2;
    int misMatch;
    int cutAdapter;
    float matchRatio;
    float lowQual;
    float qualRate;
    float nRate;
    float mean;
    int min_read_len;
    int trim[4];
    int small;
    float polyA;
    int polyAType;
    char *tile;
} filter_opt_t;


//FqInfo结构用于存储过滤的所有统计信息。
typedef struct
{
    unsigned int rawReadLength;     //raw data 读长
    unsigned int cleanReadLength;   //clean data 读长
    unsigned long rawTotalReadNum;  //raw data read个数
    unsigned long cleanTotalReadNum;//clean data read 个数
    unsigned long rawTotalBaseNum;  //raw data total base number
    unsigned long cleanTotalBaseNum; //clean data total base number
    unsigned long rawBaseA;     //raw data base A number
    unsigned long cleanBaseA;   //clean data base A number
    unsigned long rawBaseC;     //raw data base C number
    unsigned long cleanBaseC; //clean data base C number
    unsigned long rawBaseG; //raw data base G number
    unsigned long cleanBaseG; //clean data base G number
    unsigned long rawBaseT; //raw data base T number
    unsigned long cleanBaseT; //clean data base T number
    unsigned long rawBaseN; //raw data base N number
    unsigned long cleanBaseN; //clean data base N number
    unsigned long rawQ20; //rawfq文件中碱基质量>=20的碱基总数
    unsigned long cleanQ20; //cleanfq文件中碱基质量>=20的碱基总数
    unsigned long rawQ30; //rawfq文件中碱基质量>=30的碱基总数
    unsigned long cleanQ30; //cleanfq文件中碱基质量>=30的碱基总数

    unsigned long adapterNum;  //the number of read which contain adapter in raw data
    unsigned long nExceedNum;  //the number of read which n rate was exceed in raw data
    unsigned long lowQualNum;  //low qualtiy read number in raw data
    unsigned long lowMeanNum;  //low mean quality read number in raw data
    unsigned long smallInsertNum;  //samll inert number in raw data
    unsigned long polyANum;    //polyA number in raw data
    unsigned long short_length_n;    //the read is too short

    unsigned long total_short_length_n;
    unsigned long totalAdapterNum;
    unsigned long totalNExceedNum;
    unsigned long totalLowQualNum;
    unsigned long totalLowMeanNum;
    unsigned long totalSmallInsertNum;
    unsigned long totalPolyANum;

    unsigned long totalCutAdaptorNum;

    unsigned long cleanReadLengthDistribution[MAX_LENGTH];
    //base distributions by read position(Raw)
    unsigned long base[MAX_LENGTH][5]; //ACGTN
    unsigned long clean_base[MAX_LENGTH][5];

    //Distribution ofQ20+/Q30+ bases by read position(Raw)
    unsigned long q20q30[MAX_LENGTH][2];  //Q20 Q30
    unsigned long clean_q20q30[MAX_LENGTH][2];

    //Basequality value distribution by read position(Raw)
    unsigned long qual[MAX_LENGTH][MAX_QUALITY + 1];
    unsigned long clean_qual[MAX_LENGTH][MAX_QUALITY + 1];

    int maxQualityValue;  //记录fq文件中碱基的最大质量值.

} FqInfo;

typedef struct {
    filter_opt_t *filter_opt;
    bseq1_t *seqs;
    int64_t n_processed;
    int is_clean;   // 0: should be filtered; 1: good
    FqInfo *fqInfo[2];
} filter_worker_t;

filter_opt_t* filter_opt_init();

//void soapnuke_filter(bseq1_t *reads, filter_opt_t opt, int n_seq);

int statistics_se(bseq1_t *read1, filter_opt_t *opt, FqInfo *info);

int statistics_pe(bseq1_t *read1, bseq1_t *read2, filter_opt_t *opt, FqInfo *info[2]);

void seq_stat(bseq1_t *read, filter_opt_t *opt, int head_trim_n, int tail_trim_n, FqInfo *info, StatisInfo *si);

/**
 * 用于adapter为序列的情况
 * -1: no adapter  -2: filter adapter due to the adapter is too long  >0: adapter index to trim
 */
int adapter_align(bseq1_t *read, const char *adapter, filter_opt_t *opt) ;

/**
 * 用于adapter list情况 todo
 */
//int hasAdapter(set<string> &readsName, const char *seqName);

/**
 * 统计碱基的分布
 */
void calculate_base_distribution(bseq1_t* read, FqInfo *info);

#endif //BWA_FILTER_H
