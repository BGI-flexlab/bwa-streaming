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
    int skip_filter;
    int is_phred64;
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
    int max_quality_value;  //记录fq文件中碱基的最大质量值.
    int max_raw_read_len;     //raw data 读长
    int max_clean_read_len;   //clean data 读长

    unsigned long raw_base_num;  //raw data total base number
    unsigned long clean_base_num; //clean data total base number
    unsigned long raw_base_A;     //raw data base A number
    unsigned long clean_base_A;   //clean data base A number
    unsigned long raw_base_C;     //raw data base C number
    unsigned long clean_base_C; //clean data base C number
    unsigned long raw_base_G; //raw data base G number
    unsigned long clean_base_G; //clean data base G number
    unsigned long raw_base_T; //raw data base T number
    unsigned long clean_base_T; //clean data base T number
    unsigned long raw_base_N; //raw data base N number
    unsigned long clean_base_N; //clean data base N number
    unsigned long raw_q20; //rawfq文件中碱基质量>=20的碱基总数
    unsigned long clean_q20; //cleanfq文件中碱基质量>=20的碱基总数
    unsigned long raw_q30; //rawfq文件中碱基质量>=30的碱基总数
    unsigned long clean_q30; //cleanfq文件中碱基质量>=30的碱基总数

    unsigned long adapter_num;  //the number of read which contain adapter in raw data
    unsigned long n_exceed_num;  //the number of read which n rate was exceed in raw data
    unsigned long low_qual_num;  //low qualtiy read number in raw data
    unsigned long low_mean_num;  //low mean quality read number in raw data
    unsigned long small_insert_num;  //samll inert number in raw data  todo
    unsigned long polyA_num;    //polyA number in raw data
    unsigned long short_length_n;    //the read is too short

    unsigned long clean_read_len_distribution[MAX_LENGTH];
    //base distributions by read position(Raw)
    unsigned long base[MAX_LENGTH][5]; //ACGTN
    unsigned long clean_base[MAX_LENGTH][5];

    //Distribution ofQ20+/Q30+ bases by read position(Raw)
    unsigned long q20q30[MAX_LENGTH][2];  //Q20 Q30
    unsigned long clean_q20q30[MAX_LENGTH][2];

    //Basequality value distribution by read position(Raw)
    unsigned long qual[MAX_LENGTH][MAX_QUALITY + 1];
    unsigned long clean_qual[MAX_LENGTH][MAX_QUALITY + 1];
} FqInfo;

typedef struct{
    unsigned long total_raw_read_num;  //raw data read个数
    unsigned long total_clean_read_num;  //clean data read个数
    unsigned long total_short_length_num;  //长度过短的 read 数目
    unsigned long total_n_exceed_num;
    unsigned long total_low_qual_num;
    unsigned long total_low_mean_num;
    unsigned long total_adapter_num;
    unsigned long total_cut_adapter_num;
    unsigned long total_small_insert_num;
    unsigned long total_polyA_num;

    FqInfo *read1_info;
    FqInfo *read2_info;

    int id;
    char* rg_id;
}read_info_t;

typedef struct {
    const filter_opt_t *filter_opt;
    bseq1_t *seqs;
    int64_t n_processed;
    read_info_t ***read_info;
} filter_worker_t;


typedef struct {
//    ktp_aux_t *aux;
    int n_seqs;
    bseq1_t *seqs;
} ktp_data_t;

filter_opt_t* filter_opt_init();

read_info_t *read_info_init(int is_pe);
void read_info_destroy(read_info_t *read_indo);

void soapnuke_filter(const filter_opt_t *opt, int64_t n_processed, int n, bseq1_t *seqs, read_info_t *read_info[]);

void remove_bad_reads(ktp_data_t *w);

int statistics_se(bseq1_t *read1, const filter_opt_t *opt, read_info_t *read_info[]);

int statistics_pe(bseq1_t *read1, bseq1_t *read2, const filter_opt_t *opt, read_info_t *read_info[]);

void seq_stat(bseq1_t *read, const filter_opt_t *opt, int head_trim_n, int tail_trim_n, FqInfo *info, StatisInfo *si);

/**
 * 用于adapter为序列的情况
 * -1: no adapter  -2: filter adapter due to the adapter is too long  >0: adapter index to trim
 */
int adapter_align(bseq1_t *read, const char *adapter, const filter_opt_t *opt) ;

/**
 * 用于adapter list情况 todo
 */
//int hasAdapter(set<string> &readsName, const char *seqName);

/**
 * 统计碱基的分布
 */
void calculate_base_distribution(bseq1_t* read, FqInfo *info);

void report_print(read_info_t *read_info);

void merge_report(read_info_t *info, read_info_t *ptr);


#endif //BWA_FILTER_H
