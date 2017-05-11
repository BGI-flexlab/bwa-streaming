//
// Created by 黄志博 on 2017/5/2.
//
#include <getopt.h>
#include "filter.h"
#include "ksort.h"
#include "utils.h"
#include "kseq.h"
#include "bwa.h"

KSEQ_DECLARE(gzFile)

#define CHUNK_SIZE 10000000

static int usage();

int main_filter(int argc, char *argv[]);

void *kopen(const char *fn, int *_fd);
int kclose(void *a);

void kt_pipeline(int n_threads, void *(*func)(void*, int, void*), void *shared_data, int n_steps);

typedef struct {
    kseq_t *ks, *ks2;
    const filter_opt_t *opt;
    int64_t n_processed;
    int copy_comment, actual_chunk_size;
    FqInfo *fq_info;
} ktp_aux_t;

static void *process(void *shared, int step, void *_data)
{
    ktp_aux_t *aux = (ktp_aux_t*)shared;
    ktp_data_t *data = (ktp_data_t*)_data;
    int i;
    if (step == 0) {
        ktp_data_t *ret;
        int64_t size = 0;
        ret = calloc(1, sizeof(ktp_data_t));
//		ret->seqs = bseq_read(aux->actual_chunk_size, &ret->n_seqs, aux->ks, aux->ks2);
        ret->seqs = bseq_read2(aux->actual_chunk_size, &ret->n_seqs, aux->ks, aux->opt->is_pe, aux->opt->is_phred64);
        if (ret->seqs == 0) {
            free(ret);
            return 0;
        }

        if (!aux->copy_comment)
            for (i = 0; i < ret->n_seqs; ++i) {
                free(ret->seqs[i].comment);
                ret->seqs[i].comment = 0;
            }
        for (i = 0; i < ret->n_seqs; ++i) size += ret->seqs[i].l_seq;
        if (bwa_verbose >= 3)
            fprintf(stderr, "[M::%s] read %d sequences (%ld bp)...\n", __func__, ret->n_seqs, (long)size);
        return ret;
    } else if (step == 1) {
        const filter_opt_t *opt = aux->opt;
        soapnuke_filter(opt, aux->n_processed, data->n_seqs, data->seqs, aux->fq_info);
        remove_bad_reads(data);
        if (bwa_verbose >= 3) {
            int64_t size = 0;
            for (i = 0; i < data->n_seqs; ++i) size += data->seqs[i].l_seq;
            fprintf(stderr, "[M::%s] pass %d clean sequences (%ld bp)...\n", __func__, data->n_seqs, (long) size);
        }
        aux->n_processed += data->n_seqs;
        return data;
    } else if (step == 2) {
        if(aux->opt->is_pe){
            int n = data->n_seqs>>1;
            for (i = 0; i < n; ++i) {
                err_fputc('@',stdout);
                err_fputs(data->seqs[i<<1|0].name, stdout);
                err_fputs("/1\n",stdout);
                err_fputs(data->seqs[i<<1|0].seq, stdout);
                err_fputc('\n',stdout);
                err_fputs("+\n", stdout);
                err_fputs(data->seqs[i<<1|0].qual, stdout);
                err_fputc('\n',stdout);
                free(data->seqs[i<<1|0].name); free(data->seqs[i<<1|0].comment);
                free(data->seqs[i<<1|0].seq); free(data->seqs[i<<1|0].qual); free(data->seqs[i<<1|0].sam);

                err_fputc('@',stdout);
                err_fputs(data->seqs[i<<1|1].name, stdout);
                err_fputs("/2\n",stdout);
                err_fputs(data->seqs[i<<1|1].seq, stdout);
                err_fputc('\n',stdout);
                err_fputs("+\n", stdout);
                err_fputs(data->seqs[i<<1|1].qual, stdout);
                err_fputc('\n',stdout);
                free(data->seqs[i<<1|1].name); free(data->seqs[i<<1|1].comment);
                free(data->seqs[i<<1|1].seq); free(data->seqs[i<<1|1].qual); free(data->seqs[i<<1|1].sam);
            }

        } else{
            int n = data->n_seqs;
            for (i = 0; i < n; ++i) {
                err_fputc('@',stdout);
                err_fputs(data->seqs[i].name, stdout);
                err_fputc('\n',stdout);
                err_fputs(data->seqs[i].seq, stdout);
                err_fputc('\n',stdout);
                err_fputs("+\n", stdout);
                err_fputs(data->seqs[i].qual, stdout);
                err_fputc('\n',stdout);
                free(data->seqs[i].name); free(data->seqs[i].comment);
                free(data->seqs[i].seq); free(data->seqs[i].qual); free(data->seqs[i].sam);
            }
        }
        free(data->seqs); free(data);
        return 0;
    }
    return 0;
}

int main_filter(int argc, char **argv)
{
    int c, fd;
    gzFile fp;
    void *ko = 0;
    filter_opt_t *filter_opt = filter_opt_init();
    ktp_aux_t aux;
    memset(&aux, 0, sizeof(ktp_aux_t));

    static struct option long_options[] = {
            {"adapter1", required_argument, 0, 0},
            {"adapter2", required_argument, 0, 0},
            {"misMatch", required_argument, 0, 0},
            {"matchRatio", required_argument, 0, 0},
            {"cutAdaptor", no_argument, 0, 0},
            {"lowQual", required_argument, 0, 0},
            {"qualRate", required_argument, 0, 0},
            {"nRate", required_argument, 0, 0},
            {"mean", required_argument, 0, 0},
            {"trim", required_argument, 0, 0},
            {"minLen", required_argument, 0, 0},
            {"small", no_argument, 0, 0},
            {"tile", required_argument, 0, 0},
            {"polyA", required_argument, 0, 0},
            {"polyAType", required_argument, 0, 0},
            {"thread", required_argument, 0, 't'},
            {"is_pe", required_argument, 0, 'p'},
            {"phred64", required_argument, 0, 'i'},
            {0, 0, 0, 0}
    };

    int option_index = 0;
    while(1) {
        c = getopt_long(argc, argv, "pt:hv:",
                            long_options, &option_index);
        if (c == -1)
            break;

        switch(c) {
            case 0:
                if(strcmp(long_options[option_index].name, "adapter1")==0) filter_opt->adp1 = optarg;
                else if(strcmp(long_options[option_index].name, "adapter2")==0) filter_opt->adp2 = optarg;
                else if(strcmp(long_options[option_index].name, "misMatch")==0) filter_opt->misMatch = atoi(optarg);
                else if(strcmp(long_options[option_index].name, "matchRatio")==0) filter_opt->matchRatio = (float) atof(optarg);
                else if(strcmp(long_options[option_index].name, "cutAdaptor")==0) filter_opt->cutAdapter = 1;
                else if(strcmp(long_options[option_index].name, "lowQual")==0) filter_opt->lowQual = atoi(optarg);
                else if(strcmp(long_options[option_index].name, "qualRate")==0) filter_opt->qualRate = (float) atof(optarg);
                else if(strcmp(long_options[option_index].name, "nRate")==0) filter_opt->nRate = (float) atof(optarg);
                else if(strcmp(long_options[option_index].name, "mean")==0) filter_opt->mean = atoi(optarg);
                else if(strcmp(long_options[option_index].name, "minLen")==0) filter_opt->min_read_len = atoi(optarg);
                else if(strcmp(long_options[option_index].name, "small")==0) filter_opt->small = atoi(optarg);
                else if(strcmp(long_options[option_index].name, "tile")==0) filter_opt->tile = optarg;
                else if(strcmp(long_options[option_index].name, "polyA")==0) filter_opt->polyA = (float) atof(optarg);
                else if(strcmp(long_options[option_index].name, "polyAType")==0) filter_opt->polyAType = atoi(optarg);
                else if(strcmp(long_options[option_index].name, "trim")==0){
                    int index = 1;
                    char *trim_opt = optarg;
                    filter_opt->trim[0] = atoi(trim_opt);
                    while (*trim_opt) {
                        if(*trim_opt == ',') {
                            filter_opt->trim[index++] = atoi(trim_opt+1);
                            if(index>=4) break;
                        }
                        trim_opt++;
                    }
                }
                break;
            case 'p':
                filter_opt->is_pe = 0;
                break;
            case 't':
                filter_opt->n_threads = atoi(optarg);
                break;
            case 'i':
                filter_opt->is_phred64 = 1;
                break;
            case 'v':
                bwa_verbose = atoi(optarg);
                break;
            case 'h':
                usage();
            default:break;
        }
    }

    if (optind + 1 > argc) {
        return usage();
    }

    ko = kopen(argv[optind], &fd);
    if (ko == 0) {
        if (bwa_verbose >= 1) fprintf(stderr, "[E::%s] fail to open file `%s'.\n", __func__, argv[optind + 1]);
        return 1;
    }
    fp = gzdopen(fd, "r");
    aux.ks = kseq_init(fp);
    aux.opt = filter_opt;
    aux.actual_chunk_size = CHUNK_SIZE;
    aux.fq_info = fq_info_init();

    kt_pipeline(2, process, &aux, 3);

    report_print(aux.fq_info, aux.fq_info+1);

    free(filter_opt);
    free(aux.fq_info);
    kseq_destroy(aux.ks);
    err_gzclose(fp); kclose(ko);
    return 0;
}

int usage()
{
    fprintf(stderr, "\n");
    fprintf(stderr, "Usage: bwa filter [options] -\n");
    fprintf(stderr, "\nFilter options:\n\n");
    fprintf(stderr, "       --adapter1    STR      3' adapter sequence of fq1 file  [null]\n");
    fprintf(stderr, "       --adapter2    STR      5' adapter sequence of fq2 file (only for PE reads)  [null]\n");
    fprintf(stderr, "       --misMatch    INT      the max mismatch number when match the adapter [1]\n");
    fprintf(stderr, "       --matchRatio  FLOAT    adapter's shortest match ratio [0.5]\n");
    fprintf(stderr, "       --cutAdaptor           cut adaptor sequence according to --adapter1/--adapter2, just filter if is off [off]\n");
    fprintf(stderr, "       --lowQual     INT      low quality threshold  [5]\n");
    fprintf(stderr, "       --qualRate    FLOAT    low quality rate  [0.5]\n");
    fprintf(stderr, "       --nRate       FLOAT    N rate threshold  [0.05]\n");
    fprintf(stderr, "       --mean        FLOAT    filter reads with low average quality, (<)  [0]\n");
    fprintf(stderr, "       --minLen      INT      min length of reads [20]\n");
    fprintf(stderr, "       --trim        INT,INT,INT,INT\n");
    fprintf(stderr, "                              trim some bp of the read's head and tail, they means: (read1's head and tail and read2's head and tail) [0,0,0,0]\n");
    fprintf(stderr, "       --small                remove small insert. [off]\n");
    fprintf(stderr, "       --tile        STR      tile number to ignore reads, such as 1101-1104,1205. [null]\n");
    fprintf(stderr, "       --polyA       FLOAT    filter poly A, percent of A, 0 means do not filter. [0]\n");
    fprintf(stderr, "       --polyAType   INT      filter poly A type, 0->both two reads are poly a, 1->at least one reads is poly a, then filter. [0]\n");
    fprintf(stderr, "       --is_se       INT      single end read [off]\n");
    fprintf(stderr, "    -t,--thread      INT      number of threads [1]\n");
    fprintf(stderr, "    -i,--phred64     INT      set quality system as phred64, default is phred33 [off]\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Note: Please read the man page for detailed description of the command line and options.\n");
    fprintf(stderr, "\n");
    exit(1);
}
