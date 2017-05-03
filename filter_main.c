//
// Created by 黄志博 on 2017/5/2.
//
#include <getopt.h>
#include "filter.h"
#include "ksort.h"

static int usage();

int main_filter(int argc, char *argv[]);

int main_filter(int argc, char **argv) {
    int c;
    filter_opt_t *filter_opt = filter_opt_init();

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
            {0, 0, 0, 0}
    };

    int option_index = 0;
    while(1) {
        c = getopt_long(argc, argv, "pt:h:",
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
                filter_opt->is_pe = 1;
                break;
            case 't':
                filter_opt->n_threads = atoi(optarg);
                break;
            default:break;
        }
    }

    if (optind + 1 >= argc) {
        return usage();
    }

    return 0;
}

int usage() {
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
    fprintf(stderr, "\n");
    fprintf(stderr, "Note: Please read the man page for detailed description of the command line and options.\n");
    fprintf(stderr, "\n");
    return 1;
}
