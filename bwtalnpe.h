#include "bwtaln.h"
#include "bwtgap.h"
#include "utils.h"
#include "bwase.h"

int split_sample_list(char *str,char *delim, smaple_list *sl[]);

int read_sample_list(char * file, smaple_list *sl[]);

void bwa_print_sam_SQs(const bntseq_t *bns, smaple_list *sl[], int rg_number);
