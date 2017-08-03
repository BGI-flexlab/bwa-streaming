#include <stdio.h>
#include <unistd.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <stdint.h>
#include "bwtaln.h"
#include "bntseq.h"
#include "utils.h"
#include "bwape.h"
#include "readLine.h"
#include "bwtalnpe.h"
#ifdef HAVE_PTHREAD
#include <pthread.h>
#endif
#ifdef USE_MALLOC_WRAPPERS
#  include "malloc_wrap.h"
#endif


#define MIN_HASH_WIDTH 1000
#define b128_eq(a, b) ((a).x == (b).x && (a).y == (b).y)
#define b128_hash(a) ((uint32_t)(a).x)
#include "khash.h"
KHASH_INIT(b128, pair64_t, poslist_t, 1, b128_hash, b128_eq)
extern int g_log_n[256]; // in bwase.c
static kh_b128_t *g_hash;

#ifdef HAVE_PTHREAD
typedef struct {
	int tid;
	bwt_t *bwt;
	int n_seqs;
	bwa_seq_t *seqs;
	const gap_opt_t *opt;
} thread_aux_t;

static void *worker(void *data)
{
	thread_aux_t *d = (thread_aux_t*)data;
	bwa_cal_sa_reg_gap2(d->tid, d->bwt, d->n_seqs, d->seqs, d->opt);
	return 0;
}
#endif

int split_sample_list(char *str,char *delim, smaple_list *sl[])
{
	char *token[2];
	int i = 0, id;
	while(*str)
	{
		token[i] = str;
		if(i == 2) 
			break;
		while(*str != *delim && *str != '\0')
		{
			str++;
		}
		*str = '\0';
		str++;
		i++;
	}
	
	if(token[0] != NULL)
	{
		id = atoi(token[0]);
	} else {
		printf("null token!\n");
		return -1;
	}
	sl[id] = (smaple_list *) malloc(1*sizeof(smaple_list));
	sl[id]->id = id;
	if(token[1] != NULL)
	{
		sl[id]->rg = bwa_set_rg(token[1]);
		if(!sl[id]->rg)
		{
			return -1;
		}
		
		char *p, *q, *r;
		p = strstr(sl[id]->rg, "\tID:");
		if (p == 0) return -1;
		p += 4;
		for (q = p; *q && *q != '\t' && *q != '\n'; ++q);
		sl[id]->rg_id = calloc(q - p + 1, 1);
		for (q = p, r = sl[id]->rg_id; *q && *q != '\t' && *q != '\n'; ++q)
			*r++ = *q;
	    return 0;
	} else {
		return -1;
	}
	//fprintf(stderr, "ID:%d\tRG:%s\n", id, sl[id]->rg);
}

int read_sample_list(char * file, smaple_list *sl[]) 
{
	int list_fd;
	
	list_fd = open(file, O_RDONLY);
	if(list_fd < 0)
	{
		return -1;
	}
	char line[4096];
	int len;
	int i = 0;
	while((len = Readline(list_fd, line, 4096)) > 0)
	{
		if(len == 1)
			continue;
		if(split_sample_list(line, "\t", sl) == -1)
		{
			fprintf(stderr, "wrong line:%s", line);
			return -1;
		}
		i++;
	}
	close(list_fd);
	return i;
}

void bwa_alnpe_core(const char *prefix, const char *fn_fa, const gap_opt_t *opt, pe_opt_t *popt, const char *rg_line)
{
	extern bwa_seqio_t *bwa_open_reads(int mode, const char *fn_fa);
	int i, j, n_seqs, tot_seqs = 0;
	bwa_seq_t *seqs[2];
	bwa_seqio_t *ks;
	clock_t t;
	bntseq_t *bns;
	khint_t iter;
	isize_info_t last_ii; // this is for the last batch of reads
	bwt_t *bwt;
	uint8_t *pac;

	// initialization
	bwase_initialize(); // initialize g_log_n[] in bwase.c
	pac = 0; bwt = 0;
	for (i = 1; i != 256; ++i) g_log_n[i] = (int)(4.343 * log(i) + 0.5);
	bns = bns_restore(prefix);
	srand48(bns->seed);
	//fp_sa[0] = xopen(fn_sa[0], "r");
	//fp_sa[1] = xopen(fn_sa[1], "r");
	g_hash = kh_init(b128);
	last_ii.avg = -1.0;
	
	//open fa
	ks = bwa_open_reads(opt->mode, fn_fa);

	{ // load BWT
		char *str = (char*)calloc(strlen(prefix) + 10, 1);
		strcpy(str, prefix); strcat(str, ".bwt");  bwt = bwt_restore_bwt(str);
		strcpy(str, prefix); strcat(str, ".sa"); bwt_restore_sa(str, bwt);
		free(str);
		// for Illumina alignment only
		if (popt->is_preload) {
			pac = (ubyte_t*)calloc(bns->l_pac/4+1, 1);
			err_rewind(bns->fp_pac);
			err_fread_noeof(pac, 1, bns->l_pac/4+1, bns->fp_pac);
		}
	}

	//header
	int rg_number = -1;
	smaple_list **sl = NULL;
	if(popt->sample_list != NULL) 
	{
		sl = (smaple_list**)malloc(10000*sizeof(smaple_list *));
		int s;
		for (s =0; s < 10000; s++)
		{
			sl[s] = NULL;
		}
		rg_number = read_sample_list(popt->sample_list, sl);
		fprintf(stderr, "rg number : %d\n", rg_number);
		free(popt->sample_list);
	}
	if(rg_number < 0){
		char *hdr_line = 0;
		if (rg_line) {
			hdr_line = bwa_insert_header(rg_line, hdr_line);
		}
		bwa_print_sam_hdr(bns, hdr_line);
	}
	else{
		bwa_print_sam_hdr2(bns, sl, rg_number);
	}	
	
	// core loop
	while ((bwa_read_seq_pe(ks, 0x40000, &n_seqs, opt->mode, opt->trim_qual, seqs)) != 0) {
		int cnt_chg;
		isize_info_t ii;
		ubyte_t *pacseq;
		
		tot_seqs += 2*n_seqs;
		t = clock();
		
		fprintf(stderr, "[bwa_aln_core] calculate SA coordinate... ");
		for (i = 0; i < 2; ++i)
		{
#ifdef HAVE_PTHREAD
			if (opt->n_threads <= 1) { // no multi-threading at all
				//fprintf(stderr, "call fq %d\n", i);
				bwa_cal_sa_reg_gap2(0, bwt, n_seqs, seqs[i], opt);
			} else {
				pthread_t *tid;
				pthread_attr_t attr;
				thread_aux_t *data;
				int j;
				pthread_attr_init(&attr);
				pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
				data = (thread_aux_t*)calloc(opt->n_threads, sizeof(thread_aux_t));
				tid = (pthread_t*)calloc(opt->n_threads, sizeof(pthread_t));
				for (j = 0; j < opt->n_threads; ++j) {
					data[j].tid = j; data[j].bwt = bwt;
					data[j].n_seqs = n_seqs; data[j].seqs = seqs[i]; data[j].opt = opt;
					pthread_create(&tid[j], &attr, worker, data + j);
				}
				for (j = 0; j < opt->n_threads; ++j) pthread_join(tid[j], 0);
				free(data); free(tid);
			}
#else
			bwa_cal_sa_reg_gap2(0, bwt, n_seqs, seqs[i], opt);
#endif
		}
		fprintf(stderr, "%.2f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC);

		t = clock();
		//PE
		fprintf(stderr, "[bwa_sai2sam_pe_core] convert to sequence coordinate... \n");
		cnt_chg = bwa_cal_pac_pos_pe2(bns, prefix, bwt, n_seqs, seqs, &ii, popt, opt, &last_ii);
		fprintf(stderr, "[bwa_sai2sam_pe_core] time elapses: %.2f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC); t = clock();
		fprintf(stderr, "[bwa_sai2sam_pe_core] changing coordinates of %d alignments.\n", cnt_chg);

		fprintf(stderr, "[bwa_sai2sam_pe_core] align unmapped mate...\n");
		pacseq = bwa_paired_sw(bns, pac, n_seqs, seqs, popt, &ii);
		fprintf(stderr, "[bwa_sai2sam_pe_core] time elapses: %.2f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC); t = clock();

		fprintf(stderr, "[bwa_sai2sam_pe_core] refine gapped alignments... ");
		for (j = 0; j < 2; ++j)
			bwa_refine_gapped(bns, n_seqs, seqs[j], pacseq);
		fprintf(stderr, "%.2f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC); t = clock();
		if (pac == 0) free(pacseq);

		fprintf(stderr, "[bwa_sai2sam_pe_core] print alignments... ");
		for (i = 0; i < n_seqs; ++i) {
			bwa_seq_t *p[2];
			p[0] = seqs[0] + i; p[1] = seqs[1] + i;
			if (p[0]->bc[0] || p[1]->bc[0]) {
				strcat(p[0]->bc, p[1]->bc);
				strcpy(p[1]->bc, p[0]->bc);
			}
			if(rg_number < 0) 
			{
				bwa_print_sam1(bns, p[0], p[1], opt->mode, opt->max_top2);
				bwa_print_sam1(bns, p[1], p[0], opt->mode, opt->max_top2);
			}
			else
			{
				//fprintf(stderr, "sample id:%d-%d\trg:%s\n", p[0]->sample_id, p[1]->sample_id, sl[p[0]->sample_id]->rg_id);
				char *rg_id = sl[p[0]->sample_id]->rg_id;
				bwa_print_sam2(bns, p[0], p[1], opt->mode, opt->max_top2, rg_id);
				rg_id = sl[p[1]->sample_id]->rg_id;
				bwa_print_sam2(bns, p[1], p[0], opt->mode, opt->max_top2, rg_id);
			}
			if (strcmp(p[0]->name, p[1]->name) != 0) err_fatal(__func__, "paired reads have different names: \"%s\", \"%s\"\n", p[0]->name, p[1]->name);
		}
		fprintf(stderr, "%.2f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC); t = clock();

		for (j = 0; j < 2; ++j)
			bwa_free_read_seq(n_seqs, seqs[j]);
		fprintf(stderr, "[bwa_sai2sam_pe_core] %d sequences have been processed.\n", tot_seqs);
		last_ii = ii;
	}

	// destroy
	bns_destroy(bns);
	bwa_seq_close(ks);
	bwt_destroy(bwt);
	for (iter = kh_begin(g_hash); iter != kh_end(g_hash); ++iter)
		if (kh_exist(g_hash, iter)) free(kh_val(g_hash, iter).a);
	kh_destroy(b128, g_hash);
	if (pac) {
		free(pac);
	}
	if(sl != NULL)
	{
		int j;
		for(j = 0; j < 10000; j++)
		{
			if(sl[j] != NULL)
			{
				free(sl[j]->rg);
				free(sl[j]->rg_id);
				free(sl[j]);
			}
		}
		free(sl);
	}
}	


int bwa_alnpe(int argc, char *argv[])
{
	int c, opte = -1;
	gap_opt_t *opt;
	pe_opt_t *popt;
	char *prefix, *rg_line = 0;

	opt = gap_init_opt();
	popt = bwa_init_pe_opt();
	while ((c = getopt(argc, argv, "n:o:e:i:d:l:k:LR:m:t:NM:O:E:q:f:b012IYB:a:g:sPh:H:p:f:Ar:T:")) >= 0) {
		switch (c) {
		case 'n':
			if (strstr(optarg, ".")) opt->fnr = atof(optarg), opt->max_diff = -1;
			else opt->max_diff = atoi(optarg), opt->fnr = -1.0;
			break;
		case 'o': opt->max_gapo = atoi(optarg); break;
		case 'e': opte = atoi(optarg); break;
		case 'M': opt->s_mm = atoi(optarg); break;
		case 'O': opt->s_gapo = atoi(optarg); break;
		case 'E': opt->s_gape = atoi(optarg); break;
		case 'd': opt->max_del_occ = atoi(optarg); break;
		case 'i': opt->indel_end_skip = atoi(optarg); break;
		case 'l': opt->seed_len = atoi(optarg); break;
		case 'k': opt->max_seed_diff = atoi(optarg); break;
		case 'm': opt->max_entries = atoi(optarg); break;
		case 't': opt->n_threads = atoi(optarg); break;
		case 'L': opt->mode |= BWA_MODE_LOGGAP; break;
		case 'R': opt->max_top2 = atoi(optarg); break;
		case 'q': opt->trim_qual = atoi(optarg); break;
		case 'N': opt->mode |= BWA_MODE_NONSTOP; opt->max_top2 = 0x7fffffff; break;
		case 'b': opt->mode |= BWA_MODE_BAM; break;
		case '0': opt->mode |= BWA_MODE_BAM_SE; break;
		case '1': opt->mode |= BWA_MODE_BAM_READ1; break;
		case '2': opt->mode |= BWA_MODE_BAM_READ2; break;
		case 'I': opt->mode |= BWA_MODE_IL13; break;
		case 'Y': opt->mode |= BWA_MODE_CFY; break;
		case 'B': opt->mode |= atoi(optarg) << 24; break;
		case 'T': popt->sample_list = strdup(optarg);break;
		case 'r':
			if ((rg_line = bwa_set_rg(optarg)) == 0) return 1;
			break;
		case 'a': popt->max_isize = atoi(optarg); break;
		case 'g': popt->max_occ = atoi(optarg); break;
		case 's': popt->is_sw = 0; break;
		case 'P': popt->is_preload = 1; break;
		case 'h': popt->n_multi = atoi(optarg); break;
		case 'H': popt->N_multi = atoi(optarg); break;
		case 'p': popt->ap_prior = atof(optarg); break;
		case 'f': xreopen(optarg, "w", stdout); break;
		case 'A': popt->force_isize = 1; break;
		default: return 1;
		}
	}
	if (opte > 0) {
		opt->max_gape = opte;
		opt->mode &= ~BWA_MODE_GAPE;
	}

	if (optind + 2 > argc) {
		fprintf(stderr, "\n");
		fprintf(stderr, "Usage:   bwa alnpe [options] <prefix> <1.fq 2.fq cat file>\n\n");
		fprintf(stderr, "aln Options: -n NUM    max #diff (int) or missing prob under %.2f err rate (float) [%.2f]\n",
				BWA_AVG_ERR, opt->fnr);
		fprintf(stderr, "         -o INT    maximum number or fraction of gap opens [%d]\n", opt->max_gapo);
		fprintf(stderr, "         -e INT    maximum number of gap extensions, -1 for disabling long gaps [-1]\n");
		fprintf(stderr, "         -i INT    do not put an indel within INT bp towards the ends [%d]\n", opt->indel_end_skip);
		fprintf(stderr, "         -d INT    maximum occurrences for extending a long deletion [%d]\n", opt->max_del_occ);
		fprintf(stderr, "         -l INT    seed length [%d]\n", opt->seed_len);
		fprintf(stderr, "         -k INT    maximum differences in the seed [%d]\n", opt->max_seed_diff);
		fprintf(stderr, "         -m INT    maximum entries in the queue [%d]\n", opt->max_entries);
		fprintf(stderr, "         -t INT    number of threads [%d]\n", opt->n_threads);
		fprintf(stderr, "         -M INT    mismatch penalty [%d]\n", opt->s_mm);
		fprintf(stderr, "         -O INT    gap open penalty [%d]\n", opt->s_gapo);
		fprintf(stderr, "         -E INT    gap extension penalty [%d]\n", opt->s_gape);
		fprintf(stderr, "         -R INT    stop searching when there are >INT equally best hits [%d]\n", opt->max_top2);
		fprintf(stderr, "         -q INT    quality threshold for read trimming down to %dbp [%d]\n", BWA_MIN_RDLEN, opt->trim_qual);
		fprintf(stderr, "         -B INT    length of barcode\n");
		fprintf(stderr, "         -L        log-scaled gap penalty for long deletions\n");
		fprintf(stderr, "         -N        non-iterative mode: search for all n-difference hits (slooow)\n");
		fprintf(stderr, "         -I        the input is in the Illumina 1.3+ FASTQ-like format\n");
		fprintf(stderr, "         -b        the input read file is in the BAM format\n");
		fprintf(stderr, "         -0        use single-end reads only (effective with -b)\n");
		fprintf(stderr, "         -1        use the 1st read in a pair (effective with -b)\n");
		fprintf(stderr, "         -2        use the 2nd read in a pair (effective with -b)\n");
		fprintf(stderr, "         -Y        filter Casava-filtered sequences\n");
		fprintf(stderr, "         -T STR    multi-sample list of Gaea\n");
		fprintf(stderr, "\n");
		fprintf(stderr, "sampe Options: -a INT   maximum insert size [%d]\n", popt->max_isize);
		fprintf(stderr, "         -g INT   maximum occurrences for one end [%d]\n", popt->max_occ);
		fprintf(stderr, "         -h INT   maximum hits to output for paired reads [%d]\n", popt->n_multi);
		fprintf(stderr, "         -H INT   maximum hits to output for discordant pairs [%d]\n", popt->N_multi);
		fprintf(stderr, "         -p FLOAT prior of chimeric rate (lower bound) [%.1le]\n", popt->ap_prior);
		fprintf(stderr, "         -r STR   read group header line such as `@RG\\tID:foo\\tSM:bar' [null]\n");
		fprintf(stderr, "         -P       preload index into memory (for base-space reads only)\n");
		fprintf(stderr, "         -s       disable Smith-Waterman for the unmapped mate\n");
		fprintf(stderr, "         -A       disable insert size estimate (force -s)\n\n");
		fprintf(stderr, "Notes: 1. For SOLiD reads, <in1.fq> corresponds R3 reads and <in2.fq> to F3.\n");
		fprintf(stderr, "       2. For reads shorter than 30bp, applying a smaller -o is recommended to\n");
		fprintf(stderr, "          to get a sensible speed at the cost of pairing accuracy.\n");
		fprintf(stderr, "\n");
		return 1;
	}
	if (opt->fnr > 0.0) {
		int i, k;
		for (i = 17, k = 0; i <= 250; ++i) {
			int l = bwa_cal_maxdiff(i, BWA_AVG_ERR, opt->fnr);
			if (l != k) fprintf(stderr, "[bwa_aln] %dbp reads: max_diff = %d\n", i, l);
			k = l;
		}
	}
	if ((prefix = bwa_idx_infer_prefix(argv[optind])) == 0) {
		fprintf(stderr, "[%s] fail to locate the index\n", __func__);
		free(opt);
		free(popt);
		return 1;
	}
	
	//bwa_aln_core(prefix, argv[optind+1], opt);
	bwa_alnpe_core(prefix, argv[optind+1], opt, popt, rg_line);
	free(opt); free(popt); free(prefix);
	return 0;
}
