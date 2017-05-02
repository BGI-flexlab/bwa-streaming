#include <time.h>
#include <stdint.h>
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#include "bwa.h"
#include "bwtaln.h"
#include "bwtgap.h"
#include "utils.h"
#include "bwase.h"
#include "bwtalnpe.h"
#include "bntseq.h"
#include "kstring.h"
#include "ksw.h"

#ifdef HAVE_PTHREAD
#include <pthread.h>
#endif

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


void bwa_alnse_core(const char *prefix, const char *fn_fa, const gap_opt_t *opt, int n_occ, const char *rg_line, char *sample_list_file)
{
	extern bwa_seqio_t *bwa_open_reads(int mode, const char *fn_fa);
	bwt_t *bwt;
	int i, n_seqs, tot_seqs = 0;
	bwt_aln1_t *aln = 0;
	bwa_seq_t *seqs;
	bwa_seqio_t *ks;
	clock_t t;
	bntseq_t *bns;

	// initialization
	ks = bwa_open_reads(opt->mode, fn_fa);
	bwase_initialize();
	bns = bns_restore(prefix);
	srand48(bns->seed);

	{ // load BWT
		char *str = (char*)calloc(strlen(prefix) + 10, 1);
		strcpy(str, prefix); strcat(str, ".bwt");  bwt = bwt_restore_bwt(str);
		strcpy(str, prefix); strcat(str, ".sa"); bwt_restore_sa(str, bwt);
		free(str);
	}
	
	//header
	int rg_number = -1;
	smaple_list **sl = NULL;
	if(sample_list_file != NULL) 
	{
		sl = (smaple_list**)malloc(10000*sizeof(smaple_list *));
		int s;
		for (s =0; s < 10000; s++)
		{
			sl[s] = NULL;
		}
		rg_number = read_sample_list(sample_list_file, sl);
		free(sample_list_file);
	}
	if(rg_number < 0)
		bwa_print_sam_hdr(bns, rg_line);
	else
		bwa_print_sam_hdr2(bns, sl, rg_number);

	// core loop
	while ((seqs = bwa_read_seq_se(ks, 0x40000, &n_seqs, opt->mode, opt->trim_qual)) != 0) {
		tot_seqs += n_seqs;
		t = clock();

		fprintf(stderr, "[bwa_aln_core] calculate SA coordinate... ");

#ifdef HAVE_PTHREAD
		if (opt->n_threads <= 1) { // no multi-threading at all
			bwa_cal_sa_reg_gap2(0, bwt, n_seqs, seqs, opt);
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
				data[j].n_seqs = n_seqs; data[j].seqs = seqs; data[j].opt = opt;
				pthread_create(&tid[j], &attr, worker, data + j);
			}
			for (j = 0; j < opt->n_threads; ++j) pthread_join(tid[j], 0);
			free(data); free(tid);
		}
#else
		bwa_cal_sa_reg_gap2(0, bwt, n_seqs, seqs, opt);
#endif

		fprintf(stderr, "%.2f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC);
		
		//read alignment
		for (i = 0; i < n_seqs; ++i) {
			bwa_seq_t *p = seqs + i;
			int n_aln = p->n_aln;
			aln = p->aln;
			bwa_aln2seq_core(n_aln, aln, p, 1, n_occ);
		}

		t = clock();
		fprintf(stderr, "[bwa_aln_core] convert to sequence coordinate... ");
		bwa_cal_pac_pos2(bns, prefix, bwt, n_seqs, seqs, opt->max_diff, opt->fnr); // forward bwt will be destroyed here
		fprintf(stderr, "%.2f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC); t = clock();

		fprintf(stderr, "[bwa_aln_core] refine gapped alignments... ");
		bwa_refine_gapped(bns, n_seqs, seqs, 0);
		fprintf(stderr, "%.2f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC); t = clock();

		fprintf(stderr, "[bwa_aln_core] print alignments... ");
		for (i = 0; i < n_seqs; ++i)
		{
			if(rg_number < 0) 
			{
				bwa_print_sam1(bns,  seqs + i, 0, opt->mode, opt->max_top2);
			} else {
				char *rg_id = sl[(seqs + i)->sample_id]->rg_id;
				//fprintf(stderr, "sampleID:%d\nRG:%s\n", (seqs + i)->sample_id, rg_id);
				bwa_print_sam2(bns,  seqs + i, 0, opt->mode, opt->max_top2, rg_id);
			}
		}
		fprintf(stderr, "%.2f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC);

		bwa_free_read_seq(n_seqs, seqs);
		fprintf(stderr, "[bwa_aln_core] %d sequences have been processed.\n", tot_seqs);
		fprintf(stderr, "%.2f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC);
	}

	// destroy
	bwa_seq_close(ks);
	bns_destroy(bns);
	bwt_destroy(bwt);
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

int bwa_alnse(int argc, char *argv[])
{
	int c, opte = -1, n_occ = 3;
	gap_opt_t *opt;
	char *prefix, *rg_line = NULL, *sample_list = NULL;

	opt = gap_init_opt();
	while ((c = getopt(argc, argv, "n:o:e:i:d:l:k:LR:m:t:NM:O:E:q:f:b012IYB:hn:g:r:T:")) >= 0) {
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
		case 'T': sample_list = strdup(optarg);break;
		case 'h': break;
		case 'r':
			if ((rg_line = bwa_set_rg(optarg)) == 0) return 1;
			break;
		case 'g': n_occ = atoi(optarg); break;
		case 'f': xreopen(optarg, "w", stdout); break;
		default: return 1;
		}
	}
	if (opte > 0) {
		opt->max_gape = opte;
		opt->mode &= ~BWA_MODE_GAPE;
	}

	if (optind + 2 > argc) {
		fprintf(stderr, "\n");
		fprintf(stderr, "Usage:   bwa aln [options] <prefix> <in.fq>\n\n");
		fprintf(stderr, "Options: -n NUM    max #diff (int) or missing prob under %.2f err rate (float) [%.2f]\n",
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
		fprintf(stderr, "samse Options:\n");
		fprintf(stderr, "         -g INT   maximum occurrences for one end [%d]\n", n_occ);
		fprintf(stderr, "         -r STR   read group header line such as `@RG\\tID:foo\\tSM:bar' [null]\n");
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
		return 1;
	}
	bwa_alnse_core(prefix, argv[optind+1], opt, n_occ, rg_line, sample_list);
	free(opt); free(prefix);
	return 0;
}