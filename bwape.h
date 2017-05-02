#ifndef BWAPE_H
#define BWAPE_H

#include <unistd.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <stdio.h>
#include <string.h>
#include "bwtaln.h"
#include "kvec.h"
#include "bntseq.h"
#include "utils.h"
#include "bwase.h"
#include "bwa.h"
#include "ksw.h"

#ifdef __cplusplus
extern "C" {
#endif
	#ifdef USE_MALLOC_WRAPPERS
#  include "malloc_wrap.h"
#endif
	typedef struct {
		int n;
		bwtint_t *a;
	} poslist_t;

	typedef struct {
		double avg, std, ap_prior;
		bwtint_t low, high, high_bayesian;
	} isize_info_t;

	typedef struct {
		pair64_v arr;
		pair64_v pos[2];
		kvec_t(bwt_aln1_t) aln[2];
	} pe_data_t;
	
	typedef struct {
		kvec_t(bwt_aln1_t) aln;
	} aln_buf_t;

	void bwa_aln2seq_core(int n_aln, const bwt_aln1_t *aln, bwa_seq_t *s, int set_main, int n_multi);
	void bwa_aln2seq(int n_aln, const bwt_aln1_t *aln, bwa_seq_t *s);
	int bwa_approx_mapQ(const bwa_seq_t *p, int mm);
	void bwa_print_sam1(const bntseq_t *bns, bwa_seq_t *p, const bwa_seq_t *mate, int mode, int max_top2);
	bntseq_t *bwa_open_nt(const char *prefix);
	void bwa_print_sam_SQ(const bntseq_t *bns);
	
	pe_opt_t *bwa_init_pe_opt();
	int bwa_cal_pac_pos_pe(const bntseq_t *bns, const char *prefix, bwt_t *const _bwt, int n_seqs, bwa_seq_t *seqs[2], FILE *fp_sa[2], isize_info_t *ii,
					   const pe_opt_t *opt, const gap_opt_t *gopt, const isize_info_t *last_ii);
	int bwa_cal_pac_pos_pe2(const bntseq_t *bns, const char *prefix, bwt_t *const _bwt, int n_seqs, bwa_seq_t *seqs[2], isize_info_t *ii,
					   const pe_opt_t *opt, const gap_opt_t *gopt, const isize_info_t *last_ii);
	ubyte_t *bwa_paired_sw(const bntseq_t *bns, const ubyte_t *_pacseq, int n_seqs, bwa_seq_t *seqs[2], const pe_opt_t *popt, const isize_info_t *ii);
	
#ifdef __cplusplus
}
#endif

#endif // BWAPE_H