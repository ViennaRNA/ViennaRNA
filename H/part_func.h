#ifndef __PART_FUNC__
#define __PART_FUNC__

#define FLT_OR_DBL double
/* functions from part_func.c */
extern float  pf_fold(char *sequence, char *structure);
/* calculate partition function and base pair probabilities */
extern void   init_pf_fold(int length);    /* allocate space for pf_fold() */
extern void   free_pf_arrays(void);        /* free arrays from pf_fold() */
extern void   update_pf_params(int length); /*recalculate energy parameters */
extern char   bppm_symbol(const float *x);  /* string representation of structure */
extern double mean_bp_dist(int length); /* mean pair distance of ensemble */
extern char  *centroid(int length, double *dist);     /* mean pair distance of ensemble */
extern int get_pf_arrays(short **S_p, short **S1_p, char **ptype_p, FLT_OR_DBL **qb_p, FLT_OR_DBL **qm_p, FLT_OR_DBL **q1k_p, FLT_OR_DBL **qln_p);

extern	char	*pbacktrack(char *sequence);
extern	float	pf_circ_fold(char *sequence, char *structure);
extern	char	*pbacktrack_circ(char *seq);

#endif
