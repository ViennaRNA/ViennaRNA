#ifndef __VIENNA_RNA_PACKAGE_PART_FUNC_H__
#define __VIENNA_RNA_PACKAGE_PART_FUNC_H__

#define FLT_OR_DBL double
/* functions from part_func.c */
float   pf_fold(char *sequence, char *structure);
/* calculate partition function and base pair probabilities */
void    init_pf_fold(int length);    /* allocate space for pf_fold() */
void    free_pf_arrays(void);        /* free arrays from pf_fold() */
void    update_pf_params(int length); /*recalculate energy parameters */
char    bppm_symbol(float *x);  /* string representation of structure */
double  mean_bp_dist(int length); /* mean pair distance of ensemble */
char    *centroid(int length, double *dist);     /* mean pair distance of ensemble */
int     get_pf_arrays(short **S_p, short **S1_p, char **ptype_p, FLT_OR_DBL **qb_p, FLT_OR_DBL **qm_p, FLT_OR_DBL **q1k_p, FLT_OR_DBL **qln_p);

char    *pbacktrack(char *sequence);
float   pf_circ_fold(char *sequence, char *structure);
char    *pbacktrack_circ(char *seq);

extern  int st_back;

/* deprecated, use exp_E_IntLoop() from loop_energies.h instead */
double  expLoopEnergy(int u1, int u2, int type, int type2, short si1, short sj1, short sp1, short sq1) __attribute__ ((deprecated));
/* deprecated, use exp_E_Hairpin() from loop_energies.h instead */
double  expHairpinEnergy(int u, int type, short si1, short sj1, const char *string) __attribute__ ((deprecated));

#endif
