/* functions from part_func.c */
extern float  pf_fold(char *sequence, char *structure);
/* calculate partition function and base pair probabilities */
extern void   init_pf_fold(int length);    /* allocate space for pf_fold() */
extern void   free_pf_arrays(void);        /* free arrays from pf_fold() */
extern void   update_pf_params(int length); /*recalculate energy parameters */
extern char   bppm_symbol(float *x);    /* string representation of structure */
