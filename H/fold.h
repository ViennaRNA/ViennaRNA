/* function from fold.c */
extern float  fold(const char *sequence, char *structure); 
/* calculate mfe-structure of sequence */
extern float  energy_of_struct(const char *string, const char *structure);
/* calculate energy of string on structure */
extern void   free_arrays(void);           /* free arrays for mfe folding */
extern void   initialize_fold(int length); /* allocate arrays for folding */
extern void   update_fold_params(void);    /* recalculate parameters */
