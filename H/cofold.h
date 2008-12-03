/* function from fold.c */
#include "subopt.h"
extern float  cofold(char *sequence, char *structure); 
/* calculate energy of string on structure */
extern void   free_co_arrays(void);          /* free arrays for mfe folding */
extern void   initialize_cofold(int length); /* allocate arrays for folding */
extern void   update_cofold_params(void);    /* recalculate parameters */
extern SOLUTION *zukersubopt(const char *string);
