#ifndef __VIENNA_RNA_PACKAGE_COFOLD_H__
#define __VIENNA_RNA_PACKAGE_COFOLD_H__

/* function from fold.c */
#include "subopt.h"
float     cofold(const char *sequence, char *structure); 
void      free_co_arrays(void);          /* free arrays for mfe folding */
void      initialize_cofold(int length); /* allocate arrays for folding */
void      update_cofold_params(void);    /* recalculate parameters */
SOLUTION  *zukersubopt(const char *string);
float     *get_monomer_mfes();


#endif
