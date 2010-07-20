#ifndef __VIENNA_RNA_PACKAGE_COFOLD_H__
#define __VIENNA_RNA_PACKAGE_COFOLD_H__

#include "data_structures.h"

#ifdef __GNUC__
#define DEPRECATED(func) func __attribute__ ((deprecated))
#else
#define DEPRECATED(func) func
#endif

/**
*** \file cofold.h
**/

float     cofold(const char *sequence, char *structure); 
void      free_co_arrays(void);          /* free arrays for mfe folding */
void      initialize_cofold(int length); /* allocate arrays for folding */
void      update_cofold_params(void);    /* recalculate parameters */
SOLUTION  *zukersubopt(const char *string);
float     *get_monomer_mfes();

/**
*** allocate arrays for folding
*** \deprecated{This function is obsolete and will be removed soon!}
**/
void      DEPRECATED(initialize_cofold(int length));
#endif
