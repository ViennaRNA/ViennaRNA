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
/**
*** free arrays for mfe folding
**/
void      free_co_arrays(void);
/**
*** allocate arrays for folding
**/
void      initialize_cofold(int length);
/**
*** recalculate parameters
**/
void      update_cofold_params(void);
/**
*** Compute Suboptimal structures according Zuker
**/
SOLUTION  *zukersubopt(const char *string);

void get_monomere_mfes(float *e1, float *e2);

void export_cofold_arrays(int **f5_p,
                          int **c_p,
                          int **fML_p,
                          int **fM1_p,
                          int **fc_p,
                          int **indx_p,
                          char **ptype_p);

/**
*** allocate arrays for folding
*** \deprecated{This function is obsolete and will be removed soon!}
**/
DEPRECATED(void initialize_cofold(int length));
#endif
