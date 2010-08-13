/*
      minimum free energy
      RNA secondary structure with
      basepair distance d to reference structure prediction
      
*/
#ifndef __VIENNA_RNA_PACKAGE_TWO_D_FOLD_H__
#define __VIENNA_RNA_PACKAGE_TWO_D_FOLD_H__

#include "data_structures.h"

#ifdef __GNUC__
#define DEPRECATED(func) func __attribute__ ((deprecated))
#else
#define DEPRECATED(func) func
#endif

TwoDfold_vars     *get_TwoDfold_variables(const char *seq, const char *structure1, const char *structure2, int circ);
void              destroy_TwoDfold_variables(TwoDfold_vars *our_variables);

TwoDfold_solution **TwoDfold(const char *string, char *structure1, char *structure2);
TwoDfold_solution **TwoDfold_bound(TwoDfold_vars *our_variables, int distance1, int distance2);

TwoDfold_solution **TwoDfold_circ(const char *string, char *structure1, char *structure2);
TwoDfold_solution **TwoDfold_circ_bound(TwoDfold_vars *our_variables, int maxDistance1, int maxDistance2);

/* function for Pathfinder.c */
char *TwoDfold_backtrack_f5(unsigned int j, unsigned int k, unsigned int l, TwoDfold_vars *vars);

/* deprecated */
DEPRECATED(TwoDfold_solution **TwoDfold_bounded(const char *string, char *structure1, char *structure2, int maxDistance1, int maxDistance2));
/* deprecated */
DEPRECATED(TwoDfold_solution **TwoDfold_circ_bounded(const char *string, char *structure1, char *structure2, int maxDistance1, int maxDistance2));

#endif
