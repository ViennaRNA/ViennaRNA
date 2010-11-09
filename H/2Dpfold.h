/*
      minimum free energy
      RNA secondary structure with
      basepair distance d to reference structure prediction
      
*/
#ifndef __VIENNA_RNA_PACKAGE_TWO_D_PF_FOLD_H__
#define __VIENNA_RNA_PACKAGE_TWO_D_PF_FOLD_H__

#include "data_structures.h"

#ifdef __GNUC__
#define DEPRECATED(func) func __attribute__ ((deprecated))
#else
#define DEPRECATED(func) func
#endif

TwoDpfold_vars  *get_TwoDpfold_variables(const char *seq, const char *structure1, char *structure2, int circ);
TwoDpfold_vars  *get_TwoDpfold_variables_from_MFE(TwoDfold_vars *mfe_vars);
void            destroy_TwoDpfold_variables(TwoDpfold_vars *our_variables);

FLT_OR_DBL      **TwoDpfold(TwoDpfold_vars *our_variables, int maxDistance1, int maxDistance2);
FLT_OR_DBL      **TwoDpfold_circ(TwoDpfold_vars *our_variables, int maxDistance1, int maxDistance2);

char            *TwoDpfold_pbacktrack(TwoDpfold_vars *vars, unsigned int d1, unsigned int d2);
char            *TwoDpfold_pbacktrack_f5(TwoDpfold_vars *vars, unsigned int d1, unsigned int d2, unsigned int length);

#endif
