#ifndef __VIENNA_RNA_PACKAGE_LPFOLD_H__
#define __VIENNA_RNA_PACKAGE_LPFOLD_H__

#include "data_structures.h"

#ifdef __GNUC__
#define DEPRECATED(func) func __attribute__ ((deprecated))
#else
#define DEPRECATED(func) func
#endif

void    update_pf_paramsLP(int length);
struct  plist *pfl_fold(char *sequence, int winSize, int pairSize, float cutoffb, double **pU, struct plist **dpp2, FILE *pUfp, FILE *spup);
void    putoutpU_prob(double **pU,int length, int ulength, FILE *fp, int energies);

/**
*** Dunno if this function was ever used by external programs linking to RNAlib, but it
*** was declared PUBLIC before.
*** Anyway, never use this function as it will be removed soon and does nothing at all
**/
DEPRECATED(void init_pf_foldLP(int length));

#endif
