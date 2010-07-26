/* functions from part_func.c */
/* calculate partition function and base pair probabilities */
#ifndef __VIENNA_RNA_PACKAGE_PART_FUNC_CO_H__
#define __VIENNA_RNA_PACKAGE_PART_FUNC_CO_H__

#include "data_structures.h"

#ifdef __GNUC__
#define DEPRECATED(func) func __attribute__ ((deprecated))
#else
#define DEPRECATED(func) func
#endif

/**
*** \file part_func_co.h
**/

extern int    mirnatog;     /*toggles no intrabp in 2nd mol*/
extern double F_monomer[2]; /* free energies of the two monomers */

cofoldF co_pf_fold(char *sequence, char *structure); /* calculate partition function and base pair probabilities */
void    free_co_pf_arrays(void);
void    update_co_pf_params(int length); /*recalculate energy parameters */
void    compute_probabilities(double FAB, double FEA, double FEB,
                              struct plist  *prAB,
                              struct plist  *prA, struct plist  *prB,
                              int Alength);
ConcEnt *get_concentrations(double FEAB, double FEAA, double FEBB, double FEA, double FEB, double * startconc);
int     make_probsum(int length, char *name);

/**
*** DO NOT USE THIS FUNCTION ANYMORE
*** \deprecated{ This function is deprecated and will be removed soon!}
*** use \ref assign_plist_from_pr() instead!
**/
DEPRECATED(plist  *get_plist(struct plist *pl, int length, double cut_off));
DEPRECATED(void   init_co_pf_fold(int length));

#endif
