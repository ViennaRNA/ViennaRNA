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

/**
*** toggles no intrabp in 2nd mol
**/
extern int    mirnatog;

/**
*** free energies of the two monomers
**/
extern double F_monomer[2];

/**
*** calculate partition function and base pair probabilities
**/
cofoldF co_pf_fold(char *sequence, char *structure);

void    free_co_pf_arrays(void);

/**
*** recalculate energy parameters
**/
void    update_co_pf_params(int length);

void    compute_probabilities(double FAB, double FEA, double FEB,
                              struct plist  *prAB,
                              struct plist  *prA, struct plist  *prB,
                              int Alength);
ConcEnt *get_concentrations(double FEAB, double FEAA, double FEBB, double FEA, double FEB, double * startconc);

/*
#################################################
# DEPRECATED FUNCTIONS                          #
#################################################
*/

/**
*** DO NOT USE THIS FUNCTION ANYMORE
*** \deprecated{ This function is deprecated and will be removed soon!}
*** use \ref assign_plist_from_pr() instead!
**/
DEPRECATED(plist  *get_plist(struct plist *pl, int length, double cut_off));
/**
*** DO NOT USE THIS FUNCTION ANYMORE
*** \deprecated{ This function is deprecated and will be removed soon!}
**/
DEPRECATED(void   init_co_pf_fold(int length));

#endif
