/* functions from part_func.c */
/* calculate partition function and base pair probabilities */
#ifndef __VIENNA_RNA_PACKAGE_PART_FUNC_CO_H__
#define __VIENNA_RNA_PACKAGE_PART_FUNC_CO_H__

#include "data_structures.h"

extern int    mirnatog;     /*toggles no intrabp in 2nd mol*/
extern double F_monomer[2]; /* free energies of the two monomers */

cofoldF co_pf_fold(char *sequence, char *structure); /* calculate partition function and base pair probabilities */
void    init_co_pf_fold(int length);
void    free_co_pf_arrays(void);
void    update_co_pf_params(int length); /*recalculate energy parameters */
char    co_bppm_symbol(float *x);    /* string representation of structure */
void    compute_probabilities(double FAB, double FEA, double FEB,
                              struct plist  *prAB,
                              struct plist  *prA, struct plist  *prB,
                              int Alength);
ConcEnt *get_concentrations(double FEAB, double FEAA, double FEBB, double FEA, double FEB, double * startconc);
plist   *get_plist(struct plist *pl, int length, double cut_off);
int     make_probsum(int length, char *name);

#endif
