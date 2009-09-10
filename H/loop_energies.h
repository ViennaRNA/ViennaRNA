#ifndef __VIENNA_RNA_PACKAGE_LOOP_ENERGIES_H__
#define __VIENNA_RNA_PACKAGE_LOOP_ENERGIES_H__

#include "params.h"

int   E_IntLoop(int n1, int n2, int type, int type_2, int si1, int sj1, int sp1, int sq1, paramT *P);
int   E_Hairpin(int size, int type, int si1, int sj1, const char *string, paramT *P);

/* compute Boltzmann weight of a hairpin loop, multiply by scale[u+2] */
double  exp_E_Hairpin(int u, int type, short si1, short sj1, const char *string, pf_paramT *P);
/* compute Boltzmann weight of interior loop, multiply by scale[u1+u2+2] for scaling */
double  exp_E_IntLoop(int u1, int u2, int type, int type2, short si1, short sj1, short sp1, short sq1, pf_paramT *P);

#endif
