#ifndef __VIENNA_RNA_PACKAGE_PART_FUNC_UP_H__
#define __VIENNA_RNA_PACKAGE_PART_FUNC_UP_H__

#include "data_structures.h"

#define 		RNA_UP_MODE_1 	1U
#define 		RNA_UP_MODE_2 	2U
#define 		RNA_UP_MODE_3 	4U

/**
*** \file part_func_up.h
*** \brief Partition Function Cofolding as stepwise process
***
*** In this approach to cofolding the interaction between two RNA molecules is
*** seen as a stepwise process. In a first step, the target molecule has to
*** adopt a structure in which a binding site is accessible. In a second step,
*** the ligand molecule will hybridize with a region accessible to an
*** interaction. Consequently the algorithm is designed as a two step process:
*** The first step is the calculation of the probability
*** that a region within the target is unpaired, or equivalently, the
*** calculation of the free energy needed to expose a region. In the second step
*** we compute the free energy of an interaction for every possible binding site.
**/

/**
*** \brief
**/
pu_contrib *pf_unstru(char *sequence, int max_w);

/**
*** \brief
**/
interact *pf_interact(const char *s1, const char *s2, pu_contrib *p_c, pu_contrib *p_c2, int max_w, char *cstruc, int incr3, int incr5);

/**
*** \brief
**/
void free_interact(interact *pin);

/**
*** \brief
**/
int         Up_plot(pu_contrib *p_c, pu_contrib *p_c_sh, interact *pint, char *ofile, int **unpaired_values, char *select_contrib, char *head, unsigned int mode);

/**
*** \brief
**/
pu_contrib  *get_pu_contrib_struct(unsigned int n, unsigned int w);

/**
*** \brief
**/
void        free_pu_contrib_struct(pu_contrib *pu);

#endif
