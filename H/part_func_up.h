#ifndef __VIENNA_RNA_PACKAGE_PART_FUNC_UP_H__
#define __VIENNA_RNA_PACKAGE_PART_FUNC_UP_H__

#include "data_structures.h"

#define 		RNA_UP_MODE_1 	1U
#define 		RNA_UP_MODE_2 	2U
#define 		RNA_UP_MODE_3 	4U

/* prob of unpaired region of length w */
pu_contrib *pf_unstru(char *sequence, int max_w);
/* prob. of intermolecular interaction between two sequences of maximal length w*  prob unpaired from pf_unpaired */
interact *pf_interact(const char *s1, const char *s2, pu_contrib *p_c, pu_contrib *p_c2, int max_w, char *cstruc, int incr3, int incr5);
void free_interact(interact *pin);

/**
***
**/
int         Up_plot(pu_contrib *p_c, pu_contrib *p_c_sh, interact *pint, char *ofile, int **unpaired_values, char *select_contrib, char *head, unsigned int mode);

pu_contrib  *get_pu_contrib_struct(unsigned int n, unsigned int w);
void        free_pu_contrib_struct(pu_contrib *pu);

#endif
