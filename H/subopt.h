/* subopt.h */
#ifndef __VIENNA_RNA_PACKAGE_SUBOPT_H__
#define __VIENNA_RNA_PACKAGE_SUBOPT_H__

#include "data_structures.h"

SOLUTION *subopt (char *seq, char *sequence, int delta, FILE *fp);
SOLUTION *subopt_circ (char *seq, char *sequence, int delta, FILE *fp);
  /* returns list of subopt structures or writes to fp */

extern  int     subopt_sorted;                /* sort output by energy */
extern  int     density_of_states[MAXDOS+1];
extern  double  print_energy;                 /* printing threshold for use with logML */

#endif
/* End of file */
