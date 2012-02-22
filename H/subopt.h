/* subopt.h */
#ifndef __VIENNA_RNA_PACKAGE_SUBOPT_H__
#define __VIENNA_RNA_PACKAGE_SUBOPT_H__

#include "data_structures.h"

#define MAXDOS 1000

/**
 *  \file subopt.h
 *  \brief RNAsubopt and density of states declarations
 */

/**
 *  \brief Returns list of subopt structures or writes to fp
 * 
 *  This function produces <b>all</b> suboptimal secondary structures within
 *  'delta' * 0.01 kcal/mol of the optimum. The results are either
 *  directly written to a 'fp' (if 'fp' is not NULL), or
 *  (fp==NULL) returned in a #SOLUTION * list terminated
 *  by an entry were the 'structure' pointer is NULL.
 * 
 *  \param  seq
 *  \param  sequence
 *  \param  delta
 *  \param  fp
 *  \return
 */
SOLUTION *subopt (char *seq,
                  char *sequence,
                  int delta,
                  FILE *fp);

SOLUTION *subopt_par( char *seq,
                      char *structure,
                      paramT *parameters,
                      int delta,
                      int is_constrained,
                      int is_circular,
                      FILE *fp);

/**
 *  \brief Returns list of circular subopt structures or writes to fp
 * 
 *  This function is similar to subopt() but calculates secondary structures
 *  assuming the RNA sequence to be circular instead of linear
 * 
 *  \param  seq
 *  \param  sequence
 *  \param  delta
 *  \param  fp
 *  \return
 */
SOLUTION *subopt_circ ( char *seq,
                        char *sequence,
                        int delta,
                        FILE *fp);

/**
 *  \brief Sort output by energy
 */
extern  int     subopt_sorted;

extern  int     density_of_states[MAXDOS+1];

/**
 *  \brief printing threshold for use with logML
 */
extern  double  print_energy;

#endif
/* End of file */
