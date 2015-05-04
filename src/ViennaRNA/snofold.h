/* function from fold.c */
#ifndef VIENNA_RNA_PACKAGE_SNOFOLD_H
#define VIENNA_RNA_PACKAGE_SNOFOLD_H

#include <ViennaRNA/data_structures.h>

/* Normal fold */

/**
*** snofold is the stem folding array for RNAsnoop
**/
int  snofold( const char *sequence,
		char *structure,
                const int max_assym,
                const int threshold, 
                const int min_s2,
                const int max_s2,
                const int half_stem,
                const int max_half_stem);
/**
*** Free arrays and structure related to snofold
**/

void   snofree_arrays(const int length);  /* free arrays for mfe folding */
void   snoinitialize_fold(int length);    /* allocate arrays for folding */
void   snoupdate_fold_params(void);       /* recalculate parameters */
int    snoloop_energy(short *ptable,
                      short *s,
                      short *s1,
                      int i);
void   snoexport_fold_arrays( int **indx_p,
                              int **mLoop_p,
                              int **cLoop,
                              folden ***fold_p,
                              folden ***fold_p_XS);
char * snobacktrack_fold_from_pair( const char *sequence,
                                    int i,
                                    int j);
/* alifold */
float alisnofold( const char **strings,
                  const int max_assym,
                  const int threshloop, 
                  const int min_s2,
                  const int max_s2,
                  const int half_stem,
                  const int max_half_stem);
void  alisnofree_arrays(const int length);
char  *alisnobacktrack_fold_from_pair(const char **sequence,
                                      int i,
                                      int j,
                                      int *cov);
extern double cv_fact /* =1 */;
extern double nc_fact /* =1 */;

/* max number of mismatch >>>>>..((   )).>>>> */
#define MISMATCH 3

#endif
