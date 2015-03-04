#ifndef VIENNA_RNA_PACKAGE_PLEX_H
#define VIENNA_RNA_PACKAGE_PLEX_H

#include <ViennaRNA/data_structures.h>


extern int subopt_sorted;

/**
*** Lduplexfold Computes duplexes between two single sequences
**/
duplexT** Lduplexfold(const char *s1,
                      const char *s2,
                      const int threshold,
                      const int extension_cost,
                      const int alignment_length,
                      const int delta,
                      const int fast,
                      const int il_a,
                      const int il_b,
                      const int b_a,
                      const int b_b);

/**
*** Lduplexfold_XS Computes duplexes between two single sequences with accessibility
**/
duplexT** Lduplexfold_XS( const char*s1,
                          const char* s2,
                          const int **access_s1,
                          const int **access_s2,
                          const int threshold,
                          const int delta,
                          const int alignment_length,
                          const int fast,
                          const int il_a,
                          const int il_b,
                          const int b_a,
                          const int b_b);/* , const int target_dead, const int query_dead); */

/**
*** Lduplexfold_C Computes duplexes between two single sequences and takes constraint into account
**/
duplexT** Lduplexfold_C(const char *s1,
                        const char *s2,
                        const int threshold,
                        const int extension_cost,
                        const int alignment_length,
                        const int delta,
                        const int fast,
                        const char* structure,
                        const int il_a,
                        const int il_b,
                        const int b_a,
                        const int b_b);

/**
*** Lduplexfold_CXS Computes duplexes between two single sequences and takes constraint as well as accessibility into account
**/

duplexT** Lduplexfold_CXS(const char*s1,
                          const char* s2,
                          const int **access_s1,
                          const int **access_s2,
                          const int threshold,
                          const int delta,
                          const int alignment_length,
                          const int fast,
                          const char* structure,
                          const int il_a,
                          const int il_b,
                          const int b_a,
                          const int b_b); /*, const int target_dead, const int query_dead); */




int      arraySize(duplexT** array);
void     freeDuplexT(duplexT** array);

#endif
