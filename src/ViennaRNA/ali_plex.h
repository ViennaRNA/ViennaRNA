#ifndef VIENNA_RNA_PACKAGE_ALI_PLEX_H
#define VIENNA_RNA_PACKAGE_ALI_PLEX_H

#include <ViennaRNA/data_structures.h>
/**
*** aliLduplexfold computes the duplexes between two alignments
**/
duplexT** aliLduplexfold( const char *s1[],
                          const char *s2[],
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
*** aliLduplexfold computes the duplexes between two alignments. It also takes the average accessibility into account
**/
duplexT** aliLduplexfold_XS(const char* s1[],
                            const char* s2[],
                            const int **access_s1,
                            const int **access_s2, 
                            const int threshold,
                            const int alignment_length,
                            const int delta,
                            const int fast,
                            const int il_a,
                            const int il_b,
                            const int b_a,
                            const int b_b);

/*
extern duplexT aliduplexfold(const char *s1[], const char *s2[], const int extension_cost);
extern duplexT aliduplexfold_XS(const char *s1[], const char *s2[],const int **access_s1, 
const int **access_s2, const int i_pos, const int j_pos, const int threshold);
*/
#endif
