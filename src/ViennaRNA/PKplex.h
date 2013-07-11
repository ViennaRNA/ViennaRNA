#ifndef PKPLEX_H
#define PKPLEX_H

#include <ViennaRNA/data_structures.h>

extern dupVar *PlexHits;
extern int    PlexHitsArrayLength;
extern int    NumberOfHits;
extern int    verbose;


/**
 *  \brief 
 */
dupVar  **PKLduplexfold_XS( const char *s1,
                            int **access_s1,
                            const int threshold,
                            const int alignment_length,
                            const int delta);

int     arraySize(duplexT **array);

void    freeDuplexT(duplexT **array);

#endif
