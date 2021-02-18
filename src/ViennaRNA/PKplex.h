#ifndef PKPLEX_H
#define PKPLEX_H

#include <ViennaRNA/datastructures/basic.h>

extern int    verbose;


/**
 *  \brief 
 */
vrna_pkplex_t *
PKLduplexfold_XS(const char *s1,
                 const int  **access_s1,
                 int  penalty,
                 int  max_interaction_length,
                 int  delta);

int
arraySize(duplexT **array);

void
freeDuplexT(duplexT **array);

#endif
