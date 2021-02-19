#ifndef PKPLEX_H
#define PKPLEX_H

#include <ViennaRNA/datastructures/basic.h>

/**
 *  \brief 
 */
vrna_pkplex_t *
PKLduplexfold_XS(const char *s1,
                 const int  **access_s1,
                 int  penalty,
                 int  max_interaction_length,
                 int  delta);

int **
vrna_PKplex_accessibilities(const char    *sequence,
                            unsigned int  unpaired,
                            double        cutoff);

#endif
