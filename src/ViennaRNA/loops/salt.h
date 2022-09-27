#ifndef VIENNA_RNA_PACKAGE_LOOPS_SALT_H
#define VIENNA_RNA_PACKAGE_LOOPS_SALT_H


#include <math.h>
#include "ViennaRNA/utils/basic.h"

#ifdef __GNUC__
# define INLINE inline
#else
# define INLINE
#endif


double
vrna_salt_loop(int m, double salt, double T);

int
vrna_salt_stack(double salt, double T);

void
vrna_salt_ml(double saltLoop[], int lower, int upper, int *m, int *b);

#endif
