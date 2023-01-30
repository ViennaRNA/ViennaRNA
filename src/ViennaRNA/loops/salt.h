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
vrna_salt_loop(int m, double salt, double T, double standard);

int
vrna_salt_stack(double salt, double T, double standard);

void
vrna_salt_ml(double saltLoop[], int lower, int upper, int *m, int *b);

int
vrna_salt_duplex_init(double salt, double standard);

#endif
