#ifndef RNA_PUZZLER_H
#define RNA_PUZZLER_H

#include "ViennaRNA/plotting/RNApuzzler/definitions.h"

/**
 * Compute layout using RNApuzzler algorithm
 */
int layout_RNApuzzler(short const *const  pair_table,
                      float               *x,
                      float               *y,
                      double              *arc_coords,
                      puzzlerOptions      *puzzler);


#endif
