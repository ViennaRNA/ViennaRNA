#ifndef RNATURTLE_H
#define RNATURTLE_H

#include <ViennaRNA/plotting/RNApuzzler/dataTypes/tBaseInformation_struct.h>

/**
 * Compute layout using RNAturtle algorithm
 */
int
layout_RNAturtle(short const *const pair_table,
                 float              *x,
                 float              *y,
                 double             *arc_coords);


#endif
