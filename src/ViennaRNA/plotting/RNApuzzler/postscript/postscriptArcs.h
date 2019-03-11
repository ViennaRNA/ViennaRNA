#ifndef RNAPUZZLER_POSTSCRIPT_ARCS
#define RNAPUZZLER_POSTSCRIPT_ARCS

#include "ViennaRNA/plotting/RNApuzzler/dataTypes/tBaseInformation_struct.h"

/**
 * Compute arcs instead of lines for postscript loops.
 */
PRIVATE void computeAnglesAndCentersForPS(short const *const      pair_table,
                                  double *const           x,
                                  double *const           y,
                                  const tBaseInformation  *baseInformation,
                                  double                  *arcCoords);

#include "postscriptArcs.c"

#endif
