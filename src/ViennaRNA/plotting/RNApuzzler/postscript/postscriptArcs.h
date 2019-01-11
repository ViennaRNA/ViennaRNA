#include "ViennaRNA/RNApuzzler/dataTypes/tBaseInformation_struct.h"

/**
 * Compute arcs instead of lines for postscript loops.
 */
void computeAnglesAndCentersForPS(
    short const * const pair_table,
    double * const x,
    double * const y,
    const tBaseInformation* baseInformation,
    double *arcCoords
);

