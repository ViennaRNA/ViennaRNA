#ifndef RNA_TURTLE_H
#define RNA_TURTLE_H

#include <ViennaRNA/plotting/RNApuzzler/dataTypes/tBaseInformation_struct.h>

/**
 Compute affin angles for all loops
 */
void computeAffineCoordinates(
    short const * const pair_table,
    const double paired,
    const double unpaired,
    tBaseInformation* const baseInformation
);

/**
 Calculate the coordinates for the drawing with the given affin angles
 */
void affineToCartesianCoordinates(
    tBaseInformation* const baseInformation,
    unsigned short const length,
    double * const x,
    double * const y
);


/**
 * Compute layout using RNAturtle algorithm
 */
int layout_RNAturtle(
    short const * const pair_table,
    float *x,
    float *y,
    double *arc_coords
);

#endif
