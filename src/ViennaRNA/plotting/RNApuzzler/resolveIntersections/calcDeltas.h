#ifndef CALC_DELTAS_H
#define CALC_DELTAS_H

#include "ViennaRNA/RNApuzzler/definitions.h"
#include "ViennaRNA/RNApuzzler/dataTypes/configtree_struct.h"

/**
 * @brief calcDeltas
 *      The area between stems indexLeft and indexRight
 *      (by traversing the loop clockwise starting at indexLeft-stem)
 *      will be enlarged in degree as given via deltaAngle.
 *      All other areas will be used to compensate that increase
 *      (i.e. by decreasing those area's angles).
 * @param node
 * @param recursiveEnd
 * @param indexLeft
 * @param indexRight
 * @param deltaAngle
 * @param deltas
 * @return the amount of change (in positive degree) that can be accomplished with calculated deltas
 */
double calcDeltas(
        const treeNode* node,
        const treeNode* recursiveEnd,
        const int indexLeft,
        const int indexRight,
        const double deltaAngle,
        puzzlerOptions* puzzler,
        double* deltas
);


#endif
