#ifndef BOUNDING_WEDGE_H
#define BOUNDING_WEDGE_H

#include "ViennaRNA/RNApuzzler/dataTypes/configtree_struct.h"

/**
 * @brief getBoundingWedge
 *      calculates the bounding wedge for root's i-th child
 * @param root
 *      center node of the bounding wedge
 * @param childIndex
 *      index of root's child for which the bounding wedge is calculated
 * @param minAngle
 *      pointer used for returning the minimum angle
 * @param maxAngle
 *      pointer used for returning the maximum angle
 */
void getBoundingWedge(
    const treeNode* root,
    int childIndex,
    double* minAngle,
    double* maxAngle
);

#endif
