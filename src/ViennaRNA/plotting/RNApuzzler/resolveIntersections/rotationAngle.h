#ifndef RNAPUZZLER_ROTATION_ANGLE_H
#define RNAPUZZLER_ROTATION_ANGLE_H

#include "ViennaRNA/plotting/RNApuzzler/dataTypes/configtree_struct.h"
#include "ViennaRNA/plotting/RNApuzzler/resolveIntersections/intersectionType.h"

PRIVATE double getRotationAngle(const treeNode          *rootNode,
                        const treeNode          *centerNode,
                        const treeNode          *intersectorNode,
                        const intersectionType  it,
                        short                   rotationSign);


#include "rotationAngle.inc"


#endif
