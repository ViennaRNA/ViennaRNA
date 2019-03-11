#ifndef ROTATION_ANGLE_H
#define ROTATION_ANGLE_H

#include "ViennaRNA/plotting/RNApuzzler/dataTypes/configtree_struct.h"
#include "ViennaRNA/plotting/RNApuzzler/resolveIntersections/intersectionType.h"

double getRotationAngle(const treeNode          *rootNode,
                        const treeNode          *centerNode,
                        const treeNode          *intersectorNode,
                        const intersectionType  it,
                        short                   rotationSign);


#endif
