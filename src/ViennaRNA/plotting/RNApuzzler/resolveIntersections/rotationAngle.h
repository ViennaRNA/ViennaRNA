#ifndef ROTATION_ANGLE_H
#define ROTATION_ANGLE_H

#include "ViennaRNA/RNApuzzler/dataTypes/configtree_struct.h"
#include "ViennaRNA/RNApuzzler/resolveIntersections/intersectionType.h"

double getRotationAngle(
        const treeNode* rootNode,
        const treeNode* centerNode,
        const treeNode* intersectorNode,
        const intersectionType it,
        short rotationSign
);

#endif
