#ifndef RESOLVE_INTERSECTIONS_H
#define RESOLVE_INTERSECTIONS_H

#include "ViennaRNA/RNApuzzler/definitions.h"
#include "ViennaRNA/RNApuzzler/dataTypes/configtree_struct.h"

treeNode *checkAndFixIntersections(
        treeNode* node,
        const int recursionDepth,
        puzzlerOptions* puzzler
);

#endif
