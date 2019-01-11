#ifndef HANDLE_ANCESTOR_INTERSECTIONS_H
#define HANDLE_ANCESTOR_INTERSECTIONS_H

#include "ViennaRNA/RNApuzzler/dataTypes/configtree_struct.h"
#include "ViennaRNA/RNApuzzler/definitions.h"

treeNode *checkNodeAgainstAncestors(
        treeNode* node,
        puzzlerOptions* puzzler
);

#endif
