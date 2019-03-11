#ifndef RESOLVE_INTERSECTIONS_H
#define RESOLVE_INTERSECTIONS_H

#include "ViennaRNA/plotting/RNApuzzler/definitions.h"
#include "ViennaRNA/plotting/RNApuzzler/dataTypes/configtree_struct.h"

PRIVATE treeNode *checkAndFixIntersections(treeNode       *node,
                                   const int      recursionDepth,
                                   puzzlerOptions *puzzler);


#include "resolveIntersections.c"

#endif
