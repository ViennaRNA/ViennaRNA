#ifndef RNAPUZZLER_HANDLE_ANCESTOR_INTERSECTIONS_H
#define RNAPUZZLER_HANDLE_ANCESTOR_INTERSECTIONS_H

#include "ViennaRNA/plotting/RNApuzzler/dataTypes/configtree_struct.h"
#include "ViennaRNA/plotting/RNApuzzler/definitions.h"

PRIVATE treeNode *
checkNodeAgainstAncestors(treeNode        *node,
                          puzzlerOptions  *puzzler);


#include "handleAncestorIntersections.inc"

#endif
