#ifndef HANDLE_ANCESTOR_INTERSECTIONS_H
#define HANDLE_ANCESTOR_INTERSECTIONS_H

#include "ViennaRNA/plotting/RNApuzzler/dataTypes/configtree_struct.h"
#include "ViennaRNA/plotting/RNApuzzler/definitions.h"

PRIVATE treeNode *checkNodeAgainstAncestors(treeNode        *node,
                                    puzzlerOptions  *puzzler);


#include "handleAncestorIntersections.c"

#endif
