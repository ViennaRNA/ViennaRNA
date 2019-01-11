#ifndef CONFIGTREE_DEBUG_H
#define CONFIGTREE_DEBUG_H

#include "ViennaRNA/RNApuzzler/dataTypes/configtree_struct.h"
#include "ViennaRNA/RNApuzzler/definitions.h"

typedef enum {
    RGB_BLACK = 0,
    RGB_GREY = 1,
    RGB_RED = 2,
    RGB_BLUE = 3,
    RGB_GREEN = 4,
    RGB_ORANGE = 5,
    RGB_PURPLE = 6,
} RGB_SELECTOR;

void printBoxesDebug(treeNode* tree);

void printBoxesDebugToFile(treeNode* tree, const puzzlerOptions* puzzler, const int changeNumber);

void PS_printTree(
        treeNode* root,
        const puzzlerOptions* puzzler
);

void PS_printPath(
        treeNode* ancestor,
        treeNode* intersector,
        const puzzlerOptions* puzzler
);

void PS_printFancyTree(
        treeNode* node,
        puzzlerOptions* puzzler
);

void PS_printFancyPath(
        treeNode* start,
        treeNode* end,
        treeNode* feature,
        puzzlerOptions* puzzler
);

void PS_printFancySiblings(
        treeNode *node,
        treeNode *left,
        treeNode *right,
        puzzlerOptions* puzzler
);

#endif
