#ifndef OPTIMIZE_H
#define OPTIMIZE_H

#include "ViennaRNA/RNApuzzler/dataTypes/configtree_struct.h"
#include "ViennaRNA/RNApuzzler/definitions.h"

typedef enum {
  INCREASE_ALL_OTHER,
  INCREASE_LEFT_NEIGHBOR,
  INCREASE_RIGHT_NEIGHBOR,
  INCREASE_BOTH_NEIGHBORS
} increaseStrategy;

typedef enum {
  BINARY_SEARCH,
  LINEAR_SEARCH
} searchStrategy;

typedef enum {
  DISTRIBUTE_EQUALLY,
  DISTRIBUTE_PROPORTIONALLY
} distributionStrategy;

/**
 * @brief optimizeTree
 *      Performs an optimization operation to all tree nodes
 * @param root the root node of your tree
 * @return product of all partial radius ratios (0 < ratio <= 1)
 */
double optimizeTree(
    treeNode* node,
    puzzlerOptions* puzzler
);

#endif
