#ifndef RNAPUZZLER_CONFIGTREE_STRUCT_H
#define RNAPUZZLER_CONFIGTREE_STRUCT_H

#include "AABB_struct.h"
#include "config_struct.h"

/**
 * @brief The configtree struct.
 * This struct keeps your configurations (config) in a tree structure.
 * The nodes of that tree are configurations for every loop in the given RNA.
 * The main purpose of this struct is its usage for intersection detection.
 * The whole RNA has a tree-like structure.
 * The intention for intersection detection is to clean up intersections
 * beginning at small subtrees (e.g. leafs) and going up in that tree
 * level by level until the root is reached and any intersection is gone.
 */
typedef struct configtree {
  // node name
  int                     id;

  // tree information
  struct configtree       *parent;
  struct configtree       **children;
  int                     childCount;

  // RNA data
  config                  *cfg;
  int                     loop_start;
  int                     stem_start;

  // for intersection handling
  struct boundingboxLoop  *lBox;  // bounding box for this loop         first base at loop_start
  struct boundingboxStem  *sBox;  // bounding box for the prior stem    first base at stem_start

  // AABB
  AABB                    aabb;
} treeNode;

#endif
