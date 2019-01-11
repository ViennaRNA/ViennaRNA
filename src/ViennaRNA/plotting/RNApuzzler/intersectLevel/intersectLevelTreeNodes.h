#ifndef INTERSECTLEVELTREENODES_H
#define INTERSECTLEVELTREENODES_H

#include "ViennaRNA/RNApuzzler/definitions.h"
#include "ViennaRNA/RNApuzzler/dataTypes/configtree_struct.h"
#include "ViennaRNA/RNApuzzler/resolveIntersections/intersectionType.h"

/**
 * @brief intersectNodeExterior
 * @param node
 * @return
 */
short intersectNodeExterior(const treeNode* node, const puzzlerOptions* puzzler);

/**
 * @brief intersectNodeNode
 *      - intersects two nodes of a configtree. each node containts a loopbox and a stembox.
 *        intersection detection for all combinations is done in this method. The parameters
 *        LxL, LxS, SxL and SxS are used as return values to state which kind of intersection
 *        occurs.
 * @param node1
 * @param node2
 * @param bulge1 index of intersecting bulge of node1 if Bx*, -1 otherwise
 * @param bulge2 index of intersecting bulge of node2 if *xB, -1 otherwise
 * @return intersection type
 */
intersectionType intersectNodeNode(
        const treeNode* node1,
        const treeNode* node2
);

/**
 * @brief intersectNodeTree
 *      - intersects a single node of a configtree with whole subtree.
 * @param node
 *      - tree node ptr for node
 * @param tree
 *      - tree node ptr for root of subtree
 * @param intersectorNode
 *      - ptr to intersecting node in tree, undefined if not intersecting
 * @return
 *      - 1 if an intersection occured, 0 otherwise
 */
short intersectNodeTree(const treeNode* node, treeNode* tree, treeNode** intersectorNode);

/*
 * starting method for detection of intersections between trees
 * basically this one iterates over both subtrees and does the intersection check
 * for each pair from tree1 and tree2
 * this is done recursively...
 * iterate over tree1 and for each node we iterate over tree2 for intersection calcultion
 */
short intersectTrees(treeNode* tree1, treeNode* tree2);

short intersectNodeLists(
        treeNode* const * list1,
        const int size1,
        treeNode* const * list2,
        const int size2,
        const puzzlerOptions* const puzzler
);

#endif
