#ifndef RNAPUZZLER_INTERSECTLEVELTREENODES_H
#define RNAPUZZLER_INTERSECTLEVELTREENODES_H

/*
 *      RNApuzzler intersect tree nodes
 *
 *      c  Daniel Wiegreffe, Daniel Alexander, Dirk Zeckzer
 *      ViennaRNA package
 */

#include <stdio.h>

#include "ViennaRNA/utils/basic.h"

#include "definitions.inc"
#include "drawingconfig.inc"
#include "configtree.inc"
#include "boundingBoxes.inc"
#include "intersectLevelBoundingBoxes.inc"
#include "intersectionType.inc"
#include "../headers/configtree_struct.h"

/**
 * @brief intersectNodeExterior
 * @param node
 * @return
 */
PRIVATE short
intersectNodeExterior(const treeNode                    *node,
                      const vrna_plot_options_puzzler_t *puzzler);


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
PRIVATE intersectionType
intersectNodeNode(const treeNode  *node1,
                  const treeNode  *node2);


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
PRIVATE short
intersectNodeTree(const treeNode  *node,
                  treeNode        *tree,
                  treeNode        **intersectorNode);


/*
 * starting method for detection of intersections between trees
 * basically this one iterates over both subtrees and does the intersection check
 * for each pair from tree1 and tree2
 * this is done recursively...
 * iterate over tree1 and for each node we iterate over tree2 for intersection calcultion
 */
PRIVATE short
intersectTrees(treeNode *tree1,
               treeNode *tree2);


PRIVATE short
intersectNodeLists(treeNode *const                          *list1,
                   const int                                size1,
                   treeNode *const                          *list2,
                   const int                                size2,
                   const vrna_plot_options_puzzler_t *const puzzler);


PRIVATE short
intersectNodeExterior(const treeNode                    *node,
                      const vrna_plot_options_puzzler_t *puzzler)
{
  if (isExterior(node))

    return 0;

  if (isExterior(getParent(node)))

    return 0;

  double  cy  = node->lBox->c[1];
  double  r   = node->lBox->r + epsilonRecognize;

  if (puzzler->checkExteriorIntersections)
    return (cy - r) <= EXTERIOR_Y;
  else
    return 0;
}


PRIVATE short
checkBounds(const double  l1,
            const double  l2,
            const double  l3,
            const double  l4,
            const double  l5,
            const double  h1,
            const double  h2,
            const double  h3,
            const double  h4,
            const double  h5)
{
  return
    l1 < h1 && l1 < h2 && l1 < h3 && l1 < h4 && l1 < h5

    && l2 < h1 && l2 < h2 && l2 < h3 && l2 < h4 && l2 < h5

    && l3 < h1 && l3 < h2 && l3 < h3 && l3 < h4 && l3 < h5

    && l4 < h1 && l4 < h2 && l4 < h3 && l4 < h4 && l4 < h5

    && l5 < h1 && l5 < h2 && l5 < h3 && l5 < h4 && l5 < h5
  ;
}


PRIVATE short
intersectNodesBoundingBoxes(const AABB    *aabb1,
                            const AABB    *aabb2,
                            const stemBox *stem1,
                            const stemBox *stem2)
{
  const char  *fnName       = "intersectNodesBoundingBoxes";
  double      extraDistance = 0;

  extraDistance += epsilonRecognize;
  int         count = 0;

  if (stem1->bulgeDist > 0.0)

    count++;

  if (stem2->bulgeDist > 0.0)

    count++;

  if (count > 0)

    extraDistance += (1.0 / count) * (stem1->bulgeDist + stem2->bulgeDist);

  if (
    aabb1->max[0] < aabb2->min[0] - extraDistance
    ||
    aabb2->max[0] < aabb1->min[0] - extraDistance
    ||
    aabb1->max[1] < aabb2->min[1] - extraDistance
    ||
    aabb2->max[1] < aabb1->min[1] - extraDistance
    )
    return 0;
  else
    return 1;
}


PRIVATE intersectionType
intersectNodeNode(const treeNode  *node1,
                  const treeNode  *node2)
{
  const char  *fnName = "intersectNodeNode";
  int         bulge1  = -1;
  int         bulge2  = -1;

  if (node1 == node2)

    return noIntersection;

  stemBox *sBox_node1 = node1->sBox;
  loopBox *lBox_node1 = node1->lBox;
  stemBox *sBox_node2 = node2->sBox;
  loopBox *lBox_node2 = node2->lBox;

  short   intersect = intersectNodesBoundingBoxes(&(node1->aabb),
                                                  &(node2->aabb),
                                                  sBox_node1,
                                                  sBox_node2);

  if (!intersect)

    return noIntersection;

  treeNode  *parentOfNode1        = getParent(node1);
  treeNode  *parentOfNode2        = getParent(node2);
  short     node1IsParentOfNode2  = (node1 == parentOfNode2);
  short     node2IsParentOfNode1  = (node2 == parentOfNode1);
  short     nodesHaveCommonParent = (parentOfNode1 == parentOfNode2);

  /* SxS */
  if (!node1IsParentOfNode2
      && !node2IsParentOfNode1
      && !nodesHaveCommonParent
      && intersectStemStem(sBox_node1, sBox_node2))
    /* successive stems never intersect while config is not broken */
    return SxS;

  /* LxL */
  if (!node1IsParentOfNode2
      && !node2IsParentOfNode1
      && intersectLoopLoop(lBox_node1, lBox_node2))
    /* successive loops do never intersect */
    return LxL;

  /* SxL */
  if (!node2IsParentOfNode1
      && intersectStemLoop(sBox_node1, lBox_node2)
      )

    return SxL;

  /* LxS */
  if (!node1IsParentOfNode2
      && intersectStemLoop(sBox_node2, lBox_node1)
      )

    return LxS;

  /* LxB */
  if (!node1IsParentOfNode2) {
    if (intersectLoopBulges(lBox_node1, sBox_node2, &bulge2))

      return LxB;
  }

  /* BxL */
  if (!node2IsParentOfNode1) {
    if (intersectLoopBulges(lBox_node2, sBox_node1, &bulge1))

      return BxL;
  }

  /* SxB */
  if (intersectStemBulges(sBox_node1, sBox_node2, &bulge2))

    return SxB;

  /* BxS */
  if (intersectStemBulges(sBox_node2, sBox_node1, &bulge1))

    return BxS;

  /* BxB */
  if (intersectBulgesBulges(sBox_node1, sBox_node2, &bulge1, &bulge2))

    return BxB;

  return noIntersection;
}


PRIVATE short
intersectNodeTree(const treeNode  *node,
                  treeNode        *tree,
                  treeNode        **intersectorNode)
{
  intersectionType intersecting = intersectNodeNode(node, tree);

  if (intersecting != noIntersection) {
    *intersectorNode = tree;
    return 1;
  } else {
    int childCount = tree->childCount;

    for (int i = 0; i < childCount; i++) {
      if (intersectNodeTree(node, getChild(tree, i), intersectorNode))

        return 1;
    }
  }

  return 0;
}


PRIVATE short
intersect_iterateTree(treeNode  *tree1,
                      treeNode  *tree2,
                      treeNode  **intersectorNode1,
                      treeNode  **intersectorNode2)
{
  if (intersectNodeTree(tree1, tree2, intersectorNode2)) {
    *intersectorNode1 = tree1;
    return 1;
  } else {
    int childCount = tree1->childCount;

    for (int i = 0; i < childCount; i++) {
      if (intersect_iterateTree(getChild(tree1, i), tree2, intersectorNode1, intersectorNode2))

        return 1;
    }
  }

  return 0;
}


/*
 * starting method for detection of intersections between trees
 * basically this one iterates over both subtrees and does the intersection check
 * for each pair from tree1 and tree2
 * this is done recursively.
 * iterate over tree1 and for each node we iterate over tree2 for intersection calcultion
 */
PRIVATE short
intersectTrees(treeNode *tree1,
               treeNode *tree2)
{
  treeNode  *intersectorNode1 = NULL;
  treeNode  *intersectorNode2 = NULL;
  short     intersecting      = intersect_iterateTree(tree1,
                                                      tree2,
                                                      &intersectorNode1,
                                                      &intersectorNode2);

  return intersecting;
}


PRIVATE short
intersectNodeLists(treeNode *const                          *list1,
                   const int                                size1,
                   treeNode *const                          *list2,
                   const int                                size2,
                   const vrna_plot_options_puzzler_t *const puzzler)
{
  for (int index1 = 0; index1 < size1; index1++) {
    const treeNode  *node1      = list1[index1];
    short           isExterior1 = isExterior(node1);

    for (int index2 = 0; index2 < size2; index2++) {
      const treeNode *node2 = list2[index2];

      if (isExterior1) {
        if (intersectNodeExterior(node2, puzzler))
          return 1;
      } else

      if (isExterior(node2)) {
        if (intersectNodeExterior(node1, puzzler))
          return 1;
      } else {
        if (noIntersection != intersectNodeNode(node1, node2))

          return 1;
      }
    }
  }

  return 0;
}


#endif
