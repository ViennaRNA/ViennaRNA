#ifndef RNAPUZZLER_CONFIGTREE_H
#define RNAPUZZLER_CONFIGTREE_H

#include "ViennaRNA/plotting/RNApuzzler/dataTypes/tBaseInformation_struct.h"
#include "ViennaRNA/plotting/RNApuzzler/dataTypes/config_struct.h"
#include "ViennaRNA/plotting/RNApuzzler/dataTypes/AABB_struct.h"
#include "ViennaRNA/plotting/RNApuzzler/dataTypes/boundingBoxes_struct.h"
#include "ViennaRNA/plotting/RNApuzzler/dataTypes/configtree_struct.h"
#include "ViennaRNA/plotting/RNApuzzler/definitions.h"

/*--------------------------------------------------------------------------*/

PRIVATE void
updateAABB(AABB     *aabb,
           stemBox  *sBox,
           loopBox  *lBox);


PRIVATE short
isExterior(const treeNode *node);


PRIVATE treeNode *
getExterior(treeNode *node);


PRIVATE short
isInteriorLoop(const treeNode *node);


PRIVATE short
isMultiLoop(const treeNode *node);


/**
 * @brief buildConfigtree
 * @param pair_table
 *      - definition of the RNAs pairing
 * @param baseInformation
 *      - array of tBaseInformation that grants information about the RNA
 *        (only config is needed here)
 * @param x
 * @param y
 * @param bulge
 *      - distance between regular stem and bulge base
 * @return
 *      - a configtree build from the given RNA pairing.
 *        gives information about all configurations (config) applied to the RNA.
 */
PRIVATE treeNode *
buildConfigtree(const short *const      pair_table,
                const tBaseInformation  *baseInformation,
                const double            *x,
                const double            *y,
                const double            bulge);


PRIVATE void
updateBoundingBoxes(treeNode              *node,
                    const puzzlerOptions  *puzzler);


/**
 * @brief applyChangesToConfigAndBoundingBoxes
 *      Function to apply a set of config changes to the given treenode.
 *
 *      By passing -1.0 to radiusNew you will apply the minimum possible radius for this loop while not allowing to shrink the loop.
 *      In case it would actually shrink a default relative increase is applied to the radius.
 *
 *      By passing 0.0 to radiusNew you set the minimum possible radius no matter if there would be some shrinking or not.
 * @param tree the node to change
 * @param deltaCfg array of angles to change (in degree)
 * @param radiusNew desired radius to set while applying deltas
 * @param puzzler
 * @return 1 if the changes were applied, 0 otherwise
 */
PRIVATE void
applyChangesToConfigAndBoundingBoxes(treeNode             *tree,
                                     const double         *deltaCfg,
                                     const double         radiusNew,
                                     const puzzlerOptions *puzzler);


/**
 * @brief translateBoundingBoxes
 *      - Performs a translation of a whole branch by a given vector.
 *        Used to apply changes in config to the tree and its bounding boxes.
 * @param tree
 *      - tree that is being translated
 * @param vector
 *      - translation vector as double array[2]
 */
PRIVATE void
translateBoundingBoxes(treeNode     *tree,
                       const double *vector);


/**
 * @brief getChildIndex
 *      - gets the index of child node where to find the node with given ID.
 * @param tree
 *      - configtree you want to search in.
 * @param childID
 *      - ID of childnode you are looking for.
 * @return
 *      - child index or -1 if tree does not contain such a childnode.
 */
PRIVATE int
getChildIndex(const treeNode  *tree,
              const int       childID);


/**
 * @brief getChildAngle
 *      - Calculates the angle of a given child node (given by ID) at a loop.
 *        This child does not need to be a direct child node but can be any sort of grand child.
 *        If the given node is a child of the loop, the rotation angle of its center node
 *        will be calculated.
 * @param root
 *      - tree node acting as parent loop
 * @param node
 *      - child node of which you want to know the angle
 * @return
 *      - angle of child node
 *        the resulting angle might be smaller than 0° or greater than 360°
 */
PRIVATE double
getChildAngle(const treeNode  *root,
              const treeNode  *child);


/**
 * @brief getChildAngleByIndex
 *      - Calculates the clockwise angle of a given child node (given by name) at a loop.
 *        This child needs to be a direct child node
 *        The rotation angle of its center node will be calculated.
 * @param parentNode
 *      - tree node acting as parent loop
 * @param childIndex
 *      - index of child node of which you want to know the angle
 * @return
 *      - angle of child node
 *        the resulting angle might be smaller than 0° or greater than 360°
 */
PRIVATE double
getChildAngleByIndex(const treeNode *parentNode,
                     const int      childIndex);


/**
 * @brief getLoopCenter
 *      - Getter for the center coordinates of a tree node's loop.
 * @param node
 *      - your tree node
 * @param p
 *      - double[2] return value for the loop's center coordinates
 */
PRIVATE void
getLoopCenter(const treeNode  *node,
              double          p[2]);


/**
 * @brief getStemCenter
 *      - Getter for the center coordinates of a tree node's stem.
 * @param node
 *      - your tree node
 * @param p
 *      - double[2] return value for the stem's center coordinates
 */
PRIVATE void
getStemCenter(const treeNode  *node,
              double          p[2]);


/**
 * @brief getChildNode
 *      - searches for the node with given childID in a configtree.
 * @param tree
 *      - configtree you want to search in.
 * @param childID
 *      - ID of childnode you are looking for.
 * @return
 *      - ptr to node or NULL if tree does not contain such a childnode.
 */
PRIVATE treeNode *
getChildNode(treeNode   *tree,
             const int  childID);


PRIVATE treeNode *
getParent(const treeNode *node);


PRIVATE treeNode *
getChild(const treeNode *node,
         const int      index);


PRIVATE char
getNodeName(const treeNode *node);


PRIVATE int
getNodeID(const treeNode *node);


PRIVATE void
printTree(const treeNode  *node,
          const int       level);


PRIVATE int
countSubtreeNodes(const treeNode *node);


PRIVATE int
collectSubtreeNodes(treeNode          *node,
                    treeNode **const  allNodes,
                    const int         index);


PRIVATE int
countAncestorNodes(const treeNode *node);


PRIVATE void
collectAncestorNodes(const treeNode   *node,
                     treeNode **const ancestorList);


PRIVATE double
getPairedAngle(const treeNode *node);


PRIVATE void
freeTree(treeNode *node);


#include "configtree.inc"

#endif
