#include "ViennaRNA/RNApuzzler/intersectLevel/intersectLevelTreeNodes.h"
#include "ViennaRNA/RNApuzzler/definitions.h"
#include "ViennaRNA/RNApuzzler/data/cfg_reader.h"
#include "ViennaRNA/RNApuzzler/data/configtree.h"
#include "ViennaRNA/RNApuzzler/data/boundingBoxes.h"
#include "ViennaRNA/RNApuzzler/intersectLevel/intersectLevelBoundingBoxes.h"
#include "ViennaRNA/RNApuzzler/resolveIntersections/intersectionType.h"
#include "ViennaRNA/RNApuzzler/output/output.h"

#include <stdio.h>

short intersectNodeExterior(
        const treeNode* node,
        const puzzlerOptions* puzzler
) {
    if (isExterior(node)) {
        return 0;
    }
    if (isExterior(getParent(node))) {
        return 0;
    }

    double cy = node->lBox->c[1];
    double r = node->lBox->r + epsilonRecognize;
    if (puzzler->checkExteriorIntersections) {
        return (cy - r) <= EXTERIOR_Y;
    } else {
        return 0;
    }
}

short checkBounds(
    const double l1, const double l2, const double l3, const double l4, const double l5,
    const double h1, const double h2, const double h3, const double h4, const double h5
) {
    return
        (   l1 < h1 && l1 < h2 && l1 < h3 && l1 < h4 && l1 < h5

         && l2 < h1 && l2 < h2 && l2 < h3 && l2 < h4 && l2 < h5

         && l3 < h1 && l3 < h2 && l3 < h3 && l3 < h4 && l3 < h5

         && l4 < h1 && l4 < h2 && l4 < h3 && l4 < h4 && l4 < h5

         && l5 < h1 && l5 < h2 && l5 < h3 && l5 < h4 && l5 < h5
        );
}

short OLD_intersectNodesBoundingBoxes(
    const stemBox *stem1,
    const loopBox *loop1,
    const stemBox *stem2,
    const loopBox *loop2
) {
    const char *fnName = "OLD_intersectNodesBoundingBoxes";

    const double stem1_ea[2] = { stem1->e[0] * stem1->a[0],
                                 stem1->e[0] * stem1->a[1] };
    const double stem1_eb[2] = { stem1->e[1] * stem1->b[0],
                                 stem1->e[1] * stem1->b[1] };
    const double A1[2] = { stem1->c[0] - stem1_ea[0] + stem1_eb[0],
                           stem1->c[1] - stem1_ea[1] + stem1_eb[1] };
    const double B1[2] = { stem1->c[0] + stem1_ea[0] + stem1_eb[0],
                           stem1->c[1] + stem1_ea[1] + stem1_eb[1] };
    const double C1[2] = { stem1->c[0] + stem1_ea[0] - stem1_eb[0],
                           stem1->c[1] + stem1_ea[1] - stem1_eb[1] };
    const double D1[2] = { stem1->c[0] - stem1_ea[0] - stem1_eb[0],
                           stem1->c[1] - stem1_ea[1] - stem1_eb[1] };

    const double loop1L[2] = { loop1->c[0] - loop1->r,
                               loop1->c[1] - loop1->r };
    const double loop1H[2] = { loop1->c[0] + loop1->r,
                               loop1->c[1] + loop1->r };

    const double stem2_ea[2] = { stem2->e[0] * stem2->a[0],
                                 stem2->e[0] * stem2->a[1] };
    const double stem2_eb[2] = { stem2->e[1] * stem2->b[0],
                                 stem2->e[1] * stem2->b[1] };
    const double A2[2] = { stem2->c[0] - stem2_ea[0] + stem2_eb[0],
                           stem2->c[1] - stem2_ea[1] + stem2_eb[1] };
    const double B2[2] = { stem2->c[0] + stem2_ea[0] + stem2_eb[0],
                           stem2->c[1] + stem2_ea[1] + stem2_eb[1] };
    const double C2[2] = { stem2->c[0] + stem2_ea[0] - stem2_eb[0],
                           stem2->c[1] + stem2_ea[1] - stem2_eb[1] };
    const double D2[2] = { stem2->c[0] - stem2_ea[0] - stem2_eb[0],
                           stem2->c[1] - stem2_ea[1] - stem2_eb[1] };

    const double loop2L[2] = { loop2->c[0] - loop2->r,
                               loop2->c[1] - loop2->r };
    const double loop2H[2] = { loop2->c[0] + loop2->r,
                               loop2->c[1] + loop2->r };

    double extraDistance = 0;
    extraDistance += epsilonRecognize;
    int count = 0;
    if (stem1->bulgeDist > 0.0) { count++; }
    if (stem2->bulgeDist > 0.0) { count++; }
    if (count > 0) {
      extraDistance += (1.0 / count) * (stem1->bulgeDist + stem2->bulgeDist);
    }

/*
    printInformation(fnName, "1 x %12.8lf %12.8lf %12.8lf %12.8lf %12.8lf %12.8lf\n", A1[0], B1[0],  C1[0], D1[0], loop1L[0], loop1H[0]);
    printInformation(fnName, "1 y %12.8lf %12.8lf %12.8lf %12.8lf %12.8lf %12.8lf\n", A1[1], B1[1],  C1[1], D1[1], loop1L[1], loop1H[1]);
    printInformation(fnName, "2 x %12.8lf %12.8lf %12.8lf %12.8lf %12.8lf %12.8lf\n", A2[0], B2[0],  C2[0], D2[0], loop2L[0], loop2H[0]);
    printInformation(fnName, "2 y %12.8lf %12.8lf %12.8lf %12.8lf %12.8lf %12.8lf\n", A2[1], B2[1],  C2[1], D2[1], loop2L[1], loop2H[1]);
*/

    if (
        checkBounds(A1[0], B1[0], C1[0], D1[0], loop1H[0],
                    A2[0] - extraDistance,
                    B2[0] - extraDistance,
                    C2[0] - extraDistance,
                    D2[0] - extraDistance,
                    loop2L[0] - extraDistance
                   )
        ||
        checkBounds(A1[1], B1[1], C1[1], D1[1], loop1H[1],
                    A2[1] - extraDistance,
                    B2[1] - extraDistance,
                    C2[1] - extraDistance,
                    D2[1] - extraDistance,
                    loop2L[1] - extraDistance
                   )
        ||
        checkBounds(A2[0] + extraDistance,
                    B2[0] + extraDistance,
                    C2[0] + extraDistance,
                    D2[0] + extraDistance,
                    loop2H[0] + extraDistance,
                    A1[0], B1[0], C1[0], D1[0], loop1L[0]
                   )
        ||
        checkBounds(A2[1] + extraDistance,
                    B2[1] + extraDistance,
                    C2[1] + extraDistance,
                    D2[1] + extraDistance,
                    loop2H[1] + extraDistance,
                    A1[1], B1[1], C1[1], D1[1], loop1L[1]
                   )
    ) {
      return 0;
    }

    return 1;
}

short intersectNodesBoundingBoxes(
    const AABB    *aabb1,
    const AABB    *aabb2,
    const stemBox *stem1,
    const stemBox *stem2
) {
    const char *fnName = "intersectNodesBoundingBoxes";

    double extraDistance = 0;
    extraDistance += epsilonRecognize;
    int count = 0;
    if (stem1->bulgeDist > 0.0) { count++; }
    if (stem2->bulgeDist > 0.0) { count++; }
    if (count > 0) {
      extraDistance += (1.0 / count) * (stem1->bulgeDist + stem2->bulgeDist);
    }

    //printInformation(fnName, "aabb1 min-max: %12.8lf %12.8lf -- %12.8lf %12.8lf\n", aabb1->min[0], aabb1->min[1], aabb1->max[0], aabb1->max[1]);
    //printInformation(fnName, "aabb2 min-max: %12.8lf %12.8lf -- %12.8lf %12.8lf\n", aabb2->min[0], aabb2->min[1], aabb2->max[0], aabb2->max[1]);

    if (
        aabb1->max[0] < aabb2->min[0] - extraDistance
        ||
        aabb2->max[0] < aabb1->min[0] - extraDistance
        ||
        aabb1->max[1] < aabb2->min[1] - extraDistance
        ||
        aabb2->max[1] < aabb1->min[1] - extraDistance
        ) {
        return 0;
    } else {
        return 1;
    }
}

intersectionType intersectNodeNode(
        const treeNode* node1,
        const treeNode* node2
) {
    const char *fnName = "intersectNodeNode";

    int bulge1 = -1;
    int bulge2 = -1;

    if (node1 == node2) {
        return noIntersection;
    }

    stemBox* sBox_node1 = node1->sBox;
    loopBox* lBox_node1 = node1->lBox;
    stemBox* sBox_node2 = node2->sBox;
    loopBox* lBox_node2 = node2->lBox;

    // short intersectOld = OLD_intersectNodesBoundingBoxes(sBox_node1, lBox_node1, sBox_node2, lBox_node2);
    short intersect = intersectNodesBoundingBoxes(&(node1->aabb), &(node2->aabb), sBox_node1, sBox_node2);
    /*
    if (intersectOld != intersect) {
        printInformation(fnName, "aabb1 %d min-max: %12.8lf %12.8lf -- %12.8lf %12.8lf\n", getNodeID(node1), node1->aabb.min[0], node1->aabb.min[1], node1->aabb.max[0], node1->aabb.max[1]);
        printInformation(fnName, "aabb2 %d min-max: %12.8lf %12.8lf -- %12.8lf %12.8lf\n", getNodeID(node2), node2->aabb.min[0], node2->aabb.min[1], node2->aabb.max[0], node2->aabb.max[1]);
        printInformation(fnName, "result: %d\n", intersect);
    }
    */

    if (!intersect) {
        return noIntersection;
    }

    treeNode* parentOfNode1 = getParent(node1);
    treeNode* parentOfNode2 = getParent(node2);
    short node1IsParentOfNode2 = (node1 == parentOfNode2);
    short node2IsParentOfNode1 = (node2 == parentOfNode1);
    short nodesHaveCommonParent = (parentOfNode1 == parentOfNode2);

    /// SxS
    if (!node1IsParentOfNode2
        && !node2IsParentOfNode1
        && !nodesHaveCommonParent
        && intersectStemStem(sBox_node1, sBox_node2)) {
        /// successive stems never intersect while config is not broken
        return SxS;
    }

    /// LxL
    if (!node1IsParentOfNode2
        && !node2IsParentOfNode1
        && intersectLoopLoop(lBox_node1, lBox_node2)) {
        /// successive loops do never intersect
        return LxL;
    }

    /// SxL
    if (!node2IsParentOfNode1
        && intersectStemLoop(sBox_node1, lBox_node2)
       ) {
        return SxL;
    }

    /// LxS
    if (!node1IsParentOfNode2
        && intersectStemLoop(sBox_node2, lBox_node1)
       ) {
        return LxS;
    }

    /// LxB
    if (!node1IsParentOfNode2) {
        if (intersectLoopBulges(lBox_node1, sBox_node2, &bulge2)) {
            return LxB;
        }
    }

    /// BxL
    if (!node2IsParentOfNode1) {
        if (intersectLoopBulges(lBox_node2, sBox_node1, &bulge1)) {
            return BxL;
        }
    }

    /// SxB
    if (intersectStemBulges(sBox_node1, sBox_node2, &bulge2)) {
        return SxB;
    }

    /// BxS
    if (intersectStemBulges(sBox_node2, sBox_node1, &bulge1)) {
        return BxS;
    }

    /// BxB
    if (intersectBulgesBulges(sBox_node1, sBox_node2, &bulge1, &bulge2)) {
        return BxB;
    }

    return noIntersection;
}

short intersectNodeTree(
    const treeNode* node,
    treeNode* tree,
    treeNode** intersectorNode
) {
    intersectionType intersecting = intersectNodeNode(node, tree);

    if (intersecting != noIntersection) {
        *intersectorNode = tree;
        return 1;
    } else {
        int childCount = tree->childCount;
        for (int i = 0; i < childCount; i++) {
            if (intersectNodeTree(node, getChild(tree, i), intersectorNode)) {
                return 1;
            }
        }
    }

    return 0;
}

short intersect_iterateTree(
        treeNode* tree1,
        treeNode* tree2,
        treeNode** intersectorNode1,
        treeNode** intersectorNode2
) {
    if (intersectNodeTree(tree1, tree2, intersectorNode2)) {
        *intersectorNode1 = tree1;
        return 1;
    } else {
        int childCount = tree1->childCount;
        for (int i = 0; i < childCount; i++) {
            if (intersect_iterateTree(getChild(tree1, i), tree2, intersectorNode1, intersectorNode2)) {
                return 1;
            }
        }
    }

    return 0;
}

/*
 * starting method for detection of intersections between trees
 * basically this one iterates over both subtrees and does the intersection check
 * for each pair from tree1 and tree2
 * this is done recursively...
 * iterate over tree1 and for each node we iterate over tree2 for intersection calcultion
 */
short intersectTrees(
        treeNode* tree1,
        treeNode* tree2
) {
    treeNode* intersectorNode1 = NULL;
    treeNode* intersectorNode2 = NULL;
    short intersecting = intersect_iterateTree(tree1, tree2, &intersectorNode1, &intersectorNode2);

    return intersecting;
}

short intersectNodeLists(
        treeNode* const * list1,
        const int size1,
        treeNode* const * list2,
        const int size2,
        const puzzlerOptions* const puzzler
) {
    for (int index1 = 0; index1 < size1; index1++) {

        const treeNode* node1 = list1[index1];
        short isExterior1 = isExterior(node1);

        for (int index2 = 0; index2 < size2; index2++) {

            const treeNode* node2 = list2[index2];

            if (isExterior1) {
                if (intersectNodeExterior(node2, puzzler)) {
                    return 1;
                }
            } else
            if (isExterior(node2)) {
                if (intersectNodeExterior(node1, puzzler)) {
                    return 1;
                }
            } else {
                if (noIntersection != intersectNodeNode(node1, node2)) {
                    return 1;
                }
            }

        }

    }

    return 0;
}
