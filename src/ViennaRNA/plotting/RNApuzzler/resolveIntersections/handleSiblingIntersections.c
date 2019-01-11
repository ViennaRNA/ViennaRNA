#include "ViennaRNA/RNApuzzler/resolveIntersections/handleSiblingIntersections.h"
#include "ViennaRNA/RNApuzzler/dataTypes/configtree_struct.h"
#include "ViennaRNA/RNApuzzler/data/configtree.h"
#include "ViennaRNA/RNApuzzler/definitions.h"
#include "ViennaRNA/RNApuzzler/vector_math.h"
#include "ViennaRNA/RNApuzzler/intersectLevel/intersectLevelTreeNodes.h"
#include "ViennaRNA/RNApuzzler/resolveIntersections/boundingWedge.h"
#include "ViennaRNA/RNApuzzler/resolveIntersections/calcDeltas.h"
#include "ViennaRNA/RNApuzzler/resolveIntersections/handleConfigChanges.h"
#include "ViennaRNA/RNApuzzler/output/output.h"
#include "ViennaRNA/RNApuzzler/output/configtree_debug.h"

#include "ViennaRNA/utils.h"

#include <stdlib.h>
#include <math.h>

/**
 * @brief fixIntersectionOfSiblings
 *      Try to fix intersections of sibling subtrees at their common ancestor.
 * @param tree
 *      common ancestor of intersecting subtrees
 * @param left
 * @param right
 * @param deltaCfg
 * @param puzzler
 * @return
 *      1 if something was changed, 0 otherwise
 */
short fixIntersectionOfSiblings(
    treeNode* tree,
    const int left,
    const int right,
    double *deltaCfg,
    puzzlerOptions* puzzler
) {
    double wedgeMin, wedgeMax;
    getBoundingWedge(tree, right, &wedgeMin, &wedgeMax);
    double minAngle = wedgeMin;
    getBoundingWedge(tree, left, &wedgeMin, &wedgeMax);
    double maxAngle = wedgeMax;
    double targetAngle = minAngle - maxAngle;

    short changed = 0;
    if (targetAngle < 0) {
        targetAngle = fmax(targetAngle, -MATH_PI_HALF); // limit the angle to avoid malformed structures
        double changedAngle = calcDeltas(tree, getParent(tree), left, right, (-1) * targetAngle, puzzler, deltaCfg);

        if (changedAngle != 0.0) {

            treeNode* leftNode  = getChild(tree, left);
            treeNode* rightNode = getChild(tree, right);
            if (FANCY_PS) {
                PS_printFancySiblings(tree, leftNode, rightNode, puzzler);
            }

            // apply all changes
            changed = checkAndApplyConfigChanges(tree, deltaCfg, siblings, puzzler);

            if (FANCY_PS) {
                PS_printFancySiblings(tree, leftNode, rightNode, puzzler);
            }

            //printf("[%s] changed: %d\n", fnName, changed);
        }
    }

    return changed;
}

/**
 * @brief handleIntersectionOfSiblings
 *      Try to fix intersections of sibling subtrees at their common ancestor.
 * @param tree
 *      common ancestor of intersecting subtrees
 * @param listOfIntersections
 *      list of pairs of indices of subtrees that are intersecting.
 *      [numberOfIntersections, [i1, j1], [i2, j2], ...]
 * @return
 *      1 if something was changed, 0 otherwise
 */
short handleIntersectionOfSiblings(
        treeNode* tree,
        const int* listOfIntersections,
        puzzlerOptions* puzzler
) {
    char* fnName = "FIX INTERSECTION OF SIBLINGS";

    /// idea:
    /// - measure each intersection by calculating an overlap angle
    /// - increase the spaces between intersectors and decrease spaces that are not between them
    /// - distribute the overlap angle equally to all participating spaces
    /// - check for new or remaining intersections (at the end)

    if (puzzler->numberOfChangesAppliedToConfig > puzzler->maximumNumberOfConfigChangesAllowed) {
        printError(fnName, "Reached maximum number of changes. Abort.\n");
        return -1;
    }

    short changed = 0;
    int intersectionCount = listOfIntersections[0];

        /*
        printf("[%s] Summary: [", fnName);
        for (int i = 0; i < intersectionCount; i++) {
            int left = listOfIntersections[2*i+1];
            int right = listOfIntersections[2*i+2];
            treeNode* childLeft  = getChild(tree, left);
            treeNode* childRight = getChild(tree, right);
            if (i > 0) {
                printf("\n[%s]           ", fnName);
            }
            printf(" %d[%d]=%d vs. %d[%d]=%d"
                   , getNodeID(tree), left, getID(childLeft)
                   , getNodeID(tree), right, getID(childRight)
                   );
        }
        printf(" ]\n");
        */

    int childCount = tree->childCount;
    int configSize = childCount + 1;

    double* deltaCfg   = (double*) vrna_alloc(configSize * sizeof (double));

    /// init deltas with zero
    for (int i = 0; i < configSize; i++) {
        deltaCfg[i] = 0.0;
    }

    /// fix intersections of siblings
    for (int k = 0; k < intersectionCount; k++) { // for all intersections - start
        int left  = listOfIntersections[2 * k + 1];
        int right = listOfIntersections[2 * k + 2];

        /*
        treeNode* childLeft  = getChild(tree, left);
        treeNode* childRight = getChild(tree, right);
        printf("[%s] %d[%d]=%d vs. %d[%d]=%d\n"
               , fnName
               , getNodeID(tree), left, getNodeID(childLeft)
               , getNodeID(tree), right, getNodeID(childRight)
               );
        */

        changed = fixIntersectionOfSiblings(tree, left, right, deltaCfg, puzzler);
        if (changed) {
            break;
        }
    }

    free(deltaCfg);

    return changed;
}

short checkSiblings(
        treeNode* node,
        puzzlerOptions* puzzler
) {
    char* fnName = "CHECK SIBLINGS";
    short ret = _false;

    int childCount = node->childCount;
    /// create array to store all information about overlapping neighbors
    int* intersectorsBranches = (int*) vrna_alloc((childCount * childCount) * sizeof (int));
    for (int i = 0; i < childCount * childCount; i++) {
        intersectorsBranches[i] = -1;
    }

    /// actually check for those intersections ...
    for (int i = 0; i < childCount; i++) {
        int intersectorsCount = 0;
        for (int j = i+1; j < childCount; j++) {
            treeNode* childI = getChild(node, i);
            treeNode* childJ = getChild(node, j);
            if (intersectTrees(childI, childJ)) {
                intersectorsBranches[i * childCount + intersectorsCount] = j;
                intersectorsCount++;
            }
        }
    }

    /// ... and count them
    int intersectionCount = 0;
    for (int i = 0; i < (childCount * childCount); i++) {
        if (intersectorsBranches[i] != -1) {
            intersectionCount++;
        }
    }

    if (intersectionCount > 0) {
        ret |= _intersect;

        /// transform intersection information into format
        /// [ count, [intersector_a, intersector_b], [intersector_a, intersector_b], ... ]
        /// where count states how many pairs of a/b are there
        /// the i-th intersection has index a=2*i+1; b=2*i+2
        int* listOfIntersections = (int*) vrna_alloc((1 + 2 * intersectionCount) * sizeof (int));
        listOfIntersections[0] = intersectionCount;

        int counter = 0;
        for (int i = 0; i < (childCount * childCount); i++) {
            if (intersectorsBranches[i] != -1) {
                listOfIntersections[2 * counter + 1] = i / childCount;
                listOfIntersections[2 * counter + 2] = intersectorsBranches[i];
                counter++;
            }
        }

        /// resolve all of those intersections for this node's subtrees
        short retFix = handleIntersectionOfSiblings(node, listOfIntersections, puzzler);
        if (retFix < 0) {
            ret = retFix;
        } else if (retFix) {
            ret |= _changed;
        }

        free(listOfIntersections);
    }

    free(intersectorsBranches);

    return ret;
}

