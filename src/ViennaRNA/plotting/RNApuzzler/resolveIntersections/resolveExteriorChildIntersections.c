#include "ViennaRNA/RNApuzzler/resolveIntersections/resolveExteriorChildIntersections.h"
#include "ViennaRNA/RNApuzzler/definitions.h"
#include "ViennaRNA/RNApuzzler/dataTypes/configtree_struct.h"
#include "ViennaRNA/RNApuzzler/data/configtree.h"
#include "ViennaRNA/RNApuzzler/data/boundingBoxes.h"
#include "ViennaRNA/RNApuzzler/intersectLevel/intersectLevelTreeNodes.h"
#include "ViennaRNA/utils.h"

#include <stdlib.h>

//-----------------------------------------------------------------------------

void getSimpleBoundingBox(
    treeNode* node,
    double* const bounds,
    const int recursionDepth
) {
    double loopMin = node->lBox->c[0] - node->lBox->r;
    double loopMax = node->lBox->c[0] + node->lBox->r;
    if (recursionDepth == 0) {
        bounds[0] = loopMin;
        bounds[1] = loopMax;
    }

    for (int i = 0; i < node->childCount; i++) {
        treeNode* child = getChild(node, i);
        getSimpleBoundingBox(child, bounds, recursionDepth + 1);
    }

    if (loopMin < bounds[0]) { bounds[0] = loopMin; }
    if (loopMax > bounds[1]) { bounds[1] = loopMax; }

    for (int i = 0; i < node->sBox->bulgeCount; i++) {
        double pPrev[2];
        double pThis[2];
        double pNext[2];
        getBulgeCoordinates(node->sBox, i, pPrev, pThis, pNext);

        if (pThis[0] < bounds[0]) { bounds[0] = pThis[0]; }
        if (pThis[0] > bounds[1]) { bounds[1] = pThis[0]; }
    }
}

/**
 * Resolve the intersections of the children of the exterior loop
 */
void resolveExteriorChildrenIntersectionXY(
        treeNode* exteriorNode,
        short const * const pair_table,
        const double unpaired,
        const short allowFlipping,
        double *myX,
        double *myY
) {
    // number of subtrees
    int subtreeCount = exteriorNode->childCount;

    if (subtreeCount < 2) {
        return;
    }

    /// for each exterior child: get first node
    treeNode** childTreeNode = (treeNode**) vrna_alloc(subtreeCount * sizeof(treeNode*));
    for (int subtree = 0; subtree < subtreeCount; subtree++) {
        childTreeNode[subtree] = getChild(exteriorNode, subtree);
    }

    /*
    /// for each exterior child: prepare bounding box
    double** bounds = (double**) vrna_alloc(subtreeCount * sizeof(double*));
    for (int subtree = 0; subtree < subtreeCount; subtree++) {
        bounds[subtree] = (double*) vrna_alloc(2 * sizeof(double));
        bounds[subtree][0] = 0.0;
        bounds[subtree][1] = 0.0;
    }

    /// get bounding box of first child
    getSimpleBoundingBox(childTreeNode[0], bounds[0], 0);
    */

    /// for each subtree
    /// - compute number of its first non-exterior base
    /// - compute number of nucleotides before the subtree
    /// - distance between nucleotides before the subtree
    int* firstBase = (int*) vrna_alloc(subtreeCount * sizeof(int));
    int* backbone = (int*) vrna_alloc(subtreeCount * sizeof(int));
    double* distance = (double*) vrna_alloc(subtreeCount * sizeof(double));
    for (int subtree = 0; subtree < subtreeCount; subtree++) {
        backbone[subtree] = 0;
        distance[subtree] = 0.0;
    }
    int subtree = 0;
    int base = 1;
    while (base < pair_table[0] && subtree < subtreeCount) {
        if (pair_table[base] > base) {
            firstBase[subtree] = base;
            subtree++;
            base = pair_table[base];
        } else {
            base++;
            backbone[subtree]++;
        }
    }

    // store upper and lower subtrees
    int* upper = (int*) vrna_alloc((subtreeCount + 1) * sizeof(int));
    int* lower = (int*) vrna_alloc((subtreeCount + 1) * sizeof(int));
    upper[0] = 0;
    lower[0] = 0;

    /// set first subtree to upper side
    upper[0]++;
    upper[upper[0]] = 0;

    // accumulated offset of children
    double offset = 0.0;
    // accumulated translation of children
    double accumulatedTranslation = 0.0;

    /// for all subtrees
    for (int subtree = 1; subtree < subtreeCount; subtree++) {
        /// translate current subtree by accumulated offset
        if (offset > 0.0) {
            double translate[2] = { offset, 0.0 };
            translateBoundingBoxes(childTreeNode[subtree], translate);
        }
        // getSimpleBoundingBox(childTreeNode[subtree], bounds[subtree], 0);

        /// as long as the current child gets translated
        short changed = 1;
        short intersectUpper = 0;
        short intersectLower = 0;
        double fixOverlap = 0.0;
        while (changed) {
            // printf("Handling subtree %d\n", subtree);
            changed = 0;
            intersectUpper = 0;
            intersectLower = 0;

            /// check intersection of current subtree with previous upper subtrees
            for (int u = 1; u <= upper[0]; u++) {
                int upperStem = upper[u];
                // printf("Check intersection between %d and %d\n", subtree, upperStem);
                intersectUpper = intersectTrees(childTreeNode[subtree], childTreeNode[upperStem]);
                if (intersectUpper) {
                    // printf("Intersection between %d and %d\n", subtree, upperStem);
                    break;
                }
                //printf("%d vs %d: upperOverlap:%f (boundsOverlap:%f)\n", subtree, upperStem, upperOverlap, boundsOverlap);
            }

            if (allowFlipping) {
                /// if flipping is allowed:
                /// check intersection of current subtree with previous lower subtrees
                for (int l = 1; l <= lower[0]; l++) {
                    int lowerStem = lower[l];
                    intersectLower = intersectTrees(childTreeNode[subtree], childTreeNode[lowerStem]);
                    if (intersectLower) {
                        break;
                    }
                    //printf("%d vs %d: lowerOverlap:%f (boundsOverlap:%f)\n", subtree, lowerStem, lowerOverlap, boundsOverlap);
                }
            }

            if ((!allowFlipping && intersectUpper) ||
                (allowFlipping && intersectUpper && intersectLower)) {
                /// if intersections can not be resolved by flipping
                /// increase distance by constant amount per exterior base
                distance[subtree] += unpaired; // minOverlap / backbone[subtree];
                fixOverlap = unpaired * backbone[subtree];
                // printf("Increase distance by %12.8lf\n", fixOverlap);

                double translate[2] = { fixOverlap, 0.0 };
                translateBoundingBoxes(childTreeNode[subtree], translate);

                /*
                bounds[subtree][0] += fixOverlap;
                bounds[subtree][1] += fixOverlap;
                */

                offset += fixOverlap;
                // printf("Total offset: %12.8lf\n", offset);

                changed = 1;
            } else {
                if (allowFlipping && intersectUpper) {
                    /// if flipping is allowed and sufficient for resolving the intersection:
                    /// (intersection is on the upper side)
                    lower[0]++;
                    lower[lower[0]] = subtree;
                } else {
                    upper[0]++;
                    upper[upper[0]] = subtree;
                }
            }
        } // end while(changed)

        /// translate exterior bases between previous and current subtree
        int currentBase = 1;
        for (int base = pair_table[firstBase[subtree-1]];
             base < firstBase[subtree];
             base++, ++currentBase) {
            myX[base] += currentBase * distance[subtree] + accumulatedTranslation;
        }
        accumulatedTranslation += distance[subtree] * backbone[subtree];
    }

    /// Last part of the exterior loop
    for (int base = pair_table[firstBase[subtreeCount-1]];
         base < pair_table[0];
         ++base) {
        myX[base] += accumulatedTranslation;
    }

    /// modify x- and y-coordinates for all subtrees
    int currentLower = 1;
    double translation = 0.0;
    for (int subtree = 1; subtree < subtreeCount; subtree++) {
        /// translate all bases of current subtree
        translation += distance[subtree] * backbone[subtree];
        // printf("Translate subtree %d by %12.8lf\n", subtree, translation);
        for (int base = firstBase[subtree];
             base < pair_table[firstBase[subtree]];
             ++base) {
            myX[base] += translation;
        }

        if (subtree == lower[currentLower]) {
            /// flip subtrees
            double exteriorY = myY[1];
            for (int base = firstBase[subtree];
                 base < pair_table[firstBase[subtree]];
                 ++base) {
                myY[base] = 2 * exteriorY - myY[base];
            }
            ++currentLower;
        }
    }
    // processing end

    free(upper);
    free(lower);
    free(backbone);
    free(distance);
    free(firstBase);
    /*
    for (int subtree = 0; subtree < subtreeCount; subtree++) {
        free(bounds[subtree]);
    }
    free(bounds);
    */
    free(childTreeNode);
}

/**
 * Resolve the intersections of the children of the exterior loop
 */
void resolveExteriorChildrenIntersectionAffin(
        treeNode* exteriorNode,
        short const * const pair_table,
        tBaseInformation* const baseInformation,
        const double unpaired,
        const short allowFlipping
) {
    // number of subtrees
    int subtreeCount = exteriorNode->childCount;

    if (subtreeCount < 2) {
        return;
    }

    /// for each exterior child: get first node
    treeNode** childTreeNode = (treeNode**) vrna_alloc(subtreeCount * sizeof(treeNode*));
    for (int subtree = 0; subtree < subtreeCount; subtree++) {
        childTreeNode[subtree] = getChild(exteriorNode, subtree);
    }

    /*
    /// for each exterior child: prepare bounding box
    double** bounds = (double**) vrna_alloc(subtreeCount * sizeof(double*));
    for (int subtree = 0; subtree < subtreeCount; subtree++) {
        bounds[subtree] = (double*) vrna_alloc(2 * sizeof(double));
        bounds[subtree][0] = 0.0;
        bounds[subtree][1] = 0.0;
    }

    /// get bounding box of first child
    getSimpleBoundingBox(childTreeNode[0], bounds[0], 0);
    */

    /// for each subtree
    /// - compute number of its first non-exterior base
    /// - compute number of nucleotides before the subtree
    int* firstBase = (int*) vrna_alloc(subtreeCount * sizeof(int));
    int* backbone = (int*) vrna_alloc(subtreeCount * sizeof(int));
    for (int subtree = 0; subtree < subtreeCount; subtree++) {
        backbone[subtree] = 0;
    }
    int subtree = 0;
    int base = 1;
    while (base < pair_table[0] && subtree < subtreeCount) {
        if (pair_table[base] > base) {
            firstBase[subtree] = base;
            subtree++;
            base = pair_table[base];
        } else {
            base++;
            backbone[subtree]++;
        }
    }

    // store upper and lower stems
    int* upper = (int*) vrna_alloc((subtreeCount + 1) * sizeof(int));
    int* lower = (int*) vrna_alloc((subtreeCount + 1) * sizeof(int));
    upper[0] = 0;
    lower[0] = 0;

    /// set first subtree to upper side
    upper[0]++;
    upper[upper[0]] = 0;

    // accumulated offset of children
    double offset = 0.0;

    /// for all stems
    for (int subtree = 1; subtree < subtreeCount; subtree++) {
        /// translate current subtree by accumulated offset
        if (offset > 0.0) {
            double translate[2] = { offset, 0.0 };
            translateBoundingBoxes(childTreeNode[subtree], translate);
        }
        // getSimpleBoundingBox(childTreeNode[subtree], bounds[subtree], 0);

        /// as long as the current child gets translated
        short changed = 1;
        short intersectUpper = 0;
        short intersectLower = 0;
        double fixOverlap = 0.0;
        while (changed) {
            changed = 0;
            intersectUpper = 0;
            intersectLower = 0;

            /// check intersection of current subtree with previous upper stems
            for (int u = 1; u <= upper[0]; u++) {
                int upperStem = upper[u];
                intersectUpper = intersectTrees(childTreeNode[subtree], childTreeNode[upperStem]);
                if (intersectUpper) {
                    break;
                }
                //printf("%d vs %d: upperOverlap:%f (boundsOverlap:%f)\n", subtree, upperStem, upperOverlap, boundsOverlap);
            }

            if (allowFlipping) {
                /// if flipping is allowed:
                /// check intersection of current subtree with previous lower stems
                for (int l = 1; l <= lower[0]; l++) {
                    int lowerStem = lower[l];
                    intersectLower = intersectTrees(childTreeNode[subtree], childTreeNode[lowerStem]);
                    if (intersectLower) {
                        break;
                    }
                    //printf("%d vs %d: lowerOverlap:%f (boundsOverlap:%f)\n", subtree, lowerStem, lowerOverlap, boundsOverlap);
                }
            }

            if ((!allowFlipping && intersectUpper) ||
                (allowFlipping && intersectUpper && intersectLower)) {
                /// if intersections can not be resolved by flipping
                /// increase distance by constant amount per exterior base
                double deltaPerDistance = unpaired; // minOverlap / backbone[subtree];
                fixOverlap = deltaPerDistance * backbone[subtree];

                for (int base = pair_table[firstBase[subtree-1]]; base < firstBase[subtree]; base++) {
                    baseInformation[base].distance += deltaPerDistance;
                }

                double translate[2] = { fixOverlap, 0.0 };
                translateBoundingBoxes(childTreeNode[subtree], translate);

                /*
                bounds[subtree][0] += fixOverlap;
                bounds[subtree][1] += fixOverlap;
                */

                offset += fixOverlap;

                changed = 1;
            } else {
                if (allowFlipping && intersectUpper) {
                    /// if flipping is allowed and sufficient for resolving the intersection:
                    /// (intersection is on the upper side)
                    for (int base = firstBase[subtree] + 1; base <= pair_table[firstBase[subtree]] + 1; base++) {
                        if (base > pair_table[0]) {
                            break;
                        }
                        baseInformation[base].angle *= -1;
                    }

                    lower[0]++;
                    lower[lower[0]] = subtree;
                } else {
                    upper[0]++;
                    upper[upper[0]] = subtree;
                }
            }
        } // end while(changed)
    }
    // processing end

    free(upper);
    free(lower);
    free(backbone);
    free(firstBase);
    /*
    for (int subtree = 0; subtree < subtreeCount; subtree++) {
        free(bounds[subtree]);
    }
    free(bounds);
    */
    free(childTreeNode);
}

void resolveExteriorChildIntersections(
        treeNode* exteriorNode,
        short const * const pair_table,
        tBaseInformation* const baseInformation,
        const double unpaired,
        const short allowFlipping
) {
    int stemCount = exteriorNode->childCount;

    if (stemCount < 2) {
        return;
    }

    /// for each exterior child: get first node
    treeNode** node = (treeNode**) vrna_alloc(stemCount * sizeof(treeNode*));
    for (int stem = 0; stem < stemCount; stem++) {
        node[stem] = getChild(exteriorNode, stem);
    }

    /// for each exterior child: prepare bounding box
    double** bounds = (double**) vrna_alloc(stemCount * sizeof(double*));
    for (int stem = 0; stem < stemCount; stem++) {
        bounds[stem] = (double*) vrna_alloc(2 * sizeof(double));
        bounds[stem][0] = 0.0;
        bounds[stem][1] = 0.0;
    }
    getSimpleBoundingBox(node[0], bounds[0], 0);

    /// for each stem
    /// - compute number of its first base
    /// - compute number of nucleotides before the stem
    int* firstBase = (int*) vrna_alloc(stemCount * sizeof(int));
    int* backbone = (int*) vrna_alloc(stemCount * sizeof(int));
    for (int stem = 0; stem < stemCount; stem++) {
        backbone[stem] = 0;
    }
    int stem = 0;
    int base = 1;
    while (base < pair_table[0] && stem < stemCount) {
        if (pair_table[base] > base) {
            firstBase[stem] = base;
            stem++;
            base = pair_table[base];
        } else {
            base++;
            backbone[stem]++;
        }
    }

    // store upper and lower stems
    int* upper = (int*) vrna_alloc((stemCount + 1) * sizeof(int));
    int* lower = (int*) vrna_alloc((stemCount + 1) * sizeof(int));
    upper[0] = 0;
    lower[0] = 0;

    /// set first stem to upper side
    upper[0]++;
    upper[upper[0]] = 0;

    // accumulated offset of children
    double offset = 0.0;

    // processing start
    for (int stem = 1; stem < stemCount; stem++) {
        // initial setting
        if (offset > 0.0) {
            double translate[2] = { offset, 0.0 };
            translateBoundingBoxes(node[stem], translate);
        }
        getSimpleBoundingBox(node[stem], bounds[stem], 0);

        /// as long as the current child gets translated
        short changed = 1;
        while (changed) {
            changed = 0;

            // get overlap with upper and lower side
            double upperOverlap = 0.0;
            for (int u = 1; u <= upper[0]; u++) {
                int upperStem = upper[u];
                double boundsOverlap = (bounds[upperStem][1] + unpaired) - bounds[stem][0];
                if (boundsOverlap > upperOverlap && intersectTrees(node[stem], node[upperStem])) {
                    upperOverlap = boundsOverlap;
                }
                //printf("%d vs %d: upperOverlap:%f (boundsOverlap:%f)\n", stem, upperStem, upperOverlap, boundsOverlap);
            }
            double lowerOverlap = 0.0;
            for (int l = 1; l <= lower[0]; l++) {
                int lowerStem = lower[l];
                double boundsOverlap = (bounds[lowerStem][1] + unpaired) - bounds[stem][0];
                if (boundsOverlap > lowerOverlap && intersectTrees(node[stem], node[lowerStem])) {
                    lowerOverlap = boundsOverlap;
                }
                //printf("%d vs %d: lowerOverlap:%f (boundsOverlap:%f)\n", stem, lowerStem, lowerOverlap, boundsOverlap);
            }

            // fix minimum of both overlaps
            double minOverlap = (allowFlipping && (lowerOverlap < upperOverlap)) ? lowerOverlap : upperOverlap;
            if (minOverlap > 0.0) {
                double deltaPerDistance = unpaired; // minOverlap / backbone[stem];
                minOverlap = deltaPerDistance * backbone[stem];
                for (int base = pair_table[firstBase[stem-1]]; base < firstBase[stem]; base++) {
                    baseInformation[base].distance += deltaPerDistance;
                }

                double translate[2] = { minOverlap, 0.0 };
                translateBoundingBoxes(node[stem], translate);

                bounds[stem][0] += minOverlap;
                bounds[stem][1] += minOverlap;

                offset += minOverlap;

                changed = 1;
            } else {
                // flip if necessary
                if (lowerOverlap < upperOverlap) {

                    for (int base = firstBase[stem] + 1; base <= pair_table[firstBase[stem]] + 1; base++) {
                        if (base > pair_table[0]) {
                            break;
                        }
                        baseInformation[base].angle *= -1;
                    }

                    lower[0]++;
                    lower[lower[0]] = stem;
                } else {
                    upper[0]++;
                    upper[upper[0]] = stem;
                }
            }
        } // end while(changed)
    }
    // processing end

    free(upper);
    free(lower);
    free(backbone);
    free(firstBase);
    for (int stem = 0; stem < stemCount; stem++) {
        free(bounds[stem]);
    }
    free(bounds);
    free(node);
}

