#include "ViennaRNA/RNApuzzler/resolveIntersections/handleAncestorIntersections.h"
#include "ViennaRNA/RNApuzzler/dataTypes/configtree_struct.h"
#include "ViennaRNA/RNApuzzler/data/configtree.h"
#include "ViennaRNA/RNApuzzler/data/boundingBoxes.h"
#include "ViennaRNA/RNApuzzler/definitions.h"
#include "ViennaRNA/RNApuzzler/resolveIntersections/intersectionType.h"
#include "ViennaRNA/RNApuzzler/resolveIntersections/handleConfigChanges.h"
#include "ViennaRNA/RNApuzzler/resolveIntersections/calcDeltas.h"
#include "ViennaRNA/RNApuzzler/resolveIntersections/rotationAngle.h"
#include "ViennaRNA/RNApuzzler/intersectLevel/intersectLevelTreeNodes.h"
#include "ViennaRNA/RNApuzzler/intersectLevel/intersectLevelLines.h"
#include "ViennaRNA/RNApuzzler/vector_math.h"
#include "ViennaRNA/RNApuzzler/output/output.h"
#include "ViennaRNA/RNApuzzler/output/configtree_debug.h"

#include "ViennaRNA/utils.h"

#include <stdlib.h>
#include <math.h>

/**
 * - -1 if rotation is counter-clockwise
 * - 1 if rotation is clockwise
 * - 0 else
 */
short TENTATIVE3_getRotationSign(
        treeNode** const       path,
        const int              pathLength,
        const intersectionType it
) {
    char* fnName = "getRotationSign";
    short result = 0;

    if (pathLength < 2) {
        printError(fnName, "[CRITICAL] path length: %d\n", pathLength);
        return result;
    }

    /// Compute rotation angle
    double angle = 0.0;
    treeNode *currentNode = path[0];
    treeNode *childNode;
    for (int i = 1; i < pathLength; i++) {
        childNode = path[i];
        angle += getChildAngle(currentNode, childNode);
//        printInformation(fnName, "child angle %d:  %12.8lf\n", i, angle);
        angle -= MATH_PI;
//        printInformation(fnName, "turtle angle %d: %12.8lf\n", i, angle);
        currentNode = childNode;
    }
//    printInformation(fnName, "angle: %12.8lf\n", angle);

    /// add bulge angle to angle sum
    switch(it) {
    case BxL:
        if (isToTheRightPointPoint(path[0]->sBox->c, path[0]->lBox->c, path[pathLength-1]->sBox->c)) {
            angle += MATH_PI_HALF;
        } else {
            angle -= MATH_PI_HALF;
        }
        break;
    case LxB:
        if (isToTheRightPointPoint(path[pathLength-1]->sBox->c, path[pathLength-1]->lBox->c, path[0]->sBox->c)) {
            angle += MATH_PI_HALF;
        } else {
            angle -= MATH_PI_HALF;
        }
        break;
    }

    /// Determine return value
    if (angle < 0.0) {
        result = 1;
    } else if (angle > 0.0) {
        result = -1;
    } else {
        result = 0;
    }

//    printDebug(fnName, "return: %d (%12.8lf°)\n", result, angle);
    return result;
}

/**
 * - -1 if rotation is counter-clockwise
 * - 1 if rotation is clockwise
 * - 0 else
 */
short TENTATIVE2_getRotationSign(
        treeNode** const path,
        const int        pathLength
) {
    char* fnName = "getRotationSign";
    short result = 0;

    if (pathLength < 2) {
        printError(fnName, "[CRITICAL] path length: %d\n", pathLength);
        return result;
    }

    /// Compute rotation angle
    double angle = 0.0;
    treeNode *currentNode = path[0];
    treeNode *childNode;
    for (int i = 1; i < pathLength; i++) {
        childNode = path[i];
        angle += getChildAngle(currentNode, childNode);
//        printInformation(fnName, "child angle %d:  %12.8lf\n", i, angle);
        angle -= MATH_PI;
//        printInformation(fnName, "turtle angle %d: %12.8lf\n", i, angle);
        currentNode = childNode;
    }
//    printInformation(fnName, "angle: %12.8lf\n", angle);

    /// Determine return value
    if (angle < 0.0) {
        result = 1;
    } else if (angle > 0.0) {
        result = -1;
    } else {
        result = 0;
    }

//    printDebug(fnName, "return: %d (%12.8lf°)\n", result, angle);
    return result;
}

/**
 * - -1 if rotation is counter-clockwise
 * - 1 if rotation is clockwise
 * - 0 else
 */
short getRotationSign(
        const treeNode** const path,
        const int pathLength
) {
    char* fnName = "getRotationSign";

    short result = 0;

    if (pathLength < 2) {
        printError(fnName, "[CRITICAL] path length: %d\n", pathLength);
        return result;
    }

    double angle = 0.0;

    if (pathLength == 2) {
        /// path has length 2
        /// -> only one loop and intersector
        /// -> get angle of child of this loop
//        printDebug(fnName, "[PATH LENTH 2] [%s %5d] %d elements\n", puzzler->filename, puzzler->numberOfChangesAppliedToConfig, pathLength);
//        printPath(fnName, path, pathLength, -1);
//        PS_printPath(path[0], path[pathLength - 1], puzzler);
        angle = getChildAngle(path[0], path[1]);
        angle -= MATH_PI;
    } else {
        /// path has more than two usable loops

        /// center of the angle calculation
        double center[2];
        getLoopCenter(path[1], center);

        /// reference point for the axis center->refPoint
        double refPoint[2];
        refPoint[0] = center[0] + 1000 * path[1]->sBox->a[0];
        refPoint[1] = center[1] + 1000 * path[1]->sBox->a[1];

        /// initial axis equals 0°
        /// left  rotations get negative angles
        /// right rotations get positive angles

        for (int i = 2; i < pathLength; i++) {
            /// get the next node's center ...
            double point[2];
            getLoopCenter(path[i], point);

            /// update angle
            double diffAngle = anglePtPtPt2D(refPoint, center, point);
            if (!isToTheRightPointPoint(center, refPoint, point)) {
                diffAngle *= -1;
            }
            angle += diffAngle;

            /// update the reference axis
            refPoint[0] = point[0];
            refPoint[1] = point[1];
        }
    }

//    printDebug(fnName, "angle: %+7.2f\n", angle);

    if (angle > 0.0) {
        result = -1;
    }
    if (angle < 0.0) {
        result = 1;
    }
//    printDebug(fnName, "return: %d\n", result);

    if (result == 0) {
        printError(fnName, "[CRITICAL] there should always be a non-zero rotation sign. Maybe non-intersecting call?\n");
    }

    return result;
}

treeNode *fixIntersectionWithAncestor(
        treeNode* ancestor,
        treeNode* rotationNode,
        treeNode* intersector,
        const int rotationIndex,
        const short rotationSign,
        const intersectionType it,
        puzzlerOptions* puzzler
) {

    char* fnName = "fixIntersectionWithAncestor";

    if (rotationNode == ancestor && (it == LxL || it == LxS || it == LxB)) {
//        printf("[%s] [%d %d %d] no-op at ancestor with type %s\n", fnName, getNodeID(ancestor), getNodeID(rotationNode), getNodeID(intersector), intersectionTypeToString(it));
        return 0;
    }

    double interiorChildAngle;
    if (isInteriorLoop(rotationNode)) {
        /// prevent interior loops from increasing the distance to their "straight" state
        /// We only allow rotations towards straight.
        interiorChildAngle = getChildAngleByIndex(rotationNode, 0);
        short allowedRotationSign = 0;
        if (interiorChildAngle > MATH_PI) {
            allowedRotationSign = -1;
        } else if (interiorChildAngle < MATH_PI) {
            allowedRotationSign = 1;
        }
        if (rotationSign != allowedRotationSign) {
            return 0;
        }
    }

    /// get the rotation angle depending on the intersection type
    double rotationAngle = getRotationAngle(ancestor, rotationNode, intersector, it, rotationSign);

    if (isInteriorLoop(rotationNode)) {
        /// prevent interior loops from rotating over their "straight" state
        /// If necessary we limit the rotation so that this interior loop becomes straight.
        double diffToStraight = MATH_PI - interiorChildAngle;
        if (fabs(rotationAngle) > fabs(diffToStraight)) {
            rotationAngle = diffToStraight;
        }
    }

    short changed = 0;
    if (rotationAngle != 0.0) {
        /// compute the deltas for changing the configuration of this loop
        double* deltas = (double*) vrna_alloc((rotationNode->childCount + 1) * sizeof(double));
        double deltaAngle = fabs(rotationAngle);
        int indexLeft = -2;
        int indexRight = -2;
        if (rotationAngle > 0.0) {
            indexLeft = -1;
            indexRight = rotationIndex;
        } else {
            indexLeft = rotationIndex;
            indexRight = -1;
        }

        calcDeltas(rotationNode, ancestor, indexLeft, indexRight, deltaAngle, puzzler, deltas);

        /// check and apply the computed changes to the configuration of this loop
        intersectionType itLog = (isExterior(ancestor) ? exterior : it);
        changed = checkAndApplyConfigChanges(rotationNode, deltas, itLog, puzzler);

        if (FANCY_PS) {
            PS_printFancyPath(ancestor, intersector, rotationNode, puzzler);
        }

        free(deltas);
    }

    if (changed) {
        return rotationNode;
    } else {
        return NULL;
    }
}

/**
 * check if loop is interior and straight
 */
short isStraightInteriorLoop(
        treeNode *node
) {
    return
      (isInteriorLoop(node)
       && (getChildAngleByIndex(node, 0) == MATH_PI));
}

/* DEPRECATED:
 * construct complete path between ancestor and intersector including both
 */
treeNode** constructIntersectionPath(
        treeNode *ancestor,
        treeNode *intersector,
        intersectionType it,
        int *pathLength
) {
    /// - compute path length
    *pathLength = 1;
    treeNode* node = intersector;
    while (node != ancestor) {
        node = getParent(node);
        ++(*pathLength);
    }

    switch (it) {
        case LxL:
        case LxS:
        case LxB:
            /// - start at the child of the ancestor node for Lx? intersections
            if (!isExterior(ancestor)) {
                --(*pathLength);
            }
            break;
        default:
            /// - include ancestor node for all other intersections (if not straight interior loop)
            ;
    }

    /// - construct path
    treeNode** path = (treeNode**) vrna_alloc((*pathLength) * sizeof(treeNode*));
    node = intersector;
    for (int i = (*pathLength) - 1; i >= 0; node = getParent(node)) {
        path[i] = node;
        i--;
    }

    return path;
}

/**
 * compute intersection path from ancestor to intersector
 * skipping straight interior loops that can not be used for rotations
 * -> improves rotation angle and rotation loop computations
 */
treeNode** constructReducedIntersectionPath(
        treeNode *ancestor,
        treeNode *intersector,
        intersectionType it,
        int *pathLength
) {
    /// - compute path length
    *pathLength = 1;
    treeNode* node = intersector;
    while (node != ancestor) {
        node = getParent(node);
        if (!isStraightInteriorLoop(node)) {
            /// skip straight interior loops
            ++(*pathLength);
        }
    }

    switch (it) {
        case LxL:
        case LxS:
        case LxB:
            /// - start at the child of the ancestor node for Lx? intersections
            if (!isStraightInteriorLoop(ancestor)) {
                /// exclude ancestor from reduced path
                --(*pathLength);
            }
            break;
        default:
            /// - include ancestor node for all other intersections (if not straight interior loop)
            ;
    }

    /// - construct path
    treeNode** path = (treeNode**) vrna_alloc((*pathLength) * sizeof(treeNode*));
    node = intersector;
    for (int i = (*pathLength) - 1; i >= 0; node = getParent(node)) {
        if ((i == *pathLength - 1) || !isStraightInteriorLoop(node)) {
            path[i] = node;
            i--;
        }
    }

    return path;
}

/**
 */
treeNode* handleIntersectionWithAncestor(
        treeNode* ancestor,
        treeNode* intersector,
        const int recursionDepth,
        puzzlerOptions* puzzler
) {
    char* fnName = "handleIntersectionWithAncestor";
//    printf("[%s] %d vs. %d [START]\n", fnName, getNodeID(ancestor), getNodeID(intersector));

    /// Determine intersection type
    intersectionType it = intersectNodeNode(ancestor, intersector);
    if (it == noIntersection) {
        /// Early termination for no intersection.
        printError(fnName,
                   "wrong input. there is no intersection for %d and %d. [%s_DEBUG_PATH_%05d_%04d_vs_%04d.ps]\n",
                   getNodeID(ancestor),
                   getNodeID(intersector),
                   puzzler->filename,
                   puzzler->numberOfChangesAppliedToConfig,
                   getNodeID(ancestor),
                   getNodeID(intersector));
        PS_printPath(ancestor, intersector, puzzler);
        return NULL;
    }

    /// construct path from ancestor to intersector
    int pathLength = 0;
    //treeNode** path = constructIntersectionPath(ancestor, intersector, it, &pathLength);
    treeNode** path = constructReducedIntersectionPath(ancestor, intersector, it, &pathLength);

    /// at each node determine child to follow
    int* childIndex = (int*) vrna_alloc((pathLength-1) * sizeof(int));
    for (int i = 0; i < pathLength - 1; i++) {
        childIndex[i] = getChildIndex(path[i], getNodeID(path[i+1]));
    }

    /// node changed by trying to fix intersection
    treeNode *changedNode = NULL;

    /// Compute orientation of reduced path
    //short rotationSign = getRotationSign(path, pathLength);
    short rotationSign = TENTATIVE2_getRotationSign(path, pathLength);
    //short tentative_rotationSign = TENTATIVE2_getRotationSign(path, pathLength, it);

    /* DZ: check
    if (rotationSign != tentative_rotationSign && !isExterior(ancestor)) {
        printWarning(fnName,
                     "Old rotation sign != new rotation sign: %d != %d; child : %d\n",
                     rotationSign,
                     tentative_rotationSign,
                     getChildIndex(ancestor, getNodeID(intersector))
                    );
        PS_printFancyPath(ancestor, intersector, NULL, puzzler);
    }
    */

    if (rotationSign == 0) {
        /// path is straight -> no intersection
        printError(fnName, "[FAILED] invalid rotation sign (zero) for %04d and %04d. [%s_DEBUG_PATH_%05d_%04d_vs_%04d.ps]\n",
                   getNodeID(ancestor),
                   getNodeID(intersector),
                   puzzler->filename,
                   puzzler->numberOfChangesAppliedToConfig,
                   getNodeID(ancestor),
                   getNodeID(intersector));
//        PS_printPath(ancestor, intersector, puzzler);
    } else {

//    printf("[%s] Path of %d vs. %d ...", fnName, getNodeID(ancestor), getNodeID(intersector));
//    for (int i = 0; i < pathLength; i++) {
//        printf(" %d(%d)", getNodeID(path[i]), path[i]->childCount);
//    }
//    printf("\n");

        if (FANCY_PS) {
            PS_printFancyPath(ancestor, intersector, NULL, puzzler);
        }

        /// run from intersector to ancestor twice:
        /// - in the first run we only choose interior loops for rotations
        int nodeNumber = pathLength - 2; // skip intersector, start with its first ancestor
        while ((changedNode == NULL) && (nodeNumber >= 0)) {
            if (isInteriorLoop(path[nodeNumber])) {
                if (FANCY_PS) {
                    PS_printFancyPath(ancestor, intersector, path[nodeNumber], puzzler);
                }

                /// interior loop: try to improve
                changedNode = fixIntersectionWithAncestor(ancestor, path[nodeNumber], intersector, childIndex[nodeNumber], rotationSign, it, puzzler);
            }

            /// go to non-straight ancestor
            nodeNumber--;
        }

        /// - in the second run we only choose multi loop for rotations
        nodeNumber = pathLength - 2; // skip intersector, start with its first ancestor
        while ((changedNode == NULL) && (nodeNumber >= 0)) {
            if (isMultiLoop(path[nodeNumber])) {
                if (FANCY_PS) {
                    PS_printFancyPath(ancestor, intersector, path[nodeNumber], puzzler);
                }

                /// multi-loop: try to improve
                changedNode = fixIntersectionWithAncestor(ancestor, path[nodeNumber], intersector, childIndex[nodeNumber], rotationSign, it, puzzler);
            }

            /// go to non-straight ancestor
            nodeNumber--;
        }
    }

    if (changedNode == NULL) {
        printError(fnName, "[FAILED] to resolve %04d vs. %04d (%s) [%s_DEBUG_PATH_%05d_%04d_vs_%04d.ps]\n",
                   getNodeID(ancestor),
                   getNodeID(intersector),
                   (changedNode != NULL ? "changed" : "not changed"),
                   puzzler->filename,
                   puzzler->numberOfChangesAppliedToConfig,
                   getNodeID(ancestor),
                   getNodeID(intersector)
                  );
        printPath(fnName, path, pathLength, -1);
        PS_printPath(ancestor, intersector, puzzler);
    }

    free(path);
    free(childIndex);

//    printf("[%s] %d vs. %d [RETURN %d]\n", fnName, getNodeID(ancestor), getNodeID(intersector), ret);
    return changedNode;
}

void TENTATIVE2_updateExteriorBoundingBoxes(
        treeNode* exterior,
        loopBox* loop,
        double stemNorthX,
        double stemSouthX,
        double stemWestY,
        double stemEastY
) {
    double s[2]  = {stemSouthX, stemWestY}; // south / west corner
    double e[2]  = {stemNorthX, stemWestY}; // north / west corner
    double sp[2] = {stemSouthX, stemEastY}; // south / east corner
    stemBox* stem = createStemBox(s, e, sp);

    if (exterior->lBox) { free(exterior->lBox); }
    if (exterior->sBox) { free(exterior->sBox); }
    exterior->lBox = loop;
    exterior->sBox = stem;
    loop->parent = exterior;
    stem->parent = exterior;

    updateAABB(&(exterior->aabb), stem, loop);
}

void TENTATIVE3_setupExteriorBoundingBoxes(
        treeNode* exterior,
        const treeNode* topLevelAncestor,
        const treeNode* intersector,
        const puzzlerOptions *puzzler
) {
    char* fnName = "setupExteriorBoundingBoxes";

    /// determine loop
    double upperY = EXTERIOR_Y;
    // double lowerY = 0.0;
    double lowerY = upperY - puzzler->paired;

    double loopX = topLevelAncestor->lBox->c[0];

    double radius = 0.5 * (upperY - lowerY);
    double center[2] = {loopX, upperY - radius};

    /// determine stem
    double stemNorthX = loopX;
    double stemSouthX = loopX;

    const AABB *intersectorAABB = &(intersector->aabb);
    if (intersectorAABB->max[0] < loopX) {
        /// intersector is left of topLevelAncestor
        /// use distance aabb->min[0] .. loopX for stem setup
        double stemNorthX = loopX;
        double stemSouthX = intersectorAABB->min[0];
        double stemWestY = upperY;
        double StemEastY = lowerY;

        loopBox* loop = createLoopBox(center, radius);
        TENTATIVE2_updateExteriorBoundingBoxes(
            exterior,
            loop,
            stemNorthX,
            stemSouthX,
            stemWestY,
            StemEastY
        );
    } else if (loopX < intersectorAABB->min[0]) {
        /// intersector is right of topLevelAncestor
        /// use distance loopX .. aabb->max[0] for stem setup
        double stemNorthX = loopX;
        double stemSouthX = intersectorAABB->max[0];
        double stemWestY = lowerY;
        double StemEastY = upperY;

        loopBox* loop = createLoopBox(center, radius);
        TENTATIVE2_updateExteriorBoundingBoxes(
            exterior,
            loop,
            stemNorthX,
            stemSouthX,
            stemWestY,
            StemEastY
        );
    } else {
        /// intersector shares some space in x direction with topLevelAncestor
        /// use distance aabb->min[0] .. loopX for stem setup
        /// then check of intersection
        /// if noIntersection
        /// use distance loopX .. aabb->max[0] for stem setup
        double stemNorthX = loopX;
        double stemSouthX = intersectorAABB->min[0];
        double stemWestY = upperY;
        double StemEastY = lowerY;

        loopBox* loop = createLoopBox(center, radius);
        TENTATIVE2_updateExteriorBoundingBoxes(
            exterior,
            loop,
            stemNorthX,
            stemSouthX,
            stemWestY,
            StemEastY
        );

        if (noIntersection == intersectNodeNode(intersector, exterior)) {
            double stemNorthX = loopX;
            double stemSouthX = intersectorAABB->max[0];
            double stemWestY = lowerY;
            double StemEastY = upperY;

            loopBox* loop = createLoopBox(center, radius);
            TENTATIVE2_updateExteriorBoundingBoxes(
                exterior,
                loop,
                stemNorthX,
                stemSouthX,
                stemWestY,
                StemEastY
            );
        }
    }
}

void TENTATIVE_updateExteriorBoundingBoxes(
        treeNode* exterior,
        loopBox* loop,
        double stemBottom,
        double stemTop,
        double stemLeft,
        double stemRight
) {
    double s[2]  = {stemBottom, stemRight};
    double e[2]  = {stemTop, stemRight};
    double sp[2] = {stemBottom, stemLeft};
    stemBox* stem = createStemBox(s, e, sp);

    if (exterior->lBox) { free(exterior->lBox); }
    if (exterior->sBox) { free(exterior->sBox); }
    exterior->lBox = loop;
    exterior->sBox = stem;
    loop->parent = exterior;
    stem->parent = exterior;

    updateAABB(&(exterior->aabb), stem, loop);
}

void TENTATIVE2_setupExteriorBoundingBoxes(
        treeNode* exterior,
        const treeNode* topLevelAncestor,
        const treeNode* intersector,
        const puzzlerOptions *puzzler
) {
    char* fnName = "setupExteriorBoundingBoxes";

    /// determine loop
    double upperY = EXTERIOR_Y;
    // double lowerY = 0.0;
    double lowerY = upperY - puzzler->paired;

    double loopX = topLevelAncestor->lBox->c[0];

    double radius = 0.5 * (upperY - lowerY);
    double center[2] = {loopX, upperY - radius};

    /// determine stem
    double stemMinX = loopX;
    double stemMaxX = loopX;

    const AABB *intersectorAABB = &(intersector->aabb);
    intersectionType it = noIntersection;
    if ((intersectorAABB->min)[0] < stemMinX) {
        // setup left side of ancestor
        stemMinX = (intersectorAABB->min)[0];

        loopBox* loop = createLoopBox(center, radius);
        TENTATIVE_updateExteriorBoundingBoxes(
            exterior,
            loop,
            stemMinX,
            stemMaxX,
            upperY,
            lowerY
        );

        // check left side for intersection
        it = intersectNodeNode(intersector, exterior);
    }

    if (it == noIntersection) {
        // left side of ancestor does not intersect
        // -> setup right side of ancestor
        stemMinX = loopX;
        if ((intersectorAABB->max)[0] > stemMaxX) {
            stemMaxX = (intersectorAABB->max)[0];
        }

        // recreate loopBox (previous one will be destroyed)
        loopBox *loop = createLoopBox(center, radius);
        TENTATIVE_updateExteriorBoundingBoxes(
            exterior,
            loop,
            stemMaxX,
            stemMinX,
            lowerY,
            upperY
        );
    }
}

void TENTATIVE_setupExteriorBoundingBoxes(
        treeNode* exterior,
        const treeNode* topLevelAncestor,
        const treeNode* intersector,
        const puzzlerOptions *puzzler
) {
    char* fnName = "setupExteriorBoundingBoxes";

    /// determine loop
    double upperY = EXTERIOR_Y;
    // double lowerY = 0.0;
    double lowerY = upperY - puzzler->paired;

    double loopX = topLevelAncestor->lBox->c[0];

    double radius = 0.5 * (upperY - lowerY);
    double center[2] = {loopX, upperY - radius};

    /// determine stem
    double stemMinX = loopX;
    double stemMaxX = loopX;
    loopBox* loop = createLoopBox(center, radius);

    const AABB *intersectorAABB = &(intersector->aabb);
    intersectionType it = noIntersection;
    if ((intersectorAABB->min)[0] < stemMinX) {
        stemMinX = (intersectorAABB->min)[0];
    }
    if ((intersectorAABB->max)[0] > stemMaxX) {
        stemMaxX = (intersectorAABB->max)[0];
    }

    double s[2]  = {stemMinX, upperY};
    double e[2]  = {stemMaxX, upperY};
    double sp[2] = {stemMinX, lowerY};
    stemBox* stem = createStemBox(s, e, sp);

    if (exterior->lBox) { free(exterior->lBox); }
    if (exterior->sBox) { free(exterior->sBox); }
    exterior->lBox = loop;
    exterior->sBox = stem;
    loop->parent = exterior;
    stem->parent = exterior;

    updateAABB(&(exterior->aabb), stem, loop);
}

void setupExteriorBoundingBoxes(
        treeNode* exterior,
        const treeNode* topLevelAncestor,
        const treeNode* intersector,
        const puzzlerOptions *puzzler
) {
    char* fnName = "setupExteriorBoundingBoxes";

    /// get bounding box for exterior
    treeNode* parentOfIntersector = getParent(intersector);
    stemBox* refStem = intersector->sBox;
    loopBox* refLoop = intersector->lBox;
    double minX = refStem->c[0];
    double maxX = minX;
    double x[4];
    x[0] = parentOfIntersector->lBox->c[0] + parentOfIntersector->lBox->r;
    x[1] = parentOfIntersector->lBox->c[0] - parentOfIntersector->lBox->r;
    x[2] = refLoop->c[0] + refLoop->r;
    x[3] = refLoop->c[0] - refLoop->r;
    for (int i = 0; i < 4; i++) {
        minX = fmin(minX, x[i]);
        maxX = fmax(maxX, x[i]);
    }
    minX -= epsilonFix;
    maxX += epsilonFix;

    double upperY = EXTERIOR_Y;
    double lowerY = 0.0;
    // double lowerY = upperY - puzzler->paired;

    double topLevelX = topLevelAncestor->lBox->c[0];
    double s[2], e[2], sp[2];

    /// (1) iParent above direct child
    /// (1.1) iNode left of direct child
    /// (1.2) iNode right of direct child
    /// (2) iParent and iNode on same side of direct child
    /// (2.1) both left of direct child
    /// (2.1) both right of direct child
    /// (3) iParent and iNode on different sides of direct child
    /// (3.1) exterior cut point is between iParent and direct child
    /// (3.1) exterior cut point is between iNode and direct child

    if (fabs(parentOfIntersector->lBox->c[0] - topLevelX) < epsilon0) {
        /// (1)
        if (refLoop->c[0] < topLevelX) {
            /// (1.1)
//            printf("[%s] case 1.1\n", fnName);
            s[0]  = x[3];         s[1]  = lowerY;
            e[0]  = topLevelX;    e[1]  = lowerY;
            sp[0] = x[3];         sp[1] = upperY;
        } else
        if (refLoop->c[0] > topLevelX) {
            /// (1.2)
//            printf("[%s] case 1.2\n", fnName);
            s[0]  = x[2];         s[1]  = lowerY;
            e[0]  = topLevelX;    e[1]  = lowerY;
            sp[0] = x[2];         sp[1] = upperY;
        } else {
            /// this would mean rootNode and intersectorNode coincide...
            /// since this would cause an interior intersection
            /// we just ignore this case here while giving the following hint.
            printError(fnName, "unrecognized case (1.#).\n");
        }
    } else {
        short sameSide = (refLoop->c[0] - topLevelX < 0) == (parentOfIntersector->lBox->c[0] - topLevelX < 0);
        if (sameSide) {
            /// (2)
            if (parentOfIntersector->lBox->c[0] < topLevelX) {
                /// (2.1)
//                printf("[%s] case 2.1\n", fnName);
                // take minX
                s[0]  = minX;         s[1]  = upperY;
                e[0]  = topLevelX;    e[1]  = upperY;
                sp[0] = minX;         sp[1] = lowerY;
            } else
            if (parentOfIntersector->lBox->c[0] > topLevelX) {
                /// (2.2)
//                printf("[%s] case 2.2\n", fnName);
                // take maxX
                s[0]  = maxX;         s[1]  = lowerY;
                e[0]  = topLevelX;    e[1]  = lowerY;
                sp[0] = maxX;         sp[1] = upperY;
            } else {
                printError(fnName, "unrecognized case (2.#).\n");
            }
        } else {
            /// (3)
            double l1p1[2], l1p2[2], l2p1[2], l2p2[2];
            double pParent[2];
            getLoopCenter(parentOfIntersector, pParent);
            double pIntersector[2];
            getLoopCenter(intersector, pIntersector);
            if (pIntersector[1] < upperY) {
                l1p2[0] = pIntersector[0];
                l1p2[1] = pIntersector[1];
            } else {
                double dx = pIntersector[0] - pParent[0];
                double dy = pIntersector[1] - pParent[1];
                double scale = (lowerY - pParent[1]) / dy;
                l1p2[0] = pParent[0] + scale * dx;
                l1p2[1] = pParent[1] + scale * dy;
            }

            l1p1[0] = pParent[0];
            l1p1[1] = pParent[1];

            l2p1[0] = pParent[0];
            l2p1[1] = upperY;
            l2p2[0] = topLevelX;
            l2p2[1] = upperY;

//            GEOGEBRA_printLinePointPoint("lParent", l1p1, l1p2);
//            GEOGEBRA_printLinePointPoint("exterior", l2p1, l2p2);
            short intersectBetweenRootxAndParent = intersectLineSegments(l1p1, l1p2, l2p1, l2p2, NULL);
            if (intersectBetweenRootxAndParent) {
                /// (3.1)
//                printf("[%s] case 3.1\n", fnName);
                if (parentOfIntersector->lBox->c[0] < topLevelX) {
                    s[0]  = x[1];         s[1]  = upperY;
                    e[0]  = topLevelX;    e[1]  = upperY;
                    sp[0] = x[1];         sp[1] = lowerY;
                } else
                if (parentOfIntersector->lBox->c[0] > topLevelX) {
                    s[0]  = x[0];         s[1]  = lowerY;
                    e[0]  = topLevelX;    e[1]  = lowerY;
                    sp[0] = x[0];         sp[1] = upperY;
                } else {
                    printError(fnName, "unrecognized case (3.1.#).\n");
                }
            } else {
                /// (3.2)
//                printf("[%s] case 3.2\n", fnName);
                if (refLoop->c[0] < topLevelX) {
                    s[0]  = x[3];         s[1]  = lowerY;
                    e[0]  = topLevelX;    e[1]  = lowerY;
                    sp[0] = x[3];         sp[1] = upperY;
                } else
                if (refLoop->c[0] > topLevelX) {
                    s[0]  = x[2];         s[1]  = lowerY;
                    e[0]  = topLevelX;    e[1]  = lowerY;
                    sp[0] = x[2];         sp[1] = upperY;
                } else {
                    printError(fnName, "unrecognized case (3.2.#).\n");
                }
            }
        }

    }

    stemBox* stem = createStemBox(s, e, sp);
    stem->bulgeCount = 0;
    stem->bulgeDist = 0.0;

    double radius = 0.5 * (upperY - lowerY);
    double center[2] = {topLevelX, upperY - radius};
    loopBox* loop = createLoopBox(center, radius);

    if (exterior->lBox) { free(exterior->lBox); }
    if (exterior->sBox) { free(exterior->sBox); }
    exterior->lBox = loop;
    exterior->sBox = stem;
    loop->parent = exterior;
    stem->parent = exterior;

    updateAABB(&(exterior->aabb), stem, loop);
}

treeNode *checkNodeAgainstAncestors(
        treeNode* node,
        puzzlerOptions* puzzler
) {
    char* fnName = "checkNodeAgainstAncestors";

    treeNode *changedNode = NULL;

    treeNode* ancestor = getParent(node);
    treeNode* topLevelAncestor = node;

    /// Move towards root checking and fixing ancestor intersections
    while (!isExterior(ancestor)) {
        topLevelAncestor = ancestor;
        intersectionType it = intersectNodeNode(node, ancestor);
        if (it != noIntersection) {
            changedNode = handleIntersectionWithAncestor(ancestor, node, 0, puzzler);
            if (changedNode != NULL) {
                return changedNode;
            }
        }
        ancestor = getParent(ancestor);
        //printf("[ANCESTOR] %d vs %d : it:%d changed:%d\n", getNodeID(ancestor), getNodeID(node), it, (result & _changed));
    }

    /// Check and fix ancestor intersections for exterior
    if (puzzler->checkExteriorIntersections) {
        if (intersectNodeExterior(node, puzzler)) { // simple check

            // prepare complex check
            treeNode* exterior = getParent(topLevelAncestor);
            TENTATIVE3_setupExteriorBoundingBoxes(exterior, topLevelAncestor, node, puzzler);
            changedNode = handleIntersectionWithAncestor(exterior, node, 0, puzzler);
            /*
            // setupExteriorBoundingBoxes(exterior, topLevelAncestor, node, puzzler);
            TENTATIVE_setupExteriorBoundingBoxes(exterior, topLevelAncestor, node, puzzler);

            // complex check
            intersectionType it = intersectNodeNode(node, exterior);
            if (it != noIntersection) {
                changedNode = handleIntersectionWithAncestor(exterior, node, 0, puzzler);
            }
            */
        }
    }

    return changedNode;
}

