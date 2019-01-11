#include "ViennaRNA/RNApuzzler/resolveIntersections/calcDeltas.h"
#include "ViennaRNA/RNApuzzler/resolveIntersections/boundingWedge.h"
#include "ViennaRNA/RNApuzzler/dataTypes/configtree_struct.h"
#include "ViennaRNA/RNApuzzler/data/configtree.h"
#include "ViennaRNA/RNApuzzler/data/cfg_reader.h"
#include "ViennaRNA/RNApuzzler/vector_math.h"
#include "ViennaRNA/RNApuzzler/output/output.h"

#include "ViennaRNA/utils.h"

#include <stdlib.h>
#include <math.h>

/*
double DEPRECATED_getMinOuterAngle(
    const treeNode *node
) {
    stemBox* sBox = node->sBox;
    double pStemTopCorner[2];
    pStemTopCorner[0] = sBox->c[0] + sBox->e[0] * sBox->a[0] + sBox->e[1] * sBox->b[0];
    pStemTopCorner[1] = sBox->c[1] + sBox->e[0] * sBox->a[1] + sBox->e[1] * sBox->b[1];
    double pLoopCenter[2];
    getLoopCenter(node, pLoopCenter);
    double vLoopCenterToStemTopCorner[2];
    vector(pLoopCenter, pStemTopCorner, vLoopCenterToStemTopCorner);
    double vLoopCenterToStemCenter[2] = { (-1) * sBox->a[0], (-1) * sBox->a[1] };
    double minOuterAngle = angleBetweenVectors2D(vLoopCenterToStemCenter, vLoopCenterToStemTopCorner);

    return minOuterAngle;
}
*/

void calcDeltasEquidistantIncrease(
        const double targetAngleIn,
        const int configSize,
        short* increase,
        double* deltaCfg
) {
    char* fnName = "CALC DELTAS EQUIDISTANT INCREASE";
    double targetAngle = targetAngleIn;

    int increaseCount = 0;
    for (int i = 0; i < configSize; i++) {
        if (increase[i]) {
            increaseCount++;
        }
    }
    double deltaPerIncrease = targetAngle / increaseCount;

    for (int i = 0; i < configSize; i++) {
//        printDebug(fnName, "deltaCfg[%d]: %5.2lf°", i, deltaCfg[i]);
        if (increase[i]) {
            deltaCfg[i] += deltaPerIncrease;
        }
//        printDebug(NULL, " -> %5.2lf°\n", deltaCfg[i]);
    }
}

double calcDeltasMaximumFirstDecrease(
        const double targetAngleIn,
        const int indexLeft,
        const int indexRight,
        const int configSize,
        double* deltaCfg,
        double* currentAngles,
        const double minAngleHalf
) {
    char* fnName = "CALC DELTAS MAXIMUM FIRST DECREASE";
    double targetAngle = targetAngleIn;
//    printf("[%s] iLeft: %d iRight: %d\n", fnName, indexLeft , indexRight);

    int i;

    short doLoop = 1;

    while (doLoop) {
        double maxSpace = 0.0;
        int maxSpaceIndex = -1;

        if (indexLeft == -1) {
//            printf("[%s] behavior: iterate right\n", fnName);

            double sumAngles = 0.0;
            i = -1;
            while (i != indexRight) {
                i++;
                double cfg = currentAngles[i] + deltaCfg[i] - 2 * minAngleHalf;
                sumAngles += cfg;
            }
            while (i != configSize-1) {
                i++;
                double cfg = currentAngles[i] + deltaCfg[i] - 2 * minAngleHalf;
                if (sumAngles < MATH_PI) {
                    if (cfg > maxSpace) {
                        maxSpace = cfg;
                        maxSpaceIndex = i;
                    }
                } else {
                    break;
                }
                /// sum increase happens afterwards to allow bending the arc containing 180°
                sumAngles += cfg;
            }

        } else
        if (indexRight == -1) {
//            printf("[%s] behavior: iterate left\n", fnName);

            double sumAngles = 0.0;
            i = configSize - 1;
            while (i != indexLeft) {
                double cfg = currentAngles[i] + deltaCfg[i] - 2 * minAngleHalf;
                sumAngles += cfg;
                i--;
            }
            while (i != -1) {
//                printf("[%s] i: %d size: %d\n", fnName, i, configSize);
                double cfg = currentAngles[i] + deltaCfg[i] - 2 * minAngleHalf;
                if (sumAngles < MATH_PI) {
                    if (cfg > maxSpace) {
                        maxSpace = cfg;
                        maxSpaceIndex = i;
                    }
                } else {
                    break;
                }
                /// sum increase happens afterwards to allow bending the arc containing 180°
                sumAngles += cfg;
                i--;
            }

        } else {
            // default behavior
//            printf("[%s] behavior: default\n", fnName);

            i = indexRight;
            if (i == configSize-1) { i = -1; }
            while (i != indexLeft) {
//                printf("[%s] i: %d size: %d\n", fnName, i+1, configSize);
                double cfg = currentAngles[i+1] + deltaCfg[i+1] - 2 * minAngleHalf;
                if (cfg > maxSpace) {
                    maxSpace = cfg;
                    maxSpaceIndex = i+1;
//                    printf("[%s] newMax: %+7.2f°[%d] (cfg: %+7.2f° delta: %+7.2f° space: %+7.2f°)\n", fnName, maxSpace, maxSpaceIndex, toDegree(currentAngles[maxSpaceIndex]), deltaCfg[maxSpaceIndex], (-2) * minAngleHalf);
                }
                i++;
                if (i == configSize-1) { i = -1; }
            }
        }

        /// using spaces
//        for (i = 0; i < configSize; i++) {
//            if (decrease[i]) {
//                if (currentAngles[i] > maxSpace) {
//                    maxSpace = currentAngles[i];
//                    maxSpaceIndex = i;
//                }
//            }
//        }

//        maxSpace = toDegree(maxSpace);

        double diff = 0.0;
        if (maxSpaceIndex != -1) {
            double factor = (targetAngle < 0.1 * targetAngleIn) ? 1.0 : 0.5;
            diff = (-1) * fmin(factor * maxSpace, targetAngle);
//            printf("[%s] diff: %+7.2f = (-1) * fmin(%+7.2f, %+7.2f)\n", fnName, diff, (factor * maxSpace), targetAngle);
//            printf("[%s] cfg[%d]: %+7.2f° (+ %+7.2f°) diff: %+7.2f°\n", fnName, maxSpaceIndex, toDegree(currentAngles[maxSpaceIndex]), deltaCfg[maxSpaceIndex], diff);
            deltaCfg[maxSpaceIndex] += diff;
            targetAngle += diff;
        }

        doLoop = targetAngle > 0.0 && fabs(diff) > epsilon3;
    }



    return targetAngle;
}

double calcDeltasNearestNeighborsFirstDecrease(
        const double targetAngleIn,
        const int indexLeft,
        const int indexRight,
        const int configSize,
        short* decrease,
        double* space,
        double* deltaCfg
) {
    char* fnName = "CALC DELTAS NEAREST NEIGHBOR FIRST DECREASE";
    double targetAngle = targetAngleIn;

    /// count the number of possible iteration steps
    int startIndex = indexRight + 1;
    if (startIndex == configSize) { startIndex = -1; }
    int stopIndex = indexLeft + 1;
    if (stopIndex == configSize) { stopIndex = -1; }
//    printf("[%s] start: %d stop: %d\n", fnName, startIndex, stopIndex);
    int steps = 0;
    int stemIt = indexRight;
    while (stemIt != indexLeft) {
        stemIt++;
        if (stemIt == configSize) { stemIt = -1; }

        steps++;
//        printf("[%s] stemIt: %d steps: %d\n", fnName, stemIt, steps);
    }
    int numIt = steps / 2; // implicit floor() operation
//    printf("[%s] indexL: %d indexR: %d steps: %d numItFloat: %f numItInt: %d\n", fnName, indexLeft, indexRight, steps, 0.5 * steps, numIt);

    int* index = (int*) vrna_alloc(steps * sizeof(int));
    short changed = 1;
    while (changed) {
        changed = 0;
        int count = 0;

        int iL = indexLeft;
        if (iL == -1) { iL = configSize - 1; }
        int iR = indexRight + 1;
        if (iR == configSize) { iR = 0; }

        for (int i = 0; i < numIt; i++) {
            if (decrease[iL]) {
                index[count] = iL;
                count++;
            }
            if (decrease[iR]) {
                index[count] = iR;
                count++;
            }
            iL--;
            if (iL == -1) { iL = configSize - 1; }
            iR++;
            if (iR == configSize) { iR = 0; }
        }
        if (numIt < 0.5 * steps) {
            index[count] = iL;
            count++;
            iL--;
            if (iL == -1) { iL = configSize - 1; }
        }

//        printf("[%s] index queue:", fnName);
//        for (i = 0; i < count; i++) {
//            printf(" %d", index[i]);
//        }
//        printf("\n");

        if (count > 0) {
            double partAngle = targetAngle / count;
            for (int k = 0; k < count; k++) {
                int j = index[k];
                if (decrease[j]) {
                    double diff = (-1) * fmin(space[j] + deltaCfg[j], partAngle);
                    deltaCfg[j] += diff;
                    targetAngle += diff;
                    changed = changed || (diff != 0.0);
                }
            }
        }

    }
    free(index);

//    printf("[%s] return: %+7.2f°\n", fnName, targetAngle);

    return targetAngle;
}

/**
 * @brief calcDeltas
 *      The area between stems indexLeft and indexRight
 *      (by traversing the loop clockwise starting at indexLeft-stem)
 *      will be enlarged in degree as given via deltaAngle.
 *      All other areas will be used to compensate that increase
 *      (i.e. by decreasing those area's angles).
 * @param node
 * @param recursiveEnd
 * @param indexLeft
 * @param indexRight
 * @param deltaAngle
 * @param deltas
 * @return the amount of change (in positive degree) that can be accomplished with calculated deltas
 */
double calcDeltas(
        const treeNode* node,
        const treeNode* recursiveEnd,
        const int indexLeft,
        const int indexRight,
        const double deltaAngle,
        puzzlerOptions* puzzler,
        double* deltas
) {
    char* fnName = "CALC DELTAS";

    /// Check: valid angle >= 0.0
    if (deltaAngle < 0.0) {
        printError(fnName, "cannot handle negative angles! grant proper input! (deltaAngle: %+7.2f°)\n", deltaAngle);
        return 0.0;
    }

    /// Check: valid range
    if (indexLeft == indexRight) {
        printError(fnName, "non-sense input. indices have to be different. (left: %d right%d)\n", indexLeft, indexRight);
    }

    int childCount = node->childCount;
    int configSize = childCount + 1;

    /// get the current node's stem's bounding wedge
    //double minOuterAngle = getMinOuterAngle(node);
    //double minOuterAngle = 0.5 * getPairedAngle(node);
    double minOuterAngle = asin(puzzler->paired / (2 * node->cfg->radius));
//    printDebug(fnName, "%d minOuterAngle:%f\n", getNodeID(node), minOuterAngle);

    /// allocate memory for stuff used in calculation
    double* anglesMin     = (double*) vrna_alloc(childCount * sizeof (double));
    double* anglesMax     = (double*) vrna_alloc(childCount * sizeof (double));
    double* space         = (double*) vrna_alloc(configSize * sizeof (double));
    double* deltaCfg      = (double*) vrna_alloc(configSize * sizeof (double));
    short*  increase      = (short*)  vrna_alloc(configSize * sizeof (short));
    short*  decrease      = (short*)  vrna_alloc(configSize * sizeof (short));
    double* currentAngles = (double*) vrna_alloc(configSize * sizeof (double));

    /// Initialization currentAngles
    config *cfg = node->cfg;
    for (int currentArc = 0; currentArc < cfg->numberOfArcs; ++currentArc) {
        currentAngles[currentArc] = getArcAngle(cfg, currentArc);
    }

    /// get all bounding wedges (minAngle, maxAngle)
    double min, max;
    for (int currentChild = 0; currentChild < childCount; currentChild++) {
        getBoundingWedge(node, currentChild, &min, &max);
//        printDebug(fnName, "wedge[%d]  min:%5.2lf  max:%5.2lf\n",
//                   currentChild, min, max);

        anglesMin[currentChild] = min;
        anglesMax[currentChild] = max;
    }

    /// convert bounding wedges to "free" areas that can be used for compensation of changes
    space[0] = anglesMin[0] - (0 + minOuterAngle);
//    printDebug(fnName, "space[%d] = %+11.6lf : %+11.6lf - %+11.6lf\n", 0, space[0], anglesMin[0], minOuterAngle);
    for (int i = 1; i < (configSize - 1); i++) {
        space[i] = anglesMin[i] - anglesMax[i-1];
//        printDebug(fnName, "space[%d] = %+11.6lf : %+11.6lf - %+11.6lf\n", i, space[i], anglesMin[i], anglesMax[i]);
    }
    space[configSize - 1] = (MATH_TWO_PI - minOuterAngle) - anglesMax[configSize - 2];
//    printDebug(fnName, "space[%d] = %+11.6lf : %+11.6lf - %+11.6lf\n", configSize - 1, space[configSize - 1], (MATH_TWO_PI - minOuterAngle), anglesMax[configSize - 2]);

    // fix too big spaces (may become bigger than config for very large loops)
    for (int i = 0; i < configSize; i++) {
        space[i] = fmin(space[i], getArcAngle(node->cfg, i) - 2 * minOuterAngle);
    }

//    // debug
//    for (int i = 0; i < configSize; i++) {
//        printDebug(fnName, "space[%d] = %+11.6lf\n", i, space[i]);
//    }

    /// Initialization: calculation values (deltaCfg, increase, decrease)
    for (int i = 0; i < configSize; i++) {
        deltaCfg[i] = 0.0;
        increase[i] = -1;
        decrease[i] = -1;
    }

//    printf("[%s] iLeft: %d iRight: %d angle: %+7.2f°\n", fnName, indexLeft, indexRight, deltaAngle);

    /// Mark increase and decrease areas
    int currentIndex = indexLeft; // stemIndex
    while (currentIndex != indexRight) {
//        printf("[%s] (1) currentIndex: %d size: %d\n", fnName, currentIndex+1, configSize);
        increase[currentIndex+1] = 1;
        decrease[currentIndex+1] = 0;
        currentIndex++;
        if (currentIndex == configSize-1) { currentIndex = -1; }
    }
    while (currentIndex != indexLeft) {
//        printf("[%s] (2) currentIndex: %d size: %d\n", fnName, currentIndex+1, configSize);
        increase[currentIndex+1] = 0;
        decrease[currentIndex+1] = (space[currentIndex+1] > 0.0);
        currentIndex++;
        if (currentIndex == configSize-1) { currentIndex = -1; }
    }

    /// ------------------------
    /// --- ^              ^ ---
    /// --- | preparations | ---
    /// --- |              | ---
    /// ------------------------

    /// ------------------------
    /// --- |              | ---
    /// --- | calculation  | ---
    /// --- v              v ---
    /// ------------------------

//    printf("[%s] [SIGNS] %d", fnName, getNodeID(node));
//    for (int i = 0; i < configSize; i++) {
//        if (increase[i]) {
//            printf(" %d:+", i);
//        }
//        if (decrease[i]) {
//            printf(" %d:-", i);
//        }
//        if (!increase[i] && !decrease[i]) {
//            printf(" %d:#", i);
//        }
//    }
//    printf("\n");

    double targetAngle = deltaAngle;

    /// Step 1: equidistant increase
    calcDeltasEquidistantIncrease(targetAngle, configSize, increase, deltaCfg);

    /// Step 2: nearest neighbor first decrease
    targetAngle = calcDeltasNearestNeighborsFirstDecrease(targetAngle, indexLeft, indexRight, configSize, decrease, space, deltaCfg);

    /// Step 3: check if intersections are fixed
    short notFixedYet = (targetAngle != 0.0);
    if (notFixedYet) {
        /// if the intersection is not yet fixed

//        printDebug(fnName, "remaining target: %+7.2f°\n", targetAngle);

        /// check if there is a loop on a higher level that can be bend instead of this one
        /// if this is the case we can apply methods using spaces
        /// otherwise we need to use drastical measures (e.g. maximumFirstDecrease using cfg instead of spaces)

        treeNode* parent = getParent(node);
        short canGoHigher = 0;
        while (parent != recursiveEnd && !isExterior(parent)) {
//            printf("[%s] parent: %d\n", fnName, getNodeID(parent));
            short parentIsMultiLoop = isMultiLoop(parent);
            if (parentIsMultiLoop) {
//                printf("[%s] is multi loop\n", fnName);
                canGoHigher = 1;
                break;
            } else {
                // parent is interior loop
                double childAngle = getArcAngle(parent->cfg, 0);
                if (fabs(childAngle - MATH_PI) < epsilon3) {
                    // no op
                } else
                if (childAngle > MATH_PI) {
                    if (indexLeft == 0) {
//                        printf("[%s] can bend left\n", fnName);
                        canGoHigher = 1;
                        break;
                    }
                } else
                if (childAngle < MATH_PI) {
                    if (indexLeft == -1) {
//                        printf("[%s] can bend right\n", fnName);
                        canGoHigher = 1;
                        break;
                    }
                }
            }

            /// if current parent node can not be adapted check its parent
            parent = getParent(parent);
        }

//        printf("[%s] canGoHigher: %d\n", fnName, canGoHigher);
        if (!canGoHigher) {
            targetAngle = calcDeltasMaximumFirstDecrease(targetAngle, indexLeft, indexRight, configSize, deltaCfg, currentAngles, minOuterAngle);
        }
    }

    /// Step 4: equidistant increase with negative remaining target angle
    calcDeltasEquidistantIncrease((-1) * targetAngle, configSize, increase, deltaCfg);

    if (!cfgIsValid(cfg, deltaCfg)) {
        printError(fnName, "Deltas invalid 3\n");
        // printConfigError(fnName, node, deltaCfg);
//    } else {
//        printDebug(fnName, "Deltas valid 3\n");
    }

    /// Fix deltas if changes are too small.
    /// This is necessary because sometimes the calculation results in micro changes.
    /// These micro changes go along with an increase of the loops radius which causes
    /// new problems as the changes being too small to get enough distance to the
    /// changed loop and the intersector being stuck in collision (again).
    ///
    /// multiplying by factor 2.0 we always get a resulting angle between 0.1° and 0.2°
    /// don't use factor 10 as the impact of doing so is way too strong and often causes crashes
    /// in term of applicability of the changes
    short fixTooSmallChanges = 0;
    if (fixTooSmallChanges) {
        for (int cntr = 0; cntr < 100; cntr++) {
            short valid = 0;
            for (int currentArc = 0; currentArc < configSize; currentArc++) {
                if (fabs(deltaCfg[currentArc]) >= epsilon3) {
                    valid = 1;
                    break;
                }
            }
            if (valid) {
                break;
            } else {
                for (int currentArc = 0; currentArc < configSize; currentArc++) {
                    deltaCfg[currentArc] = 2.0 * deltaCfg[currentArc];
                }
            }
        }
//        if (LOG_FLAG && cntr > 0) {
//            printf("[ LOG ] fixing... (%d)\n", cntr);
//        }
    }

    /// transfer calculated deltas to return area
    for (int currentArc = 0; currentArc < configSize; currentArc++) {
        deltas[currentArc] = deltaCfg[currentArc];
//        printf("[%s] delta[%d]: %+7.2f° (space[%d]: %+7.2f°)\n", fnName, currentArc, deltas[currentArc], currentArc, space[currentArc]);
    }

    /// free allocated memory
    free(anglesMin);
    free(anglesMax);
    free(space);
    free(deltaCfg);
    free(increase);
    free(decrease);
    free(currentAngles);

    /// check if all deltas sum up to zero
    double checkSum = 0.0;
    for (int currentArc = 0; currentArc < configSize; currentArc++) {
        checkSum += deltas[currentArc];
    }
    if (fabs(checkSum) > epsilon3) {
        printf("[%s] config broke ... abort and reset\n", fnName);
        for (int currentArc = 0; currentArc < configSize; currentArc++) {
//            printf("[%s] delta[%d]: %+7.2f° (space[%d]: %+7.2f°)\n", fnName, currentArc, deltas[currentArc], currentArc, space[currentArc]);
            deltas[currentArc] = 0.0;
        }
        targetAngle = deltaAngle;
    }

    if (!cfgIsValid(cfg, deltas)) {
        printConfigError(fnName, node, deltas);

        for (int currentArc = 0; currentArc < configSize; currentArc++) {
//            printf("[%s] delta[%d]: %+7.2f° (space[%d]: %+7.2f°)\n", fnName, currentArc, deltas[currentArc], currentArc, space[currentArc]);
            deltas[currentArc] = 0.0;
        }
        targetAngle = deltaAngle;
//    } else {
//        printDebug(fnName, "Deltas valid\n");
    }

//    for (int currentArc = 0; currentArc < configSize; ++currentArc) {
//        printDebug(fnName, "delta[%d] = %05.2f°\n", currentArc, deltas[currentArc]);
//    }

    /// return the difference that can be accomplished using these deltas
    double changedAngle = deltaAngle - targetAngle;
//    printf("[%s] %+7.2f° = %+7.2f° - %+7.2f°\n", fnName, changedAngle, deltaAngle, targetAngle);
    return changedAngle;
}

