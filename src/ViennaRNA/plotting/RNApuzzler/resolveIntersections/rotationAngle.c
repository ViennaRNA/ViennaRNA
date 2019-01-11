#include "ViennaRNA/RNApuzzler/resolveIntersections/rotationAngle.h"
#include "ViennaRNA/RNApuzzler/dataTypes/configtree_struct.h"
#include "ViennaRNA/RNApuzzler/data/boundingBoxes.h"
#include "ViennaRNA/RNApuzzler/data/configtree.h"
#include "ViennaRNA/RNApuzzler/vector_math.h"
#include "ViennaRNA/RNApuzzler/intersectLevel/intersectLevelBoundingBoxes.h"
#include "ViennaRNA/RNApuzzler/output/output.h"
#include "ViennaRNA/RNApuzzler/output/GeoGebra_output.h"
#include "ViennaRNA/RNApuzzler/resolveIntersections/intersectionType.h"

#include <stdio.h>
#include <string.h>
#include <math.h>

double pointToAngle(
    const double center[2],
    const double vRef[2],
    const short rotationSign,
    const double point[2]
) {
    char* fnName = "POINT TO ANGLE";
    double vCenterToPoint[2];
    vector(center, point, vCenterToPoint);
    double angle = angleBetweenVectors2D(vRef, vCenterToPoint);
    short cw = isToTheRightPointVector(center, vRef, point);

    if (rotationSign > 0 &&  cw) {
//        angle = angle;
    } else
    if (rotationSign > 0 && !cw) {
        angle = MATH_TWO_PI - angle;
    } else
    if (rotationSign < 0 &&  cw) {
        angle = -MATH_TWO_PI + angle;
    } else
    if (rotationSign < 0 && !cw) {
        angle = -angle;
    }
//    printf("[%s] angle: %3.2f° - sign=%3.2f C=(%3.2f, %3.2f) v=(%3.2f, %3.2f) P=(%3.2f, %3.2f)\n"
//           , fnName, angle, rotationSign
//           , center[0], center[1]
//           , vRef[0], vRef[1]
//           , point[0], point[1]
//           );
    return angle;
}

double fixIntersectionOfRectangleAndCircle(
        const double staticRectCenter[2],
        const double staticRectVecA[2],
        const double staticRectVecB[2],
        const double staticRectLengthA, // used for debug only
        const double staticRectLengthB,
        const double mobileCircCenter[2],
        const double mobileCircRadius,
        const double rotationCenter[2],
        const short rotationSign
) {
    char* fnName = "fixIntersectionOfCircleAndRectangle";

    /// emergency exit
    if (rotationSign == 0) {
        printError(fnName, "invalid rotation sign\n");
        return 0.0;
    }

    /// some extra distance after intersection resolution
    const double distance = epsilonFix + mobileCircRadius;

    /// circle definition (Center = centerNode->loop->center; sourceNode->loop->center on the circles periphery)
    double vRotationCenterToInPoint[2];
    vector(rotationCenter, mobileCircCenter, vRotationCenterToInPoint);
    double rotationRadius = vectorLength2D(vRotationCenterToInPoint);

    /// line definition
    double axisOffset = staticRectLengthB + distance;
    double axisDirection[2] = { staticRectVecA[0], staticRectVecA[1] };
    double axisAnchorPositive[2] = { staticRectCenter[0] + axisOffset * staticRectVecB[0],
                                     staticRectCenter[1] + axisOffset * staticRectVecB[1]
                                   };
    double axisAnchorNegative[2] = { staticRectCenter[0] - axisOffset * staticRectVecB[0],
                                     staticRectCenter[1] - axisOffset * staticRectVecB[1]
                                   };


    /// cut point computation
    int count = 0;
    double cut[4][2];

    short numCutPointsPositive = getCutPointsOfCircleAndLine(rotationCenter, rotationRadius,
                                                     axisAnchorPositive, axisDirection,
                                                     cut[count+0], cut[count+1]);
    count += numCutPointsPositive;

    short numCutPointsNegative = getCutPointsOfCircleAndLine(rotationCenter, rotationRadius,
                                                     axisAnchorNegative, axisDirection,
                                                     cut[count+0], cut[count+1]);
    count += numCutPointsNegative;

    if (count == 0) {
        /// in case we found no cut points
        /// (example szenario: RF00100 >AANU01122137.1_2893-3214 Macaca mulatta 1099214708557, whole genome shotgun sequence.)
        /// we just take those points on the circle that are closest to the lines
        double axisNormal[2];
        normal(axisDirection, axisNormal);
        cut[count][0] = rotationCenter[0] + rotationRadius * axisNormal[0];
        cut[count][1] = rotationCenter[1] + rotationRadius * axisNormal[1];
        count++;
        cut[count][0] = rotationCenter[0] - rotationRadius * axisNormal[0];
        cut[count][1] = rotationCenter[1] - rotationRadius * axisNormal[1];
        count++;

    }

    /// transform the calculated points into rotation angles
    double angles[4];
    for (int i = 0; i < count; i++) {
        angles[i] = pointToAngle(rotationCenter, vRotationCenterToInPoint, rotationSign, cut[i]);
    }

    // fix underflows
    for (int i = 0; i < count; i++) {
        if (angles[i] == 0.0) {
            angles[i] = signbit(angles[i]) ? MIN_NEGATIVE_ANGLE : MIN_POSITIVE_ANGLE;
        }
    }

    double angle = rotationSign * MATH_TWO_PI;
    for (int i = 0; i < count; i++) {
//        printDebug(fnName, "sign: %+2.0f | old: %7.2f° | angle: %7.2f°", rotationSign, angle, angles[i]);
        if (rotationSign > 0.0 && angles[i] > 0.0) {
            angle = fmin(angle, angles[i]);
        }
        if (rotationSign < 0.0 && angles[i] < 0.0) {
            angle = fmax(angle, angles[i]);
        }
//        printDebug(NULL, "| new: %7.2f°\n", angle);
    }

    if (fabs(angle) == 0.0 || fabs(angle) == MATH_TWO_PI) {
        angle = 0.0;
        printWarning(fnName, "no valid rotation here.\n");

        short showDetails = 0;
        if (showDetails) {
            // rectangle
            double staticRectExt[2] = { staticRectLengthA, staticRectLengthB };
            printError(fnName, ""); GEOGEBRA_printRectangle("rect", staticRectVecA, staticRectVecB, staticRectCenter, staticRectExt);
            // circle
            printError(fnName, ""); GEOGEBRA_printCircle("circ", rotationCenter, rotationRadius);
            // lines
            printError(fnName, ""); GEOGEBRA_printLinePointDir("line_{pos}", axisAnchorPositive, axisDirection);
            printError(fnName, ""); GEOGEBRA_printLinePointDir("line_{neg}", axisAnchorNegative, axisDirection);
            // cut points
            for (int i = 0; i < count; i++) {
                printError(fnName, "Cut%d = (%f, %f)\n", i, cut[i][0], cut[i][1]);
            }
            // angles
            for (int i = 0; i < count; i++) {
                printError(fnName, "angle%d = %f°\n", angles[i]);
            }
        }
    }
//    printDebug(fnName, "return: %+7.2f°\n", angle);
    return angle;

}

double fixIntersectionOfCircles(
        const double staticCircleCenter[2],
        const double staticCircleRadius,
        const double mobileCircleCenter[2],
        const double mobileCircleRadius,
        const double rotationCenter[2],
        const short rotationSign
) {
    char* fnName = "fixIntersectionOfCircles";

    /// emergency exit
    if (rotationSign == 0) {
        printError(fnName, "invalid rotation sign\n");
        return 0.0;
    }

    /// some extra distance after intersection resolution
    double distance = epsilonFix;

    /// circle around centerNode
    double vRotationCenterToCircleLoopCenter[2];
    vector(rotationCenter, mobileCircleCenter, vRotationCenterToCircleLoopCenter);
    double rotationRadius = vectorLength2D(vRotationCenterToCircleLoopCenter);

    /// circle around rootNode (extended by intersectorNode's radius and some extra distance)
    double extendedStaticCircleRadius = staticCircleRadius + mobileCircleRadius + distance;

    /// cut point computation
    double cut1[2], cut2[2];
    short numCutPoints = getCutPointsOfCircles(rotationCenter, rotationRadius,
                                               staticCircleCenter, extendedStaticCircleRadius,
                                               cut1, cut2);

    /// emergency exit ... should never happen
    if (numCutPoints == 0) {
        printError(fnName, "calculated cut points: 0 expected: 1+\n");
        printError(fnName, ""); GEOGEBRA_printCircle("static", staticCircleCenter, staticCircleRadius);
        printError(fnName, ""); GEOGEBRA_printCircle("mobile", mobileCircleCenter, mobileCircleRadius);
        printError(fnName, ""); GEOGEBRA_printPoint("rotationcenter", rotationCenter);
        printError(fnName, "distance = %f\n", distance);
        printError(fnName, "rotationsign = %d\n", rotationSign);
        return 0.0;
    }

    /// get rotation angles from cut points
    double angle1 = 0.0;
    double angle2 = 0.0;
    {
        // get angle1
        double vCircleCenterToCut1[2];
        vector(rotationCenter, cut1, vCircleCenterToCut1);
        angle1 = angleBetweenVectors2D(vRotationCenterToCircleLoopCenter, vCircleCenterToCut1);
        short isCW1 = isToTheRightPointVector(rotationCenter, vRotationCenterToCircleLoopCenter, cut1);
        if (!isCW1) {
            angle1 *= -1;
        }
        // fix underflow
        if (angle1 == 0.0) {
            angle1 = signbit(angle1) ? MIN_NEGATIVE_ANGLE : MIN_POSITIVE_ANGLE;
        }

        // get angle2
        double vCircleCenterToCut2[2];
        vector(rotationCenter, cut2, vCircleCenterToCut2);
        angle2 = angleBetweenVectors2D(vRotationCenterToCircleLoopCenter, vCircleCenterToCut2);
        short isCW2 = isToTheRightPointVector(rotationCenter, vRotationCenterToCircleLoopCenter, cut2);
        if (!isCW2) {
            angle2 *= -1;
        }
        // fix underflow
        if (angle2 == 0.0) {
            angle2 = signbit(angle2) ? MIN_NEGATIVE_ANGLE : MIN_POSITIVE_ANGLE;
        }

        if (isCW1 == isCW2) {
            if (fabs(angle1) < fabs(angle2)) {
                // evaluate angle2 from the other side
                if (isCW2) {
                    angle2 = angle2 - MATH_TWO_PI;
                } else {
                    angle2 = MATH_TWO_PI - angle2;
                }
            } else {
                // evaluate angle1 from the other side
                if (isCW1) {
                    angle1 = angle1 - MATH_TWO_PI;
                } else {
                    angle1 = MATH_TWO_PI - angle1;
                }
            }
        }
    }

    double rotationAngle = 0.0;

    if (rotationSign == 1) {
        rotationAngle = fmax(angle1, angle2);
    } else if (rotationSign == -1) {
        rotationAngle = fmin(angle1, angle2);
    }

    if (rotationAngle == 0.0) {
        printWarning(fnName, "no valid rotation here.\n");
    }
    return rotationAngle;
}

/*----------------------------------------------------------------------*/

double getRotationAngleLxL(
        const treeNode* ancestor,
        const treeNode* rotationNode,
        const treeNode* intersector,
        const short rotationSign
) {
    char* fnName = "getRotationAngleLxL";

    loopBox* staticLoop = ancestor->lBox;
    loopBox* rotationLoop = rotationNode->lBox;
    loopBox* mobileLoop = intersector->lBox;

    double staticCircleCenter[2];
    getLBoxCenter(staticLoop, staticCircleCenter);
    double staticCircleRadius = staticLoop->r;

    double mobileCircleCenter[2];
    getLBoxCenter(mobileLoop, mobileCircleCenter);
    double mobileCircleRadius = mobileLoop->r;

    double rotationCenter[2];
    getLBoxCenter(rotationLoop, rotationCenter);

    double rotationAngle = fixIntersectionOfCircles(staticCircleCenter, staticCircleRadius, mobileCircleCenter, mobileCircleRadius, rotationCenter, rotationSign);
    if (rotationAngle == 0.0) {
        printWarning(fnName, "[%c %c %c] (promoted)\n", getNodeName(ancestor), getNodeName(rotationNode), getNodeName(intersector));
    }
    return rotationAngle;
}

double getRotationAngleLxS(
        const treeNode* ancestor,
        const treeNode* rotationNode,
        const treeNode* intersector,
        const short rotationSign
) {
    char* fnName = "getRotationAngleLxS";

    stemBox* staticRect = intersector->sBox;
    loopBox* mobileCirc = ancestor->lBox;
    loopBox* rotationLoop = rotationNode->lBox;

    short inverseRotationSign = (-1) * rotationSign;
    double inverseRotationAngle = fixIntersectionOfRectangleAndCircle(staticRect->c, staticRect->a, staticRect->b, staticRect->e[0], staticRect->e[1], mobileCirc->c, mobileCirc->r, rotationLoop->c, inverseRotationSign);
    double rotationAngle = (-1) * inverseRotationAngle;
    if (rotationAngle == 0.0) {
        printWarning(fnName, "[%c %c %c] (promoted)\n", getNodeName(ancestor), getNodeName(rotationNode), getNodeName(intersector));
    }
    return rotationAngle;
}

double getRotationAngleSxL(
        const treeNode* ancestor,
        const treeNode* rotationNode,
        const treeNode* intersector,
        const short rotationSign
) {
    char* fnName = "getRotationAngleSxL";

    stemBox* staticRect = ancestor->sBox;
    loopBox* mobileCirc = intersector->lBox;
    loopBox* rotationLoop = rotationNode->lBox;

    double rotationAngle = fixIntersectionOfRectangleAndCircle(staticRect->c, staticRect->a, staticRect->b, staticRect->e[0], staticRect->e[1], mobileCirc->c, mobileCirc->r, rotationLoop->c, rotationSign);
    if (rotationAngle == 0.0) {
        printWarning(fnName, "[%c %c %c] (promoted)\n", getNodeName(ancestor), getNodeName(rotationNode), getNodeName(intersector));
    }
    return rotationAngle;
}

double getRotationAngleLxB(
        const treeNode* ancestor,
        const treeNode* rotationNode,
        const treeNode* intersector,
        const short rotationSign
) {
    char* fnName = "getRotationAngleLxB";

    /// idea: construct circles around the intersecting loop and bulge and resolve their intersection

    loopBox* staticLoop = ancestor->lBox;
    stemBox* mobileStem = intersector->sBox;

    // ### static circle
    // --- grab the intersector loop as static circle
    double staticCircleCenter[2];
    getLBoxCenter(staticLoop, staticCircleCenter);
    double staticCircleRadius = staticLoop->r;

    // ### mobile circle
    // --- get bulge indices
    int mobileBulgeIndex = -1;
    short intersect = intersectLoopBulges(staticLoop, mobileStem, &mobileBulgeIndex);

    // --- define mobile circle from mobile bulge
    double mobileBulge[3][2];
    getBulgeCoordinates(mobileStem, mobileBulgeIndex, mobileBulge[0], mobileBulge[1], mobileBulge[2]);

    double mobileCircleCenter[2];
    double mobileCircleRadius = 1.0;
    circle(mobileBulge[0], mobileBulge[1], mobileBulge[2], mobileCircleCenter, &mobileCircleRadius);

    // ### rotation center
    // --- define rotation center from rotation loop
    loopBox* rotationLoop = rotationNode->lBox;
    double rotationCenter[2];
    getLBoxCenter(rotationLoop, rotationCenter);

    // ### resolve
    // --- fix intersection of circles
    double rotationAngle = fixIntersectionOfCircles(staticCircleCenter, staticCircleRadius, mobileCircleCenter, mobileCircleRadius, rotationCenter, rotationSign);
    if (rotationAngle == 0.0) {
        printWarning(fnName, "[%c %c %c] (promoted)\n", getNodeName(ancestor), getNodeName(rotationNode), getNodeName(intersector));
    }
    return rotationAngle;
}

double getRotationAngleBxL(
        const treeNode* ancestor,
        const treeNode* rotationNode,
        const treeNode* intersector,
        const short rotationSign
) {
    char* fnName = "getRotationAngleBxL";

    /// idea: construct circles around the intersecting bulge and loop and resolve their intersection

    stemBox* staticStem = ancestor->sBox;
    loopBox* mobileLoop = intersector->lBox;

    // ### static circle
    // --- get bulge indices
    int staticBulgeIndex = -1;
    short intersect = intersectLoopBulges(mobileLoop, staticStem, &staticBulgeIndex);

    // --- define static circle from static bulge
    double staticBulge[3][2];
    getBulgeCoordinates(staticStem, staticBulgeIndex, staticBulge[0], staticBulge[1], staticBulge[2]);

    double staticCircleCenter[2];
    double staticCircleRadius = 1.0;
    circle(staticBulge[0], staticBulge[1], staticBulge[2], staticCircleCenter, &staticCircleRadius);

    // ### mobile circle
    // --- grab the intersector loop as mobile circle
    double mobileCircleCenter[2];
    getLBoxCenter(mobileLoop, mobileCircleCenter);
    double mobileCircleRadius = mobileLoop->r;

    // ### rotation center
    // --- define rotation center from rotation loop
    loopBox* rotationLoop = rotationNode->lBox;
    double rotationCenter[2];
    getLBoxCenter(rotationLoop, rotationCenter);

    // ### resolve
    // --- fix intersection of circles
    double rotationAngle = fixIntersectionOfCircles(staticCircleCenter, staticCircleRadius, mobileCircleCenter, mobileCircleRadius, rotationCenter, rotationSign);
    if (rotationAngle == 0.0) {
        printWarning(fnName, "[%c %c %c] (promoted)\n", getNodeName(ancestor), getNodeName(rotationNode), getNodeName(intersector));
    }
    return rotationAngle;
}

double getRotationAngleSxS(
        const treeNode* ancestor,
        const treeNode* rotationNode,
        const treeNode* intersector,
        const short rotatioSign
) {
    char *fnName = "getRotationAngleSxS";

    double rotationAngle = getRotationAngleSxL(ancestor, rotationNode, intersector, rotatioSign);
    if (rotationAngle == 0.0) {
        printWarning(fnName, "[%c %c %c] (promoted)\n", getNodeName(ancestor), getNodeName(rotationNode), getNodeName(intersector));
    }
    return rotationAngle;
}

double getRotationAngleSxB(
        const treeNode* ancestor,
        const treeNode* rotationNode,
        const treeNode* intersector,
        const short rotationSign
) {
    char* fnName = "getRotationAngleSxB";

    stemBox* staticStem = ancestor->sBox;
    stemBox* mobileStem = intersector->sBox;
    loopBox* rotationLoop = rotationNode->lBox;

    int mobileBulgeIndex;
    short intersect = intersectStemBulges(staticStem, mobileStem, &mobileBulgeIndex);

    // ### mobile circle
    // --- define mobile circle from mobile bulge
    double mobileBulge[3][2];
    getBulgeCoordinates(mobileStem, mobileBulgeIndex, mobileBulge[0], mobileBulge[1], mobileBulge[2]);

    double mobileCircleCenter[2];
    double mobileCircleRadius = 1.0;
    circle(mobileBulge[0], mobileBulge[1], mobileBulge[2], mobileCircleCenter, &mobileCircleRadius);

    double rotationAngle = fixIntersectionOfRectangleAndCircle(staticStem->c, staticStem->a, staticStem->b, staticStem->e[0], staticStem->e[1], mobileCircleCenter, mobileCircleRadius, rotationLoop->c, rotationSign);
    if (rotationAngle == 0.0) {
        printWarning(fnName, "[%c %c %c] (promoted)\n", getNodeName(ancestor), getNodeName(rotationNode), getNodeName(intersector));
    }
    return rotationAngle;
}

double getRotationAngleBxS(
        const treeNode* ancestor,
        const treeNode* rotationNode,
        const treeNode* intersector,
        const short rotationSign
) {
    char* fnName = "getRotationAngleBxS";

    stemBox* staticStem = intersector->sBox;
    stemBox* mobileStem = ancestor->sBox;
    loopBox* rotationLoop = rotationNode->lBox;

    int mobileBulgeIndex;
    short intersect = intersectStemBulges(staticStem, mobileStem, &mobileBulgeIndex);

    // ### mobile circle
    // --- define mobile circle from mobile bulge
    double mobileBulge[3][2];
    getBulgeCoordinates(mobileStem, mobileBulgeIndex, mobileBulge[0], mobileBulge[1], mobileBulge[2]);

    double mobileCircleCenter[2];
    double mobileCircleRadius = 1.0;
    circle(mobileBulge[0], mobileBulge[1], mobileBulge[2], mobileCircleCenter, &mobileCircleRadius);

    double rotationAngle = fixIntersectionOfRectangleAndCircle(staticStem->c, staticStem->a, staticStem->b, staticStem->e[0], staticStem->e[1], mobileCircleCenter, mobileCircleRadius, rotationLoop->c, rotationSign);
    if (rotationAngle == 0.0) {
        printWarning(fnName, "[%c %c %c] (promoted)\n", getNodeName(ancestor), getNodeName(rotationNode), getNodeName(intersector));
    }
    return rotationAngle;
}

double getRotationAngleBxB(
        const treeNode* ancestor,
        const treeNode* rotationNode,
        const treeNode* intersector,
        const short rotationSign
) {
    char* fnName = "getRotationAngleBxB";

    /// idea: construct circles around both bulges and resolve their intersection

    stemBox* staticStem = ancestor->sBox;
    stemBox* mobileStem = intersector->sBox;

    // --- get bulge indices
    int staticBulgeIndex = -1;
    int mobileBulgeIndex = -1;
    short intersect = intersectBulgesBulges(staticStem, mobileStem, &staticBulgeIndex, &mobileBulgeIndex);

    // ### static circle
    // --- define static circle from static bulge
    double staticBulge[3][2];
    getBulgeCoordinates(staticStem, staticBulgeIndex, staticBulge[0], staticBulge[1], staticBulge[2]);

    double staticCircleCenter[2];
    double staticCircleRadius = 1.0;
    circle(staticBulge[0], staticBulge[1], staticBulge[2], staticCircleCenter, &staticCircleRadius);

    // ### mobile circle
    // --- define mobile circle from mobile bulge
    double mobileBulge[3][2];
    getBulgeCoordinates(mobileStem, mobileBulgeIndex, mobileBulge[0], mobileBulge[1], mobileBulge[2]);

    double mobileCircleCenter[2];
    double mobileCircleRadius = 1.0;
    circle(mobileBulge[0], mobileBulge[1], mobileBulge[2], mobileCircleCenter, &mobileCircleRadius);

    // ### rotation center
    // --- define rotation center from rotation loop
    loopBox* rotationLoop = rotationNode->lBox;
    double rotationCenter[2];
    getLBoxCenter(rotationLoop, rotationCenter);

    // ### resolve
    // --- fix intersection of circles
    double rotationAngle = fixIntersectionOfCircles(staticCircleCenter, staticCircleRadius, mobileCircleCenter, mobileCircleRadius, rotationCenter, rotationSign);
    if (rotationAngle == 0.0) {
        printWarning(fnName, "[%c %c %c] (promoted)\n", getNodeName(ancestor), getNodeName(rotationNode), getNodeName(intersector));
    }
    return rotationAngle;
}

double getRotationAngle(
        const treeNode* rootNode,
        const treeNode* centerNode,
        const treeNode* intersectorNode,
        const intersectionType it,
        short rotationSign
) {
    char* fnName = "getRotationAngle";
    /// performs the appropriate calculation method for the given intersection type

    double rotationAngle = 0.0;

    switch (it) {
        case LxL:
            rotationAngle = getRotationAngleLxL(rootNode, centerNode, intersectorNode, rotationSign);
            break;

        case LxS:
            rotationAngle = getRotationAngleLxS(rootNode, centerNode, intersectorNode, rotationSign);
            break;

        case LxB:
            rotationAngle = getRotationAngleLxB(rootNode, centerNode, intersectorNode, rotationSign);
            break;

        case SxL:
            rotationAngle = getRotationAngleSxL(rootNode, centerNode, intersectorNode, rotationSign);
            break;

        case SxS:
            rotationAngle = getRotationAngleSxS(rootNode, centerNode, intersectorNode, rotationSign);
            break;

        case SxB:
            rotationAngle = getRotationAngleSxB(rootNode, centerNode, intersectorNode, rotationSign);
            break;

        case BxL:
            rotationAngle = getRotationAngleBxL(rootNode, centerNode, intersectorNode, rotationSign);
            break;

        case BxS:
            rotationAngle = getRotationAngleBxS(rootNode, centerNode, intersectorNode, rotationSign);
            break;

        case BxB:
            rotationAngle = getRotationAngleBxB(rootNode, centerNode, intersectorNode, rotationSign);
            break;

        default:
            printf(fnName, "no computation for given intersection type\n");
    }

    return rotationAngle;
}
