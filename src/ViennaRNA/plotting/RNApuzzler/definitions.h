#ifndef DEFINITIONS_H
#define DEFINITIONS_H

#define FANCY_PS 0

#define EXTERIOR_Y 100.0

#define epsilonRecognize 14 // font size
#define epsilonFix 19       // font size + 5 @ resolveIntersections.c
#define EPSILON_0 1.0
#define EPSILON_3 1e-3
#define EPSILON_7 1e-7
#define MIN_POSITIVE_ANGLE +0.0000000001
#define MIN_NEGATIVE_ANGLE -0.0000000001

#define _false     0x0000
#define _intersect 0x0001
#define _changed   0x0002



PRIVATE void bubblesort(const int           numValues,
                const double *const valuesLevel1,
                const double *const valuesLevel2,
                int *const          indices);


/**
 * @brief
 *      Given a circle's radius and a distance between two points on the circle
 *      this function calculates the angle between those points.
 *      Note that the resulting angle will always be smaller than or equal to 180°.
 *      If knowing the wanted angle being greater than 180° just subtract the result from 360°.
 * @param radius the circle's radius
 * @param distance the distance between two points on the circle
 * @return angle in degree
 */
PRIVATE double distanceToAngle(const double radius,
                       const double distance);


PRIVATE double angleToDistance(const double radius,
                       const double degreeAngle);


#include "definitions.c"

#endif
