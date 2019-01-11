#ifndef VECTOR_MATH_H
#define VECTOR_MATH_H

#define MATH_PI      3.141592653589793
#define MATH_PI_HALF (MATH_PI / 2.0)
#define MATH_TWO_PI  (2.0 * MATH_PI)
#define TO_DEGREE    (180.0 / MATH_PI)
#define TO_RAD       (MATH_PI / 180.0)

double toDegree(double angle);
double toRad(double angle);

/**
 * @brief vectorLength2D
 *      - calculates the length of a 2D vector
 * @param vector
 *      - double array[2] with vector coordinates x,y
 * @return
 *      - length of input vector
 */
double vectorLength2D(const double vector[2]);

/**
 * @brief vectorLength2DSquared
 *      - calculates the squared length of a 2D vector
 * @param vector
 *      - double array[2] with vector coordinates x,y
 * @return
 *      - squared length of input vector
 */
double vectorLength2DSquared(const double vector[2]);

/**
 * @brief scalarProduct2D
 *      - calculates the scalar product (or dot product) of two given 2D vectors
 *        v1 * v2 = ( v1.x * v2.x + v1.y * v2.y ) / ( |v1| * |v2| )
 * @param vector1
 *      - double array[2] with vector coordinates x,y
 * @param vector2
 *      - double array[2] with vector coordinates x,y
 * @return
 *      - scalar product of both given input vectors
 */
double scalarProduct2D(const double* vector1, const double* vector2);

/**
 * @brief normalize
 *      - normalizes the vector to length 1, same direction
 * @param vector
 *      - double array[2] with vector coordinates x,y
 */
void normalize(double* vector);

/**
 * @brief isToTheRightPointPoint determines if a point is left or right to a line.
 * @param lineStart[2]
 * @param lineEnd[2]
 * @param point[2]
 * @return 1 if point is right to the vector from lineStart to lineEnd, 0 otherwise
 */
short isToTheRightPointPoint(
    const double lineStart[2],
    const double lineEnd[2],
    const double point[2]
);

/**
 * @brief isToTheRightPointVector
 *      - same as isToTheRight but works with vectors.
 * @param lineStart
 * @param lineVector
 * @param point
 * @return 1 if point is right to the vector starting at lineStart, 0 otherwise
 */
short isToTheRightPointVector(
    const double lineStart[2],
    const double lineVector[2],
    const double point[2]
);

/**
 * @brief angleBetweenVectors2D
 *      - calculates the angle between both given vectors.
 * @param vector1
 *      - double array[2] with vector coordinates x,y
 * @param vector2
 *      - double array[2] with vector coordinates x,y
 * @return
 *      - angle between vector1 and vector2 in degree
 */
double angleBetweenVectors2D(const double vector1[2], const double vector2[2]);

/**
 * @brief anglePtPtPt2D
 *      - calculates angle defined by three points.
 * @param p2
 *      - some point
 * @param p2
 *      - center point of that angle
 * @param p3
 *      - some point
 * @return
 *      - angle defined by points in degree
 */
double anglePtPtPt2D(const double* p1, const double* p2, const double* p3);

/**
 * @brief normalizeAngle
 *      - Does nothing if angle is inside [0,2*PI] or [0,360°]. Otherwise transforms the angle into this interval.
 * @param angle
 * @param useDegree
 * @return angle in interval [0,2*PI] or [0,360°]
 */
double normalizeAngle(const double angle, short useDegree);

/**
 * @brief rotatePointAroundPoint
 *      - Performs a rotation of a point around a rotation center by a given angle.
 *        A positive angle results in a clockwise rotation while
 *        a negative angle results in counterclockwise rotation.
 *        Returns to rotated point.
 * @param point
 *      - double array[2] with point coordinates x,y
 * @param rotationCenter
 *      - double array[2] with point coordinates x,y
 * @param angle
 *      - angle in radiant format
 * @param ret
 *      - insert double array[2] here to act as return value
 */
void rotatePointAroundPoint(const double* point, const double* rotationCenter, const double angle, double* ret);

/**
 * @brief rotateVectorByAngle
 *      - Performs a rotation of a vector by a given angle.
 *        A positive angle results in a clockwise rotation while
 *        a negative angle results in counterclockwise rotation.
 *        Returns to rotated vector.
 * @param vector
 *      - double array[2] with vector coordinates x,y
 * @param angle
 *      - angle in radiant format
 * @param ret
 *      - insert double array[2] here to act as return value
 */
void rotateVectorByAngle(const double* vector, const double angle, double* ret);

/**
 * @brief translatePointByVector
 *      - Performs a translation of a point by a given vector.
 *        Returns to rotated vector.
 * @param point
 *      - double array[2] with point coordinates x,y
 * @param trans
 *      - translation vector as double array[2] with vector coordinates x,y
 * @param ret
 *      - insert double array[2] here to act as return value
 */
void translatePointByVector(const double* point, const double* trans, double* ret);

/**
 * @brief solveSquareEquation
 *      - Solves a square equation with ABC method. Writes results to solution1 and solution2
 *        Input: a*x² + b*x + c = 0
 *        Output: solution1, solution2 if they exist
 * @param a
 * @param b
 * @param c
 * @param sol1
 *      - solution1 pointer
 * @param sol2
 *      - solution2 pointer
 * @return
 *      - count of solutions (0, 1 or 2)
 */
short solveSquareEquation(const double a, const double b, const double c, double* sol1, double* sol2);

/**
 * @brief getCutPointsOfCircles
 *      - Calculates the common points (cut points) of two circles given by their center
 *        and radius. The resulting cut points are returned in variables ret1 and ret2 while
 *        the methods return value states how many cut points exist.
 * @param c1
 *      - center of circle1. coordinates given as double array[2].
 * @param r1
 *      - radius of circle1.
 * @param c2
 *      - center of circle2. coordinates given as double array[2].
 * @param r2
 *      - radius of circle2.
 * @param ret1
 *      - insert double array[2] here. values are set if there is at least 1 cut point.
 * @param ret2
 *      - insert double array[2] here. values are set if there are 2 cut points.
 * @return
 *      - number of cut points: 0, 1 or 2. -1 if circles match (infinite cut points)
 */
short getCutPointsOfCircles(const double* c1, const double r1, const double* c2, const double r2, double* ret1, double* ret2);

/**
 * @brief getCutPointsOfCircleAndLine
 *      - Calculates the common points (cut points) of a given circle and a given line.
 *        The resulting cut points are returned in the variables ret1 and ret2 while
 *        the function also returns the number of cut points.
 * @param center
 *      - center of the circle. coordinates given as double array[2].
 * @param radius
 *      - radius of the circle.
 * @param anchor
 *      - anchor point of the line. coordinates given as double array[2].
 * @param direction
 *      - direction vector of the line. coordinates given as double array[2].
 * @param ret1
 *      - insert double array[2] here. values are set if there is at least 1 cut point.
 * @param ret2
 *      - insert double array[2] here. values are set if there are 2 cut point.
 * @return
 *      - number of cut points: 0, 1 or 2.
 */
short getCutPointsOfCircleAndLine(const double* center, const double radius, const double* anchor, const double* direction, double* ret1, double* ret2);

/**
 * @brief vector
 *      - creates a vector from pStart to pEnd.
 * @param pStart
 *      - double array[2] with coordinates of point
 * @param pEnd
 *      - double array[2] with coordinates of point
 * @param v
 *      - double array[2] used to return the vector
 */
void vector(const double pStart[2], const double pEnd[2], double v[2]);

/**
 * @brief normal
 *      - creates the normal vector to the given input vector.
 *        The output normal vector has length 1.
 * @param v
 *      - double array[2] with vector coordinates
 * @param n
 *      - double array[2] used to return the normal vector
 */
void normal(const double v[2], double n[2]);

void unit(const double v[2], double u[2]);

/**
 * @brief min
 *      - the typical minimum function. calculates the minimum of two given numbers.
 * @param number1
 * @param number2
 * @return
 *      - the smaller number
 */
double min(const double number1, const double number2);

/**
 * @brief sign
 *      - calculates the sign of a given number.
 * @param number
 * @return
 *      - the sign of the number
 */
double sign(const double number);

void circle(const double A[2], const double B[2], const double C[2], double center[2], double* radius);

#endif
