#ifndef RNAPUZZLER_VECTOR_MATH_H
#define RNAPUZZLER_VECTOR_MATH_H

/*
 *      RNApuzzler math helper
 *
 *      c  Daniel Wiegreffe, Daniel Alexander, Dirk Zeckzer
 *      ViennaRNA package
 */

#include <math.h>
#include <stdio.h>

#include "definitions.inc"

#define MATH_PI      3.141592653589793
#define MATH_PI_HALF (MATH_PI / 2.0)
#define MATH_TWO_PI  (2.0 * MATH_PI)
#define TO_DEGREE    (180.0 / MATH_PI)
#define TO_RAD       (MATH_PI / 180.0)

PRIVATE double
toDegree(double angle);


PRIVATE double
toRad(double angle);


/**
 * @brief vectorLength2D
 *      - calculates the length of a 2D vector
 * @param vector
 *      - double array[2] with vector coordinates x,y
 * @return
 *      - length of input vector
 */
PRIVATE double
vectorLength2D(const double vector[2]);


/**
 * @brief vectorLength2DSquared
 *      - calculates the squared length of a 2D vector
 * @param vector
 *      - double array[2] with vector coordinates x,y
 * @return
 *      - squared length of input vector
 */
PRIVATE double
vectorLength2DSquared(const double vector[2]);


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
PRIVATE double
scalarProduct2D(const double  *vector1,
                const double  *vector2);


/**
 * @brief normalize
 *      - normalizes the vector to length 1, same direction
 * @param vector
 *      - double array[2] with vector coordinates x,y
 */
PRIVATE void
normalize(double *vector);


/**
 * @brief isToTheRightPointPoint determines if a point is left or right to a line.
 * @param lineStart[2]
 * @param lineEnd[2]
 * @param point[2]
 * @return 1 if point is right to the vector from lineStart to lineEnd, 0 otherwise
 */
PRIVATE short
isToTheRightPointPoint(const double lineStart[2],
                       const double lineEnd[2],
                       const double point[2]);


/**
 * @brief isToTheRightPointVector
 *      - same as isToTheRight but works with vectors.
 * @param lineStart
 * @param lineVector
 * @param point
 * @return 1 if point is right to the vector starting at lineStart, 0 otherwise
 */
PRIVATE short
isToTheRightPointVector(const double  lineStart[2],
                        const double  lineVector[2],
                        const double  point[2]);


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
PRIVATE double
angleBetweenVectors2D(const double  vector1[2],
                      const double  vector2[2]);


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
PRIVATE double
anglePtPtPt2D(const double  *p1,
              const double  *p2,
              const double  *p3);


/**
 * @brief normalizeAngle
 *      - Does nothing if angle is inside [0,2*PI] or [0,360°]. Otherwise transforms the angle into this interval.
 * @param angle
 * @param useDegree
 * @return angle in interval [0,2*PI] or [0,360°]
 */
PRIVATE double
normalizeAngle(const double angle,
               short        useDegree);


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
PRIVATE void
rotatePointAroundPoint(const double *point,
                       const double *rotationCenter,
                       const double angle,
                       double       *ret);


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
PRIVATE void
rotateVectorByAngle(const double  *vector,
                    const double  angle,
                    double        *ret);


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
PRIVATE void
translatePointByVector(const double *point,
                       const double *trans,
                       double       *ret);


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
PRIVATE short
solveSquareEquation(const double  a,
                    const double  b,
                    const double  c,
                    double        *sol1,
                    double        *sol2);


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
PRIVATE short
getCutPointsOfCircles(const double  *c1,
                      const double  r1,
                      const double  *c2,
                      const double  r2,
                      double        *ret1,
                      double        *ret2);


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
PRIVATE short
getCutPointsOfCircleAndLine(const double  *center,
                            const double  radius,
                            const double  *anchor,
                            const double  *direction,
                            double        *ret1,
                            double        *ret2);


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
PRIVATE void
vector(const double pStart[2],
       const double pEnd[2],
       double       v[2]);


/**
 * @brief normal
 *      - creates the normal vector to the given input vector.
 *        The output normal vector has length 1.
 * @param v
 *      - double array[2] with vector coordinates
 * @param n
 *      - double array[2] used to return the normal vector
 */
PRIVATE void
normal(const double v[2],
       double       n[2]);


PRIVATE void
unit(const double v[2],
     double       u[2]);


PRIVATE void
circle(const double A[2],
       const double B[2],
       const double C[2],
       double       center[2],
       double       *radius);


PRIVATE double
toDegree(double angle)
{
  return angle * TO_DEGREE;
}


PRIVATE double
toRad(double angle)
{
  return angle * TO_RAD;
}


PRIVATE double
vectorLength2D(const double vector[2])
{
  double  x = vector[0];
  double  y = vector[1];

  return sqrt(x * x + y * y);
}


PRIVATE double
vectorLength2DSquared(const double vector[2])
{
  double  x = vector[0];
  double  y = vector[1];

  return x * x + y * y;
}


PRIVATE double
scalarProduct2D(const double  *vector1,
                const double  *vector2)
{
  double  x1            = vector1[0];
  double  y1            = vector1[1];
  double  x2            = vector2[0];
  double  y2            = vector2[1];
  double  scalarProduct = (x1 * x2 + y1 * y2);

  return scalarProduct;
}


PRIVATE void
normalize(double *vector)
{
  double length = vectorLength2D(vector);

  vector[0] /= length;
  vector[1] /= length;
}


PRIVATE short
isToTheRightPointPoint(const double lineStart[2],
                       const double lineEnd[2],
                       const double point[2])
{
  /*
   * implicite knowledge:
   * a normal to a vector is always directed to the right
   *
   * idea:
   * add the normal vector of the line to any point of the line  -> the resulting point is to the right of the line
   * add the negative normal vector to the same point            -> the resulting point is to the left  of the line
   * now get the distances of these points to the point of interest
   * if that point is closer to the right point than to the left ->   that point itself is to the right of the line
   */

  double  vline[2] = {
    lineEnd[0] - lineStart[0],
    lineEnd[1] - lineStart[1]
  };
  double  normal[2] = {
    vline[1], -vline[0]
  };

  double  right[2] = {
    lineEnd[0] + normal[0],
    lineEnd[1] + normal[1]
  };
  double  left[2] = {
    lineEnd[0] - normal[0],
    lineEnd[1] - normal[1]
  };

  double  vright[2] = {
    point[0] - right[0], point[1] - right[1]
  };
  double  vleft[2] = {
    point[0] - left[0], point[1] - left[1]
  };

  /*
   * for comparing lengths of vectors there is no need to actually compute their length
   * comparing their squares of their respective lengths grant the same results
   * while saving some computation time (spare us the sqrt operation)
   */
  double  squaredDistanceRight  = scalarProduct2D(vright, vright);
  double  squaredDistanceLeft   = scalarProduct2D(vleft, vleft);
  short   ret                   = squaredDistanceRight < squaredDistanceLeft;

  return ret;
}


PRIVATE short
isToTheRightPointVector(const double  *lineStart,
                        const double  *lineVector,
                        const double  *point)
{
  double lineEnd[2] = {
    lineStart[0] + lineVector[0], lineStart[1] + lineVector[1]
  };

  return isToTheRightPointPoint(lineStart, lineEnd, point);
}


PRIVATE double
angleBetweenVectors2D(const double  vector1[2],
                      const double  vector2[2])
{
  double  vectorNormalized1[2] = {
    vector1[0], vector1[1]
  };
  double  vectorNormalized2[2] = {
    vector2[0], vector2[1]
  };

  normalize(vectorNormalized1);
  normalize(vectorNormalized2);
  double  cosAngle  = scalarProduct2D(vectorNormalized1, vectorNormalized2);
  double  angle     = 0.0;
  if (fabs(cosAngle - -1.00) < EPSILON_7) {
    /* cosAngle == -1 -> rad: PI deg: 180° */

    angle = MATH_PI;
  } else if (fabs(cosAngle - 1.00) < EPSILON_7) {
    /* cosAngle == +1 -> rad: 0 deg: 0° */

    angle = 0;
  } else {
    angle = acos(cosAngle);
  }

  return angle;
}


PRIVATE double
anglePtPtPt2D(const double  *p1,
              const double  *p2,
              const double  *p3)
{
  double  v1[2] = {
    p1[0] - p2[0], p1[1] - p2[1]
  };
  double  v2[2] = {
    p3[0] - p2[0], p3[1] - p2[1]
  };

  return angleBetweenVectors2D(v1, v2);
}


PRIVATE double
normalizeAngle(const double angle,
               short        useDegree)
{
  double  min   = 0.0;
  double  step  = useDegree ? 360.0 : MATH_TWO_PI;
  double  max   = min + step;

  int     viciousCircleCounter = 0;

  double  ret = angle;

  while (ret < min) {
    ret += step;
    viciousCircleCounter++;
    if (viciousCircleCounter > 1000000)

      break;
  }

  while (ret > max) {
    ret -= step;
    viciousCircleCounter++;
    if (viciousCircleCounter > 1000000)

      break;
  }

  return ret;
}


PRIVATE void
rotatePointAroundPoint(const double *point,
                       const double *rotationCenter,
                       const double angle,
                       double       *ret)
{
  double  x   = point[0];
  double  y   = point[1];
  double  x0  = rotationCenter[0];
  double  y0  = rotationCenter[1];
  double  phi = -angle;  /* negative value because we rotate clockwise for positive input */

  ret[0]  = x0 + (x - x0) * cos(phi) - (y - y0) * sin(phi);
  ret[1]  = y0 + (x - x0) * sin(phi) + (y - y0) * cos(phi);
}


PRIVATE void
rotateVectorByAngle(const double  *vector,
                    const double  angle,
                    double        *ret)
{
  double c[2] = {
    0, 0
  };

  rotatePointAroundPoint(vector, c, angle, ret);
}


PRIVATE void
translatePointByVector(const double *point,
                       const double *trans,
                       double       *ret)
{
  ret[0]  = point[0] + trans[0];
  ret[1]  = point[1] + trans[1];
}


PRIVATE short
solveSquareEquation(const double  a,
                    const double  b,
                    const double  c,
                    double        *sol1,
                    double        *sol2)
{
  short   ret   = 0;
  double  discr = b * b - 4 * a * c;

  if (discr < 0.0) {
    ret = 0;
    return ret;
  }

  if (discr == 0.0)
    ret = 1;
  else
    ret = 2;

  double  answer1 = (-b + sqrt(discr)) / (2 * a);
  double  answer2 = (-b - sqrt(discr)) / (2 * a);

  *sol1 = answer1;
  *sol2 = answer2;

  return ret;
}


PRIVATE short
getCutPointsOfCircles(const double  *c1,
                      const double  r1,
                      const double  *c2,
                      const double  r2,
                      double        *ret1,
                      double        *ret2)
{
  short   answers = -2;
  double  c1x     = c1[0];
  double  c1y     = c1[1];
  double  c2x     = c2[0];
  double  c2y     = c2[1];

  double  dx = c1x - c2x;

  dx = dx < 0 ? -dx : dx;
  double  dy = c1y - c2y;
  dy = dy < 0 ? -dy : dy;
  double  dr = r1 - r2;
  dr = dr < 0 ? -dr : dr;

  double  eps = 1.0;
  /*
   * if any delta is smaller than this epsilon this delta will be considered zero
   * i.e. the values that are being compared are treated as equal
   */

  short smallDX = dx < eps;
  short smallDY = dy < eps;
  short smallDR = dr < eps;

  if (smallDX && smallDY) {
    if (smallDR)
      /* circles coincide */
      answers = -1;
    else
      /* circles are concentric but with different radius */
      answers = 0;
  } else

  if (!smallDY) {
    /*
     * (smallDX || !smallDX) && !smallDY
     * EQ1: circle1: (r1)² = (c1x - x)² + (c1y - y)²
     * EQ2: circle2: (r2)² = (c2x - x)² + (c2y - y)²
     */

    /*
     * EQ3: EQ1 - EQ2, get y
     * EQ3: y = (x * k + l) / m = (k / m) * x + (l / m)
     */
    double  k = -2 * c1x + 2 * c2x;
    double  l = c1x * c1x - c2x * c2x + c1y * c1y - c2y * c2y - r1 * r1 + r2 * r2;
    double  m = (-1) * (-2 * c1y + 2 * c2y);

    /*
     * EQ4: replace y in EQ1 with EQ3
     * transform equation into ax²+bx+c=0 shape
     */
    double  p = c1y - (l / m);
    double  q = k / m;
    double  a = q * q + 1;
    double  b = -2 * c1x - 2 * p * q;
    double  c = c1x * c1x + p * p - r1 * r1;

    double  sol1;
    double  sol2;
    answers = solveSquareEquation(a, b, c, &sol1, &sol2);

    if (answers > 0) {
      ret1[0] = sol1;
      ret1[1] = (sol1 * k + l) / m;
    }

    if (answers > 1) {
      ret2[0] = sol2;
      ret2[1] = (sol2 * k + l) / m;
    }
  } else {
    /* smallDY && !smallDX */

    double  k = -2 * c1y + 2 * c2y;
    double  l = (c1x * c1x - c2x * c2x) + (c1y * c1y - c2y * c2y) + (r2 * r2 - r1 * r1);
    double  m = (-1) * (-2 * c1x + 2 * c2x);

    double  p = c1x - (l / m);
    double  q = k / m;
    double  a = q * q + 1;
    double  b = -2 * c1y - 2 * p * q;
    double  c = c1y * c1y + p * p - r1 * r1;

    double  sol1;
    double  sol2;

    answers = solveSquareEquation(a, b, c, &sol1, &sol2);

    if (answers == 0)
      printf("no solution 2: %3.2lf %3.2lf %3.2lf\n", a, b, c);

    if (answers > 0) {
      ret1[1] = sol1;
      ret1[0] = (sol1 * k + l) / m;
    }

    if (answers > 1) {
      ret2[1] = sol2;
      ret2[0] = (sol2 * k + l) / m;
    }
  }

  return answers;
}


PRIVATE short
getCutPointsOfCircleAndLine(const double  *center,
                            const double  radius,
                            const double  *anchor,
                            const double  *direction,
                            double        *ret1,
                            double        *ret2)
{
  double  a = direction[0] * direction[0] + direction[1] * direction[1];
  double  b = 2 * direction[0] * (anchor[0] - center[0]) + 2 * direction[1] *
              (anchor[1] - center[1]);
  double  c = (anchor[0] - center[0]) * (anchor[0] - center[0]) + (anchor[1] - center[1]) *
              (anchor[1] - center[1]) - radius * radius;

  double  solution1, solution2;
  short   answers = solveSquareEquation(a, b, c, &solution1, &solution2);

  if (answers > 0) {
    ret1[0] = anchor[0] + solution1 * direction[0];
    ret1[1] = anchor[1] + solution1 * direction[1];
  }

  if (answers > 1) {
    ret2[0] = anchor[0] + solution2 * direction[0];
    ret2[1] = anchor[1] + solution2 * direction[1];
  }

  return answers;
}


PRIVATE void
vector(const double pStart[2],
       const double pEnd[2],
       double       v[2])
{
  v[0]  = pEnd[0] - pStart[0];
  v[1]  = pEnd[1] - pStart[1];
}


PRIVATE void
normal(const double v[2],
       double       n[2])
{
  double  vNormal[2];

  vNormal[0]  = v[1];
  vNormal[1]  = -v[0];

  double  vNormalUnit[2];
  unit(vNormal, vNormalUnit);

  n[0]  = vNormalUnit[0];
  n[1]  = vNormalUnit[1];
}


PRIVATE void
unit(const double v[2],
     double       u[2])
{
  double length = vectorLength2D(v);

  u[0]  = v[0] / length;
  u[1]  = v[1] / length;
}


PRIVATE void
circle(const double P1[2],
       const double P2[2],
       const double P3[2],
       double       center[2],
       double       *radius)
{
  /* initialize coefficients of the system of equations */
  double alpha[3];

  alpha[0]  = 1;
  alpha[1]  = 1;
  alpha[2]  = 1;

  double  beta[3];
  beta[0] = -P1[0];
  beta[1] = -P2[0];
  beta[2] = -P3[0];

  double  gamma[3];
  gamma[0]  = -P1[1];
  gamma[1]  = -P2[1];
  gamma[2]  = -P3[1];

  double  r[3];
  r[0]  = -(P1[0] * P1[0] + P1[1] * P1[1]);
  r[1]  = -(P2[0] * P2[0] + P2[1] * P2[1]);
  r[2]  = -(P3[0] * P3[0] + P3[1] * P3[1]);

  /* eliminate alpha 1 und 2 */
  beta[1]   -= beta[0];
  gamma[1]  -= gamma[0];
  beta[2]   -= beta[0];
  gamma[2]  -= gamma[0];
  r[1]      -= r[0];
  r[2]      -= r[0];

  double  A;
  double  B;
  double  C;

  /* solve linear equations */
  if ((fabs(beta[1]) < EPSILON_7)
      && (fabs(gamma[1]) > EPSILON_7)) {
    C = r[1] / gamma[1];
    B = (r[2] - gamma[2] * C) / beta[2];
    /* printf("beta[1] == %15.12lf\n",beta [1]); */
  } else if ((fabs(beta[2]) < EPSILON_7)
             && (fabs(gamma[2]) > EPSILON_7)) {
    C = r[2] / gamma[2];
    B = (r[1] - gamma[1] * C) / beta[1];
    /* printf("beta[2] == %15.12lf\n", beta[2]); */
  } else {
    if (fabs(gamma[1]) < EPSILON_7) {
      B = r[1] / beta[1];
      C = (r[2] - beta[2] * B) / gamma[2];
      /* printf("gamma[1] == %15.12lf\n", gamma[1]); */
    } else if (fabs(gamma[2]) < EPSILON_7) {
      B = r[2] / beta[2];
      C = (r[1] - beta[1] * B) / gamma[1];
      /* printf("gamma[2] == %15.12lf\n", gamma[2]); */
    } else {
      /*
       * printf("  GL1: %15.12lf * B + %15.12lf * C = %15.12lf\n", beta[1],
       * gamma[1], r[1]);
       * printf("  GL2: %15.12lf * B + %15.12lf * C = %15.12lf\n", beta[2],
       * gamma[2], r[2]);
       */

      gamma[2]  = gamma[2] * beta[1] - gamma[1] * beta[2];
      r[2]      = r[2] * beta[1] - r[1] * beta[2];

      /* printf("  GL C: %15.12lf * C = %15.12lf\n", gamma[2], r[2]); */
      C = r[2] / gamma[2];
      B = (r[1] - gamma[1] * C) / beta[1];
    }
  }

  A = r[0] - beta[0] * B - gamma[0] * C;

  /* calculate center and radius */
  center[0] = B / 2;
  center[1] = C / 2;
  *radius   = sqrt(center[0] * center[0] + center[1] * center[1] - A);
}


#endif
