#ifndef INTERSECTLEVELLINES_H
#define INTERSECTLEVELLINES_H
/**
 * @brief intersectLineArc
 *      checks for intersection of a circle arc and straight line
 * @param point_1 first point of line
 * @param point_2 second point of line
 * @param arc arc
 * @return 1 if intersecting, 0 otherwise
 */
short intersectLineArc(
    const double point_1[2],
    const double point_2[2],
    const double arc[6]
);

short intersectArcArc(
    const double arc1[6],
    const double arc2[6]
);


short intersectCircleCircle(
    const double c1[2], const double c1r,
    const double c2[2], const double c2r
);

short matchLinePoint(
    const double pLine[2],
    const double dirLine[2],
    const double p[2]
);

/**
 * @brief intersectLineSegments
 *      Determines if two line segments defined by the points A,B and X,Y each intersect.
 *      If they do while both lines are not parallel the coordinates of their cut point are returned in P.
 * @param A start point of first line segment
 * @param B end point of first line segment
 * @param X start point of second line segment
 * @param Y end point of second line segment
 * @param P return value for cut point of both line segments
 * @return 1 if intersecting , 0 otherwise
 */
short intersectLineSegments(
        const double A[2], const double B[2],
        const double X[2], const double Y[2],
        double P[2]
);

#endif
