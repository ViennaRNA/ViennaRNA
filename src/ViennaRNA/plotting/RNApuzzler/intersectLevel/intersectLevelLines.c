#include "ViennaRNA/RNApuzzler/intersectLevel/intersectLevelLines.h"
#include "ViennaRNA/RNApuzzler/definitions.h"
#include "ViennaRNA/RNApuzzler/vector_math.h"

#include <math.h>
#include <stdlib.h>
#include <stdio.h>

/**
 * @brief matchPointArc
 *      checks if the given point is situated on the given wedge
 * @param point double[2] array with point coordinates
 * @param arc double[6] array with arc coordinates [circleX, circleY, circleR, arcFrom, arcTo, clockwise]
 * @return 1 if point is situated on the wedge, 0 otherwise
 */
short matchPointArc(
    const double point[2],
    const double arc[6]
) {
    /// anpassen !!!
    double center[2] = {arc[0], arc[1]};
    double from = toRad(arc[3]);
    double to = toRad(arc[4]);
    short clockwise = (arc[5] > 0.5);

    double v_center_point[2];
    vector(center, point, v_center_point);

    // 0° is vector (1,0)
    // angle_cut <= 180° <-> cut[1] >= arc_y
    // berechne winkel zwischen (1,0) und v_center_cut
    // setze winkel korrekt (ggf. angle = 360° - angle)
    // prüfe, ob zwischen from und to
    double zero_degree[2] = { 1, 0 };
    double angle = angleBetweenVectors2D(v_center_point, zero_degree);
    if (point[1] < center[1]) {
        angle = MATH_TWO_PI - angle;
    }

    short ret = 0;
    if (clockwise) {
        if (from > to) {
            // normal case
            // check: from < angle < to
            ret = (from >= angle && angle >= to);
        } else {
            // interval surpasses 360° border
            // check: from < angle < 360 or 0 < angle < to
            ret = (from >= angle && angle >= 0) || (MATH_TWO_PI >= angle && angle >= to);
        }
    } else {
        if (from < to) {
            // normal case
            // check: to < angle < from
            ret = (from <= angle && angle <= to);
        } else {
            // interval surpasses 360° border
            // check: from > angle > 0 || MATH_TWO_PI > angle > to
            ret = (from <= angle && angle <= MATH_TWO_PI) || (0 <= angle && angle <= to);
        }
    }

    return ret;
}

short matchLinePoint(const double pLine[2], const double dirLine[2], const double p[2]) {
    /// line...
    ///     pX = pLineX + t * dirLineX
    ///     pY = pLineY + t * dirLineY
    ///
    ///     t = (pX - pLineX) / dirLineX    if dirLineX != 0
    ///     t = (pY - pLineY) / dirLineY    if dirLineY != 0

    double t = -1.0; // init as not matching
    if (fabs(dirLine[0]) > 0.0001) {    // != 0.0
        t = (p[0] - pLine[0]) / dirLine[0];
    }
    else if (fabs(dirLine[1]) > 0.0001) {    // != 0.0
        t = (p[1] - pLine[1]) / dirLine[1];
    }
    else {
        return 0; // this is not even a line since dirLine == (0.0, 0.0)
    }

    return (0.0 <= t && t <= 1.0);
//    /// check if the reference point matches the input point
//    double refP[2] = { pLine[0] + t * dirLine[0]
//                     , pLine[1] + t * dirLine[1] };

//    double vRef[2];
//    vector(p, refP, vRef);
//    double delta = vectorLength2D(vRef);
//    short ret = (0.0 <= t && t <= 1.0) && (delta < 0.1);
//    return ret;
}

short intersectCircleCircle(
    const double c1[2], const double c1r,
    const double c2[2], const double c2r
) {
    double v_c1_c2[2];
    vector(c1, c2, v_c1_c2);
    double distance = vectorLength2D( v_c1_c2 );

    short intersect = (distance < (c1r + c2r));
//    printf("C1=(%3.2f, %3.2f) r1=%3.2f ... C2=(%3.2f, %3.2f) r2=%3.2f\n", c1x, c1y, c1r, c2x, c2y, c2r);

    return intersect;
}

short intersectLineSegments(
        const double A[2], const double B[2],
        const double X[2], const double Y[2],
        double P[2]
) {
    if ((   X[0] < A[0] - epsilon7 && X[0] < B[0] - epsilon7
         && Y[0] < A[0] - epsilon7 && Y[0] < B[0] - epsilon7)
        ||
        (   X[0] > A[0] + epsilon7 && X[0] > B[0] + epsilon7
         && Y[0] > A[0] + epsilon7 && Y[0] > B[0] + epsilon7)
        ) {
        /// Check if the x-coordinates of X and Y are smaller than
        /// the x-coordinates of A and B -> lines can not intersect
        return 0;
    }

    if ((   X[1] < A[1] - epsilon7 && X[1] < B[1] - epsilon7
         && Y[1] < A[1] - epsilon7 && Y[1] < B[1] - epsilon7)
        ||
        (   X[1] > A[1] + epsilon7 && X[1] > B[1] + epsilon7
         && Y[1] > A[1] + epsilon7 && Y[1] > B[1] + epsilon7)
        ) {
        /// Check if the y-coordinates of X and Y are smaller than
        /// the y-coordinates of A and B -> lines can not intersect
        return 0;
    }

    double denominator = (B[0] - A[0]) * (X[1] - Y[1]) - (B[1] - A[1]) * (X[0] - Y[0]);
    if (fabs(denominator) < epsilon7) {
        /// lines are parallel

        /// check if X is situated on AB line
        double sX, sY;

        double dx = B[0] - A[0];
        double dy = B[1] - A[1];

        if (fabs(dx) > epsilon7) {
            sX = (X[0] - A[0]) / (dx);
            double refXy = A[1] + sX * (dy);
            if (fabs(refXy - X[1]) > epsilon7) {
                /// AB and XY are not part of the same line
                return 0;
            }
            sY = (Y[0] - A[0]) / (dx);
        } else {
            sX = (X[1] - A[1]) / (dy);
            double refXx = A[0] + sX * (dx);
            if (fabs(refXx - X[0]) > epsilon7) {
                /// AB and XY are not part of the same line
                return 0;
            }
            sY = (Y[1] - A[1]) / (dy);
        }

        /// check if X or Y are situated directly on AB
        if (   (0.0 <= sX && sX <= 1.0)
            || (0.0 <= sY && sY <= 1.0)) {
            return 1;
        }

        /// check if XY encloses AB
        if (   (sX < 0.0 && 1.0 < sY)
            || (sY < 0.0 && 1.0 < sX)) {
            return 1;
        }

    } else {
        /// lines are not parallel and might intersect
        /// (default case)

        double nominatorS = (X[0] - Y[0]) * (A[1] - X[1]) - (X[1] - Y[1]) * (A[0] - X[0]);
        double nominatorT = (A[0] - X[0]) * (B[1] - A[1]) - (A[1] - X[1]) * (B[0] - A[0]);
        double s = nominatorS / denominator;
        double t = nominatorT / denominator;

        if (0.0 <= s && s <= 1.0 && 0.0 <= t && t <= 1.0) {

            double Ps[2];
            Ps[0] = A[0] + s * (B[0] - A[0]);
            Ps[1] = A[1] + s * (B[1] - A[1]);

            double Pt[2];
            Pt[0] = X[0] + t * (Y[0] - X[0]);
            Pt[1] = X[1] + t * (Y[1] - X[1]);

            if (fabs(Ps[0] - Pt[0]) < epsilon7 && fabs(Ps[1] - Pt[1]) < epsilon7) {
                if (P != NULL) {
                  P[0] = Ps[0];
                  P[1] = Ps[1];
                }
                return 1;
            } else {
                // real difference
                printf("[DZ] intersectLineSegments: real difference = %15lf\n", fabs(Ps[0] - Pt[0]));
            }
        }
    }

    return 0;
}

short intersectLineArc(
    const double point_1[2],
    const double point_2[2],
    const double arc[6]
) {
    char *fnName = "intersectLineArc";

    //printf("L_A\n");
    //printf("from: %3.2f° to: %3.2f° arc[5]:%d angle:%3.2f°");

    double cut[2][2];
    //printf("Circle=Circle(%3.2f, %3.2f, %3.2f)\n", arc[0], arc[1], arc[2]);

    double center[2] = { arc[0], arc[1] };
    double radius = arc[2];
    double anchor[2] = { point_1[0], point_1[1] };
    double direction[2];
    vector(point_1, point_2, direction);
    short num_points = getCutPointsOfCircleAndLine(center, radius, anchor, direction, cut[0], cut[1]);

    short ret = 0;
    for (int i = 0; i < num_points; i++) {
        /// check if the computed intersection point is situated on the line segment
        double A[2] = { point_1[0], point_1[1] };
        double B[2] = { point_2[0], point_2[1] };
        double line[2];
        vector(A, B, line);
        double length = vectorLength2D(line);
        double v_A_cut[2];
        double v_B_cut[2];
        vector(A, cut[i], v_A_cut);
        vector(B, cut[i], v_B_cut);
//        printf("Cut%d=(%3.2f, %3.2f) ", i, cut[i][0], cut[i][1]);
        if (fabs(length - vectorLength2D(v_A_cut) - vectorLength2D(v_B_cut)) > 0.01) {
            /// consider double errors... should test for !=0
//            printf("skip\n");
            continue;
        }
//        printf("test\n");

        /// check if the computed intersection point is situated between angle_from and angle_to
        ret = ret || matchPointArc(cut[i], arc);

        if (ret) { break; }
    }

    return ret;
}

short intersectArcArc(
    const double arc1[6],
    const double arc2[6]
) {
    char *fnName = "intersectArcArc";

    double c1[2] = { arc1[0], arc1[1] };
    double r1 = arc1[2];
    double c2[2] = { arc2[0], arc2[1] };
    double r2 = arc2[2];
    if (!intersectCircleCircle(c1, r1, c2, r2)) {
        return 0;
    }

//    printf("Arc v Arc\n");

    /// get intersection points of corresponding circles
    /// then check if they match the interval of the other arc respectively
    //printf("[IMPLEMENT ME] intersectArcArc\n");

    double cut[2][2];
    short num_points = getCutPointsOfCircles(c1, r1, c2, r2, cut[0], cut[1]);
//    if (num_points <= 0) {
//            printf("[ERROR] [%s] (%12.8lf %12.8lf %12.8lf) -- (%12.8lf %12.8lf %12.8lf)\n",
//                   fnName,
//                   c1[0], c1[1], r1,
//                   c2[0], c2[1], r2
//            );
//    }
    //    printf("cut points: %d\n", num_points);
    // 0 - no common points; -1 - circles are equal
    // only test for 1 and 2 cut points
    // in other cases return false

    short ret = 0;
    for (int i = 0; i < num_points; i++) {
        short hit1 = matchPointArc(cut[i], arc1);
        short hit2 = matchPointArc(cut[i], arc2);

        ret = ret || (hit1 && hit2);

        if (ret) {
            printf("[INFO] [%s] Cut[%d]=(%3.2f, %3.2f)\n", fnName, i, cut[i][0], cut[i][1]);
            printf("[INFO] [%s] (%12.8lf %12.8lf %12.8lf) -- (%12.8lf %12.8lf %12.8lf)\n",
                   fnName,
                   c1[0], c1[1], r1,
                   c2[0], c2[1], r2
            );
            //break;
        }
    }

    return ret;
}
