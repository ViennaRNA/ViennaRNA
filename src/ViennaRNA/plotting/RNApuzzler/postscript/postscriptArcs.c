#include "ViennaRNA/RNApuzzler/postscript/postscriptArcs.h"
#include "ViennaRNA/RNApuzzler/vector_math.h"
#include "ViennaRNA/utils.h"

#include <stdlib.h>

void calcArc(
        const double center[2],
        const double radius,
        const short goClockwise,
        const int i,
        double * const x,
        double * const y,
        double *arcCoords
) {

    int from = i-1;
    int to   = i;
    double pFrom[2]   = { x[from],  y[from]  };
    double pTo[2]     = { x[to],    y[to]    };
    double pCenter[2] = { center[0], center[1] };

    double vCenterFrom[2], vCenterTo[2];
    vector(pCenter, pFrom, vCenterFrom);
    vector(pCenter, pTo  , vCenterTo  );
    double v_1_0[2] = { 1.0, 0.0 }; // zero degree axis

    double angleFrom = toDegree(angleBetweenVectors2D(v_1_0, vCenterFrom));
    double angleTo   = toDegree(angleBetweenVectors2D(v_1_0, vCenterTo));

    if (pFrom[1] < pCenter[1]) {
        angleFrom = 360.0 - angleFrom;
    }
    if (pTo[1] < pCenter[1]) {
        angleTo = 360.0 - angleTo;
    }

    arcCoords[6*i + 0] = pCenter[0];
    arcCoords[6*i + 1] = pCenter[1];
    arcCoords[6*i + 2] = radius;
    arcCoords[6*i + 3] = angleFrom;
    arcCoords[6*i + 4] = angleTo;
    arcCoords[6*i + 5] = goClockwise;
}

void calcArcsHandleStem(
        int start,
        short const * const pair_table,
        double * const x,
        double * const y,
        const tBaseInformation* baseInformation,
        double *arcCoords
);

void calcArcsHandleLoop(
    int start,
    short const * const pair_table,
    double * const x,
    double * const y,
    const tBaseInformation* baseInformation,
    double *arcCoords
) {
    int end = pair_table[start];

    /// ---------------------------------------------------
    /// count the number of points / bases on the loop
    /// ---------------------------------------------------
    int numPoints = 1;
    int i = start + 1;
    while (i < end) {
        if (pair_table[i] == 0) {
            i++;
        } else
        if (pair_table[i] > i) {
            i = pair_table[i];
        } else {
            i++;
        }
        numPoints++;
    }

    /// ---------------------------------------------------
    /// collect the list of all points / bases on the loop
    /// ---------------------------------------------------
    double** points = (double**) vrna_alloc(numPoints * sizeof(double*));
    for (int k = 0; k < numPoints; k++) {
        double* point = (double*) vrna_alloc(2 * sizeof(double));
        points[k] = point;
    }

    int k = 0;
    i = start + 1;
    while (i < end) {
        double* point = points[k];
        point[0] = x[i - 1];
        point[1] = y[i - 1];
        k++;
        if (pair_table[i] == 0) {
            i++;
        } else if (pair_table[i] > i) {
            // ... and meanwhile handle all stems
            calcArcsHandleStem(i, pair_table, x, y, baseInformation, arcCoords);
            i = pair_table[i];
        } else {
            i++;
        }
    }
    double* point = points[k];
    point[0] = x[i - 1];
    point[1] = y[i - 1];

    /// take the line from the loop's end to its start
    /// and an arbitrary point other point at the loop
    /// (we take the point which is the center value in the points list)
    /// if the chosen point is to the right of the line
    /// we can safely state the circle is traversed clockwise
    /// (counter clockwise only happens if the drawing is
    ///  allowed to flip neighboring exterior children)
    short goClockwise = isToTheRightPointPoint(points[numPoints-1],
                                               points[0],
                                               points[numPoints / 2]);

    double center[2];
    double rad;
    circle(points[0 * numPoints / 3], points[1 * numPoints / 3], points[2 * numPoints / 3], center, &rad);

    // free points array as it is no longer needed
    for (int k = 0; k < numPoints; k++) {
        free(points[k]);
    }
    free(points);

    /// ----------------------------------
    /// finally calculate arcs
    /// ----------------------------------
    i = start + 1;
    while (i < end) {
        if (pair_table[i] == 0) {
            calcArc(center, rad, goClockwise, i-1, x, y, arcCoords);
            i++;
        } else if (pair_table[i] > i) {
            calcArc(center, rad, goClockwise, i-1, x, y, arcCoords);
            i = pair_table[i];
        } else {
            i++;
        }
    }
    calcArc(center, rad, goClockwise, end-1, x, y, arcCoords);
}

void calcArcsHandleStem(
    int start,
    short const * const pair_table,
    double * const x,
    double * const y,
    const tBaseInformation* baseInformation,
    double *arcCoords
) {
    int i = start;
    config *cfg = baseInformation[i].config;
    while (cfg == NULL) {
        i++;
        cfg = baseInformation[i].config;
    }
    calcArcsHandleLoop(i, pair_table, x, y, baseInformation, arcCoords);
}

//------------------------------------------------------------------------------

void computeAnglesAndCentersForPS(
    short const * const pair_table,
    double * const x,
    double * const y,
    const tBaseInformation* baseInformation,
    double *arcCoords
) {
    int end = pair_table[0];

    // initialize (again??)
    for (int j = 0; j < end; j++) {
        arcCoords[6*j+0] = -1.;
        arcCoords[6*j+1] = -1.;
        arcCoords[6*j+2] = -1.;
        arcCoords[6*j+3] = -1.;
        arcCoords[6*j+4] = -1.;
        arcCoords[6*j+5] = -1.;
    }

    int i = 1;
    while (i < end) {
        if (pair_table[i] == 0) {
            i++;
        } else if (pair_table[i] > i) {
            calcArcsHandleStem(i, pair_table, x, y, baseInformation, arcCoords);
            i = pair_table[i];
        } else {
            i++;
        }
    }
}

