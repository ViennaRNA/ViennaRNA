#include "ViennaRNA/RNApuzzler/data/boundingBoxes.h"
#include "ViennaRNA/RNApuzzler/dataTypes/config_struct.h"
#include "ViennaRNA/RNApuzzler/data/cfg_reader.h"
#include "ViennaRNA/RNApuzzler/vector_math.h"

#include "ViennaRNA/utils.h"

#include <math.h>
#include <stdio.h>

void getBulgeXY(
    const stemBox* stem,
    const int index,
    double *x,
    double *y
) {
    double* bulge = stem->bulges[index];
    *x = stem->c[0] + bulge[2] * stem->a[0] + bulge[0] * stem->b[0] * (stem->e[1] + stem->bulgeDist);
    *y = stem->c[1] + bulge[2] * stem->a[1] + bulge[0] * stem->b[1] * (stem->e[1] + stem->bulgeDist);
}

void getBulgeCoordinatesExtraDistance(
    const stemBox* stem,
    const int index,
    const double extraDistance,
    double pPrev[2],
    double pThis[2],
    double pNext[2]
) {
    double* bulge = stem->bulges[index];
    pPrev[0] = stem->c[0] + bulge[1] * stem->a[0] + bulge[0] * stem->b[0] *  stem->e[1];
    pPrev[1] = stem->c[1] + bulge[1] * stem->a[1] + bulge[0] * stem->b[1] *  stem->e[1];
    pThis[0] = stem->c[0] + bulge[2] * stem->a[0] + bulge[0] * stem->b[0] * (stem->e[1] + extraDistance + stem->bulgeDist);
    pThis[1] = stem->c[1] + bulge[2] * stem->a[1] + bulge[0] * stem->b[1] * (stem->e[1] + extraDistance + stem->bulgeDist);
    pNext[0] = stem->c[0] + bulge[3] * stem->a[0] + bulge[0] * stem->b[0] *  stem->e[1];
    pNext[1] = stem->c[1] + bulge[3] * stem->a[1] + bulge[0] * stem->b[1] *  stem->e[1];
}

void getBulgeCoordinates(
    const stemBox* stem,
    const int index,
    double pPrev[2],
    double pThis[2],
    double pNext[2]
) {
    getBulgeCoordinatesExtraDistance(stem, index, 0, pPrev, pThis, pNext);
}

void translateLoopBox(loopBox* box, const double* vector) {
    double center[2] = { box->c[0], box->c[1] };
    double newCenter[2];
    translatePointByVector(center, vector, newCenter);
    box->c[0] = newCenter[0];
    box->c[1] = newCenter[1];
}

void rotateLoopBox(loopBox* box, const double* point, const double angle) {
    double loopCenter[2] = { box->c[0], box->c[1] };
    double newLoopCenter[2];

    rotatePointAroundPoint(loopCenter, point, angle, newLoopCenter);

    box->c[0] = newLoopCenter[0];
    box->c[1] = newLoopCenter[1];
}

void translateStemBox(stemBox* box, const double* vector) {
    double center[2] = { box->c[0], box->c[1] };
    double newCenter[2];
    translatePointByVector(center, vector, newCenter);
    box->c[0] = newCenter[0];
    box->c[1] = newCenter[1];
}

void rotateStemBox(stemBox* box, const double* point, const double angle) {
    double stemCenter[2] = { box->c[0], box->c[1] };
    double stemDirA[2] =   { box->a[0], box->a[1] };
    double stemDirB[2] =   { box->b[0], box->b[1] };
    double newStemCenter[2];
    double newStemDirA[2];
    double newStemDirB[2];

    rotatePointAroundPoint(stemCenter, point, angle, newStemCenter);
    rotateVectorByAngle(stemDirA, angle, newStemDirA);
    rotateVectorByAngle(stemDirB, angle, newStemDirB);

    box->c[0] = newStemCenter[0];
    box->c[1] = newStemCenter[1];
    box->a[0] = newStemDirA[0];
    box->a[1] = newStemDirA[1];
    box->b[0] = newStemDirB[0];
    box->b[1] = newStemDirB[1];
}

void getLoopData(
    double* center,
    double* radius,
    const int start,
    const short* const pair_table,
    const tBaseInformation* baseInformation,
    const double* x,
    const double* y
) {
    // copy/pasted from plot_layouts.c ... consider moving this one to utils or such a thing.

    int i = start;
    int end = pair_table[start];
    config* cfg = baseInformation[i].config;
    double r = cfg->radius;

    // Reminder to offset -1
    // x, y and arc_coords as well have their entry for the first base (i=1) at index 0
    // this results in a offset -1 at every reading access to x, y
    // consider the stem just before the loop
    // if i's partner is to the right of the loop's first line
    // then this loop is directed clockwise, counter clockwise otherwise
    double current[2] = {x[start-1], y[start-1]};
    double next[2] = {x[(start+1)-1], y[(start+1)-1]};
    double last[2] = {x[end-1], y[end-1]};
    short go_clockwise = isToTheRightPointPoint(current, next, last);

    double v_pair[2];
    vector(last, current, v_pair);

    double v_normal[2];
    normal(v_pair, v_normal);

    double pair_length = vectorLength2D(v_pair); // = paired

    double center_dist = sqrt(r*r - 0.25*pair_length*pair_length);

    // for clockwise go to the right... and to the left for counter clockwise loop
    short dir = go_clockwise ? 1 : -1;
    center[0] = last[0] + 0.5 * v_pair[0] + dir * center_dist * v_normal[0];
    center[1] = last[1] + 0.5 * v_pair[1] + dir * center_dist * v_normal[1];
    *radius = r;
}

loopBox* createLoopBox(
    const double center[2],
    const double radius
) {
    loopBox* box = (loopBox*) vrna_alloc(1*sizeof(loopBox));

    box->c[0] = center[0];
    box->c[1] = center[1];
    box->r    = radius;

    return box;
}

loopBox* buildLoopBox(
    const int start,
    const short* const pair_table,
    const tBaseInformation* baseInformation,
    const double* x,
    const double* y
) {
    double center[2];
    double radius;

    // calculate center coords and radius for the loop
    getLoopData(center, &radius, start, pair_table, baseInformation, x, y);

    //printf("lBox [%3.2f %3.2f] r:%3.2f\n", center_x, center_y, radius);

    loopBox* box = createLoopBox(center, radius);
    return box;
}

stemBox* createStemBox(
    const double s[2],
    const double e[2],
    const double sp[2]
) {
    stemBox* box = (stemBox*) vrna_alloc(1*sizeof(stemBox));

    double a[2] = { 0.5 * (e[0] -  s[0]), 0.5 * (e[1] -  s[1]) };
    double b[2] = { 0.5 * (s[0] - sp[0]), 0.5 * (s[1] - sp[1]) };

    double length_a = vectorLength2D( a );
    double length_b = vectorLength2D( b );

    if (length_a == 0) {
        // solve this using b's normal vector
        normal(b, a);
        // make a have length 0.1 to create a proper bounding box
        length_a = 0.1;
        a[0] = a[0] * length_a;
        a[1] = a[1] * length_a;
    }

    box->a[0] = a[0] / length_a;
    box->a[1] = a[1] / length_a;
    box->b[0] = b[0] / length_b;
    box->b[1] = b[1] / length_b;
    box->c[0] = s[0] + a[0] - b[0];
    box->c[1] = s[1] + a[1] - b[1];
    box->e[0] = length_a;
    box->e[1] = length_b;

    return box;
}

int countBulges(
    const short* const pair_table,
    const int start,
    const int end
) {
    int bulgeCount = 0;

    for (int i = start; i < end; i++) {
        if (pair_table[i] == 0) {
            bulgeCount++;
        }
    }

    for (int i = pair_table[end]; i < pair_table[start]; i++) {
        if (pair_table[i] == 0) {
            bulgeCount++;
        }
    }

    return bulgeCount;
}

double getA(
    const stemBox* box,
    const double x,
    const double y
) {
    double a[2] = { box->a[0], box->a[1] };
    double b[2] = { box->b[0], box->b[1] };
    double c[2] = { box->c[0], box->c[1] };
    double p[2] = { x - c[0], y - c[1] };

    double ret = 0.0;
    if (b[0] == 0.0) {
        ret = p[0] / a[0];
    }
    else if (b[1] == 0.0) {
        ret = p[1] / a[1];
    }
    else {
        ret = ( p[0] * b[1] - p[1] * b[0] ) / ( a[0] * b[1] - a[1] * b[0] );
    }
    return ret;
}

double * createBulge(
    const stemBox* box,
    const double* x,
    const double* y,
    const int i,
    double bSign
) {
    double* bulge = (double*) vrna_alloc(4 * sizeof(double));

    // remember -1 offset between
    double aPrev = getA(box, x[(i-1) -1], y[(i-1) -1]);
    double aThis = getA(box, x[(i-1) +0], y[(i-1) +0]);
    double aNext = getA(box, x[(i-1) +1], y[(i-1) +1]);

    bulge[0] = bSign;
    bulge[1] = aPrev;
    bulge[2] = aThis;
    bulge[3] = aNext;

    return bulge;
}

void setBulges(
    stemBox* box,
    const short* const pair_table,
    const int start,
    const int end,
    const double* x,
    const double* y,
    const int bulgeCount,
    const double bulgeDist
) {
    if (bulgeCount <= 0) {
        box->bulges = NULL;
        box->bulgeCount = 0;
        box->bulgeDist = bulgeDist;
        return;
    }

    double** bulges = (double**) vrna_alloc(bulgeCount * sizeof(double*));
    int currentBulge = 0;
    for (int i = start; i < end; i++) {
        if (pair_table[i] == 0) {
            double bSign = 1.0;
            double* bulge = createBulge(box, x, y, i, bSign);
            bulges[currentBulge] = bulge;
            currentBulge++;
        }
    }

    for (int i = pair_table[end]; i < pair_table[start]; i++) {
        if (pair_table[i] == 0) {
            double bSign = -1.0;
            double* bulge = createBulge(box, x, y, i, bSign);
            bulges[currentBulge] = bulge;
            currentBulge++;
        }
    }

    box->bulgeCount = bulgeCount;
    box->bulgeDist = bulgeDist;
    box->bulges = bulges;
}

stemBox* buildStemBox(
    const int start,
    const int end,
    const short* const pair_table,
    const double* x,
    const double* y,
    const double bulgeDist
) {
    int i_s = start;
    int i_e = end;
    int i_sp = pair_table[start];

    /// get coordinates for rectangle corners
    // -1 for offset pair_table vs. x/y
    double s[2]  = {x[i_s -1], y[i_s -1]};
    double e[2]  = {x[i_e -1], y[i_e -1]};
    double sp[2] = {x[i_sp-1], y[i_sp-1]};

    /// finally create the box
    stemBox* box = createStemBox(s, e, sp);

    /// get coordinates for bulges
    int bulgeCount = countBulges(pair_table, i_s, i_e);
    setBulges(box, pair_table, i_s, i_e, x, y, bulgeCount, bulgeDist);

    return box;
}

void printLBox(const loopBox* loop) {
    printf("lBox [%3.2f, %3.2f] r:%3.2f\n",
            loop->c[0],
            loop->c[1],
            loop->r);
}

void printSBox(const stemBox* stem) {
    const short type = 1;
    if (type == 0) {
        double a[2] = { stem->a[0] * stem->e[0], stem->a[1] * stem->e[0] };
        double b[2] = { stem->b[0] * stem->e[1], stem->b[1] * stem->e[1] };

        // s  = c - a + b
        // e  = c + a + b
        // ep = c + a - b
        // sp = c - a - b
        double s[2]  = {  stem->c[0] - a[0] + b[0], stem->c[1] - a[1] + b[1] };
        double e[2]  = {  stem->c[0] + a[0] + b[0], stem->c[1] + a[1] + b[1] };
        double ep[2] = {  stem->c[0] + a[0] - b[0], stem->c[1] + a[1] - b[1] };
        double sp[2] = {  stem->c[0] - a[0] - b[0], stem->c[1] - a[1] - b[1] };

        printf("sBox [%3.2f, %3.2f] -> [%3.2f, %3.2f] -> [%3.2f, %3.2f] -> [%3.2f, %3.2f]\n",
                       s[0],  s[1],      e[0],  e[1],     ep[0], ep[1],     sp[0], sp[1]);
    }
    if (type == 1) {
        printf("sBox a=(%3.2f, %3.2f) b=(%3.2f, %3.2f) c=(%3.2f, %3.2f)\n"
                , stem->a[0]
                , stem->a[1]
                , stem->b[0]
                , stem->b[1]
                , stem->e[0]
                , stem->e[1]
                );
    }
}

void getLBoxCenter(const loopBox* box, double c[2]) {
    c[0] = box->c[0];
    c[1] = box->c[1];
}

void getSBoxCenter(const stemBox* box, double c[2]) {
    c[0] = box->c[0];
    c[1] = box->c[1];
}

