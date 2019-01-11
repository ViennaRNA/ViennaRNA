#ifndef BOUNDING_BOXES_H
#define BOUNDING_BOXES_H

#include <ViennaRNA/RNApuzzler/dataTypes/boundingBoxes_struct.h>
#include <ViennaRNA/RNApuzzler/dataTypes/tBaseInformation_struct.h>

void translateLoopBox(loopBox* box, const double* vector);
void translateStemBox(stemBox* box, const double* vector);
void rotateLoopBox(loopBox* box, const double* point, const double angle);
void rotateStemBox(stemBox* box, const double* point, const double angle);

loopBox* buildLoopBox(
    const int start,
    const short* const pair_table,
    const tBaseInformation* baseInformation,
    const double* x,
    const double* y
);

stemBox* buildStemBox(
    const int start,
    const int end,
    const short* const pair_table,
    const double* x,
    const double* y,
    const double bulgeDist
);

loopBox* createLoopBox(
    const double center[2],
    const double radius
);

stemBox* createStemBox(
    const double s[2],
    const double e[2],
    const double sp[2]
);

void getLBoxCenter(const loopBox* box, double c[2]);
void getSBoxCenter(const stemBox* box, double c[2]);

void getBulgeXY(
    const stemBox* stem,
    const int index,
    double *x,
    double *y
);

/**
 * @brief getBulgeCoordinates
 *      Calculates the coordinates of a bulge (given by index) and returns them in pPrev, pThis and pNext.
 *      The point's coordinates are given in world space.
 * @param stem
 * @param index
 * @param pPrev
 *      point prior to the bulges peak
 * @param pThis
 *      point of the bulges peak
 * @param pNext
 *      point successive to the bulges peak
 */
void getBulgeCoordinates(
    const stemBox* stem,
    const int index,
    double pPrev[2],
    double pThis[2],
    double pNext[2]
);

/**
 * @brief getBulgeCoordinatesExtraDistance
 *      Calculates the coordinates of a bulge (given by index) and returns them in pPrev, pThis and pNext.
 *      The point's coordinates are given in world space.
 *      This function adds some extra distance in the stem's b-direction.
 *      Useful for bounding box calculation.
 * @param stem
 * @param index
 * @param extraDistance
 * @param pPrev
 *      point prior to the bulges peak
 * @param pThis
 *      point of the bulges peak
 * @param pNext
 *      point successive to the bulges peak
 */
void getBulgeCoordinatesExtraDistance(
    const stemBox* stem,
    const int index,
    const double extraDistance,
    double pPrev[2],
    double pThis[2],
    double pNext[2]
);

void printLBox(const loopBox* loop);
void printSBox(const stemBox* stem);

#endif
