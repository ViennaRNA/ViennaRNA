#include "ViennaRNA/RNApuzzler/definitions.h"
#include "ViennaRNA/RNApuzzler/vector_math.h"

#include "ViennaRNA/utils.h"

#include <stdlib.h>
#include <math.h>

puzzlerOptions* createPuzzlerOptions() {
    puzzlerOptions* puzzler = (puzzlerOptions*) vrna_alloc(sizeof(puzzlerOptions));

    // drawing behavior
    puzzler->drawArcs = 1;
    puzzler->paired = 35.0;
    puzzler->unpaired = 25.0;

    // intersection resolution behavior
    puzzler->checkExteriorIntersections = 1;
    puzzler->checkSiblingIntersections = 1;
    puzzler->checkAncestorIntersections = 1;
    puzzler->optimize = 1;

    // import behavior - unused for now
    puzzler->config = NULL;

    // other stuff
    puzzler->filename = NULL;

    puzzler->numberOfChangesAppliedToConfig = 0;
    puzzler->psNumber = 0;

    return puzzler;
}

void destroyPuzzlerOptions(puzzlerOptions* puzzler) {
    free(puzzler);
}

void bubblesort(
        const int numValues,
        const double* const valuesLevel1,
        const double* const valuesLevel2,
        int* const indices
) {
    for (int i = 0; i < numValues; i++) {
        indices[i] = i;
    }

    double thisValue = 0.0;
    double nextValue = 0.0;
    short swap = 0;
    for (int i = 0; i < numValues-1; i++) {
        for (int j = 0; j < numValues-i-1; j++) {
            thisValue = valuesLevel1[indices[j+0]];
            nextValue = valuesLevel1[indices[j+1]];
            swap = 0;
            if (nextValue - thisValue > epsilon7) {
                swap = 1;
            } else if (fabs(nextValue - thisValue) < epsilon7) {
                thisValue = valuesLevel2[indices[j+0]];
                nextValue = valuesLevel2[indices[j+1]];
                if (nextValue - thisValue > epsilon7) {
                    swap = 1;
                }
            }

            if (swap) {
                int tmp = indices[j+0];
                indices[j+0] = indices[j+1];
                indices[j+1] = tmp;
            }
        }
    }
}

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
double distanceToAngle(const double radius, const double distance) {
    return (2.0 * asin(distance / (2.0 * radius)));
}

double angleToDistance(const double radius, const double angle) {
    return (2.0 * radius * sin(angle / 2.0));
}
