#ifndef INTERSECTLEVELBOUNDINGBOXES_H
#define INTERSECTLEVELBOUNDINGBOXES_H

#include "ViennaRNA/plotting/RNApuzzler/definitions.h"
#include "ViennaRNA/plotting/RNApuzzler/data/boundingBoxes.h"

/**
 * @brief intersectStemStem
 * Checks for intersection of two given stems.
 * @param stem1
 * @param stem2
 * @return 1 if intersecting, 0 otherwise
 */
PRIVATE short intersectStemStem(const stemBox *stem1,
                        const stemBox *stem2);


/**
 * @brief intersectLoopLoop
 * Checks for intersection of two given loops.
 * @param loop1
 * @param loop2
 * @return 1 if intersecting, 0 otherwise
 */
PRIVATE short intersectLoopLoop(const loopBox *loop1,
                        const loopBox *loop2);


/**
 * @brief intersectStemLoop
 * Check for intersection of given stem and given loop.
 * @param stem
 * @param loop
 * @return 1 if intersecting, 0 otherwise
 */
PRIVATE short intersectStemLoop(const stemBox *stem,
                        const loopBox *loop);


/**
 * @brief intersectLoopBulges
 * Check for intersection of given loop and given stem's bulges.
 * @param loop
 * @param stem
 * @param bulge
 * @return 1 if intersecting, 0 otherwise
 */
PRIVATE short intersectLoopBulges(const loopBox *loop,
                          const stemBox *stem,
                          int           *bulge);


/**
 * @brief intersectBulgesBulges
 * @return
 */
PRIVATE short intersectBulgesBulges(const stemBox *stem1,
                            const stemBox *stem2,
                            int           *bulge1,
                            int           *bulge2);


/**
 * @brief intersectStemBulges
 * @param stem1
 * @param stem2
 * @param bulge2
 * @return
 */
PRIVATE short intersectStemBulges(const stemBox *stem1,
                          const stemBox *stem2,
                          int           *bulge2);


#include "intersectLevelBoundingBoxes.c"


#endif
