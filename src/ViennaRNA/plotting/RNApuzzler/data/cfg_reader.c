#include "ViennaRNA/RNApuzzler/data/cfg_reader.h"
#include "ViennaRNA/RNApuzzler/definitions.h"
#include "ViennaRNA/RNApuzzler/vector_math.h"

#include "ViennaRNA/utils.h"

#include <stdlib.h>
#include <math.h>

/*--------------------------------------------------------------------------*/
/*---   create, copy, and free config   ------------------------------------*/
/*--------------------------------------------------------------------------*/

/**
 * @brief cfgCreateConfig
 *      - constructor-like method for creating a config
 * @param radius
 *      - radius used for drawing that loop
 * @return
 *      - an initialized config struct
 */
config * cfgCreateConfig(
    const double radius
) {
    config * cfg = (config*) vrna_alloc(1*sizeof(config));

    cfg->radius = radius;
    cfg->minRadius = radius;
    cfg->defaultRadius = radius;

    cfg->cfgArcs = NULL;
    cfg->numberOfArcs = 0;

    return cfg;
}

/**
 * @brief cfgCreateConfigArc
 *      - constructor-like method for adding a new config entry to a given config
 * @param angle
 *      - angle (radiant) that can be found between stems 'from' and 'to' later on
 * @param numberOfArcSegments
 *      - number of arc segments between stems 'from' and 'to'
 * @return
 *      - an initialized configArc struct
 */
configArc cfgCreateConfigArc(
    const double angle,
    const int numberOfArcSegments
) {
    configArc newConfigArc;

    newConfigArc.numberOfArcSegments = numberOfArcSegments;

    newConfigArc.arcAngle = angle;

    return newConfigArc;
}

config *cfgCloneConfig(const config *cfg) {
    config *clonedCfg = (config *) vrna_alloc(sizeof(config));
    clonedCfg->radius = cfg->radius;
    clonedCfg->minRadius = cfg->minRadius;
    clonedCfg->defaultRadius = cfg->defaultRadius;
    clonedCfg->numberOfArcs = cfg->numberOfArcs;

    const int numberOfArcs = cfg->numberOfArcs;
    clonedCfg->cfgArcs = (configArc *) vrna_alloc(numberOfArcs * sizeof(configArc));
    for (int currentArc = 0; currentArc < numberOfArcs; ++currentArc) {
        clonedCfg->cfgArcs[currentArc].numberOfArcSegments = cfg->cfgArcs[currentArc].numberOfArcSegments;
        clonedCfg->cfgArcs[currentArc].arcAngle = cfg->cfgArcs[currentArc].arcAngle;
    }

    return clonedCfg;
}

void cfgFreeConfig(config* cfg) {
    free(cfg->cfgArcs);
    free(cfg);
}

// documentation at header file
void cfgPrintConfig(config *cfg) {
    short useDegree = 1;

    for (int currentArc = 0; currentArc < cfg->numberOfArcs; ++currentArc) {
        double angle = getArcAngle(cfg, currentArc);
        if (useDegree) { angle = toDegree(angle); }
        printf("\t%6.2f%s %3.2f\n",
                angle,
                (useDegree ? "°" : ""),
                cfg->radius);
    }
}

/*--------------------------------------------------------------------------*/
/*---   access to config elements   ----------------------------------------*/
/*--------------------------------------------------------------------------*/

double getArcAngle(
        const config *cfg,
        const int currentArc
) {
    return cfg->cfgArcs[currentArc].arcAngle;
}

double getArcAngleDegree(
        const config *cfg,
        const int currentArc
) {
    return toDegree(getArcAngle(cfg, currentArc));
}

/*--------------------------------------------------------------------------*/
/*---   radius computation   -----------------------------------------------*/
/*--------------------------------------------------------------------------*/

/**
 * Approximate the radius of a circle required to draw m base pairs
 * with a distance of 'a' to each other, and n (unpaired) consecutive
 * nucleotides with a distance of b over a specified angle.
 *
 * Uses a Newton iteration to approximate solution
 *
 * @param a paired
 * @param b unpaired
 * @param m #stems
 * @param n #backbones
 * @param angle angle
 */
double approximateConfigArcRadius(
        const double a,
        const double b,
        const short m,
        const short n,
        double angle
) {

    const short MAX_ITERATIONS = 1000;

    /// calculation:
    ///
    /// be s the length of a line at the circle (paired or unpaired / a or b)
    /// the angle over such a single line be alpha
    ///     alpha = angle / ( m + n )
    ///
    /// for such a single line the following equation holds (where r is the radius of the circle)
    ///     sin( alpha / 2 ) = ( s / 2 ) / r
    ///     r = ( s / 2 ) / sin( alpha / 2 )
    ///     r = ( s / 2 ) / sin( ( angle / ( m + n ) ) / 2 )
    ///
    /// now we replace s with a or b to get the upper or lower bound for the radius interval
    double lowerBound = ( b / 2 ) / sin( ( angle / ( m + n ) ) / 2 );
    double upperBound = ( a / 2 ) / sin( ( angle / ( m + n ) ) / 2 );
//    printf("new: lower %f upper %f\n", lowerBound, upperBound);
    double rtn = 0.5 * (lowerBound + upperBound);

    /// there is a minimum valid radius!
    /// if rtn is smaller than 0.5*a or 0.5*b than the result will become nan
    rtn = fmax(rtn, 0.5 * a);
    rtn = fmax(rtn, 0.5 * b);

    //printf("lower %f upper %f %f\n", lowerBound, upperBound, rtn);
    //printf("Angle %f\n", angle);

//    short isERROR = 0;

    int j = 0;
    for (j = 0; j < MAX_ITERATIONS; j++) {
        double dx = 2 * (m * asin(a / (2 * rtn)) + n * asin(b / (2 * rtn)) - (angle / 2))
                      / (-(a * m / (rtn * sqrt(rtn * rtn - a * a / 4)) + b * n / (rtn * sqrt(rtn * rtn - b * b / 4))));
        rtn -= dx;
//        if ((lowerBound - rtn) * (rtn - upperBound) < 0.0) {
//            if (isERROR) {
//                // print the prior iteration's state if the error still exists
//                printf("[WARNING] [GET RADIUS] jumped out of the brackets ( lower:%f upper:%f rtn:%f | a:%3.2f b:%3.2f m:%d n:%d )\n", lowerBound, upperBound, rtn, a, b, m, n);
//            }
//
//            isERROR = 1;
//        } else {
//            isERROR = 0;
//        }
        if (fabs(dx) < epsilon3) {
            break;
        }
        //printf("rtn: %f (%d)\n", rtn, j+1);
    }

    if (rtn < lowerBound) {
        printf("[WARNING] [GET RADIUS] result too small: %12.8lf < %12.8lf -> reset\n", rtn, lowerBound);
        rtn = lowerBound;
    } else if (rtn > upperBound) {
        printf("[WARNING] [GET RADIUS] result too large: %12.8lf > %12.8lf -> reset\n", rtn, upperBound);
        rtn = upperBound;
    }

    if (j >= MAX_ITERATIONS) {
        printf("[WARNING] [GET RADIUS] iterarion limit reached (%d)\n", MAX_ITERATIONS);
    }

    //printf("New radius is %f, iterations: %d\n", rtn, j);
    return rtn;
}

/**
*/
double approximateConfigRadius(
    const config *cfg,
    const double unpaired,
    const double paired
) {
    // calculate a fitting radius for each arc without compressing or stretching arc segments
    // return the maximum of those values as the radius fitting for the loop
    //printf("=========================================\n");
    //printf("radius calculation with configArc struct:\n");
    double r = 0;
    for (int currentArc = 0; currentArc < cfg->numberOfArcs; ++currentArc) {
        int stems = 1;
        int numberOfArcSegments = (cfg->cfgArcs[currentArc]).numberOfArcSegments;
        double angle = getArcAngle(cfg, currentArc);

        double tempR = approximateConfigArcRadius(paired, unpaired, stems, numberOfArcSegments, angle);
        //printf("m: %d n: %d angle: %.10f -> radius: %.10f\n", stems, numberOfArcSegments, angle, tempR);

        if (tempR > r) {
            r = tempR;
        }
    }
    //printf("radius config new: %f\n", r);
    return r;
}

//--------------------------------------------------------------------------------------------------------------------------

/**
 * @brief cfgGenerateDefaultConfig
 *      - generates a config that resembles a drawing without any
 *        constraints given by config input for the given loop
 * @param pair_table
 *      - the RNA's pairing information
 * @param start
 *      - index of the loop's first base
 * @param unpaired
 *      - default length of backbone lines
 * @param paired
 *      - default distance between paired bases
 * @param radius
 *      - radius for the given loop
 * @return
 *      - complete config for that loop
 */
config * cfgGenerateDefaultConfig(
    const short *const pair_table,
    const int  start,
    const int  unpaired,
    const int  paired,
    const double radius
) {
    /// create loop configuration
    config    *cfg = cfgCreateConfig(radius);

    /// compute angles for paired and unpaired bases
    double anglePaired   = 2 * asin(paired   / (2 * radius));    // angle over paired
    double angleUnpaired = 2 * asin(unpaired / (2 * radius));    // angle over unpaired

    /// initialize values for first arc
    int arcUnpaired = 0;
    double angleArc; // alpha + numBackbones * beta

    /// pointer to first arc
    //configArc **currentArc = NULL;
    //currentArc = &(cfg->first);

    /// start with first base after parent stem
    int i = start + 1;
    while (i <= pair_table[start]) {
        /// until last base at parent stem
        if (pair_table[i] == 0) {
            /// arc
            i++;
        } else {
            /// increment number of arcs
            ++(cfg->numberOfArcs);

            if (i != pair_table[start]) {
                /// skip subtree at stem
                i = pair_table[i] + 1;
            } else {
                /// parent stem -> finish
                break;
            }
        }
    }

    cfg->cfgArcs = (configArc *) vrna_alloc(cfg->numberOfArcs * sizeof(configArc));

    /// start with first base after parent stem
    i = start + 1;
    int currentArc = 0;
    int numberOfArcSegments = 0;
    while (i <= pair_table[start]) {
        /// until last base at parent stem
        if (pair_table[i] == 0) {
            /// arc
            arcUnpaired++;
            i++;
        } else {
            /// stem: create arc
            angleArc = anglePaired + (arcUnpaired + 1) * angleUnpaired;
            numberOfArcSegments = arcUnpaired + 1;
            cfg->cfgArcs[currentArc] = cfgCreateConfigArc(angleArc, numberOfArcSegments);
            ++currentArc;
            //currentArc = &((*currentArc)->next);
            //++(cfg->numberOfArcs);

            if (i != pair_table[start]) {
                /// initialize values for next arc
                arcUnpaired = 0;

                /// skip subtree at stem
                i = pair_table[i] + 1;
            } else {
                /// parent stem -> finish
                break;
            }
        }
    }
    /* DEPRECATED
    int i = start + 1;
    while (i <= pair_table[start]) {
        /// until last base at parent stem
        if (pair_table[i] == 0) {
            /// arc
            arcUnpaired++;
            i++;
        } else {
            /// stem: create arc
            angleArc = anglePaired + (arcUnpaired + 1) * angleUnpaired;
            numberOfArcSegments = arcUnpaired + 1;
            *currentArc = cfgCreateConfigArc(angleArc, numberOfArcSegments);
            currentArc = &((*currentArc)->next);
            ++(cfg->numberOfArcs);

            if (i != pairTable[start]) {
                /// initialize values for next arc
                arcUnpaired = 0;

                /// skip subtree at stem
                i = pair_table[i] + 1;
            } else {
                /// parent stem -> finish
                break;
            }
        }
    }
    */

    return cfg;
}

void cfgGenHandleStem(
    int baseNr,
    const short *const pair_table,
    tBaseInformation *baseInformation,
    const double unpaired,
    const double paired
);

/**
 * @brief cfgGenHandleLoop
 *      - recursively iterates through the RNA and generates default configurations.
 *        Alternates with corresponding handleStem method.
 * @param baseNr
 *      - index of the loop's first base
 * @param pair_table
 *      - the RNA's pairing information
 * @param baseInformation
 *      - array of tBaseInformation annotations (to save config)
 */
void cfgGenHandleLoop(
    int baseNr,
    const short *const pair_table,
    tBaseInformation *baseInformation,
    const double unpaired,
    const double paired
) {
    int start = baseNr;
    int end = pair_table[baseNr];

    int unpairedCount = 0;
    int stemCount = 1;

    // count stems and unpaired bases to use for bulge detection
    int i = start + 1;
    while ( i < end ) {
        if (pair_table[i] == 0) {
            // unpaired base
            unpairedCount++;
            i++;
        } else if (pair_table[i] > i) {
            // found new stem
            stemCount++;
            i = pair_table[i];
        } else {
            // returned from stem
            i++;
        }
    }

    short isBulge = (stemCount == 2 && unpairedCount == 1);
    if (isBulge) {
        if (pair_table[start + 1] == 0) {
            // unpaired on left strand
            cfgGenHandleStem(start + 2, pair_table, baseInformation, unpaired, paired);
        } else {
            // unpaired on the right strand
            cfgGenHandleStem(start + 1, pair_table, baseInformation, unpaired, paired);
        }
    } else {
        int m = stemCount;                     // compare RNApuzzler.c -> f_handle_loop
        int n = unpairedCount + stemCount;    // compare RNApuzzler.c -> f_handle_loop
        double defaultRadius = approximateConfigArcRadius(paired, unpaired, m, n, MATH_TWO_PI);
        config *cfgLoop = cfgGenerateDefaultConfig(pair_table, start, unpaired, paired, defaultRadius);
        baseInformation[start].config = cfgLoop;

        int i = start + 1;
        while ( i < end ) {
            if (pair_table[i] == 0) {
                // unpaired base
                i++;
            } else if (pair_table[i] > i) {
                // found new stem
                cfgGenHandleStem(i, pair_table, baseInformation, unpaired, paired);
                i = pair_table[i];
            } else {
                // returned from stem
                i++;
            }
        }
    }
}

/**
 * @brief cfgGenHandleStem
 *      - recursively iterates through the RNA and generates default configurations.
 *        Alternates with corresponding handleLoop method.
 * @param baseNr
 *      - index of the stem's first base
 * @param pair_table
 *      - the RNA's pairing information
 * @param baseInformation
 *      - array of tBaseInformation annotations (to save config)
 */
void cfgGenHandleStem(
    int baseNr,
    const short *const pair_table,
    tBaseInformation *baseInformation,
    const double unpaired,
    const double paired
) {
    // does nothing but iterating over the stem and calling cfgGenHandleLoop as soon as a loop is found

    short continueStem = 1;
    int i = baseNr;
    while ( continueStem ) {
        if (pair_table[i+1] == pair_table[i] - 1) {
            i++;
        } else {
            // found unpaired above stem
            cfgGenHandleLoop(i, pair_table, baseInformation, unpaired, paired);
            continueStem = 0;
        }
    }
}

// documentation at header file
void cfgGenerateConfig(
    const short *const pair_table,
    tBaseInformation *baseInformation,
    const double unpaired,
    const double paired
) {
    //printf("config generation: started\n");

    // global distance values
    //distUnpaired = unpaired;
    //distPaired = paired;

    // iterate over RNA
    // for each loop generate a default config

    int length = pair_table[0];
    int i = 1;
    while (i < length) {
        if (pair_table[i] == 0) {
            // unpaired at exterior loop
            i++;
        } else if (pair_table[i] > i) {
            // found stem
            cfgGenHandleStem(i, pair_table, baseInformation, unpaired, paired);
            i = pair_table[i];
        } else {
            // returned from stem
            i++;
        }
    }

    //printf("config generation: finished\n");
}

/*--------------------------------------------------------------------------*/
/*---   set and update config elements   -----------------------------------*/
/*--------------------------------------------------------------------------*/

// documentation at header file
/**
 * @brief cfgSetRadius
 *      - changes the value of radius for a config to the given value
 * @param config
 *      - config that is being altered
 * @param radius
 *      - new radius
 */
void cfgSetRadius(
    config *cfg,
    const double radius
) {
    cfg->radius = radius;
}

/**
 * @brief cfgUpdateMinRadius
 *      - updates the minimum possible radius for the given config
 * @param config
 * @param unpaired
 * @param paired
 */
void cfgUpdateMinRadius(
    config* cfg,
    const double unpaired,
    const double paired
) {
    double minRadius = approximateConfigRadius(cfg, unpaired, paired);
    cfg->minRadius = minRadius;
}

double cfgApplyChanges(
    config* cfg,
    const char loopName,
    const double* deltaCfg,
    const double radiusNew,
    const puzzlerOptions* puzzler
) {
    /// - start with adjusting config angles; if applicable
    if (deltaCfg != NULL) {
        for (int currentArc = 0; currentArc < cfg->numberOfArcs; currentArc++) {
            (cfg->cfgArcs[currentArc]).arcAngle += deltaCfg[currentArc];
        }
    }

    /// - then, adjust config radius
    double oldRadius = cfg->radius;
    double newRadius = -1.0;
    if (radiusNew > 0.0) {
        /// in case the input is a positive value
        /// we set the minimum of valid and input as new radius
        cfgUpdateMinRadius(cfg, puzzler->unpaired, puzzler->paired);
        newRadius = fmax(radiusNew, cfg->minRadius);
        cfgSetRadius(cfg, newRadius);
    } else
    if (radiusNew == 0.0) {
        /// set the minRadius as new value
        /// (this allows to shrink a loop)
        cfgUpdateMinRadius(cfg, puzzler->unpaired, puzzler->paired);
        newRadius = cfg->minRadius;
        cfgSetRadius(cfg, newRadius);
    } else
    if (radiusNew == -1.0) {
        /// set the minRadius as new value
        /// (this forbidds to shrink a loop)
        cfgUpdateMinRadius(cfg, puzzler->unpaired, puzzler->paired);
        if (cfg->minRadius - epsilon0 > oldRadius) {
            newRadius = cfg->minRadius;
        } else {
//            double defaultIncrease = epsilonFix;
//            printf("[ LOG ] increase radius by %5.2f\n", defaultIncrease);
//            newRadius = oldRadius + defaultIncrease;
            double defaultIncrease = 1.05;
            //printf("[ LOG ] increase radius by x%5.2f\n", defaultIncrease);
            newRadius = oldRadius * defaultIncrease;
        }
        cfgSetRadius(cfg, newRadius);
    } else {
        /// all unhandled inputs result in errors
        printf("[ERROR] set %c's new radius to -1.0 because of invalid input %10.8lf.\n", loopName, radiusNew);
        newRadius = -1.0;
    }

    return newRadius;
}

short cfgIsValid(
    config* cfg,
    const double* deltaCfg
) {
    if (deltaCfg == NULL) {
        return 0;
    }

    double sumAngles = 0.0;
    short validSingleAngles = 1;
    for (int currentArc = 0; currentArc < cfg->numberOfArcs; currentArc++) {
        double angle = getArcAngle(cfg, currentArc) + deltaCfg[currentArc];
        sumAngles += angle;

        short validAngle = 0.0 < angle && angle < MATH_TWO_PI;
        validSingleAngles = validSingleAngles && validAngle;
    }

    short validSumAngles = (fabs(sumAngles - MATH_TWO_PI) < epsilon3);

    return validSingleAngles && validSumAngles;
}

/*
void DEPRECATED_getUnpairedAngles(
    config* cfg,
    const double pairedAngle,
    double* const unpairedAngles
) {
    configArc* cfgArc = cfg->first;
    int i = 0;

    while (cfgArc != NULL) {
        double cfgAngle = toDegree(cfgArc->angle);
        double sumAlphas = cfgAngle - pairedAngle;
        int numberOfArcSegments = cfgArc->numberOfArcSegments;
        double alpha = sumAlphas / numberOfArcSegments;
        unpairedAngles[i] = alpha;

        cfgArc = cfgArc->next;
        i++;
    }
}
*/
