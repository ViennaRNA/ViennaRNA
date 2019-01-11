#include "ViennaRNA/RNApuzzler/RNAturtle.h"
#include "ViennaRNA/RNApuzzler/postscript/postscriptArcs.h"
#include "ViennaRNA/RNApuzzler/vector_math.h"
#include "ViennaRNA/RNApuzzler/data/cfg_reader.h"
#include "ViennaRNA/utils.h"

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#define printInitialConfig 0

//--------------------------------------------------------------------------------------------------------------------------

/**
 * Handle bases of the exterior loop
 *
 * @returns first paired base
 */
short handleExteriorBases(
        const short *const pair_table,
        short currentBase,
        tBaseInformation* baseInformation,
        const int direction
) {
    // const double UNPAIRED_BEND = 0.; // TODO set good value
    //		printf("unpaired %f %f\n",UNPAIRED_BEND, baseInformation[currentBase].angle );
    short const length = pair_table[0];
    //angles[currentBase] = angles[currentBase-1] + UNPAIRED_BEND; // Bend first
    if (currentBase > 1) {
        //angles[currentBase] = angles[currentBase] + UNPAIRED_BEND; // Bend first
        //angles[currentBase] -= direction * MATH_PI_HALF;
        baseInformation[currentBase].angle = baseInformation[currentBase].angle + direction * MATH_PI_HALF;
        baseInformation[currentBase].baseType = TYPE_EXTERIOR;
    }

    while (currentBase < length && pair_table[currentBase] <= 0) {
        //angles[currentBase+1] = angles[currentBase] + UNPAIRED_BEND;
        baseInformation[currentBase + 1].angle = 0.0;
        baseInformation[currentBase].baseType = TYPE_EXTERIOR;

        currentBase++;
        //		printf("FIT? %f\t%f\n",angles[currentBase-1] + UNPAIRED_BEND,angles[currentBase]);
    }
    // For following stem, since it cannot alter the first two bases
    // which may contain loop exit angles.
    if ((currentBase + 1) <= (length)) {
        //	  printf("Changed direction at %d to %d\n", currentBase, direction);
        //angles[currentBase+1]=angles[currentBase] - direction * MATH_PI_HALF;
        baseInformation[currentBase + 1].angle = direction * MATH_PI_HALF;
        baseInformation[currentBase].baseType = TYPE_EXTERIOR;
    } else {
        baseInformation[currentBase].baseType = TYPE_EXTERIOR;
    }
    return currentBase;
}

// Count no. of base pairs and consecutive pairs on a loop
// i - last base of first stem half
// m - no. of base pairs; n - no. of consecutive pairs

void countLoopPairs(short*const m, short *const n, short i, const short *const pair_table) {
    short const end = pair_table[i++]; // end is 1st base of 2nd stem half
    //		printf( "\tCounting loop pairs from %d to %d.\n", i, end);
    *m = *n = 1; // Base pair of enclosing stem, consecutive pair from last stem base to first loop base
    while (i < end) {
        if (pair_table[i] <= 0 || pair_table[i] < i) {
            (*n)++;
            i++; // Go to consecutive base
        } else // base pair ahead
        {
            (*m)++;
            i = pair_table[i]; // Go to base pair partner
        }
    }

    return;
}

//--------------------------------------------------------------------------------------------------------------------------
void handleStem(
    const short *const pair_table,
    short currentBase,
    const double paired,
    const double unpaired,
    tBaseInformation* const baseInformation,
    int direction
);

//--------------------------------------------------------------------------------------------------------------------------

/*---------------------------------------------------------------------------*/

/**
 * Detect bulges (internal loops with only one unpaired base)
 */
int detectBulge(
    short start,
    const short *const pair_table
) {
    int bulge = 0;
    int end = pair_table[start];
    int iterate = 1;
    int old = 0;
    int i = start + 1;
    int side = 0;

    do {
        //  printf("Iterate: i %d iterate  %d old %d pair_table[i]%d side %d end %d\n", i, iterate, old, pair_table[i], side,end);
        if (pair_table[i] > 0) {
            //There was a basepair before
            if (iterate > 0) {

                //Stem entry
                if (pair_table[i] == old) {
                    side++;
                    i++;
                }//Bulge found
                else {
                    //			printf("start %d pair %d\n",start, pair_table[i]);

                    //It's only a bulge if the paired base is the start base or the bulge is on the other strand
                    if (start == pair_table[i] || pair_table[i] == end - 2) {
                        //				  printf("Bulge found: %d\n", pair_table[i]);
                        bulge = pair_table[i];
                        break;
                    } else {

                        break;
                    }
                }
            }//Found first basepair
            else {
                iterate++;
                old = i;
                i = pair_table[i];
            }
        }//Consecutive Base found
        else {
            //There was a basepair before
            if (iterate > 0) {
                iterate = 0;
                i++;
            }//Just the next consecutive base
            else {
                i++;
            }
        }
        //printf("Iterate: i %d iterate  %d old %d pair_table[i]%d side %d end %d\n", i, iterate, old, pair_table[i], side,end);
    } while (i > start);

    //side =1 bulge on the strand of i
    //side =0 bulge on the opposite strand
    // printf("Bulge? %d %d\n", bulge, side);
    return bulge;
}

// i - last base of first stem half

void handleLoop(
        short i,
        const short *const pair_table,
        const double paired,
        const double unpaired,
        tBaseInformation* baseInformation,
        int direction
) {

    int start = i;
    //    printf("Angle i %d\n",i);

    short const end = pair_table[i]; // end is 1st base of 2nd stem half
    short m, n; // m - no. of base pairs; n - no. of consecutive pairs

    countLoopPairs(&m, &n, i, pair_table);


    int bulge = detectBulge(i, pair_table);
    //	printf("Bulge  found at %d\n", bulge);

    if (bulge > 0 && n - m == 1) {
        // detected bulge with a single unpaired base

        //	printf("Bulge %d %d\n", bulge, n - m);
        int length = (unpaired * (n - m + 1)) / 2;

        double alpha = acos(unpaired / (2 * length));
        //alpha = alpha;

        //	printf("Length: %d\t Angle a:%f\n", length,alpha);
        //	printf("Pair check: %d %d\n", i + 1, pair_table[i + 1]);
        //	printf("Angle: %f\n", angles[i]);
        if (pair_table[i + 1] == 0) {
            //      printf("Loop on this strand\n");
            //angles[i+1] = angles[i+1]  - direction *alpha ;
            baseInformation[i + 1].angle = baseInformation[i + 1].angle + direction*alpha;
            baseInformation[i].baseType = TYPE_BULGE;
            baseInformation[pair_table[i]].baseType = TYPE_BULGE;
            i++;
            //	printf("Angle: %f\n", angles[i]);
            //angles[i+1] = angles[i] + direction *alpha*2;
            baseInformation[i + 1].angle = -direction * alpha * 2;
            baseInformation[i].baseType = TYPE_BULGE;
            i++;
            //	printf("Angle: %f\n", angles[i]);
            //angles[i+1] = angles[i] - direction * alpha;
            baseInformation[i + 1].angle = direction * alpha;
            baseInformation[i].baseType = TYPE_BULGE;
            baseInformation[pair_table[i]].baseType = TYPE_BULGE;

            handleStem(pair_table, i, paired, unpaired, baseInformation, direction);

            i = pair_table[i];

            //angles[i+1] = angles[i+1];

        } else {

            //      printf("Loop on the opposite strand\n");
            //angles[i+1] = angles[i+1];
            baseInformation[i + 1].angle = baseInformation[i + 1].angle + 0.0;
            baseInformation[i].baseType = TYPE_BULGE;
            i++;
            //angles[i+1] = angles[i];
            baseInformation[i + 1].angle = baseInformation[i + 1].angle + 0.0;
            baseInformation[i + 1].baseType = TYPE_BULGE;
            //angles[i+2] = angles[i];
            baseInformation[i + 2].angle = baseInformation[i + 2].angle + 0.0;
            baseInformation[i + 1].baseType = TYPE_BULGE;

            handleStem(pair_table, i, paired, unpaired, baseInformation, direction);

            i = pair_table[i];
            //angles[i+1] = angles[i+1] - direction * alpha;
            baseInformation[i + 1].angle = baseInformation[i + 1].angle + direction * alpha;
            baseInformation[i].baseType = TYPE_BULGE;
            i++;
            //angles[i+1] = angles[i] + direction * alpha*2;
            baseInformation[i + 1].angle = -direction * alpha * 2;
            baseInformation[i].baseType = TYPE_BULGE;
            i++;
            //angles[i+1] = angles[i] - direction *alpha;
            baseInformation[i + 1].angle = direction * alpha;
            baseInformation[i].baseType = TYPE_BULGE;

        }
    } else {

        // loop is not bulge with a single unpaired base

        // ******************************************************************
        // ***************** Loop Drawing Algorithm - Start *****************
        // ******************************************************************


        // variables for config drawing...
        double angle_over_paired;            // alpha      at calculation
        double current_angle;                // phi        at calculation
        double current_bb_angle;             // beta       at calculation
        double current_distance;             // b          at calculation
        double current_delta_ab;             // delta_ab   at calculation
        double current_delta_bb;             // delta_bb   at calculation
        int   current_stem_count;

        //printf("--------------------------------------\n");
        //printf("no of unpaired base %d , no. of stems %d, %d\n", n - m, m, n);

        config *cfg = baseInformation[start].config;
        int currentArc = 0;
        //printf("** using the following config **\n");
        //cfg_printConfig(cfg);
        //printf("** **\n");

        //                           a     , b
        double r = cfg->radius;

        //printf("Radius: %f\n", r);

        //printf("Counted %d paired and %d consecutive pairs.\n", m,n);

        angle_over_paired = 2 * asin(paired / (2 * r));
        //printf("angle_over_paired: %.10f\n", angle_over_paired);

        current_angle = getArcAngle(cfg, currentArc);
        current_bb_angle = ( current_angle - angle_over_paired ) / (cfg->cfgArcs[currentArc]).numberOfArcSegments;
        current_distance = sqrt( 2 * r * r * ( 1 - cos( current_bb_angle ) ) );
        current_delta_ab = 0.5 * ( MATH_PI + angle_over_paired + current_bb_angle );
        current_delta_bb = MATH_PI + current_bb_angle;
        ++currentArc;

        // consider the direction that was given before (e.g. because of a bulge)
        baseInformation[i + 1].angle = baseInformation[i + 1].angle + direction * (MATH_PI - current_delta_ab);

        baseInformation[i].distance = current_distance;

        //printf("##### [%d] angle: %.5f alpha: %.5f beta: %.5f paired: %d unpaired: %.2f d_ab: %.5f d_bb: %.5f\n", config_arc_it, current_angle, angle_over_paired, current_bb_angle, paired, current_distance, current_delta_ab, current_delta_bb);

        current_stem_count = 0;

        //Double Loop Check: Loop is followed directly by a loop
        if (baseInformation[i].baseType == TYPE_LOOP1) {
            baseInformation[i].baseType = TYPE_LOOP2;
        } else {
            baseInformation[i].baseType = TYPE_LOOP1;
        }
        i++;

        while (i < end) {
            if (pair_table[i] <= 0) { // unpaired

                // handle backbone element
                baseInformation[i + 1].angle = -direction * (current_delta_bb - MATH_PI);
                baseInformation[i].distance = current_distance;
                //printf("-> [%d] unpaired -> print backbone (baseInformation[%d]: %f)\n", i, i + 1, angle[i + 1].angle);

                baseInformation[i].baseType = TYPE_LOOP1;
                i++;

            } else if (pair_table[i] > i) { // base pair ahead

                // i is the beginning of a stem and therefore the end of the current arc
                baseInformation[i + 1].angle = direction * (MATH_PI - current_delta_ab);
                //printf("goto stem with angle %3.2f°\n", toDegree(baseInformation[i+1].angle));
                //printf("-> [%d] paired -> continue loop at [%d]\n", i, pair_table[i]);
                current_stem_count++;

                baseInformation[i].baseType = TYPE_LOOP1;

                handleStem(pair_table, i, paired, unpaired, baseInformation, direction); // Recurse multiloop

                i = pair_table[i]; // Go to base pair partner

            } else { // base pair behind


                // i is the end of a stem that returned to the loop
                // recalculation if angles and distances is nessecary because the next arc has been reached
                //printf("-> [%d] end of stem and back to loop\n", i);
                if ( current_stem_count == 1) {
                    //printf("-> new config area ... recalculate current arc [%d]->[%d]\n", config_arc_it, config_arc_it + 1);

                    current_stem_count = 0;
                    current_angle = getArcAngle(cfg, currentArc);
                    current_bb_angle = ( current_angle - angle_over_paired ) / (cfg->cfgArcs[currentArc]).numberOfArcSegments;
                    current_distance = sqrt( 2 * r * r * ( 1 - cos( current_bb_angle ) ) );
                    current_delta_ab = 0.5 * ( MATH_PI + angle_over_paired + current_bb_angle );
                    current_delta_bb = MATH_PI + current_bb_angle;
                    ++currentArc;

                    //printf("##### [%d] angle: %.5f alpha: %.5f beta: %.5f a: %d unpaired: %.2f d_ab: %.5f d_bb: %.5f\n", config_arc_it, current_angle, angle_over_paired, current_bb_angle, paired, current_distance, current_delta_ab, current_delta_bb);
                }

                // consider the direction that was given before (e.g. because of a bulge)
                baseInformation[i + 1].angle = baseInformation[i + 1].angle + direction * (MATH_PI - current_delta_ab);
                //printf("from stem with angle %3.2f°\n", toDegree(baseInformation[i+1].angle));
                baseInformation[i].distance = current_distance;

                baseInformation[i].baseType = TYPE_LOOP1;

                i++;
            }
        }

        if ( (i + 1) <= pair_table[0] ) { baseInformation[i + 1].angle = direction * (MATH_PI - current_delta_ab); }

        baseInformation[i].baseType = TYPE_LOOP1;
        //printf ("Loop done\n");


        // ******************************************************************
        // ****************** Loop Drawing Algorithm - End ******************
        // ******************************************************************
    }

    return;
}

//--------------------------------------------------------------------------------------------------------------------------

void handleStem(
        const short *const pair_table,
        short i,
        const double paired,
        const double unpaired,
        tBaseInformation* const baseInformation,
        int direction
) {
    short end = pair_table[i] + 1; // first base after stem
    //	printf( "\tHandling stem from %d to %d.\n", i, end-1);

    // SKip first two bases of stem, angle already written by
    // either loop (+exit angle) or handle_unpaired
    //i+=2;
    baseInformation[i].baseType = TYPE_STEM;
    i++;

    while (pair_table[i] > 0 && // i is base paired
            (pair_table[i] == end - 1 || // First position of stem, continue!
            pair_table[i] + 1 == pair_table[i - 1] // Detect bulges on opposite strand
            )
            ) {
        //angles[i +1] = angles[i];
        baseInformation[i + 1].angle = 0.0;
        baseInformation[i].baseType = TYPE_STEM;
        i++;

    }
    if (pair_table[i] == end - 1) {
        //        printf("Kaboom");
    } else {

        //  printf("Bulge or not?: %d %d\n",pair_table[i]+1 , pair_table[i-1]);
        //      printf("Loop at %d\n",i);
        handleLoop(--i, pair_table, paired, unpaired, baseInformation, direction); // i - last base of first stem half
    }
    i = pair_table[i]; // set i as base pairing partner of last stem base
    baseInformation[i].baseType = TYPE_STEM;
    //const double newAngle = angles[pair_table[i]]+MATH_PI; // rotate old angle by 180 deg
    i++;
    while (i < end && i < pair_table[0]) {
        //angles[i++ +1] = newAngle;
        //angles[i +1] = angles[i];
        baseInformation[i].baseType = TYPE_STEM;
        i++;
    }

    return;
}

//------------------------------------------------------------------------------

/**
 Compute angle angles for all loops
 */
void computeAffineCoordinates(
        short const * const pair_table,
        const double paired,
        const double unpaired,
        tBaseInformation* const baseInformation
) {
    /// Initialization
    short const length = pair_table[0];
    short currentBase = 1;
    const int direction = -1;
    baseInformation[0].angle = 0.0;

    // angles[0]=2*MATH_PI;//0.75f*MATH_PI;	// TODO use to rotate entire drawing! //SL wir fangen mit 90 grad an
    // for( i=1; i<length; i++)
    // {
    // 	angles[i] = i*(MATH_PI/150.f);
    // 	//printf( "angle %d: %f\t", i, angles[i]);
    // }

    /// For first stem, since handle_stem cannot alter the first two bases
    /// which may contain loop exit angles (in other stems).
    if (2 <= length) {
        //angles[1]=angles[0];
        //angles[2]=angles[1];
        baseInformation[1].angle = baseInformation[0].angle;
        baseInformation[2].angle = baseInformation[1].angle;
    }

    //  int idx=0;
    /// ABRRUCH NEU DEFINIEREN
    int dangle_count = 0;
    while (currentBase < length) {
        //		printf( "Iterating get_angles while, currentBase=%d\n", currentBase);
        if (pair_table[currentBase] <= 0) {
            // currentBase unpaired
            //printf("dangle %d\n", currentBase);

            if (currentBase > 1) {
                baseInformation[currentBase - 1].baseType = TYPE_EXTERIOR;
            }
            currentBase = handleExteriorBases(pair_table, currentBase, baseInformation, direction);
            // returns first paired base currentBase
            dangle_count++;
            //printf("Dangle End %d\n", currentBase);
        }

        //		printf("first paired base %d\n", currentBase);
        if (currentBase < length) {
            // it was not the dangling end
            if (pair_table[currentBase] - pair_table[currentBase - 1] != 1
                && pair_table[currentBase] != 0
                && pair_table[currentBase - 1] != 0) {
                if (currentBase == 1) {
                    if (dangle_count < 1) {
                        //                                        printf("No Dangling End\n");
                        baseInformation[0].angle = baseInformation[1].angle = baseInformation[2].angle = -MATH_PI_HALF;
                        baseInformation[currentBase].baseType = TYPE_EXTERIOR;
                    }
                    handleStem(pair_table, currentBase, paired, unpaired, baseInformation, direction);
                    currentBase = pair_table[currentBase] + 1;

                    //Check if there is a minor dangling end at the end
                    if (currentBase == length) {
                        baseInformation[currentBase - 1].baseType = TYPE_EXTERIOR;
                        baseInformation[currentBase].baseType = TYPE_EXTERIOR;
                        printf("Found Mini Dangling End\n");
                        baseInformation[currentBase].angle = -MATH_PI_HALF;
                    }
                    continue;
                } else {
                    printf("Mini Exterior Loop found\n");

                    //		printf("Angle_ext: %f\n", angles[currentBase-1]);
                    //		printf("Angle_ext: %f\n", angles[currentBase]);
                    //		printf("Angle_ext: %f\n", angles[currentBase+1]);
                    //		printf("Angle_ext: %f\n", angles[currentBase+2]);

                    baseInformation[currentBase].angle = baseInformation[currentBase].angle + direction*MATH_PI_HALF;
                    baseInformation[currentBase + 1].distance = unpaired;
                    baseInformation[currentBase - 1].baseType = TYPE_EXTERIOR;
                    baseInformation[currentBase + 1].angle = baseInformation[currentBase + 1].angle + direction*MATH_PI_HALF;
                    baseInformation[currentBase].baseType = TYPE_EXTERIOR;

                    //angles[currentBase+1] = angles[currentBase] -3.0*MATH_PI_HALF;
                    //angles[currentBase+2] = angles[currentBase+1] -MATH_PI_HALF;

                    //printf("Angle_ext: %f\n", angles[currentBase-1]);
                    //printf("Angle_ext: %f\n", angles[currentBase]);
                    //printf("Angle_ext: %f\n", angles[currentBase+1]);
                    //printf("Angle_ext: %f\n", angles[currentBase+2]);
                    //		printf("Changed direction at %d to %d\n", currentBase, direction);
                    dangle_count++;
                }
            }

            handleStem(pair_table, currentBase, paired, unpaired, baseInformation, direction);

            currentBase = pair_table[currentBase] + 1; // currentBase is next base after stem
            //printf("I %d\n", currentBase);

            if (currentBase == length) {
                baseInformation[currentBase - 1].baseType = TYPE_EXTERIOR;
                currentBase = handleExteriorBases(pair_table, currentBase, baseInformation, direction);
            }
        }
    }
    baseInformation[length].baseType = TYPE_EXTERIOR;
    //for(currentBase =0; currentBase< length;currentBase++)
    //{
    //  	printf("angle %d  %f\n",currentBase,baseInformation[currentBase+1].angle);
    //}
}

//------------------------------------------------------------------------------

/**
 Calculate the coordinates for the drawing with the given angle angles
 */
void affineToCartesianCoordinates(
    tBaseInformation* const baseInformation,
    unsigned short const length,
    double * const x,
    double * const y
) {
    if (length < 1) {
        return;
    }

    double angle = 0.0;
    x[0] = y[0] = EXTERIOR_Y;
    for (int i = 1; i < length; i++) {
        angle = angle - baseInformation[i + 1].angle;

        x[i] = x[i - 1] + baseInformation[i].distance * cos(angle);
        y[i] = y[i - 1] + baseInformation[i].distance * sin(angle);

        // }
//        if (baseInformation[i].baseType == TYPE_EXTERIOR
//            && baseInformation[i + 1].baseType == TYPE_EXTERIOR) {
//            x[i] = x[i - 1] + baseInformation[i].distance * cos(angle);
//            y[i] = y[i - 1] + baseInformation[i].distance * sin(angle);
//        }

    }
    return;
}

//------------------------------------------------------------------------------

int layout_RNAturtle(
        short const * const pair_table,
        float *x,
        float *y,
        double *arc_coords
) {
    const char *fnName = "layout_RNAturtle";

    const short drawArcs = 1;
    const short paired = 35.0;
    const short unpaired = 25.0;

    const int length = pair_table[0];
    //printf("RNA length: %d\n", length);

    /// turtle base information
    tBaseInformation* baseInformation = vrna_alloc((length + 1) * sizeof(tBaseInformation));
    for (int i = 0; i <= length; i++) {
        baseInformation[i].baseType = TYPE_BASE_NONE;
        baseInformation[i].distance = unpaired;
        baseInformation[i].angle = 0.0;
        baseInformation[i].config = NULL;
    }

    /// generate default configuration for each loop
    cfgGenerateConfig(pair_table, baseInformation, unpaired, paired);

    if (printInitialConfig) {
        printf("** print initial config **\n");
        for (int i = 0; i <= length; i++) {
            if (baseInformation[i].config != NULL) {
                cfgPrintConfig(baseInformation[i].config);
            }
        }
    }

    /// compute loop angles
    computeAffineCoordinates(pair_table, paired, unpaired, baseInformation);

    /// transform affine coordinates into cartesian coordinates
    double* myX = (double*) vrna_alloc(length * sizeof(double));
    double* myY = (double*) vrna_alloc(length * sizeof(double));
    affineToCartesianCoordinates(baseInformation, length, myX, myY);

    if (drawArcs) {
        // compute postscript arcs instead of lines
        computeAnglesAndCentersForPS(pair_table, myX, myY, baseInformation, arc_coords);
    }

    for (int i = 0; i < length; i++) {
        x[i] = myX[i];
        y[i] = myY[i];
    }

    free(myX);
    free(myY);

    // free turtle struct
    free(baseInformation);
//    printf("baseInformation\n");

    return length;
}

