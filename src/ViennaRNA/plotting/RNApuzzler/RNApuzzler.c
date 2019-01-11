#include "ViennaRNA/RNApuzzler/RNApuzzler.h"
#include "ViennaRNA/RNApuzzler/RNAturtle.h"
#include "ViennaRNA/RNApuzzler/postscript/postscriptArcs.h"
#include "ViennaRNA/RNApuzzler/vector_math.h"
#include "ViennaRNA/RNApuzzler/definitions.h"
#include "ViennaRNA/RNApuzzler/data/cfg_reader.h"
#include "ViennaRNA/RNApuzzler/data/configtree.h"
#include "ViennaRNA/RNApuzzler/intersectLevel/intersectLevelLines.h"
#include "ViennaRNA/RNApuzzler/resolveIntersections/resolveIntersections.h"
#include "ViennaRNA/RNApuzzler/resolveIntersections/resolveExteriorChildIntersections.h"
#include "ViennaRNA/RNApuzzler/dataTypes/boundingBoxes_struct.h"
#include "ViennaRNA/RNApuzzler/data/boundingBoxes.h"

#include "ViennaRNA/RNApuzzler/output/output.h"

#include "ViennaRNA/RNApuzzler/output/GeoGebra_output.h"
#include "ViennaRNA/RNApuzzler/output/configtree_debug.h"

#include "ViennaRNA/utils.h"

#include <stdlib.h>
#include <math.h>

#define printInitialConfig 0

//----------------------------------------------------------------------------

short checkRemainingIntersections(
    double *x,
    double *y,
    double *arcCoords,
    const short printDetails,
    const tBaseInformation* baseInformation,
    const int length
) {
    char *fnName = "checkRemainingIntersections";

    if (printDetails) {
        printf("\n\n--- checking for remaining intersections ---\n");
    }

    const short skipExterior = 0;

    short intersect = 0;
    double arc_i[6];
    short i_is_arc = 0;
    double arc_j[6];
    short j_is_arc = 0;

    for (int i = 3; i < length; i++) {
        arc_i[0] = arcCoords[6*i+0];
        arc_i[1] = arcCoords[6*i+1];
        arc_i[2] = arcCoords[6*i+2];
        arc_i[3] = arcCoords[6*i+3];
        arc_i[4] = arcCoords[6*i+4];
        arc_i[5] = arcCoords[6*i+5];

        i_is_arc = (arc_i[0] != -1);
        double i0[2] = { x[i-1], y[i-1] };
        double i1[2] = { x[i-0], y[i-0] };

        if (skipExterior
            && ((i0[1] <= EXTERIOR_Y) || (i1[1] <= EXTERIOR_Y))) {
//            printf("skip i:%3d\n", i);
            continue;
        }
        short intersectIExterior = 0;
        if (baseInformation[i+0].baseType == TYPE_EXTERIOR
            && baseInformation[i+1].baseType == TYPE_EXTERIOR) {
            if (i_is_arc) {
                /// exterior line
                double xmin = fmin(i0[0], i1[0]);
                double xmax = fmax(i0[0], i1[0]);
                double p1[2] = { xmin, EXTERIOR_Y };
                double p2[2] = { xmax, EXTERIOR_Y };

                intersectIExterior = intersectLineArc( p1, p2, arc_i);
            } else {
                intersectIExterior = ((i0[1] <= EXTERIOR_Y) != (i1[1] <= EXTERIOR_Y));
            }
        }
        intersect = intersect || intersectIExterior;

        for (int j = 1; j < i-1; j++) {
//            if (j != 135 && j != 134) { continue; }

//            if (print_details) { printf("------\n%d-%d vs. %d-%d\n", i-1, i, j-1, j); }

            arc_j[0] = arcCoords[6*j+0];
            arc_j[1] = arcCoords[6*j+1];
            arc_j[2] = arcCoords[6*j+2];
            arc_j[3] = arcCoords[6*j+3];
            arc_j[4] = arcCoords[6*j+4];
            arc_j[5] = arcCoords[6*j+5];
            j_is_arc = (arc_j[0] != -1);
            double j0[2] = { x[j-1], y[j-1] };
            double j1[2] = { x[j-0], y[j-0] };

            if (skipExterior && ((j0[1] <= EXTERIOR_Y) || (j1[1] <= EXTERIOR_Y))) {
//                printf("skip i:%3d j:%3d\n", i, j);
                continue;
            }

            short intersect_ij = 0;

            if ( i_is_arc &&  j_is_arc) {
                if (arc_i[0] == arc_j[0]
                    && arc_i[1] == arc_j[1]
                    && arc_i[2] == arc_j[2]) {
                    /// two arcs of the same circle: no intersection
                    intersect_ij = 0;
                } else {
                    intersect_ij = intersectArcArc(arc_i, arc_j);
                    if (intersect_ij) {
                        printf("A=(%3.2f, %3.2f)\nB=(%3.2f, %3.2f)\n", i0[0], i0[1], i1[0], i1[1]);
                        printf("C=(%3.2f, %3.2f)\nD=(%3.2f, %3.2f)\n", j0[0], j0[1], j1[0], j1[1]);
                        printf("[INFO] [%s: %s] (%d %d) (%12.8lf %12.8lf %12.8lf) -- (%12.8lf %12.8lf %12.8lf)\n",
                                fnName,
                                "ArcArc",
                                i, j,
                                arc_i[0], arc_i[1], arc_i[2],
                                arc_j[0], arc_j[1], arc_j[2]
                              );
                    }
                }
            } else if (!i_is_arc &&  j_is_arc) {
                intersect_ij = intersectLineArc( i0, i1, arc_j);
                    if (intersect_ij) {
                        printf("[INFO] [%s: %s] (%d %d)\n",
                                fnName,
                                "LineArc",
                                i, j
                              );
                        printf("%12.8lf %12.8lf -- %12.8lf %12.8lf\n", i0[0], i0[1], i1[0], i1[1]);
                        printf("%12.8lf %12.8lf %12.8lf %12.8lf\n", arc_j[0], arc_j[1], arc_j[2], arc_j[3]);
                    }
            } else if ( i_is_arc && !j_is_arc) {
                intersect_ij = intersectLineArc( j0, j1, arc_i);
                    if (intersect_ij) {
                        printf("[INFO] [%s: %s] (%d %d)\n",
                                fnName,
                                "ArcLine",
                                i, j
                              );
                        printf("%12.8lf %12.8lf %12.8lf %12.8lf\n", arc_i[0], arc_i[1], arc_i[2], arc_i[3]);
                        printf("%12.8lf %12.8lf -- %12.8lf %12.8lf\n", j0[0], j0[1], j1[0], j1[1]);
                    }
            } else if (!i_is_arc && !j_is_arc) {
                intersect_ij = intersectLineSegments(i0, i1, j0, j1, NULL);
                    if (intersect_ij) {
                        printf("[INFO] [%s: %s] (%d %d)\n",
                                fnName,
                                "LineSegments",
                                i, j
                              );
                    }
            }

            intersect = intersect || intersect_ij;
        }
    }

    return intersect;
}

//------------------------------------------------------------------------------

void DEPRECATED_compareConfigAngleAndChildAngle(
        treeNode* node,
        puzzlerOptions* puzzler
) {

    if (isExterior(node)) {
        printf("[COMPARE] node: %c (exterior)\n", getNodeName(node));
        PS_printTree(node, puzzler);
    } else {
        printf("[COMPARE] node: %c\n", getNodeName(node));
        printf("[COMPARE]  i cfg           child         diff\n");
        double cfg = 0.0;
        for (int i = 0; i < node->childCount; i++) {
            cfg += getArcAngle(node->cfg, i);
            double cfgDegree = toDegree(cfg);
            double childAngle = getChildAngleByIndex(node, i);
            double diff = cfgDegree - childAngle;
            printf("[COMPARE] %2d %12.8lf° %12.8lf° %+.8lf°\n", i, cfgDegree, childAngle, diff);
        }
        cfg += getArcAngle(node->cfg, node->childCount);
        printf("[COMPARE] %2d %12.8lf° ------------- -------------\n", node->childCount, toDegree(cfg));
    }

    for (int i = 0; i < node->childCount; i++) {
        DEPRECATED_compareConfigAngleAndChildAngle(getChild(node, i), puzzler);
    }
}

//------------------------------------------------------------------------------

/**
 Calculate the coordinates for the drawing with the given angle angles
 */
void determineNucleotideCoordinates(
    treeNode* const node,
    short const * const pair_table,
    unsigned short const length,
    const double unpairedDistance,
    const double pairedDistance,
    double * const x,
    double * const y
) {
    if (length < 1) {
        return;
    }

    /// Handle stem of current node
    /// TODO: bulges!
    if (node->stem_start >= 1) {
        stemBox *sBox = node->sBox;

        /// - prepare bulge information
        int leftBulges = 0;
        int rightBulges = 0;
        int currentBulge = 0;
        for (int bulge = 0; bulge < sBox->bulgeCount; ++bulge) {
            if (sBox->bulges[bulge][0] < 0.0) {
                ++rightBulges;
            } else {
                ++leftBulges;
            }
        }

        /// - left side
        int ntStart = node->stem_start;
        int ntEnd   = node->loop_start;
        int ntSegments = ntEnd - ntStart - leftBulges;
        double pStart[2] = {sBox->c[0] - sBox->e[0] * sBox->a[0] + sBox->e[1] * sBox->b[0],
                            sBox->c[1] - sBox->e[0] * sBox->a[1] + sBox->e[1] * sBox->b[1],
                           };
        double pEnd[2] = {sBox->c[0] + sBox->e[0] * sBox->a[0] + sBox->e[1] * sBox->b[0],
                          sBox->c[1] + sBox->e[0] * sBox->a[1] + sBox->e[1] * sBox->b[1],
                         };

        for (int nt = ntStart; nt < ntEnd; ++nt) {
            if (pair_table[nt] == 0) {
                // bulge
                getBulgeXY(sBox, currentBulge, &(x[nt-1]), &(y[nt-1]));
                ++currentBulge;
            } else {
                x[nt-1] = pStart[0] + (nt - ntStart - currentBulge) * (pEnd[0] - pStart[0]) / ntSegments;
                y[nt-1] = pStart[1] + (nt - ntStart - currentBulge) * (pEnd[1] - pStart[1]) / ntSegments;
            }
        }
        x[ntEnd-1] = pEnd[0];
        y[ntEnd-1] = pEnd[1];

        /// - right side
        ntStart = pair_table[node->loop_start];
        ntEnd   = pair_table[node->stem_start];
        ntSegments = ntEnd - ntStart - rightBulges;
        pStart[0] = sBox->c[0] + sBox->e[0] * sBox->a[0] - sBox->e[1] * sBox->b[0];
        pStart[1] = sBox->c[1] + sBox->e[0] * sBox->a[1] - sBox->e[1] * sBox->b[1];
        pEnd[0] = sBox->c[0] - sBox->e[0] * sBox->a[0] - sBox->e[1] * sBox->b[0];
        pEnd[1] = sBox->c[1] - sBox->e[0] * sBox->a[1] - sBox->e[1] * sBox->b[1];

        for (int nt = ntStart; nt < ntEnd; ++nt) {
            if (pair_table[nt] == 0) {
                // bulge
                getBulgeXY(sBox, currentBulge, &(x[nt-1]), &(y[nt-1]));
                ++currentBulge;
            } else {
                x[nt-1] = pStart[0] + (nt - ntStart - currentBulge + leftBulges) * (pEnd[0] - pStart[0]) / ntSegments;
                y[nt-1] = pStart[1] + (nt - ntStart - currentBulge + leftBulges) * (pEnd[1] - pStart[1]) / ntSegments;
            }
        }
        x[ntEnd-1] = pEnd[0];
        y[ntEnd-1] = pEnd[1];
    }

    /// loop
    config *cfg = node->cfg;
    if (cfg != NULL) {
        double center[2] = {node->lBox->c[0], node->lBox->c[1]};
        double radius = cfg->radius;
        double pairedAngle = distanceToAngle(radius, pairedDistance);

        /// - determine angle from loop to parent stem
        double startAngle = 0.0;
        stemBox *sBox = node->sBox;
        /*
        if (abs(center[0] - sBox->c[0]) < 1e-10) {
            if (center[1] > sBox->c[1]) {
                startAngle = MATH_PI + MATH_PI_HALF;
            } else {
                startAngle = MATH_PI_HALF;
            }
        } else if (abs(center[1] - sBox->c[1]) < 1e-10) {
            if (center[0] > sBox->c[0]) {
                startAngle = MATH_PI;
            } else {
                startAngle = 0.0;
            }
        } else {
        */
            startAngle = atan2((sBox->c[1] - center[1]), (sBox->c[0] - center[0]));
        //}
        startAngle -= pairedAngle / 2.0;

        /// - for all loop arcs
        configArc *cfgArc = NULL;
        int nt = node->loop_start;
        double angle = startAngle;
        double arcAngle = 0.0;
        int numberOfArcSegments = 0;
        for (int arc = 0; arc < cfg->numberOfArcs; ++arc) {
            cfgArc = &(cfg->cfgArcs[arc]);
            numberOfArcSegments = cfgArc->numberOfArcSegments;
            arcAngle = cfgArc->arcAngle;
            for (int arcSegment = 1; arcSegment < numberOfArcSegments; ++arcSegment) {
                angle = startAngle - arcSegment * ((arcAngle - pairedAngle) / numberOfArcSegments);
                x[nt] = center[0] + radius * cos(angle);
                y[nt] = center[1] + radius * sin(angle);
                ++nt;
            }
            nt = pair_table[nt+1];
            startAngle -= arcAngle;
        }
    }

    /// children
    for (int child = 0; child < node->childCount; ++child) {
        determineNucleotideCoordinates(node->children[child],
                                       pair_table, length,
                                       unpairedDistance, pairedDistance,
                                       x, y);
    }

    /// exterior
    x[0] = EXTERIOR_Y;
    y[0] = EXTERIOR_Y;
    int start = 1;
    if (pair_table[1] != 0) {
        start = pair_table[1] + 1;
    } else {
        start = 2;
    }
    for (int nt = start; nt <= length; ++nt) {
        if (pair_table[nt] == 0) {
            x[nt-1] = x[nt-2] + unpairedDistance;
            y[nt-1] = EXTERIOR_Y;
        } else {
            nt = pair_table[nt];
        }
    }

    return;
}

//------------------------------------------------------------------------------

int layout_RNApuzzler(
        short const * const pair_table,
        float *x,
        float *y,
        double *arc_coords,
        puzzlerOptions *puzzler
) {
    const char *fnName = "layout_RNApuzzler";

//    printf("[Puzzler Settings]\n");
//    printf("exterior: %d\n", puzzler->checkExteriorIntersections);
//    printf("ancestor: %d\n", puzzler->checkAncestorIntersections);
//    printf("siblings: %d\n", puzzler->checkSiblingIntersections);
//    printf("optimize: %d\n", puzzler->optimize);
//    printf("\n");

    if (puzzler->paired / puzzler->unpaired > 2.0) {
        printWarning(fnName, "paired:unpaired > 2.0 -> layout might be destroyed!\n");
    }

    int length = pair_table[0];
    //printf("RNA length: %d\n", length);

    /// turtle base information
    tBaseInformation* baseInformation = vrna_alloc((length + 1) * sizeof(tBaseInformation));
    for (int i = 0; i <= length; i++) {
        baseInformation[i].baseType = TYPE_BASE_NONE;
        baseInformation[i].distance = puzzler->unpaired;
        baseInformation[i].angle = 0.0;
        baseInformation[i].config = NULL;
    }

    /// generate default configuration for each loop
    cfgGenerateConfig(pair_table, baseInformation, puzzler->unpaired, puzzler->paired);

    if (printInitialConfig) {
        printf("** print initial config **\n");
        for (int i = 0; i <= length; i++) {
            if (baseInformation[i].config != NULL) {
                cfgPrintConfig(baseInformation[i].config);
            }
        }
    }

    /// RNAturtle
    computeAffineCoordinates(pair_table, puzzler->paired, puzzler->unpaired, baseInformation);

    /// Transform affine coordinates into cartesian coordinates
    double* myX = (double*) vrna_alloc(length * sizeof(double));
    double* myY = (double*) vrna_alloc(length * sizeof(double));
    affineToCartesianCoordinates(baseInformation, length, myX, myY);

    /// Build RNApuzzler configuration tree from cartesian coordinates
    double distBulge = sqrt(puzzler->unpaired * puzzler->unpaired - 0.25 * puzzler->unpaired * puzzler->unpaired);
    treeNode* tree = buildConfigtree(pair_table, baseInformation, myX, myY, distBulge);

    /// current and maximal number of changes applied to config
    puzzler->numberOfChangesAppliedToConfig = 0;

    /// DZ: should be dependent on the RNA length * 10 ???
    puzzler->maximumNumberOfConfigChangesAllowed = 25000;

    /// reset angle coordinates
    /*
    for (int i = 0; i < length+1; i++) {
        baseInformation[i].distance = puzzler->unpaired;
        baseInformation[i].angle = 0.0;
    }
    */

    /// RNApuzzler
    if (puzzler->checkExteriorIntersections || puzzler->checkSiblingIntersections || puzzler->checkAncestorIntersections) {
        /// - One execution of checkAndFixIntersections should always be sufficient
        updateBoundingBoxes(tree, puzzler);

        if (FANCY_PS) {
            PS_printFancyTree(tree, puzzler);
        }
        checkAndFixIntersections(tree, 0, puzzler);
        printf("\n");
        printInformation("CHANGE COUNT", "%d %s\n\n", puzzler->numberOfChangesAppliedToConfig, puzzler->filename);
    }

    /// use configuration created by RNApuzzler for RNAturtle
    //computeAffineCoordinates(pair_table, puzzler->paired, puzzler->unpaired, baseInformation);
    //affineToCartesianCoordinates(baseInformation, length, myX, myY);

    /// determine x and y coordinates from RNApuzzler result
    /*
    for (int i = 0; i < length; i++) {
      myX[i] = 50 + i;
      myY[i] = 100 - i;
    }
    */
    determineNucleotideCoordinates(tree,
                                   pair_table, length,
                                   puzzler->unpaired, puzzler->paired,
                                   myX, myY);

    /// this section is for finding and resolving intersections
    /// of branches of the exterior loop against each other
    /// stretch backbones of the exterior until the overlap is gone
    /// may result in wide pictures
    short checkIntersectionsOfExteriorBranches = 1;
    if (checkIntersectionsOfExteriorBranches) {
        // resolveExteriorChildrenIntersectionAffin(tree, pair_table, baseInformation, puzzler->unpaired, puzzler->allowFlipping);
        // resolveExteriorChildIntersections(tree, pair_table, baseInformation, puzzler->unpaired, puzzler->allowFlipping);
        // affineToCartesianCoordinates(baseInformation, length, myX, myY);

        resolveExteriorChildrenIntersectionXY(tree, pair_table, puzzler->unpaired, puzzler->allowFlipping, myX, myY);
        if (FANCY_PS) {

            if (tree->lBox) {
                free(tree->lBox);
                tree->lBox = NULL;
            }
            if (tree->sBox) {
                free(tree->sBox);
                tree->sBox = NULL;
            }
            PS_printFancyTree(tree, puzzler);
        }
    }

    /// for all loops: compute postscript arcs instead of lines
    if (puzzler->drawArcs) {
        computeAnglesAndCentersForPS(pair_table, myX, myY, baseInformation, arc_coords);
    }

    /*
    for (int i = 0; i < length+1; i++) {
        printf("baseInformation[%d]: %d\n", i, baseInformation[i].baseType);
    }
    */

    /// final check based on line segments and arc segments
    short printDetails = 0;
    short intersect = checkRemainingIntersections(myX, myY, arc_coords, printDetails, baseInformation, length);
    printInformation("RESULT FINAL", "%s %s\n\n", (intersect ? "FAIL   " : "SUCCESS"), puzzler->filename);

    freeTree(tree);
//    printf("tree\n");
    free(baseInformation);
//    printf("baseInformation\n");

    for (int i = 0; i < length; i++) {
        x[i] = myX[i];
        y[i] = myY[i];
    }

    free(myX);
    free(myY);

    return length;
}

