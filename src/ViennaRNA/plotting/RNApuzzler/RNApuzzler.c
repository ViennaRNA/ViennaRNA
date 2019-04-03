/*
 *      RNApuzzler main algorithm
 *
 *      c  Daniel Wiegreffe, Daniel Alexander, Dirk Zeckzer
 *      ViennaRNA package
 */

#include <stdlib.h>
#include <math.h>

#include "ViennaRNA/utils/basic.h"

#include "includes/coordinates.inc"
#include "includes/postscriptArcs.inc"
#include "includes/vector_math.inc"
#include "includes/definitions.inc"
#include "includes/drawingconfig.inc"
#include "includes/configtree.inc"
#include "includes/intersectLevelLines.inc"
#include "includes/resolveIntersections.inc"
#include "includes/resolveExteriorChildIntersections.inc"
#include "headers/boundingBoxes_struct.h"
#include "includes/boundingBoxes.inc"

#include "ViennaRNA/plotting/RNApuzzler/RNApuzzler.h"

#define DEBUG   0

/* ---------------------------------------------------------------------------- */
PUBLIC vrna_plot_options_puzzler_t *
vrna_plot_options_puzzler()
{
  vrna_plot_options_puzzler_t *puzzler =
    (vrna_plot_options_puzzler_t *)vrna_alloc(sizeof(vrna_plot_options_puzzler_t));

  /* drawing behavior */
  puzzler->drawArcs = 1;
  puzzler->paired   = 35.0;
  puzzler->unpaired = 25.0;

  /* intersection resolution behavior */
  puzzler->checkExteriorIntersections = 1;
  puzzler->checkSiblingIntersections  = 1;
  puzzler->checkAncestorIntersections = 1;
  puzzler->optimize                   = 1;

  /* import behavior - unused for now */
  puzzler->config = NULL;

  /* other stuff */
  puzzler->filename = NULL;

  puzzler->numberOfChangesAppliedToConfig = 0;
  puzzler->psNumber                       = 0;

  return puzzler;
}


PUBLIC void
vrna_plot_options_puzzler_free(vrna_plot_options_puzzler_t *puzzler)
{
  free(puzzler);
}


PRIVATE short
checkRemainingIntersections(double                  *x,
                            double                  *y,
                            double                  *arcCoords,
                            const short             printDetails,
                            const tBaseInformation  *baseInformation,
                            const int               length)
{
  char        *fnName       = "checkRemainingIntersections";
  const short skipExterior  = 0;
  short       intersect     = 0;
  double      arc_i[6];
  short       i_is_arc = 0;
  double      arc_j[6];
  short       j_is_arc = 0;

  for (int i = 3; i < length; i++) {
    arc_i[0]  = arcCoords[6 * i + 0];
    arc_i[1]  = arcCoords[6 * i + 1];
    arc_i[2]  = arcCoords[6 * i + 2];
    arc_i[3]  = arcCoords[6 * i + 3];
    arc_i[4]  = arcCoords[6 * i + 4];
    arc_i[5]  = arcCoords[6 * i + 5];

    i_is_arc = (arc_i[0] != -1);
    double  i0[2] = {
      x[i - 1], y[i - 1]
    };
    double  i1[2] = {
      x[i - 0], y[i - 0]
    };

    if (skipExterior
        && ((i0[1] <= EXTERIOR_Y) || (i1[1] <= EXTERIOR_Y)))
      continue;

    short intersectIExterior = 0;
    if (baseInformation[i + 0].baseType == TYPE_EXTERIOR
        && baseInformation[i + 1].baseType == TYPE_EXTERIOR) {
      if (i_is_arc) {
        /* exterior line */
        double  xmin  = fmin(i0[0], i1[0]);
        double  xmax  = fmax(i0[0], i1[0]);
        double  p1[2] = {
          xmin, EXTERIOR_Y
        };
        double  p2[2] = {
          xmax, EXTERIOR_Y
        };

        intersectIExterior = intersectLineArc(p1, p2, arc_i);
      } else {
        intersectIExterior = ((i0[1] <= EXTERIOR_Y) != (i1[1] <= EXTERIOR_Y));
      }
    }

    intersect = intersect || intersectIExterior;

    for (int j = 1; j < i - 1; j++) {
      arc_j[0]  = arcCoords[6 * j + 0];
      arc_j[1]  = arcCoords[6 * j + 1];
      arc_j[2]  = arcCoords[6 * j + 2];
      arc_j[3]  = arcCoords[6 * j + 3];
      arc_j[4]  = arcCoords[6 * j + 4];
      arc_j[5]  = arcCoords[6 * j + 5];
      j_is_arc  = (arc_j[0] != -1);
      double  j0[2] = {
        x[j - 1], y[j - 1]
      };
      double  j1[2] = {
        x[j - 0], y[j - 0]
      };

      if (skipExterior && ((j0[1] <= EXTERIOR_Y) || (j1[1] <= EXTERIOR_Y)))

        continue;

      short intersect_ij = 0;

      if (i_is_arc && j_is_arc) {
        if (arc_i[0] == arc_j[0]
            && arc_i[1] == arc_j[1]
            && arc_i[2] == arc_j[2])
          /* two arcs of the same circle: no intersection */
          intersect_ij = 0;
        else
          intersect_ij = intersectArcArc(arc_i, arc_j);
      } else if (!i_is_arc && j_is_arc) {
        intersect_ij = intersectLineArc(i0, i1, arc_j);
      } else if (i_is_arc && !j_is_arc) {
        intersect_ij = intersectLineArc(j0, j1, arc_i);
      } else if (!i_is_arc && !j_is_arc) {
        intersect_ij = intersectLineSegments(i0, i1, j0, j1, NULL);
      }

      intersect = intersect || intersect_ij;
    }
  }

  return intersect;
}


/* ------------------------------------------------------------------------------ */

/**
 * Calculate the coordinates for the drawing with the given angle angles
 */
PRIVATE void
determineNucleotideCoordinates(treeNode *const      node,
                               short const *const   pair_table,
                               unsigned short const length,
                               const double         unpairedDistance,
                               const double         pairedDistance,
                               double *const        x,
                               double *const        y)
{
  if (length < 1)
    return;

  /* Handle stem of current node */
  if (node->stem_start >= 1) {
    stemBox *sBox = node->sBox;

    /* prepare bulge information */
    int     leftBulges    = 0;
    int     rightBulges   = 0;
    int     currentBulge  = 0;

    for (int bulge = 0; bulge < sBox->bulgeCount; ++bulge) {
      if (sBox->bulges[bulge][0] < 0.0)
        ++rightBulges;
      else
        ++leftBulges;
    }

    /* left side */
    int     ntStart     = node->stem_start;
    int     ntEnd       = node->loop_start;
    int     ntSegments  = ntEnd - ntStart - leftBulges;
    double  pStart[2]   = {
      sBox->c[0] - sBox->e[0] * sBox->a[0] + sBox->e[1] * sBox->b[0],
      sBox->c[1] - sBox->e[0] * sBox->a[1] + sBox->e[1] * sBox->b[1],
    };
    double  pEnd[2] = {
      sBox->c[0] + sBox->e[0] * sBox->a[0] + sBox->e[1] * sBox->b[0],
      sBox->c[1] + sBox->e[0] * sBox->a[1] + sBox->e[1] * sBox->b[1],
    };

    for (int nt = ntStart; nt < ntEnd; ++nt) {
      if (pair_table[nt] == 0) {
        /* bulge */
        getBulgeXY(sBox, currentBulge, &(x[nt - 1]), &(y[nt - 1]));
        ++currentBulge;
      } else {
        x[nt - 1] = pStart[0] + (nt - ntStart - currentBulge) * (pEnd[0] - pStart[0]) / ntSegments;
        y[nt - 1] = pStart[1] + (nt - ntStart - currentBulge) * (pEnd[1] - pStart[1]) / ntSegments;
      }
    }
    x[ntEnd - 1]  = pEnd[0];
    y[ntEnd - 1]  = pEnd[1];

    /* right side */
    ntStart     = pair_table[node->loop_start];
    ntEnd       = pair_table[node->stem_start];
    ntSegments  = ntEnd - ntStart - rightBulges;
    pStart[0]   = sBox->c[0] + sBox->e[0] * sBox->a[0] - sBox->e[1] * sBox->b[0];
    pStart[1]   = sBox->c[1] + sBox->e[0] * sBox->a[1] - sBox->e[1] * sBox->b[1];
    pEnd[0]     = sBox->c[0] - sBox->e[0] * sBox->a[0] - sBox->e[1] * sBox->b[0];
    pEnd[1]     = sBox->c[1] - sBox->e[0] * sBox->a[1] - sBox->e[1] * sBox->b[1];

    for (int nt = ntStart; nt < ntEnd; ++nt) {
      if (pair_table[nt] == 0) {
        /* bulge */
        getBulgeXY(sBox, currentBulge, &(x[nt - 1]), &(y[nt - 1]));
        ++currentBulge;
      } else {
        x[nt - 1] = pStart[0] + (nt - ntStart - currentBulge + leftBulges) * (pEnd[0] - pStart[0]) /
                    ntSegments;
        y[nt - 1] = pStart[1] + (nt - ntStart - currentBulge + leftBulges) * (pEnd[1] - pStart[1]) /
                    ntSegments;
      }
    }
    x[ntEnd - 1]  = pEnd[0];
    y[ntEnd - 1]  = pEnd[1];
  }

  /* loop */
  config *cfg = node->cfg;
  if (cfg != NULL) {
    double    center[2] = {
      node->lBox->c[0], node->lBox->c[1]
    };
    double    radius      = cfg->radius;
    double    pairedAngle = distanceToAngle(radius, pairedDistance);

    /* determine angle from loop to parent stem */
    double    startAngle  = 0.0;
    stemBox   *sBox       = node->sBox;
    startAngle  = atan2((sBox->c[1] - center[1]), (sBox->c[0] - center[0]));
    startAngle  -= pairedAngle / 2.0;

    /* for all loop arcs */
    configArc *cfgArc             = NULL;
    int       nt                  = node->loop_start;
    double    angle               = startAngle;
    double    arcAngle            = 0.0;
    int       numberOfArcSegments = 0;
    for (int arc = 0; arc < cfg->numberOfArcs; ++arc) {
      cfgArc              = &(cfg->cfgArcs[arc]);
      numberOfArcSegments = cfgArc->numberOfArcSegments;
      arcAngle            = cfgArc->arcAngle;

      for (int arcSegment = 1; arcSegment < numberOfArcSegments; ++arcSegment) {
        angle = startAngle - arcSegment * ((arcAngle - pairedAngle) / numberOfArcSegments);
        x[nt] = center[0] + radius * cos(angle);
        y[nt] = center[1] + radius * sin(angle);
        ++nt;
      }
      nt          = pair_table[nt + 1];
      startAngle  -= arcAngle;
    }
  }

  /* children */
  for (int child = 0; child < node->childCount; ++child) {
    determineNucleotideCoordinates(node->children[child],
                                   pair_table, length,
                                   unpairedDistance, pairedDistance,
                                   x, y);
  }

  /* exterior */
  x[0]  = EXTERIOR_Y;
  y[0]  = EXTERIOR_Y;
  int start = 1;
  if (pair_table[1] != 0)
    start = pair_table[1] + 1;
  else
    start = 2;

  for (int nt = start; nt <= length; ++nt) {
    if (pair_table[nt] == 0) {
      x[nt - 1] = x[nt - 2] + unpairedDistance;
      y[nt - 1] = EXTERIOR_Y;
    } else {
      nt = pair_table[nt];
    }
  }

  return;
}


#if DEBUG
/* ------------------------------------------------------------------------------ */

/* debug method for intersection regression test */
PRIVATE int
printInformation(const char *fnName,
                 const char *format,
                 ...)
{
  FILE    *stream = stdout;

  va_list args;

  if (fnName != NULL)
    fprintf(stream, "[INFORMATION] [%s] ", fnName);

  va_start(args, format);
  int ret = vfprintf(stream, format, args);
  va_end(args);
  fflush(stream);
  return ret;
}


#endif

PUBLIC int
vrna_plot_coords_puzzler(const char                   *structure,
                         float                        **x,
                         float                        **y,
                         double                       **arc_coords,
                         vrna_plot_options_puzzler_t  *options)
{
  if (structure) {
    int   ret = 0;
    short *pt = vrna_ptable(structure);

    ret = vrna_plot_coords_puzzler_pt(pt, x, y, arc_coords, options);

    free(pt);

    return ret;
  }

  if (x)
    *x = NULL;

  if (y)
    *y = NULL;

  if (arc_coords)
    *arc_coords = NULL;

  return 0;
}


PUBLIC int
vrna_plot_coords_puzzler_pt(short const *const          pair_table,
                            float                       **x,
                            float                       **y,
                            double                      **arc_coords,
                            vrna_plot_options_puzzler_t *options)
{
  int                         length = pair_table[0];
  vrna_plot_options_puzzler_t *puzzler;

  if ((pair_table) && (x) && (y)) {
    *x  = (float *)vrna_alloc(sizeof(float) * (length + 1));
    *y  = (float *)vrna_alloc(sizeof(float) * (length + 1));

    /* apply default options if none provided */
    if (options) {
      puzzler = options;
    } else {
      puzzler           = vrna_plot_options_puzzler();
      puzzler->filename = NULL;
      puzzler->drawArcs = (arc_coords) ? 1 : 0;

      puzzler->checkAncestorIntersections = 1;
      puzzler->checkSiblingIntersections  = 1;
      puzzler->checkExteriorIntersections = 1;
      puzzler->allowFlipping              = 0;
      puzzler->optimize                   = 1;
    }

    /* turtle base information */
    tBaseInformation *baseInformation = vrna_alloc((length + 1) * sizeof(tBaseInformation));

    for (int i = 0; i <= length; i++) {
      baseInformation[i].baseType = TYPE_BASE_NONE;
      baseInformation[i].distance = puzzler->unpaired;
      baseInformation[i].angle    = 0.0;
      baseInformation[i].config   = NULL;
    }

    /* generate default configuration for each loop */
    cfgGenerateConfig(pair_table, baseInformation, puzzler->unpaired, puzzler->paired);

    /* RNAturtle */
    computeAffineCoordinates(pair_table, puzzler->paired, puzzler->unpaired, baseInformation);

    /* Transform affine coordinates into cartesian coordinates */
    double  *myX  = (double *)vrna_alloc(length * sizeof(double));
    double  *myY  = (double *)vrna_alloc(length * sizeof(double));
    affineToCartesianCoordinates(baseInformation, length, myX, myY);

    /* Build RNApuzzler configuration tree from cartesian coordinates */
    double  distBulge = sqrt(puzzler->unpaired *
                             puzzler->unpaired -
                             0.25 *
                             puzzler->unpaired *
                             puzzler->unpaired);

    treeNode *tree = buildConfigtree(pair_table, baseInformation, myX, myY, distBulge);

    /* current and maximal number of changes applied to config */
    puzzler->numberOfChangesAppliedToConfig = 0;

    puzzler->maximumNumberOfConfigChangesAllowed = 25000;

#if 0
    /* reset angle coordinates */
    for (int i = 0; i < length + 1; i++) {
      baseInformation[i].distance = puzzler->unpaired;
      baseInformation[i].angle    = 0.0;
    }
#endif

    /* RNApuzzler */
    if (puzzler->checkExteriorIntersections ||
        puzzler->checkSiblingIntersections ||
        puzzler->checkAncestorIntersections) {
      /* - One execution of checkAndFixIntersections should always be sufficient */
      updateBoundingBoxes(tree, puzzler);
      checkAndFixIntersections(tree, 0, puzzler);
    }

    /*
     * use configuration created by RNApuzzler for RNAturtle
     * computeAffineCoordinates(pair_table, puzzler->paired, puzzler->unpaired, baseInformation);
     * affineToCartesianCoordinates(baseInformation, length, myX, myY);
     */

    determineNucleotideCoordinates(tree,
                                   pair_table, length,
                                   puzzler->unpaired, puzzler->paired,
                                   myX, myY);

    /*
     * this section is for finding and resolving intersections
     * of branches of the exterior loop against each other
     * stretch backbones of the exterior until the overlap is gone
     * may result in wide pictures
     */

    short checkIntersectionsOfExteriorBranches = 1;

    if (checkIntersectionsOfExteriorBranches) {
      resolveExteriorChildrenIntersectionXY(tree,
                                            pair_table,
                                            puzzler->unpaired,
                                            puzzler->allowFlipping,
                                            myX,
                                            myY);
    }

    /* for all loops: compute postscript arcs instead of lines */
    if ((puzzler->drawArcs) && (arc_coords)) {
      *arc_coords = (double *)vrna_alloc(sizeof(double) * 6 * length);

      for (int i = 0; i < length; i++) {
        (*arc_coords)[6 * i + 0]  = -1;
        (*arc_coords)[6 * i + 1]  = -1.;
        (*arc_coords)[6 * i + 2]  = -1.;
        (*arc_coords)[6 * i + 3]  = -1.;
        (*arc_coords)[6 * i + 4]  = -1.;
        (*arc_coords)[6 * i + 5]  = -1.;
      }

      computeAnglesAndCentersForPS(pair_table, myX, myY, baseInformation, *arc_coords);

      /* final check based on line segments and arc segments */
      short printDetails  = 0;
      short intersect     = checkRemainingIntersections(myX,
                                                        myY,
                                                        *arc_coords,
                                                        printDetails,
                                                        baseInformation,
                                                        length);

#if DEBUG
      /* Debug call for intersection regression test, uncomment for output */
      printInformation("RESULT FINAL",
                       "%s %s\n\n",
                       (intersect ? "FAIL   " : "SUCCESS"),
                       puzzler->filename);
#endif
    } else if (arc_coords) {
      *arc_coords = NULL;
    }

    freeTree(tree);
    free(baseInformation);

    for (int i = 0; i < length; i++) {
      (*x)[i] = myX[i];
      (*y)[i] = myY[i];
    }

    free(myX);
    free(myY);

    if (!options)
      vrna_plot_options_puzzler_free(puzzler);

    return length;
  }

  if (x)
    *x = NULL;

  if (y)
    *y = NULL;

  if (arc_coords)
    *arc_coords = NULL;

  return 0;
}
