/*
 *      RNAturtle algorithm
 *
 *      c  Daniel Wiegreffe, Daniel Alexander, Dirk Zeckzer
 *      ViennaRNA package
 */


#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "ViennaRNA/utils/basic.h"

#include "includes/coordinates.inc"
#include "includes/postscriptArcs.inc"
#include "includes/vector_math.inc"
#include "includes/drawingconfig.inc"
#include "headers/tBaseInformation_struct.h"

#include "ViennaRNA/plotting/RNApuzzler/RNAturtle.h"

/* -------------------------------------------------------------------------------------------------------------------------- */

PUBLIC int
vrna_plot_coords_turtle(const char  *structure,
                        float       **x,
                        float       **y,
                        double      **arc_coords)
{
  if (structure) {
    int   ret = 0;
    short *pt = vrna_ptable(structure);

    ret = vrna_plot_coords_turtle_pt(pt, x, y, arc_coords);

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
vrna_plot_coords_turtle_pt(short const *const pair_table,
                           float              **x,
                           float              **y,
                           double             **arc_coords)
{
  const short drawArcs  = 1;
  const short paired    = 35.0;
  const short unpaired  = 25.0;
  const int   length    = pair_table[0];

  if ((pair_table) && (x) && (y)) {
    *x  = (float *)vrna_alloc(sizeof(float) * (length + 1));
    *y  = (float *)vrna_alloc(sizeof(float) * (length + 1));

    /* turtle base information */
    tBaseInformation *baseInformation = vrna_alloc((length + 1) * sizeof(tBaseInformation));

    for (int i = 0; i <= length; i++) {
      baseInformation[i].baseType = TYPE_BASE_NONE;
      baseInformation[i].distance = unpaired;
      baseInformation[i].angle    = 0.0;
      baseInformation[i].config   = NULL;
    }

    /* generate default configuration for each loop */
    cfgGenerateConfig(pair_table, baseInformation, unpaired, paired);

    /* compute loop angles */
    computeAffineCoordinates(pair_table, paired, unpaired, baseInformation);

    /* transform affine coordinates into cartesian coordinates */
    double  *myX  = (double *)vrna_alloc(length * sizeof(double));
    double  *myY  = (double *)vrna_alloc(length * sizeof(double));
    affineToCartesianCoordinates(baseInformation, length, myX, myY);

    if ((drawArcs) && (arc_coords)) {
      *arc_coords = (double *)vrna_alloc(sizeof(double) * 6 * length);

      for (int i = 0; i < length; i++) {
        (*arc_coords)[6 * i + 0]  = -1;
        (*arc_coords)[6 * i + 1]  = -1.;
        (*arc_coords)[6 * i + 2]  = -1.;
        (*arc_coords)[6 * i + 3]  = -1.;
        (*arc_coords)[6 * i + 4]  = -1.;
        (*arc_coords)[6 * i + 5]  = -1.;
      }

      /* compute postscript arcs instead of lines */
      computeAnglesAndCentersForPS(pair_table, myX, myY, baseInformation, *arc_coords);
    } else if (arc_coords) {
      *arc_coords = NULL;
    }

    for (int i = 0; i < length; i++) {
      (*x)[i] = myX[i];
      (*y)[i] = myY[i];
    }

    free(myX);
    free(myY);

    /* free turtle struct */
    free(baseInformation);

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
