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

#include "coordinates.inc"

#include "ViennaRNA/plotting/RNApuzzler/postscript/postscriptArcs.h"
#include "ViennaRNA/plotting/RNApuzzler/vector_math.h"
#include "ViennaRNA/plotting/RNApuzzler/data/drawingconfig.h"

#include "ViennaRNA/plotting/RNApuzzler/RNAturtle.h"

//--------------------------------------------------------------------------------------------------------------------------


PUBLIC int
layout_RNAturtle(short const *const pair_table,
                 float              *x,
                 float              *y,
                 double             *arc_coords)
{
  const char        *fnName   = "layout_RNAturtle";
  const short       drawArcs  = 1;
  const short       paired    = 35.0;
  const short       unpaired  = 25.0;
  const int         length    = pair_table[0];

  // turtle base information
  tBaseInformation  *baseInformation = vrna_alloc((length + 1) * sizeof(tBaseInformation));

  for (int i = 0; i <= length; i++) {
    baseInformation[i].baseType = TYPE_BASE_NONE;
    baseInformation[i].distance = unpaired;
    baseInformation[i].angle    = 0.0;
    baseInformation[i].config   = NULL;
  }

  // generate default configuration for each loop
  cfgGenerateConfig(pair_table, baseInformation, unpaired, paired);

  // compute loop angles
  computeAffineCoordinates(pair_table, paired, unpaired, baseInformation);

  // transform affine coordinates into cartesian coordinates
  double  *myX  = (double *)vrna_alloc(length * sizeof(double));
  double  *myY  = (double *)vrna_alloc(length * sizeof(double));
  affineToCartesianCoordinates(baseInformation, length, myX, myY);

  if (drawArcs)
    // compute postscript arcs instead of lines
    computeAnglesAndCentersForPS(pair_table, myX, myY, baseInformation, arc_coords);

  for (int i = 0; i < length; i++) {
    x[i]  = myX[i];
    y[i]  = myY[i];
  }

  free(myX);
  free(myY);

  // free turtle struct
  free(baseInformation);

  return length;
}
