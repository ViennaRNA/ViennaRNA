#ifndef RNAPUZZLER_COORDINATE_TRANSFORM
#define RNAPUZZLER_COORDINATE_TRANSFORM

/*
 *      RNAturtle algorithm
 *
 *      c  Daniel Wiegreffe, Daniel Alexander, Dirk Zeckzer
 *      ViennaRNA package
 */

#include "../headers/tBaseInformation_struct.h"
#include "vector_math.inc"
#include "drawingconfig.inc"

/**
 * Compute affin angles for all loops
 */
PRIVATE void
computeAffineCoordinates(short const *const       pair_table,
                         const double             paired,
                         const double             unpaired,
                         tBaseInformation *const  baseInformation);


/**
 * Calculate the coordinates for the drawing with the given affin angles
 */
PRIVATE void
affineToCartesianCoordinates(tBaseInformation *const  baseInformation,
                             unsigned short const     length,
                             double *const            x,
                             double *const            y);


/* -------------------------------------------------------------------------------------------------------------------------- */

/**
 * Handle bases of the exterior loop
 *
 * @returns first paired base
 */
PRIVATE short
handleExteriorBases(const short *const  pair_table,
                    short               currentBase,
                    tBaseInformation    *baseInformation,
                    const int           direction)
{
  short const length = pair_table[0];

  if (currentBase > 1) {
    baseInformation[currentBase].angle = baseInformation[currentBase].angle + direction *
                                         MATH_PI_HALF;
    baseInformation[currentBase].baseType = TYPE_EXTERIOR;
  }

  while (currentBase < length && pair_table[currentBase] <= 0) {
    baseInformation[currentBase + 1].angle  = 0.0;
    baseInformation[currentBase].baseType   = TYPE_EXTERIOR;
    currentBase++;
  }
  /*
   * For following stem, since it cannot alter the first two bases
   * which may contain loop exit angles.
   */
  if ((currentBase + 1) <= (length)) {
    baseInformation[currentBase + 1].angle  = direction * MATH_PI_HALF;
    baseInformation[currentBase].baseType   = TYPE_EXTERIOR;
  } else {
    baseInformation[currentBase].baseType = TYPE_EXTERIOR;
  }

  return currentBase;
}


PRIVATE void
countLoopPairs(short *const       m,
               short *const       n,
               short              i,
               const short *const pair_table)
{
  short const end = pair_table[i++];  /* end is 1st base of 2nd stem half */

  *m = *n = 1;                        /* Base pair of enclosing stem, consecutive pair from last stem base to first loop base */
  while (i < end) {
    if (pair_table[i] <= 0 || pair_table[i] < i) {
      (*n)++;
      i++;       /* Go to consecutive base */
    } else {
      /* base pair ahead */
      (*m)++;
      i = pair_table[i];       /* Go to base pair partner */
    }
  }

  return;
}


/* -------------------------------------------------------------------------------------------------------------------------- */
PRIVATE void
handleStem(const short *const       pair_table,
           short                    currentBase,
           const double             paired,
           const double             unpaired,
           tBaseInformation *const  baseInformation,
           int                      direction);


/* -------------------------------------------------------------------------------------------------------------------------- */

/*---------------------------------------------------------------------------*/

/**
 * Detect bulges (internal loops with only one unpaired base)
 */
PRIVATE int
detectBulge(short               start,
            const short *const  pair_table)
{
  int bulge   = 0;
  int end     = pair_table[start];
  int iterate = 1;
  int old     = 0;
  int i       = start + 1;
  int side    = 0;

  do {
    if (pair_table[i] > 0) {
      /* There was a basepair before */
      if (iterate > 0) {
        /* Stem entry */
        if (pair_table[i] == old) {
          side++;
          i++;
        }        /* Bulge found */
        else {
          /* It's only a bulge if the paired base is the start base or the bulge is on the other strand */
          if (start == pair_table[i] || pair_table[i] == end - 2) {
            bulge = pair_table[i];
            break;
          } else {
            break;
          }
        }
      }      /* Found first basepair */
      else {
        iterate++;
        old = i;
        i   = pair_table[i];
      }
    }    /* Consecutive Base found */
    else {
      /* There was a basepair before */
      if (iterate > 0) {
        iterate = 0;
        i++;
      }      /* Just the next consecutive base */
      else {
        i++;
      }
    }
  } while (i > start);

  return bulge;
}


PRIVATE void
handleLoop(short              i,
           const short *const pair_table,
           const double       paired,
           const double       unpaired,
           tBaseInformation   *baseInformation,
           int                direction)
{
  int         start = i;
  short const end = pair_table[i];  /* end is 1st base of 2nd stem half */
  short       m, n;                 /* m - no. of base pairs; n - no. of consecutive pairs */

  countLoopPairs(&m, &n, i, pair_table);


  int         bulge = detectBulge(i, pair_table);

  if (bulge > 0 && n - m == 1) {
    /* detected bulge with a single unpaired base */
    int     length = (unpaired * (n - m + 1)) / 2;

    double  alpha = acos(unpaired / (2 * length));

    if (pair_table[i + 1] == 0) {
      baseInformation[i + 1].angle            = baseInformation[i + 1].angle + direction * alpha;
      baseInformation[i].baseType             = TYPE_BULGE;
      baseInformation[pair_table[i]].baseType = TYPE_BULGE;
      i++;

      baseInformation[i + 1].angle  = -direction * alpha * 2;
      baseInformation[i].baseType   = TYPE_BULGE;
      i++;
      if ((i + 1) <= pair_table[0])
      {	
      	baseInformation[i + 1].angle            = direction * alpha;
      }
      baseInformation[i].baseType             = TYPE_BULGE;
      baseInformation[pair_table[i]].baseType = TYPE_BULGE;

      handleStem(pair_table, i, paired, unpaired, baseInformation, direction);

      i = pair_table[i];
    } else {
      baseInformation[i + 1].angle  = baseInformation[i + 1].angle + 0.0;
      baseInformation[i].baseType   = TYPE_BULGE;
      i++;

      baseInformation[i + 1].angle    = baseInformation[i + 1].angle + 0.0;
      baseInformation[i + 1].baseType = TYPE_BULGE;

      baseInformation[i + 2].angle    = baseInformation[i + 2].angle + 0.0;
      baseInformation[i + 1].baseType = TYPE_BULGE;

      handleStem(pair_table, i, paired, unpaired, baseInformation, direction);

      i = pair_table[i];

      baseInformation[i + 1].angle  = baseInformation[i + 1].angle + direction * alpha;
      baseInformation[i].baseType   = TYPE_BULGE;
      i++;

      baseInformation[i + 1].angle  = -direction * alpha * 2;
      baseInformation[i].baseType   = TYPE_BULGE;
      i++;
      
      if ((i + 1) <= pair_table[0])
      {	
      	baseInformation[i + 1].angle  = direction * alpha;
      }
      baseInformation[i].baseType   = TYPE_BULGE;
    }
  } else {
    /* loop is not bulge with a single unpaired base */

    /*
     * ******************************************************************
     * ***************** Loop Drawing Algorithm - Start *****************
     * ******************************************************************
     */


    /* variables for config drawing */
    double  angle_over_paired;                /* alpha      at calculation */
    double  current_angle;                    /* phi        at calculation */
    double  current_bb_angle;                 /* beta       at calculation */
    double  current_distance;                 /* b          at calculation */
    double  current_delta_ab;                 /* delta_ab   at calculation */
    double  current_delta_bb;                 /* delta_bb   at calculation */
    int     current_stem_count;

    config  *cfg        = baseInformation[start].config;
    int     currentArc  = 0;
    double  r           = cfg->radius;

    angle_over_paired = 2 * asin(paired / (2 * r));

    current_angle     = getArcAngle(cfg, currentArc);
    current_bb_angle  = (current_angle - angle_over_paired) /
                        (cfg->cfgArcs[currentArc]).numberOfArcSegments;
    current_distance  = sqrt(2 * r * r * (1 - cos(current_bb_angle)));
    current_delta_ab  = 0.5 * (MATH_PI + angle_over_paired + current_bb_angle);
    current_delta_bb  = MATH_PI + current_bb_angle;
    ++currentArc;

    /* consider the direction that was given before (e.g. because of a bulge) */
    baseInformation[i + 1].angle = baseInformation[i + 1].angle + direction *
                                   (MATH_PI - current_delta_ab);

    baseInformation[i].distance = current_distance;

    current_stem_count = 0;

    /* Double Loop Check: Loop is followed directly by a loop */
    if (baseInformation[i].baseType == TYPE_LOOP1)
      baseInformation[i].baseType = TYPE_LOOP2;
    else
      baseInformation[i].baseType = TYPE_LOOP1;

    i++;

    while (i < end) {
      if (pair_table[i] <= 0) {
        /* unpaired */

        /* handle backbone element */
        baseInformation[i + 1].angle  = -direction * (current_delta_bb - MATH_PI);
        baseInformation[i].distance   = current_distance;

        baseInformation[i].baseType = TYPE_LOOP1;
        i++;
      } else if (pair_table[i] > i) {
        /* base pair ahead */

        /* i is the beginning of a stem and therefore the end of the current arc */
        baseInformation[i + 1].angle = direction * (MATH_PI - current_delta_ab);
        current_stem_count++;
        baseInformation[i].baseType = TYPE_LOOP1;
        handleStem(pair_table, i, paired, unpaired, baseInformation, direction);  /* Recurse multiloop */
        i = pair_table[i];                                                        /* Go to base pair partner */
      } else {
        /* base pair behind */


        /*
         * i is the end of a stem that returned to the loop
         * recalculation if angles and distances is nessecary because the next arc has been reached
         */
        if (current_stem_count == 1) {
          current_stem_count  = 0;
          current_angle       = getArcAngle(cfg, currentArc);
          current_bb_angle    = (current_angle - angle_over_paired) /
                                (cfg->cfgArcs[currentArc]).numberOfArcSegments;
          current_distance  = sqrt(2 * r * r * (1 - cos(current_bb_angle)));
          current_delta_ab  = 0.5 * (MATH_PI + angle_over_paired + current_bb_angle);
          current_delta_bb  = MATH_PI + current_bb_angle;
          ++currentArc;
        }

        /* consider the direction that was given before (e.g. because of a bulge) */
        baseInformation[i + 1].angle = baseInformation[i + 1].angle + direction *
                                       (MATH_PI - current_delta_ab);
        baseInformation[i].distance = current_distance;
        baseInformation[i].baseType = TYPE_LOOP1;

        i++;
      }
    }

    if ((i + 1) <= pair_table[0])
      baseInformation[i + 1].angle = direction * (MATH_PI - current_delta_ab);

    baseInformation[i].baseType = TYPE_LOOP1;


    /*
     * ******************************************************************
     * ****************** Loop Drawing Algorithm - End ******************
     * ******************************************************************
     */
  }

  return;
}


/* -------------------------------------------------------------------------------------------------------------------------- */

PRIVATE void
handleStem(const short *const       pair_table,
           short                    i,
           const double             paired,
           const double             unpaired,
           tBaseInformation *const  baseInformation,
           int                      direction)
{
  short end = pair_table[i] + 1;   /* first base after stem */

  /*
   * Skip first two bases of stem, angle already written by
   * either loop (+exit angle) or handle_unpaired
   * i+=2;
   */
  baseInformation[i].baseType = TYPE_STEM;
  i++;

  while (pair_table[i] > 0 &&                     /* i is base paired */
         (pair_table[i] == end - 1 ||             /* First position of stem, continue! */
          pair_table[i] + 1 == pair_table[i - 1]  /* Detect bulges on opposite strand */
         )
         ) {
    baseInformation[i + 1].angle  = 0.0;
    baseInformation[i].baseType   = TYPE_STEM;
    i++;
  }
  if (pair_table[i] == end - 1) {
  } else {
    handleLoop(--i, pair_table, paired, unpaired, baseInformation, direction);     /* i - last base of first stem half */
  }

  i                           = pair_table[i]; /* set i as base pairing partner of last stem base */
  baseInformation[i].baseType = TYPE_STEM;
  i++;

  while (i < end && i < pair_table[0]) {
    baseInformation[i].baseType = TYPE_STEM;
    i++;
  }

  return;
}


/* ------------------------------------------------------------------------------ */

/**
 * Compute angle angles for all loops
 */
PRIVATE void
computeAffineCoordinates(short const *const       pair_table,
                         const double             paired,
                         const double             unpaired,
                         tBaseInformation *const  baseInformation)
{
  /* Initialization */
  short const length      = pair_table[0];
  short       currentBase = 1;
  const int   direction   = -1;

  baseInformation[0].angle = 0.0;

  /*
   * For first stem, since handle_stem cannot alter the first two bases
   * which may contain loop exit angles (in other stems).
   */
  if (2 <= length) {
    baseInformation[1].angle  = baseInformation[0].angle;
    baseInformation[2].angle  = baseInformation[1].angle;
  }

  /* reset dangle_count */
  int dangle_count = 0;
  while (currentBase < length) {
    if (pair_table[currentBase] <= 0) {
      if (currentBase > 1)

        baseInformation[currentBase - 1].baseType = TYPE_EXTERIOR;

      currentBase = handleExteriorBases(pair_table, currentBase, baseInformation, direction);
      /* returns first paired base currentBase */
      dangle_count++;
    }

    if (currentBase < length) {
      /* it was not the dangling end */
      if (pair_table[currentBase] - pair_table[currentBase - 1] != 1
          && pair_table[currentBase] != 0
          && pair_table[currentBase - 1] != 0) {
        if (currentBase == 1) {
          if (dangle_count < 1) {
            baseInformation[0].angle              =
              baseInformation[1].angle            =
                baseInformation[2].angle          = -MATH_PI_HALF;
            baseInformation[currentBase].baseType = TYPE_EXTERIOR;
          }

          handleStem(pair_table, currentBase, paired, unpaired, baseInformation, direction);
          currentBase = pair_table[currentBase] + 1;

          /* Check, if there is a minor dangling end at the end */
          if (currentBase == length) {
            baseInformation[currentBase - 1].baseType = TYPE_EXTERIOR;
            baseInformation[currentBase].baseType     = TYPE_EXTERIOR;
            baseInformation[currentBase].angle        = -MATH_PI_HALF;
          }

          continue;
        } else {
          baseInformation[currentBase].angle = baseInformation[currentBase].angle +
                                               direction * MATH_PI_HALF;
          baseInformation[currentBase + 1].distance = unpaired;
          baseInformation[currentBase - 1].baseType = TYPE_EXTERIOR;
          baseInformation[currentBase + 1].angle    = baseInformation[currentBase + 1].angle +
                                                      direction * MATH_PI_HALF;
          baseInformation[currentBase].baseType = TYPE_EXTERIOR;

          dangle_count++;
        }
      }

      handleStem(pair_table, currentBase, paired, unpaired, baseInformation, direction);

      currentBase = pair_table[currentBase] + 1;       /* currentBase is next base after stem */

      if (currentBase == length) {
        baseInformation[currentBase - 1].baseType = TYPE_EXTERIOR;
        currentBase                               = handleExteriorBases(pair_table,
                                                                        currentBase,
                                                                        baseInformation,
                                                                        direction);
      }
    }
  }
  baseInformation[length].baseType = TYPE_EXTERIOR;
}


/* ------------------------------------------------------------------------------ */

/**
 * Calculate the coordinates for the drawing with the given angle angles
 */
PRIVATE void
affineToCartesianCoordinates(tBaseInformation *const  baseInformation,
                             unsigned short const     length,
                             double *const            x,
                             double *const            y)
{
  if (length < 1)
    return;

  double angle = 0.0;
  x[0] = y[0] = EXTERIOR_Y;
  for (int i = 1; i < length; i++) {
    angle = angle - baseInformation[i + 1].angle;

    x[i]  = x[i - 1] + baseInformation[i].distance * cos(angle);
    y[i]  = y[i - 1] + baseInformation[i].distance * sin(angle);
  }
  return;
}


#endif
