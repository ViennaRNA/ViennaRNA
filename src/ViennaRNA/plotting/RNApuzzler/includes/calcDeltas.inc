#ifndef RNAPUZZLER_CALC_DELTAS_H
#define RNAPUZZLER_CALC_DELTAS_H

/*
 *      RNApuzzler calc deltas
 *
 *      c  Daniel Wiegreffe, Daniel Alexander, Dirk Zeckzer
 *      ViennaRNA package
 */
#include <stdlib.h>
#include <math.h>

#include "ViennaRNA/utils/basic.h"

#include "definitions.inc"
#include "boundingWedge.inc"
#include "../headers/configtree_struct.h"
#include "configtree.inc"
#include "drawingconfig.inc"
#include "vector_math.inc"

/**
 * @brief calcDeltas
 *      The area between stems indexLeft and indexRight
 *      (by traversing the loop clockwise starting at indexLeft-stem)
 *      will be enlarged in degree as given via deltaAngle.
 *      All other areas will be used to compensate that increase
 *      (i.e. by decreasing those area's angles).
 * @param node
 * @param recursiveEnd
 * @param indexLeft
 * @param indexRight
 * @param deltaAngle
 * @param deltas
 * @return the amount of change (in positive degree) that can be accomplished with calculated deltas
 */
PRIVATE double
calcDeltas(const treeNode               *node,
           const treeNode               *recursiveEnd,
           const int                    indexLeft,
           const int                    indexRight,
           const double                 deltaAngle,
           vrna_plot_options_puzzler_t  *puzzler,
           double                       *deltas);


PRIVATE void
calcDeltasEquidistantIncrease(const double  targetAngleIn,
                              const int     configSize,
                              short         *increase,
                              double        *deltaCfg)
{
  char    *fnName     = "CALC DELTAS EQUIDISTANT INCREASE";
  double  targetAngle = targetAngleIn;

  int     increaseCount = 0;

  for (int i = 0; i < configSize; i++) {
    if (increase[i])

      increaseCount++;
  }
  double deltaPerIncrease = targetAngle / increaseCount;

  for (int i = 0; i < configSize; i++) {
    if (increase[i])

      deltaCfg[i] += deltaPerIncrease;
  }
}


PRIVATE double
calcDeltasMaximumFirstDecrease(const double targetAngleIn,
                               const int    indexLeft,
                               const int    indexRight,
                               const int    configSize,
                               double       *deltaCfg,
                               double       *currentAngles,
                               const double minAngleHalf)
{
  char    *fnName     = "CALC DELTAS MAXIMUM FIRST DECREASE";
  double  targetAngle = targetAngleIn;

  int     i;

  short   doLoop = 1;

  while (doLoop) {
    double  maxSpace      = 0.0;
    int     maxSpaceIndex = -1;

    if (indexLeft == -1) {
      double sumAngles = 0.0;
      i = -1;

      while (i != indexRight) {
        i++;
        double cfg = currentAngles[i] + deltaCfg[i] - 2 * minAngleHalf;
        sumAngles += cfg;
      }

      while (i != configSize - 1) {
        i++;
        double cfg = currentAngles[i] + deltaCfg[i] - 2 * minAngleHalf;

        if (sumAngles < MATH_PI) {
          if (cfg > maxSpace) {
            maxSpace      = cfg;
            maxSpaceIndex = i;
          }
        } else {
          break;
        }

        /* sum increase happens afterwards to allow bending the arc containing 180° */
        sumAngles += cfg;
      }
    } else if (indexRight == -1) {
      double sumAngles = 0.0;
      i = configSize - 1;

      while (i != indexLeft) {
        double cfg = currentAngles[i] + deltaCfg[i] - 2 * minAngleHalf;
        sumAngles += cfg;
        i--;
      }

      while (i != -1) {
        double cfg = currentAngles[i] + deltaCfg[i] - 2 * minAngleHalf;

        if (sumAngles < MATH_PI) {
          if (cfg > maxSpace) {
            maxSpace      = cfg;
            maxSpaceIndex = i;
          }
        } else {
          break;
        }

        /* sum increase happens afterwards to allow bending the arc containing 180° */
        sumAngles += cfg;
        i--;
      }
    } else {
      i = indexRight;
      if (i == configSize - 1)
        i = -1;

      while (i != indexLeft) {
        double cfg = currentAngles[i + 1] + deltaCfg[i + 1] - 2 * minAngleHalf;

        if (cfg > maxSpace) {
          maxSpace      = cfg;
          maxSpaceIndex = i + 1;
        }

        i++;

        if (i == configSize - 1)
          i = -1;
      }
    }

    double diff = 0.0;
    if (maxSpaceIndex != -1) {
      double factor = (targetAngle < 0.1 * targetAngleIn) ? 1.0 : 0.5;
      diff = (-1) * fmin(factor * maxSpace, targetAngle);

      deltaCfg[maxSpaceIndex] += diff;
      targetAngle             += diff;
    }

    doLoop = targetAngle > 0.0 && fabs(diff) > EPSILON_3;
  }


  return targetAngle;
}


PRIVATE double
calcDeltasNearestNeighborsFirstDecrease(const double  targetAngleIn,
                                        const int     indexLeft,
                                        const int     indexRight,
                                        const int     configSize,
                                        short         *decrease,
                                        double        *space,
                                        double        *deltaCfg)
{
  char    *fnName     = "CALC DELTAS NEAREST NEIGHBOR FIRST DECREASE";
  double  targetAngle = targetAngleIn;

  /* count the number of possible iteration steps */
  int     startIndex = indexRight + 1;

  if (startIndex == configSize)
    startIndex = -1;

  int     stopIndex = indexLeft + 1;
  if (stopIndex == configSize)
    stopIndex = -1;

  int     steps   = 0;
  int     stemIt  = indexRight;
  while (stemIt != indexLeft) {
    stemIt++;
    if (stemIt == configSize)

      stemIt = -1;

    steps++;
  }
  int   numIt = steps / 2; /* implicit floor() operation */

  int   *index  = (int *)vrna_alloc(steps * sizeof(int));
  short changed = 1;
  while (changed) {
    changed = 0;
    int count = 0;

    int iL = indexLeft;

    if (iL == -1)

      iL = configSize - 1;

    int iR = indexRight + 1;

    if (iR == configSize)

      iR = 0;

    for (int i = 0; i < numIt; i++) {
      if (decrease[iL]) {
        index[count] = iL;
        count++;
      }

      if (decrease[iR]) {
        index[count] = iR;
        count++;
      }

      iL--;

      if (iL == -1)

        iL = configSize - 1;

      iR++;

      if (iR == configSize)
        iR = 0;
    }

    if (numIt < 0.5 * steps) {
      index[count] = iL;
      count++;
      iL--;

      if (iL == -1)
        iL = configSize - 1;
    }

    if (count > 0) {
      double partAngle = targetAngle / count;
      for (int k = 0; k < count; k++) {
        int j = index[k];

        if (decrease[j]) {
          double diff = (-1) * fmin(space[j] + deltaCfg[j], partAngle);
          deltaCfg[j] += diff;
          targetAngle += diff;
          changed     = changed || (diff != 0.0);
        }
      }
    }
  }
  free(index);

  return targetAngle;
}


/**
 * @brief calcDeltas
 *      The area between stems indexLeft and indexRight
 *      (by traversing the loop clockwise starting at indexLeft-stem)
 *      will be enlarged in degree as given via deltaAngle.
 *      All other areas will be used to compensate that increase
 *      (i.e. by decreasing those area's angles).
 * @param node
 * @param recursiveEnd
 * @param indexLeft
 * @param indexRight
 * @param deltaAngle
 * @param deltas
 * @return the amount of change (in positive degree) that can be accomplished with calculated deltas
 */
PRIVATE double
calcDeltas(const treeNode               *node,
           const treeNode               *recursiveEnd,
           const int                    indexLeft,
           const int                    indexRight,
           const double                 deltaAngle,
           vrna_plot_options_puzzler_t  *puzzler,
           double                       *deltas)
{
  char *fnName = "CALC DELTAS";

  /* Check: valid angle >= 0.0 */
  if (deltaAngle < 0.0)

    return 0.0;

  int     childCount  = node->childCount;
  int     configSize  = childCount + 1;

  /* get the current node's stem's bounding wedge */
  double  minOuterAngle = asin(puzzler->paired / (2 * node->cfg->radius));

  /* allocate memory for stuff used in calculation */
  double  *anglesMin      = (double *)vrna_alloc(childCount * sizeof(double));
  double  *anglesMax      = (double *)vrna_alloc(childCount * sizeof(double));
  double  *space          = (double *)vrna_alloc(configSize * sizeof(double));
  double  *deltaCfg       = (double *)vrna_alloc(configSize * sizeof(double));
  short   *increase       = (short *)vrna_alloc(configSize * sizeof(short));
  short   *decrease       = (short *)vrna_alloc(configSize * sizeof(short));
  double  *currentAngles  = (double *)vrna_alloc(configSize * sizeof(double));

  /* Initialization currentAngles */
  config  *cfg = node->cfg;

  for (int currentArc = 0; currentArc < cfg->numberOfArcs; ++currentArc)

    currentAngles[currentArc] = getArcAngle(cfg, currentArc);

  /* get all bounding wedges (minAngle, maxAngle) */
  double min, max;

  for (int currentChild = 0; currentChild < childCount; currentChild++) {
    getBoundingWedge(node, currentChild, &min, &max);
    anglesMin[currentChild] = min;
    anglesMax[currentChild] = max;
  }

  /* convert bounding wedges to "free" areas that can be used for compensation of changes */
  space[0] = anglesMin[0] - (0 + minOuterAngle);

  for (int i = 1; i < (configSize - 1); i++)

    space[i] = anglesMin[i] - anglesMax[i - 1];

  space[configSize - 1] = (MATH_TWO_PI - minOuterAngle) - anglesMax[configSize - 2];

  /* fix too big spaces (may become bigger than config for very large loops) */
  for (int i = 0; i < configSize; i++)

    space[i] = fmin(space[i], getArcAngle(node->cfg, i) - 2 * minOuterAngle);

  /* Initialization: calculation values (deltaCfg, increase, decrease) */
  for (int i = 0; i < configSize; i++) {
    deltaCfg[i] = 0.0;
    increase[i] = -1;
    decrease[i] = -1;
  }

  /* Mark increase and decrease areas */
  int currentIndex = indexLeft;   /* stemIndex */
  while (currentIndex != indexRight) {
    increase[currentIndex + 1]  = 1;
    decrease[currentIndex + 1]  = 0;
    currentIndex++;

    if (currentIndex == configSize - 1)
      currentIndex = -1;
  }

  while (currentIndex != indexLeft) {
    increase[currentIndex + 1]  = 0;
    decrease[currentIndex + 1]  = (space[currentIndex + 1] > 0.0);
    currentIndex++;

    if (currentIndex == configSize - 1)

      currentIndex = -1;
  }

  double targetAngle = deltaAngle;

  /* Step 1: equidistant increase */
  calcDeltasEquidistantIncrease(targetAngle, configSize, increase, deltaCfg);

  /* Step 2: nearest neighbor first decrease */
  targetAngle = calcDeltasNearestNeighborsFirstDecrease(targetAngle,
                                                        indexLeft,
                                                        indexRight,
                                                        configSize,
                                                        decrease,
                                                        space,
                                                        deltaCfg);

  /* Step 3: check if intersections are fixed */
  short notFixedYet = (targetAngle != 0.0);

  if (notFixedYet) {
    /*
     * if the intersection is not yet fixed
     * check if there is a loop on a higher level that can be bend instead of this one
     * if this is the case we can apply methods using spaces
     * otherwise we need to use drastical measures (e.g. maximumFirstDecrease using cfg instead of spaces)
     */

    treeNode  *parent     = getParent(node);
    short     canGoHigher = 0;

    while (parent != recursiveEnd && !isExterior(parent)) {
      short parentIsMultiLoop = isMultiLoop(parent);
      if (parentIsMultiLoop) {
        canGoHigher = 1;
        break;
      } else {
        /* parent is internal loop */
        double childAngle = getArcAngle(parent->cfg, 0);

        if (fabs(childAngle - MATH_PI) < EPSILON_3) {
          /* no op */
        } else if (childAngle > MATH_PI) {
          if (indexLeft == 0) {
            canGoHigher = 1;
            break;
          }
        } else if (childAngle < MATH_PI) {
          if (indexLeft == -1) {
            canGoHigher = 1;
            break;
          }
        }
      }

      /* if current parent node can not be adapted check its parent */
      parent = getParent(parent);
    }

    if (!canGoHigher) {
      targetAngle = calcDeltasMaximumFirstDecrease(targetAngle,
                                                   indexLeft,
                                                   indexRight,
                                                   configSize,
                                                   deltaCfg,
                                                   currentAngles,
                                                   minOuterAngle);
    }
  }

  /* Step 4: equidistant increase with negative remaining target angle */
  calcDeltasEquidistantIncrease((-1) * targetAngle, configSize, increase, deltaCfg);

  /*
   * Fix deltas if changes are too small.
   * This is necessary because sometimes the calculation results in micro changes.
   * These micro changes go along with an increase of the loops radius which causes
   * new problems as the changes being too small to get enough distance to the
   * changed loop and the intersector being stuck in collision (again).
   *
   * multiplying by factor 2.0 we always get a resulting angle between 0.1° and 0.2°
   * don't use factor 10 as the impact of doing so is way too strong and often causes crashes
   * in term of applicability of the changes
   */
  short fixTooSmallChanges = 0;
  if (fixTooSmallChanges) {
    for (int cntr = 0; cntr < 100; cntr++) {
      short valid = 0;

      for (int currentArc = 0; currentArc < configSize; currentArc++) {
        if (fabs(deltaCfg[currentArc]) >= EPSILON_3) {
          valid = 1;
          break;
        }
      }
      if (valid) {
        break;
      } else {
        for (int currentArc = 0; currentArc < configSize; currentArc++)

          deltaCfg[currentArc] = 2.0 * deltaCfg[currentArc];
      }
    }
  }

  /* transfer calculated deltas to return area */
  for (int currentArc = 0; currentArc < configSize; currentArc++)

    deltas[currentArc] = deltaCfg[currentArc];

  /* free allocated memory */
  free(anglesMin);
  free(anglesMax);
  free(space);
  free(deltaCfg);
  free(increase);
  free(decrease);
  free(currentAngles);

  /* check if all deltas sum up to zero */
  double checkSum = 0.0;

  for (int currentArc = 0; currentArc < configSize; currentArc++)

    checkSum += deltas[currentArc];

  if (fabs(checkSum) > EPSILON_3) {
    for (int currentArc = 0; currentArc < configSize; currentArc++)

      deltas[currentArc] = 0.0;
    targetAngle = deltaAngle;
  }

  if (!cfgIsValid(cfg, deltas)) {
    for (int currentArc = 0; currentArc < configSize; currentArc++)

      deltas[currentArc] = 0.0;
    targetAngle = deltaAngle;
  }

  /* return the difference that can be accomplished using these deltas */
  double changedAngle = deltaAngle - targetAngle;
  return changedAngle;
}


#endif
