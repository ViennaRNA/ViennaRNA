#ifndef RNAPUZZLER_CONFIG_H
#define RNAPUZZLER_CONFIG_H

/*
 *      RNApuzzler template config
 *
 *      c  Daniel Wiegreffe, Daniel Alexander, Dirk Zeckzer
 *      ViennaRNA package
 */
#include <stdlib.h>
#include <math.h>

#include "ViennaRNA/utils/basic.h"

#include "definitions.inc"
#include "vector_math.inc"
#include "../headers/configtree_struct.h"
#include "../headers/tBaseInformation_struct.h"
#include "ViennaRNA/plotting/RNApuzzler/RNApuzzler.h"

PRIVATE double
getArcAngle(const config  *cfg,
            const int     currentArc);


PRIVATE double
getArcAngleDegree(const config  *cfg,
                  const int     currentArc);


/**
 * @brief cfgPrintConfig
 *      - prints the given config to the terminal
 * @param config
 *      - config to print
 */
PRIVATE void
cfgPrintConfig(config *config);


/**
 * @brief cfgGenerateConfig
 *      - generates configurations for each loop
 * @param pair_table
 *      - the RNA's pairing information
 * @param baseInformation
 *      - array of tBaseInformation annotations (for saving configs)
 * @param unpaired
 *      - default length of backbone lines
 * @param paired
 *      - default distance between paired bases
 */
PRIVATE void
cfgGenerateConfig(const short *const  pair_table,
                  tBaseInformation    *baseInformation,
                  const double        unpaired,
                  const double        paired);


/**
 * @brief cfgCloneConfig
 *      - prints the given config to the terminal
 * @param cfg
 *      - config to clone
 * @return clone of cfg
 */
PRIVATE config *
cfgCloneConfig(const config *cfg);


/**
 * @brief cfgFreeConfig
 *      - prints the given config to the terminal
 * @param cfg
 *      - config to free
 */
PRIVATE void
cfgFreeConfig(config *cfg);


/**
 * Function to apply a set of config changes.
 *
 * By passing -1.0 to radiusNew you will apply the minimum possible radius for this loop while not allowing to shrink the loop.
 * In case it would actually shrink a default relative increase is applied to the radius.
 *
 * By passing 0.0 to radiusNew you set the minimum possible radius no matter if there would be some shrinking or not.
 * @param loopName the node to change
 * @param deltaCfg array of angles to change (in degree)
 * @param radiusNew desired radius to set while applying deltas
 * @param puzzler
 */
PRIVATE double
cfgApplyChanges(config                            *cfg,
                const char                        loopName,
                const double                      *deltaCfg,
                const double                      radiusNew,
                const vrna_plot_options_puzzler_t *puzzler);


/**
 * @brief cfgIsValid
 *      - check if config is valid
 * @param config
 *      - the config to check
 * @param deltaCfg
 *      - array of config changes.
 *        contains diff-values for each config angle.
 *        in degree format
 * @return true iff config is valid
 */
PRIVATE short
cfgIsValid(config       *config,
           const double *deltaCfg);


/**
 * @brief intToMotiv
 *      - converts a given number into a readable format
 * @param _int
 *      - motiv name as number (int)
 * @return
 *      - motiv name as char
 */
PRIVATE char
intToMotiv(const int _int);


/*
 *--------------------------------------------------------------------------
 *---   create, copy, and free config   ------------------------------------
 *--------------------------------------------------------------------------
 */

/**
 * @brief cfgCreateConfig
 *      - constructor-like method for creating a config
 * @param radius
 *      - radius used for drawing that loop
 * @return
 *      - an initialized config struct
 */
PRIVATE config *
cfgCreateConfig(const double radius)
{
  config *cfg = (config *)vrna_alloc(1 * sizeof(config));

  cfg->radius         = radius;
  cfg->minRadius      = radius;
  cfg->defaultRadius  = radius;

  cfg->cfgArcs      = NULL;
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
PRIVATE configArc
cfgCreateConfigArc(const double angle,
                   const int    numberOfArcSegments)
{
  configArc newConfigArc;

  newConfigArc.numberOfArcSegments = numberOfArcSegments;

  newConfigArc.arcAngle = angle;

  return newConfigArc;
}


PRIVATE config *
cfgCloneConfig(const config *cfg)
{
  config *clonedCfg = (config *)vrna_alloc(sizeof(config));

  clonedCfg->radius         = cfg->radius;
  clonedCfg->minRadius      = cfg->minRadius;
  clonedCfg->defaultRadius  = cfg->defaultRadius;
  clonedCfg->numberOfArcs   = cfg->numberOfArcs;

  const int numberOfArcs = cfg->numberOfArcs;
  clonedCfg->cfgArcs = (configArc *)vrna_alloc(numberOfArcs * sizeof(configArc));

  for (int currentArc = 0; currentArc < numberOfArcs; ++currentArc) {
    clonedCfg->cfgArcs[currentArc].numberOfArcSegments =
      cfg->cfgArcs[currentArc].numberOfArcSegments;
    clonedCfg->cfgArcs[currentArc].arcAngle = cfg->cfgArcs[currentArc].arcAngle;
  }

  return clonedCfg;
}


PRIVATE void
cfgFreeConfig(config *cfg)
{
  free(cfg->cfgArcs);
  free(cfg);
}


/* documentation at header file */
PRIVATE void
cfgPrintConfig(config *cfg)
{
  short useDegree = 1;

  for (int currentArc = 0; currentArc < cfg->numberOfArcs; ++currentArc) {
    double angle = getArcAngle(cfg, currentArc);

    if (useDegree)
      angle = toDegree(angle);
  }
}


/*
 *--------------------------------------------------------------------------
 *---   access to config elements   ----------------------------------------
 *--------------------------------------------------------------------------
 */
PRIVATE double
getArcAngle(const config  *cfg,
            const int     currentArc)
{
  return cfg->cfgArcs[currentArc].arcAngle;
}


PRIVATE double
getArcAngleDegree(const config  *cfg,
                  const int     currentArc)
{
  return toDegree(getArcAngle(cfg, currentArc));
}


/*
 *--------------------------------------------------------------------------
 *---   radius computation   -----------------------------------------------
 *--------------------------------------------------------------------------
 */

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
PRIVATE double
approximateConfigArcRadius(const double a,
                           const double b,
                           const short  m,
                           const short  n,
                           double       angle)
{
  const short MAX_ITERATIONS = 1000;

  /*
   * calculation:
   *
   * be s the length of a line at the circle (paired or unpaired / a or b)
   * the angle over such a single line be alpha
   * alpha = angle / ( m + n )
   *
   * for such a single line the following equation holds (where r is the radius of the circle)
   * sin( alpha / 2 ) = ( s / 2 ) / r
   * r = ( s / 2 ) / sin( alpha / 2 )
   * r = ( s / 2 ) / sin( ( angle / ( m + n ) ) / 2 )
   *
   * now we replace s with a or b to get the upper or lower bound for the radius interval
   */
  double  lowerBound  = (b / 2) / sin((angle / (m + n)) / 2);
  double  upperBound  = (a / 2) / sin((angle / (m + n)) / 2);
  /* printf("new: lower %f upper %f\n", lowerBound, upperBound); */
  double  rtn = 0.5 * (lowerBound + upperBound);

  /*
   * there is a minimum valid radius!
   * if rtn is smaller than 0.5*a or 0.5*b than the result will become nan
   */
  rtn = fmax(rtn, 0.5 * a);
  rtn = fmax(rtn, 0.5 * b);

  /*
   * printf("lower %f upper %f %f\n", lowerBound, upperBound, rtn);
   * printf("Angle %f\n", angle);
   */

  /* short isERROR = 0; */

  int j = 0;

  for (j = 0; j < MAX_ITERATIONS; j++) {
    double dx = 2 * (m * asin(a / (2 * rtn)) + n * asin(b / (2 * rtn)) - (angle / 2))
                / (-(a * m / (rtn * sqrt(rtn * rtn - a * a / 4)) + b * n /
                     (rtn * sqrt(rtn * rtn - b * b / 4))));
    rtn -= dx;

    if (fabs(dx) < EPSILON_3)

      break;
  }

  if (rtn < lowerBound)
    rtn = lowerBound;
  else if (rtn > upperBound)
    rtn = upperBound;

  return rtn;
}


PRIVATE double
approximateConfigRadius(const config  *cfg,
                        const double  unpaired,
                        const double  paired)
{
  /*
   * calculate a fitting radius for each arc without compressing or stretching arc segments
   * return the maximum of those values as the radius fitting for the loop
   */
  double r = 0;

  for (int currentArc = 0; currentArc < cfg->numberOfArcs; ++currentArc) {
    int     stems               = 1;
    int     numberOfArcSegments = (cfg->cfgArcs[currentArc]).numberOfArcSegments;
    double  angle               = getArcAngle(cfg, currentArc);

    double  tempR = approximateConfigArcRadius(paired, unpaired, stems, numberOfArcSegments, angle);

    if (tempR > r)

      r = tempR;
  }

  return r;
}


/* -------------------------------------------------------------------------------------------------------------------------- */

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
PRIVATE config *
cfgGenerateDefaultConfig(const short *const pair_table,
                         const int          start,
                         const int          unpaired,
                         const int          paired,
                         const double       radius)
{
  /* create loop configuration */
  config  *cfg = cfgCreateConfig(radius);

  /* compute angles for paired and unpaired bases */
  double  anglePaired   = 2 * asin(paired / (2 * radius));        /* angle over paired */
  double  angleUnpaired = 2 * asin(unpaired / (2 * radius));      /* angle over unpaired */

  /* initialize values for first arc */
  int     arcUnpaired = 0;
  double  angleArc;  /* alpha + numBackbones * beta */

  /* start with first base after parent stem */
  int     i = start + 1;

  while (i <= pair_table[start]) {
    /* until last base at parent stem */
    if (pair_table[i] == 0) {
      /* arc */
      i++;
    } else {
      /* increment number of arcs */
      ++(cfg->numberOfArcs);

      if (i != pair_table[start])
        /* skip subtree at stem */
        i = pair_table[i] + 1;
      else
        /* parent stem -> finish */
        break;
    }
  }

  cfg->cfgArcs = (configArc *)vrna_alloc(cfg->numberOfArcs * sizeof(configArc));

  /* start with first base after parent stem */
  i = start + 1;
  int currentArc          = 0;
  int numberOfArcSegments = 0;

  while (i <= pair_table[start]) {
    /* until last base at parent stem */
    if (pair_table[i] == 0) {
      /* arc */
      arcUnpaired++;
      i++;
    } else {
      /* stem: create arc */
      angleArc                  = anglePaired + (arcUnpaired + 1) * angleUnpaired;
      numberOfArcSegments       = arcUnpaired + 1;
      cfg->cfgArcs[currentArc]  = cfgCreateConfigArc(angleArc, numberOfArcSegments);
      ++currentArc;

      if (i != pair_table[start]) {
        /* initialize values for next arc */
        arcUnpaired = 0;

        /* skip subtree at stem */
        i = pair_table[i] + 1;
      } else {
        /* parent stem -> finish */
        break;
      }
    }
  }

  return cfg;
}


PRIVATE void
cfgGenHandleStem(int                baseNr,
                 const short *const pair_table,
                 tBaseInformation   *baseInformation,
                 const double       unpaired,
                 const double       paired);


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
PRIVATE void
cfgGenHandleLoop(int                baseNr,
                 const short *const pair_table,
                 tBaseInformation   *baseInformation,
                 const double       unpaired,
                 const double       paired)
{
  int start = baseNr;
  int end   = pair_table[baseNr];

  int unpairedCount = 0;
  int stemCount     = 1;

  /* count stems and unpaired bases to use for bulge detection */
  int i = start + 1;

  while (i < end) {
    if (pair_table[i] == 0) {
      /* unpaired base */
      unpairedCount++;
      i++;
    } else if (pair_table[i] > i) {
      /* found new stem */
      stemCount++;
      i = pair_table[i];
    } else {
      /* returned from stem */
      i++;
    }
  }

  short isBulge = (stemCount == 2 && unpairedCount == 1);
  if (isBulge) {
    if (pair_table[start + 1] == 0)
      /* unpaired on left strand */
      cfgGenHandleStem(start + 2, pair_table, baseInformation, unpaired, paired);
    else
      /* unpaired on the right strand */
      cfgGenHandleStem(start + 1, pair_table, baseInformation, unpaired, paired);
  } else {
    int     m             = stemCount;                  /* compare RNApuzzler.c -> f_handle_loop */
    int     n             = unpairedCount + stemCount;  /* compare RNApuzzler.c -> f_handle_loop */
    double  defaultRadius = approximateConfigArcRadius(paired, unpaired, m, n, MATH_TWO_PI);
    config  *cfgLoop      = cfgGenerateDefaultConfig(pair_table,
                                                     start,
                                                     unpaired,
                                                     paired,
                                                     defaultRadius);
    baseInformation[start].config = cfgLoop;

    int i = start + 1;
    while (i < end) {
      if (pair_table[i] == 0) {
        /* unpaired base */
        i++;
      } else if (pair_table[i] > i) {
        /* found new stem */
        cfgGenHandleStem(i, pair_table, baseInformation, unpaired, paired);
        i = pair_table[i];
      } else {
        /* returned from stem */
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
PRIVATE void
cfgGenHandleStem(int                baseNr,
                 const short *const pair_table,
                 tBaseInformation   *baseInformation,
                 const double       unpaired,
                 const double       paired)
{
  /* iterating over the stem and calling cfgGenHandleLoop as soon as a loop is found */
  short continueStem  = 1;
  int   i             = baseNr;

  while (continueStem) {
    if (pair_table[i + 1] == pair_table[i] - 1) {
      i++;
    } else {
      /* found unpaired above stem */
      cfgGenHandleLoop(i, pair_table, baseInformation, unpaired, paired);
      continueStem = 0;
    }
  }
}


/* documentation at header file */
PRIVATE void
cfgGenerateConfig(const short *const  pair_table,
                  tBaseInformation    *baseInformation,
                  const double        unpaired,
                  const double        paired)
{
  /*
   * iterate over RNA
   * for each loop generate a default config
   */
  int length  = pair_table[0];
  int i       = 1;

  while (i < length) {
    if (pair_table[i] == 0) {
      /* unpaired at exterior loop */
      i++;
    } else if (pair_table[i] > i) {
      /* found stem */
      cfgGenHandleStem(i, pair_table, baseInformation, unpaired, paired);
      i = pair_table[i];
    } else {
      /* returned from stem */
      i++;
    }
  }
}


/*
 *--------------------------------------------------------------------------
 *---   set and update config elements   -----------------------------------
 *--------------------------------------------------------------------------
 */

/* documentation at header file */
/**
 * @brief cfgSetRadius
 *      - changes the value of radius for a config to the given value
 * @param config
 *      - config that is being altered
 * @param radius
 *      - new radius
 */
PRIVATE void
cfgSetRadius(config       *cfg,
             const double radius)
{
  cfg->radius = radius;
}


/**
 * @brief cfgUpdateMinRadius
 *      - updates the minimum possible radius for the given config
 * @param config
 * @param unpaired
 * @param paired
 */
PRIVATE void
cfgUpdateMinRadius(config       *cfg,
                   const double unpaired,
                   const double paired)
{
  double minRadius = approximateConfigRadius(cfg, unpaired, paired);

  cfg->minRadius = minRadius;
}


PRIVATE double
cfgApplyChanges(config                            *cfg,
                const char                        loopName,
                const double                      *deltaCfg,
                const double                      radiusNew,
                const vrna_plot_options_puzzler_t *puzzler)
{
  /* start with adjusting config angles; if applicable */
  if (deltaCfg != NULL) {
    for (int currentArc = 0; currentArc < cfg->numberOfArcs; currentArc++)

      (cfg->cfgArcs[currentArc]).arcAngle += deltaCfg[currentArc];
  }

  /* then, adjust config radius */
  double  oldRadius = cfg->radius;
  double  newRadius = -1.0;
  if (radiusNew > 0.0) {
    /*
     * in case the input is a positive value
     * we set the minimum of valid and input as new radius
     */
    cfgUpdateMinRadius(cfg, puzzler->unpaired, puzzler->paired);
    newRadius = fmax(radiusNew, cfg->minRadius);
    cfgSetRadius(cfg, newRadius);
  } else if (radiusNew == 0.0) {
    /*
     * set the minRadius as new value
     * (this allows to shrink a loop)
     */
    cfgUpdateMinRadius(cfg, puzzler->unpaired, puzzler->paired);
    newRadius = cfg->minRadius;
    cfgSetRadius(cfg, newRadius);
  } else if (radiusNew == -1.0) {
    /*
     * set the minRadius as new value
     * (this forbidds to shrink a loop)
     */
    cfgUpdateMinRadius(cfg, puzzler->unpaired, puzzler->paired);

    if (cfg->minRadius - EPSILON_0 > oldRadius) {
      newRadius = cfg->minRadius;
    } else {
      double defaultIncrease = 1.05;
      newRadius = oldRadius * defaultIncrease;
    }

    cfgSetRadius(cfg, newRadius);
  } else {
    /* all unhandled inputs result in errors */
    newRadius = -1.0;
  }

  return newRadius;
}


PRIVATE short
cfgIsValid(config       *cfg,
           const double *deltaCfg)
{
  if (deltaCfg == NULL)

    return 0;

  double  sumAngles         = 0.0;
  short   validSingleAngles = 1;

  for (int currentArc = 0; currentArc < cfg->numberOfArcs; currentArc++) {
    double  angle = getArcAngle(cfg, currentArc) + deltaCfg[currentArc];
    sumAngles += angle;

    short   validAngle = 0.0 < angle && angle < MATH_TWO_PI;
    validSingleAngles = validSingleAngles && validAngle;
  }

  short validSumAngles = (fabs(sumAngles - MATH_TWO_PI) < EPSILON_3);

  return validSingleAngles && validSumAngles;
}


#endif
