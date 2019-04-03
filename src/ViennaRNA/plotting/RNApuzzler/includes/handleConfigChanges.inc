#ifndef RNAPUZZLER_HANDLE_CONFIG_CHANGES_H
#define RNAPUZZLER_HANDLE_CONFIG_CHANGES_H

/*
 *      RNApuzzler handle config changes
 *
 *      c  Daniel Wiegreffe, Daniel Alexander, Dirk Zeckzer
 *      ViennaRNA package
 */

#include <stdlib.h>
#include <math.h>

#include "ViennaRNA/utils/basic.h"

#include "../headers/configtree_struct.h"
#include "drawingconfig.inc"
#include "configtree.inc"
#include "definitions.inc"
#include "intersectionType.inc"

/**
 * @brief checkAndApplyConfigChanges
 *      - Method for performing of config.
 *        Alters config as well as all corresponding boundingboxes.
 *        Determines the new radius that fits best.
 * @param tree
 *      - tree node where the config is changed.
 * @param deltaCfg
 *      - array of config changes.
 *        contains diff-values for each config angle.
 *        in degree format
 * @return 1 if something changed, 0 otherwise
 */
PRIVATE short
checkAndApplyConfigChanges(treeNode                     *tree,
                           double                       *deltaCfg,
                           const intersectionType       it,
                           vrna_plot_options_puzzler_t  *puzzler);


/**
 * @brief checkAndApplyConfigChanges
 *      - Method for performing of config.
 *        Alters config as well as all corresponding boundingboxes.
 *        Determines the new radius that fits best.
 * @param tree
 *      - tree node where the config is changed.
 * @param deltaCfg
 *      - array of config changes.
 *        contains diff-values for each config angle.
 *        in degree format
 * @return 1 if something changed, 0 otherwise
 */
PRIVATE short
checkAndApplyConfigChanges(treeNode                     *tree,
                           double                       *deltaCfg,
                           const intersectionType       it,
                           vrna_plot_options_puzzler_t  *puzzler)
{
  char    *fnName = "checkAndApplyConfigChanges";
  config  *cfg    = tree->cfg;

  /*
   * fix deltas if changes are too small
   *
   * this is necessary because sometimes the calculation results in micro changes.
   * These micro changes go along with an increase of the loops radius which causes
   * new problems as the changes being too small to get enough distance to the
   * changed loop and the intersector being stuck in collision (again).
   *
   * multiplying by factor 2.0 we always get a resulting angle between 0.1° and 0.2°
   * don't use factor 10 as the impact of doing so is way too strong and often causes crashes
   * in term of applicability of the changes
   */
  short fixTooSmallChanges = 1;

  if (fixTooSmallChanges
      && deltaCfg != NULL) {
    for (int cntr = 0; cntr < 100; cntr++) {
      short valid = 0;

      for (int currentArc = 0; currentArc < cfg->numberOfArcs; currentArc++) {
        if (fabs(deltaCfg[currentArc]) >= EPSILON_3) {
          valid = 1;
          break;
        }
      }

      if (valid) {
        break;
      } else {
        for (int currentArc = 0; currentArc < cfg->numberOfArcs; currentArc++)

          deltaCfg[currentArc] = 2.0 * deltaCfg[currentArc];
      }
    }
  }

  char *logTag = intersectionTypeToString(it);

  if (cfgIsValid(cfg, deltaCfg)) {
    (puzzler->numberOfChangesAppliedToConfig)++;
    double  oldRadius = cfg->radius;

    double  radiusNew = -1.0;    /* == unknown | calculate optimal radius */
    applyChangesToConfigAndBoundingBoxes(tree, deltaCfg, radiusNew, puzzler);

    double  newRadius = cfg->radius;

    return 1;
  } else {
    /*
     * changes result in angles outside 0° to 360°
     * for not ending up in infinite calculations without being able to apply any changes
     * we increase the counter for changes per default
     * infinite calculations occurred while testing with RNA families
     */
    (puzzler->numberOfChangesAppliedToConfig)++;
    return 0;
  }
}


#endif
