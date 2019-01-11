#include "ViennaRNA/RNApuzzler/resolveIntersections/handleConfigChanges.h"
#include "ViennaRNA/RNApuzzler/dataTypes/configtree_struct.h"
#include "ViennaRNA/RNApuzzler/data/cfg_reader.h"
#include "ViennaRNA/RNApuzzler/data/configtree.h"
#include "ViennaRNA/RNApuzzler/definitions.h"
#include "ViennaRNA/RNApuzzler/resolveIntersections/intersectionType.h"
#include "ViennaRNA/RNApuzzler/output/output.h"

#include <stdlib.h>
#include <math.h>

void logConfigChanges(
    const int id,
    const config *cfg,
    const double* deltaCfg,
    const double oldRadius,
    const double newRadius,
    const char* logTag,
    const puzzlerOptions* puzzler
) {
    const char *fnName = "configChanges";

    printInformation(fnName, "Change #%5d  %s  %5d",
                     puzzler->numberOfChangesAppliedToConfig,
                     logTag,
                     id);
    if (newRadius - oldRadius != 0.0) {
        printInformation(NULL, "  radius: %9.2lf->%9.2lf (%+7.2lf%%) diff: %+.2le",
                 oldRadius, newRadius,
                 100.0 * (newRadius / oldRadius) - 100.0,
                 newRadius - oldRadius);
    } else {
        printInformation(NULL, "  radius did not change");
    }

    if (deltaCfg != NULL) {
        for (int i = 0; i < cfg->numberOfArcs; i++) {
            if (deltaCfg[i] != 0.0) {
                printInformation(NULL, " %d: %+.2le°", i, deltaCfg[i]);
            }
        }
    } else {
        printInformation(NULL, "  EMPTY");
    }
    printInformation(NULL, "\n");
}

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
short checkAndApplyConfigChanges(
        treeNode* tree,
        double* deltaCfg,
        const intersectionType it,
        puzzlerOptions* puzzler
) {
    char* fnName = "checkAndApplyConfigChanges";
    config *cfg = tree->cfg;

    /// fix deltas if changes are too small
    ///
    /// this is necessary because sometimes the calculation results in micro changes.
    /// These micro changes go along with an increase of the loops radius which causes
    /// new problems as the changes being too small to get enough distance to the
    /// changed loop and the intersector being stuck in collision (again).
    ///
    /// multiplying by factor 2.0 we always get a resulting angle between 0.1° and 0.2°
    /// don't use factor 10 as the impact of doing so is way too strong and often causes crashes
    /// in term of applicability of the changes
    short fixTooSmallChanges = 1;
    if (fixTooSmallChanges
        && deltaCfg != NULL) {
        for (int cntr = 0; cntr < 100; cntr++) {
            short valid = 0;
            for (int currentArc = 0; currentArc < cfg->numberOfArcs; currentArc++) {
                if (fabs(deltaCfg[currentArc]) >= epsilon3) {
                    valid = 1;
                    break;
                }
            }
            if (valid) {
                break;
            } else {
                for (int currentArc = 0; currentArc < cfg->numberOfArcs; currentArc++) {
                    deltaCfg[currentArc] = 2.0 * deltaCfg[currentArc];
                }
            }
//        if (LOG_FLAG && cntr > 0) {
//            printf("[ LOG ] fixing... (%d)\n", cntr);
//        }
        }
    }

    /*
    printDebug(fnName, "\t- config old -\n");
    cfgPrintConfig(cfg);
    */

    char* logTag = intersectionTypeToString(it);

    if (cfgIsValid(cfg, deltaCfg)) {

        (puzzler->numberOfChangesAppliedToConfig)++;
        double oldRadius = cfg->radius;

        double radiusNew = -1.0; // == unknown | calculate optimal radius
        applyChangesToConfigAndBoundingBoxes(tree, deltaCfg, radiusNew, puzzler);

        double newRadius = cfg->radius;
        logConfigChanges(getNodeID(tree), cfg, deltaCfg, oldRadius, newRadius, logTag, puzzler);

        return 1;
    } else {
        /// changes result in angles outside 0° to 360°
        printError(fnName, "%s cannot apply changes to config. Invalid changes.\n", logTag);

        /*
        /// print erronious changes
        printConfigError(fnName, tree, deltaCfg);
        */

        /// for not ending up in infinite calculations without being able to apply any changes
        /// we increase the counter for changes per default
        /// infinite calculations occurred while testing with RNA families
        (puzzler->numberOfChangesAppliedToConfig)++;
        return 0;
    }
}

