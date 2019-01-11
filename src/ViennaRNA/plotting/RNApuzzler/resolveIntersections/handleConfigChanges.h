#ifndef HANDLE_CONFIG_CHANGES_H
#define HANDLE_CONFIG_CHANGES_H

#include "ViennaRNA/RNApuzzler/dataTypes/configtree_struct.h"
#include "ViennaRNA/RNApuzzler/definitions.h"
#include "ViennaRNA/RNApuzzler/resolveIntersections/intersectionType.h"

void logConfigChanges(
    const int id,
    const config *cfg,
    const double* deltaCfg,
    const double oldRadius,
    const double newRadius,
    const char* logTag,
    const puzzlerOptions* puzzler
);

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
);

#endif
