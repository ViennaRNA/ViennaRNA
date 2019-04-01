#ifndef RNAPUZZLER_HANDLE_CONFIG_CHANGES_H
#define RNAPUZZLER_HANDLE_CONFIG_CHANGES_H

#include "ViennaRNA/plotting/RNApuzzler/dataTypes/configtree_struct.h"
#include "ViennaRNA/plotting/RNApuzzler/definitions.h"
#include "ViennaRNA/plotting/RNApuzzler/resolveIntersections/intersectionType.h"


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
PRIVATE short checkAndApplyConfigChanges(treeNode               *tree,
                                 double                 *deltaCfg,
                                 const intersectionType it,
                                 puzzlerOptions         *puzzler);


#include "handleConfigChanges.inc"

#endif
