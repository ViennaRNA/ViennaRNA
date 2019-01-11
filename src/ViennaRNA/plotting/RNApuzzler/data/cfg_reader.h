#ifndef CFG_READER_H
#define CFG_READER_H

#include <ViennaRNA/plotting/RNApuzzler/dataTypes/configtree_struct.h>
#include <ViennaRNA/plotting/RNApuzzler/dataTypes/tBaseInformation_struct.h>
#include <ViennaRNA/plotting/RNApuzzler/definitions.h>

double getArcAngle(const config *cfg, const int currentArc);
double getArcAngleDegree(const config *cfg, const int currentArc);

/**
 * @brief cfgPrintConfig
 *      - prints the given config to the terminal
 * @param config
 *      - config to print
 */
void cfgPrintConfig(config *config);

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
void cfgGenerateConfig(
    const short *const pair_table,
    tBaseInformation *baseInformation,
    const double unpaired,
    const double paired
);

/**
 * @brief cfgCloneConfig
 *      - prints the given config to the terminal
 * @param cfg
 *      - config to clone
 * @return clone of cfg
 */
config *cfgCloneConfig(const config *cfg);

/**
 * @brief cfgFreeConfig
 *      - prints the given config to the terminal
 * @param cfg
 *      - config to free
 */
void cfgFreeConfig(config *cfg);

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
double cfgApplyChanges(
    config* cfg,
    const char loopName,
    const double* deltaCfg,
    const double radiusNew,
    const puzzlerOptions* puzzler
);

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
short cfgIsValid(
    config* config,
    const double* deltaCfg
);

/**
 * @brief intToMotiv
 *      - converts a given number into a readable format
 * @param _int
 *      - motiv name as number (int)
 * @return
 *      - motiv name as char
 */
char intToMotiv(const int _int);

#endif
