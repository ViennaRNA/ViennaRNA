#ifndef VIENNA_RNA_PACKAGE_PROBING_XLINK_H
#define VIENNA_RNA_PACKAGE_PROBING_XLINK_H

#include <ViennaRNA/fold_compound.h>
#include <ViennaRNA/probing/basic.h>

/**
 *  @file     ViennaRNA/probing/crosslink.h
 *  @ingroup  probing_data_strategy
 *  @brief  This module provides function to incorporate RNA structure
 *          probing data derived from cross-linking experiments into the folding recursions
 *          by means of soft constraints and grammar extensions
 */

#define VRNA_PROBING_METHOD_CROSSLINK_LIMIT_PROBABILITY     (1 << 0)
#define VRNA_PROBING_METHOD_CROSSLINK_INPUT_ENERGY          (1 << 1)
#define VRNA_PROBING_METHOD_CROSSLINK_DEFAULT_TEMPERATURE   37.0
#define VRNA_PROBING_METHOD_CROSSLINK_DEFAULT_SCALING       1.0
#define VRNA_PROBING_METHOD_CROSSLINK_DEFAULT_LIMIT         -5.0
#define VRNA_PROBING_METHOD_CROSSLINK_DEFAULT_OPTIONS       0U

/**
 *  @brief  Add crossinking data from an input file
 */
int
vrna_probing_data_xlink_file(vrna_fold_compound_t  *fc,
                             const char            *filename,
                             double                temperature,
                             double                scaling,
                             double                limit,
                             unsigned int          options);

#endif
