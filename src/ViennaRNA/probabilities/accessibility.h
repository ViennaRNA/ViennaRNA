#ifndef VIENNA_RNA_PACKAGE_EQUILIBRIUM_PROBS_ACCESSIBILITY_H
#define VIENNA_RNA_PACKAGE_EQUILIBRIUM_PROBS_ACCESSIBILITY_H

#include <ViennaRNA/datastructures/basic.h>
#include <ViennaRNA/datastructures/array.h>
#include <ViennaRNA/fold_compound.h>
#include <ViennaRNA/structures/problist.h>
#include <ViennaRNA/params/basic.h>

/**
 *  @file     ViennaRNA/probabilities/accessibility.h
 *  @ingroup  thermodynamics
 *  @brief    Accessibility implementations
 *
 */

/*
 ##############################################
 # ACCESSIBILITY RELATED FUNCTIONS            #
 ##############################################
 */


/**
 *  @addtogroup  thermodynamics
 *  @{
 */


/**
 *  @name Accessibility and derived computations
 *  @{
 */


/**
 *  @brief  Compute Accessibility
 *
 *  @warning  To avoid memory leaks, use the #vrna_array_free() macro to release the memory
 *            for the returned list of element probabilties!
 *
 *  @param  fc      The fold compound to operate on
 *  @param  min_u   The minimum sequence length for unpaired stretches
 *  @param  max_u   The maximum sequence length for unpaired stretches
 *  @param  cutoff  A probability threshold to filter the output
 *  @return         A list of accessibility values in #vrna_array() format
 */
vrna_array(vrna_ep_t)
vrna_eqm_probs_up(vrna_fold_compound_t *fc,
                  unsigned int         min_u,
                  unsigned int         max_u,
                  double               cutoff);

/* End accessibility related functions */
/**@}*/


/* End thermodynamics */
/**@}*/

#endif
