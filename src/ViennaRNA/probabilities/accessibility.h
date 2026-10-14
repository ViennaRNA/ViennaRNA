#ifndef VIENNA_RNA_PACKAGE_EQUILIBRIUM_PROBS_ACCESSIBILITY_H
#define VIENNA_RNA_PACKAGE_EQUILIBRIUM_PROBS_ACCESSIBILITY_H

#include <ViennaRNA/datastructures/basic.h>
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
