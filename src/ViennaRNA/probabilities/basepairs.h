#ifndef VIENNA_RNA_PACKAGE_EQUILIBRIUM_PROBS_BASEPAIRS_H
#define VIENNA_RNA_PACKAGE_EQUILIBRIUM_PROBS_BASEPAIRS_H

#include <ViennaRNA/datastructures/basic.h>
#include <ViennaRNA/fold_compound.h>
#include <ViennaRNA/structures/problist.h>
#include <ViennaRNA/params/basic.h>

/**
 *  @file     ViennaRNA/probabilities/basepairs.h
 *  @ingroup  thermodynamics
 *  @brief    Equilibrium Probability implementations
 *
 *  This file includes various implementations for equilibrium
 *  probability computations based on the partition function
 *  of an RNA sequence, two concatenated sequences, or a sequence
 *  alignment.
 */

/*
 #################################################
 # BASE PAIR PROBABILITY RELATED FUNCTIONS       #
 #################################################
 */


/**
 *  @addtogroup  thermodynamics
 *  @{
 */


/**
 *  @name Base pair probabilities and derived computations
 *  @{
 */

int
vrna_pairing_probs(vrna_fold_compound_t *fc,
                   char                 *structure);


/**
 *  @brief  Get a list of base pair probabiltities
 *
 *  @param  fc      The fold compound to operate on
 *  @param  cutoff  A probability threshold to filter the output
 *  @return         A list of base pair probabiltities
 */
vrna_ep_array_t
vrna_eqm_probs_bp(vrna_fold_compound_t *fc,
                  double               cutoff);


/* End base pair related functions */
/**@}*/


/* End thermodynamics */
/**@}*/

#endif
