#ifndef VIENNA_RNA_PACKAGE_EQUILIBRIUM_PROBS_HELICES_H
#define VIENNA_RNA_PACKAGE_EQUILIBRIUM_PROBS_HELICES_H

#include <ViennaRNA/datastructures/basic.h>
#include <ViennaRNA/fold_compound.h>
#include <ViennaRNA/structures/problist.h>
#include <ViennaRNA/params/basic.h>

/**
 *  @file     ViennaRNA/probabilities/hrelices.h
 *  @ingroup  thermodynamics
 *  @brief    Equilibrium Helix Probability implementations
 *
 */

/*
 ################################################
 # HELIX PROBABILITY RELATED FUNCTIONS          #
 ################################################
 */


/**
 *  @addtogroup  thermodynamics
 *  @{
 */


/**
 *  @name Helix probabilities and derived computations
 *  @{
 */

vrna_array(vrna_ep_t)
vrna_eqm_probs_hx(vrna_fold_compound_t *fc,
                  unsigned int         min_size,
                  unsigned int         max_size,
                  double               cutoff);


/**
 *  @brief  Compute stacking probabilities
 *
 *  For each possible base pair @f$(i,j)@f$, compute the probability of a stack
 *  @f$(i,j)@f$, @f$(i+1, j-1)@f$.
 *
 *  @param  fc      The fold compound data structure with precomputed base pair probabilities
 *  @param  cutoff  A cutoff value that limits the output to stacks with @f$ p > \textrm{cutoff} @f$.
 *  @return         A list of stacks with enclosing base pair @f$(i,j)@f$ and probabiltiy @f$ p @f$
 */
vrna_ep_t *
vrna_stack_prob(vrna_fold_compound_t  *fc,
                double                cutoff);


/* End helix probability related functions */
/**@}*/


/* End thermodynamics */
/**@}*/

#endif
