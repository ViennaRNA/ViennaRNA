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


/* End base pair related functions */
/**@}*/

/**
 *  @name Multimer probabilities computations
 *  @{
 */

/**
 *  @brief Compute Boltzmann probabilities of dimerization without homodimers
 *
 *  Given the pair probabilities and free energies (in the null model) for a
 *  dimer AB and the two constituent monomers A and B, compute the conditional pair
 *  probabilities given that a dimer AB actually forms.
 *  Null model pair probabilities are given as a list as produced by
 *  vrna_plist_from_probs(), the dimer probabilities 'prAB' are modified in place.
 *
 *  @param FAB        free energy of dimer AB
 *  @param FA         free energy of monomer A
 *  @param FB         free energy of monomer B
 *  @param prAB       pair probabilities for dimer
 *  @param prA        pair probabilities monomer
 *  @param prB        pair probabilities monomer
 *  @param Alength    Length of molecule A
 *  @param exp_params The precomputed Boltzmann factors
 */
void
vrna_pf_dimer_probs(double                  FAB,
                    double                  FA,
                    double                  FB,
                    vrna_ep_t               *prAB,
                    const vrna_ep_t         *prA,
                    const vrna_ep_t         *prB,
                    int                     Alength,
                    const vrna_exp_param_t  *exp_params);


/* End multimer probability related functions */
/**@}*/


/* End thermodynamics */
/**@}*/

#endif
