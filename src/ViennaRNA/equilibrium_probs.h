#ifndef VIENNA_RNA_PACKAGE_EQUILIBRIUM_PROBS_H
#define VIENNA_RNA_PACKAGE_EQUILIBRIUM_PROBS_H

/* make this interface backward compatible with RNAlib < 2.2.0 */
#define VRNA_BACKWARD_COMPAT

#include <ViennaRNA/data_structures.h>
#include <ViennaRNA/params.h>

#ifdef VRNA_WARN_DEPRECATED
# ifdef __GNUC__
#  define DEPRECATED(func) func __attribute__ ((deprecated))
# else
#  define DEPRECATED(func) func
# endif
#else
# define DEPRECATED(func) func
#endif

/**
 *  @file     equilibrium_probs.h
 *  @ingroup  pf_fold
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

void  vrna_pairing_probs(vrna_fold_compound_t *vc, char *structure);

/**
 *  @brief Get the mean base pair distance in the thermodynamic ensemble from a probability matrix
 * 
 *  @f$<d> = \sum_{a,b} p_a p_b d(S_a,S_b)@f$\n
 *  this can be computed from the pair probs @f$p_ij@f$ as\n
 *  @f$<d> = \sum_{ij} p_{ij}(1-p_{ij})@f$
 * 
 *  @ingroup pf_fold
 *
 *  @param length The length of the sequence
 *  @param pr     The matrix containing the base pair probabilities
 *  @return       The mean pair distance of the structure ensemble
 */
double vrna_mean_bp_distance_pr(int length, FLT_OR_DBL *pr);

/**
 *  @brief Get the mean base pair distance in the thermodynamic ensemble
 * 
 *  @f$<d> = \sum_{a,b} p_a p_b d(S_a,S_b)@f$\n
 *  this can be computed from the pair probs @f$p_ij@f$ as\n
 *  @f$<d> = \sum_{ij} p_{ij}(1-p_{ij})@f$
 * 
 *  @ingroup pf_fold
 *
 *  @param vc     The fold compound data structure
 *  @return       The mean pair distance of the structure ensemble
 */
double vrna_mean_bp_distance(vrna_fold_compound_t *vc);

/**
 *  @brief  Compute stacking probabilities
 *
 *  For each possible base pair @f$(i,j)@f$, compute the probability of a stack
 *  @f$(i,j)@f$, @f$(i+1, j-1)@f$.
 *
 *  @ingroup pf_fold
 *
 *  @param  vc      The fold compound data structure with precomputed base pair probabilities
 *  @param  cutoff  A cutoff value that limits the output to stacks with @f$ p > \textrm{cutoff} @f$.
 *  @return         A list of stacks with enclosing base pair @f$(i,j)@f$ and probabiltiy @f$ p @f$
 */
vrna_ep_t *vrna_stack_prob(vrna_fold_compound_t *vc, double cutoff);

/**
 *  @brief Compute Boltzmann probabilities of dimerization without homodimers
 *
 *  Given the pair probabilities and free energies (in the null model) for a
 *  dimer AB and the two constituent monomers A and B, compute the conditional pair
 *  probabilities given that a dimer AB actually forms.
 *  Null model pair probabilities are given as a list as produced by
 *  vrna_plist_from_probs(), the dimer probabilities 'prAB' are modified in place.
 *
 *  @ingroup pf_cofold
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
void  vrna_pf_dimer_probs(double                  FAB,
                          double                  FA,
                          double                  FB,
                          vrna_ep_t               *prAB,
                          const vrna_ep_t         *prA,
                          const vrna_ep_t         *prB,
                          int                     Alength,
                          const vrna_exp_param_t  *exp_params);


#endif
