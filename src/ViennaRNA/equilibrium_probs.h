#ifndef VIENNA_RNA_PACKAGE_EQUILIBRIUM_PROBS_H
#define VIENNA_RNA_PACKAGE_EQUILIBRIUM_PROBS_H

#include <ViennaRNA/datastructures/basic.h>
#include <ViennaRNA/fold_compound.h>
#include <ViennaRNA/utils/structures.h>
#include <ViennaRNA/params/basic.h>

/**
 *  @file     equilibrium_probs.h
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
 *  @brief      Compute various thermodynamic properties using the partition function
 *
 *  Many thermodynamic properties can be derived from the partition function
 *  @f[ Q = \sum_{s \in \omega} e^{\frac{-E(s)}{kT}}. @f]
 *  In particular, for nucleic acids in equilibrium the probabilty @f$ p(F) @f$ of
 *  a particular structural feature @f$ F @f$ follows Boltzmanns law, i.e.
 *  @f[ p(F) \propto \sum_{s \mid F \in s} e^{\frac{-E(s)}{kT}}. @f]
 *  The actual probabilities can then be obtained from the ratio of those
 *  structures containing @f$ F @f$ and @em all structures, i.e.
 *  @f[ p(F) = \frac{1}{Q} \sum_{s \mid F \in s} e^{\frac{-E(s)}{kT}}. @f]
 *
 *  Consequently, a particular secondary structure @f$ s @f$ has equilibrium
 *  probability
 *  @f[ p(s) = \frac{1}{Q} e^{\frac{-E(s)}{kT}} @f]
 *  which can be easily computed once @f$ Q @f$ and @f$ E(s) @f$ are known.
 *
 *  On the other hand, efficient dynamic programming algorithms exist to
 *  compute the equilibrium probabilities
 *  @f[ p_{ij} = \frac{1}{Q} \sum_{s \mid (i,j) \in s} e^{\frac{-E(s)}{kT}} @f]
 *  of base pairs @f$ (i,j) @f$ without the need for exhaustive enumeration
 *  of @f$ s @f$.
 *
 *  This interface provides the functions for all thermodynamic property
 *  computations implemented in @em RNAlib.
 */


/**
 *  @name Base pair probabilities and derived computations
 *  @{
 */

int
vrna_pairing_probs(vrna_fold_compound_t *vc,
                   char                 *structure);


/**
 *  @brief Get the mean base pair distance in the thermodynamic ensemble from a probability matrix
 *
 *  @f[
 *  <d> = \sum_{a,b} p_a p_b d(S_a,S_b)
 *  @f]
 *  this can be computed from the pair probs @f$ p_{ij} @f$ as
 *  @f[
 *  <d> = \sum_{ij} p_{ij}(1-p_{ij})
 *  @f]
 *
 *  @param length The length of the sequence
 *  @param pr     The matrix containing the base pair probabilities
 *  @return       The mean pair distance of the structure ensemble
 */
double
vrna_mean_bp_distance_pr(int        length,
                         FLT_OR_DBL *pr);


/**
 *  @brief Get the mean base pair distance in the thermodynamic ensemble
 *
 *  @f[
 *  <d> = \sum_{a,b} p_a p_b d(S_a,S_b)
 *  @f]
 *  this can be computed from the pair probs @f$p_{ij}@f$ as
 *  @f[
 *  <d> = \sum_{ij} p_{ij}(1-p_{ij})
 *  @f]
 *
 *  @param vc     The fold compound data structure
 *  @return       The mean pair distance of the structure ensemble
 */
double
vrna_mean_bp_distance(vrna_fold_compound_t *vc);


/**
 *  @brief  Compute the Ensemble Defect for a given target structure provided as a @b vrna_ptable
 *
 *  Given a target structure @f$s@f$, compute the average dissimilarity of a randomly
 *  drawn structure from the ensemble, i.e.:
 *  @f[
 *    ED(s) = 1 - \frac{1}{n} \sum_{ij, (i,j) \in s} p_{ij} - \frac{1}{n} \sum_{i}(1 - s_i)q_i
 *  @f]
 *  with sequence length @f$n@f$, the probability @f$p_{ij}@f$ of a base pair @f$(i,j)@f$,
 *  the probability @f$q_i = 1 - \sum_j p_{ij}@f$ of nucleotide @f$i@f$ being unpaired, and
 *  the indicator variable @f$s_i = 1@f$ if @f$\exists (i,j) \in s@f$, and @f$s_i = 0@f$ otherwise.
 *
 *  @pre  The #vrna_fold_compound_t input parameter @p fc must contain a valid base pair
 *        probability matrix. This means that partition function and base pair probabilities
 *        must have been computed using @p fc before execution of this function!
 *
 *  @see vrna_pf(), vrna_pairing_probs(), vrna_ensemble_defect()
 *
 *  @param  fc          A fold_compound with pre-computed base pair probabilities
 *  @param  pt          A pair table representing a target structure
 *  @return             The ensemble defect with respect to the target structure, or -1. upon failure, e.g. pre-conditions are not met
 */
double
vrna_ensemble_defect_pt(vrna_fold_compound_t *fc,
                        const short          *pt);


/**
 *  @brief  Compute the Ensemble Defect for a given target structure
 *
 *  This is a wrapper around @b vrna_ensemble_defect_pt().
 *  Given a target structure @f$s@f$, compute the average dissimilarity of a randomly
 *  drawn structure from the ensemble, i.e.:
 *  @f[
 *    ED(s) = 1 - \frac{1}{n} \sum_{ij, (i,j) \in s} p_{ij} - \frac{1}{n} \sum_{i}(1 - s_i)q_i
 *  @f]
 *  with sequence length @f$n@f$, the probability @f$p_{ij}@f$ of a base pair @f$(i,j)@f$,
 *  the probability @f$q_i = 1 - \sum_j p_{ij}@f$ of nucleotide @f$i@f$ being unpaired, and
 *  the indicator variable @f$s_i = 1@f$ if @f$\exists (i,j) \in s@f$, and @f$s_i = 0@f$ otherwise.
 *
 *  @pre  The #vrna_fold_compound_t input parameter @p fc must contain a valid base pair
 *        probability matrix. This means that partition function and base pair probabilities
 *        must have been computed using @p fc before execution of this function!
 *
 *  @see vrna_pf(), vrna_pairing_probs(), vrna_ensemble_defect_pt()
 *
 *  @param  fc          A fold_compound with pre-computed base pair probabilities
 *  @param  structure   A target structure in dot-bracket notation
 *  @return             The ensemble defect with respect to the target structure, or -1. upon failure, e.g. pre-conditions are not met
 */
double
vrna_ensemble_defect(vrna_fold_compound_t *fc,
                     const char           *structure);


/**
 *  @brief  Compute a vector of positional entropies
 *
 *  This function computes the positional entropies from base pair probabilities
 *  as
 *  @f[
 *  S(i) = - \sum_j p_{ij} \log(p_{ij}) - q_i \log(q_i)
 *  @f]
 *  with unpaired probabilities @f$ q_i = 1 - \sum_j p_{ij} @f$.
 *
 *  Low entropy regions have little
 *  structural flexibility and the reliability of the predicted structure is
 *  high. High entropy implies many structural alternatives. While these
 *  alternatives may be functionally important, they make structure prediction
 *  more difficult and thus less reliable.
 *
 *  @pre  This function requires pre-computed base pair probabilities! Thus,
 *        vrna_pf() must be called beforehand.
 *
 *  @param  fc          A fold_compound with pre-computed base pair probabilities
 *  @return             A 1-based vector of positional entropies @f$ S(i) @f$. (position 0 contains the sequence length)
 */
double *
vrna_positional_entropy(vrna_fold_compound_t *fc);


/**
 *  @brief  Compute stacking probabilities
 *
 *  For each possible base pair @f$(i,j)@f$, compute the probability of a stack
 *  @f$(i,j)@f$, @f$(i+1, j-1)@f$.
 *
 *  @param  vc      The fold compound data structure with precomputed base pair probabilities
 *  @param  cutoff  A cutoff value that limits the output to stacks with @f$ p > \textrm{cutoff} @f$.
 *  @return         A list of stacks with enclosing base pair @f$(i,j)@f$ and probabiltiy @f$ p @f$
 */
vrna_ep_t *
vrna_stack_prob(vrna_fold_compound_t  *vc,
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

/**
 *  @name Structure probability computations
 *  @{
 */

/**
 *  @brief Compute the equilibrium probability of a particular secondary structure
 *
 *  The probability @f$p(s)@f$ of a particular secondary structure @f$s@f$ can
 *  be computed as
 *  @f[
 *    p(s) = \frac{exp(-\beta E(s)}{Z}
 *  @f]
 *  from the structures free energy @f$E(s)@f$ and the partition function
 *  @f[
 *    Z = \sum_s exp(-\beta E(s)),\quad\mathrm{with}\quad\beta = \frac{1}{RT}
 *  @f]
 *  where @f$R@f$ is the gas constant and @f$T@f$ the thermodynamic temperature.
 *
 *  @pre  The fold compound @p fc must have went through a call to vrna_pf() to
 *        fill the dynamic programming matrices with the corresponding partition
 *        function.
 *
 *  @param  fc          The fold compound data structure with precomputed partition function
 *  @param  structure   The secondary structure to compute the probability for in dot-bracket notation
 *  @return             The probability of the input structure (range @f$[0:1]@f$)
 */
double
vrna_pr_structure(vrna_fold_compound_t  *fc,
                  const char            *structure);


double
vrna_pr_energy(vrna_fold_compound_t *vc,
               double               e);


/* End structure probability related functions */
/**@}*/

/* End thermodynamics */
/**@}*/

#endif
