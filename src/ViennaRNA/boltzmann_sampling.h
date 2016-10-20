#ifndef VIENNA_RNA_PACKAGE_BOLTZMANN_SAMPLING_H
#define VIENNA_RNA_PACKAGE_BOLTZMANN_SAMPLING_H

#include <ViennaRNA/data_structures.h>

/**
 *  @file boltzmann_sampling.h
 *  @ingroup subopt_and_representatives
 *  @brief Boltzmann Sampling of secondary structures from the ensemble
 *
 *  A.k.a. Stochastic backtracking
 */

/**
 *  @brief Sample a secondary structure of a subsequence from the Boltzmann ensemble according its probability
 *
 *  @ingroup subopt_stochbt
 *  @pre    The fold compound has to be obtained using the #VRNA_OPTION_HYBRID option in vrna_fold_compound()
 *  @pre    vrna_pf() has to be called first to fill the partition function matrices
 *
 *  @param  vc      The fold compound data structure
 *  @param  length  The length of the subsequence to consider (starting with 5' end)
 *  @return         A sampled secondary structure in dot-bracket notation
 */
char    *vrna_pbacktrack5(vrna_fold_compound_t *vc, int length);

/**
 *  @brief Sample a secondary structure (consensus structure) from the Boltzmann ensemble according its probability
 *
 *  @ingroup subopt_stochbt
 *  @pre    The dynamic programming (DP) matrices have to allow for unique multibranch loop decomposition, i.e.
 *          the vrna_md_t.uniq_ML flag has to be non-zero before calling vrna_fold_compound()
 *  @pre    vrna_pf() has to be called first to fill the partition function matrices
 *
 *  @note This function is polymorphic. It accepts #vrna_fold_compound_t of type
 *        #VRNA_FC_TYPE_SINGLE, and #VRNA_FC_TYPE_COMPARATIVE.
 *
 *  @note The function will automagically detect cicular RNAs based on the model_details in exp_params as
 *        provided via the #vrna_fold_compound_t
 *
 *  @param  vc      The fold compound data structure
 *  @return         A sampled secondary structure in dot-bracket notation
 */
char    *vrna_pbacktrack(vrna_fold_compound_t *vc);

#endif
