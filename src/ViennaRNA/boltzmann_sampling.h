#ifndef VIENNA_RNA_PACKAGE_BOLTZMANN_SAMPLING_H
#define VIENNA_RNA_PACKAGE_BOLTZMANN_SAMPLING_H

#include <ViennaRNA/data_structures.h>

/**
 *  @addtogroup pf_fold
 *  @{
 *    @file boltzmann_sampling.h
 *    @brief Boltzmann Sampling of secondary structures from the ensemble
 *
 *    A.k.a. Stochastic backtracking
 *  @}
 */

/**
 *  @brief Sample a secondary structure of a subsequence from the Boltzmann ensemble according its probability
 *
 *  @ingroup subopt_stochbt
 *  @pre    The fold compound has to be obtained using the #VRNA_OPTION_HYBRID option in vrna_fold_compound()
 *  @pre    vrna_pf() hasto be called first to fill the partition function matrices
 *
 *  @param  vc      The fold compound data structure
 *  @param  length  The length of the subsequence to consider (starting with 5' end)
 *  @return         A sampled secondary structure in dot-bracket notation
 */
char    *vrna_pbacktrack5(vrna_fold_compound_t *vc, int length);

/**
 *  @brief Sample a secondary structure from the Boltzmann ensemble according its probability
 *
 *  @ingroup subopt_stochbt
 *  @pre    The fold compound has to be obtained using the #VRNA_OPTION_HYBRID option in vrna_fold_compound()
 *  @pre    vrna_pf() hasto be called first to fill the partition function matrices
 *
 *  @note The function will automagically detect cicular RNAs based on the model_details in exp_params as
 *        provided via the #vrna_fold_compound_t
 *
 *  @param  vc      The fold compound data structure
 *  @param  length  The length of the subsequence to consider (starting with 5' end)
 *  @return         A sampled secondary structure in dot-bracket notation
 */
char    *vrna_pbacktrack(vrna_fold_compound_t *vc);

/**
 *  @addtogroup consensus_stochbt
 *  @{
 *
 *  @}
 */

/**
 *  @brief Sample a consensus secondary structure from the Boltzmann ensemble according its probability
 * 
 *  @ingroup consensus_stochbt
 *
 *  @see vrna_pf() for precomputing the partition function matrices, and
 *
 *  @param  vc    The #vrna_fold_compound_t of type #VRNA_VC_TYPE_ALIGNMENT with precomputed partition function matrices
 *  @param  prob  to be described (berni)
 *  @return       A sampled consensus secondary structure in dot-bracket notation
 */
char *vrna_pbacktrack_comparative(vrna_fold_compound_t *vc, double *prob);

#endif
