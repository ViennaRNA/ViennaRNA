#ifndef VIENNA_RNA_PACKAGE_BOLTZMANN_SAMPLING_H
#define VIENNA_RNA_PACKAGE_BOLTZMANN_SAMPLING_H

#include <ViennaRNA/datastructures/basic.h>

/**
 *  @file boltzmann_sampling.h
 *  @ingroup subopt_and_representatives
 *  @brief Boltzmann Sampling of secondary structures from the ensemble
 *
 *  A.k.a. Stochastic backtracking
 */


/**
 *  @addtogroup subopt_stochbt
 *  @{
 *  @brief  Functions to draw random structure samples from the ensemble according to their
 *          equilibrium probability
 */

/**
 *  @brief Sample a secondary structure of a subsequence from the Boltzmann ensemble according its probability
 *
 *  @pre    Unique multiloop decomposition has to be active upon creation of @p vc with vrna_fold_compound()
 *          or similar. This can be done easily by passing vrna_fold_compound() a model details parameter
 *          with vrna_md_t.uniq_ML = 1.
 *  @pre    vrna_pf() has to be called first to fill the partition function matrices
 *
 *  @param  vc      The fold compound data structure
 *  @param  length  The length of the subsequence to consider (starting with 5' end)
 *  @return         A sampled secondary structure in dot-bracket notation (or NULL on error)
 */
char *vrna_pbacktrack5(vrna_fold_compound_t *vc,
                       int                  length);
 
                      
/**
 *  @brief Samples multiple secondary structures non-redundantly from the Boltzmann ensemble according its probability
 *
 *  @ingroup subopt_stochbt
 *  @pre    The fold compound has to be obtained using the #VRNA_OPTION_HYBRID option in vrna_fold_compound()
 *  @pre    vrna_pf() has to be called first to fill the partition function matrices
 *
 *  @param  vc      	The fold compound data structure
 *  @param  num_samples The number of desired non-redundant samples
 *  @return         	A list of sampled secondary structures in dot-bracket notation
 */
char **vrna_non_redundant_pbacktrack(vrna_fold_compound_t *vc,
									 int				  num_samples);
									 
/**
 *  @brief Sample a secondary structure (consensus structure) from the Boltzmann ensemble according its probability
 *
 *  @pre    Unique multiloop decomposition has to be active upon creation of @p vc with vrna_fold_compound()
 *          or similar. This can be done easily by passing vrna_fold_compound() a model details parameter
 *          with vrna_md_t.uniq_ML = 1.
 *  @pre    vrna_pf() has to be called first to fill the partition function matrices
 *
 *  @note This function is polymorphic. It accepts #vrna_fold_compound_t of type
 *        #VRNA_FC_TYPE_SINGLE, and #VRNA_FC_TYPE_COMPARATIVE.
 *
 *  @note The function will automagically detect cicular RNAs based on the model_details in exp_params as
 *        provided via the #vrna_fold_compound_t
 *
 *  @param  vc      The fold compound data structure
 *  @return         A sampled secondary structure in dot-bracket notation (or NULL on error)
 */
char *vrna_pbacktrack(vrna_fold_compound_t *vc);


/**@}*/


#endif
