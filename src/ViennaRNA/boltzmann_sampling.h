#ifndef VIENNA_RNA_PACKAGE_BOLTZMANN_SAMPLING_H
#define VIENNA_RNA_PACKAGE_BOLTZMANN_SAMPLING_H

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
typedef void (vrna_boltzmann_sampling_callback)(const char  *stucture,
                                                void        *data);

typedef struct vrna_nr_memory_s *vrna_nr_memory_t;

#include <ViennaRNA/fold_compound.h>

/**
 *  @brief Sample a secondary structure of a subsequence from the Boltzmann ensemble according its probability
 *
 *  @pre    Unique multiloop decomposition has to be active upon creation of @p fc with vrna_fold_compound()
 *          or similar. This can be done easily by passing vrna_fold_compound() a model details parameter
 *          with vrna_md_t.uniq_ML = 1.
 *  @pre    vrna_pf() has to be called first to fill the partition function matrices
 *
 *  @note This function is polymorphic. It accepts #vrna_fold_compound_t of type
 *        #VRNA_FC_TYPE_SINGLE, and #VRNA_FC_TYPE_COMPARATIVE.
 *
 *  @see vrna_pbacktrack5_num(), vrna_pbacktrack5_num_cb(), vrna_pbacktrack(), vrna_pbacktrack_nr()
 *
 *  @param  fc      The fold compound data structure
 *  @param  length  The length of the subsequence to consider (starting with 5' end)
 *  @return         A sampled secondary structure in dot-bracket notation (or NULL on error)
 */
char *
vrna_pbacktrack5(vrna_fold_compound_t *fc,
                 unsigned int         length);


/**
 *  @brief Obtain a set of secondary structure samples for a subsequence from the Boltzmann ensemble according their probability
 *
 *  @pre    Unique multiloop decomposition has to be active upon creation of @p fc with vrna_fold_compound()
 *          or similar. This can be done easily by passing vrna_fold_compound() a model details parameter
 *          with vrna_md_t.uniq_ML = 1.
 *  @pre    vrna_pf() has to be called first to fill the partition function matrices
 *
 *  @note This function is polymorphic. It accepts #vrna_fold_compound_t of type
 *        #VRNA_FC_TYPE_SINGLE, and #VRNA_FC_TYPE_COMPARATIVE.
 *
 *  @see vrna_pbacktrack5(), vrna_pbacktrack5_num_cb(), vrna_pbacktrack_num(), vrna_pbacktrack_nr()
 *
 *  @param  fc            The fold compound data structure
 *  @param  length        The length of the subsequence to consider (starting with 5' end)
 *  @param  num_samples   The size of the sample set, i.e. number of structures
 *  @return               A set of secondary structure samples in dot-bracket notation (or NULL on error)
 */
char **
vrna_pbacktrack5_num(vrna_fold_compound_t *fc,
                     unsigned int         length,
                     unsigned int         num_samples);


unsigned int
vrna_pbacktrack5_num_cb(vrna_fold_compound_t              *fc,
                        unsigned int                      length,
                        unsigned int                      num_samples,
                        vrna_boltzmann_sampling_callback  *bs_cb,
                        void                              *data);


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


char **
vrna_pbacktrack_num(vrna_fold_compound_t  *fc,
                    unsigned int          num_samples);


unsigned int
vrna_pbacktrack_num_cb(vrna_fold_compound_t             *fc,
                       unsigned int                     num_samples,
                       vrna_boltzmann_sampling_callback *bs_cb,
                       void                             *data);


/**
 *  @brief Samples multiple secondary structures non-redundantly from the Boltzmann ensemble according its probability
 *
 *  @ingroup subopt_stochbt
 *  @pre    The fold compound has to be obtained using the #VRNA_OPTION_HYBRID option in vrna_fold_compound()
 *  @pre    vrna_pf() has to be called first to fill the partition function matrices
 *
 *  @note   In some cases, this function does not return the number of requested samples but a smaller number.
 *          This may happen if a) the number of requested structures is larger than the total number of structures
 *          in the ensemble, or b) numeric instabilities prevent the backtracking function to enumerate structures
 *          with very high free energies.
 *
 *  @param  vc        The fold compound data structure
 *  @param  num_samples The number of desired non-redundant samples
 *  @return           A list of sampled secondary structures in dot-bracket notation, terminated by @em NULL
 */
char **
vrna_pbacktrack_nr(vrna_fold_compound_t *vc,
                   unsigned int         num_samples);


unsigned int
vrna_pbacktrack_nr_cb(vrna_fold_compound_t              *vc,
                      unsigned int                      num_samples,
                      vrna_boltzmann_sampling_callback  *cb,
                      void                              *data);


char **
vrna_pbacktrack_nr_resume(vrna_fold_compound_t  *vc,
                          unsigned int          num_samples,
                          vrna_nr_memory_t      *nr_mem);


unsigned int
vrna_pbacktrack_nr_resume_cb(vrna_fold_compound_t             *vc,
                             unsigned int                     num_samples,
                             vrna_boltzmann_sampling_callback *bs_cb,
                             void                             *data,
                             vrna_nr_memory_t                 *nr_mem);


void
vrna_pbacktrack_nr_free(vrna_nr_memory_t s);


/**@}*/


#endif
