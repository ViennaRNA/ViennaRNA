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

/**
 *  @brief  Callback for Boltzmann sampling
 *
 * @callback
 * @parblock
 * This function will be called for each secondary structure that has been successfully backtraced
 * from the partition function DP matrices.
 * @endparblock
 *
 * @see vrna_pbacktrack5_num_cb(), vrna_pbacktrack_num_cb(), vrna_pbacktrack_nr_cb()
 *
 * @param structure The secondary structure in dot-bracket notation
 * @param data      Some arbitrary, auxiliary data address as provided to the calling function
 */
typedef void (vrna_boltzmann_sampling_callback)(const char  *stucture,
                                                void        *data);

/**
 *  @brief  Non-redundancy memory data structure
 *
 *  Initialize with @p NULL and pass its address to the corresponding
 *  functions vrna_pbacktrack_nr_resume(), etc.
 *
 *  @note Do not forget to release memory occupied by this data structure with vrna_pbacktrack_nr_free()
 *        after using it in any of the related functions!
 *
 *  @see vrna_pbacktrack_nr_resume(), vrna_pbacktrack_nr_resume_cb(), vrna_pbacktrack_nr_free()
 */
typedef struct vrna_nr_memory_s *vrna_nr_memory_t;

#include <ViennaRNA/fold_compound.h>

/**
 *  @brief Sample a secondary structure of a subsequence from the Boltzmann ensemble according its probability
 *
 *  Perform a probabilistic (stochastic) backtracing in the partition function DP arrays
 *  to obtain a secondary structure. The parameter @p length specifies the length
 *  of the substructure starting from the 5' end.
 *  The structure @f$ s @f$ is picked from the Boltzmann distributed ensemble according to its probability
 *  @f[
 *  p(s) = \frac{exp(-E(s) / kT)}{Z}
 *  @f]
 *  with partition function @f$ Z = \sum_s exp(-E(s) / kT) @f$, free energy @f$ E(s) @f$, Boltzmann constant
 *  @f$ k @f$ and thermodynamic temperature @f$ T @f$.
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
 *  Perform a probabilistic (stochastic) backtracing in the partition function DP arrays
 *  to obtain a set of secondary structures. The parameter @p length specifies the length
 *  of the substructure starting from the 5' end.
 *  Each structure @f$ s @f$ in the sample set is picked from the Boltzmann distributed ensemble
 *  according to its probability
 *  @f[
 *  p(s) = \frac{exp(-E(s) / kT)}{Z}
 *  @f]
 *  with partition function @f$ Z = \sum_s exp(-E(s) / kT) @f$, free energy @f$ E(s) @f$, Boltzmann constant
 *  @f$ k @f$ and thermodynamic temperature @f$ T @f$.
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


/**
 *  @brief Generate a set of secondary structure samples for a subsequence from the Boltzmann ensemble according their probability
 *
 *  Perform a probabilistic (stochastic) backtracing in the partition function DP arrays
 *  to obtain a set of secondary structures. The parameter @p length specifies the length
 *  of the substructure starting from the 5' end.
 *  Each structure @f$ s @f$ in the sample set is picked from the Boltzmann distributed ensemble
 *  according to its probability
 *  @f[
 *  p(s) = \frac{exp(-E(s) / kT)}{Z}
 *  @f]
 *  with partition function @f$ Z = \sum_s exp(-E(s) / kT) @f$, free energy @f$ E(s) @f$, Boltzmann constant
 *  @f$ k @f$ and thermodynamic temperature @f$ T @f$.
 *
 *  In contrast to vrna_pbacktrack5() and vrna_pbacktrack5_num() this function yields the
 *  structure samples through a callback mechanism.
 *
 *  @pre    Unique multiloop decomposition has to be active upon creation of @p fc with vrna_fold_compound()
 *          or similar. This can be done easily by passing vrna_fold_compound() a model details parameter
 *          with vrna_md_t.uniq_ML = 1.
 *  @pre    vrna_pf() has to be called first to fill the partition function matrices
 *
 *  @note This function is polymorphic. It accepts #vrna_fold_compound_t of type
 *        #VRNA_FC_TYPE_SINGLE, and #VRNA_FC_TYPE_COMPARATIVE.
 *
 *  @see vrna_pbacktrack5(), vrna_pbacktrack5_num(), vrna_pbacktrack_num_cb(), vrna_pbacktrack_nr_cb()
 *
 *  @param  fc            The fold compound data structure
 *  @param  length        The length of the subsequence to consider (starting with 5' end)
 *  @param  num_samples   The size of the sample set, i.e. number of structures
 *  @param  cb            The callback that receives the sampled structure
 *  @param  data          A data structure passed through to the callback @p bs_cb
 *  @return               The number of structures actually backtraced
 */
unsigned int
vrna_pbacktrack5_num_cb(vrna_fold_compound_t              *fc,
                        unsigned int                      length,
                        unsigned int                      num_samples,
                        vrna_boltzmann_sampling_callback  *cb,
                        void                              *data);


/**
 *  @brief Sample a secondary structure (consensus structure) from the Boltzmann ensemble according its probability
 *
 *  Perform a probabilistic (stochastic) backtracing in the partition function DP arrays
 *  to obtain a set of secondary structures. The structure @f$ s @f$ is picked from the
 *  Boltzmann distributed ensemble according to its probability
 *  @f[
 *  p(s) = \frac{exp(-E(s) / kT)}{Z}
 *  @f]
 *  with partition function @f$ Z = \sum_s exp(-E(s) / kT) @f$, free energy @f$ E(s) @f$, Boltzmann constant
 *  @f$ k @f$ and thermodynamic temperature @f$ T @f$.
 *
 *  @pre    Unique multiloop decomposition has to be active upon creation of @p fc with vrna_fold_compound()
 *          or similar. This can be done easily by passing vrna_fold_compound() a model details parameter
 *          with vrna_md_t.uniq_ML = 1.
 *  @pre    vrna_pf() has to be called first to fill the partition function matrices
 *
 *  @note This function is polymorphic. It accepts #vrna_fold_compound_t of type
 *        #VRNA_FC_TYPE_SINGLE, and #VRNA_FC_TYPE_COMPARATIVE.
 *
 *  @see vrna_pbacktrack_num(), vrna_pbacktrack_num_cb(), vrna_pbacktrack5(), vrna_pbacktrack_nr()
 *
 *  @param  fc      The fold compound data structure
 *  @return         A sampled secondary structure in dot-bracket notation (or NULL on error)
 */
char *
vrna_pbacktrack(vrna_fold_compound_t *fc);


/**
 *  @brief Obtain a set of secondary structure samples from the Boltzmann ensemble according their probability
 *
 *  Perform a probabilistic (stochastic) backtracing in the partition function DP arrays
 *  to obtain a set of secondary structures. Each structure @f$ s @f$ in the sample set is
 *  picked from the Boltzmann distributed ensemble according to its probability
 *  @f[
 *  p(s) = \frac{exp(-E(s) / kT)}{Z}
 *  @f]
 *  with partition function @f$ Z = \sum_s exp(-E(s) / kT) @f$, free energy @f$ E(s) @f$, Boltzmann constant
 *  @f$ k @f$ and thermodynamic temperature @f$ T @f$.
 *
 *  @pre    Unique multiloop decomposition has to be active upon creation of @p fc with vrna_fold_compound()
 *          or similar. This can be done easily by passing vrna_fold_compound() a model details parameter
 *          with vrna_md_t.uniq_ML = 1.
 *  @pre    vrna_pf() has to be called first to fill the partition function matrices
 *
 *  @note This function is polymorphic. It accepts #vrna_fold_compound_t of type
 *        #VRNA_FC_TYPE_SINGLE, and #VRNA_FC_TYPE_COMPARATIVE.
 *
 *  @see vrna_pbacktrack(), vrna_pbacktrack_num_cb(), vrna_pbacktrack5_num(), vrna_pbacktrack_nr()
 *
 *  @param  fc            The fold compound data structure
 *  @param  num_samples   The size of the sample set, i.e. number of structures
 *  @return               A set of secondary structure samples in dot-bracket notation (or NULL on error)
 */
char **
vrna_pbacktrack_num(vrna_fold_compound_t  *fc,
                    unsigned int          num_samples);


/**
 *  @brief Obtain a set of secondary structure samples from the Boltzmann ensemble according their probability
 *
 *  Perform a probabilistic (stochastic) backtracing in the partition function DP arrays
 *  to obtain a set of secondary structures. Each structure @f$ s @f$ in the sample set is
 *  picked from the Boltzmann distributed ensemble according to its probability
 *  @f[
 *  p(s) = \frac{exp(-E(s) / kT)}{Z}
 *  @f]
 *  with partition function @f$ Z = \sum_s exp(-E(s) / kT) @f$, free energy @f$ E(s) @f$, Boltzmann constant
 *  @f$ k @f$ and thermodynamic temperature @f$ T @f$.
 *
 *  In contrast to vrna_pbacktrack() and vrna_pbacktrack_num() this function yields the
 *  structure samples through a callback mechanism.
 *
 *  @pre    Unique multiloop decomposition has to be active upon creation of @p fc with vrna_fold_compound()
 *          or similar. This can be done easily by passing vrna_fold_compound() a model details parameter
 *          with vrna_md_t.uniq_ML = 1.
 *  @pre    vrna_pf() has to be called first to fill the partition function matrices
 *
 *  @note This function is polymorphic. It accepts #vrna_fold_compound_t of type
 *        #VRNA_FC_TYPE_SINGLE, and #VRNA_FC_TYPE_COMPARATIVE.
 *
 *  @see vrna_pbacktrack(), vrna_pbacktrack_num_cb(), vrna_pbacktrack5_num(), vrna_pbacktrack_nr()
 *
 *  @param  fc            The fold compound data structure
 *  @param  num_samples   The size of the sample set, i.e. number of structures
 *  @param  cb            The callback that receives the sampled structure
 *  @param  data          A data structure passed through to the callback @p bs_cb
 *  @return               The number of structures actually backtraced
 */
unsigned int
vrna_pbacktrack_num_cb(vrna_fold_compound_t             *fc,
                       unsigned int                     num_samples,
                       vrna_boltzmann_sampling_callback *cb,
                       void                             *data);


/**
 *  @brief Samples multiple secondary structures non-redundantly from the Boltzmann ensemble according its probability
 *
 *  @pre    Unique multiloop decomposition has to be active upon creation of @p fc with vrna_fold_compound()
 *          or similar. This can be done easily by passing vrna_fold_compound() a model details parameter
 *          with vrna_md_t.uniq_ML = 1.
 *  @pre    vrna_pf() has to be called first to fill the partition function matrices
 *
 *  @note This function is polymorphic. It accepts #vrna_fold_compound_t of type
 *        #VRNA_FC_TYPE_SINGLE, and #VRNA_FC_TYPE_COMPARATIVE.
 *
 *  @note   In some cases, this function does not return the number of requested samples but a smaller number.
 *          This may happen if a) the number of requested structures is larger than the total number of structures
 *          in the ensemble, b) numeric instabilities prevent the backtracking function to enumerate structures
 *          with very high free energies, or c) any other error occurs.
 *
 *  @see vrna_pbacktrack_nr_cb(), vrna_pbacktrack_nr_resume(), vrna_pbacktrack_num()
 *
 *  @param  fc          The fold compound data structure
 *  @param  num_samples The number of desired non-redundant samples
 *  @return             A list of sampled secondary structures in dot-bracket notation, terminated by @em NULL
 */
char **
vrna_pbacktrack_nr(vrna_fold_compound_t *fc,
                   unsigned int         num_samples);


/**
 *  @brief Samples multiple secondary structures non-redundantly from the Boltzmann ensemble according its probability
 *
 *  In contrast to vrna_pbacktrack_nr() this function yields the structure samples through a callback mechanism.
 *
 *  @pre    Unique multiloop decomposition has to be active upon creation of @p fc with vrna_fold_compound()
 *          or similar. This can be done easily by passing vrna_fold_compound() a model details parameter
 *          with vrna_md_t.uniq_ML = 1.
 *  @pre    vrna_pf() has to be called first to fill the partition function matrices
 *
 *  @note This function is polymorphic. It accepts #vrna_fold_compound_t of type
 *        #VRNA_FC_TYPE_SINGLE, and #VRNA_FC_TYPE_COMPARATIVE.
 *
 *  @note   In some cases, this function does not return the number of requested samples but a smaller number.
 *          This may happen if a) the number of requested structures is larger than the total number of structures
 *          in the ensemble, b) numeric instabilities prevent the backtracking function to enumerate structures
 *          with very high free energies, or c) any other error occurs.
 *
 *  @see vrna_pbacktrack_nr(), vrna_pbacktrack_nr_resume_cb(), vrna_pbacktrack_num_cb()
 *
 *  @param  fc            The fold compound data structure
 *  @param  num_samples   The number of desired non-redundant samples
 *  @param  cb            The callback that receives the sampled structure
 *  @param  data          A data structure passed through to the callback @p bs_cb
 *  @return               The number of structures actually backtraced
 */
unsigned int
vrna_pbacktrack_nr_cb(vrna_fold_compound_t              *fc,
                      unsigned int                      num_samples,
                      vrna_boltzmann_sampling_callback  *cb,
                      void                              *data);


/**
 *  @brief Samples multiple secondary structures non-redundantly from the Boltzmann ensemble according its probability
 *
 *  In contrast to vrna_pbacktrack_nr() where successive calls may yield redundant structures, this
 *  function provides an additional parameter that is used to store the non-redundancy memory. Consequently,
 *  this function may be called multiple times amd the total set of structures obtained is still
 *  non-redundant.
 *
 *  To achieve this, the user must provide a pointer to the memory data structure as last argument.
 *  This memory data structure must be initialized with @p NULL and will be re-set automatically to
 *  the corresponding non-redundancy memory data structure in the first call to vrna_pbacktrack_nr_resume().
 *  A successive sample call to this function may look like:
 * @code{.c}
 * vrna_nr_memory_t nonredundant_memory = NULL;
 *
 * // sample the first 100 structures
 * vrna_pbacktrack_nr_resume(fc, 100, &nonredundant_memory);
 *
 * // sample another 500 structures
 * vrna_pbacktrack_nr_resume(fc, 500, &nonredundant_memory);
 *
 * // release memory occupied by the non-redundant memory data structure
 * vrna_pbacktrack_nr_free(nonredundant_memory);
 * @endcode
 *  @pre    Unique multiloop decomposition has to be active upon creation of @p fc with vrna_fold_compound()
 *          or similar. This can be done easily by passing vrna_fold_compound() a model details parameter
 *          with vrna_md_t.uniq_ML = 1.
 *  @pre    vrna_pf() has to be called first to fill the partition function matrices
 *
 *  @note This function is polymorphic. It accepts #vrna_fold_compound_t of type
 *        #VRNA_FC_TYPE_SINGLE, and #VRNA_FC_TYPE_COMPARATIVE.
 *
 *  @note   In some cases, this function does not return the number of requested samples but a smaller number.
 *          This may happen if a) the number of requested structures is larger than the total number of structures
 *          in the ensemble, b) numeric instabilities prevent the backtracking function to enumerate structures
 *          with very high free energies, or c) any other error occurs.
 *
 *  @see    vrna_pbacktrack_nr_resume_cb(), vrna_pbacktrack_nr(), vrna_pbacktrack_num(),
 *          #vrna_nr_memory_t, vrna_pbacktrack_nr_free()
 *
 *  @param  fc            The fold compound data structure
 *  @param  num_samples   The number of desired non-redundant samples
 *  @param  nr_mem        A pointer to the non-redundancy memory data structure
 *  @return               A list of sampled secondary structures in dot-bracket notation, terminated by @em NULL
 */
char **
vrna_pbacktrack_nr_resume(vrna_fold_compound_t  *fc,
                          unsigned int          num_samples,
                          vrna_nr_memory_t      *nr_mem);


/**
 *  @brief Samples multiple secondary structures non-redundantly from the Boltzmann ensemble according its probability
 *
 *  Same as vrna_pbacktrack_nr_resume() but structure samples are obtained through a callback mechanism.
 *
 *  A successive sample call to this function may look like:
 * @code{.c}
 * vrna_nr_memory_t nonredundant_memory = NULL;
 *
 * // sample the first 100 structures
 * vrna_pbacktrack_nr_resume_cb(fc, 100, &callback_function, (void *)&callback_data, &nonredundant_memory);
 *
 * // sample another 500 structures
 * vrna_pbacktrack_nr_resume_cb(fc, 500, &callback_function, (void *)&callback_data, &nonredundant_memory);
 *
 * // release memory occupied by the non-redundant memory data structure
 * vrna_pbacktrack_nr_free(nonredundant_memory);
 * @endcode
 *  @pre    Unique multiloop decomposition has to be active upon creation of @p fc with vrna_fold_compound()
 *          or similar. This can be done easily by passing vrna_fold_compound() a model details parameter
 *          with vrna_md_t.uniq_ML = 1.
 *  @pre    vrna_pf() has to be called first to fill the partition function matrices
 *
 *  @note This function is polymorphic. It accepts #vrna_fold_compound_t of type
 *        #VRNA_FC_TYPE_SINGLE, and #VRNA_FC_TYPE_COMPARATIVE.
 *
 *  @note   In some cases, this function does not return the number of requested samples but a smaller number.
 *          This may happen if a) the number of requested structures is larger than the total number of structures
 *          in the ensemble, b) numeric instabilities prevent the backtracking function to enumerate structures
 *          with very high free energies, or c) any other error occurs.
 *
 *  @see    vrna_pbacktrack_nr_resume(), vrna_pbacktrack_nr_cb(), vrna_pbacktrack_num_cb(),
 *          #vrna_nr_memory_t, vrna_pbacktrack_nr_free()
 *
 *  @param  fc            The fold compound data structure
 *  @param  num_samples   The number of desired non-redundant samples
 *  @param  cb            The callback that receives the sampled structure
 *  @param  data          A data structure passed through to the callback @p bs_cb
 *  @param  nr_mem        A pointer to the non-redundancy memory data structure
 *  @return               The number of structures actually backtraced
 */
unsigned int
vrna_pbacktrack_nr_resume_cb(vrna_fold_compound_t             *fc,
                             unsigned int                     num_samples,
                             vrna_boltzmann_sampling_callback *cb,
                             void                             *data,
                             vrna_nr_memory_t                 *nr_mem);


/**
 *  @brief  Release memory occupied by a non-redundancy memory data structure
 *
 *  @see  #vrna_nr_memory_t, vrna_pbacktrack_nr_resume(), vrna_pbacktrack_nr_resume_cb()
 *
 *  @param  s   The non-redundancy memory data structure
 */
void
vrna_pbacktrack_nr_free(vrna_nr_memory_t s);


/**@}*/


#endif
