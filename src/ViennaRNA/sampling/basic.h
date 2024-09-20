#ifndef VIENNA_RNA_PACKAGE_BOLTZMANN_SAMPLING_H
#define VIENNA_RNA_PACKAGE_BOLTZMANN_SAMPLING_H

#ifdef VRNA_WARN_DEPRECATED
# if defined(DEPRECATED)
#   undef DEPRECATED
# endif
# if defined(__clang__)
#  define DEPRECATED(func, msg) func __attribute__ ((deprecated("", msg)))
# elif defined(__GNUC__)
#  define DEPRECATED(func, msg) func __attribute__ ((deprecated(msg)))
# else
#  define DEPRECATED(func, msg) func
# endif
#else
# define DEPRECATED(func, msg) func
#endif

/**
 *  @file ViennaRNA/sampling/basic.h
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
 *  @brief  Boltzmann sampling flag indicating default backtracing mode
 *
 *  @see    vrna_pbacktrack5_num(), vrna_pbacktrack5_cb(), vrna_pbacktrack5_resume(),
 *          vrna_pbacktrack5_resume_cb(), vrna_pbacktrack_num(), vrna_pbacktrack_cb(),
 *          vrna_pbacktrack_resume(), vrna_pbacktrack_resume_cb()
 */
#define VRNA_PBACKTRACK_DEFAULT         0

/**
 *  @brief  Boltzmann sampling flag indicating non-redundant backtracing mode
 *
 *  This flag will turn the Boltzmann sampling into non-redundant backtracing
 *  mode along the lines of @rstinline :cite:t:`michalik:2017` @endrst
 *
 *  @see    vrna_pbacktrack5_num(), vrna_pbacktrack5_cb(), vrna_pbacktrack5_resume(),
 *          vrna_pbacktrack5_resume_cb(), vrna_pbacktrack_num(), vrna_pbacktrack_cb(),
 *          vrna_pbacktrack_resume(), vrna_pbacktrack_resume_cb()
 */
#define VRNA_PBACKTRACK_NON_REDUNDANT   1

/**
 *  @brief  Callback for Boltzmann sampling
 *
 * @callback
 * @parblock
 * This function will be called for each secondary structure that has been successfully backtraced
 * from the partition function DP matrices.
 * @endparblock
 *
 * @see vrna_pbacktrack5_cb(), vrna_pbacktrack_cb(), vrna_pbacktrack5_resume_cb(),
 *      vrna_pbacktrack_resume_cb()
 *
 * @param structure The secondary structure in dot-bracket notation
 * @param data      Some arbitrary, auxiliary data address as provided to the calling function
 */
typedef void (*vrna_bs_result_f)(const char  *structure,
                                                void        *data);

DEPRECATED(typedef void (vrna_boltzmann_sampling_callback)(const char  *structure,
                                                void        *data),
           "Use vrna_bs_result_f instead!");


/**
 *  @brief  Boltzmann sampling memory data structure
 *
 *  This structure is required for properly resuming a previous sampling round in
 *  specialized Boltzmann sampling, such as non-redundant backtracking.
 *
 *  Initialize with @p NULL and pass its address to the corresponding
 *  functions vrna_pbacktrack5_resume(), etc.
 *
 *  @note Do not forget to release memory occupied by this data structure before
 *        losing its context! Use vrna_pbacktrack_mem_free().
 *
 *  @see  vrna_pbacktrack5_resume(), vrna_pbacktrack_resume(), vrna_pbacktrack5_resume_cb(),
 *        vrna_pbacktrack_resume_cb(), vrna_pbacktrack_mem_free()
 */
typedef struct vrna_pbacktrack_memory_s *vrna_pbacktrack_mem_t;

#include <ViennaRNA/fold_compound.h>

/**
 *  @brief Sample a secondary structure of a subsequence from the Boltzmann ensemble according its probability
 *
 *  Perform a probabilistic (stochastic) backtracing in the partition function DP arrays
 *  to obtain a secondary structure. The parameter @p length specifies the length
 *  of the substructure starting from the 5' end.
 *
 *  The structure @f$ s @f$ with free energy @f$ E(s) @f$ is picked from the Boltzmann distributed
 *  ensemble according to its probability
 *
 *  @f[
 *  p(s) = \frac{exp(-E(s) / kT)}{Z}
 *  @f]
 *
 *  with partition function @f$ Z = \sum_s exp(-E(s) / kT) @f$, Boltzmann constant @f$ k @f$ and
 *  thermodynamic temperature @f$ T @f$.
 *
 *  @pre    Unique multiloop decomposition has to be active upon creation of @p fc with vrna_fold_compound()
 *          or similar. This can be done easily by passing vrna_fold_compound() a model details parameter
 *          with vrna_md_t.uniq_ML = 1.<br>
 *          vrna_pf() has to be called first to fill the partition function matrices
 *
 *  @note This function is polymorphic. It accepts #vrna_fold_compound_t of type
 *        #VRNA_FC_TYPE_SINGLE, and #VRNA_FC_TYPE_COMPARATIVE.
 *
 *  @see vrna_pbacktrack5_num(), vrna_pbacktrack5_cb(), vrna_pbacktrack()
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
 *  to obtain a set of @p num_samples secondary structures. The parameter @p length specifies
 *  the length of the substructure starting from the 5' end.
 *
 *  Any structure @f$ s @f$ with free energy @f$ E(s) @f$ is picked from the Boltzmann distributed
 *  ensemble according to its probability
 *
 *  @f[
 *  p(s) = \frac{exp(-E(s) / kT)}{Z}
 *  @f]
 *
 *  with partition function @f$ Z = \sum_s exp(-E(s) / kT) @f$, Boltzmann constant @f$ k @f$ and
 *  thermodynamic temperature @f$ T @f$.
 *
 *  Using the @p options flag one can switch between regular (#VRNA_PBACKTRACK_DEFAULT) backtracing
 *  mode, and non-redundant sampling (#VRNA_PBACKTRACK_NON_REDUNDANT) along the lines of
 *  @rstinline :cite:t:`michalik:2017` @endrst.
 *
 *  @pre    Unique multiloop decomposition has to be active upon creation of @p fc with vrna_fold_compound()
 *          or similar. This can be done easily by passing vrna_fold_compound() a model details parameter
 *          with vrna_md_t.uniq_ML = 1.<br>
 *          vrna_pf() has to be called first to fill the partition function matrices
 *
 *  @note This function is polymorphic. It accepts #vrna_fold_compound_t of type
 *        #VRNA_FC_TYPE_SINGLE, and #VRNA_FC_TYPE_COMPARATIVE.
 *
 *  @warning  In non-redundant sampling mode (#VRNA_PBACKTRACK_NON_REDUNDANT), this function may
 *            not yield the full number of requested samples. This may happen if
 *            a)  the number of requested structures is larger than the total number
 *                of structuresin the ensemble,
 *            b)  numeric instabilities prevent the backtracking function to enumerate
 *                structures with high free energies, or
 *            c)  any other error occurs.
 *
 *  @see  vrna_pbacktrack5(), vrna_pbacktrack5_cb(), vrna_pbacktrack_num(),
 *        #VRNA_PBACKTRACK_DEFAULT, #VRNA_PBACKTRACK_NON_REDUNDANT
 *
 *  @param  fc            The fold compound data structure
 *  @param  num_samples   The size of the sample set, i.e. number of structures
 *  @param  length        The length of the subsequence to consider (starting with 5' end)
 *  @param  options       A bitwise OR-flag indicating the backtracing mode.
 *  @return               A set of secondary structure samples in dot-bracket notation terminated by NULL (or NULL on error)
 */
char **
vrna_pbacktrack5_num(vrna_fold_compound_t *fc,
                     unsigned int         num_samples,
                     unsigned int         length,
                     unsigned int         options);


/**
 *  @brief Obtain a set of secondary structure samples for a subsequence from the Boltzmann ensemble according their probability
 *
 *  Perform a probabilistic (stochastic) backtracing in the partition function DP arrays
 *  to obtain a set of @p num_samples secondary structures. The parameter @p length specifies
 *  the length of the substructure starting from the 5' end.
 *
 *  Any structure @f$ s @f$ with free energy @f$ E(s) @f$ is picked from the Boltzmann distributed
 *  ensemble according to its probability
 *
 *  @f[
 *  p(s) = \frac{exp(-E(s) / kT)}{Z}
 *  @f]
 *
 *  with partition function @f$ Z = \sum_s exp(-E(s) / kT) @f$, Boltzmann constant @f$ k @f$ and
 *  thermodynamic temperature @f$ T @f$.
 *
 *  Using the @p options flag one can switch between regular (#VRNA_PBACKTRACK_DEFAULT) backtracing
 *  mode, and non-redundant sampling (#VRNA_PBACKTRACK_NON_REDUNDANT) along the lines of
 *  @rstinline :cite:t:`michalik:2017` @endrst.
 *
 *  In contrast to vrna_pbacktrack5() and vrna_pbacktrack5_num() this function yields the
 *  structure samples through a callback mechanism.
 *
 *  @pre    Unique multiloop decomposition has to be active upon creation of @p fc with vrna_fold_compound()
 *          or similar. This can be done easily by passing vrna_fold_compound() a model details parameter
 *          with vrna_md_t.uniq_ML = 1.<br>
 *          vrna_pf() has to be called first to fill the partition function matrices
 *
 *  @note This function is polymorphic. It accepts #vrna_fold_compound_t of type
 *        #VRNA_FC_TYPE_SINGLE, and #VRNA_FC_TYPE_COMPARATIVE.
 *
 *  @warning  In non-redundant sampling mode (#VRNA_PBACKTRACK_NON_REDUNDANT), this function may
 *            not yield the full number of requested samples. This may happen if
 *            a)  the number of requested structures is larger than the total number
 *                of structuresin the ensemble,
 *            b)  numeric instabilities prevent the backtracking function to enumerate
 *                structures with high free energies, or
 *            c)  any other error occurs.
 *
 *  @see  vrna_pbacktrack5(), vrna_pbacktrack5_num(), vrna_pbacktrack_cb(),
 *        #VRNA_PBACKTRACK_DEFAULT, #VRNA_PBACKTRACK_NON_REDUNDANT
 *
 *  @param  fc            The fold compound data structure
 *  @param  num_samples   The size of the sample set, i.e. number of structures
 *  @param  length        The length of the subsequence to consider (starting with 5' end)
 *  @param  cb            The callback that receives the sampled structure
 *  @param  data          A data structure passed through to the callback @p cb
 *  @param  options       A bitwise OR-flag indicating the backtracing mode.
 *  @return               The number of structures actually backtraced
 */
unsigned int
vrna_pbacktrack5_cb(vrna_fold_compound_t              *fc,
                    unsigned int                      num_samples,
                    unsigned int                      length,
                    vrna_bs_result_f  cb,
                    void                              *data,
                    unsigned int                      options);


/**
 *  @brief Obtain a set of secondary structure samples for a subsequence from the Boltzmann ensemble according their probability
 *
 *  Perform a probabilistic (stochastic) backtracing in the partition function DP arrays
 *  to obtain a set of @p num_samples secondary structures. The parameter @p length specifies
 *  the length of the substructure starting from the 5' end.
 *
 *  Any structure @f$ s @f$ with free energy @f$ E(s) @f$ is picked from the Boltzmann distributed
 *  ensemble according to its probability
 *
 *  @f[
 *  p(s) = \frac{exp(-E(s) / kT)}{Z}
 *  @f]
 *
 *  with partition function @f$ Z = \sum_s exp(-E(s) / kT) @f$, Boltzmann constant @f$ k @f$ and
 *  thermodynamic temperature @f$ T @f$.
 *
 *  Using the @p options flag one can switch between regular (#VRNA_PBACKTRACK_DEFAULT) backtracing
 *  mode, and non-redundant sampling (#VRNA_PBACKTRACK_NON_REDUNDANT) along the lines of
 *  @rstinline :cite:t:`michalik:2017` @endrst.
 *
 *  In contrast to vrna_pbacktrack5_cb() this function allows for resuming a previous
 *  sampling round in specialized Boltzmann sampling, such as non-redundant backtracking.
 *  For that purpose, the user passes the address of a Boltzmann sampling data structure
 *  (#vrna_pbacktrack_mem_t) which will be re-used in each round of sampling, i.e. each
 *  successive call to vrna_pbacktrack5_resume_cb() or vrna_pbacktrack5_resume().
 *
 *  A successive sample call to this function may look like:
 * @code{.c}
 * vrna_pbacktrack_mem_t nonredundant_memory = NULL;
 *
 * // sample the first 100 structures
 * vrna_pbacktrack5_resume(fc,
 *                         100,
 *                         fc->length,
 *                         &nonredundant_memory,
 *                         options);
 *
 * // sample another 500 structures
 * vrna_pbacktrack5_resume(fc,
 *                         500,
 *                         fc->length,
 *                         &nonredundant_memory,
 *                         options);
 *
 * // release memory occupied by the non-redundant memory data structure
 * vrna_pbacktrack_mem_free(nonredundant_memory);
 * @endcode
 *
 *  @pre    Unique multiloop decomposition has to be active upon creation of @p fc with vrna_fold_compound()
 *          or similar. This can be done easily by passing vrna_fold_compound() a model details parameter
 *          with vrna_md_t.uniq_ML = 1.<br>
 *          vrna_pf() has to be called first to fill the partition function matrices
 *
 *  @note This function is polymorphic. It accepts #vrna_fold_compound_t of type
 *        #VRNA_FC_TYPE_SINGLE, and #VRNA_FC_TYPE_COMPARATIVE.
 *
 *  @warning  In non-redundant sampling mode (#VRNA_PBACKTRACK_NON_REDUNDANT), this function may
 *            not yield the full number of requested samples. This may happen if
 *            a)  the number of requested structures is larger than the total number
 *                of structuresin the ensemble,
 *            b)  numeric instabilities prevent the backtracking function to enumerate
 *                structures with high free energies, or
 *            c)  any other error occurs.
 *
 *  @see  vrna_pbacktrack5_resume_cb(), vrna_pbacktrack5_cb(), vrna_pbacktrack_resume(),
 *        #vrna_pbacktrack_mem_t, #VRNA_PBACKTRACK_DEFAULT, #VRNA_PBACKTRACK_NON_REDUNDANT,
 *        vrna_pbacktrack_mem_free
 *
 *  @param  fc            The fold compound data structure
 *  @param  num_samples   The size of the sample set, i.e. number of structures
 *  @param  length        The length of the subsequence to consider (starting with 5' end)
 *  @param  nr_mem        The address of the Boltzmann sampling memory data structure
 *  @param  options       A bitwise OR-flag indicating the backtracing mode.
 *  @return               A set of secondary structure samples in dot-bracket notation terminated by NULL (or NULL on error)
 */
char **
vrna_pbacktrack5_resume(vrna_fold_compound_t  *fc,
                        unsigned int          num_samples,
                        unsigned int          length,
                        vrna_pbacktrack_mem_t *nr_mem,
                        unsigned int          options);


/**
 *  @brief Obtain a set of secondary structure samples for a subsequence from the Boltzmann ensemble according their probability
 *
 *  Perform a probabilistic (stochastic) backtracing in the partition function DP arrays
 *  to obtain a set of @p num_samples secondary structures. The parameter @p length specifies
 *  the length of the substructure starting from the 5' end.
 *
 *  Any structure @f$ s @f$ with free energy @f$ E(s) @f$ is picked from the Boltzmann distributed
 *  ensemble according to its probability
 *
 *  @f[
 *  p(s) = \frac{exp(-E(s) / kT)}{Z}
 *  @f]
 *
 *  with partition function @f$ Z = \sum_s exp(-E(s) / kT) @f$, Boltzmann constant @f$ k @f$ and
 *  thermodynamic temperature @f$ T @f$.
 *
 *  Using the @p options flag one can switch between regular (#VRNA_PBACKTRACK_DEFAULT) backtracing
 *  mode, and non-redundant sampling (#VRNA_PBACKTRACK_NON_REDUNDANT) along the lines of
 *  @rstinline :cite:t:`michalik:2017` @endrst.
 *
 *  In contrast to vrna_pbacktrack5_resume() this function yields the structure samples
 *  through a callback mechanism.
 *
 *  A successive sample call to this function may look like:
 * @code{.c}
 * vrna_pbacktrack_mem_t nonredundant_memory = NULL;
 *
 * // sample the first 100 structures
 * vrna_pbacktrack5_resume_cb(fc,
 *                            100,
 *                            fc->length,
 *                            &callback_function,
 *                            (void *)&callback_data,
 *                            &nonredundant_memory,
 *                            options);
 *
 * // sample another 500 structures
 * vrna_pbacktrack5_resume_cb(fc,
 *                            500,
 *                            fc->length,
 *                            &callback_function,
 *                            (void *)&callback_data,
 *                            &nonredundant_memory,
 *                            options);
 *
 * // release memory occupied by the non-redundant memory data structure
 * vrna_pbacktrack_mem_free(nonredundant_memory);
 * @endcode
 *
 *  @pre    Unique multiloop decomposition has to be active upon creation of @p fc with vrna_fold_compound()
 *          or similar. This can be done easily by passing vrna_fold_compound() a model details parameter
 *          with vrna_md_t.uniq_ML = 1.<br>
 *          vrna_pf() has to be called first to fill the partition function matrices
 *
 *  @note This function is polymorphic. It accepts #vrna_fold_compound_t of type
 *        #VRNA_FC_TYPE_SINGLE, and #VRNA_FC_TYPE_COMPARATIVE.
 *
 *  @warning  In non-redundant sampling mode (#VRNA_PBACKTRACK_NON_REDUNDANT), this function may
 *            not yield the full number of requested samples. This may happen if
 *            a)  the number of requested structures is larger than the total number
 *                of structuresin the ensemble,
 *            b)  numeric instabilities prevent the backtracking function to enumerate
 *                structures with high free energies, or
 *            c)  any other error occurs.
 *
 *  @see  vrna_pbacktrack5_resume(), vrna_pbacktrack5_cb(), vrna_pbacktrack_resume_cb(),
 *        #vrna_pbacktrack_mem_t, #VRNA_PBACKTRACK_DEFAULT, #VRNA_PBACKTRACK_NON_REDUNDANT,
 *        vrna_pbacktrack_mem_free
 *
 *  @param  fc            The fold compound data structure
 *  @param  num_samples   The size of the sample set, i.e. number of structures
 *  @param  length        The length of the subsequence to consider (starting with 5' end)
 *  @param  cb            The callback that receives the sampled structure
 *  @param  data          A data structure passed through to the callback @p cb
 *  @param  nr_mem        The address of the Boltzmann sampling memory data structure
 *  @param  options       A bitwise OR-flag indicating the backtracing mode.
 *  @return               The number of structures actually backtraced
 */
unsigned int
vrna_pbacktrack5_resume_cb(vrna_fold_compound_t             *fc,
                           unsigned int                     num_samples,
                           unsigned int                     length,
                           vrna_bs_result_f cb,
                           void                             *data,
                           vrna_pbacktrack_mem_t            *nr_mem,
                           unsigned int                     options);


/**
 *  @brief Sample a secondary structure from the Boltzmann ensemble according its probability
 *
 *  Perform a probabilistic (stochastic) backtracing in the partition function DP arrays
 *  to obtain a secondary structure.
 *
 *  The structure @f$ s @f$ with free energy @f$ E(s) @f$ is picked from the Boltzmann distributed
 *  ensemble according to its probability
 *
 *  @f[
 *  p(s) = \frac{exp(-E(s) / kT)}{Z}
 *  @f]
 *
 *  with partition function @f$ Z = \sum_s exp(-E(s) / kT) @f$, Boltzmann constant @f$ k @f$ and
 *  thermodynamic temperature @f$ T @f$.
 *
 *  @pre    Unique multiloop decomposition has to be active upon creation of @p fc with vrna_fold_compound()
 *          or similar. This can be done easily by passing vrna_fold_compound() a model details parameter
 *          with vrna_md_t.uniq_ML = 1.<br>
 *          vrna_pf() has to be called first to fill the partition function matrices
 *
 *  @note This function is polymorphic. It accepts #vrna_fold_compound_t of type
 *        #VRNA_FC_TYPE_SINGLE, and #VRNA_FC_TYPE_COMPARATIVE.
 *
 *  @see vrna_pbacktrack5(), vrna_pbacktrack_num, vrna_pbacktrack_cb()
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
 *  to obtain a set of @p num_samples secondary structures.
 *
 *  Any structure @f$ s @f$ with free energy @f$ E(s) @f$ is picked from the Boltzmann distributed
 *  ensemble according to its probability
 *
 *  @f[
 *  p(s) = \frac{exp(-E(s) / kT)}{Z}
 *  @f]
 *
 *  with partition function @f$ Z = \sum_s exp(-E(s) / kT) @f$, Boltzmann constant @f$ k @f$ and
 *  thermodynamic temperature @f$ T @f$.
 *
 *  Using the @p options flag one can switch between regular (#VRNA_PBACKTRACK_DEFAULT) backtracing
 *  mode, and non-redundant sampling (#VRNA_PBACKTRACK_NON_REDUNDANT) along the lines of
 *  @rstinline :cite:t:`michalik:2017` @endrst.
 *
 *  @pre    Unique multiloop decomposition has to be active upon creation of @p fc with vrna_fold_compound()
 *          or similar. This can be done easily by passing vrna_fold_compound() a model details parameter
 *          with vrna_md_t.uniq_ML = 1.<br>
 *          vrna_pf() has to be called first to fill the partition function matrices
 *
 *  @note This function is polymorphic. It accepts #vrna_fold_compound_t of type
 *        #VRNA_FC_TYPE_SINGLE, and #VRNA_FC_TYPE_COMPARATIVE.
 *
 *  @warning  In non-redundant sampling mode (#VRNA_PBACKTRACK_NON_REDUNDANT), this function may
 *            not yield the full number of requested samples. This may happen if
 *            a)  the number of requested structures is larger than the total number
 *                of structuresin the ensemble,
 *            b)  numeric instabilities prevent the backtracking function to enumerate
 *                structures with high free energies, or
 *            c)  any other error occurs.
 *
 *  @see  vrna_pbacktrack(), vrna_pbacktrack_cb(), vrna_pbacktrack5_num(),
 *        #VRNA_PBACKTRACK_DEFAULT, #VRNA_PBACKTRACK_NON_REDUNDANT
 *
 *  @param  fc            The fold compound data structure
 *  @param  num_samples   The size of the sample set, i.e. number of structures
 *  @param  options       A bitwise OR-flag indicating the backtracing mode.
 *  @return               A set of secondary structure samples in dot-bracket notation terminated by NULL (or NULL on error)
 */
char **
vrna_pbacktrack_num(vrna_fold_compound_t  *fc,
                    unsigned int          num_samples,
                    unsigned int          options);


/**
 *  @brief Obtain a set of secondary structure samples from the Boltzmann ensemble according their probability
 *
 *  Perform a probabilistic (stochastic) backtracing in the partition function DP arrays
 *  to obtain a set of @p num_samples secondary structures.
 *
 *  Any structure @f$ s @f$ with free energy @f$ E(s) @f$ is picked from the Boltzmann distributed
 *  ensemble according to its probability
 *
 *  @f[
 *  p(s) = \frac{exp(-E(s) / kT)}{Z}
 *  @f]
 *
 *  with partition function @f$ Z = \sum_s exp(-E(s) / kT) @f$, Boltzmann constant @f$ k @f$ and
 *  thermodynamic temperature @f$ T @f$.
 *
 *  Using the @p options flag one can switch between regular (#VRNA_PBACKTRACK_DEFAULT) backtracing
 *  mode, and non-redundant sampling (#VRNA_PBACKTRACK_NON_REDUNDANT) along the lines of
 *  @rstinline :cite:t:`michalik:2017` @endrst.
 *
 *  In contrast to vrna_pbacktrack() and vrna_pbacktrack_num() this function yields the
 *  structure samples through a callback mechanism.
 *
 *  @pre    Unique multiloop decomposition has to be active upon creation of @p fc with vrna_fold_compound()
 *          or similar. This can be done easily by passing vrna_fold_compound() a model details parameter
 *          with vrna_md_t.uniq_ML = 1.<br>
 *          vrna_pf() has to be called first to fill the partition function matrices
 *
 *  @note This function is polymorphic. It accepts #vrna_fold_compound_t of type
 *        #VRNA_FC_TYPE_SINGLE, and #VRNA_FC_TYPE_COMPARATIVE.
 *
 *  @warning  In non-redundant sampling mode (#VRNA_PBACKTRACK_NON_REDUNDANT), this function may
 *            not yield the full number of requested samples. This may happen if
 *            a)  the number of requested structures is larger than the total number
 *                of structuresin the ensemble,
 *            b)  numeric instabilities prevent the backtracking function to enumerate
 *                structures with high free energies, or
 *            c)  any other error occurs.
 *
 *  @see  vrna_pbacktrack(), vrna_pbacktrack_num(), vrna_pbacktrack5_cb(),
 *        #VRNA_PBACKTRACK_DEFAULT, #VRNA_PBACKTRACK_NON_REDUNDANT
 *
 *  @param  fc            The fold compound data structure
 *  @param  num_samples   The size of the sample set, i.e. number of structures
 *  @param  cb            The callback that receives the sampled structure
 *  @param  data          A data structure passed through to the callback @p cb
 *  @param  options       A bitwise OR-flag indicating the backtracing mode.
 *  @return               The number of structures actually backtraced
 */
unsigned int
vrna_pbacktrack_cb(vrna_fold_compound_t             *fc,
                   unsigned int                     num_samples,
                   vrna_bs_result_f cb,
                   void                             *data,
                   unsigned int                     options);


/**
 *  @brief Obtain a set of secondary structure samples from the Boltzmann ensemble according their probability
 *
 *  Perform a probabilistic (stochastic) backtracing in the partition function DP arrays
 *  to obtain a set of @p num_samples secondary structures.
 *
 *  Any structure @f$ s @f$ with free energy @f$ E(s) @f$ is picked from the Boltzmann distributed
 *  ensemble according to its probability
 *
 *  @f[
 *  p(s) = \frac{exp(-E(s) / kT)}{Z}
 *  @f]
 *
 *  with partition function @f$ Z = \sum_s exp(-E(s) / kT) @f$, Boltzmann constant @f$ k @f$ and
 *  thermodynamic temperature @f$ T @f$.
 *
 *  Using the @p options flag one can switch between regular (#VRNA_PBACKTRACK_DEFAULT) backtracing
 *  mode, and non-redundant sampling (#VRNA_PBACKTRACK_NON_REDUNDANT) along the lines of
 *  @rstinline :cite:t:`michalik:2017` @endrst.
 *
 *  In contrast to vrna_pbacktrack_cb() this function allows for resuming a previous
 *  sampling round in specialized Boltzmann sampling, such as non-redundant backtracking.
 *  For that purpose, the user passes the address of a Boltzmann sampling data structure
 *  (#vrna_pbacktrack_mem_t) which will be re-used in each round of sampling, i.e. each
 *  successive call to vrna_pbacktrack_resume_cb() or vrna_pbacktrack_resume().
 *
 *  A successive sample call to this function may look like:
 * @code{.c}
 * vrna_pbacktrack_mem_t nonredundant_memory = NULL;
 *
 * // sample the first 100 structures
 * vrna_pbacktrack_resume(fc,
 *                        100,
 *                        &nonredundant_memory,
 *                        options);
 *
 * // sample another 500 structures
 * vrna_pbacktrack_resume(fc,
 *                        500,
 *                        &nonredundant_memory,
 *                        options);
 *
 * // release memory occupied by the non-redundant memory data structure
 * vrna_pbacktrack_mem_free(nonredundant_memory);
 * @endcode
 *
 *  @pre    Unique multiloop decomposition has to be active upon creation of @p fc with vrna_fold_compound()
 *          or similar. This can be done easily by passing vrna_fold_compound() a model details parameter
 *          with vrna_md_t.uniq_ML = 1.<br>
 *          vrna_pf() has to be called first to fill the partition function matrices
 *
 *  @note This function is polymorphic. It accepts #vrna_fold_compound_t of type
 *        #VRNA_FC_TYPE_SINGLE, and #VRNA_FC_TYPE_COMPARATIVE.
 *
 *  @warning  In non-redundant sampling mode (#VRNA_PBACKTRACK_NON_REDUNDANT), this function may
 *            not yield the full number of requested samples. This may happen if
 *            a)  the number of requested structures is larger than the total number
 *                of structuresin the ensemble,
 *            b)  numeric instabilities prevent the backtracking function to enumerate
 *                structures with high free energies, or
 *            c)  any other error occurs.
 *
 *  @see  vrna_pbacktrack_resume_cb(), vrna_pbacktrack_cb(), vrna_pbacktrack5_resume(),
 *        #vrna_pbacktrack_mem_t, #VRNA_PBACKTRACK_DEFAULT, #VRNA_PBACKTRACK_NON_REDUNDANT,
 *        vrna_pbacktrack_mem_free
 *
 *  @param  fc            The fold compound data structure
 *  @param  num_samples   The size of the sample set, i.e. number of structures
 *  @param  nr_mem        The address of the Boltzmann sampling memory data structure
 *  @param  options       A bitwise OR-flag indicating the backtracing mode.
 *  @return               A set of secondary structure samples in dot-bracket notation terminated by NULL (or NULL on error)
 */
char **
vrna_pbacktrack_resume(vrna_fold_compound_t   *fc,
                       unsigned int           num_samples,
                       vrna_pbacktrack_mem_t  *nr_mem,
                       unsigned int           options);


/**
 *  @brief Obtain a set of secondary structure samples from the Boltzmann ensemble according their probability
 *
 *  Perform a probabilistic (stochastic) backtracing in the partition function DP arrays
 *  to obtain a set of @p num_samples secondary structures.
 *
 *  Any structure @f$ s @f$ with free energy @f$ E(s) @f$ is picked from the Boltzmann distributed
 *  ensemble according to its probability
 *
 *  @f[
 *  p(s) = \frac{exp(-E(s) / kT)}{Z}
 *  @f]
 *
 *  with partition function @f$ Z = \sum_s exp(-E(s) / kT) @f$, Boltzmann constant @f$ k @f$ and
 *  thermodynamic temperature @f$ T @f$.
 *
 *  Using the @p options flag one can switch between regular (#VRNA_PBACKTRACK_DEFAULT) backtracing
 *  mode, and non-redundant sampling (#VRNA_PBACKTRACK_NON_REDUNDANT) along the lines of
 *  @rstinline :cite:t:`michalik:2017` @endrst.
 *
 *  In contrast to vrna_pbacktrack5_resume() this function yields the structure samples
 *  through a callback mechanism.
 *
 *  A successive sample call to this function may look like:
 * @code{.c}
 * vrna_pbacktrack_mem_t nonredundant_memory = NULL;
 *
 * // sample the first 100 structures
 * vrna_pbacktrack5_resume_cb(fc,
 *                            100,
 *                            &callback_function,
 *                            (void *)&callback_data,
 *                            &nonredundant_memory,
 *                            options);
 *
 * // sample another 500 structures
 * vrna_pbacktrack5_resume_cb(fc,
 *                            500,
 *                            &callback_function,
 *                            (void *)&callback_data,
 *                            &nonredundant_memory,
 *                            options);
 *
 * // release memory occupied by the non-redundant memory data structure
 * vrna_pbacktrack_mem_free(nonredundant_memory);
 * @endcode
 *
 *  @pre    Unique multiloop decomposition has to be active upon creation of @p fc with vrna_fold_compound()
 *          or similar. This can be done easily by passing vrna_fold_compound() a model details parameter
 *          with vrna_md_t.uniq_ML = 1.<br>
 *          vrna_pf() has to be called first to fill the partition function matrices
 *
 *  @note This function is polymorphic. It accepts #vrna_fold_compound_t of type
 *        #VRNA_FC_TYPE_SINGLE, and #VRNA_FC_TYPE_COMPARATIVE.
 *
 *  @warning  In non-redundant sampling mode (#VRNA_PBACKTRACK_NON_REDUNDANT), this function may
 *            not yield the full number of requested samples. This may happen if
 *            a)  the number of requested structures is larger than the total number
 *                of structuresin the ensemble,
 *            b)  numeric instabilities prevent the backtracking function to enumerate
 *                structures with high free energies, or
 *            c)  any other error occurs.
 *
 *  @see  vrna_pbacktrack_resume(), vrna_pbacktrack_cb(), vrna_pbacktrack5_resume_cb(),
 *        #vrna_pbacktrack_mem_t, #VRNA_PBACKTRACK_DEFAULT, #VRNA_PBACKTRACK_NON_REDUNDANT,
 *        vrna_pbacktrack_mem_free
 *
 *  @param  fc            The fold compound data structure
 *  @param  num_samples   The size of the sample set, i.e. number of structures
 *  @param  cb            The callback that receives the sampled structure
 *  @param  data          A data structure passed through to the callback @p cb
 *  @param  nr_mem        The address of the Boltzmann sampling memory data structure
 *  @param  options       A bitwise OR-flag indicating the backtracing mode.
 *  @return               The number of structures actually backtraced
 */
unsigned int
vrna_pbacktrack_resume_cb(vrna_fold_compound_t              *fc,
                          unsigned int                      num_samples,
                          vrna_bs_result_f  cb,
                          void                              *data,
                          vrna_pbacktrack_mem_t             *nr_mem,
                          unsigned int                      options);








/**
 *  @brief Sample a secondary structure of a subsequence from the Boltzmann ensemble according its probability
 *
 *  Perform a probabilistic (stochastic) backtracing in the partition function DP arrays
 *  to obtain a secondary structure. The parameters @p start and @p end specify the interval
 *  @f$ [start:end] @f$ of the subsequence with @f$ 1 \leq start < end \leq n@f$ for sequence length @f$ n @f$,
 *  the structure @f$ s_{start,end} @f$ should be drawn from.
 *
 *  The resulting substructure @f$ s_{start,end} @f$ with free energy @f$ E(s_{start, end}) @f$ is picked from
 *  the Boltzmann distributed sub ensemble of all structures within the interval @f$ [start:end] @f$ according
 *  to its probability
 *
 *  @f[
 *  p(s_{start,end}) = \frac{exp(-E(s_{start,end}) / kT)}{Z_{start,end}}
 *  @f]
 *
 *  with partition function @f$ Z_{start,end} = \sum_{s_{start,end}} exp(-E(s_{start,end}) / kT) @f$,
 *  Boltzmann constant @f$ k @f$ and thermodynamic temperature @f$ T @f$.
 *
 *  @pre    Unique multiloop decomposition has to be active upon creation of @p fc with vrna_fold_compound()
 *          or similar. This can be done easily by passing vrna_fold_compound() a model details parameter
 *          with vrna_md_t.uniq_ML = 1.<br>
 *          vrna_pf() has to be called first to fill the partition function matrices
 *
 *  @note This function is polymorphic. It accepts #vrna_fold_compound_t of type
 *        #VRNA_FC_TYPE_SINGLE, and #VRNA_FC_TYPE_COMPARATIVE.
 *
 *  @see vrna_pbacktrack_sub_num(), vrna_pbacktrack_sub_cb(), vrna_pbacktrack()
 *
 *  @param  fc      The fold compound data structure
 *  @param  start   The start of  the subsequence to consider, i.e. 5'-end position(1-based)
 *  @param  end     The end of the subsequence to consider, i.e. 3'-end position (1-based)
 *  @return         A sampled secondary structure in dot-bracket notation (or NULL on error)
 */
char *
vrna_pbacktrack_sub(vrna_fold_compound_t *fc,
                    unsigned int         start,
                    unsigned int         end);


/**
 *  @brief Obtain a set of secondary structure samples for a subsequence from the Boltzmann ensemble according their probability
 *
 *  Perform a probabilistic (stochastic) backtracing in the partition function DP arrays
 *  to obtain a set of @p num_samples secondary structures. The parameter @p length specifies
 *  the length of the substructure starting from the 5' end.
 *
 *  Any structure @f$ s @f$ with free energy @f$ E(s) @f$ is picked from the Boltzmann distributed
 *  ensemble according to its probability
 *
 *  @f[
 *  p(s) = \frac{exp(-E(s) / kT)}{Z}
 *  @f]
 *
 *  with partition function @f$ Z = \sum_s exp(-E(s) / kT) @f$, Boltzmann constant @f$ k @f$ and
 *  thermodynamic temperature @f$ T @f$.
 *
 *  Using the @p options flag one can switch between regular (#VRNA_PBACKTRACK_DEFAULT) backtracing
 *  mode, and non-redundant sampling (#VRNA_PBACKTRACK_NON_REDUNDANT) along the lines of
 *  @rstinline :cite:t:`michalik:2017` @endrst.
 *
 *  @pre    Unique multiloop decomposition has to be active upon creation of @p fc with vrna_fold_compound()
 *          or similar. This can be done easily by passing vrna_fold_compound() a model details parameter
 *          with vrna_md_t.uniq_ML = 1.<br>
 *          vrna_pf() has to be called first to fill the partition function matrices
 *
 *  @note This function is polymorphic. It accepts #vrna_fold_compound_t of type
 *        #VRNA_FC_TYPE_SINGLE, and #VRNA_FC_TYPE_COMPARATIVE.
 *
 *  @warning  In non-redundant sampling mode (#VRNA_PBACKTRACK_NON_REDUNDANT), this function may
 *            not yield the full number of requested samples. This may happen if
 *            a)  the number of requested structures is larger than the total number
 *                of structuresin the ensemble,
 *            b)  numeric instabilities prevent the backtracking function to enumerate
 *                structures with high free energies, or
 *            c)  any other error occurs.
 *
 *  @see  vrna_pbacktrack_sub(), vrna_pbacktrack_sub_cb(), vrna_pbacktrack_num(),
 *        #VRNA_PBACKTRACK_DEFAULT, #VRNA_PBACKTRACK_NON_REDUNDANT
 *
 *  @param  fc            The fold compound data structure
 *  @param  num_samples   The size of the sample set, i.e. number of structures
 *  @param  start         The start of  the subsequence to consider, i.e. 5'-end position(1-based)
 *  @param  end           The end of the subsequence to consider, i.e. 3'-end position (1-based)
 *  @param  options       A bitwise OR-flag indicating the backtracing mode.
 *  @return               A set of secondary structure samples in dot-bracket notation terminated by NULL (or NULL on error)
 */
char **
vrna_pbacktrack_sub_num(vrna_fold_compound_t *fc,
                        unsigned int         num_samples,
                        unsigned int         start,
                        unsigned int         end,
                        unsigned int         options);


/**
 *  @brief Obtain a set of secondary structure samples for a subsequence from the Boltzmann ensemble according their probability
 *
 *  Perform a probabilistic (stochastic) backtracing in the partition function DP arrays
 *  to obtain a set of @p num_samples secondary structures. The parameter @p length specifies
 *  the length of the substructure starting from the 5' end.
 *
 *  Any structure @f$ s @f$ with free energy @f$ E(s) @f$ is picked from the Boltzmann distributed
 *  ensemble according to its probability
 *
 *  @f[
 *  p(s) = \frac{exp(-E(s) / kT)}{Z}
 *  @f]
 *
 *  with partition function @f$ Z = \sum_s exp(-E(s) / kT) @f$, Boltzmann constant @f$ k @f$ and
 *  thermodynamic temperature @f$ T @f$.
 *
 *  Using the @p options flag one can switch between regular (#VRNA_PBACKTRACK_DEFAULT) backtracing
 *  mode, and non-redundant sampling (#VRNA_PBACKTRACK_NON_REDUNDANT) along the lines of
 *  @rstinline :cite:t:`michalik:2017` @endrst.
 *
 *  In contrast to vrna_pbacktrack5() and vrna_pbacktrack5_num() this function yields the
 *  structure samples through a callback mechanism.
 *
 *  @pre    Unique multiloop decomposition has to be active upon creation of @p fc with vrna_fold_compound()
 *          or similar. This can be done easily by passing vrna_fold_compound() a model details parameter
 *          with vrna_md_t.uniq_ML = 1.<br>
 *          vrna_pf() has to be called first to fill the partition function matrices
 *
 *  @note This function is polymorphic. It accepts #vrna_fold_compound_t of type
 *        #VRNA_FC_TYPE_SINGLE, and #VRNA_FC_TYPE_COMPARATIVE.
 *
 *  @warning  In non-redundant sampling mode (#VRNA_PBACKTRACK_NON_REDUNDANT), this function may
 *            not yield the full number of requested samples. This may happen if
 *            a)  the number of requested structures is larger than the total number
 *                of structuresin the ensemble,
 *            b)  numeric instabilities prevent the backtracking function to enumerate
 *                structures with high free energies, or
 *            c)  any other error occurs.
 *
 *  @see  vrna_pbacktrack5(), vrna_pbacktrack5_num(), vrna_pbacktrack_cb(),
 *        #VRNA_PBACKTRACK_DEFAULT, #VRNA_PBACKTRACK_NON_REDUNDANT
 *
 *  @param  fc            The fold compound data structure
 *  @param  num_samples   The size of the sample set, i.e. number of structures
 *  @param  start         The start of  the subsequence to consider, i.e. 5'-end position(1-based)
 *  @param  end           The end of the subsequence to consider, i.e. 3'-end position (1-based)
 *  @param  cb            The callback that receives the sampled structure
 *  @param  data          A data structure passed through to the callback @p cb
 *  @param  options       A bitwise OR-flag indicating the backtracing mode.
 *  @return               The number of structures actually backtraced
 */
unsigned int
vrna_pbacktrack_sub_cb(vrna_fold_compound_t              *fc,
                       unsigned int                      num_samples,
                       unsigned int                      start,
                       unsigned int                      end,
                       vrna_bs_result_f  cb,
                       void                              *data,
                       unsigned int                      options);


/**
 *  @brief Obtain a set of secondary structure samples for a subsequence from the Boltzmann ensemble according their probability
 *
 *  Perform a probabilistic (stochastic) backtracing in the partition function DP arrays
 *  to obtain a set of @p num_samples secondary structures. The parameter @p length specifies
 *  the length of the substructure starting from the 5' end.
 *
 *  Any structure @f$ s @f$ with free energy @f$ E(s) @f$ is picked from the Boltzmann distributed
 *  ensemble according to its probability
 *
 *  @f[
 *  p(s) = \frac{exp(-E(s) / kT)}{Z}
 *  @f]
 *
 *  with partition function @f$ Z = \sum_s exp(-E(s) / kT) @f$, Boltzmann constant @f$ k @f$ and
 *  thermodynamic temperature @f$ T @f$.
 *
 *  Using the @p options flag one can switch between regular (#VRNA_PBACKTRACK_DEFAULT) backtracing
 *  mode, and non-redundant sampling (#VRNA_PBACKTRACK_NON_REDUNDANT) along the lines of
 *  @rstinline :cite:t:`michalik:2017` @endrst.
 *
 *  In contrast to vrna_pbacktrack5_cb() this function allows for resuming a previous
 *  sampling round in specialized Boltzmann sampling, such as non-redundant backtracking.
 *  For that purpose, the user passes the address of a Boltzmann sampling data structure
 *  (#vrna_pbacktrack_mem_t) which will be re-used in each round of sampling, i.e. each
 *  successive call to vrna_pbacktrack5_resume_cb() or vrna_pbacktrack5_resume().
 *
 *  A successive sample call to this function may look like:
 * @code{.c}
 * vrna_pbacktrack_mem_t nonredundant_memory = NULL;
 *
 * // sample the first 100 structures
 * vrna_pbacktrack5_resume(fc,
 *                         100,
 *                         fc->length,
 *                         &nonredundant_memory,
 *                         options);
 *
 * // sample another 500 structures
 * vrna_pbacktrack5_resume(fc,
 *                         500,
 *                         fc->length,
 *                         &nonredundant_memory,
 *                         options);
 *
 * // release memory occupied by the non-redundant memory data structure
 * vrna_pbacktrack_mem_free(nonredundant_memory);
 * @endcode
 *
 *  @pre    Unique multiloop decomposition has to be active upon creation of @p fc with vrna_fold_compound()
 *          or similar. This can be done easily by passing vrna_fold_compound() a model details parameter
 *          with vrna_md_t.uniq_ML = 1.<br>
 *          vrna_pf() has to be called first to fill the partition function matrices
 *
 *  @note This function is polymorphic. It accepts #vrna_fold_compound_t of type
 *        #VRNA_FC_TYPE_SINGLE, and #VRNA_FC_TYPE_COMPARATIVE.
 *
 *  @warning  In non-redundant sampling mode (#VRNA_PBACKTRACK_NON_REDUNDANT), this function may
 *            not yield the full number of requested samples. This may happen if
 *            a)  the number of requested structures is larger than the total number
 *                of structuresin the ensemble,
 *            b)  numeric instabilities prevent the backtracking function to enumerate
 *                structures with high free energies, or
 *            c)  any other error occurs.
 *
 *  @see  vrna_pbacktrack5_resume_cb(), vrna_pbacktrack5_cb(), vrna_pbacktrack_resume(),
 *        #vrna_pbacktrack_mem_t, #VRNA_PBACKTRACK_DEFAULT, #VRNA_PBACKTRACK_NON_REDUNDANT,
 *        vrna_pbacktrack_mem_free
 *
 *  @param  fc            The fold compound data structure
 *  @param  num_samples   The size of the sample set, i.e. number of structures
 *  @param  start         The start of  the subsequence to consider, i.e. 5'-end position(1-based)
 *  @param  end           The end of the subsequence to consider, i.e. 3'-end position (1-based)
 *  @param  nr_mem        The address of the Boltzmann sampling memory data structure
 *  @param  options       A bitwise OR-flag indicating the backtracing mode.
 *  @return               A set of secondary structure samples in dot-bracket notation terminated by NULL (or NULL on error)
 */
char **
vrna_pbacktrack_sub_resume(vrna_fold_compound_t  *fc,
                           unsigned int          num_samples,
                           unsigned int          start,
                           unsigned int          end,
                           vrna_pbacktrack_mem_t *nr_mem,
                           unsigned int          options);


/**
 *  @brief Obtain a set of secondary structure samples for a subsequence from the Boltzmann ensemble according their probability
 *
 *  Perform a probabilistic (stochastic) backtracing in the partition function DP arrays
 *  to obtain a set of @p num_samples secondary structures. The parameter @p length specifies
 *  the length of the substructure starting from the 5' end.
 *
 *  Any structure @f$ s @f$ with free energy @f$ E(s) @f$ is picked from the Boltzmann distributed
 *  ensemble according to its probability
 *
 *  @f[
 *  p(s) = \frac{exp(-E(s) / kT)}{Z}
 *  @f]
 *
 *  with partition function @f$ Z = \sum_s exp(-E(s) / kT) @f$, Boltzmann constant @f$ k @f$ and
 *  thermodynamic temperature @f$ T @f$.
 *
 *  Using the @p options flag one can switch between regular (#VRNA_PBACKTRACK_DEFAULT) backtracing
 *  mode, and non-redundant sampling (#VRNA_PBACKTRACK_NON_REDUNDANT) along the lines of
 *  @rstinline :cite:t:`michalik:2017` @endrst.
 *
 *  In contrast to vrna_pbacktrack5_resume() this function yields the structure samples
 *  through a callback mechanism.
 *
 *  A successive sample call to this function may look like:
 * @code{.c}
 * vrna_pbacktrack_mem_t nonredundant_memory = NULL;
 *
 * // sample the first 100 structures
 * vrna_pbacktrack5_resume_cb(fc,
 *                            100,
 *                            fc->length,
 *                            &callback_function,
 *                            (void *)&callback_data,
 *                            &nonredundant_memory,
 *                            options);
 *
 * // sample another 500 structures
 * vrna_pbacktrack5_resume_cb(fc,
 *                            500,
 *                            fc->length,
 *                            &callback_function,
 *                            (void *)&callback_data,
 *                            &nonredundant_memory,
 *                            options);
 *
 * // release memory occupied by the non-redundant memory data structure
 * vrna_pbacktrack_mem_free(nonredundant_memory);
 * @endcode
 *
 *  @pre    Unique multiloop decomposition has to be active upon creation of @p fc with vrna_fold_compound()
 *          or similar. This can be done easily by passing vrna_fold_compound() a model details parameter
 *          with vrna_md_t.uniq_ML = 1.<br>
 *          vrna_pf() has to be called first to fill the partition function matrices
 *
 *  @note This function is polymorphic. It accepts #vrna_fold_compound_t of type
 *        #VRNA_FC_TYPE_SINGLE, and #VRNA_FC_TYPE_COMPARATIVE.
 *
 *  @warning  In non-redundant sampling mode (#VRNA_PBACKTRACK_NON_REDUNDANT), this function may
 *            not yield the full number of requested samples. This may happen if
 *            a)  the number of requested structures is larger than the total number
 *                of structuresin the ensemble,
 *            b)  numeric instabilities prevent the backtracking function to enumerate
 *                structures with high free energies, or
 *            c)  any other error occurs.
 *
 *  @see  vrna_pbacktrack5_resume(), vrna_pbacktrack5_cb(), vrna_pbacktrack_resume_cb(),
 *        #vrna_pbacktrack_mem_t, #VRNA_PBACKTRACK_DEFAULT, #VRNA_PBACKTRACK_NON_REDUNDANT,
 *        vrna_pbacktrack_mem_free
 *
 *  @param  fc            The fold compound data structure
 *  @param  num_samples   The size of the sample set, i.e. number of structures
 *  @param  start         The start of  the subsequence to consider, i.e. 5'-end position(1-based)
 *  @param  end           The end of the subsequence to consider, i.e. 3'-end position (1-based)
 *  @param  cb            The callback that receives the sampled structure
 *  @param  data          A data structure passed through to the callback @p cb
 *  @param  nr_mem        The address of the Boltzmann sampling memory data structure
 *  @param  options       A bitwise OR-flag indicating the backtracing mode.
 *  @return               The number of structures actually backtraced
 */
unsigned int
vrna_pbacktrack_sub_resume_cb(vrna_fold_compound_t             *fc,
                              unsigned int                     num_samples,
                              unsigned int                     start,
                              unsigned int                     end,
                              vrna_bs_result_f cb,
                              void                             *data,
                              vrna_pbacktrack_mem_t            *nr_mem,
                              unsigned int                     options);


/**
 *  @brief  Release memory occupied by a Boltzmann sampling memory data structure
 *
 *  @see  #vrna_pbacktrack_mem_t, vrna_pbacktrack5_resume(), vrna_pbacktrack5_resume_cb(),
 *        vrna_pbacktrack_resume(), vrna_pbacktrack_resume_cb()
 *
 *  @param  s   The non-redundancy memory data structure
 */
void
vrna_pbacktrack_mem_free(vrna_pbacktrack_mem_t s);


/**@}*/


#endif
