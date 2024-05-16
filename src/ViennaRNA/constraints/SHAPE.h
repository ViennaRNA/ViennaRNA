#ifndef VIENNA_RNA_PACKAGE_CONSTRAINTS_SHAPE_H
#define VIENNA_RNA_PACKAGE_CONSTRAINTS_SHAPE_H

#include <ViennaRNA/fold_compound.h>

/**
 *  @file constraints/SHAPE.h
 *  @ingroup SHAPE_reactivities
 *  @brief  This module provides function to incorporate RNA structure
 *          probing data, e.g. SHAPE reactivities, into the folding recursions
 *          by means of soft constraints
 */

/**
 *  @addtogroup probing_data
 *  @{
 */


/**
 *  @brief  A data structure that contains RNA structure probing data and specifies how this data
 *          is to be integrated into structure predictions
 *
 *  @see    vrna_probing_data_Deigan2009(), vrna_probing_data_Deigan2009_comparative(),
 *          vrna_probing_data_Zarringhalam2012(), vrna_probing_data_Zarringhalam2012_comparative(),
 *          vrna_sc_probing(), vrna_probing_data_free()
 */
typedef struct vrna_probing_data_s *vrna_probing_data_t;


/**
 *  @brief  A flag indicating probing data conversion method of @rstinline :cite:t:`deigan:2009` @endrst
 */
#define VRNA_PROBING_METHOD_DEIGAN2009                            1U

/**
 *  @brief  Default parameter for slope `m` as used in method of @rstinline :cite:t:`deigan:2009` @endrst
 *
 *  @see    vrna_probing_data_Deigan2009(), vrna_probing_data_Deigan2009_comparative(),
 *          #VRNA_PROBING_METHOD_DEIGAN2009_DEFAULT_b
 */
#define VRNA_PROBING_METHOD_DEIGAN2009_DEFAULT_m                  1.8


/**
 *  @brief  Default parameter for intercept `b` as used in method of @rstinline :cite:t:`deigan:2009` @endrst
 *
 *  @see    vrna_probing_data_Deigan2009(), vrna_probing_data_Deigan2009_comparative(),
 *          #VRNA_PROBING_METHOD_DEIGAN2009_DEFAULT_m
 */
#define VRNA_PROBING_METHOD_DEIGAN2009_DEFAULT_b                  -0.6


/**
 *  @brief  A flag indicating probing data conversion method of @rstinline :cite:t:`zarringhalam:2012` @endrst
 */
#define VRNA_PROBING_METHOD_ZARRINGHALAM2012                      2U


/**
 *  @brief  Default parameter `beta` as used in method of @rstinline :cite:t:`zarringhalam:2012` @endrst
 *
 *  @see    vrna_probing_data_Zarringhalam2012(), vrna_probing_data_Zarringhalam2012_comparative(),
 *          #VRNA_PROBING_METHOD_ZARRINGHALAM2012_DEFAULT_conversion, #VRNA_PROBING_METHOD_ZARRINGHALAM2012_DEFAULT_probability
 */
#define VRNA_PROBING_METHOD_ZARRINGHALAM2012_DEFAULT_beta         0.89


/**
 *  @brief  Default conversion method of probing data into probabilities as used in method of @rstinline :cite:t:`zarringhalam:2012` @endrst
 *
 *  @see    vrna_probing_data_Zarringhalam2012(), vrna_probing_data_Zarringhalam2012_comparative(),
 *          #VRNA_PROBING_METHOD_ZARRINGHALAM2012_DEFAULT_beta, #VRNA_PROBING_METHOD_ZARRINGHALAM2012_DEFAULT_probability
 */
#define VRNA_PROBING_METHOD_ZARRINGHALAM2012_DEFAULT_conversion   "Os1.6i-2.29"


/**
 *  @brief  Default probability value for missing data in method of @rstinline :cite:t:`zarringhalam:2012` @endrst
 *
 *  @see    vrna_probing_data_Zarringhalam2012(), vrna_probing_data_Zarringhalam2012_comparative(),
 *          #VRNA_PROBING_METHOD_ZARRINGHALAM2012_DEFAULT_beta, #VRNA_PROBING_METHOD_ZARRINGHALAM2012_DEFAULT_conversion
 */
#define VRNA_PROBING_METHOD_ZARRINGHALAM2012_DEFAULT_probability  0.5

/**
 *  @brief  A flag indicating probing data conversion method of @rstinline :cite:t:`washietl:2012` @endrst
 */
#define VRNA_PROBING_METHOD_WASHIETL2012                          3U


/**
 *  @brief  A flag indicating probing data conversion method of @rstinline :cite:t:`eddy:2014` @endrst
 */
#define VRNA_PROBING_METHOD_EDDY2014                              4U


/**
 *  @brief probing data conversion flag for comparative structure predictions indicating no parameter to be sequence specific
 *
 *  @see  vrna_probing_data_Deigan2009_comparative(), vrna_probing_data_Zarringhalam2012_comparative(),
 *        #VRNA_PROBING_METHOD_MULTI_PARAMS_1, #VRNA_PROBING_METHOD_MULTI_PARAMS_2, #VRNA_PROBING_METHOD_MULTI_PARAMS_3,
 *        #VRNA_PROBING_METHOD_MULTI_PARAMS_DEFAULT
 */
#define VRNA_PROBING_METHOD_MULTI_PARAMS_0                        0U


/**
 *  @brief probing data conversion flag for comparative structure predictions indicating 1st parameter to be sequence specific
 *
 *  @see  vrna_probing_data_Deigan2009_comparative(), vrna_probing_data_Zarringhalam2012_comparative(),
 *        #VRNA_PROBING_METHOD_MULTI_PARAMS_0, #VRNA_PROBING_METHOD_MULTI_PARAMS_2, #VRNA_PROBING_METHOD_MULTI_PARAMS_3,
 *        #VRNA_PROBING_METHOD_MULTI_PARAMS_DEFAULT
 */
#define VRNA_PROBING_METHOD_MULTI_PARAMS_1                        1U


/**
 *  @brief probing data conversion flag for comparative structure predictions indicating 2nd parameter to be sequence specific
 *
 *  @see  vrna_probing_data_Deigan2009_comparative(), vrna_probing_data_Zarringhalam2012_comparative(),
 *        #VRNA_PROBING_METHOD_MULTI_PARAMS_0, #VRNA_PROBING_METHOD_MULTI_PARAMS_1, #VRNA_PROBING_METHOD_MULTI_PARAMS_3,
 *        #VRNA_PROBING_METHOD_MULTI_PARAMS_DEFAULT
 */
#define VRNA_PROBING_METHOD_MULTI_PARAMS_2                        2U


/**
 *  @brief probing data conversion flag for comparative structure predictions indicating 3rd parameter to be sequence specific
 *
 *  @see  vrna_probing_data_Deigan2009_comparative(), vrna_probing_data_Zarringhalam2012_comparative(),
 *        #VRNA_PROBING_METHOD_MULTI_PARAMS_0, #VRNA_PROBING_METHOD_MULTI_PARAMS_1, #VRNA_PROBING_METHOD_MULTI_PARAMS_2,
 *        #VRNA_PROBING_METHOD_MULTI_PARAMS_DEFAULT
 */
#define VRNA_PROBING_METHOD_MULTI_PARAMS_3                        4U


/**
 *  @brief probing data conversion flag for comparative structure predictions indicating default parameter settings
 *
 *  Essentially, this setting indicates that all probing data is to be converted using the same
 *  parameters. Use any combination of #VRNA_PROBING_METHOD_MULTI_PARAMS_1, #VRNA_PROBING_METHOD_MULTI_PARAMS_2,
 *  #VRNA_PROBING_METHOD_MULTI_PARAMS_3, and so on to indicate that the first, second, third, or other
 *  parameter is sequence specific.
 *
 *  @see vrna_probing_data_Deigan2009_comparative(), vrna_probing_data_Zarringhalam2012_comparative()
 */
#define VRNA_PROBING_METHOD_MULTI_PARAMS_DEFAULT                  VRNA_PROBING_METHOD_MULTI_PARAMS_0


/**
 *  @brief  Apply probing data (e.g. SHAPE) to guide the structure prediction
 *
 *  @see  #vrna_probing_data_t, vrna_probing_data_Deigan2009(), vrna_probing_data_Deigan2009_comparative(),
 *        vrna_probing_data_Zarringhalam2012(), vrna_probing_data_Zarringhalam2012_comparative()
 *        
 *
 */
int
vrna_sc_probing(vrna_fold_compound_t *fc,
                vrna_probing_data_t    data);


/**
 *  @brief  Prepare probing data according to Deigan et al. 2009 method
 *
 *  Prepares a data structure to be used with vrna_sc_probing() to directed RNA
 *  folding using the simple linear ansatz
 *
 *  @f[ \Delta G_{\text{SHAPE}}(i) = m \ln(\text{SHAPE reactivity}(i)+1)+ b @f]
 *
 *  to convert probing data, e.g. SHAPE reactivity values, to pseudo energies whenever a
 *  nucleotide @f$ i @f$ contributes to a stacked pair. A positive slope @f$ m @f$
 *  penalizes high reactivities in paired regions, while a negative intercept @f$ b @f$
 *  results in a confirmatory ``bonus'' free energy for correctly predicted base pairs.
 *  Since the energy evaluation of a base pair stack involves two pairs, the pseudo
 *  energies are added for all four contributing nucleotides. Consequently, the
 *  energy term is applied twice for pairs inside a helix and only once for pairs
 *  adjacent to other structures. For all other loop types the energy model remains
 *  unchanged even when the experimental data highly disagrees with a certain motif.
 *
 *  @note For further details, we refer to @rstinline :cite:t:`deigan:2009` @endrst.
 *
 *  @see  #vrna_probing_data_t, vrna_sc_probing(), vrna_probing_data_Deigan2009_comparative(),
 *        vrna_probing_data_Zarringhalam2012(), vrna_probing_data_Zarringhalam2012_comparative()
 *
 *  @param  reactivities  1-based array of per-nucleotide probing data, e.g. SHAPE reactivities
 *  @param  n             The length of the @p reactivities list
 *  @param  m             The slope used for the probing data to soft constraints conversion strategy
 *  @param  b             The intercept used for the probing data to soft constraints conversion strategy
 *  @return               A data structure containing the probing data and any preparations necessary to use
 *                        it in vrna_sc_probing() according to the method of @rstinline :cite:t:`deigan:2009` @endrst.
 */
vrna_probing_data_t
vrna_probing_data_Deigan2009(const double *reactivities,
                           unsigned int n,
                           double       m,
                           double       b);


/**
 *  @brief  Prepare (multiple) probing data according to Deigan et al. 2009 method for comparative structure predictions
 *
 *  Similar to vrna_probing_data_Deigan2009(), this function prepares a data structure to be used
 *  with vrna_sc_probing() to directed RNA folding using the simple linear ansatz
 *
 *  @f[ \Delta G_{\text{SHAPE}}(i) = m \ln(\text{SHAPE reactivity}(i)+1)+ b @f]
 *
 *  to convert probing data, e.g. SHAPE reactivity values, to pseudo energies whenever a
 *  nucleotide @f$ i @f$ contributes to a stacked pair.
 *  This functions purpose is to allow for adding multiple probing data as required for
 *  comparative structure predictions over multiple sequence alignments (MSA) with @p n_seq
 *  sequences. For that purpose, @p reactivities can be provided for any of the sequences
 *  in the MSA. Individual probing data is always expected to be speficied in sequence coordinates,
 *  i.e. without considering gaps in the MSA. Therefore, each set of @p reactivities may have
 *  a different length as specified the parameter $p n.
 *  In addition, each set of probing data may undergo the conversion using different parameters
 *  @f$ m @f$ and @f$ b @f$. Whether or not multiple sets of conversion parameters are provided
 *  must be speficied using the @p multi_params flag parameter. Use #VRNA_PROBING_METHOD_MULTI_PARAMS_1
 *  to indicate that @p ms points to an array of slopes for each sequence. Along with that,
 *  VRNA_PROBING_METHOD_MULTI_PARAMS_2 indicates that @p bs is pointing to an array of intercepts
 *  for each sequence. Bitwise-OR of the two values renders both parameters to be sequence specific.
 *
 *  @note For further details, we refer to @rstinline :cite:t:`deigan:2009` @endrst.
 *
 *  @see  #vrna_probing_data_t, vrna_sc_probing(), vrna_probing_data_Deigan2009(),
 *        vrna_probing_data_Zarringhalam2012(), vrna_probing_data_Zarringhalam2012_comparative(),
 *        VRNA_PROBING_METHOD_MULTI_PARAMS_0, VRNA_PROBING_METHOD_MULTI_PARAMS_1,
 *        VRNA_PROBING_METHOD_MULTI_PARAMS_2, VRNA_PROBING_METHOD_MULTI_PARAMS_DEFAULT
 *
 *  @param  reactivities  0-based array of 1-based arrays of per-nucleotide probing data, e.g. SHAPE reactivities
 *  @param  n             0-based array of lengths of the @p reactivities lists
 *  @params n_seq         The number of sequences in the MSA
 *  @param  ms            0-based array of the slopes used for the probing data to soft constraints conversion strategy or the address of a single slope value to be applied for all data
 *  @param  bs            0-based array of the intercepts used for the probing data to soft constraints conversion strategy or the address of a single intercept value to be applied for all data
 *  @param  multi_params  A flag indicating what is passed through parameters @p ms and @p bs
 *  @return               A data structure containing the probing data and any preparations necessary to use
 *                        it in vrna_sc_probing() according to the method of @rstinline :cite:t:`deigan:2009` @endrst.
 */
vrna_probing_data_t
vrna_probing_data_Deigan2009_comparative(const double       **reactivities,
                                       const unsigned int *n,
                                       unsigned int       n_seq,
                                       double             *ms,
                                       double             *bs,
                                       unsigned int       multi_params);

vrna_probing_data_t
vrna_probing_data_Zarringhalam2012(const double *reactivities,
                                 unsigned int n,
                                 double       beta,
                                 const char   *pr_conversion,
                                 double       pr_default);

vrna_probing_data_t
vrna_probing_data_Zarringhalam2012_comparative(const double **reactivities,
                                             unsigned int *n,
                                             unsigned int n_seq,
                                             double       *betas,
                                             const char   **pr_conversions,
                                             double       *pr_defaults,
                                             unsigned int multi_params);

void
vrna_probing_data_free(vrna_probing_data_t d);


/* End group probing_data */
/**@}*/

/**
 *  @ingroup SHAPE_reactivities
 */
void
vrna_constraints_add_SHAPE(vrna_fold_compound_t *fc,
                           const char           *shape_file,
                           const char           *shape_method,
                           const char           *shape_conversion,
                           int                  verbose,
                           unsigned int         constraint_type);


/**
 *  @ingroup SHAPE_reactivities
 */
void
vrna_constraints_add_SHAPE_ali(vrna_fold_compound_t *fc,
                               const char           *shape_method,
                               const char           **shape_files,
                               const int            *shape_file_association,
                               int                  verbose,
                               unsigned int         constraint_type);


/**
 *  @brief  Add SHAPE reactivity data as soft constraints (Deigan et al. method)
 *
 *  This approach of SHAPE directed RNA folding uses the simple linear ansatz
 *
 *  @f[ \Delta G_{\text{SHAPE}}(i) = m \ln(\text{SHAPE reactivity}(i)+1)+ b @f]
 *
 *  to convert SHAPE reactivity values to pseudo energies whenever a
 *  nucleotide @f$ i @f$ contributes to a stacked pair. A positive slope @f$ m @f$
 *  penalizes high reactivities in paired regions, while a negative intercept @f$ b @f$
 *  results in a confirmatory ``bonus'' free energy for correctly predicted base pairs.
 *  Since the energy evaluation of a base pair stack involves two pairs, the pseudo
 *  energies are added for all four contributing nucleotides. Consequently, the
 *  energy term is applied twice for pairs inside a helix and only once for pairs
 *  adjacent to other structures. For all other loop types the energy model remains
 *  unchanged even when the experimental data highly disagrees with a certain motif.
 *
 *  @note For further details, we refer to @rstinline :cite:t:`deigan:2009` @endrst.
 *
 *  @see  vrna_sc_remove(), vrna_sc_add_SHAPE_zarringhalam(), vrna_sc_minimize_pertubation()
 *
 *  @ingroup SHAPE_reactivities
 *  @param  fc            The #vrna_fold_compound_t the soft constraints are associated with
 *  @param  reactivities  A vector of normalized SHAPE reactivities
 *  @param  m             The slope of the conversion function
 *  @param  b             The intercept of the conversion function
 *  @param  options       The options flag indicating how/where to store the soft constraints
 *  @return               1 on successful extraction of the method, 0 on errors
 */
int
vrna_sc_add_SHAPE_deigan(vrna_fold_compound_t *fc,
                         const double         *reactivities,
                         double               m,
                         double               b,
                         unsigned int         options);


/**
 *  @brief  Add SHAPE reactivity data from files as soft constraints for consensus structure prediction (Deigan et al. method)
 *
 *  @ingroup SHAPE_reactivities
 *  @param  fc            The #vrna_fold_compound_t the soft constraints are associated with
 *  @param  shape_files   A set of filenames that contain normalized SHAPE reactivity data
 *  @param  shape_file_association  An array of integers that associate the files with sequences in the alignment
 *  @param  m             The slope of the conversion function
 *  @param  b             The intercept of the conversion function
 *  @param  options       The options flag indicating how/where to store the soft constraints
 *  @return               1 on successful extraction of the method, 0 on errors
 */
int
vrna_sc_add_SHAPE_deigan_ali(vrna_fold_compound_t *fc,
                             const char           **shape_files,
                             const int            *shape_file_association,
                             double               m,
                             double               b,
                             unsigned int         options);


/**
 *  @brief  Add SHAPE reactivity data as soft constraints (Zarringhalam et al. method)
 *
 *  This method first converts the observed SHAPE reactivity of nucleotide @f$ i @f$ into a
 *  probability @f$ q_i @f$ that position @f$ i @f$ is unpaired by means of a non-linear map.
 *  Then pseudo-energies of the form
 *
 *  @f[ \Delta G_{\text{SHAPE}}(x,i) = \beta\ |x_i - q_i| @f]
 *
 *  are computed, where @f$ x_i=0 @f$ if position @f$ i @f$ is unpaired and @f$ x_i=1 @f$
 *  if @f$ i @f$ is paired in a given secondary structure. The parameter @f$ \beta @f$ serves as
 *  scaling factor. The magnitude of discrepancy between prediction and experimental observation
 *  is represented by @f$ |x_i - q_i| @f$.
 *
 *  @note For further details, we refer to @rstinline :cite:t:`zarringhalam:2012` @endrst
 *
 *  @see  vrna_sc_remove(), vrna_sc_add_SHAPE_deigan(), vrna_sc_minimize_pertubation()
 *
 *  @ingroup SHAPE_reactivities
 *  @param  fc                The #vrna_fold_compound_t the soft constraints are associated with
 *  @param  reactivities      A vector of normalized SHAPE reactivities
 *  @param  b                 The scaling factor @f$ \beta @f$ of the conversion function
 *  @param  default_value     The default value for a nucleotide where reactivity data is missing for
 *  @param  shape_conversion  A flag that specifies how to convert reactivities to probabilities
 *  @param  options           The options flag indicating how/where to store the soft constraints
 *  @return                   1 on successful extraction of the method, 0 on errors
 */
int
vrna_sc_add_SHAPE_zarringhalam(vrna_fold_compound_t *fc,
                               const double         *reactivities,
                               double               b,
                               double               default_value,
                               const char           *shape_conversion,
                               unsigned int         options);


/**
 *  @brief  Parse a character string and extract the encoded SHAPE reactivity conversion
 *          method and possibly the parameters for conversion into pseudo free energies
 *
 *  @ingroup soft_cosntraints
 *
 *  @param  method_string   The string that contains the encoded SHAPE reactivity conversion method
 *  @param  method          A pointer to the memory location where the method character will be stored
 *  @param  param_1         A pointer to the memory location where the first parameter of the corresponding method will be stored
 *  @param  param_2         A pointer to the memory location where the second parameter of the corresponding method will be stored
 *  @return                 1 on successful extraction of the method, 0 on errors
 */
int
vrna_sc_SHAPE_parse_method(const char *method_string,
                           char       *method,
                           float      *param_1,
                           float      *param_2);


/**
 *  @brief Convert SHAPE reactivity values to probabilities for being unpaired
 *
 *  This function parses the informations from a given file and stores the result
 *  in the preallocated string sequence and the #FLT_OR_DBL array values.
 *
 *  @ingroup SHAPE_reactivities
 *
 *  @see vrna_file_SHAPE_read()
 *
 *  @param shape_conversion String definining the method used for the conversion process
 *  @param values           Pointer to an array of SHAPE reactivities
 *  @param length           Length of the array of SHAPE reactivities
 *  @param default_value    Result used for position with invalid/missing reactivity values
 */
int
vrna_sc_SHAPE_to_pr(const char  *shape_conversion,
                    double      *values,
                    int         length,
                    double      default_value);


/**
 *  @brief  Add SHAPE reactivity data as soft constraints (Eddy/RNAprob-2 method)
 *
 *  This approach of SHAPE directed RNA folding uses the probability framework proposed by Eddy
 *
 *  @f[ \Delta G_{\text{SHAPE}}(i) = - RT\ln(\mathbb{P}(\text{SHAPE reactivity}(i)\mid x_i\pi_i))+ @f]
 *
 *  to convert SHAPE reactivity values to pseudo energies for given nucleotide @f$ x_i @f$ and
 *  pairedness @f$ \pi_i @f$ at position @f$ i @f$. The reactivity distribution is computed
 *  using Gaussian kernel density estimation (KDE) with bandwidth @f$ h @f$ computed using
 *  Scott factor
 * 
 * @f[ h = n^{-\frac{1}{5}} @f]
 *
 * where @f$ n @f$ is the number of data points.
 *
 *
 *  @ingroup SHAPE_reactivities
 *  @param  fc            The #vrna_fold_compound_t the soft constraints are associated with
 *  @param  reactivities  A vector of normalized SHAPE reactivities
 *  @param  unpaired_nb   Length of the array of unpaired SHAPE reactivities
 *  @param  unpaired_data Pointer to an array of unpaired SHAPE reactivities
 *  @param  paired_nb     Length of the array of paired SHAPE reactivities
 *  @param  paired_data   Pointer to an array of paired SHAPE reactivities
 *  @return               1 on successful extraction of the method, 0 on errors
 */
int
vrna_sc_add_SHAPE_eddy_2(vrna_fold_compound_t *fc,
                         const double         *reactivities,
                         int                  unpaired_nb,
                         const double         *unpaired_data,
                         int                  paired_nb,
                         const double         *paired_data);

#endif
