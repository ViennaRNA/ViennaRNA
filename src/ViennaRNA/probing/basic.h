#ifndef VIENNA_RNA_PACKAGE_PROBING_H
#define VIENNA_RNA_PACKAGE_PROBING_H

#include <ViennaRNA/fold_compound.h>

/**
 *  @file     ViennaRNA/probing/basic.h
 *  @ingroup  probing_data
 *  @brief    This module provides function to incorporate RNA structure
 *            probing data, e.g. SHAPE reactivities, into the folding recursions
 *            by means of soft constraints
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
 *
 *  This flag indicates to use an implementation that distinguishes two classes of structural context, in particular
 *  paired and unpaired positions.
 */
#define VRNA_PROBING_METHOD_EDDY2014_2                            4U


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


#define VRNA_PROBING_DATA_CHECK_SEQUENCE                          1U


/**
 *  @brief  Apply probing data (e.g. SHAPE) to guide the structure prediction
 *
 *  @see  #vrna_probing_data_t, vrna_probing_data_free(),
 *        vrna_probing_data_Deigan2009(), vrna_probing_data_Deigan2009_comparative(),
 *        vrna_probing_data_Zarringhalam2012(), vrna_probing_data_Zarringhalam2012_comparative(),
 *        vrna_probing_data_Eddy2014_2(), vrna_probing_data_Eddy2014_2_comparative()
 *
 *  @param  fc      The #vrna_fold_compound_t the probing data should be applied to in subsequent computations
 *  @param  data    The prepared probing data and probing data integration strategy
 *  @return The number of probing data sets applied, 0 upon any error
 */
int
vrna_sc_probing(vrna_fold_compound_t  *fc,
                vrna_probing_data_t   data);


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
 *  @see  #vrna_probing_data_t, vrna_probing_data_free(), vrna_sc_probing(),
 *        vrna_probing_data_Deigan2009_comparative(),
 *        vrna_probing_data_Zarringhalam2012(), vrna_probing_data_Zarringhalam2012_comparative(),
 *        vrna_probing_data_Eddy2014_2(), vrna_probing_data_Eddy2014_2_comparative()
 *
 *  @param  reactivities  1-based array of per-nucleotide probing data, e.g. SHAPE reactivities
 *  @param  n             The length of the @p reactivities list
 *  @param  m             The slope used for the probing data to soft constraints conversion strategy
 *  @param  b             The intercept used for the probing data to soft constraints conversion strategy
 *  @return               A pointer to a data structure containing the probing data and any preparations
 *                        necessary to use it in vrna_sc_probing() according to the method of
 *                        @rstinline :cite:t:`deigan:2009` @endrst or @b NULL on any error.
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
 *  in the MSA. Individual probing data is always expected to be specified in sequence coordinates,
 *  i.e. without considering gaps in the MSA. Therefore, each set of @p reactivities may have
 *  a different length as specified the parameter @p n.
 *  In addition, each set of probing data may undergo the conversion using different parameters
 *  @f$ m @f$ and @f$ b @f$. Whether or not multiple sets of conversion parameters are provided
 *  must be specified using the @p multi_params flag parameter. Use #VRNA_PROBING_METHOD_MULTI_PARAMS_1
 *  to indicate that @p ms points to an array of slopes for each sequence. Along with that,
 *  #VRNA_PROBING_METHOD_MULTI_PARAMS_2 indicates that @p bs is pointing to an array of intercepts
 *  for each sequence. Bitwise-OR of the two values renders both parameters to be sequence specific.
 *
 *  @note For further details, we refer to @rstinline :cite:t:`deigan:2009` @endrst.
 *
 *  @see  #vrna_probing_data_t, vrna_probing_data_free(), vrna_sc_probing(),
 *        vrna_probing_data_Deigan2009(),
 *        vrna_probing_data_Zarringhalam2012(), vrna_probing_data_Zarringhalam2012_comparative(),
 *        vrna_probing_data_Eddy2014_2(), vrna_probing_data_Eddy2014_2_comparative(),
 *        #VRNA_PROBING_METHOD_MULTI_PARAMS_0, #VRNA_PROBING_METHOD_MULTI_PARAMS_1,
 *        #VRNA_PROBING_METHOD_MULTI_PARAMS_2, #VRNA_PROBING_METHOD_MULTI_PARAMS_DEFAULT
 *
 *  @param  reactivities  0-based array of 1-based arrays of per-nucleotide probing data, e.g. SHAPE reactivities
 *  @param  n             0-based array of lengths of the @p reactivities lists
 *  @param  n_seq         The number of sequences in the MSA
 *  @param  ms            0-based array of the slopes used for the probing data to soft constraints conversion strategy or the address of a single slope value to be applied for all data
 *  @param  bs            0-based array of the intercepts used for the probing data to soft constraints conversion strategy or the address of a single intercept value to be applied for all data
 *  @param  multi_params  A flag indicating what is passed through parameters @p ms and @p bs
 *  @return               A pointer to a data structure containing the probing data and any preparations
 *                        necessary to use it in vrna_sc_probing() according to the method of
 *                        @rstinline :cite:t:`deigan:2009` @endrst or @b NULL on any error.
 */
vrna_probing_data_t
vrna_probing_data_Deigan2009_comparative(const double       **reactivities,
                                         const unsigned int *n,
                                         unsigned int       n_seq,
                                         double             *ms,
                                         double             *bs,
                                         unsigned int       multi_params);


/**
 *  @brief  Prepare probing data according to Zarringhalam et al. 2012 method
 *
 *  Prepares a data structure to be used with vrna_sc_probing() to directed RNA
 *  folding using the method of @rstinline :cite:t:`zarringhalam:2012` @endrst.
 *
 *  This method first converts the observed probing data of nucleotide @f$ i @f$ into a
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
 *  @see  #vrna_probing_data_t, vrna_probing_data_free(), vrna_sc_probing(),
 *        vrna_probing_data_Zarringhalam2012_comparative(),
 *        vrna_probing_data_Deigan2009(), vrna_probing_data_Deigan2009_comparative(),
 *        vrna_probing_data_Eddy2014_2(), vrna_probing_data_Eddy2014_2_comparative()
 *
 *  @param  reactivities      1-based array of per-nucleotide probing data, e.g. SHAPE reactivities
 *  @param  n                 The length of the @p reactivities list
 *  @param  beta              The scaling factor @f$ \beta @f$ of the conversion function
 *  @param  pr_conversion     A flag that specifies how to convert reactivities to probabilities
 *  @param  pr_default        The default probability for a nucleotide where reactivity data is missing for
 *  @return                   A pointer to a data structure containing the probing data and any preparations
 *                            necessary to use it in vrna_sc_probing() according to the method of
 *                            @rstinline :cite:t:`zarringhalam:2012` @endrst or @b NULL on any error.
 */
vrna_probing_data_t
vrna_probing_data_Zarringhalam2012(const double *reactivities,
                                   unsigned int n,
                                   double       beta,
                                   const char   *pr_conversion,
                                   double       pr_default);


/**
 *  @brief  Prepare probing data according to Zarringhalam et al. 2012 method for comparative structure predictions
 *
 *  Similar to vrna_probing_data_Zarringhalam2012(), this function prepares a data structure to be used
 *  with vrna_sc_probing() to guide RNA folding using the method of @rstinline :cite:t:`zarringhalam:2012` @endrst.
 *
 *  This functions purpose is to allow for adding multiple probing data as required for
 *  comparative structure predictions over multiple sequence alignments (MSA) with @p n_seq
 *  sequences. For that purpose, @p reactivities can be provided for any of the sequences
 *  in the MSA. Individual probing data is always expected to be specified in sequence coordinates,
 *  i.e. without considering gaps in the MSA. Therefore, each set of @p reactivities may have
 *  a different length as specified the parameter @p n.
 *  In addition, each set of probing data may undergo the conversion using different parameters
 *  @f$ beta @f$. Additionally, the probing data to probability conversions strategy and default
 *  values for missing data can be specified in a sequence-based manner. Whether or not multiple conversion
 *  parameters are provided must be specified using the @p multi_params flag parameter. Use
 *  #VRNA_PROBING_METHOD_MULTI_PARAMS_1 to indicate that @p betas points to an array of @f$ beta @f$ values
 *  for each sequence. #VRNA_PROBING_METHOD_MULTI_PARAMS_2 indicates that @p pr_conversions is pointing to
 *  an array of probing data to probability conversion strategies, and #VRNA_PROBING_METHOD_MULTI_PARAMS_3
 *  indicates multiple default probabilities for missing data. Bitwise-OR of the three values renders all of
 *  them to be sequence specific.
 *
 *  @note For further details, we refer to @rstinline :cite:t:`zarringhalam:2012` @endrst
 *
 *  @see  #vrna_probing_data_t, vrna_probing_data_free(), vrna_sc_probing(),
 *        vrna_probing_data_Zarringhalam2012_comparative(),
 *        vrna_probing_data_Deigan2009(), vrna_probing_data_Deigan2009_comparative(),
 *        vrna_probing_data_Eddy2014_2(), vrna_probing_data_Eddy2014_2_comparative(),
 *        #VRNA_PROBING_METHOD_MULTI_PARAMS_0, #VRNA_PROBING_METHOD_MULTI_PARAMS_1,
 *        #VRNA_PROBING_METHOD_MULTI_PARAMS_2, #VRNA_PROBING_METHOD_MULTI_PARAMS_3,
 *        #VRNA_PROBING_METHOD_MULTI_PARAMS_DEFAULT
 *
 *  @param  reactivities      0-based array of 1-based arrays of per-nucleotide probing data, e.g. SHAPE reactivities
 *  @param  n                 0-based array of lengths of the @p reactivities lists
 *  @param  n_seq             The number of sequences in the MSA
 *  @param  betas             0-based array with scaling factors @f$ \beta @f$ of the conversion function or the address of a scaling factor to be applied for all data
 *  @param  pr_conversions    0-based array of flags that specifies how to convert reactivities to probabilities or the address of a conversion strategy to be applied for all data
 *  @param  pr_defaults       0-based array of default probabilities for a nucleotide where reactivity data is missing for or the address of a single default probability to be applied for all data
 *  @param  multi_params      A flag indicating what is passed through parameters @p betas, @p pr_conversions, and @p pr_defaults
 *  @return                   A pointer to a data structure containing the probing data and any preparations
 *                            necessary to use it in vrna_sc_probing() according to the method of
 *                            @rstinline :cite:t:`zarringhalam:2012` @endrst or @b NULL on any error.
 */
vrna_probing_data_t
vrna_probing_data_Zarringhalam2012_comparative(const double **reactivities,
                                               unsigned int *n,
                                               unsigned int n_seq,
                                               double       *betas,
                                               const char   **pr_conversions,
                                               double       *pr_defaults,
                                               unsigned int multi_params);


/**
 *  @brief  Add probing data as soft constraints (Eddy/RNAprob-2 method)
 *
 *  This approach of probing data directed RNA folding uses the probability framework
 *  proposed by @rstinline :cite:t:`eddy:2014` @endrst:
 *
 *  @f[ \Delta G_{\text{data}}(i) = - RT\ln(\mathbb{P}(\text{data}(i)\mid x_i\pi_i)) @f]
 *
 *  to convert probing data to pseudo energies for given nucleotide @f$ x_i @f$ and
 *  class probability @f$ \pi_i @f$ at position @f$ i @f$. The conditional probability
 *  is taken from a prior-distribution of probing data for the respective classes.
 *
 *  Here, the method distinguishes exactly two different classes of structural context,
 *  (i) unpaired and (ii) paired positions, following the lines of the RNAprob-2
 *  method of @rstinline :cite:t:`deng:2016` @endrst. The reactivity distribution is computed
 *  using Gaussian kernel density estimation (KDE) with bandwidth @f$ h @f$ computed using
 *  Scott factor
 *
 *  @f[ h = n^{-\frac{1}{5}} @f]
 *
 *  where @f$ n @f$ is the number of data points of the prior distribution.
 *
 *  @note For further details, we refer to @rstinline :cite:t:`eddy:2014` @endrst and
 *        @rstinline :cite:t:`deng:2016` @endrst.
 *
 *  @see  #vrna_probing_data_t, vrna_probing_data_free(), vrna_sc_probing(),
 *        vrna_probing_data_Eddy2014_2_comparative(),
 *        vrna_probing_data_Deigan2009(), vrna_probing_data_Deigan2009_comparative(),
 *        vrna_probing_data_Zarringhalam2012(), vrna_probing_data_Zarringhalam2012_comparative(),
 *
 *  @param  reactivities  A 1-based vector of probing data, e.g. normalized SHAPE reactivities
 *  @param  n             Length of @p reactivities
 *  @param  unpaired_data Pointer to an array of probing data for unpaired nucleotides
 *  @param  unpaired_len  Length of @p unpaired_data
 *  @param  paired_data   Pointer to an array of probing data for paired nucleotides
 *  @param  paired_len    Length of @p paired_data
 *  @return               A pointer to a data structure containing the probing data and any preparations
 *                        necessary to use it in vrna_sc_probing() according to the method of
 *                        @rstinline :cite:t:`eddy:2014` @endrst or @b NULL on any error.
 */
vrna_probing_data_t
vrna_probing_data_Eddy2014_2(const double *reactivities,
                             unsigned int n,
                             const double *unpaired_data,
                             unsigned int unpaired_len,
                             const double *paired_data,
                             unsigned int paired_len);


/**
 *  @brief  Add probing data as soft constraints (Eddy/RNAprob-2 method) for comparative structure predictions
 *
 *  Similar to vrna_probing_data_Eddy2014_2(), this function prepares a data structure for
 *  probing data directed RNA folding. It uses the probability framework proposed by
 *  @rstinline :cite:t:`eddy:2014` @endrst:
 *
 *  @f[ \Delta G_{\text{data}}(i) = - RT\ln(\mathbb{P}(\text{data}(i)\mid x_i\pi_i)) @f]
 *
 *  to convert probing data to pseudo energies for given nucleotide @f$ x_i @f$ and
 *  class probability @f$ \pi_i @f$ at position @f$ i @f$. The conditional probability
 *  is taken from a prior-distribution of probing data for the respective classes.
 *
 *  This functions purpose is to allow for adding multiple probing data as required for
 *  comparative structure predictions over multiple sequence alignments (MSA) with @p n_seq
 *  sequences. For that purpose, @p reactivities can be provided for any of the sequences
 *  in the MSA. Individual probing data is always expected to be specified in sequence coordinates,
 *  i.e. without considering gaps in the MSA. Therefore, each set of @p reactivities may have
 *  a different length as specified the parameter @p n.
 *  In addition, each set of probing data may undergo the conversion using different prior distributions
 *  for unpaired and paired nucleotides. Whether or not multiple sets of conversion priors are provided
 *  must be specified using the @p multi_params flag parameter. Use #VRNA_PROBING_METHOD_MULTI_PARAMS_1
 *  to indicate that @p unpaired_datas points to an array of unpaired probing data for each sequence.
 *  Similarly, #VRNA_PROBING_METHOD_MULTI_PARAMS_2 indicates that @p paired_datas is pointing to an array
 *  paired probing data for each sequence. Bitwise-OR of the two values renders both parameters to be
 *  sequence specific.
 *
 *  @note For further details, we refer to @rstinline :cite:t:`eddy:2014` @endrst and
 *        @rstinline :cite:t:`deng:2016` @endrst.
 *
 *  @see  #vrna_probing_data_t, vrna_probing_data_free(), vrna_sc_probing(),
 *        vrna_probing_data_Eddy2014_2(),
 *        vrna_probing_data_Deigan2009(), vrna_probing_data_Deigan2009_comparative(),
 *        vrna_probing_data_Zarringhalam2012(), vrna_probing_data_Zarringhalam2012_comparative(),
 *        #VRNA_PROBING_METHOD_MULTI_PARAMS_0, #VRNA_PROBING_METHOD_MULTI_PARAMS_1,
 *        #VRNA_PROBING_METHOD_MULTI_PARAMS_2, #VRNA_PROBING_METHOD_MULTI_PARAMS_DEFAULT
 *
 *  @param  reactivities    0-based array of 1-based arrays of per-nucleotide probing data, e.g. SHAPE reactivities
 *  @param  n               0-based array of lengths of the @p reactivities lists
 *  @param  n_seq           The number of sequences in the MSA
 *  @param  unpaired_datas  0-based array of 0-based arrays with probing data for unpaired nucleotides or address of a single array of such data
 *  @param  unpaired_lens   0-based array of lengths for each probing data array in @p unpaired_datas
 *  @param  paired_datas    0-based array of 0-based arrays with probing data for paired nucleotides or address of a single array of such data
 *  @param  paired_lens     0-based array of lengths for each probing data array in @p paired_data
 *  @param  multi_params    A flag indicating what is passed through parameters @p unpaired_datas and @p paired_datas
 *  @return                 A pointer to a data structure containing the probing data and any preparations
 *                          necessary to use it in vrna_sc_probing() according to the method of
 *                          @rstinline :cite:t:`eddy:2014` @endrst or @b NULL on any error.
 */
vrna_probing_data_t
vrna_probing_data_Eddy2014_2_comparative(const double **reactivities,
                                         unsigned int *n,
                                         unsigned int n_seq,
                                         const double **unpaired_datas,
                                         unsigned int *unpaired_lens,
                                         const double **paired_datas,
                                         unsigned int *paired_lens,
                                         unsigned int multi_params);


/**
 *  @brief  Free memory occupied by the (prepared) probing data
 *
 *  @see    #vrna_probing_data_t, vrna_sc_probing(),
 *          vrna_probing_data_Deigan2009(), vrna_probing_data_Deigan2009_comparative(),
 *          vrna_probing_data_Zarringhalam2012(), vrna_probing_data_Zarringhalam2012_comparative(),
 *          vrna_probing_data_Eddy2014_2(), vrna_probing_data_Eddy2014_2_comparative()
 */
void
vrna_probing_data_free(vrna_probing_data_t d);


/**
 *  @brief Convert SHAPE reactivity values to probabilities for being unpaired
 *
 *  This function parses the informations from a given file and stores the result
 *  in the pre-allocated string sequence and the #FLT_OR_DBL array values.
 *
 *  @ingroup SHAPE_reactivities
 *
 *  @see vrna_file_SHAPE_read()
 *
 *  @param shape_conversion String defining the method used for the conversion process
 *  @param values           Pointer to an array of SHAPE reactivities
 *  @param length           Length of the array of SHAPE reactivities
 *  @param default_value    Result used for position with invalid/missing reactivity values
 */
int
vrna_sc_SHAPE_to_pr(const char  *shape_conversion,
                    double      *values,
                    int         length,
                    double      default_value);


double **
vrna_probing_data_load_n_distribute(unsigned int  n_seq,
                                    unsigned int  *ns,
                                    const char    **sequences,
                                    const char    **file_names,
                                    const int     *file_name_association,
                                    unsigned int  options);


/**
 *  @}
 */

#endif
