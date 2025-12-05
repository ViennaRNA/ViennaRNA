#ifndef VIENNA_RNA_PACKAGE_PROBING_STRATEGY_ZARRINGHALAM_H
#define VIENNA_RNA_PACKAGE_PROBING_STRATEGY_ZARRINGHALAM_H


#include <ViennaRNA/probing/basic.h>
#include <ViennaRNA/math/functions.h>

/**
 *  @file     ViennaRNA/probing/strategy_zarringhalam.h
 *  @ingroup  probing_data_strategy
 *  @brief    This module provides the API for the Zarringhalam et al. 20212 strategy to convert structure probing data
 */

/**
 *  @addtogroup probing_data_strategy_zarringhalam
 *  @{
 */


/**
 *  @brief  Default parameter `beta` as used in method of @rstinline :cite:t:`zarringhalam:2012` @endrst
 *
 *  @see    vrna_probing_data_zarringhalam(), vrna_probing_data_zarringhalam_comparative(),
 *          #VRNA_PROBING_METHOD_ZARRINGHALAM2012_DEFAULT_conversion, #VRNA_PROBING_METHOD_ZARRINGHALAM2012_DEFAULT_probability
 */
#define VRNA_PROBING_METHOD_ZARRINGHALAM2012_DEFAULT_beta         0.89


/**
 *  @brief  Default conversion method of probing data into probabilities as used in method of @rstinline :cite:t:`zarringhalam:2012` @endrst
 *
 *  @see    vrna_probing_data_zarringhalam(), vrna_probing_data_zarringhalam_comparative(),
 *          #VRNA_PROBING_METHOD_ZARRINGHALAM2012_DEFAULT_beta, #VRNA_PROBING_METHOD_ZARRINGHALAM2012_DEFAULT_probability
 */
#define VRNA_PROBING_METHOD_ZARRINGHALAM2012_DEFAULT_conversion   "Os1.6i-2.29"


/**
 *  @brief  Default probability value for missing data in method of @rstinline :cite:t:`zarringhalam:2012` @endrst
 *
 *  @see    vrna_probing_data_zarringhalam(), vrna_probing_data_zarringhalam_comparative(),
 *          #VRNA_PROBING_METHOD_ZARRINGHALAM2012_DEFAULT_beta, #VRNA_PROBING_METHOD_ZARRINGHALAM2012_DEFAULT_conversion
 */
#define VRNA_PROBING_METHOD_ZARRINGHALAM2012_DEFAULT_probability  0.5

double *
vrna_probing_strategy_zarringhalam(vrna_fold_compound_t *fc,
                                   const double         *data,
                                   size_t               data_size,
                                   unsigned int         target,
                                   void                 *options);


void *
vrna_probing_strategy_zarringhalam_options(double                       beta,
                                           double                       default_probability,
                                           double                       max_value,
                                           vrna_math_fun_dbl_f          cb_preprocess,
                                           vrna_math_fun_dbl_opt_t      cb_preprocess_opt,
                                           vrna_math_fun_dbl_opt_free_f cb_preprocess_opt_free);


void
vrna_probing_strategy_zarringhalam_options_free(void *option);


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
 *        vrna_probing_data_zarringhalam_comparative(),
 *        vrna_probing_data_deigan(), vrna_probing_data_deigan_comparative(),
 *        vrna_probing_data_eddy(), vrna_probing_data_eddy_comparative()
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
vrna_probing_data_zarringhalam(const double *reactivities,
                               unsigned int n,
                               double       beta,
                               const char   *pr_conversion,
                               double       pr_default);


vrna_probing_data_t
vrna_probing_data_zarringhalam_trans(const double                 *reactivities,
                                     unsigned int                 n,
                                     double                       beta,
                                     double                       pr_default,
                                     vrna_math_fun_dbl_f          cb_preprocess,
                                     vrna_math_fun_dbl_opt_t      cb_preprocess_opt,
                                     vrna_math_fun_dbl_opt_free_f cb_preprocess_opt_free);


/**
 *  @brief  Prepare probing data according to Zarringhalam et al. 2012 method for comparative structure predictions
 *
 *  Similar to vrna_probing_data_zarringhalam(), this function prepares a data structure to be used
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
 *        vrna_probing_data_zarringhalam_comparative(),
 *        vrna_probing_data_deigan(), vrna_probing_data_deigan_comparative(),
 *        vrna_probing_data_eddy(), vrna_probing_data_eddy_comparative(),
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
vrna_probing_data_zarringhalam_comparative(const double **reactivities,
                                           unsigned int *n,
                                           unsigned int n_seq,
                                           double       *betas,
                                           const char   **pr_conversions,
                                           double       *pr_defaults,
                                           unsigned int multi_params);


vrna_probing_data_t
vrna_probing_data_zarringhalam_trans_comparative(const double             **reactivities,
                                                 unsigned int             *n,
                                                 unsigned int             n_seq,
                                                 double                   *betas,
                                                 double                   *pr_defaults,
                                                 unsigned int             multi_params,
                                                 vrna_math_fun_dbl_f      *cb_preprocess,
                                                 vrna_math_fun_dbl_opt_t  *cb_preprocess_opt,
                                                 vrna_math_fun_dbl_opt_free_f *
                                                 cb_preprocess_opt_free);


/**
 *  @}
 */

#endif
