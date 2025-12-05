#ifndef VIENNA_RNA_PACKAGE_PROBING_STRATEGY_EDDY_H
#define VIENNA_RNA_PACKAGE_PROBING_STRATEGY_EDDY_H


#include <ViennaRNA/probing/basic.h>
#include <ViennaRNA/math/functions.h>

/**
 *  @file     ViennaRNA/probing/strategy_eddy.h
 *  @ingroup  probing_data_strategy
 *  @brief    This module provides the API for the Eddy 2014 strategy to convert structure probing data
 */

/**
 *  @addtogroup probing_data_strategy_eddy
 *  @{
 */


/**
 *  @brief  Default options for the @rstinline :cite:t:`eddy:2014` @endrst: probing data conversion strategy
 *
 *  @see  vrna_probing_strategy_eddy_options(), vrna_probing_data_eddy(), vrna_probing_data_eddy_trans(),
 *        vrna_probing_data_eddy_comparative(), vrna_probing_data_eddy_trans_comparative()
 */
#define VRNA_PROBING_STRATEGY_EDDY_OPTIONS_DEFAULT              0


/**
 *  @brief  Prevent temperature dependent energy rescaling in @rstinline :cite:t:`eddy:2014` @endrst strategy
 *
 *  This option flag forces the probing data conversion strategy to always use the same thermodynamic
 *  temperature @f$ T @f$, no matter what temperature the predictions are made for.
 *
 *  @see  vrna_probing_strategy_eddy_options(), vrna_probing_data_eddy(), vrna_probing_data_eddy_trans(),
 *        vrna_probing_data_eddy_comparative(), vrna_probing_data_eddy_trans_comparative()
 */
#define VRNA_PROBING_STRATEGY_EDDY_NO_TEMPERATURE_RESCALING     (1 << 5)


double *
vrna_probing_strategy_eddy(vrna_fold_compound_t *fc,
                           const double         *data,
                           size_t               data_size,
                           unsigned int         target,
                           void                 *options);


void *
vrna_probing_strategy_eddy_options(double                       temperature,
                                   unsigned char                options,
                                   const double                 *prior_unpaired,
                                   size_t                       prior_unpaired_size,
                                   const double                 *prior_paired,
                                   size_t                       prior_paired_size,
                                   vrna_math_fun_dbl_f          cb_preprocess,
                                   vrna_math_fun_dbl_opt_t      cb_preprocess_opt,
                                   vrna_math_fun_dbl_opt_free_f cb_preprocess_opt_free);


void
vrna_probing_strategy_eddy_options_free(void *options);


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
 *        vrna_probing_data_eddy_comparative(),
 *        vrna_probing_data_deigan(), vrna_probing_data_deigan_comparative(),
 *        vrna_probing_data_zarringhalam(), vrna_probing_data_zarringhalam_comparative(),
 *        #VRNA_PROBING_STRATEGY_EDDY_NO_TEMPERATURE_RESCALING,
 *        #VRNA_PROBING_STRATEGY_EDDY_OPTIONS_DEFAULT
 *
 *  @param  reactivities  A 1-based vector of probing data, e.g. normalized SHAPE reactivities
 *  @param  n             Length of @p reactivities
 *  @param  temperature   The thermodynamic temperature @f$ T @f$
 *  @param  options       Options bit flags to change the behavior of the strategy
 *  @param  unpaired_data Pointer to an array of probing data for unpaired nucleotides
 *  @param  unpaired_len  Length of @p unpaired_data
 *  @param  paired_data   Pointer to an array of probing data for paired nucleotides
 *  @param  paired_len    Length of @p paired_data
 *  @return               A pointer to a data structure containing the probing data and any preparations
 *                        necessary to use it in vrna_sc_probing() according to the method of
 *                        @rstinline :cite:t:`eddy:2014` @endrst or @b NULL on any error.
 */
struct vrna_probing_data_s *
vrna_probing_data_eddy(const double   *reactivities,
                       unsigned int   n,
                       double         temperature,
                       unsigned char  options,
                       const double   *unpaired_data,
                       unsigned int   unpaired_len,
                       const double   *paired_data,
                       unsigned int   paired_len);


struct vrna_probing_data_s *
vrna_probing_data_eddy_trans(const double                 *reactivities,
                             unsigned int                 n,
                             double                       temperature,
                             unsigned char                options,
                             const double                 *unpaired_data,
                             unsigned int                 unpaired_len,
                             const double                 *paired_data,
                             unsigned int                 paired_len,
                             vrna_math_fun_dbl_f          cb_preprocess,
                             vrna_math_fun_dbl_opt_t      cb_preprocess_opt,
                             vrna_math_fun_dbl_opt_free_f cb_preprocess_opt_free);


/**
 *  @brief  Add probing data as soft constraints (Eddy/RNAprob-2 method) for comparative structure predictions
 *
 *  Similar to vrna_probing_data_eddy(), this function prepares a data structure for
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
 *        vrna_probing_data_eddy(),
 *        vrna_probing_data_deigan(), vrna_probing_data_deigan_comparative(),
 *        vrna_probing_data_zarringhalam(), vrna_probing_data_zarringhalam_comparative(),
 *        #VRNA_PROBING_METHOD_MULTI_PARAMS_0, #VRNA_PROBING_METHOD_MULTI_PARAMS_1,
 *        #VRNA_PROBING_METHOD_MULTI_PARAMS_2, #VRNA_PROBING_METHOD_MULTI_PARAMS_DEFAULT
 *
 *  @param  reactivities    0-based array of 1-based arrays of per-nucleotide probing data, e.g. SHAPE reactivities
 *  @param  n               0-based array of lengths of the @p reactivities lists
 *  @param  n_seq           The number of sequences in the MSA
 *  @param  temperature   The thermodynamic temperature @f$ T @f$
 *  @param  options       Options bit flags to change the behavior of the strategy
 *  @param  unpaired_datas  0-based array of 0-based arrays with probing data for unpaired nucleotides or address of a single array of such data
 *  @param  unpaired_lens   0-based array of lengths for each probing data array in @p unpaired_datas
 *  @param  paired_datas    0-based array of 0-based arrays with probing data for paired nucleotides or address of a single array of such data
 *  @param  paired_lens     0-based array of lengths for each probing data array in @p paired_data
 *  @param  multi_params    A flag indicating what is passed through parameters @p unpaired_datas and @p paired_datas
 *  @return                 A pointer to a data structure containing the probing data and any preparations
 *                          necessary to use it in vrna_sc_probing() according to the method of
 *                          @rstinline :cite:t:`eddy:2014` @endrst or @b NULL on any error.
 */
struct vrna_probing_data_s *
vrna_probing_data_eddy_comparative(const double       **reactivities,
                                   const unsigned int *n,
                                   unsigned int       n_seq,
                                   double             temperature,
                                   unsigned char      options,
                                   const double       **unpaired_datas,
                                   const unsigned int *unpaired_lens,
                                   const double       **paired_datas,
                                   const unsigned int *paired_lens,
                                   unsigned int       multi_params);


struct vrna_probing_data_s *
vrna_probing_data_eddy_trans_comparative(const double                 **reactivities,
                                         const unsigned int           *n,
                                         unsigned int                 n_seq,
                                         double                       temperature,
                                         unsigned char                options,
                                         const double                 **unpaired_datas,
                                         const unsigned int           *unpaired_lens,
                                         const double                 **paired_datas,
                                         const unsigned int           *paired_lens,
                                         unsigned int                 multi_params,
                                         vrna_math_fun_dbl_f          *cb_preprocess,
                                         vrna_math_fun_dbl_opt_t      *cb_preprocess_opt,
                                         vrna_math_fun_dbl_opt_free_f *cb_preprocess_opt_free);


/**
 *  @}
 */

#endif
