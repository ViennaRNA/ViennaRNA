#ifndef VIENNA_RNA_PACKAGE_PROBING_STRATEGY_DEIGAN_H
#define VIENNA_RNA_PACKAGE_PROBING_STRATEGY_DEIGAN_H


#include <ViennaRNA/probing/basic.h>
#include <ViennaRNA/math/functions.h>

/**
 *  @file     ViennaRNA/probing/strategy_deigan.h
 *  @ingroup  probing_data_strategy
 *  @brief    This module provides the API for the Deigan et al. 2009 strategy to convert structure probing data
 */

/**
 *  @addtogroup probing_data_strategy_deigan
 *  @{
 */


/**
 *  @brief  Default parameter for slope `m` as used in method of @rstinline :cite:t:`deigan:2009` @endrst
 *
 *  @see    vrna_probing_data_deigan(), vrna_probing_data_deigan_comparative(),
 *          #VRNA_PROBING_METHOD_DEIGAN2009_DEFAULT_b
 */
#define VRNA_PROBING_METHOD_DEIGAN2009_DEFAULT_m                  1.8


/**
 *  @brief  Default parameter for intercept `b` as used in method of @rstinline :cite:t:`deigan:2009` @endrst
 *
 *  @see    vrna_probing_data_deigan(), vrna_probing_data_deigan_comparative(),
 *          #VRNA_PROBING_METHOD_DEIGAN2009_DEFAULT_m
 */
#define VRNA_PROBING_METHOD_DEIGAN2009_DEFAULT_b                  -0.6


double *
vrna_probing_strategy_deigan(vrna_fold_compound_t *fc,
                             const double         *data,
                             size_t               data_size,
                             unsigned int         target,
                             void                 *options);


void *
vrna_probing_strategy_deigan_options(double                       m,
                                     double                       b,
                                     double                       max_value,
                                     vrna_math_fun_dbl_f          cb_preprocess,
                                     vrna_math_fun_dbl_opt_t      cb_preprocess_opt,
                                     vrna_math_fun_dbl_opt_free_f cb_preprocess_opt_free);


void
vrna_probing_strategy_deigan_options_free(void *options);


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
 *  results in a confirmatory *bonus* free energy for correctly predicted base pairs.
 *  Since the energy evaluation of a base pair stack involves two pairs, the pseudo
 *  energies are added for all four contributing nucleotides. Consequently, the
 *  energy term is applied twice for pairs inside a helix and only once for pairs
 *  adjacent to other structures. For all other loop types the energy model remains
 *  unchanged even when the experimental data highly disagrees with a certain motif.
 *
 *  @note For further details, we refer to @rstinline :cite:t:`deigan:2009` @endrst.
 *
 *  @see  #vrna_probing_data_t, vrna_probing_data_free(), vrna_sc_probing(),
 *        vrna_probing_data_deigan_comparative(),
 *        vrna_probing_data_zarringhalam(), vrna_probing_data_zarringhalam_comparative(),
 *        vrna_probing_data_eddy(), vrna_probing_data_eddy_comparative()
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
vrna_probing_data_deigan(const double *reactivities,
                         unsigned int n,
                         double       m,
                         double       b);


vrna_probing_data_t
vrna_probing_data_deigan_trans(const double                 *reactivities,
                               unsigned int                 n,
                               double                       m,
                               double                       b,
                               vrna_math_fun_dbl_f          cb_preprocess,
                               vrna_math_fun_dbl_opt_t      cb_preprocess_opt,
                               vrna_math_fun_dbl_opt_free_f cb_preprocess_opt_free);


/**
 *  @brief  Prepare (multiple) probing data according to Deigan et al. 2009 method for comparative structure predictions
 *
 *  Similar to vrna_probing_data_deigan(), this function prepares a data structure to be used
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
 *        vrna_probing_data_deigan(),
 *        vrna_probing_data_zarringhalam(), vrna_probing_data_zarringhalam_comparative(),
 *        vrna_probing_data_eddy(), vrna_probing_data_eddy_comparative(),
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
vrna_probing_data_deigan_comparative(const double       **reactivities,
                                     const unsigned int *n,
                                     unsigned int       n_seq,
                                     double             *ms,
                                     double             *bs,
                                     unsigned int       multi_params);


vrna_probing_data_t
vrna_probing_data_deigan_trans_comparative(const double                 **reactivities,
                                           const unsigned int           *n,
                                           unsigned int                 n_seq,
                                           double                       *ms,
                                           double                       *bs,
                                           unsigned int                 multi_params,
                                           vrna_math_fun_dbl_f          *cb_preprocess,
                                           vrna_math_fun_dbl_opt_t      *cb_preprocess_opt,
                                           vrna_math_fun_dbl_opt_free_f *cb_preprocess_opt_free);


/**
 *  @}
 */

#endif
