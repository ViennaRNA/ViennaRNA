#ifndef VIENNA_RNA_PACKAGE_PROBING_STRATEGIES_H
#define VIENNA_RNA_PACKAGE_PROBING_STRATEGIES_H


typedef double *(*vrna_probing_strategy_f)(const double *data,
                                           size_t       data_size,
                                           void         *options);

#include <ViennaRNA/probing/basic.h>
#include <ViennaRNA/probing/transform.h>

/**
 *  @file     ViennaRNA/probing/strategies.h
 *  @ingroup  probing_data
 *  @brief    This module provides different strategies for the incorporation of RNA structure
 *            probing data, e.g. SHAPE reactivities, into the folding recursions
 *            by means of soft constraints
 */

/**
 *  @addtogroup probing_data
 *  @{
 */


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


double *
vrna_probing_strategy_deigan(const double *data,
                             size_t       data_size,
                             void         *options);


void *
vrna_probing_strategy_deigan_options(double                   m,
                                     double                   b,
                                     vrna_probing_transform_f cb_preprocess,
                                     void                     *cb_preprocess_opt,
                                     vrna_auxdata_free_f      cb_preprocess_opt_free);

void
vrna_probing_strategy_deigan_options_free(void *options);


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

double *
vrna_probing_strategy_zarringhalam_up(const double *data,
                                      size_t       data_size,
                                      void         *options);


void *
vrna_probing_strategy_zarringhalam_options(double                   beta,
                                           double                   default_probability,
                                           double                   max_value,
                                           vrna_probing_transform_f cb_preprocess,
                                           void                     *cb_preprocess_opt,
                                           vrna_auxdata_free_f      cb_preprocess_opt_free);


void
vrna_probing_strategy_zarringhalam_options_free(void *option);


/**
 *  @}
 */

#endif
