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
 *  @see    vrna_probing_data_deigan(), vrna_probing_data_deigan_comparative(),
 *          vrna_probing_data_zarringhalam(), vrna_probing_data_zarringhalam_comparative(),
 *          vrna_sc_probing(), vrna_probing_data_free()
 */
typedef struct vrna_probing_data_s *vrna_probing_data_t;


#define VRNA_REACTIVITY_MISSING                                   -999.



/**
 *  @}
 */


/**
 *  @addtogroup probing_data_strategy
 *  @{
 */

#define VRNA_PROBING_DATA_DEFAULT                                 0
#define VRNA_PROBING_DATA_SINGLE_STRATEGY                         (1 << 0)
#define VRNA_PROBING_DATA_SINGLE_WEIGHT                           (1 << 1)

#define VRNA_PROBING_DATA_LINEAR_TARGET_STACK                     (1 << 4)
#define VRNA_PROBING_DATA_LINEAR_TARGET_UP                        (1 << 5)
#define VRNA_PROBING_DATA_LINEAR_TARGET_BP                        (1 << 6)


/**
 *  @brief  Prototype of a strategy to derive pseudo energies from linear structure probing data
 *
 *  @param  fc        The fold compound the probing data will be applied to
 *  @param  data      The structure probing data (1-based array of reactivity values)
 *  @param  data_size The size of @p data, i.e. the total number of reactivity values
 *  @param  target    The structural context for which pseudo energies are requested from the strategy
 *  @param  options   An arbitrary data structure the strategy requires for working on the data
 *  @return           A pointer to an array of pseudo energies ready to be included as soft constraints
 */
typedef double *(*vrna_probing_strategy_f)(vrna_fold_compound_t *fc,
                                           const double         *data,
                                           size_t               data_size,
                                           unsigned int         target,
                                           void                 *options);

/**
 *  @brief probing data conversion flag for comparative structure predictions indicating no parameter to be sequence specific
 *
 *  @see  vrna_probing_data_deigan_comparative(), vrna_probing_data_zarringhalam_comparative(),
 *        #VRNA_PROBING_METHOD_MULTI_PARAMS_1, #VRNA_PROBING_METHOD_MULTI_PARAMS_2, #VRNA_PROBING_METHOD_MULTI_PARAMS_3,
 *        #VRNA_PROBING_METHOD_MULTI_PARAMS_DEFAULT
 */
#define VRNA_PROBING_METHOD_MULTI_PARAMS_0                        0U


/**
 *  @brief probing data conversion flag for comparative structure predictions indicating 1st parameter to be sequence specific
 *
 *  @see  vrna_probing_data_deigan_comparative(), vrna_probing_data_zarringhalam_comparative(),
 *        #VRNA_PROBING_METHOD_MULTI_PARAMS_0, #VRNA_PROBING_METHOD_MULTI_PARAMS_2, #VRNA_PROBING_METHOD_MULTI_PARAMS_3,
 *        #VRNA_PROBING_METHOD_MULTI_PARAMS_DEFAULT
 */
#define VRNA_PROBING_METHOD_MULTI_PARAMS_1                        1U


/**
 *  @brief probing data conversion flag for comparative structure predictions indicating 2nd parameter to be sequence specific
 *
 *  @see  vrna_probing_data_deigan_comparative(), vrna_probing_data_zarringhalam_comparative(),
 *        #VRNA_PROBING_METHOD_MULTI_PARAMS_0, #VRNA_PROBING_METHOD_MULTI_PARAMS_1, #VRNA_PROBING_METHOD_MULTI_PARAMS_3,
 *        #VRNA_PROBING_METHOD_MULTI_PARAMS_DEFAULT
 */
#define VRNA_PROBING_METHOD_MULTI_PARAMS_2                        2U


/**
 *  @brief probing data conversion flag for comparative structure predictions indicating 3rd parameter to be sequence specific
 *
 *  @see  vrna_probing_data_deigan_comparative(), vrna_probing_data_zarringhalam_comparative(),
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
 *  @see vrna_probing_data_deigan_comparative(), vrna_probing_data_zarringhalam_comparative()
 */
#define VRNA_PROBING_METHOD_MULTI_PARAMS_DEFAULT                  VRNA_PROBING_METHOD_MULTI_PARAMS_0


/**
 *  @}
 */



/**
 *  @addtogroup probing_data
 *  @{
 */


/**
 *  @brief  Apply probing data (e.g. SHAPE) to guide the structure prediction
 *
 *  @see  #vrna_probing_data_t, vrna_probing_data_free(),
 *        vrna_probing_data_deigan(), vrna_probing_data_deigan_comparative(),
 *        vrna_probing_data_zarringhalam(), vrna_probing_data_zarringhalam_comparative(),
 *        vrna_probing_data_eddy(), vrna_probing_data_eddy_comparative()
 *
 *  @param  fc      The #vrna_fold_compound_t the probing data should be applied to in subsequent computations
 *  @param  data    The prepared probing data and probing data integration strategy
 *  @return The number of probing data sets applied, 0 upon any error
 */
int
vrna_sc_probing(vrna_fold_compound_t  *fc,
                vrna_probing_data_t   data);


vrna_probing_data_t
vrna_probing_data_linear(const double              *data,
                         unsigned int              data_length,
                         double                    data_weight,
                         vrna_probing_strategy_f   strategy_cb,
                         void                      *strategy_cb_options,
                         vrna_auxdata_free_f       strategy_cb_options_free);


vrna_probing_data_t
vrna_probing_data_linear_multi(const double              **data,
                               unsigned int              data_size,
                               const unsigned int        *data_lengths,
                               const double              *data_weights,
                               vrna_probing_strategy_f   *strategy_cbs,
                               void                      **strategy_cbs_options,
                               vrna_auxdata_free_f       *strategy_cbs_options_free,
                               unsigned int              options);




/**
 *  @brief  Free memory occupied by the (prepared) probing data
 *
 *  @see    #vrna_probing_data_t, vrna_sc_probing(),
 *          vrna_probing_data_deigan(), vrna_probing_data_deigan_comparative(),
 *          vrna_probing_data_zarringhalam(), vrna_probing_data_zarringhalam_comparative(),
 *          vrna_probing_data_eddy(), vrna_probing_data_eddy_comparative()
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


/**
 *  @}
 */

#endif
