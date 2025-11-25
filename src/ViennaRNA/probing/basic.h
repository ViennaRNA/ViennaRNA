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


#include <ViennaRNA/probing/strategies.h>



#define VRNA_REACTIVITY_MISSING                                   -999.

#define VRNA_PROBING_DATA_DEFAULT                                 0
#define VRNA_PROBING_DATA_SINGLE_STRATEGY                         (1 << 0)
#define VRNA_PROBING_DATA_SINGLE_WEIGHT                           (1 << 1)

#define VRNA_PROBING_DATA_LINEAR_TARGET_STACK                     (1 << 4)
#define VRNA_PROBING_DATA_LINEAR_TARGET_UP                        (1 << 5)
#define VRNA_PROBING_DATA_LINEAR_TARGET_BP                        (1 << 6)

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


/**
 *  @}
 */

#endif
