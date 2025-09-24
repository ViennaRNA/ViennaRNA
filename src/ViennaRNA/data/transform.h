#ifndef VIENNA_RNA_PACKAGE_PROBING_TRANSFORM_H
#define VIENNA_RNA_PACKAGE_PROBING_TRANSFORM_H

#include <ViennaRNA/fold_compound.h>

/**
 *  @file     ViennaRNA/data/transform.h
 *  @ingroup  data_utils
 *  @brief    This module provides function to transform linear data
 */

/**
 *  @addtogroup data_utils  Handling and Manipulation of Data
 *  @{
 */


/**
 *  @brief  Linear data transform callback option
 */
typedef   void                  *vrna_data_lin_trans_opt_t;

/**
 *  @brief  Callback function to release any memory occupied by the linear data transform callback option #vrna_data_lin_trans_opt_t
 */
typedef   vrna_auxdata_free_f   vrna_data_lin_trans_opt_free_f;


/**
 * @brief Callback function prototype to transform a single value of linear data
 *
 * This is the prototype for callback functions used by the linear data transformation function
 * vrna_data_lin_transform(). During data transform, this callback will be executed for each value
 * of the linear data provided to vrna_data_lin_transform(). The callback then is responsible to
 * actually transform a provided data value. It will be provided by the calling function with two
 * parameters: @p value which is the value that is about to be transformed and @p options, which
 * is a pointer to a data structure that holds any means required to ensure that the transformation
 * callback can transform the provide data. The callback function then has to return the transformed
 * value.
 *
 * @callback
 * @parblock
 * This callback handles the transformation of a single value of linear data.
 * @endparblock
 *
 * @see vrna_data_lin_transform(),
 *      vrna_data_transform_method_bin(), vrna_data_transform_method_lm(), vrna_data_transform_method_log()
 *
 * @param value     The value that is to be transformed
 * @param options   A pointer to a transformation specific data structure
 * @return          The transformed value
 */
typedef double (*vrna_data_lin_trans_f) (double                    value,
                                         vrna_data_lin_trans_opt_t options);


/**
 *  @brief  Transform an array of linear data
 *
 *  This function transforms an array of linear data into another array of linear data
 *  of the same size. For that purpose, it utilizes a callback mechanism that transforms
 *  each single value.
 *
 *  @note When @p data is @c NULL or @p data_size equals @c 0, the function returns @c NULL.
 *        Moreover, the function simply provides a copy of the linear data if the transformation
 *        callback @p transform_cb is not provided, i.e. if it is @c NULL.
 *
 *  @see  #vrna_data_lin_trans_f, #vrna_data_lin_trans_opt_t,
 *        vrna_data_transform_method_bin(), vrna_data_transform_method_lm(), vrna_data_transform_method_log()
 *
 *  @param  data            A pointer to an array of linear data
 *  @param  data_size       The size of the array @p data is pointing to
 *  @param  transform_cb    The data transformation callback (maybe @c NULL)
 *  @param  transform_opt   The options that need to be passed through to the transformation callback @p transform_cb (maybe @c NULL)
 */
double *
vrna_data_lin_transform(const double              *data,
                        size_t                    data_size,
                        vrna_data_lin_trans_f     transform_cb,
                        vrna_data_lin_trans_opt_t transform_opt);


#define VRNA_TRANSFORM_BIN_OPTION_PROJECT               (1 << 0)
#define VRNA_TRANSFORM_BIN_OPTION_MAP_OUTOF_UPPERBOUND  (1 << 1)
#define VRNA_TRANSFORM_BIN_OPTION_MAP_OUTOF_LOWERBOUND  (1 << 2)
#define VRNA_TRANSFORM_BIN_OPTION_DEFAULT               0

#define VRNA_TRANSFORM_LM_OPTION_LOG                    (1 << 0)
#define VRNA_TRANSFORM_LM_OPTION_CLIP_SOURCE_LOW        (1 << 1)
#define VRNA_TRANSFORM_LM_OPTION_CLIP_SOURCE_HIGH       (1 << 2)
#define VRNA_TRANSFORM_LM_OPTION_CLIP_TARGET_LOW        (1 << 3)
#define VRNA_TRANSFORM_LM_OPTION_CLIP_TARGET_HIGH       (1 << 4)
#define VRNA_TRANSFORM_LM_OPTION_CLIP_SOURCE            (VRNA_TRANSFORM_LM_OPTION_CLIP_SOURCE_LOW | VRNA_TRANSFORM_LM_OPTION_CLIP_SOURCE_HIGH)
#define VRNA_TRANSFORM_LM_OPTION_CLIP_TARGET            (VRNA_TRANSFORM_LM_OPTION_CLIP_TARGET_LOW | VRNA_TRANSFORM_LM_OPTION_CLIP_TARGET_HIGH)
#define VRNA_TRANSFORM_LM_OPTION_CLIP                   (VRNA_TRANSFORM_LM_OPTION_CLIP_SOURCE_LOW | VRNA_TRANSFORM_LM_OPTION_CLIP_SOURCE_HIGH | VRNA_TRANSFORM_LM_OPTION_CLIP_TARGET_LOW | VRNA_TRANSFORM_LM_OPTION_CLIP_TARGET_HIGH)
#define VRNA_TRANSFORM_LM_OPTION_DEFAULT                (VRNA_TRANSFORM_LM_OPTION_CLIP)


vrna_data_lin_trans_f
vrna_data_transform_method_bin(double                         (*thresholds)[2],
                               unsigned int                   thresholds_num,
                               double                         oolb_value,
                               double                         ooub_value,
                               unsigned char                  options,
                               vrna_data_lin_trans_opt_t      *transform_options_p,
                               vrna_data_lin_trans_opt_free_f *transform_options_free);


vrna_data_lin_trans_f
vrna_data_transform_method_lm(double                          slope,
                              double                          intercept,
                              double                          domain[4],
                              double                          oob_value,
                              unsigned char                   options,
                              vrna_data_lin_trans_opt_t       *transform_options_p,
                              vrna_data_lin_trans_opt_free_f  *transform_options_free);


vrna_data_lin_trans_f
vrna_data_transform_method_log(double                         value_shift,
                               double                         oob_value,
                               vrna_data_lin_trans_opt_t      *transform_options_p,
                               vrna_data_lin_trans_opt_free_f *transform_options_free);


/**
 *  @}
 */

#endif
