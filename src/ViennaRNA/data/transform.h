#ifndef VIENNA_RNA_PACKAGE_DATA_TRANSFORM_H
#define VIENNA_RNA_PACKAGE_DATA_TRANSFORM_H

#include <ViennaRNA/fold_compound.h>
#include <ViennaRNA/math/functions.h>

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
 *  @brief  Options flag for transforming linear data to enforce source domain limits
 *  @see vrna_data_lin_transform()
 */
#define VRNA_TRANSFORM_ENFORCE_DOMAIN_SOURCE  (1 << 1)


/**
 *  @brief  Options flag for transforming linear data to enforce target domain limits
 *  @see vrna_data_lin_transform()
 */
#define VRNA_TRANSFORM_ENFORCE_DOMAIN_TARGET  (1 << 2)


/**
 *  @brief  Options flag for transforming linear data to enforce source and target domain limits
 *  @see vrna_data_lin_transform()
 */
#define VRNA_TRANSFORM_ENFORCE_DOMAINS        (VRNA_TRANSFORM_ENFORCE_DOMAIN_SOURCE | VRNA_TRANSFORM_ENFORCE_DOMAIN_TARGET)


/**
 *  @brief  Options flag for transforming linear data to map source values below the domain limit to the lower domain limit
 *  @see vrna_data_lin_transform()
 */
#define VRNA_TRANSFORM_MAP_SOURCE_LOW         (1 << 3)


/**
 *  @brief  Options flag for transforming linear data to map source values above the domain limit to the upper domain limit
 *  @see vrna_data_lin_transform()
 */
#define VRNA_TRANSFORM_MAP_SOURCE_HIGH        (1 << 4)


/**
 *  @brief  Options flag for transforming linear data to map source values below and above the domain limits to the domain limits
 *  @see vrna_data_lin_transform()
 */
#define VRNA_TRANSFORM_MAP_SOURCE             (VRNA_TRANSFORM_MAP_SOURCE_LOW | VRNA_TRANSFORM_MAP_SOURCE_HIGH)


/**
 *  @brief  Options flag for transforming linear data to map target values below the domain limit to the lower domain limit
 *  @see vrna_data_lin_transform()
 */
#define VRNA_TRANSFORM_MAP_TARGET_LOW         (1 << 5)


/**
 *  @brief  Options flag for transforming linear data to map target values above the domain limit to the upper domain limit
 *  @see vrna_data_lin_transform()
 */
#define VRNA_TRANSFORM_MAP_TARGET_HIGH        (1 << 6)


/**
 *  @brief  Options flag for transforming linear data to map target values below and above the domain limits to the domain limits
 *  @see vrna_data_lin_transform()
 */
#define VRNA_TRANSFORM_MAP_TARGET             (VRNA_TRANSFORM_MAP_TARGET_LOW | VRNA_TRANSFORM_MAP_TARGET_HIGH)


/**
 *  @brief  Options flag for transforming linear data to map source and target values below and above the domain limits to the respective domain limits
 *  @see vrna_data_lin_transform()
 */
#define VRNA_TRANSFORM_MAP                    (VRNA_TRANSFORM_MAP_SOURCE_LOW | VRNA_TRANSFORM_MAP_SOURCE_HIGH | VRNA_TRANSFORM_MAP_TARGET_LOW | VRNA_TRANSFORM_MAP_TARGET_HIGH)


/**
 *  @brief  Options flag for transforming linear data that indicates default settings
 *  @see vrna_data_lin_transform()
 */
#define VRNA_TRANSFORM_DEFAULT                (VRNA_TRANSFORM_ENFORCE_DOMAINS | VRNA_TRANSFORM_MAP)


/**
 *  @brief  Transform an array of linear data
 *
 *  This function transforms an array of linear data (@c source) into another array of
 *  linear data of the same size (@c target). For that purpose, it utilizes a callback
 *  mechanism that transforms each single value.
 *
 *  During data transform, the callback function (@p transform_cb) will be executed for
 *  each value of the linear data (@p data) and is responsible to actually transform the
 *  data value. The @p transform_opt parameter will be simply passed-through to the callback
 *  function as second argument. The callback function has to return the transformed
 *  value.
 *
 *  This transformation function may enforce source and target domain limits if the
 *  corresponding flag (#VRNA_TRANSFORM_ENFORCE_DOMAINS) is provided to the @p options
 *  argument. This means that one has complete control over the accepted values of the
 *  source and target domains. The limits can be provided to this function through the
 *  @p domain argument. Here, the order of the limits follows:
 *  @f[ x_\text{min}, x_\text{max}, y_\text{min}, y_\text{max}. @f]
 *
 *  If @c NULL is passed instead of an actual domain array, domain enforcement is
 *  deactivated and any value will pass. In case a domain enforcing is active and a
 *  source or target value doesn't meet the respective limits, the function assigns
 *  it the out-of-bounds value @p oob_value. This behavior can be changed to a mapping
 *  of out-of-bounds values to the respetive domain limits. For that purpose, the
 *  @p options argument requires the flag #VRNA_TRANSFORM_MAP.
 *
 *  @note   Individual control for mapping out-of-bounds values to the four domain
 *          limits, i.e. low and high values of source and target can be gained by
 *          providing the @p options argument the #VRNA_TRANSFORM_MAP_SOURCE_LOW,
 *          #VRNA_TRANSFORM_MAP_SOURCE_HIGH, #VRNA_TRANSFORM_MAP_TARGET_LOW,
 *          and #VRNA_TRANSFORM_MAP_TARGET_HIGH flags.
 *  @note   When @p data is @c NULL or @p data_size equals @c 0, the function returns
 *          @c NULL. Moreover, the function simply provides a copy of the linear data
 *          if the transformation callback @p transform_cb is not provided, i.e. if
 *          it is @c NULL. In this case, domain limits will still be enforced if the
 *          corresponding options are set.
 *
 *  @see    #vrna_math_fun_dbl_f, #vrna_math_fun_dbl_opt_t,
 *          #VRNA_TRANSFORM_DEFAULT,
 *          #VRNA_TRANSFORM_ENFORCE_DOMAINS,
 *          #VRNA_TRANSFORM_ENFORCE_DOMAIN_SOURCE, #VRNA_TRANSFORM_ENFORCE_DOMAIN_TARGET,
 *          #VRNA_TRANSFORM_MAP, #VRNA_TRANSFORM_MAP_SOURCE, #VRNA_TRANSFORM_MAP_TARGET,
 *          #VRNA_TRANSFORM_MAP_SOURCE_LOW, #VRNA_TRANSFORM_MAP_SOURCE_HIGH,
 *          #VRNA_TRANSFORM_MAP_TARGET_LOW, #VRNA_TRANSFORM_MAP_TARGET_HIGH,
 *          vrna_math_fun_dbl_bin_opt(), vrna_math_fun_dbl_linear_opt(), vrna_math_fun_dbl_log_opt(),
 *          vrna_math_fun_dbl_logistic_opt(), vrna_math_fun_dbl_gaussian_opt(), vrna_math_fun_dbl_kde_opt()
 *
 *  @param  data            A pointer to an array of linear data
 *  @param  data_size       The size of the array @p data is pointing to
 *  @param  transform_cb    The data transformation callback (maybe @c NULL)
 *  @param  transform_opt   The options that need to be passed through to the transformation callback @p transform_cb (maybe @c NULL)
 *  @param  domain          The domain limits (maybe @c NULL)
 *  @param  oob_value       Out-of-bound value
 *  @param  options         Additional options that change the behavior of the transformation function
 *  @return                 A pointer to an array of transformed linear data (or @c NULL on any error)
 */
double *
vrna_data_lin_transform(const double        *data,
                        size_t              data_size,
                        vrna_math_fun_dbl_f     transform_cb,
                        vrna_math_fun_dbl_opt_t transform_opt,
                        double              domain[4],
                        double              oob_value,
                        unsigned int        options);


/**
 *  @}
 */

#endif
