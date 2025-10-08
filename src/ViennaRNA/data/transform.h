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
 *  The transformation function may enforce source and target domain limits if the
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
 *  @see    #vrna_data_lin_trans_f, #vrna_data_lin_trans_opt_t,
 *          #VRNA_TRANSFORM_DEFAULT, #VRNA_TRANSFORM_LM_OPTION_MAP, #VRNA_TRANSFORM_ENFORCE_DOMAINS,
 *          #VRNA_TRANSFORM_ENFORCE_DOMAIN_SOURCE, #VRNA_TRANSFORM_ENFORCE_DOMAIN_TARGET,
 *          #VRNA_TRANSFORM_MAP_SOURCE_LOW, #VRNA_TRANSFORM_MAP_SOURCE_HIGH, #VRNA_TRANSFORM_MAP_SOURCE,
 *          #VRNA_TRANSFORM_MAP_TARGET_LOW, #VRNA_TRANSFORM_MAP_TARGET_HIGH, #VRNA_TRANSFORM_MAP_TARGET,
 *          vrna_data_transform_method_bin(), vrna_data_transform_method_lm(), vrna_data_transform_method_log()
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
vrna_data_lin_transform(const double                *data,
                        size_t                      data_size,
                        vrna_data_lin_trans_f       transform_cb,
                        vrna_data_lin_trans_opt_t   transform_opt,
                        double                      domain[4],
                        double                      oob_value,
                        unsigned int                options);


/**
 *  @brief  Options flag for linear data binning to activate data projection (mapping) instead of actual binning
 *  @see vrna_data_transform_method_bin()
 */
#define VRNA_TRANSFORM_BIN_OPTION_PROJECT               (1 << 0)


/**
 *  @brief  Options flag for linear data binning to indicate that values out-of-upper-bound are to be mapped to the respective domain limit
 *  @see vrna_data_transform_method_bin()
 */
#define VRNA_TRANSFORM_BIN_OPTION_MAP_OUTOF_UPPERBOUND  (1 << 1)


/**
 *  @brief  Options flag for linear data binning to indicate that values out-of-lower-bound are to be mapped to the respective domain limit
 *  @see vrna_data_transform_method_bin()
 */
#define VRNA_TRANSFORM_BIN_OPTION_MAP_OUTOF_LOWERBOUND  (1 << 2)


/**
 *  @brief  Options flag for linear data binning to indicate default settings
 *  @see vrna_data_transform_method_bin()
 */
#define VRNA_TRANSFORM_BIN_OPTION_DEFAULT               0

/**
 *  @brief  Retrieve a linear data transform callback that performs (discrete) binning and more
 *
 *  This function yields a linear data transform callback and the associated transform options
 *  data structure suitable for usage in vrna_data_lin_transform(). The transformation this
 *  callback performs can be described as (discrete) binning or data bucketing.
 *
 *  In more detail, the transformation callback can be used to map continuous data from
 *  one domain (source) into discrete or continuous data of another domain (target).
 *  User-defined domain boundaries are provided by the @p thresholds argument, which boils
 *  down to a list of pairs of values, where the first value is the boundary in the source
 *  domain. The second value is the corresponding boundary of the target domain. The first and
 *  last pair of domain boundaries denote the lower and the upper domain limits, respectively.
 *  Any value outside of these limits can either be discarded by assigning them a special
 *  out-of-bounds value (@p oolb_value, @p ooub_value), or they can be mapped directly to the
 *  domain limits. Control over this mapping is available through the @p options argument
 *  by using the binary flags #VRNA_TRANSFORM_BIN_OPTION_MAP_OUTOF_UPPERBOUND and
 *  #VRNA_TRANSFORM_BIN_OPTION_MAP_OUTOF_LOWERBOUND.
 *
 *  By default, the transformation maps any data within a source interval as specified
 *  by two consecutive entries in the @p thresholds argument to the exact value of the
 *  second target boundary. This process is also called @c binning. Alternatively, this
 *  implementation also allows for a continuous mapping (projection) of the source intervals
 *  into the target intervals. To activate this behavior, the #VRNA_TRANSFORM_BIN_OPTION_PROJECT
 *  flag must be provided to the @p options argument.
 *
 *  The @p transform_options_p and @p transform_options_free pointers are used as additional
 *  output to obtain the addresses of the transformation option data structure that has to
 *  be provided to the vrna_data_lin_transform() function and a function pointer to release
 *  the memory of the option data structure once it is not required anymore.
 *  
 *  @see  vrna_data_lin_transform(),
 *        #VRNA_TRANSFORM_BIN_OPTION_DEFAULT, #VRNA_TRANSFORM_BIN_OPTION_PROJECT,
 *        #VRNA_TRANSFORM_BIN_OPTION_MAP_OUTOF_UPPERBOUND, #VRNA_TRANSFORM_BIN_OPTION_MAP_OUTOF_LOWERBOUND,
 *        vrna_data_transform_method_lm(), vrna_data_transform_method_log()
 *
 *  @param  thresholds              A pointer to an array of data pairs holding the source and target domain boundaries
 *  @param  thresholds_num          The number of domain boundary pairs available in @p thresholds
 *  @param  oolb_value              Out-of-lower-bound value
 *  @param  ooub_value              Out-of-upper-bound value
 *  @param  options                 Additional options that change the behavior of the callback function
 *  @param  transform_options_p     A pointer to store the address of the options data structure
 *  @param  transform_options_free  A pointer to store the address of the @c free function that releases the memory of the options data structure
 *  @return                         A callback function that performs (discrete) data binning
 */
vrna_data_lin_trans_f
vrna_data_transform_method_bin(double                         (*thresholds)[2],
                               unsigned int                   thresholds_num,
                               double                         oolb_value,
                               double                         ooub_value,
                               unsigned int                   options,
                               vrna_data_lin_trans_opt_t      *transform_options_p,
                               vrna_data_lin_trans_opt_free_f *transform_options_free);



/**
 *  @brief  Options flag for transforming linear data using a linear model that enables log-transform of the source value
 *  @see vrna_data_transform_method_lm()
 */
#define VRNA_TRANSFORM_LM_OPTION_LOG                    (1 << 0)

#define VRNA_TRANSFORM_LM_OPTION_DEFAULT                0U

/**
 *  @brief  Retrieve a linear data transform callback that applies a linear model
 *
 *  This function yields a linear data transform callback and the associated transform options
 *  data structure suitable for usage in vrna_data_lin_transform(). The callback applies
 *  a linear model of the form
 *  @f[ y = a + b \cdot f(x) @f]
 *  where @f$ x @f$ is the input value (@c source), @f$ a @f$ and @f$ b @f$ are intercept
 *  and slope, and @f$ y @f$ is the output value (@c target). By default, the function
 *  @f$ f(x) = x @f$ is the identity function. However, the callback can be instructed
 *  to apply a log-transform, @f$ f(x) = log x @f$, instead. Control over the behavior
 *  of @f$ f(x) @f$ can be gained by providing a corresponding option flag to the
 *  @p options argument, e.g. #VRNA_TRANSFORM_LM_OPTION_LOG.
 *
 *  The @p transform_options_p and @p transform_options_free pointers are used as additional
 *  output to obtain the addresses of the transformation option data structure that has to
 *  be provided to the vrna_data_lin_transform() function and a function pointer to release
 *  the memory of the option data structure once it is not required anymore.
 *  
 *  @see  vrna_data_lin_transform(), #vrna_data_lin_trans_f, #vrna_data_lin_trans_opt_t, #vrna_data_lin_trans_opt_free_f,
 *        #VRNA_TRANSFORM_LM_OPTION_DEFAULT, #VRNA_TRANSFORM_LM_OPTION_LOG,
 *        vrna_data_transform_method_bin(), vrna_data_transform_method_log()
 *
 *  @param  slope                   The slope of the linear function
 *  @param  intercept               The intercept of the linear function
 *  @param  options                 Additional options that change the behavior of the callback function
 *  @param  transform_options_p     A pointer to store the address of the options data structure
 *  @param  transform_options_free  A pointer to store the address of the @c free function that releases the memory of the options data structure
 *  @return                         A callback function that performs transformation through a linear model
 */
vrna_data_lin_trans_f
vrna_data_transform_method_lm(double                          slope,
                              double                          intercept,
                              unsigned int                    options,
                              vrna_data_lin_trans_opt_t       *transform_options_p,
                              vrna_data_lin_trans_opt_free_f  *transform_options_free);


/**
 *  @brief  Options flag for transforming linear data using log transform to use a non-default base
 *  @see vrna_data_transform_method_log()
 */
#define VRNA_TRANSFORM_LOG_OPTION_NONDEFAULT_BASE        (1 << 0)


/**
 *  @brief  Options flag for transforming linear data using log transform that indicates default settings
 *  @see vrna_data_transform_method_log()
 */
#define VRNA_TRANSFORM_LOG_OPTION_DEFAULT                0


/**
 *  @brief  Retrieve a linear data transform callback that applies a log transformation
 *
 *  This function yields a linear data transform callback and the associated transform options
 *  data structure suitable for usage in vrna_data_lin_transform(). The callback applies
 *  a log transform of the form
 *  @f[ y = \log_b (x + c) @f]
 *  where @f$ x @f$ is the input value (@c source), @f$ b @f$ is the base, @f$ c @f$
 *  is a value to allow for shifting the source data, and @f$ y @f$ is the output value
 *  (@c target). By default, the natural logarithm is used, i.e. @f$ b = \mathrm{e} @f$. The
 *  @p base argument in combination with the #VRNA_TRANSFORM_LOG_OPTION_NONDEFAULT_BASE
 *  flag supplied to the @p options argument can be used to change the base to any other
 *  number.
 *
 *  The @p value_shift argument corresponds to @f$ c @f$ in the above formula and in
 *  most cases should be @c 0. The callback returns the out-of-bounds value @p oob_value
 *  if @f$ x + c \le 0 @f$.
 *
 *  The @p transform_options_p and @p transform_options_free pointers are used as additional
 *  output to obtain the addresses of the transformation option data structure that has to
 *  be provided to the vrna_data_lin_transform() function and a function pointer to release
 *  the memory of the option data structure once it is not required anymore.
 *  
 *  @see  vrna_data_lin_transform(), #vrna_data_lin_trans_f, #vrna_data_lin_trans_opt_t, #vrna_data_lin_trans_opt_free_f,
 *        #VRNA_TRANSFORM_LOG_OPTION_DEFAULT, #VRNA_TRANSFORM_LOG_OPTION_NONDEFAULT_BASE,
 *        vrna_data_transform_method_bin(), vrna_data_transform_method_lm()
 *
 *  @param  value_shift             The shift value @f$ c @f$
 *  @param  base                    The base @f$ b @f$ of the logarithm (only used if #VRNA_TRANSFORM_LOG_OPTION_NONDEFAULT_BASE is passed to @p options
 *  @param  oob_value               Out-of-bound value
 *  @param  options                 Additional options that change the behavior of the callback function
 *  @param  transform_options_p     A pointer to store the address of the options data structure
 *  @param  transform_options_free  A pointer to store the address of the @c free function that releases the memory of the options data structure
 *  @return                         A callback function that performs transformation through a linear model
 */
vrna_data_lin_trans_f
vrna_data_transform_method_log(double                         value_shift,
                               double                         base,
                               double                         oob_value,
                               unsigned int                   options,
                               vrna_data_lin_trans_opt_t      *transform_options_p,
                               vrna_data_lin_trans_opt_free_f *transform_options_free);


/**
 *  @brief  Options flag for transforming linear data using logistic function that indicates default settings
 *  @see vrna_data_transform_method_loogistic()
 */
#define VRNA_TRANSFORM_LOGISTIC_OPTION_DEFAULT         0U


/**
 *  @brief  Retrieve a linear data transform callback that applies a logistic function
 *
 *  This function yields a linear data transform callback and the associated transform options
 *  data structure suitable for usage in vrna_data_lin_transform(). The callback applies
 *  a logistic function of the form
 *  @f[ y = \frac{L}{1 + e^{-k \cdot (x - x_0)}}  @f]
 *  where @f$ x @f$ is the input value (@c source), @f$ x_0 @f$ is the mid point (@p mid_point)
 *  of the function, @f$ k @f$ is the logistic growth rate (@p growth_rate), and @f$ L @f$
 *  is the supremum (or carrying capacity) of the function (@p supremum). The standard logistic
 *  function is defined as @f$ L = 1 @f$, @f$ k = 1 @f$. and @f$ x_0 = 0 @f$, i.e.
 *  @f[ y = \frac{1}{1 + e^{-x}} @f]
 *
 *  The @p transform_options_p and @p transform_options_free pointers are used as additional
 *  output to obtain the addresses of the transformation option data structure that has to
 *  be provided to the vrna_data_lin_transform() function and a function pointer to release
 *  the memory of the option data structure once it is not required anymore.
 *  
 *  @see  vrna_data_lin_transform(), #vrna_data_lin_trans_f, #vrna_data_lin_trans_opt_t, #vrna_data_lin_trans_opt_free_f,
 *        #VRNA_TRANSFORM_LOGISTIC_OPTION_DEFAULT,
 *        vrna_data_transform_method_bin(), vrna_data_transform_method_log(),
 *        vrna_data_transform_method_lm()
 *
 *  @param  mid_point               The midpoint of the function (default = 0)
 *  @param  supremum                The supremum of the function (default = 1)
 *  @param  growth_rate             The growth rate of the function (default = 1)
 *  @param  options                 Additional options that change the behavior of the callback function
 *  @param  transform_options_p     A pointer to store the address of the options data structure
 *  @param  transform_options_free  A pointer to store the address of the @c free function that releases the memory of the options data structure
 *  @return                         A callback function that performs transformation through a logistic function
 */
vrna_data_lin_trans_f
vrna_data_transform_method_logistic(double                         mid_point,
                                    double                         supremum,
                                    double                         growth_rate,
                                    unsigned int                   options,
                                    vrna_data_lin_trans_opt_t      *transform_options_p,
                                    vrna_data_lin_trans_opt_free_f *transform_options_free);


/**
 *  @brief  Options flag for transforming linear data using a kernel density estimate (KDE) that indicates default settings
 *  @see vrna_data_transform_method_lm()
 */
#define VRNA_TRANSFORM_KDE_OPTION_DEFAULT   0U


/**
 *  @brief  Retrieve a linear data transform callback that applies a kernel density estimation (KDE)
 *
 *  This function yields a linear data transform callback and the associated transform options
 *  data structure suitable for usage in vrna_data_lin_transform(). The callback applies
 *  a kernel density estimation of the form
 *  @f[ y = \frac{1}{n \cdot h} \sum_{i = 1}^{n} K \big( \frac{x - x_i}{h} \big) @f]
 *  over a set of @f$ n > 0 @f$ samples @f$ \{ x_1, x_2, \ldots, x_n \} @f$ (@p samples), where
 *  @f$ K @f$ (@p kernel) is a non-negative kernel function, @f$ h > 0 @f$ is the bandwidth
 *  (@p bandwidth), and @f$ x @f$ is the  input value (@c source).
 *
 *  @note   If no kernel is provided, i.e. @c NULL is passed to @p kernel, the KDE uses a
 *          standard Gaussian kernel by providing vrna_data_transform_method_Gaussian()
 *          with @f$ a = \frac{1}{\sqrt{2 \pi}}, b = 0, c = 1 @f$. In this case, the
 *          bandwidth (@p bandwidth) will be automatically replaced by a rule-of-thumb bandwidth
 *          estimator @f$ h = ( \frac{4}{3n} )^{\frac{1}{5}} @f$
 *  @note   Whenever the @p transform_options_free callback will be called, this function
 *          also releases the memory occupied by @p kernel_data through the @p kernel_data_free
 *          function if it is provided, i.e. not @c NULL.
 *
 *  The @p transform_options_p and @p transform_options_free pointers are used as additional
 *  output to obtain the addresses of the transformation option data structure that has to
 *  be provided to the vrna_data_lin_transform() function and a function pointer to release
 *  the memory of the option data structure once it is not required anymore.
 *  
 *  @see  vrna_data_lin_transform(), #vrna_data_lin_trans_f, #vrna_data_lin_trans_opt_t, #vrna_data_lin_trans_opt_free_f,
 *        #VRNA_TRANSFORM_KDE_OPTION_DEFAULT, 
 *        vrna_data_transform_method_gaussian()
 *
 *  @param  samples                 The samples of the distribution to estimate @f$ \{ x_1, x_2, \ldots, x_n \} @f$
 *  @param  num_samples             The number of samples (@f$ n @f$) provided by pointer @p samples
 *  @param  kernel                  The kernel function @f$ K @f$ (maybe @c NULL)
 *  @param  kernel_data             The data that should be passed-through to the kernel function (maybe @c NULL)
 *  @param  kernel_data_free        A function to free the memory occupied by the data for the kernel function (maybe @c NULL)
 *  @param  bandwidth               The bandwidth @f$ h @f$
 *  @param  options                 Additional options that change the behavior of the callback function
 *  @param  transform_options_p     A pointer to store the address of the options data structure
 *  @param  transform_options_free  A pointer to store the address of the @c free function that releases the memory of the options data structure
 *  @return                         A callback function that performs transformation through kernel density estimation
 */
vrna_data_lin_trans_f
vrna_data_transform_method_kde(double                         *samples,
                               size_t                         num_samples,
                               vrna_data_lin_trans_f          kernel,
                               vrna_data_lin_trans_opt_t      kernel_data,
                               vrna_data_lin_trans_opt_free_f kernel_data_free,
                               double                         bandwidth,
                               unsigned int                   options,
                               vrna_data_lin_trans_opt_t      *transform_options_p,
                               vrna_data_lin_trans_opt_free_f *transform_options_free);


/**
 *  @brief  Options flag for transforming linear data using a Gaussian function that indicates default settings
 *  @see vrna_data_transform_method_gaussian()
 */
#define VRNA_TRANSFORM_GAUSSIAN_OPTION_DEFAULT  0U


/**
 *  @brief  Retrieve a linear data transform callback that applies Gaussian function
 *
 *  This function yields a linear data transform callback and the associated transform options
 *  data structure suitable for usage in vrna_data_lin_transform(). The callback applies
 *  a Gaussian function of the form
 *  @f[ y = a \cdot \exp(- \frac{(x - b)^2}{2 c^2}) @f]
 *  with arbitrary constants @f$ a, b @f$ and non-zero @f$ c @f$ and @f$ x @f$ being the
 *  input value (@c source).
 *
 *  The @p transform_options_p and @p transform_options_free pointers are used as additional
 *  output to obtain the addresses of the transformation option data structure that has to
 *  be provided to the vrna_data_lin_transform() function and a function pointer to release
 *  the memory of the option data structure once it is not required anymore.
 *  
 *  @see  vrna_data_lin_transform(), #vrna_data_lin_trans_f, #vrna_data_lin_trans_opt_t, #vrna_data_lin_trans_opt_free_f,
 *        #VRNA_TRANSFORM_GAUSSIAN_OPTION_DEFAULT
 *
 *  @param  a                       The arbitrary constant @f$ a @f$
 *  @param  b                       The arbitrary constant @f$ b @f$
 *  @param  c                       The arbitrary non-zero constant @f$ c @f$
 *  @param  options                 Additional options that change the behavior of the callback function
 *  @param  transform_options_p     A pointer to store the address of the options data structure
 *  @param  transform_options_free  A pointer to store the address of the @c free function that releases the memory of the options data structure
 *  @return                         A callback function that performs transformation through a Gaussian function
 */
vrna_data_lin_trans_f
vrna_data_transform_method_gaussian(double                          a,
                                    double                          b,
                                    double                          c,
                                    unsigned int                    options,
                                    vrna_data_lin_trans_opt_t       *transform_options_p,
                                    vrna_data_lin_trans_opt_free_f  *transform_options_free);


/**
 *  @}
 */

#endif
