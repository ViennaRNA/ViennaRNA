#ifndef VIENNA_RNA_PACKAGE_MATH_FUNCTIONS_H
#define VIENNA_RNA_PACKAGE_MATH_FUNCTIONS_H


/**
 *  @file     ViennaRNA/math/functions.h
 *  @ingroup  math
 *  @brief    This module provides mathematical functions
 */

/**
 *  @addtogroup math Mathematical Functions
 *  @{
 */


/**
 *  @brief  Option data structure type for mathematical functions
 */
typedef   void                  *vrna_math_fun_opt_t;

/**
 *  @brief  Callback function to release any memory occupied by the mathematical function option #vrna_math_fun_opt_t
 */
typedef   vrna_auxdata_free_f   vrna_math_fun_opt_free_f;


/**
 *  @brief Callback function prototype for scalar mathematical functions of the form @f$ y = f(x) @f$
 *
 *  This is the prototype for callback functions used for (scalar) mathematical functions. Its
 *  main purpose is to encapsulate the options (arguments, such as constansts, etc.) for the
 *  mathematical function and can therefore be used in a more generalized way wherever different
 *  mathematical functions may be applied in the same manner. Functions of this type have to be
 *  provided two parameters: @p value which is the value @f$ x @f$ that is about to be evaluated
 *  @p options, which is a pointer to a data structure that holds any means required to ensure
 *  that the callback can evaluate the provide data. The callback function then returns the
 *  transformed value @f$ y @f$.
 *
 *  @callback
 *  @parblock
 *  This callback handles mathematical functions working on scalar values.
 *  @endparblock
 *
 *  @see  vrna_data_lin_transform(),
 *        vrna_math_fun_bin(), vrna_math_fun_linear(), vrna_math_fun_log(), vrna_math_fun_gaussian(),
 *        vrna_math_fun_logistic(), vrna_math_fun_kde()
 *
 * @param value     The value @f$ x @f$ that is to be evaluated by the function @f$ f(x) @f$
 * @param options   A pointer to a data structure specific for the mathematical function
 * @return          The evaluated value @f$ y @f$
 */
typedef double (*vrna_math_fun_f) (double               value,
                                   vrna_math_fun_opt_t  options);


/**
 *  @brief  Options flag for linear data binning to activate data projection (mapping) instead of actual binning
 *  @see vrna_math_fun_bin(), vrna_math_fun_bin_opt()
 */
#define VRNA_MATH_FUN_BIN_OPTION_PROJECT              (1 << 0)


/**
 *  @brief  Options flag for linear data binning to indicate that values out-of-upper-bound are to be mapped to the respective domain limit
 *  @see vrna_math_fun_bin(), vrna_math_fun_bin_opt()
 */
#define VRNA_MATH_FUN_BIN_OPTION_MAP_OUTOF_UPPERBOUND (1 << 1)


/**
 *  @brief  Options flag for linear data binning to indicate that values out-of-lower-bound are to be mapped to the respective domain limit
 *  @see vrna_math_fun_bin(), vrna_math_fun_bin_opt()
 */
#define VRNA_MATH_FUN_BIN_OPTION_MAP_OUTOF_LOWERBOUND (1 << 2)


/**
 *  @brief  Options flag for linear data binning to indicate default settings
 *  @see vrna_math_fun_bin(), vrna_math_fun_bin_opt()
 */
#define VRNA_MATH_FUN_BIN_OPTION_DEFAULT              0U


double
vrna_math_fun_bin(double              v,
                  double              (*thresholds)[2],
                  size_t              num_thresholds,
                  double              oolb_value,
                  double              ooub_value,
                  unsigned int        options);


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
 *  by using the binary flags #VRNA_MATH_FUN_BIN_OPTION_MAP_OUTOF_UPPERBOUND and
 *  #VRNA_MATH_FUN_BIN_OPTION_MAP_OUTOF_LOWERBOUND.
 *
 *  By default, the transformation maps any data within a source interval as specified
 *  by two consecutive entries in the @p thresholds argument to the exact value of the
 *  second target boundary. This process is also called @c binning. Alternatively, this
 *  implementation also allows for a continuous mapping (projection) of the source intervals
 *  into the target intervals. To activate this behavior, the #VRNA_MATH_FUN_BIN_OPTION_PROJECT
 *  flag must be provided to the @p options argument.
 *
 *  The @p transform_options_p and @p transform_options_free pointers are used as additional
 *  output to obtain the addresses of the transformation option data structure that has to
 *  be provided to the vrna_data_lin_transform() function and a function pointer to release
 *  the memory of the option data structure once it is not required anymore.
 *  
 *  @see  vrna_data_lin_transform(),
 *        #VRNA_MATH_FUN_BIN_OPTION_DEFAULT, #VRNA_MATH_FUN_BIN_OPTION_PROJECT,
 *        #VRNA_MATH_FUN_BIN_OPTION_MAP_OUTOF_UPPERBOUND, #VRNA_MATH_FUN_BIN_OPTION_MAP_OUTOF_LOWERBOUND,
 *        vrna_math_fun_linear(), vrna_math_fun_log(), vrna_math_fun_kde()
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
vrna_math_fun_f
vrna_math_fun_bin_opt(double                    (*thresholds)[2],
                      unsigned int              thresholds_num,
                      double                    oolb_value,
                      double                    ooub_value,
                      unsigned int              options,
                      vrna_math_fun_opt_t       *transform_options_p,
                      vrna_math_fun_opt_free_f  *transform_options_free);



/**
 *  @brief  Options flag for transforming linear data using a linear function that enables log-transform of the source value
 *  @see vrna_math_fun_linear()
 */
#define VRNA_MATH_FUN_LINEAR_OPTION_LOG         (1 << 0)

/**
 *  @brief  Options flag for transforming linear data using linear function that indicates default settings
 *  @see vrna_math_fun_linear()
 */
#define VRNA_MATH_FUN_LINEAR_OPTION_DEFAULT     0U


double
vrna_math_fun_linear(double       v,
                     double       slope,
                     double       intercept,
                     unsigned int options);


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
 *  @p options argument, e.g. #VRNA_MATH_FUN_LINEAR_OPTION_LOG.
 *
 *  The @p transform_options_p and @p transform_options_free pointers are used as additional
 *  output to obtain the addresses of the transformation option data structure that has to
 *  be provided to the vrna_data_lin_transform() function and a function pointer to release
 *  the memory of the option data structure once it is not required anymore.
 *  
 *  @see  vrna_data_lin_transform(), #vrna_math_fun_f, #vrna_math_fun_opt_t, #vrna_math_fun_opt_free_f,
 *        #VRNA_MATH_FUN_LINEAR_OPTION_DEFAULT, #VRNA_MATH_FUN_LINEAR_OPTION_LOG,
 *        vrna_math_fun_bin(), vrna_math_fun_log()
 *
 *  @param  slope                   The slope of the linear function
 *  @param  intercept               The intercept of the linear function
 *  @param  options                 Additional options that change the behavior of the callback function
 *  @param  transform_options_p     A pointer to store the address of the options data structure
 *  @param  transform_options_free  A pointer to store the address of the @c free function that releases the memory of the options data structure
 *  @return                         A callback function that performs transformation through a linear model
 */
vrna_math_fun_f
vrna_math_fun_linear_opt(double                   slope,
                         double                   intercept,
                         unsigned int             options,
                         vrna_math_fun_opt_t      *transform_options_p,
                         vrna_math_fun_opt_free_f *transform_options_free);


/**
 *  @brief  Options flag for transforming linear data using log transform to use a non-default base
 *  @see vrna_math_fun_log()
 */
#define VRNA_MATH_FUN_LOG_OPTION_NONDEFAULT_BASE    (1 << 0)


/**
 *  @brief  Options flag for transforming linear data using log transform that indicates default settings
 *  @see vrna_math_fun_log()
 */
#define VRNA_MATH_FUN_LOG_OPTION_DEFAULT            0


double
vrna_math_fun_log(double        value,
                  double        base,
                  double        oob_value,
                  unsigned int  options);


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
 *  @p base argument in combination with the #VRNA_MATH_FUN_LOG_OPTION_NONDEFAULT_BASE
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
 *  @see  vrna_data_lin_transform(), #vrna_math_fun_f, #vrna_math_fun_opt_t, #vrna_math_fun_opt_free_f,
 *        #VRNA_MATH_FUN_LOG_OPTION_DEFAULT, #VRNA_MATH_FUN_LOG_OPTION_NONDEFAULT_BASE,
 *        vrna_math_fun_bin(), vrna_math_fun_linear()
 *
 *  @param  value_shift             The shift value @f$ c @f$
 *  @param  base                    The base @f$ b @f$ of the logarithm (only used if #VRNA_MATH_FUN_LOG_OPTION_NONDEFAULT_BASE is passed to @p options
 *  @param  oob_value               Out-of-bound value
 *  @param  options                 Additional options that change the behavior of the callback function
 *  @param  transform_options_p     A pointer to store the address of the options data structure
 *  @param  transform_options_free  A pointer to store the address of the @c free function that releases the memory of the options data structure
 *  @return                         A callback function that performs transformation through a linear model
 */
vrna_math_fun_f
vrna_math_fun_log_opt(double                    value_shift,
                      double                    base,
                      double                    oob_value,
                      unsigned int              options,
                      vrna_math_fun_opt_t       *transform_options_p,
                      vrna_math_fun_opt_free_f  *transform_options_free);


/**
 *  @brief  Options flag for transforming linear data using logistic function that indicates default settings
 *  @see vrna_math_fun_loogistic()
 */
#define VRNA_MATH_FUN_LOGISTIC_OPTION_DEFAULT     0U


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
 *  @see  vrna_data_lin_transform(), #vrna_math_fun_f, #vrna_math_fun_opt_t, #vrna_math_fun_opt_free_f,
 *        #VRNA_MATH_FUN_LOGISTIC_OPTION_DEFAULT,
 *        vrna_math_fun_bin(), vrna_math_fun_log(),
 *        vrna_math_fun_linear()
 *
 *  @param  mid_point               The midpoint of the function (default = 0)
 *  @param  supremum                The supremum of the function (default = 1)
 *  @param  growth_rate             The growth rate of the function (default = 1)
 *  @param  options                 Additional options that change the behavior of the callback function
 *  @param  transform_options_p     A pointer to store the address of the options data structure
 *  @param  transform_options_free  A pointer to store the address of the @c free function that releases the memory of the options data structure
 *  @return                         A callback function that performs transformation through a logistic function
 */
vrna_math_fun_f
vrna_math_fun_logistic_opt(double                   mid_point,
                           double                   supremum,
                           double                   growth_rate,
                           unsigned int             options,
                           vrna_math_fun_opt_t      *transform_options_p,
                           vrna_math_fun_opt_free_f *transform_options_free);


/**
 *  @brief  Options flag for transforming linear data using a kernel density estimate (KDE) that indicates default settings
 *  @see vrna_math_fun_linear()
 */
#define VRNA_MATH_FUN_KDE_OPTION_DEFAULT      0U


double
vrna_math_fun_kde(double                    value,
                  double                    *samples,
                  size_t                    num_samples,
                  vrna_math_fun_f           kernel,
                  vrna_math_fun_opt_t       kernel_data,
                  double                    bandwidth,
                  unsigned int              options);


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
 *          standard Gaussian kernel by providing vrna_math_fun_Gaussian()
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
 *  @see  vrna_data_lin_transform(), #vrna_math_fun_f, #vrna_math_fun_opt_t, #vrna_math_fun_opt_free_f,
 *        #VRNA_MATH_FUN_KDE_OPTION_DEFAULT, 
 *        vrna_math_fun_gaussian()
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
vrna_math_fun_f
vrna_math_fun_kde_opt(double                    *samples,
                      size_t                    num_samples,
                      vrna_math_fun_f           kernel,
                      vrna_math_fun_opt_t       kernel_data,
                      vrna_math_fun_opt_free_f  kernel_data_free,
                      double                    bandwidth,
                      unsigned int              options,
                      vrna_math_fun_opt_t       *transform_options_p,
                      vrna_math_fun_opt_free_f  *transform_options_free);


/**
 *  @brief  Options flag for transforming linear data using a Gaussian function that indicates default settings
 *  @see vrna_math_fun_gaussian()
 */
#define VRNA_MATH_FUN_GAUSSIAN_OPTION_DEFAULT     0U


/**
 *  @brief  Evaluate a Gaussian function for a scalar input value
 *
 *  This function applies a Gaussian function of the form
 *  @f[ y = f(x) = a \cdot \exp(- \frac{(x - b)^2}{2 c^2}) @f]
 *  with arbitrary constants @f$ a, b @f$ and non-zero @f$ c @f$ to a scalar input
 *  value @f$ x @f$.
 *
 *  @see  vrna_math_fun_gaussian_opt()
 *
 *  @param  a                       The arbitrary constant @f$ a @f$
 *  @param  a                       The arbitrary constant @f$ a @f$
 *  @param  b                       The arbitrary constant @f$ b @f$
 *  @param  c                       The arbitrary non-zero constant @f$ c @f$
 *  @return                         The evaluation for the Gaussian function @f$ y = f(x) @f$
 */
double
vrna_math_fun_gaussian(double x,
                       double a,
                       double b,
                       double c);

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
 *  @see  vrna_math_fun_gaussian(),
 *        vrna_data_lin_transform(), #vrna_math_fun_f, #vrna_math_fun_opt_t, #vrna_math_fun_opt_free_f,
 *        #VRNA_MATH_FUN_GAUSSIAN_OPTION_DEFAULT
 *
 *  @param  a             The arbitrary constant @f$ a @f$
 *  @param  b             The arbitrary constant @f$ b @f$
 *  @param  c             The arbitrary non-zero constant @f$ c @f$
 *  @param  options       Additional options that change the behavior of the callback function
 *  @param  options_p     A pointer to store the address of the options data structure
 *  @param  options_free  A pointer to store the address of the @c free function that releases the memory of the options data structure
 *  @return               A callback function that performs transformation through a Gaussian function
 */
vrna_math_fun_f
vrna_math_fun_gaussian_opt(double                   a,
                           double                   b,
                           double                   c,
                           unsigned int             options,
                           vrna_math_fun_opt_t      *options_p,
                           vrna_math_fun_opt_free_f *options_free);


/**
 *  @}
 */

#endif
