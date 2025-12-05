#ifndef VIENNA_RNA_PACKAGE_MATH_FUNCTIONS_H
#define VIENNA_RNA_PACKAGE_MATH_FUNCTIONS_H


/**
 *  @file     ViennaRNA/math/functions.h
 *  @ingroup  math_scalar
 *  @brief    This module provides mathematical functions
 */

/**
 *  @addtogroup math_scalar Mathematical Functions operating on Scalar Values
 *  @{
 */


/**
 *  @brief  Option data structure type for mathematical functions
 */
typedef   void *vrna_math_fun_dbl_opt_t;

/**
 *  @brief  Callback function to release any memory occupied by the mathematical function option #vrna_math_fun_dbl_opt_t
 */
typedef   vrna_auxdata_free_f vrna_math_fun_dbl_opt_free_f;


/**
 *  @brief Callback function prototype for scalar mathematical functions of the form @f$ y = f(x) @f$
 *
 *  This is the prototype for callback functions used for scalar mathematical functions. Its
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
 *  This callback handles mathematical functions working on scalar *double* values.
 *  @endparblock
 *
 *  @see  vrna_data_lin_transform(),
 *        vrna_math_fun_dbl_bin_opt(), vrna_math_fun_dbl_linear_opt(), vrna_math_fun_dbl_log_opt(),
 *        vrna_math_fun_dbl_gaussian_opt(), vrna_math_fun_dbl_logistic_opt(), vrna_math_fun_dbl_kde_opt()
 *
 * @param value     The value @f$ x @f$ that is to be evaluated by the function @f$ f(x) @f$
 * @param options   A pointer to a data structure specific for the mathematical function
 * @return          The evaluated value @f$ y @f$
 */
typedef double (*vrna_math_fun_dbl_f) (double                   value,
                                       vrna_math_fun_dbl_opt_t  options);


/**
 *  @brief  Options flag for linear data binning to activate data projection (mapping) instead of actual binning
 *  @see vrna_math_fun_dbl_bin(), vrna_math_fun_dbl_bin_opt()
 */
#define VRNA_MATH_FUN_BIN_OPTION_PROJECT              (1 << 0)


/**
 *  @brief  Options flag for linear data binning to indicate that values out-of-upper-bound are to be mapped to the respective domain limit
 *  @see vrna_math_fun_dbl_bin(), vrna_math_fun_dbl_bin_opt()
 */
#define VRNA_MATH_FUN_BIN_OPTION_MAP_OUTOF_UPPERBOUND (1 << 1)


/**
 *  @brief  Options flag for linear data binning to indicate that values out-of-lower-bound are to be mapped to the respective domain limit
 *  @see vrna_math_fun_dbl_bin(), vrna_math_fun_dbl_bin_opt()
 */
#define VRNA_MATH_FUN_BIN_OPTION_MAP_OUTOF_LOWERBOUND (1 << 2)


/**
 *  @brief  Options flag for linear data binning to indicate default settings
 *  @see vrna_math_fun_dbl_bin(), vrna_math_fun_dbl_bin_opt()
 */
#define VRNA_MATH_FUN_BIN_OPTION_DEFAULT              0U


/**
 *  @brief  Bin a scalar input value or project it from one domain into another
 *
 *  This function can be used to map continuous data from one domain (*source*) into
 *  discrete or continuous data of another domain (*target*). User-defined domain boundaries
 *  are provided by the @p thresholds argument, which boils down to a list of pairs of values,
 *  where the first value is the boundary in the source domain. The second value is the
 *  corresponding boundary of the target domain. The first and last pair of domain boundaries
 *  denote the lower and the upper domain limits, respectively.
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
 *  @see  vrna_data_lin_transform_opt(),
 *        #VRNA_MATH_FUN_BIN_OPTION_DEFAULT, #VRNA_MATH_FUN_BIN_OPTION_PROJECT,
 *        #VRNA_MATH_FUN_BIN_OPTION_MAP_OUTOF_UPPERBOUND, #VRNA_MATH_FUN_BIN_OPTION_MAP_OUTOF_LOWERBOUND,
 *        vrna_math_fun_dbl_linear(), vrna_math_fun_dbl_log(), vrna_math_fun_dbl_kde(),
 *        vrna_math_fun_dbl_gaussian()
 *
 *  @param  value             The input value
 *  @param  thresholds        A pointer to an array of data pairs holding the source and target domain boundaries
 *  @param  thresholds_num    The number of domain boundary pairs available in @p thresholds
 *  @param  oolb_value        Out-of-lower-bound value
 *  @param  ooub_value        Out-of-upper-bound value
 *  @param  options           Additional options that change the behavior of this function
 *  @return                   The projection of the input value into the target domain
 */
double
vrna_math_fun_dbl_bin(double        value,
                      double              (*thresholds)[2],
                      size_t        thresholds_num,
                      double        oolb_value,
                      double        ooub_value,
                      unsigned int  options);


/**
 *  @brief  Retrieve a function pointer that performs (discrete) binning and more
 *
 *  This function yields a linear data transform callback and the associated transform options
 *  data structure suitable for usage in vrna_data_lin_transform(). The transformation this
 *  callback performs can be described as (discrete) binning or data bucketing and essentially
 *  encapsulates the vrna_math_fun_dbl_bin() function.
 *
 *  The @p fun_options_p and @p fun_options_free pointers are used as additional
 *  output to obtain the addresses of the transformation option data structure that has to
 *  be provided to the vrna_data_lin_transform() function and a function pointer to release
 *  the memory of the option data structure once it is not required anymore.
 *
 *  @see  vrna_math_fun_dbl_bin(), vrna_data_lin_transform(),
 *        #vrna_math_fun_dbl_f, #vrna_math_fun_dbl_opt_t, #vrna_math_fun_dbl_opt_free_f,
 *        #VRNA_MATH_FUN_BIN_OPTION_DEFAULT, #VRNA_MATH_FUN_BIN_OPTION_PROJECT,
 *        #VRNA_MATH_FUN_BIN_OPTION_MAP_OUTOF_UPPERBOUND, #VRNA_MATH_FUN_BIN_OPTION_MAP_OUTOF_LOWERBOUND,
 *        vrna_math_fun_dbl_linear_opt(), vrna_math_fun_dbl_log_opt(), vrna_math_fun_dbl_kde_opt(),
 *        vrna_math_fun_dbl_gaussian_opt()
 *
 *  @param  thresholds          A pointer to an array of data pairs holding the source and target domain boundaries
 *  @param  thresholds_num      The number of domain boundary pairs available in @p thresholds
 *  @param  oolb_value          Out-of-lower-bound value
 *  @param  ooub_value          Out-of-upper-bound value
 *  @param  options             Additional options that change the behavior of the callback function
 *  @param  fun_options_p       A pointer to store the address of the options data structure
 *  @param  fun_options_free    A pointer to store the address of the @c free function that releases the memory of the options data structure
 *  @return                     A callback function that performs (discrete) data binning
 */
vrna_math_fun_dbl_f
vrna_math_fun_dbl_bin_opt(double                        (*thresholds)[2],
                          unsigned int                  thresholds_num,
                          double                        oolb_value,
                          double                        ooub_value,
                          unsigned int                  options,
                          vrna_math_fun_dbl_opt_t       *fun_options_p,
                          vrna_math_fun_dbl_opt_free_f  *fun_options_free);


/**
 *  @brief  Options flag for transforming linear data using a linear function that enables log-transform of the source value
 *  @see vrna_math_fun_dbl_linear(), vrna_math_fun_dbl_linear_opt()
 */
#define VRNA_MATH_FUN_LINEAR_OPTION_LOG         (1 << 0)

/**
 *  @brief  Options flag for transforming linear data using linear function that indicates default settings
 *  @see vrna_math_fun_dbl_linear(), vrna_math_fun_dbl_linear_opt()
 */
#define VRNA_MATH_FUN_LINEAR_OPTION_DEFAULT     0U


/**
 *  @brief  Transform a scalar input value through a linear function
 *
 *  This function can be used to transform an input value (parameter @p x) into
 *  an output value by applying a linear function of the form
 *  @f[ y = a + b \cdot f(x) @f]
 *  where @f$ x @f$ is the input value (@c source), @f$ a @f$ and @f$ b @f$ are *intercept*
 *  and *slope*, and @f$ y @f$ is the output value (@c target). By default, the function
 *  @f$ f(x) = x @f$ is the identity function. However, the callback can be instructed
 *  to apply a log-transform, @f$ f(x) = log x @f$, instead. Control over the behavior
 *  of @f$ f(x) @f$ can be gained by providing a corresponding option flag to the
 *  @p options argument, e.g. #VRNA_MATH_FUN_LINEAR_OPTION_LOG.
 *
 *  @see  vrna_math_fun_dbl_linear_opt(),
 *        #VRNA_MATH_FUN_LINEAR_OPTION_DEFAULT, #VRNA_MATH_FUN_LINEAR_OPTION_LOG,
 *        vrna_math_fun_dbl_linear(), vrna_math_fun_dbl_log(), vrna_math_fun_dbl_kde(),
 *        vrna_math_fun_dbl_gaussian()
 *
 *  @param  x           The input value
 *  @param  slope       The slope of the linear function
 *  @param  intercept   The intercept of the linear function
 *  @param  options     Additional options that change the behavior of the callback function
 *  @return             The transformed output value @f$ y @f$
 */
double
vrna_math_fun_dbl_linear(double       x,
                         double       slope,
                         double       intercept,
                         unsigned int options);


/**
 *  @brief  Retrieve a function pointer that applies a linear function
 *
 *  This function yields a linear data transform callback and the associated transform options
 *  data structure suitable for usage in vrna_data_lin_transform(). The callback applies
 *  a linear function @f$ y = a + b \cdot f(x) @f$ and essentially encapsulates the
 *  vrna_math_fun_dbl_linear() function.
 *
 *  The @p fun_options_p and @p fun_options_free pointers are used as additional
 *  output to obtain the addresses of the transformation option data structure that has to
 *  be provided to the vrna_data_lin_transform() function and a function pointer to release
 *  the memory of the option data structure once it is not required anymore.
 *
 *  @see  vrna_math_fun_dbl_linear(), vrna_data_lin_transform(),
 *        #vrna_math_fun_dbl_f, #vrna_math_fun_dbl_opt_t, #vrna_math_fun_dbl_opt_free_f,
 *        #VRNA_MATH_FUN_LINEAR_OPTION_DEFAULT, #VRNA_MATH_FUN_LINEAR_OPTION_LOG,
 *        vrna_math_fun_dbl_bin_opt(), vrna_math_fun_dbl_log_opt(), vrna_math_fun_dbl_kde_opt(),
 *        vrna_math_fun_dbl_gaussian_opt()
 *
 *  @param  slope               The slope of the linear function
 *  @param  intercept           The intercept of the linear function
 *  @param  options             Additional options that change the behavior of the callback function
 *  @param  fun_options_p       A pointer to store the address of the options data structure
 *  @param  fun_options_free    A pointer to store the address of the @c free function that releases the memory of the options data structure
 *  @return                     A callback function that performs transformation through a linear model
 */
vrna_math_fun_dbl_f
vrna_math_fun_dbl_linear_opt(double                       slope,
                             double                       intercept,
                             unsigned int                 options,
                             vrna_math_fun_dbl_opt_t      *fun_options_p,
                             vrna_math_fun_dbl_opt_free_f *fun_options_free);


/**
 *  @brief  Options flag for transforming linear data using log transform to use a non-default base
 *  @see vrna_math_fun_dbl_log(), vrna_math_fun_dbl_log_opt()
 */
#define VRNA_MATH_FUN_LOG_OPTION_NONDEFAULT_BASE    (1 << 0)


/**
 *  @brief  Options flag for transforming linear data using log transform that indicates default settings
 *  @see vrna_math_fun_dbl_log(), vrna_math_fun_dbl_log_opt()
 */
#define VRNA_MATH_FUN_LOG_OPTION_DEFAULT            0


/**
 *  @brief  Transform a scalar input value through a logarithmic function
 *
 *  This function transforms a scalar input value by computing the logarithm function
 *  @f[ y = \log_b x @f]
 *  where @f$ x @f$ is the input value (@c source), @f$ b @f$ is the base, and @f$ y @f$
 *  is the output value (@c target). By default, the *natural* logarithm is used, i.e.
 *  @f$ b = \mathrm{e} @f$. The @p base argument in combination with the
 *  #VRNA_MATH_FUN_LOG_OPTION_NONDEFAULT_BASE flag supplied to the @p options argument
 *  can be used to change the base to any other number.
 *
 *  The function returns the user-specified out-of-bounds value @p oob_value
 *  if @f$ x \le 0 @f$.
 *
 *  @see  vrna_math_fun_dbl_log_opt(),
 *        #VRNA_MATH_FUN_LOG_OPTION_DEFAULT, #VRNA_MATH_FUN_LOG_OPTION_NONDEFAULT_BASE,
 *        vrna_math_fun_dbl_linear(), vrna_math_fun_dbl_log(), vrna_math_fun_dbl_kde(),
 *        vrna_math_fun_dbl_gaussian()
 *
 *  @param  x           The input value
 *  @param  base        The base @f$ b @f$ of the logarithm (only used if #VRNA_MATH_FUN_LOG_OPTION_NONDEFAULT_BASE is passed to @p options
 *  @param  oob_value   Out-of-bound value
 *  @param  options     Additional options that change the behavior of the callback function
 *  @return             The transformed *target* value @f$ y @f$
 */
double
vrna_math_fun_dbl_log(double        x,
                      double        base,
                      double        oob_value,
                      unsigned int  options);


/**
 *  @brief  Retrieve a function pointer a logarithmic function
 *
 *  This function yields a linear data transform callback and the associated transform options
 *  data structure suitable for usage in vrna_data_lin_transform(). The callback applies
 *  a logarithmic function @f$ y = \log_b (x + c) @f$ and essentially encapsulates the
 *  vrna_math_fun_dbl_log() function. Note, that here, the @p value_shift argument corresponds
 *  to @f$ c @f$ in the above formula and in most cases should be @c 0. The callback then
 *  returns the out-of-bounds value @p oob_value if @f$ x + c \le 0 @f$.
 *
 *  The @p fun_options_p and @p fun_options_free pointers are used as additional
 *  output to obtain the addresses of the transformation option data structure that has to
 *  be provided to the vrna_data_lin_transform() function and a function pointer to release
 *  the memory of the option data structure once it is not required anymore.
 *
 *  @see  vrna_math_fun_dbl_log(), vrna_data_lin_transform(),
 *        #vrna_math_fun_dbl_f, #vrna_math_fun_dbl_opt_t, #vrna_math_fun_dbl_opt_free_f,
 *        #VRNA_MATH_FUN_LOG_OPTION_DEFAULT, #VRNA_MATH_FUN_LOG_OPTION_NONDEFAULT_BASE,
 *        vrna_math_fun_dbl_bin_opt(), vrna_math_fun_dbl_linear_opt(), vrna_math_fun_dbl_kde_opt(),
 *        vrna_math_fun_dbl_gaussian_opt()
 *
 *  @param  value_shift         The shift value @f$ c @f$
 *  @param  base                The base @f$ b @f$ of the logarithm (only used if #VRNA_MATH_FUN_LOG_OPTION_NONDEFAULT_BASE is passed to @p options
 *  @param  oob_value           Out-of-bound value
 *  @param  options             Additional options that change the behavior of the callback function
 *  @param  fun_options_p       A pointer to store the address of the options data structure
 *  @param  fun_options_free    A pointer to store the address of the @c free function that releases the memory of the options data structure
 *  @return                     A callback function that performs transformation through a linear model
 */
vrna_math_fun_dbl_f
vrna_math_fun_dbl_log_opt(double                        value_shift,
                          double                        base,
                          double                        oob_value,
                          unsigned int                  options,
                          vrna_math_fun_dbl_opt_t       *fun_options_p,
                          vrna_math_fun_dbl_opt_free_f  *fun_options_free);


/**
 *  @brief  Options flag for transforming linear data using logistic function that indicates default settings
 *  @see    vrna_math_fun_dbl_logistic(), vrna_math_fun_dbl_logistic_opt()
 */
#define VRNA_MATH_FUN_LOGISTIC_OPTION_DEFAULT     0U


/**
 *  @brief  Transform a scalar input value through a logistic function
 *
 *  This function transforms a scalar input value by applying a logistic function of the form
 *  @f[ y = \frac{L}{1 + e^{-k \cdot (x - x_0)}}  @f]
 *  where @f$ x @f$ is the input value (@c source), @f$ x_0 @f$ is the mid point (@p mid_point)
 *  of the function, @f$ k @f$ is the logistic growth rate (@p growth_rate), and @f$ L @f$
 *  is the supremum (or carrying capacity) of the function (@p supremum). The standard logistic
 *  function is defined as @f$ L = 1 @f$, @f$ k = 1 @f$. and @f$ x_0 = 0 @f$, i.e.
 *  @f[ y = \frac{1}{1 + e^{-x}} @f]
 *
 *  @see  vrna_math_fun_dbl_logistic_opt(),
 *        #VRNA_MATH_FUN_LOGISTIC_OPTION_DEFAULT,
 *        vrna_math_fun_dbl_bin(), vrna_math_fun_dbl_linear(), vrna_math_fun_dbl_log(),
 *        vrna_math_fun_dbl_kde(), vrna_math_fun_dbl_gaussian()
 *
 *  @param  x             The input value
 *  @param  mid_point     The midpoint of the function
 *  @param  supremum      The supremum of the function
 *  @param  growth_rate   The growth rate of the function
 *  @param  options       Additional options that change the behavior of the callback function
 *  @return               The transformed *target* value @f$ y @f$
 */
double
vrna_math_fun_dbl_logistic(double       x,
                           double       mid_point,
                           double       supremum,
                           double       growth_rate,
                           unsigned int options);


/**
 *  @brief  Retrieve a function pointer that applies a logistic function
 *
 *  This function yields a linear data transform callback and the associated transform options
 *  data structure suitable for usage in vrna_data_lin_transform(). The callback applies
 *  a logistic function of the form
 *  @f[ y = \frac{L}{1 + e^{-k \cdot (x - x_0)}}  @f]
 *  and essentially encapsulates the vrna_math_fun_dbl_logistic() function.
 *
 *  The @p fun_options_p and @p fun_options_free pointers are used as additional
 *  output to obtain the addresses of the transformation option data structure that has to
 *  be provided to the vrna_data_lin_transform() function and a function pointer to release
 *  the memory of the option data structure once it is not required anymore.
 *
 *  @see  vrna_math_fun_dbl_logistic(), vrna_data_lin_transform(),
 *        #vrna_math_fun_dbl_f, #vrna_math_fun_dbl_opt_t, #vrna_math_fun_dbl_opt_free_f,
 *        #VRNA_MATH_FUN_LOGISTIC_OPTION_DEFAULT,
 *        vrna_math_fun_dbl_bin_opt(), vrna_math_fun_dbl_linear_opt(), vrna_math_fun_dbl_kde_opt(),
 *        vrna_math_fun_dbl_gaussian_opt(), vrna_math_fun_dbl_log_opt()
 *
 *  @param  mid_point           The midpoint of the function
 *  @param  supremum            The supremum of the function
 *  @param  growth_rate         The growth rate of the function
 *  @param  options             Additional options that change the behavior of the callback function
 *  @param  fun_options_p       A pointer to store the address of the options data structure
 *  @param  fun_options_free    A pointer to store the address of the @c free function that releases the memory of the options data structure
 *  @return                     A callback function that performs transformation through a logistic function
 */
vrna_math_fun_dbl_f
vrna_math_fun_dbl_logistic_opt(double                       mid_point,
                               double                       supremum,
                               double                       growth_rate,
                               unsigned int                 options,
                               vrna_math_fun_dbl_opt_t      *fun_options_p,
                               vrna_math_fun_dbl_opt_free_f *fun_options_free);


/**
 *  @brief  Options flag for transforming linear data using a kernel density estimate (KDE) that indicates default settings
 *  @see vrna_math_fun_dbl_kde(), vrna_math_fun_dbl_kde_opt()
 */
#define VRNA_MATH_FUN_KDE_OPTION_DEFAULT      0U


/**
 *  @brief  Estimate the probability density for a scalar input value through kernel density estimation (KDE)
 *
 *  This function estimates the probability density function of a random input variable
 *  given a set of samples and a smoothing kernel. The estiamte of the probability
 *  density for a random variable @f$ x @f$ is then of the form
 *  @f[ y = \hat{f}_h(x) = \frac{1}{n \cdot h} \sum_{i = 1}^{n} K \big( \frac{x - x_i}{h} \big) @f]
 *  with a set of @f$ n > 0 @f$ samples @f$ \{ x_1, x_2, \ldots, x_n \} @f$ (@p samples),
 *  @f$ K @f$ (@p kernel) as a non-negative kernel function, and @f$ h > 0 @f$ is the bandwidth
 *  (@p bandwidth).
 *
 *  @note   If no kernel is provided, i.e. @c NULL is passed to @p kernel, the KDE uses a
 *          standard Gaussian kernel by providing vrna_math_fun_dbl_gaussian()
 *          with @f$ a = \frac{1}{\sqrt{2 \pi}}, b = 0, c = 1 @f$. In this case, the
 *          bandwidth (@p bandwidth) will be automatically replaced by a rule-of-thumb bandwidth
 *          estimator @f$ h = ( \frac{4}{3n} )^{\frac{1}{5}} @f$
 *
 *  @see  vrna_math_fun_dbl_kde_opt(),
 *        #VRNA_MATH_FUN_KDE_OPTION_DEFAULT,
 *        vrna_math_fun_dbl_bin(), vrna_math_fun_dbl_linear(), vrna_math_fun_dbl_log(),
 *        vrna_math_fun_dbl_logistic(), vrna_math_fun_dbl_gaussian()
 *
 *  @param  x             The input value
 *  @param  samples       The samples of the distribution to estimate @f$ \{ x_1, x_2, \ldots, x_n \} @f$
 *  @param  num_samples   The number of samples (@f$ n @f$) provided by pointer @p samples
 *  @param  kernel        The kernel function @f$ K @f$ (maybe @c NULL)
 *  @param  kernel_data   The data that should be passed-through to the kernel function (maybe @c NULL)
 *  @param  bandwidth     The bandwidth @f$ h @f$
 *  @param  options       Additional options that change the behavior of the callback function
 *  @return               The estimated probability density @f$ \hat{f}_h(x) @f$
 */
double
vrna_math_fun_dbl_kde(double                  x,
                      double                  *samples,
                      size_t                  num_samples,
                      vrna_math_fun_dbl_f     kernel,
                      vrna_math_fun_dbl_opt_t kernel_data,
                      double                  bandwidth,
                      unsigned int            options);


/**
 *  @brief  Retrieve a callback function that applies a kernel density estimation (KDE)
 *
 *  This function yields a linear data transform callback and the associated transform options
 *  data structure suitable for usage in vrna_data_lin_transform(). The callback applies
 *  a kernel density estimation (KDE) and essentially encapsulates the vrna_math_fun_dbl_kde()
 *  function.
 *
 *  @note   Whenever the @p fun_options_free callback will be called, this function
 *          also releases the memory occupied by @p kernel_data through the @p kernel_data_free
 *          function if it is provided, i.e. not @c NULL.
 *
 *  The @p fun_options_p and @p fun_options_free pointers are used as additional
 *  output to obtain the addresses of the transformation option data structure that has to
 *  be provided to the vrna_data_lin_transform() function and a function pointer to release
 *  the memory of the option data structure once it is not required anymore.
 *
 *  @see  vrna_math_fun_dbl_logistic(), vrna_data_lin_transform(),
 *        #vrna_math_fun_dbl_f, #vrna_math_fun_dbl_opt_t, #vrna_math_fun_dbl_opt_free_f,
 *        #VRNA_MATH_FUN_KDE_OPTION_DEFAULT,
 *        vrna_math_fun_dbl_gaussian_opt(), vrna_math_fun_dbl_bin_opt(), vrna_math_fun_dbl_linear_opt(),
 *        vrna_math_fun_dbl_log_opt(), vrna_math_fun_dbl_logistic_opt
 *
 *  @param  samples           The samples of the distribution to estimate @f$ \{ x_1, x_2, \ldots, x_n \} @f$
 *  @param  num_samples       The number of samples (@f$ n @f$) provided by pointer @p samples
 *  @param  kernel            The kernel function @f$ K @f$ (maybe @c NULL)
 *  @param  kernel_data       The data that should be passed-through to the kernel function (maybe @c NULL)
 *  @param  kernel_data_free  A function to free the memory occupied by the data for the kernel function (maybe @c NULL)
 *  @param  bandwidth         The bandwidth @f$ h @f$
 *  @param  options           Additional options that change the behavior of the callback function
 *  @param  fun_options_p     A pointer to store the address of the options data structure
 *  @param  fun_options_free  A pointer to store the address of the @c free function that releases the memory of the options data structure
 *  @return                   A callback function that performs transformation through kernel density estimation
 */
vrna_math_fun_dbl_f
vrna_math_fun_dbl_kde_opt(double                        *samples,
                          size_t                        num_samples,
                          vrna_math_fun_dbl_f           kernel,
                          vrna_math_fun_dbl_opt_t       kernel_data,
                          vrna_math_fun_dbl_opt_free_f  kernel_data_free,
                          double                        bandwidth,
                          unsigned int                  options,
                          vrna_math_fun_dbl_opt_t       *fun_options_p,
                          vrna_math_fun_dbl_opt_free_f  *fun_options_free);


/**
 *  @brief  Options flag for transforming linear data using a Gaussian function that indicates default settings
 *  @see vrna_math_fun_dbl_gaussian(), vrna_math_fun_dbl_gaussian_opt()
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
 *  @see  vrna_math_fun_dbl_gaussian_opt()
 *
 *  @param  x   The input value @f$ x @f$
 *  @param  a   The arbitrary constant @f$ a @f$
 *  @param  b   The arbitrary constant @f$ b @f$
 *  @param  c   The arbitrary non-zero constant @f$ c @f$
 *  @return     The evaluation for the Gaussian function @f$ y = f(x) @f$
 */
double
vrna_math_fun_dbl_gaussian(double x,
                           double a,
                           double b,
                           double c);


/**
 *  @brief  Retrieve a function pointer that applies Gaussian function
 *
 *  This function yields a linear data transform callback and the associated transform options
 *  data structure suitable for usage in vrna_data_lin_transform(). The callback applies
 *  a Gaussian function of the form
 *  @f[ y = a \cdot \exp(- \frac{(x - b)^2}{2 c^2}) @f]
 *  with arbitrary constants @f$ a, b @f$ and non-zero @f$ c @f$ and @f$ x @f$ being the
 *  input value (@c source).
 *
 *  The @p fun_options_p and @p fun_options_free pointers are used as additional
 *  output to obtain the addresses of the transformation option data structure that has to
 *  be provided to the vrna_data_lin_transform() function and a function pointer to release
 *  the memory of the option data structure once it is not required anymore.
 *
 *  @see  vrna_math_fun_dbl_gaussian(),
 *        vrna_data_lin_transform(), #vrna_math_fun_dbl_f, #vrna_math_fun_dbl_opt_t, #vrna_math_fun_dbl_opt_free_f,
 *        #VRNA_MATH_FUN_GAUSSIAN_OPTION_DEFAULT
 *
 *  @param  a                   The arbitrary constant @f$ a @f$
 *  @param  b                   The arbitrary constant @f$ b @f$
 *  @param  c                   The arbitrary non-zero constant @f$ c @f$
 *  @param  options             Additional options that change the behavior of the callback function
 *  @param  fun_options_p       A pointer to store the address of the options data structure
 *  @param  fun_options_free    A pointer to store the address of the @c free function that releases the memory of the options data structure
 *  @return                     A callback function that performs transformation through a Gaussian function
 */
vrna_math_fun_dbl_f
vrna_math_fun_dbl_gaussian_opt(double                       a,
                               double                       b,
                               double                       c,
                               unsigned int                 options,
                               vrna_math_fun_dbl_opt_t      *fun_options_p,
                               vrna_math_fun_dbl_opt_free_f *fun_options_free);


/**
 *  @}
 */

#endif
