/* Data transform - Kernel density estimation */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "ViennaRNA/utils/basic.h"
#include "ViennaRNA/math/functions.h"

#ifdef __GNUC__
# define INLINE inline
#else
# define INLINE
#endif


typedef struct {
  vrna_array(double)        samples;
  vrna_math_fun_dbl_f           kernel;
  vrna_math_fun_dbl_opt_t       kernel_data;
  vrna_math_fun_dbl_opt_free_f  kernel_data_free;
  double                        bandwidth;
  unsigned int                  options;
} fun_kde_opt_t;


/*
 #################################
 # GLOBAL VARIABLES              #
 #################################
 */

/*
 #################################
 # PRIVATE VARIABLES             #
 #################################
 */

/*
 #################################
 # PRIVATE FUNCTION DECLARATIONS #
 #################################
 */
PRIVATE double
fun_kde_opt(double                  v,
            vrna_math_fun_dbl_opt_t options);


INLINE PRIVATE double
fun_kde(double                  v,
        double                  *samples,
        size_t                  num_samples,
        vrna_math_fun_dbl_f     kernel,
        vrna_math_fun_dbl_opt_t kernel_data,
        double                  bandwidth);


PRIVATE void
kde_option_free(vrna_math_fun_dbl_opt_t options);


PRIVATE double
gaussian_bandwidth_default(const double *data,
                           size_t       n);


/*
 #################################
 # BEGIN OF FUNCTION DEFINITIONS #
 #################################
 */
PUBLIC vrna_math_fun_dbl_f
vrna_math_fun_dbl_kde_opt(double                        *samples,
                          size_t                        num_samples,
                          vrna_math_fun_dbl_f           kernel,
                          vrna_math_fun_dbl_opt_t       kernel_data,
                          vrna_math_fun_dbl_opt_free_f  kernel_data_free,
                          double                        bandwidth,
                          unsigned int                  options,
                          vrna_math_fun_dbl_opt_t       *fun_options_p,
                          vrna_math_fun_dbl_opt_free_f  *fun_options_free)
{
  vrna_math_fun_dbl_f cb = NULL;

  if ((fun_options_p) &&
      (fun_options_free) &&
      (samples) &&
      (num_samples > 0)) {
    fun_kde_opt_t *o = (fun_kde_opt_t *)vrna_alloc(sizeof(fun_kde_opt_t));

    o->options = options;

    /* store the samples */
    vrna_array_init_size(o->samples, num_samples);

    for (size_t i = 0; i < num_samples; i++)
      vrna_array_append(o->samples, samples[i]);

    if (kernel) {
      o->kernel           = kernel;
      o->kernel_data      = kernel_data;
      o->kernel_data_free = kernel_data_free;
      o->bandwidth        = bandwidth;
    } else {
      o->kernel = vrna_math_fun_dbl_gaussian_opt(1. / sqrt(2 * PI),
                                                 0,
                                                 1.,
                                                 VRNA_MATH_FUN_GAUSSIAN_OPTION_DEFAULT,
                                                 &(o->kernel_data),
                                                 &(o->kernel_data_free));
      o->bandwidth = gaussian_bandwidth_default(o->samples,
                                                vrna_array_size(o->samples));
    }

    cb                = fun_kde_opt;
    *fun_options_p    = (vrna_math_fun_dbl_opt_t)o;
    *fun_options_free = kde_option_free;
  }

  return cb;
}


PUBLIC double
vrna_math_fun_dbl_kde(double                  value,
                      double                  *samples,
                      size_t                  num_samples,
                      vrna_math_fun_dbl_f     kernel,
                      vrna_math_fun_dbl_opt_t kernel_data,
                      double                  bandwidth,
                      unsigned int            options)
{
  double result = 0.;

  if ((samples) &&
      (num_samples > 0)) {
    vrna_math_fun_dbl_f           k;
    vrna_math_fun_dbl_opt_t       kd;
    vrna_math_fun_dbl_opt_free_f  kd_free;

    if (kernel) {
      k   = kernel;
      kd  = kernel_data;
    } else {
      k = vrna_math_fun_dbl_gaussian_opt(1. / sqrt(2 * PI),
                                         0,
                                         1.,
                                         VRNA_MATH_FUN_GAUSSIAN_OPTION_DEFAULT,
                                         &kd,
                                         &kd_free);
      bandwidth = gaussian_bandwidth_default(samples, num_samples);
    }

    result = fun_kde(value,
                     samples,
                     num_samples,
                     k,
                     kd,
                     bandwidth);
  } else {
    vrna_log_warning("Missing samples in KDE! Returning input value instead");
  }

  return value;
}


/*
 #####################################
 # BEGIN OF STATIC HELPER FUNCTIONS  #
 #####################################
 */
PRIVATE void
kde_option_free(vrna_math_fun_dbl_opt_t options)
{
  fun_kde_opt_t *o = (fun_kde_opt_t *)options;

  vrna_array_free(o->samples);

  if (o->kernel_data_free)
    o->kernel_data_free(o->kernel_data);

  free(options);
}


PRIVATE double
fun_kde_opt(double                  value,
            vrna_math_fun_dbl_opt_t options)
{
  fun_kde_opt_t *o = (fun_kde_opt_t *)options;

  return fun_kde(value,
                 o->samples,
                 vrna_array_size(o->samples),
                 o->kernel,
                 o->kernel_data,
                 o->bandwidth);
}


INLINE PRIVATE double
fun_kde(double                  value,
        double                  *samples,
        size_t                  num_samples,
        vrna_math_fun_dbl_f     kernel,
        vrna_math_fun_dbl_opt_t kernel_data,
        double                  bandwidth)
{
  double t = 0;

  for (size_t i = 0; i < num_samples; i++)
    t += kernel((value - samples[i]) / bandwidth,
                kernel_data);

  t /= (num_samples * bandwidth);

  return t;
}


/*  assuming Gaussian basis function, this is an optimal
 *  rule-of-thumb bandwidth minimizing the mean integrated
 *  squared error
 */
PRIVATE double
gaussian_bandwidth_default(const double *data,
                           size_t       n)
{
  double  factor, mu, std;
  size_t  i, count;

  mu    = 0.;
  std   = 0.;
  count = 0;

  factor = (pow((double)n, -1. / 5)) *
           (pow(4. / 3., 1. / 5.));

  /* compute mean */
  for (i = 0; i < n; i++)
    mu += data[i];

  mu /= n;

  /* compute standard deviation */
  for (i = 0; i < n; i++)
    std += (data[i] - mu) * (data[i] - mu);

  std = sqrt(std / ((double)n - 1));

  return factor * std;
}
