/* Data transform - Kernel density estimation */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "ViennaRNA/utils/basic.h"
#include "ViennaRNA/data/transform.h"



typedef struct {
  vrna_array(double)              samples;
  vrna_data_lin_trans_f           kernel;
  vrna_data_lin_trans_opt_t       kernel_data;
  vrna_data_lin_trans_opt_free_f  kernel_data_free;
  double                          bandwidth;
  unsigned int                    options;
  double                          out_of_bounds_value;
} kde_transform_param_t;


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
transform_kde(double                    r,
              vrna_data_lin_trans_opt_t options);


PRIVATE void
transform_kde_option_free(vrna_data_lin_trans_opt_t options);


PRIVATE double
gaussian_bandwidth_default(const vrna_array(double) data);

/*
 #################################
 # BEGIN OF FUNCTION DEFINITIONS #
 #################################
 */
PUBLIC vrna_data_lin_trans_f
vrna_data_transform_method_kde(double                         *samples,
                               size_t                         num_samples,
                               vrna_data_lin_trans_f          kernel,
                               vrna_data_lin_trans_opt_t      kernel_data,
                               vrna_data_lin_trans_opt_free_f kernel_data_free,
                               double                         bandwidth,
                               unsigned int                   options,
                               vrna_data_lin_trans_opt_t      *transform_options_p,
                               vrna_data_lin_trans_opt_free_f *transform_options_free)
{
  vrna_data_lin_trans_f  cb = NULL;

  if ((transform_options_p) &&
      (transform_options_free) &&
      (samples) &&
      (num_samples > 0)) {
    kde_transform_param_t *o = (kde_transform_param_t *)vrna_alloc(sizeof(kde_transform_param_t));

    o->options              = options;

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
      o->kernel           = vrna_data_transform_method_gaussian(1. / sqrt(2 * PI),
                                                                0,
                                                                1.,
                                                                VRNA_TRANSFORM_GAUSSIAN_OPTION_DEFAULT,
                                                                &(o->kernel_data),
                                                                &(o->kernel_data_free));
      o->bandwidth        = gaussian_bandwidth_default(o->samples);
    }

    cb                      = transform_kde;
    *transform_options_p    = (vrna_data_lin_trans_opt_t)o;
    *transform_options_free = transform_kde_option_free;
  }

  return cb;
}


/*
 #####################################
 # BEGIN OF STATIC HELPER FUNCTIONS  #
 #####################################
 */
PRIVATE void
transform_kde_option_free(vrna_data_lin_trans_opt_t options)
{
  kde_transform_param_t *o = (kde_transform_param_t *)options;

  vrna_array_free(o->samples);

  if (o->kernel_data_free)
    o->kernel_data_free(o->kernel_data);

  free(options);
}


PRIVATE double
transform_kde(double                    value,
                   vrna_data_lin_trans_opt_t options)
{
  size_t                n;
  double                t, bandwidth;
  kde_transform_param_t *o  = (kde_transform_param_t *)options;

  n         = vrna_array_size(o->samples);
  t         = 0;
  bandwidth = o->bandwidth;

  for (size_t i = 0; i < n; i++)
    t += o->kernel((value - o->samples[i]) / bandwidth,
                    o->kernel_data);

  t /= (n * bandwidth);

  return t;
}


/*  assuming Gaussian basis function, this is an optimal
 *  rule-of-thumb bandwidth minimizing the mean integrated
 *  squared error
 */
PRIVATE double
gaussian_bandwidth_default(const vrna_array(double) data)
{
  double        factor, mu, std;
  unsigned int  i, count, n;

  mu    = 0.;
  std   = 0.;
  count = 0;
  n     = vrna_array_size(data);

  factor = (pow((double)n, -1. / 5)) *
            (pow(4. / 3., 1. / 5.));

  /* compute mean */
  for (i = 0; i < n; i++)
    mu += data[i];

  mu /= n;

  /* compute standard deviation */
  for (i = 0; i < n; i++)
    std += (data[i] - mu) * (data[i] - mu);

  std = sqrt(std / (n - 1));

  return (double)(factor * std);
}

