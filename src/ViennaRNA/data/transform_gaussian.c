/* Data transform - Gaussian function transform */

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
  double        a;
  double        b;
  double        c;
  double        domain[2];
  unsigned int  options;
  double        out_of_bounds_value;
} gaussian_transform_param_t;


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
transform_gaussian(double                    r,
                   vrna_data_lin_trans_opt_t options);


PRIVATE void
transform_gaussian_option_free(vrna_data_lin_trans_opt_t options);


/*
 #################################
 # BEGIN OF FUNCTION DEFINITIONS #
 #################################
 */
PUBLIC vrna_data_lin_trans_f
vrna_data_transform_method_gaussian(double                          a,
                                    double                          b,
                                    double                          c,
                                    double                          domain[4],
                                    double                          oob_value,
                                    unsigned int                    options,
                                    vrna_data_lin_trans_opt_t       *transform_options_p,
                                    vrna_data_lin_trans_opt_free_f  *transform_options_free)
{
  vrna_data_lin_trans_f  cb = NULL;

  if ((transform_options_p) &&
      (transform_options_free)) {
    gaussian_transform_param_t *o = (gaussian_transform_param_t *)vrna_alloc(sizeof(gaussian_transform_param_t));

    o->a                    = a;
    o->b                    = b;
    o->c                    = c;
    o->options              = options;
    o->out_of_bounds_value  = oob_value;

    if (domain) {
      (void)memcpy(&(o->domain[0]), &(domain[0]), sizeof(double) * 4);
    } else {
      o->domain[0] = o->domain[1] = o->domain[4] = o->domain[3] = 0.;
      o->options &= ~VRNA_TRANSFORM_GAUSSIAN_OPTION_ENFORCE_DOMAINS;
    }

    cb                      = transform_gaussian;
    *transform_options_p    = (vrna_data_lin_trans_opt_t)o;
    *transform_options_free = transform_gaussian_option_free;
  }

  return cb;
}


/*
 #####################################
 # BEGIN OF STATIC HELPER FUNCTIONS  #
 #####################################
 */
PRIVATE void
transform_gaussian_option_free(vrna_data_lin_trans_opt_t options)
{
  free(options);
}


PRIVATE double
transform_gaussian(double                    value,
                   vrna_data_lin_trans_opt_t options)
{
  double                      t;
  gaussian_transform_param_t  *o;

  o = (gaussian_transform_param_t *)options;

  if (o->options & VRNA_TRANSFORM_GAUSSIAN_OPTION_ENFORCE_DOMAIN_SOURCE) {
    if (value < o->domain[0]) {
      if (o->options & VRNA_TRANSFORM_GAUSSIAN_OPTION_MAP_SOURCE_LOW)
        value = o->domain[0];
      else
        return o->out_of_bounds_value;
    } else if (value > o->domain[1]) {
      if (o->options & VRNA_TRANSFORM_GAUSSIAN_OPTION_MAP_SOURCE_HIGH)
        value = o->domain[1];
      else
        return o->out_of_bounds_value;
    }
  }

  t = o->a * exp(-(value - o->b) * (value - o->b) / (2. * o->c * o->c));

  if (o->options & VRNA_TRANSFORM_GAUSSIAN_OPTION_ENFORCE_DOMAIN_TARGET) {
    if (t < o->domain[2]) {
      if (o->options & VRNA_TRANSFORM_GAUSSIAN_OPTION_MAP_TARGET_LOW)
        t = o->domain[2];
      else
        return o->out_of_bounds_value;
    } else if (t > o->domain[3]) {
      if (o->options & VRNA_TRANSFORM_GAUSSIAN_OPTION_MAP_TARGET_HIGH)
        t = o->domain[3];
      else
        return o->out_of_bounds_value;
    }
  }

  return t;
}
