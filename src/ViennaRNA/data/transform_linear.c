/* Data transform - Linear model transform */

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
  double        slope;
  double        intercept;
  double        domain[4];
  double        out_of_bounds_value;
  unsigned int  options;
} linear_transform_param_t;

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
transform_linear(double                     value,
                 vrna_data_lin_trans_opt_t  options);


PRIVATE void
transform_lm_option_free(vrna_data_lin_trans_opt_t options);


/*
 #################################
 # BEGIN OF FUNCTION DEFINITIONS #
 #################################
 */
PUBLIC vrna_data_lin_trans_f
vrna_data_transform_method_lm(double                          slope,
                              double                          intercept,
                              double                          domain[4],
                              double                          oob_value,
                              unsigned int                    options,
                              vrna_data_lin_trans_opt_t       *transform_options_p,
                              vrna_data_lin_trans_opt_free_f  *transform_options_free)
{
  vrna_data_lin_trans_f  cb = NULL;

  if ((transform_options_p) &&
      (transform_options_free)) {
    linear_transform_param_t *o = (linear_transform_param_t *)vrna_alloc(sizeof(linear_transform_param_t));

    o->slope                = slope;
    o->intercept            = intercept;
    o->out_of_bounds_value  = oob_value;
    o->options              = options;

    if (domain) {
      (void)memcpy(&(o->domain[0]), &(domain[0]), sizeof(double) * 4);
    } else {
      o->domain[0] = o->domain[1] = o->domain[2] = o->domain[3] = 0.;
      o->options &= ~VRNA_TRANSFORM_LM_OPTION_ENFORCE_DOMAINS;
    }

    cb                      = transform_linear;
    *transform_options_p    = (vrna_data_lin_trans_opt_t)o;
    *transform_options_free = transform_lm_option_free;
  }

  return cb;
}


/*
 #####################################
 # BEGIN OF STATIC HELPER FUNCTIONS  #
 #####################################
 */
PRIVATE void
transform_lm_option_free(vrna_data_lin_trans_opt_t options)
{
  free(options);
}


PRIVATE double
transform_linear(double                     value,
                 vrna_data_lin_trans_opt_t  options)
{
  double                    t;
  linear_transform_param_t  *o = (linear_transform_param_t *)options;

  t = value;

  if (o->options & VRNA_TRANSFORM_LM_OPTION_ENFORCE_DOMAIN_SOURCE) {
    if (t < o->domain[0]) {
      if (o->options & VRNA_TRANSFORM_LM_OPTION_MAP_SOURCE_LOW)
        t = o->domain[0];
      else
        return o->out_of_bounds_value;
    } else if (t > o->domain[1]) {
      if (o->options & VRNA_TRANSFORM_LM_OPTION_MAP_SOURCE_HIGH)
        t = o->domain[1];
      else
        return o->out_of_bounds_value;
    }
  }

  if (o->options & VRNA_TRANSFORM_LM_OPTION_LOG)
    t = o->intercept + o->slope * log(t);
  else
    t = o->intercept + o->slope * t;

  if (o->options & VRNA_TRANSFORM_LM_OPTION_ENFORCE_DOMAIN_TARGET) {
    if (t < o->domain[2]) {
      return (o->options & VRNA_TRANSFORM_LM_OPTION_MAP_TARGET_LOW) ?
                o->domain[2] :
                o->out_of_bounds_value;
    } else if (t > o->domain[3]) {
      return (o->options & VRNA_TRANSFORM_LM_OPTION_MAP_TARGET_HIGH) ?
                o->domain[3] :
                o->out_of_bounds_value;
    }
  }

  return t;
}
