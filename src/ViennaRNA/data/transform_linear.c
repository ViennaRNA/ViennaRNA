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
    o->options              = options;

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
  linear_transform_param_t  *o = (linear_transform_param_t *)options;

  return (o->options & VRNA_TRANSFORM_LM_OPTION_LOG) ?
            o->intercept + o->slope * log(value) :
            o->intercept + o->slope * value;
}
