/* Data transform - Linear model transform */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "ViennaRNA/utils/basic.h"
#include "ViennaRNA/probing/basic.h"
#include "ViennaRNA/probing/transform.h"


typedef struct {
  double        slope;
  double        intercept;
  double        oob;
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
transform_linear(double value,
                 void   *options);


PRIVATE double
transform_linear_log(double value,
                     void   *options);


PRIVATE void
transform_lm_option_free(void *options);


/*
 #################################
 # BEGIN OF FUNCTION DEFINITIONS #
 #################################
 */
PUBLIC vrna_probing_transform_f
vrna_data_transform_method_lm(double               slope,
                              double               intercept,
                              double               oob_value,
                              unsigned char        options,
                              void                 **transform_options_p,
                              vrna_auxdata_free_f  *transform_options_free)
{
  vrna_probing_transform_f  cb = NULL;

  if ((transform_options_p) &&
      (transform_options_free)) {
    linear_transform_param_t *o = (linear_transform_param_t *)vrna_alloc(sizeof(linear_transform_param_t));

    o->slope      = slope;
    o->intercept  = intercept;
    o->oob        = oob_value;

    cb                      = (options & VRNA_TRANSFORM_LM_OPTION_LOG) ? transform_linear_log : transform_linear;
    *transform_options_p    = (void *)o;
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
transform_lm_option_free(void *options)
{
  free(options);
}


PRIVATE double
transform_linear(double value,
                 void   *options)
{
  linear_transform_param_t *o = (linear_transform_param_t *)options;

  return o->intercept + o->slope * value;
}


PRIVATE double
transform_linear_log(double value,
                     void   *options)
{
  linear_transform_param_t *o = (linear_transform_param_t *)options;

  return (value > 0.) ? o->intercept + o->slope * log(value) : o->oob;
}
