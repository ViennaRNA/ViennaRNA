/* Data transform - Logistic function transform */

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
  double        mid_point;
  double        supremum;
  double        growth_rate;
  double        domain[2];
  unsigned int  options;
  double        out_of_bounds_value;
} logistic_transform_param_t;

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
/*
 * Log(1+x) for x > -1, ignore otherwise
 */
PRIVATE double
transform_logistic(double                    r,
                   vrna_data_lin_trans_opt_t options);


PRIVATE void
transform_logistic_option_free(vrna_data_lin_trans_opt_t options);


/*
 #################################
 # BEGIN OF FUNCTION DEFINITIONS #
 #################################
 */
PUBLIC vrna_data_lin_trans_f
vrna_data_transform_method_logistic(double                         mid_point,
                                    double                         supremum,
                                    double                         growth_rate,
                                    double                         domain[2],
                                    double                         oob_value,
                                    unsigned int                   options,
                                    vrna_data_lin_trans_opt_t      *transform_options_p,
                                    vrna_data_lin_trans_opt_free_f *transform_options_free)
{
  vrna_data_lin_trans_f  cb = NULL;

  if ((transform_options_p) &&
      (transform_options_free)) {
    logistic_transform_param_t *o = (logistic_transform_param_t *)vrna_alloc(sizeof(logistic_transform_param_t));

    o->mid_point                = mid_point;
    o->supremum                 = supremum;
    o->growth_rate                 = growth_rate;
    o->options              = options;
    o->out_of_bounds_value  = oob_value;

    if (domain) {
      (void)memcpy(&(o->domain[0]), &(domain[0]), sizeof(double) * 2);
    } else {
      o->domain[0] = o->domain[1] = 0.;
      o->options &= ~VRNA_TRANSFORM_LOGISTIC_OPTION_ENFORCE_DOMAIN;
    }

    cb                      = transform_logistic;
    *transform_options_p    = (vrna_data_lin_trans_opt_t)o;
    *transform_options_free = transform_logistic_option_free;
  }

  return cb;
}


/*
 #####################################
 # BEGIN OF STATIC HELPER FUNCTIONS  #
 #####################################
 */
PRIVATE void
transform_logistic_option_free(vrna_data_lin_trans_opt_t options)
{
  free(options);
}


PRIVATE double
transform_logistic(double                    value,
                   vrna_data_lin_trans_opt_t options)
{
  double t;
  logistic_transform_param_t *o = (logistic_transform_param_t *)options;

  t = value - o->mid_point;

  if (o->options & VRNA_TRANSFORM_LOGISTIC_OPTION_ENFORCE_DOMAIN) {
    if (t < o->domain[0]) {
      if (o->options & VRNA_TRANSFORM_LOGISTIC_OPTION_MAP_LOW)
        t = o->domain[0];
      else
        return o->out_of_bounds_value;
    } else if (t > o->domain[1]) {
      if (o->options & VRNA_TRANSFORM_LOGISTIC_OPTION_MAP_HIGH)
        t = o->domain[1];
      else
        return o->out_of_bounds_value;
    }
  }

  t = 1. + exp(-o->growth_rate * t);
  t = o->supremum / t;

  return t;
}
