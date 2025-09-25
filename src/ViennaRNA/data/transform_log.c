/* Data transform - Log transform */

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
  double        shift;
  double        base;
  unsigned int  options;
  double        out_of_bounds_value;
} log_transform_param_t;

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
transform_log(double                    r,
              vrna_data_lin_trans_opt_t options);


PRIVATE void
transform_log_option_free(vrna_data_lin_trans_opt_t options);


/*
 #################################
 # BEGIN OF FUNCTION DEFINITIONS #
 #################################
 */
PUBLIC vrna_data_lin_trans_f
vrna_data_transform_method_log(double                         value_shift,
                               double                         base,
                               double                         oob_value,
                               unsigned int                   options,
                               vrna_data_lin_trans_opt_t      *transform_options_p,
                               vrna_data_lin_trans_opt_free_f *transform_options_free)
{
  vrna_data_lin_trans_f  cb = NULL;

  if ((transform_options_p) &&
      (transform_options_free)) {
    log_transform_param_t *o = (log_transform_param_t *)vrna_alloc(sizeof(log_transform_param_t));

    o->shift                = value_shift;
    o->base                 = base;
    o->options              = options;
    o->out_of_bounds_value  = oob_value;

    cb                      = transform_log;
    *transform_options_p    = (vrna_data_lin_trans_opt_t)o;
    *transform_options_free = transform_log_option_free;
  }

  return cb;
}


/*
 #####################################
 # BEGIN OF STATIC HELPER FUNCTIONS  #
 #####################################
 */
PRIVATE void
transform_log_option_free(vrna_data_lin_trans_opt_t options)
{
  free(options);
}


PRIVATE double
transform_log(double                    value,
              vrna_data_lin_trans_opt_t options)
{
  log_transform_param_t *o = (log_transform_param_t *)options;

  value += o->shift;

  if (value > 0) {
    value = log(value);

    if (o->options & VRNA_TRANSFORM_LOG_OPTION_NONDEFAULT_BASE)
      value /= log(o->base);

    return value;
  }

  return o->out_of_bounds_value;
}
