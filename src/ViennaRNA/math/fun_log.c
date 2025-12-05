/* Data transform - Log transform */

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
  double        shift;
  double        base;
  unsigned int  options;
  double        out_of_bounds_value;
} fun_log_opt_t;

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
fun_log_opt(double                  v,
            vrna_math_fun_dbl_opt_t options);


INLINE PRIVATE double
fun_log(double        v,
        double        base,
        double        oob_value,
        unsigned int  options);


PRIVATE void
log_option_free(vrna_math_fun_dbl_opt_t options);


/*
 #################################
 # BEGIN OF FUNCTION DEFINITIONS #
 #################################
 */
PUBLIC vrna_math_fun_dbl_f
vrna_math_fun_dbl_log_opt(double                        value_shift,
                          double                        base,
                          double                        oob_value,
                          unsigned int                  options,
                          vrna_math_fun_dbl_opt_t       *fun_options_p,
                          vrna_math_fun_dbl_opt_free_f  *fun_options_free)
{
  vrna_math_fun_dbl_f cb = NULL;

  if ((fun_options_p) &&
      (fun_options_free)) {
    fun_log_opt_t *o = (fun_log_opt_t *)vrna_alloc(sizeof(fun_log_opt_t));

    o->shift                = value_shift;
    o->base                 = base;
    o->options              = options;
    o->out_of_bounds_value  = oob_value;

    cb                = fun_log_opt;
    *fun_options_p    = (vrna_math_fun_dbl_opt_t)o;
    *fun_options_free = log_option_free;
  }

  return cb;
}


PUBLIC double
vrna_math_fun_dbl_log(double        value,
                      double        base,
                      double        oob_value,
                      unsigned int  options)
{
  return fun_log(value,
                 base,
                 oob_value,
                 options);
}


/*
 #####################################
 # BEGIN OF STATIC HELPER FUNCTIONS  #
 #####################################
 */
PRIVATE void
log_option_free(vrna_math_fun_dbl_opt_t options)
{
  free(options);
}


PRIVATE double
fun_log_opt(double                  value,
            vrna_math_fun_dbl_opt_t options)
{
  fun_log_opt_t *o = (fun_log_opt_t *)options;

  return fun_log(value + o->shift,
                 o->base,
                 o->out_of_bounds_value,
                 o->options);
}


INLINE PRIVATE double
fun_log(double        v,
        double        base,
        double        oob_value,
        unsigned int  options)
{
  if (v > 0) {
    v = log(v);

    if (options & VRNA_MATH_FUN_LOG_OPTION_NONDEFAULT_BASE)
      v /= log(base);

    return v;
  }

  return oob_value;
}
