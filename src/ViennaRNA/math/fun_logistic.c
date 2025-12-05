/* Data transform - Logistic function transform */

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
  double        mid_point;
  double        supremum;
  double        growth_rate;
  unsigned int  options;
} fun_logistic_opt_t;

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
fun_logistic_opt(double                   v,
                 vrna_math_fun_dbl_opt_t  options);


INLINE PRIVATE double
fun_logistic(double       v,
             double       mid_point,
             double       supremum,
             double       growth_rate,
             unsigned int options);


PRIVATE void
logistic_option_free(vrna_math_fun_dbl_opt_t options);


/*
 #################################
 # BEGIN OF FUNCTION DEFINITIONS #
 #################################
 */
PUBLIC vrna_math_fun_dbl_f
vrna_math_fun_dbl_logistic_opt(double                       mid_point,
                               double                       supremum,
                               double                       growth_rate,
                               unsigned int                 options,
                               vrna_math_fun_dbl_opt_t      *fun_options_p,
                               vrna_math_fun_dbl_opt_free_f *fun_options_free)
{
  vrna_math_fun_dbl_f cb = NULL;

  if ((fun_options_p) &&
      (fun_options_free)) {
    fun_logistic_opt_t *o = (fun_logistic_opt_t *)vrna_alloc(sizeof(fun_logistic_opt_t));

    o->mid_point    = mid_point;
    o->supremum     = supremum;
    o->growth_rate  = growth_rate;
    o->options      = options;

    cb                = fun_logistic_opt;
    *fun_options_p    = (vrna_math_fun_dbl_opt_t)o;
    *fun_options_free = logistic_option_free;
  }

  return cb;
}


PUBLIC double
vrna_math_fun_dbl_logistic(double       v,
                           double       mid_point,
                           double       supremum,
                           double       growth_rate,
                           unsigned int options)
{
  return fun_logistic(v,
                      mid_point,
                      supremum,
                      growth_rate,
                      options);
}


/*
 #####################################
 # BEGIN OF STATIC HELPER FUNCTIONS  #
 #####################################
 */
PRIVATE void
logistic_option_free(vrna_math_fun_dbl_opt_t options)
{
  free(options);
}


PRIVATE double
fun_logistic_opt(double                   value,
                 vrna_math_fun_dbl_opt_t  options)
{
  fun_logistic_opt_t *o = (fun_logistic_opt_t *)options;

  return fun_logistic(value,
                      o->mid_point,
                      o->supremum,
                      o->growth_rate,
                      o->options);
}


INLINE PRIVATE double
fun_logistic(double       v,
             double       mid_point,
             double       supremum,
             double       growth_rate,
             unsigned int options)
{
  v = 1. + exp(-growth_rate * (v - mid_point));

  return supremum / v;
}
