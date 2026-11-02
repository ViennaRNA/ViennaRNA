/* Data transform - Logit function transform (Log-odds function) */

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
  double        min_p;
  double        max_p;
  double        oobl_value;
  double        oobu_value;
  unsigned int  options;
} fun_logit_opt_t;

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
fun_logit_opt(double                   p,
              vrna_math_fun_dbl_opt_t  options);


INLINE PRIVATE double
fun_logit(double p,
          double min_p,
          double max_p,
          double oobl_value,
          double oobu_value,
          unsigned int options);


PRIVATE void
logit_option_free(vrna_math_fun_dbl_opt_t options);


/*
 #################################
 # BEGIN OF FUNCTION DEFINITIONS #
 #################################
 */
PUBLIC vrna_math_fun_dbl_f
vrna_math_fun_dbl_logit_opt(const double                  bounds[2],
                            const double                  oob_values[2],
                            unsigned int                  options,
                            vrna_math_fun_dbl_opt_t       *fun_options_p,
                            vrna_math_fun_dbl_opt_free_f  *fun_options_free)
{
  vrna_math_fun_dbl_f cb = NULL;

  if ((fun_options_p) &&
      (fun_options_free)) {
    fun_logit_opt_t *o = (fun_logit_opt_t *)vrna_alloc(sizeof(fun_logit_opt_t));

    if (bounds == NULL) {
      /* default bounds imply open bounds interval, i.e. without 0 or 1 as p */
      o->min_p  = 0;
      o->max_p  = 1;
      options   &= ~VRNA_MATH_FUN_LOGIT_OPTION_BOUNDS_INCLUSIVE;
    } else {
      o->min_p  = ((bounds[0] >= 0.) && (bounds[0] <= 1.)) ? bounds[0] : 0.;
      o->max_p  = ((bounds[1] >= 0.) && (bounds[1] <= 1.)) ? bounds[1] : 1.;

      /* sanity checks */
      if ((o->min_p == 0.) &&
          (options & VRNA_MATH_FUN_LOGIT_OPTION_BOUNDS_INCLUSIVE_LOWER))
        options &= ~VRNA_MATH_FUN_LOGIT_OPTION_BOUNDS_INCLUSIVE_LOWER;

      if ((o->max_p == 1.) &&
          (options & VRNA_MATH_FUN_LOGIT_OPTION_BOUNDS_INCLUSIVE_UPPER))
        options &= ~VRNA_MATH_FUN_LOGIT_OPTION_BOUNDS_INCLUSIVE_UPPER;
    }

    if (oob_values) {
      o->oobl_value   = oob_values[0];
      o->oobu_value   = oob_values[1];
    } else {
      oob_values = VRNA_MATH_FUN_LOGIT_OOB_VALUES_DEFAULT;
      o->oobl_value   = oob_values[0];
      o->oobu_value   = oob_values[1];
    }

    o->options      = options;

    cb                = fun_logit_opt;
    *fun_options_p    = (vrna_math_fun_dbl_opt_t)o;
    *fun_options_free = logit_option_free;
  }

  return cb;
}


PUBLIC double
vrna_math_fun_dbl_logit(double        p,
                        const double  bounds[2],
                        const double  oob_values[2],
                        unsigned int  options)
{
  double min_p, max_p, oob[2];

  if (bounds == NULL) {
    /* default bounds imply open bounds interval, i.e. without 0 or 1 as p */
    min_p  = 0;
    max_p  = 1;
    options   &= ~VRNA_MATH_FUN_LOGIT_OPTION_BOUNDS_INCLUSIVE;
  } else {
    min_p  = ((bounds[0] >= 0.) && (bounds[0] <= 1.)) ? bounds[0] : 0.;
    max_p  = ((bounds[1] >= 0.) && (bounds[1] <= 1.)) ? bounds[1] : 1.;

    /* sanity checks */
    if ((min_p == 0.) &&
        (options & VRNA_MATH_FUN_LOGIT_OPTION_BOUNDS_INCLUSIVE_LOWER))
      options &= ~VRNA_MATH_FUN_LOGIT_OPTION_BOUNDS_INCLUSIVE_LOWER;

    if ((max_p == 1.) &&
        (options & VRNA_MATH_FUN_LOGIT_OPTION_BOUNDS_INCLUSIVE_UPPER))
      options &= ~VRNA_MATH_FUN_LOGIT_OPTION_BOUNDS_INCLUSIVE_UPPER;
  }

  if (oob_values) {
    oob[0] = oob_values[0];
    oob[1] = oob_values[1];
  } else {
    oob[0] = VRNA_MATH_FUN_LOGIT_OOB_VALUES_DEFAULT[0];
    oob[1] = VRNA_MATH_FUN_LOGIT_OOB_VALUES_DEFAULT[1];
  }

  return fun_logit(p, min_p, max_p, oob[0], oob[1], options);
}


/*
 #####################################
 # BEGIN OF STATIC HELPER FUNCTIONS  #
 #####################################
 */
PRIVATE void
logit_option_free(vrna_math_fun_dbl_opt_t options)
{
  free(options);
}


PRIVATE double
fun_logit_opt(double                   p,
              vrna_math_fun_dbl_opt_t  options)
{
  fun_logit_opt_t *o = (fun_logit_opt_t *)options;

  return fun_logit(p,
                   o->min_p,
                   o->max_p,
                   o->oobl_value,
                   o->oobu_value,
                   o->options);
}


INLINE PRIVATE double
fun_logit(double        p,
          double        min_p,
          double        max_p,
          double        oobl_value,
          double        oobu_value,
          unsigned int  options)
{
  if (((options & VRNA_MATH_FUN_LOGIT_OPTION_BOUNDS_INCLUSIVE_LOWER) && (p < min_p)) ||
      (p <= min_p))
    return oobl_value;
  else if (((options & VRNA_MATH_FUN_LOGIT_OPTION_BOUNDS_INCLUSIVE_UPPER) && (p > max_p)) ||
           (p >= max_p))
    return oobu_value;
  return log(p) - log(1. - p);
}
