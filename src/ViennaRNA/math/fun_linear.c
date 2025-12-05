/* Data transform - Linear model transform */

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
  double        slope;
  double        intercept;
  unsigned int  options;
} fun_linear_opt_t;

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
fun_linear_opt(double                   value,
               vrna_math_fun_dbl_opt_t  options);


INLINE PRIVATE double
fun_linear(double       value,
           double       slope,
           double       intercept,
           unsigned int options);


PRIVATE void
linear_option_free(vrna_math_fun_dbl_opt_t options);


/*
 #################################
 # BEGIN OF FUNCTION DEFINITIONS #
 #################################
 */
PUBLIC vrna_math_fun_dbl_f
vrna_math_fun_dbl_linear_opt(double                       slope,
                             double                       intercept,
                             unsigned int                 options,
                             vrna_math_fun_dbl_opt_t      *fun_options_p,
                             vrna_math_fun_dbl_opt_free_f *fun_options_free)
{
  vrna_math_fun_dbl_f cb = NULL;

  if ((fun_options_p) &&
      (fun_options_free)) {
    fun_linear_opt_t *o = (fun_linear_opt_t *)vrna_alloc(sizeof(fun_linear_opt_t));

    o->slope      = slope;
    o->intercept  = intercept;
    o->options    = options;

    cb                = fun_linear_opt;
    *fun_options_p    = (vrna_math_fun_dbl_opt_t)o;
    *fun_options_free = linear_option_free;
  }

  return cb;
}


PUBLIC double
vrna_math_fun_dbl_linear(double       v,
                         double       slope,
                         double       intercept,
                         unsigned int options)
{
  return fun_linear(v, slope, intercept, options);
}


/*
 #####################################
 # BEGIN OF STATIC HELPER FUNCTIONS  #
 #####################################
 */
PRIVATE void
linear_option_free(vrna_math_fun_dbl_opt_t options)
{
  free(options);
}


PRIVATE double
fun_linear_opt(double                   value,
               vrna_math_fun_dbl_opt_t  options)
{
  fun_linear_opt_t *o = (fun_linear_opt_t *)options;

  return fun_linear(value,
                    o->slope,
                    o->intercept,
                    o->options);
}


INLINE PRIVATE double
fun_linear(double       value,
           double       slope,
           double       intercept,
           unsigned int options)
{
  return (options & VRNA_MATH_FUN_LINEAR_OPTION_LOG) ?
         intercept + slope * log(value) :
         intercept + slope * value;
}
