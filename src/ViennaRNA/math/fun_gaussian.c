/* Data transform - Gaussian function transform */

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
  double        a;
  double        b;
  double        c;
  unsigned int  options;
} fun_gaussian_opt_t;


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
INLINE PRIVATE double
fun_gaussian(double x,
             double a,
             double b,
             double c);


PRIVATE double
fun_gaussian_opt(double                   v,
                 vrna_math_fun_dbl_opt_t  options);


PRIVATE void
gaussian_option_free(vrna_math_fun_dbl_opt_t options);


/*
 #################################
 # BEGIN OF FUNCTION DEFINITIONS #
 #################################
 */
PUBLIC vrna_math_fun_dbl_f
vrna_math_fun_dbl_gaussian_opt(double                       a,
                               double                       b,
                               double                       c,
                               unsigned int                 options,
                               vrna_math_fun_dbl_opt_t      *fun_options_p,
                               vrna_math_fun_dbl_opt_free_f *fun_options_free)
{
  vrna_math_fun_dbl_f cb = NULL;

  if ((fun_options_p) &&
      (fun_options_free)) {
    fun_gaussian_opt_t *o = (fun_gaussian_opt_t *)vrna_alloc(sizeof(fun_gaussian_opt_t));

    o->a        = a;
    o->b        = b;
    o->c        = c;
    o->options  = options;

    cb                = fun_gaussian_opt;
    *fun_options_p    = (vrna_math_fun_dbl_opt_t)o;
    *fun_options_free = gaussian_option_free;
  }

  return cb;
}


PUBLIC double
vrna_math_fun_dbl_gaussian(double x,
                           double a,
                           double b,
                           double c)
{
  return fun_gaussian(x, a, b, c);
}


/*
 #####################################
 # BEGIN OF STATIC HELPER FUNCTIONS  #
 #####################################
 */
PRIVATE void
gaussian_option_free(vrna_math_fun_dbl_opt_t options)
{
  free(options);
}


PRIVATE double
fun_gaussian_opt(double                   value,
                 vrna_math_fun_dbl_opt_t  options)
{
  fun_gaussian_opt_t *o = (fun_gaussian_opt_t *)options;

  return fun_gaussian(value, o->a, o->b, o->c);
}


INLINE PRIVATE double
fun_gaussian(double x,
             double a,
             double b,
             double c)
{
  return a * exp(-(x - b) * (x - b) / (2. * c * c));
}
