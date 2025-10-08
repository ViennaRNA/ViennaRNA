/* Data transform - Gaussian function transform */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "ViennaRNA/utils/basic.h"
#include "ViennaRNA/data/transform.h"

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
} gaussian_transform_param_t;


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
gaussian(double x,
         double a,
         double b,
         double c);

PRIVATE double
gaussian_opt(double                     r,
             vrna_data_lin_trans_opt_t  options);


PRIVATE void
transform_gaussian_option_free(vrna_data_lin_trans_opt_t options);


/*
 #################################
 # BEGIN OF FUNCTION DEFINITIONS #
 #################################
 */
PUBLIC vrna_data_lin_trans_f
vrna_data_transform_method_gaussian(double                          a,
                                    double                          b,
                                    double                          c,
                                    unsigned int                    options,
                                    vrna_data_lin_trans_opt_t       *transform_options_p,
                                    vrna_data_lin_trans_opt_free_f  *transform_options_free)
{
  vrna_data_lin_trans_f  cb = NULL;

  if ((transform_options_p) &&
      (transform_options_free)) {
    gaussian_transform_param_t *o = (gaussian_transform_param_t *)vrna_alloc(sizeof(gaussian_transform_param_t));

    o->a                    = a;
    o->b                    = b;
    o->c                    = c;
    o->options              = options;

    cb                      = gaussian_opt;
    *transform_options_p    = (vrna_data_lin_trans_opt_t)o;
    *transform_options_free = transform_gaussian_option_free;
  }

  return cb;
}


PUBLIC double
vrna_math_fun_gaussian(double x,
                       double a,
                       double b,
                       double c)
{
  return gaussian(x, a, b, c);
}

/*
 #####################################
 # BEGIN OF STATIC HELPER FUNCTIONS  #
 #####################################
 */
PRIVATE void
transform_gaussian_option_free(vrna_data_lin_trans_opt_t options)
{
  free(options);
}


PRIVATE double
gaussian_opt(double                     value,
             vrna_data_lin_trans_opt_t  options)
{
  gaussian_transform_param_t  *o = (gaussian_transform_param_t *)options;

  return gaussian(value, o->a, o->b, o->c);
}


INLINE PRIVATE double
gaussian(double x,
         double a,
         double b,
         double c)
{
  return a * exp(-(x - b) * (x - b) / (2. * c * c));
}
