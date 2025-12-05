/* Data transform - Binning and mapping */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "ViennaRNA/datastructures/array.h"
#include "ViennaRNA/utils/basic.h"

#include "ViennaRNA/math/functions.h"

#ifdef __GNUC__
# define INLINE inline
#else
# define INLINE
#endif

typedef struct {
  double (*thresholds)[2];
  size_t        num_thresholds;
  unsigned int  options;
  double        v_min;
  double        v_max;
} fun_bin_opt_t;


/*
 #################################
 # PRIVATE FUNCTION DECLARATIONS #
 #################################
 */
PRIVATE void
bin_option_free(vrna_math_fun_dbl_opt_t options);


PRIVATE double
fun_bin_opt(double                  v,
            vrna_math_fun_dbl_opt_t options);


INLINE PRIVATE double
fun_bin(double        v,
        double        (*thresholds)[2],
        size_t        num_thresholds,
        double        oolb_value,
        double        ooub_value,
        unsigned int  options);


/*
 #################################
 # BEGIN OF FUNCTION DEFINITIONS #
 #################################
 */
PUBLIC vrna_math_fun_dbl_f
vrna_math_fun_dbl_bin_opt(double                    (*thresholds)[2],
                          unsigned int                  thresholds_num,
                          double                        oolb_value,
                          double                        ooub_value,
                          unsigned int                  options,
                          vrna_math_fun_dbl_opt_t       *fun_options_p,
                          vrna_math_fun_dbl_opt_free_f  *fun_options_free)
{
  vrna_math_fun_dbl_f cb = NULL;

  if ((thresholds) &&
      (thresholds_num > 1) &&
      (fun_options_p) &&
      (fun_options_free)) {
    fun_bin_opt_t *o = (fun_bin_opt_t *)vrna_alloc(sizeof(fun_bin_opt_t));

    o->options  = options;
    o->v_min    = oolb_value;
    o->v_max    = ooub_value;

    o->num_thresholds = thresholds_num;
    o->thresholds     = (double (*)[2])vrna_alloc(thresholds_num * sizeof(*(o->thresholds)));
    o->thresholds     = (double (*)[2])memcpy(o->thresholds, thresholds,
                                              thresholds_num * sizeof(*(o->thresholds)));

    cb                = fun_bin_opt;
    *fun_options_p    = (vrna_math_fun_dbl_opt_t)o;
    *fun_options_free = bin_option_free;
  }

  return cb;
}


PUBLIC double
vrna_math_fun_dbl_bin(double        v,
                      double        (*thresholds)[2],
                      size_t        num_thresholds,
                      double        oolb_value,
                      double        ooub_value,
                      unsigned int  options)
{
  if ((thresholds) &&
      (num_thresholds > 1)) {
    return fun_bin(v,
                   thresholds,
                   num_thresholds,
                   oolb_value,
                   ooub_value,
                   options);
  } else {
    vrna_log_warning("Missing thresholds in bin! Returning input value instead");
  }

  return v;
}


/*
 #####################################
 # BEGIN OF STATIC HELPER FUNCTIONS  #
 #####################################
 */
PRIVATE void
bin_option_free(vrna_math_fun_dbl_opt_t options)
{
  fun_bin_opt_t *o = (fun_bin_opt_t *)options;

  free(o->thresholds);
  free(o);
}


PRIVATE double
fun_bin_opt(double                  v,
            vrna_math_fun_dbl_opt_t options)
{
  fun_bin_opt_t *o = (fun_bin_opt_t *)options;

  return fun_bin(v,
                 o->thresholds,
                 o->num_thresholds,
                 o->v_min,
                 o->v_max,
                 o->options);
}


INLINE PRIVATE double
fun_bin(double        v,
        double        (*thresholds)[2],
        size_t        num_thresholds,
        double        oolb_value,
        double        ooub_value,
        unsigned int  options)
{
  if (v < thresholds[0][0])
    return (options &
            VRNA_MATH_FUN_BIN_OPTION_MAP_OUTOF_LOWERBOUND) ? thresholds[0][1] : oolb_value;

  size_t start, end, mid, last;

  last  = num_thresholds - 1;
  start = 1;
  end   = last;

  /* binary search for bin the value will be mapped to */
  while (start <= end) {
    mid = start + (end - start) / 2;

    if ((v >= thresholds[mid - 1][0]) &&
        (v < thresholds[mid][0])) {
      if (options & VRNA_MATH_FUN_BIN_OPTION_PROJECT) {
        double  diff_source = thresholds[mid][0] - thresholds[mid - 1][0];
        double  diff_target = thresholds[mid][1] - thresholds[mid - 1][1];
        return (v - thresholds[mid - 1][0]) /
               diff_source *
               diff_target +
               thresholds[mid - 1][1];
      } else {
        return thresholds[mid - 1][1];
      }
    }

    if (v >= thresholds[mid][0])
      start = mid + 1;
    else
      end = mid - 1;
  }

  if (v > thresholds[last][0])
    return (options &
            VRNA_MATH_FUN_BIN_OPTION_MAP_OUTOF_UPPERBOUND) ? thresholds[last][1] : ooub_value;
  else
    return thresholds[last][1];
}
