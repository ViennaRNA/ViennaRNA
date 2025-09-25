/* Data transform - Binning and mapping */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdio.h>
#include <stdlib.h>

#include "ViennaRNA/datastructures/array.h"
#include "ViennaRNA/utils/basic.h"

#include "ViennaRNA/data/transform.h"

typedef struct {
  vrna_array(double)  threshold_source;
  vrna_array(double)  threshold_target;
  unsigned int        options;
  double              v_min;
  double              v_max;
} bin_transform_param_t;


/*
 #################################
 # PRIVATE FUNCTION DECLARATIONS #
 #################################
 */
PRIVATE void
transform_bin_option_free(vrna_data_lin_trans_opt_t options);


PRIVATE double
transform_bin(double                    v,
              vrna_data_lin_trans_opt_t options);


/*
 #################################
 # BEGIN OF FUNCTION DEFINITIONS #
 #################################
 */
PUBLIC vrna_data_lin_trans_f
vrna_data_transform_method_bin(double                         (*thresholds)[2],
                               unsigned int                   thresholds_num,
                               double                         oolb_value,
                               double                         ooub_value,
                               unsigned int                   options,
                               vrna_data_lin_trans_opt_t      *transform_options_p,
                               vrna_data_lin_trans_opt_free_f *transform_options_free)
{
  vrna_data_lin_trans_f  cb = NULL;

  if ((thresholds) &&
      (thresholds_num > 1) &&
      (transform_options_p) &&
      (transform_options_free)) {
    bin_transform_param_t *o = (bin_transform_param_t *)vrna_alloc(sizeof(bin_transform_param_t));

    o->options  = options;
    o->v_min    = oolb_value;
    o->v_max    = ooub_value;

    vrna_array_init_size(o->threshold_source, thresholds_num);
    vrna_array_init_size(o->threshold_target, thresholds_num);

    for (size_t i = 0; i < thresholds_num; i++) {
      vrna_array_append(o->threshold_source, thresholds[i][0]);
      vrna_array_append(o->threshold_target, thresholds[i][1]);
    }

    cb                      = transform_bin;
    *transform_options_p    = (vrna_data_lin_trans_opt_t)o;
    *transform_options_free = transform_bin_option_free;
  }

  return cb;
}


/*
 #####################################
 # BEGIN OF STATIC HELPER FUNCTIONS  #
 #####################################
 */
PRIVATE void
transform_bin_option_free(vrna_data_lin_trans_opt_t options)
{
  bin_transform_param_t *o = (bin_transform_param_t *)options;

  vrna_array_free(o->threshold_source);
  vrna_array_free(o->threshold_target);
  free(o);
}


PRIVATE double
transform_bin(double                    v,
              vrna_data_lin_trans_opt_t options)
{
  bin_transform_param_t *o = (bin_transform_param_t *)options;

  if (v < o->threshold_source[0])
    return (o->options & VRNA_TRANSFORM_BIN_OPTION_MAP_OUTOF_LOWERBOUND) ? o->threshold_target[0] : o->v_min;

  size_t  start, end, mid, last;

  last  = vrna_array_size(o->threshold_source) - 1;
  start = 1;
  end   = last;

  /* binary search for bin the value will be mapped to */
  while (start <= end) {
    mid = start + (end - start) / 2;

    if ((v >= o->threshold_source[mid - 1]) &&
        (v < o->threshold_source[mid])) {
      if (o->options & VRNA_TRANSFORM_BIN_OPTION_PROJECT) {
        double diff_source = o->threshold_source[mid] - o->threshold_source[mid - 1];
        double diff_target = o->threshold_target[mid] - o->threshold_target[mid - 1];
        return  (v - o->threshold_source[mid - 1]) /
                diff_source *
                diff_target +
                o->threshold_target[mid - 1];
      } else {
        return o->threshold_target[mid - 1];
      }
    }

    if (v >= o->threshold_source[mid])
      start = mid + 1;
    else
      end = mid - 1;
  }

  if (v > o->threshold_source[last])
    return (o->options & VRNA_TRANSFORM_BIN_OPTION_MAP_OUTOF_UPPERBOUND) ? o->threshold_target[last] : o->v_max;
  else
    return o->threshold_target[last];
}
