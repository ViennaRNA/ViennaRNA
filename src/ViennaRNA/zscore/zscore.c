#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdio.h>
#include <stdlib.h>

#include <svm.h>

#ifdef __GNUC__
# define INLINE inline
#else
# define INLINE
#endif

#include "ViennaRNA/utils/basic.h"
#include "ViennaRNA/utils/svm.h"
#include "ViennaRNA/zscore/basic.h"

#include "ViennaRNA/intern/zscore_dat.h"


PRIVATE INLINE double
get_zscore(vrna_fold_compound_t *fc,
           int                  i,
           int                  j,
           int                  e,
           double               *avg,
           double               *sd);


PUBLIC int
vrna_zsc_filter_init(vrna_fold_compound_t *fc,
                     double               min_z,
                     unsigned int         options)
{
  if (fc) {
    vrna_zsc_filter_free(fc);

    fc->zscore_data                   = (vrna_zsc_dat_t)vrna_alloc(sizeof(struct vrna_zsc_dat_s));
    fc->zscore_data->filter_on        = (options & VRNA_ZSCORE_FILTER_ON) ? 1 : 0;
    fc->zscore_data->pre_filter       = (options & VRNA_ZSCORE_PRE_FILTER) ? 1 : 0;
    fc->zscore_data->report_subsumed  = (options & VRNA_ZSCORE_REPORT_SUBSUMED) ? 1 : 0;
    fc->zscore_data->min_z            = min_z;
    fc->zscore_data->avg_model        = svm_load_model_string(avg_model_string);
    fc->zscore_data->sd_model         = svm_load_model_string(sd_model_string);

    if (fc->zscore_data->pre_filter)
      fc->zscore_data->current_z = (double *)vrna_alloc(sizeof(double) * (fc->window_size + 2));
    else
      fc->zscore_data->current_z = NULL;

    fc->zscore_data->current_i = 0;

    return 1;
  }

  return 0;
}


PUBLIC int
vrna_zsc_filter_update(vrna_fold_compound_t *fc,
                       double               min_z,
                       unsigned int         options)
{
  if (fc) {
    if (!fc->zscore_data) {
      vrna_zsc_filter_init(fc, min_z, options);
    } else {
      fc->zscore_data->min_z = min_z;

      if (!(options & VRNA_ZSCORE_OPTIONS_NONE)) {
        fc->zscore_data->filter_on        = (options & VRNA_ZSCORE_FILTER_ON) ? 1 : 0;
        fc->zscore_data->pre_filter       = (options & VRNA_ZSCORE_PRE_FILTER) ? 1 : 0;
        fc->zscore_data->report_subsumed  = (options & VRNA_ZSCORE_REPORT_SUBSUMED) ? 1 : 0;
      }

      if (fc->zscore_data->pre_filter) {
        if (fc->zscore_data->current_z) {
          fc->zscore_data->current_z += fc->zscore_data->current_i;
          free(fc->zscore_data->current_z);
        }

        fc->zscore_data->current_z  = (double *)vrna_alloc(sizeof(double) * (fc->window_size + 2));
        fc->zscore_data->current_i  = 0;
      } else if (fc->zscore_data->current_z) {
        fc->zscore_data->current_z += fc->zscore_data->current_i;
        free(fc->zscore_data->current_z);
        fc->zscore_data->current_z  = NULL;
        fc->zscore_data->current_i  = 0;
      }
    }

    return 1;
  }

  return 0;
}


PUBLIC void
vrna_zsc_filter_free(vrna_fold_compound_t *fc)
{
  if ((fc) && (fc->zscore_data)) {
    vrna_zsc_dat_t zsc_data = fc->zscore_data;

    zsc_data->current_z += zsc_data->current_i;
    free(zsc_data->current_z);
    svm_free_model_content(zsc_data->avg_model);
    svm_free_model_content(zsc_data->sd_model);
    free(zsc_data);

    fc->zscore_data = NULL;
  }
}


PUBLIC int
vrna_zsc_filter_on(vrna_fold_compound_t *fc)
{
  if ((fc) &&
      (fc->zscore_data) &&
      (fc->zscore_data->filter_on))
    return 1;

  return 0;
}


PUBLIC double
vrna_zsc_filter_threshold(vrna_fold_compound_t *fc)
{
  if ((fc) &&
      (fc->zscore_data) &&
      (fc->zscore_data->filter_on))
    return fc->zscore_data->min_z;

  return (double)INF;
}


PUBLIC double
vrna_zsc_compute(vrna_fold_compound_t *fc,
                 unsigned int         i,
                 unsigned int         j,
                 int                  e)
{
  if ((fc) &&
      (fc->zscore_data) &&
      (fc->zscore_data->filter_on))
    return get_zscore(fc, i, j, e, NULL, NULL);

  return (double)INF;
}


PUBLIC double
vrna_zsc_compute_raw(vrna_fold_compound_t *fc,
                     unsigned int         i,
                     unsigned int         j,
                     int                  e,
                     double               *avg,
                     double               *sd)
{
  if ((fc) &&
      (fc->zscore_data) &&
      (fc->zscore_data->filter_on))
    return get_zscore(fc, i, j, e, avg, sd);

  return (double)INF;
}


PRIVATE INLINE double
get_zscore(vrna_fold_compound_t *fc,
           int                  i,
           int                  j,
           int                  e,
           double               *avg,
           double               *sd)
{
  short           *S;
  int             info_avg, start, end, dangle_model, length;
  double          average_free_energy;
  double          sd_free_energy;
  double          z;
  vrna_zsc_dat_t  d;

  length        = fc->length;
  S             = fc->sequence_encoding2;
  dangle_model  = fc->params->model_details.dangles;
  z             = (double)INF;
  d             = fc->zscore_data;

  if (avg)
    *avg = (double)INF;

  if (sd)
    *sd = (double)INF;

  start = (dangle_model) ? MAX2(1, i - 1) : i;
  end   = (dangle_model) ? MIN2(length, j + 1) : j;

  int *AUGC = get_seq_composition(S, start, end, length);

  /*\svm*/
  average_free_energy = avg_regression(AUGC[0],
                                       AUGC[1],
                                       AUGC[2],
                                       AUGC[3],
                                       AUGC[4],
                                       d->avg_model,
                                       &info_avg);

  if (info_avg == 0) {
    double  min_sd      = minimal_sd(AUGC[0], AUGC[1], AUGC[2], AUGC[3], AUGC[4]);
    double  difference  = ((double)e / 100.) - average_free_energy;

    if (difference - (d->min_z * min_sd) <= 0.0001) {
      sd_free_energy  = sd_regression(AUGC[0], AUGC[1], AUGC[2], AUGC[3], AUGC[4], d->sd_model);
      z               = difference / sd_free_energy;
      if (avg)
        *avg = average_free_energy;

      if (sd)
        *sd = sd_free_energy;
    }
  }

  free(AUGC);

  return z;
}
