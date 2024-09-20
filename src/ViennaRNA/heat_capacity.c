#include  <stdio.h>
#include  <stdlib.h>
#include  <math.h>

#include  "ViennaRNA/utils/basic.h"
#include  "ViennaRNA/params/constants.h"
#include  "ViennaRNA/mfe/global.h"
#include  "ViennaRNA/partfunc/global.h"
#include  "ViennaRNA/heat_capacity.h"

#define MAXWIDTH  201

struct data_collector {
  struct vrna_heat_capacity_s *data;
  size_t                      num_entries;
  size_t                      allocated_memory;
};


PRIVATE float
ddiff(float f[],
      float h,
      int   m);


PRIVATE void
store_results_cb(float  t,
                 float  hc,
                 void   *data);


PUBLIC struct vrna_heat_capacity_s *
vrna_heat_capacity_simple(const char    *sequence,
                          float         T_min,
                          float         T_max,
                          float         h,
                          unsigned int  m)
{
  vrna_fold_compound_t        *fc;
  struct vrna_heat_capacity_s *result;

  result = NULL;

  if (sequence) {
    fc = vrna_fold_compound(sequence,
                            NULL,
                            VRNA_OPTION_DEFAULT);

    result = vrna_heat_capacity(fc, T_min, T_max, h, m);

    vrna_fold_compound_free(fc);
  }

  return result;
}


PUBLIC struct vrna_heat_capacity_s *
vrna_heat_capacity(vrna_fold_compound_t *fc,
                   float                T_min,
                   float                T_max,
                   float                h,
                   unsigned int         m)
{
  struct vrna_heat_capacity_s *results = NULL;

  if (fc) {
    struct data_collector d;

    d.num_entries       = 0;
    d.allocated_memory  = 127;
    d.data              =
      (struct vrna_heat_capacity_s *)vrna_alloc(sizeof(struct vrna_heat_capacity_s) *
                                                d.allocated_memory);

    (void)vrna_heat_capacity_cb(fc,
                                T_min,
                                T_max,
                                h,
                                m,
                                store_results_cb,
                                (void *)&d);

    results = (struct vrna_heat_capacity_s *)vrna_realloc(d.data,
                                                          sizeof(struct vrna_heat_capacity_s) *
                                                          (d.num_entries + 1));
    results[d.num_entries].temperature    = -K0 - 1.;
    results[d.num_entries].heat_capacity  = -K0 - 1.;
  }

  return results;
}


PUBLIC int
vrna_heat_capacity_cb(vrna_fold_compound_t        *fc,
                      float                       T_min,
                      float                       T_max,
                      float                       h,
                      unsigned int                m,
                      vrna_heat_capacity_f cb,
                      void                        *data)
{
  unsigned int  i, n;
  int           ret;
  float         hc, F[MAXWIDTH];
  double        min_en;
  vrna_md_t     md, md_init;

  ret = 0;

  if ((fc) && (cb)) {
    /* sanity checks first */
    if (m > 100)
      m = 100;
    else if (m == 0)
      m = 1;

    if (T_min > T_max) {
      hc    = T_min;
      T_min = T_max;
      T_max = hc;
    }

    if (T_min <= -K0)
      T_min = -K0;

    if (h > (T_max - T_min))
      h = T_max - T_min;

    /* now for the actual algorithm */
    n       = fc->length;
    md_init = md = fc->params->model_details;

    /* required for vrna_exp_param_rescale() in subsequent calls */
    md.sfact        = 1.;
    md.backtrack    = 0;
    md.compute_bpp  = 0;

    md.temperature = T_min - m * h;
    vrna_params_reset(fc, &md);

    min_en = (double)vrna_mfe(fc, NULL);

    vrna_exp_params_rescale(fc, &min_en);

    for (i = 0; i < 2 * m + 1; i++) {
      F[i] = vrna_pf(fc, NULL);
      /* increase temperature */
      md.temperature += h;
      /* reset all energy parameters according to temperature changes */
      vrna_params_reset(fc, &md);

      min_en = F[i] + h * 0.00727 * n;

      vrna_exp_params_rescale(fc, &min_en);
    }

    while (md.temperature <= (T_max + m * h + h)) {
      hc = -ddiff(F, h, m) * (md.temperature + K0 - m * h - h);

      /* return results */
      cb((md.temperature - (float)m * h - h), hc, data);

      for (i = 0; i < 2 * m; i++)
        F[i] = F[i + 1];

      F[2 * m] = vrna_pf(fc, NULL);

      /*       printf("%g\n", F[2*m]);*/
      md.temperature += h;

      vrna_params_reset(fc, &md);

      min_en = F[i] + h * 0.00727 * n;

      vrna_exp_params_rescale(fc, &min_en);
    }

    /* restore original state of (the model of) the fold_compound */
    vrna_params_reset(fc, &md_init);

    ret = 1;
  }

  return ret;
}


PRIVATE float
ddiff(float f[],
      float h,
      int   m)
{
  int   i;
  float fp, A, B;

  A = (float)(m * (m + 1) * (2 * m + 1) / 3);                                     /* 2*sum(x^2) */
  B = (float)(m * (m + 1) * (2 * m + 1)) * (float)(3 * m * m + 3 * m - 1) / 15.;  /* 2*sum(x^4) */

  fp = 0.;
  for (i = 0; i < 2 * m + 1; i++)
    fp += f[i] * (A - (float)((2 * m + 1) * (i - m) * (i - m)));

  fp /= ((A * A - B * ((float)(2 * m + 1))) * h * h / 2.);
  return (float)fp;
}


PRIVATE void
store_results_cb(float  t,
                 float  hc,
                 void   *data)
{
  struct data_collector *d = (struct data_collector *)data;

  /* check whether we need to resize the collector size */
  if (d->num_entries == d->allocated_memory) {
    d->allocated_memory *= 1.4;
    d->data             = (struct vrna_heat_capacity_s *)vrna_realloc(d->data,
                                                                      sizeof(struct
                                                                             vrna_heat_capacity_s) *
                                                                      d->allocated_memory);
  }

  d->data[d->num_entries].temperature   = t;
  d->data[d->num_entries].heat_capacity = hc;
  d->num_entries++;
}
