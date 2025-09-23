/* SHAPE reactivity data handling */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <string.h>
#include <limits.h>

#include "ViennaRNA/params/default.h"
#include "ViennaRNA/params/constants.h" /* defines MINPSCORE */
#include "ViennaRNA/datastructures/array.h"
#include "ViennaRNA/fold_vars.h"
#include "ViennaRNA/utils/basic.h"
#include "ViennaRNA/utils/strings.h"
#include "ViennaRNA/sequences/alignments.h"
#include "ViennaRNA/utils/log.h"
#include "ViennaRNA/io/utils.h"
#include "ViennaRNA/io/file_formats.h"
#include "ViennaRNA/params/basic.h"
#include "ViennaRNA/constraints/soft.h"
#include "ViennaRNA/probing/basic.h"
#include "ViennaRNA/probing/transform.h"
#include "ViennaRNA/probing/strategies.h"


#define gaussian(u) (1 / (sqrt(2 * PI)) * exp(-u * u / 2))


struct vrna_probing_data_s {
  unsigned int              method;
  vrna_array(double *)                  data_linear;   /* actual data */
  vrna_array(double)                    data_linear_weight;    /* weight for each data set */
  vrna_array(vrna_probing_strategy_f)   cbs_linear;
  vrna_array(void *)                    cbs_linear_options;
  vrna_array(vrna_auxdata_free_f)       cbs_linear_options_free;

  vrna_array(double)    params1;
  vrna_array(double)    params2;
  vrna_array(double *)  reactivities;
  vrna_array(double *)  transformeds;
  vrna_array(double *)  datas1;
  vrna_array(double *)  datas2;
};


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
PRIVATE void
nullify_probing_data_s(struct vrna_probing_data_s *data);


PRIVATE int
apply_probing_data(vrna_fold_compound_t        *fc,
                   struct vrna_probing_data_s  *data);

PRIVATE int
apply_Deigan2009_method(vrna_fold_compound_t        *fc,
                        struct vrna_probing_data_s  *data);


PRIVATE int
apply_Zarringhalam2012_method(vrna_fold_compound_t        *fc,
                              struct vrna_probing_data_s  *data);


PRIVATE int
apply_Washietl2012_method(vrna_fold_compound_t        *fc,
                          struct vrna_probing_data_s  *data);


PRIVATE int
apply_Eddy2014_method(vrna_fold_compound_t        *fc,
                      struct vrna_probing_data_s  *data);


PRIVATE double
get_msa_weight(const double **datasets,
               unsigned int data_size);


PRIVATE FLT_OR_DBL
conversion_deigan(double  reactivity,
                  double  m,
                  double  b);


PRIVATE FLT_OR_DBL
conversion_zarringhalam_up(double       beta,
                           double       *pr,
                           unsigned int i);


PRIVATE FLT_OR_DBL
conversion_zarringhalam_bp(double       beta,
                           double       *pr,
                           unsigned int i,
                           unsigned int j);


PRIVATE FLT_OR_DBL
conversion_eddy_up(double       kT,
                   double       *data,
                   unsigned int i);


PRIVATE FLT_OR_DBL
conversion_eddy_bp(double       kT,
                   double       *data,
                   unsigned int i,
                   unsigned int j);


/* PDF of x using Gaussian KDE
 * n is data number
 * h is bandwidth
 */
PRIVATE FLT_OR_DBL
gaussian_kde_pdf(double       x,
                 unsigned int n,
                 float        h,
                 const double *data);


PRIVATE FLT_OR_DBL
exp_pdf(double  x,
        double  lambda) VRNA_UNUSED;


/* We use same format as scitpy */
PRIVATE FLT_OR_DBL
gev_pdf(double  x,
        double  c,
        double  loc,
        double  scale) VRNA_UNUSED;


/*
 * Bandwidth for univariate KDE with Scott factor as in scipy
 * bandwidth = Scott facter * std with ddof = 1
 */
PRIVATE FLT_OR_DBL
bandwidth(unsigned int  n,
          const double  *data);


PRIVATE void
sc_parse_parameters(const char  *string,
                    char        c1,
                    char        c2,
                    float       *v1,
                    float       *v2);


PRIVATE vrna_probing_strategy_f
get_cb_stack_default(void);


PRIVATE void *
get_cb_stack_options_default(void);


PRIVATE vrna_auxdata_free_f
get_cb_stack_options_free_default(void);


PRIVATE vrna_probing_strategy_f
get_cb_up_default(void);


PRIVATE void *
get_cb_up_options_default(void);


PRIVATE vrna_auxdata_free_f
get_cb_up_options_free_default(void);


/*
 #################################
 # BEGIN OF FUNCTION DEFINITIONS #
 #################################
 */
PUBLIC int
vrna_sc_probing(vrna_fold_compound_t  *fc,
                vrna_probing_data_t   data)
{
  int ret = 0;

  if ((fc) && (data))
    return apply_probing_data(fc, data);

  return ret;
}


PUBLIC vrna_probing_data_t
vrna_probing_data_linear(const double              *data,
                         unsigned int              data_length,
                         double                    data_weight,
                         vrna_probing_strategy_f   strategy_cb,
                         void                      *strategy_cb_options,
                         vrna_auxdata_free_f       strategy_cb_options_free)
{
  struct vrna_probing_data_s  *d = NULL;

  if (data) {
    return vrna_probing_data_linear_multi(&data,
                                          1,
                                          &data_length,
                                          &data_weight,
                                          &strategy_cb,
                                          &strategy_cb_options,
                                          &strategy_cb_options_free,
                                          VRNA_PROBING_DATA_DEFAULT);
  }

  return d;
}


PUBLIC vrna_probing_data_t
vrna_probing_data_linear_multi(const double              **data,
                               unsigned int              data_size,
                               const unsigned int        *data_lengths,
                               const double              *data_weights,
                               vrna_probing_strategy_f   *strategy_cbs,
                               void                      **strategy_cbs_options,
                               vrna_auxdata_free_f       *strategy_cbs_options_free,
                               unsigned int              options)
{
  double                      weight = 1.0;
  struct vrna_probing_data_s  *d = NULL;
  vrna_probing_strategy_f   cb;
  void                      *cb_options;
  vrna_auxdata_free_f       cb_options_free;

  if ((data) &&
      (data_lengths) &&
      (data_size)) {

    d = (struct vrna_probing_data_s *)vrna_alloc(sizeof(struct vrna_probing_data_s));

    nullify_probing_data_s(d);

    vrna_array_init_size(d->data_linear, data_size);
    vrna_array_init_size(d->data_linear_weight,  data_size);
    vrna_array_init_size(d->cbs_linear, data_size);
    vrna_array_init_size(d->cbs_linear_options, data_size);
    vrna_array_init_size(d->cbs_linear_options_free, data_size);

    if (options & VRNA_PROBING_DATA_SINGLE_STRATEGY) {
      if (strategy_cbs) {
        cb              = strategy_cbs[0];
        cb_options      = (strategy_cbs_options) ? strategy_cbs_options[0] : NULL;
        cb_options_free = (strategy_cbs_options_free) ? strategy_cbs_options_free[0] : NULL;
      } else {
        cb              = get_cb_stack_default();
        cb_options      = get_cb_stack_options_default();
        cb_options_free = get_cb_stack_options_free_default();
      }
    }

    if (options & VRNA_PROBING_DATA_SINGLE_WEIGHT) {
      if (data_weights) {
        weight = data_weights[0];
      } else {
        weight = get_msa_weight(data, data_size);
      }
    } else if (!(data_weights)) {
      weight = get_msa_weight(data, data_size);
    }

    for (size_t i = 0; i < data_size; i++) {
      if ((data[i]) &&
          (data_lengths[i])) {
        /* init and store raw probing data */
        vrna_array(double)  a;
        vrna_array_init_size(a, data_lengths[i] + 1);
        for (size_t j = 0; j <= data_lengths[i]; j++)
          vrna_array_append(a, data[i][j]);

        vrna_array_append(d->data_linear, a);

        /* store weight for this data set */
        if (options & VRNA_PROBING_DATA_SINGLE_WEIGHT) {
          vrna_array_append(d->data_linear_weight, weight);
        } else {
          vrna_array_append(d->data_linear_weight, (data_weights) ? data_weights[i] : weight);
        }

        /* set corresponding conversion strategy */
        if (options & VRNA_PROBING_DATA_SINGLE_STRATEGY) {
          vrna_array_append(d->cbs_linear, cb);
          vrna_array_append(d->cbs_linear_options, cb_options);
          vrna_array_append(d->cbs_linear_options_free, cb_options_free);
        } else if ((strategy_cbs) &&
                   (strategy_cbs[i])) {
          vrna_array_append(d->cbs_linear, strategy_cbs[i]);

          if (strategy_cbs_options)
            vrna_array_append(d->cbs_linear_options, strategy_cbs_options[i]);
          else
            vrna_array_append(d->cbs_linear_options, NULL);

          if (strategy_cbs_options_free)
            vrna_array_append(d->cbs_linear_options_free, strategy_cbs_options_free[i]);
          else
            vrna_array_append(d->cbs_linear_options_free, NULL);
        } else {
          /* use default strategy */
          vrna_array_append(d->cbs_linear, get_cb_stack_default());
          vrna_array_append(d->cbs_linear_options, get_cb_stack_options_default());
          vrna_array_append(d->cbs_linear_options_free, get_cb_stack_options_free_default());
        }
      } else {
        vrna_array_append(d->data_linear, NULL);
        vrna_array_append(d->data_linear_weight, 0.);
        vrna_array_append(d->cbs_linear, NULL);
        vrna_array_append(d->cbs_linear_options, NULL);
        vrna_array_append(d->cbs_linear_options_free, NULL);
      }
    }
  }

  return d;
}


PUBLIC struct vrna_probing_data_s *
vrna_probing_data_Eddy2014_2(const double *reactivities,
                             unsigned int n,
                             const double *unpaired_data,
                             unsigned int unpaired_len,
                             const double *paired_data,
                             unsigned int paired_len,
                             double       (*trans) (double, void*),
                             void         *options)
{
  struct vrna_probing_data_s *d = NULL;

  if (reactivities)
    d = vrna_probing_data_Eddy2014_2_comparative(&reactivities,
                                                 &n,
                                                 1,
                                                 &unpaired_data,
                                                 &unpaired_len,
                                                 &paired_data,
                                                 &paired_len,
                                                 VRNA_PROBING_METHOD_MULTI_PARAMS_0,
                                                 trans,
                                                 options);

  return d;
}


PUBLIC struct vrna_probing_data_s *
vrna_probing_data_Eddy2014_2_comparative(const double **reactivities,
                                         unsigned int *n,
                                         unsigned int n_seq,
                                         const double **unpaired_datas,
                                         unsigned int *unpaired_lens,
                                         const double **paired_datas,
                                         unsigned int *paired_lens,
                                         unsigned int multi_params,
                                         double       (*trans) (double, void*),
                                         void         *options)
{
  struct vrna_probing_data_s  *d = NULL;
  double                      unpaired_h, paired_h;

  if ((reactivities) &&
      (unpaired_datas) &&
      (paired_datas) &&
      (unpaired_datas[0]) &&
      (paired_datas[0])) {
    d = (struct vrna_probing_data_s *)vrna_alloc(sizeof(struct vrna_probing_data_s));

    nullify_probing_data_s(d);

    d->method = VRNA_PROBING_METHOD_EDDY2014_2;
    vrna_array_init(d->params1);
    vrna_array_init(d->params2);
    vrna_array_init_size(d->reactivities, n_seq);
    vrna_array_init_size(d->transformeds, n_seq);
    vrna_array_init_size(d->datas1, n_seq);
    vrna_array_init_size(d->datas2, n_seq);

    if (trans == NULL)
      trans = vrna_reactivity_trans_default(VRNA_PROBING_METHOD_EDDY2014_2);

    /* prepare first probabilities */
    unpaired_h  = bandwidth(unpaired_lens[0],
        vrna_reactivity_transform(unpaired_lens[0], unpaired_datas[0], trans, options));
    paired_h    = bandwidth(paired_lens[0],
        vrna_reactivity_transform(paired_lens[0], paired_datas[0], trans, options));

    for (unsigned int i = 0; i < n_seq; i++) {
      if (reactivities[i]) {
        /* init and store reactivity data */
        vrna_array(FLT_OR_DBL)  a;
        vrna_array_init_size(a, n[i] + 1);
        for (unsigned int j = 0; j <= n[i]; j++)
          vrna_array_append(a, (FLT_OR_DBL)reactivities[i][j]);

        vrna_array_append(d->reactivities, a);
        vrna_array_append(d->transformeds, vrna_reactivity_transform(n[i], reactivities[i], trans, options));

        /* kernel-density probability computations */
        vrna_array(FLT_OR_DBL) unpaired;
        vrna_array(FLT_OR_DBL) paired;

        vrna_array_init_size(unpaired, n[i] + 1);
        vrna_array_init_size(paired, n[i] + 1);

        /* Transform unpaired and paired reactivity for distribution */
        double *unpaired_tr = vrna_reactivity_transform(unpaired_lens[i], unpaired_datas[i], trans, options);
        double *paired_tr = vrna_reactivity_transform(paired_lens[i], paired_datas[i], trans, options);


        /* Compute bandwidth? */
        if (multi_params & VRNA_PROBING_METHOD_MULTI_PARAMS_1)
          unpaired_h = bandwidth(unpaired_lens[i], unpaired_tr);

        if (multi_params & VRNA_PROBING_METHOD_MULTI_PARAMS_2)
          paired_h = bandwidth(paired_lens[i], paired_tr);

        /* convert and add */
        vrna_array_append(unpaired, 0);
        vrna_array_append(paired, 0.);

        for (unsigned int j = 1; j <= n[i]; ++j) {
          if (d->transformeds[i][j] == VRNA_REACTIVITY_MISSING) {
            vrna_array_append(unpaired, 0);
            vrna_array_append(paired, 0);
          } else {
            vrna_array_append(unpaired,
                              log(gaussian_kde_pdf(d->transformeds[i][j], unpaired_lens[i], unpaired_h,
                                                   unpaired_tr)));
            vrna_array_append(paired,
                              log(gaussian_kde_pdf(d->transformeds[i][j], paired_lens[i], paired_h,
                                                   paired_tr)));
          }
        }

        vrna_array_free(unpaired_tr);
        vrna_array_free(paired_tr);

        vrna_array_append(d->datas1, unpaired);
        vrna_array_append(d->datas2, paired);

        for (unsigned int j = 0; j <= n[i]; j++)
          printf("Before: %f. After: %f. Unpaired: %f. Paired: %f\n", d->reactivities[i][j], d->transformeds[i][j], d->datas1[i][j], d->datas2[i][j]);
      } else {
        vrna_array_append(d->reactivities, NULL);
        vrna_array_append(d->datas1, NULL);
        vrna_array_append(d->datas2, NULL);
      }
    }
  }

  return d;
}


PUBLIC void
vrna_probing_data_free(struct vrna_probing_data_s *d)
{
  if (d) {
    /* free all reactivity data */
    if (d->reactivities)
      for (unsigned int i = 0; i < vrna_array_size(d->reactivities); i++)
        vrna_array_free(d->reactivities[i]);

    vrna_array_free(d->reactivities);

    if (d->transformeds)
      for (unsigned int i = 0; i < vrna_array_size(d->transformeds); i++)
          printf("%u, %ld\n", i, d->transformeds[i]);
//        vrna_array_free(d->transformeds[i]);

    vrna_array_free(d->transformeds);

    /* free parameters */
    vrna_array_free(d->params1);
    vrna_array_free(d->params2);

    /* free auxiliary data */
    if (d->datas1) {
      for (unsigned int i = 0; i < vrna_array_size(d->datas1); i++)
        vrna_array_free(d->datas1[i]);

      vrna_array_free(d->datas1);
    }

    if (d->datas2) {
      for (unsigned int i = 0; i < vrna_array_size(d->datas2); i++)
        vrna_array_free(d->datas2[i]);

      vrna_array_free(d->datas2);
    }

    free(d);
  }
}


PUBLIC int
vrna_sc_SHAPE_to_pr(const char  *shape_conversion,
                    double      *values,
                    int         length,
                    double      default_value)
{
  int *indices;
  int i, j;
  int index;
  int ret = 1;

  if (!shape_conversion || !(*shape_conversion) || length <= 0)
    return 0;

  if (*shape_conversion == 'S')
    return 1;

  indices = vrna_alloc(sizeof(int) * (length + 1));
  for (i = 1, j = 0; i <= length; ++i) {
    if (values[i] == VRNA_REACTIVITY_MISSING)
      values[i] = default_value;
    else
      indices[j++] = i;
  }

  if (*shape_conversion == 'M') {
    double  max;
    double  map_info[4][2] = { {
                                 0.25, 0.35
                               },
                               {
                                 0.30, 0.55
                               },
                               {
                                 0.70, 0.85
                               },
                               {
                                 0, 1
                               } };

    max = values[1];
    for (i = 2; i <= length; ++i)
      max = MAX2(max, values[i]);
    map_info[3][0] = max;

    for (i = 0; indices[i]; ++i) {
      double  lower_source  = 0;
      double  lower_target  = 0;

      index = indices[i];

      if (values[index] == 0)
        continue;

      for (j = 0; j < 4; ++j) {
        if (values[index] > lower_source && values[index] <= map_info[j][0]) {
          double  diff_source = map_info[j][0] - lower_source;
          double  diff_target = map_info[j][1] - lower_target;
          values[index] = (values[index] - lower_source) / diff_source * diff_target + lower_target;
          break;
        }

        lower_source  = map_info[j][0];
        lower_target  = map_info[j][1];
      }
    }
  } else if (*shape_conversion == 'C') {
    float cutoff = 0.25;
    int   i;

    sscanf(shape_conversion + 1, "%f", &cutoff);

    for (i = 0; indices[i]; ++i) {
      index         = indices[i];
      values[index] = values[index] < cutoff ? 0 : 1;
    }
  } else if (*shape_conversion == 'L' || *shape_conversion == 'O') {
    int   i;
    float slope     = (*shape_conversion == 'L') ? 0.68 : 1.6;
    float intercept = (*shape_conversion == 'L') ? 0.2 : -2.29;

    sc_parse_parameters(shape_conversion + 1, 's', 'i', &slope, &intercept);

    for (i = 0; indices[i]; ++i) {
      double v;
      index = indices[i];









      v             = (*shape_conversion == 'L') ? values[index] : log(values[index]);
      values[index] = MAX2(MIN2((v - intercept) / slope, 1), 0);
    }
  } else {
    ret = 0;
  }

  free(indices);

  return ret;
}




/*
 #####################################
 # BEGIN OF STATIC HELPER FUNCTIONS  #
 #####################################
 */
PRIVATE int
apply_probing_data(vrna_fold_compound_t        *fc,
                   struct vrna_probing_data_s  *data)
{
  unsigned int  i, j, s, n, **a2s;
  int           ret, num_data;
  FLT_OR_DBL    *vs, **cvs, weight;
  double        *e;

  ret       = 0;
  num_data  = 0;
  n         = fc->length;
#if 1
  if ((data->data_linear) &&
      (vrna_array_size(data->data_linear) > 0)) {
    ret = 1;

    switch (fc->type) {
      case VRNA_FC_TYPE_SINGLE:
        if (vrna_array_size(data->data_linear) > 1)
          vrna_log_warning("Multiple probing data (%u) for single sequence",
                           vrna_array_size(data->data_linear));

        if (vrna_array_size(data->data_linear[0]) > 0) {
          if ((size_t)vrna_array_size(data->data_linear[0]) != (size_t)(fc->length + 1))
            vrna_log_warning("Length of probing data (%u) doesn't match length of sequence (%u)",
                             vrna_array_size(data->data_linear[0]),
                             fc->length);

          /* convert for nucleotides within a stack */
          e = data->cbs_linear[0](data->data_linear[0],
                                 vrna_array_size(data->data_linear[0]),
                                 VRNA_PROBING_DATA_LINEAR_TARGET_STACK,
                                 data->cbs_linear_options[0]);

          if (e) {
            n = MIN2(n, vrna_array_size(data->data_linear[0]) - 1);

            for (i = 1; i <= n; ++i)
              ret &= vrna_sc_add_stack(fc,
                                       i,
                                       e[i] * data->data_linear_weight[0],
                                       VRNA_OPTION_DEFAULT);

            free(e);
          }


          e = data->cbs_linear[0](data->data_linear[0],
                                 vrna_array_size(data->data_linear[0]),
                                 VRNA_PROBING_DATA_LINEAR_TARGET_UP,
                                 data->cbs_linear_options[0]);

          if (e) {
            n = MIN2(n, vrna_array_size(data->data_linear[0]) - 1);

            for (i = 1; i <= n; ++i) {
              ret &= vrna_sc_add_up(fc,
                                    i,
                                    e[i] * data->data_linear_weight[0],
                                    VRNA_OPTION_DEFAULT);
            }

            free(e);
          }

          e = data->cbs_linear[0](data->data_linear[0],
                                 vrna_array_size(data->data_linear[0]),
                                 VRNA_PROBING_DATA_LINEAR_TARGET_BP,
                                 data->cbs_linear_options[0]);

          if (e) {
            n = MIN2(n, vrna_array_size(data->data_linear[0]) - 1);

            for (i = 1; i <= n; ++i)
              for (j = i + 1; j <= n; ++j) {
                ret &= vrna_sc_add_bp(fc,
                                      i,
                                      j,
                                      (e[i] + e[j]) * data->data_linear_weight[0],
                                      VRNA_OPTION_DEFAULT);
              }

            free(e);
          }

        } else {
          vrna_log_warning("Zero-length probing data vector");
        }

        break;

      case VRNA_FC_TYPE_COMPARATIVE:
        if (vrna_array_size(data->data_linear) != fc->n_seq)
          vrna_log_warning("Number of probing data (%u) doesn't match number of sequences in the alignment (%u)",
                           vrna_array_size(data->data_linear), fc->n_seq);

        for (s = 0; s < vrna_array_size(data->data_linear); s++) {
          if (vrna_array_size(data->data_linear[s]) > 0) {
            if ((size_t)vrna_array_size(data->data_linear[s]) != (size_t)(fc->alignment->gapfree_size[s] + 1))
              vrna_log_warning("Length of probing data (%u) for sequence no. \"%u\" doesn't match length of gap-free sequence (%u)",
                               vrna_array_size(data->data_linear[s]) - 1,
                               s,
                               fc->alignment->gapfree_size[s]);

            e = data->cbs_linear[s](data->data_linear[s],
                                   vrna_array_size(data->data_linear[s]) - 1,
                                   VRNA_PROBING_DATA_LINEAR_TARGET_STACK,
                                   data->cbs_linear_options[s]);

            n = MIN2(vrna_array_size(data->data_linear[s]), fc->alignment->gapfree_size[s]);

            for (i = 1; i <= n; ++i)
              ret &= vrna_sc_add_stack_comparative_seq(fc,
                                                       s,
                                                       i,
                                                       e[i] * data->data_linear_weight[s],
                                                       VRNA_OPTION_DEFAULT);

            free(e);
          }
        }

        break;
    }
  }
#else
  switch (fc->type) {
    case VRNA_FC_TYPE_SINGLE:
      if ((vrna_array_size(data->reactivities) > 0) &&
          (n <= vrna_array_size(data->reactivities[0]))) {
        ret = 1;

        /* first convert the values according to provided slope and intercept values */
        for (i = 1; i <= n; ++i)
          ret &= vrna_sc_add_stack(fc,
                                   i,
                                   conversion_deigan(data->transformeds[0][i],
                                                     data->params1[0],
                                                     data->params2[0]),
                                   VRNA_OPTION_DEFAULT);
      }

      break;

    case VRNA_FC_TYPE_COMPARATIVE:
      if (vrna_array_size(data->reactivities) >= fc->n_seq) {
        a2s       = fc->a2s;
        weight    = get_msa_weight(fc->n_seq, data->reactivities);
        ret       = 1;

        for (s = 0; s < fc->n_seq; s++) {
          if (data->reactivities[s] != NULL) {
            n = a2s[s][fc->length];
            num_data++;

            for (i = 1; i <= n; i++) {
              ret &= vrna_sc_add_stack_comparative_seq(fc,
                                                       s,
                                                       i,
                                                       conversion_deigan(data->transformeds[s][i],
                                                                         data->params1[s],
                                                                         data->params2[s]) *
                                                       weight,
                                                       VRNA_OPTION_DEFAULT);
            }
          }
        }
      }

      if (ret)
        ret = num_data;

      break;
  }
#endif

  return ret;
}
PRIVATE int
apply_Deigan2009_method(vrna_fold_compound_t        *fc,
                        struct vrna_probing_data_s  *data)
{
  unsigned int  i, s, n, **a2s;
  int           ret, num_data;
  FLT_OR_DBL    *vs, **cvs, weight;

  ret       = 0;
  num_data  = 0;
  n         = fc->length;

  switch (fc->type) {
    case VRNA_FC_TYPE_SINGLE:
      if ((vrna_array_size(data->reactivities) > 0) &&
          (n <= vrna_array_size(data->reactivities[0]))) {
        ret = 1;

        /* first convert the values according to provided slope and intercept values */
        for (i = 1; i <= n; ++i)
          ret &= vrna_sc_add_stack(fc,
                                   i,
                                   conversion_deigan(data->transformeds[0][i],
                                                     data->params1[0],
                                                     data->params2[0]),
                                   VRNA_OPTION_DEFAULT);
      }

      break;

    case VRNA_FC_TYPE_COMPARATIVE:
      if (vrna_array_size(data->reactivities) >= fc->n_seq) {
        a2s       = fc->a2s;
        weight    = get_msa_weight((const double **)data->reactivities, fc->n_seq);
        ret       = 1;

        for (s = 0; s < fc->n_seq; s++) {
          if (data->reactivities[s] != NULL) {
            n = a2s[s][fc->length];
            num_data++;

            for (i = 1; i <= n; i++) {
              ret &= vrna_sc_add_stack_comparative_seq(fc,
                                                       s,
                                                       i,
                                                       conversion_deigan(data->transformeds[s][i],
                                                                         data->params1[s],
                                                                         data->params2[s]) *
                                                       weight,
                                                       VRNA_OPTION_DEFAULT);
            }
          }
        }
      }

      if (ret)
        ret = num_data;

      break;
  }

  return ret;
}


PRIVATE int
apply_Zarringhalam2012_method(vrna_fold_compound_t        *fc,
                              struct vrna_probing_data_s  *data)
{
  unsigned int  i, j, n, s, **a2s;
  int           ret, num_data;
  FLT_OR_DBL    *up, **bp, **ups, ***bps, weight;
  double        *pr, b;

  n         = fc->length;
  ret       = 0;
  num_data  = 0;

  switch (fc->type) {
    case VRNA_FC_TYPE_SINGLE:
      if ((vrna_array_size(data->datas1) > 0) &&
          (n <= vrna_array_size(data->datas1[0]))) {
        /*  now, convert probabilities into pseudo free energies for unpaired, and
         *  paired nucleotides
         */
        pr  = data->datas1[0];
        b   = data->params1[0];
        ret = 1;

        /* add the pseudo energies as soft constraints */
        for (i = 1; i <= n; ++i) {
          ret &= vrna_sc_add_up(fc,
                                i,
                                conversion_zarringhalam_up(b, pr, i),
                                VRNA_OPTION_DEFAULT);

          for (j = i + 1; j <= n; ++j)
            ret &= vrna_sc_add_bp(fc,
                                  i,
                                  j,
                                  conversion_zarringhalam_bp(b, pr, i, j),
                                  VRNA_OPTION_DEFAULT);
        }
      }

      break;

    case VRNA_FC_TYPE_COMPARATIVE:
      if (vrna_array_size(data->reactivities) >= fc->n_seq) {
        a2s = fc->a2s;

        /* compute weight for each set of probing data */
        weight = get_msa_weight((const double **)data->reactivities, fc->n_seq);

        ret = 1;

        for (s = 0; s < fc->n_seq; s++) {
          if (data->reactivities[s]) {
            /*  now, convert probabilities into pseudo free energies for unpaired, and
             *  paired nucleotides
             */
            n   = a2s[s][fc->length];
            pr  = data->datas1[s];
            b   = data->params1[s] * weight;
            num_data++;

            /* add the pseudo energies as soft constraints */
            for (i = 1; i <= n; ++i)
              ret &= vrna_sc_add_up_comparative_seq(fc,
                                                    s,
                                                    i,
                                                    conversion_zarringhalam_up(b, pr, i),
                                                    VRNA_OPTION_DEFAULT);

            for (i = 1; i <= n; ++i)
              for (j = i + 1; j <= n; ++j)
                ret &= vrna_sc_add_bp_comparative_seq(fc,
                                                      s,
                                                      i,
                                                      j,
                                                      conversion_zarringhalam_bp(b, pr, i, j),
                                                      VRNA_OPTION_DEFAULT);
          }
        }
      }

      if (ret)
        ret = num_data;

      break;
  }

  return ret;
}


PRIVATE int
apply_Washietl2012_method(vrna_fold_compound_t        *fc,
                          struct vrna_probing_data_s  *data VRNA_UNUSED)
{
  int ret;

  ret = 0;
  switch (fc->type) {
    case VRNA_FC_TYPE_SINGLE:
      break;

    case VRNA_FC_TYPE_COMPARATIVE:
      break;
  }

  return ret;
}


PRIVATE int
apply_Eddy2014_method(vrna_fold_compound_t        *fc,
                      struct vrna_probing_data_s  *data)
{
  unsigned int  i, j, n, s, **a2s;
  int           ret, num_data;
  FLT_OR_DBL    weight;
  double        kT;

  n         = fc->length;
  ret       = 0;
  num_data  = 0;

  kT = GASCONST * ((fc->params)->temperature + K0) / 1000; /* in kcal/mol */

  switch (fc->type) {
    case VRNA_FC_TYPE_SINGLE:
      if ((vrna_array_size(data->reactivities) > 0) &&
          (n <= vrna_array_size(data->reactivities[0])) &&
          (vrna_array_size(data->datas1) > 0) &&
          (n <= vrna_array_size(data->datas1[0])) &&
          (vrna_array_size(data->datas2) > 0) &&
          (n <= vrna_array_size(data->datas2[0]))) {

        ret = 1;

        /* add for unpaired position */
        for (i = 1; i <= n; ++i)
          ret &= vrna_sc_add_up(fc,
                                i,
                                conversion_eddy_up(kT, data->datas1[0], i),
                                VRNA_OPTION_DEFAULT);

        /* add for paired position */
        for (i = 1; i <= n; ++i)
          for (j = i + 1; j <= n; ++j)
            ret &= vrna_sc_add_bp(fc,
                                  i,
                                  j,
                                  conversion_eddy_bp(kT, data->datas2[0], i, j),
                                  VRNA_OPTION_DEFAULT);
      }

      break;

    case VRNA_FC_TYPE_COMPARATIVE:
      if (vrna_array_size(data->reactivities) >= fc->n_seq) {
        a2s = fc->a2s;
        ret = 1;

        /* compute weight for each set of probing data */
        weight = get_msa_weight((const double **)data->reactivities, fc->n_seq);

        kT *= weight;

        for (s = 0; s < fc->n_seq; s++) {
          if (data->reactivities[s]) {
            /*  now, convert probabilities into pseudo free energies for unpaired, and
             *  paired nucleotides
             */
            n = a2s[s][fc->length];
            num_data++;

            /* add for unpaired position */
            for (i = 1; i <= n; ++i)
              ret &= vrna_sc_add_up_comparative_seq(fc,
                                                    s,
                                                    i,
                                                    conversion_eddy_up(kT,
                                                                       data->datas1[s],
                                                                       i),
                                                    VRNA_OPTION_DEFAULT);

            /*
             * add for paired position
             * currently, paired probabilities have to be stored in alignment coordinates
             */
            for (i = 1; i <= n; ++i)
              for (j = i + 1; j <= n; ++j)
                ret &= vrna_sc_add_bp_comparative_seq(fc,
                                                      s,
                                                      i,
                                                      j,
                                                      conversion_eddy_bp(kT,
                                                                         data->datas2[s],
                                                                         i,
                                                                         j),
                                                      VRNA_OPTION_DEFAULT);
          }
        }
      }

      if (ret)
        ret = num_data;

      break;
  }

  return ret;
}


/*
 *  computes the individual weight for each set of probing data
 *  by simply spreading-out the probing data contributions to
 *  all sequences where no probing data is provided for.
 */
PRIVATE double
get_msa_weight(const double **datasets,
               unsigned int data_size)
{
  unsigned int  s, n_data;
  double        weight = 1.0;

  for (n_data = s = 0; s < data_size; s++)
    if (datasets[s] != NULL)
      n_data++;

  if (n_data != 0)
    weight = ((double)data_size / (double)n_data);

  return weight;
}


PRIVATE FLT_OR_DBL
conversion_deigan(double  reactivity,
                  double  m,
                  double  b)
{
  return reactivity == VRNA_REACTIVITY_MISSING ? 0. : (FLT_OR_DBL)(m * log(reactivity + 1) + b);
}


PRIVATE FLT_OR_DBL
conversion_zarringhalam_up(double       beta,
                           double       *pr,
                           unsigned int i)
{
  return beta * fabs(pr[i] - 1.);
}


PRIVATE FLT_OR_DBL
conversion_zarringhalam_bp(double       beta,
                           double       *pr,
                           unsigned int i,
                           unsigned int j)
{
  return beta * (pr[i] + pr[j]);
}


PRIVATE FLT_OR_DBL
conversion_eddy_up(double       kT,
                   double       *data,
                   unsigned int i)
{
  return -kT * data[i];
}


PRIVATE FLT_OR_DBL
conversion_eddy_bp(double       kT,
                   double       *data,
                   unsigned int i,
                   unsigned int j)
{
  return -kT * (data[i] + data[j]);
}


/*
 * ****************
 *      Eddy
 *****************
 */
PRIVATE FLT_OR_DBL
gaussian_kde_pdf(double       x,
                 unsigned int n,
                 float        h,
                 const double *data)
{
  FLT_OR_DBL total;
  int count;

  total = 0.;
  count = 0;
  for (unsigned int i = 0; i < n; i++) {
    if (data[i] != VRNA_REACTIVITY_MISSING) {
      total += gaussian((x - data[i]) / h);
      count++;
    }
  }
  return total / (count * h);
}


PRIVATE FLT_OR_DBL
exp_pdf(double  x,
        double  lambda)
{
  return x < 0 ? 0. : (FLT_OR_DBL)(lambda * exp(-lambda * x));
}


PRIVATE FLT_OR_DBL
gev_pdf(double  x,
        double  c,
        double  loc,
        double  scale)
{
  FLT_OR_DBL s, t;

  s = (FLT_OR_DBL)((x - loc) / scale);
  t = c * s;

  if (c == 0) {
    return exp(-s) * exp(-exp(-s));
  } else if (t < 1) {
    return pow(1 - t, 1 / c - 1) * exp(-pow(1 - t, 1 / c));
  } else {
    return 0;
  }
}


PRIVATE FLT_OR_DBL
bandwidth(unsigned int  n,
          const double  *data)
{
  double factor, mu, std;
  int count;

  mu  = 0.;
  std = 0.;
  count = 0;

  factor = (pow(n, -1. / 5));

  for (unsigned int i = 0; i < n; i++) {
    if (data[i] != VRNA_REACTIVITY_MISSING) {
      mu += data[i];
      count++;
    }
  }
  mu /= count;

  for (unsigned int i = 0; i < n; i++) {
    if (data[i] != VRNA_REACTIVITY_MISSING)
      std += (data[i] - mu) * (data[i] - mu);
  }
  std = sqrt(std / (count - 1));
  return (FLT_OR_DBL)(factor * std);
}


PRIVATE void
sc_parse_parameters(const char  *string,
                    char        c1,
                    char        c2,
                    float       *v1,
                    float       *v2)
{
  char        *fmt;
  const char  warning[] = "SHAPE method parameters not recognized! Using default parameters!";
  int         r;

  assert(c1);
  assert(v1);

  if (!string || !(*string))
    return;

  if (c2 == 0 || v2 == NULL) {
    fmt = vrna_strdup_printf("%c%%f", c1);
    r   = sscanf(string, fmt, v1);

    if (!r)
      vrna_log_warning(warning);

    free(fmt);

    return;
  }

  fmt = vrna_strdup_printf("%c%%f%c%%f", c1, c2);
  r   = sscanf(string, fmt, v1, v2);

  if (r != 2) {
    free(fmt);
    fmt = vrna_strdup_printf("%c%%f", c1);
    r   = sscanf(string, fmt, v1);

    if (!r) {
      free(fmt);
      fmt = vrna_strdup_printf("%c%%f", c2);
      r   = sscanf(string, fmt, v2);

      if (!r)
        vrna_log_warning(warning);
    }
  }

  free(fmt);
}


PRIVATE vrna_probing_strategy_f
get_cb_stack_default(void)
{
  return vrna_probing_strategy_deigan;
}

PRIVATE void *
get_cb_stack_options_default(void)
{

  return NULL;
}

PRIVATE vrna_auxdata_free_f
get_cb_stack_options_free_default(void)
{
  return NULL;
}


PRIVATE vrna_probing_strategy_f
get_cb_up_default(void)
{
  return NULL;
}

PRIVATE void *
get_cb_up_options_default(void)
{

  return NULL;
}

PRIVATE vrna_auxdata_free_f
get_cb_up_options_free_default(void)
{
  return NULL;
}


PRIVATE void
nullify_probing_data_s(struct vrna_probing_data_s *data)
{
  data->method = 0;

  data->data_linear              = NULL;
  data->data_linear_weight       = NULL;
  data->cbs_linear               = NULL;
  data->cbs_linear_options       = NULL;
  data->cbs_linear_options_free  = NULL;

  data->params1 = NULL;
  data->params2 = NULL;
  data->reactivities = NULL;
  data->transformeds = NULL;
  data->datas1 = NULL;
  data->datas2 = NULL;
}
