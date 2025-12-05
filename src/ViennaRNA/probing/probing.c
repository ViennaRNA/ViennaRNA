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
#include "ViennaRNA/utils/basic.h"
#include "ViennaRNA/utils/strings.h"
#include "ViennaRNA/sequences/alignments.h"
#include "ViennaRNA/utils/log.h"
#include "ViennaRNA/io/utils.h"
#include "ViennaRNA/io/file_formats.h"
#include "ViennaRNA/params/basic.h"
#include "ViennaRNA/constraints/soft.h"
#include "ViennaRNA/probing/basic.h"
#include "ViennaRNA/probing/strategy_deigan.h"

#ifdef __GNUC__
# define INLINE inline
#else
# define INLINE
#endif


struct vrna_probing_data_s {
  vrna_array(double *)                  data_linear;          /* actual data */
  vrna_array(double *)                  data_linear_weight;   /* weight for each data set */
  vrna_array(vrna_probing_strategy_f)   cbs_linear;
  vrna_array(void *)                    cbs_linear_options;
  vrna_array(vrna_auxdata_free_f)       cbs_linear_options_free;
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


INLINE PRIVATE double *
prepare_linear_data(vrna_fold_compound_t        *fc,
                    struct vrna_probing_data_s  *data,
                    size_t                      index,
                    unsigned int                target);


PRIVATE int
apply_probing_data(vrna_fold_compound_t       *fc,
                   struct vrna_probing_data_s *data);


PRIVATE double
get_msa_weight(const double **datasets,
               unsigned int data_size);


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
vrna_probing_data_linear(const double             *data,
                         unsigned int             data_length,
                         const double             *data_weights,
                         vrna_probing_strategy_f  strategy_cb,
                         void                     *strategy_cb_options,
                         vrna_auxdata_free_f      strategy_cb_options_free,
                         unsigned int             options)
{
  struct vrna_probing_data_s *d = NULL;

  if (data) {
    return vrna_probing_data_linear_multi(&data,
                                          1,
                                          &data_length,
                                          (data_weights) ? &data_weights : NULL,
                                          (strategy_cb) ? &strategy_cb : NULL,
                                          (strategy_cb_options) ? &strategy_cb_options : NULL,
                                          (strategy_cb_options_free) ? &strategy_cb_options_free :
                                          NULL,
                                          options);
  }

  return d;
}


PUBLIC vrna_probing_data_t
vrna_probing_data_linear_multi(const double             **data,
                               unsigned int             data_size,
                               const unsigned int       *data_lengths,
                               const double             **data_weights,
                               vrna_probing_strategy_f  *strategy_cbs,
                               void                     **strategy_cbs_options,
                               vrna_auxdata_free_f      *strategy_cbs_options_free,
                               unsigned int             options)
{
  double                      weight = 1.0;

  struct vrna_probing_data_s  *d = NULL;
  vrna_probing_strategy_f     cb;
  void                        *cb_options;
  vrna_auxdata_free_f         cb_options_free;

  if ((data) &&
      (data_lengths) &&
      (data_size)) {
    d = (struct vrna_probing_data_s *)vrna_alloc(sizeof(struct vrna_probing_data_s));

    nullify_probing_data_s(d);

    cb              = NULL;
    cb_options      = NULL;
    cb_options_free = NULL;

    vrna_array_init_size(d->data_linear, data_size);
    vrna_array_init_size(d->data_linear_weight, data_size);
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

    if (options & VRNA_PROBING_DATA_WEIGHT_POSITION_WISE) {
      if (!data_weights) {
        vrna_log_error("No weight data given despite request for position wise weights");
        vrna_log_warning("deactivating position wise weighting");

        options &= ~VRNA_PROBING_DATA_WEIGHT_POSITION_WISE;
      } else {
        /*  check for consistency, i.e. a weight vector must be present for
         *  each data set
         */
        for (size_t i = 0; i < data_size; i++) {
          if ((data[i]) &&
              (data_lengths[i]) &&
              (data_weights[i] == NULL)) {
            vrna_log_error(
              "Missing weight data for data set %ld despite request for position wise weights",
              i);
            vrna_log_warning("deactivating position wise weighting");
            options &= ~VRNA_PROBING_DATA_WEIGHT_POSITION_WISE;
            break;
          }
        }
      }
    }

    if ((options & VRNA_PROBING_DATA_SINGLE_WEIGHT) &&
        (!(options & VRNA_PROBING_DATA_WEIGHT_POSITION_WISE))) {
      if (data_weights) {
        weight = data_weights[0][0];
      } else {
        weight = get_msa_weight(data, data_size);
      }
    } else if (!(data_weights)) {
      weight = get_msa_weight(data, data_size);
    }

    for (size_t i = 0; i < data_size; i++) {
      if ((data[i]) &&
          (data_lengths[i])) {
        vrna_array(double)  a;

        /* init and store raw probing data */
        vrna_array_init_size(a, data_lengths[i] + 1);
        for (size_t j = 0; j <= data_lengths[i]; j++)
          vrna_array_append(a, data[i][j]);

        vrna_array_append(d->data_linear, a);
        a = NULL;

        /* prepare and store weights for this data set */
        vrna_array_init_size(a, data_lengths[i] + 1);

        if (options & VRNA_PROBING_DATA_WEIGHT_POSITION_WISE) {
          if (options & VRNA_PROBING_DATA_SINGLE_WEIGHT) {
            /* use weights from first weight vector for all data sets */
            for (size_t j = 0; j <= data_lengths[i]; j++)
              vrna_array_append(a, data_weights[0][j]);
          } else {
            for (size_t j = 0; j <= data_lengths[i]; j++)
              vrna_array_append(a, data_weights[i][j]);
          }
        } else {
          if (options & VRNA_PROBING_DATA_SINGLE_WEIGHT) {
            /* use the same weight for all positions in all data sets */
            for (size_t j = 0; j <= data_lengths[i]; j++)
              vrna_array_append(a, weight);
          } else {
            for (size_t j = 0; j <= data_lengths[i]; j++)
              vrna_array_append(a,
                                ((data_weights) &&
                                 (data_weights[i])) ? data_weights[i][0] : weight);
          }
        }

        vrna_array_append(d->data_linear_weight, a);
        a = NULL;

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
        vrna_array_append(d->data_linear_weight, NULL);
        vrna_array_append(d->cbs_linear, NULL);
        vrna_array_append(d->cbs_linear_options, NULL);
        vrna_array_append(d->cbs_linear_options_free, NULL);
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
    if (d->data_linear)
      for (size_t i = 0; i < vrna_array_size(d->data_linear); i++)
        vrna_array_free(d->data_linear[i]);

    vrna_array_free(d->data_linear);

    /* free weights */
    if (d->data_linear_weight)
      for (size_t i = 0; i < vrna_array_size(d->data_linear_weight); i++)
        vrna_array_free(d->data_linear_weight[i]);

    vrna_array_free(d->data_linear_weight);

    /* free probing strategy callback vector */
    vrna_array_free(d->cbs_linear);

    /* free probing strategy options vector */
    if ((d->cbs_linear_options_free) &&
        (d->cbs_linear_options)) {
      for (size_t i = 0; i < vrna_array_size(d->cbs_linear_options_free); i++)
        if (d->cbs_linear_options_free[i])
          d->cbs_linear_options_free[i](d->cbs_linear_options[i]);
    }

    vrna_array_free(d->cbs_linear_options);
    vrna_array_free(d->cbs_linear_options_free);

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
    if (values[i] < 0)
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
      index         = indices[i];
      v             = (*shape_conversion == 'L') ? values[index] : log(values[index]);
      values[index] = MAX2(MIN2((v - intercept) / slope, 1), 0);
    }
  } else {
    ret = 0;
  }

  free(indices);

  return ret;
}


PUBLIC unsigned int
vrna_probing_data_linear_num(struct vrna_probing_data_s *data)
{
  return (data) ? vrna_array_size(data->data_linear) : 0;
}


PUBLIC double *
vrna_probing_data_linear_raw(struct vrna_probing_data_s *data,
                             unsigned int               pos,
                             unsigned int               *data_size)
{
  double *raw;

  raw         = NULL;
  *data_size  = 0;

  if ((data) &&
      (vrna_array_size(data->data_linear) >= pos) &&
      (data->data_linear[pos])) {
    *data_size  = vrna_array_size(data->data_linear[pos]);
    raw         = vrna_alloc(sizeof(double) * *data_size);
    raw         = memcpy(raw, data->data_linear[pos], sizeof(double) * *data_size);
  }

  return raw;
}


PUBLIC double *
vrna_probing_data_linear_weight(struct vrna_probing_data_s  *data,
                                unsigned int                pos,
                                unsigned int                *data_size)
{
  double *w;

  w           = NULL;
  *data_size  = 0;

  if ((data) &&
      (vrna_array_size(data->data_linear_weight) >= pos) &&
      (data->data_linear_weight[pos])) {
    *data_size  = vrna_array_size(data->data_linear_weight[pos]);
    w           = vrna_alloc(sizeof(double) * *data_size);
    w           = memcpy(w, data->data_linear_weight[pos], sizeof(double) * *data_size);
  }

  return w;
}


PUBLIC double *
vrna_probing_data_linear_energies(struct vrna_probing_data_s  *data,
                                  unsigned int                pos,
                                  vrna_fold_compound_t        *fc,
                                  unsigned int                target,
                                  unsigned int                *data_size)
{
  double *e;

  e           = NULL;
  *data_size  = 0;

  if ((data) &&
      (vrna_array_size(data->data_linear) > pos) &&
      (data->data_linear[pos])) {
    *data_size  = vrna_array_size(data->data_linear[pos]);
    e           = prepare_linear_data(fc, data, pos, target);
  }

  return e;
}


/*
 #####################################
 # BEGIN OF STATIC HELPER FUNCTIONS  #
 #####################################
 */
INLINE PRIVATE double *
prepare_linear_data(vrna_fold_compound_t        *fc,
                    struct vrna_probing_data_s  *data,
                    size_t                      index,
                    unsigned int                target)
{
  return data->cbs_linear[index](fc,
                                 data->data_linear[index],
                                 vrna_array_size(data->data_linear[index]),
                                 target,
                                 data->cbs_linear_options[index]);
}


PRIVATE int
apply_probing_data(vrna_fold_compound_t       *fc,
                   struct vrna_probing_data_s *data)
{
  unsigned int  i, j, s, n;
  int           ret, num_data;
  double        *e;

  ret       = 0;
  num_data  = 0;
  n         = fc->length;

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
          e = prepare_linear_data(fc, data, 0, VRNA_PROBING_DATA_LINEAR_TARGET_STACK);

          if (e) {
            n = MIN2(n, vrna_array_size(data->data_linear[0]) - 1);

            for (i = 1; i <= n; ++i)
              ret &= vrna_sc_add_stack(fc,
                                       i,
                                       e[i] * data->data_linear_weight[0][i],
                                       VRNA_OPTION_DEFAULT);

            free(e);
          }

          e = prepare_linear_data(fc, data, 0, VRNA_PROBING_DATA_LINEAR_TARGET_UP);

          if (e) {
            n = MIN2(n, vrna_array_size(data->data_linear[0]) - 1);

            for (i = 1; i <= n; ++i) {
              ret &= vrna_sc_add_up(fc,
                                    i,
                                    e[i] * data->data_linear_weight[0][i],
                                    VRNA_OPTION_DEFAULT);
            }

            free(e);
          }

          e = prepare_linear_data(fc, data, 0, VRNA_PROBING_DATA_LINEAR_TARGET_BP);

          if (e) {
            n = MIN2(n, vrna_array_size(data->data_linear[0]) - 1);

            for (i = 1; i <= n; ++i)
              for (j = i + 1; j <= n; ++j) {
                ret &= vrna_sc_add_bp(fc,
                                      i,
                                      j,
                                      (e[i] * data->data_linear_weight[0][i]) +
                                      (e[j] * data->data_linear_weight[0][j]),
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
          vrna_log_warning(
            "Number of probing data (%u) doesn't match number of sequences in the alignment (%u)",
            vrna_array_size(data->data_linear),
            fc->n_seq);

        for (s = 0; s < vrna_array_size(data->data_linear); s++) {
          if (vrna_array_size(data->data_linear[s]) > 0) {
            n = MIN2(vrna_array_size(data->data_linear[s]), fc->alignment->gapfree_size[s]);

            num_data++;

            if (vrna_array_size(data->data_linear[s] - 1) != n)
              vrna_log_warning(
                "Length of probing data (%u) for sequence no. \"%u\" doesn't match length of gap-free sequence (%u)",
                vrna_array_size(data->data_linear[s]) - 1,
                s,
                n);

            e = prepare_linear_data(fc, data, s, VRNA_PROBING_DATA_LINEAR_TARGET_STACK);

            if (e) {
              for (i = 1; i <= n; ++i)
                ret &= vrna_sc_add_stack_comparative_seq(fc,
                                                         s,
                                                         i,
                                                         e[i] * data->data_linear_weight[s][i],
                                                         VRNA_OPTION_DEFAULT);

              free(e);
            }

            e = prepare_linear_data(fc, data, s, VRNA_PROBING_DATA_LINEAR_TARGET_UP);

            if (e) {
              for (i = 1; i <= n; ++i)
                ret &= vrna_sc_add_up_comparative_seq(fc,
                                                      s,
                                                      i,
                                                      e[i] * data->data_linear_weight[s][i],
                                                      VRNA_OPTION_DEFAULT);

              free(e);
            }

            e = prepare_linear_data(fc, data, s, VRNA_PROBING_DATA_LINEAR_TARGET_BP);

            if (e) {
              for (i = 1; i <= n; ++i)
                for (j = i + 1; j <= n; ++j) {
                  ret &= vrna_sc_add_bp_comparative_seq(fc,
                                                        s,
                                                        i,
                                                        j,
                                                        (e[i] * data->data_linear_weight[s][i]) +
                                                        (e[j] * data->data_linear_weight[s][j]),
                                                        VRNA_OPTION_DEFAULT);
                }

              free(e);
            }
          }
        }

        if (ret)
          ret = num_data;

        break;
    }
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
  data->data_linear             = NULL;
  data->data_linear_weight      = NULL;
  data->cbs_linear              = NULL;
  data->cbs_linear_options      = NULL;
  data->cbs_linear_options_free = NULL;
}
