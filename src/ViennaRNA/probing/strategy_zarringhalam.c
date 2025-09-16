/* SHAPE reactivity data handling */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "ViennaRNA/utils/basic.h"

#include "ViennaRNA/probing/strategies.h"

typedef struct {
  double beta;
  vrna_probing_transform_f  cb_preprocess;
  void                      *cb_preprocess_opt;
  vrna_auxdata_free_f       cb_preprocess_opt_free;
} zarringhalam_options_t;

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

PRIVATE FLT_OR_DBL
conversion_zarringhalam_up(double       pr,
                           double       beta);


PRIVATE FLT_OR_DBL
conversion_zarringhalam_bp(double       pr_i,
                           double       pr_j,
                           double       beta);


PRIVATE vrna_probing_transform_f
set_mapping_strategy(const char          *conversion_string,
                     double              default_value,
                     void                **transform_data,
                     vrna_auxdata_free_f *transform_data_free);


/*
 #################################
 # BEGIN OF FUNCTION DEFINITIONS #
 #################################
 */


PUBLIC double *
vrna_probing_strategy_zarringhalam_up(const double *data,
                                      size_t       data_size,
                                      void         *options)
{
  double                  *pseudo_energies;
  zarringhalam_options_t  *opt;

  pseudo_energies = NULL;

  if ((data) &&
      (data_size > 0)) {

    /* preprare (default) options */
    if (options)
      opt = (zarringhalam_options_t *)options;
    else
      opt = vrna_probing_strategy_zarringhalam_options(VRNA_PROBING_METHOD_ZARRINGHALAM2012_DEFAULT_beta,
                                                       NULL,
                                                       NULL,
                                                       NULL);

    /* pre-process data */
    pseudo_energies = vrna_reactivity_transform(data_size,
                                                data,
                                                opt->cb_preprocess,
                                                opt->cb_preprocess_opt);

    /* transform data into actual pseudo-energies */
    for (size_t i = 0; i <= data_size; i++)
      pseudo_energies[i] = conversion_zarringhalam_up(pseudo_energies[i],
                                                      opt->beta);

    /* release memory for default options */
    if (opt != (zarringhalam_options_t *)options)
      free(opt);
  }

  return pseudo_energies;
}


PUBLIC void *
vrna_probing_strategy_zarringhalam_options(double                   beta,
                                           vrna_probing_transform_f cb_preprocess,
                                           void                     *cb_preprocess_opt,
                                           vrna_auxdata_free_f      cb_preprocess_opt_free)
{
  zarringhalam_options_t  *opt = (zarringhalam_options_t *)vrna_alloc(sizeof(zarringhalam_options_t));

  opt->beta = beta;

  if (cb_preprocess) {
    opt->cb_preprocess = (cb_preprocess) ? cb_preprocess : vrna_reactivity_trans_method(VRNA_REACTIVITY_TRANS_NEG_IGNORE);
    opt->cb_preprocess_opt      = cb_preprocess_opt;
    opt->cb_preprocess_opt_free = cb_preprocess_opt_free;
  } else {
    /* default preprocessing of the probing data */
    opt->cb_preprocess = set_mapping_strategy(VRNA_PROBING_METHOD_ZARRINGHALAM2012_DEFAULT_conversion,
                                              VRNA_PROBING_METHOD_ZARRINGHALAM2012_DEFAULT_probability,
                                              &(opt->cb_preprocess_opt),
                                              &(opt->cb_preprocess_opt_free));
  }

  return (void *)opt;
}


PUBLIC void
vrna_probing_strategy_zarringhalam_options_free(void *options)
{
  zarringhalam_options_t  *opt = (zarringhalam_options_t *)options;

  if (opt->cb_preprocess_opt_free)
    opt->cb_preprocess_opt_free(opt->cb_preprocess_opt);

  free(opt);
}


/*
 #####################################
 # BEGIN OF STATIC HELPER FUNCTIONS  #
 #####################################
 */

PRIVATE FLT_OR_DBL
conversion_zarringhalam_up(double       pr,
                           double       beta)
{
  return beta * fabs(pr - 1.);
}


PRIVATE FLT_OR_DBL
conversion_zarringhalam_bp(double       pr_i,
                           double       pr_j,
                           double       beta)
{
  return beta * (pr_i + pr_j);
}


PRIVATE vrna_probing_transform_f
set_mapping_strategy(const char          *conversion_string,
                     double              default_value,
                     void                **transform_data,
                     vrna_auxdata_free_f *transform_data_free)
{
  int *indices;
  int i, j;
  int index;
  int ret = 1;
#if 1
  vrna_probing_transform_f  transform_function;

  transform_function    = NULL;
  *transform_data       = NULL;
  *transform_data_free  = NULL;

  switch (*conversion_string) {
    case 'S':
      transform_function = vrna_reactivity_trans_method(VRNA_REACTIVITY_TRANS_IDEN);
      break;

    case 'M':
      break;

    case 'C':
      break;

    case 'L':
      /* fall through */
    case 'O':
      transform_function = NULL;
      break;

    default:
      break;
  }

#else
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
#endif

  return NULL;
}


typedef struct {
  vrna_probing_transform_f cutoff_cb;
  void                     *cutoff_cb_data;
} zh_trans_options_t;


PRIVATE double
zh_trans_identity(double  value,
                  void    *options)
{
  zh_trans_options_t *o = (zh_trans_options_t *)options;

  return o->cutoff_cb(value, o->cutoff_cb_data);
}
