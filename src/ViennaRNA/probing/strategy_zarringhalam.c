/* SHAPE reactivity data handling */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "ViennaRNA/utils/basic.h"
#include "ViennaRNA/data/transform.h"

#include "ViennaRNA/probing/strategy_zarringhalam.h"

typedef struct {
  double                        beta;
  double                        default_probability;
  vrna_math_fun_dbl_f           cb_preprocess;
  vrna_math_fun_dbl_opt_t       cb_preprocess_opt;
  vrna_math_fun_dbl_opt_free_f  cb_preprocess_opt_free;
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
conversion_zarringhalam_up(double pr,
                           double beta);


PRIVATE FLT_OR_DBL
conversion_zarringhalam_bp(double pr,
                           double beta);


PRIVATE vrna_math_fun_dbl_f
set_mapping_strategy(const char                   *conversion_string,
                     double                       max_value,
                     vrna_math_fun_dbl_opt_t      *transform_data,
                     vrna_math_fun_dbl_opt_free_f *transform_data_free);


/*
 #################################
 # BEGIN OF FUNCTION DEFINITIONS #
 #################################
 */
PUBLIC double *
vrna_probing_strategy_zarringhalam(vrna_fold_compound_t *fc,
                                   const double         *data,
                                   size_t               data_size,
                                   unsigned int         target,
                                   void                 *options)
{
  double                  *pseudo_energies;
  zarringhalam_options_t  *opt;

  pseudo_energies = NULL;

  if (((target == VRNA_PROBING_DATA_LINEAR_TARGET_UP) ||
       (target == VRNA_PROBING_DATA_LINEAR_TARGET_BP)) &&
      (data) &&
      (data_size > 0)) {
    double max_v = data[0];
    for (size_t i = 1; i < data_size; ++i)
      max_v = MAX2(max_v, data[i]);

    /* preprare (default) options */
    if (options) {
      opt = (zarringhalam_options_t *)options;
    } else {
      opt =
        vrna_probing_strategy_zarringhalam_options(
          VRNA_PROBING_METHOD_ZARRINGHALAM2012_DEFAULT_beta,
          VRNA_PROBING_METHOD_ZARRINGHALAM2012_DEFAULT_probability,
          max_v,
          NULL,
          NULL,
          NULL);
    }

    /* pre-process data */
    double domain[4] = {
      0., max_v, 0., 1.
    };

    pseudo_energies = vrna_data_lin_transform(data,
                                              data_size,
                                              opt->cb_preprocess,
                                              opt->cb_preprocess_opt,
                                              domain,
                                              VRNA_REACTIVITY_MISSING,
                                              VRNA_TRANSFORM_ENFORCE_DOMAINS |
                                              VRNA_TRANSFORM_MAP_TARGET);

    /* transform data into actual pseudo-energies */
    for (size_t i = 0; i < data_size; i++)
      if (pseudo_energies[i] == VRNA_REACTIVITY_MISSING)
        pseudo_energies[i] = opt->default_probability;

    /* transform data into actual pseudo-energies */
    if (target == VRNA_PROBING_DATA_LINEAR_TARGET_UP)
      for (size_t i = 0; i < data_size; i++)
        pseudo_energies[i] = conversion_zarringhalam_up(pseudo_energies[i], opt->beta);
    else
      for (size_t i = 0; i < data_size; i++)
        pseudo_energies[i] = conversion_zarringhalam_bp(pseudo_energies[i], opt->beta);

    /* release memory for default options */
    if (opt != (zarringhalam_options_t *)options)
      vrna_probing_strategy_zarringhalam_options_free(options);
  }

  return pseudo_energies;
}


PUBLIC void *
vrna_probing_strategy_zarringhalam_options(double                       beta,
                                           double                       default_probability,
                                           double                       max_value,
                                           vrna_math_fun_dbl_f          cb_preprocess,
                                           vrna_math_fun_dbl_opt_t      cb_preprocess_opt,
                                           vrna_math_fun_dbl_opt_free_f cb_preprocess_opt_free)
{
  zarringhalam_options_t *opt =
    (zarringhalam_options_t *)vrna_alloc(sizeof(zarringhalam_options_t));

  opt->beta                 = beta;
  opt->default_probability  = ((default_probability < 0) || (default_probability > 1.)) ?
                              VRNA_PROBING_METHOD_ZARRINGHALAM2012_DEFAULT_probability :
                              default_probability;

  if (cb_preprocess) {
    opt->cb_preprocess          = cb_preprocess;
    opt->cb_preprocess_opt      = cb_preprocess_opt;
    opt->cb_preprocess_opt_free = cb_preprocess_opt_free;
  } else {
    /* default preprocessing of the probing data */
    opt->cb_preprocess =
      set_mapping_strategy(VRNA_PROBING_METHOD_ZARRINGHALAM2012_DEFAULT_conversion,
                           max_value,
                           &(opt->cb_preprocess_opt),
                           &(opt->cb_preprocess_opt_free));
  }

  return (void *)opt;
}


PUBLIC void
vrna_probing_strategy_zarringhalam_options_free(void *options)
{
  zarringhalam_options_t *opt = (zarringhalam_options_t *)options;

  if (opt->cb_preprocess_opt_free)
    opt->cb_preprocess_opt_free(opt->cb_preprocess_opt);

  free(opt);
}


PUBLIC struct vrna_probing_data_s *
vrna_probing_data_zarringhalam(const double *reactivities,
                               unsigned int n,
                               double       beta,
                               const char   *pr_conversion,
                               double       pr_default)
{
  vrna_math_fun_dbl_f           trans;
  vrna_math_fun_dbl_opt_t       trans_options;
  vrna_math_fun_dbl_opt_free_f  trans_options_free;

  if (reactivities) {
    double max = reactivities[0];
    for (size_t i = 1; i <= n; i++)
      max = MAX2(max, reactivities[i]);

    trans = set_mapping_strategy(((pr_conversion) && (*pr_conversion)) ?
                                 pr_conversion :
                                 VRNA_PROBING_METHOD_ZARRINGHALAM2012_DEFAULT_conversion,
                                 max,
                                 &(trans_options),
                                 &(trans_options_free));

    return vrna_probing_data_linear(reactivities,
                                    n,
                                    NULL,
                                    vrna_probing_strategy_zarringhalam,
                                    vrna_probing_strategy_zarringhalam_options(beta,
                                                                               pr_default,
                                                                               max,
                                                                               trans,
                                                                               trans_options,
                                                                               trans_options_free),
                                    vrna_probing_strategy_zarringhalam_options_free,
                                    VRNA_PROBING_DATA_DEFAULT);
  }

  return NULL;
}


PUBLIC struct vrna_probing_data_s *
vrna_probing_data_zarringhalam_trans(const double                 *reactivities,
                                     unsigned int                 n,
                                     double                       beta,
                                     double                       pr_default,
                                     vrna_math_fun_dbl_f          trans,
                                     vrna_math_fun_dbl_opt_t      trans_options,
                                     vrna_math_fun_dbl_opt_free_f trans_options_free)
{
  if (reactivities) {
    double max = reactivities[0];
    for (size_t i = 1; i <= n; i++)
      max = MAX2(max, reactivities[i]);

    return vrna_probing_data_linear(reactivities,
                                    n,
                                    NULL,
                                    vrna_probing_strategy_zarringhalam,
                                    vrna_probing_strategy_zarringhalam_options(beta,
                                                                               pr_default,
                                                                               max,
                                                                               trans,
                                                                               trans_options,
                                                                               trans_options_free),
                                    vrna_probing_strategy_zarringhalam_options_free,
                                    VRNA_PROBING_DATA_DEFAULT);
  }

  return NULL;
}


PUBLIC struct vrna_probing_data_s *
vrna_probing_data_zarringhalam_comparative(const double **reactivities,
                                           unsigned int *n,
                                           unsigned int n_seq,
                                           double       *betas,
                                           const char   **pr_conversions,
                                           double       *pr_defaults,
                                           unsigned int multi_params)
{
  vrna_array(vrna_math_fun_dbl_f)           trans;
  vrna_array(vrna_math_fun_dbl_opt_t)       trans_options;
  vrna_array(vrna_math_fun_dbl_opt_free_f)  trans_options_free;
  struct vrna_probing_data_s *d;

  if ((reactivities) &&
      (n) &&
      (n_seq > 0)) {
    vrna_array_init_size(trans, n_seq);
    vrna_array_init_size(trans_options, n_seq);
    vrna_array_init_size(trans_options_free, n_seq);

    for (size_t s = 0; s < n_seq; s++) {
      if (reactivities[s]) {
        vrna_math_fun_dbl_f           cb_trans;
        vrna_math_fun_dbl_opt_t       cb_trans_options;
        vrna_math_fun_dbl_opt_free_f  cb_trans_options_free;
        double                        max = reactivities[s][0];
        for (size_t i = 1; i <= n[s]; i++)
          max = MAX2(max, reactivities[s][i]);

        cb_trans = set_mapping_strategy(((pr_conversions) && (pr_conversions[s]) &&
                                         (*pr_conversions[s])) ?
                                        pr_conversions[s] :
                                        VRNA_PROBING_METHOD_ZARRINGHALAM2012_DEFAULT_conversion,
                                        max,
                                        &(cb_trans_options),
                                        &(cb_trans_options_free));

        vrna_array_append(trans, cb_trans);
        vrna_array_append(trans_options, cb_trans_options);
        vrna_array_append(trans_options_free, cb_trans_options_free);
      } else {
        vrna_array_append(trans, NULL);
        vrna_array_append(trans_options, NULL);
        vrna_array_append(trans_options_free, NULL);
      }
    }

    d = vrna_probing_data_zarringhalam_trans_comparative(reactivities,
                                                         n,
                                                         n_seq,
                                                         betas,
                                                         pr_defaults,
                                                         multi_params,
                                                         trans,
                                                         trans_options,
                                                         trans_options_free);

    vrna_array_free(trans);
    vrna_array_free(trans_options);
    vrna_array_free(trans_options_free);

    return d;
  }

  return NULL;
}


PUBLIC struct vrna_probing_data_s *
vrna_probing_data_zarringhalam_trans_comparative(const double                 **reactivities,
                                                 unsigned int                 *n,
                                                 unsigned int                 n_seq,
                                                 double                       *betas,
                                                 double                       *pr_defaults,
                                                 unsigned int                 multi_params,
                                                 vrna_math_fun_dbl_f          *trans,
                                                 vrna_math_fun_dbl_opt_t      *trans_options,
                                                 vrna_math_fun_dbl_opt_free_f *trans_options_free)
{
  struct vrna_probing_data_s    *d = NULL;
  double                        beta;
  double                        pr_default;

  vrna_math_fun_dbl_f           cb_trans;
  vrna_math_fun_dbl_opt_t       cb_trans_options;
  vrna_math_fun_dbl_opt_free_f  cb_trans_options_free;

  vrna_array(vrna_probing_strategy_f)   cbs_linear;
  vrna_array(void *)                    cbs_linear_options;
  vrna_array(vrna_auxdata_free_f)       cbs_linear_options_free;

  if ((reactivities) &&
      (n)) {
    beta        = (betas) ? *betas : VRNA_PROBING_METHOD_ZARRINGHALAM2012_DEFAULT_beta;
    pr_default  =
      (pr_defaults) ? *pr_defaults : VRNA_PROBING_METHOD_ZARRINGHALAM2012_DEFAULT_probability;
    cb_trans              = (trans) ? *trans : NULL;
    cb_trans_options      = (trans_options) ? *trans_options : NULL;
    cb_trans_options_free = (trans_options_free) ? *trans_options_free : NULL;

    if (((betas == NULL) && (multi_params & VRNA_PROBING_METHOD_MULTI_PARAMS_1)) ||
        ((pr_defaults == NULL) && (multi_params & VRNA_PROBING_METHOD_MULTI_PARAMS_2)) ||
        ((trans == NULL) && (multi_params & VRNA_PROBING_METHOD_MULTI_PARAMS_3))) {
      return d;
    }

    /* prepare callback arrays */
    vrna_array_init_size(cbs_linear, n_seq);
    vrna_array_init_size(cbs_linear_options, n_seq);
    vrna_array_init_size(cbs_linear_options_free, n_seq);

    /* fill callback vectors */
    for (size_t s = 0; s < n_seq; s++) {
      if (reactivities[s]) {
        double max = reactivities[s][0];
        for (size_t i = 1; i <= n[s]; i++)
          max = MAX2(max, reactivities[s][i]);

        if (multi_params & VRNA_PROBING_METHOD_MULTI_PARAMS_1)
          beta = betas[s];

        if (multi_params & VRNA_PROBING_METHOD_MULTI_PARAMS_2)
          pr_default = pr_defaults[s];

        if (multi_params & VRNA_PROBING_METHOD_MULTI_PARAMS_3) {
          cb_trans = trans[s];
          if (trans_options)
            cb_trans_options = trans_options[s];

          if (trans_options_free)
            cb_trans_options_free = trans_options_free[s];
        }

        vrna_array_append(cbs_linear, vrna_probing_strategy_zarringhalam);
        vrna_array_append(cbs_linear_options, vrna_probing_strategy_zarringhalam_options(beta,
                                                                                         pr_default,
                                                                                         max,
                                                                                         cb_trans,
                                                                                         cb_trans_options,
                                                                                         cb_trans_options_free));
        vrna_array_append(cbs_linear_options_free, vrna_probing_strategy_zarringhalam_options_free);
      } else {
        vrna_array_append(cbs_linear, NULL);
        vrna_array_append(cbs_linear_options, NULL);
        vrna_array_append(cbs_linear_options_free, NULL);
      }
    }

    d = vrna_probing_data_linear_multi(reactivities,
                                       n_seq,
                                       n,
                                       NULL,
                                       cbs_linear,
                                       cbs_linear_options,
                                       cbs_linear_options_free,
                                       VRNA_PROBING_DATA_DEFAULT);

    vrna_array_free(cbs_linear);
    vrna_array_free(cbs_linear_options);
    vrna_array_free(cbs_linear_options_free);
  }

  return d;
}


/*
 #####################################
 # BEGIN OF STATIC HELPER FUNCTIONS  #
 #####################################
 */
PRIVATE FLT_OR_DBL
conversion_zarringhalam_up(double pr,
                           double beta)
{
  return beta * fabs(pr - 1.);
}


PRIVATE FLT_OR_DBL
conversion_zarringhalam_bp(double pr,
                           double beta)
{
  return beta * pr;
}


PRIVATE vrna_math_fun_dbl_f
set_mapping_strategy(const char                   *conversion_string,
                     double                       max_value,
                     vrna_math_fun_dbl_opt_t      *transform_data,
                     vrna_math_fun_dbl_opt_free_f *transform_data_free)
{
  vrna_math_fun_dbl_f transform_function;

  transform_function    = NULL;
  *transform_data       = NULL;
  *transform_data_free  = NULL;

  switch (*conversion_string) {
    case 'S':
      /* simple identity transformation where negative values are cut-off */
    {
      double map[2][2] = { {
                             0., 0.
                           },
                           {
                             0., 1.0
                           } };
      map[1][0] = map[1][1] = max_value;

      transform_function = vrna_math_fun_dbl_bin_opt(&(map[0]),
                                                     2,
                                                     VRNA_REACTIVITY_MISSING,
                                                     VRNA_REACTIVITY_MISSING,
                                                     VRNA_MATH_FUN_BIN_OPTION_PROJECT,
                                                     transform_data,
                                                     transform_data_free);
    }
    break;

    case 'M':
      /* default mapping from Zarringhalam et al. 2012 */
    {
      double map[5][2] = { {
                             0., 0.
                           },
                           {
                             0.25, 0.35
                           },
                           {
                             0.3, 0.55
                           },
                           {
                             0.7, 0.85
                           },
                           {
                             0., 1.0
                           } };
      map[4][0]           = max_value;
      transform_function  = vrna_math_fun_dbl_bin_opt(&(map[0]),
                                                      5,
                                                      VRNA_REACTIVITY_MISSING,
                                                      VRNA_REACTIVITY_MISSING,
                                                      VRNA_MATH_FUN_BIN_OPTION_PROJECT,
                                                      transform_data,
                                                      transform_data_free);
    }
    break;

    case 'C':
      /* cut-off mapping into two bins with values 0 and 1, also negative values are treated as missing */
    {
      float   cutoff = 0.25;
      sscanf(conversion_string + 1, "%f", &cutoff);
      double  map[2][2] = { {
                              0., 0.
                            },
                            {
                              0., 1.0
                            } };
      map[1][0]           = cutoff;
      transform_function  = vrna_math_fun_dbl_bin_opt(&(map[0]),
                                                      2,
                                                      VRNA_REACTIVITY_MISSING,
                                                      VRNA_REACTIVITY_MISSING,
                                                      VRNA_MATH_FUN_BIN_OPTION_MAP_OUTOF_UPPERBOUND,
                                                      transform_data,
                                                      transform_data_free);
    }
    break;

    case 'L':
    /* fall through */
    case 'O':
    {
      float slope     = (*conversion_string == 'L') ? 0.68 : 1.6;
      float intercept = (*conversion_string == 'L') ? 0.2 : -2.29;

      /* try parsing other parameter settings */
      int   r;

      if (conversion_string[1]) {
        r = sscanf(conversion_string + 1, "s%fi%f", &slope, &intercept);

        if (r < 2) {
          r = sscanf(conversion_string + 1, "s%f", &slope);

          if (r < 1) {
            r = sscanf(conversion_string + 1,
                       "i%f",
                       &intercept);
            if (r < 1)
              vrna_log_warning(
                "Probing method parameters not recognized! Using default parameters!");
          }
        }
      }

      /*
       * weird hack to maintain backward compatibility where the formula y = (f(x) - i) / s was
       * misinterpreted and s being used as slope and i asintercept.
       */
      slope     = 1. / slope;
      intercept *= -slope;

      transform_function = vrna_math_fun_dbl_linear_opt(slope,
                                                        intercept,
                                                        (*conversion_string == 'L') ?
                                                        VRNA_MATH_FUN_LINEAR_OPTION_DEFAULT :
                                                        VRNA_MATH_FUN_LINEAR_OPTION_LOG,
                                                        transform_data,
                                                        transform_data_free);
    }
    break;

    default:
      break;
  }

  return transform_function;
}
