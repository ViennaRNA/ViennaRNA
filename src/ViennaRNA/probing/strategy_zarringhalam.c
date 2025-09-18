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
  double                    beta;
  double                    default_probability;
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
                     double              max_value,
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
    if (options) {
      opt = (zarringhalam_options_t *)options;
    } else {
      double max_v = data[0];
      for (size_t i = 1; i < data_size; ++i)
        max_v = MAX2(max_v, data[i]);

      opt = vrna_probing_strategy_zarringhalam_options(VRNA_PROBING_METHOD_ZARRINGHALAM2012_DEFAULT_beta,
                                                       VRNA_PROBING_METHOD_ZARRINGHALAM2012_DEFAULT_probability,
                                                       max_v,
                                                       NULL,
                                                       NULL,
                                                       NULL);
    }

    /* pre-process data */
    pseudo_energies = vrna_reactivity_transform(data_size,
                                                data,
                                                opt->cb_preprocess,
                                                opt->cb_preprocess_opt);

    /* transform data into actual pseudo-energies */
    for (size_t i = 0; i <= data_size; i++)
      if (pseudo_energies[i] == VRNA_REACTIVITY_MISSING)
        pseudo_energies[i] = opt->default_probability;
      else
        pseudo_energies[i] = conversion_zarringhalam_up(pseudo_energies[i],
                                                        opt->beta);

    /* release memory for default options */
    if (opt != (zarringhalam_options_t *)options)
      vrna_probing_strategy_zarringhalam_options_free(options);
  }

  return pseudo_energies;
}


PUBLIC void *
vrna_probing_strategy_zarringhalam_options(double                   beta,
                                           double                   default_probability,
                                           double                   max_value,
                                           vrna_probing_transform_f cb_preprocess,
                                           void                     *cb_preprocess_opt,
                                           vrna_auxdata_free_f      cb_preprocess_opt_free)
{
  zarringhalam_options_t  *opt = (zarringhalam_options_t *)vrna_alloc(sizeof(zarringhalam_options_t));

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
    opt->cb_preprocess = set_mapping_strategy(VRNA_PROBING_METHOD_ZARRINGHALAM2012_DEFAULT_conversion,
                                              max_value,
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
                     double              max_value,
                     void                **transform_data,
                     vrna_auxdata_free_f *transform_data_free)
{
  vrna_probing_transform_f  transform_function;

  transform_function    = NULL;
  *transform_data       = NULL;
  *transform_data_free  = NULL;

  switch (*conversion_string) {
    case 'S':
      /* simple identity transformation where negative values are cut-off */
      transform_function = vrna_reactivity_trans_method(VRNA_REACTIVITY_TRANS_NEG_IGNORE);
      break;

    case 'M':
      /* default mapping from Zarringhalam et al. 2012 */
      {
        double map[5][2] = { { 0., 0.},
                             { 0.25, 0.35 },
                             { 0.3, 0.55 },
                             { 0.7, 0.85 },
                             { 0., 1.0}
                            };
        map[4][0] = max_value;
        transform_function = vrna_data_transform_method_bin((const double **)&(map[0][0]),
                                                            5,
                                                            VRNA_REACTIVITY_MISSING,
                                                            VRNA_REACTIVITY_MISSING,
                                                            VRNA_TRANSFORM_BIN_OPTION_PROJECT,
                                                            transform_data,
                                                            transform_data_free);
      }
      break;

    case 'C':
      /* cut-off mapping into two bins with values 0 and 1, also negative values are treated as missing */
      {
        float cutoff = 0.25;
        sscanf(conversion_string + 1, "%f", &cutoff);
        double map[2][2] = { { 0., 0.},
                             { 0., 1.0}
                            };
        map[1][0] = cutoff;
        transform_function = vrna_data_transform_method_bin((const double **)&(map[0][0]),
                                                            2,
                                                            VRNA_REACTIVITY_MISSING,
                                                            VRNA_REACTIVITY_MISSING,
                                                            VRNA_TRANSFORM_BIN_OPTION_MAP_OUTOF_UPPERBOUND,
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
        float s, i;

        /* try parsing other parameter settings */
        int r;

        r = sscanf(conversion_string + 1, "s%fi%f", &s, &i);

        if (r < 2) {
          s = slope;
          i = intercept;

          r = sscanf(conversion_string + 1, "s%f", &s);

          if (r < 1) {
            r = sscanf(conversion_string + 1, "i%f", &i);
            if (r < 1)
              vrna_log_warning("Probing method parameters not recognized! Using default parameters!");
          }
        }

        transform_function = vrna_data_transform_method_lm(s,
                                                           i,
                                                           VRNA_REACTIVITY_MISSING,
                                                           (*conversion_string == 'L') ? 0 : VRNA_TRANSFORM_LM_OPTION_LOG,
                                                           transform_data,
                                                           transform_data_free);
      }
      break;

    default:
      break;
  }

  return NULL;
}
