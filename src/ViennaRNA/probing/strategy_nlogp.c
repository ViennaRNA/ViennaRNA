/* Probing data integration - negative log probability approach */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "ViennaRNA/utils/basic.h"
#include "ViennaRNA/utils/strings.h"
#include "ViennaRNA/datastructures/array.h"
#include "ViennaRNA/static/probing_data_priors.h"
#include "ViennaRNA/data/transform.h"

#include "ViennaRNA/probing/strategy_nlogp.h"


typedef struct {
  unsigned int                  options;
  double                        base;
  double                        factor;
  double                        factor_t;
  double                        min_prob;
  vrna_math_fun_dbl_f           cb_preprocess;
  vrna_math_fun_dbl_opt_t       cb_preprocess_opt;
  vrna_math_fun_dbl_opt_free_f  cb_preprocess_opt_free;
} nlogp_options_t;


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
PRIVATE double
conversion_nlogp(double           value,
                 nlogp_options_t  *options);


/*
 #################################
 # BEGIN OF FUNCTION DEFINITIONS #
 #################################
 */
PUBLIC double *
vrna_probing_strategy_nlogp(vrna_fold_compound_t  *fc,
                            const double          *data,
                            size_t                data_size,
                            unsigned int          target,
                            void                  *options)
{
  double          *pseudo_energies, t;
  nlogp_options_t *opt;

  pseudo_energies = NULL;

  if ((data) &&
      (data_size > 0)) {
    /* preprare (default) options */
    if (options) {
      opt = (nlogp_options_t *)options;
    } else {
      /* use default options */
      opt = vrna_probing_strategy_nlogp_options(0.,
                                                1.,
                                                VRNA_PROBING_STRATEGY_NLOGP_OPTIONS_DEFAULT,
                                                NULL,
                                                NULL,
                                                NULL);
    }

    if (opt->options & target) {
      /* pre-process data */

      /*  we assume to obtain probabilities here, so let us make sure that
       *  the result behaves nicely and doesn't exceed the [0,1] range of
       *  values. We do so by clipping any data exceeding the domain limits
       */
      double domain[4] = {
        0., 0., 0., 1.
      };

      pseudo_energies = vrna_data_lin_transform(data,
                                                data_size,
                                                opt->cb_preprocess,
                                                opt->cb_preprocess_opt,
                                                domain,
                                                VRNA_REACTIVITY_MISSING,
                                                VRNA_TRANSFORM_ENFORCE_DOMAIN_TARGET);

      if (opt->options & VRNA_PROBING_STRATEGY_NLOGP_FACTOR_KELVIN)
        opt->factor_t = fc->params->temperature + K0;
      else if (opt->options & VRNA_PROBING_STRATEGY_NLOGP_FACTOR_CELCIUS)
        opt->factor_t = fc->params->temperature;
      else
        opt->factor_t = 1.;

      /* transform data into actual pseudo-energies */
      for (size_t i = 0; i < data_size; i++)
        pseudo_energies[i] = conversion_nlogp(pseudo_energies[i], opt);
    }

    /* release memory for default options */
    if (opt != (nlogp_options_t *)options)
      vrna_probing_strategy_nlogp_options_free(options);
  }

  return pseudo_energies;
}


PUBLIC void *
vrna_probing_strategy_nlogp_options(double                        base,
                                    double                        factor,
                                    unsigned int                  options,
                                    vrna_math_fun_dbl_f           cb_preprocess,
                                    vrna_math_fun_dbl_opt_t       cb_preprocess_opt,
                                    vrna_math_fun_dbl_opt_free_f  cb_preprocess_opt_free)
{
  const double    *ptr;
  double          *ptr_transformed;
  size_t          ptr_size;

  nlogp_options_t *opt = (nlogp_options_t *)vrna_alloc(sizeof(nlogp_options_t));

  opt->base     = base;
  opt->factor   = factor;
  opt->factor_t = 1.;
  opt->options  = options;

  opt->cb_preprocess          = cb_preprocess;
  opt->cb_preprocess_opt      = cb_preprocess_opt;
  opt->cb_preprocess_opt_free = cb_preprocess_opt_free;

  return (void *)opt;
}


PUBLIC void
vrna_probing_strategy_nlogp_options_free(void *options)
{
  nlogp_options_t *opt = (nlogp_options_t *)options;

  if (opt->cb_preprocess_opt_free)
    opt->cb_preprocess_opt_free(opt->cb_preprocess_opt);

  free(opt);
}


PUBLIC struct vrna_probing_data_s *
vrna_probing_data_nlogp(const double  *reactivities,
                        unsigned int  n,
                        double        base,
                        double        factor,
                        unsigned int  options)
{
  return vrna_probing_data_nlogp_trans(reactivities,
                                       n,
                                       base,
                                       factor,
                                       options,
                                       NULL,
                                       NULL,
                                       NULL);
}


PUBLIC struct vrna_probing_data_s *
vrna_probing_data_nlogp_trans(const double                  *reactivities,
                              unsigned int                  n,
                              double                        base,
                              double                        factor,
                              unsigned int                  options,
                              vrna_math_fun_dbl_f           trans,
                              vrna_math_fun_dbl_opt_t       trans_options,
                              vrna_math_fun_dbl_opt_free_f  trans_options_free)
{
  if (reactivities) {
    return vrna_probing_data_linear(reactivities,
                                    n,
                                    NULL,
                                    vrna_probing_strategy_nlogp,
                                    vrna_probing_strategy_nlogp_options(base,
                                                                        factor,
                                                                        options,
                                                                        trans,
                                                                        trans_options,
                                                                        trans_options_free),
                                    vrna_probing_strategy_nlogp_options_free,
                                    VRNA_PROBING_DATA_DEFAULT);
  }

  return NULL;
}


PUBLIC struct vrna_probing_data_s *
vrna_probing_data_nlogp_comparative(const double        **reactivities,
                                    const unsigned int  *n,
                                    unsigned int        n_seq,
                                    double              base,
                                    double              factor,
                                    unsigned int        options,
                                    unsigned int        multi_params)
{
  return vrna_probing_data_nlogp_trans_comparative(reactivities,
                                                   n,
                                                   n_seq,
                                                   base,
                                                   factor,
                                                   options,
                                                   multi_params,
                                                   NULL,
                                                   NULL,
                                                   NULL);
}


PUBLIC struct vrna_probing_data_s *
vrna_probing_data_nlogp_trans_comparative(const double                  **reactivities,
                                          const unsigned int            *n,
                                          unsigned int                  n_seq,
                                          double                        base,
                                          double                        factor,
                                          unsigned int                  options,
                                          unsigned int                  multi_params,
                                          vrna_math_fun_dbl_f           *trans,
                                          vrna_math_fun_dbl_opt_t       *trans_options,
                                          vrna_math_fun_dbl_opt_free_f  *trans_options_free)
{
  struct vrna_probing_data_s    *d = NULL;
  vrna_math_fun_dbl_f           cb_trans;
  vrna_math_fun_dbl_opt_t       cb_trans_options;
  vrna_math_fun_dbl_opt_free_f  cb_trans_options_free;

  vrna_array(vrna_probing_strategy_f)   cbs_linear;
  vrna_array(void *)                    cbs_linear_options;
  vrna_array(vrna_auxdata_free_f)       cbs_linear_options_free;

  if ((reactivities) &&
      (n)) {
    cb_trans              = (trans) ? *trans : NULL;
    cb_trans_options      = (trans_options) ? *trans_options : NULL;
    cb_trans_options_free = (trans_options_free) ? *trans_options_free : NULL;

    if ((trans == NULL) && (multi_params & VRNA_PROBING_METHOD_MULTI_PARAMS_1)) {
      return d;
    }

    /* prepare callback arrays */
    vrna_array_init_size(cbs_linear, n_seq);
    vrna_array_init_size(cbs_linear_options, n_seq);
    vrna_array_init_size(cbs_linear_options_free, n_seq);

    /* fill callback vectors */
    for (size_t s = 0; s < n_seq; s++) {
      if (reactivities[s]) {
        if (multi_params & VRNA_PROBING_METHOD_MULTI_PARAMS_1) {
          cb_trans = trans[s];
          if (trans_options)
            cb_trans_options = trans_options[s];

          if (trans_options_free)
            cb_trans_options_free = trans_options_free[s];
        }

        vrna_array_append(cbs_linear, vrna_probing_strategy_nlogp);
        vrna_array_append(cbs_linear_options, vrna_probing_strategy_nlogp_options(base,
                                                                                  factor,
                                                                                  options,
                                                                                  cb_trans,
                                                                                  cb_trans_options,
                                                                                  cb_trans_options_free));
        vrna_array_append(cbs_linear_options_free, vrna_probing_strategy_nlogp_options_free);
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
PRIVATE double
conversion_nlogp(double           value,
                 nlogp_options_t  *options)
{
  double p;

  if (value == VRNA_REACTIVITY_MISSING) {
    return 0;
  } else {
    if (value < options->min_prob)
      value = options->min_prob;

    p = -options->factor * log(value);

    if (options->options & VRNA_PROBING_STRATEGY_NLOGP_NONSTD_BASE)
      p /= log(options->base);

    return p;
  }
}
