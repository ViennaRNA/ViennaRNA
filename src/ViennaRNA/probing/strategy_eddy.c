/* Probing data integration - Eddy 2014 approach */

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

#include "ViennaRNA/probing/strategy_eddy.h"

#define gaussian(u) (1 / (sqrt(2 * PI)) * exp(-u * u / 2))

typedef struct {
  unsigned char             options;
  double                    temperature;
  double                    kT;
  vrna_array(double)        priors_unpaired;
  vrna_array(double)        priors_paired;
  double                    unpaired_h;
  double                    paired_h;
  vrna_math_fun_dbl_f           cb_preprocess;
  vrna_math_fun_dbl_opt_t       cb_preprocess_opt;
  vrna_math_fun_dbl_opt_free_f  cb_preprocess_opt_free;
} eddy_options_t;


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
conversion_eddy_up(double         value,
                   eddy_options_t *options);


PRIVATE double
conversion_eddy_bp(double         value,
                   eddy_options_t *options);


/* PDF of x using Gaussian KDE
 * n is data number
 * h is bandwidth
 */
PRIVATE double
gaussian_kde_pdf(double       x,
                 unsigned int n,
                 float        h,
                 const double *data);


PRIVATE double
exp_pdf(double  x,
        double  lambda) VRNA_UNUSED;


/* We use same format as scitpy */
PRIVATE double
gev_pdf(double  x,
        double  c,
        double  loc,
        double  scale) VRNA_UNUSED;


/*
 * Bandwidth for univariate KDE with Scott factor as in scipy
 * bandwidth = Scott facter * std with ddof = 1
 */
PRIVATE double
bandwidth(unsigned int  n,
          const double  *data);


/*
 #################################
 # BEGIN OF FUNCTION DEFINITIONS #
 #################################
 */

PUBLIC double *
vrna_probing_strategy_eddy(vrna_fold_compound_t *fc,
                           const double         *data,
                           size_t               data_size,
                           unsigned int         target,
                           void                 *options)
{
  double          *pseudo_energies;
  eddy_options_t  *opt;

  pseudo_energies = NULL;

  if (((target == VRNA_PROBING_DATA_LINEAR_TARGET_UP) || (target == VRNA_PROBING_DATA_LINEAR_TARGET_BP)) &&
      (data) &&
      (data_size > 0)) {

    /* preprare (default) options */
    if (options) {
      opt = (eddy_options_t *)options;
    } else {
      opt = vrna_probing_strategy_eddy_options(fc->params->temperature,
                                               VRNA_PROBING_STRATEGY_EDDY_OPTIONS_DEFAULT,
                                               NULL,
                                               0,
                                               NULL,
                                               0,
                                               NULL,
                                               NULL,
                                               NULL);
    }

    /* pre-process data */
    pseudo_energies = vrna_data_lin_transform(data,
                                              data_size,
                                              opt->cb_preprocess,
                                              opt->cb_preprocess_opt,
                                              NULL,
                                              VRNA_REACTIVITY_MISSING,
                                              VRNA_TRANSFORM_DEFAULT);

    if (!(opt->options & VRNA_PROBING_STRATEGY_EDDY_NO_TEMPERATURE_RESCALING)) {
      opt->temperature = fc->params->temperature;
      opt->kT          = GASCONST * (opt->temperature + K0) / 1000; /* in kcal/mol */
    }

    /* transform data into actual pseudo-energies */
    if (target == VRNA_PROBING_DATA_LINEAR_TARGET_UP)
      for (size_t i = 0; i < data_size; i++)
        pseudo_energies[i] = conversion_eddy_up(pseudo_energies[i], opt);
    else
      for (size_t i = 0; i < data_size; i++)
        pseudo_energies[i] = conversion_eddy_bp(pseudo_energies[i], opt);

    /* release memory for default options */
    if (opt != (eddy_options_t *)options)
      vrna_probing_strategy_eddy_options_free(options);
  }

  return pseudo_energies;
}


PUBLIC void *
vrna_probing_strategy_eddy_options(double                         temperature,
                                   unsigned char                  options,
                                   const double                   *prior_unpaired,
                                   size_t                         prior_unpaired_size,
                                   const double                   *prior_paired,
                                   size_t                         prior_paired_size,
                                   vrna_math_fun_dbl_f          cb_preprocess,
                                   vrna_math_fun_dbl_opt_t      cb_preprocess_opt,
                                   vrna_math_fun_dbl_opt_free_f cb_preprocess_opt_free)
{
  const double  *ptr;
  double        *ptr_transformed;
  size_t        ptr_size;

  eddy_options_t  *opt = (eddy_options_t *)vrna_alloc(sizeof(eddy_options_t));

  opt->temperature    = temperature;
  opt->kT             = GASCONST * (temperature + K0) / 1000; /* in kcal/mol */
  opt->options        = options;

  if (cb_preprocess) {
    opt->cb_preprocess          = cb_preprocess;
    opt->cb_preprocess_opt      = cb_preprocess_opt;
    opt->cb_preprocess_opt_free = cb_preprocess_opt_free;
  } else {
    /* default preprocessing of the probing data */
    opt->cb_preprocess          = vrna_math_fun_dbl_log_opt(1,
                                                        0,
                                                        VRNA_REACTIVITY_MISSING,
                                                        VRNA_MATH_FUN_LOG_OPTION_DEFAULT,
                                                        &(opt->cb_preprocess_opt),
                                                        &(opt->cb_preprocess_opt_free));
  }

  /* append transformed priors */
  if ((prior_unpaired) &&
      (prior_unpaired_size > 0)) {
    ptr       = prior_unpaired;
    ptr_size  = prior_unpaired_size;
  } else {
    /* parse build-in priors */
    ptr       = vrna_str_to_dbl_array(probing_data_prior_probing_1M7_unpaired,
                                      "\n",
                                      0);
    ptr_size  = vrna_array_size(ptr);
  }

  ptr_transformed = vrna_data_lin_transform(ptr,
                                            ptr_size,
                                            opt->cb_preprocess,
                                            opt->cb_preprocess_opt,
                                            NULL,
                                            VRNA_REACTIVITY_MISSING,
                                            VRNA_TRANSFORM_DEFAULT);

  vrna_array_init_size(opt->priors_unpaired, ptr_size);

  for (size_t i = 0; i < ptr_size; i++)
    vrna_array_append(opt->priors_unpaired, ptr_transformed[i]);

  opt->unpaired_h = bandwidth(ptr_size, ptr_transformed);

  if (ptr != prior_unpaired)
    vrna_array_free(ptr);

  free(ptr_transformed);

  if ((prior_paired) &&
      (prior_paired_size > 0)) {
    ptr       = prior_paired;
    ptr_size  = prior_paired_size;
  } else {
    /* parse build-in priors */
    ptr       = vrna_str_to_dbl_array(probing_data_prior_probing_1M7_paired,
                                      "\n",
                                      0);
    ptr_size  = vrna_array_size(ptr);
  }

  ptr_transformed = vrna_data_lin_transform(ptr,
                                            ptr_size,
                                            opt->cb_preprocess,
                                            opt->cb_preprocess_opt,
                                            NULL,
                                            VRNA_REACTIVITY_MISSING,
                                            VRNA_TRANSFORM_DEFAULT);

  vrna_array_init_size(opt->priors_paired, ptr_size);

  for (size_t i = 0; i < ptr_size; i++)
    vrna_array_append(opt->priors_paired, ptr_transformed[i]);

  opt->paired_h = bandwidth(ptr_size, ptr_transformed);

  if (ptr != prior_paired)
    vrna_array_free(ptr);

  free(ptr_transformed);


  return (void *)opt;
}


PUBLIC void
vrna_probing_strategy_eddy_options_free(void *options)
{
  eddy_options_t  *opt = (eddy_options_t *)options;

  if (opt->cb_preprocess_opt_free)
    opt->cb_preprocess_opt_free(opt->cb_preprocess_opt);

  vrna_array_free(opt->priors_unpaired);
  vrna_array_free(opt->priors_paired);

  free(opt);
}


PUBLIC struct vrna_probing_data_s *
vrna_probing_data_eddy(const double   *reactivities,
                       unsigned int   n,
                       double         temperature,
                       unsigned char  options,
                       const double   *unpaired_data,
                       unsigned int   unpaired_len,
                       const double   *paired_data,
                       unsigned int   paired_len)
{
  return vrna_probing_data_eddy_trans(reactivities,
                                      n,
                                      temperature,
                                      options,
                                      unpaired_data,
                                      unpaired_len,
                                      paired_data,
                                      paired_len,
                                      NULL,
                                      NULL,
                                      NULL);
}


PUBLIC struct vrna_probing_data_s *
vrna_probing_data_eddy_trans(const double             *reactivities,
                             unsigned int             n,
                             double                   temperature,
                             unsigned char            options,
                             const double             *unpaired_data,
                             unsigned int             unpaired_len,
                             const double             *paired_data,
                             unsigned int             paired_len,
                             vrna_math_fun_dbl_f          trans,
                             vrna_math_fun_dbl_opt_t      trans_options,
                             vrna_math_fun_dbl_opt_free_f trans_options_free)
{
  if (reactivities) {
    return vrna_probing_data_linear(reactivities,
                                    n,
                                    1.,
                                    vrna_probing_strategy_eddy,
                                    vrna_probing_strategy_eddy_options(temperature,
                                                                       options,
                                                                       unpaired_data,
                                                                       unpaired_len,
                                                                       paired_data,
                                                                       paired_len,
                                                                       trans,
                                                                       trans_options,
                                                                       trans_options_free),
                                    vrna_probing_strategy_eddy_options_free);
  }

  return NULL;
}



PUBLIC struct vrna_probing_data_s *
vrna_probing_data_eddy_comparative(const double       **reactivities,
                                   const unsigned int *n,
                                   unsigned int       n_seq,
                                   double             temperature,
                                   unsigned char      options,
                                   const double       **unpaired_data,
                                   const unsigned int *unpaired_len,
                                   const double       **paired_data,
                                   const unsigned int *paired_len,
                                   unsigned int       multi_params)
{
  return vrna_probing_data_eddy_trans_comparative(reactivities,
                                                  n,
                                                  n_seq,
                                                  temperature,
                                                  options,
                                                  unpaired_data,
                                                  unpaired_len,
                                                  paired_data,
                                                  paired_len,
                                                  multi_params,
                                                  NULL,
                                                  NULL,
                                                  NULL);
}


PUBLIC struct vrna_probing_data_s *
vrna_probing_data_eddy_trans_comparative(const double **reactivities,
                                         const unsigned int *n,
                                         unsigned int n_seq,
                                         double       temperature,
                                         unsigned char  options,
                                         const double **unpaired_datas,
                                         const unsigned int *unpaired_lens,
                                         const double **paired_datas,
                                         const unsigned int *paired_lens,
                                         unsigned int multi_params,
                                         vrna_math_fun_dbl_f          *trans,
                                         vrna_math_fun_dbl_opt_t      *trans_options,
                                         vrna_math_fun_dbl_opt_free_f *trans_options_free)
{
  const double                          *priors_unpaired, *priors_paired;
  double                                priors_unpaired_size, priors_paired_size;
  struct vrna_probing_data_s            *d = NULL;
  vrna_math_fun_dbl_f                 cb_trans;
  vrna_math_fun_dbl_opt_t             cb_trans_options;
  vrna_math_fun_dbl_opt_free_f        cb_trans_options_free;
  vrna_array(vrna_probing_strategy_f)   cbs_linear;
  vrna_array(void *)                    cbs_linear_options;
  vrna_array(vrna_auxdata_free_f)       cbs_linear_options_free;

  if ((reactivities) &&
      (n)) {
    priors_unpaired       = (unpaired_datas) ? *unpaired_datas : NULL;
    priors_paired         = (paired_datas) ? *paired_datas : NULL;
    priors_unpaired_size  = (unpaired_datas) ? *unpaired_lens : 0;
    priors_paired_size    = (paired_datas) ? *paired_lens : 0;
    cb_trans              = (trans) ? *trans : NULL;
    cb_trans_options      = (trans_options) ? *trans_options : NULL;
    cb_trans_options_free = (trans_options_free) ? *trans_options_free : NULL;

    if (((unpaired_datas == NULL) && (multi_params & VRNA_PROBING_METHOD_MULTI_PARAMS_1)) ||
        ((paired_datas == NULL) && (multi_params & VRNA_PROBING_METHOD_MULTI_PARAMS_2)) ||
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
        if (multi_params & VRNA_PROBING_METHOD_MULTI_PARAMS_1) {
          priors_unpaired       = unpaired_datas[s];
          priors_unpaired_size  = unpaired_lens[s];
        }

        if (multi_params & VRNA_PROBING_METHOD_MULTI_PARAMS_2) {
          priors_paired       = paired_datas[s];
          priors_paired_size  = paired_lens[s];
        }

        if (multi_params & VRNA_PROBING_METHOD_MULTI_PARAMS_3) {
          cb_trans = trans[s];
          if (trans_options)
            cb_trans_options = trans_options[s];
          if (trans_options_free)
            cb_trans_options_free = trans_options_free[s];
        }

        vrna_array_append(cbs_linear, vrna_probing_strategy_eddy);
        vrna_array_append(cbs_linear_options, vrna_probing_strategy_eddy_options(temperature,
                                                                                 options,
                                                                                 priors_unpaired,
                                                                                 priors_unpaired_size,
                                                                                 priors_paired,
                                                                                 priors_paired_size,
                                                                                 cb_trans,
                                                                                 cb_trans_options,
                                                                                 cb_trans_options_free));
        vrna_array_append(cbs_linear_options_free, vrna_probing_strategy_eddy_options_free);
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
conversion_eddy_up(double       value,
                   eddy_options_t *options)
{
  FLT_OR_DBL p;
  if (value == VRNA_REACTIVITY_MISSING) {
    return 0;
  } else {
    p = gaussian_kde_pdf(value,
                         vrna_array_size(options->priors_unpaired),
                         options->unpaired_h,
                         options->priors_unpaired);
    p = (p > 1e-10) ? p : 1e-10;
    return  -options->kT * log(p);
  }
}


PRIVATE double
conversion_eddy_bp(double       value,
                   eddy_options_t *options)
{
  FLT_OR_DBL p;
  if (value == VRNA_REACTIVITY_MISSING) {
    return 0;
  } else {
    p = gaussian_kde_pdf(value,
                         vrna_array_size(options->priors_paired),
                         options->paired_h,
                         options->priors_paired);
    p = (p > 1e-10) ? p : 1e-10;
    return  -options->kT * log(p);
  }
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


