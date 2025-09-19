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
  double                    m;
  double                    b;
  vrna_probing_transform_f  cb_preprocess;
  void                      *cb_preprocess_opt;
  vrna_auxdata_free_f       cb_preprocess_opt_free;
} deigan_options_t;

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
conversion_deigan(double  reactivity,
                  double  m,
                  double  b);


/*
 #################################
 # BEGIN OF FUNCTION DEFINITIONS #
 #################################
 */


PUBLIC double *
vrna_probing_strategy_deigan(const double *data,
                             size_t       data_size,
                             unsigned int target,
                             void         *options)
{
  double            *pseudo_energies;
  deigan_options_t  *opt;

  pseudo_energies = NULL;

  if ((target == VRNA_PROBING_LINEAR_TARGET_STACK) &&
      (data) &&
      (data_size > 0)) {

    /* preprare (default) options */
    if (options)
      opt = (deigan_options_t *)options;
    else
      opt = vrna_probing_strategy_deigan_options(VRNA_PROBING_METHOD_DEIGAN2009_DEFAULT_m,
                                                 VRNA_PROBING_METHOD_DEIGAN2009_DEFAULT_b,
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
      pseudo_energies[i] = conversion_deigan(pseudo_energies[i],
                                             opt->m,
                                             opt->b);

    /* release memory for default options */
    if (opt != (deigan_options_t *)options)
      free(opt);
  }

  return pseudo_energies;
}


PUBLIC void *
vrna_probing_strategy_deigan_options(double                   m,
                                     double                   b,
                                     vrna_probing_transform_f cb_preprocess,
                                     void                     *cb_preprocess_opt,
                                     vrna_auxdata_free_f      cb_preprocess_opt_free)
{
  deigan_options_t  *opt = (deigan_options_t *)vrna_alloc(sizeof(deigan_options_t));

  opt->m = m;
  opt->b = b;

  if (cb_preprocess) {
    opt->cb_preprocess          = cb_preprocess;
    opt->cb_preprocess_opt      = cb_preprocess_opt;
    opt->cb_preprocess_opt_free = cb_preprocess_opt_free;
  } else {
    opt->cb_preprocess          = vrna_reactivity_trans_method(VRNA_REACTIVITY_TRANS_NEG_IGNORE);
    opt->cb_preprocess_opt      = NULL;
    opt->cb_preprocess_opt_free = NULL;
  }

  return (void *)opt;
}


PUBLIC void
vrna_probing_strategy_deigan_options_free(void *options)
{
  deigan_options_t  *opt = (deigan_options_t *)options;

  if (opt->cb_preprocess_opt_free)
    opt->cb_preprocess_opt_free(opt->cb_preprocess_opt);

  free(opt);
}


PUBLIC struct vrna_probing_data_s *
vrna_probing_data_deigan(const double             *reactivities,
                         unsigned int             n,
                         double                   m,
                         double                   b)
{
  return vrna_probing_data_deigan_trans(reactivities,
                                        n,
                                        m,
                                        b,
                                        NULL,
                                        NULL,
                                        NULL);
}


PUBLIC struct vrna_probing_data_s *
vrna_probing_data_deigan_trans(const double             *reactivities,
                               unsigned int             n,
                               double                   m,
                               double                   b,
                               vrna_probing_transform_f trans,
                               void                     *trans_options,
                               vrna_auxdata_free_f      trans_options_free)
{
  if (reactivities) {
    return vrna_probing_data_linear(reactivities,
                                    n,
                                    1.,
                                    vrna_probing_strategy_deigan,
                                    vrna_probing_strategy_deigan_options(m, b, trans, trans_options, trans_options_free),
                                    vrna_probing_strategy_deigan_options_free);
  }

  return NULL;
}


PUBLIC struct vrna_probing_data_s *
vrna_probing_data_deigan_comparative(const double       **reactivities,
                                         const unsigned int *n,
                                         unsigned int       n_seq,
                                         double             *ms,
                                         double             *bs,
                                         unsigned int       multi_params)
{
  return vrna_probing_data_deigan_trans_comparative(reactivities,
                                                    n,
                                                    n_seq,
                                                    ms,
                                                    bs,
                                                    multi_params,
                                                    NULL,
                                                    NULL,
                                                    NULL);

}

PUBLIC struct vrna_probing_data_s *
vrna_probing_data_deigan_trans_comparative(const double       **reactivities,
                                          const unsigned int *n,
                                         unsigned int       n_seq,
                                         double             *ms,
                                         double             *bs,
                                         unsigned int       multi_params,
                                         vrna_probing_transform_f trans,
                                         void                     *trans_options,
                                         vrna_auxdata_free_f      trans_options_free)
{
  struct vrna_probing_data_s  *d = NULL;
  double                      m, b;
  vrna_array(vrna_probing_strategy_f)   cbs_linear;
  vrna_array(void *)                    cbs_linear_options;
  vrna_array(vrna_auxdata_free_f)       cbs_linear_options_free;

  if ((reactivities) &&
      (n)) {

    m = (ms) ? *ms : VRNA_PROBING_METHOD_DEIGAN2009_DEFAULT_m;
    b = (bs) ? *bs : VRNA_PROBING_METHOD_DEIGAN2009_DEFAULT_b;

    if (((ms == NULL) && (multi_params & VRNA_PROBING_METHOD_MULTI_PARAMS_1)) ||
        ((bs == NULL) && (multi_params & VRNA_PROBING_METHOD_MULTI_PARAMS_2))) {
      /* if multi_params != 0, either ms or bs or both must be provided! */
      return d;
    }

    /* prepare callback arrays */
    vrna_array_init_size(cbs_linear, n_seq);
    vrna_array_init_size(cbs_linear_options, n_seq);
    vrna_array_init_size(cbs_linear_options_free, n_seq);

    /* fill callback vectors */
    for (size_t i = 0; i < n_seq; i++) {
      if (reactivities[i]) {
        if (multi_params & VRNA_PROBING_METHOD_MULTI_PARAMS_1)
          m = ms[i];

        if (multi_params & VRNA_PROBING_METHOD_MULTI_PARAMS_2)
          b = bs[i];

        vrna_array_append(cbs_linear, vrna_probing_strategy_deigan);
        vrna_array_append(cbs_linear_options, vrna_probing_strategy_deigan_options(m, b, trans, trans_options, trans_options_free));
        vrna_array_append(cbs_linear_options_free, free);
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
conversion_deigan(double  reactivity,
                  double  m,
                  double  b)
{
  return reactivity == VRNA_REACTIVITY_MISSING ? 0. : (double)(m * log(reactivity + 1) + b);
}
