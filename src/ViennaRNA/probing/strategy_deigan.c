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
                             void         *options)
{
  double            *pseudo_energies;
  deigan_options_t  *opt;

  pseudo_energies = NULL;

  if ((data) &&
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
