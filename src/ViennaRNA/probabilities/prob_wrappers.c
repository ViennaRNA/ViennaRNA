#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "ViennaRNA/utils/basic.h"
#include "ViennaRNA/probabilities/basic.h"
#include "ViennaRNA/probabilities/basepairs.h"

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
typedef struct {
  vrna_array(vrna_ep_t) probs;
  double                cutoff;
} up_storage_default;


PRIVATE void
accessibility_cb_default(vrna_fold_compound_t   *fc,
                         FLT_OR_DBL             *pr,
                         unsigned int           i,
                         unsigned int           a,
                         unsigned int           b,
                         unsigned int           type,
                         vrna_eqm_probs_data_t  data);



/*
 #################################
 # BEGIN OF FUNCTION DEFINITIONS #
 #################################
 */
PUBLIC vrna_array(vrna_ep_t)
vrna_eqm_probs_up(vrna_fold_compound_t *fc,
                  unsigned int         max_u,
                  double               cutoff)
{
  vrna_ep_t *result = NULL;

  if (fc) {
    if (max_u == 0)
      max_u = fc->length;

    vrna_eqm_probs_opt_t o = {
      .mode       = VRNA_EQM_PROBS_UP,
      .up_length  = max_u
    };

    up_storage_default payload = {
      .probs = NULL,
      .cutoff = cutoff
    };

    vrna_array_init_size(payload.probs, fc->length);

    if (vrna_eqm_probabilities(fc,
                               accessibility_cb_default,
                               (void *)&payload,
                               o)) {
      result = payload.probs;
    } else {
      vrna_array_free(payload.probs);
    }
  }

  return result;
}


/*
 #################################
 # STATIC helper functions below #
 #################################
 */
PRIVATE void
accessibility_cb_default(vrna_fold_compound_t   *fc,
                         FLT_OR_DBL             *pr,
                         unsigned int           i,
                         unsigned int           a,
                         unsigned int           b,
                         unsigned int           type,
                         vrna_eqm_probs_data_t  data)
{
  if ((type & VRNA_EQM_PROBS_UP) &&
      (type & VRNA_LOOP_TYPE_ANY)) {
    for (unsigned int k = a; k <= b; k++)
      if (pr[k] >= ((up_storage_default *)data)->cutoff)
        vrna_array_append(((up_storage_default *)data)->probs,
                          ((vrna_ep_t){
                            .i = i,
                            .j = k,
                            .p = pr[k],
                            .type = VRNA_PLIST_TYPE_UNPAIRED
                          }));
  }
}


