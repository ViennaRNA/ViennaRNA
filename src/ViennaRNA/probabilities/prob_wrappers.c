#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>    /* #defines FLT_MAX ... */

#include "ViennaRNA/utils/basic.h"
#include "ViennaRNA/eval/internal.h"
#include "ViennaRNA/probabilities/basic.h"
#include "ViennaRNA/probabilities/basepairs.h"
#include "ViennaRNA/probabilities/helices.h"
#include "ViennaRNA/probabilities/accessibility.h"

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
  unsigned int          min_size;
} storage_default;


PRIVATE void
bpp_cb_default(vrna_fold_compound_t   *fc,
               FLT_OR_DBL             *pr,
               unsigned int           i,
               unsigned int           a,
               unsigned int           b,
               unsigned int           type,
               vrna_eqm_probs_data_t  data);


PRIVATE void
accessibility_cb_default(vrna_fold_compound_t   *fc,
                         FLT_OR_DBL             *pr,
                         unsigned int           i,
                         unsigned int           a,
                         unsigned int           b,
                         unsigned int           type,
                         vrna_eqm_probs_data_t  data);


PRIVATE void
hx_prob_cb_default(vrna_fold_compound_t   *fc,
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
vrna_eqm_probs_bp(vrna_fold_compound_t *fc,
                  double               cutoff)
{
  vrna_ep_t *result = NULL;

  if (fc) {
    vrna_eqm_probs_opt_t o = {
      .mode = VRNA_EQM_PROBS_BP
    };

    storage_default payload = {
      .probs  = NULL,
      .cutoff = cutoff
    };

    vrna_array_init_size(payload.probs, fc->length);

    if (vrna_eqm_probabilities(fc,
                               bpp_cb_default,
                               (void *)&payload,
                               o)) {
      result = payload.probs;
    } else {
      vrna_array_free(payload.probs);
    }
  }

  return result;
}


PUBLIC vrna_array(vrna_ep_t)
vrna_eqm_probs_up(vrna_fold_compound_t *fc,
                  unsigned int         min_u,
                  unsigned int         max_u,
                  double               cutoff)
{
  vrna_ep_t *result = NULL;

  if (fc) {
    if (max_u == 0)
      max_u = fc->length;

    if (min_u > max_u) {
      vrna_log_warning("Reducing min_u = %u to %u",
                       min_u,
                       max_u);
      min_u = max_u;
    }

    vrna_eqm_probs_opt_t o = {
      .mode       = VRNA_EQM_PROBS_UP,
      .up_length  = max_u
    };

    storage_default payload = {
      .probs  = NULL,
      .cutoff = cutoff,
      .min_size = (min_u) ? min_u : 1
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


PUBLIC vrna_array(vrna_ep_t)
vrna_eqm_probs_hx(vrna_fold_compound_t *fc,
                  unsigned int         min_size,
                  unsigned int         max_size,
                  double               cutoff)
{
  vrna_ep_t *result = NULL;

  if (fc) {
    if (max_size == 0)
      max_size = fc->length / 2;

    if (min_size > max_size) {
      vrna_log_warning("Reducing min_size = %u to %u",
                       min_size,
                       max_size);
      min_size = max_size;
    }

    vrna_eqm_probs_opt_t o = {
      .mode       = VRNA_EQM_PROBS_HX,
      .hx_length  = max_size
    };

    storage_default payload = {
      .probs    = NULL,
      .cutoff   = cutoff,
      .min_size = (min_size) ? min_size : 1
    };

    vrna_array_init_size(payload.probs, fc->length);

    if (vrna_eqm_probabilities(fc,
                               hx_prob_cb_default,
                               (void *)&payload,
                               o)) {
      result = payload.probs;
    } else {
      vrna_array_free(payload.probs);
    }
  }

  return result;
}


PUBLIC vrna_ep_t *
vrna_stack_prob(vrna_fold_compound_t  *vc,
                double                cutoff)
{
  vrna_ep_t         *pl;
  int               i, j, plsize, length, *index, *jindx, *rtype, num;
  char              *ptype;
  FLT_OR_DBL        *qb, *probs, *scale, p;
  vrna_exp_param_t  *pf_params;
  vrna_mx_pf_t      *matrices;

  plsize  = 256;
  pl      = NULL;
  num     = 0;

  if (vc) {
    pf_params = vc->exp_params;
    length    = vc->length;
    index     = vc->iindx;
    jindx     = vc->jindx;
    rtype     = &(pf_params->model_details.rtype[0]);
    ptype     = vc->ptype;
    matrices  = vc->exp_matrices;
    qb        = matrices->qb;
    probs     = matrices->probs;
    scale     = matrices->scale;

    pl = (vrna_ep_t *)vrna_alloc(plsize * sizeof(vrna_ep_t));

    for (i = 1; i < length; i++)
      for (j = i + 3; j <= length; j++) {
        if ((p = probs[index[i] - j]) < cutoff)
          continue;

        if (qb[index[i + 1] - (j - 1)] < FLT_MIN)
          continue;

        p *= qb[index[i + 1] - (j - 1)] /
             qb[index[i] - j];
        p *= vrna_exp_eval_stack(vc, i, j, VRNA_EVAL_LOOP_DEFAULT);

        if (p > cutoff) {
          pl[num].i     = i;
          pl[num].j     = j;
          pl[num].type  = VRNA_PLIST_TYPE_STACK;
          pl[num++].p   = p;
          if (num >= plsize) {
            plsize  *= 2;
            pl      = vrna_realloc(pl, plsize * sizeof(vrna_ep_t));
          }
        }
      }

    pl[num].i = 0;
  }

  return pl;
}


/*
 #################################
 # STATIC helper functions below #
 #################################
 */
PRIVATE void
bpp_cb_default(vrna_fold_compound_t   *fc,
               FLT_OR_DBL             *pr,
               unsigned int           j,
               unsigned int           start,
               unsigned int           end,
               unsigned int           type,
               vrna_eqm_probs_data_t  data)
{
  if ((type & VRNA_EQM_PROBS_BP) &&
      (type & VRNA_LOOP_TYPE_ANY)) {
    storage_default *d = (storage_default *)data;

    for (unsigned int i = start; i <= end; i++)
      if (pr[i] >= d->cutoff)
        vrna_array_append(d->probs,
                          ((vrna_ep_t){
                            .i    = i,
                            .j    = j,
                            .size = 1,
                            .p    = pr[i],
                            .type = VRNA_PLIST_TYPE_BASEPAIR
                          }));
  }
}


PRIVATE void
accessibility_cb_default(vrna_fold_compound_t   *fc,
                         FLT_OR_DBL             *pr,
                         unsigned int           j,
                         unsigned int           start,
                         unsigned int           end,
                         unsigned int           type,
                         vrna_eqm_probs_data_t  data)
{
  if ((type & VRNA_EQM_PROBS_UP) &&
      (type & VRNA_LOOP_TYPE_ANY)) {
    storage_default *d = (storage_default *)data;

    for (unsigned int u = d->min_size; u <= end; u++)
      if (pr[u] >= d->cutoff)
        vrna_array_append(d->probs,
                          ((vrna_ep_t){
                            .i    = j + 1 - u,
                            .j    = j,
                            .size = u,
                            .p    = pr[u],
                            .type = VRNA_PLIST_TYPE_UNPAIRED
                          }));
  }
}


PRIVATE void
hx_prob_cb_default(vrna_fold_compound_t   *fc,
                   FLT_OR_DBL             *pr,
                   unsigned int           j,
                   unsigned int           i,
                   unsigned int           max_size,
                   unsigned int           type,
                   vrna_eqm_probs_data_t  data)
{
  if (type & VRNA_EQM_PROBS_HX) {
    storage_default *d = (storage_default *)data;

    if ((d->min_size == 1) &&
        (pr[1] >= d->cutoff))
      vrna_array_append(d->probs,
                        ((vrna_ep_t){
                          .i    = i,
                          .j    = j,
                          .size = 1,
                          .p    = pr[1],
                          .type = VRNA_PLIST_TYPE_BASEPAIR
                        }));

    if ((d->min_size <= 2) &&
        (pr[2] >= d->cutoff))
      vrna_array_append(d->probs,
                        ((vrna_ep_t){
                          .i    = i,
                          .j    = j,
                          .size = 2,
                          .p    = pr[2],
                          .type = VRNA_PLIST_TYPE_STACK
                        }));

    for (unsigned int s = MAX2(3, d->min_size); s <= max_size; s++)
      if (pr[s] >= d->cutoff)
        vrna_array_append(d->probs,
                          ((vrna_ep_t){
                            .i    = i,
                            .j    = j,
                            .size = s,
                            .p    = pr[s],
                            .type = VRNA_PLIST_TYPE_HELIX
                          }));
  }
}


