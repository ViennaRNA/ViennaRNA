/*
 *                partiton function for RNA secondary structures
 *
 *                Ivo L Hofacker + Ronny Lorenz
 *                Vienna RNA package
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>    /* #defines FLT_MAX ... */
#include <limits.h>

#include "ViennaRNA/utils/basic.h"
#include "ViennaRNA/utils/log.h"
#include "ViennaRNA/params/default.h"
#include "ViennaRNA/fold_vars.h"
#include "ViennaRNA/eval/exterior.h"
#include "ViennaRNA/eval/hairpin.h"
#include "ViennaRNA/eval/internal.h"
#include "ViennaRNA/eval/multibranch.h"
#include "ViennaRNA/eval/gquad.h"
#include "ViennaRNA/constraints/hard.h"
#include "ViennaRNA/constraints/soft.h"
#include "ViennaRNA/mfe/global.h"
#include "ViennaRNA/eval/structures.h"
#include "ViennaRNA/sequences/alphabet.h"
#include "ViennaRNA/partfunc/global.h"
#include "ViennaRNA/probabilities/basic.h"
#include "ViennaRNA/probabilities/basepairs.h"

#include "ViennaRNA/constraints/exterior_sc_pf.inc"
#include "ViennaRNA/constraints/hairpin_sc_pf.inc"
#include "ViennaRNA/constraints/internal_sc_pf.inc"
#include "ViennaRNA/constraints/multibranch_sc_pf.inc"

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
typedef struct helper_arrays *helper_arrays_t;

typedef void  (*outside_f)(vrna_fold_compound_t *fc,
                           unsigned int         l,
                           FLT_OR_DBL           *Qmax,
                           int                  *ov,
                           helper_arrays_t      dp_helpers);

typedef void  (*component_f)(vrna_fold_compound_t *fc,
                             unsigned int         l,
                             FLT_OR_DBL           *Qmax,
                             int                  *ov,
                             void                 *data);

typedef void (*post_f)(vrna_fold_compound_t *fc,
                       unsigned int         l,
                       void                 *data);

typedef void (*pre_f)(vrna_fold_compound_t *fc,
                       void                 *data);


struct helper_arrays {
  unsigned int  n;
  unsigned int  strands;

  struct sc_ext_exp_dat     sc_wrapper_ext;
  struct sc_hp_exp_dat      sc_wrapper_hp;
  struct sc_int_exp_dat     sc_wrapper_int;
  struct sc_mb_exp_dat      sc_wrapper_mb;

  outside_f     run_outside_external;
  outside_f     run_outside_internal;
  outside_f     run_outside_multibranch;

  FLT_OR_DBL    *probs_ext;
  FLT_OR_DBL    *probs_hp;
  FLT_OR_DBL    *probs_int;
  FLT_OR_DBL    *probs_mb;

  FLT_OR_DBL    *prm_l;
  FLT_OR_DBL    *prm_l1;
  FLT_OR_DBL    *prml;

  FLT_OR_DBL    *Y5;
  FLT_OR_DBL    **Y5p;
  FLT_OR_DBL    **Y3;
  FLT_OR_DBL    **Y3p;

  FLT_OR_DBL    *qbh;
  FLT_OR_DBL    *access_ext;
  FLT_OR_DBL    *access_hp;
  FLT_OR_DBL    *access_int;
  FLT_OR_DBL    *access_mb;

  FLT_OR_DBL    *hxp;

  unsigned int  ud_max_size;
  FLT_OR_DBL    **pmlu;
  FLT_OR_DBL    *prm_MLbu;
  FLT_OR_DBL    *qm2_aux;

  vrna_eqm_probs_f       cb;
  vrna_eqm_probs_data_t  cb_data;
  vrna_eqm_probs_opt_t   options;

  vrna_array(vrna_ep_t)  bp_correction;

  vrna_array(component_f)         component_cbs;
  vrna_array(void *)              component_cbs_data;
  vrna_array(vrna_auxdata_free_f) component_cbs_data_free;
  vrna_array(pre_f)               component_cbs_init;
  vrna_array(post_f)              component_cbs_post;

  vrna_array(pre_f)       init_cbs;
  vrna_array(post_f)      post_cbs;
  vrna_array(void *)      init_cbs_data;
  vrna_array(void *)      post_cbs_data;
};


/* default options for pf_create_bppm() */
static vrna_eqm_probs_opt_t default_options = {
  .mode = VRNA_EQM_PROBS_BPP
};


PRIVATE void
init_probs_dp_mx(vrna_fold_compound_t *fc,
                 void                 *data);


PRIVATE void
init_bpp(vrna_fold_compound_t *fc,
         void                 *data);


PRIVATE void
report_n_clear_bpp(vrna_fold_compound_t *fc,
                   unsigned int         l,
                   void                 *data);


PRIVATE void
init_accessibility(vrna_fold_compound_t *fc,
                   void                 *data);


PRIVATE void
report_n_clear_accessibility(vrna_fold_compound_t *fc,
                             unsigned int         l,
                             void                 *data);


PRIVATE int
pf_create_bppm(vrna_fold_compound_t   *vc,
               vrna_eqm_probs_f       cb,
               vrna_eqm_probs_data_t  cb_data,
               vrna_eqm_probs_opt_t   options);


PRIVATE void
outside_bp_circ(vrna_fold_compound_t  *fc,
                unsigned int          l,
                FLT_OR_DBL            *Qmax,
                int                   *ov,
                helper_arrays_t       dp_helpers);


PRIVATE void
outside_bp_external(vrna_fold_compound_t  *fc,
                     unsigned int          l,
                     FLT_OR_DBL            *Qmax,
                     int                   *ov,
                     helper_arrays_t       dp_helpers);


PRIVATE void
outside_bp_enclosed(vrna_fold_compound_t  *vc,
                    unsigned int          l,
                    FLT_OR_DBL            *Qmax,
                    int                   *ov,
                    helper_arrays_t       dp_helpers);


PRIVATE void
compute_accessibility(vrna_fold_compound_t  *fc,
                      unsigned int          l,
                      helper_arrays_t       dp_helpers);


PRIVATE void
compute_hx_prob(vrna_fold_compound_t  *fc,
                unsigned int          l,
                helper_arrays_t       dp_helpers);


PRIVATE INLINE void
ud_outside_ext_loops(vrna_fold_compound_t *vc);


PRIVATE INLINE void
ud_outside_hp_loops(vrna_fold_compound_t *vc);


PRIVATE INLINE void
ud_outside_int_loops(vrna_fold_compound_t *vc);


PRIVATE INLINE void
ud_outside_mb_loops(vrna_fold_compound_t *vc);


PRIVATE INLINE void
ud_outside_hp_loops2(vrna_fold_compound_t *vc);


PRIVATE INLINE void
ud_outside_int_loops2(vrna_fold_compound_t *vc);


PRIVATE INLINE void
ud_outside_mb_loops2(vrna_fold_compound_t *vc);


PRIVATE FLT_OR_DBL
numerator_single(vrna_fold_compound_t *vc,
                 int                  i,
                 int                  j);


PRIVATE FLT_OR_DBL
numerator_comparative(vrna_fold_compound_t  *vc,
                      int                   i,
                      int                   j);


PRIVATE helper_arrays_t
get_dp_helper_arrays(vrna_fold_compound_t *fc,
                     vrna_eqm_probs_f       cb,
                     vrna_eqm_probs_data_t  cb_data,
                     vrna_eqm_probs_opt_t   options);


PRIVATE void
free_dp_helper_arrays(helper_arrays_t dp_helpers);


PRIVATE void
rotate_dp_helper_arrays_outer(helper_arrays_t dp_helpers);


PRIVATE void
rotate_dp_helper_arrays_inner(helper_arrays_t dp_helpers);


PRIVATE void
compute_up_ext(vrna_fold_compound_t *fc,
               unsigned int         l,
               helper_arrays_t      dp_helpers);


PRIVATE void
compute_up_hp(vrna_fold_compound_t *fc,
              unsigned int         l,
              helper_arrays_t      dp_helpers);


PRIVATE void
compute_bpp_internal(vrna_fold_compound_t *fc,
                     unsigned int         l,
                     FLT_OR_DBL           *Qmax,
                     int                  *ov,
                     helper_arrays_t      dp_helpers);


PRIVATE void
compute_bpp_internal_comparative(vrna_fold_compound_t *fc,
                                 unsigned int         l,
                                 FLT_OR_DBL           *Qmax,
                                 int                  *ov,
                                 helper_arrays_t      dp_helpers);

PRIVATE void
compute_up_int(vrna_fold_compound_t *fc,
               unsigned int         l,
               helper_arrays_t      dp_helpers);


PRIVATE void
compute_up_mb(vrna_fold_compound_t *fc,
              unsigned int         l,
              helper_arrays_t      dp_helpers);


PRIVATE void
compute_gquad_prob_internal(vrna_fold_compound_t  *fc,
                            unsigned int          l);


PRIVATE void
compute_gquad_prob_internal_comparative(vrna_fold_compound_t  *fc,
                                        unsigned int          l);


PRIVATE void
compute_bpp_multibranch(vrna_fold_compound_t  *fc,
                        unsigned int          l,
                        FLT_OR_DBL            *Qmax,
                        int                   *ov,
                        helper_arrays_t      dp_helpers);


PRIVATE void
compute_bpp_multibranch_comparative(vrna_fold_compound_t  *fc,
                                    unsigned int          l,
                                    FLT_OR_DBL            *Qmax,
                                    int                   *ov,
                                    helper_arrays_t      dp_helpers);


PRIVATE FLT_OR_DBL
contrib_ext_pair(vrna_fold_compound_t *fc,
                 unsigned int         i,
                 unsigned int         j,
                 helper_arrays_t      dp_helpers);


PRIVATE FLT_OR_DBL
contrib_ext_pair_comparative(vrna_fold_compound_t *fc,
                             unsigned int         i,
                             unsigned int         j,
                             helper_arrays_t      dp_helpers);


PRIVATE FLT_OR_DBL
contrib_ext_up(vrna_fold_compound_t *fc,
               unsigned int         i,
               unsigned int         j,
               helper_arrays_t      dp_helpers);


PRIVATE FLT_OR_DBL
contrib_ext_up_comparative(vrna_fold_compound_t *fc,
                           unsigned int         i,
                           unsigned int         j,
                           helper_arrays_t      dp_helpers);


PRIVATE void
multistrand_update_Y5(vrna_fold_compound_t  *fc,
                      unsigned int          l,
                      helper_arrays_t      dp_helpers);


PRIVATE void
multistrand_update_Y3(vrna_fold_compound_t  *fc,
                      unsigned int          l,
                      helper_arrays_t      dp_helpers);


PRIVATE void
multistrand_contrib(vrna_fold_compound_t  *fc,
                    unsigned int          l,
                    helper_arrays_t      dp_helpers,
                    FLT_OR_DBL            *Qmax,
                    int                   *ov);


/*
 #################################
 # BEGIN OF FUNCTION DEFINITIONS #
 #################################
 */
PUBLIC int
vrna_eqm_probabilities(vrna_fold_compound_t   *fc,
                       vrna_eqm_probs_f       cb,
                       vrna_eqm_probs_data_t  cb_data,
                       vrna_eqm_probs_opt_t   options)
{
  if ((fc) &&
      (cb)) {

    if (options.mode == 0) /* default mode */
      options = default_options;

    return pf_create_bppm(fc,
                          cb,
                          cb_data,
                          options);
  }

  return 0;
}


PUBLIC int
vrna_pairing_probs(vrna_fold_compound_t *vc,
                   char                 *structure)
{
  int ret = 0;

  if (vc) {
    ret = pf_create_bppm(vc,
                         NULL,
                         NULL,
                         default_options);

    if (structure != NULL) {
      char *s = vrna_pairing_tendency(vc);
      memcpy(structure, s, vc->length);
      structure[vc->length] = '\0';
      free(s);
    }
  }

  return ret;
}


PUBLIC void
vrna_pf_dimer_probs(double                  FAB,
                    double                  FA,
                    double                  FB,
                    vrna_ep_t               *prAB,
                    const vrna_ep_t         *prA,
                    const vrna_ep_t         *prB,
                    int                     Alength,
                    const vrna_exp_param_t  *exp_params)
{
  /* computes binding probabilities and dimer free energies */
  int             i, j;
  double          pAB;
  double          mykT;
  const vrna_ep_t *lp2;
  vrna_ep_t       *lp1;
  int             offset;

  mykT = exp_params->kT / 1000.;

  /* pair probabilities in pr are relative to the null model (without DuplexInit) */

  /* Compute probabilities pAB, pAA, pBB */

  pAB = 1. - exp((1 / mykT) * (FAB - FA - FB));

  /*
   * compute pair probabilities given that it is a dimer
   * AB dimer
   */
  offset  = 0;
  lp2     = prA;
  if (pAB > 0) {
    for (lp1 = prAB; lp1->j > 0; lp1++) {
      float pp = 0;
      i = lp1->i;
      j = lp1->j;
      while ((offset + lp2->i < i) &&
             (lp2->i > 0))
        lp2++;
      if (offset + lp2->i == i)
        while (((offset + lp2->j) < j) &&
               (lp2->j > 0))
          lp2++;

      if (lp2->j == 0) {
        lp2     = prB;
        offset  = Alength;
      }                                          /* jump to next list */

      if ((offset + lp2->i == i) &&
          (offset + lp2->j == j)) {
        pp = lp2->p;
        lp2++;
      }

      lp1->p = (lp1->p - (1 - pAB) * pp) / pAB;
      if (lp1->p < 0.) {
        vrna_log_warning(
          "vrna_co_pf_probs: numeric instability detected, probability below zero!");
        lp1->p = 0.;
      }
    }
  }

  return;
}

/*
 #################################
 # STATIC helper functions below #
 #################################
 */
PRIVATE void
compute_accessibility(vrna_fold_compound_t  *fc,
                      unsigned int          l,
                      helper_arrays_t      dp_helpers)
{
  /* actual computation for all subsequences of up to size u ending at position l */
  compute_up_ext(fc, l, dp_helpers);
  compute_up_hp(fc, l, dp_helpers);
  compute_up_int(fc, l, dp_helpers);
  compute_up_mb(fc, l, dp_helpers);
}


/*
 *  Compute helix probabilities for all helices ending at l
 *  This function assumes probs[kl] contains the outside probabilities for 
 *  basepair (k,l)
 *  Helix probabilities are returned by callback(fc, pr, l, k, s, VRNA_EQM_PROBS_HX, data)
 *  for all helices of size s ending at l and starting at k, if any.
 */
PRIVATE void
compute_hx_prob(vrna_fold_compound_t  *fc,
                unsigned int          l,
                helper_arrays_t       dp_helpers)
{
  vrna_eqm_probs_f       cb       = dp_helpers->cb;
  vrna_eqm_probs_data_t  cb_data  = dp_helpers->cb_data;

  unsigned int      k, n, s, start, end;
  int               *my_iindx, kl;
  FLT_OR_DBL        *qb, *probs, p, *hxp, qbt1;
  vrna_mx_pf_t      *matrices;

  n         = fc->length;
  my_iindx  = fc->iindx;
  matrices  = fc->exp_matrices;
  qb        = matrices->qb;
  probs     = matrices->probs;
  hxp       = dp_helpers->hxp;

  for (k = 1; k < l; k++) {
    /* reset probabilities */
    hxp = memset(hxp, 0, sizeof(FLT_OR_DBL) * (dp_helpers->options.hx_length + 1));

    kl = my_iindx[k] - l;
    s  = 1;

    /* get initial probs, i.e. base pair probabilities */
    p = probs[kl];

    hxp[s++] = p *
               qb[kl];

    if (p > FLT_MIN) {
      /* loop over all helix sizes */
      for (s = 2; s <= dp_helpers->options.hx_length; s++) {
        start = k + s - 2;
        end   = l + 2 - s;

        if (start + 2 >= end)
          break;

        qbt1 = qb[my_iindx[start + 1] - (end - 1)];

        if (qbt1 < FLT_MIN)
          break;

        p *= vrna_exp_eval_stack(fc, start, end, VRNA_EVAL_LOOP_DEFAULT);

        hxp[s] = qbt1 *
                 p;
      }
    }

    cb(fc, hxp, l, k, s - 1, VRNA_EQM_PROBS_HX, cb_data);
  }
}


PRIVATE void
init_probs_dp_mx(vrna_fold_compound_t *fc,
                 void                 *data)
{
  unsigned int  i, j, n;
  int           *my_iindx;
  FLT_OR_DBL    *probs;

  n         = fc->length;
  my_iindx  = fc->iindx;
  probs     = fc->exp_matrices->probs;

  /* init probability matrix */
  for (i = 1; i <= n; i++)
    for (j = i; j <= n; j++)
      probs[my_iindx[i] - j] = 0.;
}


PRIVATE void
init_bpp(vrna_fold_compound_t *fc,
         void                 *data)
{
  unsigned int    n = fc->length;
  helper_arrays_t d = (helper_arrays_t)data;

  /* reserve memory for probability output */
  d->probs_ext = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * (n + 1));
  d->probs_hp  = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * (n + 1));
  d->probs_int = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * (n + 1));
  d->probs_mb  = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * (n + 1));
}


PRIVATE void
report_n_clear_bpp(vrna_fold_compound_t *fc,
                   unsigned int         l,
                   void                 *data)
{
  unsigned int          n       = fc->length;
  helper_arrays_t       d       = (helper_arrays_t)data;
  vrna_eqm_probs_f      cb      = d->cb;
  vrna_eqm_probs_data_t cb_data = d->cb_data;

  /* loop-type specific output */
  if ((d->cb) &&
      (d->options.mode & VRNA_EQM_PROBS_BPP)) {
    if (d->options.mode & VRNA_EQM_PROBS_SPLIT) {
      cb(fc, d->probs_ext, l, 1, l - 1, VRNA_EQM_PROBS_BPP | VRNA_LOOP_TYPE_EXT, cb_data);
      cb(fc, d->probs_int, l, 1, l - 1, VRNA_EQM_PROBS_BPP | VRNA_LOOP_TYPE_INT, cb_data);
      cb(fc, d->probs_mb, l, 1, l - 1, VRNA_EQM_PROBS_BPP | VRNA_LOOP_TYPE_MB, cb_data);
    }

    /* combined probabilities */
    for (unsigned int k = 1; k <= l - 1; k++)
      d->probs_ext[k] += d->probs_int[k] +
                         d->probs_mb[k];

    cb(fc, d->probs_ext, l, 1, l - 1, VRNA_EQM_PROBS_BPP | VRNA_LOOP_TYPE_ANY, cb_data);
  }

  memset(d->probs_ext, 0, sizeof(FLT_OR_DBL) * (n + 1));
  memset(d->probs_hp, 0, sizeof(FLT_OR_DBL) * (n + 1));
  memset(d->probs_int, 0, sizeof(FLT_OR_DBL) * (n + 1));
  memset(d->probs_mb, 0, sizeof(FLT_OR_DBL) * (n + 1));
}


PRIVATE void
init_accessibility(vrna_fold_compound_t *fc,
                   void                 *data)
{
  unsigned int    n = fc->length;
  helper_arrays_t d = (helper_arrays_t)data;

  if (d->options.up_length == 0)
    d->options.up_length = n;

  d->qbh = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * (n + 2));

  /* reserve memory for accessibility output */
  d->access_ext  = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * (d->options.up_length + 1));
  d->access_hp   = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * (d->options.up_length + 1));
  d->access_int  = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * (d->options.up_length + 1));
  d->access_mb   = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * (d->options.up_length + 1));

  /* recompute qm2 if necessary */
  if (fc->exp_matrices->qm2_real == NULL) {
    short int             *S, *S1;
    unsigned int          turn, type, i, j, u;
    int                   *my_iindx, ii, ij;
    FLT_OR_DBL            *qm, *qb, temp, *expMLbase;
    vrna_exp_param_t      *params;
    vrna_md_t             *md;
    vrna_hc_t             *hc;
    struct sc_mb_exp_dat  *sc_wrapper;

    S           = fc->sequence_encoding2;
    S1          = fc->sequence_encoding;
    hc          = fc->hc;
    qm          = fc->exp_matrices->qm;
    qb          = fc->exp_matrices->qb;
    expMLbase   = fc->exp_matrices->expMLbase;
    my_iindx    = fc->iindx;
    params      = fc->exp_params;
    md          = &(params->model_details);
    turn        = md->min_loop_size;
    my_iindx    = fc->iindx;
    sc_wrapper  = &(d->sc_wrapper_mb);

    d->qm2_aux = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * ((n + 1) * (n + 2) / 2 + 1));

    for (i = 1; i < n; i++) {
      ii = my_iindx[i];
      for (j = i + 1; j <= n; j++) {
        ij = ii - j;

        /* j is unpaired */
        if (hc->eval_mb(i, j, i, j - 1, VRNA_DECOMP_ML_ML, hc)) {
          temp = d->qm2_aux[ij + 1] *
                 expMLbase[1];

          if (sc_wrapper->red_ml)
            temp *= sc_wrapper->red_ml(i, j, i, j - 1, sc_wrapper);

          d->qm2_aux[ij] = temp;
        }

        /* j pairs with some u */
        for (u = i + 1; u + turn < j; u++) {
          if (hc->eval_mb(i, j, u - 1, u, VRNA_DECOMP_ML_ML_ML, hc)) {
            type = vrna_get_ptype_md(S[u], S[j], md);

            temp = qm[my_iindx[i] - (u - 1)] *
                   qb[my_iindx[u] - j] *
                   vrna_exp_E_multibranch_stem(type,
                                               S1[u - 1],
                                               S1[j + 1],
                                               params);

            if (sc_wrapper->decomp_ml)
              temp *= sc_wrapper->decomp_ml(i, j, u - 1, u, sc_wrapper);

            d->qm2_aux[ij] += temp;
          }
        }
      }
    }
  }
}


PRIVATE void
report_n_clear_accessibility(vrna_fold_compound_t *fc,
                             unsigned int         l,
                             void                 *data)
{
  helper_arrays_t       d       = (helper_arrays_t)data;
  vrna_eqm_probs_f      cb      = d->cb;
  vrna_eqm_probs_data_t cb_data = d->cb_data;
  unsigned int          end     = (l < d->options.up_length) ? l : d->options.up_length;

  /* loop-type specific output */
  if ((d->cb) &&
      (d->options.mode & VRNA_EQM_PROBS_UP)) {
    if (d->options.mode & VRNA_EQM_PROBS_SPLIT) {
      cb(fc, d->access_ext, l, 1, end, VRNA_EQM_PROBS_UP | VRNA_LOOP_TYPE_EXT, cb_data);
      cb(fc, d->access_hp, l, 1, end, VRNA_EQM_PROBS_UP | VRNA_LOOP_TYPE_HP, cb_data);
      cb(fc, d->access_int, l, 1, end, VRNA_EQM_PROBS_UP | VRNA_LOOP_TYPE_INT, cb_data);
      cb(fc, d->access_mb, l, 1, end, VRNA_EQM_PROBS_UP | VRNA_LOOP_TYPE_MB, cb_data);
    }

    /* combined probabilities */
    for (unsigned int u = 1; u <= end; u++)
      d->access_ext[u] += d->access_hp[u] +
                                   d->access_int[u] +
                                   d->access_mb[u];

    cb(fc, d->access_ext, l, 1, end, VRNA_EQM_PROBS_UP | VRNA_LOOP_TYPE_ANY, cb_data);
  }

  /* reset data */
  memset(d->access_ext, 0, sizeof(FLT_OR_DBL) * (d->options.up_length + 1));
  memset(d->access_hp, 0, sizeof(FLT_OR_DBL) * (d->options.up_length + 1));
  memset(d->access_int, 0, sizeof(FLT_OR_DBL) * (d->options.up_length + 1));
  memset(d->access_mb, 0, sizeof(FLT_OR_DBL) * (d->options.up_length + 1));
}


PRIVATE void
init_hxp(vrna_fold_compound_t *fc,
         void                 *data)
{
  unsigned int    n = fc->length;
  helper_arrays_t d = (helper_arrays_t)data;

  if (d->options.hx_length == 0)
    d->options.hx_length = (n + 1) / 2;

  /* reserve memory for probability output */
  d->hxp = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * (d->options.hx_length + 1));
}


PRIVATE void
init_multistrand(vrna_fold_compound_t *fc,
                 void                 *data)
{
  unsigned int    strands = fc->strands;
  unsigned int    n       = fc->length;
  helper_arrays_t d       = (helper_arrays_t)data;

  /* reserve memory for probability output */
  d->Y5  = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * strands);
  d->Y5p = (FLT_OR_DBL **)vrna_alloc(sizeof(FLT_OR_DBL *) * strands);
  d->Y3  = (FLT_OR_DBL **)vrna_alloc(sizeof(FLT_OR_DBL *) * strands);
  d->Y3p = (FLT_OR_DBL **)vrna_alloc(sizeof(FLT_OR_DBL *) * strands);

  for (unsigned int s = 0; s < strands; s++) {
    d->Y5p[s] = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * (n + 1));
    d->Y3[s]  = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * (n + 1));
    d->Y3p[s] = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * (n + 1));
  }
}


PRIVATE void
update_multistrand(vrna_fold_compound_t *fc,
                   unsigned int         l,
                   void                 *data)
{
  helper_arrays_t d = (helper_arrays_t)data;

  if (l > 1) {
    multistrand_update_Y5(fc, l - 1, d);
    multistrand_update_Y3(fc, l - 1, d);
  }
}


PRIVATE void
register_outside_component(vrna_fold_compound_t *fc,
                           helper_arrays_t      dp_helpers,
                           component_f          component_cb,
                           void *               component_cb_data,
                           vrna_auxdata_free_f  component_cb_data_free,
                           pre_f                component_cb_init,
                           post_f               component_cb_post)
{
  if (component_cb) {
    vrna_array_append(dp_helpers->component_cbs, component_cb);
    vrna_array_append(dp_helpers->component_cbs_data, component_cb_data);
    vrna_array_append(dp_helpers->component_cbs_data_free, component_cb_data_free);

    if (component_cb_init)
      vrna_array_append(dp_helpers->component_cbs_init, component_cb_init);

    if (component_cb_post)
      vrna_array_append(dp_helpers->component_cbs_post, component_cb_post);
  }
}



/* compute outside partition function and corresponding probabilities */
PRIVATE int
pf_create_bppm(vrna_fold_compound_t   *vc,
               vrna_eqm_probs_f       cb,
               vrna_eqm_probs_data_t  cb_data,
               vrna_eqm_probs_opt_t   options)
{
  int               n, i, j, l, ij, *pscore, *jindx, ov = 0;
  FLT_OR_DBL        Qmax = 0;
  FLT_OR_DBL        *qb, *probs, q_g;
  FLT_OR_DBL        *q1k, *qln;

  int               with_gquad;
  vrna_hc_t         *hc;
  vrna_sc_t         *sc;
  int               *my_iindx;
  int               circular, with_ud, with_ud_outside;
  vrna_exp_param_t  *pf_params;
  vrna_mx_pf_t      *matrices;
  vrna_md_t         *md;
  vrna_ud_t         *domains_up;
  vrna_smx_csr(FLT_OR_DBL) *q_gq;
  helper_arrays_t   dp_helpers;

  n           = vc->length;
  pscore      = (vc->type == VRNA_FC_TYPE_COMPARATIVE) ? vc->pscore : NULL;
  pf_params   = vc->exp_params;
  md          = &(pf_params->model_details);
  my_iindx    = vc->iindx;
  jindx       = vc->jindx;
  circular    = md->circ;
  with_gquad  = md->gquad;

  hc  = vc->hc;
  sc  = vc->sc;

  domains_up  = vc->domains_up;
  matrices    = vc->exp_matrices;

  qb    = matrices->qb;
  q_gq  = matrices->q_gq;
  probs = matrices->probs;
  q1k   = matrices->q1k;
  qln   = matrices->qln;

  with_ud         = (domains_up && domains_up->exp_energy_cb) ? 1 : 0;
  with_ud_outside = (with_ud && domains_up->probs_add) ? 1 : 0;

  /*
   * the following is a crude check whether the partition function forward recursion
   * has already been taken place
   */
  if ((qb) &&
      (probs) &&
      (circular ? matrices->qm2_real != NULL : (q1k != NULL && qln != NULL))) {
    with_gquad = pf_params->model_details.gquad;

    double      kTn = pf_params->kT / 10.;               /* kT in cal/mol  */


    dp_helpers  = get_dp_helper_arrays(vc, cb, cb_data, options);

    vrna_array(vrna_ep_t) bp_correction = dp_helpers->bp_correction;

    Qmax = 0;

    /* call registered init callbacks */
    for (size_t cnt = 0; cnt < vrna_array_size(dp_helpers->init_cbs); cnt++)
      dp_helpers->init_cbs[cnt](vc, dp_helpers->init_cbs_data[cnt]);

    /* Begin computing outside probabilities */

    for (l = n; l > 0; l--) {
      /* 1. external loop pairs, i.e. pairs not enclosed by any other pair (or external loop for circular RNAs) */
      dp_helpers->run_outside_external(vc, l, &Qmax, &ov, dp_helpers);

      /* 2. all cases where base pair (k,l) is enclosed by another pair (i,j) */
      outside_bp_enclosed(vc, l, &Qmax, &ov, dp_helpers);

      if (options.mode & VRNA_EQM_PROBS_UP)
        compute_accessibility(vc, l, dp_helpers);

      if (options.mode & VRNA_EQM_PROBS_HX)
        compute_hx_prob(vc, l, dp_helpers);

      for (size_t cnt = 0; cnt < vrna_array_size(dp_helpers->component_cbs); cnt++)
        dp_helpers->component_cbs[cnt](vc,
                                       l,
                                       &Qmax,
                                       &ov,
                                       dp_helpers->component_cbs_data[cnt]);

      /* call registered post callbacks */
      for (size_t cnt = 0; cnt < vrna_array_size(dp_helpers->post_cbs); cnt++)
        dp_helpers->post_cbs[cnt](vc, l, dp_helpers->post_cbs_data[cnt]);
    }

    if (vc->type == VRNA_FC_TYPE_SINGLE) {
      if (with_ud_outside) {
        /*
         *  The above recursions only deal with base pairs, and how they might be
         *  enclosed by other pairs. However, for unstructrued domains, we have
         *  unpaired stretches, and require information about how these are enclosed
         *  by base pairs.
         */

        /* 1. Exterior loops */
        ud_outside_ext_loops(vc);

        /* 2. Hairpin loops */
        ud_outside_hp_loops(vc);

        /* 3. Interior loops */
        ud_outside_int_loops(vc);

        /* 4. Multi branch loops */
        ud_outside_mb_loops(vc);
      }

      if ((sc) &&
          (sc->f) &&
          (sc->bt)) {
        for (i = 1; i <= n; i++)
          for (j = i + 1; j <= n; j++) {
            ij = my_iindx[i] - j;
            /*  search for possible auxiliary base pairs in hairpin loop motifs to store
             *  the corresponding probability corrections
             */
            if (hc->mx[i * n + j] & VRNA_CONSTRAINT_CONTEXT_HP_LOOP) {
              vrna_basepair_t *ptr, *aux_bps;
              aux_bps = sc->bt(i, j, i, j, VRNA_DECOMP_PAIR_HP, sc->data);
              if (aux_bps) {
                FLT_OR_DBL qhp = vrna_exp_eval_hairpin(vc, i, j, VRNA_EVAL_LOOP_DEFAULT);
                for (ptr = aux_bps; (ptr) && (ptr->i != 0); ptr++) {
                  vrna_array_append(bp_correction,
                                    ((vrna_ep_t){
                                      .i = ptr->i,
                                      .j = ptr->j,
                                      .p = probs[ij] * qhp
                                    }));
                }
              }

              free(aux_bps);
            }
          }

        /*  correct pairing probabilities for auxiliary base pairs from hairpin-, or internal loop motifs
         *  as augmented by the generalized soft constraints feature
         */
        for (i = 0; i < vrna_array_size(bp_correction); i++) {
          ij = my_iindx[bp_correction[i].i] - bp_correction[i].j;
          /* printf("correcting pair %d, %d by %f\n", bp_correction[i].i, bp_correction[i].j, bp_correction[i].p); */
          probs[ij] += bp_correction[i].p / qb[ij];
        }
      }
    }

    /* done computing the outside partition function */

    /* add inside contributions to transform probs[ij] into base pair probabilities */
    if (!(options.mode & VRNA_EQM_PROBS_STORE_OUTSIDE)) {
      for (i = 1; i <= n; i++) {
        for (j = i + 1; j <= n; j++) {
          ij = my_iindx[i] - j;

          if (with_gquad) {
            if (qb[ij] > 0.) {
              probs[ij] *= qb[ij];
              if (vc->type == VRNA_FC_TYPE_COMPARATIVE)
                probs[ij] *= exp(-pscore[jindx[j] + i] / kTn);
            } else {
#ifndef VRNA_DISABLE_C11_FEATURES
              q_g = vrna_smx_csr_get(q_gq, i, j, 0.);
#else
              q_g = vrna_smx_csr_FLT_OR_DBL_get(q_gq, i, j, 0.);
#endif
              if (q_g != 0.) {
                if (!md->circ)
                  probs[ij] += q1k[i - 1] *
                               qln[j + 1] / q1k[n];
  
                probs[ij] *= q_g;
              }
            }
          } else {
            if (qb[ij] > 0.) {
              probs[ij] *= qb[ij];

              if (vc->type == VRNA_FC_TYPE_COMPARATIVE)
                probs[ij] *= exp(-pscore[jindx[j] + i] / kTn);
            }
          }
        }
      }
    }


    if (ov > 0)
      vrna_log_warning("%d overflows occurred while backtracking;\n"
                           "you might try a smaller pf_scale than %g\n",
                           ov, pf_params->pf_scale);

    /* clean up */
    free_dp_helper_arrays(dp_helpers);
  } /* end if 'check for forward recursion' */
  else {
    vrna_log_warning("bppm calculations have to be done after calling forward recursion");
    return 0;
  }

  return 1;
}


PRIVATE void
outside_bp_enclosed(vrna_fold_compound_t  *fc,
                    unsigned int          l,
                    FLT_OR_DBL            *Qmax,
                    int                   *ov,
                    helper_arrays_t      dp_helpers)
{
  vrna_exp_param_t  *pf_params = fc->exp_params;

  /* actual computation for all base pairs (k,l) enclosed by another base pair (i,j) */

  /* 2. all cases where base pair (k,l) is enclosed by another pair (i,j) */
  dp_helpers->run_outside_internal(fc,
                                   l,
                                   Qmax,
                                   ov,
                                   dp_helpers);

  dp_helpers->run_outside_multibranch(fc,
                                      l,
                                      Qmax,
                                      ov,
                                      dp_helpers);

  if (fc->strands > 1)
    multistrand_contrib(fc,
                        l,
                        dp_helpers,
                        Qmax,
                        ov);
}


PRIVATE helper_arrays_t
get_dp_helper_arrays(vrna_fold_compound_t   *fc,
                     vrna_eqm_probs_f       cb,
                     vrna_eqm_probs_data_t  cb_data,
                     vrna_eqm_probs_opt_t   options)
{
  unsigned int  n, i, j, u, s;
  int           with_ud, *my_iindx, ii, ij;
  vrna_ud_t     *domains_up;
  struct helper_arrays *dp_helpers;

  n           = fc->length;
  domains_up  = fc->domains_up;
  with_ud     = (domains_up && domains_up->exp_energy_cb) ? 1 : 0;

  dp_helpers = (struct helper_arrays *)vrna_alloc(sizeof(struct helper_arrays));

  dp_helpers->n       = n;
  dp_helpers->strands = fc->strands;
  dp_helpers->prm_l   = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * (n + 2));
  dp_helpers->prm_l1  = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * (n + 2));
  dp_helpers->prml    = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * (n + 2));

  dp_helpers->probs_ext = NULL;
  dp_helpers->probs_hp  = NULL;
  dp_helpers->probs_int = NULL;
  dp_helpers->probs_mb  = NULL;

  dp_helpers->Y5  = NULL;
  dp_helpers->Y5p = NULL;
  dp_helpers->Y3  = NULL;
  dp_helpers->Y3p = NULL;

  dp_helpers->ud_max_size = 0;
  dp_helpers->pmlu        = NULL;
  dp_helpers->prm_MLbu    = NULL;

  dp_helpers->qm2_aux     = NULL;
  dp_helpers->qbh         = NULL;
  dp_helpers->access_ext  = NULL;
  dp_helpers->access_hp   = NULL;
  dp_helpers->access_int  = NULL;
  dp_helpers->access_mb   = NULL;

  dp_helpers->hxp         = NULL;

  dp_helpers->cb = cb;
  dp_helpers->cb_data = cb_data;
  dp_helpers->options = options;

  dp_helpers->bp_correction = NULL;

  vrna_array_init(dp_helpers->component_cbs);
  vrna_array_init(dp_helpers->component_cbs_data);
  vrna_array_init(dp_helpers->component_cbs_data_free);
  vrna_array_init(dp_helpers->component_cbs_init);
  vrna_array_init(dp_helpers->component_cbs_post);

  vrna_array_init(dp_helpers->init_cbs);
  vrna_array_init(dp_helpers->init_cbs_data);
  vrna_array_init(dp_helpers->post_cbs);
  vrna_array_init(dp_helpers->post_cbs_data);

  init_sc_ext_exp(fc, &(dp_helpers->sc_wrapper_ext));
  init_sc_hp_exp(fc, &(dp_helpers->sc_wrapper_hp));
  init_sc_int_exp(fc, &(dp_helpers->sc_wrapper_int));
  init_sc_mb_exp(fc, &(dp_helpers->sc_wrapper_mb));

  /* Register callbacks */

  /* always initialize the probs DP matrix */
  vrna_array_append(dp_helpers->init_cbs, init_probs_dp_mx);
  vrna_array_append(dp_helpers->init_cbs_data, NULL);

  if (options.mode & VRNA_EQM_PROBS_BPP) {
    vrna_array_append(dp_helpers->init_cbs, init_bpp);
    vrna_array_append(dp_helpers->init_cbs_data, (void *)dp_helpers);
    vrna_array_append(dp_helpers->post_cbs, report_n_clear_bpp);
    vrna_array_append(dp_helpers->post_cbs_data, (void *)dp_helpers);
  }

  if (options.mode & VRNA_EQM_PROBS_UP) {
    vrna_array_append(dp_helpers->init_cbs, init_accessibility);
    vrna_array_append(dp_helpers->init_cbs_data, (void *)dp_helpers);
    vrna_array_append(dp_helpers->post_cbs, report_n_clear_accessibility);
    vrna_array_append(dp_helpers->post_cbs_data, (void *)dp_helpers);
  }

  if (options.mode & VRNA_EQM_PROBS_HX) {
    vrna_array_append(dp_helpers->init_cbs, init_hxp);
    vrna_array_append(dp_helpers->init_cbs_data, (void *)dp_helpers);
  }

  if (fc->strands > 1) {
    vrna_array_append(dp_helpers->init_cbs, init_multistrand);
    vrna_array_append(dp_helpers->init_cbs_data, (void *)dp_helpers);
    vrna_array_append(dp_helpers->post_cbs, update_multistrand);
    vrna_array_append(dp_helpers->post_cbs_data, (void *)dp_helpers);
  }

  if ((fc->type == VRNA_FC_TYPE_SINGLE) &&
      (fc->sc) &&
      (fc->sc->bt))
    vrna_array_init(dp_helpers->bp_correction);


  if (fc->exp_params->model_details.gquad) {
    fc->exp_matrices->p_gq = vrna_smx_csr_FLT_OR_DBL_init(n + 1);
  }


  if (with_ud) {
    /* find out maximum size of any unstructured domain */
    for (u = 0; u < (unsigned int)domains_up->uniq_motif_count; u++)
      if (dp_helpers->ud_max_size < domains_up->uniq_motif_size[u])
        dp_helpers->ud_max_size = domains_up->uniq_motif_size[u];

    dp_helpers->pmlu = (FLT_OR_DBL **)vrna_alloc(
      sizeof(FLT_OR_DBL *) * (dp_helpers->ud_max_size + 1)); /* maximum motif size */

    for (u = 0; u <= dp_helpers->ud_max_size; u++)
      dp_helpers->pmlu[u] = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * (n + 2));

    dp_helpers->prm_MLbu =
      (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * (dp_helpers->ud_max_size + 1));

    for (u = 0; u <= dp_helpers->ud_max_size; u++)
      dp_helpers->prm_MLbu[u] = 0.;
  }

  if (fc->exp_params->model_details.circ)
    dp_helpers->run_outside_external = &outside_bp_circ;
  else
    dp_helpers->run_outside_external = &outside_bp_external;

  if (fc->type == VRNA_FC_TYPE_SINGLE) {
    dp_helpers->run_outside_internal    = &compute_bpp_internal;
    dp_helpers->run_outside_multibranch = &compute_bpp_multibranch;
  } else {
    dp_helpers->run_outside_internal    = &compute_bpp_internal_comparative;
    dp_helpers->run_outside_multibranch = &compute_bpp_multibranch_comparative;
  }

  return dp_helpers;
}


PRIVATE void
free_dp_helper_arrays(helper_arrays_t dp_helpers)
{
  unsigned int u;

  free(dp_helpers->prm_l);
  free(dp_helpers->prm_l1);
  free(dp_helpers->prml);

  free(dp_helpers->Y5);

  free_sc_ext_exp(&(dp_helpers->sc_wrapper_ext));
  free_sc_hp_exp(&(dp_helpers->sc_wrapper_hp));
  free_sc_int_exp(&(dp_helpers->sc_wrapper_int));
  free_sc_mb_exp(&(dp_helpers->sc_wrapper_mb));

  vrna_array_free(dp_helpers->bp_correction);

  /* release component callback related memory */
  for (size_t cnt = 0; cnt < vrna_array_size(dp_helpers->component_cbs_data_free); cnt++)
    if (dp_helpers->component_cbs_data_free[cnt])
      dp_helpers->component_cbs_data_free[cnt](dp_helpers->component_cbs_data[cnt]);

  vrna_array_free(dp_helpers->component_cbs);
  vrna_array_free(dp_helpers->component_cbs_data);
  vrna_array_free(dp_helpers->component_cbs_data_free);
  vrna_array_free(dp_helpers->component_cbs_init);
  vrna_array_free(dp_helpers->component_cbs_post);

  vrna_array_free(dp_helpers->init_cbs);
  vrna_array_free(dp_helpers->init_cbs_data);
  vrna_array_free(dp_helpers->post_cbs);
  vrna_array_free(dp_helpers->post_cbs_data);

  free(dp_helpers->probs_ext);
  free(dp_helpers->probs_hp);
  free(dp_helpers->probs_int);
  free(dp_helpers->probs_mb);

  if (dp_helpers->Y5p)
    for (unsigned int s = 0; s < dp_helpers->strands; s++)
      free(dp_helpers->Y5p[s]);

  free(dp_helpers->Y5p);

  if (dp_helpers->Y3)
    for (unsigned int s = 0; s < dp_helpers->strands; s++)
      free(dp_helpers->Y3[s]);

  free(dp_helpers->Y3);

  if (dp_helpers->Y3p)
    for (unsigned int s = 0; s < dp_helpers->strands; s++)
      free(dp_helpers->Y3p[s]);

  free(dp_helpers->Y3p);

  if (dp_helpers->pmlu) {
    for (u = 0; u <= dp_helpers->ud_max_size; u++)
      free(dp_helpers->pmlu[u]);
    free(dp_helpers->pmlu);
  }

  free(dp_helpers->qbh);
  free(dp_helpers->qm2_aux);

  free(dp_helpers->access_ext);
  free(dp_helpers->access_hp);
  free(dp_helpers->access_int);
  free(dp_helpers->access_mb);

  free(dp_helpers->hxp);

  free(dp_helpers->prm_MLbu);
  free(dp_helpers);
}


PRIVATE void
rotate_dp_helper_arrays_outer(helper_arrays_t dp_helpers)
{
  unsigned int  u;
  FLT_OR_DBL    *tmp;

  /* rotate prm_l and prm_l1 arrays */
  tmp                 = dp_helpers->prm_l1;
  dp_helpers->prm_l1  = dp_helpers->prm_l;
  dp_helpers->prm_l   = tmp;

  /* rotate pmlu entries required for unstructured domain feature */
  if (dp_helpers->pmlu) {
    tmp = dp_helpers->pmlu[dp_helpers->ud_max_size];

    for (u = dp_helpers->ud_max_size; u > 0; u--)
      dp_helpers->pmlu[u] = dp_helpers->pmlu[u - 1];

    dp_helpers->pmlu[0] = tmp;

    for (u = 0; u <= dp_helpers->ud_max_size; u++)
      dp_helpers->prm_MLbu[u] = 0.;
  }
}


PRIVATE void
rotate_dp_helper_arrays_inner(helper_arrays_t dp_helpers)
{
  unsigned int u;

  /* rotate prm_MLbu entries required for unstructured domain feature */
  if (dp_helpers->prm_MLbu)
    for (u = dp_helpers->ud_max_size; u > 0; u--)
      dp_helpers->prm_MLbu[u] = dp_helpers->prm_MLbu[u - 1];
}


PRIVATE void
outside_bp_external(vrna_fold_compound_t  *fc,
                     unsigned int          l,
                     FLT_OR_DBL            *Qmax,
                     int                   *ov,
                     helper_arrays_t       dp_helpers)
{
  unsigned int              k, n;
  int                       *my_iindx, kl;
  FLT_OR_DBL                *probs, *q1k, *qln, *qb, *probs_ext, q;
  vrna_mx_pf_t              *matrices;
  vrna_hc_t                 *hc;

  FLT_OR_DBL                (*contrib_f)(vrna_fold_compound_t *,
                                         unsigned int,
                                         unsigned int,
                                         helper_arrays_t dp_helpers);

  n         = fc->length;
  my_iindx  = fc->iindx;
  matrices  = fc->exp_matrices;
  qb        = matrices->qb;
  probs     = matrices->probs;
  q1k       = matrices->q1k;
  qln       = matrices->qln;
  hc        = fc->hc;
  probs_ext = dp_helpers->probs_ext;

  contrib_f = (fc->type == VRNA_FC_TYPE_SINGLE) ? &contrib_ext_pair : &contrib_ext_pair_comparative;

  for (k = 1; k < l; k++) {
    kl = my_iindx[k] - l;
    probs[kl] = q = 0.;

    if ((hc->eval_ext(1, n, k, l, VRNA_DECOMP_EXT_STEM_OUTSIDE, hc)) &&
        (qb[kl] > 0.)) {
      q = q1k[k - 1] *
          qln[l + 1] /
          q1k[n];
      q *= contrib_f(fc, k, l, dp_helpers);

      probs_ext[k] = q *
                     qb[kl];

      probs[kl] = q;
    }
  }
}


PRIVATE void
compute_up_ext(vrna_fold_compound_t *fc,
               unsigned int         l,
               helper_arrays_t      dp_helpers)
{
  unsigned int          i, u, max_u, n;
  FLT_OR_DBL            *q1k, *qln, temp, *acc_ext;
  vrna_mx_pf_t          *matrices;
  vrna_hc_t             *hc;
  vrna_eqm_probs_f      cb;
  vrna_eqm_probs_data_t cb_data;
  FLT_OR_DBL            (*contrib_up)(vrna_fold_compound_t *,
                                      unsigned int,
                                      unsigned int,
                                      helper_arrays_t);

  n         = fc->length;
  matrices  = fc->exp_matrices;
  q1k       = matrices->q1k;
  qln       = matrices->qln;
  hc        = fc->hc;
  max_u     = dp_helpers->options.up_length;
  acc_ext   = dp_helpers->access_ext;
  cb        = dp_helpers->cb;
  cb_data   = dp_helpers->cb_data;

  contrib_up = (fc->type == VRNA_FC_TYPE_SINGLE) ? &contrib_ext_up : &contrib_ext_up_comparative;

  /* compute unpaired probabilities for stretches of size u ending at position l */
  if (l < n) {
    for (u = 1; u <= max_u; u++) {
      i = l + 1 - u;
      if (i > 1) {
        if (hc->eval_ext(1, l, 1, i - 1, VRNA_DECOMP_EXT_EXT, hc)) {
          temp = contrib_up(fc, i, l, dp_helpers) *
                 q1k[i - 1] *
                 qln[l + 1] /
                 q1k[n];

          acc_ext[u] += temp;
        }
      } else {
        if (hc->eval_ext(1, l, 1, l, VRNA_DECOMP_EXT_UP, hc)) {
          temp = contrib_up(fc, 1, l, dp_helpers) *
                 qln[l + 1] /
                 q1k[n];

          acc_ext[u] += temp;
        }
        break;
      }
    }
  } else {
    for (u = 1; u <= max_u; u++) {
      i = l + 1 - u;
      if (i > 1) {
        if (hc->eval_ext(1, l, 1, i - 1, VRNA_DECOMP_EXT_EXT, hc)) {
          temp = contrib_up(fc, i, l, dp_helpers) *
                 q1k[i - 1] /
                 q1k[n];

          acc_ext[u] += temp;
        }
      } else {
        if (hc->eval_ext(1, l, 1, l, VRNA_DECOMP_EXT_UP, hc)) {
          temp = contrib_up(fc, 1, l, dp_helpers) /
                 q1k[n];

          acc_ext[u] += temp;
        }
        break;
      }
    }
  }
}


PRIVATE FLT_OR_DBL
contrib_ext_pair(vrna_fold_compound_t *fc,
                 unsigned int         i,
                 unsigned int         j,
                 helper_arrays_t      dp_helpers)
{
  unsigned char     type;
  char              *ptype;
  short             *S1, s5, s3;
  unsigned int      *sn, n;
  int               *jindx;
  FLT_OR_DBL        contribution;
  vrna_exp_param_t  *pf_params;
  vrna_sc_t         *sc;

  n         = fc->length;
  pf_params = fc->exp_params;
  S1        = fc->sequence_encoding;
  sn        = fc->strand_number;
  ptype     = fc->ptype;
  jindx     = fc->jindx;
  sc        = fc->sc;

  type  = vrna_get_ptype(jindx[j] + i, ptype);
  s5    = ((i > 1) && (sn[i] == sn[i - 1])) ? S1[i - 1] : -1;
  s3    = ((j < n) && (sn[j + 1] == sn[j])) ? S1[j + 1] : -1;

  contribution = vrna_exp_E_exterior_stem(type, s5, s3, pf_params);

  if ((sc) &&
      (sc->exp_f))
    contribution *= sc->exp_f(1, n, i, j, VRNA_DECOMP_EXT_STEM_OUTSIDE, sc->data);

  return contribution;
}


PRIVATE FLT_OR_DBL
contrib_ext_pair_comparative(vrna_fold_compound_t *fc,
                             unsigned int         i,
                             unsigned int         j,
                             helper_arrays_t      dp_helpers)
{
  unsigned char     type;
  short             **S, **S5, **S3, s5, s3;
  unsigned int      **a2s, n, s, n_seq;
  int               *jindx, *pscore;
  FLT_OR_DBL        contribution;
  double            kTn;
  vrna_exp_param_t  *pf_params;
  vrna_md_t         *md;
  vrna_sc_t         **scs;

  n         = fc->length;
  n_seq     = fc->n_seq;
  jindx     = fc->jindx;
  pf_params = fc->exp_params;
  md        = &(pf_params->model_details);
  S         = fc->S;
  S5        = fc->S5;
  S3        = fc->S3;
  a2s       = fc->a2s;
  pscore    = fc->pscore;           /* precomputed array of pair types */
  scs       = fc->scs;
  kTn       = pf_params->kT / 10.;  /* kT in cal/mol  */

  contribution = exp(pscore[jindx[j] + i] / kTn);

  for (s = 0; s < n_seq; s++) {
    type  = vrna_get_ptype_md(S[s][i], S[s][j], md);
    s5    = (a2s[s][i] > 1) ? S5[s][i] : -1;
    s3    = (a2s[s][j] < a2s[s][n]) ? S3[s][j] : -1;

    contribution *= vrna_exp_E_exterior_stem(type, s5, s3, pf_params);
  }

  if (scs) {
    for (s = 0; s < n_seq; s++)
      if (scs[s]->exp_f)
        contribution *= scs[s]->exp_f(1, n, i, j, VRNA_DECOMP_EXT_STEM_OUTSIDE, scs[s]->data);
  }

  return contribution;
}


PRIVATE FLT_OR_DBL
contrib_ext_up(vrna_fold_compound_t *fc,
                 unsigned int         i,
                 unsigned int         j,
                 helper_arrays_t      dp_helpers)
{
  unsigned int      n;
  FLT_OR_DBL        contribution, *scale;
  vrna_sc_t         *sc;

  n             = fc->length;
  scale         = fc->exp_matrices->scale;
  sc            = fc->sc;
  contribution  = scale[j - i + 1];

  if ((sc) &&
      (sc->exp_f)) {
    if (i > 1) {
      if (j < n) {
        contribution *= sc->exp_f(1, n, i - 1, j + 1, VRNA_DECOMP_EXT_EXT_EXT, sc->data);
      } else {
        contribution *= sc->exp_f(1, n, 1, i - 1, VRNA_DECOMP_EXT_EXT, sc->data);
      }
    } else if (j < n) {
      contribution *= sc->exp_f(1, n, j + 1, n, VRNA_DECOMP_EXT_EXT, sc->data);
    } else {
      contribution *= sc->exp_f(1, n, 1, n, VRNA_DECOMP_EXT_UP, sc->data);
    }
  }

  return contribution;
}


PRIVATE FLT_OR_DBL
contrib_ext_up_comparative(vrna_fold_compound_t *fc,
                             unsigned int         i,
                             unsigned int         j,
                             helper_arrays_t      dp_helpers)
{
  unsigned int      n, s, n_seq;
  FLT_OR_DBL        contribution, *scale;
  vrna_sc_t         **scs;

  n         = fc->length;
  n_seq     = fc->n_seq;
  scale     = fc->exp_matrices->scale;
  scs       = fc->scs;

  contribution = scale[j - i + 1];

  if (scs) {
    if (i > 1) {
      if (j < n) {
        for (s = 0; s < n_seq; s++)
          if (scs[s]->exp_f)
            contribution *= scs[s]->exp_f(1, n, i - 1, j + 1, VRNA_DECOMP_EXT_EXT_EXT, scs[s]->data);
      } else {
        for (s = 0; s < n_seq; s++)
          if (scs[s]->exp_f)
            contribution *= scs[s]->exp_f(1, n, 1, i - 1, VRNA_DECOMP_EXT_EXT, scs[s]->data);
      }
    } else if (j < n) {
      for (s = 0; s < n_seq; s++)
        if (scs[s]->exp_f)
          contribution *= scs[s]->exp_f(1, n, j + 1, n, VRNA_DECOMP_EXT_EXT, scs[s]->data);
    } else {
      for (s = 0; s < n_seq; s++)
        if (scs[s]->exp_f)
          contribution *= scs[s]->exp_f(1, n, 1, n, VRNA_DECOMP_EXT_UP, scs[s]->data);
    }
  }

  return contribution;
}


PRIVATE void
compute_bpp_internal(vrna_fold_compound_t *fc,
                     unsigned int         l,
                     FLT_OR_DBL           *Qmax,
                     int                  *ov,
                     helper_arrays_t      dp_helpers)
{
  unsigned char         type, type_2;
  char                  *ptype;
  short                 *S1;
  unsigned int          i, j, k, n, u1, u2, with_ud, *hc_up_int, mini, maxj;
  int                   ij, kl, *my_iindx, *jindx, *rtype;
  FLT_OR_DBL            temp, tmp2, *qb, *probs, *scale, *probs_int;
  double                max_real;
  vrna_exp_param_t      *pf_params;
  vrna_md_t             *md;
  vrna_hc_t             *hc;
  vrna_sc_t             *sc;
  vrna_ud_t             *domains_up;
  struct sc_int_exp_dat *sc_wrapper_int;

  vrna_array(vrna_ep_t) bp_correction;

  sc_wrapper_int  = &(dp_helpers->sc_wrapper_int);

  n           = (int)fc->length;
  ptype       = fc->ptype;
  S1          = fc->sequence_encoding;
  my_iindx    = fc->iindx;
  jindx       = fc->jindx;
  pf_params   = fc->exp_params;
  md          = &(pf_params->model_details);
  rtype       = &(md->rtype[0]);
  hc          = fc->hc;
  sc          = fc->sc;
  domains_up  = fc->domains_up;
  with_ud     = (domains_up && domains_up->exp_energy_cb) ? 1 : 0;
  hc_up_int   = hc->up_int;

  qb    = fc->exp_matrices->qb;
  probs = fc->exp_matrices->probs;
  scale = fc->exp_matrices->scale;

  bp_correction = dp_helpers->bp_correction;
  probs_int     = dp_helpers->probs_int;

  max_real = (sizeof(FLT_OR_DBL) == sizeof(float)) ? FLT_MAX : DBL_MAX;

  /* 2. bonding k,l as substem of 2:loop enclosed by i,j */
  for (k = 1; k < l; k++) {
    kl = my_iindx[k] - l;

    if (qb[kl] == 0.)
      continue;

    if (hc->mx[l * n + k] & VRNA_CONSTRAINT_CONTEXT_INT_LOOP_ENC) {
      type_2 = rtype[vrna_get_ptype(jindx[l] + k, ptype)];

      mini = 1;
      if (k > MAXLOOP + 2)
        mini = k - MAXLOOP - 1;

      for (i = mini; i <= k - 1; i++) {
        u1 = k - i - 1;
        if (hc_up_int[i + 1] < u1)
          continue;

        maxj = l + 1 + MAXLOOP - u1;

        if (maxj > n)
          maxj = n;

        if (maxj > l + 1 + hc_up_int[l + 1])
          maxj = l + 1 + hc_up_int[l + 1];

        u2 = 0;

        for (j = l + 1; j <= maxj; j++, u2++) {
          ij = my_iindx[i] - j;

          if (probs[ij] == 0.)
            continue;

          if (hc->eval_int(i, j, k, l, hc)) {
            int jij = jindx[j] + i;
            type  = vrna_get_ptype(jij, ptype);
            tmp2  = probs[ij] *
                    vrna_exp_E_internal(u1,
                                  u2,
                                  type,
                                  type_2,
                                  S1[i + 1],
                                  S1[j - 1],
                                  S1[k - 1],
                                  S1[l + 1],
                                  pf_params) *
                    scale[u1 + u2 + 2];

            if (sc_wrapper_int->pair)
              tmp2 *= sc_wrapper_int->pair(i, j, k, l, sc_wrapper_int);

            if (with_ud) {
              FLT_OR_DBL qql, qqr;

              qql = qqr = 0.;

              if (u1 > 0) {
                qql = domains_up->exp_energy_cb(fc,
                                                i + 1, k - 1,
                                                VRNA_UNSTRUCTURED_DOMAIN_INT_LOOP,
                                                domains_up->data);
              }

              if (u2 > 0) {
                qqr = domains_up->exp_energy_cb(fc,
                                                l + 1, j - 1,
                                                VRNA_UNSTRUCTURED_DOMAIN_INT_LOOP,
                                                domains_up->data);
              }

              temp  = tmp2;
              tmp2  += temp * qql;
              tmp2  += temp * qqr;
              tmp2  += temp * qql * qqr;
            }

            if ((sc) &&
                (sc->exp_f) &&
                (sc->bt)) {
              /* store probability correction for auxiliary pairs in internal loop motif */
              vrna_basepair_t *ptr, *aux_bps;
              aux_bps = sc->bt(i, j, k, l, VRNA_DECOMP_PAIR_IL, sc->data);
              for (ptr = aux_bps; ptr && ptr->i != 0; ptr++) {
                vrna_array_append(bp_correction,
                                  ((vrna_ep_t){
                                    .i  = ptr->i,
                                    .j  = ptr->j,
                                    .p  = tmp2 * qb[kl]
                                  }));
              }
              free(aux_bps);
            }

            probs[kl] += tmp2;

            probs_int[k] += tmp2;
          }
        }
      }
    }

    /* store actual probabilties for callback reporting */
    probs_int[k] *= qb[kl];

    if (probs[kl] > (*Qmax)) {
      (*Qmax) = probs[kl];
      if ((*Qmax) > max_real / 10.)
        vrna_log_warning("P close to overflow: %d %d %g %g\n",
                             k, l, probs[kl], qb[kl]);
    }

    if (probs[kl] >= max_real) {
      (*ov)++;
      probs[kl] = FLT_MAX;
    }
  }

  if (md->gquad)
    compute_gquad_prob_internal(fc, l);
}


PRIVATE void
compute_bpp_internal_comparative(vrna_fold_compound_t *fc,
                                 unsigned int         l,
                                 FLT_OR_DBL           *Qmax,
                                 int                  *ov,
                                 helper_arrays_t      dp_helpers)
{
  short                 **SS, **S5, **S3;
  unsigned int          i, j, k, n, mini, maxj, type, u1, u2, *tt, s, n_seq, **a2s, *hc_up_int;
  int                   ij, kl, *my_iindx, *jindx, *pscore;
  FLT_OR_DBL            tmp2, *qb, *probs, *scale, psc_exp, *probs_int;
  double                max_real, kTn;
  vrna_exp_param_t      *pf_params;
  vrna_md_t             *md;
  vrna_hc_t             *hc;
  struct sc_int_exp_dat *sc_wrapper_int;

  sc_wrapper_int  = &(dp_helpers->sc_wrapper_int);

  n         = (int)fc->length;
  n_seq     = fc->n_seq;
  pscore    = fc->pscore;
  SS        = fc->S;
  S5        = fc->S5;
  S3        = fc->S3;
  a2s       = fc->a2s;
  my_iindx  = fc->iindx;
  jindx     = fc->jindx;
  pf_params = fc->exp_params;
  md        = &(pf_params->model_details);
  hc        = fc->hc;
  hc_up_int = hc->up_int;

  qb    = fc->exp_matrices->qb;
  probs = fc->exp_matrices->probs;
  scale = fc->exp_matrices->scale;
  probs_int = dp_helpers->probs_int;

  kTn       = pf_params->kT / 10.;   /* kT in cal/mol  */
  max_real  = (sizeof(FLT_OR_DBL) == sizeof(float)) ? FLT_MAX : DBL_MAX;
  tt        = (unsigned int *)vrna_alloc(sizeof(unsigned int) * n_seq);

  /* 2. bonding k,l as substem of 2:loop enclosed by i,j */
  for (k = 1; k < l; k++) {
    kl = my_iindx[k] - l;

    if (qb[kl] == 0.)
      continue;

    if (hc->mx[l * n + k] & VRNA_CONSTRAINT_CONTEXT_INT_LOOP_ENC) {
      psc_exp = exp(pscore[jindx[l] + k] / kTn);

      for (s = 0; s < n_seq; s++)
        tt[s] = vrna_get_ptype_md(SS[s][l], SS[s][k], md);

      mini = 1;
      if (k > MAXLOOP + 2)
        mini = k - MAXLOOP - 1;

      for (i = mini; i <= k - 1; i++) {
        u1 = k - i - 1;
        if (hc_up_int[i + 1] < u1)
          continue;

        maxj = l + MAXLOOP + i + 2 - k;
        if (maxj > n)
          maxj = n;

        for (j = l + 1; j <= maxj; j++) {
          ij = my_iindx[i] - j;

          if (probs[ij] == 0.)
            continue;

          u2 = j - l - 1;

          if (hc_up_int[l + 1] < u2)
            break;

          if (hc->eval_int(i, j, k, l, hc)) {
            tmp2 = probs[ij] *
                   scale[u1 + u2 + 2] *
                   psc_exp;

            for (s = 0; s < n_seq; s++) {
              unsigned int u1_loc  = a2s[s][k - 1] - a2s[s][i];
              unsigned int u2_loc  = a2s[s][j - 1] - a2s[s][l];
              type  = vrna_get_ptype_md(SS[s][i], SS[s][j], md);
              tmp2  *= vrna_exp_E_internal(u1_loc,
                                     u2_loc,
                                     type,
                                     tt[s],
                                     S3[s][i],
                                     S5[s][j],
                                     S5[s][k],
                                     S3[s][l],
                                     pf_params);
            }

            if (sc_wrapper_int->pair)
              tmp2 *= sc_wrapper_int->pair(i, j, k, l, sc_wrapper_int);

            probs[kl] += tmp2;

            /* store actual probabilties for callback reporting */
            probs_int[k] += tmp2;
          }
        }
      }
    }

    /* store actual probabilties for callback reporting */
    probs_int[k] *= qb[kl];

    if (probs[kl] > (*Qmax)) {
      (*Qmax) = probs[kl];
      if ((*Qmax) > max_real / 10.)
        vrna_log_warning("P close to overflow: %d %d %g %g\n",
                             k, l, probs[kl], qb[kl]);
    }

    if (probs[kl] >= max_real) {
      (*ov)++;
      probs[kl] = FLT_MAX;
    }
  }

  free(tt);

  if (md->gquad)
    compute_gquad_prob_internal_comparative(fc, l);
}


PRIVATE void
compute_up_hp(vrna_fold_compound_t *fc,
              unsigned int         l,
              helper_arrays_t      dp_helpers)
{
  char                  **Ss;
  unsigned int          **a2s;
  short                 *S, *S2, **SS, **S5, **S3;
  int                   u, type, n_seq, s, max_u, ij;
  FLT_OR_DBL            q, qbt1, *probs, *scale, *puh, *qbh;
  vrna_exp_param_t      *P;
  vrna_md_t             *md;
  vrna_ud_t             *domains_up;
  unsigned int          i, j, n;
  int                   *my_iindx;
  vrna_hc_t             *hc;
  struct sc_hp_exp_dat  *sc_wrapper_hp;

  n           = fc->length;
  P           = fc->exp_params;
  md          = &(P->model_details);
  scale       = fc->exp_matrices->scale;
  domains_up  = fc->domains_up;
  probs       = fc->exp_matrices->probs;
  my_iindx  = fc->iindx;
  max_u     = dp_helpers->options.up_length;
  puh       = dp_helpers->access_hp;
  qbh       = dp_helpers->qbh;
  hc        = fc->hc;

  sc_wrapper_hp = &(dp_helpers->sc_wrapper_hp);

  if (l == n)
    return;

  j = l + 1;
  for (i = 1; i < l; i++) {
    ij = my_iindx[i] - j;
  
    if (probs[ij] < FLT_MIN)
      continue;

    if (hc->eval_hp(i, j, hc)) {
      q = 0.;
      switch (fc->type) {
        case VRNA_FC_TYPE_SINGLE:
          S     = fc->sequence_encoding;
          S2    = fc->sequence_encoding2;
          u     = j - i - 1;
          type  = vrna_get_ptype_md(S2[i], S2[j], md);

          q = exp_E_Hairpin(u, type, S[i + 1], S[j - 1], fc->sequence + i - 1, P);
          break;

        case VRNA_FC_TYPE_COMPARATIVE:
          SS    = fc->S;
          S5    = fc->S5;   /* S5[s][i] holds next base 5' of i in sequence s */
          S3    = fc->S3;   /* Sl[s][i] holds next base 3' of i in sequence s */
          Ss    = fc->Ss;
          a2s   = fc->a2s;
          n_seq = fc->n_seq;
          qbt1  = 1.;

          for (s = 0; s < n_seq; s++) {
            u = a2s[s][j - 1] - a2s[s][i];
            if (a2s[s][i] < 1)
              continue;

            type  = vrna_get_ptype_md(SS[s][i], SS[s][j], md);
            qbt1  *= exp_E_Hairpin(u, type, S3[s][i], S5[s][j], Ss[s] + a2s[s][i] - 1, P);
          }

          q = qbt1;
          break;

        default:
          break;
      }

      /* add soft constraints */
      if (sc_wrapper_hp->pair)
        q *= sc_wrapper_hp->pair(i, j, sc_wrapper_hp);

      if (domains_up && domains_up->exp_energy_cb) {
        /* we always consider both, bound and unbound state */
        q += q * domains_up->exp_energy_cb(fc,
                                           i + 1, j - 1,
                                           VRNA_UNSTRUCTURED_DOMAIN_HP_LOOP,
                                           domains_up->data);
      }

      q *= scale[j - i + 1];
      q *= probs[ij];
      qbh[i] += q;
    }
  }

  /* now, for all unpaired stretches of size u ending at l */
  for (u = 1; u <= MIN2(max_u, l - 1); u++)
    for (i = 1; i < l + 1 - u; i++)
      puh[u] += qbh[i];
}


PRIVATE void
compute_up_int(vrna_fold_compound_t *fc,
               unsigned int         l,
               helper_arrays_t      dp_helpers)
{
  short                 *S, *S1;
  unsigned int          n, u, max_u, k, i, j, p, u1, u2, turn, type, type2, with_ud;
  int                   ij, *my_iindx;
  FLT_OR_DBL            *pui, *probs, *scale, *qb, tmp;
  vrna_exp_param_t      *params;
  vrna_md_t             *md;
  vrna_hc_t             *hc;
  vrna_ud_t             *domains_up;
  struct sc_int_exp_dat *sc_wrapper_int;


  n           = fc->length;
  S           = fc->sequence_encoding2;
  S1          = fc->sequence_encoding;
  my_iindx    = fc->iindx;
  qb          = fc->exp_matrices->qb;
  probs       = fc->exp_matrices->probs;
  scale       = fc->exp_matrices->scale;
  max_u     = dp_helpers->options.up_length;
  pui       = dp_helpers->access_int;
  params    = fc->exp_params;
  md        = &(params->model_details);
  turn      = md->min_loop_size;
  domains_up  = fc->domains_up;
  with_ud     = (domains_up && domains_up->exp_energy_cb) ? 1 : 0;
  hc        = fc->hc;
  sc_wrapper_int  = &(dp_helpers->sc_wrapper_int);

  if (l == n)
    return;

  FLT_OR_DBL *q3 = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * (n + 1));
  FLT_OR_DBL *q5 = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * (n + 1));

  for (k = l - 1; k + MAXLOOP >= l; k--) {
    if (k == 0)
      break;

    for (j = l + 1; k + MAXLOOP + 1 >= j; j++) {
      if (j > n)
        break;

      u2 = j - k - 1;
      for (p = 2; p + turn < k; p++) {
        type2 = vrna_get_ptype_md(S[k], S[p], md);

        for (i = p - 1; i + MAXLOOP + 1 >= p + u2; i--) {
          if (i == 0)
            break;

          ij = my_iindx[i] - j;

          if (probs[ij] == 0.)
            continue;


          if (hc->eval_int(i, j, p, k, hc)) {
            u1   = p - i - 1;
            type = vrna_get_ptype_md(S[i], S[j], md);
            tmp  = probs[ij] *
                   vrna_exp_E_internal(u1,
                                  u2,
                                  type,
                                  type2,
                                  S1[i + 1],
                                  S1[j - 1],
                                  S1[p - 1],
                                  S1[k + 1],
                                  params) *
                    scale[u1 + u2 + 2];

            if (sc_wrapper_int->pair)
              tmp *= sc_wrapper_int->pair(i, j, p, k, sc_wrapper_int);

            if (with_ud) {
              /* here some unstructured domains stuff */
            }

            q3[k] += tmp *
                     qb[my_iindx[p] - k];
          } 
        }
      }
    }
  }

  for (i = l - 1; i + MAXLOOP >= l; i--) {
    if (i == 0)
      break;

    for (j = l + turn + 3; j <= n; j++) {
      ij = my_iindx[i] - j;

      if (probs[ij] < FLT_MIN)
        continue;

      type = vrna_get_ptype_md(S[i], S[j], md);

      for (k = j - 1; k + MAXLOOP >= j; k--) {
        if (l + turn + 2 > k)
          break;

        u2 = j - k - 1;

        for (p = l + 1; p + turn < k; p++) {
          u1 = p - i - 1;

          if (u1 + u2 > MAXLOOP)
            break;

          if (hc->eval_int(i, j, p, k, hc)) {
            type2 = vrna_get_ptype_md(S[k], S[p], md);
            tmp  = probs[ij] *
                   vrna_exp_E_internal(u1,
                                  u2,
                                  type,
                                  type2,
                                  S1[i + 1],
                                  S1[j - 1],
                                  S1[p - 1],
                                  S1[k + 1],
                                  params) *
                    scale[u1 + u2 + 2];

            if (sc_wrapper_int->pair)
              tmp *= sc_wrapper_int->pair(i, j, p, k, sc_wrapper_int);

            if (with_ud) {
              /* here some unstructured domains stuff */
            }

            q5[i] += tmp *
                     qb[my_iindx[p] - k];
          } 
        }
      }
    }  
  }

  for (u = 1; u <= max_u; u++) {
    unsigned int mink = (l > MAXLOOP) ? l - MAXLOOP - 1: 1;
    unsigned int maxk = (l > u) ? l - u : 0;

    for (k = mink; k <= maxk; k++)
      pui[u] += q3[k];

    unsigned int mini = (l > MAXLOOP) ? l - MAXLOOP - 1 : 1;
    unsigned int maxi = (l > u) ? l - u : 0;
    for (i = mini; i <= maxi; i++)
      pui[u] += q5[i];
      
  }
   
  free(q3);
  free(q5);
}


PRIVATE void
compute_up_mb(vrna_fold_compound_t *fc,
              unsigned int         l,
              helper_arrays_t      dp_helpers)
{
  short                 *S, *S1;
  unsigned int          n, u, max_u, i, j, turn, start, type;
  int                   *my_iindx, ij;
  FLT_OR_DBL            *pum, *qm, *qm2, *probs, *X, *Y, *Z, temp, up, qqm, *expMLbase, expMLclosing, *scale;
  vrna_exp_param_t      *params;
  vrna_md_t             *md;
  vrna_hc_t             *hc;
  struct sc_mb_exp_dat  *sc_wrapper;

  n         = fc->length;
  S         = fc->sequence_encoding2;
  S1        = fc->sequence_encoding;
  max_u     = dp_helpers->options.up_length;
  pum       = dp_helpers->access_mb;
  X         = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * (n + 1));
  Y         = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * (n + 1));
  Z         = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * (n + 1));
  probs     = fc->exp_matrices->probs;
  qm        = fc->exp_matrices->qm;
  qm2           = (fc->exp_matrices->qm2_real) ? fc->exp_matrices->qm2_real : dp_helpers->qm2_aux;
  params    = fc->exp_params;
  md        = &(params->model_details);
  turn          = md->min_loop_size;
  expMLbase = fc->exp_matrices->expMLbase;
  scale         = fc->exp_matrices->scale;
  expMLclosing  = params->expMLclosing;
  my_iindx  = fc->iindx;
  hc        = fc->hc;
  sc_wrapper  = &(dp_helpers->sc_wrapper_mb);

  if (l == n)
    return;

  /* fill X, Y, Z */
  for (i = 1; i < l; i++) {
    X[i] = Y[i] = Z[i] = 0.;

    for (j = l + 1; j <= n; j++) {
      if (hc->eval_mb(i + 1, j - 1, i + 1, l, VRNA_DECOMP_ML_ML, hc)) {
        ij    = my_iindx[i] - j;
        type  = vrna_get_ptype_md(S[j], S[i], md);
        temp  = probs[ij] *
                expMLbase[j - l - 1] *
                vrna_exp_E_multibranch_stem(type, S1[j - 1], S1[i + 1], params);

        if (sc_wrapper->red_ml)
          temp *= sc_wrapper->red_ml(i + 1, j - 1, i + 1, l, sc_wrapper);

        X[i] += temp;
      }
    }

    for (j = l + turn + 1; j <= n; j++) {
      if (hc->eval_mb(i + 1, j - 1, l, l + 1, VRNA_DECOMP_ML_ML_ML, hc)) {
        ij    = my_iindx[i] - j;
        type  = vrna_get_ptype_md(S[j], S[i], md);
        temp  = probs[ij] *
                qm[my_iindx[l + 1] - j + 1] *
                vrna_exp_E_multibranch_stem(type, S1[j - 1], S1[i + 1], params);

        if (sc_wrapper->decomp_ml)
          temp *= sc_wrapper->decomp_ml(i + 1, j - 1, l, l + 1, sc_wrapper);

        Y[i] += temp;
      }
    }

    for (j = l + 2 * turn + 1; j <= n; j++) {
      if (hc->eval_mb(i + 1, j - 1, l + 1, j - 1, VRNA_DECOMP_ML_ML, hc)) {
        ij    = my_iindx[i] - j;
        type  = vrna_get_ptype_md(S[j], S[i], md);
        temp  = probs[ij] *
                qm2[my_iindx[l + 1] - j + 1] *
                vrna_exp_E_multibranch_stem(type, S1[j - 1], S1[i + 1], params);

        if (sc_wrapper->red_ml)
          temp *= sc_wrapper->red_ml(i + 1, j - 1, l + 1, j - 1, sc_wrapper);

        Z[i] += temp;
      }
    }
  }

  for (u = 1; u <= max_u; u++) {
    /* process unpaired stretches of size u ending at l */
    if (l >= u) {
      start = l + 1 - u;
      temp = 0.;
      up = expMLbase[u];

      for (i = 1; i + turn < start; i++) {
        if (hc->eval_mb(i + 1, l, i + 1, start - 1, VRNA_DECOMP_ML_ML, hc)) {
          /* unpaired stretch is somewhere in the middle */
          qqm = qm[my_iindx[i + 1] - start + 1] *
                Y[i] *
                up;

          if (sc_wrapper->red_ml)
            qqm *= sc_wrapper->red_ml(i + 1, l, i + 1, start - 1, sc_wrapper);

          temp += qqm;

          /* unpaired stretch is right-most */
          qqm = qm2[my_iindx[i + 1] - start + 1] *
                X[i] *
                up;

          if (sc_wrapper->red_ml)
            qqm *= sc_wrapper->red_ml(i + 1, l, i + 1, start - 1, sc_wrapper);

          temp += qqm;
        }
      }

      for (i = 1; i < start; i++) {
        /* unpaired stretch is left-most */
        qqm = expMLbase[l - i] *
              Z[i];

        if (sc_wrapper->red_ml)
          qqm *= sc_wrapper->red_ml(i + 1, l, start, 1, sc_wrapper);

        temp += qqm;
      }

      pum[u] += temp *
                scale[2] *
                expMLclosing;
    }
  }

  free(X);
  free(Y);
  free(Z);
}

PRIVATE void
compute_bpp_multibranch(vrna_fold_compound_t  *fc,
                        unsigned int          l,
                        FLT_OR_DBL            *Qmax,
                        int                   *ov,
                        helper_arrays_t       dp_helpers)
{
  unsigned char             tt;
  char                      *ptype;
  short                     *S, *S1, s5, s3;
  unsigned int              *sn;
  int                       cnt, i, j, k, n, u, ii, ij, kl, lj, *my_iindx, *jindx,
                            *rtype, with_gquad, with_ud;
  FLT_OR_DBL                temp, ppp, prm_MLb, prmt, prmt1, *qb, *probs, *qm, *scale,
                            *expMLbase, expMLclosing, expMLstem, *probs_mb;
  double                    max_real;
  vrna_exp_param_t          *pf_params;
  vrna_md_t                 *md;
  vrna_ud_t                 *domains_up;
  vrna_hc_t                 *hc;
  struct sc_mb_exp_dat      *sc_wrapper;
  vrna_smx_csr(FLT_OR_DBL)  *q_gq;

  if ((l <= 1) ||
      (l == fc->length))
    return;

  n             = (int)fc->length;
  sn            = fc->strand_number;
  S             = fc->sequence_encoding2;
  S1            = fc->sequence_encoding;
  my_iindx      = fc->iindx;
  jindx         = fc->jindx;
  pf_params     = fc->exp_params;
  md            = &(pf_params->model_details);
  rtype         = &(md->rtype[0]);
  ptype         = fc->ptype;
  qb            = fc->exp_matrices->qb;
  qm            = fc->exp_matrices->qm;
  q_gq          = fc->exp_matrices->q_gq;
  probs         = fc->exp_matrices->probs;
  scale         = fc->exp_matrices->scale;
  expMLbase     = fc->exp_matrices->expMLbase;
  expMLclosing  = pf_params->expMLclosing;
  domains_up    = fc->domains_up;
  with_ud       = (domains_up && domains_up->exp_energy_cb) ? 1 : 0;
  with_gquad    = md->gquad;
  expMLstem     = (with_gquad) ? vrna_exp_E_multibranch_stem(0, -1, -1, pf_params) : 0;

  hc          = fc->hc;
  sc_wrapper  = &(dp_helpers->sc_wrapper_mb);
  probs_mb    = dp_helpers->probs_mb;

  prm_MLb   = 0.;
  max_real  = (sizeof(FLT_OR_DBL) == sizeof(float)) ? FLT_MAX : DBL_MAX;

  if (sn[l + 1] != sn[l]) {
    /* set prm_l to 0 to get prm_l1 in the next round to be 0 */
    for (i = 0; i <= n; i++)
      dp_helpers->prm_l[i] = 0;
  } else {
    for (k = 2; k < l; k++) {
      kl    = my_iindx[k] - l;
      i     = k - 1;
      prmt  = prmt1 = 0.0;

      ij  = my_iindx[i] - (l + 2);
      lj  = my_iindx[l + 1] - (l + 1);
      s3  = S1[i + 1];
      if (sn[k] == sn[i]) {
        for (j = l + 2; j <= n; j++, ij--, lj--) {
          if (hc->eval_mb(i, j, i + 1, j - 1, VRNA_DECOMP_PAIR_ML, hc)) {
            tt = vrna_get_ptype_md(S[j], S[i], md);

            /* which decomposition is covered here? =>
             * i + 1 = k < l < j:
             * (i,j)       -> enclosing pair
             * (k, l)      -> enclosed pair
             * (l+1, j-1)  -> multiloop part with at least one stem
             * a.k.a. (k,l) is left-most stem in multiloop closed by (k-1, j)
             */
            ppp = probs[ij] *
                  vrna_exp_E_multibranch_stem(tt, S1[j - 1], s3, pf_params) *
                  qm[lj];

            if (sc_wrapper->pair)
              ppp *= sc_wrapper->pair(i, j, sc_wrapper);

            prmt += ppp;
          }
        }

        ii  = my_iindx[i];  /* ii-j=[i,j]     */
        tt  = vrna_get_ptype(jindx[l + 1] + i, ptype);
        tt  = rtype[tt];
        if (hc->eval_mb(i, l + 1, i + 1, l, VRNA_DECOMP_PAIR_ML, hc)) {
          prmt1 = probs[ii - (l + 1)] *
                  vrna_exp_E_multibranch_stem(tt,
                               S1[l],
                               S1[i + 1],
                               pf_params) *
                  expMLclosing;

          if (sc_wrapper->pair)
            prmt1 *= sc_wrapper->pair(i, l + 1, sc_wrapper);
        }
      }

      prmt *= expMLclosing;

      dp_helpers->prml[i] = prmt;

      /* l+1 is unpaired */
      if (hc->eval_mb(k, l + 1, k, l, VRNA_DECOMP_ML_ML, hc)) {
        ppp = dp_helpers->prm_l1[i] *
              expMLbase[1];

        if (sc_wrapper->red_ml)
          ppp *= sc_wrapper->red_ml(k, l + 1, k, l, sc_wrapper);

        /* add contributions of MB loops where any unstructured domain starts at l+1 */
        if (with_ud) {
          for (cnt = 0; cnt < domains_up->uniq_motif_count; cnt++) {
            u = domains_up->uniq_motif_size[cnt];
            if (l + u < n) {
              if (hc->eval_mb(k, l + u, k, l, VRNA_DECOMP_ML_ML, hc)) {
                temp = domains_up->exp_energy_cb(fc,
                                                 l + 1,
                                                 l + u,
                                                 VRNA_UNSTRUCTURED_DOMAIN_MB_LOOP | VRNA_UNSTRUCTURED_DOMAIN_MOTIF,
                                                 domains_up->data) *
                       dp_helpers->pmlu[u][i] *
                       expMLbase[u];

                if (sc_wrapper->red_ml)
                  temp *= sc_wrapper->red_ml(k, l + u, k, l, sc_wrapper);

                ppp += temp;
              }
            }
          }
          dp_helpers->pmlu[0][i] = ppp + prmt1;
        }

        dp_helpers->prm_l[i] = ppp + prmt1;
      } else {
        /* skip configuration where l+1 is unpaired */
        dp_helpers->prm_l[i] = prmt1;

        if (with_ud)
          dp_helpers->pmlu[0][i] = prmt1;
      }

      if (hc->eval_mb(i, l, i + 1, l, VRNA_DECOMP_ML_ML, hc)) {
        ppp = prm_MLb *
              expMLbase[1];

        if (sc_wrapper->red_ml)
          ppp *= sc_wrapper->red_ml(i, l, i + 1, l, sc_wrapper);

        if (with_ud) {
          for (cnt = 0; cnt < domains_up->uniq_motif_count; cnt++) {
            u = domains_up->uniq_motif_size[cnt];
            if (1 + u <= i) {
              if (hc->eval_mb(i - u + 1, l, i + 1, l, VRNA_DECOMP_ML_ML, hc)) {
                temp = domains_up->exp_energy_cb(fc,
                                                 i - u + 1,
                                                 i,
                                                 VRNA_UNSTRUCTURED_DOMAIN_MB_LOOP | VRNA_UNSTRUCTURED_DOMAIN_MOTIF,
                                                 domains_up->data) *
                       dp_helpers->prm_MLbu[u] *
                       expMLbase[u];

                if (sc_wrapper->red_ml)
                  temp *= sc_wrapper->red_ml(i - u + 1, l, i + 1, l, sc_wrapper);

                ppp += temp;
              }
            }
          }
          dp_helpers->prm_MLbu[0] = ppp + dp_helpers->prml[i];
        }

        prm_MLb = ppp + dp_helpers->prml[i];
        /* same as:    prm_MLb = 0;
         * for (i=1; i<=k-1; i++) prm_MLb += prml[i]*expMLbase[k-i-1]; */
      } else {
        /* skip all configurations where i is unpaired */
        prm_MLb = dp_helpers->prml[i];

        if (with_ud)
          dp_helpers->prm_MLbu[0] = dp_helpers->prml[i];
      }

      dp_helpers->prml[i] = dp_helpers->prml[i] + dp_helpers->prm_l[i];

      tt = ptype[jindx[l] + k];

      if (with_gquad) {
        if ((!tt) &&
#ifndef VRNA_DISABLE_C11_FEATURES
            (vrna_smx_csr_get(q_gq, k, l, 0.) == 0.))
#else
            (vrna_smx_csr_FLT_OR_DBL_get(q_gq, k, l, 0.) == 0.))
#endif
          continue;
      } else {
        if (qb[kl] == 0.)
          continue;
      }

      temp = prm_MLb;

      if (sn[k] == sn[k - 1]) {
        if (sc_wrapper->decomp_ml) {
          for (i = 1; i <= k - 2; i++)
            temp += dp_helpers->prml[i] *
                    qm[my_iindx[i + 1] - (k - 1)] *
                    sc_wrapper->decomp_ml(i + 1, l, k - 1, k, sc_wrapper);
        } else {
          for (i = 1; i <= k - 2; i++)
            temp += dp_helpers->prml[i] *
                    qm[my_iindx[i + 1] - (k - 1)];
        }
      }

      s5  = ((k > 1) && (sn[k] == sn[k - 1])) ? S1[k - 1] : -1;
      s3  = ((l < n) && (sn[l + 1] == sn[l])) ? S1[l + 1] : -1;

      if ((with_gquad) &&
          (qb[kl] == 0.)) {
        temp *= expMLstem;
      } else if (hc->eval_mb(k, l, k, l, VRNA_DECOMP_ML_STEM, hc)) {
        if (tt == 0)
          tt = 7;

        temp *= vrna_exp_E_multibranch_stem(tt, s5, s3, pf_params);
      }

      if (sc_wrapper->red_stem)
        temp *= sc_wrapper->red_stem(k, l, k, l, sc_wrapper);

      probs[kl] += temp *
                   scale[2];

      /* store actual probabilties for callback reporting */
      probs_mb[k] = temp *
                    scale[2] *
                    qb[kl];

      if (probs[kl] > (*Qmax)) {
        (*Qmax) = probs[kl];
        if ((*Qmax) > max_real / 10.)
          vrna_log_warning("P close to overflow: %d %d %g %g\n",
                               k, l, probs[kl], qb[kl]);
      }

      if (probs[kl] >= max_real) {
        (*ov)++;
        probs[kl] = FLT_MAX;
      }

      /* rotate prm_MLbu entries required for unstructured domain feature */
      rotate_dp_helper_arrays_inner(dp_helpers);
    } /* end for (k=..) */
  }

  rotate_dp_helper_arrays_outer(dp_helpers);
}


PRIVATE void
compute_bpp_multibranch_comparative(vrna_fold_compound_t  *fc,
                                    unsigned int          l,
                                    FLT_OR_DBL            *Qmax,
                                    int                   *ov,
                                    helper_arrays_t       dp_helpers)
{
  unsigned char     tt;
  short             **S, **S5, **S3;
  unsigned int      **a2s, s, n_seq, *sn;
  int               i, j, k, n, ii, kl, ij, lj, *my_iindx, *jindx, *pscore, with_gquad;
  FLT_OR_DBL        temp, ppp, prm_MLb, prmt, prmt1, *qb, *probs, *qm, *scale,
                    *expMLbase, expMLclosing, expMLstem, *probs_mb;
  double            max_real, kTn;
  vrna_exp_param_t  *pf_params;
  vrna_md_t         *md;
  vrna_hc_t         *hc;
  vrna_sc_t         **scs;
  vrna_smx_csr(FLT_OR_DBL) *q_gq;

  if ((l <= 1) ||
      (l == fc->length))
    return;

  n             = (int)fc->length;
  n_seq         = fc->n_seq;
  S             = fc->S;
  S5            = fc->S5;
  S3            = fc->S3;
  a2s           = fc->a2s;
  sn            = fc->strand_number;
  pscore        = fc->pscore;
  my_iindx      = fc->iindx;
  jindx         = fc->jindx;
  pf_params     = fc->exp_params;
  md            = &(pf_params->model_details);
  qb            = fc->exp_matrices->qb;
  qm            = fc->exp_matrices->qm;
  q_gq          = fc->exp_matrices->q_gq;
  probs         = fc->exp_matrices->probs;
  scale         = fc->exp_matrices->scale;
  expMLbase     = fc->exp_matrices->expMLbase;
  expMLclosing  = pf_params->expMLclosing;
  with_gquad    = md->gquad;
  hc            = fc->hc;
  scs           = fc->scs;
  expMLstem     =
    (with_gquad) ? (FLT_OR_DBL)pow(vrna_exp_E_multibranch_stem(0, -1, -1, pf_params), (double)n_seq) : 0;

  kTn       = pf_params->kT / 10.;   /* kT in cal/mol  */
  prm_MLb   = 0.;
  max_real  = (sizeof(FLT_OR_DBL) == sizeof(float)) ? FLT_MAX : DBL_MAX;
  probs_mb  = dp_helpers->probs_mb;

  if (sn[l + 1] != sn[l]) {
    /* set prm_l to 0 to get prm_l1 in the next round to be 0 */
    for (i = 0; i <= n; i++)
      dp_helpers->prm_l[i] = 0;
  } else {
    for (k = 2; k < l; k++) {
      kl    = my_iindx[k] - l;
      i     = k - 1;
      prmt  = prmt1 = 0.0;

      ij  = my_iindx[i] - (l + 2);
      lj  = my_iindx[l + 1] - (l + 1);

      if (sn[k] == sn[i]) {
        for (j = l + 2; j <= n; j++, ij--, lj--) {
          if ((hc->mx[i * n + j] & VRNA_CONSTRAINT_CONTEXT_MB_LOOP) &&
              (sn[j] == sn[j - 1])) {
            /* which decomposition is covered here? =>
             * i + 1 = k < l < j:
             * (i,j)       -> enclosing pair
             * (k, l)      -> enclosed pair
             * (l+1, j-1)  -> multiloop part with at least one stem
             * a.k.a. (k,l) is left-most stem in multiloop closed by (k-1, j)
             */
            ppp = probs[ij] *
                  qm[lj];

            for (s = 0; s < n_seq; s++) {
              tt  = vrna_get_ptype_md(S[s][j], S[s][i], md);
              ppp *= vrna_exp_E_multibranch_stem(tt, S5[s][j], S3[s][i], pf_params);
            }

            if (scs) {
              for (s = 0; s < n_seq; s++) {
                if (scs[s]) {
                  if (scs[s]->exp_energy_bp)
                    ppp *= scs[s]->exp_energy_bp[jindx[j] + i];

                  /*
                   *  if(scs[s]->exp_f)
                   *    ppp *= scs[s]->exp_f(i, j, l+1, j-1, , scs[s]->data);
                   */
                }
              }
            }

            prmt += ppp;
          }
        }

        ii = my_iindx[i];   /* ii-j=[i,j]     */

        if (hc->mx[(l + 1) * n + i] & VRNA_CONSTRAINT_CONTEXT_MB_LOOP) {
          prmt1 = probs[ii - (l + 1)] *
                  (FLT_OR_DBL)pow(expMLclosing, (double)n_seq);

          for (s = 0; s < n_seq; s++) {
            tt    = vrna_get_ptype_md(S[s][l + 1], S[s][i], md);
            prmt1 *= vrna_exp_E_multibranch_stem(tt, S5[s][l + 1], S3[s][i], pf_params);
          }

          if (scs) {
            /* which decompositions are covered here? => (i, l+1) -> enclosing pair */
            for (s = 0; s < n_seq; s++) {
              if (scs[s]) {
                if (scs[s]->exp_energy_bp)
                  prmt1 *= scs[s]->exp_energy_bp[jindx[l + 1] + i];

                /*
                 *      if(sc->exp_f)
                 *        prmt1 *= sc->exp_f(i, l+1, k, l, , sc->data);
                 */
              }
            }
          }
        }
      }

      prmt *= (FLT_OR_DBL)pow(expMLclosing, (double)n_seq);

      dp_helpers->prml[i] = prmt;

      /* l+1 is unpaired */
      if (hc->up_ml[l + 1]) {
        ppp = dp_helpers->prm_l1[i] *
              expMLbase[1];

        if (scs) {
          for (s = 0; s < n_seq; s++) {
            if (scs[s]) {
              if (scs[s]->exp_energy_up)
                ppp *= scs[s]->exp_energy_up[a2s[s][l + 1]][1];

              /*
               *  if(scs[s]->exp_f)
               *    ppp *= scs[s]->exp_f(, scs[s]->data);
               */
            }
          }
        }

        dp_helpers->prm_l[i] = ppp + prmt1;
      } else {
        /* skip configuration where l+1 is unpaired */
        dp_helpers->prm_l[i] = prmt1;
      }

      /* i is unpaired */
      if (hc->up_ml[i]) {
        ppp = prm_MLb * expMLbase[1];

        if (scs) {
          for (s = 0; s < n_seq; s++) {
            if (scs[s]) {
              if (scs[s]->exp_energy_up)
                ppp *= scs[s]->exp_energy_up[a2s[s][i]][1];

              /*
               *  if(scs[s]->exp_f)
               *    ppp *= scs[s]->exp_f(, scs[s]->data);
               */
            }
          }
        }

        prm_MLb = ppp + dp_helpers->prml[i];
        /* same as:    prm_MLb = 0;
         * for (i=1; i<=k-1; i++) prm_MLb += prml[i]*expMLbase[k-i-1]; */
      } else {
        /* skip all configurations where i is unpaired */
        prm_MLb = dp_helpers->prml[i];
      }

      dp_helpers->prml[i] = dp_helpers->prml[i] + dp_helpers->prm_l[i];

      if (with_gquad) {
        if ((qb[kl] == 0.) &&
#ifndef VRNA_DISABLE_C11_FEATURES
            (vrna_smx_csr_get(q_gq, k, l, 0.) == 0.))
#else
            (vrna_smx_csr_FLT_OR_DBL_get(q_gq, k, l, 0.) == 0.))
#endif
          continue;
      } else {
        if (qb[kl] == 0.)
          continue;
      }

      temp = prm_MLb;

      if (sn[k] == sn[k - 1]) {
        for (i = 1; i <= k - 2; i++)
          if (sn[i + 1] == sn[i])
            temp += dp_helpers->prml[i] *
                    qm[my_iindx[i + 1] - (k - 1)];
      }

      if ((with_gquad) &&
          (qb[kl] == 0.)) {
        temp *= expMLstem;
      } else {
        if (hc->mx[l * n + k] & VRNA_CONSTRAINT_CONTEXT_MB_LOOP_ENC) {
          for (s = 0; s < n_seq; s++) {
            tt    = vrna_get_ptype_md(S[s][k], S[s][l], md);
            temp  *= vrna_exp_E_multibranch_stem(tt, S5[s][k], S3[s][l], pf_params);
          }
        }
      }

      probs[kl] += temp *
                   scale[2] *
                   exp(pscore[jindx[l] + k] / kTn);

      /* store actual probabilties for callback reporting */
      probs_mb[k] = temp *
                    scale[2] *
                    exp(pscore[jindx[l] + k] / kTn) *
                    qb[kl];

      if (probs[kl] > (*Qmax)) {
        (*Qmax) = probs[kl];
        if ((*Qmax) > max_real / 10.)
          vrna_log_warning("P close to overflow: %d %d %g %g\n",
                               k, l, probs[kl], qb[kl]);
      }

      if (probs[kl] >= max_real) {
        (*ov)++;
        probs[kl] = FLT_MAX;
      }

      /* rotate prm_MLbu entries required for unstructured domain feature */
      rotate_dp_helper_arrays_inner(dp_helpers);
    } /* end for (k=2..) */
  }

  rotate_dp_helper_arrays_outer(dp_helpers);
}


PRIVATE void
compute_gquad_prob_internal(vrna_fold_compound_t  *fc,
                            unsigned int          l)
{
  unsigned char     type;
  char              *ptype;
  short             *S1;
  unsigned int      i, j, k, n, u1, u2;
  int               ij, kl, *my_iindx, *jindx;
  FLT_OR_DBL        tmp2, qe, q_g, *probs, *scale;
  vrna_exp_param_t  *pf_params;
  vrna_md_t         *md;
  vrna_smx_csr(FLT_OR_DBL) *q_gq;

  n         = fc->length;
  S1        = fc->sequence_encoding;
  ptype     = fc->ptype;
  my_iindx  = fc->iindx;
  jindx     = fc->jindx;
  pf_params = fc->exp_params;
  md        = &(pf_params->model_details);
  q_gq      = fc->exp_matrices->q_gq;
  probs     = fc->exp_matrices->probs;
  scale     = fc->exp_matrices->scale;

  /* 2.5. bonding k,l as gquad enclosed by i,j */
  double *expintern = &(pf_params->expinternal[0]);

  if (l + 3 < n) {
    for (k = 2;
         k + VRNA_GQUAD_MIN_BOX_SIZE - 1 <= l;
         k++) {
      kl = my_iindx[k] - l;

#ifndef VRNA_DISABLE_C11_FEATURES
      q_g = vrna_smx_csr_get(q_gq, k, l, 0.);
#else
      q_g = vrna_smx_csr_FLT_OR_DBL_get(q_gq, k, l, 0.);
#endif

      if (q_g == 0.)
        continue;

      tmp2  = 0.;
      i     = k - 1;
      for (j = MIN2(l + MAXLOOP + 1, n);
           j > l + 3;
           j--) {
        ij    = my_iindx[i] - j;
        type  = (unsigned char)ptype[jindx[j] + i];
        if (!type)
          continue;

        u1    = j - l - 1;
        qe    = (type > 2) ? pf_params->expTermAU : 1.;

        if (md->dangles == 2)
          qe *= pf_params->expmismatchI[type][S1[i + 1]][S1[j - 1]];

        tmp2  += probs[ij] *
                 qe *
                 (FLT_OR_DBL)expintern[u1] *
                 scale[u1 + 2];

      }
      probs[kl] += tmp2;
    }
  }

  if (l + 1 < n) {
    for (k = (3 + VRNA_GQUAD_MAX_BOX_SIZE - 1 < l) ? l - VRNA_GQUAD_MAX_BOX_SIZE + 1 : 3;
         k + VRNA_GQUAD_MIN_BOX_SIZE - 1 <= l;
         k++) {
      kl = my_iindx[k] - l;

#ifndef VRNA_DISABLE_C11_FEATURES
      q_g = vrna_smx_csr_get(q_gq, k, l, 0.);
#else
      q_g = vrna_smx_csr_FLT_OR_DBL_get(q_gq, k, l, 0.);
#endif

      if (q_g == 0.)
        continue;

      tmp2 = 0.;
      for (i = (k > MAXLOOP + 1) ? k - MAXLOOP - 1 : 1;
           i + 1 < k;
           i++) {
        u1 = k - i - 1;
        for (j = l + 2;
             j <= MIN2(l + MAXLOOP - u1 + 1, n);
             j++) {
          ij    = my_iindx[i] - j;
          type  = (unsigned char)ptype[jindx[j] + i];
          if (!type)
            continue;

          u2    = j - l - 1;
          qe    = (type > 2) ? pf_params->expTermAU : 1.;
          if (md->dangles == 2)
            qe *= pf_params->expmismatchI[type][S1[i + 1]][S1[j - 1]];

          tmp2  += probs[ij] *
                   qe *
                   (FLT_OR_DBL)expintern[u1 + u2] *
                   scale[u1 + u2 + 2];
        }
      }
      probs[kl] += tmp2;
    }
  }

  if (l < n) {
    for (k = 4;
         k + VRNA_GQUAD_MIN_BOX_SIZE - 1 <= l;
         k++) {
      kl = my_iindx[k] - l;

#ifndef VRNA_DISABLE_C11_FEATURES
      q_g = vrna_smx_csr_get(q_gq, k, l, 0.);
#else
      q_g = vrna_smx_csr_FLT_OR_DBL_get(q_gq, k, l, 0.);
#endif

      if (q_g == 0.)
        continue;

      tmp2  = 0.;
      j     = l + 1;
      for (i = (k > MAXLOOP + 1) ? k - MAXLOOP - 1 : 1;
           i + 3 < k;
           i++) {
        ij    = my_iindx[i] - j;
        type  = (unsigned char)ptype[jindx[j] + i];
        if (!type)
          continue;

        u2    = k - i - 1;
        qe    = (type > 2) ? pf_params->expTermAU : 1.;
        if (md->dangles == 2)
          qe *= pf_params->expmismatchI[type][S1[i + 1]][S1[j - 1]];

        tmp2  += probs[ij] *
                 qe *
                 (FLT_OR_DBL)expintern[u2] *
                 scale[u2 + 2];

      }
      probs[kl] += tmp2;
    }
  }
}


PRIVATE void
compute_gquad_prob_internal_comparative(vrna_fold_compound_t  *fc,
                                        unsigned int          l)
{
  unsigned char     type;
  short             **S, **S5, **S3;
  unsigned int      i, j, k, n, u1, u2, u1_local, u2_local, **a2s, s, n_seq;
  int               ij, kl, *my_iindx;
  FLT_OR_DBL        tmp2, qe, q_g, *qb, *probs, *scale;
  vrna_exp_param_t  *pf_params;
  vrna_md_t         *md;
  vrna_smx_csr(FLT_OR_DBL) *q_gq;

  n         = (int)fc->length;
  n_seq     = fc->n_seq;
  S         = fc->S;
  S5        = fc->S5;
  S3        = fc->S3;
  a2s       = fc->a2s;
  my_iindx  = fc->iindx;
  pf_params = fc->exp_params;
  q_gq      = fc->exp_matrices->q_gq;
  qb        = fc->exp_matrices->qb;
  probs     = fc->exp_matrices->probs;
  scale     = fc->exp_matrices->scale;
  md        = &(pf_params->model_details);

  /* 2.5. bonding k,l as gquad enclosed by i,j */
  double *expintern = &(pf_params->expinternal[0]);

  if (l + 3 < n) {
    for (k = 2; k + VRNA_GQUAD_MIN_BOX_SIZE <= l + 1; k++) {
      kl = my_iindx[k] - l;

#ifndef VRNA_DISABLE_C11_FEATURES
      q_g = vrna_smx_csr_get(q_gq, k, l, 0.);
#else
      q_g = vrna_smx_csr_FLT_OR_DBL_get(q_gq, k, l, 0.);
#endif

      if (q_g == 0.)
        continue;

      tmp2  = 0.;
      i     = k - 1;
      for (j = MIN2(l + MAXLOOP + 1, n); j > l + 3; j--) {
        ij = my_iindx[i] - j;
        if (qb[ij] == 0.)
          continue;

        qe  = 1.;
        u1  = j - l - 1;

        for (s = 0; s < n_seq; s++) {
          type      = vrna_get_ptype_md(S[s][i], S[s][j], md);
          u1_local  = a2s[s][j - 1] - a2s[s][l];
          qe        *= (FLT_OR_DBL)expintern[u1_local];

          if (md->dangles == 2)
            qe *= (FLT_OR_DBL)pf_params->expmismatchI[type][S3[s][i]][S5[s][j]];

          if (type > 2)
            qe *= (FLT_OR_DBL)pf_params->expTermAU;
        }

        tmp2 += probs[ij] *
                qe *
                scale[u1 + 2];
      }
      probs[kl] += tmp2;
    }
  }

  if (l + 1 < n) {
    for (k = 3; k + VRNA_GQUAD_MIN_BOX_SIZE <= l + 1; k++) {
      kl = my_iindx[k] - l;

#ifndef VRNA_DISABLE_C11_FEATURES
      q_g = vrna_smx_csr_get(q_gq, k, l, 0.);
#else
      q_g = vrna_smx_csr_FLT_OR_DBL_get(q_gq, k, l, 0.);
#endif

      if (q_g == 0.)
        continue;

      tmp2 = 0.;
      unsigned int mini = 1;
      if (k > MAXLOOP + 2)
        mini = k - MAXLOOP - 1;

      for (i = mini; i + 2 <= k; i++) {
        u1 = k - i - 1;
        for (j = l + 2; j <= MIN2(l + MAXLOOP - u1 + 1, n); j++) {
          ij = my_iindx[i] - j;
          if (qb[ij] == 0.)
            continue;

          qe  = 1.;
          u2  = j - l - 1;

          for (s = 0; s < n_seq; s++) {
            type      = vrna_get_ptype_md(S[s][i], S[s][j], md);
            u1_local  = a2s[s][k - 1] - a2s[s][i];
            u2_local  = a2s[s][j - 1] - a2s[s][l];
            qe        *= (FLT_OR_DBL)expintern[u1_local + u2_local];

            if (md->dangles == 2)
              qe *= (FLT_OR_DBL)pf_params->expmismatchI[type][S3[s][i]][S5[s][j]];

            if (type > 2)
              qe *= (FLT_OR_DBL)pf_params->expTermAU;
          }

          tmp2 += probs[ij] *
                  qe *
                  scale[u1 + u2 + 2];
        }
      }
      probs[kl] += tmp2;
    }
  }

  if (l < n) {
    for (k = 4; k + VRNA_GQUAD_MIN_BOX_SIZE <= l + 1; k++) {
      kl = my_iindx[k] - l;
#ifndef VRNA_DISABLE_C11_FEATURES
      q_g = vrna_smx_csr_get(q_gq, k, l, 0.);
#else
      q_g = vrna_smx_csr_FLT_OR_DBL_get(q_gq, k, l, 0.);
#endif

      if (q_g == 0.)
        continue;

      tmp2  = 0.;
      j     = l + 1;
      unsigned int mini = 1;
      if (k > MAXLOOP + 2)
        mini = k - MAXLOOP - 1;

      for (i = mini; i < k - 3; i++) {
        ij = my_iindx[i] - j;
        if (qb[ij] == 0.)
          continue;

        qe  = 1.;
        u2  = k - i - 1;

        for (s = 0; s < n_seq; s++) {
          type      = vrna_get_ptype_md(S[s][i], S[s][j], md);
          u2_local  = a2s[s][k - 1] - a2s[s][i];
          qe        *= (FLT_OR_DBL)expintern[u2_local];

          if (md->dangles == 2)
            qe *= (FLT_OR_DBL)pf_params->expmismatchI[type][S3[s][i]][S5[s][j]];

          if (type > 2)
            qe *= (FLT_OR_DBL)pf_params->expTermAU;
        }

        tmp2 += probs[ij] *
                qe *
                scale[u2 + 2];
      }
      probs[kl] += tmp2;
    }
  }
}


PRIVATE void
multistrand_update_Y5(vrna_fold_compound_t  *fc,
                      unsigned int          l,
                      helper_arrays_t       dp_helpers)
{
  short             *S, *S1;
  unsigned int      i, s, *se, *sn, type, n, end;
  int               *my_iindx;
  FLT_OR_DBL        *probs, *q, *scale, qtmp;
  FLT_OR_DBL        *Y5;
  FLT_OR_DBL        **Y5p;
  vrna_md_t         *md;
  vrna_exp_param_t  *pf_params;
  struct sc_ext_exp_dat   *sc_wrapper;
  sc_ext_exp_cb           sc_red_stem;
  sc_ext_exp_split        sc_split;

  n         = fc->length;
  sn        = fc->strand_number;
  se        = fc->strand_end;
  my_iindx  = fc->iindx;
  q         = fc->exp_matrices->q;
  probs     = fc->exp_matrices->probs;
  scale     = fc->exp_matrices->scale;
  pf_params = fc->exp_params;
  md        = &(pf_params->model_details);
  S         = fc->sequence_encoding2;
  S1        = fc->sequence_encoding;
  sc_wrapper  = &(dp_helpers->sc_wrapper_ext);
  sc_red_stem = sc_wrapper->red_stem;
  sc_split    = sc_wrapper->split;
  Y5          = dp_helpers->Y5;
  Y5p         = dp_helpers->Y5p;

  /* compute Y5 for all strands */
  for (s = 0; s < fc->strands; s++) {
    unsigned int j;

    Y5[s] = 0;

    if ((se[s] < l) &&
        (sn[l] == sn[l + 1])) {
      /* pre-compute newly available Y5p[s][j] with j == l + 1 */
      end = se[s];
      j   = l + 1;

      Y5p[s][j] = 0.;

      if (probs[my_iindx[end] - j] > 0) {
        type  = vrna_get_ptype_md(S[j], S[end], md);
        qtmp  = probs[my_iindx[end] - j] *
                vrna_exp_E_exterior_stem(type,
                                    S1[j - 1],
                                    -1,
                                    pf_params) *
                scale[2];

        if (sc_red_stem)
          qtmp *= sc_red_stem(j, end, j, end, sc_wrapper);

        Y5p[s][j] += qtmp;
      }

      for (i = 1; i < end; i++) {
        if ((probs[my_iindx[i] - j] > 0) &&
            (sn[i] == sn[i + 1])) {
          type  = vrna_get_ptype_md(S[j], S[i], md);
          qtmp  = probs[my_iindx[i] - j] *
                  vrna_exp_E_exterior_stem(type,
                  S1[j - 1],
                  S1[i + 1],
                  pf_params) *
                  q[my_iindx[i + 1] - end] *
                  scale[2];

          if (sc_red_stem)
            qtmp *= sc_red_stem(j, i, j, i, sc_wrapper);
          if (sc_split)
            qtmp *= sc_split(i, end, i + 1, sc_wrapper);

          Y5p[s][j] += qtmp;
        }
      }

      if ((probs[my_iindx[i] - j] > 0) &&
          (sn[i] == sn[i + 1])) {
        type  = vrna_get_ptype_md(S[j], S[i], md);
        qtmp  = probs[my_iindx[i] - j] *
                vrna_exp_E_exterior_stem(type,
                                    S1[j - 1],
                                    S1[i + 1],
                                    pf_params) *
                scale[2];

        if (sc_red_stem)
          qtmp *= sc_red_stem(j, i, j, i, sc_wrapper);

        Y5p[s][j] += qtmp;
      }

      /* recompute Y5[s] */
      Y5[s] += Y5p[s][l + 1];
      for (j = l + 2; j <= n; j++) {
        qtmp = q[my_iindx[l + 1] - (j - 1)] *
               Y5p[s][j];

        if (sc_split)
          qtmp *= sc_split(l + 1, j, j, sc_wrapper);

        Y5[s] += qtmp;
      }
    }
  }
}


PRIVATE void
multistrand_update_Y3(vrna_fold_compound_t  *fc,
                      unsigned int          l,
                      helper_arrays_t       dp_helpers)
{
  short             *S, *S1;
  unsigned int      i, j, k, s, n, start, type, *ss, *sn;
  int               *my_iindx;
  FLT_OR_DBL        *q, *probs, *scale, qtmp;
  FLT_OR_DBL        **Y3;
  FLT_OR_DBL        **Y3p;
  vrna_md_t         *md;
  vrna_exp_param_t  *pf_params;
  struct sc_ext_exp_dat   *sc_wrapper;
  sc_ext_exp_cb           sc_red_stem;
  sc_ext_exp_split        sc_split;

  n         = fc->length;
  sn        = fc->strand_number;
  ss        = fc->strand_start;
  S         = fc->sequence_encoding2;
  S1        = fc->sequence_encoding;
  my_iindx  = fc->iindx;
  q         = fc->exp_matrices->q;
  probs     = fc->exp_matrices->probs;
  scale     = fc->exp_matrices->scale;
  pf_params = fc->exp_params;
  md        = &(pf_params->model_details);
  sc_wrapper  = &(dp_helpers->sc_wrapper_ext);
  sc_red_stem = sc_wrapper->red_stem;
  sc_split    = sc_wrapper->split;
  Y3          = dp_helpers->Y3;
  Y3p         = dp_helpers->Y3p;

  for (s = 0; s < fc->strands; s++) {
    start = ss[s];
    if (start == l + 1) {
      for (i = 1; i < start; i++) {
        Y3p[s][i] = 0;

        if (sn[i] == sn[i + 1]) {
          if (probs[my_iindx[i] - start] > 0) {
            type = vrna_get_ptype_md(S[start], S[i], md);
            qtmp = probs[my_iindx[i] - start] *
                         vrna_exp_E_exterior_stem(type,
                                             -1,
                                             S1[i + 1],
                                             pf_params) *
                         scale[2];

            if (sc_red_stem)
              qtmp *= sc_red_stem(start, i, start, i, sc_wrapper);

            Y3p[s][i] += qtmp;
          }

          for (j = start + 1; j <= n; j++) {
            if ((probs[my_iindx[i] - j] > 0) &&
                (sn[j - 1] == sn[j])) {
              type = vrna_get_ptype_md(S[j], S[i], md);
              qtmp = probs[my_iindx[i] - j] *
                           vrna_exp_E_exterior_stem(type,
                                               S1[j - 1],
                                               S1[i + 1],
                                               pf_params) *
                           q[my_iindx[start] - (j - 1)] *
                           scale[2];

              if (sc_red_stem)
                qtmp *= sc_red_stem(j, i, j, i, sc_wrapper);
              if (sc_split)
                qtmp *= sc_split(start, j, j, sc_wrapper);

              Y3p[s][i] += qtmp;
            }
          }
        }
      }

      for (k = 1; k < start; k++) {
        Y3[s][k] = 0.;

        if (sn[k - 1] == sn[k]) {
          for (i = 1; i < k - 1; i++) {
            if (sn[i] == sn[i + 1]) {
              qtmp = q[my_iindx[i + 1] - (k - 1)] *
                     Y3p[s][i];

              if (sc_split)
                qtmp *= sc_split(i, k - 1, i + 1, sc_wrapper);

              Y3[s][k] += qtmp;
            }
          }

          Y3[s][k] += Y3p[s][k - 1];
        }
      }
    }
  }
}


PRIVATE void
multistrand_contrib(vrna_fold_compound_t  *fc,
                    unsigned int          l,
                    helper_arrays_t       dp_helpers,
                    FLT_OR_DBL            *Qmax,
                    int                   *ov)
{
  short             *S, *S1, s5, s3;
  unsigned int      k, s, *sn, *se, *ss, end, start, type;
  int               *my_iindx, kl;
  FLT_OR_DBL        *q, *qb, *probs, tmp, qtmp;
  FLT_OR_DBL        *Y5;
  FLT_OR_DBL        **Y3;
  vrna_md_t         *md;
  vrna_exp_param_t  *pf_params;
  struct sc_ext_exp_dat   *sc_wrapper;
  sc_ext_exp_cb           sc_red_stem;
  sc_ext_exp_split        sc_split;

  if ((l <= 1) ||
      (l == fc->length))
    return;

  sn        = fc->strand_number;
  ss        = fc->strand_start;
  se        = fc->strand_end;
  S         = fc->sequence_encoding2;
  S1        = fc->sequence_encoding;
  pf_params = fc->exp_params;
  md        = &(pf_params->model_details);
  my_iindx  = fc->iindx;
  q         = fc->exp_matrices->q;
  qb        = fc->exp_matrices->qb;
  probs     = fc->exp_matrices->probs;
  sc_wrapper  = &(dp_helpers->sc_wrapper_ext);
  sc_red_stem = sc_wrapper->red_stem;
  sc_split    = sc_wrapper->split;
  Y5          = dp_helpers->Y5;
  Y3          = dp_helpers->Y3;

  for (k = l - 1; k > 1; k--) {
    kl = my_iindx[k] - l;
    if (qb[kl] > 0) {
      tmp = 0.;
      for (s = 0; s < fc->strands; s++) {
        end   = se[s];
        start = ss[s];
        if (end == k - 1)
          tmp += Y5[s];
        else if ((end < k - 1) &&
                 (sn[k - 1] == sn[k]))
          tmp += Y5[s] *
                 q[my_iindx[end + 1] - (k - 1)];
        else if (start == l + 1)
          tmp += Y3[s][k];
        else if ((start > l + 1) &&
                 (sn[l] == sn[l + 1]))
          tmp += Y3[s][k] *
                 q[my_iindx[l + 1] - (start - 1)];
      }

      type      = vrna_get_ptype_md(S[k], S[l], md);
      s5        = (sn[k - 1] == sn[k]) ? S1[k - 1] : -1;
      s3        = (sn[l] == sn[l + 1]) ? S1[l + 1] : -1;
      qtmp      = vrna_exp_E_exterior_stem(type,
                                      s5,
                                      s3,
                                      pf_params);

      if (sc_red_stem)
        qtmp *= sc_red_stem(k, l, k, l, sc_wrapper);

      probs[kl] += tmp * qtmp;
    }
  }
}


PRIVATE INLINE void
ud_outside_ext_loops(vrna_fold_compound_t *vc)
{
  unsigned int  i, j, n, cnt, *hc_up;
  int           u, *motif_list;
  FLT_OR_DBL    *q1k, *qln, temp, *scale;
  vrna_sc_t     *sc;
  vrna_ud_t     *domains_up;

  n           = vc->length;
  q1k         = vc->exp_matrices->q1k;
  qln         = vc->exp_matrices->qln;
  scale       = vc->exp_matrices->scale;
  hc_up       = vc->hc->up_ext;
  domains_up  = vc->domains_up;
  sc          = vc->sc;

  for (i = 1; i <= n; i++) {
    motif_list = vrna_ud_get_motif_size_at(vc, i, VRNA_UNSTRUCTURED_DOMAIN_EXT_LOOP);

    /* 1. Exterior loops */
    if (motif_list) {
      cnt = 0;
      while (-1 != (u = motif_list[cnt])) {
        j = i + (unsigned int)u - 1;
        if (j <= n) {
          if (hc_up[i] >= (unsigned int)u) {
            temp  = q1k[i - 1] * qln[j + 1] / q1k[n];
            temp  *= domains_up->exp_energy_cb(vc,
                                               i,
                                               j,
                                               VRNA_UNSTRUCTURED_DOMAIN_EXT_LOOP | VRNA_UNSTRUCTURED_DOMAIN_MOTIF,
                                               domains_up->data);

            if (sc)
              if (sc->exp_energy_up)
                temp *= sc->exp_energy_up[i][u];

            temp *= scale[u];

            if (temp > 0.) {
              domains_up->probs_add(vc,
                                    i,
                                    j,
                                    VRNA_UNSTRUCTURED_DOMAIN_EXT_LOOP | VRNA_UNSTRUCTURED_DOMAIN_MOTIF,
                                    temp,
                                    domains_up->data);
            }
          }
        }

        cnt++;
      }
    }

    free(motif_list);
  }
}


PRIVATE INLINE void
ud_outside_hp_loops(vrna_fold_compound_t *vc)
{
  unsigned int  i, j, k, l, n, cnt, *hc_up;
  int           kl, *my_iindx, u, *motif_list;
  FLT_OR_DBL    temp, outside, exp_motif_en, *probs, q1, q2;

  vrna_ud_t     *domains_up, *ud_bak;

  n           = vc->length;
  my_iindx    = vc->iindx;
  probs       = vc->exp_matrices->probs;
  hc_up       = vc->hc->up_hp;
  domains_up  = vc->domains_up;

  for (i = 1; i <= n; i++) {
    motif_list = vrna_ud_get_motif_size_at(vc, i, VRNA_UNSTRUCTURED_DOMAIN_HP_LOOP);

    /* 2. Hairpin loops */
    if (motif_list) {
      cnt = 0;
      while (-1 != (u = motif_list[cnt])) {
        outside = 0.;
        j       = i + u - 1;
        if (j < n) {
          if (hc_up[i] >= (unsigned int)u) {
            exp_motif_en = domains_up->exp_energy_cb(vc,
                                                     i,
                                                     j,
                                                     VRNA_UNSTRUCTURED_DOMAIN_HP_LOOP | VRNA_UNSTRUCTURED_DOMAIN_MOTIF,
                                                     domains_up->data);

            /*
             * compute the contribution of all hairpins with
             * bound motif
             */
            for (k = 1; k < i; k++)
              for (l = j + 1; l <= n; l++) {
                kl = my_iindx[k] - l;
                if (probs[kl] > 0.) {
                  ud_bak          = vc->domains_up;
                  vc->domains_up  = NULL;
                  temp            = vrna_exp_eval_hairpin(vc, k, l, VRNA_EVAL_LOOP_DEFAULT);
                  vc->domains_up  = ud_bak;

                  /* add contribution of motif */
                  if (temp > 0.) {
                    temp *= exp_motif_en * probs[kl];

                    q1 = q2 = 0.;
                    /* add contributions of other motifs in remaining unpaired segments */
                    if ((i - k - 1) > 0) {
                      q1 = domains_up->exp_energy_cb(vc,
                                                     k + 1, i - 1,
                                                     VRNA_UNSTRUCTURED_DOMAIN_HP_LOOP,
                                                     domains_up->data);
                    }

                    if ((l - j - 1) > 0) {
                      q2 = domains_up->exp_energy_cb(vc,
                                                     j + 1, l - 1,
                                                     VRNA_UNSTRUCTURED_DOMAIN_HP_LOOP,
                                                     domains_up->data);
                    }

                    outside += temp;
                    outside += temp * q1;
                    outside += temp * q1 * q2;
                    outside += temp * q2;
                  }
                }
              }
          }
        }

        if (outside > 0.) {
          domains_up->probs_add(vc,
                                i, j,
                                VRNA_UNSTRUCTURED_DOMAIN_HP_LOOP | VRNA_UNSTRUCTURED_DOMAIN_MOTIF,
                                outside,
                                domains_up->data);
        }

        cnt++;
      }
    }

    free(motif_list);
  }
}


PRIVATE INLINE void
ud_outside_hp_loops2(vrna_fold_compound_t *vc)
{
  unsigned int  i, j, k, l, u1, u2, n, m;
  int           kl, *my_iindx, u;
  FLT_OR_DBL    temp, outside, exp_motif_en, *probs, q1, q2, **qq_ud, **pp_ud;

  vrna_ud_t   *domains_up, *ud_bak;

  n           = vc->length;
  my_iindx    = vc->iindx;
  probs       = vc->exp_matrices->probs;
  domains_up  = vc->domains_up;

  qq_ud = (FLT_OR_DBL **)vrna_alloc(sizeof(FLT_OR_DBL *) * (n + 1));
  for (k = 0; k < (unsigned int)domains_up->uniq_motif_count; k++) {
    qq_ud[k]  = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * (n + 1));
    u         = domains_up->uniq_motif_size[k];
    for (i = 1; i <= n - u + 1; i++) {
      qq_ud[k][i] = domains_up->exp_energy_cb(vc,
                                              i,
                                              i + u - 1,
                                              VRNA_UNSTRUCTURED_DOMAIN_HP_LOOP | VRNA_UNSTRUCTURED_DOMAIN_MOTIF,
                                              domains_up->data);
    }
  }

  pp_ud = (FLT_OR_DBL **)vrna_alloc(sizeof(FLT_OR_DBL *) * (n + 1));
  for (k = 0; k < (unsigned int)domains_up->uniq_motif_count; k++)
    pp_ud[k] = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * (n + 1));

  for (k = 1; k < n; k++) {
    for (l = k + 1; l <= n; l++) {
      kl = my_iindx[k] - l;
      if (probs[kl] > 0.) {
        ud_bak          = vc->domains_up;
        vc->domains_up  = NULL;
        temp            = vrna_exp_eval_hairpin(vc, k, l, VRNA_EVAL_LOOP_DEFAULT);
        vc->domains_up  = ud_bak;
        temp            *= probs[kl];

        if (temp > 0.) {
          for (m = 0; m < (unsigned int)domains_up->uniq_motif_count; m++) {
            u = domains_up->uniq_motif_size[m];
            for (u1 = 0, u2 = l - k - u - 1, i = k + 1, j = k + u; j < l; i++, j++, u1++, u2--) {
              exp_motif_en  = qq_ud[m][i];
              q1            = q2 = 1.;
              if (u1 > 0) {
                q1 += domains_up->exp_energy_cb(vc,
                                                k + 1, i - 1,
                                                VRNA_UNSTRUCTURED_DOMAIN_HP_LOOP,
                                                domains_up->data);
              }

              if (u2 > 0) {
                q2 += domains_up->exp_energy_cb(vc,
                                                j + 1, l - 1,
                                                VRNA_UNSTRUCTURED_DOMAIN_HP_LOOP,
                                                domains_up->data);
              }

              outside     = temp * q1 * q2 * exp_motif_en;
              pp_ud[m][i] += outside;
            }
          }
        }
      }
    }
  }

  for (k = 0; k < (unsigned int)domains_up->uniq_motif_count; k++) {
    u = domains_up->uniq_motif_size[k];
    /* actually store the results */
    for (i = 1; i <= n - u + 1; i++) {
      if (pp_ud[k][i] > 0.) {
        domains_up->probs_add(vc,
                              i, i + u - 1,
                              VRNA_UNSTRUCTURED_DOMAIN_HP_LOOP | VRNA_UNSTRUCTURED_DOMAIN_MOTIF,
                              pp_ud[k][i],
                              domains_up->data);
      }
    }
    free(qq_ud[k]);
    free(pp_ud[k]);
  }
  free(qq_ud);
  free(pp_ud);
}


PRIVATE INLINE void
ud_outside_int_loops(vrna_fold_compound_t *vc)
{
  unsigned int  i, j, k, l, p, q, n, cnt, *hc_up, kmin,
                pmax, qmin, lmax;
  int           pq, kl, u, *motif_list, *my_iindx;
  FLT_OR_DBL    temp, q1, q2, q3, exp_motif_en, outside,
                *probs, *qb;
  vrna_ud_t     *domains_up, *ud_bak;

  n           = vc->length;
  my_iindx    = vc->iindx;
  qb          = vc->exp_matrices->qb;
  probs       = vc->exp_matrices->probs;
  hc_up       = vc->hc->up_int;
  domains_up  = vc->domains_up;

  for (i = 2; i <= n; i++) {
    motif_list = vrna_ud_get_motif_size_at(vc, i, VRNA_UNSTRUCTURED_DOMAIN_INT_LOOP);

    /* 3. Interior loops */
    if (motif_list) {
      cnt = 0;
      while (-1 != (u = motif_list[cnt])) {
        outside = 0.;
        j       = i + (unsigned int)u - 1;

        if (j < n) {
          if (hc_up[i] >= (unsigned int)u) {
            exp_motif_en = domains_up->exp_energy_cb(vc,
                                                     i,
                                                     j,
                                                     VRNA_UNSTRUCTURED_DOMAIN_INT_LOOP | VRNA_UNSTRUCTURED_DOMAIN_MOTIF,
                                                     domains_up->data);

            /* 3.1 motif is within 5' loop */
            kmin = 1;

            if (j > MAXLOOP + 2)
              kmin  = j - MAXLOOP - 1;

            for (k = kmin; k < i; k++) {
              pmax  = k + MAXLOOP + 1;
              pmax  = MIN2(pmax, n);
              for (p = j + 1; p < n; p++)
                for (q = p + 1; q < n; q++) {
                  pq = my_iindx[p] - q;
                  if (qb[pq] == 0)
                    continue;

                  lmax  = k + MAXLOOP + q - p + 2;
                  lmax  = MIN2(lmax, n);
                  for (l = q + 1; l <= lmax; l++) {
                    kl = my_iindx[k] - l;
                    if (probs[kl] > 0.) {
                      ud_bak          = vc->domains_up;
                      vc->domains_up  = NULL;
                      temp            = vrna_exp_eval_internal(vc, k, l, p, q, VRNA_EVAL_LOOP_DEFAULT);
                      vc->domains_up  = ud_bak;

                      if (temp > 0.) {
                        temp *= probs[kl] * qb[pq] * exp_motif_en;

                        q1 = q2 = q3 = 0.;
                        if ((l - q - 1) > 0) {
                          q1 = domains_up->exp_energy_cb(vc,
                                                         q + 1, l - 1,
                                                         VRNA_UNSTRUCTURED_DOMAIN_INT_LOOP,
                                                         domains_up->data);
                        }

                        if ((i - k - 1) > 0) {
                          q2 = domains_up->exp_energy_cb(vc,
                                                         k + 1, i - 1,
                                                         VRNA_UNSTRUCTURED_DOMAIN_INT_LOOP,
                                                         domains_up->data);
                        }

                        if ((p - j - 1) > 0) {
                          q3 = domains_up->exp_energy_cb(vc,
                                                         j + 1, p - 1,
                                                         VRNA_UNSTRUCTURED_DOMAIN_INT_LOOP,
                                                         domains_up->data);
                        }

                        outside += temp;
                        outside += temp * q1;
                        outside += temp * q1 * q2;
                        outside += temp * q1 * q2 * q3;
                        outside += temp * q2;
                        outside += temp * q2 * q3;
                        outside += temp * q3;
                      }
                    }
                  }
                }
            }

            /* 3.2 motif is within 3' loop */
            for (k = 1; k < i - 2; k++) {
              pmax  = k + i + MAXLOOP - j;
              pmax  = MIN2(pmax, n);
              for (p = k + 1; p <= pmax; p++) {
                qmin = p + 1;

                if (j > k + MAXLOOP + 2)
                  qmin  = p + j - k - MAXLOOP - 1;

                for (q = i - 1; q >= qmin; q--) {
                  pq = my_iindx[p] - q;
                  if (qb[pq] == 0.)
                    continue;

                  lmax  = k + q + MAXLOOP + 2 - p;
                  lmax  = MIN2(lmax, n);
                  for (l = j + 1; l < lmax; l++) {
                    kl = my_iindx[k] - l;
                    if (probs[kl] > 0.) {
                      ud_bak          = vc->domains_up;
                      vc->domains_up  = NULL;
                      temp            = vrna_exp_eval_internal(vc, k, l, p, q, VRNA_EVAL_LOOP_DEFAULT);
                      vc->domains_up  = ud_bak;

                      if (temp > 0.) {
                        FLT_OR_DBL q1, q2, q3;
                        temp *= probs[kl] * qb[pq] * exp_motif_en;

                        q1 = q2 = q3 = 0.;
                        if ((l - j - 1) > 0) {
                          q1 = domains_up->exp_energy_cb(vc,
                                                         j + 1, l - 1,
                                                         VRNA_UNSTRUCTURED_DOMAIN_INT_LOOP,
                                                         domains_up->data);
                        }

                        if ((i - q - 1) > 0) {
                          q2 = domains_up->exp_energy_cb(vc,
                                                         q + 1,
                                                         i - 1,
                                                         VRNA_UNSTRUCTURED_DOMAIN_INT_LOOP,
                                                         domains_up->data);
                        }

                        if ((p - k - 1) > 0) {
                          q3 = domains_up->exp_energy_cb(vc,
                                                         k + 1,
                                                         p - 1,
                                                         VRNA_UNSTRUCTURED_DOMAIN_INT_LOOP,
                                                         domains_up->data);
                        }

                        outside += temp;
                        outside += temp * q1;
                        outside += temp * q1 * q2;
                        outside += temp * q1 * q2 * q3;
                        outside += temp * q2;
                        outside += temp * q2 * q3;
                        outside += temp * q3;
                      }
                    }
                  }
                }
              }
            }
          }
        }

        if (outside > 0.) {
          domains_up->probs_add(vc,
                                i, j,
                                VRNA_UNSTRUCTURED_DOMAIN_INT_LOOP | VRNA_UNSTRUCTURED_DOMAIN_MOTIF,
                                outside,
                                domains_up->data);
        }

        cnt++;
      }
    }

    free(motif_list);
  }
}


PRIVATE INLINE void
ud_outside_int_loops2(vrna_fold_compound_t *vc)
{
  unsigned char *hard_constraints;
  unsigned int  i, j, k, l, p, q, u, n, pmax, qmin,
                u1, u2, uu1, uu2, u2_max, m;
  int           pq, kl, *my_iindx;
  FLT_OR_DBL    temp, q5, q3, exp_motif_en, outside,
                *probs, *qb, qq1, qq2, *qqk, *qql, *qqp, **qq_ud, **pp_ud, temp5,
                temp3;
  vrna_ud_t     *domains_up, *ud_bak;

  n                 = vc->length;
  my_iindx          = vc->iindx;
  qb                = vc->exp_matrices->qb;
  probs             = vc->exp_matrices->probs;
  hard_constraints  = vc->hc->mx;
  domains_up        = vc->domains_up;

  qqk = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * (n + 1));
  qql = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * (n + 1));
  qqp = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * (n + 1));

  qq_ud = (FLT_OR_DBL **)vrna_alloc(sizeof(FLT_OR_DBL *) * (n + 1));
  for (k = 0; k < (unsigned int)domains_up->uniq_motif_count; k++) {
    qq_ud[k]  = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * (n + 1));
    u         = domains_up->uniq_motif_size[k];
    for (i = 1; i <= n - u + 1; i++) {
      qq_ud[k][i] = domains_up->exp_energy_cb(vc,
                                              i,
                                              i + u - 1,
                                              VRNA_UNSTRUCTURED_DOMAIN_INT_LOOP | VRNA_UNSTRUCTURED_DOMAIN_MOTIF,
                                              domains_up->data);
    }
  }

  pp_ud = (FLT_OR_DBL **)vrna_alloc(sizeof(FLT_OR_DBL *) * (n + 1));
  for (k = 0; k < (unsigned int)domains_up->uniq_motif_count; k++)
    pp_ud[k] = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * (n + 1));

  for (k = 1; k < n; k++) {
    for (l = k + 1; l <= MIN2(k + MAXLOOP, n); l++) {
      qqk[l] = domains_up->exp_energy_cb(vc,
                                         k + 1, l,
                                         VRNA_UNSTRUCTURED_DOMAIN_INT_LOOP,
                                         domains_up->data);
    }
    for (l = k + 1 + 3; l <= n; l++) {
      kl = my_iindx[k] - l;
      if (probs[kl] == 0.)
        continue;

      if (hard_constraints[k * n + l] & VRNA_CONSTRAINT_CONTEXT_INT_LOOP) {
        unsigned int imin = k;
        if (l > MAXLOOP + 2)
          imin = l - MAXLOOP - 1;

        for (i = l - 1; i > imin; i--) {
          qql[i] = domains_up->exp_energy_cb(vc,
                                             i, l - 1,
                                             VRNA_UNSTRUCTURED_DOMAIN_INT_LOOP,
                                             domains_up->data);
        }

        pmax  = k + MAXLOOP + 1;
        pmax  = MIN2(pmax, l);
        for (p = k + 1; p < pmax; p++) {
          u1      = p - k - 1;
          u2_max  = MAXLOOP - u1;
          qmin    = p + 1;

          if (l > p + 2 + u2_max)
            qmin = l - 1 - u2_max;

          for (i = p - 1; i > k; i--) {
            qqp[i] = domains_up->exp_energy_cb(vc,
                                               i, p - 1,
                                               VRNA_UNSTRUCTURED_DOMAIN_INT_LOOP,
                                               domains_up->data);
          }

          q5 = 1.;
          if (u1 > 0)
            q5 += qqk[p - 1];

          for (q = qmin; q < l; q++) {
            pq = my_iindx[p] - q;

            if (hard_constraints[p * n + q] & VRNA_CONSTRAINT_CONTEXT_INT_LOOP_ENC) {
              u2              = l - q - 1;
              ud_bak          = vc->domains_up;
              vc->domains_up  = NULL;
              temp            = vrna_exp_eval_internal(vc, k, l, p, q, VRNA_EVAL_LOOP_DEFAULT);
              vc->domains_up  = ud_bak;
              temp            *= probs[kl] * qb[pq];

              q3 = 1.;
              if (u2 > 0)
                q3 += qql[q + 1];

              temp5 = temp * q3;
              temp3 = temp * q5;

              /* loop over all available motifs */
              for (m = 0; m < (unsigned int)domains_up->uniq_motif_count; m++) {
                u = domains_up->uniq_motif_size[m];
                for (i = k + 1, j = k + u; j < p; i++, j++) {
                  /* ud in 5' loop */
                  exp_motif_en  = qq_ud[m][i];
                  uu1           = i - k - 1;
                  uu2           = p - j - 1;
                  qq1           = 1.;
                  qq2           = 1.;
                  if (uu1 > 0)
                    qq1 += qqk[i - 1];

                  if (uu2 > 0)
                    qq2 += qqp[j + 1];

                  outside     = temp5 * qq1 * qq2 * exp_motif_en;
                  pp_ud[m][i] += outside;
                }

                for (i = q + 1, j = q + u; j < l; i++, j++) {
                  /* ud in 3' loop */
                  exp_motif_en  = qq_ud[m][i];
                  uu1           = i - q - 1;
                  uu2           = l - j - 1;
                  qq1           = 1.;
                  qq2           = 1.;
                  if (uu1 > 0) {
                    qq1 += domains_up->exp_energy_cb(vc,
                                                     q + 1, i - 1,
                                                     VRNA_UNSTRUCTURED_DOMAIN_INT_LOOP,
                                                     domains_up->data);
                  }

                  if (uu2 > 0)
                    qq2 += qql[j + 1];

                  outside     = temp3 * qq1 * qq2 * exp_motif_en;
                  pp_ud[m][i] += outside;
                }
              }
            }
          }
        }
      }
    }
  }

  free(qqk);
  free(qql);
  free(qqp);
  for (k = 0; k < (unsigned int)domains_up->uniq_motif_count; k++) {
    u = domains_up->uniq_motif_size[k];
    /* actually store the results */
    for (i = 1; i <= n - u + 1; i++) {
      if (pp_ud[k][i] > 0.) {
        domains_up->probs_add(vc,
                              i, i + u - 1,
                              VRNA_UNSTRUCTURED_DOMAIN_INT_LOOP | VRNA_UNSTRUCTURED_DOMAIN_MOTIF,
                              pp_ud[k][i],
                              domains_up->data);
      }
    }
    free(qq_ud[k]);
    free(pp_ud[k]);
  }
  free(qq_ud);
  free(pp_ud);
}


PRIVATE INLINE void
ud_outside_mb_loops(vrna_fold_compound_t *vc)
{
  unsigned char     *hc;
  char              *ptype;
  short             *S;
  unsigned int      i, j, k, l, n, cnt, *hc_up, tt, ud_max_size, up;
  int               kl, jkl, *my_iindx, u, *motif_list, *jindx, *rtype;
  FLT_OR_DBL        temp, *scale, outside, exp_motif_en, *probs, *qb, *qm, expMLclosing,
                    *expMLbase, *qmli, exp_motif_ml_left, exp_motif_ml_right;
  vrna_exp_param_t  *pf_params;
  vrna_md_t         *md;
  vrna_sc_t         *sc;
  vrna_ud_t         *domains_up;

  n             = vc->length;
  S             = vc->sequence_encoding;
  my_iindx      = vc->iindx;
  jindx         = vc->jindx;
  ptype         = vc->ptype;
  pf_params     = vc->exp_params;
  md            = &(vc->exp_params->model_details);
  qb            = vc->exp_matrices->qb;
  qm            = vc->exp_matrices->qm;
  probs         = vc->exp_matrices->probs;
  scale         = vc->exp_matrices->scale;
  hc_up         = vc->hc->up_ml;
  hc            = vc->hc->mx;
  domains_up    = vc->domains_up;
  sc            = vc->sc;
  rtype         = &(md->rtype[0]);
  expMLbase     = vc->exp_matrices->expMLbase;
  expMLclosing  = pf_params->expMLclosing;

  for (ud_max_size = 0, u = 0; u < domains_up->uniq_motif_count; u++)
    if (ud_max_size < domains_up->uniq_motif_size[u])
      ud_max_size = domains_up->uniq_motif_size[u];

  for (i = 1; i <= n; i++) {
    motif_list = vrna_ud_get_motif_size_at(vc, i, VRNA_UNSTRUCTURED_DOMAIN_MB_LOOP);

    /* 4. Multibranch loops */
    if (motif_list) {
      cnt = 0;
      while (-1 != (u = motif_list[cnt])) {
        outside = 0.;
        j       = i + u - 1;
        if (j < n) {
          if (hc_up[i] >= (unsigned int)u) {
            exp_motif_en = domains_up->exp_energy_cb(vc,
                                                     i,
                                                     j,
                                                     VRNA_UNSTRUCTURED_DOMAIN_MB_LOOP | VRNA_UNSTRUCTURED_DOMAIN_MOTIF,
                                                     domains_up->data);

            exp_motif_en *= expMLbase[u];

            if (sc)
              if (sc->exp_energy_up)
                exp_motif_en *= sc->exp_energy_up[i][u];

            temp = 0;

            /* 4.1 Motif [i:j] is somewhere in between branching stems */
            for (l = j + 1; l <= n; l++) {
              for (k = i - 1; k > 0; k--) {
                kl = my_iindx[k] - l;
                if (probs[kl] > 0.) {
                  if (hc[n * l + k] & VRNA_CONSTRAINT_CONTEXT_MB_LOOP) {
                    /* respect hard constraints */
                    FLT_OR_DBL qqq;
                    jkl = jindx[l] + k;
                    tt  = rtype[vrna_get_ptype(jkl, ptype)];
                    qqq = probs[kl] *
                          qm[my_iindx[k + 1] - (i - 1)] *
                          qm[my_iindx[j + 1] - (l - 1)] *
                          vrna_exp_E_multibranch_stem(tt, S[l - 1], S[k + 1], pf_params) *
                          expMLclosing *
                          scale[2];

                    if (sc)
                      if (sc->exp_energy_bp)
                        qqq *= sc->exp_energy_bp[jkl];

                    temp += qqq;
                  }
                }
              }
            }

            outside += temp
                       * exp_motif_en;

            /* 4.2 Motif is in left-most unpaired stretch of multiloop */
            FLT_OR_DBL **qm1ui =
              (FLT_OR_DBL **)vrna_alloc(sizeof(FLT_OR_DBL *) * (ud_max_size + 1));

            for (l = 0; l <= ud_max_size; l++)
              qm1ui[l] = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * (n + 2));

            exp_motif_ml_left = 0.;
            for (l = j + 1; l <= n; l++) {
              FLT_OR_DBL  lqq = 0.;
              FLT_OR_DBL  rqq = 0.;
              for (k = i - 1; k > 0; k--) {
                up  = i - k - 1;
                kl  = my_iindx[k] - l;
                if ((hc[n * l + k] & VRNA_CONSTRAINT_CONTEXT_MB_LOOP) &&
                    (probs[kl] > 0.) &&
                    (hc_up[k + 1] >= up)) {
                  int jkl = jindx[l] + k;
                  tt    = rtype[vrna_get_ptype(jkl, ptype)];
                  temp  = probs[kl] *
                          expMLbase[up] *
                          vrna_exp_E_multibranch_stem(tt, S[l - 1], S[k + 1], pf_params) *
                          expMLclosing *
                          scale[2];
                  if (sc) {
                    if (sc->exp_energy_bp)
                      temp *= sc->exp_energy_bp[jkl];

                    if (sc->exp_energy_up)
                      temp *= sc->exp_energy_up[k + 1][up];
                  }

                  lqq += temp;
                  lqq += temp *
                         domains_up->exp_energy_cb(vc,
                                                   k + 1, i - 1,
                                                   VRNA_UNSTRUCTURED_DOMAIN_MB_LOOP,
                                                   domains_up->data);
                }
              }

              for (u = j + 1; (unsigned int)u < l; u++) {
                /* 1st, l-1 is unpaired */
                if (hc_up[l - 1]) {
                  temp = qm1ui[1][u] * expMLbase[1];
                  if (sc)
                    if (sc->exp_energy_up)
                      temp *= sc->exp_energy_up[l - 1][1];

                  qm1ui[0][u] = temp;
                } else {
                  qm1ui[0][u] = 0.;
                }

                /* 2nd, l-1 is the final position of another motif [p:l-1] */
                for (cnt = 0; cnt < (unsigned int)domains_up->uniq_motif_count; cnt++) {
                  unsigned int size = (unsigned int)domains_up->uniq_motif_size[cnt];
                  if (((unsigned int)u < l - size) &&
                      (hc_up[l - size] >= size)) {
                    temp = qm1ui[size][u] *
                           expMLbase[size] *
                           domains_up->exp_energy_cb(vc,
                                                     l - size,
                                                     l - 1,
                                                     VRNA_UNSTRUCTURED_DOMAIN_MB_LOOP | VRNA_UNSTRUCTURED_DOMAIN_MOTIF,
                                                     domains_up->data);
                    if (sc)
                      if (sc->exp_energy_up)
                        temp *= sc->exp_energy_up[l - size][size];

                    qm1ui[0][u] += temp;
                  }
                }

                /* 3rd, l - 1 pairs with u */
                if (hc[n * (l - 1) + u] & VRNA_CONSTRAINT_CONTEXT_MB_LOOP_ENC) {
                  tt    = vrna_get_ptype(jindx[l - 1] + u, ptype);
                  temp  = qb[my_iindx[u] - (l - 1)] *
                          vrna_exp_E_multibranch_stem(tt, S[u - 1], S[l], pf_params);

                  qm1ui[0][u] += temp;
                }

                rqq += qm[my_iindx[j + 1] - (u - 1)] *
                       qm1ui[0][u];
              }

              /* finally, compose contribution */
              exp_motif_ml_left += lqq * rqq;

              /* rotate auxiliary arrays */
              FLT_OR_DBL *tmp = qm1ui[ud_max_size];
              for (cnt = ud_max_size; cnt > 0; cnt--)
                qm1ui[cnt] = qm1ui[cnt - 1];
              qm1ui[0] = tmp;
            }

            /* cleanup memory */
            for (l = 0; l <= ud_max_size; l++)
              free(qm1ui[l]);
            free(qm1ui);

            outside += exp_motif_ml_left
                       * exp_motif_en;

            /* 4.3 Motif is in right-most unpaired stretch of multiloop */
            qmli                = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * n);
            exp_motif_ml_right  = 0.;
            for (k = i - 1; k > 0; k--) {
              FLT_OR_DBL  lqq = 0.;
              FLT_OR_DBL  rqq = 0;

              /* update qmli[k] = qm1[k,i-1] */
              for (qmli[k] = 0., u = k + 1; (unsigned int)u < i; u++) {
                /* respect hard constraints */
                if (hc[n * u + k] & VRNA_CONSTRAINT_CONTEXT_MB_LOOP_ENC) {
                  up = (i - 1) - (u + 1) + 1;
                  if (hc_up[u + 1] >= up) {
                    temp = qb[my_iindx[k] - u] *
                           expMLbase[up];

                    /* add soft constraints */
                    if (sc)
                      if (sc->exp_energy_up)
                        temp *= sc->exp_energy_up[u + 1][up];

                    qmli[k] += temp;

                    /* add contributions of other motifs within [u+1:i-1] */
                    qmli[k] += temp *
                               domains_up->exp_energy_cb(vc,
                                                         u + 1, i - 1,
                                                         VRNA_UNSTRUCTURED_DOMAIN_MB_LOOP,
                                                         domains_up->data);
                  }
                }
              }

              for (u = k; (unsigned int)u < i; u++)
                lqq += qm[my_iindx[k + 1] - (u - 1)] *
                       qmli[u];

              for (l = j + 1; l <= n; l++) {
                kl = my_iindx[k] - l;
                if (hc[n * k + l] & VRNA_CONSTRAINT_CONTEXT_MB_LOOP) {
                  int jkl;
                  jkl = jindx[l] + k;
                  tt  = rtype[vrna_get_ptype(jkl, ptype)];
                  up  = l - j - 1;
                  if (hc_up[j + 1] >= up) {
                    temp = probs[kl] *
                           vrna_exp_E_multibranch_stem(tt, S[l - 1], S[k + 1], pf_params) *
                           expMLclosing *
                           scale[2] *
                           expMLbase[up];

                    if (sc) {
                      if (sc->exp_energy_bp)
                        temp *= sc->exp_energy_bp[jkl];

                      if (sc->exp_energy_up)
                        temp *= sc->exp_energy_up[j + 1][up];
                    }

                    rqq += temp;

                    /* add contributions of other motifs within [j+1:l-1] */
                    rqq += temp *
                           domains_up->exp_energy_cb(vc,
                                                     j + 1, l - 1,
                                                     VRNA_UNSTRUCTURED_DOMAIN_MB_LOOP,
                                                     domains_up->data);
                  }
                }
              }
              exp_motif_ml_right += rqq * lqq;
            }
            free(qmli);
            qmli = NULL;

            outside += exp_motif_ml_right *
                       exp_motif_en;
          }
        }

        if (outside > 0.) {
          domains_up->probs_add(vc,
                                i, j,
                                VRNA_UNSTRUCTURED_DOMAIN_MB_LOOP | VRNA_UNSTRUCTURED_DOMAIN_MOTIF,
                                outside,
                                domains_up->data);
        }

        cnt++;
      }
    }

    free(motif_list);
  }
}


PRIVATE INLINE void
ud_outside_mb_loops2(vrna_fold_compound_t *vc)
{
  unsigned char     *hc;
  char              *ptype;
  short             *S;
  unsigned int      i, j, k, l, n, cnt, *hc_up, tt, up, ud_max_size;
  int               kl, jkl, *my_iindx, u, *motif_list, *jindx, *rtype;
  FLT_OR_DBL        temp, *scale, outside, exp_motif_en, *probs, *qb, *qm,
                    expMLclosing, *expMLbase, *qmli, exp_motif_ml_left,
                    exp_motif_ml_right, *qqi, *qqj, *qqmi, *qqmj;
  vrna_exp_param_t  *pf_params;
  vrna_md_t         *md;
  vrna_sc_t         *sc;
  vrna_ud_t         *domains_up;

  n             = vc->length;
  S             = vc->sequence_encoding;
  my_iindx      = vc->iindx;
  jindx         = vc->jindx;
  ptype         = vc->ptype;
  pf_params     = vc->exp_params;
  md            = &(vc->exp_params->model_details);
  qb            = vc->exp_matrices->qb;
  qm            = vc->exp_matrices->qm;
  probs         = vc->exp_matrices->probs;
  scale         = vc->exp_matrices->scale;
  hc_up         = vc->hc->up_ml;
  hc            = vc->hc->mx;
  domains_up    = vc->domains_up;
  sc            = vc->sc;
  rtype         = &(md->rtype[0]);
  expMLbase     = vc->exp_matrices->expMLbase;
  expMLclosing  = pf_params->expMLclosing;

  for (ud_max_size = u = 0; u < domains_up->uniq_motif_count; u++)
    if (ud_max_size < domains_up->uniq_motif_size[u])
      ud_max_size = domains_up->uniq_motif_size[u];

  qqi   = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * (n + 1));
  qqj   = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * (n + 1));
  qqmi  = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * (n + 1));
  qqmj  = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * (n + 1));

  for (i = 1; i <= n; i++) {
    motif_list = vrna_ud_get_motif_size_at(vc, i, VRNA_UNSTRUCTURED_DOMAIN_MB_LOOP);

    /* 4. Multibranch loops */
    if (motif_list) {
      cnt = 0;
      while (-1 != (u = motif_list[cnt])) {
        outside = 0.;
        j       = i + u - 1;
        if (j < n) {
          if (hc_up[i] >= (unsigned int)u) {
            exp_motif_en = domains_up->exp_energy_cb(vc,
                                                     i,
                                                     j,
                                                     VRNA_UNSTRUCTURED_DOMAIN_MB_LOOP | VRNA_UNSTRUCTURED_DOMAIN_MOTIF,
                                                     domains_up->data);
            for (k = 1; k < i; k++) {
              qqi[k] = domains_up->exp_energy_cb(vc,
                                                 k, i - 1,
                                                 VRNA_UNSTRUCTURED_DOMAIN_MB_LOOP,
                                                 domains_up->data);
              qqmi[k] = qm[my_iindx[k] - (i - 1)];
            }
            for (l = j + 1; l <= n; l++) {
              qqj[l] = domains_up->exp_energy_cb(vc,
                                                 j + 1, l,
                                                 VRNA_UNSTRUCTURED_DOMAIN_MB_LOOP,
                                                 domains_up->data);
              qqmj[k] = qm[my_iindx[j + 1] - l];
            }
            exp_motif_en *= expMLbase[u];

            if (sc)
              if (sc->exp_energy_up)
                exp_motif_en *= sc->exp_energy_up[i][u];

            temp = 0;

            /* 4.1 Motif [i:j] is somewhere in between branching stems */
            for (l = j + 1; l <= n; l++) {
              for (k = i - 1; k > 0; k--) {
                kl = my_iindx[k] - l;
                if (probs[kl] > 0.) {
                  jkl = jindx[l] + k;
                  if (hc[l * n + k] & VRNA_CONSTRAINT_CONTEXT_MB_LOOP) {
                    /* respect hard constraints */
                    FLT_OR_DBL qqq;
                    tt  = rtype[vrna_get_ptype(jkl, ptype)];
                    qqq = probs[kl] *
                          qqmi[k + 1] *
                          qqmj[l - 1] *
                          vrna_exp_E_multibranch_stem(tt, S[l - 1], S[k + 1], pf_params) *
                          expMLclosing *
                          scale[2];

                    if (sc)
                      if (sc->exp_energy_bp)
                        qqq *= sc->exp_energy_bp[jkl];

                    temp += qqq;
                  }
                }
              }
            }

            outside += temp *
                       exp_motif_en;

            /* 4.2 Motif is in left-most unpaired stretch of multiloop */
            FLT_OR_DBL **qm1ui =
              (FLT_OR_DBL **)vrna_alloc(sizeof(FLT_OR_DBL *) * (ud_max_size + 1));

            for (l = 0; l <= ud_max_size; l++)
              qm1ui[l] = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * (n + 2));

            exp_motif_ml_left = 0.;
            for (l = j + 1; l <= n; l++) {
              FLT_OR_DBL  lqq = 0.;
              FLT_OR_DBL  rqq = 0.;
              for (k = i - 1; k > 0; k--) {
                up  = i - k - 1;
                kl  = my_iindx[k] - l;
                if ((hc[l * n + k] & VRNA_CONSTRAINT_CONTEXT_MB_LOOP) &&
                    (probs[kl] > 0.) &&
                    (hc_up[k + 1] >= up)) {
                  int jkl = jindx[l] + k;
                  tt    = rtype[vrna_get_ptype(jkl, ptype)];
                  temp  = probs[kl] *
                          expMLbase[up] *
                          vrna_exp_E_multibranch_stem(tt, S[l - 1], S[k + 1], pf_params) *
                          expMLclosing *
                          scale[2];
                  if (sc) {
                    if (sc->exp_energy_bp)
                      temp *= sc->exp_energy_bp[jkl];

                    if (sc->exp_energy_up)
                      temp *= sc->exp_energy_up[k + 1][up];
                  }

                  lqq += temp;
                  lqq += temp *
                         qqi[k + 1];
                }
              }

              for (u = j + 1; (unsigned int)u < l; u++) {
                /* 1st, l-1 is unpaired */
                if (hc_up[l - 1]) {
                  temp = qm1ui[1][u] * expMLbase[1];
                  if (sc)
                    if (sc->exp_energy_up)
                      temp *= sc->exp_energy_up[l - 1][1];

                  qm1ui[0][u] = temp;
                } else {
                  qm1ui[0][u] = 0.;
                }

                /* 2nd, l-1 is the final position of another motif [p:l-1] */
                for (cnt = 0; cnt < (unsigned int)domains_up->uniq_motif_count; cnt++) {
                  unsigned int size = (unsigned int)domains_up->uniq_motif_size[cnt];
                  if (((unsigned int)u < l - size) &&
                      (hc_up[l - size] >= size)) {
                    temp = qm1ui[size][u] *
                           expMLbase[size] *
                           domains_up->exp_energy_cb(vc,
                                                     l - size,
                                                     l - 1,
                                                     VRNA_UNSTRUCTURED_DOMAIN_MB_LOOP | VRNA_UNSTRUCTURED_DOMAIN_MOTIF,
                                                     domains_up->data);
                    if (sc)
                      if (sc->exp_energy_up)
                        temp *= sc->exp_energy_up[l - size][size];

                    qm1ui[0][u] += temp;
                  }
                }

                /* 3rd, l - 1 pairs with u */
                int ul = my_iindx[u] - (l - 1);
                if (hc[(l - 1) * n + u] & VRNA_CONSTRAINT_CONTEXT_MB_LOOP_ENC) {
                  tt    = vrna_get_ptype(jindx[l - 1] + u, ptype);
                  temp  = qb[ul] *
                          vrna_exp_E_multibranch_stem(tt, S[u - 1], S[l], pf_params);

                  qm1ui[0][u] += temp;
                }

                rqq += qqmj[u - 1] * qm1ui[0][u];
              }

              /* finally, compose contribution */
              exp_motif_ml_left += lqq * rqq;

              /* rotate auxiliary arrays */
              FLT_OR_DBL *tmp = qm1ui[ud_max_size];
              for (cnt = ud_max_size; cnt > 0; cnt--)
                qm1ui[cnt] = qm1ui[cnt - 1];
              qm1ui[0] = tmp;
            }

            /* cleanup memory */
            for (l = 0; l <= ud_max_size; l++)
              free(qm1ui[l]);
            free(qm1ui);

            outside += exp_motif_ml_left *
                       exp_motif_en;

            /* 4.3 Motif is in right-most unpaired stretch of multiloop */
            qmli                = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * n);
            exp_motif_ml_right  = 0.;
            for (k = i - 1; k > 0; k--) {
              FLT_OR_DBL  lqq = 0.;
              FLT_OR_DBL  rqq = 0;

              /* update qmli[k] = qm1[k,i-1] */
              for (qmli[k] = 0., u = k + 1; (unsigned int)u < i; u++) {
                int ku = my_iindx[k] - u;
                /* respect hard constraints */
                if (hc[k * n + u] & VRNA_CONSTRAINT_CONTEXT_MB_LOOP_ENC) {
                  up = (i - 1) - (u + 1) + 1;
                  if (hc_up[u + 1] >= up) {
                    temp = qb[ku] *
                           expMLbase[up];

                    /* add soft constraints */
                    if (sc)
                      if (sc->exp_energy_up)
                        temp *= sc->exp_energy_up[u + 1][up];

                    qmli[k] += temp;

                    /* add contributions of other motifs within [u+1:i-1] */
                    qmli[k] += temp *
                               qqi[u + 1];
                  }
                }
              }

              for (u = k; (unsigned int)u < i; u++)
                lqq += qm[my_iindx[k + 1] - (u - 1)] *
                       qmli[u];

              for (l = j + 1; l <= n; l++) {
                kl = my_iindx[k] - l;
                if (hc[k * n + l] & VRNA_CONSTRAINT_CONTEXT_MB_LOOP) {
                  int jkl;
                  jkl = jindx[l] + k;
                  tt  = rtype[vrna_get_ptype(jkl, ptype)];
                  up  = l - j - 1;
                  if (hc_up[j + 1] >= up) {
                    temp = probs[kl] *
                           vrna_exp_E_multibranch_stem(tt, S[l - 1], S[k + 1], pf_params) *
                           expMLclosing *
                           scale[2] *
                           expMLbase[up];

                    if (sc) {
                      if (sc->exp_energy_bp)
                        temp *= sc->exp_energy_bp[jkl];

                      if (sc->exp_energy_up)
                        temp *= sc->exp_energy_up[j + 1][up];
                    }

                    rqq += temp;

                    /* add contributions of other motifs within [j+1:l-1] */
                    rqq += temp *
                           qqj[l - 1];
                  }
                }
              }
              exp_motif_ml_right += rqq * lqq;
            }
            free(qmli);
            qmli = NULL;

            outside += exp_motif_ml_right *
                       exp_motif_en;
          }
        }

        if (outside > 0.) {
          domains_up->probs_add(vc,
                                i, j,
                                VRNA_UNSTRUCTURED_DOMAIN_MB_LOOP | VRNA_UNSTRUCTURED_DOMAIN_MOTIF,
                                outside,
                                domains_up->data);
        }

        cnt++;
      }
    }

    free(motif_list);
  }

  free(qqi);
  free(qqj);
  free(qqmi);
  free(qqmj);
}


PRIVATE FLT_OR_DBL
numerator_single(vrna_fold_compound_t *vc,
                 int                  i,
                 int                  j)
{
  return 1.;
}


PRIVATE FLT_OR_DBL
numerator_comparative(vrna_fold_compound_t  *vc,
                      int                   i,
                      int                   j)
{
  int     *pscore = vc->pscore;               /* precomputed array of pair types */
  double  kTn     = vc->exp_params->kT / 10.; /* kT in cal/mol  */
  int     *jindx  = vc->jindx;

  return exp(pscore[jindx[j] + i] / kTn);
}


#include "outside_circ.inc"
