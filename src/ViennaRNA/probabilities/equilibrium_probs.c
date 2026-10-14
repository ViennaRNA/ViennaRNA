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

typedef void  (*component_f)(vrna_fold_compound_t *fc,
                             unsigned int         l,
                             helper_arrays_t      eqm_prob_helpers,
                             void                 *data);

typedef void (*post_f)(vrna_fold_compound_t *fc,
                       unsigned int         l,
                       helper_arrays_t      eqm_prob_helpers,
                       void                 *data);

typedef void (*pre_f)(vrna_fold_compound_t  *fc,
                      helper_arrays_t       eqm_prob_helpers,
                      void                  *data);


struct helper_arrays {
  unsigned int  n;
  unsigned int  strands;

  /* overflow counter */
  FLT_OR_DBL    Qmax;
  int           ov;

  /* soft constraints wrappers */
  struct sc_ext_exp_dat     sc_wrapper_ext;
  struct sc_hp_exp_dat      sc_wrapper_hp;
  struct sc_int_exp_dat     sc_wrapper_int;
  struct sc_mb_exp_dat      sc_wrapper_mb;

  FLT_OR_DBL    *probs_ext;
  FLT_OR_DBL    *probs_hp;
  FLT_OR_DBL    *probs_int;
  FLT_OR_DBL    *probs_mb;

  FLT_OR_DBL    *qbh;
  FLT_OR_DBL    *access_ext;
  FLT_OR_DBL    *access_hp;
  FLT_OR_DBL    *access_int;
  FLT_OR_DBL    *access_mb;

  FLT_OR_DBL    *hxp;

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
};


typedef struct {
  unsigned int  n;
  FLT_OR_DBL    *prm_l;
  FLT_OR_DBL    *prm_l1;
  FLT_OR_DBL    *prml;

  unsigned int  ud_max_size;
  FLT_OR_DBL    **pmlu;
  FLT_OR_DBL    *prm_MLbu;
} mb_helpers_t;


/* default options for pf_create_bppm() */
static vrna_eqm_probs_opt_t default_options = {
  .mode = VRNA_EQM_PROBS_BP
};


PRIVATE void
report_n_clear_bpp(vrna_fold_compound_t *fc,
                   unsigned int         l,
                   helper_arrays_t      eqm_prob_helpers,
                   void                 *data);


PRIVATE mb_helpers_t *
get_mb_helpers(unsigned int n);


PRIVATE void
init_mb(vrna_fold_compound_t *fc,
        helper_arrays_t      d,
        void                 *data);


PRIVATE void
free_mb_helpers(void *data);


PRIVATE void
rotate_mb_helpers(vrna_fold_compound_t  *fc,
                  unsigned int          l,
                  helper_arrays_t       dp_helpers,
                  void                  *data);


PRIVATE int
pf_create_bppm(vrna_fold_compound_t   *vc,
               vrna_eqm_probs_f       cb,
               vrna_eqm_probs_data_t  cb_data,
               vrna_eqm_probs_opt_t   options);


PRIVATE void
outside_bp_external(vrna_fold_compound_t  *fc,
                     unsigned int          l,
                     helper_arrays_t       dp_helpers,
                     void                  *data);


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
rotate_dp_helper_arrays_inner(mb_helpers_t *dp_helpers);


PRIVATE void
compute_bpp_internal(vrna_fold_compound_t *fc,
                     unsigned int         l,
                     helper_arrays_t      dp_helpers,
                     void                 *data);


PRIVATE void
compute_bpp_internal_comparative(vrna_fold_compound_t *fc,
                                 unsigned int         l,
                                 helper_arrays_t      dp_helpers,
                                 void                 *data);


PRIVATE void
compute_gquad_prob_internal(vrna_fold_compound_t  *fc,
                            unsigned int          l,
                            helper_arrays_t       dp_helpers,
                            void                  *data);


PRIVATE void
compute_gquad_prob_internal_comparative(vrna_fold_compound_t  *fc,
                                        unsigned int          l,
                                        helper_arrays_t       dp_helpers,
                                        void                  *data);


PRIVATE void
compute_bpp_multibranch(vrna_fold_compound_t  *fc,
                        unsigned int           l,
                        helper_arrays_t       dp_helpers,
                        void                  *data);


PRIVATE void
compute_bpp_multibranch_comparative(vrna_fold_compound_t  *fc,
                                    unsigned int          l,
                                    helper_arrays_t       dp_helpers,
                                    void                  *data);


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


#include "outside_circ.inc"
#include "accessibility.inc"
#include "helices.inc"
#include "multistrand.inc"
#include "unstructured_domains.inc"


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
init_outside(vrna_fold_compound_t *fc,
             helper_arrays_t      eqm_prob_helpers,
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
init_outside_bpp(vrna_fold_compound_t *fc,
                 helper_arrays_t      d,
                 void                 *data)
{
  unsigned int    n = fc->length;

  init_outside(fc, d, data);

  /* reserve memory for probability output */
  d->probs_ext = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * (n + 1));
  d->probs_hp  = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * (n + 1));
  d->probs_int = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * (n + 1));
  d->probs_mb  = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * (n + 1));
}


PRIVATE void
report_n_clear_bpp(vrna_fold_compound_t *fc,
                   unsigned int         l,
                   helper_arrays_t      d,
                   void                 *data)
{

  unsigned int          n         = fc->length;
  int                   *my_iindx = fc->iindx;
  int                   jindx     = fc->jindx[l];
  int                   *pscore   = (fc->type == VRNA_FC_TYPE_COMPARATIVE) ? fc->pscore + jindx : NULL;
  double                kTn       = fc->exp_params->kT / 10.; /* kT in cal/mol  */
  FLT_OR_DBL            qbt, *qb  = fc->exp_matrices->qb - l;
  vrna_eqm_probs_f      cb        = d->cb;
  vrna_eqm_probs_data_t cb_data   = d->cb_data;

  /* loop-type specific output */
  if ((d->cb) &&
      (d->options.mode & VRNA_EQM_PROBS_BP)) {
    /* add inside contribution for base pair probability */
    for (unsigned int k = 1; k < l; k++) {
      qbt = qb[my_iindx[k]];

      if (fc->type == VRNA_FC_TYPE_COMPARATIVE)
        qbt *= exp(-pscore[k] / kTn);

      d->probs_ext[k] *= qbt;
      d->probs_hp[k]  *= qbt;
      d->probs_int[k] *= qbt;
      d->probs_mb[k]  *= qbt;
    }

    if (d->options.mode & VRNA_EQM_PROBS_SPLIT) {
      cb(fc, d->probs_ext, l, 1, l - 1, VRNA_EQM_PROBS_BP | VRNA_LOOP_TYPE_EXT, cb_data);
      cb(fc, d->probs_hp, l, 1, l - 1, VRNA_EQM_PROBS_BP | VRNA_LOOP_TYPE_HP, cb_data);
      cb(fc, d->probs_int, l, 1, l - 1, VRNA_EQM_PROBS_BP | VRNA_LOOP_TYPE_INT, cb_data);
      cb(fc, d->probs_mb, l, 1, l - 1, VRNA_EQM_PROBS_BP | VRNA_LOOP_TYPE_MB, cb_data);
    }

    /* combined probabilities */
    for (unsigned int k = 1; k <= l - 1; k++)
      d->probs_ext[k] += d->probs_hp[k] +
                         d->probs_int[k] +
                         d->probs_mb[k];

    cb(fc, d->probs_ext, l, 1, l - 1, VRNA_EQM_PROBS_BP | VRNA_LOOP_TYPE_ANY, cb_data);
  }

  memset(d->probs_ext, 0, sizeof(FLT_OR_DBL) * (n + 1));
  memset(d->probs_hp, 0, sizeof(FLT_OR_DBL) * (n + 1));
  memset(d->probs_int, 0, sizeof(FLT_OR_DBL) * (n + 1));
  memset(d->probs_mb, 0, sizeof(FLT_OR_DBL) * (n + 1));
}


PRIVATE mb_helpers_t *
get_mb_helpers(unsigned int n)
{
  return (mb_helpers_t *)vrna_alloc(sizeof(mb_helpers_t));
}


PRIVATE void
init_mb(vrna_fold_compound_t *fc,
        helper_arrays_t      d,
        void                 *data)
{
  unsigned int  u, n;
  mb_helpers_t  *dd;

  n           = fc->length;
  dd          = (mb_helpers_t *)data;
  dd->n       = n;
  dd->prm_l   = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * (n + 2));
  dd->prm_l1  = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * (n + 2));
  dd->prml    = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * (n + 2));

  if ((fc->domains_up) &&
      (fc->domains_up->exp_energy_cb)) {
    /* find out maximum size of any unstructured domain */
    for (u = 0; u < (unsigned int)fc->domains_up->uniq_motif_count; u++)
      if (dd->ud_max_size < fc->domains_up->uniq_motif_size[u])
        dd->ud_max_size = fc->domains_up->uniq_motif_size[u];

    dd->pmlu = (FLT_OR_DBL **)vrna_alloc(sizeof(FLT_OR_DBL *) *
                                         (dd->ud_max_size + 1));

    for (u = 0; u <= dd->ud_max_size; u++)
      dd->pmlu[u] = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) *
                                             (n + 2));

    dd->prm_MLbu = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) *
                                            (dd->ud_max_size + 1));

    for (u = 0; u <= dd->ud_max_size; u++)
      dd->prm_MLbu[u] = 0.;
  }
}


PRIVATE void
free_mb_helpers(void *data)
{
  unsigned int  u;
  mb_helpers_t  *dd = (mb_helpers_t *)data;

  free(dd->prm_l);
  free(dd->prm_l1);
  free(dd->prml);

  if (dd->pmlu) {
    for (u = 0; u <= dd->ud_max_size; u++)
      free(dd->pmlu[u]);

    free(dd->pmlu);
  }

  free(dd->prm_MLbu);

  free(dd);
}


PRIVATE void
rotate_mb_helpers(vrna_fold_compound_t  *fc,
                  unsigned int          l,
                  helper_arrays_t       dp_helpers,
                  void                  *data)
{
  unsigned int  u;
  FLT_OR_DBL    *tmp;
  mb_helpers_t  *dd = (mb_helpers_t *)data;

  /* rotate prm_l and prm_l1 arrays */
  tmp                 = dd->prm_l1;
  dd->prm_l1  = dd->prm_l;
  dd->prm_l   = tmp;

  /* rotate pmlu entries required for unstructured domain feature */
  if (dd->pmlu) {
    tmp = dd->pmlu[dd->ud_max_size];

    for (u = dd->ud_max_size; u > 0; u--)
      dd->pmlu[u] = dd->pmlu[u - 1];

    dd->pmlu[0] = tmp;

    for (u = 0; u <= dd->ud_max_size; u++)
      dd->prm_MLbu[u] = 0.;
  }
}


PRIVATE void
register_component(helper_arrays_t      dp_helpers,
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
    vrna_array_append(dp_helpers->component_cbs_init, component_cb_init);
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
  int               n, i, j, l, ij, *pscore, *jindx;
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
    dp_helpers  = get_dp_helper_arrays(vc, cb, cb_data, options);
    with_gquad  = pf_params->model_details.gquad;
    double      kTn = pf_params->kT / 10.;               /* kT in cal/mol  */
    vrna_array(vrna_ep_t) bp_correction = dp_helpers->bp_correction;

    /* call registered init callbacks */
    for (size_t cnt = 0; cnt < vrna_array_size(dp_helpers->component_cbs_init); cnt++)
      if (dp_helpers->component_cbs_init[cnt])
        dp_helpers->component_cbs_init[cnt](vc,
                                            dp_helpers,
                                            dp_helpers->component_cbs_data[cnt]);

    /* Begin computing (outside) probabilities */
    for (l = n; l > 0; l--) {
      /* call registered component callbacks */
      for (size_t cnt = 0;
           cnt < vrna_array_size(dp_helpers->component_cbs);
           cnt++)
        dp_helpers->component_cbs[cnt](vc,
                                       l,
                                       dp_helpers,
                                       dp_helpers->component_cbs_data[cnt]);

      /* call registered post callbacks */
      for (size_t cnt = 0;
           cnt < vrna_array_size(dp_helpers->component_cbs_post);
           cnt++)
        if (dp_helpers->component_cbs_post[cnt])
          dp_helpers->component_cbs_post[cnt](vc,
                                              l,
                                              dp_helpers,
                                              dp_helpers->component_cbs_data[cnt]);
    }

    if (vc->type == VRNA_FC_TYPE_SINGLE) {
      if (with_ud_outside) {
        /*
         *  The above recursions only deal with base pairs, and how they might be
         *  enclosed by other pairs. However, for unstructrued domains, we have
         *  unpaired stretches, and require information about how these are enclosed
         *  by base pairs.
         *
         *  Note, this part should actually be rewritten similar to the accessibility
         *  feature above. Afterall, 'regular' accessibility computations are almost
         *  identical to the unstructured domains feature here...
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


    if (dp_helpers->ov > 0)
      vrna_log_warning("%d overflows occurred while backtracking;\n"
                           "you might try a smaller pf_scale than %g\n",
                           dp_helpers->ov, pf_params->pf_scale);

    /* clean up */
    free_dp_helper_arrays(dp_helpers);
  } /* end if 'check for forward recursion' */
  else {
    vrna_log_warning("bppm calculations have to be done after calling forward recursion");
    return 0;
  }

  return 1;
}


PRIVATE helper_arrays_t
get_dp_helper_arrays(vrna_fold_compound_t   *fc,
                     vrna_eqm_probs_f       cb,
                     vrna_eqm_probs_data_t  cb_data,
                     vrna_eqm_probs_opt_t   options)
{
  unsigned int  n;
  vrna_md_t     *md;
  struct helper_arrays *dp_helpers;

  n           = fc->length;
  md          = &(fc->exp_params->model_details);

  dp_helpers = (struct helper_arrays *)vrna_alloc(sizeof(struct helper_arrays));

  dp_helpers->n       = n;
  dp_helpers->strands = fc->strands;
  dp_helpers->Qmax    = 0.;
  dp_helpers->ov      = 0;

  dp_helpers->probs_ext = NULL;
  dp_helpers->probs_hp  = NULL;
  dp_helpers->probs_int = NULL;
  dp_helpers->probs_mb  = NULL;

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

  init_sc_ext_exp(fc, &(dp_helpers->sc_wrapper_ext));
  init_sc_hp_exp(fc, &(dp_helpers->sc_wrapper_hp));
  init_sc_int_exp(fc, &(dp_helpers->sc_wrapper_int));
  init_sc_mb_exp(fc, &(dp_helpers->sc_wrapper_mb));

  /* Register component callbacks */

  /* 1. unenclosed, i.e. outside probabilities for a base pair that is not enclosed by another one */
  register_component(dp_helpers,
                     (md->circ) ? outside_bp_circ : outside_bp_external,
                     NULL,
                     NULL,
                     (options.mode & VRNA_EQM_PROBS_BP) ? init_outside_bpp : init_outside,
                     (options.mode & VRNA_EQM_PROBS_BP) ? report_n_clear_bpp : NULL);

  /* 2. enclosed, i.e. outside probabilities for a base pair that is part of an internal or multibranch loop */
  register_component(dp_helpers,
                     (fc->type == VRNA_FC_TYPE_SINGLE) ? compute_bpp_internal : compute_bpp_internal_comparative,
                     NULL,
                     NULL,
                     NULL,
                     NULL);

  register_component(dp_helpers,
                     (fc->type == VRNA_FC_TYPE_SINGLE) ? compute_bpp_multibranch : compute_bpp_multibranch_comparative,
                     (void *)get_mb_helpers(n),
                     free_mb_helpers,
                     init_mb,
                     rotate_mb_helpers);

  /* 3. register multistrand extension if necessary */
  if (fc->strands > 1) {
    register_component(dp_helpers,
                       multistrand_contrib,
                       (void *)get_multistrand_helpers(fc->strands),
                       free_multistrand_helpers,
                       init_multistrand,
                       update_multistrand);
  }

  /* register G-quadruplex probability computations */
  if (md->gquad) {
    fc->exp_matrices->p_gq = vrna_smx_csr_FLT_OR_DBL_init(n + 1);
    register_component(dp_helpers,
                       (fc->type == VRNA_FC_TYPE_SINGLE) ? compute_gquad_prob_internal : compute_gquad_prob_internal_comparative,
                       NULL,
                       NULL,
                       NULL,
                       NULL);
  }

  /* Register accessibility computations */
  if (options.mode & VRNA_EQM_PROBS_UP) {
    register_component(dp_helpers,
                       compute_accessibility,
                       NULL,
                       NULL,
                       init_accessibility,
                       report_n_clear_accessibility);
  }

  /* register helix probability computations */
  if (options.mode & VRNA_EQM_PROBS_HX) {
    register_component(dp_helpers,
                       compute_hx_prob,
                       NULL,
                       NULL,
                       init_hxp,
                       NULL);
  }

  if ((fc->type == VRNA_FC_TYPE_SINGLE) &&
      (fc->sc) &&
      (fc->sc->bt))
    vrna_array_init(dp_helpers->bp_correction);

  return dp_helpers;
}


PRIVATE void
free_dp_helper_arrays(helper_arrays_t dp_helpers)
{
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

  free(dp_helpers->probs_ext);
  free(dp_helpers->probs_hp);
  free(dp_helpers->probs_int);
  free(dp_helpers->probs_mb);

  free(dp_helpers->qbh);
  free(dp_helpers->qm2_aux);

  free(dp_helpers->access_ext);
  free(dp_helpers->access_hp);
  free(dp_helpers->access_int);
  free(dp_helpers->access_mb);

  free(dp_helpers->hxp);

  free(dp_helpers);
}


PRIVATE void
rotate_dp_helper_arrays_inner(mb_helpers_t *dp_helpers)
{
  unsigned int u;

  /* rotate prm_MLbu entries required for unstructured domain feature */
  if (dp_helpers->prm_MLbu)
    for (u = dp_helpers->ud_max_size; u > 0; u--)
      dp_helpers->prm_MLbu[u] = dp_helpers->prm_MLbu[u - 1];
}


PRIVATE void
outside_bp_external(vrna_fold_compound_t  *fc,
                     unsigned int         l,
                     helper_arrays_t      dp_helpers,
                     void                 *data)
{
  unsigned char             report_bpp;
  unsigned int              k, n;
  int                       *my_iindx, kl;
  FLT_OR_DBL                *probs, *q1k, *qln, *qb, *probs_ext, q;
  vrna_mx_pf_t              *matrices;
  vrna_hc_t                 *hc;

  FLT_OR_DBL                (*contrib_f)(vrna_fold_compound_t *,
                                         unsigned int,
                                         unsigned int,
                                         helper_arrays_t dp_helpers);

  n           = fc->length;
  my_iindx    = fc->iindx;
  matrices    = fc->exp_matrices;
  qb          = matrices->qb;
  probs       = matrices->probs;
  q1k         = matrices->q1k;
  qln         = matrices->qln;
  hc          = fc->hc;
  probs_ext   = dp_helpers->probs_ext;
  contrib_f   = (fc->type == VRNA_FC_TYPE_SINGLE) ? &contrib_ext_pair : &contrib_ext_pair_comparative;
  report_bpp  = ((dp_helpers->options.mode & VRNA_EQM_PROBS_BP) && (probs_ext)) ? 1 : 0;

  for (k = 1; k < l; k++) {
    kl = my_iindx[k] - l;
    probs[kl] = q = 0.;

    if ((hc->eval_ext(1, n, k, l, VRNA_DECOMP_EXT_STEM_OUTSIDE, hc)) &&
        (qb[kl] > 0.)) {
      q = q1k[k - 1] *
          qln[l + 1] /
          q1k[n];
      q *= contrib_f(fc, k, l, dp_helpers);

      probs[kl] += q;

      if (report_bpp)
        probs_ext[k] += q;
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


PRIVATE void
compute_bpp_internal(vrna_fold_compound_t *fc,
                     unsigned int         l,
                     helper_arrays_t      dp_helpers,
                     void                 *data)
{
  unsigned char         report_bpp;
  char                  *ptype;
  short                 *S1;
  unsigned int          i, j, k, n, u1, u2, with_ud, *hc_up_int, mini, maxj;
  int                   ij, kl, *my_iindx, *jindx, type, type_2, *rtype;
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
  report_bpp    = ((dp_helpers->options.mode & VRNA_EQM_PROBS_BP) && (probs_int)) ? 1 : 0;

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

            if (report_bpp)
              probs_int[k] += tmp2;
          }
        }
      }
    }

    if (probs[kl] > dp_helpers->Qmax) {
      dp_helpers->Qmax = probs[kl];
      if (dp_helpers->Qmax > max_real / 10.)
        vrna_log_warning("P close to overflow: %d %d %g %g\n",
                             k, l, probs[kl], qb[kl]);
    }

    if (probs[kl] >= max_real) {
      dp_helpers->ov++;
      probs[kl] = FLT_MAX;
    }
  }
}


PRIVATE void
compute_bpp_internal_comparative(vrna_fold_compound_t *fc,
                                 unsigned int         l,
                                 helper_arrays_t      dp_helpers,
                                 void                 *data)
{
  unsigned char         report_bpp;
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

  qb          = fc->exp_matrices->qb;
  probs       = fc->exp_matrices->probs;
  scale       = fc->exp_matrices->scale;
  probs_int   = dp_helpers->probs_int;
  report_bpp  = ((dp_helpers->options.mode & VRNA_EQM_PROBS_BP) && (probs_int)) ? 1 : 0;

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
            if (report_bpp)
              probs_int[k] += tmp2;
          }
        }
      }
    }

    if (probs[kl] > dp_helpers->Qmax) {
      dp_helpers->Qmax = probs[kl];
      if (dp_helpers->Qmax > max_real / 10.)
        vrna_log_warning("P close to overflow: %d %d %g %g\n",
                             k, l, probs[kl], qb[kl]);
    }

    if (probs[kl] >= max_real) {
      dp_helpers->ov++;
      probs[kl] = FLT_MAX;
    }
  }

  free(tt);
}


PRIVATE void
compute_bpp_multibranch(vrna_fold_compound_t  *fc,
                        unsigned int          l,
                        helper_arrays_t       dp_helpers,
                        void                  *data)
{
  unsigned char             report_bpp;
  char                      *ptype;
  short                     *S, *S1, s5, s3;
  unsigned int              *sn, tt;
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
  mb_helpers_t              *dd;

  if ((l <= 1) ||
      (l == fc->length))
    return;

  dd            = (mb_helpers_t *)data;
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
  report_bpp  = ((dp_helpers->options.mode & VRNA_EQM_PROBS_BP) && (probs_mb)) ? 1 : 0;

  prm_MLb   = 0.;
  max_real  = (sizeof(FLT_OR_DBL) == sizeof(float)) ? FLT_MAX : DBL_MAX;

  if (sn[l + 1] != sn[l]) {
    /* set prm_l to 0 to get prm_l1 in the next round to be 0 */
    for (i = 0; i <= n; i++)
      dd->prm_l[i] = 0;
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

      dd->prml[i] = prmt;

      /* l+1 is unpaired */
      if (hc->eval_mb(k, l + 1, k, l, VRNA_DECOMP_ML_ML, hc)) {
        ppp = dd->prm_l1[i] *
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
                       dd->pmlu[u][i] *
                       expMLbase[u];

                if (sc_wrapper->red_ml)
                  temp *= sc_wrapper->red_ml(k, l + u, k, l, sc_wrapper);

                ppp += temp;
              }
            }
          }
          dd->pmlu[0][i] = ppp + prmt1;
        }

        dd->prm_l[i] = ppp + prmt1;
      } else {
        /* skip configuration where l+1 is unpaired */
        dd->prm_l[i] = prmt1;

        if (with_ud)
          dd->pmlu[0][i] = prmt1;
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
                       dd->prm_MLbu[u] *
                       expMLbase[u];

                if (sc_wrapper->red_ml)
                  temp *= sc_wrapper->red_ml(i - u + 1, l, i + 1, l, sc_wrapper);

                ppp += temp;
              }
            }
          }
          dd->prm_MLbu[0] = ppp + dd->prml[i];
        }

        prm_MLb = ppp + dd->prml[i];
        /* same as:    prm_MLb = 0;
         * for (i=1; i<=k-1; i++) prm_MLb += prml[i]*expMLbase[k-i-1]; */
      } else {
        /* skip all configurations where i is unpaired */
        prm_MLb = dd->prml[i];

        if (with_ud)
          dd->prm_MLbu[0] = dd->prml[i];
      }

      dd->prml[i] = dd->prml[i] + dd->prm_l[i];

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
            temp += dd->prml[i] *
                    qm[my_iindx[i + 1] - (k - 1)] *
                    sc_wrapper->decomp_ml(i + 1, l, k - 1, k, sc_wrapper);
        } else {
          for (i = 1; i <= k - 2; i++)
            temp += dd->prml[i] *
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
      if (report_bpp)
        probs_mb[k] += temp *
                       scale[2];

      if (probs[kl] > dp_helpers->Qmax) {
        dp_helpers->Qmax = probs[kl];
        if (dp_helpers->Qmax > max_real / 10.)
          vrna_log_warning("P close to overflow: %d %d %g %g\n",
                               k, l, probs[kl], qb[kl]);
      }

      if (probs[kl] >= max_real) {
        dp_helpers->ov++;
        probs[kl] = FLT_MAX;
      }

      /* rotate prm_MLbu entries required for unstructured domain feature */
      rotate_dp_helper_arrays_inner(dd);
    } /* end for (k=..) */
  }
}


PRIVATE void
compute_bpp_multibranch_comparative(vrna_fold_compound_t  *fc,
                                    unsigned int          l,
                                    helper_arrays_t       dp_helpers,
                                    void                  *data)
{
  unsigned char     report_bpp;
  short             **S, **S5, **S3;
  unsigned int      **a2s, s, n_seq, *sn, tt;
  int               i, j, k, n, ii, kl, ij, lj, *my_iindx, *jindx, *pscore, with_gquad;
  FLT_OR_DBL        temp, ppp, prm_MLb, prmt, prmt1, *qb, *probs, *qm, *scale,
                    *expMLbase, expMLclosing, expMLstem, *probs_mb;
  double            max_real, kTn;
  vrna_exp_param_t  *pf_params;
  vrna_md_t         *md;
  vrna_hc_t         *hc;
  vrna_sc_t         **scs;
  vrna_smx_csr(FLT_OR_DBL) *q_gq;
  mb_helpers_t      *dd;

  if ((l <= 1) ||
      (l == fc->length))
    return;

  dd            = (mb_helpers_t *)data;
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

  kTn         = pf_params->kT / 10.;   /* kT in cal/mol  */
  prm_MLb     = 0.;
  max_real    = (sizeof(FLT_OR_DBL) == sizeof(float)) ? FLT_MAX : DBL_MAX;
  probs_mb    = dp_helpers->probs_mb;
  report_bpp  = ((dp_helpers->options.mode & VRNA_EQM_PROBS_BP) && (probs_mb)) ? 1 : 0;

  if (sn[l + 1] != sn[l]) {
    /* set prm_l to 0 to get prm_l1 in the next round to be 0 */
    for (i = 0; i <= n; i++)
      dd->prm_l[i] = 0;
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

      dd->prml[i] = prmt;

      /* l+1 is unpaired */
      if (hc->up_ml[l + 1]) {
        ppp = dd->prm_l1[i] *
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

        dd->prm_l[i] = ppp + prmt1;
      } else {
        /* skip configuration where l+1 is unpaired */
        dd->prm_l[i] = prmt1;
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

        prm_MLb = ppp + dd->prml[i];
        /* same as:    prm_MLb = 0;
         * for (i=1; i<=k-1; i++) prm_MLb += prml[i]*expMLbase[k-i-1]; */
      } else {
        /* skip all configurations where i is unpaired */
        prm_MLb = dd->prml[i];
      }

      dd->prml[i] = dd->prml[i] + dd->prm_l[i];

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
            temp += dd->prml[i] *
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
      if (report_bpp)
        probs_mb[k] += temp *
                       scale[2] *
                       exp(pscore[jindx[l] + k] / kTn);

      if (probs[kl] > dp_helpers->Qmax) {
        dp_helpers->Qmax = probs[kl];
        if (dp_helpers->Qmax > max_real / 10.)
          vrna_log_warning("P close to overflow: %d %d %g %g\n",
                               k, l, probs[kl], qb[kl]);
      }

      if (probs[kl] >= max_real) {
        dp_helpers->ov++;
        probs[kl] = FLT_MAX;
      }

      /* rotate prm_MLbu entries required for unstructured domain feature */
      rotate_dp_helper_arrays_inner(dd);
    } /* end for (k=2..) */
  }
}


PRIVATE void
compute_gquad_prob_internal(vrna_fold_compound_t  *fc,
                            unsigned int          l,
                            helper_arrays_t       dp_helpers,
                            void                  *data)
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
                                        unsigned int          l,
                                        helper_arrays_t       dp_helpers,
                                        void                  *data)
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
