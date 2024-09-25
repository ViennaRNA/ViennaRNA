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
#include "ViennaRNA/probabilities/basepairs.h"

#include "ViennaRNA/constraints/exterior_hc.inc"
#include "ViennaRNA/constraints/hairpin_hc.inc"
#include "ViennaRNA/constraints/internal_hc.inc"
#include "ViennaRNA/constraints/multibranch_hc.inc"

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
typedef struct {
  FLT_OR_DBL    *prm_l;
  FLT_OR_DBL    *prm_l1;
  FLT_OR_DBL    *prml;

  unsigned int  ud_max_size;
  FLT_OR_DBL    **pmlu;
  FLT_OR_DBL    *prm_MLbu;
} helper_arrays;


typedef struct {
  struct hc_ext_def_dat     hc_dat_ext;
  vrna_hc_eval_f hc_eval_ext;

  struct hc_hp_def_dat      hc_dat_hp;
  vrna_hc_eval_f hc_eval_hp;

  struct hc_int_def_dat     hc_dat_int;
  eval_hc                   hc_eval_int;

  struct hc_mb_def_dat      hc_dat_mb;
  vrna_hc_eval_f hc_eval_mb;

  struct sc_ext_exp_dat     sc_wrapper_ext;
  struct sc_hp_exp_dat      sc_wrapper_hp;
  struct sc_int_exp_dat     sc_wrapper_int;
  struct sc_mb_exp_dat      sc_wrapper_mb;
} constraints_helper;


PRIVATE int
pf_create_bppm(vrna_fold_compound_t *vc,
               char                 *structure);


PRIVATE INLINE void
bppm_circ(vrna_fold_compound_t  *fc,
          constraints_helper    *constraints);


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


PRIVATE helper_arrays *
get_ml_helper_arrays(vrna_fold_compound_t *fc);


PRIVATE void
free_ml_helper_arrays(helper_arrays *ml_helpers);


PRIVATE void
rotate_ml_helper_arrays_outer(helper_arrays *ml_helpers);


PRIVATE void
rotate_ml_helper_arrays_inner(helper_arrays *ml_helpers);


PRIVATE constraints_helper *
get_constraints_helper(vrna_fold_compound_t *fc);


PRIVATE void
free_constraints_helper(constraints_helper *helper);


PRIVATE void
compute_bpp_external(vrna_fold_compound_t *fc,
                     constraints_helper   *constraints);


PRIVATE void
compute_bpp_internal(vrna_fold_compound_t *fc,
                     unsigned int         l,
                     vrna_ep_t            **bp_correction,
                     int                  *corr_cnt,
                     int                  *corr_size,
                     FLT_OR_DBL           *Qmax,
                     int                  *ov,
                     constraints_helper   *constraints);


PRIVATE void
compute_bpp_internal_comparative(vrna_fold_compound_t *fc,
                                 unsigned int         l,
                                 vrna_ep_t            **bp_correction,
                                 int                  *corr_cnt,
                                 int                  *corr_size,
                                 FLT_OR_DBL           *Qmax,
                                 int                  *ov,
                                 constraints_helper   *constraints);


PRIVATE void
compute_gquad_prob_internal(vrna_fold_compound_t  *fc,
                            unsigned int          l);


PRIVATE void
compute_gquad_prob_internal_comparative(vrna_fold_compound_t  *fc,
                                        unsigned int          l);


PRIVATE void
compute_bpp_multibranch(vrna_fold_compound_t  *fc,
                        int                   l,
                        helper_arrays         *ml_helpers,
                        FLT_OR_DBL            *Qmax,
                        int                   *ov,
                        constraints_helper    *constraints);


PRIVATE void
compute_bpp_multibranch_comparative(vrna_fold_compound_t  *fc,
                                    int                   l,
                                    helper_arrays         *ml_helpers,
                                    FLT_OR_DBL            *Qmax,
                                    int                   *ov,
                                    constraints_helper    *constraints);


PRIVATE FLT_OR_DBL
contrib_ext_pair(vrna_fold_compound_t *fc,
                 unsigned int         i,
                 unsigned int         j,
                 constraints_helper   *constraints);


PRIVATE FLT_OR_DBL
contrib_ext_pair_comparative(vrna_fold_compound_t *fc,
                             unsigned int         i,
                             unsigned int         j,
                             constraints_helper   *constraints);


PRIVATE void
multistrand_update_Y5(vrna_fold_compound_t  *fc,
                      unsigned int          l,
                      FLT_OR_DBL            *Y5,
                      FLT_OR_DBL            **Y5p,
                      constraints_helper    *constraints);


PRIVATE void
multistrand_update_Y3(vrna_fold_compound_t  *fc,
                      unsigned int          l,
                      FLT_OR_DBL            **Y3,
                      FLT_OR_DBL            **Y3p,
                      constraints_helper    *constraints);


PRIVATE void
multistrand_contrib(vrna_fold_compound_t  *fc,
                    unsigned int          l,
                    FLT_OR_DBL            *Y5,
                    FLT_OR_DBL            **Y3,
                    constraints_helper    *constraints,
                    FLT_OR_DBL            *Qmax,
                    int                   *ov);


/*
 #################################
 # BEGIN OF FUNCTION DEFINITIONS #
 #################################
 */
PUBLIC int
vrna_pairing_probs(vrna_fold_compound_t *vc,
                   char                 *structure)
{
  if (vc)
    return pf_create_bppm(vc, structure);

  return 0;
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
        p *= vrna_exp_E_internal(0,
                           0,
                           vrna_get_ptype(jindx[j] + i, ptype),
                           rtype[vrna_get_ptype(jindx[j - 1] + i + 1, ptype)],
                           0,
                           0,
                           0,
                           0,
                           pf_params) *
             scale[2];

        if (p > cutoff) {
          pl[num].i     = i;
          pl[num].j     = j;
          pl[num].type  = 0;
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


/* calculate base pairing probs */
PRIVATE int
pf_create_bppm(vrna_fold_compound_t *vc,
               char                 *structure)
{
  unsigned int      s;
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
    int         corr_size = 5;
    int         corr_cnt = 0;
    vrna_ep_t   *bp_correction = vrna_alloc(sizeof(vrna_ep_t) * corr_size);
    FLT_OR_DBL  *Y5, **Y5p, **Y3, **Y3p;

    Y5  = NULL;
    Y5p = NULL;
    Y3  = NULL;
    Y3p = NULL;

    helper_arrays       *ml_helpers;
    constraints_helper  *constraints;

    ml_helpers  = get_ml_helper_arrays(vc);
    constraints = get_constraints_helper(vc);

    void                (*compute_bpp_int)(vrna_fold_compound_t *fc,
                                           unsigned int         l,
                                           vrna_ep_t            **bp_correction,
                                           int                  *corr_cnt,
                                           int                  *corr_size,
                                           FLT_OR_DBL           *Qmax,
                                           int                  *ov,
                                           constraints_helper   *constraints);

    void (*compute_bpp_mul)(vrna_fold_compound_t  *fc,
                            int                   l,
                            helper_arrays         *ml_helpers,
                            FLT_OR_DBL            *Qmax,
                            int                   *ov,
                            constraints_helper    *constraints);

    if (vc->type == VRNA_FC_TYPE_SINGLE) {
      compute_bpp_int = &compute_bpp_internal;
      compute_bpp_mul = &compute_bpp_multibranch;
    } else {
      compute_bpp_int = &compute_bpp_internal_comparative;
      compute_bpp_mul = &compute_bpp_multibranch_comparative;
    }

    Qmax = 0;

    if (vc->strands > 1) {
      Y5  = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * vc->strands);
      Y5p = (FLT_OR_DBL **)vrna_alloc(sizeof(FLT_OR_DBL *) * vc->strands);
      for (s = 0; s < vc->strands; s++)
        Y5p[s] = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * (n + 1));

      Y3  = (FLT_OR_DBL **)vrna_alloc(sizeof(FLT_OR_DBL *) * vc->strands);
      Y3p = (FLT_OR_DBL **)vrna_alloc(sizeof(FLT_OR_DBL *) * vc->strands);
      for (s = 0; s < vc->strands; s++) {
        Y3[s]   = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * (n + 1));
        Y3p[s]  = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * (n + 1));
      }
    }

    /* init diagonal entries unable to pair in pr matrix */
    for (i = 1; i <= n; i++)
      probs[my_iindx[i] - i] = 0.;

    /* 1. external loop pairs, i.e. pairs not enclosed by any other pair (or external loop for circular RNAs) */
    if (circular)
      bppm_circ(vc, constraints);
    else
      compute_bpp_external(vc, constraints);

    /* 2. all cases where base pair (k,l) is enclosed by another pair (i,j) */
    l = n;
    compute_bpp_int(vc,
                    l,
                    &bp_correction,
                    &corr_cnt,
                    &corr_size,
                    &Qmax,
                    &ov,
                    constraints);

    for (l = n - 1; l > 1; l--) {
      compute_bpp_int(vc,
                      l,
                      &bp_correction,
                      &corr_cnt,
                      &corr_size,
                      &Qmax,
                      &ov,
                      constraints);

      compute_bpp_mul(vc,
                      l,
                      ml_helpers,
                      &Qmax,
                      &ov,
                      constraints);

      if (vc->strands > 1) {
        multistrand_update_Y5(vc, l, Y5, Y5p, constraints);
        multistrand_update_Y3(vc, l, Y3, Y3p, constraints);
        multistrand_contrib(vc,
                            l,
                            Y5,
                            Y3,
                            constraints,
                            &Qmax,
                            &ov);
      }
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
                  bp_correction[corr_cnt].i   = ptr->i;
                  bp_correction[corr_cnt].j   = ptr->j;
                  bp_correction[corr_cnt++].p = probs[ij] * qhp;
                  if (corr_cnt == corr_size) {
                    corr_size     += 5;
                    bp_correction = vrna_realloc(bp_correction, sizeof(vrna_ep_t) * corr_size);
                  }
                }
              }

              free(aux_bps);
            }
          }

        /*  correct pairing probabilities for auxiliary base pairs from hairpin-, or internal loop motifs
         *  as augmented by the generalized soft constraints feature
         */
        for (i = 0; i < corr_cnt; i++) {
          ij = my_iindx[bp_correction[i].i] - bp_correction[i].j;
          /* printf("correcting pair %d, %d by %f\n", bp_correction[i].i, bp_correction[i].j, bp_correction[i].p); */
          probs[ij] += bp_correction[i].p / qb[ij];
        }
      }
    }

    for (i = 1; i <= n; i++)
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

    if (structure != NULL) {
      char *s = vrna_pairing_tendency(vc);
      memcpy(structure, s, n);
      structure[n] = '\0';
      free(s);
    }

    if (ov > 0)
      vrna_log_warning("%d overflows occurred while backtracking;\n"
                           "you might try a smaller pf_scale than %g\n",
                           ov, pf_params->pf_scale);

    /* clean up */
    free_ml_helper_arrays(ml_helpers);

    free_constraints_helper(constraints);

    free(bp_correction);
    free(Y5);
    if (Y5p)
      for (unsigned int s = 0; s < vc->strands; s++)
        free(Y5p[s]);

    free(Y5p);

    if (Y3)
      for (unsigned int s = 0; s < vc->strands; s++)
        free(Y3[s]);

    free(Y3);

    if (Y3p)
      for (unsigned int s = 0; s < vc->strands; s++)
        free(Y3p[s]);

    free(Y3p);
  } /* end if 'check for forward recursion' */
  else {
    vrna_log_warning("bppm calculations have to be done after calling forward recursion");
    return 0;
  }

  return 1;
}


PRIVATE helper_arrays *
get_ml_helper_arrays(vrna_fold_compound_t *fc)
{
  unsigned int  n, u;
  int           with_ud;
  vrna_ud_t     *domains_up;
  helper_arrays *ml_helpers;

  n           = fc->length;
  domains_up  = fc->domains_up;
  with_ud     = (domains_up && domains_up->exp_energy_cb) ? 1 : 0;

  ml_helpers = (helper_arrays *)vrna_alloc(sizeof(helper_arrays));

  ml_helpers->prm_l   = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * (n + 2));
  ml_helpers->prm_l1  = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * (n + 2));
  ml_helpers->prml    = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * (n + 2));

  ml_helpers->ud_max_size = 0;
  ml_helpers->pmlu        = NULL;
  ml_helpers->prm_MLbu    = NULL;

  if (with_ud) {
    /* find out maximum size of any unstructured domain */
    for (u = 0; u < (unsigned int)domains_up->uniq_motif_count; u++)
      if (ml_helpers->ud_max_size < domains_up->uniq_motif_size[u])
        ml_helpers->ud_max_size = domains_up->uniq_motif_size[u];

    ml_helpers->pmlu = (FLT_OR_DBL **)vrna_alloc(
      sizeof(FLT_OR_DBL *) * (ml_helpers->ud_max_size + 1)); /* maximum motif size */

    for (u = 0; u <= ml_helpers->ud_max_size; u++)
      ml_helpers->pmlu[u] = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * (n + 2));

    ml_helpers->prm_MLbu =
      (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * (ml_helpers->ud_max_size + 1));

    for (u = 0; u <= ml_helpers->ud_max_size; u++)
      ml_helpers->prm_MLbu[u] = 0.;
  }

  return ml_helpers;
}


PRIVATE void
free_ml_helper_arrays(helper_arrays *ml_helpers)
{
  unsigned int u;

  free(ml_helpers->prm_l);
  free(ml_helpers->prm_l1);
  free(ml_helpers->prml);

  if (ml_helpers->pmlu) {
    for (u = 0; u <= ml_helpers->ud_max_size; u++)
      free(ml_helpers->pmlu[u]);
    free(ml_helpers->pmlu);
  }

  free(ml_helpers->prm_MLbu);
  free(ml_helpers);
}


PRIVATE void
rotate_ml_helper_arrays_outer(helper_arrays *ml_helpers)
{
  unsigned int  u;
  FLT_OR_DBL    *tmp;

  /* rotate prm_l and prm_l1 arrays */
  tmp                 = ml_helpers->prm_l1;
  ml_helpers->prm_l1  = ml_helpers->prm_l;
  ml_helpers->prm_l   = tmp;

  /* rotate pmlu entries required for unstructured domain feature */
  if (ml_helpers->pmlu) {
    tmp = ml_helpers->pmlu[ml_helpers->ud_max_size];

    for (u = ml_helpers->ud_max_size; u > 0; u--)
      ml_helpers->pmlu[u] = ml_helpers->pmlu[u - 1];

    ml_helpers->pmlu[0] = tmp;

    for (u = 0; u <= ml_helpers->ud_max_size; u++)
      ml_helpers->prm_MLbu[u] = 0.;
  }
}


PRIVATE void
rotate_ml_helper_arrays_inner(helper_arrays *ml_helpers)
{
  unsigned int u;

  /* rotate prm_MLbu entries required for unstructured domain feature */
  if (ml_helpers->prm_MLbu)
    for (u = ml_helpers->ud_max_size; u > 0; u--)
      ml_helpers->prm_MLbu[u] = ml_helpers->prm_MLbu[u - 1];
}


PRIVATE constraints_helper *
get_constraints_helper(vrna_fold_compound_t *fc)
{
  constraints_helper *helpers;

  helpers = (constraints_helper *)vrna_alloc(sizeof(constraints_helper));

  helpers->hc_eval_ext  = prepare_hc_ext_def(fc, &(helpers->hc_dat_ext));
  helpers->hc_eval_hp   = prepare_hc_hp_def(fc, &(helpers->hc_dat_hp));
  helpers->hc_eval_int  = prepare_hc_int_def(fc, &(helpers->hc_dat_int));
  helpers->hc_eval_mb   = prepare_hc_mb_def(fc, &(helpers->hc_dat_mb));

  init_sc_ext_exp(fc, &(helpers->sc_wrapper_ext));
  init_sc_hp_exp(fc, &(helpers->sc_wrapper_hp));
  init_sc_int_exp(fc, &(helpers->sc_wrapper_int));
  init_sc_mb_exp(fc, &(helpers->sc_wrapper_mb));

  return helpers;
}


PRIVATE void
free_constraints_helper(constraints_helper *helper)
{
  free_sc_ext_exp(&(helper->sc_wrapper_ext));
  free_sc_hp_exp(&(helper->sc_wrapper_hp));
  free_sc_int_exp(&(helper->sc_wrapper_int));
  free_sc_mb_exp(&(helper->sc_wrapper_mb));

  free(helper);
}


PRIVATE void
compute_bpp_external(vrna_fold_compound_t *fc,
                     constraints_helper   *constraints)
{
  unsigned int              i, j, n;
  int                       *my_iindx, ij;
  FLT_OR_DBL                *probs, *q1k, *qln, *qb;
  vrna_mx_pf_t              *matrices;
  struct hc_ext_def_dat     *hc_dat;
  vrna_hc_eval_f evaluate;

  FLT_OR_DBL                (*contrib_f)(vrna_fold_compound_t *,
                                         unsigned int,
                                         unsigned int,
                                         constraints_helper *);

  n         = fc->length;
  my_iindx  = fc->iindx;
  matrices  = fc->exp_matrices;
  qb        = matrices->qb;
  probs     = matrices->probs;
  q1k       = matrices->q1k;
  qln       = matrices->qln;
  hc_dat    = &(constraints->hc_dat_ext);
  evaluate  = constraints->hc_eval_ext;

  contrib_f = (fc->type == VRNA_FC_TYPE_SINGLE) ? &contrib_ext_pair : &contrib_ext_pair_comparative;

  for (i = 1; i <= n; i++) {
    for (j = i + 1; j <= n; j++) {
      ij        = my_iindx[i] - j;
      probs[ij] = 0.;

      if ((evaluate(1, n, i, j, VRNA_DECOMP_EXT_STEM_OUTSIDE, hc_dat)) &&
          (qb[ij] > 0.)) {
        probs[ij] = q1k[i - 1] *
                    qln[j + 1] /
                    q1k[n];
        probs[ij] *= contrib_f(fc, i, j, constraints);
      }
    }
  }
}


PRIVATE FLT_OR_DBL
contrib_ext_pair(vrna_fold_compound_t *fc,
                 unsigned int         i,
                 unsigned int         j,
                 constraints_helper   *constraints)
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
                             constraints_helper   *constraints)
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
                     vrna_ep_t            **bp_correction,
                     int                  *corr_cnt,
                     int                  *corr_size,
                     FLT_OR_DBL           *Qmax,
                     int                  *ov,
                     constraints_helper   *constraints)
{
  unsigned char         type, type_2;
  char                  *ptype;
  short                 *S1;
  unsigned int          i, j, k, n, u1, u2, with_ud, *hc_up_int, mini, maxj;
  int                   ij, kl, *my_iindx, *jindx, *rtype;
  FLT_OR_DBL            temp, tmp2, *qb, *probs, *scale;
  double                max_real;
  vrna_exp_param_t      *pf_params;
  vrna_md_t             *md;
  vrna_hc_t             *hc;
  vrna_sc_t             *sc;
  vrna_ud_t             *domains_up;
  struct hc_int_def_dat *hc_dat_local;
  eval_hc               hc_eval;
  struct sc_int_exp_dat *sc_wrapper_int;

  hc_dat_local    = &(constraints->hc_dat_int);
  hc_eval         = constraints->hc_eval_int;
  sc_wrapper_int  = &(constraints->sc_wrapper_int);

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

          if (hc_eval(i, j, k, l, hc_dat_local)) {
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
                (*bp_correction)[*corr_cnt].i     = ptr->i;
                (*bp_correction)[*corr_cnt].j     = ptr->j;
                (*bp_correction)[(*corr_cnt)++].p = tmp2 * qb[kl];
                if ((*corr_cnt) == (*corr_size)) {
                  (*corr_size)      += 5;
                  (*bp_correction)  = vrna_realloc(*bp_correction,
                                                   sizeof(vrna_ep_t) * (*corr_size));
                }
              }
              free(aux_bps);
            }

            probs[kl] += tmp2;
          }
        }
      }
    }

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
                                 vrna_ep_t            **bp_correction,
                                 int                  *corr_cnt,
                                 int                  *corr_size,
                                 FLT_OR_DBL           *Qmax,
                                 int                  *ov,
                                 constraints_helper   *constraints)
{
  short                 **SS, **S5, **S3;
  unsigned int          i, j, k, n, mini, maxj, type, u1, u2, *tt, s, n_seq, **a2s, *hc_up_int;
  int                   ij, kl, *my_iindx, *jindx, *pscore;
  FLT_OR_DBL            tmp2, *qb, *probs, *scale, psc_exp;
  double                max_real, kTn;
  vrna_exp_param_t      *pf_params;
  vrna_md_t             *md;
  vrna_hc_t             *hc;
  eval_hc               hc_eval;
  struct hc_int_def_dat *hc_dat_local;
  struct sc_int_exp_dat *sc_wrapper_int;

  hc_eval         = constraints->hc_eval_int;
  hc_dat_local    = &(constraints->hc_dat_int);
  sc_wrapper_int  = &(constraints->sc_wrapper_int);

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

          if (hc_eval(i, j, k, l, hc_dat_local)) {
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
          }
        }
      }
    }

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
compute_bpp_multibranch(vrna_fold_compound_t  *fc,
                        int                   l,
                        helper_arrays         *ml_helpers,
                        FLT_OR_DBL            *Qmax,
                        int                   *ov,
                        constraints_helper    *constraints)
{
  unsigned char             tt;
  char                      *ptype;
  short                     *S, *S1, s5, s3;
  unsigned int              *sn;
  int                       cnt, i, j, k, n, u, ii, ij, kl, lj, *my_iindx, *jindx,
                            *rtype, with_gquad, with_ud;
  FLT_OR_DBL                temp, ppp, prm_MLb, prmt, prmt1, *qb, *probs, *qm, *scale,
                            *expMLbase, expMLclosing, expMLstem;
  double                    max_real;
  vrna_exp_param_t          *pf_params;
  vrna_md_t                 *md;
  vrna_ud_t                 *domains_up;
  struct hc_mb_def_dat      *hc_dat;
  vrna_hc_eval_f hc_eval;
  struct sc_mb_exp_dat      *sc_wrapper;
  vrna_smx_csr(FLT_OR_DBL)  *q_gq;


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

  hc_dat      = &(constraints->hc_dat_mb);
  hc_eval     = constraints->hc_eval_mb;
  sc_wrapper  = &(constraints->sc_wrapper_mb);

  prm_MLb   = 0.;
  max_real  = (sizeof(FLT_OR_DBL) == sizeof(float)) ? FLT_MAX : DBL_MAX;

  if (sn[l + 1] != sn[l]) {
    /* set prm_l to 0 to get prm_l1 in the next round to be 0 */
    for (i = 0; i <= n; i++)
      ml_helpers->prm_l[i] = 0;
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
          if (hc_eval(i, j, i + 1, j - 1, VRNA_DECOMP_PAIR_ML, hc_dat)) {
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
        if (hc_eval(i, l + 1, i + 1, l, VRNA_DECOMP_PAIR_ML, hc_dat)) {
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

      ml_helpers->prml[i] = prmt;

      /* l+1 is unpaired */
      if (hc_eval(k, l + 1, k, l, VRNA_DECOMP_ML_ML, hc_dat)) {
        ppp = ml_helpers->prm_l1[i] *
              expMLbase[1];

        if (sc_wrapper->red_ml)
          ppp *= sc_wrapper->red_ml(k, l + 1, k, l, sc_wrapper);

        /* add contributions of MB loops where any unstructured domain starts at l+1 */
        if (with_ud) {
          for (cnt = 0; cnt < domains_up->uniq_motif_count; cnt++) {
            u = domains_up->uniq_motif_size[cnt];
            if (l + u < n) {
              if (hc_eval(k, l + u, k, l, VRNA_DECOMP_ML_ML, hc_dat)) {
                temp = domains_up->exp_energy_cb(fc,
                                                 l + 1,
                                                 l + u,
                                                 VRNA_UNSTRUCTURED_DOMAIN_MB_LOOP | VRNA_UNSTRUCTURED_DOMAIN_MOTIF,
                                                 domains_up->data) *
                       ml_helpers->pmlu[u][i] *
                       expMLbase[u];

                if (sc_wrapper->red_ml)
                  temp *= sc_wrapper->red_ml(k, l + u, k, l, sc_wrapper);

                ppp += temp;
              }
            }
          }
          ml_helpers->pmlu[0][i] = ppp + prmt1;
        }

        ml_helpers->prm_l[i] = ppp + prmt1;
      } else {
        /* skip configuration where l+1 is unpaired */
        ml_helpers->prm_l[i] = prmt1;

        if (with_ud)
          ml_helpers->pmlu[0][i] = prmt1;
      }

      if (hc_eval(i, l, i + 1, l, VRNA_DECOMP_ML_ML, hc_dat)) {
        ppp = prm_MLb *
              expMLbase[1];

        if (sc_wrapper->red_ml)
          ppp *= sc_wrapper->red_ml(i, l, i + 1, l, sc_wrapper);

        if (with_ud) {
          for (cnt = 0; cnt < domains_up->uniq_motif_count; cnt++) {
            u = domains_up->uniq_motif_size[cnt];
            if (1 + u <= i) {
              if (hc_eval(i - u + 1, l, i + 1, l, VRNA_DECOMP_ML_ML, hc_dat)) {
                temp = domains_up->exp_energy_cb(fc,
                                                 i - u + 1,
                                                 i,
                                                 VRNA_UNSTRUCTURED_DOMAIN_MB_LOOP | VRNA_UNSTRUCTURED_DOMAIN_MOTIF,
                                                 domains_up->data) *
                       ml_helpers->prm_MLbu[u] *
                       expMLbase[u];

                if (sc_wrapper->red_ml)
                  temp *= sc_wrapper->red_ml(i - u + 1, l, i + 1, l, sc_wrapper);

                ppp += temp;
              }
            }
          }
          ml_helpers->prm_MLbu[0] = ppp + ml_helpers->prml[i];
        }

        prm_MLb = ppp + ml_helpers->prml[i];
        /* same as:    prm_MLb = 0;
         * for (i=1; i<=k-1; i++) prm_MLb += prml[i]*expMLbase[k-i-1]; */
      } else {
        /* skip all configurations where i is unpaired */
        prm_MLb = ml_helpers->prml[i];

        if (with_ud)
          ml_helpers->prm_MLbu[0] = ml_helpers->prml[i];
      }

      ml_helpers->prml[i] = ml_helpers->prml[i] + ml_helpers->prm_l[i];

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
            temp += ml_helpers->prml[i] *
                    qm[my_iindx[i + 1] - (k - 1)] *
                    sc_wrapper->decomp_ml(i + 1, l, k - 1, k, sc_wrapper);
        } else {
          for (i = 1; i <= k - 2; i++)
            temp += ml_helpers->prml[i] *
                    qm[my_iindx[i + 1] - (k - 1)];
        }
      }

      s5  = ((k > 1) && (sn[k] == sn[k - 1])) ? S1[k - 1] : -1;
      s3  = ((l < n) && (sn[l + 1] == sn[l])) ? S1[l + 1] : -1;

      if ((with_gquad) &&
          (qb[kl] == 0.)) {
        temp *= expMLstem;
      } else if (hc_eval(k, l, k, l, VRNA_DECOMP_ML_STEM, hc_dat)) {
        if (tt == 0)
          tt = 7;

        temp *= vrna_exp_E_multibranch_stem(tt, s5, s3, pf_params);
      }

      if (sc_wrapper->red_stem)
        temp *= sc_wrapper->red_stem(k, l, k, l, sc_wrapper);

      probs[kl] += temp *
                   scale[2];

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
      rotate_ml_helper_arrays_inner(ml_helpers);
    } /* end for (k=..) */
  }

  rotate_ml_helper_arrays_outer(ml_helpers);
}


PRIVATE void
compute_bpp_multibranch_comparative(vrna_fold_compound_t  *fc,
                                    int                   l,
                                    helper_arrays         *ml_helpers,
                                    FLT_OR_DBL            *Qmax,
                                    int                   *ov,
                                    constraints_helper    *constraints)
{
  unsigned char     tt;
  short             **S, **S5, **S3;
  unsigned int      **a2s, s, n_seq, *sn;
  int               i, j, k, n, ii, kl, ij, lj, *my_iindx, *jindx, *pscore, with_gquad;
  FLT_OR_DBL        temp, ppp, prm_MLb, prmt, prmt1, *qb, *probs, *qm, *scale,
                    *expMLbase, expMLclosing, expMLstem;
  double            max_real, kTn;
  vrna_exp_param_t  *pf_params;
  vrna_md_t         *md;
  vrna_hc_t         *hc;
  vrna_sc_t         **scs;
  vrna_smx_csr(FLT_OR_DBL) *q_gq;

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

  if (sn[l + 1] != sn[l]) {
    /* set prm_l to 0 to get prm_l1 in the next round to be 0 */
    for (i = 0; i <= n; i++)
      ml_helpers->prm_l[i] = 0;
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

      ml_helpers->prml[i] = prmt;

      /* l+1 is unpaired */
      if (hc->up_ml[l + 1]) {
        ppp = ml_helpers->prm_l1[i] *
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

        ml_helpers->prm_l[i] = ppp + prmt1;
      } else {
        /* skip configuration where l+1 is unpaired */
        ml_helpers->prm_l[i] = prmt1;
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

        prm_MLb = ppp + ml_helpers->prml[i];
        /* same as:    prm_MLb = 0;
         * for (i=1; i<=k-1; i++) prm_MLb += prml[i]*expMLbase[k-i-1]; */
      } else {
        /* skip all configurations where i is unpaired */
        prm_MLb = ml_helpers->prml[i];
      }

      ml_helpers->prml[i] = ml_helpers->prml[i] + ml_helpers->prm_l[i];

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
            temp += ml_helpers->prml[i] *
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
      rotate_ml_helper_arrays_inner(ml_helpers);
    } /* end for (k=2..) */
  }

  rotate_ml_helper_arrays_outer(ml_helpers);
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

  if (l < n - 3) {
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

  if (l < n - 1) {
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

      for (i = mini; i <= k - 2; i++) {
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
                      FLT_OR_DBL            *Y5,
                      FLT_OR_DBL            **Y5p,
                      constraints_helper    *constraints)
{
  short             *S, *S1;
  unsigned int      i, s, *se, *sn, type, n, end;
  int               *my_iindx;
  FLT_OR_DBL        *probs, *q, *scale, qtmp;
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
  sc_wrapper  = &(constraints->sc_wrapper_ext);
  sc_red_stem = sc_wrapper->red_stem;
  sc_split    = sc_wrapper->split;

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
                      FLT_OR_DBL            **Y3,
                      FLT_OR_DBL            **Y3p,
                      constraints_helper    *constraints)
{
  short             *S, *S1;
  unsigned int      i, j, k, s, n, start, type, *ss, *sn;
  int               *my_iindx;
  FLT_OR_DBL        *q, *probs, *scale, qtmp;
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
  sc_wrapper  = &(constraints->sc_wrapper_ext);
  sc_red_stem = sc_wrapper->red_stem;
  sc_split    = sc_wrapper->split;

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
                    FLT_OR_DBL            *Y5,
                    FLT_OR_DBL            **Y3,
                    constraints_helper    *constraints,
                    FLT_OR_DBL            *Qmax,
                    int                   *ov)
{
  short             *S, *S1, s5, s3;
  unsigned int      k, s, *sn, *se, *ss, end, start, type;
  int               *my_iindx, kl;
  FLT_OR_DBL        *q, *qb, *probs, tmp, qtmp;
  vrna_md_t         *md;
  vrna_exp_param_t  *pf_params;
  struct sc_ext_exp_dat   *sc_wrapper;
  sc_ext_exp_cb           sc_red_stem;
  sc_ext_exp_split        sc_split;

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
  sc_wrapper  = &(constraints->sc_wrapper_ext);
  sc_red_stem = sc_wrapper->red_stem;
  sc_split    = sc_wrapper->split;

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


/* calculate base pairing probs */
PRIVATE INLINE void
bppm_circ(vrna_fold_compound_t  *fc,
          constraints_helper    *constraints)
{
  char                      *ptype;
  unsigned char             *hard_constraints, eval, with_gquad;
  short                     *S, *S1, **SS, **S5, **S3, si, sj;
  unsigned int              s, n_seq, type, rt, *tt, **a2s, turn, u1, u2, u3, us1, us2, us3,
                            imin, imax, jmin, jmax, kmin;
  int                       n, i, j, k, l, ij, *rtype, *my_iindx, *jindx;
  FLT_OR_DBL                tmp, tmp2, tmp3, expMLclosing, *qb, *qm, *probs, *scale, *expMLbase, qo,
                            *qm2_real, q_g, qbt1;
  double                    *expintern;
  vrna_hc_t                 *hc;
  vrna_exp_param_t          *pf_params;
  vrna_mx_pf_t              *matrices;
  vrna_md_t                 *md;
  struct hc_mb_def_dat      *hc_dat_mb;
  vrna_hc_eval_f            hc_eval_mb;

  vrna_smx_csr(FLT_OR_DBL)  *q_gq = fc->exp_matrices->q_gq;

  FLT_OR_DBL                (*numerator_f)(vrna_fold_compound_t *fc,
                                           int                  i,
                                           int                  j);

  n                 = (int)fc->length;
  n_seq             = (fc->type == VRNA_FC_TYPE_SINGLE) ? 1 : fc->n_seq;
  SS                = (fc->type == VRNA_FC_TYPE_SINGLE) ? NULL : fc->S;
  S5                = (fc->type == VRNA_FC_TYPE_SINGLE) ? NULL : fc->S5;
  S3                = (fc->type == VRNA_FC_TYPE_SINGLE) ? NULL : fc->S3;
  a2s               = (fc->type == VRNA_FC_TYPE_SINGLE) ? NULL : fc->a2s;
  pf_params         = fc->exp_params;
  md                = &(pf_params->model_details);
  S                 = (fc->type == VRNA_FC_TYPE_SINGLE) ? fc->sequence_encoding2 : NULL;
  S1                = (fc->type == VRNA_FC_TYPE_SINGLE) ? fc->sequence_encoding : NULL;
  my_iindx          = fc->iindx;
  jindx             = fc->jindx;
  ptype             = (fc->type == VRNA_FC_TYPE_SINGLE) ? fc->ptype : NULL;
  turn              = md->min_loop_size;
  type              = 0;
  rt                = 0;
  tt                = NULL;
  hc                = fc->hc;
  matrices          = fc->exp_matrices;
  qb                = matrices->qb;
  qm                = matrices->qm;
  qm2_real          = matrices->qm2_real;
  probs             = matrices->probs;
  scale             = matrices->scale;
  expMLbase         = matrices->expMLbase;
  qo                = matrices->qo;
  hard_constraints  = hc->mx;
  hc_dat_mb         = &(constraints->hc_dat_mb);
  hc_eval_mb        = constraints->hc_eval_mb;

  expMLclosing  = pf_params->expMLclosing;
  expintern     = &(pf_params->expinternal[0]);
  rtype         = &(pf_params->model_details.rtype[0]);

  struct sc_mb_exp_dat  sc_mb_wrapper;
  struct sc_hp_exp_dat  sc_hp_wrapper;
  struct sc_ext_exp_dat sc_ext_wrapper;
  struct sc_int_exp_dat sc_int_wrapper;
  init_sc_mb_exp(fc, &sc_mb_wrapper);
  init_sc_hp_exp(fc, &sc_hp_wrapper);
  init_sc_ext_exp(fc, &sc_ext_wrapper);
  init_sc_int_exp(fc, &sc_int_wrapper);

  with_gquad    = (unsigned char)md->gquad;

  switch (fc->type) {
    case  VRNA_FC_TYPE_SINGLE:
      numerator_f = numerator_single;
      tt          = NULL;
      break;
    case  VRNA_FC_TYPE_COMPARATIVE:
      numerator_f = numerator_comparative;
      tt          = (unsigned int *)vrna_alloc(sizeof(unsigned int) * n_seq);
      break;
    default:
      numerator_f = NULL;
      break;
  }

  if (with_gquad) {
    matrices->p_gq = vrna_smx_csr_FLT_OR_DBL_init(n + 1);
  }

  /* 1. exterior pair i,j */
  for (i = 1; i <= n; i++) {
    for (j = i; j <= n; j++) {
      ij = my_iindx[i] - j;
      probs[ij] = 0;
    }
  }

  for (i = 1; i <= n; i++) {
    for (j = i + turn + 1; j <= n; j++) {
      ij = my_iindx[i] - j;

      if (qb[ij] > 0.) {
        probs[ij] = numerator_f(fc, i, j) / qo;

        if (fc->type == VRNA_FC_TYPE_SINGLE) {
          type  = vrna_get_ptype_md(S[i], S[j], md);
          rt    = vrna_get_ptype_md(S[j], S[i], md);
        } else {
          for (s = 0; s < n_seq; s++)
            tt[s] = vrna_get_ptype_md(SS[s][j], SS[s][i], md);
        }

        /* 1.1. Exterior Hairpin Contribution */
        tmp2 = vrna_exp_eval_hairpin(fc, j, i, VRNA_EVAL_LOOP_DEFAULT);

        if (hard_constraints[i * n + j] & VRNA_CONSTRAINT_CONTEXT_INT_LOOP) {
          /*
           * 1.2. Exterior Interior Loop Contribution
           * 1.2.1. i,j  delimtis the "left" part of the internal loop
           * (j,i) is "outer pair"
           */
          for (k = 1; k < i - 1; k++) {
            int ln1, ln3, lstart;
            ln1 = k - 1;
            ln3 = n - j;

            if (hc->up_int[j + 1] < (unsigned int)(ln1 + ln3))
              break;

            if ((ln1 + ln3) > MAXLOOP)
              break;

            lstart = (ln1 + ln3) + i - 1 - MAXLOOP;
            if (lstart < k + 1)
              lstart = k + 1;

            for (l = lstart; l < i; l++) {
              int ln2, ln2a, ln1a, type_2;
              ln2 = i - l - 1;

              if ((ln1 + ln2 + ln3) > MAXLOOP)
                continue;

              if (hc->up_int[l + 1] < (unsigned int)ln2)
                continue;

              if (qb[my_iindx[k] - l] == 0.)
                continue;

              eval = (hard_constraints[k * n + l] & VRNA_CONSTRAINT_CONTEXT_INT_LOOP) ? 1 : 0;
              if (hc->f)
                eval = hc->f(k, l, i, j, VRNA_DECOMP_PAIR_IL, hc->data);

              if (eval) {
                tmp = qb[my_iindx[k] - l];

                if (fc->type == VRNA_FC_TYPE_SINGLE) {
                  type_2 = vrna_get_ptype(jindx[l] + k, ptype);

                  tmp *= vrna_exp_E_internal(ln1 + ln3,
                                       ln2,
                                       rt,
                                       rtype[type_2],
                                       S1[j + 1],
                                       S1[i - 1],
                                       S1[k - 1],
                                       S1[l + 1],
                                       pf_params);
                } else {
                  for (s = 0; s < n_seq; s++) {
                    ln2a    = a2s[s][i - 1];
                    ln2a    -= a2s[s][l];
                    ln1a    = a2s[s][n] - a2s[s][j];
                    ln1a    += a2s[s][k - 1];
                    type_2  = vrna_get_ptype_md(SS[s][l], SS[s][k], md);
                    tmp     *= vrna_exp_E_internal(ln1a, ln2a, tt[s], type_2,
                                             S3[s][j],
                                             S5[s][i],
                                             S5[s][k],
                                             S3[s][l],
                                             pf_params);
                  }
                }

                if (sc_int_wrapper.pair_ext)
                  tmp *= sc_int_wrapper.pair_ext(k, l, i, j, &sc_int_wrapper);

                tmp2 += tmp *
                        scale[ln1 + ln2 + ln3];
              }
            }
          }

          /* 1.2.2. i,j  delimtis the "right" part of the internal loop  */
          for (k = j + 1; k < n; k++) {
            int ln1, lstart;
            ln1 = k - j - 1;

            if (hc->up_int[j + 1] < (unsigned int)ln1)
              break;

            if ((ln1 + i - 1) > MAXLOOP)
              break;

            lstart = ln1 + i - 1 + n - MAXLOOP;
            if (lstart < k + 1)
              lstart = k + 1;

            for (l = lstart; l <= n; l++) {
              int ln2, ln3, ln2a, ln1a, type_2;
              ln2 = i - 1;
              ln3 = n - l;

              if ((ln1 + ln2 + ln3) > MAXLOOP)
                continue;

              if (hc->up_int[l + 1] < (unsigned int)(ln2 + ln3))
                continue;

              if (qb[my_iindx[k] - l] == 0.)
                continue;

              eval = (hard_constraints[k * n + l] & VRNA_CONSTRAINT_CONTEXT_INT_LOOP) ? 1 : 0;
              if (hc->f)
                eval = hc->f(i, j, k, l, VRNA_DECOMP_PAIR_IL, hc->data) ? eval : 0;

              if (eval) {
                tmp = qb[my_iindx[k] - l];

                if (fc->type == VRNA_FC_TYPE_SINGLE) {
                  type_2  = vrna_get_ptype(jindx[l] + k, ptype);
                  tmp     *= vrna_exp_E_internal(ln2 + ln3,
                                           ln1,
                                           rtype[type_2],
                                           rt,
                                           S1[l + 1],
                                           S1[k - 1],
                                           S1[i - 1],
                                           S1[j + 1],
                                           pf_params);
                } else {
                  for (s = 0; s < n_seq; s++) {
                    ln1a    = a2s[s][k] - a2s[s][j + 1];
                    ln2a    = a2s[s][i - 1] + a2s[s][n] - a2s[s][l];
                    type_2  = vrna_get_ptype_md(SS[s][l], SS[s][k], md);
                    tmp     *= vrna_exp_E_internal(ln2a, ln1a, type_2, tt[s],
                                             S3[s][l],
                                             S5[s][k],
                                             S5[s][i],
                                             S3[s][j],
                                             pf_params);
                  }
                }

                if (sc_int_wrapper.pair_ext)
                  tmp *= sc_int_wrapper.pair_ext(i, j, k, l, &sc_int_wrapper);

                tmp2 += tmp *
                        scale[ln1 + ln2 + ln3];
              }
            }
          }
        }

        /* 1.3 Exterior multiloop decomposition */
        if (hc_eval_mb(i, j, i - 1, j + 1, VRNA_DECOMP_PAIR_ML_EXT, hc_dat_mb)) {
          FLT_OR_DBL sc_contrib = 1.;

          if (sc_mb_wrapper.pair_ext)
            sc_contrib = sc_mb_wrapper.pair_ext(i, j, &sc_mb_wrapper);

          /* 1.3.1 Middle part                    */
          if ((i > 2) &&
              (j < n - 1)) {
            tmp = qm[my_iindx[1] - i + 1] *
                  qm[my_iindx[j + 1] - n];

            if (fc->type == VRNA_FC_TYPE_SINGLE) {
              tmp *= vrna_exp_E_multibranch_stem(type,
                                  S1[i - 1],
                                  S1[j + 1],
                                  pf_params) *
                     expMLclosing;
            } else {
              for (s = 0; s < n_seq; s++)
                tmp *= vrna_exp_E_multibranch_stem(rtype[tt[s]],
                                    S5[s][i],
                                    S3[s][j],
                                    pf_params);

              tmp *= pow(expMLclosing, n_seq);
            }

            tmp2 += tmp *
                    sc_contrib;
          }

          /* 1.3.2 right-most part  */
          if (hc->up_ml[j + 1] >= (unsigned int)(n - j)) {
            if (i > 1) {
              if ((hc_eval_mb(i, n, i, j, VRNA_DECOMP_ML_ML, hc_dat_mb)) &&
                  (hc_eval_mb(1, j, i - 1, i, VRNA_DECOMP_ML_ML_ML, hc_dat_mb))) {
                tmp = qm2_real[my_iindx[1] - i + 1] *
                      expMLbase[n - j];

                if (fc->type == VRNA_FC_TYPE_SINGLE) {
                  tmp *= vrna_exp_E_multibranch_stem(type,
                                      S1[i - 1],
                                      S1[j + 1],
                                      pf_params) *
                         expMLclosing;
                } else {
                  for (s = 0; s < n_seq; s++)
                    tmp *= vrna_exp_E_multibranch_stem(rtype[tt[s]],
                                        S5[s][i],
                                        S3[s][j],
                                        pf_params);

                  tmp *= pow(expMLclosing, n_seq);
                }

                if (sc_mb_wrapper.red_ml)
                  tmp *= sc_mb_wrapper.red_ml(i, n, i, j, &sc_mb_wrapper);

                if (sc_mb_wrapper.decomp_ml)
                  tmp *= sc_mb_wrapper.decomp_ml(1, j, i - 1, i, &sc_mb_wrapper);

                tmp2 += tmp *
                        sc_contrib;

              }
            }
          }

          /* 1.3.3 left-most part */
          if (hc->up_ml[1] >= (unsigned int)(i - 1)) {
            if (j + 1 < n) {
              if ((hc_eval_mb(1, j, i, j, VRNA_DECOMP_ML_ML, hc_dat_mb)) &&
                  (hc_eval_mb(i, n, j, j + 1, VRNA_DECOMP_ML_ML_ML, hc_dat_mb))) {
                tmp = qm2_real[my_iindx[j + 1] - n] *
                      expMLbase[i - 1];

                if (fc->type == VRNA_FC_TYPE_SINGLE) {
                  tmp *= vrna_exp_E_multibranch_stem(type,
                                      S1[i - 1],
                                      S1[j + 1],
                                      pf_params) *
                         expMLclosing;
                } else {
                  for (s = 0; s < n_seq; s++)
                    tmp *= vrna_exp_E_multibranch_stem(rtype[tt[s]],
                                        S5[s][i],
                                        S3[s][j],
                                        pf_params);

                  tmp *= pow(expMLclosing, n_seq);
                }

                if (sc_mb_wrapper.red_ml)
                  tmp *= sc_mb_wrapper.red_ml(1, j, i, j, &sc_mb_wrapper);

                if (sc_mb_wrapper.decomp_ml)
                  tmp *= sc_mb_wrapper.decomp_ml(i, n, j, j + 1, &sc_mb_wrapper);

                tmp2 += tmp *
                        sc_contrib;
              }
            }
          }

          /* all exterior loop decompositions for pair i,j done  */
        }

        if (with_gquad) {
          unsigned int u1, u2, u3,  us1, us2, us3;
          /* consider all the other possibilities where at least one gquad is enclosed in the same loop as (i,j) */
          /* 1. internal-loop like structure */
          if (hard_constraints[i * n + j] & VRNA_CONSTRAINT_CONTEXT_INT_LOOP) {
            /* 1.1 gquad is 5' of base pair (i,j) */

            FLT_OR_DBL qint = 1.;
            short si, sj;

            switch (fc->type) {
              case VRNA_FC_TYPE_SINGLE:
                sj = S1[j + 1];
                si = S1[i - 1];

                if (md->dangles == 2)
                  qint *= pf_params->expmismatchI[rt][sj][si];

                if (rt > 2)
                  qint *= pf_params->expTermAU;

                break;

              case VRNA_FC_TYPE_COMPARATIVE:
                for (s = 0; s < n_seq; s++) {
                  if (md->dangles == 2)
                    qint *= pf_params->expmismatchI[tt[s]][S3[s][j]][S5[s][i]];

                  if (tt[s] > 2)
                    qint *= pf_params->expTermAU;
                }
                break;
            }

            /* int loop with [k,l] gquad followed by (i,j) base pair */
            if (n - j <= MAXLOOP) {
              unsigned int u1   = n - j;
              unsigned int kmax = 1 + MAXLOOP - u1;
              for (k = 1; (unsigned int)k <= kmax; k++) {
                u2 = k - 1;
                unsigned int lmax = i - 1;
                if ((unsigned int)k + VRNA_GQUAD_MAX_BOX_SIZE - 1 < lmax)
                  lmax = k + VRNA_GQUAD_MAX_BOX_SIZE - 1;

                unsigned int lmin = k + VRNA_GQUAD_MIN_BOX_SIZE - 1;
                if (lmin >= (unsigned int)i)
                  lmin = i;
                else if ((i - lmin) + u1 + u2 > MAXLOOP + 1)
                  lmin = i - (MAXLOOP - u1 - u2) - 1;
                   
                for (l = lmin; (unsigned int)l <= lmax; l++) {
                  if (((u1 + u2 == 0) && (i - l - 1 < 3)) ||
                      ((u1 + u2 < 3) && (i - l - 1 == 0)))
                    continue;

#ifndef VRNA_DISABLE_C11_FEATURES
                  q_g = vrna_smx_csr_get(q_gq, k, l, 0.);
#else
                  q_g = vrna_smx_csr_FLT_OR_DBL_get(q_gq, k, l, 0.);
#endif
                  if (q_g != 0.) {
                    tmp = scale[u1 + u2 + (i - l - 1)];

                    switch (fc->type) {
                      case VRNA_FC_TYPE_SINGLE:
                        tmp *= (FLT_OR_DBL)expintern[u1 + u2 +  (i - l - 1)];
                        break;

                      case VRNA_FC_TYPE_COMPARATIVE:
                        for (s = 0; s < n_seq; s++) {
                          us1   = (k > 1) ? a2s[s][k - 1] - a2s[s][1] : 0;
                          us2   = a2s[s][i - 1] - a2s[s][l];
                          us3   = a2s[s][n] - a2s[s][j];
                          tmp  *= (FLT_OR_DBL)expintern[us1 + us2 + us3];
                        }
                        break;
                    }

                    if (sc_int_wrapper.pair_ext)
                      tmp *= sc_int_wrapper.pair_ext(k, l, i, j, &sc_int_wrapper);

                    /* store outside part for base pair */
                    tmp2 += qint *
                            q_g *
                            tmp;
                    /* store outside part for gquad */
                    probs[my_iindx[k] - l] += qint * tmp * qb[my_iindx[i] - j] / qo;
                  }
                }
              }
            }

            /* 1.2 gquad is 3' of base pair (i,j) */
            u1 = i - 1;
            unsigned int kmax = j + MAXLOOP - u1 + 1;
            if (j + VRNA_GQUAD_MIN_BOX_SIZE <= n)
              kmax = MIN2(kmax, (unsigned int)n - VRNA_GQUAD_MIN_BOX_SIZE + 1);
            else
              kmax = j;

            for (k = j + 1; (unsigned int)k <= kmax; k++) {
              u2 = k - j - 1;
              unsigned int lmin = k + VRNA_GQUAD_MIN_BOX_SIZE - 1;
              if (lmin > (unsigned int)n)
                break;
              if (lmin + MAXLOOP - u1 - u2 < (unsigned int)n)
                lmin = n + u1 + u2 - MAXLOOP;

              unsigned int lmax = k + VRNA_GQUAD_MAX_BOX_SIZE - 1;
              if (lmax > (unsigned int)n)
                lmax = n;

              for (l = lmin; (unsigned int)l <= lmax; l++) {
                u3  = n - l;
                if (((u1 + u3 == 0) && (u2 < 3)) ||
                    ((u1 + u3 < 3) && (u2 == 0)))
                  continue;
              
#ifndef VRNA_DISABLE_C11_FEATURES
                q_g = vrna_smx_csr_get(q_gq, k, l, 0.);
#else
                q_g = vrna_smx_csr_FLT_OR_DBL_get(q_gq, k, l, 0.);
#endif
                if (q_g != 0.) {
                  tmp = scale[u1 + u2 + u3];

                  switch (fc->type) {
                    case VRNA_FC_TYPE_SINGLE:
                      tmp *= (FLT_OR_DBL)expintern[u1 + u2 + u3];
                      break;

                    case VRNA_FC_TYPE_COMPARATIVE:
                      for (s = 0; s < n_seq; s++) {
                        us1   = a2s[s][i - 1] - a2s[s][1];
                        us2   = a2s[s][k - 1] - a2s[s][j];
                        us3   = a2s[s][n] - a2s[s][l];
                        tmp  *= (FLT_OR_DBL)expintern[us1 + us2 + us3];
                      }
                      break;
                  }

                  if (sc_int_wrapper.pair_ext)
                    tmp *= sc_int_wrapper.pair_ext(i, j, k, l, &sc_int_wrapper);

                  /* store outside part for base pair */
                  tmp2 +=  qint *
                          q_g *
                          tmp;
                  /* store outside part for gquad */
                  probs[my_iindx[k] - l] += tmp * qint * qb[my_iindx[i] - j] / qo;
                }
              }
            }
          }

          if (hard_constraints[i * n + j] & VRNA_CONSTRAINT_CONTEXT_MB_LOOP) {
            /* 3. all base pairs (i,j) where G-Quad spans across n,1 boundary */
            /* 3.2 multibranch-loop like configurations */
            FLT_OR_DBL qmb = 1.;
          
            switch (fc->type) {
              case VRNA_FC_TYPE_SINGLE:
                qmb *=  vrna_exp_E_multibranch_stem(type, S1[i - 1], S1[j + 1], pf_params);
                break;
              case VRNA_FC_TYPE_COMPARATIVE:
                for (s = 0; s < n_seq; s++)
                  qmb *= vrna_exp_E_multibranch_stem(rtype[tt[s]],
                                      S5[s][i],
                                      S3[s][j],
                                      pf_params);
                break;
            }

            qmb *= pow(expMLclosing, (double)n_seq);

            FLT_OR_DBL expGQstem = pow(vrna_exp_E_multibranch_stem(0, -1, -1, pf_params), (double)n_seq);

            unsigned int kmin = j + 1;
            if (kmin + VRNA_GQUAD_MAX_BOX_SIZE - 1 <= (unsigned int)n)
              kmin = n - VRNA_GQUAD_MAX_BOX_SIZE + 1;

            for (k = kmin; k <= n; k++) {

              unsigned int lmin = 1;
              unsigned int lmax = VRNA_GQUAD_MAX_BOX_SIZE - 1 - (n - k);
              if (lmax >= (unsigned int)i)
                lmax = i - 1;
              for (l = lmin; (unsigned int)l <= lmax; l++) {
#ifndef VRNA_DISABLE_C11_FEATURES
                if ((q_g = vrna_smx_csr_get(q_gq, k, l, 0.)) != 0.) {
#else
                if ((q_g = vrna_smx_csr_FLT_OR_DBL_get(q_gq, k, l, 0.)) != 0.) {
#endif
                  /* 3.2.1 (i,j) last branch before gquad, ie. [k,l] + ML + (i,j) */
                  if (l + turn < (unsigned int)i) {
                    u1 = k - j - 1;
                    if (hc->up_ml[j + 1] >= u1) {
                      tmp = qmb *
                            qm[my_iindx[l + 1] - i + 1] *
                            expGQstem *
                            expMLbase[u1];

                      if (sc_mb_wrapper.pair_ext)
                        tmp *= sc_mb_wrapper.pair_ext(i, j, &sc_mb_wrapper);

                      /* store outside part for base pair (i,j) */
                      tmp2 += q_g *
                              tmp;
                    }
                  }

                  /* 3.2.2 (i,j) first branch after gquad, ie. [k,l] + (i,j) + ML */
                  if (j + turn < (unsigned int)k) {
                    u2 = i - l - 1;
                    if (hc->up_ml[l + 1] >= u2) {
                      tmp = qmb *
                            qm[my_iindx[j + 1] - k + 1] *
                            expGQstem *
                            expMLbase[u2];
                    
                      if (sc_mb_wrapper.pair_ext)
                        tmp *= sc_mb_wrapper.pair_ext(i, j, &sc_mb_wrapper);

                      /* store outside part for base pair (i,j) */
                      tmp2 += q_g *
                              tmp;
                    }
                  }            

                  /* 3.2.3 (i,j) somewhere in between, ie. [k,l] + ML + (i,j) + ML */
                  if ((l + turn < (unsigned int)i) &&
                      (j + turn < (unsigned int)k)) {
                    tmp = qmb *
                          qm[my_iindx[l + 1] - i + 1] *
                          qm[my_iindx[j + 1] - k + 1] *
                          expGQstem;

                    if (sc_mb_wrapper.pair_ext)
                      tmp *= sc_mb_wrapper.pair_ext(i, j, &sc_mb_wrapper);

                      /* store outside part for base pair (i,j) */
                    tmp2 += q_g *
                            tmp;
                  }
                }
              }
            }
          }
        }

        probs[ij] *= tmp2;
      }
    } /* end for j */
  } /* end for i */

  if (with_gquad) {
    /* first, for all G-Quadruplexes not spanning the n,1 junction */
    for (i = 1; i + VRNA_GQUAD_MIN_BOX_SIZE - 1 <= n; i++) {
      unsigned int maxj = MIN2(i + VRNA_GQUAD_MAX_BOX_SIZE - 1, n);
      for (j = i + VRNA_GQUAD_MIN_BOX_SIZE - 1; (unsigned int)j <= maxj; j++) {
        ij = my_iindx[i] - j;

#ifndef VRNA_DISABLE_C11_FEATURES
        if ((q_g = vrna_smx_csr_get(q_gq, i, j, 0.)) != 0.) {
#else
        if ((q_g = vrna_smx_csr_FLT_OR_DBL_get(q_gq, i, j, 0.)) != 0.) {
#endif
          tmp2 = 0.;

          /* 1. all possible G-quads spanning from i to j */
          
          /* 1.1 hairpin-loop like configurations */
          if ((i - 1) + (n - j) >= 3) {
            eval = (hc->up_ext[1] >= (unsigned int)(i - 1)) ? 1 : 0;
            if (j < n)
              eval = (hc->up_ext[j + 1] >= (unsigned int)(n - j)) ? eval : 0;
            if (hc->f) {
              if (i > 1)
                eval = (hc->f(1, i - 1, 1, i - 1, VRNA_DECOMP_EXT_UP, hc->data)) ? eval : 0;
              if (j < n)
                eval = (hc->f(j + 1, n, j + 1, n, VRNA_DECOMP_EXT_UP, hc->data)) ? eval : 0;
            }
            
            if (eval) {
              tmp = scale[(i - 1) + (n - j)]; 

              if (md->circ_penalty)
                tmp *= pow(vrna_exp_E_hairpin((i - 1) + (n - j), 0, -1, -1, NULL, pf_params), (double)n_seq);

              if (sc_ext_wrapper.red_up)
                tmp *= sc_ext_wrapper.red_up(1, i - 1, &sc_ext_wrapper) *
                       sc_ext_wrapper.red_up(j + 1, n, &sc_ext_wrapper);
 
              tmp2 += tmp;
            }
          }

          /* 1.2 internal-loop like configurations */
          /* 1.2.3 [i,j] [k,l] everything else already handled above */
          u1 = i - 1;
          if ((j + VRNA_GQUAD_MIN_BOX_SIZE <= n) &&
              (u1 <= MAXLOOP) &&
              (hc->up_int[1] >= u1)) {
            unsigned int kmax = j + 1 + MAXLOOP - u1;
            if (kmax + VRNA_GQUAD_MIN_BOX_SIZE - 1 > (unsigned int)n)
              kmax = n - VRNA_GQUAD_MIN_BOX_SIZE + 1;

            for (k = j + 1; (unsigned int)k <= kmax; k++) {
              u2 = k - j - 1;
              if (hc->up_int[j + 1] < u2)
                break;

              unsigned int lmin = k + VRNA_GQUAD_MIN_BOX_SIZE - 1;
              if (lmin + MAXLOOP < n + u1 + u2)
                lmin = n + u1 + u2 - MAXLOOP;

              for (u3 = 0, l = n; (unsigned int)l >= lmin; l--, u3++) {
                if (hc->up_int[l + 1] < u3)
                  break;

                if ((((u1 + u3) == 0) && (u2 < 3)) ||
                    (((u1 + u3) < 3) && (u2 == 0)))
                  continue;

#ifndef VRNA_DISABLE_C11_FEATURES
                FLT_OR_DBL tmp3 = vrna_smx_csr_get(q_gq, k, l, 0.);
#else
                FLT_OR_DBL tmp3 = vrna_smx_csr_FLT_OR_DBL_get(q_gq, k, l, 0.);
#endif
                
                if (tmp3 != 0.) {
                    tmp = scale[u1 + u2 + u3];

                  switch (fc->type) {
                    case VRNA_FC_TYPE_SINGLE:
                      tmp *= (FLT_OR_DBL)expintern[u1 + u2 + u3];
                      break;

                    case VRNA_FC_TYPE_COMPARATIVE:
                      for (s = 0; s < n_seq; s++) {
                        us1   = (i > 1) ? a2s[s][i - 1] - a2s[s][1] : 0;
                        us2   = a2s[s][k - 1] - a2s[s][j];
                        us3   = a2s[s][n] - a2s[s][l];
                        tmp  *= (FLT_OR_DBL)expintern[us1 + us2 + us3];
                      }
                      break;
                  }

                  /* soft constraints for at least unpaired bases in the internal-loop like conformation?! */

                  tmp2 += tmp * tmp3;
                  probs[my_iindx[k] - l] += tmp * q_g / qo; /* add contribution for 2nd gquad here instead of extra block below */
                }
              }
            }
          }

          /* 1.3 multibranch-loop like configurations with [i,j] gquad */
          FLT_OR_DBL expGQstem = pow(vrna_exp_E_multibranch_stem(0, -1, -1, pf_params) *
                                 expMLclosing, (double)n_seq);

          /* 1.3.1 no element spanning n,1 junction */
          /* 1.3.1.1 Middle part                    */
          if ((i > 2) &&
              (j < n - 1)) {
            tmp = qm[my_iindx[1] - i + 1] *
                  qm[my_iindx[j + 1] - n];


            tmp2 += tmp * expGQstem;
          }

          /* 1.3.1,2 right-most part  */
          if (hc->up_ml[j + 1] >= (unsigned int)(n - j)) {
            if (i > 1) {
              if ((hc_eval_mb(i, n, i, j, VRNA_DECOMP_ML_ML, hc_dat_mb)) &&
                  (hc_eval_mb(1, j, i - 1, i, VRNA_DECOMP_ML_ML_ML, hc_dat_mb))) {
                tmp = qm2_real[my_iindx[1] - i + 1] *
                      expMLbase[n - j];

                if (sc_mb_wrapper.red_ml)
                  tmp *= sc_mb_wrapper.red_ml(i, n, i, j, &sc_mb_wrapper);

                if (sc_mb_wrapper.decomp_ml)
                  tmp *= sc_mb_wrapper.decomp_ml(1, j, i - 1, i, &sc_mb_wrapper);

                tmp2 += tmp *
                        expGQstem;
              }
            }
          }

          /* 1.3.1.3 left-most part */
          if (hc->up_ml[1] >= (unsigned int)(i - 1)) {
            if (j + 1 < n) {
              if ((hc_eval_mb(1, j, i, j, VRNA_DECOMP_ML_ML, hc_dat_mb)) &&
                  (hc_eval_mb(i, n, j, j + 1, VRNA_DECOMP_ML_ML_ML, hc_dat_mb))) {
                tmp = qm2_real[my_iindx[j + 1] - n] *
                      expMLbase[i - 1];

                if (sc_mb_wrapper.red_ml)
                  tmp *= sc_mb_wrapper.red_ml(1, j, i, j, &sc_mb_wrapper);

                if (sc_mb_wrapper.decomp_ml)
                  tmp *= sc_mb_wrapper.decomp_ml(i, n, j, j + 1, &sc_mb_wrapper);

                tmp2 += tmp *
                        expGQstem;
              }
            }
          }

          /* 1.3.2 one element spanning n,1 junction */
          FLT_OR_DBL elem = pow(vrna_exp_E_multibranch_stem(0, -1, -1, pf_params), (double)n_seq);

          unsigned int kmin = j + 1;
          if (kmin + VRNA_GQUAD_MAX_BOX_SIZE - 1 <= (unsigned int)n)
            kmin = n - VRNA_GQUAD_MAX_BOX_SIZE + 1;

          for (k = kmin; k <= n; k++) {

            unsigned int lmin = 1;
            unsigned int lmax = VRNA_GQUAD_MAX_BOX_SIZE - 1 - (n - k);
            if (lmax >= (unsigned int)i)
              lmax = i - 1;
            for (l = lmin; (unsigned int)l <= lmax; l++) {
#ifndef VRNA_DISABLE_C11_FEATURES
              if ((q_g = vrna_smx_csr_get(q_gq, k, l, 0.)) != 0.) {
#else
              if ((q_g = vrna_smx_csr_FLT_OR_DBL_get(q_gq, k, l, 0.)) != 0.) {
#endif
                /* 3.2.1 (i,j) last branch before gquad, ie. [k,l] + ML + (i,j) */
                if (l + turn < (unsigned int)i) {
                  u1 = k - j - 1;
                  if (hc->up_ml[j + 1] >= u1) {
                    tmp = expGQstem *
                          qm[my_iindx[l + 1] - i + 1] *
                          elem *
                          expMLbase[u1];

                    if (sc_mb_wrapper.pair_ext)
                      tmp *= sc_mb_wrapper.pair_ext(i, j, &sc_mb_wrapper);

                    /* store outside part for base pair (i,j) */
                    tmp2 += q_g *
                            tmp;
                  }
                }

                /* 3.2.2 (i,j) first branch after gquad, ie. [k,l] + (i,j) + ML */
                if (j + turn < (unsigned int)k) {
                  u2 = i - l - 1;
                  if (hc->up_ml[l + 1] >= u2) {
                    tmp = expGQstem *
                          qm[my_iindx[j + 1] - k + 1] *
                          elem *
                          expMLbase[u2];
                  
                    if (sc_mb_wrapper.pair_ext)
                      tmp *= sc_mb_wrapper.pair_ext(i, j, &sc_mb_wrapper);

                    /* store outside part for base pair (i,j) */
                    tmp2 += q_g *
                            tmp;
                  }
                }            

                /* 3.2.3 (i,j) somewhere in between, ie. [k,l] + ML + (i,j) + ML */
                if ((l + turn < (unsigned int)i) &&
                    (j + turn < (unsigned int)k)) {
                  tmp = expGQstem *
                        qm[my_iindx[l + 1] - i + 1] *
                        qm[my_iindx[j + 1] - k + 1] *
                        elem;

                  if (sc_mb_wrapper.pair_ext)
                    tmp *= sc_mb_wrapper.pair_ext(i, j, &sc_mb_wrapper);

                    /* store outside part for base pair (i,j) */
                  tmp2 += q_g *
                          tmp;
                }
              }
            }
          }

          /* store total contribution */
          probs[ij] += tmp2 / qo;
        }
      }
    }

    /* finally, probabilities for all Gquads spanning the n,1 junction */
    kmin = 2;
    if (kmin + VRNA_GQUAD_MAX_BOX_SIZE - 1 <= (unsigned int)n)
      kmin = n - VRNA_GQUAD_MAX_BOX_SIZE + 1;

    for (k = kmin; k <= n; k++) {
      unsigned int lmin = 1;
      unsigned int lmax = VRNA_GQUAD_MAX_BOX_SIZE - 1 - (n - k);
      if (n - k + lmin < VRNA_GQUAD_MIN_BOX_SIZE)
        lmin = VRNA_GQUAD_MIN_BOX_SIZE + k - n - 1;

      if (lmax >= (unsigned int)n)
        lmax = k - 1;

      for (l = lmin; (unsigned int)l <= lmax; l++) {
#ifndef VRNA_DISABLE_C11_FEATURES
         if ((q_g = vrna_smx_csr_get(q_gq, k, l, 0.)) != 0.) {
#else
         if ((q_g = vrna_smx_csr_FLT_OR_DBL_get(q_gq, k, l, 0.)) != 0.) {
#endif
            tmp2 = 0;

            /* 1. hairpin-loop like case, i.e. single gquad, rest unpaired */
            u1    = k - l - 1;
            eval  = (hc->up_hp[l + 1] >= u1) ? 1 : 0;

            if ((eval) &&
                (u1 >= 3)) {
              qbt1 = scale[u1];

              if (md->circ_penalty)
                qbt1 *= pow(vrna_exp_E_hairpin(u1, 0, -1, -1, NULL, pf_params), (double)n_seq);

              if (sc_hp_wrapper.pair_ext)
                qbt1 *= sc_hp_wrapper.pair_ext(k, l, &sc_hp_wrapper);

              tmp2 += qbt1;
            }

            /* 2. internal-loop like case, i.e. single gquad plus base pair or second gquad */
            imin = l + 1;
            imax = l + MAXLOOP + 1;
            if (imax >= (unsigned int)k)
              imax = k - 1;

            for (i = imin; (unsigned int)i <= imax; i++) {
              u1 = i - l - 1;
              jmax = k - 1;
              jmin = i + turn + 1;
              if ((MAXLOOP <= n) && 
                  (jmin + MAXLOOP < k + u1))
                jmin = k + u1 - MAXLOOP - 1;

              for (j = jmax; (unsigned int)j >= jmin; j--) {
                u2 = k - j - 1;
                if (((u1 == 0) && (u2 < 3)) ||
                    ((u1 < 3) && (u2 == 0)))
                  continue;

                eval = (hc->up_int[l + 1] >= u1) ? 1 : 0;
                eval = (hc->up_int[j + 1] >= u2) ? eval : 0;
                /* 2.1 second element is a base pair */
                if ((eval) &&
                    (hc->mx[n * i + j] & (VRNA_CONSTRAINT_CONTEXT_INT_LOOP | VRNA_CONSTRAINT_CONTEXT_INT_LOOP_ENC)) && 
                    (qb[my_iindx[i] - j] > 0.)) {

                  qbt1 = scale[u1 + u2];
                  FLT_OR_DBL qint = 1.;

                  switch (fc->type) {
                    case VRNA_FC_TYPE_SINGLE:
                      sj = S1[j + 1];
                      si = S1[i - 1];

                      type    = vrna_get_ptype_md(S[j], S[i], md);
                      if (md->dangles == 2)
                        qint *= pf_params->expmismatchI[type][sj][si];

                      if (type > 2)
                        qint *= pf_params->expTermAU;

                      qbt1 *= (FLT_OR_DBL)expintern[u1 + u2];
                      break;

                    case VRNA_FC_TYPE_COMPARATIVE:
                      for (s = 0; s < n_seq; s++) {
                        type = vrna_get_ptype_md(SS[s][j], SS[s][i], md);
                        if (md->dangles == 2)
                          qbt1 *= pf_params->expmismatchI[type][S3[s][j]][S5[s][i]];

                        if (type > 2)
                          qbt1 *= pf_params->expTermAU;

                        us1   = a2s[s][i - 1] - a2s[s][l];
                        us2   = a2s[s][k - 1] - a2s[s][j];
                        qbt1  *= (FLT_OR_DBL)expintern[us1 + us2];
                      }
                      break;
                  }

                  if (sc_int_wrapper.pair_ext)
                    qbt1 *= sc_int_wrapper.pair_ext(i, j, k, l, &sc_int_wrapper);

                  /* store quad probs */
                  tmp2 += qbt1 * qint *
                          qb[my_iindx[i] - j];

                  /* store base pair probs */
                  probs[my_iindx[i] - j] += qbt1 * qint *
                                            q_g / qo;
                }

                /* 2.2 second element is a gquad */
                if ((eval) &&
                    (j - i + 1 >= VRNA_GQUAD_MIN_BOX_SIZE) &&
                    (j - i + 1 <= VRNA_GQUAD_MAX_BOX_SIZE) &&
#ifndef VRNA_DISABLE_C11_FEATURES
                    ((tmp3 = vrna_smx_csr_get(q_gq, i, j, 0.)) != 0.)) {
#else
                    ((tmp3 = vrna_smx_csr_FLT_OR_DBL_get(q_gq, i, j, 0.)) != 0.)) {
#endif
                  qbt1 = scale[u1 + u2];

                  switch (fc->type) {
                    case VRNA_FC_TYPE_SINGLE:
                      qbt1 *= (FLT_OR_DBL)expintern[u1 + u2];
                      break;

                    case VRNA_FC_TYPE_COMPARATIVE:
                      for (s = 0; s < n_seq; s++) {
                        us1   = a2s[s][i - 1] - a2s[s][l];
                        us2   = a2s[s][k - 1] - a2s[s][j];
                        qbt1  *= (FLT_OR_DBL)expintern[us1 + us2];
                      }
                      break;
                  }

                  if (sc_int_wrapper.pair_ext)
                    qbt1 *= sc_int_wrapper.pair_ext(i, j, k, l, &sc_int_wrapper);

                  /* store quad probs */
                  tmp2 += tmp3 *
                          qbt1;

                  /* store other gquad probs */
                  probs[my_iindx[i] - j] += qbt1 *
                                            q_g / qo;
                }
              }
            }

            /* 3. multibranch-loop like case, i.e. gquad plus at least two more stems */
            qbt1 = qm2_real[my_iindx[l + 1] - k + 1] *
                   pow(vrna_exp_E_multibranch_stem(0, -1, -1, pf_params) * expMLclosing, (double)n_seq);

            tmp2 += qbt1;

            /* 4. store gquad outside pf */
#ifndef VRNA_DISABLE_C11_FEATURES
            vrna_smx_csr_insert(matrices->p_gq, k, l, tmp2 * q_g / qo);

#else
            vrna_smx_csr_FLT_OR_DBL_insert(matrices->p_gq, k, l, tmp2 * q_g / qo);
#endif
         }
       }
    }
  }

  free_sc_mb_exp(&sc_mb_wrapper);
  free_sc_hp_exp(&sc_hp_wrapper);
  free_sc_ext_exp(&sc_ext_wrapper);
  free_sc_int_exp(&sc_int_wrapper);
  free(tt);
}
