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
#include "ViennaRNA/params/default.h"
#include "ViennaRNA/fold_vars.h"
#include "ViennaRNA/loops/all.h"
#include "ViennaRNA/gquad.h"
#include "ViennaRNA/constraints/hard.h"
#include "ViennaRNA/constraints/soft.h"
#include "ViennaRNA/mfe.h"
#include "ViennaRNA/eval.h"
#include "ViennaRNA/alphabet.h"
#include "ViennaRNA/part_func.h"
#include "ViennaRNA/equilibrium_probs.h"

#include "ViennaRNA/loops/external_hc.inc"
#include "ViennaRNA/loops/internal_hc.inc"

#include "ViennaRNA/loops/internal_sc_pf.inc"

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
PRIVATE int
pf_create_bppm(vrna_fold_compound_t *vc,
               char                 *structure);


PRIVATE int
pf_co_bppm(vrna_fold_compound_t *vc,
           char                 *structure);


PRIVATE INLINE void
bppm_circ(vrna_fold_compound_t *vc);


PRIVATE INLINE void
bppm_circ_comparative(vrna_fold_compound_t *vc);


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


PRIVATE double
wrap_mean_bp_distance(FLT_OR_DBL  *p,
                      int         length,
                      int         *index,
                      int         turn);


typedef struct {
  FLT_OR_DBL  *prm_l;
  FLT_OR_DBL  *prm_l1;
  FLT_OR_DBL  *prml;

  int         ud_max_size;
  FLT_OR_DBL  **pmlu;
  FLT_OR_DBL  *prm_MLbu;
} helper_arrays;


typedef struct {
  eval_hc               *hc_eval_int;
  struct hc_int_def_dat hc_dat_int;

  struct sc_int_exp_dat sc_wrapper_int;
} constraints_helper;


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
compute_bpp_external(vrna_fold_compound_t *fc);


PRIVATE void
compute_bpp_internal(vrna_fold_compound_t   *fc,
                     int                    l,
                     vrna_ep_t              **bp_correction,
                     int                    *corr_cnt,
                     int                    *corr_size,
                     FLT_OR_DBL             *Qmax,
                     int                    *ov,
                     constraints_helper     *constraints);


PRIVATE void
compute_bpp_internal_comparative(vrna_fold_compound_t   *fc,
                                 int                    l,
                                 vrna_ep_t              **bp_correction,
                                 int                    *corr_cnt,
                                 int                    *corr_size,
                                 FLT_OR_DBL             *Qmax,
                                 int                    *ov,
                                 constraints_helper     *constraints);


PRIVATE void
compute_gquad_prob_internal(vrna_fold_compound_t  *fc,
                            int                   l);


PRIVATE void
compute_gquad_prob_internal_comparative(vrna_fold_compound_t  *fc,
                                        int                   l);


PRIVATE void
compute_bpp_multibranch(vrna_fold_compound_t  *fc,
                        int                   l,
                        helper_arrays         *ml_helpers,
                        FLT_OR_DBL            *Qmax,
                        int                   *ov);


PRIVATE void
compute_bpp_multibranch_comparative(vrna_fold_compound_t  *fc,
                                    int                   l,
                                    helper_arrays         *ml_helpers,
                                    FLT_OR_DBL            *Qmax,
                                    int                   *ov);


PRIVATE FLT_OR_DBL
contrib_ext_pair(vrna_fold_compound_t *fc,
                 unsigned int         i,
                 unsigned int         j);


PRIVATE FLT_OR_DBL
contrib_ext_pair_comparative(vrna_fold_compound_t *fc,
                             unsigned int         i,
                             unsigned int         j);


/*
 #################################
 # BEGIN OF FUNCTION DEFINITIONS #
 #################################
 */
PUBLIC double
vrna_pr_structure(vrna_fold_compound_t  *fc,
                  const char            *structure)
{
  if ((fc) &&
      (fc->exp_params) &&
      (fc->exp_matrices) &&
      (fc->exp_matrices->q)) {
    unsigned int      n;
    double            e, kT, Q, dG, p;
    vrna_exp_param_t  *params = fc->exp_params;
    n = fc->length;

    if (fc->params->model_details.dangles % 2) {
      /* only compute probabilities with dangles = 2 || 0 */
      int dang_bak = fc->params->model_details.dangles;
      fc->params->model_details.dangles = 2;
      e                                 = (double)vrna_eval_structure(fc, structure);
      fc->params->model_details.dangles = dang_bak;
    } else {
      e = (double)vrna_eval_structure(fc, structure);
    }

    kT  = params->kT / 1000.;
    Q   = params->model_details.circ ? fc->exp_matrices->qo : fc->exp_matrices->q[fc->iindx[1] - n];

    dG = (-log(Q) - n * log(params->pf_scale)) * kT;

    if (fc->type == VRNA_FC_TYPE_COMPARATIVE) {
      /* add covariance term */
      e -= vrna_eval_covar_structure(fc, structure);

      /* divide ensemble free energy by number of sequences */
      dG /= fc->n_seq;
    }

    p = exp((dG - e) / kT);

    return p;
  }

  return -1.;
}


PUBLIC double
vrna_pr_energy(vrna_fold_compound_t *fc,
               double               e)
{
  if ((fc) &&
      (fc->exp_params) &&
      (fc->exp_matrices) &&
      (fc->exp_matrices->q)) {
    unsigned int      n;
    double            kT, Q, dG, p;
    vrna_exp_param_t  *params = fc->exp_params;
    n = fc->length;

    kT  = params->kT / 1000.;
    Q   = params->model_details.circ ? fc->exp_matrices->qo : fc->exp_matrices->q[fc->iindx[1] - n];

    dG = (-log(Q) - n * log(params->pf_scale)) * kT;

    if (fc->type == VRNA_FC_TYPE_COMPARATIVE)
      dG /= fc->n_seq;

    p = exp((dG - e) / kT);

    return p;
  }

  return -1.;
}


PUBLIC int
vrna_pairing_probs(vrna_fold_compound_t *vc,
                   char                 *structure)
{
  int ret = 0;

  if (vc) {
    if (vc->strands > 1)
      ret = pf_co_bppm(vc, structure);
    else
      ret = pf_create_bppm(vc, structure);
  }

  return ret;
}


PUBLIC double
vrna_mean_bp_distance_pr(int        length,
                         FLT_OR_DBL *p)
{
  int     *index = vrna_idx_row_wise((unsigned int)length);
  double  d;

  if (!p) {
    vrna_message_warning("vrna_mean_bp_distance_pr: "
                         "p == NULL. "
                         "You need to supply a valid probability matrix");
    return (double)INF / 100.;
  }

  d = wrap_mean_bp_distance(p, length, index, TURN);

  free(index);
  return d;
}


PUBLIC double
vrna_mean_bp_distance(vrna_fold_compound_t *vc)
{
  if (!vc) {
    vrna_message_warning("vrna_mean_bp_distance: run vrna_pf_fold first!");
  } else if (!vc->exp_matrices) {
    vrna_message_warning("vrna_mean_bp_distance: exp_matrices == NULL!");
  } else if (!vc->exp_matrices->probs) {
    vrna_message_warning("vrna_mean_bp_distance: probs==NULL!");
  } else {
    return wrap_mean_bp_distance(vc->exp_matrices->probs,
                                 vc->length,
                                 vc->iindx,
                                 vc->exp_params->model_details.min_loop_size);
  }

  return (double)INF / 100.;
}


PUBLIC double
vrna_ensemble_defect_pt(vrna_fold_compound_t *fc,
                        const short          *pt)
{
  unsigned int  i, j, n;
  int           ii;
  double        ed = -1.;

  if ((fc) &&
      (pt) &&
      (pt[0] == fc->length) &&
      (fc->exp_matrices) &&
      (fc->exp_matrices->probs)) {
    n = fc->length;

    FLT_OR_DBL  *probs  = fc->exp_matrices->probs;
    int         *idx    = fc->iindx;
    ed = 0.;

    for (i = 1; i < n; i++) {
      ii = idx[i];
      double pi;

      /* compute probability to be paired */
      for (pi = 0., j = 1; j < i; j++)
        pi += probs[idx[j] - i];

      for (j = i + 1; j <= n; j++)
        pi += probs[ii - j];

      if (pt[i] == 0)
        ed += pi;
      else if (pt[i] > i)
        ed += 1 - probs[ii - pt[i]];
      else
        ed += 1 - probs[idx[pt[i]] - i];
    }

    ed /= (double)n;

  }

  return ed;
}


PUBLIC double
vrna_ensemble_defect(vrna_fold_compound_t *fc,
                     const char           *structure)
{
  double ed = -1.;
  short *pt = vrna_ptable(structure);

  ed = vrna_ensemble_defect_pt(fc, pt);

  free(pt);

  return ed;
}


PUBLIC double *
vrna_positional_entropy(vrna_fold_compound_t *fc)
{
  double *pos_ent = NULL;

  if ((fc) &&
      (fc->exp_matrices) &&
      (fc->exp_matrices->probs)) {
    unsigned int  i, j, n;
    int           *my_iindx, ii, turn;
    FLT_OR_DBL    *probs;
    double        log2, a, p, *pp;

    log2      = log(2.);
    n         = fc->length;
    my_iindx  = fc->iindx;
    probs     = fc->exp_matrices->probs;
    turn      = fc->exp_params->model_details.min_loop_size;
    pos_ent   = (double *)vrna_alloc(sizeof(double) * (n + 1));
    pp        = (double *)vrna_alloc(sizeof(double) * (n + 1));

    pos_ent[0] = (double)n;

    for (i = 1; i <= n; i++) {
      ii = my_iindx[i];
      for (j = i + turn + 1; j <= n; j++) {
        p = (double)probs[ii - j];
        a = (p > 0.) ? p * log(p) : 0.;
        pos_ent[i] += a;
        pos_ent[j] += a;
        pp[i]      += p;
        pp[j]      += p;
      }
    }

    for (i = 1; i <= n; i++) {
      pos_ent[i] += (pp[i] < 1.) ? (1 - pp[i]) * log(1 - pp[i]) : 0;
      pos_ent[i] /= -log2;
    }

    free(pp);
  }

  return pos_ent;
}


PUBLIC vrna_ep_t *
vrna_stack_prob(vrna_fold_compound_t  *vc,
                double                cutoff)
{
  vrna_ep_t         *pl;
  int               i, j, plsize, turn, length, *index, *jindx, *rtype, num;
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
    turn      = pf_params->model_details.min_loop_size;

    pl = (vrna_ep_t *)vrna_alloc(plsize * sizeof(vrna_ep_t));

    for (i = 1; i < length; i++)
      for (j = i + turn + 3; j <= length; j++) {
        if ((p = probs[index[i] - j]) < cutoff)
          continue;

        if (qb[index[i + 1] - (j - 1)] < FLT_MIN)
          continue;

        p *= qb[index[i + 1] - (j - 1)] / qb[index[i] - j];
        p *=
          exp_E_IntLoop(0, 0, vrna_get_ptype(jindx[j] + i, ptype),
                        rtype[vrna_get_ptype(jindx[j - 1] + i + 1, ptype)],
                        0, 0, 0, 0, pf_params) * scale[2];/* add *scale[u1+u2+2] */
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

  /* compute pair probabilities given that it is a dimer */
  /* AB dimer */
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
        vrna_message_warning(
          "vrna_co_pf_probs: numeric instability detected, probability below zero!");
        lp1->p = 0.;
      }
    }
  }

  return;
}


/* calculate base pairing probs */
PRIVATE int
pf_create_bppm(vrna_fold_compound_t *vc,
               char                 *structure)
{
  int               n, i, j, l, ij, *pscore, *jindx, ov = 0;
  FLT_OR_DBL        Qmax = 0;
  FLT_OR_DBL        *qb, *G, *probs;
  FLT_OR_DBL        *q1k, *qln;

  int               with_gquad;
  vrna_hc_t         *hc;
  vrna_sc_t         *sc;
  int               *my_iindx;
  int               circular, turn, with_ud, with_ud_outside;
  vrna_exp_param_t  *pf_params;
  vrna_mx_pf_t      *matrices;
  vrna_md_t         *md;
  vrna_ud_t         *domains_up;

  n           = vc->length;
  pscore      = (vc->type == VRNA_FC_TYPE_COMPARATIVE) ? vc->pscore : NULL;
  pf_params   = vc->exp_params;
  md          = &(pf_params->model_details);
  my_iindx    = vc->iindx;
  jindx       = vc->jindx;
  circular    = md->circ;
  with_gquad  = md->gquad;
  turn        = md->min_loop_size;

  hc  = vc->hc;
  sc  = vc->sc;

  domains_up  = vc->domains_up;
  matrices    = vc->exp_matrices;

  qb    = matrices->qb;
  G     = matrices->G;
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
      (circular ? matrices->qm2 != NULL : (q1k != NULL && qln != NULL))) {
    with_gquad = pf_params->model_details.gquad;

    double        kTn             = pf_params->kT / 10.; /* kT in cal/mol  */
    int           corr_size       = 5;
    int           corr_cnt        = 0;
    vrna_ep_t     *bp_correction  = vrna_alloc(sizeof(vrna_ep_t) * corr_size);

    helper_arrays       *ml_helpers;
    constraints_helper  *constraints;

    ml_helpers  = get_ml_helper_arrays(vc);
    constraints = get_constraints_helper(vc);

    void          (*compute_bpp_int)(vrna_fold_compound_t *fc,
                                     int                  l,
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
                            int                   *ov);

    if (vc->type == VRNA_FC_TYPE_SINGLE) {
      compute_bpp_int = &compute_bpp_internal;
      compute_bpp_mul = &compute_bpp_multibranch;
    } else {
      compute_bpp_int = &compute_bpp_internal_comparative;
      compute_bpp_mul = &compute_bpp_multibranch_comparative;
    }

    Qmax = 0;

    /* init diagonal entries unable to pair in pr matrix */
    for (i = 1; i <= n; i++)
      for (j = i; j <= MIN2(i + turn, n); j++)
        probs[my_iindx[i] - j] = 0.;

    /* 1. external loop pairs, i.e. pairs not enclosed by any other pair (or external loop for circular RNAs) */
    compute_bpp_external(vc);

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

    for (l = n - 1; l > turn + 1; l--) {
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
                      &ov);
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
          for (j = i + turn + 1; j <= n; j++) {
            ij = my_iindx[i] - j;
            /*  search for possible auxiliary base pairs in hairpin loop motifs to store
             *  the corresponding probability corrections
             */
            if (hc->mx[i * n + j] & VRNA_CONSTRAINT_CONTEXT_HP_LOOP) {
              vrna_basepair_t *ptr, *aux_bps;
              aux_bps = sc->bt(i, j, i, j, VRNA_DECOMP_PAIR_HP, sc->data);
              if (aux_bps) {
                FLT_OR_DBL qhp = vrna_exp_E_hp_loop(vc, i, j);
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

        /*  correct pairing probabilities for auxiliary base pairs from hairpin-, or interior loop motifs
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
      for (j = i + turn + 1; j <= n; j++) {
        ij = my_iindx[i] - j;

        if (with_gquad) {
          if (qb[ij] > 0.) {
            probs[ij] *= qb[ij];
            if (vc->type == VRNA_FC_TYPE_COMPARATIVE)
              probs[ij] *= exp(-pscore[jindx[j] + i] / kTn);
          } else if (G[ij] > 0.) {
            probs[ij] +=  q1k[i - 1] *
                          G[ij] *
                          qln[j + 1] /
                          q1k[n];
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
      char *s = vrna_db_from_probs(probs, (unsigned int)n);
      memcpy(structure, s, n);
      structure[n] = '\0';
      free(s);
    }

    if (ov > 0)
      vrna_message_warning("%d overflows occurred while backtracking;\n"
                           "you might try a smaller pf_scale than %g\n",
                           ov, pf_params->pf_scale);

    /* clean up */
    free_ml_helper_arrays(ml_helpers);

    free_constraints_helper(constraints);

    free(bp_correction);
  } /* end if 'check for forward recursion' */
  else {
    vrna_message_warning("bppm calculations have to be done after calling forward recursion");
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
    for (u = 0; u < domains_up->uniq_motif_count; u++)
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
get_constraints_helper(vrna_fold_compound_t *fc){
  constraints_helper *helpers;

  helpers = (constraints_helper *)vrna_alloc(sizeof(constraints_helper));

  helpers->hc_eval_int = prepare_hc_int_def(fc, &(helpers->hc_dat_int));
  init_sc_int_exp(fc, &(helpers->sc_wrapper_int));

  return helpers;
}


PRIVATE void
free_constraints_helper(constraints_helper *helper)
{

  free_sc_int_exp(&(helper->sc_wrapper_int));

  free(helper);
}


PRIVATE void
compute_bpp_external(vrna_fold_compound_t *fc)
{
  unsigned int      i, j, n, turn;
  int               circular, *my_iindx, ij;
  FLT_OR_DBL        *probs, *q1k, *qln, *qb;
  vrna_exp_param_t  *pf_params;
  vrna_md_t         *md;
  vrna_mx_pf_t      *matrices;

  FLT_OR_DBL        (*contrib_f)(vrna_fold_compound_t *,
                                 unsigned int,
                                 unsigned int);

  n         = fc->length;
  my_iindx  = fc->iindx;
  pf_params = fc->exp_params;
  md        = &(pf_params->model_details);
  matrices  = fc->exp_matrices;
  qb        = matrices->qb;
  probs     = matrices->probs;
  q1k       = matrices->q1k;
  qln       = matrices->qln;
  circular  = md->circ;
  turn      = md->min_loop_size;

  contrib_f = (fc->type == VRNA_FC_TYPE_SINGLE) ? &contrib_ext_pair : &contrib_ext_pair_comparative;

  if (circular) {
    (fc->type == VRNA_FC_TYPE_SINGLE) ? bppm_circ(fc) : bppm_circ_comparative(fc);
  } else {
    struct hc_ext_def_dat     hc_dat_local;
    vrna_callback_hc_evaluate *evaluate = prepare_hc_ext_def(fc, &hc_dat_local);

    for (i = 1; i <= n; i++) {
      for (j = i + turn + 1; j <= n; j++) {
        ij        = my_iindx[i] - j;
        probs[ij] = 0.;

        if ((evaluate(1, n, i, j, VRNA_DECOMP_EXT_STEM_OUTSIDE, &hc_dat_local)) &&
            (qb[ij] > 0.)) {
          probs[ij] = q1k[i - 1] *
                      qln[j + 1] /
                      q1k[n];
          probs[ij] *= contrib_f(fc, i, j);
        }
      }
    }
  }
}


PRIVATE FLT_OR_DBL
contrib_ext_pair(vrna_fold_compound_t *fc,
                 unsigned int         i,
                 unsigned int         j)
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

  contribution = vrna_exp_E_ext_stem(type, s5, s3, pf_params);

  if ((sc) &&
      (sc->exp_f))
    contribution *= sc->exp_f(1, n, i, j, VRNA_DECOMP_EXT_STEM_OUTSIDE, sc->data);

  return contribution;
}


PRIVATE FLT_OR_DBL
contrib_ext_pair_comparative(vrna_fold_compound_t *fc,
                             unsigned int         i,
                             unsigned int         j)
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

    contribution *= vrna_exp_E_ext_stem(type, s5, s3, pf_params);
  }

  if (scs) {
    for (s = 0; s < n_seq; s++)
      if (scs[s]->exp_f)
        contribution *= scs[s]->exp_f(1, n, i, j, VRNA_DECOMP_EXT_STEM_OUTSIDE, scs[s]->data);
  }

  return contribution;
}


PRIVATE void
compute_bpp_internal(vrna_fold_compound_t   *fc,
                     int                    l,
                     vrna_ep_t              **bp_correction,
                     int                    *corr_cnt,
                     int                    *corr_size,
                     FLT_OR_DBL             *Qmax,
                     int                    *ov,
                     constraints_helper     *constraints)
{
  unsigned char     type, type_2;
  char              *ptype;
  short             *S1;
  unsigned int      *sn;
  int               i, j, k, n, ij, kl, u1, u2, *my_iindx, *jindx, *rtype,
                    turn, with_ud, *hc_up_int;
  FLT_OR_DBL        temp, tmp2, *qb, *probs, *scale;
  double            max_real;
  vrna_exp_param_t  *pf_params;
  vrna_md_t         *md;
  vrna_hc_t         *hc;
  vrna_sc_t         *sc;
  vrna_ud_t         *domains_up;
  eval_hc                *hc_eval;
  struct hc_int_def_dat  *hc_dat_local;
  struct sc_int_exp_dat  *sc_wrapper_int;

  hc_eval         = constraints->hc_eval_int;
  hc_dat_local    = &(constraints->hc_dat_int);
  sc_wrapper_int  = &(constraints->sc_wrapper_int);

  n           = (int)fc->length;
  ptype       = fc->ptype;
  S1          = fc->sequence_encoding;
  sn          = fc->strand_number;
  my_iindx    = fc->iindx;
  jindx       = fc->jindx;
  pf_params   = fc->exp_params;
  md          = &(pf_params->model_details);
  turn        = md->min_loop_size;
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
  for (k = 1; k < l - turn; k++) {
    kl = my_iindx[k] - l;

    if (qb[kl] == 0.)
      continue;

    if (hc->mx[l * n + k] & VRNA_CONSTRAINT_CONTEXT_INT_LOOP_ENC) {
      type_2 = rtype[vrna_get_ptype(jindx[l] + k, ptype)];

      for (i = MAX2(1, k - MAXLOOP - 1); i <= k - 1; i++) {
        u1 = k - i - 1;
        if (hc_up_int[i + 1] < u1)
          continue;

        int max_j = l + 1 + MAXLOOP - u1;

        if (max_j > n)
          max_j = n;

        if (max_j > l + 1 + hc_up_int[l + 1])
          max_j = l + 1 + hc_up_int[l + 1];

        u2 = 0;

        for (j = l + 1; j <= max_j; j++, u2++) {
          ij = my_iindx[i] - j;

          if (probs[ij] == 0.)
            continue;

#if 0
          if (hc_eval(i, j, k, l, hc_dat_local)) {
            {
#else
          if (hc->mx[i * n + j] & VRNA_CONSTRAINT_CONTEXT_INT_LOOP) {
            if ((sn[k] == sn[i]) &&
                (sn[j] == sn[l])) {
#endif
              int jij = jindx[j] + i;
              type = vrna_get_ptype(jij, ptype);
              tmp2 = probs[ij]
                     * scale[u1 + u2 + 2]
                     * exp_E_IntLoop(u1,
                                     u2,
                                     type,
                                     type_2,
                                     S1[i + 1],
                                     S1[j - 1],
                                     S1[k - 1],
                                     S1[l + 1],
                                     pf_params);

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
                /* store probability correction for auxiliary pairs in interior loop motif */
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
    }

    if (probs[kl] > (*Qmax)) {
      (*Qmax) = probs[kl];
      if ((*Qmax) > max_real / 10.)
        vrna_message_warning("P close to overflow: %d %d %g %g\n",
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
compute_bpp_internal_comparative(vrna_fold_compound_t   *fc,
                                 int                    l,
                                 vrna_ep_t              **bp_correction,
                                 int                    *corr_cnt,
                                 int                    *corr_size,
                                 FLT_OR_DBL             *Qmax,
                                 int                    *ov,
                                 constraints_helper     *constraints)
{
  short             **SS, **S5, **S3;
  unsigned int      type, *tt, s, n_seq, **a2s, *sn;
  int               i, j, k, n, ij, kl, u1, u2, *my_iindx, *jindx, turn, *pscore, *hc_up_int;
  FLT_OR_DBL        tmp2, *qb, *probs, *scale, psc_exp;
  double            max_real, kTn;
  vrna_exp_param_t  *pf_params;
  vrna_md_t         *md;
  vrna_hc_t         *hc;
  vrna_sc_t         **scs;
  eval_hc                *hc_eval;
  struct hc_int_def_dat  *hc_dat_local;
  struct sc_int_exp_dat  *sc_wrapper_int;

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
  sn        = fc->strand_number;
  my_iindx  = fc->iindx;
  jindx     = fc->jindx;
  pf_params = fc->exp_params;
  md        = &(pf_params->model_details);
  turn      = md->min_loop_size;
  hc        = fc->hc;
  scs       = fc->scs;
  hc_up_int = hc->up_int;

  qb    = fc->exp_matrices->qb;
  probs = fc->exp_matrices->probs;
  scale = fc->exp_matrices->scale;

  kTn       = pf_params->kT / 10.;   /* kT in cal/mol  */
  max_real  = (sizeof(FLT_OR_DBL) == sizeof(float)) ? FLT_MAX : DBL_MAX;
  tt        = (unsigned int *)vrna_alloc(sizeof(unsigned int) * n_seq);

  /* 2. bonding k,l as substem of 2:loop enclosed by i,j */
  for (k = 1; k < l - turn; k++) {
    kl = my_iindx[k] - l;

    if (qb[kl] == 0.)
      continue;

    if (hc->mx[l * n + k] & VRNA_CONSTRAINT_CONTEXT_INT_LOOP_ENC) {
      psc_exp = exp(pscore[jindx[l] + k] / kTn);

      for (s = 0; s < n_seq; s++)
        tt[s] = vrna_get_ptype_md(SS[s][l], SS[s][k], md);

      for (i = MAX2(1, k - MAXLOOP - 1); i <= k - 1; i++) {
        u1 = k - i - 1;
        if (hc_up_int[i + 1] < u1)
          continue;

        for (j = l + 1; j <= MIN2(l + MAXLOOP - k + i + 2, n); j++) {
          ij = my_iindx[i] - j;

          if (probs[ij] == 0.)
            continue;

          u2 = j - l - 1;

          if (hc_up_int[l + 1] < u2)
            break;

          if (hc->mx[i * n + j] & VRNA_CONSTRAINT_CONTEXT_INT_LOOP) {
            if ((sn[k] == sn[i]) &&
                (sn[j] == sn[l])) {
              tmp2 = probs[ij] *
                     scale[u1 + u2 + 2] *
                     psc_exp;

              for (s = 0; s < n_seq; s++) {
                int u1_loc  = a2s[s][k - 1] - a2s[s][i];
                int u2_loc  = a2s[s][j - 1] - a2s[s][l];
                type        = vrna_get_ptype_md(SS[s][i], SS[s][j], md);
                tmp2       *= exp_E_IntLoop(u1_loc,
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
    }

    if (probs[kl] > (*Qmax)) {
      (*Qmax) = probs[kl];
      if ((*Qmax) > max_real / 10.)
        vrna_message_warning("P close to overflow: %d %d %g %g\n",
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
                        int                   *ov)
{
  unsigned char     tt;
  char              *ptype;
  short             *S, *S1, s5, s3;
  unsigned int      *sn;
  int               cnt, i, j, k, n, u, ii, ij, kl, lj, turn, *my_iindx, *jindx,
                    *rtype, with_gquad, with_ud;
  FLT_OR_DBL        temp, ppp, prm_MLb, prmt, prmt1, *qb, *probs, *qm, *G, *scale,
                    *expMLbase, expMLclosing, expMLstem;
  double            max_real;
  vrna_exp_param_t  *pf_params;
  vrna_md_t         *md;
  vrna_hc_t         *hc;
  vrna_sc_t         *sc;
  vrna_ud_t         *domains_up;

  n             = (int)fc->length;
  S             = fc->sequence_encoding2;
  S1            = fc->sequence_encoding;
  sn            = fc->strand_number;
  my_iindx      = fc->iindx;
  jindx         = fc->jindx;
  pf_params     = fc->exp_params;
  md            = &(pf_params->model_details);
  turn          = md->min_loop_size;
  rtype         = &(md->rtype[0]);
  ptype         = fc->ptype;
  qb            = fc->exp_matrices->qb;
  qm            = fc->exp_matrices->qm;
  G             = fc->exp_matrices->G;
  probs         = fc->exp_matrices->probs;
  scale         = fc->exp_matrices->scale;
  expMLbase     = fc->exp_matrices->expMLbase;
  expMLclosing  = pf_params->expMLclosing;
  hc            = fc->hc;
  sc            = fc->sc;
  domains_up    = fc->domains_up;
  with_ud       = (domains_up && domains_up->exp_energy_cb) ? 1 : 0;
  with_gquad    = md->gquad;
  expMLstem     = (with_gquad) ? exp_E_MLstem(0, -1, -1, pf_params) : 0;

  prm_MLb   = 0.;
  max_real  = (sizeof(FLT_OR_DBL) == sizeof(float)) ? FLT_MAX : DBL_MAX;

  if (sn[l + 1] != sn[l]) {
    /* set prm_l to 0 to get prm_l1 in the next round to be 0 */
    for (i = 0; i <= n; i++)
      ml_helpers->prm_l[i] = 0;
  } else {
    for (k = 2; k < l - turn; k++) {
      kl    = my_iindx[k] - l;
      i     = k - 1;
      prmt  = prmt1 = 0.0;

      ij  = my_iindx[i] - (l + 2);
      lj  = my_iindx[l + 1] - (l + 1);
      s3  = S1[i + 1];
      if (sn[k] == sn[i]) {
        for (j = l + 2; j <= n; j++, ij--, lj--) {
          if ((hc->mx[i * n + j] & VRNA_CONSTRAINT_CONTEXT_MB_LOOP) &&
              (sn[j] == sn[j - 1])) {
            tt = vrna_get_ptype_md(S[j], S[i], md);

            /* which decomposition is covered here? =>
             * i + 1 = k < l < j:
             * (i,j)       -> enclosing pair
             * (k, l)      -> enclosed pair
             * (l+1, j-1)  -> multiloop part with at least one stem
             * a.k.a. (k,l) is left-most stem in multiloop closed by (k-1, j)
             */
            ppp = probs[ij]
                  * exp_E_MLstem(tt, S1[j - 1], s3, pf_params)
                  * qm[lj];

            if (sc) {
              if (sc->exp_energy_bp)
                ppp *= sc->exp_energy_bp[jindx[j] + i];

              /*
               *        if(sc->exp_f)
               *          ppp *= sc->exp_f(i, j, l+1, j-1, , sc->data);
               */
            }

            prmt += ppp;
          }
        }

        ii  = my_iindx[i];  /* ii-j=[i,j]     */
        tt  = vrna_get_ptype(jindx[l + 1] + i, ptype);
        tt  = rtype[tt];
        if (hc->mx[(l + 1) * n + i] & VRNA_CONSTRAINT_CONTEXT_MB_LOOP) {
          prmt1 = probs[ii - (l + 1)]
                  *expMLclosing
                  *exp_E_MLstem(tt,
                                S1[l],
                                S1[i + 1],
                                pf_params);


          if (sc) {
            /* which decompositions are covered here? => (i, l+1) -> enclosing pair */
            if (sc->exp_energy_bp)
              prmt1 *= sc->exp_energy_bp[jindx[l + 1] + i];

            /*
             *      if(sc->exp_f)
             *        prmt1 *= sc->exp_f(i, l+1, k, l, , sc->data);
             */
          }
        }
      }

      prmt *= expMLclosing;

      ml_helpers->prml[i] = prmt;

      /* l+1 is unpaired */
      if (hc->up_ml[l + 1]) {
        ppp = ml_helpers->prm_l1[i] * expMLbase[1];
        if (sc) {
          if (sc->exp_energy_up)
            ppp *= sc->exp_energy_up[l + 1][1];

          /*
           *      if(sc_exp_f)
           *        ppp *= sc->exp_f(, sc->data);
           */
        }

        /* add contributions of MB loops where any unstructured domain starts at l+1 */
        if (with_ud) {
          for (cnt = 0; cnt < domains_up->uniq_motif_count; cnt++) {
            u = domains_up->uniq_motif_size[cnt];
            if (hc->up_ml[l + 1] >= u) {
              if (l + u < n) {
                temp = domains_up->exp_energy_cb(fc,
                                                 l + 1,
                                                 l + u,
                                                 VRNA_UNSTRUCTURED_DOMAIN_MB_LOOP | VRNA_UNSTRUCTURED_DOMAIN_MOTIF,
                                                 domains_up->data)
                       * ml_helpers->pmlu[u][i]
                       * expMLbase[u];

                if (sc)
                  if (sc->exp_energy_up)
                    temp *= sc->exp_energy_up[l + 1][u];

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

      /* i is unpaired */
      if (hc->up_ml[i]) {
        ppp = prm_MLb * expMLbase[1];
        if (sc) {
          if (sc->exp_energy_up)
            ppp *= sc->exp_energy_up[i][1];

          /*
           *      if(sc->exp_f)
           *        ppp *= sc->exp_f(, sc->data);
           */
        }

        if (with_ud) {
          for (cnt = 0; cnt < domains_up->uniq_motif_count; cnt++) {
            u = domains_up->uniq_motif_size[cnt];
            if (hc->up_ml[i] >= u) {
              temp = ml_helpers->prm_MLbu[u]
                     * expMLbase[u]
                     * domains_up->exp_energy_cb(fc,
                                                 i,
                                                 i + u,
                                                 VRNA_UNSTRUCTURED_DOMAIN_MB_LOOP | VRNA_UNSTRUCTURED_DOMAIN_MOTIF,
                                                 domains_up->data);

              if (sc)
                if (sc->exp_energy_up)
                  temp *= sc->exp_energy_up[i][u];

              ppp += temp;
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
            (G[kl] == 0.))
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

      s5  = ((k > 1) && (sn[k] == sn[k - 1])) ? S1[k - 1] : -1;
      s3  = ((l < n) && (sn[l + 1] == sn[l])) ? S1[l + 1] : -1;

      if ((with_gquad) &&
          (qb[kl] == 0.)) {
        temp *= G[kl] *
                expMLstem;
      } else if (hc->mx[l * n + k] & VRNA_CONSTRAINT_CONTEXT_MB_LOOP_ENC) {
        if (tt == 0)
          tt = 7;

        temp *= exp_E_MLstem(tt, s5, s3, pf_params);
      }

      probs[kl] += temp *
                   scale[2];

      if (probs[kl] > (*Qmax)) {
        (*Qmax) = probs[kl];
        if ((*Qmax) > max_real / 10.)
          vrna_message_warning("P close to overflow: %d %d %g %g\n",
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
                                    int                   *ov)
{
  unsigned char     tt;
  short             **S, **S5, **S3;
  unsigned int      **a2s, s, n_seq, *sn;
  int               i, j, k, n, ii, kl, ij, lj, turn, *my_iindx, *jindx, *pscore, with_gquad;
  FLT_OR_DBL        temp, ppp, prm_MLb, prmt, prmt1, *qb, *probs, *qm, *G, *scale,
                    *expMLbase, expMLclosing, expMLstem;
  double            max_real, kTn;
  vrna_exp_param_t  *pf_params;
  vrna_md_t         *md;
  vrna_hc_t         *hc;
  vrna_sc_t         **scs;

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
  turn          = md->min_loop_size;
  qb            = fc->exp_matrices->qb;
  qm            = fc->exp_matrices->qm;
  G             = fc->exp_matrices->G;
  probs         = fc->exp_matrices->probs;
  scale         = fc->exp_matrices->scale;
  expMLbase     = fc->exp_matrices->expMLbase;
  expMLclosing  = pf_params->expMLclosing;
  with_gquad    = md->gquad;
  hc            = fc->hc;
  scs           = fc->scs;
  expMLstem     = (with_gquad) ? (FLT_OR_DBL)pow(exp_E_MLstem(0, -1, -1, pf_params), (double)n_seq) : 0;

  kTn       = pf_params->kT / 10.;   /* kT in cal/mol  */
  prm_MLb   = 0.;
  max_real  = (sizeof(FLT_OR_DBL) == sizeof(float)) ? FLT_MAX : DBL_MAX;

  if (sn[l + 1] != sn[l]) {
    /* set prm_l to 0 to get prm_l1 in the next round to be 0 */
    for (i = 0; i <= n; i++)
      ml_helpers->prm_l[i] = 0;
  } else {
    for (k = 2; k < l - turn; k++) {
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
            ppp = probs[ij]
                  * qm[lj];

            for (s = 0; s < n_seq; s++) {
              tt  = vrna_get_ptype_md(S[s][j], S[s][i], md);
              ppp  *= exp_E_MLstem(tt, S5[s][j], S3[s][i], pf_params);
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

        ii  = my_iindx[i];  /* ii-j=[i,j]     */

        if (hc->mx[(l + 1) * n + i] & VRNA_CONSTRAINT_CONTEXT_MB_LOOP) {
          prmt1 = probs[ii - (l + 1)] *
                  (FLT_OR_DBL)pow(expMLclosing, (double)n_seq);

          for (s = 0; s < n_seq; s++) {
            tt    = vrna_get_ptype_md(S[s][l + 1], S[s][i], md);
            prmt1 *= exp_E_MLstem(tt, S5[s][l + 1], S3[s][i], pf_params);
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
        ppp = ml_helpers->prm_l1[i] * expMLbase[1];

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
            (G[kl] == 0.))
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
        temp *= G[kl] *
                expMLstem;
      } else {
        if (hc->mx[l * n + k] & VRNA_CONSTRAINT_CONTEXT_MB_LOOP_ENC) {
          for (s = 0; s < n_seq; s++) {
            tt    = vrna_get_ptype_md(S[s][k], S[s][l], md);
            temp  *= exp_E_MLstem(tt, S5[s][k], S3[s][l], pf_params);
          }
        }
      }

      probs[kl] += temp *
                   scale[2] *
                   exp(pscore[jindx[l] + k] / kTn);

      if (probs[kl] > (*Qmax)) {
        (*Qmax) = probs[kl];
        if ((*Qmax) > max_real / 10.)
          vrna_message_warning("P close to overflow: %d %d %g %g\n",
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
                            int                   l)
{
  unsigned char     type;
  char              *ptype;
  short             *S1;
  int               i, j, k, n, ij, kl, u1, u2, *my_iindx, *jindx;
  FLT_OR_DBL        tmp2, qe, *G, *probs, *scale;
  vrna_exp_param_t  *pf_params;

  n         = (int)fc->length;
  S1        = fc->sequence_encoding;
  ptype     = fc->ptype;
  my_iindx  = fc->iindx;
  jindx     = fc->jindx;
  pf_params = fc->exp_params;
  G         = fc->exp_matrices->G;
  probs     = fc->exp_matrices->probs;
  scale     = fc->exp_matrices->scale;

  /* 2.5. bonding k,l as gquad enclosed by i,j */
  double *expintern = &(pf_params->expinternal[0]);

  if (l < n - 3) {
    for (k = 2; k <= l - VRNA_GQUAD_MIN_BOX_SIZE + 1; k++) {
      kl = my_iindx[k] - l;
      if (G[kl] == 0.)
        continue;

      tmp2  = 0.;
      i     = k - 1;
      for (j = MIN2(l + MAXLOOP + 1, n); j > l + 3; j--) {
        ij    = my_iindx[i] - j;
        type  = (unsigned char)ptype[jindx[j] + i];
        if (!type)
          continue;

        u1    = j - l - 1;
        qe    = (type > 2) ? pf_params->expTermAU : 1.;
        tmp2  += probs[ij]
                 * qe
                 * (FLT_OR_DBL)expintern[u1]
                 * pf_params->expmismatchI[type][S1[i + 1]][S1[j - 1]]
                 * scale[u1 + 2];
      }
      probs[kl] += tmp2 * G[kl];
    }
  }

  if (l < n - 1) {
    for (k = 3; k <= l - VRNA_GQUAD_MIN_BOX_SIZE + 1; k++) {
      kl = my_iindx[k] - l;
      if (G[kl] == 0.)
        continue;

      tmp2 = 0.;
      for (i = MAX2(1, k - MAXLOOP - 1); i <= k - 2; i++) {
        u1 = k - i - 1;
        for (j = l + 2; j <= MIN2(l + MAXLOOP - u1 + 1, n); j++) {
          ij    = my_iindx[i] - j;
          type  = (unsigned char)ptype[jindx[j] + i];
          if (!type)
            continue;

          u2    = j - l - 1;
          qe    = (type > 2) ? pf_params->expTermAU : 1.;
          tmp2  += probs[ij]
                   * qe
                   * (FLT_OR_DBL)expintern[u1 + u2]
                   * pf_params->expmismatchI[type][S1[i + 1]][S1[j - 1]]
                   * scale[u1 + u2 + 2];
        }
      }
      probs[kl] += tmp2 * G[kl];
    }
  }

  if (l < n) {
    for (k = 4; k <= l - VRNA_GQUAD_MIN_BOX_SIZE + 1; k++) {
      kl = my_iindx[k] - l;
      if (G[kl] == 0.)
        continue;

      tmp2  = 0.;
      j     = l + 1;
      for (i = MAX2(1, k - MAXLOOP - 1); i < k - 3; i++) {
        ij    = my_iindx[i] - j;
        type  = (unsigned char)ptype[jindx[j] + i];
        if (!type)
          continue;

        u2    = k - i - 1;
        qe    = (type > 2) ? pf_params->expTermAU : 1.;
        tmp2  += probs[ij]
                 * qe
                 * (FLT_OR_DBL)expintern[u2]
                 * pf_params->expmismatchI[type][S1[i + 1]][S1[j - 1]]
                 * scale[u2 + 2];
      }
      probs[kl] += tmp2 * G[kl];
    }
  }
}


PRIVATE void
compute_gquad_prob_internal_comparative(vrna_fold_compound_t  *fc,
                                        int                   l)
{
  unsigned char     type;
  short             **S, **S5, **S3;
  unsigned int      **a2s, s, n_seq;
  int               i, j, k, n, ij, kl, u1, u2, u1_local, u2_local, *my_iindx;
  FLT_OR_DBL        tmp2, qe, *G, *qb, *probs, *scale;
  vrna_exp_param_t  *pf_params;
  vrna_md_t         *md;

  n         = (int)fc->length;
  n_seq     = fc->n_seq;
  S         = fc->S;
  S5        = fc->S5;
  S3        = fc->S3;
  a2s       = fc->a2s;
  my_iindx  = fc->iindx;
  pf_params = fc->exp_params;
  G         = fc->exp_matrices->G;
  qb        = fc->exp_matrices->qb;
  probs     = fc->exp_matrices->probs;
  scale     = fc->exp_matrices->scale;
  md        = &(pf_params->model_details);

  /* 2.5. bonding k,l as gquad enclosed by i,j */
  double *expintern = &(pf_params->expinternal[0]);

  if (l < n - 3) {
    for (k = 2; k <= l - VRNA_GQUAD_MIN_BOX_SIZE + 1; k++) {
      kl = my_iindx[k] - l;
      if (G[kl] == 0.)
        continue;

      tmp2  = 0.;
      i     = k - 1;
      for (j = MIN2(l + MAXLOOP + 1, n); j > l + 3; j--) {
        ij    = my_iindx[i] - j;
        if (qb[ij] == 0.)
          continue;

        qe  = 1.;
        u1  = j - l - 1;

        for (s = 0; s < n_seq; s++) {
          type      = vrna_get_ptype_md(S[s][i], S[s][j], md);
          u1_local  = a2s[s][j - 1] - a2s[s][l];
          qe        *=  (FLT_OR_DBL)expintern[u1_local];

          if (md->dangles == 2)
            qe *=  (FLT_OR_DBL)pf_params->expmismatchI[type][S3[s][i]][S5[s][j]];

          if (type > 2)
            qe *= (FLT_OR_DBL)pf_params->expTermAU;
        }

        tmp2  += probs[ij]
                 * qe
                 * scale[u1 + 2];
      }
      probs[kl] += tmp2 * G[kl];
    }
  }

  if (l < n - 1) {
    for (k = 3; k <= l - VRNA_GQUAD_MIN_BOX_SIZE + 1; k++) {
      kl = my_iindx[k] - l;
      if (G[kl] == 0.)
        continue;

      tmp2 = 0.;
      for (i = MAX2(1, k - MAXLOOP - 1); i <= k - 2; i++) {
        u1 = k - i - 1;
        for (j = l + 2; j <= MIN2(l + MAXLOOP - u1 + 1, n); j++) {
          ij    = my_iindx[i] - j;
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

          tmp2  += probs[ij]
                   * qe
                   * scale[u1 + u2 + 2];
        }
      }
      probs[kl] += tmp2 * G[kl];
    }
  }

  if (l < n) {
    for (k = 4; k <= l - VRNA_GQUAD_MIN_BOX_SIZE + 1; k++) {
      kl = my_iindx[k] - l;
      if (G[kl] == 0.)
        continue;

      tmp2  = 0.;
      j     = l + 1;
      for (i = MAX2(1, k - MAXLOOP - 1); i < k - 3; i++) {
        ij    = my_iindx[i] - j;
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

        tmp2  += probs[ij]
                 * qe
                 * scale[u2 + 2];
      }
      probs[kl] += tmp2 * G[kl];
    }
  }
}


/* outside recursion of pf cofolding */
PRIVATE int
pf_co_bppm(vrna_fold_compound_t *vc,
           char                 *structure)
{
  unsigned int      *sn;
  int               n, i, j, k, l, ij, kl, type, turn, ov = 0,
                    *my_iindx, *jindx, cp;
  FLT_OR_DBL        temp, Qmax = 0;
  FLT_OR_DBL        *probs, *q1k, *qln, *q, *qb, *scale;
  vrna_exp_param_t  *pf_params;
  vrna_md_t         *md;
  short             *S, *S1;
  char              *ptype;
  vrna_mx_pf_t      *matrices;
  int               *rtype;
  helper_arrays       *ml_helpers;
  constraints_helper  *constraints;

  n         = vc->length;
  cp        = vc->cutpoint;
  pf_params = vc->exp_params;
  md        = &(pf_params->model_details);
  S         = vc->sequence_encoding2;
  S1        = vc->sequence_encoding;
  sn        = vc->strand_number;
  jindx     = vc->jindx;
  my_iindx  = vc->iindx;
  ptype     = vc->ptype;
  rtype     = &(md->rtype[0]);
  turn      = md->min_loop_size;

  matrices  = vc->exp_matrices;
  probs     = matrices->probs;
  scale     = matrices->scale;
  q1k       = matrices->q1k;
  qln       = matrices->qln;
  q         = matrices->q;
  qb        = matrices->qb;

  /* hard code min_loop_size to 0, since we can not be sure yet that this is already the case */
  turn = 0;

  /* backtracking to construct binding probabilities of pairs*/
  if ((S != NULL) && (S1 != NULL)) {
    FLT_OR_DBL    *Qlout, *Qrout;

    int           corr_size       = 5;
    int           corr_cnt        = 0;
    vrna_ep_t     *bp_correction  = vrna_alloc(sizeof(vrna_ep_t) * corr_size);

    ml_helpers  = get_ml_helper_arrays(vc);
    constraints = get_constraints_helper(vc);

    Qmax  = 0;
    Qrout = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * (n + 2));
    Qlout = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * (cp + 2));

    for (k = 1; k <= n; k++) {
      q1k[k]  = q[my_iindx[1] - k];
      qln[k]  = q[my_iindx[k] - n];
    }
    q1k[0]      = 1.0;
    qln[n + 1]  = 1.0;

    /* init diagonal entries unable to pair in pr matrix */
    for (i = 1; i <= n; i++)
      for (j = i; j <= MIN2(i + turn, n); j++)
        probs[my_iindx[i] - j] = 0.;

    /* 1. external loop pairs, i.e. pairs not enclosed by any other pair (or external loop for circular RNAs) */
    compute_bpp_external(vc);

    /* 2. all cases where base pair (k,l) is enclosed by another pair (i,j) */
    l = n;
    compute_bpp_internal(vc,
                         l,
                         &bp_correction,
                         &corr_cnt,
                         &corr_size,
                         &Qmax,
                         &ov,
                         constraints);

    for (l = n - 1; l > turn + 1; l--) {
      compute_bpp_internal(vc,
                           l,
                           &bp_correction,
                           &corr_cnt,
                           &corr_size,
                           &Qmax,
                           &ov,
                           constraints);

      compute_bpp_multibranch(vc,
                              l,
                              ml_helpers,
                              &Qmax,
                              &ov);

      /* computation of .(..(...)..&..). type features? */
      if (vc->strands <= 1)
        continue;                     /* no .(..(...)..&..). type features */

      if (l <= 2)
        continue;                     /* no .(..(...)..&..). type features */

      /*new version with O(n^3)??*/
      if (l > cp) {
        int t, kt;
        for (t = n; t > l; t--) {
          for (k = 1; k < cp; k++) {
            int samestrand;
            kt = my_iindx[k] - t;

            samestrand  = (sn[k + 1] == sn[k]) ? 1 : 0;
            type        = rtype[vrna_get_ptype(jindx[t] + k, ptype)];

            temp = probs[kt]
                   * vrna_exp_E_ext_stem(type, S1[t - 1], samestrand ? S1[k + 1] : -1, pf_params)
                   * scale[2];

            if (l + 1 < t)
              temp *= q[my_iindx[l + 1] - (t - 1)];

            if (samestrand)
              temp *= q[my_iindx[k + 1] - (cp - 1)];

            Qrout[l] += temp;
          }
        }

        for (k = l - 1; k >= cp; k--) {
          if (qb[my_iindx[k] - l]) {
            kl    = my_iindx[k] - l;
            type  = vrna_get_ptype(jindx[l] + k, ptype);
            temp  = Qrout[l];

            temp *= vrna_exp_E_ext_stem(type,
                                  (k > cp) ? S1[k - 1] : -1,
                                  (l < n) ? S1[l + 1] : -1,
                                  pf_params);
            if (k > cp)
              temp *= q[my_iindx[cp] - (k - 1)];

            probs[kl] += temp;
          }
        }
      } else if (l == cp) {
        int t, sk, s;
        for (t = 2; t < cp; t++) {
          for (s = 1; s < t; s++) {
            for (k = cp; k <= n; k++) {
              sk = my_iindx[s] - k;
              if (qb[sk]) {
                int samestrand;
                samestrand  = (sn[k] == sn[k - 1]) ? 1 : 0;
                type        = rtype[vrna_get_ptype(jindx[k] + s, ptype)];

                temp = probs[sk]
                       * vrna_exp_E_ext_stem(type, samestrand ? S1[k - 1] : -1, S1[s + 1], pf_params)
                       * scale[2];
                if (s + 1 < t)
                  temp *= q[my_iindx[s + 1] - (t - 1)];

                if (samestrand)
                  temp *= q[my_iindx[cp] - (k - 1)];

                Qlout[t] += temp;
              }
            }
          }
        }
      } else if (l < cp) {
        for (k = 1; k < l; k++) {
          if (qb[my_iindx[k] - l]) {
            type  = vrna_get_ptype(jindx[l] + k, ptype);
            temp  = Qlout[k];

            temp *= vrna_exp_E_ext_stem(type,
                                  (k > 1) ? S1[k - 1] : -1,
                                  (l < (cp - 1)) ? S1[l + 1] : -1,
                                  pf_params);
            if (l + 1 < cp)
              temp *= q[my_iindx[l + 1] - (cp - 1)];

            probs[my_iindx[k] - l] += temp;
          }
        }
      }
    }  /* end for (l=..)   */
    free(Qlout);
    free(Qrout);
    for (i = 1; i <= n; i++)
      for (j = i + turn + 1; j <= n; j++) {
        ij        = my_iindx[i] - j;
        probs[ij] *= qb[ij];
      }

    if (structure != NULL) {
      char *s = vrna_db_from_probs(probs, (unsigned int)n);
      memcpy(structure, s, n);
      structure[n] = '\0';
      free(s);
    }

    /* clean up */
    free_ml_helper_arrays(ml_helpers);

    free_constraints_helper(constraints);

    free(bp_correction);
  }   /* end if (do_backtrack) */
  else {
    vrna_message_warning("bppm calculations have to be done after calling forward recursion");
    return 0;
  }

  if (ov > 0)
    vrna_message_warning("%d overflows occurred while backtracking;\n"
                         "you might try a smaller pf_scale than %g\n",
                         ov, pf_params->pf_scale);

  return 1;
}


PRIVATE INLINE void
ud_outside_ext_loops(vrna_fold_compound_t *vc)
{
  int         i, j, u, n, cnt, *motif_list, *hc_up;
  FLT_OR_DBL  *q1k, *qln, temp, *scale;
  vrna_sc_t   *sc;
  vrna_ud_t   *domains_up;

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
        j = i + u - 1;
        if (j <= n) {
          if (hc_up[i] >= u) {
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
  int         i, j, k, l, kl, *my_iindx, u, n, cnt, *motif_list, *hc_up;
  FLT_OR_DBL  temp, outside, exp_motif_en, *probs, q1, q2;

  vrna_ud_t   *domains_up, *ud_bak;

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
          if (hc_up[i] >= u) {
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
                  temp            = vrna_exp_E_hp_loop(vc, k, l);
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
  int         i, j, k, l, kl, *my_iindx, u, u1, u2, n, turn, m;
  FLT_OR_DBL  temp, outside, exp_motif_en, *probs, q1, q2, **qq_ud, **pp_ud;

  vrna_ud_t   *domains_up, *ud_bak;

  n           = vc->length;
  my_iindx    = vc->iindx;
  probs       = vc->exp_matrices->probs;
  domains_up  = vc->domains_up;
  turn        = vc->exp_params->model_details.min_loop_size;

  qq_ud = (FLT_OR_DBL **)vrna_alloc(sizeof(FLT_OR_DBL *) * (n + 1));
  for (k = 0; k < domains_up->uniq_motif_count; k++) {
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
  for (k = 0; k < domains_up->uniq_motif_count; k++)
    pp_ud[k] = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * (n + 1));

  for (k = 1; k < n; k++) {
    for (l = k + turn + 1; l <= n; l++) {
      kl = my_iindx[k] - l;
      if (probs[kl] > 0.) {
        ud_bak          = vc->domains_up;
        vc->domains_up  = NULL;
        temp            = vrna_exp_E_hp_loop(vc, k, l);
        vc->domains_up  = ud_bak;
        temp            *= probs[kl];

        if (temp > 0.) {
          for (m = 0; m < domains_up->uniq_motif_count; m++) {
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

  for (k = 0; k < domains_up->uniq_motif_count; k++) {
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
  int         i, j, k, l, p, q, pq, kl, u, n, cnt, *motif_list, *my_iindx,
              *hc_up, kmin, pmax, qmin, lmax, turn;
  FLT_OR_DBL  temp, q1, q2, q3, exp_motif_en, outside,
              *probs, *qb;
  vrna_ud_t   *domains_up, *ud_bak;

  n           = vc->length;
  my_iindx    = vc->iindx;
  qb          = vc->exp_matrices->qb;
  probs       = vc->exp_matrices->probs;
  hc_up       = vc->hc->up_int;
  domains_up  = vc->domains_up;
  turn        = vc->exp_params->model_details.min_loop_size;

  for (i = 2; i <= n; i++) {
    motif_list = vrna_ud_get_motif_size_at(vc, i, VRNA_UNSTRUCTURED_DOMAIN_INT_LOOP);

    /* 3. Interior loops */
    if (motif_list) {
      cnt = 0;
      while (-1 != (u = motif_list[cnt])) {
        outside = 0.;
        j       = i + u - 1;

        if (j < n) {
          if (hc_up[i] >= u) {
            exp_motif_en = domains_up->exp_energy_cb(vc,
                                                     i,
                                                     j,
                                                     VRNA_UNSTRUCTURED_DOMAIN_INT_LOOP | VRNA_UNSTRUCTURED_DOMAIN_MOTIF,
                                                     domains_up->data);

            /* 3.1 motif is within 5' loop */
            kmin  = j - MAXLOOP - 1;
            kmin  = MAX2(kmin, 1);
            for (k = kmin; k < i; k++) {
              pmax  = k + MAXLOOP + 1;
              pmax  = MIN2(pmax, n);
              for (p = j + 1; p < n; p++)
                for (q = p + turn + 1; q < n; q++) {
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
                      temp            = vrna_exp_E_interior_loop(vc, k, l, p, q);
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
            for (k = 1; k < i - turn - 2; k++) {
              pmax  = k + i + MAXLOOP - j;
              pmax  = MIN2(pmax, n);
              for (p = k + 1; p <= pmax; p++) {
                qmin  = p + j - k - MAXLOOP - 1;
                qmin  = MAX2(qmin, p + turn + 1);
                for (q = i - 1; q >= qmin; q--) {
                  pq = my_iindx[p] - q;
                  if (qb[pq] == 0.)
                    continue;

                  lmax  = k + q - p + MAXLOOP + 2;
                  lmax  = MIN2(lmax, n);
                  for (l = j + 1; l < lmax; l++) {
                    kl = my_iindx[k] - l;
                    if (probs[kl] > 0.) {
                      ud_bak          = vc->domains_up;
                      vc->domains_up  = NULL;
                      temp            = vrna_exp_E_interior_loop(vc, k, l, p, q);
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
  int           i, j, k, l, p, q, pq, kl, u, n, *my_iindx, pmax, qmin, turn,
                u1, u2, uu1, uu2, u2_max, m;
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
  turn              = vc->exp_params->model_details.min_loop_size;

  qqk = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * (n + 1));
  qql = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * (n + 1));
  qqp = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * (n + 1));

  qq_ud = (FLT_OR_DBL **)vrna_alloc(sizeof(FLT_OR_DBL *) * (n + 1));
  for (k = 0; k < domains_up->uniq_motif_count; k++) {
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
  for (k = 0; k < domains_up->uniq_motif_count; k++)
    pp_ud[k] = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * (n + 1));

  for (k = 1; k < n; k++) {
    for (l = k + 1; l <= MIN2(k + MAXLOOP, n); l++) {
      qqk[l] = domains_up->exp_energy_cb(vc,
                                         k + 1, l,
                                         VRNA_UNSTRUCTURED_DOMAIN_INT_LOOP,
                                         domains_up->data);
    }
    for (l = k + turn + 1 + 3; l <= n; l++) {
      kl = my_iindx[k] - l;
      if (probs[kl] == 0.)
        continue;

      if (hard_constraints[k * n + l] & VRNA_CONSTRAINT_CONTEXT_INT_LOOP) {
        for (i = l - 1; i > MAX2(k, l - MAXLOOP - 1); i--) {
          qql[i] = domains_up->exp_energy_cb(vc,
                                             i, l - 1,
                                             VRNA_UNSTRUCTURED_DOMAIN_INT_LOOP,
                                             domains_up->data);
        }

        pmax  = k + MAXLOOP + 1;
        pmax  = MIN2(pmax, l - turn);
        for (p = k + 1; p < pmax; p++) {
          u1      = p - k - 1;
          u2_max  = MAXLOOP - u1;
          qmin    = l - 1 - u2_max;
          qmin    = MAX2(qmin, p + turn + 1);
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
              temp            = vrna_exp_E_interior_loop(vc, k, l, p, q);
              vc->domains_up  = ud_bak;
              temp            *= probs[kl] * qb[pq];

              q3 = 1.;
              if (u2 > 0)
                q3 += qql[q + 1];

              temp5 = temp * q3;
              temp3 = temp * q5;

              /* loop over all available motifs */
              for (m = 0; m < domains_up->uniq_motif_count; m++) {
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
  for (k = 0; k < domains_up->uniq_motif_count; k++) {
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
  int               i, j, k, l, kl, jkl, *my_iindx, u, n, cnt, *motif_list,
                    *hc_up, turn, tt, *jindx, *rtype, up, ud_max_size;
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
  turn          = md->min_loop_size;
  rtype         = &(md->rtype[0]);
  expMLbase     = vc->exp_matrices->expMLbase;
  expMLclosing  = pf_params->expMLclosing;

  for (ud_max_size = u = 0; u < domains_up->uniq_motif_count; u++)
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
          if (hc_up[i] >= u) {
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
            for (l = j + turn + 1; l <= n; l++) {
              for (k = i - turn - 1; k > 0; k--) {
                kl = my_iindx[k] - l;
                if (probs[kl] > 0.) {
                  if (hc[n * l + k] & VRNA_CONSTRAINT_CONTEXT_MB_LOOP) {
                    /* respect hard constraints */
                    FLT_OR_DBL qqq;
                    jkl = jindx[l] + k;
                    tt  = rtype[vrna_get_ptype(jkl, ptype)];
                    qqq = probs[kl]
                          * qm[my_iindx[k + 1] - (i - 1)]
                          * qm[my_iindx[j + 1] - (l - 1)]
                          * exp_E_MLstem(tt, S[l - 1], S[k + 1], pf_params)
                          * expMLclosing
                          * scale[2];

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
            for (l = j + turn + 1; l <= n; l++) {
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
                  temp  = probs[kl]
                          * expMLbase[up]
                          * exp_E_MLstem(tt, S[l - 1], S[k + 1], pf_params)
                          * expMLclosing
                          * scale[2];
                  if (sc) {
                    if (sc->exp_energy_bp)
                      temp *= sc->exp_energy_bp[jkl];

                    if (sc->exp_energy_up)
                      temp *= sc->exp_energy_up[k + 1][up];
                  }

                  lqq += temp;
                  lqq += temp
                         * domains_up->exp_energy_cb(vc,
                                                     k + 1, i - 1,
                                                     VRNA_UNSTRUCTURED_DOMAIN_MB_LOOP,
                                                     domains_up->data);
                }
              }

              for (u = j + turn + 1; u < l - turn; u++) {
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
                for (cnt = 0; cnt < domains_up->uniq_motif_count; cnt++) {
                  int size = domains_up->uniq_motif_size[cnt];
                  if ((u < l - size) &&
                      (hc_up[l - size] >= size)) {
                    temp = qm1ui[size][u]
                           * expMLbase[size]
                           * domains_up->exp_energy_cb(vc,
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
                  temp  = qb[my_iindx[u] - (l - 1)]
                          * exp_E_MLstem(tt, S[u - 1], S[l], pf_params);

                  qm1ui[0][u] += temp;
                }

                rqq += qm[my_iindx[j + 1] - (u - 1)] * qm1ui[0][u];
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
            for (k = i - turn - 1; k > 0; k--) {
              FLT_OR_DBL  lqq = 0.;
              FLT_OR_DBL  rqq = 0;

              /* update qmli[k] = qm1[k,i-1] */
              for (qmli[k] = 0., u = k + turn + 1; u < i; u++) {
                /* respect hard constraints */
                if (hc[n * u + k] & VRNA_CONSTRAINT_CONTEXT_MB_LOOP_ENC) {
                  up = (i - 1) - (u + 1) + 1;
                  if (hc_up[u + 1] >= up) {
                    temp = qb[my_iindx[k] - u]
                           * expMLbase[up];

                    /* add soft constraints */
                    if (sc)
                      if (sc->exp_energy_up)
                        temp *= sc->exp_energy_up[u + 1][up];

                    qmli[k] += temp;

                    /* add contributions of other motifs within [u+1:i-1] */
                    qmli[k] += temp
                               * domains_up->exp_energy_cb(vc,
                                                           u + 1, i - 1,
                                                           VRNA_UNSTRUCTURED_DOMAIN_MB_LOOP,
                                                           domains_up->data);
                  }
                }
              }

              for (u = k + turn; u < i - turn; u++)
                lqq += qm[my_iindx[k + 1] - (u - 1)]
                       * qmli[u];

              for (l = j + 1; l <= n; l++) {
                kl = my_iindx[k] - l;
                if (hc[n * k + l] & VRNA_CONSTRAINT_CONTEXT_MB_LOOP) {
                  int up, jkl;
                  jkl = jindx[l] + k;
                  tt  = rtype[vrna_get_ptype(jkl, ptype)];
                  up  = l - j - 1;
                  if (hc_up[j + 1] >= up) {
                    temp = probs[kl]
                           * exp_E_MLstem(tt, S[l - 1], S[k + 1], pf_params)
                           * expMLclosing
                           * scale[2]
                           * expMLbase[up];

                    if (sc) {
                      if (sc->exp_energy_bp)
                        temp *= sc->exp_energy_bp[jkl];

                      if (sc->exp_energy_up)
                        temp *= sc->exp_energy_up[j + 1][up];
                    }

                    rqq += temp;

                    /* add contributions of other motifs within [j+1:l-1] */
                    rqq += temp
                           * domains_up->exp_energy_cb(vc,
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

            outside += exp_motif_ml_right
                       * exp_motif_en;
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
  int               i, j, k, l, kl, jkl, *my_iindx, u, n, cnt, *motif_list,
                    *hc_up, turn, tt, *jindx, *rtype, up, ud_max_size;
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
  turn          = md->min_loop_size;
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
          if (hc_up[i] >= u) {
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
            for (l = j + turn + 1; l <= n; l++) {
              for (k = i - turn - 1; k > 0; k--) {
                kl = my_iindx[k] - l;
                if (probs[kl] > 0.) {
                  jkl = jindx[l] + k;
                  if (hc[l * n + k] & VRNA_CONSTRAINT_CONTEXT_MB_LOOP) {
                    /* respect hard constraints */
                    FLT_OR_DBL qqq;
                    tt  = rtype[vrna_get_ptype(jkl, ptype)];
                    qqq = probs[kl]
                          * qqmi[k + 1]
                          * qqmj[l - 1]
                          * exp_E_MLstem(tt, S[l - 1], S[k + 1], pf_params)
                          * expMLclosing
                          * scale[2];

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
            for (l = j + turn + 1; l <= n; l++) {
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
                  temp  = probs[kl]
                          * expMLbase[up]
                          * exp_E_MLstem(tt, S[l - 1], S[k + 1], pf_params)
                          * expMLclosing
                          * scale[2];
                  if (sc) {
                    if (sc->exp_energy_bp)
                      temp *= sc->exp_energy_bp[jkl];

                    if (sc->exp_energy_up)
                      temp *= sc->exp_energy_up[k + 1][up];
                  }

                  lqq += temp;
                  lqq += temp
                         * qqi[k + 1];
                }
              }

              for (u = j + turn + 1; u < l - turn; u++) {
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
                for (cnt = 0; cnt < domains_up->uniq_motif_count; cnt++) {
                  int size = domains_up->uniq_motif_size[cnt];
                  if ((u < l - size) &&
                      (hc_up[l - size] >= size)) {
                    temp = qm1ui[size][u]
                           * expMLbase[size]
                           * domains_up->exp_energy_cb(vc,
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
                  temp  = qb[ul]
                          * exp_E_MLstem(tt, S[u - 1], S[l], pf_params);

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

            outside += exp_motif_ml_left
                       * exp_motif_en;

            /* 4.3 Motif is in right-most unpaired stretch of multiloop */
            qmli                = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * n);
            exp_motif_ml_right  = 0.;
            for (k = i - turn - 1; k > 0; k--) {
              FLT_OR_DBL  lqq = 0.;
              FLT_OR_DBL  rqq = 0;

              /* update qmli[k] = qm1[k,i-1] */
              for (qmli[k] = 0., u = k + turn + 1; u < i; u++) {
                int ku = my_iindx[k] - u;
                /* respect hard constraints */
                if (hc[k * n + u] & VRNA_CONSTRAINT_CONTEXT_MB_LOOP_ENC) {
                  up = (i - 1) - (u + 1) + 1;
                  if (hc_up[u + 1] >= up) {
                    temp = qb[ku]
                           * expMLbase[up];

                    /* add soft constraints */
                    if (sc)
                      if (sc->exp_energy_up)
                        temp *= sc->exp_energy_up[u + 1][up];

                    qmli[k] += temp;

                    /* add contributions of other motifs within [u+1:i-1] */
                    qmli[k] += temp
                               * qqi[u + 1];
                  }
                }
              }

              for (u = k + turn; u < i - turn; u++)
                lqq += qm[my_iindx[k + 1] - (u - 1)]
                       * qmli[u];

              for (l = j + 1; l <= n; l++) {
                kl = my_iindx[k] - l;
                if (hc[k * n + l] & VRNA_CONSTRAINT_CONTEXT_MB_LOOP) {
                  int up, jkl;
                  jkl = jindx[l] + k;
                  tt  = rtype[vrna_get_ptype(jkl, ptype)];
                  up  = l - j - 1;
                  if (hc_up[j + 1] >= up) {
                    temp = probs[kl]
                           * exp_E_MLstem(tt, S[l - 1], S[k + 1], pf_params)
                           * expMLclosing
                           * scale[2]
                           * expMLbase[up];

                    if (sc) {
                      if (sc->exp_energy_bp)
                        temp *= sc->exp_energy_bp[jkl];

                      if (sc->exp_energy_up)
                        temp *= sc->exp_energy_up[j + 1][up];
                    }

                    rqq += temp;

                    /* add contributions of other motifs within [j+1:l-1] */
                    rqq += temp
                           * qqj[l - 1];
                  }
                }
              }
              exp_motif_ml_right += rqq * lqq;
            }
            free(qmli);
            qmli = NULL;

            outside += exp_motif_ml_right
                       * exp_motif_en;
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
bppm_circ(vrna_fold_compound_t *vc)
{
  unsigned char     type;
  char              *ptype;
  unsigned char     *hard_constraints, eval;
  short             *S, *S1;
  int               n, i, j, k, l, ij, *rtype, *my_iindx, *jindx, turn;
  FLT_OR_DBL        tmp, tmp2, expMLclosing, *qb, *qm, *qm1, *probs, *scale, *expMLbase, qo;
  vrna_hc_t         *hc;
  vrna_sc_t         *sc;
  vrna_exp_param_t  *pf_params;
  vrna_mx_pf_t      *matrices;
  vrna_md_t         *md;

  FLT_OR_DBL        (*numerator_f)(vrna_fold_compound_t *vc,
                                   int                  i,
                                   int                  j);

  pf_params         = vc->exp_params;
  md                = &(pf_params->model_details);
  S                 = vc->sequence_encoding2;
  S1                = vc->sequence_encoding;
  my_iindx          = vc->iindx;
  jindx             = vc->jindx;
  ptype             = vc->ptype;
  turn              = md->min_loop_size;
  hc                = vc->hc;
  sc                = vc->sc;
  matrices          = vc->exp_matrices;
  qb                = matrices->qb;
  qm                = matrices->qm;
  qm1               = matrices->qm1;
  probs             = matrices->probs;
  scale             = matrices->scale;
  expMLbase         = matrices->expMLbase;
  qo                = matrices->qo;
  hard_constraints  = hc->mx;

  expMLclosing  = pf_params->expMLclosing;
  rtype         = &(pf_params->model_details.rtype[0]);
  n             = S[0];

  switch (vc->type) {
    case  VRNA_FC_TYPE_SINGLE:
      numerator_f = numerator_single;
      break;
    case  VRNA_FC_TYPE_COMPARATIVE:
      numerator_f = numerator_comparative;
      break;
    default:
      numerator_f = NULL;
      break;
  }

  /* 1. exterior pair i,j */
  for (i = 1; i <= n; i++) {
    for (j = i; j <= MIN2(i + turn, n); j++)
      probs[my_iindx[i] - j] = 0;
    for (j = i + turn + 1; j <= n; j++) {
      ij = my_iindx[i] - j;
      if (qb[ij] > 0.) {
        probs[ij] = numerator_f(vc, i, j) / qo;

        type = vrna_get_ptype(jindx[j] + i, ptype);

        unsigned char rt = rtype[type];

        /* 1.1. Exterior Hairpin Contribution */
        tmp2 = vrna_exp_E_hp_loop(vc, j, i);

        if (hard_constraints[i * n + j] & VRNA_CONSTRAINT_CONTEXT_INT_LOOP) {
          /* 1.2. Exterior Interior Loop Contribution                     */
          /* 1.2.1. i,j  delimtis the "left" part of the interior loop    */
          /* (j,i) is "outer pair"                                        */
          for (k = 1; k < i - turn - 1; k++) {
            int ln1, ln3, lstart;
            ln1 = k - 1;
            ln3 = n - j;

            if (hc->up_int[j + 1] < (ln1 + ln3))
              break;

            if ((ln1 + ln3) > MAXLOOP)
              break;

            lstart = (ln1 + ln3) + i - 1 - MAXLOOP;
            if (lstart < k + turn + 1)
              lstart = k + turn + 1;

            for (l = lstart; l < i; l++) {
              int ln2, type_2;
              ln2 = i - l - 1;

              if (hc->up_int[l + 1] < ln2)
                continue;

              if ((ln1 + ln2 + ln3) > MAXLOOP)
                continue;

              eval = (hard_constraints[k * n + l] & VRNA_CONSTRAINT_CONTEXT_INT_LOOP) ? 1 : 0;
              if (hc->f)
                eval = hc->f(k, l, i, j, VRNA_DECOMP_PAIR_IL, hc->data);

              if (eval) {
                type_2 = vrna_get_ptype(jindx[l] + k, ptype);

                tmp = qb[my_iindx[k] - l]
                      * exp_E_IntLoop(ln1 + ln3,
                                      ln2,
                                      rt,
                                      rtype[type_2],
                                      S1[j + 1],
                                      S1[i - 1],
                                      S1[k - 1],
                                      S1[l + 1],
                                      pf_params)
                      * scale[ln1 + ln2 + ln3];

                if (sc) {
                  if (sc->exp_energy_up)
                    tmp *= sc->exp_energy_up[l + 1][ln2] *
                           sc->exp_energy_up[j + 1][ln3] *
                           sc->exp_energy_up[1][ln1];

                  if (sc->exp_f)
                    tmp *= sc->exp_f(k, l, i, j, VRNA_DECOMP_PAIR_IL, sc->data);

                  if (((ln1 + ln2 + ln3) == 0) &&
                      (sc->exp_energy_stack)) {
                    tmp *= sc->exp_energy_stack[k] *
                           sc->exp_energy_stack[l] *
                           sc->exp_energy_stack[i] *
                           sc->exp_energy_stack[j];
                  }
                }

                tmp2 += tmp;
              }
            }
          }

          /* 1.2.2. i,j  delimtis the "right" part of the interior loop  */
          for (k = j + 1; k < n - turn; k++) {
            int ln1, lstart;
            ln1 = k - j - 1;

            if (hc->up_int[j + 1] < ln1)
              break;

            if ((ln1 + i - 1) > MAXLOOP)
              break;

            lstart = ln1 + i - 1 + n - MAXLOOP;
            if (lstart < k + turn + 1)
              lstart = k + turn + 1;

            for (l = lstart; l <= n; l++) {
              int ln2, ln3, type_2;
              ln2 = i - 1;
              ln3 = n - l;

              if (hc->up_int[l + 1] < (ln2 + ln3))
                continue;

              if ((ln1 + ln2 + ln3) > MAXLOOP)
                continue;

              eval = (hard_constraints[k * n + l] & VRNA_CONSTRAINT_CONTEXT_INT_LOOP) ? 1 : 0;
              if (hc->f)
                eval = hc->f(i, j, k, l, VRNA_DECOMP_PAIR_IL, hc->data) ? eval : 0;

              if (eval) {
                type_2  = vrna_get_ptype(jindx[l] + k, ptype);
                tmp     = qb[my_iindx[k] - l] *
                          exp_E_IntLoop(ln2 + ln3,
                                        ln1,
                                        rtype[type_2],
                                        rt,
                                        S1[l + 1],
                                        S1[k - 1],
                                        S1[i - 1],
                                        S1[j + 1],
                                        pf_params) *
                          scale[ln1 + ln2 + ln3];

                if (sc) {
                  if (sc->exp_energy_up)
                    tmp *= sc->exp_energy_up[j + 1][ln1] *
                           sc->exp_energy_up[l + 1][ln3] *
                           sc->exp_energy_up[1][ln2];

                  if (sc->exp_f)
                    tmp *= sc->exp_f(i, j, k, l, VRNA_DECOMP_PAIR_IL, sc->data);

                  if (((ln1 + ln2 + ln3) == 0) &&
                      (sc->exp_energy_stack)) {
                    tmp *= sc->exp_energy_stack[i] *
                           sc->exp_energy_stack[j] *
                           sc->exp_energy_stack[k] *
                           sc->exp_energy_stack[l];
                  }
                }

                tmp2 += tmp;
              }
            }
          }
        }

        /* 1.3 Exterior multiloop decomposition */
        if (hard_constraints[i * n + j] & VRNA_CONSTRAINT_CONTEXT_MB_LOOP) {
          /* 1.3.1 Middle part                    */
          if ((i > turn + 2) &&
              (j < n - turn - 1)) {
            tmp = 0;
            if ((sc) &&
                (sc->exp_f)) {
              if (hc->f) {
                if (hc->f(i, j, i - 1, j + 1, VRNA_DECOMP_PAIR_ML, hc->data) &&
                    hc->f(j + 1, i - 1, n, 1, VRNA_DECOMP_ML_ML_ML, hc->data)) {
                  tmp = qm[my_iindx[1] - i + 1] *
                        qm[my_iindx[j + 1] - n] *
                        expMLclosing *
                        exp_E_MLstem(type,
                                     S1[i - 1],
                                     S1[j + 1],
                                     pf_params) *
                        sc->exp_f(i, j, i - 1, j + 1, VRNA_DECOMP_PAIR_ML, sc->data) *
                        sc->exp_f(j + 1, i - 1, n, 1, VRNA_DECOMP_ML_ML_ML, sc->data);
                }
              } else {
                tmp = qm[my_iindx[1] - i + 1] *
                      qm[my_iindx[j + 1] - n] *
                      expMLclosing *
                      exp_E_MLstem(type,
                                   S1[i - 1],
                                   S1[j + 1],
                                   pf_params) *
                      sc->exp_f(i, j, i - 1, j + 1, VRNA_DECOMP_PAIR_ML, sc->data) *
                      sc->exp_f(j + 1, i - 1, n, 1, VRNA_DECOMP_ML_ML_ML, sc->data);
              }
            } else {
              if (hc->f) {
                if (hc->f(i, j, i - 1, j + 1, VRNA_DECOMP_PAIR_ML, hc->data) &&
                    hc->f(j + 1, i - 1, n, 1, VRNA_DECOMP_ML_ML_ML, hc->data)) {
                  tmp = qm[my_iindx[1] - i + 1] *
                        qm[my_iindx[j + 1] - n] *
                        expMLclosing *
                        exp_E_MLstem(type,
                                     S1[i - 1],
                                     S1[j + 1],
                                     pf_params);
                }
              } else {
                tmp = qm[my_iindx[1] - i + 1] *
                      qm[my_iindx[j + 1] - n] *
                      expMLclosing *
                      exp_E_MLstem(type,
                                   S1[i - 1],
                                   S1[j + 1],
                                   pf_params);
              }
            }

            tmp2 += tmp;
          }

          /* 1.3.2 Left part  */
          if (hc->up_ml[j + 1] >= (n - j)) {
            if (sc) {
              if (hc->f) {
                for (k = turn + 2; k < i - turn - 2; k++) {
                  if (hc->f(i, j, i - 1, j + 1, VRNA_DECOMP_PAIR_ML, hc->data) &&
                      hc->f(i, n, i, j, VRNA_DECOMP_ML_ML, hc->data) &&
                      hc->f(1, i - 1, k, k + 1, VRNA_DECOMP_ML_ML_ML, hc->data)) {
                    tmp = qm[my_iindx[1] - k] *
                          qm1[jindx[i - 1] + k + 1] *
                          expMLbase[n - j] *
                          expMLclosing *
                          exp_E_MLstem(type,
                                       S1[i - 1],
                                       S1[j + 1],
                                       pf_params);


                    if (sc->exp_energy_up)
                      tmp *= sc->exp_energy_up[j + 1][n - j];

                    if (sc->exp_f)
                      tmp *= sc->exp_f(i, j, i - 1, j + 1, VRNA_DECOMP_PAIR_ML, sc->data) *
                             sc->exp_f(i, n, i, j, VRNA_DECOMP_ML_ML, sc->data) *
                             sc->exp_f(1, i - 1, k, k + 1, VRNA_DECOMP_ML_ML_ML, sc->data);

                    tmp2 += tmp;
                  }
                }
              } else {
                for (k = turn + 2; k < i - turn - 2; k++) {
                  tmp = qm[my_iindx[1] - k] *
                        qm1[jindx[i - 1] + k + 1] *
                        expMLbase[n - j] *
                        expMLclosing *
                        exp_E_MLstem(type,
                                     S1[i - 1],
                                     S1[j + 1],
                                     pf_params);


                  if (sc->exp_energy_up)
                    tmp *= sc->exp_energy_up[j + 1][n - j];

                  if (sc->exp_f)
                    tmp *= sc->exp_f(i, j, i - 1, j + 1, VRNA_DECOMP_PAIR_ML, sc->data) *
                           sc->exp_f(i, n, i, j, VRNA_DECOMP_ML_ML, sc->data) *
                           sc->exp_f(1, i - 1, k, k + 1, VRNA_DECOMP_ML_ML_ML, sc->data);

                  tmp2 += tmp;
                }
              }
            } else {
              if (hc->f) {
                for (k = turn + 2; k < i - turn - 2; k++) {
                  if (hc->f(i, j, i - 1, j + 1, VRNA_DECOMP_PAIR_ML, hc->data) &&
                      hc->f(i, n, i, j, VRNA_DECOMP_ML_ML, hc->data) &&
                      hc->f(1, i - 1, k, k + 1, VRNA_DECOMP_ML_ML_ML, hc->data)) {
                    tmp2 += qm[my_iindx[1] - k] *
                            qm1[jindx[i - 1] + k + 1] *
                            expMLbase[n - j] *
                            expMLclosing *
                            exp_E_MLstem(type,
                                         S1[i - 1],
                                         S1[j + 1],
                                         pf_params);
                  }
                }
              } else {
                for (k = turn + 2; k < i - turn - 2; k++) {
                  tmp2 += qm[my_iindx[1] - k] *
                          qm1[jindx[i - 1] + k + 1] *
                          expMLbase[n - j] *
                          expMLclosing *
                          exp_E_MLstem(type,
                                       S1[i - 1],
                                       S1[j + 1],
                                       pf_params);
                }
              }
            }
          }

          /* 1.3.3 Right part */
          if (hc->up_ml[1] >= (i - 1)) {
            if (sc) {
              if (hc->f) {
                for (k = j + turn + 2; k < n - turn - 1; k++) {
                  if (hc->f(i, j, i - 1, j + 1, VRNA_DECOMP_PAIR_ML, hc->data) &&
                      hc->f(1, j, i, j, VRNA_DECOMP_ML_ML, hc->data) &&
                      hc->f(j + 1, n, k, k + 1, VRNA_DECOMP_ML_ML_ML, hc->data)) {
                    tmp = qm[my_iindx[j + 1] - k] *
                          qm1[jindx[n] + k + 1] *
                          expMLbase[i - 1] *
                          expMLclosing *
                          exp_E_MLstem(type,
                                       S1[i - 1],
                                       S1[j + 1],
                                       pf_params);


                    if (sc->exp_energy_up)
                      tmp *= sc->exp_energy_up[1][i - 1];

                    if (sc->exp_f)
                      tmp *= sc->exp_f(i, j, i - 1, j + 1, VRNA_DECOMP_PAIR_ML, sc->data) *
                             sc->exp_f(1, j, i, j, VRNA_DECOMP_ML_ML, sc->data) *
                             sc->exp_f(j + 1, n, k, k + 1, VRNA_DECOMP_ML_ML_ML, sc->data);

                    tmp2 += tmp;
                  }
                }
              } else {
                for (k = j + turn + 2; k < n - turn - 1; k++) {
                  tmp = qm[my_iindx[j + 1] - k] *
                        qm1[jindx[n] + k + 1] *
                        expMLbase[i - 1] *
                        expMLclosing *
                        exp_E_MLstem(type,
                                     S1[i - 1],
                                     S1[j + 1],
                                     pf_params);


                  if (sc->exp_energy_up)
                    tmp *= sc->exp_energy_up[1][i - 1];

                  if (sc->exp_f)
                    tmp *= sc->exp_f(i, j, i - 1, j + 1, VRNA_DECOMP_PAIR_ML, sc->data) *
                           sc->exp_f(1, j, i, j, VRNA_DECOMP_ML_ML, sc->data) *
                           sc->exp_f(j + 1, n, k, k + 1, VRNA_DECOMP_ML_ML_ML, sc->data);

                  tmp2 += tmp;
                }
              }
            } else {
              if (hc->f) {
                for (k = j + turn + 2; k < n - turn - 1; k++) {
                  if (hc->f(i, j, i - 1, j + 1, VRNA_DECOMP_PAIR_ML, hc->data) &&
                      hc->f(1, j, i, j, VRNA_DECOMP_ML_ML, hc->data) &&
                      hc->f(j + 1, n, k, k + 1, VRNA_DECOMP_ML_ML_ML, hc->data)) {
                    tmp2 += qm[my_iindx[j + 1] - k] *
                            qm1[jindx[n] + k + 1] *
                            expMLbase[i - 1] *
                            expMLclosing *
                            exp_E_MLstem(type,
                                         S1[i - 1],
                                         S1[j + 1],
                                         pf_params);
                  }
                }
              } else {
                for (k = j + turn + 2; k < n - turn - 1; k++) {
                  tmp2 += qm[my_iindx[j + 1] - k] *
                          qm1[jindx[n] + k + 1] *
                          expMLbase[i - 1] *
                          expMLclosing *
                          exp_E_MLstem(type,
                                       S1[i - 1],
                                       S1[j + 1],
                                       pf_params);
                }
              }
            }
          }

          /* all exterior loop decompositions for pair i,j done  */
        }

        probs[ij] *= tmp2;
      } else {
        probs[ij] = 0;
      }
    }
  }
}


PRIVATE INLINE void
bppm_circ_comparative(vrna_fold_compound_t *vc)
{
  unsigned char     *hard_constraints;
  short             **S, **S5, **S3;
  unsigned int      s, n_seq, *type, **a2s;
  int               i, j, k, l, n, ij, turn, *my_iindx, *jindx, *pscore, *rtype;
  FLT_OR_DBL        *qb, *qm, *qm1, *probs, qo, *scale, *expMLbase, expMLclosing, tmp2, tmp3;
  double            kTn;
  vrna_exp_param_t  *pf_params;
  vrna_md_t         *md;
  vrna_hc_t         *hc;
  vrna_sc_t         **scs;

  n                 = (int)vc->length;
  n_seq             = vc->n_seq;
  S                 = vc->S;
  S5                = vc->S5;
  S3                = vc->S3;
  a2s               = vc->a2s;
  my_iindx          = vc->iindx;
  jindx             = vc->jindx;
  pf_params         = vc->exp_params;
  md                = &(pf_params->model_details);
  qb                = vc->exp_matrices->qb;
  qm                = vc->exp_matrices->qm;
  qm1               = vc->exp_matrices->qm1;
  qo                = vc->exp_matrices->qo;
  probs             = vc->exp_matrices->probs;
  scale             = vc->exp_matrices->scale;
  pscore            = vc->pscore;
  expMLbase         = vc->exp_matrices->expMLbase;
  expMLclosing      = pf_params->expMLclosing;
  turn              = md->min_loop_size;
  rtype             = &(md->rtype[0]);
  kTn               = pf_params->kT / 10.;   /* kT in cal/mol  */
  hc                = vc->hc;
  scs               = vc->scs;
  hard_constraints  = hc->mx;

  type = (unsigned int *)vrna_alloc(sizeof(unsigned int) * n_seq);

  for (i = 1; i <= n; i++) {
    for (j = i + turn + 1; j <= n; j++) {
      ij = my_iindx[i] - j;
      if (qb[ij] > 0.) {
        probs[ij] = exp(pscore[jindx[j] + i] / kTn) / qo;

        /* get pair types  */
        for (s = 0; s < n_seq; s++)
          type[s] = vrna_get_ptype_md(S[s][j], S[s][i], md);

        tmp2 = 0.;

        /* 1.1. Exterior Hairpin Contribution */
        tmp2 += vrna_exp_E_hp_loop(vc, j, i);
        /* 1.2. Exterior Interior Loop Contribution */
        /* recycling of k and l... */
        if (hard_constraints[i * n + j] & VRNA_CONSTRAINT_CONTEXT_INT_LOOP) {
          /* 1.2.1. first we calc exterior loop energy with constraint, that i,j  */
          /* delimtis the "right" part of the interior loop                       */
          /* (l,k) is "outer pair"                                                */
          for (k = 1; k < i - turn - 1; k++) {
            /* so first, lets calc the length of loop between j and k */
            int ln1, lstart;
            ln1 = k + n - j - 1;
            if (ln1 > MAXLOOP)
              break;

            if (hc->up_int[j + 1] < ln1)
              break;

            lstart = ln1 + i - 1 - MAXLOOP;
            if (lstart < k + turn + 1)
              lstart = k + turn + 1;

            for (l = lstart; l < i; l++) {
              int ln2, ln2a, ln1a, type_2;
              ln2 = i - l - 1;
              if (ln1 + ln2 > MAXLOOP)
                continue;

              if (hc->up_int[l + 1] < ln2)
                continue;

              if (!(hard_constraints[k * n + l] & VRNA_CONSTRAINT_CONTEXT_INT_LOOP))
                continue;

              FLT_OR_DBL qloop = 1.;
              if (qb[my_iindx[k] - l] == 0.) {
                qloop = 0.;
                continue;
              }

              for (s = 0; s < n_seq; s++) {
                ln2a    = a2s[s][i - 1];
                ln2a    -= a2s[s][l];
                ln1a    = a2s[s][n] - a2s[s][j];
                ln1a    += a2s[s][k - 1];
                type_2  = vrna_get_ptype_md(S[s][l], S[s][k], md);
                qloop   *= exp_E_IntLoop(ln1a, ln2a, type[s], type_2,
                                         S3[s][j],
                                         S5[s][i],
                                         S5[s][k],
                                         S3[s][l], pf_params);
              }
              if (scs) {
                for (s = 0; s < n_seq; s++) {
                  if (scs[s]) {
                    ln2a  = a2s[s][i - 1];
                    ln2a  -= a2s[s][l];
                    ln1a  = a2s[s][n] - a2s[s][j];
                    ln1a  += a2s[s][k - 1];

                    if (scs[s]->exp_energy_up) {
                      qloop *= scs[s]->exp_energy_up[a2s[s][l] + 1][ln2a]
                               * ((j <
                                   n) ? scs[s]->exp_energy_up[a2s[s][j] + 1][a2s[s][n] -
                                                                             a2s[s][j]] : 1.)
                               * ((k > 1) ? scs[s]->exp_energy_up[1][a2s[s][k] - 1] : 1.);
                    }

                    if ((ln1a + ln2a == 0) &&
                        (scs[s]->exp_energy_stack)) {
                      if ((S[s][i]) &&
                          (S[s][j]) &&
                          (S[s][k]) &&
                          (S[s][l])) {
                        /* don't allow gaps in stack */
                        qloop *= scs[s]->exp_energy_stack[a2s[s][k]]
                                 * scs[s]->exp_energy_stack[a2s[s][l]]
                                 * scs[s]->exp_energy_stack[a2s[s][i]]
                                 * scs[s]->exp_energy_stack[a2s[s][j]];
                      }
                    }
                  }
                }
              }

              tmp2 += qb[my_iindx[k] - l] * qloop * scale[ln1 + ln2];
            }
          }

          /* 1.2.2. second we calc exterior loop energy with constraint, that i,j */
          /* delimtis the "left" part of the interior loop                        */
          /* (j,i) is "outer pair"                                                */
          for (k = j + 1; k < n - turn; k++) {
            /* so first, lets calc the length of loop between l and i */
            int ln1, lstart;
            ln1 = k - j - 1;
            if ((ln1 + i - 1) > MAXLOOP)
              break;

            if (hc->up_int[j + 1] < ln1)
              break;

            lstart = ln1 + i - 1 + n - MAXLOOP;
            if (lstart < k + turn + 1)
              lstart = k + turn + 1;

            for (l = lstart; l <= n; l++) {
              int ln2, ln2a, ln1a, type_2;
              ln2 = i - 1 + n - l;
              if (ln1 + ln2 > MAXLOOP)
                continue;

              if (hc->up_int[l + 1] < ln2)
                continue;

              if (!(hard_constraints[k * n + l] & VRNA_CONSTRAINT_CONTEXT_INT_LOOP))
                continue;

              FLT_OR_DBL qloop = 1.;
              if (qb[my_iindx[k] - l] == 0.) {
                qloop = 0.;
                continue;
              }

              for (s = 0; s < n_seq; s++) {
                ln1a    = a2s[s][k] - a2s[s][j + 1];
                ln2a    = a2s[s][i - 1] + a2s[s][n] - a2s[s][l];
                type_2  = vrna_get_ptype_md(S[s][l], S[s][k], md);
                qloop   *= exp_E_IntLoop(ln2a, ln1a, type_2, type[s],
                                         S3[s][l],
                                         S5[s][k],
                                         S5[s][i],
                                         S3[s][j], pf_params);
              }
              if (scs) {
                for (s = 0; s < n_seq; s++) {
                  if (scs[s]) {
                    ln1a  = a2s[s][k] - a2s[s][j + 1];
                    ln2a  = a2s[s][i - 1] + a2s[s][n] - a2s[s][l];

                    if (scs[s]->exp_energy_up) {
                      qloop *= scs[s]->exp_energy_up[a2s[s][j] + 1][ln1a]
                               * ((l <
                                   n) ? scs[s]->exp_energy_up[a2s[s][l] + 1][a2s[s][n] -
                                                                             a2s[s][l]] : 1.)
                               * ((i > 1) ? scs[s]->exp_energy_up[1][a2s[s][i] - 1] : 1.);
                    }

                    if ((ln1a + ln2a == 0) &&
                        (scs[s]->exp_energy_stack)) {
                      if ((S[s][i]) &&
                          (S[s][j]) &&
                          (S[s][k]) &&
                          (S[s][l])) {
                        /* don't allow gaps in stack */
                        qloop *= scs[s]->exp_energy_stack[a2s[s][k]]
                                 * scs[s]->exp_energy_stack[a2s[s][l]]
                                 * scs[s]->exp_energy_stack[a2s[s][i]]
                                 * scs[s]->exp_energy_stack[a2s[s][j]];
                      }
                    }
                  }
                }
              }

              tmp2 += qb[my_iindx[k] - l] * qloop * scale[ln1 + ln2];
            }
          }
        }

        /* 1.3 Exterior multiloop decomposition */
        if (hard_constraints[i * n + j] & VRNA_CONSTRAINT_CONTEXT_MB_LOOP) {
          /* 1.3.1 Middle part                    */
          if ((i > turn + 2) &&
              (j < n - turn - 1)) {
            for (tmp3 = 1, s = 0; s < n_seq; s++)
              tmp3 *= exp_E_MLstem(rtype[type[s]], S5[s][i], S3[s][j], pf_params);
            tmp2 += qm[my_iindx[1] - i + 1] * qm[my_iindx[j + 1] - n] *tmp3 *pow(expMLclosing,
                                                                                 n_seq);
          }

          /* 1.3.2 Left part    */
          for (k = turn + 2; k < i - turn - 2; k++) {
            if (hc->up_ml[j + 1] < n - j)
              break;

            for (tmp3 = 1, s = 0; s < n_seq; s++)
              tmp3 *= exp_E_MLstem(rtype[type[s]], S5[s][i], S3[s][j], pf_params);

            if (scs) {
              for (s = 0; s < n_seq; s++) {
                if (scs[s]) {
                  if (scs[s]->exp_energy_bp)
                    tmp3 *= scs[s]->exp_energy_bp[jindx[j] + i];

                  if (scs[s]->exp_energy_up)
                    tmp3 *= scs[s]->exp_energy_up[a2s[s][j] + 1][a2s[s][n] - a2s[s][j]];
                }
              }
            }

            tmp2 += tmp3 *
                    qm[my_iindx[1] - k] *
                    qm1[jindx[i - 1] + k + 1] *
                    expMLbase[n - j] *
                    pow(expMLclosing, n_seq);
          }
          /* 1.3.3 Right part    */
          for (k = j + turn + 2; k < n - turn - 1; k++) {
            if (hc->up_ml[1] < i - 1)
              break;

            for (tmp3 = 1, s = 0; s < n_seq; s++)
              tmp3 *= exp_E_MLstem(rtype[type[s]], S5[s][i], S3[s][j], pf_params);

            if (scs) {
              for (s = 0; s < n_seq; s++) {
                if (scs[s]) {
                  if (scs[s]->exp_energy_bp)
                    tmp3 *= scs[s]->exp_energy_bp[jindx[j] + i];

                  if (scs[s]->exp_energy_up)
                    tmp3 *= scs[s]->exp_energy_up[a2s[s][1]][a2s[s][i] - a2s[s][1]];
                }
              }
            }

            tmp2 += qm[my_iindx[j + 1] - k] * qm1[jindx[n] + k + 1] * tmp3 * expMLbase[i - 1] *
                    pow(expMLclosing, n_seq);
          }
        }

        probs[ij] *= tmp2;
      } else {
        probs[ij] = 0;
      }
    } /* end for j=..*/
  }   /* end or i=...  */

  free(type);
}


PRIVATE double
wrap_mean_bp_distance(FLT_OR_DBL  *p,
                      int         length,
                      int         *index,
                      int         turn)
{
  int     i, j;
  double  d = 0.;

  /* compute the mean base pair distance in the thermodynamic ensemble */
  /* <d> = \sum_{a,b} p_a p_b d(S_a,S_b)
   * this can be computed from the pair probs p_ij as
   * <d> = \sum_{ij} p_{ij}(1-p_{ij}) */

  for (i = 1; i <= length; i++)
    for (j = i + turn + 1; j <= length; j++)
      d += p[index[i] - j] * (1 - p[index[i] - j]);

  return 2 * d;
}
