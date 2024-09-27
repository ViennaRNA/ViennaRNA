#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <string.h>
#include <limits.h>

#include "ViennaRNA/fold_vars.h"
#include "ViennaRNA/params/default.h"
#include "ViennaRNA/utils/basic.h"
#include "ViennaRNA/utils/log.h"
#include "ViennaRNA/sequences/alphabet.h"
#include "ViennaRNA/constraints/hard.h"
#include "ViennaRNA/constraints/soft.h"
#include "ViennaRNA/eval/gquad.h"
#include "ViennaRNA/structured_domains.h"
#include "ViennaRNA/unstructured_domains.h"
#include "ViennaRNA/eval/exterior.h"
#include "ViennaRNA/partfunc/exterior.h"

#include "ViennaRNA/intern/grammar_dat.h"

#ifdef __GNUC__
# define INLINE inline
#else
# define INLINE
#endif

#define SPEEDUP_HC  0

#include "ViennaRNA/constraints/exterior_hc.inc"
#include "ViennaRNA/constraints/exterior_sc_pf.inc"

struct vrna_mx_pf_aux_el_s {
  FLT_OR_DBL  *qq;
  FLT_OR_DBL  *qq1;

  int         qqu_size;
  FLT_OR_DBL  **qqu;
};

/*
 #################################
 # PRIVATE FUNCTION DECLARATIONS #
 #################################
 */

PRIVATE INLINE FLT_OR_DBL
reduce_ext_ext_fast(vrna_fold_compound_t        *fc,
                    unsigned int                i,
                    unsigned int                j,
                    struct vrna_mx_pf_aux_el_s  *aux_mx,
                    vrna_hc_eval_f   evaluate,
                    struct hc_ext_def_dat       *hc_dat_local,
                    struct sc_ext_exp_dat       *sc_wrapper);


PRIVATE INLINE FLT_OR_DBL
reduce_ext_stem_fast(vrna_fold_compound_t       *fc,
                     unsigned int               i,
                     unsigned int               j,
                     vrna_hc_eval_f             evaluate,
                     struct hc_ext_def_dat      *hc_dat_local,
                     struct sc_ext_exp_dat      *sc_wrapper);


PRIVATE INLINE FLT_OR_DBL
reduce_ext_up_fast(vrna_fold_compound_t       *fc,
                   unsigned int               i,
                   unsigned int               j,
                   vrna_hc_eval_f             evaluate,
                   struct hc_ext_def_dat      *hc_dat_local,
                   struct sc_ext_exp_dat      *sc_wrapper);


PRIVATE INLINE FLT_OR_DBL
split_ext_fast(vrna_fold_compound_t       *fc,
               unsigned int               i,
               unsigned int               j,
               struct vrna_mx_pf_aux_el_s *aux_mx,
               vrna_hc_eval_f  evaluate,
               struct hc_ext_def_dat      *hc_dat_local,
               struct sc_ext_exp_dat      *sc_wrapper);


PRIVATE FLT_OR_DBL
exp_E_ext_fast(vrna_fold_compound_t       *fc,
               unsigned int               i,
               unsigned int               j,
               struct vrna_mx_pf_aux_el_s *aux_mx);


/*
 #################################
 # BEGIN OF FUNCTION DEFINITIONS #
 #################################
 */
PUBLIC struct vrna_mx_pf_aux_el_s *
vrna_exp_E_ext_fast_init(vrna_fold_compound_t *fc)
{
  struct vrna_mx_pf_aux_el_s *aux_mx = NULL;

  if (fc) {
    unsigned int              u, i, j, max_j, d, n, turn, with_ud;
    int                       ij, *iidx;
    FLT_OR_DBL                *q, **q_local;
    vrna_hc_eval_f evaluate;
    struct hc_ext_def_dat     hc_dat_local;
    struct sc_ext_exp_dat     sc_wrapper;
    vrna_ud_t                 *domains_up;

    n           = fc->length;
    iidx        = fc->iindx;
    turn        = fc->exp_params->model_details.min_loop_size;
    domains_up  = fc->domains_up;
    with_ud     = (domains_up && domains_up->exp_energy_cb);

    if (fc->hc->type == VRNA_HC_WINDOW)
      evaluate = prepare_hc_ext_def_window(fc, &hc_dat_local);
    else
      evaluate = prepare_hc_ext_def(fc, &hc_dat_local);

    init_sc_ext_exp(fc, &sc_wrapper);

    /* allocate memory for helper arrays */
    aux_mx            = (struct vrna_mx_pf_aux_el_s *)vrna_alloc(sizeof(struct vrna_mx_pf_aux_el_s));
    aux_mx->qq        = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * (n + 2));
    aux_mx->qq1       = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * (n + 2));
    aux_mx->qqu_size  = 0;
    aux_mx->qqu       = NULL;

    /* pre-processing ligand binding production rule(s) and auxiliary memory */
    if (with_ud) {
      unsigned int ud_max_size = 0;
      for (u = 0; u < (unsigned int)domains_up->uniq_motif_count; u++)
        if (ud_max_size < domains_up->uniq_motif_size[u])
          ud_max_size = domains_up->uniq_motif_size[u];

      aux_mx->qqu_size  = ud_max_size;
      aux_mx->qqu       = (FLT_OR_DBL **)vrna_alloc(sizeof(FLT_OR_DBL *) * (ud_max_size + 1));

      for (u = 0; u <= ud_max_size; u++)
        aux_mx->qqu[u] = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * (n + 2));
    }

    if (fc->hc->type == VRNA_HC_WINDOW) {
      q_local = fc->exp_matrices->q_local;
      max_j   = MIN2(turn + 1, (unsigned int)fc->window_size);
      max_j   = MIN2(max_j, n);
      for (j = 1; j <= max_j; j++)
        for (i = 1; i <= j; i++)
          q_local[i][j] = reduce_ext_up_fast(fc,
                                             i,
                                             j,
                                             evaluate,
                                             &hc_dat_local,
                                             &sc_wrapper);
    } else {
      q = fc->exp_matrices->q;
      for (d = 0; d <= MIN2(turn, n); d++)
        for (i = 1; i <= n - d; i++) {
          j   = i + d;
          ij  = iidx[i] - j;

          q[ij] = reduce_ext_up_fast(fc,
                                     i,
                                     j,
                                     evaluate,
                                     &hc_dat_local,
                                     &sc_wrapper);
        }

      if (fc->aux_grammar) {
        for (size_t c = 0; c < vrna_array_size(fc->aux_grammar->exp_f); c++) {
          if (fc->aux_grammar->exp_f[c].cb) {
            for (d = 0; d <= MIN2(turn, n); d++)
              for (i = 1; i <= n - d; i++) {
                j   = i + d;
                ij  = iidx[i] - j;

                q[ij] += fc->aux_grammar->exp_f[c].cb(fc,
                                                      i,
                                                      j,
                                                      fc->aux_grammar->exp_f[c].data);
              }
          }
        }
      }
    }
  }

  return aux_mx;
}


PUBLIC void
vrna_exp_E_ext_fast_rotate(struct vrna_mx_pf_aux_el_s *aux_mx)
{
  if (aux_mx) {
    unsigned int  u;
    FLT_OR_DBL    *tmp;

    tmp         = aux_mx->qq1;
    aux_mx->qq1 = aux_mx->qq;
    aux_mx->qq  = tmp;

    /* rotate auxiliary arrays for unstructured domains */
    if (aux_mx->qqu) {
      tmp = aux_mx->qqu[aux_mx->qqu_size];
      for (u = aux_mx->qqu_size; u > 0; u--)
        aux_mx->qqu[u] = aux_mx->qqu[u - 1];
      aux_mx->qqu[0] = tmp;
    }
  }
}


PUBLIC void
vrna_exp_E_ext_fast_free(struct vrna_mx_pf_aux_el_s *aux_mx)
{
  if (aux_mx) {
    int u;

    free(aux_mx->qq);
    free(aux_mx->qq1);

    if (aux_mx->qqu) {
      for (u = 0; u <= aux_mx->qqu_size; u++)
        free(aux_mx->qqu[u]);

      free(aux_mx->qqu);
    }

    free(aux_mx);
  }
}


PUBLIC FLT_OR_DBL
vrna_exp_E_ext_fast(vrna_fold_compound_t        *fc,
                    int                         i,
                    int                         j,
                    struct vrna_mx_pf_aux_el_s  *aux_mx)
{
  if (fc) {
    if (j < i) {
      int t = j;
      vrna_log_warning(
        "vrna_exp_E_ext_fast: i (%d) larger than j (%d)! Swapping coordinates...",
        i,
        j);
      j = i;
      i = t;
    } else if ((j < 1) || (i < 1)) {
      vrna_log_warning(
        "vrna_exp_E_ext_fast: Indices too small [i = %d, j = %d]! "
        "Refusing to compute anything...",
        i,
        j);
      return 0.;
    } else if (j > (int)fc->length) {
      vrna_log_warning(
        "vrna_exp_E_ext_fast: Indices exceed sequence length (%d) [i = %d, j = %d]! "
        "Refusing to compute anything...",
        fc->length,
        i,
        j);
      return 0.;
    }

    return exp_E_ext_fast(fc, i, j, aux_mx);
  }

  return 0.;
}


PUBLIC void
vrna_exp_E_ext_fast_update(vrna_fold_compound_t       *fc,
                           int                        j,
                           struct vrna_mx_pf_aux_el_s *aux_mx VRNA_UNUSED)
{
  int                       k;
  FLT_OR_DBL                **q;
  vrna_hc_eval_f evaluate;
  struct hc_ext_def_dat     hc_dat_local;
  struct sc_ext_exp_dat     sc_wrapper;

  /*
   *  init exterior loop contributions for small segments [i, j]
   *  that can only be unpaired. This needs to be done for each
   *  j just before any contributions for [i,j] will be computed
   */
  if ((fc) && (fc->hc->type == VRNA_HC_WINDOW)) {
    q         = fc->exp_matrices->q_local;
    evaluate  = prepare_hc_ext_def_window(fc, &hc_dat_local);
    init_sc_ext_exp(fc, &sc_wrapper);


    for (k = j; k >= MAX2(1, j); k--)
      q[k][j] = reduce_ext_up_fast(fc, k, j, evaluate, &hc_dat_local, &sc_wrapper);
  }
}


PRIVATE INLINE FLT_OR_DBL
reduce_ext_ext_fast(vrna_fold_compound_t        *fc,
                    unsigned int                i,
                    unsigned int                j,
                    struct vrna_mx_pf_aux_el_s  *aux_mx,
                    vrna_hc_eval_f   evaluate,
                    struct hc_ext_def_dat       *hc_dat_local,
                    struct sc_ext_exp_dat       *sc_wrapper)
{
  unsigned int  u, cnt;
  FLT_OR_DBL    q_temp, q_temp2, q, *qq1, **qqu, *scale;
  vrna_ud_t     *domains_up;
  sc_ext_exp_cb sc_red_ext;

  domains_up  = fc->domains_up;
  qq1         = aux_mx->qq1;
  qqu         = aux_mx->qqu;
  scale       = fc->exp_matrices->scale;
  sc_red_ext  = sc_wrapper->red_ext;

  q = 0.;

  if (evaluate(i, j, i, j - 1, VRNA_DECOMP_EXT_EXT, hc_dat_local)) {
    q_temp = qq1[i] * scale[1];

    if (sc_red_ext)
      q_temp *= sc_red_ext(i, j, i, j - 1, sc_wrapper);

    if ((domains_up) && (domains_up->exp_energy_cb)) {
      for (cnt = 0; cnt < (unsigned int)domains_up->uniq_motif_count; cnt++) {
        u = domains_up->uniq_motif_size[cnt];
        if (j - u >= i) {
          if (evaluate(i, j, i, j - u, VRNA_DECOMP_EXT_EXT, hc_dat_local)) {
            q_temp2 = qqu[u][i] *
                      domains_up->exp_energy_cb(fc,
                                                j - u + 1,
                                                j,
                                                VRNA_UNSTRUCTURED_DOMAIN_EXT_LOOP | VRNA_UNSTRUCTURED_DOMAIN_MOTIF,
                                                domains_up->data) *
                      scale[u];

            if (sc_red_ext)
              q_temp2 *= sc_red_ext(i, j, i, j - u, sc_wrapper);

            q_temp += q_temp2;
          }
        }
      }
    }

    q = q_temp;
  }

  return q;
}


PRIVATE INLINE FLT_OR_DBL
reduce_ext_stem_fast(vrna_fold_compound_t       *fc,
                     unsigned int               i,
                     unsigned int               j,
                     vrna_hc_eval_f             evaluate,
                     struct hc_ext_def_dat      *hc_dat_local,
                     struct sc_ext_exp_dat      *sc_wrapper)
{
  short             **S, **S5, **S3, *S1, *S2, s5, s3;
  unsigned int      type, *sn, n, s, n_seq, **a2s, circular;
  int               *idx;
  FLT_OR_DBL        qbt, q_temp, qb;
  vrna_exp_param_t  *pf_params;
  vrna_md_t         *md;
  sc_ext_exp_cb     sc_red_stem;

  sc_red_stem = sc_wrapper->red_stem;
  n           = fc->length;
  sn          = fc->strand_number;
  pf_params   = fc->exp_params;
  md          = &(pf_params->model_details);
  circular    = md->circ;
  idx         = fc->iindx;
  qb          = (fc->hc->type == VRNA_HC_WINDOW) ?
                fc->exp_matrices->qb_local[i][j] :
                fc->exp_matrices->qb[idx[i] - j];
  qbt = 0.;

  /* exterior loop part with stem (i, j) */
  if (evaluate(i, j, i, j, VRNA_DECOMP_EXT_STEM, hc_dat_local)) {
    q_temp = qb;

    switch (fc->type) {
      case VRNA_FC_TYPE_SINGLE:
        S1      = fc->sequence_encoding;
        S2      = fc->sequence_encoding2;
        type    = vrna_get_ptype_md(S2[i], S2[j], md);
        s5      = (((i > 1) || circular) && (sn[i] == sn[i - 1])) ? S1[i - 1] : -1;
        s3      = (((j < n) || circular) && (sn[j + 1] == sn[j])) ? S1[j + 1] : -1;
        q_temp  *= vrna_exp_E_exterior_stem(type, s5, s3, pf_params);
        break;

      case VRNA_FC_TYPE_COMPARATIVE:
        n_seq = fc->n_seq;
        S     = fc->S;
        S5    = fc->S5;
        S3    = fc->S3;
        a2s   = fc->a2s;
        for (s = 0; s < n_seq; s++) {
          type    = vrna_get_ptype_md(S[s][i], S[s][j], md);
          q_temp  *= vrna_exp_E_exterior_stem(type,
                                         ((a2s[s][i] > 1) || circular) ? S5[s][i] : -1,
                                         ((a2s[s][j] < a2s[s][n]) || circular) ? S3[s][j] : -1,
                                         pf_params);
        }
        break;
    }

    if (sc_red_stem)
      q_temp *= sc_red_stem(i, j, i, j, sc_wrapper);

    qbt += q_temp;
  }

  return qbt;
}


PRIVATE INLINE FLT_OR_DBL
reduce_ext_up_fast(vrna_fold_compound_t       *fc,
                   unsigned int               i,
                   unsigned int               j,
                   vrna_hc_eval_f             evaluate,
                   struct hc_ext_def_dat      *hc_dat_local,
                   struct sc_ext_exp_dat      *sc_wrapper)
{
  unsigned int      u;
  FLT_OR_DBL        qbt, q_temp, *scale;
  vrna_ud_t         *domains_up;
  sc_ext_exp_red_up sc_red_up;

  sc_red_up = sc_wrapper->red_up;

  scale       = fc->exp_matrices->scale;
  domains_up  = fc->domains_up;
  qbt         = 0.;

  if (evaluate(i, j, i, j, VRNA_DECOMP_EXT_UP, hc_dat_local)) {
    u       = j - i + 1;
    q_temp  = scale[u];

    if (sc_red_up)
      q_temp *= sc_red_up(i, j, sc_wrapper);

    qbt += q_temp;

    if ((domains_up) && (domains_up->exp_energy_cb)) {
      qbt += q_temp *
             domains_up->exp_energy_cb(fc,
                                       i, j,
                                       VRNA_UNSTRUCTURED_DOMAIN_EXT_LOOP,
                                       domains_up->data);
    }
  }

  return qbt;
}


PRIVATE INLINE FLT_OR_DBL
split_ext_fast(vrna_fold_compound_t       *fc,
               unsigned int               i,
               unsigned int               j,
               struct vrna_mx_pf_aux_el_s *aux_mx,
               vrna_hc_eval_f  evaluate,
               struct hc_ext_def_dat      *hc_dat_local,
               struct sc_ext_exp_dat      *sc_wrapper)
{
  unsigned int      k;
  int               *idx, ij1;
  FLT_OR_DBL        qbt, *q, *qq, *qqq;
  sc_ext_exp_split  sc_split;

  sc_split = sc_wrapper->split;

  idx = fc->iindx;
  q   = (fc->hc->type == VRNA_HC_WINDOW) ?
        fc->exp_matrices->q_local[i] :
        fc->exp_matrices->q + idx[i];
  qq  = aux_mx->qq;
  qbt = 0.;

  /*
   *  use an auxiliary array qqq that already contains soft constraint
   *  contributions when we split the exterior loops
   */
  if (sc_split) {
    /* pre-process qq array */
    qqq = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * (j - i + 1));
    qqq -= i;

    for (k = j; k > i; k--)
      qqq[k] = qq[k] * sc_split(i, j, k, sc_wrapper);
  } else {
    /* pre-process qq array */
    qqq = qq;
  }

  /*
   *  the index for access in q array now is:
   *
   *  (a) q[- (j - 1)] for global, and
   *  (b) q[j - 1] for local structure prediction
   *
   *  Thus, we may use a single variable to address both cases:
   *  k = j:
   *    ij1 = (j - 1) * factor = j * factor - factor
   *  k = j - 1:
   *    ij1 = (j - 2) * factor = j * factor - 2 * factor = j * factor - factor - factor
   *  ...
   */
  int factor = (fc->hc->type == VRNA_HC_WINDOW) ? 1 : -1;

  ij1 = factor * ((int)j - 1);

  /* do actual decomposition (skip hard constraint checks if we use default settings) */
#if SPEEDUP_HC
  /*
   *  checking whether we actually are provided with hard constraints and
   *  otherwise not evaluating the default ones within the loop drastically
   *  increases speed. However, once we check for the split point between
   *  strands in hard constraints, we have to think of something else...
   */
  if ((evaluate == &hc_ext_cb_def) || (evaluate == &hc_ext_cb_def_window)) {
    for (k = j; k > i; k--) {
      qbt += q[ij1] *
             qqq[k];
      ij1 -= factor;
    }
  } else {
    for (k = j; k > i; k--) {
      if (evaluate(i, j, k - 1, k, VRNA_DECOMP_EXT_EXT_EXT, hc_dat_local))
        qbt += q[ij1] *
               qqq[k];

      ij1 -= factor;
    }
  }

#else
  for (k = j; k > i; k--) {
    if (evaluate(i, j, k - 1, k, VRNA_DECOMP_EXT_EXT_EXT, hc_dat_local))
      qbt += q[ij1] *
             qqq[k];

    ij1 -= factor;
  }
#endif

  if (qqq != qq) {
    qqq += i;
    free(qqq);
  }

  return qbt;
}


PRIVATE FLT_OR_DBL
exp_E_ext_fast(vrna_fold_compound_t       *fc,
               unsigned int               i,
               unsigned int               j,
               struct vrna_mx_pf_aux_el_s *aux_mx)
{
  unsigned int              with_ud, with_gquad;
  FLT_OR_DBL                qbt1, *qq, **qqu, **G_local;
  vrna_md_t                 *md;
  vrna_exp_param_t          *pf_params;
  vrna_ud_t                 *domains_up;
  vrna_hc_eval_f evaluate;
  struct hc_ext_def_dat     hc_dat_local;
  struct sc_ext_exp_dat     sc_wrapper;
  vrna_smx_csr(FLT_OR_DBL)  *q_gq;

  qq          = aux_mx->qq;
  qqu         = aux_mx->qqu;
  pf_params   = fc->exp_params;
  md          = &(pf_params->model_details);
  domains_up  = fc->domains_up;
  with_gquad  = md->gquad;
  with_ud     = (domains_up && domains_up->exp_energy_cb);

  if (fc->hc->type == VRNA_HC_WINDOW)
    evaluate = prepare_hc_ext_def_window(fc, &hc_dat_local);
  else
    evaluate = prepare_hc_ext_def(fc, &hc_dat_local);

  init_sc_ext_exp(fc, &sc_wrapper);

  qbt1 = 0.;

  /* all exterior loop parts [i, j] with exactly one stem (i, u) i < u < j */
  qbt1 += reduce_ext_ext_fast(fc, i, j, aux_mx, evaluate, &hc_dat_local, &sc_wrapper);
  /* exterior loop part with stem (i, j) */
  qbt1 += reduce_ext_stem_fast(fc, i, j, evaluate, &hc_dat_local, &sc_wrapper);

  if (with_gquad) {
    if (fc->hc->type == VRNA_HC_WINDOW) {
      G_local = fc->exp_matrices->G_local;
      qbt1    += G_local[i][j];
    } else {
      q_gq  = fc->exp_matrices->q_gq;
#ifndef VRNA_DISABLE_C11_FEATURES
      qbt1  += vrna_smx_csr_get(q_gq, i, j, 0.);
#else
      qbt1  += vrna_smx_csr_FLT_OR_DBL_get(q_gq, i, j, 0.);
#endif
    }
  }

  qq[i] = qbt1;

  if (with_ud)
    qqu[0][i] = qbt1;

  /* the entire stretch [i,j] is unpaired */
  qbt1 += reduce_ext_up_fast(fc, i, j, evaluate, &hc_dat_local, &sc_wrapper);

  qbt1 += split_ext_fast(fc, i, j, aux_mx, evaluate, &hc_dat_local, &sc_wrapper);

  /* apply auxiliary grammar rule for exterior loop case */
  if (fc->aux_grammar) {
    for (size_t c = 0; c < vrna_array_size(fc->aux_grammar->exp_f); c++) {
      if (fc->aux_grammar->exp_f[c].cb)
        qbt1 += fc->aux_grammar->exp_f[c].cb(fc, i, j, fc->aux_grammar->exp_f[c].data);
    }
  }

  free_sc_ext_exp(&sc_wrapper);

  return qbt1;
}
