/* Crosslinking data handling */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h> /* for uint8_t */
#include <math.h>
#include <ctype.h>
#include <string.h>
#include <limits.h>

#include "ViennaRNA/params/default.h"
#include "ViennaRNA/params/constants.h" /* defines MINPSCORE */
#include "ViennaRNA/datastructures/array.h"
#include "ViennaRNA/fold_vars.h"
#include "ViennaRNA/utils/basic.h"
#include "ViennaRNA/utils/strings.h"
#include "ViennaRNA/sequences/alignments.h"
#include "ViennaRNA/sequences/alphabet.h"
#include "ViennaRNA/utils/log.h"
#include "ViennaRNA/io/utils.h"
#include "ViennaRNA/io/file_formats.h"
#include "ViennaRNA/params/basic.h"
#include "ViennaRNA/constraints/soft.h"
#include "ViennaRNA/probing/basic.h"
#include "ViennaRNA/probabilities/basepairs.h"


#include "ViennaRNA/eval/multibranch.h"
#include "ViennaRNA/backtrack/multibranch.h"

#ifdef __GNUC__
# define INLINE inline
#else
# define INLINE
#endif

#include "ViennaRNA/constraints/exterior_sc_pf.inc"
#include "ViennaRNA/constraints/hairpin_sc_pf.inc"
#include "ViennaRNA/constraints/internal_sc_pf.inc"
#include "ViennaRNA/constraints/multibranch_sc.inc"
#include "ViennaRNA/constraints/multibranch_sc_pf.inc"

#include "ViennaRNA/probabilities/eqm_probs_shared.h"

#include "sc_iip_shared.h"

#include "ViennaRNA/constraints/soft_iip.h"

/* re-definition of vrna_mx_mfe_aux_ml_s to use its internals, see also mfe_multibrnch.c  */
struct vrna_mx_mfe_aux_ml_s {
  unsigned int  n;
  int           *Fmi;
  int           *DMLi;
  int           *DMLi1;
  int           *DMLi2;
};


struct vrna_mx_pf_aux_ml_s {
  FLT_OR_DBL    *qqm;
  FLT_OR_DBL    *qqm1;
  FLT_OR_DBL    *qqm2;
  FLT_OR_DBL    *qqm21;

  unsigned int  qqmu_size;
  FLT_OR_DBL    **qqmu;
};


typedef struct {
  FLT_OR_DBL  *vl;  /* V[l,j] */
  FLT_OR_DBL  *wl;  /* W[l,v] */
  FLT_OR_DBL  *yl;  /* Y[l,k] */
  FLT_OR_DBL  *tl;  /* T[l,k] */
  FLT_OR_DBL  *ol;
  FLT_OR_DBL  *ol2;
  FLT_OR_DBL  *pl;
  FLT_OR_DBL  *pl2;
  FLT_OR_DBL  *rl;
  FLT_OR_DBL  *sl;
  FLT_OR_DBL  *prm_lmm;
  FLT_OR_DBL  *prm_lmm2;
  FLT_OR_DBL  *prm_lmu;
  FLT_OR_DBL  *prm_lmu2;
  FLT_OR_DBL  *prm_l1mu;
  FLT_OR_DBL  *prm_l1mu2;
} xl_ml_helper_t;


PRIVATE INLINE int
ml_pair_d2(vrna_fold_compound_t         *fc,
           unsigned int                 i,
           unsigned int                 j,
           interact_region_s            *r,
           struct vrna_mx_mfe_aux_ml_s  *helpers,
           struct sc_mb_dat             *sc_wrapper);


PRIVATE INLINE int
ml_pair_d0(vrna_fold_compound_t         *fc,
           unsigned int                 i,
           unsigned int                 j,
           interact_region_s            *r,
           struct vrna_mx_mfe_aux_ml_s  *helpers,
           struct sc_mb_dat             *sc_wrapper);


PRIVATE INLINE int
mfe_M2_xlink(vrna_fold_compound_t         *fc,
             unsigned int                 i,
             unsigned int                 j,
             interact_region_s            *r,
             struct vrna_mx_mfe_aux_ml_s  *helpers,
             struct sc_mb_dat             *sc_wrapper);


PRIVATE void
fill_helper_arrays(vrna_fold_compound_t *fc,
                   unsigned int         i,
                   unsigned int         j,
                   struct sc_mb_dat     *sc_wrapper);


PRIVATE unsigned int
bt_mb_loop_xl(vrna_fold_compound_t  *fc,
              unsigned int          i,
              unsigned int          j,
              int                   en,
              vrna_bps_t bp_stack   VRNA_UNUSED,
              vrna_bts_t            bt_stack);


PRIVATE unsigned int
bt_m2(vrna_fold_compound_t  *fc,
      unsigned int          i,
      unsigned int          j,
      int                   en,
      vrna_bps_t bp_stack   VRNA_UNUSED,
      vrna_bts_t            bt_stack);


PRIVATE unsigned int
bt_A(vrna_fold_compound_t *fc,
     int                  *A,
     unsigned int         i,
     unsigned int         j,
     vrna_bps_t bp_stack  VRNA_UNUSED,
     vrna_bts_t           bt_stack);


PRIVATE unsigned int
bt_B(vrna_fold_compound_t *fc,
     int                  *B,
     unsigned int         i,
     unsigned int         j,
     vrna_bps_t bp_stack  VRNA_UNUSED,
     vrna_bts_t           bt_stack);



PRIVATE INLINE FLT_OR_DBL
ml_pair_exp(vrna_fold_compound_t        *fc,
            unsigned int                i,
            unsigned int                j,
            interact_region_s           *r,
            struct vrna_mx_pf_aux_ml_s  *helpers);


PRIVATE INLINE FLT_OR_DBL
pf_M2_xlink(vrna_fold_compound_t        *fc,
            unsigned int                i,
            unsigned int                j,
            interact_region_s           *r,
            struct vrna_mx_pf_aux_ml_s  *helpers,
            struct sc_mb_exp_dat        *sc_wrapper);


PRIVATE void
xl_fill_prms(vrna_fold_compound_t *fc,
             unsigned int         l,
             helper_arrays_t      eqm_prob_helpers,
             xl_ml_helper_t       *helpers);


PRIVATE void
xl_fill_rs(vrna_fold_compound_t *fc,
           unsigned int         l,
           helper_arrays_t      eqm_prob_helpers,
           xl_ml_helper_t       *helpers);


PRIVATE void
xl_fill_s(vrna_fold_compound_t *fc,
          unsigned int         l,
          helper_arrays_t      eqm_prob_helpers,
          xl_ml_helper_t       *helpers);


PRIVATE void
compute_bpp_multibranch_xlink(vrna_fold_compound_t *fc,
                              unsigned int         l,
                              helper_arrays_t      eqm_prob_helpers,
                              void                 *data);


PRIVATE xl_ml_helper_t *
get_xl_ml_outside_helpers(vrna_fold_compound_t *fc);


PRIVATE void
free_xl_ml_outside_helpers(void *helpers);


PUBLIC int
vrna_sc_add_iip(vrna_fold_compound_t  *fc,
                unsigned int          i,
                unsigned int          j,
                unsigned int          k,
                unsigned int          l,
                double                energy,
                unsigned int          options)
{
  int                           ret;
  FILE                          *fp;
  struct vrna_crosslink_data_s  *data;

  ret = 0;

  if ((fc) &&
      (fc->type == VRNA_FC_TYPE_SINGLE)) {

    if (fc->sc == NULL)
      vrna_sc_init(fc);

    if (fc->sc->iip == NULL) {
      data  = (struct vrna_crosslink_data_s *)vrna_alloc(sizeof(struct vrna_crosslink_data_s));
      data->overlap5 = NULL;
      data->overlap3 = NULL;
      vrna_array_init_size(data->regions, 42);

      fc->sc->iip = data;

      fc->params->model_details.uniq_ML = 1;

      if (fc->exp_params)
        fc->exp_params->model_details.uniq_ML = 1;
    }

    ret   = 1;
    data  = fc->sc->iip;

    vrna_array_append(data->regions, ((interact_region_s){
      .a      = i,
      .b      = j,
      .c      = k,
      .d      = l,
      .score  = (int)(energy * 100.),
      .energy = energy
    }));
  }

  return ret;
}


PUBLIC int
vrna_sc_iip_prepare(vrna_fold_compound_t *fc)
{
  interact_region_s *r, *other;

  if ((fc) &&
      (fc->sc) &&
      (fc->sc->iip)) {
    struct vrna_crosslink_data_s *dat = fc->sc->iip;

    /* free previously allocated overlap entries */
    if (dat->overlap5) {
      for (size_t cnt = 0; cnt < vrna_array_size(dat->overlap5); cnt++)
        vrna_array_free(dat->overlap5[cnt]);
      vrna_array_free(dat->overlap5);
    }

    if (dat->overlap3) {
      for (size_t cnt = 0; cnt < vrna_array_size(dat->overlap3); cnt++)
        vrna_array_free(dat->overlap3[cnt]);
      vrna_array_free(dat->overlap3);
    }

    /* initialize memory for overlap entries */
    vrna_array_init_size(dat->overlap5, fc->length + 1);
    vrna_array_init_size(dat->overlap3, fc->length + 1);

    for (unsigned int i = 0; i <= fc->length; i++)
      dat->overlap5[i] = NULL;

    for (unsigned int i = 0; i <= fc->length; i++)
      dat->overlap3[i] = NULL;

    /* ref-fill overlap entries */
    for (size_t cnt = 0; cnt < vrna_array_size(dat->regions); cnt++) {
      unsigned char failure = 0;

      r = &(dat->regions[cnt]);

      /* check whether the current interval pair supports a base pair already supported by a previous pair */
      for (unsigned int i = r->a; i <= r->b; i++) {
        if (dat->overlap5[i]) {
          for (size_t cnt = 0; cnt < vrna_array_size(dat->overlap5[i]); cnt++) {
            other = dat->overlap5[i][cnt];

            if (((r->c >= other->c) && (r->c <= other->d)) ||
                ((r->d >= other->c) && (r->d <= other->d))) {
              failure = 1;
              break;
            }
          }

          if (failure)
            break;
        }
      }

      if (failure) {
        vrna_log_error("Crossing overlaps detected! Omitting interval pair\n"
                       "\t[%u:%u] - [%u:%u], score=%d since it overlaps (in part) with\n"
                       "\t[%u:%u] - [%u:%u], score=%d",
                       r->a, r->b, r->c, r->d, r->score,
                       other->a, other->b, other->c, other->d, other->score);
      } else {
        for (unsigned int i = r->a; i <= r->b; i++) {
          if (dat->overlap5[i] == NULL)
            vrna_array_init(dat->overlap5[i])

          vrna_array_append(dat->overlap5[i], r);
        }

        for (unsigned int j = r->c; j <= r->d; j++) {
          if (dat->overlap3[j] == NULL)
            vrna_array_init(dat->overlap3[j])

          vrna_array_append(dat->overlap3[j], r);
        }

        if (fc->exp_params)
          r->exp_score = exp(-(FLT_OR_DBL)r->score * 10. / fc->exp_params->kT);
        else
          r->exp_score = 1.;
      }
    }

    dat->A  = (int *)vrna_alloc(sizeof(int) * (fc->length + 1));
    dat->B  = (int *)vrna_alloc(sizeof(int) * (fc->length + 1));

    for (unsigned int i = 0; i <= fc->length; i++)
      dat->A[i] = dat->B[i] = INF;

    dat->qA = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * (fc->length + 1));
    dat->qB = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * (fc->length + 1));

    for (unsigned int i = 0; i <= fc->length; i++)
      dat->qA[i] = dat->qB[i] = 0.;
  }

  return 0;
}



PUBLIC int
vrna_mfe_multibranch_loop_xlink(vrna_fold_compound_t        *fc,
                                unsigned int                i,
                                unsigned int                j,
                                struct vrna_mx_mfe_aux_ml_s *helpers)
{
  if ((fc) &&
      (helpers)) {
    unsigned char                 connects_intervals;
    int                           dangle_model, decomp, e;
    vrna_param_t                  *P;
    vrna_md_t                     *md;
    struct sc_mb_dat              sc_wrapper;
    struct vrna_crosslink_data_s  *dat;

    dat                 = fc->sc->iip;
    connects_intervals  = 0;
    P                   = fc->params;
    md                  = &(P->model_details);
    dangle_model        = md->dangles;
    e                   = decomp = INF;

    /* hack using old API to prevent calling vrna_mfe_multibranch_loop_xlink in a loop */
    decomp = vrna_E_mb_loop_fast(fc,
                                 i,
                                 j,
                                 helpers->DMLi1,
                                 helpers->DMLi2);

    /* first, check whether (i,j) connects two intervals of a segment */
    if (dat->overlap5[i]) {
      vrna_log_debug("pre_ml = %d", decomp);

      /* find supporting interval pair, if any */
      interact_region_s *r = NULL;
      for (size_t cnt = 0; cnt < vrna_array_size(dat->overlap5[i]); cnt++)
        if ((dat->overlap5[i][cnt]->c <= j) &&
            (dat->overlap5[i][cnt]->d >= j)) {
          r = dat->overlap5[i][cnt];
          break;
        }

      if (r) {
        init_sc_mb(fc, &sc_wrapper);

        switch (dangle_model) {
          case 0:
            decomp = ml_pair_d0(fc, i, j, r, helpers, &sc_wrapper);
            break;

          default:
            decomp = ml_pair_d2(fc, i, j, r, helpers, &sc_wrapper);
            break;
        }

        vrna_log_debug("mb_xl => %d [%u:%u] - [%u:%u]", decomp, r->a, r->b, r->c, r->d);
        free_sc_mb(&sc_wrapper);
      }
    }

    return decomp;
  }

  return INF;
}


PUBLIC unsigned int
vrna_bt_multibranch_loop_xlink(vrna_fold_compound_t *fc,
                               unsigned int         i,
                               unsigned int         j,
                               int                  en,
                               vrna_bps_t           bp_stack,
                               vrna_bts_t           bt_stack)
{
  if ((fc) &&
      (bp_stack) &&
      (bt_stack))
    return bt_mb_loop_xl(fc, i, j, en, bp_stack, bt_stack);

  return 0;
}


PRIVATE unsigned int
bt_mb_loop_xl(vrna_fold_compound_t  *fc,
              unsigned int          i,
              unsigned int          j,
              int                   en,
              vrna_bps_t bp_stack   VRNA_UNUSED,
              vrna_bts_t            bt_stack)
{
  unsigned char                 sliding_window;
  char                          *ptype, **ptype_local;
  short                         s5, s3, *S1, **SS, **S5, **S3;
  unsigned int                  *sn, n_seq, s, *tt, c1, c2, p, q, k, l, dangle_model,
                                type, type_2, turn, connects_intervals;
  int                           ij, e, tmp_en, *idx, *my_c, *my_fML, *my_fM2, *rtype, **c_local,
                                **fML_local;
  vrna_param_t                  *P;
  vrna_md_t                     *md;
  vrna_hc_t                     *hc;
  vrna_sc_t                     *sc;
  vrna_hc_eval_loop_f           evaluate;
  struct sc_mb_dat              sc_wrapper;
  struct vrna_crosslink_data_s  *dat;
  interact_region_s             *r;

  sc = (fc->type == VRNA_FC_TYPE_SINGLE) ? fc->sc : NULL;

  r = NULL;
  unsigned char                 use_default_bt = 1;

  if ((sc) &&
      (sc->iip) &&
      (sc->iip->overlap5[i])) {
    dat = sc->iip;
    
    /* find interval pair supported by (i, j), if any */
    for (size_t cnt = 0; cnt < vrna_array_size(dat->overlap5[i]); cnt++)
      if ((dat->overlap5[i][cnt]->c <= j) &&
          (dat->overlap5[i][cnt]->d >= j)) {
        r = dat->overlap5[i][cnt];
        use_default_bt = 0;
        break;
      }
  }

  vrna_log_debug("bt xlink (%u, %u) (default = %u)", i, j, use_default_bt);

  sliding_window  = (fc->hc->type == VRNA_HC_WINDOW) ? 1 : 0;
  n_seq           = (fc->type == VRNA_FC_TYPE_SINGLE) ? 1 : fc->n_seq;
  idx             = (sliding_window) ? NULL : fc->jindx;
  ij              = (sliding_window) ? 0 : idx[j] + i;
  S1              = (fc->type == VRNA_FC_TYPE_SINGLE) ? fc->sequence_encoding : NULL;
  SS              = (fc->type == VRNA_FC_TYPE_SINGLE) ? NULL : fc->S;
  S5              = (fc->type == VRNA_FC_TYPE_SINGLE) ? NULL : fc->S5;
  S3              = (fc->type == VRNA_FC_TYPE_SINGLE) ? NULL : fc->S3;
  P               = fc->params;
  md              = &(P->model_details);
  sn              = fc->strand_number;
  my_c            = (sliding_window) ? NULL : fc->matrices->c;
  my_fML          = (sliding_window) ? NULL : fc->matrices->fML;
  my_fM2          = (sliding_window) ? NULL : fc->matrices->fM2_real;
  c_local         = (sliding_window) ? fc->matrices->c_local : NULL;
  fML_local       = (sliding_window) ? fc->matrices->fML_local : NULL;
  ptype           = (fc->type ==
                     VRNA_FC_TYPE_SINGLE) ? (sliding_window ? NULL : fc->ptype) : NULL;
  ptype_local         = (sliding_window) ? fc->ptype_local : NULL;
  rtype               = &(md->rtype[0]);
  tt                  = NULL;
  type                = 0;
  dangle_model        = md->dangles;
  turn                = md->min_loop_size;
  hc                  = fc->hc;
  evaluate            = hc->eval_mb;
  connects_intervals  = 0;

  init_sc_mb(fc, &sc_wrapper);

  if (fc->type == VRNA_FC_TYPE_SINGLE)
    type = (sliding_window ?
            rtype[vrna_get_ptype_window(i, j, ptype_local)] :
            rtype[vrna_get_ptype(ij, ptype)]);

  p = i + 1;
  q = j - 1;

  k = q - 1;
  l = p + 1;

  s5  = -1;
  s3  = -1;

  c1 = c2 = VRNA_MX_FLAG_M;

  /* do two branches actually fit inside (i,j)? */
  if (i + 2 * (turn + 2) >= j)
    return 0;

  if (fc->type == VRNA_FC_TYPE_COMPARATIVE) {
    tt = (unsigned int *)vrna_alloc(sizeof(unsigned int) * n_seq);
    for (s = 0; s < n_seq; s++)
      tt[s] = vrna_get_ptype_md(SS[s][j], SS[s][i], md);
  } else {
    if (sn[q] == sn[j])
      s5 = S1[q];

    if (sn[i] == sn[p])
      s3 = S1[p];
  }

  if (evaluate(i, j, p, q, VRNA_DECOMP_PAIR_ML, hc)) {
    e = en -
        P->MLclosing *
        n_seq;

    if (dangle_model == 2) {
      switch (fc->type) {
        case VRNA_FC_TYPE_SINGLE:
          e -= vrna_E_multibranch_stem(type, s5, s3, P);
          break;

        case VRNA_FC_TYPE_COMPARATIVE:
          for (s = 0; s < n_seq; s++)
            e -= vrna_E_multibranch_stem(tt[s], S5[s][j], S3[s][i], P);
          break;
      }
    } else {
      switch (fc->type) {
        case VRNA_FC_TYPE_SINGLE:
          e -= vrna_E_multibranch_stem(type, -1, -1, P);
          break;

        case VRNA_FC_TYPE_COMPARATIVE:
          for (s = 0; s < n_seq; s++)
            e -= vrna_E_multibranch_stem(tt[s], -1, -1, P);
          break;
      }
    }

    if (sc_wrapper.pair)
      e -= sc_wrapper.pair(i, j, &sc_wrapper);

    /* hack to use default implementation if no crosslinking data involved */
    if ((r == NULL) ||
        (use_default_bt)) {
      for (k = p + 1; k < q; ++k) {
        if (evaluate(p, q, k - 1, k, VRNA_DECOMP_ML_ML_ML, hc)) {
          tmp_en = my_fML[idx[k - 1] + p] +
                   my_fML[idx[q] + k];

          if (sc_wrapper.decomp_ml)
            tmp_en += sc_wrapper.decomp_ml(p, q, k - 1, k, &sc_wrapper);

          if (e == tmp_en)
            goto odd_dangles_exit;
        }
      }
    } else {
      k = r->c;

      /* first 3 (obvious) cases: split point separates the intervals */
      if ((k > p + 1) &&
          (k < q) &&
          (evaluate(p, q, k - 1, k, VRNA_DECOMP_ML_ML_ML, hc))) {
        tmp_en = my_fML[idx[k - 1] + p] +
                 my_fML[idx[q] + k] +
                 r->score;

        if (sc_wrapper.decomp_ml)
          tmp_en += sc_wrapper.decomp_ml(p, q, k - 1, k, &sc_wrapper);

        if (e == tmp_en)
          goto odd_dangles_exit;
      }

      if ((k > p + 1) &&
          (evaluate(p, q, p, k - 1, VRNA_DECOMP_ML_ML, hc))) {
        tmp_en = my_fM2[idx[k - 1] + p] +
                 (q - k + 1) * P->MLbase +
                 r->score;

        if (sc_wrapper.red_ml)
          tmp_en += sc_wrapper.red_ml(p, q, p, k - 1, &sc_wrapper);

        if (e == tmp_en) {
          l = p;
          goto odd_dangles_exit;
        }
      }

      if ((k < q) &&
          (evaluate(p, q, k, q, VRNA_DECOMP_ML_ML, hc))) {
        tmp_en = my_fM2[idx[q] + k] +
                 (k - p) * P->MLbase +
                 r->score;

        if (sc_wrapper.red_ml)
          tmp_en += sc_wrapper.red_ml(p, q, k, q, &sc_wrapper);

        if (e == tmp_en) {
          l = q;
          goto odd_dangles_exit;
        }
      }

      /* nothing found so far, let's fill A and B matrices to allow for search of non-obious configurations */
      int *A = NULL, *B = NULL;

      {
        short         *S, *S2, si1, sj1;
        unsigned int  strands, *sn, tt;
        int           stem, decomp, score;

        A = (int *)vrna_alloc(sizeof(int) * (j - i));
        B = (int *)vrna_alloc(sizeof(int) * (j - i));

        A -= p; /* pointer magic */
        B -= p; /* pointer magic */

        score = r->score;

        for (l = p; l <= q; l++)
          /* reset A and B */
          A[l] = B[l] = INF;

        for (l = r->c; l <= q; l++) {
          /* got through all pairing partners k */
          for (unsigned int k = p; k < r->c; k++) {
            if ((hc->eval_mb(k, l, k, l, VRNA_DECOMP_ML_STEM, hc)) &&
                (my_c[idx[l] + k] != INF)) {
              strands = fc->strands;
              sn      = fc->strand_number;
              S       = fc->sequence_encoding;
              S2      = fc->sequence_encoding2;
              tt      = vrna_get_ptype_md(S2[k], S2[l], md);

              switch (dangle_model) {
                case 0:
                  si1 = sj1 = -1;
                  break;

                default:
                  si1 = ((strands == 1) || (sn[k - 1] == sn[k])) ? S[k - 1] : -1;
                  sj1 = ((strands == 1) || (sn[l] == sn[l + 1])) ? S[l + 1] : -1;
                  break;
              }

              stem = vrna_E_multibranch_stem(tt, si1, sj1, P) +
                     my_c[idx[l] + k];

              decomp = stem;

              if (k > r->b)
                decomp += score;

              /* fill A */
              if ((k > p) &&
                  (hc->eval_mb(p, l, k - 1, k, VRNA_DECOMP_ML_ML_STEM, hc))) {
                if (my_fML[idx[k - 1] + p] != INF) {
                  en = decomp +
                       my_fML[idx[k - 1] + p];

                  if (sc_wrapper.decomp_ml)
                    en += sc_wrapper.decomp_ml(p, l, k - 1, k, &sc_wrapper);

                  A[l] = MIN2(A[l], en);
                }
              }

              /* fill B */
              if (hc->eval_mb(p, l, k, l, VRNA_DECOMP_ML_ML, hc)) {
                en = decomp +
                     (k - p) * P->MLbase;

                if (sc_wrapper.red_ml)
                  en += sc_wrapper.red_ml(p, l, k, l, &sc_wrapper);

                B[l] = MIN2(B[l], en);
              }
            }
          }
        }
      }


      /* now, for the cases where we can't tell yet whether (i,j) is the innermost pair connecting I and I' */
      tmp_en = A[q];

      if (e == tmp_en) {
        free(tt);
        free_sc_mb(&sc_wrapper);
        return bt_A(fc, A, i, q, bp_stack, bt_stack);
      }

      for (k = r->c; k < q; k++) {
        if ((evaluate(p, q, p, k, VRNA_DECOMP_ML_ML, hc)) &&
            (A[k] != INF)) {
          tmp_en = A[k] +
                   (q - k) * P->MLbase;

          if (sc_wrapper.red_ml)
            tmp_en += sc_wrapper.red_ml(p, q, p, k, &sc_wrapper);

          if (e == tmp_en) {
            free(tt);
            free_sc_mb(&sc_wrapper);
            return bt_A(fc, A, i, k, bp_stack, bt_stack);
          }
        }

        if ((evaluate(p, q, k, k + 1, VRNA_DECOMP_ML_ML_ML, hc)) &&
            (my_fML[idx[q] + k + 1] != INF)) {
          tmp_en = my_fML[idx[q] + k + 1];

          if (sc_wrapper.decomp_ml)
            tmp_en += sc_wrapper.decomp_ml(p, q, k, k + 1, &sc_wrapper);

          if (e == A[k] + tmp_en) {
            free(tt);
            free_sc_mb(&sc_wrapper);
            vrna_bts_push(bt_stack,
                          ((vrna_sect_t){
              .i = k + 1,
              .j = q,
              .ml = c1
            }));
            return bt_A(fc, A, i, k, bp_stack, bt_stack);
          }

          if (e == B[k] + tmp_en) {
            free(tt);
            free_sc_mb(&sc_wrapper);
            vrna_bts_push(bt_stack,
                          ((vrna_sect_t){
              .i = k + 1,
              .j = q,
              .ml = c1
            }));
            return bt_B(fc, B, i, k, bp_stack, bt_stack);
          }
        }
      }
    }
  }

odd_dangles_exit:

  free(tt);
  free_sc_mb(&sc_wrapper);

  /* find split point for M2 entry */
  if (l == p) {
    return bt_m2(fc, p, k - 1, my_fM2[idx[k - 1] + p], bp_stack, bt_stack);
  } else if (l == q) {
    return bt_m2(fc, k, q, my_fM2[idx[q] + k], bp_stack, bt_stack);
  } else if (k + 2 < j) {
    /* already found split point */
    vrna_bts_push(bt_stack,
                  ((vrna_sect_t){
      .i = p,
      .j = k - 1,
      .ml = c1
    }));
    vrna_bts_push(bt_stack,
                  ((vrna_sect_t){
      .i = k,
      .j = q,
      .ml = c2
    }));
    return 1;
  }

  vrna_log_error("bt M failed k = %d, l = %d, i  %d, j = %d", k, l, i, j);

  return 0;
}


/* find split point u for base pair (u,v) */
PRIVATE unsigned int
bt_A(vrna_fold_compound_t *fc,
     int                  *A,
     unsigned int         i,
     unsigned int         v,
     vrna_bps_t bp_stack  VRNA_UNUSED,
     vrna_bts_t           bt_stack)
{
  short                         *S, *S2, si1, sj1;
  unsigned int                  u, *sn, strands, tt;
  int                           *my_c, *my_fML, *idx, stem, decomp, dangle_model, score;
  vrna_param_t                  *P;
  vrna_md_t                     *md;
  vrna_hc_t                     *hc;
  vrna_sc_t                     *sc;
  vrna_hc_eval_loop_f           evaluate;
  struct sc_mb_dat              sc_wrapper;
  struct vrna_crosslink_data_s  *dat;

  P             = fc->params;
  md            = &(P->model_details);
  dangle_model  = md->dangles;
  idx           = fc->jindx;
  my_c          = fc->matrices->c;
  my_fML        = fc->matrices->fML;
  hc            = fc->hc;
  evaluate      = hc->eval_mb;
  sc            = fc->sc;
  dat           = fc->sc->iip;

  init_sc_mb(fc, &sc_wrapper);

  /* find interval pair (I,I') where i \in I and v \in I' */
  interact_region_s *r = NULL;
  
  for (size_t cnt = 0; cnt < vrna_array_size(dat->overlap5[i]); cnt++)
    if ((dat->overlap5[i][cnt]->c <= v) &&
        (dat->overlap5[i][cnt]->d > v)) {
      r = dat->overlap5[i][cnt];
      break;
    }

  if (r) {
    score = r->score;

    for (u = i + 2; u < r->c; u++) {
      if ((evaluate(i + 1, v, u - 1, u, VRNA_DECOMP_ML_ML_STEM, hc)) &&
          (my_c[idx[v] + u] != INF)) {
        strands = fc->strands;
        sn      = fc->strand_number;
        S       = fc->sequence_encoding;
        S2      = fc->sequence_encoding2;
        tt      = vrna_get_ptype_md(S2[u], S2[v], md);

        switch (dangle_model) {
          case 0:
            si1 = sj1 = -1;
            break;

          default:
            si1 = ((strands == 1) || (sn[u - 1] == sn[u])) ? S[u - 1] : -1;
            sj1 = ((strands == 1) || (sn[v] == sn[v + 1])) ? S[v + 1] : -1;
            break;
        }

        stem = vrna_E_multibranch_stem(tt, si1, sj1, P) +
               my_c[idx[v] + u] +
               my_fML[idx[u - 1] + i + 1];

        decomp = 0;

        /* add score if pairing partner u doesn't connect both intervals */
        if (u > r->b)
          decomp += score;

        decomp += stem;

        if (sc_wrapper.decomp_ml)
          decomp += sc_wrapper.decomp_ml(i + 1, v, u - 1, u, &sc_wrapper);

        if (A[v] == decomp)
          goto exit_bt_A;
      }
    }
  }

exit_bt_A:

  A += i + 1;
  free(A);

  free_sc_mb(&sc_wrapper);

  if ((r) &&
      (u < v)) {
    vrna_bts_push(bt_stack,
                  ((vrna_sect_t){
      .i = i + 1,
      .j = u - 1,
      .ml = VRNA_MX_FLAG_M
    }));
    vrna_bts_push(bt_stack,
                  ((vrna_sect_t){
      .i = u,
      .j = v,
      .ml = VRNA_MX_FLAG_C
    }));

    return 1;
  }

  vrna_log_error("bt A failed");

  return 0;
}


PRIVATE unsigned int
bt_B(vrna_fold_compound_t *fc,
     int                  *B,
     unsigned int         i,
     unsigned int         v,
     vrna_bps_t bp_stack  VRNA_UNUSED,
     vrna_bts_t           bt_stack)
{
  short                         *S, *S2, si1, sj1;
  unsigned int                  u, *sn, strands, tt;
  int                           *my_c, *idx, stem, decomp, dangle_model, score;
  vrna_param_t                  *P;
  vrna_md_t                     *md;
  vrna_hc_t                     *hc;
  vrna_sc_t                     *sc;
  vrna_hc_eval_loop_f           evaluate;
  struct sc_mb_dat              sc_wrapper;
  struct vrna_crosslink_data_s  *dat;

  P             = fc->params;
  md            = &(P->model_details);
  dangle_model  = md->dangles;
  idx           = fc->jindx;
  my_c          = fc->matrices->c;
  hc            = fc->hc;
  evaluate      = hc->eval_mb;
  sc            = fc->sc;
  dat           = fc->sc->iip;

  init_sc_mb(fc, &sc_wrapper);

  /* find interval pair (I,I') where i \in I and v \in I' */
  interact_region_s *r = NULL;
  
  for (size_t cnt = 0; cnt < vrna_array_size(dat->overlap5[i]); cnt++)
    if ((dat->overlap5[i][cnt]->c <= v) &&
        (dat->overlap5[i][cnt]->d > v)) {
      r = dat->overlap5[i][cnt];
      break;
    }

  if (r) {
    score = r->score;

    for (u = i + 1; u < r->c; u++) {
      if ((evaluate(i + 1, v, u, v, VRNA_DECOMP_ML_ML, hc)) &&
          (my_c[idx[v] + u] != INF)) {
        strands = fc->strands;
        sn      = fc->strand_number;
        S       = fc->sequence_encoding;
        S2      = fc->sequence_encoding2;
        tt      = vrna_get_ptype_md(S2[u], S2[v], md);

        switch (dangle_model) {
          case 0:
            si1 = sj1 = -1;
            break;

          default:
            si1 = ((strands == 1) || (sn[u - 1] == sn[u])) ? S[u - 1] : -1;
            sj1 = ((strands == 1) || (sn[v] == sn[v + 1])) ? S[v + 1] : -1;
            break;
        }

        stem = vrna_E_multibranch_stem(tt, si1, sj1, P) +
               my_c[idx[v] + u] +
               (u - i - 1) * P->MLbase;

        decomp = 0;

        /* add score if pairing partner u doesn't connect both intervals */
        if (u > r->b)
          decomp += score;

        decomp += stem;

        if (sc_wrapper.red_ml)
          decomp += sc_wrapper.red_ml(i + 1, v, u, v, &sc_wrapper);

        if (B[v] == decomp)
          goto exit_bt_B;
      }
    }
  }

exit_bt_B:

  B += i + 1;
  free(B);

  free_sc_mb(&sc_wrapper);

  if ((r) &&
      (u < v)) {
    vrna_bts_push(bt_stack,
                  ((vrna_sect_t){
      .i = u,
      .j = v,
      .ml = VRNA_MX_FLAG_C
    }));

    return 1;
  }

  vrna_log_error("bt B failed");
  return 0;
}


PRIVATE unsigned int
bt_m2(vrna_fold_compound_t  *fc,
      unsigned int          i,
      unsigned int          j,
      int                   en,
      vrna_bps_t bp_stack   VRNA_UNUSED,
      vrna_bts_t            bt_stack)
{
  unsigned int        k;
  int                 *idx, *my_fML, tmp_en;
  vrna_hc_t           *hc;
  vrna_sc_t           *sc;
  vrna_hc_eval_loop_f evaluate;
  struct sc_mb_dat    sc_wrapper;

  idx       = fc->jindx;
  my_fML    = fc->matrices->fML;
  hc        = fc->hc;
  evaluate  = hc->eval_mb;

  init_sc_mb(fc, &sc_wrapper);

  for (k = i + 1; k < j - 1; ++k) {
    if (evaluate(i, j, k, k + 1, VRNA_DECOMP_ML_ML_ML, hc)) {
      tmp_en = my_fML[idx[k] + i] +
               my_fML[idx[j] + k + 1];

      if (sc_wrapper.decomp_ml)
        tmp_en += sc_wrapper.decomp_ml(i, j, k, k + 1, &sc_wrapper);

      if (en == tmp_en)
        goto odd_dangles_exit2;
    }
  }

odd_dangles_exit2:

  free_sc_mb(&sc_wrapper);

  if (k + 1 < j) {
    /* find split point for M2 entry */
    vrna_bts_push(bt_stack,
                  ((vrna_sect_t){
      .i = i,
      .j = k,
      .ml = VRNA_MX_FLAG_M
    }));
    vrna_bts_push(bt_stack,
                  ((vrna_sect_t){
      .i = k + 1,
      .j = j,
      .ml = VRNA_MX_FLAG_M
    }));
    return 1;
  }

  vrna_log_error("bt fm2 failed, segment [%u, %u]", i, j);
  return 0;
}


/*
 *  at this point we know that there is an interval pair (I,I') with
 *  I.a <= i <= I.b and I'.c <= j <= I'.d
 */
PRIVATE INLINE int
ml_pair_d2(vrna_fold_compound_t         *fc,
           unsigned int                 i,
           unsigned int                 j,
           interact_region_s            *r,
           struct vrna_mx_mfe_aux_ml_s  *helpers,
           struct sc_mb_dat             *sc_wrapper)
{
  short         *S, *S2, **SS, **S5, **S3, si1, sj1;
  unsigned int  tt, strands, *sn, p, q, s, n_seq;
  int           e;
  vrna_param_t  *P;
  vrna_md_t     *md;
  vrna_hc_t     *hc;

  e   = INF;
  hc  = fc->hc;
  p   = i + 1;
  q   = j - 1;

  if (hc->eval_mb(i, j, p, q, VRNA_DECOMP_PAIR_ML, hc)) {
    /* get contribution for enclosed part */
    e = mfe_M2_xlink(fc, i, j, r, helpers, sc_wrapper);

    if (e != INF) {
      P   = fc->params;
      md  = &(P->model_details);

      switch (fc->type) {
        case VRNA_FC_TYPE_SINGLE:
          strands = fc->strands;
          sn      = fc->strand_number;
          S       = fc->sequence_encoding;
          S2      = fc->sequence_encoding2;
          tt      = vrna_get_ptype_md(S2[j], S2[i], md);

          if (md->noGUclosure && ((tt == 3) || (tt == 4)))
            return INF; /* not allowed */

          si1 = ((strands == 1) || (sn[i] == sn[p])) ? S[p] : -1;
          sj1 = ((strands == 1) || (sn[q] == sn[j])) ? S[q] : -1;

          e += vrna_E_multibranch_stem(tt, sj1, si1, P) +
               P->MLclosing;
          break;

        case VRNA_FC_TYPE_COMPARATIVE:
          n_seq = fc->n_seq;
          SS    = fc->S;
          S5    = fc->S5;
          S3    = fc->S3;

          for (s = 0; s < n_seq; s++) {
            tt  = vrna_get_ptype_md(SS[s][j], SS[s][i], md);
            e   += vrna_E_multibranch_stem(tt, S5[s][j], S3[s][i], P);
          }

          e += n_seq * P->MLclosing;
          break;
      }

      if (sc_wrapper->pair)
        e += sc_wrapper->pair(i, j, sc_wrapper);
    }
  }

  return e;
}


PRIVATE INLINE int
ml_pair_d0(vrna_fold_compound_t         *fc,
           unsigned int                 i,
           unsigned int                 j,
           interact_region_s            *r,
           struct vrna_mx_mfe_aux_ml_s  *helpers,
           struct sc_mb_dat             *sc_wrapper)
{
  short         *S, *S2, **SS, **S5, **S3, si1, sj1;
  unsigned int  tt, strands, *sn, p, q, s, n_seq;
  int           e;
  vrna_param_t  *P;
  vrna_md_t     *md;
  vrna_hc_t     *hc;

  e   = INF;
  hc  = fc->hc;
  p   = i + 1;
  q   = j - 1;

  if (hc->eval_mb(i, j, p, q, VRNA_DECOMP_PAIR_ML, hc)) {
    /* get contribution for enclosed part */
    e = mfe_M2_xlink(fc, i, j, r, helpers, sc_wrapper);

    if (e != INF) {
      P   = fc->params;
      md  = &(P->model_details);

      switch (fc->type) {
        case VRNA_FC_TYPE_SINGLE:
          strands = fc->strands;
          sn      = fc->strand_number;
          S       = fc->sequence_encoding;
          S2      = fc->sequence_encoding2;
          tt      = vrna_get_ptype_md(S2[j], S2[i], md);

          if (md->noGUclosure && ((tt == 3) || (tt == 4)))
            return INF; /* not allowed */

          e += vrna_E_multibranch_stem(tt, -1, -1, P) +
               P->MLclosing;
          vrna_log_debug("e_closing = %d", vrna_E_multibranch_stem(tt, -1, -1, P) + P->MLclosing);
          break;

        case VRNA_FC_TYPE_COMPARATIVE:
          n_seq = fc->n_seq;
          SS    = fc->S;
          S5    = fc->S5;
          S3    = fc->S3;

          for (s = 0; s < n_seq; s++) {
            tt  = vrna_get_ptype_md(SS[s][j], SS[s][i], md);
            e   += vrna_E_multibranch_stem(tt, -1, -1, P);
          }

          e += n_seq * P->MLclosing;
          break;
      }

      if (sc_wrapper->pair)
        e += sc_wrapper->pair(i, j, sc_wrapper);
    }
  }

  return e;
}


/*  compute the MFE contribution for the enclosed part of a multibranch loop
 *  [i + 1, j - 1] with closing pair (i,j), given that (i,j) connects an
 *  interval pair
 */
PRIVATE INLINE int
mfe_M2_xlink(vrna_fold_compound_t         *fc,
             unsigned int                 i,
             unsigned int                 j,
             interact_region_s            *r,
             struct vrna_mx_mfe_aux_ml_s  *helpers,
             struct sc_mb_dat             *sc_wrapper)
{
  unsigned int                  v, k, p, q;
  int                           e, decomp, *fML, *fM2, *A, *B, *idx;
  vrna_param_t                  *P;
  vrna_hc_t                     *hc;
  struct vrna_crosslink_data_s  *dat;

  dat = fc->sc->iip;
  e   = INF;
  P   = fc->params;
  hc  = fc->hc;
  idx = fc->jindx;
  fML = fc->matrices->fML;
  fM2 = fc->matrices->fM2_real;
  A   = dat->A;
  B   = dat->B;
  p   = i + 1;
  q   = j - 1;

  if (r) {
    k = r->c;

    /* cases where with no base pair (u, v) with u < c and c <= v < d */
    if ((k > p + 1) &&
        (k < q) &&
        (hc->eval_mb(p, q, k - 1, k, VRNA_DECOMP_ML_ML_ML, hc))) {
      if ((fML[idx[k - 1] + p] != INF) &&
          (fML[idx[q] + k] != INF)) {
        decomp = fML[idx[k - 1] + p] +
                 fML[idx[q] + k] +
                 r->score;

        if (sc_wrapper->decomp_ml)
          decomp += sc_wrapper->decomp_ml(p, q, k - 1, k, sc_wrapper);

        e = MIN2(e, decomp);
        vrna_log_debug("e1 = %d", decomp);
      }
    }

    if ((k > p + 1) &&
        (hc->eval_mb(p, q, p, k - 1, VRNA_DECOMP_ML_ML, hc))) {
      if (fM2[idx[k - 1] + p] != INF) {
        decomp = fM2[idx[k - 1] + p] +
                 (q - k + 1) * P->MLbase +
                 r->score;

        if (sc_wrapper->red_ml)
          decomp += sc_wrapper->red_ml(p, q, p, k - 1, sc_wrapper);

        e = MIN2(e, decomp);
        vrna_log_debug("e2 = %d", decomp);
      }
    }

    if ((k < q) &&
        (hc->eval_mb(p, q, k, q, VRNA_DECOMP_ML_ML, hc))) {
      if (fM2[idx[q] + k] != INF) {
        decomp = fM2[idx[q] + k] +
                 (k - p) * P->MLbase +
                 r->score;

        if (sc_wrapper->red_ml)
          decomp += sc_wrapper->red_ml(p, q, k, q, sc_wrapper);

        e = MIN2(e, decomp);
        vrna_log_debug("e3 = %d", decomp);
      }
    }

    /* now, all cases where we allow for base pairs that cross the second interval boundary */
    decomp = A[q];
    vrna_log_debug("e4 = %d", decomp);
    e = MIN2(e, decomp);

    for (v = r->c; v < q; v++) {
      vrna_log_debug("v=%u, A[v] = %d, B[v] = %d", v, A[v], B[v]);
      if (hc->eval_mb(p, q, p, v, VRNA_DECOMP_ML_ML, hc)) {
        decomp = A[v];
        if (decomp != INF) {
          decomp += (q - v) * P->MLbase;

          if (sc_wrapper->red_ml)
            decomp += sc_wrapper->red_ml(p, q, p, v, sc_wrapper);

          e = MIN2(e, decomp);
          vrna_log_debug("e5 = %d", decomp);
        }
      }

      if (hc->eval_mb(p, q, v, v + 1, VRNA_DECOMP_ML_ML_ML, hc)) {
        if (fML[idx[q] + v + 1] != INF) {
          decomp = A[v];

          if (decomp != INF) {
            decomp += fML[idx[q] + v + 1];

            if (sc_wrapper->decomp_ml)
              decomp += sc_wrapper->decomp_ml(p, q, v, v + 1, sc_wrapper);

            e = MIN2(e, decomp);
            vrna_log_debug("e6 = %d", decomp);
          }

          decomp = B[v];

          if (decomp != INF) {
            decomp += fML[idx[q] + v + 1];

            if (sc_wrapper->decomp_ml)
              decomp += sc_wrapper->decomp_ml(p, q, v, v + 1, sc_wrapper);

            e = MIN2(e, decomp);
            vrna_log_debug("e7 = %d", decomp);
          }
        }
      }
    }
  }

  return e;
}


PUBLIC void
vrna_sc_iip_update_mfe(vrna_fold_compound_t *fc,
                       unsigned int         i)
{
  short                         *S, *S2, si1, sj1;
  unsigned int                  n, *sn, strands, tt, turn;
  int                           *A, *B, *fML, *c, *idx, dangle_model, score, stem, decomp, en;
  vrna_param_t                  *P;
  vrna_md_t                     *md;
  vrna_hc_t                     *hc;
  struct vrna_crosslink_data_s  *dat;
  struct sc_mb_dat              sc_wrapper;
  interact_region_s             *r;

  dat           = fc->sc->iip;
  n             = fc->length;
  P             = fc->params;
  md            = &(P->model_details);
  dangle_model  = md->dangles;
  turn          = md->min_loop_size;
  idx           = fc->jindx;
  c             = fc->matrices->c;
  fML           = fc->matrices->fML;
  hc            = fc->hc;
  A             = dat->A;
  B             = dat->B;
  r             = NULL;

  init_sc_mb(fc, &sc_wrapper);

  /* i overlaps with I, add all ML parts where a base pair starts in I and pairs outside of I */

  for (unsigned int v = i; v <= n; v++)
    /* reset A and B */
    A[v] = B[v] = INF;


  for (unsigned int v = i + turn + 1; v <= n; v++) {
    r = NULL;
    /* check whether there is an interval pair (I,I') with i \in I and v \in I' */
    if (dat->overlap5[i])
      for (size_t cnt = 0; cnt < vrna_array_size(dat->overlap5[i]); cnt++)
        if ((dat->overlap5[i][cnt]->c <= v) &&
            (dat->overlap5[i][cnt]->d > v)) {
          r = dat->overlap5[i][cnt];
          break;
        }

    if (r) {
      score = r->score;
      vrna_log_debug("score = %d I=([%u:%u], [%u:%u])", score, r->a, r->b, r->c, r->d);

      /*
       *  collect energies of ML parts where a base pair (u, v) crosses I', i.e.
       *  u < c && c <= v < d
       */
      for (unsigned int u = i + 1; u < r->c; u++) {
        if ((hc->eval_mb(u, v, u, v, VRNA_DECOMP_ML_STEM, hc)) &&
            (c[idx[v] + u] != INF)) {
          vrna_log_debug("c[%u,%u] = %d", u, v, c[idx[v] + u]);
          strands = fc->strands;
          sn      = fc->strand_number;
          S       = fc->sequence_encoding;
          S2      = fc->sequence_encoding2;
          tt      = vrna_get_ptype_md(S2[u], S2[v], md);

          switch (dangle_model) {
            case 0:
              si1 = sj1 = -1;
              break;

            default:
              si1 = ((strands == 1) || (sn[u - 1] == sn[u])) ? S[u - 1] : -1;
              sj1 = ((strands == 1) || (sn[v] == sn[v + 1])) ? S[v + 1] : -1;
          }

          stem = vrna_E_multibranch_stem(tt, si1, sj1, P) +
                 c[idx[v] + u];

          vrna_log_debug("stem(%u,%u) = %d", u, v, vrna_E_multibranch_stem(tt, si1, sj1, P));
          decomp = stem;

          if (u > r->b)
            decomp += score;

          /* fill A */
          if ((i + 1 < u) &&
              (hc->eval_mb(i + 1, v, u - 1, u, VRNA_DECOMP_ML_ML_STEM, hc))) {
            if (fML[idx[u - 1] + i + 1] != INF) {
              en = decomp +
                   fML[idx[u - 1] + i + 1];

              if (sc_wrapper.decomp_ml)
                en += sc_wrapper.decomp_ml(i + 1, v, u - 1, u, &sc_wrapper);

              vrna_log_debug("A[%u] = MIN2(%d, %d)", v, A[v], en);
              A[v] = MIN2(A[v], en);
            }
          }

          /* fill B */
          if (hc->eval_mb(i + 1, v, u, v, VRNA_DECOMP_ML_ML, hc)) {
            en = decomp +
                 (u - 1 - i) * P->MLbase;

            if (sc_wrapper.red_ml)
              en += sc_wrapper.red_ml(i + 1, v, u, v, &sc_wrapper);

            vrna_log_debug("B[%u] = MIN2(%d, %d)", v, B[v], en);
            B[v] = MIN2(B[v], en);
          }
        }
      }
    }
  }

  free_sc_mb(&sc_wrapper);
}


PUBLIC FLT_OR_DBL
vrna_pf_multibranch_loop_xlink(vrna_fold_compound_t       *fc,
                               unsigned int               i,
                               unsigned int               j,
                               struct vrna_mx_pf_aux_ml_s *helpers)
{
  if ((fc) &&
      (helpers)) {
    FLT_OR_DBL                    q = 0.;
    struct vrna_crosslink_data_s  *dat;
    interact_region_s             *r;

    dat = fc->sc->iip;

    /* first, check whether (i,j) connects two intervals of a segment */
    if (dat->overlap3[j]) {
      r = NULL;

      for (size_t cnt = 0; cnt < vrna_array_size(dat->overlap3[j]); cnt++)
        if ((dat->overlap3[j][cnt]->a <= i) &&
            (dat->overlap3[j][cnt]->b >= i)) {
          r = dat->overlap3[j][cnt];
          break;
        }

      if (r) {
        q = ml_pair_exp(fc, i, j, r, helpers);
        vrna_log_debug("mb_xl => %g [%u:%u] - [%u:%u]", q, r->a, r->b, r->c, r->d);
      } else {
        /* get default contribution */
        q = vrna_exp_E_mb_loop_fast(fc,
                                    i,
                                    j,
                                    helpers);
      }
    } else {
      /* get default contribution */
      q = vrna_exp_E_mb_loop_fast(fc,
                                  i,
                                  j,
                                  helpers);
    }

    return q;
  }

  return 0.;
}


/*
 *  at this point we know that there is an interval pair (I,I') with
 *  I.a <= i <= I.b and I'.c <= j <= I'.d
 */
PRIVATE INLINE FLT_OR_DBL
ml_pair_exp(vrna_fold_compound_t        *fc,
            unsigned int                i,
            unsigned int                j,
            interact_region_s           *r,
            struct vrna_mx_pf_aux_ml_s  *helpers)
{
  unsigned char         sliding_window;
  char                  *ptype, **ptype_local;
  short                 *S1, **SS, **S5, **S3;
  unsigned int          *sn, n_seq, s, *se, tt;
  int                   ij, *my_iindx, *jindx, *rtype;
  FLT_OR_DBL            qbt1, temp, qqqmmm, *qm, *qm2, **qm_local, *scale, expMLclosing, *qqm1,
                        *qqm21;
  vrna_hc_t             *hc;
  vrna_exp_param_t      *pf_params;
  vrna_md_t             *md;
  struct sc_mb_exp_dat  sc_wrapper;

  sliding_window  = (fc->hc->type == VRNA_HC_WINDOW) ? 1 : 0;
  n_seq           = (fc->type == VRNA_FC_TYPE_SINGLE) ? 1 : fc->n_seq;
  se              = fc->strand_end;
  my_iindx        = (sliding_window) ? NULL : fc->iindx;
  jindx           = (sliding_window) ? NULL : fc->jindx;
  ptype           = (fc->type == VRNA_FC_TYPE_SINGLE) ? (sliding_window ? NULL : fc->ptype) : NULL;
  ptype_local     = (sliding_window) ? fc->ptype_local : NULL;
  S1              = (fc->type == VRNA_FC_TYPE_SINGLE) ? fc->sequence_encoding : NULL;
  SS              = (fc->type == VRNA_FC_TYPE_SINGLE) ? NULL : fc->S;
  S5              = (fc->type == VRNA_FC_TYPE_SINGLE) ? NULL : fc->S5;
  S3              = (fc->type == VRNA_FC_TYPE_SINGLE) ? NULL : fc->S3;
  qm              = (sliding_window) ? NULL : fc->exp_matrices->qm;
  qm2             = (sliding_window) ? NULL : fc->exp_matrices->qm2_real;
  qm_local        = (sliding_window) ? fc->exp_matrices->qm_local : NULL;
  scale           = fc->exp_matrices->scale;
  pf_params       = fc->exp_params;
  md              = &(pf_params->model_details);
  ij              = (sliding_window) ? 0 : jindx[j] + i;
  sn              = fc->strand_number;
  hc              = fc->hc;
  expMLclosing    = pf_params->expMLclosing;
  qbt1            = 0.;
  rtype           = &(md->rtype[0]);

  if (hc->eval_mb(i, j, i + 1, j - 1, VRNA_DECOMP_PAIR_ML, hc)) {
    init_sc_mb_exp(fc, &sc_wrapper);

    qbt1 = pow(expMLclosing, (double)n_seq) *
           scale[2];

    switch (fc->type) {
      case VRNA_FC_TYPE_SINGLE:
        tt = (sliding_window) ?
             rtype[vrna_get_ptype_window(i, j + i, ptype_local)] :
             rtype[vrna_get_ptype(ij, ptype)];
        qbt1 *= vrna_exp_E_multibranch_stem(tt, S1[j - 1], S1[i + 1], pf_params);

        break;


      case VRNA_FC_TYPE_COMPARATIVE:
        for (s = 0; s < n_seq; s++) {
          tt    = vrna_get_ptype_md(SS[s][j], SS[s][i], md);
          qbt1  *= vrna_exp_E_multibranch_stem(tt, S5[s][j], S3[s][i], pf_params);
        }
        break;
    }

    if (sc_wrapper.pair)
      qbt1 *= sc_wrapper.pair(i, j, &sc_wrapper);

    qqqmmm = qbt1;

    temp = pf_M2_xlink(fc, i, j, r, helpers, &sc_wrapper);
    /* get contribution for enclosed part */
    qbt1 *= temp;

    vrna_log_debug("qb^M[%u,%u] = %g * %g = %g", i, j, temp, qqqmmm, temp * qqqmmm);

    free_sc_mb_exp(&sc_wrapper);
  }

  return qbt1;
}


PRIVATE INLINE FLT_OR_DBL
pf_M2_xlink(vrna_fold_compound_t        *fc,
            unsigned int                i,
            unsigned int                j,
            interact_region_s           *r,
            struct vrna_mx_pf_aux_ml_s  *helpers,
            struct sc_mb_exp_dat        *sc_wrapper)
{
  unsigned int                  v, u, k, l;
  int                           *idx;
  FLT_OR_DBL                    q, qbt1, *qm, *qm2, *expMLbase, *qA, *qB;
  vrna_hc_t                     *hc;
  struct vrna_crosslink_data_s  *dat;

  dat       = fc->sc->iip;
  q         = 0.;
  hc        = fc->hc;
  idx       = fc->iindx;
  qm        = fc->exp_matrices->qm;
  qm2       = fc->exp_matrices->qm2_real;
  expMLbase = fc->exp_matrices->expMLbase;

  qA  = dat->qA;
  qB  = dat->qB;
  k   = i + 1;
  l   = j - 1;

  u = r->b;

  /* cases where with no base pair (u, v) with u < c and c <= v < d */
  if ((i == r->b) ||
      (r->c == j)) {
    qbt1 = qm2[idx[k] - l] *
           r->exp_score;

    q += qbt1;
    vrna_log_debug("q0 = %g", qbt1);
  } else {
    if ((u > k) &&
        (u < l) &&
        (hc->eval_mb(k, l, u, u + 1, VRNA_DECOMP_ML_ML_ML, hc))) {
      qbt1 = qm[idx[k] - u] *
             qm[idx[u + 1] - l] *
             r->exp_score;

      if (sc_wrapper->decomp_ml)
        qbt1 *= sc_wrapper->decomp_ml(k, l, u, u + 1, sc_wrapper);

      q += qbt1;
      vrna_log_debug("q1 = %g", qbt1);
    }

    if ((u > k) &&
        (hc->eval_mb(k, l, k, u, VRNA_DECOMP_ML_ML, hc))) {
      qbt1 = qm2[idx[k] - u] *
             expMLbase[l - u] *
             r->exp_score;

      if (sc_wrapper->red_ml)
        qbt1 *= sc_wrapper->red_ml(k, l, k, u, sc_wrapper);

      q += qbt1;
      vrna_log_debug("q2 = %g", qbt1);
    }

    if ((u + 1 < l) &&
        (hc->eval_mb(k, l, u + 1, l, VRNA_DECOMP_ML_ML, hc))) {
      qbt1 = qm2[idx[u + 1] - l] *
             expMLbase[u + 1 - k] *
             r->exp_score;

      if (sc_wrapper->red_ml)
        qbt1 *= sc_wrapper->red_ml(k, l, u + 1, l, sc_wrapper);

      q += qbt1;
      vrna_log_debug("q3 = %g", qbt1);
    }

  /* now, all cases where we allow for base pairs that cross the first interval boundary */
  qbt1 = qA[k];
  vrna_log_debug("q4 = %g", qbt1);

  q += qbt1;

  for (u = k + 1; u <= r->b; u++) {
    vrna_log_debug("u=%u, qA[u] = %g, qB[u] = %g", u, qA[u], qB[u]);
    if (hc->eval_mb(k, l, u, l, VRNA_DECOMP_ML_ML, hc)) {
      qbt1 = qA[u] *
             expMLbase[u - k];

      if (sc_wrapper->red_ml)
        qbt1 *= sc_wrapper->red_ml(k, l, u, l, sc_wrapper);

      q += qbt1;
      vrna_log_debug("q5 = %g", qbt1);
    }

    if ((k + 1 < u) &&
        (hc->eval_mb(k, l, u - 1, u, VRNA_DECOMP_ML_ML_ML, hc))) {
      qbt1 = qA[u] *
             qm[idx[k] - (u - 1)];

      if (sc_wrapper->decomp_ml)
        qbt1 *= sc_wrapper->decomp_ml(k, l, u - 1, u, sc_wrapper);

      q += qbt1;
      vrna_log_debug("q6 = %g", qbt1);

      qbt1 = qB[u] *
             qm[idx[k] - u + 1];

      if (sc_wrapper->decomp_ml)
        qbt1 *= sc_wrapper->decomp_ml(k, l, u - 1, u, sc_wrapper);

      q += qbt1;
      vrna_log_debug("q7 = %g", qbt1);
    }
  }
  }

  return q;
}


PUBLIC void
vrna_sc_iip_update_pf(vrna_fold_compound_t *fc,
                     unsigned int         j)
{
  short                         *S, *S2, si1, sj1;
  unsigned int                  n, *sn, strands, tt, turn;
  FLT_OR_DBL                    *qA, *qB, *qm, *qb, *expMLbase, score, stem, qbt1, qq;
  int                           *idx;
  vrna_exp_param_t              *pf_params;
  vrna_md_t                     *md;
  vrna_hc_t                     *hc;
  struct vrna_crosslink_data_s  *dat;
  struct sc_mb_exp_dat          sc_wrapper;
  interact_region_s             *r;

  dat       = fc->sc->iip;
  n         = fc->length;
  pf_params = fc->exp_params;
  md        = &(fc->params->model_details);
  turn      = md->min_loop_size;
  idx       = fc->iindx;
  qb        = fc->exp_matrices->qb;
  qm        = fc->exp_matrices->qm;
  expMLbase = fc->exp_matrices->expMLbase;
  hc        = fc->hc;
  qA        = dat->qA;
  qB        = dat->qB;

  init_sc_mb_exp(fc, &sc_wrapper);

  /* i overlaps with I, add all ML parts where a base pair starts in I and pairs outside of I */

  for (unsigned int v = 1; v <= j; v++)
    /* reset A and B */
    qA[v] = qB[v] = 0.;

  for (unsigned int u = 1; u < j; u++) {
    r = NULL;

    if (dat->overlap3[j])
      for (size_t cnt = 0; cnt < vrna_array_size(dat->overlap3[j]); cnt++)
        if ((dat->overlap3[j][cnt]->a < u) &&
            (dat->overlap3[j][cnt]->b >= u)) {
          r = dat->overlap3[j][cnt];
          break;
        }

    if (r) {
      score = r->exp_score;
      vrna_log_debug("exp_score = %g I=([%u:%u], [%u:%u])", score, r->a, r->b, r->c, r->d);

      /*
       *  collect energies of ML parts where a base pair (u, v) crosses I, i.e.
       *  a < u <= b && b < v < d
       */
      for (unsigned int v = r->b + 1; v < j; v++) {
        if (hc->eval_mb(u, v, u, v, VRNA_DECOMP_ML_STEM, hc)) {
          vrna_log_debug("qb[%u,%u] = %g", u, v, qb[idx[u] - v]);
          strands = fc->strands;
          sn      = fc->strand_number;
          S       = fc->sequence_encoding;
          S2      = fc->sequence_encoding2;
          tt      = vrna_get_ptype_md(S2[u], S2[v], md);

          si1 = ((strands == 1) || (sn[u - 1] == sn[u])) ? S[u - 1] : -1;
          sj1 = ((strands == 1) || (sn[v] == sn[v + 1])) ? S[v + 1] : -1;

          qbt1 = vrna_exp_E_multibranch_stem(tt, si1, sj1, pf_params) *
                 qb[idx[u] - v];

          vrna_log_debug("stem(%u,%u) = %g", u, v, vrna_exp_E_multibranch_stem(tt,
                                                                               si1,
                                                                               sj1,
                                                                               pf_params));

          if (v < r->c)
            qbt1 *= score;

          /* fill A */
          if ((v + 1 < j) &&
              (hc->eval_mb(u, j - 1, v, v + 1, VRNA_DECOMP_ML_STEM_ML, hc))) {
            qq = qbt1 *
                 qm[idx[v + 1] - j + 1];

            if (sc_wrapper.decomp_ml)
              qq *= sc_wrapper.decomp_ml(u, j - 1, v, v + 1, &sc_wrapper);

            vrna_log_debug("qA[%u] = %g + %g", u, qA[u], qq);
            qA[u] += qq;
          }

          /* fill B */
          if (hc->eval_mb(u, j - 1, u, v, VRNA_DECOMP_ML_ML, hc)) {
            qq = qbt1 *
                 expMLbase[j - v - 1];

            if (sc_wrapper.red_ml)
              qq *= sc_wrapper.red_ml(u, j - 1, u, v, &sc_wrapper);

            vrna_log_debug("qB[%u] = %g + %g", u, qB[u], qq);
            qB[u] += qq;
          }
        }
      }
    }
  }

  free_sc_mb_exp(&sc_wrapper);
}


PUBLIC vrna_eqm_probs_component_t
vrna_sc_iip_get_outside_cb(vrna_fold_compound_t  *fc)
{
  if (fc) {
    vrna_eqm_probs_component_t c = (vrna_eqm_probs_component_t)vrna_alloc(sizeof(struct vrna_eqm_probs_component_s));
    c->cb           = compute_bpp_multibranch_xlink;
    c->cb_data      = (void *)get_xl_ml_outside_helpers(fc);
    c->cb_data_free = free_xl_ml_outside_helpers;
    c->cb_init      = NULL;
    c->cb_post      = NULL;

    return c;
  }

  return NULL;
}



/* Outside partition function, i.e. base pair probabilities below */
PRIVATE INLINE FLT_OR_DBL
contrib_enclosing_pair(unsigned int     i,
                       unsigned int     j,
                       vrna_md_t        *md,
                       vrna_exp_param_t *pf_params,
                       short            *S,
                       short            *S1,
                       helper_arrays_t  eqm_probs_helpers,
                       xl_ml_helper_t   *helpers)
{
  FLT_OR_DBL    q;
  unsigned int  tt;

  tt = vrna_get_ptype_md(S[j], S[i], md);

  q = vrna_exp_E_multibranch_stem(tt, S1[j - 1], S1[i + 1], pf_params);

  if (eqm_probs_helpers->sc_wrapper_mb.pair)
    q *= eqm_probs_helpers->sc_wrapper_mb.pair(i, j, &(eqm_probs_helpers->sc_wrapper_mb));

  return q;
}


PRIVATE INLINE FLT_OR_DBL
contrib_enclosed_pair(unsigned int     i,
                      unsigned int     j,
                      vrna_md_t        *md,
                      vrna_exp_param_t *pf_params,
                      short            *S,
                      short            *S1,
                      helper_arrays_t  eqm_probs_helpers,
                      xl_ml_helper_t   *helpers)
{
  FLT_OR_DBL    q;
  unsigned int  tt;

  tt = vrna_get_ptype_md(S[i], S[j], md);

  q = vrna_exp_E_multibranch_stem(tt, S1[i - 1], S1[j + 1], pf_params);

  if (eqm_probs_helpers->sc_wrapper_mb.red_stem)
    q *= eqm_probs_helpers->sc_wrapper_mb.red_stem(i, j, i, j, &(eqm_probs_helpers->sc_wrapper_mb));

  return q;
}


PRIVATE xl_ml_helper_t *
get_xl_ml_outside_helpers(vrna_fold_compound_t *fc)
{
  unsigned int    n;
  xl_ml_helper_t  *helpers;

  n       = fc->length;
  helpers = (xl_ml_helper_t *)vrna_alloc(sizeof(xl_ml_helper_t));

  helpers->vl       = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * (n + 2));
  helpers->wl       = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * (n + 2));
  helpers->yl       = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * (n + 2));
  helpers->tl       = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * (n + 2));
  helpers->rl       = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * (n + 2));
  helpers->sl       = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * (n + 2));
  helpers->ol       = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * (n + 2));
  helpers->ol2      = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * (n + 2));
  helpers->pl       = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * (n + 2));
  helpers->pl2      = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * (n + 2));
  helpers->prm_lmm  = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * (n + 2));
  helpers->prm_lmm2 = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * (n + 2));
  helpers->prm_lmu  = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * (n + 2));
  helpers->prm_lmu2 = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * (n + 2));
  helpers->prm_l1mu  = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * (n + 2));
  helpers->prm_l1mu2 = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * (n + 2));

  return helpers;
}


PRIVATE void
free_xl_ml_outside_helpers(void *helpers)
{
  if (helpers) {
    xl_ml_helper_t *h = (xl_ml_helper_t *)helpers;

    free(h->vl);
    free(h->wl);
    free(h->yl);
    free(h->tl);
    free(h->rl);
    free(h->sl);
    free(h->ol);
    free(h->ol2);
    free(h->pl);
    free(h->pl2);
    free(h->prm_lmm);
    free(h->prm_lmm2);
    free(h->prm_lmu);
    free(h->prm_lmu2);
    free(h->prm_l1mu);
    free(h->prm_l1mu2);

    free(helpers);
  }
}


PRIVATE void
xl_fill_prms(vrna_fold_compound_t *fc,
             unsigned int         l,
             helper_arrays_t      eqm_probs_helpers,
             xl_ml_helper_t       *helpers)
{
  short                         *S, *S1;
  unsigned int                  i, j, n, u, v, tt, overlap, turn;
  int                           *my_iindx, ij, lj1;
  FLT_OR_DBL                    qbt1, qbt2, temp, temp2, *probs, *qb, *qm, *scale,
                                *expMLbase, expMLclosing, *prm_lmm, *prm_lmm2,
                                *prm_lmu, *prm_lmu2, *prm_l1mu, *prm_l1mu2, 
                                *vl, *wl, *swp;
  vrna_exp_param_t              *pf_params;
  vrna_md_t                     *md;
  vrna_hc_t                     *hc;
  struct vrna_crosslink_data_s  *dat;
  interact_region_s             *r;
  struct sc_mb_exp_dat          *sc_wrap;


  n             = fc->length;
  S             = fc->sequence_encoding2;
  S1            = fc->sequence_encoding;
  my_iindx      = fc->iindx;
  pf_params     = fc->exp_params;
  md            = &(pf_params->model_details);
  turn          = md->min_loop_size;
  qb            = fc->exp_matrices->qb;
  qm            = fc->exp_matrices->qm;
  probs         = fc->exp_matrices->probs;
  scale         = fc->exp_matrices->scale;
  expMLbase     = fc->exp_matrices->expMLbase;
  expMLclosing  = pf_params->expMLclosing;
  hc            = fc->hc;
  dat           = fc->sc->iip;
  sc_wrap       = &(eqm_probs_helpers->sc_wrapper_mb);

  /* swap prm_lmu's from previous l */
  swp                 = helpers->prm_lmu;
  helpers->prm_lmu    = helpers->prm_l1mu;
  helpers->prm_l1mu   = swp;
  swp                 = helpers->prm_lmu2;
  helpers->prm_lmu2   = helpers->prm_l1mu2;
  helpers->prm_l1mu2  = swp;

  vl        = helpers->vl;
  wl        = helpers->wl;
  prm_lmm   = helpers->prm_lmm;
  prm_lmm2  = helpers->prm_lmm2;
  prm_lmu   = helpers->prm_lmu;
  prm_lmu2  = helpers->prm_lmu2;
  prm_l1mu  = helpers->prm_l1mu;
  prm_l1mu2 = helpers->prm_l1mu2;

  /* reset (helper) arrays */
  vl        = memset(vl, 0, sizeof(FLT_OR_DBL) * (n + 2));
  wl        = memset(wl, 0, sizeof(FLT_OR_DBL) * (n + 2));
  prm_lmm   = memset(prm_lmm, 0, sizeof(FLT_OR_DBL) * (n + 2));
  prm_lmm2  = memset(prm_lmm2, 0, sizeof(FLT_OR_DBL) * (n + 2));
  prm_lmu   = memset(prm_lmu, 0, sizeof(FLT_OR_DBL) * (n + 2));
  prm_lmu2  = memset(prm_lmu2, 0, sizeof(FLT_OR_DBL) * (n + 2));

  /* precompute W[l,v] for V[l,j] */
  if ((l < n) &&
      (dat->overlap5[l]) &&
      (dat->overlap5[l + 1])) {
    for (v = l + turn + 2; v < n; v++) {

      r = NULL;

      /* find out interval pair (I,I') where l - turn - 2, ... l + 1 \in I and v \in I' */
      for (size_t cnt = 0; cnt < vrna_array_size(dat->overlap5[l]); cnt++)
        if ((dat->overlap5[l][cnt]->b > l) &&
            (dat->overlap5[l][cnt]->c <= v) &&
            (dat->overlap5[l][cnt]->d > v)) {
          r = dat->overlap5[l][cnt];
          break;
        }

      if (r) {
        u = l + 1;

        if (hc->eval_mb(u, v, u, v, VRNA_DECOMP_ML_STEM, hc)) {
          qbt1 = qb[my_iindx[u] - v] *
                 contrib_enclosed_pair(u, v, md, pf_params, S, S1, eqm_probs_helpers, helpers);

          if (sc_wrap->red_stem)
            qbt1 *= sc_wrap->red_stem(u, v, u, v, sc_wrap);

          /* no score, since l < r->b, u = l + 1 <= r->b */
          wl[v] += qbt1;
        }

        for (u = l + 2; (u <= r->b) && (u + turn < v); u++) {
          temp = qb[my_iindx[u] - v] *
                 contrib_enclosed_pair(u, v, md, pf_params, S, S1, eqm_probs_helpers, helpers);

          if (hc->eval_mb(l + 1, v, u - 1, u, VRNA_DECOMP_ML_ML_STEM, hc)) {
            qbt1 = qm[my_iindx[l + 1] - u + 1];

            if (sc_wrap->decomp_ml)
              qbt1 *= sc_wrap->decomp_ml(l + 1, v, u - 1, u, sc_wrap);

            wl[v] += temp *
                    qbt1;
          }

          if (hc->eval_mb(l + 1, v, u, v, VRNA_DECOMP_ML_STEM, hc)) {
            qbt1 = expMLbase[u - l - 1];

            if (sc_wrap->red_stem)
              qbt1 *= sc_wrap->red_stem(l + 1, v, u, v, sc_wrap);

            wl[v] += temp *
                     qbt1;
          }
        }

        for (; (u < r->c) && (u + turn < v); u++) {
          temp = qb[my_iindx[u] - v] *
                 contrib_enclosed_pair(u, v, md, pf_params, S, S1, eqm_probs_helpers, helpers) *
                 r->exp_score;

          if (hc->eval_mb(l + 1, v, u - 1, u, VRNA_DECOMP_ML_ML_STEM, hc)) {
            qbt1 = qm[my_iindx[l + 1] - u + 1];

            if (sc_wrap->decomp_ml)
              qbt1 *= sc_wrap->decomp_ml(l + 1, v, u - 1, u, sc_wrap);

            wl[v] += temp *
                     qbt1;
          }

          if (hc->eval_mb(l + 1, v, u, v, VRNA_DECOMP_ML_STEM, hc)) {
            qbt1 = expMLbase[u - l - 1];

            if (sc_wrap->red_stem)
              qbt1 *= sc_wrap->red_stem(l + 1, v, u, v, sc_wrap);

            wl[v] += temp *
                     qbt1;
          }
        }
      }
    }
  }

  /* pre-compute V[l,j] for prm_lmm and prm_lmm2 */
  for (j = l + turn + 3; j <= n; j++) {
    if ((dat->overlap5[l]) &&
        (dat->overlap3[j])) {
      /* find out interval pair (I,I') with l \in I and j \in I' if any */
      r = NULL;

      for (size_t cnt = 0; cnt < vrna_array_size(dat->overlap5[l]); cnt++)
        if ((dat->overlap5[l][cnt]->b > l) && /* l + 1 must be in I */
            (dat->overlap5[l][cnt]->c < j) &&
            (dat->overlap5[l][cnt]->d >= j)) {
          r = dat->overlap5[l][cnt];
          break;
        }

      if (r) {
        for (v = r->c; v < j; v++) {
          temp = wl[v];

          if ((v + 1 < j) &&
              (hc->eval_mb(l + 1, j - 1, v, v + 1, VRNA_DECOMP_ML_ML_ML, hc))) {
            qbt1 = qm[my_iindx[v + 1] - j + 1];

            if (sc_wrap->decomp_ml)
              qbt1 *= sc_wrap->decomp_ml(l + 1, j - 1, v, v + 1, sc_wrap);

            vl[j] += qbt1 *
                     temp;
          }

          if (hc->eval_mb(l + 1, j - 1, l + 1, v, VRNA_DECOMP_ML_ML, hc)) {
            qbt1 = expMLbase[j - v - 1];

            if (sc_wrap->red_ml)
              qbt1 *= sc_wrap->red_ml(l + 1, j - 1, l + 1, v, sc_wrap);

            vl[j] += qbt1 *
                     temp;
          }
        }
      }
    }
  }

  /* finally, compute prm_lmm and prm_lmm2 */
  for (i = 1; i + turn + 1 < l; i++) {
    for (j = l + turn + 3; j <= n; j++) {
      ij = my_iindx[i] - j;

      if (hc->eval_mb(i, j, i + 1, j - 1, VRNA_DECOMP_PAIR_ML, hc)) {
        /* compute closing pair contribution */
        temp = probs[ij] *
               contrib_enclosing_pair(i, j, md, pf_params, S, S1, eqm_probs_helpers, helpers);

        r = NULL;

        /* check whether there is any interval pair (I,I') supported by (i,j) */
        if ((dat->overlap5[i]) &&
            (dat->overlap3[j])) {
          for (size_t cnt = 0; cnt < vrna_array_size(dat->overlap5[i]); cnt++)
            if ((dat->overlap5[i][cnt]->c <= j) &&
                (dat->overlap5[i][cnt]->d >= j)) {
              r = dat->overlap5[i][cnt];
              break;
            }
        }

        qbt1 = qbt2 = 0.;

        if (r == NULL) {
          if (hc->eval_mb(i + 1, j - 1, l, l + 1, VRNA_DECOMP_ML_ML_ML, hc)) {
            qbt1 = qm[my_iindx[l + 1] - j + 1];

            if (sc_wrap->decomp_ml)
              qbt1 *= sc_wrap->decomp_ml(i + 1, j - 1, l, l + 1, sc_wrap);

            qbt2 = qbt1; /* for prm_lmm2 */
          }
        } else if ((r->b <= l) || (r->c <= l + 1) || (r->c == j)) {
          if (hc->eval_mb(i + 1, j - 1, l, l + 1, VRNA_DECOMP_ML_ML_ML, hc)) {
            qbt1 = qm[my_iindx[l + 1] - j + 1];

            if (sc_wrap->decomp_ml)
              qbt1 *= sc_wrap->decomp_ml(i + 1, j - 1, l, l + 1, sc_wrap);

            qbt2 = qbt1; /* for prm_lmm2 */

            if (r->c > l)
              qbt2 *= r->exp_score;

            qbt1 *= r->exp_score;
          }
        } else {
          if (hc->eval_mb(i + 1, j - 1, l, l + 1, VRNA_DECOMP_ML_ML_ML, hc)) {
            qbt1 = vl[j];

            FLT_OR_DBL ulj = 0.;

            if (hc->eval_mb(l + 1, j - 1, l + 1, r->c - 1, VRNA_DECOMP_ML_ML, hc)) {
              temp2 = qm[my_iindx[l + 1] - r->c + 1] *
                      expMLbase[j - r->c] *
                      r->exp_score;

              if (sc_wrap->red_ml)
                temp2 *= sc_wrap->red_ml(l + 1, j - 1, l + 1, r->c - 1, sc_wrap);

              ulj += temp2;
            }

            if (hc->eval_mb(l + 1, j - 1, r->c, j - 1, VRNA_DECOMP_ML_ML, hc)) {
              temp2 = qm[my_iindx[r->c] - j + 1] *
                      expMLbase[r->c - l - 1] *
                      r->exp_score;

              if (sc_wrap->red_ml)
                temp2 *= sc_wrap->red_ml(l + 1, j - 1, r->c, j - 1, sc_wrap);

              ulj += temp2;
            }

            if (hc->eval_mb(l + 1, j - 1, r->c - 1, r->c, VRNA_DECOMP_ML_ML_ML, hc)) {
              temp2 = qm[my_iindx[l + 1] - r->c + 1] *
                      qm[my_iindx[r->c] - j + 1] *
                      r->exp_score;

              if (sc_wrap->decomp_ml)
                temp2 *= sc_wrap->decomp_ml(l + 1, j - 1, r->c - 1, r->c, sc_wrap);

              ulj += temp2;
            }

            
            qbt1 += ulj;

            qbt2 = qbt1; /* for prm_lmm2 */
          }
        }

        prm_lmm[i]  += qbt1 *
                       temp;
        prm_lmm2[i] += qbt2 *
                       temp;
      }
    }

    prm_lmm[i] *= expMLclosing;
    prm_lmm2[i] *= expMLclosing;
  }

  /* and the simpler cases, prm_lmu and prm_lmu2 */
  if (l < n) {
    for (i = 1; i + turn + 1 < l; i++) {
      /* all cases with l + 1 < j <= n */
      if (hc->eval_mb(i + 1, l + 1, i + 1, l, VRNA_DECOMP_ML_ML, hc)) {
        qbt1 = prm_l1mu[i] *
               expMLbase[1];
        qbt2 = prm_l1mu2[i] *
               expMLbase[1];

        if (sc_wrap->red_ml) {
          temp = sc_wrap->red_ml(i + 1, l + 1, i + 1, l, sc_wrap);
          qbt1 *= temp;
          qbt2 *= temp;
        }

        prm_lmu[i] += qbt1;
        prm_lmu2[i] += qbt2;
      }

      j = l + 1;
      ij = my_iindx[i] - j;

      if (hc->eval_mb(i, j, i + 1, j - 1, VRNA_DECOMP_PAIR_ML, hc)) {
        temp = probs[ij] *
               contrib_enclosing_pair(i, j, md, pf_params, S, S1, eqm_probs_helpers, helpers) *
               expMLclosing;

        r = NULL;

        if ((dat->overlap5[i]) &&
            (dat->overlap3[j])) {
          for (size_t cnt = 0; cnt < vrna_array_size(dat->overlap5[i]); cnt++)
            if ((dat->overlap5[i][cnt]->c <= j) &&
                (dat->overlap5[i][cnt]->d >= j)) {
              r = dat->overlap5[i][cnt];
              break;
            }
        }

        if (hc->eval_mb(i + 1, j - 1, i + 1, l, VRNA_DECOMP_ML_ML, hc)) {
          if (r) {
            prm_lmu[i] += temp *
                          r->exp_score;
            if (r->c > l) {
              prm_lmu2[i] += temp *
                             r->exp_score;
            } else {
              prm_lmu2[i] += temp;
            }
          } else {
            prm_lmu[i] += temp;
            prm_lmu2[i] += temp;
          }
        }
      }
    }
  }
}


PRIVATE void
xl_fill_y(vrna_fold_compound_t *fc,
          unsigned int         l,
          helper_arrays_t      eqm_probs_helpers,
          xl_ml_helper_t       *helpers)
{
  short                         *S, *S1;
  unsigned int                  i, k, u, v, turn, n;
  int                           *my_iindx;
  FLT_OR_DBL                    *ol, *ol2, *rl, temp, qbt1, qbt2, *qm, *qb,
                                *expMLbase, *scale, *prm_lmu, *prm_lmu2, *yl,
                                olu, ol2u, ylk;
  vrna_exp_param_t              *pf_params;
  vrna_md_t                     *md;
  vrna_hc_t                     *hc;
  struct vrna_crosslink_data_s  *dat;
  interact_region_s             *r;
  struct sc_mb_exp_dat          *sc_wrap;



  n         = fc->length;
  S             = fc->sequence_encoding2;
  S1            = fc->sequence_encoding;
  my_iindx  = fc->iindx;
  pf_params = fc->exp_params;
  md        = &(pf_params->model_details);
  turn      = md->min_loop_size;
  qb        = fc->exp_matrices->qb;
  qm        = fc->exp_matrices->qm;
  expMLbase = fc->exp_matrices->expMLbase;
  scale     = fc->exp_matrices->scale;
  hc        = fc->hc;
  dat       = fc->sc->iip;

  ol        = helpers->ol;;
  ol2       = helpers->ol2;
  rl        = helpers->rl;
  yl        = helpers->yl;
  prm_lmu   = helpers->prm_lmu;
  prm_lmu2  = helpers->prm_lmu2;
  sc_wrap   = &(eqm_probs_helpers->sc_wrapper_mb);

  yl = memset(yl, 0, sizeof(FLT_OR_DBL) * (n + 2));

  if (dat->overlap3[l]) {
    /* reset helper arrays */
    ol  = memset(ol, 0, sizeof(FLT_OR_DBL) * (n + 2));
    ol2 = memset(ol2, 0, sizeof(FLT_OR_DBL) * (n + 2));
    rl  = memset(rl, 0, sizeof(FLT_OR_DBL) * (n + 2));

    /* fill O[l,u] and O2[l,u] */
    for (u = 2; u + turn + 2 < l; u++) {
      /* find interval pair (I,I') with u \in I and l \in I', if any */
      r = NULL;

      for (size_t cnt = 0; cnt < vrna_array_size(dat->overlap3[l]); cnt++)
        if ((dat->overlap3[l][cnt]->a < u) &&
            (dat->overlap3[l][cnt]->b >= u) &&
            (dat->overlap3[l][cnt]->d > l)) {
          r = dat->overlap3[l][cnt];
          break;
        }

      if (r) {
        olu = ol2u = 0.;

        for (i = r->a; i < u; i++) {
          qbt1 = prm_lmu[i];
          qbt2 = prm_lmu2[i];

          if ((i + 1 < u) &&
              (hc->eval_mb(i + 1, l, u - 1, u, VRNA_DECOMP_ML_ML_ML, hc))) {
            temp = qm[my_iindx[i + 1] - u + 1];

            if (sc_wrap->decomp_ml)
              temp *= sc_wrap->decomp_ml(i + 1, l, u - 1, u, sc_wrap);

            olu += qbt1 *
                   temp;
            ol2u += qbt2 *
                    temp;
          }

          if (hc->eval_mb(i + 1, l, u, l, VRNA_DECOMP_ML_ML, hc)) {
            temp = expMLbase[u - i - 1];

            if (sc_wrap->red_ml)
              temp *= sc_wrap->red_ml(i + 1, l, u, l, sc_wrap);

            olu += qbt1 *
                   temp;
            ol2u += qbt2 *
                    temp;
          }
        }

        ol[u]   += olu;
        ol2[u]  += ol2u;
      }
    }

    /*
     *  fill R[l,v]
     */
    for (u = 2; u + 2 * turn + 2 < l; u++) {
      /* find out, whether there is any interval pair (I,I') with u \in I and l \in I' */
      r = NULL;

      for (size_t cnt = 0; cnt < vrna_array_size(dat->overlap3[l]); cnt++)
        if ((dat->overlap3[l][cnt]->a < u) &&
            (dat->overlap3[l][cnt]->b >= u) &&
            (dat->overlap3[l][cnt]->d > l)) {
          r = dat->overlap3[l][cnt];
          break;
        }

      if (r) {
        olu   = ol[u];
        ol2u  = ol2[u];

        for (v = u + turn + 1; v + turn + 1 < l; v++) {
          if ((r->b < v) &&
              (hc->eval_mb(u, v, u, v, VRNA_DECOMP_ML_STEM, hc))) {
            temp = qb[my_iindx[u] - v] *
                   contrib_enclosed_pair(u, v, md, pf_params, S, S1, eqm_probs_helpers, helpers);

            if (r->c <= v)
              temp *= ol2u;
            else
              temp *= olu;

            if (sc_wrap->red_stem)
              temp *= sc_wrap->red_stem(u, v, u, v, sc_wrap);

            rl[v] += temp;
          }
        }
      }
    }

    /* Finally, compute Y[l,k] */
    if ((l < n) &&
        (dat->overlap3[l + 1])) {
      for (k = turn + 4; k + turn < l; k++) {
        if (!(hc->eval_mb(k, l, k, l, VRNA_DECOMP_ML_STEM, hc)))
          continue;

        /* find out whether there is any interval I' with k \in I' and l \in I' */
        r = NULL;

        for (size_t cnt = 0; cnt < vrna_array_size(dat->overlap3[l + 1]); cnt++)
          if ((dat->overlap3[l + 1][cnt]->b + 1 < k) &&
              (dat->overlap3[l + 1][cnt]->d > k)) {
            r = dat->overlap3[l + 1][cnt];
            break;
          }

        if (r == NULL)
          continue;

        ylk = 0.;

        for (v = turn + 3; v < k; v++) {
          temp = rl[v];

          if (temp != 0.) {
            if ((v + 1 < k) &&
                (hc->eval_mb(v + 1, l, k - 1, k, VRNA_DECOMP_ML_ML_ML, hc))) {
              qbt1 = qm[my_iindx[v + 1] - k + 1];

              if (sc_wrap->decomp_ml)
                qbt1 *= sc_wrap->decomp_ml(v + 1, l, k - 1, k, sc_wrap);

              ylk += qbt1 *
                     temp;
            }

            if (hc->eval_mb(v + 1, l, k, l, VRNA_DECOMP_ML_ML, hc)) {
              qbt1 = expMLbase[k - v - 1];

              if (sc_wrap->red_ml)
                qbt1 *= sc_wrap->red_ml(v + 1, l, k, l, sc_wrap);

              ylk += qbt1 *
                     temp;
            }
          }
        }

        yl[k] += ylk;
      }
    }
  }
}


PRIVATE void
xl_fill_t(vrna_fold_compound_t *fc,
          unsigned int         l,
          helper_arrays_t      eqm_probs_helpers,
          xl_ml_helper_t       *helpers)
{
  short                         *S, *S1;
  unsigned int                  i, k, u, v, turn, n;
  int                           *my_iindx;
  FLT_OR_DBL                    *pl, *pl2, *sl, temp, qbt1, qbt2, *qm, *qb,
                                *expMLbase, *scale, *tl, *prm_lmm, *prm_lmm2,
                                plu, pl2u, tlk;
  vrna_exp_param_t              *pf_params;
  vrna_md_t                     *md;
  vrna_hc_t                     *hc;
  struct vrna_crosslink_data_s  *dat;
  interact_region_s             *r;
  struct sc_mb_exp_dat          *sc_wrap;



  n         = fc->length;
  S         = fc->sequence_encoding2;
  S1        = fc->sequence_encoding;
  my_iindx  = fc->iindx;
  pf_params = fc->exp_params;
  md        = &(pf_params->model_details);
  turn      = md->min_loop_size;
  qb        = fc->exp_matrices->qb;
  qm        = fc->exp_matrices->qm;
  expMLbase = fc->exp_matrices->expMLbase;
  scale     = fc->exp_matrices->scale;
  hc        = fc->hc;
  dat       = fc->sc->iip;

  pl        = helpers->pl;
  pl2       = helpers->pl2;
  sl        = helpers->sl;
  tl        = helpers->tl;
  prm_lmm   = helpers->prm_lmm;
  prm_lmm2  = helpers->prm_lmm2;
  sc_wrap   = &(eqm_probs_helpers->sc_wrapper_mb);

  tl  = memset(tl, 0, sizeof(FLT_OR_DBL) * (n + 2));

  if (dat->overlap3[l]) {
    /* reset arrays */
    pl  = memset(pl, 0, sizeof(FLT_OR_DBL) * (n + 2));
    pl2 = memset(pl2, 0, sizeof(FLT_OR_DBL) * (n + 2));
    sl  = memset(sl, 0, sizeof(FLT_OR_DBL) * (n + 2));

    /* fill O[l,u] and O2[l,u] */
    for (u = 2; u + turn + 2 < l; u++) {
      /* find interval pair (I,I') with u \in I and l \in I', if any */
      r = NULL;

      for (size_t cnt = 0; cnt < vrna_array_size(dat->overlap3[l]); cnt++)
        if ((dat->overlap3[l][cnt]->a < u) &&
            (dat->overlap3[l][cnt]->b >= u) &&
            (dat->overlap3[l][cnt]->d > l)) {
          r = dat->overlap3[l][cnt];
          break;
        }

      if (r) {
        plu = pl2u = 0.;

        for (i = r->a; i < u; i++) {
          qbt1 = prm_lmm[i];
          qbt2 = prm_lmm2[i];

          if ((i + 1 < u) &&
              (hc->eval_mb(i + 1, l, u - 1, u, VRNA_DECOMP_ML_ML_ML, hc))) {
            temp = qm[my_iindx[i + 1] - u + 1];

            if (sc_wrap->decomp_ml)
              temp *= sc_wrap->decomp_ml(i + 1, l, u - 1, u, sc_wrap);

            plu += qbt1 *
                   temp;
            pl2u += qbt2 *
                    temp;
          }

          if (hc->eval_mb(i + 1, l, u, l, VRNA_DECOMP_ML_ML, hc)) {
            temp = expMLbase[u - i - 1];

            if (sc_wrap->red_ml)
              temp *= sc_wrap->red_ml(i + 1, l, u, l, sc_wrap);

            plu += qbt1 *
                   temp;
            pl2u += qbt2 *
                    temp;
          }
        }

        pl[u]   += plu;
        pl2[u]  += pl2u;
      }
    }

    /*
     *  fill S[l,v]
     */
    for (u = 2; u + 2 * turn + 2 < l; u++) {
      /* find out, whether there is any interval pair (I,I') with u \in I and l \in I' */
      r = NULL;

      for (size_t cnt = 0; cnt < vrna_array_size(dat->overlap3[l]); cnt++)
        if ((dat->overlap3[l][cnt]->d > l) &&
            (dat->overlap3[l][cnt]->a < u) &&
            (dat->overlap3[l][cnt]->b >= u)) {
          r = dat->overlap3[l][cnt];
          break;
        }

      if (r) {
        plu   = pl[u];
        pl2u  = pl2[u];

        for (v = u + turn + 1; v + turn + 1 < l; v++) {
          if ((r->b < v) &&
              (hc->eval_mb(u, v, u, v, VRNA_DECOMP_ML_STEM, hc))) {
            temp = qb[my_iindx[u] - v] *
                   contrib_enclosed_pair(u, v, md, pf_params, S, S1, eqm_probs_helpers, helpers);

            if (r->c <= v)
              temp *= pl2u;
            else
              temp *= plu;

            if (sc_wrap->red_stem)
              temp *= sc_wrap->red_stem(u, v, u, v, sc_wrap);

            sl[v] += temp;
          }
        }
      }
    }

    /* Finally, compute T[l,k] */
    if ((l < n) &&
        (dat->overlap3[l + 1])) {
      for (k = turn + 4; k + turn < l; k++) {
        if (!(hc->eval_mb(k, l, k, l, VRNA_DECOMP_ML_STEM, hc)))
          continue;

        /* find out whether there is any interval I' with k \in I' and l \in I' */
        r = NULL;

        for (size_t cnt = 0; cnt < vrna_array_size(dat->overlap3[l + 1]); cnt++)
          if ((dat->overlap3[l + 1][cnt]->b + 1 < k) &&
              (dat->overlap3[l + 1][cnt]->d > k)) {
            r = dat->overlap3[l + 1][cnt];
            break;
          }

        if (r == NULL)
          continue;

        tlk = 0;

        for (v = turn + 3; v < k; v++) {
          temp = sl[v];

          if (temp != 0.) {
            if ((v + 1 < k) &&
                (hc->eval_mb(v + 1, l, k - 1, k, VRNA_DECOMP_ML_ML_ML, hc))) {
              qbt1 = qm[my_iindx[v + 1] - k + 1];

              if (sc_wrap->decomp_ml)
                qbt1 *= sc_wrap->decomp_ml(v + 1, l, k - 1, k, sc_wrap);

              tlk += qbt1 *
                     temp;
            }

            if (hc->eval_mb(v + 1, l, k, l, VRNA_DECOMP_ML_ML, hc)) {
              qbt1 = expMLbase[k - v - 1];

              if (sc_wrap->red_ml)
                qbt1 *= sc_wrap->red_ml(v + 1, l, k, l, sc_wrap);

              tlk += qbt1 *
                     temp;
            }
          }
        }

        tl[k] += tlk;
      }
    }
  }
}


PRIVATE void
compute_bpp_multibranch_xlink(vrna_fold_compound_t *fc,
                              unsigned int         l,
                              helper_arrays_t      eqm_probs_helpers,
                              void                 *data)
{
  short                         *S, *S1;
  unsigned int                  i, k, n, turn;
  int                           *my_iindx;
  FLT_OR_DBL                    pklm, qbt1, temp, *qm, *expMLbase, *probs, *scale,
                                *prm_lmm, *prm_lmm2, *prm_lmu, *prm_lmu2, *yl, *tl;
  vrna_exp_param_t              *pf_params;
  vrna_md_t                     *md;
  vrna_hc_t                     *hc;
  struct vrna_crosslink_data_s  *dat;
  interact_region_s             *r;
  xl_ml_helper_t                *h;
  struct sc_mb_exp_dat          *sc_wrap;

  n         = fc->length;
  S         = fc->sequence_encoding2;
  S1        = fc->sequence_encoding;
  pf_params = fc->exp_params;
  md        = &(pf_params->model_details);
  turn      = md->min_loop_size;
  my_iindx  = fc->iindx;
  qm        = fc->exp_matrices->qm;
  probs     = fc->exp_matrices->probs;
  expMLbase = fc->exp_matrices->expMLbase;
  scale     = fc->exp_matrices->scale;
  hc        = fc->hc;
  h         = (xl_ml_helper_t *)data;
  dat       = fc->sc->iip;
  sc_wrap   = &(eqm_probs_helpers->sc_wrapper_mb);

  xl_fill_prms(fc, l, eqm_probs_helpers, h);
  xl_fill_y(fc, l, eqm_probs_helpers, h);
  xl_fill_t(fc, l, eqm_probs_helpers, h);

  prm_lmm   = h->prm_lmm;
  prm_lmm2  = h->prm_lmm2;
  prm_lmu   = h->prm_lmu;
  prm_lmu2  = h->prm_lmu2;
  yl        = h->yl;
  tl        = h->tl;

  for (k = 2; k + turn < l; k++) {
    if (!hc->eval_mb(k, l, k, l, VRNA_DECOMP_ML_STEM, hc))
      continue;

    pklm = 0.;

    /* 1. (k,l) is left-most stem in multibranch loop */
    for (i = 1; i < k; i++) {
      temp = prm_lmm[i];

      if (dat->overlap3[l]) {
        for (size_t cnt = 0; cnt < vrna_array_size(dat->overlap3[l]); cnt++)
          if ((dat->overlap3[l][cnt]->a <= i) &&
              (dat->overlap3[l][cnt]->b > i) &&
              (dat->overlap3[l][cnt]->b >= k))
            temp = prm_lmm2[i];
      }

      if (hc->eval_mb(i + 1, l, k, l, VRNA_DECOMP_ML_STEM, hc)) {
        pklm += expMLbase[k - i - 1] *
                temp;
      }
    }

    /* 2./3. (k,l) is right-most stem or somewhere in the middle of a multibranch loop */
    pklm += yl[k];
    pklm += tl[k];

    for (i = 1; i + turn + 2 < k; i++) {
      if (hc->eval_mb(i + 1, l, k - 1, k, VRNA_DECOMP_ML_ML_STEM, hc)) {
        /* find out whether exists an interval pair (I,I') with i \in I and l \in I' */
        if (dat->overlap3[l]) {
          r = NULL;
          for (size_t cnt = 0; cnt < vrna_array_size(dat->overlap3[l]); cnt++)
            if ((dat->overlap3[l][cnt]->a <= i) &&
                (dat->overlap3[l][cnt]->b >= i) &&
                (dat->overlap3[l][cnt]->d > l) &&
                (dat->overlap3[l][cnt]->b < k)) {
              r = dat->overlap3[l][cnt];
              break;
            }

          if (r) {
            if (i < r->b) {
              temp = 0.;
              
              if (hc->eval_mb(i + 1, k - 1, i + 1, r->b, VRNA_DECOMP_ML_ML, hc)) {
                qbt1 = qm[my_iindx[i + 1] - r->b] *
                       expMLbase[k - r->b - 1];

                if (sc_wrap->red_ml)
                  qbt1 *= sc_wrap->red_ml(i + 1, k - 1, i + 1, r->b, sc_wrap);

                temp += qbt1;
              }

              if (r->b + 1 < k) {
                if (hc->eval_mb(i + 1, k - 1, r->b + 1, k - 1, VRNA_DECOMP_ML_ML, hc)) {
                  qbt1 = qm[my_iindx[r->b + 1] - k + 1] *
                         expMLbase[r->b - i];

                  if (sc_wrap->red_ml)
                    qbt1 *= sc_wrap->red_ml(i + 1, k - 1, r->b + 1, k - 1, sc_wrap);

                  temp += qbt1;
                }

                if (hc->eval_mb(i + 1, k - 1, r->b, r->b + 1, VRNA_DECOMP_ML_ML_ML, hc)) {
                  qbt1 = qm[my_iindx[i + 1] - r->b] *
                         qm[my_iindx[r->b + 1] - k + 1];

                  if (sc_wrap->decomp_ml)
                    qbt1 *= sc_wrap->decomp_ml(i + 1, k - 1, r->b, r->b + 1, sc_wrap);
  
                  temp += qbt1;
                }
              }

              pklm += temp *
                      (prm_lmu[i] + prm_lmm[i]);

              continue;
            } else if (i == r->b) {
              temp = (prm_lmu[i] + prm_lmm[i]) *
                     qm[my_iindx[i + 1] - k + 1];

              if (sc_wrap->decomp_ml)
                temp *= sc_wrap->decomp_ml(i + 1, l, k - 1, k, sc_wrap);

              pklm += temp;

              continue;
            }
          }
        }

        temp = (prm_lmu2[i] + prm_lmm2[i]) *
               qm[my_iindx[i + 1] - k + 1];

        if (sc_wrap->decomp_ml)
          temp *= sc_wrap->decomp_ml(i + 1, l, k - 1, k, sc_wrap);

        pklm += temp;
      }
    }

    pklm *= contrib_enclosed_pair(k, l, md, pf_params, S, S1, eqm_probs_helpers, h) *
            scale[2];

    probs[my_iindx[k] - l] += pklm;


    vrna_log_debug("p_mb[%u,%u] = %g", k, l, pklm);

  }
}
