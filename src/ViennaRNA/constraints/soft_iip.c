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

#include "ViennaRNA/constraints/multibranch_sc.inc"
#include "ViennaRNA/constraints/multibranch_sc_pf.inc"


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
  unsigned int  a;          /* start of 5' interval */
  unsigned int  b;          /* end of 5' interval */
  unsigned int  c;          /* start of 3' interval */
  unsigned int  d;          /* end of 3' interval */
  int           score;      /* interval score */
  FLT_OR_DBL    exp_score;  /* interval score as Boltzmann factor */
} interact_region_s;


struct vrna_crosslink_data_s {
  vrna_array(interact_region_s)   regions;  /* list of interacting regions */
  vrna_array(interact_region_s **) overlap5;  /* for each nucleotide position i we store the corresponding interval pair (I, I') where I overlaps with i */
  vrna_array(interact_region_s **) overlap3;  /* for each nucleotide position i we store the corresponding interval pair (I, I') where I' overlaps with i */
  int         *A;
  int         *B;
  FLT_OR_DBL  *qA;
  FLT_OR_DBL  *qB;
};


#if 1
typedef struct {
  FLT_OR_DBL  *vl;  /* V[l,j] */
  FLT_OR_DBL  *wl;  /* W[l,v] */
  FLT_OR_DBL  *yl;  /* Y[l,k] */
  FLT_OR_DBL  *tl;  /* T[l,k] */
  FLT_OR_DBL  *prm_lmm;
  FLT_OR_DBL  *prm_lmm2;
  FLT_OR_DBL  *prm_lmu;
  FLT_OR_DBL  *prm_lmu2;
  struct sc_mb_exp_dat  sc_wrapper_mb;
} xl_ml_helper_t;

#else

typedef struct {
  FLT_OR_DBL  *prm_lmu;
  FLT_OR_DBL  *prm_lmu2;
  FLT_OR_DBL  *prm_lmm;
  FLT_OR_DBL  *prm_lmm2;
  FLT_OR_DBL  *r_lmu;
  FLT_OR_DBL  *r_lmm;
  FLT_OR_DBL  *s_lmm;
  struct sc_mb_exp_dat  sc_wrapper_mb;
} xl_ml_helper_t;

#endif

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


PRIVATE interact_region_s *
parse_crosslink_data(FILE *fp);


INLINE PRIVATE int
xlink_hp_en(vrna_fold_compound_t  *fc,
            int                   i,
            int                   j,
            int                   k,
            int                   l,
            void                  *data);


INLINE PRIVATE int
xlink_il_en(vrna_fold_compound_t  *fc,
            int                   i,
            int                   j,
            int                   k,
            int                   l,
            void                  *data);


INLINE PRIVATE FLT_OR_DBL
xlink_hp_exp_en(vrna_fold_compound_t  *fc,
                int                   i,
                int                   j,
                int                   k,
                int                   l,
                void                  *data);


INLINE PRIVATE FLT_OR_DBL
xlink_il_exp_en(vrna_fold_compound_t  *fc,
                int                   i,
                int                   j,
                int                   k,
                int                   l,
                void                  *data);


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
             xl_ml_helper_t       *helpers);


PRIVATE void
xl_fill_rs(vrna_fold_compound_t  *fc,
           unsigned int          l,
           xl_ml_helper_t        *helpers);


PRIVATE void
xl_fill_s(vrna_fold_compound_t *fc,
          unsigned int         l,
          xl_ml_helper_t       *helpers);


PRIVATE void
compute_bpp_multibranch_xlink(vrna_fold_compound_t  *fc,
                              unsigned int          l,
                              FLT_OR_DBL            *Qmax,
                              size_t                *ov,
                              void                  *helpers);

PRIVATE xl_ml_helper_t *
get_xl_ml_outside_helpers(vrna_fold_compound_t *fc);


PRIVATE void
free_xl_ml_outside_helpers(void *helpers);


PUBLIC int
vrna_sc_iip(vrna_fold_compound_t *fc)
{
  int                           ret;
  FILE                          *fp;
  struct vrna_crosslink_data_s  *data;

  ret = 0;

  if (fc) {
    ret   = 1;
    data  = (struct vrna_crosslink_data_s *)vrna_alloc(sizeof(struct vrna_crosslink_data_s));

    data->overlap5 = NULL;
    data->overlap3 = NULL;

    vrna_array_init_size(data->regions, 42);

    if (fc->sc == NULL)
      vrna_sc_init(fc);

    fc->sc->iip = data;
    fc->params->model_details.uniq_ML = 1;

    ret &= vrna_sc_multi_cb_add(fc,
                                xlink_hp_en,
                                xlink_hp_exp_en,
                                (void *)data,
                                NULL,
                                NULL,
                                VRNA_DECOMP_PAIR_HP);
    ret &= vrna_sc_multi_cb_add(fc,
                                xlink_il_en,
                                xlink_il_exp_en,
                                (void *)data,
                                NULL,
                                NULL,
                                VRNA_DECOMP_PAIR_IL);
  }

  return ret;
}


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
      (fc->sc) &&
      (fc->sc->iip)) {
    ret   = 1;
    data  = fc->sc->iip;

    vrna_array_append(data->regions, ((interact_region_s){
      .a      = i,
      .b      = j,
      .c      = k,
      .d      = l,
      .score  = (int)(energy * 100.)
    }));
  }

  return ret;
}


PUBLIC int
vrna_sc_add_iip_file(vrna_fold_compound_t  *fc,
                     const char            *filename)
{
  int   ret = 0;
  FILE  *fp;

  if ((fc) &&
      (filename)) {
    if (!(fp = fopen(filename, "r"))) {
      vrna_log_warning("Can't open crosslink-data file \"%s\"!",
                       filename);
      return 0;
    }

    vrna_array(interact_region_s) regions = parse_crosslink_data(fp);

    fclose(fp);

    if (regions) {
      if ((fc->sc == NULL) ||
          (fc->sc->iip == NULL))
        ret = vrna_sc_iip(fc);
      else
        ret = 1;

      for (size_t i = 0; i < vrna_array_size(regions); i++)
        ret &= vrna_sc_add_iip(fc,
                               regions[i].a,
                               regions[i].b,
                               regions[i].c,
                               regions[i].d,
                               (double)regions[i].score / 100.,
                               0);
    }

    vrna_array_free(regions);
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


INLINE PRIVATE int
xlink_hp_en(vrna_fold_compound_t  *fc,
            int                   i,
            int                   j,
            int                   k,
            int                   l,
            void                  *data)
{
  struct vrna_crosslink_data_s  *dat = (struct vrna_crosslink_data_s *)data;
  interact_region_s *r;

  /* go through each interval that overlaps with position i to check whether j connects both segments */
  if (dat->overlap5[i])
    for (size_t cnt = 0; cnt < vrna_array_size(dat->overlap5[i]); cnt++) {
      r = dat->overlap5[i][cnt];
      if ((r->c <= j) &&
          (r->d >= j))
        return r->score;
    }

  return 0;
}


INLINE PRIVATE FLT_OR_DBL
xlink_hp_exp_en(vrna_fold_compound_t  *fc,
                int                   i,
                int                   j,
                int                   k,
                int                   l,
                void                  *data)
{
  struct vrna_crosslink_data_s  *dat = (struct vrna_crosslink_data_s *)data;
  interact_region_s *r;

  /* go through each interval that overlaps with position i to check whether j connects both segments */
  if (dat->overlap5[i]) {
    for (size_t cnt = 0; cnt < vrna_array_size(dat->overlap5[i]); cnt++) {
      r = dat->overlap5[i][cnt];
      if ((r->c <= j) &&
          (r->d >= j))
        return r->exp_score;
    }
  }

  return 1.0;
}


INLINE PRIVATE int
xlink_il_en(vrna_fold_compound_t  *fc,
            int                   i,
            int                   j,
            int                   k,
            int                   l,
            void                  *data)
{
  struct vrna_crosslink_data_s  *dat    = (struct vrna_crosslink_data_s *)data;
  interact_region_s *r;

  /*  go through each interval that overlaps with position i to check whether
   *  k does not overlap with the 5' interval and j connects both segments */
  if (dat->overlap5[i])
    for (size_t cnt = 0; cnt < vrna_array_size(dat->overlap5[i]); cnt++) {
      r = dat->overlap5[i][cnt];

      if ((j >= r->c) &&
          (j <= r->d) && /* done checking whether (i,j) connects interval pair */
          ((l < r->c) || (k > r->b))) {
        return r->score;
      }
    }

  return 0;
}


INLINE PRIVATE FLT_OR_DBL
xlink_il_exp_en(vrna_fold_compound_t  *fc,
                int                   i,
                int                   j,
                int                   k,
                int                   l,
                void                  *data)
{
  struct vrna_crosslink_data_s  *dat    = (struct vrna_crosslink_data_s *)data;
  interact_region_s *r;

  /*  go through each interval that overlaps with position i to check whether
   *  k does not overlap with the 5' interval and j connects both segments */
  if (dat->overlap5[i])
    for (size_t cnt = 0; cnt < vrna_array_size(dat->overlap5[i]); cnt++) {
      r = dat->overlap5[i][cnt];

      if ((j >= r->c) &&
          (j <= r->d) && /* done checking whether (i,j) connects interval pair */
          ((l < r->c) || (k > r->b))) {
        return r->exp_score;
      }
    }

  return 1.0;
}


PRIVATE interact_region_s *
parse_crosslink_data(FILE *fp)
{
  char    *line;
  int     a, b, c, d;
  double  score;

  vrna_array(interact_region_s) regions;

  vrna_array_init_size(regions, 42);

  while ((line = vrna_read_line(fp))) {
    if (sscanf(line, "%d %d %d %d %lf", &a, &b, &c, &d, &score) == 5) {
      interact_region_s r = {
        .a      = a,
        .b      = b,
        .c      = c,
        .d      = d,
        .score  = (int)(score * 100.)
      };

      vrna_array_append(regions, r);
    }

    free(line);
  }

  if (vrna_array_size(regions) == 0) {
    vrna_array_free(regions);
    regions = NULL;
  }

  return regions;
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


PUBLIC vrna_bp_probs_f
vrna_sc_iip_get_outside_cb(vrna_fold_compound_t  *fc,
                          void                  **cb_data,
                          vrna_auxdata_free_f   *cb_data_free)
{
  if ((fc) &&
      (cb_data) &&
      (cb_data_free)) {
    *cb_data      = (void *)get_xl_ml_outside_helpers(fc);
    *cb_data_free = free_xl_ml_outside_helpers;
    return compute_bpp_multibranch_xlink;
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
                       xl_ml_helper_t   *helpers)
{
  FLT_OR_DBL    q;
  unsigned int  tt;

  tt = vrna_get_ptype_md(S[j], S[i], md);

  q = vrna_exp_E_multibranch_stem(tt, S1[j - 1], S1[i + 1], pf_params);

  if (helpers->sc_wrapper_mb.pair)
    q *= helpers->sc_wrapper_mb.pair(i, j, &(helpers->sc_wrapper_mb));

  return q;
}


PRIVATE INLINE FLT_OR_DBL
contrib_enclosed_pair(unsigned int     i,
                      unsigned int     j,
                      vrna_md_t        *md,
                      vrna_exp_param_t *pf_params,
                      short            *S,
                      short            *S1,
                      xl_ml_helper_t   *helpers)
{
  FLT_OR_DBL    q;
  unsigned int  tt;

  tt = vrna_get_ptype_md(S[i], S[j], md);

  q = vrna_exp_E_multibranch_stem(tt, S1[i - 1], S1[j + 1], pf_params);

  if (helpers->sc_wrapper_mb.red_stem)
    q *= helpers->sc_wrapper_mb.red_stem(i, j, i, j, &(helpers->sc_wrapper_mb));

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
  helpers->prm_lmm  = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * (n + 2));
  helpers->prm_lmm2 = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * (n + 2));
  helpers->prm_lmu  = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * (n + 2));
  helpers->prm_lmu2 = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * (n + 2));

  init_sc_mb_exp(fc, &(helpers->sc_wrapper_mb));

  return helpers;
}


PRIVATE void
free_xl_ml_outside_helpers(void *helpers)
{
  if (helpers) {
    xl_ml_helper_t *h = (xl_ml_helper_t *)helpers;
    free_sc_mb_exp(&(h->sc_wrapper_mb));

    free(h->vl);
    free(h->wl);
    free(h->yl);
    free(h->tl);
    free(h->prm_lmm);
    free(h->prm_lmm2);
    free(h->prm_lmu);
    free(h->prm_lmu2);

    free(helpers);
  }
}


PRIVATE void
xl_fill_prms(vrna_fold_compound_t *fc,
             unsigned int         l,
             xl_ml_helper_t       *helpers)
{
  short         *S, *S1, s3;
  unsigned int  i, j, n, u, v, tt, overlap, turn;
  int           *my_iindx, ij, lj1;
  FLT_OR_DBL        qbt1, qbt2, temp, temp2, *probs, *qb, *qm, *scale, *expMLbase, expMLclosing, prm_lmm, prm_lmu;
  vrna_exp_param_t  *pf_params;
  vrna_md_t         *md;
  vrna_hc_t         *hc;
  struct vrna_crosslink_data_s  *dat;
  interact_region_s             *r;


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

  /* reset prms */
  for (i = 1; i <= n; i++) {
    helpers->prm_lmm[i]   = 0.;
    helpers->prm_lmm2[i]  = 0.;
    helpers->prm_lmu[i]   = 0.;
    helpers->prm_lmu2[i]  = 0.;
    helpers->vl[i]        = 0.;
    helpers->wl[i]        = 0.;
  }

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
                 contrib_enclosed_pair(u, v, md, pf_params, S, S1, helpers);

          if (helpers->sc_wrapper_mb.red_stem)
            qbt1 *= helpers->sc_wrapper_mb.red_stem(u, v, u, v, &(helpers->sc_wrapper_mb));

          /* no score, since l < r->b, u = l + 1 <= r->b */
          helpers->wl[v] += qbt1;
        }

        for (u = l + 2; (u <= r->b) && (u + turn < v); u++) {
          temp = qb[my_iindx[u] - v] *
                 contrib_enclosed_pair(u, v, md, pf_params, S, S1, helpers);

          if (hc->eval_mb(l + 1, v, u - 1, u, VRNA_DECOMP_ML_ML_STEM, hc)) {
            qbt1 = qm[my_iindx[l + 1] - u + 1];

            if (helpers->sc_wrapper_mb.decomp_ml)
              qbt1 *= helpers->sc_wrapper_mb.decomp_ml(l + 1, v, u - 1, u, &(helpers->sc_wrapper_mb));

            helpers->wl[v] += temp *
                              qbt1;
          }

          if (hc->eval_mb(l + 1, v, u, v, VRNA_DECOMP_ML_STEM, hc)) {
            qbt1 = expMLbase[u - l - 1];

            if (helpers->sc_wrapper_mb.red_stem)
              qbt1 *= helpers->sc_wrapper_mb.red_stem(l + 1, v, u, v, &(helpers->sc_wrapper_mb));

            helpers->wl[v] += temp *
                              qbt1;
          }
        }

        for (; (u < r->c) && (u + turn < v); u++) {
          temp = qb[my_iindx[u] - v] *
                 contrib_enclosed_pair(u, v, md, pf_params, S, S1, helpers) *
                 r->exp_score;

          if (hc->eval_mb(l + 1, v, u - 1, u, VRNA_DECOMP_ML_ML_STEM, hc)) {
            qbt1 = qm[my_iindx[l + 1] - u + 1];

            if (helpers->sc_wrapper_mb.decomp_ml)
              qbt1 *= helpers->sc_wrapper_mb.decomp_ml(l + 1, v, u - 1, u, &(helpers->sc_wrapper_mb));

            helpers->wl[v] += temp *
                              qbt1;
          }

          if (hc->eval_mb(l + 1, v, u, v, VRNA_DECOMP_ML_STEM, hc)) {
            qbt1 = expMLbase[u - l - 1];

            if (helpers->sc_wrapper_mb.red_stem)
              qbt1 *= helpers->sc_wrapper_mb.red_stem(l + 1, v, u, v, &(helpers->sc_wrapper_mb));

            helpers->wl[v] += temp *
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
          temp = helpers->wl[v];

          if ((v + 1 < j) &&
              (hc->eval_mb(l + 1, j - 1, v, v + 1, VRNA_DECOMP_ML_ML_ML, hc))) {
            qbt1 = qm[my_iindx[v + 1] - j + 1];

            if (helpers->sc_wrapper_mb.decomp_ml)
              qbt1 *= helpers->sc_wrapper_mb.decomp_ml(l + 1, j - 1, v, v + 1, &(helpers->sc_wrapper_mb));

            helpers->vl[j] += qbt1 *
                              temp;
          }

          if (hc->eval_mb(l + 1, j - 1, l + 1, v, VRNA_DECOMP_ML_ML, hc)) {
            qbt1 = expMLbase[j - v - 1];

            if (helpers->sc_wrapper_mb.red_ml)
              qbt1 *= helpers->sc_wrapper_mb.red_ml(l + 1, j - 1, l + 1, v, &(helpers->sc_wrapper_mb));

            helpers->vl[j] += qbt1 *
                              temp;
          }
        }
      }
    }
  }

  /* finally, compute prm_lmm and prm_lmm2 */
  for (i = 1; i + turn + 1 < l; i++) {
    s3                    = S1[i + 1];
    helpers->prm_lmm[i]   = 0.;
    helpers->prm_lmm2[i]  = 0.;

    for (j = l + turn + 3; j <= n; j++) {
      ij = my_iindx[i] - j;

      if (hc->eval_mb(i, j, i + 1, j - 1, VRNA_DECOMP_PAIR_ML, hc)) {
        /* compute closing pair contribution */
        temp = probs[ij] *
               contrib_enclosing_pair(i, j, md, pf_params, S, S1, helpers);

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

            if (helpers->sc_wrapper_mb.decomp_ml)
              qbt1 *= helpers->sc_wrapper_mb.decomp_ml(i + 1, j - 1, l, l + 1, &(helpers->sc_wrapper_mb));

            qbt2 = qbt1; /* for prm_lmm2 */
          }
        } else if ((r->b <= l) || (r->c <= l + 1) || (r->c == j)) {
          if (hc->eval_mb(i + 1, j - 1, l, l + 1, VRNA_DECOMP_ML_ML_ML, hc)) {
            qbt1 = qm[my_iindx[l + 1] - j + 1];

            if (helpers->sc_wrapper_mb.decomp_ml)
              qbt1 *= helpers->sc_wrapper_mb.decomp_ml(i + 1, j - 1, l, l + 1, &(helpers->sc_wrapper_mb));

            qbt2 = qbt1; /* for prm_lmm2 */

            qbt1 *= r->exp_score;
          }
        } else {
          if (hc->eval_mb(i + 1, j - 1, l, l + 1, VRNA_DECOMP_ML_ML_ML, hc)) {
            qbt1 = helpers->vl[j];

            FLT_OR_DBL ulj = 0.;

            if (hc->eval_mb(l + 1, j - 1, l + 1, r->c - 1, VRNA_DECOMP_ML_ML, hc)) {
              temp2 = qm[my_iindx[l + 1] - r->c + 1] *
                      expMLbase[j - r->c] *
                      r->exp_score;

              if (helpers->sc_wrapper_mb.red_ml)
                temp2 *= helpers->sc_wrapper_mb.red_ml(l + 1, j - 1, l + 1, r->c - 1, &(helpers->sc_wrapper_mb));

              ulj += temp2;
            }

            if (hc->eval_mb(l + 1, j - 1, r->c, j - 1, VRNA_DECOMP_ML_ML, hc)) {
              temp2 = qm[my_iindx[r->c] - j + 1] *
                      expMLbase[r->c - l - 1] *
                      r->exp_score;

              if (helpers->sc_wrapper_mb.red_ml)
                temp2 *= helpers->sc_wrapper_mb.red_ml(l + 1, j - 1, r->c, j - 1, &(helpers->sc_wrapper_mb));

              ulj += temp2;
            }

            if (hc->eval_mb(l + 1, j - 1, r->c - 1, r->c, VRNA_DECOMP_ML_ML_ML, hc)) {
              temp2 = qm[my_iindx[l + 1] - r->c + 1] *
                      qm[my_iindx[r->c] - j + 1] *
                      r->exp_score;

              if (helpers->sc_wrapper_mb.decomp_ml)
                temp2 *= helpers->sc_wrapper_mb.decomp_ml(l + 1, j - 1, r->c - 1, r->c, &(helpers->sc_wrapper_mb));

              ulj += temp2;
            }

            
            qbt1 += ulj;

            qbt2 = qbt1; /* for prm_lmm2 */
          }
        }

        helpers->prm_lmm[i]   += qbt1 *
                                 temp;
        helpers->prm_lmm2[i]  += qbt2 *
                                 temp;
      }
    }
  }

  /* and the simpler cases, prm_lmu and prm_lmu2 */
  for (i = 1; i + turn + 1 < l; i++) {
    s3                    = S1[i + 1];
    helpers->prm_lmu[i]   = 0.;
    helpers->prm_lmu2[i]  = 0.;

    for (j = l + 1; j <= n; j++) {
      ij = my_iindx[i] - j;

      if (hc->eval_mb(i, j, i + 1, j - 1, VRNA_DECOMP_PAIR_ML, hc)) {
        /* compute closing pair contribution */
        qbt1 = 0.;
        
        temp = probs[ij] *
               contrib_enclosing_pair(i, j, md, pf_params, S, S1, helpers);

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

        if (hc->eval_mb(i + 1, j - 1, i + 1, l, VRNA_DECOMP_ML_ML, hc)) {
          qbt1 = expMLbase[j - l - 1];

          if (helpers->sc_wrapper_mb.red_ml)
            qbt1 *= helpers->sc_wrapper_mb.red_ml(i + 1, j - 1, i + 1, l, &(helpers->sc_wrapper_mb));

          if (r) {
            helpers->prm_lmu[i] += qbt1 *
                                   temp *
                                   r->exp_score;
            if (r->c > l) {
              helpers->prm_lmu2[i] += qbt1 *
                                      temp *
                                      r->exp_score;
            } else {
              helpers->prm_lmu2[i] += qbt1 *
                                      temp;
            }
          } else {
            helpers->prm_lmu[i] += qbt1 *
                                   temp;
            helpers->prm_lmu2[i] += qbt1 *
                                   temp;
          }
        }
      }
    }

    helpers->prm_lmm[i] *= expMLclosing;
    helpers->prm_lmm2[i] *= expMLclosing;
    helpers->prm_lmu[i] *= expMLclosing;
    helpers->prm_lmu2[i] *= expMLclosing;

    vrna_log_debug("%u: prml[%u] = %g", l, i, helpers->prm_lmu[i] + helpers->prm_lmm[i]);

  }
}


PRIVATE void
xl_fill_y(vrna_fold_compound_t *fc,
          unsigned int         l,
          xl_ml_helper_t       *helpers)
{
  short                         *S, *S1;
  unsigned int                  i, k, u, v, turn, n;
  int                           *my_iindx;
  FLT_OR_DBL                    *ol, *ol2, *rl, temp, qbt1, qbt2, *qm, *qb, *expMLbase, *scale;
  vrna_exp_param_t              *pf_params;
  vrna_md_t                     *md;
  vrna_hc_t                     *hc;
  struct vrna_crosslink_data_s  *dat;
  interact_region_s             *r;



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

  ol  = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * (n + 2));
  ol2 = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * (n + 2));
  rl  = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * (n + 2));

  /* reset Y[lk] */
  for (k = 1; k < l; k++)
    helpers->yl[k] = 0.;

  if (dat->overlap3[l]) {
    /* fill O[l,u] and O2[l,u] */
    for (u = 2; u + turn + 2 < l; u++) {
      /* find interval pair (I,I') with u \in I and l \in I', if any */
      r = NULL;

      for (size_t cnt = 0; cnt < vrna_array_size(dat->overlap3[l]); cnt++)
        if ((dat->overlap3[l][cnt]->a < u) &&
            (dat->overlap3[l][cnt]->b >= u)) {
          r = dat->overlap3[l][cnt];
          break;
        }

      if ((r) &&
          (r->d > l)) {
        for (i = r->a; i < u; i++) {
          qbt1 = helpers->prm_lmu[i];
          qbt2 = helpers->prm_lmu2[i];

          if ((i + 1 < u) &&
              (hc->eval_mb(i + 1, l, u - 1, u, VRNA_DECOMP_ML_ML_ML, hc))) {
            temp = qm[my_iindx[i + 1] - u + 1];

            if (helpers->sc_wrapper_mb.decomp_ml)
              temp *= helpers->sc_wrapper_mb.decomp_ml(i + 1, l, u - 1, u, &(helpers->sc_wrapper_mb));

            ol[u] += qbt1 *
                     temp;
            ol2[u] += qbt2 *
                     temp;
          }

          if (hc->eval_mb(i + 1, l, u, l, VRNA_DECOMP_ML_ML, hc)) {
            temp = expMLbase[u - i - 1];

            if (helpers->sc_wrapper_mb.red_ml)
              temp *= helpers->sc_wrapper_mb.red_ml(i + 1, l, u, l, &(helpers->sc_wrapper_mb));

            ol[u] += qbt1 *
                     temp;
            ol2[u] += qbt2 *
                     temp;
          }
        }
      }
    }

    /*
     *  fill R[l,v]
     */
    for (v = turn + 3; v + turn + 1 < l; v++) {
      for (u = 2; u + turn < v; u++) {
        if (hc->eval_mb(u, v, u, v, VRNA_DECOMP_ML_STEM, hc)) {
          /* find out, whether there is any interval pair (I,I') with u \in I and l \in I' */
          r = NULL;

          for (size_t cnt = 0; cnt < vrna_array_size(dat->overlap3[l]); cnt++)
            if ((dat->overlap3[l][cnt]->a < u) &&
                (dat->overlap3[l][cnt]->b >= u)) {
              r = dat->overlap3[l][cnt];
              break;
            }

          if ((r) &&
              (r->d > l) &&
              (r->b < v)) {
            temp = qb[my_iindx[u] - v] *
                   contrib_enclosed_pair(u, v, md, pf_params, S, S1, helpers);

            if (r->c <= v)
              temp *= ol2[u];
            else
              temp *= ol[u];

            if (helpers->sc_wrapper_mb.red_stem)
              temp *= helpers->sc_wrapper_mb.red_stem(u, v, u, v, &(helpers->sc_wrapper_mb));

            rl[v] += temp;
          }
        }
      }
    }

    /* Finally, compute Y[l,k] */
    if ((l < n) &&
        (dat->overlap3[l + 1])) {
      for (k = turn + 4; k + turn < l; k++) {
        /* find out whether there is any interval I' with k \in I' and l \in I' */
        r = NULL;

        for (size_t cnt = 0; cnt < vrna_array_size(dat->overlap3[l + 1]); cnt++)
          if ((dat->overlap3[l + 1][cnt]->c < k) &&
              (dat->overlap3[l + 1][cnt]->d > k)) {
            r = dat->overlap3[l + 1][cnt];
            break;
          }

        if (r == NULL)
          continue;

        for (v = turn + 3; v < k; v++) {
          temp = rl[v];

          if (temp != 0.) {
            if ((v + 1 < k) &&
                (hc->eval_mb(v + 1, l, k - 1, k, VRNA_DECOMP_ML_ML_ML, hc))) {
              qbt1 = qm[my_iindx[v + 1] - k + 1];

              if (helpers->sc_wrapper_mb.decomp_ml)
                qbt1 *= helpers->sc_wrapper_mb.decomp_ml(v + 1, l, k - 1, k, &(helpers->sc_wrapper_mb));

              helpers->yl[k] += qbt1 *
                                temp;
            }

            if (hc->eval_mb(v + 1, l, k, l, VRNA_DECOMP_ML_ML, hc)) {
              qbt1 = expMLbase[k - v - 1];

              if (helpers->sc_wrapper_mb.red_ml)
                qbt1 *= helpers->sc_wrapper_mb.red_ml(v + 1, l, k, l, &(helpers->sc_wrapper_mb));

              helpers->yl[k] += qbt1 *
                                temp;
            }
          }
        }
      }
    }
  }

  free(ol);
  free(ol2);
  free(rl);
}


PRIVATE void
xl_fill_t(vrna_fold_compound_t *fc,
          unsigned int         l,
          xl_ml_helper_t       *helpers)
{
  short                         *S, *S1;
  unsigned int                  i, k, u, v, turn, n;
  int                           *my_iindx;
  FLT_OR_DBL                    *pl, *pl2, *sl, temp, qbt1, qbt2, *qm, *qb, *expMLbase, *scale;
  vrna_exp_param_t              *pf_params;
  vrna_md_t                     *md;
  vrna_hc_t                     *hc;
  struct vrna_crosslink_data_s  *dat;
  interact_region_s             *r;



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

  pl  = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * (n + 2));
  pl2 = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * (n + 2));
  sl  = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * (n + 2));

  /* reset Y[lk] */
  for (k = 1; k < l; k++)
    helpers->tl[k] = 0.;

  if (dat->overlap3[l]) {
    /* fill O[l,u] and O2[l,u] */
    for (u = 2; u + turn + 2 < l; u++) {
      /* find interval pair (I,I') with u \in I and l \in I', if any */
      r = NULL;

      for (size_t cnt = 0; cnt < vrna_array_size(dat->overlap3[l]); cnt++)
        if ((dat->overlap3[l][cnt]->a < u) &&
            (dat->overlap3[l][cnt]->b >= u)) {
          r = dat->overlap3[l][cnt];
          break;
        }

      if ((r) &&
          (r->d > l)) {
        for (i = r->a; i < u; i++) {
          qbt1 = helpers->prm_lmm[i];
          qbt2 = helpers->prm_lmm2[i];

          if ((i + 1 < u) &&
              (hc->eval_mb(i + 1, l, u - 1, u, VRNA_DECOMP_ML_ML_ML, hc))) {
            temp = qm[my_iindx[i + 1] - u + 1];

            if (helpers->sc_wrapper_mb.decomp_ml)
              temp *= helpers->sc_wrapper_mb.decomp_ml(i + 1, l, u - 1, u, &(helpers->sc_wrapper_mb));

            pl[u] += qbt1 *
                     temp;
            pl2[u] += qbt2 *
                     temp;
          }

          if (hc->eval_mb(i + 1, l, u, l, VRNA_DECOMP_ML_ML, hc)) {
            temp = expMLbase[u - i - 1];

            if (helpers->sc_wrapper_mb.red_ml)
              temp *= helpers->sc_wrapper_mb.red_ml(i + 1, l, u, l, &(helpers->sc_wrapper_mb));

            pl[u] += qbt1 *
                     temp;
            pl2[u] += qbt2 *
                     temp;
          }
        }
      }
    }

    /*
     *  fill S[l,v]
     */
    for (v = turn + 3; v + turn + 1 < l; v++) {
      for (u = 2; u + turn < v; u++) {
        if (hc->eval_mb(u, v, u, v, VRNA_DECOMP_ML_STEM, hc)) {
          /* find out, whether there is any interval pair (I,I') with u \in I and l \in I' */
          r = NULL;

          for (size_t cnt = 0; cnt < vrna_array_size(dat->overlap3[l]); cnt++)
            if ((dat->overlap3[l][cnt]->a < u) &&
                (dat->overlap3[l][cnt]->b >= u)) {
              r = dat->overlap3[l][cnt];
              break;
            }

          if ((r) &&
              (r->d > l) &&
              (r->b < v)) {
            temp = qb[my_iindx[u] - v] *
                   contrib_enclosed_pair(u, v, md, pf_params, S, S1, helpers);

            if (r->c <= v)
              temp *= pl2[u];
            else
              temp *= pl[u];

            if (helpers->sc_wrapper_mb.red_stem)
              temp *= helpers->sc_wrapper_mb.red_stem(u, v, u, v, &(helpers->sc_wrapper_mb));

            sl[v] += temp;
          }
        }
      }
    }

    /* Finally, compute T[l,k] */
    if ((l < n) &&
        (dat->overlap3[l + 1])) {
      for (k = turn + 4; k + turn < l; k++) {
        /* find out whether there is any interval I' with k \in I' and l \in I' */
        r = NULL;

        for (size_t cnt = 0; cnt < vrna_array_size(dat->overlap3[l + 1]); cnt++)
          if ((dat->overlap3[l + 1][cnt]->c < k) &&
              (dat->overlap3[l + 1][cnt]->d > k)) {
            r = dat->overlap3[l + 1][cnt];
            break;
          }

        if (r == NULL)
          continue;

        for (v = turn + 3; v < k; v++) {
          temp = sl[v];

          if (temp != 0.) {
            if ((v + 1 < k) &&
                (hc->eval_mb(v + 1, l, k - 1, k, VRNA_DECOMP_ML_ML_ML, hc))) {
              qbt1 = qm[my_iindx[v + 1] - k + 1];

              if (helpers->sc_wrapper_mb.decomp_ml)
                qbt1 *= helpers->sc_wrapper_mb.decomp_ml(v + 1, l, k - 1, k, &(helpers->sc_wrapper_mb));

              helpers->yl[k] += qbt1 *
                                temp;
            }

            if (hc->eval_mb(v + 1, l, k, l, VRNA_DECOMP_ML_ML, hc)) {
              qbt1 = expMLbase[k - v - 1];

              if (helpers->sc_wrapper_mb.red_ml)
                qbt1 *= helpers->sc_wrapper_mb.red_ml(v + 1, l, k, l, &(helpers->sc_wrapper_mb));

              helpers->tl[k] += qbt1 *
                                temp;
            }
          }
        }
      }
    }
  }

  free(pl);
  free(pl2);
  free(sl);
}


PRIVATE void
compute_bpp_multibranch_xlink(vrna_fold_compound_t  *fc,
                              unsigned int          l,
                              FLT_OR_DBL            *Qmax,
                              size_t                *ov,
                              void                  *helpers)
{
  short                         *S, *S1;
  unsigned int                  i, k, n, turn;
  int                           *my_iindx;
  FLT_OR_DBL                    pklm, qbt1, temp, *qm, *expMLbase, *probs, *scale;
  vrna_exp_param_t              *pf_params;
  vrna_md_t                     *md;
  vrna_hc_t                     *hc;
  struct vrna_crosslink_data_s  *dat;
  interact_region_s             *r;
  xl_ml_helper_t                *h;

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
  h         = (xl_ml_helper_t *)helpers;
  dat       = fc->sc->iip;

  xl_fill_prms(fc, l, h);
  xl_fill_y(fc, l, h);
  xl_fill_t(fc, l, h);

  for (k = 2; k + turn < l; k++) {
    if (!hc->eval_mb(k, l, k, l, VRNA_DECOMP_ML_STEM, hc))
      continue;

    pklm = 0.;

    /* 1. (k,l) is left-most stem in multibranch loop */
    for (i = 1; i < k; i++) {
      temp = h->prm_lmm[i];

      if (dat->overlap3[l]) {
        for (size_t cnt = 0; cnt < vrna_array_size(dat->overlap3[l]); cnt++)
          if ((dat->overlap3[l][cnt]->a <= i) &&
              (dat->overlap3[l][cnt]->b > i) &&
              (dat->overlap3[l][cnt]->b >= k))
            temp = h->prm_lmm2[i];
      }

      if (hc->eval_mb(i + 1, l, k, l, VRNA_DECOMP_ML_STEM, hc)) {
        pklm += expMLbase[k - i - 1] *
                temp;
      }
    }

    /* 2./3. (k,l) is right-most stem or somewhere in the middle of a multibranch loop */
    pklm += h->yl[k];
    pklm += h->tl[k];

    for (i = 1; i + turn + 2 < k; i++) {
      if (hc->eval_mb(i + 1, l, k - 1, k, VRNA_DECOMP_ML_ML_ML, hc)) {
        /* find out whether exists an interval pair (I,I') with i \in I and l \in I' */
        if (dat->overlap3[l]) {
          r = NULL;
          for (size_t cnt = 0; cnt < vrna_array_size(dat->overlap3[l]); cnt++)
            if ((dat->overlap3[l][cnt]->a <= i) &&
                (dat->overlap3[l][cnt]->b >= i)) {
              r = dat->overlap3[l][cnt];
              break;
            }

          if (r) {
            if ((i < r->b) &&
                (l + 1 <= r->d) &&
                (r->c + 1 <= k)) {
              temp = 0.;
              
              if (hc->eval_mb(i + 1, k - 1, i + 1, r->b, VRNA_DECOMP_ML_ML, hc)) {
                qbt1 = qm[my_iindx[i + 1] - r->b] *
                       expMLbase[k - r->b - 1];

                if (h->sc_wrapper_mb.red_ml)
                  qbt1 *= h->sc_wrapper_mb.red_ml(i + 1, k - 1, i + 1, r->b, &(h->sc_wrapper_mb));

                temp += qbt1;
              }

              if (r->b < k) {
                if (hc->eval_mb(i + 1, k - 1, r->b + 1, k - 1, VRNA_DECOMP_ML_ML, hc)) {
                  qbt1 = qm[my_iindx[r->b + 1] - k + 1] *
                         expMLbase[r->b - i];

                  if (h->sc_wrapper_mb.red_ml)
                    qbt1 *= h->sc_wrapper_mb.red_ml(i + 1, k - 1, r->b + 1, k - 1, &(h->sc_wrapper_mb));

                  temp += qbt1;
                }

                if (hc->eval_mb(i + 1, k - 1, r->b, r->b + 1, VRNA_DECOMP_ML_ML_ML, hc)) {
                  qbt1 = qm[my_iindx[i + 1] - r->b] *
                         qm[my_iindx[r->b + 1] - k + 1];

                  if (h->sc_wrapper_mb.decomp_ml)
                    qbt1 *= h->sc_wrapper_mb.decomp_ml(i + 1, k - 1, r->b, r->b + 1, &(h->sc_wrapper_mb));
  
                  temp += qbt1;
                }
              }

              pklm += temp *
                      (h->prm_lmu[i] + h->prm_lmm[i]);

              continue;
            } else if ((i == r->b) ||
                       (r->c >= k)) {
              temp = (h->prm_lmu[i] + h->prm_lmm[i]) *
                     qm[my_iindx[i + 1] - k + 1];

              if (h->sc_wrapper_mb.decomp_ml)
                temp *= h->sc_wrapper_mb.decomp_ml(i + 1, l, k - 1, k, &(h->sc_wrapper_mb));

              pklm += temp;

              continue;
            }
          }
        }

        temp = (h->prm_lmu2[i] + h->prm_lmm2[i]) *
               qm[my_iindx[i + 1] - k + 1];

        if (h->sc_wrapper_mb.decomp_ml)
          temp *= h->sc_wrapper_mb.decomp_ml(i + 1, l, k - 1, k, &(h->sc_wrapper_mb));

        pklm += temp;
      }
    }

    pklm *= contrib_enclosed_pair(k, l, md, pf_params, S, S1, helpers) *
            scale[2];

    probs[my_iindx[k] - l] += pklm;


    vrna_log_debug("p_mb[%u,%u] = %g", k, l, pklm);

  }
}

#if 0

PRIVATE void
compute_bpp_multibranch_xlink(vrna_fold_compound_t  *fc,
                              unsigned int          l,
                              FLT_OR_DBL            *Qmax,
                              size_t                *ov,
                              void                  *helpers)
{
  short                         *S1, *S;
  unsigned int                  i, k, n, overlap_5, overlap_3, turn, tt, with_gquad;
  int                           *my_iindx, kl;
  FLT_OR_DBL                    qbt1, temp, expMLstemGQ, *qm, *qb, *expMLbase, *prm_lmu, *prm_lmu2, *prm_lmm, *prm_lmm2, *r_lmu, *r_lmm, *s_lmm, *scale, *probs;
  vrna_exp_param_t              *pf_params;
  vrna_md_t                     *md;
  vrna_hc_t         *hc;
  vrna_smx_csr(FLT_OR_DBL)  *q_gq;
  struct vrna_crosslink_data_s  *dat;
  interact_region_s             *r;

  xl_ml_helper_t *h = (xl_ml_helper_t *)helpers;

  xl_fill_prms(fc, l, h);
  xl_fill_rs(fc, l, h);
  xl_fill_s(fc, l, h);

  dat = fc->sc->iip;

  r = dat->overlap[l];

  overlap_5 = 0;
  overlap_3 = 0;

  if ((r) &&
      (r->c <= l) &&
      (l <= r->d))
    overlap_3 = 1;

  if ((r) &&
      (r->a <= l) &&
      (l <= r->b))
    overlap_5 = 1;

  n         = fc->length;
  S             = fc->sequence_encoding2;
  S1            = fc->sequence_encoding;
  pf_params = fc->exp_params;
  md        = &(pf_params->model_details);
  turn      = md->min_loop_size;
  my_iindx      = fc->iindx;
  qb            = fc->exp_matrices->qb;
  qm            = fc->exp_matrices->qm;
  expMLbase     = fc->exp_matrices->expMLbase;
  scale         = fc->exp_matrices->scale;
  probs         = fc->exp_matrices->probs;
  q_gq          = fc->exp_matrices->q_gq;
  prm_lmu       = h->prm_lmu;
  prm_lmu2      = h->prm_lmu2;
  prm_lmm       = h->prm_lmm;
  prm_lmm2      = h->prm_lmm2;
  r_lmu         = h->r_lmu;
  r_lmm         = h->r_lmm;
  s_lmm         = h->s_lmm;
  hc            = fc->hc;
  with_gquad    = md->gquad ? 1 : 0;
  expMLstemGQ   = (with_gquad) ? vrna_exp_E_multibranch_stem(0, -1, -1, pf_params) : 0;

  if ((r == NULL) ||
      (r->a == l) ||
      (r->b == l)) {
    /* all enclosing pairs (i,j) must receive crosslink-derived contribution if any */
    for (k = 2; k + turn < l; k++) {
      kl = my_iindx[k] - l;

      if (hc->eval_mb(k, l, k, l, VRNA_DECOMP_ML_STEM, hc)) {
        temp = 0.;

        /* contributions from (k,l) being the left-most stem */
        for (i = 1; i < k; i++) {
          if (hc->eval_mb(i + 1, l, k, l, VRNA_DECOMP_ML_ML, hc)) {
            qbt1 = expMLbase[k - i - 1] *
                   prm_lmm2[i];

            if (h->sc_wrapper_mb.red_ml)
              qbt1 *= h->sc_wrapper_mb.red_ml(i + 1, l, k, l, &(h->sc_wrapper_mb));

            temp += qbt1;
          }
        }

        /* contributions from (k,l) being the right-most stem or somewhere in the middle */
        for (i = 1; i + 1 < k; i++) {
          if (hc->eval_mb(i + 1, l, k - 1, k, VRNA_DECOMP_ML_ML_ML, hc)) {
            qbt1 = qm[my_iindx[i + 1] - (k - 1)] *
                   (prm_lmu2[i] + prm_lmm2[i]);

            if (h->sc_wrapper_mb.decomp_ml)
              qbt1 *= h->sc_wrapper_mb.decomp_ml(i + 1, l, k - 1, k, &(h->sc_wrapper_mb));

            temp += qbt1;
          }
        }

        /* finally, add contribution of the pair (k,l) itself and store in probs */
        temp  *= contrib_enclosed_pair(k, l, md, pf_params, S, S1, helpers);

        probs[kl] += temp *
                     scale[2];

      }
      /* what about (k,l) being a G-quad here? */
    }
  } else {
    if ((overlap_5) &&
        (r->a < l) &&
        (r->b > l)) {
      /* add remaining contributions for (k,l) being right-most stem */
      for (k = 2; k + turn < l; k++) {
        kl = my_iindx[k] - l;

        if (hc->eval_mb(k, l, k, l, VRNA_DECOMP_ML_STEM, hc)) {
          temp = 0.;
      
          for (i = 1; i + 1 < k; i++) {
            if (hc->eval_mb(i + 1, l, k - 1, k, VRNA_DECOMP_ML_ML_ML, hc)) {
              qbt1 = qm[my_iindx[i + 1] - (k - 1)] *
                     prm_lmu2[i];

              if (h->sc_wrapper_mb.decomp_ml)
                qbt1 *= h->sc_wrapper_mb.decomp_ml(i + 1, l, k - 1, k, &(h->sc_wrapper_mb));

              temp += qbt1;
            }
          }

          /* finally, add contribution of the pair (k,l) itself and store in probs */
          temp  *= contrib_enclosed_pair(k, l, md, pf_params, S, S1, helpers);

          probs[kl] += temp *
                       scale[2];
        }

        /* G-quad cases require some attention here */
      }
    }

    if (overlap_3) { /* l is located within an interval I' */
      /* cases where k is located before interval I, so no other base pair can connected the same interval pair as the enclosing pair (i,j) */
      for (k = 2; k < r->a; k++) {
        kl = my_iindx[k] - l;
        tt    = vrna_get_ptype_md(S[k], S[l], md);
        temp  = 0.;

        /* skip this if (k,l) doesn't close a base pair or delimits a gquad */
        if (with_gquad) {
          if ((!tt) &&
#ifndef VRNA_DISABLE_C11_FEATURES
            (vrna_smx_csr_get(q_gq, k, l, 0.) == 0.))
#else
            (vrna_smx_csr_FLT_OR_DBL_get(q_gq, k, l, 0.) == 0.))
#endif
          continue;
        } else if ((qb[my_iindx[k] - l] == 0.) ||
                   (!hc->eval_mb(k, l, k, l, VRNA_DECOMP_ML_STEM, hc))) {
          continue;
        }

        /* contributions from (k,l) being the left-most stem */
        for (i = 1; i < k; i++) {
          if (hc->eval_mb(i + 1, l, k, l, VRNA_DECOMP_ML_ML, hc)) {
            qbt1 = expMLbase[k - i - 1] *
                   prm_lmm2[i];

            if (h->sc_wrapper_mb.red_ml)
              qbt1 *= h->sc_wrapper_mb.red_ml(i + 1, l, k, l, &(h->sc_wrapper_mb));

            temp += qbt1;
          }
        }

        /* (k,l) is the right-most stem or somewhere in the middle */
        for (i = 1; i + 1 < k; i++) {
          if (hc->eval_mb(i + 1, l, k - 1, k, VRNA_DECOMP_ML_ML_ML, hc)) {
            qbt1 = qm[my_iindx[i + 1] - (k - 1)] *
                   (prm_lmu2[i] + prm_lmm2[i]);

            if (h->sc_wrapper_mb.decomp_ml)
              qbt1 *= h->sc_wrapper_mb.decomp_ml(i + 1, l, k - 1, k, &(h->sc_wrapper_mb));

            temp += qbt1;
          }
        }

        if ((with_gquad) &&
            (qb[my_iindx[k] - l] == 0.)) {
          temp *= expMLstemGQ;
        } else {
          temp  *= contrib_enclosed_pair(k, l, md, pf_params, S, S1, helpers);
        }

        probs[kl] += temp *
                     scale[2];
      }


      /* cases where k is located within interval I */
      for (k = r->a; k <= r->b; k++) {
        kl = my_iindx[k] - l;
        tt    = vrna_get_ptype_md(S[k], S[l], md);
        temp  = 0.;
      
        /* skip this if (k,l) doesn't close a base pair or delimits a gquad */
        if (with_gquad) {
          if ((!tt) &&
#ifndef VRNA_DISABLE_C11_FEATURES
            (vrna_smx_csr_get(q_gq, k, l, 0.) == 0.))
#else
            (vrna_smx_csr_FLT_OR_DBL_get(q_gq, k, l, 0.) == 0.))
#endif
          continue;
        } else if ((qb[my_iindx[k] - l] == 0.) ||
                   (!hc->eval_mb(k, l, k, l, VRNA_DECOMP_ML_STEM, hc))) {
          continue;
        }

        for (i = 1; i < r->a; i++) {
          /* (k,l) is left-most stem */
          if (hc->eval_mb(i + 1, l, k, l, VRNA_DECOMP_ML_ML, hc)) {
            qbt1 = expMLbase[k - i - 1] *
                   prm_lmm2[i];

            if (h->sc_wrapper_mb.red_ml)
              qbt1 *= h->sc_wrapper_mb.red_ml(i + 1, l, k, l, &(h->sc_wrapper_mb));

            temp += qbt1;
          }

          /* (k,l) is right-most stem or in the middle */
          if ((i + 1 < k) &&
              (hc->eval_mb(i + 1, l, k - 1, k, VRNA_DECOMP_ML_ML_ML, hc))) {
            qbt1 = qm[my_iindx[i + 1] - (k - 1)] *
                   (prm_lmu2[i] + prm_lmm2[i]);

            if (h->sc_wrapper_mb.decomp_ml)
              qbt1 *= h->sc_wrapper_mb.decomp_ml(i + 1, l, k - 1, k, &(h->sc_wrapper_mb));

            temp += qbt1;
          }
        }

        for (i = r->a; i < k; i++) {
          /* (k,l) is left-most stem */
          if (hc->eval_mb(i + 1, l, k, l, VRNA_DECOMP_ML_ML, hc)) {
            qbt1 = expMLbase[k - i - 1] *
                   prm_lmm[i];

            if (h->sc_wrapper_mb.red_ml)
              qbt1 *= h->sc_wrapper_mb.red_ml(i + 1, l, k, l, &(h->sc_wrapper_mb));

            temp += qbt1;
          }

          /* (k,l) is right-most stem or in the middle */
          if ((i + 1 < k) &&
              (hc->eval_mb(i + 1, l, k - 1, k, VRNA_DECOMP_ML_ML_ML, hc))) {
            qbt1 = qm[my_iindx[i + 1] - (k - 1)] *
                   (prm_lmu[i] + prm_lmm[i]);

            if (h->sc_wrapper_mb.decomp_ml)
              qbt1 *= h->sc_wrapper_mb.decomp_ml(i + 1, l, k - 1, k, &(h->sc_wrapper_mb));

            temp += qbt1;
          }
        }

        if ((with_gquad) &&
            (qb[my_iindx[k] - l] == 0.)) {
          temp *= expMLstemGQ;
        } else {
          temp  *= contrib_enclosed_pair(k, l, md, pf_params, S, S1, helpers);
        }

        probs[kl] += temp *
                     scale[2];
      }

      /* cases where k is located in between interval I and I' */
      for (k = r->b + 1; (k < l) && (k <= r->c); k++) {
        kl = my_iindx[k] - l;
        tt    = vrna_get_ptype_md(S[k], S[l], md);
        temp  = 0.;
      
        /* skip this if (k,l) doesn't close a base pair or delimits a gquad */
        if (with_gquad) {
          if ((!tt) &&
#ifndef VRNA_DISABLE_C11_FEATURES
            (vrna_smx_csr_get(q_gq, k, l, 0.) == 0.))
#else
            (vrna_smx_csr_FLT_OR_DBL_get(q_gq, k, l, 0.) == 0.))
#endif
          continue;
        } else if ((qb[my_iindx[k] - l] == 0.) ||
                   (!hc->eval_mb(k, l, k, l, VRNA_DECOMP_ML_STEM, hc))) {
          continue;
        }

        /* contributions from (k,l) being the left-most stem */
        for (i = 1; i < k; i++) {
          if (hc->eval_mb(i + 1, l, k, l, VRNA_DECOMP_ML_ML, hc)) {
            qbt1 = expMLbase[k - i - 1] *
                   prm_lmm2[i];

            if (h->sc_wrapper_mb.red_ml)
              qbt1 *= h->sc_wrapper_mb.red_ml(i + 1, l, k, l, &(h->sc_wrapper_mb));

            temp += qbt1;
          }
        }

        /* (k,l) is the right-most stem or somewhere in the middle */
        for (i = 1; i + 1 < k; i++) {
          if (hc->eval_mb(i + 1, l, k - 1, k, VRNA_DECOMP_ML_ML_ML, hc)) {
            qbt1 = qm[my_iindx[i + 1] - (k - 1)] *
                   (prm_lmu2[i] + prm_lmm2[i]);

            if (h->sc_wrapper_mb.decomp_ml)
              qbt1 *= h->sc_wrapper_mb.decomp_ml(i + 1, l, k - 1, k, &(h->sc_wrapper_mb));

            temp += qbt1;
          }
        }

        if ((with_gquad) &&
            (qb[my_iindx[k] - l] == 0.)) {
          temp *= expMLstemGQ;
        } else {
          temp  *= contrib_enclosed_pair(k, l, md, pf_params, S, S1, helpers);
        }

        probs[kl] += temp *
                     scale[2];
      }

      /* cases where k is located on interval I' */
      for (k = r->c + 1; k < l; k++) {
        kl = my_iindx[k] - l;
        tt    = vrna_get_ptype_md(S[k], S[l], md);
        temp  = 0.;
      
        /* skip this if (k,l) doesn't close a base pair or delimits a gquad */
        if (with_gquad) {
          if ((!tt) &&
#ifndef VRNA_DISABLE_C11_FEATURES
            (vrna_smx_csr_get(q_gq, k, l, 0.) == 0.))
#else
            (vrna_smx_csr_FLT_OR_DBL_get(q_gq, k, l, 0.) == 0.))
#endif
          continue;
        } else if ((qb[my_iindx[k] - l] == 0.) ||
                   (!hc->eval_mb(k, l, k, l, VRNA_DECOMP_ML_STEM, hc))) {
          continue;
        }

        for (i = 1; i < r->a; i++) {
          /* contributions from (k,l) being the left-most stem */
          if (hc->eval_mb(i + 1, l, k, l, VRNA_DECOMP_ML_ML, hc)) {
            qbt1 = expMLbase[k - i - 1] *
                   prm_lmm2[i];

            if (h->sc_wrapper_mb.red_ml)
              qbt1 *= h->sc_wrapper_mb.red_ml(i + 1, l, k, l, &(h->sc_wrapper_mb));

            temp += qbt1;
          }

          /* (k,l) is the right-most stem or somewhere in the middle */
          if (hc->eval_mb(i + 1, l, k - 1, k, VRNA_DECOMP_ML_ML_ML, hc)) {
            qbt1 = qm[my_iindx[i + 1] - (k - 1)] *
                   (prm_lmu2[i] + prm_lmm2[i]);

            if (h->sc_wrapper_mb.decomp_ml)
              qbt1 *= h->sc_wrapper_mb.decomp_ml(i + 1, l, k - 1, k, &(h->sc_wrapper_mb));
  
            temp += qbt1;
          }
        }

        /* special cases for (k,l) being right-most stem */
        temp += r_lmu[k];

        /* special cases for (k,l) being somewhere in the middle */
        temp += r_lmm[k];

        /* remaining cases for (k,l) being left-most stem */
        for (i = r->a; i < k; i++) {
          if (hc->eval_mb(i + 1, l, k, l, VRNA_DECOMP_ML_ML, hc)) {
            qbt1 = expMLbase[k - i - 1] *
                   prm_lmm2[i];

            if (h->sc_wrapper_mb.red_ml)
              qbt1 *= h->sc_wrapper_mb.red_ml(i + 1, l, k, l, &(h->sc_wrapper_mb));

            temp += qbt1;
          }
        }

        if ((with_gquad) &&
            (qb[my_iindx[k] - l] == 0.)) {
          temp *= expMLstemGQ;
        } else {
          temp  *= contrib_enclosed_pair(k, l, md, pf_params, S, S1, helpers);
        }

        probs[kl] += temp *
                     scale[2];
      }
    } else if (overlap_5) {
      /*
       * here all the cases where l overlaps with some I
       * In particular, a < l < b, since l == a and l == b are
       * already covered above...
       * The following cases are only required for (k,l) being the
       * left-most stem, or somewhere in the middle
       */
      for (k = 2; k <= r->a; k++) {
        kl = my_iindx[k] - l;
        tt    = vrna_get_ptype_md(S[k], S[l], md);
        temp  = 0.;

        /* skip this if (k,l) doesn't close a base pair or delimits a gquad */
        if (with_gquad) {
          if ((!tt) &&
#ifndef VRNA_DISABLE_C11_FEATURES
            (vrna_smx_csr_get(q_gq, k, l, 0.) == 0.))
#else
            (vrna_smx_csr_FLT_OR_DBL_get(q_gq, k, l, 0.) == 0.))
#endif
          continue;
        } else if ((qb[my_iindx[k] - l] == 0.) ||
                   (!hc->eval_mb(k, l, k, l, VRNA_DECOMP_ML_STEM, hc))) {
          continue;
        }

        for (i = 1; i < k; i++) {
          /* contributions from (k,l) being the left-most stem */
          if (hc->eval_mb(i + 1, l, k, l, VRNA_DECOMP_ML_ML, hc)) {
            qbt1 = expMLbase[k - i - 1] *
                   prm_lmm2[i];

            if (h->sc_wrapper_mb.red_ml)
              qbt1 *= h->sc_wrapper_mb.red_ml(i + 1, l, k, l, &(h->sc_wrapper_mb));

            temp += qbt1;
          }

          /* (k,l) is somewhere in the middle */
          if ((i + 1 < k) &&
              (hc->eval_mb(i + 1, l, k - 1, k, VRNA_DECOMP_ML_ML_ML, hc))) {
            qbt1 = qm[my_iindx[i + 1] - (k - 1)] *
                   prm_lmm2[i];

            if (h->sc_wrapper_mb.decomp_ml)
              qbt1 *= h->sc_wrapper_mb.decomp_ml(i + 1, l, k - 1, k, &(h->sc_wrapper_mb));

            temp += qbt1;
          }
        }

        if ((with_gquad) &&
            (qb[my_iindx[k] - l] == 0.)) {
          temp *= expMLstemGQ;
        } else {
          temp  *= contrib_enclosed_pair(k, l, md, pf_params, S, S1, helpers);
        }

        probs[kl] += temp *
                     scale[2];
      }
      
      for (k = r->a + 1; k < l; k++) {
        kl = my_iindx[k] - l;
        tt    = vrna_get_ptype_md(S[k], S[l], md);
        temp  = 0.;

        /* skip this if (k,l) doesn't close a base pair or delimits a gquad */
        if (with_gquad) {
          if ((!tt) &&
#ifndef VRNA_DISABLE_C11_FEATURES
            (vrna_smx_csr_get(q_gq, k, l, 0.) == 0.))
#else
            (vrna_smx_csr_FLT_OR_DBL_get(q_gq, k, l, 0.) == 0.))
#endif
          continue;
        } else if ((qb[my_iindx[k] - l] == 0.) ||
                   (!hc->eval_mb(k, l, k, l, VRNA_DECOMP_ML_STEM, hc))) {
          continue;
        }

        for (i = 1; i < r->a; i++) {
          /* contributions from (k,l) being the left-most stem */
          if (hc->eval_mb(i + 1, l, k, l, VRNA_DECOMP_ML_ML, hc)) {
            qbt1 = expMLbase[k - i - 1] *
                   prm_lmm2[i];

            if (h->sc_wrapper_mb.red_ml)
              qbt1 *= h->sc_wrapper_mb.red_ml(i + 1, l, k, l, &(h->sc_wrapper_mb));

            temp += qbt1;
          }

          /* (k,l) is somewhere in the middle */
          if ((i + 1 < k) &&
              (hc->eval_mb(i + 1, l, k - 1, k, VRNA_DECOMP_ML_ML_ML, hc))) {
            qbt1 = qm[my_iindx[i + 1] - (k - 1)] *
                   prm_lmm2[i];

            if (h->sc_wrapper_mb.decomp_ml)
              qbt1 *= h->sc_wrapper_mb.decomp_ml(i + 1, l, k - 1, k, &(h->sc_wrapper_mb));

            temp += qbt1;
          }
        }

        /* remaining cases for (k,l) being left-most stem or somewhere in the middle */
        for (i = r->a; i < k; i++) {
          if (hc->eval_mb(i + 1, l, k, l, VRNA_DECOMP_ML_ML, hc)) {
            qbt1 = expMLbase[k - i - 1] *
                   s_lmm[i];

            if (h->sc_wrapper_mb.red_ml)
              qbt1 *= h->sc_wrapper_mb.red_ml(i + 1, l, k, l, &(h->sc_wrapper_mb));

            temp += qbt1;
          }

          /* (k,l) is somewhere in the middle */
          if ((i + 1 < k) &&
              (hc->eval_mb(i + 1, l, k - 1, k, VRNA_DECOMP_ML_ML_ML, hc))) {
            qbt1 = qm[my_iindx[i + 1] - (k - 1)] *
                   s_lmm[i];

            if (h->sc_wrapper_mb.decomp_ml)
              qbt1 *= h->sc_wrapper_mb.decomp_ml(i + 1, l, k - 1, k, &(h->sc_wrapper_mb));

            temp += qbt1;
          }
        }

        if ((with_gquad) &&
            (qb[my_iindx[k] - l] == 0.)) {
          temp *= expMLstemGQ;
        } else {
          temp  *= contrib_enclosed_pair(k, l, md, pf_params, S, S1, helpers);
        }

        probs[kl] += temp *
                     scale[2];
      }
    }
  }
}


PRIVATE xl_ml_helper_t *
get_xl_ml_outside_helpers(vrna_fold_compound_t *fc)
{
  unsigned int    n;
  xl_ml_helper_t  *helpers;

  n                 = fc->length;
  helpers           = (xl_ml_helper_t *)vrna_alloc(sizeof(xl_ml_helper_t));
  helpers->prm_lmu  = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * (n + 2));
  helpers->prm_lmu2 = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * (n + 2));
  helpers->prm_lmm  = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * (n + 2));
  helpers->prm_lmm2 = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * (n + 2));
  helpers->r_lmu    = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * (n + 2));
  helpers->r_lmm    = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * (n + 2));
  helpers->s_lmm    = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * (n + 2));

  init_sc_mb_exp(fc, &(helpers->sc_wrapper_mb));

  return helpers;
}


PRIVATE void
free_xl_ml_outside_helpers(void *helpers)
{
  if (helpers) {
    xl_ml_helper_t *h = (xl_ml_helper_t *)helpers;
    free(h->prm_lmu);
    free(h->prm_lmu2);
    free(h->prm_lmm);
    free(h->prm_lmm2);
    free(h->r_lmu);
    free(h->r_lmm);
    free(h->s_lmm);
    free_sc_mb_exp(&(h->sc_wrapper_mb));

    free(helpers);
  }
}


PRIVATE void
xl_fill_prms(vrna_fold_compound_t *fc,
             unsigned int         l,
             xl_ml_helper_t       *helpers)
{
  short         *S, *S1, s3;
  unsigned int  i, j, n, u, tt, overlap, turn;
  int           *my_iindx, ij, lj1;
  FLT_OR_DBL    temp, *probs, *qm, *scale, *expMLbase, expMLclosing, prm_lmm, prm_lmu;
  vrna_exp_param_t  *pf_params;
  vrna_md_t         *md;
  vrna_hc_t         *hc;
  struct vrna_crosslink_data_s  *dat;
  interact_region_s             *r;


  n             = fc->length;
  S             = fc->sequence_encoding2;
  S1            = fc->sequence_encoding;
  my_iindx      = fc->iindx;
  pf_params     = fc->exp_params;
  md            = &(pf_params->model_details);
  turn          = md->min_loop_size;
  qm            = fc->exp_matrices->qm;
  probs         = fc->exp_matrices->probs;
  scale         = fc->exp_matrices->scale;
  expMLbase     = fc->exp_matrices->expMLbase;
  expMLclosing  = pf_params->expMLclosing;
  hc            = fc->hc;
  dat           = fc->sc->iip;

//  vrna_log_error("l=%u", l);

  for (i = 1; i + turn + 1 < l; i++) {
    r   = dat->overlap[i];
    s3  = S1[i + 1];
    ij  = my_iindx[i] - (l + 1);
    lj1 = my_iindx[l + 1] - l;
    helpers->prm_lmu[i] = 0.;
    helpers->prm_lmm[i] = 0.;
    helpers->prm_lmu2[i] = 0.;
    helpers->prm_lmm2[i] = 0.;
    u       = 0;
    overlap = 0;

    if ((r) &&
        (r->a <= i) &&
        (i <= r->b))
      overlap = 1;

    for (j = l + 1; j <= n; j++, u++, ij--, lj1--) {
      if (hc->eval_mb(i, j, i + 1, j - 1, VRNA_DECOMP_PAIR_ML, hc)) {
        prm_lmu = prm_lmm = 0.;

        /* compute closing pair contribution */
        temp = probs[ij] *
               contrib_enclosing_pair(i, j, md, pf_params, S, S1, helpers) *
               expMLclosing;

        /* compute total helper contributions */
        if (hc->eval_mb(i + 1, j - 1, i + 1, l, VRNA_DECOMP_ML_ML, hc))
          prm_lmu = expMLbase[u] *
                    temp;

        if ((l + turn + 1 < j) &&
            (hc->eval_mb(i + 1, j - 1, l, l + 1, VRNA_DECOMP_ML_ML_ML, hc)))
          prm_lmm = qm[my_iindx[l + 1] - (j - 1)] *
                    temp;

        /* add total helper contributions */
        helpers->prm_lmu[i] += prm_lmu;
        helpers->prm_lmm[i] += prm_lmm;
        
        if ((overlap) &&
            (r->c <= j) &&
            (j <= r->d)) {
          prm_lmu *= r->exp_score;
          prm_lmm *= r->exp_score;
        }

        helpers->prm_lmu2[i] += prm_lmu;
        helpers->prm_lmm2[i] += prm_lmm;
      }
    }

    //vrna_log_error("prm_lmu[%u] = %g", i, helpers->prm_lmu[i]);
    //vrna_log_error("prm_lmm[%u] = %g", i, helpers->prm_lmm[i]);
    //vrna_log_error("prm_lm[%u] = %g", i, helpers->prm_lmu[i] + helpers->prm_lmm[i]);
    //vrna_log_error("prm_lmu2[%u] = %g", i, helpers->prm_lmu2[i]);
    //vrna_log_error("prm_lmm2[%u] = %g", i, helpers->prm_lmm2[i]);
    //vrna_log_error("prm_lm2[%u] = %g", i, helpers->prm_lmu2[i] + helpers->prm_lmm2[i]);
  }
}


PRIVATE void
xl_fill_rs(vrna_fold_compound_t  *fc,
           unsigned int          l,
           xl_ml_helper_t        *helpers)
{
  short         *S, *S1, s3;
  unsigned int  n, a, i, k, u, v, tt, turn;
  int           *my_iindx, ij, lj1;
  FLT_OR_DBL    temp, qbt1, *r_lmu, *r_lmm, *x, *x2, *x3, *x4, *y, *y2, *y3, *y4, *prm_lmu, *prm_lmu2, *prm_lmm, *prm_lmm2, *qm, *qb, *expMLbase;
  vrna_exp_param_t  *pf_params;
  vrna_md_t         *md;
  vrna_hc_t         *hc;
  struct vrna_crosslink_data_s  *dat;
  interact_region_s             *r;

  dat = fc->sc->iip;
  r   = dat->overlap[l];

  /* only compute this contribution if l overlaps with some I' */
  if ((r) &&
      (r->c <= l) &&
      (l <= r->d)) {
    n             = fc->length;
    S             = fc->sequence_encoding2;
    S1            = fc->sequence_encoding;
    my_iindx      = fc->iindx;
    pf_params     = fc->exp_params;
    md            = &(pf_params->model_details);
    turn          = md->min_loop_size;
    qb            = fc->exp_matrices->qb;
    qm            = fc->exp_matrices->qm;
    expMLbase     = fc->exp_matrices->expMLbase;
    hc            = fc->hc;

    r_lmu     = helpers->r_lmu;
    r_lmm     = helpers->r_lmm;
    prm_lmu   = helpers->prm_lmu;
    prm_lmm   = helpers->prm_lmm;
    prm_lmu2  = helpers->prm_lmu2;
    prm_lmm2  = helpers->prm_lmm2;

    x   = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * (n + 2));
    x2  = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * (n + 2));
    x3  = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * (n + 2));
    x4  = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * (n + 2));
    y   = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * (n + 2));
    y2  = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * (n + 2));
    y3  = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * (n + 2));
    y4  = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * (n + 2));

    /* fill x's */
    for (u = r->a + 1; u <= r->b; u++) {
      x[u]  = 0.;
      x2[u] = 0.;
      x3[u] = 0.;
      x4[u] = 0.;

      for (i = r->a; i < u; i++) {
        temp = 0.;

        /* [i + 1 : u - 1] harbors at least one more stem */
        if ((i + 1 < u) &&
            (hc->eval_mb(i + 1, l, u - 1, u, VRNA_DECOMP_ML_ML_ML, hc))) {
          qbt1 = qm[my_iindx[i + 1] - (u - 1)];

          if (helpers->sc_wrapper_mb.decomp_ml)
            qbt1 *= helpers->sc_wrapper_mb.decomp_ml(i + 1, l, u - 1, u, &(helpers->sc_wrapper_mb));

          temp += qbt1;
        }

        /* [i + 1 : u - 1] is unpaired */
        if (hc->eval_mb(i + 1, l, u, l, VRNA_DECOMP_ML_ML, hc)) {
          qbt1 = expMLbase[u - i - 1];

          if (helpers->sc_wrapper_mb.red_ml)
            qbt1 *= helpers->sc_wrapper_mb.red_ml(i + 1, l, u, l, &(helpers->sc_wrapper_mb));

          temp += qbt1;
        }

        x[u] += temp *
                prm_lmu[i];
        
        x2[u] += temp *
                 prm_lmu2[i];
        
        x3[u] += temp *
                 prm_lmm[i];
        
        x4[u] += temp *
                 prm_lmm2[i];
        
      }
    }

    /* fill y's */
    for (v = r->b + 1; v < l; v++) {
      y[v] = 0.;
      y2[v] = 0.;
      y3[v] = 0.;
      y4[v] = 0.;

      for (u = r->a + 1; u <= r->b; u++) {
        if (hc->eval_mb(u, v, u, v, VRNA_DECOMP_ML_STEM, hc)) {
          temp = qb[my_iindx[u] - v] *
                 contrib_enclosed_pair(u, v, md, pf_params, S, S1, helpers);

          y[v] += temp *
                  x[u];

          y2[v] += temp *
                   x2[u];

          y3[v] += temp *
                   x3[u];

          y4[v] += temp *
                   x4[u];
        }
      }
    }

    /* finally, fill r */
    for (k = r->c + 1; k < l; k++) {
      helpers->r_lmu[k] = 0.;
      helpers->r_lmm[k] = 0.;

      if (hc->eval_mb(k, l, k, l, VRNA_DECOMP_ML_STEM, hc)) {
        for (i = r->a; i <= r->b; i++) {
          temp = 0.;

          /* first, obvious contributions where no base pair crosses b */
          if (hc->eval_mb(i + 1, k - 1, r->b + 1, k - 1, VRNA_DECOMP_ML_ML, hc)) {
            qbt1 = expMLbase[r->b - i] *
                   qm[my_iindx[r->b + 1] - (k - 1)];

            if (helpers->sc_wrapper_mb.red_ml)
              qbt1 *= helpers->sc_wrapper_mb.red_ml(i + 1, k - 1, r->b + 1, k - 1, &(helpers->sc_wrapper_mb));

            temp += qbt1;
          }

          if ((i < r->b) &&
              (hc->eval_mb(i + 1, k - 1, i + 1, r->b, VRNA_DECOMP_ML_ML, hc))) {
            qbt1 = qm[my_iindx[i + 1] - r->b] *
                   expMLbase[k - r->b - 1];

            if (helpers->sc_wrapper_mb.red_ml)
              qbt1 *= helpers->sc_wrapper_mb.red_ml(i + 1, k - 1, i + 1, r->b, &(helpers->sc_wrapper_mb));

            temp += qbt1;
          }

          if ((i < r->b) &&
              (hc->eval_mb(i + 1, k - 1, r->b, r->b + 1, VRNA_DECOMP_ML_ML_ML, hc))) {
            qbt1 = qm[my_iindx[i + 1] - r->b] *
                   qm[my_iindx[r->b + 1] - (k - 1)];

            if (helpers->sc_wrapper_mb.decomp_ml)
              qbt1 *= helpers->sc_wrapper_mb.decomp_ml(i + 1, k - 1, r->b, r->b + 1, &(helpers->sc_wrapper_mb));

            temp += qbt1;
          }

          helpers->r_lmu[k] += temp *
                               prm_lmu2[i];

          helpers->r_lmm[k] += temp *
                               prm_lmm2[i];
        }

        for (i = r->b + 1; i + turn < k; i++) {
          temp = 0.;

          if (hc->eval_mb(i + 1, l, k - 1, k, VRNA_DECOMP_ML_ML_ML, hc)) {
            qbt1 = qm[my_iindx[i + 1] - k + 1];

            if (helpers->sc_wrapper_mb.decomp_ml)
              qbt1 *= helpers->sc_wrapper_mb.decomp_ml(i + 1, l, k - 1, k, &(helpers->sc_wrapper_mb));

            temp += qbt1;
          }

          helpers->r_lmu[k] += temp *
                               prm_lmu2[i];

          helpers->r_lmm[k] += temp *
                               prm_lmm2[i];
        }

        /* now the other ones where a base pair crosses b */
        for (v = r->b + 1; v < r->c; v++) {
          temp = 0.;

          if (hc->eval_mb(v + 1, l, k, l, VRNA_DECOMP_ML_ML, hc)) {
            qbt1 = expMLbase[k - v - 1];

            if (helpers->sc_wrapper_mb.red_ml)
              qbt1 *= helpers->sc_wrapper_mb.red_ml(v + 1, l, k, l, &(helpers->sc_wrapper_mb));

            temp += qbt1;
          }

          if  (hc->eval_mb(v + 1, l, k - 1, k, VRNA_DECOMP_ML_ML_ML, hc)) {
            qbt1 = qm[my_iindx[v + 1] - (k - 1)];

            if (helpers->sc_wrapper_mb.decomp_ml)
              qbt1 *= helpers->sc_wrapper_mb.decomp_ml(v + 1, l, k - 1, k, &(helpers->sc_wrapper_mb));

            temp += qbt1;
          }

          helpers->r_lmu[k] += temp *
                               y2[v];

          helpers->r_lmm[k] += temp *
                               y4[v];
        }

        for (v = r->c; v < k; v++) {
          temp = 0.;

          if (hc->eval_mb(v + 1, l, k, l, VRNA_DECOMP_ML_ML, hc)) {
            qbt1 = expMLbase[k - v - 1];

            if (helpers->sc_wrapper_mb.red_ml)
              qbt1 *= helpers->sc_wrapper_mb.red_ml(v + 1, l, k, l, &(helpers->sc_wrapper_mb));

            temp += qbt1;
          }

          if ((v + 1 < k) &&
              (hc->eval_mb(v + 1, l, k - 1, k, VRNA_DECOMP_ML_ML_ML, hc))) {
            qbt1 = qm[my_iindx[v + 1] - (k - 1)];

            if (helpers->sc_wrapper_mb.decomp_ml)
              qbt1 *= helpers->sc_wrapper_mb.decomp_ml(v + 1, l, k - 1, k, &(helpers->sc_wrapper_mb));

            temp += qbt1;
          }

          helpers->r_lmu[k] += temp *
                               y[v];

          helpers->r_lmm[k] += temp *
                               y3[v];
        }

        //helpers->r_lmu[k] *= contrib_enclosed_pair(k, l, md, pf_params, S, S1, helpers);
        //helpers->r_lmm[k] *= contrib_enclosed_pair(k, l, md, pf_params, S, S1, helpers);
      }
    }

    free(x);
    free(x2);
    free(x3);
    free(x4);
    free(y);
    free(y2);
    free(y3);
    free(y4);
  }
}


PRIVATE void
xl_fill_s(vrna_fold_compound_t *fc,
          unsigned int         l,
          xl_ml_helper_t       *helpers)
{
  short         *S, *S1, s3;
  unsigned int  i, j, n, u, v, tt, overlap, overlap2;
  int           *my_iindx, ij, lj1;
  FLT_OR_DBL    temp, qbt1, qbt2, *probs, *qm, *qb, *scale, *expMLbase, expMLclosing, *L, *M, *M2, *s_lmm;
  vrna_exp_param_t  *pf_params;
  vrna_md_t         *md;
  vrna_hc_t         *hc;
  struct vrna_crosslink_data_s  *dat;
  interact_region_s             *r, *r2;


  dat = fc->sc->iip;
  r   = dat->overlap[l];

  /* we only compute this part if a < l < b */
  if ((r) &&
      (r->a < l) &&
      (l < r->b)) {
    n             = fc->length;
    S             = fc->sequence_encoding2;
    S1            = fc->sequence_encoding;
    my_iindx      = fc->iindx;
    pf_params     = fc->exp_params;
    md            = &(pf_params->model_details);
    qb            = fc->exp_matrices->qb;
    qm            = fc->exp_matrices->qm;
    probs         = fc->exp_matrices->probs;
    scale         = fc->exp_matrices->scale;
    expMLbase     = fc->exp_matrices->expMLbase;
    expMLclosing  = pf_params->expMLclosing;
    hc            = fc->hc;

    L   = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * (n + 2));
    M   = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * (n + 2));
    M2  = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * (n + 2));
    s_lmm = helpers->s_lmm;

    /* fill L */
    for (v = r->b + 1; v < n; v++) {
      L[v] = 0.;

      for (u = l + 1; u <= r->b; u++) {
        if (hc->eval_mb(u, v, u, v, VRNA_DECOMP_ML_STEM, hc)) {
          /* compute closing pair contribution */
          temp  = qb[my_iindx[u] - v] *
                  contrib_enclosed_pair(u, v, md, pf_params, S, S1, helpers);

          if ((l + 1 < u) &&
              (hc->eval_mb(l + 1, v, u - 1, u, VRNA_DECOMP_ML_ML_ML, hc))) {
            qbt1 = qm[my_iindx[l + 1] - (u - 1)];

            if (helpers->sc_wrapper_mb.decomp_ml)
              qbt1 *= helpers->sc_wrapper_mb.decomp_ml(l + 1, v, u - 1, u, &(helpers->sc_wrapper_mb));

            L[v] += temp *
                    qbt1;
          }

          if (hc->eval_mb(l + 1, v, u, v, VRNA_DECOMP_ML_ML, hc)) {
            qbt1 = expMLbase[u - l - 1];

            if (helpers->sc_wrapper_mb.red_ml)
              qbt1 *= helpers->sc_wrapper_mb.red_ml(l + 1, v, u, v, &(helpers->sc_wrapper_mb));

            L[v] += temp *
                    qbt1;
          }
        }
      }
    }

    /* fill Ms */
    for (j = l + 1; j <= r->d; j++) {
      M[j] = M2[j] = 0.;

      for (v = r->b + 1; v < r->c; v++) {
        if ((v + 1 < j) &&
            (hc->eval_mb(l + 1, j - 1, v, v + 1, VRNA_DECOMP_ML_ML_ML, hc))) {
          qbt1 = qm[my_iindx[v + 1] - (j - 1)];

          if (helpers->sc_wrapper_mb.decomp_ml)
            qbt1 *= helpers->sc_wrapper_mb.decomp_ml(l + 1, j - 1, v, v + 1, &(helpers->sc_wrapper_mb));

          M2[j] += L[v] *
                   qbt1;
        }

        if (hc->eval_mb(l + 1, j - 1, l + 1, v, VRNA_DECOMP_ML_ML, hc)) {
          qbt1 = expMLbase[j - v - 1];

          if (helpers->sc_wrapper_mb.red_ml)
            qbt1 *= helpers->sc_wrapper_mb.red_ml(l + 1, j - 1, l + 1, v, &(helpers->sc_wrapper_mb));

          M2[j] += L[v] *
                   qbt1;          
        }
      }

      for (v = r->c; v < j; v++) {
        if ((v + 1 < j) &&
            (hc->eval_mb(l + 1, j - 1, v, v + 1, VRNA_DECOMP_ML_ML_ML, hc))) {
          qbt1 = qm[my_iindx[v + 1] - (j - 1)];

          if (helpers->sc_wrapper_mb.decomp_ml)
            qbt1 *= helpers->sc_wrapper_mb.decomp_ml(l + 1, j - 1, v, v + 1, &(helpers->sc_wrapper_mb));

          M[j] += L[v] *
                  qbt1;
        }

        if (hc->eval_mb(l + 1, j - 1, l + 1, v, VRNA_DECOMP_ML_ML, hc)) {
          qbt1 = expMLbase[j - v - 1];

          if (helpers->sc_wrapper_mb.red_ml)
            qbt1 *= helpers->sc_wrapper_mb.red_ml(l + 1, j - 1, l + 1, v, &(helpers->sc_wrapper_mb));

          M[j] += L[v] *
                  qbt1;          
        }
      }
    }

    /* finally, compute S[i] */
    for (i = 1; i < l; i++) {
      s_lmm[i] = 0.;
      overlap  = 0.;

      if (r->a <= i)
        overlap = 1;

      for (j = l + 1; j <= r->d; j++) {
        if (hc->eval_mb(i, j, i + 1, j - 1, VRNA_DECOMP_PAIR_ML, hc)) {
          /* compute closing pair contribution */
          temp = probs[my_iindx[i] - j] *
                 contrib_enclosing_pair(i, j, md, pf_params, S, S1, helpers) *
                 expMLclosing;

          /* first, the obvious cases where no base pair crosses b */
          qbt2 = 0.;

          if (j <= r->b) {
            if ((l + 1 < j) &&
                (hc->eval_mb(i + 1, j - 1, l, l + 1, VRNA_DECOMP_ML_ML_ML, hc))) {
              qbt1 = qm[my_iindx[l + 1] - (j - 1)];

              if (helpers->sc_wrapper_mb.decomp_ml)
                qbt1 *= helpers->sc_wrapper_mb.decomp_ml(i + 1, j - 1, l, l + 1, &(helpers->sc_wrapper_mb));

              qbt2 += qbt1;
            }
          } else {
            if ((r->b + 1 < j) &&
                (hc->eval_mb(l + 1, j - 1, r->b, r->b + 1, VRNA_DECOMP_ML_ML_ML, hc))) {
              qbt1 = qm[my_iindx[r->b + 1] - (j - 1)] *
                     qm[my_iindx[l + 1] - r->b];

              if (helpers->sc_wrapper_mb.decomp_ml)
                qbt1 *= helpers->sc_wrapper_mb.decomp_ml(l + 1, j - 1, r->b, r->b + 1, &(helpers->sc_wrapper_mb));

              qbt2 += qbt1;
            }

            if ((r->b + 1 < j) &&
                (hc->eval_mb(l + 1, j - 1, r->b + 1, j - 1, VRNA_DECOMP_ML_ML, hc))) {
              qbt1 = expMLbase[r->b - l] *
                     qm[my_iindx[r->b + 1] - (j - 1)];

              if (helpers->sc_wrapper_mb.red_ml)
                qbt1 *= helpers->sc_wrapper_mb.red_ml(l + 1, j - 1, r->b + 1, j - 1, &(helpers->sc_wrapper_mb));

              qbt2 += qbt1;
            }

            if ((r->b < j) &&
                (hc->eval_mb(l + 1, j - 1, l + 1, r->b, VRNA_DECOMP_ML_ML, hc))) {
              qbt1 = qm[my_iindx[l + 1] - r->b] *
                     expMLbase[j - r->b - 1];

              if (helpers->sc_wrapper_mb.red_ml)
                qbt1 *= helpers->sc_wrapper_mb.red_ml(l + 1, j - 1, l + 1, r->b, &(helpers->sc_wrapper_mb));

              qbt2 += qbt1;
            }

            if (r->c <= j)
              qbt2 *= r->exp_score;
          }

          /* store contributions from obvious cases */
          s_lmm[i] += temp *
                      qbt2;

          /* next, the cases where there might be a base pair spanning b */
          qbt2 = M[j];

          qbt1 = M2[j];

          if ((overlap) &&
              (r->c <= j))
            qbt1 *= r->exp_score;

          qbt2 += qbt1;

          s_lmm[i] += temp *
                      qbt2;
        }
      }

      /* now the remaining cases where d < j */
      overlap2  = 0;
      r2        = dat->overlap[i];

      if ((r2) &&
          (r2->a <= i) &&
          (i <= r2->b))
        overlap2 = 1;
      
      for (j = r->d + 1; j <= n; j++) {
        if (hc->eval_mb(i, j, i + 1, j - 1, VRNA_DECOMP_PAIR_ML, hc)) {
          /* compute closing pair contribution */
          temp = probs[my_iindx[i] - j] *
                 contrib_enclosing_pair(i, j, md, pf_params, S, S1, helpers) *
                 expMLclosing;

          if ((overlap2) &&
              (r2->c <= j) &&
              (j <= r2->d))
            temp *= r2->exp_score;

          s_lmm[i] += temp *
                      qm[my_iindx[l + 1] - (j - 1)];
        }
      }
    }

    free(L);
    free(M);
    free(M2);
  }
}
#endif
