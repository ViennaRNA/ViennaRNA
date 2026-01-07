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


#include "ViennaRNA/eval/multibranch.h"
#include "ViennaRNA/backtrack/multibranch.h"

#ifdef __GNUC__
# define INLINE inline
#else
# define INLINE
#endif

#include "ViennaRNA/constraints/multibranch_sc.inc"
#include "ViennaRNA/constraints/multibranch_sc_pf.inc"


#include "ViennaRNA/probing/crosslink.h"

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
  unsigned int  a;      /* start of 5' interval */
  unsigned int  b;      /* end of 5' interval */
  unsigned int  c;      /* start of 3' interval */
  unsigned int  d;      /* end of 3' interval */
  int           score;  /* interval score */
  FLT_OR_DBL    exp_score; /* interval score as Boltzmann factor */
} interact_region_s;


struct vrna_crosslink_data_s {
  vrna_array(interact_region_s)   regions;  /* list of interacting regions */
  vrna_array(interact_region_s *) overlap;  /* for each nucleotide position i we store the corresponding interval pair that overlaps with i either with I or I' */
  int *A;
  int *B;
  FLT_OR_DBL  *qA;
  FLT_OR_DBL  *qB;
};



PRIVATE INLINE int
ml_pair_d2(vrna_fold_compound_t         *fc,
           unsigned int                 i,
           unsigned int                 j,
           struct vrna_mx_mfe_aux_ml_s  *helpers,
           struct sc_mb_dat             *sc_wrapper);


PRIVATE INLINE int
ml_pair_d0(vrna_fold_compound_t         *fc,
           unsigned int                 i,
           unsigned int                 j,
           struct vrna_mx_mfe_aux_ml_s  *helpers,
           struct sc_mb_dat             *sc_wrapper);


PRIVATE INLINE int
mfe_M2_xlink(vrna_fold_compound_t         *fc,
             unsigned int                 i,
             unsigned int                 j,
             struct vrna_mx_mfe_aux_ml_s  *helpers,
             struct sc_mb_dat             *sc_wrapper);


PRIVATE void
fill_helper_arrays(vrna_fold_compound_t *fc,
                   unsigned int         i,
                   unsigned int         j,
                   struct sc_mb_dat     *sc_wrapper);


PRIVATE unsigned int
bt_mb_loop_xl(vrna_fold_compound_t *fc,
              unsigned int         i,                     
              unsigned int         j,                     
              int                  en,                   
              vrna_bps_t           bp_stack VRNA_UNUSED,  
              vrna_bts_t           bt_stack);


PRIVATE unsigned int
bt_m2(vrna_fold_compound_t *fc,
      unsigned int         i,                     
      unsigned int         j,                     
      int                  en,                   
      vrna_bps_t           bp_stack VRNA_UNUSED,  
      vrna_bts_t           bt_stack);


PRIVATE unsigned int
bt_A(vrna_fold_compound_t *fc,
     int                  *A,
     unsigned int         i,                     
     unsigned int         j,                     
     vrna_bps_t           bp_stack VRNA_UNUSED,  
     vrna_bts_t           bt_stack);

PRIVATE unsigned int
bt_B(vrna_fold_compound_t *fc,
     int                  *B,
     unsigned int         i,                     
     unsigned int         j,                     
     vrna_bps_t           bp_stack VRNA_UNUSED,  
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
            struct vrna_mx_pf_aux_ml_s  *helpers);


PRIVATE INLINE FLT_OR_DBL
pf_M2_xlink(vrna_fold_compound_t        *fc,
            unsigned int                i,
            unsigned int                j,
            struct vrna_mx_pf_aux_ml_s  *helpers,
            struct sc_mb_exp_dat        *sc_wrapper);


PUBLIC int
vrna_sc_xl_load(vrna_fold_compound_t  *fc,
                const char            *filename)
{
  int                           ret;
  FILE                          *fp;
  struct vrna_crosslink_data_s  *data;
  vrna_array(interact_region_s) regions;

  ret = 0;

  if ((fc) &&
      (filename)) {
    if (!(fp = fopen(filename, "r"))) {
      vrna_log_warning("Can't open crosslink-data file \"%s\"!",
                           filename);
      return 0;
    }

    regions = parse_crosslink_data(fp);

    if (regions) {
      vrna_sc_init(fc);
      fc->sc->xlink = 1;
      ret = 1;

      data = (struct vrna_crosslink_data_s *)vrna_alloc(sizeof(struct vrna_crosslink_data_s));

      data->regions       = regions;
      data->overlap       = NULL;

      fc->sc->xlink_data  = (void *)data;
      fc->params->model_details.uniq_ML = 1;

      ret &= vrna_sc_multi_cb_add(fc, xlink_hp_en, xlink_hp_exp_en, (void *)data, NULL, NULL, VRNA_DECOMP_PAIR_HP);
      ret &= vrna_sc_multi_cb_add(fc, xlink_il_en, xlink_il_exp_en, (void *)data, NULL, NULL, VRNA_DECOMP_PAIR_IL);
    }

    fclose(fp);
  }
  
  return ret;
}


PUBLIC int
vrna_sc_xl_prepare(vrna_fold_compound_t  *fc)
{
  if ((fc) &&
      (fc->sc) &&
      (fc->sc->xlink) &&
      (fc->sc->xlink_data)) {
    /* here we assume that sc->xlink_data actually points to a vrna_crosslink_data_s */
    
    struct vrna_crosslink_data_s  *dat;

    dat = (struct vrna_crosslink_data_s *)fc->sc->xlink_data;

    /* free previously allocated overlap entries */
    vrna_array_free(dat->overlap);

    /* initialize memory for overlap entries */
    vrna_array_init_size(dat->overlap, fc->length + 1);

    for (unsigned int i = 0; i <= fc->length; i++)
      dat->overlap[i]     = NULL;

    /* ref-fill overlap entries */
    for (size_t cnt = 0; cnt < vrna_array_size(dat->regions); cnt++) {
      interact_region_s *r = &(dat->regions[cnt]);

      for (unsigned int i = r->a; i <= r->b; i++)
        dat->overlap[i] = r;

      for (unsigned int j = r->c; j <= r->d; j++)
        dat->overlap[j] = r;

      if (fc->exp_params)
        r->exp_score = exp(-(FLT_OR_DBL)r->score * 10. / fc->exp_params->kT);
      else
        r->exp_score = 1.;
    }

    dat->A = (int *)vrna_alloc(sizeof(int) * (fc->length + 1));
    dat->B = (int *)vrna_alloc(sizeof(int) * (fc->length + 1));

    for (unsigned int i = 0; i <= fc->length; i++) {
      dat->A[i] = dat->B[i] = INF;
    }

    dat->qA = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * (fc->length + 1));
    dat->qB = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * (fc->length + 1));

    for (unsigned int i = 0; i <= fc->length; i++) {
      dat->qA[i] = dat->qB[i] = 0.;
    }
  }

  return 0;
}



PUBLIC int
vrna_mfe_multibranch_loop_xlink(vrna_fold_compound_t         *fc,
                                unsigned int                 i,
                                unsigned int                 j,
                                struct vrna_mx_mfe_aux_ml_s  *helpers)
{
  if ((fc) &&
      (helpers)) {

    unsigned char                 connects_intervals;
    int                           dangle_model, decomp, e;
    vrna_param_t                  *P;
    vrna_md_t                     *md;
    struct sc_mb_dat              sc_wrapper;
    struct vrna_crosslink_data_s  *dat;

    dat                 = (struct vrna_crosslink_data_s *)fc->sc->xlink_data;
    connects_intervals  = 0;
    P                   = fc->params;
    md                  = &(P->model_details);
    dangle_model        = md->dangles;
    e = decomp          = INF;

    /* hack using old API to prevent calling vrna_mfe_multibranch_loop_xlink in a loop */
    decomp = vrna_E_mb_loop_fast(fc,
                                 i,
                                 j,
                                 helpers->DMLi1,
                                 helpers->DMLi2);

    /* first, check whether (i,j) connects two intervals of a segment */
    if (dat->overlap[i]) {
      vrna_log_debug("pre_ml = %d", decomp);

      interact_region_s *r = dat->overlap[i];
      if ((r->a <= i) &&
          (r->b >= i) &&
          (r->c <= j) &&
          (r->d >= j)) {

        init_sc_mb(fc, &sc_wrapper);

        switch (dangle_model) {
          case 0:
            decomp = ml_pair_d0(fc, i, j, helpers, &sc_wrapper);
            break;

          default:
            decomp = ml_pair_d2(fc, i, j, helpers, &sc_wrapper);
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
vrna_bt_multibranch_loop_xlink(vrna_fold_compound_t  *fc,
                               unsigned int          i,         
                               unsigned int          j,         
                               int                   en,        
                               vrna_bps_t            bp_stack,  
                               vrna_bts_t            bt_stack)  
{
  if ((fc) &&
      (bp_stack) &&
      (bt_stack))
    return bt_mb_loop_xl(fc, i, j, en, bp_stack, bt_stack);

  return 0;
}


PRIVATE unsigned int
bt_mb_loop_xl(vrna_fold_compound_t *fc,
              unsigned int         i,                     
              unsigned int         j,                     
              int                  en,                   
              vrna_bps_t           bp_stack VRNA_UNUSED,  
              vrna_bts_t           bt_stack)
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
  unsigned char use_default_bt = 1;

  if ((sc) &&
      (sc->xlink)) {
    dat = (struct vrna_crosslink_data_s *)sc->xlink_data;
    r   = dat->overlap[i];

    if ((r) &&
        (r->a <= i) &&
        (r->b >= i) &&
        (r->c <= j) &&
        (r->d >= j))
      use_default_bt = 0;
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
  ptype           = (fc->type == VRNA_FC_TYPE_SINGLE) ? (sliding_window ? NULL : fc->ptype) : NULL;
  ptype_local     = (sliding_window) ? fc->ptype_local : NULL;
  rtype           = &(md->rtype[0]);
  type            = (fc->type == VRNA_FC_TYPE_SINGLE) ?
                    (sliding_window ? rtype[vrna_get_ptype_window(i, j,
                                                                  ptype_local)] : rtype[
                       vrna_get_ptype(ij, ptype)]) :
                    0;
  tt            = NULL;
  dangle_model  = md->dangles;
  turn          = md->min_loop_size;
  hc            = fc->hc;
  evaluate      = hc->eval_mb;
  connects_intervals  = 0;

  init_sc_mb(fc, &sc_wrapper);

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
    
      short *S, *S2, si1, sj1;
      unsigned int strands, *sn, tt;
      int stem, decomp, score;
    
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
                  .ml = c1}));
            return bt_A(fc, A, i, k, bp_stack, bt_stack);
          }
          
          if (e == B[k] + tmp_en) {
            free(tt);
            free_sc_mb(&sc_wrapper);
            vrna_bts_push(bt_stack,
                ((vrna_sect_t){
                  .i = k + 1,
                  .j = q,
                  .ml = c1}));
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
                    .ml = c1}));
    vrna_bts_push(bt_stack,
                  ((vrna_sect_t){
                    .i = k,
                    .j = q,
                    .ml = c2}));
    return 1;
  }

  vrna_log_error("bt M failed k = %d, l = %d, i  %d, j = %d", k, l, i, j);

  return 0;
}


PRIVATE unsigned int
bt_A(vrna_fold_compound_t *fc,
     int                  *A,
     unsigned int         i,                     
     unsigned int         v,                     
     vrna_bps_t           bp_stack VRNA_UNUSED,  
     vrna_bts_t           bt_stack)
{
  short                 *S, *S2, si1, sj1;
  unsigned int          u, *sn, strands, tt;
  int                   *my_c, *my_fML, *idx, stem, decomp, dangle_model, score;
  vrna_param_t          *P;
  vrna_md_t             *md;
  vrna_hc_t             *hc;
  vrna_sc_t             *sc;
  vrna_hc_eval_loop_f   evaluate;
  struct sc_mb_dat      sc_wrapper;
  struct vrna_crosslink_data_s  *dat;

  P         = fc->params;
  md        = &(P->model_details);
  dangle_model = md->dangles;
  idx       = fc->jindx;
  my_c      = fc->matrices->c;
  my_fML    = fc->matrices->fML;
  hc        = fc->hc;
  evaluate  = hc->eval_mb;
  sc        = fc->sc;
  dat       = (struct vrna_crosslink_data_s *)fc->sc->xlink_data;

  init_sc_mb(fc, &sc_wrapper);

  interact_region_s *r = dat->overlap[i];

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

exit_bt_A:

  A += i + 1;
  free(A);

  free_sc_mb(&sc_wrapper);

  if (u < v) {
    vrna_bts_push(bt_stack,
                  ((vrna_sect_t){
                    .i = i + 1,
                    .j = u - 1,
                    .ml = VRNA_MX_FLAG_M}));
    vrna_bts_push(bt_stack,
                  ((vrna_sect_t){
                    .i = u,
                    .j = v,
                    .ml = VRNA_MX_FLAG_C}));

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
     vrna_bps_t           bp_stack VRNA_UNUSED,  
     vrna_bts_t           bt_stack)
{
  short                 *S, *S2, si1, sj1;
  unsigned int          u, *sn, strands, tt;
  int                   *my_c, *idx, stem, decomp, dangle_model, score;
  vrna_param_t          *P;
  vrna_md_t             *md;
  vrna_hc_t             *hc;
  vrna_sc_t             *sc;
  vrna_hc_eval_loop_f   evaluate;
  struct sc_mb_dat      sc_wrapper;
  struct vrna_crosslink_data_s  *dat;

  P         = fc->params;
  md        = &(P->model_details);
  dangle_model = md->dangles;
  idx       = fc->jindx;
  my_c      = fc->matrices->c;
  hc        = fc->hc;
  evaluate  = hc->eval_mb;
  sc        = fc->sc;
  dat       = (struct vrna_crosslink_data_s *)fc->sc->xlink_data;

  init_sc_mb(fc, &sc_wrapper);

  interact_region_s *r = dat->overlap[i];
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

exit_bt_B:

  B += i + 1;
  free(B);

  free_sc_mb(&sc_wrapper);

  if (u < v) {
    vrna_bts_push(bt_stack,
                  ((vrna_sect_t){
                    .i = u,
                    .j = v,
                    .ml = VRNA_MX_FLAG_C}));

    return 1;
  }

  vrna_log_error("bt B failed");
  return 0;
}


PRIVATE unsigned int
bt_m2(vrna_fold_compound_t *fc,
      unsigned int         i,                     
      unsigned int         j,                     
      int                  en,                   
      vrna_bps_t           bp_stack VRNA_UNUSED,  
      vrna_bts_t           bt_stack)
{
  unsigned int          k;
  int                   *idx, *my_fML, tmp_en;
  vrna_hc_t             *hc;
  vrna_sc_t             *sc;
  vrna_hc_eval_loop_f   evaluate;
  struct sc_mb_dat      sc_wrapper;

  idx           = fc->jindx;
  my_fML        = fc->matrices->fML;
  hc            = fc->hc;
  evaluate      = hc->eval_mb;

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
                    .ml = VRNA_MX_FLAG_M}));
    vrna_bts_push(bt_stack,
                  ((vrna_sect_t){
                    .i = k + 1,
                    .j = j,
                    .ml = VRNA_MX_FLAG_M}));
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
    e = mfe_M2_xlink(fc, i, j, helpers, sc_wrapper);

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
    e = mfe_M2_xlink(fc, i, j, helpers, sc_wrapper);

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
    [i + 1, j - 1] with closing pair (i,j), given that (i,j) connects an
    interval pair
*/
PRIVATE INLINE int
mfe_M2_xlink(vrna_fold_compound_t         *fc,
             unsigned int                 i,
             unsigned int                 j,
             struct vrna_mx_mfe_aux_ml_s  *helpers,
             struct sc_mb_dat             *sc_wrapper)
{
  unsigned int                  v, k, p, q;
  int                           e, decomp, *fML, *fM2, *A, *B, *idx;
  vrna_param_t                  *P;
  vrna_hc_t                     *hc;
  struct vrna_crosslink_data_s  *dat;

  dat = (struct vrna_crosslink_data_s *)fc->sc->xlink_data;
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

  interact_region_s *r = dat->overlap[i];

  if ((r->a <= i) &&
      (r->b >= i) &&
      (r->c <= j) &&
      (r->d >= j)) {

    k = r->c;

    /* cases where with no base pair (u, v) with u < c and c <= v < d */
    if ((k > p + 1) &&
        (k < q) &&
        (hc->eval_mb(p, q, k - 1, k, VRNA_DECOMP_ML_ML_ML, hc))) {
      if ((fML[idx[k - 1] + p] != INF) &&
          (fML[idx[q] + k] != INF)) {
        decomp =  fML[idx[k - 1] + p] +
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
        decomp =  fM2[idx[k - 1] + p] +
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
        decomp =  fM2[idx[q] + k] +
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
vrna_sc_xl_update_mfe(vrna_fold_compound_t *fc,
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

  dat = (struct vrna_crosslink_data_s *)fc->sc->xlink_data;
  n   = fc->length;
  P   = fc->params;
  md  = &(P->model_details);
  dangle_model = md->dangles;
  turn  = md->min_loop_size;
  idx = fc->jindx;
  c   = fc->matrices->c;
  fML = fc->matrices->fML;
  hc  = fc->hc;
  A   = dat->A;
  B   = dat->B;

  init_sc_mb(fc, &sc_wrapper);

  /* i overlaps with I, add all ML parts where a base pair starts in I and pairs outside of I */

  for (unsigned int v = i; v <= n; v++)
    /* reset A and B */
    A[v] = B[v] = INF;

  interact_region_s *r = dat->overlap[i];

  if ((r) &&
      (i >= r->a) &&
      (i <= r->b)) {
    score = r->score;
    vrna_log_debug("score = %d", score);
    /*
        collect energies of ML parts where a base pair (u, v) crosses I', i.e.
        u < c && c <= v < d
     */
    for (unsigned int v = r->c; v < r->d; v++) {
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
  int                           energy = 0;
  struct vrna_crosslink_data_s  *dat = (struct vrna_crosslink_data_s *)data;

  /* go through each interval that overlaps with position i to check whether j connects both segments */
  if (dat->overlap[i]) {
    interact_region_s *r = dat->overlap[i];
    if ((r->a <= i) &&
        (r->b >= i) &&
        (r->b < j) &&
        (r->c <= j) &&
        (r->d >= j))
      energy += r->score;
  }

  return energy;

}


INLINE PRIVATE FLT_OR_DBL
xlink_hp_exp_en(vrna_fold_compound_t  *fc,
                int                   i,
                int                   j,
                int                   k,
                int                   l,
                void                  *data)
{
  return exp(-(FLT_OR_DBL)xlink_hp_en(fc, i, j, k, l, data) * 10. / fc->exp_params->kT);
}


INLINE PRIVATE int
xlink_il_en(vrna_fold_compound_t  *fc,
            int                   i,
            int                   j,
            int                   k,
            int                   l,
            void                  *data)
{
  int                           energy = 0;
  struct vrna_crosslink_data_s  *dat = (struct vrna_crosslink_data_s *)data;

  /*  go through each interval that overlaps with position i to check whether
   *  k does not overlap with the 5' interval and j connects both segments */
  if (dat->overlap[i]) {
    interact_region_s *r = dat->overlap[i];
    if ((r->a <= i) &&
        (r->b >= i) &&
        (r->b < k) &&
        (r->c <= j) &&
        (r->d >= j))
      energy += r->score;
  }

  return energy;

}


INLINE PRIVATE FLT_OR_DBL
xlink_il_exp_en(vrna_fold_compound_t  *fc,
                int                   i,
                int                   j,
                int                   k,
                int                   l,
                void                  *data)
{
  return exp(-(FLT_OR_DBL)xlink_il_en(fc, i, j, k, l, data) * 10. / fc->exp_params->kT);
}


PRIVATE interact_region_s *
parse_crosslink_data(FILE *fp)
{
  char                          *line;
  int                           a, b, c, d;
  double                        score;
  vrna_array(interact_region_s) regions;

  vrna_array_init_size(regions, 42);

  while ((line = vrna_read_line(fp))) {
    if (sscanf(line, "%d %d %d %d %lf", &a, &b, &c, &d, &score) == 5) {
      interact_region_s r = {
        .a = a,
        .b = b,
        .c = c,
        .d = d,
        .score = (int)(score * 100.)
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

    dat = (struct vrna_crosslink_data_s *)fc->sc->xlink_data;

    /* first, check whether (i,j) connects two intervals of a segment */
    if ((r = dat->overlap[j]) &&
        (r->a <= i) &&
        (r->b >= i) &&
        (r->c <= j) &&
        (r->d >= j)) {
      q = ml_pair_exp(fc, i, j, helpers);
      vrna_log_debug("mb_xl => %g [%u:%u] - [%u:%u]", q, r->a, r->b, r->c, r->d);
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
            struct vrna_mx_pf_aux_ml_s  *helpers)
{
  unsigned char             sliding_window;
  char                      *ptype, **ptype_local;
  short                     *S1, **SS, **S5, **S3;
  unsigned int              *sn, n_seq, s, *se, tt;
  int                       ij, *my_iindx, *jindx, *rtype;
  FLT_OR_DBL                qbt1, temp, qqqmmm, *qm, *qm2, **qm_local, *scale, expMLclosing, *qqm1,
                            *qqm21;
  vrna_hc_t                 *hc;
  vrna_exp_param_t          *pf_params;
  vrna_md_t                 *md;
  struct sc_mb_exp_dat      sc_wrapper;

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
        tt    = (sliding_window) ?
                rtype[vrna_get_ptype_window(i, j + i, ptype_local)] :
                rtype[vrna_get_ptype(ij, ptype)];
        qbt1  *= vrna_exp_E_multibranch_stem(tt, S1[j - 1], S1[i + 1], pf_params);

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

    /* get contribution for enclosed part */
    qbt1 *= pf_M2_xlink(fc, i, j, helpers, &sc_wrapper);

    free_sc_mb_exp(&sc_wrapper);
  }

  return qbt1;
}


PRIVATE INLINE FLT_OR_DBL
pf_M2_xlink(vrna_fold_compound_t        *fc,
            unsigned int                i,
            unsigned int                j,
            struct vrna_mx_pf_aux_ml_s  *helpers,
            struct sc_mb_exp_dat        *sc_wrapper)
{
  unsigned int                  v, u, k, l;
  int                           *idx;
  FLT_OR_DBL                    q, qbt1, *qm, *qm2, *expMLbase, *qA, *qB;
  vrna_hc_t                     *hc;
  struct vrna_crosslink_data_s  *dat;

  dat       = (struct vrna_crosslink_data_s *)fc->sc->xlink_data;
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

  interact_region_s *r = dat->overlap[i];

  u = r->b;

  /* cases where with no base pair (u, v) with u < c and c <= v < d */
  if ((u > k) &&
      (u < l) &&
      (hc->eval_mb(k, l, u, u + 1, VRNA_DECOMP_ML_ML_ML, hc))) {
    qbt1 =  qm[idx[k] - u] *
            qm[idx[u + 1] - l] *
            r->exp_score;

    if (sc_wrapper->decomp_ml)
      qbt1 *= sc_wrapper->decomp_ml(k, l, u, u + 1, sc_wrapper);

    q += qbt1;
    vrna_log_debug("q1 = %g", qbt1);
  }

  if ((u > k) &&
      (hc->eval_mb(k, l, k, u, VRNA_DECOMP_ML_ML, hc))) {
    qbt1 =  qm2[idx[k] - u] *
            expMLbase[l - u] *
            r->exp_score;

    if (sc_wrapper->red_ml)
      qbt1 *= sc_wrapper->red_ml(k, l, k, u, sc_wrapper);

    q += qbt1;
    vrna_log_debug("q2 = %g", qbt1);
  }

  if ((u < l) &&
      (hc->eval_mb(k, l, u + 1, l, VRNA_DECOMP_ML_ML, hc))) {
    qbt1 =  qm2[idx[u + 1] - l] *
            expMLbase[u - k + 1] *
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

    if (hc->eval_mb(k, l, u - 1, u, VRNA_DECOMP_ML_ML_ML, hc)) {
      qbt1 = qA[u] *
             qm[idx[k] - u + 1];

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

  return q;
}


PUBLIC void
vrna_sc_xl_update_pf(vrna_fold_compound_t *fc,
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

  dat           = (struct vrna_crosslink_data_s *)fc->sc->xlink_data;
  n             = fc->length;
  pf_params     = fc->exp_params;
  md            = &(fc->params->model_details);
  turn          = md->min_loop_size;
  idx           = fc->iindx;
  qb            = fc->exp_matrices->qb;
  qm            = fc->exp_matrices->qm;
  expMLbase     = fc->exp_matrices->expMLbase;
  hc            = fc->hc;
  qA            = dat->qA;
  qB            = dat->qB;

  init_sc_mb_exp(fc, &sc_wrapper);

  /* i overlaps with I, add all ML parts where a base pair starts in I and pairs outside of I */

  for (unsigned int v = 1; v <= j; v++)
    /* reset A and B */
    qA[v] = qB[v] = 0.;

  interact_region_s *r = dat->overlap[j];

  if ((r) &&
      (j >= r->c) &&
      (j <= r->d)) {
    score = r->exp_score;
    vrna_log_debug("exp_score = %g", score);
    /*
        collect energies of ML parts where a base pair (u, v) crosses I, i.e.
        a < u <= b && b < v < d
     */
    for (unsigned int u = r->a + 1; u <= r->b; u++) {
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

          vrna_log_debug("stem(%u,%u) = %g", u, v, vrna_exp_E_multibranch_stem(tt, si1, sj1, pf_params));

          if (v < r->c)
            qbt1 *= score;

          /* fill A */
          if ((v < j - 1) &&
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


