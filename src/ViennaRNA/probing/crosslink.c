/* Crosslinking data handling */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdio.h>
#include <stdlib.h>
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

#include "ViennaRNA/probing/crosslink.h"


/* re-definition of vrna_mx_mfe_aux_ml_s to use its internals, see also mfe_multibrnch.c  */
struct vrna_mx_mfe_aux_ml_s {
  unsigned int  n;
  int           *Fmi;
  int           *DMLi;
  int           *DMLi1;
  int           *DMLi2;
};


typedef struct {
  unsigned int  a;      /* start of 5' interval */
  unsigned int  b;      /* end of 5' interval */
  unsigned int  c;      /* start of 3' interval */
  unsigned int  d;      /* end of 3' interval */
  int           score;  /* interval score */
} interact_region_s;


struct vrna_crosslink_data_s {
  vrna_array(interact_region_s)   regions;  /* list of interacting regions */
  vrna_array(interact_region_s **) overlap5;
  vrna_array(interact_region_s **) overlap3;

  int *A;
  int *B;
};



PRIVATE INLINE int
ml_pair_d2(vrna_fold_compound_t         *fc,
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


PRIVATE int
xlink_hp_en(vrna_fold_compound_t  *fc,
            int                   i,
            int                   j,
            int                   k,
            int                   l,
            void                  *data);


PRIVATE int
xlink_il_en(vrna_fold_compound_t  *fc,
            int                   i,
            int                   j,
            int                   k,
            int                   l,
            void                  *data);


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
      data->overlap5      = NULL;
      data->overlap3      = NULL;
      fc->sc->xlink_data  = (void *)data;
      fc->params->model_details.uniq_ML = 1;

      ret &= vrna_sc_multi_cb_add(fc, xlink_hp_en, NULL, (void *)data, NULL, NULL, VRNA_DECOMP_PAIR_HP);
      ret &= vrna_sc_multi_cb_add(fc, xlink_il_en, NULL, (void *)data, NULL, NULL, VRNA_DECOMP_PAIR_IL);
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
    if (dat->overlap5) {
      for (unsigned int i = 1; i <= fc->length; i++)
        vrna_array_free(dat->overlap5[i]);

      vrna_array_free(dat->overlap5);
    }

    if (dat->overlap3) {
      for (unsigned int i = 1; i <= fc->length; i++)
        vrna_array_free(dat->overlap3[i]);

      vrna_array_free(dat->overlap3);
    }

    /* initialize memory for overlap entries */
    vrna_array_init_size(dat->overlap5, fc->length + 1);
    vrna_array_init_size(dat->overlap3, fc->length + 1);

    for (unsigned int i = 0; i <= fc->length; i++) {
      dat->overlap5[i] = NULL;
      dat->overlap3[i] = NULL;
    }

    /* ref-fill overlap entries */
    for (size_t cnt = 0; cnt < vrna_array_size(dat->regions); cnt++) {
      interact_region_s *r = &(dat->regions[cnt]);
      
      for (unsigned int i = r->a; i <= r->b; i++) {
        if (dat->overlap5[i] == NULL)
          vrna_array_init(dat->overlap5[i]);

        vrna_array_append(dat->overlap5[i], r);
      }

      for (unsigned int j = r->c; j <= r->d; j++) {
        if (dat->overlap3[j] == NULL)
          vrna_array_init(dat->overlap3[j]);

        vrna_array_append(dat->overlap3[j], r);
      }
    }

    dat->A = (int *)vrna_alloc(sizeof(int) * (fc->length + 1));
    dat->B = (int *)vrna_alloc(sizeof(int) * (fc->length + 1));

    for (unsigned int i = 0; i <= fc->length; i++) {
      dat->A[i] = dat->B[i] = INF;
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

    /* first, check whether (i,j) connects two intervals of a segment */
    if (dat->overlap5[i]) {
      for (size_t cnt = 0; cnt < vrna_array_size(dat->overlap5[i]); cnt++) {
        interact_region_s *r = dat->overlap5[i][cnt];
        if ((r->c <= j) &&
            (r->d >= j)) {
          connects_intervals = 1;
          break;
        }
      }

      if (connects_intervals) {
        init_sc_mb(fc, &sc_wrapper);

        switch (dangle_model) {
          /* currently, only dangles = 2 implemented */
          default:
            decomp = ml_pair_d2(fc, i, j, helpers, &sc_wrapper);
            break;
        }

        free_sc_mb(&sc_wrapper);

        return decomp;
      }
    }

    /* hack using old API to prevent calling vrna_mfe_multibranch_loop_xlink in a loop */
    return vrna_E_mb_loop_fast(fc,
                               i,
                               j,
                               helpers->DMLi1,
                               helpers->DMLi2);

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
  unsigned char         sliding_window;
  char                  *ptype, **ptype_local;
  short                 s5, s3, *S1, **SS, **S5, **S3;
  unsigned int          *sn, n_seq, s, *tt, c1, c2, p, q, k, l, dangle_model,
                        type, type_2, turn, connects_intervals;
  int                   ij, e, tmp_en, *idx, *my_c, *my_fML, *my_fM2, *rtype, **c_local,
                        **fML_local;
  vrna_param_t          *P;
  vrna_md_t             *md;
  vrna_hc_t             *hc;
  vrna_sc_t             *sc;
  vrna_hc_eval_loop_f   evaluate;
  struct sc_mb_dat      sc_wrapper;
  struct vrna_crosslink_data_s  *dat;

  sc = (fc->type == VRNA_FC_TYPE_SINGLE) ? fc->sc : NULL;

  /* sanity check to test wether we actually need this dedicated backtracking implementation */
  if ((sc) &&
      (!(sc->xlink)))
    return vrna_bt_multibranch_loop(fc, i, j, en, bp_stack, bt_stack);


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
  dat           = (struct vrna_crosslink_data_s *)sc->xlink_data;

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

    /* check, whether (i,j) can actually connect any pair of intervals */
    if (dat->overlap5[i]) {
      for (size_t cnt = 0; cnt < vrna_array_size(dat->overlap5[i]); cnt++) {
        interact_region_s *r = dat->overlap5[i][cnt];
        if ((r->c <= j) &&
            (r->d >= j)) {
          connects_intervals = 1;
          break;
        }
      }
    }

    if (connects_intervals) {
      /*  now, check for split point between the enclosed branches */
      for (size_t cnt = 0; cnt < vrna_array_size(dat->overlap5[i]); cnt++) {
        interact_region_s *r = dat->overlap5[i][cnt];

        /* only consider interval pairs (I,I') where I' overlaps with j */
        if ((r->c <= q) &&
            (r->d > q)) {

          k = r->b;

          /* first 3 (obvious) cases: split point separates the intervals */
          if ((k > p) &&
              (evaluate(p, q, k, k + 1, VRNA_DECOMP_ML_ML_ML, hc))) {
            tmp_en = my_fML[idx[k] + p] +
                     my_fML[idx[q] + k + 1] +
                     r->score;

            if (sc_wrapper.decomp_ml)
              tmp_en += sc_wrapper.decomp_ml(p, q, k, k + 1, &sc_wrapper);

            if (e == tmp_en)
              goto odd_dangles_exit;
          }

          if (evaluate(p, q, p, k, VRNA_DECOMP_ML_ML, hc)) {
            l = p;
            tmp_en = my_fM2[idx[k] + p] +
                     (q - k) * P->MLbase +
                     r->score;

            if (sc_wrapper.red_ml)
              tmp_en += sc_wrapper.red_ml(p, q, p, k, &sc_wrapper);

            if (e == tmp_en)
              goto odd_dangles_exit;
          }

          if (evaluate(p, q, k + 1, q, VRNA_DECOMP_ML_ML, hc)) {
            l = q;
            tmp_en = my_fM2[idx[q] + k + 1] +
                     (k + 1 - p) * P->MLbase +
                     r->score;

            if (sc_wrapper.red_ml)
              tmp_en += sc_wrapper.red_ml(p, q, k + 1, q, &sc_wrapper);

            if (e == tmp_en)
              goto odd_dangles_exit;
          }
        }
      }

      /* nothing found so far, let's fill A and B matrices to allow for search of non-obious configurations */
      int *A = NULL, *B = NULL;
      
      {
      
      short *S, *S2, si1, sj1;
      unsigned int strands, *sn, tt;
      int stem, decomp;
      
      A = (int *)vrna_alloc(sizeof(int) * (j - i));
      B = (int *)vrna_alloc(sizeof(int) * (j - i));

      A -= p; /* pointer magic */
      B -= p; /* pointer magic */

      for (l = p; l <= q; l++) {
        /* reset A and B */
        A[l] = B[l] = INF;

        /* got through all pairing partners k */
        for (unsigned int k = p; k < l; k++) {
          if ((hc->eval_mb(k, l, k, l, VRNA_DECOMP_ML_STEM, hc)) &&
              (my_c[idx[l] + k] != INF)) {
            strands = fc->strands;
            sn      = fc->strand_number;
            S       = fc->sequence_encoding;
            S2      = fc->sequence_encoding2;
            tt      = vrna_get_ptype_md(S2[k], S2[l], md);

            si1 = ((strands == 1) || (sn[k - 1] == sn[k])) ? S[k - 1] : -1;
            sj1 = ((strands == 1) || (sn[l] == sn[l + 1])) ? S[l + 1] : -1;

            stem = vrna_E_multibranch_stem(tt, si1, sj1, P) +
                   my_c[idx[l] + k];

            decomp = 0;

            /* go through all intervals I' that cover l */
            for (size_t cnt = 0; cnt < vrna_array_size(dat->overlap3[l]); cnt++) {
              interact_region_s *r = dat->overlap3[l][cnt];

              /* add up scores if pairing partner k doesn't connect both intervals */
              if ((k < r->c) &&
                  (k > r->b))
                decomp += r->score;
            }

            decomp += stem;

            /* fill A */
            if (hc->eval_mb(p, l, k - 1, k, VRNA_DECOMP_ML_ML_STEM, hc)) {
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



      for (size_t cnt = 0; cnt < vrna_array_size(dat->overlap5[i]); cnt++) {
        interact_region_s *r = dat->overlap5[i][cnt];

        /* only consider interval pairs (I,I') where I' overlaps with j */
        if ((r->c <= q) &&
            (r->d > q)) {

          /* now, for the cases where we can't tell yet whether (i,j) is the innermost pair connecting I and I' */
          tmp_en = A[q];

          if (e == tmp_en) {
            free(tt);
            free_sc_mb(&sc_wrapper);
            return bt_A(fc, A, p, q, bp_stack, bt_stack);
          }

          for (k = r->b + 1; k < q; k++) {
            if ((evaluate(p, q, p, k, VRNA_DECOMP_ML_ML, hc)) &&
                (A[k] != INF)) {
              tmp_en = A[k] +
                       (q - k) * P->MLbase;

              if (sc_wrapper.red_ml)
                tmp_en += sc_wrapper.red_ml(p, q, p, k, &sc_wrapper);

              if (e == tmp_en) {
                free(tt);
                free_sc_mb(&sc_wrapper);
                return bt_A(fc, A, p, k, bp_stack, bt_stack);
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
                return bt_A(fc, A, p, k, bp_stack, bt_stack);
              }
              
              if (e == B[k] + tmp_en) {
                free(tt);
                free_sc_mb(&sc_wrapper);
                vrna_bts_push(bt_stack,
                    ((vrna_sect_t){
                      .i = k + 1,
                      .j = q,
                      .ml = c1}));
                return bt_B(fc, B, p, k, bp_stack, bt_stack);
              }
            }
          }
        }
      }
    }

    /* regular mb-loop decomposition for dangles % 2 == 0 */
    for (k = i + 2; k < j - 2; ++k) {
      if (evaluate(p, q, k, k + 1, VRNA_DECOMP_ML_ML_ML, hc)) {
        if (sliding_window)
          tmp_en = fML_local[p][k - p] +
                   fML_local[k + 1][q - (k + 1)];
        else
          tmp_en = my_fML[idx[k] + p] +
                   my_fML[idx[q] + k + 1];

        if (sc_wrapper.decomp_ml)
          tmp_en += sc_wrapper.decomp_ml(p, q, k, k + 1, &sc_wrapper);

        if (e == tmp_en)
          goto odd_dangles_exit;
      }
    }
    
  }


odd_dangles_exit:

  free(tt);
  free_sc_mb(&sc_wrapper);

  if (k + 2 < j) {
    /* find split point for M2 entry */
    if (l == p) {
      return bt_m2(fc, p, k, my_fM2[idx[k] + p], bp_stack, bt_stack);
    } else if (l == q) {
      return bt_m2(fc, k + 1, q, my_fM2[idx[q] + k + 1], bp_stack, bt_stack);
    } else {
      /* already found split point */
      vrna_bts_push(bt_stack,
                    ((vrna_sect_t){
                      .i = p,
                      .j = k,
                      .ml = c1}));
      vrna_bts_push(bt_stack,
                    ((vrna_sect_t){
                      .i = k + 1,
                      .j = q,
                      .ml = c2}));
    }

    return 1;
  }

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
  int                   *my_c, *my_fML, *idx, stem, decomp;
  vrna_param_t          *P;
  vrna_md_t             *md;
  vrna_hc_t             *hc;
  vrna_sc_t             *sc;
  vrna_hc_eval_loop_f   evaluate;
  struct sc_mb_dat      sc_wrapper;
  struct vrna_crosslink_data_s  *dat;

  P         = fc->params;
  md        = &(P->model_details);
  idx       = fc->jindx;
  my_c      = fc->matrices->c;
  my_fML    = fc->matrices->fML;
  hc        = fc->hc;
  evaluate  = hc->eval_mb;
  sc        = fc->sc;
  dat       = (struct vrna_crosslink_data_s *)fc->sc->xlink_data;

  init_sc_mb(fc, &sc_wrapper);


  for (u = i + 1; u < v; u++) {
    if ((evaluate(i, v, u - 1, u, VRNA_DECOMP_ML_ML_STEM, hc)) &&
        (my_c[idx[v] + u] != INF)) {
      strands = fc->strands;
      sn      = fc->strand_number;
      S       = fc->sequence_encoding;
      S2      = fc->sequence_encoding2;
      tt      = vrna_get_ptype_md(S2[u], S2[v], md);

      si1 = ((strands == 1) || (sn[u - 1] == sn[u])) ? S[u - 1] : -1;
      sj1 = ((strands == 1) || (sn[v] == sn[v + 1])) ? S[v + 1] : -1;

      stem = vrna_E_multibranch_stem(tt, si1, sj1, P) +
             my_c[idx[v] + u] +
             my_fML[idx[u - 1] + i];

      decomp = 0;

      /* go through all intervals I' that cover l */
      for (size_t cnt = 0; cnt < vrna_array_size(dat->overlap3[v]); cnt++) {
        interact_region_s *r = dat->overlap3[v][cnt];

        /* add up scores if pairing partner k doesn't connect both intervals */
        if ((u < r->c) &&
            (u > r->b))
          decomp += r->score;
      }

      decomp += stem;

      if (sc_wrapper.decomp_ml)
        decomp += sc_wrapper.decomp_ml(i, v, u - 1, u, &sc_wrapper);

      if (A[v] == decomp)
        goto exit_bt_A;
    }
  }

exit_bt_A:

  A += i;
  free(A);

  free_sc_mb(&sc_wrapper);

  if (u < v) {
    vrna_bts_push(bt_stack,
                  ((vrna_sect_t){
                    .i = i,
                    .j = u - 1,
                    .ml = VRNA_MX_FLAG_M}));
    vrna_bts_push(bt_stack,
                  ((vrna_sect_t){
                    .i = u,
                    .j = v,
                    .ml = VRNA_MX_FLAG_C}));

    return 1;
  }

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
  int                   *my_c, *idx, stem, decomp;
  vrna_param_t          *P;
  vrna_md_t             *md;
  vrna_hc_t             *hc;
  vrna_sc_t             *sc;
  vrna_hc_eval_loop_f   evaluate;
  struct sc_mb_dat      sc_wrapper;
  struct vrna_crosslink_data_s  *dat;

  P         = fc->params;
  md        = &(P->model_details);
  idx       = fc->jindx;
  my_c      = fc->matrices->c;
  hc        = fc->hc;
  evaluate  = hc->eval_mb;
  sc        = fc->sc;
  dat       = (struct vrna_crosslink_data_s *)fc->sc->xlink_data;

  init_sc_mb(fc, &sc_wrapper);


  for (u = i; u < v; u++) {
    if ((evaluate(i, v, u, v, VRNA_DECOMP_ML_ML, hc)) &&
        (my_c[idx[v] + u] != INF)) {
      strands = fc->strands;
      sn      = fc->strand_number;
      S       = fc->sequence_encoding;
      S2      = fc->sequence_encoding2;
      tt      = vrna_get_ptype_md(S2[u], S2[v], md);

      si1 = ((strands == 1) || (sn[u - 1] == sn[u])) ? S[u - 1] : -1;
      sj1 = ((strands == 1) || (sn[v] == sn[v + 1])) ? S[v + 1] : -1;

      stem = vrna_E_multibranch_stem(tt, si1, sj1, P) +
             my_c[idx[v] + u] +
             (u - i) * P->MLbase;

      decomp = 0;

      /* go through all intervals I' that cover l */
      for (size_t cnt = 0; cnt < vrna_array_size(dat->overlap3[v]); cnt++) {
        interact_region_s *r = dat->overlap3[v][cnt];

        /* add up scores if pairing partner k doesn't connect both intervals */
        if ((u < r->c) &&
            (u > r->b))
          decomp += r->score;
      }

      decomp += stem;

      if (sc_wrapper.red_ml)
        decomp += sc_wrapper.red_ml(i, v, u, v, &sc_wrapper);

      if (B[v] == decomp)
        goto exit_bt_B;
    }
  }

exit_bt_B:

  B += i;
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

  return 0;
}


PRIVATE INLINE int
ml_pair_d2(vrna_fold_compound_t         *fc,
           unsigned int                 i,
           unsigned int                 j,
           struct vrna_mx_mfe_aux_ml_s  *helpers,
           struct sc_mb_dat             *sc_wrapper)
{
  short         *S, *S2, **SS, **S5, **S3, si1, sj1;
  unsigned int  tt, strands, *sn, s, n_seq;
  int           e;
  vrna_param_t  *P;
  vrna_md_t     *md;
  vrna_hc_t     *hc;

  e   = INF;
  hc  = fc->hc;

  if (hc->eval_mb(i, j, i + 1, j - 1, VRNA_DECOMP_PAIR_ML, hc)) {

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

          si1 = ((strands == 1) || (sn[i] == sn[i + 1])) ? S[i + 1] : -1;
          sj1 = ((strands == 1) || (sn[j - 1] == sn[j])) ? S[j - 1] : -1;

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
  unsigned int                  v;
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

  for (size_t cnt = 0; cnt < vrna_array_size(dat->overlap5[i]); cnt++) {
    interact_region_s *r = dat->overlap5[i][cnt];
    if ((r->c <= j) &&
        (r->d >= j)) {

      if (hc->eval_mb(i + 1, j - 1, r->b, r->b + 1, VRNA_DECOMP_ML_ML_ML, hc)) {
        if ((fML[idx[r->b] + i + 1] != INF) &&
            (fML[idx[j - 1] + r->b + 1] != INF)) {
          decomp =  fML[idx[r->b] + i + 1] +
                    fML[idx[j - 1] + r->b + 1] +
                    r->score;

          if (sc_wrapper->decomp_ml)
            decomp += sc_wrapper->decomp_ml(i + 1, j - 1, r->b, r->b + 1, sc_wrapper);

          e = MIN2(e, decomp);
        }
      }

      if (hc->eval_mb(i + 1, j - 1, i + 1, r->b, VRNA_DECOMP_ML_ML, hc)) {
        if (fM2[idx[r->b] + i + 1] != INF) {
          decomp =  fM2[idx[r->b] + i + 1] +
                    (j - 1 - r->b) * P->MLbase +
                    r->score;

          if (sc_wrapper->red_ml)
            decomp += sc_wrapper->red_ml(i + 1, j - 1, i + 1, r->b, sc_wrapper);

          e = MIN2(e, decomp);
        }
      }

      if (hc->eval_mb(i + 1, j - 1, r->b + 1, j - 1, VRNA_DECOMP_ML_ML, hc)) {
        if (fM2[idx[j - 1] + r->b + 1] != INF) {
          decomp =  fM2[idx[j - 1] + r->b + 1] +
                    (r->b - i) * P->MLbase +
                    r->score;

          if (sc_wrapper->red_ml)
            decomp += sc_wrapper->red_ml(i + 1, j - 1, r->b + 1, j - 1, sc_wrapper);

          e = MIN2(e, decomp);
        }
      }

      /* now, all cases where we allow for base pairs that cross the first interval boundary */
      decomp = A[j - 1];
      e = MIN2(e, decomp);

      for (v = r->c; v + 1 < j; v++) {
        if (hc->eval_mb(i + 1, j - 1, v, v + 1, VRNA_DECOMP_ML_ML, hc)) {
          decomp = A[v];
          if (decomp != INF) {
            decomp += (j - 1 - v) * P->MLbase;

            if (sc_wrapper->red_ml)
              decomp += sc_wrapper->red_ml(i + 1, j - 1, v, v + 1, sc_wrapper);

            e = MIN2(e, decomp);
          }
        }

        if (hc->eval_mb(i + 1, j - 1, v, v + 1, VRNA_DECOMP_ML_ML_ML, hc)) {
          decomp = A[v];

          if (decomp != INF) {
            decomp += fML[idx[j - 1] + v + 1];

            if (sc_wrapper->decomp_ml)
              decomp += sc_wrapper->decomp_ml(i + 1, j - 1, v, v + 1, sc_wrapper);

            e = MIN2(e, decomp);
          }

          decomp = B[v];

          if (decomp != INF) {
            decomp += fML[idx[j - 1] + v + 1];

            if (sc_wrapper->decomp_ml)
              decomp += sc_wrapper->decomp_ml(i + 1, j - 1, v, v + 1, sc_wrapper);

            e = MIN2(e, decomp);
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
  unsigned int                  n, stem, decomp, en, *sn, strands, tt;
  int                           *A, *B, *fML, *c, *idx;
  vrna_param_t                  *P;
  vrna_md_t                     *md;
  vrna_hc_t                     *hc;
  struct vrna_crosslink_data_s  *dat;
  struct sc_mb_dat              sc_wrapper;

  dat = (struct vrna_crosslink_data_s *)fc->sc->xlink_data;
  n   = fc->length;
  P   = fc->params;
  md  = &(P->model_details);
  idx = fc->jindx;
  c   = fc->matrices->c;
  fML = fc->matrices->fML;
  hc  = fc->hc;
  A   = dat->A;
  B   = dat->B;

  init_sc_mb(fc, &sc_wrapper);

  for (unsigned int v = i; v <= n; v++) {
    /* reset A and B */
    A[v] = B[v] = INF;

    /* got through all pairing partners u */
    for (unsigned int u = i + 1; u < v; u++) {
      if ((hc->eval_mb(u, v, u, v, VRNA_DECOMP_ML_STEM, hc)) &&
          (c[idx[v] + u] != INF)) {
        strands = fc->strands;
        sn      = fc->strand_number;
        S       = fc->sequence_encoding;
        S2      = fc->sequence_encoding2;
        tt      = vrna_get_ptype_md(S2[u], S2[v], md);

        si1 = ((strands == 1) || (sn[u - 1] == sn[u])) ? S[u - 1] : -1;
        sj1 = ((strands == 1) || (sn[v] == sn[v + 1])) ? S[v + 1] : -1;

        stem = vrna_E_multibranch_stem(tt, si1, sj1, P) +
               c[idx[v] + u];

        decomp = 0;

        /* go through all intervals I' that cover v */
        if (dat->overlap3[v]) {
          for (size_t cnt = 0; cnt < vrna_array_size(dat->overlap3[v]); cnt++) {
            interact_region_s *r = dat->overlap3[v][cnt];

            /* add up scores if pairing partner u doesn't connect both intervals */
            if ((u < r->c) &&
                (u > r->b))
              decomp += r->score;
          }
        }

        decomp += stem;

        /* fill A */
        if (hc->eval_mb(i + 1, v, u - 1, u, VRNA_DECOMP_ML_ML_STEM, hc)) {
          if (fML[idx[u - 1] + i + 1] != INF) {
            en = decomp +
                 fML[idx[u - 1] + i + 1];

            if (sc_wrapper.decomp_ml)
              en += sc_wrapper.decomp_ml(i + 1, v, u - 1, u, &sc_wrapper);

            A[v] = MIN2(A[v], en);
          }
        }

        /* fill B */
        if (hc->eval_mb(i + 1, v, u, v, VRNA_DECOMP_ML_ML, hc)) {
          en = decomp +
               (u - 1 - i) * P->MLbase;
        
          if (sc_wrapper.red_ml)
            en += sc_wrapper.red_ml(i + 1, v, u, v, &sc_wrapper);

          B[v] = MIN2(B[v], en);
        }
      }
    }
  }

  free_sc_mb(&sc_wrapper);
}


PRIVATE int
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
  if (dat->overlap5[i])
    for (size_t cnt = 0; cnt < vrna_array_size(dat->overlap5[i]); cnt++) {
      interact_region_s *r = dat->overlap5[i][cnt];
      if ((r->b < j) &&
          (r->c <= j) &&
          (r->d >= j))
      energy += r->score;
    }

  return energy;

}


PRIVATE int
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
  if (dat->overlap5[i])
    for (size_t cnt = 0; cnt < vrna_array_size(dat->overlap5[i]); cnt++) {
      interact_region_s *r = dat->overlap5[i][cnt];
      if ((r->b < k) &&
          (r->c <= j) &&
          (r->d >= j))
      energy += r->score;
    }

  return energy;

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
