/*
 * WBL 24 Aug 2018 Add AVX512 based on sources_034_578/modular_decomposition_id3.c
 * WBL 22 Aug 2018 by hand d3c17fd3e04e2419c147a1e097d3c4d2c5a6f11d lines 1355-1357
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <string.h>
#include "ViennaRNA/utils/basic.h"
#include "ViennaRNA/fold_vars.h"
#include "ViennaRNA/alphabet.h"
#include "ViennaRNA/params/default.h"
#include "ViennaRNA/constraints/hard.h"
#include "ViennaRNA/constraints/soft.h"
#include "ViennaRNA/loops/external.h"
#include "ViennaRNA/gquad.h"
#include "ViennaRNA/structured_domains.h"
#include "ViennaRNA/unstructured_domains.h"
#include "ViennaRNA/loops/multibranch.h"
#include "ViennaRNA/utils/higher_order_functions.h"

#ifdef __GNUC__
# define INLINE inline
#else
# define INLINE
#endif

#include "multibranch_hc.inc"
#include "multibranch_sc.inc"

/*
 #################################
 # PRIVATE FUNCTION DECLARATIONS #
 #################################
 */

PRIVATE int
E_mb_loop_fast(vrna_fold_compound_t *fc,
               int                  i,
               int                  j,
               int                  *dmli1,
               int                  *dmli2);


PRIVATE int
E_ml_stems_fast(vrna_fold_compound_t  *fc,
                int                   i,
                int                   j,
                int                   *fmi,
                int                   *dmli);


PRIVATE int
extend_fm_3p(int                        i,
             int                        j,
             int                        *fm,
             vrna_fold_compound_t       *fc,
             vrna_callback_hc_evaluate  *evaluate,
             struct hc_mb_def_dat       *hc_data,
             struct sc_mb_dat           *sc_wrapper);


PRIVATE int
E_mb_loop_stack(vrna_fold_compound_t  *fc,
                int                   i,
                int                   j);


PRIVATE INLINE int
ml_pair5(vrna_fold_compound_t       *fc,
         int                        i,
         int                        j,
         int                        *dmli2,
         vrna_callback_hc_evaluate  *evaluate,
         struct hc_mb_def_dat       *hc_wrapper,
         struct sc_mb_dat           *sc_wrapper);


PRIVATE INLINE int
ml_pair3(vrna_fold_compound_t       *fc,
         int                        i,
         int                        j,
         int                        *dmli1,
         vrna_callback_hc_evaluate  *evaluate,
         struct hc_mb_def_dat       *hc_wrapper,
         struct sc_mb_dat           *sc_wrapper);


PRIVATE INLINE int
ml_pair53(vrna_fold_compound_t      *fc,
          int                       i,
          int                       j,
          int                       *dmli1,
          int                       *dmli2,
          vrna_callback_hc_evaluate *evaluate,
          struct hc_mb_def_dat      *hc_wrapper,
          struct sc_mb_dat          *sc_wrapper);


/*
 #################################
 # BEGIN OF FUNCTION DEFINITIONS #
 #################################
 */
PUBLIC int
vrna_E_mb_loop_fast(vrna_fold_compound_t  *fc,
                    int                   i,
                    int                   j,
                    int                   *dmli1,
                    int                   *dmli2)
{
  int e = INF;

  if (fc)
    e = E_mb_loop_fast(fc, i, j, dmli1, dmli2);

  return e;
}


PUBLIC int
vrna_E_ml_stems_fast(vrna_fold_compound_t *fc,
                     int                  i,
                     int                  j,
                     int                  *fmi,
                     int                  *dmli)
{
  int e = INF;

  if (fc)
    e = E_ml_stems_fast(fc, i, j, fmi, dmli);

  return e;
}


PUBLIC int
vrna_E_mb_loop_stack(vrna_fold_compound_t *fc,
                     int                  i,
                     int                  j)
{
  int e = INF;

  if (fc)
    e = E_mb_loop_stack(fc, i, j);

  return e;
}


PUBLIC int
E_ml_rightmost_stem(int                   i,
                    int                   j,
                    vrna_fold_compound_t  *fc)
{
  int e;

  e = INF;

  if ((fc) && (fc->matrices) && (fc->matrices->fM1)) {
    struct hc_mb_def_dat      hc_dat_local;
    struct sc_mb_dat          sc_wrapper;
    vrna_callback_hc_evaluate *evaluate;

    evaluate = prepare_hc_mb_def(fc, &hc_dat_local);
    init_sc_mb(fc, &sc_wrapper);

    e = extend_fm_3p(i, j, fc->matrices->fM1, fc, evaluate, &hc_dat_local, &sc_wrapper);

    if ((fc->aux_grammar) && (fc->aux_grammar->cb_aux_m1)) {
      int ee = fc->aux_grammar->cb_aux_m1(fc, i, j, fc->aux_grammar->data);
      e = MIN2(e, ee);
    }

    free_sc_mb(&sc_wrapper);
  }

  return e;
}


/*
 #####################################
 # BEGIN OF STATIC HELPER FUNCTIONS  #
 #####################################
 */
PRIVATE INLINE int
ml_pair_d0(vrna_fold_compound_t       *fc,
           int                        i,
           int                        j,
           int                        *dmli1,
           vrna_callback_hc_evaluate  *evaluate,
           struct hc_mb_def_dat       *hc_wrapper,
           struct sc_mb_dat           *sc_wrapper)
{
  short         *S, **SS;
  unsigned int  tt, s, n_seq;
  int           e;
  vrna_param_t  *P;
  vrna_md_t     *md;

  e = INF;

  if (evaluate(i, j, i + 1, j - 1, VRNA_DECOMP_PAIR_ML, hc_wrapper)) {
    e = dmli1[j - 1];

    if (e != INF) {
      P   = fc->params;
      md  = &(P->model_details);

      switch (fc->type) {
        case VRNA_FC_TYPE_SINGLE:
          S   = fc->sequence_encoding2;
          tt  = vrna_get_ptype_md(S[j], S[i], md);

          if (md->noGUclosure && ((tt == 3) || (tt == 4)))
            return INF; /* not allowed */

          e += E_MLstem(tt, -1, -1, P) +
               P->MLclosing;
          break;

        case VRNA_FC_TYPE_COMPARATIVE:
          n_seq = fc->n_seq;
          SS    = fc->S;
          for (s = 0; s < n_seq; s++) {
            tt  = vrna_get_ptype_md(SS[s][j], SS[s][i], md);
            e   += E_MLstem(tt, -1, -1, P);
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
ml_pair_d1(vrna_fold_compound_t       *fc,
           int                        i,
           int                        j,
           int                        *dmli1,
           int                        *dmli2,
           vrna_callback_hc_evaluate  *evaluate,
           struct hc_mb_def_dat       *hc_wrapper,
           struct sc_mb_dat           *sc_wrapper)
{
  int e, en;

  /* new closing pair (i,j) with mb part [i+1,j-1] */
  e = ml_pair_d0(fc, i, j, dmli1, evaluate, hc_wrapper, sc_wrapper);

  /* new closing pair (i,j) with mb part [i+2,j-1] */
  en  = ml_pair5(fc, i, j, dmli2, evaluate, hc_wrapper, sc_wrapper);
  e   = MIN2(e, en);

  /* new closing pair (i,j) with mb part [i+1, j-2] */
  en  = ml_pair3(fc, i, j, dmli1, evaluate, hc_wrapper, sc_wrapper);
  e   = MIN2(e, en);

  /* new closing pair (i,j) with mb part [i+2.j-2] */
  en  = ml_pair53(fc, i, j, dmli1, dmli2, evaluate, hc_wrapper, sc_wrapper);
  e   = MIN2(e, en);

  return e;
}


PRIVATE INLINE int
ml_pair_d2(vrna_fold_compound_t       *fc,
           int                        i,
           int                        j,
           int                        *dmli1,
           vrna_callback_hc_evaluate  *evaluate,
           struct hc_mb_def_dat       *hc_wrapper,
           struct sc_mb_dat           *sc_wrapper)
{
  short         *S, *S2, **SS, **S5, **S3, si1, sj1;
  unsigned int  tt, strands, *sn, s, n_seq;
  int           e;
  vrna_param_t  *P;
  vrna_md_t     *md;

  e = INF;

  if (evaluate(i, j, i + 1, j - 1, VRNA_DECOMP_PAIR_ML, hc_wrapper)) {
    e = dmli1[j - 1];

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

          e += E_MLstem(tt, sj1, si1, P) +
               P->MLclosing;
          break;

        case VRNA_FC_TYPE_COMPARATIVE:
          n_seq = fc->n_seq;
          SS    = fc->S;
          S5    = fc->S5;
          S3    = fc->S3;

          for (s = 0; s < n_seq; s++) {
            tt  = vrna_get_ptype_md(SS[s][j], SS[s][i], md);
            e   += E_MLstem(tt, S5[s][j], S3[s][i], P);
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
ml_pair5(vrna_fold_compound_t       *fc,
         int                        i,
         int                        j,
         int                        *dmli2,
         vrna_callback_hc_evaluate  *evaluate,
         struct hc_mb_def_dat       *hc_wrapper,
         struct sc_mb_dat           *sc_wrapper)
{
  short         *S, *S2, **SS, **S3, si1;
  unsigned int  tt, strands, *sn, n_seq, s;
  int           e;
  vrna_param_t  *P;
  vrna_md_t     *md;

  e = INF;

  if (evaluate(i, j, i + 2, j - 1, VRNA_DECOMP_PAIR_ML, hc_wrapper)) {
    e = dmli2[j - 1];

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

          e += E_MLstem(tt, -1, si1, P) +
               P->MLclosing +
               P->MLbase;
          break;

        case VRNA_FC_TYPE_COMPARATIVE:
          n_seq = fc->n_seq;
          SS    = fc->S;
          S3    = fc->S3;

          for (s = 0; s < n_seq; s++) {
            tt  = vrna_get_ptype_md(SS[s][j], SS[s][i], md);
            e   += E_MLstem(tt, -1, S3[s][i], P);
          }

          e += (P->MLclosing + P->MLbase) *
               n_seq;
          break;
      }

      if (sc_wrapper->pair5)
        e += sc_wrapper->pair5(i, j, sc_wrapper);
    }
  }

  return e;
}


PRIVATE INLINE int
ml_pair3(vrna_fold_compound_t       *fc,
         int                        i,
         int                        j,
         int                        *dmli1,
         vrna_callback_hc_evaluate  *evaluate,
         struct hc_mb_def_dat       *hc_wrapper,
         struct sc_mb_dat           *sc_wrapper)
{
  short         *S, *S2, **SS, **S5, sj1;
  unsigned int  tt, strands, *sn, n_seq, s;
  int           e;
  vrna_param_t  *P;
  vrna_md_t     *md;

  e = INF;

  if (evaluate(i, j, i + 1, j - 2, VRNA_DECOMP_PAIR_ML, hc_wrapper)) {
    e = dmli1[j - 2];

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

          sj1 = ((strands == 1) || (sn[j - 1] == sn[j])) ? S[j - 1] : -1;

          e += E_MLstem(tt, sj1, -1, P) +
               P->MLclosing +
               P->MLbase;
          break;

        case VRNA_FC_TYPE_COMPARATIVE:
          n_seq = fc->n_seq;
          SS    = fc->S;
          S5    = fc->S5;

          for (s = 0; s < n_seq; s++) {
            tt  = vrna_get_ptype_md(SS[s][j], SS[s][i], md);
            e   += E_MLstem(tt, S5[s][j], -1, P);
          }

          e += (P->MLclosing + P->MLbase) *
               n_seq;
          break;
      }

      if (sc_wrapper->pair3)
        e += sc_wrapper->pair3(i, j, sc_wrapper);
    }
  }

  return e;
}


PRIVATE INLINE int
ml_pair53(vrna_fold_compound_t      *fc,
          int                       i,
          int                       j,
          int                       *dmli1,
          int                       *dmli2,
          vrna_callback_hc_evaluate *evaluate,
          struct hc_mb_def_dat      *hc_wrapper,
          struct sc_mb_dat          *sc_wrapper)
{
  short         *S, *S2, **SS, **S3, **S5, si1, sj1;
  unsigned int  tt, strands, *sn, n_seq, s;
  int           e;
  vrna_param_t  *P;
  vrna_md_t     *md;

  e = INF;

  if (evaluate(i, j, i + 2, j - 2, VRNA_DECOMP_PAIR_ML, hc_wrapper)) {
    e = dmli2[j - 2];

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

          e += E_MLstem(tt, sj1, si1, P) +
               P->MLclosing +
               2 * P->MLbase;
          break;

        case VRNA_FC_TYPE_COMPARATIVE:
          n_seq = fc->n_seq;
          SS    = fc->S;
          S5    = fc->S5;
          S3    = fc->S3;

          for (s = 0; s < n_seq; s++) {
            tt  = vrna_get_ptype_md(SS[s][j], SS[s][i], md);
            e   += E_MLstem(tt, S5[s][j], S3[s][i], P);
          }

          e += (P->MLclosing + 2 * P->MLbase) *
               n_seq;
          break;
      }

      if (sc_wrapper->pair53)
        e += sc_wrapper->pair53(i, j, sc_wrapper);
    }
  }

  return e;
}


PRIVATE int
E_mb_loop_fake(vrna_fold_compound_t *fc,
               int                  i,
               int                  j)
{
  short                     S_i1, S_j1, *S, *S2;
  unsigned int              strands, *sn;
  int                       decomp, en, e, *fC, dangle_model, tt, noGUclosure;
  vrna_param_t              *P;
  vrna_md_t                 *md;
  vrna_callback_hc_evaluate *evaluate;
  struct hc_mb_def_dat      hc_dat_local;

  S             = fc->sequence_encoding;
  S2            = fc->sequence_encoding2;
  strands       = fc->strands;
  sn            = fc->strand_number;
  fC            = fc->matrices->fc;
  P             = fc->params;
  md            = &(P->model_details);
  noGUclosure   = md->noGUclosure;
  dangle_model  = md->dangles;

  /* init values */
  e       = INF;
  decomp  = INF;

  evaluate = prepare_hc_mb_def_ext(fc, &hc_dat_local);

  tt = vrna_get_ptype_md(S2[j], S2[i], md);

  if (noGUclosure && ((tt == 3) || (tt == 4)))
    return e;

  if (strands == 1) {
    S_i1  = S[i + 1];
    S_j1  = S[j - 1];
  } else {
    S_i1  = (sn[i] == sn[i + 1]) ? S[i + 1] : -1;
    S_j1  = (sn[j - 1] == sn[j]) ? S[j - 1] : -1;
  }

  /* multibranch like cofold structure with cut somewhere between i and j */
  if (evaluate(i, j, i, j, VRNA_DECOMP_EXT_STEM, &hc_dat_local)) {
    if ((fC[i + 1] != INF) && (fC[j - 1] != INF)) {
      decomp = fC[i + 1] +
               fC[j - 1];
      switch (dangle_model) {
        case 0:
          decomp += vrna_E_ext_stem(tt, -1, -1, P);
          break;

        case 2:
          decomp += vrna_E_ext_stem(tt, S_j1, S_i1, P);
          break;

        default:
          decomp += vrna_E_ext_stem(tt, -1, -1, P);
          break;
      }
    }
  }

  if (dangle_model % 2) {
    /* dangles == 1 || dangles == 3 */
    if (evaluate(i + 1, j - 1, i + 2, j - 1, VRNA_DECOMP_EXT_EXT, &hc_dat_local)) {
      if ((fC[i + 2] != INF) && (fC[j - 1] != INF)) {
        en = fC[i + 2] +
             fC[j - 1] +
             vrna_E_ext_stem(tt, -1, S_i1, P);
        decomp = MIN2(decomp, en);
      }
    }

    if (evaluate(i + 1, j - 1, i + 1, j - 2, VRNA_DECOMP_EXT_EXT, &hc_dat_local)) {
      if ((fC[i + 1] != INF) && (fC[j - 2] != INF)) {
        en = fC[i + 1] +
             fC[j - 2] +
             vrna_E_ext_stem(tt, S_j1, -1, P);
        decomp = MIN2(decomp, en);
      }
    }

    if (evaluate(i + 1, j - 1, i + 2, j - 2, VRNA_DECOMP_EXT_EXT, &hc_dat_local)) {
      if ((fC[i + 2] != INF) && (fC[j - 2] != INF)) {
        en = fC[i + 2] +
             fC[j - 2] +
             vrna_E_ext_stem(tt, S_j1, S_i1, P);
        decomp = MIN2(decomp, en);
      }
    }
  }

  e = MIN2(e, decomp);

  return e;
}


PRIVATE int
E_mb_loop_fast(vrna_fold_compound_t *fc,
               int                  i,
               int                  j,
               int                  *dmli1,
               int                  *dmli2)
{
  unsigned int              *sn;
  int                       decomp, e, dangle_model;
  vrna_param_t              *P;
  vrna_md_t                 *md;
  vrna_callback_hc_evaluate *evaluate;
  struct hc_mb_def_dat      hc_dat_local;
  struct sc_mb_dat          sc_wrapper;

  sn            = fc->strand_number;
  P             = fc->params;
  md            = &(P->model_details);
  dangle_model  = md->dangles;

  /* init values */
  e       = INF;
  decomp  = INF;

  evaluate = prepare_hc_mb_def(fc, &hc_dat_local);

  init_sc_mb(fc, &sc_wrapper);

  /* do pointer magic for sliding window implementation */
  if (fc->hc->type == VRNA_HC_WINDOW) {
    dmli1 -= i + 1;
    if (dmli2)
      dmli2 -= i + 2;
  }

  switch (dangle_model) {
    /* no dangles */
    case 0:
      decomp = ml_pair_d0(fc, i, j, dmli1, evaluate, &hc_dat_local, &sc_wrapper);
      break;

    /* double dangles */
    case 2:
      decomp = ml_pair_d2(fc, i, j, dmli1, evaluate, &hc_dat_local, &sc_wrapper);
      break;

    /* normal dangles, aka dangles = 1 || 3 */
    default:
      decomp = ml_pair_d1(fc, i, j, dmli1, dmli2, evaluate, &hc_dat_local, &sc_wrapper);
      break;
  }

  free_sc_mb(&sc_wrapper);

  e = MIN2(e, decomp);

  /* add additional cases for possible strand nicks between i and j */
  if ((fc->type == VRNA_FC_TYPE_SINGLE) && (sn[i] != sn[j])) {
    decomp  = E_mb_loop_fake(fc, i, j);
    e       = MIN2(e, decomp);
  }

  return e;
}


PRIVATE int
E_mb_loop_stack(vrna_fold_compound_t  *fc,
                int                   i,
                int                   j)
{
  char                      *ptype, **ptype_local;
  short                     **SS;
  unsigned int              n_seq, s, *tt, sliding_window;
  int                       *c, *fML, e, decomp, en, i1k, k1j1, ij, k, *indx, turn,
                            type, type_2, *rtype, **c_local, **fML_local;
  vrna_param_t              *P;
  vrna_md_t                 *md;
  vrna_callback_hc_evaluate *evaluate;
  struct hc_mb_def_dat      hc_dat_local;
  struct sc_mb_dat          sc_wrapper;

  sliding_window = (fc->hc->type == VRNA_HC_WINDOW) ? 1 : 0;

  n_seq       = (fc->type == VRNA_FC_TYPE_SINGLE) ? 1 : fc->n_seq;
  SS          = (fc->type == VRNA_FC_TYPE_SINGLE) ? NULL : fc->S;
  indx        = fc->jindx;
  P           = fc->params;
  md          = &(P->model_details);
  turn        = md->min_loop_size;
  c           = (sliding_window) ? NULL : fc->matrices->c;
  fML         = (sliding_window) ? NULL : fc->matrices->fML;
  c_local     = (sliding_window) ? fc->matrices->c_local : NULL;
  fML_local   = (sliding_window) ? fc->matrices->fML_local : NULL;
  ptype       = (fc->type == VRNA_FC_TYPE_SINGLE) ? (sliding_window ? NULL : fc->ptype) : NULL;
  ptype_local =
    (fc->type == VRNA_FC_TYPE_SINGLE) ? (sliding_window ? fc->ptype_local : NULL) : NULL;
  rtype = (fc->type == VRNA_FC_TYPE_SINGLE) ? &(md->rtype[0]) : NULL;
  tt    = NULL;
  type  = 0;
  ij    = (sliding_window) ? 0 : indx[j] + i;
  e     = INF;

  evaluate = prepare_hc_mb_def(fc, &hc_dat_local);
  init_sc_mb(fc, &sc_wrapper);

  /* prepare type(s) for enclosing pair (i, j) */
  if (fc->type == VRNA_FC_TYPE_COMPARATIVE) {
    tt = (unsigned int *)vrna_alloc(sizeof(unsigned int) * n_seq);
    for (s = 0; s < n_seq; s++)
      tt[s] = vrna_get_ptype_md(SS[s][i], SS[s][j], md);
  } else if (sliding_window) {
    type = vrna_get_ptype_window(i, j, ptype_local);
  } else {
    type = vrna_get_ptype(ij, ptype);
  }

  if (evaluate(i, j, i + 1, j - 1, VRNA_DECOMP_PAIR_ML, &hc_dat_local)) {
    decomp = INF;
    if (sliding_window) {
      for (k = i + 2 + turn; k < j - 2 - turn; k++) {
        if (evaluate(i, j, i + 1, k, VRNA_DECOMP_ML_COAXIAL, &hc_dat_local)) {
          en = c_local[i + 1][k - i - 1] +
               fML_local[k + 1][j - 1 - k - 1];

          switch (fc->type) {
            case VRNA_FC_TYPE_SINGLE:
              type_2  = rtype[vrna_get_ptype_window(i + 1, k, ptype_local)];
              en      += P->stack[type][type_2];
              break;

            case VRNA_FC_TYPE_COMPARATIVE:
              for (s = 0; s < n_seq; s++) {
                type_2  = vrna_get_ptype_md(SS[s][k], SS[s][i + 1], md);
                en      += P->stack[tt[s]][type_2];
              }
              break;
          }

          if (sc_wrapper.coaxial_cls)
            en += sc_wrapper.coaxial_cls(i, j, i + 1, k, &sc_wrapper);

          decomp = MIN2(decomp, en);
        }

        if (evaluate(i, j, k + 1, j - 1, VRNA_DECOMP_ML_COAXIAL, &hc_dat_local)) {
          en = c_local[k + 1][j - 1 - k - 1] +
               fML_local[i + 1][k - i - 1];

          switch (fc->type) {
            case VRNA_FC_TYPE_SINGLE:
              type_2  = rtype[vrna_get_ptype_window(k + 1, j - 1, ptype_local)];
              en      += P->stack[type][type_2];

              break;

            case VRNA_FC_TYPE_COMPARATIVE:
              for (s = 0; s < n_seq; s++) {
                type_2  = vrna_get_ptype_md(SS[s][j - 1], SS[s][k + 1], md);
                en      += P->stack[tt[s]][type_2];
              }

              break;
          }

          if (sc_wrapper.coaxial_cls)
            en += sc_wrapper.coaxial_cls(i, j, k + 1, j - 1, &sc_wrapper);

          decomp = MIN2(decomp, en);
        }
      }
    } else {
      k1j1 = indx[j - 1] + i + 2 + turn + 1;
      for (k = i + 2 + turn; k < j - 2 - turn; k++, k1j1++) {
        i1k = indx[k] + i + 1;

        if (evaluate(i, j, i + 1, k, VRNA_DECOMP_ML_COAXIAL, &hc_dat_local)) {
          en = c[i1k] +
               fML[k1j1];

          switch (fc->type) {
            case VRNA_FC_TYPE_SINGLE:
              type_2  = rtype[vrna_get_ptype(i1k, ptype)];
              en      += P->stack[type][type_2];

              break;

            case VRNA_FC_TYPE_COMPARATIVE:
              for (s = 0; s < n_seq; s++) {
                type_2  = vrna_get_ptype_md(SS[s][k], SS[s][i + 1], md);
                en      += P->stack[tt[s]][type_2];
              }

              break;
          }

          if (sc_wrapper.coaxial_cls)
            en += sc_wrapper.coaxial_cls(i, j, i + 1, k, &sc_wrapper);

          decomp = MIN2(decomp, en);
        }

        if (evaluate(i, j, k + 1, j - 1, VRNA_DECOMP_ML_COAXIAL, &hc_dat_local)) {
          en = c[k1j1] +
               fML[i1k];

          switch (fc->type) {
            case VRNA_FC_TYPE_SINGLE:
              type_2  = rtype[vrna_get_ptype(k1j1, ptype)];
              en      += P->stack[type][type_2];

              break;

            case VRNA_FC_TYPE_COMPARATIVE:
              for (s = 0; s < n_seq; s++) {
                type_2  = vrna_get_ptype_md(SS[s][j - 1], SS[s][k + 1], md);
                en      += P->stack[tt[s]][type_2];
              }

              break;
          }

          if (sc_wrapper.coaxial_cls)
            en += sc_wrapper.coaxial_cls(i, j, k + 1, j - 1, &sc_wrapper);

          decomp = MIN2(decomp, en);
        }
      }
    }

    /* no TermAU penalty if coax stack */
    decomp += (2 * P->MLintern[1] + P->MLclosing) *
              n_seq;

    if (sc_wrapper.pair)
      decomp += sc_wrapper.pair(i, j, &sc_wrapper);

    e = decomp;
  }

  free_sc_mb(&sc_wrapper);
  free(tt);

  return e;
}


/*
 * compose a multibranch loop part fm[i:j]
 * by either c[i,j]/ggg[i,j] or fm[i:j-1]
 *
 * This function can be used for fM and fM1
 */
PRIVATE int
extend_fm_3p(int                        i,
             int                        j,
             int                        *fm,
             vrna_fold_compound_t       *fc,
             vrna_callback_hc_evaluate  *evaluate,
             struct hc_mb_def_dat       *hc_dat_local,
             struct sc_mb_dat           *sc_wrapper)
{
  short         *S, **SS, **S5, **S3;
  unsigned int  *sn, n_seq, s, sliding_window;
  int           en, en2, length, *indx, *c, **c_local, **fm_local, *ggg, **ggg_local, ij, type,
                dangle_model, with_gquad, e, u, k, cnt, with_ud;
  vrna_param_t  *P;
  vrna_md_t     *md;
  vrna_ud_t     *domains_up;

  sliding_window  = (fc->hc->type == VRNA_HC_WINDOW) ? 1 : 0;
  n_seq           = (fc->type == VRNA_FC_TYPE_SINGLE) ? 1 : fc->n_seq;
  length          = fc->length;
  S               = (fc->type == VRNA_FC_TYPE_SINGLE) ? fc->sequence_encoding : NULL;
  SS              = (fc->type == VRNA_FC_TYPE_SINGLE) ? NULL : fc->S;
  S5              = (fc->type == VRNA_FC_TYPE_SINGLE) ? NULL : fc->S5;
  S3              = (fc->type == VRNA_FC_TYPE_SINGLE) ? NULL : fc->S3;
  indx            = (sliding_window) ? NULL : fc->jindx;
  sn              = fc->strand_number;
  c               = (sliding_window) ? NULL : fc->matrices->c;
  ggg             = (sliding_window) ? NULL : fc->matrices->ggg;
  c_local         = (sliding_window) ? fc->matrices->c_local : NULL;
  fm_local        = (sliding_window) ? fc->matrices->fML_local : NULL;
  ggg_local       = (sliding_window) ? fc->matrices->ggg_local : NULL;
  ij              = (sliding_window) ? 0 : indx[j] + i;
  P               = fc->params;
  md              = &(P->model_details);
  dangle_model    = md->dangles;
  with_gquad      = md->gquad;
  domains_up      = fc->domains_up;
  with_ud         = (domains_up && domains_up->energy_cb) ? 1 : 0;
  e               = INF;

  if (fm == NULL) {
    if (sliding_window)
      fm_local = fc->matrices->fML_local;
    else
      fm = fc->matrices->fML;
  }

  if (evaluate(i, j, i, j, VRNA_DECOMP_ML_STEM, hc_dat_local)) {
    en = (sliding_window) ? c_local[i][j - i] : c[ij];
    if (en != INF) {
      switch (fc->type) {
        case VRNA_FC_TYPE_SINGLE:
          type = (sliding_window) ? vrna_get_ptype_window(i, j, fc->ptype_local) : vrna_get_ptype(
            ij,
            fc->ptype);
          if (dangle_model == 2)
            en += E_MLstem(type, (i == 1) ? S[length] : S[i - 1], S[j + 1], P);
          else
            en += E_MLstem(type, -1, -1, P);

          break;

        case VRNA_FC_TYPE_COMPARATIVE:
          if (dangle_model == 2) {
            for (s = 0; s < n_seq; s++) {
              type  = vrna_get_ptype_md(SS[s][i], SS[s][j], md);
              en    += E_MLstem(type, S5[s][i], S3[s][j], P);
            }
          } else {
            for (s = 0; s < n_seq; s++) {
              type  = vrna_get_ptype_md(SS[s][i], SS[s][j], md);
              en    += E_MLstem(type, -1, -1, P);
            }
          }

          break;
      }

      if (sc_wrapper->red_stem)
        en += sc_wrapper->red_stem(i, j, i, j, sc_wrapper);

      e = MIN2(e, en);
    }
  }

  if (with_gquad) {
    if (sn[i] == sn[j]) {
      en  = (sliding_window) ? ggg_local[i][j - i] : ggg[ij];
      en  += E_MLstem(0, -1, -1, P) *
             n_seq;

      e = MIN2(e, en);
    }
  }

  if (evaluate(i, j, i, j - 1, VRNA_DECOMP_ML_ML, hc_dat_local)) {
    en = (sliding_window) ? fm_local[i][j - 1 - i] : fm[indx[j - 1] + i];
    if (en != INF) {
      en += P->MLbase *
            n_seq;

      if (sc_wrapper->red_ml)
        en += sc_wrapper->red_ml(i, j, i, j - 1, sc_wrapper);

      e = MIN2(e, en);
    }
  }

  if (with_ud) {
    for (cnt = 0; cnt < domains_up->uniq_motif_count; cnt++) {
      u = domains_up->uniq_motif_size[cnt];
      k = j - u + 1;
      if ((k > i) && (evaluate(i, j, i, k - 1, VRNA_DECOMP_ML_ML, hc_dat_local))) {
        en = (sliding_window) ? fm_local[i][k - 1 - i] : fm[indx[k - 1] + i];
        if (en != INF) {
          en += u * P->MLbase *
                n_seq;

          en2 = domains_up->energy_cb(fc,
                                      k,
                                      j,
                                      VRNA_UNSTRUCTURED_DOMAIN_MB_LOOP | VRNA_UNSTRUCTURED_DOMAIN_MOTIF,
                                      domains_up->data);
          if (en2 != INF) {
            en += en2;

            if (sc_wrapper->red_ml)
              en += sc_wrapper->red_ml(i, j, i, k - 1, sc_wrapper);

            e = MIN2(e, en);
          }
        }
      }
    }
  }

  return e;
}


PRIVATE int
E_ml_stems_fast(vrna_fold_compound_t  *fc,
                int                   i,
                int                   j,
                int                   *fmi,
                int                   *dmli)
{
  char                      *ptype, **ptype_local;
  short                     *S, **SS, **S5, **S3;
  unsigned int              *sn, *se, n_seq, s;
  int                       k, en, decomp, mm5, mm3, type_2, k1j, length, *indx,
                            *c, *fm, ij, dangle_model, turn, type, *rtype, circular, e, u,
                            cnt, with_ud, sliding_window, **c_local, **fm_local;
  vrna_hc_t                 *hc;
  vrna_sc_t                 *sc;
  vrna_param_t              *P;
  vrna_md_t                 *md;
  vrna_ud_t                 *domains_up;
  vrna_callback_hc_evaluate *evaluate;
  struct hc_mb_def_dat      hc_dat_local;
  struct sc_mb_dat          sc_wrapper;

  sliding_window = (fc->hc->type == VRNA_HC_WINDOW) ? 1 : 0;

  length      = (int)fc->length;
  ptype       = (fc->type == VRNA_FC_TYPE_SINGLE) ? ((sliding_window) ? NULL : fc->ptype) : NULL;
  ptype_local =
    (fc->type == VRNA_FC_TYPE_SINGLE) ? ((sliding_window) ? fc->ptype_local : NULL) : NULL;
  S             = (fc->type == VRNA_FC_TYPE_SINGLE) ? fc->sequence_encoding : NULL;
  SS            = (fc->type == VRNA_FC_TYPE_SINGLE) ? NULL : fc->S;
  S5            = (fc->type == VRNA_FC_TYPE_SINGLE) ? NULL : fc->S5;
  S3            = (fc->type == VRNA_FC_TYPE_SINGLE) ? NULL : fc->S3;
  indx          = (sliding_window) ? NULL : fc->jindx;
  n_seq         = (fc->type == VRNA_FC_TYPE_SINGLE) ? 1 : fc->n_seq;
  sn            = fc->strand_number;
  se            = fc->strand_end;
  hc            = fc->hc;
  sc            = fc->sc;
  c             = (sliding_window) ? NULL : fc->matrices->c;
  c_local       = (sliding_window) ? fc->matrices->c_local : NULL;
  fm            = (sliding_window) ? NULL : fc->matrices->fML;
  fm_local      = (sliding_window) ? fc->matrices->fML_local : NULL;
  P             = fc->params;
  md            = &(P->model_details);
  ij            = (sliding_window) ? 0 : indx[j] + i;
  dangle_model  = md->dangles;
  turn          = md->min_loop_size;
  rtype         = &(md->rtype[0]);
  circular      = md->circ;
  domains_up    = fc->domains_up;
  with_ud       = (domains_up && domains_up->energy_cb) ? 1 : 0;
  e             = INF;
  evaluate      = prepare_hc_mb_def(fc, &hc_dat_local);

  init_sc_mb(fc, &sc_wrapper);

  /*
   *  extension with one unpaired nucleotide at the right (3' site)
   *  or full branch of (i,j)
   */
  e = extend_fm_3p(i, j, NULL, fc, evaluate, &hc_dat_local, &sc_wrapper);

  /*
   *  extension with one unpaired nucleotide at 5' site
   *  and all other variants which are needed for odd
   *  dangle models
   */
  if (evaluate(i, j, i + 1, j, VRNA_DECOMP_ML_ML, &hc_dat_local)) {
    en = (sliding_window) ? fm_local[i + 1][j - i - 1] : fm[ij + 1];
    if (en != INF) {
      en += P->MLbase *
            n_seq;

      if (sc_wrapper.red_ml)
        en += sc_wrapper.red_ml(i, j, i + 1, j, &sc_wrapper);

      e = MIN2(e, en);
    }
  }

  /* extension with bound ligand on 5'site */
  if (with_ud) {
    for (cnt = 0; cnt < domains_up->uniq_motif_count; cnt++) {
      u = domains_up->uniq_motif_size[cnt];
      k = i + u - 1;
      if ((k < j) && (evaluate(i, j, k + 1, j, VRNA_DECOMP_ML_ML, &hc_dat_local))) {
        decomp = (sliding_window) ? fm_local[i + u][j - (i + u)] : fm[ij + u];
        if (decomp != INF) {
          decomp += u * P->MLbase *
                    n_seq;

          en = domains_up->energy_cb(fc,
                                     i,
                                     k,
                                     VRNA_UNSTRUCTURED_DOMAIN_MB_LOOP | VRNA_UNSTRUCTURED_DOMAIN_MOTIF,
                                     domains_up->data);
          if (en != INF) {
            decomp += en;

            if (sc_wrapper.red_ml)
              decomp += sc_wrapper.red_ml(i, j, k + 1, j, &sc_wrapper);

            e = MIN2(e, decomp);
          }
        }
      }
    }
  }

  if (dangle_model % 2) {
    /* dangle_model = 1 || 3 */
    mm5 = mm3 = -1;

    if (fc->type == VRNA_FC_TYPE_SINGLE) {
      if ((i > 1) || circular)
        mm5 = S[i];

      if ((j < length) || circular)
        mm3 = S[j];
    }

    if (evaluate(i, j, i + 1, j, VRNA_DECOMP_ML_STEM, &hc_dat_local)) {
      en = (sliding_window) ? c_local[i + 1][j - (i + 1)] : c[ij + 1];
      if (en != INF) {
        en += P->MLbase *
              n_seq;

        switch (fc->type) {
          case VRNA_FC_TYPE_SINGLE:
            type =
              (sliding_window) ? vrna_get_ptype_window(i + 1, j, ptype_local) : vrna_get_ptype(
                ij + 1,
                ptype);
            en += E_MLstem(type, mm5, -1, P);
            break;

          case VRNA_FC_TYPE_COMPARATIVE:
            for (s = 0; s < n_seq; s++) {
              type  = vrna_get_ptype_md(SS[s][i + 1], SS[s][j], md);
              en    += E_MLstem(type, S5[s][i + 1], -1, P);
            }
            break;
        }

        if (sc_wrapper.red_ml)
          en += sc_wrapper.red_ml(i, j, i + 1, j, &sc_wrapper);

        e = MIN2(e, en);
      }
    }

    if (evaluate(i, j, i, j - 1, VRNA_DECOMP_ML_STEM, &hc_dat_local)) {
      en = (sliding_window) ? c_local[i][j - 1 - i] : c[indx[j - 1] + i];
      if (en != INF) {
        en += P->MLbase *
              n_seq;

        switch (fc->type) {
          case VRNA_FC_TYPE_SINGLE:
            type =
              (sliding_window) ? vrna_get_ptype_window(i, j - 1,
                                                       ptype_local) : vrna_get_ptype(
                indx[j - 1] + i,
                ptype);
            en += E_MLstem(type, -1, mm3, P);
            break;

          case VRNA_FC_TYPE_COMPARATIVE:
            for (s = 0; s < n_seq; s++) {
              type  = vrna_get_ptype_md(SS[s][i], SS[s][j - 1], md);
              en    += E_MLstem(type, -1, S3[s][j - 1], P);
            }
            break;
        }

        if (sc_wrapper.red_ml)
          en += sc_wrapper.red_ml(i, j, i, j - 1, &sc_wrapper);

        e = MIN2(e, en);
      }
    }

    if (evaluate(i, j, i + 1, j - 1, VRNA_DECOMP_ML_STEM, &hc_dat_local)) {
      en = (sliding_window) ? c_local[i + 1][j - 1 - (i + 1)] : c[indx[j - 1] + i + 1];
      if (en != INF) {
        en += 2 * P->MLbase *
              n_seq;

        switch (fc->type) {
          case VRNA_FC_TYPE_SINGLE:
            type =
              (sliding_window) ? vrna_get_ptype_window(i + 1, j - 1,
                                                       ptype_local) : vrna_get_ptype(
                indx[j - 1] + i + 1,
                ptype);
            en += E_MLstem(type, mm5, mm3, P);
            break;

          case VRNA_FC_TYPE_COMPARATIVE:
            for (s = 0; s < n_seq; s++) {
              type  = vrna_get_ptype_md(SS[s][i + 1], SS[s][j - 1], md);
              en    += E_MLstem(type, S5[s][i + 1], S3[s][j - 1], P);
            }
            break;
        }

        if (sc_wrapper.red_ml)
          en += sc_wrapper.red_ml(i, j, i + 1, j - 1, &sc_wrapper);

        e = MIN2(e, en);
      }
    } /* end special cases for dangles == 1 || dangles == 3 */
  }

  /* modular decomposition -------------------------------*/
  if (sliding_window) {
    fmi   -= i;
    dmli  -= i;
  }

  /* use fmi pointer that we may extend to include hard/soft constraints if necessary */
  int *fmi_tmp = fmi;

  if (hc->f) {
    fmi_tmp = (int *)vrna_alloc(sizeof(int) * (j - i + 2));
    fmi_tmp -= i;

    /* copy data */
    for (k = i + 1 + turn; k <= j - 2 - turn; k++)
      fmi_tmp[k] = fmi[k];

    /* mask unavailable decompositions */
    for (k = i + 1 + turn; k <= j - 2 - turn; k++)
      if (!hc->f(i, j, k, k + 1, VRNA_DECOMP_ML_ML_ML, hc->data))
        fmi_tmp[k] = INF;
  }

  if (sc_wrapper.decomp_ml) {
    if (fmi_tmp == fmi) {
      fmi_tmp = (int *)vrna_alloc(sizeof(int) * (j - i + 2));
      fmi_tmp -= i;

      /* copy data */
      for (k = i + 1 + turn; k <= j - 2 - turn; k++)
        fmi_tmp[k] = fmi[k];
    }

    for (k = i + 1 + turn; k <= j - 2 - turn; k++)
      if (fmi_tmp[k] != INF)
        fmi_tmp[k] += sc_wrapper.decomp_ml(i, j, k, k + 1, &sc_wrapper);
  }

  /* modular decomposition -------------------------------*/
  if (sliding_window) {
    for (decomp = INF, k = i + 1 + turn; k <= j - 2 - turn; k++) {
      if ((fmi_tmp[k] != INF) && (fm_local[k + 1][j - (k + 1)] != INF)) {
        en      = fmi_tmp[k] + fm_local[k + 1][j - (k + 1)];
        decomp  = MIN2(decomp, en);
      }
    }
  } else {
    decomp  = INF;
    k       = i + turn + 1;
    if (k >= j)
      k = j - 1;

    k1j = indx[j] + k + 1;

    /*
     *  loop over entire range but skip decompositions with in-between strand nick,
     *  this should be faster than evaluating hard constraints callback for each
     *  decomposition
     */
    while (1) {
      int last_nt = se[sn[k]];  /* go to last nucleotide of current strand */
      if (last_nt > j - turn - 2)
        last_nt = j - turn - 2;     /* at most go to last possible decomposition split before reaching j */

      if (last_nt < i)
        last_nt = i; /* do not start before i */

      const int count = last_nt - k;

      en      = vrna_fun_zip_add_min(fmi_tmp + k, fm + k1j, count);
      decomp  = MIN2(decomp, en);

      /* advance counters by processed subsegment and add 1 for the split point between strands */
      k   += count + 1;
      k1j += count + 1;

      if (k > j - turn - 2)
        break;
    }
  }

  /* end modular decomposition -------------------------------*/

  if (fmi_tmp != fmi) {
    fmi_tmp += i;
    free(fmi_tmp);
  }

  dmli[j] = decomp;               /* store for use in fast ML decompositon */

  e = MIN2(e, decomp);

  /* coaxial stacking */
  if (dangle_model == 3) {
    /* additional ML decomposition as two coaxially stacked helices */
    if (sliding_window) {
      /* additional ML decomposition as two coaxially stacked helices */
      for (decomp = INF, k = i + 1 + turn; k <= j - 2 - turn; k++) {
        if (evaluate(i, k, k + 1, j, VRNA_DECOMP_ML_COAXIAL_ENC, &hc_dat_local)) {
          type    = rtype[vrna_get_ptype_window(i, k, ptype_local)];
          type_2  = rtype[vrna_get_ptype_window(k + 1, j, ptype_local)];

          en = c_local[i][k - i] +
               c_local[k + 1][j - k - 1] +
               P->stack[type][type_2];

          if (sc)
            if (sc->f)
              en += sc->f(i, k, k + 1, j, VRNA_DECOMP_ML_COAXIAL_ENC, sc->data);

          decomp = MIN2(decomp, en);
        }
      }
    } else {
      int ik;
      k1j     = indx[j] + i + turn + 2;
      decomp  = INF;
      k       = i + turn + 1;
      if (k >= j)
        k = j - 1;

      /*
       *  loop over entire range but skip decompositions with in-between strand nick,
       *  this should be faster than evaluating hard constraints callback for each
       *  decomposition
       */
      while (1) {
        int last_nt = se[sn[k - 1]];  /* go to last nucleotide of current strand */
        if (last_nt > j - turn - 2)
          last_nt = j - turn - 2;     /* at most go to last possible decomposition split before reaching j */

        if (last_nt < i)
          last_nt = i; /* do not start before i */

        const int stop = last_nt;
        for (; k <= stop; k++, k1j++) {
          ik = indx[k] + i;
          if (evaluate(i, k, k + 1, j, VRNA_DECOMP_ML_COAXIAL_ENC, &hc_dat_local)) {
            en = c[ik] +
                 c[k1j];

            switch (fc->type) {
              case VRNA_FC_TYPE_SINGLE:
                type    = rtype[vrna_get_ptype(ik, ptype)];
                type_2  = rtype[vrna_get_ptype(k1j, ptype)];

                en += P->stack[type][type_2];

                break;

              case VRNA_FC_TYPE_COMPARATIVE:
                for (s = 0; s < n_seq; s++) {
                  type    = vrna_get_ptype_md(SS[s][k], SS[s][i], md);
                  type_2  = vrna_get_ptype_md(SS[s][j], SS[s][k + 1], md);

                  en += P->stack[type][type_2];
                }

                break;
            }

            if (sc_wrapper.coaxial_enc)
              en += sc_wrapper.coaxial_enc(i, k, k + 1, j, &sc_wrapper);

            decomp = MIN2(decomp, en);
          }
        }

        k++;
        k1j++;

        if (k > j - turn - 2)
          break;
      }
    }

    /* no TermAU penalty if coax stack */
    decomp += 2 * P->MLintern[1] *
              n_seq;
#if 0
    /*
     * This is needed for Y shaped ML loops with coax stacking of
     * interior pairts, but backtracking will fail if activated
     */
    DMLi[j] = MIN2(DMLi[j], decomp);
    DMLi[j] = MIN2(DMLi[j], DMLi[j - 1] + P->MLbase);
    DMLi[j] = MIN2(DMLi[j], DMLi1[j] + P->MLbase);
    new_fML = MIN2(new_fML, DMLi[j]);
#endif
    e = MIN2(e, decomp);
  }

  if ((fc->aux_grammar) && (fc->aux_grammar->cb_aux_m)) {
    en  = fc->aux_grammar->cb_aux_m(fc, i, j, fc->aux_grammar->data);
    e   = MIN2(e, en);
  }

  fmi[j] = e;

  free_sc_mb(&sc_wrapper);

  return e;
}
