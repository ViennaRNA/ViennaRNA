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
#include "ViennaRNA/sequences/alphabet.h"
#include "ViennaRNA/params/default.h"
#include "ViennaRNA/constraints/hard.h"
#include "ViennaRNA/constraints/soft.h"
#include "ViennaRNA/eval/exterior.h"
#include "ViennaRNA/eval/gquad.h"
#include "ViennaRNA/eval/multibranch.h"
#include "ViennaRNA/structured_domains.h"
#include "ViennaRNA/unstructured_domains.h"
#include "ViennaRNA/utils/higher_order_functions.h"

#include "ViennaRNA/intern/grammar_dat.h"

#ifdef __GNUC__
# define INLINE inline
#else
# define INLINE
#endif

#include "ViennaRNA/constraints/multibranch_hc.inc"
#include "ViennaRNA/constraints/multibranch_sc.inc"

#include "ViennaRNA/mfe/multibranch.h"

struct vrna_mx_mfe_aux_ml_s {
  unsigned int  n;
  int           *Fmi;
  int           *DMLi;
  int           *DMLi1;
  int           *DMLi2;
};


/*
 #################################
 # PRIVATE FUNCTION DECLARATIONS #
 #################################
 */

PRIVATE struct vrna_mx_mfe_aux_ml_s *
get_aux_arrays(unsigned int length);


PRIVATE void
rotate_aux_arrays(struct vrna_mx_mfe_aux_ml_s *aux);


PRIVATE void
free_aux_arrays(struct vrna_mx_mfe_aux_ml_s *aux);


PRIVATE int
mfe_multibranch_m2_fast(vrna_fold_compound_t        *fc,
                        unsigned int                i,
                        unsigned int                j,
                        struct vrna_mx_mfe_aux_ml_s *helpers);


PRIVATE int
E_mb_loop_fast(vrna_fold_compound_t         *fc,
               unsigned int                 i,
               unsigned int                 j,
               struct vrna_mx_mfe_aux_ml_s  *helpers);


PRIVATE int
E_ml_stems_fast(vrna_fold_compound_t        *fc,
                unsigned int                i,
                unsigned int                j,
                struct vrna_mx_mfe_aux_ml_s *helpers);


PRIVATE int
extend_fm_3p(unsigned int         i,
             unsigned int         j,
             int                  *fm,
             vrna_fold_compound_t *fc,
             vrna_hc_eval_f       evaluate,
             struct hc_mb_def_dat *hc_data,
             struct sc_mb_dat     *sc_wrapper);


PRIVATE int
E_mb_loop_stack(vrna_fold_compound_t  *fc,
                unsigned int          i,
                unsigned int          j);


PRIVATE INLINE int
ml_pair5(vrna_fold_compound_t *fc,
         unsigned int         i,
         unsigned int         j,
         int                  *dmli2,
         vrna_hc_eval_f       evaluate,
         struct hc_mb_def_dat *hc_wrapper,
         struct sc_mb_dat     *sc_wrapper);


PRIVATE INLINE int
ml_pair3(vrna_fold_compound_t *fc,
         unsigned int         i,
         unsigned int         j,
         int                  *dmli1,
         vrna_hc_eval_f       evaluate,
         struct hc_mb_def_dat *hc_wrapper,
         struct sc_mb_dat     *sc_wrapper);


PRIVATE INLINE int
ml_pair53(vrna_fold_compound_t  *fc,
          unsigned int          i,
          unsigned int          j,
          int                   *dmli1,
          int                   *dmli2,
          vrna_hc_eval_f        evaluate,
          struct hc_mb_def_dat  *hc_wrapper,
          struct sc_mb_dat      *sc_wrapper);


/*
 #################################
 # BEGIN OF FUNCTION DEFINITIONS #
 #################################
 */
PUBLIC struct vrna_mx_mfe_aux_ml_s *
vrna_mfe_multibranch_fast_init(unsigned int length)
{
  return get_aux_arrays(length);
}


PUBLIC void
vrna_mfe_multibranch_fast_rotate(struct vrna_mx_mfe_aux_ml_s *aux)
{
  if (aux)
    rotate_aux_arrays(aux);
}


PUBLIC void
vrna_mfe_multibranch_fast_free(struct vrna_mx_mfe_aux_ml_s *aux)
{
  if (aux)
    free_aux_arrays(aux);
}


PUBLIC int
vrna_mfe_multibranch_loop_fast(vrna_fold_compound_t         *fc,
                               unsigned int                 i,
                               unsigned int                 j,
                               struct vrna_mx_mfe_aux_ml_s  *helpers)
{
  if ((fc) &&
      (helpers))
    return E_mb_loop_fast(fc,
                          i,
                          j,
                          helpers);

  return INF;
}


PUBLIC int
vrna_mfe_multibranch_stems_fast(vrna_fold_compound_t        *fc,
                                unsigned int                i,
                                unsigned int                j,
                                struct vrna_mx_mfe_aux_ml_s *helpers)
{
  if ((fc) &&
      (helpers))
    return E_ml_stems_fast(fc,
                           i,
                           j,
                           helpers);

  return INF;
}


PUBLIC int
vrna_mfe_multibranch_m2_fast(vrna_fold_compound_t         *fc,
                             unsigned int                 i,
                             unsigned int                 j,
                             struct vrna_mx_mfe_aux_ml_s  *helpers)
{
  if ((fc) &&
      (helpers))
    return mfe_multibranch_m2_fast(fc, i, j, helpers);

  return INF;
}


PUBLIC int
vrna_mfe_multibranch_loop_stack(vrna_fold_compound_t  *fc,
                                unsigned int          i,
                                unsigned int          j)
{
  if (fc)
    return E_mb_loop_stack(fc, i, j);

  return INF;
}


PUBLIC int
vrna_mfe_multibranch_m1(vrna_fold_compound_t  *fc,
                        unsigned int          i,
                        unsigned int          j)
{
  int e;

  e = INF;

  if ((fc) &&
      (fc->matrices) &&
      (fc->matrices->fM1)) {
    struct hc_mb_def_dat  hc_dat_local;
    struct sc_mb_dat      sc_wrapper;
    vrna_hc_eval_f        evaluate;

    evaluate = prepare_hc_mb_def(fc, &hc_dat_local);
    init_sc_mb(fc, &sc_wrapper);

    e = extend_fm_3p(i, j, fc->matrices->fM1, fc, evaluate, &hc_dat_local, &sc_wrapper);

    if (fc->aux_grammar) {
      for (size_t c = 0; c < vrna_array_size(fc->aux_grammar->m1); c++) {
        if (fc->aux_grammar->m1[c].cb) {
          int ee = fc->aux_grammar->m1[c].cb(fc, i, j, fc->aux_grammar->m1[c].data);
          e = MIN2(e, ee);
        }
      }
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
PRIVATE struct vrna_mx_mfe_aux_ml_s *
get_aux_arrays(unsigned int length)
{
  unsigned int                j;
  struct vrna_mx_mfe_aux_ml_s *aux;

  aux         = (struct vrna_mx_mfe_aux_ml_s *)vrna_alloc(sizeof(struct vrna_mx_mfe_aux_ml_s));
  aux->n      = length;
  aux->Fmi    = (int *)vrna_alloc(sizeof(int) * (length + 1));  /* holds row i of fML (avoids jumps in memory)  */
  aux->DMLi   = (int *)vrna_alloc(sizeof(int) * (length + 1));  /* DMLi[j] holds  MIN(fML[i,k]+fML[k+1,j])      */
  aux->DMLi1  = (int *)vrna_alloc(sizeof(int) * (length + 1));  /*                MIN(fML[i+1,k]+fML[k+1,j])    */
  aux->DMLi2  = (int *)vrna_alloc(sizeof(int) * (length + 1));  /*                MIN(fML[i+2,k]+fML[k+1,j])    */

  /* prefill helper arrays */
  for (j = 0; j <= length; j++)
    aux->Fmi[j] = aux->DMLi[j] = aux->DMLi1[j] = aux->DMLi2[j] = INF;

  return aux;
}


PRIVATE void
rotate_aux_arrays(struct vrna_mx_mfe_aux_ml_s *aux)
{
  unsigned int  j;
  int           *FF;

  FF          = aux->DMLi2;
  aux->DMLi2  = aux->DMLi1;
  aux->DMLi1  = aux->DMLi;
  aux->DMLi   = FF;

  for (j = 1; j <= aux->n; j++)
    aux->Fmi[j] = aux->DMLi[j] = INF;
}


PRIVATE void
free_aux_arrays(struct vrna_mx_mfe_aux_ml_s *aux)
{
  free(aux->Fmi);
  free(aux->DMLi);
  free(aux->DMLi1);
  free(aux->DMLi2);
  free(aux);
}


PRIVATE INLINE int
ml_pair_d0(vrna_fold_compound_t *fc,
           unsigned int         i,
           unsigned int         j,
           int                  *dmli1,
           vrna_hc_eval_f       evaluate,
           struct hc_mb_def_dat *hc_wrapper,
           struct sc_mb_dat     *sc_wrapper)
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

          e += vrna_E_multibranch_stem(tt, -1, -1, P) +
               P->MLclosing;
          break;

        case VRNA_FC_TYPE_COMPARATIVE:
          n_seq = fc->n_seq;
          SS    = fc->S;
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


PRIVATE INLINE int
ml_pair_d1(vrna_fold_compound_t *fc,
           unsigned int         i,
           unsigned int         j,
           int                  *dmli1,
           int                  *dmli2,
           vrna_hc_eval_f       evaluate,
           struct hc_mb_def_dat *hc_wrapper,
           struct sc_mb_dat     *sc_wrapper)
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
ml_pair_d2(vrna_fold_compound_t *fc,
           unsigned int         i,
           unsigned int         j,
           int                  *dmli1,
           vrna_hc_eval_f       evaluate,
           struct hc_mb_def_dat *hc_wrapper,
           struct sc_mb_dat     *sc_wrapper)
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
ml_pair5(vrna_fold_compound_t *fc,
         unsigned int         i,
         unsigned int         j,
         int                  *dmli2,
         vrna_hc_eval_f       evaluate,
         struct hc_mb_def_dat *hc_wrapper,
         struct sc_mb_dat     *sc_wrapper)
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

          si1 = ((strands == 1) || (sn[i] == sn[i + 2])) ? S[i + 1] : -1;

          e += vrna_E_multibranch_stem(tt, -1, si1, P) +
               P->MLclosing +
               P->MLbase;
          break;

        case VRNA_FC_TYPE_COMPARATIVE:
          n_seq = fc->n_seq;
          SS    = fc->S;
          S3    = fc->S3;

          for (s = 0; s < n_seq; s++) {
            tt  = vrna_get_ptype_md(SS[s][j], SS[s][i], md);
            e   += vrna_E_multibranch_stem(tt, -1, S3[s][i], P);
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
ml_pair3(vrna_fold_compound_t *fc,
         unsigned int         i,
         unsigned int         j,
         int                  *dmli1,
         vrna_hc_eval_f       evaluate,
         struct hc_mb_def_dat *hc_wrapper,
         struct sc_mb_dat     *sc_wrapper)
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

          sj1 = ((strands == 1) || (sn[j - 2] == sn[j])) ? S[j - 1] : -1;

          e += vrna_E_multibranch_stem(tt, sj1, -1, P) +
               P->MLclosing +
               P->MLbase;
          break;

        case VRNA_FC_TYPE_COMPARATIVE:
          n_seq = fc->n_seq;
          SS    = fc->S;
          S5    = fc->S5;

          for (s = 0; s < n_seq; s++) {
            tt  = vrna_get_ptype_md(SS[s][j], SS[s][i], md);
            e   += vrna_E_multibranch_stem(tt, S5[s][j], -1, P);
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
ml_pair53(vrna_fold_compound_t  *fc,
          unsigned int          i,
          unsigned int          j,
          int                   *dmli1 VRNA_UNUSED,
          int                   *dmli2,
          vrna_hc_eval_f        evaluate,
          struct hc_mb_def_dat  *hc_wrapper,
          struct sc_mb_dat      *sc_wrapper)
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

          si1 = ((strands == 1) || (sn[i] == sn[i + 2])) ? S[i + 1] : -1;
          sj1 = ((strands == 1) || (sn[j - 2] == sn[j])) ? S[j - 1] : -1;

          e += vrna_E_multibranch_stem(tt, sj1, si1, P) +
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
            e   += vrna_E_multibranch_stem(tt, S5[s][j], S3[s][i], P);
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
E_mb_loop_fast(vrna_fold_compound_t         *fc,
               unsigned int                 i,
               unsigned int                 j,
               struct vrna_mx_mfe_aux_ml_s  *helpers)
{
  unsigned int          dangle_model;
  int                   decomp, *dmli1, *dmli2;
  vrna_param_t          *P;
  vrna_md_t             *md;
  vrna_hc_eval_f        evaluate;
  struct hc_mb_def_dat  hc_dat_local;
  struct sc_mb_dat      sc_wrapper;

  P             = fc->params;
  md            = &(P->model_details);
  dangle_model  = md->dangles;
  dmli1         = helpers->DMLi1;
  dmli2         = helpers->DMLi2;

  /* init values */
  decomp = INF;

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

  return decomp;
}


PRIVATE int
E_mb_loop_stack(vrna_fold_compound_t  *fc,
                unsigned int          i,
                unsigned int          j)
{
  char                  *ptype, **ptype_local;
  short                 **SS;
  unsigned int          n_seq, s, type, type_2, *tt, sliding_window, k;
  int                   *c, *fML, e, decomp, en, i1k, k1j1, ij, *indx,
                        *rtype, **c_local, **fML_local;
  vrna_param_t          *P;
  vrna_md_t             *md;
  vrna_hc_eval_f        evaluate;
  struct hc_mb_def_dat  hc_dat_local;
  struct sc_mb_dat      sc_wrapper;

  sliding_window = (fc->hc->type == VRNA_HC_WINDOW) ? 1 : 0;

  n_seq       = (fc->type == VRNA_FC_TYPE_SINGLE) ? 1 : fc->n_seq;
  SS          = (fc->type == VRNA_FC_TYPE_SINGLE) ? NULL : fc->S;
  indx        = fc->jindx;
  P           = fc->params;
  md          = &(P->model_details);
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
      for (k = i + 2; k < j - 2; k++) {
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
      k1j1 = indx[j - 1] + i + 2 + 1;
      for (k = i + 2; k < j - 2; k++, k1j1++) {
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
extend_fm_3p(unsigned int         i,
             unsigned int         j,
             int                  *fm,
             vrna_fold_compound_t *fc,
             vrna_hc_eval_f       evaluate,
             struct hc_mb_def_dat *hc_dat_local,
             struct sc_mb_dat     *sc_wrapper)
{
  short         *S, **SS, **S5, **S3;
  unsigned int  length, *sn, n_seq, s, sliding_window, type,
                dangle_model, with_gquad, u, k, cnt, with_ud;
  int           en, en2, *indx, *c, **c_local, **fm_local, **ggg_local, ij, e;
  vrna_param_t  *P;
  vrna_md_t     *md;
  vrna_ud_t     *domains_up;

  vrna_smx_csr(int) * c_gq;

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
  c_gq            = (sliding_window) ? NULL : fc->matrices->c_gq;
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
            en += vrna_E_multibranch_stem(type, (i == 1) ? S[length] : S[i - 1], S[j + 1], P);
          else
            en += vrna_E_multibranch_stem(type, -1, -1, P);

          break;

        case VRNA_FC_TYPE_COMPARATIVE:
          if (dangle_model == 2) {
            for (s = 0; s < n_seq; s++) {
              type  = vrna_get_ptype_md(SS[s][i], SS[s][j], md);
              en    += vrna_E_multibranch_stem(type, S5[s][i], S3[s][j], P);
            }
          } else {
            for (s = 0; s < n_seq; s++) {
              type  = vrna_get_ptype_md(SS[s][i], SS[s][j], md);
              en    += vrna_E_multibranch_stem(type, -1, -1, P);
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
#ifndef VRNA_DISABLE_C11_FEATURES
      en = (sliding_window) ? ggg_local[i][j - i] : vrna_smx_csr_get(c_gq, i, j, INF);
#else
      en = (sliding_window) ? ggg_local[i][j - i] : vrna_smx_csr_int_get(c_gq, i, j, INF);
#endif
      if (en != INF) {
        en += vrna_E_multibranch_stem(0, -1, -1, P) *
              n_seq;

        e = MIN2(e, en);
      }
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

  if (dangle_model % 2) {
    if ((i + 1 < j) &&
        (evaluate(i, j, i + 1, j, VRNA_DECOMP_ML_STEM, hc_dat_local))) {
      en = (sliding_window) ? c_local[i + 1][j - i - 1] : c[indx[j] + i + 1];
      if (en != INF) {
        switch (fc->type) {
          case VRNA_FC_TYPE_SINGLE:
            type = (sliding_window) ?
                   vrna_get_ptype_window(i + 1, j, fc->ptype_local) :
                   vrna_get_ptype(indx[j] + i + 1, fc->ptype);

            en += vrna_E_multibranch_stem(type, S[i], -1, P);
            break;

          case VRNA_FC_TYPE_COMPARATIVE:
            for (s = 0; s < n_seq; s++) {
              type  = vrna_get_ptype_md(SS[s][i + 1], SS[s][j], md);
              en    += vrna_E_multibranch_stem(type, S5[s][i + 1], -1, P);
            }
            break;
        }

        if (sc_wrapper->red_stem)
          en += sc_wrapper->red_stem(i, j, i + 1, j, sc_wrapper);

        en += P->MLbase *
              n_seq;

        e = MIN2(e, en);
      }
    }

    if ((i + 1 < j) &&
        (evaluate(i, j, i, j - 1, VRNA_DECOMP_ML_STEM, hc_dat_local))) {
      en = (sliding_window) ? c_local[i][j - 1 - i] : c[indx[j - 1] + i];
      if (en != INF) {
        switch (fc->type) {
          case VRNA_FC_TYPE_SINGLE:
            type = (sliding_window) ?
                   vrna_get_ptype_window(i, j - 1, fc->ptype_local) :
                   vrna_get_ptype(indx[j - 1] + i, fc->ptype);

            en += vrna_E_multibranch_stem(type, -1, S[j], P);
            break;

          case VRNA_FC_TYPE_COMPARATIVE:
            for (s = 0; s < n_seq; s++) {
              type  = vrna_get_ptype_md(SS[s][i], SS[s][j - 1], md);
              en    += vrna_E_multibranch_stem(type, -1, S3[s][j], P);
            }
            break;
        }

        if (sc_wrapper->red_stem)
          en += sc_wrapper->red_stem(i, j, i, j - 1, sc_wrapper);

        en += P->MLbase *
              n_seq;

        e = MIN2(e, en);
      }
    }

    if ((i + 2 < j) &&
        (evaluate(i, j, i + 1, j - 1, VRNA_DECOMP_ML_STEM, hc_dat_local))) {
      en = (sliding_window) ? c_local[i + 1][j - 1 - i - 1] : c[indx[j - 1] + i + 1];
      if (en != INF) {
        switch (fc->type) {
          case VRNA_FC_TYPE_SINGLE:
            type = (sliding_window) ?
                   vrna_get_ptype_window(i + 1, j - 1, fc->ptype_local) :
                   vrna_get_ptype(indx[j - 1] + i + 1, fc->ptype);

            en += vrna_E_multibranch_stem(type, S[i], S[j], P);
            break;

          case VRNA_FC_TYPE_COMPARATIVE:
            for (s = 0; s < n_seq; s++) {
              type  = vrna_get_ptype_md(SS[s][i + 1], SS[s][j - 1], md);
              en    += vrna_E_multibranch_stem(type, S5[s][i], S3[s][j], P);
            }
            break;
        }

        if (sc_wrapper->red_stem)
          en += sc_wrapper->red_stem(i, j, i + 1, j - 1, sc_wrapper);

        en += 2 * P->MLbase *
              n_seq;

        e = MIN2(e, en);
      }
    }
  }

  if (with_ud) {
    for (cnt = 0; cnt < (unsigned int)domains_up->uniq_motif_count; cnt++) {
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
mfe_multibranch_m2_fast(vrna_fold_compound_t        *fc,
                        unsigned int                i,
                        unsigned int                j,
                        struct vrna_mx_mfe_aux_ml_s *helpers)
{
  unsigned int      k, *sn, *se, sliding_window;
  int               en, decomp, k1j, *indx, *fm, **fm_local, *fmi, *dmli;
  vrna_hc_t         *hc;
  struct sc_mb_dat  sc_wrapper;

  sliding_window = (fc->hc->type == VRNA_HC_WINDOW) ? 1 : 0;

  indx      = (sliding_window) ? NULL : fc->jindx;
  sn        = fc->strand_number;
  se        = fc->strand_end;
  hc        = fc->hc;
  fm        = (sliding_window) ? NULL : fc->matrices->fML;
  fm_local  = (sliding_window) ? fc->matrices->fML_local : NULL;
  fmi       = helpers->Fmi;
  dmli      = helpers->DMLi;
  decomp    = INF;

  init_sc_mb(fc, &sc_wrapper);

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
    for (k = i + 1; k <= j - 2; k++)
      fmi_tmp[k] = fmi[k];

    /* mask unavailable decompositions */
    for (k = i + 1; k <= j - 2; k++)
      if (!hc->f(i, j, k, k + 1, VRNA_DECOMP_ML_ML_ML, hc->data))
        fmi_tmp[k] = INF;
  }

  if (sc_wrapper.decomp_ml) {
    if (fmi_tmp == fmi) {
      fmi_tmp = (int *)vrna_alloc(sizeof(int) * (j - i + 2));
      fmi_tmp -= i;

      /* copy data */
      for (k = i + 1; k <= j - 2; k++)
        fmi_tmp[k] = fmi[k];
    }

    for (k = i + 1; k <= j - 2; k++)
      if (fmi_tmp[k] != INF)
        fmi_tmp[k] += sc_wrapper.decomp_ml(i, j, k, k + 1, &sc_wrapper);
  }

  /* modular decomposition -------------------------------*/
  if (sliding_window) {
    for (decomp = INF, k = i + 1; k <= j - 2; k++) {
      if ((fmi_tmp[k] != INF) && (fm_local[k + 1][j - (k + 1)] != INF)) {
        en      = fmi_tmp[k] + fm_local[k + 1][j - (k + 1)];
        decomp  = MIN2(decomp, en);
      }
    }
  } else {
    decomp  = INF;
    k       = i + 1;
    if (k >= j)
      k = j - 1;

    k1j = indx[j] + k + 1;

    /*
     *  loop over entire range but skip decompositions with in-between strand nick,
     *  this should be faster than evaluating hard constraints callback for each
     *  decomposition
     */
    while (1) {
      unsigned int last_nt = se[sn[k]]; /* go to last nucleotide of current strand */
      if (last_nt > j - 2)
        last_nt = j - 2;                /* at most go to last possible decomposition split before reaching j */

      if (last_nt < i)
        last_nt = i; /* do not start before i */

      const unsigned int count = last_nt - k;

      en      = vrna_fun_zip_add_min(fmi_tmp + k, fm + k1j, count);
      decomp  = MIN2(decomp, en);

      /* advance counters by processed subsegment and add 1 for the split point between strands */
      k   += count + 1;
      k1j += count + 1;

      if (k > j - 2)
        break;
    }
  }

  /* end modular decomposition -------------------------------*/

  if (fmi_tmp != fmi) {
    fmi_tmp += i;
    free(fmi_tmp);
  }

  /* apply auxiliary grammar rule for multibranch loop case */
  if (fc->aux_grammar) {
    for (size_t c = 0; c < vrna_array_size(fc->aux_grammar->m2); c++)
      if (fc->aux_grammar->m2[c].cb)
        decomp = MIN2(decomp, fc->aux_grammar->m2[c].cb(fc, i, j, fc->aux_grammar->m2[c].data));
  }

  free_sc_mb(&sc_wrapper);

  return decomp;               /* store for use in fast ML decompositon */
}


PRIVATE int
E_ml_stems_fast(vrna_fold_compound_t        *fc,
                unsigned int                i,
                unsigned int                j,
                struct vrna_mx_mfe_aux_ml_s *helpers)
{
  char                  *ptype, **ptype_local;
  short                 *S, **SS, **S5, **S3;
  unsigned int          length, k, u, *sn, *se, n_seq, s, type_2, dangle_model, type,
                        cnt, with_ud, sliding_window, circular;
  int                   en, decomp, mm5, mm3, k1j, *indx, *c, *fm, *fm2,
                        ij, *rtype, e, **c_local, **fm_local, *fmi, *dmli;
  vrna_param_t          *P;
  vrna_md_t             *md;
  vrna_ud_t             *domains_up;
  vrna_hc_eval_f        evaluate;
  struct hc_mb_def_dat  hc_dat_local;
  struct sc_mb_dat      sc_wrapper;

  sliding_window = (fc->hc->type == VRNA_HC_WINDOW) ? 1 : 0;

  length      = fc->length;
  ptype       = (fc->type == VRNA_FC_TYPE_SINGLE) ? ((sliding_window) ? NULL : fc->ptype) : NULL;
  ptype_local = (fc->type == VRNA_FC_TYPE_SINGLE) ? ((sliding_window) ? fc->ptype_local : NULL) : NULL;
  S             = (fc->type == VRNA_FC_TYPE_SINGLE) ? fc->sequence_encoding : NULL;
  SS            = (fc->type == VRNA_FC_TYPE_SINGLE) ? NULL : fc->S;
  S5            = (fc->type == VRNA_FC_TYPE_SINGLE) ? NULL : fc->S5;
  S3            = (fc->type == VRNA_FC_TYPE_SINGLE) ? NULL : fc->S3;
  indx          = (sliding_window) ? NULL : fc->jindx;
  n_seq         = (fc->type == VRNA_FC_TYPE_SINGLE) ? 1 : fc->n_seq;
  sn            = fc->strand_number;
  se            = fc->strand_end;
  c             = (sliding_window) ? NULL : fc->matrices->c;
  c_local       = (sliding_window) ? fc->matrices->c_local : NULL;
  fm            = (sliding_window) ? NULL : fc->matrices->fML;
  fm_local      = (sliding_window) ? fc->matrices->fML_local : NULL;
  fm2           = (sliding_window) ? NULL : fc->matrices->fM2_real;
  P             = fc->params;
  md            = &(P->model_details);
  ij            = (sliding_window) ? 0 : indx[j] + i;
  dangle_model  = md->dangles;
  rtype         = &(md->rtype[0]);
  circular      = md->circ;
  domains_up    = fc->domains_up;
  with_ud       = (domains_up && domains_up->energy_cb) ? 1 : 0;
  e             = INF;
  fmi           = helpers->Fmi;
  dmli          = helpers->DMLi;

  evaluate = prepare_hc_mb_def(fc, &hc_dat_local);

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
  if ((i + 1 < j) &&
      (evaluate(i, j, i + 1, j, VRNA_DECOMP_ML_ML, &hc_dat_local))) {
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
    for (cnt = 0; cnt < (unsigned int)domains_up->uniq_motif_count; cnt++) {
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

    if ((i + 1 < j) &&
        (evaluate(i, j, i + 1, j, VRNA_DECOMP_ML_STEM, &hc_dat_local))) {
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
            en += vrna_E_multibranch_stem(type, mm5, -1, P);
            break;

          case VRNA_FC_TYPE_COMPARATIVE:
            for (s = 0; s < n_seq; s++) {
              type  = vrna_get_ptype_md(SS[s][i + 1], SS[s][j], md);
              en    += vrna_E_multibranch_stem(type, S5[s][i + 1], -1, P);
            }
            break;
        }

        if (sc_wrapper.red_ml)
          en += sc_wrapper.red_ml(i, j, i + 1, j, &sc_wrapper);

        e = MIN2(e, en);
      }
    }

    if ((i + 1 < j) &&
        (evaluate(i, j, i, j - 1, VRNA_DECOMP_ML_STEM, &hc_dat_local))) {
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
            en += vrna_E_multibranch_stem(type, -1, mm3, P);
            break;

          case VRNA_FC_TYPE_COMPARATIVE:
            for (s = 0; s < n_seq; s++) {
              type  = vrna_get_ptype_md(SS[s][i], SS[s][j - 1], md);
              en    += vrna_E_multibranch_stem(type, -1, S3[s][j - 1], P);
            }
            break;
        }

        if (sc_wrapper.red_ml)
          en += sc_wrapper.red_ml(i, j, i, j - 1, &sc_wrapper);

        e = MIN2(e, en);
      }
    }

    if ((i + 2 < j) &&
        (evaluate(i, j, i + 1, j - 1, VRNA_DECOMP_ML_STEM, &hc_dat_local))) {
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
            en += vrna_E_multibranch_stem(type, mm5, mm3, P);
            break;

          case VRNA_FC_TYPE_COMPARATIVE:
            for (s = 0; s < n_seq; s++) {
              type  = vrna_get_ptype_md(SS[s][i + 1], SS[s][j - 1], md);
              en    += vrna_E_multibranch_stem(type, S5[s][i + 1], S3[s][j - 1], P);
            }
            break;
        }

        if (sc_wrapper.red_ml)
          en += sc_wrapper.red_ml(i, j, i + 1, j - 1, &sc_wrapper);

        e = MIN2(e, en);
      }
    } /* end special cases for dangles == 1 || dangles == 3 */
  }

  /*
   * modular decomposition -------------------------------
   * get contribution with at least 2 stems
   * here, we assume that qm2[i,j] has been filled already
   */
  decomp = (fm2) ? fm2[ij] : mfe_multibranch_m2_fast(fc, i, j, helpers);

  if (sliding_window) {
    fmi   -= i;
    dmli  -= i;
  }

  dmli[j] = decomp;               /* store for use in fast ML decompositon */

  e = MIN2(e, decomp);

  /* coaxial stacking */
  if (dangle_model == 3) {
    /* additional ML decomposition as two coaxially stacked helices */
    if (sliding_window) {
      /* additional ML decomposition as two coaxially stacked helices */
      for (decomp = INF, k = i + 1; k <= j - 2; k++) {
        if (evaluate(i, k, k + 1, j, VRNA_DECOMP_ML_COAXIAL_ENC, &hc_dat_local)) {
          type    = rtype[vrna_get_ptype_window(i, k, ptype_local)];
          type_2  = rtype[vrna_get_ptype_window(k + 1, j, ptype_local)];

          en = c_local[i][k - i] +
               c_local[k + 1][j - k - 1] +
               P->stack[type][type_2];

          if (sc_wrapper.coaxial_enc)
            en += sc_wrapper.coaxial_enc(i, k, k + 1, j, &sc_wrapper);

          decomp = MIN2(decomp, en);
        }
      }
    } else {
      int ik;
      k1j     = indx[j] + i + 2;
      decomp  = INF;
      k       = i + 1;
      if (k >= j)
        k = j - 1;

      /*
       *  loop over entire range but skip decompositions with in-between strand nick,
       *  this should be faster than evaluating hard constraints callback for each
       *  decomposition
       */
      while (1) {
        unsigned int last_nt = se[sn[k - 1]]; /* go to last nucleotide of current strand */
        if (last_nt > j - 2)
          last_nt = j - 2;                    /* at most go to last possible decomposition split before reaching j */

        if (last_nt < i)
          last_nt = i; /* do not start before i */

        const unsigned int stop = last_nt;
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

        if (k > j - 2)
          break;
      }
    }

    /* no TermAU penalty if coax stack */
    decomp += 2 * P->MLintern[1] *
              n_seq;
#if 0
    /*
     * This is needed for Y shaped ML loops with coax stacking of
     * enclosed parts, but backtracking will fail if activated
     */
    DMLi[j] = MIN2(DMLi[j], decomp);
    DMLi[j] = MIN2(DMLi[j], DMLi[j - 1] + P->MLbase);
    DMLi[j] = MIN2(DMLi[j], DMLi1[j] + P->MLbase);
    new_fML = MIN2(new_fML, DMLi[j]);
#endif
    e = MIN2(e, decomp);
  }

  if (fc->aux_grammar) {
    for (size_t cnt = 0; cnt < vrna_array_size(fc->aux_grammar->m); cnt++) {
      if (fc->aux_grammar->m[cnt].cb) {
        en  = fc->aux_grammar->m[cnt].cb(fc, i, j, fc->aux_grammar->m[cnt].data);
        e   = MIN2(e, en);
      }
    }
  }

  fmi[j] = e;

  free_sc_mb(&sc_wrapper);

  return e;
}


/*
 * ###########################################
 * # deprecated functions below              #
 *###########################################
 */

#ifndef VRNA_DISABLE_BACKWARD_COMPATIBILITY

PUBLIC int
vrna_E_mb_loop_fast(vrna_fold_compound_t  *fc,
                    int                   i,
                    int                   j,
                    int                   *dmli1,
                    int                   *dmli2)
{
  if ((fc) &&
      (dmli1) &&
      (dmli2)) {
    struct vrna_mx_mfe_aux_ml_s tmp;

    tmp.n     = fc->length;
    tmp.DMLi1 = dmli1;
    tmp.DMLi2 = dmli2;

    return E_mb_loop_fast(fc, i, j, &tmp);
  }

  return INF;
}


PUBLIC int
vrna_E_ml_stems_fast(vrna_fold_compound_t *fc,
                     int                  i,
                     int                  j,
                     int                  *fmi,
                     int                  *dmli)
{
  if (fc) {
    struct vrna_mx_mfe_aux_ml_s tmp;

    tmp.n     = fc->length;
    tmp.Fmi   = fmi;
    tmp.DMLi  = dmli;

    return E_ml_stems_fast(fc, i, j, &tmp);
  }

  return INF;
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
  return vrna_mfe_multibranch_m1(fc,
                                 (unsigned int)i,
                                 (unsigned int)j);
}


#endif
