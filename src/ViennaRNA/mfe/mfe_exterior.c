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
#include "ViennaRNA/sequences/alphabet.h"
#include "ViennaRNA/constraints/hard.h"
#include "ViennaRNA/constraints/soft.h"
#include "ViennaRNA/eval/gquad.h"
#include "ViennaRNA/structured_domains.h"
#include "ViennaRNA/unstructured_domains.h"
#include "ViennaRNA/eval/structures.h"
#include "ViennaRNA/eval/exterior.h"
#include "ViennaRNA/utils/higher_order_functions.h"

#include "ViennaRNA/intern/grammar_dat.h"

#ifdef __GNUC__
# define INLINE inline
#else
# define INLINE
#endif

#include "ViennaRNA/constraints/exterior_hc.inc"
#include "ViennaRNA/constraints/exterior_sc.inc"

#include "ViennaRNA/mfe/exterior.h"

/*
 #################################
 # PRIVATE FUNCTION DECLARATIONS #
 #################################
 */
PRIVATE INLINE int
reduce_f5_up(vrna_fold_compound_t   *fc,
             unsigned int           j,
             vrna_hc_eval_f         evaluate,
             struct hc_ext_def_dat  *hc_dat_local,
             struct sc_f5_dat       *sc_wrapper);


PRIVATE INLINE int *
get_stem_contributions_d0(vrna_fold_compound_t  *fc,
                          unsigned int          j,
                          vrna_hc_eval_f        evaluate,
                          struct hc_ext_def_dat *hc_dat_local,
                          struct sc_f5_dat      *sc_wrapper);


PRIVATE INLINE int *
get_stem_contributions_d2(vrna_fold_compound_t  *fc,
                          unsigned int          j,
                          vrna_hc_eval_f        evaluate,
                          struct hc_ext_def_dat *hc_dat_local,
                          struct sc_f5_dat      *sc_wrapper);


PRIVATE INLINE int *
f5_get_stem_contributions_d5(vrna_fold_compound_t   *fc,
                             unsigned int           j,
                             vrna_hc_eval_f         evaluate,
                             struct hc_ext_def_dat  *hc_dat_local,
                             struct sc_f5_dat       *sc_wrapper);


PRIVATE INLINE int *
f5_get_stem_contributions_d3(vrna_fold_compound_t   *fc,
                             unsigned int           j,
                             vrna_hc_eval_f         evaluate,
                             struct hc_ext_def_dat  *hc_dat_local,
                             struct sc_f5_dat       *sc_wrapper);


PRIVATE INLINE int *
f5_get_stem_contributions_d53(vrna_fold_compound_t  *fc,
                              unsigned int          j,
                              vrna_hc_eval_f        evaluate,
                              struct hc_ext_def_dat *hc_dat_local,
                              struct sc_f5_dat      *sc_wrapper);


PRIVATE INLINE int
decompose_f5_ext_stem(vrna_fold_compound_t  *fc,
                      unsigned int          j,
                      int                   *stems);


PRIVATE INLINE int
decompose_f5_ext_stem_d0(vrna_fold_compound_t   *fc,
                         unsigned int           j,
                         vrna_hc_eval_f         evaluate,
                         struct hc_ext_def_dat  *hc_dat_local,
                         struct sc_f5_dat       *sc_wrapper);


PRIVATE INLINE int
decompose_f5_ext_stem_d2(vrna_fold_compound_t   *fc,
                         unsigned int           j,
                         vrna_hc_eval_f         evaluate,
                         struct hc_ext_def_dat  *hc_dat_local,
                         struct sc_f5_dat       *sc_wrapper);


PRIVATE INLINE int
decompose_f5_ext_stem_d1(vrna_fold_compound_t   *fc,
                         unsigned int           j,
                         vrna_hc_eval_f         evaluate,
                         struct hc_ext_def_dat  *hc_dat_local,
                         struct sc_f5_dat       *sc_wrapper);


PRIVATE INLINE int
add_f5_gquad(vrna_fold_compound_t   *fc,
             unsigned int           j,
             vrna_hc_eval_f         evaluate,
             struct hc_ext_def_dat  *hc_dat_local,
             struct sc_f5_dat       *sc_wrapper);


/*
 #################################
 # BEGIN OF FUNCTION DEFINITIONS #
 #################################
 */
PUBLIC int
vrna_mfe_exterior_f5(vrna_fold_compound_t *fc)
{
  if (fc) {
    unsigned int          j, length, dangle_model, with_gquad;
    int                   en, *f5;
    vrna_param_t          *P;
    vrna_hc_eval_f        evaluate;
    struct hc_ext_def_dat hc_dat_local;
    struct sc_f5_dat      sc_wrapper;
    vrna_gr_aux_t         grammar;

    length        = fc->length;
    f5            = fc->matrices->f5;
    P             = fc->params;
    dangle_model  = P->model_details.dangles;
    with_gquad    = P->model_details.gquad;
    grammar       = fc->aux_grammar;
    evaluate      = prepare_hc_ext_def(fc, &hc_dat_local);

    init_sc_f5(fc, &sc_wrapper);

    f5[0] = 0;
    f5[1] = reduce_f5_up(fc, 1, evaluate, &hc_dat_local, &sc_wrapper);

    if (grammar) {
      for (size_t c = 0; c < vrna_array_size(grammar->f); c++) {
        if (grammar->f[c].cb) {
          en    = grammar->f[c].cb(fc, 1, 1, grammar->f[c].data);
          f5[1] = MIN2(f5[1], en);
        }
      }
    }

    /*
     *  duplicated code may be faster than conditions inside loop or even
     *  using a function pointer ;)
     */
    switch (dangle_model) {
      case 2:
        for (j = 2; j <= length; j++) {
          /* extend previous solution(s) by adding an unpaired region */
          f5[j] = reduce_f5_up(fc, j, evaluate, &hc_dat_local, &sc_wrapper);

          /* decompose into exterior loop part followed by a stem */
          en    = decompose_f5_ext_stem_d2(fc, j, evaluate, &hc_dat_local, &sc_wrapper);
          f5[j] = MIN2(f5[j], en);

          if (with_gquad) {
            en    = add_f5_gquad(fc, j, evaluate, &hc_dat_local, &sc_wrapper);
            f5[j] = MIN2(f5[j], en);
          }

          if (grammar) {
            for (size_t c = 0; c < vrna_array_size(grammar->f); c++) {
              if (grammar->f[c].cb) {
                en    = grammar->f[c].cb(fc, 1, j, grammar->f[c].data);
                f5[j] = MIN2(f5[j], en);
              }
            }
          }
        }
        break;

      case 0:
        for (j = 2; j <= length; j++) {
          /* extend previous solution(s) by adding an unpaired region */
          f5[j] = reduce_f5_up(fc, j, evaluate, &hc_dat_local, &sc_wrapper);

          /* decompose into exterior loop part followed by a stem */
          en    = decompose_f5_ext_stem_d0(fc, j, evaluate, &hc_dat_local, &sc_wrapper);
          f5[j] = MIN2(f5[j], en);

          if (with_gquad) {
            en    = add_f5_gquad(fc, j, evaluate, &hc_dat_local, &sc_wrapper);
            f5[j] = MIN2(f5[j], en);
          }

          if (grammar) {
            for (size_t c = 0; c < vrna_array_size(grammar->f); c++) {
              if (grammar->f[c].cb) {
                en    = grammar->f[c].cb(fc, 1, j, grammar->f[c].data);
                f5[j] = MIN2(f5[j], en);
              }
            }
          }
        }
        break;

      default:
        for (j = 2; j <= length; j++) {
          /* extend previous solution(s) by adding an unpaired region */
          f5[j] = reduce_f5_up(fc, j, evaluate, &hc_dat_local, &sc_wrapper);

          en    = decompose_f5_ext_stem_d1(fc, j, evaluate, &hc_dat_local, &sc_wrapper);
          f5[j] = MIN2(f5[j], en);

          if (with_gquad) {
            en    = add_f5_gquad(fc, j, evaluate, &hc_dat_local, &sc_wrapper);
            f5[j] = MIN2(f5[j], en);
          }

          if (grammar) {
            for (size_t c = 0; c < vrna_array_size(grammar->f); c++) {
              if (grammar->f[c].cb) {
                en    = grammar->f[c].cb(fc, 1, j, grammar->f[c].data);
                f5[j] = MIN2(f5[j], en);
              }
            }
          }
        }
        break;
    }

    free_sc_f5(&sc_wrapper);

    return f5[length];
  }

  return INF;
}


/*
 #####################################
 # BEGIN OF STATIC HELPER FUNCTIONS  #
 #####################################
 */
PRIVATE INLINE int
decompose_f5_ext_stem_d0(vrna_fold_compound_t   *fc,
                         unsigned int           j,
                         vrna_hc_eval_f         evaluate,
                         struct hc_ext_def_dat  *hc_dat_local,
                         struct sc_f5_dat       *sc_wrapper)
{
  int e, *stems;

  stems = get_stem_contributions_d0(fc, j, evaluate, hc_dat_local, sc_wrapper);

  /* 1st case, actual decompostion */
  e = decompose_f5_ext_stem(fc, j, stems);

  /* 2nd case, reduce to single stem */
  e = MIN2(e, stems[1]);

  free(stems);

  return e;
}


PRIVATE INLINE int
decompose_f5_ext_stem_d2(vrna_fold_compound_t   *fc,
                         unsigned int           j,
                         vrna_hc_eval_f         evaluate,
                         struct hc_ext_def_dat  *hc_dat_local,
                         struct sc_f5_dat       *sc_wrapper)
{
  int e, *stems;

  stems = get_stem_contributions_d2(fc, j, evaluate, hc_dat_local, sc_wrapper);

  /* 1st case, actual decompostion */
  e = decompose_f5_ext_stem(fc, j, stems);

  /* 2nd case, reduce to single stem */
  e = MIN2(e, stems[1]);

  free(stems);

  return e;
}


PRIVATE INLINE int
decompose_f5_ext_stem_d1(vrna_fold_compound_t   *fc,
                         unsigned int           j,
                         vrna_hc_eval_f         evaluate,
                         struct hc_ext_def_dat  *hc_dat_local,
                         struct sc_f5_dat       *sc_wrapper)
{
  int e, ee, *stems;

  e = INF;

  /* A) without dangling end contributions */

  /* 1st case, actual decompostion */
  stems = get_stem_contributions_d0(fc, j, evaluate, hc_dat_local, sc_wrapper);

  ee = decompose_f5_ext_stem(fc, j, stems);

  /* 2nd case, reduce to single stem */
  ee = MIN2(ee, stems[1]);

  free(stems);

  e = MIN2(e, ee);

  /* B) with dangling end contribution on 5' side of stem */
  stems = f5_get_stem_contributions_d5(fc, j, evaluate, hc_dat_local, sc_wrapper);

  /* 1st case, actual decompostion */
  ee = decompose_f5_ext_stem(fc, j, stems);

  /* 2nd case, reduce to single stem */
  ee = MIN2(ee, stems[1]);

  free(stems);

  e = MIN2(e, ee);

  /* C) with dangling end contribution on 3' side of stem */
  stems = f5_get_stem_contributions_d3(fc, j, evaluate, hc_dat_local, sc_wrapper);

  /* 1st case, actual decompostion */
  ee = decompose_f5_ext_stem(fc, j, stems);

  /* 2nd case, reduce to single stem */
  ee = MIN2(ee, stems[1]);

  free(stems);

  e = MIN2(e, ee);

  /* D) with dangling end contribution on both sides of stem */
  stems = f5_get_stem_contributions_d53(fc, j, evaluate, hc_dat_local, sc_wrapper);

  /* 1st case, actual decompostion */
  ee = decompose_f5_ext_stem(fc, j, stems);

  /* 2nd case, reduce to single stem */
  ee = MIN2(ee, stems[1]);

  free(stems);

  e = MIN2(e, ee);

  return e;
}


/*
 *  extend f5 by adding an unpaired nucleotide or an unstructured domain
 *  to the 3' end
 */
PRIVATE INLINE int
reduce_f5_up(vrna_fold_compound_t   *fc,
             unsigned int           j,
             vrna_hc_eval_f         evaluate,
             struct hc_ext_def_dat  *hc_dat_local,
             struct sc_f5_dat       *sc_wrapper)
{
  unsigned int  u, k;
  int           e, en, *f5;
  vrna_ud_t     *domains_up;
  sc_f5_cb      sc_red_ext;

  f5          = fc->matrices->f5;
  domains_up  = fc->domains_up;
  e           = INF;

  sc_red_ext = sc_wrapper->red_ext5;

  /* check for 3' extension with one unpaired nucleotide */
  if (f5[j - 1] != INF) {
    if (evaluate(1, j, 1, j - 1, VRNA_DECOMP_EXT_EXT, hc_dat_local)) {
      e = f5[j - 1];

      if (sc_red_ext)
        e += sc_red_ext(j, 1, j - 1, sc_wrapper);
    }
  }

  if ((domains_up) && (domains_up->energy_cb)) {
    for (k = 0; k < (unsigned int)domains_up->uniq_motif_count; k++) {
      u = domains_up->uniq_motif_size[k];
      if ((j >= u) &&
          (f5[j - u] != INF)) {
        if (evaluate(1, j, 1, j - u, VRNA_DECOMP_EXT_EXT, hc_dat_local)) {
          en = f5[j - u] +
               domains_up->energy_cb(fc,
                                     j - u + 1,
                                     j,
                                     VRNA_UNSTRUCTURED_DOMAIN_EXT_LOOP | VRNA_UNSTRUCTURED_DOMAIN_MOTIF,
                                     domains_up->data);

          if (sc_red_ext)
            en += sc_red_ext(j, 1, j - u, sc_wrapper);

          e = MIN2(e, en);
        }
      }
    }
  }

  return e;
}


PRIVATE INLINE int *
get_stem_contributions_d0(vrna_fold_compound_t  *fc,
                          unsigned int          j,
                          vrna_hc_eval_f        evaluate,
                          struct hc_ext_def_dat *hc_dat_local,
                          struct sc_f5_dat      *sc_wrapper)
{
  char          *ptype;
  short         **S;
  unsigned int  i, s, n_seq, type;
  int           ij, *indx, *c, *stems;
  vrna_param_t  *P;
  vrna_md_t     *md;

  sc_f5_cb      sc_spl_stem;
  sc_f5_cb      sc_red_stem;

  stems = (int *)vrna_alloc(sizeof(int) * j);

  P     = fc->params;
  md    = &(P->model_details);
  indx  = fc->jindx;
  c     = fc->matrices->c;
  ij    = indx[j] + j - 1;
  ptype = (fc->type == VRNA_FC_TYPE_SINGLE) ? fc->ptype : NULL;
  n_seq = (fc->type == VRNA_FC_TYPE_SINGLE) ? 1 : fc->n_seq;
  S     = (fc->type == VRNA_FC_TYPE_SINGLE) ? NULL : fc->S;

  sc_spl_stem = sc_wrapper->decomp_stem5;
  sc_red_stem = sc_wrapper->red_stem5;

  switch (fc->type) {
    case VRNA_FC_TYPE_SINGLE:
      for (i = j - 1; i > 1; i--, ij--) {
        stems[i] = INF;

        if ((c[ij] != INF) &&
            (evaluate(1, j, i - 1, i, VRNA_DECOMP_EXT_EXT_STEM, hc_dat_local))) {
          stems[i]  = c[ij];
          type      = vrna_get_ptype(ij, ptype);
          stems[i]  += vrna_E_exterior_stem(type, -1, -1, P);
        }
      }
      break;

    case VRNA_FC_TYPE_COMPARATIVE:
      for (i = j - 1; i > 1; i--, ij--) {
        stems[i] = INF;
        if ((c[ij] != INF) &&
            (evaluate(1, j, i - 1, i, VRNA_DECOMP_EXT_EXT_STEM, hc_dat_local))) {
          stems[i] = c[ij];

          for (s = 0; s < n_seq; s++) {
            type      = vrna_get_ptype_md(S[s][i], S[s][j], md);
            stems[i]  += vrna_E_exterior_stem(type, -1, -1, P);
          }
        }
      }
      break;
  }

  if (sc_spl_stem)
    for (i = j - 1; i > 1; i--)
      if (stems[i] != INF)
        stems[i] += sc_spl_stem(j, i - 1, i, sc_wrapper);

  stems[1]  = INF;
  ij        = indx[j] + 1;

  if ((c[ij] != INF) &&
      (evaluate(1, j, 1, j, VRNA_DECOMP_EXT_STEM, hc_dat_local))) {
    stems[1] = c[ij];

    switch (fc->type) {
      case VRNA_FC_TYPE_SINGLE:
        type      = vrna_get_ptype(ij, ptype);
        stems[1]  += vrna_E_exterior_stem(type, -1, -1, P);
        break;

      case VRNA_FC_TYPE_COMPARATIVE:
        for (s = 0; s < n_seq; s++) {
          type      = vrna_get_ptype_md(S[s][1], S[s][j], md);
          stems[1]  += vrna_E_exterior_stem(type, -1, -1, P);
        }
        break;
    }

    if (sc_red_stem)
      stems[1] += sc_red_stem(j, 1, j, sc_wrapper);
  }

  return stems;
}


PRIVATE INLINE int *
get_stem_contributions_d2(vrna_fold_compound_t  *fc,
                          unsigned int          j,
                          vrna_hc_eval_f        evaluate,
                          struct hc_ext_def_dat *hc_dat_local,
                          struct sc_f5_dat      *sc_wrapper)
{
  char          *ptype;
  short         *S, sj1, *si1, **SS, **S5, **S3, *s3j, *sj;
  unsigned int  n, i, s, n_seq, **a2s, type, *sn;
  int           ij, *indx, *c, *stems, mm5;
  vrna_param_t  *P;
  vrna_md_t     *md;

  sc_f5_cb      sc_spl_stem;
  sc_f5_cb      sc_red_stem;

  stems = (int *)vrna_alloc(sizeof(int) * j);

  n     = fc->length;
  sn    = fc->strand_number;
  P     = fc->params;
  md    = &(P->model_details);
  indx  = fc->jindx;
  c     = fc->matrices->c;
  ij    = indx[j] + j - 1;

  sc_spl_stem = sc_wrapper->decomp_stem5;
  sc_red_stem = sc_wrapper->red_stem5;

  switch (fc->type) {
    case VRNA_FC_TYPE_SINGLE:
      S     = fc->sequence_encoding;
      ptype = fc->ptype;
      si1   = S + j - 2;
      sj1   = ((j < n) && (sn[j] == sn[j + 1])) ? S[j + 1] : -1;

      for (i = j - 1; i > 1; i--, ij--, si1--) {
        stems[i] = INF;
        if ((c[ij] != INF) &&
            (evaluate(1, j, i - 1, i, VRNA_DECOMP_EXT_EXT_STEM, hc_dat_local))) {
          type      = vrna_get_ptype(ij, ptype);
          stems[i]  = c[ij] +
                      vrna_E_exterior_stem(type, *si1, sj1, P);
        }
      }

      if (sc_spl_stem)
        for (i = j - 1; i > 1; i--)
          if (stems[i] != INF)
            stems[i] += sc_spl_stem(j, i - 1, i, sc_wrapper);

      stems[1]  = INF;
      ij        = indx[j] + 1;

      if ((c[ij] != INF) && (evaluate(1, j, 1, j, VRNA_DECOMP_EXT_STEM, hc_dat_local))) {
        type      = vrna_get_ptype(ij, ptype);
        stems[1]  = c[ij] +
                    vrna_E_exterior_stem(type, -1, sj1, P);

        if (sc_red_stem)
          stems[1] += sc_red_stem(j, 1, j, sc_wrapper);
      }

      break;

    case VRNA_FC_TYPE_COMPARATIVE:
      n_seq = fc->n_seq;
      SS    = fc->S;
      S5    = fc->S5;
      S3    = fc->S3;
      a2s   = fc->a2s;

      /* pre-compute S3[s][j - 1] */
      s3j = (short *)vrna_alloc(sizeof(short) * n_seq);
      sj  = (short *)vrna_alloc(sizeof(short) * n_seq);
      for (s = 0; s < n_seq; s++) {
        s3j[s]  = (a2s[s][j] < a2s[s][n]) ? S3[s][j] : -1;
        sj[s]   = SS[s][j];
      }

      for (i = j - 1; i > 1; i--, ij--) {
        stems[i] = INF;
        if ((c[ij] != INF) &&
            (evaluate(1, j, i - 1, i, VRNA_DECOMP_EXT_EXT_STEM, hc_dat_local))) {
          stems[i] = c[ij];
          for (s = 0; s < n_seq; s++) {
            type      = vrna_get_ptype_md(SS[s][i], sj[s], md);
            mm5       = (a2s[s][i] > 1) ? S5[s][i] : -1;
            stems[i]  += vrna_E_exterior_stem(type, mm5, s3j[s], P);
          }
        }
      }

      if (sc_spl_stem)
        for (i = j - 1; i > 1; i--)
          if (stems[i] != INF)
            stems[i] += sc_spl_stem(j, i - 1, i, sc_wrapper);

      stems[1]  = INF;
      ij        = indx[j] + 1;

      if ((c[ij] != INF) && (evaluate(1, j, 1, j, VRNA_DECOMP_EXT_STEM, hc_dat_local))) {
        stems[1] = c[ij];

        for (s = 0; s < n_seq; s++) {
          type      = vrna_get_ptype_md(SS[s][1], sj[s], md);
          stems[1]  += vrna_E_exterior_stem(type, -1, s3j[s], P);
        }

        if (sc_red_stem)
          stems[1] += sc_red_stem(j, 1, j, sc_wrapper);
      }

      free(s3j);
      free(sj);

      break;
  }

  return stems;
}


PRIVATE INLINE int *
f5_get_stem_contributions_d5(vrna_fold_compound_t   *fc,
                             unsigned int           j,
                             vrna_hc_eval_f         evaluate,
                             struct hc_ext_def_dat  *hc_dat_local,
                             struct sc_f5_dat       *sc_wrapper)
{
  char          *ptype;
  short         *S, *si1, **SS, **S5, *sj;
  unsigned int  i, s, n_seq, **a2s, type;
  int           ij, *indx, *c, *stems, mm5;
  vrna_param_t  *P;
  vrna_md_t     *md;

  sc_f5_cb      sc_spl_stem;
  sc_f5_cb      sc_red_stem;

  stems = (int *)vrna_alloc(sizeof(int) * j);

  P     = fc->params;
  md    = &(P->model_details);
  indx  = fc->jindx;
  c     = fc->matrices->c;
  ij    = indx[j] + j;

  sc_spl_stem = sc_wrapper->decomp_stem5;
  sc_red_stem = sc_wrapper->red_stem5;

  switch (fc->type) {
    case VRNA_FC_TYPE_SINGLE:
      S     = fc->sequence_encoding;
      ptype = fc->ptype;
      si1   = S + j - 1;

      for (i = j - 1; i > 1; i--, ij--, si1--) {
        stems[i] = INF;
        if ((c[ij] != INF) &&
            (evaluate(1, j, i - 1, i + 1, VRNA_DECOMP_EXT_EXT_STEM, hc_dat_local))) {
          type      = vrna_get_ptype(ij, ptype);
          stems[i]  = c[ij] +
                      vrna_E_exterior_stem(type, *si1, -1, P);
        }
      }

      if (sc_spl_stem)
        for (i = j - 1; i > 1; i--)
          if (stems[i] != INF)
            stems[i] += sc_spl_stem(j, i - 1, i + 1, sc_wrapper);

      stems[1] = INF;
      if (2 < j) {
        ij = indx[j] + 2;

        if ((c[ij] != INF) && (evaluate(1, j, 2, j, VRNA_DECOMP_EXT_STEM, hc_dat_local))) {
          type      = vrna_get_ptype(ij, ptype);
          stems[1]  = c[ij] +
                      vrna_E_exterior_stem(type, S[1], -1, P);

          if (sc_red_stem)
            stems[1] += sc_red_stem(j, 2, j, sc_wrapper);
        }
      }

      break;

    case VRNA_FC_TYPE_COMPARATIVE:
      n_seq = fc->n_seq;
      SS    = fc->S;
      S5    = fc->S5;
      a2s   = fc->a2s;

      sj = (short *)vrna_alloc(sizeof(short) * n_seq);
      for (s = 0; s < n_seq; s++)
        sj[s] = SS[s][j];

      for (i = j - 1; i > 1; i--, ij--) {
        stems[i] = INF;
        if ((c[ij] != INF) &&
            (evaluate(1, j, i - 1, i + 1, VRNA_DECOMP_EXT_EXT_STEM, hc_dat_local))) {
          stems[i] = c[ij];
          for (s = 0; s < n_seq; s++) {
            type      = vrna_get_ptype_md(SS[s][i + 1], sj[s], md);
            mm5       = (a2s[s][i + 1] > 1) ? S5[s][i + 1] : -1;
            stems[i]  = vrna_E_exterior_stem(type, mm5, -1, P);
          }
        }
      }

      if (sc_spl_stem)
        for (i = j - 1; i > 1; i--)
          if (stems[i] != INF)
            stems[i] += sc_spl_stem(j, i - 1, i + 1, sc_wrapper);

      stems[1] = INF;

      if (2 < j) {
        ij = indx[j] + 2;

        if ((c[ij] != INF) && (evaluate(1, j, 2, j, VRNA_DECOMP_EXT_STEM, hc_dat_local))) {
          stems[1] = c[ij];
          for (s = 0; s < n_seq; s++) {
            type      = vrna_get_ptype_md(SS[s][2], sj[s], md);
            mm5       = (a2s[s][2] > 1) ? S5[s][2] : -1;
            stems[i]  = vrna_E_exterior_stem(type, mm5, -1, P);
          }

          if (sc_red_stem)
            stems[1] += sc_red_stem(j, 2, j, sc_wrapper);
        }

        free(sj);
      }

      break;
  }

  return stems;
}


PRIVATE INLINE int *
f5_get_stem_contributions_d3(vrna_fold_compound_t   *fc,
                             unsigned int           j,
                             vrna_hc_eval_f         evaluate,
                             struct hc_ext_def_dat  *hc_dat_local,
                             struct sc_f5_dat       *sc_wrapper)
{
  char          *ptype;
  short         *S, sj1, **SS, **S3, *s3j1, *ssj1;
  unsigned int  i, n, s, n_seq, **a2s, type;
  int           ij, *indx, *c, *stems;
  vrna_param_t  *P;
  vrna_md_t     *md;

  sc_f5_cb      sc_spl_stem;
  sc_f5_cb      sc_red_stem;

  stems = (int *)vrna_alloc(sizeof(int) * j);

  n     = fc->length;
  P     = fc->params;
  md    = &(P->model_details);
  indx  = fc->jindx;
  c     = fc->matrices->c;
  ij    = indx[j - 1] + j - 1;

  sc_spl_stem = sc_wrapper->decomp_stem51;
  sc_red_stem = sc_wrapper->red_stem5;

  switch (fc->type) {
    case VRNA_FC_TYPE_SINGLE:
      S     = fc->sequence_encoding;
      ptype = fc->ptype;
      sj1   = S[j];

      for (i = j - 1; i > 1; i--, ij--) {
        stems[i] = INF;
        if ((i + 1 < j) &&
            (c[ij] != INF) &&
            (evaluate(1, j, i - 1, i, VRNA_DECOMP_EXT_EXT_STEM1, hc_dat_local))) {
          type      = vrna_get_ptype(ij, ptype);
          stems[i]  = c[ij] +
                      vrna_E_exterior_stem(type, -1, sj1, P);
        }
      }

      if (sc_spl_stem)
        for (i = j - 1; i > 1; i--)
          if (stems[i] != INF)
            stems[i] += sc_spl_stem(j, i - 1, i, sc_wrapper);

      stems[1] = INF;

      if (2 < j) {
        ij = indx[j - 1] + 1;

        if ((c[ij] != INF) && (evaluate(1, j, 1, j - 1, VRNA_DECOMP_EXT_STEM, hc_dat_local))) {
          type      = vrna_get_ptype(ij, ptype);
          stems[1]  = c[ij] +
                      vrna_E_exterior_stem(type, -1, sj1, P);

          if (sc_red_stem)
            stems[1] += sc_red_stem(j, 1, j - 1, sc_wrapper);
        }
      }

      break;

    case VRNA_FC_TYPE_COMPARATIVE:
      n_seq = fc->n_seq;
      SS    = fc->S;
      S3    = fc->S3;
      a2s   = fc->a2s;

      /* pre-compute S3[s][j - 1] */
      s3j1  = (short *)vrna_alloc(sizeof(short) * n_seq);
      ssj1  = (short *)vrna_alloc(sizeof(short) * n_seq);
      for (s = 0; s < n_seq; s++) {
        s3j1[s] = (a2s[s][j - 1] < a2s[s][n]) ? S3[s][j - 1] : -1;
        ssj1[s] = SS[s][j - 1];
      }

      for (i = j - 1; i > 1; i--, ij--) {
        stems[i] = INF;
        if ((i + 1 < j) &&
            (c[ij] != INF) &&
            (evaluate(1, j, i - 1, i, VRNA_DECOMP_EXT_EXT_STEM1, hc_dat_local))) {
          stems[i] = c[ij];
          for (s = 0; s < n_seq; s++) {
            type      = vrna_get_ptype_md(SS[s][i], ssj1[s], md);
            stems[i]  += vrna_E_exterior_stem(type, -1, s3j1[s], P);
          }
        }
      }

      if (sc_spl_stem)
        for (i = j - 1; i > 1; i--)
          if (stems[i] != INF)
            stems[i] += sc_spl_stem(j, i - 1, i, sc_wrapper);

      stems[1] = INF;

      if (2 < j) {
        ij = indx[j - 1] + 1;

        if ((c[ij] != INF) && (evaluate(1, j, 1, j - 1, VRNA_DECOMP_EXT_STEM, hc_dat_local))) {
          stems[1] = c[ij];

          for (s = 0; s < n_seq; s++) {
            type      = vrna_get_ptype_md(SS[s][1], ssj1[s], md);
            stems[1]  += vrna_E_exterior_stem(type, -1, s3j1[s], P);
          }

          if (sc_red_stem)
            stems[1] += sc_red_stem(j, 1, j - 1, sc_wrapper);
        }
      }

      free(s3j1);
      free(ssj1);

      break;
  }

  return stems;
}


PRIVATE INLINE int *
f5_get_stem_contributions_d53(vrna_fold_compound_t  *fc,
                              unsigned int          j,
                              vrna_hc_eval_f        evaluate,
                              struct hc_ext_def_dat *hc_dat_local,
                              struct sc_f5_dat      *sc_wrapper)
{
  char          *ptype;
  short         *S, *si1, sj1, **SS, **S5, **S3, *s3j1, *ssj1;
  unsigned int  i, n, s, n_seq, **a2s, type;
  int           ij, *indx, *c, *stems;
  vrna_param_t  *P;
  vrna_md_t     *md;

  sc_f5_cb      sc_spl_stem;
  sc_f5_cb      sc_red_stem;

  stems = (int *)vrna_alloc(sizeof(int) * j);

  n     = fc->length;
  P     = fc->params;
  md    = &(P->model_details);
  indx  = fc->jindx;
  c     = fc->matrices->c;
  ij    = indx[j - 1] + j;

  sc_spl_stem = sc_wrapper->decomp_stem51;
  sc_red_stem = sc_wrapper->red_stem5;

  switch (fc->type) {
    case VRNA_FC_TYPE_SINGLE:
      S     = fc->sequence_encoding;
      ptype = fc->ptype;
      sj1   = S[j];
      si1   = S + j - 1;

      for (i = j - 1; i > 1; i--, ij--, si1--) {
        stems[i] = INF;
        if ((i + 2 < j) &&
            (c[ij] != INF) &&
            (evaluate(1, j, i - 1, i + 1, VRNA_DECOMP_EXT_EXT_STEM1, hc_dat_local))) {
          type      = vrna_get_ptype(ij, ptype);
          stems[i]  = c[ij] +
                      vrna_E_exterior_stem(type, *si1, sj1, P);
        }
      }

      if (sc_spl_stem)
        for (i = j - 1; i > 1; i--)
          if (stems[i] != INF)
            stems[i] += sc_spl_stem(j, i - 1, i + 1, sc_wrapper);

      stems[1] = INF;

      if (3 < j) {
        ij = indx[j - 1] + 2;

        if ((c[ij] != INF) && (evaluate(1, j, 2, j - 1, VRNA_DECOMP_EXT_STEM, hc_dat_local))) {
          type      = vrna_get_ptype(ij, ptype);
          stems[1]  = c[ij] +
                      vrna_E_exterior_stem(type, S[1], sj1, P);

          if (sc_red_stem)
            stems[1] += sc_red_stem(j, 2, j - 1, sc_wrapper);
        }
      }

      break;

    case VRNA_FC_TYPE_COMPARATIVE:
      n_seq = fc->n_seq;
      SS    = fc->S;
      S5    = fc->S5;
      S3    = fc->S3;
      a2s   = fc->a2s;

      /* pre-compute S3[s][j - 1] */
      s3j1  = (short *)vrna_alloc(sizeof(short) * n_seq);
      ssj1  = (short *)vrna_alloc(sizeof(short) * n_seq);
      for (s = 0; s < n_seq; s++) {
        s3j1[s] = (a2s[s][j - 1] < a2s[s][n]) ? S3[s][j - 1] : -1;
        ssj1[s] = SS[s][j - 1];
      }

      for (i = j - 1; i > 1; i--, ij--) {
        stems[i] = INF;
        if ((i + 1 < j) &&
            (c[ij] != INF) &&
            (evaluate(1, j, i - 1, i + 1, VRNA_DECOMP_EXT_EXT_STEM1, hc_dat_local))) {
          stems[i] = c[ij];
          for (s = 0; s < n_seq; s++) {
            type      = vrna_get_ptype_md(SS[s][i + 1], ssj1[s], md);
            stems[i]  += vrna_E_exterior_stem(type,
                                              (a2s[s][i + 1] > 1) ? S5[s][i + 1] : -1,
                                              s3j1[s],
                                              P);
          }
        }
      }

      if (sc_spl_stem)
        for (i = j - 1; i > 1; i--)
          if (stems[i] != INF)
            stems[i] += sc_spl_stem(j, i - 1, i + 1, sc_wrapper);

      stems[1] = INF;
      if (3 < j) {
        ij = indx[j - 1] + 2;

        if ((c[ij] != INF) && (evaluate(1, j, 2, j - 1, VRNA_DECOMP_EXT_STEM, hc_dat_local))) {
          stems[1] = c[ij];
          for (s = 0; s < n_seq; s++) {
            type      = vrna_get_ptype_md(SS[s][2], ssj1[s], md);
            stems[1]  += vrna_E_exterior_stem(type, (a2s[s][2] > 1) ? S5[s][2] : -1, s3j1[s], P);
          }

          if (sc_red_stem)
            stems[1] += sc_red_stem(j, 2, j - 1, sc_wrapper);
        }
      }

      free(s3j1);
      free(ssj1);

      break;
  }

  return stems;
}


PRIVATE INLINE int
add_f5_gquad(vrna_fold_compound_t   *fc,
             unsigned int           j,
             vrna_hc_eval_f         evaluate VRNA_UNUSED,
             struct hc_ext_def_dat  *hc_dat_local VRNA_UNUSED,
             struct sc_f5_dat       *sc_wrapper VRNA_UNUSED)
{
  unsigned int      i;
  int               e, e_gq, *f5;
  vrna_smx_csr(int) *c_gq;

  f5    = fc->matrices->f5;
  c_gq  = fc->matrices->c_gq;
  e     = INF;

  for (i = j - 1; i > 1; i--) {
#ifndef VRNA_DISABLE_C11_FEATURES
    e_gq = vrna_smx_csr_get(c_gq, i, j, INF);
#else
    e_gq = vrna_smx_csr_int_get(c_gq, i, j, INF);
#endif
    if ((f5[i - 1] != INF) &&
        (e_gq != INF))
      e = MIN2(e, f5[i - 1] + e_gq);
  }

#ifndef VRNA_DISABLE_C11_FEATURES
  e_gq = vrna_smx_csr_get(c_gq, 1, j, INF);
#else
  e_gq = vrna_smx_csr_int_get(c_gq, 1, j, INF);
#endif
  if (e_gq != INF)
    e = MIN2(e, e_gq);

  return e;
}


PRIVATE INLINE int
decompose_f5_ext_stem(vrna_fold_compound_t  *fc,
                      unsigned int          j,
                      int                   *stems)
{
  int       e, *f5;

  f5  = fc->matrices->f5;
  e   = INF;

  const int count = j;

  e = vrna_fun_zip_add_min(f5 + 1, stems + 2, count - 2);

  return e;
}
