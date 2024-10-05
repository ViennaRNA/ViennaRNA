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

#ifdef VRNA_WITH_SVM
#include "ViennaRNA/zscore/basic.h" /* for VRNA_ZSCORE_* macros */
#endif


#include "ViennaRNA/constraints/exterior_hc.inc"
#include "ViennaRNA/constraints/exterior_sc.inc"

#ifdef VRNA_WITH_SVM
#include "ViennaRNA/intern/zscore_dat.h"
#endif

#include "ViennaRNA/mfe/exterior.h"

/*
 #################################
 # PRIVATE FUNCTION DECLARATIONS #
 #################################
 */
PRIVATE INLINE int
reduce_f3_up(vrna_fold_compound_t   *fc,
             unsigned int           i,
             vrna_hc_eval_f         evaluate,
             struct hc_ext_def_dat  *hc_dat_local,
             struct sc_f3_dat       *sc_wrapper);


PRIVATE INLINE int *
f3_get_stem_contributions_d0(vrna_fold_compound_t   *fc,
                             unsigned int           i,
                             vrna_hc_eval_f         evaluate,
                             struct hc_ext_def_dat  *hc_dat_local,
                             struct sc_f3_dat       *sc_wrapper);


PRIVATE INLINE int *
f3_get_stem_contributions_d2(vrna_fold_compound_t   *fc,
                             unsigned int           i,
                             vrna_hc_eval_f         evaluate,
                             struct hc_ext_def_dat  *hc_dat_local,
                             struct sc_f3_dat       *sc_wrapper);


PRIVATE INLINE int *
f3_get_stem_contributions_d3(vrna_fold_compound_t   *fc,
                             unsigned int           i,
                             vrna_hc_eval_f         evaluate,
                             struct hc_ext_def_dat  *hc_dat_local,
                             struct sc_f3_dat       *sc_wrapper);


PRIVATE INLINE int *
f3_get_stem_contributions_d5(vrna_fold_compound_t   *fc,
                             unsigned int           i,
                             vrna_hc_eval_f         evaluate,
                             struct hc_ext_def_dat  *hc_dat_local,
                             struct sc_f3_dat       *sc_wrapper);


PRIVATE INLINE int *
f3_get_stem_contributions_d53(vrna_fold_compound_t  *fc,
                              unsigned int          i,
                              vrna_hc_eval_f        evaluate,
                              struct hc_ext_def_dat *hc_dat_local,
                              struct sc_f3_dat      *sc_wrapper);


PRIVATE INLINE int
decompose_f3_ext_stem(vrna_fold_compound_t  *fc,
                      unsigned int          i,
                      unsigned int          max_j,
                      int                   *stems);


PRIVATE INLINE int
decompose_f3_ext_stem_d0(vrna_fold_compound_t   *fc,
                         unsigned int           i,
                         vrna_hc_eval_f         evaluate,
                         struct hc_ext_def_dat  *hc_dat_local,
                         struct sc_f3_dat       *sc_wrapper);


PRIVATE INLINE int
decompose_f3_ext_stem_d2(vrna_fold_compound_t   *fc,
                         unsigned int           i,
                         vrna_hc_eval_f         evaluate,
                         struct hc_ext_def_dat  *hc_dat_local,
                         struct sc_f3_dat       *sc_wrapper);


PRIVATE INLINE int
decompose_f3_ext_stem_d1(vrna_fold_compound_t   *fc,
                         unsigned int           i,
                         vrna_hc_eval_f         evaluate,
                         struct hc_ext_def_dat  *hc_dat_local,
                         struct sc_f3_dat       *sc_wrapper);


PRIVATE INLINE int
add_f3_gquad(vrna_fold_compound_t   *fc,
             unsigned int           i,
             vrna_hc_eval_f         evaluate,
             struct hc_ext_def_dat  *hc_dat_local,
             struct sc_f3_dat       *sc_wrapper);


/*
 #################################
 # BEGIN OF FUNCTION DEFINITIONS #
 #################################
 */
PUBLIC int
vrna_mfe_exterior_f3(vrna_fold_compound_t *fc,
                     unsigned int         i)
{
  if (fc) {
    unsigned int          dangle_model, with_gquad;
    int                   e, en;
    vrna_param_t          *P;
    vrna_md_t             *md;
    vrna_hc_eval_f        evaluate;
    struct hc_ext_def_dat hc_dat_local;
    struct sc_f3_dat      sc_wrapper;

    e = INF;

    P             = fc->params;
    md            = &(P->model_details);
    dangle_model  = md->dangles;
    with_gquad    = md->gquad;
    evaluate      = prepare_hc_ext_def_window(fc, &hc_dat_local);

    init_sc_f3(fc, i, &sc_wrapper);

    /* first case: i stays unpaired */
    e = reduce_f3_up(fc, i, evaluate, &hc_dat_local, &sc_wrapper);

    /* decompose into stem followed by exterior loop part */
    switch (dangle_model) {
      case 0:
        en  = decompose_f3_ext_stem_d0(fc, i, evaluate, &hc_dat_local, &sc_wrapper);
        e   = MIN2(e, en);
        break;

      case 2:
        en  = decompose_f3_ext_stem_d2(fc, i, evaluate, &hc_dat_local, &sc_wrapper);
        e   = MIN2(e, en);
        break;

      default:
        en  = decompose_f3_ext_stem_d1(fc, i, evaluate, &hc_dat_local, &sc_wrapper);
        e   = MIN2(e, en);
        break;
    }

    if (with_gquad) {
      en  = add_f3_gquad(fc, i, evaluate, &hc_dat_local, &sc_wrapper);
      e   = MIN2(e, en);
    }

    free_sc_f3(&sc_wrapper);

    return e;
  }

  return INF;
}


/*
 #####################################
 # BEGIN OF STATIC HELPER FUNCTIONS  #
 #####################################
 */
PRIVATE INLINE int
decompose_f3_ext_stem_d0(vrna_fold_compound_t   *fc,
                         unsigned int           i,
                         vrna_hc_eval_f         evaluate,
                         struct hc_ext_def_dat  *hc_dat_local,
                         struct sc_f3_dat       *sc_wrapper)
{
  unsigned int  maxdist, length;
  int           e, *stems;

  length  = fc->length;
  maxdist = fc->window_size;

  stems = f3_get_stem_contributions_d0(fc, i, evaluate, hc_dat_local, sc_wrapper);

  stems -= i;

  /* 1st case, actual decompostion */
  e = decompose_f3_ext_stem(fc, i, MIN2(length - 1, i + maxdist), stems);

  /* 2nd case, reduce to single stem */
  if (length <= i + maxdist)
    e = MIN2(e, stems[length]);

  stems += i;
  free(stems);

  return e;
}


PRIVATE INLINE int
decompose_f3_ext_stem_d2(vrna_fold_compound_t   *fc,
                         unsigned int           i,
                         vrna_hc_eval_f         evaluate,
                         struct hc_ext_def_dat  *hc_dat_local,
                         struct sc_f3_dat       *sc_wrapper)
{
  unsigned int  maxdist, length;
  int           e, *stems;

  length  = fc->length;
  maxdist = fc->window_size;
  stems   = f3_get_stem_contributions_d2(fc, i, evaluate, hc_dat_local, sc_wrapper);

  stems -= i;

  /* 1st case, actual decompostion */
  e = decompose_f3_ext_stem(fc, i, MIN2(length - 1, i + maxdist), stems);

  /* 2nd case, reduce to single stem */
  if (length <= i + maxdist)
    e = MIN2(e, stems[length]);

  stems += i;
  free(stems);

  return e;
}


PRIVATE INLINE int
decompose_f3_ext_stem_d1(vrna_fold_compound_t   *fc,
                         unsigned int           i,
                         vrna_hc_eval_f         evaluate,
                         struct hc_ext_def_dat  *hc_dat_local,
                         struct sc_f3_dat       *sc_wrapper)
{
  unsigned int  length, maxdist;
  int           e, ee, *stems;

  length  = fc->length;
  maxdist = fc->window_size;
  e       = INF;

  /* A) without dangling end contributions */

  /* 1st case, actual decompostion */
  stems = f3_get_stem_contributions_d0(fc, i, evaluate, hc_dat_local, sc_wrapper);

  stems -= i;

  ee = decompose_f3_ext_stem(fc, i, MIN2(length - 1, i + maxdist), stems);

  /* 2nd case, reduce to single stem */
  if (length <= i + maxdist)
    ee = MIN2(ee, stems[length]);

  stems += i;
  free(stems);

  e = MIN2(e, ee);

  /* B) with dangling end contribution on 3' side of stem */
  stems = f3_get_stem_contributions_d3(fc, i, evaluate, hc_dat_local, sc_wrapper);

  stems -= i;

  /* 1st case, actual decompostion */
  ee = decompose_f3_ext_stem(fc, i, MIN2(length - 1, i + maxdist + 1), stems);

  /* 2nd case, reduce to single stem */
  if (length <= i + maxdist)
    ee = MIN2(ee, stems[length]);

  stems += i;
  free(stems);

  e = MIN2(e, ee);

  /* C) with dangling end contribution on 5' side of stem */
  stems = f3_get_stem_contributions_d5(fc, i, evaluate, hc_dat_local, sc_wrapper);

  stems -= i;

  /* 1st case, actual decompostion */
  ee = decompose_f3_ext_stem(fc, i, MIN2(length - 1, i + maxdist + 1), stems);

  /* 2nd case, reduce to single stem */
  if (length <= i + maxdist)
    ee = MIN2(ee, stems[length]);

  stems += i;
  free(stems);

  e = MIN2(e, ee);

  /* D) with dangling end contribution on both sides of stem */
  stems = f3_get_stem_contributions_d53(fc, i, evaluate, hc_dat_local, sc_wrapper);

  stems -= i;

  /* 1st case, actual decompostion */
  ee = decompose_f3_ext_stem(fc, i, MIN2(length - 1, i + maxdist + 1), stems);

  /* 2nd case, reduce to single stem */
  if (length <= i + maxdist)
    ee = MIN2(ee, stems[length]);

  stems += i;
  free(stems);

  e = MIN2(e, ee);

  return e;
}


/*
 *  extend f3 by adding an unpaired nucleotide or an unstructured domain
 *  to the 5' end
 */
PRIVATE INLINE int
reduce_f3_up(vrna_fold_compound_t   *fc,
             unsigned int           i,
             vrna_hc_eval_f         evaluate,
             struct hc_ext_def_dat  *hc_dat_local,
             struct sc_f3_dat       *sc_wrapper)
{
  unsigned int  u, k, length;
  int           e, en, *f3;
  vrna_ud_t     *domains_up;
  sc_f3_cb      sc_red_ext;

  length      = fc->length;
  f3          = fc->matrices->f3_local;
  domains_up  = fc->domains_up;
  e           = INF;

  sc_red_ext = sc_wrapper->red_ext;

  /* check for 3' extension with one unpaired nucleotide */
  if (f3[i + 1] != INF) {
    if (evaluate(i, length, i + 1, length, VRNA_DECOMP_EXT_EXT, hc_dat_local)) {
      e = f3[i + 1];

      if (sc_red_ext)
        e += sc_red_ext(i, i + 1, length, sc_wrapper);
    }
  }

  if ((domains_up) && (domains_up->energy_cb)) {
    for (k = 0; k < (unsigned int)domains_up->uniq_motif_count; k++) {
      u = domains_up->uniq_motif_size[k];
      if ((i + u - 1 <= length) && (f3[i + u] != INF)) {
        if (evaluate(i, length, i + u - 1, length, VRNA_DECOMP_EXT_EXT, hc_dat_local)) {
          en = f3[i + u] +
               domains_up->energy_cb(fc,
                                     i,
                                     i + u - 1,
                                     VRNA_UNSTRUCTURED_DOMAIN_EXT_LOOP | VRNA_UNSTRUCTURED_DOMAIN_MOTIF,
                                     domains_up->data);

          if (sc_red_ext)
            en += sc_red_ext(i, i + u, length, sc_wrapper);

          e = MIN2(e, en);
        }
      }
    }
  }

  return e;
}


PRIVATE INLINE int *
f3_get_stem_contributions_d0(vrna_fold_compound_t   *fc,
                             unsigned int           i,
                             vrna_hc_eval_f         evaluate,
                             struct hc_ext_def_dat  *hc_dat_local,
                             struct sc_f3_dat       *sc_wrapper)
{
  char            **ptype;
  short           **S, *si;
  unsigned int    s, n_seq, type, length, j, max_j, maxdist;
  int             energy, *c, *stems;
  vrna_param_t    *P;
  vrna_md_t       *md;

#ifdef VRNA_WITH_SVM
  int             zsc_pre_filter;
  vrna_zsc_dat_t  zsc_data;
#endif

  sc_f3_cb        sc_spl_stem;
  sc_f3_cb        sc_red_stem;

  length  = fc->length;
  maxdist = fc->window_size;
  P       = fc->params;
  md      = &(P->model_details);
  c       = fc->matrices->c_local[i];
  c       -= i;
  si      = NULL;
  ptype   = (fc->type == VRNA_FC_TYPE_SINGLE) ? fc->ptype_local : NULL;
  n_seq   = (fc->type == VRNA_FC_TYPE_SINGLE) ? 1 : fc->n_seq;
  S       = (fc->type == VRNA_FC_TYPE_SINGLE) ? NULL : fc->S;
#ifdef VRNA_WITH_SVM
  zsc_data        = fc->zscore_data;
  zsc_pre_filter  = ((zsc_data) &&
                     (zsc_data->filter_on) &&
                     (zsc_data->pre_filter)) ? 1 : 0;
#endif
  stems = (int *)vrna_alloc(sizeof(int) * (maxdist + 6));
  stems -= i;

  sc_spl_stem = sc_wrapper->decomp_stem;
  sc_red_stem = sc_wrapper->red_stem;

  max_j = MIN2(length - 1, i + maxdist);

#ifdef VRNA_WITH_SVM
  /* re-set pointer */
  if (zsc_pre_filter) {
    zsc_data->current_z += zsc_data->current_i;
    /* initialize */
    memset(zsc_data->current_z, 0, sizeof(double) * (maxdist + 2));
    /* move pointer for convenience */
    zsc_data->current_i = i;
    zsc_data->current_z -= zsc_data->current_i;
  }

#endif

  switch (fc->type) {
    case VRNA_FC_TYPE_SINGLE:
      for (j = i + 1; j <= max_j; j++) {
        stems[j] = INF;
        if ((c[j] != INF) &&
            (evaluate(i, length, j, j + 1, VRNA_DECOMP_EXT_STEM_EXT, hc_dat_local))) {
          type      = vrna_get_ptype_window(i, j, ptype);
          stems[j]  = c[j] +
                      vrna_E_exterior_stem(type, -1, -1, P);
        }
      }

#ifdef VRNA_WITH_SVM
      /* if necessary, remove those stems where the z-score threshold is not satisfied */
      if (zsc_pre_filter) {
        for (j = i + 1; j <= max_j; j++) {
          if (stems[j] != INF) {
            zsc_data->current_z[j] = vrna_zsc_compute(fc, i, j, stems[j]);
            if (zsc_data->current_z[j] > zsc_data->min_z)
              stems[j] = INF;
          }
        }
      }

#endif

      break;

    case VRNA_FC_TYPE_COMPARATIVE:
      si = (short *)vrna_alloc(sizeof(short) * n_seq);

      for (s = 0; s < n_seq; s++)
        si[s] = S[s][i];

      for (j = i + 1; j <= max_j; j++) {
        stems[j] = INF;
        if ((c[j] != INF) &&
            (evaluate(i, length, j, j + 1, VRNA_DECOMP_EXT_STEM_EXT, hc_dat_local))) {
          energy = c[j];
          for (s = 0; s < n_seq; s++) {
            type    = vrna_get_ptype_md(si[s], S[s][j], md);
            energy  += vrna_E_exterior_stem(type, -1, -1, P);
          }
          stems[j] = energy;
        }
      }
      break;
  }


  if (sc_spl_stem)
    for (j = i + 1; j <= max_j; j++)
      if (stems[j] != INF)
        stems[j] += sc_spl_stem(i, j, j + 1, sc_wrapper);

  if (length <= i + maxdist) {
    j         = length;
    stems[j]  = INF;

    if ((c[j] != INF) &&
        (evaluate(i, j, i, j, VRNA_DECOMP_EXT_STEM, hc_dat_local))) {
      energy = c[j];

      switch (fc->type) {
        case VRNA_FC_TYPE_SINGLE:
          type    = vrna_get_ptype_window(i, j, ptype);
          energy  += vrna_E_exterior_stem(type, -1, -1, P);

          break;

        case VRNA_FC_TYPE_COMPARATIVE:
          for (s = 0; s < n_seq; s++) {
            type    = vrna_get_ptype_md(si[s], S[s][j], md);
            energy  += vrna_E_exterior_stem(type, -1, -1, P);
          }

          break;
      }

#ifdef VRNA_WITH_SVM
      /* if necessary, remove those stems where the z-score threshold is not satisfied */
      if ((fc->type == VRNA_FC_TYPE_SINGLE) &&
          (zsc_pre_filter) &&
          (energy != INF)) {
        zsc_data->current_z[j] = vrna_zsc_compute(fc, i, j, stems[j]);
        if (zsc_data->current_z[j] > zsc_data->min_z)
          energy = INF;
      }

#endif

      if ((sc_red_stem) &&
          (energy != INF))
        energy += sc_red_stem(i, i, j, sc_wrapper);

      stems[j] = energy;
    }
  } else {
    /*
     * make sure we do not take (i + maxdist + 1) into account when
     * decomposing for odd dangle models
     */
    stems[i + maxdist + 1] = INF;
  }

  free(si);

  stems += i;

  return stems;
}


PRIVATE INLINE int *
f3_get_stem_contributions_d2(vrna_fold_compound_t   *fc,
                             unsigned int           i,
                             vrna_hc_eval_f         evaluate,
                             struct hc_ext_def_dat  *hc_dat_local,
                             struct sc_f3_dat       *sc_wrapper)
{
  char            **ptype;
  short           **S, **S5, **S3, *S1, si1, sj1, *s5i1, *si;
  unsigned int    s, n_seq, type, length, **a2s, j, max_j, maxdist;
  int             energy, *c, *stems;
  vrna_param_t    *P;
  vrna_md_t       *md;

#ifdef VRNA_WITH_SVM
  int             zsc_pre_filter;
  vrna_zsc_dat_t  zsc_data;
#endif

  sc_f3_cb        sc_spl_stem;
  sc_f3_cb        sc_red_stem;

  length  = fc->length;
  maxdist = fc->window_size;
  P       = fc->params;
  md      = &(P->model_details);
  c       = fc->matrices->c_local[i];
  c       -= i;
#ifdef VRNA_WITH_SVM
  zsc_data        = fc->zscore_data;
  zsc_pre_filter  = ((zsc_data) &&
                     (zsc_data->filter_on) &&
                     (zsc_data->pre_filter)) ? 1 : 0;
#endif

  stems = (int *)vrna_alloc(sizeof(int) * (maxdist + 6));
  stems -= i;

  sc_spl_stem = sc_wrapper->decomp_stem;
  sc_red_stem = sc_wrapper->red_stem;

#ifdef VRNA_WITH_SVM
  /* re-set pointer */
  if (zsc_pre_filter) {
    zsc_data->current_z += zsc_data->current_i;
    /* initialize */
    memset(zsc_data->current_z, 0, sizeof(double) * (maxdist + 2));
    /* move pointer for convenience */
    zsc_data->current_i = i;
    zsc_data->current_z -= zsc_data->current_i;
  }

#endif

  switch (fc->type) {
    case VRNA_FC_TYPE_SINGLE:
      ptype = fc->ptype_local;
      S1    = fc->sequence_encoding;
      si1   = (i > 1) ? S1[i - 1] : -1;
      max_j = MIN2(length - 1, i + maxdist);

      for (j = i + 1; j <= max_j; j++) {
        stems[j] = INF;
        if ((c[j] != INF) &&
            (evaluate(i, length, j, j + 1, VRNA_DECOMP_EXT_STEM_EXT, hc_dat_local))) {
          type      = vrna_get_ptype_window(i, j, ptype);
          sj1       = S1[j + 1];
          stems[j]  = c[j] +
                      vrna_E_exterior_stem(type, si1, sj1, P);
        }
      }

#ifdef VRNA_WITH_SVM
      /* if necessary, remove those stems where the z-score threshold is not satisfied */
      if (zsc_pre_filter) {
        for (j = i + 1; j <= max_j; j++) {
          if (stems[j] != INF) {
            zsc_data->current_z[j] = vrna_zsc_compute(fc, i, j, stems[j]);
            if (zsc_data->current_z[j] > zsc_data->min_z)
              stems[j] = INF;
          }
        }
      }

#endif

      if (sc_spl_stem)
        for (j = i + 1; j <= max_j; j++)
          if (stems[j] != INF)
            stems[j] += sc_spl_stem(i, j, j + 1, sc_wrapper);

      if (length <= (i + maxdist)) {
        j         = length;
        stems[j]  = INF;

        if ((c[j] != INF) && (evaluate(i, j, i, j, VRNA_DECOMP_EXT_STEM, hc_dat_local))) {
          type = vrna_get_ptype_window(i, j, ptype);

          stems[j] = c[j] +
                     vrna_E_exterior_stem(type, si1, -1, P);

#ifdef VRNA_WITH_SVM
          if ((zsc_pre_filter) &&
              (stems[j] != INF)) {
            zsc_data->current_z[j] = vrna_zsc_compute(fc, i, j, stems[j]);
            if (zsc_data->current_z[j] > zsc_data->min_z)
              stems[j] = INF;
          }

#endif

          if ((sc_red_stem) &&
              (stems[j] != INF))
            stems[j] += sc_red_stem(i, i, j, sc_wrapper);
        }
      }

      break;

    case VRNA_FC_TYPE_COMPARATIVE:
      n_seq = fc->n_seq;
      S     = fc->S;
      S5    = fc->S5;
      S3    = fc->S3;
      a2s   = fc->a2s;
      max_j = MIN2(length - 1, i + maxdist);

      s5i1  = (short *)vrna_alloc(sizeof(short) * n_seq);
      si    = (short *)vrna_alloc(sizeof(short) * n_seq);
      for (s = 0; s < n_seq; s++) {
        s5i1[s] = (a2s[s][i] > 1) ? S5[s][i] : -1;
        si[s]   = S[s][i];
      }

      for (j = i + 1; j <= max_j; j++) {
        stems[j] = INF;
        if ((c[j] != INF) &&
            (evaluate(i, length, j, j + 1, VRNA_DECOMP_EXT_STEM_EXT, hc_dat_local))) {
          energy = c[j];
          for (s = 0; s < n_seq; s++) {
            type    = vrna_get_ptype_md(si[s], S[s][j], md);
            sj1     = (a2s[s][j] < a2s[s][length]) ? S3[s][j] : -1;
            energy  += vrna_E_exterior_stem(type, s5i1[s], sj1, P);
          }
          stems[j] = energy;
        }
      }

      if (sc_spl_stem)
        for (j = i + 1; j <= max_j; j++)
          if (stems[j] != INF)
            stems[j] += sc_spl_stem(i, j, j + 1, sc_wrapper);

      if (length <= (i + maxdist)) {
        j         = length;
        stems[j]  = INF;

        if ((c[j] != INF) && (evaluate(i, j, i, j, VRNA_DECOMP_EXT_STEM, hc_dat_local))) {
          energy = c[j];
          for (s = 0; s < n_seq; s++) {
            type    = vrna_get_ptype_md(si[s], S[s][j], md);
            energy  += vrna_E_exterior_stem(type, s5i1[s], -1, P);
          }

          if (sc_red_stem)
            energy += sc_red_stem(i, i, j, sc_wrapper);

          stems[j] = energy;
        }
      }

      free(s5i1);
      free(si);

      break;
  }

  stems += i;

  return stems;
}


PRIVATE INLINE int *
f3_get_stem_contributions_d3(vrna_fold_compound_t   *fc,
                             unsigned int           i,
                             vrna_hc_eval_f         evaluate,
                             struct hc_ext_def_dat  *hc_dat_local,
                             struct sc_f3_dat       *sc_wrapper)
{
  char            **ptype;
  short           *S1, **S, **S3, sj1, *si;
  unsigned int    s, n_seq, **a2s, type, j, max_j, length, maxdist;
  int             energy, *c, *stems;
  vrna_param_t    *P;
  vrna_md_t       *md;

#ifdef VRNA_WITH_SVM
  int             zsc_pre_filter;
  vrna_zsc_dat_t  zsc_data;
#endif

  sc_f3_cb        sc_spl_stem;
  sc_f3_cb        sc_red_stem;

  length  = fc->length;
  maxdist = fc->window_size;
  P       = fc->params;
  md      = &(P->model_details);
  c       = fc->matrices->c_local[i];
  c       -= i;
#ifdef VRNA_WITH_SVM
  zsc_data        = fc->zscore_data;
  zsc_pre_filter  = ((zsc_data) &&
                     (zsc_data->filter_on) &&
                     (zsc_data->pre_filter)) ? 1 : 0;
#endif

  stems = (int *)vrna_alloc(sizeof(int) * (maxdist + 6));
  stems -= i;

  sc_spl_stem = sc_wrapper->decomp_stem;
  sc_red_stem = sc_wrapper->red_stem;

#ifdef VRNA_WITH_SVM
  /* re-set pointer */
  if (zsc_pre_filter) {
    zsc_data->current_z += zsc_data->current_i;
    /* initialize */
    memset(zsc_data->current_z, 0, sizeof(double) * (maxdist + 2));
    /* move pointer for convenience */
    zsc_data->current_i = i;
    zsc_data->current_z -= zsc_data->current_i;
  }

#endif

  switch (fc->type) {
    case VRNA_FC_TYPE_SINGLE:
      S1    = fc->sequence_encoding;
      ptype = fc->ptype_local;
      max_j = MIN2(length - 1, i + maxdist + 1);

      for (j = i + 1; j <= max_j; j++) {
        stems[j] = INF;
        if ((c[j - 1] != INF) &&
            (evaluate(i, length, j - 1, j + 1, VRNA_DECOMP_EXT_STEM_EXT, hc_dat_local))) {
          type      = vrna_get_ptype_window(i, j - 1, ptype);
          stems[j]  = c[j - 1] +
                      vrna_E_exterior_stem(type, -1, S1[j], P);
        }
      }

#ifdef VRNA_WITH_SVM
      /* if necessary, remove those stems where the z-score threshold is not satisfied */
      if (zsc_pre_filter) {
        for (j = i + 1; j <= max_j; j++) {
          if (stems[j] != INF) {
            zsc_data->current_z[j] = vrna_zsc_compute(fc, i, j, stems[j]);
            if (zsc_data->current_z[j] > zsc_data->min_z)
              stems[j] = INF;
          }
        }
      }

#endif

      if (sc_spl_stem)
        for (j = i + 1; j <= max_j; j++)
          if (stems[j] != INF)
            stems[j] += sc_spl_stem(i, j - 1, j + 1, sc_wrapper);

      if (length <= i + maxdist) {
        j = length;

        if ((c[j - 1] != INF) &&
            (evaluate(i, j, i, j - 1, VRNA_DECOMP_EXT_STEM, hc_dat_local))) {
          type      = vrna_get_ptype_window(i, j - 1, ptype);
          stems[j]  = c[j - 1] +
                      vrna_E_exterior_stem(type, -1, S1[j], P);

#ifdef VRNA_WITH_SVM
          /* if necessary, remove those stems where the z-score threshold is not satisfied */
          if ((zsc_pre_filter) &&
              (stems[j] != INF)) {
            zsc_data->current_z[j] = vrna_zsc_compute(fc, i, j, stems[j]);
            if (zsc_data->current_z[j] > zsc_data->min_z)
              stems[j] = INF;
          }

#endif
          if ((sc_red_stem) &&
              (stems[j] != INF))
            stems[j] += sc_red_stem(i, i, j - 1, sc_wrapper);
        }
      }

      break;

    case VRNA_FC_TYPE_COMPARATIVE:
      n_seq = fc->n_seq;
      S     = fc->S;
      S3    = fc->S3;
      a2s   = fc->a2s;
      max_j = MIN2(length - 1, i + maxdist + 1);

      si = (short *)vrna_alloc(sizeof(short) * n_seq);
      for (s = 0; s < n_seq; s++)
        si[s] = S[s][i];

      for (j = i + 1; j <= max_j; j++) {
        stems[j] = INF;
        if ((c[j - 1] != INF) &&
            (evaluate(i, length, j - 1, j + 1, VRNA_DECOMP_EXT_STEM_EXT, hc_dat_local))) {
          energy = c[j - 1];
          for (s = 0; s < n_seq; s++) {
            type    = vrna_get_ptype_md(si[s], S[s][j - 1], md);
            sj1     = (a2s[s][j - 1] < a2s[s][length]) ? S3[s][j - 1] : -1;
            energy  += vrna_E_exterior_stem(type, -1, sj1, P);
          }
          stems[j] = energy;
        }
      }

      if (sc_spl_stem)
        for (j = i + 1; j <= max_j; j++)
          if (stems[j] != INF)
            stems[j] += sc_spl_stem(i, j - 1, j + 1, sc_wrapper);

      if (length <= i + maxdist) {
        j = length;

        if ((c[j - 1] != INF) &&
            (evaluate(i, j, i, j - 1, VRNA_DECOMP_EXT_STEM, hc_dat_local))) {
          energy = c[j - 1];
          for (s = 0; s < n_seq; s++) {
            type    = vrna_get_ptype_md(si[s], S[s][j - 1], md);
            sj1     = (a2s[s][j - 1] < a2s[s][length]) ? S3[s][j - 1] : -1;
            energy  += vrna_E_exterior_stem(type, -1, sj1, P);
          }

          if (sc_red_stem)
            energy += sc_red_stem(i, i, j - 1, sc_wrapper);

          stems[j] = energy;
        }
      }

      free(si);

      break;
  }

  stems += i;

  return stems;
}


PRIVATE INLINE int *
f3_get_stem_contributions_d5(vrna_fold_compound_t   *fc,
                             unsigned int           i,
                             vrna_hc_eval_f         evaluate,
                             struct hc_ext_def_dat  *hc_dat_local,
                             struct sc_f3_dat       *sc_wrapper)
{
  char            **ptype;
  short           *S1, **S, **S5, *s5i1, si, *si1;
  unsigned int    s, n_seq, **a2s, type, j, max_j, length, maxdist;
  int             energy, *c, *stems;
  vrna_param_t    *P;
  vrna_md_t       *md;

#ifdef VRNA_WITH_SVM
  int             zsc_pre_filter;
  vrna_zsc_dat_t  zsc_data;
#endif

  sc_f3_cb        sc_spl_stem;
  sc_f3_cb        sc_red_stem;

  length  = fc->length;
  maxdist = fc->window_size;
  P       = fc->params;
  md      = &(P->model_details);
  c       = fc->matrices->c_local[i + 1];
  c       -= i + 1;
#ifdef VRNA_WITH_SVM
  zsc_data        = fc->zscore_data;
  zsc_pre_filter  = ((zsc_data) &&
                     (zsc_data->filter_on) &&
                     (zsc_data->pre_filter)) ? 1 : 0;
#endif

  stems = (int *)vrna_alloc(sizeof(int) * (maxdist + 6));
  stems -= i;

  sc_spl_stem = sc_wrapper->decomp_stem1;
  sc_red_stem = sc_wrapper->red_stem;

#ifdef VRNA_WITH_SVM
  /* re-set pointer */
  if (zsc_pre_filter) {
    zsc_data->current_z += zsc_data->current_i;
    /* initialize */
    memset(zsc_data->current_z, 0, sizeof(double) * (maxdist + 2));
    /* move pointer for convenience */
    zsc_data->current_i = i;
    zsc_data->current_z -= zsc_data->current_i;
  }

#endif

  switch (fc->type) {
    case VRNA_FC_TYPE_SINGLE:
      S1    = fc->sequence_encoding;
      ptype = fc->ptype_local;
      si    = S1[i];
      max_j = MIN2(length - 1, i + maxdist + 1);

      for (j = i + 1; j <= max_j; j++) {
        stems[j] = INF;
        if ((c[j] != INF) &&
            (evaluate(i, length, j, j + 1, VRNA_DECOMP_EXT_STEM_EXT1, hc_dat_local))) {
          type      = vrna_get_ptype_window(i + 1, j, ptype);
          stems[j]  = c[j] +
                      vrna_E_exterior_stem(type, si, -1, P);
        }
      }

#ifdef VRNA_WITH_SVM
      /* if necessary, remove those stems where the z-score threshold is not satisfied */
      if (zsc_pre_filter) {
        for (j = i + 1; j <= max_j; j++) {
          if (stems[j] != INF) {
            zsc_data->current_z[j] = vrna_zsc_compute(fc, i, j, stems[j]);
            if (zsc_data->current_z[j] > zsc_data->min_z)
              stems[j] = INF;
          }
        }
      }

#endif

      if (sc_spl_stem)
        for (j = i + 1; j <= max_j; j++)
          if (stems[j] != INF)
            stems[j] += sc_spl_stem(i, j, j + 1, sc_wrapper);

      if (length <= i + maxdist) {
        j         = length;
        stems[j]  = INF;

        if ((c[j] != INF) &&
            (evaluate(i, j, i + 1, j, VRNA_DECOMP_EXT_STEM, hc_dat_local))) {
          type      = vrna_get_ptype_window(i + 1, j, ptype);
          stems[j]  = c[j] +
                      vrna_E_exterior_stem(type, si, -1, P);

#ifdef VRNA_WITH_SVM
          /* if necessary, remove those stems where the z-score threshold is not satisfied */
          if ((zsc_pre_filter) &&
              (stems[j] != INF)) {
            zsc_data->current_z[j] = vrna_zsc_compute(fc, i, j, stems[j]);
            if (zsc_data->current_z[j] > zsc_data->min_z)
              stems[j] = INF;
          }

#endif

          if ((sc_red_stem) &&
              (stems[j] != INF))
            stems[j] += sc_red_stem(i, i + 1, j, sc_wrapper);
        }
      }

      break;

    case VRNA_FC_TYPE_COMPARATIVE:
      n_seq = fc->n_seq;
      S     = fc->S;
      S5    = fc->S5;
      a2s   = fc->a2s;
      max_j = MIN2(length - 1, i + maxdist + 1);

      /* pre-compute S5[s][i] */
      s5i1  = (short *)vrna_alloc(sizeof(short) * n_seq);
      si1   = (short *)vrna_alloc(sizeof(short) * n_seq);
      for (s = 0; s < n_seq; s++) {
        s5i1[s] = (a2s[s][i + 1] > 1) ? S5[s][i + 1] : -1;
        si1[s]  = S[s][i + 1];
      }

      for (j = i + 1; j <= max_j; j++) {
        stems[j] = INF;
        if ((c[j] != INF) &&
            (evaluate(i, length, j, j + 1, VRNA_DECOMP_EXT_STEM_EXT1, hc_dat_local))) {
          energy = c[j];
          for (s = 0; s < n_seq; s++) {
            type    = vrna_get_ptype_md(si1[s], S[s][j], md);
            energy  += vrna_E_exterior_stem(type, s5i1[s], -1, P);
          }
          stems[j] = energy;
        }
      }

      if (sc_spl_stem)
        for (j = i + 1; j <= max_j; j++)
          if (stems[j] != INF)
            stems[j] += sc_spl_stem(i, j, j + 1, sc_wrapper);

      if (length <= i + maxdist) {
        j         = length;
        stems[j]  = INF;

        if ((c[j] != INF) &&
            (evaluate(i, j, i + 1, j, VRNA_DECOMP_EXT_STEM, hc_dat_local))) {
          energy = c[j];
          for (s = 0; s < n_seq; s++) {
            type    = vrna_get_ptype_md(si1[s], S[s][j], md);
            energy  += vrna_E_exterior_stem(type, s5i1[s], -1, P);
          }

          if (sc_red_stem)
            energy += sc_red_stem(i, i + 1, j, sc_wrapper);

          stems[j] = energy;
        }
      }

      free(s5i1);
      free(si1);

      break;
  }

  stems += i;

  return stems;
}


PRIVATE INLINE int *
f3_get_stem_contributions_d53(vrna_fold_compound_t  *fc,
                              unsigned int          i,
                              vrna_hc_eval_f        evaluate,
                              struct hc_ext_def_dat *hc_dat_local,
                              struct sc_f3_dat      *sc_wrapper)
{
  char            **ptype;
  short           *S1, **S, **S5, **S3, *s5i1, si1, sj1, *ssi1;
  unsigned int    s, n_seq, **a2s, type, j, max_j, length, maxdist;
  int             energy, *c, *stems;
  vrna_param_t    *P;
  vrna_md_t       *md;

#ifdef VRNA_WITH_SVM
  int             zsc_pre_filter;
  vrna_zsc_dat_t  zsc_data;
#endif

  sc_f3_cb        sc_spl_stem;
  sc_f3_cb        sc_red_stem;

  length  = fc->length;
  maxdist = fc->window_size;
  P       = fc->params;
  md      = &(P->model_details);
  c       = fc->matrices->c_local[i + 1];
  c       -= i + 1;
#ifdef VRNA_WITH_SVM
  zsc_data        = fc->zscore_data;
  zsc_pre_filter  = ((zsc_data) &&
                     (zsc_data->filter_on) &&
                     (zsc_data->pre_filter)) ? 1 : 0;
#endif

  stems = (int *)vrna_alloc(sizeof(int) * (maxdist + 6));
  stems -= i;

  sc_spl_stem = sc_wrapper->decomp_stem1;
  sc_red_stem = sc_wrapper->red_stem;

#ifdef VRNA_WITH_SVM
  /* re-set pointer */
  if (zsc_pre_filter) {
    zsc_data->current_z += zsc_data->current_i;
    /* initialize */
    memset(zsc_data->current_z, 0, sizeof(double) * (maxdist + 2));
    /* move pointer for convenience */
    zsc_data->current_i = i;
    zsc_data->current_z -= zsc_data->current_i;
  }

#endif

  switch (fc->type) {
    case VRNA_FC_TYPE_SINGLE:
      S1    = fc->sequence_encoding;
      ptype = fc->ptype_local;
      si1   = S1[i];
      max_j = MIN2(length - 1, i + maxdist + 1);

      for (j = i + 1; j <= max_j; j++) {
        stems[j] = INF;
        if ((c[j - 1] != INF) &&
            (evaluate(i, length, j - 1, j + 1, VRNA_DECOMP_EXT_STEM_EXT1, hc_dat_local))) {
          type      = vrna_get_ptype_window(i + 1, j - 1, ptype);
          stems[j]  = c[j - 1] +
                      vrna_E_exterior_stem(type, si1, S1[j], P);
        }
      }

#ifdef VRNA_WITH_SVM
      /* if necessary, remove those stems where the z-score threshold is not satisfied */
      if (zsc_pre_filter) {
        for (j = i + 1; j <= max_j; j++) {
          if (stems[j] != INF) {
            zsc_data->current_z[j] = vrna_zsc_compute(fc, i, j, stems[j]);
            if (zsc_data->current_z[j] > zsc_data->min_z)
              stems[j] = INF;
          }
        }
      }

#endif

      if (sc_spl_stem)
        for (j = i + 1; j <= max_j; j++)
          if (stems[j] != INF)
            stems[j] += sc_spl_stem(i, j - 1, j + 1, sc_wrapper);

      if (length <= i + maxdist) {
        j = length;
        if ((c[j - 1] != INF) &&
            (evaluate(i, length, i + 1, j - 1, VRNA_DECOMP_EXT_STEM, hc_dat_local))) {
          type      = vrna_get_ptype_window(i + 1, j - 1, ptype);
          stems[j]  = c[j - 1] +
                      vrna_E_exterior_stem(type, si1, S1[j], P);

#ifdef VRNA_WITH_SVM
          /* if necessary, remove those stems where the z-score threshold is not satisfied */
          if ((zsc_pre_filter) &&
              (stems[j] != INF)) {
            zsc_data->current_z[j] = vrna_zsc_compute(fc, i, j, stems[j]);
            if (zsc_data->current_z[j] > zsc_data->min_z)
              stems[j] = INF;
          }

#endif

          if ((sc_red_stem) &&
              (stems[j] != INF))
            stems[j] += sc_red_stem(i, i + 1, j - 1, sc_wrapper);
        }
      }

      break;

    case VRNA_FC_TYPE_COMPARATIVE:
      n_seq = fc->n_seq;
      S     = fc->S;
      S5    = fc->S5;
      S3    = fc->S3;
      a2s   = fc->a2s;
      max_j = MIN2(length - 1, i + maxdist + 1);

      s5i1  = (short *)vrna_alloc(sizeof(short) * n_seq);
      ssi1  = (short *)vrna_alloc(sizeof(short) * n_seq);
      for (s = 0; s < n_seq; s++) {
        s5i1[s] = (a2s[s][i + 1] > 1) ? S5[s][i + 1] : -1;
        ssi1[s] = S[s][i + 1];
      }

      for (j = i + 1; j <= max_j; j++) {
        stems[j] = INF;
        if ((c[j - 1] != INF) &&
            (evaluate(i, length, j - 1, j + 1, VRNA_DECOMP_EXT_STEM_EXT1, hc_dat_local))) {
          energy = c[j - 1];
          for (s = 0; s < n_seq; s++) {
            type    = vrna_get_ptype_md(ssi1[s], S[s][j - 1], md);
            sj1     = (a2s[s][j - 1] < a2s[s][length]) ? S3[s][j - 1] : -1;
            energy  += vrna_E_exterior_stem(type, s5i1[s], sj1, P);
          }
          stems[j] = energy;
        }
      }

      if (sc_spl_stem)
        for (j = i + 1; j <= max_j; j++)
          if (stems[j] != INF)
            stems[j] += sc_spl_stem(i, j - 1, j + 1, sc_wrapper);

      if (length <= i + maxdist) {
        j = length;
        if ((c[j - 1] != INF) &&
            (evaluate(i, length, i + 1, j - 1, VRNA_DECOMP_EXT_STEM, hc_dat_local))) {
          energy = c[j - 1];
          for (s = 0; s < n_seq; s++) {
            type    = vrna_get_ptype_md(ssi1[s], S[s][j - 1], md);
            sj1     = (a2s[s][j - 1] < a2s[s][length]) ? S3[s][j - 1] : -1;
            energy  += vrna_E_exterior_stem(type, s5i1[s], sj1, P);
          }

          if (sc_red_stem)
            energy += sc_red_stem(i, i + 1, j - 1, sc_wrapper);

          stems[j] = energy;
        }
      }

      free(s5i1);
      free(ssi1);

      break;
  }

  stems += i;

  return stems;
}


PRIVATE INLINE int
add_f3_gquad(vrna_fold_compound_t   *fc,
             unsigned int           i,
             vrna_hc_eval_f         evaluate VRNA_UNUSED,
             struct hc_ext_def_dat  *hc_dat_local VRNA_UNUSED,
             struct sc_f3_dat       *sc_wrapper VRNA_UNUSED)
{
  unsigned int  j, length, maxdist;
  int           e, *f3, *ggg;

  length  = fc->length;
  maxdist = fc->window_size;
  f3      = fc->matrices->f3_local;
  ggg     = fc->matrices->ggg_local[i];
  e       = INF;

  for (j = i + 1; (j < length) && (j <= i + maxdist); j++)
    if ((f3[j + 1] != INF) && (ggg[j - i] != INF))
      e = MIN2(e, f3[j + 1] + ggg[j - i]);

  if (length <= i + maxdist)
    e = MIN2(e, ggg[length - i]);

  return e;
}


PRIVATE INLINE int
decompose_f3_ext_stem(vrna_fold_compound_t  *fc,
                      unsigned int          i,
                      unsigned int          max_j,
                      int                   *stems)
{
  unsigned int  count;
  int           e, *f3;

  e = INF;
  if (max_j >= i) {
    f3    = fc->matrices->f3_local;
    count = max_j - i;

    /* modular decomposition */
    e = vrna_fun_zip_add_min(stems + i + 1, f3 + i + 2, count);
  }

  return e;
}
