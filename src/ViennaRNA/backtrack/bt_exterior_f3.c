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
#include "ViennaRNA/eval/exterior.h"
#include "ViennaRNA/backtrack/gquad.h"
#include "ViennaRNA/backtrack/exterior.h"

#include "ViennaRNA/intern/grammar_dat.h"

#ifdef __GNUC__
# define INLINE inline
#else
# define INLINE
#endif

#ifdef VRNA_WITH_SVM
#include "ViennaRNA/intern/zscore_dat.h"
#endif

#include "ViennaRNA/constraints/exterior_hc.inc"
#include "ViennaRNA/constraints/exterior_sc.inc"

/*
 #################################
 # PRIVATE FUNCTION DECLARATIONS #
 #################################
 */
PRIVATE unsigned int
bt_ext_loop_f3(vrna_fold_compound_t *fc,
               unsigned int         i,
               unsigned int         j,
               vrna_bps_t           bp_stack,
               vrna_bts_t           bt_stack);


PRIVATE unsigned int
bt_ext_loop_f3_comparative(vrna_fold_compound_t *fc,
                           unsigned int         i,
                           unsigned int         j,
                           vrna_bps_t           bp_stack,
                           vrna_bts_t           bt_stack);


PRIVATE unsigned int
bt_ext_loop_f3_pp(vrna_fold_compound_t  *fc,
                  unsigned int          *i,
                  unsigned int          maxj);


PRIVATE unsigned int
bt_ext_loop_f3_pp_comparative(vrna_fold_compound_t  *fc,
                              unsigned int          *i,
                              unsigned int          maxj);


/*
 #################################
 # BEGIN OF FUNCTION DEFINITIONS #
 #################################
 */
PUBLIC unsigned int
vrna_bt_exterior_f3(vrna_fold_compound_t  *fc,
                    unsigned int          i,
                    unsigned int          j,
                    vrna_bps_t            bp_stack,
                    vrna_bts_t            bt_stack)
{
  if ((fc) &&
      (bp_stack) &&
      (bt_stack)) {
    switch (fc->type) {
      case VRNA_FC_TYPE_SINGLE:
        return bt_ext_loop_f3(fc, i, j, bp_stack, bt_stack);
        break;

      case VRNA_FC_TYPE_COMPARATIVE:
        return bt_ext_loop_f3_comparative(fc, i, j, bp_stack, bt_stack);
        break;
    }
  }

  return 0;
}


PUBLIC unsigned int
vrna_bt_exterior_f3_pp(vrna_fold_compound_t *fc,
                       unsigned int         *i,
                       unsigned int         maxj)
{
  if (fc) {
    switch (fc->type) {
      case VRNA_FC_TYPE_SINGLE:
        return bt_ext_loop_f3_pp(fc, i, maxj);
        break;

      case VRNA_FC_TYPE_COMPARATIVE:
        return bt_ext_loop_f3_pp_comparative(fc, i, maxj);
        break;
    }
  }

  return 0;
}


PRIVATE unsigned int
bt_ext_loop_f3(vrna_fold_compound_t *fc,
               unsigned int         i,
               unsigned int         j,
               vrna_bps_t           bp_stack VRNA_UNUSED,
               vrna_bts_t           bt_stack)
{
  char                  **ptype;
  short                 mm5, mm3, *S1;
  unsigned int          type, length, ii, u, dangle_model, with_gquad;
  int                   fij, fj, *f3, **c, **ggg, en;
  vrna_param_t          *P;
  vrna_md_t             *md;
  vrna_sc_t             *sc;
  vrna_hc_eval_f        evaluate;
  struct hc_ext_def_dat hc_dat_local;

  length        = fc->length;
  P             = fc->params;
  md            = &(P->model_details);
  sc            = fc->sc;
  f3            = fc->matrices->f3_local;
  c             = fc->matrices->c_local;
  ggg           = fc->matrices->ggg_local;
  ptype         = fc->ptype_local;
  S1            = fc->sequence_encoding;
  dangle_model  = md->dangles;
  with_gquad    = md->gquad;
  evaluate      = prepare_hc_ext_def_window(fc, &hc_dat_local);

  ii = i;

  /* nibble off unpaired 5' bases */
  do {
    fij = f3[ii];
    fj  = INF;

    if (evaluate(ii, length, ii + 1, length, VRNA_DECOMP_EXT_EXT, &hc_dat_local)) {
      fj = f3[ii + 1];
      if (sc) {
        if (sc->energy_up)
          fj += sc->energy_up[ii][1];

        if (sc->f)
          fj += sc->f(ii, length, ii + 1, length, VRNA_DECOMP_EXT_EXT, sc->data);
      }
    }

    if (++ii > j)
      break;
  } while (fij == fj);
  ii--;

  if (ii >= j) {
    /* no more pairs */
    return 1;
  }

  /*
   *  must have found a decomposition
   *  i is paired. Find pairing partner
   */
  switch (dangle_model) {
    /* no dangles */
    case 0:
      for (u = j; u > ii; u--) {
        if (with_gquad) {
          if (fij == ggg[ii][u - ii] + f3[u + 1]) {
            vrna_bts_push(bt_stack,
                          ((vrna_sect_t){
                            .i = u + 1,
                            .j = j,
                            .ml = VRNA_MX_FLAG_F3}));

            vrna_bts_push(bt_stack,
                          ((vrna_sect_t){
                            .i = ii,
                            .j = u,
                            .ml = VRNA_MX_FLAG_G}));
            return 1;
          }
        }

        if (evaluate(ii, length, u, u + 1, VRNA_DECOMP_EXT_STEM_EXT, &hc_dat_local)) {
          type = vrna_get_ptype_window(ii, u, ptype);

          en = c[ii][u - ii] +
               vrna_E_exterior_stem(type, -1, -1, P) +
               f3[u + 1];

          if ((sc) && (sc->f))
            en += sc->f(ii, length, u, u + 1, VRNA_DECOMP_EXT_STEM_EXT, sc->data);

          if (fij == en) {
            vrna_bts_push(bt_stack,
                          ((vrna_sect_t){
                            .i = u + 1,
                            .j = j,
                            .ml = VRNA_MX_FLAG_F3}));

            vrna_bts_push(bt_stack,
                          ((vrna_sect_t){
                            .i = ii,
                            .j = u,
                            .ml = VRNA_MX_FLAG_C}));

            return 1;
          }
        }
      }
      break;

    /* dangles on both sides */
    case 2:
      mm5 = (ii > 1) ? S1[ii - 1] : -1;
      for (u = j; u > ii; u--) {
        if (with_gquad) {
          if (fij == ggg[ii][u - ii] + f3[u + 1]) {
            vrna_bts_push(bt_stack,
                          ((vrna_sect_t){
                            .i = u + 1,
                            .j = j,
                            .ml = VRNA_MX_FLAG_F3}));

            vrna_bts_push(bt_stack,
                          ((vrna_sect_t){
                            .i = ii,
                            .j = u,
                            .ml = VRNA_MX_FLAG_G}));
            return 1;
          }
        }

        if (evaluate(ii, length, u, u + 1, VRNA_DECOMP_EXT_STEM_EXT, &hc_dat_local)) {
          mm3   = (u < length) ? S1[u + 1] : -1;
          type  = vrna_get_ptype_window(ii, u, ptype);

          en = c[ii][u - ii] +
               vrna_E_exterior_stem(type, mm5, mm3, P) +
               f3[u + 1];

          if ((sc) && (sc->f))
            en += sc->f(ii, length, u, u + 1, VRNA_DECOMP_EXT_STEM_EXT, sc->data);

          if (fij == en) {
            vrna_bts_push(bt_stack,
                          ((vrna_sect_t){
                            .i = u + 1,
                            .j = j,
                            .ml = VRNA_MX_FLAG_F3}));

            vrna_bts_push(bt_stack,
                          ((vrna_sect_t){
                            .i = ii,
                            .j = u,
                            .ml = VRNA_MX_FLAG_C}));

            return 1;
          }
        }
      }
      break;

    default:
      mm5 = S1[ii];
      for (u = j; u > ii; u--) {
        if (with_gquad) {
          if (fij == ggg[ii][u - ii] + f3[u + 1]) {
            vrna_bts_push(bt_stack,
                          ((vrna_sect_t){
                            .i = u + 1,
                            .j = j,
                            .ml = VRNA_MX_FLAG_F3}));

            vrna_bts_push(bt_stack,
                          ((vrna_sect_t){
                            .i = ii,
                            .j = u,
                            .ml = VRNA_MX_FLAG_G}));
            return 1;
          }
        }

        if (u + 2 <= length) {
          if (evaluate(ii, length, u, u + 2, VRNA_DECOMP_EXT_STEM_EXT1, &hc_dat_local)) {
            mm3   = S1[u + 1];
            type  = vrna_get_ptype_window(ii + 1, u, ptype);

            en = c[ii + 1][u - ii - 1] + vrna_E_exterior_stem(type, mm5, mm3, P) + f3[u + 2];

            if (sc) {
              if (sc->energy_up)
                en += sc->energy_up[u + 1][1] +
                      sc->energy_up[ii][1];

              if (sc->f)
                en += sc->f(ii, length, u, u + 2, VRNA_DECOMP_EXT_STEM_EXT1, sc->data);
            }

            if (fij == en) {
              vrna_bts_push(bt_stack,
                            ((vrna_sect_t){
                              .i = u + 2,
                              .j = j,
                              .ml = VRNA_MX_FLAG_F3}));

              vrna_bts_push(bt_stack,
                            ((vrna_sect_t){
                              .i = ii + 1,
                              .j = u,
                              .ml = VRNA_MX_FLAG_C}));

              return 1;
            }
          }
        } else {
          if (evaluate(ii, length, ii + 1, u, VRNA_DECOMP_EXT_STEM, &hc_dat_local)) {
            mm3   = (u < length) ? S1[u + 1] : -1;
            type  = vrna_get_ptype_window(ii + 1, u, ptype);

            en = c[ii + 1][u - ii - 1] +
                 vrna_E_exterior_stem(type, mm5, mm3, P);

            if (sc) {
              if (sc->energy_up) {
                en += sc->energy_up[ii][1];
                if (u < length)
                  en += sc->energy_up[u + 1][1];
              }

              if (sc->f)
                en += sc->f(ii, length, ii + 1, u, VRNA_DECOMP_EXT_STEM, sc->data);
            }

            if (fij == en) {
              vrna_bts_push(bt_stack,
                            ((vrna_sect_t){
                              .i = u + 2,
                              .j = j,
                              .ml = VRNA_MX_FLAG_F3}));

              vrna_bts_push(bt_stack,
                            ((vrna_sect_t){
                              .i = ii + 1,
                              .j = u,
                              .ml = VRNA_MX_FLAG_C}));

              return 1;
            }
          }
        }

        if (evaluate(ii, length, u, u + 1, VRNA_DECOMP_EXT_STEM_EXT1, &hc_dat_local)) {
          type = vrna_get_ptype_window(ii + 1, u, ptype);

          en = c[ii + 1][u - ii - 1] +
               vrna_E_exterior_stem(type, mm5, -1, P) +
               f3[u + 1];

          if (sc) {
            if (sc->energy_up)
              en += sc->energy_up[ii][1];

            if (sc->f)
              en += sc->f(ii, length, u, u + 1, VRNA_DECOMP_EXT_STEM_EXT1, sc->data);
          }

          if (fij == en) {
            vrna_bts_push(bt_stack,
                          ((vrna_sect_t){
                            .i = u + 1,
                            .j = j,
                            .ml = VRNA_MX_FLAG_F3}));

            vrna_bts_push(bt_stack,
                          ((vrna_sect_t){
                            .i = ii + 1,
                            .j = u,
                            .ml = VRNA_MX_FLAG_C}));

            return 1;
          }
        }

        if (u + 2 <= length) {
          if (evaluate(ii, length, u, u + 2, VRNA_DECOMP_EXT_STEM_EXT, &hc_dat_local)) {
            mm3   = S1[u + 1];
            type  = vrna_get_ptype_window(ii, u, ptype);

            en = c[ii][u - ii] +
                 vrna_E_exterior_stem(type, -1, mm3, P) +
                 f3[u + 2];

            if (sc) {
              if (sc->energy_up)
                en += sc->energy_up[u + 1][1];

              if (sc->f)
                en += sc->f(ii, length, u, u + 2, VRNA_DECOMP_EXT_STEM_EXT, sc->data);
            }

            if (fij == en) {
              vrna_bts_push(bt_stack,
                            ((vrna_sect_t){
                              .i = u + 2,
                              .j = j,
                              .ml = VRNA_MX_FLAG_F3}));

              vrna_bts_push(bt_stack,
                            ((vrna_sect_t){
                              .i = ii,
                              .j = u,
                              .ml = VRNA_MX_FLAG_C}));

              return 1;
            }
          }
        } else {
          if (evaluate(ii, length, ii, u, VRNA_DECOMP_EXT_STEM, &hc_dat_local)) {
            mm3   = (u < length) ? S1[u + 1] : -1;
            type  = vrna_get_ptype_window(ii, u, ptype);

            en = c[ii][u - ii] +
                 vrna_E_exterior_stem(type, -1, mm3, P);

            if (sc) {
              if ((sc->energy_up) && (u < length))
                en += sc->energy_up[u + 1][1];

              if (sc->f)
                en += sc->f(ii, length, ii, u, VRNA_DECOMP_EXT_STEM, sc->data);
            }

            if (fij == en) {
              vrna_bts_push(bt_stack,
                            ((vrna_sect_t){
                              .i = u + 2,
                              .j = j,
                              .ml = VRNA_MX_FLAG_F3}));

              vrna_bts_push(bt_stack,
                            ((vrna_sect_t){
                              .i = ii,
                              .j = u,
                              .ml = VRNA_MX_FLAG_C}));
              return 1;
            }
          }
        }

        if (evaluate(ii, length, u, u + 1, VRNA_DECOMP_EXT_STEM_EXT, &hc_dat_local)) {
          type = vrna_get_ptype_window(ii, u, ptype);

          en = c[ii][u - ii] +
               vrna_E_exterior_stem(type, -1, -1, P) + f3[u + 1];

          if (sc)
            if (sc->f)
              en += sc->f(ii, length, u, u + 1, VRNA_DECOMP_EXT_STEM_EXT, sc->data);

          if (fij == en) {
            vrna_bts_push(bt_stack,
                          ((vrna_sect_t){
                            .i = u + 1,
                            .j = j,
                            .ml = VRNA_MX_FLAG_F3}));

            vrna_bts_push(bt_stack,
                          ((vrna_sect_t){
                            .i = ii,
                            .j = u,
                            .ml = VRNA_MX_FLAG_C}));

            return 1;
          }
        }
      }

      break;
  }

  return 0;
}


PRIVATE unsigned int
bt_ext_loop_f3_comparative(vrna_fold_compound_t *fc,
                           unsigned int         i,
                           unsigned int         j,
                           vrna_bps_t           bp_stack VRNA_UNUSED,
                           vrna_bts_t           bt_stack)
{
  short                 **S, **S5, **S3;
  unsigned int          n, type, ss, n_seq, **a2s, ii, u, dangle_model, with_gquad;
  int                   fij, cc, fj, *f3, **c, **ggg;
  vrna_param_t          *P;
  vrna_md_t             *md;
  vrna_sc_t             **scs;
  vrna_hc_eval_f        evaluate;
  struct hc_ext_def_dat hc_dat_local;

  n             = fc->length;
  n_seq         = fc->n_seq;
  S             = fc->S;
  S5            = fc->S5;   /* S5[s][i] holds next base 5' of i in sequence s */
  S3            = fc->S3;   /* Sl[s][i] holds next base 3' of i in sequence s */
  a2s           = fc->a2s;
  P             = fc->params;
  md            = &(P->model_details);
  scs           = fc->scs;
  f3            = fc->matrices->f3_local;
  c             = fc->matrices->c_local;
  ggg           = fc->matrices->ggg_local;
  dangle_model  = md->dangles;
  with_gquad    = md->gquad;
  evaluate      = prepare_hc_ext_def_window(fc, &hc_dat_local);

  ii = i;

  /* nibble off unpaired 5' bases */
  do {
    fij = f3[ii];
    fj  = INF;

    if (evaluate(ii, n, ii + 1, n, VRNA_DECOMP_EXT_EXT, &hc_dat_local)) {
      fj = f3[ii + 1];
      if (scs) {
        for (ss = 0; ss < n_seq; ss++)
          if (scs[ss]) {
            if (scs[ss]->energy_up)
              fj += scs[ss]->energy_up[ii][1];

            if (scs[ss]->f)
              fj += scs[ss]->f(ii, n, ii + 1, n, VRNA_DECOMP_EXT_EXT, scs[ss]->data);
          }
      }
    }

    if (++ii > j)
      break;
  } while (fij == fj);
  ii--;

  if (ii >= j) {
    /* no more pairs */
    return 1;
  }

  /*
   *  must have found a decomposition
   *  i is paired. Find pairing partner
   */
  switch (dangle_model) {
    /* no dangles */
    case 0:
      for (u = j; u > ii; u--) {
        if (with_gquad) {
          if (fij == ggg[ii][u - ii] + f3[u + 1]) {
            vrna_bts_push(bt_stack,
                          ((vrna_sect_t){
                            .i = u + 1,
                            .j = j,
                            .ml = VRNA_MX_FLAG_F3}));

            vrna_bts_push(bt_stack,
                          ((vrna_sect_t){
                            .i = ii,
                            .j = u,
                            .ml = VRNA_MX_FLAG_G}));
            return 1;
          }
        }

        if (evaluate(ii, n, u, u + 1, VRNA_DECOMP_EXT_STEM_EXT, &hc_dat_local)) {
          cc = c[ii][u - ii];
          for (ss = 0; ss < n_seq; ss++) {
            type  = vrna_get_ptype_md(S[ss][ii], S[ss][u], md);
            cc    += vrna_E_exterior_stem(type, -1, -1, P);
          }

          if (fij == cc + f3[u + 1]) {
            vrna_bts_push(bt_stack,
                          ((vrna_sect_t){
                            .i = u + 1,
                            .j = j,
                            .ml = VRNA_MX_FLAG_F3}));

            vrna_bts_push(bt_stack,
                          ((vrna_sect_t){
                            .i = ii,
                            .j = u,
                            .ml = VRNA_MX_FLAG_C}));

            return 1;
          }
        }
      }
      break;

    /* dangles on both sides */
    case 2:
      for (u = j; u > ii; u--) {
        if (with_gquad) {
          if (fij == ggg[ii][u - ii] + f3[u + 1]) {
            vrna_bts_push(bt_stack,
                          ((vrna_sect_t){
                            .i = u + 1,
                            .j = j,
                            .ml = VRNA_MX_FLAG_F3}));

            vrna_bts_push(bt_stack,
                          ((vrna_sect_t){
                            .i = ii,
                            .j = u,
                            .ml = VRNA_MX_FLAG_G}));
            return 1;
          }
        }

        if (evaluate(ii, n, u, u + 1, VRNA_DECOMP_EXT_STEM_EXT, &hc_dat_local)) {
          cc = c[ii][u - ii];
          for (ss = 0; ss < n_seq; ss++) {
            type  = vrna_get_ptype_md(S[ss][ii], S[ss][u], md);
            cc    +=
              vrna_E_exterior_stem(type, (a2s[ss][ii] > 1) ? S5[ss][ii] : -1,
                                   (a2s[ss][u] < a2s[ss][n]) ? S3[ss][u] : -1, P);
          }

          if (fij == cc + f3[u + 1]) {
            vrna_bts_push(bt_stack,
                          ((vrna_sect_t){
                            .i = u + 1,
                            .j = j,
                            .ml = VRNA_MX_FLAG_F3}));

            vrna_bts_push(bt_stack,
                          ((vrna_sect_t){
                            .i = ii,
                            .j = u,
                            .ml = VRNA_MX_FLAG_C}));

            return 1;
          }
        }
      }
      break;
  }

  return 0;
}


PRIVATE unsigned int
bt_ext_loop_f3_pp(vrna_fold_compound_t  *fc,
                  unsigned int          *i,
                  unsigned int          maxj)
{
  unsigned int j, start;

  j     = 0;
  start = *i;

  if (fc) {
    char                  **ptype;
    short                 *S1;
    unsigned int          ii, length, type, traced2, dangle_model, with_gquad, maxdist;
    int                   cc, **c, **ggg, *f3, fij;
    vrna_param_t          *P;
    vrna_md_t             *md;
    vrna_sc_t             *sc;
    vrna_hc_eval_f        evaluate;
    struct hc_ext_def_dat hc_dat_local;
#ifdef VRNA_WITH_SVM
    int                   zsc_pre_filter;
    vrna_zsc_dat_t        zsc_data;
#endif

    length        = fc->length;
    S1            = fc->sequence_encoding;
    ptype         = fc->ptype_local;
    f3            = fc->matrices->f3_local;
    c             = fc->matrices->c_local;
    ggg           = fc->matrices->ggg_local;
    sc            = fc->sc;
    P             = fc->params;
    md            = &(P->model_details);
    dangle_model  = md->dangles;
    with_gquad    = md->gquad;
    maxdist       = MIN2((unsigned int)fc->window_size, maxj);
    traced2       = 0;
    ii            = start;
    evaluate      = prepare_hc_ext_def_window(fc, &hc_dat_local);
#ifdef VRNA_WITH_SVM
    zsc_data        = fc->zscore_data;
    zsc_pre_filter  = ((zsc_data) &&
                       (zsc_data->filter_on) &&
                       (zsc_data->pre_filter) &&
                       ((unsigned int)zsc_data->current_i == ii)) ? 1 : 0;
#endif

    fij = f3[start];

    /* try to nibble off unpaired 5' bases */
    if ((sc) && (evaluate(start, length, start + 1, length, VRNA_DECOMP_EXT_EXT, &hc_dat_local))) {
      cc = f3[start + 1];

      if (sc->energy_up)
        cc += sc->energy_up[start][1];

      if (sc->f)
        cc += sc->f(start, length, start + 1, length, VRNA_DECOMP_EXT_EXT, sc->data);

      if (fij == cc)
        /* simple 5' unpaired extensions, so we skip this hit */
        return 0;
    }

    /* get pairing partner j */
    switch (dangle_model) {
      case 0:
        for (j = start + 1; j <= ii + maxdist; j++) {
          if (evaluate(start, length, j, j + 1, VRNA_DECOMP_EXT_STEM_EXT, &hc_dat_local)) {
#ifdef VRNA_WITH_SVM
            if ((zsc_pre_filter) &&
                (zsc_data->current_z[j] > zsc_data->min_z))
              continue;

#endif
            type = vrna_get_ptype_window(start, j, ptype);

            cc = c[start][j - start] +
                 vrna_E_exterior_stem(type, -1, -1, P) +
                 f3[j + 1];

            if ((sc) && (sc->f))
              cc += sc->f(start, length, j, j + 1, VRNA_DECOMP_EXT_STEM_EXT, sc->data);

            if (fij == cc) {
              traced2 = 1;
              break;
            }
          }

          if (with_gquad) {
            cc = ggg[start][j - start] + f3[j + 1];
            if (fij == cc) {
              traced2 = 1;
              break;
            }
          }
        }
        break;

      case 2:
        for (j = start + 1; j <= ii + maxdist; j++) {
          if (evaluate(start, length, j, j + 1, VRNA_DECOMP_EXT_STEM_EXT, &hc_dat_local)) {
#ifdef VRNA_WITH_SVM
            if ((zsc_pre_filter) &&
                (zsc_data->current_z[j] > zsc_data->min_z))
              continue;

#endif
            type = vrna_get_ptype_window(start, j, ptype);

            cc = c[start][j - start] +
                 vrna_E_exterior_stem(type, (start > 1) ? S1[start - 1] : -1,
                                      (j < length) ? S1[j + 1] : -1, P) +
                 f3[j + 1];

            if ((sc) && (sc->f))
              cc += sc->f(start, length, j, j + 1, VRNA_DECOMP_EXT_STEM_EXT, sc->data);

            if (fij == cc) {
              traced2 = 1;
              break;
            }
          }

          if (with_gquad) {
            cc = ggg[start][j - start] +
                 f3[j + 1];

            if (fij == cc) {
              traced2 = 1;
              break;
            }
          }
        }
        break;

      default:
        for (j = start + 1; j <= ii + maxdist; j++) {
          if (evaluate(start, length, j, j + 1, VRNA_DECOMP_EXT_STEM_EXT, &hc_dat_local)) {
            type = vrna_get_ptype_window(start, j, ptype);

            cc = c[start][j - start] +
                 vrna_E_exterior_stem(type, -1, -1, P) +
                 f3[j + 1];

            if ((sc) && (sc->f))
              cc += sc->f(start, length, j, j + 1, VRNA_DECOMP_EXT_STEM_EXT, sc->data);

            if (fij == cc) {
              traced2 = 1;
              break;
            }
          }

          if (j < length) {
            if (evaluate(start, length, j, j + 2, VRNA_DECOMP_EXT_STEM_EXT, &hc_dat_local)) {
              type  = vrna_get_ptype_window(start, j, ptype);
              cc    = c[start][j - start] +
                      vrna_E_exterior_stem(type, -1, S1[j + 1], P) +
                      f3[j + 2];

              if (sc) {
                if (sc->energy_up)
                  cc += sc->energy_up[j + 1][1];

                if (sc->f)
                  cc += sc->f(start, length, j, j + 2, VRNA_DECOMP_EXT_STEM_EXT, sc->data);
              }

              if (fij == cc) {
                traced2 = 1;
                break;
              }
            }
          }

          if (with_gquad) {
            cc = ggg[start][j - start] +
                 f3[j + 1];

            if (fij == cc) {
              traced2 = 1;
              break;
            }
          }

          if (evaluate(start, length, j, j + 1, VRNA_DECOMP_EXT_STEM_EXT1, &hc_dat_local)) {
            type = vrna_get_ptype_window(start + 1, j, ptype);

            cc = c[start + 1][j - (start + 1)] +
                 vrna_E_exterior_stem(type, S1[start], -1, P) +
                 f3[j + 1];

            if (sc) {
              if (sc->energy_up)
                cc += sc->energy_up[start][1];

              if (sc->f)
                cc += sc->f(start, length, j, j + 1, VRNA_DECOMP_EXT_STEM_EXT1, sc->data);
            }

            if (fij == cc) {
              traced2 = 1;
              break;
            }
          }

          if (j < length) {
            if (evaluate(start, length, j, j + 2, VRNA_DECOMP_EXT_STEM_EXT1, &hc_dat_local)) {
              type = vrna_get_ptype_window(start + 1, j, ptype);

              cc = c[start + 1][j - (start + 1)] +
                   vrna_E_exterior_stem(type, S1[start], S1[j + 1], P) +
                   f3[j + 2];

              if (sc) {
                if (sc->energy_up)
                  cc += sc->energy_up[start][1] +
                        sc->energy_up[j + 1][1];

                if (sc->f)
                  cc += sc->f(start, length, j, j + 2, VRNA_DECOMP_EXT_STEM_EXT1, sc->data);
              }

              if (fij == cc) {
                traced2 = 1;
                break;
              }
            }
          }
        }
        break;
    }

    if (!traced2) {
      j = 0;
      vrna_log_error("backtrack failed in short backtrack for position %u", *i);
    }
  }

  *i = start;
  return j;
}


PRIVATE unsigned int
bt_ext_loop_f3_pp_comparative(vrna_fold_compound_t  *fc,
                              unsigned int          *i,
                              unsigned int          maxj)
{
  unsigned int j, start;

  j = 0;

  if (fc) {
    short                 **S, **S5, **S3;
    unsigned int          tt, s, n_seq, **a2s, traced2, length, dangle_model, with_gquad, maxdist;
    int                   cc, **c, **ggg, *f3, fij;
    vrna_param_t          *P;
    vrna_md_t             *md;
    vrna_sc_t             **scs;
    vrna_hc_eval_f        evaluate;
    struct hc_ext_def_dat hc_dat_local;

    length        = fc->length;
    n_seq         = fc->n_seq;
    S             = fc->S;
    S5            = fc->S5;   /* S5[s][start] holds next base 5' of start in sequence s */
    S3            = fc->S3;   /* Sl[s][start] holds next base 3' of start in sequence s */
    a2s           = fc->a2s;
    f3            = fc->matrices->f3_local;
    c             = fc->matrices->c_local;
    ggg           = fc->matrices->ggg_local;
    scs           = fc->scs;
    P             = fc->params;
    md            = &(P->model_details);
    dangle_model  = md->dangles;
    with_gquad    = md->gquad;
    maxdist       = MIN2((unsigned int)fc->window_size, maxj);
    traced2       = 0;
    start         = *i;
    evaluate      = prepare_hc_ext_def_window(fc, &hc_dat_local);

    fij = f3[start];

    /* try to nibble off unpaired 5' bases */
    if ((scs) && (evaluate(start, length, start + 1, length, VRNA_DECOMP_EXT_EXT, &hc_dat_local))) {
      cc = f3[start + 1];

      for (s = 0; s < n_seq; s++)
        if (scs[s]) {
          if (scs[s]->energy_up)
            cc += scs[s]->energy_up[start][1];

          if (scs[s]->f)
            cc += scs[s]->f(start,
                            length,
                            start + 1,
                            length,
                            VRNA_DECOMP_EXT_EXT,
                            scs[s]->data);
        }

      if (fij == cc)
        /* simple 5' unpaired extensions, so we skip this hit */
        return 0;
    }

    /* get pairing partner j */
    switch (dangle_model) {
      case 0:
        for (j = start + 1; j <= MIN2(start + maxdist, length - 1); j++) {
          if (evaluate(start, length, j, j + 1, VRNA_DECOMP_EXT_STEM_EXT, &hc_dat_local)) {
            cc = c[start][j - start] +
                 f3[j + 1];

            for (s = 0; s < n_seq; s++) {
              tt  = vrna_get_ptype_md(S[s][start], S[s][j], md);
              cc  += vrna_E_exterior_stem(tt, -1, -1, P);
            }

            if (scs) {
              for (s = 0; s < n_seq; s++)
                if ((scs[s]) && (scs[s]->f))
                  cc += scs[s]->f(start, length, j, j + 1, VRNA_DECOMP_EXT_STEM_EXT, scs[s]->data);
            }

            if (fij == cc) {
              traced2 = 1;
              break;
            }
          }

          if (with_gquad) {
            cc = ggg[start][j - start] +
                 f3[j + 1];

            if (fij == cc) {
              traced2 = 1;
              break;
            }
          }
        }

        if ((!traced2) && (length <= start + maxdist)) {
          j = length;
          if (evaluate(start, length, start, j, VRNA_DECOMP_EXT_STEM, &hc_dat_local)) {
            cc = c[start][j - start];

            for (s = 0; s < n_seq; s++) {
              tt  = vrna_get_ptype_md(S[s][start], S[s][j], md);
              cc  += vrna_E_exterior_stem(tt, -1, -1, P);
            }

            if (scs) {
              for (s = 0; s < n_seq; s++)
                if ((scs[s]) && (scs[s]->f))
                  cc += scs[s]->f(start, length, start, j, VRNA_DECOMP_EXT_STEM, scs[s]->data);
            }

            if (fij == cc) {
              traced2 = 1;
              break;
            }
          }

          if (with_gquad) {
            cc = ggg[start][j - start];

            if (fij == cc) {
              traced2 = 1;
              break;
            }
          }
        }

        break;

      case 2:
        for (j = start + 1; j <= MIN2(start + maxdist, length - 1); j++) {
          if (evaluate(start, length, j, j + 1, VRNA_DECOMP_EXT_STEM_EXT, &hc_dat_local)) {
            cc = c[start][j - start] +
                 f3[j + 1];

            for (s = 0; s < n_seq; s++) {
              tt  = vrna_get_ptype_md(S[s][start], S[s][j], md);
              cc  +=
                vrna_E_exterior_stem(tt,
                                     (a2s[s][start] > 1) ? S5[s][start] : -1,
                                     (a2s[s][j] < a2s[s][length]) ? S3[s][j] : -1,
                                     P);
            }

            if (scs) {
              for (s = 0; s < n_seq; s++)
                if ((scs[s]) && (scs[s]->f))
                  cc += scs[s]->f(start, length, j, j + 1, VRNA_DECOMP_EXT_STEM_EXT, scs[s]->data);
            }

            if (fij == cc) {
              traced2 = 1;
              break;
            }
          }

          if (with_gquad) {
            cc = ggg[start][j - start] +
                 f3[j + 1];

            if (fij == cc) {
              traced2 = 1;
              break;
            }
          }
        }

        if ((!traced2) && (length <= start + maxdist)) {
          j = length;
          if (evaluate(start, length, start, j, VRNA_DECOMP_EXT_STEM, &hc_dat_local)) {
            cc = c[start][j - start];
            for (s = 0; s < n_seq; s++) {
              tt  = vrna_get_ptype_md(S[s][start], S[s][j], md);
              cc  += vrna_E_exterior_stem(tt, (a2s[s][start] > 1) ? S5[s][start] :  -1, -1, P);
            }

            if (scs) {
              for (s = 0; s < n_seq; s++)
                if ((scs[s]) && (scs[s]->f))
                  cc += scs[s]->f(start, length, start, j, VRNA_DECOMP_EXT_STEM, scs[s]->data);
            }

            if (fij == cc) {
              traced2 = 1;
              break;
            }
          }

          if (with_gquad) {
            cc = ggg[start][j - start];

            if (fij == cc) {
              traced2 = 1;
              break;
            }
          }
        }

        break;
    }

    if (!traced2) {
      j = 0;
      vrna_log_error("backtrack failed in short backtrack for position %u", *i);
    }
  }

  return j;
}


#ifndef VRNA_DISABLE_BACKWARD_COMPATIBILITY

PUBLIC int
vrna_BT_ext_loop_f3_pp(vrna_fold_compound_t *fc,
                       int                  *i,
                       int                  maxj)
{
  int           r;
  unsigned int  ii = (unsigned int)(*i);

  r = (int)vrna_bt_exterior_f3_pp(fc, &ii, (unsigned int)maxj);

  *i = (int)ii;

  return r;
}


PUBLIC int
vrna_BT_ext_loop_f3(vrna_fold_compound_t  *fc,
                    int                   *k,
                    int                   maxdist,
                    int                   *i,
                    int                   *j,
                    vrna_bp_stack_t       *bp_stack,
                    unsigned int          *stack_count)
{
  int r = 0;

  if ((fc) &&
      (k) &&
      (i) &&
      (j) &&
      (bp_stack) &&
      (stack_count)) {
    vrna_bps_t  bps = vrna_bps_init(0);
    vrna_bts_t  bts = vrna_bts_init(0);

    switch (fc->type) {
      case VRNA_FC_TYPE_SINGLE:
        r = (int)bt_ext_loop_f3(fc, (unsigned int)(*k), (unsigned int)maxdist, bps, bts);
        break;

      case VRNA_FC_TYPE_COMPARATIVE:
        r =
          (int)bt_ext_loop_f3_comparative(fc, (unsigned int)(*k), (unsigned int)maxdist, bps, bts);
        break;
    }

    while (vrna_bts_size(bts) > 0) {
      vrna_sect_t s = vrna_bts_pop(bts);
      if (s.ml == VRNA_MX_FLAG_F3) {
        *k = (int)s.i;
      } else if (s.ml == VRNA_MX_FLAG_C) {
        bp_stack[++(*stack_count)].i  = s.i;
        bp_stack[*stack_count].j      = s.j;
        *i                            = (int)s.i;
        *j                            = (int)s.j;
      } else if (s.ml == VRNA_MX_FLAG_G) {
        r = vrna_bt_gquad_mfe(fc, s.i, s.j, bps);
      }
    }

    while (vrna_bps_size(bps) > 0) {
      vrna_bp_t bp = vrna_bps_pop(bps);
      bp_stack[++(*stack_count)].i  = bp.i;
      bp_stack[*stack_count].j      = bp.j;
    }

    vrna_bps_free(bps);
    vrna_bts_free(bts);
  }

  return r;
}


#endif
