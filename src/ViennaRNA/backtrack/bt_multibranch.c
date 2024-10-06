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
#include "ViennaRNA/structured_domains.h"
#include "ViennaRNA/unstructured_domains.h"
#include "ViennaRNA/eval/multibranch.h"
#include "ViennaRNA/backtrack/gquad.h"
#include "ViennaRNA/backtrack/multibranch.h"

#include "ViennaRNA/intern/grammar_dat.h"

#ifdef __GNUC__
# define INLINE inline
#else
# define INLINE
#endif

#include "ViennaRNA/constraints/multibranch_hc.inc"
#include "ViennaRNA/constraints/multibranch_sc.inc"

/*
 #################################
 # PRIVATE FUNCTION DECLARATIONS #
 #################################
 */

PRIVATE unsigned int
bt_mb_loop_split(vrna_fold_compound_t *fc,
                 unsigned int         i,
                 unsigned int         j,
                 vrna_bps_t           bp_stack,
                 vrna_bts_t           bt_stack);


PRIVATE unsigned int
bt_mb_loop(vrna_fold_compound_t *fc,
           unsigned int         i,
           unsigned int         j,
           int                  en,
           vrna_bps_t           bp_stack,
           vrna_bts_t           bt_stack);


/*
 #################################
 # BEGIN OF FUNCTION DEFINITIONS #
 #################################
 */
PUBLIC unsigned int
vrna_bt_m(vrna_fold_compound_t  *fc,
          unsigned int          i,
          unsigned int          j,
          vrna_bps_t            bp_stack,
          vrna_bts_t            bt_stack)
{
  unsigned int  ret = 0;
  int           e, *idx;

  if ((fc) &&
      (bt_stack) &&
      (bp_stack) &&
      (fc->matrices)) {
    idx = fc->jindx;
    e   = (fc->matrices->type == VRNA_MX_WINDOW) ?
          fc->matrices->fML_local[i][j - i] :
          fc->matrices->fML[idx[j] + i];

    if ((ret = bt_mb_loop_split(fc, i, j, bp_stack, bt_stack))) {
      ret = 1;
    } else if (fc->aux_grammar) {
      for (size_t c = 0; c < vrna_array_size(fc->aux_grammar->m); c++) {
        if ((fc->aux_grammar->m[c].cb_bt) &&
            (ret =
               fc->aux_grammar->m[c].cb_bt(fc, i, j, e, bp_stack, bt_stack,
                                           fc->aux_grammar->m[c].data)))
          break;
      }
    }
  }

  return ret;
}


PUBLIC unsigned int
vrna_bt_multibranch_split(vrna_fold_compound_t  *fc,
                          unsigned int          i,
                          unsigned int          j,
                          vrna_bps_t            bp_stack,
                          vrna_bts_t            bt_stack)
{
  if ((fc) &&
      (bp_stack) &&
      (bt_stack))
    return bt_mb_loop_split(fc, i, j, bp_stack, bt_stack);

  return 0;
}


PUBLIC unsigned int
vrna_bt_multibranch_loop(vrna_fold_compound_t  *fc,
                         unsigned int          i,
                         unsigned int          j,
                         int                   en,
                         vrna_bps_t            bp_stack,
                         vrna_bts_t            bt_stack)
{
  if ((fc) &&
      (bp_stack) &&
      (bt_stack))
    return bt_mb_loop(fc, i, j, en, bp_stack, bt_stack);

  return 0;
}


PRIVATE unsigned int
bt_mb_loop_split(vrna_fold_compound_t *fc,
                 unsigned int         i,
                 unsigned int         j,
                 vrna_bps_t           bp_stack VRNA_UNUSED,
                 vrna_bts_t           bt_stack)
{
  unsigned char         sliding_window;
  char                  *ptype, **ptype_local;
  short                 *S1, **SS, **S5, **S3;
  unsigned int          n_seq, s, with_gquad, dangle_model, u, kk, cnt,
                        with_ud, type, type_2;
  int                   ij, fij, fi, en, *my_c, *my_fML, *idx, *rtype,
                        en2, **c_local, **fML_local, **ggg_local;
  vrna_param_t          *P;
  vrna_md_t             *md;
  vrna_ud_t             *domains_up;
  vrna_hc_eval_f        evaluate;
  struct hc_mb_def_dat  hc_dat_local;
  struct sc_mb_dat      sc_wrapper;
  vrna_smx_csr(int)     *c_gq;

  sliding_window  = (fc->hc->type == VRNA_HC_WINDOW) ? 1 : 0;
  n_seq           = (fc->type == VRNA_FC_TYPE_SINGLE) ? 1 : fc->n_seq;
  P               = fc->params;
  md              = &(P->model_details);
  idx             = (sliding_window) ? NULL : fc->jindx;
  ptype           = (fc->type == VRNA_FC_TYPE_SINGLE) ? (sliding_window ? NULL : fc->ptype) : NULL;
  ptype_local     = (sliding_window) ? fc->ptype_local : NULL;
  rtype           = &(md->rtype[0]);
  S1              = (fc->type == VRNA_FC_TYPE_SINGLE) ? fc->sequence_encoding : NULL;
  SS              = (fc->type == VRNA_FC_TYPE_SINGLE) ? NULL : fc->S;
  S5              = (fc->type == VRNA_FC_TYPE_SINGLE) ? NULL : fc->S5;
  S3              = (fc->type == VRNA_FC_TYPE_SINGLE) ? NULL : fc->S3;
  domains_up      = fc->domains_up;

  my_c      = (sliding_window) ? NULL : fc->matrices->c;
  my_fML    = (sliding_window) ? NULL : fc->matrices->fML;
  c_gq      = (sliding_window) ? NULL : fc->matrices->c_gq;
  c_local   = (sliding_window) ? fc->matrices->c_local : NULL;
  fML_local = (sliding_window) ? fc->matrices->fML_local : NULL;
  ggg_local = (sliding_window) ? fc->matrices->ggg_local : NULL;

  with_gquad    = md->gquad;
  with_ud       = (domains_up && domains_up->energy_cb) ? 1 : 0;
  dangle_model  = md->dangles;
  evaluate      = prepare_hc_mb_def(fc, &hc_dat_local);

  init_sc_mb(fc, &sc_wrapper);

  if (with_ud) {
    /* nibble off unpaired stretches at 3' site */
    do {
      fij = (sliding_window) ? fML_local[i][j - i] : my_fML[idx[j] + i];
      fi  = INF;

      /* process regular unpaired nucleotides (unbound by ligand) first */
      if (evaluate(i, j, i, j - 1, VRNA_DECOMP_ML_ML, &hc_dat_local)) {
        fi = P->MLbase *
             n_seq;
        fi += (sliding_window) ? fML_local[i][j - 1 - i] : my_fML[idx[j - 1] + i];

        if (sc_wrapper.red_ml)
          fi += sc_wrapper.red_ml(i, j, i, j - 1, &sc_wrapper);

        if (j == i)
          return 0; /* no more pairs */

        if (fij == fi) {
          j--;

          if (j < i)
            return 0; /* no more pairs */

          continue;
        }
      }

      /* next try to nibble off ligand */
      for (cnt = 0; cnt < (unsigned int)domains_up->uniq_motif_count; cnt++) {
        u   = domains_up->uniq_motif_size[cnt];
        kk  = j - u + 1;
        if ((kk >= i) && evaluate(i, j, i, j - u, VRNA_DECOMP_ML_ML, &hc_dat_local)) {
          en = domains_up->energy_cb(fc,
                                     kk,
                                     j,
                                     VRNA_UNSTRUCTURED_DOMAIN_MB_LOOP | VRNA_UNSTRUCTURED_DOMAIN_MOTIF,
                                     domains_up->data);

          if (sc_wrapper.red_ml)
            en += sc_wrapper.red_ml(i, j, i, j - u, &sc_wrapper);

          fi = en +
               u * P->MLbase *
               n_seq;
          fi += (sliding_window) ? fML_local[i][kk - 1 - i] : my_fML[idx[kk - 1] + i];

          if (fij == fi) {
            /* skip remaining motifs after first hit */
            j = kk - 1;
            break;
          }
        }
      }

      if (j < i)
        return 0; /* no more pairs */
    } while (fij == fi);

    /* nibble off unpaired stretches at 5' site */
    do {
      fij = (sliding_window) ? fML_local[i][j - i] : my_fML[idx[j] + i];
      fi  = INF;

      /* again, process regular unpaired nucleotides (unbound by ligand) first */
      if (evaluate(i, j, i + 1, j, VRNA_DECOMP_ML_ML, &hc_dat_local)) {
        fi = P->MLbase *
             n_seq;
        fi += (sliding_window) ? fML_local[i + 1][j - (i + 1)] : my_fML[idx[j] + i + 1];

        if (sc_wrapper.red_ml)
          fi += sc_wrapper.red_ml(i, j, i + 1, j, &sc_wrapper);

        if (i + 1 >= j)
          return 0; /* no more pairs */

        if (fij == fi) {
          i++;
          continue;
        }
      }

      /* next try to nibble off ligand again */
      for (cnt = 0; cnt < (unsigned int)domains_up->uniq_motif_count; cnt++) {
        u   = domains_up->uniq_motif_size[cnt];
        kk  = i + u - 1;
        if ((kk <= j) && evaluate(i, j, i + u, j, VRNA_DECOMP_ML_ML, &hc_dat_local)) {
          en = domains_up->energy_cb(fc,
                                     i,
                                     kk,
                                     VRNA_UNSTRUCTURED_DOMAIN_MB_LOOP | VRNA_UNSTRUCTURED_DOMAIN_MOTIF,
                                     domains_up->data);

          if (sc_wrapper.red_ml)
            en += sc_wrapper.red_ml(i, j, i + u, j, &sc_wrapper);

          fi = en +
               u * P->MLbase *
               n_seq;

          fi += (sliding_window) ? fML_local[kk + 1][j - (kk + 1)] : my_fML[idx[j] + kk + 1];

          if (fij == fi) {
            /* skip remaining motifs after first hit */
            i = kk + 1;
            break;
          }
        }
      }

      if (i > j)
        return 0; /* no more pairs */
    } while (fij == fi);
  } else {
    /* nibble off unpaired 3' bases */
    do {
      fij = (sliding_window) ? fML_local[i][j - i] : my_fML[idx[j] + i];
      fi  = INF;

      if (evaluate(i, j, i, j - 1, VRNA_DECOMP_ML_ML, &hc_dat_local)) {
        fi = P->MLbase *
             n_seq;
        fi += (sliding_window) ? fML_local[i][j - 1 - i] : my_fML[idx[j - 1] + i];

        if (sc_wrapper.red_ml)
          fi += sc_wrapper.red_ml(i, j, i, j - 1, &sc_wrapper);
      }

      if (--j == 0)
        break;
    } while (fij == fi);
    j++;

    /* nibble off unpaired 5' bases */
    do {
      fij = (sliding_window) ? fML_local[i][j - i] : my_fML[idx[j] + i];
      fi  = INF;

      if (evaluate(i, j, i + 1, j, VRNA_DECOMP_ML_ML, &hc_dat_local)) {
        fi = P->MLbase *
             n_seq;
        fi += (sliding_window) ? fML_local[i + 1][j - (i + 1)] : my_fML[idx[j] + i + 1];

        if (sc_wrapper.red_ml)
          fi += sc_wrapper.red_ml(i, j, i + 1, j, &sc_wrapper);
      }

      if (++i >= j)
        break;
    } while (fij == fi);
    i--;

    if (j <= i) /* no more pairs */
      return 0;
  }

  ij = (sliding_window) ? 0 : idx[j] + i;

  /* 1. test for single component */

  if (with_gquad) {
#ifndef VRNA_DISABLE_C11_FEATURES
    en = (sliding_window) ? ggg_local[i][j - i] : vrna_smx_csr_get(c_gq, i, j, INF);
#else
    en = (sliding_window) ? ggg_local[i][j - i] : vrna_smx_csr_int_get(c_gq, i, j, INF);
#endif

    if (en != INF) {
      en += vrna_E_multibranch_stem(0, -1, -1, P) *
            n_seq;

      if (fij == en) {
        vrna_bts_push(bt_stack,
                      ((vrna_sect_t){
                        .i = i,
                        .j = j,
                        .ml = VRNA_MX_FLAG_G}));

        return 1;
      }
    }
  }

  en = (sliding_window) ? c_local[i][j - i] : my_c[ij];

  if (sc_wrapper.red_stem)
    en += sc_wrapper.red_stem(i, j, i, j, &sc_wrapper);

  switch (dangle_model) {
    case 0:
      if (evaluate(i, j, i, j, VRNA_DECOMP_ML_STEM, &hc_dat_local)) {
        en2 = 0;
        switch (fc->type) {
          case VRNA_FC_TYPE_SINGLE:
            type = (sliding_window) ? vrna_get_ptype_window(i, j, ptype_local) : vrna_get_ptype(
              ij,
              ptype);
            en2 = vrna_E_multibranch_stem(type, -1, -1, P);
            break;

          case VRNA_FC_TYPE_COMPARATIVE:
            for (s = 0; s < n_seq; s++) {
              type  = vrna_get_ptype_md(SS[s][i], SS[s][j], md);
              en2   += vrna_E_multibranch_stem(type, -1, -1, P);
            }
            break;
        }

        if (fij == en + en2) {
          vrna_bts_push(bt_stack,
                        ((vrna_sect_t){
                          .i = i,
                          .j = j,
                          .ml = VRNA_MX_FLAG_C}));

          return 1;
        }
      }

      break;

    case 2:
      if (evaluate(i, j, i, j, VRNA_DECOMP_ML_STEM, &hc_dat_local)) {
        en2 = 0;
        switch (fc->type) {
          case VRNA_FC_TYPE_SINGLE:
            type = (sliding_window) ? vrna_get_ptype_window(i, j, ptype_local) : vrna_get_ptype(
              ij,
              ptype);
            en2 = vrna_E_multibranch_stem(type, S1[i - 1], S1[j + 1], P);
            break;

          case VRNA_FC_TYPE_COMPARATIVE:
            for (s = 0; s < n_seq; s++) {
              type  = vrna_get_ptype_md(SS[s][i], SS[s][j], md);
              en2   += vrna_E_multibranch_stem(type, S5[s][i], S3[s][j], P);
            }
            break;
        }

        if (fij == en + en2) {
          vrna_bts_push(bt_stack,
                        ((vrna_sect_t){
                          .i = i,
                          .j = j,
                          .ml = VRNA_MX_FLAG_C}));

          return 1;
        }
      }

      break;

    default:
      if (evaluate(i, j, i, j, VRNA_DECOMP_ML_STEM, &hc_dat_local)) {
        en2 = 0;
        switch (fc->type) {
          case VRNA_FC_TYPE_SINGLE:
            type = (sliding_window) ? vrna_get_ptype_window(i, j, ptype_local) : vrna_get_ptype(
              ij,
              ptype);
            en2 = vrna_E_multibranch_stem(type, -1, -1, P);
            break;

          case VRNA_FC_TYPE_COMPARATIVE:
            for (s = 0; s < n_seq; s++) {
              type  = vrna_get_ptype_md(SS[s][i], SS[s][j], md);
              en2   += vrna_E_multibranch_stem(type, -1, -1, P);
            }
            break;
        }

        if (fij == en + en2) {
          vrna_bts_push(bt_stack,
                        ((vrna_sect_t){
                          .i = i,
                          .j = j,
                          .ml = VRNA_MX_FLAG_C}));

          return 1;
        }
      }

      if (evaluate(i, j, i + 1, j, VRNA_DECOMP_ML_STEM, &hc_dat_local)) {
        en2 = P->MLbase *
              n_seq;
        en2 += (sliding_window) ? c_local[i + 1][j - (i + 1)] : my_c[ij + 1];

        switch (fc->type) {
          case VRNA_FC_TYPE_SINGLE:
            type =
              (sliding_window) ? vrna_get_ptype_window(i + 1, j, ptype_local) : vrna_get_ptype(
                ij + 1,
                ptype);
            en2 += vrna_E_multibranch_stem(type, S1[i], -1, P);
            break;

          case VRNA_FC_TYPE_COMPARATIVE:
            for (s = 0; s < n_seq; s++) {
              type  = vrna_get_ptype_md(SS[s][i], SS[s][j], md);
              en2   += vrna_E_multibranch_stem(type, S5[s][i], -1, P);
            }
            break;
        }

        if (sc_wrapper.red_stem)
          en2 += sc_wrapper.red_stem(i, j, i + 1, j, &sc_wrapper);

        if (fij == en2) {
          vrna_bts_push(bt_stack,
                        ((vrna_sect_t){
                          .i = i + 1,
                          .j = j,
                          .ml = VRNA_MX_FLAG_C}));

          return 1;
        }
      }

      if (evaluate(i, j, i, j - 1, VRNA_DECOMP_ML_STEM, &hc_dat_local)) {
        en2 = P->MLbase *
              n_seq;
        en2 += (sliding_window) ? c_local[i][j - 1 - i] : my_c[idx[j - 1] + i];

        switch (fc->type) {
          case VRNA_FC_TYPE_SINGLE:
            type =
              (sliding_window) ? vrna_get_ptype_window(i, j - 1,
                                                       ptype_local) : vrna_get_ptype(
                idx[j - 1] + i,
                ptype);
            en2 += vrna_E_multibranch_stem(type, -1, S1[j], P);
            break;

          case VRNA_FC_TYPE_COMPARATIVE:
            for (s = 0; s < n_seq; s++) {
              type  = vrna_get_ptype_md(SS[s][i], SS[s][j], md);
              en2   += vrna_E_multibranch_stem(type, -1, S3[s][j], P);
            }
            break;
        }

        if (sc_wrapper.red_stem)
          en2 += sc_wrapper.red_stem(i, j, i, j - 1, &sc_wrapper);

        if (fij == en2) {
          vrna_bts_push(bt_stack,
                        ((vrna_sect_t){
                          .i = i,
                          .j = j - 1,
                          .ml = VRNA_MX_FLAG_C}));

          return 1;
        }
      }

      if (evaluate(i, j, i + 1, j - 1, VRNA_DECOMP_ML_STEM, &hc_dat_local)) {
        en2 = 2 * P->MLbase *
              n_seq;
        en2 += (sliding_window) ? c_local[i + 1][j - 1 - (i + 1)] : my_c[idx[j - 1] + i + 1];

        switch (fc->type) {
          case VRNA_FC_TYPE_SINGLE:
            type =
              (sliding_window) ? vrna_get_ptype_window(i + 1, j - 1,
                                                       ptype_local) : vrna_get_ptype(
                idx[j - 1] + i + 1,
                ptype);
            en2 += vrna_E_multibranch_stem(type, S1[i], S1[j], P);
            break;

          case VRNA_FC_TYPE_COMPARATIVE:
            for (s = 0; s < n_seq; s++) {
              type  = vrna_get_ptype_md(SS[s][i], SS[s][j], md);
              en2   += vrna_E_multibranch_stem(type,
                                               S5[s][i],
                                               S3[s][j],
                                               P);
            }
            break;
        }

        if (sc_wrapper.red_stem)
          en2 += sc_wrapper.red_stem(i, j, i + 1, j - 1, &sc_wrapper);

        if (fij == en2) {
          vrna_bts_push(bt_stack,
                        ((vrna_sect_t){
                          .i = i + 1,
                          .j = j - 1,
                          .ml = VRNA_MX_FLAG_C}));

          return 1;
        }
      }

      break;
  }

  /* 2. Test for possible split point */
  for (u = i + 1; u + 1 < j; u++) {
    if (evaluate(i, j, u, u + 1, VRNA_DECOMP_ML_ML_ML, &hc_dat_local)) {
      if (sliding_window)
        en = fML_local[i][u - i] +
             fML_local[u + 1][j - (u + 1)];
      else
        en = my_fML[idx[u] + i] +
             my_fML[idx[j] + u + 1];

      if (sc_wrapper.decomp_ml)
        en += sc_wrapper.decomp_ml(i, j, u, u + 1, &sc_wrapper);

      if (fij == en) {
        vrna_bts_push(bt_stack,
                      ((vrna_sect_t){
                        .i = i,
                        .j = u,
                        .ml = VRNA_MX_FLAG_M}));
        vrna_bts_push(bt_stack,
                      ((vrna_sect_t){
                        .i = u + 1,
                        .j = j,
                        .ml = VRNA_MX_FLAG_M}));

        return 1;
      }
    }
  }

  /* 3. last chance! Maybe coax stack */
  if (dangle_model == 3) {
    int ik, k1j;
    k1j = (sliding_window) ? 0 : idx[j] + i + 2;
    for (u = i + 1; u <= j - 2; u++, k1j++) {
      ik = (sliding_window) ? 0 : idx[u] + i;
      if (evaluate(i, u, u + 1, j, VRNA_DECOMP_ML_COAXIAL_ENC, &hc_dat_local)) {
        en = 2 * P->MLintern[1] *
             n_seq;
        en += (sliding_window) ? c_local[i][u - i] + c_local[u + 1][j - (u + 1)] : my_c[ik] +
              my_c[k1j];

        switch (fc->type) {
          case VRNA_FC_TYPE_SINGLE:
            type =
              (sliding_window) ? rtype[vrna_get_ptype_window(i, u,
                                                             ptype_local)] : rtype[vrna_get_ptype(
                                                                                     ik,
                                                                                     ptype)
              ];
            type_2 =
              (sliding_window) ? rtype[vrna_get_ptype_window(u + 1, j,
                                                             ptype_local)] : rtype[vrna_get_ptype(
                                                                                     k1j, ptype)];
            en += P->stack[type][type_2];
            break;

          case VRNA_FC_TYPE_COMPARATIVE:
            for (s = 0; s < n_seq; s++) {
              type    = vrna_get_ptype_md(SS[s][u], SS[s][i], md);
              type_2  = vrna_get_ptype_md(SS[s][j], SS[s][u + 1], md);
              en      += P->stack[type][type_2];
            }
            break;
        }

        if (sc_wrapper.coaxial_enc)
          en += sc_wrapper.coaxial_enc(i, u, u + 1, j, &sc_wrapper);

        if (fij == en) {
          vrna_bts_push(bt_stack,
                        ((vrna_sect_t){
                          .i = i,
                          .j = u,
                          .ml = VRNA_MX_FLAG_C}));
          vrna_bts_push(bt_stack,
                        ((vrna_sect_t){
                          .i = u + 1,
                          .j = j,
                          .ml = VRNA_MX_FLAG_C}));

          return 1;
        }
      }
    }
  }

  return 0;
}


PRIVATE unsigned int
bt_mb_loop(vrna_fold_compound_t *fc,
           unsigned int         i,
           unsigned int         j,
           int                  en,
           vrna_bps_t           bp_stack VRNA_UNUSED,
           vrna_bts_t           bt_stack)
{
  unsigned char         sliding_window;
  char                  *ptype, **ptype_local;
  short                 s5, s3, *S1, **SS, **S5, **S3;
  unsigned int          *sn, n_seq, s, *tt, c1, c2, p, q, r, dangle_model, type, type_2, turn;
  int                   ij, e, tmp_en, *idx, *my_c, *my_fML, *rtype, **c_local, **fML_local;
  vrna_param_t          *P;
  vrna_md_t             *md;
  vrna_hc_eval_f        evaluate;
  struct hc_mb_def_dat  hc_dat_local;
  struct sc_mb_dat      sc_wrapper;

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
  evaluate      = prepare_hc_mb_def(fc, &hc_dat_local);

  init_sc_mb(fc, &sc_wrapper);

  p = i + 1;
  q = j - 1;

  r = q - 1;

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

  if (evaluate(i, j, p, q, VRNA_DECOMP_PAIR_ML, &hc_dat_local)) {
    e = en -
        P->MLclosing *
        n_seq;

    if (dangles == 2) {
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

    for (r = i + 2; r < j - 2; ++r) {
      if (evaluate(p, q, r, r + 1, VRNA_DECOMP_ML_ML_ML, &hc_dat_local)) {
        if (sliding_window)
          tmp_en = fML_local[p][r - p] +
                   fML_local[r + 1][q - (r + 1)];
        else
          tmp_en = my_fML[idx[r] + p] +
                   my_fML[idx[q] + r + 1];

        if (sc_wrapper.decomp_ml)
          tmp_en += sc_wrapper.decomp_ml(p, q, r, r + 1, &sc_wrapper);

        if (e == tmp_en)
          goto odd_dangles_exit;
      }
    }
  }

  if (dangles % 2) {
    /* odd dangles need more special treatment */
    if (evaluate(i, j, p + 1, q, VRNA_DECOMP_PAIR_ML, &hc_dat_local)) {
      e = en -
          (P->MLclosing + P->MLbase) *
          n_seq;

      if (sc_wrapper.pair5)
        e -= sc_wrapper.pair5(i, j, &sc_wrapper);

      switch (fc->type) {
        case VRNA_FC_TYPE_SINGLE:
          e -= vrna_E_multibranch_stem(type, -1, s3, P);
          break;

        case VRNA_FC_TYPE_COMPARATIVE:
          for (s = 0; s < n_seq; s++)
            e -= vrna_E_multibranch_stem(tt[s], -1, S3[s][i], P);
          break;
      }

      for (r = p + 1; r < q - 1; ++r) {
        if (evaluate(p + 1, q, r, r + 1, VRNA_DECOMP_ML_ML_ML, &hc_dat_local)) {
          if (sliding_window)
            tmp_en = fML_local[p + 1][r - (p + 1)] +
                     fML_local[r + 1][q - (r + 1)];
          else
            tmp_en = my_fML[idx[r] + p + 1] +
                     my_fML[idx[q] + r + 1];

          if (sc_wrapper.decomp_ml)
            tmp_en += sc_wrapper.decomp_ml(p + 1, q, r, r + 1, &sc_wrapper);

          if (e == tmp_en) {
            p += 1;
            goto odd_dangles_exit;
          }
        }
      }
    }

    if (evaluate(i, j, p, q - 1, VRNA_DECOMP_PAIR_ML, &hc_dat_local)) {
      e = en -
          (P->MLclosing + P->MLbase) *
          n_seq;

      if (sc_wrapper.pair3)
        e -= sc_wrapper.pair3(i, j, &sc_wrapper);

      switch (fc->type) {
        case VRNA_FC_TYPE_SINGLE:
          e -= vrna_E_multibranch_stem(type, s5, -1, P);
          break;

        case VRNA_FC_TYPE_COMPARATIVE:
          for (s = 0; s < n_seq; s++)
            e -= vrna_E_multibranch_stem(tt[s], S5[s][j], -1, P);
          break;
      }

      for (r = p + 1; r < q - 1; ++r) {
        if (evaluate(p, q - 1, r, r + 1, VRNA_DECOMP_ML_ML_ML, &hc_dat_local)) {
          if (sliding_window)
            tmp_en = fML_local[p][r - p] +
                     fML_local[r + 1][q - 1 - (r + 1)];
          else
            tmp_en = my_fML[idx[r] + p] +
                     my_fML[idx[q - 1] + r + 1];

          if (sc_wrapper.decomp_ml)
            tmp_en += sc_wrapper.decomp_ml(p, q - 1, r, r + 1, &sc_wrapper);

          ;

          if (e == tmp_en) {
            q -= 1;
            goto odd_dangles_exit;
          }
        }
      }
    }

    if (evaluate(i, j, p + 1, q - 1, VRNA_DECOMP_PAIR_ML, &hc_dat_local)) {
      e = en -
          (P->MLclosing + 2 * P->MLbase) *
          n_seq;

      if (sc_wrapper.pair53)
        e -= sc_wrapper.pair53(i, j, &sc_wrapper);

      switch (fc->type) {
        case VRNA_FC_TYPE_SINGLE:
          e -= vrna_E_multibranch_stem(type, s5, s3, P);
          break;

        case VRNA_FC_TYPE_COMPARATIVE:
          for (s = 0; s < n_seq; s++)
            e -= vrna_E_multibranch_stem(tt[s], S5[s][j], S3[s][i], P);
          break;
      }

      for (r = p + 1; r < q - 1; ++r) {
        if (evaluate(p + 1, q - 1, r, r + 1, VRNA_DECOMP_ML_ML_ML, &hc_dat_local)) {
          if (sliding_window)
            tmp_en = fML_local[p + 1][r - (p + 1)] +
                     fML_local[r + 1][q - 1 - (r + 1)];
          else
            tmp_en = my_fML[idx[r] + p + 1] +
                     my_fML[idx[q - 1] + r + 1];

          if (sc_wrapper.decomp_ml)
            tmp_en += sc_wrapper.decomp_ml(p + 1, q - 1, r, r + 1, &sc_wrapper);

          if (e == tmp_en) {
            p += 1;
            q -= 1;
            goto odd_dangles_exit;
          }
        }
      }
    }

    /*
     * coaxial stacking of (i.j) with (i+1.r) or (r.j-1)
     * use MLintern[1] since coax stacked pairs don't get TerminalAU
     */
    if (dangle_model == 3) {
      e = en -
          (P->MLclosing + 2 * P->MLintern[1]) *
          n_seq;

      if (sc_wrapper.pair)
        e -= sc_wrapper.pair(i, j, &sc_wrapper);

      if (fc->type == VRNA_FC_TYPE_SINGLE)
        type = rtype[type];

      for (r = p + 1; r < q - 1; ++r) {
        if (evaluate(i, j, p, r, VRNA_DECOMP_ML_COAXIAL, &hc_dat_local)) {
          if (sliding_window)
            tmp_en = c_local[p][r - p] +
                     fML_local[r + 1][q - (r + 1)];
          else
            tmp_en = my_c[idx[r] + p] +
                     my_fML[idx[q] + r + 1];

          switch (fc->type) {
            case VRNA_FC_TYPE_SINGLE:
              if (sliding_window)
                type_2 = rtype[vrna_get_ptype_window(p, r, ptype_local)];
              else
                type_2 = rtype[vrna_get_ptype(idx[r] + p, ptype)];

              tmp_en += P->stack[type][type_2];
              break;

            case VRNA_FC_TYPE_COMPARATIVE:
              for (s = 0; s < n_seq; s++) {
                type_2  = vrna_get_ptype_md(SS[s][r], SS[s][p], md);
                tmp_en  += P->stack[tt[s]][type_2];
              }
              break;
          }

          if (sc_wrapper.coaxial_cls)
            tmp_en += sc_wrapper.coaxial_cls(i, j, p, r, &sc_wrapper);

          if (e == tmp_en) {
            c1 = VRNA_MX_FLAG_C;
            goto odd_dangles_exit;
          }
        }

        if (evaluate(i, j, r + 1, q, VRNA_DECOMP_ML_COAXIAL, &hc_dat_local)) {
          if (sliding_window)
            tmp_en = c_local[r + 1][q - (r + 1)] +
                     fML_local[p][r - p];
          else
            tmp_en = my_c[idx[q] + r + 1] +
                     my_fML[idx[r] + p];

          switch (fc->type) {
            case VRNA_FC_TYPE_SINGLE:
              if (sliding_window)
                type_2 = rtype[vrna_get_ptype_window(r + 1, q, ptype_local)];
              else
                type_2 = rtype[vrna_get_ptype(idx[q] + r + 1, ptype)];

              tmp_en += P->stack[type][type_2];
              break;

            case VRNA_FC_TYPE_COMPARATIVE:
              for (s = 0; s < n_seq; s++) {
                type_2  = vrna_get_ptype_md(SS[s][q], SS[s][r + 1], md);
                tmp_en  += P->stack[tt[s]][type_2];
              }
              break;
          }

          if (sc_wrapper.coaxial_cls)
            tmp_en += sc_wrapper.coaxial_cls(i, j, r + 1, q, &sc_wrapper);

          if (e == tmp_en) {
            c2 = VRNA_MX_FLAG_C;
            goto odd_dangles_exit;
          }
        }
      }
    }

odd_dangles_exit:

    free(tt);
  }

  if (r + 2 < j) {
    vrna_bts_push(bt_stack,
                  ((vrna_sect_t){
                    .i = p,
                    .j = r,
                    .ml = c1}));
    vrna_bts_push(bt_stack,
                  ((vrna_sect_t){
                    .i = r + 1,
                    .j = q,
                    .ml = c2}));

    return 1;
  } else {
#if 0
    /* Y shaped ML loops fon't work yet */
    if (dangle_model == 3) {
      d5  = P->dangle5[tt][S1[j - 1]];
      d3  = P->dangle3[tt][S1[i + 1]];
      /* (i,j) must close a Y shaped ML loop with coax stacking */
      if (cij == fML[indx[j - 2] + i + 2] + mm + d3 + d5 + P->MLbase + P->MLbase) {
        i1  = i + 2;
        j1  = j - 2;
      } else if (cij == fML[indx[j - 2] + i + 1] + mm + d5 + P->MLbase) {
        j1 = j - 2;
      } else if (cij == fML[indx[j - 1] + i + 2] + mm + d3 + P->MLbase) {
        i1 = i + 2;
      } else /* last chance */
      if (cij != fML[indx[j - 1] + i + 1] + mm + P->MLbase) {
        fprintf(stderr, "backtracking failed in repeat");
      }

      /* if we arrive here we can express cij via fML[i1,j1]+dangles */
      bt_stack[++s].i = i1;
      bt_stack[s].j   = j1;
    }

#endif
  }

  return 0;
}


#ifndef VRNA_DISABLE_BACKWARD_COMPATIBILITY

PUBLIC int
vrna_BT_mb_loop_split(vrna_fold_compound_t  *fc,
                      int                   *i,
                      int                   *j,
                      int                   *k,
                      int                   *l,
                      unsigned int          *c1,
                      unsigned int          *c2,
                      vrna_bp_stack_t       *bp_stack,
                      unsigned int          *stack_count)
{
  int r = 0;

  if ((fc) &&
      (i) &&
      (j) &&
      (k) &&
      (l) &&
      (c1) &&
      (c2) &&
      (bp_stack) &&
      (stack_count)) {
    vrna_sect_t s;
    vrna_bps_t  bps = vrna_bps_init(0);
    vrna_bts_t  bts = vrna_bts_init(0);

    r = (int)bt_mb_loop_split(fc, *i, *j, bps, bts);

    *i = *j = *k = *l = -1;

    if (vrna_bts_size(bts) > 0) {
      s = vrna_bts_pop(bts);
      if (s.ml == VRNA_MX_FLAG_G) {
        r = vrna_bt_gquad_mfe(fc, s.i, s.j, bps);
      } else {
        *i = s.i;
        *j = s.j;
        *c1 = s.ml;
      }
    }

    if (vrna_bts_size(bts) > 0) {
      s = vrna_bts_pop(bts);
      if (s.ml == VRNA_MX_FLAG_G) {
        r = vrna_bt_gquad_mfe(fc, s.i, s.j, bps);
      } else {
        *k  = s.i;
        *l  = s.j;
        *c2 = s.ml;
      }
    }

    while (vrna_array_size(bps) > 0) {
      vrna_bp_t bp = vrna_bps_pop(bps);
      bp_stack[++(*stack_count)].i  = bp.i;
      bp_stack[*stack_count].j      = bp.j;
    }

    vrna_bps_free(bps);
    vrna_bts_free(bts);
  }

  return r;
}


PUBLIC int
vrna_BT_mb_loop(vrna_fold_compound_t  *fc,
                int                   *i,
                int                   *j,
                int                   *k,
                int                   en,
                unsigned int          *c1,
                unsigned int          *c2)
{
  int r = 0;

  if ((fc) &&
      (i) &&
      (j) &&
      (k) &&
      (c1) &&
      (c2)) {
    vrna_sect_t s;
    vrna_bps_t  bps = vrna_bps_init(0);
    vrna_bts_t  bts = vrna_bts_init(0);

    r = (int)bt_mb_loop(fc, *i, *j, en, bps, bts);

    *i = *j = *k = -1;

    if (vrna_bts_size(bts) > 0) {
      s = vrna_bts_pop(bts);
      if (s.ml == VRNA_MX_FLAG_G) {
        r = vrna_bt_gquad_mfe(fc, s.i, s.j, bps);
      } else {
        *i  = s.i;
        *k  = s.j;
        *c1 = s.ml;
      }
    }

    if (vrna_bts_size(bts) > 0) {
      s = vrna_bts_pop(bts);
      if (s.ml == VRNA_MX_FLAG_G) {
        r = vrna_bt_gquad_mfe(fc, s.i, s.j, bps);
      } else {
        *k  = s.i;
        *j  = s.j;
        *c2 = s.ml;
      }
    }

    vrna_bps_free(bps);
    vrna_bts_free(bts);
  }

  return r;
}


#endif
