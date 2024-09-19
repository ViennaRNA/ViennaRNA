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
#include "ViennaRNA/alphabet.h"
#include "ViennaRNA/constraints/hard.h"
#include "ViennaRNA/constraints/soft.h"
#include "ViennaRNA/loops/gquad.h"
#include "ViennaRNA/structured_domains.h"
#include "ViennaRNA/unstructured_domains.h"
#include "ViennaRNA/eval/exterior.h"
#include "ViennaRNA/backtrack/gquad.h"
#include "ViennaRNA/backtrack/exterior.h"

#include "ViennaRNA/grammar.inc"

#ifdef __GNUC__
# define INLINE inline
#else
# define INLINE
#endif

#ifdef VRNA_WITH_SVM
#include "ViennaRNA/zscore_dat.inc"
#endif

#include "ViennaRNA/constraints/exterior_hc.inc"
#include "ViennaRNA/constraints/exterior_sc.inc"

/*
 #################################
 # PRIVATE FUNCTION DECLARATIONS #
 #################################
 */

PRIVATE unsigned int
bt_ext_loop_f5(vrna_fold_compound_t *fc,
               unsigned int         j,
               vrna_bps_t           bp_stack,
               vrna_bts_t           bt_stack);


PRIVATE unsigned int
bt_ext_loop_f5_comparative(vrna_fold_compound_t *fc,
                           unsigned int         j,
                           vrna_bps_t           bp_stack,
                           vrna_bts_t           bt_stack);



PRIVATE unsigned int
bt_ext_loop_f3(vrna_fold_compound_t *fc,
               unsigned int                  *k,
               unsigned int                  maxdist,
               unsigned int                  *i,
               unsigned int                  *j,
               vrna_bps_t           bp_stack);


PRIVATE unsigned int
bt_ext_loop_f3_comparative(vrna_fold_compound_t *fc,
                           unsigned int                  *k,
                           unsigned int                  maxdist,
                           unsigned int                  *i,
                           unsigned int                  *j,
                           vrna_bps_t           bp_stack);


PRIVATE int
bt_ext_loop_f3_pp(vrna_fold_compound_t  *fc,
                  unsigned int                   *i,
                  unsigned int                   maxj);


PRIVATE int
bt_ext_loop_f3_pp_comparative(vrna_fold_compound_t  *fc,
                              unsigned int                   *i,
                              unsigned int                   maxj);


/*
 #################################
 # BEGIN OF FUNCTION DEFINITIONS #
 #################################
 */

PUBLIC unsigned int
vrna_bt_f(vrna_fold_compound_t *fc,
          unsigned int         i,
          unsigned int         j,
          vrna_bps_t           bp_stack,
          vrna_bts_t           bt_stack)
{
  unsigned int  n, ret = 0;
  int           e;

  e   = INF;

  if ((fc) &&
      (bp_stack) &&
      (bt_stack)) {
    n = fc->length;

    if ((i == 1) &&
        (fc->matrices) &&
        (fc->matrices->f5)) {
      e   = fc->matrices->f5[j];
      ret = vrna_bt_ext_loop_f5(fc, j, bp_stack, bt_stack);
    }

    if ((!ret) &&
        (fc->aux_grammar)) {
      for (size_t c = 0; c < vrna_array_size(fc->aux_grammar->f); c++) {
        if ((fc->aux_grammar->f[c].cb_bt) &&
            ((ret = fc->aux_grammar->f[c].cb_bt(fc, i, j, e, bp_stack, bt_stack, fc->aux_grammar->f[c].data)) != 0))
          break;
      }
    }
  }

  return ret;
}


PUBLIC unsigned int
vrna_bt_ext_loop_f5(vrna_fold_compound_t  *fc,
                    unsigned int          j,
                    vrna_bps_t            bp_stack,
                    vrna_bts_t            bt_stack)
{
  if ((fc) &&
      (bp_stack) &&
      (bt_stack)) {
    switch (fc->type) {
      case VRNA_FC_TYPE_SINGLE:
        return bt_ext_loop_f5(fc, j, bp_stack, bt_stack);
        break;

      case VRNA_FC_TYPE_COMPARATIVE:
        return bt_ext_loop_f5_comparative(fc, j, bp_stack, bt_stack);
        break;
    }
  }

  return 0;
}


PUBLIC unsigned int
vrna_bt_ext_loop_f3(vrna_fold_compound_t  *fc,
                    unsigned int                   *k,
                    unsigned int                   maxdist,
                    unsigned int                   *i,
                    unsigned int                   *j,
                    vrna_bps_t            bp_stack)
{
  if (fc) {
    switch (fc->type) {
      case VRNA_FC_TYPE_SINGLE:
        return bt_ext_loop_f3(fc, k, maxdist, i, j, bp_stack);
        break;

      case VRNA_FC_TYPE_COMPARATIVE:
        return bt_ext_loop_f3_comparative(fc, k, maxdist, i, j, bp_stack);
        break;
    }
  }

  return 0;
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
    unsigned int kk, ii, jj;
    kk = *k;
    ii = *i;
    jj = *j;

    r = (int)vrna_bt_ext_loop_f3(fc, &kk, (unsigned int)maxdist, &ii, &jj, bps);

    while (vrna_bps_size(bps) > 0) {
      vrna_bp_t bp = vrna_bps_pop(bps);
      bp_stack[++(*stack_count)].i = bp.i;
      bp_stack[*stack_count].j = bp.j;
    }

    vrna_bps_free(bps);

    *k = (int)kk;
    *i = (int)ii;
    *j = (int)jj;
  }

  return r;
}


PUBLIC int
vrna_bt_ext_loop_f3_pp(vrna_fold_compound_t *fc,
                       unsigned int                  *i,
                       unsigned int                  maxj)
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

  return -1;
}


PRIVATE unsigned int
bt_ext_loop_f5(vrna_fold_compound_t *fc,
               unsigned int         j,
               vrna_bps_t           bp_stack,
               vrna_bts_t           bt_stack)
{
  char                      *ptype;
  short                     mm5, mm3, *S1;
  unsigned int              length, *sn, type, u, dangle_model, with_gquad, cnt,
                            ii, with_ud;
  int                       fij, fi, en, e, *my_f5, *my_c, *idx, e_gq;
  vrna_param_t              *P;
  vrna_md_t                 *md;
  vrna_sc_t                 *sc;
  vrna_ud_t                 *domains_up;
  vrna_hc_eval_f  evaluate;
  struct hc_ext_def_dat     hc_dat_local;
  vrna_smx_csr(int)         *c_gq;

  length        = fc->length;
  P             = fc->params;
  md            = &(P->model_details);
  sn            = fc->strand_number;
  sc            = fc->sc;
  my_f5         = fc->matrices->f5;
  my_c          = fc->matrices->c;
  c_gq          = fc->matrices->c_gq;
  domains_up    = fc->domains_up;
  idx           = fc->jindx;
  ptype         = fc->ptype;
  S1            = fc->sequence_encoding;
  dangle_model  = md->dangles;
  with_gquad    = md->gquad;
  with_ud       = (domains_up && domains_up->energy_cb) ? 1 : 0;
  evaluate      = prepare_hc_ext_def(fc, &hc_dat_local);

  /* nibble off unpaired 3' stretches harboring bound ligands (interspersed with unpaired nucleotides) */
  if (with_ud) {
    do {
      fij = my_f5[j];
      fi  = INF;

      /* try nibble off one unpaired nucleotide first */
      if (evaluate(1, j, 1, j - 1, VRNA_DECOMP_EXT_EXT, &hc_dat_local)) {
        fi = my_f5[j - 1];

        if (sc) {
          if (sc->energy_up)
            fi += sc->energy_up[j][1];

          if (sc->f)
            fi += sc->f(1, j, 1, j - 1, VRNA_DECOMP_EXT_EXT, sc->data);
        }

        if (j == 1) {
          /* no more pairs */
          return 1;
        }

        if (fij == fi) {
          j--;
          continue;
        }
      }

      /* next, try nibble off a ligand */
      for (cnt = 0; cnt < domains_up->uniq_motif_count; cnt++) {
        u   = domains_up->uniq_motif_size[cnt];
        if ((j >= u) &&
            evaluate(1, j, 1, j - u, VRNA_DECOMP_EXT_EXT, &hc_dat_local)) {
          ii  = j - u + 1;
          en = domains_up->energy_cb(fc,
                                     ii,
                                     j,
                                     VRNA_UNSTRUCTURED_DOMAIN_EXT_LOOP | VRNA_UNSTRUCTURED_DOMAIN_MOTIF,
                                     domains_up->data);
          if (sc) {
            if (sc->energy_up)
              en += sc->energy_up[ii][u];

            if (sc->f)
              en += sc->f(1, j, 1, j - u, VRNA_DECOMP_EXT_EXT, sc->data);
          }

          fi  = my_f5[ii - 1];
          fi  += en;

          if (fij == fi) {
            /* skip remaining motifs after first hit */
            j = ii - 1;
            break;
          }
        }
      }

      if (j == 0) {
        /* no more pairs */
        return 1;
      }
    } while (fij == fi);
  } else {
    /* nibble off unpaired 3' bases */
    do {
      fij = my_f5[j];
      fi  = INF;

      if (evaluate(1, j, 1, j - 1, VRNA_DECOMP_EXT_EXT, &hc_dat_local)) {
        fi = my_f5[j - 1];

        if (sc) {
          if (sc->energy_up)
            fi += sc->energy_up[j][1];

          if (sc->f)
            fi += sc->f(1, j, 1, j - 1, VRNA_DECOMP_EXT_EXT, sc->data);
        }
      }

      if (--j == 0)
        break;
    } while (fij == fi);
    j++;
  }

  if (j < 2) {
    /* no more pairs */
    return 1;
  }

  /* must have found a decomposition */
  switch (dangle_model) {
    case 0:   /* j is paired. Find pairing partner */
      for (u = j - 1; u >= 1; u--) {
        if (with_gquad) {
#ifndef VRNA_DISABLE_C11_FEATURES
          e_gq =  vrna_smx_csr_get(c_gq, u, j, INF);
#else
          e_gq =  vrna_smx_csr_int_get(c_gq, u, j, INF);
#endif
          if ((e_gq != INF) &&
              (fij == my_f5[u - 1] + e_gq)) {
            vrna_bts_push(bt_stack,
                          ((vrna_sect_t){
                            .i = 1,
                            .j = u - 1,
                            .ml = VRNA_MX_FLAG_F5}));

            vrna_bts_push(bt_stack,
                          ((vrna_sect_t){
                            .i = u,
                            .j = j,
                            .ml = VRNA_MX_FLAG_G}));

            return 1;
          }
        }

        if (evaluate(1, j, u - 1, u, VRNA_DECOMP_EXT_EXT_STEM, &hc_dat_local)) {
          type = vrna_get_ptype(idx[j] + u, ptype);

          en = my_c[idx[j] + u];
          if (sc)
            if (sc->f)
              en += sc->f(1, j, u - 1, u, VRNA_DECOMP_EXT_EXT_STEM, sc->data);

#if 0
          if (sn[j] != sn[u])
            en += P->DuplexInit;
#endif

          if (fij == vrna_E_exterior_stem(type, -1, -1, P) + en + my_f5[u - 1]) {
            vrna_bts_push(bt_stack,
                          ((vrna_sect_t){
                            .i = 1,
                            .j = u - 1,
                            .ml = VRNA_MX_FLAG_F5}));

            vrna_bts_push(bt_stack,
                          ((vrna_sect_t){
                            .i = u,
                            .j = j,
                            .ml = VRNA_MX_FLAG_C}));

            return 1;
          }
        }
      }
      break;

    case 2:
      mm3 = ((j < length) && (sn[j + 1] == sn[j])) ? S1[j + 1] : -1;
      for (u = j - 1; u >= 1; u--) {
        if (with_gquad) {
#ifndef VRNA_DISABLE_C11_FEATURES
          e_gq =  vrna_smx_csr_get(c_gq, u, j, INF);
#else
          e_gq =  vrna_smx_csr_int_get(c_gq, u, j, INF);
#endif
          if ((e_gq != INF) &&
              (fij == my_f5[u - 1] + e_gq)) {
            vrna_bts_push(bt_stack,
                          ((vrna_sect_t){
                            .i = 1,
                            .j = u - 1,
                            .ml = VRNA_MX_FLAG_F5}));

            vrna_bts_push(bt_stack,
                          ((vrna_sect_t){
                            .i = u,
                            .j = j,
                            .ml = VRNA_MX_FLAG_G}));

            return 1;
          }
        }

        if (evaluate(1, j, u - 1, u, VRNA_DECOMP_EXT_EXT_STEM, &hc_dat_local)) {
          mm5   = ((u > 1) && (sn[u] == sn[u - 1])) ? S1[u - 1] : -1;
          type  = vrna_get_ptype(idx[j] + u, ptype);

          en = my_c[idx[j] + u];
          if (sc)
            if (sc->f)
              en += sc->f(1, j, u - 1, u, VRNA_DECOMP_EXT_EXT_STEM, sc->data);

#if 0
          if (sn[j] != sn[u])
            en += P->DuplexInit;
#endif

          if (fij == vrna_E_exterior_stem(type, mm5, mm3, P) + en + my_f5[u - 1]) {
            vrna_bts_push(bt_stack,
                          ((vrna_sect_t){
                            .i = 1,
                            .j = u - 1,
                            .ml = VRNA_MX_FLAG_F5}));

            vrna_bts_push(bt_stack,
                          ((vrna_sect_t){
                            .i = u,
                            .j = j,
                            .ml = VRNA_MX_FLAG_C}));

            return 1;
          }
        }
      }
      break;

    default:
      if (with_gquad) {
#ifndef VRNA_DISABLE_C11_FEATURES
        e_gq =  vrna_smx_csr_get(c_gq, 1, j, INF);
#else
        e_gq =  vrna_smx_csr_int_get(c_gq, 1, j, INF);
#endif
        if ((e_gq != INF) &&
            (fij == e_gq)) {

          vrna_bts_push(bt_stack,
                        ((vrna_sect_t){
                          .i = 1,
                          .j = j,
                          .ml = VRNA_MX_FLAG_G}));

          return 1;
        }
      }

      if (evaluate(1, j, 1, j, VRNA_DECOMP_EXT_STEM, &hc_dat_local)) {
        type = vrna_get_ptype(idx[j] + 1, ptype);

        en = my_c[idx[j] + 1];
        if (sc)
          if (sc->f)
            en += sc->f(1, j, 1, j, VRNA_DECOMP_EXT_STEM, sc->data);

        if (fij == en + vrna_E_exterior_stem(type, -1, -1, P)) {
          vrna_bts_push(bt_stack,
                        ((vrna_sect_t){
                          .i = 1,
                          .j = j,
                          .ml = VRNA_MX_FLAG_C}));

          return 1;
        }
      }

      if (evaluate(1, j, 1, j - 1, VRNA_DECOMP_EXT_STEM, &hc_dat_local)) {
        if (sn[j] == sn[j - 1]) {
          mm3   = S1[j];
          type  = vrna_get_ptype(idx[j - 1] + 1, ptype);

          en = my_c[idx[j - 1] + 1];
          if (sc) {
            if (sc->energy_up)
              en += sc->energy_up[j][1];

            if (sc->f)
              en += sc->f(1, j, 1, j - 1, VRNA_DECOMP_EXT_STEM, sc->data);
          }

          if (fij == en + vrna_E_exterior_stem(type, -1, mm3, P)) {
            vrna_bts_push(bt_stack,
                          ((vrna_sect_t){
                            .i = 1,
                            .j = j - 1,
                            .ml = VRNA_MX_FLAG_C}));

            return 1;
          }
        }
      }

      for (u = j - 1; u > 1; u--) {
        if (with_gquad) {
#ifndef VRNA_DISABLE_C11_FEATURES
          e_gq =  vrna_smx_csr_get(c_gq, u, j, INF);
#else
          e_gq =  vrna_smx_csr_int_get(c_gq, u, j, INF);
#endif
          if ((e_gq != INF) &&
              (fij == my_f5[u - 1] + e_gq)) {
            vrna_bts_push(bt_stack,
                          ((vrna_sect_t){
                            .i = 1,
                            .j = u - 1,
                            .ml = VRNA_MX_FLAG_F5}));

            vrna_bts_push(bt_stack,
                          ((vrna_sect_t){
                            .i = u,
                            .j = j,
                            .ml = VRNA_MX_FLAG_G}));

            return 1;
          }
        }

        type = vrna_get_ptype(idx[j] + u, ptype);

        en = my_c[idx[j] + u];
#if 0
        if (sn[j] != sn[u])
          en += P->DuplexInit;
#endif

        if (evaluate(1, j, u - 1, u, VRNA_DECOMP_EXT_EXT_STEM, &hc_dat_local)) {
          e = my_f5[u - 1] +
              en +
              vrna_E_exterior_stem(type, -1, -1, P);

          if (sc)
            if (sc->f)
              e += sc->f(1, j, u - 1, u, VRNA_DECOMP_EXT_EXT_STEM, sc->data);

          if (fij == e) {
            vrna_bts_push(bt_stack,
                          ((vrna_sect_t){
                            .i = 1,
                            .j = u - 1,
                            .ml = VRNA_MX_FLAG_F5}));

            vrna_bts_push(bt_stack,
                          ((vrna_sect_t){
                            .i = u,
                            .j = j,
                            .ml = VRNA_MX_FLAG_C}));

            return 1;
          }
        }

        if (evaluate(1, j, u - 2, u, VRNA_DECOMP_EXT_EXT_STEM, &hc_dat_local)) {
          if (sn[u] == sn[u - 1]) {
            mm5 = S1[u - 1];
            e   = my_f5[u - 2] +
                  en +
                  vrna_E_exterior_stem(type, mm5, -1, P);

            if (sc) {
              if (sc->energy_up)
                e += sc->energy_up[u - 1][1];

              if (sc->f)
                e += sc->f(1, j, u - 2, u, VRNA_DECOMP_EXT_EXT_STEM, sc->data);
            }

            if (fij == e) {
              vrna_bts_push(bt_stack,
                            ((vrna_sect_t){
                              .i = 1,
                              .j = u - 2,
                              .ml = VRNA_MX_FLAG_F5}));

              vrna_bts_push(bt_stack,
                            ((vrna_sect_t){
                              .i = u,
                              .j = j,
                              .ml = VRNA_MX_FLAG_C}));

              return 1;
            }
          }
        }

        type = vrna_get_ptype(idx[j - 1] + u, ptype);

        en = my_c[idx[j - 1] + u];

#if 0
        if (sn[j - 1] != sn[u])
          en += P->DuplexInit;
#endif

        mm5 = (sn[u] == sn[u - 1]) ? S1[u - 1] : -1;
        mm3 = (sn[j] == sn[j - 1]) ? S1[j] : -1;

        if (evaluate(1, j, u - 1, u, VRNA_DECOMP_EXT_EXT_STEM1, &hc_dat_local)) {
          e = my_f5[u - 1] +
              en +
              vrna_E_exterior_stem(type, -1, mm3, P);

          if (sc) {
            if (sc->energy_up)
              e += sc->energy_up[j][1];

            if (sc->f)
              e += sc->f(1, j, u - 1, u, VRNA_DECOMP_EXT_EXT_STEM1, sc->data);
          }

          if (fij == e) {
            vrna_bts_push(bt_stack,
                          ((vrna_sect_t){
                            .i = 1,
                            .j = u - 1,
                            .ml = VRNA_MX_FLAG_F5}));

            vrna_bts_push(bt_stack,
                          ((vrna_sect_t){
                            .i = u,
                            .j = j - 1,
                            .ml = VRNA_MX_FLAG_C}));

            return 1;
          }
        }

        if (evaluate(1, j, u - 2, u, VRNA_DECOMP_EXT_EXT_STEM1, &hc_dat_local)) {
          e = my_f5[u - 2] + en + vrna_E_exterior_stem(type, mm5, mm3, P);
          if (sc) {
            if (sc->energy_up)
              e += sc->energy_up[j][1] +
                   sc->energy_up[u - 1][1];

            if (sc->f)
              e += sc->f(1, j, u - 2, u, VRNA_DECOMP_EXT_EXT_STEM1, sc->data);
          }

          if (fij == e) {
            vrna_bts_push(bt_stack,
                          ((vrna_sect_t){
                            .i = 1,
                            .j = u - 2,
                            .ml = VRNA_MX_FLAG_F5}));

            vrna_bts_push(bt_stack,
                          ((vrna_sect_t){
                            .i = u,
                            .j = j - 1,
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
bt_ext_loop_f5_comparative(vrna_fold_compound_t *fc,
                           unsigned int         j,
                           vrna_bps_t           bp_stack,
                           vrna_bts_t           bt_stack)
{
  unsigned int              **a2s, n;
  short                     **S, **S5, **S3;
  unsigned int              tt, u, dangle_model, with_gquad, n_seq, ss;
  int                       fij, fi, en, e_gq, *my_f5, *my_c, *idx, mm5, mm3;
  vrna_param_t              *P;
  vrna_md_t                 *md;
  vrna_sc_t                 **scs;
  vrna_hc_eval_f evaluate;
  struct hc_ext_def_dat     hc_dat_local;
  vrna_smx_csr(int)         *c_gq;

  n_seq         = fc->n_seq;
  n             = fc->length;
  S             = fc->S;
  S5            = fc->S5;
  S3            = fc->S3;
  a2s           = fc->a2s;
  P             = fc->params;
  md            = &(P->model_details);
  scs           = fc->scs;
  my_f5         = fc->matrices->f5;
  my_c          = fc->matrices->c;
  c_gq          = fc->matrices->c_gq;
  idx           = fc->jindx;
  dangle_model  = md->dangles;
  with_gquad    = md->gquad;
  evaluate      = prepare_hc_ext_def(fc, &hc_dat_local);

  /* nibble off unpaired 3' bases */
  do {
    fij = my_f5[j];
    fi  = INF;

    if (evaluate(1, j, 1, j - 1, VRNA_DECOMP_EXT_EXT, &hc_dat_local)) {
      fi = my_f5[j - 1];

      if (scs) {
        for (ss = 0; ss < n_seq; ss++)
          if (scs[ss]) {
            if (scs[ss]->energy_up)
              fi += scs[ss]->energy_up[a2s[ss][j]][1];

            if (scs[ss]->f)
              fi += scs[ss]->f(1, j, 1, j - 1, VRNA_DECOMP_EXT_EXT, scs[ss]->data);
          }
      }
    }

    if (--j == 0)
      break;
  } while (fij == fi);
  j++;

  if (j < 2) {
    return 1;
  }

  /* must have found a decomposition */
  switch (dangle_model) {
    case 0:   /* j is paired. Find pairing partner */
      for (u = j - 1; u >= 1; u--) {
        if (with_gquad) {
#ifndef VRNA_DISABLE_C11_FEATURES
          e_gq = vrna_smx_csr_get(c_gq, u, j, INF);
#else
          e_gq = vrna_smx_csr_int_get(c_gq, u, j, INF);
#endif
          if ((e_gq != INF) &&
              (fij == my_f5[u - 1] + e_gq)) {
            vrna_bts_push(bt_stack,
                          ((vrna_sect_t){
                            .i = 1,
                            .j = u - 1,
                            .ml = VRNA_MX_FLAG_F5}));

            vrna_bts_push(bt_stack,
                          ((vrna_sect_t){
                            .i = u,
                            .j = j,
                            .ml = VRNA_MX_FLAG_G}));

            return 1;
          }
        }

        if (evaluate(1, j, u - 1, u, VRNA_DECOMP_EXT_EXT_STEM, &hc_dat_local)) {
          en = my_c[idx[j] + u] +
               my_f5[u - 1];

          for (ss = 0; ss < n_seq; ss++) {
            tt  = vrna_get_ptype_md(S[ss][u], S[ss][j], md);
            en  += vrna_E_exterior_stem(tt, -1, -1, P);
          }

          if (scs) {
            for (ss = 0; ss < n_seq; ss++)
              if (scs[ss] && scs[ss]->f)
                en += scs[ss]->f(1, j, u - 1, u, VRNA_DECOMP_EXT_EXT_STEM, scs[ss]->data);
          }

          if (fij == en) {
            vrna_bts_push(bt_stack,
                          ((vrna_sect_t){
                            .i = 1,
                            .j = u - 1,
                            .ml = VRNA_MX_FLAG_F5}));

            vrna_bts_push(bt_stack,
                          ((vrna_sect_t){
                            .i = u,
                            .j = j,
                            .ml = VRNA_MX_FLAG_C}));

            return 1;
          }
        }
      }
      break;

    case 2:
      for (u = j - 1; u >= 1; u--) {
        if (with_gquad) {
#ifndef VRNA_DISABLE_C11_FEATURES
          e_gq = vrna_smx_csr_get(c_gq, u, j, INF);
#else
          e_gq = vrna_smx_csr_int_get(c_gq, u, j, INF);
#endif
          if ((e_gq != INF) &&
              (fij == my_f5[u - 1] + e_gq)) {
            vrna_bts_push(bt_stack,
                          ((vrna_sect_t){
                            .i = 1,
                            .j = u - 1,
                            .ml = VRNA_MX_FLAG_F5}));

            vrna_bts_push(bt_stack,
                          ((vrna_sect_t){
                            .i = u,
                            .j = j,
                            .ml = VRNA_MX_FLAG_G}));

            return 1;
          }
        }

        if (evaluate(1, j, u - 1, u, VRNA_DECOMP_EXT_EXT_STEM, &hc_dat_local)) {
          en = my_c[idx[j] + u] +
               my_f5[u - 1];

          for (ss = 0; ss < n_seq; ss++) {
            tt  = vrna_get_ptype_md(S[ss][u], S[ss][j], md);
            mm5 = (a2s[ss][u] > 1) ? S5[ss][u] : -1;
            mm3 = (a2s[ss][j] < a2s[ss][n]) ? S3[ss][j] : -1;
            en  += vrna_E_exterior_stem(tt, mm5, mm3, P);
          }

          if (scs) {
            for (ss = 0; ss < n_seq; ss++)
              if (scs[ss] && scs[ss]->f)
                en += scs[ss]->f(1, j, u - 1, u, VRNA_DECOMP_EXT_EXT_STEM, scs[ss]->data);
          }

          if (fij == en) {
            vrna_bts_push(bt_stack,
                          ((vrna_sect_t){
                            .i = 1,
                            .j = u - 1,
                            .ml = VRNA_MX_FLAG_F5}));

            vrna_bts_push(bt_stack,
                          ((vrna_sect_t){
                            .i = u,
                            .j = j,
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
bt_ext_loop_f3(vrna_fold_compound_t *fc,
               unsigned int         *k,
               unsigned int         maxdist,
               unsigned int         *i,
               unsigned int         *j,
               vrna_bps_t           bp_stack)
{
  char                      **ptype;
  short                     mm5, mm3, *S1;
  unsigned int              type, length, ii, u, dangle_model, with_gquad;
  int                       fij, fj, *f3, **c, **ggg, en;
  vrna_param_t              *P;
  vrna_md_t                 *md;
  vrna_sc_t                 *sc;
  vrna_hc_eval_f evaluate;
  struct hc_ext_def_dat     hc_dat_local;

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

  ii = *k;

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

    if (++ii > maxdist)
      break;
  } while (fij == fj);
  ii--;

  if (ii >= maxdist) {
    /* no more pairs */
    *i  = *j = 0;
    *k  = 0;
    return 1;
  }

  /*
   *  must have found a decomposition
   *  i is paired. Find pairing partner
   */
  switch (dangle_model) {
    /* no dangles */
    case 0:
      for (u = maxdist; u > ii; u--) {
        if (with_gquad) {
          if (fij == ggg[ii][u - ii] + f3[u + 1]) {
            *i  = *j = 0;
            *k  = u + 1;
            return vrna_bt_gquad_mfe(fc, ii, u, bp_stack);
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
            *i                            = ii;
            *j                            = u;
            *k                            = u + 1;
            vrna_bps_push(bp_stack,
                         (vrna_bp_t){
                          .i = ii,
                          .j = u
                         });
            return 1;
          }
        }
      }
      break;

    /* dangles on both sides */
    case 2:
      mm5 = (ii > 1) ? S1[ii - 1] : -1;
      for (u = maxdist; u > ii; u--) {
        if (with_gquad) {
          if (fij == ggg[ii][u - ii] + f3[u + 1]) {
            *i  = *j = 0;
            *k  = u + 1;
            return vrna_bt_gquad_mfe(fc, ii, u, bp_stack);
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
            *i                            = ii;
            *j                            = u;
            *k                            = u + 1;
            vrna_bps_push(bp_stack,
                          (vrna_bp_t){
                            .i = ii,
                            .j = u
                          });
            return 1;
          }
        }
      }
      break;

    default:
      mm5 = S1[ii];
      for (u = maxdist; u > ii; u--) {
        if (with_gquad) {
          if (fij == ggg[ii][u - ii] + f3[u + 1]) {
            *i  = *j = 0;
            *k  = u + 1;
            return vrna_bt_gquad_mfe(fc, ii, u, bp_stack);
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
              *i                            = ii + 1;
              *j                            = u;
              *k                            = u + 2;
              vrna_bps_push(bp_stack,
                            (vrna_bp_t){
                              .i = ii + 1,
                              .j = u
                            });
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
              *i                            = ii + 1;
              *j                            = u;
              *k                            = u + 2;
              vrna_bps_push(bp_stack,
                            (vrna_bp_t){
                              .i = ii + 1,
                              .j = u
                            });
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
            *i                            = ii + 1;
            *j                            = u;
            *k                            = u + 1;
            vrna_bps_push(bp_stack,
                          (vrna_bp_t){
                            .i = ii + 1,
                            .j = u
                          });
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
              *i                            = ii;
              *j                            = u;
              *k                            = u + 2;
              vrna_bps_push(bp_stack,
                            (vrna_bp_t){
                              .i = ii,
                              .j = u
                            });
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
              *i                            = ii;
              *j                            = u;
              *k                            = u + 2;
              vrna_bps_push(bp_stack,
                            (vrna_bp_t){
                              .i = ii,
                              .j = u
                            });
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
            *i                            = ii;
            *j                            = u;
            *k                            = u + 1;
              vrna_bps_push(bp_stack,
                            (vrna_bp_t){
                              .i = ii,
                              .j = u
                            });
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
                           unsigned int                  *k,
                           unsigned int                  maxdist,
                           unsigned int                  *i,
                           unsigned int                  *j,
                           vrna_bps_t           bp_stack)
{
  short                     **S, **S5, **S3;
  unsigned int              n, type, ss, n_seq, **a2s, ii, u, dangle_model, with_gquad;
  int                       fij, cc, fj, *f3, **c, **ggg;
  vrna_param_t              *P;
  vrna_md_t                 *md;
  vrna_sc_t                 **scs;
  vrna_hc_eval_f evaluate;
  struct hc_ext_def_dat     hc_dat_local;

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

  ii = *k;

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

    if (++ii > maxdist)
      break;
  } while (fij == fj);
  ii--;

  if (ii >= maxdist) {
    /* no more pairs */
    *i  = *j = 0;
    *k  = 0;
    return 1;
  }

  /*
   *  must have found a decomposition
   *  i is paired. Find pairing partner
   */
  switch (dangle_model) {
    /* no dangles */
    case 0:
      for (u = maxdist; u > ii; u--) {
        if (with_gquad) {
          if (fij == ggg[ii][u - ii] + f3[u + 1]) {
            *i  = *j = 0;
            *k  = u + 1;
            return vrna_bt_gquad_mfe(fc, ii, u, bp_stack);
          }
        }

        if (evaluate(ii, n, u, u + 1, VRNA_DECOMP_EXT_STEM_EXT, &hc_dat_local)) {
          cc = c[ii][u - ii];
          for (ss = 0; ss < n_seq; ss++) {
            type  = vrna_get_ptype_md(S[ss][ii], S[ss][u], md);
            cc    += vrna_E_exterior_stem(type, -1, -1, P);
          }

          if (fij == cc + f3[u + 1]) {
            *i                            = ii;
            *j                            = u;
            *k                            = u + 1;
            vrna_bps_push(bp_stack,
                          (vrna_bp_t){
                            .i = ii,
                            .j = u
                          });
            return 1;
          }
        }
      }
      break;

    /* dangles on both sides */
    case 2:
      for (u = maxdist; u > ii; u--) {
        if (with_gquad) {
          if (fij == ggg[ii][u - ii] + f3[u + 1]) {
            *i  = *j = 0;
            *k  = u + 1;
            return vrna_bt_gquad_mfe(fc, ii, u, bp_stack);
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
            *i                            = ii;
            *j                            = u;
            *k                            = u + 1;
            vrna_bps_push(bp_stack,
                          (vrna_bp_t){
                            .i = ii,
                            .j = u
                          });
            return 1;
          }
        }
      }
      break;
  }

  return 0;
}


PRIVATE int
bt_ext_loop_f3_pp(vrna_fold_compound_t  *fc,
                  unsigned int                   *i,
                  unsigned int                   maxj)
{
  int j;
  unsigned int start;

  j     = -1;
  start = *i;

  if (fc) {
    char                      **ptype;
    short                     *S1;
    unsigned int              ii, length, type, traced2, dangle_model, with_gquad, maxdist;
    int                       cc, **c, **ggg, *f3, fij;
    vrna_param_t              *P;
    vrna_md_t                 *md;
    vrna_sc_t                 *sc;
    vrna_hc_eval_f evaluate;
    struct hc_ext_def_dat     hc_dat_local;
#ifdef VRNA_WITH_SVM
    int                       zsc_pre_filter;
    vrna_zsc_dat_t            zsc_data;
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

    if (!traced2)
      j = -1;
  }

  *i = start;
  return j;
}


PRIVATE int
bt_ext_loop_f3_pp_comparative(vrna_fold_compound_t  *fc,
                              unsigned int                   *i,
                              unsigned int                   maxj)
{
  int j;
  unsigned int start;

  j = -1;

  if (fc) {
    short                     **S, **S5, **S3;
    unsigned int              tt, s, n_seq, **a2s, traced2, length, dangle_model, with_gquad, maxdist;
    int                       cc, **c, **ggg, *f3, fij;
    vrna_param_t              *P;
    vrna_md_t                 *md;
    vrna_sc_t                 **scs;
    vrna_hc_eval_f evaluate;
    struct hc_ext_def_dat     hc_dat_local;

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
            cc +=scs[s]->f(start,
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
              for (s = 0; s < n_seq; s++) {
                if ((scs[s]) && (scs[s]->f))
                  cc += scs[s]->f(start, length, j, j + 1, VRNA_DECOMP_EXT_STEM_EXT, scs[s]->data);
              }
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
              for (s = 0; s < n_seq; s++) {
                if ((scs[s]) && (scs[s]->f))
                  cc += scs[s]->f(start, length, start, j, VRNA_DECOMP_EXT_STEM, scs[s]->data);
              }
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
              for (s = 0; s < n_seq; s++) {
                if ((scs[s]) && (scs[s]->f))
                  cc += scs[s]->f(start, length, j, j + 1, VRNA_DECOMP_EXT_STEM_EXT, scs[s]->data);
              }
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
              for (s = 0; s < n_seq; s++) {
                if ((scs[s]) && (scs[s]->f))
                  cc += scs[s]->f(start, length, start, j, VRNA_DECOMP_EXT_STEM, scs[s]->data);
              }
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

    if (!traced2)
      j = -1;
  }

  return j;
}


#ifndef VRNA_DISABLE_BACKWARD_COMPATIBILITY


PUBLIC int
vrna_BT_ext_loop_f5(vrna_fold_compound_t  *fc,
                    int                   *k,
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
        r = (int)bt_ext_loop_f5(fc, (unsigned int)(*k), bps, bts);
        break;

      case VRNA_FC_TYPE_COMPARATIVE:
        r = (int)bt_ext_loop_f5_comparative(fc, (unsigned int)(*k), bps, bts);
        break;
    }

    while (vrna_bts_size(bts) > 0) {
      vrna_sect_t s = vrna_bts_pop(bts);
      if (s.ml == VRNA_MX_FLAG_F5) {
        *k = (int)s.j;
      } else if (s.ml == VRNA_MX_FLAG_C) {
        bp_stack[++(*stack_count)].i = s.i;
        bp_stack[*stack_count].j = s.j;
        *i = (int)s.i;
        *j = (int)s.j;
      } else if (s.ml == VRNA_MX_FLAG_G) {
        r = vrna_bt_gquad_mfe(fc, s.i, s.j, bps);
      }
    }

    while (vrna_bps_size(bps) > 0) {
      vrna_bp_t bp = vrna_bps_pop(bps);
      bp_stack[++(*stack_count)].i = bp.i;
      bp_stack[*stack_count].j = bp.j;
    }

    vrna_bps_free(bps);
    vrna_bts_free(bts);
  }

  return r;
}


PUBLIC int
vrna_BT_ext_loop_f3_pp(vrna_fold_compound_t *fc,
                       int                  *i,
                       int                  maxj)
{
  int r;

  unsigned int ii = (unsigned int)(*i);

  r = vrna_bt_ext_loop_f3_pp(fc, &ii, (unsigned int)maxj);

  *i = (int)ii;

  return r;
}


#endif
