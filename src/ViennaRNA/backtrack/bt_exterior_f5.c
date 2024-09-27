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


/*
 #################################
 # BEGIN OF FUNCTION DEFINITIONS #
 #################################
 */
PUBLIC unsigned int
vrna_bt_exterior_f5(vrna_fold_compound_t  *fc,
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


PRIVATE unsigned int
bt_ext_loop_f5(vrna_fold_compound_t *fc,
               unsigned int         j,
               vrna_bps_t           bp_stack VRNA_UNUSED,
               vrna_bts_t           bt_stack)
{
  char                  *ptype;
  short                 mm5, mm3, *S1;
  unsigned int          length, *sn, type, u, dangle_model, with_gquad, cnt,
                        ii, with_ud;
  int                   fij, fi, en, e, *my_f5, *my_c, *idx, e_gq;
  vrna_param_t          *P;
  vrna_md_t             *md;
  vrna_sc_t             *sc;
  vrna_ud_t             *domains_up;
  vrna_hc_eval_f        evaluate;
  struct hc_ext_def_dat hc_dat_local;
  vrna_smx_csr(int)     *c_gq;

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
      for (cnt = 0; cnt < (unsigned int)domains_up->uniq_motif_count; cnt++) {
        u = domains_up->uniq_motif_size[cnt];
        if ((j >= u) &&
            evaluate(1, j, 1, j - u, VRNA_DECOMP_EXT_EXT, &hc_dat_local)) {
          ii  = j - u + 1;
          en  = domains_up->energy_cb(fc,
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
        e_gq = vrna_smx_csr_get(c_gq, 1, j, INF);
#else
        e_gq = vrna_smx_csr_int_get(c_gq, 1, j, INF);
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
                           vrna_bps_t           bp_stack VRNA_UNUSED,
                           vrna_bts_t           bt_stack)
{
  unsigned int          **a2s, n;
  short                 **S, **S5, **S3;
  unsigned int          tt, u, dangle_model, with_gquad, n_seq, ss;
  int                   fij, fi, en, e_gq, *my_f5, *my_c, *idx, mm5, mm3;
  vrna_param_t          *P;
  vrna_md_t             *md;
  vrna_sc_t             **scs;
  vrna_hc_eval_f        evaluate;
  struct hc_ext_def_dat hc_dat_local;

  vrna_smx_csr(int) * c_gq;

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
