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
BT_mb_loop(vrna_fold_compound_t *fc,
           int                  *i,
           int                  *j,
           int                  *k,
           int                  en,
           int                  *component1,
           int                  *component2);


PRIVATE int
BT_mb_loop_split(vrna_fold_compound_t *fc,
                 int                  *i,
                 int                  *j,
                 int                  *k,
                 int                  *l,
                 int                  *component1,
                 int                  *component2,
                 vrna_bp_stack_t      *bp_stack,
                 int                  *stack_count);


/*
 #################################
 # BEGIN OF FUNCTION DEFINITIONS #
 #################################
 */
PUBLIC int
vrna_BT_mb_loop_split(vrna_fold_compound_t  *fc,
                      int                   *i,
                      int                   *j,
                      int                   *k,
                      int                   *l,
                      int                   *c1,
                      int                   *c2,
                      vrna_bp_stack_t       *bp_stack,
                      int                   *stack_count)
{
  int r = 0;

  if (fc)
    r = BT_mb_loop_split(fc, i, j, k, l, c1, c2, bp_stack, stack_count);

  return r;
}


PUBLIC int
vrna_BT_mb_loop(vrna_fold_compound_t  *fc,
                int                   *i,
                int                   *j,
                int                   *k,
                int                   en,
                int                   *c1,
                int                   *c2)
{
  int r = 0;

  if (fc)
    r = BT_mb_loop(fc, i, j, k, en, c1, c2);

  return r;
}


PUBLIC int
vrna_BT_mb_loop_fake(vrna_fold_compound_t *fc,
                     int                  *u,
                     int                  *i,
                     int                  *j,
                     vrna_bp_stack_t      *bp_stack,
                     int                  *stack_count)
{
  char                      *ptype;
  short                     mm5, mm3, *S1;
  unsigned int              *sn, *se;
  int                       length, ii, jj, k, en, fij, fi, *my_c, *my_fc, *my_ggg,
                            *idx, with_gquad, dangle_model, turn, type;
  vrna_param_t              *P;
  vrna_md_t                 *md;
  vrna_sc_t                 *sc;
  vrna_callback_hc_evaluate *evaluate;
  struct hc_mb_def_dat      hc_dat_local;

  length        = fc->length;
  P             = fc->params;
  md            = &(P->model_details);
  sn            = fc->strand_number;
  se            = fc->strand_end;
  sc            = fc->sc;
  S1            = fc->sequence_encoding;
  ptype         = fc->ptype;
  idx           = fc->jindx;
  my_c          = fc->matrices->c;
  my_fc         = fc->matrices->fc;
  my_ggg        = fc->matrices->ggg;
  turn          = md->min_loop_size;
  with_gquad    = md->gquad;
  dangle_model  = md->dangles;
  evaluate      = prepare_hc_mb_def_ext(fc, &hc_dat_local);

  ii  = *i;
  jj  = *j;

  if (ii <= se[0]) {
    /* 'lower' part (fc[i<cut,j=cut-1]) */

    /* nibble off unpaired 5' bases */
    do {
      fij = my_fc[ii];
      fi  = INF;

      if (evaluate(ii, jj, ii + 1, jj, VRNA_DECOMP_EXT_EXT, &hc_dat_local)) {
        fi = my_fc[ii + 1];

        if (sc)
          if (sc->energy_up)
            fi += sc->energy_up[ii][1];
      }

      if (++ii == jj)
        break;
    } while (fij == fi);
    ii--;

    if (jj < ii + turn + 2) {
      /* no more pairs */
      *u = *i = *j = -1;
      return 1;
    }

    mm5 = (ii > 1 && (sn[ii - 1] == sn[ii])) ? S1[ii - 1] : -1;

    /* i or i+1 is paired. Find pairing partner */
    switch (dangle_model) {
      case 0:
        for (k = ii + turn + 1; k <= jj; k++) {
          if (evaluate(ii, jj, k, k + 1, VRNA_DECOMP_EXT_STEM_EXT, &hc_dat_local)) {
            type = vrna_get_ptype(idx[k] + ii, ptype);

            if (fij == my_fc[k + 1] + my_c[idx[k] + ii] + vrna_E_ext_stem(type, -1, -1, P)) {
              bp_stack[++(*stack_count)].i  = ii;
              bp_stack[(*stack_count)].j    = k;
              *u                            = k + 1;
              *i                            = ii;
              *j                            = k;
              return 1;
            }
          }

          if (with_gquad) {
            if (fij == my_fc[k + 1] + my_ggg[idx[k] + ii]) {
              *u  = k + 1;
              *i  = *j = -1;
              return vrna_BT_gquad_mfe(fc, ii, k, bp_stack, stack_count);
            }
          }
        }
        break;

      case 2:
        for (k = ii + turn + 1; k <= jj; k++) {
          if (evaluate(ii, jj, k, k + 1, VRNA_DECOMP_EXT_STEM_EXT, &hc_dat_local)) {
            mm3   = (sn[k] == sn[k + 1]) ? S1[k + 1] : -1;
            type  = vrna_get_ptype(idx[k] + ii, ptype);

            if (fij == my_fc[k + 1] + my_c[idx[k] + ii] + vrna_E_ext_stem(type, mm5, mm3, P)) {
              bp_stack[++(*stack_count)].i  = ii;
              bp_stack[(*stack_count)].j    = k;
              *u                            = k + 1;
              *i                            = ii;
              *j                            = k;
              return 1;
            }
          }

          if (with_gquad) {
            if (fij == my_fc[k + 1] + my_ggg[idx[k] + ii]) {
              *u  = k + 1;
              *i  = *j = -1;
              return vrna_BT_gquad_mfe(fc, ii, k, bp_stack, stack_count);
            }
          }
        }
        break;

      default:
        for (k = ii + turn + 1; k <= jj; k++) {
          type = vrna_get_ptype(idx[k] + ii, ptype);
          if (evaluate(ii, jj, k, k + 1, VRNA_DECOMP_EXT_STEM_EXT, &hc_dat_local)) {
            if (fij == my_fc[k + 1] + my_c[idx[k] + ii] + vrna_E_ext_stem(type, -1, -1, P)) {
              bp_stack[++(*stack_count)].i  = ii;
              bp_stack[(*stack_count)].j    = k;
              *u                            = k + 1;
              *i                            = ii;
              *j                            = k;
              return 1;
            }
          }

          if (evaluate(ii, jj, k, k + 2, VRNA_DECOMP_EXT_STEM_EXT, &hc_dat_local)) {
            mm3 = (sn[k] == sn[k + 1]) ? S1[k + 1] : -1;
            en  = my_c[idx[k] + ii];
            if (sc)
              if (sc->energy_up)
                en += sc->energy_up[k + 1][1];

            if (fij == my_fc[k + 2] + en + vrna_E_ext_stem(type, -1, mm3, P)) {
              bp_stack[++(*stack_count)].i  = ii;
              bp_stack[(*stack_count)].j    = k;
              *u                            = k + 2;
              *i                            = ii;
              *j                            = k;
              return 1;
            }
          }

          if (with_gquad) {
            if (fij == my_fc[k + 1] + my_ggg[idx[k] + ii]) {
              *u  = k + 1;
              *i  = *j = -1;
              return vrna_BT_gquad_mfe(fc, ii, k, bp_stack, stack_count);
            }
          }

          if (evaluate(ii, jj, k, k + 1, VRNA_DECOMP_EXT_STEM_EXT1, &hc_dat_local)) {
            mm5   = (sn[ii] == sn[ii + 1]) ? S1[ii] : -1;
            mm3   = (sn[k] == sn[k + 1]) ? S1[k + 1] : -1;
            type  = vrna_get_ptype(idx[k] + ii + 1, ptype);

            en = my_c[idx[k] + ii + 1];
            if (sc)
              if (sc->energy_up)
                en += sc->energy_up[ii][1];

            if (fij == en + my_fc[k + 1] + vrna_E_ext_stem(type, mm5, -1, P)) {
              bp_stack[++(*stack_count)].i  = ii + 1;
              bp_stack[(*stack_count)].j    = k;
              *u                            = k + 1;
              *i                            = ii + 1;
              *j                            = k;
              return 1;
            }
          }

          if ((k < jj) && (evaluate(ii, jj, k, k + 2, VRNA_DECOMP_EXT_STEM_EXT1, &hc_dat_local))) {
            mm5   = (sn[ii] == sn[ii + 1]) ? S1[ii] : -1;
            mm3   = (sn[k] == sn[k + 1]) ? S1[k + 1] : -1;
            type  = vrna_get_ptype(idx[k] + ii + 1, ptype);

            en = my_c[idx[k] + ii + 1];
            if (sc)
              if (sc->energy_up)
                en += sc->energy_up[k + 1][1];

            if (fij == en + my_fc[k + 2] + vrna_E_ext_stem(type, mm5, mm3, P)) {
              bp_stack[++(*stack_count)].i  = ii + 1;
              bp_stack[(*stack_count)].j    = k;
              *u                            = k + 2;
              *i                            = ii + 1;
              *j                            = k;
              return 1;
            }
          }
        }
        break;
    }
  } else {
    /* 'upper' part (fc[i=cut,j>cut]) */

    /* nibble off unpaired 3' bases */
    do {
      fij = my_fc[jj];
      fi  = INF;

      if (evaluate(ii, jj, ii, jj - 1, VRNA_DECOMP_EXT_EXT, &hc_dat_local)) {
        fi = my_fc[jj - 1];

        if (sc)
          if (sc->energy_up)
            fi += sc->energy_up[jj][1];
      }

      if (--jj == ii)
        break;
    } while (fij == fi);
    jj++;

    if (jj < ii + turn + 2) {
      /* no more pairs */
      *u = *i = *j = -1;
      return 1;
    }

    /* j or j-1 is paired. Find pairing partner */
    mm3 = ((jj < length) && (sn[jj] == sn[jj + 1])) ? S1[jj + 1] : -1;
    switch (dangle_model) {
      case 0:
        for (k = jj - turn - 1; k >= ii; k--) {
          if (with_gquad) {
            if (fij == my_fc[k - 1] + my_ggg[idx[jj] + k]) {
              *u  = k - 1;
              *i  = *j = -1;
              return vrna_BT_gquad_mfe(fc, k, jj, bp_stack, stack_count);
            }
          }

          if (evaluate(ii, jj, k - 1, k, VRNA_DECOMP_EXT_EXT_STEM, &hc_dat_local)) {
            type = vrna_get_ptype(idx[jj] + k, ptype);

            en = my_c[idx[jj] + k];
            if (sn[k] != sn[jj])
              en += P->DuplexInit;

            if (fij == my_fc[k - 1] + en + vrna_E_ext_stem(type, -1, -1, P)) {
              bp_stack[++(*stack_count)].i  = k;
              bp_stack[(*stack_count)].j    = jj;
              *u                            = k - 1;
              *i                            = k;
              *j                            = jj;
              return 1;
            }
          }
        }
        break;

      case 2:
        for (k = jj - turn - 1; k >= ii; k--) {
          if (with_gquad) {
            if (fij == my_fc[k - 1] + my_ggg[idx[jj] + k]) {
              *u  = k - 1;
              *i  = *j = -1;
              return vrna_BT_gquad_mfe(fc, k, jj, bp_stack, stack_count);
            }
          }

          if (evaluate(ii, jj, k - 1, k, VRNA_DECOMP_EXT_EXT_STEM, &hc_dat_local)) {
            mm5   = ((k > 1) && (sn[k - 1] == sn[k])) ? S1[k - 1] : -1;
            type  = vrna_get_ptype(idx[jj] + k, ptype);

            en = my_c[idx[jj] + k];
            if (sn[k] != sn[jj])
              en += P->DuplexInit;

            if (fij == my_fc[k - 1] + en + vrna_E_ext_stem(type, mm5, mm3, P)) {
              bp_stack[++(*stack_count)].i  = k;
              bp_stack[(*stack_count)].j    = jj;
              *u                            = k - 1;
              *i                            = k;
              *j                            = jj;
              return 1;
            }
          }
        }
        break;

      default:
        for (k = jj - turn - 1; k >= ii; k--) {
          if (with_gquad) {
            if (fij == my_fc[k - 1] + my_ggg[idx[jj] + k]) {
              *u  = k - 1;
              *i  = *j = -1;
              return vrna_BT_gquad_mfe(fc, k, jj, bp_stack, stack_count);
            }
          }

          if (evaluate(ii, jj, k - 1, k, VRNA_DECOMP_EXT_EXT_STEM, &hc_dat_local)) {
            type = vrna_get_ptype(idx[jj] + k, ptype);

            en = my_c[idx[jj] + k];
            if (sn[k] != sn[jj])
              en += P->DuplexInit;

            if (fij == my_fc[k - 1] + en + vrna_E_ext_stem(type, -1, -1, P)) {
              bp_stack[++(*stack_count)].i  = k;
              bp_stack[(*stack_count)].j    = jj;
              *u                            = k - 1;
              *i                            = k;
              *j                            = jj;
              return 1;
            }
          }

          if ((k > 1) && (sn[k - 1] == sn[k]) &&
              evaluate(ii, jj, k - 2, k, VRNA_DECOMP_EXT_EXT_STEM, &hc_dat_local)) {
            type = vrna_get_ptype(idx[jj] + k, ptype);

            en = my_c[idx[jj] + k];
            if (sn[k] != sn[jj])
              en += P->DuplexInit;

            mm5 = S1[k - 1];

            if (sc)
              if (sc->energy_up)
                en += sc->energy_up[k - 1][1];

            if (fij == my_fc[k - 2] + en + vrna_E_ext_stem(type, mm5, -1, P)) {
              bp_stack[++(*stack_count)].i  = k;
              bp_stack[(*stack_count)].j    = jj;
              *u                            = k - 2;
              *i                            = k;
              *j                            = jj;
              return 1;
            }
          }

          if ((sn[jj - 1] == sn[jj]) &&
              evaluate(ii, jj, k - 1, k, VRNA_DECOMP_EXT_EXT_STEM1, &hc_dat_local)) {
            type = vrna_get_ptype(idx[jj - 1] + k, ptype);

            mm3 = S1[jj];
            en  = my_c[idx[jj - 1] + k];
            if (sn[k] != sn[jj - 1])
              en += P->DuplexInit;         /* ??? */

            if (sc)
              if (sc->energy_up)
                en += sc->energy_up[jj][1];

            if (fij == en + my_fc[k - 1] + vrna_E_ext_stem(type, -1, mm3, P)) {
              bp_stack[++(*stack_count)].i  = k;
              bp_stack[(*stack_count)].j    = jj - 1;
              *u                            = k - 1;
              *i                            = k;
              *j                            = jj - 1;
              return 1;
            }
          }

          if ((k > ii) && (sn[jj - 1] == sn[jj]) &&
              evaluate(ii, jj, k - 2, k, VRNA_DECOMP_EXT_EXT_STEM1, &hc_dat_local)) {
            type = vrna_get_ptype(idx[jj - 1] + k, ptype);

            mm3 = S1[jj];
            en  = my_c[idx[jj - 1] + k];
            if (sn[k] != sn[jj - 1])
              en += P->DuplexInit;         /* ??? */

            mm5 = (sn[k - 1] == sn[k]) ? S1[k - 1] : -1;
            if (sc)
              if (sc->energy_up)
                en += sc->energy_up[k - 1][1];

            if (fij == my_fc[k - 2] + en + vrna_E_ext_stem(type, mm5, mm3, P)) {
              bp_stack[++(*stack_count)].i  = k;
              bp_stack[(*stack_count)].j    = jj - 1;
              *u                            = k - 2;
              *i                            = k;
              *j                            = jj - 1;
              return 1;
            }
          }
        }
        break;
    }
  }

  return 0;
}


PRIVATE int
BT_mb_loop_split(vrna_fold_compound_t *fc,
                 int                  *i,
                 int                  *j,
                 int                  *k,
                 int                  *l,
                 int                  *component1,
                 int                  *component2,
                 vrna_bp_stack_t      *bp_stack,
                 int                  *stack_count)
{
  unsigned char             sliding_window;
  char                      *ptype, **ptype_local;
  short                     *S1, **SS, **S5, **S3;
  unsigned int              n_seq, s;
  int                       ij, ii, jj, fij, fi, u, en, *my_c, *my_fML, *my_ggg,
                            turn, *idx, with_gquad, dangle_model, *rtype, kk, cnt,
                            with_ud, type, type_2, en2, **c_local, **fML_local, **ggg_local;
  vrna_param_t              *P;
  vrna_md_t                 *md;
  vrna_ud_t                 *domains_up;
  vrna_callback_hc_evaluate *evaluate;
  struct hc_mb_def_dat      hc_dat_local;
  struct sc_mb_dat          sc_wrapper;

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
  my_ggg    = (sliding_window) ? NULL : fc->matrices->ggg;
  c_local   = (sliding_window) ? fc->matrices->c_local : NULL;
  fML_local = (sliding_window) ? fc->matrices->fML_local : NULL;
  ggg_local = (sliding_window) ? fc->matrices->ggg_local : NULL;


  turn          = md->min_loop_size;
  with_gquad    = md->gquad;
  with_ud       = (domains_up && domains_up->energy_cb) ? 1 : 0;
  dangle_model  = md->dangles;
  evaluate      = prepare_hc_mb_def(fc, &hc_dat_local);

  init_sc_mb(fc, &sc_wrapper);

  ii  = *i;
  jj  = *j;

  if (with_ud) {
    /* nibble off unpaired stretches at 3' site */
    do {
      fij = (sliding_window) ? fML_local[ii][jj - ii] : my_fML[idx[jj] + ii];
      fi  = INF;

      /* process regular unpaired nucleotides (unbound by ligand) first */
      if (evaluate(ii, jj, ii, jj - 1, VRNA_DECOMP_ML_ML, &hc_dat_local)) {
        fi = P->MLbase *
             n_seq;
        fi += (sliding_window) ? fML_local[ii][jj - 1 - ii] : my_fML[idx[jj - 1] + ii];

        if (sc_wrapper.red_ml)
          fi += sc_wrapper.red_ml(ii, jj, ii, jj - 1, &sc_wrapper);

        if (jj == ii)
          return 0; /* no more pairs */

        if (fij == fi) {
          jj--;
          continue;
        }
      }

      /* next try to nibble off ligand */
      for (cnt = 0; cnt < domains_up->uniq_motif_count; cnt++) {
        u   = domains_up->uniq_motif_size[cnt];
        kk  = jj - u + 1;
        if ((kk >= ii) && evaluate(ii, jj, ii, jj - u, VRNA_DECOMP_ML_ML, &hc_dat_local)) {
          en = domains_up->energy_cb(fc,
                                     kk,
                                     jj,
                                     VRNA_UNSTRUCTURED_DOMAIN_MB_LOOP | VRNA_UNSTRUCTURED_DOMAIN_MOTIF,
                                     domains_up->data);

          if (sc_wrapper.red_ml)
            en += sc_wrapper.red_ml(ii, jj, ii, jj - u, &sc_wrapper);

          fi = en +
               u * P->MLbase *
               n_seq;
          fi += (sliding_window) ? fML_local[ii][kk - 1 - ii] : my_fML[idx[kk - 1] + ii];

          if (fij == fi) {
            /* skip remaining motifs after first hit */
            jj = kk - 1;
            break;
          }
        }
      }

      if (jj < ii)
        return 0; /* no more pairs */
    } while (fij == fi);

    /* nibble off unpaired stretches at 5' site */
    do {
      fij = (sliding_window) ? fML_local[ii][jj - ii] : my_fML[idx[jj] + ii];
      fi  = INF;

      /* again, process regular unpaired nucleotides (unbound by ligand) first */
      if (evaluate(ii, jj, ii + 1, jj, VRNA_DECOMP_ML_ML, &hc_dat_local)) {
        fi = P->MLbase *
             n_seq;
        fi += (sliding_window) ? fML_local[ii + 1][jj - (ii + 1)] : my_fML[idx[jj] + ii + 1];

        if (sc_wrapper.red_ml)
          fi += sc_wrapper.red_ml(ii, jj, ii + 1, jj, &sc_wrapper);

        if (ii + 1 == jj)
          return 0; /* no more pairs */

        if (fij == fi) {
          ii++;
          continue;
        }
      }

      /* next try to nibble off ligand again */
      for (cnt = 0; cnt < domains_up->uniq_motif_count; cnt++) {
        u   = domains_up->uniq_motif_size[cnt];
        kk  = ii + u - 1;
        if ((kk <= jj) && evaluate(ii, jj, ii + u, jj, VRNA_DECOMP_ML_ML, &hc_dat_local)) {
          en = domains_up->energy_cb(fc,
                                     ii,
                                     kk,
                                     VRNA_UNSTRUCTURED_DOMAIN_MB_LOOP | VRNA_UNSTRUCTURED_DOMAIN_MOTIF,
                                     domains_up->data);

          if (sc_wrapper.red_ml)
            en += sc_wrapper.red_ml(ii, jj, ii + u, jj, &sc_wrapper);

          fi = en +
               u * P->MLbase *
               n_seq;

          fi += (sliding_window) ? fML_local[kk + 1][jj - (kk + 1)] : my_fML[idx[jj] + kk + 1];

          if (fij == fi) {
            /* skip remaining motifs after first hit */
            ii = kk + 1;
            break;
          }
        }
      }

      if (ii > jj)
        return 0; /* no more pairs */
    } while (fij == fi);
  } else {
    /* nibble off unpaired 3' bases */
    do {
      fij = (sliding_window) ? fML_local[ii][jj - ii] : my_fML[idx[jj] + ii];
      fi  = INF;

      if (evaluate(ii, jj, ii, jj - 1, VRNA_DECOMP_ML_ML, &hc_dat_local)) {
        fi = P->MLbase *
             n_seq;
        fi += (sliding_window) ? fML_local[ii][jj - 1 - ii] : my_fML[idx[jj - 1] + ii];

        if (sc_wrapper.red_ml)
          fi += sc_wrapper.red_ml(ii, jj, ii, jj - 1, &sc_wrapper);
      }

      if (--jj == 0)
        break;
    } while (fij == fi);
    jj++;

    /* nibble off unpaired 5' bases */
    do {
      fij = (sliding_window) ? fML_local[ii][jj - ii] : my_fML[idx[jj] + ii];
      fi  = INF;

      if (evaluate(ii, jj, ii + 1, jj, VRNA_DECOMP_ML_ML, &hc_dat_local)) {
        fi = P->MLbase *
             n_seq;
        fi += (sliding_window) ? fML_local[ii + 1][jj - (ii + 1)] : my_fML[idx[jj] + ii + 1];

        if (sc_wrapper.red_ml)
          fi += sc_wrapper.red_ml(ii, jj, ii + 1, jj, &sc_wrapper);
      }

      if (++ii == jj)
        break;
    } while (fij == fi);
    ii--;

    if (jj < ii + turn + 1) /* no more pairs */
      return 0;
  }

  ij = (sliding_window) ? 0 : idx[jj] + ii;

  *component1 = *component2 = 1; /* split into two multi loop parts by default */

  /* 1. test for single component */

  if (with_gquad) {
    en = E_MLstem(0, -1, -1, P) *
         n_seq;
    en += (sliding_window) ? ggg_local[ii][jj - ii] : my_ggg[ij];

    if (fij == en) {
      *i  = *j = -1;
      *k  = *l = -1;
      return vrna_BT_gquad_mfe(fc, ii, jj, bp_stack, stack_count);
    }
  }

  en = (sliding_window) ? c_local[ii][jj - ii] : my_c[ij];

  if (sc_wrapper.red_stem)
    en += sc_wrapper.red_stem(ii, jj, ii, jj, &sc_wrapper);

  switch (dangle_model) {
    case 0:
      if (evaluate(ii, jj, ii, jj, VRNA_DECOMP_ML_STEM, &hc_dat_local)) {
        en2 = 0;
        switch (fc->type) {
          case VRNA_FC_TYPE_SINGLE:
            type = (sliding_window) ? vrna_get_ptype_window(ii, jj, ptype_local) : vrna_get_ptype(
              ij,
              ptype);
            en2 = E_MLstem(type, -1, -1, P);
            break;

          case VRNA_FC_TYPE_COMPARATIVE:
            for (s = 0; s < n_seq; s++) {
              type  = vrna_get_ptype_md(SS[s][ii], SS[s][jj], md);
              en2   += E_MLstem(type, -1, -1, P);
            }
            break;
        }

        if (fij == en + en2) {
          *i          = *j = -1;
          *k          = ii;
          *l          = jj;
          *component2 = 2;          /* 2nd part is structure enclosed by base pair */
          return 1;
        }
      }

      break;

    case 2:
      if (evaluate(ii, jj, ii, jj, VRNA_DECOMP_ML_STEM, &hc_dat_local)) {
        en2 = 0;
        switch (fc->type) {
          case VRNA_FC_TYPE_SINGLE:
            type = (sliding_window) ? vrna_get_ptype_window(ii, jj, ptype_local) : vrna_get_ptype(
              ij,
              ptype);
            en2 = E_MLstem(type, S1[ii - 1], S1[jj + 1], P);
            break;

          case VRNA_FC_TYPE_COMPARATIVE:
            for (s = 0; s < n_seq; s++) {
              type  = vrna_get_ptype_md(SS[s][ii], SS[s][jj], md);
              en2   += E_MLstem(type, S5[s][ii], S3[s][jj], P);
            }
            break;
        }

        if (fij == en + en2) {
          *i          = *j = -1;
          *k          = ii;
          *l          = jj;
          *component2 = 2;
          return 1;
        }
      }

      break;

    default:
      if (evaluate(ii, jj, ii, jj, VRNA_DECOMP_ML_STEM, &hc_dat_local)) {
        en2 = 0;
        switch (fc->type) {
          case VRNA_FC_TYPE_SINGLE:
            type = (sliding_window) ? vrna_get_ptype_window(ii, jj, ptype_local) : vrna_get_ptype(
              ij,
              ptype);
            en2 = E_MLstem(type, -1, -1, P);
            break;

          case VRNA_FC_TYPE_COMPARATIVE:
            for (s = 0; s < n_seq; s++) {
              type  = vrna_get_ptype_md(SS[s][ii], SS[s][jj], md);
              en2   += E_MLstem(type, -1, -1, P);
            }
            break;
        }

        if (fij == en + en2) {
          *i          = *j = -1;
          *k          = ii;
          *l          = jj;
          *component2 = 2;
          return 1;
        }
      }

      if (evaluate(ii, jj, ii + 1, jj, VRNA_DECOMP_ML_STEM, &hc_dat_local)) {
        en2 = P->MLbase *
              n_seq;
        en2 += (sliding_window) ? c_local[ii + 1][jj - (ii + 1)] : my_c[ij + 1];

        switch (fc->type) {
          case VRNA_FC_TYPE_SINGLE:
            type =
              (sliding_window) ? vrna_get_ptype_window(ii + 1, jj, ptype_local) : vrna_get_ptype(
                ij + 1,
                ptype);
            en2 += E_MLstem(type, S1[ii], -1, P);
            break;

          case VRNA_FC_TYPE_COMPARATIVE:
            for (s = 0; s < n_seq; s++) {
              type  = vrna_get_ptype_md(SS[s][ii], SS[s][jj], md);
              en2   += E_MLstem(type, S5[s][ii], -1, P);
            }
            break;
        }

        if (sc_wrapper.red_stem)
          en2 += sc_wrapper.red_stem(ii, jj, ii + 1, jj, &sc_wrapper);

        if (fij == en2) {
          *i          = *j = -1;
          *k          = ii + 1;
          *l          = jj;
          *component2 = 2;
          return 1;
        }
      }

      if (evaluate(ii, jj, ii, jj - 1, VRNA_DECOMP_ML_STEM, &hc_dat_local)) {
        en2 = P->MLbase *
              n_seq;
        en2 += (sliding_window) ? c_local[ii][jj - 1 - ii] : my_c[idx[jj - 1] + ii];

        switch (fc->type) {
          case VRNA_FC_TYPE_SINGLE:
            type =
              (sliding_window) ? vrna_get_ptype_window(ii, jj - 1,
                                                       ptype_local) : vrna_get_ptype(
                idx[jj - 1] + ii,
                ptype);
            en2 += E_MLstem(type, -1, S1[jj], P);
            break;

          case VRNA_FC_TYPE_COMPARATIVE:
            for (s = 0; s < n_seq; s++) {
              type  = vrna_get_ptype_md(SS[s][ii], SS[s][jj], md);
              en2   += E_MLstem(type, -1, S3[s][jj], P);
            }
            break;
        }

        if (sc_wrapper.red_stem)
          en2 += sc_wrapper.red_stem(ii, jj, ii, jj - 1, &sc_wrapper);

        if (fij == en2) {
          *i          = *j = -1;
          *k          = ii;
          *l          = jj - 1;
          *component2 = 2;
          return 1;
        }
      }

      if (evaluate(ii, jj, ii + 1, jj - 1, VRNA_DECOMP_ML_STEM, &hc_dat_local)) {
        en2 = 2 * P->MLbase *
              n_seq;
        en2 += (sliding_window) ? c_local[ii + 1][jj - 1 - (ii + 1)] : my_c[idx[jj - 1] + ii + 1];

        switch (fc->type) {
          case VRNA_FC_TYPE_SINGLE:
            type =
              (sliding_window) ? vrna_get_ptype_window(ii + 1, jj - 1,
                                                       ptype_local) : vrna_get_ptype(
                idx[jj - 1] + ii + 1,
                ptype);
            en2 += E_MLstem(type, S1[ii], S1[jj], P);
            break;

          case VRNA_FC_TYPE_COMPARATIVE:
            for (s = 0; s < n_seq; s++) {
              type  = vrna_get_ptype_md(SS[s][ii], SS[s][jj], md);
              en2   += E_MLstem(type,
                                S5[s][ii],
                                S3[s][jj],
                                P);
            }
            break;
        }

        if (sc_wrapper.red_stem)
          en2 += sc_wrapper.red_stem(ii, jj, ii + 1, jj - 1, &sc_wrapper);

        if (fij == en2) {
          *i          = *j = -1;
          *k          = ii + 1;
          *l          = jj - 1;
          *component2 = 2;
          return 1;
        }
      }

      break;
  }

  /* 2. Test for possible split point */
  for (u = ii + 1 + turn; u <= jj - 2 - turn; u++) {
    if (evaluate(ii, jj, u, u + 1, VRNA_DECOMP_ML_ML_ML, &hc_dat_local)) {
      if (sliding_window)
        en = fML_local[ii][u - ii] +
             fML_local[u + 1][jj - (u + 1)];
      else
        en = my_fML[idx[u] + ii] +
             my_fML[idx[jj] + u + 1];

      if (sc_wrapper.decomp_ml)
        en += sc_wrapper.decomp_ml(ii, jj, u, u + 1, &sc_wrapper);

      if (fij == en) {
        *i  = ii;
        *j  = u;
        *k  = u + 1;
        *l  = jj;
        return 1;
      }
    }
  }

  /* 3. last chance! Maybe coax stack */
  if (dangle_model == 3) {
    int ik, k1j;
    k1j = (sliding_window) ? 0 : idx[jj] + ii + turn + 2;
    for (u = ii + 1 + turn; u <= jj - 2 - turn; u++, k1j++) {
      ik = (sliding_window) ? 0 : idx[u] + ii;
      if (evaluate(ii, u, u + 1, jj, VRNA_DECOMP_ML_COAXIAL_ENC, &hc_dat_local)) {
        en = 2 * P->MLintern[1] *
             n_seq;
        en += (sliding_window) ? c_local[ii][u - ii] + c_local[u + 1][jj - (u + 1)] : my_c[ik] +
              my_c[k1j];

        switch (fc->type) {
          case VRNA_FC_TYPE_SINGLE:
            type =
              (sliding_window) ? rtype[vrna_get_ptype_window(ii, u,
                                                             ptype_local)] : rtype[vrna_get_ptype(
                                                                                     ik,
                                                                                     ptype)
              ];
            type_2 =
              (sliding_window) ? rtype[vrna_get_ptype_window(u + 1, jj,
                                                             ptype_local)] : rtype[vrna_get_ptype(
                                                                                     k1j, ptype)];
            en += P->stack[type][type_2];
            break;

          case VRNA_FC_TYPE_COMPARATIVE:
            for (s = 0; s < n_seq; s++) {
              type    = vrna_get_ptype_md(SS[s][u], SS[s][ii], md);
              type_2  = vrna_get_ptype_md(SS[s][jj], SS[s][u + 1], md);
              en      += P->stack[type][type_2];
            }
            break;
        }

        if (sc_wrapper.coaxial_enc)
          en += sc_wrapper.coaxial_enc(ii, u, u + 1, jj, &sc_wrapper);

        if (fij == en) {
          *i          = ii;
          *j          = u;
          *k          = u + 1;
          *l          = jj;
          *component1 = *component2 = 2;
          return 1;
        }
      }
    }
  }

  return 0;
}


PRIVATE int
BT_mb_loop(vrna_fold_compound_t *fc,
           int                  *i,
           int                  *j,
           int                  *k,
           int                  en,
           int                  *component1,
           int                  *component2)
{
  unsigned char             sliding_window;
  char                      *ptype, **ptype_local;
  short                     s5, s3, *S1, **SS, **S5, **S3;
  unsigned int              *sn, *se, n_seq, s, *tt;
  int                       ij, p, q, r, e, tmp_en, *idx, turn, dangle_model,
                            *my_c, *my_fML, *my_fc, *rtype, type, type_2, **c_local, **fML_local;
  vrna_param_t              *P;
  vrna_md_t                 *md;
  vrna_sc_t                 *sc;
  vrna_callback_hc_evaluate *evaluate;
  struct hc_mb_def_dat      hc_dat_local;
  vrna_callback_hc_evaluate *evaluate_ext;
  struct hc_mb_def_dat      hc_dat_local_ext;
  struct sc_mb_dat          sc_wrapper;

  sliding_window  = (fc->hc->type == VRNA_HC_WINDOW) ? 1 : 0;
  n_seq           = (fc->type == VRNA_FC_TYPE_SINGLE) ? 1 : fc->n_seq;
  idx             = (sliding_window) ? NULL : fc->jindx;
  ij              = (sliding_window) ? 0 : idx[*j] + *i;
  S1              = (fc->type == VRNA_FC_TYPE_SINGLE) ? fc->sequence_encoding : NULL;
  SS              = (fc->type == VRNA_FC_TYPE_SINGLE) ? NULL : fc->S;
  S5              = (fc->type == VRNA_FC_TYPE_SINGLE) ? NULL : fc->S5;
  S3              = (fc->type == VRNA_FC_TYPE_SINGLE) ? NULL : fc->S3;
  P               = fc->params;
  md              = &(P->model_details);
  sn              = fc->strand_number;
  se              = fc->strand_end;
  sc              = fc->sc;
  my_c            = (sliding_window) ? NULL : fc->matrices->c;
  my_fML          = (sliding_window) ? NULL : fc->matrices->fML;
  c_local         = (sliding_window) ? fc->matrices->c_local : NULL;
  fML_local       = (sliding_window) ? fc->matrices->fML_local : NULL;
  my_fc           = (sliding_window) ? NULL : fc->matrices->fc;
  turn            = md->min_loop_size;
  ptype           = (fc->type == VRNA_FC_TYPE_SINGLE) ? (sliding_window ? NULL : fc->ptype) : NULL;
  ptype_local     = (sliding_window) ? fc->ptype_local : NULL;
  rtype           = &(md->rtype[0]);
  type            = (fc->type == VRNA_FC_TYPE_SINGLE) ?
                    (sliding_window ? rtype[vrna_get_ptype_window(*i, *j,
                                                                  ptype_local)] : rtype[
                       vrna_get_ptype(ij, ptype)]) :
                    0;
  tt            = NULL;
  dangle_model  = md->dangles;
  evaluate      = prepare_hc_mb_def(fc, &hc_dat_local);
  evaluate_ext  = prepare_hc_mb_def_ext(fc, &hc_dat_local_ext);

  init_sc_mb(fc, &sc_wrapper);

  p = *i + 1;
  q = *j - 1;

  r = q - turn - 1;

  /* is it a fake multi-loop? */
  if (sn[*i] != sn[*j]) {
    if (evaluate_ext(*i, *j, *i, *j, VRNA_DECOMP_EXT_STEM, &hc_dat_local)) {
      int ii, jj;
      ii  = jj = 0;
      e   = my_fc[p] + my_fc[q];
      if (sc)
        if (sc->energy_bp)
          e += sc->energy_bp[ij];

      s5  = (sn[q] == sn[*j]) ? S1[q] : -1;
      s3  = (sn[*i] == sn[p]) ? S1[p] : -1;

      switch (dangle_model) {
        case 0:
          if (en == e + vrna_E_ext_stem(type, -1, -1, P))
            ii = p, jj = q;

          break;

        case 2:
          if (en == e + vrna_E_ext_stem(type, s5, s3, P))
            ii = p, jj = q;

          break;

        default:
          if (en == e + vrna_E_ext_stem(type, -1, -1, P)) {
            ii = p, jj = q;
            break;
          }

          if (evaluate_ext(p, q, p + 1, q, VRNA_DECOMP_EXT_EXT, &hc_dat_local_ext)) {
            e = my_fc[p + 1] + my_fc[q];
            if (sc) {
              if (sc->energy_up)
                e += sc->energy_up[p][1];

              if (sc->energy_bp)
                e += sc->energy_bp[ij];
            }

            if (en == e + vrna_E_ext_stem(type, -1, s3, P)) {
              ii  = p + 1;
              jj  = q;
              break;
            }
          }

          if (evaluate_ext(p, q, p, q - 1, VRNA_DECOMP_EXT_EXT, &hc_dat_local_ext)) {
            e = my_fc[p] + my_fc[q - 1];
            if (sc) {
              if (sc->energy_up)
                e += sc->energy_up[q][1];

              if (sc->energy_bp)
                e += sc->energy_bp[ij];
            }

            if (en == e + vrna_E_ext_stem(type, s5, -1, P)) {
              ii  = p;
              jj  = q - 1;
              break;
            }
          }

          if (evaluate_ext(p, q, p + 1, q - 1, VRNA_DECOMP_EXT_EXT, &hc_dat_local_ext)) {
            e = my_fc[p + 1] + my_fc[q - 1];
            if (sc) {
              if (sc->energy_up)
                e += sc->energy_up[p][1] + sc->energy_up[q][1];

              if (sc->energy_bp)
                e += sc->energy_bp[ij];
            }

            if (en == e + vrna_E_ext_stem(type, s5, s3, P)) {
              ii  = p + 1;
              jj  = q - 1;
              break;
            }
          }

          break;
      }

      if (ii) {
        /* found a decomposition */
        *component1 = 3;
        *i          = ii;
        *k          = se[0];
        *j          = jj;
        *component2 = 4;
        return 1;
      }
    }
  }

  /* true multi loop? */
  *component1 = *component2 = 1;  /* both components are MB loop parts by default */

  s5  = -1;
  s3  = -1;

  if (fc->type == VRNA_FC_TYPE_COMPARATIVE) {
    tt = (unsigned int *)vrna_alloc(sizeof(unsigned int) * n_seq);
    for (s = 0; s < n_seq; s++)
      tt[s] = vrna_get_ptype_md(SS[s][*j], SS[s][*i], md);
  } else {
    if (sn[q] == sn[*j])
      s5 = S1[q];

    if (sn[*i] == sn[p])
      s3 = S1[p];
  }

  if (evaluate(*i, *j, p, q, VRNA_DECOMP_PAIR_ML, &hc_dat_local)) {
    e = en -
        P->MLclosing *
        n_seq;

    if (dangles == 2) {
      switch (fc->type) {
        case VRNA_FC_TYPE_SINGLE:
          e -= E_MLstem(type, s5, s3, P);
          break;

        case VRNA_FC_TYPE_COMPARATIVE:
          for (s = 0; s < n_seq; s++)
            e -= E_MLstem(tt[s], S5[s][*j], S3[s][*i], P);
          break;
      }
    } else {
      switch (fc->type) {
        case VRNA_FC_TYPE_SINGLE:
          e -= E_MLstem(type, -1, -1, P);
          break;

        case VRNA_FC_TYPE_COMPARATIVE:
          for (s = 0; s < n_seq; s++)
            e -= E_MLstem(tt[s], -1, -1, P);
          break;
      }
    }

    if (sc_wrapper.pair)
      e -= sc_wrapper.pair(*i, *j, &sc_wrapper);

    for (r = *i + 2 + turn; r < *j - 2 - turn; ++r) {
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
    if (evaluate(*i, *j, p + 1, q, VRNA_DECOMP_PAIR_ML, &hc_dat_local)) {
      e = en -
          (P->MLclosing + P->MLbase) *
          n_seq;

      if (sc_wrapper.pair5)
        e -= sc_wrapper.pair5(*i, *j, &sc_wrapper);

      switch (fc->type) {
        case VRNA_FC_TYPE_SINGLE:
          e -= E_MLstem(type, -1, s3, P);
          break;

        case VRNA_FC_TYPE_COMPARATIVE:
          for (s = 0; s < n_seq; s++)
            e -= E_MLstem(tt[s], -1, S3[s][*i], P);
          break;
      }

      for (r = p + turn + 1; r < q - turn - 1; ++r) {
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

    if (evaluate(*i, *j, p, q - 1, VRNA_DECOMP_PAIR_ML, &hc_dat_local)) {
      e = en -
          (P->MLclosing + P->MLbase) *
          n_seq;

      if (sc_wrapper.pair3)
        e -= sc_wrapper.pair3(*i, *j, &sc_wrapper);

      switch (fc->type) {
        case VRNA_FC_TYPE_SINGLE:
          e -= E_MLstem(type, s5, -1, P);
          break;

        case VRNA_FC_TYPE_COMPARATIVE:
          for (s = 0; s < n_seq; s++)
            e -= E_MLstem(tt[s], S5[s][*j], -1, P);
          break;
      }

      for (r = p + turn + 1; r < q - turn - 1; ++r) {
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

    if (evaluate(*i, *j, p + 1, q - 1, VRNA_DECOMP_PAIR_ML, &hc_dat_local)) {
      e = en -
          (P->MLclosing + 2 * P->MLbase) *
          n_seq;

      if (sc_wrapper.pair53)
        e -= sc_wrapper.pair53(*i, *j, &sc_wrapper);

      switch (fc->type) {
        case VRNA_FC_TYPE_SINGLE:
          e -= E_MLstem(type, s5, s3, P);
          break;

        case VRNA_FC_TYPE_COMPARATIVE:
          for (s = 0; s < n_seq; s++)
            e -= E_MLstem(tt[s], S5[s][*j], S3[s][*i], P);
          break;
      }

      for (r = p + turn + 1; r < q - turn - 1; ++r) {
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
        e -= sc_wrapper.pair(*i, *j, &sc_wrapper);

      if (fc->type == VRNA_FC_TYPE_SINGLE)
        type = rtype[type];

      for (r = p + turn + 1; r < q - turn - 1; ++r) {
        if (evaluate(*i, *j, p, r, VRNA_DECOMP_ML_COAXIAL, &hc_dat_local)) {
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
            tmp_en += sc_wrapper.coaxial_cls(*i, *j, p, r, &sc_wrapper);

          if (e == tmp_en) {
            *component1 = 2;
            goto odd_dangles_exit;
          }
        }

        if (evaluate(*i, *j, r + 1, q, VRNA_DECOMP_ML_COAXIAL, &hc_dat_local)) {
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
            tmp_en += sc_wrapper.coaxial_cls(*i, *j, r + 1, q, &sc_wrapper);

          if (e == tmp_en) {
            *component2 = 2;
            goto odd_dangles_exit;
          }
        }
      }
    }

odd_dangles_exit:

    free(tt);
  }

  if (r <= *j - turn - 3) {
    *i  = p;
    *k  = r;
    *j  = q;
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
