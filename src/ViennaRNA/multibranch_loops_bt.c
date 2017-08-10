#ifdef HAVE_CONFIG_H
#include "config.h"
#endif


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <string.h>
#include "ViennaRNA/utils.h"
#include "ViennaRNA/fold_vars.h"
#include "ViennaRNA/energy_par.h"
#include "ViennaRNA/constraints.h"
#include "ViennaRNA/exterior_loops.h"
#include "ViennaRNA/gquad.h"
#include "ViennaRNA/structured_domains.h"
#include "ViennaRNA/unstructured_domains.h"
#include "ViennaRNA/multibranch_loops.h"

#ifdef __GNUC__
# define INLINE inline
#else
# define INLINE
#endif

#include "multibranch_loops.inc"

/*
 #################################
 # PRIVATE FUNCTION DECLARATIONS #
 #################################
 */
PRIVATE int
BT_mb_loop(vrna_fold_compound_t *vc,
           int                  *i,
           int                  *j,
           int                  *k,
           int                  en,
           int                  *component1,
           int                  *component2);


PRIVATE int
BT_mb_loop_comparative(vrna_fold_compound_t *vc,
                       int                  *i,
                       int                  *j,
                       int                  *k,
                       int                  en,
                       int                  *component1,
                       int                  *component2);


PRIVATE int
BT_mb_loop_window(vrna_fold_compound_t  *vc,
                  int                   *i,
                  int                   *j,
                  int                   *k,
                  int                   en,
                  int                   *component1,
                  int                   *component2);


PRIVATE int
BT_mb_loop_window_comparative(vrna_fold_compound_t  *vc,
                              int                   *i,
                              int                   *j,
                              int                   *k,
                              int                   en,
                              int                   *component1,
                              int                   *component2);


PRIVATE int
BT_mb_loop_split(vrna_fold_compound_t *vc,
                 int                  *i,
                 int                  *j,
                 int                  *k,
                 int                  *l,
                 int                  *component1,
                 int                  *component2,
                 vrna_bp_stack_t      *bp_stack,
                 int                  *stack_count);


PRIVATE int
BT_mb_loop_split_comparative(vrna_fold_compound_t *vc,
                             int                  *i,
                             int                  *j,
                             int                  *k,
                             int                  *l,
                             int                  *component1,
                             int                  *component2,
                             vrna_bp_stack_t      *bp_stack,
                             int                  *stack_count);


PRIVATE int
BT_mb_loop_split_window(vrna_fold_compound_t  *vc,
                        int                   *i,
                        int                   *j,
                        int                   *k,
                        int                   *l,
                        int                   *component1,
                        int                   *component2,
                        vrna_bp_stack_t       *bp_stack,
                        int                   *stack_count);


PRIVATE int
BT_mb_loop_split_window_comparative(vrna_fold_compound_t  *vc,
                                    int                   *i,
                                    int                   *j,
                                    int                   *k,
                                    int                   *l,
                                    int                   *component1,
                                    int                   *component2,
                                    vrna_bp_stack_t       *bp_stack,
                                    int                   *stack_count);


/*
 #################################
 # BEGIN OF FUNCTION DEFINITIONS #
 #################################
 */
PUBLIC int
vrna_BT_mb_loop_split(vrna_fold_compound_t  *vc,
                      int                   *i,
                      int                   *j,
                      int                   *k,
                      int                   *l,
                      int                   *c1,
                      int                   *c2,
                      vrna_bp_stack_t       *bp_stack,
                      int                   *stack_count)
{
  if (vc) {
    switch (vc->type) {
      case VRNA_FC_TYPE_SINGLE:
        if (vc->hc->type == VRNA_HC_WINDOW)
          return BT_mb_loop_split_window(vc, i, j, k, l, c1, c2, bp_stack, stack_count);
        else
          return BT_mb_loop_split(vc, i, j, k, l, c1, c2, bp_stack, stack_count);

        break;

      case VRNA_FC_TYPE_COMPARATIVE:
        if (vc->hc->type == VRNA_HC_WINDOW)
          return BT_mb_loop_split_window_comparative(vc, i, j, k, l, c1, c2, bp_stack, stack_count);

        else
          return BT_mb_loop_split_comparative(vc, i, j, k, l, c1, c2, bp_stack, stack_count);

        break;
    }
  }

  return 0;
}


PUBLIC int
vrna_BT_mb_loop(vrna_fold_compound_t  *vc,
                int                   *i,
                int                   *j,
                int                   *k,
                int                   en,
                int                   *c1,
                int                   *c2)
{
  if (vc) {
    switch (vc->type) {
      case VRNA_FC_TYPE_SINGLE:
        if (vc->hc->type == VRNA_HC_WINDOW)
          return BT_mb_loop_window(vc, i, j, k, en, c1, c2);
        else
          return BT_mb_loop(vc, i, j, k, en, c1, c2);

        break;

      case VRNA_FC_TYPE_COMPARATIVE:
        if (vc->hc->type == VRNA_HC_WINDOW)
          return BT_mb_loop_window_comparative(vc, i, j, k, en, c1, c2);
        else
          return BT_mb_loop_comparative(vc, i, j, k, en, c1, c2);

        break;
    }
  }

  return 0;
}


PUBLIC int
vrna_BT_mb_loop_fake(vrna_fold_compound_t *vc,
                     int                  *u,
                     int                  *i,
                     int                  *j,
                     vrna_bp_stack_t      *bp_stack,
                     int                  *stack_count)
{
  char                      *ptype;
  short                     mm5, mm3, *S1;
  unsigned int              strands, *sn, *so, *ss, *se;
  int                       length, ii, jj, k, en, fij, fi, *my_c, *my_fc, *my_ggg,
                            *idx, with_gquad, dangle_model, turn, type;
  vrna_param_t              *P;
  vrna_md_t                 *md;
  vrna_hc_t                 *hc;
  vrna_sc_t                 *sc;
  vrna_callback_hc_evaluate *evaluate;
  struct default_data       hc_dat_local;

  length        = vc->length;
  P             = vc->params;
  md            = &(P->model_details);
  strands       = vc->strands;
  sn            = vc->strand_number;
  so            = vc->strand_order;
  ss            = vc->strand_start;
  se            = vc->strand_end;
  hc            = vc->hc;
  sc            = vc->sc;
  S1            = vc->sequence_encoding;
  ptype         = vc->ptype;
  idx           = vc->jindx;
  my_c          = vc->matrices->c;
  my_fc         = vc->matrices->fc;
  my_ggg        = vc->matrices->ggg;
  turn          = md->min_loop_size;
  with_gquad    = md->gquad;
  dangle_model  = md->dangles;

  hc_dat_local.idx    = vc->jindx;
  hc_dat_local.mx     = hc->matrix;
  hc_dat_local.hc_up  = hc->up_ext;
  hc_dat_local.sn     = sn;

  if (hc->f) {
    evaluate            = &hc_default_user_ext;
    hc_dat_local.hc_f   = hc->f;
    hc_dat_local.hc_dat = hc->data;
  } else {
    evaluate = &hc_default_ext;
  }

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
            type = get_pair_type(idx[k] + ii, ptype);

            if (fij == my_fc[k + 1] + my_c[idx[k] + ii] + E_ExtLoop(type, -1, -1, P)) {
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
              return vrna_BT_gquad_mfe(vc, ii, k, bp_stack, stack_count);
            }
          }
        }
        break;

      case 2:
        for (k = ii + turn + 1; k <= jj; k++) {
          if (evaluate(ii, jj, k, k + 1, VRNA_DECOMP_EXT_STEM_EXT, &hc_dat_local)) {
            mm3   = (sn[k] == sn[k + 1]) ? S1[k + 1] : -1;
            type  = get_pair_type(idx[k] + ii, ptype);

            if (fij == my_fc[k + 1] + my_c[idx[k] + ii] + E_ExtLoop(type, mm5, mm3, P)) {
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
              return vrna_BT_gquad_mfe(vc, ii, k, bp_stack, stack_count);
            }
          }
        }
        break;

      default:
        for (k = ii + turn + 1; k <= jj; k++) {
          if (evaluate(ii, jj, k, k + 1, VRNA_DECOMP_EXT_STEM_EXT, &hc_dat_local)) {
            type = get_pair_type(idx[k] + ii, ptype);

            if (fij == my_fc[k + 1] + my_c[idx[k] + ii] + E_ExtLoop(type, -1, -1, P)) {
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

            if (fij == my_fc[k + 2] + en + E_ExtLoop(type, -1, mm3, P)) {
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
              return vrna_BT_gquad_mfe(vc, ii, k, bp_stack, stack_count);
            }
          }

          if (evaluate(ii, jj, k, k + 1, VRNA_DECOMP_EXT_STEM_EXT1, &hc_dat_local)) {
            mm5   = (sn[ii] == sn[ii + 1]) ? S1[ii] : -1;
            mm3   = (sn[k] == sn[k + 1]) ? S1[k + 1] : -1;
            type  = get_pair_type(idx[k] + ii + 1, ptype);

            en = my_c[idx[k] + ii + 1];
            if (sc)
              if (sc->energy_up)
                en += sc->energy_up[ii][1];

            if (fij == en + my_fc[k + 1] + E_ExtLoop(type, mm5, -1, P)) {
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
            type  = get_pair_type(idx[k] + ii + 1, ptype);

            en = my_c[idx[k] + ii + 1];
            if (sc)
              if (sc->energy_up)
                en += sc->energy_up[k + 1][1];

            if (fij == en + my_fc[k + 2] + E_ExtLoop(type, mm5, mm3, P)) {
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
              return vrna_BT_gquad_mfe(vc, k, jj, bp_stack, stack_count);
            }
          }

          if (evaluate(ii, jj, k - 1, k, VRNA_DECOMP_EXT_EXT_STEM, &hc_dat_local)) {
            type = get_pair_type(idx[jj] + k, ptype);

            en = my_c[idx[jj] + k];
            if (sn[k] != sn[jj])
              en += P->DuplexInit;

            if (fij == my_fc[k - 1] + en + E_ExtLoop(type, -1, -1, P)) {
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
              return vrna_BT_gquad_mfe(vc, k, jj, bp_stack, stack_count);
            }
          }

          if (evaluate(ii, jj, k - 1, k, VRNA_DECOMP_EXT_EXT_STEM, &hc_dat_local)) {
            mm5   = ((k > 1) && (sn[k - 1] == sn[k])) ? S1[k - 1] : -1;
            type  = get_pair_type(idx[jj] + k, ptype);

            en = my_c[idx[jj] + k];
            if (sn[k] != sn[jj])
              en += P->DuplexInit;

            if (fij == my_fc[k - 1] + en + E_ExtLoop(type, mm5, mm3, P)) {
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
              return vrna_BT_gquad_mfe(vc, k, jj, bp_stack, stack_count);
            }
          }

          if (evaluate(ii, jj, k - 1, k, VRNA_DECOMP_EXT_EXT_STEM, &hc_dat_local)) {
            type = get_pair_type(idx[jj] + k, ptype);

            en = my_c[idx[jj] + k];
            if (sn[k] != sn[jj])
              en += P->DuplexInit;

            if (fij == my_fc[k - 1] + en + E_ExtLoop(type, -1, -1, P)) {
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
            type = get_pair_type(idx[jj] + k, ptype);

            en = my_c[idx[jj] + k];
            if (sn[k] != sn[jj])
              en += P->DuplexInit;

            mm5 = S1[k - 1];

            if (sc)
              if (sc->energy_up)
                en += sc->energy_up[k - 1][1];

            if (fij == my_fc[k - 2] + en + E_ExtLoop(type, mm5, -1, P)) {
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
            type = get_pair_type(idx[jj - 1] + k, ptype);

            mm3 = S1[jj];
            en  = my_c[idx[jj - 1] + k];
            if (sn[k] != sn[jj - 1])
              en += P->DuplexInit;         /* ??? */

            if (sc)
              if (sc->energy_up)
                en += sc->energy_up[jj][1];

            if (fij == en + my_fc[k - 1] + E_ExtLoop(type, -1, mm3, P)) {
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
            type = get_pair_type(idx[jj - 1] + k, ptype);

            mm3 = S1[jj];
            en  = my_c[idx[jj - 1] + k];
            if (sn[k] != sn[jj - 1])
              en += P->DuplexInit;         /* ??? */

            mm5 = (sn[k - 1] == sn[k]) ? S1[k - 1] : -1;
            if (sc)
              if (sc->energy_up)
                en += sc->energy_up[k - 1][1];

            if (fij == my_fc[k - 2] + en + E_ExtLoop(type, mm5, mm3, P)) {
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
BT_mb_loop_split(vrna_fold_compound_t *vc,
                 int                  *i,
                 int                  *j,
                 int                  *k,
                 int                  *l,
                 int                  *component1,
                 int                  *component2,
                 vrna_bp_stack_t      *bp_stack,
                 int                  *stack_count)
{
  char                      *ptype;
  short                     *S1;
  int                       ij, ii, jj, fij, fi, u, en, *my_c, *my_fML, *my_ggg,
                            turn, *idx, with_gquad, dangle_model, *rtype, kk, cnt,
                            with_ud, type, type_2;
  vrna_param_t              *P;
  vrna_md_t                 *md;
  vrna_hc_t                 *hc;
  vrna_sc_t                 *sc;
  vrna_ud_t                 *domains_up;
  vrna_callback_hc_evaluate *evaluate;
  struct default_data       hc_dat_local;

  P           = vc->params;
  md          = &(P->model_details);
  hc          = vc->hc;
  sc          = vc->sc;
  idx         = vc->jindx;
  ptype       = vc->ptype;
  rtype       = &(md->rtype[0]);
  S1          = vc->sequence_encoding;
  domains_up  = vc->domains_up;

  my_c          = vc->matrices->c;
  my_fML        = vc->matrices->fML;
  my_ggg        = vc->matrices->ggg;
  turn          = md->min_loop_size;
  with_gquad    = md->gquad;
  with_ud       = (domains_up && domains_up->energy_cb) ? 1 : 0;
  dangle_model  = md->dangles;

  hc_dat_local.idx    = vc->jindx;
  hc_dat_local.mx     = hc->matrix;
  hc_dat_local.hc_up  = hc->up_ml;
  hc_dat_local.sn     = vc->strand_number;

  if (hc->f) {
    evaluate            = &hc_default_user;
    hc_dat_local.hc_f   = hc->f;
    hc_dat_local.hc_dat = hc->data;
  } else {
    evaluate = &hc_default;
  }

  ii  = *i;
  jj  = *j;

  if (with_ud) {
    /* nibble off unpaired stretches at 3' site */
    do {
      fij = my_fML[idx[jj] + ii];
      fi  = INF;

      /* process regular unpaired nucleotides (unbound by ligand) first */
      if (evaluate(ii, jj, ii, jj - 1, VRNA_DECOMP_ML_ML, &hc_dat_local)) {
        fi = my_fML[idx[jj - 1] + ii] +
             P->MLbase;

        if (sc) {
          if (sc->energy_up)
            fi += sc->energy_up[jj][1];

          if (sc->f)
            fi += sc->f(ii, jj, ii, jj - 1, VRNA_DECOMP_ML_ML, sc->data);
        }

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
          en = domains_up->energy_cb(vc,
                                     kk,
                                     jj,
                                     VRNA_UNSTRUCTURED_DOMAIN_MB_LOOP | VRNA_UNSTRUCTURED_DOMAIN_MOTIF,
                                     domains_up->data);

          if (sc) {
            if (sc->energy_up)
              en += sc->energy_up[kk][u];

            if (sc->f)
              en += sc->f(ii, jj, ii, jj - u, VRNA_DECOMP_ML_ML, sc->data);
          }

          fi = my_fML[idx[kk - 1] + ii] +
               u * P->MLbase;

          fi += en;

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
      fij = my_fML[idx[jj] + ii];
      fi  = INF;

      /* again, process regular unpaired nucleotides (unbound by ligand) first */
      if (evaluate(ii, jj, ii + 1, jj, VRNA_DECOMP_ML_ML, &hc_dat_local)) {
        fi = my_fML[idx[jj] + ii + 1] +
             P->MLbase;

        if (sc) {
          if (sc->energy_up)
            fi += sc->energy_up[ii][1];

          if (sc->f)
            fi += sc->f(ii, jj, ii + 1, jj, VRNA_DECOMP_ML_ML, sc->data);
        }

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
          en = domains_up->energy_cb(vc,
                                     ii,
                                     kk,
                                     VRNA_UNSTRUCTURED_DOMAIN_MB_LOOP | VRNA_UNSTRUCTURED_DOMAIN_MOTIF,
                                     domains_up->data);

          if (sc) {
            if (sc->energy_up)
              en += sc->energy_up[ii][u];

            if (sc->f)
              en += sc->f(ii, jj, ii + u, jj, VRNA_DECOMP_ML_ML, sc->data);
          }

          fi = my_fML[idx[jj] + kk + 1] +
               u * P->MLbase;
          fi += en;

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
      fij = my_fML[idx[jj] + ii];
      fi  = INF;

      if (evaluate(ii, jj, ii, jj - 1, VRNA_DECOMP_ML_ML, &hc_dat_local)) {
        fi = my_fML[idx[jj - 1] + ii] +
             P->MLbase;

        if (sc) {
          if (sc->energy_up)
            fi += sc->energy_up[jj][1];

          if (sc->f)
            fi += sc->f(ii, jj, ii, jj - 1, VRNA_DECOMP_ML_ML, sc->data);
        }
      }

      if (--jj == 0)
        break;
    } while (fij == fi);
    jj++;

    /* nibble off unpaired 5' bases */
    do {
      fij = my_fML[idx[jj] + ii];
      fi  = INF;

      if (evaluate(ii, jj, ii + 1, jj, VRNA_DECOMP_ML_ML, &hc_dat_local)) {
        fi = my_fML[idx[jj] + ii + 1] +
             P->MLbase;

        if (sc) {
          if (sc->energy_up)
            fi += sc->energy_up[ii][1];

          if (sc->f)
            fi += sc->f(ii, jj, ii + 1, jj, VRNA_DECOMP_ML_ML, sc->data);
        }
      }

      if (++ii == jj)
        break;
    } while (fij == fi);
    ii--;

    if (jj < ii + turn + 1) /* no more pairs */
      return 0;
  }

  ij = idx[jj] + ii;

  *component1 = *component2 = 1; /* split into two multi loop parts by default */

  /* 1. test for single component */

  if (with_gquad) {
    if (fij == my_ggg[ij] + E_MLstem(0, -1, -1, P)) {
      *i  = *j = -1;
      *k  = *l = -1;
      return vrna_BT_gquad_mfe(vc, ii, jj, bp_stack, stack_count);
    }
  }

  type  = get_pair_type(ij, ptype);
  en    = my_c[ij];

  if (sc)
    if (sc->f)
      en += sc->f(ii, jj, ii, jj, VRNA_DECOMP_ML_STEM, sc->data);

  switch (dangle_model) {
    case 0:
      if (evaluate(ii, jj, ii, jj, VRNA_DECOMP_ML_STEM, &hc_dat_local)) {
        if (fij == en + E_MLstem(type, -1, -1, P)) {
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
        if (fij == en + E_MLstem(type, S1[ii - 1], S1[jj + 1], P)) {
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
        if (fij == en + E_MLstem(type, -1, -1, P)) {
          *i          = *j = -1;
          *k          = ii;
          *l          = jj;
          *component2 = 2;
          return 1;
        }
      }

      if (evaluate(ii, jj, ii + 1, jj, VRNA_DECOMP_ML_STEM, &hc_dat_local)) {
        int tmp_en = fij;
        if (sc) {
          if (sc->energy_up)
            tmp_en -= sc->energy_up[ii][1];

          if (sc->f)
            tmp_en -= sc->f(ii, jj, ii + 1, jj, VRNA_DECOMP_ML_STEM, sc->data);
        }

        type = get_pair_type(ij + 1, ptype);

        if (tmp_en == my_c[ij + 1] + E_MLstem(type, S1[ii], -1, P) + P->MLbase) {
          *i          = *j = -1;
          *k          = ii + 1;
          *l          = jj;
          *component2 = 2;
          return 1;
        }
      }

      if (evaluate(ii, jj, ii, jj - 1, VRNA_DECOMP_ML_STEM, &hc_dat_local)) {
        int tmp_en = fij;
        if (sc) {
          if (sc->energy_up)
            tmp_en -= sc->energy_up[jj][1];

          if (sc->f)
            tmp_en -= sc->f(ii, jj, ii, jj - 1, VRNA_DECOMP_ML_STEM, sc->data);
        }

        type = get_pair_type(idx[jj - 1] + ii, ptype);

        if (tmp_en == my_c[idx[jj - 1] + ii] + E_MLstem(type, -1, S1[jj], P) + P->MLbase) {
          *i          = *j = -1;
          *k          = ii;
          *l          = jj - 1;
          *component2 = 2;
          return 1;
        }
      }

      if (evaluate(ii, jj, ii + 1, jj - 1, VRNA_DECOMP_ML_STEM, &hc_dat_local)) {
        int tmp_en = fij;
        if (sc) {
          if (sc->energy_up)
            tmp_en -= sc->energy_up[ii][1] + sc->energy_up[jj][1];

          if (sc->f)
            tmp_en -= sc->f(ii, jj, ii + 1, jj - 1, VRNA_DECOMP_ML_STEM, sc->data);
        }

        type = get_pair_type(idx[jj - 1] + ii + 1, ptype);

        if (tmp_en ==
            my_c[idx[jj - 1] + ii + 1] + E_MLstem(type, S1[ii], S1[jj], P) + 2 * P->MLbase) {
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
  if (hc->f) {
    for (u = ii + 1 + turn; u <= jj - 2 - turn; u++) {
      if (hc->f(ii, jj, u, u + 1, VRNA_DECOMP_ML_ML_ML, hc->data)) {
        en = my_fML[idx[u] + ii] + my_fML[idx[jj] + u + 1];
        if (sc)
          if (sc->f)
            en += sc->f(ii, jj, u, u + 1, VRNA_DECOMP_ML_ML_ML, sc->data);

        if (fij == en) {
          *i  = ii;
          *j  = u;
          *k  = u + 1;
          *l  = jj;
          return 1;
        }
      }
    }
  } else {
    for (u = ii + 1 + turn; u <= jj - 2 - turn; u++) {
      en = my_fML[idx[u] + ii] + my_fML[idx[jj] + u + 1];
      if (sc)
        if (sc->f)
          en += sc->f(ii, jj, u, u + 1, VRNA_DECOMP_ML_ML_ML, sc->data);

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
    int ik, k1j, tmp_en;
    for (k1j = idx[jj] + ii + turn + 2, u = ii + 1 + turn; u <= jj - 2 - turn; u++, k1j++) {
      ik = idx[u] + ii;
      if (evaluate(ii, u, u + 1, jj, VRNA_DECOMP_ML_COAXIAL_ENC, &hc_dat_local)) {
        type    = rtype[get_pair_type(ik, ptype)];
        type_2  = rtype[get_pair_type(k1j, ptype)];

        tmp_en = my_c[ik] +
                 my_c[k1j] +
                 P->stack[type][type_2] +
                 2 * P->MLintern[1];

        if (sc)
          if (sc->f)
            tmp_en += sc->f(ii, u, u + 1, jj, VRNA_DECOMP_ML_COAXIAL, sc->data);

        if (fij == tmp_en) {
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
BT_mb_loop_split_comparative(vrna_fold_compound_t *vc,
                             int                  *i,
                             int                  *j,
                             int                  *k,
                             int                  *l,
                             int                  *component1,
                             int                  *component2,
                             vrna_bp_stack_t      *bp_stack,
                             int                  *stack_count)
{
  unsigned int              **a2s;
  short                     **S, **S5, **S3;
  int                       ij, ii, jj, fij, fi, u, en, *my_c, *my_fML, *my_ggg,
                            turn, *idx, with_gquad, dangle_model, ss, n_seq, tt;
  vrna_param_t              *P;
  vrna_md_t                 *md;
  vrna_hc_t                 *hc;
  vrna_sc_t                 **scs;
  vrna_callback_hc_evaluate *evaluate;
  struct default_data       hc_dat_local;

  n_seq = vc->n_seq;
  S     = vc->S;
  S5    = vc->S5;
  S3    = vc->S3;
  a2s   = vc->a2s;
  P     = vc->params;
  md    = &(P->model_details);
  hc    = vc->hc;
  scs   = vc->scs;
  idx   = vc->jindx;

  my_c          = vc->matrices->c;
  my_fML        = vc->matrices->fML;
  my_ggg        = vc->matrices->ggg;
  turn          = md->min_loop_size;
  with_gquad    = md->gquad;
  dangle_model  = md->dangles;

  hc_dat_local.idx    = vc->jindx;
  hc_dat_local.mx     = hc->matrix;
  hc_dat_local.hc_up  = hc->up_ml;
  hc_dat_local.sn     = vc->strand_number;

  if (hc->f) {
    evaluate            = &hc_default_user;
    hc_dat_local.hc_f   = hc->f;
    hc_dat_local.hc_dat = hc->data;
  } else {
    evaluate = &hc_default;
  }

  ii  = *i;
  jj  = *j;

  /* nibble off unpaired 3' bases */
  do {
    fij = my_fML[idx[jj] + ii];
    fi  = INF;

    if (evaluate(ii, jj, ii, jj - 1, VRNA_DECOMP_ML_ML, &hc_dat_local)) {
      fi = my_fML[idx[jj - 1] + ii] +
           n_seq * P->MLbase;

      if (scs) {
        for (ss = 0; ss < n_seq; ss++) {
          if (scs[ss]) {
            if (scs[ss]->energy_up)
              fi += scs[ss]->energy_up[a2s[ss][jj]][1];

            if (scs[ss]->f)
              fi += scs[ss]->f(ii, jj, ii, jj - 1, VRNA_DECOMP_ML_ML, scs[ss]->data);
          }
        }
      }
    }

    if (--jj == 0)
      break;
  } while (fij == fi);
  jj++;

  /* nibble off unpaired 5' bases */
  do {
    fij = my_fML[idx[jj] + ii];
    fi  = INF;

    if (evaluate(ii, jj, ii + 1, jj, VRNA_DECOMP_ML_ML, &hc_dat_local)) {
      fi = my_fML[idx[jj] + ii + 1] +
           n_seq * P->MLbase;

      if (scs) {
        for (ss = 0; ss < n_seq; ss++) {
          if (scs[ss]) {
            if (scs[ss]->energy_up)
              fi += scs[ss]->energy_up[a2s[ss][ii]][1];

            if (scs[ss]->f)
              fi += scs[ss]->f(ii, jj, ii + 1, jj, VRNA_DECOMP_ML_ML, scs[ss]->data);
          }
        }
      }
    }

    if (++ii == jj)
      break;
  } while (fij == fi);
  ii--;

  if (jj < ii + turn + 1) /* no more pairs */
    return 0;

  ij = idx[jj] + ii;

  *component1 = *component2 = 1; /* split into two multi loop parts by default */

  /* 1. test for single component */

  if (with_gquad) {
    if (fij == my_ggg[ij] + n_seq * E_MLstem(0, -1, -1, P)) {
      *i  = *j = -1;
      *k  = *l = -1;
      return vrna_BT_gquad_mfe(vc, ii, jj, bp_stack, stack_count);
    }
  }

  en = my_c[ij];

  if (scs) {
    for (ss = 0; ss < n_seq; ss++) {
      if (scs[ss])
        if (scs[ss]->f)
          en += scs[ss]->f(ii, jj, ii, jj, VRNA_DECOMP_ML_STEM, scs[ss]->data);
    }
  }

  switch (dangle_model) {
    case 0:
      if (evaluate(ii, jj, ii, jj, VRNA_DECOMP_ML_STEM, &hc_dat_local)) {
        for (ss = 0; ss < n_seq; ss++) {
          tt  = get_pair_type_md(S[ss][ii], S[ss][jj], md);
          en  += E_MLstem(tt, -1, -1, P);
        }
      }

      break;

    case 2:
      if (evaluate(ii, jj, ii, jj, VRNA_DECOMP_ML_STEM, &hc_dat_local)) {
        for (ss = 0; ss < n_seq; ss++) {
          tt  = get_pair_type_md(S[ss][ii], S[ss][jj], md);
          en  += E_MLstem(tt, S5[ss][ii], S3[ss][jj], P);
        }
      }

      break;
  }

  if (fij == en) {
    *i          = *j = -1;
    *k          = ii;
    *l          = jj;
    *component2 = 2;          /* 2nd part is structure enclosed by base pair */
    return 1;
  }

  /* 2. Test for possible split point */
  if (hc->f) {
    for (u = ii + 1 + turn; u <= jj - 2 - turn; u++) {
      if (hc->f(ii, jj, u, u + 1, VRNA_DECOMP_ML_ML_ML, hc->data)) {
        en = my_fML[idx[u] + ii] +
             my_fML[idx[jj] + u + 1];

        if (scs) {
          for (ss = 0; ss < n_seq; ss++)
            if (scs[ss] && scs[ss]->f)
              en += scs[ss]->f(ii, jj, u, u + 1, VRNA_DECOMP_ML_ML_ML, scs[ss]->data);
        }

        if (fij == en) {
          *i  = ii;
          *j  = u;
          *k  = u + 1;
          *l  = jj;
          return 1;
        }
      }
    }
  } else {
    for (u = ii + 1 + turn; u <= jj - 2 - turn; u++) {
      en = my_fML[idx[u] + ii] +
           my_fML[idx[jj] + u + 1];

      if (scs) {
        for (ss = 0; ss < n_seq; ss++)
          if (scs[ss] && scs[ss]->f)
            en += scs[ss]->f(ii, jj, u, u + 1, VRNA_DECOMP_ML_ML_ML, scs[ss]->data);
      }

      if (fij == en) {
        *i  = ii;
        *j  = u;
        *k  = u + 1;
        *l  = jj;
        return 1;
      }
    }
  }

  return 0;
}


PRIVATE int
BT_mb_loop_split_window(vrna_fold_compound_t  *vc,
                        int                   *i,
                        int                   *j,
                        int                   *k,
                        int                   *l,
                        int                   *component1,
                        int                   *component2,
                        vrna_bp_stack_t       *bp_stack,
                        int                   *stack_count)
{
  char                      **ptype;
  short                     *S1;
  int                       ii, jj, fij, fi, u, en, **c, **fML, **ggg, turn,
                            with_gquad, dangle_model, *rtype, type, type_2;
  vrna_param_t              *P;
  vrna_md_t                 *md;
  vrna_hc_t                 *hc;
  vrna_sc_t                 *sc;
  vrna_callback_hc_evaluate *evaluate;
  struct default_data       hc_dat_local;

  P     = vc->params;
  md    = &(P->model_details);
  hc    = vc->hc;
  sc    = vc->sc;
  ptype = vc->ptype_local;
  rtype = &(md->rtype[0]);
  S1    = vc->sequence_encoding;

  c             = vc->matrices->c_local;
  fML           = vc->matrices->fML_local;
  ggg           = vc->matrices->ggg_local;
  turn          = md->min_loop_size;
  with_gquad    = md->gquad;
  dangle_model  = md->dangles;

  hc_dat_local.idx        = vc->jindx;
  hc_dat_local.mx         = hc->matrix;
  hc_dat_local.mx_window  = hc->matrix_local;
  hc_dat_local.hc_up      = hc->up_ml;
  hc_dat_local.sn         = vc->strand_number;

  if (hc->f) {
    evaluate            = &hc_default_user_window;
    hc_dat_local.hc_f   = hc->f;
    hc_dat_local.hc_dat = hc->data;
  } else {
    evaluate = &hc_default_window;
  }

  ii  = *i;
  jj  = *j;

  /* nibble off unpaired 3' bases */
  do {
    fij = fML[ii][jj - ii];
    fi  = INF;

    if (evaluate(ii, jj, ii, jj - 1, VRNA_DECOMP_ML_ML, &hc_dat_local)) {
      fi = fML[ii][jj - 1 - ii] +
           P->MLbase;

      if (sc) {
        if (sc->energy_up)
          fi += sc->energy_up[jj][1];

        if (sc->f)
          fi += sc->f(ii, jj, ii, jj - 1, VRNA_DECOMP_ML_ML, sc->data);
      }
    }

    if (--jj == 0)
      break;
  } while (fij == fi);
  jj++;

  /* nibble off unpaired 5' bases */
  do {
    fij = fML[ii][jj - ii];
    fi  = INF;

    if (evaluate(ii, jj, ii + 1, jj, VRNA_DECOMP_ML_ML, &hc_dat_local)) {
      fi = fML[ii + 1][jj - (ii + 1)] +
           P->MLbase;

      if (sc) {
        if (sc->energy_up)
          fi += sc->energy_up[ii][1];

        if (sc->f)
          fi += sc->f(ii, jj, ii + 1, jj, VRNA_DECOMP_ML_ML, sc->data);
      }
    }

    if (++ii == jj)
      break;
  } while (fij == fi);
  ii--;

  if (jj < ii + turn + 1) /* no more pairs */
    return 0;

  *component1 = *component2 = 1; /* split into two multi loop parts by default */

  /* 1. test for single component */

  if (with_gquad) {
    if (fij == ggg[ii][jj - ii] + E_MLstem(0, -1, -1, P)) {
      *i  = *j = -1;
      *k  = *l = -1;
      return vrna_BT_gquad_mfe(vc, ii, jj, bp_stack, stack_count);
    }
  }

  type  = get_pair_type_window(ii, jj, ptype);
  en    = c[ii][jj - ii];

  if (sc)
    if (sc->f)
      en += sc->f(ii, jj, ii, jj, VRNA_DECOMP_ML_STEM, sc->data);

  switch (dangle_model) {
    case 0:
      if (evaluate(ii, jj, ii, jj, VRNA_DECOMP_ML_STEM, &hc_dat_local)) {
        if (fij == en + E_MLstem(type, -1, -1, P)) {
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
        if (fij == en + E_MLstem(type, S1[ii - 1], S1[jj + 1], P)) {
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
        if (fij == en + E_MLstem(type, -1, -1, P)) {
          *i          = *j = -1;
          *k          = ii;
          *l          = jj;
          *component2 = 2;
          return 1;
        }
      }

      if (evaluate(ii, jj, ii + 1, jj, VRNA_DECOMP_ML_STEM, &hc_dat_local)) {
        int tmp_en = fij;
        if (sc) {
          if (sc->energy_up)
            tmp_en -= sc->energy_up[ii][1];

          if (sc->f)
            tmp_en -= sc->f(ii, jj, ii + 1, jj, VRNA_DECOMP_ML_STEM, sc->data);
        }

        type = get_pair_type_window(ii + 1, jj, ptype);

        if (tmp_en == c[ii + 1][jj - (ii + 1)] + E_MLstem(type, S1[ii], -1, P) + P->MLbase) {
          *i          = *j = -1;
          *k          = ii + 1;
          *l          = jj;
          *component2 = 2;
          return 1;
        }
      }

      if (evaluate(ii, jj, ii, jj - 1, VRNA_DECOMP_ML_STEM, &hc_dat_local)) {
        int tmp_en = fij;
        if (sc) {
          if (sc->energy_up)
            tmp_en -= sc->energy_up[jj][1];

          if (sc->f)
            tmp_en -= sc->f(ii, jj, ii, jj - 1, VRNA_DECOMP_ML_STEM, sc->data);
        }

        type = get_pair_type_window(ii, jj - 1, ptype);

        if (tmp_en == c[ii][jj - 1 - ii] + E_MLstem(type, -1, S1[jj], P) + P->MLbase) {
          *i          = *j = -1;
          *k          = ii;
          *l          = jj - 1;
          *component2 = 2;
          return 1;
        }
      }

      if (evaluate(ii, jj, ii + 1, jj - 1, VRNA_DECOMP_ML_STEM, &hc_dat_local)) {
        int tmp_en = fij;
        if (sc) {
          if (sc->energy_up)
            tmp_en -= sc->energy_up[ii][1] + sc->energy_up[jj][1];

          if (sc->f)
            tmp_en -= sc->f(ii, jj, ii + 1, jj - 1, VRNA_DECOMP_ML_STEM, sc->data);
        }

        type = get_pair_type_window(ii + 1, jj - 1, ptype);

        if (tmp_en ==
            c[ii + 1][jj - 1 - (ii + 1)] + E_MLstem(type, S1[ii], S1[jj], P) + 2 * P->MLbase) {
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
  if (hc->f) {
    for (u = ii + 1 + turn; u <= jj - 2 - turn; u++) {
      if (hc->f(ii, jj, u, u + 1, VRNA_DECOMP_ML_ML_ML, hc->data)) {
        en = fML[ii][u - ii] + fML[u + 1][jj - (u + 1)];
        if (sc)
          if (sc->f)
            en += sc->f(ii, jj, u, u + 1, VRNA_DECOMP_ML_ML_ML, sc->data);

        if (fij == en) {
          *i  = ii;
          *j  = u;
          *k  = u + 1;
          *l  = jj;
          return 1;
        }
      }
    }
  } else {
    for (u = ii + 1 + turn; u <= jj - 2 - turn; u++) {
      en = fML[ii][u - ii] + fML[u + 1][jj - (u + 1)];
      if (sc)
        if (sc->f)
          en += sc->f(ii, jj, u, u + 1, VRNA_DECOMP_ML_ML_ML, sc->data);

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
    int tmp_en;
    for (u = ii + 1 + turn; u <= jj - 2 - turn; u++) {
      if (evaluate(ii, u, u + 1, jj, VRNA_DECOMP_ML_COAXIAL_ENC, &hc_dat_local)) {
        type    = rtype[get_pair_type_window(ii, u, ptype)];
        type_2  = rtype[get_pair_type_window(u + 1, jj, ptype)];

        tmp_en = c[ii][u - ii] +
                 c[u + 1][jj - (u + 1)] +
                 P->stack[type][type_2] +
                 2 * P->MLintern[1];

        if (sc)
          if (sc->f)
            tmp_en += sc->f(ii, u, u + 1, jj, VRNA_DECOMP_ML_COAXIAL, sc->data);

        if (fij == tmp_en) {
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
BT_mb_loop_split_window_comparative(vrna_fold_compound_t  *vc,
                                    int                   *i,
                                    int                   *j,
                                    int                   *k,
                                    int                   *l,
                                    int                   *component1,
                                    int                   *component2,
                                    vrna_bp_stack_t       *bp_stack,
                                    int                   *stack_count)
{
  short                     **S, **S5, **S3;
  int                       ii, jj, fij, fi, u, en, **c, **fML, **ggg, n_seq,
                            turn, with_gquad, dangle_model, cc, ss, length, type;
  vrna_param_t              *P;
  vrna_md_t                 *md;
  vrna_hc_t                 *hc;
  vrna_sc_t                 **scs;
  vrna_callback_hc_evaluate *evaluate;
  struct default_data       hc_dat_local;

  length  = vc->length;
  n_seq   = vc->n_seq;
  S       = vc->S;
  S5      = vc->S5;   /* S5[s][i] holds next base 5' of i in sequence s */
  S3      = vc->S3;   /* Sl[s][i] holds next base 3' of i in sequence s */
  P       = vc->params;
  md      = &(P->model_details);
  hc      = vc->hc;
  scs     = vc->scs;

  c             = vc->matrices->c_local;
  fML           = vc->matrices->fML_local;
  ggg           = vc->matrices->ggg_local;
  turn          = md->min_loop_size;
  with_gquad    = md->gquad;
  dangle_model  = md->dangles;

  hc_dat_local.mx         = hc->matrix;
  hc_dat_local.mx_window  = hc->matrix_local;
  hc_dat_local.hc_up      = hc->up_ml;
  hc_dat_local.sn         = vc->strand_number;

  if (hc->f) {
    evaluate            = &hc_default_user_window;
    hc_dat_local.hc_f   = hc->f;
    hc_dat_local.hc_dat = hc->data;
  } else {
    evaluate = &hc_default_window;
  }

  ii  = *i;
  jj  = *j;

  /* nibble off unpaired 3' bases */
  do {
    fij = fML[ii][jj - ii];
    fi  = INF;

    if (evaluate(ii, jj, ii, jj - 1, VRNA_DECOMP_ML_ML, &hc_dat_local)) {
      fi = fML[ii][jj - 1 - ii] +
           n_seq * P->MLbase;

      if (scs) {
        for (ss = 0; ss < n_seq; ss++)
          if (scs[ss]) {
            if (scs[ss]->energy_up)
              fi += scs[ss]->energy_up[jj][1];

            if (scs[ss]->f)
              fi += scs[ss]->f(ii, jj, ii, jj - 1, VRNA_DECOMP_ML_ML, scs[ss]->data);
          }
      }
    }

    if (--jj == 0)
      break;
  } while (fij == fi);
  jj++;

  /* nibble off unpaired 5' bases */
  do {
    fij = fML[ii][jj - ii];
    fi  = INF;

    if (evaluate(ii, jj, ii + 1, jj, VRNA_DECOMP_ML_ML, &hc_dat_local)) {
      fi = fML[ii + 1][jj - (ii + 1)] +
           n_seq * P->MLbase;

      if (scs) {
        for (ss = 0; ss < n_seq; ss++)
          if (scs[ss]) {
            if (scs[ss]->energy_up)
              fi += scs[ss]->energy_up[ii][1];

            if (scs[ss]->f)
              fi += scs[ss]->f(ii, jj, ii + 1, jj, VRNA_DECOMP_ML_ML, scs[ss]->data);
          }
      }
    }

    if (++ii == jj)
      break;
  } while (fij == fi);
  ii--;

  if (jj < ii + turn + 1) /* no more pairs */
    return 0;

  *component1 = *component2 = 1; /* split into two multi loop parts by default */

  /* 1. test for single component */

  if (with_gquad) {
    if (fij == ggg[ii][jj - ii] + n_seq * E_MLstem(0, -1, -1, P)) {
      *i  = *j = -1;
      *k  = *l = -1;
      return vrna_BT_gquad_mfe(vc, ii, jj, bp_stack, stack_count);
    }
  }

  en = c[ii][jj - ii];

  if (scs) {
    for (ss = 0; ss < n_seq; ss++)
      if (scs[ss])
        if (scs[ss]->f)
          en += scs[ss]->f(ii, jj, ii, jj, VRNA_DECOMP_ML_STEM, scs[ss]->data);
  }

  switch (dangle_model) {
    case 0:
      if (evaluate(ii, jj, ii, jj, VRNA_DECOMP_ML_STEM, &hc_dat_local)) {
        cc = 0;
        for (ss = 0; ss < n_seq; ss++) {
          type  = get_pair_type_md(S[ss][ii], S[ss][jj], md);
          cc    += E_MLstem(type, -1, -1, P);
        }

        if (fij == en + cc) {
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
        cc = 0;
        for (ss = 0; ss < n_seq; ss++) {
          type  = get_pair_type_md(S[ss][ii], S[ss][jj], md);
          cc    += E_MLstem(type, (ii > 1) ? S5[ss][ii] : -1, (jj < length) ? S3[ss][jj] : -1, P);
        }

        if (fij == en + cc) {
          *i          = *j = -1;
          *k          = ii;
          *l          = jj;
          *component2 = 2;
          return 1;
        }
      }

      break;
  }

  /* 2. Test for possible split point */
  if (hc->f) {
    for (u = ii + 1 + turn; u <= jj - 2 - turn; u++) {
      if (hc->f(ii, jj, u, u + 1, VRNA_DECOMP_ML_ML_ML, hc->data)) {
        en = fML[ii][u - ii] + fML[u + 1][jj - (u + 1)];
        if (scs) {
          for (ss = 0; ss < n_seq; ss++)
            if (scs[ss])
              if (scs[ss]->f)
                en += scs[ss]->f(ii, jj, u, u + 1, VRNA_DECOMP_ML_ML_ML, scs[ss]->data);
        }

        if (fij == en) {
          *i  = ii;
          *j  = u;
          *k  = u + 1;
          *l  = jj;
          return 1;
        }
      }
    }
  } else {
    for (u = ii + 1 + turn; u <= jj - 2 - turn; u++) {
      en = fML[ii][u - ii] + fML[u + 1][jj - (u + 1)];
      if (scs) {
        for (ss = 0; ss < n_seq; ss++)
          if (scs[ss])
            if (scs[ss]->f)
              en += scs[ss]->f(ii, jj, u, u + 1, VRNA_DECOMP_ML_ML_ML, scs[ss]->data);
      }

      if (fij == en) {
        *i  = ii;
        *j  = u;
        *k  = u + 1;
        *l  = jj;
        return 1;
      }
    }
  }

  return 0;
}


PRIVATE int
BT_mb_loop(vrna_fold_compound_t *vc,
           int                  *i,
           int                  *j,
           int                  *k,
           int                  en,
           int                  *component1,
           int                  *component2)
{
  char                      *ptype;
  short                     s5, s3, *S1;
  unsigned int              strands, *sn, *so, *ss, *se;
  int                       ij, p, q, r, e, tmp_en, *idx, turn, dangle_model,
                            *my_c, *my_fML, *my_fc, *rtype, type, type_2, tt;
  vrna_param_t              *P;
  vrna_md_t                 *md;
  vrna_hc_t                 *hc;
  vrna_sc_t                 *sc;
  vrna_callback_hc_evaluate *evaluate;
  struct default_data       hc_dat_local;
  vrna_callback_hc_evaluate *evaluate_ext;
  struct default_data       hc_dat_local_ext;

  idx           = vc->jindx;
  ij            = idx[*j] + *i;
  S1            = vc->sequence_encoding;
  P             = vc->params;
  md            = &(P->model_details);
  strands       = vc->strands;
  sn            = vc->strand_number;
  so            = vc->strand_order;
  ss            = vc->strand_start;
  se            = vc->strand_end;
  hc            = vc->hc;
  sc            = vc->sc;
  my_c          = vc->matrices->c;
  my_fML        = vc->matrices->fML;
  my_fc         = vc->matrices->fc;
  turn          = md->min_loop_size;
  ptype         = vc->ptype;
  rtype         = &(md->rtype[0]);
  type          = get_pair_type(ij, ptype);
  tt            = type;
  type          = rtype[type];
  dangle_model  = md->dangles;

  hc_dat_local.idx    = vc->jindx;
  hc_dat_local.mx     = hc->matrix;
  hc_dat_local.hc_up  = hc->up_ml;
  hc_dat_local.sn     = vc->strand_number;

  if (hc->f) {
    evaluate            = &hc_default_user;
    hc_dat_local.hc_f   = hc->f;
    hc_dat_local.hc_dat = hc->data;
  } else {
    evaluate = &hc_default;
  }

  hc_dat_local_ext.idx    = vc->jindx;
  hc_dat_local_ext.mx     = hc->matrix;
  hc_dat_local_ext.hc_up  = hc->up_ext;
  hc_dat_local_ext.sn     = sn;

  if (hc->f) {
    evaluate_ext            = &hc_default_user_ext;
    hc_dat_local_ext.hc_f   = hc->f;
    hc_dat_local_ext.hc_dat = hc->data;
  } else {
    evaluate_ext = &hc_default_ext;
  }

  p = *i + 1;
  q = *j - 1;

  r = q - turn - 1;

  if (evaluate(*i, *j, p, q, VRNA_DECOMP_PAIR_ML, &hc_dat_local)) {
    /* is it a fake multi-loop? */
    /* NOTE: do we really want to evaluate it hard-constraint-wise as a multibranch loop? */
    if (sn[*i] != sn[*j]) {
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
          if (en == e + E_ExtLoop(type, -1, -1, P))
            ii = p, jj = q;

          break;

        case 2:
          if (en == e + E_ExtLoop(type, s5, s3, P))
            ii = p, jj = q;

          break;

        default:
          if (en == e + E_ExtLoop(type, -1, -1, P)) {
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

            if (en == e + E_ExtLoop(type, -1, s3, P)) {
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

            if (en == e + E_ExtLoop(type, s5, -1, P)) {
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

            if (en == e + E_ExtLoop(type, s5, s3, P)) {
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

    /* true multi loop? */
    *component1 = *component2 = 1;  /* both components are MB loop parts by default */

    s5  = (sn[q] == sn[*j]) ? S1[q] : -1;
    s3  = (sn[*i] == sn[p]) ? S1[p] : -1;

    switch (dangle_model) {
      case 0:
        e = en - E_MLstem(type, -1, -1, P) - P->MLclosing;
        if (sc) {
          if (sc->energy_bp)
            e -= sc->energy_bp[ij];

          if (sc->f)
            e -= sc->f(*i, *j, p, q, VRNA_DECOMP_PAIR_ML, sc->data);
        }

        for (r = *i + 2 + turn; r < *j - 2 - turn; ++r) {
          if (evaluate(p, q, r, r + 1, VRNA_DECOMP_ML_ML_ML, &hc_dat_local)) {
            tmp_en = my_fML[idx[r] + p] + my_fML[idx[q] + r + 1];
            if (sc)
              if (sc->f)
                tmp_en += sc->f(p, q, r, r + 1, VRNA_DECOMP_ML_ML_ML, sc->data);

            if (e == tmp_en)
              break;
          }
        }
        break;

      case 2:
        e = en - E_MLstem(type, s5, s3, P) - P->MLclosing;
        if (sc) {
          if (sc->energy_bp)
            e -= sc->energy_bp[ij];

          if (sc->f)
            e -= sc->f(*i, *j, p, q, VRNA_DECOMP_PAIR_ML, sc->data);
        }

        for (r = p + turn + 1; r < q - turn - 1; ++r) {
          if (evaluate(p, q, r, r + 1, VRNA_DECOMP_ML_ML_ML, &hc_dat_local)) {
            tmp_en = my_fML[idx[r] + p] + my_fML[idx[q] + r + 1];
            if (sc)
              if (sc->f)
                tmp_en += sc->f(p, q, r, r + 1, VRNA_DECOMP_ML_ML_ML, sc->data);

            if (e == tmp_en)
              break;
          }
        }
        break;

      default:
        e = en - P->MLclosing;
        if (sc)
          if (sc->energy_bp)
            e -= sc->energy_bp[ij];

        for (r = p + turn + 1; r < q - turn - 1; ++r) {
          if (evaluate(p, q, r, r + 1, VRNA_DECOMP_ML_ML_ML, &hc_dat_local)) {
            tmp_en = my_fML[idx[r] + p] + my_fML[idx[q] + r + 1] + E_MLstem(type, -1, -1, P);
            if (sc) {
              if (sc->f) {
                tmp_en  += sc->f(*i, *j, p, q, VRNA_DECOMP_PAIR_ML, sc->data);
                tmp_en  += sc->f(p, q, r, r + 1, VRNA_DECOMP_ML_ML_ML, sc->data);
              }
            }

            if (e == tmp_en)
              break;
          }

          if (evaluate(*i, *j, p + 1, q, VRNA_DECOMP_PAIR_ML, &hc_dat_local)) {
            if (evaluate(p + 1, q, r, r + 1, VRNA_DECOMP_ML_ML_ML, &hc_dat_local)) {
              tmp_en = e;
              if (sc) {
                if (sc->energy_up)
                  tmp_en -= sc->energy_up[p][1];

                if (sc->f) {
                  tmp_en  -= sc->f(*i, *j, p + 1, q, VRNA_DECOMP_PAIR_ML, sc->data);
                  tmp_en  -= sc->f(p + 1, q, r, r + 1, VRNA_DECOMP_ML_ML_ML, sc->data);
                }
              }

              if (tmp_en ==
                  my_fML[idx[r] + p + 1] + my_fML[idx[q] + r + 1] +
                  E_MLstem(type, -1, s3, P) + P->MLbase) {
                p += 1;
                break;
              }
            }
          }

          if (evaluate(*i, *j, p, q - 1, VRNA_DECOMP_PAIR_ML, &hc_dat_local)) {
            if (evaluate(p, q - 1, r, r + 1, VRNA_DECOMP_ML_ML_ML, &hc_dat_local)) {
              tmp_en = e;
              if (sc) {
                if (sc->energy_up)
                  tmp_en -= sc->energy_up[q][1];

                if (sc->f) {
                  tmp_en  -= sc->f(*i, *j, p, q - 1, VRNA_DECOMP_PAIR_ML, sc->data);
                  tmp_en  -= sc->f(p, q - 1, r, r + 1, VRNA_DECOMP_ML_ML_ML, sc->data);
                }
              }

              if (tmp_en ==
                  my_fML[idx[r] + p] + my_fML[idx[q - 1] + r + 1] +
                  E_MLstem(type, s5, -1, P) + P->MLbase) {
                q -= 1;
                break;
              }
            }
          }

          if (evaluate(*i, *j, p + 1, q - 1, VRNA_DECOMP_PAIR_ML, &hc_dat_local)) {
            if (evaluate(p + 1, q - 1, r, r + 1, VRNA_DECOMP_ML_ML_ML, &hc_dat_local)) {
              tmp_en = e;
              if (sc) {
                if (sc->energy_up)
                  tmp_en -= sc->energy_up[p][1] + sc->energy_up[q][1];

                if (sc->f) {
                  tmp_en  -= sc->f(*i, *j, p + 1, q - 1, VRNA_DECOMP_PAIR_ML, sc->data);
                  tmp_en  -= sc->f(p + 1, q - 1, r, r + 1, VRNA_DECOMP_ML_ML_ML, sc->data);
                }
              }

              if (tmp_en ==
                  my_fML[idx[r] + p + 1] + my_fML[idx[q - 1] + r + 1] +
                  E_MLstem(type, s5, s3, P) + 2 * P->MLbase) {
                p += 1;
                q -= 1;
                break;
              }
            }
          }

          /* coaxial stacking of (i.j) with (i+1.r) or (r.j-1) */
          /* use MLintern[1] since coax stacked pairs don't get TerminalAU */
          if (dangle_model == 3) {
            tmp_en = e;
            if (evaluate(*i, *j, p, r, VRNA_DECOMP_ML_COAXIAL, &hc_dat_local)) {
              type_2 = rtype[get_pair_type(idx[r] + p, ptype)];

              tmp_en = my_c[idx[r] + p] +
                       P->stack[tt][type_2] +
                       my_fML[idx[q] + r + 1];

              if (sc) {
                if (sc->f) {
                  tmp_en  += sc->f(*i, *j, p, q, VRNA_DECOMP_PAIR_ML, sc->data);
                  tmp_en  += sc->f(*i, *j, p, r, VRNA_DECOMP_ML_COAXIAL, sc->data);
                }
              }

              if (e == tmp_en + 2 * P->MLintern[1]) {
                *component1 = 2;
                break;
              }
            }

            if (evaluate(*i, *j, r + 1, q, VRNA_DECOMP_ML_COAXIAL, &hc_dat_local)) {
              type_2 = rtype[get_pair_type(idx[q] + r + 1, ptype)];

              tmp_en = my_c[idx[q] + r + 1] +
                       P->stack[tt][type_2] +
                       my_fML[idx[r] + p];

              if (sc) {
                if (sc->f) {
                  tmp_en  += sc->f(*i, *j, p, q, VRNA_DECOMP_PAIR_ML, sc->data);
                  tmp_en  += sc->f(*i, *j, r + 1, q, VRNA_DECOMP_ML_COAXIAL, sc->data);
                }
              }

              if (e == tmp_en + 2 * P->MLintern[1]) {
                *component2 = 2;
                break;
              }
            }
          }
        }
        break;
    }
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


PRIVATE int
BT_mb_loop_comparative(vrna_fold_compound_t *vc,
                       int                  *i,
                       int                  *j,
                       int                  *k,
                       int                  en,
                       int                  *component1,
                       int                  *component2)
{
  short                     **S, **S5, **S3;
  int                       ij, p, q, r, e, tmp_en, *idx, turn, dangle_model,
                            *my_fML, ss, n_seq, tt;
  vrna_param_t              *P;
  vrna_md_t                 *md;
  vrna_hc_t                 *hc;
  vrna_sc_t                 **scs;
  vrna_callback_hc_evaluate *evaluate;
  struct default_data       hc_dat_local;

  n_seq         = vc->n_seq;
  S             = vc->S;
  S5            = vc->S5;
  S3            = vc->S3;
  idx           = vc->jindx;
  ij            = idx[*j] + *i;
  P             = vc->params;
  md            = &(P->model_details);
  hc            = vc->hc;
  scs           = vc->scs;
  my_fML        = vc->matrices->fML;
  turn          = md->min_loop_size;
  dangle_model  = md->dangles;

  hc_dat_local.idx    = vc->jindx;
  hc_dat_local.mx     = hc->matrix;
  hc_dat_local.hc_up  = hc->up_ml;
  hc_dat_local.sn     = sn;

  if (hc->f) {
    evaluate            = &hc_default_user;
    hc_dat_local.hc_f   = hc->f;
    hc_dat_local.hc_dat = hc->data;
  } else {
    evaluate = &hc_default;
  }

  p = *i + 1;
  q = *j - 1;

  r = q - turn - 1;

  if (evaluate(*i, *j, p, q, VRNA_DECOMP_PAIR_ML, &hc_dat_local)) {
    *component1 = *component2 = 1;  /* both components are MB loop parts by default */

    e = en - n_seq * P->MLclosing;

    switch (dangle_model) {
      case 0:
        for (ss = 0; ss < n_seq; ss++) {
          tt  = get_pair_type_md(S[ss][*j], S[ss][*i], md);
          e   -= E_MLstem(tt, -1, -1, P);
        }
        break;

      case 2:
        for (ss = 0; ss < n_seq; ss++) {
          tt  = get_pair_type_md(S[ss][*j], S[ss][*i], md);
          e   -= E_MLstem(tt, S5[ss][*j], S3[ss][*i], P);
        }
        break;
    }

    if (scs) {
      for (ss = 0; ss < n_seq; ss++) {
        if (scs[ss]) {
          if (scs[ss]->energy_bp)
            e -= scs[ss]->energy_bp[ij];

          if (scs[ss]->f)
            e -= scs[ss]->f(*i, *j, p, q, VRNA_DECOMP_PAIR_ML, scs[ss]->data);
        }
      }
    }

    for (r = p + turn + 1; r < q - turn - 1; ++r) {
      if (evaluate(p, q, r, r + 1, VRNA_DECOMP_ML_ML_ML, &hc_dat_local)) {
        tmp_en = my_fML[idx[r] + p] + my_fML[idx[q] + r + 1];
        if (scs) {
          for (ss = 0; ss < n_seq; ss++) {
            if (scs[ss])
              if (scs[ss]->f)
                tmp_en += scs[ss]->f(p, q, r, r + 1, VRNA_DECOMP_ML_ML_ML, scs[ss]->data);
          }
        }

        if (e == tmp_en)
          break;
      }
    }
  }

  if (r <= *j - turn - 3) {
    *i  = p;
    *k  = r;
    *j  = q;
    return 1;
  }

  return 0;
}


PRIVATE int
BT_mb_loop_window(vrna_fold_compound_t  *vc,
                  int                   *i,
                  int                   *j,
                  int                   *k,
                  int                   en,
                  int                   *component1,
                  int                   *component2)
{
  char                      **ptype;
  short                     s5, s3, *S1;
  unsigned int              *sn;
  int                       p, q, r, e, tmp_en, turn, dangle_model,
                            **c, **fML, *rtype, type, type_2, tt;
  vrna_param_t              *P;
  vrna_md_t                 *md;
  vrna_hc_t                 *hc;
  vrna_sc_t                 *sc;
  vrna_callback_hc_evaluate *evaluate;
  struct default_data       hc_dat_local;

  S1            = vc->sequence_encoding;
  P             = vc->params;
  md            = &(P->model_details);
  sn            = vc->strand_number;
  hc            = vc->hc;
  sc            = vc->sc;
  c             = vc->matrices->c_local;
  fML           = vc->matrices->fML_local;
  turn          = md->min_loop_size;
  ptype         = vc->ptype_local;
  rtype         = &(md->rtype[0]);
  type          = get_pair_type_window(*i, *j, ptype);
  tt            = type;
  type          = rtype[type];
  dangle_model  = md->dangles;

  hc_dat_local.idx        = vc->jindx;
  hc_dat_local.mx         = hc->matrix;
  hc_dat_local.mx_window  = hc->matrix_local;
  hc_dat_local.hc_up      = hc->up_ml;
  hc_dat_local.sn         = vc->strand_number;

  if (hc->f) {
    evaluate            = &hc_default_user_window;
    hc_dat_local.hc_f   = hc->f;
    hc_dat_local.hc_dat = hc->data;
  } else {
    evaluate = &hc_default_window;
  }

  p = *i + 1;
  q = *j - 1;

  r = q - turn - 1;

  if (evaluate(*i, *j, p, q, VRNA_DECOMP_PAIR_ML, &hc_dat_local)) {
    *component1 = *component2 = 1;  /* both components are MB loop parts by default */

    s5  = (sn[q] == sn[*j]) ? S1[q] : -1;
    s3  = (sn[*i] == sn[p]) ? S1[p] : -1;

    switch (dangle_model) {
      case 0:
        e = en - E_MLstem(type, -1, -1, P) - P->MLclosing;
        if (sc) {
          if (sc->energy_bp_local)
            e -= sc->energy_bp_local[*i][*j - *i];

          if (sc->f)
            e -= sc->f(*i, *j, p, q, VRNA_DECOMP_PAIR_ML, sc->data);
        }

        for (r = *i + 2 + turn; r < *j - 2 - turn; ++r) {
          if (evaluate(p, q, r, r + 1, VRNA_DECOMP_ML_ML_ML, &hc_dat_local)) {
            tmp_en = fML[p][r - p] +
                     fML[r + 1][q - (r + 1)];

            if (sc)
              if (sc->f)
                tmp_en += sc->f(p, q, r, r + 1, VRNA_DECOMP_ML_ML_ML, sc->data);

            if (e == tmp_en)
              break;
          }
        }
        break;

      case 2:
        e = en - E_MLstem(type, s5, s3, P) - P->MLclosing;
        if (sc) {
          if (sc->energy_bp_local)
            e -= sc->energy_bp_local[*i][*j - *i];

          if (sc->f)
            e -= sc->f(*i, *j, p, q, VRNA_DECOMP_PAIR_ML, sc->data);
        }

        for (r = p + turn + 1; r < q - turn - 1; ++r) {
          if (evaluate(p, q, r, r + 1, VRNA_DECOMP_ML_ML_ML, &hc_dat_local)) {
            tmp_en = fML[p][r - p] +
                     fML[r + 1][q - (r + 1)];

            if (sc)
              if (sc->f)
                tmp_en += sc->f(p, q, r, r + 1, VRNA_DECOMP_ML_ML_ML, sc->data);

            if (e == tmp_en)
              break;
          }
        }
        break;

      default:
        e = en - P->MLclosing;

        for (r = p + turn + 1; r < q - turn - 1; ++r) {
          if (evaluate(p, q, r, r + 1, VRNA_DECOMP_ML_ML_ML, &hc_dat_local)) {
            tmp_en = fML[p][r - p] +
                     fML[r + 1][q - (r + 1)] +
                     E_MLstem(type, -1, -1, P);

            if (sc) {
              if (sc->energy_bp)
                tmp_en += sc->energy_bp_local[*i][*j - *i];

              if (sc->f) {
                tmp_en  += sc->f(*i, *j, p, q, VRNA_DECOMP_PAIR_ML, sc->data);
                tmp_en  += sc->f(p, q, r, r + 1, VRNA_DECOMP_ML_ML_ML, sc->data);
              }
            }

            if (e == tmp_en)
              break;
          }

          if (evaluate(*i, *j, p + 1, q, VRNA_DECOMP_PAIR_ML, &hc_dat_local)) {
            if (evaluate(p + 1, q, r, r + 1, VRNA_DECOMP_ML_ML_ML, &hc_dat_local)) {
              tmp_en = fML[p + 1][r - (p + 1)] +
                       fML[r + 1][q - (r + 1)] +
                       E_MLstem(type, -1, s3, P) +
                       P->MLbase;

              if (sc) {
                if (sc->energy_up)
                  tmp_en += sc->energy_up[p][1];

                if (sc->energy_bp)
                  tmp_en += sc->energy_bp_local[*i][*j - *i];

                if (sc->f) {
                  tmp_en  += sc->f(*i, *j, p + 1, q, VRNA_DECOMP_PAIR_ML, sc->data);
                  tmp_en  += sc->f(p + 1, q, r, r + 1, VRNA_DECOMP_ML_ML_ML, sc->data);
                }
              }

              if (e == tmp_en) {
                p += 1;
                break;
              }
            }
          }

          if (evaluate(*i, *j, p, q - 1, VRNA_DECOMP_PAIR_ML, &hc_dat_local)) {
            if (evaluate(p, q - 1, r, r + 1, VRNA_DECOMP_ML_ML_ML, &hc_dat_local)) {
              tmp_en = fML[p][r - p] +
                       fML[r + 1][q - 1 - (r + 1)] +
                       E_MLstem(type, s5, -1, P) +
                       P->MLbase;

              if (sc) {
                if (sc->energy_up)
                  tmp_en += sc->energy_up[q][1];

                if (sc->energy_bp)
                  tmp_en += sc->energy_bp_local[*i][*j - *i];

                if (sc->f) {
                  tmp_en  += sc->f(*i, *j, p, q - 1, VRNA_DECOMP_PAIR_ML, sc->data);
                  tmp_en  += sc->f(p, q - 1, r, r + 1, VRNA_DECOMP_ML_ML_ML, sc->data);
                }
              }

              if (e == tmp_en) {
                q -= 1;
                break;
              }
            }
          }

          if (evaluate(*i, *j, p + 1, q - 1, VRNA_DECOMP_PAIR_ML, &hc_dat_local)) {
            if (evaluate(p + 1, q - 1, r, r + 1, VRNA_DECOMP_ML_ML_ML, &hc_dat_local)) {
              tmp_en = fML[p + 1][r - (p + 1)] +
                       fML[r + 1][q - 1 - (r + 1)] +
                       E_MLstem(type, s5, s3, P) +
                       2 * P->MLbase;

              if (sc) {
                if (sc->energy_up)
                  tmp_en += sc->energy_up[p][1] +
                            sc->energy_up[q][1];

                if (sc->energy_bp)
                  tmp_en += sc->energy_bp_local[*i][*j - *i];

                if (sc->f) {
                  tmp_en  += sc->f(*i, *j, p + 1, q - 1, VRNA_DECOMP_PAIR_ML, sc->data);
                  tmp_en  += sc->f(p + 1, q - 1, r, r + 1, VRNA_DECOMP_ML_ML_ML, sc->data);
                }
              }

              if (e == tmp_en) {
                p += 1;
                q -= 1;
                break;
              }
            }
          }

          /* coaxial stacking of (i.j) with (i+1.r) or (r.j-1) */
          /* use MLintern[1] since coax stacked pairs don't get TerminalAU */
          if (dangle_model == 3) {
            tmp_en = e;
            if (evaluate(*i, *j, p, r, VRNA_DECOMP_ML_COAXIAL, &hc_dat_local)) {
              type_2 = rtype[get_pair_type_window(p, r, ptype)];

              tmp_en = c[p][r - p] +
                       P->stack[tt][type_2] +
                       fML[r + 1][q - (r + 1)];

              if (sc) {
                if (sc->energy_bp)
                  tmp_en += sc->energy_bp_local[*i][*j - *i];

                if (sc->f) {
                  tmp_en  += sc->f(*i, *j, p, q, VRNA_DECOMP_PAIR_ML, sc->data);
                  tmp_en  += sc->f(*i, *j, p, r, VRNA_DECOMP_ML_COAXIAL, sc->data);
                }
              }

              if (e == tmp_en + 2 * P->MLintern[1]) {
                *component1 = 2;
                break;
              }
            }

            if (evaluate(*i, *j, r + 1, q, VRNA_DECOMP_ML_COAXIAL, &hc_dat_local)) {
              type_2 = rtype[get_pair_type_window(r + 1, q, ptype)];

              tmp_en = c[r + 1][q - (r + 1)] +
                       P->stack[tt][type_2] +
                       fML[p][r - p];

              if (sc) {
                if (sc->energy_bp)
                  tmp_en += sc->energy_bp_local[*i][*j - *i];

                if (sc->f) {
                  tmp_en  += sc->f(*i, *j, p, q, VRNA_DECOMP_PAIR_ML, sc->data);
                  tmp_en  += sc->f(*i, *j, r + 1, q, VRNA_DECOMP_ML_COAXIAL, sc->data);
                }
              }

              if (e == tmp_en + 2 * P->MLintern[1]) {
                *component2 = 2;
                break;
              }
            }
          }
        }
        break;
    }
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


PRIVATE int
BT_mb_loop_window_comparative(vrna_fold_compound_t  *vc,
                              int                   *i,
                              int                   *j,
                              int                   *k,
                              int                   en,
                              int                   *component1,
                              int                   *component2)
{
  short                     **S, **S5, **S3;
  int                       tt, p, q, r, e, tmp_en, turn, dangle_model,
                            **fML, mm, n_seq, ss;
  vrna_param_t              *P;
  vrna_md_t                 *md;
  vrna_hc_t                 *hc;
  vrna_sc_t                 **scs;
  vrna_callback_hc_evaluate *evaluate;
  struct default_data       hc_dat_local;

  n_seq         = vc->n_seq;
  S             = vc->S;
  S5            = vc->S5;   /* S5[s][i] holds next base 5' of i in sequence s */
  S3            = vc->S3;   /* Sl[s][i] holds next base 3' of i in sequence s */
  P             = vc->params;
  md            = &(P->model_details);
  hc            = vc->hc;
  scs           = vc->scs;
  fML           = vc->matrices->fML_local;
  turn          = md->min_loop_size;
  dangle_model  = md->dangles;

  hc_dat_local.idx        = vc->jindx;
  hc_dat_local.mx         = hc->matrix;
  hc_dat_local.mx_window  = hc->matrix_local;
  hc_dat_local.hc_up      = hc->up_ml;
  hc_dat_local.sn         = vc->strand_number;

  if (hc->f) {
    evaluate            = &hc_default_user_window;
    hc_dat_local.hc_f   = hc->f;
    hc_dat_local.hc_dat = hc->data;
  } else {
    evaluate = &hc_default_window;
  }

  p = *i + 1;
  q = *j - 1;

  r = q - turn - 1;

  if (evaluate(*i, *j, p, q, VRNA_DECOMP_PAIR_ML, &hc_dat_local)) {
    mm = n_seq * P->MLclosing;

    *component1 = *component2 = 1;  /* both components are MB loop parts by default */

    switch (dangle_model) {
      case 0:
        for (ss = 0; ss < n_seq; ss++) {
          tt  = get_pair_type_md(S[ss][*j], S[ss][*i], md);
          mm  += E_MLstem(tt, -1, -1, P);
        }

        e = en - mm;
        if (scs) {
          for (ss = 0; ss < n_seq; ss++)
            if (scs[ss])
              if (scs[ss]->f)
                e -= scs[ss]->f(*i, *j, p, q, VRNA_DECOMP_PAIR_ML, scs[ss]->data);
        }

        for (r = *i + 2 + turn; r < *j - 2 - turn; ++r) {
          if (evaluate(p, q, r, r + 1, VRNA_DECOMP_ML_ML_ML, &hc_dat_local)) {
            tmp_en = fML[p][r - p] +
                     fML[r + 1][q - (r + 1)];

            if (scs) {
              for (ss = 0; ss < n_seq; ss++)
                if (scs[ss])
                  if (scs[ss]->f)
                    tmp_en += scs[ss]->f(p, q, r, r + 1, VRNA_DECOMP_ML_ML_ML, scs[ss]->data);
            }

            if (e == tmp_en)
              break;
          }
        }
        break;

      case 2:
        for (ss = 0; ss < n_seq; ss++) {
          tt  = get_pair_type_md(S[ss][*j], S[ss][*i], md);
          mm  += E_MLstem(tt, S5[ss][*j], S3[ss][*i], P);
        }

        e = en - mm;
        if (scs) {
          for (ss = 0; ss < n_seq; ss++)
            if (scs[ss])
              if (scs[ss]->f)
                e -= scs[ss]->f(*i, *j, p, q, VRNA_DECOMP_PAIR_ML, scs[ss]->data);
        }

        for (r = p + turn + 1; r < q - turn - 1; ++r) {
          if (evaluate(p, q, r, r + 1, VRNA_DECOMP_ML_ML_ML, &hc_dat_local)) {
            tmp_en = fML[p][r - p] +
                     fML[r + 1][q - (r + 1)];

            if (scs) {
              for (ss = 0; ss < n_seq; ss++)
                if (scs[ss])
                  if (scs[ss]->f)
                    tmp_en += scs[ss]->f(p, q, r, r + 1, VRNA_DECOMP_ML_ML_ML, scs[ss]->data);
            }

            if (e == tmp_en)
              break;
          }
        }
        break;
    }
  }

  if (r <= *j - turn - 3) {
    *i  = p;
    *k  = r;
    *j  = q;
    return 1;
  }

  return 0;
}
