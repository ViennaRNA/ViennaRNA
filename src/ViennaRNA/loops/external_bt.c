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
#include "ViennaRNA/gquad.h"
#include "ViennaRNA/structured_domains.h"
#include "ViennaRNA/unstructured_domains.h"
#include "ViennaRNA/loops/external.h"

#ifdef __GNUC__
# define INLINE inline
#else
# define INLINE
#endif

#include "external_hc.inc"
#include "external_sc.inc"

/*
 #################################
 # PRIVATE FUNCTION DECLARATIONS #
 #################################
 */

PRIVATE int
BT_ext_loop_f5(vrna_fold_compound_t *fc,
               int                  *k,
               int                  *i,
               int                  *j,
               vrna_bp_stack_t      *bp_stack,
               int                  *stack_count);


PRIVATE int
BT_ext_loop_f5_comparative(vrna_fold_compound_t *fc,
                           int                  *k,
                           int                  *i,
                           int                  *j,
                           vrna_bp_stack_t      *bp_stack,
                           int                  *stack_count);


PRIVATE int
BT_ext_loop_f3(vrna_fold_compound_t *fc,
               int                  *k,
               int                  maxdist,
               int                  *i,
               int                  *j,
               vrna_bp_stack_t      *bp_stack,
               int                  *stack_count);


PRIVATE int
BT_ext_loop_f3_comparative(vrna_fold_compound_t *fc,
                           int                  *k,
                           int                  maxdist,
                           int                  *i,
                           int                  *j,
                           vrna_bp_stack_t      *bp_stack,
                           int                  *stack_count);


PRIVATE int
BT_ext_loop_f3_pp(vrna_fold_compound_t  *fc,
                  int                   *i,
                  int                   maxj);


PRIVATE int
BT_ext_loop_f3_pp_comparative(vrna_fold_compound_t  *fc,
                              int                   *i,
                              int                   maxj);


/*
 #################################
 # BEGIN OF FUNCTION DEFINITIONS #
 #################################
 */
PUBLIC int
vrna_BT_ext_loop_f5(vrna_fold_compound_t  *fc,
                    int                   *k,
                    int                   *i,
                    int                   *j,
                    vrna_bp_stack_t       *bp_stack,
                    int                   *stack_count)
{
  if (fc) {
    switch (fc->type) {
      case VRNA_FC_TYPE_SINGLE:
        return BT_ext_loop_f5(fc, k, i, j, bp_stack, stack_count);
        break;

      case VRNA_FC_TYPE_COMPARATIVE:
        return BT_ext_loop_f5_comparative(fc, k, i, j, bp_stack, stack_count);
        break;
    }
  }

  return -1;
}


PUBLIC int
vrna_BT_ext_loop_f3(vrna_fold_compound_t  *fc,
                    int                   *k,
                    int                   maxdist,
                    int                   *i,
                    int                   *j,
                    vrna_bp_stack_t       *bp_stack,
                    int                   *stack_count)
{
  if (fc) {
    switch (fc->type) {
      case VRNA_FC_TYPE_SINGLE:
        return BT_ext_loop_f3(fc, k, maxdist, i, j, bp_stack, stack_count);
        break;

      case VRNA_FC_TYPE_COMPARATIVE:
        return BT_ext_loop_f3_comparative(fc, k, maxdist, i, j, bp_stack, stack_count);
        break;
    }
  }

  return -1;
}


PUBLIC int
vrna_BT_ext_loop_f3_pp(vrna_fold_compound_t *fc,
                       int                  *i,
                       int                  maxj)
{
  if (fc) {
    switch (fc->type) {
      case VRNA_FC_TYPE_SINGLE:
        return BT_ext_loop_f3_pp(fc, i, maxj);
        break;

      case VRNA_FC_TYPE_COMPARATIVE:
        return BT_ext_loop_f3_pp_comparative(fc, i, maxj);
        break;
    }
  }

  return -1;
}


PRIVATE int
BT_ext_loop_f5(vrna_fold_compound_t *fc,
               int                  *k,
               int                  *i,
               int                  *j,
               vrna_bp_stack_t      *bp_stack,
               int                  *stack_count)
{
  char                      *ptype;
  short                     mm5, mm3, *S1;
  unsigned int              *sn, type;
  int                       length, fij, fi, jj, u, en, e, *my_f5, *my_c, *my_ggg, *idx,
                            dangle_model, turn, with_gquad, cnt, ii, with_ud;
  vrna_param_t              *P;
  vrna_md_t                 *md;
  vrna_sc_t                 *sc;
  vrna_ud_t                 *domains_up;
  vrna_callback_hc_evaluate *evaluate;
  struct default_data       hc_dat_local;

  length        = fc->length;
  P             = fc->params;
  md            = &(P->model_details);
  sn            = fc->strand_number;
  sc            = fc->sc;
  my_f5         = fc->matrices->f5;
  my_c          = fc->matrices->c;
  my_ggg        = fc->matrices->ggg;
  domains_up    = fc->domains_up;
  idx           = fc->jindx;
  ptype         = fc->ptype;
  S1            = fc->sequence_encoding;
  dangle_model  = md->dangles;
  turn          = md->min_loop_size;
  with_gquad    = md->gquad;
  with_ud       = (domains_up && domains_up->energy_cb) ? 1 : 0;
  evaluate      = prepare_hc_default(fc, &hc_dat_local);

  jj = *k;

  /* nibble off unpaired 3' stretches harboring bound ligands (interspersed with unpaired nucleotides) */
  if (with_ud) {
    do {
      fij = my_f5[jj];
      fi  = INF;

      /* try nibble off one unpaired nucleotide first */
      if (evaluate(1, jj, 1, jj - 1, VRNA_DECOMP_EXT_EXT, &hc_dat_local)) {
        fi = my_f5[jj - 1];

        if (sc) {
          if (sc->energy_up)
            fi += sc->energy_up[jj][1];

          if (sc->f)
            fi += sc->f(1, jj, 1, jj - 1, VRNA_DECOMP_EXT_EXT, sc->data);
        }

        if (jj == 1) {
          /* no more pairs */
          *i  = *j = -1;
          *k  = 0;
          return 1;
        }

        if (fij == fi) {
          jj--;
          continue;
        }
      }

      /* next, try nibble off a ligand */
      for (cnt = 0; cnt < domains_up->uniq_motif_count; cnt++) {
        u   = domains_up->uniq_motif_size[cnt];
        ii  = jj - u + 1;
        if ((ii > 0) && evaluate(1, jj, 1, jj - u, VRNA_DECOMP_EXT_EXT, &hc_dat_local)) {
          en = domains_up->energy_cb(fc,
                                     ii,
                                     jj,
                                     VRNA_UNSTRUCTURED_DOMAIN_EXT_LOOP | VRNA_UNSTRUCTURED_DOMAIN_MOTIF,
                                     domains_up->data);
          if (sc) {
            if (sc->energy_up)
              en += sc->energy_up[ii][u];

            if (sc->f)
              en += sc->f(1, jj, 1, jj - u, VRNA_DECOMP_EXT_EXT, sc->data);
          }

          fi  = my_f5[ii - 1];
          fi  += en;

          if (fij == fi) {
            /* skip remaining motifs after first hit */
            jj = ii - 1;
            break;
          }
        }
      }

      if (jj == 0) {
        /* no more pairs */
        *i  = *j = -1;
        *k  = 0;
        return 1;
      }
    } while (fij == fi);
  } else {
    /* nibble off unpaired 3' bases */
    do {
      fij = my_f5[jj];
      fi  = INF;

      if (evaluate(1, jj, 1, jj - 1, VRNA_DECOMP_EXT_EXT, &hc_dat_local)) {
        fi = my_f5[jj - 1];

        if (sc) {
          if (sc->energy_up)
            fi += sc->energy_up[jj][1];

          if (sc->f)
            fi += sc->f(1, jj, 1, jj - 1, VRNA_DECOMP_EXT_EXT, sc->data);
        }
      }

      if (--jj == 0)
        break;
    } while (fij == fi);
    jj++;
  }

  if (jj < turn + 2) {
    /* no more pairs */
    *i  = *j = -1;
    *k  = 0;
    return 1;
  }

  /* must have found a decomposition */
  switch (dangle_model) {
    case 0:   /* j is paired. Find pairing partner */
      for (u = jj - turn - 1; u >= 1; u--) {
        if (with_gquad) {
          if (fij == my_f5[u - 1] + my_ggg[idx[jj] + u]) {
            *i  = *j = -1;
            *k  = u - 1;
            return vrna_BT_gquad_mfe(fc, u, jj, bp_stack, stack_count);
          }
        }

        if (evaluate(1, jj, u - 1, u, VRNA_DECOMP_EXT_EXT_STEM, &hc_dat_local)) {
          type = vrna_get_ptype(idx[jj] + u, ptype);

          en = my_c[idx[jj] + u];
          if (sc)
            if (sc->f)
              en += sc->f(1, jj, u - 1, u, VRNA_DECOMP_EXT_EXT_STEM, sc->data);

          if (sn[jj] != sn[u])
            en += P->DuplexInit;

          if (fij == vrna_E_ext_stem(type, -1, -1, P) + en + my_f5[u - 1]) {
            *i                            = u;
            *j                            = jj;
            *k                            = u - 1;
            bp_stack[++(*stack_count)].i  = u;
            bp_stack[(*stack_count)].j    = jj;
            return 1;
          }
        }
      }
      break;

    case 2:
      mm3 = ((jj < length) && (sn[jj + 1] == sn[jj])) ? S1[jj + 1] : -1;
      for (u = jj - turn - 1; u >= 1; u--) {
        if (with_gquad) {
          if (fij == my_f5[u - 1] + my_ggg[idx[jj] + u]) {
            *i  = *j = -1;
            *k  = u - 1;
            return vrna_BT_gquad_mfe(fc, u, jj, bp_stack, stack_count);
          }
        }

        if (evaluate(1, jj, u - 1, u, VRNA_DECOMP_EXT_EXT_STEM, &hc_dat_local)) {
          mm5   = ((u > 1) && (sn[u] == sn[u - 1])) ? S1[u - 1] : -1;
          type  = vrna_get_ptype(idx[jj] + u, ptype);

          en = my_c[idx[jj] + u];
          if (sc)
            if (sc->f)
              en += sc->f(1, jj, u - 1, u, VRNA_DECOMP_EXT_EXT_STEM, sc->data);

          if (sn[jj] != sn[u])
            en += P->DuplexInit;

          if (fij == vrna_E_ext_stem(type, mm5, mm3, P) + en + my_f5[u - 1]) {
            *i                            = u;
            *j                            = jj;
            *k                            = u - 1;
            bp_stack[++(*stack_count)].i  = u;
            bp_stack[(*stack_count)].j    = jj;
            return 1;
          }
        }
      }
      break;

    default:
      if (with_gquad) {
        if (fij == my_ggg[idx[jj] + 1]) {
          *i  = *j = -1;
          *k  = 0;
          return vrna_BT_gquad_mfe(fc, 1, jj, bp_stack, stack_count);
        }
      }

      if (evaluate(1, jj, 1, jj, VRNA_DECOMP_EXT_STEM, &hc_dat_local)) {
        type = vrna_get_ptype(idx[jj] + 1, ptype);

        en = my_c[idx[jj] + 1];
        if (sc)
          if (sc->f)
            en += sc->f(1, jj, 1, jj, VRNA_DECOMP_EXT_STEM, sc->data);

        if (sn[jj] != sn[1])
          en += P->DuplexInit;

        if (fij == en + vrna_E_ext_stem(type, -1, -1, P)) {
          *i                            = 1;
          *j                            = jj;
          *k                            = 0;
          bp_stack[++(*stack_count)].i  = 1;
          bp_stack[(*stack_count)].j    = jj;
          return 1;
        }
      }

      if (evaluate(1, jj, 1, jj - 1, VRNA_DECOMP_EXT_STEM, &hc_dat_local)) {
        if (sn[jj] == sn[jj - 1]) {
          mm3   = S1[jj];
          type  = vrna_get_ptype(idx[jj - 1] + 1, ptype);

          en = my_c[idx[jj - 1] + 1];
          if (sc) {
            if (sc->energy_up)
              en += sc->energy_up[jj][1];

            if (sc->f)
              en += sc->f(1, jj, 1, jj - 1, VRNA_DECOMP_EXT_STEM, sc->data);
          }

          if (sn[jj - 1] != sn[1])
            en += P->DuplexInit;

          if (fij == en + vrna_E_ext_stem(type, -1, mm3, P)) {
            *i                            = 1;
            *j                            = jj - 1;
            *k                            = 0;
            bp_stack[++(*stack_count)].i  = 1;
            bp_stack[(*stack_count)].j    = jj - 1;
            return 1;
          }
        }
      }

      for (u = jj - turn - 1; u > 1; u--) {
        if (with_gquad) {
          if (fij == my_f5[u - 1] + my_ggg[idx[jj] + u]) {
            *i  = *j = -1;
            *k  = u - 1;
            return vrna_BT_gquad_mfe(fc, u, jj, bp_stack, stack_count);
          }
        }

        type = vrna_get_ptype(idx[jj] + u, ptype);

        en = my_c[idx[jj] + u];
        if (sn[jj] != sn[u])
          en += P->DuplexInit;

        if (evaluate(1, jj, u - 1, u, VRNA_DECOMP_EXT_EXT_STEM, &hc_dat_local)) {
          e = my_f5[u - 1] +
              en +
              vrna_E_ext_stem(type, -1, -1, P);

          if (sc)
            if (sc->f)
              e += sc->f(1, jj, u - 1, u, VRNA_DECOMP_EXT_EXT_STEM, sc->data);

          if (fij == e) {
            *i                            = u;
            *j                            = jj;
            *k                            = u - 1;
            bp_stack[++(*stack_count)].i  = u;
            bp_stack[(*stack_count)].j    = jj;
            return 1;
          }
        }

        if (evaluate(1, jj, u - 2, u, VRNA_DECOMP_EXT_EXT_STEM, &hc_dat_local)) {
          if (sn[u] == sn[u - 1]) {
            mm5 = S1[u - 1];
            e   = my_f5[u - 2] +
                  en +
                  vrna_E_ext_stem(type, mm5, -1, P);

            if (sc) {
              if (sc->energy_up)
                e += sc->energy_up[u - 1][1];

              if (sc->f)
                e += sc->f(1, jj, u - 2, u, VRNA_DECOMP_EXT_EXT_STEM, sc->data);
            }

            if (fij == e) {
              *i                            = u;
              *j                            = jj;
              *k                            = u - 2;
              bp_stack[++(*stack_count)].i  = u;
              bp_stack[(*stack_count)].j    = jj;
              return 1;
            }
          }
        }

        type = vrna_get_ptype(idx[jj - 1] + u, ptype);

        en = my_c[idx[jj - 1] + u];
        if (sn[jj - 1] != sn[u])
          en += P->DuplexInit;

        mm5 = (sn[u] == sn[u - 1]) ? S1[u - 1] : -1;
        mm3 = (sn[jj] == sn[jj - 1]) ? S1[jj] : -1;

        if (evaluate(1, jj, u - 1, u, VRNA_DECOMP_EXT_EXT_STEM1, &hc_dat_local)) {
          e = my_f5[u - 1] +
              en +
              vrna_E_ext_stem(type, -1, mm3, P);

          if (sc) {
            if (sc->energy_up)
              e += sc->energy_up[jj][1];

            if (sc->f)
              e += sc->f(1, jj, u - 1, u, VRNA_DECOMP_EXT_EXT_STEM1, sc->data);
          }

          if (fij == e) {
            *i                            = u;
            *j                            = jj - 1;
            *k                            = u - 1;
            bp_stack[++(*stack_count)].i  = u;
            bp_stack[(*stack_count)].j    = jj - 1;
            return 1;
          }
        }

        if (evaluate(1, jj, u - 2, u, VRNA_DECOMP_EXT_EXT_STEM1, &hc_dat_local)) {
          e = my_f5[u - 2] + en + vrna_E_ext_stem(type, mm5, mm3, P);
          if (sc) {
            if (sc->energy_up)
              e += sc->energy_up[jj][1] +
                   sc->energy_up[u - 1][1];

            if (sc->f)
              e += sc->f(1, jj, u - 2, u, VRNA_DECOMP_EXT_EXT_STEM1, sc->data);
          }

          if (fij == e) {
            *i                            = u;
            *j                            = jj - 1;
            *k                            = u - 2;
            bp_stack[++(*stack_count)].i  = u;
            bp_stack[(*stack_count)].j    = jj - 1;
            return 1;
          }
        }
      }
      break;
  }

  return 0;
}


PRIVATE int
BT_ext_loop_f5_comparative(vrna_fold_compound_t *fc,
                           int                  *k,
                           int                  *i,
                           int                  *j,
                           vrna_bp_stack_t      *bp_stack,
                           int                  *stack_count)
{
  unsigned int              **a2s;
  short                     **S, **S5, **S3;
  unsigned int              tt;
  int                       fij, fi, jj, u, en, *my_f5, *my_c, *my_ggg, *idx,
                            dangle_model, turn, with_gquad, n_seq, ss, mm5, mm3;
  vrna_param_t              *P;
  vrna_md_t                 *md;
  vrna_sc_t                 **scs;
  vrna_callback_hc_evaluate *evaluate;
  struct default_data       hc_dat_local;

  n_seq         = fc->n_seq;
  S             = fc->S;
  S5            = fc->S5;
  S3            = fc->S3;
  a2s           = fc->a2s;
  P             = fc->params;
  md            = &(P->model_details);
  scs           = fc->scs;
  my_f5         = fc->matrices->f5;
  my_c          = fc->matrices->c;
  my_ggg        = fc->matrices->ggg;
  idx           = fc->jindx;
  dangle_model  = md->dangles;
  turn          = md->min_loop_size;
  with_gquad    = md->gquad;
  evaluate      = prepare_hc_default(fc, &hc_dat_local);

  jj = *k;

  /* nibble off unpaired 3' bases */
  do {
    fij = my_f5[jj];
    fi  = INF;

    if (evaluate(1, jj, 1, jj - 1, VRNA_DECOMP_EXT_EXT, &hc_dat_local)) {
      fi = my_f5[jj - 1];

      if (scs) {
        for (ss = 0; ss < n_seq; ss++)
          if (scs[ss]) {
            if (scs[ss]->energy_up)
              fi += scs[ss]->energy_up[a2s[ss][jj]][1];

            if (scs[ss]->f)
              fi += scs[ss]->f(1, jj, 1, jj - 1, VRNA_DECOMP_EXT_EXT, scs[ss]->data);
          }
      }
    }

    if (--jj == 0)
      break;
  } while (fij == fi);
  jj++;

  if (jj < turn + 2) {
    /* no more pairs */
    *i  = *j = -1;
    *k  = 0;
    return 1;
  }

  /* must have found a decomposition */
  switch (dangle_model) {
    case 0:   /* j is paired. Find pairing partner */
      for (u = jj - turn - 1; u >= 1; u--) {
        if (with_gquad) {
          if (fij == my_f5[u - 1] + my_ggg[idx[jj] + u]) {
            *i  = *j = -1;
            *k  = u - 1;
            return vrna_BT_gquad_mfe(fc, u, jj, bp_stack, stack_count);
          }
        }

        if (evaluate(1, jj, u - 1, u, VRNA_DECOMP_EXT_EXT_STEM, &hc_dat_local)) {
          en = my_c[idx[jj] + u] +
               my_f5[u - 1];

          for (ss = 0; ss < n_seq; ss++) {
            tt  = vrna_get_ptype_md(S[ss][u], S[ss][jj], md);
            en  += vrna_E_ext_stem(tt, -1, -1, P);
          }

          if (scs) {
            for (ss = 0; ss < n_seq; ss++)
              if (scs[ss] && scs[ss]->f)
                en += scs[ss]->f(1, jj, u - 1, u, VRNA_DECOMP_EXT_EXT_STEM, scs[ss]->data);
          }

          if (fij == en) {
            *i                            = u;
            *j                            = jj;
            *k                            = u - 1;
            bp_stack[++(*stack_count)].i  = u;
            bp_stack[(*stack_count)].j    = jj;
            return 1;
          }
        }
      }
      break;

    case 2:
      for (u = jj - turn - 1; u >= 1; u--) {
        if (with_gquad) {
          if (fij == my_f5[u - 1] + my_ggg[idx[jj] + u]) {
            *i  = *j = -1;
            *k  = u - 1;
            return vrna_BT_gquad_mfe(fc, u, jj, bp_stack, stack_count);
          }
        }

        if (evaluate(1, jj, u - 1, u, VRNA_DECOMP_EXT_EXT_STEM, &hc_dat_local)) {
          en = my_c[idx[jj] + u] +
               my_f5[u - 1];

          for (ss = 0; ss < n_seq; ss++) {
            tt  = vrna_get_ptype_md(S[ss][u], S[ss][jj], md);
            mm5 = (a2s[ss][u] > 1) ? S5[ss][u] : -1;
            mm3 = (a2s[ss][jj] < a2s[ss][S[0][0]]) ? S3[ss][jj] : -1;      /* why S[0][0] ??? */
            en  += vrna_E_ext_stem(tt, mm5, mm3, P);
          }

          if (scs) {
            for (ss = 0; ss < n_seq; ss++)
              if (scs[ss] && scs[ss]->f)
                en += scs[ss]->f(1, jj, u - 1, u, VRNA_DECOMP_EXT_EXT_STEM, scs[ss]->data);
          }

          if (fij == en) {
            *i                            = u;
            *j                            = jj;
            *k                            = u - 1;
            bp_stack[++(*stack_count)].i  = u;
            bp_stack[(*stack_count)].j    = jj;
            return 1;
          }
        }
      }
      break;
  }

  return 0;
}


PRIVATE int
BT_ext_loop_f3(vrna_fold_compound_t *fc,
               int                  *k,
               int                  maxdist,
               int                  *i,
               int                  *j,
               vrna_bp_stack_t      *bp_stack,
               int                  *stack_count)
{
  char                      **ptype;
  short                     mm5, mm3, *S1;
  unsigned int              type;
  int                       length, fij, fj, ii, u, *f3, **c, **ggg,
                            dangle_model, turn, with_gquad, en;
  vrna_param_t              *P;
  vrna_md_t                 *md;
  vrna_sc_t                 *sc;
  vrna_callback_hc_evaluate *evaluate;
  struct default_data       hc_dat_local;

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
  turn          = md->min_loop_size;
  with_gquad    = md->gquad;
  evaluate      = prepare_hc_default_window(fc, &hc_dat_local);

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

  if (ii > maxdist - turn + 1) {
    /* no more pairs */
    *i  = *j = -1;
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
      for (u = maxdist; u > ii + turn; u--) {
        if (with_gquad) {
          if (fij == ggg[ii][u - ii] + f3[u + 1]) {
            *i  = *j = -1;
            *k  = u + 1;
            return vrna_BT_gquad_mfe(fc, ii, u, bp_stack, stack_count);
          }
        }

        if (evaluate(ii, length, u, u + 1, VRNA_DECOMP_EXT_STEM_EXT, &hc_dat_local)) {
          type = vrna_get_ptype_window(ii, u, ptype);

          en = c[ii][u - ii] +
               vrna_E_ext_stem(type, -1, -1, P) +
               f3[u + 1];

          if ((sc) && (sc->f))
            en += sc->f(ii, length, u, u + 1, VRNA_DECOMP_EXT_STEM_EXT, sc->data);

          if (fij == en) {
            *i                            = ii;
            *j                            = u;
            *k                            = u + 1;
            bp_stack[++(*stack_count)].i  = ii;
            bp_stack[(*stack_count)].j    = u;
            return 1;
          }
        }
      }
      break;

    /* dangles on both sides */
    case 2:
      mm5 = (ii > 1) ? S1[ii - 1] : -1;
      for (u = maxdist; u > ii + turn; u--) {
        if (with_gquad) {
          if (fij == ggg[ii][u - ii] + f3[u + 1]) {
            *i  = *j = -1;
            *k  = u + 1;
            return vrna_BT_gquad_mfe(fc, ii, u, bp_stack, stack_count);
          }
        }

        if (evaluate(ii, length, u, u + 1, VRNA_DECOMP_EXT_STEM_EXT, &hc_dat_local)) {
          mm3   = (u < length) ? S1[u + 1] : -1;
          type  = vrna_get_ptype_window(ii, u, ptype);

          en = c[ii][u - ii] +
               vrna_E_ext_stem(type, mm5, mm3, P) +
               f3[u + 1];

          if ((sc) && (sc->f))
            en += sc->f(ii, length, u, u + 1, VRNA_DECOMP_EXT_STEM_EXT, sc->data);

          if (fij == en) {
            *i                            = ii;
            *j                            = u;
            *k                            = u + 1;
            bp_stack[++(*stack_count)].i  = ii;
            bp_stack[(*stack_count)].j    = u;
            return 1;
          }
        }
      }
      break;

    default:
      mm5 = S1[ii];
      for (u = maxdist; u > ii + turn; u--) {
        if (with_gquad) {
          if (fij == ggg[ii][u - ii] + f3[u + 1]) {
            *i  = *j = -1;
            *k  = u + 1;
            return vrna_BT_gquad_mfe(fc, ii, u, bp_stack, stack_count);
          }
        }

        if (u + 2 <= length) {
          if (evaluate(ii, length, u, u + 2, VRNA_DECOMP_EXT_STEM_EXT1, &hc_dat_local)) {
            mm3   = S1[u + 1];
            type  = vrna_get_ptype_window(ii + 1, u, ptype);

            en = c[ii + 1][u - ii - 1] + vrna_E_ext_stem(type, mm5, mm3, P) + f3[u + 2];

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
              bp_stack[++(*stack_count)].i  = ii + 1;
              bp_stack[(*stack_count)].j    = u;
              return 1;
            }
          }
        } else {
          if (evaluate(ii, length, ii + 1, u, VRNA_DECOMP_EXT_STEM, &hc_dat_local)) {
            mm3   = (u < length) ? S1[u + 1] : -1;
            type  = vrna_get_ptype_window(ii + 1, u, ptype);

            en = c[ii + 1][u - ii - 1] +
                 vrna_E_ext_stem(type, mm5, mm3, P);

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
              bp_stack[++(*stack_count)].i  = ii + 1;
              bp_stack[(*stack_count)].j    = u;
              return 1;
            }
          }
        }

        if (evaluate(ii, length, u, u + 1, VRNA_DECOMP_EXT_STEM_EXT1, &hc_dat_local)) {
          type = vrna_get_ptype_window(ii + 1, u, ptype);

          en = c[ii + 1][u - ii - 1] +
               vrna_E_ext_stem(type, mm5, -1, P) +
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
            bp_stack[++(*stack_count)].i  = ii + 1;
            bp_stack[(*stack_count)].j    = u;
            return 1;
          }
        }

        if (u + 2 <= length) {
          if (evaluate(ii, length, u, u + 2, VRNA_DECOMP_EXT_STEM_EXT, &hc_dat_local)) {
            mm3   = S1[u + 1];
            type  = vrna_get_ptype_window(ii, u, ptype);

            en = c[ii][u - ii] +
                 vrna_E_ext_stem(type, -1, mm3, P) +
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
              bp_stack[++(*stack_count)].i  = ii;
              bp_stack[(*stack_count)].j    = u;
              return 1;
            }
          }
        } else {
          if (evaluate(ii, length, ii, u, VRNA_DECOMP_EXT_STEM, &hc_dat_local)) {
            mm3   = (u < length) ? S1[u + 1] : -1;
            type  = vrna_get_ptype_window(ii, u, ptype);

            en = c[ii][u - ii] +
                 vrna_E_ext_stem(type, -1, mm3, P);

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
              bp_stack[++(*stack_count)].i  = ii;
              bp_stack[(*stack_count)].j    = u;
              return 1;
            }
          }
        }

        if (evaluate(ii, length, u, u + 1, VRNA_DECOMP_EXT_STEM_EXT, &hc_dat_local)) {
          type = vrna_get_ptype_window(ii, u, ptype);

          en = c[ii][u - ii] +
               vrna_E_ext_stem(type, -1, -1, P) + f3[u + 1];

          if (sc)
            if (sc->f)
              en += sc->f(ii, length, u, u + 1, VRNA_DECOMP_EXT_STEM_EXT, sc->data);

          if (fij == en) {
            *i                            = ii;
            *j                            = u;
            *k                            = u + 1;
            bp_stack[++(*stack_count)].i  = ii;
            bp_stack[(*stack_count)].j    = u;
            return 1;
          }
        }
      }

      break;
  }

  return 0;
}


PRIVATE int
BT_ext_loop_f3_comparative(vrna_fold_compound_t *fc,
                           int                  *k,
                           int                  maxdist,
                           int                  *i,
                           int                  *j,
                           vrna_bp_stack_t      *bp_stack,
                           int                  *stack_count)
{
  short                     **S, **S5, **S3;
  unsigned int              type, ss, n_seq, **a2s;
  int                       fij, cc, fj, ii, u, *f3, **c, **ggg,
                            dangle_model, turn, with_gquad;
  vrna_param_t              *P;
  vrna_md_t                 *md;
  vrna_sc_t                 **scs;
  vrna_callback_hc_evaluate *evaluate;
  struct default_data       hc_dat_local;

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
  turn          = md->min_loop_size;
  with_gquad    = md->gquad;
  evaluate      = prepare_hc_default_window(fc, &hc_dat_local);

  ii = *k;

  /* nibble off unpaired 5' bases */
  do {
    fij = f3[ii];
    fj  = INF;

    if (evaluate(ii, maxdist, ii + 1, maxdist, VRNA_DECOMP_EXT_EXT, &hc_dat_local)) {
      fj = f3[ii + 1];
      if (scs) {
        for (ss = 0; ss < n_seq; ss++)
          if (scs[ss]) {
            if (scs[ss]->energy_up)
              fj += scs[ss]->energy_up[ii][1];

            if (scs[ss]->f)
              fj += scs[ss]->f(ii, maxdist, ii + 1, maxdist, VRNA_DECOMP_EXT_EXT, scs[ss]->data);
          }
      }
    }

    if (++ii > maxdist)
      break;
  } while (fij == fj);
  ii--;

  if (ii > maxdist - turn + 1) {
    /* no more pairs */
    *i  = *j = -1;
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
      for (u = maxdist; u > ii + turn; u--) {
        if (with_gquad) {
          if (fij == ggg[ii][u - ii] + f3[u + 1]) {
            *i  = *j = -1;
            *k  = u + 1;
            return vrna_BT_gquad_mfe(fc, ii, u, bp_stack, stack_count);
          }
        }

        if (evaluate(ii, maxdist, u, u + 1, VRNA_DECOMP_EXT_STEM_EXT, &hc_dat_local)) {
          cc = c[ii][u - ii];
          for (ss = 0; ss < n_seq; ss++) {
            type  = vrna_get_ptype_md(S[ss][ii], S[ss][u], md);
            cc    += vrna_E_ext_stem(type, -1, -1, P);
          }

          if (fij == cc + f3[u + 1]) {
            *i                            = ii;
            *j                            = u;
            *k                            = u + 1;
            bp_stack[++(*stack_count)].i  = ii;
            bp_stack[(*stack_count)].j    = u;
            return 1;
          }
        }
      }
      break;

    /* dangles on both sides */
    case 2:
      for (u = maxdist; u > ii + turn; u--) {
        if (with_gquad) {
          if (fij == ggg[ii][u - ii] + f3[u + 1]) {
            *i  = *j = -1;
            *k  = u + 1;
            return vrna_BT_gquad_mfe(fc, ii, u, bp_stack, stack_count);
          }
        }

        if (evaluate(ii, maxdist, u, u + 1, VRNA_DECOMP_EXT_STEM_EXT, &hc_dat_local)) {
          cc = c[ii][u - ii];
          for (ss = 0; ss < n_seq; ss++) {
            type  = vrna_get_ptype_md(S[ss][ii], S[ss][u], md);
            cc    +=
              vrna_E_ext_stem(type, (a2s[ss][ii] > 1) ? S5[ss][ii] : -1,
                              (a2s[ss][u] < a2s[ss][S[0][0]]) ? S3[ss][u] : -1, P);
          }

          if (fij == cc + f3[u + 1]) {
            *i                            = ii;
            *j                            = u;
            *k                            = u + 1;
            bp_stack[++(*stack_count)].i  = ii;
            bp_stack[(*stack_count)].j    = u;
            return 1;
          }
        }
      }
      break;
  }

  return 0;
}


PRIVATE int
BT_ext_loop_f3_pp(vrna_fold_compound_t  *fc,
                  int                   *i,
                  int                   maxj)
{
  int j, start;

  j     = -1;
  start = *i;

  if (fc) {
    char                      **ptype;
    short                     *S1;
    unsigned int              type;
    int                       traced2, length, turn, dangle_model, with_gquad, maxdist, cc,
                              **c, **ggg, *f3, fij, ii;
    vrna_param_t              *P;
    vrna_md_t                 *md;
    vrna_sc_t                 *sc;
    vrna_callback_hc_evaluate *evaluate;
    struct default_data       hc_dat_local;

    length        = fc->length;
    S1            = fc->sequence_encoding;
    ptype         = fc->ptype_local;
    f3            = fc->matrices->f3_local;
    c             = fc->matrices->c_local;
    ggg           = fc->matrices->ggg_local;
    sc            = fc->sc;
    P             = fc->params;
    md            = &(P->model_details);
    turn          = md->min_loop_size;
    dangle_model  = md->dangles;
    with_gquad    = md->gquad;
    maxdist       = MIN2(fc->window_size, maxj);
    traced2       = 0;
    ii            = start;
    evaluate      = prepare_hc_default_window(fc, &hc_dat_local);

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
        for (j = start + turn; j <= ii + maxdist; j++) {
          if (evaluate(start, length, j, j + 1, VRNA_DECOMP_EXT_STEM_EXT, &hc_dat_local)) {
            type = vrna_get_ptype_window(start, j, ptype);

            cc = c[start][j - start] +
                 vrna_E_ext_stem(type, -1, -1, P) +
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
        for (j = start + turn; j <= ii + maxdist; j++) {
          if (evaluate(start, length, j, j + 1, VRNA_DECOMP_EXT_STEM_EXT, &hc_dat_local)) {
            type = vrna_get_ptype_window(start, j, ptype);

            cc = c[start][j - start] +
                 vrna_E_ext_stem(type, (start > 1) ? S1[start - 1] : -1,
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
        for (j = start + turn; j <= ii + maxdist; j++) {
          if (evaluate(start, length, j, j + 1, VRNA_DECOMP_EXT_STEM_EXT, &hc_dat_local)) {
            type = vrna_get_ptype_window(start, j, ptype);

            cc = c[start][j - start] +
                 vrna_E_ext_stem(type, -1, -1, P) +
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
                      vrna_E_ext_stem(type, -1, S1[j + 1], P) +
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
                 vrna_E_ext_stem(type, S1[start], -1, P) +
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
                   vrna_E_ext_stem(type, S1[start], S1[j + 1], P) +
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
BT_ext_loop_f3_pp_comparative(vrna_fold_compound_t  *fc,
                              int                   *i,
                              int                   maxj)
{
  int j, start;

  j = -1;

  if (fc) {
    short                     **S, **S5, **S3;
    unsigned int              tt, s, n_seq, **a2s;
    int                       traced2, length, turn, dangle_model, with_gquad, maxdist, cc, **c,
                              **ggg, *f3, fij;
    vrna_param_t              *P;
    vrna_md_t                 *md;
    vrna_sc_t                 **scs;
    vrna_callback_hc_evaluate *evaluate;
    struct default_data       hc_dat_local;

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
    turn          = md->min_loop_size;
    dangle_model  = md->dangles;
    with_gquad    = md->gquad;
    maxdist       = MIN2(fc->window_size, maxj);
    traced2       = 0;
    start         = *i;
    evaluate      = prepare_hc_default_window(fc, &hc_dat_local);

    fij = f3[start];

    /* try to nibble off unpaired 5' bases */
    if ((scs) && (evaluate(start, length, start + 1, length, VRNA_DECOMP_EXT_EXT, &hc_dat_local))) {
      cc = f3[start + 1];

      for (s = 0; s < n_seq; s++)
        if (scs[s]) {
          if (scs[s]->energy_up)
            cc += scs[s]->energy_up[start][1];

          if (scs[s]->f)
            cc +=
              scs[s]->f(start, length, start + 1, length, VRNA_DECOMP_EXT_EXT, scs[s]->data);
        }

      if (fij == cc)
        /* simple 5' unpaired extensions, so we skip this hit */
        return 0;
    }

    /* get pairing partner j */
    switch (dangle_model) {
      case 0:
        for (j = start + turn; j <= MIN2(start + maxdist, length - 1); j++) {
          if (evaluate(start, length, j, j + 1, VRNA_DECOMP_EXT_STEM_EXT, &hc_dat_local)) {
            cc = c[start][j - start];

            for (s = 0; s < n_seq; s++) {
              tt  = vrna_get_ptype_md(S[s][start], S[s][j], md);
              cc  += vrna_E_ext_stem(tt, -1, -1, P);
            }

            if (fij == cc + f3[j + 1]) {
              traced2 = 1;
              break;
            }
          }

          if (with_gquad) {
            cc = ggg[start][j - start];

            if (fij == cc + f3[j + 1]) {
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
              cc  += vrna_E_ext_stem(tt, -1, -1, P);
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
        for (j = start + turn; j <= MIN2(start + maxdist, length - 1); j++) {
          if (evaluate(start, length, j, j + 1, VRNA_DECOMP_EXT_STEM_EXT, &hc_dat_local)) {
            cc = c[start][j - start];

            for (s = 0; s < n_seq; s++) {
              tt  = vrna_get_ptype_md(S[s][start], S[s][j], md);
              cc  +=
                vrna_E_ext_stem(tt,
                                (a2s[s][start] > 1) ? S5[s][start] : -1,
                                (a2s[s][j] < a2s[s][S[0][0]]) ? S3[s][j] : -1,
                                P);
            }

            if (fij == cc + f3[j + 1]) {
              traced2 = 1;
              break;
            }
          }

          if (with_gquad) {
            cc = ggg[start][j - start];

            if (fij == cc + f3[j + 1]) {
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
              cc  += vrna_E_ext_stem(tt, (a2s[s][start] > 1) ? S5[s][start] :  -1, -1, P);
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
