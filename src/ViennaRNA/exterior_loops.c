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
#include "ViennaRNA/energy_par.h"
#include "ViennaRNA/utils.h"
#include "ViennaRNA/constraints.h"
#include "ViennaRNA/gquad.h"
#include "ViennaRNA/structured_domains.h"
#include "ViennaRNA/unstructured_domains.h"
#include "ViennaRNA/exterior_loops.h"

#ifdef __GNUC__
# define INLINE inline
#else
# define INLINE
#endif

#include "exterior_loops.inc"

/*
 #################################
 # PRIVATE FUNCTION DECLARATIONS #
 #################################
 */

PRIVATE int
E_ext_loop_5(vrna_fold_compound_t *vc);


PRIVATE int
E_ext_loop_5_comparative(vrna_fold_compound_t *vc);


PRIVATE int
E_ext_loop_3(vrna_fold_compound_t *fc,
             int                  i);


PRIVATE int
E_ext_loop_3_comparative(vrna_fold_compound_t *fc,
                         int                  i);


/*
 #################################
 # BEGIN OF FUNCTION DEFINITIONS #
 #################################
 */
PUBLIC int
vrna_E_ext_loop_5(vrna_fold_compound_t *fc)
{
  if (fc) {
    switch (fc->type) {
      case VRNA_FC_TYPE_SINGLE:
        return E_ext_loop_5(fc);
        break;

      case VRNA_FC_TYPE_COMPARATIVE:
        return E_ext_loop_5_comparative(fc);
        break;
    }
  }

  return 0;
}


PUBLIC int
vrna_E_ext_loop_3(vrna_fold_compound_t  *fc,
                  int                   i)
{
  if (fc) {
    switch (fc->type) {
      case VRNA_FC_TYPE_SINGLE:
        return E_ext_loop_3(fc, i);
        break;

      case VRNA_FC_TYPE_COMPARATIVE:
        return E_ext_loop_3_comparative(fc, i);
        break;
    }
  }

  return 0;
}


PUBLIC int
E_Stem(int          type,
       int          si1,
       int          sj1,
       int          extLoop,
       vrna_param_t *P)
{
  int energy  = 0;
  int d5      = (si1 >= 0) ? P->dangle5[type][si1] : 0;
  int d3      = (sj1 >= 0) ? P->dangle3[type][sj1] : 0;

  if (type > 2)
    energy += P->TerminalAU;

  if (si1 >= 0 && sj1 >= 0)
    energy += (extLoop) ? P->mismatchExt[type][si1][sj1] : P->mismatchM[type][si1][sj1];
  else
    energy += d5 + d3;

  if (!extLoop)
    energy += P->MLintern[type];

  return energy;
}


PUBLIC int
E_ExtLoop(int           type,
          int           si1,
          int           sj1,
          vrna_param_t  *P)
{
  int energy = 0;

  if (si1 >= 0 && sj1 >= 0)
    energy += P->mismatchExt[type][si1][sj1];
  else if (si1 >= 0)
    energy += P->dangle5[type][si1];
  else if (sj1 >= 0)
    energy += P->dangle3[type][sj1];

  if (type > 2)
    energy += P->TerminalAU;

  return energy;
}


PUBLIC int
E_ext_loop(int                  i,
           int                  j,
           vrna_fold_compound_t *vc)
{
  char                      *ptype;
  unsigned char             *hard_constraints;
  short                     *S;
  int                       ij, en, e, type, cp, *idx;
  vrna_param_t              *P;
  vrna_md_t                 *md;
  vrna_sc_t                 *sc;
  vrna_callback_hc_evaluate *evaluate;
  struct default_data       hc_dat_local;

  cp                = vc->cutpoint;
  S                 = vc->sequence_encoding;
  idx               = vc->jindx;
  ptype             = vc->ptype;
  P                 = vc->params;
  md                = &(P->model_details);
  hard_constraints  = vc->hc->matrix;
  sc                = vc->sc;

  hc_dat_local.idx  = idx;
  hc_dat_local.mx   = hard_constraints;
  hc_dat_local.sn   = vc->strand_number;

  if (vc->hc->f) {
    evaluate            = &hc_default_user;
    hc_dat_local.hc_f   = vc->hc->f;
    hc_dat_local.hc_dat = vc->hc->data;
  } else {
    evaluate = &hc_default;
  }

  e     = INF;
  ij    = idx[j] + i;
  type  = get_pair_type(ij, ptype);

  if ((cp < 0) || (((i) >= cp) || ((j) < cp))) {
    /* regular exterior loop */
    if (evaluate(i, j, i, j, VRNA_DECOMP_EXT_STEM, &hc_dat_local)) {
      switch (md->dangles) {
        case 2:
          e = E_ExtLoop(type, S[i - 1], S[j + 1], P);
          break;

        case 0:
        /* fall through */

        default:
          e = E_ExtLoop(type, -1, -1, P);
          break;
      }
      if (sc)
        if (sc->f)
          e += sc->f(i, j, i, j, VRNA_DECOMP_EXT_STEM, sc->data);
    }

    if (md->dangles % 2) {
      ij = idx[j - 1] + i;
      if (evaluate(i, j, i, j - 1, VRNA_DECOMP_EXT_STEM, &hc_dat_local)) {
        type = get_pair_type(ij, ptype);

        en = E_ExtLoop(type, -1, S[j], P);

        if (sc)
          if (sc->f)
            en += sc->f(i, j, i, j - 1, VRNA_DECOMP_EXT_STEM, sc->data);

        e = MIN2(e, en);
      }

      ij = idx[j] + i + 1;
      if (evaluate(i, j, i + 1, j, VRNA_DECOMP_EXT_STEM, &hc_dat_local)) {
        type = get_pair_type(ij, ptype);

        en = E_ExtLoop(type, S[i], -1, P);

        if (sc)
          if (sc->f)
            en += sc->f(i, j, i + 1, j, VRNA_DECOMP_EXT_STEM, sc->data);

        e = MIN2(e, en);
      }
    }
  }

  return e;
}


PRIVATE int
E_ext_loop_5(vrna_fold_compound_t *vc)
{
  char                      *ptype;
  unsigned char             *hc;
  short                     *S;
  int                       en, i, j, ij, type, length, *indx, *hc_up, *f5, *c, dangle_model,
                            *ggg, with_gquad, turn, k, u, with_ud;
  vrna_sc_t                 *sc;
  vrna_param_t              *P;
  vrna_ud_t                 *domains_up;
  vrna_callback_hc_evaluate *evaluate;
  struct default_data       hc_dat_local;

  length        = (int)vc->length;
  ptype         = vc->ptype;
  S             = vc->sequence_encoding;
  indx          = vc->jindx;
  hc            = vc->hc->matrix;
  hc_up         = vc->hc->up_ext;
  sc            = vc->sc;
  f5            = vc->matrices->f5;
  c             = vc->matrices->c;
  P             = vc->params;
  dangle_model  = P->model_details.dangles;
  ggg           = vc->matrices->ggg;
  with_gquad    = P->model_details.gquad;
  turn          = P->model_details.min_loop_size;
  domains_up    = vc->domains_up;
  with_ud       = (domains_up && domains_up->energy_cb) ? 1 : 0;

  hc_dat_local.idx    = indx;
  hc_dat_local.mx     = hc;
  hc_dat_local.hc_up  = hc_up;
  hc_dat_local.sn     = vc->strand_number;

  if (vc->hc->f) {
    evaluate            = &hc_default_user;
    hc_dat_local.hc_f   = vc->hc->f;
    hc_dat_local.hc_dat = vc->hc->data;
  } else {
    evaluate = &hc_default;
  }

  f5[0] = 0;
  for (i = 1; i <= turn + 1; i++) {
    if (f5[i - 1] != INF) {
      if (evaluate(1, i, 1, i - 1, VRNA_DECOMP_EXT_EXT, &hc_dat_local)) {
        f5[i] = f5[i - 1];

        if (sc) {
          if (sc->energy_up)
            f5[i] += sc->energy_up[i][1];

          if (sc->f)
            f5[i] += sc->f(1, i, 1, i - 1, VRNA_DECOMP_EXT_EXT, sc->data);
        }
      } else {
        f5[i] = INF;
      }
    } else {
      f5[i] = INF;
    }
  }

  if (with_ud) {
    /* do we include ligand binding? */
    /*  construct all possible combinations of
     *  f[i-1] + L[i,j] with j <= turn + 1
     */
    for (i = 1; i <= turn + 1; i++) {
      if (f5[i - 1] != INF) {
        for (k = 0; k < domains_up->uniq_motif_count; k++) {
          u = domains_up->uniq_motif_size[k];
          j = i + u - 1;
          if (j <= turn + 1) {
            if (evaluate(1, j, i - 1, i, VRNA_DECOMP_EXT_EXT, &hc_dat_local)) {
              en = f5[i - 1] +
                   domains_up->energy_cb(vc,
                                         i, j,
                                         VRNA_UNSTRUCTURED_DOMAIN_EXT_LOOP,
                                         domains_up->data);

              if (sc) {
                if (sc->energy_up)
                  en += sc->energy_up[i][u];

                if (sc->f)
                  en += sc->f(1, j, 1, j - u, VRNA_DECOMP_EXT_EXT, sc->data);
              }

              f5[j] = MIN2(f5[j], en);
            }
          }
        }
      }
    }
  }

  /* duplicated code may be faster than conditions inside loop ;) */
  switch (dangle_model) {
    /* dont use dangling end and mismatch contributions at all */
    case 0:
      for (j = turn + 2; j <= length; j++) {
        /* initialize with INF */
        f5[j] = INF;

        /* check for 3' extension with one unpaired nucleotide */
        if (f5[j - 1] != INF) {
          if (evaluate(1, j, 1, j - 1, VRNA_DECOMP_EXT_EXT, &hc_dat_local)) {
            f5[j] = f5[j - 1];

            if (sc) {
              if (sc->energy_up)
                f5[j] += sc->energy_up[j][1];

              if (sc->f)
                f5[j] += sc->f(1, j, 1, j - 1, VRNA_DECOMP_EXT_EXT, sc->data);
            }
          }
        }

        if (with_ud) {
          for (k = 0; k < domains_up->uniq_motif_count; k++) {
            u = domains_up->uniq_motif_size[k];
            if ((j - u >= 0) && (f5[j - u] != INF)) {
              if (evaluate(1, j, 1, j - u, VRNA_DECOMP_EXT_EXT, &hc_dat_local)) {
                en = f5[j - u] +
                     domains_up->energy_cb(vc,
                                           j - u + 1,
                                           j,
                                           VRNA_UNSTRUCTURED_DOMAIN_EXT_LOOP | VRNA_UNSTRUCTURED_DOMAIN_MOTIF,
                                           domains_up->data);

                if (sc) {
                  if (sc->energy_up)
                    en += sc->energy_up[j - u + 1][u];

                  if (sc->f)
                    en += sc->f(1, j, 1, j - u, VRNA_DECOMP_EXT_EXT, sc->data);
                }

                f5[j] = MIN2(f5[j], en);
              }
            }
          }
        }

        /* check for possible stems branching off the exterior loop */
        if (sc && sc->f) {
          for (i = j - turn - 1; i > 1; i--) {
            if (f5[i - 1] != INF) {
              ij = indx[j] + i;

              if (with_gquad)
                f5[j] = MIN2(f5[j], f5[i - 1] + ggg[ij]);

              if (c[ij] != INF) {
                if (evaluate(1, j, i - 1, i, VRNA_DECOMP_EXT_EXT_STEM, &hc_dat_local)) {
                  type = get_pair_type(ij, ptype);

                  en = f5[i - 1] +
                       c[ij] +
                       E_ExtLoop(type, -1, -1, P) +
                       sc->f(1, j, i - 1, i, VRNA_DECOMP_EXT_EXT_STEM, sc->data);

                  f5[j] = MIN2(f5[j], en);
                }
              }
            }
          }
        } else {
          for (i = j - turn - 1; i > 1; i--) {
            if (f5[i - 1] != INF) {
              ij = indx[j] + i;

              if (with_gquad)
                f5[j] = MIN2(f5[j], f5[i - 1] + ggg[ij]);

              if (c[ij] != INF) {
                if (evaluate(1, j, i - 1, i, VRNA_DECOMP_EXT_EXT_STEM, &hc_dat_local)) {
                  type = get_pair_type(ij, ptype);

                  en = f5[i - 1] +
                       c[ij] +
                       E_ExtLoop(type, -1, -1, P);

                  f5[j] = MIN2(f5[j], en);
                }
              }
            }
          }
        }

        ij = indx[j] + 1;

        if (with_gquad)
          f5[j] = MIN2(f5[j], ggg[ij]);

        if (c[ij] != INF) {
          if (evaluate(1, j, 1, j, VRNA_DECOMP_EXT_STEM, &hc_dat_local)) {
            type = get_pair_type(ij, ptype);

            en = c[ij] +
                 E_ExtLoop(type, -1, -1, P);

            if (sc)
              if (sc->f)
                en += sc->f(1, j, 1, j, VRNA_DECOMP_EXT_STEM, sc->data);

            f5[j] = MIN2(f5[j], en);
          }
        }
      }
      break;

    /* always use dangles on both sides */
    case 2:
      for (j = turn + 2; j < length; j++) {
        f5[j] = INF;

        if (f5[j - 1] != INF) {
          if (evaluate(1, j, 1, j - 1, VRNA_DECOMP_EXT_EXT, &hc_dat_local)) {
            f5[j] = f5[j - 1];

            if (sc) {
              if (sc->energy_up)
                f5[j] += sc->energy_up[j][1];

              if (sc->f)
                f5[j] += sc->f(1, j, 1, j - 1, VRNA_DECOMP_EXT_EXT, sc->data);
            }
          }
        }

        if (with_ud) {
          for (k = 0; k < domains_up->uniq_motif_count; k++) {
            u = domains_up->uniq_motif_size[k];
            if ((j - u >= 0) && (f5[j - u] != INF)) {
              if (evaluate(1, j, 1, j - u, VRNA_DECOMP_EXT_EXT, &hc_dat_local)) {
                en = f5[j - u] +
                     domains_up->energy_cb(vc,
                                           j - u + 1,
                                           j,
                                           VRNA_UNSTRUCTURED_DOMAIN_EXT_LOOP | VRNA_UNSTRUCTURED_DOMAIN_MOTIF,
                                           domains_up->data);

                if (sc) {
                  if (sc->energy_up)
                    en += sc->energy_up[j - u + 1][u];

                  if (sc->f)
                    en += sc->f(1, j, 1, j - u, VRNA_DECOMP_EXT_EXT, sc->data);
                }

                f5[j] = MIN2(f5[j], en);
              }
            }
          }
        }

        if (sc && sc->f) {
          for (i = j - turn - 1; i > 1; i--) {
            if (f5[i - 1] != INF) {
              ij = indx[j] + i;

              if (with_gquad)
                f5[j] = MIN2(f5[j], f5[i - 1] + ggg[ij]);

              if (c[ij] != INF) {
                if (evaluate(1, j, i - 1, i, VRNA_DECOMP_EXT_EXT_STEM, &hc_dat_local)) {
                  type = get_pair_type(ij, ptype);

                  en = f5[i - 1] +
                       c[ij] +
                       E_ExtLoop(type, S[i - 1], S[j + 1], P) +
                       sc->f(1, j, i - 1, i, VRNA_DECOMP_EXT_EXT_STEM, sc->data);

                  f5[j] = MIN2(f5[j], en);
                }
              }
            }
          }
        } else {
          for (i = j - turn - 1; i > 1; i--) {
            if (f5[i - 1] != INF) {
              ij = indx[j] + i;

              if (with_gquad)
                f5[j] = MIN2(f5[j], f5[i - 1] + ggg[ij]);

              if (c[ij] != INF) {
                if (evaluate(1, j, i - 1, i, VRNA_DECOMP_EXT_EXT_STEM, &hc_dat_local)) {
                  type = get_pair_type(ij, ptype);

                  en = f5[i - 1] +
                       c[ij] +
                       E_ExtLoop(type, S[i - 1], S[j + 1], P);

                  f5[j] = MIN2(f5[j], en);
                }
              }
            }
          }
        }

        ij = indx[j] + 1;

        if (with_gquad)
          f5[j] = MIN2(f5[j], ggg[ij]);

        if (c[ij] != INF) {
          if (evaluate(1, j, 1, j, VRNA_DECOMP_EXT_STEM, &hc_dat_local)) {
            type = get_pair_type(ij, ptype);

            en = c[ij] +
                 E_ExtLoop(type, -1, S[j + 1], P);

            if (sc)
              if (sc->f)
                en += sc->f(1, j, 1, j, VRNA_DECOMP_EXT_STEM, sc->data);

            f5[j] = MIN2(f5[j], en);
          }
        }
      }

      f5[length] = INF;
      if (f5[length - 1] != INF) {
        if (evaluate(1, length, 1, length - 1, VRNA_DECOMP_EXT_EXT, &hc_dat_local)) {
          f5[length] = f5[length - 1];
          if (sc) {
            if (sc->energy_up)
              f5[length] += sc->energy_up[length][1];

            if (sc->f)
              f5[length] += sc->f(1, length, 1, length - 1, VRNA_DECOMP_EXT_EXT, sc->data);
          }
        }
      }

      if (with_ud) {
        for (k = 0; k < domains_up->uniq_motif_count; k++) {
          u = domains_up->uniq_motif_size[k];
          if ((length - u >= 0) && (f5[length - u] != INF)) {
            if (evaluate(1, length, 1, length - u, VRNA_DECOMP_EXT_EXT, &hc_dat_local)) {
              en = f5[length - u] +
                   domains_up->energy_cb(vc,
                                         length - u + 1,
                                         length,
                                         VRNA_UNSTRUCTURED_DOMAIN_EXT_LOOP | VRNA_UNSTRUCTURED_DOMAIN_MOTIF,
                                         domains_up->data);

              if (sc) {
                if (sc->energy_up)
                  en += sc->energy_up[length - u + 1][u];

                if (sc->f)
                  en += sc->f(1, length, 1, length - u, VRNA_DECOMP_EXT_EXT, sc->data);
              }

              f5[length] = MIN2(f5[length], en);
            }
          }
        }
      }

      if (sc && sc->f) {
        for (i = length - turn - 1; i > 1; i--) {
          if (f5[i - 1] != INF) {
            ij = indx[length] + i;

            if (with_gquad)
              f5[length] = MIN2(f5[length], f5[i - 1] + ggg[ij]);

            if (c[ij] != INF) {
              if (evaluate(1, length, i - 1, i, VRNA_DECOMP_EXT_EXT_STEM, &hc_dat_local)) {
                type = get_pair_type(ij, ptype);

                en = f5[i - 1] +
                     c[ij] +
                     E_ExtLoop(type, S[i - 1], -1, P) +
                     sc->f(1, length, i - 1, i, VRNA_DECOMP_EXT_EXT_STEM, sc->data);

                f5[length] = MIN2(f5[length], en);
              }
            }
          }
        }
      } else {
        for (i = length - turn - 1; i > 1; i--) {
          if (f5[i - 1] != INF) {
            ij = indx[length] + i;

            if (with_gquad)
              f5[length] = MIN2(f5[length], f5[i - 1] + ggg[ij]);

            if (c[ij] != INF) {
              if (evaluate(1, length, i - 1, i, VRNA_DECOMP_EXT_EXT_STEM, &hc_dat_local)) {
                type = get_pair_type(ij, ptype);

                en = f5[i - 1] +
                     c[ij] +
                     E_ExtLoop(type, S[i - 1], -1, P);

                f5[length] = MIN2(f5[length], en);
              }
            }
          }
        }
      }

      ij = indx[length] + 1;

      if (with_gquad)
        f5[length] = MIN2(f5[length], ggg[ij]);

      if (c[ij] != INF) {
        if (evaluate(1, length, 1, length, VRNA_DECOMP_EXT_STEM, &hc_dat_local)) {
          type = get_pair_type(ij, ptype);

          en = c[ij] +
               E_ExtLoop(type, -1, -1, P);

          if (sc)
            if (sc->f)
              en += sc->f(1, length, 1, length, VRNA_DECOMP_EXT_STEM, sc->data);

          f5[length] = MIN2(f5[length], en);
        }
      }

      break;

    /* normal dangles, aka dangle_model = 1 || 3 */
    default:
      for (j = turn + 2; j <= length; j++) {
        f5[j] = INF;
        if (f5[j - 1] != INF) {
          if (evaluate(1, j, 1, j - 1, VRNA_DECOMP_EXT_EXT, &hc_dat_local)) {
            f5[j] = f5[j - 1];

            if (sc) {
              if (sc->energy_up)
                f5[j] += sc->energy_up[j][1];

              if (sc->f)
                f5[j] += sc->f(1, j, 1, j - 1, VRNA_DECOMP_EXT_EXT, sc->data);
            }
          }
        }

        if (with_ud) {
          for (k = 0; k < domains_up->uniq_motif_count; k++) {
            u = domains_up->uniq_motif_size[k];
            if ((j - u >= 0) && (f5[j - u] != INF)) {
              if (evaluate(1, j, 1, j - u, VRNA_DECOMP_EXT_EXT, &hc_dat_local)) {
                en = f5[j - u] +
                     domains_up->energy_cb(vc,
                                           j - u + 1,
                                           j,
                                           VRNA_UNSTRUCTURED_DOMAIN_EXT_LOOP | VRNA_UNSTRUCTURED_DOMAIN_MOTIF,
                                           domains_up->data);

                if (sc) {
                  if (sc->energy_up)
                    en += sc->energy_up[j - u + 1][u];

                  if (sc->f)
                    en += sc->f(1, j, 1, j - u, VRNA_DECOMP_EXT_EXT, sc->data);
                }

                f5[j] = MIN2(f5[j], en);
              }
            }
          }
        }

        for (i = j - turn - 1; i > 1; i--) {
          ij = indx[j] + i;
          if (f5[i - 1] != INF) {
            if (with_gquad)
              f5[j] = MIN2(f5[j], f5[i - 1] + ggg[ij]);

            if (c[ij] != INF) {
              if (evaluate(1, j, i - 1, i, VRNA_DECOMP_EXT_EXT_STEM, &hc_dat_local)) {
                type = get_pair_type(ij, ptype);

                en = f5[i - 1] +
                     c[ij] +
                     E_ExtLoop(type, -1, -1, P);

                if (sc)
                  if (sc->f)
                    en += sc->f(1, j, i - 1, i, VRNA_DECOMP_EXT_EXT_STEM, sc->data);

                f5[j] = MIN2(f5[j], en);
              }
            }
          }

          if ((f5[i - 2] != INF) && c[ij] != INF) {
            if (evaluate(1, j, i - 2, i, VRNA_DECOMP_EXT_EXT_STEM, &hc_dat_local)) {
              type = get_pair_type(ij, ptype);

              en = f5[i - 2] + c[ij] +
                   E_ExtLoop(type, S[i - 1], -1, P);

              if (sc) {
                if (sc->energy_up)
                  en += sc->energy_up[i - 1][1];

                if (sc->f)
                  en += sc->f(1, j, i - 2, i, VRNA_DECOMP_EXT_EXT_STEM, sc->data);
              }

              f5[j] = MIN2(f5[j], en);
            }
          }

          ij = indx[j - 1] + i;
          if (c[ij] != INF) {
            if (f5[i - 1] != INF) {
              if (evaluate(1, j, i - 1, i, VRNA_DECOMP_EXT_EXT_STEM1, &hc_dat_local)) {
                type = get_pair_type(ij, ptype);

                en = f5[i - 1] +
                     c[ij] +
                     E_ExtLoop(type, -1, S[j], P);

                if (sc) {
                  if (sc->energy_up)
                    en += sc->energy_up[j][1];

                  if (sc->f)
                    en += sc->f(1, j, i - 1, i, VRNA_DECOMP_EXT_EXT_STEM1, sc->data);
                }

                f5[j] = MIN2(f5[j], en);
              }
            }

            if (f5[i - 2] != INF) {
              if (evaluate(1, j, i - 2, i, VRNA_DECOMP_EXT_EXT_STEM1, &hc_dat_local)) {
                type = get_pair_type(ij, ptype);

                en = f5[i - 2] +
                     c[ij] +
                     E_ExtLoop(type, S[i - 1], S[j], P);

                if (sc) {
                  if (sc->energy_up)
                    en += sc->energy_up[i - 1][1] +
                          sc->energy_up[j][1];

                  if (sc->f)
                    en += sc->f(1, j, i - 2, i, VRNA_DECOMP_EXT_EXT_STEM1, sc->data);
                }

                f5[j] = MIN2(f5[j], en);
              }
            }
          }
        }

        ij = indx[j] + 1;

        if (with_gquad)
          f5[j] = MIN2(f5[j], ggg[ij]);

        if (c[ij] != INF) {
          if (evaluate(1, j, 1, j, VRNA_DECOMP_EXT_STEM, &hc_dat_local)) {
            type = get_pair_type(ij, ptype);

            en = c[ij] +
                 E_ExtLoop(type, -1, -1, P);

            if (sc)
              if (sc->f)
                en += sc->f(1, j, 1, j, VRNA_DECOMP_EXT_STEM, sc->data);

            f5[j] = MIN2(f5[j], en);
          }
        }

        ij = indx[j - 1] + 1;
        if (c[ij] != INF) {
          if (evaluate(1, j, 1, j - 1, VRNA_DECOMP_EXT_STEM, &hc_dat_local)) {
            type = get_pair_type(ij, ptype);

            en = c[ij] +
                 E_ExtLoop(type, -1, S[j], P);

            if (sc) {
              if (sc->energy_up)
                en += sc->energy_up[j][1];

              if (sc->f)
                en += sc->f(1, j, 1, j - 1, VRNA_DECOMP_EXT_STEM, sc->data);
            }

            f5[j] = MIN2(f5[j], en);
          }
        }
      }         /* end for j... */
      break;
  }

  return f5[length];
}


PRIVATE int
E_ext_loop_5_comparative(vrna_fold_compound_t *vc)
{
  unsigned char             *hc;
  unsigned int              **a2s;
  short                     **S, **S5, **S3;
  int                       en, i, j, ij, tt, length, *indx, *hc_up, *f5, *c, dangle_model,
                            *ggg, with_gquad, turn, n_seq, s, mm5, mm3;
  vrna_sc_t                 **scs;
  vrna_param_t              *P;
  vrna_md_t                 *md;
  vrna_callback_hc_evaluate *evaluate;
  struct default_data       hc_dat_local;

  length        = (int)vc->length;
  n_seq         = vc->n_seq;
  S             = vc->S;
  S5            = vc->S5;
  S3            = vc->S3;
  a2s           = vc->a2s;
  indx          = vc->jindx;
  hc            = vc->hc->matrix;
  hc_up         = vc->hc->up_ext;
  scs           = vc->scs;
  f5            = vc->matrices->f5;
  c             = vc->matrices->c;
  ggg           = vc->matrices->ggg;
  P             = vc->params;
  md            = &(P->model_details);
  dangle_model  = md->dangles;
  with_gquad    = md->gquad;
  turn          = md->min_loop_size;

  hc_dat_local.idx    = indx;
  hc_dat_local.mx     = hc;
  hc_dat_local.hc_up  = hc_up;
  hc_dat_local.sn     = vc->strand_number;

  if (vc->hc->f) {
    evaluate            = &hc_default_user;
    hc_dat_local.hc_f   = vc->hc->f;
    hc_dat_local.hc_dat = vc->hc->data;
  } else {
    evaluate = &hc_default;
  }

  f5[0] = 0;
  for (i = 1; i <= turn + 1; i++) {
    if (f5[i - 1] != INF) {
      if (evaluate(1, i, 1, i - 1, VRNA_DECOMP_EXT_EXT, &hc_dat_local)) {
        f5[i] = f5[i - 1];

        if (scs) {
          for (s = 0; s < n_seq; s++) {
            if (scs[s]) {
              if (scs[s]->energy_up)
                f5[i] += scs[s]->energy_up[a2s[s][j]][1];

              if (scs[s]->f)
                f5[i] += scs[s]->f(1, i, 1, i - 1, VRNA_DECOMP_EXT_EXT, scs[s]->data);
            }
          }
        }
      } else {
        f5[i] = INF;
      }
    } else {
      f5[i] = INF;
    }
  }

  /* duplicated code may be faster than conditions inside loop ;) */
  switch (dangle_model) {
    /* dont use dangling end and mismatch contributions at all */
    case 0:
      for (j = turn + 2; j <= length; j++) {
        /* initialize with INF */
        f5[j] = INF;

        /* check for 3' extension with one unpaired nucleotide */
        if (f5[j - 1] != INF) {
          if (evaluate(1, j, 1, j - 1, VRNA_DECOMP_EXT_EXT, &hc_dat_local)) {
            f5[j] = f5[j - 1];

            if (scs) {
              for (s = 0; s < n_seq; s++) {
                if (scs[s]) {
                  if (scs[s]->energy_up)
                    f5[j] += scs[s]->energy_up[a2s[s][j]][1];

                  if (scs[s]->f)
                    f5[j] += scs[s]->f(1, j, 1, j - 1, VRNA_DECOMP_EXT_EXT, scs[s]->data);
                }
              }
            }
          }
        }

        /* check for possible stems branching off the exterior loop */
        if (scs) {
          for (i = j - turn - 1; i > 1; i--) {
            if (f5[i - 1] != INF) {
              ij = indx[j] + i;

              if (with_gquad)
                f5[j] = MIN2(f5[j], f5[i - 1] + ggg[ij]);

              if (c[ij] != INF) {
                if (evaluate(1, j, i - 1, i, VRNA_DECOMP_EXT_EXT_STEM, &hc_dat_local)) {
                  en = f5[i - 1] + c[ij];
                  for (s = 0; s < n_seq; s++) {
                    tt  = get_pair_type_md(S[s][i], S[s][j], md);
                    en  += E_ExtLoop(tt, -1, -1, P);

                    if (scs[s] && scs[s]->f)
                      en += scs[s]->f(1, j, i - 1, i, VRNA_DECOMP_EXT_EXT_STEM, scs[s]->data);
                  }
                  f5[j] = MIN2(f5[j], en);
                }
              }
            }
          }
        } else {
          for (i = j - turn - 1; i > 1; i--) {
            if (f5[i - 1] != INF) {
              ij = indx[j] + i;

              if (with_gquad)
                f5[j] = MIN2(f5[j], f5[i - 1] + ggg[ij]);

              if (c[ij] != INF) {
                if (evaluate(1, j, i - 1, i, VRNA_DECOMP_EXT_EXT_STEM, &hc_dat_local)) {
                  en = f5[i - 1] +
                       c[ij];

                  for (s = 0; s < n_seq; s++) {
                    tt  = get_pair_type_md(S[s][i], S[s][j], md);
                    en  += E_ExtLoop(tt, -1, -1, P);
                  }

                  f5[j] = MIN2(f5[j], en);
                }
              }
            }
          }
        }

        ij = indx[j] + 1;

        if (with_gquad)
          f5[j] = MIN2(f5[j], ggg[ij]);

        if (c[ij] != INF) {
          if (evaluate(1, j, 1, j, VRNA_DECOMP_EXT_STEM, &hc_dat_local)) {
            en = c[ij];
            if (scs) {
              for (s = 0; s < n_seq; s++) {
                tt  = get_pair_type_md(S[s][i], S[s][j], md);
                en  += E_ExtLoop(tt, -1, -1, P);

                if (scs[s] && scs[s]->f)
                  en += scs[s]->f(1, j, 1, j, VRNA_DECOMP_EXT_STEM, scs[s]->data);
              }
            } else {
              for (s = 0; s < n_seq; s++) {
                tt  = get_pair_type_md(S[s][i], S[s][j], md);
                en  += E_ExtLoop(tt, -1, -1, P);
              }
            }

            f5[j] = MIN2(f5[j], en);
          }
        }
      }
      break;

    /* always use dangles on both sides */
    case 2:
      for (j = turn + 2; j < length; j++) {
        f5[j] = INF;

        if (f5[j - 1] != INF) {
          if (evaluate(1, j, 1, j - 1, VRNA_DECOMP_EXT_EXT, &hc_dat_local)) {
            f5[j] = f5[j - 1];

            if (scs) {
              for (s = 0; s < n_seq; s++) {
                if (scs[s]) {
                  if (scs[s]->energy_up)
                    f5[j] += scs[s]->energy_up[a2s[s][j]][1];

                  if (scs[s]->f)
                    f5[j] += scs[s]->f(1, j, 1, j - 1, VRNA_DECOMP_EXT_EXT, scs[s]->data);
                }
              }
            }
          }
        }

        if (scs) {
          for (i = j - turn - 1; i > 1; i--) {
            if (f5[i - 1] != INF) {
              ij = indx[j] + i;

              if (with_gquad)
                f5[j] = MIN2(f5[j], f5[i - 1] + ggg[ij]);

              if (c[ij] != INF) {
                if (evaluate(1, j, i - 1, i, VRNA_DECOMP_EXT_EXT_STEM, &hc_dat_local)) {
                  en = f5[i - 1] +
                       c[ij];

                  for (s = 0; s < n_seq; s++) {
                    tt  = get_pair_type_md(S[s][i], S[s][j], md);
                    mm5 = (a2s[s][i] > 1) ? S5[s][i] : -1;
                    mm3 = (a2s[s][j] < a2s[s][S[0][0]]) ? S3[s][j] : -1;      /* why S[0][0] ??? */
                    en  += E_ExtLoop(tt, mm5, mm3, P);

                    if (scs[s] && scs[s]->f)
                      en += scs[s]->f(1, j, i - 1, i, VRNA_DECOMP_EXT_EXT_STEM, scs[s]->data);
                  }
                  f5[j] = MIN2(f5[j], en);
                }
              }
            }
          }
        } else {
          for (i = j - turn - 1; i > 1; i--) {
            if (f5[i - 1] != INF) {
              ij = indx[j] + i;

              if (with_gquad)
                f5[j] = MIN2(f5[j], f5[i - 1] + ggg[ij]);

              if (c[ij] != INF) {
                if (evaluate(1, j, i - 1, i, VRNA_DECOMP_EXT_EXT_STEM, &hc_dat_local)) {
                  en = f5[i - 1] +
                       c[ij];

                  for (s = 0; s < n_seq; s++) {
                    tt  = get_pair_type_md(S[s][i], S[s][j], md);
                    mm5 = (a2s[s][i] > 1) ? S5[s][i] : -1;
                    mm3 = (a2s[s][j] < a2s[s][S[0][0]]) ? S3[s][j] : -1;      /* why S[0][0] ??? */
                    en  += E_ExtLoop(tt, mm5, mm3, P);
                  }
                  f5[j] = MIN2(f5[j], en);
                }
              }
            }
          }
        }

        ij = indx[j] + 1;

        if (with_gquad)
          f5[j] = MIN2(f5[j], ggg[ij]);

        if (c[ij] != INF) {
          if (evaluate(1, j, 1, j, VRNA_DECOMP_EXT_STEM, &hc_dat_local)) {
            en = c[ij];

            if (scs) {
              for (s = 0; s < n_seq; s++) {
                tt  = get_pair_type_md(S[s][i], S[s][j], md);
                mm3 = (a2s[s][j] < a2s[s][S[0][0]]) ? S3[s][j] : -1;      /* why S[0][0] ??? */
                en  += E_ExtLoop(tt, -1, mm3, P);

                if (scs[s] && scs[s]->f)
                  en += scs[s]->f(1, j, 1, j, VRNA_DECOMP_EXT_STEM, scs[s]->data);
              }
            } else {
              for (s = 0; s < n_seq; s++) {
                tt  = get_pair_type_md(S[s][i], S[s][j], md);
                mm3 = (a2s[s][j] < a2s[s][S[0][0]]) ? S3[s][j] : -1;      /* why S[0][0] ??? */
                en  += E_ExtLoop(tt, -1, mm3, P);
              }
            }

            f5[j] = MIN2(f5[j], en);
          }
        }
      }

      f5[length] = INF;
      if (f5[length - 1] != INF) {
        if (evaluate(1, length, 1, length - 1, VRNA_DECOMP_EXT_EXT, &hc_dat_local)) {
          f5[length] = f5[length - 1];
          if (scs) {
            for (s = 0; s < n_seq; s++) {
              if (scs[s]) {
                if (scs[s]->energy_up)
                  f5[length] += scs[s]->energy_up[a2s[s][length]][1];

                if (scs[s]->f) {
                  f5[length] += scs[s]->f(1,
                                          length,
                                          1,
                                          length - 1,
                                          VRNA_DECOMP_EXT_EXT,
                                          scs[s]->data);
                }
              }
            }
          }
        }
      }

      if (scs) {
        for (i = length - turn - 1; i > 1; i--) {
          if (f5[i - 1] != INF) {
            ij = indx[length] + i;

            if (with_gquad)
              f5[length] = MIN2(f5[length], f5[i - 1] + ggg[ij]);

            if (c[ij] != INF) {
              if (evaluate(1, length, i - 1, i, VRNA_DECOMP_EXT_EXT_STEM, &hc_dat_local)) {
                en = f5[i - 1] +
                     c[ij];

                for (s = 0; s < n_seq; s++) {
                  tt  = get_pair_type_md(S[s][i], S[s][j], md);
                  mm5 = (a2s[s][i] > 1) ? S5[s][i] : -1;
                  en  += E_ExtLoop(tt, mm5, -1, P);

                  if (scs[s] && scs[s]->f)
                    en += scs[s]->f(1, length, i - 1, i, VRNA_DECOMP_EXT_EXT_STEM, scs[s]->data);
                }

                f5[length] = MIN2(f5[length], en);
              }
            }
          }
        }
      } else {
        for (i = length - turn - 1; i > 1; i--) {
          if (f5[i - 1] != INF) {
            ij = indx[length] + i;

            if (with_gquad)
              f5[length] = MIN2(f5[length], f5[i - 1] + ggg[ij]);

            if (c[ij] != INF) {
              if (evaluate(1, length, i - 1, i, VRNA_DECOMP_EXT_EXT_STEM, &hc_dat_local)) {
                en = f5[i - 1] +
                     c[ij];

                for (s = 0; s < n_seq; s++) {
                  tt  = get_pair_type_md(S[s][i], S[s][j], md);
                  mm5 = (a2s[s][i] > 1) ? S5[s][i] : -1;
                  en  += E_ExtLoop(tt, mm5, -1, P);
                }

                f5[length] = MIN2(f5[length], en);
              }
            }
          }
        }
      }

      ij = indx[length] + 1;

      if (with_gquad)
        f5[length] = MIN2(f5[length], ggg[ij]);

      if (c[ij] != INF) {
        if (evaluate(1, length, 1, length, VRNA_DECOMP_EXT_STEM, &hc_dat_local)) {
          en = c[ij];

          if (scs) {
            for (s = 0; s < n_seq; s++) {
              tt  = get_pair_type_md(S[s][i], S[s][j], md);
              en  += E_ExtLoop(tt, -1, -1, P);

              if (scs[s] && scs[s]->f)
                en += scs[s]->f(1, length, 1, length, VRNA_DECOMP_EXT_STEM, scs[s]->data);
            }
          } else {
            for (s = 0; s < n_seq; s++) {
              tt  = get_pair_type_md(S[s][i], S[s][j], md);
              en  += E_ExtLoop(tt, -1, -1, P);
            }
          }

          f5[length] = MIN2(f5[length], en);
        }
      }

      break;
  }

  return f5[length];
}


PRIVATE int
E_ext_loop_3(vrna_fold_compound_t *fc,
             int                  i)
{
  char                      **ptype;
  short                     *S1;
  int                       e, dangle_model, *f3, j, turn, length, maxdist, with_gquad, **ggg,
                            energy, type, **c;
  vrna_param_t              *P;
  vrna_md_t                 *md;
  vrna_hc_t                 *hc;
  vrna_sc_t                 *sc;
  vrna_callback_hc_evaluate *evaluate;
  struct default_data       hc_dat_local;

  e = INF;

  length        = fc->length;
  maxdist       = fc->window_size;
  S1            = fc->sequence_encoding;
  ptype         = fc->ptype_local;
  P             = fc->params;
  md            = &(P->model_details);
  hc            = fc->hc;
  sc            = fc->sc;
  f3            = fc->matrices->f3_local;
  c             = fc->matrices->c_local;
  ggg           = fc->matrices->ggg_local;
  dangle_model  = md->dangles;
  turn          = md->min_loop_size;
  with_gquad    = md->gquad;

  hc_dat_local.mx_window  = hc->matrix_local;
  hc_dat_local.hc_up      = hc->up_ext;
  hc_dat_local.sn         = fc->strand_number;

  if (fc->hc->f) {
    evaluate            = &hc_default_user_window;
    hc_dat_local.hc_f   = fc->hc->f;
    hc_dat_local.hc_dat = fc->hc->data;
  } else {
    evaluate = &hc_default_window;
  }

  /* first case: i stays unpaired */
  if (evaluate(i, length, i + 1, length, VRNA_DECOMP_EXT_EXT, &hc_dat_local)) {
    e = f3[i + 1];
    if (sc) {
      if (sc->energy_up)
        e += sc->energy_up[i][1];

      if (sc->f)
        e += sc->f(i, length, i + 1, length, VRNA_DECOMP_EXT_EXT, sc->data);
    }
  }

  /* next all cases where i is paired */
  switch (dangle_model) {
    /* dont use dangling end and mismatch contributions at all */
    case 0:
      for (j = i + turn + 1; j < length && j <= i + maxdist; j++) {
        if ((with_gquad) && (f3[j + 1] != INF) && (ggg[i][j - i] != INF)) {
          energy = f3[j + 1] +
                   ggg[i][j - i];

          e = MIN2(e, energy);
        }

        if (evaluate(i, length, j, j + 1, VRNA_DECOMP_EXT_STEM_EXT, &hc_dat_local)) {
          if ((f3[j + 1] != INF) && (c[i][j - i] != INF)) {
            type = get_pair_type_window(i, j, ptype);

            energy = f3[j + 1] +
                     c[i][j - i] +
                     E_ExtLoop(type, -1, -1, P);

            if ((sc) && (sc->f))
              energy += sc->f(i, length, j, j + 1, VRNA_DECOMP_EXT_STEM_EXT, sc->data);

            e = MIN2(e, energy);
          }
        }
      }
      if (length <= i + maxdist) {
        j = length;

        if (with_gquad && (ggg[i][j - i] != INF))
          e = MIN2(e, ggg[i][j - i]);

        if (evaluate(i, length, i, j, VRNA_DECOMP_EXT_STEM, &hc_dat_local)) {
          if (c[i][j - i] != INF) {
            type = get_pair_type_window(i, j, ptype);

            energy = c[i][j - i] +
                     E_ExtLoop(type, -1, -1, P);

            if ((sc) && (sc->f))
              energy += sc->f(i, length, i, j, VRNA_DECOMP_EXT_STEM, sc->data);

            e = MIN2(e, energy);
          }
        }
      }

      break;
    /* always use dangle_model on both sides */
    case 2:
      for (j = i + turn + 1; j < length && j <= i + maxdist; j++) {
        if ((with_gquad) && (ggg[i][j - i] != INF) && (f3[j + 1] != INF)) {
          energy = f3[j + 1] +
                   ggg[i][j - i];

          e = MIN2(e, energy);
        }

        if (evaluate(i, length, j, j + 1, VRNA_DECOMP_EXT_STEM_EXT, &hc_dat_local)) {
          if ((f3[j + 1] != INF) && (c[i][j - i] != INF)) {
            type = get_pair_type_window(i, j, ptype);

            energy = f3[j + 1] +
                     c[i][j - i] +
                     E_ExtLoop(type,
                               (i > 1) ? S1[i - 1] : -1,
                               S1[j + 1],
                               P);

            if ((sc) && (sc->f))
              energy += sc->f(i, length, j, j + 1, VRNA_DECOMP_EXT_STEM_EXT, sc->data);

            e = MIN2(e, energy);
          }
        }
      }
      if (length <= i + maxdist) {
        j = length;

        if (with_gquad && (ggg[i][j - i] != INF))
          e = MIN2(e, ggg[i][j - i]);

        if (evaluate(i, j, i, j, VRNA_DECOMP_EXT_STEM, &hc_dat_local)) {
          if (c[i][j - i] != INF) {
            type = get_pair_type_window(i, j, ptype);

            energy = c[i][j - i] +
                     E_ExtLoop(type, (i > 1) ? S1[i - 1] : -1, -1, P);

            if ((sc) && (sc->f))
              energy += sc->f(i, length, i, j, VRNA_DECOMP_EXT_STEM, sc->data);

            e = MIN2(e, energy);
          }
        }
      }

      break;
    /* normal dangle_model, aka dangle_model = 1 */
    default:
      for (j = i + turn + 1; j < length && j <= i + maxdist; j++) {
        if (with_gquad && (f3[j + 1] != INF) && (ggg[i][j - i] != INF)) {
          energy = f3[j + 1] +
                   ggg[i][j - i];

          e = MIN2(e, energy);
        }

        type = get_pair_type_window(i, j, ptype);

        if (evaluate(i, length, j, j + 1, VRNA_DECOMP_EXT_STEM_EXT, &hc_dat_local)) {
          if ((f3[j + 1] != INF) && (c[i][j - i] != INF)) {
            energy = f3[j + 1] +
                     c[i][j - i] +
                     E_ExtLoop(type, -1, -1, P);

            if ((sc) && (sc->f))
              energy += sc->f(i, length, j, j + 1, VRNA_DECOMP_EXT_STEM_EXT, sc->data);

            e = MIN2(e, energy);
          }
        }

        if (j + 2 <= length) {
          if (evaluate(i, length, j, j + 2, VRNA_DECOMP_EXT_STEM_EXT, &hc_dat_local)) {
            if ((c[i][j - i] != INF) && (f3[j + 2] != INF)) {
              energy = c[i][j - i] +
                       f3[j + 2] +
                       E_ExtLoop(type, -1, S1[j + 1], P);

              if (sc) {
                if (sc->energy_up)
                  energy += sc->energy_up[j + 1][1];

                if (sc->f)
                  energy += sc->f(i, length, j, j + 2, VRNA_DECOMP_EXT_STEM_EXT, sc->data);
              }

              e = MIN2(e, energy);
            }
          }
        } else {
          if (evaluate(i, length, i, j, VRNA_DECOMP_EXT_STEM, &hc_dat_local)) {
            if (c[i][j - i] != INF) {
              energy = c[i][j - i] +
                       E_ExtLoop(type, -1, S1[j + 1], P);

              if ((sc) && (sc->f))
                energy += sc->f(i, length, i, j, VRNA_DECOMP_EXT_STEM, sc->data);

              e = MIN2(e, energy);
            }
          }
        }

        if (evaluate(i, length, j, j + 1, VRNA_DECOMP_EXT_STEM_EXT1, &hc_dat_local)) {
          type = get_pair_type_window(i + 1, j, ptype);

          if ((c[i + 1][j - i - 1] != INF) && (f3[j + 1] != INF)) {
            energy = f3[j + 1] +
                     c[i + 1][j - i - 1] +
                     E_ExtLoop(type, S1[i], -1, P);

            if (sc) {
              if (sc->energy_up)
                energy += sc->energy_up[i][1];

              if (sc->f)
                energy += sc->f(i, length, j, j + 1, VRNA_DECOMP_EXT_STEM_EXT1, sc->data);
            }

            e = MIN2(e, energy);
          }
        }

        if (j + 2 <= length) {
          if (evaluate(i, length, j, j + 2, VRNA_DECOMP_EXT_STEM_EXT1, &hc_dat_local)) {
            if ((c[i + 1][j - i - 1] != INF) && (f3[j + 2] != INF)) {
              energy = c[i + 1][j - i - 1] +
                       f3[j + 2] +
                       E_ExtLoop(type, S1[i], S1[j + 1], P);

              if (sc) {
                if (sc->energy_up)
                  energy += sc->energy_up[i][1] +
                            sc->energy_up[j + 1][1];

                if (sc->f)
                  energy += sc->f(i, length, j, j + 2, VRNA_DECOMP_EXT_STEM_EXT1, sc->data);
              }

              e = MIN2(e, energy);
            }
          }
        } else {
          if (evaluate(i, length, i + 1, j, VRNA_DECOMP_EXT_STEM, &hc_dat_local)) {
            if (c[i + 1][j - i - 1] != INF) {
              energy = c[i + 1][j - i - 1] +
                       E_ExtLoop(type, S1[i], S1[j + 1], P);

              if (sc) {
                if (sc->energy_up)
                  energy += sc->energy_up[i][1];

                if (sc->f)
                  energy += sc->f(i, length, i + 1, j, VRNA_DECOMP_EXT_STEM, sc->data);
              }

              e = MIN2(e, energy);
            }
          }
        }
      }

      if (length <= i + maxdist) {
        j = length;

        if (with_gquad && (ggg[i][j - i] != INF))
          e = MIN2(e, ggg[i][j - i]);

        if (evaluate(i, j, i, j, VRNA_DECOMP_EXT_STEM, &hc_dat_local)) {
          if (c[i][j - i] != INF) {
            type = get_pair_type_window(i, j, ptype);

            energy = c[i][j - i] +
                     E_ExtLoop(type, -1, -1, P);

            if ((sc) && (sc->f))
              energy += sc->f(i, j, i, j, VRNA_DECOMP_EXT_STEM, sc->data);

            e = MIN2(e, energy);
          }
        }

        if (evaluate(i, j, i + 1, j, VRNA_DECOMP_EXT_STEM, &hc_dat_local)) {
          if (c[i + 1][j - i - 1] != INF) {
            type = get_pair_type_window(i + 1, j, ptype);

            energy = c[i + 1][j - i - 1] +
                     E_ExtLoop(type, S1[i], -1, P);

            if (sc) {
              if (sc->energy_up)
                energy += sc->energy_up[i][1];

              if (sc->f)
                energy += sc->f(i, j, i + 1, j, VRNA_DECOMP_EXT_STEM, sc->data);
            }

            e = MIN2(e, energy);
          }
        }
      }

      break;
  } /* switch(dangle_model)... */

  return e;
}


PRIVATE int
E_ext_loop_3_comparative(vrna_fold_compound_t *fc,
                         int                  i)
{
  short                     **S, **S5, **S3;
  int                       e, dangle_model, *f3, j, turn, length, maxdist, with_gquad, **ggg,
                            energy, **c, n_seq, s, tt;
  vrna_param_t              *P;
  vrna_md_t                 *md;
  vrna_hc_t                 *hc;
  vrna_callback_hc_evaluate *evaluate;
  struct default_data       hc_dat_local;

  e = INF;

  length        = fc->length;
  n_seq         = fc->n_seq;
  S             = fc->S;
  S5            = fc->S5;     /* S5[s][i] holds next base 5' of i in sequence s */
  S3            = fc->S3;     /* Sl[s][i] holds next base 3' of i in sequence s */
  maxdist       = fc->window_size;
  P             = fc->params;
  md            = &(P->model_details);
  hc            = fc->hc;
  f3            = fc->matrices->f3_local;
  c             = fc->matrices->c_local;
  ggg           = fc->matrices->ggg_local;
  dangle_model  = md->dangles;
  turn          = md->min_loop_size;
  with_gquad    = md->gquad;

  hc_dat_local.mx_window  = hc->matrix_local;
  hc_dat_local.hc_up      = hc->up_ext;
  hc_dat_local.sn         = fc->strand_number;

  if (fc->hc->f) {
    evaluate            = &hc_default_user_window;
    hc_dat_local.hc_f   = fc->hc->f;
    hc_dat_local.hc_dat = fc->hc->data;
  } else {
    evaluate = &hc_default_window;
  }

  /* first case: i stays unpaired */
  if (evaluate(i, length, i + 1, length, VRNA_DECOMP_EXT_EXT, &hc_dat_local))
    e = f3[i + 1];

  /* next all cases where i is paired */
  switch (dangle_model) {
    /* dont use dangling end and mismatch contributions at all */
    case 0:
      for (j = i + turn + 1; j < length && j <= i + maxdist; j++) {
        if ((with_gquad) && (f3[j + 1] != INF) && (ggg[i][j - i] != INF)) {
          energy = f3[j + 1] +
                   ggg[i][j - i];

          e = MIN2(e, energy);
        }

        if (evaluate(i, length, j, j + 1, VRNA_DECOMP_EXT_STEM_EXT, &hc_dat_local)) {
          if ((f3[j + 1] != INF) && (c[i][j - i] != INF)) {
            energy = f3[j + 1] +
                     c[i][j - i];

            for (s = 0; s < n_seq; s++) {
              tt      = get_pair_type_md(S[s][i], S[s][j], md);
              energy  += E_ExtLoop(tt, -1, -1, P);
            }

            e = MIN2(e, energy);
          }
        }
      }
      if (length <= i + maxdist) {
        j = length;

        if (with_gquad && (ggg[i][j - i] != INF))
          e = MIN2(e, ggg[i][j - i]);

        if (evaluate(i, length, i, j, VRNA_DECOMP_EXT_STEM, &hc_dat_local)) {
          if (c[i][j - i] != INF) {
            energy = c[i][j - i];
            for (s = 0; s < n_seq; s++) {
              tt      = get_pair_type_md(S[s][i], S[s][j], md);
              energy  += E_ExtLoop(tt, -1, -1, P);
            }
            e = MIN2(e, energy);
          }
        }
      }

      break;
    /* always use dangle_model on both sides */
    case 2:
      for (j = i + turn + 1; j < length && j <= i + maxdist; j++) {
        if ((with_gquad) && (ggg[i][j - i] != INF) && (f3[j + 1] != INF)) {
          energy = f3[j + 1] +
                   ggg[i][j - i];

          e = MIN2(e, energy);
        }

        if (evaluate(i, length, j, j + 1, VRNA_DECOMP_EXT_STEM_EXT, &hc_dat_local)) {
          if ((f3[j + 1] != INF) && (c[i][j - i] != INF)) {
            energy = f3[j + 1] +
                     c[i][j - i];

            for (s = 0; s < n_seq; s++) {
              tt      = get_pair_type_md(S[s][i], S[s][j], md);
              energy  += E_ExtLoop(tt, (i > 1) ? S5[s][i] : -1, S3[s][j], P);
            }
            e = MIN2(e, energy);
          }
        }
      }
      if (length <= i + maxdist) {
        j = length;

        if (with_gquad && (ggg[i][j - i] != INF))
          e = MIN2(e, ggg[i][j - i]);

        if (evaluate(i, j, i, j, VRNA_DECOMP_EXT_STEM, &hc_dat_local)) {
          if (c[i][j - i] != INF) {
            energy = c[i][j - i];
            for (s = 0; s < n_seq; s++) {
              tt      = get_pair_type_md(S[s][i], S[s][j], md);
              energy  += E_ExtLoop(tt, (i > 1) ? S5[s][i] : -1, -1, P);
            }
            e = MIN2(e, energy);
          }
        }
      }

      break;
  } /* switch(dangle_model)... */

  return e;
}
