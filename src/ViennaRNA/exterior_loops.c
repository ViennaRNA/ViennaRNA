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

struct default_data {
  int                       *idx;
  char                      *mx;
  int                       cp;
  int                       *hc_up;
  void                      *hc_dat;
  vrna_callback_hc_evaluate *hc_f;
};


/*
 #################################
 # PRIVATE FUNCTION DECLARATIONS #
 #################################
 */

PRIVATE FLT_OR_DBL
exp_E_ext_fast(vrna_fold_compound_t *vc,
               int                  i,
               int                  j,
               vrna_mx_pf_aux_el_t  *aux_mx);


PRIVATE FLT_OR_DBL
exp_E_ext_fast_comparative(vrna_fold_compound_t *vc,
                           int                  i,
                           int                  j,
                           vrna_mx_pf_aux_el_t  *aux_mx);


PRIVATE char
hc_default(int  i,
           int  j,
           int  k,
           int  l,
           char d,
           void *data);


PRIVATE char
hc_default_user(int   i,
                int   j,
                int   k,
                int   l,
                char  d,
                void  *data);


/*
 #################################
 # BEGIN OF FUNCTION DEFINITIONS #
 #################################
 */
PUBLIC int
E_ext_loop(int                  i,
           int                  j,
           vrna_fold_compound_t *vc)
{
  char                      *ptype, *hard_constraints;
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
  hc_dat_local.cp   = cp;

  if (vc->hc->f) {
    evaluate            = &hc_default_user;
    hc_dat_local.hc_f   = vc->hc->f;
    hc_dat_local.hc_dat = vc->hc->data;
  } else {
    evaluate = &hc_default;
  }


  e     = INF;
  ij    = idx[j] + i;
  type  = ptype[ij];

  if ((cp < 0) || (((i) >= cp) || ((j) < cp))) {
    /* regular exterior loop */
    if (evaluate(i, j, i, j, VRNA_DECOMP_EXT_STEM, &hc_dat_local)) {
      if (type == 0)
        type = 7;

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
        type = vc->ptype[ij];

        if (type == 0)
          type = 7;

        en = E_ExtLoop(type, -1, S[j], P);
        if (sc)
          if (sc->f)
            en += sc->f(i, j, i, j - 1, VRNA_DECOMP_EXT_STEM, sc->data);
        e = MIN2(e, en);
      }

      ij = idx[j] + i + 1;
      if (evaluate(i, j, i + 1, j, VRNA_DECOMP_EXT_STEM, &hc_dat_local)) {
        type = vc->ptype[ij];

        if (type == 0)
          type = 7;

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


PUBLIC void
E_ext_loop_5(vrna_fold_compound_t *vc)
{
  char                      *ptype, *hc;
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
  hc_dat_local.cp     = vc->cutpoint;

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
              en = f5[i - 1]
                   + domains_up->energy_cb(vc,
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
                en = f5[j - u]
                     + domains_up->energy_cb(vc,
                                             j - u + 1, j,
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
                  type = ptype[ij];

                  if (type == 0)
                    type = 7;

                  en    = f5[i - 1] + c[ij] + E_ExtLoop(type, -1, -1, P);
                  en    += sc->f(1, j, i - 1, i, VRNA_DECOMP_EXT_EXT_STEM, sc->data);
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
                  type = ptype[ij];

                  if (type == 0)
                    type = 7;

                  en    = f5[i - 1] + c[ij] + E_ExtLoop(type, -1, -1, P);
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
            type = ptype[ij];

            if (type == 0)
              type = 7;

            en = c[ij] + E_ExtLoop(type, -1, -1, P);
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
                en = f5[j - u]
                     + domains_up->energy_cb(vc,
                                             j - u + 1, j,
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
                  type = ptype[ij];

                  if (type == 0)
                    type = 7;

                  en    = f5[i - 1] + c[ij] + E_ExtLoop(type, S[i - 1], S[j + 1], P);
                  en    += sc->f(1, j, i - 1, i, VRNA_DECOMP_EXT_EXT_STEM, sc->data);
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
                  type = ptype[ij];

                  if (type == 0)
                    type = 7;

                  en    = f5[i - 1] + c[ij] + E_ExtLoop(type, S[i - 1], S[j + 1], P);
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
            type = ptype[ij];

            if (type == 0)
              type = 7;

            en = c[ij] + E_ExtLoop(type, -1, S[j + 1], P);
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
              en = f5[length - u]
                   + domains_up->energy_cb(vc,
                                           length - u + 1, length,
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
                type = ptype[ij];

                if (type == 0)
                  type = 7;

                en          = f5[i - 1] + c[ij] + E_ExtLoop(type, S[i - 1], -1, P);
                en          += sc->f(1, length, i - 1, i, VRNA_DECOMP_EXT_EXT_STEM, sc->data);
                f5[length]  = MIN2(f5[length], en);
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
                type = ptype[ij];

                if (type == 0)
                  type = 7;

                en          = f5[i - 1] + c[ij] + E_ExtLoop(type, S[i - 1], -1, P);
                f5[length]  = MIN2(f5[length], en);
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
          type = ptype[ij];

          if (type == 0)
            type = 7;

          en = c[ij] + E_ExtLoop(type, -1, -1, P);
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
                en = f5[j - u]
                     + domains_up->energy_cb(vc,
                                             j - u + 1, j,
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
                type = ptype[ij];

                if (type == 0)
                  type = 7;

                en = f5[i - 1] + c[ij] + E_ExtLoop(type, -1, -1, P);
                if (sc)
                  if (sc->f)
                    en += sc->f(1, j, i - 1, i, VRNA_DECOMP_EXT_EXT_STEM, sc->data);
                f5[j] = MIN2(f5[j], en);
              }
            }
          }

          if ((f5[i - 2] != INF) && c[ij] != INF) {
            if (evaluate(1, j, i - 2, i, VRNA_DECOMP_EXT_EXT_STEM, &hc_dat_local)) {
              type = ptype[ij];

              if (type == 0)
                type = 7;

              en = f5[i - 2] + c[ij] + E_ExtLoop(type, S[i - 1], -1, P);

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
                type = ptype[ij];

                if (type == 0)
                  type = 7;

                en = f5[i - 1] + c[ij] + E_ExtLoop(type, -1, S[j], P);

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
                type = ptype[ij];

                if (type == 0)
                  type = 7;

                en = f5[i - 2] + c[ij] + E_ExtLoop(type, S[i - 1], S[j], P);

                if (sc) {
                  if (sc->energy_up)
                    en += sc->energy_up[i - 1][1] + sc->energy_up[j][1];
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
            type = ptype[ij];

            if (type == 0)
              type = 7;

            en = c[ij] + E_ExtLoop(type, -1, -1, P);
            if (sc)
              if (sc->f)
                en += sc->f(1, j, 1, j, VRNA_DECOMP_EXT_STEM, sc->data);
            f5[j] = MIN2(f5[j], en);
          }
        }
        ij = indx[j - 1] + 1;
        if (c[ij] != INF) {
          if (evaluate(1, j, 1, j - 1, VRNA_DECOMP_EXT_STEM, &hc_dat_local)) {
            type = ptype[ij];

            if (type == 0)
              type = 7;

            en = c[ij] + E_ExtLoop(type, -1, S[j], P);

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


PUBLIC FLT_OR_DBL
exp_E_Stem(int              type,
           int              si1,
           int              sj1,
           int              extLoop,
           vrna_exp_param_t *P)
{
  double  energy  = 1.0;
  double  d5      = (si1 >= 0) ? P->expdangle5[type][si1] : 1.;
  double  d3      = (sj1 >= 0) ? P->expdangle3[type][sj1] : 1.;

  if (si1 >= 0 && sj1 >= 0)
    energy = (extLoop) ? P->expmismatchExt[type][si1][sj1] : P->expmismatchM[type][si1][sj1];
  else
    energy = d5 * d3;

  if (type > 2)
    energy *= P->expTermAU;

  if (!extLoop)
    energy *= P->expMLintern[type];

  return (FLT_OR_DBL)energy;
}


PUBLIC FLT_OR_DBL
exp_E_ExtLoop(int               type,
              int               si1,
              int               sj1,
              vrna_exp_param_t  *P)
{
  double energy = 1.0;

  if (si1 >= 0 && sj1 >= 0)
    energy = P->expmismatchExt[type][si1][sj1];
  else if (si1 >= 0)
    energy = P->expdangle5[type][si1];
  else if (sj1 >= 0)
    energy = P->expdangle3[type][sj1];

  if (type > 2)
    energy *= P->expTermAU;

  return (FLT_OR_DBL)energy;
}


PUBLIC int
vrna_BT_ext_loop_f5(vrna_fold_compound_t  *vc,
                    int                   *k,
                    int                   *i,
                    int                   *j,
                    vrna_bp_stack_t       *bp_stack,
                    int                   *stack_count)
{
  unsigned char             type;
  char                      *ptype;
  short                     mm5, mm3, *S1;
  unsigned int              *sn;
  int                       length, fij, fi, jj, u, en, e, *my_f5, *my_c, *my_ggg, *idx,
                            dangle_model, turn, with_gquad, cnt, ii, with_ud;
  vrna_param_t              *P;
  vrna_md_t                 *md;
  vrna_hc_t                 *hc;
  vrna_sc_t                 *sc;
  vrna_ud_t                 *domains_up;
  vrna_callback_hc_evaluate *evaluate;
  struct default_data       hc_dat_local;

  length        = vc->length;
  P             = vc->params;
  md            = &(P->model_details);
  sn            = vc->strand_number;
  hc            = vc->hc;
  sc            = vc->sc;
  my_f5         = vc->matrices->f5;
  my_c          = vc->matrices->c;
  my_ggg        = vc->matrices->ggg;
  domains_up    = vc->domains_up;
  idx           = vc->jindx;
  ptype         = vc->ptype;
  S1            = vc->sequence_encoding;
  dangle_model  = md->dangles;
  turn          = md->min_loop_size;
  with_gquad    = md->gquad;
  with_ud       = (domains_up && domains_up->energy_cb) ? 1 : 0;

  hc_dat_local.idx    = idx;
  hc_dat_local.mx     = hc->matrix;
  hc_dat_local.hc_up  = hc->up_ext;
  hc_dat_local.cp     = vc->cutpoint;

  if (vc->hc->f) {
    evaluate            = &hc_default_user;
    hc_dat_local.hc_f   = hc->f;
    hc_dat_local.hc_dat = hc->data;
  } else {
    evaluate = &hc_default;
  }

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
          en = domains_up->energy_cb(vc,
                                     ii, jj,
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
            vrna_BT_gquad_mfe(vc, u, jj, bp_stack, stack_count);
            return 1;
          }
        }

        if (evaluate(1, jj, u - 1, u, VRNA_DECOMP_EXT_EXT_STEM, &hc_dat_local)) {
          type = (unsigned char)ptype[idx[jj] + u];

          if (type == 0)
            type = 7;

          en = my_c[idx[jj] + u];
          if (sc)
            if (sc->f)
              en += sc->f(1, jj, u - 1, u, VRNA_DECOMP_EXT_EXT_STEM, sc->data);
          if (sn[jj] != sn[u])
            en += P->DuplexInit;
          if (fij == E_ExtLoop(type, -1, -1, P) + en + my_f5[u - 1]) {
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
            vrna_BT_gquad_mfe(vc, u, jj, bp_stack, stack_count);
            return 1;
          }
        }

        if (evaluate(1, jj, u - 1, u, VRNA_DECOMP_EXT_EXT_STEM, &hc_dat_local)) {
          mm5   = ((u > 1) && (sn[u] == sn[u - 1])) ? S1[u - 1] : -1;
          type  = (unsigned char)ptype[idx[jj] + u];

          if (type == 0)
            type = 7;

          en = my_c[idx[jj] + u];
          if (sc)
            if (sc->f)
              en += sc->f(1, jj, u - 1, u, VRNA_DECOMP_EXT_EXT_STEM, sc->data);
          if (sn[jj] != sn[u])
            en += P->DuplexInit;
          if (fij == E_ExtLoop(type, mm5, mm3, P) + en + my_f5[u - 1]) {
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
          vrna_BT_gquad_mfe(vc, 1, jj, bp_stack, stack_count);
          return 1;
        }
      }

      if (evaluate(1, jj, 1, jj, VRNA_DECOMP_EXT_STEM, &hc_dat_local)) {
        type = (unsigned char)ptype[idx[jj] + 1];

        if (type == 0)
          type = 7;

        en = my_c[idx[jj] + 1];
        if (sc)
          if (sc->f)
            en += sc->f(1, jj, 1, jj, VRNA_DECOMP_EXT_STEM, sc->data);
        if (sn[jj] != sn[1])
          en += P->DuplexInit;
        if (fij == en + E_ExtLoop(type, -1, -1, P)) {
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
          type  = (unsigned char)ptype[idx[jj - 1] + 1];

          if (type == 0)
            type = 7;

          en = my_c[idx[jj - 1] + 1];
          if (sc) {
            if (sc->energy_up)
              en += sc->energy_up[jj][1];
            if (sc->f)
              en += sc->f(1, jj, 1, jj - 1, VRNA_DECOMP_EXT_STEM, sc->data);
          }
          if (sn[jj - 1] != sn[1])
            en += P->DuplexInit;

          if (fij == en + E_ExtLoop(type, -1, mm3, P)) {
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
            vrna_BT_gquad_mfe(vc, u, jj, bp_stack, stack_count);
            return 1;
          }
        }

        type = (unsigned char)ptype[idx[jj] + u];
        if (type == 0)
          type = 7;

        en = my_c[idx[jj] + u];
        if (sn[jj] != sn[u])
          en += P->DuplexInit;

        if (evaluate(1, jj, u - 1, u, VRNA_DECOMP_EXT_EXT_STEM, &hc_dat_local)) {
          e = my_f5[u - 1] + en + E_ExtLoop(type, -1, -1, P);
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
            e   = my_f5[u - 2] + en + E_ExtLoop(type, mm5, -1, P);
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

        type = (unsigned char)ptype[idx[jj - 1] + u];
        if (type == 0)
          type = 7;

        en = my_c[idx[jj - 1] + u];
        if (sn[jj - 1] != sn[u])
          en += P->DuplexInit;

        mm5 = (sn[u] == sn[u - 1]) ? S1[u - 1] : -1;
        mm3 = (sn[jj] == sn[jj - 1]) ? S1[jj] : -1;

        if (evaluate(1, jj, u - 1, u, VRNA_DECOMP_EXT_EXT_STEM1, &hc_dat_local)) {
          e = my_f5[u - 1] + en + E_ExtLoop(type, -1, mm3, P);

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
          e = my_f5[u - 2] + en + E_ExtLoop(type, mm5, mm3, P);
          if (sc) {
            if (sc->energy_up)
              e += sc->energy_up[jj][1]
                   + sc->energy_up[u - 1][1];
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


PUBLIC vrna_mx_pf_aux_el_t *
vrna_exp_E_ext_fast_init(vrna_fold_compound_t *vc)
{
  vrna_mx_pf_aux_el_t *aux_mx = NULL;

  if (vc) {
    char                      *hc;
    unsigned int              u, s;
    int                       i, j, d, n, turn, ij, *idx, *iidx, *hc_up;
    FLT_OR_DBL                *q, *scale;
    vrna_callback_hc_evaluate *evaluate;
    struct default_data       hc_dat_local;

    n     = (int)vc->length;
    idx   = vc->jindx;
    iidx  = vc->iindx;
    turn  = vc->exp_params->model_details.min_loop_size;
    q     = vc->exp_matrices->q;
    scale = vc->exp_matrices->scale;
    hc    = vc->hc->matrix;
    hc_up = vc->hc->up_ext;

    hc_dat_local.idx    = idx;
    hc_dat_local.mx     = hc;
    hc_dat_local.hc_up  = hc_up;
    hc_dat_local.cp     = vc->cutpoint;

    if (vc->hc->f) {
      evaluate            = &hc_default_user;
      hc_dat_local.hc_f   = vc->hc->f;
      hc_dat_local.hc_dat = vc->hc->data;
    } else {
      evaluate = &hc_default;
    }


    /* allocate memory for helper arrays */
    aux_mx            = (vrna_mx_pf_aux_el_t *)vrna_alloc(sizeof(vrna_mx_pf_aux_el_t));
    aux_mx->qq        = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * (n + 2));
    aux_mx->qq1       = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * (n + 2));
    aux_mx->qqu_size  = 0;
    aux_mx->qqu       = NULL;

    if (vc->type == VRNA_FC_TYPE_SINGLE) {
      vrna_sc_t *sc         = vc->sc;
      vrna_ud_t *domains_up = vc->domains_up;
      int       with_ud     = (domains_up && domains_up->exp_energy_cb);

      /* pre-processing ligand binding production rule(s) and auxiliary memory */
      if (with_ud) {
        int ud_max_size = 0;
        for (u = 0; u < domains_up->uniq_motif_count; u++)
          if (ud_max_size < domains_up->uniq_motif_size[u])
            ud_max_size = domains_up->uniq_motif_size[u];

        aux_mx->qqu_size  = ud_max_size;
        aux_mx->qqu       = (FLT_OR_DBL **)vrna_alloc(sizeof(FLT_OR_DBL *) * (ud_max_size + 1));

        for (u = 0; u <= ud_max_size; u++)
          aux_mx->qqu[u] = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * (n + 2));
      }

      for (d = 0; d <= turn; d++)
        for (i = 1; i <= n - d; i++) {
          j   = i + d;
          ij  = iidx[i] - j;

          if (j > n)
            continue;

          if (evaluate(i, j, i, j, VRNA_DECOMP_EXT_UP, &hc_dat_local)) {
            q[ij] = 1.0 * scale[d + 1];

            if (sc) {
              if (sc->exp_energy_up)
                q[ij] *= sc->exp_energy_up[i][d + 1];
              if (sc->exp_f)
                q[ij] *= sc->exp_f(i, j, i, j, VRNA_DECOMP_EXT_UP, sc->data);
            }

            if (with_ud) {
              q[ij] += q[ij] * domains_up->exp_energy_cb(vc,
                                                         i, j,
                                                         VRNA_UNSTRUCTURED_DOMAIN_EXT_LOOP,
                                                         domains_up->data);
            }
          } else {
            q[ij] = 0.;
          }
        }
    } else if (vc->type == VRNA_FC_TYPE_COMPARATIVE) {
      vrna_sc_t       **scs = vc->scs;
      unsigned short  **a2s = vc->a2s;
      for (d = 0; d <= turn; d++)
        for (i = 1; i <= n - d; i++) {
          j   = i + d;
          ij  = iidx[i] - j;
          if (evaluate(i, j, i, j, VRNA_DECOMP_EXT_UP, &hc_dat_local)) {
            q[ij] = 1.0 * scale[d + 1];

            if (scs) {
              for (s = 0; s < vc->n_seq; s++)
                if (scs[s]) {
                  u = d + 1 /* a2s[s][j] - a2s[s][i] + 1 */;
                  if (scs[s]->exp_energy_up)
                    q[ij] *= scs[s]->exp_energy_up[a2s[s][i]][u];
                }
            }
          } else {
            q[ij] = 0.;
          }
        }
    }
  }

  return aux_mx;
}


PUBLIC void
vrna_exp_E_ext_fast_rotate(vrna_fold_compound_t *vc,
                           vrna_mx_pf_aux_el_t  *aux_mx)
{
  if (vc && aux_mx) {
    int         u;
    FLT_OR_DBL  *tmp;

    tmp         = aux_mx->qq1;
    aux_mx->qq1 = aux_mx->qq;
    aux_mx->qq  = tmp;

    /* rotate auxiliary arrays for unstructured domains */
    if (aux_mx->qqu) {
      tmp = aux_mx->qqu[aux_mx->qqu_size];
      for (u = aux_mx->qqu_size; u > 0; u--)
        aux_mx->qqu[u] = aux_mx->qqu[u - 1];
      aux_mx->qqu[0] = tmp;
    }
  }
}


PUBLIC void
vrna_exp_E_ext_fast_free(vrna_fold_compound_t *vc,
                         vrna_mx_pf_aux_el_t  *aux_mx)
{
  if (vc && aux_mx) {
    int u;

    free(aux_mx->qq);
    free(aux_mx->qq1);

    if (aux_mx->qqu) {
      for (u = 0; u <= aux_mx->qqu_size; u++)
        free(aux_mx->qqu[u]);

      free(aux_mx->qqu);
    }
    free(aux_mx);
  }
}


PUBLIC FLT_OR_DBL
vrna_exp_E_ext_fast(vrna_fold_compound_t  *vc,
                    int                   i,
                    int                   j,
                    vrna_mx_pf_aux_el_t   *aux_mx)
{
  if (vc) {
    switch (vc->type) {
      case VRNA_FC_TYPE_SINGLE:
        return exp_E_ext_fast(vc, i, j, aux_mx);
        break;

      case VRNA_FC_TYPE_COMPARATIVE:
        return exp_E_ext_fast_comparative(vc, i, j, aux_mx);
        break;

      default:
        vrna_message_warning("vrna_exp_E_ext_fast@exterior_loops.c: Unknown fold_compound type");
        return 0.;
        break;
    }
  } else {
    return 0.;
  }
}


PRIVATE FLT_OR_DBL
exp_E_ext_fast(vrna_fold_compound_t *vc,
               int                  i,
               int                  j,
               vrna_mx_pf_aux_el_t  *aux_mx)
{
  short                     *S1;
  unsigned char             type;
  int                       n, *iidx, k, ij, kl, with_ud, u, circular, with_gquad;
  FLT_OR_DBL                qbt1, *q, *qb, *qq, *qq1, **qqu, q_temp, *scale, q_temp2, *G;
  vrna_md_t                 *md;
  vrna_exp_param_t          *pf_params;
  vrna_ud_t                 *domains_up;
  vrna_hc_t                 *hc;
  vrna_sc_t                 *sc;
  vrna_callback_hc_evaluate *evaluate;
  struct default_data       hc_dat_local;

  n                   = (int)vc->length;
  iidx                = vc->iindx;
  ij                  = iidx[i] - j;
  qq                  = aux_mx->qq;
  qq1                 = aux_mx->qq1;
  qqu                 = aux_mx->qqu;
  q                   = vc->exp_matrices->q;
  qb                  = vc->exp_matrices->qb;
  G                   = vc->exp_matrices->G;
  scale               = vc->exp_matrices->scale;
  pf_params           = vc->exp_params;
  md                  = &(pf_params->model_details);
  hc                  = vc->hc;
  sc                  = vc->sc;
  domains_up          = vc->domains_up;
  circular            = md->circ;
  with_gquad          = md->gquad;
  with_ud             = (domains_up && domains_up->exp_energy_cb);
  hc_dat_local.idx    = vc->jindx;
  hc_dat_local.mx     = hc->matrix;
  hc_dat_local.hc_up  = hc->up_ext;
  hc_dat_local.cp     = vc->cutpoint;

  if (hc->f) {
    evaluate            = &hc_default_user;
    hc_dat_local.hc_f   = hc->f;
    hc_dat_local.hc_dat = hc->data;
  } else {
    evaluate = &hc_default;
  }

  qbt1 = 0.;

  /* all exterior loop parts [i, j] with exactly one stem (i, u) i < u < j */
  if (evaluate(i, j, i, j - 1, VRNA_DECOMP_EXT_EXT, &hc_dat_local)) {
    q_temp = qq1[i] * scale[1];

    if (sc) {
      if (sc->exp_energy_up)
        q_temp *= sc->exp_energy_up[j][1];

      if (sc->exp_f)
        q_temp *= sc->exp_f(i, j, i, j - 1, VRNA_DECOMP_EXT_EXT, sc->data);
    }

    if (with_ud) {
      int cnt;
      for (cnt = 0; cnt < domains_up->uniq_motif_count; cnt++) {
        u = domains_up->uniq_motif_size[cnt];
        if (j - u >= i) {
          if (evaluate(i, j, i, j - u, VRNA_DECOMP_EXT_EXT, &hc_dat_local)) {
            q_temp2 = qqu[u][i]
                      * domains_up->exp_energy_cb(vc,
                                                  j - u + 1, j,
                                                  VRNA_UNSTRUCTURED_DOMAIN_EXT_LOOP | VRNA_UNSTRUCTURED_DOMAIN_MOTIF,
                                                  domains_up->data)
                      * scale[u];

            if (sc) {
              if (sc->exp_energy_up)
                q_temp2 *= sc->exp_energy_up[j - u + 1][u];
              if (sc->exp_f)
                q_temp2 *= sc->exp_f(i, j, i, j - u, VRNA_DECOMP_EXT_EXT, sc->data);
            }

            q_temp += q_temp2;
          }
        }
      }
    }

    qbt1 += q_temp;
  }

  /* exterior loop part with stem (i, j) */
  if (evaluate(i, j, i, j, VRNA_DECOMP_EXT_STEM, &hc_dat_local)) {
    S1    = vc->sequence_encoding;
    type  = md->pair[S1[i]][S1[j]];
    if (type == 0)
      type = 7;

    q_temp = qb[ij]
             * exp_E_ExtLoop(type, ((i > 1) || circular) ? S1[i - 1] : -1, ((j < n) || circular) ? S1[j + 1] : -1, pf_params);

    if (sc)
      if (sc->exp_f)
        q_temp *= sc->exp_f(i, j, i, j, VRNA_DECOMP_EXT_STEM, sc->data);
    qbt1 += q_temp;
  }

  if (with_gquad)
    qbt1 += G[ij];

  qq[i] = qbt1;

  if (with_ud)
    qqu[0][i] = qbt1;

  /* the entire stretch [i,j] is unpaired */
  if (evaluate(i, j, i, j, VRNA_DECOMP_EXT_UP, &hc_dat_local)) {
    u       = j - i + 1;
    q_temp  = 1.0 * scale[u];

    if (sc) {
      if (sc->exp_energy_up)
        q_temp *= sc->exp_energy_up[i][u];

      if (sc->exp_f)
        q_temp *= sc->exp_f(i, j, i, j, VRNA_DECOMP_EXT_UP, sc->data);
    }

    qbt1 += q_temp;

    if (with_ud) {
      qbt1 += q_temp * domains_up->exp_energy_cb(vc,
                                                 i, j,
                                                 VRNA_UNSTRUCTURED_DOMAIN_EXT_LOOP,
                                                 domains_up->data);
    }
  }

  kl = iidx[i] - j + 1;
  if (sc && sc->exp_f) {
    for (k = j; k > i; k--, kl++) {
      q_temp  = q[kl] * qq[k];
      q_temp  *= sc->exp_f(i, j, k - 1, k, VRNA_DECOMP_EXT_EXT_EXT, sc->data);
      qbt1    += q_temp;
    }
  } else {
    for (k = j; k > i; k--, kl++)
      qbt1 += q[kl] * qq[k];
  }

  return qbt1;
}


PRIVATE FLT_OR_DBL
exp_E_ext_fast_comparative(vrna_fold_compound_t *vc,
                           int                  i,
                           int                  j,
                           vrna_mx_pf_aux_el_t  *aux_mx)
{
  int                       n, s, n_seq, *iidx, k, ij, kl, u, circular, type;
  unsigned short            **a2s;
  short                     **S, **S5, **S3;
  FLT_OR_DBL                qbt1, *q, *qb, *qq, *qq1, q_temp, *scale;
  vrna_md_t                 *md;
  vrna_exp_param_t          *pf_params;
  vrna_hc_t                 *hc;
  vrna_sc_t                 **scs;
  vrna_callback_hc_evaluate *evaluate;
  struct default_data       hc_dat_local;

  n                   = (int)vc->length;
  n_seq               = vc->n_seq;
  iidx                = vc->iindx;
  ij                  = iidx[i] - j;
  S                   = vc->S;
  S5                  = vc->S5;     /* S5[s][i] holds next base 5' of i in sequence s */
  S3                  = vc->S3;     /* Sl[s][i] holds next base 3' of i in sequence s */
  a2s                 = vc->a2s;
  qq                  = aux_mx->qq;
  qq1                 = aux_mx->qq1;
  q                   = vc->exp_matrices->q;
  qb                  = vc->exp_matrices->qb;
  scale               = vc->exp_matrices->scale;
  pf_params           = vc->exp_params;
  md                  = &(pf_params->model_details);
  hc                  = vc->hc;
  scs                 = vc->scs;
  circular            = md->circ;
  hc_dat_local.idx    = vc->jindx;
  hc_dat_local.mx     = hc->matrix;
  hc_dat_local.hc_up  = hc->up_ext;
  hc_dat_local.cp     = vc->cutpoint;

  if (hc->f) {
    evaluate            = &hc_default_user;
    hc_dat_local.hc_f   = hc->f;
    hc_dat_local.hc_dat = hc->data;
  } else {
    evaluate = &hc_default;
  }

  qbt1 = 0.;

  /* all exterior loop parts [i, j] with exactly one stem (i, u) i < u < j */
  if (evaluate(i, j, i, j - 1, VRNA_DECOMP_EXT_EXT, &hc_dat_local)) {
    q_temp = qq1[i] * scale[1];

    if (scs) {
      for (s = 0; s < n_seq; s++) {
        if (scs[s])
          if (scs[s]->exp_energy_up)
            q_temp *= scs[s]->exp_energy_up[a2s[s][j]][1];
      }
    }

    qbt1 += q_temp;
  }

  /* exterior loop part with stem (i, j) */
  if (evaluate(i, j, i, j, VRNA_DECOMP_EXT_STEM, &hc_dat_local)) {
    q_temp = qb[ij];

    for (s = 0; s < n_seq; s++) {
      type = md->pair[S[s][i]][S[s][j]];
      if (type == 0)
        type = 7;

      q_temp *= exp_E_ExtLoop(type, ((i > 1) || circular) ? S5[s][i] : -1, ((j < n) || circular) ? S3[s][j] : -1, pf_params);
    }

    qbt1 += q_temp;
  }

  qq[i] = qbt1;

  /* the entire stretch [i,j] is unpaired */
  if (evaluate(i, j, i, j, VRNA_DECOMP_EXT_UP, &hc_dat_local)) {
    u       = j - i + 1;
    q_temp  = 1.0 * scale[u];

    if (scs) {
      for (s = 0; s < n_seq; s++) {
        if (scs[s])
          if (scs[s]->exp_energy_up)
            q_temp *= scs[s]->exp_energy_up[a2s[s][i]][a2s[s][j] - a2s[s][i] + 1];
      }
    }

    qbt1 += q_temp;
  }

  kl = iidx[i] - j + 1;
  for (k = j; k > i; k--, kl++)
    qbt1 += q[kl] * qq[k];

  return qbt1;
}


PRIVATE char
hc_default(int  i,
           int  j,
           int  k,
           int  l,
           char d,
           void *data)
{
  int                 kl, di, dj;
  char                eval;
  struct default_data *dat = (struct default_data *)data;

  eval  = (char)0;
  di    = k - i;
  dj    = j - l;
  switch (d) {
    case VRNA_DECOMP_EXT_EXT_STEM:
      kl = dat->idx[j] + l;
      if (dat->mx[kl] & VRNA_CONSTRAINT_CONTEXT_EXT_LOOP) {
        eval = (char)1;
        if (i != l) {
          /* otherwise, stem spans from i to j */
          di = l - k + 1;
          if ((di != 0) && (dat->hc_up[k + 1] < di))
            eval = (char)0;
        }
      }
      break;

    case VRNA_DECOMP_EXT_EXT_STEM1:
      kl = dat->idx[j - 1] + l;
      if (dat->mx[kl] & VRNA_CONSTRAINT_CONTEXT_EXT_LOOP) {
        eval = (char)1;
        if (i != l) {
          /* otherwise, stem spans from i to j - 1 */
          di = l - k + 1;
          if (dat->hc_up[j] == 0)
            eval = (char)0;
          if ((di != 0) && (dat->hc_up[k + 1] < di))
            eval = (char)0;
        }
      }
      break;

    case VRNA_DECOMP_EXT_STEM:
      kl = dat->idx[l] + k;
      if (dat->mx[kl] & VRNA_CONSTRAINT_CONTEXT_EXT_LOOP) {
        eval = (char)1;
        if ((di != 0) && (dat->hc_up[i] < di))
          eval = (char)0;
        if ((dj != 0) && (dat->hc_up[l + 1] < dj))
          eval = (char)0;
      }
      break;

    case VRNA_DECOMP_EXT_EXT:
      eval = (char)1;
      if ((di != 0) && (dat->hc_up[i] < di))
        eval = (char)0;
      if ((dj != 0) && (dat->hc_up[l + 1] < dj))
        eval = (char)0;
      break;

    case VRNA_DECOMP_EXT_UP:
      di    = j - i + 1;
      eval  = (dat->hc_up[i] >= di) ? (char)1 : (char)0;
      break;

    default:
      nrerror("wtf");
  }
  return eval;
}


PRIVATE char
hc_default_user(int   i,
                int   j,
                int   k,
                int   l,
                char  d,
                void  *data)
{
  char                eval;
  struct default_data *dat = (struct default_data *)data;

  eval  = hc_default(i, j, k, l, d, data);
  eval  = (dat->hc_f(i, j, k, l, d, dat->hc_dat)) ? eval : (char)0;

  return eval;
}
