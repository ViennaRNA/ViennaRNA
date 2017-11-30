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
#include "exterior_loops_sc.inc"

/*
 #################################
 # PRIVATE FUNCTION DECLARATIONS #
 #################################
 */

PRIVATE int
E_ext_loop_5(vrna_fold_compound_t *vc);


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
  if (fc)
    return E_ext_loop_5(fc);

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
  unsigned int              strands, *sn;
  int                       ij, en, e, type, *idx;
  vrna_param_t              *P;
  vrna_md_t                 *md;
  vrna_sc_t                 *sc;
  vrna_callback_hc_evaluate *evaluate;
  struct default_data       hc_dat_local;

  strands           = vc->strands;
  sn                = vc->strand_number;
  S                 = vc->sequence_encoding;
  idx               = vc->jindx;
  ptype             = vc->ptype;
  P                 = vc->params;
  md                = &(P->model_details);
  hard_constraints  = vc->hc->matrix;
  sc                = vc->sc;
  evaluate          = prepare_hc_default(vc, &hc_dat_local);

  e     = INF;
  ij    = idx[j] + i;
  type  = get_pair_type(ij, ptype);

  if ((strands == 1) || (sn[i] == sn[j])) {
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


/*
 #####################################
 # BEGIN OF STATIC HELPER FUNCTIONS  #
 #####################################
 */

/*
 *  extend f5 by adding an unpaired nucleotide or an unstructured domain
 *  to the 3' end
 */
PRIVATE INLINE int
reduce_f5_up(vrna_fold_compound_t       *vc,
             int                        j,
             vrna_callback_hc_evaluate  *evaluate,
             struct default_data        *hc_dat_local,
             struct sc_wrapper_f5       *sc_wrapper)
{
  int                 u, k, e, en, *f5;
  vrna_ud_t           *domains_up;
  sc_f5_reduce_to_ext *sc_red_ext;

  f5          = vc->matrices->f5;
  domains_up  = vc->domains_up;
  e           = INF;

  sc_red_ext = sc_wrapper->red_ext;

  /* check for 3' extension with one unpaired nucleotide */
  if (f5[j - 1] != INF) {
    if (evaluate(1, j, 1, j - 1, VRNA_DECOMP_EXT_EXT, hc_dat_local)) {
      e = f5[j - 1];

      if (sc_red_ext)
        e += sc_red_ext(j, 1, j - 1, sc_wrapper);
    }
  }

  if ((domains_up) && (domains_up->energy_cb)) {
    for (k = 0; k < domains_up->uniq_motif_count; k++) {
      u = domains_up->uniq_motif_size[k];
      if ((j - u >= 0) && (f5[j - u] != INF)) {
        if (evaluate(1, j, 1, j - u, VRNA_DECOMP_EXT_EXT, hc_dat_local)) {
          en = f5[j - u] +
               domains_up->energy_cb(vc,
                                     j - u + 1,
                                     j,
                                     VRNA_UNSTRUCTURED_DOMAIN_EXT_LOOP | VRNA_UNSTRUCTURED_DOMAIN_MOTIF,
                                     domains_up->data);

          if (sc_red_ext)
            en += sc_red_ext(j, 1, j - u, sc_wrapper);

          e = MIN2(e, en);
        }
      }
    }
  }

  return e;
}


PRIVATE INLINE int *
get_stem_contributions_d0(vrna_fold_compound_t      *vc,
                          int                       j,
                          vrna_callback_hc_evaluate *evaluate,
                          struct default_data       *hc_dat_local,
                          struct sc_wrapper_f5      *sc_wrapper)
{
  char                    *ptype;
  short                   **S;
  unsigned int            s, n_seq;
  int                     i, ij, *indx, turn, type, *c, *stems;
  vrna_param_t            *P;
  vrna_md_t               *md;

  sc_f5_split_in_ext_stem *sc_spl_stem;
  sc_f5_reduce_to_stem    *sc_red_stem;

  stems = (int *)vrna_alloc(sizeof(int) * j);

  P     = vc->params;
  md    = &(P->model_details);
  indx  = vc->jindx;
  c     = vc->matrices->c;
  turn  = md->min_loop_size;
  ij    = indx[j] + j - turn - 1;

  sc_spl_stem = sc_wrapper->decomp_stem;
  sc_red_stem = sc_wrapper->red_stem;

  switch (vc->type) {
    case VRNA_FC_TYPE_SINGLE:
      ptype = vc->ptype;

      if (sc_spl_stem) {
        for (i = j - turn - 1; i > 1; i--, ij--) {
          stems[i] = INF;
          if ((c[ij] != INF) &&
              (evaluate(1, j, i - 1, i, VRNA_DECOMP_EXT_EXT_STEM, hc_dat_local))) {
            type      = get_pair_type(ij, ptype);
            stems[i]  = c[ij] +
                        E_ExtLoop(type, -1, -1, P) +
                        sc_spl_stem(j, i - 1, i, sc_wrapper);
          }
        }
      } else {
        for (i = j - turn - 1; i > 1; i--, ij--) {
          stems[i] = INF;
          if ((c[ij] != INF) &&
              (evaluate(1, j, i - 1, i, VRNA_DECOMP_EXT_EXT_STEM, hc_dat_local))) {
            type      = get_pair_type(ij, ptype);
            stems[i]  = c[ij] +
                        E_ExtLoop(type, -1, -1, P);
          }
        }
      }

      stems[1]  = INF;
      ij        = indx[j] + 1;

      if ((c[ij] != INF) && (evaluate(1, j, 1, j, VRNA_DECOMP_EXT_STEM, hc_dat_local))) {
        type = get_pair_type(ij, ptype);

        stems[1] = c[ij] +
                   E_ExtLoop(type, -1, -1, P);

        if (sc_red_stem)
          stems[1] += sc_red_stem(j, 1, j, sc_wrapper);
      }

      break;

    case VRNA_FC_TYPE_COMPARATIVE:
      n_seq = vc->n_seq;
      S     = vc->S;

      if (sc_spl_stem) {
        for (i = j - turn - 1; i > 1; i--, ij--) {
          stems[i] = INF;
          if ((c[ij] != INF) &&
              (evaluate(1, j, i - 1, i, VRNA_DECOMP_EXT_EXT_STEM, hc_dat_local))) {
            stems[i] = c[ij] +
                       sc_spl_stem(j, i - 1, i, sc_wrapper);
            for (s = 0; s < n_seq; s++) {
              type      = get_pair_type_md(S[s][i], S[s][j], md);
              stems[i]  += E_ExtLoop(type, -1, -1, P);
            }
          }
        }
      } else {
        for (i = j - turn - 1; i > 1; i--, ij--) {
          stems[i] = INF;
          if ((c[ij] != INF) &&
              (evaluate(1, j, i - 1, i, VRNA_DECOMP_EXT_EXT_STEM, hc_dat_local))) {
            stems[i] = c[ij];
            for (s = 0; s < n_seq; s++) {
              type      = get_pair_type_md(S[s][i], S[s][j], md);
              stems[i]  += E_ExtLoop(type, -1, -1, P);
            }
          }
        }
      }

      stems[1]  = INF;
      ij        = indx[j] + 1;

      if ((c[ij] != INF) && (evaluate(1, j, 1, j, VRNA_DECOMP_EXT_STEM, hc_dat_local))) {
        stems[1] = c[ij];

        for (s = 0; s < n_seq; s++) {
          type      = get_pair_type_md(S[s][1], S[s][j], md);
          stems[1]  += E_ExtLoop(type, -1, -1, P);
        }

        if (sc_red_stem)
          stems[1] += sc_red_stem(j, 1, j, sc_wrapper);
      }

      break;
  }

  return stems;
}


PRIVATE INLINE int *
get_stem_contributions_d2(vrna_fold_compound_t      *vc,
                          int                       j,
                          vrna_callback_hc_evaluate *evaluate,
                          struct default_data       *hc_dat_local,
                          struct sc_wrapper_f5      *sc_wrapper)
{
  char                    *ptype;
  short                   *S, sj1, *si1, **SS, **S5, **S3, *s3j;
  unsigned int            s, n_seq, **a2s;
  int                     n, i, ij, *indx, turn, type, *c, *stems, mm5, mm3;
  vrna_param_t            *P;
  vrna_md_t               *md;

  sc_f5_split_in_ext_stem *sc_spl_stem;
  sc_f5_reduce_to_stem    *sc_red_stem;

  stems = (int *)vrna_alloc(sizeof(int) * j);

  n     = (int)vc->length;
  P     = vc->params;
  md    = &(P->model_details);
  indx  = vc->jindx;
  c     = vc->matrices->c;
  turn  = md->min_loop_size;
  ij    = indx[j] + j - turn - 1;

  sc_spl_stem = sc_wrapper->decomp_stem;
  sc_red_stem = sc_wrapper->red_stem;

  switch (vc->type) {
    case VRNA_FC_TYPE_SINGLE:
      S     = vc->sequence_encoding;
      ptype = vc->ptype;
      si1   = S + j - turn - 2;
      sj1   = (j < n) ? S[j + 1] : -1;

      if (sc_spl_stem) {
        for (i = j - turn - 1; i > 1; i--, ij--, si1--) {
          stems[i] = INF;
          if ((c[ij] != INF) &&
              (evaluate(1, j, i - 1, i, VRNA_DECOMP_EXT_EXT_STEM, hc_dat_local))) {
            type      = get_pair_type(ij, ptype);
            stems[i]  = c[ij] +
                        E_ExtLoop(type, *si1, sj1, P) +
                        sc_spl_stem(j, i - 1, i, sc_wrapper);
          }
        }
      } else {
        for (i = j - turn - 1; i > 1; i--, ij--, si1--) {
          stems[i] = INF;
          if ((c[ij] != INF) &&
              (evaluate(1, j, i - 1, i, VRNA_DECOMP_EXT_EXT_STEM, hc_dat_local))) {
            type      = get_pair_type(ij, ptype);
            stems[i]  = c[ij] +
                        E_ExtLoop(type, *si1, sj1, P);
          }
        }
      }

      stems[1]  = INF;
      ij        = indx[j] + 1;

      if ((c[ij] != INF) && (evaluate(1, j, 1, j, VRNA_DECOMP_EXT_STEM, hc_dat_local))) {
        type      = get_pair_type(ij, ptype);
        stems[1]  = c[ij] +
                    E_ExtLoop(type, -1, sj1, P);

        if (sc_red_stem)
          stems[1] += sc_red_stem(j, 1, j, sc_wrapper);
      }

      break;

    case VRNA_FC_TYPE_COMPARATIVE:
      n_seq = vc->n_seq;
      SS    = vc->S;
      S5    = vc->S5;
      S3    = vc->S3;
      a2s   = vc->a2s;

      /* pre-compute S3[s][j - 1] */
      s3j  = (short *)vrna_alloc(sizeof(short) * n_seq);
      for (s = 0; s < n_seq; s++)
        s3j[s] = (a2s[s][j] < a2s[s][SS[0][0]]) ? S3[s][j] : -1;

      if (sc_spl_stem) {
        for (i = j - turn - 1; i > 1; i--, ij--) {
          stems[i] = INF;
          if ((c[ij] != INF) &&
              (evaluate(1, j, i - 1, i, VRNA_DECOMP_EXT_EXT_STEM, hc_dat_local))) {
            stems[i] = c[ij] +
                       sc_spl_stem(j, i - 1, i, sc_wrapper);
            for (s = 0; s < n_seq; s++) {
              type      = get_pair_type_md(SS[s][i], SS[s][j], md);
              mm5       = (a2s[s][i] > 1) ? S5[s][i] : -1;
              stems[i]  += E_ExtLoop(type, mm5, s3j[s], P);
            }
          }
        }
      } else {
        for (i = j - turn - 1; i > 1; i--, ij--) {
          stems[i] = INF;
          if ((c[ij] != INF) &&
              (evaluate(1, j, i - 1, i, VRNA_DECOMP_EXT_EXT_STEM, hc_dat_local))) {
            stems[i] = c[ij];
            for (s = 0; s < n_seq; s++) {
              type      = get_pair_type_md(SS[s][i], SS[s][j], md);
              mm5       = (a2s[s][i] > 1) ? S5[s][i] : -1;
              stems[i]  += E_ExtLoop(type, mm5, s3j[s], P);
            }
          }
        }
      }

      stems[1]  = INF;
      ij        = indx[j] + 1;

      if ((c[ij] != INF) && (evaluate(1, j, 1, j, VRNA_DECOMP_EXT_STEM, hc_dat_local))) {
        stems[1] = c[ij];

        for (s = 0; s < n_seq; s++) {
          type      = get_pair_type_md(SS[s][1], SS[s][j], md);
          stems[1]  += E_ExtLoop(type, -1, s3j[s], P);
        }

        if (sc_red_stem)
          stems[1] += sc_red_stem(j, 1, j, sc_wrapper);
      }

      free(s3j);

      break;
  }

  return stems;
}


PRIVATE INLINE int *
get_stem_contributions_5(vrna_fold_compound_t       *vc,
                         int                        j,
                         vrna_callback_hc_evaluate  *evaluate,
                         struct default_data        *hc_dat_local,
                         struct sc_wrapper_f5       *sc_wrapper)
{
  char                    *ptype;
  short                   *S, *si1, **SS, **S5, **S3;
  unsigned int            s, n_seq, **a2s;
  int                     n, i, ij, *indx, turn, type, *c, *stems, mm5;
  vrna_param_t            *P;
  vrna_md_t               *md;

  sc_f5_split_in_ext_stem *sc_spl_stem;
  sc_f5_reduce_to_stem    *sc_red_stem;

  stems = (int *)vrna_alloc(sizeof(int) * j);

  n     = (int)vc->length;
  P     = vc->params;
  md    = &(P->model_details);
  indx  = vc->jindx;
  c     = vc->matrices->c;
  turn  = md->min_loop_size;
  ij    = indx[j] + j - turn;

  sc_spl_stem = sc_wrapper->decomp_stem;
  sc_red_stem = sc_wrapper->red_stem;

  switch (vc->type) {
    case VRNA_FC_TYPE_SINGLE:
      S     = vc->sequence_encoding;
      ptype = vc->ptype;
      si1   = S + j - turn - 1;

      if (sc_spl_stem) {
        for (i = j - turn - 1; i > 1; i--, ij--, si1--) {
          stems[i] = INF;
          if ((c[ij] != INF) &&
              (evaluate(1, j, i - 1, i + 1, VRNA_DECOMP_EXT_EXT_STEM, hc_dat_local))) {
            type      = get_pair_type(ij, ptype);
            stems[i]  = c[ij] +
                        E_ExtLoop(type, *si1, -1, P) +
                        sc_spl_stem(j, i - 1, i + 1, sc_wrapper);
          }
        }
      } else {
        for (i = j - turn - 1; i > 1; i--, ij--, si1--) {
          stems[i] = INF;
          if ((c[ij] != INF) &&
              (evaluate(1, j, i - 1, i + 1, VRNA_DECOMP_EXT_EXT_STEM, hc_dat_local))) {
            type      = get_pair_type(ij, ptype);
            stems[i]  = c[ij] +
                        E_ExtLoop(type, *si1, -1, P);
          }
        }
      }

      stems[1]  = INF;
      ij        = indx[j] + 2;

      if ((c[ij] != INF) && (evaluate(1, j, 2, j, VRNA_DECOMP_EXT_STEM, hc_dat_local))) {
        type      = get_pair_type(ij, ptype);
        stems[1]  = c[ij] +
                    E_ExtLoop(type, S[1], -1, P);

        if (sc_red_stem)
          stems[1] += sc_red_stem(j, 2, j, sc_wrapper);
      }

      break;

    case VRNA_FC_TYPE_COMPARATIVE:
      n_seq = vc->n_seq;
      SS    = vc->S;
      S5    = vc->S5;
      S3    = vc->S3;
      a2s   = vc->a2s;
      if (sc_spl_stem) {
        for (i = j - turn - 1; i > 1; i--, ij--) {
          stems[i] = INF;
          if ((c[ij] != INF) &&
              (evaluate(1, j, i - 1, i + 1, VRNA_DECOMP_EXT_EXT_STEM, hc_dat_local))) {
            stems[i] = c[ij] +
                       sc_spl_stem(j, i - 1, i + 1, sc_wrapper);
            for (s = 0; s < n_seq; s++) {
              type      = get_pair_type_md(SS[s][i + 1], SS[s][j], md);
              mm5       = (a2s[s][i + 1] > 1) ? S5[s][i + 1] : -1;
              stems[i]  = E_ExtLoop(type, mm5, -1, P);
            }
          }
        }
      } else {
        for (i = j - turn - 1; i > 1; i--, ij--) {
          stems[i] = INF;
          if ((c[ij] != INF) &&
              (evaluate(1, j, i - 1, i + 1, VRNA_DECOMP_EXT_EXT_STEM, hc_dat_local))) {
            stems[i] = c[ij];
            for (s = 0; s < n_seq; s++) {
              type      = get_pair_type_md(SS[s][i + 1], SS[s][j], md);
              mm5       = (a2s[s][i + 1] > 1) ? S5[s][i + 1] : -1;
              stems[i]  = E_ExtLoop(type, mm5, -1, P);
            }
          }
        }
      }

      stems[1]  = INF;
      ij        = indx[j] + 2;

      if ((c[ij] != INF) && (evaluate(1, j, 2, j, VRNA_DECOMP_EXT_STEM, hc_dat_local))) {
        stems[1] = c[ij];
        for (s = 0; s < n_seq; s++) {
          type      = get_pair_type_md(SS[s][2], SS[s][j], md);
          mm5       = (a2s[s][2] > 1) ? S5[s][2] : -1;
          stems[i]  = E_ExtLoop(type, mm5, -1, P);
        }

        if (sc_red_stem)
          stems[1] += sc_red_stem(j, 2, j, sc_wrapper);
      }

      break;
  }

  return stems;
}


PRIVATE INLINE int *
get_stem_contributions_3(vrna_fold_compound_t       *vc,
                         int                        j,
                         vrna_callback_hc_evaluate  *evaluate,
                         struct default_data        *hc_dat_local,
                         struct sc_wrapper_f5       *sc_wrapper)
{
  char                    *ptype;
  short                   *S, sj1, **SS, **S3, *s3j1;
  unsigned int            s, n_seq, **a2s;
  int                     n, i, ij, *indx, turn, type, *c, *stems;
  vrna_param_t            *P;
  vrna_md_t               *md;

  sc_f5_split_in_ext_stem *sc_spl_stem;
  sc_f5_reduce_to_stem    *sc_red_stem;

  stems = (int *)vrna_alloc(sizeof(int) * j);

  n     = (int)vc->length;
  P     = vc->params;
  md    = &(P->model_details);
  indx  = vc->jindx;
  c     = vc->matrices->c;
  turn  = P->model_details.min_loop_size;
  ij    = indx[j - 1] + j - turn - 1;

  sc_spl_stem = sc_wrapper->decomp_stem1;
  sc_red_stem = sc_wrapper->red_stem;

  switch (vc->type) {
    case VRNA_FC_TYPE_SINGLE:
      S     = vc->sequence_encoding;
      ptype = vc->ptype;
      sj1   = S[j];

      if (sc_spl_stem) {
        for (i = j - turn - 1; i > 1; i--, ij--) {
          stems[i] = INF;
          if ((c[ij] != INF) &&
              (evaluate(1, j, i - 1, i, VRNA_DECOMP_EXT_EXT_STEM1, hc_dat_local))) {
            type      = get_pair_type(ij, ptype);
            stems[i]  = c[ij] +
                        E_ExtLoop(type, -1, sj1, P) +
                        sc_spl_stem(j, i - 1, i, sc_wrapper);
          }
        }
      } else {
        for (i = j - turn - 1; i > 1; i--, ij--) {
          stems[i] = INF;
          if ((c[ij] != INF) &&
              (evaluate(1, j, i - 1, i, VRNA_DECOMP_EXT_EXT_STEM1, hc_dat_local))) {
            type      = get_pair_type(ij, ptype);
            stems[i]  = c[ij] +
                        E_ExtLoop(type, -1, sj1, P);
          }
        }
      }

      stems[1]  = INF;
      ij        = indx[j - 1] + 1;

      if ((c[ij] != INF) && (evaluate(1, j, 1, j - 1, VRNA_DECOMP_EXT_STEM, hc_dat_local))) {
        type      = get_pair_type(ij, ptype);
        stems[1]  = c[ij] +
                    E_ExtLoop(type, -1, sj1, P);

        if (sc_red_stem)
          stems[1] += sc_red_stem(j, 1, j - 1, sc_wrapper);
      }

      break;

    case VRNA_FC_TYPE_COMPARATIVE:
      n_seq = vc->n_seq;
      SS    = vc->S;
      S3    = vc->S3;
      a2s   = vc->a2s;

      /* pre-compute S3[s][j - 1] */
      s3j1  = (short *)vrna_alloc(sizeof(short) * n_seq);
      for (s = 0; s < n_seq; s++)
        s3j1[s] = (a2s[s][j - 1] < a2s[s][SS[0][0]]) ? S3[s][j - 1] : -1;

      if (sc_spl_stem) {
        for (i = j - turn - 1; i > 1; i--, ij--) {
          stems[i] = INF;
          if ((c[ij] != INF) &&
              (evaluate(1, j, i - 1, i, VRNA_DECOMP_EXT_EXT_STEM1, hc_dat_local))) {
            stems[i] = c[ij] +
                       sc_spl_stem(j, i - 1, i, sc_wrapper);
            for (s = 0; s < n_seq; s++) {
              type      = get_pair_type_md(SS[s][i], SS[s][j - 1], md);
              stems[i]  += E_ExtLoop(type, -1, s3j1[s], P);
            }
          }
        }
      } else {
        for (i = j - turn - 1; i > 1; i--, ij--) {
          stems[i] = INF;
          if ((c[ij] != INF) &&
              (evaluate(1, j, i - 1, i, VRNA_DECOMP_EXT_EXT_STEM1, hc_dat_local))) {
            stems[i] = c[ij];
            for (s = 0; s < n_seq; s++) {
              type      = get_pair_type_md(SS[s][i], SS[s][j - 1], md);
              stems[i]  += E_ExtLoop(type, -1, s3j1[s], P);
            }
          }
        }
      }

      stems[1]  = INF;
      ij        = indx[j - 1] + 1;

      if ((c[ij] != INF) && (evaluate(1, j, 1, j - 1, VRNA_DECOMP_EXT_STEM, hc_dat_local))) {
        stems[1] = c[ij];

        for (s = 0; s < n_seq; s++) {
          type      = get_pair_type_md(SS[s][1], SS[s][j - 1], md);
          stems[1]  += E_ExtLoop(type, -1, s3j1[s], P);
        }

        if (sc_red_stem)
          stems[1] += sc_red_stem(j, 1, j - 1, sc_wrapper);
      }

      free(s3j1);

      break;
  }

  return stems;
}


PRIVATE INLINE int *
get_stem_contributions_53(vrna_fold_compound_t      *vc,
                          int                       j,
                          vrna_callback_hc_evaluate *evaluate,
                          struct default_data       *hc_dat_local,
                          struct sc_wrapper_f5      *sc_wrapper)
{
  char                    *ptype;
  short                   *S, *si1, sj1, **SS, **S5, **S3, *s3j1;
  unsigned int            s, n_seq, **a2s;
  int                     n, i, ij, *indx, turn, type, *c, *stems, mm5, mm3;
  vrna_param_t            *P;
  vrna_md_t               *md;

  sc_f5_split_in_ext_stem *sc_spl_stem;
  sc_f5_reduce_to_stem    *sc_red_stem;

  stems = (int *)vrna_alloc(sizeof(int) * j);

  n     = (int)vc->length;
  P     = vc->params;
  md    = &(P->model_details);
  indx  = vc->jindx;
  c     = vc->matrices->c;
  turn  = md->min_loop_size;
  ij    = indx[j - 1] + j - turn;

  sc_spl_stem = sc_wrapper->decomp_stem1;
  sc_red_stem = sc_wrapper->red_stem;

  switch (vc->type) {
    case VRNA_FC_TYPE_SINGLE:
      S     = vc->sequence_encoding;
      ptype = vc->ptype;
      sj1   = S[j];
      si1   = S + j - turn - 1;

      if (sc_spl_stem) {
        for (i = j - turn - 1; i > 1; i--, ij--, si1--) {
          stems[i] = INF;
          if ((c[ij] != INF) &&
              (evaluate(1, j, i - 1, i + 1, VRNA_DECOMP_EXT_EXT_STEM1, hc_dat_local))) {
            type      = get_pair_type(ij, ptype);
            stems[i]  = c[ij] +
                        E_ExtLoop(type, *si1, sj1, P) +
                        sc_spl_stem(j, i - 1, i + 1, sc_wrapper);
          }
        }
      } else {
        for (i = j - turn - 1; i > 1; i--, ij--, si1--) {
          stems[i] = INF;
          if ((c[ij] != INF) &&
              (evaluate(1, j, i - 1, i + 1, VRNA_DECOMP_EXT_EXT_STEM1, hc_dat_local))) {
            type      = get_pair_type(ij, ptype);
            stems[i]  = c[ij] +
                        E_ExtLoop(type, *si1, sj1, P);
          }
        }
      }

      stems[1]  = INF;
      ij        = indx[j - 1] + 2;

      if ((c[ij] != INF) && (evaluate(1, j, 2, j - 1, VRNA_DECOMP_EXT_STEM, hc_dat_local))) {
        type      = get_pair_type(ij, ptype);
        stems[1]  = c[ij] +
                    E_ExtLoop(type, S[1], sj1, P);

        if (sc_red_stem)
          stems[1] += sc_red_stem(j, 2, j - 1, sc_wrapper);
      }

      break;

    case VRNA_FC_TYPE_COMPARATIVE:
      n_seq = vc->n_seq;
      SS    = vc->S;
      S5    = vc->S5;
      S3    = vc->S3;
      a2s   = vc->a2s;

      /* pre-compute S3[s][j - 1] */
      s3j1  = (short *)vrna_alloc(sizeof(short) * n_seq);
      for (s = 0; s < n_seq; s++)
        s3j1[s] = (a2s[s][j - 1] < a2s[s][SS[0][0]]) ? S3[s][j - 1] : -1;

      if (sc_spl_stem) {
        for (i = j - turn - 1; i > 1; i--, ij--) {
          stems[i] = INF;
          if ((c[ij] != INF) &&
              (evaluate(1, j, i - 1, i + 1, VRNA_DECOMP_EXT_EXT_STEM1, hc_dat_local))) {
            stems[i] = c[ij] +
                       sc_spl_stem(j, i - 1, i + 1, sc_wrapper);
            for (s = 0; s < n_seq; s++) {
              type      = get_pair_type_md(SS[s][i + 1], SS[s][j - 1], md);
              stems[i]  += E_ExtLoop(type, (a2s[s][i + 1] > 1) ? S5[s][i + 1] : -1, s3j1[s], P);
            }
          }
        }
      } else {
        for (i = j - turn - 1; i > 1; i--, ij--, si1--) {
          stems[i] = INF;
          if ((c[ij] != INF) &&
              (evaluate(1, j, i - 1, i + 1, VRNA_DECOMP_EXT_EXT_STEM1, hc_dat_local))) {
            stems[i] = c[ij];
            for (s = 0; s < n_seq; s++) {
              type      = get_pair_type_md(SS[s][i + 1], SS[s][j - 1], md);
              stems[i]  += E_ExtLoop(type, (a2s[s][i + 1] > 1) ? S5[s][i + 1] : -1, s3j1[s], P);
            }
          }
        }
      }

      stems[1]  = INF;
      ij        = indx[j - 1] + 2;

      if ((c[ij] != INF) && (evaluate(1, j, 2, j - 1, VRNA_DECOMP_EXT_STEM, hc_dat_local))) {
        stems[1] = c[ij];
        for (s = 0; s < n_seq; s++) {
          type      = get_pair_type_md(SS[s][2], SS[s][j - 1], md);
          stems[1]  += E_ExtLoop(type, (a2s[s][2] > 1) ? S5[s][2] : -1, s3j1[s], P);
        }

        if (sc_red_stem)
          stems[1] += sc_red_stem(j, 2, j - 1, sc_wrapper);
      }

      free(s3j1);

      break;
  }

  return stems;
}


#ifdef VRNA_WITH_SSE_IMPLEMENTATION
/* SSE modular decomposition -------------------------------*/
#include <emmintrin.h>
#include <smmintrin.h>

PRIVATE INLINE int
horizontal_min_Vec4i(__m128i x)
{
  __m128i min1  = _mm_shuffle_epi32(x, _MM_SHUFFLE(0, 0, 3, 2));
  __m128i min2  = _mm_min_epi32(x, min1);
  __m128i min3  = _mm_shuffle_epi32(min2, _MM_SHUFFLE(0, 0, 0, 1));
  __m128i min4  = _mm_min_epi32(min2, min3);

  return _mm_cvtsi128_si32(min4);
}


PRIVATE INLINE int
modular_decomposition(const int j,
                      const int turn,
                      const int *f5,
                      const int *stems)
{
  int       i, decomp = INF;
  __m128i   inf = _mm_set1_epi32(INF);

  const int end = j - turn;

  for (i = 2; i < end - 3; i += 4) {
    __m128i   a     = _mm_loadu_si128((__m128i *)&f5[i - 1]);
    __m128i   b     = _mm_loadu_si128((__m128i *)&stems[i]);
    __m128i   c     = _mm_add_epi32(a, b);
    __m128i   mask1 = _mm_cmplt_epi32(a, inf);
    __m128i   mask2 = _mm_cmplt_epi32(b, inf);
    __m128i   res   = _mm_or_si128(_mm_and_si128(mask1, c),
                                   _mm_andnot_si128(mask1, a));

    res = _mm_or_si128(_mm_and_si128(mask2, res),
                       _mm_andnot_si128(mask2, b));
    const int en = horizontal_min_Vec4i(res);
    decomp = MIN2(decomp, en);
  }

  for (; i < end; i++) {
    if ((f5[i - 1] != INF) && (stems[i] != INF)) {
      const int en = f5[i - 1] + stems[i];
      decomp = MIN2(decomp, en);
    }
  }

  return decomp;
}


#endif

PRIVATE INLINE int
decompose_f5_ext_stem(vrna_fold_compound_t  *vc,
                      int                   j,
                      int                   *stems)
{
  int e, i, *f5, turn;

  f5    = vc->matrices->f5;
  turn  = vc->params->model_details.min_loop_size;
  e     = INF;

  /* modular decomposition */
#if VRNA_WITH_SSE_IMPLEMENTATION
  e = modular_decomposition(j, turn, f5, stems);
#else
  for (i = 2; i < j - turn; i++)
    if ((f5[i - 1] != INF) && (stems[i] != INF)) {
      const int en = f5[i - 1] + stems[i];
      e = MIN2(e, en);
    }

#endif

  return e;
}


PRIVATE INLINE int
decompose_f5_ext_stem_d0(vrna_fold_compound_t       *vc,
                         int                        j,
                         vrna_callback_hc_evaluate  *evaluate,
                         struct default_data        *hc_dat_local,
                         struct sc_wrapper_f5       *sc_wrapper)
{
  int e, *stems;

  stems = get_stem_contributions_d0(vc, j, evaluate, hc_dat_local, sc_wrapper);

  /* 1st case, actual decompostion */
  e = decompose_f5_ext_stem(vc, j, stems);

  /* 2nd case, reduce to single stem */
  e = MIN2(e, stems[1]);

  free(stems);

  return e;
}


PRIVATE INLINE int
decompose_f5_ext_stem_d2(vrna_fold_compound_t       *vc,
                         int                        j,
                         vrna_callback_hc_evaluate  *evaluate,
                         struct default_data        *hc_dat_local,
                         struct sc_wrapper_f5       *sc_wrapper)
{
  int e, *stems;

  stems = get_stem_contributions_d2(vc, j, evaluate, hc_dat_local, sc_wrapper);

  /* 1st case, actual decompostion */
  e = decompose_f5_ext_stem(vc, j, stems);

  /* 2nd case, reduce to single stem */
  e = MIN2(e, stems[1]);

  free(stems);

  return e;
}


PRIVATE INLINE int
decompose_f5_ext_stem_d1(vrna_fold_compound_t       *vc,
                         int                        j,
                         vrna_callback_hc_evaluate  *evaluate,
                         struct default_data        *hc_dat_local,
                         struct sc_wrapper_f5       *sc_wrapper)
{
  int e, ee, *stems;

  e = INF;

  /* A) without dangling end contributions */

  /* 1st case, actual decompostion */
  stems = get_stem_contributions_d0(vc, j, evaluate, hc_dat_local, sc_wrapper);

  ee = decompose_f5_ext_stem(vc, j, stems);

  /* 2nd case, reduce to single stem */
  ee = MIN2(ee, stems[1]);

  free(stems);

  e = MIN2(e, ee);

  /* B) with dangling end contribution on 5' side of stem */
  stems = get_stem_contributions_5(vc, j, evaluate, hc_dat_local, sc_wrapper);

  /* 1st case, actual decompostion */
  ee = decompose_f5_ext_stem(vc, j, stems);

  /* 2nd case, reduce to single stem */
  ee = MIN2(ee, stems[1]);

  free(stems);

  e = MIN2(e, ee);

  /* C) with dangling end contribution on 3' side of stem */
  stems = get_stem_contributions_3(vc, j, evaluate, hc_dat_local, sc_wrapper);

  /* 1st case, actual decompostion */
  ee = decompose_f5_ext_stem(vc, j, stems);

  /* 2nd case, reduce to single stem */
  ee = MIN2(ee, stems[1]);

  free(stems);

  e = MIN2(e, ee);

  /* D) with dangling end contribution on both sides of stem */
  stems = get_stem_contributions_53(vc, j, evaluate, hc_dat_local, sc_wrapper);

  /* 1st case, actual decompostion */
  ee = decompose_f5_ext_stem(vc, j, stems);

  /* 2nd case, reduce to single stem */
  ee = MIN2(ee, stems[1]);

  free(stems);

  e = MIN2(e, ee);

  return e;
}


PRIVATE INLINE int
add_f5_gquad(vrna_fold_compound_t       *vc,
             int                        j,
             vrna_callback_hc_evaluate  *evaluate,
             struct default_data        *hc_dat_local,
             struct sc_wrapper_f5       *sc_wrapper)
{
  int e, i, ij, *indx, turn, *f5, *ggg;

  indx  = vc->jindx;
  f5    = vc->matrices->f5;
  ggg   = vc->matrices->ggg;
  turn  = vc->params->model_details.min_loop_size;
  ij    = indx[j] + j - turn - 1;
  e     = INF;

  for (i = j - turn - 1; i > 1; i--, ij--)
    if ((f5[i - 1] != INF) && (ggg[ij] != INF))
      e = MIN2(e, f5[i - 1] + ggg[ij]);

  ij  = indx[j] + 1;
  e   = MIN2(e, ggg[ij]);

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
  struct sc_wrapper_f5      sc_wrapper;

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
  evaluate      = prepare_hc_default(vc, &hc_dat_local);

  init_sc_wrapper(vc, &sc_wrapper);

  f5[0] = 0;
  for (j = 1; j <= turn + 1; j++)
    f5[j] = reduce_f5_up(vc, j, evaluate, &hc_dat_local, &sc_wrapper);

  /*
   *  duplicated code may be faster than conditions inside loop or even
   *  using a function pointer ;)
   */
  switch (dangle_model) {
    case 2:
      for (j = turn + 2; j <= length; j++) {
        /* extend previous solution(s) by adding an unpaired region */
        f5[j] = reduce_f5_up(vc, j, evaluate, &hc_dat_local, &sc_wrapper);

        /* decompose into exterior loop part followed by a stem */
        en    = decompose_f5_ext_stem_d2(vc, j, evaluate, &hc_dat_local, &sc_wrapper);
        f5[j] = MIN2(f5[j], en);

        if (with_gquad) {
          en    = add_f5_gquad(vc, j, evaluate, &hc_dat_local, &sc_wrapper);
          f5[j] = MIN2(f5[j], en);
        }
      }
      break;

    case 0:
      for (j = turn + 2; j <= length; j++) {
        /* extend previous solution(s) by adding an unpaired region */
        f5[j] = reduce_f5_up(vc, j, evaluate, &hc_dat_local, &sc_wrapper);

        /* decompose into exterior loop part followed by a stem */
        en    = decompose_f5_ext_stem_d0(vc, j, evaluate, &hc_dat_local, &sc_wrapper);
        f5[j] = MIN2(f5[j], en);

        if (with_gquad) {
          en    = add_f5_gquad(vc, j, evaluate, &hc_dat_local, &sc_wrapper);
          f5[j] = MIN2(f5[j], en);
        }
      }
      break;

    default:
      for (j = turn + 2; j <= length; j++) {
        /* extend previous solution(s) by adding an unpaired region */
        f5[j] = reduce_f5_up(vc, j, evaluate, &hc_dat_local, &sc_wrapper);

        en    = decompose_f5_ext_stem_d1(vc, j, evaluate, &hc_dat_local, &sc_wrapper);
        f5[j] = MIN2(f5[j], en);

        if (with_gquad) {
          en    = add_f5_gquad(vc, j, evaluate, &hc_dat_local, &sc_wrapper);
          f5[j] = MIN2(f5[j], en);
        }
      }
      break;
  }

  free_sc_wrapper(&sc_wrapper);

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
  evaluate      = prepare_hc_default_window(fc, &hc_dat_local);

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
  evaluate      = prepare_hc_default_window(fc, &hc_dat_local);

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
