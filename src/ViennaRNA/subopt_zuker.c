#include <stdlib.h>

#include "ViennaRNA/utils/basic.h"
#include "ViennaRNA/alphabet.h"
#include "ViennaRNA/loops/external.h"
#include "ViennaRNA/loops/internal.h"
#include "ViennaRNA/loops/multibranch.h"
#include "ViennaRNA/constraints/hard.h"
#include "ViennaRNA/subopt_zuker.h"

#define MAXSECTORS        500     /* dimension for a backtrack array */

struct ml_aux {
  int *ml;
  int *up;
};

PRIVATE void
compute_outside_int(vrna_fold_compound_t *fc,
                    unsigned int          l,
                    int                   *outside_c);

PRIVATE void
compute_outside_ml(vrna_fold_compound_t *fc,
                   unsigned int         l,
                   int                  *outside_c,
                   struct ml_aux        *ml_helpers);


PRIVATE struct ml_aux *
init_ml_helper(vrna_fold_compound_t *fc);


PRIVATE void
prepare_ml_helper(vrna_fold_compound_t  *fc,
                  unsigned int          l,
                  int                   *outside_c,
                  struct ml_aux         *ml_helper);


PRIVATE void
free_ml_helper(struct ml_aux *ml_helper);


PRIVATE char *
backtrack_enc(vrna_fold_compound_t  *fc,
              unsigned int          i,
              unsigned int          j,
              int                   e,
              int                   *outside_c);


PRIVATE char *
backtrack_ext(vrna_fold_compound_t  *fc,
              unsigned int          i,
              unsigned int          j,
              int                   e,
              int                   *f3);


PRIVATE int *
compute_f3(vrna_fold_compound_t  *fc);


PUBLIC vrna_subopt_solution_t *
vrna_subopt_zuker2(vrna_fold_compound_t *fc)
{
  char                    *s;
  short                   *S, *S1, s5, s3;
  unsigned int            k, l, n, turn, type, *sn;
  int                     e, kl, *c, *outside_c, *f5, *f3, *idx;
  struct ml_aux           *ml_helper;
  vrna_param_t            *P;
  vrna_md_t               *md;
  vrna_hc_t               *hc;
  vrna_sc_t               *sc;
  vrna_subopt_solution_t  *sol;

  sol = NULL;

  if (fc) {
    n     = fc->length;
    sn    = fc->strand_number;
    S     = fc->sequence_encoding2;
    S1    = fc->sequence_encoding;
    P     = fc->params;
    md    = &(P->model_details);
    turn  = fc->params->model_details.min_loop_size;
    idx   = fc->jindx;
    f5    = fc->matrices->f5;
    c     = fc->matrices->c;
    hc    = fc->hc;
    sc    = fc->sc;
    outside_c = (int *)vrna_alloc(sizeof(int) * ((n * (n + 1)) / 2 + 2));
    ml_helper = init_ml_helper(fc);

    f3        = compute_f3(fc);

    /* initialize outside matrix */
    for (k = 0; k < (n * (n + 1)) / 2 + 1; k++)
      outside_c[k] = INF;

    /* start filling outside_c */
    l = n;
    compute_outside_int(fc, l, outside_c);

    /* backtrack all structures with pairs (i, n) 1 <= i < n */
    for (k = l - turn - 1; k > 0; k--) {
      kl = idx[l] + k;
      if ((sn[k - 1] && sn[k]) &&
          (c[kl] != INF) &&
          (f5[k - 1] != INF)) {
        type = vrna_get_ptype_md(S[k], S[l], md);
        s5  = (k > 1) ? S1[k - 1] : -1;
        e   = c[kl] +
              f5[k - 1] +
              vrna_E_ext_stem(type, s5, -1, P);

        if (sc) {
          if (sc->f)
            e += sc->f(1, n, k - 1, k, VRNA_DECOMP_EXT_EXT_STEM, sc->data);
        }

        s = backtrack_ext(fc, k, n, e, f3);
      }
    }

    for (l = n - 1; l > turn + 1; l--) {
      compute_outside_int(fc, l, outside_c);
      compute_outside_ml(fc, l, outside_c, ml_helper);
      //compute_outside_ext(fc, l, outside_c, ext_helper);

      /* backtrack all structures with pairs (i, l) 1 <= l < n */
      for (k = l - turn - 1; k > 0; k--) {
        kl = idx[l] + k;

        if (c[kl] != INF) {
          e = INF;

          /* 1st case, (i,l) is not enclosed by another pair */
          if ((f5[k - 1] != INF) &&
              (f3[l + 1] != INF) &&
              (sn[k - 1] == sn[k]) &&
              (sn[l] == sn[l + 1])) {
            type = vrna_get_ptype_md(S[k], S[l], md);
            s5 = (k > 1) ? S1[k - 1] : -1;
            s3 = S1[l + 1];
            e = f5[k - 1] +
                f3[l + 1] +
                c[kl] +
                vrna_E_ext_stem(type, s5, s3, P);
          }

          /* other cases, (i,j) is enclosed by at least one other pair */
          if (outside_c[kl] != INF) {
            if (e < outside_c[kl] + c[kl])
              backtrack_ext(fc, k, l, e, f3);
            else
              backtrack_enc(fc, k, l, e, outside_c);
          }
        }
      }
    }

    free(f3);
    free(outside_c);
    free_ml_helper(ml_helper);
  }

  return sol;
}


PRIVATE struct ml_aux *
init_ml_helper(vrna_fold_compound_t *fc)
{
  unsigned int  i, n;
  struct ml_aux *ml_helper;

  ml_helper = (struct ml_aux *)vrna_alloc(sizeof(struct ml_aux));
  n         = fc->length;

  ml_helper->ml   = (int *)vrna_alloc(sizeof(int) * (n + 1));
  ml_helper->up   = (int *)vrna_alloc(sizeof(int) * (n + 1));

  return ml_helper;
}


PRIVATE void
prepare_ml_helpers(vrna_fold_compound_t  *fc,
                    unsigned int          l,
                    int                   *outside_c,
                    struct ml_aux         *ml_helper)
{
  short         *S, *S1;
  unsigned int  i, j, n, type, turn, u;
  int           ij, e, ee, *idx, *fML, dangle_model;
  vrna_param_t  *P;
  vrna_md_t     *md;
  vrna_hc_t     *hc;

  n     = fc->length;
  S     = fc->sequence_encoding2;
  S1    = fc->sequence_encoding;
  idx   = fc->jindx;
  P     = fc->params;
  md    = &(P->model_details);
  turn  = md->min_loop_size;
  dangle_model = md->dangles;
  hc    = fc->hc;
  fML   = fc->matrices->fML;

  /* initialize with INF */
  for (i = 0; i < l; i++) {
    ml_helper->ml[i]  = INF;
    ml_helper->up[i]  = INF;
  }

  for (j = l + 1; j <= n; j++) {
    u = j - l - 1;
    for (i = l - turn - 2; i > 0; i--) {
      if (hc->mx[j * n + i] & VRNA_CONSTRAINT_CONTEXT_MB_LOOP) {
        type = vrna_get_ptype_md(S[j], S[i], md);
        e = outside_c[idx[j] + i] +
            E_MLstem(type, S1[j - 1], S1[i + 1], P) +
            P->MLclosing;

        if (j > l + turn) {
          ee = e +
               fML[idx[j - 1] + l + 1];

          ml_helper->ml[i] = MIN2(ml_helper->ml[i], ee);
        }

        if (hc->up_ml[l + 1] >= u) {
          ee = e +
               u * P->MLbase;

          ml_helper->up[i] = MIN2(ml_helper->up[i], ee);
        }
      }
    }
  }
}


PRIVATE void
free_ml_helper(struct ml_aux *ml_helper)
{
  free(ml_helper->ml);
  free(ml_helper->up);
  free(ml_helper);
}


PRIVATE void
compute_outside_int(vrna_fold_compound_t *fc,
                    unsigned int          l,
                    int                   *outside_c)
{
  unsigned int  i, j, k, n, u1, u2, turn, min_i;
  int           ij, kl, *idx, *c, tmp;
  vrna_param_t  *P;
  vrna_md_t     *md;
  vrna_hc_t     *hc;

  n   = fc->length;
  idx = fc->jindx;
  P   = fc->params;
  md  = &(P->model_details);
  hc  = fc->hc;
  c   = fc->matrices->c;

  for(k = 2; k < l - turn; k++) {
    kl = idx[l] + k;

    if (c[kl] == INF)
      continue;

    if (hc->mx[l * n + k] & VRNA_CONSTRAINT_CONTEXT_INT_LOOP_ENC) {
      for (j = l + 1; j <= MIN2(l + MAXLOOP + 1, n); j++) {
        u2 = j - l - 1;

        if (hc->up_int[l + 1] < u2)
          break;

        min_i = (k > MAXLOOP + 1) ? k - MAXLOOP - 1 : 1;
        for (i = k - 1; i >= min_i; i--) {
          ij = idx[j] + i;

          if (c[ij] == INF)
            continue;

          u1 = k - i - 1;
          if (hc->up_int[i + 1] < u1)
            break;

          if (hc->mx[j * n + i] & VRNA_CONSTRAINT_CONTEXT_INT_LOOP) {
            tmp = outside_c[ij] +
                  vrna_eval_int_loop(fc, i, j, k, l);

            outside_c[ij] = MIN2(outside_c[ij], tmp);
          }
        }
      }
    }
  }

/*
  if (md->gquad)
    compute_gquad_outside_internal(fc, l);
*/
}

PRIVATE void
compute_outside_ml(vrna_fold_compound_t *fc,
                   unsigned int         l,
                   int                  *outside_c,
                   struct ml_aux        *ml_helpers)
{
  short         *S, *S1;
  unsigned int  i, j, k, n, u, turn, type;
  int           ij, kl, e, ppp, *fML, *idx;
  vrna_param_t  *P;
  vrna_md_t     *md;
  vrna_hc_t     *hc;

  n     = fc->length;
  S     = fc->sequence_encoding2;
  S1    = fc->sequence_encoding;
  idx   = fc->jindx;
  P     = fc->params;
  md    = &(P->model_details);
  hc    = fc->hc;
  fML   = fc->matrices->fML;
  turn  = md->min_loop_size;

  prepare_ml_helpers(fc, l, outside_c, ml_helpers);

  /* process all (k,l) enclosed as substem of multibranch loop enclosed by (i,j) */
  for (k = 2; k < l - turn; k++) {
    if (hc->mx[l * n + k] & VRNA_CONSTRAINT_CONTEXT_MB_LOOP_ENC) {
      kl    = idx[l] + k;
      type  = vrna_get_ptype_md(S[k], S[l], md);
      e     = E_MLstem(type, S1[k - 1], S1[l + 1], P);

      for (i = k - 1; i > 0; i--) {
        if (ml_helpers->ml[i] != INF) {
          u = k - i - 1;
          if (hc->up_ml[i + 1] >= u) {
            ppp = u * P->MLbase +
                  e +
                  ml_helpers->ml[i];

            outside_c[kl] = MIN2(outside_c[kl], ppp);
          }

          if (fML[idx[k - 1] + i + 1] != INF) {
            ppp = fML[idx[k - 1] + i + 1] +
                  e +
                  ml_helpers->ml[i];

            outside_c[kl] = MIN2(outside_c[kl], ppp);
          }
        }

        if ((ml_helpers->up[i] != INF) &&
            (fML[idx[k - 1] + i + 1] != INF)) {
          ppp = fML[idx[k - 1] + i + 1] +
                e +
                ml_helpers->up[i];

          outside_c[kl] = MIN2(outside_c[kl], ppp);
        }
      }
    }
  }
}

PRIVATE int *
compute_f3(vrna_fold_compound_t  *fc)
{
  short         *S, *S1;
  unsigned int  j, min_j, k, n, u, jk, turn, *sn, type;
  int           e, *c, *f3, *idx;
  vrna_param_t  *P;
  vrna_md_t     *md;
  vrna_hc_t     *hc;
  vrna_sc_t     *sc;

  n     = fc->length;
  S     = fc->sequence_encoding2;
  S1    = fc->sequence_encoding;
  sn    = fc->strand_number;
  P     = fc->params;
  md    = &(P->model_details);
  turn  = md->min_loop_size;
  c     = fc->matrices->c;
  hc    = fc->hc;
  sc    = fc->sc;

  /* init */
  f3 = (int *)vrna_alloc(sizeof(int) * (n + 2));

  f3[n + 1] = 0;

  min_j = (n > turn) ? n - turn : 1;

  f3[n] = INF;
  if (hc->up_ext[n]) {
    f3[n] = 0;
    if (sc) {
      if (sc->energy_up)
        f3[n] += sc->energy_up[n][1];

      if (sc->f)
        f3[n] += sc->f(n, n, n, n, VRNA_DECOMP_EXT_UP, sc->data);
    }
  }

  for (j = n - 1; j >= min_j; j--) {
    f3[j] = INF;
    u = n - j + 1;

    if ((hc->up_ext[j] >= u) &&
        (sn[j] == sn[j + 1])) {
      f3[j] = 0;

      if (sc) {
        if (sc->energy_up)
          f3[j] += sc->energy_up[j][1];

        if (sc->f)
          f3[j] += sc->f(j, n, j + 1, n, VRNA_DECOMP_EXT_EXT, sc->data);
      }
    }
  }

  for (j = min_j - 1; j > turn; j --) {
    /* 1st case, j is unpaired */
    if ((hc->up_ext[j]) &&
        (sn[j] == sn[j + 1])) {
      e     = f3[j + 1];

      if (sc) {
        if (sc->energy_up)
          e += sc->energy_up[j][1];

        if (sc->f)
          e += sc->f(j, n, j + 1, n, VRNA_DECOMP_EXT_EXT, sc->data);
      }

      f3[j] = MIN2(f3[j], e);
    }

    /* 2nd case, j forms pair (j,k) with j < k < n */
    for (k = j + turn + 1; k < n; k++) {
      if (hc->mx[j * n + k] & VRNA_CONSTRAINT_CONTEXT_EXT_LOOP) {
        jk = idx[k] + j;
        if ((c[jk] != INF) &&
            (f3[k + 1] != INF) &&
            (sn[k] == sn[k + 1])) {
          type = vrna_get_ptype_md(S[j], S[k], md);
          e = c[jk] +
              f3[k + 1] +
              vrna_E_ext_stem(type, S1[j - 1], S1[k + 1], P);

          if (sc) {
            if (sc->f) {
              e += sc->f(j, n, k, k + 1, VRNA_DECOMP_EXT_STEM_EXT, sc->data);
            }
          }

          f3[j] = MIN2(f3[j], e);
        }
      }
    }

    /* 3rd case, j forms pair with (j, n) */
    if (hc->mx[j * n + n] & VRNA_CONSTRAINT_CONTEXT_EXT_LOOP) {
      jk = idx[n] + j;
      if (c[jk] != INF) {
        e = c[jk] +
            vrna_E_ext_stem(type, S1[j - 1], -1, P);

        if (sc) {
          if (sc->f) {
            e += sc->f(j, n, n, k, VRNA_DECOMP_EXT_STEM, sc->data);
          }
        }

        f3[j] = MIN2(f3[j], e);
      }
    }

  }

  return f3;
}


PRIVATE char *
backtrack_enc(vrna_fold_compound_t  *fc,
              unsigned int          i,
              unsigned int          j,
              int                   e,
              int                   *outside_c)
{
  char *structure = NULL;
  
  
  return structure;
}


PRIVATE char *
backtrack_ext(vrna_fold_compound_t  *fc,
              unsigned int          i,
              unsigned int          j,
              int                   e,
              int                   *f3)
{
  char *structure = NULL;
  unsigned int n, k;
  int           s;
  sect              bt_stack[MAXSECTORS]; /* stack of partial structures for backtracking */
  vrna_bp_stack_t   *bp;

  n = fc->length;

  s = 0;

  /* add a guess of how many G's may be involved in a G quadruplex */
  bp = (vrna_bp_stack_t *)vrna_alloc(sizeof(vrna_bp_stack_t) * (4 * (1 + length / 2)));

  /* backtrack 3' part of external loop first */
  if (j < n) {
    stack_count = 0;
    k           = j + 1;
    while (k <= n) {
      if (backtrack_f3(fc, &k, &i, &j, bp, &stack_count, f3)) {
        if (i > 0) { /* store base pair for future backtracking */
          bt_stack[++s].i = i;
          bt_stack[s].j   = j;
          bt_stack[s].ml  = 2;
        }
      } else {
        vrna_message_warning("Backtracking failed in f3[%d] = %d", k, f3[k]);
      }

      if (k == 0)
        break;
    }
  }

  /* push external loop interval [1:i - 1] on bt_stack */
  if (i > 1) {
    bt_stack[++s].i = 1;
    bt_stack[s].j   = i - 1;
    bt_stack[s].ml  = 0;
  }
  /* push interval enclosed by (i,j) on bt_stack */
  bt_stack[++s].i = i;
  bt_stack[s].j   = j;
  bt_stack[s].ml  = 2;

  if (vrna_backtrack_from_intervals(fc, bp, bt_stack, s)) {

    k = j + 1;

    structure = vrna_db_from_bp_stack(bp, length);
  }

  free(bp);

  return structure;
}


PRIVATE int
backtrack_f3(vrna_fold_compound_t *fc,
             unsigned int         *k,
             unsigned int         *i,
             unsigned int         *j,
             vrna_bp_stack_t      *bp,
             int                  *stack_count,
             int                  *f3)
{
  unsigned int ii;

  ii = *k;

  /* nibble-off unpaired bases on 5' side */
  do {
    fij = f3[ii];
    fj  = INF;

    /* do hard constraints check */
    if (sn[ii] == sn[ii + 1]) {
      fj = f3[ii + 1];

      if (sc) {
        if (sc->energy_up)
          fj += sc->energy_up[ii][1];

        if (sc->f)
          fj += sc->f(ii, n, ii + 1, n, VRNA_DECOMP_EXT_EXT, sc->data);
      }
    }
    if (++ii > n)
      break;
  } while (fij == fj);
  ii--;

  if (ii > n - turn + 1) {
    /* no more pairs */
    *k = *i = *j = 0;
    return 1;
  }

  /* ii is paired, find pairing partner */
  switch (dangle_model) {
    case 0:
      for (u = ii + turn + 1; u <= n; u++) {
        if ((sn[u] == sn[u + 1]) &&
            (hc->mx[ii * n + u] & VRNA_CONSTRAINT_CONTEXT_EXT_LOOP)) {
          type  = vrna_get_ptype_md(S[ii], S[u], md);
          en    = c[idx[u] + ii];

          if (sc) {
            if (sc->f)
              en += sc->f(ii, n, u, u + 1, VRNA_DECOMP_EXT_STEM_EXT, sc->data);
          }

          if (fij == vrna_E_ext_stem(type, -1, -1, P) + en + f3[u + 1]) {
            *i = ii;
            *j = u;
            *k = u + 1;
            bp_stack[++(*stack_count)].i  = ii;
            bp_stack[(*stack_count)].j    = u;
            return 1;
          }
        }
      }
      break;

    case 2:
      for (u = ii + turn + 1; u <= n; u++) {
        if ((sn[u] == sn[u + 1]) &&
            (hc->mx[ii * n + u] & VRNA_CONSTRAINT_CONTEXT_EXT_LOOP)) {
          type  = vrna_get_ptype_md(S[ii], S[u], md);
          s5    = S1[ii - 1];
          s3    = (u < n) ? S1[u + 1] : -1;
          en    = c[idx[u] + ii];

          if (sc) {
            if (sc->f)
              en += sc->f(ii, n, u, u + 1, VRNA_DECOMP_EXT_STEM_EXT, sc->data);
          }

          if (fij == vrna_E_ext_stem(type, s5, s3, P) + en + f3[u + 1]) {
            *i = ii;
            *j = u;
            *k = u + 1;
            bp_stack[++(*stack_count)].i  = ii;
            bp_stack[(*stack_count)].j    = u;
            return 1;
          }
        }
      }
      break;

    default:
      break;
  }

  return 0;
}
