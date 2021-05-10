#include <stdlib.h>

#include "ViennaRNA/utils/basic.h"
#include "ViennaRNA/alphabet.h"
#include "ViennaRNA/loops/external.h"
#include "ViennaRNA/loops/internal.h"
#include "ViennaRNA/loops/multibranch.h"
#include "ViennaRNA/constraints/hard.h"
#include "ViennaRNA/mfe.h"
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


PRIVATE int *
compute_f3(vrna_fold_compound_t  *fc);


PRIVATE char *
backtrack_pair(vrna_fold_compound_t *fc);


PRIVATE char *
backtrack_ext5(vrna_fold_compound_t *fc,
               unsigned int         k);


PRIVATE char *
backtrack_ext3(vrna_fold_compound_t  *fc,
               unsigned int          k,
               int                   e,
               int                   *f3);


PRIVATE char *
backtrack_ext53(vrna_fold_compound_t *fc,
                unsigned int          k,
                unsigned int          l,
                int                   e,
                int                   *f3);


PRIVATE int
backtrack_f3(vrna_fold_compound_t *fc,
             unsigned int         *k,
             unsigned int         *i,
             unsigned int         *j,
             int                  *f3);


PRIVATE char *
backtrack_int(vrna_fold_compound_t  *fc,
              unsigned int          k,
              unsigned int          l,
              int                   e,
              int                   *outside_c);


PRIVATE char *
backtrack_mb(vrna_fold_compound_t  *fc,
             unsigned int          k,
             unsigned int          l,
             int                   e,
             int                   *outside_c);


PUBLIC vrna_subopt_solution_t *
vrna_subopt_zuker2(vrna_fold_compound_t *fc)
{
  char                    *s, *mfe_structure;
  short                   *S, *S1, s5, s3;
  unsigned int            i, j, k, l, n, min_i, turn, type, *sn, u1, u2, u;
  int                     e, tmp, ppp, ij, kl, *c, *outside_c, *f5, *f3, *fML,
                          *idx, dangle_model;
  float                   mfe;
  struct ml_aux           *ml_helper;
  vrna_param_t            *P;
  vrna_md_t               *md;
  vrna_hc_t               *hc;
  vrna_sc_t               *sc;
  vrna_subopt_solution_t  *sol;

  sol = NULL;

  if (fc) {
    mfe_structure = (char *)vrna_alloc(sizeof(char) * (fc->length + 1));
    mfe   = vrna_mfe(fc, mfe_structure);

    printf("%s\n%s (%6.2f)\n", fc->sequence, mfe_structure, mfe);

    n     = fc->length;
    sn    = fc->strand_number;
    S     = fc->sequence_encoding2;
    S1    = fc->sequence_encoding;
    P     = fc->params;
    md    = &(P->model_details);
    dangle_model = md->dangles;
    turn  = md->min_loop_size;
    idx   = fc->jindx;
    f5    = fc->matrices->f5;
    c     = fc->matrices->c;
    fML   = fc->matrices->fML;
    hc    = fc->hc;
    sc    = fc->sc;
    outside_c = (int *)vrna_alloc(sizeof(int) * ((n * (n + 1)) / 2 + 2));
    ml_helper = init_ml_helper(fc);

    f3        = compute_f3(fc);

    /* initialize outside matrix */
    for (k = 0; k < (n * (n + 1)) / 2 + 1; k++)
      outside_c[k] = INF;

    /* backtrack (1,n) */
    if (hc->mx[n + n] & VRNA_CONSTRAINT_CONTEXT_EXT_LOOP) {
      kl            = idx[n] + 1;
      type          = vrna_get_ptype_md(S[1], S[n], md);
      e             = vrna_E_ext_stem(type, -1, -1, P);
      outside_c[kl] = e;

      if (c[kl] != INF) {
        e += c[kl];

        /* this is the only possibility for pair (1,n), so let's backtrack */
        s = backtrack_pair(fc);
        printf("%s (%6.2f) [0]\n", s, (float)(e / 100.));
        free(s);
      }
    }

    /* backtrack all structures with pairs (k, n) 1 < k < n */
    for (k = n - turn - 1; k > 1; k--) {
      if (hc->mx[n * n + k] & VRNA_CONSTRAINT_CONTEXT_EXT_LOOP) {
        kl = idx[n] + k;
        if ((sn[k - 1] == sn[k]) &&
            (f5[k - 1] != INF)) {
          type = vrna_get_ptype_md(S[k], S[n], md);
          e   = f5[k - 1] +
                vrna_E_ext_stem(type, S1[k - 1], -1, P);

          if (sc) {
            if (sc->f)
              e += sc->f(1, n, k - 1, k, VRNA_DECOMP_EXT_EXT_STEM, sc->data);
          }

          outside_c[kl] = e;

          if (c[kl] != INF) {
            e += c[kl];
            /* this is the only possibility for pair (k,n), so let's backtrack */
            s = backtrack_ext5(fc, k);
            printf("%s (%6.2f) [1]\n", s, (float)(e / 100.));
            free(s);
          }
        }
      }
    }

    /* backtrack all structures with pairs (1, k) 1 < k < n */
    for (k = n - turn - 1; k > 1; k--) {
      if (hc->mx[n + k]  & VRNA_CONSTRAINT_CONTEXT_EXT_LOOP) {
        kl = idx[k] + 1;
        if ((sn[k] == sn[k + 1]) &&
            (f3[k + 1] != INF)) {
          type = vrna_get_ptype_md(S[1], S[k], md);
          e   = f3[k + 1] +
                vrna_E_ext_stem(type, -1, S1[k + 1], P);

          if (sc) {
            if (sc->f)
              e += sc->f(1, n, k, k + 1, VRNA_DECOMP_EXT_STEM_EXT, sc->data);
          }

          outside_c[kl] = e;

          if (c[kl] != INF) {
            e += c[kl];
            /* this is the only possibility for pair (k,n), so let's backtrack */
            s = backtrack_ext3(fc, k, e, f3);
            printf("%s (%6.2f) [2]\n", s, (float)(e / 100.));
            free(s);
          }
        }
      }
    }

    /*
        now, for all the possibilities where (k,l) may be enclosed by
        at least one other pair
    */
    for (l = n - 1; l > turn + 1; l--) {
      prepare_ml_helper(fc, l, outside_c, ml_helper);

      for (k = 2; k < l - turn; k++) {
        int e_ext, e_int, e_mb;

        type  = vrna_get_ptype_md(S[k], S[l], md);
        kl    = idx[l] + k;
        e_ext = e_int = e_mb = INF;

        /* 1. (k,l) not enclosed by any other pair */
        if (hc->mx[l * n + k] & VRNA_CONSTRAINT_CONTEXT_EXT_LOOP) {
          if ((f5[k - 1] != INF) &&
              (f3[l + 1] != INF) &&
              (sn[k - 1] == sn[k]) &&
              (sn[l] == sn[l + 1])) {
            s5 = S1[k - 1];
            s3 = S1[l + 1];
            e_ext = f5[k - 1] +
                    f3[l + 1] +
                    vrna_E_ext_stem(type, s5, s3, P);

          }
        }

        /* 2. (k,l) enclosed by a single pair forming an internal loop */
        if (hc->mx[l * n + k] & VRNA_CONSTRAINT_CONTEXT_INT_LOOP_ENC) {
          for (j = l + 1; j <= MIN2(l + MAXLOOP + 1, n); j++) {
            u2 = j - l - 1;

            if (hc->up_int[l + 1] < u2)
              break;

            min_i = (k > MAXLOOP + 1) ? k - MAXLOOP - 1 : 1;

            for (i = k - 1; i >= min_i; i--) {
              ij = idx[j] + i;

              u1 = k - i - 1;

              if (hc->up_int[i + 1] < u1)
                break;

              if (hc->mx[j * n + i] & VRNA_CONSTRAINT_CONTEXT_INT_LOOP) {
                tmp = outside_c[ij] +
                      vrna_eval_int_loop(fc, i, j, k, l);

                e_int = MIN2(e_int, tmp);
              }
            }
          }
        }

        /* 3. (k,l) enclosed as part of a multibranch loop */
        if (hc->mx[l * n + k] & VRNA_CONSTRAINT_CONTEXT_MB_LOOP_ENC) {
          e = E_MLstem(type, S1[k - 1], S1[l + 1], P);

          for (i = k - 1; i > 0; i--) {
            /* left-most or somewhere in the middle */
            if (ml_helper->ml[i] != INF) {
              u = k - i - 1;
              /* left-most */
              if (hc->up_ml[i + 1] >= u) {
                ppp = u * P->MLbase +
                      e +
                      ml_helper->ml[i];

                e_mb = MIN2(e_mb, ppp);
              }

              /* somewhere in the middle */
              if (fML[idx[k - 1] + i + 1] != INF) {
                ppp = fML[idx[k - 1] + i + 1] +
                      e +
                      ml_helper->ml[i];

                e_mb = MIN2(e_mb, ppp);
              }
            }

            /* right-most */
            if ((ml_helper->up[i] != INF) &&
                (fML[idx[k - 1] + i + 1] != INF)) {
              ppp = fML[idx[k - 1] + i + 1] +
                    e +
                    ml_helper->up[i];

              e_mb = MIN2(e_mb, ppp);
            }
          }
        }

        e = MIN2(e_ext, e_int);
        e = MIN2(e, e_mb);

        if (e != INF) {
          outside_c[kl] = e;

          if (e == e_ext) {
            e += c[kl];
            /* backtrack from external pair (k,l) */
            s = backtrack_ext53(fc, k, l, e, f3);
            printf("%s (%6.2f) [3]\n", s, (float)(e / 100.));
            free(s);
          } else if (e == e_int) {
            e += c[kl];
            /* backtrack from internal loop pair (k,l) */
            s = backtrack_int(fc, k, l, e, outside_c);
            printf("%s (%6.2f) [4]\n", s, (float)(e / 100.));
            free(s);
          } else if (e == e_mb) {
            e += c[kl];
            /* backtrack from multibranch loop pair (k,l) */
            s = backtrack_mb(fc, k, l, e, outside_c);
            printf("%s (%6.2f) [5]\n", s, (float)(e / 100.));
            free(s);
          } else {
            /* do nothing, since there seems no solution for pair (k,l) */
          }
        }
      } /* ... end for (k = 2; ...)  */
    } /* ... end for (l = n - 1; ...) */

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

  for (i = 0; i <= n; i++)
    ml_helper->up[i] = INF;

  return ml_helper;
}


PRIVATE void
prepare_ml_helper(vrna_fold_compound_t  *fc,
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
  for (i = 0; i < l; i++)
    ml_helper->ml[i]  = INF;

  if (l > turn + 2) {
    for (j = l + turn + 1; j <= n; j++) {
      for (i = l - turn - 2; i > 0; i--) {
        if (hc->mx[j * n + i] & VRNA_CONSTRAINT_CONTEXT_MB_LOOP) {
          if ((outside_c[idx[j] + i] != INF) &&
              (fML[idx[j - 1] + l + 1] != INF)) {
            type = vrna_get_ptype_md(S[j], S[i], md);
            e = outside_c[idx[j] + i] +
                fML[idx[j - 1] + l + 1] +
                E_MLstem(type, S1[j - 1], S1[i + 1], P) +
                P->MLclosing;

            ml_helper->ml[i] = MIN2(ml_helper->ml[i], e);
          }
        }
      }
    }

    /*
        given that this function is called for successively
        decreasing values of l, up[i] can be updated by a
        single loop over all i
    */
    if (hc->up_ml[l + 1]) {
      for (i = l - turn - 2; i > 0; i--) {
        if (ml_helper->up[i] != INF) {
          e = ml_helper->up[i] +
              P->MLbase;

          ml_helper->up[i] = MIN2(ml_helper->up[i], e);
        }
      }
    }

    for (i = l - turn - 2; i > 0; i--) {
      if (hc->mx[(l + 1) * n + i] & VRNA_CONSTRAINT_CONTEXT_MB_LOOP) {
        type = vrna_get_ptype_md(S[l + 1], S[i], md);
        e = outside_c[idx[l + 1] + i] +
            E_MLstem(type, S1[l], S1[i + 1], P) +
            P->MLclosing;

        ml_helper->up[i] = MIN2(ml_helper->up[i], e);
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
  idx   = fc->jindx;
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

  for (j = min_j - 1; j > turn; j--) {
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
        type = vrna_get_ptype_md(S[j], S[n], md);
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
backtrack_pair(vrna_fold_compound_t *fc)
{
  char            *structure;
  int             s;
  unsigned int    n;
  sect            bt_stack[MAXSECTORS]; /* stack of partial structures for backtracking */
  vrna_bp_stack_t *bp;

  structure = NULL;
  n         = fc->length;
  s         = 0;

  /* add a guess of how many G's may be involved in a G quadruplex */
  bp = (vrna_bp_stack_t *)vrna_alloc(sizeof(vrna_bp_stack_t) * (4 * (1 + n / 2)));

  /* push interval enclosed by (i,j) on bt_stack */
  bt_stack[++s].i = 1;
  bt_stack[s].j   = n;
  bt_stack[s].ml  = 2;

  if (vrna_backtrack_from_intervals(fc, bp, bt_stack, s))
    structure = vrna_db_from_bp_stack(bp, n);

  free(bp);

  return structure;
}


PRIVATE char *
backtrack_ext5(vrna_fold_compound_t *fc,
               unsigned int         k)
{
  char              *structure;
  int               s;
  unsigned int      n;
  sect              bt_stack[MAXSECTORS]; /* stack of partial structures for backtracking */
  vrna_bp_stack_t   *bp;

  structure = NULL;
  n         = fc->length;
  s         = 0;

  /* add a guess of how many G's may be involved in a G quadruplex */
  bp = (vrna_bp_stack_t *)vrna_alloc(sizeof(vrna_bp_stack_t) * (4 * (1 + n / 2)));

  /* push 5' external loop interval (1,k - 1) on bt_stack */
  bt_stack[++s].i = 1;
  bt_stack[s].j   = k - 1;
  bt_stack[s].ml  = 0;

  /* push interval enclosed by (i,j) on bt_stack */
  bt_stack[++s].i = k;
  bt_stack[s].j   = n;
  bt_stack[s].ml  = 2;

  if (vrna_backtrack_from_intervals(fc, bp, bt_stack, s))
    structure = vrna_db_from_bp_stack(bp, n);

  free(bp);

  return structure;
}

PRIVATE char *
backtrack_ext3(vrna_fold_compound_t *fc,
              unsigned int          k,
              int                   e,
              int                   *f3)
{
  char            *structure;
  unsigned int    i, j, n, stack_count;
  int             s;
  sect            bt_stack[MAXSECTORS]; /* stack of partial structures for backtracking */
  vrna_bp_stack_t *bp;

  structure = NULL;
  n         = fc->length;
  s         = 0;

  /* add a guess of how many G's may be involved in a G quadruplex */
  bp = (vrna_bp_stack_t *)vrna_alloc(sizeof(vrna_bp_stack_t) * (4 * (1 + n / 2)));

  /* push interval enclosed by (1, k) on bt_stack */
  bt_stack[++s].i = 1;
  bt_stack[s].j   = k;
  bt_stack[s].ml  = 2;

  /* process remaining 3' part (external loop [k + 1, n] */
  stack_count = 0;
  i           = 0;
  j           = 0;
  k++;

  while (k <= n) {
    if (backtrack_f3(fc, &k, &i, &j, f3)) {
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

  /* finally, process all remaining inside segments we've encountered so far */
  if (vrna_backtrack_from_intervals(fc, bp, bt_stack, s))
    structure = vrna_db_from_bp_stack(bp, n);

  free(bp);

  return structure;
}


PRIVATE char *
backtrack_ext53(vrna_fold_compound_t *fc,
                unsigned int          k,
                unsigned int          l,
                int                   e,
                int                   *f3)
{
  char              *structure;
  int               s;
  unsigned int      i, j, n, stack_count;
  sect              bt_stack[MAXSECTORS]; /* stack of partial structures for backtracking */
  vrna_bp_stack_t   *bp;

  structure = NULL;
  n         = fc->length;
  s         = 0;

  /* add a guess of how many G's may be involved in a G quadruplex */
  bp = (vrna_bp_stack_t *)vrna_alloc(sizeof(vrna_bp_stack_t) * (4 * (1 + n / 2)));

  /* push 5' external loop interval (1,k - 1) on bt_stack */
  bt_stack[++s].i = 1;
  bt_stack[s].j   = k - 1;
  bt_stack[s].ml  = 0;

  /* push interval enclosed by (i,j) on bt_stack */
  bt_stack[++s].i = k;
  bt_stack[s].j   = l;
  bt_stack[s].ml  = 2;

  /* process remaining 3' part (external loop [k + 1, n] */
  stack_count = 0;
  i           = 0;
  j           = 0;
  k++;

  while (k <= n) {
    if (backtrack_f3(fc, &k, &i, &j, f3)) {
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


  if (vrna_backtrack_from_intervals(fc, bp, bt_stack, s))
    structure = vrna_db_from_bp_stack(bp, n);

  free(bp);

  return structure;

}


PRIVATE char *
backtrack_int(vrna_fold_compound_t  *fc,
              unsigned int          k,
              unsigned int          l,
              int                   e,
              int                   *outside_c)
{
  char  *structure;

  structure = NULL;

  return structure;
}


PRIVATE char *
backtrack_mb(vrna_fold_compound_t  *fc,
             unsigned int          k,
             unsigned int          l,
             int                   e,
             int                   *outside_c)
{
  char  *structure;

  structure = NULL;

  return structure;
}


PRIVATE int
backtrack_f3(vrna_fold_compound_t *fc,
             unsigned int         *k,
             unsigned int         *i,
             unsigned int         *j,
             int                  *f3)
{
  short         *S, *S1, s5, s3;
  unsigned int  ii, n, u, *sn, turn, type;
  int           en, fij, fi, dangle_model, *idx, *c;
  vrna_param_t  *P;
  vrna_md_t     *md;
  vrna_hc_t     *hc;
  vrna_sc_t     *sc;

  n             = fc->length;
  sn            = fc->strand_number;
  S             = fc->sequence_encoding2;
  S1            = fc->sequence_encoding;
  idx           = fc->jindx;
  c             = fc->matrices->c;
  P             = fc->params;
  md            = &(P->model_details);
  dangle_model  = md->dangles;
  turn          = md->min_loop_size;
  hc            = fc->hc;
  sc            = fc->sc;

  ii = *k;

  /* nibble-off unpaired bases on 5' side */
  do {
    fij = f3[ii];
    fi  = INF;

    /* do hard constraints check */
    if (sn[ii] == sn[ii + 1]) {
      fi = f3[ii + 1];

      if (sc) {
        if (sc->energy_up)
          fi += sc->energy_up[ii][1];

        if (sc->f)
          fi += sc->f(ii, n, ii + 1, n, VRNA_DECOMP_EXT_EXT, sc->data);
      }
    }
    if (++ii > n)
      break;
  } while (fij == fi);
  ii--;

  if (ii + turn - 1 > n) {
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
