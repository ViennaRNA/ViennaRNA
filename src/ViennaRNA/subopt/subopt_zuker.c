#include <stdlib.h>

#include "ViennaRNA/utils/basic.h"
#include "ViennaRNA/utils/log.h"
#include "ViennaRNA/sequences/alphabet.h"
#include "ViennaRNA/eval/exterior.h"
#include "ViennaRNA/eval/hairpin.h"
#include "ViennaRNA/eval/internal.h"
#include "ViennaRNA/eval/multibranch.h"
#include "ViennaRNA/constraints/hard.h"
#include "ViennaRNA/mfe/global.h"
#include "ViennaRNA/backtrack/global.h"
#include "ViennaRNA/subopt/zuker.h"

#define MAXSECTORS        500     /* dimension for a backtrack array */


typedef struct {
  int **mb;
  int **mb_up;
  int *f3;
  int *outside_c;
} zuker_aux_mx;


PRIVATE zuker_aux_mx *
get_zuker_aux_mx(vrna_fold_compound_t *fc);


PRIVATE void
prepare_ml_helper(vrna_fold_compound_t  *fc,
                  unsigned int          l,
                  zuker_aux_mx          *aux_mx);


PRIVATE void
free_zuker_aux_mx(zuker_aux_mx  *aux_mx,
                  unsigned int  n);


PRIVATE int *
compute_f3(vrna_fold_compound_t *fc);


PRIVATE int
backtrack(vrna_fold_compound_t  *fc,
          unsigned int          k,
          unsigned int          l,
          zuker_aux_mx          *aux_mx,
          vrna_bp_stack_t       *bp);


PRIVATE int
backtrack_mb(vrna_fold_compound_t *fc,
             unsigned int         i,
             unsigned int         *k,
             unsigned int         *l,
             zuker_aux_mx         *aux_mx);


PRIVATE int
backtrack_mb_up(vrna_fold_compound_t  *fc,
                unsigned int          i,
                unsigned int          *k,
                unsigned int          *l,
                zuker_aux_mx          *aux_mx);


PRIVATE int
backtrack_f3(vrna_fold_compound_t *fc,
             unsigned int         *k,
             unsigned int         *i,
             unsigned int         *j,
             int                  *f3);


typedef struct {
  int i;
  int j;
  int e;
  int idxj;
} zuker_pair;

PRIVATE int
comp_pair(const void  *A,
          const void  *B)
{
  zuker_pair  *x, *y;
  int         ex, ey;

  x   = (zuker_pair *)A;
  y   = (zuker_pair *)B;
  ex  = x->e;
  ey  = y->e;
  if (ex > ey)
    return 1;

  if (ex < ey)
    return -1;

  return x->idxj + x->i - y->idxj + y->i;
}


PUBLIC vrna_subopt_solution_t *
vrna_subopt_zuker(vrna_fold_compound_t *fc)
{
  unsigned char           **todo;
  char                    *s;
  short                   *S, *S1;
  unsigned int            i, j, k, l, n, min_i, type, *sn, u1, u2, u,
                          num_pairs, num_struct;
  int                     e, tmp, ppp, ij, kl, *c, *outside_c, *f5, *f3, *fML,
                          *idx, dangle_model;
  vrna_param_t            *P;
  vrna_md_t               *md;
  vrna_hc_t               *hc;
  vrna_sc_t               *sc;
  vrna_subopt_solution_t  *sol;
  vrna_bp_stack_t         *bp;
  zuker_pair              *pairlist;
  zuker_aux_mx            *aux_mx;

  sol = NULL;

  if (fc) {
    (void)vrna_mfe(fc, NULL);

    n             = fc->length;
    sn            = fc->strand_number;
    S             = fc->sequence_encoding2;
    S1            = fc->sequence_encoding;
    P             = fc->params;
    md            = &(P->model_details);
    dangle_model  = md->dangles;
    idx           = fc->jindx;
    f5            = fc->matrices->f5;
    c             = fc->matrices->c;
    fML           = fc->matrices->fML;
    hc            = fc->hc;
    sc            = fc->sc;

    aux_mx = get_zuker_aux_mx(fc);

    outside_c = aux_mx->outside_c;
    f3        = aux_mx->f3;

    /* backtrack (1,n) */
    if (hc->mx[n + n] & VRNA_CONSTRAINT_CONTEXT_EXT_LOOP) {
      kl            = idx[n] + 1;
      type          = vrna_get_ptype_md(S[1], S[n], md);
      e             = vrna_E_exterior_stem(type, -1, -1, P);
      outside_c[kl] = e;
    }

    /* backtrack all structures with pairs (k, n) 1 < k < n */
    for (k = n - 1; k > 1; k--) {
      if (hc->mx[n * n + k] & VRNA_CONSTRAINT_CONTEXT_EXT_LOOP) {
        kl = idx[n] + k;
        if ((sn[k - 1] == sn[k]) &&
            (f5[k - 1] != INF)) {
          type  = vrna_get_ptype_md(S[k], S[n], md);
          e     = f5[k - 1];

          switch (dangle_model) {
            case 2:
              e += vrna_E_exterior_stem(type, S1[k - 1], -1, P);
              break;

            default:
              e += vrna_E_exterior_stem(type, -1, -1, P);
              break;
          }

          if (sc)
            if (sc->f)
              e += sc->f(1, n, k - 1, k, VRNA_DECOMP_EXT_EXT_STEM, sc->data);

          outside_c[kl] = e;
        }
      }
    }

    /* backtrack all structures with pairs (1, k) 1 < k < n */
    for (k = n - 1; k > 1; k--) {
      if (hc->mx[n + k] & VRNA_CONSTRAINT_CONTEXT_EXT_LOOP) {
        kl = idx[k] + 1;
        if ((sn[k] == sn[k + 1]) &&
            (f3[k + 1] != INF)) {
          type  = vrna_get_ptype_md(S[1], S[k], md);
          e     = f3[k + 1];

          switch (dangle_model) {
            case 2:
              e += vrna_E_exterior_stem(type, -1, S1[k + 1], P);
              break;

            default:
              e += vrna_E_exterior_stem(type, -1, -1, P);
              break;
          }

          if (sc)
            if (sc->f)
              e += sc->f(1, n, k, k + 1, VRNA_DECOMP_EXT_STEM_EXT, sc->data);

          outside_c[kl] = e;
        }
      }
    }

    /*
     *  now, for all the possibilities where (k,l) may be enclosed by
     *  at least one other pair
     */
    for (l = n - 1; l > 1; l--) {
      prepare_ml_helper(fc, l, aux_mx);

      for (k = 2; k < l; k++) {
        int e_ext, e_int, e_mb;

        type  = vrna_get_ptype_md(S[k], S[l], md);
        kl    = idx[l] + k;
        e_ext = e_int = e_mb = INF;

        /* 1. (k,l) is external pair */
        if (hc->mx[l * n + k] & VRNA_CONSTRAINT_CONTEXT_EXT_LOOP) {
          /* 1.a (k,l) not enclosed by any other pair */
          if ((f5[k - 1] != INF) &&
              (f3[l + 1] != INF) &&
              (sn[k - 1] == sn[k]) &&
              (sn[l] == sn[l + 1])) {
            e_ext = f5[k - 1] +
                    f3[l + 1];

            switch (dangle_model) {
              case 2:
                e_ext += vrna_E_exterior_stem(type, S1[k - 1], S1[l + 1], P);
                break;
              default:
                e_ext += vrna_E_exterior_stem(type, -1, -1, P);
            }

            if (sc)
              if (sc->f)
                e_ext += sc->f(1, n, k, l, VRNA_DECOMP_EXT_STEM_OUTSIDE, sc->data);
          }

          /* 1.b (k,l) is enclosed by another pair in a loop with strand nick, a.k.a. multi strand case */
          if (fc->strands > 1) {
          }
        }

        /* 2. (k,l) enclosed by a single pair forming an internal loop */
        if (hc->mx[l * n + k] & VRNA_CONSTRAINT_CONTEXT_INT_LOOP_ENC) {
          for (j = l + 1; j <= MIN2(l + MAXLOOP + 1, n); j++) {
            u2 = j - l - 1;

            if (hc->up_int[l + 1] < u2)
              break;

            min_i = (k > MAXLOOP + 1) ? k - MAXLOOP - 1 : 1;
            if (u2 + k - min_i - 1 > MAXLOOP)
              min_i = k - 1 - (MAXLOOP - u2);

            for (i = k - 1; i >= min_i; i--) {
              ij = idx[j] + i;

              u1 = k - i - 1;

              if (hc->up_int[i + 1] < u1)
                break;

              if (hc->mx[j * n + i] & VRNA_CONSTRAINT_CONTEXT_INT_LOOP) {
                tmp = outside_c[ij] +
                      vrna_eval_internal(fc, i, j, k, l, VRNA_EVAL_LOOP_NO_HC);

                e_int = MIN2(e_int, tmp);
              }
            }
          }
        }

        /* 3. (k,l) enclosed as part of a multibranch loop */
        if ((hc->mx[l * n + k] & VRNA_CONSTRAINT_CONTEXT_MB_LOOP_ENC) &&
            (sn[l] == sn[l + 1]) &&
            (sn[k - 1] == sn[k])) {
          int *aux_mb     = aux_mx->mb[l];
          int *aux_mb_up  = aux_mx->mb_up[l];

          switch (dangle_model) {
            case 2:
              e = vrna_E_multibranch_stem(type, S1[k - 1], S1[l + 1], P);
              break;
            default:
              e = vrna_E_multibranch_stem(type, -1, -1, P);
              break;
          }

          for (i = k - 1; i > 0; i--) {
            /* left-most or somewhere in the middle */
            if (aux_mb[i] != INF) {
              u = k - i - 1;
              /* left-most */
              if ((hc->up_ml[i + 1] >= u) &&
                  (sn[i] == sn[k])) {
                ppp = u * P->MLbase +
                      e +
                      aux_mb[i];

                if (sc) {
                  if (sc->energy_up)
                    ppp += sc->energy_up[i + 1][u];

                  if (sc->f)
                    ppp += sc->f(i + 1, l, k - 1, k, VRNA_DECOMP_ML_ML_STEM, sc->data) +
                           sc->f(i + 1, k - 1, i + 1, k - 1, VRNA_DECOMP_ML_UP, sc->data);
                }

                e_mb = MIN2(e_mb, ppp);
              }

              /* somewhere in the middle */
              if ((i + 2 < k) &&
                  (fML[idx[k - 1] + i + 1] != INF) &&
                  (sn[i] == sn[i + 1])) {
                ppp = fML[idx[k - 1] + i + 1] +
                      e +
                      aux_mb[i];

                if (sc)
                  if (sc->f)
                    ppp += sc->f(i + 1, l, k - 1, k, VRNA_DECOMP_ML_ML_STEM, sc->data);

                e_mb = MIN2(e_mb, ppp);
              }
            }

            /* right-most */
            if ((i + 2 < k) &&
                (aux_mb_up[i] != INF) &&
                (fML[idx[k - 1] + i + 1] != INF) &&
                (sn[i] == sn[i + 1])) {
              ppp = fML[idx[k - 1] + i + 1] +
                    e +
                    aux_mb_up[i];

              if (sc)
                if (sc->f)
                  ppp += sc->f(i + 1, l, k - 1, k, VRNA_DECOMP_ML_ML_STEM, sc->data);

              e_mb = MIN2(e_mb, ppp);
            }
          }
        }

        e             = MIN2(e_ext, e_int);
        e             = MIN2(e, e_mb);
        outside_c[kl] = e;
      } /* ... end for (k = 2; ...)  */
    }   /* ... end for (l = n - 1; ...) */

    /* now, for the actual backtracking */

    /* make todo-list of base pairs that require processing */
    todo = (unsigned char **)vrna_alloc(sizeof(unsigned char *) * n);
    for (k = 1; k < n; k++)
      todo[k] = (unsigned char *)vrna_alloc(sizeof(unsigned char) * (n + 1));

    pairlist    = (zuker_pair *)vrna_alloc(sizeof(zuker_pair) * ((n * (n + 1)) / 2 + 2));
    num_pairs   = 0;  /* number of pairs to process */
    num_struct  = 0;  /* number of Zuker suboptimal structures */

    for (l = n; l > 1; l--) {
      int idxj = idx[l];
      for (k = 1; k < l; k++) {
        if (hc->mx[l * n + k] & VRNA_CONSTRAINT_CONTEXT_ALL_LOOPS) {
          pairlist[num_pairs].i     = k;
          pairlist[num_pairs].j     = l;
          pairlist[num_pairs].e     = c[idxj + k] + outside_c[idxj + k];
          pairlist[num_pairs].idxj  = idxj;
          num_pairs++;
          todo[k][l] = 1;
        }
      }
    }

    sol = (vrna_subopt_solution_t *)vrna_alloc(sizeof(vrna_subopt_solution_t) * (num_pairs + 1));

    /* resize list to actual requirements */
    pairlist = vrna_realloc(pairlist, sizeof(zuker_pair) * (num_pairs + 1));

    pairlist[num_pairs].i = 0; /* end of list marker */

    /* sort to enable backtrack pruning */
    qsort(pairlist, num_pairs, sizeof(zuker_pair), comp_pair);

    /* go through pair list and start backtracking */
    bp = (vrna_bp_stack_t *)vrna_alloc(sizeof(vrna_bp_stack_t) * (4 * (1 + n / 2)));
    for (i = 0; pairlist[i].i; i++) {
      k = pairlist[i].i;
      l = pairlist[i].j;
      e = pairlist[i].e;

      if (todo[k][l]) {
        /* add a guess of how many G's may be involved in a G quadruplex */
        bp[0].i = 0;

        if (backtrack(fc, k, l, aux_mx, bp)) {
          s = vrna_db_from_bp_stack(bp, n);

          for (j = 1; j <= bp[0].i; j++)
            todo[bp[j].i][bp[j].j] = 0;

          sol[num_struct].energy      = (float)e / 100.;
          sol[num_struct++].structure = s;
        } else {
          vrna_log_warning("Backtracking failed for pair (%d,%d) en=%d", k, l, e);
        }
      }
    }

    /* resize solution list to actual needs */
    sol =
      (vrna_subopt_solution_t *)vrna_realloc(sol,
                                             sizeof(vrna_subopt_solution_t) * (num_struct + 1));
    sol[num_struct].structure = NULL; /* end of list marker */

    /* clean up memory */
    free(pairlist);
    free(bp);

    for (k = 1; k < n; k++)
      free(todo[k]);

    free(todo);

    free_zuker_aux_mx(aux_mx, n);
  }

  return sol;
}


PRIVATE zuker_aux_mx *
get_zuker_aux_mx(vrna_fold_compound_t *fc)
{
  unsigned int  i, n;
  zuker_aux_mx  *mx;

  n   = fc->length;
  mx  = (zuker_aux_mx *)vrna_alloc(sizeof(zuker_aux_mx));

  mx->outside_c = (int *)vrna_alloc(sizeof(int) * ((n * (n + 1)) / 2 + 2));
  mx->mb        = (int **)vrna_alloc(sizeof(int *) * (n + 1));
  mx->mb_up     = (int **)vrna_alloc(sizeof(int *) * (n + 1));

  for (i = 1; i <= n; i++) {
    mx->mb[i]     = (int *)vrna_alloc(sizeof(int) * (n + 1));
    mx->mb_up[i]  = (int *)vrna_alloc(sizeof(int) * (n + 1));
  }

  /* initialize outside matrix */
  for (i = 0; i < (n * (n + 1)) / 2 + 1; i++)
    mx->outside_c[i] = INF;

  for (i = 0; i <= n; i++)
    mx->mb_up[n][i] = INF;

  mx->f3 = compute_f3(fc);

  return mx;
}


PRIVATE void
prepare_ml_helper(vrna_fold_compound_t  *fc,
                  unsigned int          l,
                  zuker_aux_mx          *aux_mx)
{
  short         *S, *S1;
  unsigned int  i, j, n, type, *sn;
  int           e, *idx, *fML, *outside_c, *aux_mb,
                *aux_mb_up, *aux_mb_up1, dangle_model;
  vrna_param_t  *P;
  vrna_md_t     *md;
  vrna_hc_t     *hc;
  vrna_sc_t     *sc;

  n             = fc->length;
  S             = fc->sequence_encoding2;
  S1            = fc->sequence_encoding;
  sn            = fc->strand_number;
  idx           = fc->jindx;
  P             = fc->params;
  md            = &(P->model_details);
  dangle_model  = md->dangles;
  hc            = fc->hc;
  sc            = fc->sc;
  fML           = fc->matrices->fML;
  outside_c     = aux_mx->outside_c;
  aux_mb        = aux_mx->mb[l];
  aux_mb_up     = aux_mx->mb_up[l];
  aux_mb_up1    = aux_mx->mb_up[l + 1];

  /* initialize with INF */
  for (i = 0; i < l; i++) {
    aux_mb[i]     = INF;
    aux_mb_up[i]  = INF;
  }

  if ((l > 2) &&
      (sn[l] == sn[l + 1])) {
    for (j = l + 3; j <= n; j++) {
      if (sn[j] == sn[j - 1]) {
        for (i = l - 2; i > 0; i--) {
          if ((hc->mx[j * n + i] & VRNA_CONSTRAINT_CONTEXT_MB_LOOP) &&
              (sn[i] == sn[i + 1]) &&
              (outside_c[idx[j] + i] != INF) &&
              (fML[idx[j - 1] + l + 1] != INF)) {
            type  = vrna_get_ptype_md(S[j], S[i], md);
            e     = outside_c[idx[j] + i] +
                    fML[idx[j - 1] + l + 1] +
                    P->MLclosing;

            switch (dangle_model) {
              case 2:
                e += vrna_E_multibranch_stem(type, S1[j - 1], S1[i + 1], P);
                break;
              default:
                e += vrna_E_multibranch_stem(type, -1, -1, P);
                break;
            }

            if (sc) {
              if (sc->f)
                e += sc->f(i, j, i + 1, j - 1, VRNA_DECOMP_PAIR_ML, sc->data) +
                     sc->f(i + 1, j - 1, l, l + 1, VRNA_DECOMP_ML_ML_ML, sc->data);
            }

            aux_mb[i] = MIN2(aux_mb[i], e);
          }
        }
      }
    }

    /*
     *  given that this function is called for successively
     *  decreasing values of l, up[i] can be updated by a
     *  single loop over all i
     */
    if ((hc->up_ml[l + 1]) &&
        (sn[l] == sn[l + 1])) {
      for (i = l - 2; i > 0; i--) {
        if (aux_mb_up1[i] != INF) {
          e = aux_mb_up1[i] +
              P->MLbase;

          if (sc) {
            if (sc->energy_up)
              e += sc->energy_up[l + 1][1];

            if (sc->f)
              e += sc->f(i + 1, l + 1, i + 1, l, VRNA_DECOMP_ML_ML, sc->data);
          }

          aux_mb_up[i] = e;
        }
      }
    }

    for (i = l - 2; i > 0; i--) {
      if ((hc->mx[(l + 1) * n + i] & VRNA_CONSTRAINT_CONTEXT_MB_LOOP) &&
          (sn[i] == sn[i + 1])) {
        type  = vrna_get_ptype_md(S[l + 1], S[i], md);
        e     = outside_c[idx[l + 1] + i] +
                P->MLclosing;

        switch (dangle_model) {
          case 2:
            e += vrna_E_multibranch_stem(type, S1[l], S1[i + 1], P);
            break;
          default:
            e += vrna_E_multibranch_stem(type, -1, -1, P);
            break;
        }

        if (sc)
          if (sc->f)
            e += sc->f(i, l + 1, i + 1, l, VRNA_DECOMP_PAIR_ML, sc->data);

        aux_mb_up[i] = MIN2(aux_mb_up[i], e);
      }
    }
  }
}


PRIVATE void
free_zuker_aux_mx(zuker_aux_mx  *aux_mx,
                  unsigned int  n)
{
  unsigned int i;

  for (i = 1; i <= n; i++) {
    free(aux_mx->mb[i]);
    free(aux_mx->mb_up[i]);
  }

  free(aux_mx->mb);
  free(aux_mx->mb_up);

  free(aux_mx->f3);
  free(aux_mx->outside_c);

  free(aux_mx);
}


PRIVATE int *
compute_f3(vrna_fold_compound_t *fc)
{
  short         *S, *S1;
  unsigned int  j, k, n, jk, *sn, type;
  int           e, *c, *f3, *idx, dangle_model;
  vrna_param_t  *P;
  vrna_md_t     *md;
  vrna_hc_t     *hc;
  vrna_sc_t     *sc;

  n             = fc->length;
  S             = fc->sequence_encoding2;
  S1            = fc->sequence_encoding;
  sn            = fc->strand_number;
  P             = fc->params;
  md            = &(P->model_details);
  dangle_model  = md->dangles;
  idx           = fc->jindx;
  c             = fc->matrices->c;
  hc            = fc->hc;
  sc            = fc->sc;

  /* init */
  f3 = (int *)vrna_alloc(sizeof(int) * (n + 2));

  f3[n + 1] = 0;
  f3[n]     = INF;

  if ((hc->up_ext[n]) &&
      (sn[n - 1] == sn[n])) {
    f3[n] = 0;
    if (sc) {
      if (sc->energy_up)
        f3[n] += sc->energy_up[n][1];

      if (sc->f)
        f3[n] += sc->f(n, n, n, n, VRNA_DECOMP_EXT_UP, sc->data);
    }
  }

  for (j = n - 1; j >= 1; j--) {
    /* 1st case, j is unpaired */
    if ((hc->up_ext[j]) &&
        (sn[j] == sn[j + 1])) {
      e = f3[j + 1];

      if (sc) {
        if (sc->energy_up)
          e += sc->energy_up[j][1];

        if (sc->f)
          e += sc->f(j, n, j + 1, n, VRNA_DECOMP_EXT_EXT, sc->data);
      }

      f3[j] = MIN2(f3[j], e);
    }

    /* 2nd case, j forms pair (j,k) with j < k < n */
    for (k = j + 1; k < n; k++) {
      if (hc->mx[j * n + k] & VRNA_CONSTRAINT_CONTEXT_EXT_LOOP) {
        jk = idx[k] + j;
        if ((c[jk] != INF) &&
            (f3[k + 1] != INF) &&
            (sn[k] == sn[k + 1])) {
          type  = vrna_get_ptype_md(S[j], S[k], md);
          e     = c[jk] +
                  f3[k + 1];

          switch (dangle_model) {
            case 2:
              e += vrna_E_exterior_stem(type, S1[j - 1], S1[k + 1], P);
              break;
            default:
              e += vrna_E_exterior_stem(type, -1, -1, P);
              break;
          }

          if (sc)
            if (sc->f)
              e += sc->f(j, n, k, k + 1, VRNA_DECOMP_EXT_STEM_EXT, sc->data);

          f3[j] = MIN2(f3[j], e);
        }
      }
    }

    /* 3rd case, j forms pair with (j, n) */
    if (hc->mx[j * n + n] & VRNA_CONSTRAINT_CONTEXT_EXT_LOOP) {
      jk = idx[n] + j;
      if (c[jk] != INF) {
        type  = vrna_get_ptype_md(S[j], S[n], md);
        e     = c[jk];
        switch (dangle_model) {
          case 2:
            e += vrna_E_exterior_stem(type, S1[j - 1], -1, P);
            break;
          default:
            e += vrna_E_exterior_stem(type, -1, -1, P);
            break;
        }

        if (sc)
          if (sc->f)
            e += sc->f(j, n, n, k, VRNA_DECOMP_EXT_STEM, sc->data);

        f3[j] = MIN2(f3[j], e);
      }
    }
  }

  return f3;
}


PRIVATE int
backtrack(vrna_fold_compound_t  *fc,
          unsigned int          k,
          unsigned int          l,
          zuker_aux_mx          *aux_mx,
          vrna_bp_stack_t       *bp)
{
  short         *S, *S1, s5, s3;
  unsigned int  n, i, j, b, s, type, u1, u2, max_j, min_i, prev_l, prev_k,
                *sn;
  int           e, tmp, en, *f5, *idx, kl, ij, *outside_c,
                *f3, *fML, **aux_mb, **aux_mb_up, dangle_model, *mb, *mb_up;
  sect          bt_stack[MAXSECTORS];
  vrna_param_t  *P;
  vrna_md_t     *md;
  vrna_hc_t     *hc;
  vrna_sc_t     *sc;

  n             = fc->length;
  S             = fc->sequence_encoding2;
  S1            = fc->sequence_encoding;
  sn            = fc->strand_number;
  idx           = fc->jindx;
  P             = fc->params;
  md            = &(P->model_details);
  dangle_model  = md->dangles;
  f5            = fc->matrices->f5;
  hc            = fc->hc;
  sc            = fc->sc;
  s             = 0;
  b             = bp[0].i; /* number of already backtraced outside pairs */
  outside_c     = aux_mx->outside_c;
  f3            = aux_mx->f3;
  fML           = fc->matrices->fML;
  aux_mb        = aux_mx->mb;
  aux_mb_up     = aux_mx->mb_up;

  /* push interval enclosed by (i,j) on bt_stack */
  bt_stack[s].i     = k;
  bt_stack[s].j     = l;
  bt_stack[s++].ml  = VRNA_MX_FLAG_C;

backtrack_outside:

  kl    = idx[l] + k;
  e     = outside_c[kl];
  type  = vrna_get_ptype_md(S[k], S[l], md);

  /* 1st case, (k,l) enclosed by a single pair (i,j) forming an internal loop */
  if ((k > 1) &&
      (l < n) &&
      (hc->mx[k * n + l] & VRNA_CONSTRAINT_CONTEXT_INT_LOOP_ENC)) {
    min_i = (k > MAXLOOP + 1) ? k - MAXLOOP - 1 : 1;
    u1    = 0;

    for (i = k - 1; i >= min_i; i--, u1++) {
      if (hc->up_int[i + 1] < u1)
        break;

      if (sn[i] != sn[k])
        break;

      max_j = l + MAXLOOP + 1 - u1;
      if (max_j > n)
        max_j = n;

      u2 = 0;
      for (j = l + 1; j <= max_j; j++, u2++) {
        if (hc->up_int[l + 1] < u2)
          break;

        if (sn[l] != sn[j])
          break;

        if (hc->mx[j * n + i] & VRNA_CONSTRAINT_CONTEXT_INT_LOOP) {
          ij  = idx[j] + i;
          tmp = vrna_eval_internal(fc, i, j, k, l, VRNA_EVAL_LOOP_NO_HC);

          if (e == tmp + outside_c[ij]) {
            bp[++b].i = i;
            bp[b].j   = j;
            k         = i;
            l         = j;
            goto backtrack_outside;
          }
        }
      }
    }
  }

  /* 2nd case, (k,l) enclosed by a pair (i,j) forming a multibranch loop */
  if ((k > 1) &&
      (l < n) &&
      (hc->mx[k * n + l] & VRNA_CONSTRAINT_CONTEXT_MB_LOOP_ENC) &&
      (sn[k - 1] == sn[k]) &&
      (sn[l] == sn[l + 1])) {
    mb      = aux_mb[l];
    mb_up   = aux_mb_up[l];
    prev_k  = k;
    prev_l  = l;
    min_i   = 1;

    if (k > 1)
      min_i = k - 1;

    u1 = 0;

    switch (dangle_model) {
      case 2:
        tmp = vrna_E_multibranch_stem(type, S1[k - 1], S1[l + 1], P);
        break;
      default:
        tmp = vrna_E_multibranch_stem(type, -1, -1, P);
        break;
    }

    if (sc) {
      if (sc->f)
        tmp += sc->f(k, l, k, l, VRNA_DECOMP_ML_STEM, sc->data);
    }

    for (i = k - 1; i >= min_i; i--, u1++) {
      if (mb[i] == INF)
        continue;

      if (hc->up_ml[i + 1] < u1)
        break;

      if (sn[i + 1] != sn[k])
        break;

      en = tmp +
           mb[i] +
           u1 * P->MLbase;

      if (sc) {
        if (sc->energy_up)
          en += sc->energy_up[i + 1][u1];

        if (sc->f)
          en += sc->f(i + 1, l, k, l, VRNA_DECOMP_ML_ML, sc->data);
      }

      if (e == en) {
        if (backtrack_mb(fc, i, &k, &l, aux_mx)) {
          bt_stack[s].i     = prev_l + 1;
          bt_stack[s].j     = l - 1;
          bt_stack[s++].ml  = VRNA_MX_FLAG_M;

          bp[++b].i = k;
          bp[b].j   = l;

          goto backtrack_outside;
        } else {
          vrna_log_warning("backtracking failed for mb[%d] 1", i);
        }

        break;
      }
    }

    for (; i >= min_i; i--)
      u1++;

    for (; i > 0; i--, u1++) {
      if (mb[i] != INF) {
        if ((hc->up_ml[i + 1] >= u1) &&
            (sn[i + 1] == sn[k])) {
          en = tmp +
               mb[i] +
               u1 * P->MLbase;

          if (sc) {
            if (sc->energy_up)
              en += sc->energy_up[i + 1][u1];

            if (sc->f)
              en += sc->f(i + 1, l, k, l, VRNA_DECOMP_ML_ML, sc->data);
          }

          /* (k,l) is left-most pair in a multibranch loop */
          if (e == en) {
            if (backtrack_mb(fc, i, &k, &l, aux_mx)) {
              bt_stack[s].i     = prev_l + 1;
              bt_stack[s].j     = l - 1;
              bt_stack[s++].ml  = VRNA_MX_FLAG_M;

              bp[++b].i = k;
              bp[b].j   = l;

              goto backtrack_outside;
            } else {
              vrna_log_warning("backtracking failed for mb[%d] 2", i);
            }

            break;
          }
        }

        /* (k,l) is somewhere in the middle of a multibranch loop */
        if (sn[k - 1] == sn[k]) {
          en = tmp +
               mb[i] +
               fML[idx[k - 1] + i + 1];

          if (sc)
            if (sc->f)
              en += sc->f(i + 1, l, k - 1, k, VRNA_DECOMP_ML_ML_ML, sc->data);

          if (e == en) {
            if (backtrack_mb(fc, i, &k, &l, aux_mx)) {
              bt_stack[s].i     = prev_l + 1;
              bt_stack[s].j     = l - 1;
              bt_stack[s++].ml  = VRNA_MX_FLAG_M;

              bt_stack[s].i     = i + 1;
              bt_stack[s].j     = prev_k - 1;
              bt_stack[s++].ml  = VRNA_MX_FLAG_M;

              bp[++b].i = k;
              bp[b].j   = l;

              goto backtrack_outside;
            } else {
              vrna_log_warning("backtracking failed for pair (%d,%d) in mb[%d] 3", k, l, i);
            }

            break;
          }
        }
      }

      /* (k,l) is right-most pair of a multibranch loop */
      if ((mb_up[i] != INF) &&
          (sn[k - 1] == sn[k])) {
        en = tmp +
             mb_up[i] +
             fML[idx[k - 1] + i + 1];

        if (sc)
          if (sc->f)
            en += sc->f(i + 1, l, k - 1, k, VRNA_DECOMP_ML_ML_ML, sc->data);

        if (e == en) {
          if (backtrack_mb_up(fc, i, &k, &l, aux_mx)) {
            bt_stack[s].i     = i + 1;
            bt_stack[s].j     = prev_k - 1;
            bt_stack[s++].ml  = VRNA_MX_FLAG_M;

            bp[++b].i = k;
            bp[b].j   = l;

            goto backtrack_outside;
          } else {
            vrna_log_warning("backtracking failed for mb_up[%d]\n", i);
          }

          break;
        }
      }
    }
  }

  /* 3rd and last chance, (k,l) is not enclosed by any other pair */
  if ((hc->mx[k * n + l] & VRNA_CONSTRAINT_CONTEXT_EXT_LOOP) &&
      ((k == 1) || (sn[k - 1] == sn[k])) &&
      ((l == n) || (sn[l] == sn[l + 1]))) {
    switch (dangle_model) {
      case 2:
        s5  = (k > 1) ? S1[k - 1] : -1;
        s3  = (l < n) ? S1[l + 1] : -1;
        en  = vrna_E_exterior_stem(type, s5, s3, P);
        break;
      default:
        en = vrna_E_exterior_stem(type, -1, -1, P);
        break;
    }

    if (k > 1)
      en += f5[k - 1];

    if (l < n)
      en += f3[l + 1];

    if (sc)
      if (sc->f)
        en += sc->f(1, n, k, l, VRNA_DECOMP_EXT_STEM_OUTSIDE, sc->data);

    if (e == en) {
      if (k > 1) {
        /* push 5' external loop interval (1,k - 1) on bt_stack */
        bt_stack[s].i     = 1;
        bt_stack[s].j     = k - 1;
        bt_stack[s++].ml  = VRNA_MX_FLAG_F5;
      }

      if (l < n) {
        /* process remaining 3' part (external loop [k + 1, n] */
        i = 0;
        j = 0;
        k = l + 1;

        while (k <= n) {
          if (backtrack_f3(fc, &k, &i, &j, f3)) {
            if (i > 0) {
              /* store base pair for future backtracking */
              bt_stack[s].i     = i;
              bt_stack[s].j     = j;
              bt_stack[s++].ml  = VRNA_MX_FLAG_C;
            }
          } else {
            vrna_log_warning("Backtracking failed in f3[%d] = %d", k, f3[k]);
          }

          if (k == 0)
            break;
        }
      }

      goto backtrack_inside;
    }
  }

  goto backtrack_fail;

backtrack_inside:
  bp[0].i = b;

  vrna_log_error("Added %d bps", b);
  return vrna_backtrack_from_intervals(fc, bp, bt_stack, s);

backtrack_fail:

  return 0;
}


PRIVATE int
backtrack_mb(vrna_fold_compound_t *fc,
             unsigned int         i,
             unsigned int         *k,
             unsigned int         *l,
             zuker_aux_mx         *aux_mx)
{
  short         *S, *S1;
  unsigned int  j, n, type, *sn;
  int           e, en, *outside_c, *idx, *fML, dangle_model;
  vrna_param_t  *P;
  vrna_md_t     *md;
  vrna_hc_t     *hc;
  vrna_sc_t     *sc;

  n             = fc->length;
  S             = fc->sequence_encoding2;
  S1            = fc->sequence_encoding;
  sn            = fc->strand_number;
  idx           = fc->jindx;
  fML           = fc->matrices->fML;
  outside_c     = aux_mx->outside_c;
  P             = fc->params;
  md            = &(P->model_details);
  dangle_model  = md->dangles;
  hc            = fc->hc;
  sc            = fc->sc;
  e             = aux_mx->mb[*l][i];

  for (j = *l + 3; j <= n; j++) {
    if ((hc->mx[i * n + j] & VRNA_CONSTRAINT_CONTEXT_MB_LOOP) &&
        (outside_c[idx[j] + i] != INF) &&
        (fML[idx[j - 1] + *l + 1] != INF) &&
        (sn[j - 1] == sn[j])) {
      type  = vrna_get_ptype_md(S[j], S[i], md);
      en    = outside_c[idx[j] + i] +
              fML[idx[j - 1] + *l + 1] +
              P->MLclosing;

      switch (dangle_model) {
        case 2:
          en += vrna_E_multibranch_stem(type, S1[j - 1], S1[i + 1], P);
          break;
        default:
          en += vrna_E_multibranch_stem(type, -1, -1, P);
          break;
      }

      if (sc) {
        if (sc->f)
          en += sc->f(i, j, i + 1, j - 1, VRNA_DECOMP_PAIR_ML, sc->data) +
                sc->f(i + 1, j - 1, *l, *l + 1, VRNA_DECOMP_ML_ML_ML, sc->data);
      }

      if (e == en) {
        *k  = i;
        *l  = j;
        return 1;
      }
    }
  }

  return 0;
}


PRIVATE int
backtrack_mb_up(vrna_fold_compound_t  *fc,
                unsigned int          i,
                unsigned int          *k,
                unsigned int          *l,
                zuker_aux_mx          *aux_mx)
{
  short         *S, *S1;
  unsigned int  j, u, n, type, *sn;
  int           e, en, *outside_c, *idx, dangle_model;
  vrna_param_t  *P;
  vrna_md_t     *md;
  vrna_hc_t     *hc;
  vrna_sc_t     *sc;

  n             = fc->length;
  S             = fc->sequence_encoding2;
  S1            = fc->sequence_encoding;
  sn            = fc->strand_number;
  idx           = fc->jindx;
  P             = fc->params;
  md            = &(P->model_details);
  dangle_model  = md->dangles;
  outside_c     = aux_mx->outside_c;
  hc            = fc->hc;
  sc            = fc->sc;
  e             = aux_mx->mb_up[*l][i];
  u             = 0;

  /* find pairing partner j */
  for (j = *l + 1; j <= n; j++) {
    if ((hc->mx[i * n + j] & VRNA_CONSTRAINT_CONTEXT_MB_LOOP) &&
        (outside_c[idx[j] + i] != INF) &&
        (sn[*l] == sn[j])) {
      type  = vrna_get_ptype_md(S[j], S[i], md);
      en    = outside_c[idx[j] + i] +
              u * P->MLbase +
              P->MLclosing;

      switch (dangle_model) {
        case 2:
          en += vrna_E_multibranch_stem(type, S1[j - 1], S1[i + 1], P);
          break;
        default:
          en += vrna_E_multibranch_stem(type, -1, -1, P);
          break;
      }

      if (sc) {
        if (sc->energy_up)
          en += sc->energy_up[*l + 1][u];

        if (sc->f)
          en += sc->f(i, j, i + 1, j - 1, VRNA_DECOMP_PAIR_ML, sc->data) +
                sc->f(i + 1, j - 1, i + 1, *l, VRNA_DECOMP_ML_ML, sc->data);
      }

      if (e == en) {
        *k  = i;
        *l  = j;
        return 1;
      }
    }

    u++;
  }

  return 0;
}


PRIVATE int
backtrack_f3(vrna_fold_compound_t *fc,
             unsigned int         *k,
             unsigned int         *i,
             unsigned int         *j,
             int                  *f3)
{
  short         *S, *S1, s5, s3;
  unsigned int  ii, n, u, *sn, type;
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

  if (ii >= n) {
    /* no more pairs */
    *k = *i = *j = 0;
    return 1;
  }

  /* ii is paired, find pairing partner */
  switch (dangle_model) {
    case 0:
      for (u = ii + 1; u <= n; u++) {
        if ((sn[u] == sn[u + 1]) &&
            (hc->mx[ii * n + u] & VRNA_CONSTRAINT_CONTEXT_EXT_LOOP)) {
          type  = vrna_get_ptype_md(S[ii], S[u], md);
          en    = c[idx[u] + ii];

          if (sc)
            if (sc->f)
              en += sc->f(ii, n, u, u + 1, VRNA_DECOMP_EXT_STEM_EXT, sc->data);

          if (fij == vrna_E_exterior_stem(type, -1, -1, P) + en + f3[u + 1]) {
            *i  = ii;
            *j  = u;
            *k  = u + 1;
            return 1;
          }
        }
      }
      break;

    case 2:
      for (u = ii + 1; u <= n; u++) {
        if ((sn[u] == sn[u + 1]) &&
            (hc->mx[ii * n + u] & VRNA_CONSTRAINT_CONTEXT_EXT_LOOP)) {
          type  = vrna_get_ptype_md(S[ii], S[u], md);
          s5    = S1[ii - 1];
          s3    = (u < n) ? S1[u + 1] : -1;
          en    = c[idx[u] + ii];

          if (sc)
            if (sc->f)
              en += sc->f(ii, n, u, u + 1, VRNA_DECOMP_EXT_STEM_EXT, sc->data);

          if (fij == vrna_E_exterior_stem(type, s5, s3, P) + en + f3[u + 1]) {
            *i  = ii;
            *j  = u;
            *k  = u + 1;
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
