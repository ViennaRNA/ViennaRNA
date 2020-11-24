/*
 *                             MEA.c
 *               c  Ivo L Hofacker, Vienna RNA package
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>    /* #defines DBL_EPSILON */
#include <math.h>
#include "ViennaRNA/fold_vars.h"
#include "ViennaRNA/utils/basic.h"
#include "ViennaRNA/params/basic.h"
#include "ViennaRNA/gquad.h"
#include "ViennaRNA/alphabet.h"
#include "ViennaRNA/MEA.h"

/* compute an MEA structure, i.e. the structure maximising
 * EA = \sum_{(i,j) \in S} 2\gamma p_{i,j} + \sum_{i is unpaired} p^u_i
 *
 * This can be computed by a variant of the Nussinov recursion:
 * M(i,j) = min(M(i,j-1)+pu[j], min_k M(i,k-1)+C(k,j)
 * C(i,j) = 2*gamma*p_ij + M(i+1,j-1)
 *
 * Just for fun, we implement it as a sparse DP algorithm.
 * At any time we store only the current and previous row of M.
 * The C matrix is implemented as a sparse matrix:
 * For each j we store in C[j] a list of values (i, MEA([i..j])), containing
 * the MEA over all structures closed by (i,j).
 * The list is sparse since only C values where C(i,j)==M(i,j) can
 * contribute to the optimal solution.
 */

typedef struct Litem {
  int     i;
  double  A;
} Litem;

typedef struct List {
  size_t  size;  /* allocated space */
  size_t  nelem;
  Litem   *list;
} List;

/*
 #################################
 # PRIVATE FUNCTION DECLARATIONS #
 #################################
 */

PRIVATE int
comp_plist(const void *a,
           const void *b);


PRIVATE vrna_ep_t *
prune_sort(vrna_ep_t    *p,
           double       *pu,
           unsigned int n,
           double       gamma,
           short        *S,
           int          gq);


PRIVATE void
pushC(List    *c,
      int     i,
      double  a);


struct MEAdat {
  vrna_ep_t *pl;
  double    *pu;
  double    gamma;
  List      *C;
  double    *Mi;
  char      *structure;
};

PRIVATE void
mea_backtrack(const struct MEAdat *bdat,
              int                 i,
              int                 j,
              int                 paired,
              short               *S,
              vrna_exp_param_t    *pf);


PRIVATE float
compute_MEA(vrna_ep_t         *p,
            unsigned int      n,
            short             *S,
            double            gamma,
            vrna_exp_param_t  *pf,
            char              *structure);


/*
 #################################
 # BEGIN OF FUNCTION DEFINITIONS #
 #################################
 */
PUBLIC char *
vrna_MEA(vrna_fold_compound_t *fc,
         double               gamma,
         float                *mea)
{
  char      *structure;
  int       gq;
  vrna_ep_t *pl;

  structure = NULL;

  if ((fc) &&
      (mea) &&
      (fc->exp_params) &&
      (fc->exp_matrices) &&
      (fc->exp_matrices->probs)) {
    gq = fc->exp_params->model_details.gquad;

    structure = (char *)vrna_alloc(sizeof(char) * (fc->length + 1));

    fc->exp_params->model_details.gquad = 0;
    pl                                  = vrna_plist_from_probs(fc, 1e-4 / (1 + gamma));
    fc->exp_params->model_details.gquad = gq;

    *mea = compute_MEA(pl,
                       fc->length,
                       (fc->type == VRNA_FC_TYPE_SINGLE) ? fc->sequence_encoding : fc->S_cons,
                       gamma,
                       fc->exp_params,
                       structure);

    free(pl);
  }

  return structure;
}


PUBLIC char *
vrna_MEA_from_plist(vrna_ep_t   *plist,
                    const char  *sequence,
                    double      gamma,
                    vrna_md_t   *md_p,
                    float       *mea)
{
  char              *structure;
  short             *S;
  unsigned int      n;
  vrna_md_t         md;
  vrna_exp_param_t  *exp_params;

  structure = NULL;

  if ((plist) &&
      (sequence) &&
      (mea)) {
    n         = strlen(sequence);
    structure = (char *)vrna_alloc(sizeof(char) * (n + 1));

    if (md_p)
      md = *md_p;
    else
      vrna_md_set_default(&md);

    exp_params = vrna_exp_params(&md);

    S = vrna_seq_encode(sequence, &md);

    *mea = compute_MEA(plist,
                       n,
                       S,
                       gamma,
                       exp_params,
                       structure);

    free(S);
    free(exp_params);
  }

  return structure;
}


/*
 #####################################
 # BEGIN OF STATIC HELPER FUNCTIONS  #
 #####################################
 */
PRIVATE float
compute_MEA(vrna_ep_t         *p,
            unsigned int      n,
            short             *S,
            double            gamma,
            vrna_exp_param_t  *pf,
            char              *structure)
{
  unsigned int  i, j;
  int           with_gquad = 0;
  double        EA, MEA, *Mi, *Mi1, *tmp, *pu;
  vrna_ep_t     *pp, *pl;
  vrna_md_t     *md;
  Litem         *li;
  List          *C;
  struct MEAdat bdat;

  md          = &(pf->model_details);
  with_gquad  = md->gquad;

  memset(structure, '.', sizeof(char) * n);
  structure[n] = '\0';

  pu  = vrna_alloc(sizeof(double) * (n + 1));
  pp  = pl = prune_sort(p, pu, n, gamma, S, with_gquad);

  C = (List *)vrna_alloc((n + 1) * (sizeof(List)));

  Mi  = (double *)vrna_alloc((n + 1) * sizeof(double));
  Mi1 = (double *)vrna_alloc((n + 1) * sizeof(double));

  for (i = n; i > 0; i--) {
    Mi[i] = pu[i];
    for (j = i + 1; j <= n; j++) {
      Mi[j] = Mi[j - 1] + pu[j];
      for (li = C[j].list; li < C[j].list + C[j].nelem; li++) {
        EA    = li->A + Mi[(li->i) - 1];
        Mi[j] = MAX2(Mi[j], EA);
      }
      if ((unsigned int)pp->i == i && (unsigned int)pp->j == j) {
        EA = 2 * gamma * pp->p + Mi1[j - 1];
        if (Mi[j] < EA) {
          Mi[j] = EA;
          pushC(&C[j], i, EA); /* only push into C[j] list if optimal */
        }

        pp++;
      }
    }
    tmp = Mi1;
    Mi1 = Mi;
    Mi  = tmp;
  }

  MEA = Mi1[n];

  bdat.structure  = structure;
  bdat.gamma      = gamma;
  bdat.C          = C;
  bdat.Mi         = Mi1;
  bdat.pl         = pl;
  bdat.pu         = pu;
  mea_backtrack(&bdat, 1, n, 0, S, pf);
  free(Mi);
  free(Mi1);
  free(pl);
  free(pu);
  for (i = 1; i <= n; i++)
    if (C[i].list)
      free(C[i].list);

  free(C);

  return MEA;
}


PRIVATE int
comp_plist(const void *a,
           const void *b)
{
  vrna_ep_t *A, *B;
  int       di;

  A   = (vrna_ep_t *)a;
  B   = (vrna_ep_t *)b;
  di  = (B->i - A->i);
  if (di != 0)
    return di;

  return A->j - B->j;
}


PRIVATE vrna_ep_t *
prune_sort(vrna_ep_t    *p,
           double       *pu,
           unsigned int n,
           double       gamma,
           short        *S,
           int          gq)
{
  /*
   * produce a list containing all base pairs with
   * 2*gamma*p_ij > p^u_i + p^u_j
   * already sorted to be in the order we need them within the DP
   */
  unsigned int  size, i, nump = 0;
  vrna_ep_t     *pp, *pc, *pc2;

  for (i = 1; i <= n; i++)
    pu[i] = 1.;

  for (pc = p; pc->i > 0; pc++) {
    if (pc->type == VRNA_PLIST_TYPE_BASEPAIR) {
      pu[pc->i] -= pc->p;
      pu[pc->j] -= pc->p;
    }
  }

  if (gq) {
    if (!S)
      vrna_message_error("no sequence information available in MEA gquad!");

    /* remove probabilities that i or j are enclosed by a gquad */
    for (i = 1; i <= n; i++) {
      for (pc2 = p; pc2->i > 0; pc2++) {
        /* skip all non-gquads */
        if (S[pc2->i] != 3)
          continue;

        if (S[pc2->j] != 3)
          continue;

        /* remove only if i is enclosed */
        if ((pc2->i < i) && (pc2->j > i))
          pu[i] -= pc2->p;
      }
    }
  }

  size  = n + 1;
  pp    = vrna_alloc(sizeof(vrna_ep_t) * (n + 1));
  for (pc = p; pc->i > 0; pc++) {
    if (pc->i > n)
      vrna_message_error("mismatch between vrna_ep_t and structure in MEA()");

    if (pc->type == VRNA_PLIST_TYPE_BASEPAIR) {
      if (pc->p * 2 * gamma > pu[pc->i] + pu[pc->j]) {
        if (nump + 1 >= size) {
          size  += size / 2 + 1;
          pp    = vrna_realloc(pp, size * sizeof(vrna_ep_t));
        }

        pp[nump++] = *pc;
      }
    }
  }
  pp[nump].i = pp[nump].j = pp[nump].p = 0;
  qsort(pp, nump, sizeof(vrna_ep_t), comp_plist);
  return pp;
}


PRIVATE void
pushC(List    *c,
      int     i,
      double  a)
{
  if (c->nelem + 1 >= c->size) {
    c->size = MAX2(8, c->size * sqrt(2));
    c->list = vrna_realloc(c->list, sizeof(Litem) * c->size);
  }

  c->list[c->nelem].i = i;
  c->list[c->nelem].A = a;
  c->nelem++;
}


PRIVATE void
mea_backtrack(const struct MEAdat *bdat,
              int                 i,
              int                 j,
              int                 pair,
              short               *S,
              vrna_exp_param_t    *pf)
{
  /*
   * backtrack structure for the interval [i..j]
   * recursively calls itself, recomputes the necessary parts of the M matrix
   */
  int     fail, gq, k, L, l[3];
  double  *Mi, prec, *pu, EA;
  List    *C;
  Litem   *li;

  fail  = 1;
  gq    = pf->model_details.gquad;
  C     = bdat->C;
  Mi    = bdat->Mi;
  pu    = bdat->pu;

  if (pair) {
    /*
     * if pair == 1, insert pair and re-compute Mi values
     * else Mi is already filled
     */
    if (gq) {
      if ((S[i] == 3) && (S[j] == 3)) {
        get_gquad_pattern_pf(S, i, j, pf, &L, l);
        for (k = 0; k < L; k++) {
          bdat->structure[i + k - 1] \
                  = bdat->structure[i + k + L + l[0] - 1] \
                  = bdat->structure[i + k + 2 * L + l[0] + l[1] - 1] \
                  = bdat->structure[i + k + 3 * L + l[0] + l[1] + l[2] - 1] \
                  = '+';
        }
        return;
      } else {
        bdat->structure[i - 1]  = '(';
        bdat->structure[j - 1]  = ')';
        i++;
        j--;
        /* We've done this before in MEA() but didn't keep the results */
        Mi[i - 1] = 0;
        Mi[i]     = pu[i];
        for (k = i + 1; k <= j; k++) {
          Mi[k] = Mi[k - 1] + pu[k];
          for (li = C[k].list; li < C[k].list + C[k].nelem && li->i >= i; li++) {
            EA    = li->A + Mi[(li->i) - 1];
            Mi[k] = MAX2(Mi[k], EA);
          }
        }
      }
    } else {
      bdat->structure[i - 1]  = '(';
      bdat->structure[j - 1]  = ')';
      i++;
      j--;
      /* We've done this before in MEA() but didn't keep the results */
      Mi[i - 1] = 0;
      Mi[i]     = pu[i];
      for (k = i + 1; k <= j; k++) {
        Mi[k] = Mi[k - 1] + pu[k];
        for (li = C[k].list; li < C[k].list + C[k].nelem && li->i >= i; li++) {
          EA    = li->A + Mi[(li->i) - 1];
          Mi[k] = MAX2(Mi[k], EA);
        }
      }
    }
  }

  prec = DBL_EPSILON * Mi[j];
  /* Mi values are filled, do the backtrace */
  while (j > i && Mi[j] <= Mi[j - 1] + pu[j] + prec) {
    bdat->structure[j - 1] = '.';
    j--;
  }
  for (li = C[j].list; li < C[j].list + C[j].nelem && li->i >= i; li++) {
    if (Mi[j] <= li->A + Mi[(li->i) - 1] + prec) {
      if (li->i > i + 3)
        mea_backtrack(bdat, i, (li->i) - 1, 0, S, pf);

      mea_backtrack(bdat, li->i, j, 1, S, pf);
      fail = 0;
    }
  }
  if (fail && j > i)
    vrna_message_error("backtrack failed for MEA()");
}


#ifndef VRNA_DISABLE_BACKWARD_COMPATIBILITY

/*
 *###########################################
 *# deprecated functions below              #
 *###########################################
 */
PUBLIC float
MEA(vrna_ep_t *p,
    char      *structure,
    double    gamma)
{
  return MEA_seq(p, NULL, structure, gamma, NULL);
}


PUBLIC float
MEA_seq(vrna_ep_t         *p,
        const char        *sequence,
        char              *structure,
        double            gamma,
        vrna_exp_param_t  *pf)
{
  short             *S;
  double            MEA;
  vrna_exp_param_t  *exp_params;
  vrna_md_t         md;

  S = NULL;

  if (!pf) {
    /* use global variables to set model details */
    set_model_details(&md);
    exp_params = vrna_exp_params(&md);
  } else {
    exp_params = pf;
  }

  if (sequence)
    S = vrna_seq_encode(sequence, &(exp_params->model_details));

  MEA = compute_MEA(p,
                    strlen(structure),
                    S,
                    gamma,
                    exp_params,
                    structure);

  /* clean up */
  free(S);
  if (!pf)
    free(exp_params);

  return MEA;
}


#endif
