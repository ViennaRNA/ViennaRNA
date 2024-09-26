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
#include "ViennaRNA/utils/log.h"
#include "ViennaRNA/params/basic.h"
#include "ViennaRNA/partfunc/gquad.h"
#include "ViennaRNA/sequences/alphabet.h"
#include "ViennaRNA/datastructures/array.h"
#include "ViennaRNA/datastructures/sparse_mx.h"
#include "ViennaRNA/structures/mea.h"

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
  unsigned int  i;
  double        A;
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


PRIVATE int
comp_plist_circ(const void  *a,
                const void  *b);


PRIVATE vrna_ep_t *
prune_plist(vrna_ep_t     *p,
            double        *pu,
            unsigned int  *size,
            unsigned int  n,
            double        gamma,
            unsigned int  *features);


PRIVATE void
pushC(List          *c,
      unsigned int  i,
      double        a);


struct MEAdat {
  vrna_ep_t *pl;
  double    *pu;
  double    gamma;
  List      *C;
  double    *Mi;
  char      *structure;
};

PRIVATE void
mea_backtrack(vrna_fold_compound_t  *fc,
              const struct MEAdat   *bdat,
              unsigned int          i,
              unsigned int          j,
              unsigned int          paired);


PRIVATE float
compute_MEA(vrna_fold_compound_t  *fc,
            vrna_ep_t             *p,
            double                gamma,
            char                  *structure);


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
  vrna_ep_t *pl;

  structure = NULL;

  if ((fc) &&
      (mea) &&
      (fc->exp_params) &&
      (fc->exp_matrices) &&
      (fc->exp_matrices->probs)) {
    structure = (char *)vrna_alloc(sizeof(char) * (fc->length + 1));
    pl        = vrna_plist_from_probs(fc, 1e-4 / (1 + gamma));

    *mea = compute_MEA(fc,
                       pl,
                       gamma,
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
  char                  *structure;
  unsigned int          n;
  vrna_md_t             md;
  vrna_exp_param_t      *exp_params;
  vrna_fold_compound_t  *fc;

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

    fc = vrna_fold_compound(sequence, &md, VRNA_OPTION_DEFAULT | VRNA_OPTION_EVAL_ONLY);
    vrna_exp_params_subst(fc, exp_params);

    *mea = compute_MEA(fc,
                       plist,
                       gamma,
                       structure);

    vrna_fold_compound_free(fc);
    free(exp_params);
  }

  return structure;
}


/*
 #####################################
 # BEGIN OF STATIC HELPER FUNCTIONS  #
 #####################################
 */

#define HAS_CIRC_GQUAD    1U


PRIVATE float
compute_MEA(vrna_fold_compound_t  *fc,
            vrna_ep_t             *p,
            double                gamma,
            char                  *structure)
{
  unsigned int  i, j, n, size_pl, features, jmax_gq, i_gq, j_gq;
  double        EA, MEA, *Mi, *Mi1, *Mi1_gq, *tmp, *pu;
  vrna_ep_t     *pp, *pl;
  Litem         *li;
  List          *C;
  struct MEAdat bdat;

  vrna_smx_csr(double)  *gq_circ      = NULL;
  vrna_array(unsigned int)  gq_circ_imin  = NULL;
  vrna_array(unsigned int)  gq_circ_imax  = NULL;

  n       = fc->length;
  jmax_gq = 0;
  i_gq    = j_gq = 0;
  Mi1_gq  = NULL;

  memset(structure, '.', sizeof(char) * n);
  structure[n] = '\0';

  pu  = vrna_alloc(sizeof(double) * (n + 1));
  pp  = pl = prune_plist(p, pu, &size_pl, n, gamma, &features);

  /* check whether we have G-Quadruplexes that span the n,1 junction */
  if (features & HAS_CIRC_GQUAD) {
    /* sort entries so that we can safely insert them into the sparse mx */
    qsort(pl, size_pl, sizeof(vrna_ep_t), comp_plist_circ);

    /* required for backtracking */
    Mi1_gq = (double *)vrna_alloc((n + 1) * sizeof(double));

    vrna_array_init_size(gq_circ_imin, n + 1);
    for (i = 0; i <= n; i++)
      vrna_array_append(gq_circ_imin, n);

    vrna_array_init_size(gq_circ_imax, n + 1);
    memset(gq_circ_imax, 0, sizeof(unsigned int) * (n + 1));
    vrna_array_size(gq_circ_imax) = n + 1;

    for (pp = pl; pp->i > 0; pp++) {
      if ((pp->type == VRNA_PLIST_TYPE_GQUAD) &&
          (pp->i > pp->j)) {
        if (gq_circ == NULL)
          gq_circ = vrna_smx_csr_double_init(n);

        EA = (pp->j + n - pp->i + 1) * gamma * pp->p;
#ifndef VRNA_DISABLE_C11_FEATURES
        vrna_smx_csr_insert(gq_circ, pp->j, pp->i, EA);
#else
        vrna_smx_csr_double_insert(gq_circ, pp->j, pp->i, EA);
#endif

        if (pp->j > jmax_gq)
          jmax_gq = pp->j;

        if (pp->i < gq_circ_imin[pp->j])
          gq_circ_imin[pp->j] = pp->i;

        if (pp->i > gq_circ_imax[pp->j])
          gq_circ_imax[pp->j] = pp->i;
      }
    }
  }

  /* sort in the order we need it for the dp recursions below */
  qsort(pl, size_pl, sizeof(vrna_ep_t), comp_plist);

  pp = pl;

  C = (List *)vrna_alloc((n + 1) * (sizeof(List)));

  Mi  = (double *)vrna_alloc((n + 1) * sizeof(double));
  Mi1 = (double *)vrna_alloc((n + 1) * sizeof(double));

  MEA = 0.;

  for (i = n; i > 0; i--) {
    Mi[i] = pu[i];
    if ((pp->i == i) &&
        (pp->i > pp->j))
      pp++;

    for (j = i + 1; j <= n; j++) {
      Mi[j] = Mi[j - 1] + pu[j];

      for (li = C[j].list; li < C[j].list + C[j].nelem; li++) {
        EA    = li->A + Mi[(li->i) - 1];
        Mi[j] = MAX2(Mi[j], EA);
      }

      if (((unsigned int)pp->i == i) &&
          ((unsigned int)pp->j == j)) {
        EA = Mi1[j - 1];

        switch (pp->type) {
          case VRNA_PLIST_TYPE_GQUAD:
            EA += (j - i + 1) * gamma * pp->p;
            break;

          case VRNA_PLIST_TYPE_BASEPAIR:
            EA += 2 * gamma * pp->p;
            break;

          default:
            /* do nothing */
            break;
        }

        if (Mi[j] < EA) {
          Mi[j] = EA;
          pushC(&C[j], i, EA); /* only push into C[j] list if optimal */
        }

        pp++;
      }

      if ((gq_circ) &&
          (i <= jmax_gq) &&
          (j >= gq_circ_imin[i]) &&
          (j <= gq_circ_imax[i])) {
#ifndef VRNA_DISABLE_C11_FEATURES
        if ((EA = vrna_smx_csr_get(gq_circ, i, j, 0.)) != 0.) {
#else
        if ((EA = vrna_smx_csr_double_get(gq_circ, i, j, 0.)) != 0.) {
#endif
          EA += Mi1[j - 1];
          if (MEA < EA) {
            MEA   = EA;
            i_gq  = i;
            j_gq  = j;
            /* make a copy of Mi1 */
            Mi1_gq = memcpy(Mi1_gq, Mi1, sizeof(double) * (n + 1));
          }
        }
      }
    }
    tmp = Mi1;
    Mi1 = Mi;
    Mi  = tmp;
  }

  MEA = MAX2(MEA, Mi1[n]);

  bdat.structure  = structure;
  bdat.gamma      = gamma;
  bdat.pl         = pl;
  bdat.pu         = pu;
  bdat.C          = C;

  if ((i_gq > 0) &&
      (j_gq > 0)) {
    unsigned int L, l[3];
    vrna_get_gquad_pattern_pf(fc, j_gq, i_gq, &L, l);
    if (L > 0)
      vrna_db_insert_gq(structure, j_gq, L, l, n);

    bdat.Mi = Mi1_gq;
    mea_backtrack(fc, &bdat, i_gq + 1, j_gq - 1, 0);
  } else {
    bdat.Mi = Mi1;
    mea_backtrack(fc, &bdat, 1, n, 0);
  }

  vrna_array_free(gq_circ_imin);
  vrna_array_free(gq_circ_imax);
#ifndef VRNA_DISABLE_C11_FEATURES
  vrna_smx_csr_free(gq_circ);
#else
  vrna_smx_csr_double_free(gq_circ);
#endif
  free(Mi);
  free(Mi1);
  free(Mi1_gq);
  free(pl);
  free(pu);
  for (i = 1; i <= n; i++)
    if (C[i].list)
      free(C[i].list);

  free(C);

  return MEA;
}


/*
 *  sort by sequence position:
 *  1. in descending order for i
 *  2. in ascending order for j
 */
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


/*
 *  sort by sequence position:
 *  1. in ascending order for j
 *  2. in ascending order for i
 */
PRIVATE int
comp_plist_circ(const void  *a,
                const void  *b)
{
  vrna_ep_t *A, *B;
  int       di;

  A   = (vrna_ep_t *)a;
  B   = (vrna_ep_t *)b;
  di  = (A->j - B->j);
  if (di != 0)
    return di;

  return A->i - B->i;
}


PRIVATE vrna_ep_t *
prune_plist(vrna_ep_t     *p,
            double        *pu,
            unsigned int  *size_pl,
            unsigned int  n,
            double        gamma,
            unsigned int  *features)
{
  /*
   * produce a list containing all base pairs with
   * 2*gamma*p_ij > p^u_i + p^u_j
   * already sorted to be in the order we need them within the DP
   */
  unsigned int  size, i, nt;
  double        pug;
  vrna_ep_t     *pp, *pc;

  *size_pl  = 0;
  *features = 0;

  for (i = 1; i <= n; i++)
    pu[i] = 1.;

  /* collect probabilities to be unpaired */
  for (pc = p; pc->i > 0; pc++) {
    switch (pc->type) {
      case VRNA_PLIST_TYPE_GQUAD:
        if (pc->i < pc->j) {
          for (i = pc->i; i <= pc->j; i++)
            pu[i] -= pc->p;
        } else {
          for (i = pc->i; i <= n; i++)
            pu[i] -= pc->p;
          for (i = 1; i <= pc->j; i++)
            pu[i] -= pc->p;
        }

        break;

      case VRNA_PLIST_TYPE_BASEPAIR:
        pu[pc->i] -= pc->p;
        pu[pc->j] -= pc->p;
        break;

      default:
        /* do nothing */
        break;
    }
  }

  /*
   * check whether any unpaired probabilities should be explicitely
   * overwritten due to input data
   */
  for (pc = p; pc->i > 0; pc++)
    if (pc->type == VRNA_PLIST_TYPE_UNPAIRED)
      for (i = pc->i; i <= pc->j; i++)
        pu[i] = pc->p;

  /* prepare pair probability entries */
  /* here, we only keep pair/gquad probabilities that exceed
   * the sum of unpaired probabilities of the nucleotides they
   * cover, respectively.
   */
  size  = n + 1;
  pp    = vrna_alloc(sizeof(vrna_ep_t) * (n + 1));
  for (pc = p; pc->i > 0; pc++) {
    if (pc->i > n) {
      vrna_log_error("mismatch between vrna_ep_t and structure in MEA()");
      return NULL;
    }

    nt  = 0;
    pug = 0.;

    switch (pc->type) {
      case VRNA_PLIST_TYPE_GQUAD:
        if (pc->i < pc->j) {
          nt = pc->j - pc->i + 1;
          for (i = pc->i; i <= pc->j; i++)
            pug += pu[i];
        } else {
          nt = pc->j + n - pc->i + 1;
          for (i = pc->i; i <= n; i++)
            pug += pu[i];
          for (i = 1; i <= pc->j; i++)
            pug += pu[i];
        }

        break;

      case VRNA_PLIST_TYPE_BASEPAIR:
        nt  = 2;
        pug = pu[pc->i] + pu[pc->j];
        break;
    }

    if (pc->p * nt * gamma > pug) {
      if ((*size_pl) + 1 >= size) {
        size  += size / 2 + 1;
        pp    = vrna_realloc(pp, size * sizeof(vrna_ep_t));
      }

      if ((pc->type == VRNA_PLIST_TYPE_GQUAD) &&
          (pc->i > pc->j))
        *features |= HAS_CIRC_GQUAD;

      pp[(*size_pl)++] = *pc;
    }
  }
  pp[(*size_pl)].i = pp[(*size_pl)].j = pp[(*size_pl)].p = 0;
  return pp;
}


PRIVATE void
pushC(List          *c,
      unsigned int  i,
      double        a)
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
mea_backtrack(vrna_fold_compound_t  *fc,
              const struct MEAdat   *bdat,
              unsigned int          i,
              unsigned int          j,
              unsigned int          pair)
{
  /*
   * backtrack structure for the interval [i..j]
   * recursively calls itself, recomputes the necessary parts of the M matrix
   */
  short             *S;
  unsigned int      fail, L, l[3], n, gq, k;
  double            *Mi, prec, *pu, EA;
  List              *C;
  Litem             *li;
  vrna_exp_param_t  *pf;

  n     = fc->length;
  pf    = fc->exp_params;
  fail  = 1;
  gq    = pf->model_details.gquad;
  C     = bdat->C;
  Mi    = bdat->Mi;
  pu    = bdat->pu;
  S     = (fc->type == VRNA_FC_TYPE_SINGLE) ? fc->sequence_encoding : fc->S_cons;

  if (pair) {
    /*
     * if pair == 1, insert pair and re-compute Mi values
     * else Mi is already filled
     */
    if (gq) {
      if ((S[i] == 3) && (S[j] == 3)) {
        vrna_get_gquad_pattern_pf(fc, i, j, &L, l);
        if (L > 0)
          vrna_db_insert_gq(bdat->structure, i, L, l, n);
        else
          vrna_log_error("Failed to parse G-Quadruplex");

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
  while ((j > i) &&
         (Mi[j] <= (Mi[j - 1] + pu[j] + prec))) {
    bdat->structure[j - 1] = '.';
    j--;
  }

  for (li = C[j].list; li < C[j].list + C[j].nelem && li->i >= i; li++) {
    if (Mi[j] <= li->A + Mi[(li->i) - 1] + prec) {
      if (li->i > i + 3)
        mea_backtrack(fc, bdat, i, (li->i) - 1, 0);

      mea_backtrack(fc, bdat, li->i, j, 1);
      fail = 0;
    }
  }
  if (fail && j > i)
    vrna_log_error("backtrack failed for MEA()");
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
  char                  *s;
  unsigned int          n;
  double                MEA;
  vrna_fold_compound_t  *fc;

  s   = NULL;
  fc  = NULL;
  MEA = 0.;

  if ((p) &&
      (structure)) {
    n = strlen(structure);
    if (sequence) {
      if (strlen(sequence) != n) {
        vrna_log_error("sequence and structure of different length (%u vs. %u)",
                       strlen(sequence),
                       n);
        return 0.;
      }
    } else {
      s = (char *)vrna_alloc(sizeof(char) * (n + 1));
      memset(s, '.', sizeof(char) * n);
      s[n]      = '\0';
      sequence  = (const char *)s;
    }

    fc = vrna_fold_compound(sequence, NULL, VRNA_OPTION_DEFAULT | VRNA_OPTION_EVAL_ONLY);
    if (pf)
      vrna_exp_params_subst(fc, pf);

    MEA = compute_MEA(fc,
                      p,
                      gamma,
                      structure);

    /* clean up */
    vrna_fold_compound_free(fc);
    free(s);
  }

  return MEA;
}


#endif
