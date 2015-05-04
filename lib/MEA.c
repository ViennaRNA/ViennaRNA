/*
			       MEA.c
		 c  Ivo L Hofacker, Vienna RNA package
*/
/* Last changed Time-stamp: <2009-06-18 14:04:21 ivo> */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>    /* #defines DBL_EPSILON */
#include <math.h>
#include "fold_vars.h"
#include "utils.h"

/* compute an MEA structure, i.e. the structure maximising
   EA = \sum_{(i,j) \in S} 2\gamma p_{i,j} + \sum_{i is unpaired} p^u_i

   This can be computed by a variant of the Nussinov recursion:
   M(i,j) = min(M(i,j-1)+pu[j], min_k M(i,k-1)+C(k,j)
   C(i,j) = 2*gamma*p_ij + M(i+1,j-1)

   Just for fun, we implement it as a sparse DP algorithm.
   At any time we store only the current and previous row of M.
   The C matrix is implemented as a sparse matrix:
   For each j we store in C[j] a list of values (i, MEA([i..j])), containing
   the MEA over all structures closed by (i,j).
   The list is sparse since only C values where C(i,j)==M(i,j) can
   contribute to the optimal solution.
*/

typedef struct Litem {
  int i;
  double A;
} Litem;

typedef struct List {
  unsigned size;   /* allocated space */
  unsigned nelem;
  Litem *list;
} List;

PRIVATE int comp_plist(const void *a, const void *b);
PRIVATE plist *prune_sort(plist *p, double *pu, int n, double gamma);
PRIVATE void pushC(List *c, int i, double a);

struct MEAdat{
  plist *pl;
  double *pu;
  double gamma;
  List *C;
  double *Mi;
  char * structure;
};

PRIVATE void mea_backtrack(const struct MEAdat *bdat, int i, int j, int paired);

float MEA(plist *p, char *structure, double gamma) {

  int i,j,n;
  Litem *li;
  plist *pp, *pl;

  List *C;
  double MEA, *Mi, *Mi1, *tmp, *pu;
  struct MEAdat bdat;

  n = strlen(structure);
  for (i=0; i<n; i++) structure[i] = '.';

  pu = space(sizeof(double)*(n+1));
  pp = pl = prune_sort(p, pu, n, gamma);

  C = (List*) space((n+1)*(sizeof(List)));

  Mi = (double *) space((n+1)*sizeof(double));
  Mi1 = (double *) space((n+1)*sizeof(double));

  for (i=n; i>0; i--) {
    Mi[i] = pu[i];
    for (j=i+1; j<=n; j++) {
      double EA;
      Mi[j] = Mi[j-1] + pu[j];
      for (li=C[j].list; li<C[j].list+C[j].nelem; li++) {
	EA = li->A + Mi[(li->i) -1];
	Mi[j] = MAX2(Mi[j], EA);
      }
      if (pp->i == i && pp->j ==j) {
	EA = 2*gamma*pp->p +  Mi1[j-1];
	if (Mi[j]<EA) {
	  Mi[j]=EA;
	  pushC(&C[j], i, EA); /* only push into C[j] list if optimal */
	}
	pp++;
      }

    }
    tmp = Mi1; Mi1 = Mi; Mi = tmp;
  }

  MEA = Mi1[n];

  bdat.structure = structure; bdat.gamma = gamma;
  bdat.C = C;  bdat.Mi=Mi1; bdat.pl=pl; bdat.pu = pu;
  mea_backtrack(&bdat, 1, n, 0);
  free(Mi); free(Mi1); free(pl); free(pu);
  for (i=1; i<=n; i++)
    if (C[i].list) free(C[i].list);
  free(C);
  return MEA;
}

PRIVATE int comp_plist(const void *a, const void *b) {
  plist *A, *B;
  int di;
  A = (plist *)a;
  B = (plist *)b;
  di = (B->i - A->i);
  if (di!=0) return di;
  return (A->j - B->j);
}


PRIVATE plist *prune_sort(plist *p, double *pu, int n, double gamma) {
  /*
     produce a list containing all base pairs with
     2*gamma*p_ij > p^u_i + p^u_j
     already sorted to be in the order we need them within the DP
  */
  unsigned size, i, nump = 0;
  plist *pp, *pc;

  for (i=1; i<=n; i++) pu[i]=1.;

  for (pc=p; pc->i >0; pc++) {
    pu[pc->i] -= pc->p;
    pu[pc->j] -= pc->p;
  }
  size = n+1;
  pp = space(sizeof(plist)*(n+1));
  for (pc=p; pc->i >0; pc++) {
    if (pc->i > n) nrerror("mismatch between plist and structure in MEA()");
    if (pc->p*2*gamma > pu[pc->i] + pu[pc->j]) {
      if (nump+1 >= size) {
	size += size/2 + 1;
	pp = xrealloc(pp, size*sizeof(plist));
      }
      pp[nump++] = *pc;
    }
  }
  pp[nump].i = pp[nump].j = pp[nump].p = 0;
  qsort(pp, nump, sizeof(plist), comp_plist);
  return pp;
}

PRIVATE void pushC(List *c, int i, double a) {
  if (c->nelem+1>=c->size) {
    c->size = MAX2(8,c->size*sqrt(2));
    c->list = xrealloc(c->list, sizeof(Litem)*c->size);
  }
  c->list[c->nelem].i = i;
  c->list[c->nelem].A = a;
  c->nelem++;
}

PRIVATE void mea_backtrack(const struct MEAdat *bdat, int i, int j, int pair) {
  /* backtrack structure for the interval [i..j] */
  /* recursively calls itself, recomputes the necessary parts of the M matrix */
  List *C; Litem *li;
  double *Mi, prec;
  double *pu;
  int fail=1;

  C = bdat->C;
  Mi = bdat->Mi;
  pu = bdat->pu;

  if (pair) {
    int k;
    /* if pair == 1, insert pair and re-compute Mi values */
    /* else Mi is already filled */
    bdat->structure[i-1] = '(';
    bdat->structure[j-1] = ')';
    i++; j--;
    /* We've done this before in MEA() but didn't keep the results */
    Mi[i-1]=0; Mi[i]=pu[i];
    for (k=i+1; k<=j; k++) {
      Mi[k] = Mi[k-1] + pu[k];
      for (li=C[k].list; li<C[k].list+C[k].nelem && li->i >= i; li++) {
	double EA;
	EA = li->A + Mi[(li->i) -1];
	Mi[k] = MAX2(Mi[k], EA);
      }
    }
  }

  prec = DBL_EPSILON * Mi[j];
  /* Mi values are filled, do the backtrace */
  while (j>i && Mi[j] <= Mi[j-1] + pu[j] + prec) {
    bdat->structure[j-1]='.';
    j--;
  }
  for (li=C[j].list; li<C[j].list + C[j].nelem && li->i >= i; li++) {
    double EA;
    if (Mi[j] <= li->A + Mi[(li->i) -1] + prec) {
      if (li->i > i+3) mea_backtrack(bdat, i, (li->i)-1, 0);
      mea_backtrack(bdat, li->i, j, 1);
      fail = 0;
    }
  }
  if (fail && j>i) nrerror("backtrack failed for MEA()");
}
