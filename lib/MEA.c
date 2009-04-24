/*
			       MEA.c
		 c  Ivo L Hofacker, Vienna RNA package
*/
/* Last changed Time-stamp: <2009-04-21 11:02:39 ivo> */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>    /* #defines DBL_EPSILON */
#include "utils.h"
#include "PS_dot.h"  /* defines plist */

/* compute an MEA structure, i.e. the structure maximising
   EA = \sum_{(i,j) \in S} p_{i,j} + \sum_{i is unpaired} \gamma

   This can be computed by a vaiant of the Nussinov recursion:
   M(i,j) = min(M(i,j-1)+gamma, min_k M(i,k-1)+C(k,j)
   C(i,j) = p_ij + M(i+1,j-1)

   Just for fun, we implement it as a sparse DP algorithm.
   At any time we store only the current and previous row of M.
   The C matrix is implemented as a sparse matrix:
   For each j we store in C[j] a list of values (i, MEA([i..j])), containing
   the MEA over all structures closed by (i,j).
*/
#define MAX2(A, B)      ((A) > (B) ? (A) : (B))

typedef struct Litem {
  int i;
  double A;
} Litem;

static int comp_plist(const void *a, const void *b);
static plist *prune_sort(plist *p, int n, double gamma);
static void pushC(Litem *c, int i, double a);

struct MEAdat{
  plist *pl;
  double gamma;
  Litem **C;
  double *Mi;
  char * structure;
};

static void mea_backtrack(const struct MEAdat *bdat, int i, int j, int paired);

float MEA(plist *p, char *structure, double gamma) {

  int i,j,n;
  Litem *li;
  plist *pp, *pl;

  Litem **C;
  double MEA, *Mi, *Mi1, *tmp;
  struct MEAdat bdat;

  n = strlen(structure);
  for (i=0; i<n; i++) structure[i] = '.';

  pp = pl = prune_sort(p, n, gamma);

  C = (Litem **) space((n+1)*(sizeof(Litem *)));
  for (i=1; i<=n; i++)
    C[i] = (Litem *) space(sizeof(Litem)*MAX2(n, (int) (1/(2*gamma))+1));
  Mi = (double *) space((n+1)*sizeof(double));
  Mi1 = (double *) space((n+1)*sizeof(double));

  for (i=n; i>0; i--) {
    Mi[i] = gamma;
    for (j=i+1; j<=n; j++) {
      double EA;
      Mi[j] = Mi[j-1] + gamma;
      for (li=C[j]; li->i > 0; li++) {
	EA = li->A + Mi[(li->i) -1];
	Mi[j] = MAX2(Mi[j], EA);
      }
      if (pp->i == i && pp->j ==j) {
	EA = pp->p +  Mi1[j-1];
	if (Mi[j]<EA) {
	  Mi[j]=EA;
	  pushC(C[j], pp->i, EA); /* only push into C[j] list if optimal */
	}
	pp++;
      }

    }
    tmp = Mi1; Mi1 = Mi; Mi = tmp;
  }

  MEA = Mi1[n];

  bdat.structure = structure; bdat.gamma = gamma;
  bdat.C = C;  bdat.Mi=Mi1; bdat.pl=pl;
  mea_backtrack(&bdat, 1, n, 0);
  free(Mi); free(Mi1); free(pl);
  for (i=1; i<=n; i++) free(C[i]);
  free(C);
  return MEA;
}

static int comp_plist(const void *a, const void *b) {
  plist *A, *B;
  int di;
  A = (plist *)a;
  B = (plist *)b;
  di = (B->i - A->i);
  if (di!=0) return di;
  return (A->j - B->j);
}


static plist *prune_sort(plist *p, int n, double gamma) {
  unsigned size, nump = 0;
  plist *pp, *pc;

  size = n+1;
  pp = space(sizeof(plist)*(n+1));
  for (pc=p; pc->i >0; pc++) {
    if (pc->i > n) nrerror("mismatch between plist and structure in MEA()");
    if (pc->p > 2*gamma) {
      if (nump+1 >= size) {
	size += size/2 + 1;
	pp = xrealloc(pp, size*sizeof(plist));
      }
      pp[nump++] = *pc;
    }
  }
  qsort(pp, nump, sizeof(plist), comp_plist);
  return pp;
}

static void pushC(Litem *c, int i, double a) {
  Litem *ci;
  for (ci=c; ci->i >0; ci++); /* search for the end (inefficient, but list is O(1)) */
  ci->i = i;
  ci->A = a;
}

static void mea_backtrack(const struct MEAdat *bdat, int i, int j, int pair) {

  Litem **C, *li;
  double *Mi, prec;
  int fail=1;

  /* fprintf(stderr, "%4d %4d %d\n", i,j,pair); */
  C = bdat->C;
  Mi = bdat->Mi;

  if (pair) {
    int k;
    /* if pair == 1, insert pair and re-compute Mi values */
    /* else Mi is already filled */
    bdat->structure[i-1] = '(';
    bdat->structure[j-1] = ')';
    i++; j--;
    /* We've done this before in MEA() but didn't keep the results */
    Mi[i-1]=0; Mi[i]=bdat->gamma;
    for (k=i+1; k<=j; k++) {
      Mi[k] = Mi[k-1] + bdat->gamma;
      for (li=C[k]; li->i >=i; li++) {
	double EA;
	EA = li->A + Mi[(li->i) -1];
	if (Mi[k] < EA) Mi[k] = EA;
      }
    }
  }

  prec = DBL_EPSILON * Mi[j];
  /* Mi values are filled, do the backtrace */
  while (j>i && Mi[j] <= Mi[j-1] + bdat->gamma + prec) {
    bdat->structure[j-1]='.';
    j--;
  }
  for (li=C[j]; li->i >= i; li++) {
    double EA;
    if (Mi[j] <= li->A + Mi[(li->i) -1] + prec) {
      if (li->i > i+3) mea_backtrack(bdat, i, (li->i)-1, 0);
      mea_backtrack(bdat, li->i, j, 1);
      fail = 0;
    }
  }
  if (fail && j>i) nrerror("backtrack failed for MEA()");
}
