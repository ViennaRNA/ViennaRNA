/* Last changed Time-stamp: <2007-08-26 11:59:45 ivo> */
/*
           compute the duplex structure of two RNA strands,
                allowing only inter-strand base pairs.
         see cofold() for computing hybrid structures without
                             restriction.

                             Ivo Hofacker
                          Vienna RNA package
*/

#include <config.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <string.h>
#include "utils.h"
#include "energy_par.h"
#include "fold_vars.h"
#include "fold.h"
#include "pair_mat.h"
#include "params.h"
#include "alifold.h"
#include "subopt.h"
#include "loop_energies.h"
#include "duplex.h"

#ifdef _OPENMP
#include <omp.h>
#endif

/*@unused@*/
static char rcsid[] UNUSED = "$Id: duplex.c,v 1.8 2007/08/26 10:08:44 ivo Exp $";

#define STACK_BULGE1  1     /* stacking energies for bulges of size 1 */
#define NEW_NINIO     1     /* new asymetry penalty */
#define MAXSECTORS    500   /* dimension for a backtrack array */
#define LOCALITY      0.    /* locality parameter for base-pairs */
#define UNIT 100
#define MINPSCORE -2 * UNIT
#define NONE -10000         /* score for forbidden pairs */



/*
#################################
# GLOBAL VARIABLES              #
#################################
*/
/*
#################################
# PRIVATE VARIABLES             #
#################################
*/
PRIVATE paramT  *P  = NULL;
PRIVATE int     **c = NULL;                  /* energy array, given that i-j pair */
PRIVATE short   *S1 = NULL, *SS1 = NULL, *S2 = NULL, *SS2 = NULL;
PRIVATE int     n1,n2;                /* sequence lengths */

#ifdef _OPENMP

/* NOTE: all variables are assumed to be uninitialized if they are declared as threadprivate
         thus we have to initialize them before usage by a seperate function!
         OR: use copyin in the PARALLEL directive!
         e.g.:
         #pragma omp parallel for copyin(P)
*/
#pragma omp threadprivate(P, c, S1, SS1, S2, SS2, n1, n2)

#endif

/*
#################################
# PRIVATE FUNCTION DECLARATIONS #
#################################
*/
PRIVATE duplexT duplexfold_cu(const char *s1, const char *s2, int clean_up);
PRIVATE duplexT aliduplexfold_cu(const char *s1[], const char *s2[], int clean_up);
PRIVATE char    *backtrack(int i, int j);
PRIVATE char    *alibacktrack(int i, int j, const short **S1, const short **S2);
PRIVATE int     compare(const void *sub1, const void *sub2);
PRIVATE int     covscore(const int *types, int n_seq);

/*
#################################
# BEGIN OF FUNCTION DEFINITIONS #
#################################
*/

PUBLIC duplexT duplexfold(const char *s1, const char *s2){
  return duplexfold_cu(s1, s2, 1);
}

PRIVATE duplexT duplexfold_cu(const char *s1, const char *s2, int clean_up){
  int i, j, l1, Emin=INF, i_min=0, j_min=0;
  char *struc;
  duplexT mfe;

  n1 = (int) strlen(s1);
  n2 = (int) strlen(s2);

  if ((!P) || (fabs(P->temperature - temperature)>1e-6)) {
    if(P) free(P); P = scale_parameters();
    make_pair_matrix();
  }

  c = (int **) space(sizeof(int *) * (n1+1));
  for (i=1; i<=n1; i++) c[i] = (int *) space(sizeof(int) * (n2+1));

  S1  = encode_sequence(s1, 0);
  S2  = encode_sequence(s2, 0);
  SS1 = encode_sequence(s1, 1);
  SS2 = encode_sequence(s2, 1);

  for (i=1; i<=n1; i++) {
    for (j=n2; j>0; j--) {
      int type, type2, E, k,l;
      type = pair[S1[i]][S2[j]];
      c[i][j] = type ? P->DuplexInit : INF;
      if (!type) continue;
      c[i][j] += E_ExtLoop(type, (i>1) ? SS1[i-1] : -1, (j<n2) ? SS2[j+1] : -1, P);
      for (k=i-1; k>0 && k>i-MAXLOOP-2; k--) {
        for (l=j+1; l<=n2; l++) {
          if (i-k+l-j-2>MAXLOOP) break;
          type2 = pair[S1[k]][S2[l]];
          if (!type2) continue;
          E = E_IntLoop(i-k-1, l-j-1, type2, rtype[type],
                            SS1[k+1], SS2[l-1], SS1[i-1], SS2[j+1], P);
          c[i][j] = MIN2(c[i][j], c[k][l]+E);
        }
      }
      E = c[i][j];
      E += E_ExtLoop(rtype[type], (j > 1) ? SS2[j-1] : -1, (i<n1) ? SS1[i+1] : -1, P);
      if (E<Emin) {
        Emin=E; i_min=i; j_min=j;
      }
    }
  }

  struc = backtrack(i_min, j_min);
  if (i_min<n1) i_min++;
  if (j_min>1 ) j_min--;
  l1 = strchr(struc, '&')-struc;
  /*
    printf("%s %3d,%-3d : %3d,%-3d (%5.2f)\n", struc, i_min+1-l1, i_min,
       j_min, j_min+strlen(struc)-l1-2, Emin*0.01);
  */
  mfe.i = i_min;
  mfe.j = j_min;
  mfe.energy = (float) Emin/100.;
  mfe.structure = struc;
  if(clean_up) {
    for (i=1; i<=n1; i++) free(c[i]);
    free(c);
    free(S1);
    free(S2);
    free(SS1);
    free(SS2);
  }
  return mfe;
}

PUBLIC duplexT *duplex_subopt(const char *s1, const char *s2, int delta, int w) {
  int i,j, n1, n2, thresh, E, n_subopt=0, n_max;
  char *struc;
  duplexT mfe;
  duplexT *subopt;

  n_max=16;
  subopt = (duplexT *) space(n_max*sizeof(duplexT));
  mfe = duplexfold_cu(s1, s2, 0);
  free(mfe.structure);

  thresh = (int) mfe.energy*100+0.1 + delta;
  n1 = strlen(s1); n2=strlen(s2);
  for (i=n1; i>0; i--) {
    for (j=1; j<=n2; j++) {
      int type, ii,jj, Ed;
      type = pair[S2[j]][S1[i]];
      if (!type) continue;
      E = Ed = c[i][j];
      Ed += E_ExtLoop(type, (j>1) ? SS2[j-1] : -1, (i<n1) ? SS1[i+1] : -1, P);
      if (Ed>thresh) continue;
      /* too keep output small, remove hits that are dominated by a
         better one close (w) by. For simplicity we do test without
         adding dangles, which is slightly inaccurate.
      */
      for (ii=MAX2(i-w,1); (ii<=MIN2(i+w,n1)) && type; ii++) {
        for (jj=MAX2(j-w,1); jj<=MIN2(j+w,n2); jj++)
          if (c[ii][jj]<E) {type=0; break;}
      }
      if (!type) continue;

      struc = backtrack(i,j);
      fprintf(stderr, "%d %d %d\n", i,j,E);
      if (n_subopt+1>=n_max) {
        n_max *= 2;
        subopt = (duplexT *) xrealloc(subopt, n_max*sizeof(duplexT));
      }
      subopt[n_subopt].i = MIN2(i+1,n1);
      subopt[n_subopt].j = MAX2(j-1,1);
      subopt[n_subopt].energy = Ed * 0.01;
      subopt[n_subopt++].structure = struc;
    }
  }
  /* free all static globals */
  for (i=1; i<=n1; i++) free(c[i]);
  free(c);
  free(S1); free(S2); free(SS1); free(SS2);

  if (subopt_sorted) qsort(subopt, n_subopt, sizeof(duplexT), compare);
  subopt[n_subopt].i =0;
  subopt[n_subopt].j =0;
  subopt[n_subopt].structure = NULL;
  return subopt;
}

PRIVATE char *backtrack(int i, int j) {
  /* backtrack structure going backwards from i, and forwards from j
     return structure in bracket notation with & as separator */
  int k, l, type, type2, E, traced, i0, j0;
  char *st1, *st2, *struc;

  st1 = (char *) space(sizeof(char)*(n1+1));
  st2 = (char *) space(sizeof(char)*(n2+1));

  i0=MIN2(i+1,n1); j0=MAX2(j-1,1);

  while (i>0 && j<=n2) {
    E = c[i][j]; traced=0;
    st1[i-1] = '(';
    st2[j-1] = ')';
    type = pair[S1[i]][S2[j]];
    if (!type) nrerror("backtrack failed in fold duplex");
    for (k=i-1; k>0 && k>i-MAXLOOP-2; k--) {
      for (l=j+1; l<=n2; l++) {
        int LE;
        if (i-k+l-j-2>MAXLOOP) break;
        type2 = pair[S1[k]][S2[l]];
        if (!type2) continue;
        LE = E_IntLoop(i-k-1, l-j-1, type2, rtype[type],
                       SS1[k+1], SS2[l-1], SS1[i-1], SS2[j+1], P);
        if (E == c[k][l]+LE) {
          traced=1;
          i=k; j=l;
          break;
        }
      }
      if (traced) break;
    }
    if (!traced) {
      E -= E_ExtLoop(type, (i>1) ? SS1[i-1] : -1, (j<n2) ? SS2[j+1] : -1, P);
      if (E != P->DuplexInit) {
        nrerror("backtrack failed in fold duplex");
      } else break;
    }
  }
  if (i>1)  i--;
  if (j<n2) j++;

  struc = (char *) space(i0-i+1+j-j0+1+2);
  for (k=MAX2(i,1); k<=i0; k++) if (!st1[k-1]) st1[k-1] = '.';
  for (k=j0; k<=j; k++) if (!st2[k-1]) st2[k-1] = '.';
  strcpy(struc, st1+MAX2(i-1,0)); strcat(struc, "&");
  strcat(struc, st2+j0-1);

  /* printf("%s %3d,%-3d : %3d,%-3d\n", struc, i,i0,j0,j);  */
  free(st1); free(st2);

  return struc;
}

/*------------------------------------------------------------------------*/

PRIVATE int compare(const void *sub1, const void *sub2) {
  int d;
  if (((duplexT *) sub1)->energy > ((duplexT *) sub2)->energy)
    return 1;
  if (((duplexT *) sub1)->energy < ((duplexT *) sub2)->energy)
    return -1;
  d = ((duplexT *) sub1)->i - ((duplexT *) sub2)->i;
  if (d!=0) return d;
  return  ((duplexT *) sub1)->j - ((duplexT *) sub2)->j;
}

/*---------------------------------------------------------------------------*/

PUBLIC duplexT aliduplexfold(const char *s1[], const char *s2[]){
  return aliduplexfold_cu(s1, s2, 1);
}

PRIVATE duplexT aliduplexfold_cu(const char *s1[], const char *s2[], int clean_up) {
  int i, j, s, n_seq, l1, Emin=INF, i_min=0, j_min=0;
  char *struc;
  duplexT mfe;
  short **S1, **S2;
  int *type;
  n1 = (int) strlen(s1[0]);
  n2 = (int) strlen(s2[0]);

  for (s=0; s1[s]!=NULL; s++);
  n_seq = s;
  for (s=0; s2[s]!=NULL; s++);
  if (n_seq != s) nrerror("unequal number of sequences in aliduplexfold()\n");

  if ((!P) || (fabs(P->temperature - temperature)>1e-6)) {
    if(P) free(P); P = scale_parameters();
    make_pair_matrix();
  }

  c = (int **) space(sizeof(int *) * (n1+1));
  for (i=1; i<=n1; i++) c[i] = (int *) space(sizeof(int) * (n2+1));

  S1 = (short **) space((n_seq+1)*sizeof(short *));
  S2 = (short **) space((n_seq+1)*sizeof(short *));
  for (s=0; s<n_seq; s++) {
    if (strlen(s1[s]) != n1) nrerror("uneqal seqence lengths");
    if (strlen(s2[s]) != n2) nrerror("uneqal seqence lengths");
    S1[s] = encode_sequence(s1[s], 0);
    S2[s] = encode_sequence(s2[s], 0);
  }
  type = (int *) space(n_seq*sizeof(int));

  for (i=1; i<=n1; i++) {
    for (j=n2; j>0; j--) {
      int k,l,E,psc;
      for (s=0; s<n_seq; s++) {
        type[s] = pair[S1[s][i]][S2[s][j]];
      }
      psc = covscore(type, n_seq);
      for (s=0; s<n_seq; s++) if (type[s]==0) type[s]=7;
      c[i][j] = (psc>=MINPSCORE) ? (n_seq*P->DuplexInit) : INF;
      if (psc<MINPSCORE) continue;
      for(s=0; s<n_seq;s++){
        c[i][j] += E_ExtLoop(type[s], (i>1) ? S1[s][i-1] : -1, (j<n2) ? S2[s][j+1] : -1, P);
      }
      for (k=i-1; k>0 && k>i-MAXLOOP-2; k--) {
        for (l=j+1; l<=n2; l++) {
          int type2;
          if (i-k+l-j-2>MAXLOOP) break;
          if (c[k][l]>INF/2) continue;
          for (E=s=0; s<n_seq; s++) {
            type2 = pair[S1[s][k]][S2[s][l]];
            if (type2==0) type2=7;
            E += E_IntLoop(i-k-1, l-j-1, type2, rtype[type[s]],
                           S1[s][k+1], S2[s][l-1], S1[s][i-1], S2[s][j+1], P);
          }
          c[i][j] = MIN2(c[i][j], c[k][l]+E);
        }
      }
      c[i][j] -= psc;
      E = c[i][j];
      for (s=0; s<n_seq; s++) {
        E += E_ExtLoop(rtype[type[s]], (j>1) ? S2[s][j-1] : -1, (i<n1) ? S1[s][i+1] : -1, P);
      }
      if (E<Emin) {
        Emin=E; i_min=i; j_min=j;
      }
    }
  }

  struc = alibacktrack(i_min, j_min, (const short **)S1,(const short **)S2);
  if (i_min<n1) i_min++;
  if (j_min>1 ) j_min--;
  l1 = strchr(struc, '&')-struc;
  /*
    printf("%s %3d,%-3d : %3d,%-3d (%5.2f)\n", struc, i_min+1-l1, i_min,
       j_min, j_min+strlen(struc)-l1-2, Emin*0.01);
  */
  mfe.i = i_min;
  mfe.j = j_min;
  mfe.energy = (float) (Emin/(100.*n_seq));
  mfe.structure = struc;
  if (clean_up){
    for (i=1; i<=n1; i++) free(c[i]);
    free(c);
  }
  for (s=0; s<n_seq; s++) {
    free(S1[s]); free(S2[s]);
  }
  free(S1);
  free(S2);
  free(type);
  return mfe;
}

PUBLIC duplexT *aliduplex_subopt(const char *s1[], const char *s2[], int delta, int w) {
  int i,j, n1, n2, thresh, E, n_subopt=0, n_max, s, n_seq, *type;
  char *struc;
  duplexT mfe;
  duplexT *subopt;
  short **S1, **S2;

  n_max=16;
  subopt = (duplexT *) space(n_max*sizeof(duplexT));
  mfe = aliduplexfold_cu(s1, s2, 0);
  free(mfe.structure);

  for (s=0; s1[s]!=NULL; s++);
  n_seq = s;

  thresh =  (int) ((mfe.energy*100. + delta)*n_seq +0.1);
  n1 = strlen(s1[0]); n2=strlen(s2[0]);
  S1 = (short **) space((n_seq+1)*sizeof(short *));
  S2 = (short **) space((n_seq+1)*sizeof(short *));
  for (s=0; s<n_seq; s++) {
    if (strlen(s1[s]) != n1) nrerror("uneqal seqence lengths");
    if (strlen(s2[s]) != n2) nrerror("uneqal seqence lengths");
    S1[s] = encode_sequence(s1[s], 0);
    S2[s] = encode_sequence(s2[s], 0);
  }
  type = (int *) space(n_seq*sizeof(int));

  for (i=n1; i>0; i--) {
    for (j=1; j<=n2; j++) {
      int ii, jj, skip, Ed, psc;

      for (s=0; s<n_seq; s++) {
        type[s] = pair[S2[s][j]][S1[s][i]];
      }
      psc = covscore(type, n_seq);
      for (s=0; s<n_seq; s++) if (type[s]==0) type[s]=7;
      if (psc<MINPSCORE) continue;
      E = Ed = c[i][j];
      for  (s=0; s<n_seq; s++) {
        Ed += E_ExtLoop(type[s], (j>1) ? S2[s][j-1] : -1, (i<n1) ? S1[s][i+1] : -1, P);
      }
      if (Ed>thresh) continue;
      /* too keep output small, skip hits that are dominated by a
         better one close (w) by. For simplicity we don't take dangels
         into account here, thus the heuristic is somewhat inaccurate.
      */
      for (skip=0, ii=MAX2(i-w,1); (ii<=MIN2(i+w,n1)) && type; ii++) {
        for (jj=MAX2(j-w,1); jj<=MIN2(j+w,n2); jj++)
          if (c[ii][jj]<E) {skip=1; break;}
      }
      if (skip) continue;
      struc = alibacktrack(i,j,(const short **)S1, (const short **)S2);
      fprintf(stderr, "%d %d %d\n", i,j,E);
      if (n_subopt+1>=n_max) {
        n_max *= 2;
        subopt = (duplexT *) xrealloc(subopt, n_max*sizeof(duplexT));
      }
      subopt[n_subopt].i = MIN2(i+1,n1);
      subopt[n_subopt].j = MAX2(j-1,1);
      subopt[n_subopt].energy = Ed * 0.01/n_seq;
      subopt[n_subopt++].structure = struc;
    }
  }

  for (i=1; i<=n1; i++) free(c[i]);
  free(c);
  for (s=0; s<n_seq; s++) {
    free(S1[s]); free(S2[s]);
  }
  free(S1); free(S2); free(type);

  if (subopt_sorted) qsort(subopt, n_subopt, sizeof(duplexT), compare);
  subopt[n_subopt].i =0;
  subopt[n_subopt].j =0;
  subopt[n_subopt].structure = NULL;
  return subopt;
}

PRIVATE char *alibacktrack(int i, int j, const short **S1, const short **S2) {
  /* backtrack structure going backwards from i, and forwards from j
     return structure in bracket notation with & as separator */
  int k, l, *type, type2, E, traced, i0, j0, s, n_seq;
  char *st1, *st2, *struc;

  n1 = (int) S1[0][0];
  n2 = (int) S2[0][0];

  for (s=0; S1[s]!=NULL; s++);
  n_seq = s;
  for (s=0; S2[s]!=NULL; s++);
  if (n_seq != s) nrerror("unequal number of sequences in alibacktrack()\n");

  st1 = (char *) space(sizeof(char)*(n1+1));
  st2 = (char *) space(sizeof(char)*(n2+1));
  type = (int *) space(n_seq*sizeof(int));

  i0=MIN2(i+1,n1); j0=MAX2(j-1,1);

  while (i>0 && j<=n2) {
    int psc;
    E = c[i][j]; traced=0;
    st1[i-1] = '(';
    st2[j-1] = ')';
    for (s=0; s<n_seq; s++) {
      type[s] = pair[S1[s][i]][S2[s][j]];
    }
    psc = covscore(type, n_seq);
    for (s=0; s<n_seq; s++) if (type[s]==0) type[s] = 7;
    E += psc;
    for (k=i-1; k>0 && k>i-MAXLOOP-2; k--) {
      for (l=j+1; l<=n2; l++) {
        int LE;
        if (i-k+l-j-2>MAXLOOP) break;
        if (c[k][l]>INF/2) continue;
        for (s=LE=0; s<n_seq; s++) {
          type2 = pair[S1[s][k]][S2[s][l]];
          if (type2==0) type2=7;
          LE += E_IntLoop(i-k-1, l-j-1, type2, rtype[type[s]],
                           S1[s][k+1], S2[s][l-1], S1[s][i-1], S2[s][j+1], P);
        }
        if (E == c[k][l]+LE) {
          traced=1;
          i=k; j=l;
          break;
        }
      }
      if (traced) break;
    }
    if (!traced) {
      for (s=0; s<n_seq; s++) {
        E -= E_ExtLoop(type[s], (i>1) ? S1[s][i-1] : -1, (j<n2) ? S2[s][j+1] : -1, P);
      }
      if (E != n_seq*P->DuplexInit) {
        nrerror("backtrack failed in aliduplex");
      } else break;
    }
  }
  if (i>1)  i--;
  if (j<n2) j++;

  struc = (char *) space(i0-i+1+j-j0+1+2);
  for (k=MAX2(i,1); k<=i0; k++) if (!st1[k-1]) st1[k-1] = '.';
  for (k=j0; k<=j; k++) if (!st2[k-1]) st2[k-1] = '.';
  strcpy(struc, st1+MAX2(i-1,0)); strcat(struc, "&");
  strcat(struc, st2+j0-1);

  /* printf("%s %3d,%-3d : %3d,%-3d\n", struc, i,i0,j0,j);  */
  free(st1); free(st2); free(type);

  return struc;
}


PRIVATE int covscore(const int *types, int n_seq) {
  /* calculate co-variance bonus for a pair depending on  */
  /* compensatory/consistent mutations and incompatible seqs */
  /* should be 0 for conserved pairs, >0 for good pairs      */
  int k,l,s,score, pscore;
  int dm[7][7]={{0,0,0,0,0,0,0}, /* hamming distance between pairs */
                {0,0,2,2,1,2,2} /* CG */,
                {0,2,0,1,2,2,2} /* GC */,
                {0,2,1,0,2,1,2} /* GU */,
                {0,1,2,2,0,2,1} /* UG */,
                {0,2,2,1,2,0,2} /* AU */,
                {0,2,2,2,1,2,0} /* UA */};

  int pfreq[8]={0,0,0,0,0,0,0,0};
  for (s=0; s<n_seq; s++)
    pfreq[types[s]]++;

  if (pfreq[0]*2>n_seq) return NONE;
  for (k=1,score=0; k<=6; k++) /* ignore pairtype 7 (gap-gap) */
    for (l=k+1; l<=6; l++)
      /* scores for replacements between pairtypes    */
      /* consistent or compensatory mutations score 1 or 2  */
      score += pfreq[k]*pfreq[l]*dm[k][l];

  /* counter examples score -1, gap-gap scores -0.25   */
  pscore = cv_fact *
    ((UNIT*score)/n_seq - nc_fact*UNIT*(pfreq[0] + pfreq[7]*0.25));
  return pscore;
}
