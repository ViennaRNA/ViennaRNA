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
#include "duplex.h"
/*@unused@*/
static char rcsid[] UNUSED = "$Id: duplex.c,v 1.8 2007/08/26 10:08:44 ivo Exp $";

#define PUBLIC
#define PRIVATE static

#define STACK_BULGE1  1   /* stacking energies for bulges of size 1 */
#define NEW_NINIO     1   /* new asymetry penalty */

PRIVATE void  encode_seqs(const char *s1, const char *s2);
PRIVATE short *encode_seq(const char *seq);
PRIVATE char * backtrack(int i, int j);
PRIVATE char * alibacktrack(int i, int j, const short *S1[], const short *S2[]);
PRIVATE int compare(const void *sub1, const void *sub2);
PRIVATE int covscore(const int *types, int n_seq);

extern int subopt_sorted; /* from subopt.c, default 0 */

extern double cv_fact; /* from alifold.c, default 1 */
extern double nc_fact;

/*@unused@*/

#define MAXSECTORS      500     /* dimension for a backtrack array */
#define LOCALITY        0.      /* locality parameter for base-pairs */

#define MIN2(A, B)      ((A) < (B) ? (A) : (B))
#define MAX2(A, B)      ((A) > (B) ? (A) : (B))

PRIVATE paramT *P = NULL;

PRIVATE int   **c;      /* energy array, given that i-j pair */
PRIVATE short  *S1, *SS1, *S2, *SS2;
PRIVATE int   n1,n2;    /* sequence lengths */
extern  int  LoopEnergy(int n1, int n2, int type, int type_2,
                        int si1, int sj1, int sp1, int sq1);

PRIVATE int delay_free=0;
/*--------------------------------------------------------------------------*/


duplexT duplexfold(const char *s1, const char *s2) {
  int i, j, l1, Emin=INF, i_min=0, j_min=0;
  char *struc;
  duplexT mfe;
  
  n1 = (int) strlen(s1);
  n2 = (int) strlen(s2);
  
  if ((!P) || (fabs(P->temperature - temperature)>1e-6)) {
    update_fold_params();  P = scale_parameters();
    make_pair_matrix();
  }
  
  c = (int **) space(sizeof(int *) * (n1+1));
  for (i=1; i<=n1; i++) c[i] = (int *) space(sizeof(int) * (n2+1));
  
  encode_seqs(s1, s2);
  
  for (i=1; i<=n1; i++) {
    for (j=n2; j>0; j--) {
      int type, type2, E, k,l;
      type = pair[S1[i]][S2[j]];
      c[i][j] = type ? P->DuplexInit : INF;
      if (!type) continue;
      if (i>1)  c[i][j] += P->dangle5[type][SS1[i-1]];
      if (j<n2) c[i][j] += P->dangle3[type][SS2[j+1]];
      if (type>2) c[i][j] += P->TerminalAU;
      for (k=i-1; k>0 && k>i-MAXLOOP-2; k--) {
	for (l=j+1; l<=n2; l++) {
	  if (i-k+l-j-2>MAXLOOP) break;
	  type2 = pair[S1[k]][S2[l]];
	  if (!type2) continue;
	  E = LoopEnergy(i-k-1, l-j-1, type2, rtype[type],
			    SS1[k+1], SS2[l-1], SS1[i-1], SS2[j+1]);
	  c[i][j] = MIN2(c[i][j], c[k][l]+E);
	}
      }
      E = c[i][j]; 
      if (i<n1) E += P->dangle3[rtype[type]][SS1[i+1]];
      if (j>1)  E += P->dangle5[rtype[type]][SS2[j-1]];
      if (type>2) E += P->TerminalAU;
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
  if (!delay_free) {
    for (i=1; i<=n1; i++) free(c[i]);
    free(c);
    free(S1); free(S2); free(SS1); free(SS2);
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
  delay_free=1;
  mfe = duplexfold(s1, s2);
  free(mfe.structure);
  
  thresh = (int) mfe.energy*100+0.1 + delta;
  n1 = strlen(s1); n2=strlen(s2);
  for (i=n1; i>0; i--) {
    for (j=1; j<=n2; j++) {
      int type, ii,jj, Ed;
      type = pair[S2[j]][S1[i]];
      if (!type) continue;
      E = Ed = c[i][j];
      if (i<n1) Ed += P->dangle3[type][SS1[i+1]];
      if (j>1)  Ed += P->dangle5[type][SS2[j-1]];
      if (type>2) Ed += P->TerminalAU;
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
  
  for (i=1; i<=n1; i++) free(c[i]);
  free(c);
  free(S1); free(S2); free(SS1); free(SS2);
  delay_free=0;

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
	LE = LoopEnergy(i-k-1, l-j-1, type2, rtype[type],
		       SS1[k+1], SS2[l-1], SS1[i-1], SS2[j+1]);
	if (E == c[k][l]+LE) {
	  traced=1; 
	  i=k; j=l;
	  break;
	}
      }
      if (traced) break;
    }
    if (!traced) { 
      if (i>1) E -= P->dangle5[type][SS1[i-1]];
      if (j<n2) E -= P->dangle3[type][SS2[j+1]];
      if (type>2) E -= P->TerminalAU;
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

/*---------------------------------------------------------------------------*/

PRIVATE short * encode_seq(const char *sequence) {
  unsigned int i,l;
  short *S;
  l = strlen(sequence);
  S = (short *) space(sizeof(short)*(l+2));
  S[0] = (short) l;

  /* make numerical encoding of sequence */
  for (i=1; i<=l; i++)
    S[i]= (short) encode_char(toupper(sequence[i-1]));

  /* for circular folding add first base at position n+1 */
  S[l+1] = S[1];

  return S;
}

PRIVATE void encode_seqs(const char *s1, const char *s2) {
  unsigned int i,l;

  l = strlen(s1);
  S1 = encode_seq(s1);
  SS1= (short *) space(sizeof(short)*(l+1));
  /* SS1 exists only for the special X K and I bases and energy_set!=0 */
  
  for (i=1; i<=l; i++) { /* make numerical encoding of sequence */
    SS1[i] = alias[S1[i]];   /* for mismatches of nostandard bases */
  }

  l = strlen(s2);
  S2 = encode_seq(s2);
  SS2= (short *) space(sizeof(short)*(l+1));
  /* SS2 exists only for the special X K and I bases and energy_set!=0 */
  
  for (i=1; i<=l; i++) { /* make numerical encoding of sequence */
    SS2[i] = alias[S2[i]];   /* for mismatches of nostandard bases */
  }
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

#define UNIT 100
#define MINPSCORE -2 * UNIT

duplexT aliduplexfold(const char *s1[], const char *s2[]) {
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
    update_fold_params();  P = scale_parameters();
    make_pair_matrix();
  }
  
  c = (int **) space(sizeof(int *) * (n1+1));
  for (i=1; i<=n1; i++) c[i] = (int *) space(sizeof(int) * (n2+1));
  
  S1 = (short **) space((n_seq+1)*sizeof(short *));
  S2 = (short **) space((n_seq+1)*sizeof(short *));
  for (s=0; s<n_seq; s++) {
    if (strlen(s1[s]) != n1) nrerror("uneqal seqence lengths");
    if (strlen(s2[s]) != n2) nrerror("uneqal seqence lengths");
    S1[s] = encode_seq(s1[s]);
    S2[s] = encode_seq(s2[s]);
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
      if (i>1) 
	for (s=0; s<n_seq; s++) 
	  c[i][j] += P->dangle5[type[s]][S1[s][i-1]];
      if (j<n2) 
	for (s=0; s<n_seq; s++) 
	  c[i][j] += P->dangle3[type[s]][S2[s][j+1]];
      for (s=0; s<n_seq; s++) 
	if (type[s]>2) c[i][j] += P->TerminalAU;
      for (k=i-1; k>0 && k>i-MAXLOOP-2; k--) {
	for (l=j+1; l<=n2; l++) {
	  int type2;
	  if (i-k+l-j-2>MAXLOOP) break;
	  if (c[k][l]>INF/2) continue;
	  for (E=s=0; s<n_seq; s++) { 
	    type2 = pair[S1[s][k]][S2[s][l]];
	    if (type2==0) type2=7;
	    E += LoopEnergy(i-k-1, l-j-1, type2, rtype[type[s]],
			   S1[s][k+1], S2[s][l-1], S1[s][i-1], S2[s][j+1]);
	  }
	  c[i][j] = MIN2(c[i][j], c[k][l]+E);
	}
      }
      c[i][j] -= psc;
      E = c[i][j]; 
      for (s=0; s<n_seq; s++) {
	if (i<n1) E += P->dangle3[rtype[type[s]]][S1[s][i+1]];
	if (j>1)  E += P->dangle5[rtype[type[s]]][S2[s][j-1]];
	if (type[s]>2) E += P->TerminalAU;
      }
      if (E<Emin) {
	Emin=E; i_min=i; j_min=j;
      } 
    }
  }
  
  struc = alibacktrack(i_min, j_min, S1, S2);
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
  if (!delay_free) {
    for (i=1; i<=n1; i++) free(c[i]);
    free(c);
  }
  for (s=0; s<n_seq; s++) {
    free(S1[s]); free(S2[s]);
  }
  free(S1); free(S2); free(type);
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
  delay_free=1;
  mfe = aliduplexfold(s1, s2);
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
    S1[s] = encode_seq(s1[s]);
    S2[s] = encode_seq(s2[s]);
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
	if (i<n1) Ed += P->dangle3[type[s]][S1[s][i+1]];
	if (j>1)  Ed += P->dangle5[type[s]][S2[s][j-1]];
	if (type[s]>2) Ed += P->TerminalAU;
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
      struc = alibacktrack(i,j,S1,S2);
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
  delay_free=0;

  if (subopt_sorted) qsort(subopt, n_subopt, sizeof(duplexT), compare);
  subopt[n_subopt].i =0;
  subopt[n_subopt].j =0;
  subopt[n_subopt].structure = NULL;
  return subopt;
}

PRIVATE char *alibacktrack(int i, int j, const short *S1[], const short *S2[]) {
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
	  LE += LoopEnergy(i-k-1, l-j-1, type2, rtype[type[s]],
			   S1[s][k+1], S2[s][l-1], S1[s][i-1], S2[s][j+1]);
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
	if (i>1)  E -= P->dangle5[type[s]][S1[s][i-1]];
	if (j<n2) E -= P->dangle3[type[s]][S2[s][j+1]];
	if (type[s]>2) E -= P->TerminalAU;
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
#define NONE -10000 /* score for forbidden pairs */
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
