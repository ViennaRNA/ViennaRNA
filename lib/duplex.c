/* Last changed Time-stamp: <2005-07-23 16:45:18 ivo> */
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
#include "pair_mat.h"
#include "params.h"
#include "duplex.h"
/*@unused@*/
static char rcsid[] UNUSED = "$Id: duplex.c,v 1.4 2005/07/24 08:37:40 ivo Exp $";

#define PUBLIC
#define PRIVATE static

#define STACK_BULGE1  1   /* stacking energies for bulges of size 1 */
#define NEW_NINIO     1   /* new asymetry penalty */

PRIVATE void  encode_seq(const char *s1, const char *s2);
PRIVATE char * backtrack(int i, int j);
PRIVATE void update_dfold_params(void);
PRIVATE int compare(const void *sub1, const void *sub2);

extern int subopt_sorted; /* from subopt.c */

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
  
  if ((!P) || (fabs(P->temperature - temperature)>1e-6))
    update_dfold_params();
  
  c = (int **) space(sizeof(int *) * (n1+1));
  for (i=1; i<=n1; i++) c[i] = (int *) space(sizeof(int) * (n2+1));
  
  encode_seq(s1, s2);
  
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
      if (j>1) E += P->dangle5[rtype[type]][SS2[j-1]];
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
      int type;
      type = pair[S2[j]][S1[i]];
      if (!type) continue;
      E = c[i][j];
      if (i<n1) E += P->dangle3[type][SS1[i+1]];
      if (j>1)  E += P->dangle5[type][SS2[j-1]];
      if (type>2) E += P->TerminalAU;
      if (E<=thresh) {
	struc = backtrack(i,j);
#if 0
	int l1;
	l1 = strchr(struc, '&')-struc;
	printf("%s %3d,%-3d : %3d,%-3d (%5.2f)\n", struc, i+2-l1, i+1, 
	       j-1, j+strlen(struc)-l1-3, E*0.01);
#endif	
	if (n_subopt+1>=n_max) {
	  n_max *= 2;
	  subopt = (duplexT *) xrealloc(subopt, n_max*sizeof(duplexT));
	}
	subopt[n_subopt].i = MIN2(i+1,n1);
	subopt[n_subopt].j = MAX2(j-1,1);
	subopt[n_subopt].energy = E * 0.01;
	subopt[n_subopt++].structure = struc;
      }
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
  for (k=i; k<=i0; k++) if (!st1[k-1]) st1[k-1] = '.';
  for (k=j0; k<=j; k++) if (!st2[k-1]) st2[k-1] = '.';
  strcpy(struc, st1+i-1); strcat(struc, "&"); 
  strcat(struc, st2+j0-1);
  
  /* printf("%s %3d,%-3d : %3d,%-3d\n", struc, i,i0,j0,j);  */
  free(st1); free(st2);

  return struc;
}

/*---------------------------------------------------------------------------*/

PRIVATE void encode_seq(const char *s1, const char *s2) {
  unsigned int i,l;

  l = strlen(s1);
  S1 = (short *) space(sizeof(short)*(l+1));
  SS1= (short *) space(sizeof(short)*(l+1));
  /* SS1 exists only for the special X K and I bases and energy_set!=0 */
  
  for (i=1; i<=l; i++) { /* make numerical encoding of sequence */
    S1[i]= (short) encode_char(toupper(s1[i-1]));
    SS1[i] = alias[S1[i]];   /* for mismatches of nostandard bases */
  }

  l = strlen(s2);
  S2 = (short *) space(sizeof(short)*(l+1));
  SS2= (short *) space(sizeof(short)*(l+1));
  /* SS2 exists only for the special X K and I bases and energy_set!=0 */
  
  for (i=1; i<=l; i++) { /* make numerical encoding of sequence */
    S2[i]= (short) encode_char(toupper(s2[i-1]));
    SS2[i] = alias[S2[i]];   /* for mismatches of nostandard bases */
  }
}

/*---------------------------------------------------------------------------*/

PRIVATE void update_dfold_params(void)
{
  P = scale_parameters();
  make_pair_matrix();
}

/*---------------------------------------------------------------------------*/

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

