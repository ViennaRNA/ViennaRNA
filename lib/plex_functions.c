 /* Last changed Time-stamp: <2010-06-30 17:24:43 wolfgang> */
/*                
	   compute potentially pseudoknotted structures of an RNA sequence
			     Ivo Hofacker
			  Vienna RNA package
*/

/*
  library containing the function used in PKplex
  it generates pseudoknotted structures by letting the sequence form a duplex structure with itself
*/

//#include <config.h>
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
#include "PKplex.h"
#include <time.h>
#include "loop_energies.h"

//static char rcsid[] UNUSED = "$Id: plex.c,v 1.14 2007/06/12 12:50:16 htafer Exp $";

#define PUBLIC
#define PRIVATE static
#undef MAXLOOP
#define MAXLOOP 10

PRIVATE void  update_dfold_params(void);
PRIVATE void duplexfold_XS(const char *s1, int **access_s1, const int threshold, const int alignment_length);
PRIVATE char *backtrack_XS(int kk, int ll, const int ii, const int jj, const int alignment_length);
PRIVATE void encode(const char *s1);
PRIVATE void make_ptypes(const char *structure);

#define MIN2(A, B)      ((A) < (B) ? (A) : (B))
#define MAX2(A, B)      ((A) > (B) ? (A) : (B))

PRIVATE paramT *P = NULL;
PRIVATE int ***c3; /* energy array used in duplexfold */
PRIVATE short  *S1, *SS1;
PRIVATE int   n1;
extern  int  E_IntLoop(int n1, int n2, int type, int type_2, int si1, int sj1, int sp1, int sq1, paramT *P);
extern  int E_ExtLoop(int type, int si1, int sj1, paramT *P);
PRIVATE char  *ptype;   /* precomputed array of pair types */
PRIVATE int *indx; /* index for moving in the triangle matrices ptype[] */

extern dupVar *PlexHits;
extern int PlexHitsArrayLength;
extern int NumberOfHits;
extern int verbose;

/*-----------------------------------------------------------------------duplexfold_XS---------------------------------------------------------------------------*/

PRIVATE void duplexfold_XS(const char *s1, int **access_s1, const int threshold, const int alignment_length) {
  int i, j, k, l, p, q, Emin=INF, l_min=0, k_min=0, j_min=0;
  int type, type2, type3, E, tempK;
  char *struc;
  struc=NULL;

  c3 = (int ***) space(sizeof(int **) * (n1-20));
  for (i=0; i<(n1-20); i++) {
    c3[i] = (int **) space(sizeof(int *) * alignment_length);
    for (j=0; j<alignment_length; j++) {
      c3[i][j]=(int *) space(sizeof(int) * alignment_length);
    }
  }
  
  i=n1-9;
  while( i-- > 11 ){
    Emin=INF;
    j_min=0;
    l_min=0;
    k_min=0;
    
    //init all matrix elements to INF
    for (j=0; j<(n1-20); j++){
      for(k=0;k<alignment_length;k++){
	for(l=0;l<alignment_length;l++){
	  c3[j][k][l]=INF;
	}
      }
    }

    //matrix starting values for (i,j)-basepairs
    for(j=i+4; j<n1-10; j++) {     
      type=ptype[indx[j]+i];
      if (type) {
	c3[j-11][alignment_length-1][0] = P->DuplexInit;
	c3[j-11][alignment_length-1][0] += E_ExtLoop(type, SS1[i+1], SS1[j-1], P);
	//if (type>2) c3[j-11][alignment_length-1][0] += P->TerminalAU;
	//c3[j-11][alignment_length-1][0]+=P->dangle3[rtype[type]][SS1[i+1]];
	//c3[j-11][alignment_length-1][0]+=P->dangle5[rtype[type]][SS1[j-1]];
      }
    }
    
    int i_pos_begin=MAX2(9, i-alignment_length);
	
    //fill matrix
    for (k=i-1; k>i_pos_begin; k--) {
      tempK=alignment_length-i+k-1;
      for (l=i+5; l<n1-9; l++) {
	type2=ptype[indx[l]+k];
	if (!type2) continue;
	for (p=k+1; (p<=i) && (p<=k+MAXLOOP+1); p++) {
	  for (q = l-1; (q>=i+4) && (q>=l-MAXLOOP-1); q--) {
	    if (p-k+l-q-2>MAXLOOP) break;
	    type3=ptype[indx[q]+p];
	    if(!type3) continue;
	    E = E_IntLoop(p-k-1, l-q-1, type2, rtype[type3],SS1[k+1], SS1[l-1], SS1[p-1], SS1[q+1], P);
	    for (j=MAX2(i+4, l-alignment_length+1); j<=q; j++) { 
	      type=ptype[indx[j]+i];
	      if (type) {
		c3[j-11][tempK][l-j] = MIN2(c3[j-11][tempK][l-j], c3[j-11][alignment_length-i+p-1][q-j]+E);
	      }
	    }//next j
	  }//next q
	}//next p
      }//next l
    }//next k

    //read out matrix minimum
    for(j=i+4; j<n1-10; j++) { 
      type=ptype[indx[j]+i];
      if (!type) continue;
      int j_pos_end=MIN2(n1-9,j+alignment_length);
      for (k=i-1; k>i_pos_begin; k--) {
	for (l=j+1; l<j_pos_end; l++) {
	  type2=ptype[indx[l]+k];	  
	  if (!type2) continue;
	  E = c3[j-11][alignment_length-i+k-1][l-j];
	  //if (type2>2) E += P->TerminalAU;
	  E+=access_s1[i-k+1][i]+access_s1[l-j+1][l];
	  E+=E_ExtLoop(type2,((k>i_pos_begin+1)? SS1[k-1]:-1),((l<j_pos_end-1)? SS1[l+1]:-1),P);
	  //if (k>i_pos_begin+1) E += P->dangle5[type2][SS1[k-1]];
	  //if (l<j_pos_end-1) E += P->dangle3[type2][SS1[l+1]];
	  if (E<Emin) {
	    Emin=E; k_min=k; l_min=l;
	    j_min=j;
	  }
	}
      }
    }
  	
    if(Emin  < threshold){ 
      struc = backtrack_XS(k_min, l_min, i, j_min, alignment_length);
      
      //lets take care of the dangles
      //find best combination
      int dx_5, dx_3, dy_5, dy_3,dGx,dGy,bonus_x, bonus_y;
      dx_5=0; dx_3=0; dy_5=0; dy_3=0;dGx=0;dGy=0;bonus_x=0; bonus_y;
      dGx = access_s1[i-k_min+1][i];dx_3=0; dx_5=0;bonus_x=0;
      dGy = access_s1[l_min-j_min+1][l_min];
      PlexHits[NumberOfHits].tb=k_min -10 -dx_5;
      PlexHits[NumberOfHits].te=i -10 + dx_3;
      PlexHits[NumberOfHits].qb=j_min -10 - dy_5;
      PlexHits[NumberOfHits].qe=l_min -10 + dy_3;
      PlexHits[NumberOfHits].ddG=(double) Emin * 0.01;
      PlexHits[NumberOfHits].dG1=(double) dGx*0.01 ;
      PlexHits[NumberOfHits].dG2=(double) dGy*0.01 ;
      PlexHits[NumberOfHits].energy= PlexHits[NumberOfHits].ddG - PlexHits[NumberOfHits].dG1 - PlexHits[NumberOfHits].dG2;
      PlexHits[NumberOfHits].structure = struc;
	
      //output:
      if(PlexHits[NumberOfHits].energy * 100 < threshold){
	if (verbose) printf("%s %3d,%-3d : %3d,%-3d (%5.2f = %5.2f + %5.2f + %5.2f)\n", PlexHits[NumberOfHits].structure, PlexHits[NumberOfHits].tb, PlexHits[NumberOfHits].te, PlexHits[NumberOfHits].qb, PlexHits[NumberOfHits].qe, PlexHits[NumberOfHits].ddG, PlexHits[NumberOfHits].energy, PlexHits[NumberOfHits].dG1, PlexHits[NumberOfHits].dG2);
	NumberOfHits++;
	if(NumberOfHits==PlexHitsArrayLength-1) {
	  printf("Array PlexHits is full\n");
	  PlexHitsArrayLength*=2;
	  PlexHits = (dupVar *) xrealloc(PlexHits,sizeof(dupVar)*PlexHitsArrayLength);
	}
      }
    }
  }

  for (i=0; i<(n1-20); i++) {
    for (j=0; j<alignment_length; j++) {
      free(c3[i][j]);
    }
    free(c3[i]);
  }
  free(c3);
}


PRIVATE char *backtrack_XS(int k, int l, const int i, const int j, const int alignment_length) {
  /* backtrack structure going backwards from i, and forwards from j 
     return structure in bracket notation with & as separator */
  int p, q, type, type2, E, traced, i0, j0;
  char *st1, *st2, *struc;
  st1 = (char *) space(sizeof(char)*(i-k+2));
  st1[i-k+1]='\0';
  st2 = (char *) space(sizeof(char)*(l-j+2));
  st2[l-j+1]='\0';
 
  i0=k; j0=l;
  while (k<=i && l>=j) {
    E = c3[j-11][alignment_length-i+k-1][l-j]; traced=0;
    st1[k-i0] = '(';
    st2[l-j] = ')';
    
    type=ptype[indx[l]+k];
    if (!type) nrerror("backtrack failed in fold duplex bli");
    for (p=k+1; p<=i; p++) {
      for (q=l-1; q>=j; q--) {	
	int LE;
 	if (p-k+l-q-2>MAXLOOP) break;
	type2=ptype[indx[q]+p];
	if (!type2) continue;
 	LE = E_IntLoop(p-k-1, l-q-1, type, rtype[type2], SS1[k+1], SS1[l-1], SS1[p-1], SS1[q+1], P);
 	if (E == c3[j-11][alignment_length-i+p-1][q-j]+LE) {
	  traced=1; 
 	  k=p; l=q;
	  break;
	}
      }
      if (traced) break;
    }
    if (!traced) {
      E-=E_ExtLoop(type2, ((k<i)?SS1[k+1]:-1), ((l>j-1)? SS1[l-1]:-1), P); 
      //if (k<i) E -= P->dangle3[rtype[type]][SS1[k+1]];
      //if (l>j-1)  E -= P->dangle5[rtype[type]][SS1[l-1]];
      //if (type>2) E -= P->TerminalAU;
      break;
      if (E != P->DuplexInit) {
        nrerror("backtrack failed in fold duplex bal");
      } else break;
    }
  }
  struc = (char *) space(k-i0+1+j0-l+1+2);
 
  for (p=0; p<=i-i0; p++){   
    if (!st1[p]) st1[p] = '.';
  }

  for (p=0; p<=j0-j; p++) {    
    if (!st2[p]) {
      st2[p] = '.';
    }
  }
    
  strcpy(struc, st1);
  strcat(struc, "&");
  strcat(struc, st2);
  free(st1); free(st2); 
  return struc;
}

dupVar** PKLduplexfold_XS(const char *s1, int **access_s1, const int threshold, const int alignment_length, const int delta)
{
  if ((!P) || (fabs(P->temperature - temperature)>1e-6))
    update_dfold_params();

  n1 = (int) strlen(s1);
  encode(s1);

  indx = (int *) space(sizeof(int)*(n1+1));
  ptype = (char *) space(sizeof(char)*((n1*(n1+1))/2+2));
  int n;
  for (n = 1; n <= n1; n++)
    indx[n] = (n*(n-1)) >> 1;        /* n(n-1)/2 */
  make_ptypes(s1);

  P->DuplexInit=0;
  duplexfold_XS(s1,access_s1,threshold, alignment_length);
  free(S1); free(SS1);
  free(indx); free(ptype);
  return NULL;
}

/*---------------------------------UTILS------------------------------------------*/

PRIVATE void update_dfold_params(void)
{
  P = scale_parameters();
  make_pair_matrix();
}

/*---------------------------------------------------------------------------*/

PRIVATE void encode(const char *s1) {
  unsigned int i,l;

  l = strlen(s1);
  S1 = (short *) space(sizeof(short)*(l+2));
  S1[0] = (short) l;

  /* make numerical encoding of sequence */
  for (i=1; i<=l; i++)
    S1[i]= (short) encode_char(toupper(s1[i-1]));

  /* for circular folding add first base at position n+1 */
  S1[l+1] = S1[1];

  SS1= (short *) space(sizeof(short)*(l+1));
  /* SS1 exists only for the special X K and I bases and energy_set!=0 */
  
  for (i=1; i<=l; i++) { /* make numerical encoding of sequence */
    SS1[i] = alias[S1[i]];   /* for mismatches of nostandard bases */
  }
}


PRIVATE void make_ptypes(const char *structure) {
  int n,i,j,k,l;

  n=S1[0];
  for (k=1; k<n-TURN; k++)
    for (l=1; l<=2; l++) {
      int type,ntype=0,otype=0;
      i=k; j = i+TURN+l; if (j>n) continue;
      type = pair[S1[i]][S1[j]];
      while ((i>=1)&&(j<=n)) {
	if ((i>1)&&(j<n)) ntype = pair[S1[i-1]][S1[j+1]];
	if (noLonelyPairs && (!otype) && (!ntype))
	  type = 0; /* i.j can only form isolated pairs */
	ptype[indx[j]+i] = (char) type;
	otype =  type;
	type  = ntype;
	i--; j++;
      }
    }
}
