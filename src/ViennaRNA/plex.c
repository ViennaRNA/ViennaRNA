/*
           compute the duplex structure of two RNA strands,
                allowing only inter-strand base pairs.
         see cofold() for computing hybrid structures without
                             restriction.
                             Ivo Hofacker
                          Vienna RNA package

*/


/*
  library containing the function used in rnaplex
  the program rnaplex uses the following function
  Lduplexfold: finds high scoring segments
  it stores the end-position of these segments in an array
  and call then for each of these positions the duplexfold function
  which allows one to make backtracking for each of the high scoring position
  It allows one to find suboptimal partially overlapping (depends on a a parameter)
  duplexes between a long RNA and a shorter one.
  Contrarly to RNAduplex, the energy model is not in E~log(N),
  where N is the length of an interial loop but used an affine model,
  where the extension and begin parameter are fitted to the energy
  parameter used by RNAduplex. This allows one to check for duplex between a short RNA(20nt)
  and a long one at the speed of 1Mnt/s. At this speed the whole genome (3Gnt) can be analyzed for one siRNA
  in about 50 minutes.
  The algorithm is based on an idea by Durbin and Eddy:when the alginment reach a value larger than a
  given threshold this value is stored in an array. When the alignment score goes
  then under this threshold, the alignemnent begin from this value, in that way the backtracking allow us
  to find all non-overlapping high-scoring segments.
  For more information check "durbin, biological sequence analysis"
*/

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <string.h>
#include "ViennaRNA/utils.h"
#include "ViennaRNA/energy_par.h"
#include "ViennaRNA/fold_vars.h"
#include "ViennaRNA/fold.h"
#include "ViennaRNA/pair_mat.h"
#include "ViennaRNA/params.h"
#include "ViennaRNA/plex.h"
#include "ViennaRNA/ali_plex.h"
#include "ViennaRNA/loop_energies.h"
/* #################SIMD############### */

/* int subopt_sorted=0; */

#define PUBLIC
#define PRIVATE static

#define STACK_BULGE1  1   /* stacking energies for bulges of size 1 */
#define NEW_NINIO     1   /* new asymetry penalty */
#define ARRAY 32          /*array size*/
#define UNIT 100
#define MINPSCORE -2 * UNIT

/**
*** Macro that define indices for the Single Array approach defined in FLduplexfold_XS->gain of 20% in runtime
*** so that everything is done in a 1D array.
*** input is idx for i, j for j and the length of the query RNA
*** 1D is divided in 6 subarrays, one for each number of allowed state
*** The length of each subarray is 5*L. 5 the maximal stored distance on the target sequence,
*** L is the length of the query sequence
**/
#define LCI(i,j,l)      ((i     )*l + j)
#define LINI(i,j,l)     ((i +  5)*l + j)
#define LBXI(i,j,l)     ((i + 10)*l + j)
#define LBYI(i,j,l)     ((i + 15)*l + j)
#define LINIX(i,j,l)    ((i + 20)*l + j)
#define LINIY(i,j,l)    ((i + 25)*l + j)

PRIVATE void  encode_seqs(const char *s1, const char *s2);
PRIVATE short *encode_seq(const char *seq);
PRIVATE void  update_dfold_params(void);
/**
*** duplexfold(_XS)/backtrack(_XS) computes duplex interaction with standard energy and considers extension_cost
*** find_max(_XS)/plot_max(_XS) find suboptimals and MFE
*** fduplexfold(_XS) computes duplex in a plex way
**/
PRIVATE duplexT duplexfold(const char *s1, const char *s2, const int extension_cost);
PRIVATE char * backtrack(int i, int j, const int extension_cost);
PRIVATE void   find_max(const int *position, const int *position_j, const int delta, const int threshold, const int length, const char *s1, const char *s2, const int extension_cost, const int fast,const int il_a, const int il_b, const int b_a, const int b_b);
PRIVATE void   plot_max(const int max, const int max_pos, const int max_pos_j, const int alignment_length, const char *s1, const char *s2, const int extension_cost, const int fast,const int il_a, const int il_b, const int b_a, const int b_b);

/* PRIVATE duplexT duplexfold_XS(const char *s1, const char *s2,const int **access_s1, const int **access_s2, const int i_pos, const int j_pos, const int threshold); */
PRIVATE duplexT duplexfold_XS(const char *s1, const char *s2,const int **access_s1, const int **access_s2, const int i_pos, const int j_pos, const int threshold, const int i_flag, const int j_flag);
/* PRIVATE char *   backtrack_XS(int i, int j, const int** access_s1, const int** access_s2); */
PRIVATE char *backtrack_XS(int i, int j, const int **access_s1,const int **access_s2,const int i_flag, const int j_flag );
PRIVATE void find_max_XS(const int *position, const int *position_j,const int delta, const int threshold, const int alignment_length,
                         const char *s1, const char *s2, const int **access_s1, const int **access_s2, const int fast,const int il_a, const int il_b, const int b_a, const int b_b);
PRIVATE void plot_max_XS(const int max, const int max_pos, const int max_pos_j, const int alignment_length, const char *s1, const char *s2, const int **access_s1, const int **access_s2, const int fast,const int il_a, const int il_b, const int b_a, const int b_b);
PRIVATE duplexT fduplexfold(const char *s1, const char *s2, const int extension_cost, const int il_a, const int il_b, const int b_a, const int b_b);
PRIVATE char *fbacktrack(int i, int j, const int extension_cost,const int il_a, const int il_b, const int b_a, const int b_b, int *dG);
PRIVATE duplexT fduplexfold_XS(const char *s1, const char *s2, const int **access_s1, const int **access_s2, const int i_pos, const int j_pos, const int threshold,const int il_a, const int il_b, const int b_a, const int b_b);
PRIVATE char * fbacktrack_XS(int i, int j, const int **access_s1, const int **access_s2, const int i_pos, const int j_pos, const int il_a, const int il_b, const int b_a, const int b_b, int *dGe, int *dGeplex, int *dGx, int *dGy);



/*@unused@*/

#define MAXSECTORS      500     /* dimension for a backtrack array */
#define LOCALITY        0.      /* locality parameter for base-pairs */

PRIVATE vrna_param_t *P = NULL;

/**
*** energy array used in fduplexfold and fduplexfold_XS
*** We do not use the 1D array here as it is not time critical
*** It also makes the code more readable
*** c -> stack;in -> interior loop;bx/by->bulge;inx/iny->1xn loops
**/

PRIVATE int   **c=NULL, **in=NULL, **bx=NULL, **by=NULL, **inx=NULL, **iny=NULL;

/**
*** S1, SS1, ... contains the encoded sequence for target and query
*** n1, n2, n3, n4 contains target and query length
**/

PRIVATE short  *S1=NULL, *SS1=NULL, *S2=NULL, *SS2=NULL;/*contains the sequences*/
PRIVATE int   n1,n2;    /* sequence lengths */
PRIVATE int n3, n4; /*sequence length for the duplex*/;


/*-----------------------------------------------------------------------duplexfold_XS---------------------------------------------------------------------------*/

/**
*** duplexfold_XS is the pendant to the duplex function as defined in duplex.c
*** but takes the accessibility into account. It is similar to the MFE version of RNAup
*** The only approximation made is that target 3' end - query 5' end base pair is known
*** s1,s2 are the query and target sequence; access_s1, access_s2 are the accessibility
*** profiles, i_pos, j_pos are the coordinates of the closing pair.
**/



PRIVATE duplexT duplexfold_XS(const char *s1, const char *s2, const int **access_s1, const int **access_s2, const int i_pos, const int j_pos, const int threshold, const int i_flag, const int j_flag) {
  int i, j,p,q, Emin=INF, l_min=0, k_min=0;
  char *struc;
  vrna_md_t   md;

  struc=NULL;
  duplexT mfe;
  n3 = (int) strlen(s1);
  n4 = (int) strlen(s2);

  set_model_details(&md);

  if ((!P) || (fabs(P->temperature - temperature)>1e-6)) {
    update_fold_params();
    if(P)
      free(P);
    P = vrna_params(&md);
    make_pair_matrix();
  }

  c = (int **) vrna_alloc(sizeof(int *) * (n3+1));
  for (i=0; i<=n3; i++) c[i] = (int *) vrna_alloc(sizeof(int) * (n4+1));
  for (i=0; i<=n3; i++){
    for(j=0;j<=n4;j++){
      c[i][j]=INF;
    }
  }
  encode_seqs(s1, s2);
  int type, type2, type3, E, k,l;
  i=n3-i_flag; j=1+j_flag;
  type = pair[S1[i]][S2[j]];
  if(!type){
    printf("Error during initialization of the duplex in duplexfold_XS\n");
    mfe.structure=NULL;
    mfe.energy = INF;
    return mfe;
  }
  c[i][j] = P->DuplexInit;
  /**  if (type>2) c[i][j] += P->TerminalAU;
  ***  c[i][j]+=P->dangle3[rtype[type]][SS1[i+1]];
  ***  c[i][j]+=P->dangle5[rtype[type]][SS2[j-1]];
  *** The three above lines are replaced by the line below
  **/


  c[i][j] += E_ExtLoop(rtype[type], (j_flag ? SS2[j-1] : -1) , (i_flag ? SS1[i+1] : -1),  P);

/*   if(j_flag ==0 && i_flag==0){ */
/*     c[i][j] += E_ExtLoop(rtype[type], -1 , -1 , P); */
/*   }else if(j_flag ==0 && i_flag==1){ */
/*     c[i][j] += E_ExtLoop(rtype[type], -1 , SS1[i+1], P); */
/*   }else if(j_flag ==1 && i_flag==0){ */
/*     c[i][j] += E_ExtLoop(rtype[type], SS2[j-1] , -1, P); */
/*   }else { */
/*     c[i][j] += E_ExtLoop(rtype[type], SS2[j-1] , SS1[i+1], P); */
/*   } */
  /*  Just in case we have only one bp, we initialize ... */
  /*  k_min, l_min and Emin */
  k_min=i; l_min=j;Emin=c[i][j];
  for (k=i; k>1 ; k--) {
    if(k<i) c[k+1][0]=INF;
    for (l=j; l<=n4-1; l++) {
      if(!(k==i && l==j)){
        c[k][l]=INF;
      }
      type2 = pair[S1[k]][S2[l]];
      if (!type2) continue;
      for (p=k+1; p<= n3 - i_flag && p<k+MAXLOOP-1; p++) {
        for (q = l-1; q >= 1+j_flag; q--) {
          if (p-k+l-q-2>MAXLOOP) break;
          type3=pair[S1[p]][S2[q]];
          if(!type3) continue;
          E = E_IntLoop(p-k-1, l-q-1, type2, rtype[type3],SS1[k+1], SS2[l-1], SS1[p-1], SS2[q+1],P);
          c[k][l] = MIN2(c[k][l], c[p][q]+E);
        }
      }
      E = c[k][l];
      E+=access_s1[i-k+1][i_pos]+access_s2[l-1][j_pos+(l-1)-1];
      /**if (type2>2) E += P->TerminalAU;
      ***if (k>1) E += P->dangle5[type2][SS1[k-1]];
      ***if (l<n4) E += P->dangle3[type2][SS2[l+1]];
      *** Replaced by the line below
      **/
      E+=E_ExtLoop(type2, (k>1) ? SS1[k-1] : -1, (l<n4) ? SS2[l+1] : -1, P);

      if (E<Emin) {
        Emin=E; k_min=k; l_min=l;
      }
    }
  }

  if(Emin  > threshold){
    mfe.energy=INF;
    mfe.ddG=INF;
    mfe.structure=NULL;
    for (i=0; i<=n3; i++) free(c[i]);
    free(c);
    free(S1); free(S2); free(SS1); free(SS2);
    return mfe;
  } else{
    struc = backtrack_XS(k_min, l_min, access_s1, access_s2, i_flag, j_flag);
  }


  /**
  *** find best dangles combination
  **/
  int dx_5, dx_3, dy_5, dy_3,dGx,dGy,bonus_x;
  dx_5=0; dx_3=0; dy_5=0; dy_3=0;dGx=0;dGy=0;bonus_x=0;
  /* x--------x */
  /*  |||||||| */
  /* x--------x */
   dGx = access_s1[i-k_min+1][i_pos];dx_3=0; dx_5=0;bonus_x=0;
   dGy = access_s2[l_min-j+1][j_pos + (l_min-1)];

  mfe.tb=i_pos -9 - i + k_min -1 -dx_5;
  mfe.te=i_pos -9 -1 + dx_3;
  mfe.qb=j_pos -9 -1 - dy_5;
  mfe.qe=j_pos + l_min -3 -9 + dy_3;
  mfe.ddG=(double) Emin * 0.01;
  mfe.dG1=(double) dGx*0.01 ;
  mfe.dG2=(double) dGy*0.01 ;

  mfe.energy= mfe.ddG - mfe.dG1 - mfe.dG2;

  mfe.structure = struc;
  for (i=0; i<=n3; i++) free(c[i]);
  free(c);
  free(S1); free(S2); free(SS1); free(SS2);
  return mfe;
}






PRIVATE char *backtrack_XS(int i, int j, const int **access_s1,const int **access_s2, const int i_flag, const int j_flag) {
  /* backtrack structure going backwards from i, and forwards from j
     return structure in bracket notation with & as separator */
  int k, l, type, type2, E, traced, i0, j0;
  char *st1, *st2, *struc;
  st1 = (char *) vrna_alloc(sizeof(char)*(n3+1));
  st2 = (char *) vrna_alloc(sizeof(char)*(n4+1));
  i0=i;/*MAX2(i-1,1);*/j0=j;/*MIN2(j+1,n4);*/
  while (i<=n3-i_flag && j>=1+j_flag) {
    E = c[i][j]; traced=0;
    st1[i-1] = '(';
    st2[j-1] = ')';
    type = pair[S1[i]][S2[j]];
    if (!type) vrna_message_error("backtrack failed in fold duplex bli");
    for (k=i+1; k<=n3 && k>i-MAXLOOP-2; k++) {
      for (l=j-1; l>=1; l--) {
        int LE;
        if (i-k+l-j-2>MAXLOOP) break;
        type2 = pair[S1[k]][S2[l]];
        if (!type2) continue;
        LE = E_IntLoop(k-i-1, j-l-1, type, rtype[type2], SS1[i+1], SS2[j-1], SS1[k-1], SS2[l+1],P);
        if (E == c[k][l]+LE) {
          traced=1;
          i=k; j=l;
          break;
        }
      }
      if (traced) break;
    }
    if (!traced) {
#if 0
      if (i<n3) E -= P->dangle3[rtype[type]][SS1[i+1]];/* +access_s1[1][i+1]; */
      if (j>1)  E -= P->dangle5[rtype[type]][SS2[j-1]];/* +access_s2[1][j+1]; */
      if (type>2) E -= P->TerminalAU;
#endif
      E -= E_ExtLoop(rtype[type], SS2[j-1] , SS1[i+1], P);
      break;
      if (E != P->DuplexInit) {
        vrna_message_error("backtrack failed in fold duplex bal");
      } else break;
    }
  }
  /* if (i<n3)  i++; */
  /* if (j>1)   j--; */
  struc = (char *) vrna_alloc(i-i0+1+j0-j+1+2);
  for (k=MAX2(i0,1); k<=i; k++) if (!st1[k-1]) st1[k-1] = '.';
  for (k=j; k<=j0; k++) if (!st2[k-1]) st2[k-1] = '.';
  strcpy(struc, st1+MAX2(i0-1,0)); strcat(struc, "&");
  strcat(struc, st2+j-1);
  free(st1); free(st2);
  return struc;
}


/**
*** fduplexfold(_XS) computes the interaction based on the plex energy model.
*** Faster than duplex approach, but not standard model compliant
*** We use the standard matrix (c, in, etc..., because we backtrack)
**/

PRIVATE duplexT fduplexfold_XS(const char *s1, const char *s2, const int **access_s1, const int **access_s2, const int i_pos, const int j_pos, const int threshold,const int il_a, const int il_b, const int b_a, const int b_b) {
  /**
  *** i,j recursion index
  *** Emin, i_min, j_min MFE position and energy
  *** mfe struc duplex structure
  **/
  int i, j, Emin, i_min, j_min,l1;
  duplexT mfe;
  char *struc;
  /**
  *** bext=b_a bulge extension parameter for linear model
  *** iopen=il_b interior opening for linear model
  *** iext_s=2*il_a asymmetric extension for interior loop
  *** iext_ass=60+il_a symmetric extension for interior loop
  *** min_colonne=INF; max score of a row
  *** i_length;
  *** max_pos; position of best hit during recursion on target
  *** max_pos_j; position of best hit during recursion on query
  *** temp; temp variable for min_colonne
  *** min_j_colonne; position of the minimum on query in row j
  *** max=INF; absolute MFE
  *** n3,n4 length of target and query
  *** DJ contains the accessibility penalty for the query sequence
  *** maxPenalty contains the maximum penalty
  **/
  int bopen=b_b;
  int bext=b_a;
  int iopen=il_b;
  int iext_s=2*il_a;/* iext_s 2 nt nucleotide extension of interior loop, on i and j side */
  int iext_ass=50+il_a;/* iext_ass assymetric extension of interior loop, either on i or on j side. */
  int min_colonne=INF; /* enthaelt das maximum einer kolonne */
  int i_length;
  int max_pos;/* get position of the best hit */
  int max_pos_j;
  int temp=INF;
  int min_j_colonne;
  int max=INF;
  int **DJ;
  int maxPenalty[4];
  vrna_md_t   md;

  /**
  *** variable initialization
  **/
  n3 = (int) strlen(s1);
  n4 = (int) strlen(s2);

  set_model_details(&md);
  if ((!P) || (fabs(P->temperature - temperature)>1e-6)) {
    update_fold_params();
    if(P)
      free(P);
    P = vrna_params(&md);
    make_pair_matrix();
  }
  /**
  *** array initialization
  **/
  c  = (int**) vrna_alloc(sizeof(int *) * (n3+1));
  in = (int**) vrna_alloc(sizeof(int *) * (n3+1));
  bx = (int**) vrna_alloc(sizeof(int *) * (n3+1));
  by = (int**) vrna_alloc(sizeof(int *) * (n3+1));
  inx= (int**) vrna_alloc(sizeof(int *) * (n3+1));
  iny= (int**) vrna_alloc(sizeof(int *) * (n3+1));
  /* #pragma omp parallel for */
  for (i=0; i<=n3; i++){
    c[i]  = (int *) vrna_alloc(sizeof(int) * (n4+1));
    in[i] = (int *) vrna_alloc(sizeof(int) * (n4+1));
    bx[i] = (int *) vrna_alloc(sizeof(int) * (n4+1));
    by[i] = (int *) vrna_alloc(sizeof(int) * (n4+1));
    inx[i]= (int *) vrna_alloc(sizeof(int) * (n4+1));
    iny[i]= (int *) vrna_alloc(sizeof(int) * (n4+1));
  }
  for(i=0; i<n3; i++){
    for(j=0; j<n4; j++){
      in[i][j]=INF;/* no in before  1 */
      c[i][j] =INF; /* no bulge and no in before n2 */
      bx[i][j]=INF;/* no bulge before 1 */
      by[i][j]=INF;
      inx[i][j]=INF;/* no bulge before 1 */
      iny[i][j]=INF;
    }
  }
  /**
  *** sequence encoding
  **/
  encode_seqs(s1,s2);
  /**
  *** Compute max accessibility penalty for the query only once
  **/
  maxPenalty[0]=(int) -1*P->stack[2][2]/2;
  maxPenalty[1]=(int) -1*P->stack[2][2];
  maxPenalty[2]=(int) -3*P->stack[2][2]/2;
  maxPenalty[3]=(int) -2*P->stack[2][2];


  DJ=(int **)   vrna_alloc(4*sizeof(int*));
  DJ[0]=(int *) vrna_alloc((1+n4)*sizeof(int));
  DJ[1]=(int *) vrna_alloc((1+n4)*sizeof(int));
  DJ[2]=(int *) vrna_alloc((1+n4)*sizeof(int));
  DJ[3]=(int *) vrna_alloc((1+n4)*sizeof(int));

  j=n4-9;
  while(--j>9){
    int jdiff=j_pos+j-11;
    /**
    *** Depending in which direction (i:1->n vs j:m->1) the accessibility is computed we get slightly different results.
    *** We reduce the discrepancies by taking the average of d^i_k and d^j_l
    **/
    DJ[0][j] = 0.5*(access_s2[5][jdiff+4] - access_s2[4][jdiff+4] + access_s2[5][jdiff]  -access_s2[4][jdiff-1]  );
    DJ[1][j] = 0.5*(access_s2[5][jdiff+5] - access_s2[4][jdiff+5] + access_s2[5][jdiff+1]-access_s2[4][jdiff]  ) + DJ[0][j];
    DJ[2][j] = 0.5*(access_s2[5][jdiff+6] - access_s2[4][jdiff+6] + access_s2[5][jdiff+2]-access_s2[4][jdiff+1]) + DJ[1][j];
    DJ[3][j] = 0.5*(access_s2[5][jdiff+7] - access_s2[4][jdiff+7] + access_s2[5][jdiff+3]-access_s2[4][jdiff+2]) + DJ[2][j];



/*
    DJ[0][j] = access_s2[5][jdiff+4] - access_s2[4][jdiff+4]           ;
    DJ[1][j] = access_s2[5][jdiff+5] - access_s2[4][jdiff+5] + DJ[0][j];
    DJ[2][j] = access_s2[5][jdiff+6] - access_s2[4][jdiff+6] + DJ[1][j];
    DJ[3][j] = access_s2[5][jdiff+7] - access_s2[4][jdiff+7] + DJ[2][j];
    DJ[0][j] = MIN2(DJ[0][j],maxPenalty[0]);
    DJ[1][j] = MIN2(DJ[1][j],maxPenalty[1]);
    DJ[2][j] = MIN2(DJ[2][j],maxPenalty[2]);
    DJ[3][j] = MIN2(DJ[3][j],maxPenalty[3]);
*/
  }

  /**
  *** Start of the recursion
  *** first and last 10 nucleotides on target and query are dummy nucleotides
  *** allow to reduce number of if test
  **/
  i=11;
  i_length=n3-9;
  while(i < i_length) {
    int di1,di2,di3,di4;
    int idiff=i_pos-(n3-10-i);
    di1 = 0.5*(access_s1[5][idiff+4] - access_s1[4][idiff+4] + access_s1[5][idiff]   - access_s1[4][idiff-1]);
    di2 = 0.5*(access_s1[5][idiff+3] - access_s1[4][idiff+3] + access_s1[5][idiff-1] - access_s1[4][idiff-2]) + di1;
    di3 = 0.5*(access_s1[5][idiff+2] - access_s1[4][idiff+2] + access_s1[5][idiff-2] - access_s1[4][idiff-3]) + di2;
    di4 = 0.5*(access_s1[5][idiff+1] - access_s1[4][idiff+1] + access_s1[5][idiff-3] - access_s1[4][idiff-4]) + di3;
/*
    di1 =  access_s1[5][idiff]   - access_s1[4][idiff-1];
    di2 =  access_s1[5][idiff-1] - access_s1[4][idiff-2] + di1;
    di3 =  access_s1[5][idiff-2] - access_s1[4][idiff-3] + di2;
    di4 =  access_s1[5][idiff-3] - access_s1[4][idiff-4] + di3;
    di1=MIN2(di1,maxPenalty[0]);
    di2=MIN2(di2,maxPenalty[1]);
    di3=MIN2(di3,maxPenalty[2]);
    di4=MIN2(di4,maxPenalty[3]);
*/
    j=n4-9;
    min_colonne=INF;
    while (10 < --j) {
      int dj1,dj2,dj3,dj4;
      int jdiff=j_pos+j-11;
      dj1 = DJ[0][j];
      dj2 = DJ[1][j];
      dj3 = DJ[2][j];
      dj4 = DJ[3][j];
      int type, type2;
      type = pair[S1[i]][S2[j]];
      /**
      *** Start duplex
      **/
      /*
      c[i][j]=type ? P->DuplexInit + access_s1[1][idiff]+access_s2[1][jdiff] : INF;
      */
      c[i][j]=type ? P->DuplexInit: INF;
      /**
      *** update lin bx by linx liny matrix
      **/
      type2=pair[S2[j+1]][S1[i-1]];
      /**
      *** start/extend interior loop
      **/
      in[i][j]=MIN2(c[i - 1][j+1]+P->mismatchI[type2][SS2[j]][SS1[i]]+iopen+iext_s+di1+dj1,
                    in[i - 1][j]+iext_ass + di1);

      /**
      *** start/extend nx1 target
      *** use same type2 as for in
      **/
      inx[i][j]=MIN2(c[i-1][j+1]+P->mismatch1nI[type2][SS2[j]][SS1[i]]+iopen+iext_s+di1+dj1,
                     inx[i-1][j]+iext_ass+di1);
      /**
      *** start/extend 1xn target
      *** use same type2 as for in
      **/
      iny[i][j]=MIN2(c[i-1][j+1]+P->mismatch1nI[type2][SS2[j]][SS1[i]]+iopen+iext_s+di1+dj1,
                     iny[i][j+1]+iext_ass+dj1);
      /**
      *** extend interior loop
      **/
      in[i][j]=MIN2(in[i][j],in[i][j+1]+iext_ass + dj1);
      in[i][j]=MIN2(in[i][j],in[i - 1][j+1]+iext_s + di1 + dj1);
      /**
      *** start/extend bulge target
      **/
      type2=pair[S2[j]][S1[i-1]];
      bx[i][j]=MIN2(bx[i - 1][j]+bext + di1, c[i - 1][j]+bopen+bext+(type2>2?P->TerminalAU:0) + di1);
      /**
      *** start/extend bulge query
      **/
      type2=pair[S2[j+1]][S1[i]];
      by[i][j]=MIN2(by[i][j+1]+bext+dj1, c[i][j+1]+bopen+bext+(type2>2?P->TerminalAU:0)+dj1);
      /**
      ***end update recursion
      ***######################## Start stack extension##############################
      **/
      if(!type){continue;}
      c[i][j]+=E_ExtLoop(type, SS1[i-1], SS2[j+1],P);
      /**
      *** stack extension
      **/
      if((type2=pair[S1[i-1]][S2[j+1]]))
        c[i][j]=MIN2(c[i - 1][j+1]+P->stack[rtype[type]][type2]+di1+dj1, c[i][j]);
      /**
      *** 1x0 / 0x1 stack extension
      **/
      if((type2=pair[S1[i-1]][S2[j+2]]))
        c[i][j]=MIN2(c[i - 1][j+2]+P->bulge[1]+P->stack[rtype[type]][type2]+di1+dj2,c[i][j]);
      if((type2=pair[S1[i-2]][S2[j+1]]))
        c[i][j]=MIN2(c[i - 2][j+1]+P->bulge[1]+P->stack[type2][rtype[type]]+di2+dj1,c[i][j]);
      /**
      *** 1x1 / 2x2 stack extension
      **/
      if((type2=pair[S1[i-2]][S2[j+2]]))
        c[i][j]=MIN2(c[i - 2][j+2]+P->int11[type2][rtype[type]][SS1[i-1]][SS2[j+1]]+di2+dj2, c[i][j]);
      if((type2 = pair[S1[i-3]][S2[j+3]]))
        c[i][j]=MIN2(c[i - 3][j+3]+P->int22[type2][rtype[type]][SS1[i-2]][SS1[i-1]][SS2[j+1]][SS2[j+2]]+di3+dj3,c[i][j]);
      /**
      *** 1x2 / 2x1 stack extension
      *** E_IntLoop(1,2,type2, rtype[type],SS1[i-1], SS2[j+2], SS1[i-1], SS2[j+1], P) corresponds to
      *** P->int21[rtype[type]][type2][SS2[j+2]][SS1[i-1]][SS1[i-1]]
      **/
      if((type2 = pair[S1[i-3]][S2[j+2]]))
        c[i][j]=MIN2(c[i - 3][j+2]+P->int21[rtype[type]][type2][SS2[j+1]][SS1[i-2]][SS1[i-1]]+di3+dj2, c[i][j]);
      if((type2 = pair[S1[i-2]][S2[j+3]]))
        c[i][j]=MIN2(c[i - 2][j+3]+P->int21[type2][rtype[type]][SS1[i-1]][SS2[j+1]][SS2[j+2]]+di2+dj3, c[i][j]);

      /**
      *** 2x3 / 3x2 stack extension
      **/
      if((type2 = pair[S1[i-4]][S2[j+3]]))
        c[i][j]=MIN2(c[i - 4][j+3]+P->internal_loop[5]+P->ninio[2]+
                     P->mismatch23I[type2][SS1[i-3]][SS2[j+2]]+P->mismatch23I[rtype[type]][SS2[j+1]][SS1[i-1]]+di4+dj3, c[i][j]);
      if((type2 = pair[S1[i-3]][S2[j+4]]))
        c[i][j]=MIN2(c[i - 3][j+4]+P->internal_loop[5]+P->ninio[2]+
                     P->mismatch23I[type2][SS1[i-2]][SS2[j+3]]+P->mismatch23I[rtype[type]][SS2[j+1]][SS1[i-1]]+di3+dj4, c[i][j]);

      /**
      *** So now we have to handle 1x3, 3x1, 3x3, and mxn m,n > 3
      **/
      /**
      *** 3x3 or more
      **/
      c[i][j]=MIN2(in[i - 3][j+3]+P->mismatchI[rtype[type]][SS1[i-1]][SS2[j+1]]+2*iext_s+di3+dj3, c[i][j]);
      /**
      *** 2xn or more
      **/
      c[i][j]=MIN2(in[i - 4][j+2]+P->mismatchI[rtype[type]][SS1[i-1]][SS2[j+1]]+iext_s+2*iext_ass+di4+dj2, c[i][j]);
      /**
      *** nx2 or more
      **/
      c[i][j]=MIN2(in[i - 2][j+4]+P->mismatchI[rtype[type]][SS1[i-1]][SS2[j+1]]+iext_s+2*iext_ass+di2+dj4, c[i][j]);
      /**
      *** nx1 n>2
      **/
      c[i][j]=MIN2(inx[i - 3][j+1]+P->mismatch1nI[rtype[type]][SS1[i-1]][SS2[j+1]]+iext_ass+iext_ass+di3+dj1, c[i][j]);
      /**
      *** 1xn n>2
      **/
      c[i][j]=MIN2(iny[i - 1][j+3]+P->mismatch1nI[rtype[type]][SS1[i-1]][SS2[j+1]]+iext_ass+iext_ass+dj3+di1, c[i][j]);
      /**
      *** nx0 n>1
      **/
      int bAU;
      bAU=(type>2?P->TerminalAU:0);
      c[i][j]=MIN2(bx[i - 2][j+1]+di2+dj1+bext+bAU, c[i][j]);
      /**
      *** 0xn n>1
      **/
      c[i][j]=MIN2(by[i - 1][j+2]+di1+dj2+bext+bAU, c[i][j]);
      /*
      remove this line printf("%d\t",c[i][j]);
      */
      temp=min_colonne;
      min_colonne=MIN2(c[i][j]+ E_ExtLoop(rtype[type], SS2[j-1], SS1[i+1], P),min_colonne);
      if(temp>min_colonne){
        min_j_colonne=j;
      }
      /* ---------------------------------------------------------------------end update */
    }
    if(max>=min_colonne){
      max=min_colonne;
      max_pos=i;
      max_pos_j=min_j_colonne;
    }
    i++;
    /*
    remove this line printf("\n");
    */
  }
  Emin=max;
  if(Emin>threshold){
    free(S1); free(S2); free(SS1); free(SS2);
    for (i=0; i<=n3; i++) {
      free(c[i]);
      free(in[i]);
      free(bx[i]);
      free(by[i]);
      free(inx[i]);
      free(iny[i]);
    }
    for (i=0; i<=3; i++) {
      free(DJ[i]);
    }
    free(c);free(in);free(bx);free(by);free(inx);free(iny);free(DJ);
    mfe.energy=0;
    mfe.structure=NULL;
    return mfe;
  }
  i_min=max_pos;
  j_min=max_pos_j;
  int dGe, dGeplex, dGx, dGy;
  dGe=dGeplex=dGx=dGy=0;
  /* printf("MAX fduplexfold_XS %d\n",Emin); */
  struc = fbacktrack_XS(i_min, j_min, access_s1, access_s2, i_pos, j_pos,il_a, il_b,b_a,b_b,&dGe, &dGeplex, &dGx, &dGy);

  l1 = strchr(struc, '&')-struc;
  int size;
  size=strlen(struc)-1;
  int lengthx; int endx; int lengthy; int endy;
  lengthx=l1;
  lengthx-=(struc[0]=='.'?1:0);
  lengthx-=(struc[l1-1]=='.'?1:0);
  endx=(i_pos-(n3-i_min));
  lengthy=size-l1;
  lengthy-=(struc[size]=='.'?1:0);
  lengthy-=(struc[l1+1]=='.'?1:0);
  endy=j_pos+j_min+lengthy -22;
  if (i_min<n3-10) i_min++;
  if (j_min>11 ) j_min--;
  mfe.i = i_min;
  mfe.j = j_min;
  mfe.ddG = (double) Emin*0.01;
  mfe.structure = struc;
  mfe.energy_backtrack = (double) dGe * 0.01;
  mfe.energy = (double) dGeplex * 0.01;
  mfe.opening_backtrack_x = (double) dGx * 0.01;
  mfe.opening_backtrack_y = (double) dGy * 0.01;
  mfe.dG1=0;/* !remove access to complete access array (double) access_s1[lengthx][endx+10] * 0.01; */
  mfe.dG2=0;/* !remove access to complete access array (double) access_s2[lengthy][endy+10] * 0.01; */
  free(S1); free(S2); free(SS1); free(SS2);
  for (i=0; i<=n3; i++) {
    free(c[i]);
    free(in[i]);
    free(bx[i]);
    free(by[i]);
    free(inx[i]);
    free(iny[i]);
  }
  for (i=0; i<=3; i++) {
    free(DJ[i]);
  }
  free(DJ);
  free(c);free(in);free(bx);free(by);free(iny);free(inx);
  return mfe;
}

PRIVATE char *fbacktrack_XS(int i, int j, const int** access_s1, const int** access_s2, const int i_pos, const int j_pos,const int il_a, const int il_b, const int b_a, const int b_b, int *dG, int *dGplex, int *dGx, int *dGy) {
  /* backtrack structure going backwards from i, and forwards from j
     return structure in bracket notation with & as separator */
  int k, l, type, type2, E, traced, i0, j0;
  char *st1, *st2, *struc;
  int bopen=b_b;
  int bext=b_a;
  int iopen=il_b;
  int iext_s=2*il_a;/* iext_s 2 nt nucleotide extension of interior loop, on i and j side */
  int iext_ass=50+il_a;/* iext_ass assymetric extension of interior loop, either on i or on j side. */
  st1 = (char *) vrna_alloc(sizeof(char)*(n3+1));
  st2 = (char *) vrna_alloc(sizeof(char)*(n4+1));
  i0=MIN2(i+1,n3-10); j0=MAX2(j-1,11);
  int state;
  state=1; /* we start backtracking from a a pair , i.e. c-matrix */
  /* state 1 -> base pair, c
     state 2 -> interior loop, in
     state 3 -> bx loop, bx
     state 4 -> by loop, by
  */
  traced=1;
  k=i; l=j; /* stores the i,j information for subsequence usage see * */
  int idiff,jdiff;
  /**
  *** (type>2?P->TerminalAU:0)+P->dangle3[rtype[type]][SS1[i+1]]+P->dangle5[rtype[type]][SS2[j-1]];
  **/

  int maxPenalty[4];
  vrna_md_t   md;

  set_model_details(&md);

  if ((!P) || (fabs(P->temperature - temperature)>1e-6)){
    update_dfold_params();
    if(P)
      free(P);
    P = vrna_params(&md);
    make_pair_matrix();
  }
  maxPenalty[0]=(int) -1*P->stack[2][2]/2;
  maxPenalty[1]=(int) -1*P->stack[2][2];
  maxPenalty[2]=(int) -3*P->stack[2][2]/2;
  maxPenalty[3]=(int) -2*P->stack[2][2];

  type = pair[S1[i]][S2[j]];
  *dG+= E_ExtLoop(rtype[type], SS2[j-1] , SS1[i+1] , P);
  *dGplex=*dG;

  while (i>10 && j<=n4-9 && traced) {
    int di1,di2,di3,di4;
    idiff=i_pos-(n3-10-i);
    di1 = 0.5*(access_s1[5][idiff+4] - access_s1[4][idiff+4] + access_s1[5][idiff]   - access_s1[4][idiff-1]);
    di2 = 0.5*(access_s1[5][idiff+3] - access_s1[4][idiff+3] + access_s1[5][idiff-1] - access_s1[4][idiff-2]) + di1;
    di3 = 0.5*(access_s1[5][idiff+2] - access_s1[4][idiff+2] + access_s1[5][idiff-2] - access_s1[4][idiff-3]) + di2;
    di4 = 0.5*(access_s1[5][idiff+1] - access_s1[4][idiff+1] + access_s1[5][idiff-3] - access_s1[4][idiff-4]) + di3;
/*
    di1 = access_s1[5][idiff]   - access_s1[4][idiff-1];
    di2 = access_s1[5][idiff-1] - access_s1[4][idiff-2] + di1;
    di3 = access_s1[5][idiff-2] - access_s1[4][idiff-3] + di2;
    di4 = access_s1[5][idiff-3] - access_s1[4][idiff-4] + di3;
    di1=MIN2(di1,maxPenalty[0]);
    di2=MIN2(di2,maxPenalty[1]);
    di3=MIN2(di3,maxPenalty[2]);
    di4=MIN2(di4,maxPenalty[3]);
*/
    int dj1,dj2,dj3,dj4;
    jdiff=j_pos+j-11;
    dj1=0.5*(access_s2[5][jdiff+4] - access_s2[4][jdiff+4] + access_s2[5][jdiff]  -access_s2[4][jdiff-1]  );
    dj2=0.5*(access_s2[5][jdiff+5] - access_s2[4][jdiff+5] + access_s2[5][jdiff+1]-access_s2[4][jdiff]  ) + dj1;
    dj3=0.5*(access_s2[5][jdiff+6] - access_s2[4][jdiff+6] + access_s2[5][jdiff+2]-access_s2[4][jdiff+1]) + dj2;
    dj4=0.5*(access_s2[5][jdiff+7] - access_s2[4][jdiff+7] + access_s2[5][jdiff+3]-access_s2[4][jdiff+2]) + dj3;



/*
    dj1 = access_s2[5][jdiff+4] - access_s2[4][jdiff+4];
    dj2 = access_s2[5][jdiff+5] - access_s2[4][jdiff+5] + dj1;
    dj3 = access_s2[5][jdiff+6] - access_s2[4][jdiff+6] + dj2;
    dj4 = access_s2[5][jdiff+7] - access_s2[4][jdiff+7] + dj3;
    dj1=MIN2(dj1,maxPenalty[0]);
    dj2=MIN2(dj2,maxPenalty[1]);
    dj3=MIN2(dj3,maxPenalty[2]);
    dj4=MIN2(dj4,maxPenalty[3]);
*/
    traced=0;
    switch(state){
    case 1:
      type = pair[S1[i]][S2[j]];
      int bAU;
      bAU=(type>2?P->TerminalAU:0);
      if(!type) vrna_message_error("backtrack failed in fold duplex");
      type2=pair[S1[i-1]][S2[j+1]];
      if(type2 && c[i][j]== (c[i - 1][j+1]+P->stack[rtype[type]][type2]+di1+dj1)){
        k=i-1;
        l=j+1;
        (*dG)+=E_IntLoop(i-k-1, l-j-1, type2, rtype[type],SS1[k+1], SS2[l-1], SS1[i-1], SS2[j+1],P);
        *dGplex+=E_IntLoop(i-k-1, l-j-1, type2, rtype[type],SS1[k+1], SS2[l-1], SS1[i-1], SS2[j+1],P);
        *dGx+=di1;
        *dGy+=dj1;
        st1[i-1] = '(';
        st2[j-1] = ')';
        i=k;
        j=l;
        state=1;
        traced=1;
        break;
      }
      type2=pair[S1[i-1]][S2[j+2]];
      if(type2 && c[i][j]==(c[i - 1][j+2]+P->bulge[1]+P->stack[rtype[type]][type2]+di1+dj2)){
        k=i-1;
        l=j+2;
        *dG+=E_IntLoop(i-k-1, l-j-1, type2, rtype[type],SS1[k+1], SS2[l-1], SS1[i-1], SS2[j+1],P);
        *dGplex+=E_IntLoop(i-k-1, l-j-1, type2, rtype[type],SS1[k+1], SS2[l-1], SS1[i-1], SS2[j+1],P);
        *dGx+=di1;
        *dGy+=dj2;
        st1[i-1] = '(';
        st2[j-1] = ')';
        i=k;
        j=l;
        state=1;
        traced=1;
        break;
      }
      type2=pair[S1[i-2]][S2[j+1]];
      if(type2 && c[i][j]==(c[i - 2][j+1]+P->bulge[1]+P->stack[type2][rtype[type]]+di2+dj1)){
        k=i-2;
        l=j+1;
        *dG+=E_IntLoop(i-k-1, l-j-1, type2, rtype[type],SS1[k+1], SS2[l-1], SS1[i-1], SS2[j+1],P);
        *dGplex+=E_IntLoop(i-k-1, l-j-1, type2, rtype[type],SS1[k+1], SS2[l-1], SS1[i-1], SS2[j+1],P);
        *dGx+=di2;
        *dGy+=dj1;
        st1[i-1] = '(';
        st2[j-1] = ')';
        i=k;
        j=l;
        state=1;
        traced=1;
        break;
      }
      type2=pair[S1[i-2]][S2[j+2]];
      if(type2 && c[i][j]==(c[i - 2][j+2]+P->int11[type2][rtype[type]][SS1[i-1]][SS2[j+1]]+di2+dj2)){
        k=i-2;
        l=j+2;
        *dG+=E_IntLoop(i-k-1, l-j-1, type2, rtype[type],SS1[k+1], SS2[l-1], SS1[i-1], SS2[j+1],P);
        *dGplex+=E_IntLoop(i-k-1, l-j-1, type2, rtype[type],SS1[k+1], SS2[l-1], SS1[i-1], SS2[j+1],P);
        *dGx+=di2;
        *dGy+=dj2;
        st1[i-1] = '(';
        st2[j-1] = ')';
        i=k;
        j=l;
        state=1;
        traced=1;
        break;
      }
      type2 = pair[S1[i-3]][S2[j+3]];
      if(type2 && c[i][j]==(c[i - 3][j+3]+P->int22[type2][rtype[type]][SS1[i-2]][SS1[i-1]][SS2[j+1]][SS2[j+2]]+di3+dj3)){
        k=i-3;
        l=j+3;
        *dG+=E_IntLoop(i-k-1, l-j-1, type2, rtype[type],SS1[k+1], SS2[l-1], SS1[i-1], SS2[j+1],P);
        *dGplex+=E_IntLoop(i-k-1, l-j-1, type2, rtype[type],SS1[k+1], SS2[l-1], SS1[i-1], SS2[j+1],P);
        *dGx+=di3;
        *dGy+=dj3;
        st1[i-1] = '(';
        st2[j-1] = ')';
        i=k;
        j=l;
        state=1;
        traced=1;
        break;
      }
      type2 = pair[S1[i-3]][S2[j+2]];
      if(type2 && c[i][j]==(c[i - 3][j+2]+P->int21[rtype[type]][type2][SS2[j+1]][SS1[i-2]][SS1[i-1]]+di3+dj2)){
        k=i-3;
        l=j+2;
        *dG+=E_IntLoop(i-k-1, l-j-1, type2, rtype[type],SS1[k+1], SS2[l-1], SS1[i-1], SS2[j+1],P);
        *dGplex+=E_IntLoop(i-k-1, l-j-1, type2, rtype[type],SS1[k+1], SS2[l-1], SS1[i-1], SS2[j+1],P);
        *dGx+=di3;
        *dGy+=dj2;
        st1[i-1] = '(';
        st2[j-1] = ')';
        i=k;
        j=l;
        state=1;
        traced=1;
        break;
      }
      type2 = pair[S1[i-2]][S2[j+3]];
      if(type2 && c[i][j]==(c[i - 2][j+3]+P->int21[type2][rtype[type]][SS1[i-1]][SS2[j+1]][SS2[j+2]]+di2+dj3)){
        k=i-2;
        l=j+3;
        *dG+=E_IntLoop(i-k-1, l-j-1, type2, rtype[type],SS1[k+1], SS2[l-1], SS1[i-1], SS2[j+1],P);
        *dGplex+=E_IntLoop(i-k-1, l-j-1, type2, rtype[type],SS1[k+1], SS2[l-1], SS1[i-1], SS2[j+1],P);
        *dGx+=di2;
        *dGy+=dj3;
        st1[i-1] = '(';
        st2[j-1] = ')';
        i=k;
        j=l;
        state=1;
        traced=1;
        break;
      }
      type2 = pair[S1[i-4]][S2[j+3]];
      if(type2 && c[i][j]==(c[i - 4][j+3]+P->internal_loop[5]+P->ninio[2]+
                            P->mismatch23I[type2][SS1[i-3]][SS2[j+2]]+P->mismatch23I[rtype[type]][SS2[j+1]][SS1[i-1]]+di4+dj3)){
        k=i-4;
        l=j+3;
        *dG+=E_IntLoop(i-k-1, l-j-1, type2, rtype[type],SS1[k+1], SS2[l-1], SS1[i-1], SS2[j+1],P);
        *dGplex+=E_IntLoop(i-k-1, l-j-1, type2, rtype[type],SS1[k+1], SS2[l-1], SS1[i-1], SS2[j+1],P);
        *dGx+=di2;
        *dGy+=dj3;
        st1[i-1] = '(';
        st2[j-1] = ')';
        i=k;
        j=l;
        state=1;
        traced=1;
        break;
      }
      type2 = pair[S1[i-3]][S2[j+4]];
      if(type2 && c[i][j]==(c[i - 3][j+4]+P->internal_loop[5]+P->ninio[2]+
                            P->mismatch23I[type2][SS1[i-2]][SS2[j+3]]+P->mismatch23I[rtype[type]][SS2[j+1]][SS1[i-1]]+di3+dj4)){
        k=i-3;
        l=j+4;
        *dG+=E_IntLoop(i-k-1, l-j-1, type2, rtype[type],SS1[k+1], SS2[l-1], SS1[i-1], SS2[j+1],P);
        *dGplex+=E_IntLoop(i-k-1, l-j-1, type2, rtype[type],SS1[k+1], SS2[l-1], SS1[i-1], SS2[j+1],P);
        *dGx+=di2;
        *dGy+=dj3;
        st1[i-1] = '(';
        st2[j-1] = ')';
        i=k;
        j=l;
        state=1;
        traced=1;
        break;
      }
      if(c[i][j]==(in[i - 3][j+3]+P->mismatchI[rtype[type]][SS1[i-1]][SS2[j+1]]+di3+dj3+2*iext_s)){
        k=i;
        l=j;
        *dGplex+=P->mismatchI[rtype[type]][SS1[i-1]][SS2[j+1]]+2*iext_s;
        *dGx+=di3;
        *dGy+=dj3;
        st1[i-1] = '(';
        st2[j-1] = ')';
        i=i-3;
        j=j+3;
        state=2;
        traced=1;
        break;
      }
      if(c[i][j]==(in[i - 4][j+2]+P->mismatchI[rtype[type]][SS1[i-1]][SS2[j+1]]+di4+dj2+iext_s+2*iext_ass)){
        k=i;
        l=j;
        *dGplex+=P->mismatchI[rtype[type]][SS1[i-1]][SS2[j+1]]+iext_s+2*iext_ass;
        *dGx+=di4;
        *dGy+=dj2;
        st1[i-1] = '(';
        st2[j-1] = ')';
        i=i-4;
        j=j+2;
        state=2;
        traced=1;
        break;
      }
      if(c[i][j]==(in[i - 2][j+4]+P->mismatchI[rtype[type]][SS1[i-1]][SS2[j+1]]+di2+dj4+iext_s+2*iext_ass)){
        k=i;
        l=j;
        *dGplex+=P->mismatchI[rtype[type]][SS1[i-1]][SS2[j+1]]+iext_s+2*iext_ass;
        *dGx+=di2;
        *dGy+=dj4;
        st1[i-1] = '(';
        st2[j-1] = ')';
        i=i-2;
        j=j+4;
        state=2;
        traced=1;
        break;
      }
      if(c[i][j]==(inx[i - 3][j+1]+P->mismatch1nI[rtype[type]][SS1[i-1]][SS2[j+1]]+iext_ass+iext_ass+di3+dj1)){
        k=i;
        l=j;
        *dGplex+=P->mismatch1nI[rtype[type]][SS1[i-1]][SS2[j+1]]+iext_ass+iext_ass+di3+dj1;
        *dGx+=di3;
        *dGy+=dj1;
        st1[i-1] = '(';
        st2[j-1] = ')';
        i=i-3;
        j=j+1;
        state=5;
        traced=1;
        break;
      }
      if(c[i][j]==(iny[i - 1][j+3]+P->mismatch1nI[rtype[type]][SS1[i-1]][SS2[j+1]]+iext_ass+iext_ass+di1+dj3)){
        k=i;
        l=j;
        *dGplex+=P->mismatch1nI[rtype[type]][SS1[i-1]][SS2[j+1]]+iext_ass+iext_ass+di1+dj3;
        *dGx+=di1;
        *dGy+=dj3;
        st1[i-1] = '(';
        st2[j-1] = ')';
        i=i-1;
        j=j+3;
        state=6;
        traced=1;
        break;
      }
      if(c[i][j]==(bx[i - 2][j+1]+di2+dj1+bext+bAU)){
        k=i;
        l=j;
        st1[i-1] = '(';
        st2[j-1] = ')';
        *dGplex+=bext+bAU;
        *dGx+=di2;
        *dGy+=dj1;
        i=i-2;
        j=j+1;
        state=3;
        traced=1;
        break;
      }
      if(c[i][j]==(by[i - 1][j+2]+di1+dj2+bext+bAU)){
        k=i;
        l=j;
        *dGplex+=bext+bAU;
        *dGx+=di1;
        *dGy+=dj2;
        st1[i-1] = '(';
        st2[j-1] = ')';
        i=i-1;
        j=j+2;
        state=4;
        traced=1;
        break;
      }
      break;
    case 2:
      if(in[i][j]==(in[i - 1][j+1]+iext_s + di1 + dj1)){
        i--;
        j++;
        *dGplex+=iext_s;
        *dGx+=di1;
        *dGy+=dj1;
        state=2;
        traced=1;
        break;
      }
      if(in[i][j]==(in[i - 1][j]+iext_ass + di1)){
        i=i-1;
        *dGplex+=iext_ass;
        *dGx+=di1;
        state=2;
        traced=1;
        break;
      }
      if(in[i][j]==(in[i][j+1]+iext_ass + dj1)){
        j++;
        state=2;
        *dGy+=dj1;
        *dGplex+=iext_ass;
        traced=1;
        break;
      }
      type2=pair[SS2[j+1]][SS1[i-1]];
      if(type2 && in[i][j]==(c[i - 1][j+1]+P->mismatchI[type2][SS2[j]][SS1[i]]+iopen+iext_s + di1 +dj1)){
        *dGplex+=P->mismatchI[type2][SS2[j]][SS1[i]]+iopen+iext_s;
        int temp; temp=k; k=i-1; i=temp;
        temp=l; l=j+1; j=temp;
        type=pair[S1[i]][S2[j]];
        *dG+=E_IntLoop(i-k-1, l-j-1, type2, rtype[type],SS1[k+1], SS2[l-1], SS1[i-1], SS2[j+1],P);
        *dGx+=di1;
        *dGy+=dj1;
        i=k;
        j=l;
        state=1;
        traced=1;
        break;
      }
    case 3:
      if(bx[i][j]==(bx[i - 1][j]+bext+di1)){
        i--;
        *dGplex+=bext;
        *dGx+=di1;
        state=3;
        traced=1;
        break;
      }
      type2=pair[S2[j]][S1[i-1]];
      if(type2 && bx[i][j]==(c[i - 1][j]+bopen+bext+(type2>2?P->TerminalAU:0)+di1)){
        int temp; temp=k; k=i-1; i=temp;
        temp=l; l=j; j=temp;
        type=pair[S1[i]][S2[j]];
        *dG+=E_IntLoop(i-k-1, l-j-1, type2, rtype[type],SS1[k+1], SS2[l-1], SS1[i-1], SS2[j+1],P);
        *dGplex+=bopen+bext+(type2>2?P->TerminalAU:0);
        *dGx+=di1;
        i=k;
        j=l;
        state=1;
        traced=1;
        break;
      }
    case 4:
      if(by[i][j]==(by[i][j+1] + bext +dj1)){
        j++;
        *dGplex+=bext;
        state=4;
        traced=1;
        break;
      }
      type2=pair[S2[j+1]][S1[i]];
      if(type2 && by[i][j]==(c[i][j+1]+bopen+bext+(type2>2?P->TerminalAU:0) + dj1)){
        int temp; temp=k; k=i; i=temp;
        temp=l; l=j+1; j=temp;
        type=pair[S1[i]][S2[j]];
        *dG+=E_IntLoop(i-k-1, l-j-1, type2, rtype[type],SS1[k+1], SS2[l-1], SS1[i-1], SS2[j+1],P);
        *dGplex+=bopen+bext+(type2>2?P->TerminalAU:0);
        *dGy+=dj1;
        i=k;
        j=l;
        state=1;
        traced=1;
        break;
      }
    case 5:
      if(inx[i][j]==(inx[i-1][j]+iext_ass+di1)) {
        i--;
        *dGplex+=iext_ass;
        *dGx+=di1;
        state=5;
        traced=1;
        break;
      }
      type2=pair[S2[j+1]][S1[i-1]];
      if(type2 && inx[i][j]==(c[i-1][j+1]+P->mismatch1nI[type2][SS2[j]][SS1[i]]+iopen+iext_s+di1+dj1)){
        *dGplex+=P->mismatch1nI[type2][SS2[j]][SS1[i]]+iopen+iext_s;
        int temp; temp=k; k=i-1; i=temp;
        temp=l; l=j+1; j=temp;
        type=pair[S1[i]][S2[j]];
        *dG+=E_IntLoop(i-k-1, l-j-1, type2, rtype[type],SS1[k+1], SS2[l-1], SS1[i-1], SS2[j+1],P);
        *dGx+=di1;
        *dGy+=dj1;
        i=k;
        j=l;
        state=1;
        traced=1;
        break;
      }
    case 6:
      if(iny[i][j]==(iny[i][j+1]+iext_ass+dj1)) {
        j++;
        *dGplex+=iext_ass;
        *dGx+=dj1;
        state=6;
        traced=1;
        break;
      }
      type2=pair[S2[j+1]][S1[i-1]];
      if(type2 && iny[i][j]==(c[i-1][j+1]+P->mismatch1nI[type2][SS2[j]][SS1[i]]+iopen+iext_s+di1+dj1)){
        *dGplex+=P->mismatch1nI[type2][SS2[j]][SS1[i]]+iopen+iext_s;
        int temp; temp=k; k=i-1; i=temp;
        temp=l; l=j+1; j=temp;
        type=pair[S1[i]][S2[j]];
        *dG+=E_IntLoop(i-k-1, l-j-1, type2, rtype[type],SS1[k+1], SS2[l-1], SS1[i-1], SS2[j+1],P);
        *dGx+=di1;
        *dGy+=dj1;
        i=k;
        j=l;
        state=1;
        traced=1;
        break;
      }
    }
  }
  if (!traced) {
    idiff=i_pos-(n3-10-i);
    jdiff=j_pos+j-11;
    E=c[i][j];
    /**
    *** if (i>1) {E -= P->dangle5[type][SS1[i-1]]; *dG+=P->dangle5[type][SS1[i-1]];*dGplex+=P->dangle5[type][SS1[i-1]];}
    *** if (j<n4){E -= P->dangle3[type][SS2[j+1]]; *dG+=P->dangle3[type][SS2[j+1]];*dGplex+=P->dangle3[type][SS2[j+1]];}
    *** if (type>2) {E -= P->TerminalAU; *dG+=P->TerminalAU;*dGplex+=P->TerminalAU;}
    **/
    int correction;
    correction = E_ExtLoop(type, (i>1) ? SS1[i-1] : -1, (j<n4) ? SS2[j+1] : -1, P);
    *dG+=correction;
    *dGplex+=correction;
    E-=correction;

/*
 if (E != P->DuplexInit+access_s1[1][idiff]+access_s2[1][jdiff]) {
      vrna_message_error("backtrack failed in second fold duplex");
    }
*/
    if (E != P->DuplexInit) {
      vrna_message_error("backtrack failed in second fold duplex");
    }
    else{
      *dG+=P->DuplexInit;
      *dGplex+=P->DuplexInit;
      *dGx+=0;/* access_s1[1][idiff]; */
      *dGy+=0;/* access_s2[1][jdiff]; */
      st1[i-1]='(';
      st2[j-1]=')';
    }
  }
  if (i>11)  i--;
  if (j<n4-10) j++;
  struc = (char *) vrna_alloc(i0-i+1+j-j0+1+2);
  for (k=MAX2(i,1); k<=i0; k++) if (!st1[k-1]) st1[k-1] = '.';
  for (k=j0; k<=j; k++) if (!st2[k-1]) st2[k-1] = '.';
  strcpy(struc, st1+MAX2(i-1,0));
  strcat(struc, "&");
  strcat(struc, st2+j0-1);
  /* printf("%s %3d,%-3d : %3d,%-3d\n", struc, i,i0,j0,j);  */
  free(st1); free(st2);
  return struc;
}


duplexT ** Lduplexfold_XS(const char *s1, const char *s2, const int **access_s1, const int **access_s2, const int threshold, const int alignment_length, const int delta, const int fast, const int il_a, const int il_b, const int b_a, const int b_b)
{
  /**
  *** See variable definition in fduplexfold_XS
  **/
  int i, j;
  int bopen=b_b;
  int bext=b_a;
  int iopen=il_b;
  int iext_s=2*il_a;
  int iext_ass=50+il_a;
  int min_colonne=INF;
  int i_length;
  int max_pos;
  int max_pos_j;
  int min_j_colonne;
  int max=INF;
  int *position;
  int *position_j;
  int maxPenalty[4];
  int **DJ;
  /**
  *** 1D array corresponding to the standard 2d recursion matrix
  *** Makes the computation 20% faster
  **/
  int *SA;
  vrna_md_t   md;
  /**
  *** variable initialization
  **/
  n1 = (int) strlen(s1);
  n2 = (int) strlen(s2);
  /**
  *** Sequence encoding
  **/

  set_model_details(&md);

  if ((!P) || (fabs(P->temperature - temperature)>1e-6)) {
    update_dfold_params();
    if(P)
      free(P);
    P = vrna_params(&md);
    make_pair_matrix();
  }
  encode_seqs(s1,s2);
  /**
  *** Position of the high score on the target and query sequence
  **/
  position = (int *) vrna_alloc((delta+n1+3+delta) * sizeof(int));
  position_j= (int *) vrna_alloc((delta+n1+3+delta) * sizeof(int));
  /**
  *** extension penalty, computed only once, further reduce the computation time
  **/
  maxPenalty[0]=(int) -1*P->stack[2][2]/2;
  maxPenalty[1]=(int) -1*P->stack[2][2];
  maxPenalty[2]=(int) -3*P->stack[2][2]/2;
  maxPenalty[3]=(int) -2*P->stack[2][2];

  DJ=(int **) vrna_alloc(4*sizeof(int*));
  DJ[0]=(int *) vrna_alloc(n2*sizeof(int));
  DJ[1]=(int *) vrna_alloc(n2*sizeof(int));
  DJ[2]=(int *) vrna_alloc(n2*sizeof(int));
  DJ[3]=(int *) vrna_alloc(n2*sizeof(int));
  j=n2-9;
  while(--j>10){
    DJ[0][j] = 0.5*(access_s2[5][j+4] - access_s2[4][j+4] + access_s2[5][j]  -access_s2[4][j-1]  );
    DJ[1][j] = 0.5*(access_s2[5][j+5] - access_s2[4][j+5] + access_s2[5][j+1]-access_s2[4][j]  ) + DJ[0][j];
    DJ[2][j] = 0.5*(access_s2[5][j+6] - access_s2[4][j+6] + access_s2[5][j+2]-access_s2[4][j+1]) + DJ[1][j];
    DJ[3][j] = 0.5*(access_s2[5][j+7] - access_s2[4][j+7] + access_s2[5][j+3]-access_s2[4][j+2]) + DJ[2][j];
/*
    DJ[0][j] = access_s2[5][j+4] - access_s2[4][j+4]           ;
    DJ[1][j] = access_s2[5][j+5] - access_s2[4][j+5] + DJ[0][j];
    DJ[2][j] = access_s2[5][j+6] - access_s2[4][j+6] + DJ[1][j];
    DJ[3][j] = access_s2[5][j+7] - access_s2[4][j+7] + DJ[2][j];
    DJ[0][j] = MIN2(DJ[0][j],maxPenalty[0]);
    DJ[1][j] = MIN2(DJ[1][j],maxPenalty[1]);
    DJ[2][j] = MIN2(DJ[2][j],maxPenalty[2]);
    DJ[3][j] = MIN2(DJ[3][j],maxPenalty[3]);
*/
  }
  /**
  *** instead of having 4 2-dim arrays we use a unique 1-dim array
  *** The mapping 2d -> 1D is done based ont the macro
  *** LCI(i,j,l)      ((i     )*l + j)
  *** LINI(i,j,l)     ((i +  5)*l + j)
  *** LBXI(i,j,l)     ((i + 10)*l + j)
  *** LBYI(i,j,l)     ((i + 15)*l + j)
  *** LINIX(i,j,l)    ((i + 20)*l + j)
  *** LINIY(i,j,l)    ((i + 25)*l + j)
  ***
  *** SA has a length of 5 (number of columns we look back) *
  ***                  * 6 (number of structures we look at) *
  ***                  * length of the sequence
  **/

  SA=(int *) vrna_alloc(sizeof(int)*5*6*(n2+5));
  for(j=n2+4;j>=0;j--) {
    SA[(j*30)   ]=SA[(j*30)+1   ]=SA[(j*30)+2   ]=SA[(j*30)+3   ]=SA[(j*30)+4   ]=INF;
    SA[(j*30)+5 ]=SA[(j*30)+1+5 ]=SA[(j*30)+2+5 ]=SA[(j*30)+3+5 ]=SA[(j*30)+4+5 ]=INF;
    SA[(j*30)+10]=SA[(j*30)+1+10]=SA[(j*30)+2+10]=SA[(j*30)+3+10]=SA[(j*30)+4+10]=INF;
    SA[(j*30)+15]=SA[(j*30)+1+15]=SA[(j*30)+2+15]=SA[(j*30)+3+15]=SA[(j*30)+4+15]=INF;
    SA[(j*30)+20]=SA[(j*30)+1+20]=SA[(j*30)+2+20]=SA[(j*30)+3+20]=SA[(j*30)+4+20]=INF;
    SA[(j*30)+25]=SA[(j*30)+1+25]=SA[(j*30)+2+25]=SA[(j*30)+3+25]=SA[(j*30)+4+25]=INF;
  }

  i=10 ;
  i_length= n1 - 9  ;
  while(i < i_length) {
    int di1,di2,di3,di4;
    int idx=i%5;
    int idx_1=(i-1)%5;
    int idx_2=(i-2)%5;
    int idx_3=(i-3)%5;
    int idx_4=(i-4)%5;
    di1 = 0.5*(access_s1[5][i+4] - access_s1[4][i+4] + access_s1[5][i]   - access_s1[4][i-1]);
    di2 = 0.5*(access_s1[5][i+3] - access_s1[4][i+3] + access_s1[5][i-1] - access_s1[4][i-2]) + di1;
    di3 = 0.5*(access_s1[5][i+2] - access_s1[4][i+2] + access_s1[5][i-2] - access_s1[4][i-3]) + di2;
    di4 = 0.5*(access_s1[5][i+1] - access_s1[4][i+1] + access_s1[5][i-3] - access_s1[4][i-4]) + di3;
/*
    di1 = access_s1[5][i]   - access_s1[4][i-1];
    di2 = access_s1[5][i-1] - access_s1[4][i-2] + di1;
    di3 = access_s1[5][i-2] - access_s1[4][i-3] + di2;
    di4 = access_s1[5][i-3] - access_s1[4][i-4] + di3;
    di1=MIN2(di1,maxPenalty[0]);
    di2=MIN2(di2,maxPenalty[1]);
    di3=MIN2(di3,maxPenalty[2]);
    di4=MIN2(di4,maxPenalty[3]);
*/
    j=n2 - 9;
    while (--j > 9) {
      int dj1,dj2,dj3,dj4;
      dj1=DJ[0][j];
      dj2=DJ[1][j];
      dj3=DJ[2][j];
      dj4=DJ[3][j];
      int type2, type,temp;
      type  = pair[S1[i]][S2[j]];
      /**
      *** Start duplex
      **/
      /* SA[LCI(idx,j,n2)] = type ? P->DuplexInit + access_s1[1][i] + access_s2[1][j] : INF; */
      SA[LCI(idx,j,n2)] = type ? P->DuplexInit : INF;
      /**
      *** update lin bx by linx liny matrix
      **/
      type2=pair[S2[j+1]][S1[i-1]];
      /**
      *** start/extend interior loop
      **/
      SA[LINI(idx,j,n2)]=MIN2(SA[LCI(idx_1,j+1,n2)]+P->mismatchI[type2][SS2[j]][SS1[i]]+di1+dj1+iopen+iext_s,
                              SA[LINI(idx_1,j,n2)]+iext_ass + di1);

      /**
      *** start/extend nx1 target
      *** use same type2 as for in
      **/
      SA[LINIX(idx,j,n2)]=MIN2(SA[LCI(idx_1,j+1,n2)]+P->mismatch1nI[type2][SS2[j]][SS1[i]]+di1+dj1+iopen+iext_s,
                               SA[LINIX(idx_1,j,n2)]+iext_ass + di1);
      /**
      *** start/extend 1xn target
      *** use same type2 as for in
      **/
      SA[LINIY(idx,j,n2)]=MIN2(SA[LCI(idx_1,j+1,n2)]+P->mismatch1nI[type2][SS2[j]][SS1[i]]+di1+dj1+iopen+iext_s,
                               SA[LINIY(idx,j+1,n2)]+iext_ass + dj1);
      /**
      *** extend interior loop
      **/
      SA[LINI(idx,j,n2)]=MIN2(SA[LINI(idx,j,n2)],SA[LINI(idx,j+1,n2)]+iext_ass + dj1);
      SA[LINI(idx,j,n2)]=MIN2(SA[LINI(idx,j,n2)],SA[LINI(idx_1,j+1,n2)]+iext_s + di1 + dj1);
      /**
      *** start/extend bulge target
      **/
      type2=pair[S2[j]][S1[i-1]];
      SA[LBXI(idx,j,n2)]=MIN2(SA[LBXI(idx_1,j,n2)]+bext + di1, SA[LCI(idx_1,j,n2)]+bopen+bext+(type2>2?P->TerminalAU:0) + di1);
      /**
      *** start/extend bulge query
      **/
      type2=pair[S2[j+1]][S1[i]];
      SA[LBYI(idx,j,n2)]=MIN2(SA[LBYI(idx,j+1,n2)]+bext + dj1 , SA[LCI(idx,j+1,n2)]+bopen+bext+(type2>2?P->TerminalAU:0)+ dj1);
      /**
      ***end update recursion
      **/
      if(!type){continue;}
      /**
      *** stack extension
      **/
      SA[LCI(idx,j,n2)]+= E_ExtLoop(type, SS1[i-1] , SS2[j+1], P);
      /**
      *** stack extension
      **/
      if((type2=pair[S1[i-1]][S2[j+1]]))
        SA[LCI(idx,j,n2)]=MIN2(SA[LCI(idx_1,j+1,n2)]+P->stack[rtype[type]][type2]+di1+dj1, SA[LCI(idx,j,n2)]);
      /**
      *** 1x0 / 0x1 stack extension
      **/
      if((type2=pair[S1[i-1]][S2[j+2]]))
        SA[LCI(idx,j,n2)]=MIN2(SA[LCI(idx_1,j+2,n2)]+P->bulge[1]+P->stack[rtype[type]][type2]+di1+dj2, SA[LCI(idx,j,n2)]);
      if((type2=pair[S1[i-2]][S2[j+1]]))
        SA[LCI(idx,j,n2)]=MIN2(SA[LCI(idx_2,j+1,n2)]+P->bulge[1]+P->stack[type2][rtype[type]]+di2+dj1, SA[LCI(idx,j,n2)]);
      /**
      *** 1x1 / 2x2 stack extension
      **/
      if((type2=pair[S1[i-2]][S2[j+2]]))
        SA[LCI(idx,j,n2)]=MIN2(SA[LCI(idx_2,j+2,n2)]+P->int11[type2][rtype[type]][SS1[i-1]][SS2[j+1]]+di2+dj2, SA[LCI(idx,j,n2)]);
      if((type2 = pair[S1[i-3]][S2[j+3]]))
        SA[LCI(idx,j,n2)]=MIN2(SA[LCI(idx_3,j+3,n2)]+P->int22[type2][rtype[type]][SS1[i-2]][SS1[i-1]][SS2[j+1]][SS2[j+2]]+di3+dj3, SA[LCI(idx,j,n2)]);
      /**
      *** 1x2 / 2x1 stack extension
      *** E_IntLoop(1,2,type2, rtype[type],SS1[i-1], SS2[j+2], SS1[i-1], SS2[j+1], P) corresponds to
      *** P->int21[rtype[type]][type2][SS2[j+2]][SS1[i-1]][SS1[i-1]]
      **/
      if((type2 = pair[S1[i-3]][S2[j+2]]))
        SA[LCI(idx,j,n2)]=MIN2(SA[LCI(idx_3,j+2,n2)]+P->int21[rtype[type]][type2][SS2[j+1]][SS1[i-2]][SS1[i-1]]+di3+dj2, SA[LCI(idx,j,n2)]);
      if((type2 = pair[S1[i-2]][S2[j+3]]))
        SA[LCI(idx,j,n2)]=MIN2(SA[LCI(idx_2,j+3,n2)]+P->int21[type2][rtype[type]][SS1[i-1]][SS2[j+1]][SS2[j+2]]+di2+dj3, SA[LCI(idx,j,n2)]);
      /**
      *** 2x3 / 3x2 stack extension
      **/
      if((type2 = pair[S1[i-4]][S2[j+3]]))
        SA[LCI(idx,j,n2)]=MIN2(SA[LCI(idx_4,j+3,n2)]+P->internal_loop[5]+P->ninio[2]+
                               P->mismatch23I[type2][SS1[i-3]][SS2[j+2]]+P->mismatch23I[rtype[type]][SS2[j+1]][SS1[i-1]]+di4+dj3, SA[LCI(idx,j,n2)]);
      if((type2 = pair[S1[i-3]][S2[j+4]]))
        SA[LCI(idx,j,n2)]=MIN2(SA[LCI(idx_3,j+4,n2)]+P->internal_loop[5]+P->ninio[2]+
                               P->mismatch23I[type2][SS1[i-2]][SS2[j+3]]+P->mismatch23I[rtype[type]][SS2[j+1]][SS1[i-1]]+di3+dj4, SA[LCI(idx,j,n2)]);
      /**
      *** So now we have to handle 1x3, 3x1, 3x3, and mxn m,n > 3
      **/
      /**
      *** 3x3 or more
      **/
      SA[LCI(idx,j,n2)]=MIN2(SA[LINI(idx_3,j+3,n2)]+P->mismatchI[rtype[type]][SS1[i-1]][SS2[j+1]]+2*iext_s+di3+dj3,SA[LCI(idx,j,n2)]);
      /**
      *** 2xn or more
      **/
      SA[LCI(idx,j,n2)]=MIN2(SA[LINI(idx_4,j+2,n2)]+P->mismatchI[rtype[type]][SS1[i-1]][SS2[j+1]]+iext_s+2*iext_ass+di4+dj2, SA[LCI(idx,j,n2)]);
      /**
      *** nx2 or more
      **/
      SA[LCI(idx,j,n2)]=MIN2(SA[LINI(idx_2,j+4,n2)]+P->mismatchI[rtype[type]][SS1[i-1]][SS2[j+1]]+iext_s+2*iext_ass+di2+dj4, SA[LCI(idx,j,n2)]);
      /**
      *** nx1 n>2
      **/
      SA[LCI(idx,j,n2)]=MIN2(SA[LINIX(idx_3,j+1,n2)]+P->mismatch1nI[rtype[type]][SS1[i-1]][SS2[j+1]]+iext_ass+iext_ass+di3+dj1, SA[LCI(idx,j,n2)]);
      /**
      *** 1xn n>2
      **/
      SA[LCI(idx,j,n2)]=MIN2(SA[LINIY(idx_1,j+3,n2)]+P->mismatch1nI[rtype[type]][SS1[i-1]][SS2[j+1]]+iext_ass+iext_ass+dj3+di1, SA[LCI(idx,j,n2)]);
      /**
      *** nx0 n>1
      **/
      int bAU;
      bAU=(type>2?P->TerminalAU:0);
      SA[LCI(idx,j,n2)]=MIN2(SA[LBXI(idx_2,j+1,n2)]+di2+dj1+bext+bAU,SA[LCI(idx,j,n2)]);
      /**
      *** 0xn n>1
      **/
      SA[LCI(idx,j,n2)]=MIN2(SA[LBYI(idx_1,j+2,n2)]+di1+dj2+bext+bAU,SA[LCI(idx,j,n2)]);
      temp=min_colonne;
      /**
      *** (type>2?P->TerminalAU:0)+
      *** P->dangle3[rtype[type]][SS1[i+1]]+
      *** P->dangle5[rtype[type]][SS2[j-1]],
      **/
      /* remove this line printf("LCI %d:%d %d\t",i,j,SA[LCI(idx,j,n2)]); */
      /* remove this line printf("LI %d:%d %d\t",i,j, SA[LINI(idx,j,n2)]); */
      min_colonne=MIN2(SA[LCI(idx,j,n2)]+E_ExtLoop(rtype[type], SS2[j-1] , SS1[i+1] , P), min_colonne);

      if(temp>min_colonne){
        min_j_colonne=j;
      }

      /* ---------------------------------------------------------------------end update */
    }
    if(max>=min_colonne){
      max=min_colonne;
      max_pos=i;
      max_pos_j=min_j_colonne;
    }
    position[i+delta]=min_colonne;min_colonne=INF;
    position_j[i+delta]=min_j_colonne;
    /* remove this line printf("\n"); */
    i++;
  }
  /* printf("MAX: %d",max); */
  free(S1); free(S2); free(SS1); free(SS2);free(SA);
  if(max<threshold){
    find_max_XS(position, position_j, delta, threshold, alignment_length, s1, s2, access_s1, access_s2, fast,il_a, il_b,b_a, b_b);
  }
  if(max<INF){
    plot_max_XS(max, max_pos, max_pos_j, alignment_length, s1, s2, access_s1, access_s2,fast, il_a, il_b, b_a, b_b);
  }
  for (i=0; i<=3; i++) {
    free(DJ[i]);
  }
  free(DJ);
  free(position);
  free(position_j);
  return NULL;
}

PRIVATE void find_max_XS(const int *position, const int *position_j,const int delta, const int threshold, const int alignment_length, const char *s1, const char *s2, const int **access_s1, const int **access_s2, const int fast,const int il_a, const int il_b, const int b_a, const int b_b){
  int pos=n1-9;
  if(fast==1){
    while(10 < pos--){
      int temp_min=0;
      if(position[pos+delta]<(threshold)){
        int search_range;
        search_range=delta+1;
        while(--search_range){
          if(position[pos+delta-search_range]<=position[pos+delta-temp_min]){
            temp_min=search_range;
          }
        }
        pos-=temp_min;
        int max_pos_j;
        max_pos_j=position_j[pos+delta];
        int max;
        max=position[pos+delta];
        printf("target upper bound %d: query lower bound %d  (%5.2f) \n", pos-10, max_pos_j-10, ((double)max)/100);
        pos=MAX2(10,pos+temp_min-delta);
      }
    }
  }
  else if(fast==2){
    pos=n1-9;
    while(10 < pos--){
      int temp_min=0;
      if(position[pos+delta]<(threshold)){
        int search_range;
        search_range=delta+1;
        while(--search_range){
          if(position[pos+delta-search_range]<=position[pos+delta-temp_min]){
            temp_min=search_range;
          }
        }
        pos-=temp_min;
        int max_pos_j;
        max_pos_j=position_j[pos+delta];
        /* max_pos_j und pos entsprechen die realen position
            in der erweiterten sequenz.
           pos=1 -> position 1 in the sequence (and not 0 like in C)
           max_pos_j -> position 1 in the sequence ( not 0 like in C)
        */
        int alignment_length2; alignment_length2 = MIN2(n1,n2);
        int begin_t=MAX2(11, pos-alignment_length2+1);/* 10 */
        int end_t  =MIN2(n1-10, pos+1);
        int begin_q=MAX2(11, max_pos_j-1); /* 10 */
        int end_q  =MIN2(n2-10, max_pos_j+alignment_length2-1);
        char *s3 = (char*) vrna_alloc(sizeof(char)*(end_t - begin_t +2 + 20));
        char *s4 = (char*) vrna_alloc(sizeof(char)*(end_q - begin_q +2 + 20));
        strcpy(s3,"NNNNNNNNNN");strcpy(s4,"NNNNNNNNNN");
        strncat(s3, (s1+begin_t-1),  end_t - begin_t +1);
        strncat(s4, (s2+begin_q-1) , end_q - begin_q +1);
        strcat(s3,"NNNNNNNNNN");strcat(s4,"NNNNNNNNNN");
        s3[end_t -begin_t +1 +20 ]='\0';
        s4[end_q -begin_q +1 +20]='\0';
        duplexT test;
        test = fduplexfold_XS(s3, s4, access_s1, access_s2, end_t, begin_q,threshold, il_a, il_b, b_a, b_b);
        if(test.energy * 100 < threshold){
          int l1=strchr(test.structure, '&')-test.structure;
          printf(" %s %3d,%-3d : %3d,%-3d (%5.2f = %5.2f + %5.2f + %5.2f) [%5.2f] i:%d,j:%d <%5.2f>\n", test.structure,
                 begin_t-10+test.i-l1-10,
                 begin_t-10+test.i-1-10,
                 begin_q-10 + test.j-1-10 ,
                 (begin_q -11) + test.j + (int)strlen(test.structure)-l1-2-10,
                 test.ddG, test.energy, test.opening_backtrack_x, test.opening_backtrack_y, test.energy_backtrack,
                 pos-10, max_pos_j-10, ((double) position[pos+delta])/100);
          pos=MAX2(10,pos+temp_min-delta);
          free(test.structure);
        }
        free(s3);free(s4);
      }
    }
  }
  else{
    pos=n1-9;
    while( pos-- > 10 ){
      int temp_min=0;
      if(position[pos+delta]<(threshold)){
        int search_range;
        search_range=delta+1;
        while(--search_range){
          if(position[pos+delta-search_range]<=position[pos+delta-temp_min]){
            temp_min=search_range;
          }
        }
        pos-=temp_min; /* position on i */
        int max_pos_j;
        max_pos_j=position_j[pos+delta]; /* position on j */
        int begin_t=MAX2(11,pos-alignment_length);
        int end_t  =MIN2(n1-10, pos+1);
        int begin_q=MAX2(11,max_pos_j-1);
        int end_q  =MIN2(n2-10,max_pos_j+alignment_length-1);
        int i_flag;
        int j_flag;
        i_flag = (end_t   ==  pos+1?1:0);
        j_flag = (begin_q == max_pos_j-1?1:0);
        char *s3 = (char*) vrna_alloc(sizeof(char)*(end_t - begin_t +2));
        char *s4 = (char*) vrna_alloc(sizeof(char)*(end_q - begin_q +2));
        strncpy(s3, (s1+begin_t),  end_t - begin_t+1);
        strncpy(s4, (s2+begin_q) , end_q - begin_q+1);
        s3[end_t -begin_t +1 ]='\0';
        s4[end_q -begin_q +1 ]='\0';
        duplexT test;
        test = duplexfold_XS(s3,s4,access_s1,access_s2,pos, max_pos_j,threshold,i_flag,j_flag);
        if(test.energy * 100 < threshold){
          printf("%s %3d,%-3d : %3d,%-3d (%5.2f = %5.2f + %5.2f + %5.2f) i:%d,j:%d <%5.2f>\n", test.structure,
                 test.tb,test.te,test.qb,test.qe, test.ddG, test.energy, test.dG1, test.dG2, pos-10, max_pos_j-10, ((double) position[pos+delta])/100);
          pos=MAX2(10,pos+temp_min-delta);
        }
        free(s3);free(s4);
        free(test.structure);
      }
    }
  }
}

#if 0
PRIVATE int compare(const void *sub1, const void *sub2) {
  int d;
  if (((duplexT *) sub1)->ddG > ((duplexT *) sub2)->ddG)
    return 1;
  if (((duplexT *) sub1)->ddG < ((duplexT *) sub2)->ddG)
    return -1;
  d = ((duplexT *) sub1)->i - ((duplexT *) sub2)->i;
  if (d!=0) return d;
  return  ((duplexT *) sub1)->j - ((duplexT *) sub2)->j;
}
#endif

PRIVATE void plot_max_XS(const int max, const int max_pos, const int max_pos_j, const int alignment_length, const char *s1, const char *s2, const int ** access_s1, const int ** access_s2, const int fast,const int il_a, const int il_b, const int b_a, const int b_b)
{
  if(fast==1){
    printf("target upper bound %d: query lower bound %d (%5.2f)\n", max_pos-3, max_pos_j, ((double)max)/100);
  }
  else if(fast==2){
    int alignment_length2; alignment_length2 = MIN2(n1,n2);
    int begin_t=MAX2(11, max_pos-alignment_length2+1);/* 10 */
    int end_t  =MIN2(n1-10, max_pos+1);
    int begin_q=MAX2(11, max_pos_j-1); /* 10 */
    int end_q  =MIN2(n2-10, max_pos_j+alignment_length2-1);
    char *s3 = (char*) vrna_alloc(sizeof(char)*(end_t - begin_t +2 + 20));
    char *s4 = (char*) vrna_alloc(sizeof(char)*(end_q - begin_q +2 + 20));
    strcpy(s3,"NNNNNNNNNN");strcpy(s4,"NNNNNNNNNN");
    strncat(s3, (s1+begin_t-1),  end_t - begin_t +1);
    strncat(s4, (s2+begin_q-1) , end_q - begin_q +1);
    strcat(s3,"NNNNNNNNNN");strcat(s4,"NNNNNNNNNN");
    s3[end_t -begin_t +1 +20 ]='\0';
    s4[end_q -begin_q +1 +20]='\0';
    duplexT test;
    test = fduplexfold_XS(s3, s4, access_s1, access_s2, end_t, begin_q, INF, il_a, il_b, b_a, b_b);
    int l1=strchr(test.structure, '&')-test.structure;
    printf("%s %3d,%-3d : %3d,%-3d (%5.2f = %5.2f + %5.2f + %5.2f) [%5.2f] i:%d,j:%d <%5.2f>\n", test.structure,
           begin_t-10+test.i-l1-10,
           begin_t-10+test.i-1-10,
           begin_q-10 + test.j-1-10 ,
           (begin_q -11) + test.j + (int)strlen(test.structure)-l1-2-10,
           test.ddG, test.energy, test.opening_backtrack_x, test.opening_backtrack_y, test.energy_backtrack,
           max_pos-10, max_pos_j-10, (double) max/100);

    free(s3);free(s4);
    free(test.structure);
  }
  else{
    int begin_t=MAX2(11,max_pos-alignment_length);
    int end_t  =MIN2(n1-10, max_pos+1);
    int begin_q=MAX2(11, max_pos_j-1);
    int end_q  =MIN2(n2-10,max_pos_j+alignment_length-1);
    int i_flag;
    int j_flag;
    i_flag = (end_t   == max_pos+1?1:0);
    j_flag = (begin_q == max_pos_j-1?1:0);
    char *s3 = (char*) vrna_alloc(sizeof(char)*(end_t - begin_t +2)); /* +1 for \0 +1 for distance */
    char *s4 = (char*) vrna_alloc(sizeof(char)*(end_q - begin_q +2));

    strncpy(s3, (s1+begin_t-1),  end_t - begin_t+1);/* -1 to go from  */
    strncpy(s4, (s2+begin_q-1) , end_q - begin_q+1);/* -1 to go from  */
    s3[end_t -begin_t +1 ]='\0';/*  */
    s4[end_q -begin_q +1 ]='\0';
    duplexT test;
    test = duplexfold_XS(s3,s4,access_s1,access_s2,max_pos, max_pos_j,INF,i_flag,j_flag);
    printf("%s %3d,%-3d : %3d,%-3d (%5.2f = %5.2f + %5.2f + %5.2f) i:%d,j:%d <%5.2f>\n", test.structure,
           test.tb,test.te,test.qb,test.qe, test.ddG, test.energy, test.dG1, test.dG2, max_pos-10, max_pos_j - 10,(double) max/100);
    free(s3);free(s4);free(test.structure);
  }
}


/*---------------------------------------------------------duplexfold----------------------------------------------------------------------------------*/


PRIVATE duplexT duplexfold(const char *s1, const char *s2, const int extension_cost) {
  int i, j, l1, Emin=INF, i_min=0, j_min=0;
  char *struc;
  duplexT mfe;
  vrna_md_t md;

  n3 = (int) strlen(s1);
  n4 = (int) strlen(s2);

  set_model_details(&md);
  if ((!P) || (fabs(P->temperature - temperature)>1e-6)) {
    update_fold_params();
    if(P)
      free(P);
    P = vrna_params(&md);
    make_pair_matrix();
  }

  c = (int **) vrna_alloc(sizeof(int *) * (n3+1));
  for (i=0; i<=n3; i++) c[i] = (int *) vrna_alloc(sizeof(int) * (n4+1));
  encode_seqs(s1, s2);
  for (i=1; i<=n3; i++) {
    for (j=n4; j>0; j--) {
      int type, type2, E, k,l;
      type = pair[S1[i]][S2[j]];
      c[i][j] = type ? P->DuplexInit +2 * extension_cost: INF;
      if (!type) continue;
      /**
      ***       if (i>1)  c[i][j] += P->dangle5[type][SS1[i-1]]+ extension_cost;
      ***       if (j<n4) c[i][j] += P->dangle3[type][SS2[j+1]]+ extension_cost;
      ***       if (type>2) c[i][j] += P->TerminalAU;
      **/
      c[i][j] += E_ExtLoop(type, (i>1) ? SS1[i-1] : -1, (j<n4) ? SS2[j+1] : -1, P);
      for (k=i-1; k>0 && k>i-MAXLOOP-2; k--) {
        for (l=j+1; l<=n4; l++) {
          if (i-k+l-j-2>MAXLOOP) break;
          type2 = pair[S1[k]][S2[l]];
          if (!type2) continue;
          E = E_IntLoop(i-k-1, l-j-1, type2, rtype[type],
                        SS1[k+1], SS2[l-1], SS1[i-1], SS2[j+1],P)+(i-k+l-j)*extension_cost;
          c[i][j] = MIN2(c[i][j], c[k][l]+E);
        }
      }
      E = c[i][j];
      /**
      ***      if (i<n3) E += P->dangle3[rtype[type]][SS1[i+1]]+extension_cost;
      ***      if (j>1)  E += P->dangle5[rtype[type]][SS2[j-1]]+extension_cost;
      ***      if (type>2) E += P->TerminalAU;
      ***
      **/
      E += E_ExtLoop(rtype[type], (j > 1) ? SS2[j-1] : -1, (i<n3) ? SS1[i+1] : -1, P);
      if (E<Emin) {
        Emin=E; i_min=i; j_min=j;
      }
    }
  }
  struc = backtrack(i_min, j_min, extension_cost);
  if (i_min<n3) i_min++;
  if (j_min>1 ) j_min--;
  l1 = strchr(struc, '&')-struc;
  int size;
  size=strlen(struc)-1;
  Emin-= size * (extension_cost);
  mfe.i = i_min;
  mfe.j = j_min;
  mfe.energy = (double) Emin/100.;
  mfe.structure = struc;
  for (i=0; i<=n3; i++) free(c[i]);
  free(c);
  free(S1); free(S2); free(SS1); free(SS2);
  return mfe;
}

PRIVATE char *backtrack(int i, int j, const int extension_cost) {
  /* backtrack structure going backwards from i, and forwards from j
     return structure in bracket notation with & as separator */
  int k, l, type, type2, E, traced, i0, j0;
  char *st1, *st2, *struc;

  st1 = (char *) vrna_alloc(sizeof(char)*(n3+1));
  st2 = (char *) vrna_alloc(sizeof(char)*(n4+1));

  i0=MIN2(i+1,n3); j0=MAX2(j-1,1);

  while (i>0 && j<=n4) {
    E = c[i][j]; traced=0;
    st1[i-1] = '(';
    st2[j-1] = ')';
    type = pair[S1[i]][S2[j]];
    if (!type) vrna_message_error("backtrack failed in fold duplex");
    for (k=i-1; k>0 && k>i-MAXLOOP-2; k--) {
      for (l=j+1; l<=n4; l++) {
        int LE;
        if (i-k+l-j-2>MAXLOOP) break;
        type2 = pair[S1[k]][S2[l]];
        if (!type2) continue;
        LE = E_IntLoop(i-k-1, l-j-1, type2, rtype[type],
                       SS1[k+1], SS2[l-1], SS1[i-1], SS2[j+1],P)+(i-k+l-j)*extension_cost;
        if (E == c[k][l]+LE) {
          traced=1;
          i=k; j=l;
          break;
        }
      }
      if (traced) break;
    }
    if (!traced) {

      E -=  E_ExtLoop(type, (i>1) ? SS1[i-1] : -1, (j<n4) ? SS2[j+1] : -1, P);
      /**
      ***      if (i>1) E -= P->dangle5[type][SS1[i-1]]+extension_cost;
      ***      if (j<n4) E -= P->dangle3[type][SS2[j+1]]+extension_cost;
      ***      if (type>2) E -= P->TerminalAU;
      **/
      if (E != P->DuplexInit+2*extension_cost) {
        vrna_message_error("backtrack failed in fold duplex");
      } else break;
    }
  }
  if (i>1)  i--;
  if (j<n4) j++;

  struc = (char *) vrna_alloc(i0-i+1+j-j0+1+2);
  for (k=MAX2(i,1); k<=i0; k++) if (!st1[k-1]) st1[k-1] = '.';
  for (k=j0; k<=j; k++) if (!st2[k-1]) st2[k-1] = '.';
  strcpy(struc, st1+MAX2(i-1,0));
  strcat(struc, "&");
  strcat(struc, st2+j0-1);
  /* printf("%s %3d,%-3d : %3d,%-3d\n", struc, i,i0,j0,j);  */
  free(st1); free(st2);
  return struc;
}

PRIVATE duplexT fduplexfold(const char *s1, const char *s2, const int extension_cost,const int il_a, const int il_b, const int b_a, const int b_b) {
  int i, j, Emin, i_min, j_min,l1;
  duplexT mfe;
  char *struc;
  int bopen=b_b;
  int bext=b_a+extension_cost;
  int iopen=il_b;
  int iext_s=2*(il_a+extension_cost);/* iext_s 2 nt nucleotide extension of interior loop, on i and j side */
  int iext_ass=50+il_a+extension_cost;/* iext_ass assymetric extension of interior loop, either on i or on j side. */
  int min_colonne=INF; /* enthaelt das maximum einer kolonne */
  int i_length;
  int max_pos;/* get position of the best hit */
  int max_pos_j;
  int temp=INF;
  int min_j_colonne;
  int max=INF;
  vrna_md_t md;
  /* FOLLOWING NEXT 4 LINE DEFINES AN ARRAY CONTAINING POSITION OF THE SUBOPT IN S1 */

  n3 = (int) strlen(s1);
  n4 = (int) strlen(s2);
  /* delta_check is the minimal distance allowed for two hits to be accepted */
  /* if both hits are closer, reject the smaller ( in term of position)  hits  */
  /* i want to implement a function that, given a position in a long sequence and a small sequence, */
  /* duplexfold them at this position and report the result at the command line */
  /* for this i first need to rewrite backtrack in order to remove the printf functio */
  /* END OF DEFINITION FOR NEEDED SUBOPT DATA  */
  set_model_details(&md);
  if ((!P) || (fabs(P->temperature - temperature)>1e-6)) {
    update_fold_params();
    if(P)
      free(P);
    P = vrna_params(&md);
    make_pair_matrix();
  }
  /*local c array initialization---------------------------------------------*/
  c  = (int**)  vrna_alloc(sizeof(int *) * (n3+1));
  in = (int**) vrna_alloc(sizeof(int *) * (n3+1));
  bx = (int**) vrna_alloc(sizeof(int *) * (n3+1));
  by = (int**) vrna_alloc(sizeof(int *) * (n3+1));
  inx= (int**) vrna_alloc(sizeof(int *) * (n3+1));
  iny= (int**) vrna_alloc(sizeof(int *) * (n3+1));
  for (i=0; i<=n3; i++){
    c[i]  = (int *) vrna_alloc(sizeof(int) * (n4+1));
    in[i] = (int *) vrna_alloc(sizeof(int) * (n4+1));
    bx[i] = (int *) vrna_alloc(sizeof(int) * (n4+1));
    by[i] = (int *) vrna_alloc(sizeof(int) * (n4+1));
    inx[i] = (int *) vrna_alloc(sizeof(int) * (n4+1));
    iny[i] = (int *) vrna_alloc(sizeof(int) * (n4+1));
  }
  /*-------------------------------------------------------------------------*/
  /*end of array initialisation----------------------------------*/
  /*maybe int *** would be better*/
  encode_seqs(s1,s2);
  /* ------------------------------------------matrix initialisierung */
  for(i=0; i<n3; i++){
    for(j=0; j<n4; j++){
      in[i][j]=INF;/* no in before  1 */
      c[i][j] =INF; /* no bulge and no in before n2 */
      bx[i][j]=INF;/* no bulge before 1 */
      by[i][j]=INF;
      inx[i][j]=INF;/* no bulge before 1 */
      iny[i][j]=INF;
    }
  }

  /*--------------------------------------------------------local array*/


  /* -------------------------------------------------------------matrix initialisierung */
  i=11;
  i_length=n3-9;
  while(i < i_length) {
    j=n4-9;
    min_colonne=INF;
    while (10 < --j) {
      int type, type2;
      type = pair[S1[i]][S2[j]];
      /**
      *** Start duplex
      **/
      c[i][j]=type ? P->DuplexInit + 2*extension_cost : INF;
      /**
      *** update lin bx by linx liny matrix
      **/
      type2=pair[S2[j+1]][S1[i-1]];
      /**
      *** start/extend interior loop
      **/
      in[i][j]=MIN2(c[i - 1][j+1]+P->mismatchI[type2][SS2[j]][SS1[i]]+iopen+iext_s, in[i - 1][j]+iext_ass);
      /**
      *** start/extend nx1 target
      *** use same type2 as for in
      **/
      inx[i][j]=MIN2(c[i-1][j+1]+P->mismatch1nI[type2][SS2[j]][SS1[i]]+iopen+iext_s,
                     inx[i-1][j]+iext_ass);
      /**
      *** start/extend 1xn target
      *** use same type2 as for in
      **/
      iny[i][j]=MIN2(c[i-1][j+1]+P->mismatch1nI[type2][SS2[j]][SS1[i]]+iopen+iext_s,
                     iny[i][j+1]+iext_ass);
      /**
      *** extend interior loop
      **/
      in[i][j]=MIN2(in[i][j],in[i][j+1]+iext_ass);
      in[i][j]=MIN2(in[i][j],in[i - 1][j+1]+iext_s);
      /**
      *** start/extend bulge target
      **/
      type2=pair[S2[j]][S1[i-1]];
      bx[i][j]=MIN2(bx[i - 1][j]+bext, c[i - 1][j]+bopen+bext+(type2>2?P->TerminalAU:0));
      /**
      *** start/extend bulge query
      **/
      type2=pair[S2[j+1]][S1[i]];
      by[i][j]=MIN2(by[i][j+1]+bext, c[i][j+1]+bopen+bext+(type2>2?P->TerminalAU:0));
      /**
      ***end update recursion
      ***######################## Start stack extension##############################
      **/
      if(!type){continue;}
      c[i][j]+=E_ExtLoop(type, SS1[i-1] , SS2[j+1], P) + 2*extension_cost;
      /**
      *** stack extension
      **/
      if((type2=pair[S1[i-1]][S2[j+1]]))
        c[i][j]=MIN2(c[i - 1][j+1]+P->stack[rtype[type]][type2]+2*extension_cost, c[i][j]);
      /**
      *** 1x0 / 0x1 stack extension
      **/
      type2=pair[S1[i-1]][S2[j+2]];
      c[i][j]=MIN2(c[i - 1][j+2]+P->bulge[1]+P->stack[rtype[type]][type2]+3*extension_cost,c[i][j]);
      type2=pair[S1[i-2]][S2[j+1]];
      c[i][j]=MIN2(c[i - 2][j+1]+P->bulge[1]+P->stack[type2][rtype[type]]+3*extension_cost,c[i][j]);
      /**
      *** 1x1 / 2x2 stack extension
      **/
      type2=pair[S1[i-2]][S2[j+2]];
      c[i][j]=MIN2(c[i - 2][j+2]+P->int11[type2][rtype[type]][SS1[i-1]][SS2[j+1]]+4*extension_cost, c[i][j]);
      type2 = pair[S1[i-3]][S2[j+3]];
      c[i][j]=MIN2(c[i - 3][j+3]+P->int22[type2][rtype[type]][SS1[i-2]][SS1[i-1]][SS2[j+1]][SS2[j+2]]+6*extension_cost,c[i][j]);
      /**
      *** 1x2 / 2x1 stack extension
      *** E_IntLoop(1,2,type2, rtype[type],SS1[i-1], SS2[j+2], SS1[i-1], SS2[j+1], P) corresponds to
      *** P->int21[rtype[type]][type2][SS2[j+2]][SS1[i-1]][SS1[i-1]]
      **/
      type2 = pair[S1[i-3]][S2[j+2]];
      c[i][j]=MIN2(c[i - 3][j+2]+P->int21[rtype[type]][type2][SS2[j+1]][SS1[i-2]][SS1[i-1]]+5*extension_cost, c[i][j]);
      type2 = pair[S1[i-2]][S2[j+3]];
      c[i][j]=MIN2(c[i - 2][j+3]+P->int21[type2][rtype[type]][SS1[i-1]][SS2[j+1]][SS2[j+2]]+5*extension_cost, c[i][j]);

      /**
      *** 2x3 / 3x2 stack extension
      **/
      if((type2 = pair[S1[i-4]][S2[j+3]]))
        c[i][j]=MIN2(c[i - 4][j+3]+P->internal_loop[5]+P->ninio[2]+
                     P->mismatch23I[type2][SS1[i-3]][SS2[j+2]]+P->mismatch23I[rtype[type]][SS2[j+1]][SS1[i-1]]+7*extension_cost, c[i][j]);
      if((type2 = pair[S1[i-3]][S2[j+4]]))
        c[i][j]=MIN2(c[i - 3][j+4]+P->internal_loop[5]+P->ninio[2]+
                     P->mismatch23I[type2][SS1[i-2]][SS2[j+3]]+P->mismatch23I[rtype[type]][SS2[j+1]][SS1[i-1]]+7*extension_cost, c[i][j]);
      /**
      *** So now we have to handle 1x3, 3x1, 3x3, and mxn m,n > 3
      **/
      /**
      *** 3x3 or more
      **/
      c[i][j]=MIN2(in[i - 3][j+3]+P->mismatchI[rtype[type]][SS1[i-1]][SS2[j+1]]+2*iext_s+2*extension_cost, c[i][j]);
      /**
      *** 2xn or more
      **/
      c[i][j]=MIN2(in[i - 4][j+2]+P->mismatchI[rtype[type]][SS1[i-1]][SS2[j+1]]+iext_s+2*iext_ass+2*extension_cost, c[i][j]);
      /**
      *** nx2 or more
      **/
      c[i][j]=MIN2(in[i - 2][j+4]+P->mismatchI[rtype[type]][SS1[i-1]][SS2[j+1]]+iext_s+2*iext_ass+2*extension_cost, c[i][j]);
      /**
      *** nx1 n>2
      **/
      c[i][j]=MIN2(inx[i - 3][j+1]+P->mismatch1nI[rtype[type]][SS1[i-1]][SS2[j+1]]+iext_ass+iext_ass+2*extension_cost, c[i][j]);
      /**
      *** 1xn n>2
      **/
      c[i][j]=MIN2(iny[i - 1][j+3]+P->mismatch1nI[rtype[type]][SS1[i-1]][SS2[j+1]]+iext_ass+iext_ass+2*extension_cost, c[i][j]);
      /**
      *** nx0 n>1
      **/
      int bAU;
      bAU=(type>2?P->TerminalAU:0);
      c[i][j]=MIN2(bx[i - 2][j+1]+2*extension_cost+bext+bAU, c[i][j]);
      /**
      *** 0xn n>1
      **/
      c[i][j]=MIN2(by[i - 1][j+2]+2*extension_cost+bext+bAU, c[i][j]);
      temp=min_colonne;
      min_colonne=MIN2(c[i][j]+E_ExtLoop(rtype[type], SS2[j-1] , SS1[i+1] , P) + 2*extension_cost, min_colonne);
      if(temp>min_colonne){
        min_j_colonne=j;
      }
      /* ---------------------------------------------------------------------end update */
    }
    if(max>=min_colonne){
      max=min_colonne;
      max_pos=i;
      max_pos_j=min_j_colonne;
    }
    i++;
  }
  Emin=max;
  i_min=max_pos;
  j_min=max_pos_j;
  int dGe;
  dGe=0;
  struc = fbacktrack(i_min, j_min, extension_cost, il_a, il_b, b_a, b_b,&dGe);
  if (i_min<n3-10) i_min++;
  if (j_min>11 ) j_min--;
  l1 = strchr(struc, '&')-struc;
  int size;
  size=strlen(struc)-1;
  Emin-= size * (extension_cost);
  mfe.i = i_min;
  mfe.j = j_min;
  mfe.energy = (double) Emin/100.;
  mfe.energy_backtrack = (double) dGe/100.;
  mfe.structure = struc;
  free(S1); free(S2); free(SS1); free(SS2);
  for (i=0; i<=n3; i++) {
    free(c[i]);
    free(in[i]);
    free(bx[i]);
    free(by[i]);
    free(inx[i]);
    free(iny[i]);
  }
  free(c);free(in);free(bx);free(by);free(inx);free(iny);
  return mfe;
}


PRIVATE char *fbacktrack(int i, int j, const int extension_cost,const int il_a, const int il_b, const int b_a, const int b_b, int *dG) {
  /* backtrack structure going backwards from i, and forwards from j
     return structure in bracket notation with & as separator */
  int k, l, type, type2, E, traced, i0, j0;
  char *st1, *st2, *struc;
  int bopen=b_b;
  int bext=b_a+extension_cost;
  int iopen=il_b;
  int iext_s=2*(il_a+extension_cost);/* iext_s 2 nt nucleotide extension of interior loop, on i and j side */
  int iext_ass=50+il_a+extension_cost;/* iext_ass assymetric extension of interior loop, either on i or on j side. */
  st1 = (char *) vrna_alloc(sizeof(char)*(n3+1));
  st2 = (char *) vrna_alloc(sizeof(char)*(n4+1));
  i0=MIN2(i+1,n3-10); j0=MAX2(j-1,11);
  int state;
  state=1; /* we start backtracking from a a pair , i.e. c-matrix */
  /* state 1 -> base pair, c
     state 2 -> interior loop, in
     state 3 -> bx loop, bx
     state 4 -> by loop, by
  */
  traced=1;
  k=i;l=j;
  type=pair[S1[i]][S2[j]];
  *dG+=E_ExtLoop(rtype[type], SS2[j-1] , SS1[i+1] , P);
  /*     (type>2?P->TerminalAU:0)+P->dangle3[rtype[type]][SS1[i+1]]+P->dangle5[rtype[type]][SS2[j-1]]; */
  while (i>10 && j<=n4-9 && traced) {
    traced=0;
    switch(state){
    case 1:
      type = pair[S1[i]][S2[j]];
      int bAU;
      bAU=(type>2?P->TerminalAU:0);
      if(!type) vrna_message_error("backtrack failed in fold duplex");
      type2=pair[S1[i-1]][S2[j+1]];
      if(type2 && c[i][j]== (c[i - 1][j+1]+P->stack[rtype[type]][type2]+2*extension_cost)){
        k=i-1;
        l=j+1;
        (*dG)+=E_IntLoop(i-k-1, l-j-1, type2, rtype[type],SS1[k+1], SS2[l-1], SS1[i-1], SS2[j+1],P);
        st1[i-1] = '(';
        st2[j-1] = ')';
        i=k;
        j=l;
        state=1;
        traced=1;
        break;
      }
      type2=pair[S1[i-1]][S2[j+2]];
      if(type2 && c[i][j]==(c[i - 1][j+2]+P->bulge[1]+P->stack[rtype[type]][type2]+3*extension_cost)){
        k=i-1;
        l=j+2;
        *dG+=E_IntLoop(i-k-1, l-j-1, type2, rtype[type],SS1[k+1], SS2[l-1], SS1[i-1], SS2[j+1],P);
        st1[i-1] = '(';
        st2[j-1] = ')';
        i=k;
        j=l;
        state=1;
        traced=1;
        break;
      }
      type2=pair[S1[i-2]][S2[j+1]];
      if(type2 && c[i][j]==(c[i - 2][j+1]+P->bulge[1]+P->stack[type2][rtype[type]]+3*extension_cost)){
        k=i-2;
        l=j+1;
        *dG+=E_IntLoop(i-k-1, l-j-1, type2, rtype[type],SS1[k+1], SS2[l-1], SS1[i-1], SS2[j+1],P);
        st1[i-1] = '(';
        st2[j-1] = ')';
        i=k;
        j=l;
        state=1;
        traced=1;
        break;
      }
      type2=pair[S1[i-2]][S2[j+2]];
      if(type2 && c[i][j]==(c[i - 2][j+2]+P->int11[type2][rtype[type]][SS1[i-1]][SS2[j+1]]+4*extension_cost)){
        k=i-2;
        l=j+2;
        *dG+=E_IntLoop(i-k-1, l-j-1, type2, rtype[type],SS1[k+1], SS2[l-1], SS1[i-1], SS2[j+1],P);
        st1[i-1] = '(';
        st2[j-1] = ')';
        i=k;
        j=l;
        state=1;
        traced=1;
        break;
      }
      type2 = pair[S1[i-3]][S2[j+3]];
      if(type2 && c[i][j]==(c[i - 3][j+3]+P->int22[type2][rtype[type]][SS1[i-2]][SS1[i-1]][SS2[j+1]][SS2[j+2]]+6*extension_cost)){
        k=i-3;
        l=j+3;
        *dG+=E_IntLoop(i-k-1, l-j-1, type2, rtype[type],SS1[k+1], SS2[l-1], SS1[i-1], SS2[j+1],P);
        st1[i-1] = '(';
        st2[j-1] = ')';
        i=k;
        j=l;
        state=1;
        traced=1;
        break;
      }
      type2 = pair[S1[i-3]][S2[j+2]];
      if(type2 && c[i][j]==(c[i - 3][j+2]+P->int21[rtype[type]][type2][SS2[j+1]][SS1[i-2]][SS1[i-1]]+5*extension_cost)){
        k=i-3;
        l=j+2;
        *dG+=E_IntLoop(i-k-1, l-j-1, type2, rtype[type],SS1[k+1], SS2[l-1], SS1[i-1], SS2[j+1],P);
        st1[i-1] = '(';
        st2[j-1] = ')';
        i=k;
        j=l;
        state=1;
        traced=1;
        break;
      }
      type2 = pair[S1[i-2]][S2[j+3]];
      if(type2 && c[i][j]==(c[i - 2][j+3]+P->int21[type2][rtype[type]][SS1[i-1]][SS2[j+1]][SS2[j+2]]+5*extension_cost)){
        k=i-2;
        l=j+3;
        *dG+=E_IntLoop(i-k-1, l-j-1, type2, rtype[type],SS1[k+1], SS2[l-1], SS1[i-1], SS2[j+1],P);
        st1[i-1] = '(';
        st2[j-1] = ')';
        i=k;
        j=l;
        state=1;
        traced=1;
        break;
      }
      type2 = pair[S1[i-4]][S2[j+3]];
      if(type2 && c[i][j]==(c[i - 4][j+3]+P->internal_loop[5]+P->ninio[2]+
                            P->mismatch23I[type2][SS1[i-3]][SS2[j+2]]+P->mismatch23I[rtype[type]][SS2[j+1]][SS1[i-1]]+7*extension_cost)){
        k=i-4;
        l=j+3;
        *dG+=E_IntLoop(i-k-1, l-j-1, type2, rtype[type],SS1[k+1], SS2[l-1], SS1[i-1], SS2[j+1],P);
        st1[i-1] = '(';
        st2[j-1] = ')';
        i=k;
        j=l;
        state=1;
        traced=1;
        break;
      }
      type2 = pair[S1[i-3]][S2[j+4]];
      if(type2 && c[i][j]==(c[i - 3][j+4]+P->internal_loop[5]+P->ninio[2]+
                            P->mismatch23I[type2][SS1[i-2]][SS2[j+3]]+P->mismatch23I[rtype[type]][SS2[j+1]][SS1[i-1]]+7*extension_cost)){
        k=i-3;
        l=j+4;
        *dG+=E_IntLoop(i-k-1, l-j-1, type2, rtype[type],SS1[k+1], SS2[l-1], SS1[i-1], SS2[j+1],P);
        st1[i-1] = '(';
        st2[j-1] = ')';
        i=k;
        j=l;
        state=1;
        traced=1;
        break;
      }
      if(c[i][j]==(in[i - 3][j+3]+P->mismatchI[rtype[type]][SS1[i-1]][SS2[j+1]]+2*extension_cost+2*iext_s)){
        k=i;
        l=j;
        st1[i-1] = '(';
        st2[j-1] = ')';
        i=i-3;
        j=j+3;
        state=2;
        traced=1;
        break;
      }
      if(c[i][j]==(in[i - 4][j+2]+P->mismatchI[rtype[type]][SS1[i-1]][SS2[j+1]]+iext_s+2*iext_ass+2*extension_cost)){
        k=i;
        l=j;
        st1[i-1] = '(';
        st2[j-1] = ')';
        i=i-4;
        j=j+2;
        state=2;
        traced=1;
        break;
      }
      if(c[i][j]==(in[i - 2][j+4]+P->mismatchI[rtype[type]][SS1[i-1]][SS2[j+1]]+iext_s+2*iext_ass+2*extension_cost)){
        k=i;
        l=j;
        st1[i-1] = '(';
        st2[j-1] = ')';
        i=i-2;
        j=j+4;
        state=2;
        traced=1;
        break;
      }
      if(c[i][j]==(inx[i - 3][j+1]+P->mismatch1nI[rtype[type]][SS1[i-1]][SS2[j+1]]+iext_ass+iext_ass+2*extension_cost)){
        k=i;
        l=j;
        st1[i-1] = '(';
        st2[j-1] = ')';
        i=i-3;
        j=j+1;
        state=5;
        traced=1;
        break;
      }
      if(c[i][j]==(iny[i - 1][j+3]+P->mismatch1nI[rtype[type]][SS1[i-1]][SS2[j+1]]+iext_ass+iext_ass+2*extension_cost)){
        k=i;
        l=j;
        st1[i-1] = '(';
        st2[j-1] = ')';
        i=i-1;
        j=j+3;
        state=6;
        traced=1;
        break;
      }
      if(c[i][j]==(bx[i - 2][j+1]+2*extension_cost+bext+bAU)){
        k=i;
        l=j;
        st1[i-1] = '(';
        st2[j-1] = ')';
        i=i-2;
        j=j+1;
        state=3;
        traced=1;
        break;
      }
      if(c[i][j]==(by[i - 1][j+2]+2*extension_cost+bext+bAU)){
        k=i;
        l=j;
        st1[i-1] = '(';
        st2[j-1] = ')';
        i=i-1;
        j=j+2;
        state=4;
        traced=1;
        break;
      }
      break;
    case 2:
      if(in[i][j]==(in[i - 1][j+1]+iext_s)){
        i--;
        j++;
        state=2;
        traced=1;
        break;
      }
      if(in[i][j]==(in[i - 1][j]+iext_ass)){
        i=i-1;
        state=2;
        traced=1;
        break;
      }
      if(in[i][j]==(in[i][j+1]+iext_ass)){
        j++;
        state=2;
        traced=1;
        break;
      }
      type2=pair[S2[j+1]][S1[i-1]];
      if(type2 && in[i][j]==(c[i - 1][j+1]+P->mismatchI[type2][SS2[j]][SS1[i]]+iopen+iext_s)){
        int temp; temp=k; k=i-1; i=temp;
        temp=l; l=j+1; j=temp;
        type=pair[S1[i]][S2[j]];
        *dG+=E_IntLoop(i-k-1, l-j-1, type2, rtype[type],SS1[k+1], SS2[l-1], SS1[i-1], SS2[j+1],P);
        i=k;
        j=l;
        state=1;
        traced=1;
        break;
      }
    case 3:
      if(bx[i][j]==(bx[i - 1][j]+bext)){
        i--;
        state=3;
        traced=1;
        break;
      }
      type2=pair[S2[j]][S1[i-1]];
      if(type2 && bx[i][j]==(c[i - 1][j]+bopen+bext+(type2>2?P->TerminalAU:0))){
        int temp; temp=k; k=i-1; i=temp;
        temp=l; l=j; j=temp;
        type=pair[S1[i]][S2[j]];
        *dG+=E_IntLoop(i-k-1, l-j-1, type2, rtype[type],SS1[k+1], SS2[l-1], SS1[i-1], SS2[j+1],P);
        i=k;
        j=l;
        state=1;
        traced=1;
        break;
      }
    case 4:
      if(by[i][j]==(by[i][j+1] + bext)){
        j++;

        state=4;
        traced=1;
        break;
      }
      type2=pair[S2[j+1]][S1[i]];
      if(type2 && by[i][j]==(c[i][j+1]+bopen+bext+(type2>2?P->TerminalAU:0))){
        int temp; temp=k; k=i; i=temp;
        temp=l; l=j+1; j=temp;
        type=pair[S1[i]][S2[j]];
        *dG+=E_IntLoop(i-k-1, l-j-1, type2, rtype[type],SS1[k+1], SS2[l-1], SS1[i-1], SS2[j+1],P);
        i=k;
        j=l;
        state=1;
        traced=1;
        break;
      }
    case 5:
      if(inx[i][j]==(inx[i-1][j]+iext_ass)) {
        i--;
        state=5;
        traced=1;
        break;
      }
      type2=pair[S2[j+1]][S1[i-1]];
      if(type2 && inx[i][j]==(c[i-1][j+1]+P->mismatch1nI[type2][SS2[j]][SS1[i]]+iopen+iext_s)){
        int temp; temp=k; k=i-1; i=temp;
        temp=l; l=j+1; j=temp;
        type=pair[S1[i]][S2[j]];
        *dG+=E_IntLoop(i-k-1, l-j-1, type2, rtype[type],SS1[k+1], SS2[l-1], SS1[i-1], SS2[j+1],P);
        i=k;
        j=l;
        state=1;
        traced=1;
        break;
      }
    case 6:
      if(iny[i][j]==(iny[i][j+1]+iext_ass)) {
        j++;
        state=6;
        traced=1;
        break;
      }
      type2=pair[S2[j+1]][S1[i-1]];
      if(type2 && iny[i][j]==(c[i-1][j+1]+P->mismatch1nI[type2][SS2[j]][SS1[i]]+iopen+iext_s)){
        int temp; temp=k; k=i-1; i=temp;
        temp=l; l=j+1; j=temp;
        type=pair[S1[i]][S2[j]];
        *dG+=E_IntLoop(i-k-1, l-j-1, type2, rtype[type],SS1[k+1], SS2[l-1], SS1[i-1], SS2[j+1],P);
        i=k;
        j=l;
        state=1;
        traced=1;
        break;
      }
    }
  }
  if (!traced) {
    E=c[i][j];
    /**
    ***    if (i>1) {E -= P->dangle5[type][SS1[i-1]]+extension_cost; *dG+=P->dangle5[type][SS1[i-1]];}
    ***    if (j<n4){E -= P->dangle3[type][SS2[j+1]]+extension_cost; *dG+=P->dangle3[type][SS2[j+1]];}
    ***    if (type>2) {E -= P->TerminalAU; *dG+=P->TerminalAU;}
    **/
    int correction;
    correction = E_ExtLoop(type, (i>1) ? SS1[i-1] : -1, (j<n4) ? SS2[j+1] : -1, P);
    *dG+=correction;
    E-=correction+2*extension_cost;
    if (E != P->DuplexInit+2*extension_cost) {
      vrna_message_error("backtrack failed in second fold duplex");
    }
    else{
      *dG+=P->DuplexInit;
      st1[i-1]='(';
      st2[j-1]=')';
    }
  }
  if (i>11)  i--;
  if (j<n4-10) j++;
  struc = (char *) vrna_alloc(i0-i+1+j-j0+1+2);
  for (k=MAX2(i,1); k<=i0; k++) if (!st1[k-1]) st1[k-1] = '.';
  for (k=j0; k<=j; k++) if (!st2[k-1]) st2[k-1] = '.';
  strcpy(struc, st1+MAX2(i-1,0));
  strcat(struc, "&");
  strcat(struc, st2+j0-1);
  /* printf("%s %3d,%-3d : %3d,%-3d\n", struc, i,i0,j0,j);  */
  free(st1); free(st2);
  return struc;
}


duplexT ** Lduplexfold(const char *s1, const char *s2, const int threshold, const int extension_cost, const int alignment_length, const int delta, const int fast, const int il_a, const int il_b, const int b_a, const int b_b)
{
  /**
  *** See variable definition in fduplexfold_XS
  **/
  int i, j;
  int bopen=b_b;
  int bext=b_a+extension_cost;
  int iopen=il_b;
  int iext_s=2*(il_a+extension_cost);/* iext_s 2 nt nucleotide extension of interior loop, on i and j side */
  int iext_ass=50+il_a+extension_cost;/* iext_ass assymetric extension of interior loop, either on i or on j side. */
  int min_colonne=INF; /* enthaelt das maximum einer kolonne */
  int i_length;
  int max_pos;/* get position of the best hit */
  int max_pos_j;
  int temp=INF;
  int min_j_colonne;
  int max=INF;
  int *position; /* contains the position of the hits with energy > E */
  int *position_j;
  /**
  *** 1D array corresponding to the standard 2d recursion matrix
  *** Makes the computation 20% faster
  **/
  int *SA;
  vrna_md_t md;
  
  /**
  *** variable initialization
  **/
  n1 = (int) strlen(s1);
  n2 = (int) strlen(s2);
  /**
  *** Sequence encoding
  **/
  set_model_details(&md);
  if ((!P) || (fabs(P->temperature - temperature)>1e-6)) {
    update_fold_params();
    if(P)
      free(P);
    P = vrna_params(&md);
    make_pair_matrix();
  }
  encode_seqs(s1,s2);
  /**
  *** Position of the high score on the target and query sequence
  **/
  position = (int *) vrna_alloc((delta+n1+3+delta) * sizeof(int));
  position_j= (int *) vrna_alloc((delta+n1+3+delta) * sizeof(int));
  /**
  *** instead of having 4 2-dim arrays we use a unique 1-dim array
  *** The mapping 2d -> 1D is done based ont the macro
  *** LCI(i,j,l)      ((i     )*l + j)
  *** LINI(i,j,l)     ((i +  5)*l + j)
  *** LBXI(i,j,l)     ((i + 10)*l + j)
  *** LBYI(i,j,l)     ((i + 15)*l + j)
  *** LINIX(i,j,l)    ((i + 20)*l + j)
  *** LINIY(i,j,l)    ((i + 25)*l + j)
  ***
  *** SA has a length of 5 (number of columns we look back) *
  ***                  * 6 (number of structures we look at) *
  ***                  * length of the sequence
  **/
  SA=(int *) vrna_alloc(sizeof(int)*5*6*(n2+5));
  for(j=n2+4;j>=0;j--) {
    SA[(j*30)   ]=SA[(j*30)+1   ]=SA[(j*30)+2   ]=SA[(j*30)+3   ]=SA[(j*30)+4   ]=INF;
    SA[(j*30)+5 ]=SA[(j*30)+1+5 ]=SA[(j*30)+2+5 ]=SA[(j*30)+3+5 ]=SA[(j*30)+4+5 ]=INF;
    SA[(j*30)+10]=SA[(j*30)+1+10]=SA[(j*30)+2+10]=SA[(j*30)+3+10]=SA[(j*30)+4+10]=INF;
    SA[(j*30)+15]=SA[(j*30)+1+15]=SA[(j*30)+2+15]=SA[(j*30)+3+15]=SA[(j*30)+4+15]=INF;
    SA[(j*30)+20]=SA[(j*30)+1+20]=SA[(j*30)+2+20]=SA[(j*30)+3+20]=SA[(j*30)+4+20]=INF;
    SA[(j*30)+25]=SA[(j*30)+1+25]=SA[(j*30)+2+25]=SA[(j*30)+3+25]=SA[(j*30)+4+25]=INF;
  }
  i=10;
  i_length= n1 - 9 ;
  while(i < i_length) {
    int idx=i%5;
    int idx_1=(i-1)%5;
    int idx_2=(i-2)%5;
    int idx_3=(i-3)%5;
    int idx_4=(i-4)%5;
    j=n2-9;
    while (9 < --j) {
      int type, type2;
      type = pair[S1[i]][S2[j]];
      /**
      *** Start duplex
      **/
      SA[LCI(idx,j,n2)]=type ? P->DuplexInit + 2*extension_cost : INF;
      /**
      *** update lin bx by linx liny matrix
      **/
      type2=pair[S2[j+1]][S1[i-1]];
      /**
      *** start/extend interior loop
      **/
      SA[LINI(idx,j,n2)]=MIN2(SA[LCI(idx_1,j+1,n2)]+P->mismatchI[type2][SS2[j]][SS1[i]]+iopen+iext_s,
                              SA[LINI(idx_1,j,n2)]+iext_ass);
      /**
      *** start/extend nx1 target
      *** use same type2 as for in
      **/
      SA[LINIX(idx,j,n2)]=MIN2(SA[LCI(idx_1,j+1,n2)]+P->mismatch1nI[type2][SS2[j]][SS1[i]]+iopen+iext_s,
                               SA[LINIX(idx_1,j,n2)]+iext_ass);
      /**
      *** start/extend 1xn target
      *** use same type2 as for in
      **/
      SA[LINIY(idx,j,n2)]=MIN2(SA[LCI(idx_1,j+1,n2)]+P->mismatch1nI[type2][SS2[j]][SS1[i]]+iopen+iext_s,
                               SA[LINIY(idx,j+1,n2)]+iext_ass);
      /**
      *** extend interior loop
      **/
      SA[LINI(idx,j,n2)]=MIN2(SA[LINI(idx,j,n2)],SA[LINI(idx,j+1,n2)]+iext_ass);
      SA[LINI(idx,j,n2)]=MIN2(SA[LINI(idx,j,n2)],SA[LINI(idx_1,j+1,n2)]+iext_s);
      /**
      *** start/extend bulge target
      **/
      type2=pair[S2[j]][S1[i-1]];
      SA[LBXI(idx,j,n2)]=MIN2(SA[LBXI(idx_1,j,n2)]+bext, SA[LCI(idx_1,j,n2)]+bopen+bext+(type2>2?P->TerminalAU:0));
      /**
      *** start/extend bulge query
      **/
      type2=pair[S2[j+1]][S1[i]];
      SA[LBYI(idx,j,n2)]=MIN2(SA[LBYI(idx,j+1,n2)]+bext, SA[LCI(idx,j+1,n2)]+bopen+bext+(type2>2?P->TerminalAU:0));
      /**
      ***end update recursion
      ***##################### Start stack extension ######################
      **/
      if(!type){continue;}
      /**
      *** stack extension
      **/
      SA[LCI(idx,j,n2)]+= E_ExtLoop(type, SS1[i-1] , SS2[j+1], P) + 2*extension_cost;
      /**
      *** stack extension
      **/
      if((type2=pair[S1[i-1]][S2[j+1]]))
        SA[LCI(idx,j,n2)]=MIN2(SA[LCI(idx_1,j+1,n2)]+P->stack[rtype[type]][type2]+2*extension_cost, SA[LCI(idx,j,n2)]);
      /**
      *** 1x0 / 0x1 stack extension
      **/
      if((type2=pair[S1[i-1]][S2[j+2]]))
        SA[LCI(idx,j,n2)]=MIN2(SA[LCI(idx_1,j+2,n2)]+P->bulge[1]+P->stack[rtype[type]][type2]+3*extension_cost,SA[LCI(idx,j,n2)]);
      if((type2=pair[S1[i-2]][S2[j+1]]))
        SA[LCI(idx,j,n2)]=MIN2(SA[LCI(idx_2,j+1,n2)]+P->bulge[1]+P->stack[type2][rtype[type]]+3*extension_cost,SA[LCI(idx,j,n2)]);
      /**
      *** 1x1 / 2x2 stack extension
      **/
      if((type2=pair[S1[i-2]][S2[j+2]]))
        SA[LCI(idx,j,n2)]=MIN2(SA[LCI(idx_2,j+2,n2)]+P->int11[type2][rtype[type]][SS1[i-1]][SS2[j+1]]+4*extension_cost, SA[LCI(idx,j,n2)]);
      if((type2 = pair[S1[i-3]][S2[j+3]]))
        SA[LCI(idx,j,n2)]=MIN2(SA[LCI(idx_3,j+3,n2)]+P->int22[type2][rtype[type]][SS1[i-2]][SS1[i-1]][SS2[j+1]][SS2[j+2]]+6*extension_cost,SA[LCI(idx,j,n2)]);
      /**
      *** 1x2 / 2x1 stack extension
      *** E_IntLoop(1,2,type2, rtype[type],SS1[i-1], SS2[j+2], SS1[i-1], SS2[j+1], P) corresponds to
      *** P->int21[rtype[type]][type2][SS2[j+2]][SS1[i-1]][SS1[i-1]]
      **/
      if((type2 = pair[S1[i-3]][S2[j+2]]))
        SA[LCI(idx,j,n2)]=MIN2(SA[LCI(idx_3,j+2,n2)]+P->int21[rtype[type]][type2][SS2[j+1]][SS1[i-2]][SS1[i-1]]+5*extension_cost, SA[LCI(idx,j,n2)]);
      if((type2 = pair[S1[i-2]][S2[j+3]]))
        SA[LCI(idx,j,n2)]=MIN2(SA[LCI(idx_2,j+3,n2)]+P->int21[type2][rtype[type]][SS1[i-1]][SS2[j+1]][SS2[j+2]]+5*extension_cost, SA[LCI(idx,j,n2)]);
      /**
      *** 2x3 / 3x2 stack extension
      **/
      if((type2 = pair[S1[i-4]][S2[j+3]]))
        SA[LCI(idx,j,n2)]=MIN2(SA[LCI(idx_4,j+3,n2)]+P->internal_loop[5]+P->ninio[2]+
                               P->mismatch23I[type2][SS1[i-3]][SS2[j+2]]+P->mismatch23I[rtype[type]][SS2[j+1]][SS1[i-1]]+7*extension_cost, SA[LCI(idx,j,n2)]);
      if((type2 = pair[S1[i-3]][S2[j+4]]))
        SA[LCI(idx,j,n2)]=MIN2(SA[LCI(idx_3,j+4,n2)]+P->internal_loop[5]+P->ninio[2]+
                               P->mismatch23I[type2][SS1[i-2]][SS2[j+3]]+P->mismatch23I[rtype[type]][SS2[j+1]][SS1[i-1]]+7*extension_cost, SA[LCI(idx,j,n2)]);
      /**
      *** So now we have to handle 1x3, 3x1, 3x3, and mxn m,n > 3
      **/
      /**
      *** 3x3 or more
      **/
      SA[LCI(idx,j,n2)]=MIN2(SA[LINI(idx_3,j+3,n2)]+P->mismatchI[rtype[type]][SS1[i-1]][SS2[j+1]]+2*iext_s+2*extension_cost,SA[LCI(idx,j,n2)]);
      /**
      *** 2xn or more
      **/
      SA[LCI(idx,j,n2)]=MIN2(SA[LINI(idx_4,j+2,n2)]+P->mismatchI[rtype[type]][SS1[i-1]][SS2[j+1]]+iext_s+2*iext_ass+2*extension_cost, SA[LCI(idx,j,n2)]);
      /**
      *** nx2 or more
      **/
      SA[LCI(idx,j,n2)]=MIN2(SA[LINI(idx_2,j+4,n2)]+P->mismatchI[rtype[type]][SS1[i-1]][SS2[j+1]]+iext_s+2*iext_ass+2*extension_cost, SA[LCI(idx,j,n2)]);
      /**
      *** nx1 n>2
      **/
      SA[LCI(idx,j,n2)]=MIN2(SA[LINIX(idx_3,j+1,n2)]+P->mismatch1nI[rtype[type]][SS1[i-1]][SS2[j+1]]+iext_ass+iext_ass+2*extension_cost, SA[LCI(idx,j,n2)]);
      /**
      *** 1xn n>2
      **/
      SA[LCI(idx,j,n2)]=MIN2(SA[LINIY(idx_1,j+3,n2)]+P->mismatch1nI[rtype[type]][SS1[i-1]][SS2[j+1]]+iext_ass+iext_ass+2*extension_cost, SA[LCI(idx,j,n2)]);
      /**
      *** nx0 n>1
      **/
      int bAU;
      bAU=(type>2?P->TerminalAU:0);
      SA[LCI(idx,j,n2)]=MIN2(SA[LBXI(idx_2,j+1,n2)]+2*extension_cost+bext+bAU,SA[LCI(idx,j,n2)]);
      /**
      *** 0xn n>1
      **/
      SA[LCI(idx,j,n2)]=MIN2(SA[LBYI(idx_1,j+2,n2)]+2*extension_cost+bext+bAU,SA[LCI(idx,j,n2)]);
      temp=min_colonne;

      min_colonne=MIN2(SA[LCI(idx,j,n2)]+E_ExtLoop(rtype[type], SS2[j-1] , SS1[i+1] , P) + 2*extension_cost, min_colonne);
      if(temp>min_colonne){
        min_j_colonne=j;
      }
    }
    if(max>=min_colonne){
      max=min_colonne;
      max_pos=i;
      max_pos_j=min_j_colonne;
    }
    position[i+delta]=min_colonne;min_colonne=INF;
    position_j[i+delta]=min_j_colonne;
    i++;
  }
  /* printf("MAX: %d",max); */
  free(S1); free(S2); free(SS1); free(SS2);
  if(max<threshold){
    find_max(position, position_j, delta, threshold, alignment_length, s1, s2, extension_cost, fast, il_a, il_b, b_a, b_b);
  }
  if(max<INF){
    plot_max(max, max_pos, max_pos_j,alignment_length, s1, s2, extension_cost,fast, il_a, il_b, b_a, b_b);
  }
  free(SA);
  free(position);
  free(position_j);
  return NULL;
}




PRIVATE void find_max(const int *position, const int *position_j,const int delta, const int threshold, const int alignment_length, const char *s1, const char *s2, const int extension_cost, const int fast,const int il_a, const int il_b, const int b_a, const int b_b){
  int pos=n1-9;
  if(fast==1){
    while(10 < pos--){
      int temp_min=0;
      if(position[pos+delta]<(threshold)){
        int search_range;
        search_range=delta+1;
        while(--search_range){
          if(position[pos+delta-search_range]<=position[pos+delta-temp_min]){
            temp_min=search_range;
          }
        }
        pos-=temp_min;
        int max_pos_j;
        max_pos_j=position_j[pos+delta];
        int max;
        max=position[pos+delta];
        printf("target upper bound %d: query lower bound %d  (%5.2f) \n", pos-10, max_pos_j-10, ((double)max)/100);
        pos=MAX2(10,pos+temp_min-delta);
      }
    }
  }
  else if(fast==2){
    pos=n1-9;
    while(10 < pos--){
      int temp_min=0;
      if(position[pos+delta]<(threshold)){
        int search_range;
        search_range=delta+1;
        while(--search_range){
          if(position[pos+delta-search_range]<=position[pos+delta-temp_min]){
            temp_min=search_range;
          }
        }
        pos-=temp_min;
        int max_pos_j;
        max_pos_j=position_j[pos+delta];
        /* max_pos_j und pos entsprechen die realen position
           in der erweiterten sequenz.
           pos=1 -> position 1 in the sequence (and not 0 like in C)
           max_pos_j -> position 1 in the sequence ( not 0 like in C)
        */
        int alignment_length2; alignment_length2 = MIN2(n1,n2);
        int begin_t=MAX2(11, pos-alignment_length2+1);/* 10 */
        int end_t  =MIN2(n1-10, pos+1);
        int begin_q=MAX2(11, max_pos_j-1); /* 10 */
        int end_q  =MIN2(n2-10, max_pos_j+alignment_length2-1);
        char *s3 = (char*) vrna_alloc(sizeof(char)*(end_t - begin_t +2 + 20));
        char *s4 = (char*) vrna_alloc(sizeof(char)*(end_q - begin_q +2 + 20));
        strcpy(s3,"NNNNNNNNNN");strcpy(s4,"NNNNNNNNNN");
        strncat(s3, (s1+begin_t-1),  end_t - begin_t +1);
        strncat(s4, (s2+begin_q-1) , end_q - begin_q +1);
        strcat(s3,"NNNNNNNNNN");strcat(s4,"NNNNNNNNNN");
        s3[end_t -begin_t +1 +20 ]='\0';
        s4[end_q -begin_q +1 +20]='\0';
        duplexT test;
        test = fduplexfold(s3, s4, extension_cost,il_a, il_b, b_a, b_b);
        if(test.energy * 100 < threshold){
          int l1=strchr(test.structure, '&')-test.structure;
          printf("%s %3d,%-3d : %3d,%-3d (%5.2f) [%5.2f]  i:%d,j:%d <%5.2f>\n", test.structure,
                 begin_t-10+test.i-l1-10,
                 begin_t-10+test.i-1-10,
                 begin_q-10 + test.j-1-10 ,
                 (begin_q -11) + test.j + (int)strlen(test.structure)-l1-2-10,
                 test.energy,test.energy_backtrack, pos-10, max_pos_j-10, ((double) position[pos+delta])/100);
          pos=MAX2(10,pos+temp_min-delta);
        }
        free(s3);free(s4);
        free(test.structure);
      }
    }
  }
#if 0
  else if(fast==3){
    pos=n1-9;
    while(10 < pos--){
      int temp_min=0;
      if(position[pos+delta]<(threshold)){
        int search_range;
        search_range=delta+1;
        while(--search_range){
          if(position[pos+delta-search_range]<=position[pos+delta-temp_min]){
            temp_min=search_range;
          }
        }
        pos-=temp_min;
        int max_pos_j;
        max_pos_j=position_j[pos+delta];
        /* max_pos_j und pos entsprechen die realen position
           in der erweiterten sequenz.
           pos=1 -> position 1 in the sequence (and not 0 like in C)
           max_pos_j -> position 1 in the sequence ( not 0 like in C)
        */
        //Here we can start the reverse recursion for the
        //Starting from the reported pos / max_pos_j we start the recursion
        //We have to be careful with the fact that all energies are inverted.

        int alignment_length2;
        //Select the smallest interaction length in order to define the new interaction length
        alignment_length2 = MIN2(n1-pos + 1,max_pos_j - 1 + 1);
        //
        int begin_t=MAX2(11, pos-alignment_length2+1);/* 10 */
        int end_t  =MIN2(n1-10, pos+1);
        int begin_q=MAX2(11, max_pos_j-1); /* 10 */
        int end_q  =MIN2(n2-10, max_pos_j+alignment_length2-1);
        char *s3 = (char*) vrna_alloc(sizeof(char)*(end_t - begin_t +2 + 20));
        char *s4 = (char*) vrna_alloc(sizeof(char)*(end_q - begin_q +2 + 20));
        strcpy(s3,"NNNNNNNNNN");strcpy(s4,"NNNNNNNNNN");
        strncat(s3, (s1+begin_t-1),  end_t - begin_t +1);
        strncat(s4, (s2+begin_q-1) , end_q - begin_q +1);
        strcat(s3,"NNNNNNNNNN");strcat(s4,"NNNNNNNNNN");
        s3[end_t -begin_t +1 +20 ]='\0';
        s4[end_q -begin_q +1 +20]='\0';
        duplexT test;
        test = fduplexfold(s4, s3, extension_cost,il_a, il_b, b_a, b_b);
        if(test.energy * 100 < threshold){
          int structureLength=strlen(test.structure);
          int l1=strchr(test.structure, '&')-test.structure;
          int start_t,end_t,start_q,end_q;


          /*reverse structure string*/
          char *reverseStructure = (char*) vrna_alloc(sizeof(char)*(structureLength+1));
          int posStructure;
          for(posStructure=l1+1; posStructure < structureLength; posStructure++){
            if(test.structure[posStructure]==')'){
              reverseStructure[posStructure-l1-1] = '(';
            }
            else{
              reverseStructure[posStructure-l1-1] = test.structure[posStructure];
            }
          }
          reverseStructure[structureLength-1-l1]='&';
          for(posStructure=0; posStructure<l1; posStructure++){
            if(test.structure[posStructure]=='('){
              reverseStructure[structureLength+posStructure-l1] = ')';
            }
            else{
              reverseStructure[structureLength+posStructure-l1] = test.structure[posStructure];
            }
          }
          reverseStructure[structureLength]='\0';
          //          l1=strchr(reverse.structure, '&')-test.structure;


          printf("%s %3d,%-3d : %3d,%-3d (%5.2f) [%5.2f] i:%d,j:%d <%5.2f>\n", reverseStructure,
                 begin_t-10 + test.j-1-10,
                 (begin_t -11) + test.j + strlen(test.structure)-l1-2-10,
                 begin_q-10+test.i-l1-10,
                 begin_q-10+test.i-1-10,
                 test.energy,test.energy_backtrack,pos, max_pos_j, ((double) position[pos+delta])/100);
          pos=MAX2(10,pos+temp_min-delta);
        }
        free(s3);free(s4);
        free(test.structure);
      }
    }
  }
#endif
  else{
    pos=n1-9;
    while(10 < pos--){
      int temp_min=0;
      if(position[pos+delta]<(threshold)){
        int search_range;
        search_range=delta+1;
        while(--search_range){
          if(position[pos+delta-search_range]<=position[pos+delta-temp_min]){
            temp_min=search_range;
          }
        }
        pos-=temp_min;
        int max_pos_j;
        max_pos_j=position_j[pos+delta];
        /* max_pos_j und pos entsprechen die realen position
           in der erweiterten sequenz.
           pos=1 -> position 1 in the sequence (and not 0 like in C)
           max_pos_j -> position 1 in the sequence ( not 0 like in C)
        */
        int alignment_length2; alignment_length2 = MIN2(n1,n2);
        int begin_t=MAX2(11, pos-alignment_length2+1);/* 10 */
        int end_t  =MIN2(n1-10, pos+1);
        int begin_q=MAX2(11, max_pos_j-1); /* 10 */
        int end_q  =MIN2(n2-10, max_pos_j+alignment_length2-1);
        char *s3 = (char*) vrna_alloc(sizeof(char)*(end_t - begin_t +2));
        char *s4 = (char*) vrna_alloc(sizeof(char)*(end_q - begin_q +2));
        strncpy(s3, (s1+begin_t-1),  end_t - begin_t +1);
        strncpy(s4, (s2+begin_q-1) , end_q - begin_q +1);
        s3[end_t -begin_t +1 ]='\0';
        s4[end_q -begin_q +1 ]='\0';
        duplexT test;
        test = duplexfold(s3, s4, extension_cost);
        if(test.energy * 100 < threshold){
          int l1=strchr(test.structure, '&')-test.structure;
          printf("%s %3d,%-3d : %3d,%-3d (%5.2f)  i:%d,j:%d <%5.2f>\n", test.structure,
                 begin_t-10+test.i-l1,
                 begin_t-10+test.i-1,
                 begin_q-10 + test.j-1 ,
                 (begin_q -11) + test.j + (int)strlen(test.structure)-l1-2,
                 test.energy, pos-10, max_pos_j-10, ((double) position[pos+delta])/100);
          pos=MAX2(10,pos+temp_min-delta);
        }
        free(s3);free(s4);
        free(test.structure);
      }
    }
  }
}
PRIVATE void plot_max(const int max, const int max_pos, const int max_pos_j, const int alignment_length, const char *s1, const char *s2, const int extension_cost, const int fast,const int il_a, const int il_b, const int b_a, const int b_b)
{
  if(fast==1){
    printf("target upper bound %d: query lower bound %d (%5.2f)\n", max_pos-10, max_pos_j-10, ((double)max)/100);
  }
  else if(fast==2){
    int alignment_length2; alignment_length2 = MIN2(n1,n2);
    int begin_t=MAX2(11, max_pos-alignment_length2+1);/* 10 */
    int end_t  =MIN2(n1-10, max_pos+1);
    int begin_q=MAX2(11, max_pos_j-1); /* 10 */
    int end_q  =MIN2(n2-10, max_pos_j+alignment_length2-1);
    char *s3 = (char*) vrna_alloc(sizeof(char)*(end_t - begin_t +2 + 20));
    char *s4 = (char*) vrna_alloc(sizeof(char)*(end_q - begin_q +2 + 20));
    strcpy(s3,"NNNNNNNNNN");strcpy(s4,"NNNNNNNNNN");
    strncat(s3, (s1+begin_t-1),  end_t - begin_t +1);
    strncat(s4, (s2+begin_q-1) , end_q - begin_q +1);
    strcat(s3,"NNNNNNNNNN");strcat(s4,"NNNNNNNNNN");
    s3[end_t -begin_t +1 +20 ]='\0';
    s4[end_q -begin_q +1 +20]='\0';
    duplexT test;
    test = fduplexfold(s3, s4, extension_cost,il_a, il_b, b_a, b_b);
    int l1=strchr(test.structure, '&')-test.structure;
    printf("%s %3d,%-3d : %3d,%-3d (%5.2f) [%5.2f] i:%d,j:%d <%5.2f>\n", test.structure,
           begin_t-10+test.i-l1-10,
           begin_t-10+test.i-1-10,
           begin_q-10 + test.j-1-10 ,
           (begin_q -11) + test.j + (int)strlen(test.structure)-l1-2-10,
           test.energy, test.energy_backtrack,max_pos-10, max_pos_j-10,((double) max)/100);
    free(s3);free(s4);free(test.structure);
  }
  else{
    duplexT test;
    int alignment_length2; alignment_length2 = MIN2(n1,n2);
    int begin_t=MAX2(11, max_pos-alignment_length2+1);
    int end_t  =MIN2(n1-10, max_pos+1);
    int begin_q=MAX2(11, max_pos_j-1);
    int end_q  =MIN2(n2-10, max_pos_j+alignment_length2-1);
    char *s3 = (char*) vrna_alloc(sizeof(char)*(end_t - begin_t +2));
    char *s4 = (char*) vrna_alloc(sizeof(char)*(end_q - begin_q +2));
    strncpy(s3, (s1+begin_t-1),  end_t - begin_t + 1);
    strncpy(s4, (s2+begin_q-1) , end_q - begin_q +1 );
    s3[end_t -begin_t +1 ]='\0';
    s4[end_q -begin_q +1 ]='\0';
    test = duplexfold(s3, s4, extension_cost);
    int l1=strchr(test.structure, '&')-test.structure;
    printf("%s %3d,%-3d : %3d,%-3d (%5.2f) i:%d,j:%d <%5.2f>\n", test.structure,
           begin_t-10+test.i-l1,
           begin_t-10+test.i-1,
           begin_q-10 +test.j-1 ,
           (begin_q -11) + test.j + (int)strlen(test.structure)-l1-2,
           test.energy, max_pos-10, max_pos_j -10, ((double) max)/100);
    free(s3);free(s4);free(test.structure);
  }
}


PRIVATE void update_dfold_params(void)
{
  vrna_md_t md;
  if(P)
    free(P);
  set_model_details(&md);
  P = vrna_params(&md);
  make_pair_matrix();
}

PRIVATE void encode_seqs(const char *s1, const char *s2) {
  unsigned int i,l;

  l = strlen(s1);
  S1 = encode_seq(s1);
  SS1= (short *) vrna_alloc(sizeof(short)*(l+1));
  /* SS1 exists only for the special X K and I bases and energy_set!=0 */

  for (i=1; i<=l; i++) { /* make numerical encoding of sequence */
    SS1[i] = alias[S1[i]];   /* for mismatches of nostandard bases */
  }

  l = strlen(s2);
  S2 = encode_seq(s2);
  SS2= (short *) vrna_alloc(sizeof(short)*(l+1));
  /* SS2 exists only for the special X K and I bases and energy_set!=0 */

  for (i=1; i<=l; i++) { /* make numerical encoding of sequence */
    SS2[i] = alias[S2[i]];   /* for mismatches of nostandard bases */
  }
}

PRIVATE short * encode_seq(const char *sequence) {
  unsigned int i,l;
  short *S;
  l = strlen(sequence);
  S = (short *) vrna_alloc(sizeof(short)*(l+2));
  S[0] = (short) l;

  /* make numerical encoding of sequence */
  for (i=1; i<=l; i++)
    S[i]= (short) encode_char(toupper(sequence[i-1]));

  /* for circular folding add first base at position n+1 */
  S[l+1] = S[1];

  return S;
}

int arraySize(duplexT** array)
{
  int site_count=0;
  while(array[site_count]!=NULL){
    site_count++;
  }
  return site_count;
}

void freeDuplexT(duplexT** array)
{
  int size=arraySize(array);
  while(--size){
    free(array[size]->structure);
    free(array[size]);
  }
  free(array[0]->structure);
  free(array);
}
