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
#include "ViennaRNA/loop_energies.h"
#include "ViennaRNA/plex.h"
#include "ViennaRNA/ali_plex.h"


#define PUBLIC
#define PRIVATE static

#define STACK_BULGE1  1   /* stacking energies for bulges of size 1 */
#define NEW_NINIO     1   /* new asymetry penalty */
#define ARRAY 32          /*array size*/
#define UNIT 100
#define MINPSCORE -2 * UNIT
/**
*** Due to the great similarity between functions,
*** more annotation can be found in plex.c
**/

PRIVATE short *encode_seq(const char *seq);
PRIVATE void  update_dfold_params(void);
/**
*** aliduplexfold(_XS)/alibacktrack(_XS) computes duplex interaction with standard energy and considers extension_cost
*** alifind_max(_XS)/aliplot_max(_XS) find suboptimals and MFE
**/
PRIVATE duplexT aliduplexfold(const char *s1[], const char *s2[], const int extension_cost);
PRIVATE char *    alibacktrack(int i, int j, const short *s1[], const short *s2[], const int extension_cost);
PRIVATE void      alifind_max(const int *position, const int *position_j,const int delta, const int threshold,
                           const int alignment_length, const char *s1[], const char *s2[], const int extension_cost, const int fast);
PRIVATE void      aliplot_max(const int max, const int max_pos, const int max_pos_j,
                           const int alignment_length, const char *s1[], const char *s2[], const int extension_cost, const int fast);
PRIVATE duplexT aliduplexfold_XS(const char *s1[], const char *s2[],const int **access_s1,
                                 const int **access_s2, const int i_pos, const int j_pos, const int threshold,const int i_flag, const int j_flag);
PRIVATE char *    alibacktrack_XS(int i, int j, const short *s1[], const short *s2[], const int** access_s1, const int** access_s2,const int i_flag, const int j_flag);
PRIVATE void  alifind_max_XS(const int *position, const int *position_j,
                                 const int delta,  const int threshold, const int alignment_length,
                                 const char* s1[], const char* s2[],
                                 const int **access_s1, const int **access_s2, const int fast);
PRIVATE void aliplot_max_XS(const int max, const int max_pos, const int max_pos_j,
                            const int alignment_length, const char *s1[], const char* s2[],
                            const int **access_s1, const int **access_s2, const int fast);

/**
*** computes covariance score
**/

PRIVATE int covscore(const int *types, int n_seq);

extern double cv_fact; /* from alifold.c, default 1 */
extern double nc_fact;


/*@unused@*/

#define MAXSECTORS      500     /* dimension for a backtrack array */
#define LOCALITY        0.      /* locality parameter for base-pairs */

PRIVATE vrna_param_t *P = NULL;
PRIVATE int   **c = NULL;
PRIVATE int  **lc = NULL, **lin = NULL, **lbx = NULL, **lby = NULL,**linx = NULL, **liny = NULL;




PRIVATE int   n1,n2;
PRIVATE int n3, n4;
PRIVATE int delay_free=0;


/*-----------------------------------------------------------------------duplexfold_XS---------------------------------------------------------------------------*/



/*----------------------------------------------ALIDUPLEXFOLD-----------------------------------------------------------------------------------------------------------*/
PRIVATE duplexT aliduplexfold(const char *s1[], const char *s2[], const int extension_cost) {
  int i, j, s, n_seq, Emin=INF, i_min=0, j_min=0;
  char *struc;
  duplexT mfe;
  vrna_md_t   md;
  short **S1, **S2;
  int *type;
  n3 = (int) strlen(s1[0]);
  n4 = (int) strlen(s2[0]);
  for (s=0; s1[s]!=NULL; s++);
  n_seq = s;
  for (s=0; s2[s]!=NULL; s++);
  if (n_seq != s) vrna_message_error("unequal number of sequences in aliduplexfold()\n");

  set_model_details(&md);
  if ((!P) || (fabs(P->temperature - temperature)>1e-6)) {
    update_fold_params();
    if(P)
      free(P);
    P = vrna_params(&md);
    make_pair_matrix();
  }

  c = (int **) vrna_alloc(sizeof(int *) * (n3+1));
  for (i=1; i<=n3; i++) c[i] = (int *) vrna_alloc(sizeof(int) * (n4+1));

  S1 = (short **) vrna_alloc((n_seq+1)*sizeof(short *));
  S2 = (short **) vrna_alloc((n_seq+1)*sizeof(short *));
  for (s=0; s<n_seq; s++) {
    if (strlen(s1[s]) != n3) vrna_message_error("uneqal seqence lengths");
    if (strlen(s2[s]) != n4) vrna_message_error("uneqal seqence lengths");
    S1[s] = encode_seq(s1[s]);
    S2[s] = encode_seq(s2[s]);
  }
  type = (int *) vrna_alloc(n_seq*sizeof(int));

  for (i=1; i<=n3; i++) {
    for (j=n4; j>0; j--) {
      int k,l,E,psc;
      for (s=0; s<n_seq; s++) {
        type[s] = pair[S1[s][i]][S2[s][j]];
      }
      psc = covscore(type, n_seq);
      for (s=0; s<n_seq; s++) if (type[s]==0) type[s]=7;
      c[i][j] = (psc>=MINPSCORE) ? (n_seq*(P->DuplexInit + 2*extension_cost)) : INF;
      if (psc<MINPSCORE) continue;
      for (s=0; s<n_seq; s++) {
        c[i][j] += E_ExtLoop(type[s], (i>1) ? S1[s][i-1] : -1, (j<n4) ? S2[s][j+1] : -1, P) + 2*extension_cost;
      }
      for (k=i-1; k>0 && k>i-MAXLOOP-2; k--) {
        for (l=j+1; l<=n4; l++) {
          int type2;
          if (i-k+l-j-2>MAXLOOP) break;
          if (c[k][l]>INF/2) continue;
          for (E=s=0; s<n_seq; s++) {
            type2 = pair[S1[s][k]][S2[s][l]];
            if (type2==0) type2=7;
            E += E_IntLoop(i-k-1, l-j-1, type2, rtype[type[s]],
                           S1[s][k+1], S2[s][l-1], S1[s][i-1], S2[s][j+1],P) + (i-k+l-j)*extension_cost;
          }
          c[i][j] = MIN2(c[i][j], c[k][l]+E);
        }
      }
      c[i][j] -= psc;
      E = c[i][j];
      for (s=0; s<n_seq; s++) {
        E += E_ExtLoop(rtype[type[s]], (j>1) ? S2[s][j-1] : -1, (i<n3) ? S1[s][i+1] : -1, P) +2*extension_cost;
      }
      if (E<Emin) {
        Emin=E; i_min=i; j_min=j;
      }
    }
  }
  struc = alibacktrack(i_min, j_min, (const short int**) S1, (const short int**) S2 , extension_cost);
  if (i_min<n3) i_min++;
  if (j_min>1 ) j_min--;
  int size;
  size=strlen(struc)-1;
  Emin-=size * n_seq * extension_cost;
  mfe.i = i_min;
  mfe.j = j_min;
  mfe.energy = (float) (Emin/(100.*n_seq));
  mfe.structure = struc;
  if (!delay_free) {
    for (i=1; i<=n3; i++) free(c[i]);
    free(c);
  }
  for (s=0; s<n_seq; s++) {
    free(S1[s]); free(S2[s]);
  }
  free(S1); free(S2); free(type);
  return mfe;
}


PRIVATE char *alibacktrack(int i, int j, const short *S1[], const short *S2[], const int extension_cost) {
  /* backtrack structure going backwards from i, and forwards from j
     return structure in bracket notation with & as separator */
  int k, l, *type, type2, E, traced, i0, j0, s, n_seq;
  char *st1, *st2, *struc;

  n3 = (int) S1[0][0];
  n4 = (int) S2[0][0];

  for (s=0; S1[s]!=NULL; s++);
  n_seq = s;
  for (s=0; S2[s]!=NULL; s++);
  if (n_seq != s) vrna_message_error("unequal number of sequences in alibacktrack()\n");

  st1 = (char *) vrna_alloc(sizeof(char)*(n3+1));
  st2 = (char *) vrna_alloc(sizeof(char)*(n4+1));
  type = (int *) vrna_alloc(n_seq*sizeof(int));

  i0=MIN2(i+1,n3); j0=MAX2(j-1,1);

  while (i>0 && j<=n4) {
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
      for (l=j+1; l<=n4; l++) {
        int LE;
        if (i-k+l-j-2>MAXLOOP) break;
        if (c[k][l]>INF/2) continue;
        for (s=LE=0; s<n_seq; s++) {
          type2 = pair[S1[s][k]][S2[s][l]];
          if (type2==0) type2=7;
          LE += E_IntLoop(i-k-1, l-j-1, type2, rtype[type[s]],
                          S1[s][k+1], S2[s][l-1], S1[s][i-1], S2[s][j+1],P)+(i-k+l-j)*extension_cost;
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
        E -= E_ExtLoop(type[s], (i>1) ? S1[s][i-1] : -1, (j<n4) ? S2[s][j+1] : -1, P) + 2*extension_cost;
      }
      if (E != n_seq*P->DuplexInit + n_seq*2*extension_cost) {
        vrna_message_error("backtrack failed in aliduplex");
      } else break;
    }
  }
  if (i>1)  i--;
  if (j<n4) j++;

  struc = (char *) vrna_alloc(i0-i+1+j-j0+1+2);
  for (k=MAX2(i,1); k<=i0; k++) if (!st1[k-1]) st1[k-1] = '.';
  for (k=j0; k<=j; k++) if (!st2[k-1]) st2[k-1] = '.';
  strcpy(struc, st1+MAX2(i-1,0)); strcat(struc, "&");
  strcat(struc, st2+j0-1);

  /* printf("%s %3d,%-3d : %3d,%-3d\n", struc, i,i0,j0,j);  */
  free(st1); free(st2); free(type);

  return struc;
}

duplexT** aliLduplexfold(const char *s1[], const char *s2[], const int threshold, const int extension_cost, const int alignment_length, const int delta, const int fast,const int il_a, const int il_b, const int b_a, const int b_b)
{
  short **S1, **S2;
  int *type, type2;
  int i, j,s,n_seq;
  s=0;
  int bopen=b_b;
  int bext=b_a+extension_cost;
  int iopen=il_b;
  int iext_s=2*(il_a+extension_cost);/* iext_s 2 nt nucleotide extension of interior loop, on i and j side */
  int iext_ass=50+il_a+extension_cost;/* iext_ass assymetric extension of interior loop, either on i or on j side. */
  int min_colonne=INF; /* enthaelt das maximum einer kolonne */
  int i_length;
  int max_pos;/* get position of the best hit */
  int max_pos_j;
  int temp;
  int min_j_colonne;
  int max=INF;
  /* FOLLOWING NEXT 4 LINE DEFINES AN ARRAY CONTAINING POSITION OF THE SUBOPT IN S1 */
  int *position; /* contains the position of the hits with energy > E */
  int *position_j;


  n1 = (int) strlen(s1[0]);
  n2 = (int) strlen(s2[0]);
  for (s=0; s1[s]; s++);
  n_seq = s;
  for (s=0; s2[s]; s++);
  if (n_seq != s) vrna_message_error("unequal number of sequences in aliduplexfold()\n");

  position = (int *) vrna_alloc((delta+(n1)+4+delta) * sizeof(int));
  position_j= (int *) vrna_alloc((delta+(n1)+4+delta) * sizeof(int));

  if ((!P) || (fabs(P->temperature - temperature)>1e-6)){
    update_dfold_params();
  }

  lc   = (int**) vrna_alloc(sizeof(int *) * 5);
  lin  = (int**) vrna_alloc(sizeof(int *) * 5);
  lbx  = (int**) vrna_alloc(sizeof(int *) * 5);
  lby  = (int**) vrna_alloc(sizeof(int *) * 5);
  linx = (int**) vrna_alloc(sizeof(int *) * 5);
  liny = (int**) vrna_alloc(sizeof(int *) * 5);

  for (i=0; i<=4; i++){
    lc[i]  = (int *) vrna_alloc(sizeof(int) * (n2+5));
    lin[i] = (int *) vrna_alloc(sizeof(int) * (n2+5));
    lbx[i] = (int *) vrna_alloc(sizeof(int) * (n2+5));
    lby[i] = (int *) vrna_alloc(sizeof(int) * (n2+5));
    linx[i]= (int *) vrna_alloc(sizeof(int) * (n2+5));
    liny[i]= (int *) vrna_alloc(sizeof(int) * (n2+5));
  }


  S1 = (short **) vrna_alloc((n_seq+1)*sizeof(short *));
  S2 = (short **) vrna_alloc((n_seq+1)*sizeof(short *));
  for (s=0; s<n_seq; s++) {
    if (strlen(s1[s]) != n1) vrna_message_error("uneqal seqence lengths");
    if (strlen(s2[s]) != n2) vrna_message_error("uneqal seqence lengths");
    S1[s] = encode_seq(s1[s]);
    S2[s] = encode_seq(s2[s]);
  }
  type = (int *) vrna_alloc(n_seq*sizeof(int));
  /**
  *** array initialization
  **/
  for(j=n2;j>=0;j--) {
    lbx[0][j]=lbx[1][j]=lbx[2][j]=lbx[3][j]    = lbx[4][j] =INF;
    lin[0][j]=lin[1][j]=lin[2][j]=lin[3][j]    = lin[4][j] =INF;
    lc[0][j] =lc[1][j] =lc[2][j] = lc[3][j]    =  lc[4][j] =INF;
    lby[0][j]=lby[1][j]=lby[2][j]=lby[3][j]    = lby[4][j] =INF;
    liny[0][j]=liny[1][j]=liny[2][j]=liny[3][j]=liny[4][j]=INF;
    linx[0][j]=linx[1][j]=linx[2][j]=linx[3][j]=linx[4][j]=INF;
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
      int psc;
      for (s=0; s<n_seq; s++) {
        type[s] = pair[S1[s][i]][S2[s][j]];
      }
      psc = covscore(type, n_seq);
      for (s=0; s<n_seq; s++) if (type[s]==0) type[s]=7;
      lc[idx][j] = (psc>=MINPSCORE) ? (n_seq*P->DuplexInit + 2*n_seq*extension_cost) : INF;
      /**
      *** Update matrix. It is the average over all sequence of a given structure element
      *** c_stack -> stacking of c
      *** c_10, c01 -> stack from bulge
      *** c_nm -> arrives in stack from nxm loop
      *** c_in -> arrives in stack from interior loop
      *** c_bx -> arrives in stack from large bulge on target
      *** c_by -> arrives in stack from large bulge on query
      ***
      **/
      int c_stack, c_10, c_01, c_11, c_22, c_21, c_12, c_23, c_32, c_in, c_in2x, c_in2y, c_bx, c_by, c_inx, c_iny;  /* matrix c */
      int in, in_x, in_y, in_xy; /*  in begin, in_x assymetric, in_y assymetric, in_xy symetric; */
      int inx, inx_x;
      int iny, iny_y;
      int bx, bx_x;
      int by, by_y;
      in=lc[idx_1][j+1]; in_x=lin[idx_1][j]; in_y=lin[idx][j+1]; in_xy=lin[idx_1][j+1];
      inx=lc[idx_1][j+1]; inx_x=linx[idx_1][j];
      iny=lc[idx_1][j+1]; iny_y=liny[idx][j+1];
      bx=lc[idx_1][j]; bx_x=lbx[idx_1][j];
      by=lc[idx][j+1]; by_y=lby[idx][j+1];
      c_stack=lc[idx_1][j+1]; c_01=lc[idx_1][j+2];c_10=lc[idx_2][j+1];
      c_12=lc[idx_2][j+3];c_21=lc[idx_3][j+2];c_11=lc[idx_2][j+2];
      c_22=lc[idx_3][j+3];c_32=lc[idx_4][j+3];c_23=lc[idx_3][j+4];
      c_in=lin[idx_3][j+3];c_in2x=lin[idx_4][j+2];c_in2y=lin[idx_2][j+4];
      c_inx=linx[idx_3][j+1]; c_iny=liny[idx_1][j+3];
      c_bx=lbx[idx_2][j+1];c_by=lby[idx_1][j+2];
      for (s=0; s<n_seq; s++) {
        type2 = pair[S2[s][j+1]][S1[s][i-1]];
        in   +=P->mismatchI[type2][S2[s][j]][S1[s][i]]+iopen+iext_s;
        in_x +=iext_ass;
        in_y +=iext_ass;
        in_xy+=iext_s;
        inx  +=P->mismatch1nI[type2][S2[s][j]][S1[s][i]]+iopen+iext_s;
        inx_x+=iext_ass;
        iny  +=P->mismatch1nI[type2][S2[s][j]][S1[s][i]]+iopen+iext_s;
        iny_y+=iext_ass;
        type2=pair[S2[s][j]][S1[s][i-1]];
        bx   +=bopen+bext+(type2>2?P->TerminalAU:0);
        bx_x +=bext;
        type2=pair[S2[s][j+1]][S1[s][i]];
        by   +=bopen+bext+(type2>2?P->TerminalAU:0);
        by_y +=bext;
      }
      lin [idx][j]=MIN2(in, MIN2(in_x, MIN2(in_y, in_xy)));
      linx[idx][j]=MIN2(inx_x, inx);
      liny[idx][j]=MIN2(iny_y, iny);
      lby[idx][j] =MIN2(by, by_y);
      lbx[idx][j] =MIN2(bx, bx_x);

      if (psc<MINPSCORE) continue;
      for (s=0; s<n_seq; s++) {
        lc[idx][j]+=E_ExtLoop(type[s], S1[s][i-1],S2[s][j+1], P) + 2*extension_cost;
      }
      for (s=0; s<n_seq; s++) {
        type2=pair[S1[s][i-1]][S2[s][j+1]];if (type2==0) type2=7;
        c_stack+=E_IntLoop(0,0,type2, rtype[type[s]],S1[s][i], S2[s][j], S1[s][i-1], S2[s][j+1], P)+2*extension_cost;
        type2=pair[S1[s][i-1]][S2[s][j+2]];if (type2==0) type2=7;
        c_01   +=E_IntLoop(0,1,type2, rtype[type[s]],S1[s][i], S2[s][j+1], S1[s][i-1], S2[s][j+1], P)+3*extension_cost;
        type2=pair[S1[s][i-2]][S2[s][j+1]]; if (type2==0) type2=7;
        c_10   +=E_IntLoop(1,0,type2, rtype[type[s]],S1[s][i-1], S2[s][j], S1[s][i-1], S2[s][j+1], P)+3*extension_cost;
        type2=pair[S1[s][i-2]][S2[s][j+2]]; if (type2==0) type2=7;
        c_11   +=E_IntLoop(1,1,type2, rtype[type[s]],S1[s][i-1], S2[s][j+1], S1[s][i-1], S2[s][j+1], P)+4*extension_cost;
        type2 = pair[S1[s][i-3]][S2[s][j+3]];if (type2==0) type2=7;
        c_22   +=E_IntLoop(2,2,type2, rtype[type[s]],S1[s][i-2], S2[s][j+2], S1[s][i-1], S2[s][j+1], P)+6*extension_cost;
        type2 = pair[S1[s][i-3]][S2[s][j+2]];if (type2==0) type2=7;
        c_21   +=E_IntLoop(2,1,type2, rtype[type[s]],S1[s][i-2], S2[s][j+1], S1[s][i-1], S2[s][j+1], P)+5*extension_cost;
        type2 = pair[S1[s][i-2]][S2[s][j+3]];if (type2==0) type2=7;
        c_12   +=E_IntLoop(1,2,type2, rtype[type[s]],S1[s][i-1], S2[s][j+2], S1[s][i-1], S2[s][j+1], P)+5*extension_cost;
        type2 = pair[S1[s][i-4]][S2[s][j+3]];if (type2==0) type2=7;
        c_32   +=E_IntLoop(3,2,type2, rtype[type[s]],S1[s][i-3], S2[s][j+2], S1[s][i-1], S2[s][j+1], P)+7*extension_cost;
        type2 = pair[S1[s][i-3]][S2[s][j+4]];if (type2==0) type2=7;
        c_23   +=E_IntLoop(2,3,type2, rtype[type[s]],S1[s][i-2], S2[s][j+3], S1[s][i-1], S2[s][j+1], P)+7*extension_cost;
        c_in   +=P->mismatchI[rtype[type[s]]][S1[s][i-1]][S2[s][j+1]]+2*extension_cost+2*iext_s;
        c_in2x +=P->mismatchI[rtype[type[s]]][S1[s][i-1]][S2[s][j+1]]+iext_s+2*iext_ass+2*extension_cost;
        c_in2y +=P->mismatchI[rtype[type[s]]][S1[s][i-1]][S2[s][j+1]]+iext_s+2*iext_ass+2*extension_cost;
        c_inx  +=P->mismatch1nI[rtype[type[s]]][S1[s][i-1]][S2[s][j+1]]+iext_ass+iext_ass+2*extension_cost;
        c_iny  +=P->mismatch1nI[rtype[type[s]]][S1[s][i-1]][S2[s][j+1]]+iext_ass+iext_ass+2*extension_cost;
        int bAU;
        bAU=(type[s]>2?P->TerminalAU:0);
        c_bx   +=2*extension_cost+bext+bAU;
        c_by   +=2*extension_cost+bext+bAU;
      }
      lc[idx][j] =MIN2(lc[idx][j],
                       MIN2(c_stack,
                            MIN2(c_10,
                                 MIN2(c_01,
                                      MIN2(c_11,
                                           MIN2(c_21,
                                                MIN2(c_12,
                                                     MIN2(c_22,
                                                          MIN2(c_23,
                                                               MIN2(c_32,
                                                                    MIN2(c_bx,
                                                                         MIN2(c_by,
                                                                              MIN2(c_in,
                                                                                   MIN2(c_in2x,
                                                                                        MIN2(c_in2y,
                                                                                             MIN2(c_inx,c_iny)
                                                                                             )
                                                                                        )
                                                                                   )
                                                                              )
                                                                         )
                                                                    )
                                                               )
                                                          )
                                                     )
                                                )
                                           )
                                      )
                                 )
                            )
                       );
      lc[idx][j]-=psc;
      temp=lc[idx][j];
      for (s=0; s<n_seq; s++) {
        temp+=E_ExtLoop(rtype[type[s]], S2[s][j-1],S1[s][i+1],P)+2*extension_cost;
      }
      if(min_colonne > temp){
        min_colonne=temp;
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
  /* printf("MAX:%d ",max); */
  for (s=0; s<n_seq; s++) {free(S1[s]);free(S2[s]);}
  free(S1); free(S2);
  if(max<threshold){
    alifind_max(position, position_j, delta, threshold, alignment_length, s1, s2, extension_cost, fast);
  }
  aliplot_max(max, max_pos, max_pos_j,alignment_length, s1, s2, extension_cost,fast);
  for (i=0; i<=4; i++) {free(lc[i]);free(lin[i]);free(lbx[i]);free(lby[i]);free(linx[i]);free(liny[i]);}
  /* free(lc[0]);free(lin[0]);free(lbx[0]);free(lby[0]);free(linx[0]);free(liny[0]); */
  free(lc);free(lin);free(lbx);free(lby);free(linx);free(liny);
  free(position);
  free(position_j);
  free(type);
  return NULL;
}


PRIVATE void alifind_max(const int *position, const int *position_j,
                         const int delta, const int threshold, const int alignment_length,
                         const char *s1[], const char *s2[],
                         const int extension_cost, const int fast){
  int n_seq=0;
  for (n_seq=0; s1[n_seq]!=NULL; n_seq++);
  int pos=n1-9;
  if(fast==1){
    while(10<pos--){
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
        printf("target upper bound %d: query lower bound %d  (%5.2f) \n", pos-10, max_pos_j-10, ((double)max)/(n_seq*100));
        pos=MAX2(10,pos+temp_min-delta);
      }
    }
  }
  else{
    pos=n1-9;
    while(pos-- > 10) {
      /* printf("delta %d position:%d value:%d\n", delta, pos, position[pos]); */
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
        /* printf("%d %d %d\n", pos, max_pos_j,position[pos+delta]); */
        int begin_t=MAX2(11, pos-alignment_length+1);
        int end_t  =MIN2(n1-10, pos+1);
        int begin_q=MAX2(11, max_pos_j-1);
        int end_q  =MIN2(n2-10, max_pos_j+alignment_length-1);
        char **s3, **s4;
        s3 = (char**) vrna_alloc(sizeof(char*)*(n_seq+1));
        s4 = (char**) vrna_alloc(sizeof(char*)*(n_seq+1));
        int i;
        for(i=0; i<n_seq; i++){
          s3[i] = (char*) vrna_alloc(sizeof(char)*(end_t-begin_t+2));
          s4[i] = (char*) vrna_alloc(sizeof(char)*(end_q-begin_q+2));
          strncpy(s3[i], (s1[i]+begin_t-1), end_t - begin_t +1);
          strncpy(s4[i], (s2[i]+begin_q-1), end_q - begin_q +1);
          s3[i][end_t - begin_t +1]='\0';
          s4[i][end_q - begin_q +1]='\0';
        }
        duplexT test;
        test = aliduplexfold((const char**)s3, (const char**)s4, extension_cost);
        /* printf("test %d threshold %d",test.energy*100,(threshold/n_seq)); */
        if(test.energy * 100 < (int) (threshold/n_seq)){
          int l1=strchr(test.structure, '&')-test.structure;
                   printf("%s %3d,%-3d : %3d,%-3d (%5.2f)\n", test.structure,
                 begin_t -10 +test.i-l1,
                 begin_t -10 +test.i-1,
                 begin_q -10 + test.j - 1,
                 begin_q-11 + test.j + (int)strlen(test.structure) -l1 -2 , test.energy);
          pos=MAX2(10,pos+temp_min-delta);

        }
        for(i=0;i<n_seq;i++){
          free(s3[i]);free(s4[i]);
        }
        free(s3);free(s4);
        free(test.structure);
      }
    }
  }
}
PRIVATE void aliplot_max(const int max, const int max_pos, const int max_pos_j, const int alignment_length, const char *s1[], const char *s2[], const int extension_cost, const int fast)
{
  int n_seq;
  for (n_seq=0; !(s1[n_seq]==NULL); n_seq++);
  n1 = strlen(s1[0]); /* get length of alignment */
  n2 = strlen(s2[0]); /* get length of alignment */
  if(fast==1){
    printf("target upper bound %d: query lower bound %d (%5.2f)\n",
           max_pos-10, max_pos_j-10, (double) ((double)max)/(100*n_seq));
  }
  else{
    int begin_t=MAX2(11, max_pos-alignment_length+1);
    int end_t  =MIN2(n1-10, max_pos+1);
    int begin_q=MAX2(11, max_pos_j-1);
    int end_q  =MIN2(n2-10, max_pos_j+alignment_length-1);
    char **s3, **s4;
    s3 = (char**) vrna_alloc(sizeof(char*)*(n_seq+1));
    s4 = (char**) vrna_alloc(sizeof(char*)*(n_seq+1));
    int i;
    for(i=0; i<n_seq; i++){
      s3[i] = (char*) vrna_alloc(sizeof(char)*(end_t-begin_t+2));
      s4[i] = (char*) vrna_alloc(sizeof(char)*(end_q-begin_q+2));
      strncpy(s3[i], (s1[i]+begin_t-1), end_t - begin_t +1);
      strncpy(s4[i], (s2[i]+begin_q-1), end_q - begin_q +1);
      s3[i][end_t - begin_t +1]='\0';
      s4[i][end_q - begin_q +1]='\0';
    }
    duplexT test;
    s3[n_seq]=s4[n_seq]=NULL;
    test = aliduplexfold((const char**) s3,(const char**) s4, extension_cost);
    int l1=strchr(test.structure, '&')-test.structure;
    printf("%s %3d,%-3d : %3d,%-3d (%5.2f)\n",
           test.structure,
           begin_t -10 +test.i-l1,
           begin_t -10 +test.i-1,
           begin_q-10 + test.j - 1,
           begin_q -11 + test.j + (int)strlen(test.structure) - l1 - 2,
           test.energy);
    for(i=0; i<n_seq ; i++){
      free(s3[i]);free(s4[i]);
    }
    free(s3);free(s4);
    free(test.structure);
  }
}

PRIVATE duplexT aliduplexfold_XS(const char *s1[], const char *s2[],
                                 const int **access_s1, const int **access_s2,
                                 const int i_pos, const int j_pos, const int threshold,
                                 const int i_flag, const int j_flag){
  int i,j,s,p,q, Emin=INF, l_min=0, k_min=0;
  char *struc;
  short **S1,**S2;
  int *type,*type2;
  struc=NULL;
  duplexT mfe;
  vrna_md_t   md;
  int n_seq;
  n3 = (int) strlen(s1[0]);
  n4 = (int) strlen(s2[0]);
  for (s=0; s1[s]!=NULL; s++);
  n_seq = s;
  for (s=0; s2[s]!=NULL; s++);
  /* printf("%d \n",i_pos); */

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
  S1 = (short **) vrna_alloc((n_seq+1)*sizeof(short *));
  S2 = (short **) vrna_alloc((n_seq+1)*sizeof(short *));
  for (s=0; s<n_seq; s++) {
    if (strlen(s1[s]) != n3) vrna_message_error("uneqal seqence lengths");
    if (strlen(s2[s]) != n4) vrna_message_error("uneqal seqence lengths");
    S1[s] = encode_seq(s1[s]);
    S2[s] = encode_seq(s2[s]);
  }
  type =  (int *) vrna_alloc(n_seq*sizeof(int));
  type2 = (int *) vrna_alloc(n_seq*sizeof(int));
  int type3, E, k,l;
  i=n3-i_flag; j=1+j_flag;
  for (s=0; s<n_seq; s++) {
    type[s] = pair[S1[s][i]][S2[s][j]];
  }
  c[i][j] = n_seq*P->DuplexInit - covscore(type,n_seq);
  for (s=0; s<n_seq; s++) if (type[s]==0) type[s]=7;
  for (s=0; s<n_seq; s++) {
    c[i][j]+=E_ExtLoop(rtype[type[s]], (j_flag ? S2[s][j-1] : -1) , (i_flag ? S1[s][i+1] : -1),  P);
  }
  k_min=i; l_min=j; Emin=c[i][j];
  for (k=i; k>1 ; k--) {
    if(k<i) c[k+1][0]=INF;
    for (l=j; l<=n4-1; l++) {
      if(!(k==i && l==j)){
        c[k][l]=INF;
      }
      int psc2;
      for(s=0;s<n_seq;s++){
        type2[s] = pair[S1[s][k]][S2[s][l]];
      }
      psc2=covscore(type2, n_seq);
      if (psc2<MINPSCORE) continue;
      for (s=0; s<n_seq; s++) if (type2[s]==0) type2[s]=7;
      for (p=k+1; p<= n3 -i_flag && p<k+MAXLOOP-1; p++) {
        for (q = l-1; q >= 1+j_flag; q--) {
          if (p-k+l-q-2>MAXLOOP) break;
          for(E=s=0;s<n_seq;s++){
            type3=pair[S1[s][p]][S2[s][q]];
            if(type3==0) type3=7;
            E += E_IntLoop(p-k-1, l-q-1, type2[s], rtype[type3],
                           S1[s][k+1], S2[s][l-1], S1[s][p-1], S2[s][q+1],P);
          }
          c[k][l] = MIN2(c[k][l], c[p][q]+E);
        }
      }
      c[k][l]-=psc2;
      E = c[k][l];
      E+=n_seq*(access_s1[i-k+1][i_pos]+access_s2[l-1][j_pos+(l-1)-1]);
      for (s=0; s<n_seq; s++) {
        E+=E_ExtLoop(type2[s], (k>1) ? S1[s][k-1] : -1, (l<n4) ? S2[s][l+1] : -1, P);
      }
      if (E<Emin) {
        Emin=E; k_min=k; l_min=l;
      }
    }
  }
  if(Emin > threshold-1){
    mfe.structure=NULL;
    mfe.energy=INF;
    for (i=0; i<=n3; i++) free(c[i]);
    free(c);
    for(i=0; i<=n_seq;i++){
      free(S1[i]);
      free(S2[i]);
    }
    free(S1); free(S2); /* free(SS1); free(SS2); */
    free(type);free(type2);
    return mfe;
    } else{
    struc = alibacktrack_XS(k_min, l_min,(const short int**)S1,(const short int**)S2,access_s1, access_s2,i_flag,j_flag);
    }
  int dx_5, dx_3, dy_5, dy_3,dGx,dGy;
  dx_5=0; dx_3=0; dy_5=0; dy_3=0;dGx=0;dGy=0;
  dGx =n_seq*(access_s1[i-k_min+1][i_pos]);dx_3=0; dx_5=0;
  dGy =n_seq*access_s2[l_min-j+1][j_pos + (l_min-1)-1];
  mfe.tb=i_pos -9 - i + k_min -1 -dx_5;
  mfe.te=i_pos -9 -1 + dx_3;
  mfe.qb=j_pos -9 -1 - dy_5;
  mfe.qe=j_pos + l_min -3 -9 + dy_3;
  mfe.ddG=(double) (Emin *0.01);
  mfe.dG1=(double) (dGx*(0.01));
  mfe.dG2=(double) (dGy*(0.01));
  mfe.energy= mfe.ddG - mfe.dG1 - mfe.dG2;
  mfe.structure = struc;
  for (i=0; i<=n3; i++) free(c[i]);
  free(c);
  for(i=0; i<=n_seq;i++){
    free(S1[i]);
    free(S2[i]);
  }
  free(S1); free(S2); free(type);free(type2);
  return mfe;
}

PRIVATE char *alibacktrack_XS(int i, int j, const short *S1[], const short *S2[],
                              const int** access_s1, const int ** access_s2, const int i_flag, const int j_flag) {
  int n3,n4,k, l, *type, type2, E, traced, i0, j0,s,n_seq,psc;
  char *st1=NULL, *st2=NULL, *struc=NULL;

  n3 = (int) S1[0][0];
  n4 = (int) S2[0][0];
  for (s=0; S1[s]!=NULL; s++);
  n_seq = s;
  for (s=0; S2[s]!=NULL; s++);
  if (n_seq != s) vrna_message_error("unequal number of sequences in alibacktrack()\n");

  st1 = (char *) vrna_alloc(sizeof(char)*(n3+1));
  st2 = (char *) vrna_alloc(sizeof(char)*(n4+1));
  type = (int *) vrna_alloc(n_seq*sizeof(int));

  i0=i;/*MAX2(i-1,1);*/j0=j;/*MIN2(j+1,n4);*/
  while (i<=n3-i_flag && j>=1+j_flag) {
    E = c[i][j]; traced=0;
    st1[i-1] = '(';
    st2[j-1] = ')';
    for (s=0; s<n_seq; s++) {
      type[s] = pair[S1[s][i]][S2[s][j]];
    }
    psc = covscore(type,n_seq);
    for (s=0; s<n_seq; s++) if (type[s]==0) type[s] = 7;
    E += psc;
    for (k=i+1; k<=n3 && k>i-MAXLOOP-2; k++) {
      for (l=j-1; l>=1; l--) {
        int LE;
        if (i-k+l-j-2>MAXLOOP) break;
        for (s=LE=0; s<n_seq; s++) {
          type2 = pair[S1[s][k]][S2[s][l]];
          if (type2==0) type2=7;
          LE += E_IntLoop(k-i-1, j-l-1, type[s], rtype[type2],
                          S1[s][i+1], S2[s][j-1], S1[s][k-1], S2[s][l+1],P);
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
        if (type[s]>2) E -= P->TerminalAU;
      }
      break;
      if (E != n_seq*P->DuplexInit) {
        vrna_message_error("backtrack failed in fold duplex bal");
      } else break;
    }
  }
  struc = (char *) vrna_alloc(i-i0+1+j0-j+1+2);
  for (k=MAX2(i0,1); k<=i; k++) if (!st1[k-1]) st1[k-1] = '.';
  for (k=j; k<=j0; k++) if (!st2[k-1]) st2[k-1] = '.';
  strcpy(struc, st1+MAX2(i0-1,0)); strcat(struc, "&");
  strcat(struc, st2+j-1);
  free(st1);
  free(st2);
  free(type);
  return struc;
}

duplexT** aliLduplexfold_XS(const char*s1[], const char* s2[], const int **access_s1, const int **access_s2, const int threshold, const int alignment_length, const int delta,  const int fast,const int il_a, const int il_b, const int b_a, const int b_b)
{
  short **S1, **S2;
  int *type,type2;
  int i, j,s,n_seq;
  s=0;
  int bopen=b_b;
  int bext=b_a;
  int iopen=il_b;
  int iext_s=2*(il_a);/* iext_s 2 nt nucleotide extension of interior loop, on i and j side */
  int iext_ass=50+il_a;/* iext_ass assymetric extension of interior loop, either on i or on j side. */
  int min_colonne=INF; /* enthaelt das maximum einer kolonne */
  int i_length;
  int max_pos;/* get position of the best hit */
  int max_pos_j;
  int temp;
  int min_j_colonne;
  int max=INF;
  /* FOLLOWING NEXT 4 LINE DEFINES AN ARRAY CONTAINING POSITION OF THE SUBOPT IN S1 */
  int *position; /* contains the position of the hits with energy > E */
  int *position_j;
  int maxPenalty[4];

  n1 = (int) strlen(s1[0]);
  n2 = (int) strlen(s2[0]);
  for (s=0; s1[s]; s++);
  n_seq = s;
  for (s=0; s2[s]; s++);
  if (n_seq != s) vrna_message_error("unequal number of sequences in aliduplexfold()\n");

  position = (int *) vrna_alloc((delta+(n1)+4+delta) * sizeof(int));
  position_j= (int *) vrna_alloc((delta+(n1)+4+delta) * sizeof(int));

  /**
  *** extension penalty, computed only once, further reduce the computation time
  **/
  if ((!P) || (fabs(P->temperature - temperature)>1e-6)){
    update_dfold_params();
  }
  maxPenalty[0]=(int) -1*P->stack[2][2]/2;
  maxPenalty[1]=(int) -1*P->stack[2][2];
  maxPenalty[2]=(int) -3*P->stack[2][2]/2;
  maxPenalty[3]=(int) -2*P->stack[2][2];


  lc   = (int**) vrna_alloc(sizeof(int *) * 5);
  lin  = (int**) vrna_alloc(sizeof(int *) * 5);
  lbx  = (int**) vrna_alloc(sizeof(int *) * 5);
  lby  = (int**) vrna_alloc(sizeof(int *) * 5);
  linx = (int**) vrna_alloc(sizeof(int *) * 5);
  liny = (int**) vrna_alloc(sizeof(int *) * 5);

  for (i=0; i<=4; i++){
    lc[i]  = (int *) vrna_alloc(sizeof(int) * (n2+5));
    lin[i] = (int *) vrna_alloc(sizeof(int) * (n2+5));
    lbx[i] = (int *) vrna_alloc(sizeof(int) * (n2+5));
    lby[i] = (int *) vrna_alloc(sizeof(int) * (n2+5));
    linx[i]= (int *) vrna_alloc(sizeof(int) * (n2+5));
    liny[i]= (int *) vrna_alloc(sizeof(int) * (n2+5));
  }


  S1 = (short **) vrna_alloc((n_seq+1)*sizeof(short *));
  S2 = (short **) vrna_alloc((n_seq+1)*sizeof(short *));
  for (s=0; s<n_seq; s++) {
    if (strlen(s1[s]) != n1) vrna_message_error("uneqal seqence lengths");
    if (strlen(s2[s]) != n2) vrna_message_error("uneqal seqence lengths");
    S1[s] = encode_seq(s1[s]);
    S2[s] = encode_seq(s2[s]);
  }
  type = (int *) vrna_alloc(n_seq*sizeof(int));
  /**
  *** array initialization
  **/

  for(j=n2+4;j>=0;j--) {
    lbx[0][j]=lbx[1][j]=lbx[2][j]=lbx[3][j]    = lbx[4][j] =INF;
    lin[0][j]=lin[1][j]=lin[2][j]=lin[3][j]    = lin[4][j] =INF;
    lc[0][j] =lc[1][j] =lc[2][j] = lc[3][j]    =  lc[4][j] =INF;
    lby[0][j]=lby[1][j]=lby[2][j]=lby[3][j]    = lby[4][j] =INF;
    liny[0][j]=liny[1][j]=liny[2][j]=liny[3][j]=liny[4][j]=INF;
    linx[0][j]=linx[1][j]=linx[2][j]=linx[3][j]=linx[4][j]=INF;
  }
  i=10;
  i_length= n1 - 9 ;
  int di1,di2,di3,di4; /* contains accessibility penalty */
  while(i < i_length) {
    int idx=i%5;
    int idx_1=(i-1)%5;
    int idx_2=(i-2)%5;
    int idx_3=(i-3)%5;
    int idx_4=(i-4)%5;

    di1 = access_s1[5][i]   - access_s1[4][i-1];
    di2 = access_s1[5][i-1] - access_s1[4][i-2] + di1;
    di3 = access_s1[5][i-2] - access_s1[4][i-3] + di2;
    di4 = access_s1[5][i-3] - access_s1[4][i-4] + di3;
    di1=MIN2(di1,maxPenalty[0]);
    di2=MIN2(di2,maxPenalty[1]);
    di3=MIN2(di3,maxPenalty[2]);
    di4=MIN2(di3,maxPenalty[3]);
    j=n2-9;
    while (9 < --j) {
      int dj1,dj2,dj3,dj4;
      dj1 = access_s2[5][j+4] - access_s2[4][j+4];
      dj2 = access_s2[5][j+5] - access_s2[4][j+5] + dj1;
      dj3 = access_s2[5][j+6] - access_s2[4][j+6] + dj2;
      dj4 = access_s2[5][j+7] - access_s2[4][j+7] + dj3;
      dj1=MIN2(dj1,maxPenalty[0]);
      dj2=MIN2(dj2,maxPenalty[1]);
      dj3=MIN2(dj3,maxPenalty[2]);
      dj4=MIN2(dj4,maxPenalty[3]);
      int psc;
      for (s=0; s<n_seq; s++) { /* initialize type1 */
        type[s] = pair[S1[s][i]][S2[s][j]];
      }
      psc = covscore(type, n_seq);
      for (s=0; s<n_seq; s++) if (type[s]==0) type[s]=7;
      lc[idx][j] = ((psc >= MINPSCORE) ? n_seq*(P->DuplexInit) : INF);
      int c_stack, c_10, c_01, c_11, c_22, c_21, c_12, c_23, c_32, c_in, c_in2x, c_in2y, c_bx, c_by, c_inx, c_iny;  /* matrix c */
      int in,  in_x, in_y, in_xy; /*  in begin, in_x assymetric, in_y assymetric, in_xy symetric; */
      int inx, inx_x;
      int iny, iny_y;
      int bx,  bx_x;
      int by,  by_y;
      in=lc[idx_1][j+1]; in_x=lin[idx_1][j]; in_y=lin[idx][j+1]; in_xy=lin[idx_1][j+1];
      inx=lc[idx_1][j+1]; inx_x=linx[idx_1][j];
      iny=lc[idx_1][j+1]; iny_y=liny[idx][j+1];
      bx=lc[idx_1][j]; bx_x=lbx[idx_1][j];
      by=lc[idx][j+1]; by_y=lby[idx][j+1];
      c_stack=lc[idx_1][j+1]; c_01=lc[idx_1][j+2];c_10=lc[idx_2][j+1];
      c_12=lc[idx_2][j+3];c_21=lc[idx_3][j+2];c_11=lc[idx_2][j+2];
      c_22=lc[idx_3][j+3];c_32=lc[idx_4][j+3];c_23=lc[idx_3][j+4];
      c_in=lin[idx_3][j+3];c_in2x=lin[idx_4][j+2];c_in2y=lin[idx_2][j+4];
      c_inx=linx[idx_3][j+1]; c_iny=liny[idx_1][j+3];
      c_bx=lbx[idx_2][j+1];c_by=lby[idx_1][j+2];
      for (s=0; s<n_seq; s++) {
        type2 = pair[S2[s][j+1]][S1[s][i-1]];
        in   +=P->mismatchI[type2][S2[s][j]][S1[s][i]]+iopen+iext_s+di1+dj1;
        in_x +=iext_ass+di1;
        in_y +=iext_ass+dj1;
        in_xy+=iext_s+di1+dj1;
        inx  +=P->mismatch1nI[type2][S2[s][j]][S1[s][i]]+iopen+iext_s+di1+dj1;
        inx_x+=iext_ass+di1;
        iny  +=P->mismatch1nI[type2][S2[s][j]][S1[s][i]]+iopen+iext_s+di1+dj1;
        iny_y+=iext_ass+dj1;
        type2=pair[S2[s][j]][S1[s][i-1]];
        bx   +=bopen+bext+(type2>2?P->TerminalAU:0)+di1;
        bx_x +=bext+di1;
        type2=pair[S2[s][j+1]][S1[s][i]];
        by   +=bopen+bext+(type2>2?P->TerminalAU:0)+dj1;
        by_y +=bext+dj1;
      }
      lin[idx][j]=MIN2(in, MIN2(in_x, MIN2(in_y, in_xy)));
      linx[idx][j]=MIN2(inx_x, inx);
      liny[idx][j]=MIN2(iny_y, iny);
      lby[idx][j]=MIN2(by, by_y);
      lbx[idx][j]=MIN2(bx, bx_x);
      if (psc<MINPSCORE) continue;
      for (s=0; s<n_seq; s++) {
        lc[idx][j]+=E_ExtLoop(type[s], S1[s][i-1],S2[s][j+1],P);
      }
      for (s=0; s<n_seq; s++) {
        type2=pair[S1[s][i-1]][S2[s][j+1]];if (type2==0) type2=7;
        c_stack+=E_IntLoop(0,0,type2, rtype[type[s]],S1[s][i], S2[s][j], S1[s][i-1], S2[s][j+1], P)+di1+dj1;
        type2=pair[S1[s][i-1]][S2[s][j+2]];if (type2==0) type2=7;
        c_01   +=E_IntLoop(0,1,type2, rtype[type[s]],S1[s][i], S2[s][j+1], S1[s][i-1], S2[s][j+1], P)+di1+dj2;
        type2=pair[S1[s][i-2]][S2[s][j+1]];if (type2==0) type2=7;
        c_10   +=E_IntLoop(1,0,type2, rtype[type[s]],S1[s][i-1], S2[s][j], S1[s][i-1], S2[s][j+1], P)+di2+dj1;
        type2=pair[S1[s][i-2]][S2[s][j+2]];if (type2==0) type2=7;
        c_11   +=E_IntLoop(1,1,type2, rtype[type[s]],S1[s][i-1], S2[s][j+1], S1[s][i-1], S2[s][j+1], P)+di2+dj2;
        type2=pair[S1[s][i-3]][S2[s][j+3]];if (type2==0) type2=7;
        c_22   +=E_IntLoop(2,2,type2, rtype[type[s]],S1[s][i-2], S2[s][j+2], S1[s][i-1], S2[s][j+1], P)+di3+dj3;
        type2= pair[S1[s][i-3]][S2[s][j+2]];if (type2==0) type2=7;
        c_21   +=E_IntLoop(2,1,type2, rtype[type[s]],S1[s][i-2], S2[s][j+1], S1[s][i-1], S2[s][j+1], P)+di3+dj2;
        type2= pair[S1[s][i-2]][S2[s][j+3]];if (type2==0) type2=7;
        c_12   +=E_IntLoop(1,2,type2, rtype[type[s]],S1[s][i-1], S2[s][j+2], S1[s][i-1], S2[s][j+1], P)+di2+dj3;
        type2 = pair[S1[s][i-4]][S2[s][j+3]];if (type2==0) type2=7;
        c_32   +=E_IntLoop(3,2,type2, rtype[type[s]],S1[s][i-3], S2[s][j+2], S1[s][i-1], S2[s][j+1], P)+di4+dj3;
        type2 = pair[S1[s][i-3]][S2[s][j+4]];if (type2==0) type2=7;
        c_23   +=E_IntLoop(2,3,type2, rtype[type[s]],S1[s][i-2], S2[s][j+3], S1[s][i-1], S2[s][j+1], P)+dj4+di3;
        c_in   +=P->mismatchI[rtype[type[s]]][S1[s][i-1]][S2[s][j+1]]+di3+dj3+2*iext_s;
        c_in2x +=P->mismatchI[rtype[type[s]]][S1[s][i-1]][S2[s][j+1]]+iext_s+2*iext_ass+di4+dj2;
        c_in2y +=P->mismatchI[rtype[type[s]]][S1[s][i-1]][S2[s][j+1]]+iext_s+2*iext_ass+di2+dj4;
        c_inx  +=P->mismatch1nI[rtype[type[s]]][S1[s][i-1]][S2[s][j+1]]+iext_ass+iext_ass+di3+dj1;
        c_iny  +=P->mismatch1nI[rtype[type[s]]][S1[s][i-1]][S2[s][j+1]]+iext_ass+iext_ass+di1+dj3;
        int bAU;
        bAU=(type[s]>2?P->TerminalAU:0);
        c_bx   +=di2+dj1+bext+bAU;
        c_by   +=di1+dj2+bext+bAU;
      }
      lc[idx][j] =MIN2(lc[idx][j],
                       MIN2(c_stack,
                            MIN2(c_10,
                                 MIN2(c_01,
                                      MIN2(c_11,
                                           MIN2(c_21,
                                                MIN2(c_12,
                                                     MIN2(c_22,
                                                          MIN2(c_23,
                                                               MIN2(c_32,
                                                                    MIN2(c_bx,
                                                                         MIN2(c_by,
                                                                              MIN2(c_in,
                                                                                   MIN2(c_in2x,
                                                                                        MIN2(c_in2y,
                                                                                             MIN2(c_inx,c_iny)
                                                                                             )
                                                                                        )
                                                                                   )
                                                                              )
                                                                         )
                                                                    )
                                                               )
                                                          )
                                                     )
                                                )
                                           )
                                      )
                                 )
                            )
                       );
      lc[idx][j]-=psc;
      temp=lc[idx][j];

      for (s=0; s<n_seq; s++) {
        temp+=E_ExtLoop(rtype[type[s]], S2[s][j-1],S1[s][i+1],P);
      }
      if(min_colonne > temp){
        min_colonne=temp;
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
  for (s=0; s<n_seq; s++) {free(S1[s]);free(S2[s]);}
  free(S1); free(S2);
  if(max<threshold){
    alifind_max_XS(position, position_j, delta, threshold, alignment_length, s1, s2,access_s1,access_s2, fast);
  }
  aliplot_max_XS(max, max_pos, max_pos_j,alignment_length, s1, s2, access_s1, access_s2, fast);
  for (i=0; i<=4; i++) {free(lc[i]);free(lin[i]);free(lbx[i]);free(lby[i]);free(linx[i]);free(liny[i]);}
  /* free(lc[0]);free(lin[0]);free(lbx[0]);free(lby[0]);free(linx[0]);free(liny[0]); */
  free(lc);free(lin);free(lbx);free(lby);free(linx);free(liny);
  free(position);
  free(position_j);
  free(type);
  return NULL;
}




PRIVATE void alifind_max_XS(const int *position, const int *position_j,
                                const int delta, const int threshold, const int alignment_length,
                                const char *s1[], const char *s2[],
                                const int **access_s1, const int **access_s2, const int fast){

  int n_seq=0;
  for (n_seq=0; s1[n_seq]!=NULL; n_seq++);
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
        /*
        int max_pos_j;
        max_pos_j=position_j[pos+delta];
        int max;
        max=position[pos+delta];
        printf("target upper bound %d: query lower bound %d  (%5.2f) \n", pos-10, max_pos_j-10, ((double)max)/100);
        */
        pos=MAX2(10,pos+temp_min-delta);
      }
    }
  }
  else{
    pos=n1-9;
    while( pos-- > 10 ) {
      /* printf("delta %d position:%d value:%d\n", delta, pos, position[pos]); */
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
        max_pos_j=position_j[pos+delta]; /* position on j */
        int begin_t=MAX2(11, pos-alignment_length);
        int end_t  =MIN2(n1-10, pos+1);
        int begin_q=MAX2(11, max_pos_j-1);
        int end_q  =MIN2(n2-10, max_pos_j+alignment_length-1);
        int i_flag;
        int j_flag;
        i_flag = (end_t   ==  pos+1?1:0);
        j_flag = (begin_q == max_pos_j-1?1:0);
        char **s3, **s4;
        s3 = (char**) vrna_alloc(sizeof(char*)*(n_seq+1));
        s4 = (char**) vrna_alloc(sizeof(char*)*(n_seq+1));
        int i;
        for(i=0; i<n_seq; i++){
          s3[i] = (char*) vrna_alloc(sizeof(char)*(end_t-begin_t+2));
          s4[i] = (char*) vrna_alloc(sizeof(char)*(end_q-begin_q+2));
          strncpy(s3[i], (s1[i]+begin_t), end_t - begin_t +1);
          strncpy(s4[i], (s2[i]+begin_q), end_q - begin_q +1);
          s3[i][end_t - begin_t +1]='\0';
          s4[i][end_q - begin_q +1]='\0';
        }
        duplexT test;
        test = aliduplexfold_XS((const char**) s3, (const char**) s4, access_s1, access_s2,pos, max_pos_j,threshold,i_flag,j_flag);
        /*         printf("position %d approximation %d test %d threshold %d\n", pos, position[pos+delta], (int)test.energy,(int)(threshold/n_seq)); */
        if(test.energy*100  < (int) (threshold/n_seq)){
        printf("%s %3d,%-3d: %3d,%-3d (%5.2f = %5.2f + %5.2f + %5.2f)\n", test.structure,
               test.tb,test.te,test.qb,test.qe, test.ddG/n_seq, test.energy/n_seq, test.dG1/n_seq, test.dG2/n_seq);
        free(test.structure);
	pos=MAX2(10,pos+temp_min-delta);
        }
        for(i=0;i<n_seq;i++){
          free(s3[i]);free(s4[i]);
        }
        free(s3);free(s4);
        /* free(test.structure); */
      }

    }
  }
}

PRIVATE void aliplot_max_XS(const int max, const int max_pos, const int max_pos_j,const int alignment_length,
                            const char *s1[], const char* s2[],const int **access_s1, const int **access_s2,
                            const int fast){

  int n_seq;
  for (n_seq=0; !(s1[n_seq]==NULL); n_seq++);
  n1 = strlen(s1[0]); /* get length of alignment */
  n2 = strlen(s2[0]); /* get length of alignme */

  if(fast){
    printf("target upper bound %d: query lower bound %d (%5.2f)\n",
           max_pos-10, max_pos_j-10, (double) ((double)max)/(100*n_seq));
  }
  else{
    int begin_t=MAX2(11, max_pos-alignment_length); /* only get the position that binds.. */
    int end_t  =MIN2(n1-10, max_pos+1);              /* ..no dangles */
    int begin_q=MAX2(11, max_pos_j-1);
    int end_q  =MIN2(n2-10, max_pos_j+alignment_length-1);
    int i_flag;
    int j_flag;
    i_flag = (end_t   == max_pos+1?1:0);
    j_flag = (begin_q == max_pos_j-1?1:0);
    char **s3, **s4;
    s3 = (char**) vrna_alloc(sizeof(char*)*(n_seq+1));
    s4 = (char**) vrna_alloc(sizeof(char*)*(n_seq+1));
    int i;
    for(i=0; i<n_seq; i++){
      s3[i] = (char*) vrna_alloc(sizeof(char)*(end_t-begin_t+2));
      s4[i] = (char*) vrna_alloc(sizeof(char)*(end_q-begin_q+2));
      strncpy(s3[i], (s1[i]+begin_t), end_t - begin_t +1);
      strncpy(s4[i], (s2[i]+begin_q), end_q - begin_q +1);
      s3[i][end_t - begin_t +1]='\0';
      s4[i][end_q - begin_q +1]='\0';
    }
    duplexT test;
    test = aliduplexfold_XS((const char**) s3, (const char**) s4, access_s1, access_s2,max_pos, max_pos_j,INF,i_flag,j_flag);
    printf("%s %3d,%-3d: %3d,%-3d (%5.2f = %5.2f + %5.2f + %5.2f)\n", test.structure,
           test.tb,test.te,test.qb,test.qe, test.ddG/n_seq, test.energy/n_seq, test.dG1/n_seq, test.dG2/n_seq);
    free(test.structure);
    for(i=0;i<n_seq;i++){
      free(s3[i]);free(s4[i]);
    }
    free(s3);free(s4);
  }
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


PRIVATE void update_dfold_params(void)
{
  vrna_md_t md;
  if(P) free(P);
  set_model_details(&md);
  P = vrna_params(&md);
  make_pair_matrix();
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
