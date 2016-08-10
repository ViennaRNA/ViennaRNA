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

/* int subopt_sorted=0; */

#define PUBLIC
#define PRIVATE static

#define STACK_BULGE1  1   /* stacking energies for bulges of size 1 */
#define NEW_NINIO     1   /* new asymetry penalty */
#define ARRAY 32          /*array size*/
#define UNIT 100
#define MINPSCORE -2 * UNIT
PRIVATE void  encode_seqs(const char *s1, const char *s2);
PRIVATE short *encode_seq(const char *seq);
/* PRIVATE void  my_encode_seq(const char *s1, const char *s2); */
PRIVATE void  update_dfold_params(void);
/* PRIVATE int   compare(const void *sub1, const void *sub2); */
/* PRIVATE int   compare_XS(const void *sub1, const void *sub2); */
/* PRIVATE duplexT* backtrack(int threshold, const int extension_cost); */
/* static void  print_struct(duplexT const *dup); */

/* PRIVATE int   print_struct(duplexT const *dup); */
/* PRIVATE int   get_rescaled_energy(duplexT const *dup); */

PRIVATE char * backtrack_C(int i, int j, const int extension_cost, const char * structure, int *E);
PRIVATE void   find_max_C(const int *position, const int *position_j, const int delta, const int threshold, const int constthreshold, const int length, const char *s1, const char *s2, const int extension_cost, const int fast, const char* structure);
PRIVATE void   plot_max_C(const int max, const int max_pos, const int max_pos_j, const int alignment_length, const char *s1, const char *s2, const int extension_cost, const int fast, const char* structure);


PRIVATE char *   backtrack_CXS(int i, int j, const int** access_s1, const int** access_s2, const char* structure, int *E);
PRIVATE void find_max_CXS(const int *position, const int *position_j,const int delta, const int threshold, const int constthreshold, const int alignment_length, const char *s1, const char *s2, const int **access_s1, const int **access_s2, const int fast,const char* structure);
PRIVATE void     plot_max_CXS(const int max, const int max_pos, const int max_pos_j, const int alignment_length,const char *s1, const char *s2, const int **access_s1, const int **access_s2, const int fast, const char* structure);
PRIVATE duplexT duplexfold_C(const char *s1, const char *s2, const int extension_cost, const char* structure);
PRIVATE duplexT duplexfold_CXS(const char *s1, const char *s2,const int **access_s1, const int **access_s2, const int i_pos, const int j_pos, const int threshold ,const char* structure);


/*@unused@*/

#define MAXSECTORS      500     /* dimension for a backtrack array */
#define LOCALITY        0.      /* locality parameter for base-pairs */

#define MIN2(A, B)      ((A) < (B) ? (A) : (B))
#define MAX2(A, B)      ((A) > (B) ? (A) : (B))

PRIVATE vrna_param_t *P = NULL;
PRIVATE int   **c = NULL;/*, **in, **bx, **by;*/      /* energy array used in duplexfold */
/* PRIVATE int ****c_XS; */
PRIVATE int  **lc = NULL, **lin = NULL, **lbx = NULL, **lby = NULL, **linx = NULL, **liny = NULL;   /* energy array used in Lduplexfold
                                             this arrays contains only 3 columns
                                             In this way I reduce my memory use and
                                             I can make most of my computation and
                                             accession in the computer cash
                                             which is the main performance boost*/



/*PRIVATE int last_cell;                    this variable is the last_cell containing
                                            the information about the alignment
                                            useful only if there is an alignment
                                            which extends till the last nucleotide of
                                            the long sequence*/

PRIVATE short  *S1 = NULL, *SS1 = NULL, *S2 = NULL, *SS2 = NULL;/*contains the sequences*/
PRIVATE int   n1,n2;    /* sequence lengths */
PRIVATE int n3, n4; /*sequence length for the duplex*/;
PRIVATE int delay_free=0;


/*-----------------------------------------------------------------------duplexfold_XS---------------------------------------------------------------------------*/

PRIVATE duplexT duplexfold_CXS(const char *s1, const char *s2, const int **access_s1, const int **access_s2,
                               const int i_pos, const int j_pos, const int threshold, const char* structure) {
  int i, j,p,q,Emin=INF, l_min=0, k_min=0;
  char *struc;
  struc=NULL;
  duplexT mfe;
  vrna_md_t md;
  int bonus=-10000;
  n3 = (int) strlen(s1);
  n4 = (int) strlen(s2);

  int *previous_const;
  previous_const=(int *) vrna_alloc(sizeof(int) * (n4+1));
  j=0;
  previous_const[j]=1;
  int prev_temp = 1;
  while(j++<n4){
    if(structure[j-1]=='|'){
      previous_const[j]=prev_temp;
      prev_temp=j;
    }
    else{
      previous_const[j]=prev_temp;
    }
  }

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
  i=n3-1; j=2;
  type = pair[S1[i]][S2[j]];
  if(!type){
    printf("Error during initialization of the duplex in duplexfold_XS\n");
    mfe.structure=NULL;
    mfe.energy = INF;
    return mfe;
  }
  c[i][j] = P->DuplexInit + (structure[j-1]=='|' ? bonus : 0 ); /* check if first pair is constrained  */
  if(!(structure[j-2] == '|')){
    c[i][j]+=P->mismatchExt[rtype[type]][SS2[j-1]][SS1[i+1]];
  }
  else{
    c[i][j]+=P->dangle3[rtype[type]][SS1[i+1]];
  }
  if (type>2) c[i][j] += P->TerminalAU;
  for (k=i-1; k>0 ; k--) {
    c[k+1][0]=INF;
    for (l=j+1; l<=n4; l++) {
      c[k][l]=INF;
      int bonus_2 = (structure[l-1]=='|'? bonus : 0 ); /* check if position is constrained and prepare bonus accordingly */
      type2 = pair[S1[k]][S2[l]];
      if (!type2) continue;
      for (p=k+1; p< n3 && p<k+MAXLOOP-1; p++){
        for (q = l-1; q >= previous_const[l] && q > 1; q--) {
          if (p-k+l-q-2>MAXLOOP) break;
          type3=pair[S1[p]][S2[q]];
          if(!type3) continue;
          E = E_IntLoop(p-k-1, l-q-1, type2, rtype[type3],SS1[k+1], SS2[l-1], SS1[p-1], SS2[q+1],P) + bonus_2;
          c[k][l] = MIN2(c[k][l], c[p][q]+E);
        }
      }
      E = c[k][l];
      if (type2>2) E += P->TerminalAU;
      E+=access_s1[i-k+1][i_pos]+access_s2[l-1][j_pos+(l-1)-1];
      if (k>1 && l<n4 && !(structure[l]=='|') ){
        E+=P->mismatchExt[type2][SS1[k-1]][SS2[l+1]];
      }
      else if(k>1){
        E += P->dangle5[type2][SS1[k-1]];
      }
      else if(l<n4 && !(structure[l]=='|')){
        E += P->dangle3[type2][SS2[l+1]];
      }
      if (E<Emin) {
        Emin=E; k_min=k; l_min=l;
      }
    }
  }
  free(previous_const);
  if(Emin  > threshold){
    mfe.energy=INF;
    mfe.ddG=INF;
    mfe.structure=NULL;
    for (i=0; i<=n3; i++) free(c[i]);
    free(c);
    free(S1); free(S2); free(SS1); free(SS2);
    return mfe;
  } else{
    struc = backtrack_CXS(k_min, l_min, access_s1, access_s2,structure,&Emin);
  }


  /* lets take care of the dangles */
  /* find best combination  */
  int dx_5, dx_3, dy_5, dy_3,dGx,dGy,bonus_x;
  dx_5=0; dx_3=0; dy_5=0; dy_3=0;dGx=0;dGy=0;bonus_x=0;
  dGx = access_s1[i-k_min+1][i_pos];dx_3=0; dx_5=0;bonus_x=0;
  dGy = access_s2[l_min-j+1][j_pos + (l_min-1)-1];
  mfe.tb=i_pos -9 - i + k_min -1 -dx_5;
  mfe.te=i_pos -9 -1 + dx_3;
  mfe.qb=j_pos -9 -1 - dy_5;
  mfe.qe=j_pos + l_min -3 -9 + dy_3;
  mfe.ddG=(double) Emin * 0.01;
  mfe.dG1=(double) dGx*0.01 ;
  mfe.dG2=(double) dGy*0.01 ;
  /* mfe.energy += bonus_y + bonus_x; */
  mfe.energy= mfe.ddG - mfe.dG1 - mfe.dG2;

  mfe.structure = struc;
  for (i=0; i<=n3; i++) free(c[i]);
  free(c);
  free(S1); free(S2); free(SS1); free(SS2);
  return mfe;
}


PRIVATE char *backtrack_CXS (int i, int j, const int **access_s1,const int **access_s2,const char* structure, int *Emin ) {
  /* backtrack structure going backwards from i, and forwards from j
     return structure in bracket notation with & as separator */
  int k, l, type, type2, E, traced, i0, j0;
  char *st1, *st2, *struc;
  int *previous_const;
  int bonus=-10000;
  previous_const=(int *) vrna_alloc(sizeof(int) * (n4+1));
  int j_temp=0;
  previous_const[j_temp]=1;
  int prev_temp = 1;
  while(j_temp++<n4){
    if(structure[j_temp-1]=='|'){
      previous_const[j_temp]=prev_temp;
      prev_temp=j_temp;
    }
    else{
      previous_const[j_temp]=prev_temp;
    }
  }
  st1 = (char *) vrna_alloc(sizeof(char)*(n3+1));
  st2 = (char *) vrna_alloc(sizeof(char)*(n4+1));
  i0=i;/*MAX2(i-1,1);*/j0=j;/*MIN2(j+1,n4);*/
  while (i<=n3-1 && j>=2) {
    int bonus_2 = (structure[j-1]== '|'? bonus: 0);
    E = c[i][j]; traced=0;
    st1[i-1] = '(';
    st2[j-1] = ')';
    type = pair[S1[i]][S2[j]];
    if (!type) vrna_message_error("backtrack failed in fold duplex bli");
    for (k=i+1; k<=n3 && k>i-MAXLOOP-2; k++) {
      for (l=j-1; l >= previous_const[j] && l>=1; l--) {
        int LE;
        if (i-k+l-j-2>MAXLOOP) break;
        type2 = pair[S1[k]][S2[l]];
        if (!type2) continue;
        LE = E_IntLoop(k-i-1, j-l-1, type, rtype[type2], SS1[i+1], SS2[j-1], SS1[k-1], SS2[l+1],P) + bonus_2;
        if (E == c[k][l]+LE) {
          *Emin-=bonus_2;
          traced=1;
          i=k; j=l;
          break;
        }
      }
      if (traced) break;
    }
    if (!traced) {
      if(i<n3 && j>1 && !(structure[j-2]=='|')){
        E -= P->mismatchExt[rtype[type]][SS2[j-1]][SS1[i+1]];
      }
      else if (i<n3){
        E -= P->dangle3[rtype[type]][SS1[i+1]];/* +access_s1[1][i+1]; */
      }
      else if (j>1){
        E -= (!(structure[j-2]=='|') ? P->dangle5[rtype[type]][SS2[j-1]] : 0);/* +access_s2[1][j+1]; */
      }
      if (type>2) E -= P->TerminalAU;

      /* break; */
      if (E != P->DuplexInit + bonus_2)  {
        vrna_message_error("backtrack failed in fold duplex bal");
      } else {
        *Emin-=bonus_2;
        break;
      }
    }
  }
  /* if (i<n3)  i++; */
  /* if (j>1)   j--; */
  struc = (char *) vrna_alloc(i-i0+1+j0-j+1+2);
  for (k=MAX2(i0,1); k<=i; k++) if (!st1[k-1]) st1[k-1] = '.';
  for (k=j; k<=j0; k++) if (!st2[k-1]) st2[k-1] = '.';
  strcpy(struc, st1+MAX2(i0-1,0)); strcat(struc, "&");
  strcat(struc, st2+j-1);
  free(st1); free(st2);free(previous_const);
  return struc;
}


duplexT** Lduplexfold_CXS(const char *s1, const char *s2, const int **access_s1, const int **access_s2, const int threshold, const int alignment_length, const int delta, const int fast, const char* structure,const int il_a, const int il_b, const int b_a, const int b_b)/* , const int target_dead, const int query_dead) */
{

  int i, j;
  int bopen=b_b;
  int bext=b_a;
  int iopen=il_b;
  int iext_s=2*il_a;/* iext_s 2 nt nucleotide extension of interior loop, on i and j side */
  int iext_ass=50+il_a;/* iext_ass assymetric extension of interior loop, either on i or on j side. */
  int min_colonne=INF; /* enthaelt das maximum einer kolonne */
  int i_length;
  int max_pos;/* get position of the best hit */
  int max_pos_j;
  /* int temp; */
  int min_j_colonne;
  int max=INF;
  int bonus=-10000;
  int constthreshold=0; /* minimal threshold corresponding to a structure complying to all constraints */
  int maxPenalty[4];
  vrna_md_t   md;

  i=0;
  while(structure[i]!='\0'){
    if(structure[i]=='|') constthreshold+=bonus;
    i++;
  }
  int *position; /* contains the position of the hits with energy > E */
  int *position_j;
  n1 = (int) strlen(s1);
  n2 = (int) strlen(s2);
  position = (int *) vrna_alloc((delta+n1+3+delta) * sizeof(int));
  position_j= (int *) vrna_alloc((delta+n1+3+delta) * sizeof(int));

  set_model_details(&md);

  if ((!P) || (fabs(P->temperature - temperature)>1e-6)){
    update_dfold_params();
    if(P)
      free(P);
    P = vrna_params(&md);
    make_pair_matrix();
  }

  encode_seqs(s1,s2);

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
  for(j=n2;j>=0;j--) {
    lbx[0][j]=lbx[1][j]=lbx[2][j]=lbx[3][j]    = lbx[4][j] =INF;
    lin[0][j]=lin[1][j]=lin[2][j]=lin[3][j]    = lin[4][j] =INF;
    lc[0][j] =lc[1][j] =lc[2][j] = lc[3][j]    =  lc[4][j] =INF;
    lby[0][j]=lby[1][j]=lby[2][j]=lby[3][j]    = lby[4][j] =INF;
    liny[0][j]=liny[1][j]=liny[2][j]=liny[3][j]=liny[4][j]=INF;
    linx[0][j]=linx[1][j]=linx[2][j]=linx[3][j]=linx[4][j]=INF;
  }

  i=10 /*target_dead*/; /* start from 2 (        i=4) because no structure allowed to begin with a single base pair */
  i_length= n1 - 9  /*- target_dead*/ ;
  while(i < i_length) {
    int idx=i%5;
    int idx_1=(i-1)%5;
    int idx_2=(i-2)%5;
    int idx_3=(i-3)%5;
    int idx_4=(i-4)%5;
    int di1,di2,di3,di4;
    di1 = access_s1[5][i]   - access_s1[4][i-1];
    di2 = access_s1[5][i-1] - access_s1[4][i-2] + di1;
    di3 = access_s1[5][i-2] - access_s1[4][i-3] + di2;
    di4 = access_s1[5][i-3] - access_s1[4][i-4] + di3;
    di1=MIN2(di1,maxPenalty[0]);
    di2=MIN2(di2,maxPenalty[1]);
    di3=MIN2(di3,maxPenalty[2]);
    di4=MIN2(di4,maxPenalty[3]);
    j=n2 - 9 /*- query_dead*/; /* start from n2-1 because no structure allow to begin with a single base pair  */
    while (--j > 9/*query_dead - 1*/) {
      /* ----------------------------------------------------------update lin lbx lby matrix */
      int bonus_2 = (structure[j-1]=='|' ? bonus :0 );
      int dj1,dj2,dj3,dj4;
      dj1 = access_s2[5][j+4] - access_s2[4][j+4];
      dj2 = access_s2[5][j+5] - access_s2[4][j+5] + dj1;
      dj3 = access_s2[5][j+6] - access_s2[4][j+6] + dj2;
      dj4 = access_s2[5][j+7] - access_s2[4][j+7] + dj3;
      dj1=MIN2(dj1,maxPenalty[0]);
      dj2=MIN2(dj2,maxPenalty[1]);
      dj3=MIN2(dj3,maxPenalty[2]);
      dj4=MIN2(dj4,maxPenalty[3]);
      int type2, type,temp;
      type  = pair[S1[i]][S2[j]];
      lc[idx][j]= type ? P->DuplexInit + bonus_2 : INF;
      if(!bonus_2){
        type2=pair[S2[j+1]][S1[i-1]];
        lin[idx][j]=MIN2(lc[idx_1][j+1]+P->mismatchI[type2][SS2[j]][SS1[i]]+di1+dj1+iopen+iext_s,lin[idx_1][j]+iext_ass + di1);
        lin[idx][j]=MIN2(lin[idx][j],lin[idx][j+1]+iext_ass + dj1);
        lin[idx][j]=MIN2(lin[idx][j],lin[idx_1][j+1]+iext_s + di1 + dj1);
        linx[idx][j]=MIN2(lc[idx_1][j+1]+P->mismatch1nI[type2][SS2[j]][SS1[i]]+di1+dj1+iopen+iext_s,linx[idx_1][j]+iext_ass + di1);
        liny[idx][j]=MIN2(lc[idx_1][j+1]+P->mismatch1nI[type2][SS2[j]][SS1[i]]+di1+dj1+iopen+iext_s,liny[idx][j+1]+iext_ass + dj1);
        type2=pair[S2[j+1]][S1[i]];
        lby[idx][j]=MIN2(lby[idx][j+1]+bext + dj1 ,
                         lc[idx][j+1]+bopen+bext+(type2>2?P->TerminalAU:0)+dj1);
      }
      else{
        lin[idx][j] = lby[idx][j] = linx[idx][j]= liny[idx][j]=INF; /* all loop containing "|" are rejected */
      }
      type2=pair[S2[j]][S1[i-1]];
      lbx[idx][j]=MIN2(lbx[idx_1][j]+bext + di1, lc[idx_1][j]+bopen+bext+(type2>2?P->TerminalAU:0) + di1);
      /* --------------------------------------------------------------- end update recursion */
      if(!type){continue;}
      if(!(structure[j]=='|')){
        lc[idx][j]+=P->mismatchExt[type][SS1[i-1]][SS2[j+1]];
      }
      else{
        lc[idx][j]+=P->dangle5[type][SS1[i-1]];
      }
      lc[idx][j]+=(type>2?P->TerminalAU:0);
      /* type > 2 -> no GC or CG pair */
      /* ------------------------------------------------------------------update c  matrix  */
      /*  Be careful, no lc may come from a region where a "|" is in a loop, avoided in lin = lby = INF ... jedoch fuer klein loops muss man aufpassen .. */
      if((type2=pair[S1[i-1]][S2[j+1]]))
        lc[idx][j]=MIN2(lc[idx_1][j+1]+E_IntLoop(0,0,type2, rtype[type],SS1[i], SS2[j], SS1[i-1], SS2[j+1], P)+di1+dj1, lc[idx][j]); /* 0x0+1x1 */
      if((type2=pair[S1[i-2]][S2[j+1]]))
        lc[idx][j]=MIN2(lc[idx_2][j+1]+E_IntLoop(1,0,type2, rtype[type],SS1[i-1], SS2[j], SS1[i-1], SS2[j+1], P)+di2+dj1,lc[idx][j]);/* 0x1 +1x1 */
      /* kleine loops checks wird in den folgenden if test gemacht. */
      if(!(structure[j]=='|')){
        if((type2=pair[S1[i-1]][S2[j+2]]))
          lc[idx][j]=MIN2(lc[idx_1][j+2]+E_IntLoop(0,1,type2, rtype[type],SS1[i], SS2[j+1], SS1[i-1], SS2[j+1], P)+di1+dj2,lc[idx][j]);/* 1x0 + 1x1 */
        if((type2=pair[S1[i-2]][S2[j+2]]))
          lc[idx][j]=MIN2(lc[idx_2][j+2]+E_IntLoop(1,1,type2, rtype[type],SS1[i-1], SS2[j+1], SS1[i-1], SS2[j+1], P)+di2+dj2, lc[idx][j]); /*  1x1 +1x1 */
        if((type2 = pair[S1[i-3]][S2[j+2]]))
          lc[idx][j]=MIN2(lc[idx_3][j+2]+E_IntLoop(2,1,type2, rtype[type],SS1[i-2], SS2[j+1], SS1[i-1], SS2[j+1], P)+di3+dj2, lc[idx][j]); /*  2x1 +1x1 */
        if(!(structure[j+1]=='|')){
          if((type2 = pair[S1[i-3]][S2[j+3]]))
            lc[idx][j]=MIN2(lc[idx_3][j+3]+E_IntLoop(2,2,type2, rtype[type],SS1[i-2], SS2[j+2], SS1[i-1], SS2[j+1], P)+di3+dj3,lc[idx][j]);/* 2x2 + 1x1 */
          if((type2 = pair[S1[i-2]][S2[j+3]]))
            lc[idx][j]=MIN2(lc[idx_2][j+3]+E_IntLoop(1,2,type2, rtype[type],SS1[i-1], SS2[j+2], SS1[i-1], SS2[j+1], P)+di2+dj3, lc[idx][j]);/*  1x2 +1x1 */
          if((type2 = pair[S1[i-4]][S2[j+3]]))
            lc[idx][j]=MIN2(lc[idx_4][j+3]+E_IntLoop(3,2,type2, rtype[type],SS1[i-3], SS2[j+2], SS1[i-1], SS2[j+1], P)+di4+dj3, lc[idx][j]);
          if(!(structure[j+2]=='|')){
            if((type2 = pair[S1[i-3]][S2[j+4]]))
              lc[idx][j]=MIN2(lc[idx_3][j+4]+E_IntLoop(2,3,type2, rtype[type],SS1[i-2], SS2[j+3], SS1[i-1], SS2[j+1], P)+di3+dj4, lc[idx][j]);
          }
        }
      }
      /* internal->stack  */
      lc[idx][j]=MIN2(lin[idx_3][j+3]+P->mismatchI[rtype[type]][SS1[i-1]][SS2[j+1]]+di3+dj3+2*iext_s, lc[idx][j]);
      lc[idx][j]=MIN2(lin[idx_4][j+2]+P->mismatchI[rtype[type]][SS1[i-1]][SS2[j+1]]+iext_s+2*iext_ass+di4+dj2, lc[idx][j]);
      lc[idx][j]=MIN2(lin[idx_2][j+4]+P->mismatchI[rtype[type]][SS1[i-1]][SS2[j+1]]+iext_s+2*iext_ass+di2+dj4, lc[idx][j]);
      lc[idx][j]=MIN2(linx[idx_3][j+1]+P->mismatch1nI[rtype[type]][SS1[i-1]][SS2[j+1]]+iext_ass+iext_ass+di3+dj1, lc[idx][j]);
      lc[idx][j]=MIN2(liny[idx_1][j+3]+P->mismatch1nI[rtype[type]][SS1[i-1]][SS2[j+1]]+iext_ass+iext_ass+dj3+di1, lc[idx][j]);
      /* bulge -> stack */
      int bAU;
      bAU=(type>2?P->TerminalAU:0);
      lc[idx][j]=MIN2(lbx[idx_2][j+1]+di2+dj1+bext+bAU, lc[idx][j]);
      /* min2=by[i][j+1]; */
      lc[idx][j]=MIN2(lby[idx_1][j+2]+di1+dj2+bext+bAU, lc[idx][j]);
      lc[idx][j]+=bonus_2;
      /* if(j<=const5end){ */
      temp=min_colonne;
      min_colonne=MIN2(lc[idx][j]+(type>2?P->TerminalAU:0)+
                       (!(structure[j-2]=='|') ?
                        P->mismatchExt[rtype[type]][SS2[j-1]][SS1[i+1]] : P->dangle3[rtype[type]][SS1[i+1]]),
                       min_colonne);
      if(temp>min_colonne){
        min_j_colonne=j;
      }
      /* } */
      /* ---------------------------------------------------------------------end update */
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
  /* printf("MAX :%d ", max); */
  free(S1); free(S2); free(SS1); free(SS2);
  if(max<threshold+constthreshold){
    find_max_CXS(position, position_j, delta, threshold+constthreshold, constthreshold, alignment_length, s1, s2, access_s1, access_s2, fast, structure);
  }
  if(max<constthreshold){
    plot_max_CXS(max, max_pos, max_pos_j,alignment_length, s1, s2, access_s1, access_s2,fast,structure);
  }
  for (i=0; i<=4; i++) {free(lc[i]);free(lin[i]);free(lbx[i]);free(lby[i]);free(linx[i]);free(liny[i]);}
  /* free(lc[0]);free(lin[0]);free(lbx[0]);free(lby[0]); */
  free(lc);free(lin);free(lbx);free(lby);free(linx);free(liny);
  free(position);
  free(position_j);
  return NULL;
}

PRIVATE void find_max_CXS(const int *position, const int *position_j,const int delta, const int threshold, const int constthreshold, const int alignment_length, const char *s1, const char *s2, const int **access_s1, const int **access_s2, const int fast, const char* structure){
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
        /* int begin_t=MAX2(9, pos-alignment_length); */
        /* int end_t  =MIN2(n1-10, pos); */
        /* int begin_q=MAX2(9, max_pos_j-2); */
        /* int end_q  =MIN2(n2-10, max_pos_j+alignment_length-2); */
        int begin_t=MAX2(9,pos-alignment_length);
        int end_t  =pos;
        int begin_q=max_pos_j-2;
        int end_q  =MIN2(n2-9,max_pos_j+alignment_length-2);
        char *s3 = (char*) vrna_alloc(sizeof(char)*(end_t - begin_t +2));
        char *s4 = (char*) vrna_alloc(sizeof(char)*(end_q - begin_q +2));
        char *local_structure = (char*) vrna_alloc(sizeof(char) * ( end_q - begin_q +2));
        strncpy(s3, (s1+begin_t),  end_t - begin_t+1);
        strncpy(s4, (s2+begin_q) , end_q - begin_q+1 );
        strncpy(local_structure, (structure+begin_q), end_q - begin_q +1);
        s3[end_t -begin_t +1 ]='\0';
        s4[end_q -begin_q +1 ]='\0';
        local_structure[end_q - begin_q +1]='\0';
        duplexT test;
        test = duplexfold_CXS(s3,s4,access_s1,access_s2,pos, max_pos_j,threshold,local_structure);
        if(test.energy * 100 < (threshold - constthreshold)){
          int l1=strchr(test.structure, '&')-test.structure;
          int dL = strrchr(structure,'|') - strchr(structure,'|');
          dL+=1;
          if(dL <=  strlen(test.structure)-l1-1){
            printf("%s %3d,%-3d : %3d,%-3d (%5.2f = %5.2f + %5.2f + %5.2f)\n", test.structure,
                   test.tb,test.te,test.qb,test.qe, test.ddG, test.energy, test.dG1, test.dG2);
            pos=MAX2(10,pos+temp_min-delta);
          }
        }
        free(s3);free(s4);
        free(test.structure);
        free(local_structure);
      }
    }
  }
}


PRIVATE void plot_max_CXS(const int max, const int max_pos, const int max_pos_j, const int alignment_length, const char *s1, const char *s2, const int ** access_s1, const int ** access_s2, const int fast, const char* structure)
{
  if(fast==1){
    printf("target upper bound %d: query lower bound %d (%5.2f)\n", max_pos-3, max_pos_j, ((double)max)/100);
  }
  else{
    int begin_t=MAX2(9,max_pos-alignment_length);
    int end_t  =max_pos;
    int begin_q=max_pos_j-2;
    int end_q  =MIN2(n2-9,max_pos_j+alignment_length-2);
    char *s3 = (char*) vrna_alloc(sizeof(char)*(end_t - begin_t +2));
    char *s4 = (char*) vrna_alloc(sizeof(char)*(end_q - begin_q +2));
    char *local_structure = (char*) vrna_alloc(sizeof(char)*(end_q - begin_q +2));
    strncpy(s3, (s1+begin_t),  end_t - begin_t+1);
    strncpy(s4, (s2+begin_q) , end_q - begin_q+1 );
    strncpy(local_structure, (structure+begin_q) , end_q - begin_q +1 );
    s3[end_t -begin_t +1 ]='\0';
    s4[end_q -begin_q +1 ]='\0';
    local_structure[end_q - begin_q +1]='\0';
    duplexT test;
    test = duplexfold_CXS(s3,s4,access_s1,access_s2,max_pos, max_pos_j,INF,local_structure);
    int l1=  strchr(test.structure, '&')-test.structure;
    int dL = strrchr(structure,'|') - strchr(structure,'|');
    dL+=1;
    if(dL<=strlen(test.structure)-l1-1){
      printf("%s %3d,%-3d : %3d,%-3d (%5.2f = %5.2f + %5.2f + %5.2f)\n", test.structure,
             test.tb,test.te,test.qb,test.qe, test.ddG, test.energy, test.dG1, test.dG2);
    }
    free(s3);free(s4);free(test.structure);free(local_structure);

  }
}


/*---------------------------------------------------------duplexfold----------------------------------------------------------------------------------*/


PRIVATE duplexT duplexfold_C(const char *s1, const char *s2, const int extension_cost, const char* structure ) {
  int i, j, l1, Emin=INF, i_min=0, j_min=0;
  char *struc;
  duplexT mfe;
  vrna_md_t   md;
  int bonus=-10000;
  int *previous_const; /* for each "|" constraint returns the position of the next "|" constraint */

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
  previous_const=(int *) vrna_alloc(sizeof(int) * (n4+1));
  j=n4+1;
  previous_const[j-1]=n4;
  int prev_temp = n4;
  while(--j){
    if(structure[j-1]=='|'){
      previous_const[j-1]=prev_temp;
      prev_temp=j;
    }
    else{
      previous_const[j-1]=prev_temp;
    }
  }
  c = (int **) vrna_alloc(sizeof(int *) * (n3+1));
  for (i=0; i<=n3; i++) c[i] = (int *) vrna_alloc(sizeof(int) * (n4+1));
  encode_seqs(s1, s2);
  for (i=1; i<=n3; i++) {
    for (j=n4; j>0; j--) {
      int type, type2, E, k,l;
      int bonus_2 = (structure[j-1]=='|'? bonus: 0);
      type = pair[S1[i]][S2[j]];
      c[i][j] = type ? P->DuplexInit +2 * extension_cost + bonus_2: INF;
      if(!type){ continue;}
      if(j<n4 && i>1 && !(structure[j]=='|') ) {
        c[i][j]+=P->mismatchExt[type][SS1[i-1]][SS2[j+1]]+2*extension_cost;
      }
      else if(i>1){
        c[i][j] += P->dangle5[type][SS1[i-1]]+ extension_cost;
      }
      else if(j<n4 && !(structure[j]=='|')){
        c[i][j] += P->dangle3[type][SS2[j+1]]+ extension_cost;
      }
      if (type>2) c[i][j] += P->TerminalAU;
      for (k=i-1; k>0 && k>i-MAXLOOP-2; k--) {
        for (l=j+1; l<=previous_const[j]; l++) {
          if (i-k+l-j-2>MAXLOOP) break;
          type2 = pair[S1[k]][S2[l]];
          if (!type2) continue;
          E = E_IntLoop(i-k-1, l-j-1, type2, rtype[type],
                        SS1[k+1], SS2[l-1], SS1[i-1], SS2[j+1],P)+(i-k+l-j)*extension_cost + bonus_2;
          c[i][j] = MIN2(c[i][j], c[k][l]+E);
        }
      }
      E = c[i][j];
      if(i<n3 && j>1 && !(structure[j-2]=='|')){
        E+= P->mismatchExt[rtype[type]][SS2[j-1]][SS1[i+1]]+2*extension_cost;
      }
      else if (i<n3){
        E += P->dangle3[rtype[type]][SS1[i+1]]+extension_cost;
      }
      else if (j>1 && !(structure[j-2]=='|')){
        E += P->dangle5[rtype[type]][SS2[j-1]]+extension_cost;
      }
      if (type>2) E += P->TerminalAU;

      if (E<Emin) {
        Emin=E; i_min=i; j_min=j;
      }
    }
  }
  struc = backtrack_C(i_min, j_min, extension_cost,structure,&Emin);
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
  free(previous_const);
  if (!delay_free) {
    for (i=0; i<=n3; i++) free(c[i]);

    free(c);
    free(S1); free(S2); free(SS1); free(SS2);
  }
  return mfe;
}

PRIVATE char *backtrack_C(int i, int j, const int extension_cost, const char* structure, int *Emin) {
  /* backtrack structure going backwards from i, and forwards from j
     return structure in bracket notation with & as separator */
  int k, l, type, type2, E, traced, i0, j0, *previous_const;
  char *st1, *st2, *struc;
  int bonus=-10000;
  previous_const=(int *) vrna_alloc(sizeof(int) * (n4+1)); /* encodes the position of the constraints */
  int j_temp=n4+1;
  previous_const[j_temp-1]=n4;
  int prev_temp = n4;
  while(--j_temp){
    if(structure[j_temp-1]=='|'){
      previous_const[j_temp-1]=prev_temp;
      prev_temp=j_temp;
    }
    else{
      previous_const[j_temp-1]=prev_temp;
    }
  }
  st1 = (char *) vrna_alloc(sizeof(char)*(n3+1));
  st2 = (char *) vrna_alloc(sizeof(char)*(n4+1));
  i0=MIN2(i+1,n3); j0=MAX2(j-1,1);
  while (i>0 && j<=n4) {
    int bonus_2 = (structure[j-1]== '|'? bonus: 0);
    E = c[i][j]; traced=0;
    st1[i-1] = '(';
    st2[j-1] = ')';
    type = pair[S1[i]][S2[j]];
    if (!type) vrna_message_error("backtrack failed in fold duplex a");
    for (k=i-1; k>0 && k>i-MAXLOOP-2; k--) {
      for (l=j+1; l<=previous_const[j]; l++) {
        int LE;
        if (i-k+l-j-2>MAXLOOP) break;
        type2 = pair[S1[k]][S2[l]];
        if (!type2) continue;
        LE = E_IntLoop(i-k-1, l-j-1, type2, rtype[type],
                       SS1[k+1], SS2[l-1], SS1[i-1], SS2[j+1],P)+(i-k+l-j)*extension_cost + bonus_2;
        if (E == c[k][l]+LE) {
          *Emin-=bonus_2;
          traced=1;
          i=k; j=l;
          break;
        }
      }
      if (traced) break;
    }
    if (!traced) {

      if (i>1 && j<n4 && !(structure[j]=='|')){
        E -=P->mismatchExt[type][SS1[i-1]][SS2[j+1]]+2*extension_cost;
      }
      else if(i>1){
        E -= P->dangle5[type][SS1[i-1]]+extension_cost;
      }
      else if (j<n4 && !(structure[j]=='|')){
        E -= P->dangle3[type][SS2[j+1]]+ extension_cost;
      }
      /* if (j<n4) E -= P->dangle3[type][SS2[j+1]]+extension_cost; */
      if (type>2) E -= P->TerminalAU;
      if (E != P->DuplexInit+2*extension_cost + bonus_2) {
        vrna_message_error("backtrack failed in fold duplex b");
      } else {
        *Emin-=bonus_2;
        break;
      }
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
  free(st1); free(st2);
  free(previous_const);
  return struc;
}




duplexT ** Lduplexfold_C(const char *s1, const char *s2, const int threshold, const int extension_cost, const int alignment_length, const int delta, const int fast, const char* structure, const int il_a, const int il_b, const int b_a, const int b_b)
{
  /* duplexT test = duplexfold_C(s1, s2, extension_cost,structure); */

  int i, j;
  int bopen=b_b;
  int bext=b_a+extension_cost;
  int iopen=il_b;
  int iext_s=2*(il_a+extension_cost);/* iext_s 2 nt nucleotide extension of interior loop, on i and j side */
  int iext_ass=50+il_a+extension_cost;/* iext_ass assymetric extension of interior loop, either on i or on j side. */
  int min_colonne=INF; /* enthaelt das maximum einer kolonne */
  int i_length;
  int max_pos;/* get position of the best hit */
  int max_pos_j=10;
  int temp;
  int min_j_colonne=11;
  int max=INF;
  int bonus = -10000;
  int constthreshold=0; /* minimal threshold corresponding to a structure complying to all constraints */
  i=0;
  while(structure[i]!='\0'){
    if(structure[i]=='|') constthreshold+=bonus;
    i++;
  }
  /* FOLLOWING NEXT 4 LINE DEFINES AN ARRAY CONTAINING POSITION OF THE SUBOPT IN S1 */
  /* int nsubopt=10;  */ /* total number of subopt */
  int *position; /* contains the position of the hits with energy > E */
  int *position_j;
  /*   int const5end; */ /* position of the 5'most constraint. Only interaction reaching this position are taken into account. */
  /* const5end = strchr(structure,'|') - structure; */
  /* const5end++; */
  n1 = (int) strlen(s1);
  n2 = (int) strlen(s2);
  /* delta_check is the minimal distance allowed for two hits to be accepted */
  /* if both hits are closer, reject the smaller ( in term of position)  hits  */
  position = (int *) vrna_alloc((delta+n1+3+delta) * sizeof(int));
  position_j= (int *) vrna_alloc((delta+n1+3+delta) * sizeof(int));
  /* i want to implement a function that, given a position in a long sequence and a small sequence, */
  /* duplexfold them at this position and report the result at the command line */
  /* for this i first need to rewrite backtrack in order to remove the printf functio */
  /* END OF DEFINITION FOR NEEDED SUBOPT DATA  */

  if ((!P) || (fabs(P->temperature - temperature)>1e-6))
  update_dfold_params();

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
  for(j=n2;j>=0;j--) {
    lbx[0][j]=lbx[1][j]=lbx[2][j]=lbx[3][j]    = lbx[4][j] =INF;
    lin[0][j]=lin[1][j]=lin[2][j]=lin[3][j]    = lin[4][j] =INF;
    lc[0][j] =lc[1][j] =lc[2][j] = lc[3][j]    =  lc[4][j] =INF;
    lby[0][j]=lby[1][j]=lby[2][j]=lby[3][j]    = lby[4][j] =INF;
    liny[0][j]=liny[1][j]=liny[2][j]=liny[3][j]=liny[4][j]=INF;
    linx[0][j]=linx[1][j]=linx[2][j]=linx[3][j]=linx[4][j]=INF;
  }
  encode_seqs(s1,s2);
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
      int bonus_2 = (structure[j-1]=='|' ? bonus : 0) ;
      int type, type2;
      type = pair[S1[i]][S2[j]];
      lc[idx][j]=type ? P->DuplexInit + 2*extension_cost + bonus_2 : INF; /* to avoid that previous value influence result should actually not be erforderlich */
      if(!bonus_2){
        type2=pair[S2[j+1]][S1[i-1]];
        lin[idx][j]=MIN2(lc[idx_1][j+1]+P->mismatchI[type2][SS2[j]][SS1[i]]+iopen+iext_s, lin[idx_1][j]+iext_ass);
        lin[idx][j]=MIN2(lin[idx][j],lin[idx][j+1]+iext_ass);
        lin[idx][j]=MIN2(lin[idx][j],lin[idx_1][j+1]+iext_s);
        linx[idx][j]=MIN2(lc[idx_1][j+1]+P->mismatch1nI[type2][SS2[j]][SS1[i]]+iopen+iext_s,linx[idx_1][j]+iext_ass);
        liny[idx][j]=MIN2(lc[idx_1][j+1]+P->mismatch1nI[type2][SS2[j]][SS1[i]]+iopen+iext_s,liny[idx][j+1]+iext_ass);
        type2=pair[S2[j+1]][S1[i]];
        lby[idx][j]=MIN2(lby[idx][j+1]+bext, lc[idx][j+1]+bopen+bext+(type2>2?P->TerminalAU:0));
      }
      else{
        lin[idx][j] = lby[idx][j] = linx[idx][j]= liny[idx][j]=INF;
      }
      type2=pair[S2[j]][S1[i-1]];
      lbx[idx][j]=MIN2(lbx[idx_1][j]+bext, lc[idx_1][j]+bopen+bext+(type2>2?P->TerminalAU:0));
      /* --------------------------------------------------------------- end update recursion */
      if(!type){continue;}
      if(!(structure[j]=='|')){
        lc[idx][j]+=P->mismatchExt[type][SS1[i-1]][SS2[j+1]]+2*extension_cost;
      }
      else{
        lc[idx][j]+=P->dangle5[type][SS1[i-1]]+extension_cost;
      }
      lc[idx][j]+=(type>2?P->TerminalAU:0);
      /* type > 2 -> no GC or CG pair */
      /* ------------------------------------------------------------------update c  matrix  */
      /*  Be careful, no lc may come from a region where a "|" is in a loop, avoided in lin = lby = INF ... jedoch fuer klein loops muss man aufpassen .. */
      type2=pair[S1[i-1]][S2[j+1]];
      lc[idx][j]=MIN2(lc[idx_1][j+1]+E_IntLoop(0,0,type2, rtype[type],SS1[i], SS2[j], SS1[i-1], SS2[j+1], P)+2*extension_cost, lc[idx][j]);
      type2=pair[S1[i-2]][S2[j+1]];
      lc[idx][j]=MIN2(lc[idx_2][j+1]+E_IntLoop(1,0,type2, rtype[type],SS1[i-1], SS2[j], SS1[i-1], SS2[j+1], P)+3*extension_cost,lc[idx][j]);
      /* kleine loops checks wird in den folgenden if test gemacht. */
      if(!(structure[j]=='|')){
        type2=pair[S1[i-1]][S2[j+2]];
        lc[idx][j]=MIN2(lc[idx_1][j+2]+E_IntLoop(0,1,type2, rtype[type],SS1[i], SS2[j+1], SS1[i-1], SS2[j+1], P)+3*extension_cost,lc[idx][j]);
        type2=pair[S1[i-2]][S2[j+2]];
        lc[idx][j]=MIN2(lc[idx_2][j+2]+E_IntLoop(1,1,type2, rtype[type],SS1[i-1], SS2[j+1], SS1[i-1], SS2[j+1], P)+4*extension_cost, lc[idx][j]);
        type2 = pair[S1[i-3]][S2[j+2]];
        lc[idx][j]=MIN2(lc[idx_3][j+2]+E_IntLoop(2,1,type2, rtype[type],SS1[i-2], SS2[j+1], SS1[i-1], SS2[j+1], P)+5*extension_cost, lc[idx][j]);
        if(!(structure[j+1]=='|')){
          type2 = pair[S1[i-3]][S2[j+3]];
          lc[idx][j]=MIN2(lc[idx_3][j+3]+E_IntLoop(2,2,type2, rtype[type],SS1[i-2], SS2[j+2], SS1[i-1], SS2[j+1], P)+6*extension_cost,lc[idx][j]);
          type2 = pair[S1[i-2]][S2[j+3]];
          lc[idx][j]=MIN2(lc[idx_2][j+3]+E_IntLoop(1,2,type2, rtype[type],SS1[i-1], SS2[j+2], SS1[i-1], SS2[j+1], P)+5*extension_cost, lc[idx][j]);
          type2 = pair[S1[i-4]][S2[j+3]];
          lc[idx][j]=MIN2(lc[idx_4][j+3]+E_IntLoop(3,2,type2, rtype[type],SS1[i-3], SS2[j+2], SS1[i-1], SS2[j+1], P)+7*extension_cost, lc[idx][j]);
          if(!(structure[j+2]=='|')){
            type2 = pair[S1[i-3]][S2[j+4]];
            lc[idx][j]=MIN2(lc[idx_3][j+4]+E_IntLoop(2,3,type2, rtype[type],SS1[i-2], SS2[j+3], SS1[i-1], SS2[j+1], P)+7*extension_cost, lc[idx][j]);
          }
        }
      }
      /* internal->stack  */
      lc[idx][j]=MIN2(lin[idx_3][j+3]+P->mismatchI[rtype[type]][SS1[i-1]][SS2[j+1]]+2*extension_cost+2*iext_s, lc[idx][j]);
      lc[idx][j]=MIN2(lin[idx_4][j+2]+P->mismatchI[rtype[type]][SS1[i-1]][SS2[j+1]]+iext_s+2*iext_ass+2*extension_cost, lc[idx][j]);
      lc[idx][j]=MIN2(lin[idx_2][j+4]+P->mismatchI[rtype[type]][SS1[i-1]][SS2[j+1]]+iext_s+2*iext_ass+2*extension_cost, lc[idx][j]);
      lc[idx][j]=MIN2(linx[idx_3][j+1]+P->mismatch1nI[rtype[type]][SS1[i-1]][SS2[j+1]]+iext_ass+iext_ass+2*extension_cost, lc[idx][j]);
      lc[idx][j]=MIN2(liny[idx_1][j+3]+P->mismatch1nI[rtype[type]][SS1[i-1]][SS2[j+1]]+iext_ass+iext_ass+2*extension_cost, lc[idx][j]);
      /* bulge -> stack */
      int bAU;
      bAU=(type>2?P->TerminalAU:0);
      lc[idx][j]=MIN2(lbx[idx_2][j+1]+2*extension_cost+bext+bAU, lc[idx][j]);
      /* min2=by[i][j+1]; */
      lc[idx][j]=MIN2(lby[idx_1][j+2]+2*extension_cost+bext+bAU, lc[idx][j]);
      lc[idx][j]+=bonus_2;
      /*       if(j<=const5end){ */
      temp=min_colonne;
      min_colonne=MIN2(lc[idx][j]+(type>2?P->TerminalAU:0)+
                       (!(structure[j-2]=='|') ?
                        P->mismatchExt[rtype[type]][SS2[j-1]][SS1[i+1]]+2*extension_cost :
                        P->dangle3[rtype[type]][SS1[i+1]]+extension_cost),
                       min_colonne);
      if(temp>min_colonne){
        min_j_colonne=j;
        /*         } */
      }
      /* ---------------------------------------------------------------------end update       */
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
  free(S1); free(S2); free(SS1); free(SS2);
  /* printf("MAX: %d",max); */
  if(max<threshold+constthreshold){
    find_max_C(position, position_j, delta, threshold+constthreshold, constthreshold, alignment_length, s1, s2, extension_cost, fast, structure);
  }
  if(max<constthreshold){
    plot_max_C(max, max_pos, max_pos_j,alignment_length, s1, s2, extension_cost,fast,structure);
  }
  for (i=0; i<=4; i++) {free(lc[i]);free(lin[i]);free(lbx[i]);free(lby[i]);free(linx[i]);free(liny[i]);}
  /*  free(lc[0]);free(lin[0]);free(lbx[0]);free(lby[0]); */
  free(lc);free(lin);free(lbx);free(lby);free(linx);free(liny);
  free(position);
  free(position_j);
  return NULL;
}


PRIVATE void find_max_C(const int *position, const int *position_j,const int delta, const int threshold, const int constthreshold, const int alignment_length, const char *s1, const char *s2, const int extension_cost, const int fast,const char* structure){
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
        pos=MAX2(10,pos-delta);
      }
    }
  }
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
        int begin_t=MAX2(11, pos-alignment_length+1);
        int end_t  =MIN2(n1-10, pos+1);
        int begin_q=MAX2(11, max_pos_j-1);
        int end_q  =MIN2(n2-10, max_pos_j+alignment_length-2);
        char *s3 = (char*) vrna_alloc(sizeof(char)*(end_t - begin_t +2));
        char *s4 = (char*) vrna_alloc(sizeof(char)*(end_q - begin_q +2));
        char *local_structure = (char*) vrna_alloc(sizeof(char)*(end_q - begin_q +2));
        strncpy(s3, (s1+begin_t-1),  end_t - begin_t +1);
        strncpy(s4, (s2+begin_q-1) , end_q - begin_q +1);
        strncpy(local_structure, (structure+begin_q-1) , end_q - begin_q +1 );
        s3[end_t -begin_t +1 ]='\0';
        s4[end_q -begin_q +1 ]='\0';
        local_structure[end_q - begin_q +1]='\0';
        duplexT test;
        test = duplexfold_C(s3, s4, extension_cost,local_structure);
        if(test.energy * 100 < (threshold-constthreshold)){
          int l1=strchr(test.structure, '&')-test.structure;
          int dL = strrchr(structure,'|') - strchr(structure,'|');
          dL+=1;
          if(dL <=  strlen(test.structure)-l1-1){
            printf("%s %3d,%-3d : %3d,%-3d (%5.2f)\n", test.structure,
                   begin_t-10+test.i-l1,
                   begin_t-10+test.i-1,
                   begin_q-10 + test.j-1 ,
                   (begin_q -11) + test.j + (int)strlen(test.structure)-l1-2,
                   test.energy);
            pos=MAX2(10,pos-delta);
          }
        }
        free(s3);free(s4);
        free(test.structure);
        free(local_structure);
      }
    }
  }
}
PRIVATE void plot_max_C(const int max, const int max_pos, const int max_pos_j, const int alignment_length, const char *s1, const char *s2, const int extension_cost, const int fast,const char* structure)
{
  if(fast==1){
    printf("target upper bound %d: query lower bound %d (%5.2f)\n", max_pos-10, max_pos_j-10, ((double)max)/100);
  }
  else{
    duplexT test;
    int begin_t=MAX2(11, max_pos-alignment_length+1);
    int end_t  =MIN2(n1-10, max_pos+1);
    int begin_q=MAX2(11, max_pos_j-1);
    int end_q  =MIN2(n2-10, max_pos_j+alignment_length-2);
    char *s3 = (char*) vrna_alloc(sizeof(char)*(end_t - begin_t +2));
    char *s4 = (char*) vrna_alloc(sizeof(char)*(end_q - begin_q +2));
    char *local_structure = (char*) vrna_alloc(sizeof(char)*(end_q - begin_q +2));
    strncpy(s3, (s1+begin_t-1),  end_t - begin_t + 1);
    strncpy(s4, (s2+begin_q-1) , end_q - begin_q +1 );
    strncpy(local_structure, (structure+begin_q-1) , end_q - begin_q +1 );
    s3[end_t -begin_t +1 ]='\0';
    s4[end_q -begin_q +1 ]='\0';
    local_structure[end_q - begin_q +1]='\0';
    test = duplexfold_C(s3, s4, extension_cost,local_structure);
    int l1=strchr(test.structure, '&')-test.structure;
    int dL = strrchr(structure,'|') - strchr(structure,'|');
    dL+=1;
    if(dL <=  strlen(test.structure)-l1-1){
      printf("%s %3d,%-3d : %3d,%-3d (%5.2f)\n", test.structure,
             begin_t-10+test.i-l1, begin_t-10+test.i-1, begin_q-10 +test.j-1 ,
             (begin_q -11) + test.j + (int)strlen(test.structure)-l1-2  , test.energy);
      free(s3);free(s4);free(test.structure);
    }
    free(local_structure);
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

/*---------------------------------------------------------------------------*/


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


