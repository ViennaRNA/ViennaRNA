/*
           compute potentially pseudoknotted structures of an RNA sequence
                             Ivo Hofacker
                          Vienna RNA package
*/

/*
  library containing the function used in PKplex
  it generates pseudoknotted structures by letting the sequence form a duplex structure with itself
*/

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <string.h>
#include <time.h>

#include "ViennaRNA/energy_par.h"
#include "ViennaRNA/fold_vars.h"
#include "ViennaRNA/utils.h"
#include "ViennaRNA/params.h"
#include "ViennaRNA/loop_energies.h"

#include "ViennaRNA/pair_mat.h"

#include "ViennaRNA/fold.h"
#include "ViennaRNA/PKplex.h"

#undef  MAXLOOP
#define MAXLOOP 10

PRIVATE void  update_dfold_params(void);
PRIVATE void  duplexfold_XS(const char *s1, int **access_s1, const int threshold, const int max_interaction_length);
PRIVATE char  *backtrack_XS(int kk, int ll, const int ii, const int jj, const int max_interaction_length);
PRIVATE void  make_ptypes(const char *structure);

PRIVATE vrna_param_t  *P = NULL;
PRIVATE int     ***c3 = NULL;      /* energy array used in duplexfold */
PRIVATE short   *S1 = NULL, *SS1 = NULL;
PRIVATE int     n1;
PRIVATE char    *ptype = NULL;     /* precomputed array of pair types */
PRIVATE int     *indx = NULL;      /* index for moving in the triangle matrices ptype[] */

PUBLIC  dupVar  *PlexHits = NULL;
PUBLIC  int     PlexHitsArrayLength = 100;
PUBLIC  int     NumberOfHits        = 0;
PUBLIC  int     verbose             = 0;

/*-----------------------------------------------------------------------duplexfold_XS---------------------------------------------------------------------------*/

PRIVATE void  duplexfold_XS(const char *s1,
                            int **access_s1,
                            const int threshold,
                            const int max_interaction_length){

  int i, j, k, l, p, q, Emin=INF, l_min=0, k_min=0, j_min=0;
  int type, type2, type3, E, tempK;
  char *struc;
  int length = (int) strlen(s1);
  struc=NULL;

  c3 = (int ***) vrna_alloc(sizeof(int **) * (length));
  for (i=0; i<length; i++){
    c3[i] = (int **) vrna_alloc(sizeof(int *) * max_interaction_length);
    for (j=0; j<max_interaction_length; j++) {
      c3[i][j]=(int *) vrna_alloc(sizeof(int) * max_interaction_length);
    }
  }

  i = length - 9;

  while( i-- > 11 ){
    Emin=INF;
    j_min=0;
    l_min=0;
    k_min=0;

    /* init all matrix elements to INF */
    for (j=0; j < length; j++){
      for(k=0;k<max_interaction_length;k++){
        for(l=0;l<max_interaction_length;l++){
          c3[j][k][l]=INF;
        }
      }
    }
    char string[10] = {'\0'};
    /* matrix starting values for (i,j)-basepairs */
    for(j=i+4; j<n1-10; j++) {
      type=ptype[indx[j]+i];
      if (type) {
        c3[j-11][max_interaction_length-1][0] = P->DuplexInit;
        c3[j-11][max_interaction_length-1][0] += E_Hairpin(j-i-1, type,  SS1[i+1], SS1[j-1], string, P);
/*        c3[j-11][max_interaction_length-1][0] += E_ExtLoop(type, SS1[i+1], SS1[j-1], P); */
/*           c3[j-11][max_interaction_length-1][0] += E_ExtLoop(rtype[type], SS1[j-1], SS1[i+1], P); */
       }
    }

    int i_pos_begin=MAX2(9, i-max_interaction_length); /* why 9 ??? */

    /* fill matrix */
    for(k=i-1; k>i_pos_begin; k--) {
      tempK=max_interaction_length-i+k-1;
      for(l = i + 5; l < n1 - 9; l++) { /* again, why 9 less then the sequence length ? */
        type2 = ptype[indx[l] + k];
        if(!type2) continue;
        for(p = k + 1; (p <= i) && (p <= k + MAXLOOP + 1); p++){
          for(q = l - 1; (q>=i+4) && (q >= l - MAXLOOP - 1); q--){
            if (p - k + l - q - 2 > MAXLOOP) break;
            type3 = ptype[indx[q] + p];
            if(!type3) continue;
            E = E_IntLoop(p - k - 1, l - q - 1, type2, rtype[type3], SS1[k + 1], SS1[l - 1], SS1[p - 1], SS1[q + 1], P);
            for(j = MAX2(i + 4, l - max_interaction_length + 1); j <= q; j++){
              type = ptype[indx[j]+i];
              if (type){
                c3[j-11][tempK][l-j] = MIN2(c3[j-11][tempK][l-j], c3[j-11][max_interaction_length-i+p-1][q-j]+E);
              }
            }/* next j */
          }/* next q */
        }/* next p */
      }/* next l */
    }/* next k */

    /* read out matrix minimum */
    for(j=i+4; j<n1-10; j++) {
      type=ptype[indx[j]+i];
      if (!type) continue;
      int j_pos_end=MIN2(n1-9,j+max_interaction_length);
      for (k=i-1; k>i_pos_begin; k--) {
        for (l=j+1; l<j_pos_end; l++) {
          type2=ptype[indx[l]+k];
          if (!type2) continue;
          E = c3[j-11][max_interaction_length-i+k-1][l-j];
/*           printf("[%d,%d][%d,%d]\t%6.2f\t%6.2f\t%6.2f\n", i, k, l, j, E/100., access_s1[i-k+1][i]/100., access_s1[l-j+1][l]/100.); */
          E+=access_s1[i-k+1][i]+access_s1[l-j+1][l];
          E+=E_ExtLoop(type2,((k>i_pos_begin+1)? SS1[k-1]:-1),((l<j_pos_end-1)? SS1[l+1]:-1),P);
          E+=E_ExtLoop(rtype[type], SS1[j-1], SS1[i+1], P);
          if (E<Emin) {
            Emin=E; k_min=k; l_min=l;
            j_min=j;
          }
        }
      }
    }

    if(Emin  < threshold){
      struc = backtrack_XS(k_min, l_min, i, j_min, max_interaction_length);

      /* lets take care of the dangles */
      /* find best combination */
      int dx_5, dx_3, dy_5, dy_3,dGx,dGy,bonus_x, bonus_y;
      dx_5 = dx_3 = dy_5 = dy_3 = dGx = dGy = bonus_x = bonus_y = 0;
      dGx = access_s1[i-k_min+1][i];
      dGy = access_s1[l_min-j_min+1][l_min];
      PlexHits[NumberOfHits].tb=k_min -10 -dx_5;
      PlexHits[NumberOfHits].te=i -10 + dx_3;
      PlexHits[NumberOfHits].qb=j_min -10 - dy_5;
      PlexHits[NumberOfHits].qe=l_min -10 + dy_3;
      PlexHits[NumberOfHits].ddG=(double) Emin * 0.01;
      PlexHits[NumberOfHits].dG1=(double) dGx*0.01 ;
      PlexHits[NumberOfHits].dG2=(double) dGy*0.01 ;
      PlexHits[NumberOfHits].energy = PlexHits[NumberOfHits].ddG - PlexHits[NumberOfHits].dG1 - PlexHits[NumberOfHits].dG2;
      PlexHits[NumberOfHits].structure = struc;

      /* output: */
      if(PlexHits[NumberOfHits].energy * 100 < threshold){
        if (verbose) printf("%s %3d,%-3d : %3d,%-3d (%5.2f = %5.2f + %5.2f + %5.2f)\n", PlexHits[NumberOfHits].structure, PlexHits[NumberOfHits].tb, PlexHits[NumberOfHits].te, PlexHits[NumberOfHits].qb, PlexHits[NumberOfHits].qe, PlexHits[NumberOfHits].ddG, PlexHits[NumberOfHits].energy, PlexHits[NumberOfHits].dG1, PlexHits[NumberOfHits].dG2);
        NumberOfHits++;
        if(NumberOfHits==PlexHitsArrayLength-1){
          PlexHitsArrayLength*=2;
          PlexHits = (dupVar *) vrna_realloc(PlexHits,sizeof(dupVar)*PlexHitsArrayLength);
        }
      }
    }
  }

  for (i=0; i<(n1-20); i++) {
    for (j=0; j<max_interaction_length; j++) {
      free(c3[i][j]);
    }
    free(c3[i]);
  }
  free(c3);
}


PRIVATE char *backtrack_XS(int k, int l, const int i, const int j, const int max_interaction_length) {
  /* backtrack structure going backwards from i, and forwards from j
     return structure in bracket notation with & as separator */
  int p, q, type, type2, E, traced, i0, j0;
  char *st1, *st2, *struc;
  st1 = (char *) vrna_alloc(sizeof(char)*(i-k+2));
  st1[i-k+1]='\0';
  st2 = (char *) vrna_alloc(sizeof(char)*(l-j+2));
  st2[l-j+1]='\0';

  i0=k; j0=l;
  while (k<=i && l>=j) {
    E = c3[j-11][max_interaction_length-i+k-1][l-j]; traced=0;
    st1[k-i0] = '(';
    st2[l-j] = ')';

    type=ptype[indx[l]+k];
    if (!type) vrna_message_error("backtrack failed in fold duplex bli");
    for (p=k+1; p<=i; p++) {
      for (q=l-1; q>=j; q--) {
        int LE;
        if (p-k+l-q-2>MAXLOOP) break;
        type2=ptype[indx[q]+p];
        if (!type2) continue;
         LE = E_IntLoop(p-k-1, l-q-1, type, rtype[type2], SS1[k+1], SS1[l-1], SS1[p-1], SS1[q+1], P);
         if (E == c3[j-11][max_interaction_length-i+p-1][q-j]+LE) {
          traced=1;
           k=p; l=q;
          break;
        }
      }
      if (traced) break;
    }
    if (!traced) {
      E-=E_ExtLoop(type2, ((k<i)?SS1[k+1]:-1), ((l>j-1)? SS1[l-1]:-1), P);
      break;
      if (E != P->DuplexInit) {
        vrna_message_error("backtrack failed in fold duplex bal");
      } else break;
    }
  }
  struc = (char *) vrna_alloc(k-i0+1+j0-l+1+2);

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

PUBLIC  dupVar  **PKLduplexfold_XS( const char *s1,
                                    int **access_s1,
                                    const int threshold,
                                    const int max_interaction_length,
                                    const int delta){

  if ((!P) || (fabs(P->temperature - temperature)>1e-6))
    update_dfold_params();

  n1 = (int) strlen(s1);
  S1 = encode_sequence(s1, 0);
  SS1 = encode_sequence(s1, 1);

  indx  = vrna_idx_col_wise(n1);
  ptype = (char *) vrna_alloc(sizeof(char)*((n1*(n1+1))/2+2));
  make_ptypes(s1);

  P->DuplexInit=0;
  duplexfold_XS(s1,access_s1,threshold, max_interaction_length);
  free(S1); free(SS1);
  free(indx); free(ptype);
  return NULL;
}

/*---------------------------------UTILS------------------------------------------*/

PRIVATE void update_dfold_params(void){
  vrna_md_t md;
  if(P)
    free(P);
  set_model_details(&md);
  P = vrna_params(&md);
  make_pair_matrix();
}

/*---------------------------------------------------------------------------*/

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
