/*
  gquad.c

  Ronny Lorenz

  Vienna RNA package
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "../config.h"
#include "fold_vars.h"
#include "data_structures.h"
#include "energy_const.h"
#include "utils.h"
#include "aln_util.h"
#include "gquad.h"

#define MINPSCORE -2 * UNIT

PUBLIC int gquad_contribution(int L, int l1, int l2, int l3){
  int a = -1800;
  int b = 1200;
  int c = INF;

  if(L > VRNA_GQUAD_MAX_STACK_SIZE) return c;
  else if(l1 > VRNA_GQUAD_MAX_LINKER_LENGTH) return c;
  else if(l2 > VRNA_GQUAD_MAX_LINKER_LENGTH) return c;
  else if(l3 > VRNA_GQUAD_MAX_LINKER_LENGTH) return c;
  else return a*(L-1) + b*log(l1 + l2 + l3 -2);
}

PUBLIC int *get_gquad_matrix(short *S){

  int n = S[0];
  int size = (n * (n+1))/2 + 2;

  int *data = (int *)space(sizeof(int) * size);

  int *gg = (int *)space(sizeof(int)*(n+2));

  int init_size = 50;
  int actual_size = 0;
  int i, j, L, l1, l2, l3;

  /* prefill the upper triangular matrix with INF */
  for(i=0;i<size;i++)
    data[i] = INF;

  /* first make the g-island annotation */
  int cnt = 0;
  if(S[n]==3)
    gg[n] = 1;
  for(i=n-1; i > 0; i--){
    if(S[i] == 3) gg[i] = gg[i+1]+1;
  }

  int *my_index = get_indx(n);

  /* now find all quadruplexes */
  for(i = 1;
      i <= n - (4*VRNA_GQUAD_MIN_STACK_SIZE + 2);
      i++){
    for(L = MIN2(gg[i], VRNA_GQUAD_MAX_STACK_SIZE);
        L >= VRNA_GQUAD_MIN_STACK_SIZE;
        L--){
      for(l1 = VRNA_GQUAD_MIN_LINKER_LENGTH;
          l1 <= MIN2(VRNA_GQUAD_MAX_LINKER_LENGTH, 1+n-i - 2*VRNA_GQUAD_MIN_LINKER_LENGTH - 4*L);
          l1++){
        if(gg[i+L+l1] >= L)
          for(l2 = VRNA_GQUAD_MIN_LINKER_LENGTH;
              l2 <= MIN2(VRNA_GQUAD_MAX_LINKER_LENGTH, 1+n-i - l1 - VRNA_GQUAD_MIN_LINKER_LENGTH - 4*L);
              l2++){
            if(gg[i + 2*L + l1 + l2] >= L){
              for(l3 = VRNA_GQUAD_MIN_LINKER_LENGTH;
                  l3 <= MIN2(VRNA_GQUAD_MAX_LINKER_LENGTH, 1+n-i - l1 - l2 - 4*L);
                  l3++){
                if(gg[i + 3*L + l1 + l2 + l3] >= L){
                  /* insert the quadruplex into the list */
                  j = i+4*L+l1+l2+l3-1;
                  data[my_index[j]+i] = MIN2(data[my_index[j]+i], gquad_contribution(L, l1, l2, l3));
                }
              }
            }
          }
      }
    }
  }
  free(my_index);
  free(gg);
  return data;
}

PUBLIC  int gquad_ali_contribution(int i, int L, int l1, int l2, int l3, short **S, int n_seq){
  int en[2];
  gquad_ali_contribution_en(i, L, l1, l2, l3, (const short **)S, n_seq, en);
  return en[0] + en[1];
}

PUBLIC void gquad_ali_contribution_en(int i, int L, int l1, int l2, int l3, const short **S, int n_seq, int en[2]){
  /* check for compatibility in the alignment */
  int s, cnt;
  int penalty = 0;
  int gg_mismatch = 0;
  en[0] = en[1] = 0;

  /* check for compatibility in the alignment */
  for(s=0;s<n_seq;s++){
    int p0, p1, p2, p3, pen;
    p0 = p1 = p2 = p3 = pen = 0;

    /* check bottom layer */
    if(S[s][i] != 3)                      p0 = 1;
    if(S[s][i+L+l1] != 3)                 p1 = 1;
    if(S[s][i + 2*L + l1 + l2] != 3)      p2 = 1;
    if(S[s][i + 3*L + l1 + l2 + l3] != 3) p3 = 1;
    if(p0 || p1 || p2 || p3) pen += VRNA_GQUAD_MISMATCH_PENALTY; /* add 1x penalty for missing bottom layer */

    /* check top layer */
    p0 = p1 = p2 = p3 = 0;
    if(S[s][i+L-1] != 3)                      p0 = 1;
    if(S[s][i+L+l1+L-1] != 3)                 p1 = 1;
    if(S[s][i + 2*L + l1 + l2+L-1] != 3)      p2 = 1;
    if(S[s][i + 3*L + l1 + l2 + l3+L-1] != 3) p3 = 1;
    if(p0 || p1 || p2 || p3) pen += VRNA_GQUAD_MISMATCH_PENALTY; /* add 1x penalty for missing top layer */

    /* check inner layers */
    for(cnt=1;cnt<L-1;cnt++){
      if(S[s][i+cnt] != 3)                      p0 = 1;
      if(S[s][i+L+l1+cnt] != 3)                 p1 = 1;
      if(S[s][i + 2*L + l1 + l2+cnt] != 3)      p2 = 1;
      if(S[s][i + 3*L + l1 + l2 + l3+cnt] != 3) p3 = 1;
      if(p0 || p1 || p2 || p3) pen += 2*VRNA_GQUAD_MISMATCH_PENALTY; /* add 2x penalty for missing inner layer */
    }

    /* if all layers are missing, we have a complete gg mismatch */
    if(pen >= (2*VRNA_GQUAD_MISMATCH_PENALTY * (L-1)))
      gg_mismatch++;
    /* add the penalty to the score */
    penalty += pen;
  }
  /* only one ggg mismatch allowed */
  if(gg_mismatch > VRNA_GQUAD_MISMATCH_NUM_ALI){
    en[0] = n_seq * gquad_contribution(L, l1, l2, l3);
    en[1] = INF;
  } else {
    en[0] = n_seq * gquad_contribution(L, l1, l2, l3);
    en[1] = penalty;
  }
}

PUBLIC int *get_gquad_ali_matrix( short *S_cons,
                                  short **S,
                                  int n_seq,
                                  int *pscore){

  int n, size, *data, *gg, **ggg, init_size, actual_size, s, i, j, L, l1, l2, l3, cnt, *my_index, *gggg_score;


  n           = S[0][0];
  size        = (n * (n+1))/2 + 2;
  data        = (int *)space(sizeof(int) * size);
  gg          = (int *)space(sizeof(int)*(n+1));
  ggg         = (int **)space(sizeof(int*)*(n_seq));
  my_index    = get_indx(n);
  init_size   = 50;
  actual_size = 0;
  cnt         = 0;
  for(s=0;s<n_seq;s++)
    ggg[s] = (int *)space(sizeof(int) * (n+1));


  /*  make the g-island annotation for consensus sequence */
  if(S_cons[n]==3) gg[n] = 1;
  for(i = n-1; i > 0; i--)
    if(S_cons[i] == 3) gg[i] = gg[i+1]+1;

  /* make the g-island annotation for each sequence in the alignment */
  for(s=0;s<n_seq;s++){
    if(S[s][n] == 3) ggg[s][n] = 1;
    for(i=n-1;i>0;i--)
      if(S[s][i] == 3) ggg[s][i] = ggg[s][i+1]+1;
  }

  /* prefill the upper triangular matrix with INF */
  for(i=0;i<size;i++) data[i] = INF;

  /* now find all quadruplexes compatible with consensus sequence */
  for(i = 1;
      i <= n - (4*VRNA_GQUAD_MIN_STACK_SIZE + 2);
      i++){
    for(L = MIN2(gg[i], VRNA_GQUAD_MAX_STACK_SIZE);
        L >= VRNA_GQUAD_MIN_STACK_SIZE;
        L--){
      for(l1 = VRNA_GQUAD_MIN_LINKER_LENGTH;
          l1 <= MIN2(VRNA_GQUAD_MAX_LINKER_LENGTH, 1+n-i - 2*VRNA_GQUAD_MIN_LINKER_LENGTH - 4*L);
          l1++){
        if(gg[i+L+l1] >= L)
          for(l2 = VRNA_GQUAD_MIN_LINKER_LENGTH;
              l2 <= MIN2(VRNA_GQUAD_MAX_LINKER_LENGTH, 1+n-i - l1 - VRNA_GQUAD_MIN_LINKER_LENGTH - 4*L);
              l2++){
            if(gg[i + 2*L + l1 + l2] >= L){
              for(l3 = VRNA_GQUAD_MIN_LINKER_LENGTH;
                  l3 <= MIN2(VRNA_GQUAD_MAX_LINKER_LENGTH, 1+n-i - l1 - l2 - 4*L);
                  l3++){
                if(gg[i + 3*L + l1 + l2 + l3] >= L){
                  j = i+4*L+l1+l2+l3-1;
                  /* check for compatibility in the alignment */
                  int en = gquad_ali_contribution(i, L, l1, l2, l3, S, n_seq);
                  data[my_index[j]+i] = MIN2(data[my_index[j]+i], en);
                }
              }
            }
          }
      }
    }
  }

  /* clean up */
  free(my_index);
  free(gg);
  for(s=0;s<n_seq;s++)
    free(ggg[s]);
  free(ggg);
  return data;
}

PUBLIC FLT_OR_DBL *get_gquad_pf_matrix(short *S, FLT_OR_DBL *scale){

  int n = S[0];
  int size = (n * (n+1))/2 + 2;

  FLT_OR_DBL *data = (FLT_OR_DBL *)space(sizeof(FLT_OR_DBL) * size);

  double kT = (temperature+K0)*GASCONST; /* in kcal/mol */

  int *gg = (int *)space(sizeof(int)*(n+2));

  int init_size = 50;
  int actual_size = 0;
  int i, j, L, l1, l2, l3;

  /* no need for prefill since everything is initialized with 0 */

  /* first make the g-island annotation */
  int cnt = 0;
  if(S[n]==3)
    gg[n] = 1;
  for(i=n-1; i > 0; i--){
    if(S[i] == 3) gg[i] = gg[i+1]+1;
  }

  int *my_index = get_iindx(n);

  /* now find all quadruplexes */
  for(i = 1;
      i <= n - (4*VRNA_GQUAD_MIN_STACK_SIZE + 2);
      i++){
    for(L = MIN2(gg[i], VRNA_GQUAD_MAX_STACK_SIZE);
        L >= VRNA_GQUAD_MIN_STACK_SIZE;
        L--){
      for(l1 = VRNA_GQUAD_MIN_LINKER_LENGTH;
          l1 <= MIN2(VRNA_GQUAD_MAX_LINKER_LENGTH, 1+n-i - 2*VRNA_GQUAD_MIN_LINKER_LENGTH - 4*L);
          l1++){
        if(gg[i+L+l1] >= L)
          for(l2 = VRNA_GQUAD_MIN_LINKER_LENGTH;
              l2 <= MIN2(VRNA_GQUAD_MAX_LINKER_LENGTH, 1+n-i - l1 - VRNA_GQUAD_MIN_LINKER_LENGTH - 4*L);
              l2++){
            if(gg[i + 2*L + l1 + l2] >= L){
              for(l3 = VRNA_GQUAD_MIN_LINKER_LENGTH;
                  l3 <= MIN2(VRNA_GQUAD_MAX_LINKER_LENGTH, 1+n-i - l1 - l2 - 4*L);
                  l3++){
                if(gg[i + 3*L + l1 + l2 + l3] >= L){
                  /* insert the quadruplex into the list */
                  j = i+4*L+l1+l2+l3-1;
 
                  /* do we need scaling here? */
                  data[my_index[i]-j] += exp(-gquad_contribution(L, l1, l2, l3)*10./kT)*scale[j-i+1];
                }
              }
            }
          }
      }
    }
  }
  free(my_index);
  free(gg);
  return data;
}



PUBLIC plist *Gquadcomputeinnerprobability(short *S, FLT_OR_DBL *G, FLT_OR_DBL  *probs, FLT_OR_DBL *scale){ 

  int n = S[0];
  int size = (n * (n+1))/2 + 2;

  FLT_OR_DBL *data = (FLT_OR_DBL *)space(sizeof(FLT_OR_DBL) * size);
  FLT_OR_DBL *tempprobs = (FLT_OR_DBL *)space(sizeof(FLT_OR_DBL) * size);
  plist *pl=(plist *)space((S[0])*sizeof(plist));

  double kT = (temperature+K0)*GASCONST; /* in kcal/mol */

  int *gg = (int *)space(sizeof(int)*(n+2));
  int counter=0;
  int init_size = 50;
  int actual_size = 0;
  int i, j, L, l1, l2, l3;
  FLT_OR_DBL  e_con;
  /* no need for prefill since everything is initialized with 0 */

  /* first make the g-island annotation */
  int cnt = 0;
  if(S[n]==3)
    gg[n] = 1;
  for(i=n-1; i > 0; i--){
    if(S[i] == 3) gg[i] = gg[i+1]+1;
  }

  int *my_index = get_iindx(n);

  /* now find all quadruplexes */
  for(i = 1;
      i <= n - (4*VRNA_GQUAD_MIN_STACK_SIZE + 2);
      i++){
    for(L = MIN2(gg[i], VRNA_GQUAD_MAX_STACK_SIZE);
        L >= VRNA_GQUAD_MIN_STACK_SIZE;
        L--){
      for(l1 = VRNA_GQUAD_MIN_LINKER_LENGTH;
          l1 <= MIN2(VRNA_GQUAD_MAX_LINKER_LENGTH, 1+n-i - 2*VRNA_GQUAD_MIN_LINKER_LENGTH - 4*L);
          l1++){
        if(gg[i+L+l1] >= L)
          for(l2 = VRNA_GQUAD_MIN_LINKER_LENGTH;
              l2 <= MIN2(VRNA_GQUAD_MAX_LINKER_LENGTH, 1+n-i - l1 - VRNA_GQUAD_MIN_LINKER_LENGTH - 4*L);
              l2++){
            if(gg[i + 2*L + l1 + l2] >= L){
              for(l3 = VRNA_GQUAD_MIN_LINKER_LENGTH;
                  l3 <= MIN2(VRNA_GQUAD_MAX_LINKER_LENGTH, 1+n-i - l1 - l2 - 4*L);
                  l3++){
                if(gg[i + 3*L + l1 + l2 + l3] >= L){
                  /* insert the quadruplex into the list */
                    int x;
                    j = i+4*L+l1+l2+l3-1;
                    /* do we need scaling here? */
                  e_con=probs[my_index[i]-j]*exp(-gquad_contribution(L, l1, l2, l3)*10./kT)*scale[j-i+1]/G[my_index[i]-j];/*a): indexes? b: ist gesamtQ schon dividiert da?*/
                  
                  for (x=0; x<L; x++) {
                    tempprobs[my_index[i+x]-(i+x+3*L+l1+l2+l3)]+= e_con;
                    tempprobs[my_index[i+x]-(i+x+L+l1)] += e_con; /*prob??*/
                    tempprobs[my_index[i+x+L+l1]-(i+x+2*L+l1+l2)]+= e_con;
                    tempprobs[my_index[i+x+2*L+l1+l2]-(i+x+3*L+l1+l2+l3)]+= e_con;
                  }
                }
              }
            }
          }
      }
    }
  }
  for (i=0;i<n; i++) {
    for (j=i; j<=n; j++) {
      if (tempprobs[my_index[i]-j]>0.) {
        pl[counter].i=i;
        pl[counter].j=j;
        pl[counter++].p=tempprobs[my_index[i]-j];
      }
    }
    
  }
  pl[counter].i=pl[counter].j=0;
  pl[counter++].p=0.;
  /* shrink memory to actual size needed */
  pl = (plist *) xrealloc(pl, counter * sizeof(plist));

  free(my_index);
  free (tempprobs);
 return pl;
}

PUBLIC plist *get_plist_gquad_from_pr(short *S,
                                      int gi,
                                      int gj,
                                      FLT_OR_DBL *G,
                                      FLT_OR_DBL *probs,
                                      FLT_OR_DBL *scale){ 
  int L, l[3];
  return  get_plist_gquad_from_pr_max(S, gi, gj, G, probs, scale, &L, l);
}


PUBLIC plist *get_plist_gquad_from_pr_max(short *S,
                                      int gi,
                                      int gj,
                                      FLT_OR_DBL *G,
                                      FLT_OR_DBL *probs,
                                      FLT_OR_DBL *scale,
                                      int *Lmax,
                                      int lmax[3]){ 
  int n = S[0];
  int size = (n * (n+1))/2 + 2;

  FLT_OR_DBL *data = (FLT_OR_DBL *)space(sizeof(FLT_OR_DBL) * size);
  FLT_OR_DBL *tempprobs = (FLT_OR_DBL *)space(sizeof(FLT_OR_DBL) * size);
  plist *pl=(plist *)space((S[0]*S[0])*sizeof(plist));

  double kT = (temperature+K0)*GASCONST; /* in kcal/mol */

  int *gg = (int *)space(sizeof(int)*(n+2));
  int counter=0;
  int init_size = 50;
  int actual_size = 0;
  int i, j, L, l1, l2, l3;
  int cnt = 0;
  FLT_OR_DBL  e_con;
  int *my_index = get_iindx(n);
  FLT_OR_DBL ggg_max_sum = 0.;

  /* first make the g-island annotation */
  if(S[gj]==3) gg[gj] = 1;
  for(i=gj-1; i >= gi; i--)
    if(S[i] == 3) gg[i] = gg[i+1]+1;

  i = gi;
  /* now find all quadruplexes */
  for(L = MIN2(gg[i], VRNA_GQUAD_MAX_STACK_SIZE);
        L >= VRNA_GQUAD_MIN_STACK_SIZE;
        L--){
      for(l1 = VRNA_GQUAD_MIN_LINKER_LENGTH;
          l1 <= MIN2(VRNA_GQUAD_MAX_LINKER_LENGTH, 1+gj-i - 2*VRNA_GQUAD_MIN_LINKER_LENGTH - 4*L);
          l1++){
        if(gg[i+L+l1] >= L)
          for(l2 = VRNA_GQUAD_MIN_LINKER_LENGTH;
              l2 <= MIN2(VRNA_GQUAD_MAX_LINKER_LENGTH, 1+gj-i - l1 - VRNA_GQUAD_MIN_LINKER_LENGTH - 4*L);
              l2++){
            if(gg[i + 2*L + l1 + l2] >= L){
              for(l3 = VRNA_GQUAD_MIN_LINKER_LENGTH;
                  l3 <= MIN2(VRNA_GQUAD_MAX_LINKER_LENGTH, 1+gj-i - l1 - l2 - 4*L);
                  l3++){
                if(gg[i + 3*L + l1 + l2 + l3] >= L){
                  /* insert the quadruplex into the list */
                  int x;
                  FLT_OR_DBL gggm = 0;
                  j = i+4*L+l1+l2+l3-1;
                  if(j!=gj) continue;
                  /* do we need scaling here? */
                  /* is it (p | ggg(gi,gj)) * exp(ggg(i,j)) / G(i,j) ? */
                  e_con = probs[my_index[i]-j] * exp(-gquad_contribution(L, l1, l2, l3)*10./kT) * scale[j-i+1] / G[my_index[i]-j];
                  for (x=0; x<L; x++) {
                    gggm += (e_con + (1-e_con));
                    tempprobs[my_index[i+x]-(i+x+3*L+l1+l2+l3)]+= e_con;
                    tempprobs[my_index[i+x]-(i+x+L+l1)] += e_con; /*prob??*/
                    tempprobs[my_index[i+x+L+l1]-(i+x+2*L+l1+l2)]+= e_con;
                    tempprobs[my_index[i+x+2*L+l1+l2]-(i+x+3*L+l1+l2+l3)]+= e_con;
                  }
                  if(gggm > ggg_max_sum){
                    ggg_max_sum = gggm;
                    *Lmax = L; lmax[0] = l1; lmax[1] = l2; lmax[2] = l3;
                  }
                }
              }
            }
          }
      }
  }
  for (i=gi;i<gj; i++) {
    for (j=i; j<=gj; j++) {
      if (tempprobs[my_index[i]-j]>0.) {
        pl[counter].i=i;
        pl[counter].j=j;
        pl[counter++].p=tempprobs[my_index[i]-j];
      }
    }
  }
  pl[counter].i = pl[counter].j = 0;
  pl[counter++].p=0.;
  /* shrink memory to actual size needed */
  pl = (plist *) xrealloc(pl, counter * sizeof(plist));

  free(my_index);
  free (tempprobs);
  return pl;
}

PUBLIC plist *get_plist_gquad_from_db(const char *structure, float pr){
  int i, x, size, actual_size, L, cL, n;
  plist *pl;
  int   gg1, gg2, gg3, gg4;

  actual_size = 0;
  n           = 2;
  size        = strlen(structure);
  pl          = (plist *)space(n*size*sizeof(plist));

  /* seek to first occurence of '+' or end of string */
  for(i=0; i < size && structure[i] != '+'; i++);
  while(i < size){
    gg1 = gg2 = gg3 = gg4 = -1;
    /* get size of gquad */
    for(gg1 = i, L = 0; i < size && structure[i] == '+'; i++, L++);
    /* now find the rest of the quad */
    for(; i < size && structure[i] != '+'; i++);
    if(i >= size)
      nrerror("get_plist_gquad_from_db@gquad.c: misformatted dot bracket string");
    for(gg2 = i, cL = 0; i < size && structure[i] == '+'; i++, cL++);
    if(L != cL)
      nrerror("get_plist_gquad_from_db@gquad.c: misformatted dot bracket string");
    for(; i < size && structure[i] != '+'; i++);
    if(i >= size)
      nrerror("get_plist_gquad_from_db@gquad.c: misformatted dot bracket string");
    for(gg3 = i, cL = 0; i < size && structure[i] == '+'; i++, cL++);
    if(L != cL)
      nrerror("get_plist_gquad_from_db@gquad.c: misformatted dot bracket string");
    for(; i < size && structure[i] != '+'; i++);
    if(i >= size)
      nrerror("get_plist_gquad_from_db@gquad.c: misformatted dot bracket string");
    for(gg4 = i, cL = 0; i < size && structure[i] == '+'; i++, cL++);
    if(L != cL)
      nrerror("get_plist_gquad_from_db@gquad.c: misformatted dot bracket string");
    /* finally add the gquadruplex to the list */
    for (x=0; x<L; x++) {
      if (actual_size >= n * size - 5){
        n *= 2;
        pl = (plist *)xrealloc(pl, n * size * sizeof(plist));
      }
      pl[actual_size].i = gg1 + x + 1;
      pl[actual_size].j = gg4 + x + 1;
      pl[actual_size].p = pr;
      pl[actual_size++].type = 0;

      pl[actual_size].i = gg1 + x + 1;
      pl[actual_size].j = gg2 + x + 1;
      pl[actual_size].p = pr;
      pl[actual_size++].type = 0;

      pl[actual_size].i = gg2 + x + 1;
      pl[actual_size].j = gg3 + x + 1;
      pl[actual_size].p = pr;
      pl[actual_size++].type = 0;

      pl[actual_size].i = gg3 + x + 1;
      pl[actual_size].j = gg4 + x + 1;
      pl[actual_size].p = pr;
      pl[actual_size++].type = 0;
    }
    for(; i < size && structure[i] != '+'; i++);
  }

  pl[actual_size].i = pl[actual_size].j = 0;
  pl[actual_size++].p = 0;
  pl = (plist *)xrealloc(pl, actual_size * sizeof(plist));
  return pl;
}

PUBLIC  int *make_ggggscores( const short *const* S,
                              const char **AS,
                              int n_seq,
                              const char *structure){

#if 0
  /* calculate co-variance bonus for each pair depending on  */
  /* compensatory/consistent mutations and incompatible seqs */
  /* should be 0 for conserved pairs, >0 for good pairs      */
#define NONE -10000 /* score for forbidden pairs */
  int n,i,j,k,l,s, *indx;
  int *gggg_score;
  int gggg_dm[5][5] = {{2, 2, 2, 2, 2}, /* hamming distance of gg to any other nucleotide arrangements */
                      {2, 2, 2, 1, 2},  /* A */
                      {2, 2, 2, 1, 2}, /* C */
                      {1, 1, 1, 0, 1}, /* G */
                      {2, 2, 2, 1, 2}, /* U */};

  n = S[0][0];  /* length of seqs */
  indx = get_indx((unsigned int)n);

  for (i=1; i<n; i++){
    for (j=i+1; (j<i+TURN+1) && (j<=n); j++)
      gggg_score[indx[j]+i] = NONE;

    for (j=i+TURN+1; j<=n; j++) {
      int pfreq[8]={0,0,0,0,0,0,0,0};
      double score;
      for (s=0; s<n_seq; s++) {
        int type;
        if (S[s][i]==0 && S[s][j]==0) type = 7; /* gap-gap  */
        else {
          if ((AS[s][i] == '~')||(AS[s][j] == '~')) type = 7;
          else type = pair[S[s][i]][S[s][j]];
        }
        pfreq[type]++;
      }
      if (pfreq[0]*2+pfreq[7]>n_seq) { pscore[indx[j]+i] = NONE; continue;}
      for (k=1,score=0; k<=6; k++) /* ignore pairtype 7 (gap-gap) */
        for (l=k; l<=6; l++)
          score += pfreq[k]*pfreq[l]*dm[k][l];
      /* counter examples score -1, gap-gap scores -0.25   */
      pscore[indx[j]+i] = cv_fact *
        ((UNIT*score)/n_seq - nc_fact*UNIT*(pfreq[0] + pfreq[7]*0.25));
    }
  }

  if (noLonelyPairs) /* remove unwanted pairs */
    for (k=1; k<n-TURN-1; k++)
      for (l=1; l<=2; l++) {
        int type,ntype=0,otype=0;
        i=k; j = i+TURN+l;
        type = pscore[indx[j]+i];
        while ((i>=1)&&(j<=n)) {
          if ((i>1)&&(j<n)) ntype = pscore[indx[j+1]+i-1];
          if ((otype<cv_fact*MINPSCORE)&&(ntype<cv_fact*MINPSCORE))  /* too many counterexamples */
            pscore[indx[j]+i] = NONE; /* i.j can only form isolated pairs */
          otype =  type;
          type  = ntype;
          i--; j++;
        }
      }


  if (fold_constrained&&(structure!=NULL)) {
    int psij, hx, hx2, *stack, *stack2;
    stack = (int *) space(sizeof(int)*(n+1));
    stack2 = (int *) space(sizeof(int)*(n+1));

    for(hx=hx2=0, j=1; j<=n; j++) {
      switch (structure[j-1]) {
      case 'x': /* can't pair */
        for (l=1; l<j-TURN; l++) pscore[indx[j]+l] = NONE;
        for (l=j+TURN+1; l<=n; l++) pscore[indx[l]+j] = NONE;
        break;
      case '(':
        stack[hx++]=j;
        /* fallthrough */
      case '[':
        stack2[hx2++]=j;
        /* fallthrough */
      case '<': /* pairs upstream */
        for (l=1; l<j-TURN; l++) pscore[indx[j]+l] = NONE;
        break;
      case ']':
        if (hx2<=0) {
          fprintf(stderr, "%s\n", structure);
          nrerror("unbalanced brackets in constraints");
        }
        i = stack2[--hx2];
        pscore[indx[j]+i]=NONE;
        break;
      case ')':
        if (hx<=0) {
          fprintf(stderr, "%s\n", structure);
          nrerror("unbalanced brackets in constraints");
        }
        i = stack[--hx];
        psij = pscore[indx[j]+i]; /* store for later */
        for (k=j; k<=n; k++)
          for (l=i; l<=j; l++)
            pscore[indx[k]+l] = NONE;
        for (l=i; l<=j; l++)
          for (k=1; k<=i; k++)
            pscore[indx[l]+k] = NONE;
        for (k=i+1; k<j; k++)
          pscore[indx[k]+i] = pscore[indx[j]+k] = NONE;
        pscore[indx[j]+i] = (psij>0) ? psij : 0;
        /* fallthrough */
      case '>': /* pairs downstream */
        for (l=j+TURN+1; l<=n; l++) pscore[indx[l]+j] = NONE;
        break;
      }
    }
    if (hx!=0) {
      fprintf(stderr, "%s\n", structure);
      nrerror("unbalanced brackets in constraint string");
    }
    free(stack); free(stack2);
  }
  /*free dm */
  for (i=0; i<7;i++) {
    free(dm[i]);
  }
  free(dm);
#endif
}

int parse_gquad(const char *struc, int *L, int l[3]) {
  /* given a dot-bracket structure (possibly) containing gquads encoded
     by '+' signs, find first gquad, return end position or 0 if none found
     Upon return L and l[] contain the number of stacked layers, as well as
     the lengths of the linker regions.  
     To parse a string with many gquads, call parse_gquad repeatedly e.g.
     end1 = parse_gquad(struc, &L, l); ... ;
     end2 = parse_gquad(struc+end1, &L, l); end2+=end1; ... ;
     end3 = parse_gquad(struc+end2, &L, l); end3+=end2; ... ; 
  */
  int i, il, start, end, len;
  
  for (i=0; struc[i] && struc[i]!='+'; i++);
  if (struc[i] == '+') { /* start of gquad */
    for (il=0; il<=3; il++) {
      start=i; /* pos of first '+' */
      while (struc[++i] == '+');
      end=i; len=end-start; 
      if (il==0) *L=len;
      else if (len!=*L)
        nrerror("unequal stack lengths in gquad");
      if (il==3) break;
      while (struc[++i] == '.'); /* linker */
      l[il] = i-end;
      if (struc[i] != '+') nrerror("illegal character in gquad linker region");
    }
  }
  else return 0;
  /* printf("gquad at %d %d %d %d %d\n", end, *L, l[0], l[1], l[2]); */
  return end;
}

  
