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
#include "gquad.h"


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

  /* first make the g-island annotation */
  if(S[gj]==3) gg[gj] = 1;
  for(i=gj-1; i >= gi; i--)
    if(S[i] == 3) gg[i] = gg[i+1]+1;

  /* now find all quadruplexes */
  for(i = gi;
      i <= gj - (4*VRNA_GQUAD_MIN_STACK_SIZE + 2);
      i++){
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
                  j = i+4*L+l1+l2+l3-1;
                  /* do we need scaling here? */
                  /* is it (p | ggg(gi,gj)) * exp(ggg(i,j)) / G(i,j) ? */
                  e_con = probs[my_index[gi]-gj] * exp(-gquad_contribution(L, l1, l2, l3)*10./kT) * scale[j-i+1] / G[my_index[i]-j];

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
