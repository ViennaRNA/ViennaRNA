/*
  gquad.c

  Ronny Lorenz

  Vienna RNA package
*/

#include <stdio.h>
#include <stdlib.h>
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
  for(i=n-1; i >= 0; i--){
    if(S[i] == 3) gg[i] = gg[i+1]+1;
  }

  int *my_index = get_indx(n);

  /* now find all quadruplexes */
  for(i = 1;
      i < n - (4*VRNA_GQUAD_MIN_STACK_SIZE + 2);
      i++){
    for(L = MIN2(gg[i], VRNA_GQUAD_MAX_STACK_SIZE);
        L >= VRNA_GQUAD_MIN_STACK_SIZE;
        L--){
      for(l1 = VRNA_GQUAD_MIN_LINKER_LENGTH;
          l1 < MIN2(VRNA_GQUAD_MAX_LINKER_LENGTH, n-i - 2*VRNA_GQUAD_MIN_LINKER_LENGTH - 4*L);
          l1++){
        if(gg[i+L+l1] >= L)
          for(l2 = VRNA_GQUAD_MIN_LINKER_LENGTH;
              l2 < MIN2(VRNA_GQUAD_MAX_LINKER_LENGTH, n-i - l1 - VRNA_GQUAD_MIN_LINKER_LENGTH - 4*L);
              l2++){
            if(gg[i + 2*L + l1 + l2] >= L){
              for(l3 = VRNA_GQUAD_MIN_LINKER_LENGTH;
                  l3 < MIN2(VRNA_GQUAD_MAX_LINKER_LENGTH, n-i - l1 - l2 - 4*L);
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

PUBLIC FLT_OR_DBL *get_gquad_pf_matrix(short *S){

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
  for(i=n-1; i >= 0; i--){
    if(S[i] == 3) gg[i] = gg[i+1]+1;
  }

  int *my_index = get_iindx(n);

  /* now find all quadruplexes */
  for(i = 1;
      i < n - (4*VRNA_GQUAD_MIN_STACK_SIZE + 2);
      i++){
    for(L = MIN2(gg[i], VRNA_GQUAD_MAX_STACK_SIZE);
        L >= VRNA_GQUAD_MIN_STACK_SIZE;
        L--){
      for(l1 = VRNA_GQUAD_MIN_LINKER_LENGTH;
          l1 < MIN2(VRNA_GQUAD_MAX_LINKER_LENGTH, n-i - 2*VRNA_GQUAD_MIN_LINKER_LENGTH - 4*L);
          l1++){
        if(gg[i+L+l1] >= L)
          for(l2 = VRNA_GQUAD_MIN_LINKER_LENGTH;
              l2 < MIN2(VRNA_GQUAD_MAX_LINKER_LENGTH, n-i - l1 - VRNA_GQUAD_MIN_LINKER_LENGTH - 4*L);
              l2++){
            if(gg[i + 2*L + l1 + l2] >= L){
              for(l3 = VRNA_GQUAD_MIN_LINKER_LENGTH;
                  l3 < MIN2(VRNA_GQUAD_MAX_LINKER_LENGTH, n-i - l1 - l2 - 4*L);
                  l3++){
                if(gg[i + 3*L + l1 + l2 + l3] >= L){
                  /* insert the quadruplex into the list */
                  j = i+4*L+l1+l2+l3-1;
                  /* do we need scaling here? */
                  data[my_index[i]-j] += exp(-gquad_contribution(L, l1, l2, l3)*10./kT);
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
