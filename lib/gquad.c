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
#include "utils.h"
#include "gquad.h"


PUBLIC  int **annotate_gquadruplexes(short  *S){

  int n = S[0];
  int **data = (int **)space(sizeof(int *) * 5);

  int *gg = (int *)space(sizeof(int)*(n+2));

  int init_size = 50;
  int actual_size = 0;
  int i, L, l1, l2, l3;

  for(i=0;i<5;i++)
    data[i] = (int *)space(sizeof(int) * init_size);

  /* first make the g-island annotation */
  int cnt = 0;
  for(i=n-1; i >= 0; i--){
    if(S[i] == 3) gg[i] = gg[i+1]+1;
  }

  /* now find all quadruplexes */
  for(i = 1;
      i < n - (4*VRNA_GQUAD_MIN_STACK_SIZE + 2);
      i++){
    for(L = MIN2(gg[i], VRNA_GQUAD_MAX_STACK_SIZE);
        L >= VRNA_GQUAD_MIN_STACK_SIZE;
        L--){
      for(l1 = VRNA_GQUAD_MIN_LINKER_LENGTH;
          l1 < (VRNA_GQUAD_MAX_LINKER_LENGTH - 2*VRNA_GQUAD_MIN_LINKER_LENGTH);
          l1++){
        if(gg[i+L+l1] >= L)
          for(l2 = VRNA_GQUAD_MIN_LINKER_LENGTH;
              l2 < (VRNA_GQUAD_MAX_LINKER_LENGTH - l1 - VRNA_GQUAD_MIN_LINKER_LENGTH);
              l2++){
            if(gg[i + 2*L + l1 + l2] >= L){
              for(l3 = VRNA_GQUAD_MIN_LINKER_LENGTH;
                  l3 < (VRNA_GQUAD_MAX_LINKER_LENGTH - l1 - l2);
                  l3++){
                if(gg[i + 3*L + l1 + l2 + l3] >= L){
                  /* insert the quadruplex into the list */
                  printf("\n\tgg: %d, %d, %d, %d, %d\n", i, L, l1, l2, l3);
                }
              }
            }
          }
      }
    }
  }

  free(gg);
  return data;
}

PRIVATE void print_gquad(int i, int n, int L, int l1, int l2, int l3){
  int cnt;
  for(cnt=0;cnt<i-1;cnt++){
    printf(".");
  }
  for(cnt=0;cnt<L;cnt++)
    printf("#");
  for(cnt=0;cnt<l1;cnt++)
    printf(".");
  for(cnt=0;cnt<L;cnt++)
    printf("#");
  for(cnt=0;cnt<l2;cnt++)
    printf(".");
  for(cnt=0;cnt<L;cnt++)
    printf("#");
  for(cnt=0;cnt<l3;cnt++)
    printf(".");
  for(cnt=0;cnt<L;cnt++)
    printf("#");
  for(cnt=0;cnt<n-i+1-4*L-l1-l2-l3; cnt++)
    printf(".");
  printf("\tgg: %d, %d, %d, %d, %d, %d\n", i, L, l1, l2, l3, i+4*L+l1+l2+l3-1);
}

PRIVATE int gquad_contribution(int L, int l1, int l2, int l3){
  return -(4*L + l1 + l2 + l3);
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
                  //print_gquad(i, n, L, l1, l2, l3);
                  data[my_index[j]+i] = MIN2(data[my_index[j]+1], gquad_contribution(L, l1, l2, l3));
                }
              }
            }
          }
      }
    }
  }

  for(i=1;i<n;i++){
    for(j=i+(4*VRNA_GQUAD_MIN_STACK_SIZE-1)+(3*VRNA_GQUAD_MIN_LINKER_LENGTH); j<=n;j++){
      /* lets extend the precomputed quad contributions by appending unpaired nucleotides */
      data[my_index[j]+i] = MIN2(data[my_index[j]+i], data[my_index[j-1]+i]);
      /* shouldn't we also include constructions like i - ggg1 - ggg2 - j ??? */
    }
  }
  free(my_index);
  free(gg);
  return data;
}
