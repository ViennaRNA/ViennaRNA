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
