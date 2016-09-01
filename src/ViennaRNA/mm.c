/*
      Implementation of Nussinov Maximum Matching
      Ronny Lorenz
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
#include "ViennaRNA/pair_mat.h"
#include "ViennaRNA/params.h"

/* the encoded string MUST have the length of the sequence at position 0!!! */
PUBLIC unsigned int maximumMatching(const char *string){
  unsigned int i, j, l, length, max = 0;
  unsigned int *mm;            /* holds maximum matching on subsequence [i,j] */
  short *encodedString = encode_sequence(string, 0);
  int *iindx = vrna_idx_row_wise((unsigned) encodedString[0]);
  make_pair_matrix();
  length = (unsigned int)encodedString[0];
  mm = (unsigned int *) vrna_alloc(sizeof(unsigned int)*((length*(length+1))/2+2));
  for(j = 1; j<=length; j++)
    for(i=(j>TURN?(j-TURN):1); i<j; i++)
      mm[iindx[i]-j] = 0;
  for(i=length-TURN-1;i>0; i--)
    for(j=i+TURN+1; j<= length; j++){
      max = mm[iindx[i]-j+1];
      for(l=j-TURN-1; l>=i; l--)
        if(pair[encodedString[l]][encodedString[j]])
          max = MAX2(max, ((l>i) ? mm[iindx[i]-l+1] : 0) + 1 + mm[iindx[l+1]-j+1]);
       mm[iindx[i]-j] = max;
    }
  max = mm[iindx[1]-length];
  free(mm);
  free(iindx);
  free(encodedString);
  return max;
}

/* the encoded string MUST have the length of the sequence at position 0!!! */
PUBLIC unsigned int *maximumMatchingConstraint(const char *string, short *ptable){
  unsigned int i, j, l, length, max = 0;
  unsigned int *mm;            /* holds maximum matching on subsequence [i,j] */
  short *encodedString = encode_sequence(string, 0);
  int *iindx = vrna_idx_row_wise((unsigned) encodedString[0]);
  make_pair_matrix();
  length = (unsigned int)encodedString[0];
  mm = (unsigned int *) vrna_alloc(sizeof(unsigned int)*((length*(length+1))/2+2));
  for(j = 1; j<=length; j++)
    for(i=(j>TURN?(j-TURN):1); i<j; i++)
      mm[iindx[i]-j] = 0;
  for(i=length-TURN-1;i>0; i--)
    for(j=i+TURN+1; j<= length; j++){
      max = mm[iindx[i]-j+1];
      for(l=j-TURN-1; l>=i; l--){
        if(pair[encodedString[l]][encodedString[j]]){
          if(ptable[l] != j)
            max = MAX2(max, ((l>i) ? mm[iindx[i]-l+1] : 0) + 1 + mm[iindx[l+1]-j+1]);
        }
      }
      mm[iindx[i]-j] = max;
    }
  free(iindx);
  free(encodedString);
  return mm;
}

/* the encoded string MUST have the length of the sequence at position 0!!! */
PUBLIC unsigned int *maximumMatching2Constraint(const char *string, short *ptable, short *ptable2){
  unsigned int i, j, l, length, max = 0;
  unsigned int *mm;            /* holds maximum matching on subsequence [i,j] */
  short *encodedString = encode_sequence(string, 0);
  int *iindx = vrna_idx_row_wise((unsigned) encodedString[0]);
  make_pair_matrix();
  length = (unsigned int)encodedString[0];
  mm = (unsigned int *) vrna_alloc(sizeof(unsigned int)*((length*(length+1))/2+2));
  for(j = 1; j<=length; j++)
    for(i=(j>TURN?(j-TURN):1); i<j; i++)
      mm[iindx[i]-j] = 0;
  for(i=length-TURN-1;i>0; i--)
    for(j=i+TURN+1; j<= length; j++){
      max = mm[iindx[i]-j+1];
      for(l=j-TURN-1; l>=i; l--){
        if(pair[encodedString[l]][encodedString[j]]){
          if(ptable[l] != j && ptable2[l] != j)
            max = MAX2(max, ((l>i) ? mm[iindx[i]-l+1] : 0) + 1 + mm[iindx[l+1]-j+1]);
        }
      }
      mm[iindx[i]-j] = max;
    }
  free(iindx);
  free(encodedString);
  return mm;
}

