/*
 *    Implementation of Nussinov Maximum Matching
 *    Ronny Lorenz
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <string.h>
#include "ViennaRNA/utils/basic.h"
#include "ViennaRNA/params/default.h"
#include "ViennaRNA/fold_vars.h"
#include "ViennaRNA/pair_mat.h"
#include "ViennaRNA/params/basic.h"

#include "ViennaRNA/mm.h"

PUBLIC int
vrna_maximum_matching(vrna_fold_compound_t *fc)
{
  unsigned char *mx, *hc_up;
  int           i, j, l, n, turn, *mm, max, max2, max3;
  vrna_hc_t     *hc;

  n     = (int)fc->length;
  turn  = fc->params->model_details.min_loop_size;
  hc    = fc->hc;
  mx    = hc->mx;
  hc_up = (unsigned char *)vrna_alloc(sizeof(unsigned char) * n);
  mm    = (int *)vrna_alloc(sizeof(int) * (n * n));

  /* comply with hard constraints for unpaired positions */
  for (i = n - 1; i >= 0; i--)
    if (mx[n * (i + 1) + i + 1] & VRNA_CONSTRAINT_CONTEXT_ALL_LOOPS)
      hc_up[i] = 1;

  /* initialize DP matrix */
  for (j = 0; j < n; j++)
    for (i = (j >= turn ? (j - turn) : 0); i < j; i++) {
      mm[n * i + j] = hc_up[i] ? ((i > 0) ? mm[n * j + i - 1] : 0) : -1;
      mm[n * j + i] = mm[n * i + j];
    }

  /* start recursions */
  for (i = n - turn - 2; i >= 0; i--)
    for (j = i + turn + 1; j < n; j++) {
      max = -1;

      /* 1st case: i pairs with j */
      if (mx[n * (i + 1) + j + 1] & VRNA_CONSTRAINT_CONTEXT_ALL_LOOPS) {
        max2 = mm[n * (i + 1) + j - 1];

        if (max2 != -1) {
          max2 += 1;
          if (max < max2)
            max = max2;
        }
      }

      /* 2nd case: i is unpaired */
      if (hc_up[i]) {
        max2 = mm[n * (i + 1) + j];
        if (max < max2)
          max = max2;
      }

      /* 3rd case: j is unpaired */
      if (hc_up[j]) {
        max2 = mm[n * i + j - 1];
        if (max < max2)
          max = max2;
      }

      /* 4th case: split at l */
      for (l = i + 1; l < j; l++) {
        max2  = mm[n * i + l - 1];
        max3  = mm[n * j + l];
        if ((max2 != -1) && (max3 != -1)) {
          max2 += max3;
          if (max < max2)
            max = max2;
        }
      }

      mm[n * i + j] = max;
      mm[n * j + i] = max;
    }

  max = mm[n - 1];

  free(mm);
  free(hc_up);

  return max;
}


PUBLIC int
vrna_maximum_matching_simple(const char *sequence)
{
  int                   max_pairs;
  vrna_fold_compound_t  *fc;

  fc        = vrna_fold_compound(sequence, NULL, VRNA_OPTION_DEFAULT);
  max_pairs = vrna_maximum_matching(fc);

  vrna_fold_compound_free(fc);

  return max_pairs;
}


/* the encoded string MUST have the length of the sequence at position 0!!! */
PUBLIC unsigned int
maximumMatching(const char *string)
{
  unsigned int          max;

  vrna_fold_compound_t  *fc = vrna_fold_compound(string, NULL, VRNA_OPTION_DEFAULT);

  max = (unsigned int)vrna_maximum_matching(fc);
  vrna_fold_compound_free(fc);

  return max;
}


/* the encoded string MUST have the length of the sequence at position 0!!! */
PUBLIC unsigned int *
maximumMatchingConstraint(const char  *string,
                          short       *ptable)
{
  unsigned int  i, j, l, length, max = 0;
  unsigned int  *mm;           /* holds maximum matching on subsequence [i,j] */
  short         *encodedString  = encode_sequence(string, 0);
  int           *iindx          = vrna_idx_row_wise((unsigned)encodedString[0]);

  make_pair_matrix();
  length  = (unsigned int)encodedString[0];
  mm      = (unsigned int *)vrna_alloc(sizeof(unsigned int) * ((length * (length + 1)) / 2 + 2));
  for (j = 1; j <= length; j++)
    for (i = (j > TURN ? (j - TURN) : 1); i < j; i++)
      mm[iindx[i] - j] = 0;
  for (i = length - TURN - 1; i > 0; i--)
    for (j = i + TURN + 1; j <= length; j++) {
      max = mm[iindx[i] - j + 1];
      for (l = j - TURN - 1; l >= i; l--) {
        if (pair[encodedString[l]][encodedString[j]])
          if ((unsigned int)ptable[l] != j)
            max = MAX2(max, ((l > i) ? mm[iindx[i] - l + 1] : 0) + 1 + mm[iindx[l + 1] - j + 1]);
      }
      mm[iindx[i] - j] = max;
    }
  free(iindx);
  free(encodedString);
  return mm;
}


/* the encoded string MUST have the length of the sequence at position 0!!! */
PUBLIC unsigned int *
maximumMatching2Constraint(const char *string,
                           short      *ptable,
                           short      *ptable2)
{
  unsigned int  i, j, l, length, max = 0;
  unsigned int  *mm;           /* holds maximum matching on subsequence [i,j] */
  short         *encodedString  = encode_sequence(string, 0);
  int           *iindx          = vrna_idx_row_wise((unsigned)encodedString[0]);

  make_pair_matrix();
  length  = (unsigned int)encodedString[0];
  mm      = (unsigned int *)vrna_alloc(sizeof(unsigned int) * ((length * (length + 1)) / 2 + 2));
  for (j = 1; j <= length; j++)
    for (i = (j > TURN ? (j - TURN) : 1); i < j; i++)
      mm[iindx[i] - j] = 0;
  for (i = length - TURN - 1; i > 0; i--)
    for (j = i + TURN + 1; j <= length; j++) {
      max = mm[iindx[i] - j + 1];
      for (l = j - TURN - 1; l >= i; l--) {
        if (pair[encodedString[l]][encodedString[j]])
          if ((unsigned int)ptable[l] != j && ptable2[l] != j)
            max = MAX2(max, ((l > i) ? mm[iindx[i] - l + 1] : 0) + 1 + mm[iindx[l + 1] - j + 1]);
      }
      mm[iindx[i] - j] = max;
    }
  free(iindx);
  free(encodedString);
  return mm;
}
