/*
 *  ViennaRNA/utils/structures.c
 *
 *  Various functions to convert, parse, encode secondary structures
 *
 *  c  Ivo L Hofacker, Walter Fontana, Ronny Lorenz
 *              Vienna RNA package
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <limits.h>

#include "ViennaRNA/fold_vars.h"
#include "ViennaRNA/utils/basic.h"
#include "ViennaRNA/params/basic.h"
#include "ViennaRNA/structures/utils.h"

#ifdef __GNUC__
# define INLINE inline
#else
# define INLINE
#endif

/*
 #################################
 # PRIVATE FUNCTION DECLARATIONS #
 #################################
 */

/*
 #################################
 # BEGIN OF FUNCTION DEFINITIONS #
 #################################
 */

/* get a matrix containing the number of basepairs of a reference structure for each interval [i,j] with i<j
 *  access it via iindx!!!
 */
PUBLIC unsigned int *
vrna_refBPcnt_matrix(const short  *reference_pt,
                     unsigned int turn)
{
  unsigned int  i, j, k, ij, length;
  int           *iindx;
  unsigned int  *array;
  unsigned int  size;

  length  = (unsigned int)reference_pt[0];
  size    = ((length + 1) * (length + 2)) / 2;
  iindx   = vrna_idx_row_wise(length);
  array   = (unsigned int *)vrna_alloc(sizeof(unsigned int) * size);
  /* matrix containing number of basepairs of reference structure1 in interval [i,j] */;
  for (k = 0; k <= turn; k++)
    for (i = 1; i <= length - k; i++) {
      j         = i + k;
      ij        = iindx[i] - j;
      array[ij] = 0;
    }

  for (i = length - turn - 1; i >= 1; i--)
    for (j = i + turn + 1; j <= length; j++) {
      int bps;
      ij  = iindx[i] - j;
      bps = array[ij + 1];
      if ((i <= (unsigned int)reference_pt[j]) && ((unsigned int)reference_pt[j] < j))
        bps++;

      array[ij] = bps;
    }
  free(iindx);
  return array;
}


PUBLIC unsigned int *
vrna_refBPdist_matrix(const short   *pt1,
                      const short   *pt2,
                      unsigned int  turn)
{
  unsigned int  *array;
  unsigned int  n, size, i, j, ij, d;

  n     = (unsigned int)pt1[0];
  size  = ((n + 1) * (n + 2)) / 2;
  array = (unsigned int *)vrna_alloc(sizeof(unsigned int) * size);
  int           *iindx = vrna_idx_row_wise(n);

  for (i = n - turn - 1; i >= 1; i--) {
    d = 0;
    for (j = i + turn + 1; j <= n; j++) {
      ij  = iindx[i] - j;
      d   = array[ij + 1];
      if (pt1[j] != pt2[j]) {
        if (i <= (unsigned int)pt1[j] && (unsigned int)pt1[j] < j)
          /* we got an additional base pair in reference structure 1 */
          d++;

        if (i <= (unsigned int)pt2[j] && (unsigned int)pt2[j] < j)
          /* we got another base pair in reference structure 2 */
          d++;
      }

      array[ij] = d;
    }
  }
  free(iindx);
  return array;
}


/*
 #####################################
 # BEGIN OF STATIC HELPER FUNCTIONS  #
 #####################################
 */

#ifndef VRNA_DISABLE_BACKWARD_COMPATIBILITY

PUBLIC unsigned int *
make_referenceBP_array(short        *reference_pt,
                       unsigned int turn)
{
  return vrna_refBPcnt_matrix((const short *)reference_pt, turn);
}


PUBLIC unsigned int *
compute_BPdifferences(short         *pt1,
                      short         *pt2,
                      unsigned int  turn)
{
  return vrna_refBPdist_matrix((const short *)pt1, (const short *)pt2, turn);
}

#endif
