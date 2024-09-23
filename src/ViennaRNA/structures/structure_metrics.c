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

#include "ViennaRNA/utils/basic.h"
#include "ViennaRNA/params/basic.h"
#include "ViennaRNA/structures/metrics.h"

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
PUBLIC int
vrna_bp_distance_pt(const short *pt1,
                    const short *pt2)
{
  /*
   * dist = {number of base pairs in one structure but not in the other}
   * same as edit distance with pair_open pair_close as move set
   */
  int   dist;
  short i, l;

  dist = 0;

  if (pt1 && pt2) {
    l = (pt1[0] < pt2[0]) ? pt1[0] : pt2[0]; /* minimum of the two lengths */

    for (i = 1; i <= l; i++)
      if (pt1[i] != pt2[i]) {
        if (pt1[i] > i)
          dist++;

        if (pt2[i] > i)
          dist++;
      }
  }

  return dist;
}


PUBLIC int
vrna_bp_distance(const char *str1,
                 const char *str2)
{
  int   dist = 0;
  short *pt1, *pt2;

  pt1 = vrna_ptable(str1);
  pt2 = vrna_ptable(str2);

  dist = vrna_bp_distance_pt(pt1, pt2);

  free(pt1);
  free(pt2);

  return dist;
}


PUBLIC double
vrna_dist_mountain(const char   *str1,
                   const char   *str2,
                   unsigned int p)
{
  short         *pt1, *pt2;
  unsigned int  i, n;
  double        distance, w, *f1, *f2;

  distance  = -1.;
  f1        = NULL;
  f2        = NULL;

  if ((str1) && (str2)) {
    n = strlen(str1);

    if (n != strlen(str2)) {
      vrna_log_warning("vrna_dist_mountain: input structures have unequal lengths!");
      return distance;
    }

    pt1 = vrna_ptable(str1);
    pt2 = vrna_ptable(str2);
    f1  = (double *)vrna_alloc(sizeof(double) * (n + 1));
    f2  = (double *)vrna_alloc(sizeof(double) * (n + 1));

    /* count (mountain)heights for positions 1 <= i <= n */
    for (w = 0., i = 1; i <= n; i++) {
      if (pt1[i] == 0)
        continue;

      if (pt1[i] > i)
        w += 1. / (double)(pt1[i] - i);
      else
        w -= 1. / (double)(i - pt1[i]);

      f1[i] = w;
    }

    for (w = 0., i = 1; i <= n; i++) {
      if (pt2[i] == 0)
        continue;

      if (pt2[i] > i)
        w += 1. / (double)(pt2[i] - i);
      else
        w -= 1. / (double)(i - pt2[i]);

      f2[i] = w;
    }


    /* finally, compute L_p-norm */
    for (distance = 0., i = 1; i <= n; i++)
      distance += pow(fabs(f1[i] - f2[i]), (double)p);

    distance = pow(distance, 1. / (double)p);

    free(pt1);
    free(pt2);
    free(f1);
    free(f2);
  }

  return distance;
}


#ifndef VRNA_DISABLE_BACKWARD_COMPATIBILITY

/*
 * ###########################################
 * # deprecated functions below              #
 *###########################################
 */

PUBLIC int
bp_distance(const char  *str1,
            const char  *str2)
{
  return vrna_bp_distance(str1, str2);
}

#endif
