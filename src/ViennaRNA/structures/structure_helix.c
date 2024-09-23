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

#include "ViennaRNA/utils/basic.h"
#include "ViennaRNA/structures/helix.h"

/*
 #################################
 # PRIVATE FUNCTION DECLARATIONS #
 #################################
 */
PRIVATE int
hx_cmp(const void *a,
       const void *b);


/*
 #################################
 # BEGIN OF FUNCTION DEFINITIONS #
 #################################
 */
PUBLIC vrna_hx_t *
vrna_hx_from_ptable(short *pt)
{
  int       i, k, n, l, s, *stack;
  vrna_hx_t *list;

  list = NULL;

  if (pt) {
    n     = pt[0];
    l     = 0;
    s     = 1;
    list  = (vrna_hx_t *)vrna_alloc(sizeof(vrna_hx_t) * (n / 2 + 2));
    stack = (int *)vrna_alloc(sizeof(int) * (n / 2 + 2));

    stack[s] = 1;

    do {
      for (i = stack[s--]; i <= n; i++) {
        if (pt[i] > (short)i) {
          /* found a base pair */
          k = i;
          /* go through stack */
          for (; pt[k + 1] == pt[k] - 1; k++);
          list[l].start   = i;
          list[l].end     = pt[i];
          list[l].length  = k - i + 1;
          list[l].up5     = list[l].up3 = 0;
          l++;
          stack[++s]  = pt[i] + 1;
          /*
           *  only push enclosed part to stack if there is anything
           *  actually enclosed by the helix
           */
          if (k != pt[i])
            stack[++s]  = k + 1;
          break;
        } else if (pt[i]) {
          /* end of region */
          break;
        }
      }
    } while (s > 0);

    list          = vrna_realloc(list, (l + 1) * sizeof(vrna_hx_t));
    list[l].start = list[l].end = list[l].length = list[l].up5 = list[l].up3 = 0;

    free(stack);
  }

  return list;
}


PUBLIC vrna_hx_t *
vrna_hx_merge(const vrna_hx_t *list,
              int             maxdist)
{
  int       merged, i, j, s, neighbors, n;
  vrna_hx_t *merged_list;

  merged_list = NULL;

  if (list) {
    for (n = 0; list[n].length > 0; n++); /* check size of list */

    merged_list = (vrna_hx_t *)vrna_alloc(sizeof(vrna_hx_t) * (n + 1));
    memcpy(merged_list, list, sizeof(vrna_hx_t) * (n + 1));

    qsort(merged_list, n, sizeof(vrna_hx_t), hx_cmp);

    s = n + 1;

    do {
      merged = 0;
      for (i = 1; merged_list[i].length > 0; i++) {
        /*
         * GOAL: merge two consecutive helices i and i-1, if i-1
         * subsumes i, and not more than i
         */

        /* 1st, check for neighboring helices or any crossing
         * helices forming a pseudoknot
         */
        unsigned char attempt_merge = 1;
        /* go through list of remaining helices and see whether helix[i - 1]
         * (i) encloses another helix that also encloses helix[i]
         * (ii) encloses another helix apart from [i]
         * (iii) is enclosing parts of another helix forming a pseudoknot, or
         * (iv) helix[i] forms a pseudoknot with another helix not enclosed by helix[i]
         */
        unsigned int is5, ie5, is3, ie3;
        unsigned int i1s5, i1e5, i1s3, i1e3;
        i1s5 = merged_list[i - 1].start;
        i1e3 = merged_list[i - 1].end;
        i1e5 = i1s5 + merged_list[i - 1].length + merged_list[i - 1].up5 - 1;
        i1s3 = i1e3 + 1 - merged_list[i - 1].length - merged_list[i - 1].up3;
        is5 = merged_list[i].start;
        ie3 = merged_list[i].end;
        ie5 = is5 + merged_list[i].length + merged_list[i].up5 - 1;
        is3 = ie3 + 1 - merged_list[i].length - merged_list[i].up3;

        /* do not merge if helices are not nested */
        if ((is5 < i1e5) || (ie3 > i1s3))
          continue;

        /* check for conflicting helices */
        for (j = 0; j < i - 1; j++) {
          /*
           * upstream starting helices can potentially form pseudoknots into the
           * region between the two helices we are about to merge
           */
          unsigned int s5, e5, s3, e3;
          s5 = merged_list[j].start;
          e3 = merged_list[j].end;
          e5 = s5 + merged_list[j].length + merged_list[j].up5 - 1;
          s3 = e3 + 1 - merged_list[j].length - merged_list[j].up3;

          /* do not allow for cases:
           * [.(.(..).].) or [.(.].(..).)
           */
          if ((e5 < i1s5) &&
              (i1e5 < s3) &&
              ((e3 < is5) || ((e3 < i1s3) && (ie3 < s3)))) {
            attempt_merge = 0;
            break;
          }
        }

        if (!attempt_merge)
          continue;

        for (j = i + 1; merged_list[j].length > 0; j++) {
          unsigned int s5, e5, s3, e3;
          s5 = merged_list[j].start;
          e3 = merged_list[j].end;
          e5 = s5 + merged_list[j].length + merged_list[j].up5 - 1;
          s3 = e3 + 1 - merged_list[j].length - merged_list[j].up3;

          if ((i1e3 < s5) ||
              (e3 < i1s5))
            continue;

          if (e5 < i1s5) {
            if (((i1e5 < s3) && (e3 < is5)) ||
                ((ie3 < s3) && (e3 < i1s3))) {
              attempt_merge = 0;
              break;
            }
            continue;
          }

          if ((i1e5 < s5) &&
              (e5 < is5)) {
            if ((ie5 < s3) && (e3 < is3)) {
              attempt_merge = 0;
              break;
            } else if (i1e3 < s3) {
              attempt_merge = 0;
              break;
            } else if (e3 < is5) {
              attempt_merge = 0;
              break;
            } else if ((ie3 < s3) && (e3 < i1s3)) {
              attempt_merge = 0;
              break;
            }
            continue;
          }      

          if ((ie5 < s5) && (e5 < is3) && (ie3 < s3) && (e3 < i1s3)) {
            attempt_merge = 0;
            break;
          }

          if ((ie3 < s5) && (e5 < i1s3)) {
            if ((e3 < i1s3) ||
                (i1e3 < s3)) {
              attempt_merge = 0;
              break;
            }
          }
        }

        if (!attempt_merge)
          continue;

        /* check if we may merge i with i-1 */
        if ((i1e5 < is5) && (ie3 < i1s3)) {
          merged_list[i - 1].up5 += merged_list[i].up5
                                    + (is5 - i1e5 - 1);
          merged_list[i - 1].up3 += merged_list[i].up3
                                    + (i1s3 - ie3 - 1);
          merged_list[i - 1].length += merged_list[i].length;
          /* splice out helix i */
          memmove(merged_list + i, merged_list + i + 1, sizeof(vrna_hx_t) * (n - i));
          s--;
          merged = 1;
          break;
        }
      }
    } while (merged);

    merged_list = vrna_realloc(merged_list, sizeof(vrna_hx_t) * s);
  }

  return merged_list;
}


/*
 #####################################
 # BEGIN OF STATIC HELPER FUNCTIONS  #
 #####################################
 */
PRIVATE int
hx_cmp(const void *a,
       const void *b)
{
  if (((vrna_hx_t *)a)->start > ((vrna_hx_t *)b)->start)
    return 1;

  if (((vrna_hx_t *)a)->start < ((vrna_hx_t *)b)->start)
    return -1;

  return 0;
}


