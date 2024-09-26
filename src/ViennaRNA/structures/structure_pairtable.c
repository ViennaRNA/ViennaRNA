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
#include "ViennaRNA/utils/log.h"
#include "ViennaRNA/structures/mea.h"
#include "ViennaRNA/structures/pairtable.h"

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

PRIVATE INLINE int
extract_pairs(short       *pt,
              const char  *structure,
              const char  *pair);


/*
 #################################
 # BEGIN OF FUNCTION DEFINITIONS #
 #################################
 */
PUBLIC short *
vrna_ptable(const char *structure)
{
  return vrna_ptable_from_string(structure, VRNA_BRACKETS_RND);
}


PUBLIC short *
vrna_pt_pk_get(const char *structure)
{
  return vrna_ptable_from_string(structure, VRNA_BRACKETS_RND | VRNA_BRACKETS_SQR);
}


PUBLIC short *
vrna_pt_snoop_get(const char *structure)
{
  return vrna_ptable_from_string(structure, VRNA_BRACKETS_ANG);
}


PUBLIC short *
vrna_pt_ali_get(const char *structure)
{
  return vrna_ptable_from_string(structure,
                                 VRNA_BRACKETS_RND | VRNA_BRACKETS_ANG | VRNA_BRACKETS_SQR);
}


PUBLIC short *
vrna_ptable_copy(const short *pt)
{
  short *table = (short *)vrna_alloc(sizeof(short) * (pt[0] + 2));

  memcpy(table, pt, sizeof(short) * (pt[0] + 2));
  return table;
}


PUBLIC int *
vrna_loopidx_from_ptable(const short *pt)
{
  /* number each position by which loop it belongs to (positions start
   * at 1) */
  int i, hx, l, nl;
  int length;
  int *stack  = NULL;
  int *loop   = NULL;

  length  = pt[0];
  stack   = (int *)vrna_alloc(sizeof(int) * (length + 1));
  loop    = (int *)vrna_alloc(sizeof(int) * (length + 2));
  hx      = l = nl = 0;

  for (i = 1; i <= length; i++) {
    if ((pt[i] != 0) && (i < pt[i])) {
      /* ( */
      nl++;
      l           = nl;
      stack[hx++] = i;
    }

    loop[i] = l;

    if ((pt[i] != 0) && (i > pt[i])) {
      /* ) */
      --hx;
      if (hx > 0)
        l = loop[stack[hx - 1]];  /* index of enclosing loop   */
      else
        l = 0;                    /* external loop has index 0 */

      if (hx < 0) {
        vrna_log_warning("vrna_loopidx_from_ptable: "
                             "unbalanced brackets in make_pair_table");
        free(stack);
        return NULL;
      }
    }
  }
  loop[0] = nl;
  free(stack);
  return loop;
}


PUBLIC short *
vrna_pt_pk_remove(const short   *ptable,
                  unsigned int  options)
{
  short *pt = NULL;

  if (ptable) {
    char          *mea_structure;
    unsigned int  i, j, n;
    vrna_ep_t     *pairs;

    n             = (unsigned int)ptable[0];
    mea_structure = (char *)vrna_alloc(sizeof(char) * (n + 1));
    pairs         = (vrna_ep_t *)vrna_alloc(sizeof(vrna_ep_t) * (n + 1));

    /* compose list of pairs to be used in MEA() function */
    for (j = 0, i = 1; i <= n; i++)
      if (ptable[i] > i) {
        pairs[j].i    = i;
        pairs[j].j    = ptable[i];
        pairs[j].p    = 1.;
        pairs[j].type = VRNA_PLIST_TYPE_BASEPAIR;
        j++;
      }

    pairs[j].i    = 0;
    pairs[j].j    = 0;
    pairs[j].p    = 0;
    pairs[j].type = VRNA_PLIST_TYPE_BASEPAIR;

    /* use MEA() implementation to remove crossing base pairs */
    memset(mea_structure, '.', n);

    (void)MEA(pairs, mea_structure, 2.0);

    /* convert dot-bracket structure to pair table */
    pt = vrna_ptable(mea_structure);

    free(mea_structure);
    free(pairs);
  }

  return pt;
}


PUBLIC short *
vrna_ptable_from_string(const char    *string,
                        unsigned int  options)
{
  char          pairs[3];
  short         *pt;
  unsigned int  i, n;

  n = strlen(string);

  if (n > SHRT_MAX) {
    vrna_log_warning("vrna_ptable_from_string: "
                         "Structure too long to be converted to pair table (n=%d, max=%d)",
                         n,
                         SHRT_MAX);
    return NULL;
  }

  pt    = (short *)vrna_alloc(sizeof(short) * (n + 2));
  pt[0] = (short)n;


  if ((options & VRNA_BRACKETS_RND) &&
      (!extract_pairs(pt, string, "()"))) {
    free(pt);
    return NULL;
  }

  if ((options & VRNA_BRACKETS_ANG) &&
      (!extract_pairs(pt, string, "<>"))) {
    free(pt);
    return NULL;
  }

  if ((options & VRNA_BRACKETS_CLY) &&
      (!extract_pairs(pt, string, "{}"))) {
    free(pt);
    return NULL;
  }

  if ((options & VRNA_BRACKETS_SQR) &&
      (!extract_pairs(pt, string, "[]"))) {
    free(pt);
    return NULL;
  }

  if (options & VRNA_BRACKETS_ALPHA) {
    for (i = 65; i < 91; i++) {
      pairs[0]  = (char)i;
      pairs[1]  = (char)(i + 32);
      pairs[2]  = '\0';
      if (!extract_pairs(pt, string, pairs)) {
        free(pt);
        return NULL;
      }
    }
  }

  return pt;
}


/*
 #####################################
 # BEGIN OF STATIC HELPER FUNCTIONS  #
 #####################################
 */

/* requires that pt[0] already contains the length of the string! */
PRIVATE INLINE int
extract_pairs(short       *pt,
              const char  *structure,
              const char  *pair)
{
  const char    *ptr;
  char          open, close;
  short         *stack;
  unsigned int  i, j, n;
  int           hx;

  n     = (unsigned int)pt[0];
  stack = (short *)vrna_alloc(sizeof(short) * (n + 1));

  open  = pair[0];
  close = pair[1];

  for (hx = 0, i = 1, ptr = structure; (i <= n) && (*ptr != '\0'); ptr++, i++) {
    if (*ptr == open) {
      stack[hx++] = i;
    } else if (*ptr == close) {
      j = stack[--hx];

      if (hx < 0) {
        vrna_log_warning("%s\nunbalanced brackets '%2s' found while extracting base pairs",
                             structure,
                             pair);
        free(stack);
        return 0;
      }

      pt[i] = j;
      pt[j] = i;
    }
  }

  free(stack);

  if (hx != 0) {
    vrna_log_warning("%s\nunbalanced brackets '%2s' found while extracting base pairs",
                         structure,
                         pair);
    return 0;
  }

  return 1; /* success */
}


#ifndef VRNA_DISABLE_BACKWARD_COMPATIBILITY

/*
 * ###########################################
 * # deprecated functions below              #
 *###########################################
 */

PUBLIC short *
make_pair_table(const char *structure)
{
  return vrna_ptable(structure);
}


PUBLIC short *
copy_pair_table(const short *pt)
{
  return vrna_ptable_copy(pt);
}


PUBLIC short *
make_pair_table_pk(const char *structure)
{
  return vrna_pt_pk_get(structure);
}


PUBLIC short *
make_pair_table_snoop(const char *structure)
{
  return vrna_pt_snoop_get(structure);
}


PUBLIC short *
alimake_pair_table(const char *structure)
{
  return vrna_pt_ali_get(structure);
}


PUBLIC int *
make_loop_index_pt(short *pt)
{
  return vrna_loopidx_from_ptable((const short *)pt);
}


#endif
