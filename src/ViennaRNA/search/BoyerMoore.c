/*
 *  ViennaRNA/search/BoyerMoore.c
 *
 *  Variations of the Boyer-Moore search algorithms
 *
 *  (c) 2018, Ronny Lorenz
 *
 *  ViennaRNA package
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "ViennaRNA/utils/basic.h"
#include "ViennaRNA/utils/log.h"
#include "ViennaRNA/search/BoyerMoore.h"

/*
 #################################
 # GLOBAL VARIABLES              #
 #################################
 */

/*
 #################################
 # PRIVATE VARIABLES and STRUCTS #
 #################################
 */

/*
 #################################
 # PRIVATE FUNCTION DECLARATIONS #
 #################################
 */
PRIVATE size_t *
get_BM_BCT_num(const unsigned int *needle,
               size_t             needle_size,
               unsigned int       maxnum);


PRIVATE size_t *
get_BM_BCT(const char *needle,
           size_t     needle_size);


PRIVATE const unsigned int *
BoyerMooreHorspool_num(const unsigned int *needle,
                       size_t             needle_size,
                       const unsigned int *haystack,
                       size_t             haystack_size,
                       size_t             start,
                       size_t             *bad_chars,
                       unsigned char      cyclic);


PRIVATE const char *
BoyerMooreHorspool(const char     *needle,
                   size_t         needle_size,
                   const char     *haystack,
                   size_t         haystack_size,
                   size_t         start,
                   size_t         *bad_chars,
                   unsigned char  cyclic);


/*
 #################################
 # BEGIN OF FUNCTION DEFINITIONS #
 #################################
 */
PUBLIC const unsigned int *
vrna_search_BMH_num(const unsigned int  *needle,
                    size_t              needle_size,
                    const unsigned int  *haystack,
                    size_t              haystack_size,
                    size_t              start,
                    size_t              *badchars,
                    unsigned char       cyclic)
{
  const unsigned int  *hit;
  unsigned int        max;
  size_t              *bc, i;

  if ((!needle) || (!haystack) || (start > haystack_size))
    return NULL;

  bc = badchars;

  /* create bad character table in case none is supplied */
  if (!bc) {
    max = needle[0];
    for (i = 1; i < needle_size; i++)
      max = MAX2(max, needle[i]);

    for (i = 1; i < haystack_size; i++)
      max = MAX2(max, haystack[i]);

    bc = get_BM_BCT_num(needle, needle_size, max);
  }

  /* perform actual search */
  hit = BoyerMooreHorspool_num(needle,
                               needle_size,
                               haystack,
                               haystack_size,
                               start,
                               bc,
                               cyclic);

  if (bc != badchars)
    free(bc);

  return hit;
}


PUBLIC const char *
vrna_search_BMH(const char    *needle,
                size_t        needle_size,
                const char    *haystack,
                size_t        haystack_size,
                size_t        start,
                size_t        *badchars,
                unsigned char cyclic)
{
  const char  *hit;
  size_t      *bc;

  if ((!needle) || (!haystack) || (start > haystack_size))
    return NULL;

  bc = badchars;

  /* create bad character table in case none is supplied */
  if (!bc)
    bc = get_BM_BCT(needle, needle_size);

  /* perform actual search */
  hit = BoyerMooreHorspool(needle,
                           needle_size,
                           haystack,
                           haystack_size,
                           start,
                           bc,
                           cyclic);

  if (bc != badchars)
    free(bc);

  return hit;
}


PUBLIC size_t *
vrna_search_BM_BCT_num(const unsigned int *pattern,
                       size_t             pattern_size,
                       unsigned int       num_max)
{
  if (!pattern)
    return NULL;

  return get_BM_BCT_num(pattern, pattern_size, num_max);
}


PUBLIC size_t *
vrna_search_BM_BCT(const char *pattern)
{
  if (!pattern)
    return NULL;

  size_t pattern_size = strlen(pattern);

  return get_BM_BCT(pattern, pattern_size);
}


/*
 #################################
 # STATIC helper functions below #
 #################################
 */
PRIVATE size_t *
get_BM_BCT_num(const unsigned int *needle,
               size_t             needle_size,
               unsigned int       maxnum)
{
  size_t *table, i;

  table = vrna_alloc(sizeof(size_t) * (maxnum + 2));

  /* store maximum element value at position 0 */
  table[0] = maxnum;

  /* use remainder of array for actual bad character table */
  for (i = 1; i <= maxnum + 1; i++)
    table[i] = needle_size;

  for (i = 0; i < needle_size - 1; i++)
    table[needle[i] + 1] = needle_size - i - 1;

  return table;
}


PRIVATE size_t *
get_BM_BCT(const char *needle,
           size_t     needle_size)
{
  size_t *table, i;

  table = vrna_alloc(sizeof(size_t) * (127 + 2));

  /* store maximum element value at position 0 */
  table[0] = 127;

  /* use remainder of array for actual bad character table */
  for (i = 1; i <= 127 + 1; i++)
    table[i] = needle_size;

  for (i = 0; i < needle_size - 1; i++)
    table[needle[i] + 1] = needle_size - i - 1;

  return table;
}


/*
 * An implementation of the Boyer-Moore-Horspool search algorithm
 *
 * Returns a pointer to the first occurence of needle within haystack.
 * If needle can't be found the function returns NULL
 *
 * Note, the Bad Character Table BCT must contain the maximum number
 * representation of the objects in needle and haystack, i.e. the
 * alphabet size. Actual skip data then starts at T[1] for element
 * with value 0
 */
#define BMH { \
    hit     = NULL; \
    shift   = start; \
    margin  = (cyclic) ? 0 : needle_size; \
    max     = bad_chars[0]; \
    /* pop first element since the Bad Character Table starts at position 1 */ \
    bad_chars++; \
    /* main loop - go through haystack */ \
    while (shift + margin < haystack_size) { \
      /* \
       *  matching loop, note that we allow for possibly wrapping the \
       *  pattern around the haystack \
       */\
      for (i = needle_size - 1; \
           haystack[(shift + i) % haystack_size] == needle[i]; \
           i--) { \
        if (i == 0) \
        return haystack + shift; \
      } \
      val = haystack[(shift + needle_size - 1) % haystack_size]; \
      if (val > max) { \
        vrna_log_warning("vrna_search_BMH: " \
                         "haystack value %d at hit %d " \
                         "out of bad character table range [%d : %d]\n" \
                         "Aborting search...", \
                         (shift + needle_size - 1) % haystack_size, \
                         val, \
                         0, \
                         max); \
        return NULL; \
      } \
      shift += bad_chars[(size_t)val]; \
    } \
}

PRIVATE const unsigned int *
BoyerMooreHorspool_num(const unsigned int *needle,
                       size_t             needle_size,
                       const unsigned int *haystack,
                       size_t             haystack_size,
                       size_t             start,
                       size_t             *bad_chars,
                       unsigned char      cyclic)
{
  const unsigned int  *hit;
  unsigned int        val, max;
  size_t              i, shift, margin;

  /* empty pattern matches element in haystack */
  if (needle_size == 0)
    return haystack;

  /* empty haystack can't contain any pattern */
  if (haystack_size == 0)
    return NULL;

  /* haystack mustn't be shorter than needle */
  if (haystack_size < needle_size)
    return NULL;

  /* begin actual algorithm */
  BMH;

  return hit;
}


PRIVATE const char *
BoyerMooreHorspool(const char     *needle,
                   size_t         needle_size,
                   const char     *haystack,
                   size_t         haystack_size,
                   size_t         start,
                   size_t         *bad_chars,
                   unsigned char  cyclic)
{
  const char  *hit;
  char        val, max;
  size_t      i, shift, margin;

  /* empty pattern matches element in haystack */
  if (!needle)
    return haystack;

  /* empty pattern matches element in haystack */
  if (needle_size == 0)
    return haystack;

  /* empty haystack can't contain any pattern */
  if (haystack_size == 0)
    return NULL;

  /* haystack mustn't be shorter than needle */
  if (haystack_size < needle_size)
    return NULL;

  /* begin actual algorithm */
  BMH;

  return hit;
}
