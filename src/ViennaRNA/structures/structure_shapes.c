/*
 *  ViennaRNA/structures/shapes.c
 *
 *  Functions to convert secondary structures into abstract shapes
 *
 *  c  Ronny Lorenz
 *     ViennaRNA package
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdio.h>
#include <string.h>

#include "ViennaRNA/utils/basic.h"
#include "ViennaRNA/structures/shapes.h"

#ifdef __GNUC__
# define INLINE inline
#else
# define INLINE
#endif


#define SHAPES_OPEN   '['
#define SHAPES_CLOSE  ']'
#define SHAPES_UP     '_'

struct shrep {
  struct shrep  *pred;
  struct shrep  *succ;
  char          character;
};


/*
 #################################
 # PRIVATE FUNCTION DECLARATIONS #
 #################################
 */
PRIVATE INLINE struct shrep *
shrep_insert_after(struct shrep *target,
                   char         character);


PRIVATE INLINE struct shrep *
shrep_insert_before(struct shrep  *target,
                    char          character);


PRIVATE INLINE void
shrep_concat(struct shrep **target,
             struct shrep *suffix);


PRIVATE INLINE char *
db2shapes(const char    *structure,
          unsigned int  level);


PRIVATE INLINE char *
db2shapes_pt(const short  *pt,
             unsigned int n,
             unsigned int level);


PRIVATE INLINE struct shrep *
get_shrep(const short   *pt,
          unsigned int  start,
          unsigned int  end,
          unsigned int  level);


/*
 #################################
 # BEGIN OF FUNCTION DEFINITIONS #
 #################################
 */
PUBLIC char *
vrna_abstract_shapes(const char   *structure,
                     unsigned int level)
{
  if (structure) {
    if (level > 5)
      level = 5;

    return db2shapes(structure, level);
  }

  return NULL;
}


PUBLIC char *
vrna_abstract_shapes_pt(const short   *pt,
                        unsigned int  level)
{
  if (pt) {
    if (level > 5)
      level = 5;

    return db2shapes_pt(pt, (unsigned int)pt[0], level);
  }

  return NULL;
}


/*
 #####################################
 # BEGIN OF STATIC HELPER FUNCTIONS  #
 #####################################
 */
PRIVATE INLINE struct shrep *
shrep_insert_after(struct shrep *target,
                   char         character)
{
  struct shrep *entry, *ptr;

  entry             = (struct shrep *)vrna_alloc(sizeof(struct shrep));
  entry->character  = character;
  entry->pred       = NULL;
  entry->succ       = NULL;

  if (target) {
    for (ptr = target; ptr->succ; ptr = ptr->succ);
    entry->pred = ptr;
    ptr->succ   = entry;
  }

  return entry;
}


PRIVATE INLINE struct shrep *
shrep_insert_before(struct shrep  *target,
                    char          character)
{
  struct shrep *entry, *ptr;

  entry             = (struct shrep *)vrna_alloc(sizeof(struct shrep));
  entry->character  = character;
  entry->pred       = NULL;
  entry->succ       = NULL;

  if (target) {
    for (ptr = target; ptr->pred; ptr = ptr->pred);
    entry->succ = ptr;
    ptr->pred   = entry;
  }

  return entry;
}


PRIVATE INLINE void
shrep_concat(struct shrep **target,
             struct shrep *suffix)
{
  struct shrep *ptr_s, *ptr_p;

  ptr_p = *target;
  ptr_s = suffix;

  if (ptr_p)
    for (; ptr_p->succ; ptr_p = ptr_p->succ);

  if (ptr_s)
    for (; ptr_s->pred; ptr_s = ptr_s->pred);

  if (!ptr_p) {
    *target = ptr_s;
  } else if (ptr_s) {
    ptr_p->succ = ptr_s;
    ptr_s->pred = ptr_p;
  }
}


PRIVATE INLINE char *
db2shapes(const char    *structure,
          unsigned int  level)
{
  unsigned int  n;
  char          *SHAPE;
  short         *pt;

  n   = strlen(structure);
  pt  = vrna_ptable(structure);

  SHAPE = db2shapes_pt(pt, n, level);

  free(pt);

  return SHAPE;
}


PRIVATE INLINE char *
db2shapes_pt(const short  *pt,
             unsigned int n,
             unsigned int level)
{
  unsigned int  start;
  struct shrep  *ptr, *ptr2;
  char          *SHAPE;

  SHAPE = NULL;
  start = 1;

  ptr = get_shrep(pt, start, n, level);

  if (ptr) {
    SHAPE = (char *)vrna_alloc(sizeof(char) * (n + 1));

    for (; ptr->pred; ptr = ptr->pred);

    for (n = 0; ptr; n++) {
      SHAPE[n]  = ptr->character;
      ptr2      = ptr;
      ptr       = ptr->succ;
      free(ptr2);
    }

    free(ptr);

    SHAPE     = vrna_realloc(SHAPE, sizeof(char) * (n + 1));
    SHAPE[n]  = '\0';
  }

  return SHAPE;
}


PRIVATE INLINE struct shrep *
get_shrep(const short   *pt,
          unsigned int  start,
          unsigned int  end,
          unsigned int  level)
{
  struct shrep  *shape_string, *substring, *sub5, *sub3;
  unsigned int  components, i, k, l, ext, u5, u3, cnt;
  unsigned char no_stack, bulge;

  shape_string = NULL;

  /* find out if there is more than one component */
  for (i = start, components = 0; i <= end; i++) {
    if (i < pt[i]) {
      components++;
      i = pt[i];
      if (components > 1)
        break;
    }
  }

  /* find out if we process the exterior loop */
  for (i = start - 1, ext = 1; i > 0; i--) {
    if (pt[i] > end) {
      ext = 0;
      break;
    }
  }

  for (; start <= end; start++) {
    if (start < pt[start]) {
      /* we have base pair (start, pt[start]) */
      substring = sub5 = sub3 = NULL;
      k         = start + 1;
      l         = pt[start] - 1;

      while (1) {
        u5  = 0;
        u3  = 0;

        for (; pt[k] == 0; k++, u5++);  /* skip unpaired 5' side */
        for (; pt[l] == 0; l--, u3++);  /* skip unpaired 3' side */

        if (k >= l) {
          /* hairpin loop */
          if (level == 0)
            for (cnt = 0; cnt < u5; cnt++)
              substring = shrep_insert_after(substring, SHAPES_UP);

          break;
        } else {
          if (pt[k] != l) {
            /* multi branch loop or external loop */
            substring = get_shrep(pt, k, l, level);
            if (level == 1) {
              if (u5)
                sub5 = shrep_insert_after(sub5, SHAPES_UP);

              if (u3)
                sub3 = shrep_insert_before(sub3, SHAPES_UP);
            } else if (level == 0) {
              for (cnt = 0; cnt < u5; cnt++)
                sub5 = shrep_insert_after(sub5, SHAPES_UP);

              for (cnt = 0; cnt < u3; cnt++)
                sub3 = shrep_insert_before(sub3, SHAPES_UP);
            }

            break;
          } else {
            /* internal loop with enclosed pair (k,l) */
            no_stack  = (u5 + u3 > 0) ? 1 : 0;
            bulge     = ((u5 > 0) && (u3 > 0)) ? 1 : 0;

            switch (level) {
              case 4:
                if (bulge) {
                  sub5  = shrep_insert_after(sub5, SHAPES_OPEN);
                  sub3  = shrep_insert_before(sub3, SHAPES_CLOSE);
                }

                break;

              case 3:
                if (no_stack) {
                  sub5  = shrep_insert_after(sub5, SHAPES_OPEN);
                  sub3  = shrep_insert_before(sub3, SHAPES_CLOSE);
                }

                break;

              case 2: /* fall through */
              case 1:
                if (no_stack) {
                  if (u5)
                    sub5 = shrep_insert_after(sub5, SHAPES_UP);

                  if (u3)
                    sub3 = shrep_insert_before(sub3, SHAPES_UP);

                  sub5  = shrep_insert_after(sub5, SHAPES_OPEN);
                  sub3  = shrep_insert_before(sub3, SHAPES_CLOSE);
                }

                break;

              case 0:
                for (unsigned int cnt = 0; cnt < u5; cnt++)
                  sub5 = shrep_insert_after(sub5, SHAPES_UP);

                for (unsigned int cnt = 0; cnt < u3; cnt++)
                  sub3 = shrep_insert_before(sub3, SHAPES_UP);

                sub5  = shrep_insert_after(sub5, SHAPES_OPEN);
                sub3  = shrep_insert_before(sub3, SHAPES_CLOSE);
                break;

              default:
                break;
            }

            k++;
            l--;
          }
        }
      }

      sub5 = shrep_insert_before(sub5, SHAPES_OPEN);
      shrep_concat(&sub5, substring);
      shrep_concat(&sub5, sub3);
      sub5 = shrep_insert_after(sub5, SHAPES_CLOSE);
      shrep_concat(&shape_string, sub5);

      start = pt[start];
    } else if ((pt[start] == 0) && (level < 3)) {
      if (level == 0) {
        shape_string = shrep_insert_after(shape_string, SHAPES_UP);
      } else if ((level == 1) ||
                 ((components < 2) && (ext == 0))) {
        shape_string = shrep_insert_after(shape_string, SHAPES_UP);
        for (; (start <= end) && (pt[start] == 0); start++);
        start--;
      }
    }
  }

  return shape_string;
}
