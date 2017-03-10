/*
 *  combinatorics.c
 *
 *  Various implementations that deal with combinatorial aspects
 *  of RNA/DNA folding
 *
 *  (c) 2016, Ronny Lorenz
 *
 *  Vienna RNA package
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdio.h>
#include <stdlib.h>

#include "ViennaRNA/utils.h"
#include "ViennaRNA/combinatorics.h"

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
struct perm_list {
  int               value;
  struct perm_list  *next, *prev;
};

struct necklace_content {
  int value;
  int count;
};


/*
 #################################
 # PRIVATE FUNCTION DECLARATIONS #
 #################################
 */
PRIVATE struct  perm_list *perm_list_insert(struct perm_list  *before,
                                            int               value);


PRIVATE struct  perm_list *perm_list_remove_val(struct perm_list  *head,
                                                int               value);


PRIVATE struct  perm_list *perm_list_head(struct perm_list *entry);


PRIVATE void    perm_list_destroy(struct perm_list *entry);


PRIVATE int     cmpfunc(const void  *a,
                        const void  *b);


PRIVATE void    sawada_fast_finish_perm(struct necklace_content *content,
                                        unsigned int            ***results,
                                        unsigned int            *result_count,
                                        unsigned int            *result_size,
                                        unsigned int            n);


PRIVATE void    sawada_fast(unsigned int            t,
                            unsigned int            p,
                            unsigned int            s,
                            struct necklace_content *content,
                            unsigned int            k,
                            unsigned int            *r,
                            struct perm_list        *a,
                            unsigned int            n,
                            unsigned int            ***results,
                            unsigned int            *result_count,
                            unsigned int            *result_size);


/*
 #################################
 # BEGIN OF FUNCTION DEFINITIONS #
 #################################
 */

/*
 * This is an implementation of Joe Sawada's
 * "Fast algorithm to generate necklaces with fixed content"
 * Theoretical Computer Science 301 (2003) 477 - 489
 */
PUBLIC unsigned int **
vrna_enumerate_necklaces(const unsigned int *entity_counts)
{
  unsigned int            n, i, *r, **result, result_count, result_size, num_entities;
  struct necklace_content *content;
  struct  perm_list       *a;

  num_entities = 0;

  if (entity_counts)
    for (i = 0; entity_counts[i] > 0; i++)
      num_entities++;

  /* count total number of strands */
  for (n = i = 0; i < num_entities; i++)
    n += entity_counts[i];

  /* first re-order entity_counts such that k - 1 has most occurences */
  content = (struct necklace_content *)vrna_alloc(sizeof(struct necklace_content) * num_entities);
  for (i = 0; i < num_entities; i++) {
    content[i].value  = i;
    content[i].count  = entity_counts[i];
  }
  qsort(content, num_entities, sizeof(struct necklace_content), cmpfunc);

  /* create a-list (available characters) */
  a = NULL;
  for (i = 0; i < num_entities; i++)
    a = perm_list_insert(a, i); /* inserts new element before a */

  /* create r-array */
  r = (unsigned int *)vrna_alloc(sizeof(unsigned int) * (n + 1));

  /* prepare result list */
  result        = NULL;
  result_count  = 0;
  result_size   = 20;

  result = (unsigned int **)vrna_alloc(sizeof(unsigned int *) * (result_size));
  for (i = 0; i < result_size; i++)
    result[i] = (unsigned int *)vrna_alloc(sizeof(unsigned int) * (n + 1));

  /* do initial step */
  for (i = 1; i <= n; i++)
    result[result_count][i] = content[num_entities - 1].value;

  result[result_count][1] = 0;
  content[0].count        = content[0].count - 1;
  if (content[0].count == 0)
    a = perm_list_remove_val(a, 0);

  /* now we iterate to retrieve full permutation(s) */
  sawada_fast(2, 1, 2, content, num_entities, r, a, n, &result, &result_count, &result_size);

  /* resize results list to actual requirements */
  for (i = result_count; i < result_size; i++)
    free(result[i]);
  result                = (unsigned int **)vrna_realloc(result, sizeof(unsigned int *) * (result_count + 1));
  result[result_count]  = NULL;

  /* cleanup memory */
  free(r);
  free(content);
  perm_list_destroy(a);

  return result;
}


/*
 #################################
 # STATIC helper functions below #
 #################################
 */
PRIVATE struct perm_list *
perm_list_insert(struct perm_list *before,
                 int              value)
{
  struct perm_list *new;

  new = (struct perm_list *)vrna_alloc(sizeof(struct perm_list));

  new->value  = value;
  new->next   = NULL;
  new->prev   = NULL;

  if (before) {
    new->prev     = before->prev;
    new->next     = before;
    before->prev  = new;
  }

  return new;
}


PRIVATE struct perm_list *
perm_list_remove_val(struct perm_list *head,
                     int              value)
{
  struct perm_list *ptr;

  ptr = head;

  while (ptr) {
    if (ptr->value == value) {
      if (ptr->prev)
        ptr->prev->next = ptr->next;
      else
        head = ptr->next;

      if (ptr->next)
        ptr->next->prev = ptr->prev;

      free(ptr);
      break;
    }

    ptr = ptr->next;
  }

  return head;
}


PRIVATE struct perm_list *
perm_list_head(struct perm_list *entry)
{
  struct perm_list *head;

  head = entry;
  if (head)
    while (head->prev)
      head = head->prev;

  return head;
}


PRIVATE void
perm_list_destroy(struct perm_list *entry)
{
  struct perm_list *head, *ptr;

  head = perm_list_head(entry);
  while (head) {
    ptr   = head;
    head  = ptr->next;
    free(ptr);
  }
}


PRIVATE int
cmpfunc(const void  *a,
        const void  *b)
{
  return ((struct necklace_content *)a)->count - ((struct necklace_content *)b)->count;
}


PRIVATE void
sawada_fast_finish_perm(struct necklace_content *content,
                        unsigned int            ***results,
                        unsigned int            *result_count,
                        unsigned int            *result_size,
                        unsigned int            n)
{
  unsigned int i;

  /* adjust results list size if necessary */
  if ((*result_count + 1) == (*result_size)) {
    *result_size  *= 1.2;
    (*results)    = (unsigned int **)vrna_realloc(*results, sizeof(unsigned int *) * (*result_size));
    for (i = *result_count + 1; i < *result_size; i++)
      (*results)[i] = (unsigned int *)vrna_alloc(sizeof(unsigned int) * (n + 1));
  }

  /* create new (next) entry, copy-over curent result, and transliterate current one */
  for (i = 1; i <= n; i++) {
    unsigned int v = (*results)[*result_count][i];
    (*results)[(*result_count) + 1][i]  = v;
    (*results)[(*result_count)][i]      = content[v].value;
  }
  (*result_count)++;
}


PRIVATE void
sawada_fast(unsigned int            t,
            unsigned int            p,
            unsigned int            s,
            struct necklace_content *content,
            unsigned int            k,
            unsigned int            *r,
            struct perm_list        *a,
            unsigned int            n,
            unsigned int            ***results,
            unsigned int            *result_count,
            unsigned int            *result_size)
{
  if (content[k - 1].count == n - t + 1) {
    if ((content[k - 1].count == r[t - p]) && (n % p == 0))
      sawada_fast_finish_perm(content, results, result_count, result_size, n);
    else if (content[k - 1].count > r[t - p])
      sawada_fast_finish_perm(content, results, result_count, result_size, n);
  } else if (content[0].count != (n - t + 1)) {
    unsigned int      *result = (*results)[*result_count];
    struct perm_list  *ptr, *before, *after;
    ptr = perm_list_head(a);
    unsigned int      j   = ptr->value;
    unsigned int      sp  = s;
    while (j >= result[t - p]) {
      r[s]              = t - s;
      result[t]         = j;
      content[j].count  = content[j].count - 1;

      if (content[j].count == 0) {
        /* detach current element */
        if (ptr->prev) {
          before        = ptr->prev;
          before->next  = ptr->next;
        } else {
          before = NULL;
        }

        if (ptr->next) {
          after       = ptr->next;
          after->prev = ptr->prev;
        } else {
          after = NULL;
        }

        if (!before)
          a = ptr->next;
      }

      if (j != k - 1)
        sp = t + 1;

      if (j == result[t - p])
        sawada_fast(t + 1, p, sp, content, k, r, a, n, results, result_count, result_size);
      else
        sawada_fast(t + 1, t, sp, content, k, r, a, n, results, result_count, result_size);

      if (content[j].count == 0) {
        /* re-attach current element */
        if (before)
          before->next = ptr;
        else
          a = ptr;

        if (after)
          after->prev = ptr;
      }

      content[j].count = content[j].count + 1;

      result = (*results)[*result_count];

      if (ptr->next) {
        ptr = ptr->next;
        j   = ptr->value;
      } else {
        break;
      }
    }
    result[t] = k - 1;
  }
}
