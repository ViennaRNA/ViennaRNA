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
#include <string.h>

#include "ViennaRNA/utils/basic.h"
#include "ViennaRNA/utils/log.h"
#include "ViennaRNA/search/BoyerMoore.h"
#include "ViennaRNA/combinatorics/basic.h"

#ifdef __GNUC__
# define INLINE inline
#else
# define INLINE
#endif

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
  int           value;
  unsigned int  count;
};


/*
 #################################
 # PRIVATE FUNCTION DECLARATIONS #
 #################################
 */
PRIVATE struct  perm_list *
perm_list_insert(struct perm_list *before,
                 int              value);


PRIVATE struct  perm_list *
perm_list_remove_val(struct perm_list *head,
                     int              value);


PRIVATE struct  perm_list *
perm_list_head(struct perm_list *entry);


PRIVATE void
perm_list_destroy(struct perm_list *entry);


PRIVATE int
cmpfunc(const void  *a,
        const void  *b);


PRIVATE void
sawada_fast_finish_perm(struct necklace_content *content,
                        unsigned int            ***results,
                        unsigned int            *result_count,
                        unsigned int            *result_size,
                        unsigned int            n);


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
            unsigned int            *result_size);


PRIVATE void
n_choose_k(unsigned int *current,
           size_t       start,
           size_t       end,
           size_t       selected,
           size_t       k,
           unsigned int ***output,
           size_t       *output_size,
           size_t       *cnt);


PRIVATE INLINE unsigned int
boustrophedon_at(size_t start,
                 size_t end,
                 size_t pos);


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
    result[result_count][i] = num_entities - 1;

  result[result_count][1] = 0;
  content[0].count        = content[0].count - 1;
  if (content[0].count == 0)
    a = perm_list_remove_val(a, 0);

  /* now we iterate to retrieve full permutation(s) */
  sawada_fast(2, 1, 2, content, num_entities, r, a, n, &result, &result_count, &result_size);

  /* resize results list to actual requirements */
  for (i = result_count; i < result_size; i++)
    free(result[i]);
  result =
    (unsigned int **)vrna_realloc(result, sizeof(unsigned int *) * (result_count + 1));
  result[result_count] = NULL;

  /* cleanup memory */
  free(r);
  free(content);
  perm_list_destroy(a);

  return result;
}


PUBLIC unsigned int
vrna_rotational_symmetry_num(const unsigned int *string,
                             size_t             string_length)
{
  return vrna_rotational_symmetry_pos_num(string,
                                          string_length,
                                          NULL);
}


PUBLIC unsigned int
vrna_rotational_symmetry_pos_num(const unsigned int *string,
                                 size_t             string_length,
                                 unsigned int       **positions)
{
  const unsigned int  *ptr;
  unsigned int        matches, max, shifts_size;
  size_t              *badchars, shift, i;

  if ((!string) || (string_length == 0)) {
    if (positions)
      *positions = NULL;

    return 0;
  }

  /* any string is at least order 1 */
  matches = 1;

  if (positions) {
    shifts_size = 10; /* initial guess for the order of rotational symmetry */
    *positions  = vrna_alloc(sizeof(unsigned int) * shifts_size);

    /* store trivial symmetry */
    (*positions)[matches - 1] = 0;
  }

  /* strings of length 1 are order 1 */
  if (string_length == 1) {
    /* resize positions array to actual length */
    if (positions)
      *positions = vrna_realloc(*positions, sizeof(unsigned int) * matches);

    return matches;
  }

  /* determine largest number/character in string */
  max = string[0];
  for (i = 1; i < string_length; i++)
    max = MAX2(max, string[i]);

  badchars  = vrna_search_BM_BCT_num(string, string_length, max);
  shift     = 1; /* skip trivial symmetry */

  /* detect order of rotational symmetry */

  /*
   *  Note, that finding the smallest shift s of the string that
   *  results in an identity mapping of the string to itself
   *  already determines the order of rotational symmetry R, i.e.
   *  R = n / r where n is the length of the string
   */
  ptr = vrna_search_BMH_num(string,
                            string_length,
                            string,
                            string_length,
                            shift,
                            badchars,
                            1);
  if (ptr) {
    shift   = ptr - string;
    matches = string_length / shift;
    if (positions) {
      *positions = vrna_realloc(*positions, sizeof(unsigned int) * matches);
      for (i = 0; i < matches; i++)
        (*positions)[i] = i * shift;
    }
  }

  free(badchars);

  return matches;
}


PUBLIC unsigned int
vrna_rotational_symmetry(const char *string)
{
  return vrna_rotational_symmetry_pos(string,
                                      NULL);
}


PUBLIC unsigned int
vrna_rotational_symmetry_pos(const char   *string,
                             unsigned int **positions)
{
  const char    *ptr;
  unsigned int  matches, shifts_size;
  size_t        *badchars, shift, i, string_length;

  if (!string) {
    if (positions)
      *positions = NULL;

    return 0;
  }

  string_length = strlen(string);

  if (string_length == 0) {
    if (positions)
      *positions = NULL;

    return 0;
  }

  /* any string is at least order 1 */
  matches = 1;

  if (positions) {
    shifts_size = 10; /* initial guess for the order of rotational symmetry */
    *positions  = vrna_alloc(sizeof(unsigned int) * shifts_size);

    /* store trivial symmetry */
    (*positions)[matches - 1] = 0;
  }

  /* strings of length 1 are order 1 */
  if (string_length == 1) {
    /* resize positions array to actual length */
    if (positions)
      *positions = vrna_realloc(*positions, sizeof(unsigned int) * matches);

    return matches;
  }

  /* determine largest number/character in string */
  badchars = vrna_search_BM_BCT(string);

  shift = 1; /* skip trivial symmetry */

  /* detect order of rotational symmetry */

  /*
   *  Note, that finding the smallest shift s of the string that
   *  results in an identity mapping of the string to itself
   *  already determines the order of rotational symmetry R, i.e.
   *  R = n / r where n is the length of the string
   */
  ptr = vrna_search_BMH(string,
                        string_length,
                        string,
                        string_length,
                        shift,
                        badchars,
                        1);

  if (ptr) {
    shift   = ptr - string;
    matches = string_length / shift;
    if (positions) {
      *positions = vrna_realloc(*positions, sizeof(unsigned int) * matches);
      for (i = 0; i < matches; i++)
        (*positions)[i] = i * shift;
    }
  }

  free(badchars);

  return matches;
}


/**
 * @brief Compute the order of rotational symmetry for a secondary structure @p s
 *
 * This is the size of the stabilizer of @p s, i.e. the set of cyclic permutations
 * of the strand identifiers of the base pairs in @p s that map @p s onto
 * itself.
 */
PUBLIC unsigned int
vrna_rotational_symmetry_db(vrna_fold_compound_t  *fc,
                            const char            *structure)
{
  return vrna_rotational_symmetry_db_pos(fc,
                                         structure,
                                         NULL);
}


PUBLIC unsigned int
vrna_rotational_symmetry_db_pos(vrna_fold_compound_t  *fc,
                                const char            *structure,
                                unsigned int          **positions)
{
  unsigned int n, permutations;

  permutations = 0;

  if (positions)
    *positions = NULL;

  if ((fc) && (structure)) {
    n = strlen(structure);
    if (fc->length != n) {
      vrna_log_warning("vrna_rotational_symmetry_db*: "
                       "Sequence and structure have unequal lengths (%d vs. %d)",
                       fc->length,
                       n);
    } else {
      unsigned int  *shifts, s, r, i, j, ii, jj, string_permutations;
      short         *pt;

      shifts = NULL;

      /* any structure has rotational symmetry of at least order 1, i.e. identity */
      string_permutations = permutations = 1;

      if (positions) {
        *positions = vrna_alloc(sizeof(unsigned int));
        /* store trivial symmetry, i.e. identity */
        (*positions)[0] = 0;
      }

      /* single strands only exhibit rotational symmetry if the string is circular */
      if ((fc->strands == 1) && (fc->params->model_details.circ)) {
        /* compute rotational symmetry for the circular sequence */
        string_permutations = vrna_rotational_symmetry_pos(fc->sequence,
                                                           &shifts);
      } else if (fc->strands > 1) {
        /* determine rotational symmetry of the current strand permutation */
        string_permutations = vrna_rotational_symmetry_pos_num(fc->strand_order,
                                                               fc->strands,
                                                               &shifts);
      }

      /*
       *  There are no rotationally symmetric structures if the strand ordering is
       *  not rotationally symmetric, i.e. string_permutations == 1
       */
      if (string_permutations > 1) {
        /*
         *  For each cyclic permutation of the strand(s), we check if the structure
         *  is also an identify mapping. For that purpose, we simply convert the
         *  structure into a pair table, perform the permutation and check for
         *  identity with the initial structure
         */
        pt  = vrna_ptable(structure);
        s   = 0;

        for (r = 1; r < string_permutations; r++) {
          /*
           *  initial shift is given by smallest shift on sequence level
           *  Here, we distinguish circular RNAs from multi stranded ones
           *  as we've stored the shifts on sequence level and strand level,
           *  respectively.
           */
          if (fc->strands == 1) {
            /*
             *  1. Circular string, i.e. shifts[r] already contains the shifts in
             *  nucleotide positions
             */
            s += shifts[r] - shifts[r - 1];
          } else {
            /*
             * 2. For multiple strands, we have to compute the actual nucleotide shift, since
             * shifts[r] represents a shift by the first shifts[r] strands, i.e. their total
             * length in nucleotides
             */
            for (i = shifts[r - 1]; i < shifts[r]; i++)
              s += fc->nucleotides[fc->strand_order[i]].length;
          }

          /*
           *  Finally, go through the structure and decrease rot_s if the rotationally
           *  symmetric arrangement of the strand(s) still leads to an asymetry of the
           *  structure.
           */
          for (i = 1; i <= n; i++) {
            j   = pt[i]; /* pairing partner of i (1-based), or 0 if unpaired */
            ii  = (i + s);
            if (ii > n)
              ii = (ii % (n + 1)) + 1;  /* position i' in the image of rotation by s */

            jj = pt[ii];                /* pairing partner of i' (1-based), or 0 if unpaired */

            /* i is paired? */
            if (j != 0) {
              j += s;
              if (j > n)
                j = (j % (n + 1)) + 1;  /* pairing partner of i (1-based) in the image of rotation by s */
            }

            /* finally check if j is identical to jj, and stop loop in case of mismatch */
            if (j != jj)
              break;
          }

          /* check whether we had a match of the structure under current permutation */
          if (i == n + 1) {
            /*
             *  now, we know the minimal shift to retrieve an identical structure, so
             *  the order of rotational symmetry is given by multiples of that shift
             */
            permutations = fc->length / s;

            /* store permutation if structure was matches successfully */
            if (positions) {
              *positions = vrna_realloc(*positions, sizeof(unsigned int) * permutations);
              for (i = 0; i < permutations; i++)
                (*positions)[i] = i * s;
            }

            /* we are finished */
            break;
          }
        }
        free(pt);
      }

      free(shifts);
    }
  }

  return permutations;
}


PUBLIC unsigned int **
vrna_n_multichoose_k(size_t n,
                     size_t k)
{
  size_t        result_size = 2;
  unsigned int  **result    = NULL;
  unsigned int  *current    = (unsigned int *)vrna_alloc(sizeof(unsigned int) * k);

  result = (unsigned int **)vrna_alloc(sizeof(unsigned int *) * result_size);

  /* We want to enumerate n multichoose k for total strand number n and
   * interacting strands k. For that purpose, we enumerate n + k - 1 choose k
   * and decrease each index position i by i to obtain n multichoose k
   */
  size_t counter = 0;

  n_choose_k(current, 0, n + k - 2, 0, k, &result, &result_size, &counter);

  for (size_t j = 0; j < counter; j++)
    for (size_t i = 0; i < k; i++)
      result[j][i] -= i;

  /* resize to actual requirements */
  result = (unsigned int **)vrna_realloc(result, sizeof(unsigned int *) * (counter + 1));

  /* add end of list marker */
  result[counter] = NULL;

  free(current);

  return result;
}


PUBLIC unsigned int
vrna_boustrophedon_pos(size_t start,
                       size_t end,
                       size_t pos)
{
  if ((end >= start) &&
      (pos <= end - start + 1))
    return boustrophedon_at(start, end, pos);

  return 0;
}


PUBLIC unsigned int *
vrna_boustrophedon(size_t start,
                   size_t end)
{
  unsigned int *seq, pos;

  seq = NULL;

  if (end >= start) {
    seq     = (unsigned int *)vrna_alloc(sizeof(unsigned int) * (end - start + 2));
    seq[0]  = end - start + 1;

    for (pos = 1; pos <= end - start + 1; pos++)
      seq[pos] = boustrophedon_at(start, end, pos);
  }

  return seq;
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
    (*results)    =
      (unsigned int **)vrna_realloc(*results, sizeof(unsigned int *) * (*result_size));
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

    after = before = NULL;

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


PRIVATE void
n_choose_k(unsigned int *current,
           size_t       start,
           size_t       end,
           size_t       selected,
           size_t       k,
           unsigned int ***output,
           size_t       *output_size,
           size_t       *cnt)
{
  if (selected == k) {
    if (*output_size == *cnt) {
      *output_size  *= 2;
      *output       =
        (unsigned int **)vrna_realloc(*output, sizeof(unsigned int *) * (*output_size));
    }

    (*output)[(*cnt)] = (unsigned int *)vrna_alloc(sizeof(unsigned int) * k);

    for (size_t j = 0; j < k; j++)
      (*output)[(*cnt)][j] = current[j];

    (*cnt)++;
    return;
  }

  for (size_t i = start; i <= end && end - i + 1 >= k - selected; i++) {
    current[selected] = (unsigned int)i;
    n_choose_k(current, i + 1, end, selected + 1, k, output, output_size, cnt);
  }

  return;
}


PRIVATE INLINE unsigned int
boustrophedon_at(size_t start,
                 size_t end,
                 size_t pos)
{
  size_t  count   = pos - 1;
  size_t  advance = (size_t)(count / 2);

  return start +
         (end - start) * (count % 2) +
         advance -
         (2 * (count % 2)) * advance;
}
