#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "ViennaRNA/utils/basic.h"
#include "ViennaRNA/utils/structures.h"
#include "ViennaRNA/eval.h"
#include "ViennaRNA/landscape/walk.h"

#ifndef bool
#define bool int
#define true 1
#define false 0
#endif

#ifdef __GNUC__
# define INLINE inline
#else
# define INLINE
#endif

#include "local_neighbors.inc"

#define DEBUG   0

/*
 #################################
 # PRIVATE FUNCTION DECLARATIONS #
 #################################
 */
PRIVATE bool
isDeletion(vrna_move_t *m);


PRIVATE bool
isInsertion(vrna_move_t *m);


PRIVATE bool
isShift(vrna_move_t *m);


PRIVATE bool
isLexicographicallySmaller(short        *ptStructure,
                           vrna_move_t  *m,
                           vrna_move_t  *n);


/*
 #################################
 # BEGIN OF FUNCTION DEFINITIONS #
 #################################
 */
PRIVATE bool
isDeletion(vrna_move_t *m)
{
  return m->pos_5 < 0 && m->pos_3 < 0;
}


PRIVATE bool
isInsertion(vrna_move_t *m)
{
  return m->pos_5 > 0 && m->pos_3 > 0;
}


PRIVATE bool
isShift(vrna_move_t *m)
{
  return (m->pos_5 < 0 && m->pos_3 > 0) || (m->pos_5 > 0 && m->pos_3 < 0);
}


/**
 * Determines if the neighbored structure is lexicographically smaller if the given move will be performed.
 * "(" < ")" < "."
 * The comparison is easy for insertions and deletions. A deletion is always greater than an insertion on
 * a common structure. Shift moves are more complicated. Actually one has to order all position of each
 * move and compare the first change of a position. Instead of this, it is more general but slower to
 * create the structures as pair tables explicitly and compare only the first position that differs.
 */
PRIVATE bool
isLexicographicallySmaller(short        *ptStructure,
                           vrna_move_t  *m,
                           vrna_move_t  *n)
{
  if (isDeletion(m) && isDeletion(n)) {
    if (-m->pos_5 > -n->pos_5)
      return true;
    else
      return false;
  }

  if (isDeletion(m) && isInsertion(n))
    return false;

  if (isInsertion(m) && isDeletion(n))
    return true;

  if (isInsertion(m) && isInsertion(n)) {
    if (m->pos_5 < n->pos_5)
      return true;
    else if ((m->pos_5 == n->pos_5) && (m->pos_3 < n->pos_3))
      return true;
    else
      return false;
  }

  /* ...string comparison for shifts and other moves */
  short *s_m  = vrna_ptable_copy(ptStructure);
  short *s_n  = vrna_ptable_copy(ptStructure);
  vrna_move_apply(s_m, m);
  vrna_move_apply(s_n, n);
  bool  isSmaller = false;
  for (int i = 1; i < s_m[0]; i++) {
    if (s_m[i] != s_n[i]) {
      char  c_m = '.';
      char  c_n = '.';
      if (s_m[i] != 0) {
        if (s_m[i] < i)
          c_m = '(';

        if (s_m[i] > i)
          c_m = ')';
      }

      if (s_n[i] != 0) {
        if (s_n[i] < i)
          c_n = '(';
        else
          c_n = ')';
      }

      isSmaller = c_m < c_n;
      break;
    }
  }
  free(s_m);
  free(s_n);
  return isSmaller;
}


PUBLIC vrna_move_t *
vrna_path(vrna_fold_compound_t  *vc,
          short                 *ptStartAndResultStructure,
          unsigned int          steps,
          unsigned int          options)
{
  int         initialNumberOfMoves  = vc->length;
  vrna_move_t *moves                = NULL;

  if (!(options & VRNA_PATH_NO_TRANSITION_OUTPUT))
    moves = vrna_alloc(sizeof(vrna_move_t) * (initialNumberOfMoves + 1));

  int         numberOfMoves = 0;

  int         energy = vrna_eval_structure_pt(vc, ptStartAndResultStructure);

  vrna_move_t *moveset = vrna_neighbors(vc, ptStartAndResultStructure, options);

  vrna_move_t *newMoveSet = NULL;
  int         energyNeighbor;
  bool        isDeepest   = false;
  int         iterations  = steps;

  while (((options & VRNA_PATH_STEEPEST_DESCENT) && !isDeepest) ||
         ((options & VRNA_PATH_RANDOM) && iterations > 0)) {
    vrna_move_t m = {
      0, 0
    };
    if (options & VRNA_PATH_STEEPEST_DESCENT) {
      /* determine the deepest neighbor */
      int lowestEnergyIndex = -1;
      int lowestEnergy      = 0;
      int i                 = 0;
      for (vrna_move_t *moveNeighbor = moveset; moveNeighbor->pos_5 != 0; moveNeighbor++, i++) {
        energyNeighbor = vrna_eval_move_shift_pt(vc, moveNeighbor, ptStartAndResultStructure);
        if (energyNeighbor <= lowestEnergy) {
          /* make the walk unique */
          if ((energyNeighbor == lowestEnergy) &&
              !isLexicographicallySmaller(ptStartAndResultStructure, moveNeighbor, &m))
            continue;

          lowestEnergy      = energyNeighbor;
          lowestEnergyIndex = i;
          m                 = moveset[lowestEnergyIndex];
        }
      }
      if (lowestEnergyIndex == -1) {
        isDeepest = true;
        free(moveset);
        break;
      }

      energyNeighbor  = lowestEnergy;
      m               = moveset[lowestEnergyIndex];
    } else if (options & VRNA_PATH_RANDOM) {
      int length = 0;
      for (vrna_move_t *moveNeighbor = moveset; moveNeighbor->pos_5 != 0; moveNeighbor++)
        length++;
      int index = rand() % length;
      m               = moveset[index];
      energyNeighbor  = vrna_eval_move_shift_pt(vc, &m, ptStartAndResultStructure);
      iterations--;
    }

    if (!(options & VRNA_PATH_NO_TRANSITION_OUTPUT)) {
      if (numberOfMoves > initialNumberOfMoves) {
        initialNumberOfMoves  += vc->length;
        moves                 =
          vrna_realloc(moves, sizeof(vrna_move_t) * (initialNumberOfMoves + 1));
      }

      moves[numberOfMoves] = m;
      numberOfMoves++;
    }

    int newLength   = 0;
    int movesLength = 0;
    for (vrna_move_t *n = moveset; n->pos_5 != 0; n++)
      movesLength++;
    /* compute neighbors for next round */
    newMoveSet =
      vrna_neighbors_successive(vc,
                                &m,
                                ptStartAndResultStructure,
                                (const vrna_move_t *)moveset,
                                movesLength,
                                &newLength,
                                options);
    free(moveset);
    moveset = newMoveSet;

    /* adjust pt for next round */
    vrna_move_apply(ptStartAndResultStructure, &m);
    energy += energyNeighbor;

    /* alternative neighbor generation
     * newMoveSet = vrna_neighbors(vc, ptStartAndResultStructure, options);
     * free(moveset);
     * moveset = newMoveSet;
     */
  }

  if (!(options & VRNA_PATH_NO_TRANSITION_OUTPUT)) {
    vrna_move_t end = {
      0, 0
    };
    moves[numberOfMoves]  = end;
    moves                 = vrna_realloc(moves, sizeof(vrna_move_t) * (numberOfMoves + 1));
  }

  return moves;
}


PUBLIC vrna_move_t *
vrna_path_random(vrna_fold_compound_t *vc,
                 short                *pt,
                 unsigned int         steps,
                 unsigned int         options)
{
  options &= ~VRNA_PATH_STEEPEST_DESCENT;
  options |= VRNA_PATH_RANDOM;

  return vrna_path(vc, pt, steps, options);
}


struct heap {
  short         *pt;
  vrna_move_t   next_move;
  int           *values;
  vrna_move_t   *moves;
  unsigned int  *reverse_idx;
  unsigned int  *reverse_idx_remove;
  unsigned int  num_elements;
  unsigned int  mem_elements;
  void          (*free_cb)(vrna_move_t *);
};


PRIVATE struct heap *
heap_init(unsigned int n,
          void (*free_cb)(vrna_move_t *))
{
  struct heap *h      = (struct heap *)vrna_alloc(sizeof(struct heap));
  unsigned int  size  = (n * (n + 1)) / 2 + 2;

  h->pt                 = NULL;
  h->num_elements       = 0;
  h->mem_elements       = n;
  h->moves              = (vrna_move_t *)vrna_alloc(sizeof(vrna_move_t) * n);
  h->values             = (int *)vrna_alloc(sizeof(int) * n);
  h->reverse_idx        = (unsigned int *)vrna_alloc(sizeof(unsigned int) * size);
  h->reverse_idx_remove = (unsigned int *)vrna_alloc(sizeof(unsigned int) * size);
  h->free_cb            = free_cb;

  return h;
}


PRIVATE void
heap_destroy(struct heap *h)
{
  if (h) {
    free(h->moves);
    free(h->values);
    free(h->reverse_idx);
    free(h->reverse_idx_remove);
    free(h);
  }
}


PRIVATE INLINE unsigned int
heap_parent(unsigned int i)
{
  return floor(i / 2);
}


PRIVATE INLINE unsigned int
heap_left_child(unsigned int i)
{
  return 2 * i;
}


PRIVATE INLINE unsigned int
heap_right_child(unsigned int i)
{
  return 2 * i + 1;
}


PRIVATE INLINE unsigned int
rev_idx(const vrna_move_t *m)
{
  int i, j;

  j = m->pos_5;
  i = m->pos_3;

  if ((i < 0) && (j < 0)) {
    i = -i;
    j = -j;
  }

  /* implement shift moves! */

  return (unsigned int)((i * (i - 1)) / 2 + j);
}


PRIVATE INLINE void
heap_update_rev_idx(struct heap       *h,
                    const vrna_move_t *m,
                    unsigned int      pos)
{
  unsigned int *idx = (vrna_move_is_deletion(m)) ?
                      h->reverse_idx_remove :
                      h->reverse_idx;

  idx[rev_idx(m)] = pos;
}


PRIVATE INLINE unsigned int
heap_find_move(struct heap        *h,
               const vrna_move_t  *m)
{
  unsigned int *idx = (vrna_move_is_deletion(m)) ?
                      h->reverse_idx_remove :
                      h->reverse_idx;

  return idx[rev_idx(m)];
}


PRIVATE INLINE void
heap_swap(struct heap   *h,
          unsigned int  a,
          unsigned int  b)
{
  int         v;
  vrna_move_t m;

  v             = h->values[b];
  m             = h->moves[b];
  h->values[b]  = h->values[a];
  h->moves[b]   = h->moves[a];
  h->values[a]  = v;
  h->moves[a]   = m;

  /* update reverse index */
  heap_update_rev_idx(h, &m, a);
  heap_update_rev_idx(h, &(h->moves[b]), b);
}


PRIVATE int
min_heapify(struct heap   *h,
            unsigned int  i)
{
  int ret = 0;

  while (i > 1) {
    unsigned int  parent  = heap_parent(i);
    int           v       = h->values[parent];

    /*
     * stop heapify-up if heap property is fullfilled, i.e.
     * current node value is larger than that of its parent,
     * or equal but lexigraphically larger
     */
    if ((h->values[i] > v) ||
        ((h->values[i] == v) &&
         (vrna_move_compare(&(h->moves[parent]), &(h->moves[i])) < 0)))
      break;

    heap_swap(h, parent, i);

    i   = parent;
    ret = 1;
  }

  return ret;
}


PRIVATE void
min_heapify_down(struct heap  *h,
                 unsigned int pos)
{
  int           child_v, child_v2, v;
  vrna_move_t   *m, *child_m;
  unsigned int  last_pos, child_pos, child_pos2;

  last_pos = h->num_elements;

  /* nothing to do if already last element */
  if (pos == last_pos)
    return;

  v           = h->values[pos];
  m           = &(h->moves[pos]);
  child_pos   = heap_left_child(pos);
  child_pos2  = heap_right_child(pos);
  child_v     = INF;
  child_m     = NULL;

  /* compare to 1st child */
  if (child_pos <= last_pos) {
    child_v = h->values[child_pos];
    child_m = &(h->moves[child_pos]);
    if ((child_v > v) ||
        ((child_v == v) && (vrna_move_compare(m, child_m) < 0))) {
      child_pos = 0;
      child_v   = v;
      child_m   = m;
    }
  } else {
    child_pos = 0;
    child_v   = v;
    child_m   = m;
  }

  /* compare to 2nd child */
  if (child_pos2 <= last_pos) {
    v = h->values[child_pos2];
    m = &(h->moves[child_pos2]);
    if ((v < child_v) ||
        ((v == child_v) && (vrna_move_compare(m, child_m) < 0))) {
      child_pos = child_pos2;
    }
  }

  if (child_pos) {
    /* swap current node with child */
    heap_swap(h, pos, child_pos);

    min_heapify_down(h, child_pos);
  }
}


PRIVATE void
min_heap_insert(struct heap *h,
                int         value,
                vrna_move_t m)
{
  unsigned int n;

  if ((h) /* && (value <= 0) */) {
    n = ++h->num_elements;

    if (n == h->mem_elements) {
      h->mem_elements *= 1.4;
      h->values       = (int *)vrna_realloc(h->values, sizeof(int) * h->mem_elements);
      h->moves        =
        (vrna_move_t *)vrna_realloc(h->moves, sizeof(vrna_move_t) * h->mem_elements);
    }

    h->values[n]  = value;
    h->moves[n]   = m;

    heap_update_rev_idx(h, &m, n);

    min_heapify(h, n);
  }
}


PRIVATE void
min_heap_remove(struct heap       *h,
                const vrna_move_t m)
{
  if (h) {
    /* get position of last entry in heap */
    unsigned int  last_pos = h->num_elements;

    /* obtain position of element to remove */
    unsigned int  pos = heap_find_move(h, &m);

    if (!pos)
      /* vrna_message_warning("move %d=%d doesn't exist in heap!", m.pos_5, m.pos_3); */
      return;

    /* delete entry for current element */
    heap_update_rev_idx(h, &m, 0);

    h->num_elements--;

    /* we only need to do anything if we didn't remove the last element */
    if (pos != last_pos) {
      h->moves[pos]   = h->moves[last_pos];
      h->values[pos]  = h->values[last_pos];

      /* update reverse index */
      heap_update_rev_idx(h, &(h->moves[pos]), pos);

      if (!min_heapify(h, pos))
        min_heapify_down(h, pos);
    }
  }
}


PRIVATE void
gradient_descent_update_cb(vrna_fold_compound_t *fc,
                           const vrna_move_t    neighbor,
                           unsigned int         state,
                           void                 *data)
{
  struct heap   *h;
  unsigned int  pos;
  int           dG, dG_old;

  h = (struct heap *)data;

  switch (state) {
    case VRNA_NEIGHBOR_REMOVED:
      min_heap_remove(h, neighbor);
      break;

    case VRNA_NEIGHBOR_NEW:
      dG  = vrna_eval_move_pt(fc, h->pt, neighbor.pos_5, neighbor.pos_3);
      min_heap_insert(h, dG, neighbor);
      break;

    case VRNA_NEIGHBOR_CHANGED:
      pos = heap_find_move(h, &neighbor);

      if (!pos) /* insert as new if not already present? */
        return;

      dG      = vrna_eval_move_pt(fc, h->pt, neighbor.pos_5, neighbor.pos_3);
      dG_old  = h->values[pos];

      /* update heap entry value */
      h->values[pos] = dG;

      /* restore min-heap condition */
      if (dG > dG_old)
        min_heapify_down(h, pos);
      else if (dG < dG_old)
        min_heapify(h, pos);

      break;

    default:
      vrna_message_warning("unrecognized state in neighbor callback");
      break;
  }
}


vrna_move_t *
gradient_descent(vrna_fold_compound_t *fc,
                 short                *pt,
                 unsigned int         options)
{
  size_t      num_moves, mem_moves;
  vrna_move_t *moves_applied;

  num_moves     = 0;
  moves_applied = NULL;

  /* obtain initial set of moves to neighboring structures */
  vrna_move_t *neighbors = vrna_neighbors(fc, pt, options);

  /* create initial heap for fast traversal */
  struct heap *h = heap_init(2 * fc->length, NULL);

  for (int i = 0; neighbors[i].pos_5 != 0; i++) {
    int dG = vrna_eval_move_pt(fc, pt, neighbors[i].pos_5, neighbors[i].pos_3);
    min_heap_insert(h, dG, neighbors[i]);
  }

  h->pt = pt;

  if (!(options & VRNA_PATH_NO_TRANSITION_OUTPUT)) {
    mem_moves = 42;
    moves_applied = (vrna_move_t *)vrna_alloc(sizeof(vrna_move_t) * mem_moves);
  }

  /* get current energy */
  while ((h->values[1] <= 0) && (h->num_elements > 0)) {
#if DEBUG
    printf("heap entries so far:\n");
    for (int cnt = 1; cnt <= h->num_elements; cnt++) {
      printf("\t[%d] => %d, %d v=%d\n", cnt, h->moves[cnt].pos_5, h->moves[cnt].pos_3, h->values[cnt]);
    }
#endif

    vrna_move_t next_move = h->moves[1];

    /*
     *  only accept side-ways moves, i.e. dG == 0,
     *  if they result in lexicographically smaller
     *  structure
     */
    if ((h->values[1] == 0) &&
        (vrna_move_is_deletion(&next_move)))
      break;

    vrna_move_neighbor_diff_cb(fc,
                               pt,
                               &(next_move),
                               &gradient_descent_update_cb,
                               (void *)h,
                               options);

    if (moves_applied) {
      moves_applied[num_moves++] = next_move;
      if (num_moves == mem_moves) {
        mem_moves *= 1.4;
        moves_applied = (vrna_move_t *)vrna_realloc(moves_applied, sizeof(vrna_move_t) * mem_moves);
      }
    }
  }

  heap_destroy(h);
  free(neighbors);

  if (moves_applied) {
    moves_applied = (vrna_move_t *)vrna_realloc(moves_applied, sizeof(vrna_move_t) * (num_moves + 1));
    moves_applied[num_moves] = vrna_move_init(0, 0);
  }

  return moves_applied;
}


PUBLIC vrna_move_t *
vrna_path_gradient(vrna_fold_compound_t *vc,
                   short                *pt,
                   unsigned int         options)
{
  options &= ~VRNA_PATH_RANDOM;
  options |= VRNA_PATH_STEEPEST_DESCENT;

#if 0
  return gradient_descent(vc, pt, options);
#else
  return vrna_path(vc, pt, 0, options);
#endif
}
