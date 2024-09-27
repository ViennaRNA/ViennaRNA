#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "ViennaRNA/utils/basic.h"
#include "ViennaRNA/structures/pairtable.h"
#include "ViennaRNA/utils/log.h"
#include "ViennaRNA/datastructures/heap.h"
#include "ViennaRNA/eval/structures.h"
#include <ViennaRNA/landscape/neighbor.h>
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


struct heap_rev_idx {
  vrna_heap_t heap;
  short       *pt;
  size_t      *reverse_idx;
  size_t      *reverse_idx_remove;
};


struct move_en {
  vrna_move_t move;
  int         en;
};

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


PRIVATE vrna_move_t *
do_path(vrna_fold_compound_t  *vc,
        short                 *ptStartAndResultStructure,
        unsigned int          steps,
        unsigned int          options);


PRIVATE vrna_move_t *
gradient_descent(vrna_fold_compound_t *fc,
                 short                *pt,
                 unsigned int         options);


/*
 #################################
 # BEGIN OF FUNCTION DEFINITIONS #
 #################################
 */
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


PUBLIC vrna_move_t *
vrna_path_gradient(vrna_fold_compound_t *vc,
                   short                *pt,
                   unsigned int         options)
{
  options &= ~VRNA_PATH_RANDOM;
  options |= VRNA_PATH_STEEPEST_DESCENT;

  if ((options & VRNA_MOVESET_SHIFT) ||
      (options & VRNA_MOVESET_NO_LP))
    return vrna_path(vc, pt, 0, options);
  else
    return gradient_descent(vc, pt, options);
}


PUBLIC vrna_move_t *
vrna_path(vrna_fold_compound_t  *vc,
          short                 *ptStartAndResultStructure,
          unsigned int          steps,
          unsigned int          options)
{
  if ((vc) && (ptStartAndResultStructure))
    return do_path(vc, ptStartAndResultStructure, steps, options);

  return NULL;
}


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


PRIVATE vrna_move_t *
do_path(vrna_fold_compound_t  *vc,
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
    vrna_move_t m = { 0 };
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
    vrna_move_t end = { 0 };
    moves[numberOfMoves]  = end;
    moves                 = vrna_realloc(moves, sizeof(vrna_move_t) * (numberOfMoves + 1));
  }

  return moves;
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


PRIVATE struct move_en *
move_en_init(vrna_move_t  move,
             int          en)
{
  struct move_en *m = (struct move_en *)vrna_alloc(sizeof(struct move_en));

  m->move = move;
  m->en   = en;

  return m;
}


PRIVATE struct heap_rev_idx *
gradient_descent_data(size_t  n,
                      short   *pt)
{
  size_t              size = (n * (n + 1)) / 2 + 2;

  struct heap_rev_idx *d = (struct heap_rev_idx *)vrna_alloc(sizeof(struct heap_rev_idx));

  d->reverse_idx        = (size_t *)vrna_alloc(sizeof(size_t) * size);
  d->reverse_idx_remove = (size_t *)vrna_alloc(sizeof(size_t) * size);
  d->pt                 = pt;

  return d;
}


PRIVATE void
gradient_descent_data_free(struct heap_rev_idx *d)
{
  free(d->reverse_idx);
  free(d->reverse_idx_remove);
  free(d);
}


PRIVATE void
set_move_pos(const void *m,
             size_t     pos,
             void       *d)
{
  vrna_move_t         *move = &(((struct move_en *)m)->move);

  struct heap_rev_idx *lookup = (struct heap_rev_idx *)d;

  size_t              *idx = (vrna_move_is_removal(move)) ?
                             lookup->reverse_idx_remove :
                             lookup->reverse_idx;

  idx[rev_idx(move)] = pos;
}


PRIVATE size_t
get_move_pos(const void *m,
             void       *d)
{
  vrna_move_t         *move   = &(((struct move_en *)m)->move);
  struct heap_rev_idx *lookup = (struct heap_rev_idx *)d;

  size_t              *idx = (vrna_move_is_removal(move)) ?
                             lookup->reverse_idx_remove :
                             lookup->reverse_idx;

  return idx[rev_idx(move)];
}


PRIVATE int
move_en_compare(const void  *a,
                const void  *b,
                void        *data)
{
  const struct move_en  *m1 = (const struct move_en *)a;
  const struct move_en  *m2 = (const struct move_en *)b;

  if (m1->en < m2->en)
    return -1;
  else if (m1->en > m2->en)
    return 1;
  else
    return vrna_move_compare(&(m1->move),
                             &(m2->move),
                             ((struct heap_rev_idx *)data)->pt);
}


PRIVATE void
gradient_descent_update_cb(vrna_fold_compound_t *fc,
                           const vrna_move_t    neighbor,
                           unsigned int         state,
                           void                 *data)
{
  int                 dG;
  struct move_en      *mm;
  struct heap_rev_idx *lookup;
  vrna_heap_t         h;

  lookup  = (struct heap_rev_idx *)data;
  h       = lookup->heap;

  switch (state) {
    case VRNA_NEIGHBOR_INVALID:
      mm = move_en_init(neighbor, 0);
      free(vrna_heap_remove(h, mm));
      free(mm);

      break;

    case VRNA_NEIGHBOR_NEW:
      dG = vrna_eval_move_pt(fc, lookup->pt, neighbor.pos_5, neighbor.pos_3);
      if (dG <= 0) {
        mm = move_en_init(neighbor, dG);
        vrna_heap_insert(h, mm);
      }

      break;

    case VRNA_NEIGHBOR_CHANGE:
      dG = vrna_eval_move_pt(fc, lookup->pt, neighbor.pos_5, neighbor.pos_3);
      if (dG <= 0) {
        mm = move_en_init(neighbor, dG);
        free(vrna_heap_update(h, mm));
      } else {
        mm = move_en_init(neighbor, 0);
        free(vrna_heap_remove(h, mm));
        free(mm);
      }

      break;

    default:
      vrna_log_warning("unrecognized state in neighbor callback");
      break;
  }
}


PRIVATE vrna_move_t *
gradient_descent(vrna_fold_compound_t *fc,
                 short                *pt,
                 unsigned int         options)
{
  size_t                num_moves, mem_moves, i;
  int                   dG;
  const struct move_en  *next_move_en;
  struct heap_rev_idx   *lookup;
  vrna_heap_t           h;
  vrna_move_t           *neighbors, *moves_applied, next_move;
  void                  *ptr;

  num_moves     = 0;
  moves_applied = NULL;

  /* obtain initial set of moves to neighboring structures */
  neighbors = vrna_neighbors(fc, pt, options);

  /* create initial heap for fast traversal */
  lookup  = gradient_descent_data(fc->length, pt);
  h       = vrna_heap_init(2 * fc->length,
                           &move_en_compare,
                           &get_move_pos,
                           &set_move_pos,
                           (void *)lookup);
  lookup->heap = h;

  for (i = 0; neighbors[i].pos_5 != 0; i++) {
    dG = vrna_eval_move_pt(fc, pt, neighbors[i].pos_5, neighbors[i].pos_3);
    if (dG <= 0) {
      struct move_en *mm = move_en_init(neighbors[i], dG);
      vrna_heap_insert(h, mm);
    }
  }

  if (!(options & VRNA_PATH_NO_TRANSITION_OUTPUT)) {
    mem_moves     = 42;
    moves_applied = (vrna_move_t *)vrna_alloc(sizeof(vrna_move_t) * mem_moves);
  }

  /* get current energy */
  while ((next_move_en = vrna_heap_top(h))) {
    dG        = next_move_en->en;
    next_move = next_move_en->move;

    if ((dG > 0) ||
        ((dG == 0) && vrna_move_is_removal(&(next_move))))
      break;

    vrna_move_neighbor_diff_cb(fc,
                               pt,
                               next_move,
                               &gradient_descent_update_cb,
                               (void *)lookup,
                               options);

    if (moves_applied) {
      moves_applied[num_moves++] = next_move;
      if (num_moves == mem_moves) {
        mem_moves     *= 1.4;
        moves_applied = (vrna_move_t *)vrna_realloc(moves_applied,
                                                    sizeof(vrna_move_t) * mem_moves);
      }
    }
  }

  /* remove remaining entries in heap */
  while ((ptr = vrna_heap_pop(h)))
    free(ptr);

  gradient_descent_data_free(lookup);
  vrna_heap_free(h);
  free(neighbors);

  if (moves_applied) {
    moves_applied = (vrna_move_t *)vrna_realloc(moves_applied,
                                                sizeof(vrna_move_t) * (num_moves + 1));
    moves_applied[num_moves] = vrna_move_init(0, 0);
  }

  return moves_applied;
}
