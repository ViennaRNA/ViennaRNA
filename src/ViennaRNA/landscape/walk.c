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

#include "local_neighbors.inc"

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
  int           *reverse_idx;
  int           *reverse_idx_remove;
  unsigned int  num_elements;
  unsigned int  mem_elements;
  void(*free_cb)(vrna_move_t *);
};


PRIVATE struct heap *
heap_init(unsigned int n,
          void(*free_cb)(vrna_move_t *))
{
  struct heap *h = (struct heap *)vrna_alloc(sizeof(struct heap));

  h->pt           = NULL;
  h->num_elements = 0;
  h->mem_elements = n;
  h->moves        = (vrna_move_t *)vrna_alloc(sizeof(vrna_move_t) * n);
  h->values       = (int *)vrna_alloc(sizeof(int) * n);
  h->reverse_idx  = (int *)vrna_alloc(sizeof(int) * ((n * (n + 1)) / 2 + 2));
  h->reverse_idx_remove  = (int *)vrna_alloc(sizeof(int) * ((n * (n + 1)) / 2 + 2));
  h->free_cb      = free_cb;

  return h;
}


PRIVATE void
heap_destroy(struct heap *h)
{
  if (h) {
    free(h->moves);
    free(h->values);
    free(h);
  }
}


PRIVATE void
min_heapify(struct heap   *h,
            unsigned int  i)
{
  while (i != 0) {
    unsigned int parent = floor(i / 2);
    int v               = h->values[parent];
    vrna_move_t  m      = h->moves[parent];

    printf("compare to parent %d (%d <= %d ?)\n", parent, h->values[i], v);
    if (h->values[i] < v) {
      if (h->values[i] == v) {
      /* only swap if current move is lexigraphically smaller than v */

      }
      /* swap elements */
      h->values[parent] = h->values[i];
      h->values[i]      = v;

      h->moves[parent]  = h->moves[i];
      h->moves[i]       = m;

      int k, l, p, q;
      k = h->moves[parent].pos_5;
      l = h->moves[parent].pos_3;
      p = h->moves[i].pos_5;
      q = h->moves[i].pos_3;

      if ((k < 0) && (l < 0)) {
        k = -k;
        l = -l;
        h->reverse_idx_remove[(k * (k - 1)) / 2 + l] = parent;
      } else {
        h->reverse_idx[(k * (k - 1)) / 2 + l] = parent;
      }

      if ((p < 0) && (q < 0)) {
        p = -p;
        q = -q;
        h->reverse_idx_remove[(p * (p - 1)) / 2 + q] = i;
      } else {
        h->reverse_idx[(p * (p - 1)) / 2 + q] = i;
      }
    } else {
      break;
    }
    i = parent;
  }
}


PRIVATE void
min_heapify_down( struct heap   *h,
                  unsigned int  pos)
{
  int has_child1 = 0, has_child2 = 0, cnt;
  int child_v, child_v1, child_v2, v = h->values[pos];
  vrna_move_t child_m, child_m1, child_m2, m = h->moves[pos];
  unsigned int child[2];

  child[0] = 2 * pos + 1;
  child[1] = 2 * pos + 2;

  if (child[0] < h->num_elements) {
    has_child1 = 1;
    child_m1 = h->moves[child[0]];
    child_v1 = h->values[child[0]];
  }

  if (child[1] < h->num_elements) {
    has_child2 = 1;
    child_m2 = h->moves[child[1]];
    child_v2 = h->values[child[1]];
  }

  child_m = child_m2;
  child_v = child_v2;
  cnt     = 1;

  if (has_child1) {
    if ((!has_child2) || (child_v1 < child_v2)) {
      child_m = child_m1;
      child_v = child_v1;
      cnt = 0;
    }
  } else if (!has_child2) {
    return;
  }

  /* swap current node with child */
  h->values[pos]  = child_v;
  h->moves[pos]   = child_m;

  h->values[child[cnt]] = v;
  h->moves[child[cnt]]  = m;

  int i        = m.pos_5;
  int j        = m.pos_3;
  int  child_i = child_m.pos_5;
  int  child_j = child_m.pos_3;

  printf("move %d=%d v: %d at pos %d below child %d: %d=%d v: %d at pos %d\n",
         i, j, v, pos,
         cnt, child_i, child_j, h->values[pos], child[cnt]);

  if ((i < 0) && (j < 0)) {
    i = -i;
    j = -j;
    h->reverse_idx_remove[(i * (i - 1)) / 2 + j] = child[cnt];
  } else {
    h->reverse_idx[(i * (i - 1)) / 2 + j] = child[cnt];
  }

  if ((child_i < 0) && (child_j < 0)) {
    child_i = -child_i;
    child_j = -child_j;
    h->reverse_idx_remove[(child_i + (child_i - 1)) / 2 + child_j] = pos;
  } else {
    h->reverse_idx[(child_i + (child_i - 1)) / 2 + child_j] = pos;
  }

  min_heapify_down(h, child[cnt]);
}


PRIVATE void
min_heap_insert(struct heap *h,
                int         value,
                vrna_move_t m)
{
  if (h) {
    if (h->num_elements == h->mem_elements) {
      h->mem_elements *= 1.4;
      h->values = (int *)vrna_realloc(h->values, sizeof(int) * h->mem_elements);
      h->moves  = (vrna_move_t *)vrna_realloc(h->moves, sizeof(vrna_move_t) * h->mem_elements);
    }

    printf("inserting %d=%d with value %d\n", m.pos_5, m.pos_3, value);
    h->values[h->num_elements]  = value;
    h->moves[h->num_elements]   = m;
    int i = m.pos_5;
    int j = m.pos_3;
    if ((i < 0) && (j < 0)) {
      i = -i;
      j = -j;
      h->reverse_idx_remove[(i * (i - 1)) / 2 + j] = h->num_elements;
    } else {
      h->reverse_idx[(i * (i - 1)) / 2 + j] = h->num_elements;
    }

    h->num_elements++;

    min_heapify(h, h->num_elements - 1);
  }
}



PRIVATE void
min_heap_remove(struct heap *h,
                const vrna_move_t m)
{
  if (h) {
    int i = m.pos_5;
    int j = m.pos_3;
    int pos = -1;

    if ((i < 0) && (j < 0)) {
      i = -i;
      j = -j;
      pos = h->reverse_idx_remove[(i * (i - 1)) / 2 + j];
      h->reverse_idx_remove[(i * (i - 1)) / 2 + j] = -1;
    } else {
      pos = h->reverse_idx[(i * (i - 1)) / 2 + j];
      h->reverse_idx[(i * (i - 1)) / 2 + j] = -1;
    }

    h->num_elements--;

    /* obtain position of element to remove */

    if (pos < 0) {
      vrna_message_warning("move %d=%d doesn't exist in heap!", m.pos_5, m.pos_3);
      return;
    }

    printf("removing from pos %d [%d]\n", pos, h->num_elements);

    /* get data from last entry in heap */
    unsigned int  last_pos    = h->num_elements;
    int           last_value  = h->values[last_pos];
    int           last_i      = h->moves[last_pos].pos_5;
    int           last_j      = h->moves[last_pos].pos_3;

    /* finally remove the move */
    if (pos != h->num_elements) {
      /* we only need to do anything if we didn't remove the last element */
      h->moves[pos] = h->moves[h->num_elements];
      h->values[pos] = h->values[h->num_elements];
      if ((last_i < 0) && (last_j < 0)) {
        last_i = -last_i;
        last_j = -last_j;
        h->reverse_idx_remove[(last_i * (last_i - 1)) / 2 + last_j] = pos;
      } else {
        h->reverse_idx[(last_i * (last_i - 1)) / 2 + last_j] = pos;
      }
      min_heapify_down(h, pos);
    }
  }
}

void
gradient_descent_update_cb(vrna_fold_compound_t *fc,
                           const vrna_move_t          neighbor,
                           unsigned int               state,
                           void                       *data)
{
  struct heap *h = (struct heap *)data;
  int i, j, pos, dG, dG_old;

  printf("callback move %d=%d\n", neighbor.pos_5, neighbor.pos_3);

  switch (state) {
    case VRNA_NEIGHBOR_REMOVED:
      printf("removing move %d=%d\n", neighbor.pos_5, neighbor.pos_3);
      min_heap_remove(h, neighbor);
      break;

    case VRNA_NEIGHBOR_NEW:
      i   = neighbor.pos_5;
      j   = neighbor.pos_3;

      dG = vrna_eval_move_pt(fc, h->pt, neighbor.pos_5, neighbor.pos_3);
      printf("insert move %d=%d value: %d\n", neighbor.pos_5, neighbor.pos_3, dG);
      min_heap_insert(h, dG, neighbor);
      break;

    case VRNA_NEIGHBOR_CHANGED:
      i   = neighbor.pos_5;
      j   = neighbor.pos_3;

      if ((i < 0) && (j < 0)) {
        i = -i;
        j = -j;
        pos = h->reverse_idx_remove[(i * (i - 1)) / 2 + j];
      } else {
        pos = h->reverse_idx[(i * (i - 1)) / 2 + j];
      }

      /* overwrite in case sign changed */
      h->moves[pos] = neighbor;

      dG = vrna_eval_move_pt(fc, h->pt, neighbor.pos_5, neighbor.pos_3);
      dG_old = h->values[pos];

      /* update heap entry */
      h->values[pos] = dG;
      printf("changed move %d=%d at pos %d from %d to %d\n", neighbor.pos_5, neighbor.pos_3, pos, dG_old, dG);
      /* restore min-heap condition */
      if (dG > dG_old)
        min_heapify_down(h, pos);
      else
        min_heapify(h, pos);

      break;

    default:
      vrna_message_warning("unrecognized state in neighbor callback");
      break;
  }
}


vrna_move_t *
gradient_descent( vrna_fold_compound_t  *fc,
                  short                 *pt,
                  unsigned int          options)
{
  vrna_move_t *moves_applied = NULL;

  /* obtain initial set of moves to neighboring structures */
  vrna_move_t *neighbors = vrna_neighbors(fc, pt, options);

  /* create initial heap for fast traversal */
  struct heap *h  = heap_init(2 * fc->length, NULL);

  for (int i = 0; neighbors[i].pos_5 != 0; i++) {
    int dG = vrna_eval_move_pt(fc, pt, neighbors[i].pos_5, neighbors[i].pos_3);
    min_heap_insert(h, dG, neighbors[i]);
  }

  h->pt = pt;
  printf("\n\n===============\nstart with move %d=%d, value: %d\n", h->moves[0].pos_5, h->moves[0].pos_3, h->values[0]);

  int moves=0;
  /* get current energy */
  while(h->values[0] < 0) {
    vrna_move_t next_move = h->moves[0];
    printf("\n===============\nnext step: %d=%d value: %d\n", next_move.pos_5, next_move.pos_3, h->values[0]);
    vrna_move_neighbor_diff_cb(fc,
                               pt,
                               &(next_move),
                               &gradient_descent_update_cb,
                               (void *)h,
                               options);

    printf("next value: %d\n", h->values[0]);
    moves++;
    if (moves == 40)
      break;
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


