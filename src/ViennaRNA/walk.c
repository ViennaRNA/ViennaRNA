#include <stdio.h>
#include <stdlib.h>

#include "ViennaRNA/utils.h"
#include "ViennaRNA/structure_utils.h"
#include "ViennaRNA/eval.h"
#include "ViennaRNA/walk.h"

#ifndef bool
#define bool int
#define true 1
#define false 0
#endif

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
    if (MIN2(m->pos_5, m->pos_3) > MIN2(n->pos_5, n->pos_3))
      return true;
    else
      return false;
  }

  if (isDeletion(m) && isInsertion(n))
    return true;

  if (isInsertion(m) && isDeletion(n))
    return true;

  if (isInsertion(m) && isInsertion(n)) {
    if (MIN2(m->pos_5, m->pos_3) < MIN2(n->pos_5, n->pos_3))
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
vrna_path_gradient(vrna_fold_compound_t *vc,
                   short                *pt,
                   unsigned int         options)
{
  options &= ~VRNA_PATH_RANDOM;
  options |= VRNA_PATH_STEEPEST_DESCENT;

  return vrna_path(vc, pt, 0, options);
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
