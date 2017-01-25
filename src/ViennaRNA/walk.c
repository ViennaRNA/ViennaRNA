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

PRIVATE bool
isDeletion(move *m)
{
  return m->bpLeft < 0 && m->bpRight < 0;
}


PRIVATE bool
isInsertion(move *m)
{
  return m->bpLeft > 0 && m->bpRight > 0;
}


PRIVATE bool
isShift(move *m)
{
  return (m->bpLeft < 0 && m->bpRight > 0) || (m->bpLeft > 0 && m->bpRight < 0);
}


/**
 * Determines if the neighbored structure is lexicographically smaller if the given move will be performed.
 * „(“ < „)“ < „.“
 * The comparison is easy for insertions and deletions. A deletion is always greater than an insertion on a common structure.
 * Shift moves are more complicated. Actually one has to order all position of each move and compare the first change of a position.
 * Instead of this, it is more general but slower to create the structures as pair tables explicitly and compare only the first position that differs.
 */
PRIVATE bool
isLexicographicallySmaller(short  *ptStructure,
                           move   *m,
                           move   *n)
{
  if (isDeletion(m) && isDeletion(n)) {
    if (MIN2(m->bpLeft, m->bpRight) > MIN2(n->bpLeft, n->bpRight))
      return true;
    else
      return false;
  }

  if (isDeletion(m) && isInsertion(n))
    return true;

  if (isInsertion(m) && isDeletion(n))
    return true;

  if (isInsertion(m) && isInsertion(n)) {
    if (MIN2(m->bpLeft, m->bpRight) < MIN2(n->bpLeft, n->bpRight))
      return true;
    else
      return false;
  }

  /* ...string comparison for shifts and other moves */
  short *s_m  = vrna_ptable_copy(ptStructure);
  short *s_n  = vrna_ptable_copy(ptStructure);
  vrna_perform_move(m, s_m);
  vrna_perform_move(n, s_n);
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


int
vrna_eval_move_shift_pt(vrna_fold_compound_t  *vc,
                        move                  *m,
                        short                 *structure)
{
  if ((m->bpLeft < 0 && m->bpRight > 0) || (m->bpLeft > 0 && m->bpRight < 0)) {
    /* split shift move */
    int   unchangedPosition = m->bpLeft > 0 ? m->bpLeft : m->bpRight;
    int   insertedPosition  = m->bpLeft < 0 ? -m->bpLeft : -m->bpRight;
    int   d1                = -structure[unchangedPosition];
    int   d2                = -unchangedPosition;
    move  deletion;
    if (d1 < d2) {
      deletion.bpLeft   = d2;
      deletion.bpRight  = d1;
    } else {
      deletion.bpLeft   = d1;
      deletion.bpRight  = d2;
    }

    int   i1  = unchangedPosition;
    int   i2  = insertedPosition;
    move  insertion;
    if (i1 > i2) {
      insertion.bpLeft  = i2;
      insertion.bpRight = i1;
    } else {
      insertion.bpLeft  = i1;
      insertion.bpRight = i2;
    }

    int   energy  = vrna_eval_move_pt(vc, structure, deletion.bpLeft, deletion.bpRight);
    short *tmpS   = vrna_ptable_copy(structure);
    vrna_perform_move(&deletion, tmpS);
    energy += vrna_eval_move_pt(vc, tmpS, insertion.bpLeft, insertion.bpRight);
    free(tmpS);
    return energy;
  } else {
    return vrna_eval_move_pt(vc, structure, m->bpLeft, m->bpRight);
  }
}


typedef enum {
  STEEPEST_DESCENT = 1, RANDOM = 2, SHIFTS = 4, OUTPUT_PATH = 8
} options;

PRIVATE move *
vrna_walk(vrna_fold_compound_t  *vc,
          short                 *ptStartAndResultStructure,
          options               opt,
          unsigned int          steps)
{
  int   initialNumberOfMoves  = vc->length;
  move  *moves                = NULL;

  if (opt & OUTPUT_PATH)
    moves = vrna_alloc(sizeof(move) * (initialNumberOfMoves + 1));

  int   numberOfMoves = 0;

  int   energy = vrna_eval_structure_pt(vc, ptStartAndResultStructure);

  int   shifts    = opt & SHIFTS ? true : false;
  move  *moveset  = vrna_build_neighbors(vc, ptStartAndResultStructure, shifts);

  move  *newMoveSet = NULL;
  int   energyNeighbor;
  bool  isDeepest   = false;
  int   iterations  = steps;

  while (((opt & STEEPEST_DESCENT) && !isDeepest) || ((opt & RANDOM) && iterations > 0)) {
    move m = {
      0, 0
    };
    if (opt & STEEPEST_DESCENT) {
      /* determine the deepest neighbor */
      int lowestEnergyIndex = -1;
      int lowestEnergy      = 0;
      int i                 = 0;
      for (move *moveNeighbor = moveset; moveNeighbor->bpLeft != 0; moveNeighbor++, i++) {
        energyNeighbor = vrna_eval_move_shift_pt(vc, moveNeighbor, ptStartAndResultStructure);
        if (energyNeighbor <= lowestEnergy) {
          if ((energyNeighbor == lowestEnergy) && !isLexicographicallySmaller(ptStartAndResultStructure, moveNeighbor, &m))  /* make the walk unique */
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
    } else if (opt & RANDOM) {
      int length = 0;
      for (move *moveNeighbor = moveset; moveNeighbor->bpLeft != 0; moveNeighbor++)
        length++;
      int index = rand() % length;
      m               = moveset[index];
      energyNeighbor  = vrna_eval_move_shift_pt(vc, &m, ptStartAndResultStructure);
      iterations--;
    }

    if (opt & OUTPUT_PATH) {
      if (numberOfMoves > initialNumberOfMoves) {
        initialNumberOfMoves  += vc->length;
        moves                 = vrna_realloc(moves, sizeof(move) * (initialNumberOfMoves + 1));
      }

      moves[numberOfMoves] = m;
      numberOfMoves++;
    }

    int newLength   = 0;
    int movesLength = 0;
    for (move *n = moveset; n->bpLeft != 0; n++)
      movesLength++;
    /* compute neighbors for next round */
    newMoveSet = vrna_build_neighbors_iteratively(vc, &m, ptStartAndResultStructure, (const move *)moveset, movesLength, &newLength, shifts);
    free(moveset);
    moveset = newMoveSet;

    /* adjust pt for next round */
    vrna_perform_move(&m, ptStartAndResultStructure);
    energy += energyNeighbor;

    /* alternative neighbor generation
     * newMoveSet = vrna_build_neighbors(vc, ptStartAndResultStructure, shifts);
     * free(moveset);
     * moveset = newMoveSet;
     */
  }

  if (opt & OUTPUT_PATH) {
    move end = {
      0, 0
    };
    moves[numberOfMoves]  = end;
    moves                 = vrna_realloc(moves, sizeof(move) * (numberOfMoves + 1));
  }

  return moves;
}


void
vrna_walk_gradient(vrna_fold_compound_t *vc,
                   short                *ptStartAndResultStructure,
                   int                  shifts)
{
  options opt = STEEPEST_DESCENT;

  if (shifts != 0)
    opt |= SHIFTS;

  vrna_walk(vc, ptStartAndResultStructure, opt, 0);
}


move *
vrna_walk_gradient_path(vrna_fold_compound_t  *vc,
                        short                 *ptStartAndResultStructure,
                        int                   shifts)
{
  options opt = STEEPEST_DESCENT | OUTPUT_PATH;

  if (shifts != 0)
    opt |= SHIFTS;

  return vrna_walk(vc, ptStartAndResultStructure, opt, 0);
}


void
vrna_walk_random(vrna_fold_compound_t *vc,
                 short                *ptStartAndResultStructure,
                 unsigned int         steps,
                 int                  shifts)
{
  options opt = RANDOM;

  if (shifts != 0)
    opt |= SHIFTS;

  vrna_walk(vc, ptStartAndResultStructure, opt, steps);
}


move *
vrna_walk_random_path(vrna_fold_compound_t  *vc,
                      short                 *ptStartAndResultStructure,
                      unsigned int          steps,
                      int                   shifts)
{
  options opt = RANDOM | OUTPUT_PATH;

  if (shifts != 0)
    opt |= SHIFTS;

  return vrna_walk(vc, ptStartAndResultStructure, opt, steps);
}
