#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "ViennaRNA/pair_mat.h"
#include "ViennaRNA/file_utils.h"
#include "ViennaRNA/structure_utils.h"
#include "ViennaRNA/model.h"
#include "ViennaRNA/move_set.h"
#include "ViennaRNA/eval.h"

#include "ViennaRNA/neighbor.h"

#ifndef bool
#define bool int
#define true 1
#define false 0
#endif

/**************************************/
/* private helper methods             */
/**************************************/

/**
 * compatible base pair?
 * @param vc - fold compound with sequence
 * @param i - one based index of letter in the sequence string
 * @param j - one based index of letter in the sequence string
 * @return true if a pair is possible
 */
PRIVATE bool
is_compatible(const vrna_fold_compound_t  *vc,
              int                         i,
              int                         j)
{
  return vc->params->model_details.pair[vc->sequence_encoding[i]][vc->sequence_encoding[j]] != 0; /* see pair_mat.h */
}


PRIVATE bool
is_crossing(size_t  i,
            size_t  j,
            size_t  k,
            size_t  l)
{
  if ((i <= k && k <= j && j <= l) || (k <= i && i <= l && l <= j))
    return true;
  else
    return false;
}


/**
 * delete all base pairs in a secondary RNA structure and return a list of single moves.
 */
PRIVATE move *
deletions(vrna_fold_compound_t  *vc,
          const short           *pt,
          int                   *length)
{
  int   len                       = vc->length;
  int   maxStructures             = len / 2;
  move  *listOfDeletionStructures = (move *)malloc(sizeof(move) * (maxStructures + 1));
  int   count                     = 0;

  for (int i = 1; i <= len; i++) {
    if ((pt[i] != 0) & (pt[i] > i)) {
      move *m = &listOfDeletionStructures[count];
      m->bpLeft   = -i;
      m->bpRight  = -pt[i];
      count++;
    }
  }

  *length = count;
  return listOfDeletionStructures;
}


/**
 * insert all possible base pairs in an RNA secondary structure and return a list of single moves.
 */
PRIVATE move *
insertions(vrna_fold_compound_t *vc,
           const short          *pt,
           int                  *length)
{
  int   len = vc->length;
  /* estimate memory for structures. */
  int   maxStructures = (len * len) / 2;

  /* generate structures */
  move  *listOfInsertionStructures  = (move *)malloc(sizeof(move) * (maxStructures + 1));
  int   count                       = 0;
  int   mingap                      = vc->params->model_details.min_loop_size;

  for (int i = 1; i <= len; i++) {
    if (pt[i] == 0) {
      for (int j = i + 1; j <= len; j++) {
        if (pt[j] < i && pt[j] != 0)
          break; /* otherwise it would be crossing. */

        if (pt[j] > j) {
          j = pt[j]; /* skip neighbored base pairs */
          continue;
        }

        if ((j - i) <= mingap)
          continue;

        if ((pt[j] == 0) && is_compatible(vc, i, j)) {
          move *m = &listOfInsertionStructures[count];
          m->bpLeft   = i;
          m->bpRight  = j;
          count++;
        }
      }
    }
  }

  *length = count;
  return listOfInsertionStructures;
}


/**
 * creates all shift moves from position i in the interval [i,end)
 * @param vc - the fold compound with sequence length and parameters
 * @param i - the first position of the move
 * @param start - the start of the interval on the structure; has to be <= i.
 *                if it is not == i, then it is assumed that [start,i] is within the same enclosing loop index!
 * @param end - the end of the interval on the structure
 * @param structure - a secondary structure as pair table
 * @param structures - allocated memory where the move will be stored
 * @param count - the position in the memory where the move will be stored
 */
PRIVATE void
shift_bpins_to_right(const vrna_fold_compound_t *vc,
                     int                        i,
                     int                        start,
                     int                        end,
                     const short                *structure,
                     move                       *structures,
                     int                        *count)
{
  int length  = MIN2(vc->length + 1, end);
  int mingap  = vc->params->model_details.min_loop_size;

  for (int j = start + 1; j < length; j++) {
    /* skip neighbored base pairs */
    while (j < length && (structure[j] > j))
      j = structure[j] + 1;

    if (j >= length)
      break;

    /* stop if other bases pairs would cross the loop at i */
    if (structure[j] < start && structure[j] > 0)
      break;

    /* test if it is a valid pair */
    if (((j - i) > mingap) && is_compatible(vc, i, j)) {
      move *newMove = &structures[(*count)++];
      newMove->bpLeft   = i;
      newMove->bpRight  = -j;
    }
  }
}


/**
 * creates all shift moves from position i in the interval (end,i]
 * @param vc - the fold compound with sequence length and parameters
 * @param i - the first position of the move
 * @param start - the start of the interval on the structure; has to be >= i.
 *                if it is not == i, then it is assumed that [i,start] is within the same enclosing loop index!
 * @param end - the end of the interval on the structure
 * @param structure - a secondary structure as pair table
 * @param structures - allocated memory where the move will be stored
 * @param count - the position in the memory where the move will be stored
 */
PRIVATE void
shift_bpins_to_left(const vrna_fold_compound_t  *vc,
                    int                         i,
                    int                         start,
                    int                         end,
                    const short                 *structure,
                    move                        *structures,
                    int                         *count)
{
  int stop    = MAX2(0, end);
  int mingap  = vc->params->model_details.min_loop_size;

  for (int j = start - 1; j > stop; j--) {
    /* skip neighbored base pairs */
    while (j > stop && (structure[j] < j && structure[j] > 0))
      j = structure[j] - 1;

    if (j <= stop)
      break;

    /* stop if other bases pairs would cross the loop at i */
    if (structure[j] > start)
      break;

    /* test if it is a valid pair */
    if (((i - j) > mingap) && is_compatible(vc, j, i)) {
      move *newMove = &structures[(*count)++];
      newMove->bpLeft   = -j;
      newMove->bpRight  = i;
    }
  }
}



/**
 * Generate all possible shift moves.
 *
 */
PRIVATE move *
shifts(vrna_fold_compound_t *vc,
       const short          *pt,
       int                  *length)
{
  /* Maximal n/2 base pairs per structure times maximal n shift moves --> (n^2)/2 memory */
  int   len               = vc->length;
  int   maxStructures     = (len * len) / 2;
  move  *listOfShiftMoves = (move *)vrna_alloc(sizeof(move) * (maxStructures + 1));

  int   count   = 0;
  int   end = len +1;
  int   pt_i;

  for (int i = 1; i <= len; i++) {
    pt_i = pt[i];
    /* check if we have a pair. */
    if (i < pt_i) {
      /* iterate over all positions with the same loopIndex to find possible pairs. */
      /* 1. first position is fix; search in left direction. */
      shift_bpins_to_left(vc,i,i,0,pt,listOfShiftMoves,&count);
      /* first position is fix; search in right direction (skip the second position) */
      shift_bpins_to_right(vc,i,i,pt_i,pt,listOfShiftMoves,&count);
      shift_bpins_to_right(vc,i,pt_i,end,pt,listOfShiftMoves,&count);
      /* second position is fix */
      shift_bpins_to_left(vc,pt_i,pt_i,i,pt,listOfShiftMoves,&count);
      shift_bpins_to_left(vc,pt_i,i,0,pt,listOfShiftMoves,&count);
      shift_bpins_to_right(vc,pt_i,pt_i,end,pt,listOfShiftMoves,&count);
    }
  }


  *length = count;
  return listOfShiftMoves;
}


/**
 * creates all shift moves to position i from base pairs in the interval [i,end)
 * @param vc - the fold compound with sequence length and parameters
 * @param i - the first position of the move
 * @param start - the start of the interval on the structure; has to be <= i.
 *                if it is not == i, then it is assumed that [start,i] is within the same enclosing loop index!
 * @param end - the end of the interval on the structure
 * @param structure - a secondary structure as pair table
 * @param structures - allocated memory where the move will be stored
 * @param count - the position in the memory where the move will be stored
 */
PRIVATE void
shift_bpins_to_i_from_right(const vrna_fold_compound_t  *vc,
                            int                         i,
                            int                         start,
                            int                         end,
                            const short                 *structure,
                            move                        *structures,
                            int                         *count)
{
  int length  = MIN2(vc->length, end);
  int mingap  = vc->params->model_details.min_loop_size;

  for (int j = start + 1; j < length; j++) {
    /* skip neighbored base pairs */
    while (j < length && (structure[j] > j)) {
      /* test if it is a valid pair */
      if (((j - i) > mingap) && is_compatible(vc, i, j)) {
        move *newMove = &structures[(*count)++];
        newMove->bpLeft   = -i;
        newMove->bpRight  = j;
      }

      j = structure[j];
      if (structure[j] < start && structure[j] > 0)
        break;

      /* test if it is a valid pair */
      if (((j - i) > mingap) && is_compatible(vc, i, j)) {
        move *newMove = &structures[(*count)++];
        newMove->bpLeft   = -i;
        newMove->bpRight  = j;
      }
    }

    if (j > length)
      break;

    /* stop if other bases pairs would cross the loop at i */
    if (structure[j] < start && structure[j] > 0)
      break;
  }
}


/**
 * creates all shift moves to position i from base pairs in the interval (end,i]
 * @param vc - the fold compound with sequence length and parameters
 * @param i - the first position of the move
 * @param start - the start of the interval on the structure; has to be >= i.
 *                if it is not == i, then it is assumed that [i,start] is within the same enclosing loop index!
 * @param end - the end of the interval on the structure
 * @param structure - a secondary structure as pair table
 * @param structures - allocated memory where the move will be stored
 * @param count - the position in the memory where the move will be stored
 */
PRIVATE void
shift_bpins_to_i_from_left(const vrna_fold_compound_t *vc,
                           int                        i,
                           int                        start,
                           int                        end,
                           const short                *structure,
                           move                       *structures,
                           int                        *count)
{
  int stop    = MAX2(0, end);
  int mingap  = vc->params->model_details.min_loop_size;

  for (int j = start - 1; j > stop; j--) {
    /* skip neighbored base pairs */
    while (j > stop && (structure[j] < j && structure[j] > 0)) {
      /* test if it is a valid pair */
      if (((i - j) > mingap) && is_compatible(vc, j, i)) {
        move *newMove = &structures[(*count)++];
        newMove->bpLeft   = j;
        newMove->bpRight  = -i;
      }

      j = structure[j];
      if (structure[j] > start)
        break;

      /* test if it is a valid pair */
      if (((i - j) > mingap) && is_compatible(vc, j, i)) {
        move *newMove = &structures[(*count)++];
        newMove->bpLeft   = j;
        newMove->bpRight  = -i;
      }
    }

    if (j < 1)
      break;

    /* stop if other bases pairs would cross the loop at i */
    if (structure[j] > start)
      break;
  }
}


/**
 * searches for neighbored pairs on the left side next to i. The positions of these pairs are used to generate
 * shift moves with positions in the interval [start, end], which has also to be on the left side next to i.
 */
PRIVATE void
pairs_to_left_most_position_whithin_eclosing_loop_and_shifts_to_interval(const vrna_fold_compound_t *vc,
                                                                         int i,
                                                                         int start,
                                                                         int end,
                                                                         const short *structure,
                                                                         move *structures,
                                                                         int *count,
                                                                         void
                                                                         (*shiftsInInterval)(const vrna_fold_compound_t *, int, int, int, const short *, move *, int *),
                                                                         int includeBorder)
{
  int j = i - 1;

  for (; j > 0; j--) {
    /* skip neighbored base pairs */
    while (j > 0 && (structure[j] < j && structure[j] > 0)) {
      (*shiftsInInterval)(vc, j, start, end, structure, structures, count);
      j = structure[j];
      (*shiftsInInterval)(vc, j, start, end, structure, structures, count);
      continue;
    }
    /* stop if other bases pairs would cross the loop at i */
    if (structure[j] > i) {
      if (includeBorder > 0)
        (*shiftsInInterval)(vc, j, start, end, structure, structures, count);

      break;
    }
  }
}


/**
 * searches for neighbored pairs on the right side next to i. The positions of these pairs are used to generate
 * shift moves with positions in the interval [start, end], which has also to be on the right side next to i.
 */
PRIVATE void
pairs_to_right_most_position_whithin_eclosing_loop_and_shifts_to_interval(const vrna_fold_compound_t *vc,
                                                                          int i,
                                                                          int start,
                                                                          int end,
                                                                          const short *structure,
                                                                          move *structures,
                                                                          int *count,
                                                                          void
                                                                          (*shiftsInInterval)(const vrna_fold_compound_t *, int, int, int, const short *, move *, int *),
                                                                          int includeBorder)
{
  int length  = vc->length;
  int j       = i + 1;

  for (; j <= length; j++) {
    /* skip neighbored base pairs */
    while (j < length && (structure[j] > j)) {
      (*shiftsInInterval)(vc, j, start, end, structure, structures, count);
      j = structure[j];
      if (structure[j] < i && structure[j] > 0)
        break;

      (*shiftsInInterval)(vc, j, start, end, structure, structures, count);
      continue;
    }
    /* stop if other bases pairs would cross the loop at i */
    if (structure[j] < i && structure[j] > 0) {
      if (includeBorder > 0)
        (*shiftsInInterval)(vc, j, start, end, structure, structures, count);

      break;
    }
  }
}


/**
 * searches for neighbored pairs within the given interval. The positions of these pairs are used to generate
 * shift moves with positions in the interval [start, end], which has also to be on the right side next to i.
 */
PRIVATE void
pairs_from_interval_into_shifts_to_interval(const vrna_fold_compound_t *vc,
                                            int i_start,
                                            int i_end,
                                            int start,
                                            int end,
                                            const short *structure,
                                            move *structures,
                                            int *count,
                                            void
                                            (*shiftsInInterval)(const vrna_fold_compound_t *, int, int, int, const short *, move *, int *))
{
  int length  = i_end;
  int j       = i_start + 1;

  for (; j < length; j++) {
    /* skip neighbored base pairs */
    while (j < length && (structure[j] > j)) {
      (*shiftsInInterval)(vc, j, start, end, structure, structures, count);
      j = structure[j];
      (*shiftsInInterval)(vc, j, start, end, structure, structures, count);
      continue;
    }
    /* stop if other bases pairs would cross the loop at i */
    if ((structure[j] < i_start && structure[j] > 0) || (structure[j] > i_end)) {
      printf("there was a crossing shift in a previously freed interval! This is wrong if non-crossing structures are considered.\n");
      break;
    }
  }
}


typedef enum {
  INCREASED, DECREASED, SWITCHED
} intervalType;


PRIVATE move *
generateInsertionsThatWereNotPossibleBeforeThisShiftMove(const vrna_fold_compound_t *vc,
                                                         const short                *structure,
                                                         move                       *freedInterval,
                                                         intervalType               t,
                                                         int                        positivePosition,
                                                         int                        previousPairedPosition,
                                                         int                        newPairedPosition,
                                                         int                        *length)
{
  /* Maximum memory is interval size times structure length for base pairs between these intervals + structure length for bps with the prev. paired position */
  int     intervalLength  = freedInterval->bpRight - freedInterval->bpLeft + 1;
  size_t  resultSize      = intervalLength * vc->length + intervalLength;
  move    *resultList     = vrna_alloc(sizeof(move) * (resultSize + 1));

  /* moves from all pairs inside the freed interval with previous paired position. */
  int     count = 0;

  /* generate inserts between middle and left and right side */
  for (int i = freedInterval->bpLeft; i <= freedInterval->bpRight; i++) {
    /* skip inner neighbors */
    while (structure[i] > i)
      i = structure[i] + 1;
    /* end */
    if (i > freedInterval->bpRight)
      break;

    shift_bpins_to_right(vc, i, freedInterval->bpRight, vc->length + 1, structure, resultList, &count);
    shift_bpins_to_left(vc, i, freedInterval->bpLeft, 0, structure, resultList, &count);
  }

  /* inserts within interval to previousPairedPosition */
  if (previousPairedPosition == freedInterval->bpLeft)
    shift_bpins_to_right(vc, previousPairedPosition, previousPairedPosition - 1, freedInterval->bpRight + 1, structure, resultList, &count);
  else
    shift_bpins_to_left(vc, previousPairedPosition, previousPairedPosition + 1, freedInterval->bpLeft - 1, structure, resultList, &count);

  /* convert shifts to inserts */
  for (int i = 0; i < count; i++) {
    move *m = &resultList[i];
    m->bpLeft   = abs(m->bpLeft);
    m->bpRight  = abs(m->bpRight);
  }

  resultList[count].bpLeft  = 0;
  resultList[count].bpRight = 0;
  resultList                = vrna_realloc(resultList, sizeof(move) * (count + 1));
  *length                   = count;
  return resultList;
}


/**
 * @brief Compute the freed interval after a shift move. These positions can be used to
 *        generate base shifts and insertions that were not possible before a shift move.
 *
 * @structure - the structure as pair table
 * @move - the shift move that will be performed on this structure
 * @freedInterval - output of the left and right position of the freed interval [l,r]
 * @return the type of the freed interval
 */
PRIVATE intervalType
computeFreedInterval(const short  *structure,
                     const move   *m,
                     move         *freedInterval)
{
  int           positivePosition        = MAX2(m->bpLeft, m->bpRight);
  int           newPairedPosition       = abs(MIN2(m->bpLeft, m->bpRight));
  int           previousPairedPosition  = structure[positivePosition];
  intervalType  t;

  /*    |  +)..-) //+ = newPairedPos; - = prevPaired; | = positivePosition (unchanged). */
  if (positivePosition < previousPairedPosition && positivePosition < newPairedPosition) {
    if (newPairedPosition < previousPairedPosition) {
      freedInterval->bpLeft   = newPairedPosition + 1;
      freedInterval->bpRight  = previousPairedPosition;
      t                       = DECREASED;
    } else {
      /*    |  -)..+) */
      freedInterval->bpLeft   = previousPairedPosition;
      freedInterval->bpRight  = newPairedPosition - 1;
      t                       = INCREASED;
    }
  }

  /*  (+  |..-) */
  if (positivePosition < previousPairedPosition && positivePosition > newPairedPosition) {
    freedInterval->bpLeft   = positivePosition + 1;
    freedInterval->bpRight  = previousPairedPosition;
    t                       = SWITCHED;
  }

  /*  (-..|  +) */
  if (positivePosition > previousPairedPosition && positivePosition < newPairedPosition) {
    freedInterval->bpLeft   = previousPairedPosition;
    freedInterval->bpRight  = positivePosition - 1;
    t                       = SWITCHED;
  }

  /*  -(..+(  | */
  if (positivePosition > previousPairedPosition && positivePosition > newPairedPosition) {
    if (newPairedPosition > previousPairedPosition) {
      freedInterval->bpLeft   = previousPairedPosition;
      freedInterval->bpRight  = newPairedPosition - 1;
      t                       = DECREASED;
    } else {
      /*  +(..-(  | */
      freedInterval->bpLeft   = newPairedPosition + 1;
      freedInterval->bpRight  = previousPairedPosition;
      t                       = INCREASED;
    }
  }

  return t;
}


/**
 * @brief Generate all new shift and insertion moves between the freed interval of the shift move and
 *        its environment. Only moves that cross the former structure will be generated.
 *
 * For generating new shift moves that were not possible before, we have to consider 3 interval types
 * (increased, decreased, switched).
 * All in all we have to distinguish 6 cases (all interval types and the mirrored cases).
 *
 *  * Example 1:
 * increase the freed interval with a shift move
 * AAAAGACAAGAAACAAAAGAGAAACAACAAACAAGAAACAAACAAAA
 * ....(....(...)....(.(...)..)......(...)...).... // structure before the shift
 * ....(....(...)....(.(...)......)..(...)...).... // structure after the shift
 * ............................[__]............... // freed interval
 * ..................[________]................... // interval that can pair with the freed interval
 *
 * Example 2:
 * switch the freed interval with a shift move
 * AAAAGACAAGAAACAAAAGAGAAACAACAAACAAGAAACAAACAAAA
 * ....(....(...)....(.(...)..)......(...)...).... // structure before the shift
 * ....(.(..(...)....).(...).........(...)...).... // structure after the shift
 * ...................[_______]................... // freed interval
 * ....[_].....................[_____________].... // intervals that can pair with the freed interval
 *
 * Example 3:
 * decrease the freed interval with a shift move
 * AAAAGACAAGAAACAAAAGAGAAACAACAAACAAGAAACAAACAAAA
 * ....(....(...)....(.(...)......)..(...)...).... // structure before the shift
 * ....(....(...)....(.(...)..)......(...)...).... // structure after the shift
 * ............................[__]............... // freed interval
 * ....[____________].............[__________].... // intervals that can pair with the freed interval
 *
 * @param vc - the fold compound with sequence and parameters
 * @param previousStructure - the structure as pair table
 * @param currentMove - the move that can be applied to the previous structure
 * @param length - outputs the length of the output
 * @result a list with shift and insertion moves that are possible on the current structure and were not possible
 *         on the previous structure
 */
PRIVATE move *
generateShiftsAndInsertionsThatWereNotPossibleBeforeThisShiftMove(const vrna_fold_compound_t  *vc,
                                                                  const short                 *previousStructure,
                                                                  const move                  *currentMove,
                                                                  int                         *length)
{
  short         *currentStructure = vrna_ptable_copy(previousStructure);

  vrna_perform_move(currentMove, currentStructure);
  move          freedInterval = {
    0, 0
  };
  int           positivePosition        = MAX2(currentMove->bpLeft, currentMove->bpRight);
  int           newPairedPosition       = abs(MIN2(currentMove->bpLeft, currentMove->bpRight));
  int           previousPairedPosition  = previousStructure[positivePosition];
  intervalType  t                       = computeFreedInterval(previousStructure, currentMove, &freedInterval);

  move          *allNewShifts = vrna_alloc(sizeof(move) * (vc->length * vc->length));
  int           count         = 0;
  /* moves from all pairs inside the freed interval with previous paired position. */
  if (previousPairedPosition == freedInterval.bpLeft)
    shift_bpins_to_i_from_right(vc, previousPairedPosition, freedInterval.bpLeft - 1, freedInterval.bpRight + 1, currentStructure, allNewShifts,
                                &count);
  else
    shift_bpins_to_i_from_left(vc, previousPairedPosition, freedInterval.bpRight + 1, freedInterval.bpLeft - 1, currentStructure, allNewShifts,
                               &count);

  if (t == INCREASED) {
    /* compute only pairs to one interval instead of two */
    int i_left;
    int i_right;
    if (positivePosition < previousPairedPosition) {
      /* [to pair][freed]   */
      i_left  = positivePosition;
      i_right = previousPairedPosition;
      pairs_to_left_most_position_whithin_eclosing_loop_and_shifts_to_interval(vc, previousPairedPosition, freedInterval.bpLeft - 1,
                                                                               freedInterval.bpRight + 1, currentStructure, allNewShifts, &count,
                                                                               &shift_bpins_to_right, 0);
      /* and now from pairs of the freed interval to the environment */
      pairs_from_interval_into_shifts_to_interval(vc, freedInterval.bpLeft, freedInterval.bpRight, freedInterval.bpLeft, 0, currentStructure,
                                                  allNewShifts, &count, &shift_bpins_to_left);
    } else {
      /* [freed][to pair] */
      i_left  = previousPairedPosition;
      i_right = positivePosition;
      pairs_to_right_most_position_whithin_eclosing_loop_and_shifts_to_interval(vc, previousPairedPosition, freedInterval.bpRight + 1,
                                                                                freedInterval.bpLeft - 1, currentStructure, allNewShifts, &count,
                                                                                &shift_bpins_to_left, 0);

      /* and now from pairs of the freed interval to the environment */
      pairs_from_interval_into_shifts_to_interval(vc, freedInterval.bpLeft - 1, freedInterval.bpRight, freedInterval.bpRight, vc->length + 1,
                                                  currentStructure, allNewShifts, &count, &shift_bpins_to_right);
    }
  } else {
    /* consider both sides [to pair][freed][to pair] */
    /*without current base pair (was possible before) and with exterior pair */
    int startLeftDirection  = 0;
    int startRightDirection = 0;
    if (t == DECREASED) {
      if (positivePosition < newPairedPosition) {
        startLeftDirection  = positivePosition - 1;
        startRightDirection = previousPairedPosition;
      } else {
        startLeftDirection  = previousPairedPosition - 1;
        startRightDirection = positivePosition + 1;
      }
    }

    if (t == SWITCHED) {
      if (newPairedPosition < positivePosition) {
        startLeftDirection  = newPairedPosition - 1;
        startRightDirection = previousPairedPosition + 1;
      } else {
        startLeftDirection  = previousPairedPosition - 1;
        startRightDirection = newPairedPosition + 1;
      }
    }

    pairs_to_left_most_position_whithin_eclosing_loop_and_shifts_to_interval(vc, startLeftDirection + 1, freedInterval.bpLeft - 1,
                                                                             freedInterval.bpRight + 1, currentStructure, allNewShifts, &count,
                                                                             &shift_bpins_to_right, 1);

    pairs_to_right_most_position_whithin_eclosing_loop_and_shifts_to_interval(vc, startRightDirection - 1, freedInterval.bpRight + 1,
                                                                              freedInterval.bpLeft - 1, currentStructure, allNewShifts, &count,
                                                                              &shift_bpins_to_left, 1);
    /* and now from pairs of the freed interval to the environment */
    pairs_from_interval_into_shifts_to_interval(vc, freedInterval.bpLeft - 1, freedInterval.bpRight + 1, freedInterval.bpLeft, 0, currentStructure,
                                                allNewShifts, &count, &shift_bpins_to_left);
    pairs_from_interval_into_shifts_to_interval(vc, freedInterval.bpLeft - 1, freedInterval.bpRight + 1, freedInterval.bpRight, vc->length + 1,
                                                currentStructure, allNewShifts, &count, &shift_bpins_to_right);
  }

  /* and finally all shifts with the new paired position as stable position */
  if (positivePosition < newPairedPosition) {
    /* exclude previous base pair --> two calls */
    shift_bpins_to_left(vc, newPairedPosition, newPairedPosition, positivePosition, currentStructure, allNewShifts, &count);
    shift_bpins_to_left(vc, newPairedPosition, positivePosition, 0, currentStructure, allNewShifts, &count);

    shift_bpins_to_right(vc, newPairedPosition, newPairedPosition, vc->length + 1, currentStructure, allNewShifts, &count);
  } else {
    shift_bpins_to_left(vc, newPairedPosition, newPairedPosition, 0, currentStructure, allNewShifts, &count);

    shift_bpins_to_right(vc, newPairedPosition, newPairedPosition, positivePosition, currentStructure, allNewShifts, &count);
    shift_bpins_to_right(vc, newPairedPosition, positivePosition, vc->length + 1, currentStructure, allNewShifts, &count);
  }

  /* insertion moves */
  int     length_insertionMoves = 0;
  move    *insertionMoves       = generateInsertionsThatWereNotPossibleBeforeThisShiftMove(vc, currentStructure, &freedInterval, t,
                                                                                           positivePosition, previousPairedPosition, newPairedPosition,
                                                                                           &length_insertionMoves);

  /* join shift and insertion lists */
  int     resultSize      = count + length_insertionMoves;
  move    *resultList     = vrna_alloc(sizeof(move) * (resultSize + 1));
  size_t  resultPosition  = 0;
  memcpy(&resultList[resultPosition], allNewShifts, sizeof(move) * count);
  resultPosition += count;
  memcpy(&resultList[resultPosition], insertionMoves, sizeof(move) * length_insertionMoves);
  resultPosition += length_insertionMoves;
  free(allNewShifts);
  free(insertionMoves);
  free(currentStructure);
  /* add terminator */
  resultList[resultPosition].bpLeft   = 0;
  resultList[resultPosition].bpRight  = 0;
  *length                             = resultSize;
  return resultList;
}


/**
 * @brief generate all shifts on the given structure that cross the current move.
 *        The move should already be performed on the structure.
 */
PRIVATE move *
generateCrossingShifts(const vrna_fold_compound_t *vc,
                       const short                *structure,
                       const move                 *currentMove,
                       int                        *length)
{
  int     leftPosition    = MIN2(abs(currentMove->bpRight), abs(currentMove->bpLeft));
  int     rightPosition   = MAX2(abs(currentMove->bpRight), abs(currentMove->bpLeft));
  int     structureLength = vc->length;

  int     count = 0;

  /* Maximum memory */
  size_t  maxCrossingPairs  = 2 * (rightPosition - leftPosition) * (structureLength - (rightPosition - leftPosition));
  size_t  resultSize        = maxCrossingPairs * (structureLength - (rightPosition - leftPosition));

  move    *resultList = vrna_alloc(sizeof(move) * (resultSize + 1));

  /* left to middle */
  pairs_to_left_most_position_whithin_eclosing_loop_and_shifts_to_interval(vc, leftPosition, leftPosition - 1, rightPosition + 1, structure,
                                                                           resultList, &count, &shift_bpins_to_right, 1);

  /* right to middle */
  pairs_to_right_most_position_whithin_eclosing_loop_and_shifts_to_interval(vc, rightPosition, rightPosition + 1, leftPosition - 1, structure,
                                                                            resultList, &count, &shift_bpins_to_left, 1);

  /* middle to left */
  pairs_from_interval_into_shifts_to_interval(vc, leftPosition, rightPosition, leftPosition + 1, 0, structure, resultList, &count,
                                              &shift_bpins_to_left);
  /* middle to right */
  pairs_from_interval_into_shifts_to_interval(vc, leftPosition, rightPosition, rightPosition - 1, structureLength + 1, structure, resultList, &count,
                                              &shift_bpins_to_right);

  resultSize = count;

  resultList[resultSize].bpLeft   = 0;
  resultList[resultSize].bpRight  = 0;
  *length                         = resultSize;
  return resultList;
}


/**
 * @brief generate all insertions on the given structure that cross the current deletion move.
 *        The move should already be performed on the structure.
 */
PRIVATE move *
generateCrossingInserts(const vrna_fold_compound_t  *vc,
                        const short                 *structure,
                        const move                  *currentMove,
                        int                         *length)
{
  int   leftPosition  = MIN2(abs(currentMove->bpRight), abs(currentMove->bpLeft));
  int   rightPosition = MAX2(abs(currentMove->bpRight), abs(currentMove->bpLeft));

  int   len           = vc->length;
  int   allocatedSize = 2 * (rightPosition - leftPosition) * (len - (rightPosition - leftPosition));
  move  *resultList   = vrna_alloc(sizeof(move) * (allocatedSize));
  int   count         = 0;

  int   startLeft = leftPosition + 1;

  /* generate inserts between middle and left and right side */
  for (int i = leftPosition; i <= rightPosition; i++) {
    /* skip inner neighbors */
    while (structure[i] > i)
      i = structure[i] + 1;

    shift_bpins_to_right(vc, i, rightPosition - 1, len + 1, structure, resultList, &count);
    /* ensure that the deleted pair will not be inserted twice */
    if (i == rightPosition)
      startLeft--;

    shift_bpins_to_left(vc, i, startLeft, 0, structure, resultList, &count);
  }

  /* convert shifts to inserts */
  for (int i = 0; i < count; i++) {
    move *m = &resultList[i];
    m->bpLeft   = abs(m->bpLeft);
    m->bpRight  = abs(m->bpRight);
  }

  resultList                = vrna_realloc(resultList, sizeof(move) * (count + 1));
  resultList[count].bpLeft  = 0;
  resultList[count].bpRight = 0;
  *length                   = count;
  return resultList;
}


move *
buildNeighborsForDeletionMove(const vrna_fold_compound_t  *vc,
                              const move                  *currentMove,
                              const short                 *previousStructure,
                              const move                  *previousMoves,
                              int                         length,
                              int                         *newLength,
                              short                       shift)
{
  /* more moves */
  move  *newMoves = (move *)vrna_alloc(sizeof(move) * (length));
  int   newCount  = 0;

  /* insert insertion moves and deletion moves and shift moves, except currentMove and except shifts that have a common position with the deletion. */
  for (int i = 0; i < length; i++) {
    const move *pm = &previousMoves[i];
    if (!(abs(pm->bpLeft) == abs(currentMove->bpLeft) || abs(pm->bpRight) == abs(currentMove->bpRight)
          || abs(pm->bpLeft) == abs(currentMove->bpRight) || abs(pm->bpRight) == abs(currentMove->bpLeft))) {
      newMoves[newCount] = *pm;
      newCount++;
    }
  }
  short *currentStructure = vrna_ptable_copy(previousStructure);
  int   leftPosition      = MIN2(abs(currentMove->bpRight), abs(currentMove->bpLeft));
  int   rightPosition     = MAX2(abs(currentMove->bpRight), abs(currentMove->bpLeft));
  currentStructure[leftPosition]  = 0;
  currentStructure[rightPosition] = 0;

  /* create crossing insert structures in within loopIndex with the deleted bp positions. */
  int   lengthInserts     = 0;
  move  *crossingInserts  = generateCrossingInserts(vc, currentStructure, currentMove, &lengthInserts);

  int   lengthCrossingShifts  = 0;
  move  *crossingShiftMoves   = NULL;

  if (shift)
    /* generate shifts that cross the deleted base pair */
    crossingShiftMoves = generateCrossingShifts(vc, currentStructure, currentMove, &lengthCrossingShifts);

  int totalSize = newCount + lengthCrossingShifts + lengthInserts;
  newMoves = vrna_realloc(newMoves, sizeof(move) * (totalSize + 1));
  memcpy(&newMoves[newCount], crossingInserts, sizeof(move) * (lengthInserts));
  newCount += lengthInserts;
  memcpy(&newMoves[newCount], crossingShiftMoves, sizeof(move) * (lengthCrossingShifts));
  newCount += lengthCrossingShifts;

  if (shift)
    free(crossingShiftMoves);

  free(crossingInserts);
  free(currentStructure);

  *newLength                  = newCount;
  newMoves[newCount].bpLeft   = 0;
  newMoves[newCount].bpRight  = 0;
  return newMoves;
}


move *
buildNeighborsForInsertionMove(const vrna_fold_compound_t *vc,
                               const move                 *currentMove,
                               const short                *previousStructure,
                               const move                 *previousMoves,
                               int                        length,
                               int                        *newLength,
                               short                      shift)
{
  int   allocatedSize = length;
  move  *newMoves     = (move *)vrna_alloc(sizeof(move) * (length));
  int   newCount      = 0;

  /* insert current move as deletion */
  move  nm = {
    -abs(currentMove->bpLeft), -abs(currentMove->bpRight)
  };

  newMoves[newCount++] = nm;

  /* copy filtered previous moves */
  for (int i = 0; i < length; i++) {
    const move *pm = &previousMoves[i];
    if (!is_crossing(abs(pm->bpLeft), abs(pm->bpRight), abs(currentMove->bpLeft), abs(currentMove->bpRight))) {
      /* insert all moves if they are not crossing with the current move. */
      newMoves[newCount++] = *pm;
      continue;
    } else {
      if (shift) {
        /* current is insertion and previous move is crossing and previous move is insertion
         * --> convert to shift if one position is equal */
        if (pm->bpLeft > 0 && pm->bpRight > 0) {
          move  insertionAsShift  = *pm;
          bool  isConvertable     = false;
          if ((pm->bpLeft == currentMove->bpLeft) || (pm->bpLeft == currentMove->bpRight)) {
            insertionAsShift.bpRight  = -pm->bpRight;
            isConvertable             = true;
          }

          if ((pm->bpRight == currentMove->bpLeft) || (pm->bpRight == currentMove->bpRight)) {
            insertionAsShift.bpLeft = -pm->bpLeft;
            isConvertable           = true;
          }

          if ((pm->bpLeft == currentMove->bpLeft && pm->bpRight == currentMove->bpRight)
              || (pm->bpRight == currentMove->bpLeft && pm->bpLeft == currentMove->bpRight))
            continue;

          if (isConvertable) {
            if (newCount >= allocatedSize) {
              allocatedSize += vc->length;
              newMoves      = vrna_realloc(newMoves, sizeof(move) * allocatedSize);
            }

            newMoves[newCount++] = insertionAsShift;
            continue;
          }
        }
      }
    }
  }
  *newLength                  = newCount;
  newMoves                    = (move *)vrna_realloc(newMoves, sizeof(move) * (newCount + 1));
  newMoves[newCount].bpLeft   = 0;
  newMoves[newCount].bpRight  = 0;
  return newMoves;
}


move *
buildNeighborsForShiftMove(const vrna_fold_compound_t *vc,
                           const move                 *currentMove,
                           const short                *previousStructure,
                           const move                 *previousMoves,
                           int                        length,
                           int                        *newLength,
                           short                      shift)
{
  int   allocatedSize = length;
  move  *newMoves     = (move *)vrna_alloc(sizeof(move) * (allocatedSize));
  int   newCount      = 0;

  int   currentPositivePosition = currentMove->bpLeft;

  if (currentMove->bpRight > 0)
    currentPositivePosition = currentMove->bpRight;

  /* copy filtered previous moves */
  int previousPairedPosition = previousStructure[currentPositivePosition];
  for (int i = 0; i < length; i++) {
    const move  *pm     = &previousMoves[i];
    bool        isShift = (pm->bpLeft > 0 && pm->bpRight < 0) || (pm->bpLeft < 0 && pm->bpRight > 0);
    if (!is_crossing(abs(pm->bpLeft), abs(pm->bpRight), abs(currentMove->bpLeft), abs(currentMove->bpRight))) {
      /* if current move is a shift and pm is a shift */
      if (isShift) {
        /*  current negative position cannot be pm.positive position --> allow only if current pos == shift pos. */
        int positivePosition = pm->bpLeft;
        if (pm->bpRight > 0)
          positivePosition = pm->bpRight;

        if (positivePosition != previousPairedPosition) {
          /*  no conflict with this shift --> copy it */
          newMoves[newCount++] = *pm;
          continue;
        }
      } else {
        /*  if non-crossing and previous move is ins or del --> copy */
        newMoves[newCount++] = *pm;
        continue;
      }
    } else {
      /*  if current is shift and prev is crossing and prev is shift */
      if (isShift) {
        /*  if is shift */
        int positivePosition = MAX2(pm->bpLeft, pm->bpRight);
        if ((positivePosition == currentMove->bpLeft || positivePosition == currentMove->bpRight)
            && !((pm->bpLeft == currentMove->bpLeft && pm->bpRight == currentMove->bpRight)
                 || (pm->bpRight == currentMove->bpLeft && pm->bpLeft == currentMove->bpRight))) {
          /*  if one position in common, but not the same move */
          newMoves[newCount++] = *pm;
          continue;
        }
      }
    }
  }

  /* insert current move as deletion */
  move  nm = {
    -abs(currentMove->bpLeft), -abs(currentMove->bpRight)
  };
  newMoves[newCount++] = nm;

  /* insert back shift */
  int   left      = MIN2(currentPositivePosition, previousPairedPosition);
  int   right     = MAX2(currentPositivePosition, previousPairedPosition);
  move  backShift = {
    left, right
  };
  if (previousPairedPosition == backShift.bpLeft)
    backShift.bpLeft = -backShift.bpLeft;
  else
    backShift.bpRight = -backShift.bpRight;

  newMoves[newCount++] = backShift;

  /*  compute all shift moves and insertion moves that end or start in the interval, that was freed by the current shift move. */
  int   length_newFreedIntervalMoves  = 0;
  move  *newFreedIntervalMoves        = generateShiftsAndInsertionsThatWereNotPossibleBeforeThisShiftMove(vc, previousStructure, currentMove,
                                                                                                          &length_newFreedIntervalMoves);

  int   expectedSize = newCount + length_newFreedIntervalMoves;
  if (expectedSize >= allocatedSize) {
    allocatedSize = expectedSize + 1;
    newMoves      = vrna_realloc(newMoves, sizeof(move) * (allocatedSize));
  }

  memcpy(&newMoves[newCount], newFreedIntervalMoves, sizeof(move) * (length_newFreedIntervalMoves));
  newCount += length_newFreedIntervalMoves;
  free(newFreedIntervalMoves);

  *newLength                  = newCount;
  newMoves                    = (move *)vrna_realloc(newMoves, sizeof(move) * (newCount + 1));
  newMoves[newCount].bpLeft   = 0;
  newMoves[newCount].bpRight  = 0;
  return newMoves;
}


/**************************************/
/* public helper methods              */
/**************************************/

void
vrna_perform_move(const move  *m,
                  short       *ptStructure)
{
  /* deletion */
  if (m->bpLeft < 0 && m->bpRight < 0) {
    ptStructure[-m->bpLeft]   = 0;
    ptStructure[-m->bpRight]  = 0;
    return;
  }

  /* insertion */
  if (m->bpLeft > 0 && m->bpRight > 0) {
    ptStructure[m->bpLeft]  = m->bpRight;
    ptStructure[m->bpRight] = m->bpLeft;
    return;
  }

  /* shift right */
  if (m->bpLeft > 0 /*&& m->bpRight < 0*/) {
    short previousPairedPosition = ptStructure[m->bpLeft];
    ptStructure[previousPairedPosition] = 0;
    short newPairedPosition = -m->bpRight;
    ptStructure[m->bpLeft]          = newPairedPosition;
    ptStructure[newPairedPosition]  = m->bpLeft;
    return;
  }

  /* shift left */
  if (m->bpLeft < 0 /*&& m->bpRight > 0*/) {
    short previousPairedPosition = ptStructure[m->bpRight];
    ptStructure[previousPairedPosition] = 0;
    short newPairedPosition = -m->bpLeft;
    ptStructure[m->bpRight]         = newPairedPosition;
    ptStructure[newPairedPosition]  = m->bpRight;
    return;
  }
}


void
vrna_perform_move_str(const move  *m,
                      const short *ptStructure,
                      char        *structure)
{
  /* deletion */
  if (m->bpLeft < 0 && m->bpRight < 0) {
    structure[-m->bpLeft - 1]   = '.';
    structure[-m->bpRight - 1]  = '.';
    return;
  }

  /* insertion */
  if (m->bpLeft > 0 && m->bpRight > 0) {
    structure[m->bpLeft - 1]  = '(';
    structure[m->bpRight - 1] = ')';
    return;
  }

  /* shift right */
  if (m->bpLeft > 0) {
    short previousPairedPosition = ptStructure[m->bpLeft];
    structure[previousPairedPosition - 1] = '.';
    int   left  = m->bpLeft - 1;
    int   right = -m->bpRight - 1;
    structure[left]   = '(';
    structure[right]  = ')';
    return;
  }

  /* shift left */
  if (m->bpLeft < 0) {
    short previousPairedPosition = ptStructure[m->bpRight];
    structure[previousPairedPosition - 1] = '.';
    int   left  = -m->bpLeft - 1;
    int   right = m->bpRight - 1;
    structure[left]   = '(';
    structure[right]  = ')';
    return;
  }
}


void
vrna_update_loop_indices(int          length,
                         int          *loopIndices,
                         const move   *m,
                         const short  *structure)
{
  int index   = loopIndices[abs(m->bpLeft)];
  int bpLeft  = m->bpLeft;
  int bpRight = m->bpRight;

  /* left index has to be smaller for the next iterations */
  if (abs(bpLeft) > abs(bpRight)) {
    bpLeft  = m->bpRight;
    bpRight = m->bpLeft;
  }

  if (bpLeft < 0 && bpRight < 0) {
    /* deletion -> decrease loopIndex. */
    /* search for next closing bracket or unpaired loopIndex to the right. */
    int currentIndex = 0;
    for (int i = -bpRight + 1; i <= length; i++) {
      if (loopIndices[i] < index) {
        /* if unpaired or closing. an loopIndex < current. */
        if (structure[i] <= 0 || structure[i] < i) {
          currentIndex = loopIndices[i];
          break;
        }
      }
    }

    int highestIndexInside = index;
    for (int i = -bpLeft; i <= -bpRight; i++) {
      if (loopIndices[i] > highestIndexInside)
        highestIndexInside = loopIndices[i];

      if (loopIndices[i] == index)
        loopIndices[i] = currentIndex;
      else
        loopIndices[i]--;
    }
    /*  decrease all following indices. */
    for (int i = -bpRight + 1; i <= length; i++)
      if (loopIndices[i] > highestIndexInside)
        loopIndices[i]--;

    loopIndices[0]--;
    return;
  } else if (bpLeft > 0 && bpRight > 0) {
    /* insertion --> 1.Increase the index for all following loopIndices inside the interval. */
    /* 2.Then look at the index on the left side and increase it for the base pair interval (ignore bp in between). */
    /* 3.Increase the indices to the right side of the interval. */
    int currentIndex = 0;
    for (int i = bpLeft; i > 0; i--) {
      if (structure[i] > i) {
        /*  if bp opens */
        currentIndex = loopIndices[i];
        break;
      }
    }
    currentIndex++;
    for (int i = bpLeft; i <= bpRight; i++) {
      if (loopIndices[i] >= currentIndex)
        loopIndices[i]++;
      else
        loopIndices[i] = currentIndex;
    }
    for (int i = bpRight + 1; i <= length; i++)
      if (loopIndices[i] >= currentIndex)
        loopIndices[i]++;

    loopIndices[0]++;
    return;
  } else {
    /* shift move: split it into one deletion and one insertion. */
    move deletionMove;
    if (bpLeft > 0) {
      deletionMove.bpLeft   = -bpLeft;
      deletionMove.bpRight  = -structure[bpLeft];
    } else {
      deletionMove.bpRight  = -bpRight;
      deletionMove.bpLeft   = -structure[bpRight];
    }

    /* left index absolute value has to be smaller than the right one for the other methods */
    if (deletionMove.bpLeft < deletionMove.bpRight) {
      int tmp = deletionMove.bpLeft;
      deletionMove.bpLeft   = deletionMove.bpRight;
      deletionMove.bpRight  = tmp;
    }

    move  insertionMove = {
      abs(bpLeft), abs(bpRight)
    };

    vrna_update_loop_indices(length, loopIndices, &deletionMove, structure);
    short *tmpStructure = vrna_ptable_copy(structure);
    vrna_perform_move(&deletionMove, tmpStructure);
    vrna_update_loop_indices(length, loopIndices, &insertionMove, tmpStructure);
    free(tmpStructure);
    return;
  }
}


/**************************************/
/* public neighbor methods            */
/**************************************/

move *
vrna_build_neighbors(vrna_fold_compound_t *vc,
                     const short          *pt,
                     int                  shift)
{
  int   lengthDeletions;
  move  *deletionList = deletions(vc, pt, &lengthDeletions);

  int   lengthInsertions;
  move  *insertionList = insertions(vc, pt, &lengthInsertions);

  int   totalLength = lengthDeletions + lengthInsertions;
  /* merge deletion list with insertion list. */
  move  *moveSet        = (move *)vrna_realloc(insertionList, sizeof(move) * (totalLength + 1));
  move  *ptInsertionEnd = moveSet + lengthInsertions;

  memcpy(ptInsertionEnd, deletionList, lengthDeletions * sizeof(move));

  if (shift == 1) {
    /*  append shift moves */
    int   lengthShifts;
    move  *shiftList = shifts(vc, pt, &lengthShifts);
    totalLength += lengthShifts;
    moveSet     = (move *)vrna_realloc(moveSet, sizeof(move) * (totalLength + 1));
    move  *ptMoveSetEnd = moveSet + lengthInsertions + lengthDeletions;
    memcpy(ptMoveSetEnd, shiftList, lengthShifts * sizeof(move));
    free(shiftList);
  }

  /* terminate list */
  moveSet[totalLength].bpLeft   = 0;
  moveSet[totalLength].bpRight  = 0;

  free(deletionList);

  return moveSet;
}


move *
vrna_build_neighbors_iteratively(const vrna_fold_compound_t *vc,
                                 const move                 *currentMove,
                                 const short                *previousStructure,
                                 const move                 *previousMoves,
                                 int                        length,
                                 int                        *newLength,
                                 short                      shift)
{
  bool  isDeletion  = currentMove->bpLeft < 0 && currentMove->bpRight < 0 ? true : false;
  bool  isInsertion = currentMove->bpLeft > 0 && currentMove->bpRight > 0 ? true : false;
  bool  isShift     = !isDeletion && !isInsertion;

  move  *newMoves;

  if (isDeletion)
    newMoves = buildNeighborsForDeletionMove(vc, currentMove, previousStructure, previousMoves, length, newLength, shift);

  if (isInsertion)
    newMoves = buildNeighborsForInsertionMove(vc, currentMove, previousStructure, previousMoves, length, newLength, shift);

  if (isShift)
    newMoves = buildNeighborsForShiftMove(vc, currentMove, previousStructure, previousMoves, length, newLength, shift);

  return newMoves;
}
