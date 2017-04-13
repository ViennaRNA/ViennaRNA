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

typedef enum {
  INCREASED, DECREASED, SWITCHED
} intervalType;

typedef void (shiftsInInterval)(const vrna_fold_compound_t *,
                                int,
                                int,
                                int,
                                const short *,
                                vrna_move_t *,
                                int *);

PRIVATE bool
is_compatible(const vrna_fold_compound_t  *vc,
              int                         i,
              int                         j);


PRIVATE bool
is_crossing(size_t  i,
            size_t  j,
            size_t  k,
            size_t  l);


PRIVATE vrna_move_t *
deletions(vrna_fold_compound_t  *vc,
          const short           *pt,
          int                   *length);


PRIVATE vrna_move_t *
insertions(vrna_fold_compound_t *vc,
           const short          *pt,
           int                  *length);


PRIVATE void
shift_bpins_to_right(const vrna_fold_compound_t *vc,
                     int                        i,
                     int                        start,
                     int                        end,
                     const short                *structure,
                     vrna_move_t                *structures,
                     int                        *count);


PRIVATE void
shift_bpins_to_left(const vrna_fold_compound_t  *vc,
                    int                         i,
                    int                         start,
                    int                         end,
                    const short                 *structure,
                    vrna_move_t                 *structures,
                    int                         *count);


PRIVATE vrna_move_t *
shifts(vrna_fold_compound_t *vc,
       const short          *pt,
       int                  *length);


PRIVATE void
shift_bpins_to_i_from_right(const vrna_fold_compound_t  *vc,
                            int                         i,
                            int                         start,
                            int                         end,
                            const short                 *structure,
                            vrna_move_t                 *structures,
                            int                         *count);


PRIVATE void
shift_bpins_to_i_from_left(const vrna_fold_compound_t *vc,
                           int                        i,
                           int                        start,
                           int                        end,
                           const short                *structure,
                           vrna_move_t                *structures,
                           int                        *count);


PRIVATE void
pairs_to_left_most_position_whithin_eclosing_loop_and_shifts_to_interval(
  const vrna_fold_compound_t  *vc,
  int                         i,
  int                         start,
  int                         end,
  const short                 *structure,
  vrna_move_t                 *structures,
  int                         *count,
  shiftsInInterval            *interval_func,
  int                         includeBorder);


PRIVATE void
pairs_to_right_most_position_whithin_eclosing_loop_and_shifts_to_interval(
  const vrna_fold_compound_t  *vc,
  int                         i,
  int                         start,
  int                         end,
  const short                 *structure,
  vrna_move_t                 *structures,
  int                         *count,
  shiftsInInterval            *interval_func,
  int                         includeBorder);


PRIVATE void
pairs_from_interval_into_shifts_to_interval(const vrna_fold_compound_t  *vc,
                                            int                         i_start,
                                            int                         i_end,
                                            int                         start,
                                            int                         end,
                                            const short                 *structure,
                                            vrna_move_t                 *structures,
                                            int                         *count,
                                            shiftsInInterval            *interval_func);


PRIVATE vrna_move_t *
generateInsertionsThatWereNotPossibleBeforeThisShiftMove(const vrna_fold_compound_t *vc,
                                                         const short                *structure,
                                                         vrna_move_t                *freedInterval,
                                                         intervalType               t,
                                                         int                        positivePosition,
                                                         int                        previousPairedPosition,
                                                         int                        newPairedPosition,
                                                         int                        *length);


PRIVATE intervalType
computeFreedInterval(const short        *structure,
                     const vrna_move_t  *m,
                     vrna_move_t        *freedInterval);


PRIVATE vrna_move_t *
generateShiftsThatWereNotPossibleBeforeThisShiftMove(const vrna_fold_compound_t  *vc,
                                                                  const short                 *prev_pt,
                                                                  const vrna_move_t           *curr_move,
                                                                  int                         *length);


PRIVATE vrna_move_t *
generateCrossingShifts(const vrna_fold_compound_t *vc,
                       const short                *structure,
                       const vrna_move_t          *curr_move,
                       int                        *length);


PRIVATE vrna_move_t *
generateCrossingInserts(const vrna_fold_compound_t  *vc,
                        const short                 *structure,
                        const vrna_move_t           *curr_move,
                        int                         *length);


PRIVATE vrna_move_t *
buildNeighborsForDeletionMove(const vrna_fold_compound_t  *vc,
                              const vrna_move_t           *curr_move,
                              const short                 *prev_pt,
                              const vrna_move_t           *prev_neighbors,
                              int                         length,
                              int                         *size_neighbors,
                              unsigned int                options);


PRIVATE vrna_move_t *
buildNeighborsForInsertionMove(const vrna_fold_compound_t *vc,
                               const vrna_move_t          *curr_move,
                               const short                *prev_pt,
                               const vrna_move_t          *prev_neighbors,
                               int                        length,
                               int                        *size_neighbors,
                               unsigned int               options);


PRIVATE vrna_move_t *
buildNeighborsForShiftMove(const vrna_fold_compound_t *vc,
                           const vrna_move_t          *curr_move,
                           const short                *prev_pt,
                           const vrna_move_t          *prev_neighbors,
                           int                        length,
                           int                        *size_neighbors,
                           unsigned int               options);


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
PRIVATE vrna_move_t *
deletions(vrna_fold_compound_t  *vc,
          const short           *pt,
          int                   *length)
{
  int         len                       = vc->length;
  int         maxStructures             = len / 2;
  vrna_move_t *listOfDeletionStructures =
    (vrna_move_t *)malloc(sizeof(vrna_move_t) * (maxStructures + 1));
  int         count = 0;

  for (int i = 1; i <= len; i++) {
    if ((pt[i] != 0) & (pt[i] > i)) {
      vrna_move_t *m = &listOfDeletionStructures[count];
      m->pos_5  = -i;
      m->pos_3  = -pt[i];
      count++;
    }
  }

  *length = count;
  return listOfDeletionStructures;
}


/**
 * insert all possible base pairs in an RNA secondary structure and return a list of single moves.
 */
PRIVATE vrna_move_t *
insertions(vrna_fold_compound_t *vc,
           const short          *pt,
           int                  *length)
{
  int         len = vc->length;
  /* estimate memory for structures. */
  int         maxStructures = (len * len) / 2;

  /* generate structures */
  vrna_move_t *listOfInsertionStructures =
    (vrna_move_t *)malloc(sizeof(vrna_move_t) * (maxStructures + 1));
  int         count   = 0;
  int         mingap  = vc->params->model_details.min_loop_size;

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
          vrna_move_t *m = &listOfInsertionStructures[count];
          m->pos_5  = i;
          m->pos_3  = j;
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
                     vrna_move_t                *structures,
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
      vrna_move_t *newMove = &structures[(*count)++];
      newMove->pos_5  = i;
      newMove->pos_3  = -j;
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
                    vrna_move_t                 *structures,
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
      vrna_move_t *newMove = &structures[(*count)++];
      newMove->pos_5  = -j;
      newMove->pos_3  = i;
    }
  }
}


/**
 * Generate all possible shift moves.
 *
 */
PRIVATE vrna_move_t *
shifts(vrna_fold_compound_t *vc,
       const short          *pt,
       int                  *length)
{
  /* Maximal n/2 base pairs per structure times maximal n shift moves --> (n^2)/2 memory */
  int         len               = vc->length;
  int         maxStructures     = (len * len) / 2;
  vrna_move_t *listOfShiftMoves = (vrna_move_t *)vrna_alloc(
    sizeof(vrna_move_t) * (maxStructures + 1));

  int         count = 0;
  int         end   = len + 1;
  int         pt_i;

  for (int i = 1; i <= len; i++) {
    pt_i = pt[i];
    /* check if we have a pair. */
    if (i < pt_i) {
      /* iterate over all positions with the same loopIndex to find possible pairs. */
      /* 1. first position is fix; search in left direction. */
      shift_bpins_to_left(vc, i, i, 0, pt, listOfShiftMoves, &count);
      /* first position is fix; search in right direction (skip the second position) */
      shift_bpins_to_right(vc, i, i, pt_i, pt, listOfShiftMoves, &count);
      shift_bpins_to_right(vc, i, pt_i, end, pt, listOfShiftMoves, &count);
      /* second position is fix */
      shift_bpins_to_left(vc, pt_i, pt_i, i, pt, listOfShiftMoves, &count);
      shift_bpins_to_left(vc, pt_i, i, 0, pt, listOfShiftMoves, &count);
      shift_bpins_to_right(vc, pt_i, pt_i, end, pt, listOfShiftMoves, &count);
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
                            vrna_move_t                 *structures,
                            int                         *count)
{
  int length  = MIN2(vc->length, end);
  int mingap  = vc->params->model_details.min_loop_size;

  for (int j = start + 1; j < length; j++) {
    /* skip neighbored base pairs */
    while (j < length && (structure[j] > j)) {
      /* test if it is a valid pair */
      if (((j - i) > mingap) && is_compatible(vc, i, j)) {
        vrna_move_t *newMove = &structures[(*count)++];
        newMove->pos_5  = -i;
        newMove->pos_3  = j;
      }

      j = structure[j];
      if (structure[j] < start && structure[j] > 0)
        break;

      /* test if it is a valid pair */
      if (((j - i) > mingap) && is_compatible(vc, i, j)) {
        vrna_move_t *newMove = &structures[(*count)++];
        newMove->pos_5  = -i;
        newMove->pos_3  = j;
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
                           vrna_move_t                *structures,
                           int                        *count)
{
  int stop    = MAX2(0, end);
  int mingap  = vc->params->model_details.min_loop_size;

  for (int j = start - 1; j > stop; j--) {
    /* skip neighbored base pairs */
    while (j > stop && (structure[j] < j && structure[j] > 0)) {
      /* test if it is a valid pair */
      if (((i - j) > mingap) && is_compatible(vc, j, i)) {
        vrna_move_t *newMove = &structures[(*count)++];
        newMove->pos_5  = j;
        newMove->pos_3  = -i;
      }

      j = structure[j];
      if (structure[j] > start)
        break;

      /* test if it is a valid pair */
      if (((i - j) > mingap) && is_compatible(vc, j, i)) {
        vrna_move_t *newMove = &structures[(*count)++];
        newMove->pos_5  = j;
        newMove->pos_3  = -i;
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
pairs_to_left_most_position_whithin_eclosing_loop_and_shifts_to_interval(
  const vrna_fold_compound_t  *vc,
  int                         i,
  int                         start,
  int                         end,
  const short                 *structure,
  vrna_move_t                 *structures,
  int                         *count,
  shiftsInInterval            *interval_func,
  int                         includeBorder)
{
  int j = i - 1;

  for (; j > 0; j--) {
    /* skip neighbored base pairs */
    while (j > 0 && (structure[j] < j && structure[j] > 0)) {
      interval_func(vc, j, start, end, structure, structures, count);
      j = structure[j];
      interval_func(vc, j, start, end, structure, structures, count);
      continue;
    }
    /* stop if other bases pairs would cross the loop at i */
    if (structure[j] > i) {
      if (includeBorder > 0)
        interval_func(vc, j, start, end, structure, structures, count);

      break;
    }
  }
}


/**
 * searches for neighbored pairs on the right side next to i. The positions of these pairs are used to generate
 * shift moves with positions in the interval [start, end], which has also to be on the right side next to i.
 */
PRIVATE void
pairs_to_right_most_position_whithin_eclosing_loop_and_shifts_to_interval(
  const vrna_fold_compound_t  *vc,
  int                         i,
  int                         start,
  int                         end,
  const short                 *structure,
  vrna_move_t                 *structures,
  int                         *count,
  shiftsInInterval            *interval_func,
  int                         includeBorder)
{
  int length  = vc->length;
  int j       = i + 1;

  for (; j <= length; j++) {
    /* skip neighbored base pairs */
    while (j < length && (structure[j] > j)) {
      interval_func(vc, j, start, end, structure, structures, count);
      j = structure[j];
      if (structure[j] < i && structure[j] > 0)
        break;

      interval_func(vc, j, start, end, structure, structures, count);
      continue;
    }
    /* stop if other bases pairs would cross the loop at i */
    if (structure[j] < i && structure[j] > 0) {
      if (includeBorder > 0)
        interval_func(vc, j, start, end, structure, structures, count);

      break;
    }
  }
}


/**
 * searches for neighbored pairs within the given interval. The positions of these pairs are used to generate
 * shift moves with positions in the interval [start, end], which has also to be on the right side next to i.
 */
PRIVATE void
pairs_from_interval_into_shifts_to_interval(const vrna_fold_compound_t  *vc,
                                            int                         i_start,
                                            int                         i_end,
                                            int                         start,
                                            int                         end,
                                            const short                 *structure,
                                            vrna_move_t                 *structures,
                                            int                         *count,
                                            shiftsInInterval            *interval_func)
{
  int length  = i_end;
  int j       = i_start + 1;

  for (; j < length; j++) {
    /* skip neighbored base pairs */
    while (j < length && (structure[j] > j)) {
      interval_func(vc, j, start, end, structure, structures, count);
      j = structure[j];
      interval_func(vc, j, start, end, structure, structures, count);
      continue;
    }
    /* stop if other bases pairs would cross the loop at i */
    if ((structure[j] < i_start && structure[j] > 0) || (structure[j] > i_end)) {
      vrna_message_warning(
        "there was a crossing shift in a previously freed interval! This is wrong if non-crossing structures are considered.\n");
      break;
    }
  }
}


PRIVATE vrna_move_t *
generateInsertionsThatWereNotPossibleBeforeThisShiftMove(const vrna_fold_compound_t *vc,
                                                         const short                *structure,
                                                         vrna_move_t                *freedInterval,
                                                         intervalType               t,
                                                         int                        positivePosition,
                                                         int                        previousPairedPosition,
                                                         int                        newPairedPosition,
                                                         int                        *length)
{
  /* Maximum memory is interval size times structure length for base pairs between these intervals + structure length for bps with the prev. paired position */
  int         intervalLength  = freedInterval->pos_3 - freedInterval->pos_5 + 1;
  size_t      resultSize      = intervalLength * vc->length + intervalLength;
  vrna_move_t *resultList     = vrna_alloc(sizeof(vrna_move_t) * (resultSize + 1));

  /* moves from all pairs inside the freed interval with previous paired position. */
  int         count = 0;

  /* generate inserts between middle and left and right side */
  for (int i = freedInterval->pos_5; i <= freedInterval->pos_3; i++) {
    /* skip inner neighbors */
    while (structure[i] > i)
      i = structure[i] + 1;
    /* end */
    if (i > freedInterval->pos_3)
      break;

    shift_bpins_to_right(vc,
                         i,
                         freedInterval->pos_3,
                         vc->length + 1,
                         structure,
                         resultList,
                         &count);
    shift_bpins_to_left(vc, i, freedInterval->pos_5, 0, structure, resultList, &count);
  }

  /* inserts within interval to previousPairedPosition */
  if (previousPairedPosition == freedInterval->pos_5) {
    shift_bpins_to_right(vc,
                         previousPairedPosition,
                         previousPairedPosition - 1,
                         freedInterval->pos_3 + 1,
                         structure,
                         resultList,
                         &count);
  } else {
    shift_bpins_to_left(vc,
                        previousPairedPosition,
                        previousPairedPosition + 1,
                        freedInterval->pos_5 - 1,
                        structure,
                        resultList,
                        &count);
  }

  /* convert shifts to inserts */
  for (int i = 0; i < count; i++) {
    vrna_move_t *m = &resultList[i];
    m->pos_5  = abs(m->pos_5);
    m->pos_3  = abs(m->pos_3);
  }

  resultList[count].pos_5 = 0;
  resultList[count].pos_3 = 0;
  resultList              = vrna_realloc(resultList, sizeof(vrna_move_t) * (count + 1));
  *length                 = count;
  return resultList;
}


/**
 * @brief Compute the freed interval after a shift move. These positions can be used to
 *        generate base shifts and insertions that were not possible before a shift move.
 *
 * @structure - the structure as pair table
 * @m - the shift move that will be performed on this structure
 * @freedInterval - output of the left and right position of the freed interval [l,r]
 * @return the type of the freed interval
 */
PRIVATE intervalType
computeFreedInterval(const short        *structure,
                     const vrna_move_t  *m,
                     vrna_move_t        *freedInterval)
{
  int           positivePosition        = MAX2(m->pos_5, m->pos_3);
  int           newPairedPosition       = abs(MIN2(m->pos_5, m->pos_3));
  int           previousPairedPosition  = structure[positivePosition];
  intervalType  t;

  /*    |  +)..-) //+ = newPairedPos; - = prevPaired; | = positivePosition (unchanged). */
  if (positivePosition < previousPairedPosition && positivePosition < newPairedPosition) {
    if (newPairedPosition < previousPairedPosition) {
      freedInterval->pos_5  = newPairedPosition + 1;
      freedInterval->pos_3  = previousPairedPosition;
      t                     = DECREASED;
    } else {
      /*    |  -)..+) */
      freedInterval->pos_5  = previousPairedPosition;
      freedInterval->pos_3  = newPairedPosition - 1;
      t                     = INCREASED;
    }
  }

  /*  (+  |..-) */
  if (positivePosition < previousPairedPosition && positivePosition > newPairedPosition) {
    freedInterval->pos_5  = positivePosition + 1;
    freedInterval->pos_3  = previousPairedPosition;
    t                     = SWITCHED;
  }

  /*  (-..|  +) */
  if (positivePosition > previousPairedPosition && positivePosition < newPairedPosition) {
    freedInterval->pos_5  = previousPairedPosition;
    freedInterval->pos_3  = positivePosition - 1;
    t                     = SWITCHED;
  }

  /*  -(..+(  | */
  if (positivePosition > previousPairedPosition && positivePosition > newPairedPosition) {
    if (newPairedPosition > previousPairedPosition) {
      freedInterval->pos_5  = previousPairedPosition;
      freedInterval->pos_3  = newPairedPosition - 1;
      t                     = DECREASED;
    } else {
      /*  +(..-(  | */
      freedInterval->pos_5  = newPairedPosition + 1;
      freedInterval->pos_3  = previousPairedPosition;
      t                     = INCREASED;
    }
  }

  return t;
}


/**
 * @brief Generate all new shift moves between the freed interval of the shift move and
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
 * @param prev_pt - the structure as pair table
 * @param curr_move - the move that can be applied to the previous structure
 * @param length - outputs the length of the output
 * @result a list with shift and insertion moves that are possible on the current structure and were not possible
 *         on the previous structure
 */
PRIVATE vrna_move_t *
generateShiftsThatWereNotPossibleBeforeThisShiftMove(const vrna_fold_compound_t  *vc,
                                                                  const short                 *prev_pt,
                                                                  const vrna_move_t           *curr_move,
                                                                  int                         *length)
{
  short         *currentStructure = vrna_ptable_copy(prev_pt);

  vrna_move_apply(currentStructure, curr_move);
  vrna_move_t   freedInterval = {
    0, 0
  };
  int           positivePosition        = MAX2(curr_move->pos_5, curr_move->pos_3);
  int           newPairedPosition       = abs(MIN2(curr_move->pos_5, curr_move->pos_3));
  int           previousPairedPosition  = prev_pt[positivePosition];
  intervalType  t                       = computeFreedInterval(prev_pt,
                                                               curr_move,
                                                               &freedInterval);

  vrna_move_t   *allNewShifts = vrna_alloc(sizeof(vrna_move_t) * (vc->length * vc->length));
  int           count         = 0;
  /* moves from all pairs inside the freed interval with previous paired position. */
  if (previousPairedPosition == freedInterval.pos_5) {
    shift_bpins_to_i_from_right(vc,
                                previousPairedPosition,
                                freedInterval.pos_5 - 1,
                                freedInterval.pos_3 + 1,
                                currentStructure,
                                allNewShifts,
                                &count);
  } else {
    shift_bpins_to_i_from_left(vc,
                               previousPairedPosition,
                               freedInterval.pos_3 + 1,
                               freedInterval.pos_5 - 1,
                               currentStructure,
                               allNewShifts,
                               &count);
  }

  if (t == INCREASED) {
    /* compute only pairs to one interval instead of two */
    int i_left;
    int i_right;
    if (positivePosition < previousPairedPosition) {
      /* [to pair][freed]   */
      i_left  = positivePosition;
      i_right = previousPairedPosition;
      pairs_to_left_most_position_whithin_eclosing_loop_and_shifts_to_interval(vc,
                                                                               previousPairedPosition,
                                                                               freedInterval.pos_5 - 1,
                                                                               freedInterval.pos_3 + 1,
                                                                               currentStructure,
                                                                               allNewShifts,
                                                                               &count,
                                                                               &shift_bpins_to_right,
                                                                               0);
      /* and now from pairs of the freed interval to the environment */
      pairs_from_interval_into_shifts_to_interval(vc,
                                                  freedInterval.pos_5,
                                                  freedInterval.pos_3,
                                                  freedInterval.pos_5,
                                                  0,
                                                  currentStructure,
                                                  allNewShifts,
                                                  &count,
                                                  &shift_bpins_to_left);
    } else {
      /* [freed][to pair] */
      i_left  = previousPairedPosition;
      i_right = positivePosition;
      pairs_to_right_most_position_whithin_eclosing_loop_and_shifts_to_interval(vc,
                                                                                previousPairedPosition,
                                                                                freedInterval.pos_3 + 1,
                                                                                freedInterval.pos_5 - 1,
                                                                                currentStructure,
                                                                                allNewShifts,
                                                                                &count,
                                                                                &shift_bpins_to_left,
                                                                                0);

      /* and now from pairs of the freed interval to the environment */
      pairs_from_interval_into_shifts_to_interval(vc,
                                                  freedInterval.pos_5 - 1,
                                                  freedInterval.pos_3,
                                                  freedInterval.pos_3,
                                                  vc->length + 1,
                                                  currentStructure,
                                                  allNewShifts,
                                                  &count,
                                                  &shift_bpins_to_right);
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

    pairs_to_left_most_position_whithin_eclosing_loop_and_shifts_to_interval(vc,
                                                                             startLeftDirection + 1,
                                                                             freedInterval.pos_5 - 1,
                                                                             freedInterval.pos_3 + 1,
                                                                             currentStructure,
                                                                             allNewShifts,
                                                                             &count,
                                                                             &shift_bpins_to_right,
                                                                             1);

    pairs_to_right_most_position_whithin_eclosing_loop_and_shifts_to_interval(vc,
                                                                              startRightDirection - 1,
                                                                              freedInterval.pos_3 + 1,
                                                                              freedInterval.pos_5 - 1,
                                                                              currentStructure,
                                                                              allNewShifts,
                                                                              &count,
                                                                              &shift_bpins_to_left,
                                                                              1);
    /* and now from pairs of the freed interval to the environment */
    pairs_from_interval_into_shifts_to_interval(vc,
                                                freedInterval.pos_5 - 1,
                                                freedInterval.pos_3 + 1,
                                                freedInterval.pos_5,
                                                0,
                                                currentStructure,
                                                allNewShifts,
                                                &count,
                                                &shift_bpins_to_left);
    pairs_from_interval_into_shifts_to_interval(vc,
                                                freedInterval.pos_5 - 1,
                                                freedInterval.pos_3 + 1,
                                                freedInterval.pos_3,
                                                vc->length + 1,
                                                currentStructure,
                                                allNewShifts,
                                                &count,
                                                &shift_bpins_to_right);
  }

  /* and finally all shifts with the new paired position as stable position */
  if (positivePosition < newPairedPosition) {
    /* exclude previous base pair --> two calls */
    shift_bpins_to_left(vc,
                        newPairedPosition,
                        newPairedPosition,
                        positivePosition,
                        currentStructure,
                        allNewShifts,
                        &count);
    shift_bpins_to_left(vc,
                        newPairedPosition,
                        positivePosition,
                        0,
                        currentStructure,
                        allNewShifts,
                        &count);

    shift_bpins_to_right(vc,
                         newPairedPosition,
                         newPairedPosition,
                         vc->length + 1,
                         currentStructure,
                         allNewShifts,
                         &count);
  } else {
    shift_bpins_to_left(vc,
                        newPairedPosition,
                        newPairedPosition,
                        0,
                        currentStructure,
                        allNewShifts,
                        &count);

    shift_bpins_to_right(vc,
                         newPairedPosition,
                         newPairedPosition,
                         positivePosition,
                         currentStructure,
                         allNewShifts,
                         &count);
    shift_bpins_to_right(vc,
                         newPairedPosition,
                         positivePosition,
                         vc->length + 1,
                         currentStructure,
                         allNewShifts,
                         &count);
  }


  free(currentStructure);
  /* add terminator */
  allNewShifts[count].pos_5  = 0;
  allNewShifts[count].pos_3  = 0;
  *length                           = count;
  return allNewShifts;
}


/**
 * @brief generate all shifts on the given structure that cross the current move.
 *        The move should already be performed on the structure.
 */
PRIVATE vrna_move_t *
generateCrossingShifts(const vrna_fold_compound_t *vc,
                       const short                *structure,
                       const vrna_move_t          *curr_move,
                       int                        *length)
{
  int         leftPosition    = MIN2(abs(curr_move->pos_3), abs(curr_move->pos_5));
  int         rightPosition   = MAX2(abs(curr_move->pos_3), abs(curr_move->pos_5));
  int         structureLength = vc->length;

  int         count = 0;

  /* Maximum memory */
  size_t      maxCrossingPairs = 2 * (rightPosition - leftPosition) *
                                 (structureLength - (rightPosition - leftPosition));
  size_t      resultSize = maxCrossingPairs *
                           (structureLength - (rightPosition - leftPosition));

  vrna_move_t *resultList = vrna_alloc(sizeof(vrna_move_t) * (resultSize + 1));

  /* left to middle */
  pairs_to_left_most_position_whithin_eclosing_loop_and_shifts_to_interval(vc,
                                                                           leftPosition,
                                                                           leftPosition - 1,
                                                                           rightPosition + 1,
                                                                           structure,
                                                                           resultList,
                                                                           &count,
                                                                           &shift_bpins_to_right,
                                                                           1);

  /* right to middle */
  pairs_to_right_most_position_whithin_eclosing_loop_and_shifts_to_interval(vc,
                                                                            rightPosition,
                                                                            rightPosition + 1,
                                                                            leftPosition - 1,
                                                                            structure,
                                                                            resultList,
                                                                            &count,
                                                                            &shift_bpins_to_left,
                                                                            1);

  /* middle to left */
  pairs_from_interval_into_shifts_to_interval(vc,
                                              leftPosition,
                                              rightPosition,
                                              leftPosition + 1,
                                              0,
                                              structure,
                                              resultList,
                                              &count,
                                              &shift_bpins_to_left);
  /* middle to right */
  pairs_from_interval_into_shifts_to_interval(vc,
                                              leftPosition,
                                              rightPosition,
                                              rightPosition - 1,
                                              structureLength + 1,
                                              structure,
                                              resultList,
                                              &count,
                                              &shift_bpins_to_right);

  resultSize = count;

  resultList[resultSize].pos_5  = 0;
  resultList[resultSize].pos_3  = 0;
  *length                       = resultSize;
  return resultList;
}


/**
 * @brief generate all insertions on the given structure that cross the current deletion move.
 *        The move should already be performed on the structure.
 */
PRIVATE vrna_move_t *
generateCrossingInserts(const vrna_fold_compound_t  *vc,
                        const short                 *structure,
                        const vrna_move_t           *curr_move,
                        int                         *length)
{
  int         leftPosition  = MIN2(abs(curr_move->pos_3), abs(curr_move->pos_5));
  int         rightPosition = MAX2(abs(curr_move->pos_3), abs(curr_move->pos_5));

  int         len           = vc->length;
  int         allocatedSize = 2 * (rightPosition - leftPosition) *
                              (len - (rightPosition - leftPosition));
  vrna_move_t *resultList = vrna_alloc(sizeof(vrna_move_t) * (allocatedSize));
  int         count       = 0;

  int         startLeft = leftPosition + 1;

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
    vrna_move_t *m = &resultList[i];
    m->pos_5  = abs(m->pos_5);
    m->pos_3  = abs(m->pos_3);
  }

  resultList              = vrna_realloc(resultList, sizeof(vrna_move_t) * (count + 1));
  resultList[count].pos_5 = 0;
  resultList[count].pos_3 = 0;
  *length                 = count;
  return resultList;
}


PRIVATE vrna_move_t *
buildNeighborsForDeletionMove(const vrna_fold_compound_t  *vc,
                              const vrna_move_t           *curr_move,
                              const short                 *prev_pt,
                              const vrna_move_t           *prev_neighbors,
                              int                         length,
                              int                         *size_neighbors,
                              unsigned int                options)
{
  /* more moves */
  vrna_move_t *newMoves = (vrna_move_t *)vrna_alloc(sizeof(vrna_move_t) * (length));
  int         newCount  = 0;

  /* insert insertion moves and deletion moves and shift moves, except curr_move and except shifts that have a common position with the deletion. */
  for (int i = 0; i < length; i++) {
    const vrna_move_t *pm = &prev_neighbors[i];
    if (!(abs(pm->pos_5) == abs(curr_move->pos_5) ||
          abs(pm->pos_3) == abs(curr_move->pos_3)
          || abs(pm->pos_5) == abs(curr_move->pos_3) ||
          abs(pm->pos_3) == abs(curr_move->pos_5))) {
      newMoves[newCount] = *pm;
      newCount++;
    }
  }
  short       *currentStructure = vrna_ptable_copy(prev_pt);
  int         leftPosition      = MIN2(abs(curr_move->pos_3), abs(curr_move->pos_5));
  int         rightPosition     = MAX2(abs(curr_move->pos_3), abs(curr_move->pos_5));
  currentStructure[leftPosition]  = 0;
  currentStructure[rightPosition] = 0;

  /* create crossing insert structures in within loopIndex with the deleted bp positions. */
  int         lengthInserts     = 0;
  vrna_move_t *crossingInserts;

  if(options & VRNA_MOVESET_INSERTION) {
  crossingInserts = generateCrossingInserts(vc,
                                                          currentStructure,
                                                          curr_move,
                                                          &lengthInserts);
  }

  int         lengthCrossingShifts  = 0;
  vrna_move_t *crossingShiftMoves   = NULL;

  if (options & VRNA_MOVESET_SHIFT) {
    /* generate shifts that cross the deleted base pair */
    crossingShiftMoves = generateCrossingShifts(vc,
                                                currentStructure,
                                                curr_move,
                                                &lengthCrossingShifts);
  }

  int totalSize = newCount + lengthCrossingShifts + lengthInserts;
  newMoves = vrna_realloc(newMoves, sizeof(vrna_move_t) * (totalSize + 1));

  if(lengthInserts > 0){
    memcpy(&newMoves[newCount], crossingInserts, sizeof(vrna_move_t) * (lengthInserts));
    newCount += lengthInserts;
  }
  if(lengthCrossingShifts > 0){
    memcpy(&newMoves[newCount], crossingShiftMoves, sizeof(vrna_move_t) * (lengthCrossingShifts));
    newCount += lengthCrossingShifts;
  }
  if (options & VRNA_MOVESET_SHIFT)
    free(crossingShiftMoves);

  if(options & VRNA_MOVESET_INSERTION)
    free(crossingInserts);
  free(currentStructure);

  *size_neighbors           = newCount;
  newMoves[newCount].pos_5  = 0;
  newMoves[newCount].pos_3  = 0;
  return newMoves;
}


PRIVATE vrna_move_t *
buildNeighborsForInsertionMove(const vrna_fold_compound_t *vc,
                               const vrna_move_t          *curr_move,
                               const short                *prev_pt,
                               const vrna_move_t          *prev_neighbors,
                               int                        length,
                               int                        *size_neighbors,
                               unsigned int               options)
{
  int         allocatedSize = length;
  vrna_move_t *newMoves     = (vrna_move_t *)vrna_alloc(sizeof(vrna_move_t) * (length));
  int         newCount      = 0;

  /* insert current move as deletion */
  vrna_move_t nm = {
    -abs(curr_move->pos_5), -abs(curr_move->pos_3)
  };

  newMoves[newCount++] = nm;

  /* copy filtered previous moves */
  for (int i = 0; i < length; i++) {
    const vrna_move_t *pm = &prev_neighbors[i];
    if (!is_crossing(abs(pm->pos_5), abs(pm->pos_3), abs(curr_move->pos_5),
                     abs(curr_move->pos_3))) {
      /* insert all moves if they are not crossing with the current move. */
      newMoves[newCount++] = *pm;
      continue;
    } else {
      if (options & VRNA_MOVESET_SHIFT) {
        /* current is insertion and previous move is crossing and previous move is insertion
         * --> convert to shift if one position is equal */
        if (pm->pos_5 > 0 && pm->pos_3 > 0) {
          vrna_move_t insertionAsShift  = *pm;
          bool        isConvertable     = false;
          if ((pm->pos_5 == curr_move->pos_5) || (pm->pos_5 == curr_move->pos_3)) {
            insertionAsShift.pos_3  = -pm->pos_3;
            isConvertable           = true;
          }

          if ((pm->pos_3 == curr_move->pos_5) || (pm->pos_3 == curr_move->pos_3)) {
            insertionAsShift.pos_5  = -pm->pos_5;
            isConvertable           = true;
          }

          if ((pm->pos_5 == curr_move->pos_5 && pm->pos_3 == curr_move->pos_3)
              || (pm->pos_3 == curr_move->pos_5 && pm->pos_5 == curr_move->pos_3))
            continue;

          if (isConvertable) {
            if (newCount >= allocatedSize) {
              allocatedSize += vc->length;
              newMoves      = vrna_realloc(newMoves, sizeof(vrna_move_t) * allocatedSize);
            }

            newMoves[newCount++] = insertionAsShift;
            continue;
          }
        }
      }
    }
  }
  *size_neighbors = newCount;
  newMoves        =
    (vrna_move_t *)vrna_realloc(newMoves, sizeof(vrna_move_t) * (newCount + 1));
  newMoves[newCount].pos_5  = 0;
  newMoves[newCount].pos_3  = 0;
  return newMoves;
}


PRIVATE vrna_move_t *
buildNeighborsForShiftMove(const vrna_fold_compound_t *vc,
                           const vrna_move_t          *curr_move,
                           const short                *prev_pt,
                           const vrna_move_t          *prev_neighbors,
                           int                        length,
                           int                        *size_neighbors,
                           unsigned int               options)
{
  int         allocatedSize = length;
  vrna_move_t *newMoves     = (vrna_move_t *)vrna_alloc(sizeof(vrna_move_t) * (allocatedSize));
  int         newCount      = 0;

  int         currentPositivePosition = curr_move->pos_5;
  int         currentNegativePosition = curr_move->pos_5;
  if (curr_move->pos_3 > 0){
    currentPositivePosition = curr_move->pos_3;
    currentNegativePosition = curr_move->pos_5;
  }
  /* copy filtered previous moves */
  int previousPairedPosition = prev_pt[currentPositivePosition];
  for (int i = 0; i < length; i++) {
    const vrna_move_t *pm     = &prev_neighbors[i];
    bool              isShift = (pm->pos_5 > 0 && pm->pos_3 < 0) ||
                                (pm->pos_5 < 0 && pm->pos_3 > 0);
    if (!is_crossing(abs(pm->pos_5), abs(pm->pos_3), abs(curr_move->pos_5),
                     abs(curr_move->pos_3))) {
      /* if current move is a shift and pm is a shift */
      if (isShift) {
        /*  current negative position cannot be pm.positive position --> allow only if current pos == shift pos. */
        int positivePosition = pm->pos_5;
        if (pm->pos_3 > 0)
          positivePosition = pm->pos_3;

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
        int positivePosition = MAX2(pm->pos_5, pm->pos_3);
        if ((positivePosition == curr_move->pos_5 || positivePosition == curr_move->pos_3)
            && !((pm->pos_5 == curr_move->pos_5 && pm->pos_3 == curr_move->pos_3)
                 || (pm->pos_3 == curr_move->pos_5 && pm->pos_5 == curr_move->pos_3))) {
          /*  if one position in common, but not the same move */
          newMoves[newCount++] = *pm;
          continue;
        }
      }
    }
  }

  /* insert current move as deletion */
  vrna_move_t nm = {
    -abs(curr_move->pos_5), -abs(curr_move->pos_3)
  };
  newMoves[newCount++] = nm;

  /* insert back shift */
  int         left      = MIN2(currentPositivePosition, previousPairedPosition);
  int         right     = MAX2(currentPositivePosition, previousPairedPosition);
  vrna_move_t backShift = {
    left, right
  };
  if (previousPairedPosition == backShift.pos_5)
    backShift.pos_5 = -backShift.pos_5;
  else
    backShift.pos_3 = -backShift.pos_3;

  newMoves[newCount++] = backShift;



  if(options & VRNA_MOVESET_INSERTION){
    vrna_move_t freedInterval;
    intervalType t = computeFreedInterval(prev_pt,curr_move,&freedInterval);
    /* insertion moves */
    int         length_insertionMoves = 0;
    short * curr_pt = vrna_ptable_copy(prev_pt);
    vrna_move_apply(curr_pt, curr_move);
    vrna_move_t *insertionMoves;
    insertionMoves                    = generateInsertionsThatWereNotPossibleBeforeThisShiftMove(vc,
                                                                                               curr_pt,
                                                                                               &freedInterval,
                                                                                               t,
                                                                                               currentPositivePosition,
                                                                                               previousPairedPosition,
                                                                                               currentNegativePosition,
                                                                                               &length_insertionMoves);
    int expectedSize = newCount + length_insertionMoves;
    if (expectedSize >= allocatedSize) {
      allocatedSize = expectedSize + 1;
      newMoves      = vrna_realloc(newMoves, sizeof(vrna_move_t) * (allocatedSize));
    }

    memcpy(&newMoves[newCount], insertionMoves,
           sizeof(vrna_move_t) * (length_insertionMoves));
    newCount += length_insertionMoves;
    free(insertionMoves);
    free(curr_pt);
  }


  /*  compute all shift moves and insertion moves that end or start in the interval, that was freed by the current shift move. */
  int         length_newFreedIntervalMoves  = 0;
  vrna_move_t *newFreedIntervalMoves        =
    generateShiftsThatWereNotPossibleBeforeThisShiftMove(vc,
                                                                      prev_pt,
                                                                      curr_move,
                                                                      &length_newFreedIntervalMoves);

  int expectedSize = newCount + length_newFreedIntervalMoves;
  if (expectedSize >= allocatedSize) {
    allocatedSize = expectedSize + 1;
    newMoves      = vrna_realloc(newMoves, sizeof(vrna_move_t) * (allocatedSize));
  }

  memcpy(&newMoves[newCount], newFreedIntervalMoves,
         sizeof(vrna_move_t) * (length_newFreedIntervalMoves));
  newCount += length_newFreedIntervalMoves;
  free(newFreedIntervalMoves);

  *size_neighbors = newCount;
  newMoves        =
    (vrna_move_t *)vrna_realloc(newMoves, sizeof(vrna_move_t) * (newCount + 1));
  newMoves[newCount].pos_5  = 0;
  newMoves[newCount].pos_3  = 0;
  return newMoves;
}


/**************************************/
/* public helper methods              */
/**************************************/

void
vrna_move_apply(short             *pt,
                const vrna_move_t *m)
{
  /* deletion */
  if (m->pos_5 < 0 && m->pos_3 < 0) {
    pt[-m->pos_5] = 0;
    pt[-m->pos_3] = 0;
    return;
  }

  /* insertion */
  if (m->pos_5 > 0 && m->pos_3 > 0) {
    pt[m->pos_5]  = m->pos_3;
    pt[m->pos_3]  = m->pos_5;
    return;
  }

  /* shift right */
  if (m->pos_5 > 0 /*&& m->pos_3 < 0*/) {
    short previousPairedPosition = pt[m->pos_5];
    pt[previousPairedPosition] = 0;
    short newPairedPosition = -m->pos_3;
    pt[m->pos_5]          = newPairedPosition;
    pt[newPairedPosition] = m->pos_5;
    return;
  }

  /* shift left */
  if (m->pos_5 < 0 /*&& m->pos_3 > 0*/) {
    short previousPairedPosition = pt[m->pos_3];
    pt[previousPairedPosition] = 0;
    short newPairedPosition = -m->pos_5;
    pt[m->pos_3]          = newPairedPosition;
    pt[newPairedPosition] = m->pos_3;
    return;
  }
}


void
vrna_move_apply_db(char               *structure,
                   const short        *pt,
                   const vrna_move_t  *m)
{
  /* deletion */
  if (m->pos_5 < 0 && m->pos_3 < 0) {
    structure[-m->pos_5 - 1]  = '.';
    structure[-m->pos_3 - 1]  = '.';
    return;
  }

  /* insertion */
  if (m->pos_5 > 0 && m->pos_3 > 0) {
    structure[m->pos_5 - 1] = '(';
    structure[m->pos_3 - 1] = ')';
    return;
  }

  /* shift right */
  if (m->pos_5 > 0) {
    short previousPairedPosition = pt[m->pos_5];
    structure[previousPairedPosition - 1] = '.';
    int   left  = m->pos_5 - 1;
    int   right = -m->pos_3 - 1;
    structure[left]   = '(';
    structure[right]  = ')';
    return;
  }

  /* shift left */
  if (m->pos_5 < 0) {
    short previousPairedPosition = pt[m->pos_3];
    structure[previousPairedPosition - 1] = '.';
    int   left  = -m->pos_5 - 1;
    int   right = m->pos_3 - 1;
    structure[left]   = '(';
    structure[right]  = ')';
    return;
  }
}


void
vrna_loopidx_update(int               *loopIndices,
                    const short       *pt,
                    int               length,
                    const vrna_move_t *m)
{
  int index = loopIndices[abs(m->pos_5)];
  int pos_5 = m->pos_5;
  int pos_3 = m->pos_3;

  /* left index has to be smaller for the next iterations */
  if (abs(pos_5) > abs(pos_3)) {
    pos_5 = m->pos_3;
    pos_3 = m->pos_5;
  }

  if (pos_5 < 0 && pos_3 < 0) {
    /* deletion -> decrease loopIndex. */
    /* search for next closing bracket or unpaired loopIndex to the right. */
    int currentIndex = 0;
    for (int i = -pos_3 + 1; i <= length; i++) {
      if (loopIndices[i] < index) {
        /* if unpaired or closing. an loopIndex < current. */
        if (pt[i] <= 0 || pt[i] < i) {
          currentIndex = loopIndices[i];
          break;
        }
      }
    }

    int highestIndexInside = index;
    for (int i = -pos_5; i <= -pos_3; i++) {
      if (loopIndices[i] > highestIndexInside)
        highestIndexInside = loopIndices[i];

      if (loopIndices[i] == index)
        loopIndices[i] = currentIndex;
      else
        loopIndices[i]--;
    }
    /*  decrease all following indices. */
    for (int i = -pos_3 + 1; i <= length; i++)
      if (loopIndices[i] > highestIndexInside)
        loopIndices[i]--;

    loopIndices[0]--;
    return;
  } else if (pos_5 > 0 && pos_3 > 0) {
    /* insertion --> 1.Increase the index for all following loopIndices inside the interval. */
    /* 2.Then look at the index on the left side and increase it for the base pair interval (ignore bp in between). */
    /* 3.Increase the indices to the right side of the interval. */
    int currentIndex = 0;
    for (int i = pos_5; i > 0; i--) {
      if (pt[i] > i) {
        /*  if bp opens */
        currentIndex = loopIndices[i];
        break;
      }
    }
    currentIndex++;
    for (int i = pos_5; i <= pos_3; i++) {
      if (loopIndices[i] >= currentIndex)
        loopIndices[i]++;
      else
        loopIndices[i] = currentIndex;
    }
    for (int i = pos_3 + 1; i <= length; i++)
      if (loopIndices[i] >= currentIndex)
        loopIndices[i]++;

    loopIndices[0]++;
    return;
  } else {
    /* shift move: split it into one deletion and one insertion. */
    vrna_move_t deletionMove;
    if (pos_5 > 0) {
      deletionMove.pos_5  = -pos_5;
      deletionMove.pos_3  = -pt[pos_5];
    } else {
      deletionMove.pos_3  = -pos_3;
      deletionMove.pos_5  = -pt[pos_3];
    }

    /* left index absolute value has to be smaller than the right one for the other methods */
    if (deletionMove.pos_5 < deletionMove.pos_3) {
      int tmp = deletionMove.pos_5;
      deletionMove.pos_5  = deletionMove.pos_3;
      deletionMove.pos_3  = tmp;
    }

    vrna_move_t insertionMove = {
      abs(pos_5), abs(pos_3)
    };

    vrna_loopidx_update(loopIndices, pt, length, &deletionMove);
    short       *tmppt = vrna_ptable_copy(pt);
    vrna_move_apply(tmppt, &deletionMove);
    vrna_loopidx_update(loopIndices, tmppt, length, &insertionMove);
    free(tmppt);
    return;
  }
}


/**************************************/
/* public neighbor methods            */
/**************************************/

vrna_move_t *
vrna_neighbors(vrna_fold_compound_t *vc,
               const short          *pt,
               unsigned int         options)
{
  vrna_move_t *moveSet = NULL;
  int totalLength = 0;
  vrna_move_t *ptMoveSetEnd;
  int         lengthDeletions=0;
  int         lengthInsertions=0;

  if (options & VRNA_MOVESET_DELETION) {
    /*  append deletion moves */
    vrna_move_t *deletionList = deletions(vc, pt, &lengthDeletions);
    totalLength += lengthDeletions;
    moveSet     = (vrna_move_t *)vrna_realloc(moveSet, sizeof(vrna_move_t) * (totalLength + 1));
    ptMoveSetEnd = moveSet;
    memcpy(ptMoveSetEnd, deletionList, lengthDeletions * sizeof(vrna_move_t));
    free(deletionList);
  }

  if (options & VRNA_MOVESET_INSERTION) {
    /*  append insertion moves */
    vrna_move_t *insertionList = insertions(vc, pt, &lengthInsertions);
    totalLength += lengthInsertions;
    moveSet     = (vrna_move_t *)vrna_realloc(moveSet, sizeof(vrna_move_t) * (totalLength + 1));
    ptMoveSetEnd = moveSet + lengthDeletions;
    memcpy(ptMoveSetEnd, insertionList, lengthInsertions * sizeof(vrna_move_t));
    free(insertionList);
  }

  if (options & VRNA_MOVESET_SHIFT) {
    /*  append shift moves */
    int         lengthShifts;
    vrna_move_t *shiftList = shifts(vc, pt, &lengthShifts);
    totalLength += lengthShifts;
    moveSet     = (vrna_move_t *)vrna_realloc(moveSet, sizeof(vrna_move_t) * (totalLength + 1));
    ptMoveSetEnd = moveSet + lengthInsertions + lengthDeletions;
    memcpy(ptMoveSetEnd, shiftList, lengthShifts * sizeof(vrna_move_t));
    free(shiftList);
  }

  if(totalLength > 0){
    /* terminate list */
    moveSet[totalLength].pos_5  = 0;
    moveSet[totalLength].pos_3  = 0;
  }
  return moveSet;
}


vrna_move_t *
vrna_neighbors_successive(const vrna_fold_compound_t  *vc,
                          const vrna_move_t           *curr_move,
                          const short                 *prev_pt,
                          const vrna_move_t           *prev_neighbors,
                          int                         length,
                          int                         *size_neighbors,
                          unsigned int                options)
{
  bool        isDeletion  = curr_move->pos_5 < 0 && curr_move->pos_3 < 0 ? true : false;
  bool        isInsertion = curr_move->pos_5 > 0 && curr_move->pos_3 > 0 ? true : false;
  bool        isShift     = !isDeletion && !isInsertion;

  vrna_move_t *newMoves = NULL;

  if (isDeletion) {
    newMoves = buildNeighborsForDeletionMove(vc,
                                             curr_move,
                                             prev_pt,
                                             prev_neighbors,
                                             length,
                                             size_neighbors,
                                             options);
  }

  if (isInsertion) {
    newMoves = buildNeighborsForInsertionMove(vc,
                                              curr_move,
                                              prev_pt,
                                              prev_neighbors,
                                              length,
                                              size_neighbors,
                                              options);
  }

  if (isShift) {
    newMoves = buildNeighborsForShiftMove(vc,
                                          curr_move,
                                          prev_pt,
                                          prev_neighbors,
                                          length,
                                          size_neighbors,
                                          options);
  }

  return newMoves;
}
