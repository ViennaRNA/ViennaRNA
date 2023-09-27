#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "ViennaRNA/utils/basic.h"
#include "ViennaRNA/fold_compound.h"
#include "ViennaRNA/landscape/move.h"
#include "ViennaRNA/landscape/neighbor.h"

#ifdef __GNUC__
# define INLINE inline
#else
# define INLINE
#endif
/*
 * Below follows a callback-based fast version that generates a list of moves to neighboring structures
 * that appear after application of a particular move. In particular, the implementation below reports
 * two list, one for all neighbors that become available, and one for all the moves that become invalid
 * after application of the move.
 */


#define FOR_UNPAIRED(ptable, idx, start, cond, next, code)    { \
    for ((start); (cond); (next)) { \
      if (ptable[idx] > idx) { \
        idx = ptable[idx]; \
      } else { code } \
    } }

#define FOR_PAIRED(ptable, idx, start, cond, next, code)      { \
    for ((start); (cond); (next)) { \
      if (ptable[idx] > idx) { \
        { code }; \
        idx = ptable[idx]; \
      } \
    } }


#define ST_NEW    VRNA_NEIGHBOR_NEW
#define ST_CHG    VRNA_NEIGHBOR_CHANGE
#define ST_DEL    VRNA_NEIGHBOR_INVALID

PRIVATE INLINE void
enclosing_pair(const short  *pt,
               int          pair_pos_5,
               int          pair_pos_3,
               int          *enc_pos_5,
               int          *enc_pos_3);


PRIVATE void
generate_local_nb(vrna_fold_compound_t  *fc,
                  const short           *pt,
                  const vrna_move_t     *move,
                  const vrna_move_t     *affected_pair,
                  vrna_move_update_f    cb,
                  void                  *data,
                  unsigned int          options);


PRIVATE void
generate_local_nb_insertion(vrna_fold_compound_t  *fc,
                            const short           *pt,
                            const vrna_move_t     *move,
                            const vrna_move_t     *affected_pair,
                            vrna_move_update_f    cb,
                            void                  *data,
                            unsigned int          options);


PRIVATE void
generate_local_nb_deletion(vrna_fold_compound_t *fc,
                           const short          *pt,
                           const vrna_move_t    *move,
                           const vrna_move_t    *affected_pair,
                           vrna_move_update_f   cb,
                           void                 *data,
                           unsigned int         options);


PRIVATE void
generate_local_nb_shift(vrna_fold_compound_t  *fc,
                        const short           *pt,
                        const vrna_move_t     *move,
                        const vrna_move_t     *affected_pair,
                        vrna_move_update_f    cb,
                        void                  *data,
                        unsigned int          options);


PRIVATE void
generate_conflicts_local_nb(vrna_fold_compound_t  *fc,
                            const short           *pt,
                            const vrna_move_t     *move,
                            vrna_move_update_f    cb,
                            void                  *data,
                            unsigned int          options);


PRIVATE void
generate_conflicts_local_nb_insertion(vrna_fold_compound_t  *fc,
                                      const short           *pt,
                                      const vrna_move_t     *move,
                                      vrna_move_update_f    cb,
                                      void                  *data,
                                      unsigned int          options);


PRIVATE void
generate_conflicts_local_nb_deletion(vrna_fold_compound_t *fc,
                                     const short          *pt,
                                     const vrna_move_t    *move,
                                     vrna_move_update_f   cb,
                                     void                 *data,
                                     unsigned int         options);


PRIVATE void
generate_conflicts_local_nb_shift(vrna_fold_compound_t  *fc,
                                  const short           *pt,
                                  const vrna_move_t     *move,
                                  vrna_move_update_f    cb,
                                  void                  *data,
                                  unsigned int          options);


#include "landscape/local_neighbors.inc"

PRIVATE INLINE int
is_compatible(const vrna_fold_compound_t  *vc,
              int                         i,
              int                         j)
{
  /* better use hard constraints here! */
  if (i > j) {
    int k = i;
    i = j;
    j = k;
  }

  if (i + vc->params->model_details.min_loop_size < j)
    return vc->params->model_details.pair[vc->sequence_encoding2[i]][vc->sequence_encoding2[j]] !=
           0;                                                                                         /* see pair_mat.h */

  return 0;
}


PUBLIC int
vrna_move_neighbor_diff_cb(vrna_fold_compound_t *fc,
                           short                *ptable,
                           vrna_move_t          move,
                           vrna_move_update_f   cb,
                           void                 *data,
                           unsigned int         options)
{
  if ((fc) && (ptable) && (cb)) {
    /* crude check if ptable has correct size */
    if ((unsigned int)ptable[0] == fc->length) {
      /* save affected pair in case we loose it due to application of move later on */
      vrna_move_t affected_pair = vrna_move_init(move.pos_5, move.pos_3);

      if ((affected_pair.pos_5 < 0) &&
          (affected_pair.pos_3 > 0)) {
        /* shift move, 5' position changes, 3' position stays constant */
        affected_pair.pos_5 = ptable[affected_pair.pos_3];
      } else if ((affected_pair.pos_5 > 0) &&
                 (affected_pair.pos_3 < 0)) {
        /* shift move, 3' position changes, 5' position stays constant */
        affected_pair.pos_3 = ptable[affected_pair.pos_5];
      } else if (affected_pair.pos_5 < 0) {
        /* deletion move */
        affected_pair.pos_5 = -affected_pair.pos_5;
        affected_pair.pos_3 = -affected_pair.pos_3;
      }

      if (affected_pair.pos_5 > affected_pair.pos_3)
        affected_pair = vrna_move_init(affected_pair.pos_3, affected_pair.pos_5);

      /* 1. remove the neighbor we are about to change into */
      cb(fc, move, ST_DEL, data);

      /* 2. remove neighbors that become invalid after application of 'move' */
      generate_conflicts_local_nb(fc, ptable, &move, cb, data, options);

      /*
       *  we apply the move here such that novel and changed loops
       *  can be evaluated correctly in the calling context
       */
      vrna_move_apply(ptable, &move);

      /* 3. Detect novel neighbors and those that require an update after application of 'move' */
      generate_local_nb(fc, ptable, &move, &affected_pair, cb, data, options);

      /*
       *  undo the move if required
       */
      if (options & VRNA_MOVE_NO_APPLY) {
        if ((move.pos_5 < 0) &&
            (move.pos_3 > 0)) {
          /* shift move, 5' position changes, 3' position stays constant */
          ptable[-move.pos_5]         = 0;
          ptable[move.pos_3]          = affected_pair.pos_5;
          ptable[affected_pair.pos_5] = move.pos_3;
        } else if ((affected_pair.pos_5 > 0) &&
                   (affected_pair.pos_3 < 0)) {
          /* shift move, 3' position changes, 5' position stays constant */
          ptable[-move.pos_3]         = 0;
          ptable[move.pos_5]          = affected_pair.pos_3;
          ptable[affected_pair.pos_3] = move.pos_5;
        } else {
          vrna_move_t m = vrna_move_init(-move.pos_5, -move.pos_3);
          vrna_move_apply(ptable, &m);
        }
      }

      return 1; /* success */
    }
  }

  return 0;
}


PUBLIC vrna_move_t *
vrna_move_neighbor_diff(vrna_fold_compound_t  *fc,
                        short                 *ptable,
                        vrna_move_t           move,
                        vrna_move_t           **invalid_moves,
                        unsigned int          options)
{
  struct movelist *mlist;
  vrna_move_t     *valid_neighbors;

  valid_neighbors = NULL;

  if ((fc) && (ptable)) {
    mlist = init_incremental_movelist(42);

    if (invalid_moves)
      *invalid_moves = NULL;

    if (vrna_move_neighbor_diff_cb(fc,
                                   ptable,
                                   move,
                                   &add_to_incremental_move_list,
                                   (void *)mlist, options)) {
      /* prepare output */
      valid_neighbors = mlist->moves;

      /* shrink list to actually required size */
      valid_neighbors = (vrna_move_t *)vrna_realloc(valid_neighbors,
                                                    sizeof(vrna_move_t) *
                                                    (mlist->num_moves + 1));

      /* set end-of-list marker */
      valid_neighbors[mlist->num_moves] = vrna_move_init(0, 0);

      if (invalid_moves) {
        *invalid_moves = mlist->moves_invalid;

        /* shrink list to actually required size */
        *invalid_moves = (vrna_move_t *)vrna_realloc(*invalid_moves,
                                                     sizeof(vrna_move_t) *
                                                     (mlist->num_moves_invalid + 1));

        /* set end-of-list marker */
        (*invalid_moves)[mlist->num_moves_invalid] = vrna_move_init(0, 0);
      } else {
        free(mlist->moves_invalid);
      }

      mlist->moves          = NULL;
      mlist->moves_invalid  = NULL;

      free_incremental_movelist(mlist);

      return valid_neighbors;
    }

    free_incremental_movelist(mlist);
  }

  if (invalid_moves)
    *invalid_moves = NULL;

  return valid_neighbors;
}


PRIVATE void
generate_local_nb(vrna_fold_compound_t  *fc,
                  const short           *pt,
                  const vrna_move_t     *move,
                  const vrna_move_t     *affected_pair,
                  vrna_move_update_f    cb,
                  void                  *data,
                  unsigned int          options)
{
  /* check whether move is valid given the options */
  if ((move->pos_5 > 0) &&
      (move->pos_3 > 0) &&
      (options & VRNA_MOVESET_INSERTION))
    generate_local_nb_insertion(fc, pt, move, affected_pair, cb, data, options);
  else if ((move->pos_5 < 0) &&
           (move->pos_3 < 0) &&
           (options & VRNA_MOVESET_DELETION))
    generate_local_nb_deletion(fc, pt, move, affected_pair, cb, data, options);
  else if (options & VRNA_MOVESET_SHIFT)
    generate_local_nb_shift(fc, pt, move, affected_pair, cb, data, options);
}


PRIVATE void
generate_conflicts_local_nb(vrna_fold_compound_t  *fc,
                            const short           *pt,
                            const vrna_move_t     *move,
                            vrna_move_update_f    cb,
                            void                  *data,
                            unsigned int          options)
{
  /* check whether move is valid given the options */
  if ((move->pos_5 > 0) &&
      (move->pos_3 > 0) &&
      (options & VRNA_MOVESET_INSERTION)) {
    generate_conflicts_local_nb_insertion(fc, pt, move, cb, data, options);
  } else if ((move->pos_5 < 0) &&
             (move->pos_3 < 0) &&
             (options & VRNA_MOVESET_DELETION)) {
    generate_conflicts_local_nb_deletion(fc, pt, move, cb, data, options);
  } else if (options & VRNA_MOVESET_SHIFT) {
    generate_conflicts_local_nb_shift(fc, pt, move, cb, data, options);
  }
}


PRIVATE INLINE void
insertions(vrna_fold_compound_t *fc,
           const short          *pt,
           unsigned int         first_i,
           unsigned int         last_i,
           unsigned int         first_j,
           unsigned int         last_j,
           unsigned int         status,
           vrna_move_update_f   cb,
           void                 *data)
{
  unsigned int i, j;

  if (first_j == 0) {
    FOR_UNPAIRED(pt, i, i = first_i, i <= last_i, i++, {
      FOR_UNPAIRED(pt, j, j = i + 1, j <= last_j, j++, {
        if (is_compatible(fc, i, j))
          cb(fc, vrna_move_init(i, j), status, data);
      });
    });
  } else {
    FOR_UNPAIRED(pt, i, i = first_i, i <= last_i, i++, {
      FOR_UNPAIRED(pt, j, j = first_j, j <= last_j, j++, {
        if (is_compatible(fc, i, j))
          cb(fc, vrna_move_init(i, j), status, data);
      });
    });
  }
}


PRIVATE INLINE void
deletion_range(vrna_fold_compound_t *fc,
               const short          *pt,
               unsigned int         start,
               unsigned int         stop,
               unsigned int         status,
               vrna_move_update_f   cb,
               void                 *data)
{
  unsigned int i;

  FOR_PAIRED(pt, i, i = start, i <= stop, i++, {
    cb(fc, vrna_move_init(-i, -pt[i]), status, data);
  });
}


PRIVATE INLINE void
shift_pos(vrna_fold_compound_t  *fc,
          const short           *pt,
          unsigned int          i,
          unsigned int          start,
          unsigned int          end,
          unsigned int          status,
          vrna_move_update_f    cb,
          void                  *data)
{
  unsigned int k;

  if (end < i) {
    FOR_UNPAIRED(pt, k, k = start, k <= end, k++, {
      if (is_compatible(fc, k, i))
        cb(fc, vrna_move_init(-k, i), status, data);
    });
  } else {
    FOR_UNPAIRED(pt, k, k = start, k <= end, k++, {
      if (is_compatible(fc, i, k))
        cb(fc, vrna_move_init(i, -k), status, data);
    });
  }
}


PRIVATE INLINE void
shift_both(vrna_fold_compound_t *fc,
           const short          *pt,
           unsigned int         i,
           unsigned int         j,
           unsigned int         start,
           unsigned int         end,
           unsigned int         status,
           vrna_move_update_f   cb,
           void                 *data)
{
  unsigned int k;

  if (end < i) {
    FOR_UNPAIRED(pt, k, k = start, k <= end, k++, {
      if (is_compatible(fc, k, i))
        cb(fc, vrna_move_init(-k, i), status, data);

      if (is_compatible(fc, k, j))
        cb(fc, vrna_move_init(-k, j), status, data);
    });
  } else if (start < j) {
    FOR_UNPAIRED(pt, k, k = start, k <= end, k++, {
      if (is_compatible(fc, i, k))
        cb(fc, vrna_move_init(i, -k), status, data);

      if (is_compatible(fc, k, j))
        cb(fc, vrna_move_init(-k, j), status, data);
    });
  } else {
    FOR_UNPAIRED(pt, k, k = start, k <= end, k++, {
      if (is_compatible(fc, i, k))
        cb(fc, vrna_move_init(i, -k), status, data);

      if (is_compatible(fc, j, k))
        cb(fc, vrna_move_init(j, -k), status, data);
    });
  }
}


PRIVATE void
generate_local_nb_insertion(vrna_fold_compound_t  *fc,
                            const short           *pt,
                            const vrna_move_t     *move,
                            const vrna_move_t     *affected_pair,
                            vrna_move_update_f    cb,
                            void                  *data,
                            unsigned int          options)
{
  unsigned int  n = fc->length;
  int           i, j, k, mi, mj, enc5, enc3;
  int           min_loop_size = fc->params->model_details.min_loop_size;

  mi  = affected_pair->pos_5;
  mj  = affected_pair->pos_3;

  /* find out actual enclosing pair that delimits the loop affected by the current move */
  enclosing_pair(pt, mi, mj, &enc5, &enc3);

  /* 1. valid base pair deletions */
  if (options & VRNA_MOVESET_DELETION) {
    /* 1.1 removal of enclosing pair */
    if (enc5 > 0)
      cb(fc, vrna_move_init(-enc5, -enc3), ST_CHG, data);

    /* 1.2 removal of base pair that has just been inserted */
    cb(fc, vrna_move_init(-mi, -mj), ST_NEW, data);

    /* 1.3 removal of other base pairs that branch-off from the current loop */
    deletion_range(fc, pt, enc5 + 1, mi - 1, ST_CHG, cb, data);
    deletion_range(fc, pt, mi + 1, mj - 1, ST_CHG, cb, data);
    deletion_range(fc, pt, mj + 1, enc3 - 1, ST_CHG, cb, data);
  }

  /* 2. valid base pair insertions */
  if (options & VRNA_MOVESET_INSERTION) {
    /* 5' side of newly inserted base pair */

    insertions(fc, pt, enc5 + 1, mi - 1, 0, mi - 1, ST_CHG, cb, data);
    insertions(fc, pt, enc5 + 1, mi - 1, mj + 1, enc3 - 1, ST_CHG, cb, data);

    /* within newly inserted base pair */
    insertions(fc, pt, mi + 1, mj - 1, 0, mj - 1, ST_CHG, cb, data);

    /* 3' side of newly inserted base pair */
    insertions(fc, pt, mj + 1, enc3 - 1, 0, enc3 - 1, ST_CHG, cb, data);
  }

  if (options & VRNA_MOVESET_SHIFT) {
    /* 1. novel base pair shift moves due to inserted pair */

    /* 1.1 shift either mi or mj towards 5' side
     * interval [enc5 + 1:mi-1]
     */
    shift_both(fc, pt, mi, mj, enc5 + 1, mi - 1, ST_NEW, cb, data);

    /* 1.2 shift either mi or mj into the
     * enclosed interval [mi+1:mj-1]
     */
    shift_both(fc, pt, mi, mj, mi + 1, mj - 1, ST_NEW, cb, data);

    /* 1.3 shift either mi or mj towards the 3' side
     * interval [mj+1:enc3 - 1]
     */
    shift_both(fc, pt, mi, mj, mj + 1, enc3 - 1, ST_NEW, cb, data);

    /* 2. shift moves that need to be re-evaluated */

    /* 2.1 Shifts of the enclosing pair, if any */
    if (enc5 > 0) {
      shift_both(fc, pt, enc5, enc3, enc5 + 1, mi - 1, ST_CHG, cb, data);
      shift_both(fc, pt, enc5, enc3, mj + 1, enc3 - 1, ST_CHG, cb, data);
    }

    /* 2.2 shift of base pairs in same loop outside of pair inserted by move (5' side) */
    FOR_PAIRED(pt, i, i = enc5 + 1, i < mi, i++, {
      /* 2.2.1 shifts within the same loop (5' side) */
      shift_both(fc, pt, i, pt[i], enc5 + 1, i - 1, ST_CHG, cb, data);

      /* 2.2.2 shifts within the loop enclosed by (i, j) */
      shift_both(fc, pt, i, pt[i], i + 1, pt[i] - 1, ST_CHG, cb, data);

      /* 2.2.3 shifts within the same loop (3' side) */
      shift_both(fc, pt, i, pt[i], pt[i] + 1, mi - 1, ST_CHG, cb, data);
      shift_both(fc, pt, i, pt[i], mj + 1, enc3 - 1, ST_CHG, cb, data);
    });

    /* 2.3 shift of base pairs in same loop outside of pair inserted by move (3' side) */
    FOR_PAIRED(pt, i, i = mj + 1, i < enc3, i++, {
      /* 2.3.1 shifts within the same loop (5' side) */
      shift_both(fc, pt, i, pt[i], enc5 + 1, mi - 1, ST_CHG, cb, data);
      shift_both(fc, pt, i, pt[i], mj + 1, i - 1, ST_CHG, cb, data);

      /* 2.3.2 shifts within the loop enclosed by (i, j) */
      shift_both(fc, pt, i, pt[i], i + 1, pt[i] - 1, ST_CHG, cb, data);

      /* 2.3.3 shifts within the same loop (3' side) */
      shift_both(fc, pt, i, pt[i], pt[i] + 1, enc3 - 1, ST_CHG, cb, data);
    });

    /* 2.4 shifts of base pairs enclosed by novel pair introduced by move */
    FOR_PAIRED(pt, i, i = mi + 1, i < mj, i++, {
      /* 2.4.1 shifts within the loop enclosed by novel pair (5' side) */
      shift_both(fc, pt, i, pt[i], mi + 1, i - 1, ST_CHG, cb, data);

      /* 2.4.2 shifts within the loop enclosed by (i,j) */
      shift_both(fc, pt, i, pt[i], i + 1, pt[i] - 1, ST_CHG, cb, data);

      /* 2.4.3 shifts within the loop enclosed by novel pair (3' side) */
      shift_both(fc, pt, i, pt[i], pt[i] + 1, mj - 1, ST_CHG, cb, data);
    });
  }
}


PRIVATE void
generate_local_nb_deletion(vrna_fold_compound_t *fc,
                           const short          *pt,
                           const vrna_move_t    *move,
                           const vrna_move_t    *affected_pair,
                           vrna_move_update_f   cb,
                           void                 *data,
                           unsigned int         options)
{
  unsigned int  n = fc->length;
  int           i, j, k, enc5, enc3;
  int           min_loop_size = fc->params->model_details.min_loop_size;
  int           mi            = affected_pair->pos_5;
  int           mj            = affected_pair->pos_3;

  /* find out actual enclosing pair that delimits the loop affected by the current move */
  enclosing_pair(pt, mi, mj, &enc5, &enc3);

  /* 1. valid base pair deletions */
  if (options & VRNA_MOVESET_DELETION) {
    /* 1.1 removal of enclosing pair */
    if (enc5 > 0)
      cb(fc, vrna_move_init(-enc5, -enc3), ST_CHG, data);

    /* 1.2 branching base pairs 5' to base pair deleted by current move */
    deletion_range(fc, pt, enc5 + 1, mi - 1, ST_CHG, cb, data);

    /* 1.3 branching base pairs enclosed by base pair deleted by current move */
    deletion_range(fc, pt, mi + 1, mj - 1, ST_CHG, cb, data);

    /* 1.4 branching base pairs on 3' side of base pair deleted by current move */
    deletion_range(fc, pt, mj + 1, enc3 - 1, ST_CHG, cb, data);
  }

  /* 2. valid base pair insertions */
  if (options & VRNA_MOVESET_INSERTION) {
    /* 2.1 re-insertion of base pair removed by current move */
    cb(fc, vrna_move_init(mi, mj), ST_NEW, data);

    /* 2.2 insertion of novel pairs that start on 5' side of current move */

    /* 2.2.1 base pairs that end before 5' side of current move */
    insertions(fc, pt, enc5 + 1, mi - 1, 0, mi - 1, ST_CHG, cb, data);
    /* 2.2.2 base pairs that end within interval spanned by the base
     * pair we just deleted
     */
    insertions(fc, pt, enc5 + 1, mi - 1, mi, mj, ST_NEW, cb, data);
    /* 2.2.3 base pairs that end after 3' nucleotide of current move */
    insertions(fc, pt, enc5 + 1, mi - 1, mj + 1, enc3 - 1, ST_CHG, cb, data);

    /* 2.3 insertion of novel pairs that start at 5' nucleotide of current position */

    /* 2.3.1 ending within loop enclosed by base pair deleted by current move */
    insertions(fc, pt, mi, mi, 0, mj - 1, ST_NEW, cb, data);

    /* 2.3.2 ending after loop enclosed by base pair deleted by current move */
    insertions(fc, pt, mi, mi, mj + 1, enc3 - 1, ST_NEW, cb, data);

    /* 2.4 insertion of novel pairs that start within loop closed by current move */
    insertions(fc, pt, mi + 1, mj - 1, 0, mj - 1, ST_CHG, cb, data);
    insertions(fc, pt, mi + 1, mj - 1, mj, enc3 - 1, ST_NEW, cb, data);

    /* 2.5 insertion of novel pairs that start at 3' nucleotide of current move */
    insertions(fc, pt, mj, mj, 0, enc3 - 1, ST_NEW, cb, data);

    /* 2.6 insertion of novel pairs that start after 3' nucleotide of current move */
    insertions(fc, pt, mj + 1, enc3 - 1, 0, enc3 - 1, ST_CHG, cb, data);
  }

  /* 3. valid shift moves */
  if (options & VRNA_MOVESET_SHIFT) {
    /* 3.1 shifts of the enclosing pair, if any */
    if (enc5 > 0) {
      shift_both(fc, pt, enc5, enc3, enc5 + 1, mi - 1, ST_CHG, cb, data);
      shift_both(fc, pt, enc5, enc3, mi, mj, ST_NEW, cb, data);
      shift_both(fc, pt, enc5, enc3, mj + 1, enc3 - 1, ST_CHG, cb, data);
    }

    /* 3.2 shifts of base pair previously outside of removed pair 5' side */
    FOR_PAIRED(pt, i, i = enc5 + 1, i < mi, i++, {
      shift_both(fc, pt, i, pt[i], enc5 + 1, i - 1, ST_CHG, cb, data);
      shift_both(fc, pt, i, pt[i], i + 1, pt[i] - 1, ST_CHG, cb, data);
      shift_both(fc, pt, i, pt[i], pt[i] + 1, mi - 1, ST_CHG, cb, data);
      shift_both(fc, pt, i, pt[i], mi, mj, ST_NEW, cb, data);
      shift_both(fc, pt, i, pt[i], mj + 1, enc3 - 1, ST_CHG, cb, data);
    });

    /* 3.2 shifts of base pair previously outside of removed pair 3' side */
    FOR_PAIRED(pt, i, i = mj + 1, i < enc3, i++, {
      shift_both(fc, pt, i, pt[i], enc5 + 1, mi - 1, ST_CHG, cb, data);
      shift_both(fc, pt, i, pt[i], mi, mj, ST_NEW, cb, data);
      shift_both(fc, pt, i, pt[i], mj + 1, i - 1, ST_CHG, cb, data);
      shift_both(fc, pt, i, pt[i], i + 1, pt[i] - 1, ST_CHG, cb, data);
      shift_both(fc, pt, i, pt[i], pt[i] + 1, enc3 - 1, ST_CHG, cb, data);
    });

    /* 3.3. shifts of base pairs previously enclosed by removed pair */
    FOR_PAIRED(pt, i, i = mi + 1, i < mj, i++, {
      /* 3.3.1 shifts outside towards 5' side */
      shift_both(fc, pt, i, pt[i], enc5 + 1, mi, ST_NEW, cb, data);
      shift_both(fc, pt, i, pt[i], mi + 1, i - 1, ST_CHG, cb, data);

      /* 3.3.2 shifts inside (i, j) */
      shift_both(fc, pt, i, pt[i], i + 1, pt[i] - 1, ST_CHG, cb, data);

      /* 3.3.3 shifts outside towards 3' side */
      shift_both(fc, pt, i, pt[i], pt[i] + 1, mj - 1, ST_CHG, cb, data);
      shift_both(fc, pt, i, pt[i], mj, enc3 - 1, ST_NEW, cb, data);
    });
  }
}


PRIVATE void
generate_local_nb_shift(vrna_fold_compound_t  *fc,
                        const short           *pt,
                        const vrna_move_t     *move,
                        const vrna_move_t     *affected_pair,
                        vrna_move_update_f    cb,
                        void                  *data,
                        unsigned int          options)
{
  int i, j, k, mi, mj, old_i, old_j, enc5, enc3;

  old_i = affected_pair->pos_5;
  old_j = affected_pair->pos_3;

  if (move->pos_5 < 0) {
    mi  = -move->pos_5;
    mj  = move->pos_3;
  } else {
    mi  = move->pos_5;
    mj  = -move->pos_3;
  }

  if (mj < mi) {
    k   = mi;
    mi  = mj;
    mj  = k;
  }

  /* find out actual enclosing pair that delimits the loop affected by the current move */
  enclosing_pair(pt, mi, mj, &enc5, &enc3);

  /* 1. updates of insertion moves */
  if (options & VRNA_MOVESET_INSERTION) {
    /* 1.1 new insertions using the now vacant segment */
    k = old_i;
    if ((mi == k) ||
        (mj == k))
      k = old_j;

    /* base pairs still outside of (mi, mj) at 5' side */
    insertions(fc, pt, enc5 + 1, k - 1, k, mi - 1, ST_NEW, cb, data);
    insertions(fc, pt, k, k, 0, mi - 1, ST_NEW, cb, data);
    insertions(fc, pt, k, mi - 1, mj + 1, enc3 - 1, ST_NEW, cb, data);

    /* base pairs still outside of (mi, mj) at 3' side */
    insertions(fc, pt, enc5 + 1, mi - 1, mj + 1, k, ST_NEW, cb, data);
    insertions(fc, pt, mj + 1, k, k, enc3 - 1, ST_NEW, cb, data);

    /* base pairs inside of pair */
    insertions(fc, pt, mi + 1, k, k + 1, mj - 1, ST_NEW, cb, data);
    insertions(fc, pt, mi + 1, k - 1, k, mj - 1, ST_NEW, cb, data);

    /* 1.2 changed insertion outside shifted pair (i on 5' side) */
    insertions(fc, pt, enc5 + 1, MIN2(k, mi) - 1, 0, MIN2(k, mi) - 1, ST_CHG, cb, data);
    insertions(fc, pt, enc5 + 1, MIN2(k, mi) - 1, MAX2(k, mj) + 1, enc3 - 1, ST_CHG, cb, data);

    /* 1.3 insertions outside shifted pair (i on 3' side) */
    insertions(fc, pt, MAX2(k, mj) + 1, enc3 - 1, 0, enc3 - 1, ST_CHG, cb, data);

    /* 1.4 changed insertions within shifted pair. */

    /* 1.4.1 base pair insertions that are now located outside of
     * shifted pair and at 5' side
     */
    insertions(fc, pt, k + 1, mi - 1, 0, mj - 1, ST_CHG, cb, data);

    /* 1.4.2 base pair insertions that are now located outside of
     * shifted pair and at 3' side
     */
    insertions(fc, pt, mj + 1, k - 1, 0, k - 1, ST_CHG, cb, data);

    /* 1.4.3 base pairs that are now within shifted pair 5' side */
    insertions(fc, pt, mi + 1, old_i - 1, 0, old_i - 1, ST_CHG, cb, data);

    /* 1.4.4 base pairs that are now within shifted pair 3' side */
    insertions(fc, pt, old_j + 1, mj - 1, 0, mj - 1, ST_CHG, cb, data);

    /* 1.4.5 base pairs that are still within the shifted pair */
    insertions(fc,
               pt,
               MAX2(old_i, mi) + 1,
               MIN2(old_j, mj) - 1,
               0,
               MIN2(old_j, mj) - 1,
               ST_CHG,
               cb,
               data);
  }

  /* 2. updates and new deletion moves */
  if (options & VRNA_MOVESET_DELETION) {
    k = old_i;
    if ((mi == k) ||
        (mj == k))
      k = old_j;

    /* the only new deletion move should be the one that deletes
     * the base pair that arises after the shift
     */
    cb(fc, vrna_move_init(-mi, -mj), ST_NEW, data);

    /*
     * handle other deletion moves affected by the shift
     * 1. enclosing pair
     */
    if (enc5 > 0)
      cb(fc, vrna_move_init(-enc5, -enc3), ST_CHG, data);

    /* 2. all base pairs (i, j) with i starting before shifted move */
    deletion_range(fc, pt, enc5 + 1, mi - 1, ST_CHG, cb, data);

    /* 3. all base pairs within the shifted move */
    deletion_range(fc, pt, mi + 1, mj - 1, ST_CHG, cb, data);

    /* 4. all base pairs (i, j) with i starting after shifted move */
    deletion_range(fc, pt, mj + 1, enc3 - 1, ST_CHG, cb, data);
  }

  /* 3. updates and new shift moves */
  if (options & VRNA_MOVESET_SHIFT) {
    k = old_i;
    if ((mi == k) ||
        (mj == k))
      k = old_j;

    /* shift back to the previous configuration */
    if (old_i == mi)
      cb(fc, vrna_move_init(mi, -old_j), ST_NEW, data);
    else if (old_i == mj)
      cb(fc, vrna_move_init(mj, -old_j), ST_NEW, data);
    else if (old_j == mi)
      cb(fc, vrna_move_init(-old_i, mi), ST_NEW, data);
    else
      cb(fc, vrna_move_init(-old_i, mj), ST_NEW, data);

    /* 3.1 All new shift moves that arise due to the novel base pair,
     * i.e. all shifts of the position that stayed constant during the
     * shift
     */
    if ((mi == old_i) || (mi == old_j)) {
      /* mi was the constant part */
      shift_pos(fc, pt, mj, enc5 + 1, mi - 1, ST_NEW, cb, data);
      shift_pos(fc, pt, mj, mi + 1, mj - 1, ST_NEW, cb, data);
      shift_pos(fc, pt, mj, mj + 1, enc3 - 1, ST_NEW, cb, data);
    } else {
      /* mj was the constant part */
      shift_pos(fc, pt, mi, enc5 + 1, mi - 1, ST_NEW, cb, data);
      shift_pos(fc, pt, mi, mi + 1, mj - 1, ST_NEW, cb, data);
      shift_pos(fc, pt, mi, mj + 1, enc3 - 1, ST_NEW, cb, data);
    }

    /*
     * 3.2 all new shift moves that arise from moving one of the pairing partners
     * 3.2.1 shifts of the enclosing pair
     */
    if (enc5 > 0) {
      shift_both(fc, pt, enc5, enc3, k, mi - 1, ST_NEW, cb, data);
      shift_both(fc, pt, enc5, enc3, mj + 1, k, ST_NEW, cb, data);
    }

    /* 3.2.2 shifts of other base pairs 5' of the shift */
    FOR_PAIRED(pt, i, i = enc5 + 1, i < MIN2(k, mi), i++, {
      shift_both(fc, pt, i, pt[i], k, mi - 1, ST_NEW, cb, data);
      shift_both(fc, pt, i, pt[i], mj + 1, k, ST_NEW, cb, data);
    });

    /* 3.2.3 shifts of other base pairs at 3' side of the shift */
    FOR_PAIRED(pt, j, j = MAX2(k, mj) + 1, j < enc3, j++, {
      shift_both(fc, pt, j, pt[j], k, mi - 1, ST_NEW, cb, data);
      shift_both(fc, pt, j, pt[j], mj + 1, k, ST_NEW, cb, data);
    });

    /* 3.2.3 shifts of pairs now within the shifted pair. */
    FOR_PAIRED(pt, i, i = mi + 1, i < old_i, i++, {
      shift_both(fc, pt, i, pt[i], old_i, mj - 1, ST_NEW, cb, data);
    });

    FOR_PAIRED(pt, i, i = old_j + 1, i < mj, i++, {
      shift_both(fc, pt, i, pt[i], mi + 1, old_j, ST_NEW, cb, data);
    });

    if ((mi < k) && (k < mj)) {
      /* 3.2.4 novel base pairs inside the shifted pair */
      FOR_PAIRED(pt, i, i = old_i + 1, i < old_j, i++, {
        shift_both(fc, pt, i, pt[i], mi + 1, old_i, ST_NEW, cb, data);
        shift_both(fc, pt, i, pt[i], old_j, mj - 1, ST_NEW, cb, data);
      });
    } else {
      /* 3.2.5 base pairs that are now located outside the shifted base pair */
      FOR_PAIRED(pt, i, i = old_i + 1, i < mi, i++, {
        shift_both(fc, pt, i, pt[i], enc5 + 1, old_i, ST_NEW, cb, data);
        shift_both(fc, pt, i, pt[i], mj + 1, enc3 - 1, ST_NEW, cb, data);
      });

      FOR_PAIRED(pt, j, j = mj + 1, j < old_j, j++, {
        shift_both(fc, pt, j, pt[j], enc5 + 1, mi - 1, ST_NEW, cb, data);
        shift_both(fc, pt, j, pt[j], old_j, enc3 - 1, ST_NEW, cb, data);
      });
    }

    /* 3.3 now for the shift moves that changed */
    /* this should affect all other shift moves within the same loop(s),
     * ie. those enclosed by the same enclosing pair, and those enclosed
     * by the newly formed base pair
     */

    /* 3.3.1 shifts of the current base pair where the part that stayed
     * constant doesn't move
     */
    if ((mi == old_i) || (mi == old_j)) {
      /* mi was the constant part */
      shift_pos(fc, pt, mi, enc5 + 1, MIN2(old_i, mi) - 1, ST_CHG, cb, data);
      shift_pos(fc, pt, mi, mi + 1, MIN2(old_j, mj) - 1, ST_CHG, cb, data);
      shift_pos(fc, pt, mi, MAX2(old_j, mj) + 1, enc3 - 1, ST_CHG, cb, data);
    } else {
      /* mj was the constant part */
      shift_pos(fc, pt, mj, enc5 + 1, MIN2(old_i, mi) - 1, ST_CHG, cb, data);
      shift_pos(fc, pt, mj, MAX2(old_i, mi) + 1, mj - 1, ST_CHG, cb, data);
      shift_pos(fc, pt, mj, MAX2(old_j, mj) + 1, enc3 - 1, ST_CHG, cb, data);
    }

    /* 3.3.2 shifts of the enclosing base pair (if any) */
    if (enc5 > 0) {
      shift_pos(fc, pt, enc5, enc5 + 1, MIN2(k, mi) - 1, ST_CHG, cb, data);
      shift_pos(fc, pt, enc5, MAX2(k, mj) + 1, enc3 - 1, ST_CHG, cb, data);
      shift_pos(fc, pt, enc3, enc5 + 1, MIN2(k, mi) - 1, ST_CHG, cb, data);
      shift_pos(fc, pt, enc3, MAX2(k, mj) + 1, enc3 - 1, ST_CHG, cb, data);
    }

    /* shifts of all other base pairs */
    FOR_PAIRED(pt, i, i = enc5 + 1, i < MIN2(k, mi), i++, {
      shift_both(fc, pt, i, pt[i], enc5 + 1, i - 1, ST_CHG, cb, data);
      shift_both(fc, pt, i, pt[i], i + 1, pt[i] - 1, ST_CHG, cb, data);
      shift_both(fc, pt, i, pt[i], pt[i] + 1, MIN2(k, mi) - 1, ST_CHG, cb, data);
      shift_both(fc, pt, i, pt[i], MAX2(k, mj) + 1, enc3 - 1, ST_CHG, cb, data);
    });

    FOR_PAIRED(pt, i, i = MAX2(k, mi) + 1, i < MIN2(k, mj), i++, {
      shift_both(fc, pt, i, pt[i], MAX2(k, mi) + 1, i - 1, ST_CHG, cb, data);
      shift_both(fc, pt, i, pt[i], i + 1, pt[i] - 1, ST_CHG, cb, data);
      shift_both(fc, pt, i, pt[i], pt[i] + 1, MIN2(k, mj) - 1, ST_CHG, cb, data);
    });

    FOR_PAIRED(pt, j, j = MAX2(k, mj) + 1, j < enc3, j++, {
      shift_both(fc, pt, j, pt[j], enc5 + 1, MIN2(k, mi) - 1, ST_CHG, cb, data);
      shift_both(fc, pt, j, pt[j], MAX2(k, mj) + 1, j - 1, ST_CHG, cb, data);
      shift_both(fc, pt, j, pt[j], j + 1, pt[j] - 1, ST_CHG, cb, data);
      shift_both(fc, pt, j, pt[j], pt[j] + 1, enc3 - 1, ST_CHG, cb, data);
    });
  }
}


PRIVATE INLINE void
enclosing_pair(const short  *pt,
               int          pair_pos_5,
               int          pair_pos_3,
               int          *enc_pos_5,
               int          *enc_pos_3)
{
  int n = pt[0]; /* length of structure */

  *enc_pos_5  = 0;
  *enc_pos_3  = n + 1;

  for (int i = pair_pos_5 - 1; i > 0; i--) {
    if (pt[i] == 0) {
      continue;
    } else if (pt[i] < i) {
      i = pt[i]; /* hop over branching stems */
    } else if (pt[i] > i) {
      /* found enclosing pair */
      *enc_pos_5  = i;
      *enc_pos_3  = pt[i];
      break;
    }
  }
}


PRIVATE void
generate_conflicts_local_nb_insertion(vrna_fold_compound_t  *fc,
                                      const short           *pt,
                                      const vrna_move_t     *move,
                                      vrna_move_update_f    cb,
                                      void                  *data,
                                      unsigned int          options)
{
  unsigned int  n = fc->length;
  int           i, j, enc5, enc3;
  int           min_loop_size = fc->params->model_details.min_loop_size;
  int           mi            = move->pos_5;
  int           mj            = move->pos_3;

  /* find out actual enclosing pair that delimits the loop affected by the current move */
  enclosing_pair(pt, mi, mj, &enc5, &enc3);

  /*
   *  Determine all previously possible base pairs insertions
   *  that became invalid after inserting (move->pos_5, move->pos_3)
   */
  if (options & VRNA_MOVESET_INSERTION) {
    /*
     *  1. base pairs that start before move->pos_5 and end within
     *  interval [move->pos_5, move->pos_3]
     */

    /* 1.1 remove neighbors that would introduce base pair (i, mi) */
    insertions(fc, pt, enc5 + 1, mi - 1, mi, mi, ST_DEL, cb, data);
    /*
     * 1.2 remove neighbors that would introduce base pair (i, k)
     * with mi < k < mj
     */
    insertions(fc, pt, enc5 + 1, mi - 1, mi + 1, mj - 1, ST_DEL, cb, data);
    /* 1.3 remove neighbors that would introduce base pair (i, mj) */
    insertions(fc, pt, enc5 + 1, mi - 1, mj, mj, ST_DEL, cb, data);

    /*
     *  2. base pairs that start at move->pos_5 and end within interval
     *  [move->pos_5, move->pos_3]
     */
    insertions(fc, pt, mi, mi, 0, mj - 1, ST_DEL, cb, data);

    /*
     *  2.1 neighbors that would introduce base pairs (mi, k)
     *  with mj < k < enc3
     */
    insertions(fc, pt, mi, mi, mj + 1, enc3 - 1, ST_DEL, cb, data);

    /*
     *  3. base pairs that start in (move->pos_5, move->pos_3) and
     *  end after move->pos_3
     */

    /*  3.1 neighbors that introduce base pair (i, mj) */
    insertions(fc, pt, mi + 1, mj - 1, mj, mj, ST_DEL, cb, data);
    /*
     *  3.2 neighbors that introduce base pair (i, k)
     *  with mj < k < enc3
     */
    insertions(fc, pt, mi + 1, mj - 1, mj + 1, enc3 - 1, ST_DEL, cb, data);

    /*
     *  4. base pairs that start at move->pos_3 and end somewhere downstream
     */
    insertions(fc, pt, mj, mj, mj + 1, enc3 - 1, ST_DEL, cb, data);
  }

  if (options & VRNA_MOVESET_DELETION) {
    /* no deletions are affected if the next move is an insertion */
  }

  /* next all invalid shift moves */
  if (options & VRNA_MOVESET_SHIFT) {
    /* 1st, the pairing partners of the enclosing pair (if existing) must not
     * interfere with interval [move->pos_5:move->pos_3]
     */
    if (enc5 > 0)
      shift_both(fc, pt, enc5, enc3, mi, mj, ST_DEL, cb, data);

    /* 2nd, all other base pairs within the same loop and outside
     * the new base pair must not interfere with the interval
     * [move->pos_5:move->pos_3]
     */
    FOR_PAIRED(pt, i, i = enc5 + 1, i < mi, i++, {
      shift_both(fc, pt, i, pt[i], mi, mj, ST_DEL, cb, data);
    });

    FOR_PAIRED(pt, j, j = mj + 1, j < enc3, j++, {
      shift_both(fc, pt, j, pt[j], mi, mj, ST_DEL, cb, data);
    });

    /* Finally, invalidate shift moves from pairs that will be enclosed by
     * (move->pos_5, move->pos_3) outside the interval [move->pos_5,move->pos_3]
     */
    FOR_PAIRED(pt, j, j = mi + 1, j < mj, j++, {
      shift_both(fc, pt, j, pt[j], enc5 + 1, mi, ST_DEL, cb, data);
      shift_both(fc, pt, j, pt[j], mj, enc3 - 1, ST_DEL, cb, data);
    });
  }
}


PRIVATE void
generate_conflicts_local_nb_deletion(vrna_fold_compound_t *fc,
                                     const short          *pt,
                                     const vrna_move_t    *move,
                                     vrna_move_update_f   cb,
                                     void                 *data,
                                     unsigned int         options)
{
  unsigned int  n = fc->length;
  int           i, j, enc5, enc3;
  int           min_loop_size = fc->params->model_details.min_loop_size;
  int           mi            = -move->pos_5;
  int           mj            = -move->pos_3;

  /* upon deletion of a base pair, conflicts can only
   * arise for shift moves of base pair (move->pos_5, move->pos_3)
   */
  if (options & VRNA_MOVESET_SHIFT) {
    /* find out actual enclosing pair that delimits the loop affected by the current move */
    enclosing_pair(pt, mi, mj, &enc5, &enc3);

    /* moves to the region 5' of move->pos_5 */
    shift_both(fc, pt, mi, mj, enc5 + 1, mi - 1, ST_DEL, cb, data);
    shift_both(fc, pt, mi, mj, mi + 1, mj - 1, ST_DEL, cb, data);
    shift_both(fc, pt, mi, mj, mj + 1, enc3 - 1, ST_DEL, cb, data);
  }
}


PRIVATE void
generate_conflicts_local_nb_shift(vrna_fold_compound_t  *fc,
                                  const short           *pt,
                                  const vrna_move_t     *move,
                                  vrna_move_update_f    cb,
                                  void                  *data,
                                  unsigned int          options)
{
  unsigned int  n = fc->length;
  int           i, j, p, q, mi, mj, enc5, enc3;

  /* here, we first determine the base pair (p,q) that is about
   * to change one of its pairing partners
   */
  if (move->pos_5 < 0) {
    p   = pt[move->pos_3];
    q   = move->pos_3;
    mi  = -move->pos_5;
    mj  = move->pos_3;
  } else {
    p   = move->pos_5;
    q   = pt[move->pos_5];
    mi  = move->pos_5;
    mj  = -move->pos_3;
  }

  if (p > q) {
    i = p;
    p = q;
    q = i;
  }

  if (mi > mj) {
    i   = mi;
    mi  = mj;
    mj  = i;
  }

  /* find out actual enclosing pair that delimits the loop affected by the current move */
  enclosing_pair(pt, p, q, &enc5, &enc3);

  /* we now have the current base pair as (p, q) which will change
   * to (mi, mj) or (p', q') after application of the move
   * Here, we distinguish 6 cases, where either p or q stays constant
   * and the pairing partner is shifted towards 3' end of the enclosing
   * loop, the 5' end of the enclosing loop, or towards the segment
   * enclosed by (p,q).
   */
  if (options & VRNA_MOVESET_INSERTION) {
    if (mi == q) {
      /*
       * 1. former pair (p,q) shifts to (q, p') with p < q < p'
       *
       * ......(......)....
       * ......p......q....
       *
       * becomes
       *
       * .............(...)
       * .............q...p'
       */

      /* 1.1. invalidate insertions where one pairing partner involves p' */
      insertions(fc, pt, enc5 + 1, p - 1, mj, mj, ST_DEL, cb, data);
      insertions(fc, pt, q + 1, mj - 1, mj, mj, ST_DEL, cb, data);
      insertions(fc, pt, mj, mj, mj + 1, enc3 - 1, ST_DEL, cb, data);

      /* 1.2 invalidate insertion from the outside into interval [q + 1: mj] */
      insertions(fc, pt, enc5 + 1, p - 1, q + 1, mj - 1, ST_DEL, cb, data);
      insertions(fc, pt, q + 1, mj - 1, mj + 1, enc3 - 1, ST_DEL, cb, data);
    } else if (mj == p) {
      /*
       * 2. q of (p,q) moved to position q' < p forming new pair (q', p)
       *
       * ......(......)....
       * ......p......q....
       *
       * becomes
       *
       * .(....)...........
       * .q'...p...........
       */

      /* 2.1 invalidate insertions where one pairing partner involves q' */
      insertions(fc, pt, enc5 + 1, mi - 1, mi, mi, ST_DEL, cb, data);
      insertions(fc, pt, mi, mi, mi + 1, p - 1, ST_DEL, cb, data);
      insertions(fc, pt, mi, mi, q + 1, enc3 - 1, ST_DEL, cb, data);

      /* 2.2 invalidate insertions where one pairing partner is within interval [q':p-1] */
      insertions(fc, pt, enc5 + 1, mi - 1, mi + 1, p - 1, ST_DEL, cb, data);
      insertions(fc, pt, mi + 1, p - 1, q + 1, enc3 - 1, ST_DEL, cb, data);
    } else if (p - mi > 0) {
      /*
       * 3. former pair (p,q) shifts to (p', q) with p' < p < q
       *
       * ......(......)....
       * ......p......q....
       *
       * becomes
       *
       * ...(.........)....
       * ...p'........q....
       */

      /* 3.1 invalidate insertions where one pairing partner involves p' */
      insertions(fc, pt, enc5 + 1, mi - 1, mi, mi, ST_DEL, cb, data);
      insertions(fc, pt, mi, mi, mi + 1, p - 1, ST_DEL, cb, data);
      insertions(fc, pt, mi, mi, q + 1, enc3 - 1, ST_DEL, cb, data);

      /* 3.2 invalidate insertions that span into free'd interval
       * [p':p]
       */
      insertions(fc, pt, enc5 + 1, mi - 1, mi + 1, p - 1, ST_DEL, cb, data);
      insertions(fc, pt, mi + 1, p - 1, q + 1, enc3 - 1, ST_DEL, cb, data);
    } else if (mj - q > 0) {
      /*
       * 4. former pair (p,q) shifts to (p, q') with p < q < q'
       *
       * ......(......)....
       * ......p......q....
       *
       * becomes
       *
       * ......(..........)
       * ......p..........q'
       */

      /* 4.1. invalidate insertions where one pairing partner involves q' */
      insertions(fc, pt, enc5 + 1, p - 1, mj, mj, ST_DEL, cb, data);
      insertions(fc, pt, q + 1, mj - 1, mj, mj, ST_DEL, cb, data);
      insertions(fc, pt, mj, mj, mj + 1, enc3 - 1, ST_DEL, cb, data);

      /* 4.2 invalidate insertions that span into free'd interval
       * [q:q']
       */
      insertions(fc, pt, enc5 + 1, p - 1, q + 1, mj - 1, ST_DEL, cb, data);
      insertions(fc, pt, q + 1, mj - 1, mj + 1, enc3 - 1, ST_DEL, cb, data);
    } else if (mi - p > 0) {
      /*
       * 5. former pair (p,q) shifts to (p', q) with p < p' < q
       *
       * ......(......)....
       * ......p......q....
       *
       * becomes
       *
       * .........(...)....
       * .........p'..q....
       */

      /* 5.1 invalidate insertions where one pairing partner is p' */
      insertions(fc, pt, p + 1, mi, mi, q - 1, ST_DEL, cb, data);
    } else if (q - mj > 0) {
      /*
       * 6. former pair (p,q) shifts to (p, q') with p < q' < q
       *
       * ......(......)....
       * ......p......q....
       *
       * becomes
       *
       * ......(...).......
       * ......p...q'......
       */

      /* 6.1 invalidate insertions where one pairing partner is q' */
      insertions(fc, pt, p + 1, mj, mj, q - 1, ST_DEL, cb, data);
    }
  }

  if (options & VRNA_MOVESET_DELETION) {
    /* the only deletion move invalidated is the original pair that moves */
    cb(fc, vrna_move_init(-p, -q), ST_DEL, data);
  }

  if (options & VRNA_MOVESET_SHIFT) {
    if (p == mj) {
      /*
       * 1. q of (p,q) moved to position q' < p forming new pair (q', p)
       *
       * ......(......)....
       * ......p......q....
       *
       * becomes
       *
       * .(....)...........
       * .q'...p...........
       */

      /* 1.1 invalidate all shifts of p for pair (p, q) */
      shift_pos(fc, pt, q, enc5 + 1, p - 1, ST_DEL, cb, data);
      shift_pos(fc, pt, q, p + 1, q - 1, ST_DEL, cb, data);
      shift_pos(fc, pt, q, q + 1, enc3 - 1, ST_DEL, cb, data);

      /* 1.2 invalidate moves of enclosing pair positions, if any, into interval [mi:p-1] */
      if (enc5 > 0)
        shift_both(fc, pt, enc5, enc3, mi, p - 1, ST_DEL, cb, data);

      /* 1.3 invalidate moves of other pairs into interval [mi:p-1] */
      FOR_PAIRED(pt, i, i = enc5 + 1, i < mi, i++, {
        shift_both(fc, pt, i, pt[i], mi, p - 1, ST_DEL, cb, data);
      });

      FOR_PAIRED(pt, i, i = q + 1, i < enc3, i++, {
        shift_both(fc, pt, i, pt[i], mi, p - 1, ST_DEL, cb, data);
      });

      /* 1.4 invalidate moves from within the interval [mi+1:p-1] outside of the interval */
      FOR_PAIRED(pt, j, j = mi + 1, j < p, j++, {
        shift_both(fc, pt, j, pt[j], enc5 + 1, mi, ST_DEL, cb, data);
        shift_both(fc, pt, j, pt[j], q + 1, enc3 - 1, ST_DEL, cb, data);
      });
    } else if (q == mi) {
      /*
       * 2. former pair (p,q) shifts to (q, p') with p < q < p'
       *
       * ......(......)....
       * ......p......q....
       *
       * becomes
       *
       * .............(...)
       * .............q...p'
       */

      /* 2.1 invalidate all shifts of q for pair (p,q) */
      shift_pos(fc, pt, p, enc5 + 1, p - 1, ST_DEL, cb, data);
      shift_pos(fc, pt, p, p + 1, q - 1, ST_DEL, cb, data);
      shift_pos(fc, pt, p, q + 1, enc3 - 1, ST_DEL, cb, data);

      /* 2.2 invalidate moves of enclosing pair, if any, into interval [q+1:mj] */
      if (enc5 > 0)
        shift_both(fc, pt, enc5, enc3, q + 1, mj, ST_DEL, cb, data);

      /* 2.3 invalidate moves of other pairs into interval [q+1:mj] */
      FOR_PAIRED(pt, i, i = enc5 + 1, i < p, i++, {
        shift_both(fc, pt, i, pt[i], q + 1, mj, ST_DEL, cb, data);
      });

      FOR_PAIRED(pt, i, i = mj + 1, i < enc3, i++, {
        shift_both(fc, pt, i, pt[i], q + 1, mj, ST_DEL, cb, data);
      });

      /* 2.4 invalidate move of pairs within interval [q+1:mj-1] outside of the interval */
      FOR_PAIRED(pt, j, j = q + 1, j < mj, j++, {
        shift_both(fc, pt, j, pt[j], enc5 + 1, p - 1, ST_DEL, cb, data);
        shift_both(fc, pt, j, pt[j], mj, enc3 - 1, ST_DEL, cb, data);
      });
    } else if (p - mi > 0) {
      /*
       * 3. former pair (p,q) shifts to (p', q) with p' < p < q
       *
       * ......(......)....
       * ......p......q....
       *
       * becomes
       *
       * ...(.........)....
       * ...p'........q....
       */

      /* 3.1 invalidate any other shift of q for pair (p, q) */
      shift_pos(fc, pt, p, enc5 + 1, p - 1, ST_DEL, cb, data);
      shift_pos(fc, pt, p, p + 1, q - 1, ST_DEL, cb, data);
      shift_pos(fc, pt, p, q + 1, enc3 - 1, ST_DEL, cb, data);

      /* 3.2 invalidate moves of enclosing pair, if any, into interval [mi:p-1] */
      if (enc5 > 0)
        shift_both(fc, pt, enc5, enc3, mi, p - 1, ST_DEL, cb, data);

      /* 3.3 invalidate shift of any other pair into the interval [mi:p-1] */
      FOR_PAIRED(pt, i, i = enc5 + 1, i < mi, i++, {
        shift_both(fc, pt, i, pt[i], mi, p - 1, ST_DEL, cb, data);
      });

      FOR_PAIRED(pt, i, i = q + 1, i < enc3, i++, {
        shift_both(fc, pt, i, pt[i], mi, p - 1, ST_DEL, cb, data);
      });

      /* 3.4 invalidate shift moves from pairs within interval to the outside */
      FOR_PAIRED(pt, j, j = mi + 1, j < p, j++, {
        shift_both(fc, pt, j, pt[j], enc5 + 1, mi, ST_DEL, cb, data);
        shift_both(fc, pt, j, pt[j], q + 1, enc3 - 1, ST_DEL, cb, data);
      });
    } else if (mj - q > 0) {
      /*
       * 4. former pair (p,q) shifts to (p, q') with p < q < q'
       *
       * ......(......)....
       * ......p......q....
       *
       * becomes
       *
       * ......(..........)
       * ......p..........q'
       */

      /* 4.1 invalidate any shift of p for pair (p, q) */
      shift_pos(fc, pt, q, enc5 + 1, p - 1, ST_DEL, cb, data);
      shift_pos(fc, pt, q, p + 1, q - 1, ST_DEL, cb, data);
      shift_pos(fc, pt, q, q + 1, enc3 - 1, ST_DEL, cb, data);

      /* 4.2 invalidate moves of enclosing pair, if any, into interval [q+1:mj] */
      if (enc5 > 0)
        shift_both(fc, pt, enc5, enc3, q + 1, mj, ST_DEL, cb, data);

      /* 4.3 invalidate shifts from other pairs into the interval [q+1: mj] */
      FOR_PAIRED(pt, i, i = enc5 + 1, i < p, i++, {
        shift_both(fc, pt, i, pt[i], q + 1, mj, ST_DEL, cb, data);
      });

      FOR_PAIRED(pt, i, i = mj + 1, i < enc3, i++, {
        shift_both(fc, pt, i, pt[i], q + 1, mj, ST_DEL, cb, data);
      });

      /* 4.4 invalidate shifts from within the interval [q + 1: mj - 1] to the outside */
      FOR_PAIRED(pt, j, j = q + 1, j < mj, j++, {
        shift_both(fc, pt, j, pt[j], enc5 + 1, p - 1, ST_DEL, cb, data);
        shift_both(fc, pt, j, pt[j], mj, enc3 - 1, ST_DEL, cb, data);
      });
    } else if (mi - p > 0) {
      /*
       * 5. former pair (p,q) shifts to (p', q) with p < p' < q
       *
       * ......(......)....
       * ......p......q....
       *
       * becomes
       *
       * .........(...)....
       * .........p'..q....
       */

      /* 5.1 invalidate shifts of q for pair (p, q) */
      shift_pos(fc, pt, p, enc5 + 1, p - 1, ST_DEL, cb, data);
      shift_pos(fc, pt, p, p + 1, q - 1, ST_DEL, cb, data);
      shift_pos(fc, pt, p, q + 1, enc3 - 1, ST_DEL, cb, data);

      /* 5.2 invalidate moves of pairs into interval [p + 1:mi] */
      FOR_PAIRED(pt, j, j = mi + 1, j < q, j++, {
        shift_both(fc, pt, j, pt[j], p + 1, mi, ST_DEL, cb, data);
      });

      /* 5.3 invalidate moves of pairs from within interval [p + 1: mi] */
      FOR_PAIRED(pt, i, i = p + 1, i < mi, i++, {
        shift_both(fc, pt, i, pt[i], mi, q - 1, ST_DEL, cb, data);
      });
    } else if (q - mj > 0) {
      /*
       * 6. former pair (p,q) shifts to (p, q') with p < q' < q
       *
       * ......(......)....
       * ......p......q....
       *
       * becomes
       *
       * ......(...).......
       * ......p...q'......
       */

      /* 6.1 invalidate any shifts of p for pair (p, q) */
      shift_pos(fc, pt, q, enc5 + 1, p - 1, ST_DEL, cb, data);
      shift_pos(fc, pt, q, p + 1, q - 1, ST_DEL, cb, data);
      shift_pos(fc, pt, q, q + 1, enc3 - 1, ST_DEL, cb, data);

      /* 6.2 invalidate shifts of pairs into interval [mj:q-1] */
      FOR_PAIRED(pt, i, i = p + 1, i < mj, i++, {
        shift_both(fc, pt, i, pt[i], mj, q - 1, ST_DEL, cb, data);
      });

      /* 6.3 invalidate shifts from within interval [mj+1:q-1] */
      FOR_PAIRED(pt, j, j = mj + 1, j < q, j++, {
        shift_both(fc, pt, j, pt[j], p + 1, mj, ST_DEL, cb, data);
      });
    }
  }
}
