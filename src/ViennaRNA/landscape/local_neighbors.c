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
    }}

#define FOR_PAIRED(ptable, idx, start, cond, next, code)      {\
    for ((start); (cond); (next)) { \
      if (ptable[idx] > idx) { \
        { code }; \
        idx = ptable[idx]; \
      } \
    }}


PRIVATE INLINE void
enclosing_pair(const short *pt,
               int pair_pos_5,
               int pair_pos_3,
               int *enc_pos_5,
               int *enc_pos_3);


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


PRIVATE INLINE void
shift_range5_cb(vrna_fold_compound_t  *fc,
                const short           *pt,
                unsigned int          i,
                unsigned int          j,
                unsigned int          min_k,
                unsigned int          max_k,
                unsigned int          status,
                vrna_move_update_f    cb,
                void                  *data);


PRIVATE INLINE void
shift_range3_cb(vrna_fold_compound_t  *fc,
                const short           *pt,
                unsigned int          i,
                unsigned int          j,
                unsigned int          min_k,
                unsigned int          max_k,
                unsigned int          status,
                vrna_move_update_f    cb,
                void                  *data);


PRIVATE INLINE void
shift_rangeenc_cb(vrna_fold_compound_t  *fc,
                  const short           *pt,
                  unsigned int          i,
                  unsigned int          j,
                  unsigned int          min_k,
                  unsigned int          max_k,
                  unsigned int          status,
                  vrna_move_update_f    cb,
                  void                  *data);


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
    return vc->params->model_details.pair[vc->sequence_encoding2[i]][vc->sequence_encoding2[j]] != 0; /* see pair_mat.h */

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
      cb(fc, move, VRNA_NEIGHBOR_INVALID, data);

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
          ptable[-move.pos_5] = 0;
          ptable[move.pos_3] = affected_pair.pos_5;
          ptable[affected_pair.pos_5] = move.pos_3;
        } else if ((affected_pair.pos_5 > 0) &&
                   (affected_pair.pos_3 < 0)) {
          /* shift move, 3' position changes, 5' position stays constant */
          ptable[-move.pos_3] = 0;
          ptable[move.pos_5] = affected_pair.pos_3;
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
insertions_range_cb(vrna_fold_compound_t  *fc,
                    const short           *pt,
                    int                   i,
                    int                   first_j,
                    int                   last_j,
                    unsigned int          status,
                    vrna_move_update_f    cb,
                    void                  *data)
{
  int j;

  if (first_j > last_j)
    return;

  /* skip everything before first_j */
  for (j = i + 1; j < first_j; j++)
    if (pt[j] > j)
      j = pt[j];

  if (j > last_j)
    return;

  for (; j <= last_j; j++) {
    if (pt[j] > j)
      j = pt[j]; /* hop over branching stems */
    else if ((pt[j] == 0) && (is_compatible(fc, i, j)))
      cb(fc, vrna_move_init(i, j), status, data);
  }
}


PRIVATE INLINE void
deletion_range(vrna_fold_compound_t *fc,
               const short          *pt,
               unsigned int         start,
               unsigned int         stop,
               unsigned int         status,
                vrna_move_update_f  cb,
                void                *data)
{
  unsigned int i;

  FOR_PAIRED(pt, i, i = start, i <= stop, i++, {
    cb(fc, vrna_move_init(-i, -pt[i]), status, data);
  });
}


PRIVATE INLINE void
shift_range5_cb(vrna_fold_compound_t  *fc,
                const short           *pt,
                unsigned int          i,
                unsigned int          j,
                unsigned int          min_k,
                unsigned int          max_k,
                unsigned int          status,
                vrna_move_update_f    cb,
                void                  *data)
{
  unsigned int k;

  FOR_UNPAIRED(pt, k, k = min_k, k <= max_k, k++, {
    if (is_compatible(fc, k, i))
      cb(fc, vrna_move_init(-k, i), status, data);

    if (is_compatible(fc, k, j))
      cb(fc, vrna_move_init(-k, j), status, data);
  });
}


PRIVATE INLINE void
shift_range3_cb(vrna_fold_compound_t  *fc,
                const short           *pt,
                unsigned int          i,
                unsigned int          j,
                unsigned int          min_k,
                unsigned int          max_k,
                unsigned int          status,
                vrna_move_update_f    cb,
                void                  *data)
{
  unsigned int k;

  FOR_UNPAIRED(pt, k, k = min_k, k <= max_k, k++, {
    if (is_compatible(fc, i, k))
      cb(fc, vrna_move_init(i, -k), status, data);

    if (is_compatible(fc, j, k))
      cb(fc, vrna_move_init(j, -k), status, data);
  });
}


PRIVATE INLINE void
shift_rangeenc_cb(vrna_fold_compound_t  *fc,
                  const short           *pt,
                  unsigned int          i,
                  unsigned int          j,
                  unsigned int          min_k,
                  unsigned int          max_k,
                  unsigned int          status,
                  vrna_move_update_f    cb,
                  void                  *data)
{
  unsigned int k;

  FOR_UNPAIRED(pt, k, k = min_k, k <= max_k, k++, {
    if (is_compatible(fc, i, k))
      cb(fc, vrna_move_init(i, -k), status, data);

    if (is_compatible(fc, k, j))
      cb(fc, vrna_move_init(-k, j), status, data);
  });
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
  int           i, j, k, move_i, move_j, enclosing_5, enclosing_3;
  int           min_loop_size = fc->params->model_details.min_loop_size;

  move_i  = affected_pair->pos_5;
  move_j  = affected_pair->pos_3;

  /* find out actual enclosing pair that delimits the loop affected by the current move */
  enclosing_pair(pt, move_i, move_j, &enclosing_5, &enclosing_3);

  /* 1. valid base pair deletions */
  if (options & VRNA_MOVESET_DELETION) {
    /* 1.1 removal of enclosing pair */
    if (enclosing_5 > 0)
      cb(fc, vrna_move_init(-enclosing_5, -enclosing_3), VRNA_NEIGHBOR_CHANGE, data);

    /* 1.2 removal of base pair that has just been inserted */
    cb(fc, vrna_move_init(-move_i, -move_j), VRNA_NEIGHBOR_NEW, data);

    /* 1.3 removal of other base pairs that branch-off from the current loop */
    deletion_range(fc, pt, enclosing_5 + 1, move_i - 1, VRNA_NEIGHBOR_CHANGE, cb, data);
    deletion_range(fc, pt, move_i + 1, move_j - 1, VRNA_NEIGHBOR_CHANGE, cb, data);
    deletion_range(fc, pt, move_j + 1, enclosing_3 - 1, VRNA_NEIGHBOR_CHANGE, cb, data);
  }

  /* 2. valid base pair insertions */
  if (options & VRNA_MOVESET_INSERTION) {
    /* 5' side of newly inserted base pair */
    FOR_UNPAIRED(pt, i, i = enclosing_5 + 1, i < move_i, i++, {
      /*
       * determine all potential pairing partners for current nucleotide within
       * 5' side of newly inserted base pair
       */
      insertions_range_cb(fc,
                          pt,
                          i,
                          i + 1,
                          move_i - 1,
                          VRNA_NEIGHBOR_CHANGE,
                          cb,
                          data);

      /*
       * determine all potential pairing partners for current nucleotide within
       * 3' side of newly inserted base pair
       */
      insertions_range_cb(fc,
                          pt,
                          i,
                          move_j + 1,
                          enclosing_3 - 1,
                          VRNA_NEIGHBOR_CHANGE,
                          cb,
                          data);
    });

    /* within newly inserted base pair */
    FOR_UNPAIRED(pt, i, i = move_i + 1, i < move_j, i++, {
      insertions_range_cb(fc,
                          pt,
                          i,
                          i + 1,
                          move_j - 1,
                          VRNA_NEIGHBOR_CHANGE,
                          cb,
                          data);
    });

    /* 3' side of newly inserted base pair */
    FOR_UNPAIRED(pt, i, i = move_j + 1, i < enclosing_3, i++, {
      /*
       * determine all potential pairing partners for current nucleotide within
       * 5' side of newly inserted base pair
       */
      insertions_range_cb(fc,
                          pt,
                          i,
                          i + 1,
                          enclosing_3 - 1,
                          VRNA_NEIGHBOR_CHANGE,
                          cb,
                          data);
    });
  }

  if (options & VRNA_MOVESET_SHIFT) {
    /* 1. novel base pair shift moves due to inserted pair */

    /* 1.1 shift either move_i or move_j towards 5' side
     * interval [enclosing_5 + 1:move_i-1]
     */
    shift_range5_cb(fc,
                    pt,
                    (unsigned int)move_i,
                    (unsigned int)move_j,
                    (unsigned int)enclosing_5 + 1,
                    (unsigned int)move_i - 1,
                    VRNA_NEIGHBOR_NEW,
                    cb,
                    data);

    /* 1.2 shift either move_i or move_j into the
     * enclosed interval [move_i+1:move_j-1]
     */
    shift_rangeenc_cb(fc,
                 pt,
                 (unsigned int)move_i,
                 (unsigned int)move_j,
                 (unsigned int)move_i + 1,
                 (unsigned int)move_j - 1,
                 VRNA_NEIGHBOR_NEW,
                 cb,
                 data);

    /* 1.3 shift either move_i or move_j towards the 3' side
     * interval [move_j+1:enclosing_3 - 1]
     */
    shift_range3_cb(fc,
                    pt,
                    (unsigned int)move_i,
                    (unsigned int)move_j,
                    (unsigned int)move_j + 1,
                    (unsigned int)enclosing_3 - 1,
                    VRNA_NEIGHBOR_NEW,
                    cb,
                    data);

    /* 2. shift moves that need to be re-evaluated */

    /* 2.1 Shifts of the enclosing pair, if any */
    if (enclosing_5 > 0) {
      shift_rangeenc_cb(fc,
                        pt,
                        (unsigned int)enclosing_5,
                        (unsigned int)enclosing_3,
                        (unsigned int)enclosing_5 + 1,
                        (unsigned int)move_i - 1,
                        VRNA_NEIGHBOR_CHANGE,
                        cb,
                        data);

      shift_rangeenc_cb(fc,
                        pt,
                        (unsigned int)enclosing_5,
                        (unsigned int)enclosing_3,
                        (unsigned int)move_j + 1,
                        (unsigned int)enclosing_3 - 1,
                        VRNA_NEIGHBOR_CHANGE,
                        cb,
                        data);
    }

    /* 2.2 shift of base pairs in same loop outside of pair inserted by move (5' side) */
    FOR_PAIRED(pt, i, i = enclosing_5 + 1, i < move_i, i++, {
      j = pt[i];
      /* 2.2.1 shifts within the same loop (5' side) */
      shift_range5_cb(fc,
                      pt,
                      (unsigned int)i,
                      (unsigned int)j,
                      (unsigned int)enclosing_5 + 1,
                      (unsigned int)i - 1,
                      VRNA_NEIGHBOR_CHANGE,
                      cb,
                      data);

      /* 2.2.2 shifts within the loop enclosed by (i, j) */
      shift_rangeenc_cb(fc,
                   pt,
                   (unsigned int)i,
                   (unsigned int)j,
                   (unsigned int)i + 1,
                   (unsigned int)j - 1,
                   VRNA_NEIGHBOR_CHANGE,
                   cb,
                   data);

      /* 2.2.3 shifts within the same loop (3' side) */
      shift_range3_cb(fc,
                      pt,
                      (unsigned int)i,
                      (unsigned int)j,
                      (unsigned int)j + 1,
                      (unsigned int)move_i - 1,
                      VRNA_NEIGHBOR_CHANGE,
                      cb,
                      data);

      shift_range3_cb(fc,
                      pt,
                      (unsigned int)i,
                      (unsigned int)j,
                      (unsigned int)move_j + 1,
                      (unsigned int)enclosing_3 - 1,
                      VRNA_NEIGHBOR_CHANGE,
                      cb,
                      data);
    });

    /* 2.3 shift of base pairs in same loop outside of pair inserted by move (3' side) */
    FOR_PAIRED(pt, i, i = move_j + 1, i < enclosing_3, i++, {
      j = pt[i];
      /* 2.3.1 shifts within the same loop (5' side) */
      shift_range5_cb(fc,
                      pt,
                      (unsigned int)i,
                      (unsigned int)j,
                      (unsigned int)enclosing_5 + 1,
                      (unsigned int)move_i - 1,
                      VRNA_NEIGHBOR_CHANGE,
                      cb,
                      data);
      shift_range5_cb(fc,
                      pt,
                      (unsigned int)i,
                      (unsigned int)j,
                      (unsigned int)move_j + 1,
                      (unsigned int)i - 1,
                      VRNA_NEIGHBOR_CHANGE,
                      cb,
                      data);

      /* 2.3.2 shifts within the loop enclosed by (i, j) */
      shift_rangeenc_cb(fc,
                   pt,
                   (unsigned int)i,
                   (unsigned int)j,
                   (unsigned int)i + 1,
                   (unsigned int)j - 1,
                   VRNA_NEIGHBOR_CHANGE,
                   cb,
                   data);

      /* 2.3.3 shifts within the same loop (3' side) */
      shift_range3_cb(fc,
                      pt,
                      (unsigned int)i,
                      (unsigned int)j,
                      (unsigned int)j + 1,
                      (unsigned int)enclosing_3 - 1,
                      VRNA_NEIGHBOR_CHANGE,
                      cb,
                      data);
    });

    /* 2.4 shifts of base pairs enclosed by novel pair introduced by move */
    FOR_PAIRED(pt, i, i = move_i + 1, i < move_j, i++, {
      j = pt[i];

      /* 2.4.1 shifts within the loop enclosed by novel pair (5' side) */
      shift_range5_cb(fc,
                      pt,
                      (unsigned int)i,
                      (unsigned int)j,
                      (unsigned int)move_i + 1,
                      (unsigned int)i - 1,
                      VRNA_NEIGHBOR_CHANGE,
                      cb,
                      data);

      /* 2.4.2 shifts within the loop enclosed by (i,j) */
      shift_rangeenc_cb(fc,
                   pt,
                   (unsigned int)i,
                   (unsigned int)j,
                   (unsigned int)i + 1,
                   (unsigned int)j - 1,
                   VRNA_NEIGHBOR_CHANGE,
                   cb,
                   data);

      /* 2.4.3 shifts within the loop enclosed by novel pair (3' side) */
      shift_range3_cb(fc,
                      pt,
                      (unsigned int)i,
                      (unsigned int)j,
                      (unsigned int)j + 1,
                      (unsigned int)move_j - 1,
                      VRNA_NEIGHBOR_CHANGE,
                      cb,
                      data);
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
  int           i, j, k, enclosing_5, enclosing_3;
  int           min_loop_size = fc->params->model_details.min_loop_size;
  int           move_i        = affected_pair->pos_5;
  int           move_j        = affected_pair->pos_3;

  /* find out actual enclosing pair that delimits the loop affected by the current move */
  enclosing_pair(pt, move_i, move_j, &enclosing_5, &enclosing_3);

  /* 1. valid base pair deletions */
  if (options & VRNA_MOVESET_DELETION) {
    /* 1.1 removal of enclosing pair */
    if (enclosing_5 > 0)
      cb(fc, vrna_move_init(-enclosing_5, -enclosing_3), VRNA_NEIGHBOR_CHANGE, data);

    /* 1.2 branching base pairs 5' to base pair deleted by current move */
    deletion_range(fc, pt, enclosing_5 + 1, move_i - 1, VRNA_NEIGHBOR_CHANGE, cb, data);

    /* 1.3 branching base pairs enclosed by base pair deleted by current move */
    deletion_range(fc, pt, move_i + 1, move_j - 1, VRNA_NEIGHBOR_CHANGE, cb, data);

    /* 1.4 branching base pairs on 3' side of base pair deleted by current move */
    deletion_range(fc, pt, move_j + 1, enclosing_3 - 1, VRNA_NEIGHBOR_CHANGE, cb, data);
  }

  /* 2. valid base pair insertions */
  if (options & VRNA_MOVESET_INSERTION) {
    /* 2.1 re-insertion of base pair removed by current move */
    cb(fc, vrna_move_init(move_i, move_j), VRNA_NEIGHBOR_NEW, data);

    /* 2.2 insertion of novel pairs that start on 5' side of current move */
    FOR_UNPAIRED(pt, i, i = enclosing_5 + 1, i < move_i, i++, {
      /* 2.2.1 base pairs that end before 5' side of current move */
      insertions_range_cb(fc,
                          pt,
                          i,
                          i + 1,
                          move_i - 1,
                          VRNA_NEIGHBOR_CHANGE,
                          cb,
                          data);

      /* 2.2.2 base pairs that end within interval spanned by the base
       * pair we just deleted
       */
      insertions_range_cb(fc,
                          pt,
                          i,
                          move_i,
                          move_j,
                          VRNA_NEIGHBOR_NEW,
                          cb,
                          data);

      /* 2.2.3 base pairs that end after 3' nucleotide of current move */
      insertions_range_cb(fc,
                          pt,
                          i,
                          move_j + 1,
                          enclosing_3 - 1,
                          VRNA_NEIGHBOR_CHANGE,
                          cb,
                          data);
    });

    /* 2.3 insertion of novel pairs that start at 5' nucleotide of current position */

    /* 2.3.1 ending within loop enclosed by base pair deleted by current move */
    insertions_range_cb(fc,
                        pt,
                        move_i,
                        move_i + 1,
                        move_j - 1,
                        VRNA_NEIGHBOR_NEW,
                        cb,
                        data);

    /* 2.3.2 ending after loop enclosed by base pair deleted by current move */
    insertions_range_cb(fc,
                        pt,
                        move_i,
                        move_j + 1,
                        enclosing_3 - 1,
                        VRNA_NEIGHBOR_NEW,
                        cb,
                        data);

    /* 2.4 insertion of novel pairs that start within loop closed by current move */
    FOR_UNPAIRED(pt, i, i = move_i + 1, i < move_j, i++, {
      /* 2.4.1 */
      insertions_range_cb(fc,
                          pt,
                          i,
                          i + 1,
                          move_j - 1,
                          VRNA_NEIGHBOR_CHANGE,
                          cb,
                          data);

      /* 2.4.2 */
      insertions_range_cb(fc,
                          pt,
                          i,
                          move_j,
                          enclosing_3 - 1,
                          VRNA_NEIGHBOR_NEW,
                          cb,
                          data);
    });

    /* 2.5 insertion of novel pairs that start at 3' nucleotide of current move */
    insertions_range_cb(fc,
                        pt,
                        move_j,
                        move_j + 1,
                        enclosing_3 - 1,
                        VRNA_NEIGHBOR_NEW,
                        cb,
                        data);

    /* 2.6 insertion of novel pairs that start after 3' nucleotide of current move */
    FOR_UNPAIRED(pt, i, i = move_j + 1, i < enclosing_3, i++, {
      insertions_range_cb(fc,
                          pt,
                          i,
                          i + 1,
                          enclosing_3 - 1,
                          VRNA_NEIGHBOR_CHANGE,
                          cb,
                          data);
    });
  }

  /* 3. valid shift moves */
  if (options & VRNA_MOVESET_SHIFT) {
    /* 3.1 shifts of the enclosing pair, if any */
    if (enclosing_5 > 0) {
      shift_rangeenc_cb(fc,
                        pt,
                        (unsigned int)enclosing_5,
                        (unsigned int)enclosing_3,
                        (unsigned int)enclosing_5 + 1,
                        (unsigned int)move_i - 1,
                        VRNA_NEIGHBOR_CHANGE,
                        cb,
                        data);

      shift_rangeenc_cb(fc,
                        pt,
                        (unsigned int)enclosing_5,
                        (unsigned int)enclosing_3,
                        (unsigned int)move_i,
                        (unsigned int)move_j,
                        VRNA_NEIGHBOR_NEW,
                        cb,
                        data);

      shift_rangeenc_cb(fc,
                        pt,
                        (unsigned int)enclosing_5,
                        (unsigned int)enclosing_3,
                        (unsigned int)move_j + 1,
                        (unsigned int)enclosing_3 - 1,
                        VRNA_NEIGHBOR_CHANGE,
                        cb,
                        data);
    }

    /* 3.2 shifts of base pair previously outside of removed pair 5' side */
    FOR_PAIRED(pt, i, i = enclosing_5 + 1, i < move_i, i++, {
      j = pt[i];

      shift_range5_cb(fc,
                      pt,
                      (unsigned int)i,
                      (unsigned int)j,
                      (unsigned int)enclosing_5 + 1,
                      (unsigned int)i - 1,
                      VRNA_NEIGHBOR_CHANGE,
                      cb,
                      data);

      shift_rangeenc_cb(fc,
                   pt,
                   (unsigned int)i,
                   (unsigned int)j,
                   (unsigned int)i + 1,
                   (unsigned int)j - 1,
                   VRNA_NEIGHBOR_CHANGE,
                   cb,
                   data);

      shift_range3_cb(fc,
                      pt,
                      (unsigned int)i,
                      (unsigned int)j,
                      (unsigned int)j + 1,
                      (unsigned int)move_i - 1,
                      VRNA_NEIGHBOR_CHANGE,
                      cb,
                      data);

      shift_rangeenc_cb(fc,
                        pt,
                        (unsigned int)i,
                        (unsigned int)j,
                        (unsigned int)move_i,
                        (unsigned int)move_j,
                        VRNA_NEIGHBOR_NEW,
                        cb,
                        data);

      shift_range3_cb(fc,
                      pt,
                      (unsigned int)i,
                      (unsigned int)j,
                      (unsigned int)move_j + 1,
                      (unsigned int)enclosing_3 - 1,
                      VRNA_NEIGHBOR_CHANGE,
                      cb,
                      data);
    });

    /* 3.2 shifts of base pair previously outside of removed pair 3' side */
    FOR_PAIRED(pt, i, i = move_j + 1, i < enclosing_3, i++, {
      j = pt[i];

      shift_range5_cb(fc,
                      pt,
                      (unsigned int)i,
                      (unsigned int)j,
                      (unsigned int)enclosing_5 + 1,
                      (unsigned int)move_i - 1,
                      VRNA_NEIGHBOR_CHANGE,
                      cb,
                      data);

      shift_rangeenc_cb(fc,
                        pt,
                        (unsigned int)i,
                        (unsigned int)j,
                        (unsigned int)move_i,
                        (unsigned int)move_j,
                        VRNA_NEIGHBOR_NEW,
                        cb,
                        data);

      shift_range5_cb(fc,
                      pt,
                      (unsigned int)i,
                      (unsigned int)j,
                      (unsigned int)move_j + 1,
                      (unsigned int)i - 1,
                      VRNA_NEIGHBOR_CHANGE,
                      cb,
                      data);

      shift_rangeenc_cb(fc,
                   pt,
                   (unsigned int)i,
                   (unsigned int)j,
                   (unsigned int)i + 1,
                   (unsigned int)j - 1,
                   VRNA_NEIGHBOR_CHANGE,
                   cb,
                   data);

      shift_range3_cb(fc,
                      pt,
                      (unsigned int)i,
                      (unsigned int)j,
                      (unsigned int)j + 1,
                      (unsigned int)enclosing_3 - 1,
                      VRNA_NEIGHBOR_CHANGE,
                      cb,
                      data);
    });

    /* 3.3. shifts of base pairs previously enclosed by removed pair */
    FOR_PAIRED(pt, i, i = move_i + 1, i < move_j, i++, {
      j = pt[i];
      /* 3.3.1 shifts outside towards 5' side */
      shift_rangeenc_cb(fc,
                        pt,
                        (unsigned int)i,
                        (unsigned int)j,
                        (unsigned int)enclosing_5 + 1,
                        (unsigned int)move_i,
                        VRNA_NEIGHBOR_NEW,
                        cb,
                        data);

      shift_range5_cb(fc,
                      pt,
                      (unsigned int)i,
                      (unsigned int)j,
                      (unsigned int)move_i + 1,
                      (unsigned int)i - 1,
                      VRNA_NEIGHBOR_CHANGE,
                      cb,
                      data);

      /* 3.3.2 shifts inside (i, j) */
      shift_rangeenc_cb(fc,
                   pt,
                   (unsigned int)i,
                   (unsigned int)j,
                   (unsigned int)i + 1,
                   (unsigned int)j - 1,
                   VRNA_NEIGHBOR_CHANGE,
                   cb,
                   data);

      /* 3.3.3 shifts outside towards 3' side */
      shift_range3_cb(fc,
                      pt,
                      (unsigned int)i,
                      (unsigned int)j,
                      (unsigned int)j + 1,
                      (unsigned int)move_j - 1,
                      VRNA_NEIGHBOR_CHANGE,
                      cb,
                      data);

      shift_rangeenc_cb(fc,
                        pt,
                        (unsigned int)i,
                        (unsigned int)j,
                        (unsigned int)move_j,
                        (unsigned int)enclosing_3 - 1,
                        VRNA_NEIGHBOR_NEW,
                        cb,
                        data);
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
  int i, j, k, move_i, move_j, old_i, old_j, enclosing_5, enclosing_3;

  enclosing_5 = 0;
  enclosing_3 = fc->length + 1;
  old_i       = affected_pair->pos_5;
  old_j       = affected_pair->pos_3;

  if (move->pos_5 < 0) {
    move_i  = -move->pos_5;
    move_j  = move->pos_3;
  } else {
    move_i  = move->pos_5;
    move_j  = -move->pos_3;
  }

  if (move_j < move_i) {
    k       = move_i;
    move_i  = move_j;
    move_j  = k;
  }

  /* find out actual enclosing pair that delimits the loop affected by the current move */
  enclosing_pair(pt, move_i, move_j, &enclosing_5, &enclosing_3);

  /* 1. updates of insertion moves */
  if (options & VRNA_MOVESET_INSERTION) {
    /* 1.1 new insertions using the now vacant segment */
    k = old_i;
    if ((move_i == k) ||
        (move_j == k))
      k = old_j;

    /* base pairs still outside of (move_i, move_j) at 5' side */
    FOR_UNPAIRED(pt, i, i = k, i < move_i, i++, {
      FOR_UNPAIRED(pt, j, j = enclosing_5 + 1, j < i, j++, {
        if (is_compatible(fc, j, i))
          cb(fc, vrna_move_init(j, i), VRNA_NEIGHBOR_NEW, data);
      });

      FOR_UNPAIRED(pt, j, j = move_j + 1, j < enclosing_3, j++, {
        if (is_compatible(fc, i, j))
          cb(fc, vrna_move_init(i, j), VRNA_NEIGHBOR_NEW, data);
      });
    });

    /* base pairs still outside of (move_i, move_j) at 3' side */
    FOR_UNPAIRED(pt, j, j = move_j + 1, j <= k, j++, {
      FOR_UNPAIRED(pt, i, i = enclosing_5 + 1, i < move_i, i++, {
        if (is_compatible(fc, i, j))
          cb(fc, vrna_move_init(i, j), VRNA_NEIGHBOR_NEW, data);
      });

      FOR_UNPAIRED(pt, i, i = k + 1, i < enclosing_3, i++, {
        if (is_compatible(fc, j, i))
          cb(fc, vrna_move_init(j, i), VRNA_NEIGHBOR_NEW, data);
      });
    });

    /* base pairs inside of pair */
    FOR_UNPAIRED(pt, i, i = move_i + 1, i <= k, i++, {
      FOR_UNPAIRED(pt, j, j = k + 1, j < move_j, j++, {
        if (is_compatible(fc, i, j))
          cb(fc, vrna_move_init(i, j), VRNA_NEIGHBOR_NEW, data);
      });
    });

    FOR_UNPAIRED(pt, j, j = k, j < move_j, j++, {
      FOR_UNPAIRED(pt, i, i = move_i + 1, i < k, i++, {
        if (is_compatible(fc, i, j))
          cb(fc, vrna_move_init(i, j), VRNA_NEIGHBOR_NEW, data);
      });
    });

    /* 1.2 changed insertion outside shifted pair (i on 5' side) */
    FOR_UNPAIRED(pt, i, i = enclosing_5 + 1, i < MIN2(k, move_i), i++, {
      FOR_UNPAIRED(pt, j, j = i + 1, j < MIN2(k, move_i), j++, {
        if (is_compatible(fc, i, j))
          cb(fc, vrna_move_init(i, j), VRNA_NEIGHBOR_CHANGE, data);
      });

      FOR_UNPAIRED(pt, j, j = MAX2(k, move_j) + 1, j < enclosing_3, j++, {
        if (is_compatible(fc, i, j))
          cb(fc, vrna_move_init(i, j), VRNA_NEIGHBOR_CHANGE, data);
      });
    });

    /* 1.3 insertions outside shifted pair (i on 3' side) */
    FOR_UNPAIRED(pt, i, i = MAX2(k, move_j) + 1, i < enclosing_3, i++, {
      FOR_UNPAIRED(pt, j, j = i + 1, j < enclosing_3, j++, {
        if (is_compatible(fc, i, j))
          cb(fc, vrna_move_init(i, j), VRNA_NEIGHBOR_CHANGE, data);
      });
    });

    /* 1.4 changed insertions within shifted pair. */

    /* 1.4.1 base pair insertions that are now located outside of
     * shifted pair and at 5' side
     */
    FOR_UNPAIRED(pt, i, i = k + 1, i < move_i, i++, {
      FOR_UNPAIRED(pt, j, j = i + 1, j < move_j, j++, {
        if (is_compatible(fc, i, j))
          cb(fc, vrna_move_init(i, j), VRNA_NEIGHBOR_CHANGE, data);
      });
    });

    /* 1.4.2 base pair insertions that are now located outside of
     * shifted pair and at 3' side
     */
    FOR_UNPAIRED(pt, i, i = move_j + 1, i < k, i++, {
      FOR_UNPAIRED(pt, j, j = i + 1, j < k, j++, {
        if (is_compatible(fc, i, j))
          cb(fc, vrna_move_init(i, j), VRNA_NEIGHBOR_CHANGE, data);
      });
    });

    /* 1.4.3 base pairs that are now within shifted pair 5' side */
    FOR_UNPAIRED(pt, i, i = move_i + 1, i < old_i, i++, {
      FOR_UNPAIRED(pt, j, j = i + 1, j < old_i, j++, {
        if (is_compatible(fc, i, j))
          cb(fc, vrna_move_init(i, j), VRNA_NEIGHBOR_CHANGE, data);
      });
    });

    /* 1.4.4 base pairs that are now within shifted pair 3' side */
    FOR_UNPAIRED(pt, i, i = old_j + 1, i < move_j, i++, {
      FOR_UNPAIRED(pt, j, j = i + 1, j < move_j, j++, {
        if (is_compatible(fc, i, j))
          cb(fc, vrna_move_init(i, j), VRNA_NEIGHBOR_CHANGE, data);
      });
    });

    /* 1.4.5 base pairs that are still within the shifted pair */
    FOR_UNPAIRED(pt, i, i = MAX2(old_i, move_i) + 1, i < MIN2(old_j, move_j), i++, {
      FOR_UNPAIRED(pt, j, j = i + 1, j < MIN2(old_j, move_j), j++, {
        if (is_compatible(fc, i, j))
          cb(fc, vrna_move_init(i, j), VRNA_NEIGHBOR_CHANGE, data);
      });
    });
  }

  /* 2. updates and new deletion moves */
  if (options & VRNA_MOVESET_DELETION) {
    k = old_i;
    if ((move_i == k) ||
        (move_j == k))
      k = old_j;

    /* the only new deletion move should be the one that deletes
     * the base pair that arises after the shift
     */
    cb(fc, vrna_move_init(-move_i, -move_j), VRNA_NEIGHBOR_NEW, data);

    /* handle other deletion moves affected by the shift */
    /* 1. enclosing pair */
    if (enclosing_5 > 0)
      cb(fc, vrna_move_init(-enclosing_5, -enclosing_3), VRNA_NEIGHBOR_CHANGE, data);

    /* 2. all base pairs (i, j) with i starting before shifted move */
    deletion_range(fc, pt, enclosing_5 + 1, move_i - 1, VRNA_NEIGHBOR_CHANGE, cb, data);

    /* 3. all base pairs within the shifted move */
    deletion_range(fc, pt, move_i + 1, move_j - 1,VRNA_NEIGHBOR_CHANGE, cb, data);

    /* 4. all base pairs (i, j) with i starting after shifted move */
    deletion_range(fc, pt, move_j + 1, enclosing_3 - 1, VRNA_NEIGHBOR_CHANGE, cb, data);
  }

  /* 3. updates and new shift moves */
  if (options & VRNA_MOVESET_SHIFT) {
    k = old_i;
    if ((move_i == k) ||
        (move_j == k))
      k = old_j;

    /* shift back to the previous configuration */
    if (old_i == move_i)
      cb(fc, vrna_move_init(move_i, -old_j), VRNA_NEIGHBOR_NEW, data);
    else if (old_i == move_j)
      cb(fc, vrna_move_init(move_j, -old_j), VRNA_NEIGHBOR_NEW, data);
    else if (old_j == move_i)
      cb(fc, vrna_move_init(-old_i, move_i), VRNA_NEIGHBOR_NEW, data);
    else
      cb(fc, vrna_move_init(-old_i, move_j), VRNA_NEIGHBOR_NEW, data);

    /* 3.1 All new shift moves that arise due to the novel base pair,
     * i.e. all shifts of the position that stayed constant during the
     * shift
     */
    if ((move_i == old_i) || (move_i == old_j)) { /* move_i was the constant part */
      FOR_UNPAIRED(pt, i, i = enclosing_5 + 1, i < move_i, i++, {
        if (is_compatible(fc, i, move_j))
          cb(fc, vrna_move_init(-i, move_j), VRNA_NEIGHBOR_NEW, data);
      });

      FOR_UNPAIRED(pt, i, i = move_i + 1, i < move_j, i++, {
        if (is_compatible(fc, i, move_j))
          cb(fc, vrna_move_init(-i, move_j), VRNA_NEIGHBOR_NEW, data);
      });

      FOR_UNPAIRED(pt, j, j = move_j + 1, j < enclosing_3, j++, {
        if (is_compatible(fc, move_j, j))
          cb(fc, vrna_move_init(move_j, -j), VRNA_NEIGHBOR_NEW, data);
      });
    } else { /* move_j was the constant part */
      FOR_UNPAIRED(pt, i, i = enclosing_5 + 1, i < move_i, i++, {
        if (is_compatible(fc, i, move_i))
          cb(fc, vrna_move_init(-i, move_i), VRNA_NEIGHBOR_NEW, data);
      });

      FOR_UNPAIRED(pt, j, j = move_i + 1, j < move_j, j++, {
        if (is_compatible(fc, move_i, j))
          cb(fc, vrna_move_init(move_i, -j), VRNA_NEIGHBOR_NEW, data);
      });

      FOR_UNPAIRED(pt, j, j = move_j + 1, j < enclosing_3, j++, {
        if (is_compatible(fc, move_i, j))
          cb(fc, vrna_move_init(move_i, -j), VRNA_NEIGHBOR_NEW, data);
      });
    }

    /* 3.2 all new shift moves that arise from moving one of the pairing partners */
    /* 3.2.1 shifts of the enclosing pair */
    if (enclosing_5 > 0) {
      FOR_UNPAIRED(pt, i, i = k, i < MIN2(k, move_i), i++, {
        if (is_compatible(fc, enclosing_5, i))
          cb(fc, vrna_move_init(enclosing_5, -i), VRNA_NEIGHBOR_NEW, data);
        if (is_compatible(fc, i, enclosing_3))
          cb(fc, vrna_move_init(-i, enclosing_3), VRNA_NEIGHBOR_NEW, data);
      });

      FOR_UNPAIRED(pt, j, j = MIN2(k, move_j) + 1, j <= k, j++, {
        if (is_compatible(fc, enclosing_5, j))
          cb(fc, vrna_move_init(enclosing_5, -j), VRNA_NEIGHBOR_NEW, data);
        if (is_compatible(fc, j, enclosing_3))
          cb(fc, vrna_move_init(-j, enclosing_3), VRNA_NEIGHBOR_NEW, data);
      });
    }

    /* 3.2.2 shifts of other base pairs 5' of the shift */
    FOR_PAIRED(pt, i, i = enclosing_5 + 1, i < MIN2(k, move_i), i++, {
      FOR_UNPAIRED(pt, j, j = k, j < move_i, j++, {
        if (is_compatible(fc, i, j))
          cb(fc, vrna_move_init(i, -j), VRNA_NEIGHBOR_NEW, data);
        if (is_compatible(fc, pt[i], j))
          cb(fc, vrna_move_init(pt[i], -j), VRNA_NEIGHBOR_NEW, data);
      });

      FOR_UNPAIRED(pt, j, j = move_j + 1, j <= k, j++, {
        if (is_compatible(fc, i, j))
          cb(fc, vrna_move_init(i, -j), VRNA_NEIGHBOR_NEW, data);
        if (is_compatible(fc, pt[i], j))
          cb(fc, vrna_move_init(pt[i], -j), VRNA_NEIGHBOR_NEW, data);
      });
    });

    /* 3.2.3 shifts of other base pairs at 3' side of the shift */
    FOR_PAIRED(pt, j, j = MAX2(k, move_j) + 1, j < enclosing_3, j++, {
      FOR_UNPAIRED(pt, i, i = k, i < move_i, i++, {
        if (is_compatible(fc, i, j))
          cb(fc, vrna_move_init(-i, j), VRNA_NEIGHBOR_NEW, data);
        if (is_compatible(fc, i, pt[j]))
          cb(fc, vrna_move_init(-i, pt[j]), VRNA_NEIGHBOR_NEW, data);
      });

      FOR_UNPAIRED(pt, i, i = move_j + 1, i <= k, i++, {
        if (is_compatible(fc, i, j))
          cb(fc, vrna_move_init(-i, j), VRNA_NEIGHBOR_NEW, data);
        if (is_compatible(fc, i, pt[j]))
          cb(fc, vrna_move_init(-i, pt[j]), VRNA_NEIGHBOR_NEW, data);
      });
    });

    /* 3.2.3 shifts of pairs now within the shifted pair. */
    FOR_PAIRED(pt, i, i = move_i + 1, i < old_i, i++, {
      FOR_UNPAIRED(pt, j, j = old_i, j < move_j, j++, {
        if (is_compatible(fc, i, j))
          cb(fc, vrna_move_init(i, -j), VRNA_NEIGHBOR_NEW, data);
        if (is_compatible(fc, pt[i], j))
          cb(fc, vrna_move_init(pt[i], -j), VRNA_NEIGHBOR_NEW, data);
      });
    });

    FOR_PAIRED(pt, i, i = old_j + 1, i < move_j, i++, {
      FOR_UNPAIRED(pt, j, j = move_i + 1, j <= old_j, j++, {
        if (is_compatible(fc, j, i))
          cb(fc, vrna_move_init(-j, i), VRNA_NEIGHBOR_NEW, data);
        if (is_compatible(fc, j, pt[i]))
          cb(fc, vrna_move_init(-j, pt[i]), VRNA_NEIGHBOR_NEW, data);
      });
    });

    if ((move_i < k) && (k < move_j)) {
      /* 3.2.4 novel base pairs inside the shifted pair */
      FOR_PAIRED(pt, i, i = old_i + 1, i < old_j, i++, {
        FOR_UNPAIRED(pt, j, j = move_i + 1, j <= old_i, j++, {
          if (is_compatible(fc, j, i))
            cb(fc, vrna_move_init(-j, i), VRNA_NEIGHBOR_NEW, data);
          if (is_compatible(fc, j, pt[i]))
            cb(fc, vrna_move_init(-j, pt[i]), VRNA_NEIGHBOR_NEW, data);
        });

        FOR_UNPAIRED(pt, j, j = old_j, j < move_j, j++, {
          if (is_compatible(fc, i, j))
            cb(fc, vrna_move_init(i, -j), VRNA_NEIGHBOR_NEW, data);
          if (is_compatible(fc, pt[i], j))
            cb(fc, vrna_move_init(pt[i], -j), VRNA_NEIGHBOR_NEW, data);
        });
      });
    } else {
      /* 3.2.5 base pairs that are now located outside the shifted base pair */
      FOR_PAIRED(pt, i, i = old_i + 1, i < move_i, i++, {
        FOR_UNPAIRED(pt, j, j = enclosing_5 + 1, j < old_i, j++, {
          if (is_compatible(fc, j, i))
            cb(fc, vrna_move_init(-j, i), VRNA_NEIGHBOR_NEW, data);
          if (is_compatible(fc, j, pt[i]))
            cb(fc, vrna_move_init(-j, pt[i]), VRNA_NEIGHBOR_NEW, data);
        });

        FOR_UNPAIRED(pt, j, j = move_j + 1, j < enclosing_3, j++, {
          if (is_compatible(fc, i, j))
            cb(fc, vrna_move_init(i, -j), VRNA_NEIGHBOR_NEW, data);
          if (is_compatible(fc, pt[i], j))
            cb(fc, vrna_move_init(pt[i], -j), VRNA_NEIGHBOR_NEW, data);
        });
      });

      FOR_PAIRED(pt, j, j = move_j + 1, j < old_j, j++, {
        FOR_UNPAIRED(pt, i, i = enclosing_5 + 1, i < move_i, i++, {
          if (is_compatible(fc, i, j))
            cb(fc, vrna_move_init(-i, j), VRNA_NEIGHBOR_NEW, data);
          if (is_compatible(fc, i, pt[j]))
            cb(fc, vrna_move_init(-i, pt[j]), VRNA_NEIGHBOR_NEW, data);
        });

        FOR_UNPAIRED(pt, i, i = move_j + 1, i < enclosing_3, i++, {
          if (is_compatible(fc, j, i))
            cb(fc, vrna_move_init(j, -i), VRNA_NEIGHBOR_NEW, data);
          if (is_compatible(fc, pt[j], i))
            cb(fc, vrna_move_init(pt[j], -i), VRNA_NEIGHBOR_NEW, data);
        });
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
    if ((move_i == old_i) || (move_i == old_j)) {
      /* move_i was the constant part */
      FOR_UNPAIRED(pt, i, i = enclosing_5 + 1, i < MIN2(old_i, move_i), i++, {
        if (is_compatible(fc, i, move_i))
          cb(fc, vrna_move_init(-i, move_i), VRNA_NEIGHBOR_CHANGE, data);
      });

      FOR_UNPAIRED(pt, j, j = move_i + 1, j < MIN2(old_j, move_j), j++, {
        if (is_compatible(fc, move_i, j))
          cb(fc, vrna_move_init(move_i, -j), VRNA_NEIGHBOR_CHANGE, data);
      });

      FOR_UNPAIRED(pt, j, j = MAX2(old_j, move_j) + 1, j < enclosing_3, j++, {
        if (is_compatible(fc, move_i, j))
          cb(fc, vrna_move_init(move_i, -j), VRNA_NEIGHBOR_CHANGE, data);
      });
    } else {
      /* move_j was the constant part */
      FOR_UNPAIRED(pt, i, i = enclosing_5 + 1, i < MIN2(old_i, move_i), i++, {
        if (is_compatible(fc, i, move_j))
          cb(fc, vrna_move_init(-i, move_j), VRNA_NEIGHBOR_CHANGE, data);
      });

      FOR_UNPAIRED(pt, i, i = MAX2(old_i, move_i) + 1, i < move_j, i++, {
        if (is_compatible(fc, i, move_j))
          cb(fc, vrna_move_init(-i , move_j), VRNA_NEIGHBOR_CHANGE, data);
      });

      FOR_UNPAIRED(pt, j, j = MAX2(old_j, move_j) + 1, j < enclosing_3, j++, {
        if (is_compatible(fc, move_j, j))
          cb(fc, vrna_move_init(move_j, -j), VRNA_NEIGHBOR_CHANGE, data);
      });
    }

    /* 3.3.2 shifts of the enclosing base pair (if any) */
    if (enclosing_5 > 0) {
      FOR_UNPAIRED(pt, j, j = enclosing_5 + 1, j < MIN2(k, move_i), j++, {
        if (is_compatible(fc, enclosing_5, j))
          cb(fc, vrna_move_init(enclosing_5, -j), VRNA_NEIGHBOR_CHANGE, data);
      });

      FOR_UNPAIRED(pt, j, j = MAX2(k, move_j) + 1, j < enclosing_3, j++, {
        if (is_compatible(fc, enclosing_5, j))
          cb(fc, vrna_move_init(enclosing_5, -j), VRNA_NEIGHBOR_CHANGE, data);
      });

      FOR_UNPAIRED(pt, i, i = enclosing_5 + 1, j < MIN2(k, move_i), i++, {
        if (is_compatible(fc, i, enclosing_3))
          cb(fc, vrna_move_init(-i, enclosing_3), VRNA_NEIGHBOR_CHANGE, data);
      });

      FOR_UNPAIRED(pt, i, i = MAX2(k, move_j) + 1, i < enclosing_3, i++, {
        if (is_compatible(fc, i, enclosing_3))
          cb(fc, vrna_move_init(-i, enclosing_3), VRNA_NEIGHBOR_CHANGE, data);
      });
    }

    /* shifts of all other base pairs */
    FOR_PAIRED(pt, i, i = enclosing_5 + 1, i < MIN2(k, move_i), i++, {
      FOR_UNPAIRED(pt, j, j = enclosing_5 + 1, j < i, j++, {
        if (is_compatible(fc, j, i))
          cb(fc, vrna_move_init(-j, i), VRNA_NEIGHBOR_CHANGE, data);
        if (is_compatible(fc, j, pt[i]))
          cb(fc, vrna_move_init(-j, pt[i]), VRNA_NEIGHBOR_CHANGE, data);
      });

      FOR_UNPAIRED(pt, j, j = i + 1, j < pt[i], j++, {
        if (is_compatible(fc, i, j))
          cb(fc, vrna_move_init(i, -j), VRNA_NEIGHBOR_CHANGE, data);
        if (is_compatible(fc, j, pt[i]))
          cb(fc, vrna_move_init(-j, pt[i]), VRNA_NEIGHBOR_CHANGE, data);
      });

      FOR_UNPAIRED(pt, j, j = pt[i] + 1, j < MIN2(k, move_i), j++, {
        if (is_compatible(fc, i, j))
          cb(fc, vrna_move_init(i, -j), VRNA_NEIGHBOR_CHANGE, data);
        if (is_compatible(fc, pt[i], j))
          cb(fc, vrna_move_init(pt[i], -j), VRNA_NEIGHBOR_CHANGE, data);
      });

      FOR_UNPAIRED(pt, j, j = MAX2(k, move_j) + 1, j < enclosing_3, j++, {
        if (is_compatible(fc, i, j))
          cb(fc, vrna_move_init(i, -j), VRNA_NEIGHBOR_CHANGE, data);
        if (is_compatible(fc, pt[i], j))
          cb(fc, vrna_move_init(pt[i], -j), VRNA_NEIGHBOR_CHANGE, data);
      });
    });

    FOR_PAIRED(pt, i, i = MAX2(k, move_i) + 1, i < MIN2(k, move_j), i++, {
      FOR_UNPAIRED(pt, j, j = MAX2(k, move_i) + 1, j < i, j++, {
        if (is_compatible(fc, j, i))
          cb(fc, vrna_move_init(-j, i), VRNA_NEIGHBOR_CHANGE, data);
        if (is_compatible(fc, j, pt[i]))
          cb(fc, vrna_move_init(-j, pt[i]), VRNA_NEIGHBOR_CHANGE, data);
      });
      FOR_UNPAIRED(pt, j, j = i + 1, j < pt[i], j++, {
        if (is_compatible(fc, i, j))
          cb(fc, vrna_move_init(i, -j), VRNA_NEIGHBOR_CHANGE, data);
        if (is_compatible(fc, j, pt[i]))
          cb(fc, vrna_move_init(-j, pt[i]), VRNA_NEIGHBOR_CHANGE, data);
      });
      FOR_UNPAIRED(pt, j, j = pt[i] + 1, j < MIN2(k, move_j), j++, {
        if (is_compatible(fc, i, j))
          cb(fc, vrna_move_init(i, -j), VRNA_NEIGHBOR_CHANGE, data);
        if (is_compatible(fc, pt[i], j))
          cb(fc, vrna_move_init(pt[i], -j), VRNA_NEIGHBOR_CHANGE, data);
      });
    });

    FOR_PAIRED(pt, j, j = MAX2(k, move_j) + 1, j < enclosing_3, j++, {
      FOR_UNPAIRED(pt, i, i = enclosing_5 + 1, i < MIN2(k, move_i), i++, {
        if (is_compatible(fc, i, j))
          cb(fc, vrna_move_init(-i, j), VRNA_NEIGHBOR_CHANGE, data);
        if (is_compatible(fc, i, pt[j]))
          cb(fc, vrna_move_init(-i, pt[j]), VRNA_NEIGHBOR_CHANGE, data);
      });

      FOR_UNPAIRED(pt, i, i = MAX2(k, move_j) + 1, i < j, i++, {
        if (is_compatible(fc, i, j))
          cb(fc, vrna_move_init(-i, j), VRNA_NEIGHBOR_CHANGE, data);
        if (is_compatible(fc, i, pt[j]))
          cb(fc, vrna_move_init(-i, pt[j]), VRNA_NEIGHBOR_CHANGE, data);
      });

      FOR_UNPAIRED(pt, i, i = j + 1, i < pt[j], i++, {
        if (is_compatible(fc, j, i))
          cb(fc, vrna_move_init(j, -i), VRNA_NEIGHBOR_CHANGE, data);
        if (is_compatible(fc, i, pt[j]))
          cb(fc, vrna_move_init(-i, pt[j]), VRNA_NEIGHBOR_CHANGE, data);
      });

      FOR_UNPAIRED(pt, i, i = pt[j] + 1, i < enclosing_3, i++, {
        if (is_compatible(fc, j, i))
          cb(fc, vrna_move_init(j, -i), VRNA_NEIGHBOR_CHANGE, data);
        if (is_compatible(fc, pt[j], i))
          cb(fc, vrna_move_init(pt[j], -i), VRNA_NEIGHBOR_CHANGE, data);
      });
    });
  }
}


PRIVATE INLINE void
enclosing_pair(const short *pt,
               int pair_pos_5,
               int pair_pos_3,
               int *enc_pos_5,
               int *enc_pos_3)
{
  int n = pt[0]; /* length of structure */

  *enc_pos_5 = 0;
  *enc_pos_3 = n + 1;

  for (int i = pair_pos_5 - 1; i > 0; i--) {
    if (pt[i] == 0) {
      continue;
    } else if (pt[i] < i) {
      i = pt[i]; /* hop over branching stems */
    } else if (pt[i] > i) {
      /* found enclosing pair */
      *enc_pos_5 = i;
      *enc_pos_3 = pt[i];
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
  int           i, j, enclosing_5, enclosing_3;
  int           min_loop_size = fc->params->model_details.min_loop_size;
  int           move_i        = move->pos_5;
  int           move_j        = move->pos_3;

  /* find out actual enclosing pair that delimits the loop affected by the current move */
  enclosing_pair(pt, move_i, move_j, &enclosing_5, &enclosing_3);

  /*
   *  Determine all previously possible base pairs insertions
   *  that became invalid after inserting (move->pos_5, move->pos_3)
   */
  if (options & VRNA_MOVESET_INSERTION) {
    /*
     *  1. base pairs that start before move->pos_5 and end within
     *  interval [move->pos_5, move->pos_3]
     */
    FOR_UNPAIRED(pt, i, i = enclosing_5 + 1, i < move_i, i++, {
      /* 1.1 remove neighbors that would introduce base pair (i, move_i) */
      if (is_compatible(fc, i, move_i))
        cb(fc, vrna_move_init(i, move_i), VRNA_NEIGHBOR_INVALID, data);

      /*
       * 1.2 remove neighbors that would introduce base pair (i, k)
       * with move_i < k < move_j
       */
      insertions_range_cb(fc,
                          pt,
                          i,
                          move_i + 1,
                          move_j - 1,
                          VRNA_NEIGHBOR_INVALID,
                          cb,
                          data);

      /* 1.3 remove neighbors that would introduce base pair (i, move_j) */
      if (is_compatible(fc, i, move_j))
        cb(fc, vrna_move_init(i, move_j), VRNA_NEIGHBOR_INVALID, data);
    });

    /*
     *  2. base pairs that start at move->pos_5 and end within interval
     *  [move->pos_5, move->pos_3]
     */
    insertions_range_cb(fc,
                        pt,
                        move_i,
                        move_i + 1,
                        move_j - 1,
                        VRNA_NEIGHBOR_INVALID,
                        cb,
                        data);

    /*
     *  2.1 neighbors that would introduce base pairs (move_i, k)
     *  with move_j < k < enclosing_3
     */
    insertions_range_cb(fc,
                        pt,
                        move_i,
                        move_j + 1,
                        enclosing_3 - 1,
                        VRNA_NEIGHBOR_INVALID,
                        cb,
                        data);

    /*
     *  3. base pairs that start in (move->pos_5, move->pos_3) and
     *  end after move->pos_3
     */
    FOR_UNPAIRED(pt, i, i = move_i + 1, i < move_j, i++, {
      /*  3.1 neighbors that introduce base pair (i, move_j) */
      if (is_compatible(fc, i, move_j))
        cb(fc, vrna_move_init(i, move_j), VRNA_NEIGHBOR_INVALID, data);

      /*
       *  3.2 neighbors that introduce base pair (i, k)
       *  with move_j < k < enclosing_3
       */
      insertions_range_cb(fc,
                          pt,
                          i,
                          move_j + 1,
                          enclosing_3 - 1,
                          VRNA_NEIGHBOR_INVALID,
                          cb,
                          data);
    });

    /*
     *  4. base pairs that start at move->pos_3 and end somewhere downstream
     */
    insertions_range_cb(fc,
                        pt,
                        move_j,
                        move_j + 1,
                        enclosing_3 - 1,
                        VRNA_NEIGHBOR_INVALID,
                        cb,
                        data);
  }

  if (options & VRNA_MOVESET_INSERTION) {
    /* no deletions are affected if the next move is an insertion */
  }

  /* next all invalid shift moves */
  if (options & VRNA_MOVESET_SHIFT) {
    /* 1st, the pairing partners of the enclosing pair (if existing) must not
     * interfere with interval [move->pos_5:move->pos_3]
     */
    if (enclosing_5 > 0) {
      FOR_UNPAIRED(pt, j, j = move_i, j <= move_j, j++, {
        if (is_compatible(fc, enclosing_5, j))
          cb(fc, vrna_move_init(enclosing_5, -j), VRNA_NEIGHBOR_INVALID, data);

        if (is_compatible(fc, j, enclosing_3))
          cb(fc, vrna_move_init(-j, enclosing_3), VRNA_NEIGHBOR_INVALID, data);
      });
    }

    /* 2nd, all other base pairs within the same loop and outside
     * the new base pair must not interfere with the interval
     * [move->pos_5:move->pos_3]
     */
    FOR_PAIRED(pt, i, i = enclosing_5 + 1, i < move_i, i++, {
      FOR_UNPAIRED(pt, j, j = move_i, j <= move_j, j++, {
        if (is_compatible(fc, i, j))
          cb(fc, vrna_move_init(i, -j), VRNA_NEIGHBOR_INVALID, data);

        if (is_compatible(fc, pt[i], j))
          cb(fc, vrna_move_init(pt[i], -j), VRNA_NEIGHBOR_INVALID, data);
      });
    });

    FOR_PAIRED(pt, j, j = move_j + 1, j < enclosing_3, j++, {
      FOR_UNPAIRED(pt, i, i = move_i, i <= move_j, i++, {
        if (is_compatible(fc, i, j))
          cb(fc, vrna_move_init(-i, j), VRNA_NEIGHBOR_INVALID, data);

        if (is_compatible(fc, i, pt[j]))
          cb(fc, vrna_move_init(-i, pt[j]), VRNA_NEIGHBOR_INVALID, data);
      });
    });

    /* Finally, invalidate shift moves from pairs that will be enclosed by
     * (move->pos_5, move->pos_3) outside the interval [move->pos_5,move->pos_3]
     */
    FOR_PAIRED(pt, j, j = move_i + 1, j < move_j, j++, {
      FOR_UNPAIRED(pt, i, i = enclosing_5 + 1, i <= move_i, i++, {
        if (is_compatible(fc, i, j))
          cb(fc, vrna_move_init(-i, j), VRNA_NEIGHBOR_INVALID, data);

        if (is_compatible(fc, i, pt[j]))
          cb(fc, vrna_move_init(-i, pt[j]), VRNA_NEIGHBOR_INVALID, data);
      });

      FOR_UNPAIRED(pt, i, i = move_j, i < enclosing_3, i++, {
        if (is_compatible(fc, j, i))
          cb(fc, vrna_move_init(j, -i), VRNA_NEIGHBOR_INVALID, data);

        if (is_compatible(fc, pt[j], i))
          cb(fc, vrna_move_init(pt[j], -i), VRNA_NEIGHBOR_INVALID, data);
      });
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
  int           i, j, enclosing_5, enclosing_3;
  int           min_loop_size = fc->params->model_details.min_loop_size;
  int           move_i        = -move->pos_5;
  int           move_j        = -move->pos_3;

  /* upon deletion of a base pair, conflicts can only
   * arise for shift moves of base pair (move->pos_5, move->pos_3)
   */
  if (options & VRNA_MOVESET_SHIFT) {
    /* find out actual enclosing pair that delimits the loop affected by the current move */
    enclosing_pair(pt, move_i, move_j, &enclosing_5, &enclosing_3);

    /* moves to the region 5' of move->pos_5 */
    FOR_UNPAIRED(pt, i, i = enclosing_5 + 1, i < move_i, i++, {
      if (is_compatible(fc, i, move_i))
        cb(fc, vrna_move_init(-i, move_i), VRNA_NEIGHBOR_INVALID, data);

      if (is_compatible(fc, i, move_j))
        cb(fc, vrna_move_init(-i, move_j), VRNA_NEIGHBOR_INVALID, data);
    });

    /* moves within [move->pos_5:move->pos_3] */
    FOR_UNPAIRED(pt, i, i = move_i + 1, i < move_j, i++, {
      if (is_compatible(fc, move_i, i))
        cb(fc, vrna_move_init(move_i, -i), VRNA_NEIGHBOR_INVALID, data);

      if (is_compatible(fc, i, move_j))
        cb(fc, vrna_move_init(-i, move_j), VRNA_NEIGHBOR_INVALID, data);
    });

    /* moves to the region 3' of move->pos_3 */
    FOR_UNPAIRED(pt, j, j = move_j + 1, j < enclosing_3, j++, {
      if (is_compatible(fc, move_i, j))
        cb(fc, vrna_move_init(move_i, -j), VRNA_NEIGHBOR_INVALID, data);

      if (is_compatible(fc, move_j, j))
        cb(fc, vrna_move_init(move_j, -j), VRNA_NEIGHBOR_INVALID, data);
    });
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
  int           i, j, p, q, move_i, move_j, enclosing_5, enclosing_3;

  /* here, we first determine the base pair (p,q) that is about
   * to change one of its pairing partners
   */
  if (move->pos_5 < 0) {
    p       = pt[move->pos_3];
    q       = move->pos_3;
    move_i  = -move->pos_5;
    move_j  = move->pos_3;
  } else {
    p       = move->pos_5;
    q       = pt[move->pos_5];
    move_i  = move->pos_5;
    move_j  = -move->pos_3;
  }

  if (p > q) {
    i = p;
    p = q;
    q = i;
  }

  if (move_i > move_j) {
    i       = move_i;
    move_i  = move_j;
    move_j  = i;
  }

  /* find out actual enclosing pair that delimits the loop affected by the current move */
  enclosing_pair(pt, p, q, &enclosing_5, &enclosing_3);

  /* we now have the current base pair as (p, q) which will change
   * to (move_i, move_j) or (p', q') after application of the move
   * Here, we distinguish 6 cases, where either p or q stays constant
   * and the pairing partner is shifted towards 3' end of the enclosing
   * loop, the 5' end of the enclosing loop, or towards the segment
   * enclosed by (p,q).
   */
  if (options & VRNA_MOVESET_INSERTION) {
    if (move_i == q) {
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
      j = move_j;
      FOR_UNPAIRED(pt, i, i = enclosing_5 + 1, i < p, i++, {
        if (is_compatible(fc, i, j))
          cb(fc, vrna_move_init(i, j), VRNA_NEIGHBOR_INVALID, data);
      });

      FOR_UNPAIRED(pt, i, i = q + 1, i < move_j, i++, {
        if (is_compatible(fc, i, j))
          cb(fc, vrna_move_init(i, j), VRNA_NEIGHBOR_INVALID, data);
      });

      FOR_UNPAIRED(pt, i, i = move_j + 1, i < enclosing_3, i++, {
        if (is_compatible(fc, j, i))
          cb(fc, vrna_move_init(j, i), VRNA_NEIGHBOR_INVALID, data);
      });

      /* 1.2 invalidate insertion from the outside into interval [q + 1: move_j] */
      FOR_UNPAIRED(pt, j, j = q + 1, j < move_j, j++, {
        FOR_UNPAIRED(pt, i, i = enclosing_5, i < p, i++, {
          if (is_compatible(fc, i, j))
            cb(fc, vrna_move_init(i, j), VRNA_NEIGHBOR_INVALID, data);
        });

        FOR_UNPAIRED(pt, i, i = move_j, i < enclosing_3, i++, {
          if (is_compatible(fc, j, i))
            cb(fc, vrna_move_init(j, i), VRNA_NEIGHBOR_INVALID, data);
        });
      });
    } else if (move_j == p) {
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
      j = move_i;
      FOR_UNPAIRED(pt, i, i = enclosing_5 + 1, i < move_i, i++, {
        if (is_compatible(fc, i, j))
          cb(fc, vrna_move_init(i, j), VRNA_NEIGHBOR_INVALID, data);
      });

      FOR_UNPAIRED(pt, i, i = move_i + 1, i < p, i++, {
        if (is_compatible(fc, j, i))
          cb(fc, vrna_move_init(j, i), VRNA_NEIGHBOR_INVALID, data);
      });

      FOR_UNPAIRED(pt, i, i = q + 1, i < enclosing_3, i++, {
        if (is_compatible(fc, j, i))
          cb(fc, vrna_move_init(j, i), VRNA_NEIGHBOR_INVALID, data);
      });

      /* 2.2 invalidate insertions where one pairing partner is within interval [q':p-1] */
      FOR_UNPAIRED(pt, j, j = move_i + 1, j < p, j++, {
        FOR_UNPAIRED(pt, i, i = enclosing_5 + 1, i < move_i, i++, {
          if (is_compatible(fc, i, j))
            cb(fc, vrna_move_init(i, j), VRNA_NEIGHBOR_INVALID, data);
        });

        FOR_UNPAIRED(pt, i, i = q + 1, i < enclosing_3, i++, {
          if (is_compatible(fc, j, i))
            cb(fc, vrna_move_init(j, i), VRNA_NEIGHBOR_INVALID, data);
        });
      });
    } else if (p - move_i > 0) {
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
      j = move_i;
      FOR_UNPAIRED(pt, i, i = enclosing_5 + 1, i < move_i, i++, {
        if (is_compatible(fc, i, j))
          cb(fc, vrna_move_init(i, j), VRNA_NEIGHBOR_INVALID, data);
      });

      FOR_UNPAIRED(pt, i, i = move_i + 1, i < p, i++, {
        if (is_compatible(fc, j, i))
          cb(fc, vrna_move_init(j, i), VRNA_NEIGHBOR_INVALID, data);
      });

      FOR_UNPAIRED(pt, i, i = q + 1, i < enclosing_3, i++, {
        if (is_compatible(fc, j, i))
          cb(fc, vrna_move_init(j, i), VRNA_NEIGHBOR_INVALID, data);
      });

      /* 3.2 invalidate insertions that span into free'd interval
       * [p':p]
       */
      FOR_UNPAIRED(pt, j, j = move_i + 1, j < p, j++, {
        FOR_UNPAIRED(pt, i, i = enclosing_5 + 1, i < move_i, i++, {
          if (is_compatible(fc, i, j))
            cb(fc, vrna_move_init(i, j), VRNA_NEIGHBOR_INVALID, data);
        });

        FOR_UNPAIRED(pt, i, i = q + 1, i < enclosing_3, i++, {
          if (is_compatible(fc, j, i))
            cb(fc, vrna_move_init(j, i), VRNA_NEIGHBOR_INVALID, data);
        });
      });
    } else if (move_j - q > 0) {
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
      j = move_j;
      FOR_UNPAIRED(pt, i, i = enclosing_5 + 1, i < p, i++, {
        if (is_compatible(fc, i, j))
          cb(fc, vrna_move_init(i, j), VRNA_NEIGHBOR_INVALID, data);
      });

      FOR_UNPAIRED(pt, i, i = q + 1, i < move_j, i++, {
        if (is_compatible(fc, i, j))
          cb(fc, vrna_move_init(i, j), VRNA_NEIGHBOR_INVALID, data);
      });

      FOR_UNPAIRED(pt, i, i = move_j + 1, i < enclosing_3, i++, {
        if (is_compatible(fc, j, i))
          cb(fc, vrna_move_init(j, i), VRNA_NEIGHBOR_INVALID, data);
      });

      /* 4.2 invalidate insertions that span into free'd interval
       * [q:q']
       */
      FOR_UNPAIRED(pt, j, j = q + 1, j < move_j, j++, {
        FOR_UNPAIRED(pt, i, i = enclosing_5 + 1, i < p, i++, {
          if (is_compatible(fc, i, j))
            cb(fc, vrna_move_init(i, j), VRNA_NEIGHBOR_INVALID, data);
        });

        FOR_UNPAIRED(pt, i, i = move_j + 1, i < enclosing_3, i++, {
          if (is_compatible(fc, j, i))
            cb(fc, vrna_move_init(j, i), VRNA_NEIGHBOR_INVALID, data);
        });
      });
    } else if (move_i - p > 0) {
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
      FOR_UNPAIRED(pt, i, i = p + 1, i <= move_i, i++, {
        FOR_UNPAIRED(pt, j, j = move_i, j < q, j++, {
          if (is_compatible(fc, i, j))
            cb(fc, vrna_move_init(i, j), VRNA_NEIGHBOR_INVALID, data);
        });
      });
    } else if (q - move_j > 0) {
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
      FOR_UNPAIRED(pt, i, i = p + 1, i <= move_j, i++, {
        FOR_UNPAIRED(pt, j, j = move_j, j < q, j++, {
          if (is_compatible(fc, i, j))
            cb(fc, vrna_move_init(i, j), VRNA_NEIGHBOR_INVALID, data);
        });
      });
    }
  }

  if (options & VRNA_MOVESET_DELETION) {
    /* the only deletion move invalidated is the original pair that moves */
    cb(fc, vrna_move_init(-p, -q), VRNA_NEIGHBOR_INVALID, data);
  }

  if (options & VRNA_MOVESET_SHIFT) {
    if (p == move_j) {
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
      FOR_UNPAIRED(pt, i, i = enclosing_5 + 1, i < p, i++, {
        if (is_compatible(fc, i, q))
          cb(fc, vrna_move_init(-i, q), VRNA_NEIGHBOR_INVALID, data);
      });

      FOR_UNPAIRED(pt, i, i = p + 1, i < q, i++, {
        if (is_compatible(fc, i, q))
          cb(fc, vrna_move_init(-i, q), VRNA_NEIGHBOR_INVALID, data);
      });

      FOR_UNPAIRED(pt, i, i = q + 1, i < enclosing_3, i++, {
        if (is_compatible(fc, q, i))
          cb(fc, vrna_move_init(q, -i), VRNA_NEIGHBOR_INVALID, data);
      });

      /* 1.2 invalidate moves of enclosing pair positions, if any, into interval [move_i:p-1] */
      if (enclosing_5 > 0) {
        FOR_UNPAIRED(pt, j, j = move_i, j < p, j++, {
          if (is_compatible(fc, enclosing_5, j))
            cb(fc, vrna_move_init(enclosing_5, -j), VRNA_NEIGHBOR_INVALID, data);

          if (is_compatible(fc, j, enclosing_3))
            cb(fc, vrna_move_init(-j, enclosing_3), VRNA_NEIGHBOR_INVALID, data);
        });
      }

      /* 1.3 invalidate moves of other pairs into interval [move_i:p-1] */
      FOR_UNPAIRED(pt, j, j = move_i, j < p, j++, {
        FOR_PAIRED(pt, i, i = enclosing_5 + 1, i < move_i, i++, {
          if (is_compatible(fc, i, j))
            cb(fc, vrna_move_init(i, -j), VRNA_NEIGHBOR_INVALID, data);

          if (is_compatible(fc, pt[i], j))
            cb(fc, vrna_move_init(pt[i], -j), VRNA_NEIGHBOR_INVALID, data);
        });

        FOR_PAIRED(pt, i, i = q + 1, i < enclosing_3, i++, {
          if (is_compatible(fc, j, i))
            cb(fc, vrna_move_init(-j, i), VRNA_NEIGHBOR_INVALID, data);

          if (is_compatible(fc, j, pt[i]))
            cb(fc, vrna_move_init(-j, pt[i]), VRNA_NEIGHBOR_INVALID, data);
        });
      });

      /* 1.4 invalidate moves from within the interval [move_i+1:p-1] outside of the interval */
      FOR_PAIRED(pt, j, j = move_i + 1, j < p, j++, {
        FOR_UNPAIRED(pt, i, i = enclosing_5 + 1, i < move_i, i++, {
          if (is_compatible(fc, i, j))
            cb(fc, vrna_move_init(-i, j), VRNA_NEIGHBOR_INVALID, data);

          if (is_compatible(fc, i, pt[j]))
            cb(fc, vrna_move_init(-i, pt[j]), VRNA_NEIGHBOR_INVALID, data);
        });

        FOR_UNPAIRED(pt, i, i = q + 1, i < enclosing_3, i++, {
          if (is_compatible(fc, j, i))
            cb(fc, vrna_move_init(j, -i), VRNA_NEIGHBOR_INVALID, data);

          if (is_compatible(fc, pt[j], i))
            cb(fc, vrna_move_init(pt[j], -i), VRNA_NEIGHBOR_INVALID, data);
        });
      });
    } else if (q == move_i) {
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
      FOR_UNPAIRED(pt, i, i = enclosing_5 + 1, i < p, i++, {
        if (is_compatible(fc, i, p))
          cb(fc, vrna_move_init(-i, p), VRNA_NEIGHBOR_INVALID, data);
      });

      FOR_UNPAIRED(pt, i, i = p + 1, i < q, i++, {
        if (is_compatible(fc, p, i))
          cb(fc, vrna_move_init(p, -i), VRNA_NEIGHBOR_INVALID, data);
      });

      FOR_UNPAIRED(pt, i, i = q + 1, i < enclosing_3, i++, {
        if (is_compatible(fc, p, i))
          cb(fc, vrna_move_init(p, -i), VRNA_NEIGHBOR_INVALID, data);
      });

      /* 2.2 invalidate moves of enclosing pair, if any, into interval [q+1:move_j] */
      if (enclosing_5 > 0) {
        FOR_UNPAIRED(pt, j, j = q + 1, j <= move_j, j++, {
          if (is_compatible(fc, enclosing_5, j))
            cb(fc, vrna_move_init(enclosing_5, -j), VRNA_NEIGHBOR_INVALID, data);

          if (is_compatible(fc, j, enclosing_3))
            cb(fc, vrna_move_init(-j, enclosing_3), VRNA_NEIGHBOR_INVALID, data);
        });
      }

      /* 2.3 invalidate moves of other pairs into interval [q+1:move_j] */
      FOR_UNPAIRED(pt, j, j = q + 1, j <= move_j, j++, {
        FOR_PAIRED(pt, i, i = enclosing_5 + 1, i < p, i++, {
          if (is_compatible(fc, i, j))
            cb(fc, vrna_move_init(i, -j), VRNA_NEIGHBOR_INVALID, data);

          if (is_compatible(fc, pt[i], j))
            cb(fc, vrna_move_init(pt[i], -j), VRNA_NEIGHBOR_INVALID, data);

        });

        FOR_PAIRED(pt, i, i = move_j + 1, i < enclosing_3, i++, {
          if (is_compatible(fc, j, i))
            cb(fc, vrna_move_init(-j, i), VRNA_NEIGHBOR_INVALID, data);

          if (is_compatible(fc, j, pt[i]))
            cb(fc, vrna_move_init(-j, pt[i]), VRNA_NEIGHBOR_INVALID, data);
        });
      });

      /* 2.4 invalidate move of pairs within interval [q+1:move_j-1] outside of the interval */
      FOR_PAIRED(pt, j, j = q + 1, j < move_j, j++, {
        FOR_UNPAIRED(pt, i, i = enclosing_5 + 1, i < p, i++, {
          if (is_compatible(fc, i, j))
            cb(fc, vrna_move_init(-i, j), VRNA_NEIGHBOR_INVALID, data);

          if (is_compatible(fc, i, pt[j]))
            cb(fc, vrna_move_init(-i, pt[j]), VRNA_NEIGHBOR_INVALID, data);
        });

        FOR_UNPAIRED(pt, i, i = move_j + 1, i < enclosing_3, i++, {
          if (is_compatible(fc, j, i))
            cb(fc, vrna_move_init(j, -i), VRNA_NEIGHBOR_INVALID, data);

          if (is_compatible(fc, pt[j], i))
            cb(fc, vrna_move_init(pt[j], -i), VRNA_NEIGHBOR_INVALID, data);
        });
      });
    } else if (p - move_i > 0) {
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
      FOR_UNPAIRED(pt, i, i = enclosing_5 + 1, i < p, i++, {
        if (is_compatible(fc, i, p))
          cb(fc, vrna_move_init(-i, p), VRNA_NEIGHBOR_INVALID, data);
      });

      FOR_UNPAIRED(pt, i, i = p + 1, i < q, i++, {
        if (is_compatible(fc, p, i))
          cb(fc, vrna_move_init(p, -i), VRNA_NEIGHBOR_INVALID, data);
      });

      FOR_UNPAIRED(pt, j, j = q + 1, j < enclosing_3, j++, {
        if (is_compatible(fc, p, j))
          cb(fc, vrna_move_init(p, -j), VRNA_NEIGHBOR_INVALID, data);
      });

      /* 3.2 invalidate moves of enclosing pair, if any, into interval [move_i:p-1] */
      if (enclosing_5 > 0) {
        FOR_UNPAIRED(pt, i, i = move_i, i < p, i++, {
          if (is_compatible(fc, enclosing_5, i))
            cb(fc, vrna_move_init(enclosing_5, -i), VRNA_NEIGHBOR_INVALID, data);

          if (is_compatible(fc, i, enclosing_3))
            cb(fc, vrna_move_init(-i, enclosing_3), VRNA_NEIGHBOR_INVALID, data);
        });
      }

      /* 3.3 invalidate shift of any other pair into the interval [move_i:p-1] */
      FOR_UNPAIRED(pt, j, j = move_i, j < p, j++, {
        FOR_PAIRED(pt, i, i = enclosing_5 + 1, i < move_i, i++, {
          if (is_compatible(fc, i, j))
            cb(fc, vrna_move_init(i, -j), VRNA_NEIGHBOR_INVALID, data);

          if (is_compatible(fc, pt[i], j))
            cb(fc, vrna_move_init(pt[i], -j), VRNA_NEIGHBOR_INVALID, data);
        });

        FOR_PAIRED(pt, i, i = q + 1, i < enclosing_3, i++, {
          if (is_compatible(fc, j, i))
            cb(fc, vrna_move_init(-j, i), VRNA_NEIGHBOR_INVALID, data);

          if (is_compatible(fc, j, pt[i]))
            cb(fc, vrna_move_init(-j, pt[i]), VRNA_NEIGHBOR_INVALID, data);

        });
      });

      /* 3.4 invalidate shift moves from pairs within interval to the outside */
      FOR_PAIRED(pt, j, j = move_i + 1, j < p, j++, {
        FOR_UNPAIRED(pt, i, i = enclosing_5 + 1, i < move_i, i++, {
          if (is_compatible(fc, i, j))
            cb(fc, vrna_move_init(-i, j), VRNA_NEIGHBOR_INVALID, data);

          if (is_compatible(fc, i, pt[j]))
            cb(fc, vrna_move_init(-i, pt[j]), VRNA_NEIGHBOR_INVALID, data);
        });

        FOR_UNPAIRED(pt, i, i = q + 1, i < enclosing_3, i++, {
          if (is_compatible(fc, j, i))
            cb(fc, vrna_move_init(j, -i), VRNA_NEIGHBOR_INVALID, data);

          if (is_compatible(fc, pt[j], i))
            cb(fc, vrna_move_init(pt[j], -i), VRNA_NEIGHBOR_INVALID, data);
        });
      });
    } else if (move_j - q > 0) {
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
      FOR_UNPAIRED(pt, i, i = enclosing_5 + 1, i < p, i++, {
        if (is_compatible(fc, i, q))
          cb(fc, vrna_move_init(-i, q), VRNA_NEIGHBOR_INVALID, data);
      });

      FOR_UNPAIRED(pt, i, i = p + 1, i < q, i++, {
        if (is_compatible(fc, i, q))
          cb(fc, vrna_move_init(-i, q), VRNA_NEIGHBOR_INVALID, data);
      });

      FOR_UNPAIRED(pt, i, i = q + 1, i < enclosing_3, i++, {
        if (is_compatible(fc, q, i))
          cb(fc, vrna_move_init(q, -i), VRNA_NEIGHBOR_INVALID, data);
      });

      /* 4.2 invalidate moves of enclosing pair, if any, into interval [q+1:move_j] */
      if (enclosing_5 > 0) {
        FOR_UNPAIRED(pt, j, j = q + 1, j <= move_j, j++, {
          if (is_compatible(fc, enclosing_5, j))
            cb(fc, vrna_move_init(enclosing_5, -j), VRNA_NEIGHBOR_INVALID, data);

          if (is_compatible(fc, j, enclosing_3))
            cb(fc, vrna_move_init(-j, enclosing_3), VRNA_NEIGHBOR_INVALID, data);
        });
      }

      /* 4.3 invalidate shifts from other pairs into the interval [q+1: move_j] */
      FOR_UNPAIRED(pt, j, j = q + 1, j <= move_j, j++, {
        FOR_PAIRED(pt, i, i = enclosing_5 + 1, i < p, i++, {
          if (is_compatible(fc, i, j))
            cb(fc, vrna_move_init(i, -j), VRNA_NEIGHBOR_INVALID, data);

          if (is_compatible(fc, pt[i], j))
            cb(fc, vrna_move_init(pt[i], -j), VRNA_NEIGHBOR_INVALID, data);

        });

        FOR_PAIRED(pt, i, i = move_j + 1, i < enclosing_3, i++, {
          if (is_compatible(fc, j, i))
            cb(fc, vrna_move_init(-j, i), VRNA_NEIGHBOR_INVALID, data);

          if (is_compatible(fc, j, pt[i]))
            cb(fc, vrna_move_init(-j, pt[i]), VRNA_NEIGHBOR_INVALID, data);
        });
      });

      /* 4.4 invalidate shifts from within the interval [q + 1: move_j - 1] to the outside */
      FOR_PAIRED(pt, j, j = q + 1, j < move_j, j++, {
        FOR_UNPAIRED(pt, i, i = enclosing_5, i < p, i++, {
          if (is_compatible(fc, i, j))
            cb(fc, vrna_move_init(-i, j), VRNA_NEIGHBOR_INVALID, data);

          if (is_compatible(fc, i, pt[j]))
            cb(fc, vrna_move_init(-i, pt[j]), VRNA_NEIGHBOR_INVALID, data);
        });

        FOR_UNPAIRED(pt, i, i = move_j + 1, i < enclosing_3, i++, {
          if (is_compatible(fc, j, i))
            cb(fc, vrna_move_init(j, -i), VRNA_NEIGHBOR_INVALID, data);

          if (is_compatible(fc, pt[j], i))
            cb(fc, vrna_move_init(pt[j], -i), VRNA_NEIGHBOR_INVALID, data);
        });
      });
    } else if (move_i - p > 0) {
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
      FOR_UNPAIRED(pt, i, i = enclosing_5 + 1, i < p, i++, {
        if (is_compatible(fc, i, p))
          cb(fc, vrna_move_init(-i, p), VRNA_NEIGHBOR_INVALID, data);
      });

      FOR_UNPAIRED(pt, i, i = p + 1, i < q, i++, {
        if (is_compatible(fc, p, i))
          cb(fc, vrna_move_init(p, -i), VRNA_NEIGHBOR_INVALID, data);
      });

      FOR_UNPAIRED(pt, i, i = q + 1, i < enclosing_3, i++, {
        if (is_compatible(fc, p, i))
          cb(fc, vrna_move_init(p, -i), VRNA_NEIGHBOR_INVALID, data);
      });

      /* 5.2 invalidate moves of pairs into interval [p + 1:move_i] */
      FOR_UNPAIRED(pt, i, i = p + 1, i <= move_i, i++, {
        FOR_PAIRED(pt, j, j = move_i + 1, j < q, j++, {
          if (is_compatible(fc, i, j))
            cb(fc, vrna_move_init(-i, j), VRNA_NEIGHBOR_INVALID, data);

          if (is_compatible(fc, i, pt[j]))
            cb(fc, vrna_move_init(-i, pt[j]), VRNA_NEIGHBOR_INVALID, data);
        });
      });

      /* 5.3 invalidate moves of pairs from within interval [p + 1: move_i] */
      FOR_PAIRED(pt, i, i = p + 1, i < move_i, i++, {
        FOR_UNPAIRED(pt, j, j = move_i, j < q, j++, {
          if (is_compatible(fc, i, j))
            cb(fc, vrna_move_init(i, -j), VRNA_NEIGHBOR_INVALID, data);

          if (is_compatible(fc, pt[i], j))
            cb(fc, vrna_move_init(pt[i], -j), VRNA_NEIGHBOR_INVALID, data);
        });
      });
    } else if (q - move_j > 0) {
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
      FOR_UNPAIRED(pt, i, i = enclosing_5, i < p, i++, {
        if (is_compatible(fc, i, q))
          cb(fc, vrna_move_init(-i, q), VRNA_NEIGHBOR_INVALID, data);
      });

      FOR_UNPAIRED(pt, i, i = p + 1, i < q, i++, {
        if (is_compatible(fc, i, q))
          cb(fc, vrna_move_init(-i, q), VRNA_NEIGHBOR_INVALID, data);
      });

      FOR_UNPAIRED(pt, i, i = q + 1, i < enclosing_3, i++,{
        if (is_compatible(fc, q, i))
          cb(fc, vrna_move_init(q, -i), VRNA_NEIGHBOR_INVALID, data);
      });

      /* 6.2 invalidate shifts of pairs into interval [move_j:q-1] */
      FOR_UNPAIRED(pt, j, j = move_j, j < q, j++, {
        FOR_PAIRED(pt, i, i = p + 1, i < move_j, i++, {
          if (is_compatible(fc, i, j))
            cb(fc, vrna_move_init(i, -j), VRNA_NEIGHBOR_INVALID, data);

          if (is_compatible(fc, pt[i], j))
            cb(fc, vrna_move_init(pt[i], -j), VRNA_NEIGHBOR_INVALID, data);
        });
      });

      /* 6.3 invalidate shifts from within interval [move_j+1:q-1] */
      FOR_PAIRED(pt, j, j = move_j + 1, j < q, j++, {
        FOR_UNPAIRED(pt, i, i = p + 1, i <= move_j, i++, {
          if (is_compatible(fc, i, j))
            cb(fc, vrna_move_init(-i, j), VRNA_NEIGHBOR_INVALID, data);

          if (is_compatible(fc, i, pt[j]))
            cb(fc, vrna_move_init(-i, pt[j]), VRNA_NEIGHBOR_INVALID, data);
        });
      });
    }
  }
}
