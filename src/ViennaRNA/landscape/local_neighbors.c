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


PRIVATE void
generate_local_nb(vrna_fold_compound_t      *fc,
                  const short               *pt,
                  const vrna_move_t         *move,
                  vrna_move_update_f cb,
                  void                      *data,
                  unsigned int              options);


PRIVATE void
generate_local_nb_insertion(vrna_fold_compound_t      *fc,
                            const short               *pt,
                            const vrna_move_t         *move,
                            vrna_move_update_f cb,
                            void                      *data,
                            unsigned int              options);


PRIVATE void
generate_local_nb_deletion(vrna_fold_compound_t       *fc,
                           const short                *pt,
                           const vrna_move_t          *move,
                           vrna_move_update_f  cb,
                           void                       *data,
                           unsigned int               options);


PRIVATE void
generate_local_nb_shift(vrna_fold_compound_t      *fc,
                        const short               *pt,
                        const vrna_move_t         *move,
                        vrna_move_update_f cb,
                        void                      *data,
                        unsigned int              options);


PRIVATE void
generate_conflicts_local_nb(vrna_fold_compound_t      *fc,
                            const short               *pt,
                            const vrna_move_t         *move,
                            vrna_move_update_f cb,
                            void                      *data,
                            unsigned int              options);


PRIVATE void
generate_conflicts_local_nb_insertion(vrna_fold_compound_t      *fc,
                                      const short               *pt,
                                      const vrna_move_t         *move,
                                      vrna_move_update_f cb,
                                      void                      *data,
                                      unsigned int              options);

PRIVATE void
generate_conflicts_local_nb_deletion(vrna_fold_compound_t      *fc,
                                     const short               *pt,
                                     const vrna_move_t         *move,
                                     vrna_move_update_f cb,
                                     void                      *data,
                                     unsigned int              options);

PRIVATE void
generate_conflicts_local_nb_shift(vrna_fold_compound_t      *fc,
                                  const short               *pt,
                                  const vrna_move_t         *move,
                                  vrna_move_update_f cb,
                                  void                      *data,
                                  unsigned int              options);


#include "landscape/local_neighbors.inc"

PRIVATE INLINE int
is_compatible(const vrna_fold_compound_t  *vc,
              int                         i,
              int                         j)
{
  return vc->params->model_details.pair[vc->sequence_encoding2[i]][vc->sequence_encoding2[j]] != 0; /* see pair_mat.h */
}

PUBLIC int
vrna_move_neighbor_diff_cb(vrna_fold_compound_t       *fc,
                           short                      *ptable,
                           vrna_move_t                move,
                           vrna_move_update_f  cb,
                           void                       *data,
                           unsigned int               options)
{
  if ((fc) && (ptable) && (cb)) {
    /* crude check if ptable has correct size */
    if ((unsigned int)ptable[0] == fc->length) {

      /* 1. remove the neighbor we are about to change into */
      cb(fc, move, VRNA_NEIGHBOR_INVALID, data);

      /* 2. remove neighbors that become invalid after application of 'move' */
      generate_conflicts_local_nb(fc, ptable, &move, cb, data, options);

      /* 3. Detect novel neighbors and those that require an update after application of 'move' */
      generate_local_nb(fc, ptable, &move, cb, data, options);

      vrna_move_apply(ptable, &move);
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
    invalid_moves = NULL;

  return valid_neighbors;
}


PRIVATE void
generate_local_nb(vrna_fold_compound_t      *fc,
                  const short               *pt,
                  const vrna_move_t         *move,
                  vrna_move_update_f cb,
                  void                      *data,
                  unsigned int              options)
{
  /* check whether move is valid given the options */
  if ((move->pos_5 > 0) &&
      (move->pos_3 > 0) &&
      (options & VRNA_MOVESET_INSERTION))
    generate_local_nb_insertion(fc, pt, move, cb, data, options);
  else if ((move->pos_5 < 0) &&
           (move->pos_3 < 0) &&
           (options & VRNA_MOVESET_DELETION))
    generate_local_nb_deletion(fc, pt, move, cb, data, options);
  else if (options & VRNA_MOVESET_SHIFT)
    generate_local_nb_shift(fc, pt, move, cb, data, options);
}


PRIVATE void
generate_conflicts_local_nb(vrna_fold_compound_t      *fc,
                            const short               *pt,
                            const vrna_move_t         *move,
                            vrna_move_update_f cb,
                            void                      *data,
                            unsigned int              options)
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
insertions_range_cb(vrna_fold_compound_t      *fc,
                    const short               *pt,
                    int                       i,
                    int                       min_loop_size,
                    int                       last_j,
                    unsigned int              status,
                    vrna_move_update_f cb,
                    void                      *data)
{
  int j;

  for (j = i + 1; j < MIN2(i + min_loop_size, last_j) + 1; j++)
    if (pt[j] > j)
      j = pt[j];

  if ((j > last_j) || (j < i + min_loop_size + 1))
    return;

  while (pt[j] > j)
    j = pt[j];

  for (; j <= last_j; j++) {
    if (pt[j] > j)
      j = pt[j]; /* hop over branching stems */
    else if ((pt[j] == 0) && (is_compatible(fc, i, j)))
      cb(fc, vrna_move_init(i, j), status, data);
  }
}


PRIVATE INLINE void
insertions_range_cb2(vrna_fold_compound_t       *fc,
                     const short                *pt,
                     int                        i,
                     int                        min_j,
                     int                        last_j,
                     unsigned int               status,
                     vrna_move_update_f  cb,
                     void                       *data)
{
  int j;

  for (j = min_j; j <= last_j; j++) {
    if (pt[j] > j)
      j = pt[j]; /* hop over branching stems */
    else if ((pt[j] == 0) && (is_compatible(fc, i, j)))
      cb(fc, vrna_move_init(i, j), status, data);
  }
}


PRIVATE INLINE void
deletions_range_cb(vrna_fold_compound_t       *fc,
                   const short                *pt,
                   int                        i,
                   int                        first_j,
                   int                        last_j,
                   unsigned int               status,
                   vrna_move_update_f  cb,
                   void                       *data)
{
  int j;

  for (j = first_j; j <= last_j; j++) {
    if (pt[j] > j)
      j = pt[j]; /* hop over branching stems */
    else if ((pt[j] == 0) && (is_compatible(fc, i, j)))
      cb(fc, vrna_move_init(i, j), status, data);
  }
}


PRIVATE void
generate_local_nb_insertion(vrna_fold_compound_t      *fc,
                            const short               *pt,
                            const vrna_move_t         *move,
                            vrna_move_update_f cb,
                            void                      *data,
                            unsigned int              options)
{
  unsigned int  n = fc->length;
  int           i;
  int           enclosing_5   = 0;          /* exterior loop */
  int           enclosing_3   = (int)n + 1; /* exterior loop */
  int           min_loop_size = fc->params->model_details.min_loop_size;

  /* find out actual enclosing pair that delimits the loop affected by the current move */
  for (i = move->pos_5 - 1; i > 0; i--) {
    if (pt[i] == 0) {
      continue;
    } else if (pt[i] < i) {
      i = pt[i]; /* hop over branching stems */
    } else if (pt[i] > i) {
      /* found enclosing pair */
      enclosing_5 = i;
      enclosing_3 = pt[i];
      break;
    }
  }

  /* 1. valid base pair deletions */
  if (options & VRNA_MOVESET_DELETION) {
    /* 1.1 removal of enclosing pair */
    if (enclosing_5 > 0)
      cb(fc, vrna_move_init(-enclosing_5, -enclosing_3), VRNA_NEIGHBOR_CHANGE, data);

    /* 1.2 removal of base pair that has just been inserted */
    cb(fc, vrna_move_init(-move->pos_5, -move->pos_3), VRNA_NEIGHBOR_NEW, data);

    /* 1.3 removal of other base pairs that branch-off from the current loop */
    for (i = enclosing_5 + 1; i < move->pos_5; i++)
      if (pt[i] > i) {
        cb(fc, vrna_move_init(-i, -pt[i]), VRNA_NEIGHBOR_CHANGE, data);
        i = pt[i]; /* hop over branching stems */
      }

    for (i = move->pos_5 + 1; i < move->pos_3; i++)
      if (pt[i] > i) {
        cb(fc, vrna_move_init(-i, -pt[i]), VRNA_NEIGHBOR_CHANGE, data);
        i = pt[i]; /* hop over branching stems */
      }

    for (i = move->pos_3 + 1; i < enclosing_3; i++)
      if (pt[i] > i) {
        cb(fc, vrna_move_init(-i, -pt[i]), VRNA_NEIGHBOR_CHANGE, data);
        i = pt[i]; /* hop over branching stems */
      }
  }

  /* 2. valid base pair insertions */
  if (options & VRNA_MOVESET_INSERTION) {
    /* 5' side of newly inserted base pair */
    for (i = enclosing_5 + 1; i < move->pos_5; i++) {
      if (pt[i] > i) {
        i = pt[i]; /* hop over branching stems */
      } else {
        /*
         * determine all potential pairing partners for current nucleotide within
         * 5' side of newly inserted base pair
         */
        insertions_range_cb(fc,
                            pt,
                            i,
                            min_loop_size,
                            move->pos_5 - 1,
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
                            move->pos_3 - i,
                            enclosing_3 - 1,
                            VRNA_NEIGHBOR_CHANGE,
                            cb,
                            data);
      }
    }

    /* within newly inserted base pair */
    for (i = move->pos_5 + 1; i < move->pos_3; i++) {
      if (pt[i] > i) {
        i = pt[i]; /* hop over branching stems */
      } else {
        insertions_range_cb(fc,
                            pt,
                            i,
                            min_loop_size,
                            move->pos_3 - 1,
                            VRNA_NEIGHBOR_CHANGE,
                            cb,
                            data);
      }
    }

    /* 3' side of newly inserted base pair */
    for (i = move->pos_3 + 1; i < enclosing_3; i++) {
      if (pt[i] > i) {
        i = pt[i]; /* hop over branching stems */
      } else {
        /*
         * determine all potential pairing partners for current nucleotide within
         * 5' side of newly inserted base pair
         */
        insertions_range_cb(fc,
                            pt,
                            i,
                            min_loop_size,
                            enclosing_3 - 1,
                            VRNA_NEIGHBOR_CHANGE,
                            cb,
                            data);
      }
    }
  }

  if (options & VRNA_MOVESET_SHIFT) {
    /* not implemented yet */
  }
}


PRIVATE void
generate_local_nb_deletion(vrna_fold_compound_t       *fc,
                           const short                *pt,
                           const vrna_move_t          *move,
                           vrna_move_update_f  cb,
                           void                       *data,
                           unsigned int               options)
{
  unsigned int  n = fc->length;
  int           i, j;
  int           enclosing_5   = 0;          /* exterior loop */
  int           enclosing_3   = (int)n + 1; /* exterior loop */
  int           min_loop_size = fc->params->model_details.min_loop_size;
  int           move_i        = -move->pos_5;
  int           move_j        = -move->pos_3;

  /* find out actual enclosing pair that delimits the loop affected by the current move */
  for (i = move_i - 1; i > 0; i--) {
    if (pt[i] == 0) {
      continue;
    } else if (pt[i] < i) {
      i = pt[i]; /* hop over branching stems */
    } else if (pt[i] > i) {
      /* found enclosing pair */
      enclosing_5 = i;
      enclosing_3 = pt[i];
      break;
    }
  }


  /* actally generate the list of valid moves */

  /* 1. valid base pair deletions */
  if (options & VRNA_MOVESET_DELETION) {
    /* 1.1 removal of enclosing pair */
    if (enclosing_5 > 0)
      cb(fc, vrna_move_init(-enclosing_5, -enclosing_3), VRNA_NEIGHBOR_CHANGE, data);

    /* 1.2 branching base pairs 5' to base pair deleted by current move */
    for (i = enclosing_5 + 1; i < move_i; i++) {
      if (pt[i] > i) {
        cb(fc, vrna_move_init(-i, -pt[i]), VRNA_NEIGHBOR_CHANGE, data);
        i = pt[i]; /* hop over branching stems */
      }
    }

    /* 1.3 branching base pairs enclosed by base pair deleted by current move */
    for (i = move_i + 1; i < move_j; i++) {
      if (pt[i] > i) {
        cb(fc, vrna_move_init(-i, -pt[i]), VRNA_NEIGHBOR_CHANGE, data);
        i = pt[i]; /* hop over branching stems */
      }
    }

    /* 1.4 branching base pairs on 3' side of base pair deleted by current move */
    for (i = move_j + 1; i < enclosing_3; i++) {
      if (pt[i] > i) {
        cb(fc, vrna_move_init(-i, -pt[i]), VRNA_NEIGHBOR_CHANGE, data);
        i = pt[i]; /* hop over branching stems */
      }
    }
  }

  /* 2. valid base pair insertions */
  if (options & VRNA_MOVESET_INSERTION) {
    /* 2.1 re-insertion of base pair removed by current move */
    cb(fc, vrna_move_init(move_i, move_j), VRNA_NEIGHBOR_NEW, data);

    /* 2.2 insertion of novel pairs that start on 5' side of current move */
    for (i = enclosing_5 + 1; i < move_i; i++) {
      if (pt[i] > i) {
        i = pt[i]; /* hop over branching stems */
      } else {
        /* 2.2.1 base pairs that end before 5' side of current move */
        insertions_range_cb(fc,
                            pt,
                            i,
                            min_loop_size,
                            move_i - 1,
                            VRNA_NEIGHBOR_CHANGE,
                            cb,
                            data);

        /* 2.2.2 base pairs that end at 5' side of current move */
        if (is_compatible(fc, i, move_i) && (move_i - i > min_loop_size))
          cb(fc, vrna_move_init(i, move_i), VRNA_NEIGHBOR_NEW, data);

        /* 2.2.3 base pairs that end within loop that was delimited by current move */
        insertions_range_cb(fc,
                            pt,
                            i,
                            MAX2(move_i - i, min_loop_size),
                            move_j - 1,
                            VRNA_NEIGHBOR_NEW,
                            cb,
                            data);

        /* 2.2.4 base pairs that end at 3' position of current loop */
        if (is_compatible(fc, i, move_j))
          cb(fc, vrna_move_init(i, move_j), VRNA_NEIGHBOR_NEW, data);

        /* 2.2.5 base pairs that end after 3' nucleotide of current move */
        insertions_range_cb(fc,
                            pt,
                            i,
                            move_j - i,
                            enclosing_3 - 1,
                            VRNA_NEIGHBOR_CHANGE,
                            cb,
                            data);
      }
    }

    /* 2.3 insertion of novel pairs that start at 5' nucleotide of current position */
    i = move_i;

    /* 2.3.1 ending within loop enclosed by base pair deleted by current move */
    insertions_range_cb(fc,
                        pt,
                        i,
                        min_loop_size,
                        move_j - 1,
                        VRNA_NEIGHBOR_NEW,
                        cb,
                        data);

    /* 2.3.2 ending after loop enclosed by base pair deleted by current move */
    insertions_range_cb(fc,
                        pt,
                        i,
                        move_j - i,
                        enclosing_3 - 1,
                        VRNA_NEIGHBOR_NEW,
                        cb,
                        data);

    /* 2.4 insertion of novel pairs that start within loop closed by current move */
    for (i = move_i + 1; i < move_j; i++) {
      if (pt[i] > i) {
        i = pt[i]; /* hop over branching stems */
      } else {
        /* 2.4.1 */
        insertions_range_cb(fc,
                            pt,
                            i,
                            min_loop_size,
                            move_j - 1,
                            VRNA_NEIGHBOR_CHANGE,
                            cb,
                            data);

        /* 2.4.2 */
        j = move_j;

        if ((is_compatible(fc, i, j)) && (move_j - i > min_loop_size))
          cb(fc, vrna_move_init(i, j), VRNA_NEIGHBOR_NEW, data);

        /* 2.4.3 */
        insertions_range_cb(fc,
                            pt,
                            i,
                            MAX2(move_j - i, min_loop_size),
                            enclosing_3 - 1,
                            VRNA_NEIGHBOR_NEW,
                            cb,
                            data);
      }
    }

    /* 2.5 insertion of novel pairs that start at 3' nucleotide of current move */
    i = move_j;

    insertions_range_cb(fc,
                        pt,
                        i,
                        min_loop_size,
                        enclosing_3 - 1,
                        VRNA_NEIGHBOR_NEW,
                        cb,
                        data);

    /* 2.6 insertion of novel pairs that start after 3' nucleotide of current move */
    for (i = move_j + 1; i < enclosing_3; i++) {
      if (pt[i] > i) {
        i = pt[i]; /* hop over branching stems */
      } else {
        insertions_range_cb(fc,
                            pt,
                            i,
                            min_loop_size,
                            enclosing_3 - 1,
                            VRNA_NEIGHBOR_CHANGE,
                            cb,
                            data);
      }
    }
  }

  /* 3. valid shift moves */
  if (options & VRNA_MOVESET_SHIFT) {
    /* not implemented yet */
  }
}


PRIVATE void
generate_local_nb_shift(vrna_fold_compound_t      *fc,
                        const short               *pt,
                        const vrna_move_t         *move,
                        vrna_move_update_f cb,
                        void                      *data,
                        unsigned int              options)
{
}


PRIVATE void
generate_conflicts_local_nb_insertion(vrna_fold_compound_t      *fc,
                                      const short               *pt,
                                      const vrna_move_t         *move,
                                      vrna_move_update_f cb,
                                      void                      *data,
                                      unsigned int              options)
{
  unsigned int  n = fc->length;
  int           i, j;
  int           enclosing_5   = 0;          /* exterior loop */
  int           enclosing_3   = (int)n + 1; /* exterior loop */
  int           min_loop_size = fc->params->model_details.min_loop_size;
  int           move_i        = move->pos_5;
  int           move_j        = move->pos_3;

  /* find out actual enclosing pair that delimits the loop affected by the current move */
  for (i = move_i - 1; i > 0; i--) {
    if (pt[i] == 0) {
      continue;
    } else if (pt[i] < i) {
      i = pt[i]; /* hop over branching stems */
    } else if (pt[i] > i) {
      /* found enclosing pair */
      enclosing_5 = i;
      enclosing_3 = pt[i];
      break;
    }
  }

  /* no deletions are affected if the next move is an insertion */

  /*
   *  Determine all previously possible base pairs insertions
   *  that became invalid after inserting (move->pos_5, move->pos_3)
   */
  if (options & VRNA_MOVESET_INSERTION) {
    /*
     *  1. base pairs that start before move->pos_5 and end in
     *  [move->pos_5, move->pos_3]
     */
    for (i = enclosing_5 + 1; i < move_i; i++) {
      if (pt[i] > i) {
        i = pt[i]; /* hop over branching stems */
      } else if (pt[i] == 0) {
        /* 1.1 remove neighbors that would introduce base pair (i, move_i) */
        if ((is_compatible(fc, i, move_i)) && (move_i - i > min_loop_size))
          cb(fc, vrna_move_init(i, move_i), VRNA_NEIGHBOR_INVALID, data);

        /*
         * 1.2 remove neighbors that would introduce base pair (i, k)
         * with move_i < k < move_j
         */
        insertions_range_cb2(fc,
                             pt,
                             i,
                             MAX2(move_i, i + min_loop_size) + 1,
                             move_j - 1,
                             VRNA_NEIGHBOR_INVALID,
                             cb,
                             data);

        /* 1.3 remove neighbors that would introduce base pair (i, move_j) */
        if (is_compatible(fc, i, move_j))
          cb(fc, vrna_move_init(i, move_j), VRNA_NEIGHBOR_INVALID, data);
      }
    }

    i = move_i;

    /*
     *  2. base pairs that start at move->pos_5 and end within interval (move_i, move_j)
     */
    insertions_range_cb(fc,
                        pt,
                        i,
                        min_loop_size,
                        move_j - 1,
                        VRNA_NEIGHBOR_INVALID,
                        cb,
                        data);

    /*
     *  2.1 neighbors that would introduce base pairs (move_i, k)
     *  with move_j < k < enclosing_3
     */
    insertions_range_cb2(fc,
                         pt,
                         i,
                         move_j + 1,
                         enclosing_3 - 1,
                         VRNA_NEIGHBOR_INVALID,
                         cb,
                         data);

    /*
     *  3. base pairs that start in (move->pos_5, move->pos_3) and
     *  end after move->pos_3
     */
    for (i = move_i + 1; i < move_j; i++) {
      if (pt[i] > i) {
        i = pt[i]; /* hop over branching stems */
      } else if (pt[i] == 0) {
        /*  3.1 neighbors that introduce base pair (i, move_j) */
        if ((is_compatible(fc, i, move_j)) && (move_j - i > min_loop_size))
          cb(fc, vrna_move_init(i, move_j), VRNA_NEIGHBOR_INVALID, data);

        /*
         *  3.2 neighbors that introduce base pair (i, k)
         *  with move_j < k < enclosing_3
         */
        insertions_range_cb(fc,
                            pt,
                            i,
                            MAX2(move_j - i, min_loop_size),
                            enclosing_3 - 1,
                            VRNA_NEIGHBOR_INVALID,
                            cb,
                            data);
      }
    }

    /*
     *  4. base pairs that start at move->pos_3 and end somewhere downstream
     */
    i = move_j;
    insertions_range_cb(fc,
                        pt,
                        i,
                        min_loop_size,
                        enclosing_3 - 1,
                        VRNA_NEIGHBOR_INVALID,
                        cb,
                        data);
  }

  /* next all invalid shift moves */
  if (options & VRNA_MOVESET_SHIFT) {
    /* 1st, the pairing partners of the enclosing pair (if existing) must not
     * move into the range [move->pos_5:move->pos_3]
     */
    if ((enclosing_5 > 0) &&
        (enclosing_3 <= n)) {

      if (is_compatible(fc, enclosing_5, move_i))
        cb(fc, vrna_move_init(enclosing_5, -move_i), VRNA_NEIGHBOR_INVALID, data);

      if (is_compatible(fc, move_j, enclosing_3))
        cb(fc, vrna_move_init(-move_j, enclosing_3), VRNA_NEIGHBOR_INVALID, data);

      for (j = move_i + 1; j < move_j; j++) {
        if (pt[j] > j) {
          j = pt[j];
        } else {
          if (is_compatible(fc, enclosing_5, j))
            cb(fc, vrna_move_init(enclosing_5, -j), VRNA_NEIGHBOR_INVALID, data);

          if (is_compatible(fc, j, enclosing_3))
            cb(fc, vrna_move_init(-j, enclosing_3), VRNA_NEIGHBOR_INVALID, data);
        }
      }
    }

    /* 2nd, all other base pairs within the same loop must not move
     * there pairing partners into the interval [move->pos_5:move->pos_3]
     */
    for (i = enclosing_5; i < move_i; i++) {
      if (pt[i] > i) {
        if (is_compatible(fc, i, move_i))
           cb(fc, vrna_move_init(i, -move_i), VRNA_NEIGHBOR_INVALID, data);

        if (is_compatible(fc, i, move_j))
           cb(fc, vrna_move_init(i, -move_j), VRNA_NEIGHBOR_INVALID, data);

        if (is_compatible(fc, pt[i], move_i))
           cb(fc, vrna_move_init(pt[i], -move_i), VRNA_NEIGHBOR_INVALID, data);

        if (is_compatible(fc, pt[i], move_j))
           cb(fc, vrna_move_init(pt[i], -move_j), VRNA_NEIGHBOR_INVALID, data);

        for (j = move_i + 1; j < move_j; j++) {
          if (pt[j] > j) {
            j = pt[j];
          } else {
            if (is_compatible(fc, i, j))
              cb(fc, vrna_move_init(i, -j), VRNA_NEIGHBOR_INVALID, data);

            if (is_compatible(fc, pt[i], j))
              cb(fc, vrna_move_init(pt[i], -j), VRNA_NEIGHBOR_INVALID, data);
          }
        }

        /* hop over base pair */
        i = pt[i];
      }
    }

    for (j = move_j + 1; j < enclosing_3; j++) {
      if (pt[j] > j) {
        if (is_compatible(fc, move_i, j))
          cb(fc, vrna_move_init(-move_i, j), VRNA_NEIGHBOR_INVALID, data);

        if (is_compatible(fc, move_j, j))
          cb(fc, vrna_move_init(-move_j, j), VRNA_NEIGHBOR_INVALID, data);

        if (is_compatible(fc, move_i, pt[j]))
          cb(fc, vrna_move_init(-move_i, pt[j]), VRNA_NEIGHBOR_INVALID, data);

        if (is_compatible(fc, move_j, pt[j]))
          cb(fc, vrna_move_init(-move_j, pt[j]), VRNA_NEIGHBOR_INVALID, data);

        for (i = move_i + 1; i < move_j; i++) {
          if (pt[i] > i) {
            i = pt[i];
          } else {
            if (is_compatible(fc, i, j))
              cb(fc, vrna_move_init(-i, j), VRNA_NEIGHBOR_INVALID, data);

            if (is_compatible(fc, i, pt[j]))
              cb(fc, vrna_move_init(-i, pt[j]), VRNA_NEIGHBOR_INVALID, data);
          }
        }

        /* hop over base pair */
        j = pt[j];
      }
    }
    
    /* Finally, invalidate shift moves from pairs enclosed by
     * (move->pos_5, move->pos_3) into the interval outside the
     * interval [move->pos_5,move->pos_3]
     */
    for (j = move_i + 1; j < move_j; j++) {
      if (j < pt[j]) {
        for (i = enclosing_5 + 1; i < move_i; i++) {
          if (pt[i] > i) {
            i = pt[i];
          } else {
            if (is_compatible(fc, i, j))
              cb(fc, vrna_move_init(-i, j), VRNA_NEIGHBOR_INVALID, data);

            if (is_compatible(fc, i, pt[j]))
              cb(fc, vrna_move_init(-i, pt[j]), VRNA_NEIGHBOR_INVALID, data);
          }
        }

        for (i = move_j + 1; i < enclosing_3; i++) {
          if (i < pt[i]) {
            i = pt[i];
          } else {
            if (is_compatible(fc, j, i))
              cb(fc, vrna_move_init(j, -i), VRNA_NEIGHBOR_INVALID, data);

            if (is_compatible(fc, pt[j], i))
              cb(fc, vrna_move_init(pt[j], -i), VRNA_NEIGHBOR_INVALID, data);
          }
        }

        /* hop over base pair */
        j = pt[j];
      }
    }
  }
}


PRIVATE void
generate_conflicts_local_nb_deletion(vrna_fold_compound_t      *fc,
                                     const short               *pt,
                                     const vrna_move_t         *move,
                                     vrna_move_update_f cb,
                                     void                      *data,
                                     unsigned int              options)
{
  unsigned int  n = fc->length;
  int           i, j;
  int           enclosing_5   = 0;          /* exterior loop */
  int           enclosing_3   = (int)n + 1; /* exterior loop */
  int           min_loop_size = fc->params->model_details.min_loop_size;
  int           move_i        = -move->pos_5;
  int           move_j        = -move->pos_3;

  /* upon deletion of a base pair, conflicts can only
   * arise for shift moves of base pair (move->pos_5, move->pos_3)
   */
  if (options & VRNA_MOVESET_SHIFT) {
    /* find out actual enclosing pair that delimits the loop affected by the current move */
    for (i = move_i - 1; i > 0; i--) {
      if (pt[i] == 0) {
        continue;
      } else if (pt[i] < i) {
        i = pt[i]; /* hop over branching stems */
      } else if (pt[i] > i) {
        /* found enclosing pair */
        enclosing_5 = i;
        enclosing_3 = pt[i];
        break;
      }
    }

    /* moves to the region 5' of move->pos_5 */
    for (i = enclosing_5 + 1; i < move_i; i++) {
      if (pt[i] > i) {
        i = pt[i];
      } else {
        if (is_compatible(fc, i, move_i))
          cb(fc, vrna_move_init(-i, move_i), VRNA_NEIGHBOR_INVALID, data);
        if (is_compatible(fc, i, move_j))
          cb(fc, vrna_move_init(-i, move_j), VRNA_NEIGHBOR_INVALID, data);
      }
    }

    /* moves within [move->pos_5:move->pos_3] */
    for (i = move_i + 1; i < move_j; i++) {
      if (pt[i] > i) {
        i = pt[i];
      } else {
        if (is_compatible(fc, move_i, i))
          cb(fc, vrna_move_init(move_i, -i), VRNA_NEIGHBOR_INVALID, data);
        if (is_compatible(fc, i, move_j))
          cb(fc, vrna_move_init(-i, move_j), VRNA_NEIGHBOR_INVALID, data);
      }
    }

    /* moves to the region 3' of move->pos_3 */
    for (j = move_j + 1; j < enclosing_3; j++) {
      if (pt[j] > j) {
        j = pt[j];
      } else {
        if (is_compatible(fc, move_i, j))
          cb(fc, vrna_move_init(move_i, -j), VRNA_NEIGHBOR_INVALID, data);
        if (is_compatible(fc, move_j, j))
          cb(fc, vrna_move_init(move_j, -j), VRNA_NEIGHBOR_INVALID, data);
      }
    }
  }
}

PRIVATE void
generate_conflicts_local_nb_shift(vrna_fold_compound_t      *fc,
                                  const short               *pt,
                                  const vrna_move_t         *move,
                                  vrna_move_update_f cb,
                                  void                      *data,
                                  unsigned int              options)
{
  unsigned int  n = fc->length;
  int           enclosing_5   = 0;          /* exterior loop */
  int           enclosing_3   = (int)n + 1; /* exterior loop */
  int i, j, p,q, move_i, move_j;

  if (move->pos_5 < 0) {
    p = pt[move->pos_3];
    q = move->pos_3;
    move_i = -move->pos_5;
    move_j = move->pos_3;
  } else {
    p = move->pos_5;
    q = pt[move->pos_5];
    move_i = move->pos_5;
    move_j = -move->pos_3;
  }

  if (p > q) {
    i = p;
    p = q;
    q = i;
  }

  if (move_i > move_j) {
    i = move_i;
    move_i = move_j;
    move_j = i;
  }

  /* find out actual enclosing pair that delimits the loop affected by the current move */
  for (i = move_i - 1; i > 0; i--) {
    if (pt[i] == 0) {
      continue;
    } else if (pt[i] < i) {
      i = pt[i]; /* hop over branching stems */
    } else if (pt[i] > i) {
      /* found enclosing pair */
      enclosing_5 = i;
      enclosing_3 = pt[i];
      break;
    }
  }

  /* we now have the current base pair as (p, q) which will change
   * to (move_i, move_j) after application of the move
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
      for (i = enclosing_5 + 1; i < p; i++) {
        if (pt[i] > i) {
          i = pt[i];
        } else {
          if (is_compatible(fc, i, j))
            cb(fc, vrna_move_init(i, j), VRNA_NEIGHBOR_INVALID, data);
        }
      }

      for (i = q + 1; i < move_j; i++) {
        if (pt[i] > i) {
          i = pt[i];
        } else {
          if (is_compatible(fc, i, j))
            cb(fc, vrna_move_init(i, j), VRNA_NEIGHBOR_INVALID, data);
        }
      }

      for (i = move_j + 1; i < enclosing_3; i++) {
        if (pt[i] > i) {
          i = pt[i];
        } else {
          if (is_compatible(fc, j, i))
            cb(fc, vrna_move_init(j, i), VRNA_NEIGHBOR_INVALID, data);
        }
      }

      /* 1.2 invalidate insertion from the outside into interval [q + 1: move_j] */
      for (j = q + 1; j <= move_j; j++) {
        if (pt[j] > j) {
          j = pt[j];
        } else {
          for (i = enclosing_5; i < p; i++) {
            if (pt[i] > i) {
              i = pt[i];
            } else {
              if (is_compatible(fc, i, j))
                cb(fc, vrna_move_init(i, j), VRNA_NEIGHBOR_INVALID, data);
            }
          }

          for (i = move_j; i < enclosing_3; i++) {
            if (pt[i] > i) {
              i = pt[i];
            } else {
              if (is_compatible(fc, j, i))
                cb(fc, vrna_move_init(j, i), VRNA_NEIGHBOR_INVALID, data);
            }
          }
        }
      }
    } else if (p == move_j) {
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
      for (i = enclosing_5; i < move_i; i++) {
        if (pt[i] > i) {
          i = pt[i];
        } else {
          if (is_compatible(fc, i, j))
            cb(fc, vrna_move_init(i, j), VRNA_NEIGHBOR_INVALID, data);
        }
      }

      for (i = move_i + 1; i < p; i++) {
        if (pt[i] > i) {
          i = pt[i];
        } else {
          if (is_compatible(fc, j, i))
            cb(fc, vrna_move_init(j, i), VRNA_NEIGHBOR_INVALID, data);
        }
      }

      for (i = q + 1; i < enclosing_3; i++) {
        if (pt[i] > i) {
          i = pt[i];
        } else {
          if (is_compatible(fc, j, i))
            cb(fc, vrna_move_init(j, i), VRNA_NEIGHBOR_INVALID, data);
        }
      }

      /* 2.2 invalidate insertions where one pairing partner is within interval [q':p-1] */
      for (j = move_i; j < p; j++) {
        if (pt[j] > j) {
          j = pt[j];
        } else {
          for (i = enclosing_5 + 1; i <= move_i; i++) {
            if (pt[i] > i) {
              i = pt[i];
            } else {
              if (is_compatible(fc, i, j))
                cb(fc, vrna_move_init(i, j), VRNA_NEIGHBOR_INVALID, data);
            }
          }
          for (i = q + 1; i < enclosing_3; i++) {
            if (pt[i] > i) {
              i = pt[i];
            } else {
              if (is_compatible(fc, j, i))
                cb(fc, vrna_move_init(j, i), VRNA_NEIGHBOR_INVALID, data);
            }
          }
        }
      }
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

      /* 3.1 invalidate insertions where one pairing partner is p' */
      for (j = move_i; j < p; j++) {
        if (pt[j] > j) {
          j = pt[j];
        } else {
          for (i = enclosing_5 + 1; i <= move_i; i++) {
            if (pt[i] > i) {
              i = pt[i];
            } else {
              if (is_compatible(fc, i, j))
                cb(fc, vrna_move_init(i, j), VRNA_NEIGHBOR_INVALID, data);
            }
          }

          for (i = q + 1; i < enclosing_3; i++) {
            if (pt[i] > i) {
              i = pt[i];
            } else {
              if (is_compatible(fc, j, i))
                cb(fc, vrna_move_init(j, i), VRNA_NEIGHBOR_INVALID, data);
            }
          }
        }
      }
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

      /* 4.1 invalidate insertions where one pairing partner is q' */
      for (j = q + 1; j <= move_j; j++) {
        if (pt[j] > j) {
          j = pt[j];
        } else {
          for (i = enclosing_5 + 1; i < p; i++) {
            if (pt[i] > i) {
              i = pt[i];
            } else {
              if (is_compatible(fc, i, j))
                cb(fc, vrna_move_init(i, j), VRNA_NEIGHBOR_INVALID, data);
            }
          }

          for (i = move_j; i < enclosing_3; i++) {
            if (pt[i] > i) {
              i = pt[i];
            } else {
              if (is_compatible(fc, i, j))
                cb(fc, vrna_move_init(j, i), VRNA_NEIGHBOR_INVALID, data);
            }
          }
        }
      }
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
      for (i = p + 1; i <= move_i; i++) {
        if (pt[i] > i) {
          i = pt[i];
        } else {
          for (j = move_i; j < q; j++) {
            if (pt[j] > j) {
              j = pt[j];
            } else {
              if (is_compatible(fc, i, j))
                cb(fc, vrna_move_init(i, j), VRNA_NEIGHBOR_INVALID, data);
            }
          }
        }
      }
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
      for (i = p + 1; i <= move_j; i++) {
        if (pt[i] > i) {
          i = pt[i];
        } else {
          for (j = move_j; j < q; j++) {
            if (pt[j] > j) {
              j = pt[j];
            } else {
              if (is_compatible(fc, i, j))
                cb(fc, vrna_move_init(i, j), VRNA_NEIGHBOR_INVALID, data);
            }
          }
        }
      }
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
      for (i = enclosing_5 + 1; i < p; i++) {
        if (pt[i] > i) {
          i = pt[i];
        } else {
          if (is_compatible(fc, i, q))
            cb(fc, vrna_move_init(-i, q), VRNA_NEIGHBOR_INVALID, data);
        }
      }

      for (i = p + 1; i < q; i++) {
        if (pt[i] > i) {
          i = pt[i];
        } else {
          if (is_compatible(fc, i, q))
            cb(fc, vrna_move_init(-i, q), VRNA_NEIGHBOR_INVALID, data);
        }
      }

      for (i = q + 1; i < enclosing_3; i++) {
        if (pt[i] > i) {
          i = pt[i];
        } else {
          if (is_compatible(fc, q, i))
            cb(fc, vrna_move_init(q, -i), VRNA_NEIGHBOR_INVALID, data);
        }
      }

      /* 1.2 invalidate moves of enclosing pair positions, if any, into interval [move_i:p-1] */
      if (enclosing_5 > 0) {
        for (j = move_i; j < p; j++) {
          if (pt[j] > j) {
            j = pt[j];
          } else {
           if (is_compatible(fc, enclosing_5, j))
              cb(fc, vrna_move_init(enclosing_5, -j), VRNA_NEIGHBOR_INVALID, data);
            if (is_compatible(fc, j, enclosing_3))
              cb(fc, vrna_move_init(-j, enclosing_3), VRNA_NEIGHBOR_INVALID, data);
          }
        }
      }

      /* 1.3 invalidate moves of other pairs into interval [move_i:p-1] */
      for (j = move_i; j < p; j++) {
        if (pt[j] > j) {
          j = pt[j];
        } else {
          for (i = enclosing_5 + 1; i < move_i; i++) {
            if (pt[i] > i) {
              if (is_compatible(fc, i, j))
                cb(fc, vrna_move_init(i, -j), VRNA_NEIGHBOR_INVALID, data);
              if (is_compatible(fc, pt[i], j))
                cb(fc, vrna_move_init(pt[i], -j), VRNA_NEIGHBOR_INVALID, data);

              i = pt[i];
            }
          }

          for (i = q + 1; i < enclosing_3; i++) {
            if (pt[i] > i) {
              if (is_compatible(fc, j, i))
                cb(fc, vrna_move_init(-j, i), VRNA_NEIGHBOR_INVALID, data);
              if (is_compatible(fc, j, pt[i]))
                cb(fc, vrna_move_init(-j, pt[i]), VRNA_NEIGHBOR_INVALID, data);

              i = pt[i];
            }
          }
        }
      }

      /* 1.4 invalidate moves from within the interval [move_i+1:p-1] outside of the interval */
      for (j = move_i + 1; j < p; j++) {
        if (pt[j] > j) {
          for (i = enclosing_5 + 1; i < move_i; i++) {
            if (pt[i] > i) {
              i = pt[i];
            } else {
              if (is_compatible(fc, i, j))
                cb(fc, vrna_move_init(-i, j), VRNA_NEIGHBOR_INVALID, data);
              if (is_compatible(fc, i, pt[j]))
                cb(fc, vrna_move_init(-i, pt[j]), VRNA_NEIGHBOR_INVALID, data);
            }
          }

          for (i = q + 1; i < enclosing_3; i++) {
            if (pt[i] > i) {
              i = pt[i];
            } else {
              if (is_compatible(fc, j, i))
                cb(fc, vrna_move_init(j, -i), VRNA_NEIGHBOR_INVALID, data);
              if (is_compatible(fc, pt[j], i))
                cb(fc, vrna_move_init(pt[j], -i), VRNA_NEIGHBOR_INVALID, data);
            }
          }

          j = pt[j];
        }
      }
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
      for (i = enclosing_5 + 1; i < p; i++) {
        if (pt[i] > i) {
          i = pt[i];
        } else {
          if (is_compatible(fc, i, p))
            cb(fc, vrna_move_init(-i, p), VRNA_NEIGHBOR_INVALID, data);
        }
      }

      for (i = p + 1; i < q; i++) {
        if (pt[i] > i) {
          i = pt[i];
        } else {
          if (is_compatible(fc, p, i))
            cb(fc, vrna_move_init(p, -i), VRNA_NEIGHBOR_INVALID, data);
        }
      }

      for (i = q + 1; i < enclosing_3; i++) {
        if (pt[i] > i) {
          i = pt[i];
        } else {
          if (is_compatible(fc, p, i))
            cb(fc, vrna_move_init(p, -i), VRNA_NEIGHBOR_INVALID, data);
        }
      }

      /* 2.2 invalidate moves of enclosing pair, if any, into interval [q+1:move_j] */
      if (enclosing_5 > 0) {
        for (j = q + 1; j <= move_j; j++) {
          if (pt[j] > j) {
            j = pt[j];
          } else {
            if (is_compatible(fc, enclosing_5, j))
              cb(fc, vrna_move_init(enclosing_5, -j), VRNA_NEIGHBOR_INVALID, data);
            if (is_compatible(fc, j, enclosing_3))
              cb(fc, vrna_move_init(-j, enclosing_3), VRNA_NEIGHBOR_INVALID, data);
          }
        }
      }

      /* 2.3 invalidate moves of other pairs into interval [q+1:move_j] */
      for (j = q + 1; j <= move_j; j++) {
        if (pt[j] > j) {
          j = pt[j];
        } else {
          for (i = enclosing_5 + 1; i < p; i++) {
            if (pt[i] > i) {
              if (is_compatible(fc, i, j))
                cb(fc, vrna_move_init(i, -j), VRNA_NEIGHBOR_INVALID, data);
              if (is_compatible(fc, pt[i], j))
                cb(fc, vrna_move_init(pt[i], -j), VRNA_NEIGHBOR_INVALID, data);

              i = pt[i];
            }
          }

          for (i = move_j + 1; i < enclosing_3; i++) {
            if (pt[i] > i) {
              if (is_compatible(fc, j, i))
                cb(fc, vrna_move_init(-j, i), VRNA_NEIGHBOR_INVALID, data);
              if (is_compatible(fc, j, pt[i]))
                cb(fc, vrna_move_init(-j, pt[i]), VRNA_NEIGHBOR_INVALID, data);

              i = pt[i];
            }
          }
        }
      }

      /* 2.4 invalidate move of pairs within interval [q+1:move_j-1] outside of the interval */
      for (j = q + 1; j < move_j; j++) {
        if (pt[j] > j) {
          for (i = enclosing_5 + 1; i < p; i++) {
            if (pt[i] > i) {
              i = pt[i];
            } else {
              if (is_compatible(fc, i, j))
                cb(fc, vrna_move_init(-i, j), VRNA_NEIGHBOR_INVALID, data);
              if (is_compatible(fc, i, pt[j]))
                cb(fc, vrna_move_init(-i, pt[j]), VRNA_NEIGHBOR_INVALID, data);
            }
          }

          for (i = move_j + 1; i < enclosing_3; i++) {
            if (pt[i] > i) {
              i = pt[i];
            } else {
              if (is_compatible(fc, j, i))
                cb(fc, vrna_move_init(j, -i), VRNA_NEIGHBOR_INVALID, data);
              if (is_compatible(fc, pt[j], i))
                cb(fc, vrna_move_init(pt[j], -i), VRNA_NEIGHBOR_INVALID, data);
            }
          }

          j = pt[j];
        }
      }
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
      for (i = enclosing_5 + 1; i < p; i++) {
        if (pt[i] > i) {
          i = pt[i];
        } else {
          if (is_compatible(fc, i, p))
            cb(fc, vrna_move_init(-i, p), VRNA_NEIGHBOR_INVALID, data);
        }
      }

      for (i = p + 1; i < q; i++) {
        if (pt[i] > i) {
          i = pt[i];
        } else {
          if (is_compatible(fc, p, i))
            cb(fc, vrna_move_init(p, -i), VRNA_NEIGHBOR_INVALID, data);
        }
      }

      for (j = q + 1; j < enclosing_3; j++) {
        if (pt[j] > j) {
          j = pt[j];
        } else {
          if (is_compatible(fc, p, j))
            cb(fc, vrna_move_init(p, -j), VRNA_NEIGHBOR_INVALID, data);
        }
      }

      /* 3.2 invalidate moves of enclosing pair, if any, into interval [move_i:p-1] */
      if (enclosing_5 > 0) {
        for (i = move_i; i < p; i++) {
          if (pt[i] > i) {
            i = pt[i];
          } else {
            if (is_compatible(fc, enclosing_5, i))
              cb(fc, vrna_move_init(enclosing_5, -i), VRNA_NEIGHBOR_INVALID, data);
            if (is_compatible(fc, i, enclosing_3))
              cb(fc, vrna_move_init(-i, enclosing_3), VRNA_NEIGHBOR_INVALID, data);
          }
        }
      }

      /* 3.3 invalidate shift of any other pair into the interval [move_i:p-1] */
      for (j = move_i; j < p; j++) {
        if (pt[j] > j) {
          j = pt[j];
        } else {
          for (i = enclosing_5 + 1; i < move_i; i++) {
            if (pt[i] > i) {
              if (is_compatible(fc, i, j))
                cb(fc, vrna_move_init(i, -j), VRNA_NEIGHBOR_INVALID, data);
              if (is_compatible(fc, pt[i], j))
                cb(fc, vrna_move_init(pt[i], -j), VRNA_NEIGHBOR_INVALID, data);

              i = pt[i];
            }
          }

          for (i = q + 1; i < enclosing_3; i++) {
            if (pt[i] > i) {
              if (is_compatible(fc, j, i))
                cb(fc, vrna_move_init(-j, i), VRNA_NEIGHBOR_INVALID, data);
              if (is_compatible(fc, j, pt[i]))
                cb(fc, vrna_move_init(-j, pt[i]), VRNA_NEIGHBOR_INVALID, data);

              i = pt[i];
            }
          }
        }
      }

      /* 3.4 invalidate shift moves from pairs within interval to the outside */
      for (j = move_i + 1; j < p; j++) {
        if (pt[j] > j) {
          for (i = enclosing_5 + 1; i < move_i; i++) {
            if (pt[i] > i) {
              i = pt[i];
            } else {
              if (is_compatible(fc, i, j))
                cb(fc, vrna_move_init(-i, j), VRNA_NEIGHBOR_INVALID, data);
              if (is_compatible(fc, i, pt[j]))
                cb(fc, vrna_move_init(-i, pt[j]), VRNA_NEIGHBOR_INVALID, data);
            }
          }

          for (i = q + 1; i < enclosing_3; i++) {
            if (pt[i] > i) {
              i = pt[i];
            } else {
              if (is_compatible(fc, j, i))
                cb(fc, vrna_move_init(j, -i), VRNA_NEIGHBOR_INVALID, data);
              if (is_compatible(fc, pt[j], i))
                cb(fc, vrna_move_init(pt[j], -i), VRNA_NEIGHBOR_INVALID, data);
            }
          }

          j = pt[j];
        }
      }
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
      for (i = enclosing_5 + 1; i < p; i++) {
        if (pt[i] > i) {
          i = pt[i];
        } else {
          if (is_compatible(fc, i, q))
            cb(fc, vrna_move_init(-i, q), VRNA_NEIGHBOR_INVALID, data);
        }
      }

      for (i = p + 1; i < q; i++) {
        if (pt[i] > i) {
          i = pt[i];
        } else {
          if (is_compatible(fc, i, q))
            cb(fc, vrna_move_init(-i, q), VRNA_NEIGHBOR_INVALID, data);
        }
      }

      for (i = q + 1; i < enclosing_3; i++) {
        if (pt[i] > i) {
          i = pt[i];
        } else {
          if (is_compatible(fc, q, i))
            cb(fc, vrna_move_init(q, -i), VRNA_NEIGHBOR_INVALID, data);
        }
      }

      /* 4.2 invalidate moves of enclosing pair, if any, into interval [q+1:move_j] */
      if (enclosing_5 > 0) {
        for (j = q + 1; j <= move_j; j++) {
          if (pt[j] > j) {
            j = pt[j];
          } else {
            if (is_compatible(fc, enclosing_5, j))
              cb(fc, vrna_move_init(enclosing_5, -j), VRNA_NEIGHBOR_INVALID, data);
            if (is_compatible(fc, j, enclosing_3))
              cb(fc, vrna_move_init(-j, enclosing_3), VRNA_NEIGHBOR_INVALID, data);
          }
        }
      }

      /* 4.3 invalidate shifts from other pairs into the interval [q+1: move_j] */
      for (j = q + 1; j <= move_j; j++) {
        if (pt[j] > j) {
          j = pt[j];
        } else {
          for (i = enclosing_5 + 1; i < p; i++) {
            if (pt[i] > i) {
              if (is_compatible(fc, i, j))
                cb(fc, vrna_move_init(i, -j), VRNA_NEIGHBOR_INVALID, data);
              if (is_compatible(fc, pt[i], j))
                cb(fc, vrna_move_init(pt[i], -j), VRNA_NEIGHBOR_INVALID, data);

              i = pt[i];
            }
          }
          for (i = move_j + 1; i < enclosing_3; i++) {
            if (pt[i] > i) {
              if (is_compatible(fc, j, i))
                cb(fc, vrna_move_init(-j, i), VRNA_NEIGHBOR_INVALID, data);
              if (is_compatible(fc, j, pt[i]))
                cb(fc, vrna_move_init(-j, pt[i]), VRNA_NEIGHBOR_INVALID, data);

              i = pt[i];
            }
          }
        }
      }

      /* 4.4 invalidate shifts from within the interval [q + 1: move_j - 1] to the outside */
      for (j = q + 1; j < move_j; j++) {
        if (pt[j] > j) {
          for (i = enclosing_5; i < p; i++) {
            if (pt[i] > i) {
              i = pt[i];
            } else {
              if (is_compatible(fc, i, j))
                cb(fc, vrna_move_init(-i, j), VRNA_NEIGHBOR_INVALID, data);
              if (is_compatible(fc, i, pt[j]))
                cb(fc, vrna_move_init(-i, pt[j]), VRNA_NEIGHBOR_INVALID, data);
            }
          }

          for (i = move_j + 1; i < enclosing_3; i++) {
            if (pt[i] > i) {
              i = pt[i];
            } else {
              if (is_compatible(fc, j, i))
                cb(fc, vrna_move_init(j, -i), VRNA_NEIGHBOR_INVALID, data);
              if (is_compatible(fc, pt[j], i))
                cb(fc, vrna_move_init(pt[j], -i), VRNA_NEIGHBOR_INVALID, data);
            }
          }

          j = pt[j];
        }
      }
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
      for (i = enclosing_5 + 1; i < p; i++) {
        if (pt[i] > i) {
          i = pt[i];
        } else {
          if (is_compatible(fc, i, p))
            cb(fc, vrna_move_init(-i, p), VRNA_NEIGHBOR_INVALID, data);
        }
      }

      for (i = p + 1; i < q; i++) {
        if (pt[i] > i) {
          i = pt[i];
        } else {
          if (is_compatible(fc, p, i))
            cb(fc, vrna_move_init(p, -i), VRNA_NEIGHBOR_INVALID, data);
        }
      }

      for (i = q + 1; i < enclosing_3; i++) {
        if (pt[i] > i) {
          i = pt[i];
        } else {
          if (is_compatible(fc, p, i))
            cb(fc, vrna_move_init(p, -i), VRNA_NEIGHBOR_INVALID, data);
        }
      }

      /* 5.2 invalidate moves of pairs into interval [p + 1:move_i] */
      for (i = p + 1; i <= move_i; i++) {
        if (pt[i] > i) {
          i = pt[i];
        } else {
          for (j = move_i + 1; j < q; j++) {
            if (pt[j] > j) {
              if (is_compatible(fc, i, j))
                cb(fc, vrna_move_init(-i, j), VRNA_NEIGHBOR_INVALID, data);
              if (is_compatible(fc, i, pt[j]))
                cb(fc, vrna_move_init(-i, pt[j]), VRNA_NEIGHBOR_INVALID, data);
              j = pt[j];
            }
          }
        }
      }

      /* 5.3 invalidate moves of pairs from within interval [p + 1: move_i] */
      for (i = p + 1; i < move_i; i++) {
        if (pt[i] > i) {
          for (j = move_i + 1; j < q; j++) {
            if (pt[j] > j) {
              j = pt[j];
            } else {
              if (is_compatible(fc, i, j))
                cb(fc, vrna_move_init(i, -j), VRNA_NEIGHBOR_INVALID, data);
              if (is_compatible(fc, pt[i], j))
                cb(fc, vrna_move_init(pt[i], -j), VRNA_NEIGHBOR_INVALID, data);
            }
          }
          i = pt[i];
        }
      }
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
      for (i = enclosing_5; i < p; i++) {
        if (pt[i] > i) {
          i = pt[i];
        } else {
          if (is_compatible(fc, i, q))
            cb(fc, vrna_move_init(-i, q), VRNA_NEIGHBOR_INVALID, data);
        }
      }

      for (i = p + 1; i < q; i++) {
        if (pt[i] > i) {
          i = pt[i];
        } else {
          if (is_compatible(fc, i, q))
            cb(fc, vrna_move_init(-i, q), VRNA_NEIGHBOR_INVALID, data);
        }
      }

      for (i = q + 1; i < enclosing_3; i++) {
        if (pt[i] > i) {
          i = pt[i];
        } else {
          if (is_compatible(fc, q, i))
            cb(fc, vrna_move_init(q, -i), VRNA_NEIGHBOR_INVALID, data);
        }
      }

      /* 6.2 invalidate shifts of pairs into interval [move_j:q-1] */
      for (j = move_j; j < q; j++) {
        if (pt[j] > j) {
          j = pt[j];
        } else {
          for (i = p + 1; i < move_j; i++) {
            if (pt[i] > i) {
              if (is_compatible(fc, i, j))
                cb(fc, vrna_move_init(i, -j), VRNA_NEIGHBOR_INVALID, data);
              if (is_compatible(fc, pt[i], j))
                cb(fc, vrna_move_init(pt[i], -j), VRNA_NEIGHBOR_INVALID, data);

              i = pt[i];
            }
          }
        }
      }

      /* 6.3 invalidate shifts from within interval [move_j+1:q-1] */
      for (j = move_j + 1; j < q; j++) {
        if (pt[j] > j) {
          for (i = p + 1; i < move_j; i++) {
            if (pt[i] > i) {
              i = pt[i];
            } else {
              if (is_compatible(fc, i, j))
                cb(fc, vrna_move_init(-i, j), VRNA_NEIGHBOR_INVALID, data);
              if (is_compatible(fc, i, pt[j]))
                cb(fc, vrna_move_init(-i, pt[j]), VRNA_NEIGHBOR_INVALID, data);
            }
          }

          j = pt[j];
        }
      }
    }
  }
}
