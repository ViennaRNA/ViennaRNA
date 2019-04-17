/*
 * This is a default callback implementation for local neighborhood retrieval
 * after (successive) move applications
 */

struct movelist {
  vrna_move_t   *moves;
  vrna_move_t   *moves_invalid;
  unsigned int  num_moves;
  unsigned int  num_moves_invalid;
  unsigned int  mem_moves;
  unsigned int  mem_moves_invalid;
};


PRIVATE void
add_to_incremental_move_list(vrna_fold_compound_t *fc,
                             vrna_move_t          next_neighbor,
                             unsigned int         state,
                             void                 *data)
{
  /* add move to list */
  struct movelist *mlist = (struct movelist *)data;

  if ((state == VRNA_NEIGHBOR_CHANGE) || (state == VRNA_NEIGHBOR_NEW)) {
    mlist->moves[mlist->num_moves] = next_neighbor;
    mlist->num_moves++;

    /* check whether we must increase memory block for deletion moves */
    if (mlist->num_moves == mlist->mem_moves) {
      mlist->mem_moves  *= 1.4;
      mlist->moves      = (vrna_move_t *)vrna_realloc(mlist->moves,
                                                      sizeof(vrna_move_t) *
                                                      mlist->mem_moves);
    }
  } else if (state == VRNA_NEIGHBOR_INVALID) {
    mlist->moves_invalid[mlist->num_moves_invalid] = next_neighbor;
    mlist->num_moves_invalid++;

    /* check whether we must increase memory block for deletion moves */
    if (mlist->num_moves_invalid == mlist->mem_moves_invalid) {
      mlist->mem_moves_invalid  *= 1.4;
      mlist->moves_invalid      = (vrna_move_t *)vrna_realloc(mlist->moves_invalid,
                                                              sizeof(vrna_move_t) *
                                                              mlist->mem_moves_invalid);
    }
  }
}


PRIVATE struct movelist *
init_incremental_movelist(unsigned int n)
{
  struct movelist *mlist = (struct movelist *)vrna_alloc(sizeof(struct movelist));

  mlist->num_moves  = 0;
  mlist->mem_moves  = n;
  mlist->moves      = (vrna_move_t *)vrna_alloc(sizeof(vrna_move_t) *
                                                mlist->mem_moves);
  mlist->num_moves_invalid  = 0;
  mlist->mem_moves_invalid  = n;
  mlist->moves_invalid      = (vrna_move_t *)vrna_alloc(sizeof(vrna_move_t) *
                                                        mlist->mem_moves_invalid);

  return mlist;
}


PRIVATE void
free_incremental_movelist(struct movelist *m)
{
  free(m->moves);
  free(m->moves_invalid);
  free(m);
}
