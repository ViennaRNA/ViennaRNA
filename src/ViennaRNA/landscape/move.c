#include <stdio.h>
#include <stdlib.h>

#include "ViennaRNA/utils/basic.h"
#include "ViennaRNA/landscape/move.h"


PUBLIC vrna_move_t
vrna_move_init(int  pos_5,
               int  pos_3)
{
  return (vrna_move_t){
           pos_5, pos_3, NULL
  };
}


PUBLIC void
vrna_move_list_free(vrna_move_t *moves)
{
  if (moves) {
    for (vrna_move_t *move = moves; move->pos_5 != 0; move++) {
      if (move->next != NULL)
        if (move->next->pos_5 != 0)
          vrna_move_list_free(move->next);
    }
    free(moves);
  }
}


PUBLIC void
vrna_move_apply(short             *pt,
                const vrna_move_t *m)
{
  /* deletion */
  if (vrna_move_is_removal(m)) {
    pt[-m->pos_5] = 0;
    pt[-m->pos_3] = 0;
  } else if (vrna_move_is_insertion(m)) {
    pt[m->pos_5]  = m->pos_3;
    pt[m->pos_3]  = m->pos_5;
  } else
  /* shift right */
  if (m->pos_5 > 0 && m->pos_3 < 0) {
    short previousPairedPosition = pt[m->pos_5];
    pt[previousPairedPosition] = 0;
    short newPairedPosition = -m->pos_3;
    pt[m->pos_5]          = newPairedPosition;
    pt[newPairedPosition] = m->pos_5;
  } else

  /* shift left */
  if (m->pos_5 < 0 && m->pos_3 > 0) {
    short previousPairedPosition = pt[m->pos_3];
    pt[previousPairedPosition] = 0;
    short newPairedPosition = -m->pos_5;
    pt[m->pos_3]          = newPairedPosition;
    pt[newPairedPosition] = m->pos_3;
  }

  /* apply successive moves if m.next is a list */
  if (m->next != NULL)
    for (vrna_move_t *move = m->next; move->pos_5 != 0; move++)
      vrna_move_apply(pt, move);
}


PUBLIC void
vrna_move_apply_db(char               *structure,
                   const short        *pt,
                   const vrna_move_t  *m)
{
  /* deletion */
  if (vrna_move_is_removal(m)) {
    structure[-m->pos_5 - 1]  = '.';
    structure[-m->pos_3 - 1]  = '.';
    return;
  }

  /* insertion */
  if (vrna_move_is_insertion(m)) {
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


PUBLIC int
vrna_move_is_removal(const vrna_move_t *m)
{
  return (m->pos_5 < 0) && (m->pos_3 < 0);
}


PUBLIC int
vrna_move_is_insertion(const vrna_move_t *m)
{
  return (m->pos_5 > 0) && (m->pos_3 > 0);
}


PUBLIC int
vrna_move_is_shift(const vrna_move_t *m)
{
  return ((m->pos_5 < 0) && (m->pos_3 > 0)) ||
         ((m->pos_5 > 0) && (m->pos_3 < 0));
}


PUBLIC int
vrna_move_compare(const vrna_move_t *a,
                  const vrna_move_t *b,
                  const short       *ptable)
{
  /* assume, both moves a and b are compatible with current structure */

  if (vrna_move_is_removal(a)) {
    if (vrna_move_is_removal(b)) {
      if (a->pos_5 > b->pos_5)
        return 1;
      else if (a->pos_5 < b->pos_5)
        return -1;
      else
        return 0;
    } else if (vrna_move_is_insertion(b)) {
      return 1;
    } else {
      /*
       * 'b' is shift move
       * Implement me!
       */
    }
  } else if (vrna_move_is_insertion(a)) {
    if (vrna_move_is_insertion(b)) {
      if (a->pos_5 < b->pos_5)
        return -1;
      else if (a->pos_5 > b->pos_5)
        return 1;
      else if (a->pos_3 < b->pos_3)
        return -1;
      else if (a->pos_3 > b->pos_3)
        return 1;
      else
        return 0;
    } else if (vrna_move_is_removal(b)) {
      return -1;
    } else {
      /*
       * 'b' is shift move
       * Implement me!
       */
    }
  } else {
    /*
     * 'a' is shift move
     * Implement me!
     */
  }

  return 0;
}
