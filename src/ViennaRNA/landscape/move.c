#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "ViennaRNA/pair_mat.h"
#include "ViennaRNA/io/utils.h"
#include "ViennaRNA/utils/structures.h"
#include "ViennaRNA/model.h"
#include "ViennaRNA/eval.h"
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
  if (m->pos_5 < 0 && m->pos_3 < 0) {
    pt[-m->pos_5] = 0;
    pt[-m->pos_3] = 0;
  } else

  /* insertion */
  if (m->pos_5 > 0 && m->pos_3 > 0) {
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
