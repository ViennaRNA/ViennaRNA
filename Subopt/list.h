/*
  $Log: list.h,v $
  Revision 1.1  1997/08/04 21:05:32  walter
  Initial revision

*/

#ifndef	__LIST_H
#define	__LIST_H

#include "include.h"

/*---------------------- Macros and type definitions ----------------------*/

typedef struct LST_BUCKET
  {
    struct LST_BUCKET *next;
  }
LST_BUCKET;

typedef struct
  {
    int count;			/* Number of elements currently in list */
    LST_BUCKET *head;		/* Pointer to head element of list      */
    LST_BUCKET *z;		/* Pointer to last node of list         */
    LST_BUCKET hz[2];		/* Space for head and z nodes           */
  }
LIST;

/* Return a pointer to the user space given the address of the header of
 * a node.
 */

#define	LST_USERSPACE(h)	((void*)((LST_BUCKET*)(h) + 1))

/* Return a pointer to the header of a node, given the address of the
 * user space.
 */

#define	LST_HEADER(n)		((LST_BUCKET*)(n) - 1)

/* Return a pointer to the user space of the list's head node. This user
 * space does not actually exist, but it is useful to be able to address
 * it to enable insertion at the start of the list.
 */

#define	LST_HEAD(l)		LST_USERSPACE((l)->head)

/* Determine if a list is empty
 */

#define	LST_EMPTY(l)		((l)->count == 0)

/*-------------------------- Function Prototypes --------------------------*/

void *lst_newnode (int size);
void lst_freenode (void *node);
LIST *lst_init (void);
void lst_kill (LIST * l, void (*freeNode) ());
void lst_insertafter (LIST * l, void *node, void *after);
void *lst_deletenext (LIST * l, void *node);
void *lst_first (LIST * l);
void *lst_next (void *prev);
void lst_mergesort (LIST * l, int (*cmp_func) ());

#endif
