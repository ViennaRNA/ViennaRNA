/*
  $Log: list.h,v $
  Revision 1.2  2000/10/10 08:50:01  ivo
  some annotation for lclint

  Revision 1.1  1997/08/04 21:05:32  walter
  Initial revision

*/

#ifndef	__LIST_H
#define	__LIST_H

/*---------------------- Macros and type definitions ----------------------*/

typedef struct LST_BUCKET {
  struct LST_BUCKET *next;
}
LST_BUCKET;

typedef struct {
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

/*@only@*//*@out@*/ void *lst_newnode (int size);
void lst_freenode (/*@only@*/ void *node);
/*@only@*//*@out@*/  LIST *lst_init (void);
void lst_kill (LIST * l, void (*freeNode) ());
void lst_insertafter (LIST * l, /*@keep@*/ void *node, void *after);
void *lst_deletenext (/*@only@*/ LIST * l, void *node);
/*@dependent@*/ void *lst_first (LIST * l);
/*@dependent@*/ void *lst_next (void *prev);
void lst_mergesort (LIST * l, int (*cmp_func) ());

#endif
