/*
  $Log: list.c,v $
  Revision 1.5  2003/07/14 13:36:58  ivo
  use vrna_alloc() instead of malloc

  Revision 1.4  2000/10/10 08:53:52  ivo
  include dmalloc.h header if DMALLOC defined

  Revision 1.4  2000/10/10 08:04:34  ivo
  include dmalloc header id DMALLOC defined

  Revision 1.3  1998/03/30 14:24:51  ivo
  use RNA package utils.h

  Revision 1.2  1997/10/09  19:01:50  steve
  *** empty log message ***

  Revision 1.1  1997/08/04 21:05:32  walter
  Initial revision

*/
/*
   (C) 1991 Kendall Bennett.
*/

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#include "ViennaRNA/utils.h"
#include "ViennaRNA/list.h"

#define PUBLIC
PUBLIC void *
lst_newnode (int size)
/****************************************************************************
*
* Function:	lst_newnode
* Parameters:	size - Amount of memory to allocate for node
* Returns:      Pointer to the allocated node's user space.
*
* Description:	Allocates the memory required for a node, adding a small
*		header at the start of the node. We return a reference to
*		the user space of the node, as if it had been allocated via
*		malloc().
*
****************************************************************************/
{
  LST_BUCKET *node;

  node = (LST_BUCKET *) vrna_alloc(size + sizeof (LST_BUCKET));

  return LST_USERSPACE (node);	/* Return pointer to user space */
}

PUBLIC void
lst_freenode (void *node)
/****************************************************************************
*
* Function:	lst_freenode
* Parameters:	node - Node to free.
*
* Description:  Frees a node previously allocated with lst_newnode().
*
****************************************************************************/
{
  free (LST_HEADER (node));
}

PUBLIC LIST *
lst_init (void)
/****************************************************************************
*
* Function:	lst_init
* Returns:      Pointer to a newly created list.
*
* Description:	Initialises a list and returns a pointer to it.
*
****************************************************************************/
{
  LIST *l;

  if ((l = (LIST *) vrna_alloc(sizeof (LIST))) != NULL)
    {
      l->count = 0;
      l->head = &(l->hz[0]);
      l->z = &(l->hz[1]);
      l->head->next = l->z->next = l->z;
    }

  return l;
}

PUBLIC void
lst_kill (LIST * l, void (*freeNode) (void *node))
/****************************************************************************
*
* Function:	lst_kill
* Parameters:	l - List to kill
*		freeNode - Pointer to user routine to free a node
*
* Description:	Kills the list l, by deleting all of the elements contained
*		within the list one by one and then deleting the list
*		itself. Note that we call the user supplied routine
*		(*freeNode)() to free each list node. This allows the user
*		program to perform any extra processing needed to kill each
*		node (if each node contains pointers to other items on the
*		heap for example). If no extra processing is required, just
*		pass the address of lst_freenode(), ie:
*
*		lst_kill(myList,lst_freenode);
*
****************************************************************************/
{
  LST_BUCKET *n, *p;

  n = l->head->next;
  while (n != l->z)
    {				/* Free all nodes in list  */
      p = n;
      n = n->next;
      (*freeNode) (LST_USERSPACE (p));
    }
  free (l);			/* Free the list itself    */
}

PUBLIC void
lst_insertafter (LIST * l, void *node, void *after)
/****************************************************************************
*
* Function:	lst_insertafter
* Parameters:	l - List to insert node into
*		node - Pointer to user space of node to insert
*		after - Pointer to user space of node to insert node after
*
* Description:	Inserts a new node into the list after the node 'after'. To
*		insert a new node at the beginning of the list, user the
*		macro LST_HEAD in place of 'after'. ie:
*
*		lst_insertafter(mylist,node,LST_HEAD(mylist));
*
****************************************************************************/
{
  LST_BUCKET *n = LST_HEADER (node), *a = LST_HEADER (after);

  n->next = a->next;
  a->next = n;
  l->count++;
}

PUBLIC void *
lst_deletenext (LIST * l, void *node)
/****************************************************************************
*
* Function:	lst_deletenext
* Parameters:	l - List to delete node from.
*		node - Node to delete the next node from
* Returns:	Pointer to the deleted node's userspace.
*
* Description:	Removes the node AFTER 'node' from the list l.
*
****************************************************************************/
{
  LST_BUCKET *n = LST_HEADER (node);

  node = LST_USERSPACE (n->next);
  n->next = n->next->next;
  l->count--;
  return node;
}

PUBLIC void *
lst_first (LIST * l)
/****************************************************************************
*
* Function:	lst_first
* Parameters:	l - List to obtain first node from
* Returns:	Pointer to first node in list, NULL if list is empty.
*
* Description:	Returns a pointer to the user space of the first node in
*		the list. If the list is empty, we return NULL.
*
****************************************************************************/
{
  LST_BUCKET *n;

  n = l->head->next;
  return (n == l->z ? NULL : LST_USERSPACE (n));
}

PUBLIC void *
lst_next (void *prev)
/****************************************************************************
*
* Function:	lst_next
* Parameters:	prev - Previous node in list to obtain next node from
* Returns:	Pointer to the next node in the list, NULL at end of list.
*
* Description:	Returns a pointer to the user space of the next node in the
*		list given a pointer to the user space of the previous node.
*		If we have reached the end of the list, we return NULL. The
*		end of the list is detected when the next pointer of a node
*		points back to itself, as does the dummy last node's next
*		pointer. This enables us to detect the end of the list
*		without needed access to the list data structure itself.
*
*		NOTE:	We do no checking to ensure that 'prev' is NOT a
*			NULL pointer.
*
****************************************************************************/
{
  LST_BUCKET *n = LST_HEADER (prev);

  n = n->next;
  return (n == n->next ? NULL : LST_USERSPACE (n));
}

/* Static globals required by merge()   */

static LST_BUCKET *z;
static int (*cmp) (void *, void *);

static LST_BUCKET *
merge (LST_BUCKET * a, LST_BUCKET * b, LST_BUCKET ** end)
/****************************************************************************
*
* Function:	merge
* Parameters:	a,b - Sublist's to merge
* Returns:	Pointer to the merged sublists.
*
* Description:	Merges two sorted lists of nodes together into a single
*		sorted list.
*
****************************************************************************/
{
  LST_BUCKET *c;

  /* Go through the lists, merging them together in sorted order  */

  c = z;
  while (a != z && b != z)
    {
      if ((*cmp) (LST_USERSPACE (a), LST_USERSPACE (b)) <= 0)
	{
	  c->next = a;
	  c = a;
	  a = a->next;
	}
      else
	{
	  c->next = b;
	  c = b;
	  b = b->next;
	}
    };

  /* If one of the lists is not exhausted, then re-attach it to the end
   * of the newly merged list
   */

  if (a != z)
    c->next = a;
  if (b != z)
    c->next = b;

  /* Set *end to point to the end of the newly merged list        */

  while (c->next != z)
    c = c->next;
  *end = c;

  /* Determine the start of the merged lists, and reset z to point to
   * itself
   */

  c = z->next;
  z->next = z;
  return c;
}

PUBLIC void
lst_mergesort (LIST * l, int (*cmp_func) (void *, void *))
/****************************************************************************
*
* Function:	lst_mergesort
* Parameters:	l - List to merge sort
*		cmp_func - Function to compare two user spaces
*
* Description:	Mergesort's all the nodes in the list. 'cmp' must point to
*		a comparison function that can compare the user spaces of
*		two different nodes. 'cmp' should work the same as
*		strcmp(), in terms of the values it returns.
*
****************************************************************************/
{
  int i, N;
  LST_BUCKET *a, *b;		/* Pointers to sublists to merge                */
  LST_BUCKET *c;		/* Pointer to end of sorted sublists            */
  LST_BUCKET *head;		/* Pointer to dummy head node for list          */
  LST_BUCKET *todo;		/* Pointer to sublists yet to be sorted         */
  LST_BUCKET *t;		/* Temporary                                                            */

  /* Set up globals required by merge() and pointer to head       */

  z = l->z;
  cmp = cmp_func;
  head = l->head;

  for (N = 1, a = z; a != head->next; N = N + N)
    {
      todo = head->next;
      c = head;
      while (todo != z)
	{

	  /* Build first sublist to be merged, and splice from main list
	   */

	  a = t = todo;
	  for (i = 1; i < N; i++)
	    t = t->next;
	  b = t->next;
	  t->next = z;
	  t = b;

	  /* Build second sublist to be merged and splice from main list
	   */

	  for (i = 1; i < N; i++)
	    t = t->next;
	  todo = t->next;
	  t->next = z;

	  /* Merge the two sublists created, and set 'c' to point to the
	   * end of the newly merged sublists.
	   */

	  c->next = merge (a, b, &t);
	  c = t;
	}
    }
}

#ifdef LIST_TEST

/*---------------------------------------------------------------*/
/*---------------------------------------------------------------*/

/* Simple program to test the list routines */

typedef struct
{
  char name[40];
  int age;
}
REC;

/*---------------------------------------------------------------*/

int
my_cmp (REC * r1, REC * r2)
{
  return strcmp (r1->name, r2->name);
}

/*---------------------------------------------------------------*/

void
main (void)
{
  LIST *list;
  int done = 0;
  REC *rec;
  char line[80];

  list = lst_init ();

  printf ("Type a list of names and ages. Empty line quits\n\n");

  while (!done)
    {
      rec = lst_newnode (sizeof (REC));
      gets (line);
      if ((done = (line[0] == '\0')) != 1)
	{
	  strcpy (rec->name, line);
	  gets (line);
	  rec->age = atoi (line);
	  lst_insertafter (list, rec, LST_HEAD (list));
	}
    };

  printf ("\nThe list you typed in was:\n\n");

  for (rec = lst_first (list); rec; rec = lst_next (rec))
    printf ("Name: %s, Age: %d\n", rec->name, rec->age);

  printf ("\nSorting the list...\n\n");

  lst_mergesort (list, my_cmp);

  for (rec = lst_first (list); rec; rec = lst_next (rec))
    printf ("Name: %s, Age: %d\n", rec->name, rec->age);

  lst_kill (list, lst_freenode);
}

/*---------------------------------------------------------------*/

#endif
