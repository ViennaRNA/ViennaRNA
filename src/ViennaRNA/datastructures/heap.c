/*
 * This is a generic implementation of a binaray heap
 */
#include <stdlib.h>
#include <math.h>

#include "ViennaRNA/utils/basic.h"

#include  "ViennaRNA/datastructures/heap.h"

#ifdef __GNUC__
# define INLINE inline
#else
# define INLINE
#endif


struct vrna_heap_s {
  size_t                      num_entries;    /* number of entries in the heap */
  size_t                      mem_entries;    /* memory allocated for the entries */
  void                        **entries;      /* the array of entries */

  vrna_heap_cmp_f      cmp;           /* the compare function */
  vrna_heap_get_pos_f  get_entry_pos; /* position look-up function */
  vrna_heap_set_pos_f  set_entry_pos; /* position store function */

  void                        *data;          /* arbitrary user-data pointer */
};


/*
 #################################
 # PRIVATE FUNCTION DECLARATIONS #
 #################################
 */
PRIVATE INLINE unsigned int
parent_pos(unsigned int i);


PRIVATE INLINE unsigned int
child0_pos(unsigned int i);


PRIVATE INLINE unsigned int
child1_pos(unsigned int i);


PRIVATE INLINE void
swap_entries(struct vrna_heap_s *h,
             size_t             a,
             size_t             b);


PRIVATE INLINE int
heapify_up(struct vrna_heap_s *h,
           size_t             i);


PRIVATE INLINE void
heapify_down(struct vrna_heap_s *h,
             size_t             pos);


/*
 #################################
 # BEGIN OF FUNCTION DEFINITIONS #
 #################################
 */
PUBLIC struct vrna_heap_s *
vrna_heap_init(size_t                     n,
               vrna_heap_cmp_f     cmp,
               vrna_heap_get_pos_f get_entry_pos,
               vrna_heap_set_pos_f set_entry_pos,
               void                       *data)
{
  struct vrna_heap_s *h = NULL;

  if (cmp) {
    h                 = (struct vrna_heap_s *)vrna_alloc(sizeof(struct vrna_heap_s));
    h->num_entries    = 0;
    h->mem_entries    = n + 1;
    h->get_entry_pos  = NULL;
    h->set_entry_pos  = NULL;
    h->entries        = (void **)vrna_alloc(sizeof(void *) * (n + 1));
    h->cmp            = cmp;
    h->data           = data;

    /* only assign getter/setter functions if both are provided */
    if ((get_entry_pos) && (set_entry_pos)) {
      h->get_entry_pos  = get_entry_pos;
      h->set_entry_pos  = set_entry_pos;
    }
  }

  return h;
}


PUBLIC void
vrna_heap_free(struct vrna_heap_s *h)
{
  if (h) {
    free(h->entries);
    free(h);
  }
}


PUBLIC void
vrna_heap_insert(struct vrna_heap_s *h,
                 void               *v)
{
  size_t n;

  if ((h) && (v)) {
    n = ++h->num_entries;

    /* increase array of entries if necessary */
    if (n == h->mem_entries) {
      h->mem_entries  *= 1.4;
      h->entries      = (void **)vrna_realloc(h->entries, sizeof(void *) * h->mem_entries);
    }

    h->entries[n] = v;

    /* notify external storage about current position */
    if (h->set_entry_pos)
      h->set_entry_pos(v, n, h->data);

    heapify_up(h, n);
  }
}


PUBLIC const void *
vrna_heap_top(struct vrna_heap_s *h)
{
  if ((h) && (h->num_entries > 0))
    return (const void *)h->entries[1];

  return NULL;
}


PUBLIC void *
vrna_heap_pop(struct vrna_heap_s *h)
{
  if ((h) && (h->num_entries > 0)) {
    void *entry = h->entries[1];


    /* notify external storage about deletion */
    if (h->set_entry_pos)
      h->set_entry_pos(entry, 0, h->data);

    h->num_entries--;

    /* restore heap condition if there are entries left in the heap */
    if (h->num_entries) {
      swap_entries(h, 1, h->num_entries + 1);

      heapify_down(h, 1);
    }

    return entry;
  }

  return NULL;
}


PUBLIC void *
vrna_heap_remove(struct vrna_heap_s *h,
                 const void         *v)
{
  if ((h) && (h->get_entry_pos)) {
    /* get position of last entry in heap */
    size_t  last_pos = h->num_entries;

    /* obtain position of element to remove */
    size_t  pos = h->get_entry_pos(v, h->data);

    if (!pos)
      return NULL;

    void    *entry = h->entries[pos];

    /* notify external storage about deletion */
    h->set_entry_pos(v, 0, h->data);

    h->num_entries--;

    /* we only need to do anything if we didn't remove the last element */
    if (pos != last_pos) {
      h->entries[pos] = h->entries[last_pos];

      h->set_entry_pos(h->entries[pos], pos, h->data);

      if (!heapify_up(h, pos))
        heapify_down(h, pos);
    }

    return entry;
  }

  return NULL;
}


PUBLIC void *
vrna_heap_update(struct vrna_heap_s *h,
                 void               *v)
{
  if ((h) && (v) && (h->get_entry_pos)) {
    size_t pos = h->get_entry_pos(v, h->data);
    if (pos) {
      void  *old_entry = h->entries[pos];
      h->entries[pos] = v;

      int   cmp = h->cmp(v, old_entry, h->data);

      if (cmp < 0)
        heapify_up(h, pos);
      else if (cmp > 0)
        heapify_down(h, pos);

      return old_entry;
    } else {
      vrna_heap_insert(h, v);
    }
  }

  return NULL;
}


PUBLIC size_t
vrna_heap_size(struct vrna_heap_s *h)
{
  if (h)
    return h->num_entries;

  return 0;
}


/*
 #####################################
 # BEGIN OF STATIC HELPER FUNCTIONS  #
 #####################################
 */
PRIVATE INLINE unsigned int
parent_pos(unsigned int i)
{
  return floor(i / 2);
}


PRIVATE INLINE unsigned int
child0_pos(unsigned int i)
{
  return 2 * i;
}


PRIVATE INLINE unsigned int
child1_pos(unsigned int i)
{
  return 2 * i + 1;
}


PRIVATE INLINE void
swap_entries(struct vrna_heap_s *h,
             size_t             a,
             size_t             b)
{
  void *v;

  v             = h->entries[b];
  h->entries[b] = h->entries[a];
  h->entries[a] = v;

  /* notify external storage about updated positions */
  if (h->set_entry_pos) {
    h->set_entry_pos(v, a, h->data);
    h->set_entry_pos(h->entries[b], b, h->data);
  }
}


PRIVATE INLINE int
heapify_up(struct vrna_heap_s *h,
           size_t             i)
{
  int ret = 0;

  while (i > 1) {
    size_t  parent  = parent_pos(i);
    void    *v      = h->entries[parent];

    if (h->cmp(v, h->entries[i], h->data) < 0)
      break;

    swap_entries(h, parent, i);

    i   = parent;
    ret = 1;
  }

  return ret;
}


PRIVATE INLINE void
heapify_down(struct vrna_heap_s *h,
             size_t             pos)
{
  void    *child_v, *v;
  size_t  last_pos, child_pos, child_pos2;

  last_pos = h->num_entries;

  /* nothing to do if already last element */
  if (pos == last_pos)
    return;

  v           = h->entries[pos];
  child_pos   = child0_pos(pos);
  child_pos2  = child1_pos(pos);
  child_v     = NULL;

  /* compare to 1st child */
  if (child_pos <= last_pos) {
    child_v = h->entries[child_pos];
    if (h->cmp(v, child_v, h->data) < 0) {
      child_pos = 0;
      child_v   = v;
    }
  } else {
    child_pos = 0;
    child_v   = v;
  }

  /* compare to 2nd child */
  if (child_pos2 <= last_pos) {
    v = h->entries[child_pos2];
    if (h->cmp(v, child_v, h->data) < 0)
      child_pos = child_pos2;
  }

  if (child_pos) {
    /* swap current node with child */
    swap_entries(h, pos, child_pos);

    heapify_down(h, child_pos);
  }
}
