#ifndef VIENNA_RNA_PACKAGE_HEAP_H
#define VIENNA_RNA_PACKAGE_HEAP_H

#ifdef VRNA_WARN_DEPRECATED
# if defined(DEPRECATED)
#   undef DEPRECATED
# endif
# if defined(__clang__)
#  define DEPRECATED(func, msg) func __attribute__ ((deprecated("", msg)))
# elif defined(__GNUC__)
#  define DEPRECATED(func, msg) func __attribute__ ((deprecated(msg)))
# else
#  define DEPRECATED(func, msg) func
# endif
#else
# define DEPRECATED(func, msg) func
#endif

/**
 *  @file     ViennaRNA/datastructures/heap.h
 *  @ingroup  heap_utils
 *  @brief    Implementation of an abstract heap data structure
 */

/**
 *  @addtogroup heap_utils
 *  @{
 */


/**
 *  @brief  An abstract heap data structure
 *
 *  @see  vrna_heap_init(), vrna_heap_free(), vrna_heap_insert(),
 *        vrna_heap_pop(), vrna_heap_top(), vrna_heap_remove(),
 *        vrna_heap_update()
 */
typedef struct vrna_heap_s *vrna_heap_t;


/**
 *  @brief  Heap compare function prototype
 *
 *  Use this prototype to design the compare function for the heap implementation.
 *  The arbitrary data pointer @p data may be used to get access to further information
 *  required to actually compare the two values @p a and @p b.
 *
 *  @note The heap implementation acts as a <em>min-heap</em>, therefore, the minimum
 *        element will be present at the heap's root. In case a <em>max-heap</em> is
 *        required, simply reverse the logic of this compare function.
 *
 *  @param  a     The first object to compare
 *  @param  b     The second object to compare
 *  @param  data  An arbitrary data pointer passed through from the heap implementation
 *  @return       A value less than zero if @p a < @p b, a value greater than zero if @p a > @p b, and 0 otherwise
 */
typedef int (*vrna_heap_cmp_f)(const void *a,
                                     const void *b,
                                     void       *data);

DEPRECATED(typedef int (vrna_callback_heap_cmp)(const void *a,
                                     const void *b,
                                     void       *data),
           "Use vrna_heap_cmp_f instead!");


/**
 *  @brief  Retrieve the position of a particular heap entry within the heap
 *
 *  @param  a     The object to look-up within the heap
 *  @param  data  An arbitrary data pointer passed through from the heap implementation
 *  @return       The position of the element @p a within the heap, or 0 if it is not in the heap
 */
typedef size_t (*vrna_heap_get_pos_f)(const void  *a,
                                            void        *data);

DEPRECATED(typedef size_t (vrna_callback_heap_get_pos)(const void  *a,
                                            void        *data),
           "Use vrna_heap_get_pos_f instead!");


/**
 *  @brief  Store the position of a particular heap entry within the heap
 *
 *  @param  a     The object whose position shall be stored
 *  @param  pos   The current position of @p a within the heap, or 0 if a was deleted
 *  @param  data  An arbitrary data pointer passed through from the heap implementation
 */
typedef void (*vrna_heap_set_pos_f)(const void  *a,
                                          size_t      pos,
                                          void        *data);

DEPRECATED(typedef void (vrna_callback_heap_set_pos)(const void  *a,
                                          size_t      pos,
                                          void        *data),
           "USe vrna_heap_set_pos_f instead!");


/**
 *  @brief  Initialize a heap data structure
 *
 *  This function initializes a heap data structure. The implementation is based on a
 *  <em>min-heap</em>, i.e. the minimal element is located at the root of the heap.
 *  However, by reversing the logic of the compare function, one can easily transform
 *  this into a <em>max-heap</em> implementation.
 *
 *  Beside the regular operations on a heap data structure, we implement removal and
 *  update of arbitrary elements within the heap. For that purpose, however, one requires
 *  a reverse-index lookup system that, (i) for a given element stores the current position
 *  in the heap, and (ii) allows for fast lookup of an elements current position within the
 *  heap. The corresponding getter- and setter- functions may be provided through the
 *  arguments @p get_entry_pos and @p set_entry_pos, respectively.
 *
 *  Sometimes, it is difficult to simply compare two data structures without any context.
 *  Therefore, the compare function is provided with a user-defined data pointer that
 *  can hold any context required.
 *
 *  @warning  If any of the arguments @p get_entry_pos or @p set_entry_pos is NULL, the
 *            operations vrna_heap_update() and vrna_heap_remove() won't work.
 *
 *  @see  vrna_heap_free(), vrna_heap_insert(), vrna_heap_pop(), vrna_heap_top(),
 *        vrna_heap_remove(), vrna_heap_update(), #vrna_heap_t, #vrna_heap_cmp_f,
 *        #vrna_heap_get_pos_f, #vrna_heap_set_pos_f
 *
 *  @param  n             The initial size of the heap, i.e. the number of elements to store
 *  @param  cmp           The address of a compare function that will be used to fullfill the partial order requirement
 *  @param  get_entry_pos The address of a function that retrieves the position of an element within the heap (or NULL)
 *  @param  set_entry_pos The address of a function that stores the position of an element within the heap (or NULL)
 *  @param  data          An arbitrary data pointer passed through to the compare function @p cmp, and the set/get functions @p get_entry_pos / @p set_entry_pos
 *  @return               An initialized heap data structure, or NULL on error
 */
vrna_heap_t
vrna_heap_init(size_t                     n,
               vrna_heap_cmp_f     cmp,
               vrna_heap_get_pos_f get_entry_pos,
               vrna_heap_set_pos_f set_entry_pos,
               void                       *data);


/**
 *  @brief  Free memory occupied by a heap data structure
 *
 *  @see    vrna_heap_init()
 *
 *  @param  h   The heap that should be free'd
 */
void
vrna_heap_free(vrna_heap_t h);


/**
 *  @brief  Get the size of a heap data structure, i.e. the number of stored elements
 *
 *  @param  h   The heap data structure
 *  @return     The number of elements currently stored in the heap, or 0 upon any error
 */
size_t
vrna_heap_size(struct vrna_heap_s *h);


/**
 *  @brief  Insert an element into the heap
 *
 *  @see    vrna_heap_init(), vrna_heap_pop(), vrna_heap_top(), vrna_heap_free(),
 *          vrna_heap_remove(), vrna_heap_update()
 *
 *  @param  h     The heap data structure
 *  @param  v     A pointer to the object that is about to be inserted into the heap
 */
void
vrna_heap_insert(vrna_heap_t  h,
                 void         *v);


/**
 *  @brief  Pop (remove and return) the object at the root of the heap
 *
 *  This function removes the root from the heap and returns it to the caller.
 *
 *  @see    vrna_heap_init(), vrna_heap_top(), vrna_heap_insert(), vrna_heap_free()
 *          vrna_heap_remove(), vrna_heap_update()
 *
 *  @param  h   The heap data structure
 *  @return     The object at the root of the heap, i.e. the minimal element (or NULL if (a) the heap is empty or (b) any error occurred)
 */
void *
vrna_heap_pop(vrna_heap_t h);


/**
 *  @brief  Get the object at the root of the heap
 *
 *  @see    vrna_heap_init(), vrna_heap_pop(), vrna_heap_insert(), vrna_heap_free()
 *          vrna_heap_remove(), vrna_heap_update()
 *
 *  @param  h   The heap data structure
 *  @return     The object at the root of the heap, i.e. the minimal element (or NULL if (a) the heap is empty or (b) any error occurred)
 */
const void *
vrna_heap_top(vrna_heap_t h);


/**
 *  @brief  Remove an arbitrary element within the heap
 *
 *  @see    vrna_heap_init(), #vrna_heap_get_pos_f, #vrna_heap_set_pos_f,
 *          vrna_heap_pop(), vrna_heap_free()
 *
 *  @warning  This function won't work if the heap was not properly initialized with
 *            callback functions for fast reverse-index mapping!
 *
 *  @param  h   The heap data structure
 *  @param  v   The object to remove from the heap
 *  @return     The object that was removed from the heap (or NULL if (a) it wasn't found or (b) any error occurred)
 */
void *
vrna_heap_remove(vrna_heap_t  h,
                 const void   *v);


/**
 *  @brief  Update an arbitrary element within the heap
 *
 *  @note   If the object that is to be updated is not currently stored in the heap,
 *          it will be inserted. In this case, the function returns NULL.
 *
 *  @warning  This function won't work if the heap was not properly initialized with
 *            callback functions for fast reverse-index mapping!
 *
 *  @see    vrna_heap_init(), #vrna_heap_get_pos_f, #vrna_heap_set_pos_f
 *          vrna_heap_pop(), vrna_heap_remove(), vrna_heap_free()
 *
 *  @param  h   The heap data structure
 *  @param  v   The object to update
 *  @return     The 'previous' object within the heap that now got replaced by @p v (or NULL if (a) it wasn't found or (b) any error occurred)
 */
void *
vrna_heap_update(vrna_heap_t  h,
                 void         *v);


/**
 * @}
 */

#endif
