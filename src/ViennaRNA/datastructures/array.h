#ifndef VIENNA_RNA_PACKAGE_ARRAY_H
#define VIENNA_RNA_PACKAGE_ARRAY_H

#include <stddef.h>
#include "ViennaRNA/utils/log.h"

#if !defined(VRNA_NO_INLINE)
# ifdef SWIG
#   define VRNA_NO_INLINE 
# else
#   if defined(_MSC_VER)
#     define VRNA_NO_INLINE __declspec(noinline)
#   else
#     define VRNA_NO_INLINE __attribute__ ((noinline))
#   endif
# endif
#endif

/**
 *  @file     datastructures/array.h
 *  @ingroup  data_structures, array_utils
 *  @brief    A macro-based dynamic array implementation
 */

/**
 *  @addtogroup array_utils
 *  @{
 */

/**
 *  @brief  The header of an array
 */
typedef struct vrna_array_header_s {
  size_t  num;  /**< @brief The number of elements in an array */
  size_t  size; /**< @brief The actual capacity of an array */
} vrna_array_header_t;

/**
 *  @brief Define an array
 */
#define vrna_array(Type) Type *

/**
 *  @brief  Make an array @p Name of type @p Type
 */
#define vrna_array_make(Type, Name) Type * Name; vrna_array_init(Name)


#ifndef VRNA_ARRAY_GROW_FORMULA
/**
 *  @brief The default growth formula for array
 */
#define VRNA_ARRAY_GROW_FORMULA(n)                      (1.4 * (n) + 8)
#endif

/**
 *  @brief  Retrieve a pointer to the header of an array @p input
 */
#define VRNA_ARRAY_HEADER(input)                        ((vrna_array_header_t *)(input) - 1)
/**
 *  @brief  Get the number of elements of an array @p input
 */
#define vrna_array_size(input)                          (VRNA_ARRAY_HEADER(input)->num)
/**
 *  @brief  Get the size of an array @p input, i.e. its actual capacity
 */
#define vrna_array_capacity(input)                      (VRNA_ARRAY_HEADER(input)->size)

/**
 *  @brief  Explicitely set the capacity of an array @p a
 */
#define vrna_array_set_capacity(a, capacity) do { \
  if (a) { \
    void **a_ptr = (void **)&(a); \
    *a_ptr = vrna__array_set_capacity((a), (capacity), sizeof(*(a))); \
  } \
} while (0)


/**
 *  @brief Explicitely set the capacity of an array
 *
 *  @note Do not use this function. Rather resort to the #vrna_array_set_capacity macro
 */
VRNA_NO_INLINE void *
vrna__array_set_capacity(void    *array,
                         size_t  capacity,
                         size_t  element_size);

/**
 *  @brief Initialize an array @p a with a particular pre-allocated size @p init_size
 *
 */
#define vrna_array_init_size(a, init_size) do { \
  void **a_ptr = (void **)&(a); \
  size_t size = sizeof(*(a)) * (size_t)(init_size) + sizeof(vrna_array_header_t); \
  vrna_array_header_t *h = (void *)vrna_alloc(size); \
  h->num           = 0; \
  h->size          = (size_t)init_size; \
  *a_ptr           = (void *)(h + 1); \
} while (0)

/**
 *  @brief Initialize an array @p a
 */
#define vrna_array_init(a)  vrna_array_init_size(a, VRNA_ARRAY_GROW_FORMULA(0));


/**
 *  @brief  Release memory of an array @p a
 */
#define vrna_array_free(a) do { \
  if (a) { \
    vrna_array_header_t *h = VRNA_ARRAY_HEADER(a); \
    free(h); \
  } \
} while (0)


/**
 *  @brief  Safely append an item to an array @p a
 */
#define vrna_array_append(a, item) do { \
  if (vrna_array_capacity(a) < vrna_array_size(a) + 1) \
    vrna_array_grow(a, vrna_array_size(a) + 1); \
  (a)[vrna_array_size(a)++] = (item); \
} while (0)


/**
 *  @brief  Grow an array @p a to provide a minimum capacity @p min_capacity
 */
#define vrna_array_grow(a, min_capacity) do { \
  size_t new_capacity = VRNA_ARRAY_GROW_FORMULA(vrna_array_capacity(a)); \
  if (new_capacity < (min_capacity)) \
    new_capacity = (min_capacity); \
  vrna_array_set_capacity(a, new_capacity); \
} while (0)

/**
 * @}
 */


#endif
