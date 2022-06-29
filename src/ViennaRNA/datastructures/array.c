#include <stddef.h>
#include <stdlib.h>
#include <string.h>

#include "ViennaRNA/utils/basic.h"
#include "ViennaRNA/datastructures/array.h"



PUBLIC VRNA_NO_INLINE void *
vrna__array_set_capacity(void    *array,
                         size_t  capacity,
                         size_t  element_size)
{
  vrna_array_header_t *h = VRNA_ARRAY_HEADER(array);

  if (capacity == h->size)
    return array;

  /* shrink array (and remove trailing elements */
  if (capacity < h->num) {
    if (h->size < capacity) {
      size_t new_capacity = VRNA_ARRAY_GROW_FORMULA(h->size);
      if (new_capacity < capacity)
        new_capacity = capacity;
      vrna__array_set_capacity(array, new_capacity, element_size);
    }
    h->num = capacity;
  }

  /* move memory to new location */
  size_t size = sizeof(vrna_array_header_t) + element_size * capacity;
  vrna_array_header_t *nh = (vrna_array_header_t *)vrna_alloc(size);
  memmove(nh, h, sizeof(vrna_array_header_t) + element_size * h->num);
  nh->num     = h->num;
  nh->size    = capacity;
  free(h);

  return nh + 1;
}
