#include <stddef.h>
#include <stdlib.h>

#include "ViennaRNA/utils/basic.h"
#include "ViennaRNA/datastructures/array.h"

PUBLIC void *
vrna_array_init(size_t element_size,
                size_t init_size)
{
  void    *array;
  size_t  size;

  size      = element_size * init_size + sizeof(struct vrna_array_header_s);
  array     = (void *)vrna_alloc(size);
  array     += sizeof(struct vrna_array_header_s);
  struct vrna_array_header_s *header = VRNA_ARRAY_HEADER(array);
  header->num           = 0;
  header->size          = init_size;
  header->element_size  = element_size;

  return array;
}


PUBLIC void
vrna_array_free(void *array) {
  if (array) {
    array -= sizeof(struct vrna_array_header_s);
    free(array);
  }
}


PUBLIC void *
vrna_array_grow(void    *array,
                size_t  grow_size) {
  size_t                      size;
  struct vrna_array_header_s  *header;

  if (array) {
    header  = VRNA_ARRAY_HEADER(array);
    size    = header->element_size * (header->size + grow_size) + sizeof(struct vrna_array_header_s);
    array  -= sizeof(struct vrna_array_header_s);
    array   = (void *)vrna_realloc(array, size);

    if (array) {
      array += sizeof(struct vrna_array_header_s);
      header = VRNA_ARRAY_HEADER(array);
      header->size += grow_size;
    }
  }

  return array;
}
