#ifndef VIENNA_RNA_PACKAGE_ARRAY_H
#define VIENNA_RNA_PACKAGE_ARRAY_H

#include <stddef.h>

typedef struct vrna_array_header_s {
  size_t  num;
  size_t  size;
  size_t  element_size;
} vrna_array_header_t;


#define VRNA_ARRAY_HEADER(input)                        ((vrna_array_header_t *)(input) - 1)
#define VRNA_ARRAY_SIZE(input)                          (VRNA_ARRAY_HEADER(input)->num)
#define VRNA_ARRAY_RESIZE_REQUIRED(input, num_elements) (VRNA_ARRAY_HEADER(input)->num + num_elements >= VRNA_ARRAY_HEADER(input)->size)


void *
vrna_array_init(size_t element_size,
                size_t init_size);


void
vrna_array_free(void *array);


void *
vrna_array_grow(void    *array,
                size_t  grow_size);

#endif
