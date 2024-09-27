#include <stdlib.h>

#include "ViennaRNA/utils/basic.h"
#include "ViennaRNA/datastructures/string.h"

PRIVATE vrna_string_t
vrna_string_make_length(void const *init_str,
                        size_t      len);

PRIVATE vrna_string_t
vrna_string_append_length(vrna_string_t str,
                          void const    *other,
                          size_t        other_len);


size_t
vrna_string_length(vrna_string_t const str) {
  return VRNA_STRING_HEADER(str)->len;
}

size_t
vrna_string_size(vrna_string_t const str) {
  return VRNA_STRING_HEADER(str)->size;
}


PRIVATE void
set_string_length(vrna_string_t str,
                  size_t        len)
{
  VRNA_STRING_HEADER(str)->len = len;
}

PRIVATE void
set_string_size(vrna_string_t str,
                size_t        cap) {
  VRNA_STRING_HEADER(str)->size = cap;
}

PRIVATE vrna_string_t
vrna_string_make_length(void const *init_str,
                        size_t      len)
{
  vrna_string_t         str;
  vrna_string_header_t  *header;
  size_t                header_size = sizeof(vrna_string_header_t);

  void *ptr = vrna_alloc(header_size + len + 1);

  if (ptr == NULL)
    return NULL;

  if (!init_str)
    memset(ptr, 0, header_size + len + 1);

  str = (char *)ptr + header_size;
  header = VRNA_STRING_HEADER(str);
  header->len   = len;
  header->size  = len;
  if (len && init_str)
    memcpy(str, init_str, len);
  str[len] = '\0';

  return str;
}


vrna_string_t
vrna_string_make(char const *str) {
  size_t len = str ? strlen(str) : 0;
  return vrna_string_make_length(str, len);
}


void
vrna_string_free(vrna_string_t str) {
  if (str == NULL)
    return;

  free((vrna_string_header_t *)str - 1);
}


size_t
vrna_string_available_space(vrna_string_t const str) {
  vrna_string_header_t *h = VRNA_STRING_HEADER(str);

  if (h->size > h->len)
    return h->size - h->len;

  return 0;
}

vrna_string_t
vrna_string_make_space_for(vrna_string_t str,
                           size_t        add_len)
{
  size_t len      = vrna_string_length(str);
  size_t  new_len = len + add_len;
  void *ptr, *new_ptr;
  size_t available, new_size;

  available = vrna_string_available_space(str);
  if (available >= add_len) /* Return if there is enough space left */
    return str;

  ptr = (char *)str - sizeof(vrna_string_header_t);
  new_size = sizeof(vrna_string_header_t) + new_len + 1;

  new_ptr = vrna_realloc(ptr, new_size);

  if (new_ptr == NULL)
    return NULL;

  str = (char *)new_ptr + sizeof(vrna_string_header_t);

  set_string_size(str, new_len);

  return str;
}

PRIVATE vrna_string_t
vrna_string_append_length(vrna_string_t str,
                          void const    *other,
                          size_t        other_len)
{
  size_t curr_len = vrna_string_length(str);

  str = vrna_string_make_space_for(str, other_len);
  if (str == NULL)
    return NULL;

  memcpy(str + curr_len, other, other_len);
  str[curr_len + other_len] = '\0';
  set_string_length(str, curr_len + other_len);

  return str;
}

vrna_string_t
vrna_string_append(vrna_string_t str,
                   vrna_string_t const other)
{
  return vrna_string_append_length(str,
                                   other,
                                   vrna_string_length(other));
}

vrna_string_t
vrna_string_append_cstring(vrna_string_t  str,
                           char const     *other)
{
  return vrna_string_append_length(str,
                                   other,
                                   strlen(other));
}
