#ifndef VIENNA_RNA_PACKAGE_STRING_H
#define VIENNA_RNA_PACKAGE_STRING_H

#include <stddef.h>
#include <string.h>

/**
 *  @file     datastructures/array.h
 *  @ingroup  data_structures, strings
 *  @brief    A macro-based dynamic array implementation
 */

/**
 *  @addtogroup strings
 *  @{
 */

typedef char *vrna_string_t;

/**
 *  @brief  The header of an array
 */
typedef struct vrna_string_header_s {
  size_t  len;  /**< @brief The length of the string */
  size_t  size; /**< @brief The actual capacity of an array */
  size_t  shift_post;
  char    backup;
} vrna_string_header_t;


#define VRNA_STRING_HEADER(s) ((vrna_string_header_t *)s - 1)

vrna_string_t
vrna_string_make(char const *str);

void
vrna_string_free(vrna_string_t str);


size_t
vrna_string_length(vrna_string_t const str);


size_t
vrna_string_size(vrna_string_t const str);


vrna_string_t
vrna_string_append(vrna_string_t str,
                   vrna_string_t const other);

vrna_string_t
vrna_string_append_cstring(vrna_string_t  str,
                           char const     *other);

vrna_string_t
vrna_string_make_space_for(vrna_string_t str,
                           size_t        add_len);

size_t
vrna_string_available_space(vrna_string_t const str);

/**
 * @}
 */

#endif
