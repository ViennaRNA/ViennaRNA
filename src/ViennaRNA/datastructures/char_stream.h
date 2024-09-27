#ifndef VIENNA_RNA_PACKAGE_CHAR_STREAM_H
#define VIENNA_RNA_PACKAGE_CHAR_STREAM_H

/**
 *  @file     ViennaRNA/datastructures/char_stream.h
 *  @ingroup  utils, buffer_utils
 *  @brief    Implementation of a dynamic, buffered character stream
 */

/**
 *  @addtogroup   buffer_utils
 *  @{
 */

#include <stdarg.h>
#include <stdio.h>

/* below is our own implementation of a dynamic char * stream */
typedef struct vrna_cstr_s *vrna_cstr_t;

/**
 *  @brief  Create a dynamic char * stream data structure
 *
 *  @see  vrna_cstr_free(), vrna_cstr_close(), vrna_cstr_fflush(), vrna_cstr_discard(), vrna_cstr_printf()
 *
 *  @param  size    The initial size of the buffer in characters
 *  @param  output  An optional output file stream handle that is used to write the collected data to (defaults to @em stdout if @em NULL)
 */
vrna_cstr_t
vrna_cstr(size_t  size,
          FILE    *output);


/**
 *  @brief  Discard the current content of the dynamic char * stream data structure
 *
 *  @see  vrna_cstr_free(), vrna_cstr_close(), vrna_cstr_fflush(), vrna_cstr_printf()
 *
 *  @param  buf   The dynamic char * stream data structure to free
 */
void
vrna_cstr_discard(struct vrna_cstr_s *buf);


/**
 *  @brief  Free the memory occupied by a dynamic char * stream data structure
 *
 *  This function first flushes any remaining character data within the stream
 *  and then free's the memory occupied by the data structure.
 *
 *  @see vrna_cstr_close(), vrna_cstr_fflush(), vrna_cstr()
 *
 *  @param  buf   The dynamic char * stream data structure to free
 */
void
vrna_cstr_free(vrna_cstr_t buf);


/**
 *  @brief  Free the memory occupied by a dynamic char * stream and close the output stream
 *
 *  This function first flushes any remaining character data within the stream
 *  then closes the attached output file stream (if any), and finally free's the
 *  memory occupied by the data structure.
 *
 *  @see vrna_cstr_free(), vrna_cstr_fflush(), vrna_cstr()
 *
 *  @param  buf   The dynamic char * stream data structure to free
 */
void
vrna_cstr_close(vrna_cstr_t buf);


/**
 *  @brief  Flush the dynamic char * output stream
 *
 *  This function flushes the collected char * stream, either by writing
 *  to the attached file handle, or simply by writing to @em stdout if
 *  no file handle has been attached upon construction using vrna_cstr().
 *
 *  @post The stream buffer is empty after execution of this function
 *
 *  @see  vrna_cstr(), vrna_cstr_close(), vrna_cstr_free()
 *
 *  @param  buf   The dynamic char * stream data structure to flush
 */
void
vrna_cstr_fflush(struct vrna_cstr_s *buf);


const char *
vrna_cstr_string(vrna_cstr_t buf);


int
vrna_cstr_vprintf(vrna_cstr_t buf,
                  const char  *format,
                  va_list     args);


int
vrna_cstr_printf(vrna_cstr_t  buf,
                 const char   *format,
                 ...);


void
vrna_cstr_message_info(vrna_cstr_t  buf,
                       const char   *format,
                       ...);


void
vrna_cstr_message_vinfo(vrna_cstr_t buf,
                        const char  *format,
                        va_list     args);


void
vrna_cstr_message_warning(struct vrna_cstr_s  *buf,
                          const char          *format,
                          ...);


void
vrna_cstr_message_vwarning(struct vrna_cstr_s *buf,
                           const char         *format,
                           va_list            args);


void
vrna_cstr_print_fasta_header(vrna_cstr_t  buf,
                             const char   *head);


void
vrna_cstr_printf_structure(struct vrna_cstr_s *buf,
                           const char         *structure,
                           const char         *format,
                           ...);


void
vrna_cstr_vprintf_structure(struct vrna_cstr_s  *buf,
                            const char          *structure,
                            const char          *format,
                            va_list             args);


void
vrna_cstr_printf_comment(struct vrna_cstr_s *buf,
                         const char         *format,
                         ...);


void
vrna_cstr_vprintf_comment(struct vrna_cstr_s  *buf,
                          const char          *format,
                          va_list             args);


void
vrna_cstr_printf_thead(struct vrna_cstr_s *buf,
                       const char         *format,
                       ...);


void
vrna_cstr_vprintf_thead(struct vrna_cstr_s  *buf,
                        const char          *format,
                        va_list             args);


void
vrna_cstr_printf_tbody(struct vrna_cstr_s *buf,
                       const char         *format,
                       ...);


void
vrna_cstr_vprintf_tbody(struct vrna_cstr_s  *buf,
                        const char          *format,
                        va_list             args);


void
vrna_cstr_print_eval_sd_corr(struct vrna_cstr_s *buf);


void
vrna_cstr_print_eval_ext_loop(struct vrna_cstr_s  *buf,
                              int                 energy);


void
vrna_cstr_print_eval_ext_loop_revert(struct vrna_cstr_s  *buf,
                                     int                 energy);


void
vrna_cstr_print_eval_hp_loop(struct vrna_cstr_s *buf,
                             int                i,
                             int                j,
                             char               si,
                             char               sj,
                             int                energy);


void
vrna_cstr_print_eval_hp_loop_revert(struct vrna_cstr_s  *buf,
                                    int                 i,
                                    int                 j,
                                    char                si,
                                    char                sj,
                                    int                 energy);


void
vrna_cstr_print_eval_int_loop(struct vrna_cstr_s  *buf,
                              int                 i,
                              int                 j,
                              char                si,
                              char                sj,
                              int                 k,
                              int                 l,
                              char                sk,
                              char                sl,
                              int                 energy);


void
vrna_cstr_print_eval_int_loop_revert(struct vrna_cstr_s *buf,
                                     int                i,
                                     int                j,
                                     char               si,
                                     char               sj,
                                     int                k,
                                     int                l,
                                     char               sk,
                                     char               sl,
                                     int                energy);


void
vrna_cstr_print_eval_mb_loop(struct vrna_cstr_s *buf,
                             int                i,
                             int                j,
                             char               si,
                             char               sj,
                             int                energy);


void
vrna_cstr_print_eval_mb_loop_revert(struct vrna_cstr_s  *buf,
                                    int                 i,
                                    int                 j,
                                    char                si,
                                    char                sj,
                                    int                 energy);


void
vrna_cstr_print_eval_gquad(struct vrna_cstr_s *buf,
                           unsigned int       i,
                           unsigned int       j,
                           unsigned int       L,
                           unsigned int       l[3],
                           int                energy);


/**
 *  @}
 */

#endif
