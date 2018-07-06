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
 *  @brief Functions that provide dynamically buffered stream-like data structures
 */

#include <stdarg.h>

/* below is our own implementation of a dynamic char * stream */
typedef struct vrna_cstr_s *vrna_cstr_t;

vrna_cstr_t
vrna_cstr(size_t  size,
          FILE    *output);


void
vrna_cstr_free(vrna_cstr_t buf);


void
vrna_cstr_close(vrna_cstr_t buf);


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
                           int                i,
                           int                L,
                           int                l[3],
                           int                energy);


/**
 *  @}
 */

#endif
