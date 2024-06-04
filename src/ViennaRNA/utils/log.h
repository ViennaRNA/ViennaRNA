#ifndef VIENNA_RNA_PACKAGE_UTILS_LOG_H
#define VIENNA_RNA_PACKAGE_UTILS_LOG_H

#ifdef VRNA_WARN_DEPRECATED
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
 *  @file     ViennaRNA/utils/log.h
 *  @ingroup  utils_log
 *  @brief    Logging system of the @em ViennaRNA @em Package
 */

/**
 *  @addtogroup  utils_log
 *  @{
 */

#include <stdio.h>
#include <stdarg.h>

/**
 *  @addtogroup  message_utils
 *  @{
 */

/**
 *  @brief Print an error message and die
 *
 *  This function is a wrapper to @em fprintf(stderr, ...) that
 *  puts a capital <b>ERROR:</b> in front of the message and then exits
 *  the calling program.
 *
 *  @see vrna_message_verror(), vrna_message_warning(), vrna_message_info()
 *
 *  @param format The error message to be printed
 *  @param ...    Optional arguments for the formatted message string
 */
void
vrna_message_error(const char *format,
                   ...);


/**
 *  @brief Print an error message and die
 *
 *  This function is a wrapper to @em vfprintf(stderr, ...) that
 *  puts a capital <b>ERROR:</b> in front of the message and then exits
 *  the calling program.
 *
 *  @see vrna_message_error(), vrna_message_warning(), vrna_message_info()
 *
 *  @param format The error message to be printed
 *  @param args   The argument list for the formatted message string
 */
void
vrna_message_verror(const char  *format,
                    va_list     args);


/**
 *  @brief Print a warning message
 *
 *  This function is a wrapper to @em fprintf(stderr, ...) that
 *  puts a capital <b>WARNING:</b> in front of the message.
 *
 *  @see vrna_message_vwarning(), vrna_message_error(), vrna_message_info()
 *
 *  @param format The warning message to be printed
 *  @param ...    Optional arguments for the formatted message string
 */
void
vrna_message_warning(const char *format,
                     ...);


/**
 *  @brief Print a warning message
 *
 *  This function is a wrapper to @em fprintf(stderr, ...) that
 *  puts a capital <b>WARNING:</b> in front of the message.
 *
 *  @see vrna_message_vwarning(), vrna_message_error(), vrna_message_info()
 *
 *  @param format The warning message to be printed
 *  @param args   The argument list for the formatted message string
 */
void
vrna_message_vwarning(const char  *format,
                      va_list     args);


/**
 *  @brief Print an info message
 *
 *  This function is a wrapper to @em fprintf(...).
 *
 *  @see vrna_message_vinfo(), vrna_message_error(), vrna_message_warning()
 *
 *  @param fp     The file pointer where the message is printed to
 *  @param format The warning message to be printed
 *  @param ...    Optional arguments for the formatted message string
 */
void
vrna_message_info(FILE        *fp,
                  const char  *format,
                  ...);


/**
 *  @brief Print an info message
 *
 *  This function is a wrapper to @em fprintf(...).
 *
 *  @see vrna_message_vinfo(), vrna_message_error(), vrna_message_warning()
 *
 *  @param fp     The file pointer where the message is printed to
 *  @param format The info message to be printed
 *  @param args   The argument list for the formatted message string
 */
void
vrna_message_vinfo(FILE       *fp,
                   const char *format,
                   va_list    args);


/**
 *  @}
 */

#ifndef VRNA_DISABLE_BACKWARD_COMPATIBILITY

/**
 *  @brief Print a warning message
 *
 *  Print a warning message to @e stderr
 *
 *  @deprecated Use vrna_message_warning() instead! (since v2.3.0)
 */
DEPRECATED(void warn_user(const char message[]), "Use vrna_message_warning() instead");

/**
 *  @brief Die with an error message
 *
 *  @deprecated Use vrna_message_error() instead! (since v2.3.0)
 */
DEPRECATED(void nrerror(const char message[]), "Use vrna_message_error() instead()");

#endif

#endif
