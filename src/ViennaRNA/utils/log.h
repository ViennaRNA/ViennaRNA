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

#include <stdio.h>
#include <stdarg.h>
#include <time.h>

/**
 *  @addtogroup  utils_log
 *  @{
 */

typedef struct vrna_log_event_s {
  const char  *format_string;
  va_list     params;
  int         level;
  int         line_number;
  const char  *file_name;
} vrna_log_event_t;


typedef void (*vrna_log_lock_f)(int lock, void *lock_data);
typedef void (*vrna_log_cb_f)(vrna_log_event_t *event, void *log_data);

enum {
  VRNA_LOG_LEVEL_UNKNOWN  = -1,
  VRNA_LOG_LEVEL_DEBUG    = 10,
  VRNA_LOG_LEVEL_INFO     = 20,
  VRNA_LOG_LEVEL_WARNING  = 30,
  VRNA_LOG_LEVEL_ERROR    = 40,
  VRNA_LOG_LEVEL_CRITICAL = 50
};
  
#define VRNA_LOG_LEVEL_DEFAULT      VRNA_LOG_LEVEL_ERROR
#define VRNA_LOG_OPTION_QUIET       1U
#define VRNA_LOG_OPTION_TRACE_CALL  2U
#define VRNA_LOG_OPTION_TRACE_TIME  4U
#define VRNA_LOG_OPTION_DEFAULT     0U

#define vrna_log_debug(...) \
    do { \
        vrna_log(VRNA_LOG_LEVEL_DEBUG, __FILE__, __LINE__, __VA_ARGS__); \
    } while (0)

#define vrna_log_info(...) \
    do { \
        vrna_log(VRNA_LOG_LEVEL_INFO, __FILE__, __LINE__, __VA_ARGS__); \
    } while (0)

#define vrna_log_warning(...) \
    do { \
        vrna_log(VRNA_LOG_LEVEL_WARNING, __FILE__, __LINE__, __VA_ARGS__); \
    } while (0)

#define vrna_log_error(...) \
    do { \
        vrna_log(VRNA_LOG_LEVEL_ERROR, __FILE__, __LINE__, __VA_ARGS__); \
    } while (0)

#define vrna_log_critical(...) \
    do { \
        vrna_log(VRNA_LOG_LEVEL_CRITICAL, __FILE__, __LINE__, __VA_ARGS__); \
    } while (0)


void
vrna_log(int level,
         const char *file_name,
         int        line_number,
         const char *format_string,
         ...);


int vrna_log_level(void);


int vrna_log_level_set(int level);


unsigned int vrna_log_options(void);


void vrna_log_options_set(unsigned int options);

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
