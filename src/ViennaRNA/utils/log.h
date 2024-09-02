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


/**
 *  @brief  The log levels
 */
typedef enum {
  VRNA_LOG_LEVEL_UNKNOWN  = -1,   /**< Unknown log level */
  VRNA_LOG_LEVEL_DEBUG    = 10,   /**< Debug log level */
  VRNA_LOG_LEVEL_INFO     = 20,   /**< Info log level */
  VRNA_LOG_LEVEL_WARNING  = 30,   /**< Warning log level */
  VRNA_LOG_LEVEL_ERROR    = 40,   /**< Error log level */
  VRNA_LOG_LEVEL_CRITICAL = 50,   /**< Critical log level */
  VRNA_LOG_LEVEL_SILENT   = 999   /**< Silent log level */
} vrna_log_levels_e;


/**
 *  @brief  A log event
 */
typedef struct vrna_log_event_s {
  const char        *format_string; /**< The printf-like format string containing the log information */
  va_list           params;         /**< The parameters for the printf-like format string */
  vrna_log_levels_e level;          /**< The log level */
  int               line_number;    /**< The source code line number that issued the log */
  const char        *file_name;     /**< The source code file that issued the log */
} vrna_log_event_t;


/**
 *  @brief  The lock function prototype that may be passed to the logging system
 *
 *  @see  vrna_log_lock_set()
 *
 *  @param  lock      A parameter indicating whether to lock (lock != 0) or unlock (lock == 0)
 *  @param  lock_data An arbitrary user-defined data pointer for the user-defined locking system
 */
typedef void (*vrna_log_lock_f)(int   lock,
                                void  *lock_data);


/**
 *  @brief  The log callback function prototype
 *
 *  @see    vrna_log_cb_add(), vrna_log_cb_num(), vrna_log_cb_remove()
 *
 *  @param  event     The log event
 *  @param  log_data  An arbitrary user-defined data pointer for the user-define log message receiver
 */
typedef void (*vrna_log_cb_f)(vrna_log_event_t  *event,
                              void              *log_data);

typedef void (*vrna_logdata_free_f)(void *data);

/**
 *  @brief  Default log level
 *
 *  @see  vrna_log_level_set(), vrna_log_reset(),
 *        vrna_log_level(), #vrna_log_levels_e
 */
#define VRNA_LOG_LEVEL_DEFAULT      VRNA_LOG_LEVEL_ERROR


/**
 *  @brief  Log option to turn off internal logging
 *
 *  When this option is set via vrna_log_options_set()
 *  the internal logging system will be deactivated and
 *  only user-defined callbacks will be seeing any logs.
 *
 *  @see  vrna_log_options_set(), vrna_log_options(), vrna_log_reset(),
 *        #VRNA_LOG_OPTION_TRACE_CALL, #VRNA_LOG_OPTION_TRACE_TIME,
 *        #VRNA_LOG_OPTION_DEFAULT
 */
#define VRNA_LOG_OPTION_QUIET       1U


/**
 *  @brief  Log option to turn on call tracing
 *
 *  When this option is set via vrna_log_options_set()
 *  the internal logging system will include a call
 *  trace to the log output, i.e. the source code file
 *  and line numbers will be included in the log message.
 *
 *  @see  vrna_log_options_set(), vrna_log_options(), vrna_log_reset(),
 *        #VRNA_LOG_OPTION_QUIET, #VRNA_LOG_OPTION_TRACE_TIME,
 *        #VRNA_LOG_OPTION_DEFAULT
 */
#define VRNA_LOG_OPTION_TRACE_CALL  2U


/**
 *  @brief  Log option to turn on time stamp
 *
 *  When this option is set via vrna_log_options_set()
 *  the internal logging system will include a time stamp
 *  to the log output, i.e. the time when the log message
 *  was issued will be included in the log message.
 *
 *  @see  vrna_log_options_set(), vrna_log_options(), vrna_log_reset(),
 *        #VRNA_LOG_OPTION_QUIET, #VRNA_LOG_OPTION_TRACE_CALL,
 *        #VRNA_LOG_OPTION_DEFAULT
 */
#define VRNA_LOG_OPTION_TRACE_TIME  4U


/**
 *  @brief  Log option representing the default options
 *
 *  When this option is set via vrna_log_options_set()
 *  the default options will be set.
 *
 *  @see  vrna_log_options_set(), vrna_log_options(), vrna_log_reset(),
 *        #VRNA_LOG_OPTION_QUIET, #VRNA_LOG_OPTION_TRACE_CALL,
 *        #VRNA_LOG_OPTION_TRACE_TIME
 */
#define VRNA_LOG_OPTION_DEFAULT     0U


/**
 *  @brief  Issue a debug log message
 *
 *  This macro expects a printf-like format string followed by a variable list of
 *  arguments for the format string and passes this content to the log system.
 *
 *  @see  #vrna_log_info, #vrna_log_warning, #vrna_log_error, #vrna_log_critical,
 *        vrna_log(), vrna_log_level_set(), vrna_log_options_set(), vrna_log_fp_set()
 */
#ifndef VRNA_LOG_NO_DEBUG_RNALIB
# define vrna_log_debug(...) \
        do { \
          vrna_log(VRNA_LOG_LEVEL_DEBUG, __FILE__, __LINE__, __VA_ARGS__); \
        } while (0)
#else
# define vrna_log_debug(...)
#endif
/**
 *  @brief  Issue an info log message
 *
 *  This macro expects a printf-like format string followed by a variable list of
 *  arguments for the format string and passes this content to the log system.
 *
 *  @see  #vrna_log_debug, #vrna_log_warning, #vrna_log_error, #vrna_log_critical,
 *        vrna_log(), vrna_log_level_set(), vrna_log_options_set(), vrna_log_fp_set()
 */
#define vrna_log_info(...) \
        do { \
          vrna_log(VRNA_LOG_LEVEL_INFO, __FILE__, __LINE__, __VA_ARGS__); \
        } while (0)

/**
 *  @brief  Issue a warning log message
 *
 *  This macro expects a printf-like format string followed by a variable list of
 *  arguments for the format string and passes this content to the log system.
 *
 *  @see  #vrna_log_debug, #vrna_log_info, #vrna_log_error, #vrna_log_critical,
 *        vrna_log(), vrna_log_level_set(), vrna_log_options_set(), vrna_log_fp_set()
 */
#define vrna_log_warning(...) \
        do { \
          vrna_log(VRNA_LOG_LEVEL_WARNING, __FILE__, __LINE__, __VA_ARGS__); \
        } while (0)

/**
 *  @brief  Issue an error log message
 *
 *  This macro expects a printf-like format string followed by a variable list of
 *  arguments for the format string and passes this content to the log system.
 *
 *  @see  #vrna_log_debug, #vrna_log_info, #vrna_log_warning, #vrna_log_critical,
 *        vrna_log(), vrna_log_level_set(), vrna_log_options_set(), vrna_log_fp_set()
 */
#define vrna_log_error(...) \
        do { \
          vrna_log(VRNA_LOG_LEVEL_ERROR, __FILE__, __LINE__, __VA_ARGS__); \
        } while (0)

/**
 *  @brief  Issue a critical log message
 *
 *  This macro expects a printf-like format string followed by a variable list of
 *  arguments for the format string and passes this content to the log system.
 *
 *  @see  #vrna_log_debug, #vrna_log_info, #vrna_log_warning, #vrna_log_error,
 *        vrna_log(), vrna_log_level_set(), vrna_log_options_set(), vrna_log_fp_set()
 */
#define vrna_log_critical(...) \
        do { \
          vrna_log(VRNA_LOG_LEVEL_CRITICAL, __FILE__, __LINE__, __VA_ARGS__); \
        } while (0)


/**
 *  @brief  Issue a log message
 *
 *  This is the low-level log message function. Usually, you don't want to call
 *  it directly but rather call one of the following high-level macros instead:
 *
 *  - #vrna_log_debug
 *  - #vrna_log_info
 *  - #vrna_log_warning
 *  - #vrna_log_error
 *  - #vrna_log_critical
 *
 *  @see  #vrna_log_debug, #vrna_log_info, #vrna_log_warning, #vrna_log_error,
 *        #vrna_log_critical, vrna_log_level_set(), vrna_log_options_set(),
 *        vrna_log_fp_set()
 *
 *  @param  level         The log level
 *  @param  file_name     The source code file name of the file that issued the log
 *  @param  line_number   The source code line number that issued the log
 *  @param  format_string The printf-like format string containing the log message
 *  @param  ...           The variable argument list for the printf-like @p format_string
 */
void
vrna_log(vrna_log_levels_e  level,
         const char         *file_name,
         int                line_number,
         const char         *format_string,
         ...);


/**
 *  @brief  Get the current default log level
 *
 *  @see vrna_log_level_set(), #vrna_log_levels_e
 *
 *  @returns  The current default log level
 */
vrna_log_levels_e
vrna_log_level(void);


/**
 *  @brief  Set the default log level
 *
 *  Set the log level for the default log output system.
 *  Any user-defined log callback mechanism will not be affected...
 *
 *  @see  vrna_log_level(), #vrna_log_levels_e,
 *        vrna_log_cb_add(), vrna_log_reset()
 *
 *  @param    level   The new log level for the default logging system
 *  @returns          The (updated) log level of the default logging system
 */
int
vrna_log_level_set(vrna_log_levels_e level);


/**
 *  @brief  Get the current log options of the default logging system
 *
 *  @see  vrna_log_options_set(),
 *        #VRNA_LOG_OPTION_QUIET, #VRNA_LOG_OPTION_TRACE_TIME
 *        #VRNA_LOG_OPTION_TRACE_CALL, #VRNA_LOG_OPTION_DEFAULT
 *
 *  @return The current options for the default logging system
 */
unsigned int
vrna_log_options(void);


/**
 *  @brief  Set the log options for the default logging system
 *
 *  @see  vrna_log_options(),
 *        #VRNA_LOG_OPTION_QUIET, #VRNA_LOG_OPTION_TRACE_TIME
 *        #VRNA_LOG_OPTION_TRACE_CALL, #VRNA_LOG_OPTION_DEFAULT
 *
 *  @param  options   The new options for the default logging system
 */
void
vrna_log_options_set(unsigned int options);


/**
 *  @brief  Get the output file pointer for the default logging system
 *
 *  @returns  The file pointer where the default logging system will print log messages to
 */
FILE *
vrna_log_fp(void);


/**
 *  @brief  Set the output file pointer for the default logging system
 *
 *  @param  fp  The file pointer where the default logging system should print log messages to
 */
void
vrna_log_fp_set(FILE *fp);


/**
 *  @brief  Add a user-defined log message callback
 *
 *  This function will add the user-defined callback @p cb to the logging system
 *  that will receive log messages from RNAlib. The callback will be
 *  called for each issued message that has a level of at least @p level.
 *  The pointer @p data will be passed-through to the callback and may
 *  store arbitrary data required for the callback.
 *
 *  @see    #vrna_log_cb_f, vrna_log_cb_num(), vrna_log_cb_remove(), vrna_log_reset()
 *
 *  @param  cb    The callback function
 *  @param  data  The data passed through to the callback function
 *  @param  data_release  A function that releases memory occupied by @p data (maybe NULL)
 *  @param  level The log level threshold for this callback
 *  @returns      The current number of log message callbacks stored in the logging system
 */
size_t
vrna_log_cb_add(vrna_log_cb_f       cb,
                void                *data,
                vrna_logdata_free_f data_release,
                vrna_log_levels_e   level);


/**
 *  @brief  Get the current number of log message callbacks
 *
 *  @returns      The current number of log message callbacks stored in the logging system
 */
size_t
vrna_log_cb_num(void);


/**
 *  @brief Remove a log message callback
 *
 *  This function removes the log message callback @p cb from the logging system.
 *  It does so by searching through the list of known log message callbacks and
 *  comparing function (and data) addresses.
 *
 *  @warning  The first callback stored in the logging system that matches @p cb
 *            will be removed! If @p data is supplied as well, the first callback
 *            that matches both, function and data address will be removed.
 *
 *  @see    #vrna_log_cb_f, vrna_log_cb_num(), vrna_log_cb_add(), vrna_log_reset()
 *
 *  @param  cb    The callback function to remove
 *  @param  data  The data that goes along with the callback
 *  @returns      0 on any error, e.g. if the callback was not found, non-zero if it was removed
 */
size_t
vrna_log_cb_remove(vrna_log_cb_f  cb,
                   void           *data);


/**
 *  @brief  Specify a lock function to be used for the logging system
 *
 *  To prevent undefined behavior in multi-threaded calls to the log
 *  system, each log message should be issued as an atomic block. For
 *  this to happen, a locking-/unlocking mechanism is required that
 *  ensures that log messages from other threads will be blocked until
 *  the current log message has been finalized. By default, we use
 *  pthreads mutex locking. Using this function, the locking mechanism
 *  can be changed to something else, or implemented after all if
 *  the ViennaRNA Package was compiled without pthreads support.
 *
 *  @see  #vrna_log_lock_f, vrna_log_reset()
 *
 *  @param  cb    The locking-/unlocking callback
 *  @param  data  An arbitrary data pointer passed through to the callback
 */
void
vrna_log_lock_set(vrna_log_lock_f cb,
                  void            *data);


/**
 *  @brief  Reset the logging system
 *
 *  This resets the logging system and restores default settings
 */
void
vrna_log_reset(void);


/**
 *  @}
 */


#endif
