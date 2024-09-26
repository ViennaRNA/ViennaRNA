/*
 *    ViennaRNA/utils/basic.c
 *
 *               c  Ivo L Hofacker and Walter Fontana
 *                        Vienna RNA package
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <time.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <stdint.h>
#include <stdarg.h>
#include <errno.h>

/* isatty() is not available for Windows systems */
#ifndef _WIN32
#include <unistd.h>
#endif

#if VRNA_WITH_PTHREADS
# include <pthread.h>
#endif

#include "ViennaRNA/intern/color_output.h"

#include "ViennaRNA/datastructures/array.h"
#include "ViennaRNA/utils/basic.h"
#include "ViennaRNA/utils/log.h"

#define EXIT_ON_ERROR

/*
 #################################
 # GLOBAL VARIABLES              #
 #################################
 */

/*
 #################################
 # PRIVATE VARIABLES             #
 #################################
 */
typedef struct {
  vrna_log_cb_f       cb;
  void                *cb_data;
  vrna_logdata_free_f data_release;
  vrna_log_levels_e   level;
} logger_callback;

PRIVATE struct {
  FILE              *default_file;
  vrna_log_levels_e default_level;
  unsigned int      options;
  vrna_log_lock_f   lock;
  void              *lock_data;
  vrna_array(logger_callback) callbacks;
#if VRNA_WITH_PTHREADS
  pthread_mutex_t   mtx;            /* semaphore to prevent concurrent access */
#endif
} logger = {
  .default_file   = NULL,
  .default_level  = VRNA_LOG_LEVEL_DEFAULT,
  .options        = VRNA_LOG_OPTION_DEFAULT,
  .lock           = NULL,
  .lock_data      = NULL,
  .callbacks      = NULL,
#if VRNA_WITH_PTHREADS
  .mtx            = PTHREAD_MUTEX_INITIALIZER
#endif
};

/*
 #################################
 # PRIVATE FUNCTION DECLARATIONS #
 #################################
 */
PRIVATE void
log_v(vrna_log_event_t *event);


PRIVATE void
log_default(vrna_log_event_t *event);


PRIVATE void
lock(void);


PRIVATE void
unlock(void);


#ifndef VRNA_WITHOUT_TTY_COLORS
PRIVATE const char *
get_log_level_color(vrna_log_levels_e level);
#endif


PRIVATE const char *
get_log_level_string(vrna_log_levels_e level);


/*
 #################################
 # BEGIN OF FUNCTION DEFINITIONS #
 #################################
 */
PUBLIC void
vrna_log(vrna_log_levels_e  level,
         const char         *file_name,
         int                line_number,
         const char         *format_string,
         ...)
{
  vrna_log_event_t event = {
    .format_string  = format_string,
    .level          = level,
    .line_number    = line_number,
    .file_name      = file_name
  };

  va_start(event.params, format_string);
  log_v(&event);
  va_end(event.params);
}


PUBLIC vrna_log_levels_e
vrna_log_level(void)
{
  return logger.default_level;
}


PUBLIC int
vrna_log_level_set(vrna_log_levels_e level)
{
  switch (level) {
    case VRNA_LOG_LEVEL_DEBUG:
    /* fall through */
    case VRNA_LOG_LEVEL_INFO:
    /* fall through */
    case VRNA_LOG_LEVEL_WARNING:
    /* fall through */
    case VRNA_LOG_LEVEL_ERROR:
    /* fall through */
    case VRNA_LOG_LEVEL_CRITICAL:
    /* fall through */
    case VRNA_LOG_LEVEL_SILENT:
      logger.default_level = level;
      break;
    default:
      vrna_log_warning("unkown log level specified! Not doing anything");
      level = VRNA_LOG_LEVEL_UNKNOWN;
      break;
  }

  return level;
}


PUBLIC FILE *
vrna_log_fp(void)
{
  if (!logger.default_file)
    logger.default_file = stderr;

  return logger.default_file;
}


PUBLIC void
vrna_log_fp_set(FILE *fp)
{
  if (fp)
    logger.default_file = fp;
}


PUBLIC unsigned int
vrna_log_options(void)
{
  return logger.options;
}


PUBLIC void
vrna_log_options_set(unsigned int options)
{
  logger.options = options;
}


PUBLIC size_t
vrna_log_cb_add(vrna_log_cb_f       cb,
                void                *data,
                vrna_logdata_free_f data_release,
                vrna_log_levels_e   level)
{
  /* initialize the logger, if not done already */
  if (logger.callbacks == NULL) {
    /* initialize callbacks if not done so far */
    vrna_array_init(logger.callbacks);
  }

  if (cb) {
    logger_callback logger_cb = {
      .cb           = cb,
      .cb_data      = data,
      .data_release = data_release,
      .level        = level
    };

    /* append callback to callback list */
    vrna_array_append(logger.callbacks, logger_cb);
  }

  return vrna_array_size(logger.callbacks);
}


PUBLIC size_t
vrna_log_cb_num(void)
{
  return vrna_array_size(logger.callbacks);
}


PUBLIC size_t
vrna_log_cb_remove(vrna_log_cb_f  cb,
                   void           *data)
{
  /* initialize the logger, if not done already */
  if (logger.callbacks == NULL) {
    /* initialize callbacks if not done so far */
    vrna_array_init(logger.callbacks);
  }

  /* search through logger callbacks to find
   * the logger we want to remove
   */
  if ((cb) &&
      (vrna_array_size(logger.callbacks) > 0)) {
    size_t i;

    for (i = 0; i < vrna_array_size(logger.callbacks); i++) {
      if ((logger.callbacks[i].cb == cb) &&
          (logger.callbacks[i].cb_data == data))
        break;
    }

    /* did we find the callback ? */
    if (i < vrna_array_size(logger.callbacks)) {
      /* release data if necessary */
      if (logger.callbacks[i].data_release)
        logger.callbacks[i].data_release(logger.callbacks[i].cb_data);

      /* move all callbacks after the one we found */
      if (i < vrna_array_size(logger.callbacks) - 1) {
        (void)memmove(logger.callbacks + i,
                      logger.callbacks + i + 1,
                      (vrna_array_size(logger.callbacks) - i - 1) * sizeof(logger_callback));
      }

      /* reduce number of callbacks */
      VRNA_ARRAY_HEADER(logger.callbacks)->num--;

      return 1;
    }
  }

  return 0; /* unsuccessful */
}


PUBLIC void
vrna_log_lock_set(vrna_log_lock_f cb,
                  void            *data)
{
  logger.lock       = cb;
  logger.lock_data  = data;
}


PUBLIC void
vrna_log_reset(void)
{
  if (logger.callbacks) {
    for (size_t i = 0; i < vrna_array_size(logger.callbacks); i++)
      if (logger.callbacks[i].data_release)
        logger.callbacks[i].data_release(logger.callbacks[i].cb_data);

    vrna_array_free(logger.callbacks);
  }

  /* initialize everything to default settings */
  logger.default_file   = stderr;
  logger.lock           = NULL;
  logger.lock_data      = NULL;
  logger.default_level  = VRNA_LOG_LEVEL_DEFAULT;
  logger.options        = VRNA_LOG_OPTION_DEFAULT;

  vrna_array_init(logger.callbacks);
#if VRNA_WITH_PTHREADS
  (void)pthread_mutex_init(&(logger.mtx), NULL);
#endif
}


PUBLIC void
vrna_message_error(const char *format,
                   ...)
{
  va_list args;

  va_start(args, format);
  vrna_message_verror(format, args);
  va_end(args);
}


PUBLIC void
vrna_message_verror(const char  *format,
                    va_list     args)
{
  vrna_log_event_t event = {
    .format_string  = format,
    .level          = VRNA_LOG_LEVEL_ERROR,
    .line_number    = __LINE__,
    .file_name      = __FILE__
  };

  va_copy(event.params, args);
  log_v(&event);
  va_end(event.params);

#ifdef EXIT_ON_ERROR
  exit(EXIT_FAILURE);
#endif
}


PUBLIC void
vrna_message_warning(const char *format,
                     ...)
{
  va_list args;

  va_start(args, format);
  vrna_message_vwarning(format, args);
  va_end(args);
}


PUBLIC void
vrna_message_vwarning(const char  *format,
                      va_list     args)
{
  vrna_log_event_t event = {
    .format_string  = format,
    .level          = VRNA_LOG_LEVEL_WARNING,
    .line_number    = __LINE__,
    .file_name      = __FILE__
  };

  va_copy(event.params, args);
  log_v(&event);
  va_end(event.params);
}


PUBLIC void
vrna_message_info(FILE        *fp,
                  const char  *format,
                  ...)
{
  va_list args;

  va_start(args, format);
  vrna_message_vinfo(fp, format, args);
  va_end(args);
}


PUBLIC void
vrna_message_vinfo(FILE       *fp_p,
                   const char *format,
                   va_list    args)
{
  FILE *fp_bak;

  /* remember current file pointer */
  fp_bak = vrna_log_fp();

  if (!fp_p)
    fp_p = stdout;

  /* set requested file pointer */
  vrna_log_fp_set(fp_p);

  vrna_log_event_t event = {
    .format_string  = format,
    .level          = VRNA_LOG_LEVEL_INFO,
    .line_number    = __LINE__,
    .file_name      = __FILE__
  };

  va_copy(event.params, args);
  log_v(&event);
  va_end(event.params);

  /* restore previous file pointer */
  vrna_log_fp_set(fp_bak);
}


/*
 #################################
 # STATIC helper functions below #
 #################################
 */
PRIVATE void
log_default(vrna_log_event_t *event)
{
  if (!logger.default_file)
    logger.default_file = stderr;

  /* print time unless turned off explicitely */
  if (logger.options & VRNA_LOG_OPTION_TRACE_TIME) {
    char    timebuf[64];
    time_t  t = time(NULL);
    timebuf[strftime(timebuf, sizeof(timebuf), "%H:%M:%S", localtime(&t))] = '\0';
    fprintf(logger.default_file, "%s ", timebuf);
  }

  /* print log level */
#ifndef VRNA_WITHOUT_TTY_COLORS
  if (isatty(fileno(logger.default_file))) {
    fprintf(logger.default_file,
            "%s%-9s" ANSI_COLOR_RESET " ",
            get_log_level_color(event->level),
            get_log_level_string(event->level));
  } else {
#endif
  fprintf(logger.default_file,
          "%-9s ",
          get_log_level_string(event->level));
#ifndef VRNA_WITHOUT_TTY_COLORS
}


#endif

  /* print file name / line number trace unless turned off explicitely */
  if (logger.options & VRNA_LOG_OPTION_TRACE_CALL) {
#ifndef VRNA_WITHOUT_TTY_COLORS
    if (isatty(fileno(logger.default_file))) {
      fprintf(logger.default_file,
              "\x1b[90m%s:%d:" ANSI_COLOR_RESET " ",
              event->file_name,
              event->line_number);
    } else {
#endif
    fprintf(logger.default_file, "%s:%d: ",
            event->file_name,
            event->line_number);
#ifndef VRNA_WITHOUT_TTY_COLORS
  }

#endif
  }

  /* print actual message */
  vfprintf(logger.default_file, event->format_string, event->params);
  fprintf(logger.default_file, "\n");
  fflush(logger.default_file);
}


PRIVATE void
log_v(vrna_log_event_t *event)
{
  lock();

  /* initialize the logger, if not done already */
  if (logger.callbacks == NULL) {
    /* initialize callbacks if not done so far */
    vrna_array_init(logger.callbacks);
  }

  /* process log output for default implementation */
  if (!(logger.options & VRNA_LOG_OPTION_QUIET)) {
    /* print log if not in quiet mode */
    if (event->level >= logger.default_level)
      log_default(event);
  }

  /* process log for any user-defined output */
  for (size_t i = 0; i < vrna_array_size(logger.callbacks); i++) {
    logger_callback *cb = &(logger.callbacks[i]);

    if (event->level >= cb->level)
      cb->cb(event, cb->cb_data);
  }

  unlock();
}


PRIVATE void
lock(void)
{
  if (logger.lock)
#if VRNA_WITH_PTHREADS
  {
#endif
    logger.lock(1, logger.lock_data);

#if VRNA_WITH_PTHREADS
} else {
  pthread_mutex_lock(&(logger.mtx));
}
#endif
}


PRIVATE void
unlock(void)
{
  if (logger.lock)
#if VRNA_WITH_PTHREADS
  {
#endif
    logger.lock(0, logger.lock_data);

#if VRNA_WITH_PTHREADS
} else {
  pthread_mutex_unlock(&(logger.mtx));
}
#endif
}


PRIVATE const char *
get_log_level_string(vrna_log_levels_e level)
{
  switch (level) {
    case VRNA_LOG_LEVEL_DEBUG:
      return "[DEBUG]";
    case VRNA_LOG_LEVEL_INFO:
      return "[INFO]";
    case VRNA_LOG_LEVEL_WARNING:
      return "[WARNING]";
    case VRNA_LOG_LEVEL_ERROR:
      return "[ERROR]";
    case VRNA_LOG_LEVEL_CRITICAL:
      return "[FATAL]";
    default:
      return "[UNKNOWN]";
  }
}


#ifndef VRNA_WITHOUT_TTY_COLORS
PRIVATE const char *
get_log_level_color(vrna_log_levels_e level)
{
  switch (level) {
    case VRNA_LOG_LEVEL_DEBUG:
      return ANSI_COLOR_CYAN_B;
    case VRNA_LOG_LEVEL_INFO:
      return ANSI_COLOR_BLUE_B;
    case VRNA_LOG_LEVEL_WARNING:
      return ANSI_COLOR_YELLOW_B;
    case VRNA_LOG_LEVEL_ERROR:
      return ANSI_COLOR_RED_B;
    case VRNA_LOG_LEVEL_CRITICAL:
      return ANSI_COLOR_MAGENTA_B;
    default:
      return ANSI_COLOR_RESET;
  }
}
#endif


#ifndef VRNA_DISABLE_BACKWARD_COMPATIBILITY

/*
 * ###########################################
 * # deprecated functions below              #
 *###########################################
 */
PUBLIC void
warn_user(const char message[])
{
  vrna_message_warning(message);
}


PUBLIC void
nrerror(const char message[])
{
  vrna_message_error(message);
}


#endif
