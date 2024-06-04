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

#include "ViennaRNA/color_output.inc"

#include "ViennaRNA/utils/basic.h"
#include "ViennaRNA/utils/log.h"

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

/*
 #################################
 # PRIVATE FUNCTION DECLARATIONS #
 #################################
 */

/*
 #################################
 # BEGIN OF FUNCTION DEFINITIONS #
 #################################
 */

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
#ifndef VRNA_WITHOUT_TTY_COLORS
  if (isatty(fileno(stderr))) {
    fprintf(stderr, ANSI_COLOR_RED_B "ERROR: " ANSI_COLOR_RESET ANSI_COLOR_BRIGHT);
    vfprintf(stderr, format, args);
    fprintf(stderr, ANSI_COLOR_RESET "\n");
  } else {
#endif
  fprintf(stderr, "ERROR: ");
  vfprintf(stderr, format, args);
  fprintf(stderr, "\n");
#ifndef VRNA_WITHOUT_TTY_COLORS
}


#endif

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
#ifndef VRNA_WITHOUT_TTY_COLORS
  if (isatty(fileno(stderr))) {
    fprintf(stderr, ANSI_COLOR_MAGENTA_B "WARNING: " ANSI_COLOR_RESET ANSI_COLOR_BRIGHT);
    vfprintf(stderr, format, args);
    fprintf(stderr, ANSI_COLOR_RESET "\n");
  } else {
#endif
  fprintf(stderr, "WARNING: ");
  vfprintf(stderr, format, args);
  fprintf(stderr, "\n");
#ifndef VRNA_WITHOUT_TTY_COLORS
}


#endif
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
vrna_message_vinfo(FILE       *fp,
                   const char *format,
                   va_list    args)
{
  if (!fp)
    fp = stdout;

#ifndef VRNA_WITHOUT_TTY_COLORS
  if (isatty(fileno(fp))) {
    fprintf(fp, ANSI_COLOR_BLUE_B);
    vfprintf(fp, format, args);
    fprintf(fp, ANSI_COLOR_RESET "\n");
  } else {
#endif
  vfprintf(fp, format, args);
  fprintf(fp, "\n");
#ifndef VRNA_WITHOUT_TTY_COLORS
}


#endif
}


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
