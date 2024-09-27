/*
 *  Buffered Character Stream
 *
 *  (c) 2018, Ronny Lorenz, ViennaRNA Package
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef _WIN32
#ifdef __MINGW32__
#include <unistd.h>
#else
#include "ViennaRNA/intern/unistd_win.h"
#endif
#else
#include <unistd.h>
#endif

#include <stdint.h>

#include "ViennaRNA/utils/basic.h"
#include "ViennaRNA/datastructures/char_stream.h"
#include "ViennaRNA/intern/color_output.h"

#define CSTR_OVERHEAD 4096

struct vrna_cstr_s {
  char          *string;
  size_t        size;
  FILE          *output;
  unsigned char istty;
};


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
struct vrna_cstr_s *
vrna_cstr(size_t  size,
          FILE    *output)
{
  struct vrna_cstr_s *buf = NULL;

  if (size == 0)
    size = CSTR_OVERHEAD;

  buf         = (struct vrna_cstr_s *)vrna_alloc(sizeof(struct vrna_cstr_s));
  buf->string = (char *)vrna_alloc(sizeof(char) * size);
  buf->size   = size;
  buf->output = (output) ? output : stdout;
  buf->istty  = isatty(fileno(buf->output));

  if (buf->string == NULL) {
    free(buf);
    return NULL;
  }

  buf->string[0] = '\0'; /* just in case */

  return buf;
}


void
vrna_cstr_free(struct vrna_cstr_s *buf)
{
  if (buf) {
    vrna_cstr_fflush(buf);
    free(buf->string);
    free(buf);
  }
}


void
vrna_cstr_discard(struct vrna_cstr_s *buf)
{
  if (buf) {
    buf->size       = CSTR_OVERHEAD;
    buf->string     = (char *)vrna_realloc(buf->string, sizeof(char) * buf->size);
    buf->string[0]  = '\0';
  }
}


void
vrna_cstr_close(struct vrna_cstr_s *buf)
{
  if (buf) {
    vrna_cstr_fflush(buf);

    free(buf->string);

    if ((buf->output != stdout) && (buf->output != stderr))
      fclose(buf->output);

    free(buf);
  }
}


const char *
vrna_cstr_string(struct vrna_cstr_s *buf)
{
  if (buf)
    return (const char *)buf->string;

  return NULL;
}


void
vrna_cstr_fflush(struct vrna_cstr_s *buf)
{
  if (buf) {
    if ((buf->output) &&
        (buf->string) &&
        (buf->string[0] != '\0')) {
      fprintf(buf->output, "%s", buf->string);
      (void)fflush(buf->output);
    }

    buf->size       = CSTR_OVERHEAD;
    buf->string     = (char *)vrna_realloc(buf->string, sizeof(char) * buf->size);
    buf->string[0]  = '\0';
  }
}


int
vrna_cstr_printf(struct vrna_cstr_s *buf,
                 const char         *format,
                 ...)
{
  int     r;
  va_list argp;

  if ((!buf) || (!format))
    return -1;

  va_start(argp, format);
  r = vrna_cstr_vprintf(buf, format, argp);
  va_end(argp); /* Each va_start() or va_copy() needs a va_end() */

  return r;
}


int
vrna_cstr_vprintf(struct vrna_cstr_s  *buf,
                  const char          *format,
                  va_list             args)
{
  char    *ptr;
  int     r;
  size_t  size_avail, size_old, size_new, l1, l2;

  if ((!buf) && (!format))
    return -1;

  va_list copy;

  va_copy(copy, args);

  r           = -1;
  ptr         = buf->string;
  size_avail  = buf->size;
  size_old    = (ptr) ? strlen(ptr) : 0;

  /* retrieve the number of characters that the string requires */
#ifdef _WIN32
  /*
   * vsnprintf() in Windows is not ANSI compliant, although it's
   * "...included for compliance to the ANSI standard"
   * Thus, we use _vscprintf() that explicitly counts characters
   */
  size_new = _vscprintf(format, args);
#else
  size_new = vsnprintf(NULL, 0, format, args);
#endif

  /* determine longer and shorter part of new string for size_t overflow protection */
  if (size_old > size_new) {
    l1  = size_old;
    l2  = size_new;
  } else {
    l1  = size_new;
    l2  = size_old;
  }

  if ((size_new > 0) &&
      (l1 < SIZE_MAX) &&
      (SIZE_MAX > (l1 + l2))) {
    /* increase string memory if necessary */
    if ((size_old + size_new + 1) > size_avail) {
      size_avail = size_old + size_new + 1;
      if (size_avail + CSTR_OVERHEAD < SIZE_MAX)
        size_avail += CSTR_OVERHEAD;

      ptr = (char *)vrna_realloc(ptr, sizeof(char) * (size_avail));
    }

    if (ptr == NULL) {
      r = -1;
    } else if ((r = vsnprintf(ptr + size_old, size_new + 1, format, copy)) < 0) {
      free(ptr);
    } else {
      buf->string = ptr;
      buf->size   = size_avail;
      r           = (int)(size_old + size_new);
    }
  } else if (size_new == 0) {
    /* we do not treat empty format string as error */
    r = (int)size_old;
  }

  va_end(copy); /* Each va_start() or va_copy() needs a va_end() */

  return r;
}


PUBLIC void
vrna_cstr_message_info(struct vrna_cstr_s *buf,
                       const char         *format,
                       ...)
{
  va_list args;

  if ((!buf) || (!format))
    return;

  va_start(args, format);
  vrna_cstr_message_vinfo(buf, format, args);
  va_end(args); /* Each va_start() or va_copy() needs a va_end() */
}


PUBLIC void
vrna_cstr_message_vinfo(struct vrna_cstr_s  *buf,
                        const char          *format,
                        va_list             args)
{
  if ((!buf) || (!format))
    return;

#ifndef VRNA_WITHOUT_TTY_COLORS
  if (buf->istty) {
    vrna_cstr_printf(buf, ANSI_COLOR_BLUE_B);
    vrna_cstr_vprintf(buf, format, args);
    vrna_cstr_printf(buf, ANSI_COLOR_RESET "\n");
  } else {
#endif
  vrna_cstr_vprintf(buf, format, args);
  vrna_cstr_printf(buf, "\n");
#ifndef VRNA_WITHOUT_TTY_COLORS
}


#endif
}

PUBLIC void
vrna_cstr_message_warning(struct vrna_cstr_s  *buf,
                          const char          *format,
                          ...)
{
  va_list args;

  if ((!buf) || (!format))
    return;

  va_start(args, format);
  vrna_cstr_message_vwarning(buf, format, args);
  va_end(args); /* Each va_start() or va_copy() needs a va_end() */
}


PUBLIC void
vrna_cstr_message_vwarning(struct vrna_cstr_s *buf,
                           const char         *format,
                           va_list            args)
{
  if ((!buf) || (!format))
    return;

#ifndef VRNA_WITHOUT_TTY_COLORS
  if (buf->istty) {
    vrna_cstr_printf(buf, ANSI_COLOR_MAGENTA_B "WARNING: " ANSI_COLOR_RESET ANSI_COLOR_BRIGHT);
    vrna_cstr_vprintf(buf, format, args);
    vrna_cstr_printf(buf, ANSI_COLOR_RESET "\n");
  } else {
#endif
  vrna_cstr_printf(buf, "WARNING: ");
  vrna_cstr_vprintf(buf, format, args);
  vrna_cstr_printf(buf, "\n");
#ifndef VRNA_WITHOUT_TTY_COLORS
}


#endif
}

PUBLIC void
vrna_cstr_print_fasta_header(struct vrna_cstr_s *buf,
                             const char         *head)
{
  if (head) {
#ifndef VRNA_WITHOUT_TTY_COLORS
    if (buf->istty)
      vrna_cstr_printf(buf, ANSI_COLOR_YELLOW ">%s" ANSI_COLOR_RESET "\n", head);
    else
#endif
    vrna_cstr_printf(buf, ">%s\n", head);
  }
}


PUBLIC void
vrna_cstr_printf_structure(struct vrna_cstr_s *buf,
                           const char         *structure,
                           const char         *format,
                           ...)

{
  va_list args;

  if (!buf)
    return;

  va_start(args, format);
  vrna_cstr_vprintf_structure(buf, structure, format, args);
  va_end(args); /* Each va_start() or va_copy() needs a va_end() */
}


PUBLIC void
vrna_cstr_vprintf_structure(struct vrna_cstr_s  *buf,
                            const char          *structure,
                            const char          *format,
                            va_list             args)
{
  if (!buf)
    return;

  if (structure)
    vrna_cstr_printf(buf, structure);

  if ((format) && (*format != '\0')) {
#ifndef VRNA_WITHOUT_TTY_COLORS
    if (buf->istty) {
      vrna_cstr_printf(buf, ANSI_COLOR_GREEN);
      vrna_cstr_vprintf(buf, format, args);
      vrna_cstr_printf(buf, ANSI_COLOR_RESET);
    } else {
#endif
    vrna_cstr_vprintf(buf, format, args);
#ifndef VRNA_WITHOUT_TTY_COLORS
  }

#endif
  }

  if ((structure) || ((format) && (*format != '\0')))
    vrna_cstr_printf(buf, "\n");
}


PUBLIC void
vrna_cstr_printf_comment(struct vrna_cstr_s *buf,
                         const char         *format,
                         ...)

{
  va_list args;

  if (!buf)
    return;

  va_start(args, format);
  vrna_cstr_vprintf_comment(buf, format, args);
  va_end(args); /* Each va_start() or va_copy() needs a va_end() */
}


PUBLIC void
vrna_cstr_vprintf_comment(struct vrna_cstr_s  *buf,
                          const char          *format,
                          va_list             args)
{
  if (!buf)
    return;

  if ((format) && (*format != '\0')) {
#ifndef VRNA_WITHOUT_TTY_COLORS
    if (buf->istty) {
      vrna_cstr_printf(buf, ANSI_COLOR_CYAN);
      vrna_cstr_vprintf(buf, format, args);
      vrna_cstr_printf(buf, ANSI_COLOR_RESET);
    } else {
#endif
    vrna_cstr_vprintf(buf, format, args);
#ifndef VRNA_WITHOUT_TTY_COLORS
  }

#endif
  }

  if ((format) && (*format != '\0'))
    vrna_cstr_printf(buf, "\n");
}


PUBLIC void
vrna_cstr_printf_thead(struct vrna_cstr_s *buf,
                       const char         *format,
                       ...)

{
  va_list args;

  if (!buf)
    return;

  va_start(args, format);
  vrna_cstr_vprintf_thead(buf, format, args);
  va_end(args); /* Each va_start() or va_copy() needs a va_end() */
}


PUBLIC void
vrna_cstr_vprintf_thead(struct vrna_cstr_s  *buf,
                        const char          *format,
                        va_list             args)
{
  if (!buf)
    return;

  if ((format) && (*format != '\0')) {
#ifndef VRNA_WITHOUT_TTY_COLORS
    if (buf->istty) {
      vrna_cstr_printf(buf, ANSI_COLOR_UNDERLINE ANSI_COLOR_BRIGHT);
      vrna_cstr_vprintf(buf, format, args);
      vrna_cstr_printf(buf, ANSI_COLOR_RESET);
    } else {
#endif
    vrna_cstr_vprintf(buf, format, args);
#ifndef VRNA_WITHOUT_TTY_COLORS
  }

#endif
  }

  if ((format) && (*format != '\0'))
    vrna_cstr_printf(buf, "\n");
}


PUBLIC void
vrna_cstr_printf_tbody(struct vrna_cstr_s *buf,
                       const char         *format,
                       ...)

{
  va_list args;

  if (!buf)
    return;

  va_start(args, format);
  vrna_cstr_vprintf_tbody(buf, format, args);
  va_end(args); /* Each va_start() or va_copy() needs a va_end() */
}


PUBLIC void
vrna_cstr_vprintf_tbody(struct vrna_cstr_s  *buf,
                        const char          *format,
                        va_list             args)
{
  if (!buf)
    return;

  if ((format) && (*format != '\0')) {
#ifndef VRNA_WITHOUT_TTY_COLORS
    if (buf->istty) {
      vrna_cstr_printf(buf, ANSI_COLOR_GREEN);
      vrna_cstr_vprintf(buf, format, args);
      vrna_cstr_printf(buf, ANSI_COLOR_RESET);
    } else {
#endif
    vrna_cstr_vprintf(buf, format, args);
#ifndef VRNA_WITHOUT_TTY_COLORS
  }

#endif
  }

  if ((format) && (*format != '\0'))
    vrna_cstr_printf(buf, "\n");
}


PUBLIC void
vrna_cstr_print_eval_sd_corr(struct vrna_cstr_s *buf)
{
  if (!buf)
    return;

#ifndef VRNA_WITHOUT_TTY_COLORS
  if (buf->istty) {
    vrna_cstr_printf(buf,
                     ANSI_COLOR_BRIGHT "Correcting for presence of structured domains" ANSI_COLOR_RESET "\n");
  } else {
#endif
  vrna_cstr_printf(buf, "Correcting for presence of structured domains\n");
#ifndef VRNA_WITHOUT_TTY_COLORS
}


#endif
}

PUBLIC void
vrna_cstr_print_eval_ext_loop(struct vrna_cstr_s  *buf,
                              int                 energy)
{
  if (!buf)
    return;

#ifndef VRNA_WITHOUT_TTY_COLORS
  if (buf->istty) {
    vrna_cstr_printf(buf,
                     ANSI_COLOR_CYAN "External loop" ANSI_COLOR_RESET
                     "                           : "
                     ANSI_COLOR_GREEN "%5d" ANSI_COLOR_RESET "\n",
                     energy);
  } else {
#endif
  vrna_cstr_printf(buf,
                   "External loop"
                   "                           : "
                   "%5d\n",
                   energy);
#ifndef VRNA_WITHOUT_TTY_COLORS
}


#endif
}


PUBLIC void
vrna_cstr_print_eval_ext_loop_revert(struct vrna_cstr_s  *buf,
                                     int                 energy)
{
  if (!buf)
    return;

#ifndef VRNA_WITHOUT_TTY_COLORS
  if (buf->istty) {
    vrna_cstr_printf(buf,
                     ANSI_COLOR_MAGENTA "External loop" ANSI_COLOR_RESET
                     "                           : "
                     ANSI_COLOR_GREEN "%5d" ANSI_COLOR_RESET "\n",
                     -energy);
  } else {
#endif
  vrna_cstr_printf(buf,
                   "External loop"
                   "                           : "
                   "%5d\n",
                   -energy);
#ifndef VRNA_WITHOUT_TTY_COLORS
}


#endif
}


PUBLIC void
vrna_cstr_print_eval_hp_loop(struct vrna_cstr_s *buf,
                             int                i,
                             int                j,
                             char               si,
                             char               sj,
                             int                energy)
{
  if (!buf)
    return;

#ifndef VRNA_WITHOUT_TTY_COLORS
  if (buf->istty) {
    vrna_cstr_printf(buf,
                     ANSI_COLOR_CYAN "Hairpin  loop" ANSI_COLOR_RESET
                     " (%3d,%3d) "
                     ANSI_COLOR_BRIGHT "%c%c" ANSI_COLOR_RESET
                     "              : "
                     ANSI_COLOR_GREEN "%5d" ANSI_COLOR_RESET "\n",
                     i, j,
                     si, sj,
                     energy);
  } else {
#endif
  vrna_cstr_printf(buf,
                   "Hairpin  loop"
                   " (%3d,%3d) %c%c              : "
                   "%5d\n",
                   i, j,
                   si, sj,
                   energy);
#ifndef VRNA_WITHOUT_TTY_COLORS
}


#endif
}


PUBLIC void
vrna_cstr_print_eval_hp_loop_revert(struct vrna_cstr_s  *buf,
                                    int                 i,
                                    int                 j,
                                    char                si,
                                    char                sj,
                                    int                 energy)
{
  if (!buf)
    return;

#ifndef VRNA_WITHOUT_TTY_COLORS
  if (buf->istty) {
    vrna_cstr_printf(buf,
                     ANSI_COLOR_MAGENTA "Hairpin  loop" ANSI_COLOR_RESET
                     " (%3d,%3d) "
                     ANSI_COLOR_BRIGHT "%c%c" ANSI_COLOR_RESET
                     "              : "
                     ANSI_COLOR_RED "%5d" ANSI_COLOR_RESET "\n",
                     i, j,
                     si, sj,
                     -energy);
  } else {
#endif
  vrna_cstr_printf(buf,
                   "Hairpin  loop"
                   " (%3d,%3d) %c%c              : "
                   "%5d\n",
                   i, j,
                   si, sj,
                   -energy);
#ifndef VRNA_WITHOUT_TTY_COLORS
}


#endif
}


PUBLIC void
vrna_cstr_print_eval_int_loop(struct vrna_cstr_s  *buf,
                              int                 i,
                              int                 j,
                              char                si,
                              char                sj,
                              int                 k,
                              int                 l,
                              char                sk,
                              char                sl,
                              int                 energy)
{
  if (!buf)
    return;

#ifndef VRNA_WITHOUT_TTY_COLORS
  if (buf->istty) {
    vrna_cstr_printf(buf,
                     ANSI_COLOR_CYAN "Interior loop" ANSI_COLOR_RESET
                     " (%3d,%3d) "
                     ANSI_COLOR_BRIGHT "%c%c" ANSI_COLOR_RESET
                     "; (%3d,%3d) "
                     ANSI_COLOR_BRIGHT "%c%c" ANSI_COLOR_RESET
                     ": "
                     ANSI_COLOR_GREEN "%5d" ANSI_COLOR_RESET "\n",
                     i, j,
                     si, sj,
                     k, l,
                     sk, sl,
                     energy);
  } else {
#endif
  vrna_cstr_printf(buf,
                   "Interior loop"
                   " (%3d,%3d) "
                   "%c%c"
                   "; (%3d,%3d) "
                   "%c%c"
                   ": "
                   "%5d\n",
                   i, j,
                   si, sj,
                   k, l,
                   sk, sl,
                   energy);
#ifndef VRNA_WITHOUT_TTY_COLORS
}


#endif
}


PUBLIC void
vrna_cstr_print_eval_int_loop_revert(struct vrna_cstr_s *buf,
                                     int                i,
                                     int                j,
                                     char               si,
                                     char               sj,
                                     int                k,
                                     int                l,
                                     char               sk,
                                     char               sl,
                                     int                energy)
{
  if (!buf)
    return;

#ifndef VRNA_WITHOUT_TTY_COLORS
  if (buf->istty) {
    vrna_cstr_printf(buf,
                     ANSI_COLOR_MAGENTA "Interior loop" ANSI_COLOR_RESET
                     " (%3d,%3d) "
                     ANSI_COLOR_BRIGHT "%c%c" ANSI_COLOR_RESET
                     "; (%3d,%3d) "
                     ANSI_COLOR_BRIGHT "%c%c" ANSI_COLOR_RESET
                     ": "
                     ANSI_COLOR_RED "%5d" ANSI_COLOR_RESET "\n",
                     i, j,
                     si, sj,
                     k, l,
                     sk, sl,
                     -energy);
  } else {
#endif
  vrna_cstr_printf(buf,
                   "Interior loop"
                   " (%3d,%3d) "
                   "%c%c"
                   "; (%3d,%3d) "
                   "%c%c"
                   ": "
                   "%5d\n",
                   i, j,
                   si, sj,
                   k, l,
                   sk, sl,
                   -energy);
#ifndef VRNA_WITHOUT_TTY_COLORS
}


#endif
}


PUBLIC void
vrna_cstr_print_eval_mb_loop(struct vrna_cstr_s *buf,
                             int                i,
                             int                j,
                             char               si,
                             char               sj,
                             int                energy)
{
  if (!buf)
    return;

#ifndef VRNA_WITHOUT_TTY_COLORS
  if (buf->istty) {
    vrna_cstr_printf(buf,
                     ANSI_COLOR_CYAN "Multi    loop" ANSI_COLOR_RESET
                     " (%3d,%3d) "
                     ANSI_COLOR_BRIGHT "%c%c" ANSI_COLOR_RESET
                     "              : "
                     ANSI_COLOR_GREEN "%5d" ANSI_COLOR_RESET "\n",
                     i, j,
                     si, sj,
                     energy);
  } else {
#endif
  vrna_cstr_printf(buf,
                   "Multi    loop"
                   " (%3d,%3d) %c%c              : "
                   "%5d\n",
                   i, j,
                   si, sj,
                   energy);
#ifndef VRNA_WITHOUT_TTY_COLORS
}


#endif
}


PUBLIC void
vrna_cstr_print_eval_mb_loop_revert(struct vrna_cstr_s  *buf,
                                    int                 i,
                                    int                 j,
                                    char                si,
                                    char                sj,
                                    int                 energy)
{
  if (!buf)
    return;

#ifndef VRNA_WITHOUT_TTY_COLORS
  if (buf->istty) {
    vrna_cstr_printf(buf,
                     ANSI_COLOR_MAGENTA "Multi    loop" ANSI_COLOR_RESET
                     " (%3d,%3d) "
                     ANSI_COLOR_BRIGHT "%c%c" ANSI_COLOR_RESET
                     "              : "
                     ANSI_COLOR_RED "%5d" ANSI_COLOR_RESET "\n",
                     i, j,
                     si, sj,
                     -energy);
  } else {
#endif
  vrna_cstr_printf(buf,
                   "Multi    loop"
                   " (%3d,%3d) %c%c              : "
                   "%5d\n",
                   i, j,
                   si, sj,
                   -energy);
#ifndef VRNA_WITHOUT_TTY_COLORS
}


#endif
}


PUBLIC void
vrna_cstr_print_eval_gquad(struct vrna_cstr_s *buf,
                           unsigned int                i,
                           unsigned int                j,
                           unsigned int       L,
                           unsigned int       l[3],
                           int                energy)
{
  if (!buf)
    return;

#ifndef VRNA_WITHOUT_TTY_COLORS
  if (buf->istty) {
    vrna_cstr_printf(buf,
                     ANSI_COLOR_CYAN "G-Quadruplex " ANSI_COLOR_RESET
                     " (%3u,%3u) "
                     ANSI_COLOR_BRIGHT "L%u  " ANSI_COLOR_RESET
                     "(%2u,%2u,%2u)  : "
                     ANSI_COLOR_GREEN "%5d" ANSI_COLOR_RESET "\n",
                     i, j,
                     L, l[0], l[1], l[2],
                     energy);
  } else {
#endif
  vrna_cstr_printf(buf,
                   "G-Quadruplex "
                   " (%3u,%3u) "
                   "L%u  "
                   "(%2u,%2u,%2u)  : "
                   "%5d\n",
                   i, j,
                   L, l[0], l[1], l[2],
                   energy);
#ifndef VRNA_WITHOUT_TTY_COLORS
}


#endif
}
