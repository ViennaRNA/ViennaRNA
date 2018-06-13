/*
 *    ViennaRNA/utils/strings.c
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
#include <errno.h>
#include <time.h>
#include <string.h>
#include <sys/types.h>
#include <stdint.h>
#include <stdarg.h>

#include "ViennaRNA/utils/basic.h"
#include "ViennaRNA/utils/strings.h"

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

#ifndef HAVE_STRDUP
char *
strdup(const char *s)
{
  char *dup;

  dup = vrna_alloc(strlen(s) + 1);
  strcpy(dup, s);
  return dup;
}


#endif

PUBLIC char *
vrna_strdup_printf(const char *format,
                   ...)
{
  char    *result;
  va_list argp;

  va_start(argp, format);
  result = vrna_strdup_vprintf(format, argp);
  va_end(argp); /* Each va_start() or va_copy() needs a va_end() */

  return result;
}


PUBLIC char *
vrna_strdup_vprintf(const char  *format,
                    va_list     argp)
{
  char    *result;
  int     r;

  result = NULL;

#ifndef HAVE_VASPRINTF
  int     count;
  va_list copy;
  va_copy(copy, argp);

  r = -1;

  /* retrieve the number of characters that the string requires */
#ifdef _WIN32
  /*
   * vsnprintf() in Windows is not ANSI compliant, although it's
   * "...included for compliance to the ANSI standard"
   * Thus, we use _vscprintf() that explicitly counts characters
   */
  count = _vscprintf(format, argp);
#else
  count = vsnprintf(NULL, 0, format, argp);
#endif

  if ((count >= 0) && (count < INT_MAX)) {
    char *buf = (char *)vrna_alloc(sizeof(char) * (count + 1));
    if (buf == NULL)
      r = -1;
    else if ((r = vsnprintf(buf, count + 1, format, copy)) < 0)
      free(buf);
    else
      result = buf;
  }

  va_end(copy);  /* Each va_start() or va_copy() needs a va_end() */
#else
  /* the default is to use vasprintf() if available */
  r = vasprintf(&result, format, argp);
#endif

  /* check for any memory allocation error indicated by r == -1 */
  if (r == -1) {
    vrna_message_warning("vrna_strdup_printf: memory allocation failure!");
    result = NULL;
  }

  return result;
}


PUBLIC int
vrna_strcat_printf(char       **dest,
                   const char *format,
                   ...)
{
  int     r;
  va_list argp;

  va_start(argp, format);
  r = vrna_strcat_vprintf(dest, format, argp);
  va_end(argp); /* Each va_start() or va_copy() needs a va_end() */

  return r;
}


PUBLIC int
vrna_strcat_vprintf(char        **dest,
                    const char  *format,
                    va_list     args)
{
  char    *buf;
  int     r, l1, l2;
  size_t  old_count, new_count;

  if ((!dest) || (!format))
    return -1;

  va_list copy;
  va_copy(copy, args);

  r         = -1;
  buf       = *dest;
  old_count = (buf) ? strlen(buf) : 0;

  /* retrieve the number of characters that the string requires */
#ifdef _WIN32
  /*
   * vsnprintf() in Windows is not ANSI compliant, although it's
   * "...included for compliance to the ANSI standard"
   * Thus, we use _vscprintf() that explicitly counts characters
   */
  new_count = _vscprintf(format, args);
#else
  new_count = vsnprintf(NULL, 0, format, args);
#endif

  /* determine longer and shorter part of new string for INT overflow protection */
  if (old_count > new_count) {
    l1  = old_count;
    l2  = new_count;
  } else {
    l1  = new_count;
    l2  = old_count;
  }

  if ((new_count > 0) && (l1 < SIZE_MAX) && ((SIZE_MAX - l1) > l2)) {
    buf = (char *)vrna_realloc(buf, sizeof(char) * (old_count + new_count + 1));
    if (buf == NULL) {
      r = -1;
    } else if ((r = vsnprintf(buf + old_count, new_count + 1, format, copy)) < 0) {
      free(buf);
    } else {
      *dest = buf;
      r     = old_count + new_count;
    }
  } else if (new_count == 0) {
    /* we do not treat empty format string as error */
    r = (int)old_count;
  }

  va_end(copy);  /* Each va_start() or va_copy() needs a va_end() */

  /* check for any memory allocation error indicated by r == -1 */
  if (r == -1) {
    vrna_message_warning("vrna_strcat_printf: memory allocation failure!");
    *dest = NULL;
  }

  return r;
}


PUBLIC char *
vrna_random_string(int        l,
                   const char symbols[])
{
  char  *r;
  int   i, rn, base;

  base  = (int)strlen(symbols);
  r     = (char *)vrna_alloc(sizeof(char) * (l + 1));

  for (i = 0; i < l; i++) {
    rn    = (int)(vrna_urn() * base); /* [0, base-1] */
    r[i]  = symbols[rn];
  }
  r[l] = '\0';
  return r;
}


/*-----------------------------------------------------------------*/

PUBLIC int
vrna_hamming_distance(const char  *s1,
                      const char  *s2)
{
  int h = 0;

  for (; *s1 && *s2; s1++, s2++)
    if (*s1 != *s2)
      h++;

  return h;
}


PUBLIC int
vrna_hamming_distance_bound(const char  *s1,
                            const char  *s2,
                            int         boundary)
{
  int h = 0;

  for (; *s1 && *s2 && boundary; s1++, s2++, boundary--)
    if (*s1 != *s2)
      h++;

  return h;
}


PUBLIC void
vrna_seq_toRNA(char *sequence)
{
  unsigned int i;

  if (sequence) {
    for (i = 0; sequence[i]; i++) {
      if (sequence[i] == 'T')
        sequence[i] = 'U';

      if (sequence[i] == 't')
        sequence[i] = 'u';
    }
  }
}


PUBLIC void
vrna_seq_toupper(char *sequence)
{
  unsigned int i;

  if (sequence)
    for (i = 0; sequence[i]; i++)
      sequence[i] = toupper(sequence[i]);
}


PUBLIC char *
vrna_cut_point_insert(const char  *string,
                      int         cp)
{
  char  *ctmp;
  int   len;

  if (cp > 0) {
    len   = strlen(string);
    ctmp  = (char *)vrna_alloc((len + 2) * sizeof(char));
    /* first sequence */
    (void)strncpy(ctmp, string, cp - 1);
    /* spacer */
    ctmp[cp - 1] = '&';
    /* second sequence */
    (void)strcat(ctmp, string + cp - 1);
  } else {
    ctmp = strdup(string);
  }

  return ctmp;
}


PUBLIC char *
vrna_cut_point_remove(const char  *string,
                      int         *cp)
{
  char *pos, *copy = NULL;
  unsigned int len;

  *cp = -1;

  if (string) {
    len = strlen(string);
    copy = strdup(string);
    if ((pos = strchr(copy, '&'))) {
      *cp = (int)(pos - copy) + 1;
      if (*cp >= len)
        *cp = -1;

      if (strchr(pos + 1, '&'))
        vrna_message_error("more than one cut-point in input");

      for (; *pos; pos++)
        *pos = *(pos + 1);               /* splice out the & */
    }
  }

  return copy;
}


PUBLIC char **
vrna_strsplit(const char  *string,
              const char  *delimiter)
{
  char          delim[2], *ptr, *ptr2, *token, *save, **split;
  unsigned int  n;

  split = NULL;
  n     = 0;

  if (string) {
    if ((delimiter) && (*delimiter))
      delim[0] = *delimiter;
    else
      delim[0] = '&';

    delim[1] = '\0';

    /* copy string such that we can alter it via strtok() */
    ptr2 = strdup(string);

    /* count how many elements we'll extract */
    ptr = ptr2;

    while (*ptr++)
      if (*ptr == *delim)
        n++;

    /*
     * allocate (n + 1) + 1 elements in split list
     * n + 1 elements plus 1 additional element to indicate
     * the last element in split
     */
    split = (char **)vrna_alloc(sizeof(char *) * (n + 2));

    n     = 0;
    token = strtok_r(ptr2, delim, &save);

    while (token != NULL) {
      split[n++]  = vrna_strdup_printf("%s", token);
      token       = strtok_r(NULL, delim, &save);
    }

    split[n] = NULL;

    free(ptr2);
  }

  return split;
}


PUBLIC char *
vrna_seq_ungapped(const char *seq)
{
  char  *tmp_sequence, *b;
  int   i;

  tmp_sequence = NULL;

  if (seq) {
    tmp_sequence = strdup(seq);

    b = tmp_sequence;
    i = 0;
    do {
      if ((*b == '-') || (*b == '_') || (*b == '~') || (*b == '.'))
        continue;

      tmp_sequence[i] = *b;
      i++;
    } while (*(++b));

    tmp_sequence    = (char *)vrna_realloc(tmp_sequence, (i + 1) * sizeof(char));
    tmp_sequence[i] = '\0';
  }

  return tmp_sequence;
}


#ifndef VRNA_DISABLE_BACKWARD_COMPATIBILITY

/*###########################################*/
/*# deprecated functions below              #*/
/*###########################################*/

PUBLIC void
str_uppercase(char *sequence)
{
  vrna_seq_toupper(sequence);
}


PUBLIC void
str_DNA2RNA(char *sequence)
{
  vrna_seq_toRNA(sequence);
}


PUBLIC char *
random_string(int         l,
              const char  symbols[])
{
  return vrna_random_string(l, symbols);
}


PUBLIC int
hamming(const char  *s1,
        const char  *s2)
{
  return vrna_hamming_distance(s1, s2);
}


PUBLIC int
hamming_bound(const char  *s1,
              const char  *s2,
              int         boundary)
{
  return vrna_hamming_distance_bound(s1, s2, boundary);
}


#endif
