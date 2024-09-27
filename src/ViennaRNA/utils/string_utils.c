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
#include "ViennaRNA/utils/log.h"
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
  va_list copy;
  va_copy(copy, argp);

  /* the default is to use vasprintf() if available */
  r = vasprintf(&result, format, copy);

  va_end(copy); /* Each va_start() or va_copy() needs a va_end() */
#endif

  /* check for any memory allocation error indicated by r == -1 */
  if (r == -1) {
    vrna_log_warning("vrna_strdup_printf: memory allocation failure!");
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
  int     r;
  size_t  l1, l2;
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
    vrna_log_warning("vrna_strcat_printf: memory allocation failure!");
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
vrna_seq_toupper(char *sequence)
{
  unsigned int i;

  if (sequence)
    for (i = 0; sequence[i]; i++)
      sequence[i] = toupper(sequence[i]);
}


PUBLIC void
vrna_seq_reverse(char *sequence)
{
  if (sequence) {
    char *p1 = sequence;
    char *p2 = sequence + strlen(sequence) - 1;

    while (p1 < p2) {
      char tmp = *p1;
      *p1++ = *p2;
      *p2-- = tmp;
    }
  }
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
  char          *pos, *copy = NULL;
  unsigned int  len;

  *cp = -1;

  if (string) {
    len   = strlen(string);
    copy  = strdup(string);
    if ((pos = strchr(copy, '&'))) {
      *cp = (int)(pos - copy) + 1;
      if (*cp >= len)
        *cp = -1;

      if (strchr(pos + 1, '&')) {
        vrna_log_error("more than one cut-point in input");
        free(copy);
        return NULL;
      }

      for (; *pos; pos++)
        *pos = *(pos + 1);               /* splice out the & */
    }
  }

  return copy;
}


PUBLIC unsigned int
vrna_strtrim(char         *string,
             const char   *delimiters,
             unsigned int keep,
             unsigned int options)
{
  char          delim_ws[7] = {
    32, 9, 10, 11, 12, 13, 0
  };
  const char    *delim, *ptrd;
  char          *str_end, *ptr, *ptr_out, *ptr_start, *ptr_end;
  unsigned int  ret, hits;

  ret = 0;

  if (string) {
    if ((delimiters) &&
        (*delimiters))
      delim = delimiters;
    else
      delim = &(delim_ws[0]);

    /* find first non-delimiter position */
    for (ptr_start = string; *ptr_start != '\0'; ptr_start++) {
      /* check whether we've found a delimiter */
      for (ptrd = delim; *ptrd != '\0'; ptrd++)
        if (*ptrd == *ptr_start)
          break;

      /* abort if we didn't find any of the delimiter characters */
      if (*ptrd == '\0')
        break;
    }

    /* find last non-delimiter position */
    for (ptr_end = ptr = ptr_start; *ptr != '\0'; ptr++) {
      for (ptrd = delim; *ptrd != '\0'; ptrd++)
        if (*ptrd == *ptr)
          break;

      if (*ptrd == '\0')
        ptr_end = ptr;
    }

    ptr_end++;

    str_end = ptr_out = ptr;

    if (options & VRNA_TRIM_LEADING) {
      ptr = ptr_start -
            sizeof(char) * keep;

      if (ptr < string)
        ptr = string;

      /* adjust start and end pointer positions */
      ptr_start -= ptr - string;
      ptr_end   -= ptr - string;

      /* actually remove leading delimiters */
      for (ptr_out = string; ptr < ptr_start; ptr++)
        if (options & VRNA_TRIM_SUBST_BY_FIRST)
          *(ptr_out++) = delim[0];
        else
          *(ptr_out++) = *ptr;

      for (; *ptr != '\0'; ptr++)
        *(ptr_out++) = *ptr;

      *ptr_out = '\0';
    }

    if (options & VRNA_TRIM_IN_BETWEEN) {
      hits = 0;

      /* remove in-between delimiters */
      for (ptr_out = ptr = ptr_start; ptr < ptr_end; ptr++) {
        for (ptrd = delim; *ptrd != '\0'; ptrd++)
          if (*ptrd == *ptr)
            break;

        /* did we find a delimiter? */
        if (*ptrd != '\0') {
          if (hits++ < keep) {
            if (options & VRNA_TRIM_SUBST_BY_FIRST)
              *ptr_out = delim[0];
            else
              *ptr_out = *ptr;

            ptr_out++;
          }
        } else {
          hits      = 0;
          *ptr_out  = *ptr;
          ptr_out++;
        }
      }

      /* adjust end pointer position */
      ptr_end -= ptr - ptr_out;

      /* shift trailing part to correct position */
      for (; *ptr != '\0'; ptr++)
        *(ptr_out++) = *ptr;

      *ptr_out = '\0';
    }

    if (options & VRNA_TRIM_TRAILING) {
      hits = 0;

      /* seek to end */
      for (ptr_out = ptr = ptr_end; *ptr != '\0'; ptr++) {
        if (hits++ < keep) {
          if (options & VRNA_TRIM_SUBST_BY_FIRST)
            *ptr_out = delim[0];
          else
            *ptr_out = *ptr;

          ptr_out++;
        }
      }

      *ptr_out = '\0';
    }

    return (str_end - ptr_out) / sizeof(char);
  }

  return ret;
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
#ifdef _WIN32
# ifndef __MINGW32__
    token = strtok_s(ptr2, delim, &save);
# else
    token = strtok_r(ptr2, delim, &save);
# endif
#else
    token = strtok_r(ptr2, delim, &save);
#endif

    while (token != NULL) {
      split[n++]  = vrna_strdup_printf("%s", token);
#ifdef _WIN32
# ifndef __MINGW32__
      token = strtok_s(NULL, delim, &save);
# else
      token = strtok_r(NULL, delim, &save);
# endif
#else
      token = strtok_r(NULL, delim, &save);
#endif
    }

    split[n] = NULL;

    free(ptr2);
  }

  return split;
}


PUBLIC char *
vrna_strjoin(const char **strings,
             const char *delimiter)
{
  char    *s = NULL;
  size_t  n, offset, *lengths, mem_strings, total_length;

  if (strings) {
    total_length  = 0;
    mem_strings   = 32;
    lengths       = (size_t *)vrna_alloc(sizeof(size_t) * mem_strings);

    for (n = 0; strings[n]; n++) {
      if (n == mem_strings) {
        mem_strings += 32;
        lengths      = (size_t *)vrna_realloc(lengths, sizeof(size_t) * mem_strings);
      }

      lengths[n]    = strlen(strings[n]);
      total_length += lengths[n];
    }

    if ((delimiter) && (*delimiter))
      total_length += (n - 1);

    /* finally, glue the strings together */
    s = (char *)vrna_alloc(sizeof(char) * (total_length + 1));

    for (offset = 0, n = 0; strings[n]; n++) {
      memcpy(s + offset, strings[n], sizeof(char) * lengths[n]);
      offset += lengths[n];

      if ((delimiter) &&
          (*delimiter) &&
          (strings[n + 1]))
        s[offset++] = *delimiter;
    }

    s[total_length] = '\0';

    free(lengths);
  }

  return s;
}


#if 0
PUBLIC char *
vrna_strsplice(const char   *string,
               const char   *delimiter,
               unsigned int **positions,
               unsigned int options)
{
  char *result = NULL;
  
  if (string) {
    if (delimiter) {
      if (options & VRNA_STRSPLICE_IN){
        if (positions) {
          /* count how many more characters we require for the fully spliced string */
          for (size_t n = 0; positions[n] != 0; n++);

          size_t dl = strlen(delimiter);
          size_t l  = strlen(string);

          result = (char *)vrna_alloc(sizeof(char) * (l + dl * n + 1));

          /* finally, construct the spliced sequence */
          size_t start = 0;
          size_t end   = 0;
          size_t last_pos = 0;
          /* handle first case separately */
          memcpy(result, string, sizeof(char) * ((*positions)[0] - 1));
          memcpy(result + (*positions)[0] - 1, delimiter, sizeof(char) * dl);
          start += (*positions)[0] - 1;
          end   += (*positions)[0] - 1 + dl;

          for (size_t i = 1; i < n; i++) {
            memcpy(result + end, string + start, sizeof(char) * positions
          }

        } else {
          result = strdup(string);
        }
      } else if (options & VRNA_STRSPLICE_OUT) {
        
      }
    } else {
      /* no delimiter specified, so we don't need to do anything */
      result = strdup(string);
      if ((options & VRNA_STRSPLICE_OUT) &&
          (positions)) {
        *positions = (unsigned int *)vrna_alloc(sizeof(unsigned int));
        (*positions)[0] = 0;
      }
    }
  }

  return result;
}

#endif


PUBLIC size_t *
vrna_strchr(const char  *str,
            int          c,
            size_t       n)
{
  size_t  str_len, pos_num, *positions  = NULL;
  char    *ptr, *ptr2;

  if (str) {
    str_len   = strlen(str);
    pos_num   = 0;
    positions = vrna_alloc(sizeof(size_t) * (str_len + 1));
    ptr       = NULL;

    if (n == 0)
      n = str_len;

    /* find first occurence */
    ptr = strchr(str, (int)c);

    if (ptr != NULL) {
      positions[++pos_num] = ptr - str + 1;
      ptr2 = ptr + 1;

      /* find further occurrences */
      do {
        ptr = strchr((const char*)ptr2, (int)c);
        if (ptr) {
          positions[++pos_num] = ptr - str + 1;
          ptr2 = ptr + 1;
        }
      } while ((*ptr2) &&
               (ptr != NULL) &&
               (pos_num <= n));
    }

    positions[0] = pos_num;

    positions = vrna_realloc(positions, sizeof(size_t) * (pos_num + 1));
  }

  return positions;
}

#ifndef VRNA_DISABLE_BACKWARD_COMPATIBILITY

/*
 * ###########################################
 * # deprecated functions below              #
 * ###########################################
 */

PUBLIC void
str_uppercase(char *sequence)
{
  vrna_seq_toupper(sequence);
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
