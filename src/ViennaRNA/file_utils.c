/*
 *                             file_utils.c
 *
 *               c  Ronny Lorenz
 *               Vienna RNA package
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
#include <libgen.h>

#include "ViennaRNA/utils.h"
#include "ViennaRNA/file_utils.h"

#define PRIVATE  static
#define PUBLIC

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
#ifdef _WIN32
#ifdef __MINGW32__
#include <direct.h>
#endif
#define DIRSEPC '\\'
#define DIRSEPS "\\"
#else
#define DIRSEPC '/'
#define DIRSEPS "/"
#endif

/*
 #################################
 # PRIVATE FUNCTION DECLARATIONS #
 #################################
 */
PRIVATE int is_absolute_path(const char *p);


/*
 #################################
 # BEGIN OF FUNCTION DEFINITIONS #
 #################################
 */
PUBLIC void
vrna_file_copy(FILE *from,
               FILE *to)
{
  int c;

  while ((c = getc(from)) != EOF)
    (void)putc(c, to);
}


PUBLIC char *
vrna_read_line(FILE *fp)
{
  /* reads lines of arbitrary length from fp */

  char  s[512], *line, *cp;
  int   len = 0, size = 0, l;

  line = NULL;
  do {
    if (fgets(s, 512, fp) == NULL)
      break;

    cp = strchr(s, '\n');
    if (cp != NULL)
      *cp = '\0';

    l = len + (int)strlen(s);
    if (l + 1 > size) {
      size  = (int)((l + 1) * 1.2);
      line  = (char *)vrna_realloc(line, size * sizeof(char));
    }

    strcat(line + len, s);
    len = l;
  } while (cp == NULL);

  return line;
}


PUBLIC int
vrna_mkdir_p(const char *path)
{
  struct stat sb;
  char        *slash, *ptr;
  int         done = 0;

  if (!is_absolute_path(path))
    ptr = vrna_strdup_printf(".%c%s", DIRSEPC, path);
  else
    ptr = strdup(path);

  slash = ptr;

  while (!done) {
    slash += strspn(slash, DIRSEPS);
    slash += strcspn(slash, DIRSEPS);

    done    = (*slash == '\0');
    *slash  = '\0';

    if (stat(ptr, &sb)) {
#ifdef _WIN32
      if (errno != ENOENT || (_mkdir(ptr) &&
                              errno != EEXIST)) {
#else
      if (errno != ENOENT || (mkdir(ptr, 0777) &&
                              errno != EEXIST)) {
#endif
        vrna_message_warning("Can't create directory %s", ptr);
        free(ptr);
        return -1;
      }
    } else if (!S_ISDIR(sb.st_mode)) {
      vrna_message_warning("File exists but is not a directory %s: %s", ptr, strerror(ENOTDIR));
      free(ptr);
      return -1;
    }

    *slash = DIRSEPC;
  }

  free(ptr);
  return 0;
}


PUBLIC char *
vrna_basename(const char *path)
{
  char  *name, *p, *ptr;
  int   pos;

  name = NULL;

  if (path) {
    ptr = strrchr(path, DIRSEPC);

    if (ptr && (*(ptr + 1) != '\0'))
      name = strdup(ptr + 1);
    else if (!ptr)
      name = strdup(path);
  }

  return name;
}


PUBLIC char *
vrna_dirname(const char *path)
{
  char  *name, *p, *ptr;
  int   pos;

  name = NULL;

  if (path) {
    if (!is_absolute_path(path))
      ptr = vrna_strdup_printf(".%c%s", DIRSEPC, path);
    else
      ptr = strdup(path);

    pos = (int)strlen(ptr);
    p   = ptr + pos;

    do  /* remove part after last separator */
      *p = '\0';
    while ((--p > ptr) && (*p != DIRSEPC));

    if (p > ptr)
      name = ptr;
  }

  return name;
}


PUBLIC char *
vrna_filename_sanitize(const char *name,
                       const char *replacement)
{
  if (name) {
    const char    *ptr, *start, *illegal_chars;
    char          *sanitized_name;
    unsigned int  i, n;

    illegal_chars   = "\\/?%*:|\"<> ";
    sanitized_name  = (char *)vrna_alloc(sizeof(char) * (strlen(name) + 1));
    start           = name;
    i               = 0;
    while (ptr = strpbrk(start, illegal_chars)) {
      /* find illegal chars */
      strncpy(sanitized_name + i, start, ptr - start);
      i += ptr - start;
      if (replacement && (*replacement))
        sanitized_name[i++] = *replacement;

      start = ptr + 1; /* skip invalid character */
    }
    /* copy remaining part */
    if (start < (name + strlen(name))) {
      unsigned int diff = name - start + strlen(name);
      strncpy(sanitized_name + i, start, diff);
      i += diff;
    }

    /* resize the output string to actual requirements */
    sanitized_name    = (char *)vrna_realloc(sanitized_name, sizeof(char) * (i + 1));
    sanitized_name[i] = '\0';

    /* check for reserved unix file names */
    if ((!strcmp(sanitized_name, ".")) || (!strcmp(sanitized_name, ".."))) {
      sanitized_name    = (char *)vrna_realloc(sanitized_name, sizeof(char));
      sanitized_name[0] = '\0';
    }

    /* check for length restrictions */
    n = strlen(sanitized_name);
    if (n > 255) {
      char *suff = NULL;
      /* try to leave file suffix, i.e. everything after last dot '.', intact */
      if ((suff = strrchr(sanitized_name, '.')) && (sanitized_name + n - suff < 255)) {
        unsigned int n_suff = sanitized_name + n - suff;
        memmove(sanitized_name + (255 - n_suff), sanitized_name + n - n_suff, sizeof(char) * n_suff);
      }

      sanitized_name      = (char *)vrna_realloc(sanitized_name, sizeof(char) * 256);
      sanitized_name[255] = '\0';
    }

    /* finally, return the sanitized file name */
    return sanitized_name;
  } else {
    return NULL;
  }
}


#ifdef _WIN32
PRIVATE int
is_drive_char(const char c)
{
  if ((c >= 'a' && c <= 'z') || (c >= 'A' && c <= 'Z'))
    return 1;

  return 0;
}


#endif


PRIVATE int
is_absolute_path(const char *p)
{
  if (*p == DIRSEPC)
    return 1;

#ifdef _WIN32
  if (is_drive_char((const char)*p) && (strlen(p) > 3))
    if ((*(p + 1) == ':') && ((*(p + 2) == '\\') || (*(p + 2) == '/')))
      return 1;

#endif
  return 0;
}
