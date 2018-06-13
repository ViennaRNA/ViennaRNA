#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <limits.h>

#include "ViennaRNA/utils/basic.h"
#include "ViennaRNA/utils/strings.h"

#include "input_id_helpers.h"


struct id_data {
  char      *name;
  int       auto_id;
  char      *prefix;
  char      *delimiter;
  int       digits;
  long int  number;
};


void
ID_number_increase(long int   *num,
                   const char *name)
{
  if (*num == LONG_MAX) {
    vrna_message_warning("%s ID number overflow, beginning with 1 (again)!", name);
    *num = 1;
  } else {
    (*num)++;
  }
}


struct id_data *
init_id_data(const char *name,
             const char *default_prefix,
             const char *default_delimiter,
             int        default_digits,
             long       default_start)
{
  struct id_data *dat;

  dat = (struct id_data *)vrna_alloc(sizeof(struct id_data));

  dat->name       = (name) ? strdup(name) : NULL;
  dat->auto_id    = 0;
  dat->prefix     = (default_prefix) ? strdup(default_prefix) : NULL;
  dat->delimiter  = (default_delimiter) ? strdup(default_delimiter) : NULL;
  dat->digits     = default_digits;
  dat->number     = default_start;

  return dat;
}


void
free_id_data(struct id_data *dat)
{
  if (dat) {
    free(dat->name);
    free(dat->prefix);
    free(dat->delimiter);
    free(dat);
  }
}


void
set_auto_id(struct id_data  *dat,
            int             sw)
{
  if (dat)
    dat->auto_id = sw;
}


int
get_auto_id(struct id_data *dat)
{
  if (dat)
    return dat->auto_id;

  return 0;
}


void
set_id_prefix(struct id_data  *dat,
              const char      *prefix)
{
  if (dat) {
    free(dat->prefix);
    dat->prefix = strdup(prefix);
  }
}


void
set_id_delim(struct id_data *dat,
             const char     *delim)
{
  if (dat) {
    free(dat->delimiter);
    dat->delimiter = strdup(delim);
  }
}


const char *
get_id_delim(struct id_data *dat)
{
  if (dat)
    return (const char *)dat->delimiter;

  return NULL;
}


void
set_id_digits(struct id_data  *dat,
              int             digits)
{
  if (dat)
    dat->digits = digits;
}


void
set_id_start(struct id_data *dat,
             long int       start)
{
  if (dat)
    dat->number = start;
}


void
set_next_id(char            **ID_default,
            struct id_data  *dat)
{
  if (dat) {
    if (dat->number == LONG_MAX) {
      vrna_message_warning("%s ID number overflow, beginning with 1 (again)!", dat->name);
      dat->number = 1;
    }

    if ((ID_default) && (dat->auto_id)) {
      /* generate an ID */
      free(*ID_default);
      *ID_default = vrna_strdup_printf("%s%s%0*ld",
                                       dat->prefix,
                                       dat->delimiter,
                                       dat->digits,
                                       dat->number);
    }

    dat->number++;
  }
}


long int
get_current_id(struct id_data *dat)
{
  if (dat)
    return dat->number - 1;

  return -1;
}


char *
fileprefix_from_id(const char     *id,
                   struct id_data *dat,
                   int            ID_full)
{
  char *prefix = NULL;

  if ((id) && (id[0] != '\0') && (dat)) {
    if ((ID_full) || (dat->auto_id)) {
      /* use full ID */
      prefix = strdup(id);
    } else {
      /* use first word from default ID */
      prefix = (char *)vrna_alloc(sizeof(char) * (strlen(id) + 1));
      (void)sscanf(id, "%s", prefix);
      prefix = (char *)vrna_realloc(prefix, sizeof(char) * (strlen(prefix) + 1));
    }
  }

  return prefix;
}


char *
fileprefix_from_id_alifold(const char     *id,
                           struct id_data *dat,
                           int            continuous_names)
{
  char *prefix = NULL;

  if (dat) {
    if (dat->number == LONG_MAX) {
      vrna_message_warning("%s ID number overflow, beginning with 1 (again)!", dat->name);
      dat->number = 1;
    }

    if (id && (!dat->auto_id)) {
      prefix = strdup(id);
    } else if (dat->auto_id || (dat->number > 1) || continuous_names) {
      prefix = vrna_strdup_printf("%s%s%0*ld",
                                  dat->prefix,
                                  dat->delimiter,
                                  dat->digits,
                                  dat->number);
    }

    dat->number++;
  }

  return prefix;
}
