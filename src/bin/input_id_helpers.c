#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <limits.h>
#include "ViennaRNA/utils.h"

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


#if 0

char *
ID_generate(char        *ID_default,
            int         ID_auto_switch,
            const char  *ID_prefix,
            const char  *ID_delimiter,
            int         ID_digits,
            long int    ID_number,
            int         ID_full)
{
  char *ID, *tmp;

  ID = tmp = NULL;

  if (ID_auto_switch) {
    /* generate ID automatically */
    ID = vrna_strdup_printf("%s%s%0*ld", ID_prefix, ID_delimiter, ID_digits, ID_number);
  } else {
    /* compose ID if possible */
    /* Use already existing string as ID, e.g. from FASTA header */
    if ((ID_default) && (ID_default[0] != '\0')) {
      if (ID_full) {
        /* use full default ID */
        tmp = strdup(ID_default);
      } else {
        /* extract first word from default ID */
        tmp = (char *)vrna_alloc(sizeof(char) * (strlen(ID_default) + 1));
        (void)sscanf(ID_default, "%s", tmp);
        tmp = (char *)vrna_realloc(tmp, sizeof(char) * (strlen(tmp) + 1));
      }
    }

    if (ID_prefix) {
      /* prepend ID with prefix */
      if (tmp) {
        ID = vrna_strdup_printf("%s%s%s", ID_prefix, ID_delimiter, tmp);
        free(tmp);
      } else {
        ID = vrna_strdup_printf("%s", ID_prefix);
      }
    } else {
      ID = tmp;
    }
  }

  return ID;
}


#endif
