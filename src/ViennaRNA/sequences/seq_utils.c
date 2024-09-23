/*
 *    ViennaRNA/sequences/seq_utils.c
 *
 *    c  Ivo L Hofacker, Ronny Lorenz
 *       Vienna RNA package
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
#include "ViennaRNA/sequences/utils.h"


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


PUBLIC char *
vrna_DNA_complement(const char *sequence)
{
  char    *complement, *ptr;
  size_t  n;

  complement = NULL;

  if (sequence) {
    n           = strlen(sequence);
    complement  = (char *)vrna_alloc(sizeof(char) * (n + 1));
    /* copy the input string */
    complement  = memcpy(complement, sequence, sizeof(char) * n);

    /* complement characters */
    for (ptr = complement; *ptr; ptr++) {
      switch (*ptr) {
        case 'A':
          *ptr = 'T';
          break;

        case 'a':
          *ptr = 't';
          break;

        case 'C':
          *ptr = 'G';
          break;

        case 'c':
          *ptr = 'g';
          break;

        case 'G':
          *ptr = 'C';
          break;

        case 'g':
          *ptr = 'c';
          break;

        case 'T': /* fall through */
        case 'U':
          *ptr = 'A';
          break;

        case 't': /* fall through */
        case 'u':
          *ptr = 'a';
          break;

        default:
          break;
      }
    }

    complement[n] = '\0';
  }

  return complement;
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

/*
 * ###########################################
 * # deprecated functions below              #
 * ###########################################
 */

PUBLIC void
str_DNA2RNA(char *sequence)
{
  vrna_seq_toRNA(sequence);
}


#endif
