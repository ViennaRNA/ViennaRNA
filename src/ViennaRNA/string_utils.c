/*
                               string_utils.c

                 c  Ivo L Hofacker and Walter Fontana
                          Vienna RNA package
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

#include "ViennaRNA/utils.h"
#include "ViennaRNA/string_utils.h"

#ifndef HAVE_STRDUP
char *strdup(const char *s) {
  char *dup;

  dup = vrna_alloc(strlen(s)+1);
  strcpy(dup, s);
  return(dup);
}
#endif


PUBLIC char *
vrna_random_string(int l, const char symbols[]){

  char *r;
  int   i, rn, base;

  base = (int) strlen(symbols);
  r = (char *) vrna_alloc(sizeof(char)*(l+1));

  for (i = 0; i < l; i++) {
    rn = (int) (vrna_urn()*base);  /* [0, base-1] */
    r[i] = symbols[rn];
  }
  r[l] = '\0';
  return r;
}

/*-----------------------------------------------------------------*/

PUBLIC int
vrna_hamming_distance(const char *s1,
                      const char *s2){

  int h=0;

  for (; *s1 && *s2; s1++, s2++)
    if (*s1 != *s2) h++;
  return h;
}

PUBLIC int
vrna_hamming_distance_bound(const char *s1,
                            const char *s2,
                            int boundary){

  int h=0;

  for (; *s1 && *s2 && boundary; s1++, s2++, boundary--)
    if (*s1 != *s2) h++;
  return h;
}

PUBLIC  void
vrna_seq_toRNA(char *sequence){

  unsigned int i;
  if(sequence){
    for(i = 0; sequence[i]; i++){
      if(sequence[i] == 'T') sequence[i] = 'U';
      if(sequence[i] == 't') sequence[i] = 'u';
    }
  }
}

PUBLIC void
vrna_seq_toupper(char *sequence){

  unsigned int i;
  if(sequence){
    for(i=0;sequence[i];i++)
      sequence[i] = toupper(sequence[i]);
  }
}

PUBLIC char *
vrna_cut_point_insert(const char *string,
                      int cp){

  char *ctmp;
  int len;

  if(cp > 0){
    len = strlen(string);
    ctmp = (char *)vrna_alloc((len+2) * sizeof(char));
    /* first sequence */
    (void) strncpy(ctmp, string, cp-1);
    /* spacer */
    ctmp[cp-1] = '&';
    /* second sequence */
    (void) strcat(ctmp, string+cp-1);
  } else {
    ctmp = strdup(string);
  }
  return ctmp;
}

PUBLIC char *
vrna_cut_point_remove(const char *string,
                      int *cp){

  char *pos, *copy = NULL;

  *cp = -1;

  if(string){
    copy = (char *) vrna_alloc(strlen(string)+1);
    (void) sscanf(string, "%s", copy);
    pos = strchr(copy, '&');
    if (pos) {
      *cp = (int)(pos - copy) + 1;
      if (*cp >= strlen(copy)) *cp = -1;
      if (strchr(pos+1, '&')) vrna_message_error("more than one cut-point in input");
      for (;*pos;pos++) *pos = *(pos+1); /* splice out the & */
    }
  }

  return copy;
}

#ifdef  VRNA_BACKWARD_COMPAT

/*###########################################*/
/*# deprecated functions below              #*/
/*###########################################*/

PUBLIC void
str_uppercase(char *sequence){

  vrna_seq_toupper(sequence);
}

PUBLIC void
str_DNA2RNA(char *sequence){

  vrna_seq_toRNA(sequence);
}

PUBLIC char *
random_string(int l, const char symbols[]){

  return vrna_random_string(l, symbols);
}

PUBLIC int
hamming(const char *s1,
        const char *s2){

  return vrna_hamming_distance(s1, s2);
}

PUBLIC int
hamming_bound(const char *s1,
              const char *s2,
              int boundary){

  return vrna_hamming_distance_bound(s1, s2, boundary);
}

#endif
