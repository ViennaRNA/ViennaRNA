/*
                               string_utils.c

                 c  Ivo L Hofacker and Walter Fontana
                          Vienna RNA package
*/

#include <config.h>
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

PRIVATE const char Law_and_Order[] = "_ACGUTXKI";

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

PUBLIC short *
vrna_seq_encode(const char *sequence,
                vrna_md_t *md){

  unsigned int  i, l;
  short         *S = vrna_seq_encode_simple(sequence, md);

  l = (unsigned int)strlen(sequence);

  for(i=1; i<=l; i++)
    S[i] = md->alias[S[i]];

  S[l+1] = S[1];
  S[0] = S[l];

  return S;
}

PUBLIC short *
vrna_seq_encode_simple( const char *sequence,
                        vrna_md_t *md){

  unsigned int i,l = (unsigned int)strlen(sequence);
  short         *S = (short *) vrna_alloc(sizeof(short)*(l+2));

  for(i=1; i<=l; i++) /* make numerical encoding of sequence */
    S[i]= (short) vrna_nucleotide_encode(toupper(sequence[i-1]), md);

  S[l+1] = S[1];
  S[0] = (short) l;

  return S;
}

PUBLIC  int
vrna_nucleotide_encode( char c,
                        vrna_md_t *md){

  /* return numerical representation of nucleotide used e.g. in vrna_md_t.pair[][] */
  int code;
  if (md->energy_set>0) code = (int) (c-'A')+1;
  else {
    const char *pos;
    pos = strchr(Law_and_Order, c);
    if (pos==NULL) code=0;
    else code = (int) (pos-Law_and_Order);
    if (code>5) code = 0;
    if (code>4) code--; /* make T and U equivalent */
  }
  return code;
}

PUBLIC  char
vrna_nucleotide_decode( int enc,
                        vrna_md_t *md){

  if(md->energy_set > 0)
    return (char)enc + 'A' - 1;
  else
    return (char)Law_and_Order[enc];
}

PUBLIC void
vrna_aln_encode(const char *sequence,
                    short **S_p,
                    short **s5_p,
                    short **s3_p,
                    char **ss_p,
                    unsigned short **as_p,
                    vrna_md_t *md){

  unsigned  int   i,l;
  unsigned  short p;

  l     = strlen(sequence);

  (*s5_p)   = (short *)         vrna_alloc((l + 2) * sizeof(short));
  (*s3_p)   = (short *)         vrna_alloc((l + 2) * sizeof(short));
  (*as_p)  = (unsigned short *)vrna_alloc((l + 2) * sizeof(unsigned short));
  (*ss_p)   = (char *)          vrna_alloc((l + 2) * sizeof(char));

  /* make numerical encoding of sequence */
  (*S_p)    = vrna_seq_encode_simple(sequence, md);

  (*s5_p)[0] = (*s5_p)[1] = 0;

  if(md->oldAliEn){
    /* use alignment sequences in all energy evaluations */
    (*ss_p)[0]=sequence[0];
    for(i=1; i<l; i++){
      (*s5_p)[i] = (*S_p)[i-1];
      (*s3_p)[i] = (*S_p)[i+1];
      (*ss_p)[i] = sequence[i];
      (*as_p)[i] = i;
    }
    (*ss_p)[l]   = sequence[l];
    (*as_p)[l]   = l;
    (*s5_p)[l]   = (*S_p)[l-1];
    (*s3_p)[l]   = 0;
    (*S_p)[l+1]  = (*S_p)[1];
    (*s5_p)[1]   = 0;
    if(md->circ){
      (*s5_p)[1]   = (*S_p)[l];
      (*s3_p)[l]   = (*S_p)[1];
      (*ss_p)[l+1] = (*S_p)[1];
    }
  }
  else{
    if(md->circ){
      for(i=l; i>0; i--){
        char c5;
        c5 = sequence[i-1];
        if ((c5=='-')||(c5=='_')||(c5=='~')||(c5=='.')) continue;
        (*s5_p)[1] = (*S_p)[i];
        break;
      }
      for (i=1; i<=l; i++) {
        char c3;
        c3 = sequence[i-1];
        if ((c3=='-')||(c3=='_')||(c3=='~')||(c3=='.')) continue;
        (*s3_p)[l] = (*S_p)[i];
        break;
      }
    }
    else  (*s5_p)[1]=(*s3_p)[l]=0;

    for(i=1,p=0; i<=l; i++){
      char c5;
      c5 = sequence[i-1];
      if ((c5=='-')||(c5=='_')||(c5=='~')||(c5=='.'))
        (*s5_p)[i+1]=(*s5_p)[i];
      else { /* no gap */
        (*ss_p)[p++]=sequence[i-1]; /*start at 0!!*/
        (*s5_p)[i+1]=(*S_p)[i];
      }
      (*as_p)[i]=p;
    }
    for (i=l; i>=1; i--) {
      char c3;
      c3 = sequence[i-1];
      if ((c3=='-')||(c3=='_')||(c3=='~')||(c3=='.'))
        (*s3_p)[i-1]=(*s3_p)[i];
      else
        (*s3_p)[i-1]=(*S_p)[i];
    }
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
