/*
    alphabet.c
    
    Code for handling nucleotide and base pair alphabet
    
    Part of the ViennaRNA Package
*/

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <limits.h>

#include "ViennaRNA/utils.h"
#include "ViennaRNA/alphabet.h"

/*
  For now, we neglect all non-standard nucleotides in an input sequence, i.e. only
  ACGTUN is allowed.

  However, the standard nucleotide ambiguity code table would allow for many more:

  A = Adenylic acid
  C = Cytidylic acid
  G = Guanylic acid
  T = Thymidylic acid
  U = Uridylic acid
  I = Inosylic acid
  R = A or G = puRine
  Y = C or T = pYrimidine
  K = G or T = Keto
  M = A or C = aMino
  S = G or C = Strong base pair
  W = A or T = Weak base pair
  B = not A (G or C or T)
  D = not C (A or G or T)
  H = not G (A or C or T)
  V = not T/U (A or C or G)
  N = aNy base  (by convention, X is used for unknown amino acids, N for unknown nucleotides)

  For the future, we aim to accept all of the above codes.
*/

/*
#################################
# PRIVATE VARIABLES             #
#################################
*/
PRIVATE const char Law_and_Order[] = "_ACGUTXKI";

/*
#################################
# PRIVATE FUNCTION DECLARATIONS #
#################################
*/
PRIVATE char  *wrap_get_ptypes(const short *S, vrna_md_t *md);  /* provides backward compatibility for old ptypes array in pf computations */

/*
#################################
# BEGIN OF FUNCTION DEFINITIONS #
#################################
*/

PUBLIC unsigned int
vrna_sequence_length_max(unsigned int options){

  if(options & VRNA_OPTION_WINDOW)
    return (unsigned int)INT_MAX;

/*
  return (unsigned int)sqrt((double)INT_MAX);
*/
  /*
      many functions in RNAlib still rely on the sequence length
      at pos 0 in the integer encoded sequence array S. Since this
      encoding is stored in a short * array, the maximum length
      of any sequence is SHRT_MAX
  */
  return (unsigned int)SHRT_MAX;
}


PUBLIC int
vrna_nucleotide_IUPAC_identity( char nt,
                                char mask){

  char n1,n2,*p;

  p   = NULL;
  n1  = toupper(nt);
  n2  = toupper(mask);

  switch(n1){
    case 'A': p = strchr("ARMWDHVN", n2);
              break;
    case 'C': p = strchr("CYMSBHVN", n2);
              break;
    case 'G': p = strchr("GRKSBDVN", n2);
              break;
    case 'T': p = strchr("TYKWBDHN", n2);
              break;
    case 'U': p = strchr("UYKWBDHN", n2);
              break;
    case 'I': p = strchr("IN", n2);
              break;
    case 'R': p = strchr("AGR", n2);
              break;
    case 'Y': p = strchr("CTUY", n2);
              break;
    case 'K': p = strchr("GTUK", n2);
              break;
    case 'M': p = strchr("ACM", n2);
              break;
    case 'S': p = strchr("GCS", n2);
              break;
    case 'W': p = strchr("ATUW", n2);
              break;
    case 'B': p = strchr("GCTBU", n2);
              break;
    case 'D': p = strchr("AGTUD", n2);
              break;
    case 'H': p = strchr("ACTUH", n2);
              break;
    case 'V': p = strchr("ACGV", n2);
              break;
    case 'N': p = strchr("ACGTUN", n2);
              break;
  }

  return (p) ? 1 : 0;
}


PUBLIC char *
vrna_ptypes(const short *S,
                vrna_md_t *md){

  char *ptype;
  int n,i,j,k,l,*idx;
  int min_loop_size = md->min_loop_size;

  n     = S[0];

  if((unsigned int)n > vrna_sequence_length_max(VRNA_OPTION_DEFAULT)){
    vrna_message_warning("vrna_ptypes@alphabet.c: sequence length of %d exceeds addressable range", n);
    return NULL;
  }

  ptype = (char *)vrna_alloc(sizeof(char)*((n*(n+1))/2+2));
  idx   = vrna_idx_col_wise(n);

  for (k=1; k<n-min_loop_size; k++)
    for (l=1; l<=2; l++) {
      int type,ntype=0,otype=0;
      i=k; j = i+min_loop_size+l; if (j>n) continue;
      type = md->pair[S[i]][S[j]];
      while ((i>=1)&&(j<=n)) {
        if ((i>1)&&(j<n)) ntype = md->pair[S[i-1]][S[j+1]];
        if (md->noLP && (!otype) && (!ntype))
          type = 0; /* i.j can only form isolated pairs */
        ptype[idx[j]+i] = (char) type;
        otype =  type;
        type  = ntype;
        i--; j++;
      }
    }
  free(idx);
  return ptype;
}

PUBLIC short *
vrna_seq_encode(const char *sequence,
                vrna_md_t *md){

  unsigned int  i, l;
  short         *S = NULL;
  
  if(sequence && md){
    S = vrna_seq_encode_simple(sequence, md);

    l = (unsigned int)strlen(sequence);

    for(i=1; i<=l; i++)
      S[i] = md->alias[S[i]];

    S[l+1] = S[1];
    S[0] = S[l];
  }

  return S;
}

PUBLIC short *
vrna_seq_encode_simple( const char *sequence,
                        vrna_md_t *md){

  unsigned int  i, l;
  short         *S = NULL;

  if(sequence && md){
    l = (unsigned int)strlen(sequence);
    S = (short *) vrna_alloc(sizeof(short)*(l+2));

    for(i=1; i<=l; i++) /* make numerical encoding of sequence */
      S[i]= (short) vrna_nucleotide_encode(toupper(sequence[i-1]), md);

    S[l+1] = S[1];
    S[0] = (short) l;
  }

  return S;
}

PUBLIC  int
vrna_nucleotide_encode( char c,
                        vrna_md_t *md){

  /* return numerical representation of nucleotide used e.g. in vrna_md_t.pair[][] */
  int code = -1;

  if(md){
    if (md->energy_set>0) code = (int) (c-'A')+1;
    else {
      const char *pos;
      pos = strchr(Law_and_Order, c);
      if (pos==NULL) code=0;
      else code = (int) (pos-Law_and_Order);
      if (code>5) code = 0;
      if (code>4) code--; /* make T and U equivalent */
    }
  }

  return code;
}

PUBLIC  char
vrna_nucleotide_decode( int enc,
                        vrna_md_t *md){

  if(md){
    if(md->energy_set > 0)
      return (char)enc + 'A' - 1;
    else
      return (char)Law_and_Order[enc];
  } else {
    return (char)0;
  }
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

PRIVATE char *
wrap_get_ptypes(const short *S,
                vrna_md_t *md){

  char *ptype;
  int n,i,j,k,l,*idx;

  n     = S[0];
  ptype = (char *)vrna_alloc(sizeof(char)*((n*(n+1))/2+2));
  idx   = vrna_idx_row_wise(n);
  int min_loop_size = md->min_loop_size;

  for (k=1; k<n-min_loop_size; k++)
    for (l=1; l<=2; l++) {
      int type,ntype=0,otype=0;
      i=k; j = i+min_loop_size+l; if (j>n) continue;
      type = md->pair[S[i]][S[j]];
      while ((i>=1)&&(j<=n)) {
        if ((i>1)&&(j<n)) ntype = md->pair[S[i-1]][S[j+1]];
        if (md->noLP && (!otype) && (!ntype))
          type = 0; /* i.j can only form isolated pairs */
        ptype[idx[i]-j] = (char) type;
        otype =  type;
        type  = ntype;
        i--; j++;
      }
    }
  free(idx);
  return ptype;
}

#ifdef  VRNA_BACKWARD_COMPAT

/*###########################################*/
/*# deprecated functions below              #*/
/*###########################################*/

PUBLIC char *
get_ptypes( const short *S,
            vrna_md_t *md,
            unsigned int idx_type){

  if(S){
    if((unsigned int)S[0] > vrna_sequence_length_max(VRNA_OPTION_DEFAULT)){
      vrna_message_warning("get_ptypes@alphabet.c: sequence length of %d exceeds addressable range", (int)S[0]);
      return NULL;
    }

    if(idx_type)
      return wrap_get_ptypes(S, md);
    else
      return vrna_ptypes(S, md);
  } else {
    return NULL;
  }
}

#endif

