/*
                               utils.c

                 c  Ivo L Hofacker and Walter Fontana
                          Vienna RNA package
*/
/* Last changed Time-stamp: <2008-11-25 16:34:36 ivo> */

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <errno.h>
#include <time.h>
#include <string.h>
#include "../config.h"
#include "utils.h"

#ifdef WITH_DMALLOC
#include "dmalloc.h"
#endif

#define PRIVATE  static
#define PUBLIC

/*@notnull@ @only@*/
PUBLIC unsigned short xsubi[3];

PRIVATE char  scale1[] = "....,....1....,....2....,....3....,....4";
PRIVATE char  scale2[] = "....,....5....,....6....,....7....,....8";
PRIVATE char  *inbuf = NULL;

PRIVATE char  *inbuf2 = NULL;
PRIVATE unsigned int typebuf2 = 0;

/* default values for the type/rtype stuff */
PRIVATE const char Law_and_Order[] = "_ACGUTXKI";

PRIVATE int rtype[8] = {0, 2, 1, 4, 3, 6, 5, 7};
PRIVATE int BP_pair[NBASES][NBASES]=
/* _  A  C  G  U  X  K  I */
{{ 0, 0, 0, 0, 0, 0, 0, 0},
 { 0, 0, 0, 0, 5, 0, 0, 5},
 { 0, 0, 0, 1, 0, 0, 0, 0},
 { 0, 0, 2, 0, 3, 0, 0, 0},
 { 0, 6, 0, 4, 0, 0, 0, 6},
 { 0, 0, 0, 0, 0, 0, 2, 0},
 { 0, 0, 0, 0, 0, 1, 0, 0},
 { 0, 6, 0, 0, 5, 0, 0, 0}};


PRIVATE int get_char_encoding(char c, model_detailsT *md);

/*-------------------------------------------------------------------------*/

PUBLIC void *space(unsigned size) {
  void *pointer;

  if ( (pointer = (void *) calloc(1, (size_t) size)) == NULL) {
#ifdef EINVAL
    if (errno==EINVAL) {
      fprintf(stderr,"SPACE: requested size: %d\n", size);
      nrerror("SPACE allocation failure -> EINVAL");
    }
    if (errno==ENOMEM)
#endif
      nrerror("SPACE allocation failure -> no memory");
  }
  return  pointer;
}

#ifdef WITH_DMALLOC
#define space(S) calloc(1,(S))
#endif

#undef xrealloc
/* dmalloc.h #define's xrealloc */
void *xrealloc (void *p, unsigned size) {
  if (p == 0)
    return space(size);
  p = (void *) realloc(p, size);
  if (p == NULL) {
#ifdef EINVAL
    if (errno==EINVAL) {
      fprintf(stderr,"xrealloc: requested size: %d\n", size);
      nrerror("xrealloc allocation failure -> EINVAL");
    }
    if (errno==ENOMEM)
#endif
      nrerror("xrealloc allocation failure -> no memory");
  }
  return p;
}

/*------------------------------------------------------------------------*/

PUBLIC void nrerror(const char message[])       /* output message upon error */
{
  fprintf(stderr, "ERROR: %s\n", message);
  exit(EXIT_FAILURE);
}

PUBLIC void warn_user(const char message[]){
  fprintf(stderr, "WARNING: %s\n", message);
}

/*------------------------------------------------------------------------*/
PUBLIC void init_rand(void)
{
  time_t t;
  (void) time(&t);
  xsubi[0] = xsubi[1] = xsubi[2] = (unsigned short) t;  /* lower 16 bit */
  xsubi[1] += (unsigned short) ((unsigned)t >> 6);
  xsubi[2] += (unsigned short) ((unsigned)t >> 12);
#ifndef HAVE_ERAND48
  srand((unsigned int) t);
#endif
}

/*------------------------------------------------------------------------*/

PUBLIC double urn(void)
     /* uniform random number generator; urn() is in [0,1] */
     /* uses a linear congruential library routine */
     /* 48 bit arithmetic */
{
#ifdef HAVE_ERAND48
  extern double erand48(unsigned short[]);
  return erand48(xsubi);
#else
  return ((double) rand())/RAND_MAX;
#endif
}

/*------------------------------------------------------------------------*/

PUBLIC int int_urn(int from, int to)
{
  return ( ( (int) (urn()*(to-from+1)) ) + from );
}

/*------------------------------------------------------------------------*/

PUBLIC void filecopy(FILE *from, FILE *to)
{
  int c;

  while ((c = getc(from)) != EOF) (void)putc(c, to);
}

/*-----------------------------------------------------------------*/

PUBLIC char *time_stamp(void)
{
  time_t  cal_time;

  cal_time = time(NULL);
  return ( ctime(&cal_time) );
}

/*-----------------------------------------------------------------*/

PUBLIC char *random_string(int l, const char symbols[])
{
  char *r;
  int   i, rn, base;

  base = (int) strlen(symbols);
  r = (char *) space(sizeof(char)*(l+1));

  for (i = 0; i < l; i++) {
    rn = (int) (urn()*base);  /* [0, base-1] */
    r[i] = symbols[rn];
  }
  r[l] = '\0';
  return r;
}

/*-----------------------------------------------------------------*/

PUBLIC int   hamming(const char *s1, const char *s2)
{
  int h=0;

  for (; *s1 && *s2; s1++, s2++)
    if (*s1 != *s2) h++;
  return h;
}

PUBLIC int   hamming_bound(const char *s1, const char *s2, int boundary)
{
  int h=0;

  for (; *s1 && *s2 && boundary; s1++, s2++, boundary--)
    if (*s1 != *s2) h++;
  return h;
}
/*-----------------------------------------------------------------*/

PUBLIC char *get_line(FILE *fp) /* reads lines of arbitrary length from fp */
{
  char s[512], *line, *cp;
  int len=0, size=0, l;
  line=NULL;
  do {
    if (fgets(s, 512, fp)==NULL) break;
    cp = strchr(s, '\n');
    if (cp != NULL) *cp = '\0';
    l = len + (int)strlen(s);
    if (l+1>size) {
      size = (int)((l+1)*1.2);
      line = (char *) xrealloc(line, size*sizeof(char));
    }
    strcat(line+len, s);
    len=l;
  } while(cp==NULL);

  return line;
}

PUBLIC int  skip_comment_lines(char **line){
  if((*line = get_line(stdin))==NULL) return -1;

  while((**line=='*')||(**line=='\0')){
    free(*line);
    if((*line = get_line(stdin))==NULL) return -1;
  }
  return 0;
}

PUBLIC  unsigned int get_input_line(char **string, unsigned int option){
  char  *line;
  int   i, l, r;

  /*
  * read lines until informative data appears or
  * report an error if anything goes wrong
  */
  if((line = get_line(stdin))==NULL) return VRNA_INPUT_ERROR;

  if(!(option & VRNA_INPUT_NOSKIP_COMMENTS))
    while((*line=='*')||(*line=='\0')){
      free(line);
      if((line = get_line(stdin))==NULL) return VRNA_INPUT_ERROR;
    }

  l = (int) strlen(line);

  /* break on '@' sign if not disabled */
  if(*line == '@'){
    free(line);
    return VRNA_INPUT_QUIT;
  }
  /* print line read if not disabled */
  /* if(!(option & VRNA_INPUT_NOPRINT)) printf("%s\n", line); */

  /* eliminate whitespaces at the end of the line read */
  if(!(option & VRNA_INPUT_NO_TRUNCATION)){
    for(i = l-1; i >= 0; i--){
      if      (line[i] == ' ')  continue;
      else if (line[i] == '\t') continue;
      else                      break;
    }
    line[(i >= 0) ? (i+1) : 0] = '\0';
  }

  if(*line == '>'){
    /* fasta header */
    /* alloc memory for the string */
    *string = (char *) space(sizeof(char) * (strlen(line) + 1));
    r = VRNA_INPUT_FASTA_HEADER;
    i = sscanf(line, ">%s", *string);
    if(i > 0){
      i       = (int)     strlen(*string);
      *string = (char *)  xrealloc(*string, (i+1)*sizeof(char));
      free(line);
      return r;
    }
    else{
      free(line);
      free(*string);
      *string = NULL;
      return VRNA_INPUT_ERROR;
    }
  }
  else{
    *string = strdup(line);
    free(line);
  }
  return VRNA_INPUT_MISC;
}

PUBLIC  unsigned int get_multi_input_line(char **string, unsigned int option){
  char  *line;
  int   i, l;
  int   state = 0;
  int   str_length = 0;

  line = (inbuf) ? inbuf : get_line(stdin);
  inbuf = NULL;
  do{

    /*
    * read lines until informative data appears or
    * report an error if anything goes wrong
    */
    if(!line) return VRNA_INPUT_ERROR;

    l = (int)strlen(line);

    /* eliminate whitespaces at the end of the line read */
    if(!(option & VRNA_INPUT_NO_TRUNCATION)){
      for(i = l-1; i >= 0; i--){
        if      (line[i] == ' ')  continue;
        else if (line[i] == '\t') continue;
        else                      break;
      }
      line[(i >= 0) ? (i+1) : 0] = '\0';
    }

    l           = (int)strlen(line);
    str_length  = (*string) ? (int) strlen(*string) : 0;

    switch(*line){
      case  '@':    /* user abort */
                    if(state) inbuf = line;
                    else      free(line);
                    return (state==2) ? VRNA_INPUT_CONSTRAINT : (state==1) ? VRNA_INPUT_SEQUENCE : VRNA_INPUT_QUIT;

      case  '\0':   /* empty line */
                    if(option & VRNA_INPUT_NOSKIP_BLANK_LINES){
                      if(state) inbuf = line;
                      else      free(line);
                      return (state==2) ? VRNA_INPUT_CONSTRAINT : (state==1) ? VRNA_INPUT_SEQUENCE : VRNA_INPUT_BLANK_LINE;
                    }
                    break;

      case  '#': case  '%': case  ';': case  '/': case  '*': case ' ':
                    /* comments */
                    if(option & VRNA_INPUT_NOSKIP_COMMENTS){
                      if(state) inbuf   = line;
                      else      *string = line;
                      return (state == 2) ? VRNA_INPUT_CONSTRAINT : (state==1) ? VRNA_INPUT_SEQUENCE : VRNA_INPUT_COMMENT;
                    }
                    break;

      case  '>':    /* fasta header */
                    if(state) inbuf   = line;
                    else      *string = line;
                    return (state==2) ? VRNA_INPUT_CONSTRAINT : (state==1) ? VRNA_INPUT_SEQUENCE : VRNA_INPUT_FASTA_HEADER;

      case  'x': case 'e': case 'l': case '&':   /* seems to be a constraint or line starting with second sequence for dimer calculations */
                    i = 1;
                    /* lets see if this assumption holds for the complete line */
                    while((line[i] == 'x') || (line[i] == 'e') || (line[i] == 'l')) i++;
                    /* lines solely consisting of 'x's, 'e's or 'l's will be considered as structure constraint */
                    
                    if(
                            ((line[i]>64) && (line[i]<91))  /* A-Z */
                        ||  ((line[i]>96) && (line[i]<123)) /* a-z */
                      ){
                      if(option & VRNA_INPUT_FASTA_HEADER){
                        /* are we in structure mode? Then we remember this line for the next round */
                        if(state == 2){ inbuf = line; return VRNA_INPUT_CONSTRAINT;}
                        else{
                          *string = (char *)xrealloc(*string, sizeof(char) * (str_length + l + 1));
                          strcpy(*string + str_length, line);
                          state = 1;
                        }
                        break;
                      }
                      /* otherwise return line read */
                      else{ *string = line; return VRNA_INPUT_SEQUENCE;}
                    }
                    /* mmmh? it really seems to be a constraint */
                    /* fallthrough */
      case  '<': case  '.': case  '|': case  '(': case ')': case '[': case ']': case '{': case '}': case ',': case '+':
                    /* seems to be a structure or a constraint */
                    /* either we concatenate this line to one that we read previously */
                    if(option & VRNA_INPUT_FASTA_HEADER){
                      if(state == 1){
                        inbuf = line;
                        return VRNA_INPUT_SEQUENCE;
                      }
                      else{
                        *string = (char *)xrealloc(*string, sizeof(char) * (str_length + l + 1));
                        strcpy(*string + str_length, line);
                        state = 2;
                      }
                    }
                    /* or we return it as it is */
                    else{
                      *string = line;
                      return VRNA_INPUT_CONSTRAINT;
                    }
                    break;
      default:      if(option & VRNA_INPUT_FASTA_HEADER){
                      /* are we already in sequence mode? */
                      if(state == 2){
                        inbuf = line;
                        return VRNA_INPUT_CONSTRAINT;
                      }
                      else{
                        *string = (char *)xrealloc(*string, sizeof(char) * (str_length + l + 1));
                        strcpy(*string + str_length, line);
                        state = 1;
                      }
                    }
                    /* otherwise return line read */
                    else{
                      *string = line;
                      return VRNA_INPUT_SEQUENCE;
                    }
    }
    free(line);
    line = get_line(stdin);
  }while(line);

  return (state==2) ? VRNA_INPUT_CONSTRAINT : (state==1) ? VRNA_INPUT_SEQUENCE : VRNA_INPUT_ERROR;
}

PUBLIC  unsigned int read_record(char **header, char **sequence, char ***rest, unsigned int options){
  unsigned int  input_type, return_type, tmp_type;
  int           rest_count;
  char          *input_string;

  rest_count    = 0;
  return_type   = tmp_type = 0;
  input_string  = *header = *sequence = NULL;
  *rest         = (char **)space(sizeof(char *));

  /* remove unnecessary option flags from options variable... */
  options &= ~VRNA_INPUT_FASTA_HEADER;

  /* read first input or last buffered input */
  if(typebuf2){
    input_type    = typebuf2;
    input_string  = inbuf2;
    typebuf2      = 0;
    inbuf2        = NULL;
  }
  else input_type  = get_multi_input_line(&input_string, options);

  if(input_type & (VRNA_INPUT_QUIT | VRNA_INPUT_ERROR)) return input_type;

  /* skip everything until we read either a fasta header or a sequence */
  while(input_type & (VRNA_INPUT_MISC | VRNA_INPUT_CONSTRAINT | VRNA_INPUT_BLANK_LINE)){
    free(input_string); input_string = NULL;
    input_type    = get_multi_input_line(&input_string, options);
    if(input_type & (VRNA_INPUT_QUIT | VRNA_INPUT_ERROR)) return input_type;
  }

  if(input_type & VRNA_INPUT_FASTA_HEADER){
    return_type  |= VRNA_INPUT_FASTA_HEADER; /* remember that we've read a fasta header */
    *header       = input_string;
    input_string  = NULL;
    /* get next data-block with fasta support if not explicitely forbidden by VRNA_INPUT_NO_SPAN */
    input_type  = get_multi_input_line(
                    &input_string,
                    ((options & VRNA_INPUT_NO_SPAN) ? 0 : VRNA_INPUT_FASTA_HEADER) | options
                  );
    if(input_type & (VRNA_INPUT_QUIT | VRNA_INPUT_ERROR)) return (return_type | input_type);
  }

  if(input_type & VRNA_INPUT_SEQUENCE){
    return_type  |= VRNA_INPUT_SEQUENCE; /* remember that we've read a sequence */
    *sequence     = input_string;
    input_string  = NULL;
  } else nrerror("sequence input missing");

  /* read the rest until we find user abort, EOF, new sequence or new fasta header */
  if(!(options & VRNA_INPUT_NO_REST)){
    options |= VRNA_INPUT_NOSKIP_COMMENTS; /* allow commetns to appear in rest output */
    tmp_type = VRNA_INPUT_QUIT | VRNA_INPUT_ERROR | VRNA_INPUT_SEQUENCE | VRNA_INPUT_FASTA_HEADER;
    if(options & VRNA_INPUT_NOSKIP_BLANK_LINES) tmp_type |= VRNA_INPUT_BLANK_LINE;
    while(!((input_type = get_multi_input_line(&input_string, options)) & tmp_type)){
      *rest = xrealloc(*rest, sizeof(char **)*(++rest_count + 1));
      (*rest)[rest_count-1] = input_string;
      input_string = NULL;
    }
    /*
    if(input_type & (VRNA_INPUT_QUIT | VRNA_INPUT_ERROR)) return input_type;
    */

    /*  finished reading everything...
    *   we now put the last line into the buffer if necessary
    *   since it should belong to the next record
    */
    inbuf2 = input_string;
    typebuf2 = input_type;
  }
  (*rest)[rest_count] = NULL;
  return (return_type);
}


/*-----------------------------------------------------------------*/

PUBLIC char *pack_structure(const char *struc) {
  /* 5:1 compression using base 3 encoding */
  int i,j,l,pi;
  unsigned char *packed;

  l = (int) strlen(struc);
  packed = (unsigned char *) space(((l+4)/5+1)*sizeof(unsigned char));

  j=i=pi=0;
  while (i<l) {
    register int p;
    for (p=pi=0; pi<5; pi++) {
      p *= 3;
      switch (struc[i]) {
      case '(':
      case '\0':
        break;
      case '.':
        p++;
        break;
      case ')':
        p += 2;
        break;
      default: nrerror("pack_structure: illegal charcter in structure");
      }
      if (i<l) i++;
    }
    packed[j++] = (unsigned char) (p+1); /* never use 0, so we can use
                                            strcmp()  etc. */
  }
  packed[j] = '\0';      /* for str*() functions */
  return (char *) packed;
}

PUBLIC char *unpack_structure(const char *packed) {
  /* 5:1 compression using base 3 encoding */
  int i,j,l;
  char *struc;
  unsigned const char *pp;
  char code[3] = {'(', '.', ')'};

  l = (int) strlen(packed);
  pp = (const unsigned char *) packed;
  struc = (char *) space((l*5+1)*sizeof(char));   /* up to 4 byte extra */

  for (i=j=0; i<l; i++) {
    register int p, c, k;

    p = (int) pp[i] - 1;
    for (k=4; k>=0; k--) {
      c = p % 3;
      p /= 3;
      struc[j+k] = code[c];
    }
    j += 5;
  }
  struc[j--] = '\0';
  while (struc[j] == '(') /* strip trailing ( */
    struc[j--] = '\0';

  return struc;
}

/*--------------------------------------------------------------------------*/

PUBLIC short *make_pair_table(const char *structure)
{
    /* returns array representation of structure.
       table[i] is 0 if unpaired or j if (i.j) pair.  */
   short i,j,hx;
   short length;
   short *stack;
   short *table;

   length = (short) strlen(structure);
   stack = (short *) space(sizeof(short)*(length+1));
   table = (short *) space(sizeof(short)*(length+2));
   table[0] = length;

   for (hx=0, i=1; i<=length; i++) {
      switch (structure[i-1]) {
       case '(':
         stack[hx++]=i;
         break;
       case ')':
         j = stack[--hx];
         if (hx<0) {
            fprintf(stderr, "%s\n", structure);
            nrerror("unbalanced brackets in make_pair_table");
         }
         table[i]=j;
         table[j]=i;
         break;
       default:   /* unpaired base, usually '.' */
         table[i]= 0;
         break;
      }
   }
   if (hx!=0) {
      fprintf(stderr, "%s\n", structure);
      nrerror("unbalanced brackets in make_pair_table");
   }
   free(stack);
   return(table);
}

PUBLIC short *make_pair_table_pk(const char *structure){
   short i,j,hx, hx2;
   short length;
   short *stack;
   short *stack2;
   short *table;

   length = (short) strlen(structure);
   stack  = (short *) space(sizeof(short)*(length+1));
   stack2 = (short *) space(sizeof(short)*(length+1));
   table  = (short *) space(sizeof(short)*(length+2));
   table[0] = length;

   for (hx=0, hx2=0, i=1; i<=length; i++) {
      switch (structure[i-1]) {
       case '(':
         stack[hx++]=i;
         break;
       case ')':
         j = stack[--hx];
         if (hx<0) {
            fprintf(stderr, "%s\n", structure);
            nrerror("unbalanced '()' brackets in make_pair_table_pk");
         }
         table[i]=j;
         table[j]=i;
         break;
       case '[':
         stack2[hx2++]=i;
         break;
       case ']':
         j = stack2[--hx2];
         if (hx2<0) {
            fprintf(stderr, "%s\n", structure);
            nrerror("unbalanced '[]' brackets in make_pair_table_pk");
         }
         table[i]=j;
         table[j]=i;
         break;
       default:   /* unpaired base, usually '.' */
         table[i]= 0;
         break;
      }
   }
   if (hx!=0) {
      fprintf(stderr, "%s\n", structure);
      nrerror("unbalanced '()' brackets in make_pair_table_pk");
   } else if (hx2!=0) {
      fprintf(stderr, "%s\n", structure);
      nrerror("unbalanced '[]' brackets in make_pair_table_pk");
   }
   free(stack);
   free(stack2);
   return(table);
}

PUBLIC short *make_pair_table_snoop(const char *structure)
{
    /* returns array representation of structure.
       table[i] is 0 if unpaired or j if (i.j) pair.  */
   short i,j,hx;
   short length;
   short *stack;
   short *table;

   length = (short) strlen(structure);
   stack = (short *) space(sizeof(short)*(length+1));
   table = (short *) space(sizeof(short)*(length+2));
   table[0] = length;

   for (hx=0, i=1; i<=length; i++) {
     switch (structure[i-1]) {
     case '<':
       stack[hx++]=i;
       break;
     case '>':
       j = stack[--hx];
       if (hx<0) {
         fprintf(stderr, "%s\n", structure);
         nrerror("unbalanced brackets in make_pair_table");
       }
       table[i]=j;
       table[j]=i;
       break;
     default:   /* unpaired base, usually '.' */
       table[i]= table[i];
       break;
     }
   }
   if (hx!=0) {
     fprintf(stderr, "%s\n", structure);
     nrerror("unbalanced brackets in make_pair_table");
   }
   free(stack);
   return table ;
}


PUBLIC short *alimake_pair_table(const char *structure)
{
    /* returns array representation of structure.
       table[i] is 0 if unpaired or j if (i.j) pair.  */
   short i,j,hx;
   short length;
   short *stack;
   short *table;

   length = (short) strlen(structure);
   stack = (short *) space(sizeof(short)*(length+1));
   table = (short *) space(sizeof(short)*(length+2));
   table[0] = length;

   for (hx=0, i=1; i<=length; i++) {
      switch (structure[i-1]) {
       case '(':
         stack[hx++]=i;
         break;
       case ')':
         j = stack[--hx];
         if (hx<0) {
            fprintf(stderr, "%s\n", structure);
            nrerror("unbalanced brackets in make_pair_table");
         }
         table[i]=j;
         table[j]=i;
         break;
       default:   /* unpaired base, usually '.' */
         table[i]= 0;
         break;
      }
   }
   for (hx=0, i=1; i<=length; i++) {
      switch (structure[i-1]) {
       case '<':
         stack[hx++]=i;
         break;
       case '>':
         j = stack[--hx];
         if (hx<0) {
            fprintf(stderr, "%s\n", structure);
            nrerror("unbalanced brackets in make_pair_table");
         }
         table[i]=j;
         table[j]=i;
         break;
       default:   /* unpaired base, usually '.' */
         table[i]= table[i];
         break;
      }
   }
   for (hx=0, i=1; i<=length; i++) {
     switch (structure[i-1]) {
     case '[':
       stack[hx++]=i;
       break;
     case ']':
       j = stack[--hx];
       if (hx<0) {
         fprintf(stderr, "%s\n", structure);
         nrerror("unbalanced brackets in make_pair_table");
       }
       table[i]=j;
       table[j]=i;
       break;
     default:   /* unpaired base, usually '.' */
       break;
     }
   }
   if (hx!=0) {
      fprintf(stderr, "%s\n", structure);
      nrerror("unbalanced brackets in make_pair_table");
   }
   free(stack);
   return(table);
}

PUBLIC short *copy_pair_table(const short *pt){
  short *table = (short *)space(sizeof(short) * (pt[0]+2));
  memcpy(table, pt, sizeof(short)*(pt[0]+2));
  return table;
}


PUBLIC int *make_loop_index_pt(short *pt){

  /* number each position by which loop it belongs to (positions start
     at 1) */
  int i,hx,l,nl;
  int length;
  int *stack = NULL;
  int *loop = NULL;

  length = pt[0];
  stack  = (int *) space(sizeof(int)*(length+1));
  loop   = (int *) space(sizeof(int)*(length+2));
  hx=l=nl=0;

  for (i=1; i<=length; i++) {
    if ((pt[i] != 0) && (i < pt[i])) { /* ( */
      nl++; l=nl;
      stack[hx++]=i;
    }
    loop[i]=l;

    if ((pt[i] != 0) && (i > pt[i])) { /* ) */
      --hx;
      if (hx>0)
        l = loop[stack[hx-1]];  /* index of enclosing loop   */
      else l=0;                 /* external loop has index 0 */
      if (hx<0) {
        nrerror("unbalanced brackets in make_pair_table");
      }
    }
  }
  loop[0] = nl;
  free(stack);
  return (loop);
}

/*---------------------------------------------------------------------------*/

PUBLIC int bp_distance(const char *str1, const char *str2)
{
  /* dist = {number of base pairs in one structure but not in the other} */
  /* same as edit distance with pair_open pair_close as move set */
   int dist;
   short i,l;
   short *t1, *t2;

   dist = 0;
   t1 = make_pair_table(str1);
   t2 = make_pair_table(str2);

   l = (t1[0]<t2[0])?t1[0]:t2[0];    /* minimum of the two lengths */

   for (i=1; i<=l; i++)
     if (t1[i]!=t2[i]) {
       if (t1[i]>i) dist++;
       if (t2[i]>i) dist++;
     }
   free(t1); free(t2);
   return dist;
}

#ifndef HAVE_STRDUP
char *strdup(const char *s) {
  char *dup;

  dup = space(strlen(s)+1);
  strcpy(dup, s);
  return(dup);
}
#endif

PUBLIC  void  print_tty_input_seq(void){
  print_tty_input_seq_str("Input string (upper or lower case)");
}

PUBLIC  void  print_tty_input_seq_str(const char *s){
  printf("\n%s; @ to quit\n", s);
  printf("%s%s\n", scale1, scale2);
  (void) fflush(stdout);
}

PUBLIC  void  str_DNA2RNA(char *sequence){
  unsigned int l, i;
  if(sequence != NULL){
    l = strlen(sequence);
    for(i = 0; i < l; i++){
      if(sequence[i] == 'T') sequence[i] = 'U';
      if(sequence[i] == 't') sequence[i] = 'u';
    }
  }
}

PUBLIC void str_uppercase(char *sequence){
  unsigned int l, i;
  if(sequence){
    l = strlen(sequence);
    for(i=0;i<l;i++)
      sequence[i] = toupper(sequence[i]);
  }
}

PUBLIC int *get_iindx(unsigned int length){
  int i;
  int *idx = (int *)space(sizeof(int) * (length+1));
  for (i=1; i <= length; i++)
    idx[i] = (((length + 1 - i) * (length - i))>>1) + length + 1;
  return idx;
}

PUBLIC int *get_indx(unsigned int length){
  unsigned int i;
  int *idx = (int *)space(sizeof(int) * (length+1));
  for (i = 1; i <= length; i++)
    idx[i] = (i*(i-1)) >> 1;        /* i(i-1)/2 */
  return idx;
}

PUBLIC char *extract_record_rest_structure( const char **lines,
                                            unsigned int length,
                                            unsigned int option){

  char *structure = NULL;
  int r, i, l, cl, stop;
  char *c;

  if(lines){
    for(r = i = stop = 0; lines[i]; i++){
      l   = (int)strlen(lines[i]);
      c   = (char *) space(sizeof(char) * (l+1));
      (void) sscanf(lines[i], "%s", c);
      cl  = (int)strlen(c);

      /* line commented out ? */
      if((*c == '#') || (*c == '%') || (*c == ';') || (*c == '/') || (*c == '*' || (*c == '\0'))){
        /* skip leading comments only, i.e. do not allow comments inside the constraint */
        if(!r)  continue;
        else    break;
      }

      /* append the structure part to the output */
      r += cl+1;
      structure = (char *)xrealloc(structure, r*sizeof(char));
      strcat(structure, c);
      free(c);
      /* stop if the assumed structure length has been reached */
      if((length > 0) && (r-1 == length)) break;
      /* stop if not allowed to read from multiple lines */
      if(!(option & VRNA_OPTION_MULTILINE)) break;
    }
  }
  return structure;
}




/* get a matrix containing the number of basepairs of a reference structure for each interval [i,j] with i<j
*  access it via iindx!!!
*/
PUBLIC unsigned int *make_referenceBP_array(short *reference_pt, unsigned int turn){
  unsigned int i,j,k,ij,length;
  int *iindx;
  unsigned int *array;
  unsigned int size;
  length = (unsigned int)reference_pt[0];
  size  = ((length+1)*(length+2))/2;
  iindx = get_iindx(length);
  array = (unsigned int *) space(sizeof(unsigned int)*size);    /* matrix containing number of basepairs of reference structure1 in interval [i,j] */;
  for (k=0; k<=turn; k++)
    for (i=1; i<=length-k; i++) {
      j=i+k;
      ij = iindx[i]-j;
      array[ij] = 0;
    }

  for (i = length-turn-1; i >= 1; i--)
    for (j = i+turn+1; j <= length; j++){
      int bps;
      ij = iindx[i]-j;
      bps = array[ij+1];
      if((i<=(unsigned int)reference_pt[j]) && ((unsigned int)reference_pt[j] < j))
        bps++;
      array[ij] = bps;
    }
  free(iindx);
  return array;
}

PUBLIC unsigned int *compute_BPdifferences(short *pt1, short *pt2, unsigned int turn){
  unsigned int *array;
  unsigned int n, size, i, j, ij, d;
  n = (unsigned int)pt1[0];
  size = ((n+1)*(n+2))/2;
  array = (unsigned int *)space(sizeof(unsigned int) * size);
  int *iindx = get_iindx(n);
  for(i = n - turn - 1; i>=1; i--){
    d = 0;
    for(j = i+turn+1; j <= n; j++){
      ij = iindx[i]-j;
      d = array[ij+1];
      if(pt1[j] != pt2[j]){
        if(i <= (unsigned int)pt1[j] && (unsigned int)pt1[j] < j){
          /* we got an additional base pair in reference structure 1 */
          d++;
        }
        if(i <= (unsigned int)pt2[j] && (unsigned int)pt2[j] < j){
          /* we got another base pair in reference structure 2 */
          d++;
        }
      }
      array[ij] = d;

    }
  }
  free(iindx);
  return array;
}

PUBLIC  void  fill_pair_matrices(model_detailsT *md){

  int i,j;

  /* nullify everything */
  for(i = 0;i <= MAXALPHA; i++)
    memset(md->pair[i], 0, (MAXALPHA + 1) * sizeof(int));

  memset(md->alias, 0, (MAXALPHA + 1) * sizeof(short));

  /* start setting actual base pair type encodings */
  switch(md->energy_set){
    case  0:    for(i = 0; i < 5; i++)
                  md->alias[i] = (short) i;

                md->alias[5] = 3; /* X <-> G */
                md->alias[6] = 2; /* K <-> C */
                md->alias[7] = 0; /* I <-> default base '@' */

                for(i = 0; i < NBASES; i++)
                    for(j = 0; j < NBASES; j++)
                      md->pair[i][j] = BP_pair[i][j];

                if(md->noGU)
                  md->pair[3][4] = md->pair[4][3] = 0;

                if(md->nonstandards != NULL) {  /* allow nonstandard bp's (encoded by type=7) */
                   for(i = 0; i < (int)strlen(md->nonstandards); i += 2)
                      md->pair[get_char_encoding(md->nonstandards[i], md)]
                        [get_char_encoding(md->nonstandards[i+1], md)] = 7;
                }

                break;

    case 1:     for(i = 1; i < MAXALPHA;){
                  md->alias[i++] = 3;  /* A <-> G */
                  md->alias[i++] = 2;  /* B <-> C */
                }
                for(i = 1; i < MAXALPHA; i++){
                  md->pair[i][i+1] = 2;    /* AB <-> GC */
                  i++;
                  md->pair[i][i-1] = 1;    /* BA <-> CG */
                }

                break;

    case 2:     for(i = 1; i < MAXALPHA;){
                  md->alias[i++] = 1;  /* A <-> A*/
                  md->alias[i++] = 4;  /* B <-> U */
                }
                for(i = 1; i < MAXALPHA; i++){
                  md->pair[i][i+1] = 5;    /* AB <-> AU */
                  i++;
                  md->pair[i][i-1] = 6;    /* BA <-> UA */
                }

                break;

    case 3:     for(i = 1; i < MAXALPHA - 2; ){
                  md->alias[i++] = 3;  /* A <-> G */
                  md->alias[i++] = 2;  /* B <-> C */
                  md->alias[i++] = 1;  /* C <-> A */
                  md->alias[i++] = 4;  /* D <-> U */
                }
                for(i = 1; i < MAXALPHA - 2; i++){
                  md->pair[i][i+1] = 2;    /* AB <-> GC */
                  i++;
                  md->pair[i][i-1] = 1;    /* BA <-> CG */
                  i++;
                  md->pair[i][i+1] = 5;    /* CD <-> AU */
                  i++;
                  md->pair[i][i-1] = 6;    /* DC <-> UA */
                }

                break;

    default:    nrerror("Which energy_set are YOU using??");
                break;
  }

  /* set the reverse base pair types */
  for(i = 0; i <= MAXALPHA; i++){
    for(j = 0; j <= MAXALPHA; j++){
      md->rtype[md->pair[i][j]] = md->pair[j][i];
    }
  }

  /* was used for energy_set == 0
  for(i = 0; i < NBASES; i++)
      for(j = 0; j < NBASES; j++)
       md->rtype[md->pair[i][j]] = md->pair[j][i];
  */
}

PUBLIC  short *get_sequence_encoding( const char *sequence,
                                      short type,
                                      model_detailsT *md){

  unsigned int i,l = (unsigned int)strlen(sequence);
  short         *S = (short *) space(sizeof(short)*(l+2));

  switch(type){
    /* standard encoding as always used for S */
    case 0:   for(i=1; i<=l; i++) /* make numerical encoding of sequence */
                S[i]= (short) get_char_encoding(toupper(sequence[i-1]), md);
              S[l+1] = S[1];
              S[0] = (short) l;
              break;
    /* encoding for mismatches of nostandard bases (normally used for S1) */
    case 1:   for(i=1; i<=l; i++)
                S[i] = md->alias[(short)get_char_encoding(toupper(sequence[i-1]), md)];
              S[l+1] = S[1];
              S[0] = S[l];
              break;
  }

  return S;
}

PRIVATE int get_char_encoding(char c, model_detailsT *md){
  /* return numerical representation of base used e.g. in pair[][] */
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

PUBLIC  char  get_encoded_char(int enc, model_detailsT *md){
  if(md->energy_set > 0) return (char)enc + 'A' - 1;
  else return (char)Law_and_Order[enc];
}

PUBLIC  char  *get_ptypes(const short *S,
                          model_detailsT *md,
                          unsigned int idx_type){

  char *ptype;
  int n,i,j,k,l,*idx;

  n     = S[0];
  ptype = (char *)space(sizeof(char)*((n*(n+1))/2+2));
  idx   = (idx_type) ? get_iindx(n) : get_indx(n);

  if(idx_type){
    for (k=1; k<n-TURN; k++)
      for (l=1; l<=2; l++) {
        int type,ntype=0,otype=0;
        i=k; j = i+TURN+l; if (j>n) continue;
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
  } else {
    for (k=1; k<n-TURN; k++)
      for (l=1; l<=2; l++) {
        int type,ntype=0,otype=0;
        i=k; j = i+TURN+l; if (j>n) continue;
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
  }
  free(idx);
  return ptype;
}
