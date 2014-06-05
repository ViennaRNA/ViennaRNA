/*
                               utils.c

                 c  Ivo L Hofacker and Walter Fontana
                          Vienna RNA package
*/

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <errno.h>
#include <time.h>
#include <string.h>
#include <sys/types.h>
#include <stdint.h>

/* for getpid() we need some distinction between UNIX and Win systems */
#ifdef _WIN32
#include <windows.h>
#define getpid() GetCurrentProcessId() /* rename windows specific getpid function */
#else
#include <unistd.h>
#endif

#include "../config.h"
#include "ViennaRNA/utils.h"

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


PRIVATE int   get_char_encoding(char c, model_detailsT *md);
PRIVATE char  *wrap_get_ptypes(const short *S, model_detailsT *md);  /* provides backward compatibility for old ptypes array in pf computations */

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

PRIVATE uint32_t
rj_mix( uint32_t a,
        uint32_t b,
        uint32_t c){

/*
  This is Robert Jenkins' 96 bit Mix function

  we use it to produce a more diverse seed for our random number
  generators. E.g.:
  
  seed = rj_mix(clock(), time(NULL), getpid());

  original comments on that function can be found below
*/


/*
--------------------------------------------------------------------
mix -- mix 3 32-bit values reversibly.
For every delta with one or two bits set, and the deltas of all three
  high bits or all three low bits, whether the original value of a,b,c
  is almost all zero or is uniformly distributed,
* If mix() is run forward or backward, at least 32 bits in a,b,c
  have at least 1/4 probability of changing.
* If mix() is run forward, every bit of c will change between 1/3 and
  2/3 of the time.  (Well, 22/100 and 78/100 for some 2-bit deltas.)
mix() was built out of 36 single-cycle latency instructions in a 
  structure that could supported 2x parallelism, like so:
      a -= b; 
      a -= c; x = (c>>13);
      b -= c; a ^= x;
      b -= a; x = (a<<8);
      c -= a; b ^= x;
      c -= b; x = (b>>13);
      ...
  Unfortunately, superscalar Pentiums and Sparcs can't take advantage 
  of that parallelism.  They've also turned some of those single-cycle
  latency instructions into multi-cycle latency instructions.  Still,
  this is the fastest good hash I could find.  There were about 2^^68
  to choose from.  I only looked at a billion or so.
--------------------------------------------------------------------
*/

  a=a-b;  a=a-c;  a=a^(c >> 13);
  b=b-c;  b=b-a;  b=b^(a << 8); 
  c=c-a;  c=c-b;  c=c^(b >> 13);
  a=a-b;  a=a-c;  a=a^(c >> 12);
  b=b-c;  b=b-a;  b=b^(a << 16);
  c=c-a;  c=c-b;  c=c^(b >> 5);
  a=a-b;  a=a-c;  a=a^(c >> 3);
  b=b-c;  b=b-a;  b=b^(a << 10);
  c=c-a;  c=c-b;  c=c^(b >> 15);
  return c;
}

/*------------------------------------------------------------------------*/
PUBLIC void init_rand(void)
{

  uint32_t seed = rj_mix(clock(), time(NULL), getpid());

  xsubi[0] = xsubi[1] = xsubi[2] = (unsigned short) seed;  /* lower 16 bit */
  xsubi[1] += (unsigned short) ((unsigned)seed >> 6);
  xsubi[2] += (unsigned short) ((unsigned)seed >> 12);
#ifndef HAVE_ERAND48
  srand((unsigned int) seed);
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

/*--------------------------------------------------------------------------*/


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

PUBLIC void
get_sequence_encoding_gapped( const char *sequence,
                              short **S_p,
                              short **s5_p,
                              short **s3_p,
                              char **ss_p,
                              unsigned short **as_p,
                              model_detailsT *md){

  unsigned  int   i,l;
  unsigned  short p;

  l     = strlen(sequence);

  (*s5_p)   = (short *)         space((l + 2) * sizeof(short));
  (*s3_p)   = (short *)         space((l + 2) * sizeof(short));
  (*as_p)  = (unsigned short *)space((l + 2) * sizeof(unsigned short));
  (*ss_p)   = (char *)          space((l + 2) * sizeof(char));
  (*S_p)    = (short *)         space((l + 2) * sizeof(short));


  (*S_p)[0]  = (short) l;
  (*s5_p)[0] = (*s5_p)[1] = 0;

  /* make numerical encoding of sequence */
  for(i=1; i<=l; i++){
    (*S_p)[i] = (short) get_char_encoding(toupper(sequence[i-1]), md);
  }

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
vrna_get_ptypes(const short *S,
                model_detailsT *md){

  char *ptype;
  int n,i,j,k,l,*idx;
  int min_loop_size = md->min_loop_size;

  n     = S[0];
  ptype = (char *)space(sizeof(char)*((n*(n+1))/2+2));
  idx   = get_indx(n);

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

PRIVATE char *
wrap_get_ptypes(const short *S,
                model_detailsT *md){

  char *ptype;
  int n,i,j,k,l,*idx;

  n     = S[0];
  ptype = (char *)space(sizeof(char)*((n*(n+1))/2+2));
  idx   = get_iindx(n);

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
  free(idx);
  return ptype;
}


PUBLIC int *
get_pscores(const short *const* S,
            const char **AS,
            int n_seq,
            float **distance_matrix,
            model_detailsT *md){

  int *pscore, *indx;


  /* calculate co-variance bonus for each pair depending on  */
  /* compensatory/consistent mutations and incompatible seqs */
  /* should be 0 for conserved pairs, >0 for good pairs      */
#define NONE -10000 /* score for forbidden pairs */
#define UNIT 100
#define MINPSCORE -2 * UNIT
  int n,i,j,k,l,s, max_span;

  int olddm[7][7]={{0,0,0,0,0,0,0}, /* hamming distance between pairs */
                  {0,0,2,2,1,2,2} /* CG */,
                  {0,2,0,1,2,2,2} /* GC */,
                  {0,2,1,0,2,1,2} /* GU */,
                  {0,1,2,2,0,2,1} /* UG */,
                  {0,2,2,1,2,0,2} /* AU */,
                  {0,2,2,2,1,2,0} /* UA */};

  float **dm;
  n       = S[0][0];  /* length of seqs */
  pscore  = (int *) space(sizeof(int)*((n*(n+1))/2+2));
  indx    = get_indx(n);

  if(distance_matrix){
    dm = distance_matrix;
  }
  else { /*use usual matrix*/
    dm=(float **)space(7*sizeof(float*));
    for (i=0; i<7;i++) {
      dm[i]=(float *)space(7*sizeof(float));
      for (j=0; j<7; j++)
        dm[i][j] = (float) olddm[i][j];
    }
  }

  max_span = md->max_bp_span;
  if((max_span < TURN+2) || (max_span > n))
    max_span = n;
  for (i=1; i<n; i++) {
    for (j=i+1; (j<i+TURN+1) && (j<=n); j++)
      pscore[indx[j]+i] = NONE;
    for (j=i+TURN+1; j<=n; j++) {
      int pfreq[8]={0,0,0,0,0,0,0,0};
      double score;
      for (s=0; s<n_seq; s++) {
        int type;
        if (S[s][i]==0 && S[s][j]==0) type = 7; /* gap-gap  */
        else {
          if ((AS[s][i] == '~')||(AS[s][j] == '~')) type = 7;
          else type = md->pair[S[s][i]][S[s][j]];
        }
        pfreq[type]++;
      }
      if (pfreq[0]*2+pfreq[7]>n_seq) { pscore[indx[j]+i] = NONE; continue;}
      for (k=1,score=0; k<=6; k++) /* ignore pairtype 7 (gap-gap) */
        for (l=k; l<=6; l++)
          score += pfreq[k]*pfreq[l]*dm[k][l];
      /* counter examples score -1, gap-gap scores -0.25   */
      pscore[indx[j]+i] = md->cv_fact *
        ((UNIT*score)/n_seq - md->nc_fact*UNIT*(pfreq[0] + pfreq[7]*0.25));

      if((j - i + 1) > max_span){
        pscore[indx[j]+i] = NONE;
      }
    }
  }

  if (md->noLP) /* remove unwanted pairs */
    for (k=1; k<n-TURN-1; k++)
      for (l=1; l<=2; l++) {
        int type,ntype=0,otype=0;
        i=k; j = i+TURN+l;
        type = pscore[indx[j]+i];
        while ((i>=1)&&(j<=n)) {
          if ((i>1)&&(j<n)) ntype = pscore[indx[j+1]+i-1];
          if ((otype<md->cv_fact*MINPSCORE)&&(ntype<md->cv_fact*MINPSCORE))  /* too many counterexamples */
            pscore[indx[j]+i] = NONE; /* i.j can only form isolated pairs */
          otype =  type;
          type  = ntype;
          i--; j++;
        }
      }


  /*free dm */
  for (i=0; i<7;i++) {
    free(dm[i]);
  }
  free(dm);
  free(indx);
  return pscore;
}

/*###########################################*/
/*# deprecated functions below              #*/
/*###########################################*/

PUBLIC char *
get_ptypes( const short *S,
            model_detailsT *md,
            unsigned int idx_type){

  if(idx_type)
    return wrap_get_ptypes(S, md);
  else
    return vrna_get_ptypes(S, md);
}

