/*
                               utils.c

                 c  Ivo L Hofacker and Walter Fontana
                          Vienna RNA package
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

/* for getpid() we need some distinction between UNIX and Win systems */
#ifdef _WIN32
#include <windows.h>
#define getpid() GetCurrentProcessId() /* rename windows specific getpid function */
#else
#include <unistd.h>
#endif

#include "ViennaRNA/utils.h"

#ifdef WITH_DMALLOC
#include "dmalloc.h"
#endif

#define PRIVATE  static
#define PUBLIC

#define EXIT_ON_ERROR

#include "ViennaRNA/color_output.inc"

/*
#################################
# GLOBAL VARIABLES              #
#################################
*/
/*@notnull@ @only@*/
PUBLIC unsigned short xsubi[3];


/*
#################################
# PRIVATE VARIABLES             #
#################################
*/
PRIVATE char  scale1[] = "....,....1....,....2....,....3....,....4";
PRIVATE char  scale2[] = "....,....5....,....6....,....7....,....8";
PRIVATE char  *inbuf = NULL;
PRIVATE char  *inbuf2 = NULL;
PRIVATE unsigned int typebuf2 = 0;


/*
#################################
# PRIVATE FUNCTION DECLARATIONS #
#################################
*/
PRIVATE uint32_t  rj_mix( uint32_t a, uint32_t b, uint32_t c);


/*
#################################
# BEGIN OF FUNCTION DEFINITIONS #
#################################
*/

#ifndef WITH_DMALLOC
/* include the following two functions only if not including <dmalloc.h> */

PUBLIC void *
vrna_alloc(unsigned size){

  void *pointer;

  if ( (pointer = (void *) calloc(1, (size_t) size)) == NULL) {
#ifdef EINVAL
    if (errno==EINVAL) {
      fprintf(stderr,"vrna_alloc: requested size: %d\n", size);
      vrna_message_error("Memory allocation failure -> EINVAL");
    }
    if (errno==ENOMEM)
#endif
      vrna_message_error("Memory allocation failure -> no memory");
  }
  return  pointer;
}

PUBLIC void *
vrna_realloc(void *p, unsigned size){

  if (p == NULL)
    return vrna_alloc(size);
  p = (void *) realloc(p, size);
  if (p == NULL) {
#ifdef EINVAL
    if (errno==EINVAL) {
      fprintf(stderr,"vrna_realloc: requested size: %d\n", size);
      vrna_message_error("vrna_realloc allocation failure -> EINVAL");
    }
    if (errno==ENOMEM)
#endif
      vrna_message_error("vrna_realloc allocation failure -> no memory");
  }
  return p;
}

#endif

/*------------------------------------------------------------------------*/

PUBLIC void
vrna_message_error(const char *format, ...){

  va_list args;
  va_start(args, format);
  vrna_message_verror(format, args);
  va_end(args);
}


PUBLIC void
vrna_message_verror(const char *format, va_list args){

#ifndef VRNA_WITHOUT_TTY_COLORS
  if(isatty(fileno(stderr))){
    fprintf(stderr, ANSI_COLOR_RED_B "ERROR: " ANSI_COLOR_RESET ANSI_COLOR_BRIGHT);
    vfprintf(stderr, format, args);
    fprintf(stderr, ANSI_COLOR_RESET "\n");
  } else {
#endif
    fprintf(stderr, "ERROR: ");
    vfprintf(stderr, format, args);
    fprintf(stderr, "\n");
#ifndef VRNA_WITHOUT_TTY_COLORS
  }
#endif

#ifdef EXIT_ON_ERROR
  exit(EXIT_FAILURE);
#endif
}


PUBLIC void
vrna_message_warning(const char *format, ...){

  va_list args;
  va_start(args, format);
  vrna_message_vwarning(format, args);
  va_end(args);
}


PUBLIC void
vrna_message_vwarning(const char *format, va_list args){

#ifndef VRNA_WITHOUT_TTY_COLORS
  if(isatty(fileno(stderr))){
    fprintf(stderr, ANSI_COLOR_MAGENTA_B "WARNING: " ANSI_COLOR_RESET ANSI_COLOR_BRIGHT);
    vfprintf(stderr, format, args);
    fprintf(stderr, ANSI_COLOR_RESET "\n");
  } else {
#endif
    fprintf(stderr, "WARNING: ");
    vfprintf(stderr, format, args);
    fprintf(stderr, "\n");
#ifndef VRNA_WITHOUT_TTY_COLORS
  }
#endif
}


PUBLIC void
vrna_message_info(FILE *fp, const char *format, ...){

  va_list args;
  va_start(args, format);
  vrna_message_vinfo(fp, format, args);
  va_end(args);
}


PUBLIC void
vrna_message_vinfo(FILE *fp, const char *format, va_list args){

  if(!fp)
    fp = stdout;

#ifndef VRNA_WITHOUT_TTY_COLORS
  if(isatty(fileno(fp))){
    fprintf(fp, ANSI_COLOR_BLUE_B);
    vfprintf(fp, format, args);
    fprintf(fp, ANSI_COLOR_RESET "\n");
  } else {
#endif
    vfprintf(fp, format, args);
    fprintf(fp, "\n");
#ifndef VRNA_WITHOUT_TTY_COLORS
  }
#endif
}


/*------------------------------------------------------------------------*/
PUBLIC void
vrna_init_rand(void){

  uint32_t seed = rj_mix(clock(), time(NULL), getpid());

  xsubi[0] = xsubi[1] = xsubi[2] = (unsigned short) seed;  /* lower 16 bit */
  xsubi[1] += (unsigned short) ((unsigned)seed >> 6);
  xsubi[2] += (unsigned short) ((unsigned)seed >> 12);
#ifndef HAVE_ERAND48
  srand((unsigned int) seed);
#endif
}

/*------------------------------------------------------------------------*/

/* uniform random number generator; vrna_urn() is in [0,1] */
/* uses a linear congruential library routine */
/* 48 bit arithmetic */
PUBLIC double
vrna_urn(void){

#ifdef HAVE_ERAND48
  extern double erand48(unsigned short[]);
  return erand48(xsubi);
#else
  return ((double) rand())/RAND_MAX;
#endif
}

/*------------------------------------------------------------------------*/

PUBLIC int
vrna_int_urn(int from, int to){

  return ( ( (int) (vrna_urn()*(to-from+1)) ) + from );
}

/*------------------------------------------------------------------------*/

/*-----------------------------------------------------------------*/

PUBLIC char *
vrna_time_stamp(void){

  time_t  cal_time;

  cal_time = time(NULL);
  return ( ctime(&cal_time) );
}

/*-----------------------------------------------------------------*/

PUBLIC  unsigned int get_input_line(char **string, unsigned int option){
  char  *line;
  int   i, l, r;

  /*
  * read lines until informative data appears or
  * report an error if anything goes wrong
  */
  if((line = vrna_read_line(stdin))==NULL) return VRNA_INPUT_ERROR;

  if(!(option & VRNA_INPUT_NOSKIP_COMMENTS))
    while((*line=='*')||(*line=='\0')){
      free(line);
      if((line = vrna_read_line(stdin))==NULL) return VRNA_INPUT_ERROR;
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
    *string = (char *) vrna_alloc(sizeof(char) * (strlen(line) + 1));
    r = VRNA_INPUT_FASTA_HEADER;
    i = sscanf(line, ">%s", *string);
    if(i > 0){
      i       = (int)     strlen(*string);
      *string = (char *)  vrna_realloc(*string, (i+1)*sizeof(char));
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


PUBLIC void
vrna_message_input_seq_simple(void){

  vrna_message_input_seq("Input string (upper or lower case)");
}

PUBLIC  void
vrna_message_input_seq(const char *s){
#ifndef VRNA_WITHOUT_TTY_COLORS
  if(isatty(fileno(stdout))){
    printf("\n" ANSI_COLOR_CYAN "%s; @ to quit" ANSI_COLOR_RESET "\n", s);
    printf(ANSI_COLOR_BRIGHT "%s%s" ANSI_COLOR_RESET "\n", scale1, scale2);
  } else {
#endif
    printf("\n%s; @ to quit\n", s);
    printf("%s%s\n", scale1, scale2);
#ifndef VRNA_WITHOUT_TTY_COLORS
  }
#endif
  (void) fflush(stdout);
}

PUBLIC int *
vrna_idx_row_wise(unsigned int length){

  int i;
  int *idx = (int *)vrna_alloc(sizeof(int) * (length+1));
  for (i=1; i <= length; i++)
    idx[i] = (((length + 1 - i) * (length - i)) / 2) + length + 1;
  return idx;
}

PUBLIC int *
vrna_idx_col_wise(unsigned int length){

  unsigned int i;
  int *idx = (int *)vrna_alloc(sizeof(int) * (length+1));
  for (i = 1; i <= length; i++)
    idx[i] = (i*(i-1)) / 2; 
  return idx;
}


/*
#################################
# STATIC helper functions below #
#################################
*/
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

#ifdef  VRNA_BACKWARD_COMPAT

/*###########################################*/
/*# deprecated functions below              #*/
/*###########################################*/

PUBLIC int *
get_iindx(unsigned int length){

  return vrna_idx_row_wise(length);
}

PUBLIC int *
get_indx(unsigned int length){

  return vrna_idx_col_wise(length);
}

PUBLIC void
print_tty_input_seq(void){

  vrna_message_input_seq_simple();
}

PUBLIC void
print_tty_input_seq_str(const char *s){

  vrna_message_input_seq(s);
}

PUBLIC void
warn_user(const char message[]){

  vrna_message_warning(message);
}

PUBLIC void
nrerror(const char message[]){

  vrna_message_error(message);
}

PUBLIC void *space(unsigned size) {

  return vrna_alloc(size);
}

#undef  xrealloc
/* dmalloc.h #define's vrna_realloc */
PUBLIC void *xrealloc(void *p, unsigned size){

  return vrna_realloc(p, size);
}

PUBLIC void
init_rand(void){

  vrna_init_rand();
}

PUBLIC double urn(void){

  return vrna_urn();
}

PUBLIC int
int_urn(int from, int to){

  return vrna_int_urn(from, to);
}

PUBLIC void
filecopy(FILE *from, FILE *to){

  vrna_file_copy(from, to);
}

PUBLIC char *
time_stamp(void){

  return vrna_time_stamp();
}

PUBLIC char *
get_line(FILE *fp){ /* reads lines of arbitrary length from fp */

  return vrna_read_line(fp);
}
#endif
