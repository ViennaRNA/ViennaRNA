/*
                               utils.c

                 c  Ivo L Hofacker and Walter Fontana
                          Vienna RNA package
*/
/* Last changed Time-stamp: <2008-11-25 16:34:36 ivo> */

#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <time.h>
#include <string.h>
#include "../config.h"
#include "utils.h"

#ifdef WITH_DMALLOC
#include "dmalloc.h"
#endif
/*@unused@*/
static char rcsid[] = "$Id: utils.c,v 1.19 2008/12/16 22:30:30 ivo Exp $";

#define PRIVATE  static
#define PUBLIC

/*@notnull@ @only@*/
PUBLIC unsigned short xsubi[3];

static char  scale1[] = "....,....1....,....2....,....3....,....4";
static char  scale2[] = "....,....5....,....6....,....7....,....8";

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
  fprintf(stderr, "\n%s\n", message);
  exit(EXIT_FAILURE);
}

PUBLIC void warn_user(const char message[]){
  fprintf(stderr, "\nWARNING: %s\n", message);
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
      if(!(option & VRNA_INPUT_NOPRINT_COMMENTS)) printf("%s\n", line);
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
  //if(!(option & VRNA_INPUT_NOPRINT)) printf("%s\n", line);

  /* eliminate whitespaces at the end of the line read */
  if(!(option & VRNA_INPUT_NOELIM_WS_SUFFIX)){
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

PUBLIC short *copy_pair_table(const short *pt){
  short *table = (short *)space(sizeof(short) * (pt[0]+2));
  memcpy(table, pt, sizeof(short)*(pt[0]+2));
  return table;
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
}

PUBLIC  void  print_tty_constraint_full(void){
  print_tty_constraint(VRNA_CONSTRAINT_PIPE | VRNA_CONSTRAINT_DOT | VRNA_CONSTRAINT_X | VRNA_CONSTRAINT_ANG_BRACK | VRNA_CONSTRAINT_RND_BRACK);
}

PUBLIC  void  print_tty_constraint(unsigned int option){
  if(!(option & VRNA_CONSTRAINT_NO_HEADER)) printf("Input structure constraints using the following notation:\n");
  if(option & VRNA_CONSTRAINT_PIPE)       printf("| : paired with another base\n");
  if(option & VRNA_CONSTRAINT_DOT)        printf(". : no constraint at all\n");
  if(option & VRNA_CONSTRAINT_X)          printf("x : base must not pair\n");
  if(option & VRNA_CONSTRAINT_ANG_BRACK)  printf("< : base i is paired with a base j<i\n> : base i is paired with a base j>i\n");
  if(option & VRNA_CONSTRAINT_RND_BRACK)  printf("matching brackets ( ): base i pairs base j\n");
}

PUBLIC  void  str_DNA2RNA(char *sequence){
  unsigned int l, i;
  if(sequence != NULL){
    l = strlen(sequence);
    for(i = 0; i < l; i++){
      sequence[i] = toupper(sequence[i]);
      if(sequence[i] == 'T') sequence[i] = 'U';
    }
  }
}

PUBLIC  void  str_RNA2RNA(char *sequence){
  unsigned int l, i;
  if(sequence != NULL){
    l = strlen(sequence);
    for(i = 0; i < l; i++)
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

PUBLIC void constrain_ptypes(const char *constraint, char *ptype, int *BP, int min_loop_size, unsigned int idx_type){
  int n,i,j,k,l;
  int hx, *stack;
  char type;
  int *index;

  if(constraint == NULL) return;

  n = (int)strlen(constraint);

  stack = (int *) space(sizeof(int)*(n+1));

  if(!idx_type){ /* index allows access in energy matrices at pos (i,j) via index[j]+i */
    index = get_indx((unsigned)n);

    for(hx=0, j=1; j<=n; j++){
      switch(constraint[j-1]){
        case '|':   if(BP) BP[j] = -1;
                    break;
        case 'x':   /* can't pair */
                    for (l=1; l<j-min_loop_size; l++) ptype[index[j]+l] = 0;
                    for (l=j+min_loop_size+1; l<=n; l++) ptype[index[l]+j] = 0;
                    break;
        case '(':   stack[hx++]=j;
                    /* fallthrough */
        case '<':   /* pairs upstream */
                    for (l=1; l<j-min_loop_size; l++) ptype[index[j]+l] = 0;
                    break;
        case ')':   if (hx<=0) {
                      fprintf(stderr, "%s\n", constraint);
                      nrerror("unbalanced brackets in constraint");
                    }
                    i = stack[--hx];
                    type = ptype[index[j]+i];
                    for (k=i+1; k<=n; k++) ptype[index[k]+i] = 0;
                    /* don't allow pairs i<k<j<l */
                    for (l=j; l<=n; l++)
                      for (k=i+1; k<=j; k++) ptype[index[l]+k] = 0;
                    /* don't allow pairs k<i<l<j */
                    for (l=i; l<=j; l++)
                      for (k=1; k<=i; k++) ptype[index[l]+k] = 0;
                    for (k=1; k<j; k++) ptype[index[j]+k] = 0;
                    ptype[index[j]+i] = (type==0) ? 7 : type;
                    /* fallthrough */
        case '>':   /* pairs downstream */
                    for (l=j+min_loop_size+1; l<=n; l++) ptype[index[l]+j] = 0;
                    break;
      }
    }
  }
  else{ /* index allows access in energy matrices at pos (i,j) via index[i]-j */
    index = get_iindx((unsigned)n);

    for(hx=0, j=1; j<=n; j++) {
      switch (constraint[j-1]) {
        case 'x':   /* can't pair */
                    for (l=1; l<j-min_loop_size; l++) ptype[index[l]-j] = 0;
                    for (l=j+min_loop_size+1; l<=n; l++) ptype[index[j]-l] = 0;
                    break;
        case '(':   stack[hx++]=j;
                    /* fallthrough */
        case '<':   /* pairs upstream */
                    for (l=1; l<j-min_loop_size; l++) ptype[index[l]-j] = 0;
                    break;
        case ')':   if (hx<=0) {
                      fprintf(stderr, "%s\n", constraint);
                      nrerror("unbalanced brackets in constraints");
                    }
                    i = stack[--hx];
                    type = ptype[index[i]-j];
                    /* don't allow pairs i<k<j<l */
                    for (k=i; k<=j; k++)
                      for (l=j; l<=n; l++) ptype[index[k]-l] = 0;
                    /* don't allow pairs k<i<l<j */
                    for (k=1; k<=i; k++)
                      for (l=i; l<=j; l++) ptype[index[k]-l] = 0;
                    ptype[index[i]-j] = (type==0) ? 7 : type;
                    /* fallthrough */
        case '>':   /* pairs downstream */
                    for (l=j+min_loop_size+1; l<=n; l++) ptype[index[j]-l] = 0;
                    break;
      }
    }
  }
  if (hx!=0) {
    fprintf(stderr, "%s\n", constraint);
    nrerror("unbalanced brackets in constraint string");
  }
  free(index);
  free(stack);
}

/* get a matrix containing the number of basepairs of a reference structure for each interval [i,j] with i<j
*  access it via iindx!!!
*/
PUBLIC unsigned int *make_referenceBP_array(short *reference_pt, unsigned int turn){
  unsigned int i,j,k,ij,length, *iindx;
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
  unsigned int *iindx = get_iindx(n);
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
