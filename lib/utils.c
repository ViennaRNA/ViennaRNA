/*
			       utils.c

		 c  Ivo L Hofacker and Walter Fontana
			  Vienna RNA package
*/
/* Last changed Time-stamp: <1998-07-18 20:53:29 ivo> */

#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <time.h>
#include <string.h>

static char rcsid[] = "$Id: utils.c,v 1.6 1998/07/19 14:31:15 ivo Exp $";

#define PRIVATE  static
#define PUBLIC

PUBLIC void  *space(unsigned int size);
PUBLIC void   nrerror(const char message[]);
PUBLIC double urn(void);
PUBLIC int    int_urn(int from, int to);
PUBLIC void   filecopy(FILE *from, FILE *to);
PUBLIC char  *time_stamp(void);
PUBLIC char  *random_string(int l, const char symbols[]);
PUBLIC int    hamming(const char *s1, const char *s2);
PUBLIC char  *get_line(FILE *fp);

PUBLIC unsigned short xsubi[3];

/*-------------------------------------------------------------------------*/

PUBLIC void *space(unsigned size)
{
    void *pointer;
    
    if ( (pointer = (void *) calloc(1, size)) == NULL) {
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

/*------------------------------------------------------------------------*/

PUBLIC void nrerror(const char message[])       /* output message upon error */

               

{
    fprintf(stderr, "\n%s\n", message);
    exit(-1);
}

/*------------------------------------------------------------------------*/
PUBLIC void init_rand(void)
{
   time_t t;
   time(&t);
   xsubi[0] = (unsigned short) t;
   xsubi[1] = (unsigned short) (t >> 16);
   xsubi[2] = 5246;
}

/*------------------------------------------------------------------------*/
 
PUBLIC double urn(void)    
		/* uniform random number generator; urn() is in [0,1] */
                /* uses a linear congruential library routine */ 
                /* 48 bit arithmetic */
{
    extern double erand48(unsigned short[3]);

    return erand48(xsubi);
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
    
    while ((c = getc(from)) != EOF) putc(c, to);
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
  
   base = strlen(symbols);
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

   line = NULL;
   do {
      if (fgets(s, 512, fp)==NULL) break;
      cp = strchr(s, '\n');
      if (cp != NULL) *cp = '\0';
      if (line==NULL)
	 line = space(strlen(s)+1);
      else
	 line = (char *) realloc(line, strlen(s)+strlen(line)+1);
      strcat(line,s);
   } while(cp==NULL);

   return line;
}


/*-----------------------------------------------------------------*/

PUBLIC char *pack_structure(const char *struc) {
  /* 5:1 compression using base 3 encoding */
  int i,j,l,pi;
  unsigned char *packed;

  l = strlen(struc);
  packed = (unsigned char *) space(((l+4)/5+1)*sizeof(unsigned char));

  j=i=pi=0; 
  while (i<l) {
    register unsigned char p;
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
    packed[j++] = p+1;   /* never use 0, so we can use strcmp() etc. */
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

  l = strlen(packed);
  pp = (unsigned char *) packed;
  struc = (char *) space((l*5+1)*sizeof(char));   /* up to 4 byte extra */

  j=0;
  for (i=j=0; i<l; i++) {
    unsigned char p;
    int k, c;
    
    p = pp[i]-1;
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

				   
/*---------------------------------------------------------------------------*/ 

PUBLIC short *make_pair_table(const char *structure)
{
    /* returns array representation of structure.
       table[i] is 0 if unpaired or j if (i.j) pair.  */
   int i,j,hx;
   int length;
   short *stack;
   short *table;
   
   length = strlen(structure);
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

/*---------------------------------------------------------------------------*/

PUBLIC int bp_distance(const char *str1, const char *str2)
{
  /* dist = {number of base pairs in one structure but not in the other} */
  /* same as edit distance with pair_open pair_close as move set */
   int dist,i,l;
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
