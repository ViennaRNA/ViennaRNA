/*
			       utils.c

		 c  Ivo L Hofacker and Walter Fontana
			  Vienna RNA package
*/
/* Last changed Time-stamp: <97/10/27 14:20:34 ivo> */

#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <errno.h>
#include <time.h>
#include <string.h>
#ifdef DBMALLOC
#include "/usr/local/debug_include/malloc.h"
#endif

static char rcsid[] = "$Id: utils.c,v 1.2 1997/10/27 13:21:01 ivo Rel $";

#define PRIVATE  static
#define PUBLIC

PUBLIC void  *space(unsigned int size);
PUBLIC void   nrerror(const char message[]);
PUBLIC double urn(void);
PUBLIC int    int_urn(int from, int to);
PUBLIC void   filecopy(FILE *from, FILE *to);
PUBLIC char  *time_stamp(void);
PUBLIC char  *random_string(int l, const char symbols[]);
PUBLIC int    hamming(const char s1[], const char s2[]);
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

PUBLIC int   hamming(const char s1[], const char s2[])
{
   int h=0,i;
   
   for (i=0; i<strlen(s1); i++)
      if (s1[i]!=s2[i]) h++;
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



