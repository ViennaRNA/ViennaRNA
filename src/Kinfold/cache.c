/*
  Last changed Time-stamp: <2006-10-03 16:31:54 xtof>
  c  Christoph Flamm and Ivo L Hofacker
  {xtof,ivo}@tbi.univie.ac.at
  Kinfold: $Name:  $
  $Id: cache.c,v 1.3 2006/10/04 12:45:12 xtof Exp $
*/

#include "config.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#if HAVE_LIBRNA_API3
#include <ViennaRNA/utils.h>
#else
#include <utils.h>
#endif

#include "cache_util.h"

#ifdef __GNUC__
# define INLINE inline
#else
# define INLINE
#endif

/*
  modify cache_f(), cache_comp() and the typedef of cache_entry
  in cache_utils.h to suit your application
*/

/* PUBLIC FUNCTIONES */
cache_entry *lookup_cache (char *x);
int write_cache (cache_entry *x);
/*  void delete_cache (cache_entry *x); */
void kill_cache();
void initialize_cache();

/* PRIVATE FUNCTIONES */
/*  static int cache_comp(cache_entry *x, cache_entry *y); */
INLINE static unsigned cache_f (char *x);

/* #define CACHESIZE 67108864 -1 */ /* 2^26 -1   must be power of 2 -1 */
/* #define CACHESIZE 33554432 -1 */ /* 2^25 -1   must be power of 2 -1 */
/* #define CACHESIZE 16777216 -1 */ /* 2^24 -1   must be power of 2 -1 */ 
/* #define CACHESIZE  4194304 -1 */ /* 2^22 -1   must be power of 2 -1 */
#define CACHESIZE  1048576 -1  /* 2^20 -1   must be power of 2 -1 */
/* #define CACHESIZE   262144 -1 */ /* 2^18 -1   must be power of 2 -1 */
/* next is default */
/* #define CACHESIZE    65536 -1 */ /* 2^16 -1   must be power of 2 -1 */
/* #define CACHESIZE    16384 -1 */ /* 2^14 -1   must be power of 2 -1 */
/* #define CACHESIZE     4096 -1 */ /* 2^12 -1   must be power of 2 -1 */

static cache_entry *cachetab[CACHESIZE+1];
static char UNUSED rcsid[] ="$Id: cache.c,v 1.3 2006/10/04 12:45:12 xtof Exp $";
unsigned long collisions=0;

/* stolen from perl source */
char coeff[] = {
                61,59,53,47,43,41,37,31,29,23,17,13,11,7,3,1,
                61,59,53,47,43,41,37,31,29,23,17,13,11,7,3,1,
                61,59,53,47,43,41,37,31,29,23,17,13,11,7,3,1,
                61,59,53,47,43,41,37,31,29,23,17,13,11,7,3,1,
                61,59,53,47,43,41,37,31,29,23,17,13,11,7,3,1,
                61,59,53,47,43,41,37,31,29,23,17,13,11,7,3,1,
                61,59,53,47,43,41,37,31,29,23,17,13,11,7,3,1,
                61,59,53,47,43,41,37,31,29,23,17,13,11,7,3,1};

/* key must not be longer than 128 */
INLINE static unsigned cache_f(char *x) { 
  register char *s;
  register int i;
  register unsigned cache;

  s = x;

  for (i=0,    cache = 0;
       /* while */ *s;
       s++,           i++ , cache *= 5 ) {
    cache += *s * coeff[i];
  }

   /* divide through CACHESIZE for normalization */
  return ((cache) & (CACHESIZE));
}


/* returns NULL unless x is in the cache */
cache_entry *lookup_cache (char *x) {
  int cacheval;
  cache_entry *c;

  cacheval=cache_f(x);
  if ((c=cachetab[cacheval]))
    if (strcmp(c->structure,x)==0) return c;
  
  return NULL;
}

/* returns 1 if x already was in the cache */
int write_cache (cache_entry *x) {
  int cacheval;
  cache_entry *c;
  
  cacheval=cache_f(x->structure);
  if ((c=cachetab[cacheval])) {
    free(c->structure);
    free(c->neighbors);
    free(c->rates);
    free(c->energies);
    free(c);
  }
  cachetab[cacheval]=x;
  return 0;
}

/**/
void initialize_cache () { }

/**/
void kill_cache () {
  int i;
  
  for (i=0;i<CACHESIZE+1;i++) {
    if ( cachetab[i] ) {
      free (cachetab[i]->structure);
      free (cachetab[i]->neighbors);
      free (cachetab[i]->rates);
      free (cachetab[i]->energies);
      free (cachetab[i]);
    }
    cachetab[i]=NULL;
  }
}

#if 0
/**/
static int cache_comp(cache_entry *x, cache_entry *y) {
  return strcmp(((cache_entry *)x)->structure, ((cache_entry *)y)->structure);
}
#endif

/* End of file */
