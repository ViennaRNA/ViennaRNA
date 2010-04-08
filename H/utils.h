#ifndef __VIENNA_RNA_PACKAGE_UTILS_H__
#define __VIENNA_RNA_PACKAGE_UTILS_H__

/** \file
***
**/

#ifdef HAVE_CONFIG_H
#include <config.h>
#ifndef HAVE_STRDUP
char *strdup(const char *s);
#endif
#endif
#ifdef WITH_DMALLOC
/* use dmalloc library to check for memory management bugs */
#include "dmalloc.h"
#define space(S) calloc(1,(S))
#else
/*@only@*/ /*@notnull@*/
void  *space(unsigned size) /*@ensures MaxSet(result) == (size-1);@*/;
			    /* allocate space safely */
/*@only@*/ /*@notnull@*/
void  *xrealloc(/*@null@*/ /*@only@*/ /*@out@*/ /*@returned@*/ void *p, unsigned size) /*@modifies *p @*/ /*@ensures MaxSet(result) == (size-1) @*/;
#endif

/*@exits@*/ void nrerror(const char message[]);  /* die with error message */
void   init_rand(void);                /* make random number seeds */
extern unsigned short xsubi[3];               /* current 48bit random number */
double urn(void);                      /* random number from [0..1] */
int    int_urn(int from, int to);      /* random integer */
void   filecopy(FILE *from, FILE *to); /* inefficient `cp' */
/*@observer@*/ char  *time_stamp(void);               /* current date in a string */
/*@only@*/ /*@notnull@*/ char  *random_string(int l, const char symbols[]);
/* random string of length l using characters from symbols[] */
/**
*** calculate hamming distance
**/
int    hamming(const char *s1, const char *s2);
/*@only@*/ /*@null@*/ char  *get_line(const FILE *fp); /* read one (arbitrary length) line from fp */


/**
*** pack secondary secondary structure, 5:1 compression using base 3 encoding
**/
char *pack_structure(const char *struc);
/**
*** unpack sec structure packed with pack_structure()
**/
char *unpack_structure(const char *packed);
/**
*** returns a newly allocated table, such that:  table[i]=j if (i.j) pair or
*** 0 if i is unpaired, table[0] contains the length of the structure.
**/
short *make_pair_table(const char *structure);
/**
*** dist = {number of base pairs in one structure but not in the other}
*** same as edit distance with open-pair close-pair as move-set
**/
int bp_distance(const char *str1, const char *str2);



#endif
