/* Header file for utils.c */

#ifdef DMALLOC
/* use dmalloc library to check for memory management bugs */
#include "/usr/local/include/dmalloc.h"
#define space(S) calloc(1,(S))
#else
/*@only@*/ /*@notnull@*/
extern void  *space(unsigned size);           /* allocate space safely */
#endif
/*@exits@*/
extern void   nrerror(const char message[]);  /* die with error message */
extern void   init_rand(void);                /* make random number seeds */
extern unsigned short xsubi[3];               /* current 48bit random number */
extern double urn(void);                      /* random number from [0..1] */
extern int    int_urn(int from, int to);      /* random integer */
extern void   filecopy(FILE *from, FILE *to); /* inefficient `cp' */
/*@observer@*/
extern char  *time_stamp(void);               /* current date in a string */
extern char  *random_string(int l, const char symbols[]);
/* random string of length l using characters from symbols[] */
extern int    hamming(const char *s1, const char *s2);
/* calculate hamming distance */
extern char  *get_line(FILE *fp); /* read one (arbitrary length) line from fp */


extern char *pack_structure(const char *struc);
/* pack secondary secondary structure, 5:1 compression using base 3 encoding */
extern char *unpack_structure(const char *packed);
/* unpack sec structure packed with pack_structure() */
extern short *make_pair_table(const char *structure);
/* returns a newly allocated table, such that:  table[i]=j if (i.j) pair or
   0 if i is unpaired, table[0] contains the length of the structure. */

extern int bp_distance(const char *str1, const char *str2);
/* dist = {number of base pairs in one structure but not in the other} 
   same as edit distance with open-pair close-pair as move-set */
