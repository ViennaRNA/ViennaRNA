#ifndef _UTILS_H_
#define _UTILS_H_

/* Header file for utils.c */

#if defined(__cplusplus)
extern "C" {
#endif

    int   naview_xy_coordinates(short *pair_table, float *X, float *Y);
    void  *space(unsigned size);           /* allocate space safely */

    /*@exits@*/
    void   nrerror(const char message[]);  /* die with error message */
    void   init_rand(void);                /* make random number seeds */
    //unsigned short xsubi[3];               /* current 48bit random number */
    double urn(void);                      /* random number from [0..1] */
    int    int_urn(int from, int to);      /* random integer */
    void   filecopy(FILE *from, FILE *to); /* inefficient `cp' */
    /*@observer@*/
    char  *time_stamp(void);               /* current date in a string */
    char  *randostring_(int l, const char symbols[]);
    /* random string of length l using characters from symbols[] */
    int    hamming(const char *s1, const char *s2);
    /* calculate hamming distance */
    char  *get_line(FILE *fp); /* read one (arbitrary length) line from fp */


    char *pack_structure(const char *struc);
    /* pack secondary secondary structure, 5:1 compression using base 3 encoding */
    char *unpack_structure(const char *packed);
    /* unpack sec structure packed with pack_structure() */
    short *make_pair_table(const char *structure);
    /* returns a newly allocated table, such that:  table[i]=j if (i.j) pair or
       0 if i is unpaired, table[0] contains the length of the structure. */

    int bp_distance(const char *str1, const char *str2);
    /* dist = {number of base pairs in one structure but not in the other}
       same as edit distance with open-pair close-pair as move-set */

#if defined(__cplusplus)
}
#endif

#endif
