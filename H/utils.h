/* Header file for utils.c */

extern void  *space(unsigned size);           /* allocate space safely */
extern void   nrerror(const char message[]);  /* die with error message */
extern void   init_rand(void);                /* make random number seeds */
extern unsigned short xsubi[3];               /* 48bit random number */
extern double urn(void);                      /* random number from [0..1] */
extern int    int_urn(int from, int to);      /* random integer */
extern void   filecopy(FILE *from, FILE *to); 
extern char  *time_stamp(void);               /* current date in a string */
extern char  *random_string(int l, const char symbols[]);  
extern int    hamming(const char s1[], const char s2[]);  /* hamming distance */
extern char *get_line(FILE *fp);              /* read one line */

