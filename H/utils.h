/* Header file for utils.c */

extern void  *space(unsigned size);
extern void   nrerror(const char message[]);
extern void   init_rand(void);
extern double urn(void);
extern int    int_urn(int from, int to);
extern void   filecopy(FILE *from, FILE *to);
extern char  *time_stamp(void);
extern char  *random_string(int l, const char symbols[]);
extern int    hamming(const char s1[], const char s2[]);
extern char *get_line(FILE *fp);
extern unsigned short xsubi[3];

