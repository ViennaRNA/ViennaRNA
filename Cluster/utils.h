/* Header file for utils.c */

extern void  *space(unsigned size);
extern void   nrerror(char *message);
extern double urn();
extern int    int_urn(int from, int to);
extern void   filecopy(FILE *from, FILE *to);
extern char  *time_stamp();
extern char  *random_string(int l, char *symbols);
extern int    hamming(char *s1, char *s2);
extern char *get_line(FILE *fp);
extern unsigned short xsubi[3];

