#ifndef ALN_UTIL_H
#define ALN_UTIL_H

extern int read_clustal(FILE *clust, char *AlignedSeqs[], char *names[]);
extern /*@only@*/ /*@notnull@*/ char *consensus(const char *AS[]);
extern /*@only@*/ /*@notnull@*/ char *consens_mis(const char *AS[]);

#endif
