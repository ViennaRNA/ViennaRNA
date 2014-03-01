#ifndef __VIENNA_RNA_PACKAGE_ALN_UTIL_H__
#define __VIENNA_RNA_PACKAGE_ALN_UTIL_H__

int read_clustal( FILE *clust,
                  char *AlignedSeqs[],
                  char *names[]);
/*@only@*/ /*@notnull@*/ char *consensus(const char *AS[]);
/*@only@*/ /*@notnull@*/ char *consens_mis(const char *AS[]);

char *
get_ungapped_sequence(const char *seq);

#endif
