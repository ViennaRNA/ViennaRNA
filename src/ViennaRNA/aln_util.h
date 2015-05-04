#ifndef VIENNA_RNA_PACKAGE_ALN_UTIL_H
#define VIENNA_RNA_PACKAGE_ALN_UTIL_H

#ifdef __GNUC__
#define DEPRECATED(func) func __attribute__ ((deprecated))
#else
#define DEPRECATED(func) func
#endif

int read_clustal( FILE *clust,
                  char *AlignedSeqs[],
                  char *names[]);
/*@only@*/ /*@notnull@*/ char *consensus(const char *AS[]);
/*@only@*/ /*@notnull@*/ char *consens_mis(const char *AS[]);

char *
get_ungapped_sequence(const char *seq);

/**
 *  @brief Get the mean pairwise identity in steps from ?to?(ident)
 * 
 *  @ingroup consensus_fold
 * 
 *  @param Alseq
 *  @param n_seq  The number of sequences in the alignment
 *  @param length The length of the alignment
 *  @param mini
 *  @return       The mean pairwise identity
 */
int vrna_ali_get_mpi( char *Alseq[],
                      int n_seq,
                      int length,
                      int *mini);

/**
 *  @brief Get the mean pairwise identity in steps from ?to?(ident)
 * 
 *  @ingroup consensus_fold
 *
 *  @deprecated Use vrna_ali_get_mpi() as a replacement
 *
 *  @param Alseq
 *  @param n_seq  The number of sequences in the alignment
 *  @param length The length of the alignment
 *  @param mini
 *  @return       The mean pairwise identity
 */
DEPRECATED(int get_mpi(char *Alseq[], int n_seq, int length, int *mini));

/*
#############################################################
# some helper functions that might be useful in the library #
#############################################################
*/

/**
 *  @brief Get arrays with encoded sequence of the alignment
 *
 *  this function assumes that in S, S5, s3, ss and as enough
 *  space is already allocated (size must be at least sequence length+2)
 * 
 *  @ingroup consensus_fold
 * 
 *  @param sequence The gapped sequence from the alignment
 *  @param S        pointer to an array that holds encoded sequence
 *  @param s5      pointer to an array that holds the next base 5' of alignment position i
 *  @param s3      pointer to an array that holds the next base 3' of alignment position i
 *  @param ss
 *  @param as
 *  @param circ    assume the molecules to be circular instead of linear (circ=0)
 */
void encode_ali_sequence( const char *sequence,
                          short *S,
                          short *s5,
                          short *s3,
                          char *ss,
                          unsigned short *as,
                          int circ);

/**
 *  @brief Allocate memory for sequence array used to deal with aligned sequences
 * 
 *  Note that these arrays will also be initialized according to the sequence alignment given
 * 
 *  @ingroup consensus_fold
 * 
 *  @see free_sequence_arrays()
 * 
 *  @param sequences  The aligned sequences
 *  @param S          A pointer to the array of encoded sequences
 *  @param S5         A pointer to the array that contains the next 5' nucleotide of a sequence position
 *  @param S3         A pointer to the array that contains the next 3' nucleotide of a sequence position
 *  @param a2s        A pointer to the array that contains the alignment to sequence position mapping
 *  @param Ss         A pointer to the array that contains the ungapped sequence
 *  @param circ       assume the molecules to be circular instead of linear (circ=0)
 */
void  alloc_sequence_arrays(const char **sequences,
                            short ***S,
                            short ***S5,
                            short ***S3,
                            unsigned short ***a2s,
                            char ***Ss,
                            int circ);

/**
 *  @brief Free the memory of the sequence arrays used to deal with aligned sequences
 * 
 *  This function frees the memory previously allocated with alloc_sequence_arrays()
 * 
 *  @ingroup consensus_fold
 * 
 *  @see alloc_sequence_arrays()
 * 
 *  @param n_seq      The number of aligned sequences
 *  @param S          A pointer to the array of encoded sequences
 *  @param S5         A pointer to the array that contains the next 5' nucleotide of a sequence position
 *  @param S3         A pointer to the array that contains the next 3' nucleotide of a sequence position
 *  @param a2s        A pointer to the array that contains the alignment to sequence position mapping
 *  @param Ss         A pointer to the array that contains the ungapped sequence
 */
void  free_sequence_arrays( unsigned int n_seq,
                            short ***S,
                            short ***S5,
                            short ***S3,
                            unsigned short ***a2s,
                            char ***Ss);



#endif
