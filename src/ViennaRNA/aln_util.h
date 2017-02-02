#ifndef VIENNA_RNA_PACKAGE_ALN_UTIL_H
#define VIENNA_RNA_PACKAGE_ALN_UTIL_H

#ifdef VRNA_WARN_DEPRECATED
# ifdef __GNUC__
#  define DEPRECATED(func) func __attribute__ ((deprecated))
# else
#  define DEPRECATED(func) func
# endif
#else
# define DEPRECATED(func) func
#endif

/**
 *  @file aln_util.h
 *  @ingroup utils
 *  @brief Various utility- and helper-functions for sequence alignments and comparative structure prediction
 */

/**
 *  @{
 *  @ingroup   aln_utils
 *
 */

/** @brief Typename for the base pair info repesenting data structure #vrna_pinfo_s */
typedef struct vrna_pinfo_s     vrna_pinfo_t;

/* make this interface backward compatible with RNAlib < 2.2.0 */
#define VRNA_BACKWARD_COMPAT


#ifdef VRNA_BACKWARD_COMPAT

/* the following typedefs are for backward compatibility only */

/**
 *  @brief Old typename of #vrna_pinfo_s
 *  @deprecated Use #vrna_pinfo_t instead!
*/
typedef struct vrna_pinfo_s     pair_info;

#endif

/**
 *  @brief A base pair info structure
 *
 *  For each base pair (i,j) with i,j in [0, n-1] the structure lists:
 *  - its probability 'p'
 *  - an entropy-like measure for its well-definedness 'ent'
 *  - the frequency of each type of pair in 'bp[]'
 *    + 'bp[0]' contains the number of non-compatible sequences
 *    + 'bp[1]' the number of CG pairs, etc.
 */
struct vrna_pinfo_s {
   unsigned i;    /**<  @brief  nucleotide position i */
   unsigned j;    /**<  @brief  nucleotide position j */
   float p;       /**< @brief  Probability */
   float ent;     /**< @brief  Pseudo entropy for @f$ p(i,j) = S_i + S_j - p_ij*ln(p_ij) @f$ */
   short bp[8];   /**< @brief  Frequencies of pair_types */
   char comp;     /**< @brief  1 iff pair is in mfe structure */
};

int read_clustal( FILE *clust,
                  char *AlignedSeqs[],
                  char *names[]);

char *consensus(const char *AS[]);

char *consens_mis(const char *AS[]);

char *
get_ungapped_sequence(const char *seq);

/**
 *  @brief Get the mean pairwise identity in steps from ?to?(ident)
 * 
 *  @ingroup consensus_fold
 * 
 *  @param alignment  Aligned sequences
 *  @return       The mean pairwise identity
 */
int vrna_aln_mpi( const char **alignment);

/**
 *  \brief Retrieve an array of #vrna_pinfo_t structures from precomputed pair probabilities
 *
 *  This array of structures contains information about positionwise pair probabilies,
 *  base pair entropy and more
 *
 *  \see #vrna_pinfo_t, and vrna_pf()
 *
 *  \param  vc          The #vrna_fold_compound_t of type #VRNA_FC_TYPE_COMPARATIVE with precomputed partition function matrices
 *  \param  structure   An optional structure in dot-bracket notation (Maybe NULL)
 *  \param  threshold   Do not include results with pair probabilities below threshold
 *  \return             The #vrna_pinfo_t array
 */
vrna_pinfo_t *vrna_aln_pinfo(vrna_fold_compound_t *vc,
                                  const char *structure,
                                  double threshold);

int *
vrna_aln_pscore(const char  **alignment,
                vrna_md_t   *md);


/**
 *  @brief Get the mean pairwise identity in steps from ?to?(ident)
 * 
 *  @ingroup consensus_fold
 *
 *  @deprecated Use vrna_aln_mpi() as a replacement
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


/**
 * @}
 */


#endif
