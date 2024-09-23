#ifndef VIENNA_RNA_PACKAGE_SEQUENCES_MSA_H
#define VIENNA_RNA_PACKAGE_SEQUENCES_MSA_H

#ifdef VRNA_WARN_DEPRECATED
# if defined(__clang__)
#  define DEPRECATED(func, msg) func __attribute__ ((deprecated("", msg)))
# elif defined(__GNUC__)
#  define DEPRECATED(func, msg) func __attribute__ ((deprecated(msg)))
# else
#  define DEPRECATED(func, msg) func
# endif
#else
# define DEPRECATED(func, msg) func
#endif

/**
 *  @file ViennaRNA/utils/alignments.h
 *  @ingroup utils, aln_utils
 *  @brief Various utility- and helper-functions for sequence alignments and comparative structure prediction
 */

/**
 *  @addtogroup aln_utils
 *  @{
 */

/** @brief Typename for the base pair info repesenting data structure #vrna_pinfo_s */
typedef struct vrna_pinfo_s vrna_pinfo_t;


/**
 *  @brief  Use default alignment settings
 */
#define VRNA_ALN_DEFAULT      0U


/**
 *  @brief  Convert to RNA alphabet
 */
#define VRNA_ALN_RNA          1U


/**
 *  @brief  Convert to DNA alphabet
 */
#define VRNA_ALN_DNA          2U


/**
 *  @brief  Convert to uppercase nucleotide letters
 */
#define VRNA_ALN_UPPERCASE    4U


/**
 *  @brief  Convert to lowercase nucleotide letters
 */
#define VRNA_ALN_LOWERCASE    8U

/**
 *  @brief  Flag indicating Shannon Entropy measure
 *
 *  Shannon Entropy is defined as @f$ H = - \sum_c p_c \cdot \log_2 p_c @f$
 */
#define VRNA_MEASURE_SHANNON_ENTROPY  1U

#ifndef VRNA_DISABLE_BACKWARD_COMPATIBILITY

/* the following typedefs are for backward compatibility only */

/**
 *  @brief Old typename of #vrna_pinfo_s
 *  @deprecated Use #vrna_pinfo_t instead!
 *  @ingroup  aln_utils_deprecated
 */
typedef struct vrna_pinfo_s pair_info;

#endif

#include <ViennaRNA/fold_compound.h>
#include <ViennaRNA/model.h>

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
  unsigned  i;      /**<  @brief  nucleotide position i */
  unsigned  j;      /**<  @brief  nucleotide position j */
  float     p;      /**< @brief  Probability */
  float     ent;    /**< @brief  Pseudo entropy for @f$ p(i,j) = S_i + S_j - p_ij*ln(p_ij) @f$ */
  short     bp[8];  /**< @brief  Frequencies of pair_types */
  char      comp;   /**< @brief  1 iff pair is in mfe structure */
};


/**
 *  @brief Get the mean pairwise identity in steps from ?to?(ident)
 *
 *  @param alignment  Aligned sequences
 *  @return       The mean pairwise identity
 */
int
vrna_aln_mpi(const char **alignment);


/**
 *  \brief Retrieve an array of #vrna_pinfo_t structures from precomputed pair probabilities
 *
 *  This array of structures contains information about positionwise pair probabilies,
 *  base pair entropy and more
 *
 *  \see #vrna_pinfo_t, and vrna_pf()
 *
 *  \param  fc          The #vrna_fold_compound_t of type #VRNA_FC_TYPE_COMPARATIVE with precomputed partition function matrices
 *  \param  structure   An optional structure in dot-bracket notation (Maybe NULL)
 *  \param  threshold   Do not include results with pair probabilities below threshold
 *  \return             The #vrna_pinfo_t array
 */
vrna_pinfo_t *
vrna_aln_pinfo(vrna_fold_compound_t *fc,
               const char           *structure,
               double               threshold);


int *
vrna_aln_pscore(const char  **alignment,
                vrna_md_t   *md);


int
vrna_pscore(vrna_fold_compound_t  *fc,
            unsigned int          i,
            unsigned int          j);


int
vrna_pscore_freq(vrna_fold_compound_t *fc,
                 const unsigned int   *frequencies,
                 unsigned int         pairs);

/**
 *  @brief  Slice out a subalignment from a larger alignment
 *
 *  @note   The user is responsible to free the memory occupied by
 *          the returned subalignment
 *
 *  @see    vrna_aln_free()
 *
 *  @param  alignment   The input alignment
 *  @param  i           The first column of the subalignment (1-based)
 *  @param  j           The last column of the subalignment (1-based)
 *  @return             The subalignment between column @f$i@f$ and @f$j@f$
 */
char **
vrna_aln_slice(const char   **alignment,
               unsigned int i,
               unsigned int j);


/**
 *  @brief  Free memory occupied by a set of aligned sequences
 *
 *  @param  alignment   The input alignment
 */
void
vrna_aln_free(char **alignment);


/**
 *  @brief  Create a copy of an alignment with only uppercase letters in the sequences
 *
 *  @see  vrna_aln_copy
 *
 *  @param  alignment   The input sequence alignment (last entry must be @em NULL terminated)
 *  @return             A copy of the input alignment where lowercase sequence letters are replaced by uppercase letters
 */
char **
vrna_aln_uppercase(const char **alignment);


/**
 *  @brief  Create a copy of an alignment where DNA alphabet is replaced by RNA alphabet
 *
 *  @see  vrna_aln_copy
 *
 *  @param  alignment   The input sequence alignment (last entry must be @em NULL terminated)
 *  @return             A copy of the input alignment where DNA alphabet is replaced by RNA alphabet (T -> U)
 */
char **
vrna_aln_toRNA(const char **alignment);


/**
 *  @brief  Make a copy of a multiple sequence alignment
 *
 *  This function allows one to create a copy of a multiple sequence alignment. The @p options parameter
 *  additionally allows for sequence manipulation, such as converting DNA to RNA alphabet, and conversion
 *  to uppercase letters.
 *
 *  @see  vrna_aln_copy(), #VRNA_ALN_RNA, #VRNA_ALN_UPPERCASE, #VRNA_ALN_DEFAULT
 *
 *  @param  alignment   The input sequence alignment (last entry must be @em NULL terminated)
 *  @param  options     Option flags indicating whether the aligned sequences should be converted
 *  @return             A (manipulated) copy of the input alignment
 */
char **
vrna_aln_copy(const char    **alignment,
              unsigned int  options);


/**
 *  @brief Compute base pair conservation of a consensus structure
 *
 *  This function computes the base pair conservation (fraction of canonical base pairs)
 *  of a consensus structure given a multiple sequence alignment. The base pair types
 *  that are considered canonical may be specified using the #vrna_md_t.pair array.
 *  Passing @em NULL as parameter @p md results in default pairing rules, i.e. canonical
 *  Watson-Crick and GU Wobble pairs.
 *
 *  @param  alignment   The input sequence alignment (last entry must be @em NULL terminated)
 *  @param  structure   The consensus structure in dot-bracket notation
 *  @param  md          Model details that specify compatible base pairs (Maybe @em NULL)
 *  @return             A 1-based vector of base pair conservations
 */
float *
vrna_aln_conservation_struct(const char       **alignment,
                             const char       *structure,
                             const vrna_md_t  *md);


/**
 *  @brief Compute nucleotide conservation in an alignment
 *
 *  This function computes the conservation of nucleotides in alignment columns.
 *  The simples measure is Shannon Entropy and can be selected by passing the
 *  #VRNA_MEASURE_SHANNON_ENTROPY flag in the @p options parameter.
 *
 *  @note Currently, only #VRNA_MEASURE_SHANNON_ENTROPY is supported as
 *        conservation measure.
 *
 *  @see #VRNA_MEASURE_SHANNON_ENTROPY
 *
 *  @param  alignment   The input sequence alignment (last entry must be @em NULL terminated)
 *  @param  md          Model details that specify known nucleotides (Maybe @em NULL)
 *  @param  options     A flag indicating which measure of conservation should be applied
 *  @return             A 1-based vector of column conservations
 */
float *
vrna_aln_conservation_col(const char      **alignment,
                          const vrna_md_t *md_p,
                          unsigned int    options);


/**
 *  @brief  Compute the consensus sequence for a given multiple sequence alignment
 *
 *  @param  alignment   The input sequence alignment (last entry must be @em NULL terminated)
 *  @param  md_p        Model details that specify known nucleotides (Maybe @em NULL)
 *  @return             The consensus sequence of the alignment, i.e. the most frequent nucleotide for each alignment column
 */
char *
vrna_aln_consensus_sequence(const char      **alignment,
                            const vrna_md_t *md_p);

/**
 *  @brief  Compute the Most Informative Sequence (MIS) for a given multiple sequence alignment
 *
 *  The most informative sequence (MIS) @rstinline :cite:p:`freyhult:2005` @endrst displays for each alignment column
 *  the nucleotides with frequency greater than the background frequency, projected into IUPAC
 *  notation. Columns where gaps are over-represented are in lower case.
 *
 *  @param  alignment   The input sequence alignment (last entry must be @em NULL terminated)
 *  @param  md_p        Model details that specify known nucleotides (Maybe @em NULL)
 *  @return             The most informative sequence for the alignment
 */
char *
vrna_aln_consensus_mis(const char       **alignment,
                       const vrna_md_t  *md_p);


#ifndef VRNA_DISABLE_BACKWARD_COMPATIBILITY

#include <stdio.h>
/**
 *  @ingroup  aln_utils_deprecated
 */
DEPRECATED(int read_clustal(FILE  *clust,
                            char  *AlignedSeqs[],
                            char  *names[]),
          "Use vrna_file_msa_read() and vrna_file_msa_read_record() instead");


/**
 *  @ingroup  aln_utils_deprecated
 */
DEPRECATED(char *consensus(const char *AS[]),
          "Use vrna_aln_consensus_sequence() instead!");


/**
 *  @ingroup  aln_utils_deprecated
 */
DEPRECATED(char *consens_mis(const char *AS[]),
          "Use vrna_aln_consensus_mis() instead!");


/**
 *  @ingroup  aln_utils_deprecated
 */
DEPRECATED(char *get_ungapped_sequence(const char *seq),
          "Use vrna_seq_ungapped() instead!");


/**
 *  @brief Get the mean pairwise identity in steps from ?to?(ident)
 *
 *  @deprecated Use vrna_aln_mpi() as a replacement
 *  @ingroup  aln_utils_deprecated
 *  @param Alseq
 *  @param n_seq  The number of sequences in the alignment
 *  @param length The length of the alignment
 *  @param mini
 *  @return       The mean pairwise identity
 */
DEPRECATED(int get_mpi(char *Alseq[],
                       int  n_seq,
                       int  length,
                       int  *mini),
          "Use vrna_aln_mpi() instead");

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
 *  @ingroup  aln_utils_deprecated
 *  @param sequence The gapped sequence from the alignment
 *  @param S        pointer to an array that holds encoded sequence
 *  @param s5      pointer to an array that holds the next base 5' of alignment position i
 *  @param s3      pointer to an array that holds the next base 3' of alignment position i
 *  @param ss
 *  @param as
 *  @param circ    assume the molecules to be circular instead of linear (circ=0)
 */
DEPRECATED(void encode_ali_sequence(const char      *sequence,
                                    short           *S,
                                    short           *s5,
                                    short           *s3,
                                    char            *ss,
                                    unsigned short  *as,
                                    int             circ),
          "This function is obsolete");


/**
 *  @brief Allocate memory for sequence array used to deal with aligned sequences
 *
 *  Note that these arrays will also be initialized according to the sequence alignment given
 *
 *  @see free_sequence_arrays()
 *
 *  @ingroup  aln_utils_deprecated
 *  @param sequences  The aligned sequences
 *  @param S          A pointer to the array of encoded sequences
 *  @param S5         A pointer to the array that contains the next 5' nucleotide of a sequence position
 *  @param S3         A pointer to the array that contains the next 3' nucleotide of a sequence position
 *  @param a2s        A pointer to the array that contains the alignment to sequence position mapping
 *  @param Ss         A pointer to the array that contains the ungapped sequence
 *  @param circ       assume the molecules to be circular instead of linear (circ=0)
 */
DEPRECATED(void  alloc_sequence_arrays(const char     **sequences,
                                       short          ***S,
                                       short          ***S5,
                                       short          ***S3,
                                       unsigned short ***a2s,
                                       char           ***Ss,
                                       int            circ),
          "This function is obsolete");


/**
 *  @brief Free the memory of the sequence arrays used to deal with aligned sequences
 *
 *  This function frees the memory previously allocated with alloc_sequence_arrays()
 *
 *  @see alloc_sequence_arrays()
 *
 *  @ingroup  aln_utils_deprecated
 *  @param n_seq      The number of aligned sequences
 *  @param S          A pointer to the array of encoded sequences
 *  @param S5         A pointer to the array that contains the next 5' nucleotide of a sequence position
 *  @param S3         A pointer to the array that contains the next 3' nucleotide of a sequence position
 *  @param a2s        A pointer to the array that contains the alignment to sequence position mapping
 *  @param Ss         A pointer to the array that contains the ungapped sequence
 */
DEPRECATED(void  free_sequence_arrays(unsigned int    n_seq,
                                      short           ***S,
                                      short           ***S5,
                                      short           ***S3,
                                      unsigned short  ***a2s,
                                      char            ***Ss),
          "This fucntion is obsolete");

#endif

/**
 * @}
 */


#endif
