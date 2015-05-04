#ifndef __VIENNA_RNA_PACKAGE_ALIFOLD_H__
#define __VIENNA_RNA_PACKAGE_ALIFOLD_H__

#include "data_structures.h"

/**
 *  \addtogroup consensus_fold
 *  \brief compute various properties (consensus MFE structures, partition function,
 *  Boltzmann distributed stochastic samples, ...) for RNA sequence alignments
 *
 *  Consensus structures can be predicted by a modified version of the
 *  fold() algorithm that takes a set of aligned sequences instead
 *  of a single sequence. The energy function consists of the mean energy
 *  averaged over the sequences, plus a covariance term that favors pairs
 *  with consistent and compensatory mutations and penalizes pairs that
 *  cannot be formed by all structures. For details see \cite hofacker:2002 and
 *  \cite bernhart:2008.
 *  @{
 *    \file alifold.h
 *    \brief compute various properties (consensus MFE structures, partition function, Boltzmann
 *    distributed stochastic samples, ...) for RNA sequence alignments
 *
 *  @}
 */


/**
 *  \addtogroup consensus_mfe_fold
 *  \ingroup consensus_fold
 *  @{
 *
 *  @}
 */

/**
 *  \brief This variable controls the weight of the covariance term in the
 *  energy function of alignment folding algorithms
 * 
 *  \ingroup consensus_fold
 * 
 *  Default is 1.
 */
extern  double  cv_fact;
/**
 *  \brief This variable controls the magnitude of the penalty for non-compatible sequences in
 *  the covariance term of alignment folding algorithms.
 * 
 *  \ingroup consensus_fold
 * 
 *  Default is 1.
 */
extern  double  nc_fact;

/*
##############################################
# MFE VARIANTS OF THE ALIFOLD IMPLEMENTATION #
##############################################
*/

/**
 *  \brief Update the energy parameters for alifold function
 * 
 *  Call this to recalculate the pair matrix and energy parameters after a
 *  change in folding parameters like #temperature
 */
void update_alifold_params(void);


/**
 *  \brief Compute MFE and according consensus structure of an alignment of sequences
 * 
 *  This function predicts the consensus structure for the aligned 'sequences'
 *  and returns the minimum free energy; the mfe structure in bracket
 *  notation is returned in 'structure'.
 * 
 *  Sufficient space must be allocated for 'structure' before calling
 *  alifold().
 * 
 *  \ingroup consensus_mfe_fold
 * 
 *  \param strings    A pointer to a NULL terminated array of character arrays
 *  \param structure  A pointer to a character array that may contain a constraining consensus structure
 *                    (will be overwritten by a consensus structure that exhibits the MFE)
 *  \return           The free energy score in kcal/mol
 */
float  alifold( const char **strings,
                char *structure);


/**
 *  \brief Compute MFE and according structure of an alignment of sequences assuming the sequences are circular instead of linear
 * 
 *  \ingroup consensus_mfe_fold
 * 
 *  \param strings    A pointer to a NULL terminated array of character arrays
 *  \param structure  A pointer to a character array that may contain a constraining consensus structure
 *                    (will be overwritten by a consensus structure that exhibits the MFE)
 *  \return           The free energy score in kcal/mol
 */
float  circalifold( const char **strings,
                    char *structure);

/**
 *  \brief Free the memory occupied by MFE alifold functions
 * 
 *  \ingroup consensus_mfe_fold
 * 
 */
void    free_alifold_arrays(void);

/**
 *  \brief Get the mean pairwise identity in steps from ?to?(ident)
 * 
 *  \ingroup consensus_fold
 * 
 *  \param Alseq
 *  \param n_seq  The number of sequences in the alignment
 *  \param length The length of the alignment
 *  \param mini
 *  \return       The mean pairwise identity
 */
int get_mpi(char *Alseq[],
            int n_seq,
            int length,
            int *mini);

/**
 *  \brief Read a ribosum or other user-defined scoring matrix
 * 
 *  \ingroup consensus_fold
 * 
 */
float   **readribosum(char *name);

/**
 *  \brief Calculate the free energy of a consensus structure given a set of aligned sequences
 * 
 *  \ingroup consensus_fold
 * 
 *  \param  sequences   The NULL terminated array of sequences
 *  \param  structure   The consensus structure
 *  \param  n_seq       The number of sequences in the alignment
 *  \param  energy      A pointer to an array of at least two floats that will hold the free energies
 *                      (energy[0] will contain the free energy, energy[1] will be filled with the covariance energy term)
 *  \returns free energy in kcal/mol
 * 
 */
float   energy_of_alistruct(const char **sequences,
                            const char *structure,
                            int n_seq,
                            float *energy);

float   energy_of_ali_gquad_structure(const char **sequences,
                                      const char *structure,
                                      int n_seq,
                                      float *energy);

/*
#############################################################
# some helper functions that might be useful in the library #
#############################################################
*/

/**
 *  \brief Get arrays with encoded sequence of the alignment
 *
 *  this function assumes that in S, S5, s3, ss and as enough
 *  space is already allocated (size must be at least sequence length+2)
 * 
 *  \ingroup consensus_fold
 * 
 *  \param sequence The gapped sequence from the alignment
 *  \param S        pointer to an array that holds encoded sequence
 *  \param s5      pointer to an array that holds the next base 5' of alignment position i
 *  \param s3      pointer to an array that holds the next base 3' of alignment position i
 *  \param ss
 *  \param as
 *  \param circ    assume the molecules to be circular instead of linear (circ=0)
 */
void encode_ali_sequence( const char *sequence,
                          short *S,
                          short *s5,
                          short *s3,
                          char *ss,
                          unsigned short *as,
                          int circ);

/**
 *  \brief Allocate memory for sequence array used to deal with aligned sequences
 * 
 *  Note that these arrays will also be initialized according to the sequence alignment given
 * 
 *  \ingroup consensus_fold
 * 
 *  \see free_sequence_arrays()
 * 
 *  \param sequences  The aligned sequences
 *  \param S          A pointer to the array of encoded sequences
 *  \param S5         A pointer to the array that contains the next 5' nucleotide of a sequence position
 *  \param S3         A pointer to the array that contains the next 3' nucleotide of a sequence position
 *  \param a2s        A pointer to the array that contains the alignment to sequence position mapping
 *  \param Ss         A pointer to the array that contains the ungapped sequence
 *  \param circ       assume the molecules to be circular instead of linear (circ=0)
 */
void  alloc_sequence_arrays(const char **sequences,
                            short ***S,
                            short ***S5,
                            short ***S3,
                            unsigned short ***a2s,
                            char ***Ss,
                            int circ);

/**
 *  \brief Free the memory of the sequence arrays used to deal with aligned sequences
 * 
 *  This function frees the memory previously allocated with alloc_sequence_arrays()
 * 
 *  \ingroup consensus_fold
 * 
 *  \see alloc_sequence_arrays()
 * 
 *  \param n_seq      The number of aligned sequences
 *  \param S          A pointer to the array of encoded sequences
 *  \param S5         A pointer to the array that contains the next 5' nucleotide of a sequence position
 *  \param S3         A pointer to the array that contains the next 3' nucleotide of a sequence position
 *  \param a2s        A pointer to the array that contains the alignment to sequence position mapping
 *  \param Ss         A pointer to the array that contains the ungapped sequence
 */
void  free_sequence_arrays( unsigned int n_seq,
                            short ***S,
                            short ***S5,
                            short ***S3,
                            unsigned short ***a2s,
                            char ***Ss);

/*
#############################################################
# PARTITION FUNCTION VARIANTS OF THE ALIFOLD IMPLEMENTATION #
#############################################################
*/


/**
 *  \addtogroup consensus_pf_fold
 *  \ingroup consensus_fold
 *  @{
 *
 *  @}
 */

/**
 *  \brief
 * 
 *  \ingroup consensus_pf_fold
 * 
 *  \param  sequences
 *  \param  structure
 *  \param  pl
 *  \param  parameters
 *  \param  calculate_bppm
 *  \param  is_constrained
 *  \param  is_circular
 *  \return
 */
float alipf_fold_par( const char **sequences,
                      char *structure,
                      plist **pl,
                      pf_paramT *parameters,
                      int calculate_bppm,
                      int is_constrained,
                      int is_circular);

/**
 *  \brief
 * 
 *  The partition function version of alifold() works in analogy to
 *  pf_fold(). Pair probabilities and information about sequence
 *  covariations are returned via the 'pi' variable as a list of
 *  #pair_info structs. The list is terminated by the first entry with
 *  pi.i = 0.
 * 
 *  \ingroup consensus_pf_fold
 * 
 *  \param sequences
 *  \param structure
 *  \param pl
 *  \return
 */
float alipf_fold( const char **sequences,
                  char *structure,
                  plist **pl);

/**
 *  \brief
 * 
 *  \ingroup consensus_pf_fold
 * 
 *  \param sequences
 *  \param structure
 *  \param pl
 *  \return
 */
float alipf_circ_fold(const char **sequences,
                      char *structure,
                      plist **pl);


/**
 *  \brief Get a pointer to the base pair probability array
 * 
 *  Accessing the base pair probabilities for a pair (i,j) is achieved by
 *  \verbatim FLT_OR_DBL *pr = export_bppm(); pr_ij = pr[iindx[i]-j]; \endverbatim
 * 
 *  \ingroup consensus_pf_fold
 * 
 *  \see get_iindx()
 *  \return A pointer to the base pair probability array
 */
FLT_OR_DBL *export_ali_bppm(void);

/**
 *  \brief
 *
 *  \ingroup consensus_pf_fold
 * 
 */
void  free_alipf_arrays(void);

/**
 *  \addtogroup consensus_stochbt
 *  @{
 *
 *  @}
 */

/**
 *  \brief Sample a consensus secondary structure from the Boltzmann ensemble according its probability\n
 * 
 *  \ingroup consensus_stochbt
 * 
 *  \param  prob  to be described (berni)
 *  \return       A sampled consensus secondary structure in dot-bracket notation
 */
char  *alipbacktrack(double *prob);


/**
 *  \brief Get pointers to (almost) all relavant arrays used in alifold's partition function computation
 *
 *  \ingroup consensus_fold
 *
 *  \note To obtain meaningful pointers, call alipf_fold first!
 *
 *  \see pf_alifold(), alipf_circ_fold()
 *
 *  \param S_p      A pointer to the 'S' array (integer representation of nucleotides)
 *  \param S5_p     A pointer to the 'S5' array
 *  \param S3_p     A pointer to the 'S3' array
 *  \param a2s_p    A pointer to the pair type matrix
 *  \param Ss_p     A pointer to the 'Ss' array
 *  \param qb_p     A pointer to the Q<sup>B</sup> matrix
 *  \param qm_p     A pointer to the Q<sup>M</sup> matrix
 *  \param q1k_p    A pointer to the 5' slice of the Q matrix (\f$q1k(k) = Q(1, k)\f$)
 *  \param qln_p    A pointer to the 3' slice of the Q matrix (\f$qln(l) = Q(l, n)\f$)
 *  \return         Non Zero if everything went fine, 0 otherwise
 */
int get_alipf_arrays(short ***S_p,
		     short ***S5_p,
		     short ***S3_p,
		     unsigned short ***a2s_p,
		     char ***Ss_p,
		     FLT_OR_DBL **qb_p,
		     FLT_OR_DBL **qm_p,
		     FLT_OR_DBL **q1k_p,
		     FLT_OR_DBL **qln_p,
		     short **pscore);

#endif
