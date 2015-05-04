#ifndef __VIENNA_RNA_PACKAGE_ALIFOLD_H__
#define __VIENNA_RNA_PACKAGE_ALIFOLD_H__

#include "data_structures.h"

/**
 *  \file alifold.h
 *  \brief compute various properties (consensus MFE structures, partition function, Boltzmann distributed stochastic samples, ...) for 
 *  RNA sequence alignments
 */

/**
 *  \brief This variable controls the weight of the covariance term in the
 *  energy function of alignment folding algorithms
 * 
 *  Default is 1.
 */
extern  double  cv_fact;
/**
 *  \brief This variable controls the magnitude of the penalty for non-compatible sequences in
 *  the covariance term of alignment folding algorithms.
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
 *  \param strings    A pointer to a NULL terminated array of character arrays
 *  \param structure  A pointer to a character array that may contain a constraining consensus structure
 *                    (will be overwritten by a consensus structure that exhibits the MFE)
 *  \return           The free energy score in kcal/mol
 */
float  alifold(const char **strings, char *structure);


/**
 *  \brief Compute MFE and according structure of an alignment of sequences assuming the sequences are circular instead of linear
 * 
 *  \param strings    A pointer to a NULL terminated array of character arrays
 *  \param structure  A pointer to a character array that may contain a constraining consensus structure
 *                    (will be overwritten by a consensus structure that exhibits the MFE)
 *  \return           The free energy score in kcal/mol
 */
float  circalifold(const char **strings, char *structure);

/**
 *  \brief Free the memory occupied by MFE alifold functions
 */
void    free_alifold_arrays(void);

/**
 *  \brief Get the mean pairwise identity in steps from ?to?(ident)
 * 
 *  \param Alseq
 *  \param n_seq  The number of sequences in the alignment
 *  \param length The length of the alignment
 *  \param mini
 *  \return       The mean pairwise identity
 */
int get_mpi(char *Alseq[], int n_seq, int length, int *mini);

/**
 *  \brief Read a ribosum or other user-defined scoring matrix
 */
float   **readribosum(char *name);

/**
 *  \brief Calculate the free energy of a consensus structure given a set of aligned sequences
 * 
 *  \param  sequences   The NULL terminated array of sequences
 *  \param  structure   The consensus structure
 *  \param  n_seq       The number of sequences in the alignment
 *  \param  energy      A pointer to an array of at least two floats that will hold the free energies
 *                      (energy[0] will contain the free energy, energy[1] will be filled with the covariance energy term)
 *  \returns free energy in kcal/mol
 * 
 */
float   energy_of_alistruct(const char **sequences, const char *structure, int n_seq, float *energy);

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
 *  \param sequence The gapped sequence from the alignment
 *  \param S        pointer to an array that holds encoded sequence
 *  \param s5      pointer to an array that holds the next base 5' of alignment position i
 *  \param s3      pointer to an array that holds the next base 3' of alignment position i
 *  \param ss
 *  \param as
 *  \param circ    assume the molecules to be circular instead of linear (circ=0)
 */
void encode_ali_sequence(const char *sequence, short *S, short *s5, short *s3, char *ss, unsigned short *as, int circ);

/**
 *  \brief Allocate memory for sequence array used to deal with aligned sequences
 * 
 *  Note that these arrays will also be initialized according to the sequence alignment given
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
void  alloc_sequence_arrays(const char **sequences, short ***S, short ***S5, short ***S3, unsigned short ***a2s, char ***Ss, int circ);

/**
 *  \brief Free the memory of the sequence arrays used to deal with aligned sequences
 * 
 *  This function frees the memory previously allocated with alloc_sequence_arrays()
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
void  free_sequence_arrays(unsigned int n_seq, short ***S, short ***S5, short ***S3, unsigned short ***a2s, char ***Ss);

/*
#############################################################
# PARTITION FUNCTION VARIANTS OF THE ALIFOLD IMPLEMENTATION #
#############################################################
*/

/**
 *  \brief
 * 
 *  The partition function version of alifold() works in analogy to
 *  pf_fold(). Pair probabilities and information about sequence
 *  covariations are returned via the 'pi' variable as a list of
 *  #pair_info structs. The list is terminated by the first entry with
 *  pi.i = 0.
 * 
 *  \param sequences
 *  \param structure
 *  \param pl
 *  \return
 */
float alipf_fold(const char **sequences, char *structure, plist **pl);

/**
 *  \brief
 * 
 *  \param sequences
 *  \param structure
 *  \param pl
 *  \return
 */
float alipf_circ_fold(const char **sequences, char *structure, plist **pl);

/**
 *  \brief
 */
void  free_alipf_arrays(void);

/**
 *  \brief Sample a consensus secondary structure from the Boltzmann ensemble according its probability\n
 * 
 *  \param  prob  to be described (berni)
 *  \return       A sampled consensus secondary structure in dot-bracket notation
 */
char  *alipbacktrack(double *prob) ;

#endif
