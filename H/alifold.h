#ifndef __VIENNA_RNA_PACKAGE_ALIFOLD_H__
#define __VIENNA_RNA_PACKAGE_ALIFOLD_H__

#include "data_structures.h"

/**
*** \file alifold.h
**/

/**
*** covariance scaling factor (default = 1.)
**/
extern  double  cv_fact;
/** */
extern  double  nc_fact /* =1 */;

/*
##############################################
# MFE VARIANTS OF THE ALIFOLD IMPLEMENTATION #
##############################################
*/

/**
*** Update the energy parameters for alifold function
**/
void update_alifold_params(void);


/**
*** Compute MFE and according structure of an alignment of sequences
***
*** \param strings    A pointer to a NULL terminated array of character arrays
*** \param structure  A pointer to a character array that may contain a constraining consensus structure
***                   (will be overwritten by a consensus structure that exhibits the MFE)
*** \return           The minimum free energy in kcal/mol
**/
float  alifold(const char **strings, char *structure);


/**
*** Compute MFE and according structure of an alignment of sequences assuming the sequences are circular instead of linear
***
*** \param strings    A pointer to a NULL terminated array of character arrays
*** \param structure  A pointer to a character array that may contain a constraining consensus structure
***                   (will be overwritten by a consensus structure that exhibits the MFE)
*** \return           The minimum free energy in kcal/mol
**/
float  circalifold(const char **strings, char *structure);

/**
*** Free the memory occupied by MFE alifold functions
**/
void    free_alifold_arrays(void);

/**
*** Get the mean pairwise identity in steps from ?to?(ident)
*** \param Alseq
*** \param n_seq  The number of sequences in the alignment
*** \param length The length of the alignment
*** \param mini   
*** \return       The mean pairwise identity
**/
int get_mpi(char *Alseq[], int n_seq, int length, int *mini);

/**
*** how to chose a ribosum matrix:<BR>
*** ribosum matrices exist for starlike clusters of<BR>
*** X=45 55 60 65 70 75 80 85 90 95 100<BR>
*** they are further seperated by only regarding sequences with a minimal pairwise idensity of Y=25-95, step 5, not all are present.<BR>
*** now the question is, which matrix to use when.<BR>
*** the suggestion of the dr. will:<BR>
*** with a mpi of Z<BR>
*** X~Z and Y > Z ??<BR>
*** if we say the minimum of the pis is M,<BR>
*** then we may be able to use:<BR>
*** X~Z and Y ~ M?<BR>
*** I'd say we do a default matrix (e.g. 85/60) but better is try out (all 170??)<BR>
*** and then we use the best in average.<BR>
*** actually, it would be preferrable to make a very big testset and simply check it out.<BR>
*** (create a function to derive the best matrix)<BR>
*** furthermore:<BR>
*** default, function or user defined.<BR>
*** <BR>
*** ntscd:<BR>
*** fijklmn<BR>
*** pijpklpmn<BR>
***
**/
float   **readribosum(char *name);

/**
*** Calculate the free energy of a consensus structure given a set of aligned sequences
***
*** \param  sequences   The NULL terminated array of sequences
*** \param  structure   The consensus structure
*** \param  n_seq       The number of sequences in the alignment
*** \param  energy      A pointer to an array of at least two floats that will hold the free energies
***                     (energy[0] will contain the free energy, energy[1] will be filled with the covariance energy term)
*** \returns free energy in kcal/mol
***
**/
float   energy_of_alistruct(const char **sequences, const char *structure, int n_seq, float *energy);

/*
#############################################################
# some helper functions that might be useful in the library #
#############################################################
*/

/**
*** Get arrays with encoded sequence of the alignment
***
*** this function assumes that in S, S5, s3, ss and as enough
*** space is already allocated (size must be at least sequence length+2)
*** \param sequence The gapped sequence from the alignment
*** \param S        pointer to an array that holds encoded sequence
*** \param s5      pointer to an array that holds the next base 5' of alignment position i
*** \param s3      pointer to an array that holds the next base 3' of alignment position i
*** \param ss      
*** \param as      
*** \param circ    assume the molecules to be circular instead of linear (circ=0)
**/
void encode_ali_sequence(const char *sequence, short *S, short *s5, short *s3, char *ss, unsigned short *as, int circ);

/**
*** Allocate memory for sequence array used to deal with aligned sequences
*** Note that these arrays will also be initialized according to the sequence alignment given
***
*** \see free_sequence_arrays()
***
*** \param sequences  The aligned sequences
*** \param S          A pointer to the array of encoded sequences
*** \param S5         A pointer to the array that contains the next 5' nucleotide of a sequence position
*** \param S3         A pointer to the array that contains the next 3' nucleotide of a sequence position
*** \param a2s        A pointer to the array that contains the alignment to sequence position mapping
*** \param Ss         A pointer to the array that contains the ungapped sequence
*** \param circ       assume the molecules to be circular instead of linear (circ=0)
**/
void  alloc_sequence_arrays(const char **sequences, short ***S, short ***S5, short ***S3, unsigned short ***a2s, char ***Ss, int circ);

/** Free the memory of the sequence arrays used to deal with aligned sequences
*** This function frees the memory previously allocated with alloc_sequence_arrays()
***
*** \see alloc_sequence_arrays()
***
*** \param n_seq      The number of aligned sequences
*** \param S          A pointer to the array of encoded sequences
*** \param S5         A pointer to the array that contains the next 5' nucleotide of a sequence position
*** \param S3         A pointer to the array that contains the next 3' nucleotide of a sequence position
*** \param a2s        A pointer to the array that contains the alignment to sequence position mapping
*** \param Ss         A pointer to the array that contains the ungapped sequence
**/
void  free_sequence_arrays(unsigned int n_seq, short ***S, short ***S5, short ***S3, unsigned short ***a2s, char ***Ss);

/*
#############################################################
# PARTITION FUNCTION VARIANTS OF THE ALIFOLD IMPLEMENTATION #
#############################################################
*/

float alipf_fold(const char **sequences, char *structure, plist **pl);
float alipf_circ_fold(const char **sequences, char *structure, plist **pl);
void  free_alipf_arrays(void);
char  *alipbacktrack(double *prob) ;

#endif
