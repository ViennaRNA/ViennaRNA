#ifndef __VIENNA_RNA_PACKAGE_ALIFOLD_H__
#define __VIENNA_RNA_PACKAGE_ALIFOLD_H__

#include <ViennaRNA/data_structures.h>
#include <ViennaRNA/params.h>
#include <ViennaRNA/ribo.h>

#ifdef __GNUC__
#define DEPRECATED(func) func __attribute__ ((deprecated))
#else
#define DEPRECATED(func) func
#endif

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
 *  \deprecated See #vrna_md_t.cv_fact, and vrna_ali_fold() to avoid using global variables
 *
 *  Default is 1.
 */
DEPRECATED(extern  double  cv_fact);
/**
 *  \brief This variable controls the magnitude of the penalty for non-compatible sequences in
 *  the covariance term of alignment folding algorithms.
 * 
 *  \ingroup consensus_fold
 * 
 *  \deprecated See #vrna_md_t.nc_fact, and vrna_ali_fold() to avoid using global variables
 *
 *  Default is 1.
 */
DEPRECATED(extern  double  nc_fact);

/*
##############################################
# MFE VARIANTS OF THE ALIFOLD IMPLEMENTATION #
##############################################
*/

/**
 *  \brief Compute MFE and according consensus structure of an alignment of sequences
 * 
 *  This function predicts the consensus structure for the alignment stored in
 *  vc and returns the minimum free energy; the mfe structure in bracket
 *  notation is returned in 'structure'.
 *  \note vc has to be of type #VRNA_VC_TYPE_ALIGNMENT.
 * 
 *  \note Sufficient space must be allocated for 'structure' before calling
 *  vrna_ali_fold(). Passing NULL to the 'structure' or setting #vrna_md_t.backtrack to
 *  0 turns of backtracing an no structure will be returned.
 * 
 *  \ingroup consensus_mfe_fold
 *
 *  \see vrna_get_fold_compound_ali()
 * 
 *  \param vc         The fold compound structure of type #VRNA_VC_TYPE_ALIGNMENT
 *  \param structure  A pointer to a character array that will be overwritten by
 *                    a consensus structure that exhibits the MFE. (Maybe NULL)
 *  \return           The free energy score in kcal/mol
 */
float vrna_ali_fold( vrna_fold_compound *vc,
                    char *structure);

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
 *  \deprecated Usage of this function is discouraged! Use vrna_ali_fold() instead
 *  \see vrna_ali_fold()
 *
 *  \param strings    A pointer to a NULL terminated array of character arrays
 *  \param structure  A pointer to a character array that may contain a constraining consensus structure
 *                    (will be overwritten by a consensus structure that exhibits the MFE)
 *  \return           The free energy score in kcal/mol
 */
DEPRECATED(float alifold( const char **strings, char *structure));

/**
 *  \brief Compute MFE and according structure of an alignment of sequences assuming the sequences are circular instead of linear
 * 
 *  \ingroup consensus_mfe_fold
 * 
 *  \deprecated Usage of this function is discouraged! Use vrna_ali_fold() instead
 *  \see vrna_ali_fold()
 * 
 *  \param strings    A pointer to a NULL terminated array of character arrays
 *  \param structure  A pointer to a character array that may contain a constraining consensus structure
 *                    (will be overwritten by a consensus structure that exhibits the MFE)
 *  \return           The free energy score in kcal/mol
 */
DEPRECATED(float circalifold( const char **strings, char *structure));

/**
 *  \brief Free the memory occupied by MFE alifold functions
 * 
 *  \deprecated Usage of this function is discouraged! It only
 *  affects memory being free'd that was allocated by an old API
 *  function before. Release of memory occupied by the newly introduced
 *  #vrna_fold_compound is handled by vrna_vrna_free_fold_compound()
 *
 *  \see vrna_vrna_free_fold_compound()
 *
 *  \ingroup consensus_mfe_fold
 * 
 */
DEPRECATED(void free_alifold_arrays(void));

/**
 *  \brief Calculate the free energy of a consensus structure given a set of aligned sequences
 * 
 *  \ingroup consensus_fold
 *
 *  \deprecated Usage of this function is discouraged! Use vrna_eval_structure(), and vrna_eval_covar_structure() instead!
 *
 *  \param  sequences   The NULL terminated array of sequences
 *  \param  structure   The consensus structure
 *  \param  n_seq       The number of sequences in the alignment
 *  \param  energy      A pointer to an array of at least two floats that will hold the free energies
 *                      (energy[0] will contain the free energy, energy[1] will be filled with the covariance energy term)
 *  \returns free energy in kcal/mol
 * 
 */
DEPRECATED(float energy_of_alistruct(const char **sequences, const char *structure, int n_seq, float *energy));

DEPRECATED(float energy_of_ali_gquad_structure(const char **sequences, const char *structure, int n_seq, float *energy));

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
 *  \brief Compute partition function and base pair probabilities for
 *  a sequence alignment.
 * 
 *  The partition function version of vrna_ali_fold() works in analogy to
 *  vrna_pf_fold(). Pair probabilities are returned via the 'pl' variable
 *  as a list of #plist structs. The list is terminated by the first entry
 *  with pl.i = 0.
 * 
 *  \ingroup consensus_pf_fold
 *
 *  \see vrna_ali_get_pair_info() for a replacement of pl with more detailed
 *  information
 * 
 *  \param vc         The #vrna_fold_compound of type #VRNA_VC_TYPE_ALIGNMENT
 *  \param structure  A pointer to a character array of length of the alignment (Maybe NULL)
 *  \param pl         A pointer to a #plist pointer where the pair probabilities are stored (Maybe NULL)
 *  \return           Gibbs free energy of the consensus fold space
 */
float vrna_ali_pf_fold(vrna_fold_compound *vc,
                      char *structure,
                      plist **pl);

/**
 *  \brief Retrieve an array of #pair_info structures from precomputed pair probabilities
 *
 *  This array of structures contains information about positionwise pair probabilies,
 *  base pair entropy and more
 *
 *  \see #pair_info, and vrna_ali_pf_fold()
 *
 *  \param  vc          The #vrna_fold_compound of type #VRNA_VC_TYPE_ALIGNMENT with precomputed partition function matrices
 *  \param  structure   An optional structure in dot-bracket notation (Maybe NULL)
 *  \param  threshold   Do not include results with pair probabilities below threshold
 *  \return             The #pair_info array
 */
pair_info *vrna_ali_get_pair_info(vrna_fold_compound *vc,
                                  const char *structure,
                                  double threshold);

/**
 *  \brief
 * 
 *  \ingroup consensus_pf_fold
 * 
 *  \deprecated Use vrna_ali_pf_fold() instead
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
DEPRECATED(float alipf_fold_par( const char **sequences,
                      char *structure,
                      plist **pl,
                      pf_paramT *parameters,
                      int calculate_bppm,
                      int is_constrained,
                      int is_circular));

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
 *  \deprecated Use vrna_ali_pf_fold() instead
 *
 *  \param sequences
 *  \param structure
 *  \param pl
 *  \return
 */
DEPRECATED(float alipf_fold( const char **sequences, char *structure, plist **pl));

/**
 *  \brief
 * 
 *  \ingroup consensus_pf_fold
 *
 *  \deprecated Use vrna_ali_pf_fold() instead
 * 
 *  \param sequences
 *  \param structure
 *  \param pl
 *  \return
 */
DEPRECATED(float alipf_circ_fold(const char **sequences, char *structure, plist **pl));


/**
 *  \brief Get a pointer to the base pair probability array
 * 
 *  Accessing the base pair probabilities for a pair (i,j) is achieved by
 *  \verbatim FLT_OR_DBL *pr = export_bppm(); pr_ij = pr[iindx[i]-j]; \endverbatim
 * 
 *  \ingroup consensus_pf_fold
 *
 *  \deprecated Usage of this function is discouraged! The new #vrna_fold_compound
 *  allows direct access to the folding matrices, including the pair probabilities!
 *  The pair probability array returned here reflects the one of the latest call
 *  to vrna_ali_pf_fold(), or any of the old API calls for consensus structure
 *  partition function folding.
 * 
 *  \see #vrna_fold_compound, vrna_get_fold_compound_ali(), and vrna_ali_pf_fold()
 *
 *  \return A pointer to the base pair probability array
 */
DEPRECATED(FLT_OR_DBL *export_ali_bppm(void));

/**
 *  \brief Free the memory occupied by folding matrices allocated by alipf_fold, alipf_circ_fold, etc.
 *
 *  \ingroup consensus_pf_fold
 * 
 *  \deprecated Usage of this function is discouraged! This function only free's memory
 *  allocated by old API function calls. Memory allocated by any of the new API calls (starting with vrna_)
 *  will be not affected!
 *
 *  \see #vrna_fold_compound, vrna_vrna_free_fold_compound()
 *
 */
DEPRECATED(void  free_alipf_arrays(void));

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
 *  \see vrna_ali_pf_fold() for precomputing the partition function matrices, and
 *
 *  \param  vc    The #vrna_fold_compound of type #VRNA_VC_TYPE_ALIGNMENT with precomputed partition function matrices
 *  \param  prob  to be described (berni)
 *  \return       A sampled consensus secondary structure in dot-bracket notation
 */
char *vrna_ali_pbacktrack(vrna_fold_compound *vc, double *prob);

/**
 *  \brief Sample a consensus secondary structure from the Boltzmann ensemble according its probability\n
 * 
 *  \ingroup consensus_stochbt
 *
 *  \deprecated Use vrna_ali_pbacktrack() instead!
 *
 *  \param  prob  to be described (berni)
 *  \return       A sampled consensus secondary structure in dot-bracket notation
 */
DEPRECATED(char  *alipbacktrack(double *prob));

/**
 *  \brief Get pointers to (almost) all relavant arrays used in alifold's partition function computation
 *
 *  \ingroup consensus_fold
 *
 *  \note To obtain meaningful pointers, call alipf_fold first!
 *
 *  \see pf_alifold(), alipf_circ_fold()
 *
 *  \deprecated It is discouraged to use this function! The new #vrna_fold_compound allows
 *  direct access to all necessary consensus structure prediction related variables!
 *
 *  \see #vrna_fold_compound, vrna_get_fold_compound_ali(), vrna_ali_pf_fold()
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
DEPRECATED(int get_alipf_arrays(short ***S_p,
                     short ***S5_p,
                     short ***S3_p,
                     unsigned short ***a2s_p,
                     char ***Ss_p,
                     FLT_OR_DBL **qb_p,
                     FLT_OR_DBL **qm_p,
                     FLT_OR_DBL **q1k_p,
                     FLT_OR_DBL **qln_p,
                     int **pscore));


/**
 *  \brief Update the energy parameters for alifold function
 *
 *  Call this to recalculate the pair matrix and energy parameters after a
 *  change in folding parameters like #temperature
 *
 *  \ingroup  consensus_fold
 *  \deprecated Usage of this function is discouraged! The new API uses #vrna_fold_compound
 *  to lump all folding related necessities together, including the energy parameters. Use
 *  vrna_update_fold_params() to update the energy parameters within a #vrna_fold_compound.
 */
DEPRECATED(void update_alifold_params(void));

#endif
