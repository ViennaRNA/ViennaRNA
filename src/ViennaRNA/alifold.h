#ifndef VIENNA_RNA_PACKAGE_ALIFOLD_H
#define VIENNA_RNA_PACKAGE_ALIFOLD_H

#include <ViennaRNA/datastructures/basic.h>
#include <ViennaRNA/params/basic.h>
#include <ViennaRNA/params/ribosum.h>
#include <ViennaRNA/mfe/global.h>
#include <ViennaRNA/backtrack/global.h>
#include <ViennaRNA/partfunc/global.h>
#include <ViennaRNA/sequences/alignments.h>
#include <ViennaRNA/structures/problist.h>
#include <ViennaRNA/sampling/basic.h>

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
 *  @file alifold.h
 *  @ingroup mfe_global_deprecated
 *  @brief Functions for comparative structure prediction using RNA sequence alignments
 *
 */

#ifndef VRNA_DISABLE_BACKWARD_COMPATIBILITY

/*
#################################################
# DEPRECATED FUNCTIONS                          #
#################################################
*/

/**
 * @ingroup mfe_global_deprecated
 * @{
 */

/**
 *  @brief Compute MFE and according consensus structure of an alignment of sequences
 * 
 *  This function predicts the consensus structure for the aligned 'sequences'
 *  and returns the minimum free energy; the mfe structure in bracket
 *  notation is returned in 'structure'.
 * 
 *  Sufficient space must be allocated for 'structure' before calling
 *  alifold().
 * 
 *  @deprecated Usage of this function is discouraged! Use vrna_alifold(), or vrna_mfe() instead!
 *
 *  @see vrna_alifold(), vrna_mfe()
 *
 *  @param strings    A pointer to a NULL terminated array of character arrays
 *  @param structure  A pointer to a character array that may contain a constraining consensus structure
 *                    (will be overwritten by a consensus structure that exhibits the MFE)
 *  @return           The free energy score in kcal/mol
 */
DEPRECATED(float alifold( const char **strings, char *structure),
          "Use vrna_alifold() or vrna_mfe() instead");

/**
 *  @brief Compute MFE and according structure of an alignment of sequences assuming the sequences are circular instead of linear
 * 
 *  @deprecated Usage of this function is discouraged! Use vrna_alicircfold(), and vrna_mfe() instead!
 *
 *  @see vrna_alicircfold(), vrna_alifold(), vrna_mfe()
 * 
 *  @param strings    A pointer to a NULL terminated array of character arrays
 *  @param structure  A pointer to a character array that may contain a constraining consensus structure
 *                    (will be overwritten by a consensus structure that exhibits the MFE)
 *  @return           The free energy score in kcal/mol
 */
DEPRECATED(float circalifold( const char **strings, char *structure),
          "Use vrna_alicircfold() or vrna_mfe() instead");

/**
 *  @brief Free the memory occupied by MFE alifold functions
 * 
 *  @deprecated Usage of this function is discouraged! It only
 *  affects memory being free'd that was allocated by an old API
 *  function before. Release of memory occupied by the newly introduced
 *  #vrna_fold_compound_t is handled by vrna_fold_compound_free()
 *
 *  @see vrna_fold_compound_free()
 * 
 */
DEPRECATED(void free_alifold_arrays(void),
          "This function is obsolete");

/* End group mfe_global_deprecated */
/**@}*/

/**
 *  @brief Calculate the free energy of a consensus structure given a set of aligned sequences
 * 
 *  @ingroup consensus_fold
 *
 *  @deprecated Usage of this function is discouraged! Use vrna_eval_structure(), and vrna_eval_covar_structure() instead!
 *
 *  @param  sequences   The NULL terminated array of sequences
 *  @param  structure   The consensus structure
 *  @param  n_seq       The number of sequences in the alignment
 *  @param  energy      A pointer to an array of at least two floats that will hold the free energies
 *                      (energy[0] will contain the free energy, energy[1] will be filled with the covariance energy term)
 *  @returns free energy in kcal/mol
 * 
 */
DEPRECATED(float energy_of_alistruct(const char **sequences, const char *structure, int n_seq, float *energy),
           "Use vrna_eval_structure() and vrna_eval_covar_structure() instead");

DEPRECATED(float energy_of_ali_gquad_structure(const char **sequences, const char *structure, int n_seq, float *energy),
          "Use vrna_eval_structure() and vrna_eval_covar_structure() instead");

/**
 *  @brief This variable controls the weight of the covariance term in the
 *  energy function of alignment folding algorithms
 * 
 *  @ingroup consensus_fold
 *
 *  @deprecated See #vrna_md_t.cv_fact, and vrna_mfe() to avoid using global variables
 *
 *  Default is 1.
 */
DEPRECATED(extern  double  cv_fact,
          "Use the cv_fact attribute of the vrna_md_t datastructure instead");
/**
 *  @brief This variable controls the magnitude of the penalty for non-compatible sequences in
 *  the covariance term of alignment folding algorithms.
 * 
 *  @ingroup consensus_fold
 * 
 *  @deprecated See #vrna_md_t.nc_fact, and vrna_mfe() to avoid using global variables
 *
 *  Default is 1.
 */
DEPRECATED(extern  double  nc_fact,
          "Use the nc_fact attribute of the vrna_md_t datastructure instead");

/**
 * @ingroup part_func_global_deprecated
 * @{
 */

/**
 *  @brief
 * 
 *  @deprecated Use vrna_pf() instead
 *
 *  @param  sequences
 *  @param  structure
 *  @param  pl
 *  @param  parameters
 *  @param  calculate_bppm
 *  @param  is_constrained
 *  @param  is_circular
 *  @return
 */
DEPRECATED(float alipf_fold_par( const char **sequences,
                      char *structure,
                      vrna_ep_t **pl,
                      vrna_exp_param_t *parameters,
                      int calculate_bppm,
                      int is_constrained,
                      int is_circular),
          "Use vrna_pf_alifold() or vrna_pf() instead");

/**
 *  @brief
 * 
 *  The partition function version of alifold() works in analogy to
 *  pf_fold(). Pair probabilities and information about sequence
 *  covariations are returned via the 'pi' variable as a list of
 *  #vrna_pinfo_t structs. The list is terminated by the first entry with
 *  pi.i = 0.
 * 
 *  @deprecated Use vrna_pf() instead
 *
 *  @param sequences
 *  @param structure
 *  @param pl
 *  @return
 */
DEPRECATED(float alipf_fold( const char **sequences, char *structure, vrna_ep_t **pl),
          "Use vrna_pf_alifold() or vrna_pf() instead");

/**
 *  @brief
 *
 *  @deprecated Use vrna_pf() instead
 * 
 *  @param sequences
 *  @param structure
 *  @param pl
 *  @return
 */
DEPRECATED(float alipf_circ_fold(const char **sequences, char *structure, vrna_ep_t **pl),
          "Use vrna_pf_circalifold() or vrna_pf() instead");


/**
 *  @brief Get a pointer to the base pair probability array
 * 
 *  Accessing the base pair probabilities for a pair (i,j) is achieved by
 *  @verbatim FLT_OR_DBL *pr = export_bppm(); pr_ij = pr[iindx[i]-j]; @endverbatim
 *
 *  @deprecated Usage of this function is discouraged! The new #vrna_fold_compound_t
 *  allows direct access to the folding matrices, including the pair probabilities!
 *  The pair probability array returned here reflects the one of the latest call
 *  to vrna_pf(), or any of the old API calls for consensus structure
 *  partition function folding.
 * 
 *  @see #vrna_fold_compound_t, vrna_fold_compound_comparative(), and vrna_pf()
 *
 *  @return A pointer to the base pair probability array
 */
DEPRECATED(FLT_OR_DBL *export_ali_bppm(void),
          "Use the new API with vrna_fold_compound_t datastructure instead");

/**
 *  @brief Free the memory occupied by folding matrices allocated by alipf_fold, alipf_circ_fold, etc.
 *
 *  @deprecated Usage of this function is discouraged! This function only free's memory
 *  allocated by old API function calls. Memory allocated by any of the new API calls (starting with vrna_)
 *  will be not affected!
 *
 *  @see #vrna_fold_compound_t, vrna_vrna_fold_compound_free()
 *
 */
DEPRECATED(void  free_alipf_arrays(void),
          "This function is obsolete");

/**
 *  @brief Sample a consensus secondary structure from the Boltzmann ensemble according its probability
 * 
 *  @deprecated Use vrna_pbacktrack() instead!
 *
 *  @param  prob  to be described (berni)
 *  @return       A sampled consensus secondary structure in dot-bracket notation
 */
DEPRECATED(char  *alipbacktrack(double *prob),
          "Use the new API and vrna_pbacktrack() instead");

/**
 *  @brief Get pointers to (almost) all relavant arrays used in alifold's partition function computation
 *
 *  @note To obtain meaningful pointers, call alipf_fold first!
 *
 *  @see  #vrna_fold_compound_t, vrna_fold_compound_comparative(), vrna_pf(),
 *        pf_alifold(), alipf_circ_fold()
 *
 *  @deprecated It is discouraged to use this function! The new #vrna_fold_compound_t allows
 *              direct access to all necessary consensus structure prediction related variables!
 *
 *  @param S_p      A pointer to the 'S' array (integer representation of nucleotides)
 *  @param S5_p     A pointer to the 'S5' array
 *  @param S3_p     A pointer to the 'S3' array
 *  @param a2s_p    A pointer to the alignment-column to sequence position mapping array
 *  @param Ss_p     A pointer to the 'Ss' array
 *  @param qb_p     A pointer to the Q<sup>B</sup> matrix
 *  @param qm_p     A pointer to the Q<sup>M</sup> matrix
 *  @param q1k_p    A pointer to the 5' slice of the Q matrix (@f$q1k(k) = Q(1, k)@f$)
 *  @param qln_p    A pointer to the 3' slice of the Q matrix (@f$qln(l) = Q(l, n)@f$)
 *  @param pscore   A pointer to the start of a pscore list
 *  @return         Non Zero if everything went fine, 0 otherwise
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
                     short **pscore),
          "Use the new API with vrna_fold_compound_t datastructure instead");


/* End group part_func_global_deprecated */
/**@}*/

/**
 *  @brief Update the energy parameters for alifold function
 *
 *  Call this to recalculate the pair matrix and energy parameters after a
 *  change in folding parameters like #temperature
 *
 *  @ingroup  consensus_fold
 *  @deprecated Usage of this function is discouraged! The new API uses #vrna_fold_compound_t
 *  to lump all folding related necessities together, including the energy parameters. Use
 *  vrna_update_fold_params() to update the energy parameters within a #vrna_fold_compound_t.
 */
DEPRECATED(void update_alifold_params(void),
          "Use the new API with vrna_fold_compound_t datastructure instead");

#endif


#endif
