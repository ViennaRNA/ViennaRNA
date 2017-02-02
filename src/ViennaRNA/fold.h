#ifndef VIENNA_RNA_PACKAGE_FOLD_H
#define VIENNA_RNA_PACKAGE_FOLD_H

#include <ViennaRNA/data_structures.h>
#include <ViennaRNA/params.h>
#include <ViennaRNA/mfe.h>
#include <ViennaRNA/eval.h>

#ifdef VRNA_WARN_DEPRECATED
# ifdef __GNUC__
#  define DEPRECATED(func) func __attribute__ ((deprecated))
# else
#  define DEPRECATED(func) func
# endif
#else
# define DEPRECATED(func) func
#endif

/* make this interface backward compatible with RNAlib < 2.2.0 */
#define VRNA_BACKWARD_COMPAT

/**
 *  @file fold.h
 *  @ingroup  mfe_fold
 *  @brief MFE calculations for single RNA sequences
 */

/**
 *  @addtogroup mfe_fold_single
 *  @brief This module contains all functions and variables related to the calculation
 *  of global minimum free energy structures for single sequences.
 *
 *  The library provides a fast dynamic programming minimum free energy
 *  folding algorithm as described by "Zuker & Stiegler (1981)" @cite zuker:1981.
 *
 *  @{
 *  @ingroup  mfe_fold_single
 */

/**
 *  @brief Compute Minimum Free Energy (MFE), and a corresponding secondary structure for an RNA sequence
 *
 *  This simplified interface to vrna_mfe() computes the MFE and, if required, a secondary structure for an
 *  RNA sequence using default options. Memory required for dynamic programming (DP) matrices will
 *  be allocated and free'd on-the-fly. Hence, after return of this function, the recursively filled
 *  matrices are not available any more for any post-processing, e.g. suboptimal backtracking, etc.
 *
 *  @note In case you want to use the filled DP matrices for any subsequent post-processing step, or
 *  you require other conditions than specified by the default model details, use vrna_mfe(),
 *  and the data structure #vrna_fold_compound_t instead.
 *
 *  @see vrna_circfold(), vrna_mfe(), vrna_fold_compound(), #vrna_fold_compound_t
 *
 *  @param sequence   RNA sequence
 *  @param structure  A pointer to the character array where the
 *         secondary structure in dot-bracket notation will be written to
 *  @return the minimum free energy (MFE) in kcal/mol
 */
float
vrna_fold(const char *sequence,
          char *structure);

/**
 *  @brief Compute Minimum Free Energy (MFE), and a corresponding secondary structure for a circular RNA sequence
 *
 *  This simplified interface to vrna_mfe() computes the MFE and, if required, a secondary structure for a
 *  circular RNA sequence using default options. Memory required for dynamic programming (DP) matrices will
 *  be allocated and free'd on-the-fly. Hence, after return of this function, the recursively filled
 *  matrices are not available any more for any post-processing, e.g. suboptimal backtracking, etc.
 *
 *  Folding of circular RNA sequences is handled as a post-processing step of the forward
 *  recursions. See @cite hofacker:2006 for further details.
 *
 *  @note In case you want to use the filled DP matrices for any subsequent post-processing step, or
 *  you require other conditions than specified by the default model details, use vrna_mfe(),
 *  and the data structure #vrna_fold_compound_t instead.
 *
 *  @see vrna_fold(), vrna_mfe(), vrna_fold_compound(), #vrna_fold_compound_t
 *
 *  @param sequence   RNA sequence
 *  @param structure  A pointer to the character array where the
 *         secondary structure in dot-bracket notation will be written to
 *  @return the minimum free energy (MFE) in kcal/mol
 */
float
vrna_circfold(const char *sequence,
              char *structure);

#ifdef  VRNA_BACKWARD_COMPAT

/**
 *  @brief Compute minimum free energy and an appropriate secondary
 *  structure of an RNA sequence
 *
 *  The first parameter given, the RNA sequence, must be @a uppercase and should only contain
 *  an alphabet @f$\Sigma@f$ that is understood by the RNAlib\n
 *  (e.g. @f$ \Sigma = \{A,U,C,G\} @f$)\n
 *
 *  The second parameter, @a structure, must always point to an allocated
 *  block of memory with a size of at least @f$\mathrm{strlen}(\mathrm{sequence})+1@f$
 *
 *  If the third parameter is NULL, global model detail settings are assumed for the folding
 *  recursions. Otherwise, the provided parameters are used.
 *
 *  The fourth parameter indicates whether a secondary structure constraint in enhanced dot-bracket
 *  notation is passed through the structure parameter or not. If so, the characters " | x < > " are
 *  recognized to mark bases that are paired, unpaired, paired upstream, or downstream, respectively.
 *  Matching brackets " ( ) " denote base pairs, dots "." are used for unconstrained bases.
 *
 *  To indicate that the RNA sequence is circular and thus has to be post-processed, set the last
 *  parameter to non-zero
 *
 *  After a successful call of fold_par(), a backtracked secondary structure (in dot-bracket notation)
 *  that exhibits the minimum of free energy will be written to the memory @a structure is pointing to.
 *  The function returns the minimum of free energy for any fold of the sequence given.
 *
 *  @note OpenMP: Passing NULL to the 'parameters' argument involves access to several global model
 *        detail variables and thus is not to be considered threadsafe
 *
 *  @deprecated use vrna_mfe() instead!
 *
 *  @see vrna_mfe(), fold(), circfold(), #vrna_md_t, set_energy_model(), get_scaled_parameters()
 *
 *  @param sequence       RNA sequence
 *  @param structure      A pointer to the character array where the
 *                        secondary structure in dot-bracket notation will be written to
 *  @param parameters     A data structure containing the pre-scaled energy contributions
 *                        and the model details. (NULL may be passed, see OpenMP notes above)
 *  @param is_constrained Switch to indicate that a structure constraint is passed via the structure argument (0==off)
 *  @param is_circular    Switch to (de-)activate post-processing steps in case RNA sequence is circular (0==off)
 *
 *  @return the minimum free energy (MFE) in kcal/mol
 */
DEPRECATED(float
fold_par( const char *sequence,
          char *structure,
          vrna_param_t *parameters,
          int is_constrained,
          int is_circular));

/**
 *  @brief Compute minimum free energy and an appropriate secondary structure of an RNA sequence
 *
 *  This function essentially does the same thing as fold_par(). However, it takes its model details,
 *  i.e. #temperature, #dangles, #tetra_loop, #noGU, #no_closingGU, #fold_constrained, #noLonelyPairs
 *  from the current global settings within the library
 *
 *  @deprecated use vrna_fold(), or vrna_mfe() instead!
 *
 *  @see fold_par(), circfold()
 *
 *  @param sequence RNA sequence
 *  @param structure A pointer to the character array where the
 *         secondary structure in dot-bracket notation will be written to
 *  @return the minimum free energy (MFE) in kcal/mol
 */
DEPRECATED(float fold( const char *sequence, char *structure));

/**
 *  @brief Compute minimum free energy and an appropriate secondary structure of a circular RNA sequence
 *
 *  This function essentially does the same thing as fold_par(). However, it takes its model details,
 *  i.e. #temperature, #dangles, #tetra_loop, #noGU, #no_closingGU, #fold_constrained, #noLonelyPairs
 *  from the current global settings within the library
 *
 *  @deprecated Use vrna_circfold(), or vrna_mfe() instead!
 *
 *  @see fold_par(), circfold()
 *
 *  @param sequence RNA sequence
 *  @param structure A pointer to the character array where the
 *         secondary structure in dot-bracket notation will be written to
 *  @return the minimum free energy (MFE) in kcal/mol
 */
DEPRECATED(float circfold( const char *sequence, char *structure));


/**
 *  @brief Free arrays for mfe folding
 *
 *  @deprecated See vrna_fold(), vrna_circfold(), or vrna_mfe() and #vrna_fold_compound_t for the usage of the new API!
 *
 */
DEPRECATED(void free_arrays(void));



/**
 *  @brief Recalculate energy parameters
 *
 *  @deprecated For non-default model settings use the new API with vrna_params_subst() and vrna_mfe() instead!
 *
 */
DEPRECATED(void update_fold_params(void));

/**
 *  @brief Recalculate energy parameters
 *
 *  @deprecated For non-default model settings use the new API with vrna_params_subst() and vrna_mfe() instead!
 *
 */
DEPRECATED(void update_fold_params_par(vrna_param_t *parameters));

/**
 *
 *  @deprecated See vrna_mfe() and #vrna_fold_compound_t for the usage of the new API!
 *
 */
DEPRECATED(void
export_fold_arrays( int **f5_p,
                    int **c_p,
                    int **fML_p,
                    int **fM1_p,
                    int **indx_p,
                    char **ptype_p));

/**
 *
 *  @deprecated See vrna_mfe() and #vrna_fold_compound_t for the usage of the new API!
 *
 */
DEPRECATED(void
export_fold_arrays_par( int **f5_p,
                        int **c_p,
                        int **fML_p,
                        int **fM1_p,
                        int **indx_p,
                        char **ptype_p,
                        vrna_param_t **P_p));

/**
 *
 *  @deprecated See vrna_mfe() and #vrna_fold_compound_t for the usage of the new API!
 *
 */
DEPRECATED(void
export_circfold_arrays( int *Fc_p,
                        int *FcH_p,
                        int *FcI_p,
                        int *FcM_p,
                        int **fM2_p,
                        int **f5_p,
                        int **c_p,
                        int **fML_p,
                        int **fM1_p,
                        int **indx_p,
                        char **ptype_p));

/**
 *
 *  @deprecated See vrna_mfe() and #vrna_fold_compound_t for the usage of the new API!
 *
 */
DEPRECATED(void
export_circfold_arrays_par( int *Fc_p,
                            int *FcH_p,
                            int *FcI_p,
                            int *FcM_p,
                            int **fM2_p,
                            int **f5_p,
                            int **c_p,
                            int **fML_p,
                            int **fM1_p,
                            int **indx_p,
                            char **ptype_p,
                            vrna_param_t **P_p));



/* finally moved the loop energy function declarations to this header...  */
/* BUT: The functions only exist for backward compatibility reasons!      */
/* You better include "loop_energies.h" and call the functions:           */
/* E_Hairpin() and E_IntLoop() which are (almost) threadsafe as they get  */
/* a pointer to the energy parameter data structure as additional argument */

/**
 *  @deprecated {This function is deprecated and will be removed soon.
 *  Use @ref E_IntLoop() instead!}
 */
DEPRECATED(int LoopEnergy(int n1,
                          int n2,
                          int type,
                          int type_2,
                          int si1,
                          int sj1,
                          int sp1,
                          int sq1));

/**
 *  @deprecated {This function is deprecated and will be removed soon.
 *  Use @ref E_Hairpin() instead!}
 */
DEPRECATED(int HairpinE(int size,
                        int type,
                        int si1,
                        int sj1,
                        const char *string));

/**
 *  Allocate arrays for folding\n
 *  @deprecated See vrna_mfe() and #vrna_fold_compound_t for the usage of the new API!
 *
 */
DEPRECATED(void initialize_fold(int length));

/**
 *
 */
DEPRECATED(char *backtrack_fold_from_pair(char *sequence,
                                          int i,
                                          int j));


#endif

/**
 *  @}
 */

#endif
