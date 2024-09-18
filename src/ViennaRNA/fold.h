#ifndef VIENNA_RNA_PACKAGE_FOLD_H
#define VIENNA_RNA_PACKAGE_FOLD_H

#ifndef VRNA_DISABLE_BACKWARD_COMPATIBILITY

#include <ViennaRNA/datastructures/basic.h>
#include <ViennaRNA/params/basic.h>
#include <ViennaRNA/mfe/global.h>
#include <ViennaRNA/backtrack/global.h>
#include <ViennaRNA/eval/structures.h>

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
 *  @file fold.h
 *  @ingroup  mfe_global_deprecated
 *  @brief MFE calculations for single RNA sequences
 */

/**
 *  @addtogroup mfe_global_deprecated
 *  @{
 */

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
          int is_circular),
"Use the new API and vrna_mfe() instead");

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
DEPRECATED(float fold( const char *sequence, char *structure),
"Use vrna_fold() or vrna_mfe() instead");

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
DEPRECATED(float circfold( const char *sequence, char *structure),
"Use vrna_circfold() or vrna_mfe() instead");


/**
 *  @brief Free arrays for mfe folding
 *
 *  @deprecated See vrna_fold(), vrna_circfold(), or vrna_mfe() and #vrna_fold_compound_t for the usage of the new API!
 *
 */
DEPRECATED(void free_arrays(void),
"This function is obsolete");



/**
 *  @brief Recalculate energy parameters
 *
 *  @deprecated For non-default model settings use the new API with vrna_params_subst() and vrna_mfe() instead!
 *
 */
DEPRECATED(void update_fold_params(void),
"This function is obsolete");

/**
 *  @brief Recalculate energy parameters
 *
 *  @deprecated For non-default model settings use the new API with vrna_params_subst() and vrna_mfe() instead!
 *
 */
DEPRECATED(void update_fold_params_par(vrna_param_t *parameters),
"Use the new API with vrna_fold_compound_t datastructure instead");

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
                    char **ptype_p),
"Use the new API with vrna_fold_compound_t datastructure instead");

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
                        vrna_param_t **P_p),
"Use the new API with vrna_fold_compound_t datastructure instead");

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
                        char **ptype_p),
"Use the new API with vrna_fold_compound_t datastructure instead");

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
                            vrna_param_t **P_p),
"Use the new API with vrna_fold_compound_t datastructure instead");



/* finally moved the loop energy function declarations to this header...  */
/* BUT: The functions only exist for backward compatibility reasons!      */
/* You better include "loop_energies.h" and call the functions:           */
/* vrna_E_hairpin() and E_IntLoop() which are (almost) threadsafe as they get  */
/* a pointer to the energy parameter data structure as additional argument */

/**
 *  @deprecated {This function is deprecated and will be removed soon.
 *  Use @ref vrna_E_internal() instead!}
 */
DEPRECATED(int LoopEnergy(int n1,
                          int n2,
                          int type,
                          int type_2,
                          int si1,
                          int sj1,
                          int sp1,
                          int sq1),
"This function is obsolete");

/**
 *  @deprecated {This function is deprecated and will be removed soon.
 *  Use @ref vrna_E_hairpin() instead!}
 */
DEPRECATED(int HairpinE(int size,
                        int type,
                        int si1,
                        int sj1,
                        const char *string),
"Use vrna_E_hairpin() instead");

/**
 *  Allocate arrays for folding\n
 *  @deprecated See vrna_mfe() and #vrna_fold_compound_t for the usage of the new API!
 *
 */
DEPRECATED(void initialize_fold(int length),
"This function is obsolete");

/**
 *
 */
DEPRECATED(char *backtrack_fold_from_pair(char *sequence,
                                          int i,
                                          int j),
"This function is obsolete. Consider using vrna_backtrack_from_intervals() instead");


#endif

/**
 *  @}
 */

#endif
