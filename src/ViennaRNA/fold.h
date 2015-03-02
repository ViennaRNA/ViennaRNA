#ifndef __VIENNA_RNA_PACKAGE_FOLD_H__
#define __VIENNA_RNA_PACKAGE_FOLD_H__

#include <ViennaRNA/data_structures.h>
#include <ViennaRNA/params.h>
#include <ViennaRNA/eval.h>

#ifdef __GNUC__
#define DEPRECATED(func) func __attribute__ ((deprecated))
#else
#define DEPRECATED(func) func
#endif

/**
 *  \addtogroup mfe_fold
 *  \ingroup folding_routines
 *  \brief This section covers all functions and variables related to the calculation
 *  of minimum free energy (MFE) structures.
 *
 *  The library provides a fast dynamic programming minimum free energy
 *  folding algorithm as described in \cite zuker:1981.
 *  All relevant parts that directly implement the "Zuker & Stiegler" algorithm for single
 *  sequences are described in this section.
 *
 *  Folding of circular RNA sequences is handled as a post-processing step of the forward
 *  recursions. See \cite hofacker:2006 for further details.
 *
 *  Nevertheless, the RNAlib also
 *  provides interfaces for the prediction of consensus MFE structures of sequence alignments,
 *  MFE structure for two hybridized sequences, local optimal structures and many more. For
 *  those more specialized variants of MFE folding routines, please consult the appropriate
 *  subsections (Modules) as listed above.
 *  
 *  \file fold.h
 *  \brief MFE calculations for single RNA sequences
 * 
 *  This file includes (almost) all function declarations within the RNAlib that are related to
 *  MFE folding...
 */

/**
 *  \addtogroup mfe_fold
 *  @{
 *    \brief This module contains all functions and variables related to the calculation
 *    of global minimum free energy structures for single sequences.
 *
 *    The library provides a fast dynamic programming minimum free energy
 *    folding algorithm as described by \ref zuker_81 "Zuker & Stiegler (1981)".
 *  @}
 */

/**
 *  \brief Compute minimum free energy and an appropriate secondary
 *  structure of an RNA sequence
 * 
 *  The first parameter given, the RNA sequence, must be \a uppercase and should only contain
 *  an alphabet \f$\Sigma\f$ that is understood by the RNAlib\n
 *  (e.g. \f$ \Sigma = \{A,U,C,G\} \f$)\n
 *
 *  The second parameter, \a structure, must always point to an allocated
 *  block of memory with a size of at least \f$\mathrm{strlen}(\mathrm{sequence})+1\f$
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
 *  that exhibits the minimum of free energy will be written to the memory \a structure is pointing to.
 *  The function returns the minimum of free energy for any fold of the sequence given.
 *
 *  \note OpenMP: Passing NULL to the 'parameters' argument involves access to several global model
 *        detail variables and thus is not to be considered threadsafe
 *
 *  \ingroup mfe_fold
 *
 *  \deprecated use vrna_fold() instead
 *
 *  \see vrna_fold(), fold(), circfold(), #vrna_md_t, set_energy_model(), get_scaled_parameters()
 *
 *  \param sequence       RNA sequence
 *  \param structure      A pointer to the character array where the
 *                        secondary structure in dot-bracket notation will be written to
 *  \param parameters     A data structure containing the prescaled energy contributions
 *                        and the model details. (NULL may be passed, see OpenMP notes above)
 *  \param is_constrained Switch to indicate that a structure contraint is passed via the structure argument (0==off)
 *  \param is_circular    Switch to (de-)activate postprocessing steps in case RNA sequence is circular (0==off)
 *
 *  \return the minimum free energy (MFE) in kcal/mol
 */
DEPRECATED(float 
fold_par( const char *sequence,
          char *structure,
          vrna_param_t *parameters,
          int is_constrained,
          int is_circular));

/**
 *  \brief Compute minimum free energy and an appropriate secondary
 *  structure of an RNA sequence
 * 
 *  The first parameter given, the RNA sequence, must be \a uppercase and should only contain
 *  an alphabet \f$\Sigma\f$ that is understood by the RNAlib\n
 *  (e.g. \f$ \Sigma = \{A,U,C,G\} \f$)\n
 *
 *  The second parameter, \a structure, must always point to an allocated
 *  block of memory with a size of at least \f$\mathrm{strlen}(\mathrm{sequence})+1\f$
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
 *  that exhibits the minimum of free energy will be written to the memory \a structure is pointing to.
 *  The function returns the minimum of free energy for any fold of the sequence given.
 *
 *  \note OpenMP: Passing NULL to the 'parameters' argument involves access to several global model
 *        detail variables and thus is not to be considered threadsafe
 *
 *  \ingroup mfe_fold
 *
 *  \see fold(), circfold(), #vrna_md_t, set_energy_model(), get_scaled_parameters()
 *
 *  \param vc             fold compound
 *  \param structure      A pointer to the character array where the
 *                        secondary structure in dot-bracket notation will be written to
 *
 *  \return the minimum free energy (MFE) in kcal/mol
 */
float
vrna_fold(vrna_fold_compound *vc,
          char *structure);

/**
 *  \brief Compute minimum free energy and an appropriate secondary structure of an RNA sequence
 *
 *  This function essentially does the same thing as fold_par(). However, it takes its model details,
 *  i.e. #temperature, #dangles, #tetra_loop, #noGU, #no_closingGU, #fold_constrained, #noLonelyPairs
 *  from the current global settings within the library
 *
 *  Use fold_par() for a completely threadsafe variant
 *
 *  \ingroup mfe_fold
 *
 *  \deprecated use vrna_fold() instead
 *
 *  \see fold_par(), circfold()
 *
 *  \param sequence RNA sequence
 *  \param structure A pointer to the character array where the
 *         secondary structure in dot-bracket notation will be written to
 *  \return the minimum free energy (MFE) in kcal/mol
 */
DEPRECATED(float fold( const char *sequence, char *structure));

/**
 *  \brief Compute minimum free energy and an appropriate secondary structure of a circular RNA sequence
 *
 *  This function essentially does the same thing as fold_par(). However, it takes its model details,
 *  i.e. #temperature, #dangles, #tetra_loop, #noGU, #no_closingGU, #fold_constrained, #noLonelyPairs
 *  from the current global settings within the library
 *
 *  Use fold_par() for a completely threadsafe variant
 *
 *  \ingroup mfe_fold
 *
 *  \see fold_par(), circfold()
 *
 *  \param sequence RNA sequence
 *  \param structure A pointer to the character array where the
 *         secondary structure in dot-bracket notation will be written to
 *  \return the minimum free energy (MFE) in kcal/mol
 */
DEPRECATED(float circfold( const char *sequence, char *structure));


/**
 *  \brief Free arrays for mfe folding
 *
 *  \ingroup mfe_fold
 *
 */
DEPRECATED(void free_arrays(void));



/**
 *  \brief Recalculate energy parameters
 *
 *  \deprecated use vrna_update_fold_params() instead
 *
 *  \ingroup mfe_fold
 */
DEPRECATED(void update_fold_params(void));

/**
 *  \brief Recalculate energy parameters
 *
 *  \deprecated use vrna_update_fold_params() instead
 *
 *  \ingroup mfe_fold
 * 
 */
DEPRECATED(void update_fold_params_par(vrna_param_t *parameters));

/**
 *
 *  \ingroup mfe_fold
 * 
 */
void
vrna_update_fold_params(vrna_fold_compound *vc,
                        vrna_param_t *parameters);

/**
 *
 *  \ingroup mfe_fold
 * 
 */
char  *backtrack_fold_from_pair(char *sequence,
                                int i,
                                int j);

plist *
vrna_backtrack_from_intervals(vrna_fold_compound *vc,
                              bondT *bp_stack,
                              sect bt_stack[],
                              int s);

/**
 *
 *  \ingroup mfe_fold
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
 *  \ingroup mfe_fold
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
 *  \ingroup mfe_fold
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
 *  \ingroup mfe_fold
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
/* a pointer to the energy parameter datastructure as additional argument */

/**
 *  \deprecated {This function is deprecated and will be removed soon.
 *  Use \ref E_IntLoop() instead!}
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
 *  \deprecated {This function is deprecated and will be removed soon.
 *  Use \ref E_Hairpin() instead!}
 */
DEPRECATED(int HairpinE(int size,
                        int type,
                        int si1,
                        int sj1,
                        const char *string));

/**
 *  Allocate arrays for folding\n
 *  \deprecated {This function is deprecated and will be removed soon!}
 * 
 */
DEPRECATED(void initialize_fold(int length));

#endif
