#ifndef VIENNA_RNA_PACKAGE_COFOLD_H
#define VIENNA_RNA_PACKAGE_COFOLD_H

#include <ViennaRNA/data_structures.h>
#include <ViennaRNA/params.h>

#ifdef __GNUC__
#define DEPRECATED(func) func __attribute__ ((deprecated))
#else
#define DEPRECATED(func) func
#endif

/**
 *  @addtogroup cofold
 *  @brief Predict structures formed by two molecules upon hybridization.
 *
 *  The function of an RNA molecule often depends on its interaction with
 *  other RNAs. The following routines therefore allow to predict structures
 *  formed by two RNA molecules upon hybridization.\n
 *  One approach to co-folding two RNAs consists of concatenating the two
 *  sequences and keeping track of the concatenation point in all energy
 *  evaluations. Correspondingly, many of the cofold() and
 *  co_pf_fold() routines below take one sequence string as argument
 *  and use the the global variable #cut_point to mark the concatenation
 *  point. Note that while the <i>RNAcofold</i> program uses the '&' character
 *  to mark the chain break in its input, you should not use an '&' when using
 *  the library routines (set #cut_point instead).\n
 *  In a second approach to co-folding two RNAs, cofolding is seen as a
 *  stepwise process. In the first step the probability of an unpaired region
 *  is calculated and in a second step this probability of an unpaired region
 *  is  multiplied with the probability of an interaction between the two RNAs.
 *  This approach is implemented for the interaction between a long
 *  target sequence and a short ligand RNA. Function pf_unstru() calculates
 *  the partition function over all unpaired regions in the input
 *  sequence. Function pf_interact(), which calculates the
 *  partition function over all possible interactions between two
 *  sequences, needs both sequence as separate strings as input.
 *
 */

/**
 *  @addtogroup mfe_cofold
 *  @{
 *  @file cofold.h
 *
 *  @brief MFE version of cofolding routines
 *
 *  This file includes (almost) all function declarations within the <b>RNAlib</b> that are related to
 *  MFE Cofolding...
 *  This also includes the Zuker suboptimals calculations, since they are implemented using the cofold
 *  routines.
 */

/**
 *  @brief Compute the minimum free energy of two interacting RNA molecules
 *
 *  The code is analog to the vrna_mfe() function.
 *
 *  @ingroup mfe_cofold
 *
 *  @param    vc  fold compound
 *  @param    structure Will hold the barcket dot structure of the dimer molecule
 *  @return   minimum free energy of the structure
 */
float vrna_mfe_dimer( vrna_fold_compound_t *vc,
                      char *structure);

/**
 *  @brief Add a separating '&' character into a string according to cut-point position
 *
 *  If the cut-point position is less or equal to zero, this function just
 *  returns a copy of the provided string. Otherwise, the cut-point character
 *  is set at the corresponding position
 *
 *  @param  string    The original string
 *  @param  cp        The cut-point position
 *  @return           A copy of the provided string including the cut-point character
 */
char *vrna_cut_point_insert(const char *string,
                            int cp);

/**
 *  @brief  Remove a separating '&' character from a string
 *
 *  This function removes the cut-point indicating '&' character from a string
 *  and memorizes its position in a provided integer variable. If not '&' is
 *  found in the input, the integer variable is set to -1. The function returns
 *  a copy of the input string with the '&' being sliced out.
 *
 *  @param  string  The original string
 *  @param  cp      The cut-point position
 *  @return         A copy of the input string with the '&' being sliced out
 */
char *vrna_cut_point_remove(const char *string,
                            int *cp);

/**
 *  @brief Compute the minimum free energy of two interacting RNA molecules
 *
 *  The code is analog to the fold() function. If #cut_point ==-1 results
 *  should be the same as with fold().
 *
 *  @ingroup mfe_cofold
 *
 *  @deprecated use vrna_mfe_dimer() instead
 *
 *  @param    sequence  The two sequences concatenated
 *  @param    structure Will hold the barcket dot structure of the dimer molecule
 *  @return   minimum free energy of the structure
 */
DEPRECATED(float
cofold( const char *sequence,
        char *structure));

/**
 *  @brief Compute the minimum free energy of two interacting RNA molecules
 *
 *  @deprecated use vrna_mfe_dimer() instead
 *
 */
DEPRECATED(float
cofold_par( const char *string,
            char *structure,
            vrna_param_t *parameters,
            int is_constrained));

/**
 *  @brief Free memory occupied by cofold()
 *
 *  @deprecated This function will only free memory allocated by a prior call of cofold() or cofold_par().
 *  See vrna_mfe_dimer() for how to use the new API
 *
 *  @note folding matrices now reside in the fold compound, and should be free'd there
 *  @see  vrna_fc_destroy(), vrna_mfe_dimer()
 */
DEPRECATED(void free_co_arrays(void));

/**
 *  @brief Recalculate parameters
 *  @deprecated See vrna_params_update() for an alternative using the new API
 */
DEPRECATED(void update_cofold_params(void));

/**
 *  @brief Recalculate parameters
 *  @deprecated See vrna_params_update() for an alternative using the new API
 */
DEPRECATED(void update_cofold_params_par(vrna_param_t *parameters));


/**
 *  @brief Export the arrays of partition function cofold (with gquadruplex support)
 *
 *  Export the cofold arrays for use e.g. in the concentration
 *  Computations or suboptimal secondary structure backtracking
 *
 *  @deprecated folding matrices now reside within the fold compound. Thus, this function will
 *  only work in conjunction with a prior call to cofold() or cofold_par()
 *
 *  @see vrna_mfe_dimer() for the new API
 *
 *  @param  f5_p    A pointer to the 'f5' array, i.e. array conatining best free energy in interval [1,j]
 *  @param  c_p     A pointer to the 'c' array, i.e. array containing best free energy in interval [i,j] given that i pairs with j
 *  @param  fML_p   A pointer to the 'M' array, i.e. array containing best free energy in interval [i,j] for any multiloop segment with at least one stem
 *  @param  fM1_p   A pointer to the 'M1' array, i.e. array containing best free energy in interval [i,j] for multiloop segment with exactly one stem
 *  @param  fc_p    A pointer to the 'fc' array, i.e. array ...
 *  @param  ggg_p   A pointer to the 'ggg' array, i.e. array containing best free energy of a gquadruplex delimited by [i,j]
 *  @param  indx_p  A pointer to the indexing array used for accessing the energy matrices
 *  @param  ptype_p A pointer to the ptype array containing the base pair types for each possibility (i,j)
 */
DEPRECATED(void export_cofold_arrays_gq(int **f5_p,
                                        int **c_p,
                                        int **fML_p,
                                        int **fM1_p,
                                        int **fc_p,
                                        int **ggg_p,
                                        int **indx_p,
                                        char **ptype_p));

/**
 *  @brief Export the arrays of partition function cofold
 *
 *  Export the cofold arrays for use e.g. in the concentration
 *  Computations or suboptimal secondary structure backtracking
 *
 *  @deprecated folding matrices now reside within the #vrna_fold_compound_t. Thus, this function will
 *  only work in conjunction with a prior call to the deprecated functions cofold() or cofold_par()
 *
 *  @see vrna_mfe_dimer() for the new API
 *
 *  @param  f5_p    A pointer to the 'f5' array, i.e. array conatining best free energy in interval [1,j]
 *  @param  c_p     A pointer to the 'c' array, i.e. array containing best free energy in interval [i,j] given that i pairs with j
 *  @param  fML_p   A pointer to the 'M' array, i.e. array containing best free energy in interval [i,j] for any multiloop segment with at least one stem
 *  @param  fM1_p   A pointer to the 'M1' array, i.e. array containing best free energy in interval [i,j] for multiloop segment with exactly one stem
 *  @param  fc_p    A pointer to the 'fc' array, i.e. array ...
 *  @param  indx_p  A pointer to the indexing array used for accessing the energy matrices
 *  @param  ptype_p A pointer to the ptype array containing the base pair types for each possibility (i,j)
 */
DEPRECATED(void export_cofold_arrays( int **f5_p,
                                      int **c_p,
                                      int **fML_p,
                                      int **fM1_p,
                                      int **fc_p,
                                      int **indx_p,
                                      char **ptype_p));


/**
 *  @}
 */

/**
 *  @brief Compute Zuker type suboptimal structures
 *
 *  Compute Suboptimal structures according to M. Zuker, i.e. for every
 *  possible base pair the minimum energy structure containing the resp. base pair.
 *  Returns a list of these structures and their energies.
 *
 *  @ingroup subopt_zuker
 *
 *  @deprecated use vrna_zukersubopt() instead
 *
 *  @param  string  RNA sequence
 *  @return         List of zuker suboptimal structures
 */
DEPRECATED(SOLUTION  *zukersubopt(const char *string));

/**
 *  @brief Compute Zuker type suboptimal structures
 *
 *  @ingroup subopt_zuker
 *
 *  @deprecated use vrna_zukersubopt() instead
 *
 */
DEPRECATED(SOLUTION  *zukersubopt_par(const char *string, vrna_param_t *parameters));

/**
 *  @brief Compute Zuker type suboptimal structures
 *
 *  Compute Suboptimal structures according to M. Zuker, i.e. for every
 *  possible base pair the minimum energy structure containing the resp. base pair.
 *  Returns a list of these structures and their energies.
 *
 *  @ingroup subopt_zuker
 *
 *  @param  vc  fold compound
 *  @return     List of zuker suboptimal structures
 */
SOLUTION *vrna_zukersubopt(vrna_fold_compound_t *vc);

/**
 *  @brief get_monomer_free_energies
 *
 *  Export monomer free energies out of cofold arrays
 *  @deprecated{This function is obsolete and will be removed soon!}
 *
 *  @param e1 A pointer to a variable where the energy of molecule A will be written to
 *  @param e2 A pointer to a variable where the energy of molecule B will be written to
 */
DEPRECATED(void get_monomere_mfes( float *e1, float *e2));


/**
 *  allocate arrays for folding
 *  @deprecated{This function is obsolete and will be removed soon!}
 */
DEPRECATED(void initialize_cofold(int length));

#endif
