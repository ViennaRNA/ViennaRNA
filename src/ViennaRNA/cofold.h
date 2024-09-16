#ifndef VIENNA_RNA_PACKAGE_COFOLD_H
#define VIENNA_RNA_PACKAGE_COFOLD_H

#include <ViennaRNA/datastructures/basic.h>
#include <ViennaRNA/params/basic.h>
#include <ViennaRNA/mfe/global.h>

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
 *  @file     cofold.h
 *  @ingroup mfe_global_deprecated
 *  @brief    MFE implementations for RNA-RNA interaction
 */

#ifndef VRNA_DISABLE_BACKWARD_COMPATIBILITY

/**
 *  @brief Compute the minimum free energy of two interacting RNA molecules
 *
 *  The code is analog to the fold() function. If #cut_point ==-1 results
 *  should be the same as with fold().
 *
 *  @ingroup mfe_global_deprecated
 *
 *  @deprecated use vrna_mfe_dimer() instead
 *
 *  @param    sequence  The two sequences concatenated
 *  @param    structure Will hold the barcket dot structure of the dimer molecule
 *  @return   minimum free energy of the structure
 */
DEPRECATED(float
           cofold(const char  *sequence,
                  char        *structure),
           "Use vrna_cofold() instead");

/**
 *  @brief Compute the minimum free energy of two interacting RNA molecules
 *
 *  @deprecated use vrna_mfe_dimer() instead
 *
 *  @ingroup mfe_global_deprecated
 */
DEPRECATED(float
           cofold_par(const char    *string,
                      char          *structure,
                      vrna_param_t  *parameters,
                      int           is_constrained),
           "Use the new API and vrna_mfe_dimer() instead");

/**
 *  @brief Free memory occupied by cofold()
 *
 *  @deprecated This function will only free memory allocated by a prior call of cofold() or cofold_par().
 *  See vrna_mfe_dimer() for how to use the new API
 *
 *  @note folding matrices now reside in the fold compound, and should be free'd there
 *
 *  @see  vrna_fc_destroy(), vrna_mfe_dimer()
 *
 *  @ingroup mfe_global_deprecated
 */
DEPRECATED(void
           free_co_arrays(void),
           "This function is obsolete");

/**
 *  @brief Recalculate parameters
 *  @deprecated See vrna_params_subst() for an alternative using the new API
 *
 *  @ingroup mfe_global_deprecated
 */
DEPRECATED(void
           update_cofold_params(void),
           "This function is obsolete");

/**
 *  @brief Recalculate parameters
 *  @deprecated See vrna_params_subst() for an alternative using the new API
 *
 *  @ingroup mfe_global_deprecated
 */
DEPRECATED(void
           update_cofold_params_par(vrna_param_t *parameters),
           "Use the new API with vrna_fold_compound_t instead");


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
 *  @ingroup mfe_global_deprecated
 *  @param  f5_p    A pointer to the 'f5' array, i.e. array conatining best free energy in interval [1,j]
 *  @param  c_p     A pointer to the 'c' array, i.e. array containing best free energy in interval [i,j] given that i pairs with j
 *  @param  fML_p   A pointer to the 'M' array, i.e. array containing best free energy in interval [i,j] for any multiloop segment with at least one stem
 *  @param  fM1_p   A pointer to the 'M1' array, i.e. array containing best free energy in interval [i,j] for multiloop segment with exactly one stem
 *  @param  fc_p    A pointer to the 'fc' array, i.e. array ...
 *  @param  ggg_p   A pointer to the 'ggg' array, i.e. array containing best free energy of a gquadruplex delimited by [i,j]
 *  @param  indx_p  A pointer to the indexing array used for accessing the energy matrices
 *  @param  ptype_p A pointer to the ptype array containing the base pair types for each possibility (i,j)
 */
DEPRECATED(void
           export_cofold_arrays_gq(int  **f5_p,
                                   int  **c_p,
                                   int  **fML_p,
                                   int  **fM1_p,
                                   int  **fc_p,
                                   int  **ggg_p,
                                   int  **indx_p,
                                   char **ptype_p),
           "Use the new API with vrna_fold_compound_t instead");

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
 *  @ingroup mfe_global_deprecated
 *  @param  f5_p    A pointer to the 'f5' array, i.e. array conatining best free energy in interval [1,j]
 *  @param  c_p     A pointer to the 'c' array, i.e. array containing best free energy in interval [i,j] given that i pairs with j
 *  @param  fML_p   A pointer to the 'M' array, i.e. array containing best free energy in interval [i,j] for any multiloop segment with at least one stem
 *  @param  fM1_p   A pointer to the 'M1' array, i.e. array containing best free energy in interval [i,j] for multiloop segment with exactly one stem
 *  @param  fc_p    A pointer to the 'fc' array, i.e. array ...
 *  @param  indx_p  A pointer to the indexing array used for accessing the energy matrices
 *  @param  ptype_p A pointer to the ptype array containing the base pair types for each possibility (i,j)
 */
DEPRECATED(void
           export_cofold_arrays(int   **f5_p,
                                int   **c_p,
                                int   **fML_p,
                                int   **fM1_p,
                                int   **fc_p,
                                int   **indx_p,
                                char  **ptype_p),
           "Use the new API with vrna_fold_compound_t instead");


/**
 *  allocate arrays for folding
 *  @deprecated{This function is obsolete and will be removed soon!}
 *
 *  @ingroup mfe_global_deprecated
 */
DEPRECATED(void
           initialize_cofold(int length),
           "This function is obsolete");

#endif

#endif
