#ifndef VIENNA_RNA_PACKAGE_MFE_H
#define VIENNA_RNA_PACKAGE_MFE_H

#include <stdio.h>
#include <ViennaRNA/datastructures/basic.h>
#include <ViennaRNA/fold_compound.h>

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
 *
 *  @file     ViennaRNA/mfe/global.h
 *  @ingroup  mfe, mfe_global
 *  @brief    Compute (global) Minimum Free energy (MFE)
 *
 *  This file includes (almost) all function declarations within the RNAlib that are related to
 *  MFE folding...
 */


/**
 *  @addtogroup  mfe_global
 *  @{
 */

/**
 *  @name Basic global MFE prediction interface
 *  @{
 */

/**
 *  @brief Compute minimum free energy and an appropriate secondary
 *  structure of an RNA sequence, or RNA sequence alignment
 *
 *  Depending on the type of the provided #vrna_fold_compound_t, this function
 *  predicts the MFE for a single sequence (or connected component of multiple
 *  sequences), or an averaged MFE for a sequence alignment.
 *  If backtracking is activated, it also constructs the corresponding secondary
 *  structure, or consensus structure.
 *  Therefore, the second parameter, @a structure, has to point to an allocated
 *  block of memory with a size of at least @f$\mathrm{strlen}(\mathrm{sequence})+1@f$ to
 *  store the backtracked MFE structure. (For consensus structures, this is the length of
 *  the alignment + 1. If @p NULL is passed, no backtracking will be performed.
 *
 *  @note This function is polymorphic. It accepts #vrna_fold_compound_t of type
 *        #VRNA_FC_TYPE_SINGLE, and #VRNA_FC_TYPE_COMPARATIVE.
 *
 *  @see #vrna_fold_compound_t, vrna_fold_compound(), vrna_fold(), vrna_circfold(),
 *        vrna_fold_compound_comparative(), vrna_alifold(), vrna_circalifold()
 *
 *  @param fc             fold compound
 *  @param structure      A pointer to the character array where the
 *                        secondary structure in dot-bracket notation will be written to (Maybe NULL)
 *
 *  @return the minimum free energy (MFE) in kcal/mol
 */
float
vrna_mfe(vrna_fold_compound_t *fc,
         char                 *structure);


/**
 *  @brief Compute the minimum free energy of two interacting RNA molecules
 *
 *  The code is analog to the vrna_mfe() function.
 *
 *  @deprecated This function is obsolete since vrna_mfe() can handle complexes multiple
 *              sequences since v2.5.0. Use vrna_mfe() for connected component MFE instead
 *              and compute MFEs of unconnected states separately.
 *
 *  @see  vrna_mfe()
 *
 *  @param    fc  fold compound
 *  @param    structure Will hold the barcket dot structure of the dimer molecule
 *  @return   minimum free energy of the structure
 */
DEPRECATED(float
           vrna_mfe_dimer(vrna_fold_compound_t  *fc,
                          char                  *structure),
           "Use vrna_mfe() instead");


/**
 * End basic MFE interface
 * @}
 */


/**
 *  @name Simplified global MFE prediction using sequence(s) or multiple sequence alignment(s)
 *  @{
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
 *        you require other conditions than specified by the default model details, use vrna_mfe(),
 *        and the data structure #vrna_fold_compound_t instead.
 *
 *  @see vrna_circfold(), vrna_mfe()
 *
 *  @param sequence   RNA sequence
 *  @param structure  A pointer to the character array where the
 *         secondary structure in dot-bracket notation will be written to
 *  @return the minimum free energy (MFE) in kcal/mol
 */
float
vrna_fold(const char  *sequence,
          char        *structure);


/**
 *  @brief Compute Minimum Free Energy (MFE), and a corresponding secondary structure for a circular RNA sequence
 *
 *  This simplified interface to vrna_mfe() computes the MFE and, if required, a secondary structure for a
 *  circular RNA sequence using default options. Memory required for dynamic programming (DP) matrices will
 *  be allocated and free'd on-the-fly. Hence, after return of this function, the recursively filled
 *  matrices are not available any more for any post-processing, e.g. suboptimal backtracking, etc.
 *
 *  Folding of circular RNA sequences is handled as a post-processing step of the forward
 *  recursions. See @rstinline :cite:t:`hofacker:2006` @endrst for further details.
 *
 *  @note In case you want to use the filled DP matrices for any subsequent post-processing step, or
 *        you require other conditions than specified by the default model details, use vrna_mfe(),
 *        and the data structure #vrna_fold_compound_t instead.
 *
 *  @see vrna_fold(), vrna_mfe()
 *
 *  @param sequence   RNA sequence
 *  @param structure  A pointer to the character array where the
 *         secondary structure in dot-bracket notation will be written to
 *  @return the minimum free energy (MFE) in kcal/mol
 */
float
vrna_circfold(const char  *sequence,
              char        *structure);


/**
 *  @brief  Compute Minimum Free Energy (MFE), and a corresponding consensus secondary structure
 *          for an RNA sequence alignment using a comparative method
 *
 *  This simplified interface to vrna_mfe() computes the MFE and, if required, a consensus secondary
 *  structure for an RNA sequence alignment using default options. Memory required for dynamic programming
 *  (DP) matrices will be allocated and free'd on-the-fly. Hence, after return of this function, the
 *  recursively filled matrices are not available any more for any post-processing, e.g. suboptimal
 *  backtracking, etc.
 *
 *  @note In case you want to use the filled DP matrices for any subsequent post-processing step, or
 *        you require other conditions than specified by the default model details, use vrna_mfe(),
 *        and the data structure #vrna_fold_compound_t instead.
 *
 *  @see vrna_circalifold(), vrna_mfe()
 *
 *  @param sequences  RNA sequence alignment
 *  @param structure  A pointer to the character array where the
 *         secondary structure in dot-bracket notation will be written to
 *  @return the minimum free energy (MFE) in kcal/mol
 */
float
vrna_alifold(const char **sequences,
             char       *structure);


/**
 *  @brief  Compute Minimum Free Energy (MFE), and a corresponding consensus secondary structure
 *          for a sequence alignment of circular RNAs using a comparative method
 *
 *  This simplified interface to vrna_mfe() computes the MFE and, if required, a consensus secondary
 *  structure for an RNA sequence alignment using default options. Memory required for dynamic programming
 *  (DP) matrices will be allocated and free'd on-the-fly. Hence, after return of this function, the
 *  recursively filled matrices are not available any more for any post-processing, e.g. suboptimal
 *  backtracking, etc.
 *
 *  Folding of circular RNA sequences is handled as a post-processing step of the forward
 *  recursions. See @rstinline :cite:t:`hofacker:2006` @endrst for further details.
 *
 *  @note In case you want to use the filled DP matrices for any subsequent post-processing step, or
 *        you require other conditions than specified by the default model details, use vrna_mfe(),
 *        and the data structure #vrna_fold_compound_t instead.
 *
 *  @see vrna_alifold(), vrna_mfe()
 *
 *  @param sequences  Sequence alignment of circular RNAs
 *  @param structure  A pointer to the character array where the
 *         secondary structure in dot-bracket notation will be written to
 *  @return the minimum free energy (MFE) in kcal/mol
 */
float
vrna_circalifold(const char **sequences,
                 char       *structure);


/**
 *  @brief Compute Minimum Free Energy (MFE), and a corresponding secondary structure for two dimerized RNA sequences
 *
 *  This simplified interface to vrna_mfe() computes the MFE and, if required, a secondary structure for
 *  two RNA sequences upon dimerization using default options. Memory required for dynamic programming
 *  (DP) matrices will be allocated and free'd on-the-fly. Hence, after return of this function, the
 *  recursively filled matrices are not available any more for any post-processing, e.g. suboptimal
 *  backtracking, etc.
 *
 *  @note In case you want to use the filled DP matrices for any subsequent post-processing step, or
 *        you require other conditions than specified by the default model details, use vrna_mfe(),
 *        and the data structure #vrna_fold_compound_t instead.
 *
 *  @deprecated This function is obsolete since vrna_mfe()/vrna_fold() can handle complexes multiple
 *              sequences since v2.5.0. Use vrna_mfe()/vrna_fold() for connected component MFE instead
 *              and compute MFEs of unconnected states separately.
 *
 *  @see  vrna_fold(), vrna_mfe(), vrna_fold_compound(), #vrna_fold_compound_t, vrna_cut_point_insert()
 *
 *  @param sequence   two RNA sequences separated by the '&' character
 *  @param structure  A pointer to the character array where the
 *         secondary structure in dot-bracket notation will be written to
 *  @return the minimum free energy (MFE) in kcal/mol
 */
DEPRECATED(float
           vrna_cofold(const char *sequence,
                       char       *structure),
           "USe vrna_fold() instead");


/**
 * End simplified global MFE interface
 * @}
 */

/**
 * End group mfe_global
 * @}
 */


#endif
