#ifndef VIENNA_RNA_PACKAGE_MFE_H
#define VIENNA_RNA_PACKAGE_MFE_H

#include <stdio.h>
#include <ViennaRNA/datastructures/basic.h>
#include <ViennaRNA/fold_compound.h>

/**
 *
 *  @file mfe.h
 *  @ingroup  mfe, mfe_global
 *  @brief Compute Minimum Free energy (MFE) and backtrace corresponding secondary
 *         structures from RNA sequence data.
 *
 *  This file includes (almost) all function declarations within the RNAlib that are related to
 *  MFE folding...
 */

/**
 *  @addtogroup mfe
 *  @{
 *  @brief  Predicting the Minimum Free Energy (MFE) and a corresponding (consensus) secondary structure
 *
 *  In a nutshell we provide two different flavors for MFE prediction:
 *  * @ref mfe_global - to compute the MFE for the entire sequence
 *  * @ref mfe_window - to compute MFEs for each window using a sliding window approach
 *
 *  Each of these flavors, again, provides two implementations to either compute the MFE based on
 *  *  single RNA (DNA) sequence(s), or
 *  *  a comparative approach using multiple sequence alignments (MSA).
 *
 *  For the latter, a consensus secondary structure is predicted and our implementations compute
 *  an average of free energies for each sequence in the MSA plus an additional covariance
 *  pseudo-energy term.
 *
 *  The implementations for @ref mfe_backtracking are generally agnostic with respect to whether
 *  local or global structure prediction is in place.
 *  @}
 */


/**
 *  @addtogroup  mfe_global
 *  @{
 *  @brief  Variations of the global Minimum Free Energy (MFE) prediction algorithm
 *
 *  We provide implementations of the global MFE prediction algorithm for
 *  * Single sequences,
 *  * Multiple sequence alignments (MSA), and
 *  * RNA-RNA hybrids
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
 *  predicts the MFE for a single sequence, or a corresponding averaged MFE for
 *  a sequence alignment. If backtracking is activated, it also constructs the
 *  corresponding secondary structure, or consensus structure.
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
 *  @param vc             fold compound
 *  @param structure      A pointer to the character array where the
 *                        secondary structure in dot-bracket notation will be written to (Maybe NULL)
 *
 *  @return the minimum free energy (MFE) in kcal/mol
 */
float
vrna_mfe(vrna_fold_compound_t *vc,
         char                 *structure);


/**
 *  @brief Compute the minimum free energy of two interacting RNA molecules
 *
 *  The code is analog to the vrna_mfe() function.
 *
 *  @param    vc  fold compound
 *  @param    structure Will hold the barcket dot structure of the dimer molecule
 *  @return   minimum free energy of the structure
 */
float
vrna_mfe_dimer(vrna_fold_compound_t *vc,
               char                 *structure);


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
 *  you require other conditions than specified by the default model details, use vrna_mfe(),
 *  and the data structure #vrna_fold_compound_t instead.
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
 *  recursions. See @cite hofacker:2006 for further details.
 *
 *  @note In case you want to use the filled DP matrices for any subsequent post-processing step, or
 *  you require other conditions than specified by the default model details, use vrna_mfe(),
 *  and the data structure #vrna_fold_compound_t instead.
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
 *  you require other conditions than specified by the default model details, use vrna_mfe(),
 *  and the data structure #vrna_fold_compound_t instead.
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
 *  recursions. See @cite hofacker:2006 for further details.
 *
 *  @note In case you want to use the filled DP matrices for any subsequent post-processing step, or
 *  you require other conditions than specified by the default model details, use vrna_mfe(),
 *  and the data structure #vrna_fold_compound_t instead.
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
 *  you require other conditions than specified by the default model details, use vrna_mfe(),
 *  and the data structure #vrna_fold_compound_t instead.
 *
 *  @see vrna_mfe_dimer(), vrna_fold_compound(), #vrna_fold_compound_t, vrna_cut_point_insert()
 *
 *  @param sequence   two RNA sequences separated by the '&' character
 *  @param structure  A pointer to the character array where the
 *         secondary structure in dot-bracket notation will be written to
 *  @return the minimum free energy (MFE) in kcal/mol
 */
float
vrna_cofold(const char  *sequence,
            char        *structure);


/**
 * End simplified global MFE interface
 * @}
 */

/**
 * End group mfe_global
 * @}
 */

/**
 *  @addtogroup mfe_backtracking
 *  @{
 *  @brief   Backtracking related interfaces
 */

/**
 *  @brief
 */
int
vrna_backtrack_from_intervals(vrna_fold_compound_t  *vc,
                              vrna_bp_stack_t       *bp_stack,
                              sect                  bt_stack[],
                              int                   s);


/**
 *  @brief Backtrack an MFE (sub)structure
 *
 *  This function allows one to backtrack the MFE structure for a (sub)sequence
 *
 *  @note On error, the function returns #INF / 100. and stores the empty string
 *        in @p structure.
 *
 *  @pre  Requires pre-filled MFE dynamic programming matrices, i.e. one has to call vrna_mfe()
 *        prior to calling this function
 *
 *  @see vrna_mfe(), vrna_pbacktrack5()
 *
 *  @param fc             fold compound
 *  @param length         The length of the subsequence, starting from the 5' end
 *  @param structure      A pointer to the character array where the secondary structure in
 *                        dot-bracket notation will be written to. (Must have size of at least $p length + 1)
 *
 *  @return               The minimum free energy (MFE) for the specified @p length in kcal/mol and
 *                        a corresponding secondary structure in dot-bracket notation (stored in @p structure)
 */
float
vrna_backtrack5(vrna_fold_compound_t  *fc,
                unsigned int          length,
                char                  *structure);

int
vrna_backtrack_window(vrna_fold_compound_t  *fc,
                      const char            *Lfold_filename,
                      long                  file_pos,
                      char                  **structure,
                      double                mfe);

/**
 * End backtracking related interfaces
 * @}
 */


#endif
