#ifndef VIENNA_RNA_PACKAGE_BT_H
#define VIENNA_RNA_PACKAGE_BT_H

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
 *  @file     ViennaRNA/backtrack/global.h
 *  @ingroup  mfe, mfe_global
 *  @brief    Backtracking for (global) Minimum Free energy (MFE)
 *
 */

/**
 *  @addtogroup mfe_backtracking
 *  @{
 */

/**
 *  @brief  Backtrack a secondary structure with pre-evaluated structure components
 */
int
vrna_backtrack_from_intervals(vrna_fold_compound_t  *fc,
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


/**
 * End backtracking related interfaces
 * @}
 */


#endif
