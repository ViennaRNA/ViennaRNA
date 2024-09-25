#ifndef VIENNA_RNA_PACKAGE_BT_INTERNAL_H
#define VIENNA_RNA_PACKAGE_BT_INTERNAL_H

#include <ViennaRNA/utils/basic.h>
#include <ViennaRNA/datastructures/basic.h>
#include <ViennaRNA/fold_compound.h>


#ifdef VRNA_WARN_DEPRECATED
# if defined(DEPRECATED)
#   undef DEPRECATED
# endif
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
 *  @file     ViennaRNA/backtrack/internal.h
 *  @ingroup  mfe
 *  @brief    Backtracking of internal loops for MFE computations
 */

/**
 *  @addtogroup mfe_backtracking
 *  @{
 */


/**
 *  @brief Backtrack a stacked pair closed by @f$ (i,j) @f$
 *
 */
int
vrna_bt_stacked_pairs(vrna_fold_compound_t  *fc,
                      unsigned int          i,
                      unsigned int          j,
                      int                   *en,
                      vrna_bps_t            bp_stack,
                      vrna_bts_t            bt_stack);


/**
 *  @brief Backtrack an internal loop closed by @f$ (i,j) @f$
 *
 */
int
vrna_bt_internal_loop(vrna_fold_compound_t *fc,
                      unsigned int         i,
                      unsigned int         j,
                      int                  en,
                      vrna_bps_t           bp_stack,
                      vrna_bts_t           bt_stack);


/**
 * @}
 */


#ifndef VRNA_DISABLE_BACKWARD_COMPATIBILITY


/**
 *  @addtogroup global_deprecated
 *  @{
 */
DEPRECATED(int
           vrna_BT_stack(vrna_fold_compound_t *fc,
                         int                  *i,
                         int                  *j,
                         int                  *en,
                         vrna_bp_stack_t      *bp_stack,
                         unsigned int         *stack_count),
           "Use vrna_bt_stack() instead!");


DEPRECATED(int
           vrna_BT_int_loop(vrna_fold_compound_t  *fc,
                            int                   *i,
                            int                   *j,
                            int                   en,
                            vrna_bp_stack_t       *bp_stack,
                            unsigned int          *stack_count),
           "Use vrna_bt_internal_loop() instead!");

/**
 * @}
 */


#endif

#endif
