#ifndef VIENNA_RNA_PACKAGE_BT_HAIRPIN_H
#define VIENNA_RNA_PACKAGE_BT_HAIRPIN_H

#include <math.h>
#include <string.h>
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
 *
 *  @file     ViennaRNA/backtrack/hairpin.h
 *  @ingroup  eval, eval_loops, eval_loops_hp
 *  @brief    Backtracking of hairpin loops for MFE calculations
 */


/**
 *  @addtogroup mfe_backtracking
 *  @{
 */

/**
 *  @brief Backtrack a hairpin loop closed by @f$ (i,j) @f$
 *
 *  @note This function is polymorphic! The provided #vrna_fold_compound_t may be of type
 *        #VRNA_FC_TYPE_SINGLE or #VRNA_FC_TYPE_COMPARATIVE
 *
 */
int
vrna_bt_hairpin(vrna_fold_compound_t  *fc,
                unsigned int          i,
                unsigned int          j,
                int                   en,
                vrna_bps_t            bp_stack,
                vrna_bts_t            bt_stack);


/** @} */

#ifndef VRNA_DISABLE_BACKWARD_COMPATIBILITY

/**
 *  @addtogroup global_deprecated
 *  @{
 */
DEPRECATED(int
           vrna_BT_hp_loop(vrna_fold_compound_t *fc,
                           int                  i,
                           int                  j,
                           int                  en,
                           vrna_bp_stack_t      *bp_stack,
                           unsigned int         *stack_count),
           "Use vrna_bt_hairpin() instead!");


/** @} */

#endif

#endif
