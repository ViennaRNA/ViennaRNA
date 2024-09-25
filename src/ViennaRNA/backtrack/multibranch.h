#ifndef VIENNA_RNA_PACKAGE_BT_MULTIBRANCH_H
#define VIENNA_RNA_PACKAGE_BT_MULTIBRANCH_H

#include <ViennaRNA/utils/basic.h>
#include <ViennaRNA/datastructures/basic.h>
#include <ViennaRNA/fold_compound.h>
#include <ViennaRNA/params/basic.h>

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
 *  @file     ViennaRNA/backtrack/multibranch.h
 *  @ingroup  mfe
 *  @brief    Backtracking of multibranch loops for MFE
 */


/**
 *  @addtogroup mfe_backtracking
 *  @{
 */

unsigned int
vrna_bt_m(vrna_fold_compound_t  *fc,
          unsigned int          i,
          unsigned int          j,
          vrna_bps_t            bp_stack,
          vrna_bts_t            bt_stack);


/**
 *  @brief  Backtrack the decomposition of a multi branch loop closed by @f$ (i,j) @f$
 *
 *  @param    fc          The #vrna_fold_compound_t filled with all relevant data for backtracking
 *  @param    i           5' position of base pair closing the loop (will be set to 5' position
 *                        of leftmost decomposed block upon successful backtracking)
 *  @param    j           3' position of base pair closing the loop (will be set to 3' position
 *                        of rightmost decomposed block upon successful backtracking)
 *  @param    k           Split position that delimits leftmost from rightmost block, [i,k] and
 *                        [k+1, j], respectively. (Will be set upon successful backtracking)
 *  @param    en          The energy contribution of the substructure enclosed by @f$ (i,j) @f$
 *  @param    component1  Type of leftmost block (1 = ML, 2 = C)
 *  @param    component2  Type of rightmost block (1 = ML, 2 = C)
 *  @returns              1, if backtracking succeeded, 0 otherwise.
 */
unsigned int
vrna_bt_multibranch_loop(vrna_fold_compound_t  *fc,
                         unsigned int          i,
                         unsigned int          j,
                         int                   en,
                         vrna_bps_t            bp_stack,
                         vrna_bts_t            bt_stack);


unsigned int
vrna_bt_multibranch_split(vrna_fold_compound_t  *fc,
                          unsigned int          i,
                          unsigned int          j,
                          vrna_bps_t            bp_stack,
                          vrna_bts_t            bt_stack);


/**
 * @}
 */


#ifndef VRNA_DISABLE_BACKWARD_COMPATIBILITY

/**
 *  @addtogroup global_deprecated
 *  @{
 */
DEPRECATED(int
           vrna_BT_mb_loop(vrna_fold_compound_t *fc,
                           int                  *i,
                           int                  *j,
                           int                  *k,
                           int                  en,
                           unsigned int         *component1,
                           unsigned int         *component2),
           "Use vrna_bt_multibranch_loop() instead!");


DEPRECATED(int
           vrna_BT_mb_loop_split(vrna_fold_compound_t *fc,
                                 int                  *i,
                                 int                  *j,
                                 int                  *k,
                                 int                  *l,
                                 unsigned int         *component1,
                                 unsigned int         *component2,
                                 vrna_bp_stack_t      *bp_stack,
                                 unsigned int         *stack_count),
           "Use vrna_bt_multibranch_split() instead!");

/**
 * @}
 */


#endif

#endif
