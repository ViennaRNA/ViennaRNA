#ifndef VIENNA_RNA_PACKAGE_BT_EXTERIOR_H
#define VIENNA_RNA_PACKAGE_BT_EXTERIOR_H

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
 *  @file     ViennaRNA/backtrack/exterior.h
 *  @ingroup  mfe
 *  @brief    Backtracking of exterior loops for MFE computations
 */


/**
 *  @addtogroup mfe_backtracking
 *  @{
 */

unsigned int
vrna_bt_f(vrna_fold_compound_t  *fc,
          unsigned int          i,
          unsigned int          j,
          vrna_bps_t            bp_stack,
          vrna_bts_t            bt_stack);


unsigned int
vrna_bt_exterior_f5(vrna_fold_compound_t  *fc,
                    unsigned int          j,
                    vrna_bps_t            bp_stack,
                    vrna_bts_t            bt_stack);


unsigned int
vrna_bt_exterior_f3(vrna_fold_compound_t  *fc,
                    unsigned int          i,
                    unsigned int          j,
                    vrna_bps_t            bp_stack,
                    vrna_bts_t            bt_stack);


unsigned int
vrna_bt_exterior_f3_pp(vrna_fold_compound_t *fc,
                       unsigned int         *i,
                       unsigned int         maxdist);


/**
 * @}
 */

#ifndef VRNA_DISABLE_BACKWARD_COMPATIBILITY

/**
 *  @addtogroup global_deprecated
 *  @{
 */
DEPRECATED(int
           vrna_BT_ext_loop_f5(vrna_fold_compound_t *fc,
                               int                  *k,
                               int                  *i,
                               int                  *j,
                               vrna_bp_stack_t      *bp_stack,
                               unsigned int         *stack_count),
           "Use vrna_bt_exterior_f5() or vrna_bt_f() instead!");


DEPRECATED(int
           vrna_BT_ext_loop_f3(vrna_fold_compound_t *fc,
                               int                  *k,
                               int                  maxdist,
                               int                  *i,
                               int                  *j,
                               vrna_bp_stack_t      *bp_stack,
                               unsigned int         *stack_count),
           "Use vrna_bt_exterior_f3() instead!");


DEPRECATED(int
           vrna_BT_ext_loop_f3_pp(vrna_fold_compound_t  *fc,
                                  int                   *i,
                                  int                   maxdist),
           "Use vrna_bt_exterior_f3_pp() instead!");

/**
 * @}
 */


#endif

#endif
