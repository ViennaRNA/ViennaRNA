#ifndef VIENNA_RNA_PACKAGE_MFE_EXTERIOR_H
#define VIENNA_RNA_PACKAGE_MFE_EXTERIOR_H

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
 *  @file     ViennaRNA/mfe/external.h
 *  @ingroup  mfe
 *  @brief    Energy evaluation of exterior loops for MFE computations
 */

/**
 *  @addtogroup   eval_loops_ext
 *  @{
 */


/**
 *  @name Minimum Free Energy API
 *  @{
 */


int
vrna_E_ext_loop_5(vrna_fold_compound_t *fc);


int
vrna_E_ext_loop_3(vrna_fold_compound_t  *fc,
                  int                   i);


/* End basic interface */
/**@}*/


/**
 * @}
 */


/**
 *  @addtogroup mfe_backtracking
 *  @{
 */

unsigned int
vrna_bt_f(vrna_fold_compound_t *fc,
          unsigned int         i,
          unsigned int         j,
          vrna_bps_t           bp_stack,
          vrna_bts_t           bt_stack);


unsigned int
vrna_bt_ext_loop_f5(vrna_fold_compound_t  *fc,
                    unsigned int          j,
                    vrna_bps_t            bp_stack,
                    vrna_bts_t            bt_stack);


unsigned int
vrna_bt_ext_loop_f3(vrna_fold_compound_t  *fc,
                    unsigned int                   *k,
                    unsigned int                   maxdist,
                    unsigned int                   *i,
                    unsigned int                   *j,
                    vrna_bps_t            bp_stack);


int
vrna_bt_ext_loop_f3_pp(vrna_fold_compound_t *fc,
                       unsigned int                  *i,
                       unsigned int                  maxdist);


/**
 * @}
 */

#ifndef VRNA_DISABLE_BACKWARD_COMPATIBILITY

DEPRECATED(int
vrna_BT_ext_loop_f5(vrna_fold_compound_t  *fc,
                    int                   *k,
                    int                   *i,
                    int                   *j,
                    vrna_bp_stack_t       *bp_stack,
                    unsigned int          *stack_count),
            "Use vrna_bt_ext_loop_f5() or vrna_bt_f() instead!");


DEPRECATED(int
vrna_BT_ext_loop_f3(vrna_fold_compound_t  *fc,
                    int                   *k,
                    int                   maxdist,
                    int                   *i,
                    int                   *j,
                    vrna_bp_stack_t       *bp_stack,
                    unsigned int          *stack_count),
            "Use vrna_bt_ext_loop_f3() instead!");


DEPRECATED(int
vrna_BT_ext_loop_f3_pp(vrna_fold_compound_t *fc,
                       int                  *i,
                       int                  maxdist),
            "Use vrna_bt_ext_loop_f3_pp() instead!");



#endif

#endif
