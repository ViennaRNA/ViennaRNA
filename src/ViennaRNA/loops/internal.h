#ifndef VIENNA_RNA_PACKAGE_LOOPS_INTERNAL_H
#define VIENNA_RNA_PACKAGE_LOOPS_INTERNAL_H

#include <math.h>

#include <ViennaRNA/utils/basic.h>
#include <ViennaRNA/params/default.h>
#include <ViennaRNA/datastructures/basic.h>
#include <ViennaRNA/fold_compound.h>
#include <ViennaRNA/params/basic.h>
#include <ViennaRNA/constraints/hard.h>
#include <ViennaRNA/constraints/soft.h>
#include <ViennaRNA/params/salt.h>

/* include energy evaluation API for backward compatibility */
#include <ViennaRNA/eval/internal.h>


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

#ifdef __GNUC__
# define INLINE inline
#else
# define INLINE
#endif

/**
 *  @file     ViennaRNA/loops/internal.h
 *  @ingroup  eval, eval_loops, eval_loops_int
 *  @brief    Energy evaluation of interior loops for MFE and partition function calculations
 */

/**
 *  @addtogroup   eval_loops_int
 *  @{
 */


/**
 *  @name Basic free energy interface
 *  @{
 */

int
vrna_E_int_loop(vrna_fold_compound_t  *fc,
                int                   i,
                int                   j);


int
vrna_E_ext_int_loop(vrna_fold_compound_t  *fc,
                    int                   i,
                    int                   j,
                    int                   *ip,
                    int                   *iq);


int
vrna_E_stack(vrna_fold_compound_t *fc,
             int                  i,
             int                  j);


/* End basic interface */
/**@}*/


/**
 *  @name Boltzmann weight (partition function) interface
 *  @{
 */


/* j < i indicates circular folding, i.e. collect contributions for exterior int loops */
FLT_OR_DBL
vrna_exp_E_int_loop(vrna_fold_compound_t  *fc,
                    int                   i,
                    int                   j);


/* End partition function interface */
/**@}*/

/**
 * @}
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
              unsigned int                  i,
              unsigned int                  j,
              int                   *en,
              vrna_bps_t            bp_stack,
              vrna_bts_t            bt_stack);


/**
 *  @brief Backtrack an interior loop closed by @f$ (i,j) @f$
 *
 */
int
vrna_bt_int_loop(vrna_fold_compound_t *fc,
                 unsigned int         i,
                 unsigned int         j,
                 int                  en,
                 vrna_bps_t           bp_stack,
                 vrna_bts_t           bt_stack);


/**
 * @}
 */


#ifndef VRNA_DISABLE_BACKWARD_COMPATIBILITY


DEPRECATED(int
vrna_BT_stack(vrna_fold_compound_t  *fc,
              int                   *i,
              int                   *j,
              int                   *en,
              vrna_bp_stack_t       *bp_stack,
              unsigned int          *stack_count),
          "Use vrna_bt_stack() instead!");


DEPRECATED(int
vrna_BT_int_loop(vrna_fold_compound_t *fc,
                 int                  *i,
                 int                  *j,
                 int                  en,
                 vrna_bp_stack_t      *bp_stack,
                 unsigned int         *stack_count),
          "Use vrna_bt_int_loop() instead!");


#endif

#endif
