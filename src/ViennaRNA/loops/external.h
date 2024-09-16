#ifndef VIENNA_RNA_PACKAGE_LOOPS_EXTERNAL_H
#define VIENNA_RNA_PACKAGE_LOOPS_EXTERNAL_H

#include <ViennaRNA/datastructures/basic.h>
#include <ViennaRNA/datastructures/array.h>
#include <ViennaRNA/fold_compound.h>
#include <ViennaRNA/params/basic.h>

/* include energy evaluation API here for backward compatibility */
#include <ViennaRNA/eval/external.h>


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
 *  @file     ViennaRNA/loops/external.h
 *  @ingroup  eval, eval_loops, eval_loops_ext
 *  @brief    Energy evaluation of exterior loops for MFE and partition function calculations
 */

/**
 *  @addtogroup   eval_loops_ext
 *  @{
 */


/**
 *  @name Basic free energy interface
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
 *  @name Boltzmann weight (partition function) interface
 *  @{
 */


/**
 *  @brief  Auxiliary helper arrays for fast exterior loop computations
 *
 *  @see  vrna_exp_E_ext_fast_init(), vrna_exp_E_ext_fast_rotate(),
 *        vrna_exp_E_ext_fast_free(), vrna_exp_E_ext_fast()
 */
typedef struct vrna_mx_pf_aux_el_s *vrna_mx_pf_aux_el_t;


vrna_mx_pf_aux_el_t
vrna_exp_E_ext_fast_init(vrna_fold_compound_t *fc);


void
vrna_exp_E_ext_fast_rotate(vrna_mx_pf_aux_el_t aux_mx);


void
vrna_exp_E_ext_fast_free(vrna_mx_pf_aux_el_t aux_mx);


FLT_OR_DBL
vrna_exp_E_ext_fast(vrna_fold_compound_t  *fc,
                    int                   i,
                    int                   j,
                    vrna_mx_pf_aux_el_t   aux_mx);


void
vrna_exp_E_ext_fast_update(vrna_fold_compound_t *fc,
                           int                  j,
                           vrna_mx_pf_aux_el_t  aux_mx);


/* End partition function interface */
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
