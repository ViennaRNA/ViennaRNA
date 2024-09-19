#ifndef VIENNA_RNA_PACKAGE_LOOPS_EXTERNAL_H
#define VIENNA_RNA_PACKAGE_LOOPS_EXTERNAL_H

#include <ViennaRNA/datastructures/basic.h>
#include <ViennaRNA/datastructures/array.h>
#include <ViennaRNA/fold_compound.h>
#include <ViennaRNA/params/basic.h>

/* include energy evaluation API here for backward compatibility */
#include <ViennaRNA/eval/exterior.h>

/* include MFE API here for backward compatibility */
#include <ViennaRNA/mfe/exterior.h>

/* include backtrack API here for backward compatibility */
#include <ViennaRNA/backtrack/exterior.h>


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


#endif
