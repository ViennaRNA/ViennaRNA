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
 *  @file     ViennaRNA/mfe/exterior.h
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
vrna_mfe_exterior_f5(vrna_fold_compound_t *fc);


int
vrna_mfe_exterior_f3(vrna_fold_compound_t *fc,
                     unsigned int         i);


/* End basic interface */
/**@}*/

#ifndef VRNA_DISABLE_BACKWARD_COMPATIBILITY

DEPRECATED(int
           vrna_E_ext_loop_5(vrna_fold_compound_t *fc),
           "Use vrna_mfe_exterior_f5() instead");


DEPRECATED(int
           vrna_E_ext_loop_3(vrna_fold_compound_t *fc),
           "Use vrna_mfe_exterior_f3() instead");
#endif

/**
 * @}
 */

#endif
