#ifndef VIENNA_RNA_PACKAGE_PARTFUNC_INTERNAL_H
#define VIENNA_RNA_PACKAGE_PARTFUNC_INTERNAL_H

#include <math.h>

#include <ViennaRNA/utils/basic.h>
#include <ViennaRNA/params/default.h>
#include <ViennaRNA/datastructures/basic.h>
#include <ViennaRNA/fold_compound.h>
#include <ViennaRNA/params/basic.h>
#include <ViennaRNA/constraints/hard.h>
#include <ViennaRNA/constraints/soft.h>
#include <ViennaRNA/params/salt.h>


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
 *  @file     ViennaRNA/partfunc/internal.h
 *  @ingroup  eval, eval_loops, eval_loops_int
 *  @brief    Energy evaluation of internal loops for MFE and partition function calculations
 */

/**
 *  @addtogroup   eval_loops_int
 *  @{
 */



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



#endif
