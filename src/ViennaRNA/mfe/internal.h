#ifndef VIENNA_RNA_PACKAGE_MFE_INTERNAL_H
#define VIENNA_RNA_PACKAGE_MFE_INTERNAL_H

#include <math.h>

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
 *  @file     ViennaRNA/mfe/internal.h
 *  @ingroup  mfe
 *  @brief    Energy evaluation of internal loops for MFE computations
 */

/**
 *  @addtogroup   eval_loops_int
 *  @{
 */


/**
 *  @name Minimum Free Energy API
 *  @{
 */

int
vrna_mfe_internal(vrna_fold_compound_t  *fc,
                  unsigned int          i,
                  unsigned int          j);


int
vrna_mfe_internal_ext(vrna_fold_compound_t  *fc,
                      unsigned int          i,
                      unsigned int          j,
                      unsigned int          *ip,
                      unsigned int          *iq);

/* End basic interface */
/**@}*/

#ifndef VRNA_DISABLE_BACKWARD_COMPATIBILITY

DEPRECATED(int
           vrna_E_int_loop(vrna_fold_compound_t *fc,
                           int                  i,
                           int                  j),
           "Use vrna_mfe_internal() instead!");


DEPRECATED(int
           vrna_E_ext_int_loop(vrna_fold_compound_t *fc,
                               int                  i,
                               int                  j,
                               int                  *ip,
                               int                  *iq),
           "Use vrna_mfe_internal_ext() instead!");


#endif

/**
 * @}
 */


#endif
