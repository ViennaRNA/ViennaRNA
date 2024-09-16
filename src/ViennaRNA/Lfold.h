#ifndef VIENNA_RNA_PACKAGE_LFOLD_H
#define VIENNA_RNA_PACKAGE_LFOLD_H

#ifndef VRNA_DISABLE_BACKWARD_COMPATIBILITY

/**
 *  @file Lfold.h
 *  @ingroup mfe_window_deprecated
 *  @brief    Functions for locally optimal MFE structure prediction
 */

#ifdef VRNA_WARN_DEPRECATED
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

#include <ViennaRNA/mfe/local.h>

/**
 *  @brief The local analog to fold().
 *
 *  Computes the minimum free energy structure including only base pairs
 *  with a span smaller than 'maxdist'
 *
 *  @ingroup mfe_window_deprecated
 *
 *  @deprecated Use vrna_mfe_window() instead!
 */
DEPRECATED(float Lfold(const char *string,
                       const char *structure,
                       int        maxdist),
           "Use vrna_Lfold() or vrna_Lfold_cb() instead");

#ifdef VRNA_WITH_SVM
/**
 *  @brief
 *
 *  @ingroup mfe_window_deprecated
 *
 *  @deprecated Use vrna_mfe_window_zscore() instead!
 */
DEPRECATED(float Lfoldz(const char  *string,
                        const char  *structure,
                        int         maxdist,
                        int         zsc,
                        double      min_z),
           "Use vrna_Lfoldz() or vrna_Lfoldz_cb() instead");
#endif

/**
 *  @brief
 *
 *  @ingroup mfe_window_deprecated
 */
DEPRECATED(float aliLfold(const char  **AS,
                          const char  *structure,
                          int         maxdist),
           "Use vrna_aliLfold() or vrna_aliLfold_cb() instead");


/**
 *  @brief
 *
 *  @ingroup mfe_window_deprecated
 *
 */
DEPRECATED(float aliLfold_cb(const char               **AS,
                             int                      maxdist,
                             vrna_mfe_window_f cb,
                             void                     *data),
           "Use vrna_aliLfold() or vrna_aliLfold_cb() instead");


#endif

#endif
