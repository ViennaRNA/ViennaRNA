#ifndef VIENNA_RNA_PACKAGE_MFE_MULTIBRANCH_H
#define VIENNA_RNA_PACKAGE_MFE_MULTIBRANCH_H

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
 *  @file     ViennaRNA/mfe/multibranch.h
 *  @ingroup  mfe
 *  @brief    Energy evaluation of multibranch loops for MFE
 */

/**
 *  @addtogroup  eval_loops_mb
 *  @{
 */


/**
 *  @name Minimum Free Energy API
 *  @{
 */

typedef struct vrna_mx_mfe_aux_ml_s *vrna_mx_mfe_aux_ml_t;

vrna_mx_mfe_aux_ml_t
vrna_mfe_multibranch_fast_init(unsigned int length);


void
vrna_mfe_multibranch_fast_rotate(vrna_mx_mfe_aux_ml_t aux);


void
vrna_mfe_multibranch_fast_free(vrna_mx_mfe_aux_ml_t aux);


int
vrna_mfe_multibranch_loop_fast(vrna_fold_compound_t *fc,
                               unsigned int         i,
                               unsigned int         j,
                               vrna_mx_mfe_aux_ml_t helpers);


int
vrna_mfe_multibranch_stems_fast(vrna_fold_compound_t        *fc,
                                unsigned int                i,
                                unsigned int                j,
                                struct vrna_mx_mfe_aux_ml_s *helpers);


int
vrna_mfe_multibranch_m2_fast(vrna_fold_compound_t         *fc,
                             unsigned int                 i,
                             unsigned int                 j,
                             struct vrna_mx_mfe_aux_ml_s  *helpers);


/**
 *  @brief Evaluate energy of multi branch loop helices stacking onto closing pair (i,j)
 *
 *  Computes total free energy for coaxial stacking of (i.j) with (i+1.k) or (k+1.j-1)
 */
int
vrna_mfe_multibranch_loop_stack(vrna_fold_compound_t  *fc,
                                unsigned int          i,
                                unsigned int          j);


int
vrna_mfe_multibranch_m1(vrna_fold_compound_t  *fc,
                        unsigned int          i,
                        unsigned int          j);

/* End MFE interface */
/** @} */

/* End eval_loops_mb */
/** @} */


#ifndef VRNA_DISABLE_BACKWARD_COMPATIBILITY

DEPRECATED(int
           vrna_E_mb_loop_stack(vrna_fold_compound_t  *fc,
                                int                   i,
                                int                   j),
           "Use vrna_mfe_multibranch_loop_stack() instead!");


DEPRECATED(int
           vrna_E_mb_loop_fast(vrna_fold_compound_t *fc,
                               int                  i,
                               int                  j,
                               int                  *dmli1,
                               int                  *dmli2),
           "Use vrna_mfe_multibranch_loop_fast() instead!");


DEPRECATED(int
           E_ml_rightmost_stem(int                  i,
                               int                  j,
                               vrna_fold_compound_t *fc),
           "Use vrna_mfe_multibranch_m1() instead!");


DEPRECATED(int
           vrna_E_ml_stems_fast(vrna_fold_compound_t  *fc,
                                int                   i,
                                int                   j,
                                int                   *fmi,
                                int                   *dmli),
           "Use vrna_mfe_multibranch_stems_fast() instead!");

#endif

#endif
