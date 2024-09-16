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

#ifdef __GNUC__
# define INLINE inline
#else
# define INLINE
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


/**
 *  @brief Evaluate energy of multi branch loop helices stacking onto closing pair (i,j)
 *
 *  Computes total free energy for coaxial stacking of (i.j) with (i+1.k) or (k+1.j-1)
 */
int
vrna_E_mb_loop_stack(vrna_fold_compound_t *fc,
                     int                  i,
                     int                  j);


int
vrna_E_mb_loop_fast(vrna_fold_compound_t  *fc,
                    int                   i,
                    int                   j,
                    int                   *dmli1,
                    int                   *dmli2);


int
E_ml_rightmost_stem(int                   i,
                    int                   j,
                    vrna_fold_compound_t  *fc);


int
vrna_E_ml_stems_fast(vrna_fold_compound_t *fc,
                     int                  i,
                     int                  j,
                     int                  *fmi,
                     int                  *dmli);


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
vrna_bt_m(vrna_fold_compound_t    *fc,
          unsigned int            i,
          unsigned int            j,
          vrna_bps_t              bp_stack,
          vrna_bts_t              bt_stack);

/**
 *  @brief  Backtrack the decomposition of a multi branch loop closed by @f$ (i,j) @f$
 *
 *  @param    fc          The #vrna_fold_compound_t filled with all relevant data for backtracking
 *  @param    i           5' position of base pair closing the loop (will be set to 5' position
 *                        of leftmost decomposed block upon successful backtracking)
 *  @param    j           3' position of base pair closing the loop (will be set to 3' position
 *                        of rightmost decomposed block upon successful backtracking)
 *  @param    k           Split position that delimits leftmost from rightmost block, [i,k] and
 *                        [k+1, j], respectively. (Will be set upon successful backtracking)
 *  @param    en          The energy contribution of the substructure enclosed by @f$ (i,j) @f$
 *  @param    component1  Type of leftmost block (1 = ML, 2 = C)
 *  @param    component2  Type of rightmost block (1 = ML, 2 = C)
 *  @returns              1, if backtracking succeeded, 0 otherwise.
 */
unsigned int
vrna_bt_mb_loop(vrna_fold_compound_t  *fc,
                unsigned int          i,
                unsigned int          j,
                int                   en,
                vrna_bps_t            bp_stack,
                vrna_bts_t            bt_stack);


unsigned int
vrna_bt_mb_loop_split(vrna_fold_compound_t  *fc,
                      unsigned int          i,
                      unsigned int          j,
                      vrna_bps_t            bp_stack,
                      vrna_bts_t            bt_stack);


/**
 * @}
 */


#ifndef VRNA_DISABLE_BACKWARD_COMPATIBILITY

DEPRECATED(int
vrna_BT_mb_loop(vrna_fold_compound_t  *fc,
                int                   *i,
                int                   *j,
                int                   *k,
                int                   en,
                unsigned int          *component1,
                unsigned int          *component2),
          "Use vrna_bt_mb_loop() instead!");


DEPRECATED(int
vrna_BT_mb_loop_split(vrna_fold_compound_t  *fc,
                      int                   *i,
                      int                   *j,
                      int                   *k,
                      int                   *l,
                      unsigned int          *component1,
                      unsigned int          *component2,
                      vrna_bp_stack_t       *bp_stack,
                      unsigned int          *stack_count),
          "Use vrna_bt_mb_loop_split() instead!");



#endif

#endif
