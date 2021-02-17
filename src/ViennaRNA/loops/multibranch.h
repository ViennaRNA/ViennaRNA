#ifndef VIENNA_RNA_PACKAGE_LOOPS_MULTIBRANCH_H
#define VIENNA_RNA_PACKAGE_LOOPS_MULTIBRANCH_H

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
 *  @file     ViennaRNA/loops/multibranch.h
 *  @ingroup  eval, eval_loops, eval_loops_mb
 *  @brief    Energy evaluation of multibranch loops for MFE and partition function calculations
 */

/**
 *  @addtogroup  eval_loops_mb
 *  @{
 *  @brief  Functions to evaluate the free energy contributions for mutlibranch loops
 */


/**
 *  @name Basic free energy interface
 *  @{
 */


/**
 *  @brief Evaluate energy of a multi branch helices stacking onto closing pair (i,j)
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
 *  @name Boltzmann weight (partition function) interface
 *  @{
 */


/**
 *  @brief  Auxiliary helper arrays for fast exterior loop computations
 *
 *  @see vrna_exp_E_ml_fast_init(), vrna_exp_E_ml_fast_rotate(),
 *  vrna_exp_E_ml_fast_free(), vrna_exp_E_ml_fast()
 */
typedef struct vrna_mx_pf_aux_ml_s *vrna_mx_pf_aux_ml_t;


FLT_OR_DBL
vrna_exp_E_mb_loop_fast(vrna_fold_compound_t  *fc,
                        int                   i,
                        int                   j,
                        vrna_mx_pf_aux_ml_t   aux_mx);


vrna_mx_pf_aux_ml_t
vrna_exp_E_ml_fast_init(vrna_fold_compound_t *fc);


void
vrna_exp_E_ml_fast_rotate(vrna_mx_pf_aux_ml_t aux_mx);


void
vrna_exp_E_ml_fast_free(vrna_mx_pf_aux_ml_t aux_mx);


const FLT_OR_DBL *
vrna_exp_E_ml_fast_qqm(vrna_mx_pf_aux_ml_t aux_mx);


const FLT_OR_DBL *
vrna_exp_E_ml_fast_qqm1(vrna_mx_pf_aux_ml_t aux_mx);


FLT_OR_DBL
vrna_exp_E_ml_fast(vrna_fold_compound_t *fc,
                   int                  i,
                   int                  j,
                   vrna_mx_pf_aux_ml_t  aux_mx);


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
int
vrna_BT_mb_loop(vrna_fold_compound_t  *fc,
                int                   *i,
                int                   *j,
                int                   *k,
                int                   en,
                int                   *component1,
                int                   *component2);


int
vrna_BT_mb_loop_fake(vrna_fold_compound_t *fc,
                     int                  *u,
                     int                  *i,
                     int                  *j,
                     vrna_bp_stack_t      *bp_stack,
                     int                  *stack_count);


int
vrna_BT_mb_loop_split(vrna_fold_compound_t  *fc,
                      int                   *i,
                      int                   *j,
                      int                   *k,
                      int                   *l,
                      int                   *component1,
                      int                   *component2,
                      vrna_bp_stack_t       *bp_stack,
                      int                   *stack_count);


/**
 * @}
 */


#ifndef VRNA_DISABLE_BACKWARD_COMPATIBILITY

/**
 *  @addtogroup   eval_deprecated
 *  @{
 */

/**
 *  @def E_MLstem(A,B,C,D)
 *  <H2>Compute the Energy contribution of a Multiloop stem</H2>
 *  This definition is a wrapper for the E_Stem() funtion.
 *  It is substituted by an E_Stem() funtion call with argument
 *  extLoop=0, so the energy contribution returned reflects a
 *  stem introduced in a multiloop.<BR>
 *  As for the parameters B (si1) and C (sj1) of the substituted
 *  E_Stem() function, you can inhibit to take 5'-, 3'-dangles
 *  or mismatch contributions to be taken into account by passing
 *  -1 to these parameters.
 *
 *  @see    E_Stem()
 *  @param  A The pair type of the stem-closing pair
 *  @param  B The 5'-mismatching nucleotide
 *  @param  C The 3'-mismatching nucleotide
 *  @param  D The datastructure containing scaled energy parameters
 *  @return   The energy contribution of the introduced multiloop stem
 */
PRIVATE INLINE int
E_MLstem(int          type,
         int          si1,
         int          sj1,
         vrna_param_t *P)
{
  int energy = 0;

  if (si1 >= 0 && sj1 >= 0)
    energy += P->mismatchM[type][si1][sj1];
  else if (si1 >= 0)
    energy += P->dangle5[type][si1];
  else if (sj1 >= 0)
    energy += P->dangle3[type][sj1];

  if (type > 2)
    energy += P->TerminalAU;

  energy += P->MLintern[type];

  return energy;
}


/**
 *  @def exp_E_MLstem(A,B,C,D)
 *  This is the partition function variant of @ref E_MLstem()
 *  @see E_MLstem()
 *  @return The Boltzmann weighted energy contribution of the introduced multiloop stem
 */
PRIVATE INLINE FLT_OR_DBL
exp_E_MLstem(int              type,
             int              si1,
             int              sj1,
             vrna_exp_param_t *P)
{
  double energy = 1.0;

  if (si1 >= 0 && sj1 >= 0)
    energy = P->expmismatchM[type][si1][sj1];
  else if (si1 >= 0)
    energy = P->expdangle5[type][si1];
  else if (sj1 >= 0)
    energy = P->expdangle3[type][sj1];

  if (type > 2)
    energy *= P->expTermAU;

  energy *= P->expMLintern[type];
  return (FLT_OR_DBL)energy;
}


#endif

/**
 * @}
 */

#endif
