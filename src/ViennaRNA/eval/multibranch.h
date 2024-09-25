#ifndef VIENNA_RNA_PACKAGE_EVAL_MULTIBRANCH_H
#define VIENNA_RNA_PACKAGE_EVAL_MULTIBRANCH_H

#include <ViennaRNA/utils/basic.h>
#include <ViennaRNA/datastructures/basic.h>
#include <ViennaRNA/datastructures/array.h>
#include <ViennaRNA/fold_compound.h>
#include <ViennaRNA/params/basic.h>
#include <ViennaRNA/eval/basic.h>

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
 *  @file     ViennaRNA/eval/multibranch.h
 *  @ingroup  eval, eval_loops, eval_loops_mb
 *  @brief    Energy evaluation of multibranch loops for MFE and partition function calculations
 */

/**
 *  @addtogroup  eval_loops_mb
 *  @{
 */


/**
 *  @name Basic free energy interface
 *  @{
 */


/**
 *  @brief  Evaluate the free energy contribution of a stem branching off a multibranch loop
 *
 *  This function yields the free energy contribution for the terminal base
 *  pairs of a stem branching off a multibranch loop. In essence, this consists of
 *  (i) a terminal mismatch or dangling end contribution, (ii) the score for a stem
 *  according to the affine multibranch loop model, and (iii) a terminal AU/GU penalty,
 *  if applicable.
 *
 *  @note By default, terminal mismatch energies are applied that correspond to the
 *        neighboring nucleotides provided by their encodings @p n5d and @p n3d. Whenever
 *        the encodings are negative, the implementation switches to usage of dangling
 *        end energies (for the non-negative base). If both encodings are negative, no
 *        terminal mismatch contributions are added.
 *
 *  @see vrna_exp_E_multibranch_stem(), vrna_E_exterior_stem()
 *
 *  @param  type  The base pair encoding
 *  @param  si1   The encoded nucleotide directly adjacent at the 5' side of the base pair (may be -1)
 *  @param  sj1   The encoded nucleotide directly adjacent at the 3' side of the base pair (may be -1)
 *  @param  p     The pre-computed energy parameters
 *  @return       The energy contribution of the introduced mutlibranch loop stem
 */
int
vrna_E_multibranch_stem(unsigned int  type,
                        int           si1,
                        int           sj1,
                        vrna_param_t  *P);


/* End basic interface */
/**@}*/


/**
 *  @name Boltzmann weight (partition function) interface
 *  @{
 */


/**
 *  @brief  Evaluate the free energy contribution of a stem branching off a multibranch loop  (Boltzmann factor version)
 *
 *  This function yields the free energy contribution as Boltzmann factor @f$ exp(-E/kT) @f$
 *  for the terminal base pairs of a stem branching off a multibranch loop. In essence, this
 *  consists of (i) a terminal mismatch or dangling end contribution, (ii) the score for a
 *  stem according to the affine multibranch loop model, and (iii) a terminal AU/GU penalty,
 *  if applicable.
 *
 *  @note By default, terminal mismatch energies are applied that correspond to the
 *        neighboring nucleotides provided by their encodings @p n5d and @p n3d. Whenever
 *        the encodings are negative, the implementation switches to usage of dangling
 *        end energies (for the non-negative base). If both encodings are negative, no
 *        terminal mismatch contributions are added.
 *
 *  @see vrna_E_multibranch_stem(), vrna_exp_E_exterior_stem()
 *
 *  @param  type  The base pair encoding
 *  @param  si1   The encoded nucleotide directly adjacent at the 5' side of the base pair (may be -1)
 *  @param  sj1   The encoded nucleotide directly adjacent at the 3' side of the base pair (may be -1)
 *  @param  p     The pre-computed energy parameters (Boltzmann factor version)
 *  @return       The Boltzmann factor of the energy contribution for the introduced mutlibranch loop stem
 */
FLT_OR_DBL
vrna_exp_E_multibranch_stem(unsigned int      type,
                            int               si1,
                            int               sj1,
                            vrna_exp_param_t  *P);


/* End partition function interface */
/**@}*/

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
 *
 *  @param  A The pair type of the stem-closing pair
 *  @param  B The 5'-mismatching nucleotide
 *  @param  C The 3'-mismatching nucleotide
 *  @param  D The datastructure containing scaled energy parameters
 *  @return   The energy contribution of the introduced multiloop stem
 */
DEPRECATED(PRIVATE INLINE int
           E_MLstem(int           type,
                    int           si1,
                    int           sj1,
                    vrna_param_t  *P),
           "Use vrna_E_multibranch_stem() instead!");


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
 *
 *  @see E_MLstem()
 *
 *  @return The Boltzmann weighted energy contribution of the introduced multiloop stem
 */
DEPRECATED(PRIVATE INLINE FLT_OR_DBL
           exp_E_MLstem(int               type,
                        int               si1,
                        int               sj1,
                        vrna_exp_param_t  *P),
           "Use vrna_exp_E_multibranch_stem() instead!");


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
