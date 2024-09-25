#ifndef VIENNA_RNA_PACKAGE_EVAL_EXTERNAL_H
#define VIENNA_RNA_PACKAGE_EVAL_EXTERNAL_H

#include <ViennaRNA/datastructures/basic.h>
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

/**
 *  @file     ViennaRNA/eval/exterior.h
 *  @ingroup  eval, eval_loops, eval_loops_ext
 *  @brief    Energy evaluation of exterior loops
 */

/**
 *  @addtogroup   eval_loops_ext
 *  @{
 */


/**
 *  @name Basic free energy interface
 *  @{
 */


/**
 *  @brief  Evaluate a stem branching off the exterior loop
 *
 *  Given a base pair @f$(i,j)@f$ encoded by @em type, compute the energy contribution
 *  including dangling-end/terminal-mismatch contributions. Instead of returning the
 *  energy contribution per-se, this function returns the corresponding Boltzmann factor.
 *  If either of the adjacent nucleotides @f$(i - 1)@f$ and @f$(j+1)@f$ must not
 *  contribute stacking energy, the corresponding encoding must be @f$-1@f$.
 *
 *  @note By default, terminal mismatch energies are applied that correspond to the
 *        neighboring nucleotides provided by their encodings @p n5d and @p n3d. Whenever
 *        the encodings are negative, the implementation switches to usage of dangling
 *        end energies (for the non-negative base). If both encodings are negative, no
 *        terminal mismatch contributions are added.
 *
 *  @see vrna_exp_E_exterior_stem()
 *
 *  @param  type  The base pair encoding
 *  @param  n5d   The encoded nucleotide directly adjacent at the 5' side of the base pair (may be -1)
 *  @param  n3d   The encoded nucleotide directly adjacent at the 3' side of the base pair (may be -1)
 *  @param  p     The pre-computed energy parameters
 *  @return       The energy contribution of the introduced exterior-loop stem
 */
int
vrna_E_exterior_stem(unsigned int type,
                     int          n5d,
                     int          n3d,
                     vrna_param_t *p);


int
vrna_E_exterior_loop(unsigned int n,
                     vrna_md_t    *md);


/**
 *  @brief Evaluate the free energy of a base pair in the exterior loop
 *
 *  Evalue the free energy of a base pair connecting two nucleotides in the exterior
 *  loop and take hard constraints into account.
 *
 *  Typically, this is simply dangling end contributions of the adjacent
 *  nucleotides, potentially a terminal A-U mismatch penalty, and maybe some
 *  generic soft constraint contribution for that decomposition.
 *
 *  @note   For dangles == 1 || 3 this function also evaluates the three additional
 *          pairs (i + 1, j), (i, j - 1), and (i + 1, j - 1) and returns the minimum
 *          for all four possibilities in total.
 *
 *  @note   By default, all user-supplied hard- and soft constraints will be taken
 *          into account! Use the #VRNA_EVAL_LOOP_NO_HC and #VRNA_EVAL_LOOP_NO_SC
 *          bit flags as input for @p options to change the default behavior if necessary.
 *
 *  @see    vrna_E_exterior_stem(),
 *          #VRNA_EVAL_LOOP_NO_HC, #VRNA_EVAL_LOOP_NO_SC, #VRNA_EVAL_LOOP_NO_CONSTRAINTS
 *
 *  @param  fc    Fold compound to work on (defines the model and parameters)
 *  @param  i     5' position of the base pair
 *  @param  j     3' position of the base pair
 *  @param  options   A bit-field that specifies which aspects (not) to consider during evaluation
 *  @returns          Free energy for the terminal base pair of a stem branching off the exterior loop in deka-kal/mol or #INF if the pair is forbidden
 */
int
vrna_eval_exterior_stem(vrna_fold_compound_t  *fc,
                        unsigned int          i,
                        unsigned int          j,
                        unsigned int          options);


/* End basic interface */
/**@}*/


/**
 *  @name Boltzmann weight (partition function) interface
 *  @{
 */


/**
 *  @brief  Evaluate a stem branching off the exterior loop (Boltzmann factor version)
 *
 *  Given a base pair @f$(i,j)@f$ encoded by @em type, compute the energy contribution
 *  including dangling-end/terminal-mismatch contributions. Instead of returning the
 *  energy contribution per-se, this function returns the corresponding Boltzmann factor.
 *
 *  @note By default, terminal mismatch energies are applied that correspond to the
 *        neighboring nucleotides provided by their encodings @p n5d and @p n3d. Whenever
 *        the encodings are negative, the implementation switches to usage of dangling
 *        end energies (for the non-negative base). If both encodings are negative, no
 *        terminal mismatch contributions are added.
 *
 *  @see vrna_E_exterior()
 *
 *  @param  type  The base pair encoding
 *  @param  n5d   The encoded nucleotide directly adjacent at the 5' side of the base pair (may be -1)
 *  @param  n3d   The encoded nucleotide directly adjacent at the 3' side of the base pair (may be -1)
 *  @param  p     The pre-computed energy parameters (Boltzmann factor version)
 *  @return The Boltzmann weighted energy contribution of the introduced exterior-loop stem
 */
FLT_OR_DBL
vrna_exp_E_exterior_stem(unsigned int     type,
                         int              n5d,
                         int              n3d,
                         vrna_exp_param_t *p);


FLT_OR_DBL
vrna_exp_E_exterior_loop(unsigned int n,
                         vrna_md_t    *md);


/* End partition function interface */
/**@}*/


/* End eval_loops_ext */
/**
 * @}
 */


#ifndef VRNA_DISABLE_BACKWARD_COMPATIBILITY

/**
 *  @addtogroup   eval_deprecated
 *  @{
 */
DEPRECATED(int
           vrna_eval_ext_stem(vrna_fold_compound_t  *fc,
                              int                   i,
                              int                   j),
           "Use vrna_eval_exterior_stem() instead!");


/**
 *  @brief Compute the energy contribution of a stem branching off a loop-region
 *
 *  This function computes the energy contribution of a stem that branches off
 *  a loop region. This can be the case in multiloops, when a stem branching off
 *  increases the degree of the loop but also <I>immediately interior base pairs</I>
 *  of an exterior loop contribute free energy.
 *  To switch the behavior of the function according to the evaluation of a multiloop-
 *  or exterior-loop-stem, you pass the flag 'extLoop'.
 *  The returned energy contribution consists of a TerminalAU penalty if the pair type
 *  is greater than 2, dangling end contributions of mismatching nucleotides adjacent to
 *  the stem if only one of the si1, sj1 parameters is greater than 0 and mismatch energies
 *  if both mismatching nucleotides are positive values.
 *  Thus, to avoid incorporating dangling end or mismatch energies just pass a negative number,
 *  e.g. -1 to the mismatch argument.
 *
 *  This is an illustration of how the energy contribution is assembled:
 *  <PRE>
 *        3'  5'
 *        |   |
 *        X - Y
 *  5'-si1     sj1-3'
 *  </PRE>
 *
 *  Here, (X,Y) is the base pair that closes the stem that branches off a loop region.
 *  The nucleotides si1 and sj1 are the 5'- and 3'- mismatches, respectively. If the base pair
 *  type of (X,Y) is greater than 2 (i.e. an A-U or G-U pair, the TerminalAU penalty will be
 *  included in the energy contribution returned. If si1 and sj1 are both nonnegative numbers,
 *  mismatch energies will also be included. If one of si1 or sj1 is a negative value, only
 *  5' or 3' dangling end contributions are taken into account. To prohibit any of these mismatch
 *  contributions to be incorporated, just pass a negative number to both, si1 and sj1.
 *  In case the argument extLoop is 0, the returned energy contribution also includes
 *  the <I>internal-loop-penalty</I> of a multiloop stem with closing pair type.
 *
 *  @note   This function is threadsafe
 *
 *  @see    vrna_E_multibranch_stem(), _ExtLoop()
 *
 *  @deprecated     Please use one of the functions vrna_E_exterior_stem() and
 *                  vrna_E_multibranch_stem() instead! Use the former for cases where @p extLoop != 0
 *                  and the latter otherwise.
 *
 *  @param  type    The pair type of the first base pair un the stem
 *  @param  si1     The 5'-mismatching nucleotide
 *  @param  sj1     The 3'-mismatching nucleotide
 *  @param  extLoop A flag that indicates whether the contribution reflects the one of an exterior loop or not
 *  @param  P       The data structure containing scaled energy parameters
 *  @return         The Free energy of the branch off the loop in dcal/mol
 *
 */
DEPRECATED(int
           E_Stem(int           type,
                  int           si1,
                  int           sj1,
                  int           extLoop,
                  vrna_param_t  *P),
           "This function is obsolete. Use vrna_E_exterior_stem() or vrna_E_multibranch_stem() instead");


DEPRECATED(int
           E_ExtLoop(int          type,
                     int          si1,
                     int          sj1,
                     vrna_param_t *P),
           "Use vrna_E_exterior_stem() instead");


DEPRECATED(int
           vrna_E_ext_stem(unsigned int type,
                           int          n5d,
                           int          n3d,
                           vrna_param_t *p),
           "Use vrna_E_exterior_stem() instead!");

DEPRECATED(FLT_OR_DBL
           vrna_exp_E_ext_stem(unsigned int     type,
                               int              n5d,
                               int              n3d,
                               vrna_exp_param_t *p),
           "Use vrna_exp_E_exterior_stem() instead!");


/**
 *  This is the partition function variant of @ref E_ExtLoop()
 *  @deprecated Use vrna_exp_E_ext_stem() instead!
 *
 *  @see E_ExtLoop()
 *
 *  @return The Boltzmann weighted energy contribution of the introduced exterior-loop stem
 */
DEPRECATED(FLT_OR_DBL
           exp_E_ExtLoop(int              type,
                         int              si1,
                         int              sj1,
                         vrna_exp_param_t *P),
           "Use vrna_exp_E_exterior_stem() instead");


/**
 *  <H2>Compute the Boltzmann weighted energy contribution of a stem branching off a loop-region</H2>
 *  This is the partition function variant of @ref E_Stem()
 *
 *  @note This function is threadsafe
 *
 *  @see E_Stem()
 *
 *  @return The Boltzmann weighted energy contribution of the branch off the loop
 */
DEPRECATED(FLT_OR_DBL
           exp_E_Stem(int               type,
                      int               si1,
                      int               sj1,
                      int               extLoop,
                      vrna_exp_param_t  *P),
           "This function is obsolete. Use vrna_exp_E_exterior_stem() or vrna_exp_E_multibranch_stem() instead!");


/**
 * @}
 */


#endif

#endif
