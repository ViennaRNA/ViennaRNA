#ifndef VIENNA_RNA_PACKAGE_PARTFUNC_WINDOW_H
#define VIENNA_RNA_PACKAGE_PARTFUNC_WINDOW_H

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
 *  @file     ViennaRNA/partfunc/local.h
 *  @ingroup  pf_fold, part_func_window
 *  @brief    Partition function and equilibrium probability implementation for the sliding window algorithm
 *
 *  This file contains the implementation for sliding window partition function and equilibrium
 *  probabilities. It also provides the unpaired probability implementation from
 *  @rstinline :cite:t:`bernhart:2011` @endrst
 */


#include <ViennaRNA/datastructures/basic.h>

/**
 *  @addtogroup part_func_window
 *  @{
 */

/**
 * @brief Sliding window probability computation callback
 *
 * @callback
 * @parblock
 * This function will be called for each probability data set in the sliding
 * window probability computation implementation of vrna_probs_window().
 * The argument @a type specifies the type of probability that is passed to
 * this function.
 * @endparblock
 *
 * #### Types: ####
 *  * #VRNA_PROBS_WINDOW_BPP  - @copybrief #VRNA_PROBS_WINDOW_BPP
 *  * #VRNA_PROBS_WINDOW_UP   - @copybrief #VRNA_PROBS_WINDOW_UP
 *  * #VRNA_PROBS_WINDOW_PF   - @copybrief #VRNA_PROBS_WINDOW_PF
 *
 * The above types usually come exclusively. However, for unpaired
 * probabilities, the #VRNA_PROBS_WINDOW_UP flag is OR-ed together
 * with one of the loop type contexts
 *
 *  * #VRNA_EXT_LOOP  - @copybrief #VRNA_EXT_LOOP
 *  * #VRNA_HP_LOOP   - @copybrief #VRNA_HP_LOOP
 *  * #VRNA_INT_LOOP  - @copybrief #VRNA_INT_LOOP
 *  * #VRNA_MB_LOOP   - @copybrief #VRNA_MB_LOOP
 *  * #VRNA_ANY_LOOP  - @copybrief #VRNA_ANY_LOOP
 *
 * to indicate the particular type of data available through the @p pr
 * pointer.
 *
 * @see vrna_probs_window(), vrna_pfl_fold_up_cb()
 *
 * @param pr      An array of probabilities
 * @param pr_size The length of the probability array
 * @param i       The i-position (5') of the probabilities
 * @param max     The (theoretical) maximum length of the probability array
 * @param type    The type of data that is provided
 * @param data    Auxiliary data
 */
typedef void (*vrna_probs_window_f)(FLT_OR_DBL    *pr,
                                          int           pr_size,
                                          int           i,
                                          int           max,
                                          unsigned int  type,
                                          void          *data);

DEPRECATED(typedef void (vrna_probs_window_callback)(FLT_OR_DBL    *pr,
                                          int           pr_size,
                                          int           i,
                                          int           max,
                                          unsigned int  type,
                                          void          *data),
           "Use vrna_probs_window_f instead!");


#include <ViennaRNA/fold_compound.h>
#include <ViennaRNA/structures/problist.h>

/**
 *  @brief  Exterior loop
 */
#define VRNA_EXT_LOOP   1U

/**
 *  @brief  Hairpin loop
 */
#define VRNA_HP_LOOP    2U

/**
 *  @brief  Internal loop
 */
#define VRNA_INT_LOOP   4U

/**
 *  @brief  Multibranch loop
 */
#define VRNA_MB_LOOP    8U

/**
 *  @brief  Any loop
 */
#define VRNA_ANY_LOOP   (VRNA_EXT_LOOP | VRNA_HP_LOOP | VRNA_INT_LOOP | VRNA_MB_LOOP)


/**
 *  @brief  Trigger base pairing probabilities
 *
 *  Passing this flag to vrna_probs_window() activates callback execution for base
 *  pairing probabilities. In turn, the corresponding callback receives this flag
 *  through the @p type argument whenever base pairing probabilities are provided.
 *
 *  Detailed information for the algorithm to compute unpaired probabilities
 *  can be taken from @rstinline :cite:t:`bernhart:2005` @endrst.
 *
 *  @see  vrna_probs_window()
 */
#define VRNA_PROBS_WINDOW_BPP  4096U

/**
 *  @brief  Trigger unpaired probabilities
 *
 *  Passing this flag to vrna_probs_window() activates callback execution for
 *  unpaired probabilities. In turn, the corresponding callback receives this flag
 *  through the @p type argument whenever unpaired probabilities are provided.
 *
 *  Detailed information for the algorithm to compute unpaired probabilities
 *  can be taken from @rstinline :cite:t:`bernhart:2011` @endrst.
 *
 *  @see  vrna_probs_window()
 */
#define VRNA_PROBS_WINDOW_UP   8192U

/**
 *  @brief Trigger base pair stack probabilities
 *
 *  Passing this flag to vrna_probs_window() activates callback execution for
 *  stacking probabilities. In turn, the corresponding callback receives this flag
 *  through the @p type argument whenever stack probabilities are provided.
 *
 *  @bug  Currently, this flag is a placeholder doing nothing as the corresponding
 *        implementation for stack probability computation is missing.
 *
 *  @see  vrna_probs_window()
 */
#define VRNA_PROBS_WINDOW_STACKP   16384U

/**
 *  @brief  Trigger detailed unpaired probabilities split up into different loop type contexts
 *
 *  Passing this flag to vrna_probs_window() activates callback execution for
 *  unpaired probabilities. In contrast to #VRNA_PROBS_WINDOW_UP this flag requests
 *  unpaired probabilities to be split up into different loop type contexts. In turn,
 *  the corresponding callback receives the #VRNA_PROBS_WINDOW_UP flag OR-ed together
 *  with the corresponding loop type, i.e.:
 *
 *  * #VRNA_EXT_LOOP  - @copybrief #VRNA_EXT_LOOP
 *  * #VRNA_HP_LOOP   - @copybrief #VRNA_HP_LOOP
 *  * #VRNA_INT_LOOP  - @copybrief #VRNA_INT_LOOP
 *  * #VRNA_MB_LOOP   - @copybrief #VRNA_MB_LOOP
 *  * #VRNA_ANY_LOOP  - @copybrief #VRNA_ANY_LOOP
 *
 *  @see  vrna_probs_window(), #VRNA_PROBS_WINDOW_UP
 */
#define VRNA_PROBS_WINDOW_UP_SPLIT   32768U


/**
 *  @brief  Trigger partition function
 *
 *  Passing this flag to vrna_probs_window() activates callback execution for
 *  partition function. In turn, the corresponding callback receives this flag
 *  through it's @p type argument whenever partition function data is provided.
 *
 *  @note Instead of actually providing the partition function @f$Z@f$, the
 *        callback is always provided with the corresponding enemble free energy
 *        @f$\Delta G = - RT \ln Z@f$.
 *
 *  @see  vrna_probs_window()
 */
#define VRNA_PROBS_WINDOW_PF        65536U

/**
 *  @name Basic local partition function interface
 *  @{
 */

/**
 *  @brief  Compute various equilibrium probabilities under a sliding window approach
 *
 *  This function applies a sliding window scan for the sequence provided with the
 *  argument @p fc and reports back equilibrium probabilities through the callback
 *  function @p cb. The data reported to the callback depends on the @p options flag.
 *
 *  @note   The parameter @p ulength only affects computation and resulting data if unpaired
 *          probability computations are requested through the @p options flag.
 *
 *  #### Options: ####
 *  * #VRNA_PROBS_WINDOW_BPP      - @copybrief #VRNA_PROBS_WINDOW_BPP
 *  * #VRNA_PROBS_WINDOW_UP       - @copybrief #VRNA_PROBS_WINDOW_UP
 *  * #VRNA_PROBS_WINDOW_UP_SPLIT - @copybrief #VRNA_PROBS_WINDOW_UP_SPLIT
 *
 *  Options may be OR-ed together
 *
 *  @see  vrna_pfl_fold_cb(), vrna_pfl_fold_up_cb()
 *
 *  @param  fc            The fold compound with sequence data, model settings and precomputed energy parameters
 *  @param  ulength       The maximal length of an unpaired segment (only for unpaired probability computations)
 *  @param  cb            The callback function which collects the pair probability data for further processing
 *  @param  data          Some arbitrary data structure that is passed to the callback @p cb
 *  @param  options       Option flags to control the behavior of this function
 *  @return               0 on failure, non-zero on success
 */
int
vrna_probs_window(vrna_fold_compound_t        *fc,
                  int                         ulength,
                  unsigned int                options,
                  vrna_probs_window_f  cb,
                  void                        *data);

/* End basic interface */
/**@}*/

/**
 *  @name Simplified global partition function computation using sequence(s) or multiple sequence alignment(s)
 *  @{
 */

/**
 *  @brief  Compute base pair probabilities using a sliding-window approach
 *
 *  This is a simplified wrapper to vrna_probs_window() that given a nucleid acid sequence,
 *  a window size, a maximum base pair span, and a cutoff value computes the pair probabilities
 *  for any base pair in any window. The pair probabilities are returned as a list and the user
 *  has to take care to free() the memory occupied by the list.
 *
 *  @note This function uses default model settings! For custom model settings, we refer to
 *        the function vrna_probs_window().<br>
 *        In case of any computation errors, this function returns @p NULL
 *
 *  @see    vrna_probs_window(), vrna_pfl_fold_cb(), vrna_pfl_fold_up()
 *
 *  @param  sequence      The nucleic acid input sequence
 *  @param  window_size   The size of the sliding window
 *  @param  max_bp_span   The maximum distance along the backbone between two nucleotides that form a base pairs
 *  @param  cutoff        A cutoff value that omits all pairs with lower probability
 *  @return               A list of base pair probabilities, terminated by an entry with #vrna_ep_t.i and #vrna_ep_t.j set to 0
 */
vrna_ep_t *
vrna_pfl_fold(const char  *sequence,
              int         window_size,
              int         max_bp_span,
              float       cutoff);


/**
 *  @brief  Compute base pair probabilities using a sliding-window approach (callback version)
 *
 *  This is a simplified wrapper to vrna_probs_window() that given a nucleid acid sequence,
 *  a window size, a maximum base pair span, and a cutoff value computes the pair probabilities
 *  for any base pair in any window. It is similar to vrna_pfl_fold() but uses a callback mechanism
 *  to return the pair probabilities.
 *
 *  Read the details for vrna_probs_window() for details on the callback implementation!
 *
 *  @note This function uses default model settings! For custom model settings, we refer to
 *        the function vrna_probs_window().
 *
 *  @see    vrna_probs_window(), vrna_pfl_fold(), vrna_pfl_fold_up_cb()
 *
 *  @param  sequence      The nucleic acid input sequence
 *  @param  window_size   The size of the sliding window
 *  @param  max_bp_span   The maximum distance along the backbone between two nucleotides that form a base pairs
 *  @param  cb            The callback function which collects the pair probability data for further processing
 *  @param  data          Some arbitrary data structure that is passed to the callback @p cb
 *  @return               0 on failure, non-zero on success
 */
int
vrna_pfl_fold_cb(const char                 *sequence,
                 int                        window_size,
                 int                        max_bp_span,
                 vrna_probs_window_f cb,
                 void                       *data);


/**
 *  @brief  Compute probability of contiguous unpaired segments
 *
 *  This is a simplified wrapper to vrna_probs_window() that given a nucleic acid sequence,
 *  a maximum length of unpaired segments (@p ulength), a window size, and a maximum base
 *  pair span computes the equilibrium probability of any segment not exceeding @p ulength.
 *  The probabilities to be unpaired are returned as a 1-based, 2-dimensional matrix with
 *  dimensions @f$ N \times M @f$, where @f$N@f$ is the length of the sequence and @f$M@f$
 *  is the maximum segment length. As an example, the probability of a segment of size 5
 *  starting at position 100 is stored in the matrix entry @f$X[100][5]@f$.
 *
 *  It is the users responsibility to free the memory occupied by this matrix.
 *
 *  @note This function uses default model settings! For custom model settings, we refer to
 *        the function vrna_probs_window().
 *
 *  @param  sequence      The nucleic acid input sequence
 *  @param  ulength       The maximal length of an unpaired segment
 *  @param  window_size   The size of the sliding window
 *  @param  max_bp_span   The maximum distance along the backbone between two nucleotides that form a base pairs
 *  @return               The probabilities to be unpaired for any segment not exceeding @p ulength
 */
double **
vrna_pfl_fold_up(const char *sequence,
                 int        ulength,
                 int        window_size,
                 int        max_bp_span);


/**
 *  @brief  Compute probability of contiguous unpaired segments
 *
 *  This is a simplified wrapper to vrna_probs_window() that given a nucleic acid sequence,
 *  a maximum length of unpaired segments (@p ulength), a window size, and a maximum base
 *  pair span computes the equilibrium probability of any segment not exceeding @p ulength.
 *  It is similar to vrna_pfl_fold_up() but uses a callback mechanism to return the unpaired
 *  probabilities.
 *
 *  Read the details for vrna_probs_window() for details on the callback implementation!
 *
 *
 *  @note This function uses default model settings! For custom model settings, we refer to
 *        the function vrna_probs_window().
 *
 *  @param  sequence      The nucleic acid input sequence
 *  @param  ulength       The maximal length of an unpaired segment
 *  @param  window_size   The size of the sliding window
 *  @param  max_bp_span   The maximum distance along the backbone between two nucleotides that form a base pairs
 *  @param  cb            The callback function which collects the pair probability data for further processing
 *  @param  data          Some arbitrary data structure that is passed to the callback @p cb
 *  @return               0 on failure, non-zero on success
 */
int
vrna_pfl_fold_up_cb(const char                  *sequence,
                    int                         ulength,
                    int                         window_size,
                    int                         max_bp_span,
                    vrna_probs_window_f  cb,
                    void                        *data);


/* End simplified interface */
/**@}*/


/**@}*/

#endif
