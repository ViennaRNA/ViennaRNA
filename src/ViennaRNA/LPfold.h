#ifndef VIENNA_RNA_PACKAGE_LPFOLD_H
#define VIENNA_RNA_PACKAGE_LPFOLD_H

/**
 *  @file     LPfold.h
 *  @ingroup  local_fold
 *  @brief    Partition function and equilibrium probability implementation for the sliding window algorithm
 *
 *  This file contains the implementation for sliding window partition function and equilibrium
 *  probabilities. It also provides the unpaired probability implementation from Bernhart et al.
 *  2011 @cite bernhart:2011
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
 * @ingroup  local_pf_fold
 * @see vrna_probs_window(), vrna_pfl_fold_up_cb()
 *
 * @param pr      An array of probabilities
 * @param pr_size The length of the probability array
 * @param i       The i-position (5') of the probabilities
 * @param max     The (theoretical) maximum length of the probability array
 * @param type    The type of data that is provided
 * @param data    Auxiliary data
 */

typedef void (vrna_probs_window_callback)(FLT_OR_DBL    *pr,
                                          int           pr_size,
                                          int           i,
                                          int           max,
                                          unsigned int  type,
                                          void          *data);

#include <ViennaRNA/data_structures.h>
#include <ViennaRNA/params.h>

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
 *  can be taken from @cite bernhart:2005.
 *
 *  @see  vrna_probs_window()
 *  @ingroup  local_pf_fold
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
 *  can be taken from @cite bernhart:2011.
 *
 *  @see  vrna_probs_window()
 *  @ingroup  local_pf_fold
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
 *  @ingroup  local_pf_fold
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
 *  @ingroup  local_pf_fold
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
 *  @ingroup  local_pf_fold
 */
#define VRNA_PROBS_WINDOW_PF        65536U

/**
 *  @brief  Compute base pair probabilities using a sliding-window approach
 *
 *  This is a simplified wrapper to vrna_probs_window() that given a nucleid acid sequence,
 *  a window size, a maximum base pair span, and a cutoff value computes the pair probabilities
 *  for any base pair in any window. The pair probabilities are returned as a list and the user
 *  has to take care to free() the memory occupied by the list.
 *
 *  @note This function uses default model settings! For custom model settings, we refer to
 *        the function vrna_probs_window().
 *
 *  @note In case of any computation errors, this function returns @p NULL
 *
 *  @see    vrna_probs_window(), vrna_pfl_fold_cb(), vrna_pfl_fold_up()
 *
 *  @ingroup  local_pf_fold
 *  @param  sequence      The nucleic acid input sequence
 *  @param  window_size   The size of the sliding window
 *  @param  max_bp_span   The maximum distance along the backbone between two nucleotides that form a base pairs
 *  @param  cutoff        A cutoff value that omits all pairs with lower probability
 *  @return               A list of base pair probabilities, terminated by an entry with @ref vrna_ep_t.i and @ref vrna_ep_t.j set to 0
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
 *  @ingroup  local_pf_fold
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
                 vrna_probs_window_callback *cb,
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
 *  @ingroup local_pf_fold
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
 *  @ingroup local_pf_fold
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
                    vrna_probs_window_callback  *cb,
                    void                        *data);


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
 *  @ingroup local_pf_fold
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
                  vrna_probs_window_callback  *cb,
                  void                        *data);


#ifndef VRNA_DISABLE_BACKWARD_COMPATIBILITY

/**
 *  \brief
 *
 *  \ingroup local_pf_fold
 *
 *  \param  length
 */
DEPRECATED(void update_pf_paramsLP(int length),
"This function is obsolete");


/**
 *  \brief
 *
 *  \ingroup local_pf_fold
 *
 */
DEPRECATED(void update_pf_paramsLP_par(int              length,
                                       vrna_exp_param_t *parameters),
"Use the new API with vrna_folc_compound_t instead");


/**
 *  \brief Compute partition functions for locally stable secondary structures
 *
 *  pfl_fold computes partition functions for every window of size
 *  'winSize' possible in a RNA molecule, allowing only pairs with a span
 *  smaller than 'pairSize'. It returns the mean pair probabilities averaged
 *  over all windows containing the pair in 'pl'. 'winSize' should
 *  always be >= 'pairSize'. Note that in contrast to Lfold(),
 *  bases outside of the window do not influence the structure at all. Only
 *  probabilities higher than 'cutoffb' are kept.
 *
 *  If 'pU' is supplied (i.e is not the NULL pointer), pfl_fold()
 *  will also compute the mean probability that regions of length 'u' and smaller are
 *  unpaired. The parameter 'u' is supplied in 'pup[0][0]'. On return
 *  the 'pup' array will contain these probabilities, with the entry on
 *  'pup[x][y]' containing the mean probability that x and the y-1
 *  preceding bases are unpaired. The 'pU' array needs to be large
 *  enough to hold n+1 float* entries, where n is the sequence length.
 *
 *  If an array dpp2 is supplied, the probability of base pair (i,j)
 *  given that there already exists a base pair (i+1,j-1) is also
 *  computed and saved in this array. If pUfp is given (i.e. not NULL), pU
 *  is not saved but put out imediately. If spup is given (i.e. is not NULL),
 *  the pair probabilities in pl are not saved but put out imediately.
 *
 *  \ingroup local_pf_fold
 *
 *  \param  sequence  RNA sequence
 *  \param  winSize   size of the window
 *  \param  pairSize  maximum size of base pair
 *  \param  cutoffb   cutoffb for base pairs
 *  \param  pU        array holding all unpaired probabilities
 *  \param  dpp2  array of dependent pair probabilities
 *  \param  pUfp     file pointer for pU
 *  \param  spup     file pointer for pair probabilities
 *  \return list of pair probabilities
 */
DEPRECATED(vrna_ep_t *pfl_fold(char          *sequence,
                               int           winSize,
                               int           pairSize,
                               float         cutoffb,
                               double        **pU,
                               vrna_ep_t     **dpp2,
                               FILE          *pUfp,
                               FILE          *spup),
"Use vrna_pfl_fold(), vrna_pfl_fold_cb(), vrna_pfl_fold_up(), or vrna_pfl_fold_up_cb() instead");


/**
 *  \brief Compute partition functions for locally stable secondary structures
 *
 *  \ingroup local_pf_fold
 *
 */
DEPRECATED(vrna_ep_t *pfl_fold_par(char              *sequence,
                                      int               winSize,
                                      int               pairSize,
                                      float             cutoffb,
                                      double            **pU,
                                      vrna_ep_t         **dpp2,
                                      FILE              *pUfp,
                                      FILE              *spup,
                                      vrna_exp_param_t  *parameters),
"Use the new API and vrna_probs_window() instead");


DEPRECATED(void putoutpU_prob_par(double            **pU,
                                  int               length,
                                  int               ulength,
                                  FILE              *fp,
                                  int               energies,
                                  vrna_exp_param_t  *parameters),
"");


/**
 *  \brief Writes the unpaired probabilities (pU) or opening energies into a file
 *
 *  Can write either the unpaired probabilities (accessibilities) pU or
 *  the opening energies -log(pU)kT into a file
 *
 *  \ingroup local_pf_fold
 *
 *  \param  pU       pair probabilities
 *  \param  length   length of RNA sequence
 *  \param  ulength  maximum length of unpaired stretch
 *  \param  fp file pointer of destination file
 *  \param  energies  switch to put out as  opening energies
 */
DEPRECATED(void    putoutpU_prob(double **pU,
                                 int    length,
                                 int    ulength,
                                 FILE   *fp,
                                 int    energies),
"");


DEPRECATED(void putoutpU_prob_bin_par(double            **pU,
                                      int               length,
                                      int               ulength,
                                      FILE              *fp,
                                      int               energies,
                                      vrna_exp_param_t  *parameters),
"");


/**
 *  \brief Writes the unpaired probabilities (pU) or opening energies into a binary file
 *
 *  Can write either the unpaired probabilities (accessibilities) pU or
 *  the opening energies -log(pU)kT into a file
 *
 *  \ingroup local_pf_fold
 *
 *  \param  pU       pair probabilities
 *  \param  length   length of RNA sequence
 *  \param  ulength  maximum length of unpaired stretch
 *  \param  fp file pointer of destination file
 *  \param  energies  switch to put out as  opening energies
 */
DEPRECATED(void    putoutpU_prob_bin(double **pU,
                                     int    length,
                                     int    ulength,
                                     FILE   *fp,
                                     int    energies),
"");


/**
 *  Dunno if this function was ever used by external programs linking to RNAlib, but it
 *  was declared PUBLIC before.
 *  Anyway, never use this function as it will be removed soon and does nothing at all
 */
DEPRECATED(void init_pf_foldLP(int length),
"This function is obsolete");

#endif

#endif
