#ifndef VIENNA_RNA_PACKAGE_PK_PLEX_H
#define VIENNA_RNA_PACKAGE_PK_PLEX_H

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
 *  @file pk_plex.h
 *  @ingroup  pseudoknots
 *  @brief    Heuristics for two-step pseudoknot forming interaction predictions
 */

/**
 *  @addtogroup   pseudoknots
 *  @{
 */

/**
 *  @brief Pseudoknot loop scoring function prototype
 *
 *  This function is used to evaluate a formed pseudoknot (PK) interaction in vrna_pk_plex().
 *  It is supposed to take a PK-free secondary structure as input and coordinates
 *  of an additional interaction site. From this data, the energy (penalty) to score the
 *  PK loop is derived and returned in decakal/mol.
 *  Upon passing zero in any of the interaction site coordinates (@p start_5, @p end_5,
 *  @p start_3, @p end_3) or a @em NULL pointer in @p pt, the function must return a
 *  PK loop score. This minimum PK loop score is used in the first phase of the heuristic
 *  implemented in vrna_pk_plex() to assess whether a particular interaction is further
 *  taken into account in a later, more thorough evaluation step.
 *
 *  The simplest scoring function would simply return a constant score for any
 *  PK loop, no matter what type of loop is formed and how large the loop is. This
 *  is the default if vrna_pk_plex_opt_defaults() or vrna_pk_plex_opt() is used to
 *  generate options for vrna_pk_plex().
 *
 *  @see vrna_pk_plex_opt_fun(), vrna_pk_plex()
 *
 *  @param  pt        The secondary structure (without pseudoknot) in pair table notation
 *  @param  start_5   The start coordinate of the 5' site of the pseudoknot interaction
 *  @param  end_5     The end coordinate of the 5' site of the pseudoknot interaction
 *  @param  start_3   The start coordinate of the 3' site of the pseudoknot interaction
 *  @param  end_3     The end coordinate of the 3' site of the pseudoknot interaction
 *  @param  data      An arbitrary data structure passed from the calling function
 *  @return           The energy (penalty) of the resulting pseudoknot
 */
typedef int (*vrna_pk_plex_score_f)(const short *pt,
                                          int         start_5,
                                          int         end_5,
                                          int         start_3,
                                          int         end_3,
                                          void        *data);

DEPRECATED(typedef int (vrna_callback_pk_plex_score)(const short *pt,
                                          int         start_5,
                                          int         end_5,
                                          int         start_3,
                                          int         end_3,
                                          void        *data),
          "Use vrna_pk_plex_score_f instead!");


/**
 *  @brief  RNA PKplex options object
 *
 *  @see  vrna_pk_plex_opt_defaults(), vrna_pk_plex_opt(), vrna_pk_plex_opt_fun(),
 *        vrna_pk_plex(), #vrna_pk_plex_score_f
 */
typedef struct vrna_pk_plex_option_s *vrna_pk_plex_opt_t;

/**
 *  @brief  Convenience typedef for results of the RNA PKplex prediction
 *
 *  @see #vrna_pk_plex_results_s, vrna_pk_plex()
 */
typedef struct vrna_pk_plex_result_s vrna_pk_plex_t;

#include <ViennaRNA/fold_compound.h>

/**
 *  @brief A result of the RNA PKplex interaction prediction
 *
 *  @see #vrna_pk_plex_t
 */
struct vrna_pk_plex_result_s {
  char          *structure; /**< @brief Secondary Structure in dot-bracket notation */
  double        energy;     /**< @brief Net free energy in kcal/mol */
  double        dGpk;       /**< @brief Free energy of PK loop in kcal/mol */
  double        dGint;      /**< @brief Free energy of PK forming duplex interaction */
  double        dG1;        /**< @brief Opening energy for the 5' interaction site used in the heuristic */
  double        dG2;        /**< @brief Opening energy for the 3' interaction site used in the heuristic */
  unsigned int  start_5;    /**< @brief Start coordinate of the 5' interaction site */
  unsigned int  end_5;      /**< @brief End coordinate of the 5' interaction site */
  unsigned int  start_3;    /**< @brief Start coordinate of the 3' interaction site */
  unsigned int  end_3;      /**< @brief End coordinate of the 3' interaction site */
};

/**
 *  @brief Predict Pseudoknot interactions in terms of a two-step folding process
 *
 *  Computes simple pseudoknot interactions according to the PKplex algorithm. This
 *  simple heuristic first compiles a list of potential interaction sites that may
 *  form a pseudoknot. The resulting candidate interactions are then fixed and an
 *  PK-free MFE structure for the remainder of the sequence is computed.
 *
 *  The @p accessibility argument is a list of opening energies for potential
 *  interaction sites. It is used in the first step of the algorithm to identify
 *  potential interactions. Upon passing @em NULL, the opening energies are determined
 *  automatically based on the current model settings.
 *
 *  Depending on the @p options, the function can return the MFE (incl. PK loops)
 *  or suboptimal structures within an energy band around the MFE. The PK loop is
 *  internally scored by a scoring function that in the simplest cases assigns a
 *  constant value for each PK loop. More complicated scoring functions can be
 *  passed as well, see #vrna_pk_plex_score_f and vrna_pk_plex_opt_fun().
 *
 *  The function returns @em NULL on any error. Otherwise, a list of structures
 *  and interaction coordinates with corresponding energy contributions is returned.
 *  If no PK-interaction that satisfies the options is found, the list only
 *  consists of the PK-free MFE structure.
 *
 *  @param  fc  fold compound with the input sequence and model settings
 *  @param  accessibility   An array of opening energies for the implemented heuristic (maybe @em NULL)
 *  @param  options         An #vrna_pk_plex_opt_t options data structure that determines the algorithm parameters
 *  @return                 A list of potentially pseudoknotted structures (Last element in the list indicated by @em NULL value in #vrna_pk_plex_result_s.structure)
 */
vrna_pk_plex_t *
vrna_pk_plex(vrna_fold_compound_t *fc,
             const int            **accessibility,
             vrna_pk_plex_opt_t   options);


/**
 *  @brief  Obtain a list of opening energies suitable for PKplex computations
 *
 *  @see vrna_pk_plex()
 *
 *  @param  fc        fold compound with the input sequence and model settings for accessibility computations
 *  @param  unpaired  The maximum number of unpaired nucleotides, i.e. length of interaction
 *  @param  cutoff    A cutoff value for unpaired probabilities
 *  @return           Opening energies as required for vrna_pk_plex()
 */
int **
vrna_pk_plex_accessibility(vrna_fold_compound_t *fc,
                           unsigned int         unpaired,
                           double               cutoff);


/**
 *  @brief  Default options for PKplex algorithm
 *
 *  @see vrna_pk_plex(), vrna_pk_plex_opt(), vrna_pk_plex_opt_fun()
 *
 *  @return An options data structure suitabe for PKplex computations
 */
vrna_pk_plex_opt_t
vrna_pk_plex_opt_defaults(void);


/**
 *  @brief  Simple options for PKplex algorithm
 *
 *  @see vrna_pk_plex(), vrna_pk_plex_opt_defaults(), vrna_pk_plex_opt_fun()
 *
 *  @param  delta                   Size of energy band around MFE for suboptimal results in dekacal/mol
 *  @param  max_interaction_length  Maximum length of interaction
 *  @param  pk_penalty              Energy constant to score the PK forming loop
 *  @return                         An options data structure suitabe for PKplex computations
 */
vrna_pk_plex_opt_t
vrna_pk_plex_opt(unsigned int delta,
                 unsigned int max_interaction_length,
                 int          pk_penalty);


/**
 *  @brief  Simple options for PKplex algorithm
 *
 *  @see vrna_pk_plex(), vrna_pk_plex_opt_defaults(), vrna_pk_plex_opt(), #vrna_pk_plex_score_f
 *
 *  @param  delta                   Size of energy band around MFE for suboptimal results in dekacal/mol
 *  @param  max_interaction_length  Maximum length of interaction
 *  @param  scoring_function        Energy evaluating function to score the PK forming loop
 *  @param  scoring_data            An arbitrary data structure passed to the scoring function (maybe @em NUL)
 *  @return                         An options data structure suitabe for PKplex computations
 */
vrna_pk_plex_opt_t
vrna_pk_plex_opt_fun(unsigned int                 delta,
                     unsigned int                 max_interaction_length,
                     vrna_pk_plex_score_f  scoring_function,
                     void                         *scoring_data);


/**
 * @}
 */

#endif
