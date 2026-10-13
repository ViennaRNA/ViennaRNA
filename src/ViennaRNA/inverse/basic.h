#ifndef VIENNA_RNA_PACKAGE_INVERSE_H
#define VIENNA_RNA_PACKAGE_INVERSE_H

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
 *  @file     ViennaRNA/inverse/basic.h
 *  @ingroup  inverse_fold
 *  @brief    Inverse folding routines
 */

/**
 *  @addtogroup inverse_fold
 *  @{
 *  
 *  @brief RNA sequence design
 */

/**
 *  \brief This global variable points to the allowed bases, initially "AUGC".
 *  It can be used to design sequences from reduced alphabets.
 */
DEPRECATED(extern char *symbolset,
           "This variable is obsolete");

/** when to stop inverse_pf_fold() */
DEPRECATED(extern  float final_cost,
           "This variable is obsolete");
/** default 0: try to minimize structure distance even if no exact solution can be found */
DEPRECATED(extern  int   give_up,
           "This variable is obsolete");
/** print out substructure on which inverse_fold() fails */
DEPRECATED(extern  int   inv_verbose,
           "This variable is obsolete");

/**
 *  \brief Find sequences with predefined structure
 * 
 *  This function searches for a sequence with minimum free energy structure
 *  provided in the parameter 'target', starting with sequence 'start'.
 *  It returns 0 if the search was successful, otherwise a structure distance
 *  in terms of the energy difference between the search result and the actual
 *  target 'target' is returned. The found sequence is returned in 'start'.
 *  If #give_up is set to 1, the function will return as soon as it is
 *  clear that the search will be unsuccessful, this speeds up the algorithm
 *  if you are only interested in exact solutions.
 * 
 *  \param  start   The start sequence
 *  \param  target  The target secondary structure in dot-bracket notation
 *  \return         The distance to the target in case a search was unsuccessful, 0 otherwise
 */
float inverse_fold( char *start,
                    const char *target);

/**
 *  \brief Find sequence that maximizes probability of a predefined structure
 * 
 *  This function searches for a sequence with maximum probability to fold into
 *  the provided structure 'target' using the partition function algorithm.
 *  It returns \f$-kT \cdot \log(p)\f$ where \f$p\f$ is the frequency of 'target' in
 *  the ensemble of possible structures. This is usually much slower than
 *  inverse_fold().
 * 
 *  \param  start   The start sequence
 *  \param  target  The target secondary structure in dot-bracket notation
 *  \return         The distance to the target in case a search was unsuccessful, 0 otherwise
 */
float inverse_pf_fold(char *start,
                      const char *target);

/**
 *  @}
 */
#endif
