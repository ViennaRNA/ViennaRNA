#ifndef VIENNA_RNA_PACKAGE_STRUCTURES_MEA_H
#define VIENNA_RNA_PACKAGE_STRUCTURES_MEA_H

#include <ViennaRNA/datastructures/basic.h>
#include <ViennaRNA/params/basic.h>

/**
 *  @file MEA.h
 *  @ingroup subopt_and_representatives
 *  @brief Computes a MEA (maximum expected accuracy) structure.
 */

/**
 *  @brief Compute a MEA (maximum expected accuracy) structure
 *
 *  The algorithm maximizes the expected accuracy
 *
 *  @f[
 *    A(S) = \sum_{(i,j) \in S} 2 \gamma p_{ij} + \sum_{i \notin S} p^u_i
 *  @f]
 *
 *  Higher values of @f$\gamma@f$ result in more base pairs of lower
 *  probability and thus higher sensitivity. Low values of @f$\gamma@f$ result in structures
 *  containing only highly likely pairs (high specificity).
 *  The code of the MEA function also demonstrates the use of sparse dynamic
 *  programming scheme to reduce the time and memory complexity of folding.
 *
 *  @pre  vrna_pf() must be executed on input parameter @p fc
 *
 *  @ingroup  mea_fold
 *
 *  @param  fc    The fold compound data structure with pre-filled base pair probability matrix
 *  @param  gamma The weighting factor for base pairs vs. unpaired nucleotides
 *  @param  mea   A pointer to a variable where the MEA value will be written to
 *  @return       An MEA structure (or NULL on any error)
 */
char *
vrna_MEA(vrna_fold_compound_t *fc,
         double               gamma,
         float                *mea);


/**
 *  @brief Compute a MEA (maximum expected accuracy) structure from a list of probabilities
 *
 *  The algorithm maximizes the expected accuracy
 *
 *  @f[
 *    A(S) = \sum_{(i,j) \in S} 2 \gamma p_{ij} + \sum_{i \notin S} p^u_i
 *  @f]
 *
 *  Higher values of @f$\gamma@f$ result in more base pairs of lower
 *  probability and thus higher sensitivity. Low values of @f$\gamma@f$ result in structures
 *  containing only highly likely pairs (high specificity).
 *  The code of the MEA function also demonstrates the use of sparse dynamic
 *  programming scheme to reduce the time and memory complexity of folding.
 *
 *  @note The unpaired probabilities @f$p^u_i = 1 - \sum_{j \neq i} p_{ij}@f$ are usually
 *        computed from the supplied pairing probabilities @f$p_{ij}@f$ as stored in @p plist
 *        entries of type #VRNA_PLIST_TYPE_BASEPAIR. To overwrite individual @f$p^u_o@f$
 *        values simply add entries with type #VRNA_PLIST_TYPE_UNPAIRED<br>
 *        To include G-Quadruplex support, the corresponding field in @p md must be set.
 *
 *  @ingroup  mea_fold
 *
 *  @param  plist     A list of base pair probabilities the MEA structure is computed from
 *  @param  sequence  The RNA sequence that corresponds to the list of probability values
 *  @param  gamma     The weighting factor for base pairs vs. unpaired nucleotides
 *  @param  md        A model details data structure (maybe NULL)
 *  @param  mea       A pointer to a variable where the MEA value will be written to
 *  @return           An MEA structure (or NULL on any error)
 */
char *
vrna_MEA_from_plist(vrna_ep_t   *plist,
                    const char  *sequence,
                    double      gamma,
                    vrna_md_t   *md,
                    float       *mea);


#ifndef VRNA_DISABLE_BACKWARD_COMPATIBILITY

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
 *  @brief Computes a MEA (maximum expected accuracy) structure.
 *
 *  @ingroup  mea_fold
 *
 *  The algorithm maximizes the expected accuracy
 *
 *  @f[ A(S) = \sum_{(i,j) \in S} 2 \gamma p_{ij} + \sum_{i \notin S} p^u_i @f]
 *
 *  Higher values of @f$\gamma@f$ result in more base pairs of lower
 *  probability and thus higher sensitivity. Low values of @f$\gamma@f$ result in structures
 *  containing only highly likely pairs (high specificity).
 *  The code of the MEA function also demonstrates the use of sparse dynamic
 *  programming scheme to reduce the time and memory complexity of folding.
 *
 *  @deprecated Use vrna_MEA() or vrna_MEA_from_plist() instead!
 */
DEPRECATED(float
           MEA(plist  *p,
               char   *structure,
               double gamma),
           "Use vrna_MEA() or vrna_MEA_from_plist() instead!");


DEPRECATED(float
           MEA_seq(plist            *p,
                   const char       *sequence,
                   char             *structure,
                   double           gamma,
                   vrna_exp_param_t *pf),
           "Use vrna_MEA() or vrna_MEA_from_plist() instead!");


#endif

#endif
