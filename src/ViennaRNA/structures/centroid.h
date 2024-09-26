#ifndef VIENNA_RNA_PACKAGE_STRUCTURES_CENTROID_H
#define VIENNA_RNA_PACKAGE_STRUCTURES_CENTROID_H

#include <ViennaRNA/datastructures/basic.h>
#include <ViennaRNA/fold_compound.h>
#include <ViennaRNA/structures/problist.h>

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
 *  @file   centroid.h
 *  @ingroup subopt_and_representatives
 *  @brief  Centroid structure computation
 */

/**
 *  @brief Get the centroid structure of the ensemble
 *
 *  The centroid is the structure with the minimal average distance to all other structures
 *  \n @f$ <d(S)> = \sum_{(i,j) \in S} (1-p_{ij}) + \sum_{(i,j) \notin S} p_{ij} @f$ \n
 *  Thus, the centroid is simply the structure containing all pairs with @f$p_ij>0.5@f$
 *  The distance of the centroid to the ensemble is written to the memory adressed by @a dist.
 *
 *  @ingroup              centroid_fold
 *  @param[in]    fc      The fold compound data structure
 *  @param[out]   dist    A pointer to the distance variable where the centroid distance will be written to
 *  @return               The centroid structure of the ensemble in dot-bracket notation (@p NULL on error)
 */
char *
vrna_centroid(vrna_fold_compound_t  *fc,
              double                *dist);


/**
 *  @brief Get the centroid structure of the ensemble
 *
 *  This function is a threadsafe replacement for @ref centroid() with a #vrna_ep_t input
 *
 *  The centroid is the structure with the minimal average distance to all other structures
 *  \n @f$ <d(S)> = \sum_{(i,j) \in S} (1-p_{ij}) + \sum_{(i,j) \notin S} p_{ij} @f$ \n
 *  Thus, the centroid is simply the structure containing all pairs with @f$p_ij>0.5@f$
 *  The distance of the centroid to the ensemble is written to the memory adressed by @a dist.
 *
 *  @ingroup            centroid_fold
 *  @param[in]  length  The length of the sequence
 *  @param[out] dist    A pointer to the distance variable where the centroid distance will be written to
 *  @param[in]  pl      A pair list containing base pair probability information about the ensemble
 *  @return             The centroid structure of the ensemble in dot-bracket notation (@p NULL on error)
 */
char *
vrna_centroid_from_plist(int        length,
                         double     *dist,
                         vrna_ep_t  *pl);


/**
 *  @brief Get the centroid structure of the ensemble
 *
 *  This function is a threadsafe replacement for @ref centroid() with a probability array input
 *
 *  The centroid is the structure with the minimal average distance to all other structures
 *  \n @f$ <d(S)> = \sum_{(i,j) \in S} (1-p_{ij}) + \sum_{(i,j) \notin S} p_{ij} @f$ \n
 *  Thus, the centroid is simply the structure containing all pairs with @f$p_ij>0.5@f$
 *  The distance of the centroid to the ensemble is written to the memory adressed by @a dist.
 *
 *  @ingroup              centroid_fold
 *  @param[in]    length  The length of the sequence
 *  @param[out]   dist    A pointer to the distance variable where the centroid distance will be written to
 *  @param[in]    probs   An upper triangular matrix containing base pair probabilities (access via iindx @ref vrna_idx_row_wise() )
 *  @return               The centroid structure of the ensemble in dot-bracket notation (@p NULL on error)
 */
char *
vrna_centroid_from_probs(int        length,
                         double     *dist,
                         FLT_OR_DBL *probs);


#ifndef VRNA_DISABLE_BACKWARD_COMPATIBILITY

/**
 *  @brief Get the centroid structure of the ensemble
 *
 *  @deprecated This function was renamed to vrna_centroid_from_plist()
 */
DEPRECATED(char *get_centroid_struct_pl(int       length,
                                        double    *dist,
                                        vrna_ep_t *pl),
           "Use vrna_centroid_from_plist() instead");

/**
 *  @brief Get the centroid structure of the ensemble
 *
 *  @deprecated This function was renamed to vrna_centroid_from_probs()
 */
DEPRECATED(char *get_centroid_struct_pr(int         length,
                                        double      *dist,
                                        FLT_OR_DBL  *pr),
           "Use vrna_centroid_from_probs() instead");

#endif

#endif
