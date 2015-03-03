#ifndef __VIENNA_RNA_PACKAGE_CENTROID_H__
#define __VIENNA_RNA_PACKAGE_CENTROID_H__

#include <ViennaRNA/data_structures.h>

#ifdef __GNUC__
#define DEPRECATED(func) func __attribute__ ((deprecated))
#else
#define DEPRECATED(func) func
#endif

/**
 *  \file   centroid.h
 *  \brief  Centroid structure computation
 */


/**
 *  \brief Get the centroid structure of the ensemble
 * 
 *  This function is a threadsafe replacement for \ref centroid() with a 'plist' input
 * 
 *  The centroid is the structure with the minimal average distance to all other structures
 *  \n \f$ <d(S)> = \sum_{(i,j) \in S} (1-p_{ij}) + \sum_{(i,j) \notin S} p_{ij} \f$ \n
 *  Thus, the centroid is simply the structure containing all pairs with \f$p_ij>0.5\f$
 *  The distance of the centroid to the ensemble is written to the memory adressed by \a dist.
 *
 *  \ingroup            centroid_fold
 *  \param[in]  length  The length of the sequence
 *  \param[out] dist    A pointer to the distance variable where the centroid distance will be written to
 *  \param[in]  pl      A pair list containing base pair probability information about the ensemble
 *  \return             The centroid structure of the ensemble in dot-bracket notation
 */
char  *vrna_get_centroid_struct_pl( int length,
                                    double *dist,
                                    plist *pl);

/**
 *  \brief Get the centroid structure of the ensemble
 *
 *  \deprecated This function was renamed to vrna_get_centroid_struct_pl()
 */
DEPRECATED(char  *get_centroid_struct_pl(int length,
                              double *dist,
                              plist *pl));

/**
 *  \brief Get the centroid structure of the ensemble
 * 
 *  This function is a threadsafe replacement for \ref centroid() with a probability array input
 * 
 *  The centroid is the structure with the minimal average distance to all other structures
 *  \n \f$ <d(S)> = \sum_{(i,j) \in S} (1-p_{ij}) + \sum_{(i,j) \notin S} p_{ij} \f$ \n
 *  Thus, the centroid is simply the structure containing all pairs with \f$p_ij>0.5\f$
 *  The distance of the centroid to the ensemble is written to the memory adressed by \a dist.
 * 
 *  \ingroup              centroid_fold
 *  \param[in]    length  The length of the sequence
 *  \param[out]   dist    A pointer to the distance variable where the centroid distance will be written to
 *  \param[in]    pr      A upper triangular matrix containing base pair probabilities (access via iindx \ref vrna_get_iindx() )
 *  \return               The centroid structure of the ensemble in dot-bracket notation
 */
char  *vrna_get_centroid_struct_pr( int length,
                                    double *dist,
                                    FLT_OR_DBL *pr);

/**
 *  \brief Get the centroid structure of the ensemble
 *
 *  \deprecated This function was renamed to vrna_get_centroid_struct_pr()
 */
DEPRECATED(char  *get_centroid_struct_pr(int length,
                              double *dist,
                              FLT_OR_DBL *pr));


/**
 *  \brief Get the centroid structure of the ensemble
 * 
 *  The centroid is the structure with the minimal average distance to all other structures
 *  \n \f$ <d(S)> = \sum_{(i,j) \in S} (1-p_{ij}) + \sum_{(i,j) \notin S} p_{ij} \f$ \n
 *  Thus, the centroid is simply the structure containing all pairs with \f$p_ij>0.5\f$
 *  The distance of the centroid to the ensemble is written to the memory adressed by \a dist.
 * 
 *  \ingroup              centroid_fold
 *  \param[in]    vc      The fold compound data structure
 *  \param[out]   dist    A pointer to the distance variable where the centroid distance will be written to
 *  \return               The centroid structure of the ensemble in dot-bracket notation
 */
char *vrna_get_centroid_struct( vrna_fold_compound *vc,
                                double *dist);

#endif
