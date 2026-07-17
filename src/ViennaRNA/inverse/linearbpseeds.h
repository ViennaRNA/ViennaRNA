#ifndef VIENNARNA_DLIB_WRAPPER_H
#define VIENNARNA_DLIB_WRAPPER_H

#ifdef __cplusplus
extern "C" {
#endif

/**
 *  @file     ViennaRNA/inverse/linearbpseeds.h
 *  @ingroup  inverse_fold
 *  @brief    Sequence seed sampling with the LinearBP approach 
 */

  
/**
 *  @addtogroup inverse_fold
 *  @{

/**
 *  @brief  The options data structure for vrna_inverse_seed_linearbp()
 *
 *  @see vrna_inverse_seed_linearbp()
 */
typedef struct {
    unsigned int  biseparable;  /**< @brief Flag indicating if biseparanle mode is switched on (A nd Cs) in unpaired positions */
    unsigned int  maxmodulo;    /**< @brief The maximal modulus to be considered for m/mA/mC */
    double        GCweight;     /**< @brief The weight of the GC pairs (default is 0.0) */
} vrna_inverse_linearbp_opt_t;


/**
 *  @brief  Generate a set of sequence seeds compatible with a given structure according to the LinearBP design algorithm
 *
 *  @param  structure   A secondary structure in dot-bracket format
 *  @param  num_samples The target number of solutions
 *  @param  options     A pointer to the options data structure to modifiy the behavior of the algorithm
 *  @return             A pointer to an array of sequences (end of list indicated by @b NULL entry) or @b NULL on error
*/

char **
vrna_inverse_seed_linearbp(const char                   *structure,
                           size_t                       num_samples,
                           vrna_inverse_linearbp_opt_t  *options);

/**
 *  @}
 */


#ifdef __cplusplus
}
#endif

#endif
