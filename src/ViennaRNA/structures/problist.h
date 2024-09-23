#ifndef VIENNA_RNA_PACKAGE_STRUCTURES_PROBLIST_H
#define VIENNA_RNA_PACKAGE_STRUCTURES_PROBLIST_H

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
 *  @file     ViennaRNA/structures/problist.h
 *  @ingroup  struct_utils
 *  @brief    Various utility- and helper-functions for secondary structure parsing, converting, etc.
 */

/**
 *  @addtogroup struct_utils_plist
 *  @{
 */

/**
 *  @brief Convenience typedef for data structure #vrna_elem_prob_s
 *  @ingroup  struct_utils_plist
 */
typedef struct vrna_elem_prob_s vrna_ep_t;


/**
 *  @brief  A Base Pair element
 */
#define VRNA_PLIST_TYPE_BASEPAIR      0


/**
 *  @brief  A G-Quadruplex element
 */
#define VRNA_PLIST_TYPE_GQUAD         1


/**
 *  @brief  A Hairpin loop motif element
 */
#define VRNA_PLIST_TYPE_H_MOTIF       2


/**
 *  @brief  An Internal loop motif element
 */
#define VRNA_PLIST_TYPE_I_MOTIF       3


/**
 *  @brief  An Unstructured Domain motif element
 */
#define VRNA_PLIST_TYPE_UD_MOTIF      4


/**
 *  @brief  A Base Pair stack element
 */
#define VRNA_PLIST_TYPE_STACK         5


/**
 *  @brief  An unpaired base
 */
#define VRNA_PLIST_TYPE_UNPAIRED      6


/**
 *  @brief  One pair of a base triplet
 */
#define VRNA_PLIST_TYPE_TRIPLE        7


/**
 *  @brief  Data structure representing a single entry of an element probability list
 *          (e.g. list of pair probabilities)
 *
 *  @see vrna_plist(), vrna_plist_from_probs(), vrna_db_from_plist(),
 *
 *  #VRNA_PLIST_TYPE_BASEPAIR, #VRNA_PLIST_TYPE_GQUAD, #VRNA_PLIST_TYPE_H_MOTIF, #VRNA_PLIST_TYPE_I_MOTIF,
 *  #VRNA_PLIST_TYPE_UD_MOTIF, #VRNA_PLIST_TYPE_STACK
 */
struct vrna_elem_prob_s {
  int   i;    /**<  @brief  Start position (usually 5' nucleotide that starts the element, e.g. base pair) */
  int   j;    /**<  @brief  End position (usually 3' nucleotide that ends the element, e.g. base pair) */
  float p;    /**<  @brief  Probability of the element */
  int   type; /**<  @brief  Type of the element */
};

/**
 *  @brief Create a #vrna_ep_t from a dot-bracket string
 *
 *  The dot-bracket string is parsed and for each base pair an
 *  entry in the plist is created. The probability of each pair in
 *  the list is set by a function parameter.
 *
 *  The end of the plist is marked by sequence positions i as well as j
 *  equal to 0. This condition should be used to stop looping over its
 *  entries
 *
 *  @param struc  The secondary structure in dot-bracket notation
 *  @param pr     The probability for each base pair used in the plist
 *  @return       The plist array
 */
vrna_ep_t *vrna_plist(const char  *struc,
                      float       pr);


/**
 *  @brief Create a #vrna_ep_t from base pair probability matrix
 *
 *  The probability matrix provided via the #vrna_fold_compound_t is parsed
 *  and all pair probabilities above the given threshold are used to create
 *  an entry in the plist
 *
 *  The end of the plist is marked by sequence positions i as well as j
 *  equal to 0. This condition should be used to stop looping over its
 *  entries
 *
 *  @ingroup              part_func_global
 *  @param[in]  fc        The fold compound
 *  @param[in]  cut_off   The cutoff value
 *  @return               A pointer to the plist that is to be created
 */
vrna_ep_t *vrna_plist_from_probs(vrna_fold_compound_t *fc,
                                 double               cut_off);


int
vrna_plist_append(vrna_ep_t       **target,
                  const vrna_ep_t *list);


/* End pair list interface */
/** @} */


/**
 *  @addtogroup gquad_parse
 *  @{
 */
plist *
get_plist_gquad_from_db(const char  *structure,
                        float       pr);

/** @} */


#ifndef VRNA_DISABLE_BACKWARD_COMPATIBILITY

/*###########################################*/
/*# deprecated functions below              #*/
/*###########################################*/

/**
 *  @brief Create a #vrna_ep_t from a dot-bracket string
 *
 *  The dot-bracket string is parsed and for each base pair an
 *  entry in the plist is created. The probability of each pair in
 *  the list is set by a function parameter.
 *
 *  The end of the plist is marked by sequence positions i as well as j
 *  equal to 0. This condition should be used to stop looping over its
 *  entries
 *
 *  @deprecated   Use vrna_plist() instead
 *
 *  @ingroup part_func_global_deprecated
 *
 *  @param pl     A pointer to the #vrna_ep_t that is to be created
 *  @param struc  The secondary structure in dot-bracket notation
 *  @param pr     The probability for each base pair
 */
DEPRECATED(void assign_plist_from_db(vrna_ep_t  **pl,
                                     const char *struc,
                                     float      pr),
           "Use vrna_plist() instead");

/**
 *  @brief Create a vrna_ep_t from a probability matrix
 *
 *  The probability matrix given is parsed and all pair probabilities above
 *  the given threshold are used to create an entry in the plist
 *
 *  The end of the plist is marked by sequence positions i as well as j
 *  equal to 0. This condition should be used to stop looping over its
 *  entries
 *
 *  @note This function is threadsafe
 *  @deprecated Use vrna_plist_from_probs() instead!
 *
 *  @ingroup part_func_global_deprecated
 *
 *  @param[out] pl      A pointer to the vrna_ep_t that is to be created
 *  @param[in]  probs   The probability matrix used for creating the plist
 *  @param[in]  length  The length of the RNA sequence
 *  @param[in]  cutoff  The cutoff value
 */
DEPRECATED(void  assign_plist_from_pr(vrna_ep_t   **pl,
                                      FLT_OR_DBL  *probs,
                                      int         length,
                                      double      cutoff),
           "Use vrna_plist_from_probs() instead");

#endif

#endif
