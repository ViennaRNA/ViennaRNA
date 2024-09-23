#ifndef VIENNA_RNA_PACKAGE_STRUCTURES_UTILS_H
#define VIENNA_RNA_PACKAGE_STRUCTURES_UTILS_H

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
 *  @file     ViennaRNA/structures/utils.h
 *  @ingroup  struct_utils
 *  @brief    Various utility- and helper-functions for secondary structure parsing, converting, etc.
 */

/**
 *  @addtogroup   struct_utils
 *  @{
 */

/**
 *  @brief Make a reference base pair count matrix
 *
 *  Get an upper triangular matrix containing the number of basepairs of a reference
 *  structure for each interval [i,j] with i<j. Access it via iindx!!!
 */
unsigned int *
vrna_refBPcnt_matrix(const short  *reference_pt,
                     unsigned int turn);


/**
 *  @brief Make a reference base pair distance matrix
 *
 *  Get an upper triangular matrix containing the base pair distance of two
 *  reference structures for each interval [i,j] with i<j. Access it via iindx!!!
 *
 */
unsigned int *
vrna_refBPdist_matrix(const short   *pt1,
                      const short   *pt2,
                      unsigned int  turn);


/** @} */

#ifndef VRNA_DISABLE_BACKWARD_COMPATIBILITY

/*###########################################*/
/*# deprecated functions below              #*/
/*###########################################*/

/**
 *  @brief Make a reference base pair count matrix
 *
 *  Get an upper triangular matrix containing the number of basepairs of a reference
 *  structure for each interval [i,j] with i<j. Access it via iindx!!!
 *
 *  @deprecated Use vrna_refBPcnt_matrix() instead
 *  @ingroup        struct_utils_deprecated
 */
DEPRECATED(unsigned int *make_referenceBP_array(short         *reference_pt,
                                                unsigned int  turn),
           "Use vrna_refBPcnt_matrix() instead");

/**
 *  @brief Make a reference base pair distance matrix
 *
 *  Get an upper triangular matrix containing the base pair distance of two
 *  reference structures for each interval [i,j] with i<j. Access it via iindx!!!
 *
 *  @deprecated Use vrna_refBPdist_matrix() instead
 *  @ingroup        struct_utils_deprecated
 */
DEPRECATED(unsigned int *compute_BPdifferences(short        *pt1,
                                               short        *pt2,
                                               unsigned int turn),
           "Use vrna_refBPdist_matrix() instead");

#endif

#endif
