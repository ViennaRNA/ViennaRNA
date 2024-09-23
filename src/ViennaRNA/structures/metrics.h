#ifndef VIENNA_RNA_PACKAGE_STRUCTURES_METRICS_H
#define VIENNA_RNA_PACKAGE_STRUCTURES_METRICS_H

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
 *  @file     ViennaRNA/structures/metrics.h
 *  @ingroup  struct_utils
 *  @brief    Various utility- and helper-functions for secondary structure parsing, converting, etc.
 */

/**
 *  @addtogroup struct_utils_metrics
 *  @{
 */

/**
 *  @brief Compute the "base pair" distance between two pair tables pt1 and pt2 of secondary structures.
 *
 *  The pair tables should have the same length.
 *  dist = number of base pairs in one structure but not in the other
 *  same as edit distance with open-pair close-pair as move-set
 *
 *  @see vrna_bp_distance()
 *
 *  @param pt1   First structure in dot-bracket notation
 *  @param pt2   Second structure in dot-bracket notation
 *  @return       The base pair distance between pt1 and pt2
 */
int
vrna_bp_distance_pt(const short *pt1,
                    const short *pt2);

/**
 *  @brief Compute the "base pair" distance between two secondary structures s1 and s2.
 *
 *  This is a wrapper around @b vrna_bp_distance_pt().
 *  The sequences should have the same length.
 *  dist = number of base pairs in one structure but not in the other
 *  same as edit distance with open-pair close-pair as move-set
 *
 *  @see vrna_bp_distance_pt()
 *
 *  @param str1   First structure in dot-bracket notation
 *  @param str2   Second structure in dot-bracket notation
 *  @return       The base pair distance between str1 and str2
 */
int
vrna_bp_distance(const char *str1,
                 const char *str2);


double
vrna_dist_mountain(const char   *str1,
                   const char   *str2,
                   unsigned int p);


/* End metrics interface */
/** @} */


#ifndef VRNA_DISABLE_BACKWARD_COMPATIBILITY

/*###########################################*/
/*# deprecated functions below              #*/
/*###########################################*/


/**
 *  @brief Compute the "base pair" distance between two secondary structures s1 and s2.
 *
 *  The sequences should have the same length.
 *  dist = number of base pairs in one structure but not in the other
 *  same as edit distance with open-pair close-pair as move-set
 *
 *  @deprecated   Use vrna_bp_distance instead
 *  @ingroup        struct_utils_deprecated
 *  @param str1   First structure in dot-bracket notation
 *  @param str2   Second structure in dot-bracket notation
 *  @return       The base pair distance between str1 and str2
 */
DEPRECATED(int bp_distance(const char *str1,
                           const char *str2),
           "Use vrna_bp_distance() instead");

#endif

#endif
