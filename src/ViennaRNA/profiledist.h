#ifndef VIENNA_RNA_PACKAGE_PROFILEDIST_H
#define VIENNA_RNA_PACKAGE_PROFILEDIST_H

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

#include <ViennaRNA/datastructures/basic.h>

/** \file profiledist.h  */

#ifndef VRNA_DISABLE_BACKWARD_COMPATIBILITY

/**
 *  \brief Align the 2 probability profiles T1, T2\n
 *
 *  This is like a Needleman-Wunsch alignment,
 *  we should really use affine gap-costs ala Gotoh
 */
float profile_edit_distance(const float *T1,
                            const float *T2);


/**
 *  \brief condense pair probability matrix into a vector containing probabilities
 *  for unpaired, upstream paired and downstream paired.
 *
 *  This resulting probability profile is used as input for profile_edit_distance
 *
 *  \param bppm   A pointer to the base pair probability matrix
 *  \param length The length of the sequence
 *  \returns      The bp profile
 */
float *Make_bp_profile_bppm(FLT_OR_DBL  *bppm,
                            int         length);


/**
 *  \brief print string representation of probability profile
 */
void  print_bppm(const float *T);


/**
 *  \brief free space allocated in Make_bp_profile
 *
 *  Backward compatibility only. You can just use plain free()
 */
void  free_profile(float *T);


/**
 *  \note This function is NOT threadsafe
 *
 *  \see Make_bp_profile_bppm()
 *
 *  \deprecated This function is deprecated and will be removed soon! See \ref Make_bp_profile_bppm() for a replacement
 *
 */
DEPRECATED(float *Make_bp_profile(int length),
"Use Make_bp_profile_bppm() instead");

#endif

#endif
