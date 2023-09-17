#ifndef VIENNA_RNA_PACKAGE_FIND_PATH_H
#define VIENNA_RNA_PACKAGE_FIND_PATH_H

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
 *  @file     ViennaRNA/landscape/findpath.h
 *  @ingroup  paths
 *  @brief    A breadth-first search heuristic for optimal direct folding paths
 */

/**
 *  @addtogroup   paths_direct
 *  @{
 */

#include <ViennaRNA/fold_compound.h>
#include <ViennaRNA/landscape/paths.h>

/**
 *  @brief Find energy of a saddle point between 2 structures (search only direct path)
 *
 *  This function uses an inplementation of the @em findpath algorithm @rstinline :cite:p:`flamm:2001` @endrst
 *  for near-optimal direct refolding path prediction.
 *
 *  Model details, and energy parameters are used as provided via the parameter 'fc'.
 *  The #vrna_fold_compound_t does not require memory for any DP matrices,
 *  but requires all most basic init values as one would get from a call like this:
 *  @code{.c}
 * fc = vrna_fold_compound(sequence, NULL, VRNA_OPTION_DEFAULT);
 *  @endcode
 *
 *  @see vrna_path_findpath_saddle_ub(), vrna_fold_compound(), #vrna_fold_compound_t, vrna_path_findpath()
 *
 *  @param fc     The #vrna_fold_compound_t with precomputed sequence encoding and model details
 *  @param s1     The start structure in dot-bracket notation
 *  @param s2     The target structure in dot-bracket notation
 *  @param width  A number specifying how many strutures are being kept at each step during the search
 *  @returns      The saddle energy in 10cal/mol
 */
int
vrna_path_findpath_saddle(vrna_fold_compound_t  *fc,
                          const char            *s1,
                          const char            *s2,
                          int                   width);


/**
 *  @brief Find energy of a saddle point between 2 structures (search only direct path)
 *
 *  This function uses an inplementation of the @em findpath algorithm @rstinline :cite:p:`flamm:2001` @endrst
 *  for near-optimal direct refolding path prediction.
 *
 *  Model details, and energy parameters are used as provided via the parameter 'fc'.
 *  The #vrna_fold_compound_t does not require memory for any DP matrices,
 *  but requires all most basic init values as one would get from a call like this:
 *  @code{.c}
 * fc = vrna_fold_compound(sequence, NULL, VRNA_OPTION_DEFAULT);
 *  @endcode
 *
 *  @warning  The argument @p maxE (@f$E_{max}@f$) enables one to specify an upper bound, or maximum free
 *            energy for the saddle point between the two input structures. If no path
 *            with @f$E_{saddle} < E_{max}@f$ is found, the function simply returns @p maxE
 *
 *  @see  vrna_path_findpath_saddle(), vrna_fold_compound(), #vrna_fold_compound_t, vrna_path_findpath()
 *
 *  @param fc     The #vrna_fold_compound_t with precomputed sequence encoding and model details
 *  @param s1 The start structure in dot-bracket notation
 *  @param s2 The target structure in dot-bracket notation
 *  @param width  A number specifying how many strutures are being kept at each step during the search
 *  @param maxE   An upper bound for the saddle point energy in 10cal/mol
 *  @returns      The saddle energy in 10cal/mol
 */
int
vrna_path_findpath_saddle_ub(vrna_fold_compound_t *fc,
                             const char           *s1,
                             const char           *s2,
                             int                  width,
                             int                  maxE);


/**
 *  @brief Find refolding path between 2 structures (search only direct path)
 *
 *  This function uses an inplementation of the @em findpath algorithm @rstinline :cite:p:`flamm:2001` @endrst
 *  for near-optimal direct refolding path prediction.
 *
 *  Model details, and energy parameters are used as provided via the parameter 'fc'.
 *  The #vrna_fold_compound_t does not require memory for any DP matrices,
 *  but requires all most basic init values as one would get from a call like this:
 *  @code{.c}
 * fc = vrna_fold_compound(sequence, NULL, VRNA_OPTION_DEFAULT);
 *  @endcode
 *
 *  @see vrna_path_findpath_ub(), vrna_fold_compound(), #vrna_fold_compound_t, vrna_path_findpath_saddle()
 *
 *  @param fc       The #vrna_fold_compound_t with precomputed sequence encoding and model details
 *  @param s1       The start structure in dot-bracket notation
 *  @param s2       The target structure in dot-bracket notation
 *  @param width    A number specifying how many strutures are being kept at each step during the search
 *  @returns        The saddle energy in 10cal/mol
 */
vrna_path_t *
vrna_path_findpath(vrna_fold_compound_t *fc,
                   const char           *s1,
                   const char           *s2,
                   int                  width);


/**
 *  @brief Find refolding path between 2 structures
 *  (search only direct path)
 *
 *  This function uses an inplementation of the @em findpath algorithm @rstinline :cite:p:`flamm:2001` @endrst
 *  for near-optimal direct refolding path prediction.
 *
 *  Model details, and energy parameters are used as provided via the parameter 'fc'.
 *  The #vrna_fold_compound_t does not require memory for any DP matrices,
 *  but requires all most basic init values as one would get from a call like this:
 *  @code{.c}
 * fc = vrna_fold_compound(sequence, NULL, VRNA_OPTION_DEFAULT);
 *  @endcode
 *
 *  @warning  The argument @p maxE enables one to specify an upper bound, or maximum free
 *            energy for the saddle point between the two input structures. If no path
 *            with @f$E_{saddle} < E_{max}@f$ is found, the function simply returns @em NULL
 *
 *  @see vrna_path_findpath(), vrna_fold_compound(), #vrna_fold_compound_t, vrna_path_findpath_saddle()
 *
 *  @param fc       The #vrna_fold_compound_t with precomputed sequence encoding and model details
 *  @param s1       The start structure in dot-bracket notation
 *  @param s2       The target structure in dot-bracket notation
 *  @param width    A number specifying how many strutures are being kept at each step during the search
 *  @param maxE     An upper bound for the saddle point energy in 10cal/mol
 *  @returns        The saddle energy in 10cal/mol
 */
vrna_path_t *
vrna_path_findpath_ub(vrna_fold_compound_t  *fc,
                      const char            *s1,
                      const char            *s2,
                      int                   width,
                      int                   maxE);


#ifndef VRNA_DISABLE_BACKWARD_COMPATIBILITY

/**
 *  @brief Find energy of a saddle point between 2 structures
 *  (search only direct path)
 *
 *  @deprecated Use vrna_path_findpath_saddle() instead!
 *
 *  @ingroup paths_deprecated
 *
 *  @param seq RNA sequence
 *  @param s1 A pointer to the character array where the first
 *         secondary structure in dot-bracket notation will be written to
 *  @param s2 A pointer to the character array where the second
 *         secondary structure in dot-bracket notation will be written to
 *  @param width integer how many strutures are being kept during the search
 *  @returns the saddle energy in 10cal/mol
 */
DEPRECATED(int
           find_saddle(const char *seq,
                       const char *s1,
                       const char *s2,
                       int        width),
           "Use vrna_path_findpath_saddle() instead!");


/**
 *  @brief Free memory allocated by get_path() function
 *
 *  @deprecated Use vrna_path_free() instead!
 *
 *  @ingroup paths_deprecated
 *
 *  @param path pointer to memory to be freed
 */
DEPRECATED(void
           free_path(vrna_path_t *path),
           "Use vrna_path_free() instead!");


/**
 *  @brief Find refolding path between 2 structures
 *  (search only direct path)
 *
 *  @deprecated Use vrna_path_findpath() instead!
 *
 *  @ingroup paths_deprecated
 *
 *  @param seq RNA sequence
 *  @param s1 A pointer to the character array where the first
 *         secondary structure in dot-bracket notation will be written to
 *  @param s2 A pointer to the character array where the second
 *         secondary structure in dot-bracket notation will be written to
 *  @param width integer how many strutures are being kept during the search
 *  @returns direct refolding path between two structures
 */
DEPRECATED(vrna_path_t *
           get_path(const char  *seq,
                    const char  *s1,
                    const char  *s2,
                    int         width),
           "Use vrna_path_findpath() instead!");


#endif

/**
 *  @}
 */

#endif
