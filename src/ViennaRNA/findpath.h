#ifndef VIENNA_RNA_PACKAGE_FIND_PATH_H
#define VIENNA_RNA_PACKAGE_FIND_PATH_H

/**
 *  @file     findpath.h
 *  @ingroup  paths
 *  @brief    A breadth-first search heuristic for optimal direct folding paths
 */

/**
 *  @addtogroup   direct_paths
 *  @{
 *
 *  @brief Heuristics to explore direct, optimal (re-)folding paths between two secondary structures
 *
 */

/* below are several convenience typedef's we use throughout the ViennaRNA library */

/**
 *  @brief Typename for the refolding path data structure #vrna_path_s
 */
typedef struct vrna_path_s vrna_path_t;


#ifndef VRNA_DISABLE_BACKWARD_COMPATIBILITY

/* the following typedefs are for backward compatibility only */

/**
 *  @brief Old typename of #vrna_path_s
 *  @deprecated Use #vrna_path_t instead!
 */
typedef struct vrna_path_s path_t;

#endif

#include <ViennaRNA/datastructures/basic.h>
#include <ViennaRNA/fold_compound.h>

/**
 *  @brief  An element of a refolding path list
 *  @see    vrna_path_findpath()
 */
struct vrna_path_s {
  double  en; /**<  @brief  Free energy of current structure */
  char    *s; /**<  @brief  Secondary structure in dot-bracket notation */
};


/**
 *  @brief Find energy of a saddle point between 2 structures (search only direct path)
 *
 *  This function uses an inplementation of the @em findpath algorithm @cite flamm:2001
 *  for near-optimal direct refolding path prediction.
 *
 *  Model details, and energy parameters are used as provided via the parameter 'vc'.
 *  The #vrna_fold_compound_t does not require memory for any DP matrices,
 *  but requires all most basic init values as one would get from a call like this:
 *  @code{.c}
 * vc = vrna_fold_compound(sequence, NULL, VRNA_OPTION_DEFAULT);
 *  @endcode
 *
 *  @see vrna_path_findpath_saddle_ub(), vrna_fold_compound(), #vrna_fold_compound_t, vrna_path_findpath()
 *
 *  @param vc     The #vrna_fold_compound_t with precomputed sequence encoding and model details
 *  @param s1     The start structure in dot-bracket notation
 *  @param s2     The target structure in dot-bracket notation
 *  @param width  A number specifying how many strutures are being kept at each step during the search
 *  @returns      The saddle energy in 10cal/mol
 */
int
vrna_path_findpath_saddle(vrna_fold_compound_t  *vc,
                          const char            *s1,
                          const char            *s2,
                          int                   width);


/**
 *  @brief Find energy of a saddle point between 2 structures (search only direct path)
 *
 *  This function uses an inplementation of the @em findpath algorithm @cite flamm:2001
 *  for near-optimal direct refolding path prediction.
 *
 *  Model details, and energy parameters are used as provided via the parameter 'vc'.
 *  The #vrna_fold_compound_t does not require memory for any DP matrices,
 *  but requires all most basic init values as one would get from a call like this:
 *  @code{.c}
 * vc = vrna_fold_compound(sequence, NULL, VRNA_OPTION_DEFAULT);
 *  @endcode
 *
 *  @warning  The argument @p maxE (@f$E_{max}@f$) enables one to specify an upper bound, or maximum free
 *            energy for the saddle point between the two input structures. If no path
 *            with @f$E_{saddle} < E_{max}@f$ is found, the function simply returns @p maxE
 *
 *  @see  vrna_path_findpath_saddle(), vrna_fold_compound(), #vrna_fold_compound_t, vrna_path_findpath()
 *
 *  @param vc     The #vrna_fold_compound_t with precomputed sequence encoding and model details
 *  @param s1 The start structure in dot-bracket notation
 *  @param s2 The target structure in dot-bracket notation
 *  @param width  A number specifying how many strutures are being kept at each step during the search
 *  @param maxE   An upper bound for the saddle point energy in 10cal/mol
 *  @returns      The saddle energy in 10cal/mol
 */
int
vrna_path_findpath_saddle_ub(vrna_fold_compound_t *vc,
                             const char           *s1,
                             const char           *s2,
                             int                  width,
                             int                  maxE);


/**
 *  @brief Find refolding path between 2 structures (search only direct path)
 *
 *  This function uses an inplementation of the @em findpath algorithm @cite flamm:2001
 *  for near-optimal direct refolding path prediction.
 *
 *  Model details, and energy parameters are used as provided via the parameter 'vc'.
 *  The #vrna_fold_compound_t does not require memory for any DP matrices,
 *  but requires all most basic init values as one would get from a call like this:
 *  @code{.c}
 * vc = vrna_fold_compound(sequence, NULL, VRNA_OPTION_DEFAULT);
 *  @endcode
 *
 *  @see vrna_path_findpath_ub(), vrna_fold_compound(), #vrna_fold_compound_t, vrna_path_findpath_saddle()
 *
 *  @param vc       The #vrna_fold_compound_t with precomputed sequence encoding and model details
 *  @param s1       The start structure in dot-bracket notation
 *  @param s2       The target structure in dot-bracket notation
 *  @param width    A number specifying how many strutures are being kept at each step during the search
 *  @returns        The saddle energy in 10cal/mol
 */
vrna_path_t *
vrna_path_findpath(vrna_fold_compound_t *vc,
                   const char           *s1,
                   const char           *s2,
                   int                  width);


/**
 *  @brief Find refolding path between 2 structures
 *  (search only direct path)
 *
 *  This function uses an inplementation of the @em findpath algorithm @cite flamm:2001
 *  for near-optimal direct refolding path prediction.
 *
 *  Model details, and energy parameters are used as provided via the parameter 'vc'.
 *  The #vrna_fold_compound_t does not require memory for any DP matrices,
 *  but requires all most basic init values as one would get from a call like this:
 *  @code{.c}
 * vc = vrna_fold_compound(sequence, NULL, VRNA_OPTION_DEFAULT);
 *  @endcode
 *
 *  @warning  The argument @p maxE enables one to specify an upper bound, or maximum free
 *            energy for the saddle point between the two input structures. If no path
 *            with @f$E_{saddle} < E_{max}@f$ is found, the function simply returns @em NULL
 *
 *  @see vrna_path_findpath(), vrna_fold_compound(), #vrna_fold_compound_t, vrna_path_findpath_saddle()
 *
 *  @param vc       The #vrna_fold_compound_t with precomputed sequence encoding and model details
 *  @param s1       The start structure in dot-bracket notation
 *  @param s2       The target structure in dot-bracket notation
 *  @param width    A number specifying how many strutures are being kept at each step during the search
 *  @param maxE     An upper bound for the saddle point energy in 10cal/mol
 *  @returns        The saddle energy in 10cal/mol
 */
vrna_path_t *
vrna_path_findpath_ub(vrna_fold_compound_t  *vc,
                      const char            *s1,
                      const char            *s2,
                      int                   width,
                      int                   maxE);


#ifndef VRNA_DISABLE_BACKWARD_COMPATIBILITY

/**
 *  \brief Find energy of a saddle point between 2 structures
 *  (search only direct path)
 *
 *  \param seq RNA sequence
 *  \param s1 A pointer to the character array where the first
 *         secondary structure in dot-bracket notation will be written to
 *  \param s2 A pointer to the character array where the second
 *         secondary structure in dot-bracket notation will be written to
 *  \param width integer how many strutures are being kept during the search
 *  \returns the saddle energy in 10cal/mol
 */
int
find_saddle(const char  *seq,
            const char  *s1,
            const char  *s2,
            int         width);


/**
 *  \brief Free memory allocated by get_path() function
 *
 *  \param path pointer to memory to be freed
 */
void
free_path(vrna_path_t *path);


/**
 *  \brief Find refolding path between 2 structures
 *  (search only direct path)
 *
 *  \param seq RNA sequence
 *  \param s1 A pointer to the character array where the first
 *         secondary structure in dot-bracket notation will be written to
 *  \param s2 A pointer to the character array where the second
 *         secondary structure in dot-bracket notation will be written to
 *  \param width integer how many strutures are being kept during the search
 *  \returns direct refolding path between two structures
 */
vrna_path_t *
get_path(const char *seq,
         const char *s1,
         const char *s2,
         int        width);


#endif

/**
 *  @}
 */

#endif
