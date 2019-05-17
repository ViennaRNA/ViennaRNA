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


/**
 *  @brief  Options data structure for (re-)folding path implementations
 *  @ingroup paths
 */
typedef struct vrna_path_options_s *vrna_path_options_t;


#ifndef VRNA_DISABLE_BACKWARD_COMPATIBILITY

/* the following typedefs are for backward compatibility only */

/**
 *  @brief Old typename of #vrna_path_s
 *  @deprecated Use #vrna_path_t instead!
 */
DEPRECATED(typedef struct vrna_path_s path_t, "Use vrna_path_t instead!");

#endif

#include <ViennaRNA/datastructures/basic.h>
#include <ViennaRNA/fold_compound.h>
#include <ViennaRNA/landscape/move.h>

#define   VRNA_PATH_TYPE_DOT_BRACKET    1U
#define   VRNA_PATH_TYPE_MOVES          2U

#define   VRNA_PATH_DIRECT_FINDPATH     4U
#define   VRNA_PATH_DIRECT_DEFAULT      (VRNA_PATH_DIRECT_FINDPATH | VRNA_PATH_TYPE_DOT_BRACKET)


/**
 *  @brief  An element of a refolding path list
 *  @see    vrna_path_findpath()
 *  @ingroup  paths
 */
struct vrna_path_s {
  unsigned int  type;

  double        en; /**<  @brief  Free energy of current structure */
  char          *s; /**<  @brief  Secondary structure in dot-bracket notation */

  vrna_move_t   move;
};


/**
 *  @brief    Release (free) memory occupied by a (re-)folding path
 *  @ingroup  paths
 */
void
vrna_path_free(vrna_path_t *path);


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


/**
 *  @brief    Create options data structure for findpath direct (re-)folding path heuristic
 *
 *  @see      #VRNA_PATH_TYPE_DOT_BRACKET, #VRNA_PATH_TYPE_MOVES, vrna_path_options_free(),
 *            vrna_path_direct(), vrna_path_direct_ub()
 *
 *  @ingroup  paths
 *
 *  @param    width   Width of the breath-first search strategy
 *  @param    type    Setting that specifies how the return (re-)folding path should be encoded
 *  @returns          An options data structure with settings for the findpath direct path heuristic
 */
vrna_path_options_t
vrna_path_options_findpath(int          width,
                           unsigned int type);


/**
 *  @brief  Release (free) memory occupied by an options data structure for (re-)folding path implementations
 *  @see    vrna_path_options_findpath(), vrna_path_direct(), vrna_path_direct_ub()
 *
 *  @ingroup  paths
 *
 *  @param  options   The options data structure to be free'd
 */
void
vrna_path_options_free(vrna_path_options_t options);


/**
 *  @brief  Determine an optimal direct (re-)folding path between two secondary structures
 *
 *  @see    vrna_path_direct_ub(), vrna_path_options_findpath(), vrna_path_options_free(),
 *          vrna_path_free()
 *
 *  @param  fc      The #vrna_fold_compound_t with precomputed sequence encoding and model details
 *  @param  s1      The start structure in dot-bracket notation
 *  @param  s2      The target structure in dot-bracket notation
 *  @param  options An options data structure that specifies the path heuristic and corresponding settings
 *  @returns        An optimal (re-)folding path between the two input structures
 */
vrna_path_t *
vrna_path_direct(vrna_fold_compound_t *fc,
                 const char           *s1,
                 const char           *s2,
                 vrna_path_options_t  options);


/**
 *  @brief  Determine an optimal direct (re-)folding path between two secondary structures
 *
 *  @see    vrna_path_direct_ub(), vrna_path_options_findpath(), vrna_path_options_free(),
 *          vrna_path_free()
 *
 *  @param  fc      The #vrna_fold_compound_t with precomputed sequence encoding and model details
 *  @param  s1      The start structure in dot-bracket notation
 *  @param  s2      The target structure in dot-bracket notation
 *  @param  maxE    Upper bound for the saddle point along the (re-)folding path
 *  @param  options An options data structure that specifies the path heuristic and corresponding settings
 *  @returns        An optimal (re-)folding path between the two input structures
 */
vrna_path_t *
vrna_path_direct_ub(vrna_fold_compound_t  *fc,
                    const char            *s1,
                    const char            *s2,
                    int                   maxE,
                    vrna_path_options_t   options);


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
           get_path(const char *seq,
                    const char *s1,
                    const char *s2,
                    int width),
           "Use vrna_path_findpath() instead!");


#endif

/**
 *  @}
 */

#endif
