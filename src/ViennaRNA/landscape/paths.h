#ifndef VIENNA_RNA_PACKAGE_PATHS_H
#define VIENNA_RNA_PACKAGE_PATHS_H

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
 *  @file     ViennaRNA/landscape/paths.h
 *  @ingroup  paths
 *  @brief    API for computing (optimal) (re-)folding paths between secondary structures
 */

/**
 *  @addtogroup   paths
 *  @{
 */

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

/**
 *  @brief Old typename of #vrna_path_s
 *  @deprecated Use #vrna_path_t instead!
 *  @ingroup paths_deprecated
 */
DEPRECATED(typedef struct vrna_path_s path_t,
             "Use vrna_path_t instead!");

#endif

#include <ViennaRNA/fold_compound.h>
#include <ViennaRNA/landscape/move.h>

/**
 *  @brief  Flag to indicate producing a (re-)folding path as list of dot-bracket structures
 *
 *  @see    #vrna_path_t, vrna_path_options_findpath(), vrna_path_direct(), vrna_path_direct_ub()
 */
#define   VRNA_PATH_TYPE_DOT_BRACKET    1U

/**
 *  @brief  Flag to indicate producing a (re-)folding path as list of transition moves
 *
 *  @see    #vrna_path_t, vrna_path_options_findpath(), vrna_path_direct(), vrna_path_direct_ub()
 */
#define   VRNA_PATH_TYPE_MOVES          2U

/**
 *  @brief  An element of a refolding path list
 *
 *  Usually, one has to deal with an array of vrna_path_s, e.g. returned from one of
 *  the refolding-path algorithms.
 *
 *  Since in most cases the length of the list is not known in advance, such lists
 *  have an <em>end-of-list</em> marker, which is either:
 *  - a value of <em>NULL</em> for vrna_path_s::s if vrna_path_s::type = #VRNA_PATH_TYPE_DOT_BRACKET, or
 *  - a vrna_path_s::move with zero in both fields vrna_move_t::pos_5 and vrna_move_t::pos_3 if
 *    vrna_path_s::type = #VRNA_PATH_TYPE_MOVES.
 *
 *  In the following we show an example for how to cover both cases of iteration:
 *  @code{.c}
 *  vrna_path_t *ptr = path; // path was returned from one of the refolding path functions, e.g. vrna_path_direct()
 *
 *  if (ptr) {
 *    if (ptr->type == VRNA_PATH_TYPE_DOT_BRACKET) {
 *      for (; ptr->s; ptr++)
 *        printf("%s [%6.2f]\n", ptr->s, ptr->en);
 *    } else if (ptr->type == VRNA_PATH_TYPE_MOVES) {
 *      for (; ptr->move.pos_5 != 0; ptr++)
 *        printf("move %d:%d, dG = %6.2f\n", ptr->move.pos_5, ptr->move.pos_3, ptr->en);
 *    }
 *  }
 *  @endcode
 *
 *  @see    vrna_path_free()
 */
struct vrna_path_s {
  unsigned int type;  /**< @brief The type of the path element.
                       *  @details A value of #VRNA_PATH_TYPE_DOT_BRACKET indicates that
                       *  vrna_path_s::s consists of the secondary structure in dot-bracket
                       *  notation, and vrna_path_s::en the corresponding free energy.<br>
                       *  On the other hand, if the value is #VRNA_PATH_TYPE_MOVES, vrna_path_s::s
                       *  is <em>NULL</em> and vrna_path_s::move is set to the transition move
                       *  that transforms a previous structure into it's neighbor along the path.
                       *  In this case, the attribute vrna_path_s::en states the change in free
                       *  energy with respect to the structure before application of vrna_path_s::move.
                       */
  double      en;     /**< @brief Free energy of current structure */
  char        *s;     /**< @brief Secondary structure in dot-bracket notation */
  vrna_move_t move;   /**< @brief Move that transforms the previous structure into it's next neighbor along the path */
};


/**
 *  @brief    Release (free) memory occupied by a (re-)folding path
 *
 *  @see      vrna_path_direct(), vrna_path_direct_ub(), vrna_path_findpath(), vrna_path_findpath_ub()
 *
 *  @param    path    The refolding path to be free'd
 */
void
vrna_path_free(vrna_path_t *path);


/**
 *  @brief  Release (free) memory occupied by an options data structure for (re-)folding path implementations
 *
 *  @see    vrna_path_options_findpath(), vrna_path_direct(), vrna_path_direct_ub()
 *
 *  @param  options   The options data structure to be free'd
 */
void
vrna_path_options_free(vrna_path_options_t options);


/** @} */

/**
 *  @addtogroup paths_direct
 *  @{
 */


/**
 *  @brief    Create options data structure for findpath direct (re-)folding path heuristic
 *
 *  This function returns an options data structure that switches the vrna_path_direct()
 *  and vrna_path_direct_ub() API functions to use the <em>findpath</em> @rstinline :cite:p:`flamm:2001` @endrst
 *  heuristic. The parameter @p width specifies the width of the breadth-first search
 *  while the second parameter @p type allows one to set the type of the returned
 *  (re-)folding path.
 *
 *  Currently, the following return types are available:
 *  - A list of dot-bracket structures and corresponding free energy (flag: #VRNA_PATH_TYPE_DOT_BRACKET)
 *  - A list of transition moves and corresponding free energy changes (flag: #VRNA_PATH_TYPE_MOVES)
 *
 *  @see      #VRNA_PATH_TYPE_DOT_BRACKET, #VRNA_PATH_TYPE_MOVES, vrna_path_options_free(),
 *            vrna_path_direct(), vrna_path_direct_ub()
 *
 *  @param    width   Width of the breath-first search strategy
 *  @param    type    Setting that specifies how the return (re-)folding path should be encoded
 *  @returns          An options data structure with settings for the findpath direct path heuristic
 */
vrna_path_options_t
vrna_path_options_findpath(int          width,
                           unsigned int type);


/**
 *  @brief  Determine an optimal direct (re-)folding path between two secondary structures
 *
 *  This is the generic wrapper function to retrieve (an optimal) (re-)folding path
 *  between two secondary structures @p s1 and @p s2. The actual algorithm that
 *  is used to generate the (re-)folding path is determined by the settings specified
 *  in the @p options data structure. This data structure also determines the return
 *  type, which might be either:
 *  - a list of dot-bracket structures with corresponding free energy, or
 *  - a list of transition moves with corresponding free energy change
 *
 *  If the @p options parameter is passed a <em>NULL</em> pointer, this function
 *  defaults to the <em>findpath heuristic</em> @rstinline :cite:p:`flamm:2001` @endrst with a breadth-first
 *  search width of @f$ 10 @f$, and the returned path consists of dot-bracket structures
 *  with corresponding free energies.
 *
 *  @see    vrna_path_direct_ub(), vrna_path_options_findpath(), vrna_path_options_free(),
 *          vrna_path_free()
 *
 *  @param  fc      The #vrna_fold_compound_t with precomputed sequence encoding and model details
 *  @param  s1      The start structure in dot-bracket notation
 *  @param  s2      The target structure in dot-bracket notation
 *  @param  options An options data structure that specifies the path heuristic and corresponding settings (maybe <em>NULL</em>)
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
 *  This function is similar to vrna_path_direct(), but allows to specify an <em>upper-bound</em>
 *  for the saddle point energy. The underlying algorithms will stop determining an (optimal)
 *  (re-)folding path, if none can be found that has a saddle point below the specified
 *  upper-bound threshold @p maxE.
 *
 *  @warning  The argument @p maxE enables one to specify an upper bound, or maximum free
 *            energy for the saddle point between the two input structures. If no path
 *            with @f$E_{saddle} < E_{max}@f$ is found, the function simply returns @em NULL
 *
 *  @see    vrna_path_direct_ub(), vrna_path_options_findpath(), vrna_path_options_free(),
 *          vrna_path_free()
 *
 *  @param  fc      The #vrna_fold_compound_t with precomputed sequence encoding and model details
 *  @param  s1      The start structure in dot-bracket notation
 *  @param  s2      The target structure in dot-bracket notation
 *  @param  maxE    Upper bound for the saddle point along the (re-)folding path
 *  @param  options An options data structure that specifies the path heuristic and corresponding settings (maybe <em>NULL</em>)
 *  @returns        An optimal (re-)folding path between the two input structures
 */
vrna_path_t *
vrna_path_direct_ub(vrna_fold_compound_t  *fc,
                    const char            *s1,
                    const char            *s2,
                    int                   maxE,
                    vrna_path_options_t   options);


/** @} */

#endif
