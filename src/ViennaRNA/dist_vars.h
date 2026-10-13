#ifndef VIENNA_RNA_PACKAGE_DIST_VARS_H
#define VIENNA_RNA_PACKAGE_DIST_VARS_H

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
 *  @file dist_vars.h
 *  @brief Global variables for Distance-Package
 */

#ifndef VRNA_DISABLE_BACKWARD_COMPATIBILITY

/**
 *  @brief Produce an alignment of the two structures being compared by
 *  tracing the editing path giving the minimum distance.
 * 
 *  set to 1 if you want backtracking
 */
DEPRECATED(extern int   edit_backtrack,
           "This global variable is obsolete");

/**
 *  @brief Contains the two aligned structures after a call to one of the distance
 *  functions with #edit_backtrack set to 1.
 */
DEPRECATED(extern char *aligned_line[4],
           "This global variable is obsolete");

/**
 *  @brief Specify the cost matrix to be used for distance calculations
 * 
 *  if 0, use the default cost matrix (upper matrix in example), otherwise
 *  use Shapiro's costs (lower matrix).
 */
DEPRECATED(extern int  cost_matrix,
           "This global variable is obsolete");

/*  Global type defs for Distance-Package */

/**
 *  @brief Postorder data structure
 */
typedef struct {
                 int  type;
                 int  weight;
                 int  father;
                 int  sons;
                 int  leftmostleaf;
               } Postorder_list;

/**
 *  @brief  Tree data structure
 */
typedef struct {
                 Postorder_list *postorder_list;
                 int            *keyroots;
               } Tree;

/**
 *  @brief  Some other data structure
 */
typedef struct {
                 int    type;
                 int    sign;
                 float  weight;
               } swString;
#endif

#endif
