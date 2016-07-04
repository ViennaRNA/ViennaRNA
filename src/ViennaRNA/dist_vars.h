#ifndef VIENNA_RNA_PACKAGE_DIST_VARS_H
#define VIENNA_RNA_PACKAGE_DIST_VARS_H

/**
 *  @file dist_vars.h
 *  @brief Global variables for Distance-Package
 */

/**
 *  @brief Produce an alignment of the two structures being compared by
 *  tracing the editing path giving the minimum distance.
 * 
 *  set to 1 if you want backtracking
 */
extern int   edit_backtrack;

/**
 *  @brief Contains the two aligned structures after a call to one of the distance
 *  functions with #edit_backtrack set to 1.
 */
extern char *aligned_line[4];

/**
 *  @brief Specify the cost matrix to be used for distance calculations
 * 
 *  if 0, use the default cost matrix (upper matrix in example), otherwise
 *  use Shapiro's costs (lower matrix).
 */
extern int  cost_matrix;

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
