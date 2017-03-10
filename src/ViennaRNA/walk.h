#ifndef VIENNA_RNA_PACKAGE_WALK_H
#define VIENNA_RNA_PACKAGE_WALK_H

/**
 *  @file     walk.h
 *  @ingroup  paths
 *  @brief    Methods to generate particular paths such as gradient or random walks through the energy landscape of an RNA sequence
 *
 */

#include <ViennaRNA/data_structures.h>
#include <ViennaRNA/neighbor.h>

/**
 *  @addtogroup paths
 *
 *  @{
 *  @ingroup  paths
 */

/**
 * @brief Option flag to request a steepest descent / gradient path
 * @see   vrna_path()
 */
#define  VRNA_PATH_STEEPEST_DESCENT 128

/**
 * @brief Option flag to request a random walk path
 * @see   vrna_path()
 */
#define  VRNA_PATH_RANDOM           256

/**
 * @brief Option flag to omit returning the transition path
 * @see   vrna_path(), vrna_path_gradient(), vrna_path_random()
 */
#define  VRNA_PATH_NO_TRANSITION_OUTPUT          512

/**
 * @brief Option flag to request defaults (steepest descent / default move set)
 * @see   vrna_path(), #VRNA_PATH_STEEPEST_DESCENT, #VRNA_MOVESET_DEFAULT
 */

#define VRNA_PATH_DEFAULT   (VRNA_PATH_STEEPEST_DESCENT | VRNA_MOVESET_DEFAULT)

/**
 *  @brief Compute a path, store the final structure, and return a list of transition moves
 *  from the start to the final structure.
 *
 *  This function computes, given a start structure in pair table format, a transition path,
 *  updates the pair table to the final structure of the path. Finally, if not requested otherwise
 *  by using the #VRNA_PATH_NO_TRANSITION_OUTPUT flag in the @p options field, this function returns a list
 *  of individual transitions that lead from the start to the final
 *  structure if requested.
 *
 *  The currently available transition paths are
 *  - Steepest Descent / Gradient walk (flag: #VRNA_PATH_STEEPEST_DESCENT)
 *  - Random walk (flag: #VRNA_PATH_RANDOM)
 *
 *  The type of transitions must be set through the @p options parameter
 *
 *  @note   Since the result is written to the input structure you may want to use
 *          vrna_ptable_copy() before calling this function to keep the initial structure
 *
 *  @see    vrna_path_gradient(), vrna_path_random(), vrna_ptable(), vrna_ptable_copy(), vrna_fold_compound()
 *          #VRNA_PATH_STEEPEST_DESCENT, #VRNA_PATH_RANDOM, #VRNA_MOVESET_DEFAULT, #VRNA_MOVESET_SHIFT,
 *          #VRNA_PATH_NO_TRANSITION_OUTPUT
 *
 *  @param[in]      vc      A vrna_fold_compound_t containing the energy parameters and model details
 *  @param[in,out]  pt      The pair table containing the start structure. Used to update to the final structure after execution of this function
 *  @param[in]      options Options to modify the behavior of this function
 *  @return                 A list of transition moves (default), or NULL (if options & #VRNA_PATH_NO_TRANSITION_OUTPUT)
 */
vrna_move_t *
vrna_path(vrna_fold_compound_t  *vc,
          short                 *pt,
          unsigned int          steps,
          unsigned int          options);


/**
 *  @brief Compute a steepest descent / gradient path, store the final structure, and return a
 *  list of transition moves from the start to the final structure.
 *
 *  This function computes, given a start structure in pair table format, a steepest descent path,
 *  updates the pair table to the final structure of the path. Finally, if not requested otherwise
 *  by using the #VRNA_PATH_NO_TRANSITION_OUTPUT flag in the @p options field, this function returns a list
 *  of individual transitions that lead from the start to the final
 *  structure if requested.
 *
 *  @note   Since the result is written to the input structure you may want to use
 *          vrna_ptable_copy() before calling this function to keep the initial structure
 *
 *  @see    vrna_path_random(), vrna_path(), vrna_ptable(), vrna_ptable_copy(), vrna_fold_compound()
 *          #VRNA_MOVESET_DEFAULT, #VRNA_MOVESET_SHIFT, #VRNA_PATH_NO_TRANSITION_OUTPUT
 *
 *  @param[in]      vc      A vrna_fold_compound_t containing the energy parameters and model details
 *  @param[in,out]  pt      The pair table containing the start structure. Used to update to the final structure after execution of this function
 *  @param[in]      options Options to modify the behavior of this function
 *  @return                 A list of transition moves (default), or NULL (if options & #VRNA_PATH_NO_TRANSITION_OUTPUT)
 */
vrna_move_t *
vrna_path_gradient(vrna_fold_compound_t *vc,
                   short                *pt,
                   unsigned int         options);


/**
 *  @brief Generate a random walk / path of a given length, store the final structure, and return a
 *  list of transition moves from the start to the final structure.
 *
 *  This function generates, given a start structure in pair table format,  a random walk / path,
 *  updates the pair table to the final structure of the path. Finally, if not requested otherwise
 *  by using the #VRNA_PATH_NO_TRANSITION_OUTPUT flag in the @p options field, this function returns a list
 *  of individual transitions that lead from the start to the final
 *  structure if requested.
 *
 *  @note   Since the result is written to the input structure you may want to use
 *          vrna_ptable_copy() before calling this function to keep the initial structure
 *
 *  @see    vrna_path_gradient(), vrna_path(), vrna_ptable(), vrna_ptable_copy(), vrna_fold_compound()
 *          #VRNA_MOVESET_DEFAULT, #VRNA_MOVESET_SHIFT, #VRNA_PATH_NO_TRANSITION_OUTPUT
 *
 *  @param[in]      vc      A vrna_fold_compound_t containing the energy parameters and model details
 *  @param[in,out]  pt      The pair table containing the start structure. Used to update to the final structure after execution of this function
 *  @param[in]      steps   The length of the path, i.e. the total number of transitions / moves
 *  @param[in]      options Options to modify the behavior of this function
 *  @return                 A list of transition moves (default), or NULL (if options & #VRNA_PATH_NO_TRANSITION_OUTPUT)
 */
vrna_move_t *
vrna_path_random(vrna_fold_compound_t *vc,
                 short                *pt,
                 unsigned int         steps,
                 unsigned int         options);


/**
 *  @}
 */

#endif /* VIENNA_RNA_PACKAGE_WALK_H */
