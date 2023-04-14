#ifndef VIENNA_RNA_PACKAGE_NEIGHBOR_H
#define VIENNA_RNA_PACKAGE_NEIGHBOR_H

#ifdef VRNA_WARN_DEPRECATED
# if defined(DEPRECATED)
#   undef DEPRECATED
# endif
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
 *  @file     ViennaRNA/landscape/neighbor.h
 *  @ingroup  neighbors
 *  @brief    Methods to compute the neighbors of an RNA secondary structure.
 */

/**
 *  @addtogroup neighbors
 *  @{
 */

#include <ViennaRNA/fold_compound.h>
#include <ViennaRNA/landscape/move.h>

/**
 *  @brief  Prototype of the neighborhood update callback
 *
 *  @see  vrna_move_neighbor_diff_cb(), #VRNA_NEIGHBOR_CHANGE, #VRNA_NEIGHBOR_INVALID, #VRNA_NEIGHBOR_NEW
 *
 *  @param  fc        The fold compound the calling function is working on
 *  @param  neighbor  The move that generates the (changed or new) neighbor
 *  @param  state     The state of the neighbor (move) as supplied by argument @p neighbor
 *  @param  data      Some arbitrary data pointer as passed to vrna_move_neighbor_diff_cb()
 */
typedef void (*vrna_move_update_f)(vrna_fold_compound_t *fc,
                                   vrna_move_t          neighbor,
                                   unsigned int         state,
                                   void                 *data);

DEPRECATED(typedef void (vrna_callback_move_update)(vrna_fold_compound_t *fc,
                                                    vrna_move_t neighbor,
                                                    unsigned int state,
                                                    void *data),
           "Use vrna_move_update_f instead!");


/**
 *  @brief  State indicator for a neighbor that has been changed
 *
 *  @see vrna_move_neighbor_diff_cb()
 */
#define VRNA_NEIGHBOR_CHANGE   1


/**
 *  @brief  State indicator for a neighbor that has been invalidated
 *
 *  @see vrna_move_neighbor_diff_cb()
 */
#define VRNA_NEIGHBOR_INVALID   2


/**
 *  @brief  State indicator for a neighbor that has become newly available
 *
 *  @see vrna_move_neighbor_diff_cb()
 */
#define VRNA_NEIGHBOR_NEW       3


/**
 * @brief Alters the loopIndices array that was constructed with vrna_loopidx_from_ptable().
 *
 * The loopIndex of the current move will be inserted.
 * The correctness of the input will not be checked because the speed should be optimized.
 *
 * @param[in,out] loopidx   The loop index data structure that needs an update
 * @param[in]     pt        A pair table on which the move will be executed
 * @param         length    The length of the structure
 * @param[in]     m         The move that is applied to the current structure
 */
void
vrna_loopidx_update(int               *loopidx,
                    const short       *pt,
                    int               length,
                    const vrna_move_t *m);


/**
 * @brief Generate neighbors of a secondary structure
 *
 * This function allows one to generate all structural neighbors (according to a particular move set)
 * of an RNA secondary structure. The neighborhood is then returned as a list of transitions / moves
 * required to transform the current structure into the actual neighbor.
 *
 * @see vrna_neighbors_successive(), vrna_move_apply(),
 *      #VRNA_MOVESET_INSERTION, #VRNA_MOVESET_DELETION, #VRNA_MOVESET_SHIFT, #VRNA_MOVESET_DEFAULT
 *
 * @param[in] fc        A vrna_fold_compound_t containing the energy parameters and model details
 * @param[in] pt        The pair table representation of the structure
 * @param     options   Options to modify the behavior of this function, e.g. available move set
 * @return              Neighbors as a list of moves / transitions (the last element in the list has both of its fields set to 0)
 */
vrna_move_t *
vrna_neighbors(vrna_fold_compound_t *fc,
               const short          *pt,
               unsigned int         options);


/**
 * @brief Generate neighbors of a secondary structure (the fast way)
 *
 * This function implements a fast way to generate all neighbors of a secondary structure
 * that results from successive applications of individual moves. The speed-up results from
 * updating an already known list of valid neighbors before the individual move towards the
 * current structure took place. In essence, this function removes neighbors that are not
 * accessible anymore and inserts neighbors emerging after a move took place.
 *
 * @see vrna_neighbors(), vrna_move_apply(),
 *      #VRNA_MOVESET_INSERTION, #VRNA_MOVESET_DELETION, #VRNA_MOVESET_SHIFT, #VRNA_MOVESET_DEFAULT
 *
 * @param[in]   fc                  A vrna_fold_compound_t containing the energy parameters and model details
 * @param[in]   curr_move           The move that was/will be applied to @p prev_pt
 * @param[in]   prev_pt             A pair table representation of the structure before @p curr_move is/was applied
 * @param[in]   prev_neighbors      The list of neighbors of @p prev_pt
 * @param       size_prev_neighbors The size of @p prev_neighbors, i.e. the lists length
 * @param[out]  size_neighbors      A pointer to store the size / length of the new neighbor list
 * @param       options             Options to modify the behavior of this function, e.g. available move set
 * @return                          Neighbors as a list of moves / transitions (the last element in the list has both of its fields set to 0)
 */
vrna_move_t *
vrna_neighbors_successive(const vrna_fold_compound_t  *fc,
                          const vrna_move_t           *curr_move,
                          const short                 *prev_pt,
                          const vrna_move_t           *prev_neighbors,
                          int                         size_prev_neighbors,
                          int                         *size_neighbors,
                          unsigned int                options);


/**
 *  @brief  Apply a move to a secondary structure and indicate which neighbors have changed consequentially
 *
 *  This function applies a move to a secondary structure and explores the local neighborhood of the
 *  affected loop. Any changes to previously compatible neighbors that have been affected by this loop
 *  will be reported through a callback function. In particular, any of the three cases might appear:
 *  - A previously available neighbor move has changed, usually the free energy change of the move (#VRNA_NEIGHBOR_CHANGE)
 *  - A previously available neighbor move became invalid (#VRNA_NEIGHBOR_INVALID)
 *  - A new neighbor move becomes available (#VRNA_NEIGHBOR_NEW)
 *
 *  @see  vrna_move_neighbor_diff(), #VRNA_NEIGHBOR_CHANGE, #VRNA_NEIGHBOR_INVALID, #VRNA_NEIGHBOR_NEW,
 *        #vrna_move_update_f, #VRNA_MOVE_NO_APPLY
 *
 *  @param  fc        A fold compound for the RNA sequence(s) that this function operates on
 *  @param  ptable    The current structure as pair table
 *  @param  move      The move to apply
 *  @param  cb        The address of the callback function that is passed the neighborhood changes
 *  @param  data      An arbitrary data pointer that will be passed through to the callback function @p cb
 *  @param  options   Options to modify the behavior of this function, .e.g available move set
 *  @return           Non-zero on success, 0 otherwise
 */
int
vrna_move_neighbor_diff_cb(vrna_fold_compound_t *fc,
                           short                *ptable,
                           vrna_move_t          move,
                           vrna_move_update_f   cb,
                           void                 *data,
                           unsigned int         options);


/**
 *  @brief  Apply a move to a secondary structure and indicate which neighbors have changed consequentially
 *
 *  Similar to vrna_move_neighbor_diff_cb(), this function applies a move to a secondary structure and
 *  reports back the neighbors of the current structure become affected by this move. Instead of executing
 *  a callback for each of the affected neighbors, this function compiles two lists of neighbor moves, one
 *  that is returned and consists of all moves that are novel or may have changed in energy, and a second,
 *  @p invalid_moves, that consists of all the neighbor moves that become invalid, respectively.
 *
 *  @param  fc        A fold compound for the RNA sequence(s) that this function operates on
 *  @param  ptable    The current structure as pair table
 *  @param  move      The move to apply
 *  @param  invalid_moves The address of a move list where the function stores those moves that become invalid
 *  @param  options   Options to modify the behavior of this function, .e.g available move set
 *  @return           A list of moves that might have changed in energy or are novel compared to the structure before application of the move
 */
vrna_move_t *
vrna_move_neighbor_diff(vrna_fold_compound_t  *fc,
                        short                 *ptable,
                        vrna_move_t           move,
                        vrna_move_t           **invalid_moves,
                        unsigned int          options);


/**
 *  @}
 */
#endif /* VIENNA_RNA_PACKAGE_NEIGHBOR_H */
