#ifndef VIENNA_RNA_PACKAGE_NEIGHBOR_H
#define VIENNA_RNA_PACKAGE_NEIGHBOR_H

/**
 *  @file     ViennaRNA/landscape/neighbor.h
 *  @ingroup  neighbors
 *  @brief    Methods to compute the neighbors of an RNA secondary structure.
 */

/**
 *  @addtogroup neighbors
 *  @{
 *
 *  @brief Different functions to generate structural neighbors of a secondary structure according to a particular Move Set.
 *
 *  This module contains methods to compute the neighbors of an RNA secondary structure. Neighbors of a given structure
 *  are all structures that differ in exactly one base pair. That means one can insert an delete base pairs in the
 *  given structure. These insertions and deletions of base pairs are usually called moves. A third move which is
 *  considered in these methods is a shift move. A shifted base pair has one stable position and one position that
 *  changes. These moves are encoded as follows: \n
 *  - insertion: (i, j) where i,j > 0 \n
 *  - deletion: (i, j) where i,j < 0 \n
 *  shift: (i, j) where either i > 0, j < 0 or i < 0, j > 0 \n
 *
 *  The negative position of a shift indicates the position that has changed.
 *
 *  \code
 *  Example:
 *           We have given a sequence and a structure.
 *           Sequence  AAGGAAACC
 *           Structure ..(.....)
 *           Indices   123456789
 *
 *           The given base pair is (3,9) and the neighbors are the insertion (4, 8), the deletion (-3,-9), the shift (3,-8)
 *           and the shift (-4, 9).
 *           This leads to the neighbored structures:
 *           ...(....)
 *           .........
 *           ...(...).
 *           ....(...)
 *  \endcode
 *
 *  A simple method to construct all insertions is to iterate over the positions of a sequence twice. The first
 *  iteration has the index i in [1, sequence length], the second iteration has the index j in [i+1, sequence length].
 *  All pairs (i,j) with compatible letters and which are non-crossing with present base pairs are valid neighbored
 *  insertion moves.
 *  Valid deletion moves are all present base pairs with negative sign.
 *  Valid shift moves are constructed by taking all paired positions as fix position of a shift move and iterating over
 *  all positions of the sequence. If the letters of a position are compatible and if it the move is non-crossing with
 *  existing base pairs, we have a valid shift move.
 *  The method of generating shift moves can be accelerated by skipping neighbored base pairs.
 *
 *  If we need to construct all neighbors several times for subsequent moves, we can speed up the task by using
 *  the move set of the previous structure. The previous move set has to be filtered, such that all moves that would
 *  cross the next selected move are non-crossing. Next, the selected move has to be removed. Then one has to only
 *  to generate all moves that were not possible before.
 *  One move is the inverted selected move (if it was an insertion, simply make the indices negative).
 *  The generation of all other new moves is different and depends on the selected move. It is easy for an insertion move,
 *  because we have only to include all non-crossing shift moves, that are possible with the new base pair. For that we can
 *  either iterate over the sequence or we can select all crossing shift moves in the filter procedure and convert them into
 *  shifts.
 *
 *  The generation of new moves given a deletion is a little bit more complex, because we can create more moves. At first
 *  we can insert the deleted pair as insertion move. Then we generate all insertions that would have crossed the deleted
 *  base pair. Finally we construct all crossing shift moves.
 *
 *  If the given move is a shift, we can save much time by specifying the intervals for the generation of new moves.
 *  The interval which was enclosed by the positive position of the shift move and the previous paired position is the
 *  freed interval after applying the move. This freed interval includes all positions and base pairs that we need to
 *  construct new insertions and shifts. All these new moves have one position in the freed interval and the other position
 *  in the environment of the freed interval. The environment are all position which are outside the freed interval, but within
 *  the same enclosing loop of the shift move. The environment for valid base pairs can be divided into one or more intervals,
 *  depending on the shift move. The following examples describe a few scenarios to specify the intervals of the environment.
 *
 * \internal
 * Example 1:
 * increase the freed interval with a shift move
 * AAAAGACAAGAAACAAAAGAGAAACAACAAACAAGAAACAAACAAAA
 * ....(....(...)....(.(...)..)......(...)...).... // structure before the shift
 * ....(....(...)....(.(...)......)..(...)...).... // structure after the shift
 * ............................[__]............... // freed interval
 * ..................[________]................... // interval that can pair with the freed interval
 *
 * Example 2:
 * switch the freed interval with a shift move
 * AAAAGACAAGAAACAAAAGAGAAACAACAAACAAGAAACAAACAAAA
 * ....(....(...)....(.(...)..)......(...)...).... // structure before the shift
 * ....(.(..(...)....).(...).........(...)...).... // structure after the shift
 * ...................[_______]................... // freed interval
 * ....[_].....................[_____________].... // intervals that can pair with the freed interval
 *
 * Example 3:
 * decrease the freed interval with a shift move
 * AAAAGACAAGAAACAAAAGAGAAACAACAAACAAGAAACAAACAAAA
 * ....(....(...)....(.(...)......)..(...)...).... // structure before the shift
 * ....(....(...)....(.(...)..)......(...)...).... // structure after the shift
 * ............................[__]............... // freed interval
 * ....[____________].............[__________].... // intervals that can pair with the freed interval
 * \endinternal
 *
 *  @image html shift_move_intervals.svg
 *  @image latex shift_move_intervals.eps
 *
 *  Given the intervals of the environment and the freed interval, the new shift moves can be constructed quickly.
 *  One has to take all positions of pairs from the environment in order to create valid pairs with positions in the freed interval.
 *  The same procedure can be applied for the other direction. This is taking all paired positions within the freed interval
 *  in order to look for pairs with valid positions in the intervals of the environment.
 *
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
typedef void (vrna_callback_move_update)(vrna_fold_compound_t *fc,
                                         vrna_move_t          neighbor,
                                         unsigned int         state,
                                         void                 *data);


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
 * #VRNA_MOVESET_INSERTION, #VRNA_MOVESET_DELETION, #VRNA_MOVESET_SHIFT, #VRNA_MOVESET_DEFAULT
 *
 * @param[in] vc        A vrna_fold_compound_t containing the energy parameters and model details
 * @param[in] pt        The pair table representation of the structure
 * @param     options   Options to modify the behavior of this function, e.g. available move set
 * @return              Neighbors as a list of moves / transitions (the last element in the list has both of its fields set to 0)
 */
vrna_move_t *
vrna_neighbors(vrna_fold_compound_t *vc,
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
 * #VRNA_MOVESET_INSERTION, #VRNA_MOVESET_DELETION, #VRNA_MOVESET_SHIFT, #VRNA_MOVESET_DEFAULT
 *
 * @param[in]   vc                  A vrna_fold_compound_t containing the energy parameters and model details
 * @param[in]   curr_move           The move that was/will be applied to @p prev_pt
 * @param[in]   prev_pt             A pair table representation of the structure before @p curr_move is/was applied
 * @param[in]   prev_neighbors      The list of neighbors of @p prev_pt
 * @param       size_prev_neighbors The size of @p prev_neighbors, i.e. the lists length
 * @param[out]  size_neighbors      A pointer to store the size / length of the new neighbor list
 * @param       options             Options to modify the behavior of this function, e.g. available move set
 * @return                          Neighbors as a list of moves / transitions (the last element in the list has both of its fields set to 0)
 */
vrna_move_t *
vrna_neighbors_successive(const vrna_fold_compound_t  *vc,
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
 *        #vrna_callback_move_update
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
vrna_move_neighbor_diff_cb(vrna_fold_compound_t       *fc,
                           short                      *ptable,
                           vrna_move_t                move,
                           vrna_callback_move_update  *cb,
                           void                       *data,
                           unsigned int               options);


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
