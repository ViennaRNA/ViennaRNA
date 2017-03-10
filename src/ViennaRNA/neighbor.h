#ifndef VIENNA_RNA_PACKAGE_NEIGHBOR_H
#define VIENNA_RNA_PACKAGE_NEIGHBOR_H

/**
 *  @file     neighbor.h
 *  @ingroup  neighbors
 *  @brief    Methods to compute the neighbors of an RNA secondary structure.
 */

/**
 *  @addtogroup neighbors
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
 *  @{
 *  @ingroup  neighbors
 */

typedef struct vrna_move_s vrna_move_t;

#include <ViennaRNA/data_structures.h>

/**
 *  @brief Option flag indicating insertion move
 *  @see vrna_neighbors(), vrna_neighbors_successive, vrna_path()
 */
#define VRNA_MOVESET_INSERTION   4

/**
 *  @brief Option flag indicating deletion move
 *  @see vrna_neighbors(), vrna_neighbors_successive, vrna_path()
 */
#define VRNA_MOVESET_DELETION    8

/**
 *  @brief Option flag indicating shift move
 *  @see vrna_neighbors(), vrna_neighbors_successive, vrna_path()
 */
#define VRNA_MOVESET_SHIFT       16

/**
 *  @brief Option flag indicating default move set, i.e. insertions/deletion of a base pair
 *  @see vrna_neighbors(), vrna_neighbors_successive, vrna_path()
 */
#define VRNA_MOVESET_DEFAULT     (VRNA_MOVESET_INSERTION | VRNA_MOVESET_DELETION)


/**
 * @brief   An atomic representation of the transition / move from one structure to its neighbor
 *
 * An atomic transition / move may be (a) the insertion of a base pair (both fields are positive),
 * (b) the deletion of a base pair (both fields are negative), or (c) a base pair shift where one
 * position stays constant while the other is allowed to shift along the same loop it resides in
 * (one field position and the other negative, where the positive field indicates the constant
 * position and the absolute value of the negative field is the new position of the pairing partner).
 *
 * A value of 0 is either field is typically used to indicate the lists last element.
 */
struct vrna_move_s {
  int pos_5;  /**< The 5' position of a base pair, or any position of a shifted pair */
  int pos_3;  /**< The 3' position of a base pair, or any position of a shifted pair */
};

/**
 * @brief Apply a particular move / transition to a secondary structure, i.e. transform a structure
 *
 * @param[in,out] pt The pair table representation of the secondary structure
 * @param[in]     m  The move to apply
 */
void
vrna_move_apply(short             *pt,
                const vrna_move_t *m);


void
vrna_move_apply_to_db(char              *structure,
                      const short       *pt,
                      const vrna_move_t *m);


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
 *  @}
 */
#endif /* VIENNA_RNA_PACKAGE_NEIGHBOR_H */
