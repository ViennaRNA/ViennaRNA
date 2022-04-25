#ifndef VIENNA_RNA_PACKAGE_MOVE_H
#define VIENNA_RNA_PACKAGE_MOVE_H


/**
 *  @file     ViennaRNA/landscape/move.h
 *  @ingroup  neighbors
 *  @brief    Methods to operate with structural neighbors of RNA secondary structures.
 */

/**
 *  @addtogroup neighbors
 *  @{
 */


/**
 *  @brief  A single move that transforms a secondary structure into one of its neighbors
 */
typedef struct vrna_move_s vrna_move_t;

/**
 *  @brief Option flag indicating insertion move
 *
 *  @see vrna_neighbors(), vrna_neighbors_successive, vrna_path()
 */
#define VRNA_MOVESET_INSERTION   4

/**
 *  @brief Option flag indicating deletion move
 *
 *  @see vrna_neighbors(), vrna_neighbors_successive, vrna_path()
 */
#define VRNA_MOVESET_DELETION    8

/**
 *  @brief Option flag indicating shift move
 *
 *  @see vrna_neighbors(), vrna_neighbors_successive, vrna_path()
 */
#define VRNA_MOVESET_SHIFT       16
/**
 *  @brief Option flag indicating moves without lonely base pairs
 *
 *  @see vrna_neighbors(), vrna_neighbors_successive, vrna_path()
 */
#define VRNA_MOVESET_NO_LP       32

/**
 *  @brief Option flag indicating default move set, i.e. insertions/deletion of a base pair
 *
 *  @see vrna_neighbors(), vrna_neighbors_successive, vrna_path()
 */
#define VRNA_MOVESET_DEFAULT     (VRNA_MOVESET_INSERTION | VRNA_MOVESET_DELETION)

#define VRNA_MOVE_NO_APPLY       64

/**
 * @brief   An atomic representation of the transition / move from one structure to its neighbor
 *
 *  An atomic transition / move may be one of the following:
 *  - a <strong>base pair insertion</strong>,
 *  - a <strong>base pair removal</strong>, or
 *  - a <strong>base pair shift</strong> where an existing base pair changes one of its
 *    pairing partner.
 *
 *  These moves are encoded by two integer values that represent the affected 5' and 3'
 *  nucleotide positions. Furthermore, we use the following convention on the signedness
 *  of these encodings:
 *  - both values are positive for <em>insertion moves</em>
 *  - both values are negative for <em>base pair removals</em>
 *  - both values have different signedness for <em>shift moves</em>, where the positive
 *    value indicates the nucleotide that stays constant, and the others absolute value
 *    is the new pairing partner
 *
 *  @note   A value of 0 in either field is used as list-end indicator and doesn't represent
 *          any valid move.
 */
struct vrna_move_s {
  int         pos_5;  /**< @brief The (absolute value of the) 5' position of a base pair, or any position of a shifted pair */
  int         pos_3;  /**< @brief The (absolute value of the) 3' position of a base pair, or any position of a shifted pair */
  vrna_move_t *next;  /**< @brief The next base pair (if an elementary move changes more than one base pair), or @p NULL
                       *   Has to be terminated with move 0,0
                       */
};


/**
 *  @brief  Create an atomic move
 *
 *  @see #vrna_move_s
 *
 *  @param  pos_5   The 5' position of the move (positive for insertions, negative for removal, any value for shift moves)
 *  @param  pos_3   The 3' position of the move (positive for insertions, negative for removal, any value for shift moves)
 *  @return         An atomic move as specified by @p pos_5 and @p pos_3
 */
vrna_move_t
vrna_move_init(int  pos_5,
               int  pos_3);


/**
 * delete all moves in a zero terminated list.
 */
void
vrna_move_list_free(vrna_move_t *moves);


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
vrna_move_apply_db(char               *structure,
                   const short        *pt,
                   const vrna_move_t  *m);


/**
 *  @brief  Test whether a move is a base pair removal
 *
 *  @param  m   The move to test against
 *  @return     Non-zero if the move is a base pair removal, 0 otherwise
 */
int
vrna_move_is_removal(const vrna_move_t *m);


/**
 *  @brief  Test whether a move is a base pair insertion
 *
 *  @param  m   The move to test against
 *  @return     Non-zero if the move is a base pair insertion, 0 otherwise
 */
int
vrna_move_is_insertion(const vrna_move_t *m);


/**
 *  @brief  Test whether a move is a base pair shift
 *
 *  @param  m   The move to test against
 *  @return     Non-zero if the move is a base pair shift, 0 otherwise
 */
int
vrna_move_is_shift(const vrna_move_t *m);


/**
 *  @brief  Compare two moves
 *
 *  The function compares two moves @p m and @p b and returns
 *  whether move @p m is lexicographically smaller (-1), larger (1)
 *  or equal to move @p b.
 *
 *  If any of the moves @p m or @p b is a shift move, this
 *  comparison only makes sense in a structure context. Thus,
 *  the third argument with the current structure must be provided.
 *
 *  @note    This function returns 0 (equality) upon any error, e.g. missing input
 *
 *  @warning Currently, shift moves are not supported!
 *
 *  @param    m   The first move of the comparison
 *  @param    b   The second move of the comparison
 *  @param    pt  The pair table of the current structure that is compatible with both moves (maybe NULL if moves are guaranteed to be no shifts)
 *  @return       -1 if @p m < @p b, 1 if @p m > @p b, 0 otherwise
 */
int
vrna_move_compare(const vrna_move_t *m,
                  const vrna_move_t *b,
                  const short       *pt);


/**
 *  @}
 */

#endif
