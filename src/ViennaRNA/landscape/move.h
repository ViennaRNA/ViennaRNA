#ifndef VIENNA_RNA_PACKAGE_MOVE_H
#define VIENNA_RNA_PACKAGE_MOVE_H


/**
 *  @file     ViennaRNA/landscape/move.h
 *  @ingroup  neighbors
 *  @brief    Methods to operate with structural neighbors of RNA secondary structures.
 */

typedef struct vrna_move_s vrna_move_t;

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
 *  @brief Option flag indicating moves without lonely base pairs
 *  @see vrna_neighbors(), vrna_neighbors_successive, vrna_path()
 */
#define VRNA_MOVESET_NO_LP       32

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
  int         pos_5;  /**< The 5' position of a base pair, or any position of a shifted pair */
  int         pos_3;  /**< The 3' position of a base pair, or any position of a shifted pair */
  vrna_move_t *next;  /**< The next base pair (if an elementary move changes more than one base pair)
                       *   Has to be terminated with move 0,0 */
};

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


#endif
