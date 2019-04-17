#ifndef VIENNA_RNA_PACKAGE_COMBINATORICS_H
#define VIENNA_RNA_PACKAGE_COMBINATORICS_H

/**
 *  @file     combinatorics.h
 *  @ingroup  utils, combinatorics_utils
 *  @brief    Various implementations that deal with combinatorial aspects of objects
 */

#include <ViennaRNA/fold_compound.h>

/**
 *  @addtogroup   combinatorics_utils
 *  @{
 *  @brief  Implementations to solve various combinatorial aspects for strings of objects
 */

/**
 *  @brief  Enumerate all necklaces with fixed content
 *
 *  This function implements <em>A fast algorithm to generate necklaces with fixed content</em>
 *  as published by Joe Sawada in 2003 @cite sawada:2003.
 *
 *  The function receives a list of counts (the elements on the necklace) for each
 *  type of object within a necklace. The list starts at index 0 and ends with an
 *  entry that has a count of 0. The algorithm then enumerates all non-cyclic permutations
 *  of the content, returned as a list of necklaces. This list, again, is zero-terminated,
 *  i.e. the last entry of the list is a @p NULL pointer.
 *
 *  @param  type_counts A 0-terminated list of entity counts
 *  @return             A list of all non-cyclic permutations of the entities
 */
unsigned int **
vrna_enumerate_necklaces(const unsigned int *type_counts);


/**
 *  @brief Determine the order of rotational symmetry for a string of objects represented
 *         by natural numbers
 *
 *  The algorithm applies a fast search of the provided string within itself, assuming
 *  the end of the string wraps around to connect with it's start.
 *  For example, a string of the form @p 011011 has rotational symmetry
 *  of order @p 2
 *
 *  This is a simplified version of vrna_rotational_symmetry_pos_num() that may
 *  be useful if one is only interested in the degree of rotational symmetry but
 *  not the actual set of rotational symmetric strings.
 *
 *  @see vrna_rotational_symmetry_pos_num(), vrna_rotationa_symmetry()
 *
 *  @param  string        The string of elements encoded as natural numbers
 *  @param  string_length The length of the string
 *  @return               The order of rotational symmetry
 */
unsigned int
vrna_rotational_symmetry_num(const unsigned int *string,
                             size_t             string_length);


/**
 *  @brief Determine the order of rotational symmetry for a string of objects represented
 *         by natural numbers
 *
 *  The algorithm applies a fast search of the provided string within itself, assuming
 *  the end of the string wraps around to connect with it's start.
 *  For example, a string of the form @p 011011 has rotational symmetry
 *  of order @p 2
 *
 *  If the argument @p positions is not @p NULL, the function stores an array of
 *  string start positions for rotational shifts that map the string back onto itself.
 *  This array has length of order of rotational symmetry, i.e. the number returned
 *  by this function. The first element @p positions[0] always contains a shift value
 *  of @p 0 representing the trivial rotation.
 *
 *  @note Do not forget to release the memory occupied by @p positions after a successful
 *        execution of this function.
 *
 *  @see vrna_rotational_symmetry_num(), vrna_rotational_symmetry(), vrna_rotational_symmetry_pos()
 *
 *  @param  string        The string of elements encoded as natural numbers
 *  @param  string_length The length of the string
 *  @param  positions     A pointer to an (undefined) list of alternative string start
 *                        positions that lead to an identity mapping (may be NULL)
 *  @return               The order of rotational symmetry
 */
unsigned int
vrna_rotational_symmetry_pos_num(const unsigned int *string,
                                 size_t             string_length,
                                 unsigned int       **positions);


/**
 *  @brief Determine the order of rotational symmetry for a NULL-terminated string of ASCII characters
 *
 *  The algorithm applies a fast search of the provided string within itself, assuming
 *  the end of the string wraps around to connect with it's start.
 *  For example, a string of the form @p AABAAB has rotational symmetry
 *  of order @p 2
 *
 *  This is a simplified version of vrna_rotational_symmetry_pos() that may
 *  be useful if one is only interested in the degree of rotational symmetry but
 *  not the actual set of rotational symmetric strings.
 *
 *  @see vrna_rotational_symmetry_pos(), vrna_rotationa_symmetry_num()
 *
 *  @param  string        A NULL-terminated string of characters
 *  @return               The order of rotational symmetry
 */
unsigned int
vrna_rotational_symmetry(const char *string);


/**
 *  @brief Determine the order of rotational symmetry for a NULL-terminated string of ASCII characters
 *
 *  The algorithm applies a fast search of the provided string within itself, assuming
 *  the end of the string wraps around to connect with it's start.
 *  For example, a string of the form @p AABAAB has rotational symmetry
 *  of order @p 2
 *
 *  If the argument @p positions is not @p NULL, the function stores an array of
 *  string start positions for rotational shifts that map the string back onto itself.
 *  This array has length of order of rotational symmetry, i.e. the number returned
 *  by this function. The first element @p positions[0] always contains a shift value
 *  of @p 0 representing the trivial rotation.
 *
 *  @note Do not forget to release the memory occupied by @p positions after a successful
 *        execution of this function.
 *
 *  @see vrna_rotational_symmetry(), vrna_rotational_symmetry_num(), vrna_rotational_symmetry_num_pos()
 *
 *  @param  string        A NULL-terminated string of characters
 *  @param  positions     A pointer to an (undefined) list of alternative string start
 *                        positions that lead to an identity mapping (may be NULL)
 *  @return               The order of rotational symmetry
 */
unsigned int
vrna_rotational_symmetry_pos(const char   *string,
                             unsigned int **positions);


/**
 *  @brief  Determine the order of rotational symmetry for a dot-bracket structure
 *
 *  Given a (permutation of multiple) RNA strand(s) and a particular secondary structure
 *  in dot-bracket notation, compute the degree of rotational symmetry. In case there
 *  is only a single linear RNA strand, the structure always has degree 1, as there are
 *  no rotational symmetries due to the direction of the nucleic acid sequence and
 *  the fixed positions of 5' and 3' ends. However, for circular RNAs, rotational
 *  symmetries might arise if the sequence consists of a concatenation of @f$k@f$
 *  identical subsequences.
 *
 *  This is a simplified version of vrna_rotational_symmetry_db_pos() that may
 *  be useful if one is only interested in the degree of rotational symmetry but
 *  not the actual set of rotational symmetric strings.
 *
 *  @see vrna_rotational_symmetry_db_pos(), vrna_rotational_symmetry(), vrna_rotational_symmetry_num()
 *
 *  @param  fc          A fold_compound data structure containing the nucleic acid sequence(s),
 *                      their order, and model settings
 *  @param  structure   The dot-bracket structure the degree of rotational symmetry is checked for
 *  @return             The degree of rotational symmetry of the @p structure (0 in case of any errors)
 */
unsigned int
vrna_rotational_symmetry_db(vrna_fold_compound_t  *fc,
                            const char            *structure);


/**
 *  @brief  Determine the order of rotational symmetry for a dot-bracket structure
 *
 *  Given a (permutation of multiple) RNA strand(s) and a particular secondary structure
 *  in dot-bracket notation, compute the degree of rotational symmetry. In case there
 *  is only a single linear RNA strand, the structure always has degree 1, as there are
 *  no rotational symmetries due to the direction of the nucleic acid sequence and
 *  the fixed positions of 5' and 3' ends. However, for circular RNAs, rotational
 *  symmetries might arise if the sequence consists of a concatenation of @f$k@f$
 *  identical subsequences.
 *
 *  If the argument @p positions is not @p NULL, the function stores an array of
 *  string start positions for rotational shifts that map the string back onto itself.
 *  This array has length of order of rotational symmetry, i.e. the number returned
 *  by this function. The first element @p positions[0] always contains a shift value
 *  of @p 0 representing the trivial rotation.
 *
 *  @note Do not forget to release the memory occupied by @p positions after a successful
 *        execution of this function.
 *
 *  @see  vrna_rotational_symmetry_db(), vrna_rotational_symmetry_pos(),
 *        vrna_rotational_symmetry_pos_num()
 *
 *  @param  fc          A fold_compound data structure containing the nucleic acid sequence(s),
 *                      their order, and model settings
 *  @param  structure   The dot-bracket structure the degree of rotational symmetry is checked for
 *  @param  positions   A pointer to an (undefined) list of alternative string start
 *                      positions that lead to an identity mapping (may be NULL)
 *  @return             The degree of rotational symmetry of the @p structure (0 in case of any errors)
 */
unsigned int
vrna_rotational_symmetry_db_pos(vrna_fold_compound_t  *fc,
                                const char            *structure,
                                unsigned int          **positions);


/**
 *  @}
 */
#endif
