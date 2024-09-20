#ifndef VIENNA_RNA_PACKAGE_COMBINATORICS_H
#define VIENNA_RNA_PACKAGE_COMBINATORICS_H

/**
 *  @file     combinatorics/combinatorics.h
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
 *  as published by @rstinline :cite:t:`sawada:2003` @endrst.
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
 *  @brief  Obtain a list of k-combinations with repetition (n multichoose k)
 *
 *  This function compiles a list of k-combinations, or k-multicombination, i.e.
 *  a list of multisubsets of size k from a set of integer values from 0 to n - 1.
 *  For that purpose, we enumerate n + k - 1 choose k and decrease each index
 *  position i by i to obtain n multichoose k.
 *
 *  @param  n   Maximum number to choose from (interval of integers from 0 to @p n - 1)
 *  @param  k   Number of elements to choose, i.e. size of each multisubset
 *  @return     A list of lists of elements of combinations (last entry is terminated by @b NULL
 */
unsigned int **
vrna_n_multichoose_k(size_t n,
                     size_t k);


/**
 *  @brief Generate a sequence of Boustrophedon distributed numbers
 *
 *  This function generates a sequence of positive natural numbers within the
 *  interval @f$ [start, end] @f$ in a Boustrophedon fashion. That is, the
 *  numbers @f$ start, \ldots, end @f$ in the resulting list are alternating
 *  between left and right ends of the interval while progressing to the inside,
 *  i.e. the list consists of a sequence of natural numbers of the form:
 *
 *  @f[ start, end, start + 1, end - 1, start + 2, end - 2, \ldots @f]
 *
 *  The resulting list is 1-based and contains the length of the sequence
 *  of numbers at it's 0-th position.
 *
 *  Upon failure, the function returns @b NULL
 *
 *  @see vrna_boustrophedon_pos()
 *
 *  @param  start   The first number of the list (left side of the interval)
 *  @param  end     The last number of the list (right side of the interval)
 *  @return         A list of alternating numbers from the interval @f$ [start, end] @f$ (or @b NULL on error)
 */
unsigned int *
vrna_boustrophedon(size_t start,
                   size_t end);


/**
 *  @brief Obtain the i-th element in a Boustrophedon distributed interval of natural numbers
 *
 *
 *  @see vrna_boustrophedon()
 *
 *  @param  start   The first number of the list (left side of the interval)
 *  @param  end     The last number of the list (right side of the interval)
 *  @param  pos     The index of the number within the Boustrophedon distributed sequence (1-based)
 *  @return         The @p pos-th element in the Boustrophedon distributed sequence of natural numbers of the interval
 */
unsigned int
vrna_boustrophedon_pos(size_t start,
                       size_t end,
                       size_t pos);


/**
 *  @}
 */
#endif
