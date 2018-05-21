#ifndef VIENNA_RNA_PACKAGE_COMBINATORICS_H
#define VIENNA_RNA_PACKAGE_COMBINATORICS_H

/**
 *  @file     combinatorics.h
 *  @ingroup  utils
 *  @brief     Various implementations that deal with combinatorial aspects of RNA/DNA folding
 */

/**
 *  @addtogroup   utils
 *  @{
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
 *  @param  string        The string of elements encoded as natural numbers
 *  @param  string_length The length of the string
 *  @return               The order of rotational symmetry
 */
unsigned int
vrna_rotational_symmetry_num(const unsigned int *string,
                             size_t             string_length);


/**
 *  @brief Determine the order of rotational symmetry for a NULL-terminated string of ASCII characters
 *
 *  The algorithm applies a fast search of the provided string within itself, assuming
 *  the end of the string wraps around to connect with it's start.
 *  For example, a string of the form @p AABAAB has rotational symmetry
 *  of order @p 2
 *
 *  @param  string        A NULL-terminated string of characters
 *  @return               The order of rotational symmetry
 */
unsigned int
vrna_rotational_symmetry(const char *string);


/**
 *  @}
 */
#endif
