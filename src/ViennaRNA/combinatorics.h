#ifndef VIENNA_RNA_PACKAGE_COMBINATORICS_H
#define VIENNA_RNA_PACKAGE_COMBINATORICS_H

/**
 *  @file     combinatorics.h
 *  @ingroup  utils
 *  @brief     Various implementations that deal with combinatorial aspects of RNA/DNA folding
 */

/**
 *  @{
 *  @ingroup   utils
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
 *  @}
 */
#endif
