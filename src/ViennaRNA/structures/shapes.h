#ifndef VIENNA_RNA_PACKAGE_STRUCTURES_SHAPES_H
#define VIENNA_RNA_PACKAGE_STRUCTURES_SHAPES_H

/**
 *  @file     ViennaRNA/structures/shapes.h
 *  @ingroup  struct_utils
 *  @brief    Abstract shapes of secondary structures
 */

/**
 *  @addtogroup struct_utils_abstract_shapes
 *  @{
 */

/**
 *  @brief  Convert a secondary structure in dot-bracket notation to its abstract shapes representation
 *
 *  This function converts a secondary structure into its abstract shapes representation as
 *  presented by @rstinline :cite:t:`giegerich:2004` @endrst.
 *
 *  @see vrna_abstract_shapes_pt()
 *
 *  @param  structure A secondary structure in dot-bracket notation
 *  @param  level     The abstraction level (integer in the range of 0 to 5)
 *  @return           The secondary structure in abstract shapes notation
 */
char *
vrna_abstract_shapes(const char    *structure,
                     unsigned int  level);


/**
 *  @brief  Convert a secondary structure to its abstract shapes representation
 *
 *  This function converts a secondary structure into its abstract shapes representation as
 *  presented by @rstinline :cite:t:`giegerich:2004` @endrst. This function is equivalent to
 *  vrna_db_to_shapes(), but requires a pair table input instead of a dot-bracket structure.
 *
 *  @note   The length of the structure must be present at @p pt[0]!
 *
 *  @see vrna_abstract_shapes()
 *
 *  @param  pt      A secondary structure in pair table format
 *  @param  level   The abstraction level (integer in the range of 0 to 5)
 *  @return         The secondary structure in abstract shapes notation
 */
char *
vrna_abstract_shapes_pt(const short  *pt,
                        unsigned int level);


/* End abstract shapes interface */
/** @} */


#endif
