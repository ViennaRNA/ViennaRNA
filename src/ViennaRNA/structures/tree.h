#ifndef VIENNA_RNA_PACKAGE_STRUCTURES_TREE_H
#define VIENNA_RNA_PACKAGE_STRUCTURES_TREE_H

#ifdef VRNA_WARN_DEPRECATED
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
 *  @file     ViennaRNA/structures/tree.h
 *  @ingroup  struct_utils
 *  @brief    Various utility- and helper-functions for secondary structure parsing, converting, etc.
 */


/**
 *  @addtogroup struct_utils_tree
 *  @{
 */

/**
 *  @brief  Homeomorphically Irreducible Tree (HIT) representation of a secondary structure
 *  @see    vrna_db_to_tree_string()
 */
#define   VRNA_STRUCTURE_TREE_HIT             1U


/**
 *  @brief  (short) Coarse Grained representation of a secondary structure
 *  @see    vrna_db_to_tree_string()
 */
#define   VRNA_STRUCTURE_TREE_SHAPIRO_SHORT   2U


/**
 *  @brief  (full)  Coarse Grained representation of a secondary structure
 *  @see    vrna_db_to_tree_string()
 */
#define   VRNA_STRUCTURE_TREE_SHAPIRO         3U


/**
 *  @brief  (extended) Coarse Grained representation of a secondary structure
 *  @see    vrna_db_to_tree_string()
 */
#define   VRNA_STRUCTURE_TREE_SHAPIRO_EXT     4U


/**
 *  @brief  (weighted) Coarse Grained representation of a secondary structure
 *  @see    vrna_db_to_tree_string()
 */
#define   VRNA_STRUCTURE_TREE_SHAPIRO_WEIGHT  5U

/**
 *  @brief  Expanded Tree representation of a secondary structure
 *  @see    vrna_db_to_tree_string()
 */
#define   VRNA_STRUCTURE_TREE_EXPANDED        6U


/**
 *  @brief  Convert a Dot-Bracket structure string into tree string representation
 *
 *  This function allows one to convert a secondary structure in dot-bracket notation
 *  into one of the various tree representations for secondary structures. The resulting
 *  tree is then represented as a string of parenthesis and node symbols, similar to
 *  to the Newick format.
 *
 *  Currently we support conversion into the following formats, denoted by the value
 *  of parameter @p type:
 *  * #VRNA_STRUCTURE_TREE_HIT            - @copybrief #VRNA_STRUCTURE_TREE_HIT
 *                                          (See also @rstinline :cite:t:`fontana:1993b` @endrst)
 *  * #VRNA_STRUCTURE_TREE_SHAPIRO_SHORT  - @copybrief #VRNA_STRUCTURE_TREE_SHAPIRO_SHORT
 *                                          (same as @rstinline :cite:t:`shapiro:1988` @endrst, but with root node @p R and without @p S nodes for the stems)
 *  * #VRNA_STRUCTURE_TREE_SHAPIRO        - @copybrief #VRNA_STRUCTURE_TREE_SHAPIRO
 *                                          (See also @rstinline :cite:t:`shapiro:1988` @endrst)
 *  * #VRNA_STRUCTURE_TREE_SHAPIRO_EXT    - @copybrief #VRNA_STRUCTURE_TREE_SHAPIRO_EXT
 *                                          (same as @rstinline :cite:t:`shapiro:1988` @endrst, but external nodes denoted as @p E )
 *  * #VRNA_STRUCTURE_TREE_SHAPIRO_WEIGHT - @copybrief #VRNA_STRUCTURE_TREE_SHAPIRO_WEIGHT
 *                                          (same as #VRNA_STRUCTURE_TREE_SHAPIRO_EXT but with additional weights
 *                                          for number of unpaired nucleotides in loop, and number of pairs in stems)
 *  * #VRNA_STRUCTURE_TREE_EXPANDED       - @copybrief #VRNA_STRUCTURE_TREE_EXPANDED
 *
 *  @see  @ref sec_structure_representations_tree
 *
 *  @param  structure   The null-terminated dot-bracket structure string
 *  @param  type        A switch to determine the type of tree string representation
 *  @return             A tree representation of the input @p structure
 */
char *
vrna_db_to_tree_string(const char   *structure,
                       unsigned int type);


/**
 *  @brief  Remove weights from a linear string tree representation of a secondary structure
 *
 *  This function strips the weights of a linear string tree representation such as @p HIT,
 *  or Coarse Grained Tree sensu @rstinline :cite:t:`shapiro:1988` @endrst
 *
 *  @see vrna_db_to_tree_string()
 *
 *  @param  structure   A linear string tree representation of a secondary structure with weights
 *  @return             A linear string tree representation of a secondary structure without weights
 */
char *
vrna_tree_string_unweight(const char *structure);


/**
 *  @brief  Convert a linear tree string representation of a secondary structure back to Dot-Bracket notation
 *
 *  @warning  This function only accepts <em>Expanded</em> and <em>HIT</em> tree representations!
 *
 *  @see vrna_db_to_tree_string(), #VRNA_STRUCTURE_TREE_EXPANDED, #VRNA_STRUCTURE_TREE_HIT,
 *       @ref sec_structure_representations_tree
 *
 *  @param  tree  A linear tree string representation of a secondary structure
 *  @return       A dot-bracket notation of the secondary structure provided in @p tree
 */
char *
vrna_tree_string_to_db(const char *tree);


/* End tree representations */
/** @} */


#endif
