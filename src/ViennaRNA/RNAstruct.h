#ifndef VIENNA_RNA_PACKAGE_RNASTRUCT_H
#define VIENNA_RNA_PACKAGE_RNASTRUCT_H

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
 *  @addtogroup   struct_utils_deprecated
 *
 *  @{
 *
 *  @file RNAstruct.h
 *  @brief Parsing and Coarse Graining of Structures
 *
 *   Example:
 *  @verbatim
 *   .((..(((...)))..((..)))).   is the bracket or full tree
 *   becomes expanded:   - expand_Full() -
 *   ((U)(((U)(U)((((U)(U)(U)P)P)P)(U)(U)(((U)(U)P)P)P)P)(U)R)
 *   HIT:                - b2HIT() -
 *   ((U1)((U2)((U3)P3)(U2)((U2)P2)P2)(U1)R)
 *   Coarse:             - b2C() -
 *   ((H)((H)M)R)
 *   becomes expanded:   - expand_Shapiro() -
 *   (((((H)S)((H)S)M)S)R)
 *   weighted Shapiro:   - b2Shapiro() -
 *   ((((((H3)S3)((H2)S2)M4)S2)E2)R)
 *  @endverbatim
 */

#define STRUC     2000

/**
 *  @brief Converts the full structure from bracket notation to the HIT
 *  notation including root.
 *
 *  @deprecated   See vrna_db_to_tree_string() and #VRNA_STRUCTURE_TREE_HIT
 *                for a replacement
 *
 *  @param structure
 *  @return
 */
DEPRECATED(char *b2HIT(const char *structure),
           "Use vrna_db_to_tree_string() instead!");             /* Full   -> HIT    [incl. root] */


/**
 *  @brief Converts the full structure from bracket notation to the a
 *  coarse grained notation using the 'H' 'B' 'I' 'M' and 'R' identifiers.
 *
 *  @deprecated   See vrna_db_to_tree_string() and #VRNA_STRUCTURE_TREE_SHAPIRO_SHORT
 *                for a replacement
 *
 *  @param structure
 *  @return
 */
DEPRECATED(char *b2C(const char *structure),
           "Use vrna_db_to_tree_string() instead!");               /* Full   -> Coarse [incl. root] */


/**
 *  @brief Converts the full structure from bracket notation to the
 *  <i>weighted</i> coarse grained notation using the 'H' 'B' 'I' 'M' 'S' 'E' and
 *  'R' identifiers.
 *
 *  @deprecated   See vrna_db_to_tree_string() and #VRNA_STRUCTURE_TREE_SHAPIRO_WEIGHT
 *                for a replacement
 *
 *  @param structure
 *  @return
 */
DEPRECATED(char *b2Shapiro(const char *structure),
           "Use vrna_db_to_tree_string() instead!");         /* Full -> weighted Shapiro [i.r.] */


/**
 *  @brief Adds a root to an un-rooted tree in any except bracket notation.
 *
 *  @param  structure
 *  @return
 */
DEPRECATED(char *add_root(const char *structure),
           "");                   /* {Tree} -> ({Tree}R)          */


/**
 *  @brief Inserts missing 'S' identifiers in unweighted coarse grained structures
 *  as obtained from b2C().
 *
 *  @param coarse
 *  @return
 */
DEPRECATED(char *expand_Shapiro(const char *coarse),
           "Use vrna_db_to_tree_string() instead!");


/* add S for stacks to coarse struct */
/**
 *  @brief Convert the full structure from bracket notation to the
 *  expanded notation including root.
 *
 *  @param structure
 *  @return
 */
DEPRECATED(char *expand_Full(const char *structure),
           "Use vrna_db_to_tree_string() instead!");       /* Full   -> FFull         */


/**
 *  @brief Restores the bracket notation from an expanded full or HIT tree, that is
 *  any tree using only identifiers 'U' 'P' and 'R'.
 *
 *  @param ffull
 *  @return
 */
DEPRECATED(char *unexpand_Full(const char *ffull),
           "Use vrna_tree_string_to_db() instead!");         /* FFull  -> Full          */


/**
 *  @brief Strip weights from any weighted tree.
 *
 *  @param wcoarse
 *  @return
 */
DEPRECATED(char *unweight(const char *wcoarse),
           "Use vrna_tree_string_unweight() instead!");            /* remove weights from coarse struct */


/**
 *  @brief Converts two aligned structures in expanded notation.
 *
 *  Takes two aligned structures as produced by
 *  tree_edit_distance() function back to bracket notation with '_'
 *  as the gap character. The result overwrites the input.
 *
 *  @param align
 */
DEPRECATED(void   unexpand_aligned_F(char *align[2]),
           "");


/**
 *  @brief Collects a statistic of structure elements of the full structure in
 *  bracket notation.
 *
 *  The function writes to the following global variables:
 *  #loop_size, #loop_degree, #helix_size, #loops, #pairs, #unpaired
 *
 *  @param structure
 */
DEPRECATED(void   parse_structure(const char *structure),
           "");  /* make structure statistics */


/**
 *  @brief contains a list of all loop sizes. loop_size[0] contains the
 *  number of external bases.
 */
DEPRECATED(extern int loop_size[STRUC],
           "");         /* loop sizes of a structure */

/**
 *  @brief contains a list of all stack sizes.
 */
DEPRECATED(extern int helix_size[STRUC],
           "");        /* helix sizes of a structure */

/**
 *  @brief contains the corresponding list of loop degrees.
 */
DEPRECATED(extern int loop_degree[STRUC],
           "");       /* loop degrees of a structure */

/**
 *  @brief contains the number of loops ( and therefore of stacks ).
 */
DEPRECATED(extern int loops,
           "");                    /* n of loops and stacks */

/**
 *  @brief contains the number of unpaired bases.
 */
DEPRECATED(extern int unpaired,
           "");

/**
 *  @brief contains the number of base pairs in the last parsed structure.
 */
DEPRECATED(extern int pairs,
           "");          /* n of unpaired digits and pairs */

/**
 * @}
 */

#endif
