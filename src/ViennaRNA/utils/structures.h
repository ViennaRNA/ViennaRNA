#ifndef VIENNA_RNA_PACKAGE_STRUCT_UTILS_H
#define VIENNA_RNA_PACKAGE_STRUCT_UTILS_H

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
 *  @file     ViennaRNA/utils/structures.h
 *  @ingroup  struct_utils
 *  @brief    Various utility- and helper-functions for secondary structure parsing, converting, etc.
 */

/**
 *  @addtogroup   struct_utils
 *  @{
 *  @brief  Functions to create, parse, convert, manipulate, and compare secondary structure representations
 */


/**
 *  @brief Convenience typedef for data structure #vrna_hx_s
 *  @ingroup  struct_utils_helix_list
 */
typedef struct vrna_hx_s vrna_hx_t;


/**
 *  @brief Convenience typedef for data structure #vrna_elem_prob_s
 *  @ingroup  struct_utils_plist
 */
typedef struct vrna_elem_prob_s vrna_ep_t;


/**
 *  @addtogroup struct_utils_dot_bracket
 *  @{
 *  @brief  The Dot-Bracket notation as introduced already in the early times of the ViennaRNA Package
 *          denotes base pairs by matching pairs of parenthesis `()` and unpaired nucleotides by dots `.`.
 *
 *  As a simple example, consider a helix of size 4 enclosing a hairpin of size 4. In dot-bracket
 *  notation, this is annotated as
 *
 *  `((((....))))`
 *
 *  <b>Extended Dot-Bracket Notation</b>
 *
 *  A more generalized version of the original Dot-Bracket notation may use additional pairs
 *  of brackets, such as <tt><></tt>, <tt>{}</tt>, and <tt>[]</tt>, and matching pairs of
 *  uppercase/lowercase letters. This allows for anotating pseudo-knots, since different
 *  pairs of brackets are not required to be nested.
 *
 *  The follwing annotations of a simple structure with two crossing helices of size 4 are equivalent:
 *
 *  `<<<<[[[[....>>>>]]]]`<br>
 *  `((((AAAA....))))aaaa`<br>
 *  `AAAA{{{{....aaaa}}}}`
 */

/**
 *  @brief  Bitflag to indicate secondary structure notations using uppercase/lowercase letters from the latin alphabet
 *
 *  @see  vrna_ptable_from_string()
 */
#define VRNA_BRACKETS_ALPHA    4U


/**
 *  @brief  Bitflag to indicate secondary structure notations using round brackets (parenthesis), <tt>()</tt>
 *
 *  @see  vrna_ptable_from_string(), vrna_db_flatten(), vrna_db_flatten_to()
 */
#define VRNA_BRACKETS_RND      8U


/**
 *  @brief  Bitflag to indicate secondary structure notations using curly brackets, <tt>{}</tt>
 *
 *  @see  vrna_ptable_from_string(), vrna_db_flatten(), vrna_db_flatten_to()
 */
#define VRNA_BRACKETS_CLY      16U


/**
 *  @brief  Bitflag to indicate secondary structure notations using angular brackets, <tt><></tt>
 *
 *  @see  vrna_ptable_from_string(), vrna_db_flatten(), vrna_db_flatten_to()
 */
#define VRNA_BRACKETS_ANG      32U


/**
 *  @brief  Bitflag to indicate secondary structure notations using square brackets, <tt>[]</tt>
 *
 *  @see  vrna_ptable_from_string(), vrna_db_flatten(), vrna_db_flatten_to()
 */
#define VRNA_BRACKETS_SQR      64U


/**
 *  @brief  Default bitmask to indicate secondary structure notation using any pair of brackets
 *
 *  This set of matching brackets/parenthesis is always nested, i.e. pseudo-knot free, in WUSS
 *  format. However, in general different kinds of brackets are mostly used for annotating
 *  pseudo-knots. Thus special care has to be taken to remove pseudo-knots if this bitmask
 *  is used in functions that return secondary structures without pseudo-knots!
 *
 *  @see  vrna_ptable_from_string(), vrna_db_flatten(), vrna_db_flatten_to(), vrna_db_pk_remove()
 *        vrna_pt_pk_remove()
 */
#define VRNA_BRACKETS_DEFAULT  \
  (VRNA_BRACKETS_RND | \
   VRNA_BRACKETS_CLY | \
   VRNA_BRACKETS_ANG | \
   VRNA_BRACKETS_SQR)


/**
 *  @brief  Bitmask to indicate secondary structure notation using any pair of brackets or uppercase/lowercase alphabet letters
 *
 *  @see  vrna_ptable_from_string(), vrna_db_pk_remove(), vrna_db_flatten(),
 *        vrna_db_flatten_to()
 */
#define VRNA_BRACKETS_ANY \
  (VRNA_BRACKETS_RND | \
   VRNA_BRACKETS_CLY | \
   VRNA_BRACKETS_ANG | \
   VRNA_BRACKETS_SQR | \
   VRNA_BRACKETS_ALPHA)


#include <stdio.h>

#include <ViennaRNA/datastructures/basic.h>

/**
 *  @brief Pack secondary secondary structure, 5:1 compression using base 3 encoding
 *
 *  Returns a binary string encoding of the secondary structure using
 *  a 5:1 compression scheme. The string is NULL terminated and can
 *  therefore be used with standard string functions such as strcmp().
 *  Useful for programs that need to keep many structures in memory.
 *
 *  @see  vrna_db_unpack()
 *  @param struc    The secondary structure in dot-bracket notation
 *  @return         The binary encoded structure
 */
char *
vrna_db_pack(const char *struc);


/**
 *  @brief Unpack secondary structure previously packed with vrna_db_pack()
 *
 *  Translate a compressed binary string produced by vrna_db_pack() back into
 *  the familiar dot-bracket notation.
 *
 *  @see  vrna_db_pack()
 *  @param packed   The binary encoded packed secondary structure
 *  @return         The unpacked secondary structure in dot-bracket notation
 */
char *
vrna_db_unpack(const char *packed);


/**
 *  @brief Substitute pairs of brackets in a string with parenthesis
 *
 *  This function can be used to replace brackets of unusual types,
 *  such as angular brackets @p <> , to dot-bracket format.
 *  The @p options parameter is used tpo specify which types of brackets
 *  will be replaced by round parenthesis @p () .
 *
 *  @see vrna_db_flatten_to(),
 *       #VRNA_BRACKETS_RND, #VRNA_BRACKETS_ANG, #VRNA_BRACKETS_CLY, #VRNA_BRACKETS_SQR,
 *       #VRNA_BRACKETS_DEFAULT
 *
 *  @param  structure   The structure string where brackets are flattened in-place
 *  @param  options     A bitmask to specify which types of brackets should be flattened out
 */
void
vrna_db_flatten(char          *structure,
                unsigned int  options);


/**
 *  @brief Substitute pairs of brackets in a string with another type of pair characters
 *
 *  This function can be used to replace brackets in a structure annotation string,
 *  such as square brackets @p [] , to another type of pair characters,
 *  e.g. angular brackets @p <> .
 *
 *  The @p target array must contain a character for the 'pair open' annotation at
 *  position 0, and one for 'pair close' at position 1. T@p options parameter is used
 *  to specify which types of brackets will be replaced by the new pairs.
 *
 *  @see vrna_db_flatten(),
 *       #VRNA_BRACKETS_RND, #VRNA_BRACKETS_ANG, #VRNA_BRACKETS_CLY, #VRNA_BRACKETS_SQR,
 *       #VRNA_BRACKETS_DEFAULT
 *
 *  @param  string      The structure string where brackets are flattened in-place
 *  @param  target      The new pair characters the string will be flattened to
 *  @param  options     A bitmask to specify which types of brackets should be flattened out
 */
void
vrna_db_flatten_to(char         *string,
                   const char   target[3],
                   unsigned int options);


/**
 *  @brief Convert a pair table into dot-parenthesis notation
 *
 *  @param pt The pair table to be copied
 *  @return   A char pointer to the dot-bracket string
 */
char *
vrna_db_from_ptable(short *pt);


/**
 *  @brief  Convert a list of base pairs into dot-bracket notation
 *
 *  @see vrna_plist()
 *  @param  pairs   A #vrna_ep_t containing the pairs to be included in
 *                  the dot-bracket string
 *  @param  n       The length of the structure (number of nucleotides)
 *  @return         The dot-bracket string containing the provided base pairs
 */
char *
vrna_db_from_plist(vrna_ep_t    *pairs,
                   unsigned int n);


/**
 *  @brief  Convert a secondary structure in dot-bracket notation to a nucleotide annotation of loop contexts
 *
 *  @param  structure   The secondary structure in dot-bracket notation
 *  @return             A string annotating each nucleotide according to it's structural context
 */
char *
vrna_db_to_element_string(const char *structure);


/**
 *  @brief  Remove pseudo-knots from an input structure
 *
 *  This function removes pseudo-knots from an input structure
 *  by determining the minimum number of base pairs that need
 *  to be removed to make the structure pseudo-knot free.
 *
 *  To accomplish that, we use a dynamic programming algorithm
 *  similar to the Nussinov maxmimum matching approach.
 *
 *  The input structure must be in a dot-bracket string like form
 *  where crossing base pairs are denoted by the use of additional
 *  types of matching brackets, e.g. @p <>, @p {}, @p [], @p {}.
 *  Furthermore, crossing pairs may be annotated by matching
 *  uppercase/lowercase letters from the alphabet @p A-Z. For the latter,
 *  the uppercase letter must be the 5' and the lowercase letter
 *  the 3' nucleotide of the base pair. The actual type of brackets
 *  to be recognized by this function must be specifed through the
 *  @p options parameter.
 *
 *  @note Brackets in the input structure string that are not covered
 *        by the @p options bitmask will be silently ignored!
 *
 *  @see vrna_pt_pk_remove(), vrna_db_flatten(),
 *       #VRNA_BRACKETS_RND, #VRNA_BRACKETS_ANG, #VRNA_BRACKETS_CLY, #VRNA_BRACKETS_SQR,
 *       #VRNA_BRACKETS_ALPHA, #VRNA_BRACKETS_DEFAULT, #VRNA_BRACKETS_ANY
 *
 *  @param  structure   Input structure in dot-bracket format that may include pseudo-knots
 *  @param  options     A bitmask to specify which types of brackets should be processed
 *  @return             The input structure devoid of pseudo-knots in dot-bracket notation
 */
char *
vrna_db_pk_remove(const char *structure,
                  unsigned int options);

/* End dot-bracket interface */
/**@}*/

/**
 *  @addtogroup struct_utils_pair_table
 *  @{
 */

/**
 *  @brief Create a pair table from a dot-bracket notation of a secondary structure
 *
 *  Returns a newly allocated table, such that table[i]=j if (i.j) pair
 *  or 0 if i is unpaired, table[0] contains the length of the structure.
 *
 *  @see  vrna_ptable_from_string(), vrna_db_from_ptable()
 *
 *  @param  structure The secondary structure in dot-bracket notation
 *  @return           A pointer to the created pair_table
 */
short *
vrna_ptable(const char *structure);


/**
 *  @brief  Create a pair table for a secondary structure string
 *
 *  This function takes an input string of a secondary structure annotation
 *  in @ref dot-bracket-notation or @ref dot-bracket-ext-notation, and converts
 *  it into a pair table representation.
 *
 *  @note   This function also extracts crossing base pairs, i.e. pseudo-knots
 *          if more than a single matching bracket type is allowed through the
 *          bitmask @p options.
 *
 *  @see vrna_ptable(), vrna_db_from_ptable(), vrna_db_flatten_to(), vrna_pt_pk_remove()
 *       #VRNA_BRACKETS_RND, #VRNA_BRACKETS_ANG, #VRNA_BRACKETS_CLY, #VRNA_BRACKETS_SQR,
 *       VRNA_BRACKETS_ALPHA, #VRNA_BRACKETS_DEFAULT, #VRNA_BRACKETS_ANY
 *
 *  @param  string    Secondary structure in @ref dot-bracket-ext-notation
 *  @param  options   A bitmask to specify which brackets are recognized during conversion to pair table
 *  @return           A pointer to a new pair table of the provided secondary structure
 */
short *
vrna_ptable_from_string(const char    *string,
                        unsigned int  options);


/**
 *  @brief Create a pair table of a secondary structure (pseudo-knot version)
 *
 *  Returns a newly allocated table, such that table[i]=j if (i.j) pair
 *  or 0 if i is unpaired, table[0] contains the length of the structure.
 *
 *  In contrast to vrna_ptable() this function also recognizes the base pairs
 *  denoted by '[' and ']' brackets. Thus, this function behaves like
 *  @code{.c}
 *  vrna_ptable_from_string(structure, #VRNA_BRACKETS_RND | VRNA_BRACKETS_SQR)
 *  @endcode
 *
 *  @see    vrna_ptable_from_string()
 *
 *  @param  structure The secondary structure in (extended) dot-bracket notation
 *  @return           A pointer to the created pair_table
 */
short *
vrna_pt_pk_get(const char *structure);


/**
 *  @brief Get an exact copy of a pair table
 *
 *  @param pt The pair table to be copied
 *  @return   A pointer to the copy of 'pt'
 */
short *
vrna_ptable_copy(const short *pt);


/**
 * @brief Create a pair table of a secondary structure (snoop align version)
 *
 */
short *
vrna_pt_ali_get(const char *structure);


/**
 * @brief Create a pair table of a secondary structure (snoop version)
 *
 *  returns a newly allocated table, such that:  table[i]=j if (i.j) pair or
 *  0 if i is unpaired, table[0] contains the length of the structure.
 *  The special pseudoknotted H/ACA-mRNA structure is taken into account.
 */
short *
vrna_pt_snoop_get(const char *structure);


/**
 *  @brief  Remove pseudo-knots from a pair table
 *
 *  This function removes pseudo-knots from an input structure
 *  by determining the minimum number of base pairs that need
 *  to be removed to make the structure pseudo-knot free.
 *
 *  To accomplish that, we use a dynamic programming algorithm
 *  similar to the Nussinov maxmimum matching approach.
 *
 *  @see    vrna_db_pk_remove()
 *
 *  @param  ptable  Input structure that may include pseudo-knots
 *  @param  options
 *  @return         The input structure devoid of pseudo-knots
 */
short *
vrna_pt_pk_remove(const short   *ptable,
                  unsigned int  options);


/* End pair table interface */
/**@}*/


/**
 *  @addtogroup struct_utils_plist
 *  @{
 */

/**
 *  @brief  A Base Pair element
 */
#define VRNA_PLIST_TYPE_BASEPAIR      0


/**
 *  @brief  A G-Quadruplex element
 */
#define VRNA_PLIST_TYPE_GQUAD         1


/**
 *  @brief  A Hairpin loop motif element
 */
#define VRNA_PLIST_TYPE_H_MOTIF       2


/**
 *  @brief  An Internal loop motif element
 */
#define VRNA_PLIST_TYPE_I_MOTIF       3


/**
 *  @brief  An Unstructured Domain motif element
 */
#define VRNA_PLIST_TYPE_UD_MOTIF      4


/**
 *  @brief  A Base Pair stack element
 */
#define VRNA_PLIST_TYPE_STACK         5


/**
 *  @brief  Data structure representing a single entry of an element probability list
 *          (e.g. list of pair probabilities)
 *
 *  @see vrna_plist(), vrna_plist_from_probs(), vrna_db_from_plist(),
 *  #VRNA_PLIST_TYPE_BASEPAIR, #VRNA_PLIST_TYPE_GQUAD, #VRNA_PLIST_TYPE_H_MOTIF, #VRNA_PLIST_TYPE_I_MOTIF,
 *  #VRNA_PLIST_TYPE_UD_MOTIF, #VRNA_PLIST_TYPE_STACK
 */
struct vrna_elem_prob_s {
  int   i;    /**<  @brief  Start position (usually 5' nucleotide that starts the element, e.g. base pair) */
  int   j;    /**<  @brief  End position (usually 3' nucleotide that ends the element, e.g. base pair) */
  float p;    /**<  @brief  Probability of the element */
  int   type; /**<  @brief  Type of the element */
};

/**
 *  @brief Create a #vrna_ep_t from a dot-bracket string
 *
 *  The dot-bracket string is parsed and for each base pair an
 *  entry in the plist is created. The probability of each pair in
 *  the list is set by a function parameter.
 *
 *  The end of the plist is marked by sequence positions i as well as j
 *  equal to 0. This condition should be used to stop looping over its
 *  entries
 *
 *  @param struc  The secondary structure in dot-bracket notation
 *  @param pr     The probability for each base pair used in the plist
 *  @return       The plist array
 */
vrna_ep_t *vrna_plist(const char  *struc,
                      float       pr);


/**
 *  @brief Create a #vrna_ep_t from base pair probability matrix
 *
 *  The probability matrix provided via the #vrna_fold_compound_t is parsed
 *  and all pair probabilities above the given threshold are used to create
 *  an entry in the plist
 *
 *  The end of the plist is marked by sequence positions i as well as j
 *  equal to 0. This condition should be used to stop looping over its
 *  entries
 *
 *  @ingroup              part_func_global
 *  @param[in]  vc        The fold compound
 *  @param[in]  cut_off   The cutoff value
 *  @return               A pointer to the plist that is to be created
 */
vrna_ep_t *vrna_plist_from_probs(vrna_fold_compound_t *vc,
                                 double               cut_off);


/* End pair list interface */
/**@}*/


/**
 *  @addtogroup struct_utils_wuss
 *  @{
 *  @brief  The WUSS notation, as frequently used for consensus secondary structures in @ref msa-formats-stockholm.
 *
 *  This notation allows for a fine-grained annotation of base pairs and unpaired nucleotides, including pseudo-knots.
 *  Below, you'll find a list of secondary structure elements and their corresponding WUSS annotation
 *  (See also the infernal user guide at http://eddylab.org/infernal/Userguide.pdf)
 *  @parblock
 *  - <b>Base pairs</b><br>
 *    Nested base pairs are annotated by matching pairs of the symbols `<>`,
 *    `()`, `{}`, and `[]`. Each of the matching pairs
 *    of parenthesis have their special meaning, however, when used as input in our programs,
 *    e.g. structure constraint, these details are usually ignored. Furthermore, base pairs
 *    that constitute as pseudo-knot are denoted by letters from the latin alphabet and are,
 *    if not denoted otherwise, ignored entirely in our programs.
 *
 *  - <b>Hairpin loops</b><br>
 *    Unpaired nucleotides that constitute the hairpin loop are indicated by underscores, `_`.
 *
 *    Example: `<<<<<_____>>>>>`
 *
 *  - <b>Bulges and interior loops</b><br>
 *    Residues that constitute a bulge or interior loop are denoted by dashes, `-`.
 *
 *    Example: `(((--<<_____>>-)))`
 *
 *  - <b>Multibranch loops</b><br>
 *    Unpaired nucleotides in multibranch loops are indicated by commas `,`.
 *  
 *    Example: `(((,,<<_____>>,<<____>>)))`
 *
 *  - <b>External residues</b><br>
 *    Single stranded nucleotides in the exterior loop, i.e. not enclosed by any other pair are
 *    denoted by colons, `:`.
 *
 *    Example: `<<<____>>>:::`
 *
 *  - <b>Insertions</b><br>
 *    In cases where an alignment represents the consensus with a known structure, insertions relative
 *    to the known structure are denoted by periods, `.`. Regions where local structural
 *    alignment was invoked, leaving regions of both target and query sequence unaligned, are indicated
 *    by tildes, `~`.
 *    @note These symbols only appear in alignments of a known (query) structure annotation to a target
 *    sequence of unknown structure.
 *
 *  - <b>Pseudo-knots</b><br>
 *    The WUSS notation allows for annotation of pseudo-knots using pairs of upper-case/lower-case letters.
 *    @note Our programs and library functions usually ignore pseudo-knots entirely treating them as
 *    unpaired nucleotides, if not stated otherwise.
 *
 *    Example:  `<<<_AAA___>>>aaa`
 *  @endparblock
 */

/**
 *  @brief  Convert a WUSS annotation string to dot-bracket format
 *
 *  @note This function flattens all brackets, and treats pseudo-knots annotated
 *        by matching pairs of upper/lowercase letters as unpaired nucleotides
 *
 *  @param  wuss  The input string in WUSS notation
 *  @return       A dot-bracket notation of the input secondary structure
 */
char *
vrna_db_from_WUSS(const char *wuss);


/* End WUSS notation interface */
/**@}*/


/**
 *  @addtogroup struct_utils_abstract_shapes
 *  @{
 *  @brief  Abstract Shapes, introduced by Giegerich et al. in (2004) @cite giegerich:2004,
 *          collapse the secondary structure while retaining the nestedness of helices and
 *          hairpin loops.
 *
 *  The abstract shapes representation abstracts the structure from individual base pairs
 *  and their corresponding location in the sequence, while retaining the inherent nestedness
 *  of helices and hairpin loops.
 *
 *  Below is a description of what is included in the abstract shapes abstraction for each
 *  respective level together with an example structure:
 *
 *      CGUCUUAAACUCAUCACCGUGUGGAGCUGCGACCCUUCCCUAGAUUCGAAGACGAG
 *      ((((((...(((..(((...))))))...(((..((.....))..)))))))))..
 *
 *  ______
 *
 *  Shape Level | Description                     |   Result
 *  ----------- | ------------------------------- | --------
 *  1           | Most accurate - all loops and all unpaired | `[_[_[]]_[_[]_]]_`
 *  2           | Nesting pattern for all loop types and unpaired regions in external loop and multiloop | `[[_[]][_[]_]]`
 *  3           | Nesting pattern for all loop types but no unpaired regions | `[[[]][[]]]`
 *  4           | Helix nesting pattern in external loop and multiloop | `[[][[]]]`
 *  5           | Most abstract - helix nesting pattern and no unpaired regions | `[[][]]`
 *
 *  @note   Our implementations also provide the special Shape Level 0, which does not
 *          collapse any structural features but simply convert base pairs and unpaired
 *          nucleotides into their corresponding set of symbols for abstract shapes.
 */

/**
 *  @brief  Convert a secondary structure in dot-bracket notation to its abstract shapes representation
 *
 *  This function converts a secondary structure into its abstract shapes representation as
 *  presented by Giegerich et al. 2004 @cite giegerich:2004.
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
 *  presented by Giegerich et al. 2004 @cite giegerich:2004. This function is equivalent to
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
/**@}*/


/**
 *  @addtogroup struct_utils_helix_list
 *  @{
 */

/**
 *  @brief  Data structure representing an entry of a helix list
 */
struct vrna_hx_s {
  unsigned int  start;
  unsigned int  end;
  unsigned int  length;
  unsigned int  up5;
  unsigned int  up3;
};


/**
 *  @brief  Convert a pair table representation of a secondary structure into a helix list
 *
 *  @param  pt  The secondary structure in pair table representation
 *  @return     The secondary structure represented as a helix list
 */
vrna_hx_t *
vrna_hx_from_ptable(short *pt);


/**
 *  @brief  Create a merged helix list from another helix list
 */
vrna_hx_t *
vrna_hx_merge(const vrna_hx_t *list,
              int             maxdist);


/* End helix list interface */
/**@}*/


/**
 *  @brief Get a loop index representation of a structure
 */
int *
vrna_loopidx_from_ptable(const short *pt);


/**
 *  @brief Compute the "base pair" distance between two pair tables pt1 and pt2 of secondary structures.
 *
 *  The pair tables should have the same length.
 *  dist = number of base pairs in one structure but not in the other
 *  same as edit distance with open-pair close-pair as move-set
 *
 *  @see vrna_bp_distance()
 *
 *  @param pt1   First structure in dot-bracket notation
 *  @param pt2   Second structure in dot-bracket notation
 *  @return       The base pair distance between pt1 and pt2
 */
int
vrna_bp_distance_pt(const short *pt1,
                    const short *pt2);

/**
 *  @brief Compute the "base pair" distance between two secondary structures s1 and s2.
 *
 *  This is a wrapper around @b vrna_bp_distance_pt().
 *  The sequences should have the same length.
 *  dist = number of base pairs in one structure but not in the other
 *  same as edit distance with open-pair close-pair as move-set
 *
 *  @see vrna_bp_distance_pt()
 *
 *  @param str1   First structure in dot-bracket notation
 *  @param str2   Second structure in dot-bracket notation
 *  @return       The base pair distance between str1 and str2
 */
int
vrna_bp_distance(const char *str1,
                 const char *str2);


double
vrna_dist_mountain(const char   *str1,
                   const char   *str2,
                   unsigned int p);


/**
 *  @brief Make a reference base pair count matrix
 *
 *  Get an upper triangular matrix containing the number of basepairs of a reference
 *  structure for each interval [i,j] with i<j. Access it via iindx!!!
 */
unsigned int *
vrna_refBPcnt_matrix(const short  *reference_pt,
                     unsigned int turn);


/**
 *  @brief Make a reference base pair distance matrix
 *
 *  Get an upper triangular matrix containing the base pair distance of two
 *  reference structures for each interval [i,j] with i<j. Access it via iindx!!!
 *
 */
unsigned int *
vrna_refBPdist_matrix(const short   *pt1,
                      const short   *pt2,
                      unsigned int  turn);


/**
 *  @brief Create a dot-bracket like structure string from base pair probability matrix
 */
char *
vrna_db_from_probs(const FLT_OR_DBL *pr,
                   unsigned int     length);


/**
 *  @brief Get a pseudo dot bracket notation for a given probability information
 */
char
vrna_bpp_symbol(const float *x);


/**
 *  @brief Create a dot-backet/parenthesis structure from backtracking stack
 *
 *  This function is capable to create dot-bracket structures from suboptimal
 *  structure prediction sensu M. Zuker
 *
 *  @param bp     Base pair stack containing the traced base pairs
 *  @param length The length of the structure
 *  @return       The secondary structure in dot-bracket notation as
 *                provided in the input
 */
char *
vrna_db_from_bp_stack(vrna_bp_stack_t *bp,
                      unsigned int    length);


void
vrna_letter_structure(char            *structure,
                      vrna_bp_stack_t *bp,
                      unsigned int    length);


/**
 *  @addtogroup struct_utils_tree
 *  @{
 *  @brief Secondary structures can be readily represented as trees, where internal
 *  nodes represent base pairs, and leaves represent unpaired nucleotides.
 *  The dot-bracket structure string already is a tree represented by a string
 *  of parenthesis (base pairs) and dots for the leaf nodes (unpaired nucleotides).
 *
 *  Alternatively, one may find representations with two types of node labels,
 *  `P` for paired and `U` for unpaired; a dot is then replaced by `(U)`, and
 *  each closed bracket is assigned an additional identifier `P`.
 *  We call this the expanded notation. In @cite fontana:1993b a condensed
 *  representation of the secondary structure is proposed, the so-called
 *  homeomorphically irreducible tree (HIT) representation. Here a stack is
 *  represented as a single pair of matching brackets labeled `P` and
 *  weighted by the number of base pairs.  Correspondingly, a contiguous
 *  strain of unpaired bases is shown as one pair of matching brackets
 *  labeled `U` and weighted by its length.  Generally any string consisting
 *  of matching brackets and identifiers is equivalent to a plane tree with
 *  as many different types of nodes as there are identifiers.
 *  
 *  Bruce Shapiro proposed a coarse grained representation @cite shapiro:1988,
 *  which, does not retain the full information of the secondary structure. He
 *  represents the different structure elements by single matching brackets
 *  and labels them as
 *  
 *  - `H`  (hairpin loop),
 *  - `I`  (interior loop),
 *  - `B`  (bulge),
 *  - `M`  (multi-loop), and
 *  - `S`  (stack).
 *  
 *  We extend his alphabet by an extra letter for external elements `E`.
 *  Again these identifiers may be followed by a weight corresponding to the
 *  number of unpaired bases or base pairs in the structure element.  All tree
 *  representations (except for the dot-bracket form) can be encapsulated into
 *  a virtual root (labeled `R`).
 *  
 *  The following example illustrates the different linear tree representations
 *  used by the package:
 *  
 *  Consider the secondary structure represented by the dot-bracket string (full tree)
 *  `.((..(((...)))..((..)))).` which is the most convenient
 *  condensed notation used by our programs and library functions.
 *  
 *  Then, the following tree representations are equivalent:
 *  
 *  - Expanded tree:<br>
 *    `((U)(((U)(U)((((U)(U)(U)P)P)P)(U)(U)(((U)(U)P)P)P)P)(U)R)`
 *  - HIT representation (Fontana et al. 1993 @cite fontana:1993b):<br>
 *    `((U1)((U2)((U3)P3)(U2)((U2)P2)P2)(U1)R)`
 *  - Coarse Grained Tree Representation (Shapiro 1988 @cite shapiro:1988):
 *    + Short (with root node `R`, without stem nodes `S`):<br>
 *      `((H)((H)M)R)`
 *    + Full (with root node `R`):<br>
 *      `(((((H)S)((H)S)M)S)R)`
 *    + Extended (with root node `R`, with external nodes `E`):<br>
 *      `((((((H)S)((H)S)M)S)E)R)`
 *    + Weighted (with root node `R`, with external nodes `E`):<br>
 *      `((((((H3)S3)((H2)S2)M4)S2)E2)R)`
 *
 *  The Expanded tree is rather clumsy and mostly included for the sake of
 *  completeness. The different versions of Coarse Grained Tree Representations
 *  are variatios of Shapiro's linear tree notation.
 *  
 *  For the output of aligned structures from string editing, different
 *  representations are needed, where we put the label on both sides.
 *  The above examples for tree representations would then look like:
 *  
 *  @verbatim
 *  a) (UU)(P(P(P(P(UU)(UU)(P(P(P(UU)(UU)(UU)P)P)P)(UU)(UU)(P(P(UU)(U...
 *  b) (UU)(P2(P2(U2U2)(P2(U3U3)P3)(U2U2)(P2(U2U2)P2)P2)(UU)P2)(UU)
 *  c) (B(M(HH)(HH)M)B)
 *     (S(B(S(M(S(HH)S)(S(HH)S)M)S)B)S)
 *     (E(S(B(S(M(S(HH)S)(S(HH)S)M)S)B)S)E)
 *  d) (R(E2(S2(B1(S2(M4(S3(H3)S3)((H2)S2)M4)S2)B1)S2)E2)R)
 *  @endverbatim
 *
 *  Aligned structures additionally contain the gap character `_`.
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
 *                                          (See also Fontana et al. 1993 @cite fontana:1993b)
 *  * #VRNA_STRUCTURE_TREE_SHAPIRO_SHORT  - @copybrief #VRNA_STRUCTURE_TREE_SHAPIRO_SHORT
 *                                          (same as Shapiro 1988 @cite shapiro:1988, but with root node @p R and without @p S nodes for the stems)
 *  * #VRNA_STRUCTURE_TREE_SHAPIRO        - @copybrief #VRNA_STRUCTURE_TREE_SHAPIRO
 *                                          (See also Shapiro 1988 @cite shapiro:1988)
 *  * #VRNA_STRUCTURE_TREE_SHAPIRO_EXT    - @copybrief #VRNA_STRUCTURE_TREE_SHAPIRO_EXT
 *                                          (same as Shapiro 1988 @cite shapiro:1988, but external nodes denoted as @p E )
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
 *  or Coarse Grained Tree sensu Shapiro @cite shapiro:1988
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
/**@}*/

#ifndef VRNA_DISABLE_BACKWARD_COMPATIBILITY

/*###########################################*/
/*# deprecated functions below              #*/
/*###########################################*/

/**
 *  @brief Create a #vrna_ep_t from a dot-bracket string
 *
 *  The dot-bracket string is parsed and for each base pair an
 *  entry in the plist is created. The probability of each pair in
 *  the list is set by a function parameter.
 *
 *  The end of the plist is marked by sequence positions i as well as j
 *  equal to 0. This condition should be used to stop looping over its
 *  entries
 *
 *  @deprecated   Use vrna_plist() instead
 *
 *  @ingroup part_func_global_deprecated
 *
 *  @param pl     A pointer to the #vrna_ep_t that is to be created
 *  @param struc  The secondary structure in dot-bracket notation
 *  @param pr     The probability for each base pair
 */
DEPRECATED(void assign_plist_from_db(vrna_ep_t  **pl,
                                     const char *struc,
                                     float      pr),
           "Use vrna_plist() instead");

/**
 *  @brief Pack secondary secondary structure, 5:1 compression using base 3 encoding
 *
 *  Returns a binary string encoding of the secondary structure using
 *  a 5:1 compression scheme. The string is NULL terminated and can
 *  therefore be used with standard string functions such as strcmp().
 *  Useful for programs that need to keep many structures in memory.
 *
 *  @deprecated     Use vrna_db_pack() as a replacement
 *  @ingroup        struct_utils_deprecated
 *  @param struc    The secondary structure in dot-bracket notation
 *  @return         The binary encoded structure
 */
DEPRECATED(char *pack_structure(const char *struc),
           "Use vrna_db_pack() instead");

/**
 *  @brief Unpack secondary structure previously packed with pack_structure()
 *
 *  Translate a compressed binary string produced by pack_structure() back into
 *  the familiar dot-bracket notation.
 *
 *  @deprecated     Use vrna_db_unpack() as a replacement
 *  @ingroup        struct_utils_deprecated
 *  @param packed   The binary encoded packed secondary structure
 *  @return         The unpacked secondary structure in dot-bracket notation
 */
DEPRECATED(char *unpack_structure(const char *packed),
           "Use vrna_db_unpack() instead");

/**
 *  @brief Create a pair table of a secondary structure
 *
 *  Returns a newly allocated table, such that table[i]=j if (i.j) pair
 *  or 0 if i is unpaired, table[0] contains the length of the structure.
 *
 *  @deprecated Use vrna_ptable() instead
 *  @ingroup        struct_utils_deprecated
 *
 *  @param  structure The secondary structure in dot-bracket notation
 *  @return           A pointer to the created pair_table
 */
DEPRECATED(short *make_pair_table(const char *structure),
           "Use vrna_ptable() instead");

DEPRECATED(short *make_pair_table_pk(const char *structure),
           "Use vrna_ptable_from_string() instead");

/**
 *  @brief Get an exact copy of a pair table
 *
 *  @deprecated Use vrna_ptable_copy() instead
 *  @ingroup        struct_utils_deprecated
 *
 *  @param pt The pair table to be copied
 *  @return   A pointer to the copy of 'pt'
 */
DEPRECATED(short *copy_pair_table(const short *pt),
           "Use vrna_ptable_copy() instead");

/**
 *  Pair table for snoop align
 *
 *  @deprecated Use vrna_pt_ali_get() instead!
 *  @ingroup        struct_utils_deprecated
 */
DEPRECATED(short *alimake_pair_table(const char *structure),
           "Use vrna_pt_ali_get() instead");

/**
 *  returns a newly allocated table, such that:  table[i]=j if (i.j) pair or
 *  0 if i is unpaired, table[0] contains the length of the structure.
 *  The special pseudoknotted H/ACA-mRNA structure is taken into account.
 *  @deprecated Use vrna_pt_snoop_get() instead!
 *  @ingroup        struct_utils_deprecated
 */
DEPRECATED(short *make_pair_table_snoop(const char *structure),
           "Use vrna_pt_snoop_get() instead");

DEPRECATED(int *make_loop_index_pt(short *pt),
           "Use vrna_loopidx_from_ptable() instead");

/**
 *  @brief Compute the "base pair" distance between two secondary structures s1 and s2.
 *
 *  The sequences should have the same length.
 *  dist = number of base pairs in one structure but not in the other
 *  same as edit distance with open-pair close-pair as move-set
 *
 *  @deprecated   Use vrna_bp_distance instead
 *  @ingroup        struct_utils_deprecated
 *  @param str1   First structure in dot-bracket notation
 *  @param str2   Second structure in dot-bracket notation
 *  @return       The base pair distance between str1 and str2
 */
DEPRECATED(int bp_distance(const char *str1,
                           const char *str2),
           "Use vrna_bp_distance() instead");

/**
 *  @brief Make a reference base pair count matrix
 *
 *  Get an upper triangular matrix containing the number of basepairs of a reference
 *  structure for each interval [i,j] with i<j. Access it via iindx!!!
 *
 *  @deprecated Use vrna_refBPcnt_matrix() instead
 *  @ingroup        struct_utils_deprecated
 */
DEPRECATED(unsigned int *make_referenceBP_array(short         *reference_pt,
                                                unsigned int  turn),
           "Use vrna_refBPcnt_matrix() instead");

/**
 *  @brief Make a reference base pair distance matrix
 *
 *  Get an upper triangular matrix containing the base pair distance of two
 *  reference structures for each interval [i,j] with i<j. Access it via iindx!!!
 *
 *  @deprecated Use vrna_refBPdist_matrix() instead
 *  @ingroup        struct_utils_deprecated
 */
DEPRECATED(unsigned int *compute_BPdifferences(short        *pt1,
                                               short        *pt2,
                                               unsigned int turn),
           "Use vrna_refBPdist_matrix() instead");

/**
 *  @brief Create a vrna_ep_t from a probability matrix
 *
 *  The probability matrix given is parsed and all pair probabilities above
 *  the given threshold are used to create an entry in the plist
 *
 *  The end of the plist is marked by sequence positions i as well as j
 *  equal to 0. This condition should be used to stop looping over its
 *  entries
 *
 *  @note This function is threadsafe
 *  @deprecated Use vrna_plist_from_probs() instead!
 *
 *  @ingroup part_func_global_deprecated
 *
 *  @param[out] pl      A pointer to the vrna_ep_t that is to be created
 *  @param[in]  probs   The probability matrix used for creating the plist
 *  @param[in]  length  The length of the RNA sequence
 *  @param[in]  cutoff  The cutoff value
 */
DEPRECATED(void  assign_plist_from_pr(vrna_ep_t   **pl,
                                      FLT_OR_DBL  *probs,
                                      int         length,
                                      double      cutoff),
           "Use vrna_plist_from_probs() instead");

/**
 *  @brief Create a dot-backet/parenthesis structure from backtracking stack
 *
 *  @deprecated use vrna_parenthesis_structure() instead
 *  @ingroup        struct_utils_deprecated
 *
 *  @note This function is threadsafe
 */
DEPRECATED(void parenthesis_structure(char            *structure,
                                      vrna_bp_stack_t *bp,
                                      int             length),
           "Use vrna_parenthesis_structure() instead");

/**
 *  @brief Create a dot-backet/parenthesis structure from backtracking stack
 *  obtained by zuker suboptimal calculation in cofold.c
 *
 *  @deprecated use vrna_parenthesis_zuker instead
 *  @ingroup        struct_utils_deprecated
 *
 *  @note This function is threadsafe
 */
DEPRECATED(void parenthesis_zuker(char            *structure,
                                  vrna_bp_stack_t *bp,
                                  int             length),
           "Use vrna_parenthesis_zuker() instead");

DEPRECATED(void letter_structure(char             *structure,
                                 vrna_bp_stack_t  *bp,
                                 int              length),
           "Use vrna_letter_structure() instead");

/**
 *  @brief Create a dot-bracket like structure string from base pair probability matrix
 *  @deprecated Use vrna_db_from_probs() instead!
 *  @ingroup        struct_utils_deprecated
 */
DEPRECATED(void  bppm_to_structure(char         *structure,
                                   FLT_OR_DBL   *pr,
                                   unsigned int length),
           "Use vrna_db_from_probs() instead");

/**
 *  @brief Get a pseudo dot bracket notation for a given probability information
 *  @deprecated Use vrna_bpp_symbol() instead!
 *  @ingroup        struct_utils_deprecated
 */
DEPRECATED(char    bppm_symbol(const float *x),
           "Use vrna_bpp_symbol() instead");

#endif

/**
 * @}
 */

#endif
