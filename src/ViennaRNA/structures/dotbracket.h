#ifndef VIENNA_RNA_PACKAGE_STRUCTURES_DOTBRACKET_H
#define VIENNA_RNA_PACKAGE_STRUCTURES_DOTBRACKET_H

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
 *  @file     ViennaRNA/structures/dotbracket.h
 *  @ingroup  struct_utils
 *  @brief    Various functions for dot-bracket representation of secondary structures
 */

#include <ViennaRNA/structures/problist.h>

/**
 *  @addtogroup struct_utils_dot_bracket
 *  @{
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


#define   VRNA_GQUAD_DB_SYMBOL      '+'
#define   VRNA_GQUAD_DB_SYMBOL_END  '~'

/**
 *  @brief Pack secondary secondary structure, 5:1 compression using base 3 encoding
 *
 *  Returns a binary string encoding of the secondary structure using
 *  a 5:1 compression scheme. The string is NULL terminated and can
 *  therefore be used with standard string functions such as strcmp().
 *  Useful for programs that need to keep many structures in memory.
 *
 *  @see  vrna_db_unpack()
 *
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
 *
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
 *  This function also converts pair table formatted structures that contain
 *  pseudoknots. Non-nested base pairs result in additional pairs of
 *  parenthesis and brackets within the resulting dot-parenthesis string.
 *  The following pairs are awailable: (), []. {}. <>, as well as pairs of
 *  matching upper-/lower-case characters from the alphabet A-Z.
 *
 *  @note In cases where the level of non-nested base pairs exceeds the
 *        maximum number of 30 different base pair indicators (4 parenthesis/brackets,
 *        26 matching characters), a warning is printed and the remaining base pairs
 *        are left out from the conversion.
 *
 *  @param pt The pair table to be copied
 *  @return   A char pointer to the dot-bracket string
 */
char *
vrna_db_from_ptable(const short *pt);


/**
 *  @brief  Convert a list of base pairs into dot-bracket notation
 *
 *  @see vrna_plist()
 *
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
/** @} */


/**
 *  @addtogroup struct_utils_wuss
 *  @{
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
/** @} */


/**
 *  @addtogroup   struct_utils
 *  @{
 */

/**
 *  @brief Get a loop index representation of a structure
 */
int *
vrna_loopidx_from_ptable(const short *pt);

/** @} */


/**
 *  @addtogroup   struct_utils
 *  @{
 */

/**
 *  @brief Create a dot-bracket like structure string from base pair probability matrix
 */
char *
vrna_db_from_probs(const FLT_OR_DBL *pr,
                   unsigned int     length);

char *
vrna_pairing_tendency(vrna_fold_compound_t *fc);


/**
 *  @brief Get a pseudo dot bracket notation for a given probability information
 */
char
vrna_bpp_symbol(const float *x);


/**
 *  @brief Create a dot-backet/parenthesis structure from backtracking stack
 *
 *  This function is capable to create dot-bracket structures from base pairs
 *  stored within the base pair stack @p bp_stack.
 *
 *  @param bp_stack Base pair stack containing the traced base pairs
 *  @param length   The length of the structure
 *  @return         The secondary structure in dot-bracket notation as
 *                  provided in the input
 */
char *
vrna_db_from_bps(vrna_bps_t   bp_stack,
                 unsigned int length);

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


/** @} */

/**
 *  @addtogroup gquad_parse
 *  @{
 */
/**
 *  @brief  Parse a G-Quadruplex from a dot-bracket structure string
 *
 *  Given a dot-bracket structure (possibly) containing gquads encoded
 *  by '+' signs (and an optional '~' end sign, find first gquad, return
 *  end position (1-based) or 0 if none found.
 *  Upon return L and l[] contain the number of stacked layers, as well as
 *  the lengths of the linker regions.
 *
 *  @note   For circular RNAs and G-Quadruplexes spanning the n,1-junction
 *          the sum of linkers and g-runs is lower than the end position.
 *          This condition can be used to check whether or not to accept
 *          a G-Quadruplex parsed from the dot-bracket string. Also note,
 *          that such n,1-junction spanning G-Quadruplexes must end with
 *          a `~` sign, to be unambigous.
 *
 *
 *  To parse a string with many gquads, call vrna_gq_parse() repeatedly e.g.
 *
 *  @code
 *  end1 = vrna_gq_parse(struc, &L, l); ... ;
 *  end2 = vrna_gq_parse(struc+end1, &L, l); end2+=end1; ... ;
 *  end3 = vrna_gq_parse(struc+end2, &L, l); end3+=end2; ... ;
 *  @endcode
 *
 *  @param  db_string   The input structure in dot-bracket notation
 *  @param  L           A pointer to an unsigned integer to store the layer (stack) size
 *  @param  l           An array of three values to store the respective linker lenghts
 *  @return             The end position of the G-Quadruplex (1-based) or 0 if not found
 */
unsigned int
vrna_gq_parse(const char *db_string,
              unsigned int *L,
              unsigned int l[3]);

void
vrna_db_insert_gq(char          *db,
                  unsigned int  i,
                  unsigned int  L,
                  unsigned int l[3],
                  unsigned int  n);

/** @} */

#ifndef VRNA_DISABLE_BACKWARD_COMPATIBILITY

/*###########################################*/
/*# deprecated functions below              #*/
/*###########################################*/

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

/**
 *  @addtogroup gquad_deprecated
 *  @{
 */

/**
 *  @brief  Parse a G-Quadruplex from a dot-bracket structure string
 *
 *  Given a dot-bracket structure (possibly) containing gquads encoded
 *  by '+' signs, find first gquad, return end position or 0 if none found
 *  Upon return L and l[] contain the number of stacked layers, as well as
 *  the lengths of the linker regions.
 *  To parse a string with many gquads, call parse_gquad repeatedly e.g.
 *  end1 = parse_gquad(struc, &L, l); ... ;
 *  end2 = parse_gquad(struc+end1, &L, l); end2+=end1; ... ;
 *  end3 = parse_gquad(struc+end2, &L, l); end3+=end2; ... ;
 */
DEPRECATED(int
           parse_gquad(const char *struc,
                       int        *L,
                       int        l[3]),
           "Use vrna_gq_parse() instead");


/**
 * @}
 */

#endif

#endif
