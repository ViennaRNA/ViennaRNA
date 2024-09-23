#ifndef VIENNA_RNA_PACKAGE_STRUCTURES_PAIRTABLE_H
#define VIENNA_RNA_PACKAGE_STRUCTURES_PAIRTABLE_H

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
 *  @file     ViennaRNA/structures/pairtable.h
 *  @ingroup  struct_utils
 *  @brief    Various utility- and helper-functions for secondary structure parsing, converting, etc.
 */

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
 *  @param  structure Secondary structure in @ref dot-bracket-ext-notation
 *  @param  options   A bitmask to specify which brackets are recognized during conversion to pair table
 *  @return           A pointer to a new pair table of the provided secondary structure
 */
short *
vrna_ptable_from_string(const char    *structure,
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
 *  vrna_ptable_from_string(structure, VRNA_BRACKETS_RND | VRNA_BRACKETS_SQR)
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


#ifndef VRNA_DISABLE_BACKWARD_COMPATIBILITY

/*###########################################*/
/*# deprecated functions below              #*/
/*###########################################*/

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

#endif

#endif
