#ifndef __VIENNA_RNA_PACKAGE_STRUCT_UTILS_H__
#define __VIENNA_RNA_PACKAGE_STRUCT_UTILS_H__

/**
 *  \file structure_utils.h
 *  \brief Various utility- and helper-functions for secondary structure parsing, converting, etc.
 */

#include <stdio.h>

#include <ViennaRNA/data_structures.h>

#ifdef __GNUC__
#define DEPRECATED(func) func __attribute__ ((deprecated))
#else
#define DEPRECATED(func) func
#endif

/**
 *  \brief Pack secondary secondary structure, 5:1 compression using base 3 encoding
 *
 *  Returns a binary string encoding of the secondary structure using
 *  a 5:1 compression scheme. The string is NULL terminated and can
 *  therefore be used with standard string functions such as strcmp().
 *  Useful for programs that need to keep many structures in memory.
 *
 *  \param struc    The secondary structure in dot-bracket notation
 *  \return         The binary encoded structure
 */
char *pack_structure(const char *struc);

/**
 *  \brief Unpack secondary structure previously packed with pack_structure()
 *
 *  Translate a compressed binary string produced by pack_structure() back into
 *  the familiar dot-bracket notation.
 *
 *  \param packed   The binary encoded packed secondary structure
 *  \return         The unpacked secondary structure in dot-bracket notation
 */
char *unpack_structure(const char *packed);

/**
 *  \brief Create a pair table of a secondary structure
 *
 *  Returns a newly allocated table, such that table[i]=j if (i.j) pair
 *  or 0 if i is unpaired, table[0] contains the length of the structure.
 *
 *  \param  structure The secondary structure in dot-bracket notation
 *  \return           A pointer to the created pair_table
 */
short *make_pair_table(const char *structure);

short *make_pair_table_pk(const char *structure);

/**
 *  \brief Get an exact copy of a pair table
 *
 *  \param pt The pair table to be copied
 *  \return   A pointer to the copy of 'pt' 
 */
short *copy_pair_table(const short *pt);

/**
***Pair table for snoop align
***
***
**/
short *alimake_pair_table(const char *structure);

/**
*** returns a newly allocated table, such that:  table[i]=j if (i.j) pair or
*** 0 if i is unpaired, table[0] contains the length of the structure.
*** The special pseudoknotted H/ACA-mRNA structure is taken into account.
**/
short *make_pair_table_snoop(const char *structure);

/**
 *  \brief Compute the "base pair" distance between two secondary structures s1 and s2.
 * 
 *  The sequences should have the same length.
 *  dist = number of base pairs in one structure but not in the other
 *  same as edit distance with open-pair close-pair as move-set
 * 
 *  \param str1   First structure in dot-bracket notation
 *  \param str2   Second structure in dot-bracket notation
 *  \return       The base pair distance between str1 and str2
 */

int *make_loop_index_pt(short *pt);


/**
 *  \brief Convert a pair table into dot-parenthesis notation
 *
 *  \param pt The pair table to be copied
 *  \return   A char pointer to the dot-bracket string
 */
char *vrna_pt_to_db(short *pt);


int bp_distance(const char *str1,
                const char *str2);

unsigned int  *make_referenceBP_array(short *reference_pt,
                                      unsigned int turn);

unsigned int  *compute_BPdifferences( short *pt1,
                                      short *pt2,
                                      unsigned int turn);


#endif
