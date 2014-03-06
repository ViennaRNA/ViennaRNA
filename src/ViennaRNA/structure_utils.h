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
char *vrna_pack_structure(const char *struc);

/**
 *  \brief Pack secondary secondary structure, 5:1 compression using base 3 encoding
 *
 *  Returns a binary string encoding of the secondary structure using
 *  a 5:1 compression scheme. The string is NULL terminated and can
 *  therefore be used with standard string functions such as strcmp().
 *  Useful for programs that need to keep many structures in memory.
 *
 *  \deprecated     Use vrna_pack_structure() as a replacement
 *  \param struc    The secondary structure in dot-bracket notation
 *  \return         The binary encoded structure
 */
DEPRECATED(char *pack_structure(const char *struc));

/**
 *  \brief Unpack secondary structure previously packed with pack_structure()
 *
 *  Translate a compressed binary string produced by pack_structure() back into
 *  the familiar dot-bracket notation.
 *
 *  \param packed   The binary encoded packed secondary structure
 *  \return         The unpacked secondary structure in dot-bracket notation
 */
char *vrna_unpack_structure(const char *packed);

/**
 *  \brief Unpack secondary structure previously packed with pack_structure()
 *
 *  Translate a compressed binary string produced by pack_structure() back into
 *  the familiar dot-bracket notation.
 *
 *  \deprecated     Use vrna_unpack_structure() as a replacement
 *  \param packed   The binary encoded packed secondary structure
 *  \return         The unpacked secondary structure in dot-bracket notation
 */
DEPRECATED(char *unpack_structure(const char *packed));

/**
 *  \brief Create a pair table of a secondary structure
 *
 *  Returns a newly allocated table, such that table[i]=j if (i.j) pair
 *  or 0 if i is unpaired, table[0] contains the length of the structure.
 *
 *  \param  structure The secondary structure in dot-bracket notation
 *  \return           A pointer to the created pair_table
 */
short *vrna_pt_get(const char *structure);

/**
 *  \brief Create a pair table of a secondary structure
 *
 *  Returns a newly allocated table, such that table[i]=j if (i.j) pair
 *  or 0 if i is unpaired, table[0] contains the length of the structure.
 *
 *  \deprecated Use vrna_pt_get() instead
 *
 *  \param  structure The secondary structure in dot-bracket notation
 *  \return           A pointer to the created pair_table
 */
DEPRECATED(short *make_pair_table(const char *structure));


/**
 *  \brief Create a pair table of a secondary structure (pseudo-knot version)
 *
 *  Returns a newly allocated table, such that table[i]=j if (i.j) pair
 *  or 0 if i is unpaired, table[0] contains the length of the structure.
 *
 *  In contrast to vrna_pt_get() this function also recognizes the base pairs
 *  denoted by '[' and ']' brackets.
 *
 *  \param  structure The secondary structure in (extended) dot-bracket notation
 *  \return           A pointer to the created pair_table
 */
short *vrna_pt_pk_get(const char *structure);

DEPRECATED(short *make_pair_table_pk(const char *structure));

/**
 *  \brief Get an exact copy of a pair table
 *
 *  \param pt The pair table to be copied
 *  \return   A pointer to the copy of 'pt' 
 */
short *vrna_pt_copy(const short *pt);

/**
 *  \brief Get an exact copy of a pair table
 *
 *  \deprecated Use vrna_pt_copy() instead
 *
 *  \param pt The pair table to be copied
 *  \return   A pointer to the copy of 'pt' 
 */
DEPRECATED(short *copy_pair_table(const short *pt));

/**
 * \brief Create a pair table of a secondary structure (snoop align version)
 *
 */
short *vrna_pt_ali_get(const char *structure);

/**
*** Pair table for snoop align
***
*** \deprecated Use vrna_pt_ali_get() instead!
**/
DEPRECATED(short *alimake_pair_table(const char *structure));

/**
 * \brief Create a pair table of a secondary structure (snoop version)
 *
 *  returns a newly allocated table, such that:  table[i]=j if (i.j) pair or
 *  0 if i is unpaired, table[0] contains the length of the structure.
 *  The special pseudoknotted H/ACA-mRNA structure is taken into account.
 */
short *vrna_pt_snoop_get(const char *structure);

/**
*** returns a newly allocated table, such that:  table[i]=j if (i.j) pair or
*** 0 if i is unpaired, table[0] contains the length of the structure.
*** The special pseudoknotted H/ACA-mRNA structure is taken into account.
*** \deprecated Use vrna_pt_snoop_get() instead!
**/
DEPRECATED(short *make_pair_table_snoop(const char *structure));


/**
 *  \brief Get a loop index representation of a structure
 */
int *vrna_get_loop_index(const short *pt);

DEPRECATED(int *make_loop_index_pt(short *pt));

/**
 *  \brief Convert a pair table into dot-parenthesis notation
 *
 *  \param pt The pair table to be copied
 *  \return   A char pointer to the dot-bracket string
 */
char *vrna_pt_to_db(short *pt);


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
int vrna_bp_distance( const char *str1,
                      const char *str2);

/**
 *  \brief Compute the "base pair" distance between two secondary structures s1 and s2.
 * 
 *  The sequences should have the same length.
 *  dist = number of base pairs in one structure but not in the other
 *  same as edit distance with open-pair close-pair as move-set
 *
 *  \deprecated   Use vrna_bp_distance instead
 *  \param str1   First structure in dot-bracket notation
 *  \param str2   Second structure in dot-bracket notation
 *  \return       The base pair distance between str1 and str2
 */
DEPRECATED(int bp_distance(const char *str1, const char *str2));


/**
 *  \brief Make a reference base pair count matrix
 *
 *  Get an upper triangular matrix containing the number of basepairs of a reference
 *  structure for each interval [i,j] with i<j. Access it via iindx!!!
 */
unsigned int  *vrna_refBPcnt_matrix(const short *reference_pt,
                                    unsigned int turn);

/**
 *  \brief Make a reference base pair count matrix
 *
 *  Get an upper triangular matrix containing the number of basepairs of a reference
 *  structure for each interval [i,j] with i<j. Access it via iindx!!!
 *
 *  \deprecated Use vrna_refBPcnt_matrix() instead
 */
DEPRECATED(unsigned int  *make_referenceBP_array(short *reference_pt,
                                      unsigned int turn));


/**
 *  \brief Make a reference base pair distance matrix
 *
 *  Get an upper triangular matrix containing the base pair distance of two
 *  reference structures for each interval [i,j] with i<j. Access it via iindx!!!
 *
 */
unsigned int  *vrna_refBPdist_matrix( const short *pt1,
                                      const short *pt2,
                                      unsigned int turn);

/**
 *  \brief Make a reference base pair distance matrix
 *
 *  Get an upper triangular matrix containing the base pair distance of two
 *  reference structures for each interval [i,j] with i<j. Access it via iindx!!!
 *
 *  \deprecated Use vrna_refBPdist_matrix() instead
 */
DEPRECATED(unsigned int  *compute_BPdifferences( short *pt1,
                                      short *pt2,
                                      unsigned int turn));

/**
 *  \brief Create a plist from base pair probability matrix
 * 
 *  The probability matrix provided via the #vrna_fold_compound is parsed
 *  and all pair probabilities above the given threshold are used to create
 *  an entry in the plist
 * 
 *  The end of the plist is marked by sequence positions i as well as j
 *  equal to 0. This condition should be used to stop looping over its
 *  entries
 * 
 *  \ingroup            pf_fold
 *  \param[in]  vc      The fold compound
 *  \param[in]  cutoff  The cutoff value
 *  \return             A pointer to the plist that is to be created
 */
plist *vrna_get_plist_from_pr(vrna_fold_compound *vc, double cut_off);

/**
 *  \brief Create a plist from a probability matrix
 * 
 *  The probability matrix given is parsed and all pair probabilities above
 *  the given threshold are used to create an entry in the plist
 * 
 *  The end of the plist is marked by sequence positions i as well as j
 *  equal to 0. This condition should be used to stop looping over its
 *  entries
 * 
 *  \note This function is threadsafe
 *  \deprecated Use vrna_get_plist_from_pr() instead!
 *
 *  \ingroup            pf_fold
 *  \param[out] pl      A pointer to the plist that is to be created
 *  \param[in]  probs   The probability matrix used for creting the plist
 *  \param[in]  length  The length of the RNA sequence
 *  \param[in]  cutoff  The cutoff value
 */
DEPRECATED(void  assign_plist_from_pr( plist **pl,
                            FLT_OR_DBL *probs,
                            int length,
                            double cutoff));

/**
 *  \brief Create a dot-bracket like structure string from base pair probability matrix
 */
void  bppm_to_structure(char *structure,
                        FLT_OR_DBL *pr,
                        unsigned int length);

/**
 *  \brief Get a pseudo dot bracket notation for a given probability information
 */
char    bppm_symbol(const float *x);

/**
 *  \brief Create a dot-backet/parenthesis structure from backtracking stack
 * 
 */
void
vrna_parenthesis_structure( char *structure,
                            bondT *bp,
                            int length);

/**
 *  \brief Create a dot-backet/parenthesis structure from backtracking stack
 *
 *  \deprecated use vrna_parenthesis_structure() instead
 * 
 *  \note This function is threadsafe
 */
DEPRECATED(void 
parenthesis_structure(char *structure,
                      bondT *bp,
                      int length));

/**
 *  \brief Create a dot-backet/parenthesis structure from backtracking stack
 *  obtained by zuker suboptimal calculation in cofold.c
 * 
 *  \note This function is threadsafe
 */
void
vrna_parenthesis_zuker( char *structure,
                        bondT *bp,
                        int length);

/**
 *  \brief Create a dot-backet/parenthesis structure from backtracking stack
 *  obtained by zuker suboptimal calculation in cofold.c
 * 
 *  \deprecated use vrna_parenthesis_zuker instead
 * 
 *  \note This function is threadsafe
 */
DEPRECATED(void 
parenthesis_zuker(char *structure,
                  bondT *bp,
                  int length));


void
vrna_letter_structure(char *structure,
                      bondT *bp,
                      int length);

DEPRECATED(void 
letter_structure( char *structure,
                  bondT *bp,
                  int length));

/**
 *  \brief Create a plist from a dot-bracket string
 * 
 *  The dot-bracket string is parsed and for each base pair an
 *  entry in the plist is created. The probability of each pair in
 *  the list is set by a function parameter.
 * 
 *  The end of the plist is marked by sequence positions i as well as j
 *  equal to 0. This condition should be used to stop looping over its
 *  entries
 * 
 *  \param struc  The secondary structure in dot-bracket notation
 *  \param pr     The probability for each base pair used in the plist
 *  \return       The plist array
 */
plist *vrna_get_plist_from_db(const char *struc, float pr);

/**
 *  \brief Create a plist from a dot-bracket string
 * 
 *  The dot-bracket string is parsed and for each base pair an
 *  entry in the plist is created. The probability of each pair in
 *  the list is set by a function parameter.
 * 
 *  The end of the plist is marked by sequence positions i as well as j
 *  equal to 0. This condition should be used to stop looping over its
 *  entries
 * 
 *  \deprecated   Use vrna_get_plist_from_db() instead
 * 
 *  \param pl     A pointer to the plist that is to be created
 *  \param struc  The secondary structure in dot-bracket notation
 *  \param pr     The probability for each base pair
 */
DEPRECATED(void assign_plist_from_db(plist **pl, const char *struc, float pr));


#endif
