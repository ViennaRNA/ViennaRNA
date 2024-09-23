#ifndef VIENNA_RNA_PACKAGE_SEQUENCES_UTILS_H
#define VIENNA_RNA_PACKAGE_SEQUENCES_UTILS_H

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
 *  @file     ViennaRNA/sequences/utils.h
 *  @ingroup  utils, string_utils
 *  @brief    General utility- and helper-functions for RNA sequences used throughout the ViennaRNA Package
 */

/**
 *  @addtogroup   string_utils
 *  @{
 */

/**
 *  @brief Convert an input sequence (possibly containing DNA alphabet characters) to RNA alphabet
 *
 *  This function substitudes <i>T</i> and <i>t</i> with <i>U</i> and <i>u</i>, respectively
 *
 *  @param sequence The sequence to be converted
 */
void
vrna_seq_toRNA(char *sequence);


/**
 *  @brief  Retrieve a DNA sequence which resembles the complement of the input sequence
 *
 *  This function returns a mew DNA string which is the complement
 *  of the input, i.e. the nucleotide letters `A`,`C`,`G`, and `T`
 *  are substituted by their complements `T`,`G`,`C`, and `A`, respectively.
 *
 *  Any characters not belonging to the alphabet of the 4 canonical
 *  bases of DNA are not altered.
 *
 *  @note This function also handles lower-case input sequences and
 *        treats `U` of the RNA alphabet equally to `T`
 *
 *  @see vrna_seq_reverse()
 *
 *  @param  sequence  the input DNA sequence
 *  @return           The complement of the input DNA sequence
 */
char *
vrna_DNA_complement(const char *sequence);


/**
 *  @brief  Remove gap characters from a nucleotide sequence
 *
 *  @param  sequence  The original, null-terminated nucleotide sequence
 *  @return           A copy of the input sequence with all gap characters removed
 */
char *
vrna_seq_ungapped(const char *sequence);


/**
 *  @}
 */

#ifndef VRNA_DISABLE_BACKWARD_COMPATIBILITY

/**
 *  @brief Convert a DNA input sequence to RNA alphabet
 *
 *  @deprecated Use vrna_seq_toRNA() instead!
 */
DEPRECATED(void
           str_DNA2RNA(char *sequence),
           "Use vrna_seq_toRNA() instead");

#endif

#endif
