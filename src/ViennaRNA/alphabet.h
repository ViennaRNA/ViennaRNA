#ifndef VIENNA_RNA_PACKAGE_ALPHABET_H
#define VIENNA_RNA_PACKAGE_ALPHABET_H

/* make this interface backward compatible with RNAlib < 2.2.0 */
#define VRNA_BACKWARD_COMPAT

#ifdef VRNA_WARN_DEPRECATED
# ifdef __GNUC__
#  define DEPRECATED(func) func __attribute__ ((deprecated))
# else
#  define DEPRECATED(func) func
# endif
#else
# define DEPRECATED(func) func
#endif

/**
 *  @file     alphabet.h
 *  @ingroup  utils
 *  @brief    Functions to process, convert, and generally handle different nucleotide
 *            and/or base pair alphabets
 */

/**
 *  @{
 *  @ingroup utils
 */

#include <ViennaRNA/model.h>

unsigned int vrna_sequence_length_max(unsigned int options);

int vrna_nucleotide_IUPAC_identity(char a, char b);

/**
 *  @brief Get an array of the numerical encoding for each possible base pair (i,j)
 *
 *  @note This array is always indexed in column-wise order, in contrast to previously
 *  different indexing between mfe and pf variants!
 *
 *  @see  vrna_idx_col_wise(), #vrna_fold_compound_t
 *
 */
char  *vrna_ptypes( const short *S,
                    vrna_md_t *md);

/**
 *  @brief Get a numerical representation of the nucleotide sequence
 *
 */
short *vrna_seq_encode( const char *sequence,
                        vrna_md_t *md);

/**
 *  @brief Get a numerical representation of the nucleotide sequence (simple version)
 *
 */
short *vrna_seq_encode_simple(const char *sequence,
                              vrna_md_t *md);

/**
 *  @brief  Encode a nucleotide character to numerical value
 *
 *  This function encodes a nucleotide character to its numerical representation as required by many functions in RNAlib.
 *
 *  @see  vrna_nucleotide_decode(), vrna_seq_encode()
 *
 *  @param  c   The nucleotide character to encode
 *  @param  md  The model details that determine the kind of encoding
 *  @return     The encoded nucleotide
 */
int vrna_nucleotide_encode( char c,
                            vrna_md_t *md);

/**
 *  @brief  Decode a numerical representation of a nucleotide back into nucleotide alphabet
 *
 *  This function decodes a numerical representation of a nucleotide character back into nucleotide alphabet
 *
 *  @see  vrna_nucleotide_encode(), vrna_seq_encode()
 *
 *  @param  enc The encoded nucleotide
 *  @param  md  The model details that determine the kind of decoding
 *  @return     The decoded nucleotide character
 */
char vrna_nucleotide_decode(int enc,
                            vrna_md_t *md);

void vrna_aln_encode( const char *sequence,
                      short **S_p,
                      short **s5_p,
                      short **s3_p,
                      char **ss_p,
                      unsigned short **as_p,
                      vrna_md_t *md);

/**
 *  @}
 */

#ifdef  VRNA_BACKWARD_COMPAT

DEPRECATED(char  *get_ptypes(const short *S, vrna_md_t *md, unsigned int idx_type));

#endif

#endif
