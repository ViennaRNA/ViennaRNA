#ifndef VIENNA_RNA_PACKAGE_SEQUENCES_ALPHABET_H
#define VIENNA_RNA_PACKAGE_SEQUENCES_ALPHABET_H

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
 *  @file     alphabet.h
 *  @ingroup  utils, alphabet_utils
 *  @brief    Functions to process, convert, and generally handle different nucleotide
 *            and/or base pair alphabets
 */

/**
 *  @addtogroup alphabet_utils
 *  @{
 *  @brief  Functions to cope with various aspects related to the nucleotide sequence alphabet
 */

#include <ViennaRNA/fold_compound.h>
#include <ViennaRNA/model.h>

unsigned int
vrna_sequence_length_max(unsigned int options);


int
vrna_nucleotide_IUPAC_identity(char a,
                               char b);


void
vrna_ptypes_prepare(vrna_fold_compound_t  *fc,
                    unsigned int          options);


/**
 *  @brief Get an array of the numerical encoding for each possible base pair (i,j)
 *
 *  @note This array is always indexed in column-wise order, in contrast to previously
 *        different indexing between mfe and pf variants!
 *
 *  @see  vrna_idx_col_wise(), #vrna_fold_compound_t
 *
 */
char *
vrna_ptypes(const short *S,
            vrna_md_t   *md);


/**
 *  @brief Get a numerical representation of the nucleotide sequence
 *
 *  @param  sequence    The input sequence in upper-case letters
 *  @param  md          A pointer to a #vrna_md_t data structure that specifies the conversion type
 *  @return             A list of integer encodings for each sequence letter (1-based). Position 0 denotes the length of the list
 */
short *
vrna_seq_encode(const char  *sequence,
                vrna_md_t   *md);


/**
 *  @brief Get a numerical representation of the nucleotide sequence (simple version)
 *
 */
short *
vrna_seq_encode_simple(const char *sequence,
                       vrna_md_t  *md);


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
int
vrna_nucleotide_encode(char       c,
                       vrna_md_t  *md);


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
char
vrna_nucleotide_decode(int        enc,
                       vrna_md_t  *md);


void
vrna_aln_encode(const char    *sequence,
                short         **S_p,
                short         **s5_p,
                short         **s3_p,
                char          **ss_p,
                unsigned int  **as_p,
                vrna_md_t     *md);


unsigned int
vrna_get_ptype_md(int       i,
                  int       j,
                  vrna_md_t *md);


unsigned int
vrna_get_ptype(int  ij,
               char *ptype);


unsigned int
vrna_get_ptype_window(int   i,
                      int   j,
                      char  **ptype);


/**
 *  @}
 */

#ifndef VRNA_DISABLE_BACKWARD_COMPATIBILITY

DEPRECATED(char *get_ptypes(const short   *S,
                            vrna_md_t     *md,
                            unsigned int  idx_type),
           "Use vrna_pytpes() instead");

#endif

#endif
