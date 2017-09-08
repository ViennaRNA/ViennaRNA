#ifndef VIENNA_RNA_PACKAGE_SEQUENCE_H
#define VIENNA_RNA_PACKAGE_SEQUENCE_H

/**
 *  @{
 *
 *  @file sequence.h
 *  @brief  Functions and data structures related to sequence representations
 *
 *  @ingroup utils
 */


/** @brief Typename for nucleotide sequence representation data structure #vrna_sequence_s */
typedef struct vrna_sequence_s vrna_seq_t;


#define VRNA_SEQUENCE_RNA       1U

#define VRNA_SEQUENCE_DNA       2U

/**
 *  @brief  A enumerator used in #vrna_sequence_s to distinguish different nucleotide sequences
 */
typedef enum {
  VRNA_SEQ_UNKNOWN,   /**< @brief Nucleotide sequence represents an Unkown type */
  VRNA_SEQ_RNA,       /**< @brief Nucleotide sequence represents an RNA type */
  VRNA_SEQ_DNA        /**< @brief Nucleotide sequence represents a DNA type */
} vrna_seq_type_e;


/**
 *  @brief  Data structure representing a nucleotide sequence
 */
struct vrna_sequence_s {
  vrna_seq_type_e type;       /**< @brief The type of sequence */
  char            *string;    /**< @brief The string representation of the sequence */
  short           *encoding;  /**< @brief The integer representation of the sequence */
  unsigned int    length;     /**< @brief The length of the sequence */
};


vrna_seq_t *vrna_sequence(const char    *string,
                          unsigned int  options);


int           vrna_sequence_add(vrna_fold_compound_t  *vc,
                                const char            *string,
                                unsigned int          options);


int           vrna_sequence_remove(vrna_fold_compound_t *vc,
                                   unsigned int         i);


void          vrna_sequence_remove_all(vrna_fold_compound_t *vc);


void          vrna_sequence_prepare(vrna_fold_compound_t *fc);


/**
 *  @}
 */

#endif
