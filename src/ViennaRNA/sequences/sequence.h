#ifndef VIENNA_RNA_PACKAGE_SEQUENCES_SEQUENCE_H
#define VIENNA_RNA_PACKAGE_SEQUENCES_SEQUENCE_H

/**
 *  @file sequence.h
 *  @brief  Functions and data structures related to sequence representations
 *  @ingroup utils, alphabet_utils
 */

/**
 *  @addtogroup alphabet_utils
 *  @{
 */


/** @brief Typename for nucleotide sequence representation data structure #vrna_sequence_s */
typedef struct vrna_sequence_s vrna_seq_t;

typedef struct vrna_alignment_s vrna_msa_t;

#include <ViennaRNA/fold_compound.h>


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
  char            *name;
  char            *string;    /**< @brief The string representation of the sequence */
  short           *encoding;  /**< @brief The integer representation of the sequence */
  short           *encoding5;
  short           *encoding3;
  unsigned int    length;     /**< @brief The length of the sequence */
};


struct vrna_alignment_s {
  unsigned int  n_seq;
  vrna_seq_t          *sequences;
  char                **gapfree_seq;
  unsigned int        *gapfree_size;  /* for MAF alignment coordinates */
  unsigned long long  *genome_size;     /* for MAF alignment coordinates */
  unsigned long long  *start;           /* for MAF alignment coordinates */
  unsigned char       *orientation;     /* for MAF alignment coordinates */
  unsigned int        **a2s;
};


vrna_seq_t *
vrna_sequence(const char    *string,
              unsigned int  options);


int
vrna_sequence_add(vrna_fold_compound_t  *fc,
                  const char            *string,
                  unsigned int          options);


int
vrna_sequence_remove(vrna_fold_compound_t *fc,
                     unsigned int         i);


void
vrna_sequence_remove_all(vrna_fold_compound_t *fc);


void
vrna_sequence_prepare(vrna_fold_compound_t *fc);


int
vrna_sequence_order_update(vrna_fold_compound_t *fc,
                           const unsigned int   *order);


int
vrna_msa_add( vrna_fold_compound_t      *fc,
              const char                **alignment,
              const char                **names,
              const unsigned char       *orientation,
              const unsigned long long  *start,
              const unsigned long long  *genome_size,
              unsigned int              options);


/**
 *  @}
 */

#endif
