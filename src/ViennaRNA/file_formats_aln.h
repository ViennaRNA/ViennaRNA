#ifndef VIENNA_RNA_PACKAGE_FILE_FORMATS_ALN_H
#define VIENNA_RNA_PACKAGE_FILE_FORMATS_ALN_H

/**
 *  @addtogroup   file_utils
 *
 *  @{
 *
 *  @file file_formats_aln.h
 *  @brief Functions dealing with file formats for RNA sequence alignments
 *
 */

#include <stdio.h>

/**
 *  @brief  Option flag for vrna_file_alignment_read() to enable parsing of ClustalW formatted files
 */
#define VRNA_FILE_FORMAT_ALN_CLUSTAL      1U

/**
 *  @brief Option flag for vrna_file_alignment_read() to enable parsing of Stockholm 1.0 formatted files
 */
#define VRNA_FILE_FORMAT_ALN_STOCKHOLM    2U

/**
 *  @brief Option flag for vrna_file_alignment_read() to enable parsing of FASTA (Pearson) formatted files
 */
#define VRNA_FILE_FORMAT_ALN_FASTA        4U

/**
 *  @brief Option flag for vrna_file_alignment_read() to enable parsing of default file formats
 */
#define VRNA_FILE_FORMAT_ALN_DEFAULT      (   VRNA_FILE_FORMAT_ALN_CLUSTAL \
                                            | VRNA_FILE_FORMAT_ALN_STOCKHOLM \
                                            | VRNA_FILE_FORMAT_ALN_FASTA \
                                          )

/**
 *  @brief Option flag for vrna_file_alignment_read() to disable validation of the alignment
 */
#define VRNA_FILE_FORMAT_ALN_NOCHECK      4096U

/**
 *  @brief Read a multiple sequence alignment from file
 *
 *  This function reads the (first) multiple sequence alignment from
 *  an input file. The read alignment is split into the sequence id/name
 *  part and the actual sequence information and stored in memory as
 *  arrays of ids/names and sequences. If the alignment file format
 *  allows for additional information, such as an ID of the entire alignment
 *  or consensus structure information, this data is retrieved as well
 *  and made available. The @p options parameter allows to specify the
 *  set of alignment file formats that should be used to retrieve the data.
 *  If 0 is passed as option, the list of alignment file formats defaults to
 *  #VRNA_FILE_FORMAT_ALN_DEFAULT.
 *
 *  Currently, the list of parsable multiple sequence alignment file formats
 *  consists of:
 *  - @ref msa-formats-clustal
 *  - @ref msa-formats-stockholm
 *  - @ref msa-formats-fasta
 *  .
 *
 *  @note After successfully reading an alignment, this function performs
 *        a validation of the data that includes uniqueness of the sequence
 *        identifiers, and equal sequence lengths. This check can be
 *        deactivated by passing #VRNA_FILE_FORMAT_ALN_NOCHECK in the
 *        @p options parameter.
 *
 *  @see  #VRNA_FILE_FORMAT_ALN_CLUSTAL, #VRNA_FILE_FORMAT_ALN_STOCKHOLM,
 *        #VRNA_FILE_FORMAT_ALN_FASTA, #VRNA_FILE_FORMAT_ALN_DEFAULT,
 *        #VRNA_FILE_FORMAT_ALN_NOCHECK
 *
 *  @param  filename    The name of input file that contains the alignment
 *  @param  names       An address to the pointer where sequence identifiers
 *                      should be written to
 *  @param  aln         An address to the pointer where aligned sequences should
 *                      be written to
 *  @param  id          An address to the pointer where the alignment ID should
 *                      be written to (Maybe NULL)
 *  @param  structure   An address to the pointer where consensus structure
 *                      information should be written to (Maybe NULL)
 *  @param  options     Options to manipulate the behavior of this function
 *  @return             The number of sequences in the alignment
 */
int
vrna_file_alignment_read( const char *filename,
                          char ***names,
                          char ***aln,
                          char  **id,
                          char  **structure,
                          unsigned int options);

int
vrna_file_stockholm_read_record(  FILE *fp,
                                  char ***aln,
                                  char ***names,
                                  char  **id,
                                  char  **structure,
                                  int   verbosity);


/**
 * @}
 */

#endif
