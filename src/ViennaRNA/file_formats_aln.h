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

#define VRNA_FILE_FORMAT_ALN_CLUSTAL      1U
#define VRNA_FILE_FORMAT_ALN_STOCKHOLM    2U
#define VRNA_FILE_FORMAT_ALN_FASTA        4U
#define VRNA_FILE_FORMAT_ALN_DEFAULT      (   VRNA_FILE_FORMAT_ALN_CLUSTAL \
                                            | VRNA_FILE_FORMAT_ALN_STOCKHOLM \
                                            | VRNA_FILE_FORMAT_ALN_FASTA \
                                          )

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
