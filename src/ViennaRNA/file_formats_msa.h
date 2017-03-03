#ifndef VIENNA_RNA_PACKAGE_FILE_FORMATS_MSA_H
#define VIENNA_RNA_PACKAGE_FILE_FORMATS_MSA_H

/**
 *  @addtogroup   file_utils
 *
 *  @{
 *
 *  @file file_formats_msa.h
 *  @brief Functions dealing with file formats for Multiple Sequence Alignments (MSA)
 *
 */

#include <stdio.h>

/**
 *  @brief  Option flag indicating ClustalW formatted files
 *  @see vrna_file_msa_read(), vrna_file_msa_read_record(), vrna_file_msa_detect_format()
 */
#define VRNA_FILE_FORMAT_MSA_CLUSTAL      1U

/**
 *  @brief Option flag indicating Stockholm 1.0 formatted files
 *  @see vrna_file_msa_read(), vrna_file_msa_read_record(), vrna_file_msa_detect_format()
 */
#define VRNA_FILE_FORMAT_MSA_STOCKHOLM    2U

/**
 *  @brief Option flag indicating FASTA (Pearson) formatted files
 *  @see vrna_file_msa_read(), vrna_file_msa_read_record(), vrna_file_msa_detect_format()
 */
#define VRNA_FILE_FORMAT_MSA_FASTA        4U

/**
 *  @brief Option flag indicating MAF formatted files
 *  @see vrna_file_msa_read(), vrna_file_msa_read_record(), vrna_file_msa_detect_format()
 */
#define VRNA_FILE_FORMAT_MSA_MAF          8U

/**
 *  @brief Option flag indicating most informative sequence (MIS) output
 *
 *  The default reference sequence output for an alignment is simply a consensus sequence.
 *  This flag allows to write the most informative equence (MIS) instead.
 *
 *  @see vrna_file_msa_write()
 */
#define VRNA_FILE_FORMAT_MSA_MIS          16U

/**
 *  @brief Option flag indicating the set of default file formats
 *  @see vrna_file_msa_read(), vrna_file_msa_read_record(), vrna_file_msa_detect_format()
 */
#define VRNA_FILE_FORMAT_MSA_DEFAULT      ( \
    VRNA_FILE_FORMAT_MSA_CLUSTAL \
    | VRNA_FILE_FORMAT_MSA_STOCKHOLM \
    | VRNA_FILE_FORMAT_MSA_FASTA \
    | VRNA_FILE_FORMAT_MSA_MAF \
    )

/**
 *  @brief Option flag to disable validation of the alignment
 *  @see  vrna_file_msa_read(), vrna_file_msa_read_record()
 */
#define VRNA_FILE_FORMAT_MSA_NOCHECK      4096U

/**
 *  @brief Return flag of vrna_file_msa_detect_format() to indicate unknown or malformatted alignment
 *  @see vrna_file_msa_detect_format()
 */
#define VRNA_FILE_FORMAT_MSA_UNKNOWN      8192U

#define VRNA_FILE_FORMAT_MSA_APPEND       16384U

#define VRNA_FILE_FORMAT_MSA_QUIET        32768U

#define VRNA_FILE_FORMAT_MSA_SILENT       65536U

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
 *  #VRNA_FILE_FORMAT_MSA_DEFAULT.
 *
 *  Currently, the list of parsable multiple sequence alignment file formats
 *  consists of:
 *  - @ref msa-formats-clustal
 *  - @ref msa-formats-stockholm
 *  - @ref msa-formats-fasta
 *  - @ref msa-formats-maf
 *  .
 *
 *  @note After successfully reading an alignment, this function performs
 *        a validation of the data that includes uniqueness of the sequence
 *        identifiers, and equal sequence lengths. This check can be
 *        deactivated by passing #VRNA_FILE_FORMAT_MSA_NOCHECK in the
 *        @p options parameter.
 *
 *  @see  vrna_file_msa_read_record(), #VRNA_FILE_FORMAT_MSA_CLUSTAL,
 *        #VRNA_FILE_FORMAT_MSA_STOCKHOLM, #VRNA_FILE_FORMAT_MSA_FASTA,
 *        #VRNA_FILE_FORMAT_MSA_MAF, #VRNA_FILE_FORMAT_MSA_DEFAULT,
 *        #VRNA_FILE_FORMAT_MSA_NOCHECK
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
 *  @return             The number of sequences in the alignment, or -1 if
 *                      no alignment record could be found
 */
int
vrna_file_msa_read(const char   *filename,
                   char         ***names,
                   char         ***aln,
                   char         **id,
                   char         **structure,
                   unsigned int options);


/**
 *  @brief Read a multiple sequence alignment from file handle
 *
 *  Similar to vrna_file_msa_read(), this function reads a multiple
 *  sequence alignment from an input file handle. Since using a file
 *  handle, this function is not limited to the first alignment record,
 *  but allows for looping over all alignments within the input.
 *
 *  The read alignment is split into the sequence id/name
 *  part and the actual sequence information and stored in memory as
 *  arrays of ids/names and sequences. If the alignment file format
 *  allows for additional information, such as an ID of the entire alignment
 *  or consensus structure information, this data is retrieved as well
 *  and made available. The @p options parameter allows to specify the
 *  alignment file format used to retrieve the data. A single format
 *  must be specified here, see vrna_file_msa_detect_format() for helping
 *  to determine the correct MSA file format.
 *
 *  Currently, the list of parsable multiple sequence alignment file formats
 *  consists of:
 *  - @ref msa-formats-clustal
 *  - @ref msa-formats-stockholm
 *  - @ref msa-formats-fasta
 *  - @ref msa-formats-maf
 *  .
 *
 *  @note After successfully reading an alignment, this function performs
 *        a validation of the data that includes uniqueness of the sequence
 *        identifiers, and equal sequence lengths. This check can be
 *        deactivated by passing #VRNA_FILE_FORMAT_MSA_NOCHECK in the
 *        @p options parameter.
 *
 *  @see  vrna_file_msa_read(), vrna_file_msa_detect_format(),
 *        #VRNA_FILE_FORMAT_MSA_CLUSTAL, #VRNA_FILE_FORMAT_MSA_STOCKHOLM,
 *        #VRNA_FILE_FORMAT_MSA_FASTA, #VRNA_FILE_FORMAT_MSA_MAF,
 *        #VRNA_FILE_FORMAT_MSA_DEFAULT, #VRNA_FILE_FORMAT_MSA_NOCHECK
 *
 *  @param  fp          The file pointer the data will be retrieved from
 *  @param  names       An address to the pointer where sequence identifiers
 *                      should be written to
 *  @param  aln         An address to the pointer where aligned sequences should
 *                      be written to
 *  @param  id          An address to the pointer where the alignment ID should
 *                      be written to (Maybe NULL)
 *  @param  structure   An address to the pointer where consensus structure
 *                      information should be written to (Maybe NULL)
 *  @param  options     Options to manipulate the behavior of this function
 *  @return             The number of sequences in the alignment, or -1 if
 *                      no alignment record could be found
 */
int
vrna_file_msa_read_record(FILE          *fp,
                          char          ***names,
                          char          ***aln,
                          char          **id,
                          char          **structure,
                          unsigned int  options);


/**
 *  @brief Detect the format of a multiple sequence alignment file
 *
 *  This function attempts to determine the format of a file that
 *  supposedly contains a multiple sequence alignment (MSA). This is
 *  useful in cases where a MSA file contains more than a single record
 *  and therefore vrna_file_msa_read() can not be applied, since
 *  it only retrieves the first.
 *  Here, one can try to guess the correct file format using this
 *  function and then loop over the file, record by record using one
 *  of the low-level record retrieval functions for the corresponding
 *  MSA file format.
 *
 *  @note This function parses the entire first record within the
 *        specified file. As a result, it returns #VRNA_FILE_FORMAT_MSA_UNKNOWN
 *        not only if it can't detect the file's format, but also
 *        in cases where the file doesn't contain sequences!
 *
 *  @see  vrna_file_msa_read(), vrna_file_stockholm_read_record(),
 *        vrna_file_clustal_read_record(), vrna_file_fasta_read_record()
 *
 *  @param  filename  The name of input file that contains the alignment
 *  @param  options   Options to manipulate the behavior of this function
 *  @return           The MSA file format, or #VRNA_FILE_FORMAT_MSA_UNKNOWN
 */
unsigned int
vrna_file_msa_detect_format(const char    *filename,
                            unsigned int  options);


/**
 *  @brief Write multiple sequence alignment file
 *
 *  @note Currently, we only support Stockholm 1.0 formatted output
 *
 *  @param  filename  The output filename
 *  @param  names     The array of sequence names / identifies
 *  @param  aln       The array of aligned sequences
 *  @param  id        An optional ID for the alignment
 *  @param  structure An optional consensus structure
 *  @param  source    A string describing the source of the alignment
 *  @param  options   Options to manipulate the behavior of this function
 *  @return           Non-null upon successful written alignment
 */
int
vrna_file_msa_write(const char    *filename,
                    const char    **names,
                    const char    **aln,
                    const char    *id,
                    const char    *structure,
                    const char    *source,
                    unsigned int  options);


/**
 * @}
 */

#endif
