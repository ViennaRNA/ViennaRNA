#ifndef VIENNA_RNA_PACKAGE_FILE_UTILS_H
#define VIENNA_RNA_PACKAGE_FILE_UTILS_H

/**
 *  @file     file_utils.h
 *  @ingroup  file_utils
 *  @brief    Several utilities for file handling
 */

/**
 *  @addtogroup file_utils
 *  @brief      Functions dealing with file formats for RNA sequences, structures, and alignments
 *
 *  @{
 *  @ingroup  file_utils
 */

/**
 *  @brief Inefficient `cp'
 */
void vrna_file_copy(FILE  *from,
                    FILE  *to);


/**
 *  @brief Read a line of arbitrary length from a stream
 *
 *  Returns a pointer to the resulting string. The necessary memory is
 *  allocated and should be released using @e free() when the string is
 *  no longer needed.
 *
 *  @param  fp  A file pointer to the stream where the function should read from
 *  @return     A pointer to the resulting string
 */
char *vrna_read_line(FILE *fp);


/**
 *  @brief  Recursivly create a directory tree
 */
int vrna_mkdir_p(const char *path);


/**
 *  @brief  Extract the filename from a file path
 */
char *vrna_basename(const char *path);


/**
 *  @brief  Extract the directory part of a file path
 */
char *vrna_dirname(const char *path);


/**
 *  @brief  Sanitize a file name
 *
 *  Returns a new file name where all invalid characters are
 *  substituted by a replacement character. If no replacement
 *  character is supplied, invalid characters are simply removed
 *  from the filename. File names may also never exeed a length
 *  of 255 characters. Longer file names will be truncated!
 *
 *  For now we consider the following characters invalid:
 *  - backslash '\'
 *  - slash '/'
 *  - question mark '?'
 *  - percent sign '%'
 *  - asterisk '*'
 *  - colon ':'
 *  - pipe symbol '|'
 *  - double quote '"'
 *  - triangular brackets '<' and '>'
 *
 *  Furthermore, the (resulting) file name must not be a reserved
 *  file name, such as:
 *  - '.'
 *  - '..'
 *
 *  @note This function may return NULL if the input is pointing
 *        to NULL, or an empty string if the input only consists
 *        of invalid characters which are simply removed!
 *
 *  @param  name        The input file name
 *  @param  replacment  The replacement character, or NULL
 *  @return             The sanitized file name, or NULL
 */
char *vrna_filename_sanitize(const char *name,
                             const char *replacement);


/**
 *  @}
 */

#endif
