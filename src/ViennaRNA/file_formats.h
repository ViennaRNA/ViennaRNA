#ifndef VIENNA_RNA_PACKAGE_FILE_FORMATS_H
#define VIENNA_RNA_PACKAGE_FILE_FORMATS_H

#ifdef VRNA_WARN_DEPRECATED
# ifdef __GNUC__
#  define DEPRECATED(func) func __attribute__ ((deprecated))
# else
#  define DEPRECATED(func) func
# endif
#else
# define DEPRECATED(func) func
#endif

/* make this interface backward compatible with RNAlib < 2.2.0 */
#define VRNA_BACKWARD_COMPAT

/**
 *  @file     file_formats.h
 *  @ingroup  file_utils
 *  @brief    Read and write different file formats for RNA sequences, structures
 */

/**
 *  @{
 *  @ingroup  file_utils
 */

#include <stdio.h>

#include <ViennaRNA/data_structures.h>
#include <ViennaRNA/file_formats_msa.h>

/**
 *  @brief Print a secondary structure as helix list
 *
 *  @param  seq     The RNA sequence
 *  @param  db      The structure in dot-bracket format
 *  @param  energy  Free energy of the structure in kcal/mol
 *  @param  file    The file handle used to print to (print defaults to 'stdout' if(file == NULL) )
 */
void vrna_file_helixlist( const char *seq,
                          const char *db,
                          float energy,
                          FILE *file);

/**
 *  @brief Print a secondary structure as connect table
 *
 *  Connect table file format looks like this:
@verbatim
300  ENERGY = 7.0  example
  1 G       0    2   22    1
  2 G       1    3   21    2
@endverbatim
 *  where the headerline is followed by 6 columns with:
 *  1. Base number: index n
 *  2. Base (A, C, G, T, U, X)
 *  3. Index n-1  (0 if first nucleotide)
 *  4. Index n+1  (0 if last nucleotide)
 *  5. Number of the base to which n is paired. No pairing is indicated by 0 (zero).
 *  6. Natural numbering.
 *
 *  @param  seq         The RNA sequence
 *  @param  db          The structure in dot-bracket format
 *  @param  energy      The free energy of the structure
 *  @param  identifier  An optional identifier for the sequence
 *  @param  file  The file handle used to print to (print defaults to 'stdout' if(file == NULL) )
 */
void vrna_file_connect( const char *seq,
                        const char *db,
                        float energy,
                        const char *identifier,
                        FILE *file);

/**
 *  @brief Print a secondary structure in bpseq format
 *
 *  @param  seq         The RNA sequence
 *  @param  db          The structure in dot-bracket format
 *  @param  file  The file handle used to print to (print defaults to 'stdout' if(file == NULL) )
 */
void vrna_file_bpseq( const char *seq,
                      const char *db,
                      FILE *file);

#if VRNA_WITH_JSON_SUPPORT

/**
 *  @brief Print a secondary structure in jsonformat
 *
 *  @param  seq         The RNA sequence
 *  @param  db          The structure in dot-bracket format
 *  @param  energy      The free energy
 *  @param  identifier  An identifier for the sequence
 *  @param  file  The file handle used to print to (print defaults to 'stdout' if(file == NULL) )
 */
void vrna_file_json(const char *seq,
                    const char *db,
                    double energy,
                    const char *identifier,
                    FILE *file);

#endif

/**
 *  @brief  Tell a function that an input is assumed to span several lines
 *
 *  If used as input-option a function might also be returning this state telling
 *  that it has read data from multiple lines.
 *
 *  @see vrna_extract_record_rest_structure(), vrna_file_fasta_read_record()
 *
 */
#define VRNA_OPTION_MULTILINE             32U
/**
 *  @brief parse multiline constraint
 *  @deprecated see vrna_extract_record_rest_structure()
 */
#define VRNA_CONSTRAINT_MULTILINE         32U

/**
 *  @brief  Get a (fasta) data set from a file or stdin
 * 
 *  This function may be used to obtain complete datasets from a filehandle or stdin.
 *  A dataset is always defined to contain at least a sequence. If data starts with a
 *  fasta header, i.e. a line like
 *  @verbatim >some header info @endverbatim
 *  then vrna_file_fasta_read_record() will assume that the sequence that follows the header may span
 *  over several lines. To disable this behavior and to assign a single line to the argument
 *  'sequence' one can pass #VRNA_INPUT_NO_SPAN in the 'options' argument.
 *  If no fasta header is read in the beginning of a data block, a sequence must not span over
 *  multiple lines!\n
 *  Unless the options #VRNA_INPUT_NOSKIP_COMMENTS or #VRNA_INPUT_NOSKIP_BLANK_LINES are passed,
 *  a sequence may be interrupted by lines starting with a comment character or empty lines.\n
 *  A sequence is regarded as completely read if it was either assumed to not span over multiple
 *  lines, a secondary structure or structure constraint follows the sequence on the next line,
 *  or a new header marks the beginning of a new sequence...\n
 *  All lines following the sequence (this includes comments) that do not initiate a new dataset
 *  according to the above definition are available through the line-array 'rest'.
 *  Here one can usually find the structure constraint or other information belonging to the
 *  current dataset. Filling of 'rest' may be prevented by passing #VRNA_INPUT_NO_REST to the
 *  options argument.\n
 * 
 *  @note This function will exit any program with an error message if no sequence could be read!
 *  @note This function is NOT threadsafe! It uses a global variable to store information about
 *  the next data block.
 * 
 *  The main purpose of this function is to be able to easily parse blocks of data
 *  in the header of a loop where all calculations for the appropriate data is done inside the
 *  loop. The loop may be then left on certain return values, e.g.:
 *  @code
char *id, *seq, **rest;
int  i;
id = seq = NULL;
rest = NULL;
while(!(vrna_file_fasta_read_record(&id, &seq, &rest, NULL, 0) & (VRNA_INPUT_ERROR | VRNA_INPUT_QUIT))){
  if(id) printf("%s\n", id);
  printf("%s\n", seq);
  if(rest)
    for(i=0;rest[i];i++){
      printf("%s\n", rest[i]);
      free(rest[i]);
    }
  free(rest);
  free(seq);
  free(id);
}
 *  @endcode
 *  In the example above, the while loop will be terminated when vrna_file_fasta_read_record() returns
 *  either an error, EOF, or a user initiated quit request.\n
 *  As long as data is read from stdin (we are passing NULL as the file pointer), the id is
 *  printed if it is available for the current block of data. The sequence will be printed in
 *  any case and if some more lines belong to the current block of data each line will be printed
 *  as well.
 * 
 *  @note Do not forget to free the memory occupied by header, sequence and rest!
 * 
 *  @param  header    A pointer which will be set such that it points to the header of the record
 *  @param  sequence  A pointer which will be set such that it points to the sequence of the record
 *  @param  rest      A pointer which will be set such that it points to an array of lines which also belong to the record
 *  @param  file      A file handle to read from (if NULL, this function reads from stdin)
 *  @param  options   Some options which may be passed to alter the behavior of the function, use 0 for no options
 *  @return           A flag with information about what the function actually did read
 */
unsigned int vrna_file_fasta_read_record(char **header,
                                    char **sequence,
                                    char  ***rest,
                                    FILE *file,
                                    unsigned int options);

/** @brief Extract a dot-bracket structure string from (multiline)character array
 *
 * This function extracts a dot-bracket structure string from the 'rest' array as
 * returned by vrna_file_fasta_read_record() and returns it. All occurences of comments within the
 * 'lines' array will be skipped as long as they do not break the structure string.
 * If no structure could be read, this function returns NULL.
 *
 * @pre      The argument 'lines' has to be a 2-dimensional character array as obtained
 *            by vrna_file_fasta_read_record()
 * @see vrna_file_fasta_read_record()
 *
 * @param lines   The (multiline) character array to be parsed
 * @param length  The assumed length of the dot-bracket string (passing a value < 1 results in no length limit)
 * @param option  Some options which may be passed to alter the behavior of the function, use 0 for no options
 * @return        The dot-bracket string read from lines or NULL
 */
char *vrna_extract_record_rest_structure( const char **lines,
                                          unsigned int length,
                                          unsigned int option);

/**
  * @brief Read data from a given SHAPE reactivity input file
  *
  * This function parses the informations from a given file and stores the result
  * in the preallocated string sequence and the double array values.
  *
  * @param file_name     Path to the constraints file
  * @param length        Length of the sequence (file entries exceeding this limit will cause an error)
  * @param default_value Value for missing indices
  * @param sequence      Pointer to an array used for storing the sequence obtained from the SHAPE reactivity file
  * @param values        Pointer to an array used for storing the values obtained from the SHAPE reactivity file
  */
int vrna_file_SHAPE_read( const char *file_name,
                          int length,
                          double default_value,
                          char *sequence,
                          double *values);


#ifdef  VRNA_BACKWARD_COMPAT

/**
 *  @brief  Extract a hard constraint encoded as pseudo dot-bracket string
 *
 *  @deprecated     Use vrna_extract_record_rest_structure() instead!
 *  @pre      The argument 'lines' has to be a 2-dimensional character array as obtained
 *            by vrna_file_fasta_read_record()
 *  @see      vrna_file_fasta_read_record(), #VRNA_CONSTRAINT_DB_PIPE, #VRNA_CONSTRAINT_DB_DOT, #VRNA_CONSTRAINT_DB_X
 *            #VRNA_CONSTRAINT_DB_ANG_BRACK, #VRNA_CONSTRAINT_DB_RND_BRACK
 *
 *  @param  cstruc  A pointer to a character array that is used as pseudo dot-bracket
 *                  output
 *  @param  lines   A 2-dimensional character array with the extension lines from the FASTA
 *                  input
 *  @param  option  The option flags that define the behavior and recognition pattern of
 *                  this function
 */
DEPRECATED(void vrna_extract_record_rest_constraint( char **cstruc, const char **lines, unsigned int option));

/** @brief Extract a dot-bracket structure string from (multiline)character array
 *
 * @deprecated This function is deprecated! Use \fn vrna_extract_record_rest_structure() as a replacment.
 */
DEPRECATED(char *extract_record_rest_structure( const char **lines,
                                                unsigned int length,
                                                unsigned int option));

/**
 *  @brief  Get a data record from stdin
 * 
 *  @deprecated This function is deprecated! Use vrna_file_fasta_read_record() as a replacment.
 *
 */
DEPRECATED(unsigned int read_record(char **header,
                                    char **sequence,
                                    char  ***rest,
                                    unsigned int options));


DEPRECATED(unsigned int get_multi_input_line(char **string, unsigned int options));

#endif

/**
 * @}
 */

#endif
