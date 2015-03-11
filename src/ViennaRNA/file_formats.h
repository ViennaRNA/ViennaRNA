#ifndef VIENNA_RNA_PACKAGE_FILE_FORMATS_H
#define VIENNA_RNA_PACKAGE_FILE_FORMATS_H

#ifdef __GNUC__
#define DEPRECATED(func) func __attribute__ ((deprecated))
#else
#define DEPRECATED(func) func
#endif

/* make this interface backward compatible with RNAlib < 2.2.0 */
#define VRNA_BACKWARD_COMPAT

/**
 *  @file file_formats.h
 *  @brief Various functions dealing with file formats for RNA sequences, structures, and alignments
 *
 *  @section  file-formats  Input/Output file formats used in RNAlib
 *
 *  @subsection     constraint-formats  File formats for secondary structure constraints
 *  @subsubsection  constraint-formats-file Constraints Definition File
 *
 *  The RNAlib can parse and apply data from constraint definition text files, where each constraint
 *  is given as a line of whitespace delimited commands. The syntax we use extends the one used in
 *  <a href="http://mfold.rna.albany.edu/?q=mfold">mfold</a> /
 *  <a href="http://mfold.rna.albany.edu/?q=DINAMelt/software">UNAfold</a> where
 *  each line begins with a command character followed by a set of positions. Additionally, we allow
 *  an optional loop type context specifier in form of a sequence of characters, and an orientation
 *  flag enables to force a nucleotide to pair upstream, or downstream
 *
 *  @paragraph  const_file_commands Constraint commands
 *  The following set of commands is recognized:
 *  -#  @f$ F \ldots @f$ Force
 *  -#  @f$ P \ldots @f$ Prohibit
 *  -#  @f$ W \ldots @f$ Weakly enforce, i.e. remove conflicts only
 *  -#  @f$ U \ldots @f$ Soft constraint for unpaired position(s)
 *  -#  @f$ B \ldots @f$ Soft constraint for base pair(s)
 *
 *  @paragraph  const_file_loop_types Specification of the loop type context
 *  The optional loop type context specifier @f$ [WHERE] @f$ may be a combination of the following:
 *  -#  @f$ E \ldots @f$ Exterior loop
 *  -#  @f$ H \ldots @f$ Hairpin loop
 *  -#  @f$ I \ldots @f$ Interior loop (enclosing pair)
 *  -#  @f$ i \ldots @f$ Interior loop (enclosed pair)
 *  -#  @f$ M \ldots @f$ Multibranch loop (enclosing pair)
 *  -#  @f$ m \ldots @f$ Multibranch loop (enclosed pair)
 *
 *  If no @f$ [WHERE] @f$ flags are set, all contexts are considered
 *
 *  @paragraph const_file_orientation Controlling the orientation of base pairing
 *  For particular nucleotides that are forced to pair, the following @f$ [ORIENTATION] @f$ flags
 *  may be used:
 *  -#  @f$ U \ldots @f$ Upstream
 *  -#  @f$ D \ldots @f$ Downstream
 *
 *  If no @f$ [ORIENTATION] @f$ flag is set, both directions are considered.
 *
 *  @paragraph const_file_seq_coords Sequence coordinates
 *  Sequence positions of nucleotides/base pairs are @f$ 1- @f$ based and consist of three
 *  positions @f$ i @f$, @f$ j @f$, and @f$ k @f$. Alternativly, four positions may be provided
 *  as a pair of two position ranges @f$ [i:j] @f$, and @f$ [k:l] @f$ using the '-' sign as
 *  delimiter within each range, i.e. @f$ i-j @f$, and @f$ k-l @f$.
 *
 *  @paragraph  const_file_syntax Valid constraint commands
 *  Below are resulting general cases that are considered @em valid constraints:
 *
 *  -#  @b "Forcing a range of nucleotide positions to be paired":\n
 *      Syntax: @code F i 0 k [WHERE] [ORIENTATION] @endcode\n
 *      Description:\n
 *      Enforces the set of @f$ k @f$ consecutive nucleotides starting at
 *      position @f$ i @f$ to be paired. The optional loop type specifier @f$ [WHERE] @f$
 *      allows to force them to appear as closing/enclosed pairs of certain types of loops.
 *  -#  @b "Forcing a set of consecutive base pairs to form":\n
 *      Syntax: @verbatim F i j k [WHERE] @endverbatim\n
 *      Description:\n
 *      Enforces the base pairs @f$ (i,j), \ldots, (i+(k-1), j-(k-1)) @f$ to form.
 *      The optional loop type specifier @f$ [WHERE] @f$ allows to specify in which loop
 *      context the base pair must appear.
 *  -#  @b "Prohibiting a range of nucleotide positions to be paired":\n
 *      Syntax: @verbatim P i 0 k [WHERE] @endverbatim\n
 *      Description:\n
 *      Prohibit a set of @f$ k @f$ consecutive nucleotides to participate
 *      in base pairing, i.e. make these positions unpaired. The optional loop type specifier
 *      @f$ [WHERE] @f$ allows to specify in which loop type context they must not appear paired.
 *  -#  @b "Probibiting a set of consecutive base pairs to form":\n
 *      Syntax: @verbatim P i j k [WHERE] @endverbatim\n
 *      Description:\n
 *      Probibit the base pairs @f$ (i,j), \ldots, (i+(k-1), j-(k-1)) @f$ to form.
 *      The optional loop type specifier @f$ [WHERE] @f$ allows to specify the type of
 *      loop they are disallowed to be the closing or an enclosed pair of.
 *  -#  @b "Prohibiting two ranges of nucleotides to pair with each other":\n
 *      Syntax: @verbatim P i-j k-l [WHERE] @endverbatim
 *      Description:\n
 *      Prohibit any nucleotide @f$ p \in [i:j] @f$ to pair with any other nucleotide
 *      @f$ q \in [k:l] @f$. The optional loop type specifier @f$ [WHERE] @f$ allows to
 *      specify the type of loop they are disallowed to be the closing or an enclosed pair of.
 *  -#  @b "Weakly enforce a range of nucleotide positions to be unpaired":\n
 *      Syntax: @verbatim W i 0 k [WHERE] @endverbatim
 *      Description:\n
 *      This command is meant as a complement to @em prohibiting nucleotides to be paired,
 *      as described above. It also marks the corresponding nucleotides to be unpaired, however,
 *      they are not required to appear in the optional loop type context, if an energetically
 *      better structure includes them as unpaired nucleotides within another loop.
 *  -#  @b "Weakly enforce a set of consecutive base pairs":\n
 *      Syntax: @verbatim W i j k @endverbatim\n
 *      Description:\n
 *      Remove all base pairs that conflict with a set of consecutive base pairs
 *      @f$ (i,j), \ldots, (i+(k-1), j-(k-1)) @f$. Two base pairs @f$ (i,j) @f$ and @f$ (p,q) @f$
 *      conflict with each other if @f$ i < p < j < q @f$, or @f$ p < i < q < j @f$.
 */

#include <stdio.h>

#include <ViennaRNA/data_structures.h>

/**
 *  @brief Print a secondary structure as helix list
 *
 *  @param  db    The structure in dot-bracket format
 *  @param  file  The file handle used to print to (print defaults to 'stdout' if(file == NULL) )
 */
void vrna_structure_print_helix_list(const char *db, FILE *file);

/**
 *  @brief Print a secondary structure as connect table
 *
 *  @param  seq         The RNA sequence
 *  @param  db          The structure in dot-bracket format
 *  @param  energy      The free energy of the structure
 *  @param  identifier  An optional identifier for the sequence
 *  @param  file  The file handle used to print to (print defaults to 'stdout' if(file == NULL) )
 */
void vrna_structure_print_ct( const char *seq,
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
void vrna_structure_print_bpseq(const char *seq,
                                const char *db,
                                FILE *file);

#if WITH_JSON_SUPPORT

void vrna_structure_print_json( const char *seq,
                                const char *db,
                                double energy,
                                const char *identifier,
                                FILE *file);

#endif

/**
 *  @brief  Get a (fasta) data set from a file or stdin
 * 
 *  This function may be used to obtain complete datasets from a filehandle or stdin.
 *  A dataset is always defined to contain at least a sequence. If data starts with a
 *  fasta header, i.e. a line like
 *  @verbatim >some header info @endverbatim
 *  then vrna_read_fasta_record() will assume that the sequence that follows the header may span
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
while(!(vrna_read_fasta_record(&id, &seq, &rest, NULL, 0) & (VRNA_INPUT_ERROR | VRNA_INPUT_QUIT))){
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
 *  In the example above, the while loop will be terminated when vrna_read_fasta_record() returns
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
unsigned int vrna_read_fasta_record(char **header,
                                    char **sequence,
                                    char  ***rest,
                                    FILE *file,
                                    unsigned int options);

/* @brief Extract a dot-bracket structure string from (multiline)character array
 *
 * This function extracts a dot-bracket structure string from the 'rest' array as
 * returned by vrna_read_fasta_record() and returns it. All occurences of comments within the
 * 'lines' array will be skipped as long as they do not break the structure string.
 * If no structure could be read, this function returns NULL.
 *
 * @see vrna_read_fasta_record()
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
 *  @brief  Extract a hard constraint encoded as pseudo dot-bracket string
 *
 *  @pre      The argument 'lines' has to be a 2-dimensional character array as obtained
 *            by vrna_read_fasta_record()
 *  @see      vrna_read_fasta_record(), #VRNA_CONSTRAINT_DB_PIPE, #VRNA_CONSTRAINT_DB_DOT, #VRNA_CONSTRAINT_DB_X
 *            #VRNA_CONSTRAINT_DB_ANG_BRACK, #VRNA_CONSTRAINT_DB_RND_BRACK
 *
 *  @param  cstruc  A pointer to a character array that is used as pseudo dot-bracket
 *                  output
 *  @param  lines   A 2-dimensional character array with the extension lines from the FASTA
 *                  input
 *  @param  option  The option flags that define the behavior and recognition pattern of
 *                  this function
 */
void vrna_extract_record_rest_constraint( char **cstruc,
                                          const char **lines,
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
int vrna_read_SHAPE_file( const char *file_name,
                          int length,
                          double default_value,
                          char *sequence,
                          double *values);

/**
 *  @brief  Read constraints from an input file
 *
 *  This function reads constraint definitions from a file and converts them
 *  into an array of #plist data structures. The data fields of each individual
 *  returned plist entry may adopt the following configurations:
 *  - plist.i == plist.j @f$ \rightarrow @f$ single nucleotide constraint
 *  - plist.i != plist.j @f$ \rightarrow @f$ base pair constraint
 *  - plist.i == 0 @f$ \rightarrow @f$ End of list
 *
 */
plist *vrna_read_constraints_file(const char *filename,
                                  unsigned int length,
                                  unsigned int options);

#ifdef  VRNA_BACKWARD_COMPAT

/* @brief Extract a dot-bracket structure string from (multiline)character array
 *
 * @deprecated This function is deprecated! Use vrna_extract_record_rest_structure() as a replacment.
 */
DEPRECATED(char *extract_record_rest_structure( const char **lines,
                                                unsigned int length,
                                                unsigned int option));

/**
 *  @brief  Get a data record from stdin
 * 
 *  @deprecated This function is deprecated! Use vrna_read_fasta_record() as a replacment.
 *
 */
DEPRECATED(unsigned int read_record(char **header,
                                    char **sequence,
                                    char  ***rest,
                                    unsigned int options));


DEPRECATED(unsigned int get_multi_input_line(char **string, unsigned int options));

#endif

#endif
