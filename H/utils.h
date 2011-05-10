#ifndef __VIENNA_RNA_PACKAGE_UTILS_H__
#define __VIENNA_RNA_PACKAGE_UTILS_H__

/**
 *  \file utils.h
 *  \brief Various utility- and helper-functions used throughout the Vienna RNA package
 */

/**
 *  Output flag of \ref get_input_line():  "An ERROR has occured, maybe EOF"
 */
#define VRNA_INPUT_ERROR                  1U
/**
 *  Output flag of \ref get_input_line():  "the user requested quitting the program"
 */
#define VRNA_INPUT_QUIT                   2U
/**
 *  Output flag of \ref get_input_line():  "something was read"
 */
#define VRNA_INPUT_MISC                   4U
/** Input/Output flag of \ref get_input_line():\n
 *  if used as input option this tells get_input_line() that the data to be read should comply
 *  with the FASTA format
 * 
 *  the function will return this flag if a fasta header was read
 */
#define VRNA_INPUT_FASTA_HEADER           8U
/** Input flag for get_input_line():\n
 *  Tell get_input_line() that we assume to read a nucleotide sequence
 * 
 */
#define VRNA_INPUT_SEQUENCE               16U
/** Input flag for get_input_line():\n
 *  Tell get_input_line() that we assume to read a structure constraint
 * 
 */
#define VRNA_INPUT_CONSTRAINT             32U
/**
 *  Input switch for \ref get_input_line():
 *  "do not trunkate the line by eliminating white spaces at end of line"
 */
#define VRNA_INPUT_NO_TRUNCATION          256U
/**
 *  Input switch for read_record():  "do fill rest array"
 */
#define VRNA_INPUT_NO_REST                512U
/**
 *  Input switch for read_record():  "never allow data to span more than one line"
 */
#define VRNA_INPUT_NO_SPAN                1024U
/**
 *  Input switch for read_record():  "do not skip empty lines"
 */
#define VRNA_INPUT_NOSKIP_BLANK_LINES     2048U
/**
 *  Output flag for read_record():  "read an empty line"
 */
#define VRNA_INPUT_BLANK_LINE             4096U
/**
 *  Input switch for \ref get_input_line():  "do not skip comment lines"
 */
#define VRNA_INPUT_NOSKIP_COMMENTS        128U
/**
 *  Output flag for read_record():  "read a comment"
 */
#define VRNA_INPUT_COMMENT                8192U




/**
 *  pipe sign '|' switch for structure constraints (paired with another base)
 */
#define VRNA_CONSTRAINT_PIPE              1U
/**
 *  dot '.' switch for structure constraints (no constraint at all)
 */
#define VRNA_CONSTRAINT_DOT               2U
/**
 *  'x' switch for structure constraint (base must not pair)
 */
#define VRNA_CONSTRAINT_X                 4U
/**
 *  angle brackets '<', '>' switch for structure constraint (paired downstream/upstream)
 */
#define VRNA_CONSTRAINT_ANG_BRACK         8U
/**
 *  round brackets '(',')' switch for structure constraint (base i pairs base j)
 */
#define VRNA_CONSTRAINT_RND_BRACK         16U
/**
 *  constraint may span over several lines
 */
#define VRNA_CONSTRAINT_MULTILINE         32U
/**
 *  do not print the header information line
 */
#define VRNA_CONSTRAINT_NO_HEADER         64U
/**
 *  placeholder for all constraining characters
 */
#define VRNA_CONSTRAINT_ALL              128U

/**
 *  Get the minimum of two comparable values
 */
#define MIN2(A, B)      ((A) < (B) ? (A) : (B))
/**
 *  Get the maximum of two comparable values
 */
#define MAX2(A, B)      ((A) > (B) ? (A) : (B))
/**
 *  Get the minimum of three comparable values
 */
#define MIN3(A, B, C)   (MIN2(  (MIN2((A),(B))) ,(C)))
/**
 *  Get the maximum of three comparable values
 */
#define MAX3(A, B, C)   (MAX2(  (MAX2((A),(B))) ,(C)))

#ifdef HAVE_CONFIG_H
#include <config.h>
#ifndef HAVE_STRDUP
char *strdup(const char *s);
#endif
#endif
#ifdef WITH_DMALLOC
/* use dmalloc library to check for memory management bugs */
#include "dmalloc.h"
#define space(S) calloc(1,(S))
#else

/**
 *  \brief Allocate space safely
 *
 *  \param size The size of the memory to be allocated in bytes
 *  \return     A pointer to the allocated memory
 */
/*@only@*/ /*@notnull@*/
void  *space(unsigned size) /*@ensures MaxSet(result) == (size-1);@*/;

/**
 *  \brief Reallocate space safely
 *
 *  \param p    A pointer to the memory region to be reallocated
 *  \param size The size of the memory to be allocated in bytes
 *  \return     A pointer to the newly allocated memory
 */
/*@only@*/ /*@notnull@*/
void  *xrealloc(/*@null@*/ /*@only@*/ /*@out@*/ /*@returned@*/ void *p, unsigned size) /*@modifies *p @*/ /*@ensures MaxSet(result) == (size-1) @*/;
#endif

/**
 *  \brief Die with an error message
 *
 *  \see warn_user()
 *  \param message The error message to be printed before exiting with 'FAILURE'
 */
/*@exits@*/
void nrerror(const char message[]);

/**
 *  \brief Print a warning message
 *
 *  Print a warning message to \e stderr
 *
 *  \param  message   The warning message
 */
void warn_user(const char message[]);

/**
 *  \brief  Make random number seeds
 */
void   init_rand(void);

/**
 * \brief Current 48 bit random number
 *
 *  This variable is used by urn(). These should be set to some
 *  random number seeds before the first call to urn().
 *
 *  \see urn()
 */
extern unsigned short xsubi[3];

/**
 *  \brief get a random number from [0..1]
 *
 *  \note Usually implemented by calling \e erand48().
 *  \return   A random number in range [0..1]
 */
double urn(void);

/**
 *  \brief Generates a pseudo random integer in a specified range
 *
 *  \param from   The first number in range
 *  \param to     The last number in range
 *  \return       A pseudo random number in range [from, to]
 */
int    int_urn(int from, int to);

void   filecopy(FILE *from, FILE *to); /* inefficient `cp' */

/**
 *  \brief Get a timestamp
 *
 *  Returns a string containing the current date in the format
 *  \verbatim Fri Mar 19 21:10:57 1993\endverbatim
 *
 *  \return A string containing the timestamp
 */
/*@observer@*/
char  *time_stamp(void);

/**
 *  \brief Create a random string using characters from a specified symbol set
 *
 *  \param l        The length of the sequence
 *  \param symbols  The symbol set
 *  \return         A random string of length 'l' containing characters from the symbolset
 */
/*@only@*/ /*@notnull@*/
char  *random_string(int l, const char symbols[]);

/**
 *  \brief Calculate hamming distance between two sequences
 *
 *  Calculate the number of positions in which 
 *  \param s1   The first sequence
 *  \param s2   The second sequence
 *  \return     The hamming distance between s1 and s2
 */
int   hamming(const char *s1, const char *s2);

/**
 *  \brief Calculate hamming distance between two sequences up to a specified length
 *
 *  This function is similar to hamming() but instead of comparing both sequences
 *  up to their actual length only the first 'n' characters are taken into account
 *  \param s1   The first sequence
 *  \param s2   The second sequence
 *  \return     The hamming distance between s1 and s2
 */
int   hamming_bound(const char *s1, const char *s2, int n);

/**
 *  \brief Read a line of arbitrary length from a stream
 *
 *  Returns a pointer to the resulting string. The necessary memory is
 *  allocated and should be released using \e free() when the string is
 *  no longer needed.
 *
 *  \param  fp  A file pointer to the stream where the function should read from
 *  \return     A pointer to the resulting string
 */
/*@only@*/ /*@null@*/
char  *get_line(FILE *fp);

int skip_comment_lines(char **line);

/**
 *  Retrieve a line from 'stdin' savely while skipping comment characters and
 *  other features
 *  This function returns the type of input it has read if recognized.
 *  An option argument allows to switch between different reading modes.\n
 *  Currently available options are:\n
 *  #VRNA_INPUT_NOPRINT_COMMENTS, #VRNA_INPUT_NOSKIP_COMMENTS, #VRNA_INPUT_NOELIM_WS_SUFFIX
 * 
 *  pass a collection of options as one value like this:
 *  \verbatim get_input_line(string, option_1 | option_2 | option_n) \endverbatim
 * 
 *  If the function recognizes the type of input, it will report it in the return
 *  value. It also reports if a user defined 'quit' command (@-sign on 'stdin')
 *  was given. Possible return values are:\n
 *  #VRNA_INPUT_FASTA_HEADER, #VRNA_INPUT_ERROR, #VRNA_INPUT_MISC, #VRNA_INPUT_QUIT
 * 
 *  \param string   A pointer to the character array that contains the line read
 *  \param options  A collection of options for switching the functions behavior
 *  \return         A flag with information about what has been read
 */
unsigned int get_input_line(char **string, unsigned int options);

unsigned int get_multi_input_line(char **string, unsigned int options);

/**
 *  \brief  Get a data record from stdin
 * 
 *  This function may be used to obtain complete datasets from stdin. A dataset is always
 *  defined to contain at least a sequence. If data on stdin starts with a fasta header,
 *  i.e. a line like
 *  \verbatim >some header info \endverbatim
 *  then read_record() will assume that the sequence that follows the header may span
 *  over several lines. To disable this behavior and to assign a single line to the argument
 *  'sequence' one can pass VRNA_INPUT_NO_SPAN in the 'options' argument.
 *  If no fasta header is read in the beginning of a data block, a sequence must not span over
 *  multiple lines!\n
 *  Unless the options #VRNA_INPUT_NOSKIP_COMMENTS or #VRNA_INPUT_NOSKIP_BLANK_LINES are passed,
 *  a sequence may be interrupted by lines starting with a comment character or empty lines.\n
 *  A sequence is regarded as completely read if it was either assumed to not span over multiple
 *  lines, a secondary structure or structure constraint follows the sequence on the next line
 *  or a new header marks the beginning of a new sequence...\n
 *  All lines following the sequence (this includes comments) and not initiating a new dataset are
 *  available through the line-array 'rest'. Here one can usually find the structure constraint or
 *  other information belonging to the current dataset. Filling of 'rest' may be prevented by
 *  passing #VRNA_INPUT_NO_REST to the options argument.\n
 * 
 *  \note This function will exit any program with an error message if no sequence could be read!
 * 
 *  The main purpose of this function is to be able to easily parse blocks of data from stdin
 *  in the header of a loop where all calculations for the appropriate data is done inside the
 *  loop. The loop may be then left on certain return values, e.g.:
 *  \verbatim
char *id, *seq, **rest;
int  i;
while(!(read_record(&id, &seq, &rest, 0) & (VRNA_INPUT_ERROR | VRNA_INPUT_QUIT))){
  if(id) printf("%s\n", id);
  printf("%s\n", seq);
  if(rest)
    for(i=0;rest[i];i++)
      printf("%s\n", rest[i]);
} \endverbatim
 * 
 *  In the example above, the while loop will be terminated when read_record() returns either an
 *  error or a user initiated quit request.\n
 *  As long as data is read from stdin, the id is printed if it is available for the current block
 *  of data. The sequence will be printed in any case and if some more lines belong to the current
 *  block of data each line will be printed as well.
 * 
 *  \note Do not forget to free the memory occupied by header, sequence and rest!
 * 
 *  \param  header    A pointer which will be set such that it points to the header of the record
 *  \param  sequence  A pointer which will be set such that it points to the sequence of the record
 *  \param  rest      A pointer which will be set such that it points to an array of lines which also belong to the record
 *  \param  options   Some options which may be passed to alter the behavior of the function, use 0 for no options
 *  \return           A flag with information about what the function actually did read
 */
unsigned int read_record(char **header, char **sequence, char  ***rest, unsigned int options);

/**
 *  \brief Pack secondary secondary structure, 5:1 compression using base 3 encoding
 *
 *  Returns a binary string encoding of the secondary structure using
 *  a 5:1 compression scheme. The string is NULL terminated and can
 *  therefore be used with standard string functions such as strcmp().
 *  Useful for programs that need to keep many structures in memory.
 *
 *  \param struc    The secondary structure in dot-bracket notation
 *  \return         The binary encoded structure
 */
char *pack_structure(const char *struc);

/**
 *  \brief Unpack secondary structure previously packed with pack_structure()
 *
 *  Translate a compressed binary string produced by pack_structure() back into
 *  the familiar dot-bracket notation.
 *
 *  \param packed   The binary encoded packed secondary structure
 *  \return         The unpacked secondary structure in dot-bracket notation
 */
char *unpack_structure(const char *packed);

/**
 *  \brief Create a pair table of a secondary structure
 *
 *  Returns a newly allocated table, such that table[i]=j if (i.j) pair
 *  or 0 if i is unpaired, table[0] contains the length of the structure.
 *
 *  \param  structure The secondary structure in dot-bracket notation
 *  \return           A pointer to the created pair_table
 */
short *make_pair_table(const char *structure);

/**
 *  \brief Get an exact copy of a pair table
 *
 *  \param pt The pair table to be copied
 *  \return   A pointer to the copy of 'pt' 
 */
short *copy_pair_table(const short *pt);

/**
 *  \brief Compute the "base pair" distance between two secondary structures s1 and s2.
 * 
 *  The sequences should have the same length.
 *  dist = number of base pairs in one structure but not in the other
 *  same as edit distance with open-pair close-pair as move-set
 * 
 *  \param str1   First structure in dot-bracket notation
 *  \param str2   Second structure in dot-bracket notation
 *  \return       The base pair distance between str1 and str2
 */
int bp_distance(const char *str1, const char *str2);

/**
 *  \brief Print a line to \e stdout that asks for an input sequence
 *
 *  There will also be a ruler (scale line) printed that helps orientation of the sequence positions
 */
void print_tty_input_seq(void);

/**
 *  \brief Print a line with a user defined string and a ruler to stdout.
 *
 *  (usually this is used to ask for user input)
 *  There will also be a ruler (scale line) printed that helps orientation of the sequence positions
 * 
 *  \param s A user defined string that will be printed to stdout
 */
void print_tty_input_seq_str(const char *s);

/**
 *  \brief Print structure constraint characters to stdout
 *  (full constraint support)
 *
 */
void print_tty_constraint_full(void);

/**
 *  \brief Print structure constraint characters to stdout.
 *  (constraint support is specified by option parameter)
 *
 *  Currently available options are:\n
 *  #VRNA_CONSTRAINT_PIPE (paired with another base)\n
 *  #VRNA_CONSTRAINT_DOT (no constraint at all)\n
 *  #VRNA_CONSTRAINT_X (base must not pair)\n
 *  #VRNA_CONSTRAINT_ANG_BRACK (paired downstream/upstream)\n
 *  #VRNA_CONSTRAINT_RND_BRACK (base i pairs base j)\n
 * 
 *  pass a collection of options as one value like this:
 *  \verbatim print_tty_constraint(option_1 | option_2 | option_n) \endverbatim
 * 
 *  \param option Option switch that tells which constraint help will be printed
 */
void print_tty_constraint(unsigned int option);

/**
 *  \brief Convert a DNA input sequence to RNA alphabet
 *
 *  This function substitudes <i>T</i> and <i>t</i> with <i>U</i> and <i>u</i>, respectively
 * 
 *  \param sequence The sequence to be converted
 */
void str_DNA2RNA(char *sequence);

/**
 *  \brief Convert an input sequence to uppercase
 * 
 *  \param sequence The sequence to be converted
 */
void  str_uppercase(char *sequence);

/**
 *  \brief Get an index mapper array (iindx) for accessing the energy matrices, e.g. in partition function related functions.
 *
 *  Access of a position "(i,j)" is then accomplished by using \verbatim (i,j) ~ iindx[i]-j \endverbatim
 *  This function is necessary as most of the two-dimensional energy matrices are actually one-dimensional arrays throughout
 *  the ViennaRNAPackage
 * 
 *  Consult the implemented code to find out about the mapping formula ;)
 * 
 *  \see get_indx()
 *  \param length The length of the RNA sequence
 *  \return       The mapper array
 */
int   *get_iindx(unsigned int length);

/**
 *  \brief Get an index mapper array (indx) for accessing the energy matrices, e.g. in MFE related functions.
 *
 *  Access of a position "(i,j)" is then accomplished by using \verbatim (i,j) ~ indx[j]+i \endverbatim
 *  This function is necessary as most of the two-dimensional energy matrices are actually one-dimensional arrays throughout
 *  the ViennaRNAPackage
 * 
 *  Consult the implemented code to find out about the mapping formula ;)
 * 
 *  \see get_iindx()
 *  \param length The length of the RNA sequence
 *  \return       The mapper array
 * 
 */
int   *get_indx(unsigned int length);

void getConstraint(char **cstruc, const char **lines, unsigned int option);

/**
 *  \brief Insert constraining pair types according to constraint structure string
 *
 *  \see get_indx(), get_iindx()
 *
 *  \param constraint     The structure constraint string
 *  \param length         The actual length of the sequence (constraint may be shorter)
 *  \param ptype          A pointer to the basepair type array
 *  \param min_loop_size  The minimal loop size (usually \ref TURN )
 *  \param idx_type       Define the access type for base pair type array (0 = indx, 1 = iindx)
 */
void constrain_ptypes(const char *constraint, unsigned int length, char *ptype, int *BP, int min_loop_size, unsigned int idx_type);

unsigned int  *make_referenceBP_array(short *reference_pt, unsigned int turn);

unsigned int  *compute_BPdifferences(short *pt1, short *pt2, unsigned int turn);

#endif
