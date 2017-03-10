#ifndef VIENNA_RNA_PACKAGE_UTILS_H
#define VIENNA_RNA_PACKAGE_UTILS_H

/* make this interface backward compatible with RNAlib < 2.2.0 */
#define VRNA_BACKWARD_COMPAT

#ifdef VRNA_WARN_DEPRECATED
# ifdef __GNUC__
#  define DEPRECATED(func) func __attribute__ ((deprecated))
# else
#  define DEPRECATED(func) func
# endif
#else
# define DEPRECATED(func) func
#endif

/**
 *  @file     utils.h
 *  @ingroup  utils
 *  @brief    General utility- and helper-functions used throughout the @em ViennaRNA @em Package
 */

/**
 *  @{
 *  @ingroup  utils
 */

#include <stdio.h>
#include <stdarg.h>

#include <ViennaRNA/data_structures.h>
#include <ViennaRNA/string_utils.h>
#include <ViennaRNA/structure_utils.h>
#include <ViennaRNA/file_utils.h>
#include <ViennaRNA/alphabet.h>

/* two helper macros to indicate whether a function should be exported in
the library or stays hidden */
#define PUBLIC
#define PRIVATE static

/**
 *  @brief Output flag of get_input_line():  @e "An ERROR has occured, maybe EOF"
 */
#define VRNA_INPUT_ERROR                  1U
/**
 *  @brief @brief Output flag of get_input_line():  @e "the user requested quitting the program"
 */
#define VRNA_INPUT_QUIT                   2U
/**
 *  @brief Output flag of get_input_line():  @e "something was read"
 */
#define VRNA_INPUT_MISC                   4U

/**
 *  @brief  Input/Output flag of get_input_line():\n
 *  if used as input option this tells get_input_line() that the data to be read should comply
 *  with the FASTA format
 * 
 *  the function will return this flag if a fasta header was read
 */
#define VRNA_INPUT_FASTA_HEADER           8U

/*
 *  @brief  Input flag for get_input_line():\n
 *  Tell get_input_line() that we assume to read a nucleotide sequence
 * 
 */
#define VRNA_INPUT_SEQUENCE               16U

/** @brief  Input flag for get_input_line():\n
 *  Tell get_input_line() that we assume to read a structure constraint
 * 
 */
#define VRNA_INPUT_CONSTRAINT             32U

/**
 *  @brief  Input switch for get_input_line():
 *  @e "do not trunkate the line by eliminating white spaces at end of line"
 */
#define VRNA_INPUT_NO_TRUNCATION          256U

/**
 *  @brief  Input switch for vrna_file_fasta_read_record():  @e "do fill rest array"
 */
#define VRNA_INPUT_NO_REST                512U

/**
 *  @brief  Input switch for vrna_file_fasta_read_record():  @e "never allow data to span more than one line"
 */
#define VRNA_INPUT_NO_SPAN                1024U

/**
 *  @brief  Input switch for vrna_file_fasta_read_record():  @e "do not skip empty lines"
 */
#define VRNA_INPUT_NOSKIP_BLANK_LINES     2048U

/**
 *  @brief  Output flag for vrna_file_fasta_read_record():  @e "read an empty line"
 */
#define VRNA_INPUT_BLANK_LINE             4096U

/**
 *  @brief Input switch for get_input_line():  @e "do not skip comment lines"
 */
#define VRNA_INPUT_NOSKIP_COMMENTS        128U

/**
 *  @brief  Output flag for vrna_file_fasta_read_record():  @e "read a comment"
 */
#define VRNA_INPUT_COMMENT                8192U

/**
 *  @brief Get the minimum of two comparable values
 */
#define MIN2(A, B)      ((A) < (B) ? (A) : (B))

/**
 *  @brief Get the maximum of two comparable values
 */
#define MAX2(A, B)      ((A) > (B) ? (A) : (B))

/**
 *  @brief Get the minimum of three comparable values
 */
#define MIN3(A, B, C)   (MIN2(  (MIN2((A),(B))) ,(C)))

/**
 *  @brief Get the maximum of three comparable values
 */
#define MAX3(A, B, C)   (MAX2(  (MAX2((A),(B))) ,(C)))


#ifdef WITH_DMALLOC
/* use dmalloc library to check for memory management bugs */
#include "dmalloc.h"
#define vrna_alloc(S)       calloc(1,(S))
#define vrna_realloc(p, S)  xrealloc(p, S)
#else

/**
 *  @brief Allocate space safely
 *
 *  @param size The size of the memory to be allocated in bytes
 *  @return     A pointer to the allocated memory
 */
void  *vrna_alloc(unsigned size);

/**
 *  @brief Reallocate space safely
 *
 *  @param p    A pointer to the memory region to be reallocated
 *  @param size The size of the memory to be allocated in bytes
 *  @return     A pointer to the newly allocated memory
 */
void  *vrna_realloc(void *p, unsigned size);

#endif

/**
 *  @brief Print an error message and die
 *
 *  This function is a wrapper to @em fprintf(stderr, ...) that
 *  puts a capital <b>ERROR:</b> in front of the message and then exits
 *  the calling program.
 *
 *  @see vrna_message_verror(), vrna_message_warning(), vrna_message_info()
 *
 *  @param format The error message to be printed
 *  @param ...    Optional arguments for the formatted message string
 */
void vrna_message_error(const char *format, ...);


/**
 *  @brief Print an error message and die
 *
 *  This function is a wrapper to @em vfprintf(stderr, ...) that
 *  puts a capital <b>ERROR:</b> in front of the message and then exits
 *  the calling program.
 *
 *  @see vrna_message_error(), vrna_message_warning(), vrna_message_info()
 *
 *  @param format The error message to be printed
 *  @param args   The argument list for the formatted message string
 */
void vrna_message_verror(const char *format, va_list args);


/**
 *  @brief Print a warning message
 *
 *  This function is a wrapper to @em fprintf(stderr, ...) that
 *  puts a capital <b>WARNING:</b> in front of the message.
 *
 *  @see vrna_message_vwarning(), vrna_message_error(), vrna_message_info()
 *
 *  @param format The warning message to be printed
 *  @param ...    Optional arguments for the formatted message string
 */
void vrna_message_warning(const char *format, ...);


/**
 *  @brief Print a warning message
 *
 *  This function is a wrapper to @em fprintf(stderr, ...) that
 *  puts a capital <b>WARNING:</b> in front of the message.
 *
 *  @see vrna_message_vwarning(), vrna_message_error(), vrna_message_info()
 *
 *  @param format The warning message to be printed
 *  @param args   The argument list for the formatted message string
 */
void vrna_message_vwarning(const char *format, va_list args);


/**
 *  @brief Print an info message
 *
 *  This function is a wrapper to @em fprintf(...).
 *
 *  @see vrna_message_vinfo(), vrna_message_error(), vrna_message_warning()
 *
 *  @param fp     The file pointer where the message is printed to
 *  @param format The warning message to be printed
 *  @param ...    Optional arguments for the formatted message string
 */
void vrna_message_info(FILE *fp, const char *format, ...);


/**
 *  @brief Print an info message
 *
 *  This function is a wrapper to @em fprintf(...).
 *
 *  @see vrna_message_vinfo(), vrna_message_error(), vrna_message_warning()
 *
 *  @param fp     The file pointer where the message is printed to
 *  @param format The info message to be printed
 *  @param args   The argument list for the formatted message string
 */
void vrna_message_vinfo(FILE *fp, const char *format, va_list args);


/**
 *  @brief  Initialize seed for random number generator
 */
void vrna_init_rand(void);

/**
 * @brief Current 48 bit random number
 *
 *  This variable is used by vrna_urn(). These should be set to some
 *  random number seeds before the first call to vrna_urn().
 *
 *  @see vrna_urn()
 */
extern unsigned short xsubi[3];

/**
 *  @brief get a random number from [0..1]
 *
 *  @see  vrna_int_urn(), vrna_init_rand()
 *  @note Usually implemented by calling @e erand48().
 *  @return   A random number in range [0..1]
 */
double vrna_urn(void);

/**
 *  @brief Generates a pseudo random integer in a specified range
 *
 *  @see  vrna_urn(), vrna_init_rand()
 *  @param from   The first number in range
 *  @param to     The last number in range
 *  @return       A pseudo random number in range [from, to]
 */
int vrna_int_urn(int from, int to);

/**
 *  @brief Get a timestamp
 *
 *  Returns a string containing the current date in the format
 *  @verbatim Fri Mar 19 21:10:57 1993 @endverbatim
 *
 *  @return A string containing the timestamp
 */
char  *vrna_time_stamp(void);

/**
 *  Retrieve a line from 'stdin' savely while skipping comment characters and
 *  other features
 *  This function returns the type of input it has read if recognized.
 *  An option argument allows one to switch between different reading modes.\n
 *  Currently available options are:\n
 *  #VRNA_INPUT_NOPRINT_COMMENTS, #VRNA_INPUT_NOSKIP_COMMENTS, #VRNA_INPUT_NOELIM_WS_SUFFIX
 * 
 *  pass a collection of options as one value like this:
 *  @verbatim get_input_line(string, option_1 | option_2 | option_n) @endverbatim
 * 
 *  If the function recognizes the type of input, it will report it in the return
 *  value. It also reports if a user defined 'quit' command (@-sign on 'stdin')
 *  was given. Possible return values are:\n
 *  #VRNA_INPUT_FASTA_HEADER, #VRNA_INPUT_ERROR, #VRNA_INPUT_MISC, #VRNA_INPUT_QUIT
 * 
 *  @param string   A pointer to the character array that contains the line read
 *  @param options  A collection of options for switching the functions behavior
 *  @return         A flag with information about what has been read
 */
unsigned int get_input_line(char **string,
                            unsigned int options);


/**
 *  @brief Print a line to @e stdout that asks for an input sequence
 *
 *  There will also be a ruler (scale line) printed that helps orientation of the sequence positions
 */
void vrna_message_input_seq_simple(void);


/**
 *  @brief Print a line with a user defined string and a ruler to stdout.
 *
 *  (usually this is used to ask for user input)
 *  There will also be a ruler (scale line) printed that helps orientation of the sequence positions
 * 
 *  @param s A user defined string that will be printed to stdout
 */
void vrna_message_input_seq(const char *s);

/**
 *  @brief Get an index mapper array (iindx) for accessing the energy matrices, e.g. in partition function related functions.
 *
 *  Access of a position "(i,j)" is then accomplished by using @verbatim (i,j) ~ iindx[i]-j @endverbatim
 *  This function is necessary as most of the two-dimensional energy matrices are actually one-dimensional arrays throughout
 *  the ViennaRNA Package
 * 
 *  Consult the implemented code to find out about the mapping formula ;)
 * 
 *  @see vrna_idx_col_wise()
 *  @param length The length of the RNA sequence
 *  @return       The mapper array
 */
int *vrna_idx_row_wise(unsigned int length);

/**
 *  @brief Get an index mapper array (indx) for accessing the energy matrices, e.g. in MFE related functions.
 *
 *  Access of a position "(i,j)" is then accomplished by using @verbatim (i,j) ~ indx[j]+i @endverbatim
 *  This function is necessary as most of the two-dimensional energy matrices are actually one-dimensional arrays throughout
 *  the ViennaRNAPackage
 * 
 *  Consult the implemented code to find out about the mapping formula ;)
 * 
 *  @see vrna_idx_row_wise()
 *  @param length The length of the RNA sequence
 *  @return       The mapper array
 * 
 */
int *vrna_idx_col_wise(unsigned int length);

/**
 *  @}
 */

#ifdef  VRNA_BACKWARD_COMPAT

DEPRECATED(int   *get_indx(unsigned int length));

DEPRECATED(int   *get_iindx(unsigned int length));

/**
 *  @brief Read a line of arbitrary length from a stream
 *
 *  Returns a pointer to the resulting string. The necessary memory is
 *  allocated and should be released using @e free() when the string is
 *  no longer needed.
 *
 *	@deprecated	Use vrna_read_line() as a substitute!
 *
 *  @param  fp  A file pointer to the stream where the function should read from
 *  @return     A pointer to the resulting string
 */
DEPRECATED(char  *get_line(FILE *fp));

/**
 *  @brief Print a line to @e stdout that asks for an input sequence
 *
 *  There will also be a ruler (scale line) printed that helps orientation of the sequence positions
 *  @deprecated Use vrna_message_input_seq_simple() instead!
 */
DEPRECATED(void print_tty_input_seq(void));

/**
 *  @brief Print a line with a user defined string and a ruler to stdout.
 *
 *  (usually this is used to ask for user input)
 *  There will also be a ruler (scale line) printed that helps orientation of the sequence positions
 * 
 *  @deprecated Use vrna_message_input_seq() instead!
 */
DEPRECATED(void print_tty_input_seq_str(const char *s));

/**
 *  @brief Print a warning message
 *
 *  Print a warning message to @e stderr
 *
 *  @deprecated Use vrna_message_warning() instead!
 */
DEPRECATED(void warn_user(const char message[]));

/**
 *  @brief Die with an error message
 *
 *  @deprecated Use vrna_message_error() instead!
 */
DEPRECATED(void nrerror(const char message[]));

/**
 *  @brief Allocate space safely
 *
 *  @deprecated Use vrna_alloc() instead!
 */
DEPRECATED(void *space(unsigned size));

/**
 *  @brief Reallocate space safely
 *
 *  @deprecated Use vrna_realloc() instead!
 */
DEPRECATED(void *xrealloc(void *p, unsigned size));

/**
 *  @brief  Make random number seeds
 *  @deprecated Use vrna_init_rand() instead!
 */
DEPRECATED(void init_rand(void));

/**
 *  @brief get a random number from [0..1]
 *
 *  @deprecated Use vrna_urn() instead!
 */
DEPRECATED(double urn(void));

/**
 *  @brief Generates a pseudo random integer in a specified range
 *
 *  @deprecated Use vrna_int_urn() instead!
 */
DEPRECATED(int int_urn(int from, int to));

/**
 *  @brief  Inefficient `cp`
 *
 *  @deprecated Use vrna_file_copy() instead!
 */
DEPRECATED(void filecopy(FILE *from, FILE *to));

/**
 *  @brief Get a timestamp
 *
 *  @deprecated Use vrna_time_stamp() instead!
 */
DEPRECATED(char *time_stamp(void));

#endif

#endif
