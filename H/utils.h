#ifndef __VIENNA_RNA_PACKAGE_UTILS_H__
#define __VIENNA_RNA_PACKAGE_UTILS_H__

/** \file
***
**/

/** Output flag of \func get_input_line():  "An ERROR has occured, maybe EOF" **/
#define VRNA_INPUT_ERROR                  1U
/** Output flag of \func get_input_line():  "the user requested quitting the program" **/
#define VRNA_INPUT_QUIT                   2U
/** Output flag of \func get_input_line():  "something was read" **/
#define VRNA_INPUT_MISC                   4U
/** Output flag of \func get_input_line():  "a fasta header was read" **/
#define VRNA_INPUT_FASTA_HEADER           8U
/** **/
#define VRNA_INPUT_SEQUENCE               16U
/** **/
#define VRNA_INPUT_STRUCTURE              32U
/** Input switch for \func get_input_line():  "do not print the line read to stdout" **/
#define VRNA_INPUT_NOPRINT                64U
/** Input switch for \func get_input_line():  "do not skip comment lines" **/
#define VRNA_INPUT_NOSKIP_COMMENTS        128U
/** Input switch for \func get_input_line():  "do not eliminate white spaces at end of line" **/
#define VRNA_INPUT_NOELIM_WS_SUFFIX       256U

/** use the pipe sign '|' for structure constraints **/
#define VRNA_CONSTRAINT_PIPE              1U


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
/*@only@*/ /*@notnull@*/
void  *space(unsigned size) /*@ensures MaxSet(result) == (size-1);@*/;
			    /* allocate space safely */
/*@only@*/ /*@notnull@*/
void  *xrealloc(/*@null@*/ /*@only@*/ /*@out@*/ /*@returned@*/ void *p, unsigned size) /*@modifies *p @*/ /*@ensures MaxSet(result) == (size-1) @*/;
#endif

/**
*** Die with an error message
*** \see 
*** \param message The error message to be printed before exiting with 'FAILURE'
**/
/*@exits@*/
void nrerror(const char message[]);

/**
*** Print a warning message
**/
void warn_user(const char message[]);

void   init_rand(void);                /* make random number seeds */
extern unsigned short xsubi[3];               /* current 48bit random number */
double urn(void);                      /* random number from [0..1] */
int    int_urn(int from, int to);      /* random integer */
void   filecopy(FILE *from, FILE *to); /* inefficient `cp' */
/*@observer@*/ char  *time_stamp(void);               /* current date in a string */
/*@only@*/ /*@notnull@*/ char  *random_string(int l, const char symbols[]);
/* random string of length l using characters from symbols[] */
/**
*** calculate hamming distance
**/
int    hamming(const char *s1, const char *s2);
/*@only@*/ /*@null@*/
char  *get_line(FILE *fp); /* read one (arbitrary length) line from fp */

int skip_comment_lines(char **line);

/**
*** Retrieve a line from 'stdin' savely while skipping comment characters and
*** other features
*** This function returns the type of input it has read if recognized.
*** An option argument allows to switch between different reading modes.\n
*** Currently available options are:\n
*** \ref VRNA_INPUT_NOPRINT, \ref VRNA_INPUT_NOSKIP_COMMENTS, ref\ VRNA_INPUT_NOELIM_WS_SUFFIX
***
*** pass a collection of options as one value like this:
*** \verbatim get_input_line(string, option_1 | option_2 | option_n) \endverbatim
***
*** If the function recognizes the type of input, it will report it in the return
*** value. It also reports if a user defined 'quit' command (@-sign on 'stdin')
*** was given. Possible return values are:\n
*** \ref VRNA_INPUT_FASTA_HEADER, \ref VRNA_INPUT_ERROR, \ref VRNA_INPUT_MISC, \ref VRNA_INPUT_QUIT
***
*** \param string   A pointer to the character array that contains the line read
*** \param options  A collection of options for switching the functions behavior
*** \return         A flag with information about what has been read
**/
unsigned int get_intput_line(char **string, unsigned int options);


/**
*** pack secondary secondary structure, 5:1 compression using base 3 encoding
**/
char *pack_structure(const char *struc);
/**
*** unpack sec structure packed with pack_structure()
**/
char *unpack_structure(const char *packed);
/**
*** returns a newly allocated table, such that:  table[i]=j if (i.j) pair or
*** 0 if i is unpaired, table[0] contains the length of the structure.
**/
short *make_pair_table(const char *structure);
/**
*** dist = {number of base pairs in one structure but not in the other}
*** same as edit distance with open-pair close-pair as move-set
**/
int bp_distance(const char *str1, const char *str2);

/**
*** Just print a line to stdout that asks for an input sequence
*** There will also be a scale line printed that helps orientation of the sequence positions
**/
void print_tty_input_seq(void);
void print_tty_input_seq_str(const char *s);

void print_tty_constraint_str(void);

void str_DNA2RNA(char *sequence);
void str_RNA2RNA(char *sequence);


#endif
