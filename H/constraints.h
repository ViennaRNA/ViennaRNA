#ifndef __VIENNA_RNA_PACKAGE_CONSTRAINTS_H__
#define __VIENNA_RNA_PACKAGE_CONSTRAINTS_H__

#include "data_structures.h"

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
#define VRNA_CONSTRAINT_ALL               128U


#define VRNA_CONSTRAINT_DB                256U

#define VRNA_CONSTRAINT_FILE              512U

#define VRNA_CONSTRAINT_IINDX             1024U


#define IN_EXT_LOOP     (char)0x01
#define IN_HP_LOOP      (char)0x02
#define IN_INT_LOOP     (char)0x04
#define IN_MB_LOOP      (char)0x08
#define IN_INT_LOOP_ENC (char)0x10
#define IN_MB_LOOP_ENC  (char)0x20

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
 *  \brief Print structure constraint characters to stdout
 *  (full constraint support)
 *
 */
void print_tty_constraint_full(void);

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
void constrain_ptypes(const char *constraint,
                      unsigned int length,
                      char *ptype,
                      int min_loop_size,
                      unsigned int idx_type);

void getConstraint( char **cstruc,
                    const char **lines,
                    unsigned int option);

hard_constraintT  *get_hard_constraints(  const char *constraint,
                                          unsigned int n,
                                          char *ptype,
                                          unsigned int min_loop_size,
                                          unsigned int options);

void destroy_hard_constraints(hard_constraintT *hc);


#endif
