#ifndef VIENNA_RNA_PACKAGE_CONSTRAINTS_HARD_H
#define VIENNA_RNA_PACKAGE_CONSTRAINTS_HARD_H

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
 *  @file       constraints_hard.h
 *  @ingroup    hard_constraints
 *  @brief      Functions and data structures for handling of secondary structure hard constraints
 */
#include <ViennaRNA/data_structures.h>

/**
 *  @addtogroup hard_constraints
 *
 *  @brief  This module covers all functionality for hard constraints in secondary
 *          structure prediction
 */

/**
 *  @brief Typename for the hard constraints data structure #vrna_hc_s
 *  @ingroup  hard_constraints
 */
typedef struct  vrna_hc_s vrna_hc_t;

/**
 *  @brief Typename for the single nucleotide hard constraint data structure #vrna_hc_up_s
 *  @ingroup  hard_constraints
 */
typedef struct vrna_hc_up_s vrna_hc_up_t;

/**
 * @brief Callback to evaluate whether or not a particular decomposition step is contributing to the solution space
 *
 * @ingroup hard_constraints
 *
 * This is the prototype for callback functions used by the folding recursions to evaluate generic
 * hard constraints. The first four parameters passed indicate the delimiting nucleotide positions
 * of the decomposition, and the parameter @p denotes the decomposition step. The last parameter
 * @p data is the auxiliary data structure associated to the hard constraints via vrna_hc_add_data(),
 * or NULL if no auxiliary data was added.
 *
 * @see #VRNA_DECOMP_PAIR_HP, #VRNA_DECOMP_PAIR_IL, #VRNA_DECOMP_PAIR_ML, #VRNA_DECOMP_ML_ML_ML,
 *      #VRNA_DECOMP_ML_STEM, #VRNA_DECOMP_ML_ML, #VRNA_DECOMP_ML_UP, #VRNA_DECOMP_ML_ML_STEM,
 *      #VRNA_DECOMP_ML_COAXIAL, #VRNA_DECOMP_EXT_EXT, #VRNA_DECOMP_EXT_UP, #VRNA_DECOMP_EXT_STEM,
 *      #VRNA_DECOMP_EXT_EXT_EXT, #VRNA_DECOMP_EXT_STEM_EXT, #VRNA_DECOMP_EXT_EXT_STEM,
 *      #VRNA_DECOMP_EXT_EXT_STEM1, vrna_hc_add_f(), vrna_hc_add_data()
 *
 * @param i         Left (5') delimiter position of substructure
 * @param j         Right (3') delimiter position of substructure
 * @param k         Left delimiter of decomposition
 * @param l         Right delimiter of decomposition
 * @param d         Decomposition step indicator
 * @param data      Auxiliary data
 * @return          Pseudo energy contribution in deka-kalories per mol
 */
typedef char (vrna_callback_hc_evaluate)(int i, int j, int k, int l, char d, void *data);

/**
 *  @brief  do not print the header information line
 *  @deprecated   This mode is not supported anymore!
 *
 */
#define VRNA_CONSTRAINT_NO_HEADER         0

/**
 *  @brief  Flag for vrna_constraints_add() to indicate that constraint is passed in pseudo dot-bracket notation
 *
 *  @see vrna_constraints_add(), vrna_message_constraint_options(), vrna_message_constraint_options_all()
 *
 *  @ingroup  hard_constraints
 *
 */
#define VRNA_CONSTRAINT_DB                16384U

/**
 *  @brief Switch for dot-bracket structure constraint to enforce base pairs
 *
 *  This flag should be used to really enforce base pairs given in dot-bracket constraint rather than
 *  just weakly-enforcing them.
 *
 *  @ingroup  hard_constraints
 *
 *  @see  vrna_hc_add_from_db(), vrna_constraints_add(), vrna_message_constraint_options(),
 *        vrna_message_constraint_options_all()
 */
#define VRNA_CONSTRAINT_DB_ENFORCE_BP           32768U

/**
 *  @brief  Flag that is used to indicate the pipe '|' sign in pseudo dot-bracket
 *  notation of hard constraints.
 *
 *  Use this definition to indicate the pipe sign '|' (paired with another base)
 *
 *  @ingroup  hard_constraints
 *
 *  @see  vrna_hc_add_from_db(), vrna_constraints_add(), vrna_message_constraint_options(),
 *        vrna_message_constraint_options_all()
 */
#define VRNA_CONSTRAINT_DB_PIPE              65536U

/**
 *  @brief  dot '.' switch for structure constraints (no constraint at all)
 *
 *  @ingroup  hard_constraints
 *
 *  @see  vrna_hc_add_from_db(), vrna_constraints_add(), vrna_message_constraint_options(),
 *        vrna_message_constraint_options_all()
 */
#define VRNA_CONSTRAINT_DB_DOT               131072U
/**
 *  @brief  'x' switch for structure constraint (base must not pair)
 *
 *  @ingroup  hard_constraints
 *
 *  @see  vrna_hc_add_from_db(), vrna_constraints_add(), vrna_message_constraint_options(),
 *        vrna_message_constraint_options_all()
 */
#define VRNA_CONSTRAINT_DB_X                 262144U
/**
 *  @brief  angle brackets '<', '>' switch for structure constraint (paired downstream/upstream)
 *
 *  @see  vrna_hc_add_from_db(), vrna_constraints_add(), vrna_message_constraint_options(),
 *        vrna_message_constraint_options_all()
 */
#define VRNA_CONSTRAINT_DB_ANG_BRACK         524288U
/**
 *  @brief  round brackets '(',')' switch for structure constraint (base i pairs base j)
 *
 *  @ingroup  hard_constraints
 *
 *  @see  vrna_hc_add_from_db(), vrna_constraints_add(), vrna_message_constraint_options(),
 *        vrna_message_constraint_options_all()
 */
#define VRNA_CONSTRAINT_DB_RND_BRACK         1048576U

/**
 *  @brief  Flag that is used to indicate the character 'l' in pseudo dot-bracket
 *  notation of hard constraints.
 *
 *  Use this definition to indicate the usage of 'l' character (intramolecular pairs only)
 *
 *  @ingroup  hard_constraints
 *
 *  @see  vrna_hc_add_from_db(), vrna_constraints_add(), vrna_message_constraint_options(),
 *        vrna_message_constraint_options_all()
 */
#define VRNA_CONSTRAINT_DB_INTRAMOL    2097152U

/**
 *  @brief  Flag that is used to indicate the character 'e' in pseudo dot-bracket
 *  notation of hard constraints.
 *
 *  Use this definition to indicate the usage of 'e' character (intermolecular pairs only)
 *
 *  @ingroup  hard_constraints
 *
 *  @see  vrna_hc_add_from_db(), vrna_constraints_add(), vrna_message_constraint_options(),
 *        vrna_message_constraint_options_all()
 */
#define VRNA_CONSTRAINT_DB_INTERMOL    4194304U

/**
 *  @brief '+' switch for structure constraint (base is involved in a gquad)
 *
 *  @ingroup  hard_constraints
 *
 *  @see  vrna_hc_add_from_db(), vrna_constraints_add(), vrna_message_constraint_options(),
 *        vrna_message_constraint_options_all()
 *
 *  @warning  This flag is for future purposes only! No implementation recognizes it yet.
 */
#define VRNA_CONSTRAINT_DB_GQUAD                8388608U

/**
 *  @brief Switch for dot-bracket structure constraint with default symbols
 *
 *  This flag conveniently combines all possible symbols in dot-bracket notation
 *  for hard constraints and #VRNA_CONSTRAINT_DB
 *
 *  @ingroup  hard_constraints
 *
 *  @see  vrna_hc_add_from_db(), vrna_constraints_add(), vrna_message_constraint_options(),
 *        vrna_message_constraint_options_all()
 */
#define VRNA_CONSTRAINT_DB_DEFAULT \
    (   VRNA_CONSTRAINT_DB \
      | VRNA_CONSTRAINT_DB_PIPE \
      | VRNA_CONSTRAINT_DB_DOT \
      | VRNA_CONSTRAINT_DB_X \
      | VRNA_CONSTRAINT_DB_ANG_BRACK \
      | VRNA_CONSTRAINT_DB_RND_BRACK \
      | VRNA_CONSTRAINT_DB_INTRAMOL \
      | VRNA_CONSTRAINT_DB_INTERMOL \
      | VRNA_CONSTRAINT_DB_GQUAD \
    )

/**
 *  @brief  Hard constraints flag, base pair in the exterior loop
 *
 *  @ingroup  hard_constraints
 *
 */
#define VRNA_CONSTRAINT_CONTEXT_EXT_LOOP      (char)0x01

/**
 *  @brief  Hard constraints flag, base pair encloses hairpin loop
 *
 *  @ingroup  hard_constraints
 *
 */
#define VRNA_CONSTRAINT_CONTEXT_HP_LOOP       (char)0x02

/**
 *  @brief  Hard constraints flag, base pair encloses an interior loop
 *
 *  @ingroup  hard_constraints
 *
 */
#define VRNA_CONSTRAINT_CONTEXT_INT_LOOP      (char)0x04

/**
 *  @brief  Hard constraints flag, base pair encloses a multi branch loop
 *
 *  @ingroup  hard_constraints
 *
 */
#define VRNA_CONSTRAINT_CONTEXT_INT_LOOP_ENC  (char)0x08

/**
 *  @brief  Hard constraints flag, base pair is enclosed in an interior loop
 *
 *  @ingroup  hard_constraints
 *
 */
#define VRNA_CONSTRAINT_CONTEXT_MB_LOOP       (char)0x10

/**
 *  @brief  Hard constraints flag, base pair is enclosed in a multi branch loop
 *
 *  @ingroup  hard_constraints
 *
 */
#define VRNA_CONSTRAINT_CONTEXT_MB_LOOP_ENC   (char)0x20

/**
 *  @brief  Hard constraint flag to indicate enforcement of constraints
 */
#define VRNA_CONSTRAINT_CONTEXT_ENFORCE       (char)0x40

/**
 *  @brief  Hard constraint flag to indicate not to remove base pairs that conflict with a given constraint
 */
#define VRNA_CONSTRAINT_CONTEXT_NO_REMOVE     (char)0x80

/**
 * @brief  Hard constraints flag, shortcut for all base pairs
 *
 *  @ingroup  hard_constraints
 *
 */
#define VRNA_CONSTRAINT_CONTEXT_ALL_LOOPS     (char)(   VRNA_CONSTRAINT_CONTEXT_EXT_LOOP \
                                                      | VRNA_CONSTRAINT_CONTEXT_HP_LOOP \
                                                      | VRNA_CONSTRAINT_CONTEXT_INT_LOOP \
                                                      | VRNA_CONSTRAINT_CONTEXT_INT_LOOP_ENC \
                                                      | VRNA_CONSTRAINT_CONTEXT_MB_LOOP \
                                                      | VRNA_CONSTRAINT_CONTEXT_MB_LOOP_ENC)

/**
 *  @brief  The hard constraints data structure
 *
 *  The content of this data structure determines the decomposition pattern
 *  used in the folding recursions. Attribute 'matrix' is used as source for
 *  the branching pattern of the decompositions during all folding recursions.
 *  Any entry in matrix[i,j] consists of the 6 LSB that allows one to distinguish the
 *  following types of base pairs:
 *  - in the exterior loop (#VRNA_CONSTRAINT_CONTEXT_EXT_LOOP)
 *  - enclosing a hairpin (#VRNA_CONSTRAINT_CONTEXT_HP_LOOP)
 *  - enclosing an interior loop (#VRNA_CONSTRAINT_CONTEXT_INT_LOOP)  
 *  - enclosed by an exterior loop (#VRNA_CONSTRAINT_CONTEXT_INT_LOOP_ENC)
 *  - enclosing a multi branch loop (#VRNA_CONSTRAINT_CONTEXT_MB_LOOP)
 *  - enclosed by a multi branch loop (#VRNA_CONSTRAINT_CONTEXT_MB_LOOP_ENC)
 *
 *  The four linear arrays 'up_xxx' provide the number of available unpaired
 *  nucleotides (including position i) 3' of each position in the sequence.
 *
 *  @see  vrna_hc_init(), vrna_hc_free(), #VRNA_CONSTRAINT_CONTEXT_EXT_LOOP,
 *        #VRNA_CONSTRAINT_CONTEXT_HP_LOOP, #VRNA_CONSTRAINT_CONTEXT_INT_LOOP,
 *        #VRNA_CONSTRAINT_CONTEXT_MB_LOOP, #VRNA_CONSTRAINT_CONTEXT_MB_LOOP_ENC
 *
 *  @ingroup hard_constraints
 */
struct vrna_hc_s {
  char    *matrix;  /**<  @brief  Upper triangular matrix that encodes where a
                                  base pair or unpaired nucleotide is allowed
                    */
  int     *up_ext;  /**<  @brief  A linear array that holds the number of allowed
                                  unpaired nucleotides in an exterior loop
                    */
  int     *up_hp;   /**<  @brief  A linear array that holds the number of allowed
                                  unpaired nucleotides in a hairpin loop
                    */
  int     *up_int;  /**<  @brief  A linear array that holds the number of allowed
                                  unpaired nucleotides in an interior loop
                    */
  int     *up_ml;   /**<  @brief  A linear array that holds the number of allowed
                                  unpaired nucleotides in a multi branched loop
                    */

  vrna_callback_hc_evaluate *f; /**<  @brief  A function pointer that returns whether or
                                              not a certain decomposition may be evaluated
                                */

  void *data;                         /**<  @brief  A pointer to some structure where the user
                                                    may store necessary data to evaluate its
                                                    generic hard constraint function
                                      */

  vrna_callback_free_auxdata *free_data;  /**<  @brief  A pointer to a function to free memory
                                                        occupied by auxiliary data

                                                The function this pointer is pointing to will be
                                                called upon destruction of the #vrna_hc_s, and
                                                provided with the vrna_hc_s.data pointer that
                                                may hold auxiliary data. Hence, to avoid leaking
                                                memory, the user may use this pointer to free
                                                memory occupied by auxiliary data.
                                          */
};

/**
 *  @brief  A single hard constraint for a single nucleotide
 *
 *  @ingroup hard_constraints
 */
struct vrna_hc_up_s {
  int   position; /**<  @brief The sequence position (1-based)  */
  char  options;  /**<  @brief The hard constraint option       */
};

/**
 *  @brief Print a help message for pseudo dot-bracket structure constraint characters to stdout.
 *  (constraint support is specified by option parameter)
 *
 *  Currently available options are:\n
 *  #VRNA_CONSTRAINT_DB_PIPE (paired with another base)\n
 *  #VRNA_CONSTRAINT_DB_DOT (no constraint at all)\n
 *  #VRNA_CONSTRAINT_DB_X (base must not pair)\n
 *  #VRNA_CONSTRAINT_DB_ANG_BRACK (paired downstream/upstream)\n
 *  #VRNA_CONSTRAINT_DB_RND_BRACK (base i pairs base j)\n
 *
 *  pass a collection of options as one value like this:
 *  @verbatim vrna_message_constraints(option_1 | option_2 | option_n) @endverbatim
 *
 *  @ingroup  constraints
 *
 *  @see  vrna_message_constraint_options_all(), vrna_constraints_add(), #VRNA_CONSTRAINT_DB,
 *        #VRNA_CONSTRAINT_DB_PIPE, #VRNA_CONSTRAINT_DB_DOT, #VRNA_CONSTRAINT_DB_X, #VRNA_CONSTRAINT_DB_ANG_BRACK,
 *        #VRNA_CONSTRAINT_DB_RND_BRACK, #VRNA_CONSTRAINT_DB_INTERMOL, #VRNA_CONSTRAINT_DB_INTRAMOL
 *
 *  @param option Option switch that tells which constraint help will be printed
 */
void vrna_message_constraint_options(unsigned int option);

/**
 *  @brief Print structure constraint characters to stdout
 *  (full constraint support)
 *
 *  @ingroup  constraints
 *
 *  @see  vrna_message_constraint_options(), vrna_constraints_add(), #VRNA_CONSTRAINT_DB,
 *        #VRNA_CONSTRAINT_DB_PIPE, #VRNA_CONSTRAINT_DB_DOT, #VRNA_CONSTRAINT_DB_X, #VRNA_CONSTRAINT_DB_ANG_BRACK,
 *        #VRNA_CONSTRAINT_DB_RND_BRACK, #VRNA_CONSTRAINT_DB_INTERMOL, #VRNA_CONSTRAINT_DB_INTRAMOL
 */
void vrna_message_constraint_options_all(void);

/**
 *  @brief  Initialize/Reset hard constraints to default values
 *
 *  This function resets the hard constraints to their default values, i.e.
 *  all positions may be unpaired in all contexts, and base pairs are
 *  allowed in all contexts, if they resemble canonical pairs.
 *  Previously set hard constraints will be removed vefore initialization.
 *
 *  @ingroup  hard_constraints
 *
 *  @see  vrna_hc_add_bp(), vrna_hc_add_bp_nonspecific(), vrna_hc_add_up()
 *
 *  @param  vc  The fold compound
 */
void vrna_hc_init(vrna_fold_compound_t *vc);

/**
 *  @brief  Make a certain nucleotide unpaired
 *
 *  @ingroup  hard_constraints
 *
 *  @see  vrna_hc_add_bp(), vrna_hc_add_bp_nonspecific(), vrna_hc_init(),
 *        #VRNA_CONSTRAINT_CONTEXT_EXT_LOOP, #VRNA_CONSTRAINT_CONTEXT_HP_LOOP,
 *        #VRNA_CONSTRAINT_CONTEXT_INT_LOOP, #VRNA_CONSTRAINT_CONTEXT_MB_LOOP,
 *        #VRNA_CONSTRAINT_CONTEXT_ALL_LOOPS
 *
 *  @param  vc      The #vrna_fold_compound_t the hard constraints are associated with
 *  @param  i       The position that needs to stay unpaired (1-based)
 *  @param  option  The options flag indicating how/where to store the hard constraints
 */
void vrna_hc_add_up(vrna_fold_compound_t *vc,
                    int i,
                    char option);

/**
 *  @brief Apply a list of hard constraints for single nucleotides
 *
 *  @ingroup  hard_constraints
 *
 *  @param  vc          The #vrna_fold_compound_t the hard constraints are associated with
 *  @param  constraints The list off constraints to apply, last entry must have position
 *                      attribute set to 0
 */
int
vrna_hc_add_up_batch( vrna_fold_compound_t *vc,
                      vrna_hc_up_t *constraints);

/**
 *  @brief  Favorize/Enforce  a certain base pair (i,j)
 *
 *  @ingroup  hard_constraints
 *
 *  @see  vrna_hc_add_bp_nonspecific(), vrna_hc_add_up(), vrna_hc_init(),
 *        #VRNA_CONSTRAINT_CONTEXT_EXT_LOOP, #VRNA_CONSTRAINT_CONTEXT_HP_LOOP,
 *        #VRNA_CONSTRAINT_CONTEXT_INT_LOOP, #VRNA_CONSTRAINT_CONTEXT_INT_LOOP_ENC,
 *        #VRNA_CONSTRAINT_CONTEXT_MB_LOOP, #VRNA_CONSTRAINT_CONTEXT_MB_LOOP_ENC,
 *        #VRNA_CONSTRAINT_CONTEXT_ENFORCE, #VRNA_CONSTRAINT_CONTEXT_ALL_LOOPS
 *
 *  @param  vc      The #vrna_fold_compound_t the hard constraints are associated with
 *  @param  i       The 5' located nucleotide position of the base pair (1-based)
 *  @param  j       The 3' located nucleotide position of the base pair (1-based)
 *  @param  option  The options flag indicating how/where to store the hard constraints
 */
void vrna_hc_add_bp(vrna_fold_compound_t *vc,
                    int i,
                    int j,
                    char option);

/**
 *  @brief  Enforce a nucleotide to be paired (upstream/downstream)
 *
 *  @ingroup  hard_constraints
 *
 *  @see  vrna_hc_add_bp(), vrna_hc_add_up(), vrna_hc_init(),
 *        #VRNA_CONSTRAINT_CONTEXT_EXT_LOOP, #VRNA_CONSTRAINT_CONTEXT_HP_LOOP,
 *        #VRNA_CONSTRAINT_CONTEXT_INT_LOOP, #VRNA_CONSTRAINT_CONTEXT_INT_LOOP_ENC,
 *        #VRNA_CONSTRAINT_CONTEXT_MB_LOOP, #VRNA_CONSTRAINT_CONTEXT_MB_LOOP_ENC,
 *        #VRNA_CONSTRAINT_CONTEXT_ALL_LOOPS
 *
 *  @param  vc      The #vrna_fold_compound_t the hard constraints are associated with
 *  @param  i       The position that needs to stay unpaired (1-based)
 *  @param  d       The direction of base pairing (@f$ d < 0 @f$: pairs upstream,
 *                  @f$ d > 0 @f$: pairs downstream, @f$ d == 0 @f$: no direction)
 *  @param  option  The options flag indicating in which loop type context the pairs may appear
 */
void vrna_hc_add_bp_nonspecific(vrna_fold_compound_t *vc,
                                int i,
                                int d,
                                char option);

/**
 *  @brief  Free the memory allocated by a #vrna_hc_t data structure
 *
 *  Use this function to free all memory that was allocated for a data structure
 *  of type #vrna_hc_t .
 *
 *  @see get_hard_constraints(), #vrna_hc_t
 *
 *  @ingroup  hard_constraints
 *
 */
void vrna_hc_free(vrna_hc_t *hc);


/**
 *  @brief  Add a function pointer pointer for the generic hard constraint
 *          feature
 */
void vrna_hc_add_f( vrna_fold_compound_t *vc,
                    vrna_callback_hc_evaluate *f);

/**
 *  @brief Add an auxiliary data structure for the generic hard constraints callback function
 *
 *  @ingroup generic_hc
 *
 *  @see vrna_hc_add_f()
 *
 *  @param  vc        The fold compound the generic hard constraint function should be bound to
 *  @param  data      A pointer to the data structure that holds required data for function 'f'
 *  @param  free_data A pointer to a function that free's the memory occupied by @data (Maybe NULL)
 */
void vrna_hc_add_data(vrna_fold_compound_t *vc,
                      void *data,
                      vrna_callback_free_auxdata *f);


/**
 *  @brief Add hard constraints from pseudo dot-bracket notation
 *
 *  This function allows one to apply hard constraints from a pseudo dot-bracket
 *  notation. The @p options parameter controls, which characters are recognized
 *  by the parser. Use the #VRNA_CONSTRAINT_DB_DEFAULT convenience macro, if you
 *  want to allow all known characters
 *
 *  @ingroup  hard_constraints
 *
 *  @see  #VRNA_CONSTRAINT_DB_PIPE, #VRNA_CONSTRAINT_DB_DOT, #VRNA_CONSTRAINT_DB_X,
 *        #VRNA_CONSTRAINT_DB_ANG_BRACK, #VRNA_CONSTRAINT_DB_RND_BRACK, #VRNA_CONSTRAINT_DB_INTRAMOL,
 *        #VRNA_CONSTRAINT_DB_INTERMOL, #VRNA_CONSTRAINT_DB_GQUAD
 *
 *  @param  vc            The fold compound
 *  @param  constraint    A pseudo dot-bracket notation of the hard constraint.
 *  @param  options       The option flags
 */
int
vrna_hc_add_from_db(vrna_fold_compound_t *vc,
                    const char *constraint,
                    unsigned int options);


#ifdef  VRNA_BACKWARD_COMPAT

/**
 *  @brief Print structure constraint characters to stdout.
 *  (constraint support is specified by option parameter)
 *
 *  @deprecated Use vrna_message_constraints() instead!
 *  @param option Option switch that tells which constraint help will be printed
 */
DEPRECATED(void print_tty_constraint(unsigned int option));

/**
 *  @brief Print structure constraint characters to stdout
 *  (full constraint support)
 *
 *  @deprecated Use vrna_message_constraint_options_all() instead!
 */
DEPRECATED(void print_tty_constraint_full(void));

/**
 *  @brief Insert constraining pair types according to constraint structure string
 *
 *  @deprecated   Do not use this function anymore! Structure constraints are now handled through #vrna_hc_t and related functions.
 *
 *  @param constraint     The structure constraint string
 *  @param length         The actual length of the sequence (constraint may be shorter)
 *  @param ptype          A pointer to the basepair type array
 *  @param BP             (not used anymore)
 *  @param min_loop_size  The minimal loop size (usually #TURN )
 *  @param idx_type       Define the access type for base pair type array (0 = indx, 1 = iindx)
 */
DEPRECATED(void constrain_ptypes(const char *constraint, unsigned int length, char *ptype, int *BP, int min_loop_size, unsigned int idx_type));

#endif

#endif
