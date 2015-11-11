#ifndef VIENNA_RNA_PACKAGE_CONSTRAINTS_H
#define VIENNA_RNA_PACKAGE_CONSTRAINTS_H

#ifdef __GNUC__
#define DEPRECATED(func) func __attribute__ ((deprecated))
#else
#define DEPRECATED(func) func
#endif


/* make this interface backward compatible with RNAlib < 2.2.0 */
#define VRNA_BACKWARD_COMPAT

/**
 *  @addtogroup constraints
 *  @ingroup folding_routines
 *  @brief This module covers all functions and variables related to the
 *  problem of incorporating secondary structure constraints into the folding
 *  recursions.
 *
 *  This module provides general functions that allow for an easy control of
 *  constrained secondary structure prediction and evaluation.
 *  Secondary Structure constraints can be subdivided into two groups:
 *
 *  - @ref hard_constraints, and
 *  - @ref soft_constraints.
 *
 *  While Hard-Constraints directly influence the production rules used in the folding
 *  recursions by allowing, disallowing, or enforcing certain decomposition steps,
 *  Soft-constraints on the other hand are used to change position specific contributions
 *  in the recursions by adding bonuses/penalties in form of pseudo free energies
 *  to certain loop configurations.
 *
 *  @file constraints.h
 */

/** @brief Typename for the hard constraints data structure #vrna_hc_s
 *  @ingroup  constraints
 */
typedef struct  vrna_hc_s vrna_hc_t;

/** @brief Typename for the soft constraints data structure #vrna_sc_s
 *  @ingroup  constraints
 */
typedef struct  vrna_sc_s vrna_sc_t;

/**
 * @brief Callback to evaluate whether or not a particular decomposition step is contributing to the solution space
 *
 * @ingroup hard_constraints
 *
 * @param i         Left (5') delimiter position of substructure
 * @param j         Right (3') delimiter position of substructure
 * @param k
 * @param l
 * @param d         Decomposition step indicator
 * @param data      Auxiliary data
 * @return          Pseudo energy contribution in deka-kalories per mol
 */
typedef char (vrna_callback_hc_evaluate)(int i, int j, int k, int l, char d, void *data);

/**
 * @brief Callback to retrieve pseudo energy contribution for soft constraint feature
 *
 * @ingroup generalized_sc
 *
 * @param i         Left (5') delimiter position of substructure
 * @param j         Right (3') delimiter position of substructure
 * @param k
 * @param l
 * @param d         Decomposition step indicator
 * @param data      Auxiliary data
 * @return          Pseudo energy contribution in deka-kalories per mol
 */
typedef int (vrna_callback_sc_energy)(int i, int j, int k, int l, char d, void *data);

/**
 * @brief Callback to retrieve pseudo energy contribution as Boltzmann Factors for soft constraint feature
 *
 * @ingroup generalized_sc
 *
 * @param i         Left (5') delimiter position of substructure
 * @param j         Right (3') delimiter position of substructure
 * @param k
 * @param l
 * @param d         Decomposition step indicator
 * @param data      Auxiliary data
 * @return          Pseudo energy contribution in deka-kalories per mol
 */
typedef FLT_OR_DBL (vrna_callback_sc_exp_energy)(int i, int j, int k, int l, char d, void *data);

/**
 * @brief Callback to retrieve auxiliary base pairs for soft constraint feature
 *
 * @ingroup generalized_sc
 *
 * @param i         Left (5') delimiter position of substructure
 * @param j         Right (3') delimiter position of substructure
 * @param k
 * @param l
 * @param d         Decomposition step indicator
 * @param data      Auxiliary data
 * @return          List of additional base pairs
 */
typedef vrna_basepair_t *(vrna_callback_sc_backtrack)(int i, int j, int k, int l, char d, void *data);

#include <ViennaRNA/data_structures.h>

/**
 *  @brief  Flag that is used to indicate the pipe '|' sign in pseudo dot-bracket
 *  notation of hard constraints.
 *
 *  Use this definition to indicate the pipe sign '|' (paired with another base)
 *
 *  @ingroup  constraints
 *
 *  @see vrna_constraints_add(), vrna_message_constraint_options(), vrna_message_constraint_options_all()
 *
 */
#define VRNA_CONSTRAINT_DB_PIPE              1U

/**
 *  @brief  dot '.' switch for structure constraints (no constraint at all)
 *
 *  @ingroup  constraints
 *
 *  @see vrna_constraints_add(), vrna_message_constraint_options(), vrna_message_constraint_options_all()
 *
 */
#define VRNA_CONSTRAINT_DB_DOT               2U
/**
 *  @brief  'x' switch for structure constraint (base must not pair)
 *
 *  @ingroup  constraints
 *
 *  @see vrna_constraints_add(), vrna_message_constraint_options(), vrna_message_constraint_options_all()
 *
 */
#define VRNA_CONSTRAINT_DB_X                 4U
/**
 *  @brief  angle brackets '<', '>' switch for structure constraint (paired downstream/upstream)
 *
 *  @ingroup  constraints
 *
 *  @see vrna_constraints_add(), vrna_message_constraint_options(), vrna_message_constraint_options_all()
 *
 */
#define VRNA_CONSTRAINT_DB_ANG_BRACK         8U
/**
 *  @brief  round brackets '(',')' switch for structure constraint (base i pairs base j)
 *
 *  @ingroup  constraints
 *
 *  @see vrna_constraints_add(), vrna_message_constraint_options(), vrna_message_constraint_options_all()
 *
 */
#define VRNA_CONSTRAINT_DB_RND_BRACK         16U

/**
 *  @brief  Flag that is used to indicate the character 'l' in pseudo dot-bracket
 *  notation of hard constraints.
 *
 *  Use this definition to indicate the usage of 'l' character (intramolecular pairs only)
 *
 *  @ingroup  constraints
 *
 *  @see vrna_constraints_add(), vrna_message_constraint_options(), vrna_message_constraint_options_all()
 *
 */
#define VRNA_CONSTRAINT_DB_INTRAMOL    2048U

/**
 *  @brief  Flag that is used to indicate the character 'e' in pseudo dot-bracket
 *  notation of hard constraints.
 *
 *  Use this definition to indicate the usage of 'e' character (intermolecular pairs only)
 *
 *  @ingroup  constraints
 *
 *  @see vrna_constraints_add(), vrna_message_constraint_options(), vrna_message_constraint_options_all()
 *
 */
#define VRNA_CONSTRAINT_DB_INTERMOL    4096U

/**
 *  @brief '+' switch for structure constraint (base is involved in a gquad)
 *
 *  @ingroup  constraints
 *
 *  @see vrna_constraints_add(), vrna_message_constraint_options(), vrna_message_constraint_options_all()
 *  @warning  This flag is for future purposes only! No implementation recognizes it yet.
 */
#define VRNA_CONSTRAINT_DB_GQUAD                8192U

/**
 *  @brief Switch for dot-bracket structure constraint to enforce base pairs
 *
 *  This flag should be used to really enforce base pairs given in dot-bracket constraint rather than
 *  just weakly-enforcing them.
 *
 *  @ingroup  constraints
 *
 *  @see vrna_constraints_add()
 */
#define VRNA_CONSTRAINT_DB_ENFORCE_BP           16384U

/**
 *  @brief  constraint may span over several lines
 *
 *  @ingroup  hard_constraints
 *
 */
#define VRNA_CONSTRAINT_MULTILINE         32U
/**
 *  @brief  do not print the header information line
 *
 *  @ingroup  hard_constraints
 *
 */
#define VRNA_CONSTRAINT_NO_HEADER         64U
/**
 *  @brief  placeholder for all constraining characters
 *
 *  @ingroup  constraints
 *
 */
#define VRNA_CONSTRAINT_ALL               128U

/**
 *  @brief  Flag for vrna_constraints_add() to indicate that constraint is passed in pseudo dot-bracket notation
 *
 *  @see vrna_constraints_add(), vrna_message_constraint_options(), vrna_message_constraint_options_all()
 *
 *  @ingroup  constraints
 *
 */
#define VRNA_CONSTRAINT_DB                256U

/**
 *  @brief  Flag for vrna_constraints_add() to indicate that constraints are present in a text file
 *
 *  @see vrna_constraints_add()
 *
 *  @ingroup  constraints
 *
 */
#define VRNA_CONSTRAINT_FILE              512U

/**
 *  @brief  Soft constraints flag, apply constraints for MFE calculations
 *
 *  @ingroup  soft_constraints
 *
 */
#define VRNA_CONSTRAINT_SOFT_MFE          8192U

/**
 *  @brief  Soft constraints flag, apply constraints for partition function calculations
 *
 *  @ingroup  soft_constraints
 *
 */
#define VRNA_CONSTRAINT_SOFT_PF           16384U


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

#define VRNA_CONSTRAINT_CONTEXT_ENFORCE       (char)0x40

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
 *  @brief  Generalized constraint folding flag indicating hairpin loop decomposition step
 *
 *  @ingroup  generalized_sc
 *
 */
#define VRNA_DECOMP_PAIR_HP     1

/**
 *  @brief  Generalized constraint folding flag indicating interior loop decomposition step
 *
 *  @ingroup  generalized_sc
 *
 */
#define VRNA_DECOMP_PAIR_IL     2

/**
 *  @ingroup  generalized_sc
 *
 */
#define VRNA_DECOMP_PAIR_ML     3

/**
 *  @ingroup  generalized_sc
 *
 */
#define VRNA_DECOMP_ML_ML_ML    5

/**
 *  @ingroup  generalized_sc
 *
 */
#define VRNA_DECOMP_ML_STEM     4

/**
 *  @ingroup  generalized_sc
 *
 */
#define VRNA_DECOMP_ML_ML       6

/**
 *  @ingroup  generalized_sc
 *
 */
#define VRNA_DECOMP_ML_UP       11

/**
 *  @ingroup  generalized_sc
 *
 */
#define VRNA_DECOMP_ML_COAXIAL  13

/**
 *  @ingroup  generalized_sc
 *
 */
#define VRNA_DECOMP_EXT_EXT     9

/**
 *  @ingroup  generalized_sc
 *
 */
#define VRNA_DECOMP_EXT_UP_3    7

/**
 *  @ingroup  generalized_sc
 *
 */
#define VRNA_DECOMP_EXT_UP      8

/**
 *  @ingroup  generalized_sc
 *
 */
#define VRNA_DECOMP_EXT_STEM 14

/**
 *  @ingroup  generalized_sc
 *
 */
#define VRNA_DECOMP_EXT_EXT_EXT 15

/**
 *  @ingroup  generalized_sc
 *
 */
#define VRNA_DECOMP_EXT_STEM_EXT 16

/**
 *  @ingroup  generalized_sc
 *
 */
#define VRNA_DECOMP_EXT_STEM_OUTSIDE 17

/**
 *  @ingroup  generalized_sc
 *
 */
#define VRNA_DECOMP_EXT_EXT_STEM 18

/**
 *  @ingroup  generalized_sc
 *
 */
#define VRNA_DECOMP_EXT_EXT_STEM1 19

/**
 *  @ingroup  generalized_sc
 *
 */
#define VRNA_DECOMP_ML_ML_STEM 20

/**
 *  @brief  The hard constraints data structure
 *
 *  The content of this data structure determines the decomposition pattern
 *  used in the folding recursions. Attribute 'matrix' is used as source for
 *  the branching pattern of the decompositions during all folding recursions.
 *  Any entry in matrix[i,j] consists of the 6 LSB that allows to distinguish the
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
 *        #VRNA_CONSTRAINT_CONTEXT_EXT_LOOP_ENC, #VRNA_CONSTRAINT_CONTEXT_MB_LOOP,
 *        #VRNA_CONSTRAINT_CONTEXT_MB_LOOP_ENC
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
                                                    generalized hard constraint function
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
 *  @brief  The soft constraints data structure
 *
 *  @ingroup soft_constraints
 */
struct vrna_sc_s {
  int         **energy_up;            /**<  @brief Energy contribution for stretches of unpaired nucleotides */
  int         *energy_bp;             /**<  @brief Energy contribution for base pairs */
  FLT_OR_DBL  **exp_energy_up;        /**<  @brief Boltzmann Factors of the energy contributions for unpaired sequence stretches */
  FLT_OR_DBL  *exp_energy_bp;         /**<  @brief Boltzmann Factors of the energy contribution for base pairs */

  int         *energy_stack;          /**<  @brief Pseudo Energy contribution per base pair involved in a stack */
  FLT_OR_DBL  *exp_energy_stack;      /**<  @brief Boltzmann weighted pseudo energy contribution per nucleotide involved in a stack */

  /* generalized soft contraints below */
  vrna_callback_sc_energy *f;         /**<  @brief  A function pointer used for pseudo
                                                    energy contribution in MFE calculations
                                            @see    vrna_sc_add_f()
                                      */

  vrna_callback_sc_backtrack *bt;     /**<  @brief  A function pointer used to obtain backtraced
                                                    base pairs in loop regions that were altered
                                                    by soft constrained pseudo energy contributions
                                            @see    vrna_sc_add_bt()
                                      */

  vrna_callback_sc_exp_energy *exp_f; /**<  @brief  A function pointer used for pseudo energy
                                                    contribution boltzmann factors in PF
                                                    calculations
                                            @see    vrna_sc_add_exp_f()
                                      */

  void *data;                         /**<  @brief  A pointer to the data object provided for
                                                    for pseudo energy contribution functions of the
                                                    generalized soft constraints feature
                                      */
  vrna_callback_free_auxdata *free_data;
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
 *  @brief  Add constraints to a #vrna_fold_compound_t data structure
 *
 *  Use this function to add/update the hard/soft constraints
 *  The function allows for passing a string 'constraint' that can either be a
 *  filename that points to a constraints definition file or it may be a
 *  pseudo dot-bracket notation indicating hard constraints. Depending on
 *  the type of the string the user has to pass #VRNA_CONSTRAINT_FILE or
 *  #VRNA_CONSTRAINT_DB in the option parameter, respectively. If none of these
 *  to options are passed, no action is performed, other than to guarantee that
 *  at least a hard constraints data structure of type #vrna_hc_t with default
 *  values is present in 'vc'. Already existing hard constraints are not touched.
 *
 *  In case, a psuedo dot-bracket string is passed as the second argument, the
 *  user has to specify, which characters are allowed to be interpreted as
 *  constraints by passing the corresponding options via the third parameter.
 *
 *  @see      vrna_hc_init(), vrna_sc_init(), vrna_hc_add_up(), vrna_hc_add_bp(),
 *            vrna_sc_add_up(), vrna_sc_add_bp(), vrna_sc_add_SHAPE_deigan(),
 *            vrna_sc_add_SHAPE_zarringhalam(), vrna_hc_free(), vrna_sc_free(),
 *            #VRNA_CONSTRAINT_FILE, #VRNA_CONSTRAINT_DB, #VRNA_CONSTRAINT_DB_PIPE,
 *            #VRNA_CONSTRAINT_DB_DOT, #VRNA_CONSTRAINT_DB_X, #VRNA_CONSTRAINT_DB_ANG_BRACK,
 *            #VRNA_CONSTRAINT_DB_RND_BRACK, #VRNA_CONSTRAINT_DB_INTRAMOL,
 *            #VRNA_CONSTRAINT_DB_INTERMOL, #VRNA_CONSTRAINT_DB_GQUAD
 *
 *  @ingroup  constraints
 *
 *  @param  vc            The fold compound
 *  @param  constraint    A string with either the filename of the constraint definitions
 *                        or a pseudo dot-bracket notation of the hard constraint. May be NULL.
 *  @param  options       The option flags
 */
void vrna_constraints_add(vrna_fold_compound_t *vc,
                          const char *constraint,
                          unsigned int options);

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

#ifdef WITH_GEN_HC

/**
 *  @brief  Add a function pointer pointer for the generalized hard constraint
 *          feature
 */
void vrna_hc_add_f( vrna_fold_compound_t *vc,
                    vrna_callback_hc_evaluate *f);

void vrna_hc_add_data(vrna_fold_compound_t *vc,
                      void *data,
                      vrna_callback_free_auxdata *f);

#endif

/**
 *  @brief Initialize an empty soft constraints data structure within a #vrna_fold_compound_t
 *
 *  This function adds a proper soft constraints data structure
 *  to the #vrna_fold_compound_t data structure.
 *  If soft constraints already exist within the fold compound, they are removed.
 *
 *  \note Accepts vrna_fold_compound_t of type #VRNA_VC_TYPE_SINGLE and #VRNA_VC_TYPE_ALIGNMENT
 *
 *  @ingroup  soft_constraints
 *
 *  @see  vrna_sc_add_bp(), vrna_sc_add_up(), vrna_sc_add_SHAPE_deigan(),
 *        vrna_sc_add_SHAPE_zarringhalam(), vrna_sc_remove(), vrna_sc_add_f(),
 *        vrna_sc_add_exp_f(), vrna_sc_add_pre(), vrna_sc_add_post()
 *  @param  vc  The #vrna_fold_compound_t where an empty soft constraint feature is to be added to
 */
void vrna_sc_init(vrna_fold_compound_t *vc);

/**
 *  @brief  Add soft constraints for paired nucleotides
 *
 *  @ingroup  soft_constraints
 *
 *  @param  vc          The #vrna_fold_compound_t the soft constraints are associated with
 *  @param  constraints A two-dimensional array of pseudo free energies in @f$ kcal / mol @f$
 *  @param  options     The options flag indicating how/where to store the soft constraints
 */
void vrna_sc_add_bp(vrna_fold_compound_t *vc,
                    const double **constraints,
                    unsigned int options);

/**
 *  @brief  Add soft constraints for unpaired nucleotides
 *
 *  @ingroup  soft_constraints
 *
 *  @param  vc          The #vrna_fold_compound_t the soft constraints are associated with
 *  @param  constraints A vector of pseudo free energies in @f$ kcal / mol @f$
 *  @param  options     The options flag indicating how/where to store the soft constraints
 */
void vrna_sc_add_up(vrna_fold_compound_t *vc,
                    const double *constraints,
                    unsigned int options);

/**
 *  @brief  Remove soft constraints from #vrna_fold_compound_t
 *
 *  \note Accepts vrna_fold_compound_t of type #VRNA_VC_TYPE_SINGLE and #VRNA_VC_TYPE_ALIGNMENT
 *
 *  @ingroup  soft_constraints
 *
 *  @param  vc  The #vrna_fold_compound_t possibly containing soft constraints
 */
void vrna_sc_remove(vrna_fold_compound_t *vc);

/**
 *  @brief  Free memory occupied by a #vrna_sc_t data structure
 *
 *  @ingroup  soft_constraints
 *
 *  @param  sc  The data structure to free from memory
 */
void vrna_sc_free(vrna_sc_t *sc);

/**
 *  @brief  Add SHAPE reactivity data as soft constraints (Deigan et al. method)
 *
 *  This approach of SHAPE directed RNA folding uses the simple linear ansatz
 *  @f[ \Delta G_{\text{SHAPE}}(i) = m \ln(\text{SHAPE reactivity}(i)+1)+ b @f]
 *  to convert SHAPE reactivity values to pseudo energies whenever a
 *  nucleotide @f$ i @f$ contributes to a stacked pair. A positive slope @f$ m @f$
 *  penalizes high reactivities in paired regions, while a negative intercept @f$ b @f$
 *  results in a confirmatory ``bonus'' free energy for correctly predicted base pairs.
 *  Since the energy evaluation of a base pair stack involves two pairs, the pseudo
 *  energies are added for all four contributing nucleotides. Consequently, the
 *  energy term is applied twice for pairs inside a helix and only once for pairs
 *  adjacent to other structures. For all other loop types the energy model remains
 *  unchanged even when the experimental data highly disagrees with a certain motif.
 *
 *  @see  For further details, we refer to @cite deigan:2009.
 *  @see  vrna_sc_remove(), vrna_sc_add_SHAPE_zarringhalam(), vrna_sc_minimize_pertubation()
 *  @ingroup  soft_constraints
 *  @param  vc            The #vrna_fold_compound_t the soft constraints are associated with
 *  @param  reactivities  A vector of normalized SHAPE reactivities
 *  @param  m             The slope of the conversion function
 *  @param  b             The intercept of the conversion function
 *  @param  options       The options flag indicating how/where to store the soft constraints
 *  @return               1 on successful extraction of the method, 0 on errors
 */
int vrna_sc_add_SHAPE_deigan( vrna_fold_compound_t *vc,
                              const double *reactivities,
                              double m,
                              double b,
                              unsigned int options);

/**
 *  @brief  Add SHAPE reactivity data from files as soft constraints for consensus structure prediction (Deigan et al. method)
 *
 *  @ingroup  soft_constraints
 *  @param  vc            The #vrna_fold_compound_t the soft constraints are associated with
 *  @param  shape_files   A set of filenames that contain normalized SHAPE reactivity data
 *  @param  shape_file_association  An array of integers that associate the files with sequences in the alignment
 *  @param  m             The slope of the conversion function
 *  @param  b             The intercept of the conversion function
 *  @param  options       The options flag indicating how/where to store the soft constraints
 *  @return               1 on successful extraction of the method, 0 on errors
 */
int vrna_sc_add_SHAPE_deigan_ali( vrna_fold_compound_t *vc,
                                  const char **shape_files,
                                  const int *shape_file_association,
                                  double m,
                                  double b,
                                  unsigned int options);

/**
 *  @brief  Add SHAPE reactivity data as soft constraints (Zarringhalam et al. method)
 *
 *  This method first converts the observed SHAPE reactivity of nucleotide @f$ i @f$ into a
 *  probability @f$ q_i @f$ that position @f$ i @f$ is unpaired by means of a non-linear map.
 *  Then pseudo-energies of the form @f[ \Delta G_{\text{SHAPE}}(x,i) = \beta\ |x_i - q_i| @f]
 *  are computed, where @f$ x_i=0 @f$ if position @f$ i @f$ is unpaired and @f$ x_i=1 @f$
 *  if @f$ i @f$ is paired in a given secondary structure. The parameter @f$ \beta @f$ serves as
 *  scaling factor. The magnitude of discrepancy between prediction and experimental observation
 *  is represented by @f$ |x_i - q_i| @f$.
 *
 *  @see For further details, we refer to @cite zarringhalam:2012
 *  @see  vrna_sc_remove(), vrna_sc_add_SHAPE_deigan(), vrna_sc_minimize_pertubation()
 *  @ingroup  soft_constraints
 *  @param  vc            The #vrna_fold_compound_t the soft constraints are associated with
 *  @param  reactivities  A vector of normalized SHAPE reactivities
 *  @param  b             The scaling factor @f$ \beta @f$ of the conversion function
 *  @param  options       The options flag indicating how/where to store the soft constraints
 *  @return               1 on successful extraction of the method, 0 on errors
 */
int vrna_sc_add_SHAPE_zarringhalam( vrna_fold_compound_t *vc,
                                    const double *reactivities,
                                    double b,
                                    double default_value,
                                    const char *shape_conversion,
                                    unsigned int options);

/**
 *  @brief  Parse a character string and extract the encoded SHAPE reactivity conversion
 *          method and possibly the parameters for conversion into pseudo free energies
 *
 *  @ingroup soft_cosntraints
 *
 *  @param  method_string   The string that contains the encoded SHAPE reactivity conversion method
 *  @param  method          A pointer to the memory location where the method character will be stored
 *  @param  param_1         A pointer to the memory location where the first parameter of the corresponding method will be stored
 *  @param  param_2         A pointer to the memory location where the second parameter of the corresponding method will be stored
 *  @return                 1 on successful extraction of the method, 0 on errors
 */
int vrna_sc_SHAPE_parse_method( const char *method_string,
                                char *method,
                                float *param_1,
                                float *param_2);

/**
 *  @brief Convert SHAPE reactivity values to probabilities for being unpaired
 *
 *  This function parses the informations from a given file and stores the result
 *  in the preallocated string sequence and the double array values.
 *
 *  @ingroup  soft_constraints
 *
 *  @see vrna_file_SHAPE_read()
 *  @param shape_conversion String definining the method used for the conversion process
 *  @param values           Pointer to an array of SHAPE reactivities
 *  @param length           Length of the array of SHAPE reactivities
 *  @param default_value    Result used for position with invalid/missing reactivity values
 */
int vrna_sc_SHAPE_to_pr(const char *shape_conversion,
                        double *values,
                        int length,
                        double default_value);

void vrna_sc_add_data(vrna_fold_compound_t *vc,
                      void *data,
                      vrna_callback_free_auxdata *free_data);

/**
 *  @brief  Bind a function pointer for generalized soft constraint feature (MFE version)
 *
 *  This function allows to easily bind a function pointer and corresponding data structure
 *  to the soft constraint part #vrna_sc_t of the #vrna_fold_compound_t.
 *  The function for evaluating the generalized soft constraint feature has to return
 *  a pseudo free energy @f$ \hat{E} @f$ in @f$ dacal/mol @f$, where @f$ 1 dacal/mol = 10 cal/mol @f$.
 *
 *  @ingroup generalized_sc
 *
 *  @param  vc    The fold compound the generalized soft constraint function should be bound to
 *  @param  f     A pointer to the function that evaluates the generalized soft constraint feature
 *  @param  data  A pointer to the data structure that holds required data for function 'f'
 */
void vrna_sc_add_f( vrna_fold_compound_t *vc,
                    vrna_callback_sc_energy *f);

/**
 *  @brief  Bind a backtracking function pointer for generalized soft constraint feature
 *
 *  This function allows to easily bind a function pointer to the soft constraint part
 *  #vrna_sc_t of the #vrna_fold_compound_t.
 *  The provided function should be used for backtracking purposes in loop regions
 *  that were altered via the generalized soft constraint feature. It has to return
 *  an array of #vrna_basepair_t data structures, were the last element in the list is indicated
 *  by a value of -1 in it's i position.
 *
 *  @ingroup generalized_sc
 *
 *  @param  vc    The fold compound the generalized soft constraint function should be bound to
 *  @param  f     A pointer to the function that returns additional base pairs
 */
void vrna_sc_add_bt(vrna_fold_compound_t *vc,
                    vrna_callback_sc_backtrack *f);

/**
 *  @brief  Bind a function pointer for generalized soft constraint feature (PF version)
 *
 *  This function allows to easily bind a function pointer and corresponding data structure
 *  to the soft constraint part #vrna_sc_t of the #vrna_fold_compound_t.
 *  The function for evaluating the generalized soft constraint feature has to return
 *  a pseudo free energy @f$ \hat{E} @f$ as Boltzmann factor, i.e. @f$ exp(- \hat{E} / kT) @f$.
 *  The required unit for @f$ E @f$ is @f$ cal/mol @f$.
 *
 *  @ingroup generalized_sc
 *
 *  @param  vc    The fold compound the generalized soft constraint function should be bound to
 *  @param  exp_f A pointer to the function that evaluates the generalized soft constraint feature
 *  @param  data  A pointer to the data structure that holds required data for function 'f'
 */
void vrna_sc_add_exp_f( vrna_fold_compound_t *vc,
                        vrna_callback_sc_exp_energy *exp_f);

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
