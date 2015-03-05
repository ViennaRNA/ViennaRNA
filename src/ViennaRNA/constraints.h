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
 *  @brief This section covers all functions and variables related to the
 *  problem of incorporating secondary structure constraints into the folding
 *  recursions.
 *
 *  Secondary Structure constraints can further be subdivided into so called
 *  Hard-Constraints and Soft-Constraints. Hard-Constraints directly influence
 *  the production rules used in the folding recursions by allowing, disallowing
 *  or enforcing certain decomposition steps. Soft-constraints on the other hand
 *  are used to change position specific contributions in the recursions by
 *  adding a bonus/penalty to certain loop configurations.
 *
 */

#include <ViennaRNA/data_structures.h>

/**
 *  @brief  Flag that is used to indicate the pipe '|' sign in pseudo dot-bracket
 *  notation of hard constraints.
 *
 *  Use this definition to indicate the pipe sign '|' (paired with another base)
 *  
 *  @ingroup  hard_constraints
 *
 */
#define VRNA_CONSTRAINT_PIPE              1U

/**
 *  @brief  dot '.' switch for structure constraints (no constraint at all)
 *  
 *  @ingroup  hard_constraints
 *
 */
#define VRNA_CONSTRAINT_DOT               2U
/**
 *  @brief  'x' switch for structure constraint (base must not pair)
 *  
 *  @ingroup  hard_constraints
 *
 */
#define VRNA_CONSTRAINT_X                 4U
/**
 *  @brief  angle brackets '<', '>' switch for structure constraint (paired downstream/upstream)
 *  
 *  @ingroup  hard_constraints
 *
 */
#define VRNA_CONSTRAINT_ANG_BRACK         8U
/**
 *  @brief  round brackets '(',')' switch for structure constraint (base i pairs base j)
 *  
 *  @ingroup  hard_constraints
 *
 */
#define VRNA_CONSTRAINT_RND_BRACK         16U

/**
 *  @brief  Flag that is used to indicate the character 'l' in pseudo dot-bracket
 *  notation of hard constraints.
 *
 *  Use this definition to indicate the usage of 'l' character (intramolecular pairs only)
 *
 *  @ingroup  hard_constraints
 *
 */
#define VRNA_CONSTRAINT_INTRAMOLECULAR    2048U

/**
 *  @brief  Flag that is used to indicate the character 'e' in pseudo dot-bracket
 *  notation of hard constraints.
 *
 *  Use this definition to indicate the usage of 'e' character (intermolecular pairs only)
 *
 *  @ingroup  hard_constraints
 *
 */
#define VRNA_CONSTRAINT_INTERMOLECULAR    4096U

/**
 *  @brief '+' switch for structure constraint (base is involved in a gquad)
 */
#define VRNA_CONSTRAINT_G                8192U

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
 *  @ingroup  hard_constraints
 *
 */
#define VRNA_CONSTRAINT_ALL               128U

/**
 *  @brief  
 *  
 *  @ingroup  hard_constraints
 *
 */
#define VRNA_CONSTRAINT_DB                256U

/**
 *  @brief  
 *  
 *  @ingroup  hard_constraints
 *
 */
#define VRNA_CONSTRAINT_FILE              512U

/**
 *  @brief  
 *  
 *  @ingroup  hard_constraints
 *
 */
#define VRNA_CONSTRAINT_IINDX             1024U

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
 *  @brief  Soft constraints flag, apply constraints for unpaired nucleotides
 *  
 *  @ingroup  soft_constraints
 *
 */
#define VRNA_CONSTRAINT_SOFT_UP           32768U

/**
 *  @brief  Soft constraints flag, apply constraints for base pairs
 *  
 *  @ingroup  soft_constraints
 *
 */
#define VRNA_CONSTRAINT_SOFT_BP           65536U

/**
 *  @brief  Soft constraints flag, remove gaps from aligned sequence
 *  
 *  @ingroup  soft_constraints
 *
 */
#define VRNA_CONSTRAINT_UNGAP             131072U

/**
 *  @brief  Hard constraints flag, base pair in the exterior loop
 *  
 *  @ingroup  hard_constraints
 *
 */
#define VRNA_HC_CONTEXT_EXT_LOOP       (char)0x01
/**
 *  @brief  Hard constraints flag, base pair encloses hairpin loop
 *  
 *  @ingroup  hard_constraints
 *
 */
#define VRNA_HC_CONTEXT_HP_LOOP        (char)0x02
/**
 *  @brief  Hard constraints flag, base pair encloses an interior loop
 *  
 *  @ingroup  hard_constraints
 *
 */
#define VRNA_HC_CONTEXT_INT_LOOP       (char)0x04
/**
 *  @brief  Hard constraints flag, base pair is enclosed in an interior loop
 *  
 *  @ingroup  hard_constraints
 *
 */
#define VRNA_HC_CONTEXT_MB_LOOP        (char)0x08
/**
 *  @brief  Hard constraints flag, base pair encloses a multi branch loop
 *
 *  @ingroup  hard_constraints
 *
 */
#define VRNA_HC_CONTEXT_INT_LOOP_ENC  (char)0x10
/**
 *  @brief  Hard constraints flag, base pair is enclosed in a multi branch loop
 *  
 *  @ingroup  hard_constraints
 *
 */
#define VRNA_HC_CONTEXT_MB_LOOP_ENC   (char)0x20

#define VRNA_HC_CONTEXT_ALL_LOOPS     VRNA_HC_CONTEXT_EXT_LOOP \
                                      | VRNA_HC_CONTEXT_HP_LOOP \
                                      | VRNA_HC_CONTEXT_INT_LOOP \
                                      | VRNA_HC_CONTEXT_INT_LOOP_ENC \
                                      | VRNA_HC_CONTEXT_MB_LOOP \
                                      | VRNA_HC_CONTEXT_MB_LOOP_ENC

#define VRNA_HC_CONTEXT_ENFORCE       (char)0x40

/**
 *  @brief  Generalized constraint folding flag indicating hairpin loop decomposition step
 *
 *  @ingroup  constraints
 *
 */
#define VRNA_DECOMP_PAIR_HP     1

/**
 *  @brief  Generalized constraint folding flag indicating interior loop decomposition step
 *
 *  @ingroup  constraints
 *
 */
#define VRNA_DECOMP_PAIR_IL     2

/**
 *  @ingroup  constraints
 *
 */
#define VRNA_DECOMP_PAIR_ML     3

/**
 *  @ingroup  constraints
 *
 */
#define VRNA_DECOMP_ML_ML_ML    5

/**
 *  @ingroup  constraints
 *
 */
#define VRNA_DECOMP_ML_UP_3     4

/**
 *  @ingroup  constraints
 *
 */
#define VRNA_DECOMP_ML_UP_5     6

/**
 *  @ingroup  constraints
 *
 */
#define VRNA_DECOMP_ML_UP       11

/**
 *  @ingroup  constraints
 *
 */
#define VRNA_DECOMP_EXT_EXT     9

/**
 *  @ingroup  constraints
 *
 */
#define VRNA_DECOMP_EXT_UP_3    7

/**
 *  @ingroup  constraints
 *
 */
#define VRNA_DECOMP_EXT_UP_5    10

/**
 *  @ingroup  constraints
 *
 */
#define VRNA_DECOMP_EXT_UP      8

/**
 *  @ingroup  constraints
 *
 */
#define VRNA_DECOMP_EXT_STEM_UP 12


/**
 *  @brief    A flag passed to the generalized soft constraints pre-, and post- functions to
 *            indicate Minimum Free Energy (MFE) processing
 *  @details  This flag is passed as second argument to the pre-, and post- processing
 *            funtions that are bound to the #vrna_sc_t structure via vrna_sc_add_pre(), and
 *            vrna_sc_add_post(), respectively.
 *            Use it in your implementation of the pre-, and post-processing functions to
 *            determine the mode of action required for corresponding pre-, and post-
 *            processing of data available to the function.
 *  @note     This flag will be passed by calls of vrna_fold(), vrna_ali_fold(), vrna_cofold(),
 *            and vrna_subopt()
 *  @ingroup  soft_constraints
 */
#define VRNA_SC_GEN_MFE         (char)1

/**
 *  @brief    A flag passed to the generalized soft constraints pre-, and post- functions to
 *            indicate Partition function (PF) processing
 *  @details  This flag is passed as second argument to the pre-, and post- processing
 *            funtions that are bound to the #vrna_sc_t structure via vrna_sc_add_pre(), and
 *            vrna_sc_add_post(), respectively.
 *            Use it in your implementation of the pre-, and post-processing functions to
 *            determine the mode of action required for corresponding pre-, and post-
 *            processing of data available to the function.
 *  @note     This flag will be passed by calls of vrna_pf_fold(), vrna_ali_pf_fold(), and
 *            vrna_co_pf_fold().
 *  @ingroup  soft_constraints
 */
#define VRNA_SC_GEN_PF          (char)2

/**
 *  @brief  The hard constraints data structure
 *
 *  The content of this data structure determines the decomposition pattern
 *  used in the folding recursions. Attribute 'matrix' is used as source for
 *  the branching pattern of the decompositions during all folding recursions.
 *  Any entry in matrix[i,j] consists of the 6 LSB that allows to distinguish the
 *  following types of base pairs:
 *  - in the exterior loop (#VRNA_HC_CONTEXT_EXT_LOOP)
 *  - enclosing a hairpin (#VRNA_HC_CONTEXT_HP_LOOP)
 *  - enclosing an interior loop (#VRNA_HC_CONTEXT_INT_LOOP)
 *  - enclosed by an exterior loop (#VRNA_HC_CONTEXT_INT_LOOP_ENC)
 *  - enclosing a multi branch loop (#VRNA_HC_CONTEXT_MB_LOOP)
 *  - enclosed by a multi branch loop (#VRNA_HC_CONTEXT_MB_LOOP_ENC)
 *
 *  The four linear arrays 'up_xxx' provide the number of available unpaired
 *  nucleotides (including position i) 3' of each position in the sequence.
 *
 *  @see  get_hard_constraints(), vrna_hc_free(), #VRNA_HC_CONTEXT_EXT_LOOP,
 *        #VRNA_HC_CONTEXT_HP_LOOP, #VRNA_HC_CONTEXT_INT_LOOP, #VRNA_HC_CONTEXT_EXT_LOOP_ENC, #VRNA_HC_CONTEXT_MB_LOOP, #VRNA_HC_CONTEXT_MB_LOOP_ENC
 *        
 *  @ingroup hard_constraints
 */
typedef struct vrna_hc_t {
  char    *matrix;  /**<  @brief  Upper triangular matrix encoding where a
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
} vrna_hc_t;

/**
 *  @brief  The soft constraints data structure
 *
 *  @ingroup soft_constraints
 */
typedef struct vrna_sc_t {
  double      *constraints;           /**<  @brief Backup storage for energy contributions of single nucleotides */
  int         **free_energies;        /**<  @brief Energy contribution for unpaired sequence stretches */
  int         *en_basepair;           /**<  @brief Energy contribution for base pairs */
  FLT_OR_DBL  **boltzmann_factors;    /**<  @brief Boltzmann Factors of the energy contributions for unpaired sequence stretches */
  FLT_OR_DBL  *exp_en_basepair;       /**<  @brief Boltzmann Factors of the energy contribution for base pairs */

  int         *en_stack;              /**<  @brief Pseudo Energy contribution per base pair involved in a stack */
  FLT_OR_DBL  *exp_en_stack;          /**<  @brief Boltzmann weighted pseudo energy contribution per nucleotide involved in a stack */

  /* generalized soft contraints below */
  int (*f)( int,
            int,
            int,
            int,
            char,
            void *);                  /**<  @brief  A function pointer used for pseudo
                                                    energy contribution in MFE calculations
                                            @see    vrna_sc_add_f()
                                      */

  FLT_OR_DBL (*exp_f)(int,
                      int,
                      int,
                      int,
                      char,
                      void *);        /**<  @brief  A function pointer used for pseudo energy
                                                    contribution boltzmann factors in PF
                                                    calculations
                                            @see    vrna_sc_add_exp_f()
                                      */

  void (*pre)(vrna_fold_compound *,
              char);                  /**<  @brief    A function pointer to some generalized soft
                                                      constraints preprocessing function.
                                            @details  This function will be called just before
                                                      the forward recursions start
                                      */

  void (*post)( vrna_fold_compound *,
                char);                /**<  @brief    A function pointer to some generalized soft
                                                      constraints postprocessing function.
                                            @details  This function will be called right after the
                                                      forward recursions or backtracking, whatever
                                                      is last.
                                      */

  void *data;                         /**<  @brief  A pointer to the data object provided for
                                                    for pseudo energy contribution functions of the
                                                    generalized soft constraints feature
                                      */
} vrna_sc_t;

/**
 *  @brief Print a help message for pseudo dot-bracket structure constraint characters to stdout.
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
 *  @verbatim vrna_message_constraints(option_1 | option_2 | option_n) @endverbatim
 *
 *  @ingroup  hard_constraints
 *
 *  @param option Option switch that tells which constraint help will be printed
 */
void vrna_message_constraint_options(unsigned int option);

/**
 *  @brief Print structure constraint characters to stdout
 *  (full constraint support)
 *
 *  @ingroup  hard_constraints
 *
 *  @see  vrna_message_constraint_options()
 */
void vrna_message_constraints_all(void);

/**
 *  @brief  Add hard constraints to a #vrna_fold_compound data structure
 *
 *  Use this function to add/update a data structure of type #vrna_hc_t that
 *  specifies which decomposition steps are allowed/enforced during the recursions.
 *  The function allows for passing a string 'constraint' that can either be a
 *  filename that points to a hard constraints definition file or it may be a
 *  pseudo dot-bracket notation indicating the hard constraint. Depending on
 *  the type of the string the user has to pass #VRNA_CONSTRAINT_FILE or
 *  #VRNA_CONSTRAINT_DB in the option parameter, respectively. If none of these
 *  to options are passed, no further hard constraints then the ones induced by
 *  canonical base pairing (as supplied with the ptype argument) are applied.
 *
 *  @see      vrna_hc_reset(), vrna_hc_free(), #VRNA_CONSTRAINT_FILE, #VRNA_CONSTRAINT_DB
 *
 *  @ingroup  hard_constraints
 *
 *  @param  vc            The fold compound
 *  @param  constraint    A string with either the filename of the hard constraint definitions
 *                        or a pseudo dot-bracket notation of the hard constraint. May be NULL.
 *  @param  options       The option flags
 */
void vrna_hc_add( vrna_fold_compound *vc,
                  const char *constraint,
                  unsigned int options);

void vrna_hc_add_up(vrna_fold_compound *vc,
                    int i,
                    char option);

void vrna_hc_add_bp(vrna_fold_compound *vc,
                    int i,
                    int j,
                    char option);

/**
 *  @brief  Reset hard constraints to default values
 *
 *  This function resets the hard constraints to its default values, i.e.
 *  all positions may be unpaired in all contexts, and base pairs are
 *  allowed in all contexts, if they resemble canonical pairs
 *
 *  @ingroup  hard_constraints
 *
 *  @param  vc            The fold compound
 */
void vrna_hc_reset(vrna_fold_compound *vc);


/**
 *  @brief  Free the memory allocated by a #vrna_hc_t data structure
 *
 *  Use this function to free all memory that was allocated for a data structure
 *  of type #vrna_hc_t .
 *
 *  @see get_hard_constraints(), #vrna_hc_t
 *  @ingroup  hard_constraints
 *
 */
void vrna_hc_free(vrna_hc_t *hc);

/**
 *  @brief Add soft constraints to a fold_compound
 *
 *  This function adds a proper soft constraints data structure
 *  to the fold_compound data structure.
 *  Current behavior upon already existing soft constraints is
 *  to remove them first before adding the new ones.
 *  This behavior will probably change in the near future in favor
 *  of adding additional soft constraints on top of others.
 *
 *  If no constraints are passed an empty soft constraints data
 *  structure is created.
 *  Otherwise, the provided constraint values are assumed to
 *  specify energy contributions in units of kcal/mol for unpaired
 *  positions.
 *
 *  @ingroup  soft_constraints
 *
 */
void vrna_sc_add( vrna_fold_compound *vc,
                  const double *constraints,
                  unsigned int options);

void vrna_sc_add_ali( vrna_fold_compound *vc,
                      const double **constraints,
                      unsigned int options);

void vrna_sc_add_bp(vrna_fold_compound *vc,
                    const double **constraints,
                    unsigned int options);

void vrna_sc_add_bp_mfe(vrna_fold_compound *vc,
                        const double **constraints,
                        unsigned int options);

void vrna_sc_add_bp_pf( vrna_fold_compound *vc,
                        const double **constraints,
                        unsigned int options);

void vrna_sc_add_up(vrna_fold_compound *vc,
                    const double *constraints,
                    unsigned int options);

void vrna_sc_add_up_mfe(vrna_fold_compound *vc,
                        const double *constraints,
                        unsigned int options);

void vrna_sc_add_up_pf( vrna_fold_compound *vc,
                        const double *constraints,
                        unsigned int options);

void vrna_sc_remove(vrna_fold_compound *vc);

void vrna_sc_free(vrna_sc_t *sc);

int vrna_sc_SHAPE_parse_method( const char *method_string,
                                char *method,
                                float *param_1,
                                float *param_2);

int vrna_sc_SHAPE_add_deigan( vrna_fold_compound *vc,
                              const double *reactivities,
                              double m,
                              double b,
                              unsigned int options);

int vrna_sc_SHAPE_add_deigan_ali( vrna_fold_compound *vc,
                                  const char **shape_files,
                                  const int *shape_file_association,
                                  double m,
                                  double b,
                                  unsigned int options);

int vrna_sc_SHAPE_add_zarringhalam( vrna_fold_compound *vc,
                                    const double *reactivities,
                                    double b,
                                    double default_value,
                                    const char *shape_conversion,
                                    unsigned int options);

/**
 *  @brief Convert SHAPE reactivity values to probabilities for being unpaired
 *
 *  This function parses the informations from a given file and stores the result
 *  in the preallocated string sequence and the double array values.
 *
 *  @ingroup  soft_constraints
 *
 *  @see vrna_read_SHAPE_file()
 *  @param shape_conversion String definining the method used for the conversion process
 *  @param values           Pointer to an array of SHAPE reactivities
 *  @param length           Length of the array of SHAPE reactivities
 *  @param default_value    Result used for position with invalid/missing reactivity values
 */
int vrna_sc_SHAPE_to_pr(const char *shape_conversion,
                        double *values,
                        int length,
                        double default_value);

/**
 *  @brief  Bind a function pointer for generalized soft constraint feature (MFE version)
 *
 *  This function allows to easily bind a function pointer and corresponding data structure
 *  to the soft constraint part #vrna_sc_t of the #vrna_fold_compound.
 *  The function for evaluating the generalized soft constraint feature has to return
 *  a pseudo free energy @f$ \hat{E} @f$ in @f$ dacal/mol @f$, where @f$ 1 dacal/mol = 10 cal/mol @f$.
 *
 *  @ingroup soft_constraints
 *  @param  vc    The fold compound the generalized soft constraint function should be bound to
 *  @param  f     A pointer to the function that evaluates the generalized soft constraint feature
 *  @param  data  A pointer to the data structure that holds required data for function 'f'
 */
void vrna_sc_add_f( vrna_fold_compound *vc,
                    int (*f)( int, int, int, int, char, void *),
                    void *data);

/**
 *  @brief  Bind a function pointer for generalized soft constraint feature (PF version)
 *
 *  This function allows to easily bind a function pointer and corresponding data structure
 *  to the soft constraint part #vrna_sc_t of the #vrna_fold_compound.
 *  The function for evaluating the generalized soft constraint feature has to return
 *  a pseudo free energy @f$ \hat{E} @f$ as Boltzmann factor, i.e. @f$ exp(- \hat{E} / kT) @f$.
 *  The required unit for @f$ E @f$ is @f$ cal/mol @f$.
 *
 *  @ingroup soft_constraints
 *  @param  vc    The fold compound the generalized soft constraint function should be bound to
 *  @param  exp_f A pointer to the function that evaluates the generalized soft constraint feature
 *  @param  data  A pointer to the data structure that holds required data for function 'f'
 */
void vrna_sc_add_exp_f( vrna_fold_compound *vc,
                        FLT_OR_DBL (*exp_f)( int, int, int, int, char, void *),
                        void *data);

/**
 *  @brief  Add a pre-processing function for the generalized soft constraints feature
 *
 *  @note     The function pointer passed will be used by calls of vrna_fold(), vrna_pf_fold(),
 *            vrna_ali_fold(), vrna_ali_pf_fold(), vrna_cofold(), vrna_co_pf_fold(), and
 *            vrna_subopt(). Each call is provided with the current #vrna_fold_compound, and
 *            a flag to indicate whether general free energy evaluation (#VRNA_SC_GEN_MFE), or
 *            partition function computations (#VRNA_SC_GEN_PF) will take place next.
 *  @see      #VRNA_SC_GEN_MFE, #VRNA_SC_GEN_PF, #vrna_sc_t, #vrna_fold_compound
 *  @ingroup  soft_constraints
 *  @param    vc    The fold compound the generalized soft constraint function should be bound to
 *  @param    post  A pointer to the pre-processing function
 */
void vrna_sc_add_pre(vrna_fold_compound *vc, void (*pre)( vrna_fold_compound *, char));

/**
 *  @brief  Add a post-processing function for the generalized soft constraints feature
 *
 *  @note     The function pointer passed will be used by calls of vrna_fold(), vrna_pf_fold(),
 *            vrna_ali_fold(), vrna_ali_pf_fold(), vrna_cofold(), vrna_co_pf_fold(), and
 *            vrna_subopt(). Each call is provided with the current #vrna_fold_compound, and
 *            a flag to indicate whether general free energy evaluation (#VRNA_SC_GEN_MFE), or
 *            partition function computations (#VRNA_SC_GEN_PF) has taken place before.
 *  @see      #VRNA_SC_GEN_MFE, #VRNA_SC_GEN_PF, #vrna_sc_t, #vrna_fold_compound
 *  @ingroup  soft_constraints
 *  @param    vc    The fold compound the generalized soft constraint function should be bound to
 *  @param    post  A pointer to the post-processing function
 */
void vrna_sc_add_post( vrna_fold_compound *vc, void (*post)( vrna_fold_compound *, char));

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
 *  @deprecated Use vrna_message_constraints_all() instead!
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
