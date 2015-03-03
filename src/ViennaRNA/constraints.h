#ifndef __VIENNA_RNA_PACKAGE_CONSTRAINTS_H__
#define __VIENNA_RNA_PACKAGE_CONSTRAINTS_H__

#include <ViennaRNA/data_structures.h>

/**
 *  \addtogroup constraints
 *  \ingroup folding_routines
 *  \brief This section covers all functions and variables related to the
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


/**
 *  \brief  Flag that is used to indicate the pipe '|' sign in pseudo dot-bracket
 *  notation of hard constraints.
 *
 *  Use this definition to indicate the pipe sign '|' (paired with another base)
 *  
 *  \ingroup  hard_constraints
 *
 */
#define VRNA_CONSTRAINT_PIPE              1U

/**
 *  \brief  dot '.' switch for structure constraints (no constraint at all)
 *  
 *  \ingroup  hard_constraints
 *
 */
#define VRNA_CONSTRAINT_DOT               2U
/**
 *  \brief  'x' switch for structure constraint (base must not pair)
 *  
 *  \ingroup  hard_constraints
 *
 */
#define VRNA_CONSTRAINT_X                 4U
/**
 *  \brief  angle brackets '<', '>' switch for structure constraint (paired downstream/upstream)
 *  
 *  \ingroup  hard_constraints
 *
 */
#define VRNA_CONSTRAINT_ANG_BRACK         8U
/**
 *  \brief  round brackets '(',')' switch for structure constraint (base i pairs base j)
 *  
 *  \ingroup  hard_constraints
 *
 */
#define VRNA_CONSTRAINT_RND_BRACK         16U

/**
 *  \brief  Flag that is used to indicate the character 'l' in pseudo dot-bracket
 *  notation of hard constraints.
 *
 *  Use this definition to indicate the usage of 'l' character (intramolecular pairs only)
 *
 *  \ingroup  hard_constraints
 *
 */
#define VRNA_CONSTRAINT_INTRAMOLECULAR    2048U

/**
 *  \brief  Flag that is used to indicate the character 'e' in pseudo dot-bracket
 *  notation of hard constraints.
 *
 *  Use this definition to indicate the usage of 'e' character (intermolecular pairs only)
 *
 *  \ingroup  hard_constraints
 *
 */
#define VRNA_CONSTRAINT_INTERMOLECULAR    4096U

/**
 *  \brief  constraint may span over several lines
 *  
 *  \ingroup  hard_constraints
 *
 */
#define VRNA_CONSTRAINT_MULTILINE         32U
/**
 *  \brief  do not print the header information line
 *  
 *  \ingroup  hard_constraints
 *
 */
#define VRNA_CONSTRAINT_NO_HEADER         64U
/**
 *  \brief  placeholder for all constraining characters
 *  
 *  \ingroup  hard_constraints
 *
 */
#define VRNA_CONSTRAINT_ALL               128U

/**
 *  \brief  
 *  
 *  \ingroup  hard_constraints
 *
 */
#define VRNA_CONSTRAINT_DB                256U

/**
 *  \brief  
 *  
 *  \ingroup  hard_constraints
 *
 */
#define VRNA_CONSTRAINT_FILE              512U

/**
 *  \brief  
 *  
 *  \ingroup  hard_constraints
 *
 */
#define VRNA_CONSTRAINT_IINDX             1024U

/**
 *  \brief  Soft constraints flag, apply constraints for MFE calculations
 *  
 *  \ingroup  soft_constraints
 *
 */
#define VRNA_CONSTRAINT_SOFT_MFE          8192U

/**
 *  \brief  Soft constraints flag, apply constraints for partition function calculations
 *  
 *  \ingroup  soft_constraints
 *
 */
#define VRNA_CONSTRAINT_SOFT_PF           16384U


/**
 *  \brief  Soft constraints flag, apply constraints for unpaired nucleotides
 *  
 *  \ingroup  soft_constraints
 *
 */
#define VRNA_CONSTRAINT_SOFT_UP           32768U

/**
 *  \brief  Soft constraints flag, apply constraints for base pairs
 *  
 *  \ingroup  soft_constraints
 *
 */
#define VRNA_CONSTRAINT_SOFT_BP           65536U

/**
 *  \brief  Soft constraints flag, remove gaps from aligned sequence
 *  
 *  \ingroup  soft_constraints
 *
 */
#define VRNA_CONSTRAINT_UNGAP             131072U

/**
 *  \brief  Hard constraints flag, base pair in the exterior loop
 *  
 *  \ingroup  hard_constraints
 *
 */
#define VRNA_HC_CONTEXT_EXT_LOOP       (char)0x01
/**
 *  \brief  Hard constraints flag, base pair encloses hairpin loop
 *  
 *  \ingroup  hard_constraints
 *
 */
#define VRNA_HC_CONTEXT_HP_LOOP        (char)0x02
/**
 *  \brief  Hard constraints flag, base pair encloses an interior loop
 *  
 *  \ingroup  hard_constraints
 *
 */
#define VRNA_HC_CONTEXT_INT_LOOP       (char)0x04
/**
 *  \brief  Hard constraints flag, base pair is enclosed in an interior loop
 *  
 *  \ingroup  hard_constraints
 *
 */
#define VRNA_HC_CONTEXT_MB_LOOP        (char)0x08
/**
 *  \brief  Hard constraints flag, base pair encloses a multi branch loop
 *
 *  \ingroup  hard_constraints
 *
 */
#define VRNA_HC_CONTEXT_INT_LOOP_ENC  (char)0x10
/**
 *  \brief  Hard constraints flag, base pair is enclosed in a multi branch loop
 *  
 *  \ingroup  hard_constraints
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
 *  \brief  Generalized constraint folding flag indicating hairpin loop decomposition step
 *
 *  \ingroup  constraints
 *
 */
#define VRNA_DECOMP_PAIR_HP     1

/**
 *  \brief  Generalized constraint folding flag indicating interior loop decomposition step
 *
 *  \ingroup  constraints
 *
 */
#define VRNA_DECOMP_PAIR_IL     2

/**
 *  \ingroup  constraints
 *
 */
#define VRNA_DECOMP_PAIR_ML     3

/**
 *  \ingroup  constraints
 *
 */
#define VRNA_DECOMP_ML_ML_ML    5

/**
 *  \ingroup  constraints
 *
 */
#define VRNA_DECOMP_ML_UP_3     4

/**
 *  \ingroup  constraints
 *
 */
#define VRNA_DECOMP_ML_UP_5     6

/**
 *  \ingroup  constraints
 *
 */
#define VRNA_DECOMP_ML_UP       11

/**
 *  \ingroup  constraints
 *
 */
#define VRNA_DECOMP_EXT_EXT     9

/**
 *  \ingroup  constraints
 *
 */
#define VRNA_DECOMP_EXT_UP_3    7

/**
 *  \ingroup  constraints
 *
 */
#define VRNA_DECOMP_EXT_UP_5    10

/**
 *  \ingroup  constraints
 *
 */
#define VRNA_DECOMP_EXT_UP      8

/**
 *  \ingroup  constraints
 *
 */
#define VRNA_DECOMP_EXT_STEM_UP 12

/**
 *  \brief  The hard constraints data structure
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
 *  \see  get_hard_constraints(), vrna_hc_free(), #VRNA_HC_CONTEXT_EXT_LOOP,
 *        #VRNA_HC_CONTEXT_HP_LOOP, #VRNA_HC_CONTEXT_INT_LOOP, #VRNA_HC_CONTEXT_EXT_LOOP_ENC, #VRNA_HC_CONTEXT_MB_LOOP, #VRNA_HC_CONTEXT_MB_LOOP_ENC
 *        
 *  \ingroup hard_constraints
 */
typedef struct vrna_hc_t {
  char    *matrix;  /**<  \brief  Upper triangular matrix encoding where a
                                  base pair or unpaired nucleotide is allowed
                    */
  int     *up_ext;  /**<  \brief  A linear array that holds the number of allowed
                                  unpaired nucleotides in an exterior loop
                    */
  int     *up_hp;   /**<  \brief  A linear array that holds the number of allowed
                                  unpaired nucleotides in a hairpin loop
                    */
  int     *up_int;  /**<  \brief  A linear array that holds the number of allowed
                                  unpaired nucleotides in an interior loop
                    */
  int     *up_ml;   /**<  \brief  A linear array that holds the number of allowed
                                  unpaired nucleotides in a multi branched loop
                    */
} vrna_hc_t;

/**
 *  \brief  The soft constraints data structure
 *
 *  \ingroup soft_constraints
 */
typedef struct vrna_sc_t {
  double      *constraints;         /**<  \brief Backup storage for energy contributions of single nucleotides */
  int         **free_energies;      /**<  \brief Energy contribution for unpaired sequence stretches */
  int         *en_basepair;         /**<  \brief Energy contribution for base pairs */
  FLT_OR_DBL  **boltzmann_factors;  /**<  \brief Boltzmann Factors of the energy contributions for unpaired sequence stretches */
  FLT_OR_DBL  *exp_en_basepair;     /**<  \brief Boltzmann Factors of the energy contribution for base pairs */

  int         *en_stack;            /**<  \brief Pseudo Energy contribution per base pair involved in a stack */
  FLT_OR_DBL  *exp_en_stack;        /**<  \brief Boltzmann weighted pseudo energy contribution per nucleotide involved in a stack */

  /* generalized soft contraints */
  int (*f)( int,
            int,
            int,
            int,
            char,
            void *);            /**<  \brief  A function pointer used for pseudo
                                              energy contribution in MFE calculations
                                */

  FLT_OR_DBL (*exp_f)(int,
                      int,
                      int,
                      int,
                      char,
                      void *);  /**<  \brief  A function pointer used for pseudo energy
                                              contribution boltzmann factors in PF
                                              calculations
                                */

  void *data;                   /**<  \brief  A pointer to the data object necessary for
                                              for pseudo energy contribution functions
                                */
} vrna_sc_t;

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
 *
 *  \ingroup  hard_constraints
 *
 *  \param option Option switch that tells which constraint help will be printed
 */
void print_tty_constraint(unsigned int option);

/**
 *  \brief Print structure constraint characters to stdout
 *  (full constraint support)
 *
 *  \ingroup  hard_constraints
 *
 *
 */
void print_tty_constraint_full(void);

/**
 *  \brief Insert constraining pair types according to constraint structure string
 *
 *  \see vrna_get_indx(), vrna_get_iindx()
 *
 *  \ingroup  hard_constraints
 *
 *
 *  \param constraint     The structure constraint string
 *  \param length         The actual length of the sequence (constraint may be shorter)
 *  \param ptype          A pointer to the basepair type array
 *  \param BP             (not used anymore)
 *  \param min_loop_size  The minimal loop size (usually \ref TURN )
 *  \param idx_type       Define the access type for base pair type array (0 = indx, 1 = iindx)
 */
void constrain_ptypes(const char *constraint,
                      unsigned int length,
                      char *ptype,
                      int *BP,
                      int min_loop_size,
                      unsigned int idx_type);

/**
 *  \brief  Get a hard constraint from pseudo dot-bracket notation as specified
 *          in the ViennaRNA Package extension of the FASTA format.
 *
 *  \pre      The argument 'lines' has to be a 2-dimensional character array as obtained
 *            by read_record()
 *  \see      read_record(), #VRNA_CONSTRAINT_PIPE, #VRNA_CONSTRAINT_DOT, #VRNA_CONSTRAINT_X
 *            #VRNA_CONSTRAINT_ANG_BRACK, #VRNA_CONSTRAINT_RND_BRACK
 *  \ingroup  hard_constraints
 *
 *  \param  cstruc  A pointer to a character array that is used as pseudo dot-bracket
 *                  output
 *  \param  lines   A 2-dimensional character array with the extension lines from the FASTA
 *                  input
 *  \param  option  The option flags that define the behavior and recognition pattern of
 *                  this function
 */
void getConstraint( char **cstruc,
                    const char **lines,
                    unsigned int option);

/**
 *  \brief  Add hard constraints to a #vrna_fold_compound data structure
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
 *  \see      vrna_hc_reset(), vrna_hc_free(), #VRNA_CONSTRAINT_FILE, #VRNA_CONSTRAINT_DB
 *
 *  \ingroup  hard_constraints
 *
 *  \param  vc            The fold compound
 *  \param  constraint    A string with either the filename of the hard constraint definitions
 *                        or a pseudo dot-bracket notation of the hard constraint. May be NULL.
 *  \param  options       The option flags
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
 *  \brief  Reset hard constraints to default values
 *
 *  This function resets the hard constraints to its default values, i.e.
 *  all positions may be unpaired in all contexts, and base pairs are
 *  allowed in all contexts, if they resemble canonical pairs
 *
 *  \ingroup  hard_constraints
 *
 *  \param  vc            The fold compound
 */
void vrna_hc_reset(vrna_fold_compound *vc);


/**
 *  \brief  Free the memory allocated by a #vrna_hc_t data structure
 *
 *  Use this function to free all memory that was allocated for a data structure
 *  of type #vrna_hc_t .
 *
 *  \see get_hard_constraints(), #vrna_hc_t
 *  \ingroup  hard_constraints
 *
 */
void vrna_hc_free(vrna_hc_t *hc);


/**
  * \brief Parse information from a given soft constraints file
  *
  * This function parses the informations from a given file and stores the result
  * in the preallocated string sequence and the double array values.
  *
  * \param file_name     Path to the constraints file
  * \param length        Length of the sequence (file entries exceeding this limit will cause an error)
  * \param default_value Value for missing indices
  * \param sequence      Pointer to an array used for storring the sequence obtained from the constraints file
  * \param values        Pointer to an array used for storring the values obtained from the constraints file
  */
int parse_soft_constraints_file(const char *file_name,
                                int length,
                                double default_value,
                                char *sequence,
                                double *values);


/**
  * \brief Convert SHAPE reactivity values to probabilities for being unpaired
  *
  * This function parses the informations from a given file and stores the result
  * in the preallocated string sequence and the double array values.
  *
  * \param shape_conversion String definining the method used for the conversion process
  * \param values           Pointer to an array of SHAPE reactivities
  * \param length           Length of the array of SHAPE reactivities
  * \param default_value    Result used for position with invalid/missing reactivity values
  */
int convert_shape_reactivities_to_probabilities(const char *shape_conversion,
                                                double *values,
                                                int length,
                                                double default_value);

/**
 *  \brief Add soft constraints to a fold_compound
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
 */
void vrna_sc_add( vrna_fold_compound *vc,
                  const double *constraints,
                  unsigned int options);

void vrna_sc_add_ali( vrna_fold_compound *vc,
                      const double **constraints,
                      unsigned int options);

int vrna_sc_add_deigan_ali(vrna_fold_compound *vc,
                            const char **shape_files,
                            const int *shape_file_association,
                            double m,
                            double b,
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

void vrna_sc_add_sp(vrna_fold_compound *vc,
                    const double *constraints,
                    unsigned int options);

void vrna_sc_add_sp_mfe(vrna_fold_compound *vc,
                        const double *constraints,
                        unsigned int options);

void vrna_sc_add_sp_pf( vrna_fold_compound *vc,
                        const double *constraints,
                        unsigned int options);


void vrna_sc_remove(vrna_fold_compound *vc);

void vrna_sc_free(vrna_sc_t *sc);

int parse_soft_constraints_shape_method(const char *method_string,
                                        char *method,
                                        float *param_1,
                                        float *param_2);

#endif
