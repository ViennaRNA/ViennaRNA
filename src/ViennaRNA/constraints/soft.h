#ifndef VIENNA_RNA_PACKAGE_CONSTRAINTS_SOFT_H
#define VIENNA_RNA_PACKAGE_CONSTRAINTS_SOFT_H

/**
 *  @file     constraints/soft.h
 *  @ingroup  soft_constraints
 *  @brief    Functions and data structures for secondary structure soft constraints
 */


/**
 *  @brief Typename for the soft constraints data structure #vrna_sc_s
 *
 *  @ingroup soft_constraints
 */
typedef struct  vrna_sc_s vrna_sc_t;

#include <stdlib.h>

#include <ViennaRNA/datastructures/basic.h>
#include <ViennaRNA/fold_compound.h>
#include <ViennaRNA/constraints/basic.h>

/**
 * @brief Callback to retrieve pseudo energy contribution for soft constraint feature
 *
 * This is the prototype for callback functions used by the folding recursions to evaluate generic
 * soft constraints. The first four parameters passed indicate the delimiting nucleotide positions
 * of the decomposition, and the parameter @p denotes the decomposition step. The last parameter
 * @p data is the auxiliary data structure associated to the hard constraints via vrna_sc_add_data(),
 * or NULL if no auxiliary data was added.
 *
 * @callback
 * @parblock
 * This callback enables one to add (pseudo-)energy contributions to individual decompositions
 * of the secondary structure.
 * @endparblock
 *
 * @ingroup soft_constraints_generic
 *
 * @see #VRNA_DECOMP_PAIR_HP, #VRNA_DECOMP_PAIR_IL, #VRNA_DECOMP_PAIR_ML, #VRNA_DECOMP_ML_ML_ML,
 *      #VRNA_DECOMP_ML_STEM, #VRNA_DECOMP_ML_ML, #VRNA_DECOMP_ML_UP, #VRNA_DECOMP_ML_ML_STEM,
 *      #VRNA_DECOMP_ML_COAXIAL, #VRNA_DECOMP_EXT_EXT, #VRNA_DECOMP_EXT_UP, #VRNA_DECOMP_EXT_STEM,
 *      #VRNA_DECOMP_EXT_EXT_EXT, #VRNA_DECOMP_EXT_STEM_EXT, #VRNA_DECOMP_EXT_EXT_STEM,
 *      #VRNA_DECOMP_EXT_EXT_STEM1, vrna_sc_add_f(), vrna_sc_add_exp_f(), vrna_sc_add_bt(),
 *      vrna_sc_add_data()
 *
 * @param i         Left (5') delimiter position of substructure
 * @param j         Right (3') delimiter position of substructure
 * @param k         Left delimiter of decomposition
 * @param l         Right delimiter of decomposition
 * @param d         Decomposition step indicator
 * @param data      Auxiliary data
 * @return          Pseudo energy contribution in deka-kalories per mol
 */
typedef int (*vrna_sc_f)(int            i,
                         int            j,
                         int            k,
                         int            l,
                         unsigned char  d,
                         void           *data);

/**
 *
 * @ingroup soft_constraints_generic
 *
 */
DEPRECATED(typedef int (vrna_callback_sc_energy)(int i,
                                                 int j,
                                                 int k,
                                                 int l,
                                                 unsigned char d,
                                                 void *data),
           "Use vrna_sc_f instead!");


/**
 *
 * @ingroup soft_constraints_generic
 *
 */
typedef int (*vrna_sc_direct_f)(vrna_fold_compound_t  *fc,
                                int                   i,
                                int                   j,
                                int                   k,
                                int                   l,
                                void                  *data);

/**
 * @brief Callback to retrieve pseudo energy contribution as Boltzmann Factors for soft constraint feature
 *
 * This is the prototype for callback functions used by the partition function recursions to evaluate generic
 * soft constraints. The first four parameters passed indicate the delimiting nucleotide positions
 * of the decomposition, and the parameter @p denotes the decomposition step. The last parameter
 * @p data is the auxiliary data structure associated to the hard constraints via vrna_sc_add_data(),
 * or NULL if no auxiliary data was added.
 *
 * @callback
 * @parblock
 * This callback enables one to add (pseudo-)energy contributions to individual decompositions
 * of the secondary structure (Partition function variant, i.e. contributions must be returned as Boltzmann factors).
 * @endparblock
 *
 * @ingroup soft_constraints_generic
 *
 * @see #VRNA_DECOMP_PAIR_HP, #VRNA_DECOMP_PAIR_IL, #VRNA_DECOMP_PAIR_ML, #VRNA_DECOMP_ML_ML_ML,
 *      #VRNA_DECOMP_ML_STEM, #VRNA_DECOMP_ML_ML, #VRNA_DECOMP_ML_UP, #VRNA_DECOMP_ML_ML_STEM,
 *      #VRNA_DECOMP_ML_COAXIAL, #VRNA_DECOMP_EXT_EXT, #VRNA_DECOMP_EXT_UP, #VRNA_DECOMP_EXT_STEM,
 *      #VRNA_DECOMP_EXT_EXT_EXT, #VRNA_DECOMP_EXT_STEM_EXT, #VRNA_DECOMP_EXT_EXT_STEM,
 *      #VRNA_DECOMP_EXT_EXT_STEM1, vrna_sc_add_exp_f(), vrna_sc_add_f(), vrna_sc_add_bt(),
 *      vrna_sc_add_data()
 *
 * @param i         Left (5') delimiter position of substructure
 * @param j         Right (3') delimiter position of substructure
 * @param k         Left delimiter of decomposition
 * @param l         Right delimiter of decomposition
 * @param d         Decomposition step indicator
 * @param data      Auxiliary data
 * @return          Pseudo energy contribution in deka-kalories per mol
 */
typedef FLT_OR_DBL (*vrna_sc_exp_f)(int           i,
                                    int           j,
                                    int           k,
                                    int           l,
                                    unsigned char d,
                                    void          *data);

/**
 *
 * @ingroup soft_constraints_generic
 *
 */
DEPRECATED(typedef FLT_OR_DBL (vrna_callback_sc_exp_energy)(int i,
                                                            int j,
                                                            int k,
                                                            int l,
                                                            unsigned char d,
                                                            void *data),
           "Use vrna_sc_exp_f instead!");


/**
 *
 * @ingroup soft_constraints_generic
 *
 */
typedef FLT_OR_DBL (*vrna_sc_exp_direct_f)(vrna_fold_compound_t *fc,
                                           int                  i,
                                           int                  j,
                                           int                  k,
                                           int                  l,
                                           void                 *data);


/**
 * @brief Callback to retrieve auxiliary base pairs for soft constraint feature
 *
 * @callback
 * @parblock
 * This callback enables one to add auxiliary base pairs in the backtracking steps
 * of hairpin- and internal loops.
 * @endparblock
 *
 * @ingroup soft_constraints_generic
 *
 * @see #VRNA_DECOMP_PAIR_HP, #VRNA_DECOMP_PAIR_IL, #VRNA_DECOMP_PAIR_ML, #VRNA_DECOMP_ML_ML_ML,
 *      #VRNA_DECOMP_ML_STEM, #VRNA_DECOMP_ML_ML, #VRNA_DECOMP_ML_UP, #VRNA_DECOMP_ML_ML_STEM,
 *      #VRNA_DECOMP_ML_COAXIAL, #VRNA_DECOMP_EXT_EXT, #VRNA_DECOMP_EXT_UP, #VRNA_DECOMP_EXT_STEM,
 *      #VRNA_DECOMP_EXT_EXT_EXT, #VRNA_DECOMP_EXT_STEM_EXT, #VRNA_DECOMP_EXT_EXT_STEM,
 *      #VRNA_DECOMP_EXT_EXT_STEM1, vrna_sc_add_bt(), vrna_sc_add_f(), vrna_sc_add_exp_f(),
 *      vrna_sc_add_data()
 *
 * @param i         Left (5') delimiter position of substructure
 * @param j         Right (3') delimiter position of substructure
 * @param k         Left delimiter of decomposition
 * @param l         Right delimiter of decomposition
 * @param d         Decomposition step indicator
 * @param data      Auxiliary data
 * @return          List of additional base pairs
 */
typedef vrna_basepair_t *(*vrna_sc_bt_f)(int            i,
                                         int            j,
                                         int            k,
                                         int            l,
                                         unsigned char  d,
                                         void           *data);
/**
 *
 * @ingroup soft_constraints_generic
 *
 */

DEPRECATED(typedef vrna_basepair_t *(vrna_callback_sc_backtrack)(int i,
                                                                 int j,
                                                                 int k,
                                                                 int l,
                                                                 unsigned char d,
                                                                 void *data),
           "Use vrna_sc_bt_f instead");


/**
 *  @brief  The type of a soft constraint
 *
 *  @ingroup soft_constraints
 */
typedef enum {
  VRNA_SC_DEFAULT,  /**<  @brief  Default Soft Constraints */
  VRNA_SC_WINDOW    /**<  @brief  Soft Constraints suitable for local structure prediction using
                     *    window approach.
                     *    @see    vrna_mfe_window(), vrna_mfe_window_zscore(), pfl_fold()
                     */
} vrna_sc_type_e;


/**
 *  @brief  A base pair constraint
 *
 *  @ingroup soft_constraints
 */
typedef struct {
  unsigned int  interval_start;
  unsigned int  interval_end;
  int           e;
} vrna_sc_bp_storage_t;


/**
 *  @brief  The soft constraints data structure
 *
 *  @ingroup soft_constraints
 */
struct vrna_sc_s {
  /**
   *  @name Common data fields
   *  @{
   */
  const vrna_sc_type_e  type; /**< @brief Type of the soft constraints data structure */
  unsigned int          n;    /**< @brief Length of the sequence this soft constraints data structure belongs to */

  unsigned char         state; /**< @brief  Current state of the soft constraints data structure */

  int                   **energy_up;      /**<  @brief Energy contribution for stretches of unpaired nucleotides */
  FLT_OR_DBL            **exp_energy_up;  /**<  @brief Boltzmann Factors of the energy contributions for unpaired sequence stretches */

  int                   *up_storage;      /**<  @brief  Storage container for energy contributions per unpaired nucleotide */
  vrna_sc_bp_storage_t  **bp_storage;     /**<  @brief  Storage container for energy contributions per base pair */

  int           *energy_stack;                    /**<  @brief Pseudo Energy contribution per base pair involved in a stack */
  FLT_OR_DBL    *exp_energy_stack;                /**<  @brief Boltzmann weighted pseudo energy contribution per nucleotide involved in a stack */
  /**
   *  @}
   *
   */
#ifndef VRNA_DISABLE_C11_FEATURES
  /* C11 support for unnamed unions/structs */
  union {
    struct {
#endif
  /**
   *  @name Global structure prediction data fields
   *  @{
   */
  int *energy_bp;                               /**<  @brief Energy contribution for base pairs */
  FLT_OR_DBL *exp_energy_bp;                    /**<  @brief Boltzmann Factors of the energy contribution for base pairs */
  /**
   *  @}
   *
   */
#ifndef VRNA_DISABLE_C11_FEATURES
  /* C11 support for unnamed unions/structs */
};
struct {
#endif
  /**
   *  @name Local structure prediction data fields
   *  @{
   */
  int         **energy_bp_local;                        /**<  @brief Energy contribution for base pairs (sliding window approach) */
  FLT_OR_DBL  **exp_energy_bp_local;                    /**<  @brief Boltzmann Factors of the energy contribution for base pairs (sliding window approach) */
  /**
   *  @}
   *
   */
#ifndef VRNA_DISABLE_C11_FEATURES
  /* C11 support for unnamed unions/structs */
};
};
#endif

  /**
   *  @name User-defined data fields
   *  @{
   */
  /* generic soft contraints below */
  vrna_sc_f     f;            /**<  @brief  A function pointer used for pseudo
                               *            energy contribution in MFE calculations
                               *    @see    vrna_sc_add_f()
                               */

  vrna_sc_bt_f  bt;           /**<  @brief  A function pointer used to obtain backtraced
                               *            base pairs in loop regions that were altered
                               *            by soft constrained pseudo energy contributions
                               *    @see    vrna_sc_add_bt()
                               */

  vrna_sc_exp_f exp_f;        /**<  @brief  A function pointer used for pseudo energy
                               *            contribution boltzmann factors in PF
                               *            calculations
                               *    @see    vrna_sc_add_exp_f()
                               */

  void                *data;  /**<  @brief  A pointer to the data object provided for
                               *            for pseudo energy contribution functions of the
                               *            generic soft constraints feature
                               */

  vrna_auxdata_prepare_f  prepare_data;
  vrna_auxdata_free_f     free_data;

  /**
   *  @}
   */
};

/**
 *  @brief Initialize an empty soft constraints data structure within a #vrna_fold_compound_t
 *
 *  This function adds a proper soft constraints data structure
 *  to the #vrna_fold_compound_t data structure.
 *  If soft constraints already exist within the fold compound, they are removed.
 *
 *  @note Accepts vrna_fold_compound_t of type #VRNA_FC_TYPE_SINGLE and #VRNA_FC_TYPE_COMPARATIVE
 *
 *  @ingroup soft_constraints
 *
 *  @see  vrna_sc_set_bp(), vrna_sc_set_up(), vrna_sc_add_SHAPE_deigan(),
 *        vrna_sc_add_SHAPE_zarringhalam(), vrna_sc_remove(), vrna_sc_add_f(),
 *        vrna_sc_add_exp_f(), vrna_sc_add_pre(), vrna_sc_add_post()
 *
 *  @param  fc  The #vrna_fold_compound_t where an empty soft constraint feature is to be added to
 */
void
vrna_sc_init(vrna_fold_compound_t *fc);


void
vrna_sc_init_window(vrna_fold_compound_t *fc);


/**
 *  @brief  Prepare soft constraints
 *
 *  @ingroup soft_constraints
 */
int
vrna_sc_prepare(vrna_fold_compound_t  *fc,
                unsigned int          options);


/**
 *  @brief  Update/prepare soft constraints for sliding-window computations
 *
 *  @ingroup soft_constraints
 */
int
vrna_sc_update(vrna_fold_compound_t *fc,
               unsigned int         i,
               unsigned int         options);


/**
 *  @brief  Remove soft constraints from #vrna_fold_compound_t
 *
 *  @note Accepts #vrna_fold_compound_t of type #VRNA_FC_TYPE_SINGLE and #VRNA_FC_TYPE_COMPARATIVE
 *
 *  @ingroup soft_constraints
 *
 *  @param  fc  The #vrna_fold_compound_t possibly containing soft constraints
 */
void
vrna_sc_remove(vrna_fold_compound_t *fc);


/**
 *  @brief  Free memory occupied by a #vrna_sc_t data structure
 *
 *  @ingroup soft_constraints
 *
 *  @param  sc  The data structure to free from memory
 */
void
vrna_sc_free(vrna_sc_t *sc);


/**
 *  @addtogroup soft_constraints_bp
 *  @{
 */

/**
 *  @brief  Set soft constraints for paired nucleotides
 *
 *  @note     This function replaces any pre-exisitng soft constraints with the ones supplied
 *            in @p constraints.
 *
 *  @see      vrna_sc_add_bp(), vrna_sc_set_up(), vrna_sc_add_up()
 *
 *  @param  fc          The #vrna_fold_compound_t the soft constraints are associated with
 *  @param  constraints A two-dimensional array of pseudo free energies in @f$ kcal / mol @f$
 *  @param  options     The options flag indicating how/where to store the soft constraints
 *  @return             Non-zero on successful application of the constraint, 0 otherwise.
 */
int
vrna_sc_set_bp(vrna_fold_compound_t *fc,
               const FLT_OR_DBL     **constraints,
               unsigned int         options);


/**
 *  @brief  Set soft constraints for paired nucleotides in comparative structure predictions
 *
 *  Similar to vrna_sc_set_bp() this function allows to set soft constraints @f$ e_{i,j}^\mathrm{BP} @f$
 *  for all base pairs @f$ (i, j) @f$ at once using a 1-based upper-triangular matrix @f$ E_\mathrm{BP} @f$.
 *  Since this function is supposed to be used for comparative structure predictions over a
 *  multiple sequence alignment (MSA), a 0-based array of matrices must be supplied as
 *  parameter @p constraints. If no constraints are to be used for sequence @f$ s @f$ in the
 *  MSA, the corresponding entry may be set to @b NULL.
 *
 *  @note     This function replaces any pre-exisitng soft constraints with the ones supplied
 *            in @p constraints.
 *
 *  @warning  Currently, base pair constraints must be provided in alignment coordinates
 *            rather than sequence coordinates! This may change in the future!
 *
 *  @see      vrna_sc_set_bp_comparative_seq(), vrna_sc_add_bp_comparative(),
 *            vrna_sc_set_up_comparative(), vrna_sc_add_up_comparative(),
 *            vrna_sc_set_stack_comparative(), vrna_sc_add_stack_comparative(),
 *            vrna_sc_set_bp()
 *
 *  @param  fc          The #vrna_fold_compound_t the soft constraints are associated with
 *  @param  constraints A 0-based array of 1-based two-dimensional arrays with pseudo free energies in @f$ kcal / mol @f$
 *  @param  options     The options flag indicating how/where to store the soft constraints
 *  @return             The number of sequences in the MSA constraints have been applied to
 */
int
vrna_sc_set_bp_comparative(vrna_fold_compound_t *fc,
                           const FLT_OR_DBL     ***constraints,
                           unsigned int         options);


/**
 *  @brief  Set soft constraints for paired nucleotides in comparative structure predictions
 *
 *  This is a convenience wrapper for vrna_sc_set_bp_comparative() where only one particular
 *  sequence @p s is provided with constraints.
 *
 *  @note     This function replaces any pre-exisitng soft constraints with the ones supplied
 *            in @p constraints.
 *
 *  @warning  This function not only re-sets the constraints for sequence @p s
 *            in the MSA but will also remove all constraints for all other sequences!
 *
 *  @warning  Currently, base pair constraints must be provided in alignment coordinates
 *            rather than sequence coordinates! This may change in the future!
 *
 *  @see      vrna_sc_set_bp_comparative(), vrna_sc_add_bp_comparative(),
 *            vrna_sc_set_up_comparative_seq(), vrna_sc_add_up_comparative_seq(),
 *            vrna_sc_set_stack_comparative_seq(), vrna_sc_add_stack_comparative_seq(),
 *            vrna_sc_set_bp()
 *
 *  @param  fc          The #vrna_fold_compound_t the soft constraints are associated with
 *  @param  s           The 0-based number of the sequence in the alignment the constraints are provided for
 *  @param  constraints A 1-based two-dimensional array of pseudo free energies in @f$ kcal / mol @f$
 *  @param  options     The options flag indicating how/where to store the soft constraints
 *  @return             The number of sequences in the MSA constraints have been applied to
 */
int
vrna_sc_set_bp_comparative_seq(vrna_fold_compound_t *fc,
                               unsigned int         s,
                               const FLT_OR_DBL     **constraints,
                               unsigned int         options);


/**
 *  @brief  Add soft constraints for paired nucleotides
 *
 *  @see      vrna_sc_set_bp(), vrna_sc_set_up(), vrna_sc_add_up()
 *
 *  @param  fc          The #vrna_fold_compound_t the soft constraints are associated with
 *  @param  i           The 5' position of the base pair the soft constraint is added for
 *  @param  j           The 3' position of the base pair the soft constraint is added for
 *  @param  energy      The free energy (soft-constraint) in @f$ kcal / mol @f$
 *  @param  options     The options flag indicating how/where to store the soft constraints
 *  @return             Non-zero on successful application of the constraint, 0 otherwise.
 */
int
vrna_sc_add_bp(vrna_fold_compound_t *fc,
               unsigned int         i,
               unsigned int         j,
               FLT_OR_DBL           energy,
               unsigned int         options);


/**
 *  @brief  Add soft constraints for paired nucleotides in comparative structure predictions
 *
 *  Similar to vrna_sc_add_bp(), this function allows to add soft constraints @f$ e_{i,j}^\mathrm{BP} @f$
 *  for all base pairs @f$ (i, j) @f$ in the multiple sequence alignment (MSA).
 *  The actual pairing partners @f$ i @f$ and @f$ j @f$ for each sequence in the MSA are
 *  provided in the form of 0-based arrays as parameters @p is and @p js. The corresponding
 *  energy contributions are provided as 0-based array in parameter @p energies. If no constraint
 *  is provided for sequence @f$ s @f$ in the MSA, the corresponding @p is value must be set
 *  to 0.
 *
 *  @note     Consecutive calls of this function with the same @p is and @p js accumulate
 *            to corresponding @p energies values, i.e. energies are added up onto each other.
 *
 *  @warning  Currently, base pair constraints must be provided in alignment coordinates
 *            rather than sequence coordinates! This may change in the future!
 *
 *  @see      vrna_sc_add_bp_comparative_seq(), vrna_sc_set_bp_comparative(),
 *            vrna_sc_set_up_comparative(), vrna_sc_add_up_comparative(),
 *            vrna_sc_set_stack_comparative(), vrna_sc_add_stack_comparative(),
 *            vrna_sc_add_bp()
 *
 *  @param  fc          The #vrna_fold_compound_t the soft constraints are associated with
 *  @param  is          A 0-based array of 5' position of the base pairs the soft constraint is added for
 *  @param  js          A 0-based array of 3' position of the base pairs the soft constraint is added for
 *  @param  energies    A 0-based array of free energies (soft-constraint) in @f$ kcal / mol @f$
 *  @param  options     The options flag indicating how/where to store the soft constraints
 *  @return             The number of sequences in the MSA constraints have been applied to
 */
int
vrna_sc_add_bp_comparative(vrna_fold_compound_t *fc,
                           unsigned int         *is,
                           unsigned int         *js,
                           const FLT_OR_DBL     *energies,
                           unsigned int         options);


/**
 *  @brief  Add soft constraints for paired nucleotides in comparative structure predictions
 *
 *  This is a convenience wrapper for vrna_sc_add_bp_comparative() where only one particular
 *  sequence @p s is provided with constraints.
 *
 *  @note     Consecutive calls of this function with the same @p i and @p j accumulate
 *            to corresponding @p energies values, i.e. energies are added up onto each other.
 *
 *  @warning  Currently, base pair constraints must be provided in alignment coordinates
 *            rather than sequence coordinates! This may change in the future!
 *
 *  @see      vrna_sc_add_bp_comparative(), vrna_sc_set_bp_comparative_seq(),
 *            vrna_sc_set_up_comparative_seq(), vrna_sc_add_up_comparative_seq(),
 *            vrna_sc_set_stack_comparative_seq(), vrna_sc_add_stack_comparative_seq(),
 *            vrna_sc_add_bp()
 *
 *  @param  fc        The #vrna_fold_compound_t the soft constraints are associated with
 *  @param  s         The 0-based number of the sequence in the alignment the constraints are provided for
 *  @param  i         5' position of the base pairs the soft constraint is added for
 *  @param  j         3' position of the base pairs the soft constraint is added for
 *  @param  energy    Free energy (soft-constraint) in @f$ kcal / mol @f$
 *  @param  options   The options flag indicating how/where to store the soft constraints
 *  @return           The number of sequences in the MSA constraints have been applied to
 */
int
vrna_sc_add_bp_comparative_seq(vrna_fold_compound_t *fc,
                               unsigned int         s,
                               unsigned int         i,
                               unsigned int         j,
                               FLT_OR_DBL           energy,
                               unsigned int         options);


/**
 *  @}
 *  @addtogroup soft_constraints_up
 *  @{
 */

/**
 *  @brief  Set soft constraints for unpaired nucleotides
 *
 *  @note     This function replaces any pre-exisitng soft constraints with the ones supplied
 *            in @p constraints.
 *
 *  @see      vrna_sc_add_up(), vrna_sc_set_bp(), vrna_sc_add_bp()
 *
 *  @param  fc          The #vrna_fold_compound_t the soft constraints are associated with
 *  @param  constraints A vector of pseudo free energies in @f$ kcal / mol @f$
 *  @param  options     The options flag indicating how/where to store the soft constraints
 *  @return             Non-zero on successful application of the constraint, 0 otherwise.
 */
int
vrna_sc_set_up(vrna_fold_compound_t *fc,
               const FLT_OR_DBL     *constraints,
               unsigned int         options);


/**
 *  @brief  Set soft constraints for unpaired nucleotides in comparative structure predictions
 *
 *  Use this function to set soft constraints for unpaired nucleotides for each sequence
 *  in the multiple sequence alignment (MSA). The constraints are provided as 0-based array
 *  of 1-based array with the actual pseudo energies, where the first dimension corresponds
 *  to the number (0-based) of the respective sequence in the alignment. If no constraints
 *  are provided for a particular sequence @f$ s @f$ in the MSA, the corresponding entry
 *  must be set to @b NULL.
 *
 *  @note     This function replaces any pre-exisitng soft constraints with the ones supplied
 *            in @p constraints.
 *
 *  @warning  Pseudo energies for all sequences must be provided in sequence coordinates
 *            rather than alignment coordinates!
 *
 *  @see      vrna_sc_set_up_comparative_seq(), vrna_sc_add_up_comparative(),
 *            vrna_sc_set_bp_comparative(), vrna_sc_add_stack_comparative(),
 *            vrna_sc_set_stack_comparative(),
 *            vrna_sc_set_up()
 *
 *  @param  fc          The #vrna_fold_compound_t the soft constraints are associated with
 *  @param  constraints A 0-based array of 1-based arrays with pseudo free energies in @f$ kcal / mol @f$
 *  @param  options     The options flag indicating how/where to store the soft constraints
 *  @return             The number of sequences in the MSA constraints have been applied to
 */
int
vrna_sc_set_up_comparative(vrna_fold_compound_t *fc,
                           const FLT_OR_DBL     **constraints,
                           unsigned int         options);


/**
 *  @brief  Set soft constraints for unpaired nucleotides in comparative structure predictions
 *
 *  This is a convenience wrapper for vrna_sc_set_up_comparative() where only one particular
 *  sequence @p s is provided with constraints.
 *
 *  @note     This function replaces any pre-exisitng soft constraints with the ones supplied
 *            in @p constraints.
 *
 *  @warning  Pseudo energies for all sequences must be provided in sequence coordinates
 *            rather than alignment coordinates!
 *
 *  @see      vrna_sc_set_up_comparative(), vrna_sc_add_up_comparative_seq(),
 *            vrna_sc_set_bp_comparative_seq(), vrna_sc_add_stack_comparative_seq(),
 *            vrna_sc_set_stack_comparative_seq(),
 *            vrna_sc_set_up()
 *
 *  @param  fc          The #vrna_fold_compound_t the soft constraints are associated with
 *  @param  s           The 0-based number of the sequence in the alignment the constraints are provided for
 *  @param  constraints A 1-based arrays with pseudo free energies in @f$ kcal / mol @f$
 *  @param  options     The options flag indicating how/where to store the soft constraints
 *  @return             The number of sequences in the MSA constraints have been applied to
 */
int
vrna_sc_set_up_comparative_seq(vrna_fold_compound_t *fc,
                               unsigned int         s,
                               const FLT_OR_DBL     *constraints,
                               unsigned int         options);


/**
 *  @brief  Add soft constraints for unpaired nucleotides
 *
 *  @see      vrna_sc_set_up(), vrna_sc_add_bp(), vrna_sc_set_bp(),
 *            vrna_sc_add_up_comparative()
 *
 *  @param  fc          The #vrna_fold_compound_t the soft constraints are associated with
 *  @param  i           The nucleotide position the soft constraint is added for
 *  @param  energy      The free energy (soft-constraint) in @f$ kcal / mol @f$
 *  @param  options     The options flag indicating how/where to store the soft constraints
 *  @return             Non-zero on successful application of the constraint, 0 otherwise.
 */
int
vrna_sc_add_up(vrna_fold_compound_t *fc,
               unsigned int         i,
               FLT_OR_DBL           energy,
               unsigned int         options);


/**
 *  @brief  Add soft constraints for unpaired nucleotides in comparative structure predictions
 *
 *  Use this function to add soft constraints for unpaired nucleotides for each sequence
 *  in the multiple sequence alignment (MSA). The constraints are provided as 0-based array
 *  of pseudo energies, one for each sequence in the MSA. If no constraints are provided for
 *  a particular sequence, the corrsponding @p is value must be 0.
 *
 *  @warning  Pseudo energies for all sequences must be provided in sequence coordinates
 *            rather than alignment coordinates!
 *
 *  @see      vrna_sc_set_up_comparative_seq(), vrna_sc_add_up_comparative(),
 *            vrna_sc_set_bp_comparative(), vrna_sc_add_stack_comparative(),
 *            vrna_sc_set_stack_comparative(),
 *            vrna_sc_add_up()
 *
 *  @param  fc          The #vrna_fold_compound_t the soft constraints are associated with
 *  @param  is          A 0-based array of nucleotide positions the soft constraint is added for
 *  @param  energies    A 0-based array of free energies (soft-constraint) in @f$ kcal / mol @f$
 *  @param  options     The options flag indicating how/where to store the soft constraints
 *  @return             Non-zero on successful application of the constraint, 0 otherwise.
 */
int
vrna_sc_add_up_comparative(vrna_fold_compound_t *fc,
                           unsigned int         *is,
                           const FLT_OR_DBL     *energies,
                           unsigned int         options);


/**
 *  @brief  Add soft constraints for unpaired nucleotides in comparative structure predictions
 *
 *  This is a convenience wrapper for vrna_sc_add_up_comparative() where only one particular
 *  sequence @p s is provided with constraints.
 *
 *  @warning  Pseudo energies for all sequences must be provided in sequence coordinates
 *            rather than alignment coordinates!
 *
 *  @see      vrna_sc_set_up_comparative_seq(), vrna_sc_add_up_comparative(),
 *            vrna_sc_set_bp_comparative(), vrna_sc_add_stack_comparative(),
 *            vrna_sc_set_stack_comparative(),
 *            vrna_sc_add_up()
 *
 *  @param  fc          The #vrna_fold_compound_t the soft constraints are associated with
 *  @param  s           The 0-based number of the sequence in the alignment the constraints are provided for
 *  @param  i           The 1-based nucleotide position the soft constraint is added for
 *  @param  energy      The free energies (soft-constraint) in @f$ kcal / mol @f$
 *  @param  options     The options flag indicating how/where to store the soft constraints
 *  @return             Non-zero on successful application of the constraint, 0 otherwise.
 */
int
vrna_sc_add_up_comparative_seq(vrna_fold_compound_t *fc,
                               unsigned int         s,
                               unsigned int         i,
                               const FLT_OR_DBL     energy,
                               unsigned int         options);


/**
 *  @}
 *  @addtogroup soft_constraints_st
 *  @{
 */

int
vrna_sc_set_stack(vrna_fold_compound_t  *fc,
                  const FLT_OR_DBL      *constraints,
                  unsigned int          options);


int
vrna_sc_set_stack_comparative(vrna_fold_compound_t  *fc,
                              const FLT_OR_DBL      **constraints,
                              unsigned int          options);


int
vrna_sc_set_stack_comparative_seq(vrna_fold_compound_t  *fc,
                                  unsigned int          s,
                                  const FLT_OR_DBL      *constraints,
                                  unsigned int          options);

int
vrna_sc_add_stack(vrna_fold_compound_t  *fc,
                  unsigned int          i,
                  FLT_OR_DBL            energy,
                  unsigned int          options);


int
vrna_sc_add_stack_comparative(vrna_fold_compound_t  *fc,
                              unsigned int          *is,
                              const FLT_OR_DBL      *energies,
                              unsigned int          options);


int
vrna_sc_add_stack_comparative_seq(vrna_fold_compound_t  *fc,
                                  unsigned int          s,
                                  unsigned int          i,
                                  FLT_OR_DBL            energy,
                                  unsigned int          options);


/**
 *  @}
 *  @addtogroup soft_constraints_generic
 *  @{
 */


/**
 *  @brief Add an auxiliary data structure for the generic soft constraints callback function
 *
 *  @see vrna_sc_add_f(), vrna_sc_add_exp_f(), vrna_sc_add_bt()
 *
 *  @param  fc        The fold compound the generic soft constraint function should be bound to
 *  @param  data      A pointer to the data structure that holds required data for function 'f'
 *  @param  free_data A pointer to a function that free's the memory occupied by @p data (Maybe NULL)
 *  @return           Non-zero on successful binding the data (and free-function), 0 otherwise
 */
int
vrna_sc_add_data(vrna_fold_compound_t *fc,
                 void                 *data,
                 vrna_auxdata_free_f  free_data);


int
vrna_sc_add_auxdata(vrna_fold_compound_t    *fc,
                    void                    *data,
                    vrna_auxdata_prepare_f  prepare_cb,
                    vrna_auxdata_free_f     free_cb);


int
vrna_sc_set_data_comparative(vrna_fold_compound_t *fc,
                             void                 **data,
                             vrna_auxdata_free_f  *free_data,
                             unsigned int         options);


int
vrna_sc_set_data_comparative_seq(vrna_fold_compound_t *fc,
                                 unsigned int         s,
                                 void                 *data,
                                 vrna_auxdata_free_f  free_data,
                                 unsigned int         options);


int
vrna_sc_set_auxdata_comparative(vrna_fold_compound_t    *fc,
                                void                    **data,
                                vrna_auxdata_prepare_f  *prepare_cbs,
                                vrna_auxdata_free_f     *free_data,
                                unsigned int            options);


int
vrna_sc_set_auxdata_comparative_seq(vrna_fold_compound_t    *fc,
                                    unsigned int            s,
                                    void                    *data,
                                    vrna_auxdata_prepare_f  prepare_cb,
                                    vrna_auxdata_free_f     free_data,
                                    unsigned int            options);


/**
 *  @brief  Bind a function pointer for generic soft constraint feature (MFE version)
 *
 *  This function allows one to easily bind a function pointer and corresponding data structure
 *  to the soft constraint part #vrna_sc_t of the #vrna_fold_compound_t.
 *  The function for evaluating the generic soft constraint feature has to return
 *  a pseudo free energy @f$ \hat{E} @f$ in @f$ dacal/mol @f$, where @f$ 1 dacal/mol = 10 cal/mol @f$.
 *
 *  @see vrna_sc_add_data(), vrna_sc_add_bt(), vrna_sc_add_exp_f()
 *
 *  @param  fc    The fold compound the generic soft constraint function should be bound to
 *  @param  f     A pointer to the function that evaluates the generic soft constraint feature
 *  @return       Non-zero on successful binding the callback function, 0 otherwise
 */
int
vrna_sc_add_f(vrna_fold_compound_t  *fc,
              vrna_sc_f             f);


size_t
vrna_sc_multi_cb_add(vrna_fold_compound_t   *fc,
                     vrna_sc_direct_f       cb,
                     vrna_sc_exp_direct_f   cb_exp,
                     void                   *data,
                     vrna_auxdata_prepare_f prepare_cb,
                     vrna_auxdata_free_f    free_cb,
                     unsigned int           decomp_type);


size_t
vrna_sc_multi_cb_add_comparative(vrna_fold_compound_t   *fc,
                                 vrna_sc_direct_f       *cbs,
                                 vrna_sc_exp_direct_f   *cbs_exp,
                                 void                   **datas,
                                 vrna_auxdata_prepare_f *prepare_cbs,
                                 vrna_auxdata_free_f    *free_cbs,
                                 unsigned int           *ds,
                                 unsigned int           multi_params);

int
vrna_sc_set_f_comparative(vrna_fold_compound_t  *fc,
                          vrna_sc_f             *f,
                          unsigned int          options);


/**
 *  @brief  Bind a backtracking function pointer for generic soft constraint feature
 *
 *  This function allows one to easily bind a function pointer to the soft constraint part
 *  #vrna_sc_t of the #vrna_fold_compound_t.
 *  The provided function should be used for backtracking purposes in loop regions
 *  that were altered via the generic soft constraint feature. It has to return
 *  an array of #vrna_basepair_t data structures, were the last element in the list is indicated
 *  by a value of -1 in it's i position.
 *
 *  @see vrna_sc_add_data(), vrna_sc_add_f(), vrna_sc_add_exp_f()
 *
 *  @param  fc    The fold compound the generic soft constraint function should be bound to
 *  @param  f     A pointer to the function that returns additional base pairs
 *  @return       Non-zero on successful binding the callback function, 0 otherwise
 */
int
vrna_sc_add_bt(vrna_fold_compound_t *fc,
               vrna_sc_bt_f         f);


/**
 *  @brief  Bind a function pointer for generic soft constraint feature (PF version)
 *
 *  This function allows one to easily bind a function pointer and corresponding data structure
 *  to the soft constraint part #vrna_sc_t of the #vrna_fold_compound_t.
 *  The function for evaluating the generic soft constraint feature has to return
 *  a pseudo free energy @f$ \hat{E} @f$ as Boltzmann factor, i.e. @f$ exp(- \hat{E} / kT) @f$.
 *  The required unit for @f$ E @f$ is @f$ cal/mol @f$.
 *
 *  @see vrna_sc_add_bt(), vrna_sc_add_f(), vrna_sc_add_data()
 *
 *  @param  fc    The fold compound the generic soft constraint function should be bound to
 *  @param  exp_f A pointer to the function that evaluates the generic soft constraint feature
 *  @return       Non-zero on successful binding the callback function, 0 otherwise
 */
int
vrna_sc_add_exp_f(vrna_fold_compound_t  *fc,
                  vrna_sc_exp_f         exp_f);


int
vrna_sc_set_exp_f_comparative(vrna_fold_compound_t  *fc,
                              vrna_sc_exp_f         *exp_f,
                              unsigned int          options);

/**
 *  @}
 */

#endif
