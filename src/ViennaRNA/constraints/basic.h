#ifndef VIENNA_RNA_PACKAGE_CONSTRAINTS_H
#define VIENNA_RNA_PACKAGE_CONSTRAINTS_H

#include <ViennaRNA/fold_compound.h>

/**
 *  @file     constraints/basic.h
 *  @ingroup  constraints
 *  @brief    Functions and data structures for constraining secondary structure predictions and evaluation
 */

/**
 *  @addtogroup constraints
 *
 *  @brief This module provides general functions that allow for an easy control of
 *  constrained secondary structure prediction and evaluation.
 *
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
 *  Secondary structure constraints are always applied at decomposition level, i.e.
 *  in each step of the recursive structure decomposition, for instance during MFE
 *  prediction. Below is a visualization of the decomposition scheme
 *
 *  @image html   recursions.svg
 *  @image latex  recursions.eps
 *
 *  For @ref hard_constraints the following option flags may be used to constrain
 *  the pairing behavior of single, or pairs of nucleotides:
 *
 *  - #VRNA_CONSTRAINT_CONTEXT_EXT_LOOP - @copybrief #VRNA_CONSTRAINT_CONTEXT_EXT_LOOP
 *  - #VRNA_CONSTRAINT_CONTEXT_HP_LOOP  - @copybrief #VRNA_CONSTRAINT_CONTEXT_HP_LOOP
 *  - #VRNA_CONSTRAINT_CONTEXT_INT_LOOP  - @copybrief #VRNA_CONSTRAINT_CONTEXT_INT_LOOP
 *  - #VRNA_CONSTRAINT_CONTEXT_INT_LOOP_ENC - @copybrief #VRNA_CONSTRAINT_CONTEXT_INT_LOOP_ENC
 *  - #VRNA_CONSTRAINT_CONTEXT_MB_LOOP - @copybrief #VRNA_CONSTRAINT_CONTEXT_MB_LOOP
 *  - #VRNA_CONSTRAINT_CONTEXT_MB_LOOP_ENC - @copybrief #VRNA_CONSTRAINT_CONTEXT_MB_LOOP_ENC
 *  - #VRNA_CONSTRAINT_CONTEXT_ENFORCE - @copybrief #VRNA_CONSTRAINT_CONTEXT_ENFORCE
 *  - #VRNA_CONSTRAINT_CONTEXT_NO_REMOVE - @copybrief #VRNA_CONSTRAINT_CONTEXT_NO_REMOVE
 *  - #VRNA_CONSTRAINT_CONTEXT_ALL_LOOPS - @copybrief #VRNA_CONSTRAINT_CONTEXT_ALL_LOOPS
 *
 *  However, for @ref soft_constraints we do not allow for simple loop type
 *  dependent constraining. But soft constraints are equipped with generic constraint
 *  support. This enables the user to pass arbitrary callback functions that
 *  return auxiliary energy contributions for evaluation the evaluation of any
 *  decomposition.
 *
 *  The callback will then always be notified about the type of decomposition
 *  that is happening, and the corresponding delimiting sequence positions. The
 *  following decomposition steps are distinguished, and should be captured by
 *  the user's implementation of the callback:
 *
 *  - #VRNA_DECOMP_PAIR_HP - @copybrief #VRNA_DECOMP_PAIR_HP
 *  - #VRNA_DECOMP_PAIR_IL - @copybrief #VRNA_DECOMP_PAIR_IL
 *  - #VRNA_DECOMP_PAIR_ML - @copybrief #VRNA_DECOMP_PAIR_ML
 *  - #VRNA_DECOMP_ML_ML_ML - @copybrief #VRNA_DECOMP_ML_ML_ML
 *  - #VRNA_DECOMP_ML_STEM - @copybrief #VRNA_DECOMP_ML_STEM
 *  - #VRNA_DECOMP_ML_ML - @copybrief #VRNA_DECOMP_ML_ML
 *  - #VRNA_DECOMP_ML_UP - @copybrief #VRNA_DECOMP_ML_UP
 *  - #VRNA_DECOMP_ML_ML_STEM - @copybrief #VRNA_DECOMP_ML_ML_STEM
 *  - #VRNA_DECOMP_ML_COAXIAL - @copybrief #VRNA_DECOMP_ML_COAXIAL
 *  - #VRNA_DECOMP_EXT_EXT - @copybrief #VRNA_DECOMP_EXT_EXT
 *  - #VRNA_DECOMP_EXT_UP - @copybrief #VRNA_DECOMP_EXT_UP
 *  - #VRNA_DECOMP_EXT_STEM - @copybrief #VRNA_DECOMP_EXT_STEM
 *  - #VRNA_DECOMP_EXT_EXT_EXT - @copybrief #VRNA_DECOMP_EXT_EXT_EXT
 *  - #VRNA_DECOMP_EXT_STEM_EXT - @copybrief #VRNA_DECOMP_EXT_STEM_EXT
 *  - #VRNA_DECOMP_EXT_STEM_OUTSIDE - @copybrief #VRNA_DECOMP_EXT_STEM_OUTSIDE
 *  - #VRNA_DECOMP_EXT_EXT_STEM - @copybrief #VRNA_DECOMP_EXT_EXT_STEM
 *  - #VRNA_DECOMP_EXT_EXT_STEM1 - @copybrief #VRNA_DECOMP_EXT_EXT_STEM1
 *
 *  Simplified interfaces to the soft constraints framework can be obtained
 *  by the implementations in the submodules
 *
 *  - @ref SHAPE_reactivities and
 *  - @ref ligands.
 *
 *  An implementation that generates soft constraints for unpaired nucleotides
 *  by minimizing the discrepancy between their predicted and expected pairing
 *  probability is available in submodule @ref perturbation.
 *
 */


/**
 *  @brief  Flag for vrna_constraints_add() to indicate that constraints are present in a text file
 *
 *  @see vrna_constraints_add()
 *  @deprecated   Use 0 instead!
 *  @ingroup  constraints
 *
 */
#define VRNA_CONSTRAINT_FILE      0

/**
 *  @brief  Indicate generation of constraints for MFE folding
 *  @deprecated   This flag has no meaning anymore, since constraints are now always stored!
 *  @ingroup  constraints
 *
 */
#define VRNA_CONSTRAINT_SOFT_MFE  0

/**
 *  @brief  Indicate generation of constraints for partition function computation
 *  @deprecated   Use #VRNA_OPTION_PF instead!
 *  @ingroup  constraints
 *
 */
#define VRNA_CONSTRAINT_SOFT_PF   VRNA_OPTION_PF

/**
 *  @brief  Flag passed to generic softt constraints callback to indicate hairpin loop decomposition step
 *
 *  @ingroup  constraints
 *
 *  @details This flag notifies the soft or hard constraint callback function that the current
 *  decomposition step evaluates a hairpin loop enclosed by the base pair @f$(i,j)@f$.
 *
 *  @image html   decomp_hp.svg
 *  @image latex  decomp_hp.eps
 *
 */
#define VRNA_DECOMP_PAIR_HP     (unsigned char)1

/**
 *  @brief  Indicator for interior loop decomposition step
 *
 *  @ingroup  constraints
 *
 *  @details This flag notifies the soft or hard constraint callback function that the current
 *  decomposition step evaluates an interior loop enclosed by the base pair @f$(i,j)@f$,
 *  and enclosing the base pair @f$(k,l)@f$.
 *
 *  @image html   decomp_il.svg
 *  @image latex  decomp_il.eps
 *
 */
#define VRNA_DECOMP_PAIR_IL     (unsigned char)2

/**
 *  @brief  Indicator for multibranch loop decomposition step
 *
 *  @ingroup  constraints
 *
 *  @details This flag notifies the soft or hard constraint callback function that the current
 *  decomposition step evaluates a multibranch loop enclosed by the base pair @f$(i,j)@f$,
 *  and consisting of some enclosed multi loop content from k to l.
 *
 *  @image html   decomp_ml.svg
 *  @image latex  decomp_ml.eps
 *
 */
#define VRNA_DECOMP_PAIR_ML     (unsigned char)3
#define VRNA_DECOMP_PAIR_ML_EXT     (unsigned char)23

#define VRNA_DECOMP_PAIR_ML_OUTSIDE     (unsigned char)4
/**
 *  @brief  Indicator for decomposition of multibranch loop part
 *
 *  @ingroup  constraints
 *
 *  @details This flag notifies the soft or hard constraint callback function that the current
 *  decomposition step evaluates a multibranch loop part in the interval @f$[i:j]@f$,
 *  which will be decomposed into two multibranch loop parts @f$[i:k]@f$, and @f$[l:j]@f$.
 *
 *  @image html   decomp_ml_ml_ml.svg
 *  @image latex  decomp_ml_ml_ml.eps
 *
 */
#define VRNA_DECOMP_ML_ML_ML    (unsigned char)5

/**
 *  @brief  Indicator for decomposition of multibranch loop part
 *
 *  @ingroup  constraints
 *
 *  @details This flag notifies the soft or hard constraint callback function that the current
 *  decomposition step evaluates a multibranch loop part in the interval @f$[i:j]@f$,
 *  which will be considered a single stem branching off with base pair @f$(k,l)@f$.
 *
 *  @image html   decomp_ml_stem.svg
 *  @image latex  decomp_ml_stem.eps
 *
 */
#define VRNA_DECOMP_ML_STEM     (unsigned char)6

/**
 *  @brief  Indicator for decomposition of multibranch loop part
 *
 *  @ingroup  constraints
 *
 *  @details This flag notifies the soft or hard constraint callback function that the current
 *  decomposition step evaluates a multibranch loop part in the interval @f$[i:j]@f$,
 *  which will be decomposed into a (usually) smaller multibranch loop part @f$[k:l]@f$.
 *
 *  @image html   decomp_ml_ml.svg
 *  @image latex  decomp_ml_ml.eps
 *
 */
#define VRNA_DECOMP_ML_ML       (unsigned char)7

/**
 *  @brief  Indicator for decomposition of multibranch loop part
 *
 *  @ingroup  constraints
 *
 *  @details This flag notifies the soft or hard constraint callback function that the current
 *  decomposition step evaluates a multibranch loop part in the interval @f$[i:j]@f$,
 *  which will be considered a multibranch loop part that only consists of unpaired
 *  nucleotides.
 *
 *  @image html   decomp_ml_up.svg
 *  @image latex  decomp_ml_up.eps
 *
 */
#define VRNA_DECOMP_ML_UP       (unsigned char)8

/**
 *  @brief  Indicator for decomposition of multibranch loop part
 *
 *  @ingroup  constraints
 *
 *  @details This flag notifies the soft or hard constraint callback function that the current
 *  decomposition step evaluates a multibranch loop part in the interval @f$[i:j]@f$,
 *  which will decomposed into a multibranch loop part @f$[i:k]@f$, and a stem with
 *  enclosing base pair @f$(l,j)@f$.
 *
 *  @image html   decomp_ml_ml_stem.svg
 *  @image latex  decomp_ml_ml_stem.eps
 *
 */
#define VRNA_DECOMP_ML_ML_STEM (unsigned char)9

/**
 *  @brief  Indicator for decomposition of multibranch loop part
 *
 *  @ingroup  constraints
 *
 *  @details This flag notifies the soft or hard constraint callback function that the current
 *  decomposition step evaluates a multibranch loop part in the interval @f$[i:j]@f$,
 *  where two stems with enclosing pairs @f$(i,k)@f$ and @f$(l,j)@f$ are coaxially stacking
 *  onto each other.
 *
 *  @image html   decomp_ml_coaxial.svg
 *  @image latex  decomp_ml_coaxial.eps
 *
 */
#define VRNA_DECOMP_ML_COAXIAL  (unsigned char)10

/**
 *  @brief  Indicator for decomposition of multibranch loop part
 *
 *  @ingroup  constraints
 *
 *  @details This flag notifies the soft or hard constraint callback function that the current
 *  decomposition step evaluates a multibranch loop part in the interval @f$[i:j]@f$,
 *  where two stems with enclosing pairs @f$(i,k)@f$ and @f$(l,j)@f$ are coaxially stacking
 *  onto each other.
 *
 *  @image html   decomp_ml_coaxial.svg
 *  @image latex  decomp_ml_coaxial.eps
 *
 */
#define VRNA_DECOMP_ML_COAXIAL_ENC  (unsigned char)11

/**
 *  @brief  Indicator for decomposition of exterior loop part
 *
 *  @ingroup  constraints
 *
 *  @def VRNA_DECOMP_EXT_EXT
 *  @details This flag notifies the soft or hard constraint callback function that the current
 *  decomposition step evaluates an exterior loop part in the interval @f$[i:j]@f$,
 *  which will be decomposed into a (usually) smaller exterior loop part @f$[k:l]@f$.
 *
 *  @image html   decomp_ext_ext.svg
 *  @image latex  decomp_ext_ext.eps
 *
 */
#define VRNA_DECOMP_EXT_EXT     (unsigned char)12

/**
 *  @brief  Indicator for decomposition of exterior loop part
 *
 *  @ingroup  constraints
 *
 *  @details This flag notifies the soft or hard constraint callback function that the current
 *  decomposition step evaluates an exterior loop part in the interval @f$[i:j]@f$,
 *  which will be considered as an exterior loop component consisting of only unpaired
 *  nucleotides.
 *
 *  @image html   decomp_ext_up.svg
 *  @image latex  decomp_ext_up.eps
 *
 */
#define VRNA_DECOMP_EXT_UP      (unsigned char)13

/**
 *  @brief  Indicator for decomposition of exterior loop part
 *
 *  @ingroup  constraints
 *
 *  @details This flag notifies the soft or hard constraint callback function that the current
 *  decomposition step evaluates an exterior loop part in the interval @f$[i:j]@f$,
 *  which will be considered a stem with enclosing pair @f$(k,l)@f$.
 *
 *  @image html   decomp_ext_stem.svg
 *  @image latex  decomp_ext_stem.eps
 *
 */
#define VRNA_DECOMP_EXT_STEM (unsigned char)14

/**
 *  @brief  Indicator for decomposition of exterior loop part
 *
 *  @ingroup  constraints
 *
 *  @details This flag notifies the soft or hard constraint callback function that the current
 *  decomposition step evaluates an exterior loop part in the interval @f$[i:j]@f$,
 *  which will be decomposed into two exterior loop parts @f$[i:k]@f$ and @f$[l:j]@f$.
 *
 *  @image html   decomp_ext_ext_ext.svg
 *  @image latex  decomp_ext_ext_ext.eps
 *
 */
#define VRNA_DECOMP_EXT_EXT_EXT (unsigned char)15

/**
 *  @brief  Indicator for decomposition of exterior loop part
 *
 *  @ingroup  constraints
 *
 *  @details This flag notifies the soft or hard constraint callback function that the current
 *  decomposition step evaluates an exterior loop part in the interval @f$[i:j]@f$,
 *  which will be decomposed into a stem branching off with base pair @f$(i,k)@f$, and
 *  an exterior loop part @f$[l:j]@f$.
 *
 *  @image html   decomp_ext_stem_ext.svg
 *  @image latex  decomp_ext_stem_ext.eps
 *
 */
#define VRNA_DECOMP_EXT_STEM_EXT (unsigned char)16

/**
 *  @brief  Indicator for decomposition of exterior loop part
 *
 *  @ingroup  constraints
 *
 */
#define VRNA_DECOMP_EXT_STEM_OUTSIDE (unsigned char)17

/**
 *  @brief  Indicator for decomposition of exterior loop part
 *
 *  @ingroup  constraints
 *
 *  @details This flag notifies the soft or hard constraint callback function that the current
 *  decomposition step evaluates an exterior loop part in the interval @f$[i:j]@f$,
 *  which will be decomposed into an exterior loop part @f$[i:k]@f$, and a stem
 *  branching off with base pair @f$(l,j)@f$.
 *
 *  @image html   decomp_ext_ext_stem.svg
 *  @image latex  decomp_ext_ext_stem.eps
 *
 */
#define VRNA_DECOMP_EXT_EXT_STEM (unsigned char)18

/**
 *  @brief  Indicator for decomposition of exterior loop part
 *
 *  @ingroup  constraints
 *
 *  @def VRNA_DECOMP_EXT_EXT_STEM1
 *  @details This flag notifies the soft or hard constraint callback function that the current
 *  decomposition step evaluates an exterior loop part in the interval @f$[i:j]@f$,
 *  which will be decomposed into an exterior loop part @f$[i:k]@f$, and a stem
 *  branching off with base pair @f$(l,j-1)@f$.
 *
 *  @image html   decomp_ext_ext_stem1.svg
 *  @image latex  decomp_ext_ext_stem1.eps
 *
 */
#define VRNA_DECOMP_EXT_EXT_STEM1 (unsigned char)19

#define VRNA_DECOMP_EXT_STEM_EXT1 (unsigned char)20

#define VRNA_DECOMP_EXT_L         (unsigned char)21
#define VRNA_DECOMP_EXT_EXT_L     (unsigned char)22

/**
 *  @brief  Add constraints to a #vrna_fold_compound_t data structure
 *
 *  Use this function to add/update the hard/soft constraints
 *  The function allows for passing a string 'constraint' that can either be a
 *  filename that points to a constraints definition file or it may be a
 *  pseudo dot-bracket notation indicating hard constraints. For the latter, the
 *  user has to pass the #VRNA_CONSTRAINT_DB option. Also, the
 *  user has to specify, which characters are allowed to be interpreted as
 *  constraints by passing the corresponding options via the third parameter.
 *
 *  @see      vrna_hc_init(), vrna_hc_add_up(), vrna_hc_add_up_batch(), vrna_hc_add_bp(),
 *            vrna_sc_init(), vrna_sc_set_up(), vrna_sc_set_bp(),
 *            vrna_sc_add_SHAPE_deigan(),  vrna_sc_add_SHAPE_zarringhalam(),
 *            vrna_hc_free(), vrna_sc_free(),
 *            #VRNA_CONSTRAINT_DB, #VRNA_CONSTRAINT_DB_DEFAULT, #VRNA_CONSTRAINT_DB_PIPE,
 *            #VRNA_CONSTRAINT_DB_DOT, #VRNA_CONSTRAINT_DB_X, #VRNA_CONSTRAINT_DB_ANG_BRACK,
 *            #VRNA_CONSTRAINT_DB_RND_BRACK, #VRNA_CONSTRAINT_DB_INTRAMOL,
 *            #VRNA_CONSTRAINT_DB_INTERMOL, #VRNA_CONSTRAINT_DB_GQUAD
 *
 *  @ingroup  constraints
 *
 *  The following is an example for adding hard constraints given in
 *  pseudo dot-bracket notation. Here, @p vc is the #vrna_fold_compound_t object,
 *  @p structure is a char array with the hard constraint in dot-bracket notation,
 *  and @p enforceConstraints is a flag indicating whether or not constraints for
 *  base pairs should be enforced instead of just doing a removal of base pair that
 *  conflict with the constraint.
 *
 *  @snippet RNAfold.c Adding hard constraints from pseudo dot-bracket
 *
 *  In constrat to the above, constraints may also be read from file:
 *
 *  @snippet RNAfold.c Adding hard constraints from file
 *
 *  @see  vrna_hc_add_from_db(), vrna_hc_add_up(), vrna_hc_add_up_batch()
 *        vrna_hc_add_bp_unspecific(), vrna_hc_add_bp()
 *
 *  @param  vc            The fold compound
 *  @param  constraint    A string with either the filename of the constraint definitions
 *                        or a pseudo dot-bracket notation of the hard constraint. May be NULL.
 *  @param  options       The option flags
 */
void vrna_constraints_add(vrna_fold_compound_t  *vc,
                          const char            *constraint,
                          unsigned int          options);


#endif
