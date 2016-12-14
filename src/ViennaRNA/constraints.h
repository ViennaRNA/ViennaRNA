#ifndef VIENNA_RNA_PACKAGE_CONSTRAINTS_H
#define VIENNA_RNA_PACKAGE_CONSTRAINTS_H

#include <ViennaRNA/data_structures.h>

/* include all structure constraint related headers */
#include <ViennaRNA/constraints_hard.h>
#include <ViennaRNA/constraints_soft.h>
#include <ViennaRNA/constraints_SHAPE.h>
#include <ViennaRNA/perturbation_fold.h>
#include <ViennaRNA/constraints_ligand.h>

/**
 *  @file     constraints.h
 *  @brief    Functions and data structures for constraining secondary structure predictions and evaluation
 *  @ingroup  constraints
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
#define VRNA_DECOMP_PAIR_HP     1

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
#define VRNA_DECOMP_PAIR_IL     2

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
#define VRNA_DECOMP_PAIR_ML     3

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
#define VRNA_DECOMP_ML_ML_ML    5

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
#define VRNA_DECOMP_ML_STEM     4

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
#define VRNA_DECOMP_ML_ML       6

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
#define VRNA_DECOMP_ML_UP       11

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
#define VRNA_DECOMP_ML_ML_STEM 20

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
#define VRNA_DECOMP_ML_COAXIAL  13

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
#define VRNA_DECOMP_ML_COAXIAL_ENC  22

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
#define VRNA_DECOMP_EXT_EXT     9

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
#define VRNA_DECOMP_EXT_UP      8

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
#define VRNA_DECOMP_EXT_STEM 14

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
#define VRNA_DECOMP_EXT_EXT_EXT 15

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
#define VRNA_DECOMP_EXT_STEM_EXT 16

/**
 *  @brief  Indicator for decomposition of exterior loop part
 *
 *  @ingroup  constraints
 *
 */
#define VRNA_DECOMP_EXT_STEM_OUTSIDE 17

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
#define VRNA_DECOMP_EXT_EXT_STEM 18

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

 */
#define VRNA_DECOMP_EXT_EXT_STEM1 19


#define VRNA_DECOMP_EXT_L         20


#define VRNA_DECOMP_EXT_EXT_L     21

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
void vrna_constraints_add(vrna_fold_compound_t *vc,
                          const char *constraint,
                          unsigned int options);

#endif
