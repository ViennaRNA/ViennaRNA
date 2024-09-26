#ifndef VIENNA_RNA_PACKAGE_GRAMMAR_MFE_H
#define VIENNA_RNA_PACKAGE_GRAMMAR_MFE_H

/**
 *  @file     ViennaRNA/grammar/mfe.h
 *  @ingroup  grammar
 *  @brief    Implementations for the RNA folding grammar
 */

/**
 * @addtogroup grammar
 * @{
 */

#include <ViennaRNA/fold_compound.h>
#include <ViennaRNA/datastructures/basic.h>
#include <ViennaRNA/datastructures/array.h>
#include <ViennaRNA/grammar/basic.h>

/**
 *  @brief  Function prototype for auxiliary grammar rules (inside version, MFE)
 *
 *  This function will be called during the inside recursions of MFE predictions
 *  for subsequences from position @p i to @p j and is supposed to return an
 *  energy contribution in dekacal/mol.
 *
 *  @callback
 *  @parblock
 *  This callback allows for extending the MFE secondary structure decomposition with
 *  additional rules.
 *  @endparblock
 *
 *  @see  vrna_gr_add_aux_f(), vrna_gr_add_aux_c(), vrna_gr_add_aux_m(),
 *        vrna_gr_add_aux_m1(), vrna_gr_add_aux_aux(), vrna_gr_add_aux_exp_f()
 *
 *  @param  fc    The fold compound to work on
 *  @param  i     The 5' delimiter of the sequence segment
 *  @param  j     The 3' delimiter of the sequence segment
 *  @param  data  An arbitrary user-provided data pointer
 *  @return       The free energy computed by the auxiliary grammar rule in dekacal/mol
 */
typedef int (*vrna_gr_inside_f)(vrna_fold_compound_t  *fc,
                                unsigned int          i,
                                unsigned int          j,
                                void                  *data);


/**
 *  @brief  Function prototype for auxiliary grammar rules (outside version, MFE)
 *
 *  This function will be called during the backtracking stage of MFE predictions
 *  for subsequences from position @p i to @p j and is supposed to identify the
 *  the structure components that give rise to the corresponding energy contribution
 *  @p e as determined in the inside (forward) step.
 *
 *  For that purpose, the function may modify the base pair stack (@p bp_stack) by adding
 *  all base pairs identified through the additional rule. In addition, when the rule
 *  sub-divides the sequence segment into smaller pieces, these pieces can be put onto
 *  the backtracking (@p bt_stack) stack.
 *  On successful backtrack, the function returns non-zero, and exactly zero (0)
 *  when the backtracking failed.
 *
 *  @callback
 *  @parblock
 *  This callback allows for backtracking (sub)structures obtained from extending the
 *  MFE secondary structure decomposition with additional rules.
 *  @endparblock
 *
 *  @see  vrna_gr_add_aux_f(), vrna_gr_add_aux_c(), vrna_gr_add_aux_m(),
 *        vrna_gr_add_aux_m1(), vrna_gr_add_aux_aux(), vrna_gr_add_aux_exp_f()
 *
 *  @param  fc            The fold compound to work on
 *  @param  i             The 5' delimiter of the sequence segment
 *  @param  j             The 3' delimiter of the sequence segment
 *  @param  e             The energy of the segment [i:j]
 *  @param  bp_stack      The base pair stack
 *  @param  bt_stack      The backtracking stack (all segments that need to be further investigated)
 *  @param  data  An arbitrary user-provided data pointer
 *  @return       The free energy computed by the auxiliary grammar rule in dekacal/mol
 */
typedef int (*vrna_gr_outside_f)(vrna_fold_compound_t *fc,
                                 unsigned int         i,
                                 unsigned int         j,
                                 int                  e,
                                 vrna_bps_t           bp_stack,
                                 vrna_bts_t           bt_stack,
                                 void                 *data);


/**
 *  @brief  Add an auxiliary grammar rule for the F-decomposition (MFE version)
 *
 *  This function binds callback functions for auxiliary grammar rules (inside and outside)
 *  in the F-decomposition, i.e. the external loop decomposition stage.
 *
 *  While the inside rule (@p cb) computes a minimum free energy contribution for any
 *  subsequence the outside rule (@p cb_bt) is used for backtracking the corresponding
 *  structure.
 *
 *  Both callbacks will be provided with the @p data pointer that can be used to
 *  store whatever data is needed in the callback evaluations. The @p data_prepare
 *  callback may be used to prepare the @p data just before the start of the recursions. If
 *  present, it will be called prior the actual decompositions automatically. You may use
 *  the @p data_release callback to properly free the memory of @p data once it is not required
 *  anymore. Hence, it serves as a kind of destructor for @p data which will be called as
 *  soon as the grammar rules of @p fc are re-set to defaults or if the @p fc is destroyed.
 *
 *  @see  vrna_gr_add_aux_c(), vrna_gr_add_aux_m(), vrna_gr_add_aux_m1(), vrna_gr_add_aux(),
 *        vrna_gr_add_aux_exp_f(), #vrna_gr_inside_f, #vrna_gr_outside_f, #vrna_fold_compound_t,
 *        #vrna_auxdata_prepare_f, #vrna_auxdata_free_f
 *
 *  @param  fc            The fold compound that is to be extended by auxiliary grammar rules
 *  @param  cb            The auxiliary grammar callback for the inside step
 *  @param  cb_bt         The auxiliary grammar callback for the outside (backtracking) step
 *  @param  data          A pointer to some data that will be passed through to the inside and outside callbacks
 *  @param  data_prepare  A callback to prepare @p data
 *  @param  data_release  A callback to free-up memory occupied by @p data
 *  @return               The current number of auxiliary grammar rules for the MFE F-decomposition stage, or 0 on error.
 */
unsigned int
vrna_gr_add_aux_f(vrna_fold_compound_t    *fc,
                  vrna_gr_inside_f        cb,
                  vrna_gr_outside_f       cb_bt,
                  void                    *data,
                  vrna_auxdata_prepare_f  data_prepare,
                  vrna_auxdata_free_f     data_release);


/**
 *  @brief  Add an auxiliary grammar rule for the C-decomposition (MFE version)
 *
 *  This function binds callback functions for auxiliary grammar rules (inside and outside)
 *  in the C-decomposition, i.e. the base pair decomposition stage.
 *
 *  While the inside rule (@p cb) computes a minimum free energy contribution for any
 *  subsequence the outside rule (@p cb_bt) is used for backtracking the corresponding
 *  structure.
 *
 *  Both callbacks will be provided with the @p data pointer that can be used to
 *  store whatever data is needed in the callback evaluations. The @p data_prepare
 *  callback may be used to prepare the @p data just before the start of the recursions. If
 *  present, it will be called prior the actual decompositions automatically. You may use
 *  the @p data_release callback to properly free the memory of @p data once it is not required
 *  anymore. Hence, it serves as a kind of destructor for @p data which will be called as
 *  soon as the grammar rules of @p fc are re-set to defaults or if the @p fc is destroyed.
 *
 *  @see  vrna_gr_add_aux_f(), vrna_gr_add_aux_m(), vrna_gr_add_aux_m1(), vrna_gr_add_aux(),
 *        vrna_gr_add_aux_exp_c(), #vrna_gr_inside_f, #vrna_gr_outside_f, #vrna_fold_compound_t,
 *        #vrna_auxdata_prepare_f, #vrna_auxdata_free_f
 *
 *  @param  fc            The fold compound that is to be extended by auxiliary grammar rules
 *  @param  cb            The auxiliary grammar callback for the inside step
 *  @param  cb_bt         The auxiliary grammar callback for the outside (backtracking) step
 *  @param  data          A pointer to some data that will be passed through to the inside and outside callbacks
 *  @param  data_prepare  A callback to prepare @p data
 *  @param  data_release  A callback to free-up memory occupied by @p data
 *  @return               The current number of auxiliary grammar rules for the MFE C-decomposition stage, or 0 on error.
 */
unsigned int
vrna_gr_add_aux_c(vrna_fold_compound_t    *fc,
                  vrna_gr_inside_f        cb,
                  vrna_gr_outside_f       cb_bt,
                  void                    *data,
                  vrna_auxdata_prepare_f  data_prepare,
                  vrna_auxdata_free_f     data_release);


/**
 *  @brief  Add an auxiliary grammar rule for the M-decomposition (MFE version)
 *
 *  This function binds callback functions for auxiliary grammar rules (inside and outside)
 *  in the M-decomposition, i.e. the multibranch loop decomposition stage.
 *
 *  While the inside rule (@p cb) computes a minimum free energy contribution for any
 *  subsequence the outside rule (@p cb_bt) is used for backtracking the corresponding
 *  structure.
 *
 *  Both callbacks will be provided with the @p data pointer that can be used to
 *  store whatever data is needed in the callback evaluations. The @p data_prepare
 *  callback may be used to prepare the @p data just before the start of the recursions. If
 *  present, it will be called prior the actual decompositions automatically. You may use
 *  the @p data_release callback to properly free the memory of @p data once it is not required
 *  anymore. Hence, it serves as a kind of destructor for @p data which will be called as
 *  soon as the grammar rules of @p fc are re-set to defaults or if the @p fc is destroyed.
 *
 *  @see  vrna_gr_add_aux_f(), vrna_gr_add_aux_c(), vrna_gr_add_aux_m1(), vrna_gr_add_aux(),
 *        vrna_gr_add_aux_exp_m(), #vrna_gr_inside_f, #vrna_gr_outside_f, #vrna_fold_compound_t,
 *        #vrna_auxdata_prepare_f, #vrna_auxdata_free_f
 *
 *  @param  fc            The fold compound that is to be extended by auxiliary grammar rules
 *  @param  cb            The auxiliary grammar callback for the inside step
 *  @param  cb_bt         The auxiliary grammar callback for the outside (backtracking) step
 *  @param  data          A pointer to some data that will be passed through to the inside and outside callbacks
 *  @param  data_prepare  A callback to prepare @p data
 *  @param  data_release  A callback to free-up memory occupied by @p data
 *  @return               The current number of auxiliary grammar rules for the MFE M-decomposition stage, or 0 on error.
 */
unsigned int
vrna_gr_add_aux_m(vrna_fold_compound_t    *fc,
                  vrna_gr_inside_f        cb,
                  vrna_gr_outside_f       cb_bt,
                  void                    *data,
                  vrna_auxdata_prepare_f  data_prepare,
                  vrna_auxdata_free_f     data_release);


/**
 *  @brief  Add an auxiliary grammar rule for the M1-decomposition (MFE version)
 *
 *  This function binds callback functions for auxiliary grammar rules (inside and outside)
 *  in the M1-decomposition, i.e. the multibranch loop components with exactly one branch
 *  decomposition stage.
 *
 *  While the inside rule (@p cb) computes a minimum free energy contribution for any
 *  subsequence the outside rule (@p cb_bt) is used for backtracking the corresponding
 *  structure.
 *
 *  Both callbacks will be provided with the @p data pointer that can be used to
 *  store whatever data is needed in the callback evaluations. The @p data_prepare
 *  callback may be used to prepare the @p data just before the start of the recursions. If
 *  present, it will be called prior the actual decompositions automatically. You may use
 *  the @p data_release callback to properly free the memory of @p data once it is not required
 *  anymore. Hence, it serves as a kind of destructor for @p data which will be called as
 *  soon as the grammar rules of @p fc are re-set to defaults or if the @p fc is destroyed.
 *
 *  @see  vrna_gr_add_aux_f(), vrna_gr_add_aux_c(), vrna_gr_add_aux_m(), vrna_gr_add_aux(),
 *        vrna_gr_add_aux_exp_m1(), #vrna_gr_inside_f, #vrna_gr_outside_f, #vrna_fold_compound_t,
 *        #vrna_auxdata_prepare_f, #vrna_auxdata_free_f
 *
 *  @param  fc            The fold compound that is to be extended by auxiliary grammar rules
 *  @param  cb            The auxiliary grammar callback for the inside step
 *  @param  cb_bt         The auxiliary grammar callback for the outside (backtracking) step
 *  @param  data          A pointer to some data that will be passed through to the inside and outside callbacks
 *  @param  data_prepare  A callback to prepare @p data
 *  @param  data_release  A callback to free-up memory occupied by @p data
 *  @return               The current number of auxiliary grammar rules for the MFE M1-decomposition stage, or 0 on error.
 */
unsigned int
vrna_gr_add_aux_m1(vrna_fold_compound_t   *fc,
                   vrna_gr_inside_f       cb,
                   vrna_gr_outside_f      cb_bt,
                   void                   *data,
                   vrna_auxdata_prepare_f data_prepare,
                   vrna_auxdata_free_f    data_release);


/**
 *  @brief  Add an auxiliary grammar rule for the M2-decomposition (MFE version)
 *
 *  This function binds callback functions for auxiliary grammar rules (inside and outside)
 *  in the M2-decomposition, i.e. the multibranch loop components with at least two branches.
 *
 *  While the inside rule (@p cb) computes a minimum free energy contribution for any
 *  subsequence the outside rule (@p cb_bt) is used for backtracking the corresponding
 *  structure.
 *
 *  Both callbacks will be provided with the @p data pointer that can be used to
 *  store whatever data is needed in the callback evaluations. The @p data_prepare
 *  callback may be used to prepare the @p data just before the start of the recursions. If
 *  present, it will be called prior the actual decompositions automatically. You may use
 *  the @p data_release callback to properly free the memory of @p data once it is not required
 *  anymore. Hence, it serves as a kind of destructor for @p data which will be called as
 *  soon as the grammar rules of @p fc are re-set to defaults or if the @p fc is destroyed.
 *
 *  @see  vrna_gr_add_aux_f(), vrna_gr_add_aux_c(), vrna_gr_add_aux_m(), vrna_gr_add_aux(),
 *        vrna_gr_add_aux_exp_m1(), #vrna_gr_inside_f, #vrna_gr_outside_f, #vrna_fold_compound_t,
 *        #vrna_auxdata_prepare_f, #vrna_auxdata_free_f
 *
 *  @param  fc            The fold compound that is to be extended by auxiliary grammar rules
 *  @param  cb            The auxiliary grammar callback for the inside step
 *  @param  cb_bt         The auxiliary grammar callback for the outside (backtracking) step
 *  @param  data          A pointer to some data that will be passed through to the inside and outside callbacks
 *  @param  data_prepare  A callback to prepare @p data
 *  @param  data_release  A callback to free-up memory occupied by @p data
 *  @return               The current number of auxiliary grammar rules for the MFE M2-decomposition stage, or 0 on error.
 */
unsigned int
vrna_gr_add_aux_m2(vrna_fold_compound_t   *fc,
                   vrna_gr_inside_f       cb,
                   vrna_gr_outside_f      cb_bt,
                   void                   *data,
                   vrna_auxdata_prepare_f data_prepare,
                   vrna_auxdata_free_f    data_release);


/**
 *  @brief  Add an auxiliary grammar rule (MFE version)
 *
 *  This function binds callback functions for auxiliary grammar rules (inside and outside)
 *  as additional, independent decomposition steps.
 *
 *  While the inside rule (@p cb) computes a minimum free energy contribution for any
 *  subsequence the outside rule (@p cb_bt) is used for backtracking the corresponding
 *  structure.
 *
 *  Both callbacks will be provided with the @p data pointer that can be used to
 *  store whatever data is needed in the callback evaluations. The @p data_prepare
 *  callback may be used to prepare the @p data just before the start of the recursions. If
 *  present, it will be called prior the actual decompositions automatically. You may use
 *  the @p data_release callback to properly free the memory of @p data once it is not required
 *  anymore. Hence, it serves as a kind of destructor for @p data which will be called as
 *  soon as the grammar rules of @p fc are re-set to defaults or if the @p fc is destroyed.
 *
 *  @see  vrna_gr_add_aux_f(), vrna_gr_add_aux_c(), vrna_gr_add_aux_m(), vrna_gr_add_aux_m1(),
 *        vrna_gr_add_aux_exp(), #vrna_gr_inside_f, #vrna_gr_outside_f, #vrna_fold_compound_t,
 *        #vrna_auxdata_prepare_f, #vrna_auxdata_free_f
 *
 *  @param  fc            The fold compound that is to be extended by auxiliary grammar rules
 *  @param  cb            The auxiliary grammar callback for the inside step
 *  @param  cb_bt         The auxiliary grammar callback for the outside (backtracking) step
 *  @param  data          A pointer to some data that will be passed through to the inside and outside callbacks
 *  @param  data_prepare  A callback to prepare @p data
 *  @param  data_release  A callback to free-up memory occupied by @p data
 *  @return               The current number of auxiliary grammar rules for MFE predictions, or 0 on error.
 */
unsigned int
vrna_gr_add_aux(vrna_fold_compound_t    *fc,
                vrna_gr_inside_f        cb,
                vrna_gr_outside_f       cb_bt,
                void                    *data,
                vrna_auxdata_prepare_f  data_prepare,
                vrna_auxdata_free_f     data_release);


/**
 * @}
 */

#endif
