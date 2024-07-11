#ifndef VIENNA_RNA_PACKAGE_GRAMMAR_H
#define VIENNA_RNA_PACKAGE_GRAMMAR_H

/**
 *  @file     grammar.h
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

/**
 *  @brief  A pointer to the auxiliary grammar data structure
 */
typedef struct vrna_gr_aux_s *vrna_gr_aux_t;


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
                                int                   i,
                                int                   j,
                                void                  *data);


/**
 *  @brief  Function prototype for auxiliary grammar rules (inside version, partition function)
 *
 *  This function will be called during the inside recursions of partition function predictions
 *  for subsequences from position @p i to @p j and is supposed to return a Boltzmann factor
 *  with energy contribution in cal/mol.
 *
 *  @callback
 *  @parblock
 *  This callback allows for extending the partition function secondary structure decomposition
 *  with additional rules.
 *  @endparblock
 *
 *  @see  vrna_gr_add_aux_exp_f(), vrna_gr_add_aux_exp_c(), vrna_gr_add_aux_exp_m(),
 *        vrna_gr_add_aux_exp_m1(), vrna_gr_add_aux_exp(), vrna_gr_add_aux_f()
 *
 *  @param  fc    The fold compound to work on
 *  @param  i     The 5' delimiter of the sequence segment
 *  @param  j     The 3' delimiter of the sequence segment
 *  @param  data  An arbitrary user-provided data pointer
 *  @return       The partition function omputed by the auxiliary grammar rule with energies in cal/mol
 */
typedef FLT_OR_DBL (*vrna_gr_inside_exp_f)(vrna_fold_compound_t *fc,
                                           int                  i,
                                           int                  j,
                                           void                 *data);


typedef FLT_OR_DBL (*vrna_gr_outside_exp_f)(vrna_fold_compound_t  *fc,
                                            int                   i,
                                            int                   j,
                                            void                  *data);


/**
 *  @brief  Function prototype for auxiliary grammar rules (outside version, MFE)
 *
 *  This function will be called during the backtracking stage of MFE predictions
 *  for subsequences from position @p i to @p j and is supposed to identify the
 *  the structure components that give rise to the corresponding energy contribution
 *  determined in the inside (forward) step.
 *
 *  For that purpose, the function may modify the base pair (@p bp_stack) by adding
 *  all base pairs identified through the additional rule. In addition, when the rule
 *  sub-divides the sequence segment into smaller pieces, these pieces can be put onto
 *  the backtracking (@p bt_stack) stack(s). The stack sizes always have to be increased
 *  accordingly. On successful backtrack, the function returns non-zero, and exactly
 *  zero (0) when the backtracking failed.
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
 *  @param  bp_stack      The base pair stack
 *  @param  bp_stack_size A pointer to the current size of the @p bp_stack size
 *  @param  bt_stack      The backtracking stack (all segments that need to be further investigated)
 *  @param  bt_stack_size A pointer to the current size of the @p bt_stack size
 *  @param  data  An arbitrary user-provided data pointer
 *  @return       The free energy computed by the auxiliary grammar rule in dekacal/mol
 */
typedef int (*vrna_gr_outside_f)(vrna_fold_compound_t *fc,
                                 unsigned int         i,
                                 unsigned int         j,
                                 vrna_bp_stack_t      *bp_stack,
                                 int                  *bp_stack_size,
                                 vrna_sect_t          *bt_stack,
                                 int                  *bt_stack_size,
                                 void                 *data);


/**
 *  @brief  Prepare the auxiliary grammar rule data
 *
 *  @note   This function is mainly for internal use. Users of the auxiliary grammar API usually
 *          do not need to call this function except for debugging purposes.
 *
 *  @param  fc      The fold compound storing the auxiliary grammar rules
 *  @param  options Options flag(s) that denote what needs to be prepared
 *  @return         non-zero on success, 0 otherwise
 */
unsigned int
vrna_gr_prepare(vrna_fold_compound_t  *fc,
                unsigned int          options);


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
 *  @brief  Add an auxiliary grammar rule for the F-decomposition (partition function version)
 *
 *  This function binds callback functions for auxiliary grammar rules (inside and outside)
 *  in the F-decomposition, i.e. the external loop decomposition stage.
 *
 *  While the inside rule (@p cb) computes the partition function for any subsequence, the
 *  outside rule (@p cb_out) is used for (base pairing) probabilities.
 *
 *  Both callbacks will be provided with the @p data pointer that can be used to
 *  store whatever data is needed in the callback evaluations. The @p data_prepare
 *  callback may be used to prepare the @p data just before the start of the recursions. If
 *  present, it will be called prior the actual decompositions automatically. You may use
 *  the @p data_release callback to properly free the memory of @p data once it is not required
 *  anymore. Hence, it serves as a kind of destructor for @p data which will be called as
 *  soon as the grammar rules of @p fc are re-set to defaults or if the @p fc is destroyed.
 *
 *  @bug  Calling the @p cb_out callback is not implemented yet!
 *
 *  @see  vrna_gr_add_aux_exp_c(), vrna_gr_add_aux_exp_m(), vrna_gr_add_aux_exp_m1(), vrna_gr_add_aux_exp(),
 *        vrna_gr_add_aux_f(), #vrna_gr_inside_exp_f, #vrna_gr_outside_exp_f, #vrna_fold_compound_t,
 *        #vrna_auxdata_prepare_f, #vrna_auxdata_free_f
 *
 *  @param  fc            The fold compound that is to be extended by auxiliary grammar rules
 *  @param  cb            The auxiliary grammar callback for the inside step
 *  @param  cb_out        The auxiliary grammar callback for the outside (probability) step
 *  @param  data          A pointer to some data that will be passed through to the inside and outside callbacks
 *  @param  data_prepare  A callback to prepare @p data
 *  @param  data_release  A callback to free-up memory occupied by @p data
 *  @return               The current number of auxiliary grammar rules for the partition function F-decomposition stage, or 0 on error.
 */
unsigned int
vrna_gr_add_aux_exp_f(vrna_fold_compound_t    *fc,
                      vrna_gr_inside_exp_f    cb,
                      vrna_gr_outside_exp_f   cb_out,
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
 *  @brief  Add an auxiliary grammar rule for the C-decomposition (partition function version)
 *
 *  This function binds callback functions for auxiliary grammar rules (inside and outside)
 *  in the C-decomposition, i.e. the base pair decomposition stage.
 *
 *  While the inside rule (@p cb) computes the partition function for any subsequence, the
 *  outside rule (@p cb_out) is used for (base pairing) probabilities.
 *
 *  Both callbacks will be provided with the @p data pointer that can be used to
 *  store whatever data is needed in the callback evaluations. The @p data_prepare
 *  callback may be used to prepare the @p data just before the start of the recursions. If
 *  present, it will be called prior the actual decompositions automatically. You may use
 *  the @p data_release callback to properly free the memory of @p data once it is not required
 *  anymore. Hence, it serves as a kind of destructor for @p data which will be called as
 *  soon as the grammar rules of @p fc are re-set to defaults or if the @p fc is destroyed.
 *
 *  @bug  Calling the @p cb_out callback is not implemented yet!
 *
 *  @see  vrna_gr_add_aux_exp_f(), vrna_gr_add_aux_exp_m(), vrna_gr_add_aux_exp_m1(), vrna_gr_add_aux_exp(),
 *        vrna_gr_add_aux_c(), #vrna_gr_inside_exp_f, #vrna_gr_outside_exp_f, #vrna_fold_compound_t,
 *        #vrna_auxdata_prepare_f, #vrna_auxdata_free_f
 *
 *  @param  fc            The fold compound that is to be extended by auxiliary grammar rules
 *  @param  cb            The auxiliary grammar callback for the inside step
 *  @param  cb_out        The auxiliary grammar callback for the outside (probability) step
 *  @param  data          A pointer to some data that will be passed through to the inside and outside callbacks
 *  @param  data_prepare  A callback to prepare @p data
 *  @param  data_release  A callback to free-up memory occupied by @p data
 *  @return               The current number of auxiliary grammar rules for the partition function C-decomposition stage, or 0 on error.
 */
unsigned int
vrna_gr_add_aux_exp_c(vrna_fold_compound_t    *fc,
                      vrna_gr_inside_exp_f    cb,
                      vrna_gr_outside_exp_f   cb_out,
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
 *  @brief  Add an auxiliary grammar rule for the M-decomposition (partition function version)
 *
 *  This function binds callback functions for auxiliary grammar rules (inside and outside)
 *  in the M-decomposition, i.e. the multibranch loop decomposition stage.
 *
 *  While the inside rule (@p cb) computes the partition function for any subsequence, the
 *  outside rule (@p cb_out) is used for (base pairing) probabilities.
 *
 *  Both callbacks will be provided with the @p data pointer that can be used to
 *  store whatever data is needed in the callback evaluations. The @p data_prepare
 *  callback may be used to prepare the @p data just before the start of the recursions. If
 *  present, it will be called prior the actual decompositions automatically. You may use
 *  the @p data_release callback to properly free the memory of @p data once it is not required
 *  anymore. Hence, it serves as a kind of destructor for @p data which will be called as
 *  soon as the grammar rules of @p fc are re-set to defaults or if the @p fc is destroyed.
 *
 *  @bug  Calling the @p cb_out callback is not implemented yet!
 *
 *  @see  vrna_gr_add_aux_exp_f(), vrna_gr_add_aux_exp_c(), vrna_gr_add_aux_exp_m1(), vrna_gr_add_aux_exp(),
 *        vrna_gr_add_aux_m(), #vrna_gr_inside_exp_f, #vrna_gr_outside_exp_f, #vrna_fold_compound_t,
 *        #vrna_auxdata_prepare_f, #vrna_auxdata_free_f
 *
 *  @param  fc            The fold compound that is to be extended by auxiliary grammar rules
 *  @param  cb            The auxiliary grammar callback for the inside step
 *  @param  cb_out        The auxiliary grammar callback for the outside (probability) step
 *  @param  data          A pointer to some data that will be passed through to the inside and outside callbacks
 *  @param  data_prepare  A callback to prepare @p data
 *  @param  data_release  A callback to free-up memory occupied by @p data
 *  @return               The current number of auxiliary grammar rules for the partition function M-decomposition stage, or 0 on error.
 */
unsigned int
vrna_gr_add_aux_exp_m(vrna_fold_compound_t    *fc,
                      vrna_gr_inside_exp_f    cb,
                      vrna_gr_outside_exp_f   cb_out,
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
 *  @brief  Add an auxiliary grammar rule for the M1-decomposition (partition function version)
 *
 *  This function binds callback functions for auxiliary grammar rules (inside and outside)
 *  in the M1-decomposition, i.e. the multibranch loop components with exactly one branch
 *  decomposition stage.
 *
 *  While the inside rule (@p cb) computes the partition function for any subsequence, the
 *  outside rule (@p cb_out) is used for (base pairing) probabilities.
 *
 *  Both callbacks will be provided with the @p data pointer that can be used to
 *  store whatever data is needed in the callback evaluations. The @p data_prepare
 *  callback may be used to prepare the @p data just before the start of the recursions. If
 *  present, it will be called prior the actual decompositions automatically. You may use
 *  the @p data_release callback to properly free the memory of @p data once it is not required
 *  anymore. Hence, it serves as a kind of destructor for @p data which will be called as
 *  soon as the grammar rules of @p fc are re-set to defaults or if the @p fc is destroyed.
 *
 *  @bug  Calling the @p cb_out callback is not implemented yet!
 *
 *  @see  vrna_gr_add_aux_exp_f(), vrna_gr_add_aux_exp_c(), vrna_gr_add_aux_exp_m(), vrna_gr_add_aux_exp(),
 *        vrna_gr_add_aux_m1(), #vrna_gr_inside_exp_f, #vrna_gr_outside_exp_f, #vrna_fold_compound_t,
 *        #vrna_auxdata_prepare_f, #vrna_auxdata_free_f
 *
 *  @param  fc            The fold compound that is to be extended by auxiliary grammar rules
 *  @param  cb            The auxiliary grammar callback for the inside step
 *  @param  cb_out        The auxiliary grammar callback for the outside (probability) step
 *  @param  data          A pointer to some data that will be passed through to the inside and outside callbacks
 *  @param  data_prepare  A callback to prepare @p data
 *  @param  data_release  A callback to free-up memory occupied by @p data
 *  @return               The current number of auxiliary grammar rules for the partition function M1-decomposition stage, or 0 on error.
 */
unsigned int
vrna_gr_add_aux_exp_m1(vrna_fold_compound_t   *fc,
                       vrna_gr_inside_exp_f   cb,
                       vrna_gr_outside_exp_f  cb_out,
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
 *  @brief  Add an auxiliary grammar rule (partition function version)
 *
 *  This function binds callback functions for auxiliary grammar rules (inside and outside)
 *  as additional, independent decomposition steps.
 *
 *  While the inside rule (@p cb) computes the partition function for any subsequence, the
 *  outside rule (@p cb_out) is used for (base pairing) probabilities.
 *
 *  Both callbacks will be provided with the @p data pointer that can be used to
 *  store whatever data is needed in the callback evaluations. The @p data_prepare
 *  callback may be used to prepare the @p data just before the start of the recursions. If
 *  present, it will be called prior the actual decompositions automatically. You may use
 *  the @p data_release callback to properly free the memory of @p data once it is not required
 *  anymore. Hence, it serves as a kind of destructor for @p data which will be called as
 *  soon as the grammar rules of @p fc are re-set to defaults or if the @p fc is destroyed.
 *
 *  @bug  Calling the @p cb_out callback is not implemented yet!
 *
 *  @see  vrna_gr_add_aux_exp_f(), vrna_gr_add_aux_exp_c(), vrna_gr_add_aux_exp_m(), vrna_gr_add_aux_exp_m1(),
 *        vrna_gr_add_aux(), #vrna_gr_inside_exp_f, #vrna_gr_outside_exp_f, #vrna_fold_compound_t,
 *        #vrna_auxdata_prepare_f, #vrna_auxdata_free_f
 *
 *  @param  fc            The fold compound that is to be extended by auxiliary grammar rules
 *  @param  cb            The auxiliary grammar callback for the inside step
 *  @param  cb_out        The auxiliary grammar callback for the outside (probability) step
 *  @param  data          A pointer to some data that will be passed through to the inside and outside callbacks
 *  @param  data_prepare  A callback to prepare @p data
 *  @param  data_release  A callback to free-up memory occupied by @p data
 *  @return               The current number of auxiliary grammar rules for partition function predictions, or 0 on error.
 */
unsigned int
vrna_gr_add_aux_exp(vrna_fold_compound_t    *fc,
                    vrna_gr_inside_exp_f    cb,
                    vrna_gr_outside_exp_f   cb_out,
                    void                    *data,
                    vrna_auxdata_prepare_f  data_prepare,
                    vrna_auxdata_free_f     data_release);


/**
 *  @brief  Add status event based data preparation callbacks
 *
 *  This function binds additional data structures and corresponding callback
 *  functions for the auxiliary grammar extension API. This might be helpful whenever
 *  certain preparation steps need to be done prior and/or after the actual run of the
 *  prediction algorithms.
 *
 *  @param  fc      The fold compound
 *  @param  cb      The recursion status callback that performs the preparation
 *  @param  data    The data pointer the @p cb callback is working on
 *  @param  prepare_cb  A preparation callback function for parameter @p data
 *  @param  free_cb     A callback to release memory for @p data
 *  @return             The number of status function callbacks bound to the fold compound or 0 on error
 */
unsigned int
vrna_gr_add_status(vrna_fold_compound_t     *fc,
                   vrna_recursion_status_f  cb,
                   void                     *data,
                   vrna_auxdata_prepare_f   prepare_cb,
                   vrna_auxdata_free_f      free_cb);


/**
 *  @brief  Remove all auxiliary grammar rules
 *
 *  This function re-sets the fold compound to the default rules
 *  by removing all auxiliary grammar rules
 *
 *  @param  fc  The fold compound
 */
void
vrna_gr_reset(vrna_fold_compound_t *fc);


/**
 * @}
 */

#endif
