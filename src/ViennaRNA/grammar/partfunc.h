#ifndef VIENNA_RNA_PACKAGE_GRAMMAR_PARTFUNC_H
#define VIENNA_RNA_PACKAGE_GRAMMAR_PARTFUNC_H

/**
 *  @file     ViennaRNA/grammar/partfunc.h
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
                                           unsigned int         i,
                                           unsigned int         j,
                                           void                 *data);


typedef FLT_OR_DBL (*vrna_gr_outside_exp_f)(vrna_fold_compound_t  *fc,
                                            unsigned int          i,
                                            unsigned int          j,
                                            void                  *data);


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
 *  @brief  Add an auxiliary grammar rule for the M2-decomposition (partition function version)
 *
 *  This function binds callback functions for auxiliary grammar rules (inside and outside)
 *  in the M2-decomposition, i.e. the multibranch loop components with at least two branches.
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
 *  @return               The current number of auxiliary grammar rules for the partition function M2-decomposition stage, or 0 on error.
 */
unsigned int
vrna_gr_add_aux_exp_m2(vrna_fold_compound_t   *fc,
                       vrna_gr_inside_exp_f   cb,
                       vrna_gr_outside_exp_f  cb_out,
                       void                   *data,
                       vrna_auxdata_prepare_f data_prepare,
                       vrna_auxdata_free_f    data_release);


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
 * @}
 */

#endif
