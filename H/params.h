#ifndef __VIENNA_RNA_PACKAGE_PARAMS_H__
#define __VIENNA_RNA_PACKAGE_PARAMS_H__

#include "energy_const.h"
#include "data_structures.h"

#ifdef __GNUC__
#define DEPRECATED(func) func __attribute__ ((deprecated))
#else
#define DEPRECATED(func) func
#endif

/**
 *  \addtogroup energy_parameters
 *  \brief All relevant functions to retrieve and copy precalculated energy parameter sets as well as
 *  reading/writing the energy parameter set from/to file(s).
 *
 *  This module covers all relevant functions for precalculation of the energy parameters
 *  necessary for the folding routines provided by RNAlib. Furthermore, the energy parameter set
 *  in the RNAlib can be easily exchanged by a user-defined one. It is also possible to write the
 *  current energy parameter set into a text file.
 *  @{
 *
 *  \file params.h
 */

/**
 * \brief Get precomputed energy contributions for all the known loop types
 *
 *  \note OpenMP: This function relies on several global model settings variables and thus is
 *        not to be considered threadsafe. See get_scaled_parameters() for a completely threadsafe
 *        implementation.
 *
 * \return     A set of precomputed energy contributions
 */
paramT *scale_parameters(void);

/**
 * \brief Get precomputed energy contributions for all the known loop types
 *
 *  Call this function to retrieve precomputed energy contributions, i.e. scaled
 *  according to the temperature passed. Furthermore, this function assumes a
 *  data structure that contains the model details as well, such that subsequent
 *  folding recursions are able to retrieve the correct model settings
 *
 *  \see #model_detailsT, set_model_details()
 *
 *  \param temperature  The temperature in degrees Celcius
 *  \param md           The model details
 *  \return             precomputed energy contributions and model settings
 */
paramT *get_scaled_parameters(double temperature,
                              model_detailsT md);

paramT *get_parameter_copy(paramT *par);

/**
 *  get a datastructure of type \ref pf_paramT which contains
 *  the Boltzmann weights of several energy parameters scaled
 *  according to the current temperature
 *  \return The datastructure containing Boltzmann weights for use in partition function calculations
 */
pf_paramT *get_scaled_pf_parameters(void);

/**
 *  \brief Get precomputed Boltzmann factors of the loop type
 *  dependent energy contributions with independent thermodynamic
 *  temperature
 *
 *  This function returns a data structure that contains
 *  all necessary precalculated Boltzmann factors for each
 *  loop type contribution.<br>
 *  In contrast to get_scaled_pf_parameters(), this function
 *  enables setting of independent temperatures for both, the
 *  individual energy contributions as well as the thermodynamic
 *  temperature used in
 *  \f$ exp(-\Delta G / kT) \f$
 *
 *  \see get_scaled_pf_parameters(), get_boltzmann_factor_copy()
 *
 *  \param  temperature   The temperature in degrees Celcius used for (re-)scaling the energy contributions
 *  \param  betaScale     A scaling value that is used as a multiplication factor for the absolute
 *                        temperature of the system
 *  \param  md            The model details to be used
 *  \param  pf_scale      The scaling factor for the Boltzmann factors
 *  \return               A set of precomputed Boltzmann factors
 */
pf_paramT *get_boltzmann_factors( double temperature,
                                  double betaScale,
                                  model_detailsT md,
                                  double pf_scale);

/**
 *  \brief Get a copy of already precomputed Boltzmann factors
 *
 *  \see get_boltzmann_factors(), get_scaled_pf_parameters()
 *
 *  \param  parameters  The input data structure that shall be copied
 *  \return             A copy of the provided Boltzmann factor dataset
 */
pf_paramT *get_boltzmann_factor_copy(pf_paramT *parameters);

/**
 *  \brief Get precomputed Boltzmann factors of the loop type
 *  dependent energy contributions (alifold variant)
 *
 */
pf_paramT *get_scaled_alipf_parameters(unsigned int n_seq);

/**
 *  \brief Get precomputed Boltzmann factors of the loop type
 *  dependent energy contributions (alifold variant) with
 *  independent thermodynamic temperature
 *
 */
PUBLIC pf_paramT *get_boltzmann_factors_ali(unsigned int n_seq,
                                            double temperature,
                                            double betaScale,
                                            model_detailsT md,
                                            double pf_scale);

/**
 *  @}
 */

DEPRECATED(paramT     *copy_parameters(void));
DEPRECATED(paramT     *set_parameters(paramT *dest));
DEPRECATED(pf_paramT  *scale_pf_parameters(void));
DEPRECATED(pf_paramT  *copy_pf_param(void));
DEPRECATED(pf_paramT  *set_pf_param(paramT *dest));



#endif
