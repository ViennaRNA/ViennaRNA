#ifndef __VIENNA_RNA_PACKAGE_PARAMS_H__
#define __VIENNA_RNA_PACKAGE_PARAMS_H__

#include "energy_const.h"
#include "data_structures.h"

/**
 *  \file params.h
 *  \brief Several functions to obtain (pre)scaled energy parameter data containers
 */

extern paramT *scale_parameters(void);

extern paramT *copy_parameters(void);

extern paramT *set_parameters(paramT *dest);


pf_paramT *scale_pf_parameters(void);

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
 *  \see get_scaled_pf_parameters();
 *
 *  \param  dangle_model  The dangle model to be used (possible values: 0 or 2)
 *  \param  temperature   The temperature in degC used for (re-)scaling the energy contributions
 *  \param  alpha         A scaling value that is used as a multiplication factor for the absolute
 *                        temperature of the system
 *  \returns              A set of precomputed Boltzmann factors
 */
pf_paramT *get_boltzmann_factors( int dangle_model,
                                  double temperature,
                                  double alpha,
                                  double pf_scale);

pf_paramT *get_boltzmann_factor_copy(pf_paramT *parameters);

pf_paramT *get_scaled_alipf_parameters(unsigned int n_seq);

extern pf_paramT *copy_pf_param(void);

extern pf_paramT *set_pf_param(paramT *dest);

#endif
