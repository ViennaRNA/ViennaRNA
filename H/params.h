#ifndef __VIENNA_RNA_PACKAGE_PARAMS_H__
#define __VIENNA_RNA_PACKAGE_PARAMS_H__

#include "energy_const.h"
#include "data_structures.h"

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

pf_paramT *get_scaled_alipf_parameters(unsigned int n_seq);

extern pf_paramT *copy_pf_param(void);
extern pf_paramT *set_pf_param(paramT *dest);

#endif
