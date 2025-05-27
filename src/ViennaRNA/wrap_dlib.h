#ifndef VIENNARNA_DLIB_WRAPPER_H
#define VIENNARNA_DLIB_WRAPPER_H

#ifdef __cplusplus
extern "C" {
#endif


/**
 *  @brief  Compute equilibrium concentration of interacting RNA complexes
 *
 *  This function allows one to compute, given a set of initial RNA monomer concentrations,
 *  equilibrium concentrations for a set of interacting RNA complexes.
 *
 *  @ingroup  pf_cofold
 *
 *  @param  eq_constants            The equilibrium constants for the individual complexes
 *  @param  concentration_strands   The initial concentrations of the monomer species (also used for output of the equilibrium concentrations for the monomer species)
 *  @param  A                       An @p NxM matrix indicating the composition of the individual complexes (M), i.e. how many of the monomers (N) each complex consists of
 *  @param  num_strands             The number of monomer species
 *  @param  num_complexes           The number of interacting complexes
 *  @return                         The equilibrium concentrations for each complex
 */
double *
vrna_equilibrium_conc(const double        *eq_constants,
                      double              *concentration_strands,
                      const unsigned int  **A,
                      size_t              num_strands,
                      size_t              num_complexes);


#ifdef __cplusplus
}
#endif

#endif
