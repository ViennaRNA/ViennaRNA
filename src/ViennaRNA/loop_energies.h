#ifndef VIENNA_RNA_PACKAGE_LOOP_ENERGIES_H
#define VIENNA_RNA_PACKAGE_LOOP_ENERGIES_H

/**
 *  @file     loop_energies.h
 *  @ingroup  loops
 *  @brief    Energy evaluation for MFE and partition function calculations
 */

/**
 *  @{
 *  @ingroup   loops
 * 
 *  <P>
 *  This file contains functions for the calculation of the free energy @f$\Delta G@f$
 *  of a hairpin- [ E_Hairpin() ] or interior-loop [ E_IntLoop()] .<BR>
 *  The unit of the free energy returned is @f$10^{-2} * \mathrm{kcal}/\mathrm{mol}@f$
 *  </P>
 *  <P>
 *  In case of computing the partition function, this file also supplies functions
 *  which return the Boltzmann weights @f$e^{-\Delta G/kT} @f$ for a hairpin- [ exp_E_Hairpin() ]
 *  or interior-loop [ exp_E_IntLoop() ].
 *  </P>
 */


/* below we include the loop type specific energy evaluation functions */

#include <ViennaRNA/exterior_loops.h>

#include <ViennaRNA/hairpin_loops.h>

#include <ViennaRNA/interior_loops.h>

#include <ViennaRNA/multibranch_loops.h>

/**
 * @}
 */

#endif
