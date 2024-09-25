#ifndef VIENNA_RNA_PACKAGE_LOOPS_ALL_H
#define VIENNA_RNA_PACKAGE_LOOPS_ALL_H

/**
 *  @file     ViennaRNA/loops/all.h
 *  @ingroup  eval, eval_loops
 *  @brief    Energy evaluation for MFE and partition function calculations
 *
 *  <P>
 *  This file contains functions for the calculation of the free energy @f$\Delta G@f$
 *  of a hairpin- [ vrna_E_hairpin() ] or internal-loop [ vrna_E_internal()] .<BR>
 *  The unit of the free energy returned is @f$10^{-2} * \mathrm{kcal}/\mathrm{mol}@f$
 *  </P>
 *  <P>
 *  In case of computing the partition function, this file also supplies functions
 *  which return the Boltzmann weights @f$e^{-\Delta G/kT} @f$ for a hairpin- [ vrna_exp_E_hairpin() ]
 *  or internal-loop [ vrna_exp_E_internal() ].
 *  </P>
 */

/**
 *  @addtogroup   eval_loops
 *  @{
 */


/* below we include the loop type specific energy evaluation functions */

#include <ViennaRNA/loops/external.h>

#include <ViennaRNA/loops/hairpin.h>

#include <ViennaRNA/loops/internal.h>

#include <ViennaRNA/loops/multibranch.h>

/**
 * @}
 */

#endif
