#ifndef VIENNA_RNA_PACKAGE_LIGAND_H
#define VIENNA_RNA_PACKAGE_LIGAND_H

/**
 *  @file ligand.h
 *
 *  @{
 *
 *  @ingroup soft_constraints
 */
#include <ViennaRNA/data_structures.h>

/**
 *  @brief  Add soft constraints for hairpin or interior loop binding motif
 *
 *  @ingroup generalized_sc
 *
 *  Here is an example that adds a theophylline binding motif. Free energy
 *  contribution is derived from @f$k_d = 0.32 \mu mol / l @f$, taken from
 *  Jenison et al. 1994
 *  @code{.c}
vrna_sc_add_hi_motif( vc,
                      "GAUACCAG&CCCUUGGCAGC",
                      "(...((((&)...)))...)",
                      -9.22,
                      VRNA_OPTION_MFE); @endcode
 *
 *  @param  vc        The #vrna_fold_compound_t the motif is applied to
 *  @param  seq       The sequence motif (may be interspaced by '&' character
 *  @param  structure The structure motif (may be interspaced by '&' character
 *  @param  energy    The free energy of the motif (e.g. binding free energy)
 *  @param  options   Options indicating whether to use the motif in MFE prediction, and/or PF predictions
 *  @return           non-zero value if application of the motif using soft constraints was successful
 *  
 */
int
vrna_sc_add_hi_motif( vrna_fold_compound_t *vc,
                      const char *seq,
                      const char *structure,
                      FLT_OR_DBL energy,
                      unsigned int options);


/**
 * @}
 */
#endif
