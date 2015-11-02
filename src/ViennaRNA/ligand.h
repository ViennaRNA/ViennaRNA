#ifndef VIENNA_RNA_PACKAGE_LIGAND_H
#define VIENNA_RNA_PACKAGE_LIGAND_H

#include <ViennaRNA/data_structures.h>

/**
 *  @brief  Add soft constraints for hairpin or interior loop binding motif
 *
 *  @ingroup soft_constraints
 *
 *  @param  vc        The #vrna_fold_compound_t the motif is applied to
 *  @param  seq       The sequence motif (may be interspaced by '&' character
 *  @param  structure The structure motif (may be interspaced by '&' character
 *  @double energy    The free energy of the motif (e.g. binding free energy)
 *  @param  options   Options indicating whether to use the motif in MFE prediction, and/or PF predictions
 *  @return           non-zero value if application of the motif using soft constraints was successful
 *  
 */
int
vrna_sc_add_hi_motif( vrna_fold_compound_t *vc,
                      const char *seq,
                      const char *structure,
                      double energy,
                      unsigned int options);

#endif
