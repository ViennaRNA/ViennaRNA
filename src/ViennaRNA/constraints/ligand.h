#ifndef VIENNA_RNA_PACKAGE_LIGAND_H
#define VIENNA_RNA_PACKAGE_LIGAND_H

/**
 *  @file     constraints/ligand.h
 *  @ingroup  ligand_binding
 *  @brief    Functions for incorporation of ligands binding to hairpin and interior loop motifs using the soft constraints framework
 */

/**
 *  @addtogroup constraints_ligand
 *
 *  @brief  Ligand binding to specific hairpin/interior loop like motifs using the @ref soft_constraints feature
 *
 *  Here is an example that adds a theophylline binding motif. Free energy
 *  contribution is derived from @f$k_d = 0.32 \mu mol / l @f$, taken from
 *  Jenison et al. 1994
 *
 *  @image html   theo_aptamer.svg
 *  @image latex  theo_aptamer.eps
 *
 *  @code{.c}
vrna_sc_add_hi_motif(vc,
                   "GAUACCAG&CCCUUGGCAGC",
                   "(...((((&)...)))...)",
                   -9.22, VRNA_OPTION_DEFAULT); @endcode
 */


typedef struct vrna_sc_motif_s vrna_sc_motif_t;

#include <ViennaRNA/datastructures/basic.h>
#include <ViennaRNA/fold_compound.h>

struct vrna_sc_motif_s {
  int i;
  int j;
  int k;
  int l;
  int number;
};


/**
 *  @brief  Add soft constraints for hairpin or interior loop binding motif
 *
 *  @ingroup  constraints_ligand
 *
 *  @param  vc        The #vrna_fold_compound_t the motif is applied to
 *  @param  seq       The sequence motif (may be interspaced by '&' character
 *  @param  structure The structure motif (may be interspaced by '&' character
 *  @param  energy    The free energy of the motif (e.g. binding free energy)
 *  @param  options   Options
 *  @return           non-zero value if application of the motif using soft constraints was successful
 *
 */
int
vrna_sc_add_hi_motif(vrna_fold_compound_t *vc,
                     const char           *seq,
                     const char           *structure,
                     FLT_OR_DBL           energy,
                     unsigned int         options);


vrna_sc_motif_t *
vrna_sc_ligand_detect_motifs(vrna_fold_compound_t *vc,
                             const char           *structure);


vrna_sc_motif_t *
vrna_sc_ligand_get_all_motifs(vrna_fold_compound_t *vc);


#endif
