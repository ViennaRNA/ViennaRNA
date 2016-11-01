#ifndef VIENNA_RNA_PACKAGE_LIGAND_H
#define VIENNA_RNA_PACKAGE_LIGAND_H

/**
 *  @file     constraints_ligand.h
 *  @ingroup  ligand_binding
 *  @brief    Functions for incorporation of ligands binding to hairpin and interior loop motifs using the soft constraints framework
 */
typedef struct vrna_sc_motif_s  vrna_sc_motif_t;

#include <ViennaRNA/data_structures.h>

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
vrna_sc_add_hi_motif( vrna_fold_compound_t *vc,
                      const char *seq,
                      const char *structure,
                      FLT_OR_DBL energy,
                      unsigned int options);

vrna_sc_motif_t *
vrna_sc_ligand_detect_motifs( vrna_fold_compound_t *vc,
                              const char *structure);

int
vrna_sc_get_hi_motif( vrna_fold_compound_t *vc,
                      int *i,
                      int *j,
                      int *k,
                      int *l);

#endif
