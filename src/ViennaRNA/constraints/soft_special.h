#ifndef VIENNA_RNA_PACKAGE_CONSTRAINTS_SOFT_SPECIAL_H
#define VIENNA_RNA_PACKAGE_CONSTRAINTS_SOFT_SPECIAL_H

/**
 *  @file     constraints/soft_special.h
 *  @ingroup  soft_constraints
 *  @brief    Specialized implementations that utilize the soft constraint callback mechanism
 */


void
vrna_sc_m6A(vrna_fold_compound_t  *fc,
            const unsigned int    *modification_sites);

void
vrna_sc_psi(vrna_fold_compound_t  *fc,
            const unsigned int    *modification_sites);

void
vrna_sc_dihydrouridine(vrna_fold_compound_t *fc,
                       const unsigned int   *modification_sites);

#endif
