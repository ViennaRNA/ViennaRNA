#ifndef VIENNA_RNA_PACKAGE_CONSTRAINTS_SOFT_SPECIAL_H
#define VIENNA_RNA_PACKAGE_CONSTRAINTS_SOFT_SPECIAL_H

/**
 *  @file     constraints/soft_special.h
 *  @ingroup  soft_constraints
 *  @brief    Specialized implementations that utilize the soft constraint callback mechanism
 */

typedef struct vrna_sc_cb_mod_param_s *vrna_sc_cb_mod_parameters_t;

void
vrna_sc_m6A(vrna_fold_compound_t  *fc,
            const unsigned int    *modification_sites);

void
vrna_sc_psi(vrna_fold_compound_t  *fc,
            const unsigned int    *modification_sites);

void
vrna_sc_dihydrouridine(vrna_fold_compound_t *fc,
                       const unsigned int   *modification_sites);

void
vrna_sc_inosine(vrna_fold_compound_t  *fc,
                const unsigned int    *modification_sites);

vrna_sc_cb_mod_parameters_t
vrna_sc_cb_mod_read_from_json_file(const char *filename,
                                   vrna_md_t *md);


vrna_sc_cb_mod_parameters_t
vrna_sc_cb_mod_read_from_json(const char *json,
                              vrna_md_t *md);

void
vrna_sc_mod_json(vrna_fold_compound_t *fc,
                 const char *json_file,
                 unsigned int *modification_sites);

#endif
