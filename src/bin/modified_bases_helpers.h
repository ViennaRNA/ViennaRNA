#ifndef VRNA_MODIFIED_BASES_HELPERS
#define VRNA_MODIFIED_BASES_HELPERS

#include  "ViennaRNA/fold_compound.h"
#include  "ViennaRNA/constraints/soft_special.h"


#define   SPECIAL_BASES_DIHYDROURIDINE    1UL


size_t **
mod_positions_seq_prepare(char                *sequence,
                          unsigned char       mod_dihydrouridine,
                          vrna_sc_mod_param_t *params,
                          int                 verbose,
                          size_t              *param_set_num);


void
mod_bases_apply(vrna_fold_compound_t  *fc,
                size_t                param_set_num,
                size_t                **mod_positions,
                unsigned char         mod_dihydrouridine,
                vrna_sc_mod_param_t   *params);


vrna_sc_mod_param_t *
mod_params_collect_from_string(const char           *string,
                               size_t               *num_params,
                               vrna_sc_mod_param_t  *mod_params,
                               vrna_md_t            *md,
                               unsigned int         *special_bases);


vrna_sc_mod_param_t *
mod_params_collect_from_files(const char          **filenames,
                              unsigned int        file_num,
                              size_t              *num_params,
                              vrna_sc_mod_param_t *mod_params,
                              vrna_md_t           *md);


#endif
