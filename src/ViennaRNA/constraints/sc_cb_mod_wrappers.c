#include <stdlib.h>

#include "ViennaRNA/fold_compound.h"
#include "ViennaRNA/utils/basic.h"
#include "ViennaRNA/static/energy_parameter_sets.h"
#include "ViennaRNA/constraints/soft_special.h"

/*
 #################################
 # BEGIN OF FUNCTION DEFINITIONS #
 #################################
 */
PUBLIC int
vrna_sc_mod_m6A(vrna_fold_compound_t  *fc,
                const unsigned int    *modification_sites,
                unsigned int          options)
{
  return vrna_sc_mod_json(fc,
                          (const char *)parameter_set_rna_mod_m6A_parameters,
                          modification_sites,
                          options);
}


PUBLIC int
vrna_sc_mod_pseudouridine(vrna_fold_compound_t  *fc,
                          const unsigned int    *modification_sites,
                          unsigned int            options)
{
  return vrna_sc_mod_json(fc,
                          (const char *)parameter_set_rna_mod_pseudouridine_parameters,
                          modification_sites,
                          options);
}


PUBLIC int
vrna_sc_mod_inosine(vrna_fold_compound_t  *fc,
                    const unsigned int    *modification_sites,
                    unsigned int          options)
{
  return vrna_sc_mod_json(fc,
                          (const char *)parameter_set_rna_mod_inosine_parameters,
                          modification_sites,
                          options);
}


PUBLIC int
vrna_sc_mod_7DA(vrna_fold_compound_t  *fc,
                const unsigned int    *modification_sites,
                unsigned int          options)
{
  return vrna_sc_mod_json(fc,
                          (const char *)parameter_set_rna_mod_7DA_parameters,
                          modification_sites,
                          options);
}


PUBLIC int
vrna_sc_mod_purine(vrna_fold_compound_t *fc,
                   const unsigned int   *modification_sites,
                   unsigned int         options)
{
  return vrna_sc_mod_json(fc,
                          (const char *)parameter_set_rna_mod_purine_parameters,
                          modification_sites,
                          options);
}


PUBLIC int
vrna_sc_mod_dihydrouridine(vrna_fold_compound_t *fc,
                           const unsigned int   *modification_sites,
                           unsigned int         options)
{
  return vrna_sc_mod_json(fc,
                          (const char *)parameter_set_rna_mod_dihydrouridine_parameters,
                          modification_sites,
                          options);
}
