
#include "ViennaRNA/utils/basic.h"
#include "ViennaRNA/model.h"

#include "gengetopt_helpers.h"

PUBLIC void
set_geometry(vrna_md_t *md) {
  if (md) {
    md->helical_rise = vrna_md_defaults_helical_rise_get();
    md->backbone_length = vrna_md_defaults_backbone_length_get();
  }
}

PUBLIC void
set_salt_DNA(vrna_md_t *md) {
  if (md) { /* in case of DNA, also reset typical variables for salt correction */
    md->helical_rise = VRNA_MODEL_HELICAL_RISE_DNA;
    md->backbone_length = VRNA_MODEL_BACKBONE_LENGTH_DNA;
    md->saltDPXInitFact = VRNA_MODEL_SALT_DPXINIT_FACT_DNA;
  }
}
