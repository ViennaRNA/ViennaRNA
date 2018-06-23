/*
 * Various wrappers to create simplified interfaces for Minimum Free Energy
 * prediction
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdio.h>
#include <stdlib.h>

#include "ViennaRNA/fold_compound.h"
#include "ViennaRNA/model.h"
#include "ViennaRNA/utils/basic.h"
#include "ViennaRNA/mfe.h"


/* wrappers for single sequences */
PUBLIC float
vrna_fold(const char  *string,
          char        *structure)
{
  float                 mfe;
  vrna_fold_compound_t  *vc;
  vrna_md_t             md;

  vrna_md_set_default(&md);
  vc  = vrna_fold_compound(string, &md, 0);
  mfe = vrna_mfe(vc, structure);

  vrna_fold_compound_free(vc);

  return mfe;
}


PUBLIC float
vrna_circfold(const char  *string,
              char        *structure)
{
  float                 mfe;
  vrna_fold_compound_t  *vc;
  vrna_md_t             md;

  vrna_md_set_default(&md);
  md.circ = 1;
  vc      = vrna_fold_compound(string, &md, 0);
  mfe     = vrna_mfe(vc, structure);

  vrna_fold_compound_free(vc);

  return mfe;
}


/* wrappers for multiple sequence alignments */

PUBLIC float
vrna_alifold(const char **strings,
             char       *structure)
{
  float                 mfe;
  vrna_fold_compound_t  *vc;
  vrna_md_t             md;

  vrna_md_set_default(&md);

  vc  = vrna_fold_compound_comparative(strings, &md, VRNA_OPTION_DEFAULT);
  mfe = vrna_mfe(vc, structure);

  vrna_fold_compound_free(vc);

  return mfe;
}


PUBLIC float
vrna_circalifold(const char **sequences,
                 char       *structure)
{
  float                 mfe;
  vrna_fold_compound_t  *vc;
  vrna_md_t             md;

  vrna_md_set_default(&md);
  md.circ = 1;

  vc  = vrna_fold_compound_comparative(sequences, &md, VRNA_OPTION_DEFAULT);
  mfe = vrna_mfe(vc, structure);

  vrna_fold_compound_free(vc);

  return mfe;
}


/* wrappers for RNA-RNA cofolding interaction */
PUBLIC float
vrna_cofold(const char  *seq,
            char        *structure)
{
  float                 mfe;
  vrna_fold_compound_t  *vc;
  vrna_md_t             md;

  vrna_md_set_default(&md);
  md.min_loop_size = 0;  /* set min loop length to 0 */

  /* get compound structure */
  vc = vrna_fold_compound(seq, &md, 0);

  mfe = vrna_mfe_dimer(vc, structure);

  vrna_fold_compound_free(vc);

  return mfe;
}
