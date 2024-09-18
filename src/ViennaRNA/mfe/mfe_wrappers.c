/*
 * Various wrappers to create simplified interfaces for Minimum Free Energy
 * prediction
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "ViennaRNA/fold_compound.h"
#include "ViennaRNA/model.h"
#include "ViennaRNA/utils/basic.h"
#include "ViennaRNA/backtrack/global.h"
#include "ViennaRNA/mfe/global.h"


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


PUBLIC float
vrna_mfe_dimer(vrna_fold_compound_t *vc,
               char                 *structure)
{
  char          *s2, *ss1, *ss2;
  unsigned int  l1, l2;
  float         mfe, mfe1, mfe2;

  mfe = vrna_mfe(vc, structure);

  /*
   *  for backward compatibility reasons, we alson need
   *  to see whether the unconnected structure is better
   */
  if (vc->strands > 1) {
    l1  = vc->nucleotides[0].length;
    l2  = vc->nucleotides[1].length;
    s2  = vc->nucleotides[1].string;
    ss1 = (char *)vrna_alloc(sizeof(char) * (l1 + 1));
    ss2 = (char *)vrna_alloc(sizeof(char) * (l2 + 1));

    mfe1 = vrna_backtrack5(vc, l1, ss1);

    vrna_fold_compound_t *fc2 = vrna_fold_compound(s2,
                                                   &(vc->params->model_details),
                                                   VRNA_OPTION_DEFAULT);

    mfe2 = vrna_mfe(fc2, ss2);

    if (mfe1 + mfe2 < mfe) {
      mfe = mfe1 + mfe2;
      memcpy(structure, ss1, sizeof(char) * l1);
      memcpy(structure + l1, ss2, sizeof(char) * l2);
      structure[l1 + l2] = '\0';
    }

    vrna_fold_compound_free(fc2);
    free(ss1);
    free(ss2);
  }

  return mfe;
}
