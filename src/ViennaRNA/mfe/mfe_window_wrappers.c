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
#include "ViennaRNA/mfe/local.h"


/* wrappers for local MFE prediction */
PUBLIC float
vrna_Lfold(const char *string,
           int        window_size,
           FILE       *file)
{
  float                 energy;
  vrna_fold_compound_t  *vc;
  vrna_md_t             md;

  vrna_md_set_default(&md);

  md.window_size  = window_size;
  md.max_bp_span  = window_size;

  vc = vrna_fold_compound(string, &md, VRNA_OPTION_WINDOW);

  energy = vrna_mfe_window(vc, file);

  vrna_fold_compound_free(vc);

  return energy;
}


PUBLIC float
vrna_Lfold_cb(const char        *string,
              int               window_size,
              vrna_mfe_window_f cb,
              void              *data)
{
  float                 energy;
  vrna_fold_compound_t  *vc;
  vrna_md_t             md;

  vrna_md_set_default(&md);

  md.window_size  = window_size;
  md.max_bp_span  = window_size;

  vc = vrna_fold_compound(string, &md, VRNA_OPTION_WINDOW);

  energy = vrna_mfe_window_cb(vc, cb, data);

  vrna_fold_compound_free(vc);

  return energy;
}


#ifdef VRNA_WITH_SVM

PUBLIC float
vrna_Lfoldz(const char  *string,
            int         window_size,
            double      min_z,
            FILE        *file)
{
  float                 energy;
  vrna_fold_compound_t  *vc;
  vrna_md_t             md;

  vrna_md_set_default(&md);

  md.window_size  = window_size;
  md.max_bp_span  = window_size;

  vc = vrna_fold_compound(string, &md, VRNA_OPTION_WINDOW);

  energy = vrna_mfe_window_zscore(vc, min_z, file);

  vrna_fold_compound_free(vc);

  return energy;
}


PUBLIC float
vrna_Lfoldz_cb(const char               *string,
               int                      window_size,
               double                   min_z,
               vrna_mfe_window_zscore_f cb,
               void                     *data)
{
  float                 energy;
  vrna_fold_compound_t  *vc;
  vrna_md_t             md;

  vrna_md_set_default(&md);

  md.window_size  = window_size;
  md.max_bp_span  = window_size;

  vc = vrna_fold_compound(string, &md, VRNA_OPTION_WINDOW);

  energy = vrna_mfe_window_zscore_cb(vc, min_z, cb, data);

  vrna_fold_compound_free(vc);

  return energy;
}


#endif

/* wrappers for local MFE prediction using multiple sequence alignments */
PUBLIC float
vrna_aliLfold(const char  **AS,
              int         window_size,
              FILE        *file)
{
  float                 energy;
  vrna_fold_compound_t  *vc;
  vrna_md_t             md;

  vrna_md_set_default(&md);

  md.window_size  = window_size;
  md.max_bp_span  = window_size;

  vc = vrna_fold_compound_comparative(AS, &md, VRNA_OPTION_MFE | VRNA_OPTION_WINDOW);

  energy = vrna_mfe_window(vc, file);

  vrna_fold_compound_free(vc);

  return energy;
}


PUBLIC float
vrna_aliLfold_cb(const char         **AS,
                 int                window_size,
                 vrna_mfe_window_f  cb,
                 void               *data)
{
  float                 energy;
  vrna_fold_compound_t  *vc;
  vrna_md_t             md;

  vrna_md_set_default(&md);

  md.window_size  = window_size;
  md.max_bp_span  = window_size;

  vc = vrna_fold_compound_comparative(AS, &md, VRNA_OPTION_MFE | VRNA_OPTION_WINDOW);

  energy = vrna_mfe_window_cb(vc, cb, data);

  vrna_fold_compound_free(vc);

  return energy;
}
