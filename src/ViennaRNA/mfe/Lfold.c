/*
 *                minimum free energy
 *                RNA secondary structure prediction
 *                with maximum distance base pairs
 *
 *                c Ivo Hofacker, Peter Stadler
 *
 *                ViennaRNA package
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#ifndef VRNA_DISABLE_BACKWARD_COMPATIBILITY

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>

#include "ViennaRNA/utils/basic.h"
#include "ViennaRNA/datastructures/basic.h"
#include "ViennaRNA/mfe/local.h"
#include "ViennaRNA/Lfold.h"


/*
 #################################
 # BEGIN OF FUNCTION DEFINITIONS #
 #################################
 */


/*
 * ###########################################
 * # deprecated functions below              #
 *###########################################
 */
PUBLIC float
Lfold(const char  *string,
      const char  *structure VRNA_UNUSED,
      int         window_size)
{
  float                 energy;
  vrna_fold_compound_t  *vc;
  vrna_md_t             md;

  set_model_details(&md);

  md.window_size  = window_size;
  md.max_bp_span  = window_size;

  vc = vrna_fold_compound(string, &md, VRNA_OPTION_WINDOW);

  energy = vrna_mfe_window(vc, NULL);

  vrna_fold_compound_free(vc);

  return energy;
}


#ifdef VRNA_WITH_SVM

PUBLIC float
Lfoldz(const char *string,
       const char *structure VRNA_UNUSED,
       int        window_size,
       int        zsc,
       double     min_z)
{
  float                 energy;
  vrna_fold_compound_t  *vc;
  vrna_md_t             md;

  set_model_details(&md);

  md.window_size  = window_size;
  md.max_bp_span  = window_size;

  vc = vrna_fold_compound(string, &md, VRNA_OPTION_WINDOW);

  energy = (zsc) ? vrna_mfe_window_zscore(vc, min_z, NULL) : vrna_mfe_window(vc, NULL);

  vrna_fold_compound_free(vc);

  return energy;
}


#endif


PUBLIC float
aliLfold(const char *AS[],
         const char *structure VRNA_UNUSED,
         int        maxdist)
{
  float                 en;
  vrna_fold_compound_t  *fc;
  vrna_md_t             md;

  set_model_details(&md);

  md.max_bp_span = md.window_size = maxdist;

  fc = vrna_fold_compound_comparative((const char **)AS, &md, VRNA_OPTION_MFE | VRNA_OPTION_WINDOW);

  en = vrna_mfe_window(fc, NULL);

  vrna_fold_compound_free(fc);

  return en;
}


PUBLIC float
aliLfold_cb(const char        **AS,
            int               maxdist,
            vrna_mfe_window_f cb,
            void              *data)
{
  float                 en;
  vrna_fold_compound_t  *fc;
  vrna_md_t             md;

  set_model_details(&md);

  md.max_bp_span = md.window_size = maxdist;

  fc = vrna_fold_compound_comparative(AS, &md, VRNA_OPTION_MFE | VRNA_OPTION_WINDOW);

  en = vrna_mfe_window_cb(fc, cb, data);

  vrna_fold_compound_free(fc);

  return en;
}


#endif
