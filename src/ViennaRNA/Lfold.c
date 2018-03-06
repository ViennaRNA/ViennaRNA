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

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>

#include "ViennaRNA/utils.h"
#include "ViennaRNA/data_structures.h"
#include "ViennaRNA/Lfold.h"


/*
 #################################
 # BEGIN OF FUNCTION DEFINITIONS #
 #################################
 */
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
vrna_Lfold_cb(const char                *string,
              int                       window_size,
              vrna_mfe_window_callback  *cb,
              void                      *data)
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
vrna_Lfoldz_cb(const char                       *string,
               int                              window_size,
               double                           min_z,
               vrna_mfe_window_zscore_callback  *cb,
               void                             *data)
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
vrna_aliLfold_cb(const char               **AS,
                 int                      window_size,
                 vrna_mfe_window_callback *cb,
                 void                     *data)
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


/*###########################################*/
/*# deprecated functions below              #*/
/*###########################################*/

#ifndef VRNA_DISABLE_BACKWARD_COMPATIBILITY

PUBLIC float
Lfold(const char  *string,
      char        *structure,
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


PUBLIC float
Lfoldz(const char *string,
       char       *structure,
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

#ifndef VRNA_WITH_SVM
  energy = vrna_mfe_window(vc, NULL);
#else
  energy = (zsc) ? vrna_mfe_window_zscore(vc, min_z, NULL) : vrna_mfe_window(vc, NULL);
#endif

  vrna_fold_compound_free(vc);

  return energy;
}


PUBLIC float
aliLfold(const char *AS[],
         char       *structure,
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
aliLfold_cb(const char                **AS,
            int                       maxdist,
            vrna_mfe_window_callback  *cb,
            void                      *data)
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
