/*
 *                minimum free energy
 *                RNA secondary structure prediction
 *                with maximum distance base pairs
 *
 *                c Ivo Hofacker, Peter Stadler
 *
 *                Vienna RNA package
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <string.h>
#include <limits.h>
#include "ViennaRNA/utils.h"
#include "ViennaRNA/energy_par.h"
#include "ViennaRNA/energy_const.h" /* defines MINPSCORE */
#include "ViennaRNA/fold_vars.h"
#include "ViennaRNA/params.h"
#include "ViennaRNA/loop_energies.h"
#include "ViennaRNA/gquad.h"
#include "ViennaRNA/ribo.h"
#include "ViennaRNA/alifold.h"
#include "ViennaRNA/fold.h"
#include "ViennaRNA/Lfold.h"
#include "ViennaRNA/aln_util.h"


#ifdef VRNA_WITH_SVM
#include "svm.h"
#include "svm_utils.h"
#endif

#ifdef __GNUC__
# define INLINE inline
#else
# define INLINE
#endif

#define MAXSECTORS                  500   /* dimension for a backtrack array */
#define INT_CLOSE_TO_UNDERFLOW(i)   ((i) <= (INT_MIN / 16))
#define UNDERFLOW_CORRECTION        (INT_MIN / 32)

#define NONE -10000 /* score for forbidden pairs */


typedef struct {
  FILE  *output;
  int   dangle_model;
} hit_data;

/*
 #################################
 # GLOBAL VARIABLES              #
 #################################
 */

/*
 #################################
 # PRIVATE VARIABLES             #
 #################################
 */

/*
 #################################
 # PRIVATE FUNCTION DECLARATIONS #
 #################################
 */
PRIVATE float wrap_Lfold(vrna_fold_compound_t             *vc,
                         int                              with_zsc,
                         double                           min_z,
                         vrna_mfe_window_callback         *cb,
#ifdef VRNA_WITH_SVM
                         vrna_mfe_window_zscore_callback  *cb_z,
#endif
                         void                             *data);


PRIVATE void  make_ptypes(vrna_fold_compound_t  *vc,
                          int                   i);


PRIVATE char *backtrack(vrna_fold_compound_t  *vc,
                        int                   start,
                        int                   maxdist);


PRIVATE int   fill_arrays(vrna_fold_compound_t            *vc,
                          int                             with_zsc,
                          double                          min_z,
#ifdef VRNA_WITH_SVM
                          struct svm_model                *avg_model,
                          struct svm_model                *sd_model,
#endif
                          int                             *underflow,
                          vrna_mfe_window_callback        *cb,
#ifdef VRNA_WITH_SVM
                          vrna_mfe_window_zscore_callback *cb_z,
#endif
                          void                            *data);


#ifdef VRNA_WITH_SVM
PRIVATE void
default_callback_z(int        start,
                   int        end,
                   const char *structure,
                   float      en,
                   float      zscore,
                   void       *data);


#endif

PRIVATE void
default_callback(int        start,
                 int        end,
                 const char *structure,
                 float      en,
                 void       *data);


PRIVATE double
cov_score(vrna_fold_compound_t  *fc,
          int                   i,
          int                   j,
          float                 **dm);


PRIVATE void  make_pscores(vrna_fold_compound_t *fc,
                           int                  start,
                           float                **dm);


PRIVATE int   fill_arrays_comparative(vrna_fold_compound_t      *fc,
                                      vrna_mfe_window_callback  *cb,
                                      void                      *data);


PRIVATE char *backtrack_comparative(vrna_fold_compound_t  *fc,
                                    int                   start,
                                    int                   maxdist);


PRIVATE void
default_callback_comparative(int        start,
                             int        end,
                             const char *structure,
                             float      en,
                             void       *data);


PRIVATE INLINE void allocate_dp_matrices(vrna_fold_compound_t *fc);


PRIVATE INLINE void free_dp_matrices(vrna_fold_compound_t *fc);


PRIVATE INLINE void rotate_dp_matrices(vrna_fold_compound_t *fc,
                                       int                  i);


PRIVATE INLINE void init_constraints(vrna_fold_compound_t *fc,
                                     float                **dm);


PRIVATE INLINE void rotate_constraints(vrna_fold_compound_t *fc,
                                       float                **dm,
                                       int                  i);


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
  vrna_md_t md;
  hit_data  data;

  vrna_md_set_default(&md);

  data.output       = (file) ? file : stdout;
  data.dangle_model = md.dangles;

  return vrna_Lfold_cb(string, window_size, &default_callback, (void *)&data);
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

  energy = wrap_Lfold(vc, 0, 0.0,
                      cb,
#ifdef VRNA_WITH_SVM
                      NULL,
#endif
                      data);

  vrna_fold_compound_free(vc);

  return energy;
}


PUBLIC float
vrna_mfe_window(vrna_fold_compound_t  *vc,
                FILE                  *file)
{
  hit_data data;

  data.output       = (file) ? file : stdout;
  data.dangle_model = vc->params->model_details.dangles;
  if (vc->type == VRNA_FC_TYPE_COMPARATIVE)
    return vrna_mfe_window_cb(vc, &default_callback_comparative, (void *)&data);
  else
    return vrna_mfe_window_cb(vc, &default_callback, (void *)&data);
}


PUBLIC float
vrna_mfe_window_cb(vrna_fold_compound_t     *vc,
                   vrna_mfe_window_callback *cb,
                   void                     *data)
{
  return wrap_Lfold(vc, 0, 0.0,
                    cb,
#ifdef VRNA_WITH_SVM
                    NULL,
#endif
                    data);
}


#ifdef VRNA_WITH_SVM

PUBLIC float
vrna_Lfoldz(const char  *string,
            int         window_size,
            double      min_z,
            FILE        *file)
{
  vrna_md_t md;
  hit_data  data;

  vrna_md_set_default(&md);

  data.output       = (file) ? file : stdout;
  data.dangle_model = md.dangles;

  return vrna_Lfoldz_cb(string, window_size, min_z, &default_callback_z, (void *)&data);
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

  energy = wrap_Lfold(vc, 1, min_z, NULL, cb, data);

  vrna_fold_compound_free(vc);

  return energy;
}


PUBLIC float
vrna_mfe_window_zscore(vrna_fold_compound_t *vc,
                       double               min_z,
                       FILE                 *file)
{
  hit_data data;

  data.output       = (file) ? file : stdout;
  data.dangle_model = vc->params->model_details.dangles;

  return vrna_mfe_window_zscore_cb(vc, min_z, &default_callback_z, (void *)&data);
}


PUBLIC float
vrna_mfe_window_zscore_cb(vrna_fold_compound_t            *vc,
                          double                          min_z,
                          vrna_mfe_window_zscore_callback *cb,
                          void                            *data)
{
  return wrap_Lfold(vc, 1, min_z, NULL, cb, data);
}


#endif


PUBLIC float
vrna_aliLfold(const char  **AS,
              int         maxdist,
              FILE        *fp)
{
  vrna_md_t md;
  hit_data  data;

  vrna_md_set_default(&md);

  data.output       = (fp) ? fp : stdout;
  data.dangle_model = md.dangles;

  return vrna_aliLfold_cb(AS, maxdist, &default_callback_comparative, (void *)&data);
}


PUBLIC float
vrna_aliLfold_cb(const char               **AS,
                 int                      maxdist,
                 vrna_mfe_window_callback *cb,
                 void                     *data)
{
  float                 en;
  vrna_fold_compound_t  *fc;
  vrna_md_t             md;

  vrna_md_set_default(&md);

  md.max_bp_span = md.window_size = maxdist;

  fc = vrna_fold_compound_comparative((const char **)AS, &md, VRNA_OPTION_MFE | VRNA_OPTION_WINDOW);

  en = (float)fill_arrays_comparative(fc, cb, data) / 100.;

  vrna_fold_compound_free(fc);

  return en;
}


/*
 #####################################
 # BEGIN OF STATIC HELPER FUNCTIONS  #
 #####################################
 */
PRIVATE float
wrap_Lfold(vrna_fold_compound_t             *vc,
           int                              with_zsc,
           double                           min_z,
           vrna_mfe_window_callback         *cb,
#ifdef VRNA_WITH_SVM
           vrna_mfe_window_zscore_callback  *cb_z,
#endif
           void                             *data)
{
  int               energy, underflow;
  float             mfe_local;

#ifdef VRNA_WITH_SVM
  struct svm_model  *avg_model  = NULL;
  struct svm_model  *sd_model   = NULL;
#endif

  if (vc->type == VRNA_FC_TYPE_COMPARATIVE)
    return (float)fill_arrays_comparative(vc, cb, data) / 100.;

  if (!vrna_fold_compound_prepare(vc, VRNA_OPTION_MFE | VRNA_OPTION_WINDOW)) {
    vrna_message_warning("vrna_mfe_window@Lfold.c: Failed to prepare vrna_fold_compound");
    return (float)(INF / 100.);
  }

#ifdef VRNA_WITH_SVM  /*svm*/
  if (with_zsc) {
    avg_model = svm_load_model_string(avg_model_string);
    sd_model  = svm_load_model_string(sd_model_string);
  }

#endif

  /* keep track of how many times we were close to an integer underflow */
  underflow = 0;

#ifdef VRNA_WITH_SVM
  energy = fill_arrays(vc, with_zsc, min_z, avg_model, sd_model, &underflow, cb, cb_z, data);
  if (with_zsc) {
    svm_free_model_content(avg_model);
    svm_free_model_content(sd_model);
  }

#else
  energy = fill_arrays(vc, with_zsc, min_z, &underflow, cb, data);
#endif

  mfe_local = (underflow > 0) ? ((float)underflow * (float)(UNDERFLOW_CORRECTION)) / 100. : 0.;
  mfe_local += (float)energy / 100.;

  return mfe_local;
}


PRIVATE INLINE void
allocate_dp_matrices(vrna_fold_compound_t *fc)
{
  int       i, j, length, maxdist, **c, **fML;
  vrna_hc_t *hc;

  length  = fc->length;
  maxdist = MIN2(fc->window_size, length);
  hc      = fc->hc;
  c       = fc->matrices->c_local;
  fML     = fc->matrices->fML_local;

  /* reserve additional memory for j-dimension */
  for (i = length; (i > length - maxdist - 5) && (i >= 0); i--) {
    c[i]                = (int *)vrna_alloc(sizeof(int) * (maxdist + 5));
    fML[i]              = (int *)vrna_alloc(sizeof(int) * (maxdist + 5));
    hc->matrix_local[i] = (char *)vrna_alloc(sizeof(char) * (maxdist + 5));
    if (fc->type == VRNA_FC_TYPE_SINGLE)
      fc->ptype_local[i] = vrna_alloc(sizeof(char) * (maxdist + 5));
    else if (fc->type == VRNA_FC_TYPE_COMPARATIVE)
      fc->pscore_local[i] = vrna_alloc(sizeof(int) * (maxdist + 5));
  }

  if (fc->type == VRNA_FC_TYPE_SINGLE)
    for (j = length; j > length - maxdist - 4; j--)
      for (i = (length - maxdist - 4 > 0) ? length - maxdist - 4 : 1; i < j; i++)
        c[i][j - i] = fML[i][j - i] = INF;
  else if (fc->type == VRNA_FC_TYPE_COMPARATIVE)
    for (j = length; j > length - maxdist - 3; j--)
      for (i = (length - maxdist - 2 > 0) ? length - maxdist - 2 : 1; i < j; i++)
        c[i][j - i] = fML[i][j - i] = INF;
}


PRIVATE INLINE void
free_dp_matrices(vrna_fold_compound_t *fc)
{
  int       i, length, maxdist, **c, **fML, **ggg, with_gquad;
  vrna_hc_t *hc;

  length      = fc->length;
  maxdist     = MIN2(fc->window_size, length);
  hc          = fc->hc;
  c           = fc->matrices->c_local;
  fML         = fc->matrices->fML_local;
  ggg         = fc->matrices->ggg_local;
  with_gquad  = fc->params->model_details.gquad;


  /* free additional memory for j-dimension */
  for (i = 0; (i < maxdist + 5) && (i <= length); i++) {
    if (fc->type == VRNA_FC_TYPE_SINGLE) {
      free(fc->ptype_local[i]);
      fc->ptype_local[i] = NULL;
    } else if (fc->type == VRNA_FC_TYPE_COMPARATIVE) {
      free(fc->pscore_local[i]);
      fc->pscore_local[i] = NULL;
    }

    free(c[i]);
    c[i] = NULL;
    free(fML[i]);
    fML[i] = NULL;
    free(hc->matrix_local[i]);
    hc->matrix_local[i] = NULL;
  }

  if (with_gquad) {
    for (i = 0; (i <= maxdist + 5) && (i <= length); i++)
      free(ggg[i]);
    free(ggg);
    fc->matrices->ggg_local = NULL;
  }
}


PRIVATE INLINE void
rotate_dp_matrices(vrna_fold_compound_t *fc,
                   int                  i)
{
  int       ii, maxdist, length, **c, **fML;
  vrna_hc_t *hc;

  length  = fc->length;
  maxdist = fc->window_size;
  c       = fc->matrices->c_local;
  fML     = fc->matrices->fML_local;
  hc      = fc->hc;

  if (i + maxdist + 4 <= length) {
    c[i - 1]                          = c[i + maxdist + 4];
    c[i + maxdist + 4]                = NULL;
    fML[i - 1]                        = fML[i + maxdist + 4];
    fML[i + maxdist + 4]              = NULL;
    hc->matrix_local[i - 1]           = hc->matrix_local[i + maxdist + 4];
    hc->matrix_local[i + maxdist + 4] = NULL;

    if ((fc->params->model_details.gquad) && (i > 1))
      vrna_gquad_mx_local_update(fc, i - 1);

    for (ii = 0; ii < maxdist + 5; ii++) {
      c[i - 1][ii]    = INF;
      fML[i - 1][ii]  = INF;
    }
  }
}


PRIVATE INLINE void
init_constraints(vrna_fold_compound_t *fc,
                 float                **dm)
{
  int i, length, maxdist;

  length  = fc->length;
  maxdist = fc->window_size;

  switch (fc->type) {
    case VRNA_FC_TYPE_SINGLE:
      for (i = length; (i >= length - maxdist - 4) && (i > 0); i--) {
        make_ptypes(fc, i);
        vrna_hc_prepare(fc, i);
      }
      break;

    case VRNA_FC_TYPE_COMPARATIVE:
      for (i = length; (i >= length - maxdist - 4) && (i > 0); i--) {
        make_pscores(fc, i, dm);
        vrna_hc_prepare(fc, i);
      }
      break;
  }
}


PRIVATE INLINE void
rotate_constraints(vrna_fold_compound_t *fc,
                   float                **dm,
                   int                  i)
{
  int length, maxdist;

  length  = fc->length;
  maxdist = fc->window_size;

  switch (fc->type) {
    case VRNA_FC_TYPE_SINGLE:
      if (i + maxdist + 4 <= length) {
        fc->ptype_local[i - 1]            = fc->ptype_local[i + maxdist + 4];
        fc->ptype_local[i + maxdist + 4]  = NULL;
        if (i > 1) {
          make_ptypes(fc, i - 1);
          vrna_hc_prepare(fc, i - 1);
        }
      }

      break;

    case VRNA_FC_TYPE_COMPARATIVE:
      if (i + maxdist + 4 <= length) {
        fc->pscore_local[i - 1]           = fc->pscore_local[i + maxdist + 4];
        fc->pscore_local[i + maxdist + 4] = NULL;
        if (i > 1) {
          make_pscores(fc, i - 1, dm);
          vrna_hc_prepare(fc, i - 1);
        }
      }

      break;
  }
}


PRIVATE int
fill_arrays(vrna_fold_compound_t            *vc,
            int                             zsc,
            double                          min_z,
#ifdef VRNA_WITH_SVM
            struct svm_model                *avg_model,
            struct svm_model                *sd_model,
#endif
            int                             *underflow,
            vrna_mfe_window_callback        *cb,
#ifdef VRNA_WITH_SVM
            vrna_mfe_window_zscore_callback *cb_z,
#endif
            void                            *data)
{
  /* fill "c", "fML" and "f3" arrays and return  optimal energy */

  char          **ptype, *prev, hc_decompose;
  short         *S;
  int           i, j, length, energy, maxdist, **c, **fML, *f3, no_close,
                type, with_gquad, dangle_model, noLP, noGUclosure, turn, fij,
                lind, *cc, *cc1, *Fmi, *DMLi, *DMLi1, *DMLi2, do_backtrack,
                prev_i, new_c, stackEnergy;
  double        prevz;
  vrna_param_t  *P;
  vrna_md_t     *md;
  vrna_hc_t     *hc;

  length        = vc->length;
  S             = vc->sequence_encoding2;
  ptype         = vc->ptype_local;
  maxdist       = vc->window_size;
  P             = vc->params;
  md            = &(P->model_details);
  dangle_model  = md->dangles;
  with_gquad    = md->gquad;
  noLP          = md->noLP;
  noGUclosure   = md->noGUclosure;
  turn          = md->min_loop_size;
  hc            = vc->hc;
  do_backtrack  = 0;
  prev_i        = 0;
  prev          = NULL;
  prevz         = 0.;

  c   = vc->matrices->c_local;
  fML = vc->matrices->fML_local;
  f3  = vc->matrices->f3_local;

  cc    = (int *)vrna_alloc(sizeof(int) * (maxdist + 5)); /* linear array for calculating canonical structures */
  cc1   = (int *)vrna_alloc(sizeof(int) * (maxdist + 5)); /*   "     "        */
  Fmi   = (int *)vrna_alloc(sizeof(int) * (maxdist + 5)); /* holds row i of fML (avoids jumps in memory) */
  DMLi  = (int *)vrna_alloc(sizeof(int) * (maxdist + 5)); /* DMLi[j] holds MIN(fML[i,k]+fML[k+1,j])  */
  DMLi1 = (int *)vrna_alloc(sizeof(int) * (maxdist + 5)); /*             MIN(fML[i+1,k]+fML[k+1,j])  */
  DMLi2 = (int *)vrna_alloc(sizeof(int) * (maxdist + 5)); /*             MIN(fML[i+2,k]+fML[k+1,j])  */

  /* reserve additional memory for j-dimension */
  allocate_dp_matrices(vc);

  init_constraints(vc, NULL);

  for (j = 0; j < maxdist + 5; j++)
    Fmi[j] = DMLi[j] = DMLi1[j] = DMLi2[j] = INF;

  if (with_gquad)
    vrna_gquad_mx_local_update(vc, length - maxdist - 4);

  for (i = length - turn - 1; i >= 1; i--) {
    /* i,j in [1..length] */
    for (j = i + turn + 1; j <= length && j <= i + maxdist; j++) {
      hc_decompose  = hc->matrix_local[i][j - i];
      type          = ptype[i][j - i];

      no_close = (((type == 3) || (type == 4)) && noGUclosure);

      if (hc_decompose) {
        /* we have a pair */
        new_c       = INF;
        stackEnergy = INF;

        if (!no_close) {
          /* check for hairpin loop */
          energy  = vrna_E_hp_loop(vc, i, j);
          new_c   = MIN2(new_c, energy);

          /* check for multibranch loops */
          energy  = vrna_E_mb_loop_fast(vc, i, j, DMLi1, DMLi2);
          new_c   = MIN2(new_c, energy);
        }

        if (dangle_model == 3) {
          /* coaxial stacking */
          energy  = E_mb_loop_stack(i, j, vc);
          new_c   = MIN2(new_c, energy);
        }

        /* check for interior loops */
        energy  = vrna_E_int_loop(vc, i, j);
        new_c   = MIN2(new_c, energy);

        /* remember stack energy for --noLP option */
        if (noLP) {
          stackEnergy = vrna_E_stack(vc, i, j);
          new_c       = MIN2(new_c, cc1[j - 1 - (i + 1)] + stackEnergy);
          cc[j - i]   = new_c;
          c[i][j - i] = cc1[j - 1 - (i + 1)] + stackEnergy;
        } else {
          c[i][j - i] = new_c;
        }
      } /* end >> if (pair) << */
      else {
        c[i][j - i] = INF;
      }

      /* done with c[i,j], now compute fML[i,j] */
      /* free ends ? -----------------------------------------*/
      fML[i][j - i] = vrna_E_ml_stems_fast(vc, i, j, Fmi, DMLi);
    } /* for (j...) */

    /* calculate energies of 5' and 3' fragments */
    {
      char *ss = NULL;

      f3[i] = vrna_E_ext_loop_3(vc, i);

      /* backtrack partial structure */
      if (f3[i] < f3[i + 1]) {
        do_backtrack = 1;
      } else if (do_backtrack) {
        int pairpartner; /*i+1?? is paired with pairpartner*/
        fij   = f3[i + 1];
        lind  = i + 1;
        /*start "short" backtrack*/

        /*get paired base*/
        while (fij == f3[lind + 1])
          lind++;

        pairpartner = vrna_BT_ext_loop_f3_pp(vc, lind, maxdist - (lind - (i + 1)), fij);
        if (pairpartner == -1)
          vrna_message_error("backtrack failed in short backtrack 1");

        if (zsc) {
#ifdef VRNA_WITH_SVM
          int     info_avg;
          double  average_free_energy;
          double  sd_free_energy;
          double  my_z;
          int     *AUGC = get_seq_composition(S, lind - 1, MIN2((pairpartner + 1), length), length);
          /*\svm*/
          average_free_energy = avg_regression(AUGC[0],
                                               AUGC[1],
                                               AUGC[2],
                                               AUGC[3],
                                               AUGC[4],
                                               avg_model,
                                               &info_avg);
          if (info_avg == 0) {
            double  difference;
            double  min_sd = minimal_sd(AUGC[0], AUGC[1], AUGC[2], AUGC[3], AUGC[4]);
            difference = (fij - f3[pairpartner + 1]) / 100. - average_free_energy;
            if (difference - (min_z * min_sd) <= 0.0001) {
              sd_free_energy =
                sd_regression(AUGC[0], AUGC[1], AUGC[2], AUGC[3], AUGC[4], sd_model);
              my_z = difference / sd_free_energy;
              if (my_z <= min_z) {
                ss = backtrack(vc, lind, pairpartner + 1);
                if (prev) {
                  if ((i + strlen(ss) < prev_i + strlen(prev)) ||
                      strncmp(ss + prev_i - i, prev, strlen(prev))) /* ss does not contain prev */
                    cb_z(prev_i, prev_i + strlen(prev) - 1, prev,
                         (f3[prev_i] - f3[prev_i + strlen(prev) - 1]) / 100., prevz, data);

                  free(prev);
                }

                prev    = ss;
                prev_i  = lind;
                prevz   = my_z;
              }
            }
          }

          free(AUGC);
          do_backtrack = 0;
#endif
        } else {
          /* original code for Lfold*/
          ss = backtrack(vc, lind, pairpartner + 1);
          if (prev) {
            if ((i + strlen(ss) < prev_i + strlen(prev)) ||
                strncmp(ss + prev_i - i, prev, strlen(prev)))
              /* ss does not contain prev */
              cb(prev_i, prev_i + strlen(prev) - 1, prev,
                 (f3[prev_i] - f3[prev_i + strlen(prev) - 1]) / 100., data);

            free(prev);
          }

          prev          = ss;
          prev_i        = lind;
          do_backtrack  = 0;
        }
      }

      if (i == 1) {
        if (prev) {
#ifdef VRNA_WITH_SVM
          if (zsc)
            cb_z(prev_i, prev_i + strlen(prev) - 1, prev,
                 (f3[prev_i] - f3[prev_i + strlen(prev) - 1]) / 100., prevz, data);
          else
            cb(prev_i, prev_i + strlen(prev) - 1, prev, (f3[prev_i] - f3[prev_i + strlen(
                                                                           prev) - 1]) / 100.,
               data);

#else
          cb(prev_i, prev_i + strlen(prev) - 1, prev, (f3[prev_i] - f3[prev_i + strlen(
                                                                         prev) - 1]) / 100., data);
#endif
          free(prev);
          prev = NULL;
        } else if ((f3[i] < 0) && (!zsc)) {
          do_backtrack = 1;
        }

        if (do_backtrack) {
          int     pairpartner; /*i+1?? is paired with pairpartner*/
          double  average_free_energy;
          double  sd_free_energy;
          int     info_avg;
          double  my_z;
          fij   = f3[i];
          lind  = i;
          while (fij == f3[lind + 1])
            lind++;

          pairpartner = vrna_BT_ext_loop_f3_pp(vc, lind, maxdist - (lind - i), fij);
          if (pairpartner == -1)
            vrna_message_error("backtrack failed in short backtrack 2");

          if (zsc) {
#ifdef VRNA_WITH_SVM
            int *AUGC = get_seq_composition(S, lind - 1, MIN2((pairpartner + 1), length), length);
            average_free_energy = avg_regression(AUGC[0],
                                                 AUGC[1],
                                                 AUGC[2],
                                                 AUGC[3],
                                                 AUGC[4],
                                                 avg_model,
                                                 &info_avg);
            if (info_avg == 0) {
              double  difference;
              double  min_sd = minimal_sd(AUGC[0], AUGC[1], AUGC[2], AUGC[3], AUGC[4]);
              difference = (fij - f3[pairpartner + 1]) / 100. - average_free_energy;
              if (difference - (min_z * min_sd) <= 0.0001) {
                sd_free_energy = sd_regression(AUGC[0],
                                               AUGC[1],
                                               AUGC[2],
                                               AUGC[3],
                                               AUGC[4],
                                               sd_model);
                my_z = difference / sd_free_energy;
                if (my_z <= min_z) {
                  ss = backtrack(vc, lind, pairpartner + 1);
                  cb_z(lind, lind + strlen(ss) - 1, ss, (f3[lind] - f3[lind + strlen(
                                                                         ss) - 1]) / 100., my_z,
                       data);
                }
              }
            }

            free(AUGC);
#endif
          } else {
            ss = backtrack(vc, lind, pairpartner + 1);
            cb(1, lind + strlen(ss) - 1, ss, (f3[lind] - f3[lind + strlen(ss) - 1]) / 100., data);
            free(ss);
          }
        }

        do_backtrack = 0;
      }
    }
    {
      int *FF; /* rotate the auxilliary arrays */

      /* check for values close to integer underflow */
      if (INT_CLOSE_TO_UNDERFLOW(f3[i])) {
        /* correct f3 free energies and increase underflow counter */
        int cnt;
        for (cnt = i; cnt <= length && cnt <= lind + maxdist + 2; cnt++)
          f3[cnt] -= UNDERFLOW_CORRECTION;
        (*underflow)++;
      }

      FF    = DMLi2;
      DMLi2 = DMLi1;
      DMLi1 = DMLi;
      DMLi  = FF;
      FF    = cc1;
      cc1   = cc;
      cc    = FF;
      for (j = 0; j < maxdist + 5; j++)
        cc[j] = Fmi[j] = DMLi[j] = INF;

      rotate_dp_matrices(vc, i);
      rotate_constraints(vc, NULL, i);
    }
  }

  free(cc);
  free(cc1);
  free(Fmi);
  free(DMLi);
  free(DMLi1);
  free(DMLi2);

  free_dp_matrices(vc);

  return f3[1];
}


PRIVATE char *
backtrack(vrna_fold_compound_t  *vc,
          int                   start,
          int                   maxdist)
{
  /*------------------------------------------------------------------
   *  trace back through the "c", "f3" and "fML" arrays to get the
   *  base pairing list. No search for equivalent structures is done.
   *  This is fast, since only few structure elements are recalculated.
   *  ------------------------------------------------------------------*/
  sect            sector[MAXSECTORS]; /* backtracking sectors */
  char            *string, *structure, **ptype;
  int             i, j, k, length, no_close, type, s, b, bt_type, turn,
                  dangle_model, noLP, noGUclosure, **c, dangle3, ml, cij,
                  canonical, p, q, comp1, comp2, max3;
  vrna_param_t    *P;
  vrna_md_t       *md;
  vrna_bp_stack_t *bp_stack;

  string        = vc->sequence;
  length        = vc->length;
  ptype         = vc->ptype_local;
  P             = vc->params;
  md            = &(P->model_details);
  dangle_model  = md->dangles;
  noLP          = md->noLP;
  noGUclosure   = md->noGUclosure;
  bt_type       = md->backtrack_type;
  turn          = md->min_loop_size;

  s         = 0;                                                                                /* depth of backtracking stack */
  b         = 0;                                                                                /* number of base pairs */
  bp_stack  = (vrna_bp_stack_t *)vrna_alloc(sizeof(vrna_bp_stack_t) * (4 * (1 + length / 2)));  /* add a guess of how many G's may be involved in a G quadruplex */

  c = vc->matrices->c_local;

  sector[++s].i = start;
  sector[s].j   = MIN2(length, maxdist + 1);
  sector[s].ml  = (bt_type == 'M') ? 1 : ((bt_type == 'C') ? 2 : 0);

  structure = (char *)vrna_alloc((MIN2(length - start, maxdist) + 3) * sizeof(char));

  memset(structure, '.', MIN2(length - start, maxdist) + 1);

  dangle3 = 0;

  while (s > 0) {
    canonical = 1;     /* (i,j) closes a canonical structure */

    /* pop one element from stack */
    i   = sector[s].i;
    j   = sector[s].j;
    ml  = sector[s--].ml;  /* ml is a flag indicating if backtracking is to
                           * occur in the fML- (1) or in the f-array (0) */

    if (j < i + turn + 1)
      continue;                     /* no more pairs in this interval */

    switch (ml) {
      /* backtrack in f3 */
      case 0:
        if (vrna_BT_ext_loop_f3(vc, &i, j, &p, &q, bp_stack, &b)) {
          if (i > 0) {
            sector[++s].i = i;
            sector[s].j   = j;
            sector[s].ml  = 0;
          }

          if (p > 0) {
            if (((i == q + 2) || (dangle_model == 2)) && (q < length))
              dangle3 = 1;

            i = p;
            j = q;
            goto repeat1;
          }

          continue;
        } else {
          vrna_message_error("backtracking failed in f3 for sequence:\n%s\n", string);
        }

      /* trace back in fML array */
      case 1:
        if (vrna_BT_mb_loop_split(vc, &i, &j, &p, &q, &comp1, &comp2, bp_stack, &b)) {
          if (i > 0) {
            sector[++s].i = i;
            sector[s].j   = j;
            sector[s].ml  = comp1;
          }

          if (p > 0) {
            sector[++s].i = p;
            sector[s].j   = q;
            sector[s].ml  = comp2;
          }

          continue;
        } else {
          vrna_message_error("backtracking failed in fML for sequence:\n%s\n", string);
        }

        break;

      /* backtrack in c */
      case 2:
        bp_stack[++b].i = i;
        bp_stack[b].j   = j;
        goto repeat1;

      default:
        vrna_message_error("Backtracking failed due to unrecognized DP matrix!");
        break;
    }

repeat1:

    /*----- begin of "repeat:" -----*/
    if (canonical)
      cij = c[i][j - i];

    type = ptype[i][j - i];
    if (type == 0)
      type = 7;

    if (noLP) {
      if (vrna_BT_stack(vc, &i, &j, &cij, bp_stack, &b)) {
        canonical = 0;
        goto repeat1;
      }
    }

    canonical = 1;

    no_close = (((type == 3) || (type == 4)) && noGUclosure);

    if (no_close) {
      if (cij == FORBIDDEN)
        continue;
    } else {
      if (vrna_BT_hp_loop(vc, i, j, cij, bp_stack, &b))
        continue;
    }

    if (vrna_BT_int_loop(vc, &i, &j, cij, bp_stack, &b)) {
      if (i < 0)
        continue;
      else
        goto repeat1;
    }

    /* (i.j) must close a multi-loop */
    if (vrna_BT_mb_loop(vc, &i, &j, &k, cij, &comp1, &comp2)) {
      sector[++s].i = i;
      sector[s].j   = k;
      sector[s].ml  = comp1;
      sector[++s].i = k + 1;
      sector[s].j   = j;
      sector[s].ml  = comp2;
    } else {
      vrna_message_error("backtracking failed in repeat for sequence:\n%s\n", string);
    }

    /* end of repeat: --------------------------------------------------*/
  } /* end of infinite while loop */


  bp_stack[0].i = b;

  /* and now create a dot-brakcet string from the base pair stack... */
  max3 = 1;
  for (i = 1; i <= b; i++) {
    if (bp_stack[i].i == bp_stack[i].j) {
      /* Gquad bonds are marked as bp[i].i == bp[i].j */
      structure[bp_stack[i].i - start] = '+';
    } else {
      /* the following ones are regular base pairs */
      structure[bp_stack[i].i - start]  = '(';
      structure[bp_stack[i].j - start]  = ')';
    }

    if (max3 < bp_stack[i].j - start)
      max3 = bp_stack[i].j - start;
  }

  free(bp_stack);

  structure = (char *)vrna_realloc(structure,
                                   sizeof(char) * (max3 + dangle3 + 2));
  structure[max3 + dangle3 + 1] = '\0';

  return structure;
}


PRIVATE int
fill_arrays_comparative(vrna_fold_compound_t      *fc,
                        vrna_mfe_window_callback  *cb,
                        void                      *data)
{
  /* fill "c", "fML" and "f3" arrays and return  optimal energy */
  short         **S;
  char          **strings, *prev;
  int           **pscore, i, j, length, energy, turn,
                n_seq, **c, **fML, *f3, *cc, *cc1, *Fmi,
                *DMLi, *DMLi1, *DMLi2, maxdist, prev_i, prev_j,
                prev_en, do_backtrack, with_gquad, new_c, psc,
                stackEnergy, jjj, iii, eee;
  float         **dm;
  vrna_mx_mfe_t *matrices;
  vrna_param_t  *P;
  vrna_md_t     *md;
  vrna_hc_t     *hc;

  int           olddm[7][7] = { { 0, 0, 0, 0, 0, 0, 0 },   /* hamming distance between pairs PRIVATE needed??*/
                                { 0, 0, 2, 2, 1, 2, 2 } /* CG */,
                                { 0, 2, 0, 1, 2, 2, 2 } /* GC */,
                                { 0, 2, 1, 0, 2, 1, 2 } /* GU */,
                                { 0, 1, 2, 2, 0, 2, 1 } /* UG */,
                                { 0, 2, 2, 1, 2, 0, 2 } /* AU */,
                                { 0, 2, 2, 2, 1, 2, 0 } /* UA */ };

  do_backtrack  = 0;
  prev_i        = 0;
  prev_j        = 0;
  prev_en       = INF;

  dm    = NULL;
  prev  = NULL;

  strings     = fc->sequences;
  S           = fc->S;
  n_seq       = fc->n_seq;
  length      = fc->length;
  maxdist     = fc->window_size;
  matrices    = fc->matrices;
  hc          = fc->hc;
  P           = fc->params;
  md          = &(P->model_details);
  turn        = md->min_loop_size;
  with_gquad  = md->gquad;

  if (md->ribo) {
    if (RibosumFile != NULL)
      dm = readribosum(RibosumFile);
    else
      dm = get_ribosum((const char **)strings, n_seq, S[0][0]);
  } else {
    /*use usual matrix*/
    dm = (float **)vrna_alloc(7 * sizeof(float *));
    for (i = 0; i < 7; i++) {
      dm[i] = (float *)vrna_alloc(7 * sizeof(float));
      for (j = 0; j < 7; j++)
        dm[i][j] = (float)olddm[i][j];
    }
  }

  pscore = fc->pscore_local;      /* precomputed array of pair types */

  c   = fc->matrices->c_local;    /* energy array, given that i-j pair */
  fML = fc->matrices->fML_local;  /* multi-loop auxiliary energy array */
  f3  = fc->matrices->f3_local;   /* energy of 5' end */

  cc    = (int *)vrna_alloc(sizeof(int) * (maxdist + 5));
  cc1   = (int *)vrna_alloc(sizeof(int) * (maxdist + 5));
  Fmi   = (int *)vrna_alloc(sizeof(int) * (maxdist + 5));
  DMLi  = (int *)vrna_alloc(sizeof(int) * (maxdist + 5));
  DMLi1 = (int *)vrna_alloc(sizeof(int) * (maxdist + 5));
  DMLi2 = (int *)vrna_alloc(sizeof(int) * (maxdist + 5));

  /* reserve additional memory for j-dimension */
  allocate_dp_matrices(fc);

  init_constraints(fc, dm);

  for (j = 0; j < maxdist + 5; j++)
    Fmi[j] = DMLi[j] = DMLi1[j] = DMLi2[j] = INF;

  if (with_gquad)
    vrna_gquad_mx_local_update(fc, length - maxdist - 4);

  for (i = length - turn - 1; i >= 1; i--) {
    /* i,j in [1..length] */
    for (j = i + 1; j <= length && j <= i + turn; j++)
      c[i][j - i] = fML[i][j - i] = INF;
    for (j = i + turn + 1; j <= length && j <= i + maxdist; j++) {
      psc = pscore[i][j - i];

      if (hc->matrix_local[i][j - i]) {
        /* we evaluate this pair */
        new_c       = INF;
        stackEnergy = INF;

        /* hairpin ----------------------------------------------*/
        energy  = vrna_E_hp_loop(fc, i, j);
        new_c   = MIN2(new_c, energy);

        /* check for multibranch loops */
        energy  = vrna_E_mb_loop_fast(fc, i, j, DMLi1, NULL);
        new_c   = MIN2(new_c, energy);

        /* check for interior loops */
        energy  = vrna_E_int_loop(fc, i, j);
        new_c   = MIN2(new_c, energy);

        /* remember stack energy for --noLP option */
        if (md->noLP) {
          stackEnergy = vrna_E_stack(fc, i, j);
          new_c       = MIN2(new_c, cc1[j - 1 - (i + 1)] + stackEnergy);
          cc[j - i]   = new_c - psc; /* add covariance bonnus/penalty */
          c[i][j - i] = cc1[j - 1 - (i + 1)] + stackEnergy - psc;
        } else {
          c[i][j - i] = new_c - psc; /* add covariance bonnus/penalty */
        }
      } /* end >> if (pair) << */
      else {
        c[i][j - i] = INF;
      }

      /* done with c[i,j], now compute fML[i,j] */
      fML[i][j - i] = vrna_E_ml_stems_fast(fc, i, j, Fmi, DMLi);
    } /* for (j...) */

    /* calculate energies of 5' and 3' fragments */
    f3[i] = vrna_E_ext_loop_3(fc, i);

    if (f3[i] < f3[i + 1]) {
      do_backtrack = 1;
    } else if (do_backtrack) {
      eee = f3[i + 1];
      iii = i + 1;

      while (eee == f3[iii + 1])
        iii++;

      jjj = vrna_BT_ext_loop_f3_pp(fc, iii, maxdist - (iii - (i + 1)), eee);
      if (jjj == -1)
        vrna_message_error("backtrack failed in short backtrack 1 for columns %d to %d", iii, jjj);

      energy = f3[iii] - f3[jjj];

      /* do not spam output with relatively unstable structures */
      char *sss = backtrack_comparative(fc, iii, jjj - iii + 1);
      if (prev) {
        if ((jjj < prev_j) || strncmp(sss + prev_i - iii, prev, prev_j + 1 - prev_i + 1 - 1))   /* -1 because of 3' dangle unpaired position */
          cb(prev_i, prev_j + 1, prev, (prev_en) / (100. * n_seq), data);

        free(prev);
      }

      prev    = sss;
      prev_i  = iii;
      prev_j  = jjj;
      prev_en = energy;

      do_backtrack = 0;
    }

    if (i == 1) {
      if (f3[1] != f3[prev_i]) {
        /* find pairing partner, if any */
        eee = f3[1];
        iii = 1;

        while (eee == f3[iii + 1]) {
          iii++;
          if (iii > maxdist)
            break;
        }

        if (iii < maxdist) {
          /* otherwise, c-array columns don't exist anymore */
          jjj = vrna_BT_ext_loop_f3_pp(fc, iii, maxdist - (iii - 1), eee);
          if (jjj != -1) {
            energy = f3[iii] - f3[jjj];

            char *ss = backtrack_comparative(fc, iii, jjj - iii + 1);

            if (prev)
              if ((jjj < prev_j) || strncmp(ss + prev_i - iii, prev, prev_j + 1 - prev_i + 1 - 1))  /* -1 because of 3' dangle unpaired position */
                cb(prev_i, prev_j + 1, prev, (prev_en) / (100. * n_seq), data);

            /* execute callback */
            cb(iii, jjj + 1, ss, energy / (100. * n_seq), data);

            free(prev);
            prev = NULL;
            free(ss);
          }
        }
      }

      if (prev) {
        cb(prev_i, prev_j + 1, prev, (prev_en) / (100. * n_seq), data);
        free(prev);
        prev = NULL;
      }
    }

    {
      int *FF; /* rotate the auxilliary arrays */
      FF    = DMLi2;
      DMLi2 = DMLi1;
      DMLi1 = DMLi;
      DMLi  = FF;
      FF    = cc1;
      cc1   = cc;
      cc    = FF;
      for (j = 0; j < maxdist + 5; j++)
        cc[j] = Fmi[j] = DMLi[j] = INF;

      rotate_dp_matrices(fc, i);
      rotate_constraints(fc, dm, i);
    }
  }

  free(cc);
  free(cc1);
  free(Fmi);
  free(DMLi);
  free(DMLi1);
  free(DMLi2);

  for (i = 0; i < 7; i++)
    free(dm[i]);
  free(dm);

  free_dp_matrices(fc);

  return matrices->f3[1];
}


PRIVATE char *
backtrack_comparative(vrna_fold_compound_t  *fc,
                      int                   start,
                      int                   maxdist)
{
  /*------------------------------------------------------------------
   *  trace back through the "c", "f3" and "fML" arrays to get the
   *  base pairing list. No search for equivalent structures is done.
   *  This is fast, since only few structure elements are recalculated.
   *  ------------------------------------------------------------------*/
  sect            sector[MAXSECTORS]; /* backtracking sectors */

  char            *structure;
  int             **pscore, i, j, k, length, turn, dangle_model, noLP,
                  s, **c, ml, cij, p, q, canonical, b, dangle3, max3,
                  comp1, comp2;
  vrna_param_t    *P;
  vrna_md_t       *md;
  vrna_bp_stack_t *bp_stack;

  length        = fc->length;
  pscore        = fc->pscore_local; /* precomputed array of pair types */
  P             = fc->params;
  md            = &(P->model_details);
  dangle_model  = md->dangles;
  noLP          = md->noLP;

  c = fc->matrices->c_local;      /* energy array, given that i-j pair */

  turn = md->min_loop_size;

  s         = 0;                                                                                /* depth of backtracking stack */
  b         = 0;                                                                                /* number of base pairs */
  bp_stack  = (vrna_bp_stack_t *)vrna_alloc(sizeof(vrna_bp_stack_t) * (4 * (1 + length / 2)));  /* add a guess of how many G's may be involved in a G quadruplex */

  sector[++s].i = start;
  sector[s].j   = MIN2(length, start + maxdist + 1);
  sector[s].ml  = (backtrack_type == 'M') ? 1 : ((backtrack_type == 'C') ? 2 : 0);

  structure = (char *)vrna_alloc((MIN2(length - start, maxdist) + 3) * sizeof(char));

  memset(structure, '.', MIN2(length - start, maxdist) + 1);

  dangle3 = 0;

  while (s > 0) {
    canonical = 1;     /* (i,j) closes a canonical structure */

    i   = sector[s].i;
    j   = sector[s].j;
    ml  = sector[s--].ml;   /* ml is a flag indicating if backtracking is to
                            * occur in the fML- (1) or in the f-array (0) */
    if (j < i + turn + 1)
      continue;             /* no more pairs in this interval */

    switch (ml) {
      /* backtrack in f3 */
      case 0:
        if (vrna_BT_ext_loop_f3(fc, &i, j, &p, &q, bp_stack, &b)) {
          if (i > 0) {
            sector[++s].i = i;
            sector[s].j   = j;
            sector[s].ml  = 0;
          }

          if (p > 0) {
            if (((i == q + 2) || (dangle_model == 2)) && (q < length))
              dangle3 = 1;

            i = p;
            j = q;
            goto repeat1_comparative;
          }

          continue;
        } else {
          vrna_message_error("backtracking failed in f3\n");
        }

        break;

      /* trace back in fML array */
      case 1:
        if (vrna_BT_mb_loop_split(fc, &i, &j, &p, &q, &comp1, &comp2, bp_stack, &b)) {
          if (i > 0) {
            sector[++s].i = i;
            sector[s].j   = j;
            sector[s].ml  = comp1;
          }

          if (p > 0) {
            sector[++s].i = p;
            sector[s].j   = q;
            sector[s].ml  = comp2;
          }

          continue;
        } else {
          vrna_message_error("backtracking failed in fML\n");
        }

        break;

      /* backtrack in c */
      case 2:
        bp_stack[++b].i = i;
        bp_stack[b].j   = j;
        goto repeat1_comparative;
        break;

      default:
        vrna_message_error("Backtracking failed due to unrecognized DP matrix!");
        break;
    }

repeat1_comparative:

    /*----- begin of "repeat:" -----*/
    if (canonical)
      cij = c[i][j - i];

    if (noLP) {
      if (vrna_BT_stack(fc, &i, &j, &cij, bp_stack, &b)) {
        canonical = 0;
        goto repeat1_comparative;
      }
    }

    canonical = 1;
    cij       += pscore[i][j - i];

    if (vrna_BT_hp_loop(fc, i, j, cij, bp_stack, &b))
      continue;

    if (vrna_BT_int_loop(fc, &i, &j, cij, bp_stack, &b)) {
      if (i < 0)
        continue;
      else
        goto repeat1_comparative;
    }

    /* (i.j) must close a multi-loop */
    if (vrna_BT_mb_loop(fc, &i, &j, &k, cij, &comp1, &comp2)) {
      sector[++s].i = i;
      sector[s].j   = k;
      sector[s].ml  = comp1;
      sector[++s].i = k + 1;
      sector[s].j   = j;
      sector[s].ml  = comp2;
    } else {
      vrna_message_error("backtracking failed in repeat\n");
    }

    /* end of repeat: --------------------------------------------------*/
  }

  bp_stack[0].i = b;

  /* and now create a dot-brakcet string from the base pair stack... */
  max3 = 1;
  for (i = 1; i <= b; i++) {
    if (bp_stack[i].i == bp_stack[i].j) {
      /* Gquad bonds are marked as bp[i].i == bp[i].j */
      structure[bp_stack[i].i - start] = '+';
    } else {
      /* the following ones are regular base pairs */
      structure[bp_stack[i].i - start]  = '(';
      structure[bp_stack[i].j - start]  = ')';
    }

    if (max3 < bp_stack[i].j - start)
      max3 = bp_stack[i].j - start;
  }

  free(bp_stack);

  structure = (char *)vrna_realloc(structure,
                                   sizeof(char) * (max3 + dangle3 + 2));
  structure[max3 + dangle3 + 1] = '\0';

  return structure;
}


PRIVATE void
make_ptypes(vrna_fold_compound_t  *vc,
            int                   i)
{
  int       j, k, type, n, maxdist, turn, noLP;
  short     *S;
  char      **ptype;
  vrna_md_t *md;

  n       = (int)vc->length;
  S       = vc->sequence_encoding2;
  ptype   = vc->ptype_local;
  maxdist = vc->window_size;
  md      = &(vc->params->model_details);
  turn    = md->min_loop_size;
  noLP    = md->noLP;

  for (k = turn + 1; k < maxdist; k++) {
    j = i + k;
    if (j > n)
      break;

    type = md->pair[S[i]][S[j]];

    if (noLP && type) {
      if (!ptype[i + 1][j - 1 - i - 1])
        if (j == n || i == 1 || (!md->pair[S[i - 1]][S[j + 1]]))
          type = 0;
    }

    ptype[i][j - i] = type;
  }
}


PRIVATE double
cov_score(vrna_fold_compound_t  *fc,
          int                   i,
          int                   j,
          float                 **dm)
{
  char      **AS;
  short     **S;
  int       n_seq, k, l, s, type;
  double    score;
  vrna_md_t *md;
  int       pfreq[8] = {
    0, 0, 0, 0, 0, 0, 0, 0
  };

  n_seq = fc->n_seq;
  AS    = fc->sequences;
  S     = fc->S;
  md    = &(fc->params->model_details);

  for (s = 0; s < n_seq; s++) {
    if (S[s][i] == 0 && S[s][j] == 0) {
      type = 7;                             /* gap-gap  */
    } else {
      if ((AS[s][i] == '~') || (AS[s][j] == '~'))
        type = 7;
      else
        type = md->pair[S[s][i]][S[s][j]];
    }

    pfreq[type]++;
  }

  if (pfreq[0] * 2 + pfreq[7] > n_seq) {
    return NONE;
  } else {
    for (k = 1, score = 0.; k <= 6; k++) /* ignore pairtype 7 (gap-gap) */
      for (l = k; l <= 6; l++)
        /* scores for replacements between pairtypes    */
        /* consistent or compensatory mutations score 1 or 2  */
        score += pfreq[k] * pfreq[l] * dm[k][l];
  }

  /* counter examples score -1, gap-gap scores -0.25   */
  return md->cv_fact * ((UNIT * score) / n_seq - md->nc_fact * UNIT * (pfreq[0] + pfreq[7] * 0.25));
}


PRIVATE void
make_pscores(vrna_fold_compound_t *fc,
             int                  i,
             float                **dm)
{
  /* calculate co-variance bonus for each pair depending on  */
  /* compensatory/consistent mutations and incompatible seqs */
  /* should be 0 for conserved pairs, >0 for good pairs      */
  int       n, j, **pscore, maxd, turn, noLP;
  vrna_md_t *md;

  n       = (int)fc->length;
  maxd    = fc->window_size;
  pscore  = fc->pscore_local;
  md      = &(fc->params->model_details);
  turn    = md->min_loop_size;
  noLP    = md->noLP;

  /*fill pscore[start], too close*/
  for (j = i + 1; (j < i + turn + 1) && (j <= n); j++)
    pscore[i][j - i] = NONE;
  for (j = i + turn + 1; ((j <= n) && (j <= i + maxd)); j++)
    pscore[i][j - i] = cov_score(fc, i, j, dm);

  if (noLP) {
    /* remove unwanted lonely pairs */
    int otype = 0, ntype = 0;
    for (j = i + turn; ((j < n) && (j < i + maxd)); j++) {
      if ((i > 1) && (j < n))
        otype = cov_score(fc, i - 1, j + 1, dm);

      if (i < n)
        ntype = pscore[i + 1][j - 1 - (i + 1)];
      else
        ntype = NONE;

      if ((otype < -4 * UNIT) && (ntype < -4 * UNIT)) /* worse than 2 counterex */
        pscore[i][j - i] = NONE;                      /* i.j can only form isolated pairs */
    }
  }

  if ((j - i + 1) > maxd)
    pscore[i][j - i] = NONE;
}


PRIVATE void
default_callback(int        start,
                 int        end,
                 const char *structure,
                 float      en,
                 void       *data)
{
  FILE  *output       = ((hit_data *)data)->output;
  int   dangle_model  = ((hit_data *)data)->dangle_model;

  if ((dangle_model == 2) && (start > 1))
    fprintf(output, ".%s (%6.2f) %4d\n", structure, en, start - 1);
  else
    fprintf(output, "%s (%6.2f) %4d\n ", structure, en, start);
}


#ifdef VRNA_WITH_SVM
PRIVATE void
default_callback_z(int        start,
                   int        end,
                   const char *structure,
                   float      en,
                   float      zscore,
                   void       *data)
{
  FILE  *output       = ((hit_data *)data)->output;
  int   dangle_model  = ((hit_data *)data)->dangle_model;

  if ((dangle_model == 2) && (start > 1))
    fprintf(output, ".%s (%6.2f) %4d z= %.3f\n", structure, en, start - 1, zscore);
  else
    fprintf(output, "%s (%6.2f) %4d z= %.3f\n ", structure, en, start, zscore);
}


#endif


PRIVATE void
default_callback_comparative(int        start,
                             int        end,
                             const char *structure,
                             float      en,
                             void       *data)
{
  if (csv == 1)
    printf("%s ,%6.2f, %4d, %4d\n", structure, en, start, end);
  else
    printf("%s (%6.2f) %4d - %4d\n", structure, en, start, end);
}


/*###########################################*/
/*# deprecated functions below              #*/
/*###########################################*/

#ifdef  VRNA_BACKWARD_COMPAT

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
aliLfold(const char *strings[],
         char       *structure,
         int        maxdist)
{
  vrna_md_t md;
  hit_data  data;

  set_model_details(&md);

  data.output       = stdout;
  data.dangle_model = md.dangles;

  return aliLfold_cb(strings, maxdist, &default_callback_comparative, (void *)&data);
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

  fc = vrna_fold_compound_comparative((const char **)AS, &md, VRNA_OPTION_MFE | VRNA_OPTION_WINDOW);

  en = vrna_mfe_window_cb(fc, cb, data);

  vrna_fold_compound_free(fc);

  return en;
}


#endif
