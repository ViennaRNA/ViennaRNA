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

#define MAXSECTORS                  500   /* dimension for a backtrack array */
#define INT_CLOSE_TO_UNDERFLOW(i)   ((i) <= (INT_MIN / 16))
#define UNDERFLOW_CORRECTION        (INT_MIN / 32)
#define UNIT 100
#define MINPSCORE -2 * UNIT


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
  int               i, energy, underflow, n, maxdist;
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

  n       = vc->length;
  maxdist = vc->window_size;

  /* reserve additional memory for j-dimension */
  for (i = (int)vc->length; (i > (int)vc->length - vc->window_size - 5) && (i >= 0); i--) {
    vc->ptype_local[i]          = vrna_alloc(sizeof(char) * (vc->window_size + 5));
    vc->matrices->c_local[i]    = (int *)vrna_alloc(sizeof(int) * (vc->window_size + 5));
    vc->matrices->fML_local[i]  = (int *)vrna_alloc(sizeof(int) * (vc->window_size + 5));
  }

  for (i = n; (i >= (int)n - (int)maxdist - 4) && (i > 0); i--)
    make_ptypes(vc, i);

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

  /* free additional memory for j-dimension */
  for (i = 0; (i < vc->window_size + 5) && (i <= (int)vc->length); i++) {
    free(vc->ptype_local[i]);
    free(vc->matrices->c_local[i]);
    free(vc->matrices->fML_local[i]);
  }

  return mfe_local;
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

  int           i, j, k, length, energy, maxdist;
  int           **c, **fML, *f3, **ggg;
  int           decomp, new_fML;
  int           no_close, type, type_2, tt, with_gquad, dangle_model, noLP, noGUclosure, turn;
  int           *rtype;
  int           fij;
  int           lind;

  int           *cc     = NULL; /* linear array for calculating canonical structures */
  int           *cc1    = NULL; /*   "     "        */
  int           *Fmi    = NULL; /* holds row i of fML (avoids jumps in memory) */
  int           *DMLi   = NULL; /* DMLi[j] holds MIN(fML[i,k]+fML[k+1,j])  */
  int           *DMLi1  = NULL; /*             MIN(fML[i+1,k]+fML[k+1,j])  */
  int           *DMLi2  = NULL; /*             MIN(fML[i+2,k]+fML[k+1,j])  */

  short         *S, *S1;
  char          *string, **ptype, *prev;
  double        prevz         = 0.;
  int           do_backtrack  = 0, prev_i = 0;
  vrna_param_t  *P;
  vrna_md_t     *md;


  string        = vc->sequence;
  length        = vc->length;
  S             = vc->sequence_encoding2;
  S1            = vc->sequence_encoding;
  ptype         = vc->ptype_local;
  maxdist       = vc->window_size;
  P             = vc->params;
  md            = &(P->model_details);
  dangle_model  = md->dangles;
  with_gquad    = md->gquad;
  noLP          = md->noLP;
  noGUclosure   = md->noGUclosure;
  turn          = md->min_loop_size;
  rtype         = &(md->rtype[0]);

  prev = NULL;

  c   = vc->matrices->c_local;
  fML = vc->matrices->fML_local;
  f3  = vc->matrices->f3_local;
  ggg = vc->matrices->ggg_local;

  cc    = (int *)vrna_alloc(sizeof(int) * (maxdist + 5));
  cc1   = (int *)vrna_alloc(sizeof(int) * (maxdist + 5));
  Fmi   = (int *)vrna_alloc(sizeof(int) * (maxdist + 5));
  DMLi  = (int *)vrna_alloc(sizeof(int) * (maxdist + 5));
  DMLi1 = (int *)vrna_alloc(sizeof(int) * (maxdist + 5));
  DMLi2 = (int *)vrna_alloc(sizeof(int) * (maxdist + 5));


  for (j = 0; j < maxdist + 5; j++)
    Fmi[j] = DMLi[j] = DMLi1[j] = DMLi2[j] = INF;
  for (j = length; j > length - maxdist - 4; j--)
    for (i = (length - maxdist - 4 > 0) ? length - maxdist - 4 : 1; i < j; i++)
      c[i][j - i] = fML[i][j - i] = INF;

  if (with_gquad) {
    vrna_gquad_mx_local_update(vc, length - maxdist - 4);
    ggg = vc->matrices->ggg_local;
  }

  for (i = length - turn - 1; i >= 1; i--) {
    /* i,j in [1..length] */
    for (j = i + turn + 1; j <= length && j <= i + maxdist; j++) {
      int p, q;
      type = ptype[i][j - i];

      no_close = (((type == 3) || (type == 4)) && noGUclosure);

      if (type) {
        /* we have a pair */
        int new_c = 0, stackEnergy = INF;
        /* hairpin ----------------------------------------------*/

        new_c = (no_close) ? FORBIDDEN : E_Hairpin(j - i - 1,
                                                   type,
                                                   S1[i + 1],
                                                   S1[j - 1],
                                                   string + i - 1,
                                                   P);

        /*--------------------------------------------------------
         *  check for elementary structures involving more than one
         *  closing pair.
         *  --------------------------------------------------------*/

        for (p = i + 1; p <= MIN2(j - 2 - turn, i + MAXLOOP + 1); p++) {
          int minq = j - i + p - MAXLOOP - 2;
          if (minq < p + 1 + turn)
            minq = p + 1 + turn;

          for (q = minq; q < j; q++) {
            type_2 = ptype[p][q - p];

            if (type_2 == 0)
              continue;

            type_2 = rtype[type_2];

            if (noGUclosure)
              if (no_close || (type_2 == 3) || (type_2 == 4))
                if ((p > i + 1) || (q < j - 1))
                  continue;

            /* continue unless stack */

            energy = E_IntLoop(p - i - 1,
                               j - q - 1,
                               type,
                               type_2,
                               S1[i + 1],
                               S1[j - 1],
                               S1[p - 1],
                               S1[q + 1],
                               P);
            new_c = MIN2(new_c, energy + c[p][q - p]);
            if ((p == i + 1) && (j == q + 1))
              stackEnergy = energy;                       /* remember stack energy */
          } /* end q-loop */
        } /* end p-loop */

        /* multi-loop decomposition ------------------------*/
        if (!no_close) {
          decomp  = DMLi1[j - 1 - (i + 1)];
          tt      = rtype[type];
          switch (dangle_model) {
            /* no dangle_model */
            case 0:
              decomp += E_MLstem(tt, -1, -1, P);
              break;
            /* double dangle_model */
            case 2:
              decomp += E_MLstem(tt, S1[j - 1], S1[i + 1], P);
              break;
            /* normal dangle_model, aka dangle_model = 1 */
            default:
              decomp  += E_MLstem(tt, -1, -1, P);
              decomp  = MIN2(decomp, DMLi2[j - 1 - (i + 2)] + E_MLstem(tt,
                                                                       -1,
                                                                       S1[i + 1],
                                                                       P) + P->MLbase);
              decomp = MIN2(decomp, DMLi2[j - 2 - (i + 2)] + E_MLstem(tt,
                                                                      S1[j - 1],
                                                                      S1[i + 1],
                                                                      P) + 2 * P->MLbase);
              decomp = MIN2(decomp, DMLi1[j - 2 - (i + 1)] + E_MLstem(tt,
                                                                      S1[j - 1],
                                                                      -1,
                                                                      P) + P->MLbase);
              break;
          }
          new_c = MIN2(new_c, decomp + P->MLclosing);
        }

        /* coaxial stacking of (i.j) with (i+1.k) or (k+1.j-1) */

        if (dangle_model == 3) {
          decomp = INF;
          for (k = i + 2 + turn; k < j - 2 - turn; k++) {
            type_2  = ptype[i + 1][k - i - 1];
            type_2  = rtype[type_2];
            if (type_2)
              decomp = MIN2(decomp, c[i + 1][k - i - 1] + P->stack[type][type_2] +
                            fML[k + 1][j - 1 - k - 1]);

            type_2  = ptype[k + 1][j - 1 - k - 1];
            type_2  = rtype[type_2];
            if (type_2)
              decomp = MIN2(decomp, c[k + 1][j - 1 - k - 1] + P->stack[type][type_2] +
                            fML[i + 1][k - i - 1]);
          }
          /* no TermAU penalty if coax stack */
          decomp  += 2 * P->MLintern[1] + P->MLclosing;
          new_c   = MIN2(new_c, decomp);
        }

        if (with_gquad) {
          /* include all cases where a g-quadruplex may be enclosed by base pair (i,j) */
          if (!no_close) {
            tt      = rtype[type];
            energy  = E_GQuad_IntLoop_L(i, j, type, S1, ggg, maxdist, P);
            new_c   = MIN2(new_c, energy);
          }
        }

        new_c     = MIN2(new_c, cc1[j - 1 - (i + 1)] + stackEnergy);
        cc[j - i] = new_c;
        if (noLP)
          c[i][j - i] = cc1[j - 1 - (i + 1)] + stackEnergy;
        else
          c[i][j - i] = cc[j - i];
      } /* end >> if (pair) << */
      else {
        c[i][j - i] = INF;
      }

      /* done with c[i,j], now compute fML[i,j] */
      /* free ends ? -----------------------------------------*/
      new_fML = INF;
      switch (dangle_model) {
        /* no dangle_model */
        case 0:
          new_fML = fML[i + 1][j - i - 1] + P->MLbase;
          new_fML = MIN2(new_fML, fML[i][j - 1 - i] + P->MLbase);
          new_fML = MIN2(new_fML, c[i][j - i] + E_MLstem(type, -1, -1, P));
          break;
        /* double dangle_model */
        case 2:
          new_fML = fML[i + 1][j - i - 1] + P->MLbase;
          new_fML = MIN2(fML[i][j - 1 - i] + P->MLbase, new_fML);
          new_fML =
            MIN2(new_fML,
                 c[i][j - i] +
                 E_MLstem(type, (i > 1) ? S1[i - 1] : -1, (j < length) ? S1[j + 1] : -1, P));
          break;
        /* normal dangle_model, aka dangle_model = 1 */
        default:  /* i unpaired */
          new_fML = fML[i + 1][j - i - 1] + P->MLbase;
          /* j unpaired */
          new_fML = MIN2(new_fML, fML[i][j - 1 - i] + P->MLbase);
          /* i,j */
          if (type)
            new_fML = MIN2(new_fML, c[i][j - i] + E_MLstem(type, -1, -1, P));

          /* i+1,j */
          tt = ptype[i + 1][j - i - 1];
          if (tt)
            new_fML = MIN2(new_fML, c[i + 1][j - i - 1] + E_MLstem(tt, S1[i], -1, P) + P->MLbase);

          /* i, j-1 */
          tt = ptype[i][j - 1 - i];
          if (tt)
            new_fML = MIN2(new_fML, c[i][j - 1 - i] + E_MLstem(tt, -1, S1[j], P) + P->MLbase);

          /* i+1,j-1 */
          tt = ptype[i + 1][j - 1 - i - 1];
          if (tt) {
            new_fML = MIN2(new_fML, c[i + 1][j - 1 - i - 1] + E_MLstem(tt,
                                                                       S1[i],
                                                                       S1[j],
                                                                       P) + 2 * P->MLbase);
          }

          break;
      }

      if (with_gquad)
        new_fML = MIN2(new_fML, ggg[i][j - i] + E_MLstem(0, -1, -1, P));

      /* modular decomposition -------------------------------*/
      for (decomp = INF, k = i + 1 + turn; k <= j - 2 - turn; k++)
        decomp = MIN2(decomp, Fmi[k - i] + fML[k + 1][j - k - 1]);

      DMLi[j - i] = decomp;               /* store for use in ML decompositon */
      new_fML     = MIN2(new_fML, decomp);

      /* coaxial stacking */
      if (dangle_model == 3) {
        /* additional ML decomposition as two coaxially stacked helices */
        for (decomp = INF, k = i + 1 + turn; k <= j - 2 - turn; k++) {
          type    = ptype[i][k - i];
          type    = rtype[type];
          type_2  = ptype[k + 1][j - k - 1];
          type_2  = rtype[type_2];
          if (type && type_2)
            decomp = MIN2(decomp,
                          c[i][k - i] + c[k + 1][j - k - 1] + P->stack[type][type_2]);
        }

        decomp += 2 * P->MLintern[1];          /* no TermAU penalty if coax stack */
#if 0
        /* This is needed for Y shaped ML loops with coax stacking of
        * interior pairts, but backtracking will fail if activated */
        DMLi[j - i] = MIN2(DMLi[j - i], decomp);
        DMLi[j - i] = MIN2(DMLi[j - i], DMLi[j - 1 - i] + P->MLbase);
        DMLi[j - i] = MIN2(DMLi[j - i], DMLi1[j - (i + 1)] + P->MLbase);
        new_fML     = MIN2(new_fML, DMLi[j - i]);
#endif
        new_fML = MIN2(new_fML, decomp);
      }

      fML[i][j - i] = Fmi[j - i] = new_fML;     /* substring energy */
    } /* for (j...) */

    /* calculate energies of 5' and 3' fragments */
    {
      char *ss = NULL;

      /* first case: i stays unpaired */
      f3[i] = f3[i + 1];

      /* next all cases where i is paired */
      switch (dangle_model) {
        /* dont use dangling end and mismatch contributions at all */
        case 0:
          for (j = i + turn + 1; j < length && j <= i + maxdist; j++) {
            type = ptype[i][j - i];

            if (with_gquad)
              f3[i] = MIN2(f3[i], f3[j + 1] + ggg[i][j - i]);

            if (type)
              f3[i] = MIN2(f3[i], f3[j + 1] + c[i][j - i] + E_ExtLoop(type, -1, -1, P));
          }
          if (length <= i + maxdist) {
            j = length;

            if (with_gquad)
              f3[i] = MIN2(f3[i], ggg[i][j - i]);

            type = ptype[i][j - i];
            if (type)
              f3[i] = MIN2(f3[i], c[i][j - i] + E_ExtLoop(type, -1, -1, P));
          }

          break;
        /* always use dangle_model on both sides */
        case 2:
          for (j = i + turn + 1; j < length && j <= i + maxdist; j++) {
            type = ptype[i][j - i];

            if (with_gquad)
              if (ggg[i][j - i] != INF)
                f3[i] = MIN2(f3[i], f3[j + 1] + ggg[i][j - i]);

            if (type) {
              f3[i] =
                MIN2(f3[i],
                     f3[j + 1] + c[i][j - i] +
                     E_ExtLoop(type, (i > 1) ? S1[i - 1] : -1, S1[j + 1], P));
            }
          }
          if (length <= i + maxdist) {
            j = length;

            if (with_gquad)
              f3[i] = MIN2(f3[i], ggg[i][j - i]);

            type = ptype[i][j - i];
            if (type)
              f3[i] = MIN2(f3[i], c[i][j - i] + E_ExtLoop(type, (i > 1) ? S1[i - 1] : -1, -1, P));
          }

          break;
        /* normal dangle_model, aka dangle_model = 1 */
        default:
          for (j = i + turn + 1; j < length && j <= i + maxdist; j++) {
            type = ptype[i][j - i];

            if (with_gquad)
              f3[i] = MIN2(f3[i], f3[j + 1] + ggg[i][j - i]);

            if (type) {
              f3[i] = MIN2(f3[i], f3[j + 1] + c[i][j - i] + E_ExtLoop(type, -1, -1, P));
              f3[i] =
                MIN2(f3[i],
                     ((j + 2 <=
                       length) ? f3[j + 2] : 0) + c[i][j - i] + E_ExtLoop(type, -1, S1[j + 1], P));
            }

            type = ptype[i + 1][j - i - 1];
            if (type) {
              f3[i] = MIN2(f3[i], f3[j + 1] + c[i + 1][j - i - 1] + E_ExtLoop(type, S1[i], -1, P));
              f3[i] = MIN2(f3[i],
                           ((j + 1 <
                             length) ? f3[j + 2] : 0) + c[i + 1][j - i - 1] +
                           E_ExtLoop(type, S1[i], S1[j + 1], P));
            }
          }
          if (length <= i + maxdist) {
            j = length;

            if (with_gquad)
              f3[i] = MIN2(f3[i], ggg[i][j - i]);

            type = ptype[i][j - i];
            if (type)
              f3[i] = MIN2(f3[i], c[i][j - i] + E_ExtLoop(type, -1, -1, P));

            type = ptype[i + 1][j - i - 1];
            if (type)
              f3[i] = MIN2(f3[i], c[i + 1][j - i - 1] + E_ExtLoop(type, S1[i], -1, P));
          }

          break;
      } /* switch(dangle_model)... */

      /* backtrack partial structure */
      if (f3[i] < f3[i + 1]) {
        do_backtrack = 1;
      } else if (do_backtrack) {
        int pairpartner; /*i+1?? is paired with pairpartner*/
        int cc;
        int traced2 = 0;
        fij   = f3[i + 1];
        lind  = i + 1;
        /*start "short" backtrack*/

        /*get paired base*/
        while (fij == f3[lind + 1])
          lind++;

        /*get pairpartner*/
        for (pairpartner = lind + turn; pairpartner <= lind + maxdist; pairpartner++) {
          type = ptype[lind][pairpartner - lind];
          switch (dangle_model) {
            case 0:
              if (type) {
                cc = c[lind][pairpartner - lind] + E_ExtLoop(type, -1, -1, P);
                if (fij == cc + f3[pairpartner + 1])
                  traced2 = 1;
              } else if (with_gquad) {
                cc = ggg[lind][pairpartner - lind];
                if (fij == cc + f3[pairpartner + 1])
                  traced2 = 1;
              }

              break;
            case 2:
              if (type) {
                cc = c[lind][pairpartner - lind] +
                     E_ExtLoop(type, (lind > 1) ? S1[lind - 1] : -1,
                               (pairpartner < length) ? S1[pairpartner + 1] : -1, P);
                if (fij == cc + f3[pairpartner + 1])
                  traced2 = 1;
              } else if (with_gquad) {
                cc = ggg[lind][pairpartner - lind];
                if (fij == cc + f3[pairpartner + 1])
                  traced2 = 1;
              }

              break;
            default:
              if (type) {
                cc = c[lind][pairpartner - lind] + E_ExtLoop(type, -1, -1, P);
                if (fij == cc + f3[pairpartner + 1]) {
                  traced2 = 1;
                  break;
                } else if (pairpartner < length) {
                  cc = c[lind][pairpartner - lind] + E_ExtLoop(type, -1, S1[pairpartner + 1], P);
                  if (fij == cc + f3[pairpartner + 2]) {
                    traced2 = 1;
                    break;
                  }
                }
              } else if (with_gquad) {
                cc = ggg[lind][pairpartner - lind];
                if (fij == cc + f3[pairpartner + 1])
                  traced2 = 1;
              }

              type = ptype[lind + 1][pairpartner - lind - 1];
              if (type) {
                cc = c[lind + 1][pairpartner - (lind + 1)] + E_ExtLoop(type, S1[lind], -1, P);
                if (fij == cc + f3[pairpartner + 1]) {
                  traced2 = 1;
                  break;
                } else if (pairpartner < length) {
                  cc = c[lind + 1][pairpartner - (lind + 1)] +
                       E_ExtLoop(type, S1[lind], S1[pairpartner + 1], P);
                  if (fij == cc + f3[pairpartner + 2])
                    traced2 = 1;
                }
              }

              break;
          }
          if (traced2)
            break;
        }
        if (!traced2)
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
          int     cc;
          double  average_free_energy;
          double  sd_free_energy;
          int     info_avg;
          double  my_z;
          int     traced2 = 0;
          fij   = f3[i];
          lind  = i;
          while (fij == f3[lind + 1])
            lind++;
          /*get pairpartner*/
          for (pairpartner = lind + turn; pairpartner <= lind + maxdist; pairpartner++) {
            type = ptype[lind][pairpartner - lind];
            switch (dangle_model) {
              case 0:
                if (type) {
                  cc = c[lind][pairpartner - lind] + E_ExtLoop(type, -1, -1, P);
                  if (fij == cc + f3[pairpartner + 1])
                    traced2 = 1;
                } else if (with_gquad) {
                  cc = ggg[lind][pairpartner - lind];
                  if (fij == cc + f3[pairpartner + 1])
                    traced2 = 1;
                }

                break;
              case 2:
                if (type) {
                  cc = c[lind][pairpartner - lind] + E_ExtLoop(type,
                                                               (lind > 1) ? S1[lind - 1] : -1,
                                                               (pairpartner <
                                                                length) ? S1[pairpartner + 1] : -1,
                                                               P);
                  if (fij == cc + f3[pairpartner + 1])
                    traced2 = 1;
                } else if (with_gquad) {
                  cc = ggg[lind][pairpartner - lind];
                  if (fij == cc + f3[pairpartner + 1])
                    traced2 = 1;
                }

                break;
              default:
                if (type) {
                  cc = c[lind][pairpartner - lind] + E_ExtLoop(type, -1, -1, P);
                  if (fij == cc + f3[pairpartner + 1]) {
                    traced2 = 1;
                    break;
                  } else if (pairpartner < length) {
                    cc = c[lind][pairpartner - lind] + E_ExtLoop(type, -1, S1[pairpartner + 1], P);
                    if (fij == cc + f3[pairpartner + 1]) {
                      traced2 = 1;
                      break;
                    }
                  }
                } else if (with_gquad) {
                  cc = ggg[lind][pairpartner - lind];
                  if (fij == cc + f3[pairpartner + 1])
                    traced2 = 1;
                }

                type = ptype[lind + 1][pairpartner - lind - 1];
                if (type) {
                  cc = c[lind + 1][pairpartner - (lind + 1)] + E_ExtLoop(type, S1[lind], -1, P);
                  if (fij == cc + f3[pairpartner + 1]) {
                    traced2 = 1;
                    break;
                  } else if (pairpartner < length) {
                    cc = c[lind + 1][pairpartner - (lind + 1)] +
                         E_ExtLoop(type, S1[lind], S1[pairpartner + 1], P);
                    if (fij == cc + f3[pairpartner + 2]) {
                      traced2 = 1;
                      break;
                    }
                  }
                }
            }
            if (traced2)
              break;
          }
          if (!traced2)
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
      int ii, *FF; /* rotate the auxilliary arrays */

      /* check for values close to integer underflow */
      if (INT_CLOSE_TO_UNDERFLOW(f3[i])) {
        /* correct f3 free energies and increase underflow counter */
        int cnt, cnt2;
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

      /*
       * rotate the DP matrices
       * NOTE: here we rotate them only locally, i.e. their
       * actual configuration within vc remains intact
       */
      if (i + maxdist + 4 <= length) {
        c[i - 1]                = c[i + maxdist + 4];
        c[i + maxdist + 4]      = NULL;
        fML[i - 1]              = fML[i + maxdist + 4];
        fML[i + maxdist + 4]    = NULL;
        ptype[i - 1]            = ptype[i + maxdist + 4];
        ptype[i + maxdist + 4]  = NULL;
        if (i > 1) {
          make_ptypes(vc, i - 1);
          if (with_gquad) {
            vrna_gquad_mx_local_update(vc, i - 1);
            ggg = vc->matrices->ggg_local;
          }
        }

        for (ii = 0; ii < maxdist + 5; ii++) {
          c[i - 1][ii]    = INF;
          fML[i - 1][ii]  = INF;
        }
      }
    }
  }

  free(cc);
  free(cc1);
  free(Fmi);
  free(DMLi);
  free(DMLi1);
  free(DMLi2);

  if (with_gquad) {
    for (i = 0; (i <= maxdist + 5) && (i <= length); i++)
      free(ggg[i]);
    free(ggg);
    vc->matrices->ggg_local = NULL;
  }

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
  sect          sector[MAXSECTORS];   /* backtracking sectors */
  int           i, j, k, length, energy, new, no_close, type, type_2, tt, s = 0;
  int           with_gquad, bt_type, turn, dangle_model, noLP, noGUclosure, *rtype;
  int           **c, **fML, *f3, **ggg;
  char          *string, *structure, **ptype;
  short         *S, *S1;
  vrna_param_t  *P;
  vrna_md_t     *md;

  string        = vc->sequence;
  length        = vc->length;
  S             = vc->sequence_encoding2;
  S1            = vc->sequence_encoding;
  ptype         = vc->ptype_local;
  P             = vc->params;
  md            = &(P->model_details);
  dangle_model  = md->dangles;
  noLP          = md->noLP;
  noGUclosure   = md->noGUclosure;
  with_gquad    = md->gquad;
  bt_type       = md->backtrack_type;
  turn          = md->min_loop_size;
  rtype         = &(md->rtype[0]);

  c   = vc->matrices->c_local;
  fML = vc->matrices->fML_local;
  f3  = vc->matrices->f3_local;
  ggg = vc->matrices->ggg_local;

  /* length = strlen(string); */
  sector[++s].i = start;
  sector[s].j   = MIN2(length, maxdist + 1);
  sector[s].ml  = (bt_type == 'M') ? 1 : ((bt_type == 'C') ? 2 : 0);

  structure = (char *)vrna_alloc((MIN2(length - start, maxdist) + 3) * sizeof(char));
  for (i = 0; i <= MIN2(length - start, maxdist); i++)
    structure[i] = '-';

  while (s > 0) {
    int ml, fij, cij, traced, i1, j1, mm, mm5, mm3, mm53, p, q, jj = 0, gq = 0;
    int canonical = 1;     /* (i,j) closes a canonical structure */
    i   = sector[s].i;
    j   = sector[s].j;
    ml  = sector[s--].ml;  /* ml is a flag indicating if backtracking is to
                           * occur in the fML- (1) or in the f-array (0) */
    if (ml == 2) {
      structure[i - start]  = '(';
      structure[j - start]  = ')';
      goto repeat1;
    }

    if (j < i + turn + 1)
      continue;                     /* no more pairs in this interval */

    fij = (ml) ? fML[i][j - i] : f3[i];

    if (ml == 0) {
      /* backtrack in f3 */

      if (fij == f3[i + 1]) {
        sector[++s].i = i + 1;
        sector[s].j   = j;
        sector[s].ml  = ml;
        continue;
      }

      /* i or i+1 is paired. Find pairing partner */
      switch (dangle_model) {
        case 0:
          for (traced = 0, k = j; k > i + turn; k--) {
            if (with_gquad) {
              if (fij == ggg[i][k - i] + f3[k + 1]) {
                /* found the decomposition */
                traced  = i;
                jj      = k + 1;
                gq      = 1;
                break;
              }
            }

            jj    = k + 1;
            type  = ptype[i][k - i];
            if (type) {
              if (fij == c[i][k - i] + E_ExtLoop(type, -1, -1, P) + f3[k + 1]) {
                traced = i;
                break;
              }
            }
          }
          break;
        case 2:
          for (traced = 0, k = j; k > i + turn; k--) {
            if (with_gquad) {
              if (fij == ggg[i][k - i] + f3[k + 1]) {
                /* found the decomposition */
                traced  = i;
                jj      = k + 1;
                gq      = 1;
                break;
              }
            }

            jj    = k + 1;
            type  = ptype[i][k - i];
            if (type) {
              if (fij ==
                  c[i][k - i] +
                  E_ExtLoop(type, (i > 1) ? S1[i - 1] : -1, (k < length) ? S1[k + 1] : -1,
                            P) + f3[k + 1]) {
                traced = i;
                break;
              }
            }
          }
          break;
        default:
          for (traced = 0, k = j; k > i + turn; k--) {
            if (with_gquad) {
              if (fij == ggg[i][k - i] + f3[k + 1]) {
                /* found the decomposition */
                traced  = i;
                jj      = k + 1;
                gq      = 1;
                break;
              }
            }

            jj    = k + 1;
            type  = ptype[i + 1][k - (i + 1)];
            if (type) {
              if (fij == c[i + 1][k - (i + 1)] + E_ExtLoop(type, S1[i], -1, P) + f3[k + 1])
                traced = i + 1;

              if (k < length) {
                if (fij ==
                    c[i + 1][k - (i + 1)] + E_ExtLoop(type, S1[i], S1[k + 1], P) + f3[k + 2]) {
                  traced  = i + 1;
                  jj      = k + 2;
                }
              }
            }

            type = ptype[i][k - i];
            if (type) {
              if (fij == c[i][k - i] + E_ExtLoop(type, -1, -1, P) + f3[k + 1])
                traced = i;

              if (k < length) {
                if (fij == c[i][k - i] + E_ExtLoop(type, -1, S1[k + 1], P) + f3[k + 2]) {
                  traced  = i;
                  jj      = k + 2;
                }
              }
            }

            if (traced)
              break;
          }
          break;
      } /* switch(dangle_model)...*/

      if (!traced)
        vrna_message_error("backtrack failed in f3");

      if (j == length) {
        /* backtrack only one component, unless j==length */
        sector[++s].i = jj;
        sector[s].j   = j;
        sector[s].ml  = ml;
      }

      i = traced;
      j = k;

      if (with_gquad && gq)
        /* goto backtrace of gquadruplex */
        goto repeat_gquad;

      structure[i - start]  = '(';
      structure[j - start]  = ')';
      if (((jj == j + 2) || (dangle_model == 2)) && (j < length))
        structure[j + 1 - start] = '.';

      goto repeat1;
    } else {
      /* trace back in fML array */
      if (fML[i][j - 1 - i] + P->MLbase == fij) {
        /* 3' end is unpaired */
        sector[++s].i = i;
        sector[s].j   = j - 1;
        sector[s].ml  = ml;
        continue;
      }

      if (fML[i + 1][j - (i + 1)] + P->MLbase == fij) {
        /* 5' end is unpaired */
        sector[++s].i = i + 1;
        sector[s].j   = j;
        sector[s].ml  = ml;
        continue;
      }

      if (with_gquad) {
        if (fij == ggg[i][j - i] + E_MLstem(0, -1, -1, P))
          /* go to backtracing of quadruplex */
          goto repeat_gquad;
      }

      switch (dangle_model) {
        case 0:
          tt = ptype[i][j - i];
          if (fij == c[i][j - i] + E_MLstem(tt, -1, -1, P)) {
            structure[i - start]  = '(';
            structure[j - start]  = ')';
            goto repeat1;
          }

          break;
        case 2:
          tt = ptype[i][j - i];
          if (fij ==
              c[i][j - i] +
              E_MLstem(tt, (i > 1) ? S1[i - 1] : -1, (j < length) ? S1[j + 1] : -1, P)) {
            structure[i - start]  = '(';
            structure[j - start]  = ')';
            goto repeat1;
          }

          break;
        default:
          tt = ptype[i][j - i];
          if (fij == c[i][j - i] + E_MLstem(tt, -1, -1, P)) {
            structure[i - start]  = '(';
            structure[j - start]  = ')';
            goto repeat1;
          }

          tt = ptype[i + 1][j - (i + 1)];
          if (fij == c[i + 1][j - (i + 1)] + E_MLstem(tt, S1[i], -1, P) + P->MLbase) {
            structure[++i - start]  = '(';
            structure[j - start]    = ')';
            goto repeat1;
          }

          tt = ptype[i][j - 1 - i];
          if (fij == c[i][j - 1 - i] + E_MLstem(tt, -1, S1[j], P) + P->MLbase) {
            structure[i - start]    = '(';
            structure[--j - start]  = ')';
            goto repeat1;
          }

          tt = ptype[i + 1][j - 1 - (i + 1)];
          if (fij == c[i + 1][j - 1 - (i + 1)] + E_MLstem(tt, S1[i], S1[j], P) + 2 * P->MLbase) {
            structure[++i - start]  = '(';
            structure[--j - start]  = ')';
            goto repeat1;
          }

          break;
      } /* switch(dangle_model)... */

      /* modular decomposition */
      for (k = i + 1 + turn; k <= j - 2 - turn; k++)
        if (fij == (fML[i][k - i] + fML[k + 1][j - (k + 1)]))
          break;

      if ((dangle_model == 3) && (k > j - 2 - turn)) {
        /* must be coax stack */
        ml = 2;
        for (k = i + 1 + turn; k <= j - 2 - turn; k++) {
          type    = ptype[i][k - i];
          type    = rtype[type];
          type_2  = ptype[k + 1][j - (k + 1)];
          type_2  = rtype[type_2];
          if (type && type_2)
            if (fij == c[i][k - i] + c[k + 1][j - (k + 1)] + P->stack[type][type_2] +
                2 * P->MLintern[1])
              break;
        }
      }

      sector[++s].i = i;
      sector[s].j   = k;
      sector[s].ml  = ml;
      sector[++s].i = k + 1;
      sector[s].j   = j;
      sector[s].ml  = ml;

      if (k > j - 2 - turn)
        vrna_message_error("backtrack failed in fML");

      continue;
    }

repeat1:

    /*----- begin of "repeat:" -----*/
    if (canonical)
      cij = c[i][j - i];

    type = ptype[i][j - i];


    if (noLP) {
      if (cij == c[i][j - i]) {
        /* (i.j) closes canonical structures, thus
         *  (i+1.j-1) must be a pair                */
        type_2                    = ptype[i + 1][j - 1 - (i + 1)];
        type_2                    = rtype[type_2];
        cij                       -= P->stack[type][type_2];
        structure[i + 1 - start]  = '(';
        structure[j - 1 - start]  = ')';
        i++;
        j--;
        canonical = 0;
        goto repeat1;
      }
    }

    canonical = 1;


    no_close = (((type == 3) || (type == 4)) && noGUclosure);
    if (no_close) {
      if (cij == FORBIDDEN)
        continue;
    } else
    if (cij == E_Hairpin(j - i - 1, type, S1[i + 1], S1[j - 1], string + i - 1, P)) {
      continue;
    }

    for (p = i + 1; p <= MIN2(j - 2 - turn, i + MAXLOOP + 1); p++) {
      int minq;
      minq = j - i + p - MAXLOOP - 2;
      if (minq < p + 1 + turn)
        minq = p + 1 + turn;

      for (q = j - 1; q >= minq; q--) {
        type_2 = ptype[p][q - p];
        if (type_2 == 0)
          continue;

        type_2 = rtype[type_2];
        if (noGUclosure)
          if (no_close || (type_2 == 3) || (type_2 == 4))
            if ((p > i + 1) || (q < j - 1))
              continue;

        /* continue unless stack */

        /* energy = oldLoopEnergy(i, j, p, q, type, type_2); */
        energy =
          E_IntLoop(p - i - 1, j - q - 1, type, type_2, S1[i + 1], S1[j - 1], S1[p - 1], S1[q + 1],
                    P);

        new     = energy + c[p][q - p];
        traced  = (cij == new);
        if (traced) {
          structure[p - start]  = '(';
          structure[q - start]  = ')';
          i                     = p, j = q;
          goto repeat1;
        }
      }
    }

    /* end of repeat: --------------------------------------------------*/

    /* (i.j) must close a multi-loop */
    tt  = rtype[type];
    i1  = i + 1;
    j1  = j - 1;

    if (with_gquad) {
      /*
       * The case that is handled here actually resembles something like
       * an interior loop where the enclosing base pair is of regular
       * kind and the enclosed pair is not a canonical one but a g-quadruplex
       * that should then be decomposed further...
       */
      if (backtrack_GQuad_IntLoop_L(cij, i, j, type, S, ggg, maxdist, &p, &q, P)) {
        i = p;
        j = q;
        goto repeat_gquad;
      }
    }

    sector[s + 1].ml = sector[s + 2].ml = 1;

    switch (dangle_model) {
      case 0:
        mm = P->MLclosing + E_MLstem(tt, -1, -1, P);
        for (k = i + 2 + turn; k < j - 2 - turn; k++)
          if (cij == fML[i + 1][k - (i + 1)] + fML[k + 1][j - 1 - (k + 1)] + mm)
            break;

        break;
      case 2:
        mm = P->MLclosing + E_MLstem(tt, S1[j - 1], S1[i + 1], P);
        for (k = i + 2 + turn; k < j - 2 - turn; k++)
          if (cij == fML[i + 1][k - (i + 1)] + fML[k + 1][j - 1 - (k + 1)] + mm)
            break;

        break;
      default:
        mm    = P->MLclosing + E_MLstem(tt, -1, -1, P);
        mm5   = P->MLclosing + E_MLstem(tt, S1[j - 1], -1, P) + P->MLbase;
        mm3   = P->MLclosing + E_MLstem(tt, -1, S1[i + 1], P) + P->MLbase;
        mm53  = P->MLclosing + E_MLstem(tt, S1[j - 1], S1[i + 1], P) + 2 * P->MLbase;
        for (k = i + 2 + turn; k < j - 2 - turn; k++) {
          if (cij == fML[i + 1][k - (i + 1)] + fML[k + 1][j - 1 - (k + 1)] + mm) {
            break;
          } else if (cij == fML[i + 2][k - (i + 2)] + fML[k + 1][j - 1 - (k + 1)] + mm3) {
            i1 = i + 2;
            break;
          } else if (cij == fML[i + 1][k - (i + 1)] + fML[k + 1][j - 2 - (k + 1)] + mm5) {
            j1 = j - 2;
            break;
          } else if (cij == fML[i + 2][k - (i + 2)] + fML[k + 1][j - 2 - (k + 1)] + mm53) {
            i1  = i + 2;
            j1  = j - 2;
            break;
          }

          /* coaxial stacking of (i.j) with (i+1.k) or (k.j-1) */
          /* use MLintern[1] since coax stacked pairs don't get TerminalAU */
          if (dangle_model == 3) {
            int en;
            type_2  = ptype[i + 1][k - (i + 1)];
            type_2  = rtype[type_2];
            if (type_2) {
              en = c[i + 1][k - (i + 1)] + P->stack[type][type_2] + fML[k + 1][j - 1 - (k + 1)];
              if (cij == en + 2 * P->MLintern[1] + P->MLclosing) {
                ml                = 2;
                sector[s + 1].ml  = 2;
                break;
              }
            }

            type_2  = ptype[k + 1][j - 1 - (k + 1)];
            type_2  = rtype[type_2];
            if (type_2) {
              en = c[k + 1][j - 1 - (k + 1)] + P->stack[type][type_2] + fML[i + 1][k - (i + 1)];
              if (cij == en + 2 * P->MLintern[1] + P->MLclosing) {
                sector[s + 2].ml = 2;
                break;
              }
            }
          }
        }
        break;
    } /* switch(dangle_model)... */

    if (k <= j - 3 - turn) {
      /* found the decomposition */
      sector[++s].i = i1;
      sector[s].j   = k;
      sector[++s].i = k + 1;
      sector[s].j   = j1;
    } else {
#if 0
      /* Y shaped ML loops fon't work yet */
      if (dangle_model == 3) {
        /* (i,j) must close a Y shaped ML loop with coax stacking */
        if (cij == fML[i + 1][j - 2 - (i + 2)] + mm + d3 + d5 + P->MLbase + P->MLbase) {
          i1  = i + 2;
          j1  = j - 2;
        } else if (cij == fML[i + 1][j - 2 - (i + 1)] + mm + d5 + P->MLbase) {
          j1 = j - 2;
        } else if (cij == fML[i + 2][j - 1 - (i + 2)] + mm + d3 + P->MLbase) {
          i1 = i + 2;
        } else /* last chance */
        if (cij != fML[i + 1][j - 1 - (i + 1)] + mm + P->MLbase) {
          fprintf(stderr, "backtracking failed in repeat");
        }

        /* if we arrive here we can express cij via fML[i1,j1]+dangle_model */
        sector[++s].i = i1;
        sector[s].j   = j1;
      } else
#endif
      vrna_message_error("backtracking failed in repeat");
    }

    continue; /* this is a workarround to not accidentally proceed in the following block */

repeat_gquad:
    /*
     * now we do some fancy stuff to backtrace the stacksize and linker lengths
     * of the g-quadruplex that should reside within position i,j
     */
    {
      int l[3], L, a;
      L = -1;

      get_gquad_pattern_mfe(S, i, j, P, &L, l);
      if (L != -1) {
        /* fill the G's of the quadruplex into the structure string */
        for (a = 0; a < L; a++) {
          structure[i + a - start]                                  = '+';
          structure[i + L + l[0] + a - start]                       = '+';
          structure[i + L + l[0] + L + l[1] + a - start]            = '+';
          structure[i + L + l[0] + L + l[1] + L + l[2] + a - start] = '+';
        }
        goto repeat_gquad_exit;
      }

      vrna_message_error("backtracking failed in repeat_gquad");
    }
repeat_gquad_exit:
    __asm("nop");
  }

  for (i = strlen(structure) - 1; i > 0 && structure[i] == '-'; i--)
    structure[i] = '\0';
  for (; i >= 0; i--)
    if (structure[i] == '-')
      structure[i] = '.';

  return structure;
}


PRIVATE int
fill_arrays_comparative(vrna_fold_compound_t      *fc,
                        vrna_mfe_window_callback  *cb,
                        void                      *data)
{
  /* fill "c", "fML" and "f3" arrays and return  optimal energy */

  char            **strings, **Ss, *prev;
  unsigned short  **a2s;
  short           **S, **S5, **S3, *S_cons;
  int             **pscore, u, i, j, k, length, energy, turn, *rtype, decomp, new_fML, MLenergy,
                  *type, type_2, tt, s, n_seq, lastf, lastf2, thisj, lastj, **c, **fML, *f3, **ggg,
                  *cc, *cc1, *Fmi, *DMLi, *DMLi1, *DMLi2, energyout, energyprev, maxdist,
                  dangle_model,
                  noLP, do_backtrack, prev_i, with_gquad, min_q, max_q, c0;
  float           **dm;
  vrna_mx_mfe_t   *matrices;
  vrna_param_t    *P;
  vrna_md_t       *md;

  int             olddm[7][7] = { { 0, 0, 0, 0, 0, 0, 0 }, /* hamming distance between pairs PRIVATE needed??*/
                                  { 0, 0, 2, 2, 1, 2, 2 } /* CG */,
                                  { 0, 2, 0, 1, 2, 2, 2 } /* GC */,
                                  { 0, 2, 1, 0, 2, 1, 2 } /* GU */,
                                  { 0, 1, 2, 2, 0, 2, 1 } /* UG */,
                                  { 0, 2, 2, 1, 2, 0, 2 } /* AU */,
                                  { 0, 2, 2, 2, 1, 2, 0 } /* UA */ };

  do_backtrack  = 0;
  prev_i        = 0;

  lastf = lastf2 = INF;
  dm    = NULL;
  prev  = NULL;

  strings       = fc->sequences;
  S             = fc->S;
  S_cons        = fc->S_cons;
  n_seq         = fc->n_seq;
  length        = fc->length;
  maxdist       = fc->window_size;
  matrices      = fc->matrices;
  P             = fc->params;
  md            = &(P->model_details);
  turn          = md->min_loop_size;
  dangle_model  = md->dangles;
  with_gquad    = md->gquad;
  noLP          = md->noLP;
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

  /* reserve additional memory for j-dimension */
  for (i = length; (i > length - maxdist - 5) && (i >= 0); i--) {
    fc->pscore_local[i]         = vrna_alloc(sizeof(int) * (maxdist + 5));
    fc->matrices->c_local[i]    = (int *)vrna_alloc(sizeof(int) * (maxdist + 5));
    fc->matrices->fML_local[i]  = (int *)vrna_alloc(sizeof(int) * (maxdist + 5));
  }

  for (i = length; (i >= length - maxdist - 4) && (i > 0); i--)
    make_pscores(fc, i, dm);

  pscore  = fc->pscore_local; /* precomputed array of pair types */
  S       = fc->S;
  S5      = fc->S5;           /*S5[s][i] holds next base 5' of i in sequence s*/
  S3      = fc->S3;           /*Sl[s][i] holds next base 3' of i in sequence s*/
  Ss      = fc->Ss;
  a2s     = fc->a2s;
  rtype   = &(md->rtype[0]);

  c   = fc->matrices->c_local;              /* energy array, given that i-j pair */
  fML = fc->matrices->fML_local;            /* multi-loop auxiliary energy array */
  f3  = fc->matrices->f3_local;             /* energy of 5' end */
  ggg = fc->matrices->ggg_local;

  cc    = (int *)vrna_alloc(sizeof(int) * (maxdist + 5));
  cc1   = (int *)vrna_alloc(sizeof(int) * (maxdist + 5));
  Fmi   = (int *)vrna_alloc(sizeof(int) * (maxdist + 5));
  DMLi  = (int *)vrna_alloc(sizeof(int) * (maxdist + 5));
  DMLi1 = (int *)vrna_alloc(sizeof(int) * (maxdist + 5));
  DMLi2 = (int *)vrna_alloc(sizeof(int) * (maxdist + 5));

  type = (int *)vrna_alloc(n_seq * sizeof(int));
  for (j = 0; j < maxdist + 5; j++)
    Fmi[j] = DMLi[j] = DMLi1[j] = DMLi2[j] = INF;
  for (j = length; j > length - maxdist - 3; j--)
    for (i = (length - maxdist - 2 > 0) ? length - maxdist - 2 : 1; i < j; i++)
      c[i][j - i] = fML[i][j - i] = INF;

  if (with_gquad) {
    vrna_gquad_mx_local_update(fc, length - maxdist - 4);
    ggg = fc->matrices->ggg_local;
  }

  for (i = length - turn - 1; i >= 1; i--) {
    /* i,j in [1..length] */
    for (j = i + 1; j <= length && j <= i + turn; j++)
      c[i][j - i] = fML[i][j - i] = INF;
    for (j = i + turn + 1; j <= length && j <= i + maxdist; j++) {
      int p, q, psc;
      /* bonus = 0;*/
      for (s = 0; s < n_seq; s++) {
        type[s] = md->pair[S[s][i]][S[s][j]];
        if (type[s] == 0)
          type[s] = 7;
      }

      psc = pscore[i][j - i];
      if (psc >= cv_fact * MINPSCORE) {
        /* we have a pair 2 consider */
        int new_c = 0, stackEnergy = INF;
        /* hairpin ----------------------------------------------*/
        for (new_c = s = 0; s < n_seq; s++) {
          if ((a2s[s][j - 1] - a2s[s][i]) < 3) {
            new_c += 600;
          } else {
            new_c +=
              E_Hairpin(a2s[s][j - 1] - a2s[s][i],
                        type[s],
                        S3[s][i],
                        S5[s][j],
                        Ss[s] + (a2s[s][i - 1]),
                        P);
          }
        }
        /*--------------------------------------------------------
         *  check for elementary structures involving more than one
         *  closing pair.
         *  --------------------------------------------------------*/
        for (p = i + 1; p <= MIN2(j - 2 - turn, i + MAXLOOP + 1); p++) {
          int minq = j - i + p - MAXLOOP - 2;
          if (minq < p + 1 + turn)
            minq = p + 1 + turn;

          for (q = minq; q < j; q++) {
            if (pscore[p][q - p] < MINPSCORE)
              continue;

            for (energy = s = 0; s < n_seq; s++) {
              type_2 = md->pair[S[s][q]][S[s][p]]; /* q,p not p,q! */
              if (type_2 == 0)
                type_2 = 7;

              energy += E_IntLoop(a2s[s][p - 1] - a2s[s][i],
                                  a2s[s][j - 1] - a2s[s][q],
                                  type[s],
                                  type_2,
                                  S3[s][i],
                                  S5[s][j],
                                  S5[s][p],
                                  S3[s][q],
                                  P);
            }
            new_c = MIN2(energy + c[p][q - p], new_c);
            if ((p == i + 1) && (j == q + 1))
              stackEnergy = energy;                       /* remember stack energy */
          } /* end q-loop */
        } /* end p-loop */


        /* multi-loop decomposition ------------------------*/
        decomp = DMLi1[j - 1 - (i + 1)];
        if (dangle_model) {
          for (s = 0; s < n_seq; s++) {
            tt      = rtype[type[s]];
            decomp  += E_MLstem(tt, S5[s][j], S3[s][i], P);
          }
        } else {
          for (s = 0; s < n_seq; s++) {
            tt      = rtype[type[s]];
            decomp  += E_MLstem(tt, -1, -1, P);
          }
        }

        MLenergy  = decomp + n_seq * P->MLclosing;
        new_c     = MIN2(new_c, MLenergy);

        if (with_gquad) {
          /* include all cases where a g-quadruplex may be enclosed by base pair (i,j) */
          energy = 0;
          for (s = 0; s < n_seq; s++) {
            tt = type[s];
            if (dangle_model == 2)
              energy += P->mismatchI[tt][S3[s][i]][S5[s][j]];

            if (tt > 2)
              energy += P->TerminalAU;
          }
          for (p = i + 2; p < j - VRNA_GQUAD_MIN_BOX_SIZE; p++) {
            u = p - i - 1;
            if (u > MAXLOOP)
              break;

            if (S_cons[p] != 3)
              continue;

            min_q = j - i + p - MAXLOOP - 2;
            c0    = p + VRNA_GQUAD_MIN_BOX_SIZE - 1;
            min_q = MAX2(c0, min_q);
            c0    = j - 1;
            max_q = p + VRNA_GQUAD_MAX_BOX_SIZE + 1;
            max_q = MIN2(c0, max_q);
            for (q = min_q; q < max_q; q++) {
              if (S_cons[q] != 3)
                continue;

              c0    = energy + ggg[p][q - p] + n_seq * P->internal_loop[u + j - q - 1];
              new_c = MIN2(new_c, c0);
            }
          }

          p = i + 1;
          if (S_cons[p] == 3) {
            if (p < j - VRNA_GQUAD_MIN_BOX_SIZE) {
              min_q = j - i + p - MAXLOOP - 2;
              c0    = p + VRNA_GQUAD_MIN_BOX_SIZE - 1;
              min_q = MAX2(c0, min_q);
              c0    = j - 3;
              max_q = p + VRNA_GQUAD_MAX_BOX_SIZE + 1;
              max_q = MIN2(c0, max_q);
              for (q = min_q; q < max_q; q++) {
                if (S_cons[q] != 3)
                  continue;

                c0    = energy + ggg[p][q - p] + n_seq * P->internal_loop[j - q - 1];
                new_c = MIN2(new_c, c0);
              }
            }
          }

          q = j - 1;
          if (S_cons[q] == 3) {
            for (p = i + 4; p < j - VRNA_GQUAD_MIN_BOX_SIZE; p++) {
              u = p - i - 1;
              if (u > MAXLOOP)
                break;

              if (S_cons[p] != 3)
                continue;

              c0    = energy + ggg[p][q - p] + n_seq * P->internal_loop[u];
              new_c = MIN2(new_c, c0);
            }
          }
        }

        new_c     = MIN2(new_c, cc1[j - 1 - (i + 1)] + stackEnergy);
        cc[j - i] = new_c - psc; /* add covariance bonnus/penalty */
        if (noLP)
          c[i][j - i] = cc1[j - 1 - (i + 1)] + stackEnergy - psc;
        else
          c[i][j - i] = cc[j - i];
      } /* end >> if (pair) << */
      else {
        c[i][j - i] = INF;
      }

      /* done with c[i,j], now compute fML[i,j] */
      /* free ends ? -----------------------------------------*/
      new_fML = INF;
      switch (dangle_model) {
        /* no dangling end contributions */
        case 0:
          new_fML = fML[i + 1][j - i - 1] + n_seq * P->MLbase;
          new_fML = MIN2(fML[i][j - 1 - i] + n_seq * P->MLbase, new_fML);
          energy  = c[i][j - i] /*+P->MLintern[type]*/;
          if (energy < INF)
            for (s = 0; s < n_seq; s++)
              energy += E_MLstem(type[s], -1, -1, P);

          new_fML = MIN2(energy, new_fML);
          break;

        /* double dangles */
        case 2:
          new_fML = fML[i + 1][j - i - 1] + n_seq * P->MLbase;
          new_fML = MIN2(fML[i][j - 1 - i] + n_seq * P->MLbase, new_fML);
          energy  = c[i][j - i] /*+P->MLintern[type]*/;
          if (energy < INF)
            for (s = 0; s < n_seq; s++)
              energy += E_MLstem(type[s], (i > 1) ? S5[s][i] : -1, (j < length) ? S3[s][j] : -1, P);

          new_fML = MIN2(energy, new_fML);
          break;

        default:
          vrna_message_warning("dangle model %d not implemented in comparative Lfold!",
                               dangle_model);
          break;
      }

      if (with_gquad)
        new_fML = MIN2(new_fML, ggg[i][j - i] + n_seq * E_MLstem(0, -1, -1, P));

      /* modular decomposition -------------------------------*/

      for (decomp = INF, k = i + 1 + turn; k <= j - 2 - turn; k++)
        decomp = MIN2(decomp, Fmi[k - i] + fML[k + 1][j - k - 1]);

      DMLi[j - i] = decomp;               /* store for use in ML decompositon */
      new_fML     = MIN2(new_fML, decomp);


      fML[i][j - i] = Fmi[j - i] = new_fML;     /* substring energy */
    } /* for (j...) */

    /* calculate energies of 5' and 3' fragments */
    {
      char  *ss;
      int   thisf = 0;

      /* first case: i stays unpaired */
      f3[i] = f3[i + 1];

      /* next all cases where i is paired */
      switch (dangle_model) {
        /* no dangling end contributions */
        case 0:
          for (j = i + turn + 1; j < length && j <= i + maxdist; j++) {
            if (f3[j + 1] < INF) {
              if (with_gquad) {
                energy = ggg[i][j - i];

                if (energy / (j - i + 1) < thisf) {
                  thisf = energy / (j - i + 1);
                  thisj = j;
                }

                energy  += f3[j + 1];
                f3[i]   = MIN2(f3[i], energy);
              }

              if (c[i][j - i] < INF) {
                energy = c[i][j - i];
                for (s = 0; s < n_seq; s++) {
                  tt = md->pair[S[s][i]][S[s][j]];
                  if (tt == 0)
                    tt = 7;

                  energy += E_ExtLoop(tt, -1, -1, P);
                }

                if (energy / (j - i + 1) < thisf) {
                  thisf = energy / (j - i + 1);
                  thisj = j;
                }

                energy  += f3[j + 1];
                f3[i]   = MIN2(f3[i], energy);
              }
            }
          }
          if (length <= i + maxdist) {
            j = length;

            if (with_gquad) {
              energy = ggg[i][j - i];

              if (energy / (j - i + 1) < thisf) {
                thisf = energy / (j - i + 1);
                thisj = j;
              }

              f3[i] = MIN2(f3[i], energy);
            }

            if (c[i][j - i] < INF) {
              energy = c[i][j - i];
              for (s = 0; s < n_seq; s++) {
                tt = md->pair[S[s][i]][S[s][j]];
                if (tt == 0)
                  tt = 7;

                energy += E_ExtLoop(tt, -1, -1, P);
              }

              /*  thisf=MIN2(energy/(j-i+1),thisf); ???*/
              if (energy / (j - i + 1) < thisf) {
                thisf = energy / (j - i + 1);
                thisj = j;
              }

              f3[i] = MIN2(f3[i], energy);
            }
          }

          break;

        /* double dangles */
        case 2:
          for (j = i + turn + 1; j < length && j <= i + maxdist; j++) {
            if (f3[j + 1] < INF) {
              if (with_gquad) {
                energy = ggg[i][j - i];

                if (energy / (j - i + 1) < thisf) {
                  thisf = energy / (j - i + 1);
                  thisj = j;
                }

                energy  += f3[j + 1];
                f3[i]   = MIN2(f3[i], energy);
              }

              if (c[i][j - i] < INF) {
                energy = c[i][j - i];
                for (s = 0; s < n_seq; s++) {
                  tt = md->pair[S[s][i]][S[s][j]];
                  if (tt == 0)
                    tt = 7;

                  energy += E_ExtLoop(tt, (i > 1) ? S5[s][i] : -1, S3[s][j], P);
                }

                if (energy / (j - i + 1) < thisf) {
                  thisf = energy / (j - i + 1);
                  thisj = j;
                }

                energy  += f3[j + 1];
                f3[i]   = MIN2(f3[i], energy);
              }
            }
          }
          if (length <= i + maxdist) {
            j = length;

            if (with_gquad) {
              energy = ggg[i][j - i];

              if (energy / (j - i + 1) < thisf) {
                thisf = energy / (j - i + 1);
                thisj = j;
              }

              f3[i] = MIN2(f3[i], energy);
            }

            if (c[i][j - i] < INF) {
              energy = c[i][j - i];
              for (s = 0; s < n_seq; s++) {
                tt = md->pair[S[s][i]][S[s][j]];
                if (tt == 0)
                  tt = 7;

                energy += E_ExtLoop(tt, (i > 1) ? S5[s][i] : -1, -1, P);
              }

              /*  thisf=MIN2(energy/(j-i+1),thisf); ???*/
              if (energy / (j - i + 1) < thisf) {
                thisf = energy / (j - i + 1);
                thisj = j;
              }

              f3[i] = MIN2(f3[i], energy);
            }
          }

          break;

        default:
          vrna_message_warning("dangle model %d not implemented in comparative Lfold!",
                               dangle_model);
          break;
      }


      /* backtrack partial structure */
      /* if (i+maxdist<length) {*/
      if (i < length - 1) {
        if (f3[i] != f3[i + 1]) {
          do_backtrack    = 1;
          backtrack_type  = 'F';
          if (prev_i == 0) {
            free(prev);
            prev          = backtrack_comparative(fc, i, MIN2(maxdist, length - i));
            prev_i        = i;
            do_backtrack  = 0;
            lastf2        = lastf;
            energyprev    = f3[i];
          }
        } else if ((thisf < lastf) && (thisf < lastf2) && ((thisf / (n_seq * 100)) < -0.01)) {
          /*?????????*/
          do_backtrack    = 2;
          backtrack_type  = 'C';
        } else if (do_backtrack) {
          if (do_backtrack == 1) {
            ss        = backtrack_comparative(fc, i + 1, MIN2(maxdist, length - i) /*+1*/);
            energyout = f3[i] - f3[i + strlen(ss) - 1];/*??*/
          } else {
            ss        = backtrack_comparative(fc, i + 1, lastj - i - 2);
            energyout = c[i + 1][lastj - (i + 1)];
            switch (dangle_model) {
              case 0:
                for (s = 0; s < n_seq; s++) {
                  int type;
                  type = md->pair[S[s][i + 1]][S[s][lastj - i]];
                  if (type == 0)
                    type = 7;

                  energyout += E_ExtLoop(type, -1, -1, P);
                }
                break;

              case 2:
                for (s = 0; s < n_seq; s++) {
                  int type;
                  type = md->pair[S[s][i + 1]][S[s][lastj - i]];
                  if (type == 0)
                    type = 7;

                  energyout += E_ExtLoop(type, (i > 1) ? S5[s][i + 1] : -1, S3[s][lastj - i], P);
                }
                break;

              default:
                vrna_message_warning("dangle model %d not implemented in comparative Lfold!",
                                     dangle_model);
                break;
            }
          }

          if ((prev) && ((prev_i + strlen(prev) > i + 1 + strlen(ss)) || (do_backtrack == 2))) {
            char *outstr = (char *)vrna_alloc(sizeof(char) * (strlen(prev) + 1));
            strncpy(outstr, strings[0] + prev_i - 1, strlen(prev));
            outstr[strlen(prev)] = '\0';

            /* execute callback */
            cb(prev_i, prev_i + (int)strlen(prev) - 1, prev, energyprev / (100. * n_seq), data);

            free(outstr);
          }

          free(prev);
          prev            = ss;
          energyprev      = energyout;
          prev_i          = i + 1;
          do_backtrack    = 0;
          backtrack_type  = 'F';
        }
      }

      lastf2  = lastf;
      lastf   = thisf;
      lastj   = thisj;


      if (i == 1) {
        char *outstr = NULL;
        if (prev) {
          outstr = (char *)vrna_alloc(sizeof(char) * (strlen(prev) + 1));
          strncpy(outstr, strings[0] + prev_i - 1, strlen(prev));
          outstr[strlen(prev)] = '\0';

          /* execute callback */
          cb(prev_i, prev_i + (int)strlen(prev) - 1, prev, (energyprev) / (100. * n_seq), data);
        }

        if ((f3[prev_i] != f3[1]) || !prev) {
          ss = backtrack_comparative(fc, i, maxdist);
          free(outstr);

          outstr = (char *)vrna_alloc(sizeof(char) * (strlen(ss) + 1));
          strncpy(outstr, strings[0], strlen(ss));
          outstr[strlen(ss)] = '\0';

          /* execute callback */
          cb(1, (int)strlen(ss) - 1, ss, (f3[1] - f3[1 + strlen(ss) - 1]) / (100. * n_seq), data);

          free(ss);
          ss = NULL;
        }

        free(outstr);
        free(prev);
        prev = NULL;
      }
    }
    {
      int ii, *FF; /* rotate the auxilliary arrays */
      FF    = DMLi2;
      DMLi2 = DMLi1;
      DMLi1 = DMLi;
      DMLi  = FF;
      FF    = cc1;
      cc1   = cc;
      cc    = FF;
      for (j = 0; j < maxdist + 5; j++)
        cc[j] = Fmi[j] = DMLi[j] = INF;
      if (i <= length - maxdist - 4) {
        matrices->c_local[i - 1]              = matrices->c_local[i + maxdist + 4];
        matrices->c_local[i + maxdist + 4]    = NULL;
        matrices->fML_local[i - 1]            = matrices->fML_local[i + maxdist + 4];
        matrices->fML_local[i + maxdist + 4]  = NULL;
        fc->pscore_local[i - 1]               = fc->pscore_local[i + maxdist + 4];
        fc->pscore_local[i + maxdist + 4]     = NULL;
        c                                     = matrices->c_local;
        fML                                   = matrices->fML_local;

        if (i > 1) {
          make_pscores(fc, i - 1, dm);
          if (with_gquad) {
            vrna_gquad_mx_local_update(fc, i - 1);
            ggg = fc->matrices->ggg_local;
          }
        }

        pscore = fc->pscore_local;  /* precomputed array of pair types */

        for (ii = 0; ii < maxdist + 5; ii++)
          c[i - 1][ii] = fML[i - 1][ii] = INF;
      }
    }
  }

  free(type);
  free(cc);
  free(cc1);
  free(Fmi);
  free(DMLi);
  free(DMLi1);
  free(DMLi2);

  for (i = 0; i < 7; i++)
    free(dm[i]);
  free(dm);

  /* free additional memory for j-dimension */
  for (i = 0; (i < maxdist + 5) && (i <= length); i++) {
    free(fc->pscore_local[i]);
    free(fc->matrices->c_local[i]);
    free(fc->matrices->fML_local[i]);
  }

  if (with_gquad) {
    for (i = 0; (i <= maxdist + 5) && (i <= length); i++)
      free(ggg[i]);
    free(ggg);
    fc->matrices->ggg_local = NULL;
  }

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

  char            *structure, **Ss;
  unsigned short  **a2s;
  short           **S, **S5, **S3, *S_cons;
  int             **pscore, i, j, k, energy, length, turn, dangle_model, noLP, *type, type_2, tt,
                  n_seq,
                  with_gquad, *rtype, s, ss, **c, **fML, *f3, **ggg, l1, minq, maxq, c0;
  vrna_param_t    *P;
  vrna_md_t       *md;

  n_seq         = fc->n_seq;
  length        = fc->length;
  S             = fc->S;
  S_cons        = fc->S_cons;
  S5            = fc->S5; /*S5[s][i] holds next base 5' of i in sequence s*/
  S3            = fc->S3; /*Sl[s][i] holds next base 3' of i in sequence s*/
  Ss            = fc->Ss;
  a2s           = fc->a2s;
  pscore        = fc->pscore_local; /* precomputed array of pair types */
  P             = fc->params;
  md            = &(P->model_details);
  rtype         = &(md->rtype[0]);
  with_gquad    = md->gquad;
  dangle_model  = md->dangles;
  noLP          = md->noLP;

  c   = fc->matrices->c_local;    /* energy array, given that i-j pair */
  fML = fc->matrices->fML_local;  /* multi-loop auxiliary energy array */
  f3  = fc->matrices->f3_local;   /* energy of 5' end */
  ggg = fc->matrices->ggg_local;

  type          = (int *)vrna_alloc(n_seq * sizeof(int));
  s             = 0;
  turn          = md->min_loop_size;
  sector[++s].i = start;
  sector[s].j   = MIN2(length, start + maxdist + 1);
  sector[s].ml  = (backtrack_type == 'M') ? 1 : ((backtrack_type == 'C') ? 2 : 0);

  structure = (char *)vrna_alloc((MIN2(length - start, maxdist) + 3) * sizeof(char));
  for (i = 0; i <= MIN2(length - start, maxdist); i++)
    structure[i] = '.';

  while (s > 0) {
    int ml, fij, cij, traced, i1, j1, mm, p, q, jj = 0, gq = 0;
    int canonical = 1;     /* (i,j) closes a canonical structure */
    i   = sector[s].i;
    j   = sector[s].j;
    ml  = sector[s--].ml;  /* ml is a flag indicating if backtracking is to
                           * occur in the fML- (1) or in the f-array (0) */
    if (ml == 2) {
      structure[i - start]  = '(';
      structure[j - start]  = ')';
      goto repeat1_comparative;
    }

    if (j < i + turn + 1)
      continue;                 /* no more pairs in this interval */

    fij = (ml) ? fML[i][j - i] : f3[i];

    if (ml == 0) {
      /* backtrack in f3 */

      if (fij == f3[i + 1]) {
        sector[++s].i = i + 1;
        sector[s].j   = j;
        sector[s].ml  = ml;
        continue;
      }

      /* i is paired. Find pairing partner */
      switch (dangle_model) {
        case 0:
          for (k = i + turn + 1, traced = 0; k <= j; k++) {
            int cc;
            if (with_gquad) {
              if (fij == ggg[i][k - i] + f3[k + 1]) {
                traced  = i;
                jj      = k + 1;
                gq      = 1;
                break;
              }
            }

            jj  = k + 1;
            cc  = c[i][k - (i)];
            if (cc < INF) {
              for (ss = 0; ss < n_seq; ss++) {
                type[ss] = md->pair[S[ss][i]][S[ss][k]];
                if (type[ss] == 0)
                  type[ss] = 7;

                cc += E_ExtLoop(type[ss], -1, -1, P);
              }

              if (fij == cc + f3[k + 1])
                traced = i;
            }

            if (traced)
              break;
          }

          break;
        case 2:
          for (k = i + turn + 1, traced = 0; k <= j; k++) {
            int cc;
            if (with_gquad) {
              if (fij == ggg[i][k - i] + f3[k + 1]) {
                /* found the decomposition */
                traced  = i;
                jj      = k + 1;
                gq      = 1;
                break;
              }
            }

            jj  = k + 1;
            cc  = c[i][k - (i)];
            if (cc < INF) {
              for (ss = 0; ss < n_seq; ss++) {
                type[ss] = md->pair[S[ss][i]][S[ss][k]];
                if (type[ss] == 0)
                  type[ss] = 7;

                cc +=
                  E_ExtLoop(type[ss], (i > 1) ? S5[ss][i] : -1, (k < length) ? S3[ss][k] : -1, P);
              }

              if (fij == cc + f3[k + 1])
                traced = i;
            }

            if (traced)
              break;
          }
          break;
      }


      if (!traced)
        vrna_message_error("backtrack failed in f3");

      if (j == length) {
        /* backtrack only one component, unless j==length */
        sector[++s].i = jj;
        sector[s].j   = j;
        sector[s].ml  = ml;
      }

      i = traced;
      j = k;

      if (with_gquad && gq)
        /* goto backtrace of gquadruplex */
        goto repeat_gquad_comparative;

      structure[i - start]  = '(';
      structure[j - start]  = ')';
      goto repeat1_comparative;
    } else {
      /* trace back in fML array */
      if (fML[i][j - 1 - i] + n_seq * P->MLbase == fij) {
        /* 3' end is unpaired */
        sector[++s].i = i;
        sector[s].j   = j - 1;
        sector[s].ml  = ml;
        continue;
      }

      if (fML[i + 1][j - (i + 1)] + n_seq * P->MLbase == fij) {
        /* 5' end is unpaired */
        sector[++s].i = i + 1;
        sector[s].j   = j;
        sector[s].ml  = ml;
        continue;
      }

      if (with_gquad) {
        if (fij == ggg[i][j - i] + n_seq * E_MLstem(0, -1, -1, P))
          /* go to backtracing of quadruplex */
          goto repeat_gquad_comparative;
      }

      cij = c[i][j - i];
      if (dangle_model) {
        for (ss = 0; ss < n_seq; ss++) {
          tt = md->pair[S[ss][i]][S[ss][j]];
          if (tt == 0)
            tt = 7;

          cij += E_MLstem(tt, (i > 1) ? S5[ss][i] : -1, (j < length) ? S3[ss][j] : -1, P);
        }
      } else {
        for (ss = 0; ss < n_seq; ss++) {
          tt = md->pair[S[ss][i]][S[ss][j]];
          if (tt == 0)
            tt = 7;

          cij += E_MLstem(tt, -1, -1, P);
        }
      }

      if (fij == cij) {
        /* found a pair */
        structure[i - start]  = '(';
        structure[j - start]  = ')';
        goto repeat1_comparative;
      }

      for (k = i + 1 + turn; k <= j - 2 - turn; k++)
        if (fij == (fML[i][k - i] + fML[k + 1][j - (k + 1)]))
          break;

      sector[++s].i = i;
      sector[s].j   = k;
      sector[s].ml  = ml;
      sector[++s].i = k + 1;
      sector[s].j   = j;
      sector[s].ml  = ml;

      if (k > j - 2 - turn)
        vrna_message_error("backtrack failed in fML");

      continue;
    }

repeat1_comparative:

    /*----- begin of "repeat:" -----*/
    if (canonical)
      cij = c[i][j - i];

    for (ss = 0; ss < n_seq; ss++) {
      type[ss] = md->pair[S[ss][i]][S[ss][j]];
      if (type[ss] == 0)
        type[ss] = 7;
    }

    /* check for case when fill_arrays_comparative wants to backtrack in 'C' but only saw a G-quadruplex without enclosing base pairs */
    if ((cij == INF) && (with_gquad) && (ggg[i][j - i] < INF))
      goto repeat_gquad_comparative;

    if (noLP) {
      if (cij == c[i][j - i]) {
        /* (i.j) closes canonical structures, thus
         *  (i+1.j-1) must be a pair                */
        for (ss = 0; ss < n_seq; ss++) {
          type_2 = md->pair[S[ss][j - 1]][S[ss][i + 1]];  /* j,i not i,j */
          if (type_2 == 0)
            type_2 = 7;

          cij -= P->stack[type[ss]][type_2];
        }
        cij                       += pscore[i][j - i];
        structure[i + 1 - start]  = '(';
        structure[j - 1 - start]  = ')';
        i++;
        j--;
        canonical = 0;
        goto repeat1_comparative;
      }
    }

    canonical = 1;
    cij       += pscore[i][j - i];
    {
      int cc = 0;
      for (ss = 0; ss < n_seq; ss++) {
        if ((a2s[ss][j - 1] - a2s[ss][i]) < 3)
          cc += 600;
        else
          cc +=
            E_Hairpin(a2s[ss][j - 1] - a2s[ss][i], type[ss], S3[ss][i], S5[ss][j],
                      Ss[ss] + a2s[ss][i - 1], P);
      }
      if (cij == cc) /* found hairpin */
        continue;
    }

    for (p = i + 1; p <= MIN2(j - 2 - turn, i + MAXLOOP + 1); p++) {
      int minq;
      minq = j - i + p - MAXLOOP - 2;
      if (minq < p + 1 + turn)
        minq = p + 1 + turn;

      for (q = j - 1; q >= minq; q--) {
        if (c[p][q - p] >= INF)
          continue;

        for (ss = energy = 0; ss < n_seq; ss++) {
          type_2 = md->pair[S[ss][q]][S[ss][p]];  /* q,p not p,q */
          if (type_2 == 0)
            type_2 = 7;

          energy += E_IntLoop(a2s[ss][p - 1] - a2s[ss][i],
                              a2s[ss][j - 1] - a2s[ss][q],
                              type[ss],
                              type_2,
                              S3[ss][i],
                              S5[ss][j],
                              S5[ss][p],
                              S3[ss][q],
                              P);
        }
        traced = (cij == energy + c[p][q - p]);
        if (traced) {
          structure[p - start]  = '(';
          structure[q - start]  = ')';
          i                     = p, j = q;
          goto repeat1_comparative;
        }
      }
    }

    /* end of repeat: --------------------------------------------------*/
    i1  = i + 1;
    j1  = j - 1;

    if (with_gquad) {
      /*
       * The case that is handled here actually resembles something like
       * an interior loop where the enclosing base pair is of regular
       * kind and the enclosed pair is not a canonical one but a g-quadruplex
       * that should then be decomposed further...
       */
      mm = 0;
      for (ss = 0; ss < n_seq; ss++) {
        tt = type[ss];

        if (dangle_model == 2)
          mm += P->mismatchI[tt][S3[ss][i]][S5[ss][j]];

        if (tt > 2)
          mm += P->TerminalAU;
      }

      for (p = i + 2;
           p < j - VRNA_GQUAD_MIN_BOX_SIZE;
           p++) {
        if (S_cons[p] != 3)
          continue;

        l1 = p - i - 1;
        if (l1 > MAXLOOP)
          break;

        minq  = j - i + p - MAXLOOP - 2;
        c0    = p + VRNA_GQUAD_MIN_BOX_SIZE - 1;
        minq  = MAX2(c0, minq);
        c0    = j - 1;
        maxq  = p + VRNA_GQUAD_MAX_BOX_SIZE + 1;
        maxq  = MIN2(c0, maxq);
        for (q = minq; q < maxq; q++) {
          if (S_cons[q] != 3)
            continue;

          c0 = mm + ggg[p][q - p] + n_seq * P->internal_loop[l1 + j - q - 1];
          if (cij == c0) {
            i = p;
            j = q;
            goto repeat_gquad_comparative;
          }
        }
      }
      p = i1;
      if (S_cons[p] == 3) {
        if (p < j - VRNA_GQUAD_MIN_BOX_SIZE) {
          minq  = j - i + p - MAXLOOP - 2;
          c0    = p + VRNA_GQUAD_MIN_BOX_SIZE - 1;
          minq  = MAX2(c0, minq);
          c0    = j - 3;
          maxq  = p + VRNA_GQUAD_MAX_BOX_SIZE + 1;
          maxq  = MIN2(c0, maxq);
          for (q = minq; q < maxq; q++) {
            if (S_cons[q] != 3)
              continue;

            if (cij == mm + ggg[p][q - p] + n_seq * P->internal_loop[j - q - 1]) {
              i = p;
              j = q;
              goto repeat_gquad_comparative;
            }
          }
        }
      }

      q = j1;
      if (S_cons[q] == 3) {
        for (p = i1 + 3; p < j - VRNA_GQUAD_MIN_BOX_SIZE; p++) {
          l1 = p - i - 1;
          if (l1 > MAXLOOP)
            break;

          if (S_cons[p] != 3)
            continue;

          if (cij == mm + ggg[p][q - p] + n_seq * P->internal_loop[l1]) {
            i = p;
            j = q;
            goto repeat_gquad_comparative;
          }
        }
      }
    }

    /* (i.j) must close a multi-loop */
    mm = n_seq * P->MLclosing;
    if (dangle_model) {
      for (ss = 0; ss < n_seq; ss++) {
        tt  = rtype[type[ss]];
        mm  += E_MLstem(tt, S5[ss][j], S3[ss][i], P);
      }
    } else {
      for (ss = 0; ss < n_seq; ss++) {
        tt  = rtype[type[ss]];
        mm  += E_MLstem(tt, -1, -1, P);
      }
    }

    sector[s + 1].ml = sector[s + 2].ml = 1;

    for (k = i + turn + 2; k < j - turn - 2; k++)
      if (cij == fML[i + 1][k - (i + 1)] + fML[k + 1][j - 1 - (k + 1)] + mm)
        break;

    if (k <= j - 3 - turn) {
      /* found the decomposition */
      sector[++s].i = i1;
      sector[s].j   = k;
      sector[++s].i = k + 1;
      sector[s].j   = j1;
    } else {
      vrna_message_error("backtracking failed in repeat");
    }

    continue; /* this is a workarround to not accidentally proceed in the following block */

repeat_gquad_comparative:
    /*
     * now we do some fancy stuff to backtrace the stacksize and linker lengths
     * of the g-quadruplex that should reside within position i,j
     */
    {
      int cnt1, l[3], L, size;
      size = j - i + 1;

      for (L = 0; L < VRNA_GQUAD_MIN_STACK_SIZE; L++) {
        if (S_cons[i + L] != 3)
          break;

        if (S_cons[j - L] != 3)
          break;
      }

      if (L == VRNA_GQUAD_MIN_STACK_SIZE) {
        /* continue only if minimum stack size starting from i is possible */
        for (; L <= VRNA_GQUAD_MAX_STACK_SIZE; L++) {
          if (S_cons[i + L - 1] != 3)
            break;                      /* break if no more consecutive G's 5' */

          if (S_cons[j - L + 1] != 3)
            break;                      /* break if no more consecutive G'1 3' */

          for (l[0] = VRNA_GQUAD_MIN_LINKER_LENGTH;
               (l[0] <= VRNA_GQUAD_MAX_LINKER_LENGTH)
               && (size - 4 * L - 2 * VRNA_GQUAD_MIN_LINKER_LENGTH - l[0] >= 0);
               l[0]++) {
            /* check whether we find the second stretch of consecutive G's */
            for (cnt1 = 0; (cnt1 < L) && (S_cons[i + L + l[0] + cnt1] == 3); cnt1++);
            if (cnt1 < L)
              continue;

            for (l[1] = VRNA_GQUAD_MIN_LINKER_LENGTH;
                 (l[1] <= VRNA_GQUAD_MAX_LINKER_LENGTH)
                 && (size - 4 * L - VRNA_GQUAD_MIN_LINKER_LENGTH - l[0] - l[1] >= 0);
                 l[1]++) {
              /* check whether we find the third stretch of consectutive G's */
              for (cnt1 = 0; (cnt1 < L) && (S_cons[i + 2 * L + l[0] + l[1] + cnt1] == 3); cnt1++);
              if (cnt1 < L)
                continue;

              /*
               * the length of the third linker now depends on position j as well
               * as the other linker lengths... so we do not have to loop too much
               */
              l[2] = size - 4 * L - l[0] - l[1];
              if (l[2] < VRNA_GQUAD_MIN_LINKER_LENGTH)
                break;

              if (l[2] > VRNA_GQUAD_MAX_LINKER_LENGTH)
                continue;

              /* check for contribution */
              if (ggg[i][j - i] == E_gquad_ali(i, L, l, (const short **)S, n_seq, P)) {
                int a;
                /* fill the G's of the quadruplex into base pair stack */
                for (a = 0; a < L; a++) {
                  structure[i + a - start]                                  = '+';
                  structure[i + L + l[0] + a - start]                       = '+';
                  structure[i + L + l[0] + L + l[1] + a - start]            = '+';
                  structure[i + L + l[0] + L + l[1] + L + l[2] + a - start] = '+';
                }
                goto repeat_gquad_comparative_exit;
              }
            }
          }
        }
      }

      vrna_message_error("backtracking failed in repeat_gquad_comparative");
    }
repeat_gquad_comparative_exit:
    __asm("nop");
  }
  if (start + maxdist < length)
    for (i = strlen(structure); i > 0 && structure[i - 1] == '.'; i--)
      structure[i] = '\0';

  free(type);

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
