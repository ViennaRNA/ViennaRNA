/* Last changed Time-stamp: <2007-09-04 11:16:01 ivo> */
/*
                  minimum free energy
                  RNA secondary structure prediction
                  with maximum distance base pairs

                  c Ivo Hofacker, Peter Stadler

                  Vienna RNA package
*/

#include <config.h>
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
#include "ViennaRNA/Lfold.h"

#ifdef USE_SVM
#include "svm.h"
#include "svm_utils.h"
#endif

#define MAXSECTORS                  500   /* dimension for a backtrack array */
#define INT_CLOSE_TO_UNDERFLOW(i)   ((i) <= (INT_MIN/16))
#define UNDERFLOW_CORRECTION        (INT_MIN/32)

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
PRIVATE float wrap_Lfold( vrna_fold_compound_t *vc,
                          int with_zsc,
                          double min_z,
                          FILE *file);
PRIVATE void  make_ptypes(vrna_fold_compound_t *vc,
                          int i);
PRIVATE char  *backtrack( vrna_fold_compound_t *vc,
                          int start,
                          int maxdist);
PRIVATE int   fill_arrays(vrna_fold_compound_t *vc,
                          int with_zsc,
                          double min_z,
#ifdef USE_SVM
                          struct svm_model *avg_model,
                          struct svm_model *sd_model,
#endif
                          int *underflow,
                          FILE *output);

/*
#################################
# BEGIN OF FUNCTION DEFINITIONS #
#################################
*/

PUBLIC float
vrna_mfe_window( vrna_fold_compound_t *vc,
            FILE *file){

  return wrap_Lfold(vc, 0, 0.0, file);
}

#ifdef USE_SVM

PUBLIC float
vrna_mfe_window_zscore( vrna_fold_compound_t *vc,
             double min_z,
             FILE *file){

  return wrap_Lfold(vc, 1, min_z, file);
}

#endif

/*
#####################################
# BEGIN OF STATIC HELPER FUNCTIONS  #
#####################################
*/

PRIVATE float
wrap_Lfold( vrna_fold_compound_t *vc,
            int with_zsc,
            double min_z,
            FILE *file){

  int     i, energy, underflow, n, maxdist;
  float   mfe_local;
  FILE    *out;

#ifdef USE_SVM
  struct svm_model  *avg_model = NULL;
  struct svm_model  *sd_model = NULL;
#endif

  n       = vc->length;
  maxdist = vc->window_size;
  out     = (file) ? file : stdout;

  for(i = n; (i >= (int)n - (int)maxdist - 4) && (i > 0); i--)
    make_ptypes(vc, i);

#ifdef USE_SVM  /*svm*/
  if(with_zsc){
    avg_model = svm_load_model_string(avg_model_string);
    sd_model  = svm_load_model_string(sd_model_string);
  }
#endif

  /* keep track of how many times we were close to an integer underflow */
  underflow = 0;

#ifdef USE_SVM
  energy = fill_arrays(vc, with_zsc, min_z, avg_model, sd_model, &underflow, out);
  if(with_zsc){
    svm_free_model_content(avg_model);
    svm_free_model_content(sd_model);
  }
#else
  energy = fill_arrays(vc, with_zsc, min_z, &underflow, out);
#endif

  mfe_local = (underflow > 0) ? ((float)underflow * (float)(UNDERFLOW_CORRECTION)) / 100. : 0.;
  mfe_local += (float)energy/100.;

  return mfe_local;
}

PRIVATE int
fill_arrays(vrna_fold_compound_t *vc,
            int zsc,
            double min_z,
#ifdef USE_SVM
            struct svm_model *avg_model,
            struct svm_model *sd_model,
#endif
            int *underflow,
            FILE *output){

  /* fill "c", "fML" and "f3" arrays and return  optimal energy */

  int   i, j, k, length, energy, maxdist;
  int   **c, **fML, *f3, **ggg;
  int   decomp, new_fML;
  int   no_close, type, type_2, tt, with_gquad, dangle_model, noLP, noGUclosure, turn;
  int   *rtype;
  int   fij;
  int   lind;

  int   *cc = NULL;        /* linear array for calculating canonical structures */
  int   *cc1 = NULL;       /*   "     "        */
  int   *Fmi = NULL;       /* holds row i of fML (avoids jumps in memory) */
  int   *DMLi = NULL;      /* DMLi[j] holds MIN(fML[i,k]+fML[k+1,j])  */
  int   *DMLi1 = NULL;     /*             MIN(fML[i+1,k]+fML[k+1,j])  */
  int   *DMLi2 = NULL;     /*             MIN(fML[i+2,k]+fML[k+1,j])  */

  short         *S, *S1;
  char          *string, **ptype, *prev;
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

  prev          = NULL;

  c           = vc->matrices->c_local;
  fML         = vc->matrices->fML_local;
  f3          = vc->matrices->f3_local;
  ggg         = vc->matrices->ggg_local;

  cc    = (int *)   vrna_alloc(sizeof(int)   *(maxdist+5));
  cc1   = (int *)   vrna_alloc(sizeof(int)   *(maxdist+5));
  Fmi   = (int *)   vrna_alloc(sizeof(int)   *(maxdist+5));
  DMLi  = (int *)   vrna_alloc(sizeof(int)   *(maxdist+5));
  DMLi1 = (int *)   vrna_alloc(sizeof(int)   *(maxdist+5));
  DMLi2 = (int *)   vrna_alloc(sizeof(int)   *(maxdist+5));


  for (j=0; j<maxdist+5; j++)
    Fmi[j]=DMLi[j]=DMLi1[j]=DMLi2[j]=INF;
  for (j=length; j>length-maxdist-4; j--) {
    for (i=(length-maxdist-4>0)?length-maxdist-4:1 ; i<j; i++)
      c[i][j-i] = fML[i][j-i] = INF;
  }

  if(with_gquad){
    vrna_gquad_mx_local_update(vc, length - maxdist - 4);
  }

  for (i = length-turn-1; i >= 1; i--) { /* i,j in [1..length] */
    for (j = i+turn+1; j <= length && j <= i+maxdist; j++) {
      int p, q;
      type = ptype[i][j-i];

      no_close = (((type==3)||(type==4))&&noGUclosure);

      if (type) {   /* we have a pair */
        int new_c=0, stackEnergy=INF;
        /* hairpin ----------------------------------------------*/

        new_c = (no_close) ? FORBIDDEN : E_Hairpin(j-i-1, type, S1[i+1], S1[j-1], string+i-1, P);

        /*--------------------------------------------------------
          check for elementary structures involving more than one
          closing pair.
          --------------------------------------------------------*/

        for (p = i+1; p <= MIN2(j-2-turn,i+MAXLOOP+1) ; p++){
          int minq = j-i+p-MAXLOOP-2;
          if (minq<p+1+turn) minq = p+1+turn;
          for (q = minq; q < j; q++) {
            type_2 = ptype[p][q-p];

            if (type_2==0) continue;
            type_2 = rtype[type_2];

            if (noGUclosure)
              if (no_close||(type_2==3)||(type_2==4))
                if ((p>i+1)||(q<j-1)) continue;  /* continue unless stack */

            energy = E_IntLoop(p-i-1, j-q-1, type, type_2, S1[i+1], S1[j-1], S1[p-1], S1[q+1],P);
            new_c = MIN2(new_c, energy + c[p][q-p]);
            if ((p==i+1)&&(j==q+1)) stackEnergy = energy; /* remember stack energy */
          } /* end q-loop */
        } /* end p-loop */

        /* multi-loop decomposition ------------------------*/
        if (!no_close) {
          decomp  = DMLi1[j-1-(i+1)];
          tt      = rtype[type];
          switch(dangle_model){
            /* no dangle_model */
            case 0:   decomp += E_MLstem(tt, -1, -1, P);
                      break;
            /* double dangle_model */
            case 2:   decomp += E_MLstem(tt, S1[j-1], S1[i+1], P);
                      break;
            /* normal dangle_model, aka dangle_model = 1 */
            default:  decomp += E_MLstem(tt, -1, -1, P);
                      decomp = MIN2(decomp, DMLi2[j-1-(i+2)] + E_MLstem(tt, -1, S1[i+1], P) + P->MLbase);
                      decomp = MIN2(decomp, DMLi2[j-2-(i+2)] + E_MLstem(tt, S1[j-1], S1[i+1], P) + 2*P->MLbase);
                      decomp = MIN2(decomp, DMLi1[j-2-(i+1)] + E_MLstem(tt, S1[j-1], -1, P) + P->MLbase);
                      break;
          }
          new_c = MIN2(new_c, decomp + P->MLclosing);
        }

        /* coaxial stacking of (i.j) with (i+1.k) or (k+1.j-1) */

        if (dangle_model==3) {
          decomp = INF;
          for (k = i+2+turn; k < j-2-turn; k++) {
            type_2 = ptype[i+1][k-i-1]; type_2 = rtype[type_2];
            if (type_2)
              decomp = MIN2(decomp, c[i+1][k-i-1]+P->stack[type][type_2]+
                            fML[k+1][j-1-k-1]);
            type_2 = ptype[k+1][j-1-k-1]; type_2 = rtype[type_2];
            if (type_2)
              decomp = MIN2(decomp, c[k+1][j-1-k-1]+P->stack[type][type_2]+
                            fML[i+1][k-i-1]);
          }
          /* no TermAU penalty if coax stack */
          decomp += 2*P->MLintern[1] + P->MLclosing;
          new_c = MIN2(new_c, decomp);
        }

        if(with_gquad){
          /* include all cases where a g-quadruplex may be enclosed by base pair (i,j) */
          if (!no_close) {
            tt = rtype[type];
            energy = E_GQuad_IntLoop_L(i, j, type, S1, ggg, maxdist, P);
            new_c = MIN2(new_c, energy);
          }
        }

        new_c = MIN2(new_c, cc1[j-1-(i+1)]+stackEnergy);
        cc[j-i] = new_c;
        if (noLP)
          c[i][j-i] = cc1[j-1-(i+1)]+stackEnergy;
        else
          c[i][j-i] = cc[j-i];

      } /* end >> if (pair) << */

      else c[i][j-i] = INF;

      /* done with c[i,j], now compute fML[i,j] */
      /* free ends ? -----------------------------------------*/
      new_fML = INF;
      switch(dangle_model){
        /* no dangle_model */
        case 0:   new_fML = fML[i+1][j-i-1] + P->MLbase;
                  new_fML = MIN2(new_fML, fML[i][j-1-i] + P->MLbase);
                  new_fML = MIN2(new_fML, c[i][j-i] + E_MLstem(type, -1, -1, P));
                  break;
        /* double dangle_model */
        case 2:   new_fML = fML[i+1][j-i-1] + P->MLbase;
                  new_fML = MIN2(fML[i][j-1-i] + P->MLbase, new_fML);
                  new_fML = MIN2(new_fML,  c[i][j-i] + E_MLstem(type, (i>1) ? S1[i-1] : -1, (j<length) ? S1[j+1] : -1, P));
                  break;
        /* normal dangle_model, aka dangle_model = 1 */
        default:  /* i unpaired */
                  new_fML = fML[i+1][j-i-1] + P->MLbase;
                  /* j unpaired */
                  new_fML = MIN2(new_fML, fML[i][j-1-i] + P->MLbase);
                  /* i,j */
                  if(type) new_fML = MIN2(new_fML, c[i][j-i] + E_MLstem(type, -1, -1, P));
                  /* i+1,j */
                  tt = ptype[i+1][j-i-1];
                  if(tt) new_fML = MIN2(new_fML, c[i+1][j-i-1] + E_MLstem(tt, S1[i], -1, P) + P->MLbase);
                  /* i, j-1 */
                  tt = ptype[i][j-1-i];
                  if(tt) new_fML = MIN2(new_fML, c[i][j-1-i] + E_MLstem(tt, -1, S1[j], P) + P->MLbase);
                  /* i+1,j-1 */
                  tt = ptype[i+1][j-1-i-1];
                  if(tt) new_fML = MIN2(new_fML, c[i+1][j-1-i-1] + E_MLstem(tt, S1[i], S1[j], P) + 2*P->MLbase);
                  break;
      }

      if(with_gquad){
        new_fML = MIN2(new_fML, ggg[i][j - i] + E_MLstem(0, -1, -1, P));
      }

      /* modular decomposition -------------------------------*/
      for (decomp = INF, k = i+1+turn; k <= j-2-turn; k++)
        decomp = MIN2(decomp, Fmi[k-i]+fML[k+1][j-k-1]);

      DMLi[j-i] = decomp;               /* store for use in ML decompositon */
      new_fML   = MIN2(new_fML, decomp);

      /* coaxial stacking */
      if (dangle_model==3) {
        /* additional ML decomposition as two coaxially stacked helices */
        for (decomp = INF, k = i+1+turn; k <= j-2-turn; k++) {
          type = ptype[i][k-i]; type = rtype[type];
          type_2 = ptype[k+1][j-k-1]; type_2 = rtype[type_2];
          if (type && type_2)
            decomp = MIN2(decomp,
                          c[i][k-i]+c[k+1][j-k-1]+P->stack[type][type_2]);
        }

        decomp += 2*P->MLintern[1];          /* no TermAU penalty if coax stack */
#if 0
        /* This is needed for Y shaped ML loops with coax stacking of
           interior pairts, but backtracking will fail if activated */
        DMLi[j-i] = MIN2(DMLi[j-i], decomp);
        DMLi[j-i] = MIN2(DMLi[j-i], DMLi[j-1-i]+P->MLbase);
        DMLi[j-i] = MIN2(DMLi[j-i], DMLi1[j-(i+1)]+P->MLbase);
        new_fML = MIN2(new_fML, DMLi[j-i]);
#endif
        new_fML = MIN2(new_fML, decomp);
      }
      fML[i][j-i] = Fmi[j-i] = new_fML;     /* substring energy */
    } /* for (j...) */

    /* calculate energies of 5' and 3' fragments */
    {
      static int do_backtrack = 0, prev_i=0;
      char *ss=NULL;
      double prevz;

      /* first case: i stays unpaired */
      f3[i] = f3[i+1];

      /* next all cases where i is paired */
      switch(dangle_model){
        /* dont use dangling end and mismatch contributions at all */
        case 0:   for(j=i+turn+1; j<length && j<=i+maxdist; j++){
                    type = ptype[i][j-i];

                    if(with_gquad){
                      f3[i] = MIN2(f3[i], f3[j+1] + ggg[i][j-i]);
                    }

                    if(type)
                      f3[i] = MIN2(f3[i], f3[j+1] + c[i][j-i] + E_ExtLoop(type, -1, -1, P));
                  }
                  if(length<=i+maxdist){
                    j=length;

                    if(with_gquad){
                      f3[i] = MIN2(f3[i], ggg[i][j-i]);
                    }

                    type = ptype[i][j-i];
                    if(type)
                      f3[i] = MIN2(f3[i], c[i][j-i] + E_ExtLoop(type, -1, -1, P));
                  }
                  break;
        /* always use dangle_model on both sides */
        case 2:   for(j=i+turn+1; j<length && j<=i+maxdist; j++){
                    type = ptype[i][j-i];

                    if(with_gquad){
                      if(ggg[i][j-i] != INF)
                        f3[i] = MIN2(f3[i], f3[j+1] + ggg[i][j-i]);
                    }

                    if(type)
                      f3[i] = MIN2(f3[i], f3[j+1] + c[i][j-i] + E_ExtLoop(type, (i>1) ? S1[i-1] : -1, S1[j+1], P));
                  }
                  if(length<=i+maxdist){
                    j=length;

                    if(with_gquad){
                      f3[i] = MIN2(f3[i], ggg[i][j-i]);
                    }

                    type = ptype[i][j-i];
                    if(type)
                      f3[i] = MIN2(f3[i], c[i][j-i] + E_ExtLoop(type, (i>1) ? S1[i-1] : -1, -1, P));
                  }
                  break;
        /* normal dangle_model, aka dangle_model = 1 */
        default:  for(j=i+turn+1; j<length && j<=i+maxdist; j++){
                    type = ptype[i][j-i];

                    if(with_gquad){
                      f3[i] = MIN2(f3[i], f3[j+1] + ggg[i][j-i]);
                    }

                    if(type){
                      f3[i] = MIN2(f3[i], f3[j+1] + c[i][j-i] + E_ExtLoop(type, -1, -1, P));
                      f3[i] = MIN2(f3[i], ((j+2<=length) ? f3[j+2] : 0) + c[i][j-i] + E_ExtLoop(type, -1, S1[j+1], P));
                    }
                    type = ptype[i+1][j-i-1];
                    if(type){
                      f3[i] = MIN2(f3[i], f3[j+1] + c[i+1][j-i-1] + E_ExtLoop(type, S1[i], -1, P));
                      f3[i] = MIN2(f3[i], ((j + 1 < length) ? f3[j+2] : 0) + c[i+1][j-i-1] + E_ExtLoop(type, S1[i], S1[j+1], P));
                    }
                  }
                  if(length<=i+maxdist){
                    j     = length;

                    if(with_gquad){
                      f3[i] = MIN2(f3[i], ggg[i][j-i]);
                    }

                    type  = ptype[i][j-i];
                    if(type)
                      f3[i] = MIN2(f3[i], c[i][j-i] + E_ExtLoop(type, -1, -1, P));
                    type  = ptype[i+1][j-i-1];
                    if(type)
                      f3[i] = MIN2(f3[i], c[i+1][j-i-1] + E_ExtLoop(type, S1[i], -1, P));
                  }
                  break;
      } /* switch(dangle_model)... */

      /* backtrack partial structure */
      if (f3[i] < f3[i+1]){
        do_backtrack=1;
      }
      else if (do_backtrack) {
        int pairpartner; /*i+1?? is paired with pairpartner*/
        int cc;
        int traced2=0;
        fij = f3[i+1];
        lind=i+1;
        /*start "short" backtrack*/

        /*get paired base*/
        while(fij==f3[lind+1])
          lind++;

        /*get pairpartner*/
        for (pairpartner = lind + turn; pairpartner <= lind + maxdist; pairpartner++){
          type = ptype[lind][pairpartner-lind];
          switch(dangle_model){
            case 0:   if(type){
                        cc = c[lind][pairpartner-lind] + E_ExtLoop(type, -1, -1, P);
                        if(fij == cc + f3[pairpartner + 1])
                          traced2 = 1;
                      }
                      else if(with_gquad) {
                        cc = ggg[lind][pairpartner-lind];
                        if(fij == cc + f3[pairpartner + 1])
                          traced2 = 1;
                      }

                      break;
            case 2:   if(type){
                        cc = c[lind][pairpartner-lind] + E_ExtLoop(type, (lind > 1) ? S1[lind-1] : -1, (pairpartner < length) ? S1[pairpartner+1] : -1, P);
                        if(fij == cc + f3[pairpartner + 1])
                          traced2 = 1;
                      }
                      else if(with_gquad){
                        cc = ggg[lind][pairpartner-lind];
                        if(fij == cc + f3[pairpartner + 1])
                          traced2 = 1;
                      }

                      break;
            default:  if(type){
                        cc = c[lind][pairpartner-lind] + E_ExtLoop(type, -1, -1, P);
                        if(fij == cc + f3[pairpartner + 1]){
                          traced2 = 1;
                          break;
                        }
                        else if(pairpartner < length){
                          cc = c[lind][pairpartner-lind] + E_ExtLoop(type, -1, S1[pairpartner+1], P);
                          if(fij == cc + f3[pairpartner + 2]){
                            traced2 = 1;
                            break;
                          }
                        }
                      }
                      else if(with_gquad){
                        cc = ggg[lind][pairpartner-lind];
                        if(fij == cc + f3[pairpartner + 1])
                          traced2 = 1;
                      }

                      type = ptype[lind+1][pairpartner-lind-1];
                      if(type){
                        cc = c[lind+1][pairpartner-(lind+1)] + E_ExtLoop(type, S1[lind], -1, P);
                        if(fij == cc + f3[pairpartner+1]){
                          traced2 = 1;
                          break;
                        }
                        else if(pairpartner < length){
                          cc = c[lind+1][pairpartner-(lind+1)] + E_ExtLoop(type, S1[lind], S1[pairpartner+1], P);
                          if(fij == cc + f3[pairpartner+2])
                            traced2 = 1;
                        }
                      }
                      break;
          }
          if(traced2) break;
        }
        if (!traced2) vrna_message_error("backtrack failed in short backtrack 1");
        if (zsc){
#ifdef USE_SVM
          int info_avg;
          double average_free_energy;
          double sd_free_energy;
          double my_z;
          int *AUGC = get_seq_composition(S, lind-1, MIN2((pairpartner+1),length), length);
          /*\svm*/
          average_free_energy = avg_regression(AUGC[0], AUGC[1], AUGC[2], AUGC[3], AUGC[4], avg_model, &info_avg);
          if (info_avg == 0)  {
            double difference;
            double min_sd = minimal_sd(AUGC[0],AUGC[1],AUGC[2],AUGC[3],AUGC[4]);
            difference=(fij-f3[pairpartner+1])/100.-average_free_energy;
            if ( difference - ( min_z * min_sd ) <= 0.0001 ) {
              sd_free_energy = sd_regression(AUGC[0],AUGC[1],AUGC[2],AUGC[3],AUGC[4],sd_model);
              my_z=difference/sd_free_energy;
              if (my_z<=min_z){
                ss =  backtrack(vc, lind, pairpartner+1);
                if (prev) {
                  if ((i+strlen(ss)<prev_i+strlen(prev)) ||
                      strncmp(ss+prev_i-i,prev,strlen(prev))) { /* ss does not contain prev */
                    if (dangle_model==2)
                      fprintf(output, ".%s (%6.2f) %4d z= %.3f\n", prev, (f3[prev_i]-f3[prev_i+strlen(prev)-1])/100., prev_i-1, prevz);
                    else
                      fprintf(output, "%s (%6.2f) %4d z=%.3f\n ", prev, (f3[prev_i]-f3[prev_i+strlen(prev)])/100., prev_i, prevz);
                  }
                  free(prev);
                }
                prev=ss; prev_i = lind; prevz=my_z;
              }
            }

          }
          free(AUGC);
          do_backtrack=0;
#endif
        }
        else {
          /* original code for Lfold*/
          ss =  backtrack(vc, lind , pairpartner+1);
          if (prev) {
            if ((i+strlen(ss)<prev_i+strlen(prev)) || strncmp(ss+prev_i-i,prev,strlen(prev))){
              /* ss does not contain prev */
              if (dangle_model==2){
                fprintf(output, ".%s (%6.2f) %4d\n", prev, (f3[prev_i]-f3[prev_i+strlen(prev)-1])/100., prev_i-1);
              } else
                fprintf(output, "%s (%6.2f) %4d\n", prev, (f3[prev_i]-f3[prev_i+strlen(prev)])/100., prev_i);
            }
            free(prev);
          }
          prev=ss;
          prev_i = lind;
          do_backtrack=0;
        }
      }
      if (i==1) {
        if (prev) {
          if(zsc) {
            if (dangle_model==2)
              fprintf(output, ".%s (%6.2f) %4d z= %.2f\n", prev, (f3[prev_i]-f3[prev_i+strlen(prev)-1])/100., prev_i-1, prevz);
           else
              fprintf(output, "%s (%6.2f) %4dz= %.2f \n", prev, (f3[prev_i]-f3[prev_i+strlen(prev)])/100., prev_i, prevz);
          }
          else {
            if (dangle_model==2)
              fprintf(output, ".%s (%6.2f) %4d\n", prev, (f3[prev_i]-f3[prev_i+strlen(prev)-1])/100., prev_i-1);
            else
              fprintf(output, "%s (%6.2f) %4d\n", prev, (f3[prev_i]-f3[prev_i+strlen(prev)])/100., prev_i);
          }
          free(prev); prev=NULL;
        } else if ((f3[i]<0) && (!zsc)) do_backtrack=1;

        if (do_backtrack) {
          int pairpartner; /*i+1?? is paired with pairpartner*/
          int cc;
          double average_free_energy;
          double sd_free_energy;
          int info_avg;
          double my_z;
          int traced2 = 0;
          fij = f3[i];
          lind=i;
          while(fij==f3[lind+1]) lind++;
          /*get pairpartner*/
          for(pairpartner = lind + turn; pairpartner <= lind + maxdist; pairpartner++){
            type = ptype[lind][pairpartner-lind];
            switch(dangle_model){
              case 0:   if(type){
                          cc = c[lind][pairpartner-lind] + E_ExtLoop(type, -1, -1, P);
                          if(fij == cc + f3[pairpartner + 1])
                            traced2 = 1;
                        }
                        else if(with_gquad){
                          cc = ggg[lind][pairpartner-lind];
                          if(fij == cc + f3[pairpartner + 1])
                            traced2 = 1;
                        }

                        break;
              case 2:   if(type){
                          cc = c[lind][pairpartner-lind] + E_ExtLoop(type, (lind > 1) ? S1[lind-1] : -1, (pairpartner < length) ? S1[pairpartner+1] : -1, P);
                          if(fij == cc + f3[pairpartner + 1])
                            traced2 = 1;
                        }
                        else if(with_gquad){
                          cc = ggg[lind][pairpartner-lind];
                          if(fij == cc + f3[pairpartner + 1])
                            traced2 = 1;
                        }

                        break;
              default:  if(type){
                          cc = c[lind][pairpartner-lind] + E_ExtLoop(type, -1, -1, P);
                          if(fij == cc + f3[pairpartner + 1]){
                            traced2 = 1;
                            break;
                          }
                          else if(pairpartner < length){
                            cc = c[lind][pairpartner-lind] + E_ExtLoop(type, -1, S1[pairpartner + 1], P);
                            if(fij == cc + f3[pairpartner + 1]){
                              traced2 = 1;
                              break;
                            }
                          }
                        }
                        else if(with_gquad){
                          cc = ggg[lind][pairpartner-lind];
                          if(fij == cc + f3[pairpartner + 1])
                            traced2 = 1;
                        }

                        type = ptype[lind+1][pairpartner-lind-1];
                        if(type){
                          cc = c[lind+1][pairpartner-(lind+1)] + E_ExtLoop(type, S1[lind], -1, P);
                          if(fij == cc + f3[pairpartner+1]){
                            traced2 = 1;
                            break;
                          }
                          else if (pairpartner < length){
                            cc = c[lind+1][pairpartner-(lind+1)] + E_ExtLoop(type, S1[lind], S1[pairpartner+1], P);
                            if(fij == cc + f3[pairpartner + 2]){
                              traced2 =1;
                              break;
                            }
                          }
                        }
            }
            if(traced2) break;
          }
          if (!traced2) vrna_message_error("backtrack failed in short backtrack 2");

          if(zsc){
#ifdef USE_SVM
            int *AUGC = get_seq_composition(S, lind-1, MIN2((pairpartner+1),length), length);
            average_free_energy = avg_regression(AUGC[0],AUGC[1],AUGC[2],AUGC[3],AUGC[4],avg_model,&info_avg);
            if (info_avg == 0)  {
              double difference;
              double min_sd = minimal_sd(AUGC[0],AUGC[1],AUGC[2],AUGC[3],AUGC[4]);
              difference=(fij-f3[pairpartner+1])/100.-average_free_energy;
              if ( difference - ( min_z * min_sd ) <= 0.0001 ) {
                sd_free_energy = sd_regression(AUGC[0],AUGC[1],AUGC[2],AUGC[3],AUGC[4],sd_model);
                my_z=difference/sd_free_energy;
                if (my_z<=min_z){
                  ss =  backtrack(vc, lind , pairpartner+1);
                  fprintf(output, "%s (%6.2f) %4d z= %.2f\n", ss, (f3[lind]-f3[lind+strlen(ss)-1])/100., lind, my_z);
                }
              }
            }
            free(AUGC);
#endif
          }
          else {
            ss =  backtrack(vc, lind , pairpartner+1);
            if (dangle_model==2)
              fprintf(output, "%s (%6.2f) %4d\n", ss, (f3[lind]-f3[lind+strlen(ss)-1])/100., 1);
            else
              fprintf(output, "%s (%6.2f) %4d\n", ss, (f3[lind]-f3[lind+strlen(ss)])/100., 1);
            free(ss);
          }
        }
        do_backtrack=0;
      }
    }
    {
      int ii, *FF; /* rotate the auxilliary arrays */

      /* check for values close to integer underflow */
      if(INT_CLOSE_TO_UNDERFLOW(f3[i])){
        /* correct f3 free energies and increase underflow counter */
        int cnt, cnt2;
        for(cnt=i; cnt <= length && cnt <= lind + maxdist + 2; cnt++) {
          f3[cnt] -= UNDERFLOW_CORRECTION;
        }
        (*underflow)++;
      }

      FF = DMLi2; DMLi2 = DMLi1; DMLi1 = DMLi; DMLi = FF;
      FF = cc1; cc1=cc; cc=FF;
      for(j = 0; j < maxdist + 5; j++){
        cc[j] = Fmi[j] = DMLi[j] = INF;
      }

      /*
        rotate the DP matrices
        NOTE: here we rotate them only locally, i.e. their
        actual configuration within vc remains intact
      */
      if( i + maxdist + 4 <= length ){
        c[i - 1]                = c[i + maxdist + 4];
        c[i + maxdist + 4]      = NULL;
        fML[i - 1]              = fML[i + maxdist + 4];
        fML[i + maxdist + 4]    = NULL;
        ptype[i - 1]            = ptype[i + maxdist + 4];
        ptype[i + maxdist + 4]  = NULL;
        if( i > 1 ){
          make_ptypes(vc, i - 1);
          if(with_gquad){
            vrna_gquad_mx_local_update(vc, i - 1);
          }
        }
        for(ii = 0; ii < maxdist + 5; ii++){
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

  return f3[1];
}

PRIVATE char *
backtrack(vrna_fold_compound_t *vc,
          int start,
          int maxdist){

  /*------------------------------------------------------------------
    trace back through the "c", "f3" and "fML" arrays to get the
    base pairing list. No search for equivalent structures is done.
    This is fast, since only few structure elements are recalculated.
    ------------------------------------------------------------------*/
  sect          sector[MAXSECTORS];   /* backtracking sectors */
  int           i, j, k, length, energy, new, no_close, type, type_2, tt, s=0;
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
  bt_type       = md->backtrack_type;
  turn          = md->min_loop_size;
  rtype         = &(md->rtype[0]);

  c       = vc->matrices->c_local;
  fML     = vc->matrices->fML_local;
  f3      = vc->matrices->f3_local;
  ggg     = vc->matrices->ggg_local;

  /* length = strlen(string); */
  sector[++s].i = start;
  sector[s].j   = MIN2(length, maxdist+1);
  sector[s].ml  = (bt_type=='M') ? 1 : ((bt_type=='C')?2:0);

  structure = (char *) vrna_alloc((MIN2(length-start, maxdist)+3)*sizeof(char));
  for (i=0; i<=MIN2(length-start, maxdist); i++) structure[i] = '-';

  while (s>0) {
    int ml, fij, cij, traced, i1, j1, mm, mm5, mm3, mm53, p, q, jj=0, gq=0;
    int canonical = 1;     /* (i,j) closes a canonical structure */
    i  = sector[s].i;
    j  = sector[s].j;
    ml = sector[s--].ml;   /* ml is a flag indicating if backtracking is to
                              occur in the fML- (1) or in the f-array (0) */
    if (ml==2) {
      structure[i-start] = '(';
      structure[j-start] = ')';
      goto repeat1;
    }

    if (j < i + turn + 1) continue; /* no more pairs in this interval */

    fij = (ml)? fML[i][j-i] : f3[i];

    if (ml == 0) { /* backtrack in f3 */

      if (fij == f3[i+1]) {
        sector[++s].i = i+1;
        sector[s].j   = j;
        sector[s].ml  = ml;
        continue;
      }
      /* i or i+1 is paired. Find pairing partner */
      switch(dangle_model){
        case 0:   for(traced = 0, k=j; k>i+turn; k--){

                    if(with_gquad){
                      if(fij == ggg[i][k-i] + f3[k+1]){
                        /* found the decomposition */
                        traced = i; jj = k + 1; gq = 1;
                        break;
                      }
                    }

                    jj    = k+1;
                    type  = ptype[i][k-i];
                    if(type)
                      if(fij == c[i][k-i] + E_ExtLoop(type, -1, -1, P) + f3[k+1]){
                        traced = i;
                        break;
                      }
                  }
                  break;
        case 2:   for(traced = 0, k=j; k>i+turn; k--){

                    if(with_gquad){
                      if(fij == ggg[i][k-i] + f3[k+1]){
                        /* found the decomposition */
                        traced = i; jj = k + 1; gq = 1;
                        break;
                      }
                    }

                    jj    = k+1;
                    type  = ptype[i][k-i];
                    if(type)
                      if(fij == c[i][k-i] + E_ExtLoop(type, (i>1) ? S1[i-1] : -1, (k<length) ? S1[k+1] : -1, P) + f3[k+1]){
                        traced = i;
                        break;
                      }
                  }
                  break;
        default:  for(traced = 0,k=j; k>i+turn; k--){

                    if(with_gquad){
                      if(fij == ggg[i][k-i] + f3[k+1]){
                        /* found the decomposition */
                        traced = i; jj = k + 1; gq = 1;
                        break;
                      }
                    }

                    jj = k+1;
                    type = ptype[i+1][k-(i+1)];
                    if(type){
                      if(fij == c[i+1][k-(i+1)] + E_ExtLoop(type, S1[i], -1, P) + f3[k+1]){
                        traced=i+1;
                      }
                      if(k < length){
                        if(fij == c[i+1][k-(i+1)] + E_ExtLoop(type, S1[i], S1[k+1], P) + f3[k+2]){
                          traced  = i+1;
                          jj      = k+2;
                        }
                      }
                    }
                    type = ptype[i][k-i];
                    if(type){
                      if(fij == c[i][k-i] + E_ExtLoop(type, -1, -1, P) + f3[k+1]){
                        traced = i;
                      }
                      if(k<length){
                        if(fij == c[i][k-i] + E_ExtLoop(type, -1, S1[k+1], P) + f3[k+2]){
                          traced  = i;
                          jj      = k+2;
                        }
                      }
                    }
                    if(traced) break;
                  }
                  break;
      } /* switch(dangle_model)...*/

      if (!traced) vrna_message_error("backtrack failed in f3");
      if (j==length) { /* backtrack only one component, unless j==length */
        sector[++s].i = jj;
        sector[s].j   = j;
        sector[s].ml  = ml;
      }
      i=traced; j=k;

      if(with_gquad && gq){
        /* goto backtrace of gquadruplex */
        goto repeat_gquad;
      }

      structure[i-start] = '('; structure[j-start] = ')';
      if (((jj==j+2) || (dangle_model==2)) && (j < length)) structure[j+1-start] = '.';
      goto repeat1;
    }
    else { /* trace back in fML array */
      if (fML[i][j-1-i]+P->MLbase == fij) {  /* 3' end is unpaired */
        sector[++s].i = i;
        sector[s].j   = j-1;
        sector[s].ml  = ml;
        continue;
      }
      if (fML[i+1][j-(i+1)]+P->MLbase == fij) { /* 5' end is unpaired */
        sector[++s].i = i+1;
        sector[s].j   = j;
        sector[s].ml  = ml;
        continue;
      }

      if(with_gquad){
        if(fij == ggg[i][j-i] + E_MLstem(0, -1, -1, P)){
          /* go to backtracing of quadruplex */
          goto repeat_gquad;
        }
      }

      switch(dangle_model){
        case 0:   tt = ptype[i][j-i];
                  if(fij == c[i][j-i] + E_MLstem(tt, -1, -1, P)){
                    structure[i-start] = '(';
                    structure[j-start] = ')';
                    goto repeat1;
                  }
                  break;
        case 2:   tt = ptype[i][j-i];
                  if(fij == c[i][j-i] + E_MLstem(tt, (i>1) ? S1[i-1] : -1, (j < length) ? S1[j+1] : -1, P)){
                    structure[i-start] = '(';
                    structure[j-start] = ')';
                    goto repeat1;
                  }
                  break;
        default:  tt = ptype[i][j-i];
                  if(fij == c[i][j-i] + E_MLstem(tt, -1, -1, P)){
                    structure[i-start] = '(';
                    structure[j-start] = ')';
                    goto repeat1;
                  }
                  tt = ptype[i+1][j-(i+1)];
                  if(fij == c[i+1][j-(i+1)] + E_MLstem(tt, S1[i], -1, P) + P->MLbase){
                    structure[++i-start] = '(';
                    structure[j-start] = ')';
                    goto repeat1;
                  }
                  tt = ptype[i][j-1-i];
                  if(fij == c[i][j-1-i] + E_MLstem(tt, -1, S1[j], P) + P->MLbase){
                    structure[i-start] = '(';
                    structure[--j-start] = ')';
                    goto repeat1;
                  }
                  tt = ptype[i+1][j-1-(i+1)];
                  if(fij == c[i+1][j-1-(i+1)] + E_MLstem(tt, S1[i], S1[j], P) + 2*P->MLbase){
                    structure[++i-start] = '(';
                    structure[--j-start] = ')';
                    goto repeat1;
                  }
                  break;
      } /* switch(dangle_model)... */

      /* modular decomposition */
      for (k = i+1+turn; k <= j-2-turn; k++)
        if (fij == (fML[i][k-i]+fML[k+1][j-(k+1)]))
          break;

      if ((dangle_model==3)&&(k>j-2-turn)) { /* must be coax stack */
        ml = 2;
        for (k = i+1+turn; k <= j-2-turn; k++) {
          type = ptype[i][k-i];  type= rtype[type];
          type_2 = ptype[k+1][j-(k+1)]; type_2= rtype[type_2];
          if (type && type_2)
            if (fij == c[i][k-i]+c[k+1][j-(k+1)]+P->stack[type][type_2]+
                       2*P->MLintern[1])
              break;
        }
      }

      sector[++s].i = i;
      sector[s].j   = k;
      sector[s].ml  = ml;
      sector[++s].i = k+1;
      sector[s].j   = j;
      sector[s].ml  = ml;

      if (k>j-2-turn) vrna_message_error("backtrack failed in fML");
      continue;
    }

  repeat1:

    /*----- begin of "repeat:" -----*/
    if (canonical)  cij = c[i][j-i];

    type = ptype[i][j-i];


    if (noLP)
      if (cij == c[i][j-i]) {
        /* (i.j) closes canonical structures, thus
           (i+1.j-1) must be a pair                */
        type_2 = ptype[i+1][j-1-(i+1)]; type_2 = rtype[type_2];
        cij -= P->stack[type][type_2];
        structure[i+1-start] = '('; structure[j-1-start] = ')';
        i++; j--;
        canonical=0;
        goto repeat1;
      }
    canonical = 1;


    no_close = (((type==3)||(type==4))&&noGUclosure);
    if (no_close) {
      if (cij == FORBIDDEN) continue;
    } else
      if (cij == E_Hairpin(j-i-1, type, S1[i+1], S1[j-1],string+i-1, P))
        continue;

    for (p = i+1; p <= MIN2(j-2-turn,i+MAXLOOP+1); p++) {
      int minq;
      minq = j-i+p-MAXLOOP-2;
      if (minq<p+1+turn) minq = p+1+turn;
      for (q = j-1; q >= minq; q--) {

        type_2 = ptype[p][q-p];
        if (type_2==0) continue;
        type_2 = rtype[type_2];
        if (noGUclosure)
          if (no_close||(type_2==3)||(type_2==4))
            if ((p>i+1)||(q<j-1)) continue;  /* continue unless stack */

        /* energy = oldLoopEnergy(i, j, p, q, type, type_2); */
        energy = E_IntLoop(p-i-1, j-q-1, type, type_2, S1[i+1], S1[j-1], S1[p-1], S1[q+1],P);

        new = energy+c[p][q-p];
        traced = (cij == new);
        if (traced) {
          structure[p-start] = '(';
          structure[q-start] = ')';
          i = p, j = q;
          goto repeat1;
        }
      }
    }

    /* end of repeat: --------------------------------------------------*/

    /* (i.j) must close a multi-loop */
    tt = rtype[type];
    i1 = i+1; j1 = j-1;

    if(with_gquad){
      /*
        The case that is handled here actually resembles something like
        an interior loop where the enclosing base pair is of regular
        kind and the enclosed pair is not a canonical one but a g-quadruplex
        that should then be decomposed further...
      */
      if(backtrack_GQuad_IntLoop_L(cij, i, j, type, S, ggg, maxdist, &p, &q, P)){
        i = p; j = q;
        goto repeat_gquad;
      }
    }

    sector[s+1].ml  = sector[s+2].ml = 1;

    switch(dangle_model){
      case 0:   mm = P->MLclosing + E_MLstem(tt, -1, -1, P);
                for(k = i+2+turn; k < j-2-turn; k++){
                  if(cij == fML[i+1][k-(i+1)] + fML[k+1][j-1-(k+1)] + mm)
                    break;
                }
                break;
      case 2:   mm = P->MLclosing + E_MLstem(tt, S1[j-1], S1[i+1], P);
                for(k = i+2+turn; k < j-2-turn; k++){
                  if(cij == fML[i+1][k-(i+1)] + fML[k+1][j-1-(k+1)] + mm)
                    break;
                }
                break;
      default:  mm    = P->MLclosing + E_MLstem(tt, -1, -1, P);
                mm5   = P->MLclosing + E_MLstem(tt, S1[j-1], -1, P) + P->MLbase;
                mm3   = P->MLclosing + E_MLstem(tt, -1, S1[i+1], P) + P->MLbase;
                mm53  = P->MLclosing + E_MLstem(tt, S1[j-1], S1[i+1], P) + 2*P->MLbase;
                for(k = i+2+turn; k < j-2-turn; k++){
                  if(cij == fML[i+1][k-(i+1)] + fML[k+1][j-1-(k+1)] + mm)
                    break;
                  else if(cij == fML[i+2][k-(i+2)] + fML[k+1][j-1-(k+1)] + mm3){
                    i1 = i+2;
                    break;
                  }
                  else if(cij == fML[i+1][k-(i+1)] + fML[k+1][j-2-(k+1)] + mm5){
                    j1 = j-2;
                    break;
                  }
                  else if(cij == fML[i+2][k-(i+2)] + fML[k+1][j-2-(k+1)] + mm53){
                    i1 = i+2;
                    j1 = j-2;
                    break;
                  }
                  /* coaxial stacking of (i.j) with (i+1.k) or (k.j-1) */
                  /* use MLintern[1] since coax stacked pairs don't get TerminalAU */
                  if (dangle_model==3) {
                    int en;
                    type_2 = ptype[i+1][k-(i+1)]; type_2 = rtype[type_2];
                    if (type_2) {
                      en = c[i+1][k-(i+1)]+P->stack[type][type_2]+fML[k+1][j-1-(k+1)];
                      if (cij == en+2*P->MLintern[1]+P->MLclosing) {
                        ml = 2;
                        sector[s+1].ml  = 2;
                        break;
                      }
                    }
                    type_2 = ptype[k+1][j-1-(k+1)]; type_2 = rtype[type_2];
                    if (type_2) {
                      en = c[k+1][j-1-(k+1)]+P->stack[type][type_2]+fML[i+1][k-(i+1)];
                      if (cij == en+2*P->MLintern[1]+P->MLclosing) {
                        sector[s+2].ml = 2;
                        break;
                      }
                    }
                  }
                }
                break;
    } /* switch(dangle_model)... */

    if (k<=j-3-turn) { /* found the decomposition */
      sector[++s].i = i1;
      sector[s].j   = k;
      sector[++s].i = k+1;
      sector[s].j   = j1;
    } else {
#if 0
      /* Y shaped ML loops fon't work yet */
      if (dangle_model==3) {
        /* (i,j) must close a Y shaped ML loop with coax stacking */
        if (cij ==  fML[i+1][j-2-(i+2)] + mm + d3 + d5 + P->MLbase + P->MLbase) {
          i1 = i+2;
          j1 = j-2;
        } else if (cij ==  fML[i+1][j-2-(i+1)] + mm + d5 + P->MLbase)
          j1 = j-2;
        else if (cij ==  fML[i+2][j-1-(i+2)] + mm + d3 + P->MLbase)
          i1 = i+2;
        else /* last chance */
          if (cij != fML[i+1][j-1-(i+1)] + mm + P->MLbase)
            fprintf(stderr,  "backtracking failed in repeat");
        /* if we arrive here we can express cij via fML[i1,j1]+dangle_model */
        sector[++s].i = i1;
        sector[s].j   = j1;
      }
      else
#endif
        vrna_message_error("backtracking failed in repeat");
    }

    continue; /* this is a workarround to not accidentally proceed in the following block */

  repeat_gquad:
    /*
      now we do some fancy stuff to backtrace the stacksize and linker lengths
      of the g-quadruplex that should reside within position i,j
    */
    {
      int l[3], L, a;
      L = -1;

      get_gquad_pattern_mfe(S, i, j, P, &L, l);
      if(L != -1){
        /* fill the G's of the quadruplex into the structure string */
        for(a=0;a<L;a++){
          structure[i+a-start] = '+';
          structure[i+L+l[0]+a-start] = '+';
          structure[i+L+l[0]+L+l[1]+a-start] = '+';
          structure[i+L+l[0]+L+l[1]+L+l[2]+a-start] = '+';
        }
        goto repeat_gquad_exit;
      }
      vrna_message_error("backtracking failed in repeat_gquad");
    }
  repeat_gquad_exit:
    asm("nop");

  }

  for (i=strlen(structure)-1; i>0 && structure[i] == '-'; i--)
    structure[i] = '\0';
  for (;i>=0; i--)
   if (structure[i]=='-') structure[i]='.';

  return structure;
}

PRIVATE void
make_ptypes(vrna_fold_compound_t *vc, int i){

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

  for(k = turn + 1; k < maxdist; k++){
    j = i + k;
    if (j > n)
      continue;
    type = md->pair[S[i]][S[j]];

    if(noLP && type){
      if(!ptype[i + 1][j - 1 - i - 1])
        if(j == n || i == 1 || (!md->pair[S[i - 1]][S[j + 1]]))
          type = 0;
    }
    ptype[i][j - i] = type;
  }
}

/*###########################################*/
/*# deprecated functions below              #*/
/*###########################################*/

#ifdef  VRNA_BACKWARD_COMPAT

PUBLIC float Lfold( const char *string,
                    char *structure,
                    int window_size){

  float               energy;
  vrna_fold_compound_t  *vc;
  vrna_md_t           md;

  vrna_md_set_globals(&md);

  md.window_size = window_size;
  md.max_bp_span = window_size;

  vc  = vrna_fold_compound(string, &md, VRNA_OPTION_MFE | VRNA_OPTION_WINDOW);

  energy = wrap_Lfold(vc, 0, 0.0, NULL);

  vrna_fold_compound_free(vc);

  return energy;
}

PUBLIC float
Lfoldz( const char *string,
        char *structure,
        int window_size,
        int zsc,
        double min_z){

  float               energy;
  vrna_fold_compound_t  *vc;
  vrna_md_t           md;

  vrna_md_set_globals(&md);

  md.window_size = window_size;
  md.max_bp_span = window_size;

  vc  = vrna_fold_compound(string, &md, VRNA_OPTION_MFE | VRNA_OPTION_WINDOW);

#ifndef USE_SVM
  zsc = 0;  /* deactivate z-scoring if no compiled-in svm support is available */
#endif

  energy = wrap_Lfold(vc, zsc, min_z, NULL);

  vrna_fold_compound_free(vc);

  return energy;
}

#endif
