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
#include "utils.h"
#include "energy_par.h"
#include "fold_vars.h"
#include "pair_mat.h"
#include "params.h"
#include "loop_energies.h"
#include "Lfold.h"

#ifdef USE_SVM
#include "svm.h"
#include "svm_utils.h"
#endif

#ifdef _OPENMP
#include <omp.h>
#endif


/*@unused@*/
static char rcsid[] UNUSED = "$Id: Lfold.c,v 1.9 2007/09/04 09:20:12 ivo Exp $";


#define PAREN

#define STACK_BULGE1      1   /* stacking energies for bulges of size 1 */
#define NEW_NINIO         1   /* new asymetry penalty */
/*@unused@*/
#define MAXSECTORS        500     /* dimension for a backtrack array */
#define LOCALITY          0.      /* locality parameter for base-pairs */

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
PRIVATE paramT        *P = NULL;
PRIVATE int           **c = NULL;        /* energy array, given that i-j pair */
PRIVATE int           *cc = NULL;        /* linear array for calculating canonical structures */
PRIVATE int           *cc1 = NULL;       /*   "     "        */
PRIVATE int           *f3 = NULL;        /* energy of 5' end */
PRIVATE int           **fML = NULL;      /* multi-loop auxiliary energy array */
PRIVATE int           *Fmi = NULL;       /* holds row i of fML (avoids jumps in memory) */
PRIVATE int           *DMLi = NULL;      /* DMLi[j] holds MIN(fML[i,k]+fML[k+1,j])  */
PRIVATE int           *DMLi1 = NULL;     /*             MIN(fML[i+1,k]+fML[k+1,j])  */
PRIVATE int           *DMLi2 = NULL;     /*             MIN(fML[i+2,k]+fML[k+1,j])  */
PRIVATE char          **ptype = NULL;    /* precomputed array of pair types */
PRIVATE short         *S = NULL, *S1 = NULL;
PRIVATE unsigned int  length;
PRIVATE char          *prev = NULL;
#ifdef USE_SVM
PRIVATE struct svm_model  *avg_model = NULL;
PRIVATE struct svm_model  *sd_model = NULL;
#endif

#ifdef _OPENMP

/* NOTE: all variables are assumed to be uninitialized if they are declared as threadprivate
         thus we have to initialize them before usage by a seperate function!
         OR: use copyin in the PARALLEL directive!
         e.g.:
         #pragma omp parallel for copyin(P, ...)
*/

#ifdef USE_SVM
#pragma omp threadprivate(P, c, cc, cc1, f3, fML, Fmi, DMLi, DMLi1, DMLi2, ptype, S, S1, length, sd_model, avg_model)
#else
#pragma omp threadprivate(P, c, cc, cc1, f3, fML, Fmi, DMLi, DMLi1, DMLi2, ptype, S, S1, length)
#endif

#endif

/*
#################################
# PRIVATE FUNCTION DECLARATIONS #
#################################
*/
PRIVATE void  initialize_Lfold(int length, int maxdist);
PRIVATE void  update_fold_params(void);
PRIVATE void  get_arrays(unsigned int size, int maxdist);
PRIVATE void  free_arrays(int maxdist);
PRIVATE void  make_ptypes(const short *S, int i, int maxdist, int n);
PRIVATE char  *backtrack(const char *sequence, int start, int maxdist);
PRIVATE int   fill_arrays(const char *sequence, int maxdist, int zsc, double min_z);

/*
#################################
# BEGIN OF FUNCTION DEFINITIONS #
#################################
*/

/*--------------------------------------------------------------------------*/
PRIVATE void initialize_Lfold(int length, int maxdist){

  if (length<1) nrerror("initialize_Lfold: argument must be greater 0");
  get_arrays((unsigned) length, maxdist);
  update_fold_params();
}

/*--------------------------------------------------------------------------*/
PRIVATE void get_arrays(unsigned int size, int maxdist){
  int i;
  c     = (int **)  space(sizeof(int *) *(size+1));
  fML   = (int **)  space(sizeof(int *) *(size+1));
  ptype = (char **) space(sizeof(char *)*(size+1));
  f3    = (int *)   space(sizeof(int)   *(size+2));  /* has to be one longer */
  cc    = (int *)   space(sizeof(int)   *(maxdist+5));
  cc1   = (int *)   space(sizeof(int)   *(maxdist+5));
  Fmi   = (int *)   space(sizeof(int)   *(maxdist+5));
  DMLi  = (int *)   space(sizeof(int)   *(maxdist+5));
  DMLi1 = (int *)   space(sizeof(int)   *(maxdist+5));
  DMLi2 = (int *)   space(sizeof(int)   *(maxdist+5));
  for (i=size; (i>(int)size-maxdist-5) && (i>=0); i--) {
    c[i] = (int *) space(sizeof(int)*(maxdist+5));
  }
  for (i=size; (i>(int)size-maxdist-5) && (i>=0); i--) {
    fML[i] = (int *) space(sizeof(int)*(maxdist+5));
  }
  for (i=size; (i>(int)size-maxdist-5) && (i>=0); i--) {
    ptype[i] = (char *) space(sizeof(char)*(maxdist+5));
  }
}

/*--------------------------------------------------------------------------*/

PRIVATE void free_arrays(int maxdist){
  int i;
  for (i=0; (i<maxdist+5) && (i<=length); i++){
    free(c[i]);
    free(fML[i]);
    free(ptype[i]);
  }
  free(c);
  free(fML);
  free(ptype);
  free(f3);
  free(cc);
  free(cc1);
  free(Fmi);
  free(DMLi);
  free(DMLi1);
  free(DMLi2);

  f3    = cc = cc1 = Fmi = DMLi = DMLi1 = DMLi2 = NULL;
  c     = fML = NULL;
  ptype = NULL;
}

/*--------------------------------------------------------------------------*/

PUBLIC  float Lfold(const char *string, char *structure, int maxdist){
  return Lfoldz(string, structure, maxdist, 0, 0.0);
}

PUBLIC  float Lfoldz(const char *string, char *structure, int maxdist, int zsc, double min_z){
  int i, energy;

  length = (int) strlen(string);
  if (maxdist>length) maxdist = length;
  initialize_Lfold(length, maxdist);
  if (fabs(P->temperature - temperature)>1e-6) update_fold_params();

  S   = encode_sequence(string, 0);
  S1  = encode_sequence(string, 1);

  for (i=length; i>=(int)length-(int)maxdist-4 && i>0; i--)
    make_ptypes(S, i, maxdist, length);

#ifdef USE_SVM  /*svm*/
  if(zsc){
    avg_model = svm_load_model_string(avg_model_string);
    sd_model  = svm_load_model_string(sd_model_string);
  }
#endif

  energy = fill_arrays(string, maxdist, zsc, min_z);

#ifdef USE_SVM  /*svm*/
  if(zsc){
    svm_destroy_model(avg_model);
    svm_destroy_model(sd_model);
  }
#endif

  free(S); free(S1);
  free_arrays(maxdist);

  return (float) energy/100.;
}

PRIVATE int fill_arrays(const char *string, int maxdist, int zsc, double min_z) {
  /* fill "c", "fML" and "f3" arrays and return  optimal energy */

  int   i, j, k, length, energy;
  int   decomp, new_fML;
  int   no_close, type, type_2, tt;
  int   fij;
  int   lind;
  short mm5, mm3;

  length = (int) strlen(string);
  prev=NULL;
  for (j=0; j<maxdist+5; j++)
    Fmi[j]=DMLi[j]=DMLi1[j]=DMLi2[j]=INF;
  for (j=length; j>length-maxdist-4; j--) {
    for (i=(length-maxdist-4>0)?length-maxdist-4:1 ; i<j; i++)
      c[i][j-i] = fML[i][j-i] = INF;
  }

  for (i = length-TURN-1; i >= 1; i--) { /* i,j in [1..length] */
    for (j = i+TURN+1; j <= length && j <= i+maxdist; j++) {
      int p, q;
      type = ptype[i][j-i];

      no_close = (((type==3)||(type==4))&&no_closingGU);

      if (type) {   /* we have a pair */
        int new_c=0, stackEnergy=INF;
        /* hairpin ----------------------------------------------*/

        new_c = (no_close) ? FORBIDDEN : E_Hairpin(j-i-1, type, S1[i+1], S1[j-1], string+i-1, P);

        /*--------------------------------------------------------
          check for elementary structures involving more than one
          closing pair.
          --------------------------------------------------------*/

        for (p = i+1; p <= MIN2(j-2-TURN,i+MAXLOOP+1) ; p++){
          int minq = j-i+p-MAXLOOP-2;
          if (minq<p+1+TURN) minq = p+1+TURN;
          for (q = minq; q < j; q++) {
            type_2 = ptype[p][q-p];

            if (type_2==0) continue;
            type_2 = rtype[type_2];

            if (no_closingGU)
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
          switch(dangles){
            /* no dangles */
            case 0:   decomp += E_MLstem(tt, -1, -1, P);
                      break;
            /* double dangles */
            case 2:   decomp += E_MLstem(tt, S1[j-1], S1[i+1], P);
                      break;
            /* normal dangles, aka dangles = 1 */
            default:  decomp += E_MLstem(tt, -1, -1, P);
                      decomp = MIN2(decomp, DMLi2[j-1-(i+2)] + E_MLstem(tt, -1, S1[i+1], P) + P->MLbase);
                      decomp = MIN2(decomp, DMLi2[j-2-(i+2)] + E_MLstem(tt, S1[j-1], S1[i+1], P) + 2*P->MLbase);
                      decomp = MIN2(decomp, DMLi1[j-2-(i+1)] + E_MLstem(tt, S1[j-1], -1, P) + P->MLbase);
                      break;
          }
          new_c = MIN2(new_c, decomp + P->MLclosing);
        }

        /* coaxial stacking of (i.j) with (i+1.k) or (k+1.j-1) */

        if (dangles==3) {
          decomp = INF;
          for (k = i+2+TURN; k < j-2-TURN; k++) {
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

        new_c = MIN2(new_c, cc1[j-1-(i+1)]+stackEnergy);
        cc[j-i] = new_c;
        if (noLonelyPairs)
          c[i][j-i] = cc1[j-1-(i+1)]+stackEnergy;
        else
          c[i][j-i] = cc[j-i];

      } /* end >> if (pair) << */

      else c[i][j-i] = INF;


      /* done with c[i,j], now compute fML[i,j] */
      /* free ends ? -----------------------------------------*/
      new_fML = INF;
      switch(dangles){
        /* no dangles */
        case 0:   new_fML = fML[i+1][j-i-1] + P->MLbase;
                  new_fML = MIN2(new_fML, fML[i][j-1-i] + P->MLbase);
                  new_fML = MIN2(new_fML, c[i][j-i] + E_MLstem(type, -1, -1, P));
                  break;
        /* double dangles */
        case 2:   new_fML = fML[i+1][j-i-1] + P->MLbase;
                  new_fML = MIN2(fML[i][j-1-i] + P->MLbase, new_fML);
                  new_fML = MIN2(new_fML,  c[i][j-i] + E_MLstem(type, (i>1) ? S1[i-1] : -1, (j<length) ? S1[j+1] : -1, P));
                  break;
        /* normal dangles, aka dangles = 1 */
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

      /* modular decomposition -------------------------------*/
      for (decomp = INF, k = i+1+TURN; k <= j-2-TURN; k++)
        decomp = MIN2(decomp, Fmi[k-i]+fML[k+1][j-k-1]);

      DMLi[j-i] = decomp;               /* store for use in ML decompositon */
      new_fML   = MIN2(new_fML, decomp);

      /* coaxial stacking */
      if (dangles==3) {
        /* additional ML decomposition as two coaxially stacked helices */
        for (decomp = INF, k = i+1+TURN; k <= j-2-TURN; k++) {
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
      f3[i] = f3[i+1];
      switch(dangles){
        /* dont use dangling end and mismatch contributions at all */
        case 0:   for(j=i+TURN+1; j<length && j<=i+maxdist; j++){
                    type = ptype[i][j-i];
                    if(type)
                      f3[i] = MIN2(f3[i], f3[j+1] + c[i][j-i] + E_ExtLoop(type, -1, -1, P));
                  }
                  if(length<=i+maxdist){
                    j=length;
                    type = ptype[i][j-i];
                    if(type)
                      f3[i] = MIN2(f3[i], c[i][j-i] + E_ExtLoop(type, -1, -1, P));
                  }
                  break;
        /* always use dangles on both sides */
        case 2:   for(j=i+TURN+1; j<length && j<=i+maxdist; j++){
                    type = ptype[i][j-i];
                    if(type)
                      f3[i] = MIN2(f3[i], f3[j+1] + c[i][j-i] + E_ExtLoop(type, (i>1) ? S1[i-1] : -1, S1[j+1], P));
                  }
                  if(length<=i+maxdist){
                    j=length;
                    type = ptype[i][j-i];
                    if(type)
                      f3[i] = MIN2(f3[i], c[i][j-i] + E_ExtLoop(type, (i>1) ? S1[i-1] : -1, -1, P));
                  }
                  break;
        /* normal dangles, aka dangles = 1 */
        default:  for(j=i+TURN+1; j<length && j<=i+maxdist; j++){
                    type = ptype[i][j-i];
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
                    type  = ptype[i][j-i];
                    if(type)
                      f3[i] = MIN2(f3[i], c[i][j-i] + E_ExtLoop(type, -1, -1, P));
                    type  = ptype[i+1][j-i-1];
                    if(type)
                      f3[i] = MIN2(f3[i], c[i+1][j-i-1] + E_ExtLoop(type, S1[i], -1, P));
                  }
                  break;
      } /* switch(dangles)... */

      /* backtrack partial structure */
      if (f3[i] != f3[i+1] && (f3[i+1] < 0)) do_backtrack=1;
      else if (do_backtrack) {
        int pairpartner; /*i+1?? is paired with pairpartner*/
        int en, cc;
        int traced2=0;
        fij = f3[i+1];
        lind=i+1;

        /*start "short" backtrack*/

        /*get paired base*/
        while(fij==f3[lind+1]) lind++;

        /*get pairpartner*/
        for (pairpartner = lind + TURN; pairpartner <= lind + maxdist; pairpartner++){
          type = ptype[lind][pairpartner-lind];
          switch(dangles){
            case 0:   if(type){
                        cc = c[lind][pairpartner-lind] + E_ExtLoop(type, -1, -1, P);
                        if(fij == cc + f3[pairpartner + 1])
                          traced2 = 1;
                      }
                      break;
            case 2:   if(type){
                        cc = c[lind][pairpartner-lind] + E_ExtLoop(type, (lind > 1) ? S1[lind-1] : -1, (pairpartner < length) ? S1[pairpartner+1] : -1, P);
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
        if (!traced2) nrerror("backtrack failed in short backtrack");
        if (zsc){
#ifdef USE_SVM
          int info_avg;
          double average_free_energy;
          double sd_free_energy;
          double my_z;
          int *AUGC = get_seq_composition(S, lind-1, MIN2((pairpartner+1),length));
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
                ss =  backtrack(string, lind , pairpartner+1);
                if (prev) {
                  if ((i+strlen(ss)<prev_i+strlen(prev)) ||
                      strncmp(ss+prev_i-i,prev,strlen(prev))) { /* ss does not contain prev */
                    if (dangles==2)
                      printf(".%s (%6.2f) %4d z= %.3f\n", prev, (f3[prev_i]-f3[prev_i+strlen(prev)-1])/100., prev_i-1, prevz);
                    else
                      printf("%s (%6.2f) %4d z=%.3f\n ", prev, (f3[prev_i]-f3[prev_i+strlen(prev)])/100., prev_i, prevz);
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
          ss =  backtrack(string, lind , pairpartner+1);
          if (prev) {
            if ((i+strlen(ss)<prev_i+strlen(prev)) || strncmp(ss+prev_i-i,prev,strlen(prev))){
              /* ss does not contain prev */
              if (dangles==2)
                printf(".%s (%6.2f) %4d\n", prev, (f3[prev_i]-f3[prev_i+strlen(prev)-1])/100., prev_i-1);
              else
                printf("%s (%6.2f) %4d\n", prev, (f3[prev_i]-f3[prev_i+strlen(prev)])/100., prev_i);
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
            if (dangles==2)
              printf(".%s (%6.2f) %4d z= %.2f\n", prev, (f3[prev_i]-f3[prev_i+strlen(prev)-1])/100., prev_i-1, prevz);
           else
              printf("%s (%6.2f) %4dz= %.2f \n", prev, (f3[prev_i]-f3[prev_i+strlen(prev)])/100., prev_i, prevz);
            free(prev); prev=NULL;
          }
          else {
            if (dangles==2)
              printf(".%s (%6.2f) %4d\n", prev, (f3[prev_i]-f3[prev_i+strlen(prev)-1])/100., prev_i-1);
            else
              printf("%s (%6.2f) %4d\n", prev, (f3[prev_i]-f3[prev_i+strlen(prev)])/100., prev_i);
          }
        } else if (f3[i]<0)do_backtrack=1;

        if (do_backtrack) {
          int pairpartner; /*i+1?? is paired with pairpartner*/
          int en, cc;
          double average_free_energy;
          double sd_free_energy;
          int info_avg;
          double my_z;
          int traced2 = 0;
          fij = f3[i];
          lind=i;
          while(fij==f3[lind+1]) lind++;

          /*get pairpartner*/
          for(pairpartner = lind + TURN; pairpartner <= lind + maxdist; pairpartner++){
            type = ptype[lind][pairpartner-lind];
            switch(dangles){
              case 0:   if(type){
                          cc = c[lind][pairpartner-lind] + E_ExtLoop(type, -1, -1, P);
                          if(fij == cc + f3[pairpartner + 1])
                            traced2 = 1;
                        }
                        break;
              case 2:   if(type){
                          cc = c[lind][pairpartner-lind] + E_ExtLoop(type, (lind > 1) ? S1[lind-1] : -1, (pairpartner < length) ? S1[pairpartner+1] : -1, P);
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
          if (!traced2) nrerror("backtrack failed in short backtrack");

          if(zsc){
#ifdef USE_SVM
            int *AUGC = get_seq_composition(S, lind-1, MIN2((pairpartner+1),length));
            average_free_energy = avg_regression(AUGC[0],AUGC[1],AUGC[2],AUGC[3],AUGC[4],avg_model,&info_avg);
            if (info_avg == 0)  {
              double difference;
              double min_sd = minimal_sd(AUGC[0],AUGC[1],AUGC[2],AUGC[3],AUGC[4]);
              difference=(fij-f3[pairpartner+1])/100.-average_free_energy;
              if ( difference - ( min_z * min_sd ) <= 0.0001 ) {
                sd_free_energy = sd_regression(AUGC[0],AUGC[1],AUGC[2],AUGC[3],AUGC[4],sd_model);
                my_z=difference/sd_free_energy;
                if (my_z<=min_z){
                  ss =  backtrack(string, lind , pairpartner+1);
                  printf("%s (%6.2f) %4d z= %.2f\n", ss, (f3[lind]-f3[lind+strlen(ss)-1])/100., lind, my_z);
                }
              }
            }
            free(AUGC);
#endif
          }
          else {
            ss =  backtrack(string, lind , pairpartner+1);
            if (dangles==2)
              printf("%s (%6.2f) %4d\n", ss, (f3[lind]-f3[lind+strlen(ss)-1])/100., 1);
            else
              printf("%s (%6.2f) %4d\n", ss, (f3[lind]-f3[lind+strlen(ss)])/100., 1);
            free(ss);
          }
        }
        do_backtrack=0;
      }
    }
    {
      int ii, *FF; /* rotate the auxilliary arrays */
      FF = DMLi2; DMLi2 = DMLi1; DMLi1 = DMLi; DMLi = FF;
      FF = cc1; cc1=cc; cc=FF;
      for (j=0; j< maxdist+5; j++) {cc[j]=Fmi[j]=DMLi[j]=INF; }
      if (i+maxdist+4<=length) {
        c[i-1] = c[i+maxdist+4]; c[i+maxdist+4] = NULL;
        fML[i-1] = fML[i+maxdist+4]; fML[i+maxdist+4]=NULL;
        ptype[i-1] = ptype[i+maxdist+4]; ptype[i+maxdist+4] = NULL;
        if (i>1) make_ptypes(S, i-1, maxdist, length);
        for (ii=0; ii<maxdist+5; ii++) {
          c[i-1][ii] = INF;
          fML[i-1][ii] = INF;
        }
      }
    }
  }

  return f3[1];
}

PRIVATE char *backtrack(const char *string, int start, int maxdist){
  /*------------------------------------------------------------------
    trace back through the "c", "f3" and "fML" arrays to get the
    base pairing list. No search for equivalent structures is done.
    This is fast, since only few structure elements are recalculated.
    ------------------------------------------------------------------*/
  sect  sector[MAXSECTORS];   /* backtracking sectors */
  int   i, j, k, energy, new, no_close, type, type_2, tt, s=0;
  char  *structure;

  /* length = strlen(string); */
  sector[++s].i = start;
  sector[s].j   = MIN2(length, maxdist+1);
  sector[s].ml  = (backtrack_type=='M') ? 1 : ((backtrack_type=='C')?2:0);

  structure = (char *) space((MIN2(length-start, maxdist)+3)*sizeof(char));
  for (i=0; i<=MIN2(length-start, maxdist); i++) structure[i] = '-';

  while (s>0) {
    int ml, fij, cij, traced, i1, j1, d3, d5, mm, mm5, mm3, mm53, p, q, jj=0;
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

    if (j < i+TURN+1) continue; /* no more pairs in this interval */

    fij = (ml)? fML[i][j-i] : f3[i];

    if (ml == 0) { /* backtrack in f3 */

      if (fij == f3[i+1]) {
        sector[++s].i = i+1;
        sector[s].j   = j;
        sector[s].ml  = ml;
        continue;
      }
      /* i or i+1 is paired. Find pairing partner */
      switch(dangles){
        case 0:   for(traced = 0, k=j; k>i+TURN; k--){
                    jj    = k+1;
                    type  = ptype[i][k-i];
                    if(type)
                      if(fij == c[i][k-i] + E_ExtLoop(type, -1, -1, P) + f3[k+1]){
                        traced = i;
                        break;
                      }
                  }
                  break;
        case 2:   for(traced = 0, k=j; k>i+TURN; k--){
                    jj    = k+1;
                    type  = ptype[i][k-i];
                    if(type)
                      if(fij == c[i][k-i] + E_ExtLoop(type, (i>1) ? S1[i-1] : -1, (k<length) ? S1[k+1] : -1, P) + f3[k+1]){
                        traced = i;
                        break;
                      }
                  }
                  break;
        default:  for(traced = 0,k=j; k>i+TURN; k--){
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
      } /* switch(dangles)...*/

      if (!traced) nrerror("backtrack failed in f3");
      if (j==length) { /* backtrack only one component, unless j==length */
        sector[++s].i = jj;
        sector[s].j   = j;
        sector[s].ml  = ml;
      }
      i=traced; j=k;
      structure[i-start] = '('; structure[j-start] = ')';
      if (((jj==j+2) || (dangles==2)) && (j < length)) structure[j+1-start] = '.';
      goto repeat1;
    }
    else { /* trace back in fML array */
      int cij1=INF, ci1j=INF, ci1j1=INF;

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


      switch(dangles){
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
      } /* switch(dangles)... */

      /* modular decomposition */
      for (k = i+1+TURN; k <= j-2-TURN; k++)
        if (fij == (fML[i][k-i]+fML[k+1][j-(k+1)]))
          break;

      if ((dangles==3)&&(k>j-2-TURN)) { /* must be coax stack */
        ml = 2;
        for (k = i+1+TURN; k <= j-2-TURN; k++) {
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

      if (k>j-2-TURN) nrerror("backtrack failed in fML");
      continue;
    }

  repeat1:

    /*----- begin of "repeat:" -----*/
    if (canonical)  cij = c[i][j-i];

    type = ptype[i][j-i];


    if (noLonelyPairs)
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


    no_close = (((type==3)||(type==4))&&no_closingGU);
    if (no_close) {
      if (cij == FORBIDDEN) continue;
    } else
      if (cij == E_Hairpin(j-i-1, type, S1[i+1], S1[j-1],string+i-1, P))
        continue;

    for (p = i+1; p <= MIN2(j-2-TURN,i+MAXLOOP+1); p++) {
      int minq;
      minq = j-i+p-MAXLOOP-2;
      if (minq<p+1+TURN) minq = p+1+TURN;
      for (q = j-1; q >= minq; q--) {

        type_2 = ptype[p][q-p];
        if (type_2==0) continue;
        type_2 = rtype[type_2];
        if (no_closingGU)
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
    sector[s+1].ml  = sector[s+2].ml = 1;

    switch(dangles){
      case 0:   mm = P->MLclosing + E_MLstem(tt, -1, -1, P);
                for(k = i+2+TURN; k < j-2-TURN; k++){
                  if(cij == fML[i+1][k-(i+1)] + fML[k+1][j-1-(k+1)] + mm)
                    break;
                }
                break;
      case 2:   mm = P->MLclosing + E_MLstem(tt, S1[j-1], S1[i+1], P);
                for(k = i+2+TURN; k < j-2-TURN; k++){
                  if(cij == fML[i+1][k-(i+1)] + fML[k+1][j-1-(k+1)] + mm)
                    break;
                }
                break;
      default:  mm    = P->MLclosing + E_MLstem(tt, -1, -1, P);
                mm5   = P->MLclosing + E_MLstem(tt, S1[j-1], -1, P) + P->MLbase;
                mm3   = P->MLclosing + E_MLstem(tt, -1, S1[i+1], P) + P->MLbase;
                mm53  = P->MLclosing + E_MLstem(tt, S1[j-1], S1[i+1], P) + 2*P->MLbase;
                for(k = i+2+TURN; k < j-2-TURN; k++){
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
                  if (dangles==3) {
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
    } /* switch(dangles)... */

    if (k<=j-3-TURN) { /* found the decomposition */
      sector[++s].i = i1;
      sector[s].j   = k;
      sector[++s].i = k+1;
      sector[s].j   = j1;
    } else {
#if 0
      /* Y shaped ML loops fon't work yet */
      if (dangles==3) {
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
        /* if we arrive here we can express cij via fML[i1,j1]+dangles */
        sector[++s].i = i1;
        sector[s].j   = j1;
      }
      else
#endif
        nrerror("backtracking failed in repeat");
    }

  }

  for (i=strlen(structure)-1; i>0 && structure[i] == '-'; i--)
    structure[i] = '\0';
  for (;i>=0; i--)
   if (structure[i]=='-') structure[i]='.';

  return structure;
}


PRIVATE void update_fold_params(void){
  if(P) free(P);
  P = scale_parameters();
  make_pair_matrix();
}

/*---------------------------------------------------------------------------*/

PRIVATE void make_ptypes(const short *S, int i, int maxdist, int n) {
  int j,k, type;

  for (k=TURN+1; k<maxdist; k++) {
    j = i+k;
    if (j>n) continue;
    type = pair[S[i]][S[j]];
    if (noLonelyPairs && type) {
      if (!ptype[i+1][j-1-i-1])
        if (j==n || i==1 || (!pair[S[i-1]][S[j+1]])) type=0;
    }
    ptype[i][j-i]=type;
  }
}
