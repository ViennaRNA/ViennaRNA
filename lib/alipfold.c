/* Last changed Time-stamp: <2009-02-24 14:37:05 ivo> */
/*
                  partiton function and base pair probabilities
                  for RNA secvondary structures
                  of a set of aligned sequences

                  Ivo L Hofacker
                  Vienna RNA package
*/

/**
*** \file alipfold.c
**/

#include <config.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>    /* #defines FLT_MIN */
#include <limits.h>

#include "utils.h"
#include "energy_par.h"
#include "fold_vars.h"
#include "pair_mat.h"
#include "PS_dot.h"
#include "ribo.h"
#include "params.h"
#include "loop_energies.h"
#include "part_func.h"
#include "alifold.h"

#ifdef _OPENMP
#include <omp.h>
#endif

/*@unused@*/
static char rcsid[] = "$Id: alipfold.c,v 1.17 2009/02/24 14:21:33 ivo Exp $";

#define STACK_BULGE1  1   /* stacking energies for bulges of size 1 */
#define NEW_NINIO     1   /* new asymetry penalty */
#define ISOLATED  256.0
#define UNIT 100
#define MINPSCORE -2 * UNIT
#define PMIN 0.0008

/*
#################################
# PUBLIC GLOBAL VARIABLES       #
#################################
*/

/*
#################################
# PRIVATE GLOBAL VARIABLES      #
#################################
*/
PRIVATE FLT_OR_DBL      *expMLbase=NULL;
PRIVATE FLT_OR_DBL      *q=NULL, *qb=NULL, *qm=NULL, *qm1=NULL, *qqm=NULL, *qqm1=NULL, *qq=NULL, *qq1=NULL;
PRIVATE FLT_OR_DBL      *prml=NULL, *prm_l=NULL, *prm_l1=NULL, *q1k=NULL, *qln=NULL;
PRIVATE FLT_OR_DBL      *probs=NULL;
PRIVATE FLT_OR_DBL      *scale=NULL;
PRIVATE short           *pscore=NULL;   /* precomputed array of covariance bonus/malus */
/* some additional things for circfold  */
PRIVATE int             circular=0;
PRIVATE FLT_OR_DBL      qo, qho, qio, qmo, *qm2=NULL;
PRIVATE int             *jindx=NULL;
PRIVATE int             *my_iindx=NULL;
PRIVATE int             do_bppm = 1;             /* do backtracking per default */
PRIVATE short           **S=NULL;
PRIVATE short           **S5=NULL;               /*S5[s][i] holds next base 5' of i in sequence s*/
PRIVATE short           **S3=NULL;               /*S3[s][i] holds next base 3' of i in sequence s*/
PRIVATE char            **Ss=NULL;
PRIVATE unsigned short  **a2s=NULL;
PRIVATE int             N_seq = 0;
PRIVATE pf_paramT       *pf_params = NULL;
PRIVATE char            *pstruc=NULL;
PRIVATE double          alpha = 1.0;
PRIVATE int             struct_constrained = 0;

#ifdef _OPENMP

/* NOTE: all variables are assumed to be uninitialized if they are declared as threadprivate
         thus we have to initialize them before usage by a seperate function!
         OR: use copyin in the PARALLEL directive!
         e.g.:
         #pragma omp parallel for copyin(pf_params)
*/
#pragma omp threadprivate(expMLbase, q, qb, qm, qm1, qqm, qqm1, qq, qq1,\
                          probs, prml, prm_l, prm_l1, q1k, qln,\
                          scale, pscore, circular,\
                          qo, qho, qio, qmo, qm2, jindx, my_iindx,\
                          S, S5, S3, Ss, a2s, N_seq, pf_params, pstruc, alpha, struct_constrained)

#endif

/*
#################################
# PRIVATE FUNCTION DECLARATIONS #
#################################
*/

PRIVATE void      init_alipf_fold(int length, int n_seq, pf_paramT *parameters);
PRIVATE void      scale_pf_params(unsigned int length, int n_seq, pf_paramT *parameters);
PRIVATE void      get_arrays(unsigned int length);
PRIVATE void      make_pscores(const short *const *S, const char **AS, int n_seq, const char *structure);
PRIVATE pair_info *make_pairinfo(const short *const* S, const char **AS, int n_seq);
PRIVATE void      alipf_circ(const char **sequences, char *structure);
PRIVATE void      alipf_linear(const char **sequences, char *structure);
PRIVATE void      alipf_create_bppm(const char **sequences, char *structure, struct plist **pl);
PRIVATE void      backtrack(int i, int j, int n_seq, double *prob);
PRIVATE void      backtrack_qm1(int i,int j, int n_seq, double *prob);



/*
#################################
# BEGIN OF FUNCTION DEFINITIONS #
#################################
*/

PRIVATE void init_alipf_fold(int length, int n_seq, pf_paramT *parameters){
  if (length<1) nrerror("init_alipf_fold: length must be greater 0");

#ifdef _OPENMP
/* Explicitly turn off dynamic threads */
  omp_set_dynamic(0);
#endif

#ifdef SUN4
  nonstandard_arithmetic();
#else
#ifdef HP9
  fpsetfastmode(1);
#endif
#endif
  make_pair_matrix();
  free_alipf_arrays(); /* free previous allocation */
  get_arrays((unsigned) length);
  scale_pf_params((unsigned) length, n_seq, parameters);
}

/**
*** Allocate memory for all matrices and other stuff
**/
PRIVATE void get_arrays(unsigned int length){
  unsigned int size,i;

  if((length+1) >= (unsigned int)sqrt((double)INT_MAX))
    nrerror("get_arrays@alipfold.c: sequence length exceeds addressable range");

  size  = sizeof(FLT_OR_DBL) * ((length+1)*(length+2)/2);

  q     = (FLT_OR_DBL *) space(size);
  qb    = (FLT_OR_DBL *) space(size);
  qm    = (FLT_OR_DBL *) space(size);
  qm1   = (FLT_OR_DBL *) space(size);
  qm2   = (circular) ? (FLT_OR_DBL *) space(sizeof(FLT_OR_DBL)*(length+2)) : NULL;

  pscore    = (short *) space(sizeof(short)*((length+1)*(length+2)/2));
  q1k       = (FLT_OR_DBL *) space(sizeof(FLT_OR_DBL)*(length+1));
  qln       = (FLT_OR_DBL *) space(sizeof(FLT_OR_DBL)*(length+2));
  qq        = (FLT_OR_DBL *) space(sizeof(FLT_OR_DBL)*(length+2));
  qq1       = (FLT_OR_DBL *) space(sizeof(FLT_OR_DBL)*(length+2));
  qqm       = (FLT_OR_DBL *) space(sizeof(FLT_OR_DBL)*(length+2));
  qqm1      = (FLT_OR_DBL *) space(sizeof(FLT_OR_DBL)*(length+2));
  prm_l     = (FLT_OR_DBL *) space(sizeof(FLT_OR_DBL)*(length+2));
  prm_l1    = (FLT_OR_DBL *) space(sizeof(FLT_OR_DBL)*(length+2));
  prml      = (FLT_OR_DBL *) space(sizeof(FLT_OR_DBL)*(length+2));
  expMLbase = (FLT_OR_DBL *) space(sizeof(FLT_OR_DBL)*(length+1));
  scale     = (FLT_OR_DBL *) space(sizeof(FLT_OR_DBL)*(length+1));

  my_iindx  = get_iindx(length);
  iindx     = get_iindx(length); /* for backward compatibility and Perl wrapper */
  jindx     = get_indx(length);
}

/*----------------------------------------------------------------------*/


PUBLIC void free_alipf_arrays(void){
  if(q)         free(q);
  if(qb)        free(qb);
  if(qm)        free(qm);
  if(qm1)       free(qm1);
  if(qm2)       free(qm2);
  if(pscore)    free(pscore);
  if(qq)        free(qq);
  if(qq1)       free(qq1);
  if(qqm)       free(qqm);
  if(qqm1)      free(qqm1);
  if(q1k)       free(q1k);
  if(qln)       free(qln);
  if(prm_l)     free(prm_l);
  if(prm_l1)    free(prm_l1);
  if(prml)      free(prml);
  if(expMLbase) free(expMLbase);
  if(scale)     free(scale);
  if(my_iindx)  free(my_iindx);
  if(iindx)     free(iindx); /* for backward compatibility and Perl wrapper */
  if(jindx)     free(jindx);

  if(S){
    free_sequence_arrays(N_seq, &S, &S5, &S3, &a2s, &Ss);
    N_seq = 0;
    S = NULL;
  }
  pr = NULL; /* ? */
  q = probs = qb = qm = qm1 = qm2 = qq = qq1 = qqm = qqm1 = q1k = qln = prml = prm_l = prm_l1 = expMLbase = scale = NULL;
  my_iindx   = jindx = iindx = NULL;
  pscore  = NULL;

#ifdef SUN4
  standard_arithmetic();
#else
#ifdef HP9
  fpsetfastmode(0);
#endif
#endif
}

/*-----------------------------------------------------------------*/
PUBLIC float alipf_fold(const char **sequences, char *structure, plist **pl){
  return alipf_fold_par(sequences, structure, pl, NULL, do_backtrack, fold_constrained, 0);
}

PUBLIC float alipf_circ_fold(const char **sequences, char *structure, plist **pl){
  return alipf_fold_par(sequences, structure, pl, NULL, do_backtrack, fold_constrained, 1);
}

PUBLIC float alipf_fold_par(const char **sequences,
                            char *structure,
                            plist **pl,
                            pf_paramT *parameters,
                            int calculate_bppm,
                            int is_constrained,
                            int is_circular){

  int         n, s, n_seq;
  FLT_OR_DBL  Q;
  float       free_energy;

  circular            = is_circular;
  do_bppm             = calculate_bppm;
  struct_constrained  = is_constrained;

  if(circular)
    oldAliEn  = 1; /* may be removed if circular alipf fold works with gapfree stuff */

  n = (int) strlen(sequences[0]);
  for (s=0; sequences[s]!=NULL; s++);
  n_seq = N_seq = s;

  init_alipf_fold(n, n_seq, parameters);

  alloc_sequence_arrays(sequences, &S, &S5, &S3, &a2s, &Ss, circular);
  make_pscores((const short *const*)S, sequences, n_seq, structure);

  alipf_linear(sequences, structure);

  /* calculate post processing step for circular  */
  /* RNAs                                         */
  if(circular)
    alipf_circ(sequences, structure);

  if (backtrack_type=='C')      Q = qb[my_iindx[1]-n];
  else if (backtrack_type=='M') Q = qm[my_iindx[1]-n];
  else Q = (circular) ? qo : q[my_iindx[1]-n];

  /* ensemble free energy in Kcal/mol */
  if (Q<=FLT_MIN) fprintf(stderr, "pf_scale too large\n");
  free_energy = (-log(Q)-n*log(pf_params->pf_scale))*pf_params->kT/(1000.0 * n_seq);
  /* in case we abort because of floating point errors */
  if (n>1600) fprintf(stderr, "free energy = %8.2f\n", free_energy);

  /* backtracking to construct binding probabilities of pairs*/
  if(do_bppm){
    alipf_create_bppm(sequences, structure, pl);
    /*
    *  Backward compatibility:
    *  This block may be removed if deprecated functions
    *  relying on the global variable "pr" vanish from within the package!
    */
    pr = probs;
  }

  return free_energy;
}

PUBLIC FLT_OR_DBL *alipf_export_bppm(void){
  return probs;
}



PRIVATE void alipf_linear(const char **sequences, char *structure)
{
  int         s, n, n_seq, i,j,k,l, ij, u, u1, d, ii, *type, type_2, tt;
  FLT_OR_DBL  temp, Qmax=0;
  FLT_OR_DBL  qbt1, *tmp;
  double      kTn;

  FLT_OR_DBL  expMLclosing  = pf_params->expMLclosing;

  for(s=0; sequences[s]!=NULL; s++);

  n_seq = s;
  n     = (int) strlen(sequences[0]);
  kTn   = pf_params->kT/10.;   /* kT in cal/mol  */
  type  = (int *)space(sizeof(int) * n_seq);

  /* array initialization ; qb,qm,q
     qb,qm,q (i,j) are stored as ((n+1-i)*(n-i) div 2 + n+1-j */

  for (d=0; d<=TURN; d++)
    for (i=1; i<=n-d; i++) {
      j=i+d;
      ij = my_iindx[i]-j;
      q[ij]=1.0*scale[d+1];
      qb[ij]=qm[ij]=0.0;
    }

  for (i=1; i<=n; i++)
    qq[i]=qq1[i]=qqm[i]=qqm1[i]=0;

  for (j=TURN+2;j<=n; j++) {
    for (i=j-TURN-1; i>=1; i--) {
      int ij, psc;
      /* construction of partition function for segment i,j */
      /* calculate pf given that i and j pair: qb(i,j)      */
       ij = my_iindx[i]-j;

      for (s=0; s<n_seq; s++) {
        type[s] = pair[S[s][i]][S[s][j]];
        if (type[s]==0) type[s]=7;
      }
      psc = pscore[ij];
      if (psc>=cv_fact*MINPSCORE) {   /* otherwise ignore this pair */

        /* hairpin contribution */
        for (qbt1=1,s=0; s<n_seq; s++) {
          u = a2s[s][j-1]-a2s[s][i];
          if (a2s[s][i]<1) continue;
          char loopseq[10];
          if (u<7){
            strncpy(loopseq, Ss[s]+a2s[s][i]-1, 10);
          }
          qbt1 *= exp_E_Hairpin(u, type[s], S3[s][i], S5[s][j], loopseq, pf_params);

        }
        qbt1 *= scale[j-i+1];

        /* interior loops with interior pair k,l */
        for (k=i+1; k<=MIN2(i+MAXLOOP+1,j-TURN-2); k++){

          for (l=MAX2(k+TURN+1,j-1-MAXLOOP+k-i-1); l<=j-1; l++){
            double qloop=1;
            if (qb[my_iindx[k]-l]==0) {qloop=0; continue;}
            for (s=0; s<n_seq; s++) {
              u1 = a2s[s][k-1]-a2s[s][i]/*??*/;
              type_2 = pair[S[s][l]][S[s][k]]; if (type_2 == 0) type_2 = 7;
              qloop *= exp_E_IntLoop( u1, a2s[s][j-1]-a2s[s][l],
                                  type[s], type_2, S3[s][i],
                                  S5[s][j], S5[s][k], S3[s][l],
                                  pf_params
                                );
            }
            qbt1 += qb[my_iindx[k]-l] * qloop * scale[k-i+j-l];
          }
        }

        /* multi-loop loop contribution */
        ii = my_iindx[i+1]; /* ii-k=[i+1,k-1] */
        temp = 0.0;
        for (k=i+2; k<=j-1; k++) temp += qm[ii-(k-1)]*qqm1[k];
        for (s=0; s<n_seq; s++) {
          tt = rtype[type[s]];
          temp *= exp_E_MLstem(tt, S5[s][j], S3[s][i], pf_params)* expMLclosing;
        }
        temp *= scale[2] ;
        qbt1 += temp;
        qb[ij] = qbt1;
        qb[ij] *= exp(psc/kTn);
      } /* end if (type!=0) */
      else qb[ij] = 0.0;
      /* construction of qqm matrix containing final stem
         contributions to multiple loop partition function
         from segment i,j */
      qqm[i] = qqm1[i]*expMLbase[1];  /* expMLbase[1]^n_seq */
      for (qbt1=1, s=0; s<n_seq; s++) {
        qbt1 *= exp_E_MLstem(type[s], (i>1) || circular ? S5[s][i] : -1, (j<n) || circular ? S3[s][j] : -1, pf_params);
      }
      qqm[i] += qb[ij]*qbt1;
      qm1[jindx[j]+i] = qqm[i]; /* for circ folding and stochBT */

      /* construction of qm matrix containing multiple loop
         partition function contributions from segment i,j */
      temp = 0.0;
      ii = my_iindx[i];  /* ii-k=[i,k-1] */
      for (k=i+1; k<=j; k++)
        temp += (qm[ii-(k-1)]+expMLbase[k-i])*qqm[k];
      qm[ij] = (temp + qqm[i]);

      /* auxiliary matrix qq for cubic order q calculation below */
      qbt1 = qb[ij];
      if (qbt1>0)
        for (s=0; s<n_seq; s++) {
          qbt1 *= exp_E_ExtLoop(type[s], (i>1) || circular ? S5[s][i] : -1, (j<n) || circular ? S3[s][j] : -1, pf_params);
        }
      qq[i] = qq1[i]*scale[1] + qbt1;

      /* construction of partition function for segment i,j */
      temp = 1.0*scale[1+j-i] + qq[i];
      for (k=i; k<=j-1; k++) temp += q[ii-k]*qq[k+1];
      q[ij] = temp;

#ifndef LARGE_PF
      if (temp>Qmax) {
        Qmax = temp;
        if (Qmax>FLT_MAX/10.)
          fprintf(stderr, "%d %d %g\n", i,j,temp);
      }
      if (temp>FLT_MAX) {
        PRIVATE char msg[128];
        sprintf(msg, "overflow in pf_fold while calculating q[%d,%d]\n"
          "use larger pf_scale", i,j);
        nrerror(msg);
      }
#endif
    }
    tmp = qq1;  qq1 =qq;  qq =tmp;
    tmp = qqm1; qqm1=qqm; qqm=tmp;
  }

  free(type);
}

PRIVATE void alipf_create_bppm(const char **sequences, char *structure, plist **pl)
{
  int s;
  int n, n_seq, i,j,k,l, ij, kl, ii, ll, tt, *type, ov=0;
  FLT_OR_DBL temp, Qmax=0, prm_MLb;
  FLT_OR_DBL prmt,prmt1;
  FLT_OR_DBL qbt1, *tmp, tmp2, tmp3;
  FLT_OR_DBL  expMLclosing      = pf_params->expMLclosing;

  double kTn;
  n = (int) strlen(sequences[0]);
  for (s=0; sequences[s]!=NULL; s++);
  n_seq = s;
  type  = (int *)space(sizeof(int) * n_seq);

  kTn = pf_params->kT/10.;   /* kT in cal/mol  */

  for (i=0; i<=n; i++)
    prm_l[i]=prm_l1[i]=prml[i]=0;

  /* backtracking to construct binding probabilities of pairs*/
  Qmax=0;

  for (k=1; k<=n; k++) {
    q1k[k] = q[my_iindx[1] - k];
    qln[k] = q[my_iindx[k] -n];
  }
  q1k[0] = 1.0;
  qln[n+1] = 1.0;

  probs = q;     /* recycling */

  /* 1. exterior pair i,j and initialization of pr array */
  if(circular){
    for (i=1; i<=n; i++) {
      for (j=i; j<=MIN2(i+TURN,n); j++) probs[my_iindx[i]-j] = 0;
      for (j=i+TURN+1; j<=n; j++) {
        ij = my_iindx[i]-j;
        if (qb[ij]>0.) {
          probs[ij] =  exp(pscore[ij]/kTn)/qo;

          /* get pair types  */
          for (s=0; s<n_seq; s++) {
            type[s] = pair[S[s][j]][S[s][i]];
            if (type[s]==0) type[s]=7;
          }
          int rt;

          /* 1.1. Exterior Hairpin Contribution */
          int u = i + n - j -1;
          for (qbt1=1.,s=0; s<n_seq; s++) {

            char loopseq[10];
            char *ts;
            if (u<7){
              strcpy(loopseq , sequences[s]+j-1);
              strncat(loopseq, sequences[s], i);
            }
            qbt1 *= exp_E_Hairpin(u, type[s], S[s][j+1], S[s][(i>1) ? i-1 : n], loopseq, pf_params);
          }
          tmp2 = qbt1 * scale[u];

          /* 1.2. Exterior Interior Loop Contribution */
          /* recycling of k and l... */
          /* 1.2.1. first we calc exterior loop energy with constraint, that i,j  */
          /* delimtis the "left" part of the interior loop                        */
          /* (j,i) is "outer pair"                                                */
          for(k=1; k < i-TURN-1; k++){
            /* so first, lets calc the length of loop between j and k */
            int ln1, lstart;
            ln1 = k + n - j - 1;
            if(ln1>MAXLOOP) break;
            lstart = ln1+i-1-MAXLOOP;
            if(lstart<k+TURN+1) lstart = k + TURN + 1;
            for(l=lstart; l < i; l++){
              int ln2,ln2a,ln1a, type_2;
              ln2 = i - l - 1;
                if(ln1+ln2>MAXLOOP) continue;

              double qloop=1.;
              if (qb[my_iindx[k]-l]==0.){ qloop=0.; continue;}

              for (s=0; s<n_seq; s++){
                ln2a= a2s[s][i-1];
                ln2a-=a2s[s][l];
                ln1a= a2s[s][n]-a2s[s][j];
                ln1a+=a2s[s][k-1];
                type_2 = pair[S[s][l]][S[s][k]];
                if (type_2 == 0) type_2 = 7;
                qloop *= exp_E_IntLoop(ln1a, ln2a, type[s], type_2,
                            S[s][j+1],
                            S[s][i-1],
                            S[s][(k>1) ? k-1 : n],
                            S[s][l+1], pf_params);
              }
              tmp2 += qb[my_iindx[k] - l] * qloop * scale[ln1+ln2];
            }
          }

          /* 1.2.2. second we calc exterior loop energy with constraint, that i,j  */
          /* delimtis the "right" part of the interior loop                        */
          /* (l,k) is "outer pair"                                                */
          for(k=j+1; k < n-TURN; k++){
            /* so first, lets calc the length of loop between l and i */
            int ln1, lstart;
            ln1 = k - j - 1;
            if((ln1 + i - 1)>MAXLOOP) break;
            lstart = ln1+i-1+n-MAXLOOP;
            if(lstart<k+TURN+1) lstart = k + TURN + 1;
            for(l=lstart; l <= n; l++){
              int ln2, type_2;
              ln2 = i - 1 + n - l;
              if(ln1+ln2>MAXLOOP) continue;
              double qloop=1.;
              if (qb[my_iindx[k]-l]==0.){ qloop=0.; continue;}

              for (s=0; s<n_seq; s++){
                    ln1 = a2s[s][k] - a2s[s][j+1];
                    ln2 = a2s[s][i-1] + a2s[s][n] - a2s[s][l];
                    type_2 = pair[S[s][l]][S[s][k]];
                    if (type_2 == 0) type_2 = 7;
                    qloop *= exp_E_IntLoop(ln2, ln1, type_2, type[s],
                            S3[s][l],
                            S5[s][k],
                            S5[s][i],
                            S3[s][j], pf_params);
              }
              tmp2 += qb[my_iindx[k] - l] * qloop * scale[(k-j-1)+(i-1+n-l)];
            }
          }

          /* 1.3 Exterior multiloop decomposition */
          /* 1.3.1 Middle part                    */
          if((i>TURN+2) && (j<n-TURN-1)){

            for (tmp3=1, s=0; s<n_seq; s++){
              tmp3 *= exp_E_MLstem(rtype[type[s]], S5[s][i], S3[s][j], pf_params);
            }
            tmp2 += qm[my_iindx[1]-i+1] * qm[my_iindx[j+1]-n] * tmp3 * pow(expMLclosing,n_seq);
          }
          /* 1.3.2 Left part    */
          for(k=TURN+2; k < i-TURN-2; k++){

            for (tmp3=1, s=0; s<n_seq; s++){
              tmp3 *= exp_E_MLstem(rtype[type[s]], S5[s][i], S3[s][j], pf_params);
            }
            tmp2 += qm[my_iindx[1]-k] * qm1[jindx[i-1]+k+1] * tmp3 * expMLbase[n-j] * pow(expMLclosing,n_seq);
          }
          /* 1.3.3 Right part    */
          for(k=j+TURN+2; k < n-TURN-1;k++){

            for (tmp3=1, s=0; s<n_seq; s++){
              tmp3 *= exp_E_MLstem(rtype[type[s]], S5[s][i], S3[s][j], pf_params);
            }
            tmp2 += qm[my_iindx[j+1]-k] * qm1[jindx[n]+k+1] * tmp3 * expMLbase[i-1] * pow(expMLclosing,n_seq);
          }
          probs[ij] *= tmp2;
        }
        else probs[ij] = 0;
      }  /* end for j=..*/
    }  /* end or i=...  */
  } /* end if(circular)  */
  else{
    for (i=1; i<=n; i++) {
      for (j=i; j<=MIN2(i+TURN,n); j++) probs[my_iindx[i]-j] = 0;
      for (j=i+TURN+1; j<=n; j++) {
        ij = my_iindx[i]-j;
        if (qb[ij]>0.){
          probs[ij] = q1k[i-1]*qln[j+1]/q1k[n] * exp(pscore[ij]/kTn);
          for (s=0; s<n_seq; s++) {
            int typ;
            typ = pair[S[s][i]][S[s][j]]; if (typ==0) typ=7;
            probs[ij] *= exp_E_ExtLoop(typ, (i>1) ? S5[s][i] : -1, (j<n) ? S3[s][j] : -1, pf_params);
          }
        } else
          probs[ij] = 0;
      }
    }
  } /* end if(!circular)  */
  for (l=n; l>TURN+1; l--) {

    /* 2. bonding k,l as substem of 2:loop enclosed by i,j */
    for (k=1; k<l-TURN; k++) {
      double pp = 0;
      kl = my_iindx[k]-l;
      if (qb[kl]==0) continue;
      for (s=0; s<n_seq; s++) {
            type[s] = pair[S[s][l]][S[s][k]];
            if (type[s]==0) type[s]=7;
      }

      for (i=MAX2(1,k-MAXLOOP-1); i<=k-1; i++)
        for (j=l+1; j<=MIN2(l+ MAXLOOP -k+i+2,n); j++) {
          ij = my_iindx[i] - j;
          if ((probs[ij]>0.)) {
            double qloop=1;
            for (s=0; s<n_seq; s++) {
              int typ;
              typ = pair[S[s][i]][S[s][j]]; if (typ==0) typ=7;
              qloop *=  exp_E_IntLoop(a2s[s][k-1]-a2s[s][i], a2s[s][j-1]-a2s[s][l], typ, type[s], S3[s][i], S5[s][j], S5[s][k], S3[s][l], pf_params);
            }
            pp += probs[ij]*qloop*scale[k-i + j-l];
          }
        }
      probs[kl] += pp * exp(pscore[kl]/kTn);
    }
    /* 3. bonding k,l as substem of multi-loop enclosed by i,j */
    prm_MLb = 0.;
    if (l<n) for (k=2; k<l-TURN; k++) {
      i = k-1;
      prmt = prmt1 = 0.0;

      ii = my_iindx[i];     /* ii-j=[i,j]     */
      ll = my_iindx[l+1];   /* ll-j=[l+1,j-1] */
      prmt1 = probs[ii-(l+1)];
      for (s=0; s<n_seq; s++) {
        tt = pair[S[s][l+1]][S[s][i]]; if (tt==0) tt=7;
        prmt1 *= exp_E_MLstem(tt, S5[s][l+1], S3[s][i], pf_params) * expMLclosing;
      }

      for (j=l+2; j<=n; j++) {
        double pp=1;
        if (probs[ii-j]==0) continue;
        for (s=0; s<n_seq; s++) {
          tt=pair[S[s][j]][S[s][i]]; if (tt==0) tt=7;
          pp *=  exp_E_MLstem(tt, S5[s][j], S3[s][i], pf_params)*expMLclosing;
        }
        prmt +=  probs[ii-j]*pp*qm[ll-(j-1)];
      }
      kl = my_iindx[k]-l;

      prml[ i] = prmt;
      prm_l[i] = prm_l1[i]*expMLbase[1]+prmt1; /* expMLbase[1]^n_seq */

      prm_MLb = prm_MLb*expMLbase[1] + prml[i];
      /* same as:    prm_MLb = 0;
         for (i=1; i<=k-1; i++) prm_MLb += prml[i]*expMLbase[k-i-1]; */

      prml[i] = prml[ i] + prm_l[i];

      if (qb[kl] == 0.) continue;

      temp = prm_MLb;

      for (i=1;i<=k-2; i++)
        temp += prml[i]*qm[my_iindx[i+1] - (k-1)];

      for (s=0; s<n_seq; s++) {
        tt=pair[S[s][k]][S[s][l]]; if (tt==0) tt=7;
        temp *= exp_E_MLstem(tt, S5[s][k], S3[s][l], pf_params);
      }
      probs[kl] += temp * scale[2] * exp(pscore[kl]/kTn);


#ifndef LARGE_PF
      if (probs[kl]>Qmax) {
        Qmax = probs[kl];
        if (Qmax>FLT_MAX/10.)
          fprintf(stderr, "%d %d %g %g\n", i,j,probs[kl],qb[kl]);
      }
      if (probs[kl]>FLT_MAX) {
        ov++;
        probs[kl]=FLT_MAX;
      }
#endif
    } /* end for (k=2..) */
    tmp = prm_l1; prm_l1=prm_l; prm_l=tmp;

  }  /* end for (l=..)   */

  for (i=1; i<=n; i++)
    for (j=i+TURN+1; j<=n; j++) {
      ij = my_iindx[i]-j;
      probs[ij] *= qb[ij] *exp(-pscore[ij]/kTn);
    }

  /* did we get an adress where to save a pair-list? */
  if (pl != NULL)
    assign_plist_from_pr(pl, probs, n, /*cut_off:*/ 1e-6);

  if (structure!=NULL)
    bppm_to_structure(structure, probs, n);

  if (ov>0) fprintf(stderr, "%d overflows occurred while backtracking;\n"
        "you might try a smaller pf_scale than %g\n",
        ov, pf_params->pf_scale);

  free(type);
}

PRIVATE void scale_pf_params(unsigned int length, int n_seq, pf_paramT *parameters){
  unsigned int i;
  double  kT, scaling_factor;

  if(pf_params) free(pf_params);

  if(parameters){
    pf_params = get_boltzmann_factor_copy(parameters);
  } else {
    model_detailsT  md;
    set_model_details(&md);
    pf_params = get_boltzmann_factors_ali(n_seq, temperature, alpha, md, pf_scale);
  }

  scaling_factor  = pf_params->pf_scale;
  kT              = pf_params->kT / n_seq;

  /* scaling factors (to avoid overflows) */
  if (scaling_factor == -1) { /* mean energy for random sequences: 184.3*length cal */
    scaling_factor = exp(-(-185+(pf_params->temperature-37.)*7.27)/kT);
    if (scaling_factor<1) scaling_factor=1;
    pf_params->pf_scale = scaling_factor;
  }
  scale[0] = 1.;
  scale[1] = 1./scaling_factor;

  expMLbase[0] = 1;
  expMLbase[1] = pf_params->expMLbase/scaling_factor;
  for (i=2; i<=length; i++) {
    scale[i] = scale[i/2]*scale[i-(i/2)];
    expMLbase[i] = pow(pf_params->expMLbase, (double)i) * scale[i];
  }
}


/*---------------------------------------------------------------------------*/
PRIVATE int compare_pair_info(const void *pi1, const void *pi2) {
  pair_info *p1, *p2;
  int  i, nc1, nc2;
  p1 = (pair_info *)pi1;  p2 = (pair_info *)pi2;
  for (nc1=nc2=0, i=1; i<=6; i++) {
    if (p1->bp[i]>0) nc1++;
    if (p2->bp[i]>0) nc2++;
  }
  /* sort mostly by probability, add
     epsilon * comp_mutations/(non-compatible+1) to break ties */
  return (p1->p + 0.01*nc1/(p1->bp[0]+1.)) <
         (p2->p + 0.01*nc2/(p2->bp[0]+1.)) ? 1 : -1;
}

pair_info *make_pairinfo(const short *const* S, const char **AS, int n_seq) {
  int i,j,n, num_p=0, max_p = 64;
  pair_info *pi;
  double *duck, p;
  n = S[0][0];
  max_p = 64; pi = space(max_p*sizeof(pair_info));
  duck =  (double *) space((n+1)*sizeof(double));
  for (i=1; i<n; i++)
    for (j=i+TURN+1; j<=n; j++)
      if ((p=probs[my_iindx[i]-j])>0) {
        duck[i] -=  p * log(p);
        duck[j] -=  p * log(p);
      }

  for (i=1; i<n; i++)
    for (j=i+TURN+1; j<=n; j++) {
      if ((p=probs[my_iindx[i]-j])>=PMIN) {
        int type, s;
        pi[num_p].i = i;
        pi[num_p].j = j;
        pi[num_p].p = p;
        pi[num_p].ent =  duck[i]+duck[j]-p*log(p);
        for (type=0; type<8; type++) pi[num_p].bp[type]=0;
        for (s=0; s<n_seq; s++) {
          if (S[s][i]==0 && S[s][j]==0) type = 7; /* gap-gap  */
          else {
            if ((AS[s][i] == '~')||(AS[s][j] == '~')) type = 7;
            else type = pair[S[s][i]][S[s][j]];
          }
          pi[num_p].bp[type]++;
        }
        num_p++;
        if (num_p>=max_p) {
          max_p *= 2;
          pi = xrealloc(pi, max_p * sizeof(pair_info));
        }
      }
    }
  free(duck);
  pi = xrealloc(pi, (num_p+1)*sizeof(pair_info));
  pi[num_p].i=0;
  qsort(pi, num_p, sizeof(pair_info), compare_pair_info );
  return pi;
}

/*---------------------------------------------------------------------------*/
PRIVATE void make_pscores(const short *const* S, const char **AS,
                          int n_seq, const char *structure) {
  /* calculate co-variance bonus for each pair depending on  */
  /* compensatory/consistent mutations and incompatible seqs */
  /* should be 0 for conserved pairs, >0 for good pairs      */
#define NONE -10000 /* score for forbidden pairs */
  int n,i,j,k,l,s,score;

  int olddm[7][7]={{0,0,0,0,0,0,0}, /* hamming distance between pairs */
                  {0,0,2,2,1,2,2} /* CG */,
                  {0,2,0,1,2,2,2} /* GC */,
                  {0,2,1,0,2,1,2} /* GU */,
                  {0,1,2,2,0,2,1} /* UG */,
                  {0,2,2,1,2,0,2} /* AU */,
                  {0,2,2,2,1,2,0} /* UA */};

  float **dm;
  int noLP = pf_params->model_details.noLP;

  n=S[0][0];  /* length of seqs */
  if (ribo) {
    if (RibosumFile !=NULL) dm=readribosum(RibosumFile);
    else dm=get_ribosum(AS,n_seq,n);
  }
  else { /*use usual matrix*/
    dm=(float **)space(7*sizeof(float*));
    for (i=0; i<7;i++) {
      dm[i]=(float *)space(7*sizeof(float));
      for (j=0; j<7; j++)
        dm[i][j] = olddm[i][j];
    }
  }

  n=S[0][0];  /* length of seqs */
  for (i=1; i<n; i++) {
    for (j=i+1; (j<i+TURN+1) && (j<=n); j++)
      pscore[my_iindx[i]-j] = NONE;
    for (j=i+TURN+1; j<=n; j++) {
      int pfreq[8]={0,0,0,0,0,0,0,0};
      for (s=0; s<n_seq; s++) {
        int type;
        if (S[s][i]==0 && S[s][j]==0) type = 7; /* gap-gap  */
        else {
          if ((AS[s][i] == '~')||(AS[s][j] == '~')) type = 7;
          else type = pair[S[s][i]][S[s][j]];
        }
        pfreq[type]++;
      }
      if (pfreq[0]*2+pfreq[7]>n_seq) { pscore[my_iindx[i]-j] = NONE; continue;}
      for (k=1,score=0; k<=6; k++) /* ignore pairtype 7 (gap-gap) */
        for (l=k; l<=6; l++)
          /* scores for replacements between pairtypes    */
          /* consistent or compensatory mutations score 1 or 2  */
          score += pfreq[k]*pfreq[l]*dm[k][l];
      /* counter examples score -1, gap-gap scores -0.25  */
      pscore[my_iindx[i]-j] = cv_fact *
        ((UNIT*score)/n_seq - nc_fact*UNIT*(pfreq[0] + pfreq[7]*0.25));
    }
  }

  if (noLP) /* remove unwanted pairs */
    for (k=1; k<=n-TURN-1; k++)
      for (l=1; l<=2; l++) {
        int type,ntype=0,otype=0;
        i=k; j = i+TURN+l;
        type = pscore[my_iindx[i]-j];
        while ((i>=1)&&(j<=n)) {
          if ((i>1)&&(j<n)) ntype = pscore[my_iindx[i-1]-j-1];
          if ((otype<cv_fact*MINPSCORE)&&(ntype<cv_fact*MINPSCORE))
            /* too many counterexamples */
            pscore[my_iindx[i]-j] = NONE; /* i.j can only form isolated pairs */
          otype =  type;
          type  = ntype;
          i--; j++;
        }
      }


  if (struct_constrained&&(structure!=NULL)) {
    int psij, hx, *stack;
    stack = (int *) space(sizeof(int)*(n+1));

    for(hx=0, j=1; j<=n; j++) {
      switch (structure[j-1]) {
      case 'x': /* j can't pair */
        for (l=1; l<j-TURN; l++) pscore[my_iindx[l]-j] = NONE;
        for (l=j+TURN+1; l<=n; l++) pscore[my_iindx[j]-l] = NONE;
        break;
      case '(':
        stack[hx++]=j;
        /* fallthrough */
      case '<': /* j pairs upstream */
        for (l=1; l<j-TURN; l++) pscore[my_iindx[l]-j] = NONE;
        break;
      case ')': /* j pairs with i */
        if (hx<=0) {
          fprintf(stderr, "%s\n", structure);
          nrerror("unbalanced brackets in constraints");
        }
        i = stack[--hx];
        psij = pscore[my_iindx[i]-j]; /* store for later */
        for (l=i; l<=j; l++)
          for (k=j; k<=n; k++) pscore[my_iindx[l]-k] = NONE;
        for (k=1; k<=i; k++)
          for (l=i; l<=j; l++) pscore[my_iindx[k]-l] = NONE;
        for (k=i+1; k<j; k++)
          pscore[my_iindx[k]-j] = pscore[my_iindx[i]-k] = NONE;
        pscore[my_iindx[i]-j] = (psij>0) ? psij : 0;
        /* fallthrough */
      case '>': /* j pairs downstream */
        for (l=j+TURN+1; l<=n; l++) pscore[my_iindx[j]-l] = NONE;
        break;
      }
    }
    if (hx!=0) {
      fprintf(stderr, "%s\n", structure);
      nrerror("unbalanced brackets in constraint string");
    }
    free(stack);
  }
  for (i=0; i<7;i++) {
    free(dm[i]);
  }
  free(dm);
}

/* calculate partition function for circular case   */
/* NOTE: this is the postprocessing step ONLY        */
/* You have to call alipf_linear first to calculate  */
/* circular case!!!                                  */

PUBLIC void alipf_circ(const char **sequences, char *structure){

  int u, p, q, k, l, n_seq, s, *type;
  FLT_OR_DBL  expMLclosing      = pf_params->expMLclosing;

  int n = (int) strlen(sequences[0]);
  for (s=0; sequences[s]!=NULL; s++);
  n_seq = s;

  double kTn;
  FLT_OR_DBL qbt1, qot;
  kTn = pf_params->kT/10.;   /* kT in cal/mol  */
  type  = (int *)space(sizeof(int) * n_seq);

  qo = qho = qio = qmo = 0.;
  /* calculate the qm2 matrix  */
  for(k=1; k<n-TURN; k++){
    qot = 0.;
    for (u=k+TURN+1; u<n-TURN-1; u++)
      qot += qm1[jindx[u]+k]*qm1[jindx[n]+(u+1)];
    qm2[k] = qot;
  }

  for(p=1;p<n;p++){
    for(q=p+TURN+1;q<=n;q++){
      int psc;
      u = n-q + p-1;
      if (u<TURN) continue;

      psc = pscore[my_iindx[p]-q];

      if(psc<cv_fact*MINPSCORE) continue;

      /* 1. exterior hairpin contribution  */
      /* Note, that we do not scale Hairpin Energy by u+2 but by u cause the scale  */
      /* for the closing pair was already done in the forward recursion              */
      for (qbt1=1,s=0; s<n_seq; s++) {
        char loopseq[10];
            type[s] = pair[S[s][q]][S[s][p]];
            if (type[s]==0) type[s]=7;

            if (u<9){
              strcpy(loopseq , sequences[s]+q-1);
              strncat(loopseq, sequences[s], p);
            }
        qbt1 *= exp_E_Hairpin(u, type[s], S[s][q+1], S[s][(p>1) ? p-1 : n], loopseq, pf_params);
      }
      qho += qb[my_iindx[p]-q] * qbt1 * scale[u];

      /* 2. exterior interior loop contribution*/

      for(k=q+1; k < n; k++){
        int ln1, lstart;
        ln1 = k - q - 1;
        if(ln1+p-1>MAXLOOP) break;
        lstart = ln1+p-1+n-MAXLOOP;
        if(lstart<k+TURN+1) lstart = k + TURN + 1;
        for(l=lstart;l <= n; l++){
          int ln2, type_2;

              ln2 = (p - 1) + (n - l);
              if((ln1+ln2) > MAXLOOP) continue;
              double qloop=1.;
              if (qb[my_iindx[k]-l]==0.){ qloop=0.; continue;}

              for (s=0; s<n_seq; s++){
                int ln1a=a2s[s][k-1]-a2s[s][q];
                int ln2a=a2s[s][n]-a2s[s][l]+a2s[s][p-1];
                type_2 = pair[S[s][l]][S[s][k]];
                if (type_2 == 0) type_2 = 7;
                qloop *= exp_E_IntLoop(ln2a, ln1a, type_2, type[s], S3[s][l], S5[s][k], S5[s][p], S3[s][q], pf_params);
              }
              qio += qb[my_iindx[p]-q] * qb[my_iindx[k]-l] * qloop * scale[ln1+ln2];
            }
      } /* end of kl double loop */
    }
  } /* end of pq double loop */

  /* 3. exterior multiloop contribution  */
  for(k=TURN+2; k<n-2*TURN-3; k++)
    qmo += qm[my_iindx[1]-k] * qm2[k+1] * pow(expMLclosing,n_seq);

  /* add additional pf of 1.0 to take open chain into account */
  qo = qho + qio + qmo + 1.0*scale[n];

  free(type);
}


/*brauch ma nurnoch pscores!*/
PUBLIC char *alipbacktrack(double *prob) {
  double r, gr, qt, kTn;
  int k,i,j, start,s,n, n_seq;
  double probs=1;

  if (q == NULL)
    nrerror("can't backtrack without pf arrays.\n"
            "Call pf_fold() before pbacktrack()");

  n = S[0][0];
  n_seq = N_seq;
  kTn = pf_params->kT/10.;
  /*sequence = seq;*/
  if (do_bppm==0) {
    for (k=1; k<=n; k++) {
      qln[k] = q[my_iindx[k] -n];
    }
    qln[n+1] = 1.0;
  }

  pstruc = space((n+1)*sizeof(char));

  for (i=0; i<n; i++) pstruc[i] = '.';

  start = 1;
  while (start<n) {
  /* find i position of first pair */
    probs=1.;
    for (i=start; i<n; i++) {
      gr = urn() * qln[i];
      if (gr > qln[i+1]*scale[1]) {
        *prob=*prob*probs*(1-qln[i+1]*scale[1]/qln[i]);
        break; /* i is paired */
      }
      probs*=qln[i+1]*scale[1]/qln[i];
    }
    if (i>=n) {
      *prob=*prob*probs;
      break; /* no more pairs */
    }
    /* now find the pairing partner j */
    r = urn() * (qln[i] - qln[i+1]*scale[1]);
    for (qt=0, j=i+1; j<=n; j++) {
      int xtype;
      /*  type = ptype[my_iindx[i]-j];
          if (type) {*/
      double qkl;
      if (qb[my_iindx[i]-j]>0) {
        qkl = qb[my_iindx[i]-j]*qln[j+1];  /*if psc too small qb=0!*/
        for (s=0; s< n_seq; s++) {
          xtype=pair[S[s][i]][S[s][j]];
          if (xtype==0) xtype=7;
          qkl *= exp_E_ExtLoop(xtype, (i>1) ? S5[s][i] : -1, (j<n) ? S3[s][j] : -1, pf_params);
        }
        qt += qkl; /*?*exp(pscore[my_iindx[i]-j]/kTn)*/
        if (qt > r) {
          *prob=*prob*(qkl/(qln[i] - qln[i+1]*scale[1]));/*probs*=qkl;*/
          break; /* j is paired */
        }
      }
    }
    if (j==n+1) nrerror("backtracking failed in ext loop");
    start = j+1;
    backtrack(i,j, n_seq, prob); /*?*/
  }

  return pstruc;
}


PRIVATE void backtrack(int i, int j, int n_seq, double *prob) {
  /*backtrack given i,j basepair!*/
  double kTn = pf_params->kT/10.;
  double tempwert;
  int *type = (int *)space(sizeof(int) * n_seq);

  do {
    double r, qbt1;
    int k, l, u, u1,s;
    pstruc[i-1] = '('; pstruc[j-1] = ')';
    for (s=0; s<n_seq; s++) {
      type[s] = pair[S[s][i]][S[s][j]];
      if (type[s]==0) type[s]=7;
    }
    r = urn() * (qb[my_iindx[i]-j]/exp(pscore[my_iindx[i]-j]/kTn)); /*?*exp(pscore[my_iindx[i]-j]/kTn)*/

    qbt1=1.;
    for (s=0; s<n_seq; s++){
      u = a2s[s][j-1]-a2s[s][i];
      if (a2s[s][i]<1) continue;
      char loopseq[10];
      if(u < 7){
        strncpy(loopseq, Ss[s]+a2s[s][i]-1, 10);
      }
      qbt1 *= exp_E_Hairpin(u, type[s], S3[s][i], S5[s][j], loopseq, pf_params);
    }
    qbt1 *= scale[j-i+1];

    if (qbt1>r) {
      *prob=*prob*qbt1/(qb[my_iindx[i]-j]/exp(pscore[my_iindx[i]-j]/kTn));/*probs*=qbt1;*/
      free(type);
      return; /* found the hairpin we're done */
    }

    for (k=i+1; k<=MIN2(i+MAXLOOP+1,j-TURN-2); k++){

      for (l=MAX2(k+TURN+1,j-1-MAXLOOP+k-i-1); l<j; l++){
        double qloop=1;
        int type_2;
        if (qb[my_iindx[k]-l]==0) {qloop=0; continue;}
        for (s=0; s<n_seq; s++) {
          u1 = a2s[s][k-1]-a2s[s][i]/*??*/;
          type_2 = pair[S[s][l]][S[s][k]]; if (type_2 == 0) type_2 = 7;
          qloop *= exp_E_IntLoop(u1, a2s[s][j-1]-a2s[s][l], type[s], type_2,
                                 S3[s][i], S5[s][j],S5[s][k], S3[s][l], pf_params);
        }
        qbt1 += qb[my_iindx[k]-l] * qloop * scale[k-i+j-l];

        if (qbt1 > r) {
         *prob=*prob*qb[my_iindx[k]-l] * qloop * scale[k-i+j-l]/(qb[my_iindx[i]-j]/exp(pscore[my_iindx[i]-j]/kTn));
         /*
          prob*=qb[my_iindx[k]-l] * qloop * scale[k-i+j-l];
         */
          break;
        }
      }
      if (qbt1 > r) break;
    }
    if (l<j) {
      i=k; j=l;
    }
    else {
       *prob=*prob*(1-qbt1/(qb[my_iindx[i]-j]/exp(pscore[my_iindx[i]-j]/kTn)));
      break;
    }
  } while (1);

  /* backtrack in multi-loop */
  tempwert=(qb[my_iindx[i]-j]/exp(pscore[my_iindx[i]-j]/kTn));
  {
    double r, qt;
    int k, ii, jj;
    double qttemp=0;;
    i++; j--;
    /* find the first split index */
    ii = my_iindx[i]; /* ii-j=[i,j] */
    jj = jindx[j]; /* jj+i=[j,i] */
    for (qt=0., k=i+1; k<j; k++) qttemp += qm[ii-(k-1)]*qm1[jj+k];
    r = urn() * qttemp;
    for (qt=0., k=i+1; k<j; k++) {
      qt += qm[ii-(k-1)]*qm1[jj+k];
      if (qt>=r){
        *prob=*prob*qm[ii-(k-1)]*qm1[jj+k]/qttemp;/*qttemp;tempwert*/
        /*        prob*=qm[ii-(k-1)]*qm1[jj+k];*/
        break;
      }
    }
    if (k>=j) nrerror("backtrack failed, can't find split index ");

    backtrack_qm1(k, j, n_seq, prob);

    j = k-1;
    while (j>i) {
      /* now backtrack  [i ... j] in qm[] */
      jj = jindx[j];/*habides??*/
      ii = my_iindx[i];
      r = urn() * qm[ii - j];
      qt = qm1[jj+i]; k=i;
      if (qt<r)
        for (k=i+1; k<=j; k++) {
          qt += (qm[ii-(k-1)]+expMLbase[k-i]/*n_seq??*/)*qm1[jj+k];
          if (qt >= r) {
            *prob=*prob*(qm[ii-(k-1)]+expMLbase[k-i])*qm1[jj+k]/qm[ii - j];/*???*/
            /*            probs*=qt;*/
            break;
          }
        }
      else {
        *prob=*prob*qt/qm[ii - j];/*??*/
      }
      if (k>j) nrerror("backtrack failed in qm");

      backtrack_qm1(k,j, n_seq, prob);

      if (k<i+TURN) break; /* no more pairs */
      r = urn() * (qm[ii-(k-1)] + expMLbase[k-i]);
      if (expMLbase[k-i] >= r) {
        break; /* no more pairs */
        *prob=*prob*expMLbase[k-i]/(qm[ii-(k-1)] + expMLbase[k-i]);
      }
      j = k-1;
      //whatishere??
    }
  }
  free(type);
}

PRIVATE void backtrack_qm1(int i,int j, int n_seq, double *prob) {
  /* i is paired to l, i<l<j; backtrack in qm1 to find l */
  int ii, l, xtype,s;
  double qt, r, tempz;
  r = urn() * qm1[jindx[j]+i];
  ii = my_iindx[i];
  for (qt=0., l=i+TURN+1; l<=j; l++) {
    if (qb[ii-l]==0) continue;
    tempz=1.;
    for (s=0; s<n_seq; s++) {
      xtype = pair[S[s][i]][S[s][l]];
      if (xtype==0) xtype=7;
      tempz*=exp_E_MLstem(xtype, S5[s][i], S3[s][l], pf_params);
    }
    qt +=  qb[ii-l]*tempz*expMLbase[j-l];
    if (qt>=r) {
      *prob=*prob*qb[ii-l]*tempz*expMLbase[j-l]/qm1[jindx[j]+i];
      /* probs*=qb[ii-l]*tempz*expMLbase[j-l];*/
      break;
    }
  }
  if (l>j) nrerror("backtrack failed in qm1");

  backtrack(i,l, n_seq, prob);
}

PUBLIC FLT_OR_DBL *export_ali_bppm(void){
  return probs;
}
