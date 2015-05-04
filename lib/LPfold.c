/* Last changed Time-stamp: <2009-02-18 14:19:51 ivo> */
/*
  local pair probabilities for RNA secondary structures

  Stephan Bernhart, Ivo L Hofacker
  Vienna RNA package
*/
/*
  todo: compute energy z-score for each window

*/
#include <config.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>    /* #defines FLT_MAX ... */
#include "utils.h"
#include "energy_par.h"
#include "fold_vars.h"
#include "pair_mat.h"
#include "PS_dot.h"
#include "part_func.h"
#include "params.h"
#include "loop_energies.h"
#include "LPfold.h"
#include "Lfold.h"

#ifdef _OPENMP
#include <omp.h>
#endif


/*@unused@*/
PRIVATE char rcsid[] UNUSED = "$Id: LPfold.c,v 1.8 2009/02/18 20:34:38 ivo Exp $";

#define ISOLATED  256.0

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

PRIVATE float       cutoff;
PRIVATE int         num_p=0; /* for counting basepairs in pairlist pl, can actually be moved into pfl_fold */
PRIVATE FLT_OR_DBL  *expMLbase=NULL;
PRIVATE FLT_OR_DBL  **q=NULL, **qb=NULL, **qm=NULL, *qqm=NULL, *qqm1=NULL, *qq=NULL, *qq1=NULL, **pR=NULL, **qm2=NULL, **QI5=NULL,  **q2l=NULL, **qmb=NULL;/*,**QI3,*/
PRIVATE FLT_OR_DBL  *prml=NULL, *prm_l=NULL, *prm_l1=NULL, *q1k=NULL, *qln=NULL;
PRIVATE FLT_OR_DBL  *scale=NULL;
PRIVATE char        **ptype=NULL; /* precomputed array of pair types */
PRIVATE int         *jindx=NULL;
PRIVATE int         *my_iindx=NULL;
PRIVATE int         init_length = 0;  /* length in last call to init_pf_fold() */
PRIVATE pf_paramT   *pf_params=NULL;
PRIVATE short       *S=NULL, *S1=NULL;
PRIVATE int         unpaired;
PRIVATE int         ulength;
PRIVATE int         pUoutput;
PRIVATE double      alpha = 1.0;

#ifdef _OPENMP

/* NOTE: all variables are assumed to be uninitialized if they are declared as threadprivate
         thus we have to initialize them before usage by a seperate function!
         OR: use copyin in the PARALLEL directive!
         e.g.:
         #pragma omp parallel for copyin(pf_params)
*/
#pragma omp threadprivate(cutoff, num_p, scale, ptype, jindx, my_iindx, init_length, pf_params,\
                          expMLbase, q, qb, qm, qqm, qqm1, qq, qq1, pR, qm2, QI5, q2l, qmb,\
                          prml, prm_l, prm_l1, q1k, qln,\
                          S, S1, unpaired, ulength, pUoutput, alpha)

#endif

/*
#################################
# PRIVATE FUNCTION DECLARATIONS #
#################################
*/

PRIVATE void  init_partfunc_L(int length, pf_paramT *parameters);
PRIVATE void  get_arrays_L(unsigned int length);
PRIVATE void  free_pf_arrays_L(void);
PRIVATE void  scale_pf_params(unsigned int length, pf_paramT *parameters);
PRIVATE void  GetPtype(int j, int pairsize, const short *S, int n);
PRIVATE void  FreeOldArrays(int i);
PRIVATE void  GetNewArrays(int j, int winSize);
PRIVATE void  printpbar(FLT_OR_DBL **prb,int winSize, int i, int n);
PRIVATE plist *get_deppp(struct plist *pl, int start, int pairsize, int length);
PRIVATE plist *get_plistW(struct plist *pl, int length, int start, FLT_OR_DBL **Tpr, int winSize);
PRIVATE void  print_plist(int length, int start, FLT_OR_DBL **Tpr, int winSize, FILE *fp);
PRIVATE void  compute_pU(int k, int ulength, double **pU, int winSize, int n, char *sequence);
PRIVATE void  putoutpU(double **pU,int k, int ulength, FILE *fp);
/*PRIVATE void make_ptypes(const short *S, const char *structure);*/


/*
#################################
# BEGIN OF FUNCTION DEFINITIONS #
#################################
*/

PRIVATE void init_partfunc_L(int length, pf_paramT *parameters){
  if (length<1) nrerror("init_partfunc_L: length must be greater 0");
#ifdef _OPENMP
/* Explicitly turn off dynamic threads */
  omp_set_dynamic(0);
  free_pf_arrays_L(); /* free previous allocation */
#else
  if (init_length>0) free_pf_arrays_L(); /* free previous allocation */
#endif

#ifdef SUN4
  nonstandard_arithmetic();
#else
#ifdef HP9
  fpsetfastmode(1);
#endif
#endif
  make_pair_matrix();
  get_arrays_L((unsigned) length);
  scale_pf_params((unsigned) length, parameters);

#ifndef _OPENMP
  init_length = length;
#endif
}

PRIVATE void get_arrays_L(unsigned int length){
  /*arrays in 2 dimensions*/

  q         = (FLT_OR_DBL **) space(sizeof(FLT_OR_DBL *)*(length+1));
  qb        = (FLT_OR_DBL **) space(sizeof(FLT_OR_DBL *)*(length+1));
  qm        = (FLT_OR_DBL **) space(sizeof(FLT_OR_DBL *)*(length+1));
  pR        = (FLT_OR_DBL **) space(sizeof(FLT_OR_DBL *)*(length+1));
  q1k       = (FLT_OR_DBL *)  space(sizeof(FLT_OR_DBL)  *(length+1));
  qln       = (FLT_OR_DBL *)  space(sizeof(FLT_OR_DBL)  *(length+2));
  qq        = (FLT_OR_DBL *)  space(sizeof(FLT_OR_DBL)  *(length+2));
  qq1       = (FLT_OR_DBL *)  space(sizeof(FLT_OR_DBL)  *(length+2));
  qqm       = (FLT_OR_DBL *)  space(sizeof(FLT_OR_DBL)  *(length+2));
  qqm1      = (FLT_OR_DBL *)  space(sizeof(FLT_OR_DBL)  *(length+2));
  prm_l     = (FLT_OR_DBL *)  space(sizeof(FLT_OR_DBL)  *(length+2));
  prm_l1    = (FLT_OR_DBL *)  space(sizeof(FLT_OR_DBL)  *(length+2));
  prml      = (FLT_OR_DBL *)  space(sizeof(FLT_OR_DBL)  *(length+2));
  expMLbase = (FLT_OR_DBL *)  space(sizeof(FLT_OR_DBL)  *(length+1));
  scale     = (FLT_OR_DBL *)  space(sizeof(FLT_OR_DBL)  *(length+1));
  ptype     = (char **)       space(sizeof(char *)      *(length+2));

  if (ulength>0) {
    /* QI3 = (FLT_OR_DBL **) space((length+1)*sizeof(FLT_OR_DBL *));*/
    QI5 = (FLT_OR_DBL **) space((length+1)*sizeof(FLT_OR_DBL *));
    qmb = (FLT_OR_DBL **) space((length+1)*sizeof(FLT_OR_DBL *));
    qm2 = (FLT_OR_DBL **) space((length+1)*sizeof(FLT_OR_DBL *));
    q2l = (FLT_OR_DBL **) space((length+1)*sizeof(FLT_OR_DBL *));
  }
  my_iindx  = get_iindx(length);
  iindx     = get_iindx(length); /* for backward compatibility and Perl wrapper */
  jindx     = get_indx(length);
}

PRIVATE void free_pf_arrays_L(void){
  if(q)         free(q);
  if(qb)        free(qb);
  if(qm)        free(qm);
  if(pR)        free(pR);
  if(qm2)       free(qm2);
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
  if(ptype)     free(ptype);
  if(QI5)       free(QI5);
  if(qmb)       free(qmb);
  if(q2l)       free(q2l);
  if(pf_params) free(pf_params);

  q = qb = qm = pR = QI5 = qmb = qm2 = q2l = NULL;
  qq = qq1 = qqm = qqm1 = q1k = qln = prml = prm_l = prm_l1 = expMLbase = NULL;
  my_iindx = jindx = iindx = NULL;
  pf_params = NULL;
  ptype     = NULL;
  scale = NULL;

#ifdef SUN4
  standard_arithmetic();
#else
#ifdef HP9
  fpsetfastmode(0);
#endif
#endif

#ifndef _OPENMP
  init_length=0;
#endif
}

PUBLIC void update_pf_paramsLP(int length){
  update_pf_paramsLP_par(length, NULL);
}

PUBLIC void update_pf_paramsLP_par(int length, pf_paramT *parameters){
#ifdef _OPENMP
  init_partfunc_L(length, parameters);
#else
  if(parameters) init_partfunc_L(length, parameters);
  else if (length > init_length) init_partfunc_L(length, parameters);
  else {
    /*   make_pair_matrix();*/
    scale_pf_params((unsigned) length, parameters);
  }
#endif
}

PUBLIC plist *pfl_fold( char *sequence,
                        int winSize,
                        int pairSize,
                        float cutoffb,
                        double **pU,
                        struct plist **dpp2,
                        FILE *pUfp,
                        FILE *spup){
  return pfl_fold_par(sequence, winSize, pairSize, cutoffb, pU, dpp2, pUfp, spup, NULL);
}

PUBLIC plist *pfl_fold_par( char *sequence,
                            int winSize,
                            int pairSize,
                            float cutoffb,
                            double **pU,
                            struct plist **dpp2,
                            FILE *pUfp,
                            FILE *spup,
                            pf_paramT *parameters){

  int         n, m, i, j, k, l, u, u1, ii, type, type_2, tt, ov, do_dpp, simply_putout, noGUclosure;
  double      max_real;
  FLT_OR_DBL  temp, Qmax, prm_MLb, prmt, prmt1, qbt1, *tmp, expMLclosing;
  plist       *dpp, *pl;

  ov            = 0;
  Qmax          = 0;
  do_dpp        = 0;
  simply_putout = 0;
  dpp           = NULL;
  pl            = NULL;
  pUoutput      = 0;
  ulength       = 0;
  cutoff        = cutoffb;

  if(pU != NULL)  ulength       = (int)pU[0][0]+0.49;
  if(spup !=NULL) simply_putout = 1; /*can't have one without the other*/
  if(pUfp!=NULL)  pUoutput      = 1;
  else if((pUoutput)&&(ulength!=0)){
    fprintf(stderr, "There was a problem with non existing File Pointer for unpaireds, terminating process\n");
    return pl;
  }
  dpp = *dpp2;
  if(dpp !=NULL)  do_dpp=1;

  n = (int) strlen(sequence);
  if (n<TURN+2) return 0;

#ifdef _OPENMP
  /* always init everything since all global static variables are uninitialized when entering a thread */
  init_partfunc_L(n, parameters);
#else
  if(parameters) init_partfunc_L(n, parameters);
  else if (n > init_length) init_partfunc_L(n, parameters);
  else if (fabs(pf_params->temperature - temperature)>1e-6) update_pf_paramsLP_par(n, parameters);
#endif

  expMLclosing  = pf_params->expMLclosing;
  noGUclosure   = pf_params->model_details.noGUclosure;


  max_real = (sizeof(FLT_OR_DBL) == sizeof(float)) ? FLT_MAX : DBL_MAX;

  S   = encode_sequence(sequence, 0);
  S1  = encode_sequence(sequence, 1);

  /*  make_ptypes(S, structure); das machmadochlieber lokal, ey!*/

  /*here, I allocate memory for pU, if has to be saved, I allocate all in one go,
    if pU is put out and freed, I only allocate what I really need*/

  if (ulength>0){
    if (pUoutput) {
      for (i=1; i<=ulength; i++) pU[i]=(double *)space((MAX2(MAXLOOP,ulength)+2)*sizeof(double));
    }
    else {
      for (i=1; i<=n; i++) pU[i]=(double *)space((MAX2(MAXLOOP,ulength)+2)*sizeof(double));
     }
  }

  /*array initialization ; qb,qm,q
    qb,qm,q (i,j) are stored as ((n+1-i)*(n-i) div 2 + n+1-j */
  num_p = 0;
  pl    = (struct plist *)space(1000*sizeof(struct plist));


  /*ALWAYS q[i][j] => i>j!!*/
  for (j=1; j<MIN2(TURN+2,n); j++) { /*allocate start*/
    GetNewArrays(j, winSize);
    GetPtype(j,pairSize,S,n);
    for (i=1; i<=j; i++) q[i][j]=scale[(j-i+1)];
  }
  for (j=TURN+2;j<=n+winSize; j++) {
    if (j<=n) {
      GetNewArrays(j, winSize);
      GetPtype(j,pairSize,S,n);
      for (i=MAX2(1,j-winSize); i<=j/*-TURN*/; i++)
        q[i][j]=scale[(j-i+1)];
      for (i=j-TURN-1;i>=MAX2(1,(j-winSize+1)); i--) {
        /* construction of partition function of segment i,j*/
        /*firstly that given i bound to j : qb(i,j) */
        u = j-i-1;
        type = ptype[i][j];
        if (type!=0) {
          /*hairpin contribution*/
          if (((type==3)||(type==4))&&noGUclosure) qbt1 = 0;
          else
            qbt1 = exp_E_Hairpin(u, type, S1[i+1], S1[j-1], sequence+i-1, pf_params) * scale[u+2];

          /* interior loops with interior pair k,l */
          for (k=i+1; k<=MIN2(i+MAXLOOP+1,j-TURN-2); k++) {
            u1 = k-i-1;
            for (l=MAX2(k+TURN+1,j-1-MAXLOOP+u1); l<j; l++) {
              type_2 = ptype[k][l];
              if (type_2) {
                type_2 = rtype[type_2];
                qbt1 += qb[k][l] *
                  exp_E_IntLoop(u1, j-l-1, type, type_2,
                                S1[i+1], S1[j-1], S1[k-1], S1[l+1], pf_params) * scale[k-i+j-l];
              }
            }
          }
          /*multiple stem loop contribution*/
          ii = my_iindx[i+1]; /* ii-k=[i+1,k-1] */
          temp = 0.0;
          for (k=i+2; k<=j-1; k++) temp += qm[i+1][k-1]*qqm1[k];
          tt = rtype[type];
          qbt1 += temp * expMLclosing * exp_E_MLstem(tt, S1[j-1], S1[i+1], pf_params) * scale[2];

          qb[i][j] = qbt1;
        } /* end if (type!=0) */
        else qb[i][j] = 0.0;

        /* construction of qqm matrix containing final stem
           contributions to multiple loop partition function
           from segment i,j */
        qqm[i] = qqm1[i]*expMLbase[1];
        if (type) {
          qbt1 = qb[i][j] * exp_E_MLstem(type, (i>1) ? S1[i-1] : -1, (j<n) ? S1[j+1] : -1, pf_params);
          qqm[i] += qbt1;
        }

        /*construction of qm matrix containing multiple loop
          partition function contributions from segment i,j */
        temp = 0.0;
        /*ii = my_iindx[i];   ii-k=[i,k-1] */
        /*new qm2 computation done here*/
        for (k=i+1; k<=j; k++) temp += (qm[i][k-1])*qqm[k];
        if (ulength>0) qm2[i][j]=temp;/*new qm2 computation done here*/
        for (k=i+1; k<=j; k++) temp += expMLbase[k-i] * qqm[k];
        qm[i][j] = (temp + qqm[i]);

        /*auxiliary matrix qq for cubic order q calculation below */
        qbt1 = qb[i][j];
        if (type) {
          qbt1 *= exp_E_ExtLoop(type, (i>1) ? S1[i-1] : -1, (j < n) ? S1[j+1] : -1, pf_params);
        }
        qq[i] = qq1[i]*scale[1] + qbt1;

        /*construction of partition function for segment i,j */
        temp = 1.0*scale[1+j-i] + qq[i];
        for (k=i; k<=j-1; k++) temp += q[i][k]*qq[k+1];
        q[i][j] = temp;

        if (temp>Qmax) {
          Qmax = temp;
          if (Qmax>max_real/10.)
            fprintf(stderr, "Q close to overflow: %d %d %g\n", i,j,temp);
        }
        if (temp>=max_real) {
          PRIVATE char msg[128];
          snprintf(msg, 128, "overflow in pf_fold while calculating q[%d,%d]\n"
                  "use larger pf_scale", i,j);
          nrerror(msg);
        }
      } /*end for i*/
      tmp = qq1;  qq1 =qq;  qq =tmp;
      tmp = qqm1; qqm1=qqm; qqm=tmp;
    }

    /* just as a general service, I save here the free energy of the windows
       no output is generated, however,...
    */
    if ((j>=winSize) && (j<=n) && (ulength) && !(pUoutput)) {
      double Fwindow=0.;
      Fwindow=(-log(q[j-winSize+1][j])-winSize*log(pf_params->pf_scale))*pf_params->kT/1000.0;

      pU[j][0]=Fwindow;
      /*
      if (ulength>=winSize)
        pU[j][winSize]=scale[winSize]/q[j-winSize+1][j];
      */
    }
    if (j>winSize) {
      Qmax=0;
      /* i=j-winSize; */
      /* initialize multiloopfs */
      for (k=j-winSize; k<=MIN2(n,j); k++) {
        prml[k]=0;
        prm_l[k]=0;
        /*        prm_l1[k]=0;  others stay*/
      }
      prm_l1[j-winSize]=0;
      k=j-winSize;
      for (l=k+TURN+1; l<=MIN2(n,k+winSize-1); l++) {
        int a;
        pR[k][l] = 0; /* set zero at start */
        type=ptype[k][l];
        if (qb[k][l]==0) continue;

        for (a=MAX2(1,l-winSize+2); a<MIN2(k,n-winSize+2);a++)
          pR[k][l]+=q[a][k-1]*q[l+1][a+winSize-1]/q[a][a+winSize-1];

        if (l-k+1==winSize)
          pR[k][l]+=1./q[k][l];
        else {
          if (k+winSize-1<=n)          /* k outermost */
            pR[k][l]+=q[l+1][k+winSize-1]/q[k][k+winSize-1];
          if (l-winSize+1>=1)  /*l outermost*/
            pR[k][l]+=q[l-winSize+1][k-1]/q[l-winSize+1][l];
        }
        pR[k][l] *= exp_E_ExtLoop(type, (k>1) ? S1[k-1] : -1, (l<n) ? S1[l+1] : -1, pf_params);

        type_2 = ptype[k][l];
        type_2 = rtype[type_2];

        for (i=MAX2(MAX2(l-winSize+1,k-MAXLOOP-1),1); i<=k-1; i++) {
          for (m=l+1; m<=MIN2(MIN2(l+ MAXLOOP -k+i+2,i+winSize-1),n); m++) {
            type = ptype[i][m];
            if ((pR[i][m]>0))
              pR[k][l] += pR[i][m]*exp_E_IntLoop(k-i-1, m-l-1, type, type_2,
                                                 S1[i+1], S1[m-1], S1[k-1], S1[l+1], pf_params) * scale[k-i+m-l];
          }
        }
        if (ulength) { /* NOT IF WITHIN INNER LOOP */
          for (i=MAX2(MAX2(l-winSize+1,k-MAXLOOP-1),1); i<=k-1; i++) {
            for (m=l+1; m<=MIN2(MIN2(l+ MAXLOOP -k+i+2,i+winSize-1),n); m++) {
              type = ptype[i][m];
              if ((pR[i][m]>0)){
                temp=pR[i][m]*qb[k][l]*exp_E_IntLoop(k-i-1, m-l-1, type, type_2,
                                                     S1[i+1], S1[m-1], S1[k-1], S1[l+1], pf_params) * scale[k-i+m-l];
                QI5[l][m-l-1]+=temp;
                QI5[i][k-i-1]+=temp;
              }
            }
           }
        }
      }
      /* 3. bonding k,l as substem of multi-loop enclosed by i,m */
      prm_MLb = 0.;
      if(k>1) /*sonst nix!*/
        for (l=MIN2(n-1,k+winSize-2); l>=k+TURN+1; l--) { /* opposite direction */
          m=l+1;
          prmt = prmt1 = 0.0;
          tt = ptype[k-1][m]; tt=rtype[tt];
          prmt1 = pR[k-1][m] * expMLclosing * exp_E_MLstem(tt, S1[l], S1[k], pf_params);
          for (i=MAX2(1,l-winSize+2); i<k-1/*TURN*/; i++) {
            tt = ptype[i][m]; tt = rtype[tt];
            prmt += pR[i][m] * exp_E_MLstem(tt, S1[m-1], S1[i+1], pf_params) * qm[i+1][k-1];
          }
          tt = ptype[k][l];
          prmt *= expMLclosing;
          prml[ m] = prmt;
          prm_l[m] = prm_l1[m]*expMLbase[1]+prmt1;

          prm_MLb = prm_MLb*expMLbase[1] + prml[m];
          /* same as:    prm_MLb = 0;
             for (i=n; i>k; i--)  prm_MLb += prml[i]*expMLbase[k-i-1];
          */
          prml[m] = prml[ m] + prm_l[m];

          if (qb[k][l] == 0.) continue;

          temp = prm_MLb;

          if (ulength) {
            double dang;
            /* coefficient for computations of unpairedarrays */
            dang  =   qb[k][l] * exp_E_MLstem(tt, S1[k-1], S1[l+1], pf_params) * scale[2];
            for (m=MIN2(k+winSize-2,n);m>=l+2; m--){
              qmb[l][m-l-1] +=  prml[m]*dang;
              q2l[l][m-l-1] +=  (prml[m]-prm_l[m])*dang;
            }
          }

          for (m=MIN2(k+winSize-2,n);m>=l+2; m--)
            temp += prml[m]*qm[l+1][m-1];

          temp      *= exp_E_MLstem(tt, (k>1) ? S1[k-1] : -1, (l<n) ? S1[l+1] : -1, pf_params) * scale[2];
          pR[k][l]  += temp;

          if (pR[k][l]>Qmax) {
            Qmax = pR[k][l];
            if (Qmax>max_real/10.)
              fprintf(stderr, "P close to overflow: %d %d %g %g\n",
                      i, m, pR[k][l], qb[k][l]);
          }
          if (pR[k][l]>=max_real) {
            ov++;
            pR[k][l]=FLT_MAX;
          }

        } /* end for (l=..) */
      tmp = prm_l1; prm_l1=prm_l; prm_l=tmp;

      /* end for (l=..)   */
      if ((ulength)&&(k-MAXLOOP-1>0)){
        /* if (pUoutput) pU[k-MAXLOOP-1]=(double *)space((ulength+2)*sizeof(double)); */
        compute_pU(k-MAXLOOP-1,ulength,pU, winSize, n, sequence);

        /* here, we put out and free pUs not in use any more (hopefully) */
        if (pUoutput)
          putoutpU(pU,k-MAXLOOP-1, ulength, pUfp);
      }

      if (j-(2*winSize+MAXLOOP+1)>0) {
        printpbar(pR,winSize,j-(2*winSize+MAXLOOP+1),n);
        if (simply_putout) {
          print_plist(n, j-(2*winSize+MAXLOOP+1), pR, winSize, spup);
        }
        else{
          pl=get_plistW(pl, n, j-(2*winSize+MAXLOOP+1), pR, winSize);
        }
        if (do_dpp)dpp=get_deppp(dpp,j-(2*winSize-MAXLOOP),pairSize, n);
        FreeOldArrays(j-(2*winSize+MAXLOOP+1));
      }
    }   /* end if (do_backtrack)*/

  }/* end for j */

  /* finish output and free */
  for (j=MAX2(1,n-MAXLOOP); j<=n;j++) {
    /* if (pUoutput) pU[j]=(double *)space((ulength+2)*sizeof(double)); */
    if (ulength) compute_pU(j,ulength,pU, winSize, n, sequence);
    /*here, we put out and free pUs not in use any more (hopefully)*/
    if (pUoutput) putoutpU(pU,j, ulength, pUfp);
  }
  for (j=MAX2(n-winSize-MAXLOOP,1); j<=n; j++) {
    printpbar(pR,winSize,j,n);
    if (simply_putout) {
      print_plist(n, j, pR, winSize, spup);
    }
    else {
      pl=get_plistW(pl, n, j, pR, winSize);
    }
    if ((do_dpp)&&j<n) dpp=get_deppp(dpp,j,pairSize, n);
    FreeOldArrays(j);
  }
  /* free_pf_arrays_L(); */
  free(S);
  free(S1);
  S = S1 = NULL;
  if (ov>0) fprintf(stderr, "%d overflows occurred while backtracking;\n"
                    "you might try a smaller pf_scale than %g\n",
                    ov, pf_params->pf_scale);
  *dpp2=dpp;

  return pl;
}

PRIVATE void scale_pf_params(unsigned int length, pf_paramT *parameters){
  unsigned int i;
  double  kT, scaling_factor;

  if(pf_params) free(pf_params);

  if(parameters){
    pf_params = get_boltzmann_factor_copy(parameters);
  } else {
    model_detailsT  md;
    set_model_details(&md);
    pf_params = get_boltzmann_factors(temperature, alpha, md, pf_scale);
  }

  scaling_factor = pf_params->pf_scale;
  kT = pf_params->kT;   /* kT in cal/mol  */

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

PRIVATE void printpbar(FLT_OR_DBL **prb,int winSize, int i, int n) {
  int j;
  int howoften=0; /* how many samples do we have for this pair */
  int pairdist;

  for (j=i+TURN; j<MIN2(i+winSize,n+1); j++) {
    pairdist=(j-i+1);
    /*4cases*/
    howoften=MIN2(winSize-pairdist+1,i); /*pairdist,start*/
    howoften=MIN2(howoften,n-j+1);       /*end*/
    howoften=MIN2(howoften,n-winSize+1); /*windowsize*/
    prb[i][j] *= qb[i][j]/howoften;
  }
  return;
}

PRIVATE void FreeOldArrays(int i) {
  /*free arrays no longer needed*/
  free(pR[i]+i);
  free(q[i]+i);
  free(qb[i]+i);
  free(qm[i]+i);
  if (ulength!=0) {
    free(qm2[i]+i);
    free(QI5[i]);
    free(qmb[i]);
    free(q2l[i]);
  }
  free(ptype[i]+i);
  return;
}

PRIVATE void GetNewArrays(int j, int winSize) {
  /*allocate new part of arrays*/
  pR[j]=(FLT_OR_DBL *)space((winSize+1)*sizeof(FLT_OR_DBL));
  pR[j]-=j;
  q[j]=(FLT_OR_DBL *)space((winSize+1)*sizeof(FLT_OR_DBL));
  q[j]-=j;
  qb[j]=(FLT_OR_DBL *)space((winSize+1)*sizeof(FLT_OR_DBL));
  qb[j]-=j;
  qm[j]=(FLT_OR_DBL *)space((winSize+1)*sizeof(FLT_OR_DBL));
  qm[j]-=j;
  if (ulength!=0) {
    qm2[j]=(FLT_OR_DBL *)space((winSize+1)*sizeof(FLT_OR_DBL));
    qm2[j]-=j;
    QI5[j]=(FLT_OR_DBL *)space((winSize+1)*sizeof(FLT_OR_DBL));
    qmb[j]=(FLT_OR_DBL *)space((winSize+1)*sizeof(FLT_OR_DBL));
    q2l[j]=(FLT_OR_DBL *)space((winSize+1)*sizeof(FLT_OR_DBL));
  }
  ptype[j]=(char *)space((winSize+1)*sizeof(char));
  ptype[j]-=j;
  return;
}


PRIVATE void GetPtype(int i, int winSize,const short *S,int n) {
  /*make new entries in ptype array*/
  int j;
  int type;
  for (j=i; j<=MIN2(i+winSize,n); j++) {
    type = pair[S[i]][S[j]];
    ptype[i][j] = (char) type;
  }
  return;
}


PRIVATE plist *get_plistW(plist *pl, int length,
                                 int start, FLT_OR_DBL **Tpr, int winSize) {
  /* get pair probibilities out of pr array */
  int  j,  max_p;
  max_p=1000;
  while (max_p<num_p)
    max_p*=2;

  for (j=start+1; j<=MIN2(start+winSize, length); j++) {
    if (Tpr[start][j]<cutoff) continue;
    if (num_p==max_p-1) {
      max_p*=2;
      pl=(plist *)xrealloc(pl,max_p*sizeof(plist));
    }
    pl[num_p].i=start;
    pl[num_p].j=j;
    pl[num_p++].p=Tpr[start][j];
  }

  /* mark end of data with zeroes */
  pl[num_p].i=0;
  pl[num_p].j=0;
  pl[num_p].p=0.;
  /* pl=(struct plist *)xrealloc(pl,(count)*sizeof(struct plist)); */
  return pl;
}


PRIVATE plist *get_deppp(plist *pl, int start, int pairsize, int length) {
  /* compute dependent pair probabilities */
  int i, j, count=0;
  double tmp;
  plist *temp;
  temp=(plist *)space(pairsize*sizeof(plist)); /* holds temporary deppp */
  for (j=start+TURN; j<MIN2(start+pairsize,length); j++) {

    if ((qb[start][j]*qb[start-1][(j+1)])>10e-200) {
      int type=ptype[start-1][j+1];
      int type_2=rtype[ptype[start][j]];
      tmp=qb[start][j]/qb[start-1][(j+1)]*exp_E_IntLoop(0, 0, type, type_2,
                                                        S1[start], S1[j], S1[start-1], S1[j+1], pf_params) * scale[2];
       temp[count].i=start;
      temp[count].j=j;
      temp[count++].p=tmp;
    }
  }
  /* write it to list of deppps */
  for (i=0; pl[i].i!=0; i++);
  pl=(plist *)xrealloc(pl,(i+count+1)*sizeof(plist));
  for (j=0; j<count; j++) {
    pl[i+j].i=temp[j].i;
    pl[i+j].j=temp[j].j;
    pl[i+j].p=temp[j].p;
  }
  pl[i+count].i=0;
  pl[i+count].j=0;
  pl[i+count].p=0;
  free(temp);
  return pl;
}


PRIVATE void print_plist(int length,int start, FLT_OR_DBL **Tpr, int winSize, FILE *fp) {
  /* print out of pr array, do not save */
  int  j;


  for (j=start+1; j<=MIN2(start+winSize, length); j++) {
    if (Tpr[start][j]<cutoff) continue;
    fprintf(fp,"%d  %d  %g\n",start,j,Tpr[start][j]);
  }

  /* mark end of data with zeroes */

  return ;
}

PRIVATE void compute_pU(int k, int ulength, double **pU, int winSize,int n, char *sequence) {
/*  here, we try to add a function computing all unpaired probabilities starting at some i,
    going down to $unpaired, to be unpaired, i.e. a list with entries from 1 to unpaired for
    every i, with the probability of a stretch of length x, starting at i-x+1, to be unpaired
*/
  int startu;
  int i5;
  int j3, len, obp;
  double temp;
  double *QBE;
  FLT_OR_DBL  expMLclosing      = pf_params->expMLclosing;

  QBE=(double *) space((MAX2(ulength,MAXLOOP)+2)*sizeof(double));

  /* first, we will */
  /* for k<=ulength, pU[k][k]=0, because no bp can enclose it */
  if (pUoutput&&k+ulength<=n)  pU[k+ulength]=(double *)space((ulength+2)*sizeof(double));
  /*compute pu[k+ulength][ulength] */
   for (i5=MAX2(k+ulength-winSize+1,1);i5<=k;i5++) {
    for (j3=k+ulength+1; j3<=MIN2(n,i5+winSize-1); j3++) {
      /*  if (k>400) {
        printf("i%d j%d  ",i5,j3);
        fflush(stdout);
        } */
      if (ptype[i5][j3]!=0) {/**/
        /* (.. >-----|..........)
          i5  j     j+ulength  j3              */
        /*Multiloops*/
        temp = (i5<k) ? qm2[i5+1][k] * expMLbase[j3-k-1] : 0.; /* (..{}{}-----|......) */

        if(j3-1>k+ulength)
          temp  +=  qm2[k+ulength+1][j3-1] * expMLbase[k+ulength-i5]; /* (..|-----|{}{}) */

        if((i5<k)&&(j3-1>k+ulength))
          temp  +=  qm[i5+1][k] * qm[k+ulength+1][j3-1] * expMLbase[ulength]; /* ({}|-----|{}) */

        /* add dangles, multloopclosing etc. */
        temp  *=  exp_E_MLstem(rtype[ptype[i5][j3]], S1[j3-1], S1[i5+1], pf_params) * scale[2] * expMLclosing;
        /*add hairpins*/
        temp  +=  exp_E_Hairpin(j3-i5-1, ptype[i5][j3], S1[i5+1], S1[j3-1], sequence+i5-1, pf_params) * scale[j3-i5+1];
        /*add outer probability*/
        temp *= pR[i5][j3];
        pU[k+ulength][ulength] += temp;

      }
    }
   }
   /* code doubling to avoid if within loop */
#if 0
  /*initialization for interior loops,
    it is not recomended to have verysmall ulengths!!*/
  if (ulength<MAXLOOP) {
    int k5;
    int l3;
    int outype;
    /* kl bp is 5' */
    /* MAXLOOP>((l5-k5-1)+(j3-l3-1)
      k-winSize+ulength<i5<k-TURN-1;
      k+ulength<j3<=k+MAXLOOP+1
      if i then use l3, it is easier by far:
      j3-MAXLOOP<=l3<=k
      i5<k5<k-TURN k5<=i5+l3+2+MAXLOOP-j3
      k5+TURN<l3<=k
    */
    for (i5=MAX2(k+ulength-winSize,1);i5<k-TURN-1;i5++) {

      for (j3=k+ulength+1; j3<=MIN2(n,MIN2(i5+winSize-1,k+MAXLOOP+1)); j3++) {
        double temp=0;
        if (outype=ptype[i5][j3]>0) /* oder so halt */
          for (l3=MAX2(i5+TURN+1,j3-MAXLOOP-1); l3<=k; l3++){
            for (k5=i5+1; k5<=MIN2(l3-TURN-1,MAXLOOP+i5+l3+2-j3); k5++){
              if (ptype[k5][l3]) {
                temp+= qb[k5][l3]*expLoopEnergy(k5-i5-1, j3-l3-1, outype, rtype[ptype[k5][l3]], S1[i5+1], S1[j3-1], S1[k5-1], S1[l3+1]);
              }
            }
          }
        temp*=pR[i5][j3];
        pU[k+ulength][ulength]+= temp;
      }
    }
    /* kl bp is 3' */
    /*
      k+ulength-MAXLOOP<=i5<=k
      k+ulength+1+TURN<j3<i5+winSize
      k+ulength+1<=k5<i5+MAXLOOP+2 || k5<j3-TURN
      k5<l3<j3 || j3-k5-i5-2-ML<=l3<j3
    */
    for (i5=MAX2(1,MAX2(k+ulength-winSize,k+ulength-MAXLOOP));i5<=k; i5++){
      for (j3=k+ulength+TURN+2; j3<MIN2(n+1,i5+winSize); j3++) {
        double temp = 0;
        if (outype=ptype[i5][j3]>0) /* oder so halt */
          for (k5=k+ulength+1; k5<MIN2(j3-TURN-1,i5+MAXLOOP+2); k5++) {
            for (l3=MAX2(k5+TURN+1,j3+k5-i5-2-MAXLOOP); l3<j3; l3++) {
              if (ptype[k5][l3])
                temp += qb[k5][l3]*expLoopEnergy(k5-i5-1, j3-l3-1, outype, rtype[ptype[k5][l3]], S1[i5+1], S1[j3-1], S1[k5-1], S1[l3+1]);
            }
          }
        temp*=pR[i5][j3];
        pU[k+ulength][ulength]+= temp;
      }
    }
  }
  /* Add up Is QI5[l][m-l-1] QI3 */
  /* Add up Interior loop terms */
  temp=0.;

  for (len=winSize; len>=ulength; len--) temp+=QI3[k][len];
  for (;len>0; len--) {
    temp += QI3[k][len];
    QBE[len] += temp;
  }
#endif
  temp=0.;
  for (len=winSize; len>=MAX2(ulength,MAXLOOP); len--) temp+=QI5[k][len];
  for (;len>0; len--) {
    temp += QI5[k][len];
    QBE[len] += temp;  /* replace QBE with QI */
  }
  /* Add Hairpinenergy to QBE */
  temp=0.;
  for(obp = MIN2(n, k + winSize - 1); obp > k + ulength; obp--)
    if(ptype[k][obp])
      temp += pR[k][obp] * exp_E_Hairpin(obp-k-1, ptype[k][obp], S1[k+1], S1[obp-1], sequence+k-1, pf_params) * scale[obp-k+1];
  for(obp = MIN2(n, MIN2(k + winSize - 1, k + ulength)); obp > k + 1; obp--){
    if (ptype[k][obp])
      temp += pR[k][obp] * exp_E_Hairpin(obp-k-1, ptype[k][obp], S1[k+1], S1[obp-1], sequence+k-1, pf_params) * scale[obp-k+1];
    QBE[obp-k-1] += temp;  /* add hairpins to QBE (all in one array) */
  }
  /* doubling the code to get the if out of the loop */

  /* Add up Multiloopterms  qmb[l][m]+=prml[m]*dang;
    q2l[l][m]+=(prml[m]-prm_l[m])*dang; */

  temp=0.;
  for(len = winSize; len >= ulength; len--)
    temp += q2l[k][len] * expMLbase[len];
  for( ; len > 0; len--){
    temp += q2l[k][len] * expMLbase[len];
    QBE[len] += temp; /* add (()()____) type cont. to I3 */
  }
  for(len = 1; len < ulength; len++){
    for(obp = k + len + TURN; obp <= MIN2(n, k + winSize - 1); obp++){
      /* add (()___()) */
      QBE[len] += qmb[k][obp-k-1] * qm[k+len+1/*2*/][obp-1] * expMLbase[len];
    }
  }
  for (len=1; len<ulength; len++) {
    for (obp=k+len+TURN+TURN; obp<=MIN2(n,k+winSize-1); obp++) {
      if (ptype[k][obp]) {
        temp      =   exp_E_MLstem(rtype[ptype[k][obp]], S1[obp-1], S1[k+1], pf_params) * scale[2] * expMLbase[len] * expMLclosing; /* k:obp */
        QBE[len]  +=  pR[k][obp] * temp * qm2[k+len+1][obp-1]; /* add (___()()) */
      }
    }
  }
  /* After computing all these contributions in QBE[len], that k is paired
    and the unpaired stretch is AT LEAST len long, we start to add that to
    the old unpaired thingies; */
  for(len = 1; len < MIN2(MAX2(ulength, MAXLOOP), n - k); len++){
    pU[k+len][len] += pU[k+len][len+1] + QBE[len];
  }

  /* open chain */
  if ((ulength>=winSize)&&(k>=ulength)) {
    pU[k][winSize]=scale[winSize]/q[k-winSize+1][k];
  }
  /*open chain*/
  if ((ulength>=winSize)&&(k>=ulength)) {
    pU[k][winSize]=scale[winSize]/q[k-winSize+1][k];
  }
  /* now the not enclosed by any base pair terms for whatever it is we do not need anymore...
    ... which should be e.g; k, again */
  for(startu = MIN2(ulength, k); startu > 0; startu--){
    temp=0.;
    for(i5 = MAX2(1, k - winSize + 2); i5 <= MIN2(k - startu, n - winSize + 1); i5++){
      temp += q[i5][k - startu] * q[k + 1][i5 + winSize - 1] * scale[startu]/q[i5][i5 + winSize - 1];
    }
    /* the 2 Cases where the borders are on the edge of the interval */
    if((k >= winSize) && (startu + 1 <= winSize))
      temp += q[k - winSize + 1][k - startu]*scale[startu]/q[k - winSize + 1][k];
    if((k <= n - winSize+ startu) && (k - startu >= 0) && (k < n) && (startu + 1 <= winSize))
      temp += q[k + 1][k - startu + winSize] * scale[startu] / q[k - startu + 1][k - startu + winSize];

    /* Divide by number of possible windows */
    pU[k][startu] += temp;
    {
      int leftmost, rightmost;

      leftmost      = MAX2(1, k - winSize + 1);
      rightmost     = MIN2(n - winSize + 1, k - startu + 1);
      pU[k][startu] /= (rightmost - leftmost + 1);
    }
  }
  free(QBE);
  return;
}


PRIVATE void putoutpU(double **pUx, int k, int ulength, FILE *fp) {
  /*put out unpaireds for k, and free pU[k], make sure we don't need pU[k] any more!!*/
  /*could use that for hairpins, also!*/
  int i;
  fprintf(fp,"%d\t",k);
  for (i=1; i<=MIN2(ulength,k); i++) {
    fprintf(fp,"%.5g\t",pUx[k][i]);
  }
  fprintf(fp,"\n");
  free(pUx[k]);
}

PUBLIC void putoutpU_prob(double **pU,int length, int ulength, FILE *fp, int energies) {
  /*put out unpaireds */
  int i,k;
  double kT= pf_params->kT/1000.0;
  double temp;
  if (energies) fprintf(fp,"#opening energies\n #i$\tl=");
  else  fprintf(fp,"#unpaired probabilities\n #i$\tl=");
  for (i=1; i<=ulength; i++) {
    fprintf(fp,"%d\t", i);
  }
  fprintf(fp,"\n");

  for (k=1; k<=length; k++){
    fprintf(fp,"%d\t",k);
    for (i=1; i<=ulength; i++) {
      if (i>k) {
        fprintf(fp,"NA\t");
        continue;
      }
      if (energies) temp=-log(pU[k][i])*kT;
      else temp=pU[k][i];
      fprintf(fp,"%.7g\t",temp);
    }
    fprintf(fp,"\n");
    free(pU[k]);
  }
  fflush(fp);
}

PUBLIC void putoutpU_prob_bin(double **pU,int length, int ulength, FILE *fp, int energies) {
  /*put out unpaireds */
  int i,k;
  double kT= pf_params->kT/1000.0;
  double temp;
  int *p;
  p = (int*) space(sizeof(int)*1);
  //write first line
  p[0]=ulength; //u length
  fwrite(p,sizeof(int),1,fp);
  p[0]=length; //seq length
  fwrite(p,sizeof(int),1,fp);
  for (k=3; k<=(length+20); k++){ //all the other lines are set to 1000000 because we are at ulength=0
    p[0]=1000000;
    fwrite(p,sizeof(int),1,fp);
  }
  //data
  for (i=1; i<=ulength; i++) {
    for (k=1; k<=11; k++){//write first ten entries to 1000000
      p[0]=1000000;
      fwrite(p,sizeof(int),1,fp);
    }
    for (k=1; k<=length; k++){//write data now
      if (i>k) {
	p[0]=1000000;         //check if u > pos
	fwrite(p,sizeof(int),1,fp);
	continue;
      }
      else{
	p[0]= (int) rint(100 *(-log(pU[k][i])*kT));
	fwrite(p,sizeof(int),1,fp);
      }
    }
    for (k=1; k<=9; k++){//finish by writing the last 10 entries
      p[0]=1000000;
      fwrite(p,sizeof(int),1,fp);
    }
  }
  //free pU array;
  for (k=1; k<=length; k++){
    free(pU[k]);
  }
  free(p);
  fflush(fp);
}


/*###########################################*/
/*# deprecated functions below              #*/
/*###########################################*/

PUBLIC void init_pf_foldLP(int length){ /* DO NOTHING */}

/*
 Here: Space for questions...
*/
