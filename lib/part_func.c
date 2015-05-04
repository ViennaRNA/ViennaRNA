/* Last changed Time-stamp: <2007-12-05 13:52:28 ronny> */
/*
                  partiton function for RNA secondary structures

                  Ivo L Hofacker
                  Vienna RNA package
*/
/*
  $Log: part_func.c,v $
  Revision 1.29  2008/02/23 10:10:49  ivo
  list returned from StackProb was sometimes not terminated correctly

  Revision 1.28  2008/01/08 15:08:10  ivo
  circular fold would fail for open chain

  Revision 1.27  2007/12/05 13:04:04  ivo
  add various circfold variants from Ronny

  Revision 1.26  2007/09/19 12:41:56  ivo
  add computation of centroid() structure for RNAfold -p

  Revision 1.25  2007/04/30 15:12:00  ivo
  merge RNAup into package

  Revision 1.24  2007/03/03 17:57:44  ivo
  make sure entries in scale[] decrease to 0

  Revision 1.23  2006/12/01 12:40:23  ivo
  undo Ulli's accidental commit

  Revision 1.21  2006/08/04 15:39:06  ivo
  new function stackProb returns probability for stacks
  p[(i,j)(i+1,j-1)]

  Revision 1.20  2004/08/12 12:14:46  ivo
  update

  Revision 1.19  2004/05/14 16:28:05  ivo
  fix the bug in make_ptype here too (fixed in 1.27 of fold.c)

  Revision 1.18  2004/02/17 10:46:52  ivo
  make sure init_pf_fold is called before scale_parameters

  Revision 1.17  2004/02/09 18:37:59  ivo
  new mean_bp_dist() function to compute ensemble diversity

  Revision 1.16  2003/08/04 09:14:09  ivo
  finish up stochastic backtracking

  Revision 1.15  2002/03/19 16:51:12  ivo
  more on stochastic backtracking (still incomplete)

  Revision 1.14  2002/02/08 17:37:23  ivo
  set freed S,S1 pointers to NULL

  Revision 1.13  2001/11/16 17:30:04  ivo
  add stochastic backtracking (still incomplete)
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
#include "params.h"
#include "loop_energies.h"
#include "part_func.h"

#ifdef _OPENMP
#include <omp.h>
#endif

/*@unused@*/
static char rcsid[] UNUSED = "$Id: part_func.c,v 1.29 2008/02/23 10:10:49 ivo Exp $";

#define ISOLATED  256.0

/*
#################################
# GLOBAL VARIABLES              #
#################################
*/
PUBLIC  int     st_back=0;

/*
#################################
# PRIVATE VARIABLES             #
#################################
*/
PRIVATE FLT_OR_DBL  *q, *qb=NULL, *qm, *qm1, *qqm, *qqm1, *qq, *qq1;
PRIVATE FLT_OR_DBL  *probs, *prml, *prm_l, *prm_l1, *q1k, *qln;
PRIVATE FLT_OR_DBL  *scale;
PRIVATE FLT_OR_DBL  *expMLbase;
PRIVATE FLT_OR_DBL  qo, qho, qio, qmo, *qm2;
PRIVATE int         *jindx;
PRIVATE int         init_length;  /* length in last call to init_pf_fold() */
PRIVATE int         circular=0;
PRIVATE char        *pstruc;
PRIVATE char        *sequence;
PRIVATE char        *ptype;       /* precomputed array of pair types */
PRIVATE pf_paramT   *pf_params;   /* the precomputed Boltzmann weights */
PRIVATE short       *S, *S1;

#ifdef _OPENMP

/* NOTE: all variables are assumed to be uninitialized if they are declared as threadprivate
         thus we have to initialize them before usage by a seperate function!
         OR: use copyin in the PARALLEL directive!
         e.g.:
         #pragma omp parallel for copyin(pf_params)
*/
#pragma omp threadprivate(q, qb, qm, qm1, qqm, qqm1, qq, qq1, prml, prm_l, prm_l1, q1k, qln,\
                          probs, scale, expMLbase, qo, qho, qio, qmo, qm2, jindx, init_length,\
                          circular, pstruc, sequence, ptype, pf_params, S, S1)

#endif

/*
#################################
# PRIVATE FUNCTION DECLARATIONS #
#################################
*/
PRIVATE void  init_partfunc(int length);
PRIVATE void  scale_pf_params(unsigned int length);
PRIVATE void  get_arrays(unsigned int length);
PRIVATE void  make_ptypes(const short *S, const char *structure);
PRIVATE void  pf_circ(const char *sequence, char *structure);
PRIVATE void  pf_linear(const char *sequence, char *structure);
PRIVATE void  pf_create_bppm(const char *sequence, char *structure);
PRIVATE void  backtrack(int i, int j);
PRIVATE void  backtrack_qm(int i, int j);
PRIVATE void  backtrack_qm1(int i,int j);
PRIVATE void  backtrack_qm2(int u, int n);

/*
#################################
# BEGIN OF FUNCTION DEFINITIONS #
#################################
*/

PRIVATE void init_partfunc(int length){
  if (length<1) nrerror("init_pf_fold: length must be greater 0");

#ifdef _OPENMP
/* Explicitly turn off dynamic threads */
  omp_set_dynamic(0);
  free_pf_arrays(); /* free previous allocation */
#else
  if (init_length>0) free_pf_arrays(); /* free previous allocation */
#endif

#ifdef SUN4
  nonstandard_arithmetic();
#else
#ifdef HP9
  fpsetfastmode(1);
#endif
#endif
  make_pair_matrix();
  get_arrays((unsigned) length);
  scale_pf_params((unsigned) length);

  init_length = length;
}

PRIVATE void get_arrays(unsigned int length){
  unsigned int size;

  size  = sizeof(FLT_OR_DBL) * ((length+1)*(length+2)/2);

  q     = (FLT_OR_DBL *) space(size);
  qb    = (FLT_OR_DBL *) space(size);
  qm    = (FLT_OR_DBL *) space(size);
  qm1   = (st_back || circular) ? (FLT_OR_DBL *) space(size) : NULL;
  qm2   = (circular) ? (FLT_OR_DBL *) space(sizeof(FLT_OR_DBL)*(length+2)) : NULL;
  probs = (FLT_OR_DBL *) space(size);

  ptype     = (char *) space(sizeof(char)*((length+1)*(length+2)/2));
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

  iindx     = get_iindx(length);
  jindx     = get_indx(length);
}

/**
*** Allocate memory for all matrices and other stuff
**/
PUBLIC void free_pf_arrays(void){
  if(q)         free(q);
  if(qb)        free(qb);
  if(qm)        free(qm);
  if(qm1)       free(qm1);
  if(qm2)       free(qm2);
  if(ptype)     free(ptype);
  if(qq)        free(qq);
  if(qq1)       free(qq1);
  if(qqm)       free(qqm);
  if(qqm1)      free(qqm1);
  if(q1k)       free(q1k);
  if(qln)       free(qln);
  if(probs)     free(probs);
  if(prm_l)     free(prm_l);
  if(prm_l1)    free(prm_l1);
  if(prml)      free(prml);
  if(expMLbase) free(expMLbase);
  if(scale)     free(scale);
  if(iindx)     free(iindx);
  if(jindx)     free(jindx);
  if(S)         free(S);
  if(S1)        free(S1);
  if(pr)        free(pr);

  S = S1 = NULL;
  q = pr = probs = qb = qm = qm1 = qm2 = qq = qq1 = qqm = qqm1 = q1k = qln = prm_l = prm_l1 = prml = expMLbase = scale = NULL;
  iindx = jindx = NULL;
  ptype = NULL;

#ifdef SUN4
  standard_arithmetic();
#else
#ifdef HP9
  fpsetfastmode(0);
#endif
#endif

  init_length = 0;
}

/*-----------------------------------------------------------------*/
PUBLIC float pf_fold(const char *sequence, char *structure){

  FLT_OR_DBL  Q;
  double      free_energy;
  int         n = (int) strlen(sequence);

  circular = 0;

#ifdef _OPENMP
  /* always init everything since all global static variables are uninitialized when entering a thread */
  init_partfunc(n);
#else
  if (n > init_length) init_partfunc(n);
#endif
  if (fabs(pf_params->temperature - temperature)>1e-6) update_pf_params(n);

  S   = encode_sequence(sequence, 0);
  S1  = encode_sequence(sequence, 1);

  make_ptypes(S, structure);

  /* do the linear pf fold and fill all matrices  */
  pf_linear(sequence, structure);


  if (backtrack_type=='C')      Q = qb[iindx[1]-n];
  else if (backtrack_type=='M') Q = qm[iindx[1]-n];
  else Q = q[iindx[1]-n];

  /* ensemble free energy in Kcal/mol              */
  if (Q<=FLT_MIN) fprintf(stderr, "pf_scale too large\n");
  free_energy = (-log(Q)-n*log(pf_scale))*pf_params->kT/1000.0;
  /* in case we abort because of floating point errors */
  if (n>1600) fprintf(stderr, "free energy = %8.2f\n", free_energy);

  /* calculate base pairing probability matrix (bppm)  */
  if(do_backtrack){
    pf_create_bppm(sequence, structure);
    /*
    *  Backward compatibility:
    *  This block may be removed if deprecated functions
    *  relying on the global variable "pr" vanish from within the package!
    */
    {
      if(pr) free(pr);
      pr = (FLT_OR_DBL *) space(sizeof(FLT_OR_DBL) * ((n+1)*(n+2)/2));
      memcpy(pr, probs, sizeof(FLT_OR_DBL) * ((n+1)*(n+2)/2));
    }
  }
  return free_energy;
}

PUBLIC float pf_circ_fold(const char *sequence, char *structure){

  FLT_OR_DBL Q;

  double free_energy;
  int n = (int) strlen(sequence);

  circular = 1;
#ifdef _OPENMP
  /* always init everything since all global static variables are uninitialized when entering a thread */
  init_partfunc(n);
#else
  if (n >init_length) init_partfunc(n);
#endif
  if (fabs(pf_params->temperature - temperature)>1e-6) update_pf_params(n);

  S   = encode_sequence(sequence, 0);
  S1  = encode_sequence(sequence, 1);

  make_ptypes(S, structure);

  /* do the linear pf fold and fill all matrices  */
  pf_linear(sequence, structure);

  /* calculate post processing step for circular  */
  /* RNAs                                          */
  pf_circ(sequence, structure);

  if (backtrack_type=='C')      Q = qb[iindx[1]-n];
  else if (backtrack_type=='M') Q = qm[iindx[1]-n];
  else Q = qo;

  /* ensemble free energy in Kcal/mol              */
  if (Q<=FLT_MIN) fprintf(stderr, "pf_scale too large\n");
  free_energy = (-log(Q)-n*log(pf_scale))*pf_params->kT/1000.0;
  /* in case we abort because of floating point errors */
  if (n>1600) fprintf(stderr, "free energy = %8.2f\n", free_energy);

  /* calculate base pairing probability matrix (bppm)  */
  if(do_backtrack){
    pf_create_bppm(sequence, structure);
    /*
    *  Backward compatibility:
    *  This block may be removed if deprecated functions
    *  relying on the global variable "pr" vanish from within the package!
    */
    {
      if(pr) free(pr);
      pr = (FLT_OR_DBL *) space(sizeof(FLT_OR_DBL) * ((n+1)*(n+2)/2));
      memcpy(pr, probs, sizeof(FLT_OR_DBL) * ((n+1)*(n+2)/2));
    }
  }
  return free_energy;
}

PRIVATE void pf_linear(const char *sequence, char *structure){

  int n, i,j,k,l, ij, u,u1,d,ii, type, type_2, tt;
  FLT_OR_DBL temp, Qmax=0;
  FLT_OR_DBL qbt1, *tmp;

  FLT_OR_DBL  expMLclosing = pf_params->expMLclosing;
  double      max_real;

  max_real = (sizeof(FLT_OR_DBL) == sizeof(float)) ? FLT_MAX : DBL_MAX;

  n = (int) strlen(sequence);

  /*array initialization ; qb,qm,q
    qb,qm,q (i,j) are stored as ((n+1-i)*(n-i) div 2 + n+1-j */

  for (d=0; d<=TURN; d++)
    for (i=1; i<=n-d; i++) {
      j=i+d;
      ij = iindx[i]-j;
      q[ij]=1.0*scale[d+1];
      qb[ij]=qm[ij]=0.0;
    }

  for (i=1; i<=n; i++)
    qq[i]=qq1[i]=qqm[i]=qqm1[i]=0;

  for (j=TURN+2;j<=n; j++) {
    for (i=j-TURN-1; i>=1; i--) {
      /* construction of partition function of segment i,j*/
      /*firstly that given i binds j : qb(i,j) */
      u = j-i-1; ij = iindx[i]-j;
      type = ptype[ij];
      if (type!=0) {
        /*hairpin contribution*/
        if (((type==3)||(type==4))&&no_closingGU) qbt1 = 0;
        else
          qbt1 = exp_E_Hairpin(u, type, S1[i+1], S1[j-1], sequence+i-1, pf_params) * scale[u+2];
        /* interior loops with interior pair k,l */
        for (k=i+1; k<=MIN2(i+MAXLOOP+1,j-TURN-2); k++) {
          u1 = k-i-1;
          for (l=MAX2(k+TURN+1,j-1-MAXLOOP+u1); l<j; l++) {
            type_2 = ptype[iindx[k]-l];
            if (type_2) {
              type_2 = rtype[type_2];
              qbt1 += qb[iindx[k]-l] * (scale[u1+j-l+1] *
                                        exp_E_IntLoop(u1, j-l-1, type, type_2,
                                        S1[i+1], S1[j-1], S1[k-1], S1[l+1], pf_params));
            }
          }
        }
        /*multiple stem loop contribution*/
        ii = iindx[i+1]; /* ii-k=[i+1,k-1] */
        temp = 0.0;
        for (k=i+2; k<=j-1; k++) temp += qm[ii-(k-1)]*qqm1[k];
        tt = rtype[type];
        qbt1 += temp * expMLclosing * exp_E_MLstem(tt, S1[j-1], S1[i+1], pf_params) * scale[2];
        qb[ij] = qbt1;
      } /* end if (type!=0) */
      else qb[ij] = 0.0;

      /* construction of qqm matrix containing final stem
         contributions to multiple loop partition function
         from segment i,j */
      qqm[i] = qqm1[i]*expMLbase[1];
      if (type) {
        qbt1 = qb[ij] * exp_E_MLstem(type, ((i>1) || circular) ? S1[i-1] : -1, ((j<n) || circular) ? S1[j+1] : -1, pf_params);
        qqm[i] += qbt1;
      }
      if (qm1) qm1[jindx[j]+i] = qqm[i]; /* for stochastic backtracking and circfold */

      /*construction of qm matrix containing multiple loop
        partition function contributions from segment i,j */
      temp = 0.0;
      ii = iindx[i];  /* ii-k=[i,k-1] */
      for (k=j; k>i; k--) temp += (qm[ii-(k-1)] + expMLbase[k-i])*qqm[k];
      qm[ij] = (temp + qqm[i]);

      /*auxiliary matrix qq for cubic order q calculation below */
      qbt1 = qb[ij];
      if(type)
        qbt1 *= exp_E_ExtLoop(type, ((i>1) || circular) ? S1[i-1] : -1, ((j<n) || circular) ? S1[j+1] : -1, pf_params);

      qq[i] = qq1[i]*scale[1] + qbt1;

      /*construction of partition function for segment i,j */
      temp = 1.0*scale[1+j-i] + qq[i];
      for (k=i; k<=j-1; k++) temp += q[ii-k]*qq[k+1];
      q[ij] = temp;
      if (temp>Qmax) {
        Qmax = temp;
        if (Qmax>max_real/10.)
          fprintf(stderr, "Q close to overflow: %d %d %g\n", i,j,temp);
      }
      if (temp>=max_real) {
        PRIVATE char msg[128];
        sprintf(msg, "overflow in pf_fold while calculating q[%d,%d]\n"
                     "use larger pf_scale", i,j);
        nrerror(msg);
      }
    }
    tmp = qq1;  qq1 =qq;  qq =tmp;
    tmp = qqm1; qqm1=qqm; qqm=tmp;
  }
}

/* calculate partition function for circular case */
/* NOTE: this is the postprocessing step ONLY     */
/* You have to call pf_linear first to calculate  */
/* complete circular case!!!                      */
PRIVATE void pf_circ(const char *sequence, char *structure){

  int u, p, q, k, l;
  int n = (int) strlen(sequence);

  FLT_OR_DBL  qot;
  FLT_OR_DBL  expMLclosing  = pf_params->expMLclosing;

  qo = qho = qio = qmo = 0.;
  /* construct qm2 matrix from qm1 entries  */
  for(k=1; k<n-TURN-1; k++){
    qot = 0.;
    for (u=k+TURN+1; u<n-TURN-1; u++)
      qot += qm1[jindx[u]+k]*qm1[jindx[n]+(u+1)];
    qm2[k] = qot;
   }

  for(p = 1; p < n; p++){
    for(q = p + TURN + 1; q <= n; q++){
      int type;
      /* 1. get exterior hairpin contribution  */
      u = n-q + p-1;
      if (u<TURN) continue;
      type = ptype[iindx[p]-q];
      if (!type) continue;
       /* cause we want to calc the exterior loops, we need the reversed pair type from now on  */
      type=rtype[type];

      char loopseq[10];
      if (u<7){
        strcpy(loopseq , sequence+q-1);
        strncat(loopseq, sequence, p);
      }
      qho += (((type==3)||(type==4))&&no_closingGU) ? 0. : qb[iindx[p]-q] * exp_E_Hairpin(u, type, S1[q+1], S1[p-1],  loopseq, pf_params) * scale[u];

      /* 2. exterior interior loops, i "define" the (k,l) pair as "outer pair"  */
      /* so "outer type" is rtype[type[k,l]] and inner type is type[p,q]        */
      qot = 0.;
      for(k=q+1; k < n; k++){
        int ln1, lstart;
        ln1 = k - q - 1;
        if(ln1+p-1>MAXLOOP) break;
        lstart = ln1+p-1+n-MAXLOOP;
        if(lstart<k+TURN+1) lstart = k + TURN + 1;
        for(l=lstart;l <= n; l++){
          int ln2, type2;
          ln2 = (p - 1) + (n - l);

          if((ln1+ln2) > MAXLOOP) continue;

          type2 = ptype[iindx[k]-l];
          if(!type2) continue;
          qio += qb[iindx[p]-q] * qb[iindx[k]-l] * exp_E_IntLoop(ln2, ln1, rtype[type2], type, S1[l+1], S1[k-1], S1[p-1], S1[q+1], pf_params) * scale[ln1+ln2];
        }
      } /* end of kl double loop */
    }
  } /* end of pq double loop */

  /* 3. Multiloops  */
  for(k=TURN+2; k<n-2*TURN-3; k++)
    qmo += qm[iindx[1]-k] * qm2[k+1] * expMLclosing;

  /* add an additional pf of 1.0 to take the open chain into account too           */
  qo = qho + qio + qmo + 1.0*scale[n];
}

/* calculate base pairing probs */
PUBLIC void pf_create_bppm(const char *sequence, char *structure){
  int n, i,j,k,l, ij, kl, ii,ll, type, type_2, tt, ov=0;
  FLT_OR_DBL  temp, Qmax=0, prm_MLb;
  FLT_OR_DBL  prmt,prmt1;
  FLT_OR_DBL  *tmp;
  FLT_OR_DBL  tmp2;
  FLT_OR_DBL  expMLclosing = pf_params->expMLclosing;
  double      max_real;

  max_real = (sizeof(FLT_OR_DBL) == sizeof(float)) ? FLT_MAX : DBL_MAX;

  if((S != NULL) && (S1 != NULL)){
    n = S[0];
    Qmax=0;

    for (k=1; k<=n; k++) {
      q1k[k] = q[iindx[1] - k];
      qln[k] = q[iindx[k] -n];
    }
    q1k[0] = 1.0;
    qln[n+1] = 1.0;

/*  pr = q; */     /* recycling */


    /* 1. exterior pair i,j and initialization of pr array */
    if(circular){
      for (i=1; i<=n; i++) {
        for (j=i; j<=MIN2(i+TURN,n); j++)
          probs[iindx[i]-j] = 0;
        for (j=i+TURN+1; j<=n; j++) {
          ij = iindx[i]-j;
          type = ptype[ij];
          if (type&&(qb[ij]>0.)) {
            probs[ij] = 1./qo;
            int rt = rtype[type];

            /* 1.1. Exterior Hairpin Contribution */
            int u = i + n - j -1;
            /* get the loop sequence */
            char loopseq[10];
            if (u<7){
              strcpy(loopseq , sequence+j-1);
              strncat(loopseq, sequence, i);
            }
            tmp2 = exp_E_Hairpin(u, rt, S1[j+1], S1[i-1], loopseq, pf_params) * scale[u];

            /* 1.2. Exterior Interior Loop Contribution                    */
            /* 1.2.1. i,j  delimtis the "left" part of the interior loop    */
            /* (j,i) is "outer pair"                                                */
            for(k=1; k < i-TURN-1; k++){
              int ln1, lstart;
              ln1 = k + n - j - 1;
              if(ln1>MAXLOOP) break;
              lstart = ln1+i-1-MAXLOOP;
              if(lstart<k+TURN+1) lstart = k + TURN + 1;
              for(l=lstart; l < i; l++){
                int ln2, type_2;
                type_2 = ptype[iindx[k]-l];
                if (type_2==0) continue;
                ln2 = i - l - 1;
                if(ln1+ln2>MAXLOOP) continue;
                tmp2 += qb[iindx[k] - l]*exp_E_IntLoop(ln1, ln2, rt, rtype[type_2], S1[j+1], S1[i-1], S1[k-1], S1[l+1], pf_params) * scale[ln1 + ln2];
              }
            }
            /* 1.2.2. i,j  delimtis the "right" part of the interior loop  */
            for(k=j+1; k < n-TURN; k++){
              int ln1, lstart;
              ln1 = k - j - 1;
              if((ln1 + i - 1)>MAXLOOP) break;
              lstart = ln1+i-1+n-MAXLOOP;
              if(lstart<k+TURN+1) lstart = k + TURN + 1;
              for(l=lstart; l <= n; l++){
                int ln2, type_2;
                type_2 = ptype[iindx[k]-l];
                if (type_2==0) continue;
                ln2 = i - 1 + n - l;
                if(ln1+ln2>MAXLOOP) continue;
                tmp2 += qb[iindx[k] - l]*exp_E_IntLoop(ln2, ln1, rtype[type_2], rt, S1[l+1], S1[k-1], S1[i-1], S1[j+1], pf_params) * scale[ln1 + ln2];
              }
            }
            /* 1.3 Exterior multiloop decomposition */
            /* 1.3.1 Middle part                    */
            if((i>TURN+2) && (j<n-TURN-1))
              tmp2 += qm[iindx[1]-i+1] * qm[iindx[j+1]-n] * expMLclosing * exp_E_MLstem(type, S1[i-1], S1[j+1], pf_params);

            /* 1.3.2 Left part                      */
            for(k=TURN+2; k < i-TURN-2; k++)
              tmp2 += qm[iindx[1]-k] * qm1[jindx[i-1]+k+1] * expMLbase[n-j] * expMLclosing * exp_E_MLstem(type, S1[i-1], S1[j+1], pf_params);

            /* 1.3.3 Right part                      */
            for(k=j+TURN+2; k < n-TURN-1;k++)
              tmp2 += qm[iindx[j+1]-k] * qm1[jindx[n]+k+1] * expMLbase[i-1] * expMLclosing * exp_E_MLstem(type, S1[i-1], S1[j+1], pf_params);

            /* all exterior loop decompositions for pair i,j done  */
            probs[ij] *= tmp2;

          }
          else probs[ij] = 0;
        }
      }
    } /* end if(circular)  */
    else {
      for (i=1; i<=n; i++) {
        for (j=i; j<=MIN2(i+TURN,n); j++) probs[iindx[i]-j] = 0;
        for (j=i+TURN+1; j<=n; j++) {
          ij = iindx[i]-j;
          type = ptype[ij];
          if (type&&(qb[ij]>0.)) {
            probs[ij] = q1k[i-1]*qln[j+1]/q1k[n];
            probs[ij] *= exp_E_ExtLoop(type, (i>1) ? S1[i-1] : -1, (j<n) ? S1[j+1] : -1, pf_params);
          } else
            probs[ij] = 0;
        }
      }
    } /* end if(!circular)  */

    for (l=n; l>TURN+1; l--) {

      /* 2. bonding k,l as substem of 2:loop enclosed by i,j */
      for (k=1; k<l-TURN; k++) {
        kl = iindx[k]-l;
        type_2 = ptype[kl]; type_2 = rtype[type_2];
        if (qb[kl]==0) continue;

        for (i=MAX2(1,k-MAXLOOP-1); i<=k-1; i++)
          for (j=l+1; j<=MIN2(l+ MAXLOOP -k+i+2,n); j++) {
            ij = iindx[i] - j;
            type = ptype[ij];
            if ((probs[ij]>0)) {
              /* add *scale[u1+u2+2] */
              probs[kl] += probs[ij] * (scale[k-i+j-l] *
                        exp_E_IntLoop(k-i-1, j-l-1, type, type_2,
                        S1[i+1], S1[j-1], S1[k-1], S1[l+1], pf_params));
            }
          }
      }
      /* 3. bonding k,l as substem of multi-loop enclosed by i,j */
      prm_MLb = 0.;
      if (l<n) for (k=2; k<l-TURN; k++) {
        i = k-1;
        prmt = prmt1 = 0.0;

        ii = iindx[i];     /* ii-j=[i,j]     */
        ll = iindx[l+1];   /* ll-j=[l+1,j-1] */
        tt = ptype[ii-(l+1)]; tt=rtype[tt];
        prmt1 = probs[ii-(l+1)] * expMLclosing * exp_E_MLstem(tt, S1[l], S1[i+1], pf_params);

        for (j=l+2; j<=n; j++) {
          tt = ptype[ii-j]; tt = rtype[tt];
          prmt += probs[ii-j] * exp_E_MLstem(tt, S1[j-1], S1[i+1], pf_params) * qm[ll-(j-1)];
        }
        kl = iindx[k]-l;
        tt = ptype[kl];
        prmt *= expMLclosing;
        prml[ i] = prmt;
        prm_l[i] = prm_l1[i]*expMLbase[1]+prmt1;

        prm_MLb = prm_MLb*expMLbase[1] + prml[i];
        /* same as:    prm_MLb = 0;
           for (i=1; i<=k-1; i++) prm_MLb += prml[i]*expMLbase[k-i-1]; */

        prml[i] = prml[ i] + prm_l[i];

        if (qb[kl] == 0.) continue;

        temp = prm_MLb;

        for (i=1;i<=k-2; i++)
          temp += prml[i]*qm[iindx[i+1] - (k-1)];

        temp    *= exp_E_MLstem(tt, (k>1) ? S1[k-1] : -1, (l<n) ? S1[l+1] : -1, pf_params) * scale[2];
        probs[kl]  += temp;

        if (probs[kl]>Qmax) {
          Qmax = probs[kl];
          if (Qmax>max_real/10.)
            fprintf(stderr, "P close to overflow: %d %d %g %g\n",
              i, j, probs[kl], qb[kl]);
        }
        if (probs[kl]>=max_real) {
          ov++;
          probs[kl]=FLT_MAX;
        }

      } /* end for (k=..) */
      tmp = prm_l1; prm_l1=prm_l; prm_l=tmp;

    }  /* end for (l=..)   */

    for (i=1; i<=n; i++)
      for (j=i+TURN+1; j<=n; j++) {
        ij = iindx[i]-j;
        probs[ij] *= qb[ij];
      }

    if (structure!=NULL)
      bppm_to_structure(structure, probs, n);
    if (ov>0) fprintf(stderr, "%d overflows occurred while backtracking;\n"
        "you might try a smaller pf_scale than %g\n",
        ov, pf_scale);
  } /* end if((S != NULL) && (S1 != NULL))  */
  else
    nrerror("bppm calculations have to be done after calling forward recursion\n");
  return;
}

PRIVATE void scale_pf_params(unsigned int length){
  unsigned int i;
  double  kT;

  if(pf_params) free(pf_params);
  pf_params = get_scaled_pf_parameters();

  kT = pf_params->kT;   /* kT in cal/mol  */

   /* scaling factors (to avoid overflows) */
  if (pf_scale == -1) { /* mean energy for random sequences: 184.3*length cal */
    pf_scale = exp(-(-185+(pf_params->temperature-37.)*7.27)/kT);
    if (pf_scale<1) pf_scale=1;
  }
  scale[0] = 1.;
  scale[1] = 1./pf_scale;
  expMLbase[0] = 1;
  expMLbase[1] = pf_params->expMLbase/pf_scale;
  for (i=2; i<=length; i++) {
    scale[i] = scale[i/2]*scale[i-(i/2)];
    expMLbase[i] = pow(pf_params->expMLbase, (double)i) * scale[i];
  }
}

/*----------------------------------------------------------------------*/


/*----------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/

PUBLIC void update_pf_params(int length){
#ifdef _OPENMP
  make_pair_matrix(); /* is this really necessary? */
  scale_pf_params((unsigned) length);
#else
  if (length>init_length) init_partfunc(length);  /* init not update */
  else {
    make_pair_matrix();
    scale_pf_params((unsigned) length);
  }
#endif
}

/*---------------------------------------------------------------------------*/

PUBLIC char bppm_symbol(const float *x){
  if( x[0] > 0.667 )  return '.';
  if( x[1] > 0.667 )  return '(';
  if( x[2] > 0.667 )  return ')';
  if( (x[1]+x[2]) > x[0] ) {
    if( (x[1]/(x[1]+x[2])) > 0.667) return '{';
    if( (x[2]/(x[1]+x[2])) > 0.667) return '}';
    else return '|';
  }
  if( x[0] > (x[1]+x[2]) ) return ',';
  return ':';
}

PUBLIC void bppm_to_structure(char *structure, FLT_OR_DBL *pr, unsigned int length){
  int    i, j;
  int   *iindx = get_iindx(length);
  float  P[3];   /* P[][0] unpaired, P[][1] upstream p, P[][2] downstream p */

  for( j=1; j<=length; j++ ) {
    P[0] = 1.0;
    P[1] = P[2] = 0.0;
    for( i=1; i<j; i++) {
      P[2] += pr[iindx[i]-j];    /* j is paired downstream */
      P[0] -= pr[iindx[i]-j];    /* j is unpaired */
    }
    for( i=j+1; i<=length; i++ ) {
      P[1] += pr[iindx[j]-i];    /* j is paired upstream */
      P[0] -= pr[iindx[j]-i];    /* j is unpaired */
    }
    structure[j-1] = bppm_symbol(P);
  }
  structure[length] = '\0';
  free(iindx);
}


/*---------------------------------------------------------------------------*/
PRIVATE void make_ptypes(const short *S, const char *structure){
  int n,i,j,k,l;

  n=S[0];
  for (k=1; k<n-TURN; k++)
    for (l=1; l<=2; l++) {
      int type,ntype=0,otype=0;
      i=k; j = i+TURN+l; if (j>n) continue;
      type = pair[S[i]][S[j]];
      while ((i>=1)&&(j<=n)) {
        if ((i>1)&&(j<n)) ntype = pair[S[i-1]][S[j+1]];
        if (noLonelyPairs && (!otype) && (!ntype))
          type = 0; /* i.j can only form isolated pairs */
        qb[iindx[i]-j] = 0.;
        ptype[iindx[i]-j] = (char) type;
        otype =  type;
        type  = ntype;
        i--; j++;
      }
    }

  if (fold_constrained && (structure != NULL))
    constrain_ptypes(structure, (unsigned int)n, ptype, NULL, TURN, 1);
}

/*
  stochastic backtracking in pf_fold arrays
  returns random structure S with Boltzman probabilty
  p(S) = exp(-E(S)/kT)/Z
*/
char *pbacktrack(char *seq){
  double r, qt;
  int i,j,n, start;
  sequence = seq;
  n = strlen(sequence);

  if (init_length<1)
    nrerror("can't backtrack without pf arrays.\n"
            "Call pf_fold() before pbacktrack()");
  pstruc = space((n+1)*sizeof(char));

  for (i=0; i<n; i++) pstruc[i] = '.';

  start = 1;
  while (start<n) {
  /* find i position of first pair */
    for (i=start; i<n; i++) {
      r = urn() * qln[i];
      if (r > qln[i+1]*scale[1])  break; /* i is paired */
    }
    if (i>=n) break; /* no more pairs */
    /* now find the pairing partner j */
    r = urn() * (qln[i] - qln[i+1]*scale[1]);
    for (qt=0, j=i+1; j<=n; j++) {
      int type;
      type = ptype[iindx[i]-j];
      if (type) {
        double qkl;
        qkl = qb[iindx[i]-j];
        if (j<n) qkl *= qln[j+1];
        qkl *=  exp_E_ExtLoop(type, (i>1) ? S1[i-1] : -1, (j<n) ? S1[j+1] : -1, pf_params);
        qt += qkl;
        if (qt > r) break; /* j is paired */
      }
    }
    if (j==n+1) nrerror("backtracking failed in ext loop");
    start = j+1;
    backtrack(i,j);
  }

  return pstruc;
}
char *pbacktrack_circ(char *seq){
  double r, qt;
  int i, j, k, l, n;
  FLT_OR_DBL  expMLclosing      = pf_params->expMLclosing;

  sequence = seq;
  n = strlen(sequence);

  if (init_length<1)
    nrerror("can't backtrack without pf arrays.\n"
      "Call pf_circ_fold() before pbacktrack_circ()");
  pstruc = space((n+1)*sizeof(char));

  /* initialize pstruct with single bases  */
  for (i=0; i<n; i++) pstruc[i] = '.';

  qt = 1.0*scale[n];
  r = urn() * qo;

  /* open chain? */
  if(qt > r) return pstruc;

  for(i=1; (i < n); i++){
    for(j=i+TURN+1;(j<=n); j++){

      int type, u;
      /* 1. first check, wether we can do a hairpin loop  */
      u = n-j + i-1;
      if (u<TURN) continue;

      type = ptype[iindx[i]-j];
      if (!type) continue;

      type=rtype[type];

      char loopseq[10];
      if (u<7){
        strcpy(loopseq , sequence+j-1);
        strncat(loopseq, sequence, i);
      }

      qt += qb[iindx[i]-j] * exp_E_Hairpin(u, type, S1[j+1], S1[i-1],  loopseq, pf_params) * scale[u];
      /* found a hairpin? so backtrack in the enclosed part and we're done  */
      if(qt>r){ backtrack(i,j); return pstruc;}

      /* 2. search for (k,l) with which we can close an interior loop  */
      for(k=j+1; (k < n); k++){
        int ln1, lstart;
        ln1 = k - j - 1;
        if(ln1+i-1>MAXLOOP) break;

        lstart = ln1+i-1+n-MAXLOOP;
        if(lstart<k+TURN+1) lstart = k + TURN + 1;
        for(l=lstart; (l <= n); l++){
            int ln2, type2;
            ln2 = (i - 1) + (n - l);
            if((ln1+ln2) > MAXLOOP) continue;

            type2 = ptype[iindx[k]-l];
            if(!type) continue;
            type2 = rtype[type2];
            qt += qb[iindx[i]-j] * qb[iindx[k]-l] * exp_E_IntLoop(ln2, ln1, type2, type, S1[l+1], S1[k-1], S1[i-1], S1[j+1], pf_params) * scale[ln1 + ln2];
            /* found an exterior interior loop? also this time, we can go straight  */
            /* forward and backtracking the both enclosed parts and we're done      */
            if(qt>r){ backtrack(i,j); backtrack(k,l); return pstruc;}
        }
      } /* end of kl double loop */
    }
  } /* end of ij double loop  */
  {
    /* as we reach this part, we have to search for our barrier between qm and qm2  */
    qt = 0.;
    r = urn()*qmo;
    for(k=TURN+2; k<n-2*TURN-3; k++){
      qt += qm[iindx[1]-k] * qm2[k+1] * expMLclosing;
      /* backtrack in qm and qm2 if we've found a valid barrier k  */
      if(qt>r){ backtrack_qm(1,k); backtrack_qm2(k+1,n); return pstruc;}
    }
  }
  /* if we reach the actual end of this function, an error has occured  */
  /* cause we HAVE TO find an exterior loop or an open chain!!!         */
  nrerror("backtracking failed in exterior loop");
  return pstruc;
}

PRIVATE void backtrack_qm(int i, int j){
  /* divide multiloop into qm and qm1  */
  double qmt, r;
  int k;
  while(j>i){
    /* now backtrack  [i ... j] in qm[] */
    r = urn() * qm[iindx[i] - j];
    qmt = qm1[jindx[j]+i]; k=i;
    if(qmt<r)
      for(k=i+1; k<=j; k++){
        qmt += (qm[iindx[i]-(k-1)]+expMLbase[k-i])*qm1[jindx[j]+k];
        if(qmt >= r) break;
      }
    if(k>j) nrerror("backtrack failed in qm");

    backtrack_qm1(k,j);

    if(k<i+TURN) break; /* no more pairs */
    r = urn() * (qm[iindx[i]-(k-1)] + expMLbase[k-i]);
    if(expMLbase[k-i] >= r) break; /* no more pairs */
    j = k-1;
  }
}

PRIVATE void backtrack_qm1(int i,int j){
  /* i is paired to l, i<l<j; backtrack in qm1 to find l */
  int ii, l, type;
  double qt, r;
  r = urn() * qm1[jindx[j]+i];
  ii = iindx[i];
  for (qt=0., l=i+TURN+1; l<=j; l++) {
    type = ptype[ii-l];
    if (type)
      qt +=  qb[ii-l] * exp_E_MLstem(type, S1[i-1], S1[l+1], pf_params) * expMLbase[j-l];
    if (qt>=r) break;
  }
  if (l>j) nrerror("backtrack failed in qm1");
  backtrack(i,l);
}

PRIVATE void backtrack_qm2(int k, int n){
  double qom2t, r;
  int u;
  r= urn()*qm2[k];
  /* we have to search for our barrier u between qm1 and qm1  */
  for (qom2t = 0.,u=k+TURN+1; u<n-TURN-1; u++){
    qom2t += qm1[jindx[u]+k]*qm1[jindx[n]+(u+1)];
    if(qom2t > r) break;
  }
  if(u==n-TURN) nrerror("backtrack failed in qm2");
  backtrack_qm1(k,u);
  backtrack_qm1(u+1,n);
}

PRIVATE void backtrack(int i, int j){
  do {
    double r, qbt1;
    int k, l, type, u, u1;

    pstruc[i-1] = '('; pstruc[j-1] = ')';

    r = urn() * qb[iindx[i]-j];
    type = ptype[iindx[i]-j];
    u = j-i-1;
    /*hairpin contribution*/
    if (((type==3)||(type==4))&&no_closingGU) qbt1 = 0;
    else
      qbt1 = exp_E_Hairpin(u, type, S1[i+1], S1[j-1], sequence+i-1, pf_params)*
        scale[u+2]; /* add scale[u+2] */

    if (qbt1>=r) return; /* found the hairpin we're done */

    for (k=i+1; k<=MIN2(i+MAXLOOP+1,j-TURN-2); k++) {
      u1 = k-i-1;
      for (l=MAX2(k+TURN+1,j-1-MAXLOOP+u1); l<j; l++) {
        int type_2;
        type_2 = ptype[iindx[k]-l];
        if (type_2) {
          type_2 = rtype[type_2];
          /* add *scale[u1+u2+2] */
          qbt1 += qb[iindx[k]-l] * (scale[u1+j-l+1] *
            exp_E_IntLoop(u1, j-l-1, type, type_2,
                          S1[i+1], S1[j-1], S1[k-1], S1[l+1], pf_params));
        }
        if (qbt1 > r) break;
      }
      if (qbt1 > r) break;
    }
    if (l<j) {
      i=k; j=l;
    }
    else break;
  } while (1);

  /* backtrack in multi-loop */
  {
    double r, qt;
    int k, ii, jj;

    i++; j--;
    /* find the first split index */
    ii = iindx[i]; /* ii-j=[i,j] */
    jj = jindx[j]; /* jj+i=[j,i] */
    for (qt=0., k=i+1; k<j; k++) qt += qm[ii-(k-1)]*qm1[jj+k];
    r = urn() * qt;
    for (qt=0., k=i+1; k<j; k++) {
      qt += qm[ii-(k-1)]*qm1[jj+k];
      if (qt>=r) break;
    }
    if (k>=j) nrerror("backtrack failed, can't find split index ");

    backtrack_qm1(k, j);

    j = k-1;
    backtrack_qm(i,j);
  }
}

PUBLIC void assign_plist_from_pr(plist **pl, double *probs, int length, double cut_off){
  int i, j, n, count, *index;
  count = 0;
  n     = 2;

  index = get_iindx(length);

  /* first guess of the size needed for pl */
  *pl = (plist *)space(n*length*sizeof(plist));

  for (i=1; i<length; i++) {
    for (j=i+1; j<=length; j++) {
      /* skip all entries below the cutoff */
      if (probs[index[i]-j] < cut_off) continue;
      /* do we need to allocate more memory? */
      if (count == n * length - 1){
        n *= 2;
        *pl = (plist *)xrealloc(*pl, n * length * sizeof(plist));
      }
      (*pl)[count].i    = i;
      (*pl)[count].j    = j;
      (*pl)[count++].p  = probs[index[i] - j];
    }
  }
  /* mark the end of pl */
  (*pl)[count].i   = 0;
  (*pl)[count].j   = 0;
  (*pl)[count++].p = 0.;
  /* shrink memory to actual size needed */
  *pl = (plist *)xrealloc(*pl, count * sizeof(plist));
  free(index);
}

/* this function is a threadsafe replacement for centroid() */
PUBLIC char *get_centroid_struct_pl(int length, double *dist, plist *pl){
  /* compute the centroid structure of the ensemble, i.e. the strutcure
     with the minimal average distance to all other structures
     <d(S)> = \sum_{(i,j) \in S} (1-p_{ij}) + \sum_{(i,j) \notin S} p_{ij}
     Thus, the centroid is simply the structure containing all pairs with
     p_ij>0.5 */
  int i;
  char *centroid;

  if (pl==NULL)
    nrerror("get_centroid_struct: pl==NULL!");

  *dist = 0.;
  centroid = (char *) space((length+1)*sizeof(char));
  for (i=0; i<length; i++) centroid[i]='.';
  for (i=0; pl[i].i>0; i++){
    if ((pl[i].p)>0.5) {
      centroid[pl[i].i-1] = '(';
      centroid[pl[i].j-1] = ')';
      *dist += (1-pl[i].p);
    } else
      *dist += pl[i].p;
  }
  centroid[length] = '\0';
  return centroid;
}

/* this function is a threadsafe replacement for centroid() */
PUBLIC char *get_centroid_struct_pr(int length, double *dist, double *probs){
  /* compute the centroid structure of the ensemble, i.e. the strutcure
     with the minimal average distance to all other structures
     <d(S)> = \sum_{(i,j) \in S} (1-p_{ij}) + \sum_{(i,j) \notin S} p_{ij}
     Thus, the centroid is simply the structure containing all pairs with
     p_ij>0.5 */
  int i,j;
  double p;
  char  *centroid;
  int   *my_iindx = get_iindx(length);

  if (probs == NULL)
    nrerror("get_centroid_struct_pr: probs==NULL!");

  *dist = 0.;
  centroid = (char *) space((length+1)*sizeof(char));
  for (i=0; i<length; i++) centroid[i]='.';
  for (i=1; i<=length; i++)
    for (j=i+TURN+1; j<=length; j++) {
      if ((p=probs[my_iindx[i]-j])>0.5) {
        centroid[i-1] = '(';
        centroid[j-1] = ')';
        *dist += (1-p);
      } else
        *dist += p;
    }
  free(my_iindx);
  centroid[length] = '\0';
  return centroid;
}

PUBLIC plist *stackProb(double cutoff){
  plist *pl;
  int i,j,plsize=256;
  int length, num = 0;
  if (pr==NULL)
    nrerror("pr==NULL. You need to call pf_fold() before stackProb()");

  pl = (plist *) space(plsize*sizeof(plist));
  length = S[0];
  for (i=1; i<length; i++)
    for (j=i+TURN+3; j<=length; j++) {
      double p;
      if ((p=pr[iindx[i]-j])<cutoff) continue;
      if (qb[iindx[i+1]-(j-1)]<FLT_MIN) continue;
      p *= qb[iindx[i+1]-(j-1)]/qb[iindx[i]-j];
      p *= exp_E_IntLoop(0,0,ptype[iindx[i]-j],rtype[ptype[iindx[i+1]-(j-1)]],
                         0,0,0,0, pf_params)*scale[2];/* add *scale[u1+u2+2] */
      if (p>cutoff) {
        pl[num].i = i;
        pl[num].j = j;
        pl[num++].p = p;
        if (num>=plsize) {
          plsize *= 2;
          pl = xrealloc(pl, plsize*sizeof(plist));
        }
      }
    }
  pl[num].i=0;
  return pl;
}

/*-------------------------------------------------------------------------*/
/* make arrays used for pf_fold available to other routines */
PUBLIC int get_pf_arrays(short **S_p, short **S1_p, char **ptype_p, FLT_OR_DBL **qb_p, FLT_OR_DBL **qm_p, FLT_OR_DBL **q1k_p, FLT_OR_DBL **qln_p){
  if(qb == NULL) return(0); /* check if pf_fold() has been called */
  *S_p = S; *S1_p = S1; *ptype_p = ptype;
  *qb_p = qb; *qm_p = qm;
  *q1k_p = q1k; *qln_p = qln;
  return(1); /* success */
}

PUBLIC double mean_bp_distance(int length){
  return mean_bp_distance_pr(length, probs);
}

PUBLIC double mean_bp_distance_pr(int length, double *pr){
  /* compute the mean base pair distance in the thermodynamic ensemble */
  /* <d> = \sum_{a,b} p_a p_b d(S_a,S_b)
     this can be computed from the pair probs p_ij as
     <d> = \sum_{ij} p_{ij}(1-p_{ij}) */
  int i,j;
  double d=0;
  int *my_iindx = get_iindx((unsigned int) length);

  if (pr==NULL)
    nrerror("pr==NULL. You need to supply a valid probability matrix for mean_bp_distance_pr()");

  for (i=1; i<=length; i++)
    for (j=i+TURN+1; j<=length; j++)
      d += pr[my_iindx[i]-j] * (1-pr[my_iindx[i]-j]);

  free(my_iindx);
  return 2*d;
}

PUBLIC FLT_OR_DBL *export_bppm(void){
  return probs;
}

/*###########################################*/
/*# deprecated functions below              #*/
/*###########################################*/

/* this function is deprecated as it is not threadsafe */
PUBLIC char *centroid(int length, double *dist) {
  /* compute the centroid structure of the ensemble, i.e. the strutcure
     with the minimal average distance to all other structures
     <d(S)> = \sum_{(i,j) \in S} (1-p_{ij}) + \sum_{(i,j) \notin S} p_{ij}
     Thus, the centroid is simply the structure containing all pairs with
     p_ij>0.5 */
  int i,j;
  double p;
  char *centroid;

  if (pr==NULL)
    nrerror("pr==NULL. You need to call pf_fold() before centroid()");

  *dist = 0.;
  centroid = (char *) space((length+1)*sizeof(char));
  for (i=0; i<length; i++) centroid[i]='.';
  for (i=1; i<=length; i++)
    for (j=i+TURN+1; j<=length; j++) {
      if ((p=pr[iindx[i]-j])>0.5) {
        centroid[i-1] = '(';
        centroid[j-1] = ')';
        *dist += (1-p);
      } else
        *dist += p;
    }
  return centroid;
}


/* This function is deprecated as it uses the global array pr for calculations */
PUBLIC double mean_bp_dist(int length) {
  /* compute the mean base pair distance in the thermodynamic ensemble */
  /* <d> = \sum_{a,b} p_a p_b d(S_a,S_b)
     this can be computed from the pair probs p_ij as
     <d> = \sum_{ij} p_{ij}(1-p_{ij}) */
  int i,j;
  double d=0;

  if (pr==NULL)
    nrerror("pr==NULL. You need to call pf_fold() before mean_bp_dist()");

  for (i=1; i<=length; i++)
    for (j=i+TURN+1; j<=length; j++)
      d += pr[iindx[i]-j] * (1-pr[iindx[i]-j]);
  return 2*d;
}

/*----------------------------------------------------------------------*/
PUBLIC double expHairpinEnergy(int u, int type, short si1, short sj1,
                                const char *string) {
/* compute Boltzmann weight of a hairpin loop, multiply by scale[u+2] */
  double q, kT;
  kT = pf_params->kT;   /* kT in cal/mol  */
  if(u <= 30)
    q = pf_params->exphairpin[u];
  else
    q = pf_params->exphairpin[30] * exp( -(pf_params->lxc*log( u/30.))*10./kT);
  if ((tetra_loop)&&(u==4)) {
    char tl[7]={0}, *ts;
    strncpy(tl, string, 6);
    if ((ts=strstr(pf_params->Tetraloops, tl)))
      return (pf_params->exptetra[(ts-pf_params->Tetraloops)/7]);
  }
  if ((tetra_loop)&&(u==6)) {
    char tl[9]={0}, *ts;
    strncpy(tl, string, 6);
    if ((ts=strstr(pf_params->Hexaloops, tl)))
      return  (pf_params->exphex[(ts-pf_params->Hexaloops)/9]);
  }
  if (u==3) {
    char tl[6]={0}, *ts;
    strncpy(tl, string, 5);
    if ((ts=strstr(pf_params->Triloops, tl)))
      return (pf_params->exptri[(ts-pf_params->Triloops)/6]);
    if (type>2)
      q *= pf_params->expTermAU;
  }
  else /* no mismatches for tri-loops */
    q *= pf_params->expmismatchH[type][si1][sj1];

  return q;
}
PUBLIC double expLoopEnergy(int u1, int u2, int type, int type2,
                             short si1, short sj1, short sp1, short sq1) {
/* compute Boltzmann weight of interior loop,
   multiply by scale[u1+u2+2] for scaling */
  double z=0;
  int no_close = 0;

  if ((no_closingGU) && ((type2==3)||(type2==4)||(type==2)||(type==4)))
    no_close = 1;

  if ((u1==0) && (u2==0)) /* stack */
    z = pf_params->expstack[type][type2];
  else if (no_close==0) {
    if ((u1==0)||(u2==0)) { /* bulge */
      int u;
      u = (u1==0)?u2:u1;
      z = pf_params->expbulge[u];
      if (u2+u1==1) z *= pf_params->expstack[type][type2];
      else {
        if (type>2) z *= pf_params->expTermAU;
        if (type2>2) z *= pf_params->expTermAU;
      }
    }
    else {     /* interior loop */
      if (u1+u2==2) /* size 2 is special */
        z = pf_params->expint11[type][type2][si1][sj1];
      else if ((u1==1) && (u2==2))
        z = pf_params->expint21[type][type2][si1][sq1][sj1];
      else if ((u1==2) && (u2==1))
        z = pf_params->expint21[type2][type][sq1][si1][sp1];
      else if ((u1==2) && (u2==2))
        z = pf_params->expint22[type][type2][si1][sp1][sq1][sj1];
      else if (((u1==2)&&(u2==3))||((u1==3)&&(u2==2))){ /*2-3 is special*/
        z = pf_params->expinternal[5]*
          pf_params->expmismatch23I[type][si1][sj1]*
          pf_params->expmismatch23I[type2][sq1][sp1];
        z *= pf_params->expninio[2][1];
      }
      else if ((u1==1)||(u2==1)) {  /*1-n is special*/
        z = pf_params->expinternal[u1+u2]*
          pf_params->expmismatch1nI[type][si1][sj1]*
          pf_params->expmismatch1nI[type2][sq1][sp1];
        z *= pf_params->expninio[2][abs(u1-u2)];
      }
      else {
        z = pf_params->expinternal[u1+u2]*
          pf_params->expmismatchI[type][si1][sj1]*
          pf_params->expmismatchI[type2][sq1][sp1];
        z *= pf_params->expninio[2][abs(u1-u2)];
      }
    }
  }
  return z;
}

PUBLIC void init_pf_circ_fold(int length){
/* DO NOTHING */
}

PUBLIC void init_pf_fold(int length){
/* DO NOTHING */
}


