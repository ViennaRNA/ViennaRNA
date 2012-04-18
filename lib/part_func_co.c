/* Last changed Time-stamp: <2007-05-09 16:11:21 ivo> */
/*
                  partiton function for RNA secondary structures

                  Ivo L Hofacker
                  Stephan Bernhart
                  Vienna RNA package
*/
/*
  $Log: part_func_co.c,v $
  Revision 1.10  2007/05/10 17:27:01  ivo
  make sure the relative error eps is positive in newton iteration

  Revision 1.9  2006/05/10 15:12:11  ivo
  some compiler choked on  double semicolon after declaration

  Revision 1.8  2006/04/05 12:52:31  ivo
  Fix performance bug (O(n^4) loop)

  Revision 1.7  2006/01/19 11:30:04  ivo
  compute_probabilities should only look at one dimer at a time

  Revision 1.6  2006/01/18 12:55:40  ivo
  major cleanup of berni code
  fix bugs related to confusing which free energy is returned by co_pf_fold()

  Revision 1.5  2006/01/16 11:32:25  ivo
  small bug in multiloop pair probs

  Revision 1.4  2006/01/05 18:13:40  ivo
  update

  Revision 1.3  2006/01/04 15:14:29  ivo
  fix bug in concentration calculations

  Revision 1.2  2004/12/23 12:14:41  berni
  *** empty log message ***

  Revision 1.1  2004/12/22 10:46:17  berni

  Partition function Cofolding 0.9, Computation of concentrations.

  Revision 1.16  2003/08/04 09:14:09  ivo
  finish up stochastic backtracking

  Revision 1.15  2002/03/19 16:51:12  ivo
  more on stochastic backtracking (still incomplete)

  Revision 1.13  2001/11/16 17:30:04  ivo
  add stochastic backtracking (still incomplete)
*/

#include <config.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>    /* #defines FLT_MAX ... */
#include <limits.h>

#include "utils.h"
#include "energy_par.h"
#include "fold_vars.h"
#include "pair_mat.h"
#include "PS_dot.h"
#include "params.h"
#include "loop_energies.h"
#include "part_func.h"
#include "part_func_co.h"

#ifdef _OPENMP
#include <omp.h>
#endif


/*@unused@*/
PRIVATE char rcsid[] UNUSED = "$Id: part_func_co.c,v 1.10 2007/05/10 17:27:01 ivo Exp $";

#define ISOLATED  256.0
#undef TURN
#define TURN 0
#define SAME_STRAND(I,J) (((I)>=cut_point)||((J)<cut_point))

/* #define SAME_STRAND(I,J) (((J)<cut_point)||((I)>=cut_point2)||(((I)>=cut_point)&&((J)<cut_point2)))
 */

/*
#################################
# GLOBAL VARIABLES              #
#################################
*/
int     mirnatog      = 0;
double  F_monomer[2]  = {0,0}; /* free energies of the two monomers */

/*
#################################
# PRIVATE VARIABLES             #
#################################
*/
PRIVATE FLT_OR_DBL  *expMLbase=NULL;
PRIVATE FLT_OR_DBL  *q=NULL, *qb=NULL, *qm=NULL, *qm1=NULL, *qqm=NULL, *qqm1=NULL, *qq=NULL, *qq1=NULL;
PRIVATE FLT_OR_DBL  *prml=NULL, *prm_l=NULL, *prm_l1=NULL, *q1k=NULL, *qln=NULL, *probs=NULL;
PRIVATE FLT_OR_DBL  *scale=NULL;
PRIVATE pf_paramT   *pf_params = NULL;
PRIVATE char        *ptype=NULL; /* precomputed array of pair types */
PRIVATE int         *jindx=NULL;
PRIVATE int         *my_iindx=NULL;
PRIVATE int         init_length; /* length in last call to init_pf_fold() */
PRIVATE int         do_bppm = 1;             /* do backtracking per default */
PRIVATE short       *S=NULL, *S1=NULL;
PRIVATE char        *pstruc=NULL;
PRIVATE char        *sequence=NULL;
PRIVATE double      alpha = 1.0;
PRIVATE int         struct_constrained = 0;

#ifdef _OPENMP

/* NOTE: all variables are assumed to be uninitialized if they are declared as threadprivate
         thus we have to initialize them before usage by a seperate function!
         OR: use copyin in the PARALLEL directive!
         e.g.:
         #pragma omp parallel for copyin(pf_params)
*/
#pragma omp threadprivate(expMLbase, q, qb, qm, qm1, qqm, qqm1, qq, qq1, prml, prm_l, prm_l1, q1k, qln,\
                          scale, pf_params, ptype, jindx, my_iindx, init_length, S, S1, pstruc, sequence, probs, do_bppm, alpha, struct_constrained)

#endif


/*
#################################
# PRIVATE FUNCTION DECLARATIONS #
#################################
*/
PRIVATE void    init_partfunc_co(int length, pf_paramT *parameters);
PRIVATE void    pf_co(const char *sequence);
PRIVATE void    pf_co_bppm(const char *sequence, char *structure);
PRIVATE double  *Newton_Conc(double ZAB, double ZAA, double ZBB, double concA, double concB,double* ConcVec);
PRIVATE void    scale_pf_params(unsigned int length, pf_paramT *parameters);
PRIVATE void    get_arrays(unsigned int length);
PRIVATE void    make_ptypes(const short *S, const char *structure);
PRIVATE void    backtrack(int i, int j);


/*
#################################
# BEGIN OF FUNCTION DEFINITIONS #
#################################
*/

PRIVATE void init_partfunc_co(int length, pf_paramT *parameters){
  if (length<1) nrerror("init_pf_fold: length must be greater 0");

#ifdef _OPENMP
/* Explicitly turn off dynamic threads */
  omp_set_dynamic(0);
  free_co_pf_arrays(); /* free previous allocation */
#else
  if (init_length>0) free_co_pf_arrays(); /* free previous allocation */
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
  scale_pf_params((unsigned) length, parameters);
  init_length = length;
}

PRIVATE void get_arrays(unsigned int length){
  unsigned int size;

  if((length +1) >= (unsigned int)sqrt((double)INT_MAX))
    nrerror("get_arrays@part_func_co.c: sequence length exceeds addressable range");

  size      = sizeof(FLT_OR_DBL) * ((length+1)*(length+2)/2);
  q         = (FLT_OR_DBL *) space(size);
  qb        = (FLT_OR_DBL *) space(size);
  qm        = (FLT_OR_DBL *) space(size);
  probs     = (FLT_OR_DBL *) space(size);
  qm1       = (FLT_OR_DBL *) space(size);
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
  ptype     = (char *) space(sizeof(char)*((length+1)*(length+2)/2));
  my_iindx  = get_iindx(length);
  iindx     = get_iindx(length); /* for backward compatibility and Perl wrapper */
  jindx     = get_indx(length);
}

PUBLIC void free_co_pf_arrays(void){
  if(q)         free(q);
  if(qb)        free(qb);
  if(qm)        free(qm);
  if(qm1)       free(qm1);
  if(ptype)     free(ptype);
  if(qq)        free(qq);
  if(qq1)       free(qq1);
  if(qqm)       free(qqm);
  if(qqm1)      free(qqm1);
  if(q1k)       free(q1k);
  if(qln)       free(qln);
  if(prm_l)     free(prm_l);
  if(prm_l1)    free(prm_l1);
  if(prml)      free(prml);
  if(probs)     free(probs);
  if(expMLbase) free(expMLbase);
  if(scale)     free(scale);
  if(my_iindx)  free(my_iindx);
  if(iindx)     free(iindx); /* for backward compatibility and Perl wrapper */
  if(jindx)     free(jindx);
  if(S)         free(S);
  if(S1)        free(S1);

  init_length=0;
  q = qb = qm = qm1 = qq = qq1 = qqm = qqm1 = q1k = qln = prm_l = prm_l1 = prml = expMLbase = scale = probs = NULL;
  ptype = NULL;
  S = S1 = NULL;
  my_iindx = jindx = iindx = NULL;

#ifdef SUN4
  standard_arithmetic();
#else
#ifdef HP9
  fpsetfastmode(0);
#endif
#endif
}

/*-----------------------------------------------------------------*/
PUBLIC cofoldF co_pf_fold(char *sequence, char *structure){
  return co_pf_fold_par(sequence, structure, NULL, do_backtrack, fold_constrained);
}

PUBLIC cofoldF co_pf_fold_par(char *sequence,
                              char *structure,
                              pf_paramT *parameters,
                              int calculate_bppm,
                              int is_constrained){

  int         n;
  FLT_OR_DBL  Q;
  cofoldF     X;
  double      free_energy;

  n                   = (int) strlen(sequence);
  do_bppm             = calculate_bppm;
  struct_constrained  = is_constrained;

#ifdef _OPENMP
  /* always init everything since all global static variables are uninitialized when entering a thread */
  init_partfunc_co(n, parameters);
#else
  if(parameters) init_partfunc_co(n, parameters);
  else if (n > init_length) init_partfunc_co(n, parameters);
  else if (fabs(pf_params->temperature - temperature)>1e-6) update_co_pf_params_par(n, parameters);
#endif

 /* printf("mirnatog=%d\n",mirnatog); */

  if(S) free(S);
  S   = encode_sequence(sequence, 0);
  if(S1) free(S1);
  S1  = encode_sequence(sequence, 1);

  make_ptypes(S, structure);

  pf_co(sequence);

  if (backtrack_type=='C')      Q = qb[my_iindx[1]-n];
  else if (backtrack_type=='M') Q = qm[my_iindx[1]-n];
  else Q = q[my_iindx[1]-n];
  /* ensemble free energy in Kcal/mol */
  if (Q<=FLT_MIN) fprintf(stderr, "pf_scale too large\n");
  free_energy = (-log(Q)-n*log(pf_params->pf_scale))*pf_params->kT/1000.0;
  /* in case we abort because of floating point errors */
  if (n>1600) fprintf(stderr, "free energy = %8.2f\n", free_energy);
  /*probability of molecules being bound together*/


  /*Computation of "real" Partition function*/
  /*Need that for concentrations*/
  if (cut_point>0){
    double kT, pbound, QAB, QToT, Qzero;

    kT = pf_params->kT/1000.0;
    Qzero=q[my_iindx[1]-n];
    QAB=(q[my_iindx[1]-n]-q[my_iindx[1]-(cut_point-1)]*q[my_iindx[cut_point]-n])*pf_params->expDuplexInit;
    /*correction for symmetry*/
    if((n-(cut_point-1)*2)==0) {
      if ((strncmp(sequence, sequence+cut_point-1, cut_point-1))==0) {
        QAB/=2;
      }}

    QToT=q[my_iindx[1]-(cut_point-1)]*q[my_iindx[cut_point]-n]+QAB;
    pbound=1-(q[my_iindx[1]-(cut_point-1)]*q[my_iindx[cut_point]-n]/QToT);
     X.FAB  = -kT*(log(QToT)+n*log(pf_params->pf_scale));
    X.F0AB = -kT*(log(Qzero)+n*log(pf_params->pf_scale));
    X.FcAB = (QAB>1e-17) ? -kT*(log(QAB)+n*log(pf_params->pf_scale)) : 999;
    X.FA = -kT*(log(q[my_iindx[1]-(cut_point-1)]) + (cut_point-1)*log(pf_params->pf_scale));
    X.FB = -kT*(log(q[my_iindx[cut_point]-n]) + (n-cut_point+1)*log(pf_params->pf_scale));

    /* printf("QAB=%.9f\tQtot=%.9f\n",QAB/scale[n],QToT/scale[n]);*/
  }
  else {
    X.FA = X.FB = X.FAB = X.F0AB = free_energy;
    X.FcAB = 0;
  }

  /* backtracking to construct binding probabilities of pairs*/
  if(do_bppm){
    pf_co_bppm(sequence, structure);
    /*
    *  Backward compatibility:
    *  This block may be removed if deprecated functions
    *  relying on the global variable "pr" vanish from within the package!
    */
    pr = probs;
    /*
    {
      if(pr) free(pr);
      pr = (FLT_OR_DBL *) space(sizeof(FLT_OR_DBL) * ((n+1)*(n+2)/2));
      memcpy(pr, probs, sizeof(FLT_OR_DBL) * ((n+1)*(n+2)/2));
    }
    */
  }
  return X;
}

/* forward recursion of pf cofolding */
PRIVATE void pf_co(const char *sequence){
  int         n, i,j,k,l, ij, u,u1,ii, type, type_2, tt;
  FLT_OR_DBL  temp, Qmax=0;
  FLT_OR_DBL  qbt1, *tmp;
  FLT_OR_DBL  expMLclosing;
  double      max_real;
  int         noGUclosure = pf_params->model_details.noGUclosure;

  max_real = (sizeof(FLT_OR_DBL) == sizeof(float)) ? FLT_MAX : DBL_MAX;
  n = (int) strlen(sequence);

  expMLclosing = pf_params->expMLclosing;


  /*array initialization ; qb,qm,q
    qb,qm,q (i,j) are stored as ((n+1-i)*(n-i) div 2 + n+1-j */

  /* for (d=0; d<=TURN; d++) */
  for (i=1; i<=n/*-d*/; i++) {
      ij = my_iindx[i]-i;
      q[ij]=scale[1];
      qb[ij]=qm[ij]=0.0;
    }

  for (i=0; i<=n; i++)
    qq[i]=qq1[i]=qqm[i]=qqm1[i]=prm_l[i]=prm_l1[i]=prml[i]=0;

  for (j=TURN+2;j<=n; j++) {
    for (i=j-TURN-1; i>=1; i--) {
      /* construction of partition function of segment i,j*/
       /*firstly that given i bound to j : qb(i,j) */
      u = j-i-1; ij = my_iindx[i]-j;
      type = ptype[ij];
      qbt1=0;
      if (type!=0) {
        /*hairpin contribution*/
        if SAME_STRAND(i,j){
          if (((type==3)||(type==4))&&noGUclosure) qbt1 = 0;
          else
            qbt1 = exp_E_Hairpin(u, type, S1[i+1], S1[j-1], sequence+i-1, pf_params)*scale[u+2];

        }

        /* interior loops with interior pair k,l */
        for (k=i+1; k<=MIN2(i+MAXLOOP+1,j-TURN-2); k++) {
          u1 = k-i-1;
          for (l=MAX2(k+TURN+1,j-1-MAXLOOP+u1); l<j; l++) {
            if ((SAME_STRAND(i,k))&&(SAME_STRAND(l,j))){
              type_2 = ptype[my_iindx[k]-l];
              if (type_2) {
                type_2 = rtype[type_2];
                qbt1 += qb[my_iindx[k]-l] *
                  exp_E_IntLoop(u1, j-l-1, type, type_2,
                                S1[i+1], S1[j-1], S1[k-1], S1[l+1], pf_params)*scale[u1+j-l+1];
              }
            }
          }
        }
        /*multiple stem loop contribution*/
        ii = my_iindx[i+1]; /* ii-k=[i+1,k-1] */
        temp = 0.0;
        if (SAME_STRAND(i,i+1) && SAME_STRAND(j-1,j)) {
          for (k=i+2; k<=j-1; k++) {
            if (SAME_STRAND(k-1,k))
              temp += qm[ii-(k-1)]*qqm1[k];
          }
          tt = rtype[type];
          temp*=exp_E_MLstem(tt, S1[j-1], S1[i+1], pf_params)*scale[2];
          temp*=expMLclosing;
          qbt1 += temp;
        }
        /*qc contribution*/
        temp=0.0;
        if (!SAME_STRAND(i,j)){
          tt = rtype[type];
          temp=q[my_iindx[i+1]-(cut_point-1)]*q[my_iindx[cut_point]-(j-1)];
          if ((j==cut_point)&&(i==cut_point-1)) temp=scale[2];
          else if (i==cut_point-1) temp=q[my_iindx[cut_point]-(j-1)]*scale[1];
          else if (j==cut_point) temp=q[my_iindx[i+1]-(cut_point-1)]*scale[1];
          if (j>cut_point) temp*=scale[1];
          if (i<cut_point-1) temp*=scale[1];
          temp *= exp_E_ExtLoop(tt, SAME_STRAND(j-1,j) ? S1[j-1] : -1, SAME_STRAND(i,i+1) ? S1[i+1] : -1, pf_params);
          qbt1+=temp;
        }
        qb[ij] = qbt1;
      } /* end if (type!=0) */
      else qb[ij] = 0.0;
      /* construction of qqm matrix containing final stem
         contributions to multiple loop partition function
         from segment i,j */
      if (SAME_STRAND(j-1,j)) {
        qqm[i] = qqm1[i]*expMLbase[1];
      }
      else qqm[i]=0;
      if (type&&SAME_STRAND(i-1,i)&&SAME_STRAND(j,j+1)) {
        qbt1 = qb[ij];
        qbt1 *= exp_E_MLstem(type, (i>1) ? S1[i-1] : -1, (j<n) ? S1[j+1] : -1, pf_params);
        qqm[i] += qbt1;
      }

      if (qm1) qm1[jindx[j]+i] = qqm[i]; /* for stochastic backtracking */


      /*construction of qm matrix containing multiple loop
        partition function contributions from segment i,j */
      temp = 0.0;
      ii = my_iindx[i];  /* ii-k=[i,k] */

      for (k=i+1; k<=j; k++) {
        if (SAME_STRAND(k-1,k)) temp += (qm[ii-(k-1)])*qqm[k];
        if (SAME_STRAND(i,k))   temp += expMLbase[k-i]*qqm[k];

      }

      qm[ij] = (temp + qqm[i]);

      /*auxiliary matrix qq for cubic order q calculation below */
      qbt1 = qb[ij];
      if (type) {
        qbt1 *= exp_E_ExtLoop(type, ((i>1)&&(SAME_STRAND(i-1,i))) ? S1[i-1] : -1, ((j<n)&&(SAME_STRAND(j,j+1))) ? S1[j+1] : -1, pf_params);
      }
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
        snprintf(msg, 127, "overflow in co_pf_fold while calculating q[%d,%d]\n"
                "use larger pf_scale", i,j);
        nrerror(msg);
      }
    }
    tmp = qq1;  qq1 =qq;  qq =tmp;
    tmp = qqm1; qqm1=qqm; qqm=tmp;
  }
}

/* backward recursion of pf cofolding */
PRIVATE void pf_co_bppm(const char *sequence, char *structure){
  int         n, i,j,k,l, ij, kl, ii, ll, type, type_2, tt, ov=0;
  FLT_OR_DBL  temp, Qmax=0, prm_MLb;
  FLT_OR_DBL  prmt,prmt1;
  FLT_OR_DBL  *tmp;
  FLT_OR_DBL  expMLclosing;
  double      max_real;

  max_real = (sizeof(FLT_OR_DBL) == sizeof(float)) ? FLT_MAX : DBL_MAX;
  n = (int) strlen(sequence);

  expMLclosing = pf_params->expMLclosing;

  /* backtracking to construct binding probabilities of pairs*/
  if ((S != NULL) && (S1 != NULL)) {
    FLT_OR_DBL   *Qlout, *Qrout;
    Qmax=0;
    Qrout=(FLT_OR_DBL *)space(sizeof(FLT_OR_DBL) * (n+2));
    Qlout=(FLT_OR_DBL *)space(sizeof(FLT_OR_DBL) * (cut_point+2));

    for (k=1; k<=n; k++) {
      q1k[k] = q[my_iindx[1] - k];
      qln[k] = q[my_iindx[k] -n];
    }
    q1k[0] = 1.0;
    qln[n+1] = 1.0;

    /*    pr = q;     /  * recycling */

    /* 1. exterior pair i,j and initialization of pr array */
    for (i=1; i<=n; i++) {
      for (j=i; j<=MIN2(i+TURN,n); j++) probs[my_iindx[i]-j] = 0;
      for (j=i+TURN+1; j<=n; j++) {
        ij = my_iindx[i]-j;
        type = ptype[ij];
        if (type&&(qb[ij]>0.)) {
          probs[ij] = q1k[i-1]*qln[j+1]/q1k[n];
          probs[ij] *= exp_E_ExtLoop(type, ((i>1)&&(SAME_STRAND(i-1,i))) ? S1[i-1] : -1, ((j<n)&&(SAME_STRAND(j,j+1))) ? S1[j+1] : -1, pf_params);
        } else
          probs[ij] = 0;
      }
    }

    for (l=n; l>TURN+1; l--) {

      /* 2. bonding k,l as substem of 2:loop enclosed by i,j */
      for (k=1; k<l-TURN; k++) {
        kl = my_iindx[k]-l;
        type_2 = ptype[kl]; type_2 = rtype[type_2];
        if (qb[kl]==0) continue;

        for (i=MAX2(1,k-MAXLOOP-1); i<=k-1; i++)
          for (j=l+1; j<=MIN2(l+ MAXLOOP -k+i+2,n); j++) {
            if ((SAME_STRAND(i,k))&&(SAME_STRAND(l,j))){
              ij = my_iindx[i] - j;
              type = ptype[ij];
              if ((probs[ij]>0)) {
                probs[kl] += probs[ij]*exp_E_IntLoop(k-i-1, j-l-1, type, type_2,
                                               S1[i+1], S1[j-1], S1[k-1], S1[l+1], pf_params)*scale[k-i+j-l];
              }
            }
          }
      }
      /* 3. bonding k,l as substem of multi-loop enclosed by i,j */
      prm_MLb = 0.;
      if ((l<n)&&(SAME_STRAND(l,l+1)))
        for (k=2; k<l-TURN; k++) {
          i = k-1;
          prmt = prmt1 = 0.0;

          ii = my_iindx[i];     /* ii-j=[i,j]     */
          ll = my_iindx[l+1];   /* ll-j=[l+1,j] */
          tt = ptype[ii-(l+1)]; tt=rtype[tt];
          if (SAME_STRAND(i,k)){
            prmt1 = probs[ii-(l+1)]*expMLclosing;
            prmt1 *= exp_E_MLstem(tt, S1[l], S1[i+1], pf_params);
            for (j=l+2; j<=n; j++) {
              if (SAME_STRAND(j-1,j)){ /*??*/
                tt = ptype[ii-j]; tt = rtype[tt];
                prmt += probs[ii-j]*exp_E_MLstem(tt, S1[j-1], S1[i+1], pf_params)*qm[ll-(j-1)];
              }
            }
          }
          kl = my_iindx[k]-l;
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

          for (i=1;i<=k-2; i++) {
            if ((SAME_STRAND(i,i+1))&&(SAME_STRAND(k-1,k))){
              temp += prml[i]*qm[my_iindx[i+1] - (k-1)];
            }
          }
          temp *= exp_E_MLstem( tt,
                                ((k>1)&&SAME_STRAND(k-1,k)) ? S1[k-1] : -1,
                                ((l<n)&&SAME_STRAND(l,l+1)) ? S1[l+1] : -1,
                                pf_params) * scale[2];
          probs[kl] += temp;

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

        } /* end for (k=..) multloop*/
      else  /* set prm_l to 0 to get prm_l1 to be 0 */
        for (i=0; i<=n; i++) prm_l[i]=0;

      tmp = prm_l1; prm_l1=prm_l; prm_l=tmp;
      /*computation of .(..(...)..&..). type features?*/
      if (cut_point<=0) continue;  /* no .(..(...)..&..). type features*/
      if ((l==n)||(l<=2)) continue; /* no .(..(...)..&..). type features*/
      /*new version with O(n^3)??*/
      if (l>cut_point) {
        if (l<n) {
          int t,kt;
          for (t=n; t>l; t--) {
            for (k=1; k<cut_point; k++) {
              kt=my_iindx[k]-t;
              type=rtype[ptype[kt]];
              temp = probs[kt] * exp_E_ExtLoop(type, S1[t-1], (SAME_STRAND(k,k+1)) ? S1[k+1] : -1, pf_params) * scale[2];
              if (l+1<t)               temp*=q[my_iindx[l+1]-(t-1)];
              if (SAME_STRAND(k,k+1))  temp*=q[my_iindx[k+1]-(cut_point-1)];
              Qrout[l]+=temp;
            }
          }
        }
        for (k=l-1; k>=cut_point; k--) {
          if (qb[my_iindx[k]-l]) {
            kl=my_iindx[k]-l;
            type=ptype[kl];
            temp = Qrout[l];
            temp *= exp_E_ExtLoop(type, (k>cut_point) ? S1[k-1] : -1, (l < n) ? S1[l+1] : -1, pf_params);
            if (k>cut_point) temp*=q[my_iindx[cut_point]-(k-1)];
            probs[kl]+=temp;
          }
        }
      }
      else if (l==cut_point ) {
        int t, sk,s;
        for (t=2; t<cut_point;t++) {
          for (s=1; s<t; s++) {
            for (k=cut_point; k<=n; k++) {
              sk=my_iindx[s]-k;
              if (qb[sk]) {
                type=rtype[ptype[sk]];
                temp=probs[sk]*exp_E_ExtLoop(type, (SAME_STRAND(k-1,k)) ? S1[k-1] : -1, S1[s+1], pf_params)*scale[2];
                if (s+1<t)               temp*=q[my_iindx[s+1]-(t-1)];
                if (SAME_STRAND(k-1,k))  temp*=q[my_iindx[cut_point]-(k-1)];
                Qlout[t]+=temp;
              }
            }
          }
        }
      }
      else if (l<cut_point) {
        for (k=1; k<l; k++) {
          if (qb[my_iindx[k]-l]) {
            type=ptype[my_iindx[k]-l];
            temp=Qlout[k];
            temp *= exp_E_ExtLoop(type, (k>1) ? S1[k-1] : -1, (l<(cut_point-1)) ? S1[l+1] : -1, pf_params);
            if (l+1<cut_point) temp*=q[my_iindx[l+1]-(cut_point-1)];
            probs[my_iindx[k]-l]+=temp;
          }
        }
      }
    }  /* end for (l=..)   */
    free(Qlout);
    free(Qrout);
    for (i=1; i<=n; i++)
      for (j=i+TURN+1; j<=n; j++) {
        ij = my_iindx[i]-j;
        probs[ij] *= qb[ij];
      }

    if (structure!=NULL)
      bppm_to_structure(structure, probs, n);
  }   /* end if (do_backtrack)*/

  if (ov>0) fprintf(stderr, "%d overflows occurred while backtracking;\n"
                    "you might try a smaller pf_scale than %g\n",
                    ov, pf_params->pf_scale);
}


PRIVATE void scale_pf_params(unsigned int length, pf_paramT *parameters){
  unsigned int  i;
  double        kT, scaling_factor;

  if(pf_params) free(pf_params);

  if(parameters){
    pf_params = get_boltzmann_factor_copy(parameters);
  } else {
    model_detailsT md;
    set_model_details(&md);
    pf_params = get_boltzmann_factors(temperature, alpha, md, pf_scale);
  }

  scaling_factor  = pf_params->pf_scale;
  kT              = pf_params->kT;        /* kT in cal/mol  */

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

/*----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/

PUBLIC void update_co_pf_params(int length){
  update_co_pf_params_par(length, NULL);
}

PUBLIC void update_co_pf_params_par(int length, pf_paramT *parameters){
  make_pair_matrix();
  scale_pf_params((unsigned) length, parameters);
}

/*---------------------------------------------------------------------------*/
PRIVATE void make_ptypes(const short *S, const char *structure) {
  int n,i,j,k,l;
  int noLP = pf_params->model_details.noLP;

  n=S[0];
  for (k=1; k<=n-TURN-1; k++)
    for (l=1; l<=2; l++) {
      int type,ntype=0,otype=0;
      i=k; j = i+TURN+l;
      if (j>n) continue;
      type = pair[S[i]][S[j]];
      while ((i>=1)&&(j<=n)) {
        if ((i>1)&&(j<n)) ntype = pair[S[i-1]][S[j+1]];
        if (noLP && (!otype) && (!ntype))
          type = 0; /* i.j can only form isolated pairs */
        qb[my_iindx[i]-j] = 0.;
        ptype[my_iindx[i]-j] = (char) type;
        otype =  type;
        type  = ntype;
        i--; j++;
      }

    }

  if (struct_constrained&&(structure!=NULL)) {
    constrain_ptypes(structure, (unsigned int)n, ptype, NULL, TURN, 1);
    for(j=1; j<=n; j++) {
      switch (structure[j-1]) {
        case 'l': /*only intramolecular basepairing*/
                  if (j<cut_point) for (l=cut_point; l<=n; l++) ptype[my_iindx[j]-l] = 0;
                  else for (l=1; l<cut_point; l++) ptype[my_iindx[l]-j] =0;
                  break;
        case 'e': /*only intermolecular bp*/
                  if (j<cut_point) {
                    for (l=1; l<j; l++) ptype[my_iindx[l]-j] =0;
                    for (l=j+1; l<cut_point; l++) ptype[my_iindx[j]-l] = 0;
                  }
                  else {
                    for (l=cut_point; l<j; l++) ptype[my_iindx[l]-j] =0;
                    for (l=j+1; l<=n; l++) ptype[my_iindx[j]-l] = 0;
                  }
                  break;
      }
    }
  }
  if (mirnatog==1) {   /*microRNA toggle: no intramolec. bp in 2. molec*/
    for (j=cut_point; j<n; j++) {
      for (l=j+1; l<=n; l++) {
        ptype[my_iindx[j]-l] = 0;
      }
    }
  }
}

/*
  stochastic backtracking in pf_fold arrays
  returns random structure S with Boltzman probabilty
  p(S) = exp(-E(S)/kT)/Z
*/
PRIVATE void backtrack_qm1(int i,int j) {
  /* i is paired to l, i<l<j; backtrack in qm1 to find l */
  int ii, l, type;
  double qt, r;
  r = urn() * qm1[jindx[j]+i];
  ii = my_iindx[i];
  for (qt=0., l=i+TURN+1; l<=j; l++) {
    type = ptype[ii-l];
    if (type)
      qt +=  qb[ii-l]*exp_E_MLstem(type, S1[i-1], S1[l+1], pf_params) * expMLbase[j-l];
    if (qt>=r) break;
  }
  if (l>j) nrerror("backtrack failed in qm1");
  backtrack(i,l);
}

PRIVATE void backtrack(int i, int j) {
  int noGUclosure = pf_params->model_details.noGUclosure;

  do {
    double r, qbt1;
    int k, l, type, u, u1;

    pstruc[i-1] = '('; pstruc[j-1] = ')';

    r = urn() * qb[my_iindx[i]-j];
    type = ptype[my_iindx[i]-j];
    u = j-i-1;
    /*hairpin contribution*/
    if (((type==3)||(type==4))&&noGUclosure) qbt1 = 0;
    else
      qbt1 = exp_E_Hairpin(u, type, S1[i+1], S1[j-1], sequence+i-1, pf_params)*scale[u+2];

    if (qbt1>r) return; /* found the hairpin we're done */

    for (k=i+1; k<=MIN2(i+MAXLOOP+1,j-TURN-2); k++) {
      u1 = k-i-1;
      for (l=MAX2(k+TURN+1,j-1-MAXLOOP+u1); l<j; l++) {
        int type_2;
        type_2 = ptype[my_iindx[k]-l];
        if (type_2) {
          type_2 = rtype[type_2];
          qbt1 += qb[my_iindx[k]-l] *
            exp_E_IntLoop(u1, j-l-1, type, type_2,
                          S1[i+1], S1[j-1], S1[k-1], S1[l+1], pf_params)*scale[u1+j-l+1];
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
    ii = my_iindx[i]; /* ii-j=[i,j] */
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
    while (j>i) {
      /* now backtrack  [i ... j] in qm[] */
      jj = jindx[j];
      ii = my_iindx[i];
      r = urn() * qm[ii - j];
      qt = qm1[jj+i]; k=i;
      if (qt<r)
        for (k=i+1; k<=j; k++) {
          qt += (qm[ii-(k-1)]+expMLbase[k-i])*qm1[jj+k];
          if (qt >= r) break;
        }
      if (k>j) nrerror("backtrack failed in qm");

      backtrack_qm1(k,j);

      if (k<i+TURN) break; /* no more pairs */
      r = urn() * (qm[ii-(k-1)] + expMLbase[k-i]);
      if (expMLbase[k-i] >= r) break; /* no more pairs */
      j = k-1;
    }
  }
}

PUBLIC void compute_probabilities(double FAB, double FA,double FB,
                                  struct plist *prAB,
                                  struct plist *prA, struct plist *prB,
                                  int Alength) {
  /*computes binding probabilities and dimer free energies*/
  int i, j;
  double pAB;
  double mykT;
  struct plist  *lp1, *lp2;
  int offset;

  mykT=pf_params->kT/1000.;

  /* pair probabilities in pr are relative to the null model (without DuplexInit) */

  /*Compute probabilities pAB, pAA, pBB*/

  pAB=1.-exp((1/mykT)*(FAB-FA-FB));

  /* compute pair probabilities given that it is a dimer */
  /* AB dimer */
  offset=0;
  lp2=prA;
  if (pAB>0)
    for (lp1=prAB; lp1->j>0; lp1++) {
      float pp=0;
      i=lp1->i; j=lp1->j;
      while (offset+lp2->i < i && lp2->i>0) lp2++;
      if (offset+lp2->i == i)
        while ((offset+lp2->j) < j  && (lp2->j>0)) lp2++;
      if (lp2->j == 0) {lp2=prB; offset=Alength;}/* jump to next list */
      if ((offset+lp2->i==i) && (offset+lp2->j ==j)) {
        pp = lp2->p;
        lp2++;
      }
      lp1->p=(lp1->p-(1-pAB)*pp)/pAB;
      if(lp1->p < 0.){
        warn_user("part_func_co: numeric instability detected, probability below zero!");
        lp1->p = 0.;
      }
    }
  return;
}

PRIVATE double *Newton_Conc(double KAB, double KAA, double KBB, double concA, double concB,double* ConcVec) {
  double TOL, EPS, xn, yn, det, cA, cB;
  int i=0;
  /*Newton iteration for computing concentrations*/
  cA=concA;
  cB=concB;
  TOL=1e-6; /*Tolerance for convergence*/
  ConcVec=(double*)space(5*sizeof(double)); /* holds concentrations */
  do {
    det = (4.0 * KAA * cA + KAB *cB + 1.0) * (4.0 * KBB * cB + KAB *cA + 1.0) - (KAB *cB) * (KAB *cA);
    xn  = ( (2.0 * KBB * cB*cB + KAB *cA *cB + cB - concB) * (KAB *cA) -
            (2.0 * KAA * cA*cA + KAB *cA *cB + cA - concA) * (4.0 * KBB * cB + KAB *cA + 1.0) ) /det;
    yn  = ( (2.0 * KAA * cA*cA + KAB *cA *cB + cA - concA) * (KAB *cB) -
            (2.0 * KBB * cB*cB + KAB *cA *cB + cB - concB) * (4.0 * KAA * cA + KAB *cB + 1.0) ) /det;
    EPS = fabs(xn/cA) + fabs(yn/cB);
    cA += xn;
    cB += yn;
    i++;
    if (i>10000) {
      fprintf(stderr, "Newton did not converge after %d steps!!\n",i);
      break;
    }
  } while(EPS>TOL);

  ConcVec[0]= cA*cB*KAB ;/*AB concentration*/
  ConcVec[1]= cA*cA*KAA ;/*AA concentration*/
  ConcVec[2]= cB*cB*KBB ;/*BB concentration*/
  ConcVec[3]= cA;        /* A concentration*/
  ConcVec[4]= cB;        /* B concentration*/

  return ConcVec;
}

PUBLIC struct ConcEnt *get_concentrations(double FcAB, double FcAA, double FcBB, double FEA, double FEB, double *startconc)
{
  /*takes an array of start concentrations, computes equilibrium concentrations of dimers, monomers, returns array of concentrations in strucutre ConcEnt*/
  double *ConcVec;
  int i;
  struct ConcEnt *Concentration;
  double KAA, KAB, KBB, kT;

  kT=pf_params->kT/1000.;
  Concentration=(struct ConcEnt *)space(20*sizeof(struct ConcEnt));
 /* Compute equilibrium constants */
  /* again note the input free energies are not from the null model (without DuplexInit) */

  KAA = exp(( 2.0 * FEA - FcAA)/kT);
  KBB = exp(( 2.0 * FEB - FcBB)/kT);
  KAB = exp(( FEA + FEB - FcAB)/kT);
  /* printf("Kaa..%g %g %g\n", KAA, KBB, KAB); */
  for (i=0; ((startconc[i]!=0)||(startconc[i+1]!=0));i+=2) {
    ConcVec=Newton_Conc(KAB, KAA, KBB, startconc[i], startconc[i+1], ConcVec);
    Concentration[i/2].A0=startconc[i];
    Concentration[i/2].B0=startconc[i+1];
    Concentration[i/2].ABc=ConcVec[0];
    Concentration[i/2].AAc=ConcVec[1];
    Concentration[i/2].BBc=ConcVec[2];
    Concentration[i/2].Ac=ConcVec[3];
    Concentration[i/2].Bc=ConcVec[4];

   if (!(((i+2)/2)%20))  {
     Concentration=(struct ConcEnt *)xrealloc(Concentration,((i+2)/2+20)*sizeof(struct ConcEnt));
     }
    free(ConcVec);
  }

  return Concentration;
}

PUBLIC FLT_OR_DBL *export_co_bppm(void){
  return probs;
}

/*###########################################*/
/*# deprecated functions below              #*/
/*###########################################*/


PUBLIC struct plist *get_plist(struct plist *pl, int length, double cut_off) {
  int i, j,n, count;
  /*get pair probibilities out of pr array*/
  count=0;
  n=2;
  for (i=1; i<length; i++) {
    for (j=i+1; j<=length; j++) {
      if (pr[my_iindx[i]-j]<cut_off) continue;
      if (count==n*length-1) {
        n*=2;
        pl=(struct plist *)xrealloc(pl,n*length*sizeof(struct plist));
      }
      pl[count].i=i;
      pl[count].j=j;
      pl[count++].p=pr[my_iindx[i]-j];
      /*      printf("gpl: %2d %2d %.9f\n",i,j,pr[my_iindx[i]-j]);*/
    }
  }
  pl[count].i=0;
  pl[count].j=0; /*->??*/
  pl[count++].p=0.;
  pl=(struct plist *)xrealloc(pl,(count)*sizeof(struct plist));
  return pl;
}

PUBLIC void init_co_pf_fold(int length){ /* DO NOTHING */ }
