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
#include "utils.h"
#include "energy_par.h"
#include "fold_vars.h"
#include "pair_mat.h"
#include "PS_dot.h"
#include "part_func_co.h"
/*@unused@*/
static char rcsid[] UNUSED = "$Id: part_func_co.c,v 1.10 2007/05/10 17:27:01 ivo Exp $";

#define MAX(x,y) (((x)>(y)) ? (x) : (y))
#define MIN(x,y) (((x)<(y)) ? (x) : (y))
#define PUBLIC
#define PRIVATE static

PUBLIC  cofoldF co_pf_fold(char *sequence, char *structure);
PUBLIC  void  init_co_pf_fold(int length);
PUBLIC  void  free_co_pf_arrays(void);
PUBLIC  void  update_co_pf_params(int length);
PUBLIC  char  co_bppm_symbol(float *x);

PUBLIC struct ConcEnt *get_concentrations(double FEAB, double FEAA, double FEBB, double FEA, double FEB, double *startconc);
PUBLIC struct plist *get_plist(struct plist *pl, int length, double cut_off);
PUBLIC struct plist *get_mfe_plist(struct plist *pl);
PUBLIC double F_monomer[2]; /* free energies of the two monomers */
/* PUBLIC int make_probsum(int length, char *name); */
PRIVATE double *Newton_Conc(double ZAB, double ZAA, double ZBB, double concA, double concB,double* ConcVec);
PRIVATE void  sprintf_bppm(int length, char *structure);
PRIVATE void  scale_pf_params(unsigned int length);
PRIVATE void  get_arrays(unsigned int length);
PRIVATE double expLoopEnergy(int u1, int u2, int type, int type2,
			     short si1, short sj1, short sp1, short sq1);
PRIVATE double expHairpinEnergy(int u, int type, short si1, short sj1,
				const char *string);
PRIVATE void make_ptypes(const short *S, const char *structure);

PRIVATE FLT_OR_DBL expMLclosing, expMLintern[NBPAIRS+1], *expMLbase;
PRIVATE FLT_OR_DBL expTermAU;
PRIVATE FLT_OR_DBL expdangle5[NBPAIRS+1][5], expdangle3[NBPAIRS+1][5];
PRIVATE FLT_OR_DBL lxc, exptetra[40], expTriloop[40];
PRIVATE FLT_OR_DBL expstack[NBPAIRS+1][NBPAIRS+1];
PRIVATE FLT_OR_DBL expmismatchI[NBPAIRS+1][5][5],
  expmismatchH[NBPAIRS+1][5][5], expmismatchM[NBPAIRS+1][5][5];
PRIVATE FLT_OR_DBL expint11[NBPAIRS+1][NBPAIRS+1][5][5];
PRIVATE FLT_OR_DBL expint21[NBPAIRS+1][NBPAIRS+1][5][5][5];
PRIVATE FLT_OR_DBL expint22[NBPAIRS+1][NBPAIRS+1][5][5][5][5];
PRIVATE FLT_OR_DBL *exphairpin;
PRIVATE FLT_OR_DBL expbulge[MAXLOOP+1];
PRIVATE FLT_OR_DBL expInit;
PRIVATE FLT_OR_DBL expinternal[MAXLOOP+1];
PRIVATE FLT_OR_DBL expninio[5][MAXLOOP+1];
PRIVATE FLT_OR_DBL *q, *qb, *qm, *qm1, *qqm, *qqm1, *qq, *qq1;
PRIVATE FLT_OR_DBL *prml, *prm_l, *prm_l1, *q1k, *qln;
PRIVATE FLT_OR_DBL *scale;
PRIVATE char *ptype; /* precomputed array of pair types */
PRIVATE int *jindx;
PRIVATE int init_length; /* length in last call to init_pf_fold() */
#define ISOLATED  256.0
#undef TURN
#define TURN 0
#define SAME_STRAND(I,J) (((I)>=cut_point)||((J)<cut_point))
/* #define SAME_STRAND(I,J) (((J)<cut_point)||((I)>=cut_point2)||(((I)>=cut_point)&&((J)<cut_point2)))
 */
int mirnatog = 0;

/*-----------------------------------------------------------------*/
static  short *S, *S1;
PUBLIC cofoldF co_pf_fold(char *sequence, char *structure)
{

  int n, i,j,k,l, ij, kl, u,u1,ii,ll, type, type_2, tt, ov=0;
  FLT_OR_DBL temp, Q, Qmax=0, prm_MLb;
  FLT_OR_DBL prmt,prmt1;
  FLT_OR_DBL qbt1, *tmp;
  cofoldF X;
  double free_energy, max_real;

  max_real = (sizeof(FLT_OR_DBL) == sizeof(float)) ? FLT_MAX : DBL_MAX;
  n = (int) strlen(sequence);
  if (n>init_length) init_co_pf_fold(n);  /* (re)allocate space */
 /* printf("mirnatog=%d\n",mirnatog); */
  S = (short *) xrealloc(S, sizeof(short)*(n+1));
  S1= (short *) xrealloc(S1, sizeof(short)*(n+1));
  S[0] = n;
  for (l=1; l<=n; l++) {
    S[l]  = (short) encode_char(toupper(sequence[l-1]));
    S1[l] = alias[S[l]];
  }
  make_ptypes(S, structure);

  /*array initialization ; qb,qm,q
    qb,qm,q (i,j) are stored as ((n+1-i)*(n-i) div 2 + n+1-j */

  /* for (d=0; d<=TURN; d++) */
  for (i=1; i<=n/*-d*/; i++) {
      ij = iindx[i]-i;
      q[ij]=scale[1];
      qb[ij]=qm[ij]=0.0;
    }

  for (i=1; i<=n; i++)
    qq[i]=qq1[i]=qqm[i]=qqm1[i]=prm_l[i]=prm_l1[i]=prml[i]=0;

  for (j=TURN+2;j<=n; j++) {
    for (i=j-TURN-1; i>=1; i--) {
      /* construction of partition function of segment i,j*/
       /*firstly that given i bound to j : qb(i,j) */
      u = j-i-1; ij = iindx[i]-j;
      type = ptype[ij];
      qbt1=0;
      if (type!=0) {
	/*hairpin contribution*/
	if SAME_STRAND(i,j){
	  if (((type==3)||(type==4))&&no_closingGU) qbt1 = 0;
	  else
	    qbt1 = expHairpinEnergy(u, type, S1[i+1], S1[j-1], sequence+i-1);

	}

	/* interior loops with interior pair k,l */
	for (k=i+1; k<=MIN(i+MAXLOOP+1,j-TURN-2); k++) {
	  u1 = k-i-1;
	  for (l=MAX(k+TURN+1,j-1-MAXLOOP+u1); l<j; l++) {
	    if ((SAME_STRAND(i,k))&&(SAME_STRAND(l,j))){
	      type_2 = ptype[iindx[k]-l];
	      if (type_2) {
		type_2 = rtype[type_2];
		qbt1 += qb[iindx[k]-l] *
		  expLoopEnergy(u1, j-l-1, type, type_2,
				S1[i+1], S1[j-1], S1[k-1], S1[l+1]);
	      }
	    }
	  }
	}
	/*multiple stem loop contribution*/
	ii = iindx[i+1]; /* ii-k=[i+1,k-1] */
	temp = 0.0;
	if (SAME_STRAND(i,i+1) && SAME_STRAND(j-1,j)) {
	  for (k=i+2; k<=j-1; k++) {
	    if (SAME_STRAND(k-1,k))
	      temp += qm[ii-(k-1)]*qqm1[k];
	  }
	  tt = rtype[type];
	  temp*=expMLclosing*expMLintern[tt]*scale[2];
	  temp*=expdangle3[tt][S1[i+1]];
	  temp*=expdangle5[tt][S1[j-1]];
	  qbt1 += temp;
	}
	/*qc contribution*/
	temp=0.0;
	if (!SAME_STRAND(i,j)){
	  tt = rtype[type];
	  temp=q[iindx[i+1]-(cut_point-1)]*q[iindx[cut_point]-(j-1)];
	  if ((j==cut_point)&&(i==cut_point-1)) temp=scale[2];
	  else if (i==cut_point-1) temp=q[iindx[cut_point]-(j-1)]*scale[1];
	  else if (j==cut_point) temp=q[iindx[i+1]-(cut_point-1)]*scale[1];
	  if (j>cut_point) temp*=scale[1];
	  if (i<cut_point-1) temp*=scale[1];
	  if (SAME_STRAND(i,i+1)) temp*=expdangle3[tt][S1[i+1]];
	  else if (type>2) temp*= expTermAU;
	  if (SAME_STRAND(j-1,j)) temp*=expdangle5[tt][S1[j-1]];

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
	qbt1 = qb[ij]*expMLintern[type];
	if (i>1) qbt1 *= expdangle5[type][S1[i-1]];
	if (j<n) qbt1 *= expdangle3[type][S1[j+1]];
	else if (type>2) qbt1 *= expTermAU;
	qqm[i] += qbt1;
      }

      if (qm1) qm1[jindx[j]+i] = qqm[i]; /* for stochastic backtracking */


      /*construction of qm matrix containing multiple loop
	partition function contributions from segment i,j */
      temp = 0.0;
      ii = iindx[i];  /* ii-k=[i,k] */

      for (k=i+1; k<=j; k++) {
	if (SAME_STRAND(k-1,k))   temp += (qm[ii-(k-1)])*qqm[k];
	if (SAME_STRAND(i,k)) temp +=expMLbase[k-i]*qqm[k];

      }

      qm[ij] = (temp + qqm[i]);

      /*auxiliary matrix qq for cubic order q calculation below */
      qbt1 = qb[ij];
      if (type) {
	if ((i>1)&&(SAME_STRAND(i-1,i))) qbt1 *= expdangle5[type][S1[i-1]];
	if ((j<n)&&(SAME_STRAND(j,j+1))) qbt1 *= expdangle3[type][S1[j+1]];
	else if (type>2) qbt1 *= expTermAU;
      }
      qq[i] = qq1[i]*scale[1] + qbt1;
       /*construction of partition function for segment i,j */
      temp = 1.0*scale[1+j-i] + qq[i];
      for (k=i; k<=j-1; k++) {
	temp += q[ii-k]*qq[k+1];
      }
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
  if (backtrack_type=='C')      Q = qb[iindx[1]-n];
  else if (backtrack_type=='M') Q = qm[iindx[1]-n];
  else Q = q[iindx[1]-n];
  /* ensemble free energy in Kcal/mol */
  if (Q<=FLT_MIN) fprintf(stderr, "pf_scale too large\n");
  free_energy = (-log(Q)-n*log(pf_scale))*(temperature+K0)*GASCONST/1000.0;
  /* in case we abort because of floating point errors */
  if (n>1600) fprintf(stderr, "free energy = %8.2f\n", free_energy);
  /*probability of molecules being bound together*/


  /*Computation of "real" Partition function*/
  /*Need that for concentrations*/
  if (cut_point>0){
    double kT, pbound, QAB, QToT, Qzero;

    kT = (temperature+K0)*GASCONST/1000.0;
    Qzero=q[iindx[1]-n];
    QAB=(q[iindx[1]-n]-q[iindx[1]-(cut_point-1)]*q[iindx[cut_point]-n])*expInit;
    /*correction for symmetry*/
    if((n-(cut_point-1)*2)==0) {
      if ((strncmp(sequence, sequence+cut_point-1, cut_point-1))==0) {
	QAB/=2;
      }}

    QToT=q[iindx[1]-(cut_point-1)]*q[iindx[cut_point]-n]+QAB;
    pbound=1-(q[iindx[1]-(cut_point-1)]*q[iindx[cut_point]-n]/QToT);
     X.FAB  = -kT*(log(QToT)+n*log(pf_scale));
    X.F0AB = -kT*(log(Qzero)+n*log(pf_scale));
    X.FcAB = (QAB>1e-17) ? -kT*(log(QAB)+n*log(pf_scale)) : 999;
    X.FA = -kT*(log(q[iindx[1]-(cut_point-1)]) + (cut_point-1)*log(pf_scale));
    X.FB = -kT*(log(q[iindx[cut_point]-n]) + (n-cut_point+1)*log(pf_scale));

    /* printf("QAB=%.9f\tQtot=%.9f\n",QAB/scale[n],QToT/scale[n]);*/
  }
  else {
    X.FA = X.FB = X.FAB = X.F0AB = free_energy;
    X.FcAB = 0;
  }
  /* printf("freen=%.6f\n",free_energy); whereto?*/
  /* backtracking to construct binding probabilities of pairs*/
  /* printf("qis %f\n",q[iindx[1]-(cut_point-1)]);*/
  /*new: expInit added*/ /*new*/
  if (do_backtrack) {
    FLT_OR_DBL   *Qlout, *Qrout;
    Qmax=0;
    Qrout=(FLT_OR_DBL *)space(sizeof(FLT_OR_DBL) * (n+2));
    Qlout=(FLT_OR_DBL *)space(sizeof(FLT_OR_DBL) * (cut_point+2));

    for (k=1; k<=n; k++) {
      q1k[k] = q[iindx[1] - k];
      qln[k] = q[iindx[k] -n];
    }
    q1k[0] = 1.0;
    qln[n+1] = 1.0;

    /*    pr = q;     /  * recycling */

    /* 1. exterior pair i,j and initialization of pr array */
    for (i=1; i<=n; i++) {
      for (j=i; j<=MIN(i+TURN,n); j++) pr[iindx[i]-j] = 0;
      for (j=i+TURN+1; j<=n; j++) {
	ij = iindx[i]-j;
	type = ptype[ij];
	if (type&&(qb[ij]>0.)) {
	  pr[ij] = q1k[i-1]*qln[j+1]/q1k[n];
	  if ((i>1)&&(SAME_STRAND(i-1,i))) pr[ij] *= expdangle5[type][S1[i-1]];
	  if ((j<n)&&(SAME_STRAND(j,j+1))) pr[ij] *= expdangle3[type][S1[j+1]];
	  else if (type>2) pr[ij] *= expTermAU;
	} else
	  pr[ij] = 0;
      }
    }

    for (l=n; l>TURN+1; l--) {

      /* 2. bonding k,l as substem of 2:loop enclosed by i,j */
      for (k=1; k<l-TURN; k++) {
	kl = iindx[k]-l;
	type_2 = ptype[kl]; type_2 = rtype[type_2];
	if (qb[kl]==0) continue;

	for (i=MAX(1,k-MAXLOOP-1); i<=k-1; i++)
	  for (j=l+1; j<=MIN(l+ MAXLOOP -k+i+2,n); j++) {
	    if ((SAME_STRAND(i,k))&&(SAME_STRAND(l,j))){
	      ij = iindx[i] - j;
	      type = ptype[ij];
	      if ((pr[ij]>0)) {
		pr[kl] += pr[ij]*expLoopEnergy(k-i-1, j-l-1, type, type_2,
					       S1[i+1], S1[j-1], S1[k-1], S1[l+1]);
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

	  ii = iindx[i];     /* ii-j=[i,j]     */
	  ll = iindx[l+1];   /* ll-j=[l+1,j] */
	  tt = ptype[ii-(l+1)]; tt=rtype[tt];
	  if (SAME_STRAND(i,k)){
	    prmt1 = pr[ii-(l+1)]*expMLclosing*expMLintern[tt]*
	      expdangle3[tt][S1[i+1]]*expdangle5[tt][S1[l]];
	    for (j=l+2; j<=n; j++) {
	      if (SAME_STRAND(j-1,j)){ /*??*/
		tt = ptype[ii-j]; tt = rtype[tt];
		prmt += pr[ii-j]*expdangle3[tt][S1[i+1]]*
		  expdangle5[tt][S1[j-1]] *qm[ll-(j-1)];
	      }
	    }
	  }
	  kl = iindx[k]-l;
	  tt = ptype[kl];
	  prmt *= expMLclosing*expMLintern[tt];
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
	      temp += prml[i]*qm[iindx[i+1] - (k-1)];
	    }
	  }
	  temp *= expMLintern[tt]*scale[2];
	  if ((k>1)&&SAME_STRAND(k-1,k)) temp *= expdangle5[tt][S1[k-1]];
	  if ((l<n)&&SAME_STRAND(l,l+1)) temp *= expdangle3[tt][S1[l+1]];
	  pr[kl] += temp;

	  if (pr[kl]>Qmax) {
	    Qmax = pr[kl];
	    if (Qmax>max_real/10.)
	      fprintf(stderr, "P close to overflow: %d %d %g %g\n",
		      i, j, pr[kl], qb[kl]);
	  }
	  if (pr[kl]>=max_real) {
	    ov++;
	    pr[kl]=FLT_MAX;
	  }

	} /* end for (k=..) multloop*/
      else  /* set prm_l to 0 to get prm_l1 to be 0 */
	for (i=0; i<=n; i++) prm_l[i]=0;

     tmp = prm_l1; prm_l1=prm_l; prm_l=tmp;
      /*computation of .(..(...)..&..). type features?*/
     if (cut_point<=0) continue;  /* no .(..(...)..&..). type features*/
     if ((l==n)||(l<=2)) continue; /* no .(..(...)..&..). type features*/
     /*new version with O(n^3)??*/
#if 1
     if (l>cut_point) {
       if (l<n) {
	 int t,kt;
	 for (t=n; t>l; t--) {
	   for (k=1; k<cut_point; k++) {
	     kt=iindx[k]-t;
	     type=rtype[ptype[kt]];
	     temp=pr[kt]*scale[2]*expdangle5[type][S1[t-1]];
	     if (l+1<t) temp*=q[iindx[l+1]-(t-1)];
	     if (SAME_STRAND(k,k+1)) {
	       temp*=q[iindx[k+1]-(cut_point-1)];
	       temp*=expdangle3[type][S1[k+1]];
	     }
	     else if (type>2) temp*=expTermAU;
	     Qrout[l]+=temp;
	   }
	 }
       }
       for (k=l-1; k>=cut_point; k--) {
	 if (qb[iindx[k]-l]) {
	   kl=iindx[k]-l;
	   type=ptype[kl];
	   temp = Qrout[l];
	   if (k>cut_point) {
	     temp*=q[iindx[cut_point]-(k-1)]*expdangle5[type][S1[k-1]];
	   }
	   if(l<n) temp*=expdangle3[type][S1[l+1]];
	   else if (type>2 )  temp*=expTermAU;
	   pr[kl]+=temp;
	 }
       }
     }
     else if (l==cut_point ) {
       int t, sk,s;
       for (t=2; t<cut_point;t++) {
	 for (s=1; s<t; s++) {
	   for (k=cut_point; k<=n; k++) {
	     sk=iindx[s]-k;
	     if (qb[sk]) {
	       type=rtype[ptype[sk]];
	       temp=pr[sk]*expdangle3[type][S1[s+1]]*scale[2];
	       if (s+1<t) temp*=q[iindx[s+1]-(t-1)];
	       if (SAME_STRAND(k-1,k)) temp*=expdangle5[type][S1[k-1]]*q[iindx[cut_point]-(k-1)];
	       Qlout[t]+=temp;
	     }
	   }
	 }
       }
     }
     else if (l<cut_point) {
       for (k=1; k<l; k++) {
	 if (qb[iindx[k]-l]) {
	   type=ptype[iindx[k]-l];
	   temp=Qlout[k];
	   if (l+1<cut_point) temp*=q[iindx[l+1]-(cut_point-1)];
	   if (k>1) temp*=expdangle5[type][S1[k-1]];
	   if (l<(cut_point-1)) temp*=expdangle3[type][S1[l+1]];
	   else if (type>2) temp*=expTermAU;
	   pr[iindx[k]-l]+=temp;
	 }
       }
     }
#endif
#if 0
     /* old variant, p[lease delete me soon */
     if (l>cut_point) { /*right of cut*/
       /*k=cut_point*/
       kl=iindx[cut_point]-l;
       if (qb[kl]>0.) {
	 tt=ptype[kl];
	 /*i<cut_point-1*/
	 for (i=1; i<cut_point-1; i++) {
	   /*j==l+1*/
	   if (qb[iindx[i]-(l+1)]>0.) {
	     type=ptype[iindx[i]-(l+1)]; type=rtype[type];
	     pr[kl]+=q[iindx[i+1]-(cut_point-1)]*pr[iindx[i]-(l+1)]*
	       expdangle5[type][S1[l]]*expdangle3[type][S1[i+1]]*
	       expdangle3[tt][S1[l+1]]*scale[2];

	   }
	   /*j>l+1*/
	   for (j=l+2; j<=n; j++) {
	     if (qb[iindx[i]-(j)]>0.) {
	       type=ptype[iindx[i]-(j)]; type=rtype[type];
	       pr[kl]+=q[iindx[i+1]-(cut_point-1)]*pr[iindx[i]-(j)]*
		 expdangle5[type][S1[j-1]]*expdangle3[type][S1[i+1]]*
		 expdangle3[tt][S1[l+1]]*q[iindx[l+1]-(j-1)]*scale[2];

	     }
	   }
	 }
	 i=cut_point-1; /*i=cut_point-1;  */
	 /*j==l+1*/
	 if (qb[iindx[i]-(l+1)]>0.) {
	   type=ptype[iindx[i]-(l+1)]; type=rtype[type];
	   temp=pr[iindx[i]-(l+1)]*expdangle5[type][S1[l]]*
	     expdangle3[tt][S1[l+1]]*scale[2];
	   if (type>2) temp*=expTermAU;
	   pr[kl]+=temp;

	 }
	 /*j>l+1*/
	 for (j=l+2; j<=n; j++) {
	   if (qb[iindx[i]-(j)]>0.) {
	     type=ptype[iindx[i]-(j)]; type=rtype[type];
	     temp=pr[iindx[i]-(j)]*expdangle5[type][S1[j-1]]*
	       expdangle3[tt][S1[l+1]]*q[iindx[l+1]-(j-1)]*scale[2];
	     if (type>2) temp*=expTermAU;
	     pr[kl]+=temp;

	   }
	 }
       }

     for (k=cut_point+1; k<l;k++) { /*k>cut_point*/
	 kl=iindx[k]-l;
	 if (qb[kl]>0.) {
	   tt=ptype[kl];
	   /*i<cut_point-1*/
	   for (i=1; i<cut_point-1; i++) {
	     /*j==l+1*/
	     if (qb[iindx[i]-(l+1)]>0.) {
	       type=ptype[iindx[i]-(l+1)]; type=rtype[type];
	       pr[kl]+=q[iindx[cut_point]-(k-1)]*q[iindx[i+1]-(cut_point-1)]*
		 pr[iindx[i]-(l+1)]*expdangle5[type][S1[l]]*
		 expdangle3[type][S1[i+1]]*expdangle3[tt][S1[l+1]]*
		 expdangle5[tt][S1[k-1]]*scale[2];
	     }
	     /*j>l+1*/
	     for (j=l+2; j<=n; j++) {
	       if (qb[iindx[i]-(j)]>0.) {
		 type=ptype[iindx[i]-(j)]; type=rtype[type];
		 pr[kl]+=q[iindx[cut_point]-(k-1)]*q[iindx[i+1]-(cut_point-1)]*
		   pr[iindx[i]-(j)]*expdangle5[type][S1[j-1]]*
		   expdangle3[type][S1[i+1]]*expdangle3[tt][S1[l+1]]*
		   q[iindx[l+1]-(j-1)]*expdangle5[tt][S1[k-1]]*scale[2];
	       }
	     }
	   }
	   i=cut_point-1; /*i=cut_point-1;  */
	   /*j==l+1*/
	   if (qb[iindx[i]-(l+1)]>0.) {
	     type=ptype[iindx[i]-(l+1)]; type=rtype[type];
	     temp=q[iindx[cut_point]-(k-1)]*pr[iindx[i]-(l+1)]*
	       expdangle5[type][S1[l]]*expdangle3[tt][S1[l+1]]*
	       expdangle5[tt][S1[k-1]]*scale[2];
	     if (type>2) temp*=expTermAU;
	     pr[kl]+=temp;
	   }
	   /*j>l+1*/
	   for (j=l+2; j<=n; j++) {
	     if (qb[iindx[i]-(j)]>0.) {
	       type=ptype[iindx[i]-(j)]; type=rtype[type];
	       temp=q[iindx[cut_point]-(k-1)]*pr[iindx[i]-(j)]*
		 expdangle5[type][S1[j-1]]*expdangle3[tt][S1[l+1]]*
		 q[iindx[l+1]-(j-1)]*expdangle5[tt][S1[k-1]]*scale[2];
	       if (type>2) temp*=expTermAU;
	       pr[kl]+=temp;
	     }
	   }
	 }
       }
      }
      if (l<cut_point-1) { /*l bis 1 vor cut*/
     /*l left of cut*/
     for (k=2;k<l; k++) {

	 kl=iindx[k]-l;
	 if (qb[kl]>0.) {
	   tt=ptype[kl];                    /*i=k-1*/
	   if (qb[iindx[k-1]-cut_point]>0.) { /*j=cut_point*/
	     type=ptype[iindx[k-1]-cut_point];type=rtype[type];
	     pr[kl]+=pr[iindx[k-1]-cut_point]*q[iindx[l+1]-(cut_point-1)]*
	       expdangle3[tt][S1[l+1]]*expdangle5[tt][S1[k-1]]*
	       expdangle3[type][S1[k]]*scale[2];

	   }

	   for (j=cut_point+1; j<=n;j++) { /*j ab cut_point+1*/
	     if (qb[iindx[k-1]-j]>0.) {
	       type=ptype[iindx[k-1]-j];type=rtype[type];
	       pr[kl]+=pr[iindx[k-1]-j]*q[iindx[l+1]-(cut_point-1)]*
		 q[iindx[cut_point]-(j-1)]*expdangle3[tt][S1[l+1]]*
		 expdangle5[tt][S1[k-1]]*expdangle5[type][S1[j-1]]*
		 expdangle3[type][S1[k]]*scale[2];
	     }
	   }
	   /*i<k-1*/
	   for (i=1; i<k-1; i++) {
	     if (qb[iindx[i]-cut_point]>0.) { /*j=cut_point*/
	       type=ptype[iindx[i]-cut_point];type=rtype[type];
	       pr[kl]+=pr[iindx[i]-cut_point]*q[iindx[l+1]-(cut_point-1)]*
		 expdangle3[tt][S1[l+1]]*expdangle5[tt][S1[k-1]]*
		 expdangle3[type][S1[i+1]]*q[iindx[i+1]-(k-1)]*scale[2];
	     }

	     for (j=cut_point+1; j<=n;j++) { /*j ab cut_point+1*/
	       if (qb[iindx[i]-j]>0.) {
		 type=ptype[iindx[i]-j];type=rtype[type];
		 pr[kl]+=pr[iindx[i]-j]*q[iindx[l+1]-(cut_point-1)]*
		   q[iindx[cut_point]-(j-1)]*expdangle3[tt][S1[l+1]]*
		   expdangle5[tt][S1[k-1]]*expdangle5[type][S1[j-1]]*
		   expdangle3[type][S1[i+1]]*q[iindx[i+1]-(k-1)]*scale[2];
	       }
	     }
	   }
	 }
       }
      }
      /*l=cut_point-1*/
      if (l==cut_point-1) {
	for(k=2; k<l; k++) {
       kl=iindx[k]-l;
       if (qb[kl]>0.) {
	 tt=ptype[kl];                    /*i=k-1*/
	 if (qb[iindx[k-1]-cut_point]>0.) { /*j=cut_point*/
	   type=ptype[iindx[k-1]-cut_point];type=rtype[type];
	   temp=pr[iindx[k-1]-cut_point]*expdangle5[tt][S1[k-1]]*
	     expdangle3[type][S1[k]]*scale[2];
	   if (tt>2) temp*=expTermAU;
	   pr[kl]+=temp;
	 }

	 for (j=cut_point+1; j<=n;j++) { /*j ab cut_point+1*/
	   if (qb[iindx[k-1]-j]>0.) {
	     type=ptype[iindx[k-1]-j];type=rtype[type];
	     temp=pr[iindx[k-1]-j]*q[iindx[cut_point]-(j-1)]*
	       expdangle5[tt][S1[k-1]]*expdangle5[type][S1[j-1]]*
	       expdangle3[type][S1[k]]*scale[2];
	     if (tt>2) temp*=expTermAU;
	     pr[kl]+=temp;
	   }
	 }
	 /*i<k-1*/
	 for (i=1; i<k-1; i++) {
	   if (qb[iindx[i]-cut_point]>0.) { /*j=cut_point*/
	     type=ptype[iindx[i]-cut_point];type=rtype[type];
	     temp=pr[iindx[i]-cut_point]*expdangle5[tt][S1[k-1]]*
	       expdangle3[type][S1[i+1]]*q[iindx[i+1]-(k-1)]*scale[2];
	     if (tt>2) temp*=expTermAU;
	     pr[kl]+=temp;
	   }

	   for (j=cut_point+1; j<=n;j++) { /*j ab cut_point+1*/
	     if (qb[iindx[i]-j]>0.) {
	       type=ptype[iindx[i]-j];type=rtype[type];
	       temp=pr[iindx[i]-j]*q[iindx[cut_point]-(j-1)]*
		 expdangle5[tt][S1[k-1]]*expdangle5[type][S1[j-1]]*
		 expdangle3[type][S1[i+1]]*q[iindx[i+1]-(k-1)]*scale[2];
	       if (tt>2) temp*=expTermAU;
	       pr[kl]+=temp;
	     }
	   }
	 }
       }
	}
      }
#endif
    }  /* end for (l=..)   */
    free(Qlout);
    free(Qrout);
    for (i=1; i<=n; i++)
      for (j=i+TURN+1; j<=n; j++) {
	ij = iindx[i]-j;
	pr[ij] *= qb[ij];
      }

    if (structure!=NULL)
      sprintf_bppm(n, structure);
  }   /* end if (do_backtrack)*/

  if (ov>0) fprintf(stderr, "%d overflows occurred while backtracking;\n"
		    "you might try a smaller pf_scale than %g\n",
		    ov, pf_scale);
  return X;
}

/*------------------------------------------------------------------------*/
/* dangling ends should never be destabilizing, i.e. expdangle>=1         */
/* specific heat needs smooth function (2nd derivative)                   */
/* we use a*(sin(x+b)+1)^2, with a=2/(3*sqrt(3)), b=Pi/6-sqrt(3)/2,       */
/* in the interval b<x<sqrt(3)/2                                          */

#define SCALE 10
/* #define SMOOTH(X) ((X)/SCALE<-1.2283697)?0:(((X)/SCALE>0.8660254)?(X):\
  SCALE*0.38490018*(sin((X)/SCALE-0.34242663)+1)*(sin((X)/SCALE-0.34242663)+1))*/
#define SMOOTH(X) ((X)<0 ? 0 : (X))

PRIVATE void scale_pf_params(unsigned int length)
{
  /* scale energy parameters and pre-calculate Boltzmann weights */
  unsigned int i, j, k, l;
  double  kT, TT;
  double  GT;



  kT = (temperature+K0)*GASCONST;   /* kT in cal/mol  */
  TT = (temperature+K0)/(Tmeasure);

   /* scaling factors (to avoid overflows) */
  if (pf_scale==-1) { /* mean energy for random sequences: 184.3*length cal */
    pf_scale = exp(-(-185+(temperature-37.)*7.27)/kT);
    if (pf_scale<1) pf_scale=1;
  }
  scale[0] = 1.;
  for (i=1; i<=length; i++) {
    scale[i] = scale[i-1]/pf_scale;
  }

  /* loop energies: hairpins, bulges, interior, mulit-loops */
  for (i=0; i<=MIN(30,length); i++) {
    GT =  hairpin37[i]*TT;
    exphairpin[i] = exp( -GT*10./kT);
  }
  for (i=0; i<=MIN(30, MAXLOOP); i++) {
    GT =  bulge37[i]*TT;
    expbulge[i] = exp( -GT*10./kT);
    GT =  internal_loop37[i]*TT;
    expinternal[i] = exp( -GT*10./kT);
  }
  /* special case of size 2 interior loops (single mismatch) */
  if (james_rule) expinternal[2] = exp ( -80*10/kT);

  lxc = lxc37*TT;
  for (i=31; i<length; i++) {
    GT = hairpin37[30]*TT + (lxc*log( i/30.));
    exphairpin[i] = exp( -GT*10./kT);
  }
  for (i=31; i<=MAXLOOP; i++) {
    GT = bulge37[30]*TT + (lxc*log( i/30.));
    expbulge[i] = exp( -GT*10./kT);
    GT = internal_loop37[30]*TT + (lxc*log( i/30.));
    expinternal[i] = exp( -GT*10./kT);
  }

  for (i=0; i<5; i++) {
    GT = F_ninio37[i]*TT;
    for (j=0; j<=MAXLOOP; j++)
      expninio[i][j]=exp(-MIN(MAX_NINIO,j*GT)*10/kT);
  }
  for (i=0; (i*7)<strlen(Tetraloops); i++) {
    GT = TETRA_ENTH37 - (TETRA_ENTH37-TETRA_ENERGY37[i])*TT;
    exptetra[i] = exp( -GT*10./kT);
  }
  for (i=0; (i*5)<strlen(Triloops); i++)
    expTriloop[i] = exp(-Triloop_E37[i]*10/kT);

  GT =  ML_closing37*TT;
  expMLclosing = exp( -GT*10/kT);

  for (i=0; i<=NBPAIRS; i++) { /* includes AU penalty */
    GT =  ML_intern37*TT;
    /* if (i>2) GT += TerminalAU; */
    expMLintern[i] = exp( -GT*10./kT);
  }
  expTermAU = exp(-TerminalAU*10/kT);

  GT =  ML_BASE37*TT;
  for (i=0; i<length; i++) {
    expMLbase[i] = exp( -10.*i*GT/kT)*scale[i];
  }

  /* if dangles==0 just set their energy to 0,
      don't let dangle energies become > 0 (at large temps),
      but make sure go smoothly to 0                        */
  for (i=0; i<=NBPAIRS; i++)
    for (j=0; j<=4; j++) {
      if (dangles) {
	GT = dangle5_H[i][j] - (dangle5_H[i][j] - dangle5_37[i][j])*TT;
	expdangle5[i][j] = exp(SMOOTH(-GT)*10./kT);
	GT = dangle3_H[i][j] - (dangle3_H[i][j] - dangle3_37[i][j])*TT;
	expdangle3[i][j] = exp(SMOOTH(-GT)*10./kT);
      } else
	expdangle3[i][j] = expdangle5[i][j] = 1;
      if (i>2) /* add TermAU penalty into dangle3 */
	expdangle3[i][j] *= expTermAU;
    }

  /* stacking energies */
  for (i=0; i<=NBPAIRS; i++)
    for (j=0; j<=NBPAIRS; j++) {
      GT =  enthalpies[i][j] - (enthalpies[i][j] - stack37[i][j])*TT;
      expstack[i][j] = exp( -GT*10/kT);
    }

  /* mismatch energies */
  for (i=0; i<=NBPAIRS; i++)
    for (j=0; j<5; j++)
      for (k=0; k<5; k++) {
	GT = mism_H[i][j][k] - (mism_H[i][j][k] - mismatchI37[i][j][k])*TT;
	expmismatchI[i][j][k] = exp(-GT*10.0/kT);
	GT = mism_H[i][j][k] - (mism_H[i][j][k] - mismatchH37[i][j][k])*TT;
	expmismatchH[i][j][k] = exp(-GT*10.0/kT);
	GT = mism_H[i][j][k] - (mism_H[i][j][k] - mismatchM37[i][j][k])*TT;
	expmismatchM[i][j][k] = exp(-GT*10.0/kT);
      }

  /* interior lops of length 2 */
  for (i=0; i<=NBPAIRS; i++)
    for (j=0; j<=NBPAIRS; j++)
      for (k=0; k<5; k++)
	for (l=0; l<5; l++) {
	  GT = int11_H[i][j][k][l] -
	    (int11_H[i][j][k][l] - int11_37[i][j][k][l])*TT;
	  expint11[i][j][k][l] = exp(-GT*10./kT);
	}
  /* interior 2x1 loops */
  for (i=0; i<=NBPAIRS; i++)
    for (j=0; j<=NBPAIRS; j++)
      for (k=0; k<5; k++)
	for (l=0; l<5; l++) {
	  int m;
	  for (m=0; m<5; m++) {
	    GT = int21_H[i][j][k][l][m] -
	      (int21_H[i][j][k][l][m] - int21_37[i][j][k][l][m])*TT;
	    expint21[i][j][k][l][m] = exp(-GT*10./kT);
	  }
	}
  /* interior 2x2 loops */
  for (i=0; i<=NBPAIRS; i++)
    for (j=0; j<=NBPAIRS; j++)
      for (k=0; k<5; k++)
	for (l=0; l<5; l++) {
	  int m,n;
	  for (m=0; m<5; m++)
	    for (n=0; n<5; n++) {
	      GT = int22_H[i][j][k][l][m][n] -
		(int22_H[i][j][k][l][m][n]-int22_37[i][j][k][l][m][n])*TT;
	      expint22[i][j][k][l][m][n] = exp(-GT*10./kT);
	    }
	}
  /*initialization of duplex energy*/
  expInit=exp(-DuplexInit*TT*10./kT);
}

/*----------------------------------------------------------------------*/
PRIVATE double expHairpinEnergy(int u, int type, short si1, short sj1,
				const char *string) {
  double q;
  q = exphairpin[u];
  if ((tetra_loop)&&(u==4)) {
    char tl[7]={0}, *ts;
    strncpy(tl, string, 6);
    if ((ts=strstr(Tetraloops, tl)))
      q *= exptetra[(ts-Tetraloops)/7];
  }
  if (u==3) {
    char tl[6]="", *ts;
    strncpy(tl, string, 5);
    if ((ts=strstr(Triloops, tl)))
      q *= expTriloop[(ts-Triloops)/6];
    if (type>2)
      q *= expTermAU;
  }
  else /* no mismatches for tri-loops */
    q *= expmismatchH[type][si1][sj1];

  q *= scale[u+2];
  return q;
}

PRIVATE double expLoopEnergy(int u1, int u2, int type, int type2,
			     short si1, short sj1, short sp1, short sq1) {
  double z=0;
  int no_close = 0;

  if ((no_closingGU) && ((type2==3)||(type2==4)||(type==2)||(type==4)))
    no_close = 1;

  if ((u1==0) && (u2==0)) /* stack */
    z = expstack[type][type2];
  else if (no_close==0) {
    if ((u1==0)||(u2==0)) { /* bulge */
      int u;
      u = (u1==0)?u2:u1;
      z = expbulge[u];
      if (u2+u1==1) z *= expstack[type][type2];
      else {
	if (type>2) z *= expTermAU;
	if (type2>2) z *= expTermAU;
      }
    }
    else {     /* interior loop */
      if (u1+u2==2) /* size 2 is special */
	z = expint11[type][type2][si1][sj1];
      else if ((u1==1) && (u2==2))
	z = expint21[type][type2][si1][sq1][sj1];
      else if ((u1==2) && (u2==1))
	z = expint21[type2][type][sq1][si1][sp1];
      else if ((u1==2) && (u2==2))
	z = expint22[type][type2][si1][sp1][sq1][sj1];
      else {
	z = expinternal[u1+u2]*
	  expmismatchI[type][si1][sj1]*
	  expmismatchI[type2][sq1][sp1];
	z *= expninio[2][abs(u1-u2)];
      }
    }
  }
  return z*scale[u1+u2+2];
}

/*----------------------------------------------------------------------*/

PRIVATE void get_arrays(unsigned int length)
{
  unsigned int size,i;

  size = sizeof(FLT_OR_DBL) * ((length+1)*(length+2)/2);
  q   = (FLT_OR_DBL *) space(size);
  qb  = (FLT_OR_DBL *) space(size);
  qm  = (FLT_OR_DBL *) space(size);
  pr  = (FLT_OR_DBL *) space(size); /*q is needed later*/

  qm1 = (FLT_OR_DBL *) space(size);

  ptype = (char *) space(sizeof(char)*((length+1)*(length+2)/2));
  q1k = (FLT_OR_DBL *) space(sizeof(FLT_OR_DBL)*(length+1));

  qln = (FLT_OR_DBL *) space(sizeof(FLT_OR_DBL)*(length+2));
  qq  = (FLT_OR_DBL *) space(sizeof(FLT_OR_DBL)*(length+2));
  qq1 = (FLT_OR_DBL *) space(sizeof(FLT_OR_DBL)*(length+2));
  qqm  = (FLT_OR_DBL *) space(sizeof(FLT_OR_DBL)*(length+2));
  qqm1 = (FLT_OR_DBL *) space(sizeof(FLT_OR_DBL)*(length+2));
  prm_l = (FLT_OR_DBL *) space(sizeof(FLT_OR_DBL)*(length+2));
  prm_l1 =(FLT_OR_DBL *) space(sizeof(FLT_OR_DBL)*(length+2));
  prml = (FLT_OR_DBL *) space(sizeof(FLT_OR_DBL)*(length+2));
  exphairpin = (FLT_OR_DBL *) space(sizeof(FLT_OR_DBL)*(length+1));
  expMLbase  = (FLT_OR_DBL *) space(sizeof(FLT_OR_DBL)*(length+1));
  scale = (FLT_OR_DBL *) space(sizeof(FLT_OR_DBL)*(length+1));
  iindx = (int *) space(sizeof(int)*(length+1));
  jindx = (int *) space(sizeof(int)*(length+1));
  for (i=1; i<=length; i++) {
    iindx[i] = ((length+1-i)*(length-i))/2 +length+1;
    jindx[i] = (i*(i-1))/2;
  }
}

/*----------------------------------------------------------------------*/

PUBLIC void init_co_pf_fold(int length)
{
  if (length<1) nrerror("init_pf_fold: length must be greater 0");
  if (init_length>0) free_co_pf_arrays(); /* free previous allocation */
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
  init_length=length;
}

PUBLIC void free_co_pf_arrays(void)
{
  free(q);
  free(qb);
  free(qm);
  free(pr);
  if (qm1 != NULL) {free(qm1); qm1 = NULL;}
  free(ptype);
  free(qq); free(qq1);
  free(qqm); free(qqm1);
  free(q1k); free(qln);
  free(prm_l); free(prm_l1); free(prml);
  free(exphairpin);
  free(expMLbase);
    free(scale);
  free(iindx);
  free(jindx);
#ifdef SUN4
  standard_arithmetic();
#else
#ifdef HP9
  fpsetfastmode(0);
#endif
#endif
  init_length=0;
  free(S); S=NULL;
  free(S1); S1=NULL;
}
/*---------------------------------------------------------------------------*/

PUBLIC void update_co_pf_params(int length)
{
  make_pair_matrix();
  scale_pf_params((unsigned) length);
}

/*---------------------------------------------------------------------------*/

PUBLIC char co_bppm_symbol(float *x)
{
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

/*---------------------------------------------------------------------------*/
#define L 3
PRIVATE void sprintf_bppm(int length, char *structure)
{
  int    i,j;
  float  P[L];   /* P[][0] unpaired, P[][1] upstream p, P[][2] downstream p */

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
    structure[j-1] = co_bppm_symbol(P);
  }
  structure[length] = '\0';
}

/*---------------------------------------------------------------------------*/
PRIVATE void make_ptypes(const short *S, const char *structure) {
  int n,i,j,k,l;

  n=S[0];
  for (k=1; k<=n-TURN-1; k++)
    for (l=1; l<=2; l++) {
      int type,ntype=0,otype=0;
      i=k; j = i+TURN+l;
      if (j>n) continue;
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

  if (fold_constrained&&(structure!=NULL)) {
    int hx, *stack;
    char type;
    stack = (int *) space(sizeof(int)*(n+1));

    for(hx=0, j=1; j<=n; j++) {
      switch (structure[j-1]) {
      case 'x': /* can't pair */
	for (l=1; l<j-TURN; l++) ptype[iindx[l]-j] = 0;
	for (l=j+TURN+1; l<=n; l++) ptype[iindx[j]-l] = 0;
	break;
      case '(':
	stack[hx++]=j;
	/* fallthrough */
      case '<': /* pairs upstream */
	for (l=1; l<j-TURN; l++) ptype[iindx[l]-j] = 0;
	break;
      case ')':
	if (hx<=0) {
	  fprintf(stderr, "%s\n", structure);
	  nrerror("unbalanced brackets in constraints");
	}
	i = stack[--hx];
	type = ptype[iindx[i]-j];
	/* don't allow pairs i<k<j<l */
	for (k=i; k<=j; k++)
	  for (l=j; l<=n; l++) ptype[iindx[k]-l] = 0;
	/* don't allow pairs k<i<l<j */
	for (k=1; k<=i; k++)
	  for (l=i; l<=j; l++) ptype[iindx[k]-l] = 0;
	ptype[iindx[i]-j] = (type==0)?7:type;
	/* fallthrough */
      case '>': /* pairs downstream */
	for (l=j+TURN+1; l<=n; l++) ptype[iindx[j]-l] = 0;
	break;
      case 'l': /*only intramolecular basepairing*/
	if (j<cut_point) for (l=cut_point; l<=n; l++) ptype[iindx[j]-l] = 0;
	else for (l=1; l<cut_point; l++) ptype[iindx[l]-j] =0;
	break;
      case 'e': /*only intermolecular bp*/
	if (j<cut_point) {
	  for (l=1; l<j; l++) ptype[iindx[l]-j] =0;
	  for (l=j+1; l<cut_point; l++) ptype[iindx[j]-l] = 0;
	}
	else {
	  for (l=cut_point; l<j; l++) ptype[iindx[l]-j] =0;
	  for (l=j+1; l<=n; l++) ptype[iindx[j]-l] = 0;
	}
	break;
      }
    }
    if (hx!=0) {
      fprintf(stderr, "%s\n", structure);
      nrerror("unbalanced brackets in constraint string");
    }
    free(stack);
  }
  if (mirnatog==1) {   /*microRNA toggle: no intramolec. bp in 2. molec*/
    for (j=cut_point; j<n; j++) {
      for (l=j+1; l<=n; l++) {
	ptype[iindx[j]-l] = 0;
      }
    }
    }
}

/*
  stochastic backtracking in pf_fold arrays
  returns random structure S with Boltzman probabilty
  p(S) = exp(-E(S)/kT)/Z
*/
static void backtrack(int i, int j);

static char *pstruc;
static char *sequence;

static void backtrack_qm1(int i,int j) {
  /* i is paired to l, i<l<j; backtrack in qm1 to find l */
  int ii, l, type;
  double qt, r;

  r = urn() * qm1[jindx[j]+i];
  ii = iindx[i];
  for (qt=0., l=i+TURN+1; l<=j; l++) {
    type = ptype[ii-l];
    if (type)
      qt +=  qb[ii-l]*expMLintern[type]*
	expdangle5[type][S1[i-1]] * expdangle3[type][S1[l+1]] *
	expMLbase[j-l];
    if (qt>=r) break;
  }
  if (l>j) nrerror("backtrack failed in qm1");
  backtrack(i,l);
}

static void backtrack(int i, int j) {
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
      qbt1 = expHairpinEnergy(u, type, S1[i+1], S1[j-1], sequence+i-1);

    if (qbt1>r) return; /* found the hairpin we're done */

    for (k=i+1; k<=MIN(i+MAXLOOP+1,j-TURN-2); k++) {
      u1 = k-i-1;
      for (l=MAX(k+TURN+1,j-1-MAXLOOP+u1); l<j; l++) {
	int type_2;
	type_2 = ptype[iindx[k]-l];
	if (type_2) {
	  type_2 = rtype[type_2];
	  qbt1 += qb[iindx[k]-l] *
	    expLoopEnergy(u1, j-l-1, type, type_2,
			  S1[i+1], S1[j-1], S1[k-1], S1[l+1]);
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
    while (j>i) {
      /* now backtrack  [i ... j] in qm[] */
      jj = jindx[j];
      ii = iindx[i];
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

  mykT=(temperature+K0)*GASCONST/1000.;

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

  kT=(temperature+K0)*GASCONST/1000.;
  Concentration=(struct ConcEnt *)space(20*sizeof(struct ConcEnt));
 /* Compute equilibrium constants */
  /* again note the input free energies are not from the null model (without DuplexInit) */

  KAA = exp(( 2.0 * FEA - FcAA)/kT);
  KBB = exp(( 2.0 * FEB - FcBB)/kT);
  KAB = exp(( FEA + FEB - FcAB)/kT);
  printf("Kaa..%g %g %g\n", KAA, KBB, KAB);
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

PUBLIC struct plist *get_plist(struct plist *pl, int length, double cut_off) {
  int i, j,n, count;
  /*get pair probibilities out of pr array*/
  count=0;
  n=2;
  for (i=1; i<length; i++) {
    for (j=i+1; j<=length; j++) {
      if (pr[iindx[i]-j]<cut_off) continue;
      if (count==n*length-1) {
	n*=2;
	pl=(struct plist *)xrealloc(pl,n*length*sizeof(struct plist));
      }
      pl[count].i=i;
      pl[count].j=j;
      pl[count++].p=pr[iindx[i]-j];
      /*      printf("gpl: %2d %2d %.9f\n",i,j,pr[iindx[i]-j]);*/
    }
  }
  pl[count].i=0;
  pl[count].j=0; /*->??*/
  pl[count++].p=0.;
  pl=(struct plist *)xrealloc(pl,(count)*sizeof(struct plist));
  return pl;
}


#if 0
PUBLIC int make_probsum(int length, char *name) {
  /*compute probability of any base to be paired (preliminary)*/
  double *Spprob;
  double *Pprob;
  int i,j;
  FILE *fp;
  char *filename;
  filename=(char *)space((strlen(name)+10)*sizeof(char));
  sprintf(filename,"%s_pp.dat",name);
  fp=fopen(filename,"w");
  Spprob=(double *)space((length+1)*sizeof(double));
  if (cut_point>0) Pprob=(double *)space((length+1)*sizeof(double));
  for (i=1; i<length; i++) {
    for (j=i+1; j<=length; j++) {
      Spprob[i]+=pr[iindx[i]-j];
      Spprob[j]+=pr[iindx[i]-j];
      if (!SAME_STRAND(i,j)) {
	Pprob[i]+=pr[iindx[i]-j];
	Pprob[j]+=pr[iindx[i]-j];
      }
    }
  }

  /*  plot_probsum(Spprob,Pprop, seq);*/
  /*daweilamal:*/
  for (i=1; i<=length; i++) {
    fprintf(fp,"%4d %.8f\n",i,Spprob[i]);
  }
  fprintf(fp,"&\n");
  if (cut_point>0) for (i=1;i<=length; i++) {
    fprintf(fp,"%4d %.8f\n",i,Pprob[i]);
  }
  fclose(fp);
  free(Spprob);
  free(filename);
  if (cut_point>0) free(Pprob);
  return 1;
}
#endif
