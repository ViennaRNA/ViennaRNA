/* Last changed Time-stamp: <2004-05-14 18:28:57 ivo> */
/*                
		  partiton function for RNA secondary structures

		  Ivo L Hofacker
		  Vienna RNA package
*/
/*
  $Log: part_func.c,v $
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

/*@unused@*/
static char rcsid[] UNUSED = "$Id: part_func.c,v 1.20 2004/08/12 12:14:46 ivo Exp $";

#define MAX(x,y) (((x)>(y)) ? (x) : (y))
#define MIN(x,y) (((x)<(y)) ? (x) : (y))
#define PUBLIC
#define PRIVATE static

PUBLIC  float pf_fold(char *sequence, char *structure);
PUBLIC  void  init_pf_fold(int length);
PUBLIC  void  free_pf_arrays(void);
PUBLIC  void  update_pf_params(int length);
PUBLIC  char  bppm_symbol(float *x);
PUBLIC  int   st_back=0;
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
PRIVATE FLT_OR_DBL expinternal[MAXLOOP+1];
PRIVATE FLT_OR_DBL expninio[5][MAXLOOP+1];
PRIVATE FLT_OR_DBL *q, *qb, *qm, *qm1, *qqm, *qqm1, *qq, *qq1;
PRIVATE FLT_OR_DBL *prml, *prm_l, *prm_l1, *q1k, *qln;
PRIVATE FLT_OR_DBL *scale;
PRIVATE char *ptype; /* precomputed array of pair types */ 
PRIVATE int *jindx;
PRIVATE int init_length;  /* length in last call to init_pf_fold() */
PRIVATE double init_temp; /* temperature in last call to scale_pf_params */
#define ISOLATED  256.0

/*-----------------------------------------------------------------*/
static  short *S, *S1;
PUBLIC float pf_fold(char *sequence, char *structure)
{

  int n, i,j,k,l, ij, kl, u,u1,d,ii,ll, type, type_2, tt, ov=0;
  FLT_OR_DBL temp, Q, Qmax=0, prm_MLb;
  FLT_OR_DBL prmt,prmt1;
  FLT_OR_DBL qbt1, *tmp;
   
  double free_energy;
  double max_real;

  max_real = (sizeof(FLT_OR_DBL) == sizeof(float)) ? FLT_MAX : DBL_MAX;

  n = (int) strlen(sequence);
  if (n>init_length) init_pf_fold(n);  /* (re)allocate space */
  if ((init_temp - temperature)>1e-6) update_pf_params(n);

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

  for (d=0; d<=TURN; d++) 
    for (i=1; i<=n-d; i++) {
      j=i+d;
      ij = iindx[i]-j;
      q[ij]=1.0*scale[d+1];
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
      if (type!=0) {
	/*hairpin contribution*/
	if (((type==3)||(type==4))&&no_closingGU) qbt1 = 0;
	else 
	  qbt1 = expHairpinEnergy(u, type, S1[i+1], S1[j-1], sequence+i-1);
	
	/* interior loops with interior pair k,l */
	for (k=i+1; k<=MIN(i+MAXLOOP+1,j-TURN-2); k++) {
	  u1 = k-i-1;
	  for (l=MAX(k+TURN+1,j-1-MAXLOOP+u1); l<j; l++) {
	    type_2 = ptype[iindx[k]-l];
	    if (type_2) {
	      type_2 = rtype[type_2];
	      qbt1 += qb[iindx[k]-l] * 
		expLoopEnergy(u1, j-l-1, type, type_2,
			      S1[i+1], S1[j-1], S1[k-1], S1[l+1]);
	    }
	  }
	}
	/*multiple stem loop contribution*/
	ii = iindx[i+1]; /* ii-k=[i+1,k-1] */
	temp = 0.0;
	for (k=i+2; k<=j-1; k++) temp += qm[ii-(k-1)]*qqm1[k]; 
	tt = rtype[type];
	qbt1 += temp*expMLclosing*expMLintern[tt]*scale[2]*
	  expdangle3[tt][S1[i+1]]*expdangle5[tt][S1[j-1]];
	 
	qb[ij] = qbt1;
      } /* end if (type!=0) */
      else qb[ij] = 0.0;
       
      /* construction of qqm matrix containing final stem
	 contributions to multiple loop partition function
	 from segment i,j */
      qqm[i] = qqm1[i]*expMLbase[1];
      if (type) {
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
      ii = iindx[i];  /* ii-k=[i,k-1] */
      for (k=i+1; k<=j; k++) temp += (qm[ii-(k-1)]+expMLbase[k-i])*qqm[k];
      qm[ij] = (temp + qqm[i]);
      
      /*auxiliary matrix qq for cubic order q calculation below */
      qbt1 = qb[ij];
      if (type) {
	if (i>1) qbt1 *= expdangle5[type][S1[i-1]];
	if (j<n) qbt1 *= expdangle3[type][S1[j+1]];
	else if (type>2) qbt1 *= expTermAU;
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
      
  /* backtracking to construct binding probabilities of pairs*/
   
  if (do_backtrack) {
    Qmax=0;

    for (k=1; k<=n; k++) {
      q1k[k] = q[iindx[1] - k];
      qln[k] = q[iindx[k] -n];
    }
    q1k[0] = 1.0;
    qln[n+1] = 1.0;
      
    pr = q;     /* recycling */

    /* 1. exterior pair i,j and initialization of pr array */
    for (i=1; i<=n; i++) {
      for (j=i; j<=MIN(i+TURN,n); j++) pr[iindx[i]-j] = 0;
      for (j=i+TURN+1; j<=n; j++) {
	ij = iindx[i]-j;
	type = ptype[ij];
	if (type&&(qb[ij]>0.)) {
	  pr[ij] = q1k[i-1]*qln[j+1]/q1k[n];
	  if (i>1) pr[ij] *= expdangle5[type][S1[i-1]];
	  if (j<n) pr[ij] *= expdangle3[type][S1[j+1]];
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
	    ij = iindx[i] - j;
	    type = ptype[ij];
	    if ((pr[ij]>0)) {
	      pr[kl] += pr[ij]*expLoopEnergy(k-i-1, j-l-1, type, type_2,
					     S1[i+1], S1[j-1], S1[k-1], S1[l+1]);
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
	prmt1 = pr[ii-(l+1)]*expMLclosing*expMLintern[tt]*
	  expdangle3[tt][S1[i+1]]*expdangle5[tt][S1[l]];
	for (j=l+2; j<=n; j++) {
	  tt = ptype[ii-j]; tt = rtype[tt];
	  prmt += pr[ii-j]*expdangle3[tt][S1[i+1]]*
	    expdangle5[tt][S1[j-1]] *qm[ll-(j-1)];
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

	for (i=1;i<=k-2; i++) 
	  temp += prml[i]*qm[iindx[i+1] - (k-1)];

	temp *= expMLintern[tt]*scale[2];
	if (k>1) temp *= expdangle5[tt][S1[k-1]];
	if (l<n) temp *= expdangle3[tt][S1[l+1]];
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

      } /* end for (k=..) */
      tmp = prm_l1; prm_l1=prm_l; prm_l=tmp;

    }  /* end for (l=..)   */
      
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
   
  return free_energy; 
}

/*------------------------------------------------------------------------*/
/* dangling ends should never be destabilizing, i.e. expdangle>=1         */
/* specific heat needs smooth function (2nd derivative)                   */
/* we use a*(sin(x+b)+1)^2, with a=2/(3*sqrt(3)), b=Pi/6-sqrt(3)/2,       */
/* in the interval b<x<sqrt(3)/2                                          */

#define SCALE 10
#define SMOOTH(X) ((X)/SCALE<-1.2283697)?0:(((X)/SCALE>0.8660254)?(X):\
          SCALE*0.38490018*(sin((X)/SCALE-0.34242663)+1)*(sin((X)/SCALE-0.34242663)+1))
/* #define SMOOTH(X) ((X)<0 ? 0 : (X)) */

PRIVATE void scale_pf_params(unsigned int length)
{
  /* scale energy parameters and pre-calculate Boltzmann weights */
  unsigned int i, j, k, l;
  double  kT, TT;
  double  GT;

   
  init_temp = temperature;
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
	expdangle3[i][j] =  exp(SMOOTH(-GT)*10./kT);
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
    char tl[6]={0}, *ts;
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

  if (st_back) {
    qm1 = (FLT_OR_DBL *) space(size);
  }
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
   
PUBLIC void init_pf_fold(int length)
{
  if (length<1) nrerror("init_pf_fold: length must be greater 0");
  if (init_length>0) free_pf_arrays(); /* free previous allocation */
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

PUBLIC void free_pf_arrays(void)
{
  free(q); q=pr=NULL;
  free(qb);
  free(qm);
  if (qm1 != NULL) {free(qm1); qm1 = NULL;}
  free(ptype);
  free(qq); free(qq1);
  free(qqm); free(qqm1);
  free(q1k); free(qln);
  free(prm_l); free(prm_l1); free(prml);
  free(exphairpin);
  free(expMLbase);
  free(scale);
  free(iindx); free(jindx);
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

PUBLIC void update_pf_params(int length)
{
  if (length>init_length) init_pf_fold(length);  /* init not update */
  else {
    make_pair_matrix();
    scale_pf_params((unsigned) length);
  }
}

/*---------------------------------------------------------------------------*/

PUBLIC char bppm_symbol(float *x)
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
    structure[j-1] = bppm_symbol(P);
  }
  structure[length] = '\0';
}   

/*---------------------------------------------------------------------------*/
PRIVATE void make_ptypes(const short *S, const char *structure) {
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
      }
    }
    if (hx!=0) {
      fprintf(stderr, "%s\n", structure);
      nrerror("unbalanced brackets in constraint string");
    }
    free(stack);
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
char *pbacktrack(char *seq) {
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
	if (i>1) qkl *= expdangle5[type][S1[i-1]];
	if (j<n) qkl *= expdangle3[type][S1[j+1]];
	else if (type>2) qkl *= expTermAU;
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
    int k, l, ii, jj, type;

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
  return d;
}
