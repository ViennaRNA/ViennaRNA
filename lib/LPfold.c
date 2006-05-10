/* Last changed Time-stamp: <2006-05-08 15:55:54 ivo> */
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
/*@unused@*/
static char rcsid[] UNUSED = "$Id: LPfold.c,v 1.3 2006/05/10 15:10:01 ivo Exp $";

#define MAX(x,y) (((x)>(y)) ? (x) : (y))
#define MIN(x,y) (((x)<(y)) ? (x) : (y))
#define PUBLIC
#define PRIVATE static

static int num_p=0; /*for counting basepairs*/
PUBLIC  void  init_pf_foldLP(int length);
PUBLIC  void  free_pf_arraysLP(void);
PUBLIC  void  update_pf_paramsLP(int length);
/*PUBLIC  int   st_back=0;*/
PRIVATE void  scale_pf_params(unsigned int length);
PRIVATE void  get_arrays(unsigned int length);
PRIVATE double expLoopEnergy(int u1, int u2, int type, int type2,
			     short si1, short sj1, short sp1, short sq1);
PRIVATE double expHairpinEnergy(int u, int type, short si1, short sj1,
				const char *string);
/*PRIVATE void make_ptypes(const short *S, const char *structure);*/
/*new functions*/
PRIVATE void GetPtype(int j, int pairsize, const short *S, int n);
PRIVATE void FreeOldArrays(int i);
PRIVATE void GetNewArrays(int j, int winSize);
PRIVATE void printpbar(FLT_OR_DBL **prb,int winSize, int i, int n, float cutoff);
PRIVATE struct plist *get_plistW(struct plist *pl, int length, double cutoff, int start, FLT_OR_DBL **Tpr, int winSize);
/*end*/
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
PRIVATE FLT_OR_DBL **q, **qb, **qm, *qm1, *qqm, *qqm1, *qq, *qq1, **pR;
PRIVATE FLT_OR_DBL *prml, *prm_l, *prm_l1, *q1k, *qln;
PRIVATE FLT_OR_DBL *scale;
PRIVATE char **ptype; /* precomputed array of pair types */
PRIVATE int *jindx;
PRIVATE int init_length;  /* length in last call to init_pf_fold() */
PRIVATE double init_temp; /* temperature in last call to scale_pf_params */
#define ISOLATED  256.0

/*-----------------------------------------------------------------*/
static  short *S, *S1;
PUBLIC int pfl_fold(char *sequence, int winSize, int pairSize, float cutoff, struct plist **pl)
{

  int n, m,i,j,k,l, u,u1,ii, type, type_2, tt, ov=0;
  FLT_OR_DBL temp, Qmax=0, prm_MLb;
  FLT_OR_DBL prmt,prmt1;
  FLT_OR_DBL qbt1, *tmp;

  double free_energy;
  double max_real;

  max_real = (sizeof(FLT_OR_DBL) == sizeof(float)) ? FLT_MAX : DBL_MAX;

  n = (int) strlen(sequence);
  if (n<TURN+2) return 0;
  if (n>init_length) init_pf_foldLP(n);  /* (re)allocate space */
  if ((init_temp - temperature)>1e-6) update_pf_paramsLP(n);

  S = (short *) xrealloc(S, sizeof(short)*(n+1));
  S1= (short *) xrealloc(S1, sizeof(short)*(n+1));
  S[0] = n;
  for (l=1; l<=n; l++) {
    S[l]  = (short) encode_char(toupper(sequence[l-1]));
    S1[l] = alias[S[l]];
  }
  /*  make_ptypes(S, structure); das machmadochlieber lokal, ey!*/

  /*array initialization ; qb,qm,q
    qb,qm,q (i,j) are stored as ((n+1-i)*(n-i) div 2 + n+1-j */
  num_p=0;
  *pl=space(1000*sizeof(struct plist));


  /*  for (i=1; i<=n; i++)
      qq[i]=qq1[i]=qqm[i]=qqm1[i]=prm_l[i]=prm_l1[i]=prml[i]=0; oisa gaunza unnedig, wenn nicht recycelt*/
  /*ALWAYS q[i][j] => i>j!!*/
  for (j=1; j<MIN(TURN+2,n); j++) { /*allocate start*/
    GetNewArrays(j, winSize);
    GetPtype(j,pairSize,S,n);
    for (i=1; i<=j; i++) {
      q[i][j]=scale[(j-i+1)];
    }
  }
  for (j=TURN+2;j<=n+winSize; j++) {
    if (j<=n) {
      GetNewArrays(j, winSize);
      GetPtype(j,pairSize,S,n);
      for (i=MAX(1,j-winSize); i<=j/*-TURN*/; i++) {
	q[i][j]=scale[(j-i+1)];
      }
      for (i=j-TURN-1;i>=MAX(1,(j-winSize+1)); i--) {
	/* construction of partition function of segment i,j*/
	/*firstly that given i bound to j : qb(i,j) */
	u = j-i-1;
	type = ptype[i][j];
	if (type!=0) {
	  /*hairpin contribution*/
	  if (((type==3)||(type==4))&&no_closingGU) qbt1 = 0;
	  else
	    qbt1 = expHairpinEnergy(u, type, S1[i+1], S1[j-1], sequence+i-1);

	  /* interior loops with interior pair k,l */
	  for (k=i+1; k<=MIN(i+MAXLOOP+1,j-TURN-2); k++) {
	    u1 = k-i-1;
	    for (l=MAX(k+TURN+1,j-1-MAXLOOP+u1); l<j; l++) {
	      type_2 = ptype[k][l];
	      if (type_2) {
		type_2 = rtype[type_2];
		qbt1 += qb[k][l] *
		  expLoopEnergy(u1, j-l-1, type, type_2,
				S1[i+1], S1[j-1], S1[k-1], S1[l+1]);
	      }
	    }
	  }
	  /*multiple stem loop contribution*/
	  ii = iindx[i+1]; /* ii-k=[i+1,k-1] */
	  temp = 0.0;
	  for (k=i+2; k<=j-1; k++) temp += qm[i+1][k-1]*qqm1[k];
	  tt = rtype[type];
	  qbt1 += temp*expMLclosing*expMLintern[tt]*scale[2]*
	    expdangle3[tt][S1[i+1]]*expdangle5[tt][S1[j-1]];

	  qb[i][j] = qbt1;
	} /* end if (type!=0) */
	else qb[i][j] = 0.0;

	/* construction of qqm matrix containing final stem
	   contributions to multiple loop partition function
	   from segment i,j */
	qqm[i] = qqm1[i]*expMLbase[1];
	if (type) {
	  qbt1 = qb[i][j]*expMLintern[type];
	  if (i>1) qbt1 *= expdangle5[type][S1[i-1]];/*wann dangle??*/
	  if (j<n) qbt1 *= expdangle3[type][S1[j+1]];/*wann dangle??*/
	  else if (type>2) qbt1 *= expTermAU;
	  qqm[i] += qbt1;
	}
	if (qm1) qm1[jindx[j]+i] = qqm[i]; /* for stochastic backtracking */

	/*construction of qm matrix containing multiple loop
	  partition function contributions from segment i,j */
	temp = 0.0;
	/*ii = iindx[i];   ii-k=[i,k-1] */
	for (k=i+1; k<=j; k++) temp += (qm[i][k-1]+expMLbase[k-i])*qqm[k];
	qm[i][j] = (temp + qqm[i]);

	/*auxiliary matrix qq for cubic order q calculation below */
	qbt1 = qb[i][j];
	if (type) {
	  if (i>1/*MAX(1,j-winSize+1)*/) qbt1 *= expdangle5[type][S1[i-1]];
	  if (j<n/*winSize*/) qbt1 *= expdangle3[type][S1[j+1]];
	  else if (type>2) qbt1 *= expTermAU;
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
	  sprintf(msg, "overflow in pf_fold while calculating q[%d,%d]\n"
		  "use larger pf_scale", i,j);
	  nrerror(msg);
	}
      } /*end for i??*/
      tmp = qq1;  qq1 =qq;  qq =tmp;
      tmp = qqm1; qqm1=qqm; qqm=tmp;
    }

#if 0  /* no output in library functions! */

    if ((j>winSize)&&(j<=n)) {
      /*Ausgabe von freenergy fenster*/
      int middle=j-winSize/2;
      double Fwind;
      Fwind=(-log(q[j-winSize+1][j])-winSize*log(pf_scale))*(temperature+K0)*GASCONST/1000.0;
       printf("%6d\t%.7f\n",middle, Fwind);
    }
#endif
    if (j>winSize) {
      Qmax=0;
      /*i=j-winSize;*/
      /*initialize multiloopfs*/
      for (k=j-winSize+1; k<=MIN(n,j); k++) {
	prml[k]=0;
	prm_l[k]=0;
	/*	prm_l1[k]=0;  others stay??*/
      }
      k=j-winSize;
      for (l=k+TURN+1; l<=MIN(n,k+winSize-1); l++) {
	int a;
	pR[k][l] = 0; /*set zero at start*/
	type=ptype[k][l];
	if (qb[k][l]==0) continue;

	for (a=MAX(1,l-winSize+2); a<MIN(k,n-winSize+2);a++)
	  pR[k][l]+=q[a][k-1]*q[l+1][a+winSize-1]/q[a][a+winSize-1];

	if (l-k+1==winSize)
	  pR[k][l]+=1./q[k][l];
	else {
	  if (k+winSize-1<=n)          /*k outermost*/
	    pR[k][l]+=q[l+1][k+winSize-1]/q[k][k+winSize-1];
	  if (l-winSize+1>=1)  /*l outermost*/
	    pR[k][l]+=q[l-winSize+1][k-1]/q[l-winSize+1][l];
	}
	if (k>1)
	  pR[k][l]*=expdangle5[type][S1[k-1]];
	if (l<n)
	  pR[k][l]*= expdangle3[type][S1[l+1]];
	else if (type>2)
	  pR[k][l] *= expTermAU;


	/*initialize multiloopfs*/
	prml[j-winSize]=0;
	prm_l[j-winSize]=0;
	prm_l1[j-winSize]=0;

	/* pr = q; */    /* recycling */

	type_2 = ptype[k][l]; type_2 = rtype[type_2];

	for (i=MAX(MAX(l-winSize+1,k-MAXLOOP-1),1); i<=k-1; i++) {
	  for (m=l+1; m<=MIN(MIN(l+ MAXLOOP -k+i+2,i+winSize-1),n); m++) {
	    type = ptype[i][m];
	    if ((pR[i][m]>0))
	      pR[k][l] += pR[i][m]*expLoopEnergy(k-i-1, m-l-1, type, type_2,
						 S1[i+1], S1[m-1], S1[k-1], S1[l+1]);
	  }
	}
      }
      /* 3. bonding k,l as substem of multi-loop enclosed by i,m */
      prm_MLb = 0.;
      if(k>1) /*sonst nix!*/
	for (l=MIN(n-1,k+winSize-2/*??*/); l>=k+TURN+1; l--) { /*opposite direction*/
	  m=l+1;
	  prmt = prmt1 = 0.0;
	  tt = ptype[k-1][m]; tt=rtype[tt];
	  prmt1 = pR[k-1][m]*expMLclosing*expMLintern[tt]*
	    expdangle3[tt][S1[k]]*expdangle5[tt][S1[l]];
	  for (i=MAX(1,l-winSize+2/*?*/); i<k-1/*TURN*/; i++) {
	    tt = ptype[i][m]; tt = rtype[tt];
	    prmt += pR[i][m]*expdangle3[tt][S1[i+1]]*
	      expdangle5[tt][S1[m-1]] *qm[i+1][k-1];
	  }
	  tt = ptype[k][l];
	  prmt *= expMLclosing*expMLintern[tt];
	  prml[ m] = prmt;
	  prm_l[m] = prm_l1[m]*expMLbase[1]+prmt1;

	  prm_MLb = prm_MLb*expMLbase[1] + prml[m];
	  /* same as:    prm_MLb = 0;
	     for (i=1; i<=k-1; i++) prm_MLb += prml[i]*expMLbase[k-i-1]; */

	  prml[m] = prml[ m] + prm_l[m];

	  if (qb[k][l] == 0.) continue;

	  temp = prm_MLb;

	  for (m=MIN(k+winSize-2,n)/*1??*/;m>=l+2; m--)
	    temp += prml[m]*qm[l+1][m-1];

	  temp *= expMLintern[tt]*scale[2];
	  if (k>1) temp *= expdangle5[tt][S1[k-1]];
	  if (l<n) temp *= expdangle3[tt][S1[l+1]];
	  else temp *= expTermAU;
	  pR[k][l] += temp;

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

	} /* end for (k=..) */
      tmp = prm_l1; prm_l1=prm_l; prm_l=tmp;

      /* end for (l=..)   */

      if (j-2*winSize>0) {
	printpbar(pR,winSize,j-2*winSize,n,cutoff);
	*pl=get_plistW(*pl, n, cutoff, j-2*winSize, pR, winSize);
	FreeOldArrays(j-2*winSize);
      }

    }   /* end if (do_backtrack)*/


  }/*end for j */

  /*finish output and free*/
  for (j=n-winSize+1; j<=n; j++) {
    printpbar(pR,winSize,j,n,cutoff);
    *pl=get_plistW(*pl, n, cutoff, j, pR, winSize);
    FreeOldArrays(j);
  }
  if (ov>0) fprintf(stderr, "%d overflows occurred while backtracking;\n"
		    "you might try a smaller pf_scale than %g\n",
		    ov, pf_scale);

  return 1;
}

/*------------------------------------------------------------------------*/
/* dangling ends should never be destabilizing, i.e. expdangle>=1         */
/* specific heat needs smooth function (2nd derivative)                   */
/* we use a*(sin(x+b)+1)^2, with a=2/(3*sqrt(3)), b=Pi/6-sqrt(3)/2,       */
/* in the interval b<x<sqrt(3)/2                                          */

#define SCALE 10
#define SMOOTH(X) ((X)/SCALE<-1.2283697)?0:(((X)/SCALE>0.8660254)?(X):	\
					    SCALE*0.38490018*(sin((X)/SCALE-0.34242663)+1)*(sin((X)/SCALE-0.34242663)+1))
/*#define SMOOTH(X) ((X)<0 ? 0 : (X)) */

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
  }
  if (pf_scale<1) pf_scale=1;
  scale[0] = 1.;
  for (i=1; i<=length; i++)
    scale[i] = scale[i-1]/pf_scale;

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
{/*arrays in 2 dimensions*/

  q   = (FLT_OR_DBL **) space((length+1)*sizeof(FLT_OR_DBL *));
  qb  = (FLT_OR_DBL **) space((length+1)*sizeof(FLT_OR_DBL *));
  qm  = (FLT_OR_DBL **) space((length+1)*sizeof(FLT_OR_DBL *));
  pR = (FLT_OR_DBL **) space((length+1)*sizeof(FLT_OR_DBL *));
  ptype = (char **) space((length+1)*sizeof(char *));
  /*  if (st_back) {
      qm1 = (FLT_OR_DBL *) space(size);
      }*/
  /*length??*/
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
  /*  for (i=1; i<=length; i++) {
      iindx[i] = ((length+1-i)*(length-i))/2 +length+1;
      jindx[i] = (i*(i-1))/2;
      }*/
}

/*----------------------------------------------------------------------*/

PUBLIC void init_pf_foldLP(int length)
{
  if (length<1) nrerror("init_pf_fold: length must be greater 0");
  if (init_length>0) free_pf_arraysLP(); /* free previous allocation */
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

PUBLIC void free_pf_arraysLP(void)
{
  free(q);
  free(qb);
  free(qm);
  free(pR);
  q=pR=NULL;
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

PUBLIC void update_pf_paramsLP(int length)
{
  if (length>init_length) init_pf_foldLP(length);  /* init not update */
  else {
    /*   make_pair_matrix();*/
    scale_pf_params((unsigned) length);
  }
}

/*---------------------------------------------------------------------------*/


PRIVATE void printpbar(FLT_OR_DBL **prb,int winSize, int i, int n, float cutoff) {
  int j;
  int howoften=0; /* how many samples do we have for this pair */
  int pairdist;

  for (j=i+TURN; j<MIN(i+winSize,n+1); j++) {
    pairdist=(j-i+1);
    /*4cases*/
    howoften=MIN(winSize-pairdist+1,i); /*pairdist,start*/
    howoften=MIN(howoften,n-j+1);       /*end*/
    howoften=MIN(howoften,n-winSize+1); /*windowsize*/
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
  ptype[j]=(char *)space((winSize+1)*sizeof(char));
  ptype[j]-=j;
  return;
}


PRIVATE void GetPtype(int i, int winSize,const short *S,int n) {
  /*make new entries in ptype array*/
  int j;
  int type;
  for (j=i; j<=MIN(i+winSize,n); j++) {
    type = pair[S[i]][S[j]];
    /*qb[i][j]=0;  what for?*/
    ptype[i][j] = (char) type;
  }
  return;
}

PRIVATE struct plist *get_plistW(struct plist *pl, int length, double cutoff,
				 int start, FLT_OR_DBL **Tpr, int winSize) {
  /*get pair probibilities out of pr array*/
  int  j,  max_p;
  int n=2;

  max_p=1000;
  while (max_p<num_p)
    max_p*=2;

  for (j=start+1; j<=MIN(start+winSize, length); j++) {
    if (Tpr[start][j]<cutoff) continue;
    if (num_p==max_p-1) {
      max_p*=2;
      pl=(struct plist *)xrealloc(pl,max_p*sizeof(struct plist));
    }
    pl[num_p].i=start;
    pl[num_p].j=j;
    pl[num_p++].p=Tpr[start][j];
  }

  /* mark end of data with zeroes */
  pl[num_p].i=0;
  pl[num_p].j=0;
  pl[num_p].p=0.;
  /*  pl=(struct plist *)xrealloc(pl,(count)*sizeof(struct plist));*/
  return pl;
}
