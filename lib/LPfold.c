/* Last changed Time-stamp: <2009-10-12 13:24:15 ivo> */
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
#include "LPfold.h"

/*@unused@*/
static char rcsid[] UNUSED = "$Id: LPfold.c,v 1.8 2009/02/18 20:34:38 ivo Exp $";

#define MAX(x,y) (((x)>(y)) ? (x) : (y))
#define MIN(x,y) (((x)<(y)) ? (x) : (y))
#define PUBLIC
#define PRIVATE static

static float cutoff;
static int num_p=0; /*for counting basepairs*/
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
PRIVATE void printpbar(FLT_OR_DBL **prb,int winSize, int i, int n);
PRIVATE  void  init_pf_foldLP(int length);
PRIVATE  void  free_pf_arraysLP(void);

PRIVATE struct plist *get_deppp(struct plist *pl, int start, int pairsize, int length);
PRIVATE struct plist *get_plistW(struct plist *pl, int length, int start, FLT_OR_DBL **Tpr, int winSize);
PRIVATE void print_plist(int length, int start, FLT_OR_DBL **Tpr, int winSize, FILE *fp);
PRIVATE void compute_pU(int k, int ulength, double **pU, int winSize, int n, char *sequence);
PRIVATE void putoutpU(double **pU,int k, int ulength, FILE *fp);
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
PRIVATE FLT_OR_DBL **q, **qb, **qm, *qm1, *qqm, *qqm1, *qq, *qq1, **pR, **qm2, **QI5,  **q2l, **qmb;/*,**QI3,*/
PRIVATE FLT_OR_DBL *prml, *prm_l, *prm_l1, *q1k, *qln;
PRIVATE FLT_OR_DBL *scale;
PRIVATE char **ptype; /* precomputed array of pair types */
PRIVATE int *jindx;
PRIVATE int init_length;  /* length in last call to init_pf_fold() */
PRIVATE double init_temp; /* temperature in last call to scale_pf_params */
#define ISOLATED  256.0

/*-----------------------------------------------------------------*/
static  short *S, *S1;
static int unpaired;
static int ulength;
static int pUoutput;
PUBLIC struct plist *pfl_fold(char *sequence, int winSize, int pairSize, float cutoffb, double **pU, struct plist **dpp2, FILE *pUfp, FILE *spup)
{

  int n, m,i,j,k,l, u,u1,ii, type, type_2, tt, ov=0;
  FLT_OR_DBL temp, Qmax=0, prm_MLb;
  FLT_OR_DBL prmt,prmt1;
  FLT_OR_DBL qbt1, *tmp;
  double max_real;
  int do_dpp=0;
   int simply_putout=0;
  struct plist *dpp;
  struct plist *pl;
  pUoutput=0;
  ulength=0;
  cutoff=cutoffb;
  if (pU != NULL) ulength=(int) pU[0][0]+0.49;
  if (spup !=NULL) simply_putout=1; /*can't have one without the other*/
  if (pUfp!=NULL) pUoutput=1;
  else if ((pUoutput)&&(ulength!=0)) {
    fprintf(stderr, "There was a problem with non existing File Pointer for unpaireds, terminating process\n");
    return pl;
  }

  dpp=*dpp2;
  if (dpp !=NULL) do_dpp=1;


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

  /*here, I allocate memory for pU, if has to be saved, I allocate all in one go,
    if pU is put out and freed, I only allocate what I really need*/

  if (ulength>0){
    if (pUoutput) {
      for (i=1; i<=ulength; i++) pU[i]=(double *)space((MAX(MAXLOOP,ulength)+2)*sizeof(double));
    }
    else {
      for (i=1; i<=n; i++) pU[i]=(double *)space((MAX(MAXLOOP,ulength)+2)*sizeof(double));
     }
  }

  /*array initialization ; qb,qm,q
    qb,qm,q (i,j) are stored as ((n+1-i)*(n-i) div 2 + n+1-j */
  num_p=0;
  pl=(struct plist *)space(1000*sizeof(struct plist));


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
	  if (i>1) qbt1 *= expdangle5[type][S1[i-1]];
	  if (j<n) qbt1 *= expdangle3[type][S1[j+1]];
	  else if (type>2) qbt1 *= expTermAU;
	  qqm[i] += qbt1;
	}
	if (qm1) qm1[jindx[j]+i] = qqm[i]; /* for stochastic backtracking */

	/*construction of qm matrix containing multiple loop
	  partition function contributions from segment i,j */
	temp = 0.0;
	/*ii = iindx[i];   ii-k=[i,k-1] */
	/*new qm2 computation done here*/
	for (k=i+1; k<=j; k++) temp += (qm[i][k-1])*qqm[k];
	if (ulength>0) qm2[i][j]=temp;/*new qm2 computation done here*/
	for (k=i+1; k<=j; k++) temp += expMLbase[k-i]*qqm[k];
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
	  snprintf(msg, 128, "overflow in pf_fold while calculating q[%d,%d]\n"
		  "use larger pf_scale", i,j);
	  nrerror(msg);
	}
      } /*end for i*/
      tmp = qq1;  qq1 =qq;  qq =tmp;
      tmp = qqm1; qqm1=qqm; qqm=tmp;
    }

    /*just as a general service, I save here the free energy of the windows
      no output is generated, however,...     */
    if ((j>=winSize)&&(j<=n)&&(ulength)&&!(pUoutput)) {
      double Fwindow=0.;
      Fwindow=(-log(q[j-winSize+1][j])-winSize*log(pf_scale))*(temperature+K0)*GASCONST/1000.0;

      pU[j][0]=Fwindow;
      /* if (ulength>=winSize)
	 pU[j][winSize]=scale[winSize]/q[j-winSize+1][j];
      */
    }
    if (j>winSize) {
      Qmax=0;
      /* i=j-winSize; */
      /* initialize multiloopfs */
      for (k=j-winSize; k<=MIN(n,j); k++) {
	prml[k]=0;
	prm_l[k]=0;
	/*	prm_l1[k]=0;  others stay*/
      }
      prm_l1[j-winSize]=0;
      k=j-winSize;
      for (l=k+TURN+1; l<=MIN(n,k+winSize-1); l++) {
	int a;
	pR[k][l] = 0; /* set zero at start */
	type=ptype[k][l];
	if (qb[k][l]==0) continue;

	for (a=MAX(1,l-winSize+2); a<MIN(k,n-winSize+2);a++)
	  pR[k][l]+=q[a][k-1]*q[l+1][a+winSize-1]/q[a][a+winSize-1];

	if (l-k+1==winSize)
	  pR[k][l]+=1./q[k][l];
	else {
	  if (k+winSize-1<=n)          /* k outermost */
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
	type_2 = ptype[k][l]; type_2 = rtype[type_2];

	for (i=MAX(MAX(l-winSize+1,k-MAXLOOP-1),1); i<=k-1; i++) {
	  for (m=l+1; m<=MIN(MIN(l+ MAXLOOP -k+i+2,i+winSize-1),n); m++) {
	    type = ptype[i][m];
	    if ((pR[i][m]>0))
	      pR[k][l] += pR[i][m]*expLoopEnergy(k-i-1, m-l-1, type, type_2,
						 S1[i+1], S1[m-1], S1[k-1], S1[l+1]);
	  }
	}
	if (ulength) { /* NOT IF WITHIN INNER LOOP */
	  for (i=MAX(MAX(l-winSize+1,k-MAXLOOP-1),1); i<=k-1; i++) {
	    for (m=l+1; m<=MIN(MIN(l+ MAXLOOP -k+i+2,i+winSize-1),n); m++) {
	      type = ptype[i][m];
	      if ((pR[i][m]>0)){
		temp=pR[i][m]*qb[k][l]*expLoopEnergy(k-i-1, m-l-1, type, type_2,
						     S1[i+1], S1[m-1], S1[k-1], S1[l+1]);
		QI5[l][m-l-1]+=temp;
		QI5[i][k-i-1]+=temp;
	      }
	    }
	   }
	}
      }
      /* 3. bonding k,l as substem of multi-loop enclosed by i,m */
      prm_MLb = 0.;
      if (k>1) /*sonst nix!*/
	for (l=MIN(n-1,k+winSize-2); l>=k+TURN+1; l--) { /* opposite direction */
	  m=l+1;
	  prmt = prmt1 = 0.0;
	  tt = ptype[k-1][m]; tt=rtype[tt];
	  prmt1 = pR[k-1][m]*expMLclosing*expMLintern[tt]*
	    expdangle3[tt][S1[k]]*expdangle5[tt][S1[l]];
	  for (i=MAX(1,l-winSize+2); i<k-1/*TURN*/; i++) {
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
	     for (i=n; i>k; i--)  prm_MLb += prml[i]*expMLbase[k-i-1];
	  */
	  prml[m] = prml[ m] + prm_l[m];

	  if (qb[k][l] == 0.) continue;

	  temp = prm_MLb;

	  if (ulength) {
	    double dang;
	    /* coefficient for computations of unpairedarrays */
	    dang =  expMLintern[tt]*scale[2];
	    dang *= expdangle5[tt][S1[k-1]];
	    dang *= expdangle3[tt][S1[l+1]];
	    dang*=qb[k][l];
	    for (m=MIN(k+winSize-2,n);m>=l+2; m--){
	      qmb[l][m-l-1]+=prml[m]*dang;
	      q2l[l][m-l-1]+=(prml[m]-prm_l[m])*dang;
	    }
	  }

	  for (m=MIN(k+winSize-2,n);m>=l+2; m--)
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

	} /* end for (l=..) */
      tmp = prm_l1; prm_l1=prm_l; prm_l=tmp;

      /* end for (l=..)   */
      if ((ulength)&&(k-MAXLOOP-1>0)){
	compute_pU(k-MAXLOOP-1,ulength,pU, winSize, n, sequence);
	
	/* here, we print and free pUs not in use any more (hopefully)*/
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

  }/*end for j */

  /*finish output and free*/
  for (j=MAX(1,n-MAXLOOP); j<=n;j++) {
    /*   if (pUoutput) pU[j]=(double *)space((ulength+2)*sizeof(double));*/
    if (ulength) compute_pU(j,ulength,pU, winSize, n, sequence);
    /*here, we put out and free pUs not in use any more (hopefully)*/
    if (pUoutput) putoutpU(pU,j, ulength, pUfp);
  }
  for (j=MAX(n-winSize-MAXLOOP,1); j<=n; j++) {
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
  free_pf_arraysLP();
  if (ov>0) fprintf(stderr, "%d overflows occurred while backtracking;\n"
		    "you might try a smaller pf_scale than %g\n",
		    ov, pf_scale);
  *dpp2=dpp;
  
  return pl;
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
  scale[0] = 1.; scale[1] = 1./pf_scale;
  for (i=2; i<=length; i++)
    scale[i] = scale[i/2]*scale[i-(i/2)];

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
  ptype = (char **) space((length+2)*sizeof(char *));
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

  if (ulength>0) {
  QI5 = (FLT_OR_DBL **) space((length+1)*sizeof(FLT_OR_DBL *));
  /* QI3 = (FLT_OR_DBL **) space((length+1)*sizeof(FLT_OR_DBL *));*/
  qmb = (FLT_OR_DBL **) space((length+1)*sizeof(FLT_OR_DBL *));
  qm2 =(FLT_OR_DBL **) space((length+1)*sizeof(FLT_OR_DBL *));
  q2l =(FLT_OR_DBL **) space((length+1)*sizeof(FLT_OR_DBL *));
  }
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
  if (unpaired) free(qm2);
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
  if (ulength!=0) {
    free(QI5);
    free(qmb);
    free(qm2);
    free(q2l);
  }

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


PRIVATE void printpbar(FLT_OR_DBL **prb,int winSize, int i, int n) {
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
  for (j=i; j<=MIN(i+winSize,n); j++) {
    type = pair[S[i]][S[j]];
    ptype[i][j] = (char) type;
  }
  return;
}

PRIVATE struct plist *get_plistW(struct plist *pl, int length,
				 int start, FLT_OR_DBL **Tpr, int winSize) {
  /*get pair probibilities out of pr array*/
  int  j,  max_p;
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
PRIVATE struct plist *get_deppp(struct plist *pl, int start, int pairsize, int length) {
  /*compute dependent pair probabilities*/
  int i, j, count=0;
  double tmp;
  struct plist *temp;
  temp=(plist *)space(pairsize*sizeof(plist)); /*holds temporary deppp*/
  for (j=start+TURN; j<MIN(start+pairsize,length); j++) {

    if ((qb[start][j]*qb[start-1][(j+1)])>10e-200) {
      int type=ptype[start-1][j+1];
      int type_2=rtype[ptype[start][j]];
      tmp=qb[start][j]/qb[start-1][(j+1)]*expLoopEnergy(0, 0, type, type_2,
							S1[start], S1[j], S1[start-1], S1[j+1]);
       temp[count].i=start;
      temp[count].j=j;
      temp[count++].p=tmp;
    }
  }
  /*write it to list of deppps*/
  for (i=0; pl[i].i!=0; i++);
  pl=(struct plist *)xrealloc(pl,(i+count+1)*sizeof(struct plist));
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
  /*print out of pr array, do not save*/
  int  j;


  for (j=start+1; j<=MIN(start+winSize, length); j++) {
    if (Tpr[start][j]<cutoff) continue;
    fprintf(fp,"%d  %d  %g\n",start,j,Tpr[start][j]);
  }

  /* mark end of data with zeroes */

  return ;
}

PRIVATE void compute_pU(int k, int ulength, double **pU, int winSize,int n, char *sequence) {
/*here, we try to add a function computing all unpaired probabilities starting at some i, going down to $unpaired, to be unpaired, i.e. a list with entries from 1 to unpaired for every i, with the probability of a stretch of length x, starting at i-x+1, to be unpaired*/
  int startu;
  int i5;
  int j3, len, obp;
  double temp;
  double *QBE;

  QBE=(double *) space((MAX(ulength,MAXLOOP)+2)*sizeof(double));

  /*first, we will*/
 /*for k<=ulength, pU[k][k]=0, because no bp can enclose it*/
  if (pUoutput&&k+ulength<=n)  pU[k+ulength]=(double *)space((ulength+2)*sizeof(double));
  /*compute pu[k+ulength][ulength] */
   for (i5=MAX(k+ulength-winSize+1,1);i5<=k;i5++) {
    for (j3=k+ulength+1; j3<=MIN(n,i5+winSize-1); j3++) {
      /*  if (k>400) {
	printf("i%d j%d  ",i5,j3);
	fflush(stdout);
	}*/
      if (ptype[i5][j3]!=0) {/**/
	/*(.. >-----|..........)
	  i5  j     j+ulength  j3	      */
	/*Multiloops*/
	if (i5<k)temp=qm2[i5+1][k]*expMLbase[j3-k-1];/*(..{}{}-----|......)*/
	else temp=0.;
	if (j3-1>k+ulength)temp+=qm2[k+ulength+1][j3-1]*expMLbase[(k+ulength-i5)];/*(..|-----|{}{})*/
	if ((i5<k)&&(j3-1>k+ulength))temp+=qm[i5+1][k]*qm[k+ulength+1][j3-1]*expMLbase[ulength];/*({}|-----|{})*/
	/* add dangles, multloopclosing etc. */
	temp*=scale[2];
	temp*=expMLclosing*expMLintern[rtype[ptype[i5][j3]]];
	temp*=expdangle3[rtype[ptype[i5][j3]]][S1[i5+1]]*expdangle5[rtype[ptype[i5][j3]]][S1[j3-1]];
	/*add hairpins*/
	temp+=expHairpinEnergy(j3-i5-1, ptype[i5][j3], S1[i5+1], S1[j3-1], sequence+i5-1);
	/*add outer probability*/
	temp*=pR[i5][j3];
	pU[k+ulength][ulength]+=temp;

      }
    }
   }
   /*code doubling to avoid if within loop*/
#if 0
  /*initialization for interior loops,
    it is not recomended to have verysmall ulengths!!*/
  if (ulength<MAXLOOP) {
    int k5;
    int l3;
    int outype;
    /*kl bp is 5'*/
    /*MAXLOOP>((l5-k5-1)+(j3-l3-1)
      k-winSize+ulength<i5<k-TURN-1;
      k+ulength<j3<=k+MAXLOOP+1
      if i then use l3, it is easier by far:
      j3-MAXLOOP<=l3<=k
      i5<k5<k-TURN k5<=i5+l3+2+MAXLOOP-j3
      k5+TURN<l3<=k
    */
    for (i5=MAX(k+ulength-winSize,1);i5<k-TURN-1;i5++) {

      for (j3=k+ulength+1; j3<=MIN(n,MIN(i5+winSize-1,k+MAXLOOP+1)); j3++) {
	double temp=0;
	if (outype=ptype[i5][j3]>0) /*oder so halt*/
	  for (l3=MAX(i5+TURN+1,j3-MAXLOOP-1); l3<=k; l3++){
	    for (k5=i5+1; k5<=MIN(l3-TURN-1,MAXLOOP+i5+l3+2-j3); k5++){
	      if (ptype[k5][l3]) {
		temp+= qb[k5][l3]*expLoopEnergy(k5-i5-1, j3-l3-1, outype, rtype[ptype[k5][l3]], S1[i5+1], S1[j3-1], S1[k5-1], S1[l3+1]);
	      }
	    }
	  }
	temp*=pR[i5][j3];
	pU[k+ulength][ulength]+= temp;
      }
    }
    /*kl bp is 3'*/
    /*
      k+ulength-MAXLOOP<=i5<=k
      k+ulength+1+TURN<j3<i5+winSize
      k+ulength+1<=k5<i5+MAXLOOP+2 || k5<j3-TURN
      k5<l3<j3 || j3-k5-i5-2-ML<=l3<j3
    */
    for (i5=MAX(1,MAX(k+ulength-winSize,k+ulength-MAXLOOP));i5<=k; i5++){
      for (j3=k+ulength+TURN+2; j3<MIN(n+1,i5+winSize); j3++) {
	double temp = 0;
	if (outype=ptype[i5][j3]>0) /*oder so halt*/
	  for (k5=k+ulength+1; k5<MIN(j3-TURN-1,i5+MAXLOOP+2); k5++) {
	    for (l3=MAX(k5+TURN+1,j3+k5-i5-2-MAXLOOP); l3<j3; l3++) {
	      if (ptype[k5][l3])
		temp += qb[k5][l3]*expLoopEnergy(k5-i5-1, j3-l3-1, outype, rtype[ptype[k5][l3]], S1[i5+1], S1[j3-1], S1[k5-1], S1[l3+1]);
	    }
	  }
	temp*=pR[i5][j3];
	pU[k+ulength][ulength]+= temp;
      }
    }
  }
  /*Add up Is QI5[l][m-l-1] QI3*/
  /*Add up Interior loop terms*/
  temp=0.;

  for (len=winSize; len>=ulength; len--) temp+=QI3[k][len];
  for ( ;len>0; len--) {
    temp+=QI3[k][len];
    QBE[len]+=temp;
  }
#endif
  temp=0.;
  for (len=winSize; len>=MAX(ulength,MAXLOOP); len--) temp+=QI5[k][len];
  for ( ;len>0; len--) { /*grenzen?*/
    temp+=QI5[k][len];
    QBE[len]+=temp;  /*replace QBE with QI*/
  }
  /*Add Hairpinenergy to QBE*/
  temp=0.;
  for (obp=MIN(n,k+winSize-1);obp>k+ulength; obp--)  if (ptype[k][obp])  temp+=pR[k][obp]*expHairpinEnergy(obp-k-1, ptype[k][obp], S1[k+1], S1[obp-1], sequence+k-1);
  for (obp=MIN(n,MIN(k+winSize-1,k+ulength)); obp>k+1; obp--) {
    if (ptype[k][obp])  temp+=pR[k][obp]*expHairpinEnergy(obp-k-1, ptype[k][obp], S1[k+1], S1[obp-1], sequence+k-1);
    QBE[obp-k-1]+=temp;  /*add hairpins to QBE (all in one array)*/
  }
  /*doubling the code to get the if out of the loop*/

  /*Add up Multiloopterms  qmb[l][m]+=prml[m]*dang;
    q2l[l][m]+=(prml[m]-prm_l[m])*dang; */

  temp=0.;
  for (len=winSize;len>=ulength; len--)temp+=q2l[k][len]*expMLbase[len];
  for ( ;len>0; len--) {
    temp+=q2l[k][len]*expMLbase[len];
    QBE[len]+=temp;/*add (()()____) type cont. to I3*/
  }
  for (len=1; len<ulength;len++) {
    for (obp=k+len+TURN; obp<=MIN(n,k+winSize-1); obp++) {
      /*add (()___())*/
      QBE[len]+=qmb[k][obp-k-1]*qm[k+len+1/*2*/][obp-1]*expMLbase[len];
    }
  }
  for (len=1; len<ulength; len++) {
    for (obp=k+len+TURN+TURN; obp<=MIN(n,k+winSize-1); obp++) {
      if (ptype[k][obp]) {
	temp=scale[2]*expMLbase[len]; /*k:obp*/
	/**/
	temp*=expMLclosing*expMLintern[rtype[ptype[k][obp]]];
	temp*=expdangle3[rtype[ptype[k][obp]]][S1[k+1]]*expdangle5[rtype[ptype[k][obp]]][S1[obp-1]];
	QBE[len]+=pR[k][obp]*temp*qm2[k+len+1][obp-1]; /*add (___()())*/
      }
    }
  }
  /*After computing all these contributions in QBE[len], that k is paired
    and the unpaired stretch is AT LEAST len long, we start to add that to
    the old unpaired thingies;*/
  for (len=1; len<MIN(MAX(ulength,MAXLOOP),n-k); len++) {
    pU[k+len][len]+=pU[k+len][len+1]+QBE[len];
  }

  /*open chain*/
  if ((ulength>=winSize)&&(k>=ulength)) {
    pU[k][winSize]=scale[winSize]/q[k-winSize+1][k];
  }
  /* now the not enclosed by any base pair terms for whatever it is we do not need anymore...
    ... which should be e.g; k, again */
  for (startu=MIN(ulength,k); startu>0; startu--) {
    temp=0.;
    for (i5=MAX(1,k-winSize+2); i5<=MIN(k-startu,n-winSize+1); i5++) {
      temp+=q[i5][k-startu]*q[k+1][i5+winSize-1]*scale[startu]/q[i5][i5+winSize-1];
    }
    /* the 2 Cases where the borders are on the edge of the interval */
    if((k>=winSize)&&(startu+1<=winSize)) temp+=q[k-winSize+1][k-startu]*scale[startu]/q[k-winSize+1][k];
    if((k<=n-winSize+startu)&&(k-startu>=0)&&(k<n)&&(startu+1<=winSize)) temp+=q[k+1][k-startu+winSize]*scale[startu]/q[k-startu+1][k-startu+winSize];
    /*Divide by number of possible windows*/
    pU[k][startu]+=temp;
      {
      int leftmost, rightmost;

      leftmost = MAX(1,k-winSize+1);
      rightmost = MIN(n-winSize+1,k-startu+1);
      pU[k][startu]/=(rightmost-leftmost+1);
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
  for (i=1; i<=MIN(ulength,k); i++) {
    fprintf(fp,"%.5g\t",pUx[k][i]);
  }
  fprintf(fp,"\n");
  free(pUx[k]);
}

PUBLIC void putoutpU_prob(double **pU,int length, int ulength, FILE *fp, int energies) {
  /*put out unpaireds */
  int i,k;
  double kT= (temperature+K0)*GASCONST/1000.0;
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
  free(pU[0]);
  free(pU);
  fflush(fp);
}


/*
 Here: Space for questions...
*/
