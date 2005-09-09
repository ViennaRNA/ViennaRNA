/* Last changed Time-stamp: <2005-08-31 15:32:49 ivo> */
/*                
		  partiton function and base pair probabilities
		  for RNA secvondary structures 
		  of a set of aligned sequences

		  Ivo L Hofacker
		  Vienna RNA package
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>    /* #defines FLT_MIN */
#include "utils.h"
#include "energy_par.h"
#include "fold_vars.h"
#include "pair_mat.h"
#include "alifold.h"
/*@unused@*/
static char rcsid[] = "$Id: alipfold.c,v 1.10 2005/09/09 08:03:16 ivo Exp $";

#define MAX(x,y) (((x)>(y)) ? (x) : (y))
#define MIN(x,y) (((x)<(y)) ? (x) : (y))
#define PUBLIC
#define PRIVATE static
#define STACK_BULGE1  1   /* stacking energies for bulges of size 1 */
#define NEW_NINIO     1   /* new asymetry penalty */



PUBLIC  float alipf_fold(char **sequences, char *structure, pair_info **pi);
PRIVATE void  init_alipf_fold(int length, int n_seq);
PRIVATE void  free_alipf_arrays(void);
/* PRIVATE void  update_alipf_params(int length); */
PRIVATE void  sprintf_bppm(int length, char *structure);
PRIVATE void  scale_pf_params(unsigned int length, int n_seq);
PRIVATE void  get_arrays(unsigned int length);
PRIVATE double expLoopEnergy(int u1, int u2, int type, int type2,
			     short si1, short sj1, short sp1, short sq1);
PRIVATE void make_pscores(const short *const *S, const char *const* AS,
			  int n_seq, const char *structure);
PRIVATE pair_info *make_pairinfo(const short *const* S, char **AS, 
				 int n_seq);
PRIVATE short * encode_seq(const char *sequence);
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
PRIVATE FLT_OR_DBL *q, *qb, *qm, *qqm, *qqm1, *qq, *qq1;
PRIVATE FLT_OR_DBL *prml, *prm_l, *prm_l1, *q1k, *qln;
PRIVATE FLT_OR_DBL *scale;
PRIVATE short *pscore;   /* precomputed array of covariance bonus/malus */ 
PRIVATE int init_length; /* length in last call to init_pf_fold() */
#define ISOLATED  256.0

#define UNIT 100
#define MINPSCORE -2 * UNIT

extern double cv_fact /* =1 */;
extern double nc_fact /* =1 */;

/*-----------------------------------------------------------------*/
PUBLIC float alipf_fold(char **sequences, char *structure, pair_info **pi)
{
  short **S;
  int s, *type;
  int n, n_seq, i,j,k,l, ij, kl, u,u1,d,ii,ll, type_2, tt, ov=0;
  FLT_OR_DBL temp, Q, Qmax=0, prm_MLb;
  FLT_OR_DBL prmt,prmt1;
  FLT_OR_DBL qbt1, *tmp;
   
  double free_energy, kTn;

  n = (int) strlen(sequences[0]);
  for (s=0; sequences[s]!=NULL; s++); 
  n_seq = s;
  init_alipf_fold(n, n_seq);  /* (re)allocate space */
  kTn = (temperature+K0)*GASCONST*n_seq/10.;   /* kT in cal/mol  */
  
  S = (short **) space(sizeof(short *)*(n_seq+1));
  type = (int *) space(n_seq*sizeof(int));
  for (s=0; s<n_seq; s++) { 
    if (strlen(sequences[s]) != n) nrerror("uneqal seqence lengths");
    S[s] = encode_seq(sequences[s]);
  }
  make_pscores((const short *const*)S, sequences, n_seq, structure);
   
  /* array initialization ; qb,qm,q
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
      int ij, psc;
      /* construction of partition function for segment i,j */
      /* calculate pf given that i and j pair: qb(i,j)      */
      u = j-i-1; ij = iindx[i]-j;

      for (s=0; s<n_seq; s++) {
	type[s] = pair[S[s][i]][S[s][j]];
	if (type[s]==0) type[s]=7;
      }
      psc = pscore[ij];
      if (psc>=cv_fact*MINPSCORE) {   /* otherwise ignore this pair */
	
	/* hairpin contribution */
	for (qbt1=1,s=0; s<n_seq; s++) {
	  qbt1 *= exphairpin[u];
	  if ((tetra_loop)&&(u==4)) {
	    char tl[7]={0}, *ts;
	    strncpy(tl, sequences[s]+i-1, 6);
	    if ((ts=strstr(Tetraloops, tl)))
	      qbt1 *= exptetra[(ts-Tetraloops)/7];
	  } 
	  if (u==3) {
	    char tl[6]={0,0,0,0,0,0}, *ts;
	    strncpy(tl, sequences[s]+i-1, 5);
	    if ((ts=strstr(Triloops, tl))) 
	      qbt1 *= expTriloop[(ts-Triloops)/6];
	    if (type[s]>2) 
	      qbt1 *= expTermAU;
	  }
	  else /* no mismatches for tri-loops */
	    qbt1 *= expmismatchH[type[s]][S[s][i+1]][S[s][j-1]];
	}
	qbt1 *= scale[u+2];

	/* interior loops with interior pair k,l */
	for (k=i+1; k<=MIN(i+MAXLOOP+1,j-TURN-2); k++) {
	  u1 = k-i-1;
	  for (l=MAX(k+TURN+1,j-1-MAXLOOP+u1); l<=j-1; l++) {
	    double qloop=1;
	    if (qb[iindx[k]-l]==0) {qloop=0; continue;}
	    for (s=0; s<n_seq; s++) {
	      type_2 = pair[S[s][l]][S[s][k]]; if (type_2 == 0) type_2 = 7;
	      qloop *= expLoopEnergy(u1, j-l-1, type[s], type_2,
				     S[s][i+1], S[s][j-1],
				     S[s][k-1], S[s][l+1]);
	    }
	    qbt1 += qb[iindx[k]-l] * qloop * scale[u1+j-l-1+2];
	  }
	}
	/* multi-loop loop contribution */
	ii = iindx[i+1]; /* ii-k=[i+1,k-1] */
	temp = 0.0;
	for (k=i+2; k<=j-1; k++) temp += qm[ii-(k-1)]*qqm1[k];
	for (s=0; s<n_seq; s++) {
	  tt = rtype[type[s]];
	  temp *= expMLintern[tt]*expMLclosing*
	    expdangle3[tt][S[s][i+1]]*expdangle5[tt][S[s][j-1]];
	}
	qbt1 += temp*scale[2];
	qb[ij] = qbt1;
	qb[ij] *= exp(psc/kTn);
      } /* end if (type!=0) */
      else qb[ij] = 0.0;
       
      /* construction of qqm matrix containing final stem
	 contributions to multiple loop partition function
	 from segment i,j */
      qqm[i] = qqm1[i]*expMLbase[1];  /* expMLbase[1]^n_seq */
      for (qbt1=1, s=0; s<n_seq; s++) {
	qbt1 *= expMLintern[type[s]];
	if (i>1) qbt1 *= expdangle5[type[s]][S[s][i-1]];
	if (j<n) qbt1 *= expdangle3[type[s]][S[s][j+1]];
	else if (type[s]>2) qbt1 *= expTermAU;
      }
      qqm[i] += qb[ij]*qbt1;      
      /* construction of qm matrix containing multiple loop
	 partition function contributions from segment i,j */
      temp = 0.0;
      ii = iindx[i];  /* ii-k=[i,k-1] */
      for (k=i+1; k<=j; k++)
	temp += (qm[ii-(k-1)]+expMLbase[k-i])*qqm[k];
      qm[ij] = (temp + qqm[i]);
      
      /* auxiliary matrix qq for cubic order q calculation below */
      qbt1 = qb[ij];
      if (qbt1>0)
	for (s=0; s<n_seq; s++) {
	  if (i>1) qbt1 *= expdangle5[type[s]][S[s][i-1]];
	  if (j<n) qbt1 *= expdangle3[type[s]][S[s][j+1]];
	  else if (type[s]>2) qbt1 *= expTermAU;
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
	if (qb[ij]>0.) {
	  pr[ij] = q1k[i-1]*qln[j+1]/q1k[n] * exp(pscore[ij]/kTn);
	  for (s=0; s<n_seq; s++) {
	    int typ;
	    typ = pair[S[s][i]][S[s][j]]; if (typ==0) typ=7;
	    if (i>1) pr[ij] *= expdangle5[typ][S[s][i-1]];
	    if (j<n) pr[ij] *= expdangle3[typ][S[s][j+1]];
	    else if (typ>2) pr[ij] *= expTermAU;
	  }
	} else
	  pr[ij] = 0;
      }
    }
      
    for (l=n; l>TURN+1; l--) {

      /* 2. bonding k,l as substem of 2:loop enclosed by i,j */
      for (k=1; k<l-TURN; k++) {
	double pp = 0;
	kl = iindx[k]-l;
	if (qb[kl]==0) continue;
	for (s=0; s<n_seq; s++) {
          type[s] = pair[S[s][l]][S[s][k]];
          if (type[s]==0) type[s]=7;
        }
	
	for (i=MAX(1,k-MAXLOOP-1); i<=k-1; i++) 
	  for (j=l+1; j<=MIN(l+ MAXLOOP -k+i+2,n); j++) {
	    ij = iindx[i] - j;
	    if ((pr[ij]>0)) {
	      double qloop=1;
	      for (s=0; s<n_seq; s++) {
		int typ;
		typ = pair[S[s][i]][S[s][j]]; if (typ==0) typ=7;
		qloop *=  expLoopEnergy(k-i-1, j-l-1, typ, type[s],
				S[s][i+1], S[s][j-1], S[s][k-1], S[s][l+1]);
	      }
	      pp += pr[ij]*qloop*scale[k-i-1 + j-l-1 + 2];
	    } 
	  }
	pr[kl] += pp * exp(pscore[kl]/kTn);
      }
      /* 3. bonding k,l as substem of multi-loop enclosed by i,j */
      prm_MLb = 0.;
      if (l<n) for (k=2; k<l-TURN; k++) {
	i = k-1;
	prmt = prmt1 = 0.0;
	    
	ii = iindx[i];     /* ii-j=[i,j]     */
	ll = iindx[l+1];   /* ll-j=[l+1,j-1] */
	prmt1 = pr[ii-(l+1)];
	for (s=0; s<n_seq; s++) {
	  tt = pair[S[s][l+1]][S[s][i]]; if (tt==0) tt=7;
	  prmt1 *= expMLclosing*expMLintern[tt]*
	    expdangle3[tt][S[s][i+1]]*expdangle5[tt][S[s][l]];
	} 
	for (j=l+2; j<=n; j++) {
	  double pp=1;
	  if (pr[ii-j]==0) continue;
	  for (s=0; s<n_seq; s++) {
	    tt=pair[S[s][j]][S[s][i]]; if (tt==0) tt=7;
	    pp *=  expdangle3[tt][S[s][i+1]]*
	      expdangle5[tt][S[s][j-1]];
	  }
	  prmt +=  pr[ii-j]*pp*qm[ll-(j-1)];
	}
	kl = iindx[k]-l;
	for (s=0; s<n_seq; s++) {
	  int typ;
	  typ=pair[S[s][k]][S[s][l]]; if (typ==0) typ=7;
	  prmt *= expMLclosing*expMLintern[typ];
	}
	prml[ i] = prmt;
	prm_l[i] = prm_l1[i]*expMLbase[1]+prmt1; /* expMLbase[1]^n_seq */
	
	prm_MLb = prm_MLb*expMLbase[1] + prml[i];
	/* same as:    prm_MLb = 0;
	   for (i=1; i<=k-1; i++) prm_MLb += prml[i]*expMLbase[k-i-1]; */
	
	prml[i] = prml[ i] + prm_l[i];
	
	if (qb[kl] == 0.) continue; 
	    
	temp = prm_MLb;
	
	for (i=1;i<=k-2; i++) 
	  temp += prml[i]*qm[iindx[i+1] - (k-1)];

	for (s=0; s<n_seq; s++) {
	  tt=pair[S[s][k]][S[s][l]]; if (tt==0) tt=7;
	  temp *= expMLintern[tt];
	  if (k>1) temp *= expdangle5[tt][S[s][k-1]];
	  if (l<n) temp *= expdangle3[tt][S[s][l+1]];
	  else if (tt>2) temp *= expTermAU;
	}
	pr[kl] += temp * scale[2] * exp(pscore[kl]/kTn);
	
#ifndef LARGE_PF
	if (pr[kl]>Qmax) {
	  Qmax = pr[kl];
	  if (Qmax>FLT_MAX/10.)
	    fprintf(stderr, "%d %d %g %g\n", i,j,pr[kl],qb[kl]);
	}
	if (pr[kl]>FLT_MAX) {
	  ov++;
	  pr[kl]=FLT_MAX;
	}
#endif
      } /* end for (k=..) */
      tmp = prm_l1; prm_l1=prm_l; prm_l=tmp;

    }  /* end for (l=..)   */
    
    for (i=1; i<=n; i++)
      for (j=i+TURN+1; j<=n; j++) {
	ij = iindx[i]-j;
	pr[ij] *= qb[ij] *exp(-pscore[ij]/kTn);
      }

    if (pi != NULL)
      *pi = make_pairinfo((const short **)S, sequences, n_seq);
  
    if (structure!=NULL)
      sprintf_bppm(n, structure);
  }   /* end if (do_backtrack)*/

  for (i=0; i<n_seq; i++) free(S[i]);
  free(S);
  if (ov>0) fprintf(stderr, "%d overflows occurred while backtracking;\n"
		    "you might try a smaller pf_scale than %g\n",
		    ov, pf_scale);
  free(type);
  free_alipf_arrays();
  return (free_energy); 
}

/*------------------------------------------------------------------------*/
/* dangling ends should never be destabilizing, i.e. expdangle>=1         */
/* specific heat needs smooth function (2nd derivative)                   */
/* we use a*(sin(x+b)+1)^2, with a=2/(3*sqrt(3)), b=Pi/6-sqrt(3)/2,       */
/* in the interval b<x<sqrt(3)/2                                          */

#define SCALE 10
#define SMOOTH(X) ((X)/SCALE<-1.2283697)?0:(((X)/SCALE>0.8660254)?(X):\
          SCALE*0.38490018*(sin((X)/SCALE-0.34242663)+1)*(sin((X)/SCALE-0.34242663)+1))

PRIVATE void scale_pf_params(unsigned int length, int n_seq)
{
  /* scale energy parameters and pre-calculate Boltzmann weights */
  unsigned int i, j, k, l;
  double  kT, TT, kTn;
  double  GT;

   
   
  kT = (temperature+K0)*GASCONST;   /* kT in cal/mol  */
  kTn = kT*n_seq;
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
    exphairpin[i] = exp( -GT*10./kTn);
  }
  for (i=0; i<=MIN(30, MAXLOOP); i++) {
    GT =  bulge37[i]*TT;
    expbulge[i] = exp( -GT*10./kTn);
    GT =  internal_loop37[i]*TT;
    expinternal[i] = exp( -GT*10./kTn);
  }
  /* special case of size 2 interior loops (single mismatch) */
  if (james_rule) expinternal[2] = exp ( -80*10/kTn);
   
  lxc = lxc37*TT;
  for (i=31; i<length; i++) {
    GT = hairpin37[30]*TT + (lxc*log( i/30.));
    exphairpin[i] = exp( -GT*10./kTn);
  }
  for (i=31; i<=MAXLOOP; i++) {
    GT = bulge37[30]*TT + (lxc*log( i/30.));
    expbulge[i] = exp( -GT*10./kTn);
    GT = internal_loop37[30]*TT + (lxc*log( i/30.));
    expinternal[i] = exp( -GT*10./kTn);
  }

  for (i=0; i<5; i++) {
    GT = F_ninio37[i]*TT;
    for (j=0; j<=MAXLOOP; j++)
      expninio[i][j]=exp(-MIN(MAX_NINIO,j*GT)*10/kTn);
  }
  for (i=0; (i*7)<strlen(Tetraloops); i++) {
    GT = TETRA_ENTH37 - (TETRA_ENTH37-TETRA_ENERGY37[i])*TT;
    exptetra[i] = exp( -GT*10./kTn);
  }
  for (i=0; (i*5)<strlen(Triloops); i++) 
    expTriloop[i] = exp(-Triloop_E37[i]*10/kTn);

  GT =  ML_closing37*TT;
  expMLclosing = exp( -GT*10/kTn);

  for (i=0; i<=NBPAIRS; i++) { /* includes AU penalty */
    GT =  ML_intern37*TT;
    /* if (i>2) GT += TerminalAU; */
    expMLintern[i] = exp( -GT*10./kTn);
  }
  expTermAU = exp(-TerminalAU*10/kTn);

  GT =  ML_BASE37*TT;
  for (i=0; i<length; i++) {
    expMLbase[i] = exp( -10.*i*GT/kT)*scale[i];
  }

  /* if dangles==0 just set their energy to 0,
      don't let dangle energies become > 0 (at large temps),
      but make sure go smoothly to 0                        */
  for (i=0; i<=NBPAIRS; i++)
    for (j=0; j<=4; j++) {
      GT = dangle5_H[i][j] - (dangle5_H[i][j] - dangle5_37[i][j])*TT;
      expdangle5[i][j] = dangles?exp(SMOOTH(-GT)*10./kTn):1.;
      GT = dangle3_H[i][j] - (dangle3_H[i][j] - dangle3_37[i][j])*TT;
      expdangle3[i][j] =  dangles?exp(SMOOTH(-GT)*10./kTn):1.;
      if (i>2) /* add TermAU penalty into dangle3 */
	expdangle3[i][j] *= expTermAU;
    }

  /* stacking energies */
  for (i=0; i<=NBPAIRS; i++)
    for (j=0; j<=NBPAIRS; j++) {
      GT =  enthalpies[i][j] - (enthalpies[i][j] - stack37[i][j])*TT;
      expstack[i][j] = exp( -GT*10/kTn);
    }

  /* mismatch energies */
  for (i=0; i<=NBPAIRS; i++)
    for (j=0; j<5; j++)
      for (k=0; k<5; k++) {
	GT = mism_H[i][j][k] - (mism_H[i][j][k] - mismatchI37[i][j][k])*TT;
	expmismatchI[i][j][k] = exp(-GT*10.0/kTn);
	GT = mism_H[i][j][k] - (mism_H[i][j][k] - mismatchH37[i][j][k])*TT;
	expmismatchH[i][j][k] = exp(-GT*10.0/kTn);
	GT = mism_H[i][j][k] - (mism_H[i][j][k] - mismatchM37[i][j][k])*TT;
	expmismatchM[i][j][k] = exp(-GT*10.0/kTn);
      }

  /* interior lops of length 2 */
  for (i=0; i<=NBPAIRS; i++)
    for (j=0; j<=NBPAIRS; j++)
      for (k=0; k<5; k++)
	for (l=0; l<5; l++) {
	  GT = int11_H[i][j][k][l] -
	    (int11_H[i][j][k][l] - int11_37[i][j][k][l])*TT;
	  expint11[i][j][k][l] = exp(-GT*10./kTn);
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
	    expint21[i][j][k][l][m] = exp(-GT*10./kTn);
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
	      expint22[i][j][k][l][m][n] = exp(-GT*10./kTn);
	    }
	}  
}

/*----------------------------------------------------------------------*/

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
  return z;
}
 
/*----------------------------------------------------------------------*/

PRIVATE void get_arrays(unsigned int length)
{
  unsigned int size,i;
   
  size = sizeof(FLT_OR_DBL) * ((length+1)*(length+2)/2);
  q   = (FLT_OR_DBL *) space(size);
  qb  = (FLT_OR_DBL *) space(size);
  qm  = (FLT_OR_DBL *) space(size);
  pscore = (short *) space(sizeof(short)*((length+1)*(length+2)/2));
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
  for (i=1; i<=length; i++) {
    iindx[i] = ((length+1-i)*(length-i))/2 +length+1;
  }
}

/*----------------------------------------------------------------------*/
   
PUBLIC void init_alipf_fold(int length, int n_seq)
{
  if (length<1) nrerror("init_pf_fold: length must be greater 0");
  if (init_length>0) free_alipf_arrays(); /* free previous allocation */
#ifdef SUN4
  nonstandard_arithmetic();
#else
#ifdef HP9
  fpsetfastmode(1);
#endif
#endif
  make_pair_matrix();
  get_arrays((unsigned) length);
  scale_pf_params((unsigned) length, n_seq);
  init_length=length;
}

PRIVATE void free_alipf_arrays(void)
{
  free(q);
  free(qb);
  free(qm);
  free(pscore);
  free(qq); free(qq1);
  free(qqm); free(qqm1);
  free(q1k); free(qln);
  free(prm_l); free(prm_l1); free(prml);
  free(exphairpin);
  free(expMLbase);
  free(scale);
  free(iindx);
#ifdef SUN4
  standard_arithmetic();
#else
#ifdef HP9
  fpsetfastmode(0);
#endif
#endif
  init_length=0;
}
/*---------------------------------------------------------------------------*/
#define PMIN 0.0008
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

pair_info *make_pairinfo(const short *const* S, char **AS, int n_seq) {
  int i,j,n, num_p=0, max_p = 64;
  pair_info *pi;
  double *duck, p;
  n = S[0][0];
  max_p = 64; pi = space(max_p*sizeof(pair_info));
  duck =  (double *) space((n+1)*sizeof(double));
  for (i=1; i<n; i++)
    for (j=i+TURN+1; j<=n; j++) 
      if ((p=pr[iindx[i]-j])>0) {
	duck[i] -=  p * log(p);
	duck[j] -=  p * log(p);
      }

  for (i=1; i<n; i++)
    for (j=i+TURN+1; j<=n; j++) {
      if ((p=pr[iindx[i]-j])>=PMIN) {
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
	  pi = realloc(pi, max_p * sizeof(pair_info));
	  if (pi==NULL) nrerror("out of memory in alipf_fold");
	}
      }  
    }
  free(duck);
  pi = realloc(pi, (num_p+1)*sizeof(pair_info));
  pi[num_p].i=0;
  qsort(pi, num_p, sizeof(pair_info), compare_pair_info );
  return pi;
}
/*---------------------------------------------------------------------------*/

#define L 3
PRIVATE void sprintf_bppm(int length, char *structure)
{
  extern char  bppm_symbol(float *x);
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

PRIVATE short * encode_seq(const char *sequence) {
  unsigned int i,l;
  short *S;
  l = strlen(sequence);
  S = (short *) space(sizeof(short)*(l+1));
  S[0] = (short) l;
  
  /* make numerical encoding of sequence */
  for (i=1; i<=l; i++) 
    S[i]= (short) encode_char(toupper(sequence[i-1]));
  
  return S;
}

/*---------------------------------------------------------------------------*/

PRIVATE void make_pscores(const short *const *S, const char *const *AS,
			  int n_seq, const char *structure) {
  /* calculate co-variance bonus for each pair depending on  */
  /* compensatory/consistent mutations and incompatible seqs */
  /* should be 0 for conserved pairs, >0 for good pairs      */
#define NONE -10000 /* score for forbidden pairs */
  int n,i,j,k,l,s,score;
  int dm[7][7]={{0,0,0,0,0,0,0}, /* hamming distance between pairs */
                {0,0,2,2,1,2,2} /* CG */,
                {0,2,0,1,2,2,2} /* GC */,
                {0,2,1,0,2,1,2} /* GU */,
                {0,1,2,2,0,2,1} /* UG */,
                {0,2,2,1,2,0,2} /* AU */,
                {0,2,2,2,1,2,0} /* UA */};

  n=S[0][0];  /* length of seqs */
  for (i=1; i<n; i++) {
    for (j=i+1; (j<i+TURN+1) && (j<=n); j++) 
      pscore[iindx[i]-j] = NONE;
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
      if (pfreq[0]*2>n_seq) { pscore[iindx[i]-j] = NONE; continue;}
      for (k=1,score=0; k<=6; k++) /* ignore pairtype 7 (gap-gap) */
	for (l=k+1; l<=6; l++) 
	  /* scores for replacements between pairtypes    */
	  /* consistent mutations (if l==k+2) score 1     */
	  /* compensatory (all else) score 2              */
	  score += pfreq[k]*pfreq[l]*dm[k][l];
      /* counter examples score -1, gap-gap scores -0.25  */
      pscore[iindx[i]-j] = cv_fact *
	((UNIT*score)/n_seq - nc_fact*UNIT*(pfreq[0] + pfreq[7]*0.25));
    }
  }

  if (noLonelyPairs) /* remove unwanted pairs */
    for (k=1; k<=n-TURN-1; k++) 
      for (l=1; l<=2; l++) {
	int type,ntype=0,otype=0;
	i=k; j = i+TURN+l;
	type = pscore[iindx[i]-j];
	while ((i>=1)&&(j<=n)) {
	  if ((i>1)&&(j<n)) ntype = pscore[iindx[i-1]-j-1];
	  if ((otype<cv_fact*MINPSCORE)&&(ntype<cv_fact*MINPSCORE))
	    /* too many counterexamples */
	    pscore[iindx[i]-j] = NONE; /* i.j can only form isolated pairs */
	  otype =  type;
	  type  = ntype;
	  i--; j++;
	}
      }
  

  if (fold_constrained&&(structure!=NULL)) {
    int psij, hx, *stack;
    stack = (int *) space(sizeof(int)*(n+1));
    
    for(hx=0, j=1; j<=n; j++) {
      switch (structure[j-1]) {
      case 'x': /* j can't pair */ 
        for (l=1; l<j-TURN; l++) pscore[iindx[l]-j] = NONE;
        for (l=j+TURN+1; l<=n; l++) pscore[iindx[j]-l] = NONE;
        break;
      case '(':
        stack[hx++]=j;
        /* fallthrough */
      case '<': /* j pairs upstream */
        for (l=1; l<j-TURN; l++) pscore[iindx[l]-j] = NONE;
        break;
      case ')': /* j pairs with i */
        if (hx<=0) {
          fprintf(stderr, "%s\n", structure);
          nrerror("unbalanced brackets in constraints");
        }
        i = stack[--hx];
	psij = pscore[iindx[i]-j]; /* store for later */
        for (l=i; l<=j; l++) 
	  for (k=j; k<=n; k++) pscore[iindx[l]-k] = NONE;
	for (k=1; k<=i; k++) 
	  for (l=i; l<=j; l++) pscore[iindx[k]-l] = NONE;
	for (k=i+1; k<j; k++)
	  pscore[iindx[k]-j] = pscore[iindx[i]-k] = NONE;
	pscore[iindx[i]-j] = (psij>0) ? psij : 0;
        /* fallthrough */
      case '>': /* j pairs downstream */
        for (l=j+TURN+1; l<=n; l++) pscore[iindx[j]-l] = NONE;
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
