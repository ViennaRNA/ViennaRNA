/* Last changed Time-stamp: <1999-05-06 11:58:09 ivo> */
/*                
	    partiton function for RNA secondary structures

			    Ivo L Hofacker
	       Sebastian Bonhoeffer, John S McCaskill,
			 and Peter F Stadler
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

static char rcsid[] = "$Id: part_func.c,v 1.8 1999/05/06 09:58:27 ivo Exp $";

#define MAX(x,y) (((x)>(y)) ? (x) : (y))
#define MIN(x,y) (((x)<(y)) ? (x) : (y))
#define PUBLIC
#define PRIVATE static
#define STACK_BULGE1  1   /* stacking energies for bulges of size 1 */
#define NEW_NINIO     1   /* new asymetry penalty */


PUBLIC  float pf_fold(char *sequence, char *structure);
PUBLIC  void  init_pf_fold(int length);
PUBLIC  void  free_pf_arrays(void);
PUBLIC  void  update_pf_params(int length);
PUBLIC  char  bppm_symbol(float *x);
PRIVATE void  sprintf_bppm(int length, char *structure);
PRIVATE void  scale_pf_params(int length);
PRIVATE void  get_arrays(int length);
PRIVATE void  BP_calculate(char *structure, int *BP, int length);

PRIVATE FLT_OR_DBL expMLclosing, expMLintern, *expMLbase;
PRIVATE FLT_OR_DBL expdangle5[NBPAIRS+1][5], expdangle3[NBPAIRS+1][5];
PRIVATE FLT_OR_DBL lxc, exptetra[40], expTriloop[40];
PRIVATE FLT_OR_DBL expstack[NBPAIRS+1][NBPAIRS+1];
PRIVATE FLT_OR_DBL expmismatchI[NBPAIRS+1][5][5], expmismatchH[NBPAIRS+1][5][5];
PRIVATE FLT_OR_DBL expSint2[NBPAIRS+1][NBPAIRS+1][5][5];
PRIVATE FLT_OR_DBL *exphairpin;
PRIVATE FLT_OR_DBL expbulge[MAXLOOP+1];
PRIVATE FLT_OR_DBL expinternal[MAXLOOP+1];
PRIVATE FLT_OR_DBL expninio[5][MAXLOOP+1];
PRIVATE FLT_OR_DBL *q, *qb, *qm, *qqm, *qqm1, *qq, *qq1;
PRIVATE FLT_OR_DBL *prml, *prm_l, *prm_l1, *q1k, *qln;
PRIVATE FLT_OR_DBL *scale;
PRIVATE char msg[200];
PRIVATE int *jindx;
PRIVATE int init_length; /* length in last call to init_pf_fold() */
#define ISOLATED  256.0

/*-----------------------------------------------------------------*/
PUBLIC float pf_fold(char *sequence, char *structure)
{
   short *S, *S1;
   int n, i,j,k,l,  ij,kl, u,u1,u2,d,ii,ll, type, type_2, tt, ov=0;
   FLT_OR_DBL temp, Q, Qmax=0, prm_MLb;
   FLT_OR_DBL prmt,prmt1;
   FLT_OR_DBL qbt1,s1temp, *tmp;
   
   float free_energy;
   char *pos;

   int   *BP; /* contains the structure constrainsts: BP[i]
		 negative: x = base must not pair
		 positive int: base is paired with int      */
   int  cpos, cforbid;

   n = strlen(sequence);
   if (n>init_length) init_pf_fold(n);  /* (re)allocate space */
   
   if (fold_constrained) {
      if (structure==NULL) nrerror("no constraints given!");
      BP = (int *)space(sizeof(int)*(n+2));
      BP_calculate(structure,BP,n);
   }
    
   S = (short *) space(sizeof(short)*(n+1));
   S1= (short *) space(sizeof(short)*(n+1));

   for (l=1; l<=n; l++) {
      if (energy_set>0) S[l]=sequence[l-1]-'A'+1;
      else {
	 pos = strchr(Law_and_Order, sequence[l-1]);
	 if (pos==NULL) S[l]=0;
	 else S[l]= pos-Law_and_Order;
      }
      S1[l] = alias[S[l]];
   }
   
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

   if (noLonelyPairs) { /* mark isolated pairs */
     for (k=2; k<n-TURN-1; k++) 
       for (l=1; l<=2; l++) {
	 int ntype,otype=0;
	 i = k; j = i+TURN+l;
	 type = pair[S[i]][S[j]];
	 while ((i>1)&&(j<n)) {
	   ntype = pair[S[i-1]][S[j+1]];
	   if (type&&(!otype)&&(!ntype)) /* i.j can only form isolated pairs */
	     qb[iindx[i]-j] = ISOLATED;
	   else qb[iindx[i]-j] = 0.;
	   otype =  type;
	   type  = ntype;
	   i--; j++;
	 }
       }
   } 

   for (j=TURN+2;j<=n; j++) {
      if (fold_constrained) {
	 cforbid=0; cpos=n+1;
	 for (k=j; k>j-TURN-1; k--) {
	    if (BP[k]>j) cforbid=1;         /* no pair (i,j) i<k allowed */
	    if (BP[j]<0) cforbid=1;         /* j doesn't pair */
	    if ((BP[k]<k-TURN)&&(BP[k]>0))
	       cpos=((cpos<BP[k])?cpos:BP[k]);    /* no pair (i,j) i>=cpos */
	 }
      }
      for (i=j-TURN-1; i>=1; i--) {
	 /* construction of partition function of segment i,j*/
	 /*firstly that given i bound to j : qb(i,j) */
	 u = j-i-1;
	 type = pair[S[i]][S[j]];
	 if (noLonelyPairs && (qb[iindx[i]-j]==ISOLATED))
	   type=0; /* don't allow this pair */
	   
	 if (fold_constrained) {
	    if (BP[i]>j) cforbid=1;         /* no more pairs for this j */
	    if ((BP[i]<i)&&(BP[i]>0))
	      cpos=((cpos<BP[i])?cpos:BP[i]);  /* no pairs while i>=cpos */

	    if (BP[i]==j) { 
	       if (type==0) type=7;         /* constrained unusual pair */
	       if ((cforbid)||(i>cpos))
		  nrerror("conflicting constraints");
	       cforbid=1;                   /* no pairs for smaller i */
	    } else
	       if ((BP[i]<0)||cforbid||(i>=cpos)) type=0;
	 }
	 
	 if (type!=0) {
	    /*hairpin contribution*/ 
	    qbt1 = exphairpin[u];
	    if ((tetra_loop)&&(u==4)) {
	      char tl[5]={0,0,0,0,0}, *ts;
	      strncpy(tl, sequence+i, 4);
	      if ((ts=strstr(Tetraloops, tl)))
		qbt1 *= exptetra[(ts-Tetraloops)/5];
	    } 
	    if (u==3) {
	      char tl[6]={0,0,0,0,0,0}, *ts;
	      strncpy(tl, sequence+i-1, 5);
	      if ((ts=strstr(Triloops, tl))) 
		qbt1 *= expTriloop[(ts-Tetraloops)/6];
	    }
	    else /* no mismatches for tri-loops */
	      qbt1 *= expmismatchH[type][S1[i+1]][S1[j-1]];
	    
	    qbt1 *= scale[u+2];

	    /* interior loops with interior pair k,l */
	    for (k=i+1; k<=MIN(i+MAXLOOP+1,j-TURN-2); k++) {
	       u1 = k-i-1;
	       for (l=MAX(k+TURN+1,j-1-MAXLOOP+u1); l<=j-1; l++) {
		 
		 type_2 = pair[S[k]][S[l]];
		 if ((fold_constrained)&&(BP[k]==l)&&(type_2==0))
		   type_2=7; /* nonstandard */
		 if (type_2) {
		   register int rtype2;
		   rtype2  = rtype[type_2];
		   u2 = j-l-1;
		   
		   if (u1==0) {
		     if (u2==0) /* stack */
		       s1temp = expstack[type][type_2];
		     else { /* bulge */
		       s1temp = expbulge[u2];
#if STACK_BULGE1
		       if (u2==1) s1temp *= expstack[type][type_2];
#endif
		     }
		   } else {
		     if (u2==0) { /* bulge */
		       s1temp = expbulge[u1];
#if STACK_BULGE1
		       if (u1==1) s1temp *= expstack[type][type_2];
#endif
		     }
		     else { /* interior loop */
		       if ((u1+u2==2)&&james_rule) /* size 2 is special */
			 s1temp = expSint2[type][rtype2][S1[i+1]][S1[j-1]];
		       else {
			 s1temp = expinternal[u1+u2]*
			   expmismatchI[type][S1[i+1]][S1[j-1]]*
			   expmismatchI[rtype2][S1[l+1]][S1[k-1]];
#if NEW_NINIO
			 s1temp *= expninio[2][abs(u1-u2)];
#else
			 m = MIN(4,MIN(u1,u2));
			 s1temp *= expninio[m][abs(u1-u2)];
#endif
		       }
		     }
		   }
		   qbt1 += qb[iindx[k]-l]*s1temp*scale[u1+u2+2];
		 }
	       }
	    }
	    /*multiple stem loop contribution*/
	    ii = iindx[i+1]; /* ii-k=[i+1,k-1] */
	    temp = 0.0;
	    for (k=i+2; k<=j-1; k++) temp += qm[ii-(k-1)]*qqm1[k]; 
	    tt = rtype[type];
	    qbt1 += temp*expMLclosing*expMLintern*scale[2]*
	       expdangle3[tt][S1[i+1]]*expdangle5[tt][S1[j-1]];

	    qb[iindx[i]-j] = qbt1;
	 } /* end if (type!=0) */
	 else qb[iindx[i]-j] = 0.0;
	 
	 /* construction of qqm matrix containing final stem
	    contributions to multiple loop partition function
	    from segment i,j */
	 ij = iindx[i]-j;
	 qqm[i] = qqm1[i]*expMLbase[1];
	 if (type) {
	    qbt1 = qb[ij]*expMLintern;
	    if (i>1) qbt1 *= expdangle5[type][S1[i-1]];
	    if (j<n) qbt1 *= expdangle3[type][S1[j+1]];
	    qqm[i] += qbt1;
	 }
	 
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
	 }
	 qq[i] = qq1[i]*scale[1] + qbt1;
	 
	 /*construction of partition function for segment i,j */
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
	    type = pair[S[i]][S[j]];
	    ij = iindx[i]-j;
	    if (fold_constrained)
	       if ((BP[i]==j)&&(type==0)) type=7;
	    if (type&&(qb[ij]>0.)) {
	       pr[ij] = q1k[i-1]*qln[j+1]/q1k[n];
	       if (i>1) pr[ij] *= expdangle5[type][S1[i-1]];
	       if (j<n) pr[ij] *= expdangle3[type][S1[j+1]];
	    } else
	       pr[ij] = 0;
	 }
      }
      
      for (l=n; l>TURN+1; l--) {

	 /* 2. bonding k,l as substem of 2:loop enclosed by i,j */
	 for (k=1; k<l-TURN; k++) {
	    kl = iindx[k]-l;
	    type_2 = pair[S[k]][S[l]];
	    if (qb[kl]==0) continue;
	
	    for (i=MAX(1,k-MAXLOOP-1); i<=k-1; i++) 
	      for (j=l+1; j<=MIN(l+ MAXLOOP -k+i+2,n); j++) {
		ij = iindx[i] - j;
		type = pair[S[i]][S[j]];
		if ((pr[ij]>0)) {
		  u1 = k-i-1;
		  u2 = j-l-1;
		  if (u1==0) {
		    if (u2==0) /* stack */
		      s1temp = expstack[type][type_2];
		    else {      /* bulge */
		      s1temp = expbulge[u2];
#if STACK_BULGE1
		      if (u2==1) s1temp *= expstack[type][type_2];
#endif
		    }
		  } else {
		    if (u2==0) { /* bulge */
		      s1temp = expbulge[u1];
#if STACK_BULGE1
		      if (u1==1) s1temp *= expstack[type][type_2];
#endif
		    }
		    else {     /* interior loop */
		      register int rtype2;
		      rtype2  = rtype[type_2];
		      if ((u1+u2==2)&&james_rule) /* size 2 is special */
			s1temp = expSint2[type][rtype2][S1[i+1]][S1[j-1]];
		      else {
			s1temp = expinternal[u1+u2]*
			  expmismatchI[type][S1[i+1]][S1[j-1]]*
			  expmismatchI[rtype2][S1[l+1]][S1[k-1]];
#if NEW_NINIO
			s1temp *= expninio[2][abs(u1-u2)];
#else
			m = MIN(4,MIN(u1,u2));
			s1temp *= expninio[m][abs(u1-u2)];
#endif
		      }
		    }
		  }
		  pr[kl] += (pr[ij])*s1temp*scale[u1+u2+2];
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
	    tt = pair[S[l+1]][S[i]];
	    prmt1 = pr[ii-(l+1)]*expMLclosing*expMLintern*
	       expdangle3[tt][S1[i+1]]*expdangle5[tt][S1[l]];
	    for (j=l+2; j<=n; j++) {
	       tt = pair[S[j]][S[i]];
	       prmt += pr[ii-j]*expdangle3[tt][S1[i+1]]*
		 expdangle5[tt][S1[j-1]] *qm[ll-(j-1)];
	    }
	    prmt *= expMLclosing*expMLintern;
	    prml[ i] = prmt;
	    prm_l[i] = prm_l1[i]*expMLbase[1]+prmt1;

	    prm_MLb = prm_MLb*expMLbase[1] + prml[i];
	    /* same as:    prm_MLb = 0;
	    for (i=1; i<=k-1; i++) prm_MLb += prml[i]*expMLbase[k-i-1]; */

	    prml[i] = prml[ i] + prm_l[i];
	    
	    kl = iindx[k]-l;
	    ll = jindx[l];

	    if (qb[kl] == 0.) continue; 
	    
	    temp = prm_MLb;

	    for (i=1;i<=k-2; i++) 
	       temp += prml[i]*qm[iindx[i+1] - (k-1)];

	    tt = pair[S[k]][S[l]];
	    temp *= expMLintern*scale[2];
	    if (k>1) temp *= expdangle5[tt][S1[k-1]];
	    if (l<n) temp *= expdangle3[tt][S1[l+1]];
	    pr[kl] += temp;
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
	    pr[ij] *= qb[ij];
	 }
      
      if (structure!=NULL)
	 sprintf_bppm(n, structure);
   }   /* end if (do_backtrack)*/
   
   free(S);
   free(S1);
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

PRIVATE void scale_pf_params(int length)
{
  /* scale energy parameters and pre-calculate Boltzmann weights */
   int i, j, k, l;
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
      exphairpin[i] = exp( -GT*10/kT);
   }
   for (i=0; i<=MIN(30, MAXLOOP); i++) {
      GT =  bulge37[i]*TT;
      expbulge[i] = exp( -GT*10/kT);
      GT =  internal_loop37[i]*TT;
      expinternal[i] = exp( -GT*10/kT);
   }
   /* special case of size 2 interior loops (single mismatch) */
   if (james_rule) expinternal[2] = exp ( -internal2_energy*10/kT);
   
   lxc = lxc37*TT;
   for (i=31; i<length; i++) {
      GT = hairpin37[30]*TT + (lxc*log( i/30.));
      exphairpin[i] = exp( -GT*10.0/kT);
   }
   for (i=31; i<=MAXLOOP; i++) {
      GT = bulge37[30]*TT + (lxc*log( i/30.));
      expbulge[i] = exp( -GT*10/kT);
      GT = internal_loop37[30]*TT + (lxc*log( i/30.));
      expinternal[i] = exp( -GT*10/kT);
   }

   for (i=0; i<5; i++) {
      GT = F_ninio37[i]*TT;
      for (j=0; j<=MAXLOOP; j++)
	 expninio[i][j]=exp(-MIN(MAX_NINIO,j*GT)*10/kT);
   }
   for (i=0; (i*5)<strlen(Tetraloops); i++) {
     GT = TETRA_ENTH37 - (TETRA_ENTH37-TETRA_ENERGY37[i])*TT;
     exptetra[i] = exp( -GT*10/kT);
   }
   for (i=0; (i*5)<strlen(Triloops); i++) 
     expTriloop[i] = exp(-Triloop_E37[i]*10/kT);

   GT =  ML_closing37*TT;
   expMLclosing = exp( -GT*10/kT);

   GT =  ML_intern37*TT;
   expMLintern = exp( -GT*10/kT);
   
   GT =  ML_BASE37*TT;
   for (i=0; i<length; i++) {
      expMLbase[i] = exp( -i*GT*10/kT)*scale[i];
   }

   /* if dangles==0 just set their energy to 0,
      don't let dangle energies become > 0 (at large temps),
      but make sure go smoothly to 0                        */
   for (i=0; i<=NBPAIRS; i++)
      for (j=0; j<=4; j++) {
	GT = dangle5_H[i][j] - (dangle5_H[i][j] - dangle5_37[i][j])*TT;
	expdangle5[i][j] = dangles?exp(SMOOTH(-GT)*10/kT):1.;
	GT = dangle3_H[i][j] - (dangle3_H[i][j] - dangle3_37[i][j])*TT;
	expdangle3[i][j] =  dangles?exp(SMOOTH(-GT)*10/kT):1.;
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
       }

   /* interior lops of length 2 */
   for (i=0; i<=NBPAIRS; i++)
     for (j=0; j<=NBPAIRS; j++)
       for (k=0; k<5; k++)
         for (l=0; l<5; l++) {
           GT = Sint2_H[i][j][k][l] -
	     (Sint2_H[i][j][k][l] - Sint2_37[i][j][k][l])*TT;
	   expSint2[i][j][k][l] = exp(-GT*10./kT);
	 }
}

/*----------------------------------------------------------------------*/

PRIVATE void get_arrays(int length)
{
   unsigned int size,i;
   
   size = sizeof(FLT_OR_DBL) * ((length+1)*(length+2)/2);
   q   = (FLT_OR_DBL *) space(size);
   qb  = (FLT_OR_DBL *) space(size);
   qm  = (FLT_OR_DBL *) space(size);
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
   get_arrays(length);
   scale_pf_params(length);
   init_length=length;
}

PUBLIC void free_pf_arrays(void)
{
   free(q);
   free(qb);
   free(qm);
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
}
/*---------------------------------------------------------------------------*/

PUBLIC void update_pf_params(int length)
{
   make_pair_matrix();
   scale_pf_params(length);
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

PRIVATE void BP_calculate(char *structure, int *BP, int length)
{
   int i,j,ct;
   /* only ".x()" are recognized in pf_fold() */
   
   for(i=0;i<length;i++) {
      switch (structure[i]) {
       case 'x': BP[i+1] = -4; break;
       case '(': 
	 ct = 1;
	 for(j=i+1;j<length;j++) {
	    if(structure[j]==')') ct--;
	    if(structure[j]=='(') ct++;
	    if(ct==0) {
	       BP[i+1]=j+1;
	       BP[j+1]=i+1;
	       break;
	    }
	 }
	 break;
      }
   }
}
