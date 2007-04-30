/* Last changed Time-stamp: <2007-04-30 12:44:11 ulim> */
/*                
		  partiton function for RNA secondary structures

		  Ivo L Hofacker

		  Vienna RNA package
*/
/*
  $Log: part_func_up.c,v $
  Revision 1.1  2007/04/30 15:13:13  ivo
  merge RNAup into package

  Revision 1.11  2006/07/17 11:11:43  ulim
  removed all globals from fold_vars.h,c, cleaned code

  Revision 1.10  2006/07/12 09:19:29  ulim
  global variables w, incr3 and incr5 are now local

  Revision 1.9  2006/07/11 12:45:02  ulim
  remove redundancy in function pf_interact(...)

  Revision 1.8  2006/03/08 15:26:37  ulim
  modified -o[1|2], added meaningful default

  Revision 1.5  2006/01/23 11:27:04  ulim
  include file into new package version. cleaned it

  Revision 1.2  2005/07/29 15:13:37  ulim
  put the function, calculating the probability of an unpaired region in
  an RNA and the function calculating the prob. of interaction between 2 RNAs
  in a seperate file (pf_two.c)

  Revision 1.1  2005/07/26 13:27:12  ulim
  Initial revision

  Revision 1.2  2005/07/01 13:14:57  ulim
  fixed error in scaling, included new commandline options -incr5, -incr3 to
  allow a variable number of unpaired positions 5' and 3' of the site of
  interaction between the two RNAs

  Revision 1.1  2005/04/19 08:16:38  ulim
  Initial revision
*/

#include <config.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>    /* #defines FLT_MAX ... */
#include "fold.h"
#include "utils.h"
#include "energy_par.h"
#include "fold_vars.h"
#include "pair_mat.h"
#include "params.h"
#include "part_func.h"
#include "part_func_up.h"
#include "duplex.h"
/*@unused@*/
static char rcsid[] UNUSED = "$Id: part_func_up.c,v 1.1 2007/04/30 15:13:13 ivo Exp $";

#define MAX(x,y) (((x)>(y)) ? (x) : (y))
#define MIN(x,y) (((x)<(y)) ? (x) : (y))
#define EQUAL(A,B) (fabs((A)-(B)) < 1000*DBL_EPSILON)
#define PUBLIC
#define PRIVATE static
/* #define PF2_DEBUG 1 */

PUBLIC  pu_contrib *pf_unstru(char *sequence, char *structure, int w);
PUBLIC  void free_pf_unstru(void);
PUBLIC  double **pf_interact(const char *s1, const char *s2, pu_contrib *p_c, int w, int incr3, int incr5);
/* free_pf_two: first argument output of pf_unstru() !!!! */ 
PUBLIC  void free_pf_two(pu_contrib *p_con, double **p_in);
PUBLIC int Up_plot(pu_contrib *p_c, double **pint, int len, char *ofile, int w, char *select_contrib);
PUBLIC double expLoopEnergy(int u1, int u2, int type, int type2, short si1, short sj1, short sp1, short sq1);
PUBLIC double expHairpinEnergy(int u, int type, short si1, short sj1, const char *string);


PRIVATE void  init_pf_two(int length);
/* allocate/ free space for pf_unpaired() and pf_up() */
PRIVATE void scale_int(const char *s, const char *sl,double *sc_int,int incr3, int incr5);
PRIVATE void get_unpaired(int length);
PRIVATE void free_unpaired(void);
PRIVATE void encode_seq(const char *s1, const char *s2);

PRIVATE pf_paramT *Pf = NULL;/* use this structure for all the exp-arrays*/
PRIVATE FLT_OR_DBL *qb, *qm, *prpr, *puij;/* add arrays for pf_unpaired()*/
PRIVATE FLT_OR_DBL *q1k, *qln;
PRIVATE double *qqm2, *qq_1m2, *qqm, *qqm1;

PRIVATE double *scale, *expMLbase;
PRIVATE char *ptype; /* precomputed array of pair types */ 
PRIVATE int init_length;  /* length in last call to init_pf_fold()*/
PRIVATE double init_temp; /* temperature in last call to scale_pf_params */

#define ISOLATED  256.0

/*-----------------------------------------------------------------*/
static  short *S, *S1, *SS, *SS2;
PUBLIC pu_contrib *pf_unstru(char *sequence, char *structure, int w)
{
  int n, i,j,k,l, ij, kl, u,u1,d,type, type_2, tt, uu;
  int o,p,po,x,y;
  unsigned int size,nu;
  double temp, tqm2, bla; 
  double qbt1, *tmp, Zuij;
  pu_contrib *pu_stru;

  temp=0;
  n = (int) strlen(sequence);
  /* contributions to probability of being unpaired witihin a(n)
     H hairpin,
     I interior loop,
     M muliloop,
     E exterior loop*/
  size = sizeof(double *) * (n+1);
  
  pu_stru = (pu_contrib *) space (sizeof(pu_contrib)*1);
  pu_stru->length=n;
  /* pu_stru->X[i][j] where i <= j and i [1...n], j = [1...w[ */ 
  pu_stru->H = (double **) space(size);
  pu_stru->I = (double **) space(size);
  pu_stru->M = (double **) space(size);
  pu_stru->E = (double **) space(size);
  for(i=0;i<=n;i++){
    pu_stru->H[i] = (double *) space(sizeof(double)*(w+1));
    pu_stru->I[i] = (double *) space(sizeof(double)*(w+1));
    pu_stru->M[i] = (double *) space(sizeof(double)*(w+1));
    pu_stru->E[i] = (double *) space(sizeof(double)*(w+1));
  }
     
  init_pf_two(n);
   
  for (d=0; d<=TURN; d++) 
    for (i=1; i<=n-d; i++) {
      j=i+d;
      ij = iindx[i]-j;
      puij[ij]=0.0;
      if(d < w) {
	pu_stru->H[i][d]=pu_stru->I[i][d]=pu_stru->M[i][d]=pu_stru->E[i][d]=0.0;
      }
    }
  nu = (unsigned int) n;
  size = ((nu+1)*(nu+2)/2);
  for (i=0; i<size; i++){
    prpr[i]=pr[i];
  }
  
  for (i=1; i<=n; i++)
    for (j=i+TURN+1; j<=n; j++) {
      ij = iindx[i]-j;
      /* i need the part_func of all structures outside bp[ij] */
      if(qb[ij] > 0.0) prpr[ij]=pr[ij]/qb[ij];
    }
  
  for (i=1; i<=n; i++)
    {
      qqm2[i]=qq_1m2[i]=0;
    }
  
  /* no do_backtrack test anymore ! */
  for (i=1; i<=n; i++)
    {
      /* set auxillary arrays to 0, reuse qqm and qqm1*/
      qqm[i]=qqm1[i]=0;
    }
  /* 2. exterior bp (p,o) encloses unpaired region ]i-1,i+w[*/
  for (o=TURN+2;o<=n; o++) {
    for (p=o-TURN-1; p>=1; p--) {
      /* construction of partition function of segment [p,o], given that
	 an unpaired region ]i-1,i+w[ exists within [p,o]: qpq_w[i]*/
      /*firstly that given p bound to o :
	qb(p,o)) and region ]i-1,i+w[ is unpaired*/
      u = o-p-1; po = iindx[p]-o;
      type = ptype[po];
      if (type!=0) {  
	/*hairpin contribution*/
 	temp = expHairpinEnergy(u, type, S1[p+1], S1[o-1], sequence+p-1) * prpr[po] * scale[u+2]; /* add scale[u+2] */
	/* p < i <= o-2 region [i,j] ,= w */
	for(i=p+1; i<=o-1;i++) {
	  uu=0;
	  for(j=i; j < MIN(i+w,o); j++) {
	    ij=iindx[i]-j;
	    if (((type==3)||(type==4))&&no_closingGU) {
	      puij[ij]+= 0.;
	      pu_stru->H[i][j-i]+=0.;
	    }
	    else {
	      puij[ij]+= temp;
	      pu_stru->H[i][j-i]+=temp;
	    }
	  }
	}
      }
      /* interior loops with interior pair k,l and an unpaired region of
	 length w between p and k || l and o*/
      for (k=p+1; k<=MIN(p+MAXLOOP+1,o-TURN-2); k++) {
	u1 = k-p-1;
	for (l=MAX(k+TURN+1,o-1-MAXLOOP+u1); l<o; l++) {
	  kl=iindx[k]-l;
	  type_2 = ptype[kl];
	  if (type_2) {
	    type_2 = rtype[type_2];
	    /* add *scale[u1+u2+2] */ 
	    temp = qb[kl]  * (scale[u1+o-l+1] *
		    expLoopEnergy(u1, o-l-1, type, type_2,
   				  S1[p+1], S1[o-1], S1[k-1], S1[l+1])) *
	      prpr[po];
	    /* p < i <= k-2, region [i,j] <= w */
	    for (i=p+1; i <= k-1; i++ ){
	      uu=0;
	      for (j=i; j < MIN(i+w,k); j++ ){
		ij=iindx[i]-j;
		puij[ij]+=temp;
		pu_stru->I[i][j-i]+=temp;
	      }
	    }

	    /* p < k < l < i <= o-2, region [i,j] <= w */
	    for (i=l+1; i <= o-1; i++){
	      uu=0;
	      for (j=i; j < MIN(i+w,o); j++ ){
		ij=iindx[i]-j;
		puij[ij]+=temp;
		pu_stru->I[i][j-i]+=temp;
	      }
	    }	
	  }
	} /* end of l */
      } /* end of k */
      
      /* multiple stem loop contribution
	 calculate qm2[iindx[i]-j] in the course of the calculation
	 of the multiple stem loop contribution:
	 advantage: you save memory:
	 instead of a (n+1)*n array for qqm2 you only need 2*n arrays
	 disadvantage: you have to use two times the op-loop for the full
	 multiloop contribution
	 first op-loop: index o goes from 1...n and
	                index p from o-TURN-1 ... 1
	 second op-loop: index o goes from n...1 and
	                 index p from o+TURN+1 ... n !!
	 HERE index o goes from 1...n and index p o-TURN-1 ... 1 ,
	 we calculate the contributions to multiple stem loop 
	 where exp(i+w-1-p)*(qqm2 values between i+w and o-1)
	 AND qm[iindex[p+1]-(i-1)]*exp(beta*w)*qm[iindex[i+w]-(o-1)]
	 you have to recalculate of qqm matrix containing final stem
	 contributions to multiple loop partition function
	 from segment p,o */

      /* recalculate qqm[] 
	 qqm[p] := (contribution with exact one loop in region (p,o)*/
      qqm[p] = qqm1[p]*expMLbase[1];
      if (type) {
	qbt1 = qb[po]*Pf->expMLintern[type];
	if (p>1) qbt1 *= Pf->expdangle5[type][S1[p-1]];
	if (o<n) qbt1 *= Pf->expdangle3[type][S1[o+1]];
	else if (type>2) qbt1 *= Pf->expTermAU;
	qqm[p] += qbt1;
	/* revers dangles for prpr[po]*... */
	temp=0.; 
	tt=rtype[type];
	temp = prpr[po]*Pf->expdangle3[tt][S1[p+1]]*Pf->expdangle5[tt][S1[o-1]];
	temp *=Pf->expMLclosing*Pf->expMLintern[tt]*scale[2];
      }
      tqm2=0.;
      for(i=p+1; i < o; i++) {
	uu=0;
	tqm2+=qm[iindx[p]-i]*qqm[i+1];
	
	if(  type !=0 && (i <= o - 1) ) {
	  for(j=i; j< MIN(i+w,o);j++) {
	    int iwo, p1i;
	    double temp2;
	    temp2=0;
	    iwo=p1i=0;
	    ij=iindx[i]-j;
	    /* iindx[x] - y : x <= y */
	    x=j+1; y=o-1;
	    if( x < y) iwo=iindx[j+1]-(o-1);
	    x=p+1; y=i-1;
	    if( x < y ) p1i=iindx[p+1]-(i-1);
	    /* w unpaired bases left of structured region */
	    temp2 = expMLbase[j-p]*qq_1m2[j+1];
	    /* w unpaired bases between structured region */
	    temp2 += qm[p1i]*expMLbase[j-i+1]*qm[iwo];
	      
	    /* multiply with revers dangles for prpr[po]... */
	    temp2 *= temp;
	    puij[ij]+=temp2;
	    pu_stru->M[i][j-i]+=temp2;
	  } /* end of j ...  */
	}
      }/*end of for i ... */
      /* qqm2[p] contrib with at least 2 loops in region (p,o) */ 
      qqm2[p]=tqm2;
    } /* end for (p=..) */
    tmp = qqm1; qqm1=qqm; qqm=tmp;
    tmp = qqm2; qqm2=qq_1m2; qq_1m2=tmp;
  }/* end for (o=..) */

  for (i=1; i<=n; i++) {
    /* set auxillary arrays to 0 */
    qqm[i]=qqm1[i]=0;
    qqm2[i]=qq_1m2[i]=0;
  }
    
  /* 2. exterior bp (p,o) encloses unpaired region [i,j]
     HERE index o goes from n...1 and index p from o+TURN+1 ... n,
     that is, we add the one multiloop contribution that we
     could not calculate before  */
  for (o=n-TURN-1;o>=1; o--) {
    for (p=o+TURN+1; p<=n; p++) {
      po=iindx[o]-p;
      type=ptype[po];
      /* recalculate of qqm matrix containing final stem
	 contributions to multiple loop partition function
	 from segment [o,p] */
      qqm[p] = qqm1[p]*expMLbase[1];
      if (type) {
	qbt1 = qb[po]*Pf->expMLintern[type];
	if (o>1) qbt1 *= Pf->expdangle5[type][S1[o-1]];
	if (p<n) qbt1 *= Pf->expdangle3[type][S1[p+1]];
	else if (type>2) qbt1 *= Pf->expTermAU;
	qqm[p] += qbt1;
	/* revers dangles for prpr[po]...  */
	temp=0.;
	tt=rtype[type];
	temp = prpr[po]*Pf->expdangle3[tt][S1[o+1]]*Pf->expdangle5[tt][S1[p-1]];
	temp *= Pf->expMLclosing*Pf->expMLintern[tt]*scale[2];
      }
      tqm2=0.;
      for(i=o+1; i < p; i++) {
	uu=0;
	tqm2+=qqm[i]*qm[iindx[i+1]-p];
	
	if(  type !=0 && i <= p-1 ) {
	  double temp2;
	  temp2=0;
	         
	  /* w unpaired bases right of structured region */
	  temp2 = qq_1m2[i-1]*expMLbase[p-i];
	  /* multiply with revers dangles for prpr[po]*... */
	  temp2 *= temp;
	  
	  for (j=i; j<MIN(i+w,p);j++) {
	    ij=iindx[i]-j;
	    puij[ij]+=temp2;
	    pu_stru->M[i][j-i]+=temp2;
	  }
	}
      }/*end of for i ....*/
      qqm2[p]=tqm2;
    }/* end for (p=..) */
    tmp = qqm1; qqm1=qqm; qqm=tmp;
    tmp = qqm2; qqm2=qq_1m2; qq_1m2=tmp;
  }/* end for (o=..) */

  /* 1. region [i,j] exterior to all loops */
  Zuij=0.;bla=0;
  for (i=1; i<=n; i++) {
    uu=0;
    for(j=i; j<MIN(i+w,n+1);j++){
      ij=iindx[i]-j;
      temp=q1k[i-1]*1*scale[j-i+1]*qln[j+1]/q1k[n];
      puij[ij]+=temp;
      pu_stru->E[i][j-i]+=temp;
      
      Zuij+=puij[ij];/* partition function over all contributions to puij*/
    }
    bla+=puij[ij];
  }
#ifdef PF2_DEBUG
  /* different tests */
  printf("pf_unpaired Zuij=%.3f   Z=%.3f\n",Zuij,q1k[n]); 
  /* test if puij[iindx[i]-i] == 1 - p_paired von [i],
     where p_paired von [i] = sum_i pr[ij oder ji]*/
  float bla_pf;
  for (i=1; i<=n; i++) {
    uu=iindx[i]-i;
    bla=puij[uu];
    bla_pf=.0;
    for(j=1; j<=n;j++){
      if(i < j)	{
	ij=iindx[i]-j;
      }
      else if( i > j){
	ij=iindx[j]-i;
      }
      else{
	continue;
      }	
      bla_pf+=pr[ij];
    }
    printf("%3d p_unstru[i]: %12.9f\t1-p_paired[i]: %12.9f\n",i,bla,(1-bla_pf));
  }
  /* check if pu_stru->H[ij]+pu_stru->I[ij]+pu_stru->M[ij]+pu_stru->E[ij] == puij[ij]  */
  for (i=1; i<n; i++) {
    for (j=i; j < MIN((i+w),n); j++) {
      double blubb;
      ij = iindx[i]-j;
      blubb = pu_stru->H[i][j-i]+pu_stru->I[i][j-i]+pu_stru->M[i][j-i]+pu_stru->E[i][j-i];
      if(!EQUAL(puij[ij],blubb)) {
	printf("i=%d,j=%d,ij=%d\npuij = %.18f\nblubb= %.18f\n",i,j,ij,puij[ij],blubb);
      }
    }
  }
#endif
  
  puij[0]=Zuij;
  
  return pu_stru;  
}
/*------------------------------------------------------------------------*/
/* s1 is the longer seq */
PUBLIC  double **pf_interact(const char *s1, const char *s2, pu_contrib *p_c, int w, int incr3, int incr5)
{
  int i, j, k,l,n1,n2,add_i5,add_i3;
  double temptest, Ztest, temp, Z, rev_d, E, Z2,**p_c_S, int_scale;
  FLT_OR_DBL ****qint_4, **qint_ik, **p_intik;
  double Fup, temppfs, *temps; /* get free energy*/
  unsigned int size;
  PRIVATE double **pint;  /*array for pf_up() output */
  
  n1 = (int) strlen(s1); 
  n2 = (int) strlen(s2);
   
  p_c_S = (double **) space (sizeof(double *)*(n1+1));
  
  for (i=1; i<=n1; i++) {
    p_c_S[i] = (double *) space (sizeof(double)*(w+incr3+incr5+1));    
    for (j=0; j < (w+incr3+incr5); j++) {
      p_c_S[i][j] = p_c->H[i][j]+p_c->I[i][j]+p_c->M[i][j]+p_c->E[i][j];
      
    }
  }

  /*array for pf_up() output */
  pint = (double **) space(sizeof(double*)*(2));
  for(i=0; i<=1;i++){
    pint[i] = (double *) space(sizeof(double)*(n1+2));
  }
    
  /* use a different scaling for pf_interact*/
  scale_int(s2,s1,&int_scale,incr3,incr5);
  /* set the global scale array and the global variable pf_scale to the
     values used to scale the interaction, keep their former values !! */
  temppfs = pf_scale;
  pf_scale = int_scale;
  /* in order to scale expLoopEnergy correctly call*/
  init_pf_fold(n1);
  scale_stru_pf_params((unsigned) n1);
  
  qint_ik = (FLT_OR_DBL **) space(sizeof(FLT_OR_DBL *) * (n1+1));
  for (i=1; i<=n1; i++) {
    qint_ik[i] = (FLT_OR_DBL *) space(sizeof(FLT_OR_DBL) * (n1+1));
  }
  
  encode_seq(s1, s2);
  /* static  short *S~S1, *S1~SS1, *SS~S2, *SS2; */
  for (i=1; i<=n1; i++) {
    pint[0][i]=pint[1][i]=0.;
  }
  E=0.;
  Z=0.;

  /*  qint_4[i][j][k][l] contribution that region (k-i) in seq1 (l=n1)
      is paired to region (l-j) in seq 2(l=n2) that is
      a region closed by bp k-l  and bp i-j */
  
  qint_4 = (FLT_OR_DBL ****) space(sizeof(FLT_OR_DBL ***) * (n1+1));
  
  /* qint_4[i][j][k][l] */
  for (i=1; i<=n1; i++) {
    /* note: qint_4[i] will be freed before we allocate qint_4[i+1] */
    qint_4[i] = (FLT_OR_DBL ***) space(sizeof(FLT_OR_DBL **) * (n2+1));
    for (j=n2; j>0; j--) {
      qint_4[i][j] = (FLT_OR_DBL **) space(sizeof(FLT_OR_DBL*) * (w+1));
      for (k=0; k<=w; k++) {
	qint_4[i][j][k] = (FLT_OR_DBL *) space(sizeof(FLT_OR_DBL) * (w+1));
      }
    }
    for (j=n2; j>0; j--) {
      int type, type2,ii,ki;
      type = pair[S[i]][SS[j]];
      
      qint_4[i][j][0][0] = type ? Pf->expDuplexInit : 0;

      if (!type) continue;
      if (i>1)  qint_4[i][j][0][0] *= Pf->expdangle5[type][S1[i-1]];
      if (j<n2 && dangles) qint_4[i][j][0][0] *= Pf->expdangle3[type][SS2[j+1]];
      else if (type>2) qint_4[i][j][0][0] *= Pf->expTermAU;
      
      rev_d=1.;
      if (i<n1 && dangles) rev_d *= Pf->expdangle3[rtype[type]][S1[i+1]];
      else if (type>2) rev_d *= Pf->expTermAU;
      if (j>1) rev_d *= Pf->expdangle5[rtype[type]][SS2[j-1]];
      
      /* add inc5 and incr3 */
      if((i-incr5) > 0 ) add_i5=i-incr5;
      else add_i5=1;
      add_i3=incr3;

      /* only one bp (no interior loop) */
      qint_ik[i][i]+=qint_4[i][j][0][0]*rev_d*p_c_S[add_i5][add_i3]*scale[((int) w/2)];
      Z+=qint_4[i][j][0][0]*rev_d*p_c_S[add_i5][add_i3]*scale[((int) w/2)];

      temp=0.;
      for (k=i-1; k>0 && k>i-MAXLOOP-2 && k>i-w; k--) {
	for (l=j+1; l<j+w && l<=n2; l++) {
	  int a,b,ia,ib,isw;
	  double scalew, tt;
	  
	  type2 = pair[S[k]][SS[l]];
	  if (!type2) continue;
	  
	  /* to save memory keep only qint_4[i-w...i][][][] in memory 
	     use indices qint_4[i][j][a={0,1,...,w-1}][b={0,1,...,w-1}] */
	  a=i-k;/* k -> a from 1...w-1*/
	  b=l-j;/* l -> b from 1...w-1 */
	    
	  /* scale everything to w/2 */
	  isw = ((int) w/2);
	  if ((a+b) < isw ){
	    scalew = ( scale[isw - (a+b)] );
	  } else if ( (a+b) > isw ) {
	    scalew = 1/( scale[(a+b) - isw] );
	  } else {
	    scalew = 1;
	  }
	  
	  if (i-k+l-j-2<=MAXLOOP) {
            /* exLoopEnergy is scaled correctly */
	    E = expLoopEnergy(i-k-1,l-j-1, type2, rtype[type],
			      S1[k+1], SS2[l-1], S1[i-1], SS2[j+1]) *
	                      scale[i-k+l-j]; /* add *scale[u1+u2+2] */
 
	    qint_4[i][j][a][b] += ( qint_4[k][l][0][0]*E);
	    
	    /* use ia and ib to go from a....w-1 and from b....w-1  */
	    ia=ib=1;
	    while((a+ia)<w && i-(a+ia)>=1 && (b+ib)<w && (j+b+ib)<=n2) {
	      int iaa,ibb;

	      qint_4[i][j][a+ia][b+ib] += qint_4[k][l][ia][ib]*E;

	      iaa=ia+1;
	      while(a+iaa<w && i-(a+iaa)>=1) {
		qint_4[i][j][a+iaa][b+ib] += qint_4[k][l][iaa][ib]*E;
		++iaa;
	      }
	      
	      ibb=ib+1;
	      while( (b+ibb)<w && (j+b+ibb)<=n2 ) {
		qint_4[i][j][a+ia][b+ibb] += qint_4[k][l][ia][ibb]*E;
		++ibb;
	      } 
	      ++ia;
	      ++ib;
	    }
	  }
	  /* scale everything to w/2*/
	  /* qint_ik[k][i] all interactions where k and i both are paired */
	  /* substract incr5 from k */
          if(k-incr5 > 0) add_i5=k-incr5;
          else add_i5=1;
	  /* add incr3 to i */
          add_i3=i-k+incr3;
	  tt = qint_4[i][j][a][b]*p_c_S[add_i5][add_i3]*scalew*rev_d; 
	  temp+= tt;
	  qint_ik[k][i]+= tt; 
	}
      }
      Z+=temp;
    }
  
    /* free qint_4 values not needed any more */
    if(i > w) {
      int bla;
      bla=i-w;
      for (j=n2; j>0; j--) {
	for (k=0; k<=w; k++){
	  free(qint_4[bla][j][k]);
	}
	free(qint_4[bla][j]);
      }
      free(qint_4[bla]);
    }
  }

  
  Z2=0.0;
  for (i=1; i<=n1; i++) {
    for (k=i; k<=n1 && k<i+w; k++) {
      Z2+=qint_ik[i][k];
      for(l=i;l<=k;l++) {
	double delG;
	
	/* pint[0]: prob that position l is within a paired region */
	/* qint_ik[i][k] as well as Z are scaled to scale[((int) w/2) */
	pint[0][l]+=qint_ik[i][k]/Z;
	/* pint[1]: minimal delta G at position [l] */
	pint[1][l]=MIN(pint[1][l],
			( -log(qint_ik[i][k])-( ((int) w/2)*log(pf_scale)) )*
			 ((Pf->temperature+K0)*GASCONST/1000.0) );
	delG=( -log(qint_ik[i][k])-( ((int) w/2)*log(pf_scale)) )*((Pf->temperature+K0)*GASCONST/1000.0);
      }
    }
  }
  
  if(n1 > w){
    for (i=n1-w+1; i<=n1; i++) {
      for (j=n2; j>0; j--) {
	for (k=0; k<=w; k++) {
	  free(qint_4[i][j][k]);
	}
	free(qint_4[i][j]);
      }
      free(qint_4[i]);
    }
    free(qint_4);
  } else {
    for (i=1; i<=n1; i++) {
      for (j=n2; j>0; j--) {
	for (k=0; k<=w; k++) {
	  free(qint_4[i][j][k]);
	}
	free(qint_4[i][j]);
      }
      free(qint_4[i]);
    }
    free(qint_4);
  }
  
  for (i=1; i<=n1; i++) {
    free(qint_ik[i]);
  }
  free(qint_ik);
  
  Fup = (-log(Z)-n1*log(pf_scale))*(Pf->temperature+K0)*GASCONST/1000.0;

  /* reset the global variables pf_scale and scale to their original values */
  pf_scale = temppfs;/* reset pf_scale */
  scale_stru_pf_params((unsigned) n1);/* reset the scale array */
  free_pf_arrays(); /* for arrays for pf_fold(...) */
  
  if(expMLbase != NULL) {
    free(expMLbase);
    expMLbase = NULL;
  }
  if(scale != NULL) {
    free(scale);
    scale = NULL;
  }
  for (i=1; i<=n1; i++) {
    free(p_c_S[i]);      
  }
  free(p_c_S);
  
  return pint;  
}
/*------------------------------------------------------------------------*/
/* use an extra scale for pf_interact, here sl is the longer sequence */
PRIVATE void scale_int(const char *s, const char *sl, double *sc_int,int incr3, int incr5)
{
  int i,n,nl,l_scales;
  float e_duplex;
  duplexT mfe;
  double  kT;
  
  /* sc_int is similar to pf_scale: i.e. one time the scale */
  n=strlen(s);
  nl=strlen(sl);
  l_scales = (2*nl)+incr3+incr5;

  expMLbase  = (double *) space(sizeof(double)*(nl+1));
  scale = (double *) space(sizeof(double)*((nl+1)*2));
  
  /* use RNA duplex to get a realistic estimate for the best possible
     interaction energy between the short RNA s and its target sl */
  mfe = duplexfold(s,sl);
  
  kT = (Pf->temperature+K0)*GASCONST/1000.0;   /* in Kcal */

  *sc_int = 3.42;
  
  /* sc_int is similar to pf_scale: i.e. one time the scale */
  *sc_int = exp(-(mfe.energy)/kT/n);
  
  /* free the structure returned by duplexfold */
  free(mfe.structure);  
}

/*----------------------------------------------------------------------*/
/* init_pf_two(n) :gets the arrays, that you need, from part_func.c */ 
/* get_pf_arrays(&S, &S1, &ptype, &qb, &qm, &q1k, &qln); get_unpaired(n); */
/* init_pf_fold(), update_pf_params, encode_char(), make_ptypes() are called by pf_fold() */  
PRIVATE void init_pf_two(int length)
{
  
#ifdef SUN4
  nonstandard_arithmetic();
#else
#ifdef HP9
  fpsetfastmode(1);
#endif 
#endif
  expMLbase  = (double *) space(sizeof(double)*(length+1));
  scale = (double *) space(sizeof(double)*((length+1)*2));
  make_pair_matrix();
  /* gets the arrays, that you need, from part_func.c */
  /* pr: base pairing prob. matrix, global declared in fold_vars.h */
  /* iindx also global declared in fold_vars.h */
  if(get_pf_arrays(&S, &S1, &ptype, &qb, &qm, &q1k, &qln) == 0) {
    nrerror("init_pf_two:"
	    "pf_fold() has to be called before calling pf_unstru()\n");
  }
  get_unpaired(length);
  scale_stru_pf_params((unsigned) length);
  
  init_length=length;
  if(init_temp != Pf->temperature)
    nrerror("init_pf_two: inconsistency with temperature");  
}

/*---------------------------------------------------------------------------*/
PRIVATE void get_unpaired(int length)
{
  unsigned int size;
   
  size = sizeof(FLT_OR_DBL) * ((length+1)*(length+2)/2);
  prpr = (FLT_OR_DBL *) space(size);
  puij = (FLT_OR_DBL *) space(size);
  qqm2 = (double *) space(sizeof(double)*(length+2));
  qq_1m2 = (double *) space(sizeof(double)*(length+2));
  qqm = (double *) space(sizeof(double)*(length+2));
  qqm1 = (double *) space(sizeof(double)*(length+2));
}
PRIVATE void free_unpaired(void)
{
  free(puij);
  puij=NULL;
  free(prpr);
  prpr=NULL;
  free(qqm2);
  free(qq_1m2);
  free(qqm);
  free(qqm1);
}
/* free all but the output structure for pf_unstru */
PUBLIC void free_pf_unstru(void) {
  if(scale != NULL) {
    free(scale);
    scale = NULL;
  }
  if(expMLbase != NULL) {
    free(expMLbase);
    expMLbase = NULL;
  }
  free_unpaired();
}
/*---------------------------------------------------------------------------*/

/* first argument output of pf_unstru() !!!! */ 
PUBLIC void free_pf_two(pu_contrib *p_con, double **pin)
{
  int i;  
  if(p_con != NULL) {
    for(i=0;i<=p_con->length;i++) {
      free(p_con->H[i]);
      free(p_con->I[i]);
      free(p_con->M[i]);
      free(p_con->E[i]);
    }
    free(p_con->H);
    free(p_con->I);
    free(p_con->M);
    free(p_con->E);
    free(p_con);
    p_con = NULL;
  }
  if(S != NULL && pin != NULL){
    free(S);
    S=NULL;
  }
  if(S1 != NULL && pin != NULL){
    free(S1);
    S1=NULL;
  }
  if(SS != NULL){
    free(SS);
    SS=NULL;
  }
  if(SS2 != NULL){
    free(SS2);
    SS2=NULL;
  }
  if(pin != NULL){
    for(i=0; i<=1;i++){
      free(pin[i]);
    }
    free(pin);
    pin=NULL;
  }
}
/*---------------------------------------------------------------------------*/

PRIVATE void encode_seq(const char *s1, const char *s2) {
  unsigned int i,l;

  l = strlen(s1);
  /* S and S1 are freed by free_pf_arrays(); ! */
  S = (short *) space(sizeof(short)*(l+1));
  S1= (short *) space(sizeof(short)*(l+1));
  /* S1 exists only for the special X K and I bases and energy_set!=0 */
  S[0] = l;
  for (i=1; i<=l; i++) { /* make numerical encoding of sequence */
    S[i]= (short) encode_char(toupper(s1[i-1]));
    S1[i] = alias[S[i]];   /* for mismatches of nostandard bases */
  }
  l = strlen(s2);
  SS = (short *) xrealloc(SS, sizeof(short)*(l+1));
  SS2= (short *) xrealloc(SS2, sizeof(short)*(l+1));
  /* SS2 exists only for the special X K and I bases and energy_set!=0 */
  SS[0] = l;
  for (i=1; i<=l; i++) { /* make numerical encoding of sequence */
    SS[i]= (short) encode_char(toupper(s2[i-1]));
    SS2[i] = alias[SS[i]];   /* for mismatches of nostandard bases */
  }
}

/*-------------------------------------------------------------------------*/
 /* scale energy parameters and pre-calculate Boltzmann weights:
  most of this is done in structure Pf see params.c,h (function:
  scale_pf_parameters(), only arrays scale and expMLbase are handled here*/
PUBLIC void scale_stru_pf_params(unsigned int length)
{
  unsigned int i;
  double  kT;


  /* Do this only at the first call for scale_pf_parameters()
     and/or if temperature has changed*/
  if(init_temp != temperature) {
    Pf=scale_pf_parameters();
  }

  init_temp = Pf->temperature;
  kT = (Pf->temperature+K0)*GASCONST;   /* kT in cal/mol  */

  /* scaling factors (to avoid overflows) */
  if (pf_scale==-1) { /* mean energy for random sequences: 184.3*length cal */
    pf_scale = exp(-(-185+(Pf->temperature-37.)*7.27)/kT);
    if (pf_scale<1) pf_scale=1;
  }
  scale[0] = 1.;
  for (i=1; i<=(length*2); i++) {
    scale[i] = scale[i-1]/pf_scale;
  }

  for (i=0; i<length; i++) {
    expMLbase[i] = exp(i*Pf->expMLbase)*scale[i];
  }

}
/*-------------------------------------------------------------------------*/

PUBLIC int Up_plot(pu_contrib *p_c, double **pint, int len, char *ofile, int w, char *select_contrib)
{
  /* produce output */
  FILE *wastl;
  int i,j,g,uu;
  float mini,maxi;
  double **p_cont; /* p_u AND its contributions */
  char *unpaired[6]={"total probability of being unpaired","probability of being unpaired in the exterior loop","probability of being unpaired in hairpin loops","probability of being unpaired in interior loops","probability of being unpaired in multiloops","interaction probability"};
  
  if(p_c != NULL) {
    wastl = fopen(ofile,"w");
    if (wastl==NULL) {
      fprintf(stderr, "p_cont: can't open %s for Up_plot\n", ofile);
      perror(ofile);
      exit(EXIT_FAILURE);
    }
  } else if (pint != NULL) {
    wastl = fopen(ofile,"a");
    if (wastl==NULL) {
      fprintf(stderr, "p_cont: can't open %s for Up_plot\n", ofile);
      perror(ofile);
      exit(EXIT_FAILURE);
    }
  } else {
    return(0);
  }
  
  /* collect output values
     p_cont[0] sum of pr_unpaired H+I+M+E
     p_cont[1] E exterior loops
     p_cont[2] H hairpin loops
     p_cont[3] I interior loops
     p_cont[4] M multi loops
  */ 
  p_cont = NULL;
  if(p_c != NULL) {
    p_cont = (double **) space(sizeof(double*)*(5));
    
    if(strchr(select_contrib, 'S')) {
      p_cont[0] = (double *) space(sizeof(double)*(len+3));
      p_cont[0][0] = 1.;
    } else {
      p_cont[0] = (double *) space(sizeof(double)*(3));
      p_cont[0][0] = 0.;
    }
    if(strchr(select_contrib, 'E')) {
      p_cont[1] = (double *) space(sizeof(double)*(len+3));
      p_cont[1][0] = 1.;
    } else {
      p_cont[1] = (double *) space(sizeof(double)*(3));
      p_cont[1][0] = 0.;
    }
    if(strchr(select_contrib, 'H')) {
      p_cont[2] = (double *) space(sizeof(double)*(len+3));
      p_cont[2][0] = 1.;
    } else {
      p_cont[2] = (double *) space(sizeof(double)*(3));
      p_cont[2][0] = 0.;
    }
    if(strchr(select_contrib, 'I')) {
      p_cont[3] = (double *) space(sizeof(double)*(len+3));
      p_cont[3][0] = 1.;
    } else {
      p_cont[3] = (double *) space(sizeof(double)*(3));
      p_cont[3][0] = 0.;
    }
    if(strchr(select_contrib, 'M')) {
      p_cont[4] = (double *) space(sizeof(double)*(len+3));
      p_cont[4][0] = 1.;
    } else {
      p_cont[4] = (double *) space(sizeof(double)*(3));
      p_cont[4][0] = 0.;
    }

    for (i=1; i<=len+2; i++) {
      /* not necessary for calculation, stores contributions to 
	 prob. upaired[i], i gives the START of the unpaired region*/
      if((int) p_cont[0][0] ) p_cont[0][i]=0.;
      if((int) p_cont[1][0] ) p_cont[1][i]=0.;
      if((int) p_cont[2][0] ) p_cont[2][i]=0.;
      if((int) p_cont[3][0] ) p_cont[3][i]=0.;
      if((int) p_cont[4][0] ) p_cont[4][i]=0.; 
    }
    
    for (i=1; i<=len; i++) {
      for (j=i; j < MIN((i+w),len+1); j++) {
	double blubb;
	if( (j-i+1) == w && i+w-1 <= len) {
	  blubb = p_c->H[i][j-i]+p_c->I[i][j-i]+p_c->M[i][j-i]+p_c->E[i][j-i];
	  if((int) p_cont[0][0] ) p_cont[0][i+w-1]+=blubb;
	  if((int) p_cont[1][0] ) p_cont[1][i+w-1]+=p_c->E[i][j-i];
	  if((int) p_cont[2][0] ) p_cont[2][i+w-1]+=p_c->H[i][j-i];
	  if((int) p_cont[3][0] ) p_cont[3][i+w-1]+=p_c->I[i][j-i];
	  if((int) p_cont[4][0] ) p_cont[4][i+w-1]+=p_c->M[i][j-i];
	}
      }
    }
  }
     
  if(p_cont != NULL){
    for (g=0; g<=4; g++){
      if((int) p_cont[g][0]) {
	fprintf(wastl, "# %s\n",unpaired[g]);
	for (i=1; i<=len; i++){
	 
	  fprintf(wastl,"%7d\t%7.10f\n",i,p_cont[g][i]);
	}
	fprintf(wastl,"&\n");
      }
    } 
  }

  if(pint != NULL){
    g=5;
    fprintf(wastl, "# %s\n",unpaired[g]);
    for (i=1; i<=len; i++){
      fprintf(wastl,"%3d\t%17.14f\n",i,pint[0][i]);
    }
    
    g=0;
    fprintf(wastl, "&\n# free energy of interaction\n");
    for (i=1; i<=len; i++){
      fprintf(wastl,"%3d\t%17.14f\n",i,pint[1][i]);
    }
  }
  fclose(wastl);
  if(p_c != NULL) {
    for(i=0; i < 5; i++){
      free(p_cont[i]);
    }
    free(p_cont);
  }
  return(1); /*success*/
}

