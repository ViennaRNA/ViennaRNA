/* Last changed Time-stamp: <2008-03-26 17:05:48 ulim> */
/*                
		  partiton function for RNA secondary structures

		  Ivo L Hofacker

		  Vienna RNA package
*/
/*
  $Log: part_func_up.c,v $
  Revision 1.3  2008/05/08 14:11:55  ivo
  minor output changes

  Revision 1.2  2007/12/13 10:19:54  ivo
  major RNAup update from Ulli

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
#include <unistd.h>
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
static char rcsid[] UNUSED = "$Id: part_func_up.c,v 1.3 2008/05/08 14:11:55 ivo Exp $";
#define CO_TURN 0
#define MAX(x,y) (((x)>(y)) ? (x) : (y))
#define MIN(x,y) (((x)<(y)) ? (x) : (y))
#define ZERO(A) (fabs(A) < DBL_EPSILON)
#define EQUAL(A,B) (fabs((A)-(B)) < 1000*DBL_EPSILON)
#define PUBLIC
#define PRIVATE static
/* #define PF2_DEBUG 1 *//* tests for prop_unpaired */
/* #define NUMERIC 1 */

PUBLIC pu_contrib *pf_unstru(char *sequence, int w);
PUBLIC interact *pf_interact(const char *s1, const char *s2, pu_contrib *p_c, pu_contrib *p_c2, int w, char *cstruc, int incr3, int incr5);
/* free_pf_two: first argument output of pf_unstru() !!!! */ 
/* PUBLIC  void free_pf_two(pu_contrib *p_con, double **p_in); */
PUBLIC void free_pu_contrib(pu_contrib *p_con);
PUBLIC void free_interact(interact *pin);
PUBLIC int Up_plot(pu_contrib *p_c, pu_contrib *p_c_sh, interact *pint, int len, char *ofile, int w, char *select_contrib);
PUBLIC double expLoopEnergy(int u1, int u2, int type, int type2, short si1, short sj1, short sp1, short sq1);
PUBLIC double expHairpinEnergy(int u, int type, short si1, short sj1, const char *string);

typedef struct constrain { /* constrains for cofolding */
  int *indx;
  char *ptype;
} constrain;
PRIVATE void scale_stru_pf_params(unsigned int length);
PRIVATE void free_pf_unstru(void);
PRIVATE void  init_pf_two(int length);
/* allocate free space for pf_unpaired() and pf_up() */
PRIVATE void scale_int(const char *s, const char *sl,double *sc_int,int incr3, int incr5);
PRIVATE void get_unpaired(int length);
PRIVATE void free_unpaired(void);
PRIVATE void encode_seq(const char *s1, const char *s2);
PRIVATE constrain *get_ptypes(char *S, const char *structure);

PRIVATE pf_paramT *Pf = NULL;/* use this structure for all the exp-arrays*/
PRIVATE FLT_OR_DBL *qb, *qm, *prpr; /* add arrays for pf_unpaired()*/
PRIVATE FLT_OR_DBL *q1k, *qln;
PRIVATE double *qqm2, *qq_1m2, *qqm, *qqm1;

PRIVATE double *scale, *expMLbase;
PRIVATE char *ptype; /* precomputed array of pair types */ 
PRIVATE int init_length;  /* length in last call to init_pf_fold()*/
PRIVATE double init_temp; /* temperature in last call to scale_pf_params */
/* make iptypes array for intermolecular constrains (ipidx for indexing)*/



#define ISOLATED  256.0

/*-----------------------------------------------------------------*/
PRIVATE short *S, *S1, *SS, *SS2;
/* you have to call pf_fold(sequence, structure); befor pf_unstru */
PUBLIC pu_contrib *pf_unstru(char *sequence, int w)
{
  int n, i,j,v,k,l, ij, kl, u,u1,d,type, type_2, tt, uu;
  int o,p,po;
  unsigned int size, nu;
  double temp, tqm2, bla; 
  double qbt1, *tmp, Zuij, sum_l;
  double *store_H, *store_Io, **store_I2o; /* hairp., interior contribs */
  double *store_M_qm_o,*store_M_mlbase;/* multiloop contributions */
  pu_contrib *pu_test;
#ifdef PF2_DEBUG  
  pu_contrib *pu_stru;
  FLT_OR_DBL *puij;
#endif
  sum_l=0.0; 
  temp=0;
  n = (int) strlen(sequence);
  double sum_M[n+1];
  /* contributions to probability of being unpaired witihin a(n)
     H hairpin,
     I interior loop,
     M muliloop,
     E exterior loop*/
  size = sizeof(double *) * (n+1);
#ifdef PF2_DEBUG
  size = sizeof(FLT_OR_DBL) * ((n+1)*(n+2)/2); 
  puij = (FLT_OR_DBL *) space(size);
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
#endif

  pu_test = (pu_contrib *) space (sizeof(pu_contrib)*1);
  pu_test->length=n;
  /* pu_test->X[i][j] where i <= j and i [1...n], j = [1...w[ */ 
  pu_test->H = (double **) space(size);
  pu_test->I = (double **) space(size);
  pu_test->M = (double **) space(size);
  pu_test->E = (double **) space(size);
  
  for(i=0;i<=n;i++){
    pu_test->H[i] = (double *) space(sizeof(double)*(w+1));
    pu_test->I[i] = (double *) space(sizeof(double)*(w+1));
    pu_test->M[i]=  (double *) space(sizeof(double)*(w+1));
    pu_test->E[i] = (double *) space(sizeof(double)*(w+1));
    
  }
  init_pf_two(n);
   
  for (d=0; d<=TURN; d++) 
    for (i=1; i<=n-d; i++) {
      j=i+d;
      ij = iindx[i]-j;
#ifdef PF2_DEBUG      
      puij[ij]=0.0;
#endif      
      if(d < w) {
#ifdef PF2_DEBUG	
	pu_stru->H[i][d]=pu_stru->I[i][d]=pu_stru->M[i][d]=pu_stru->E[i][d]=0.;
#endif
	pu_test->H[i][d]=pu_test->I[i][d]=pu_test->M[i][d]=pu_test->E[i][d]=0.;
      }
    }
  nu = (unsigned int) n;
  size = ((nu+1)*(nu+2)/2);
  for (i=0; i<size; i++){
    prpr[i]= pr[i];
  }
  
  for (i=1; i<=n; i++)
    for (j=i+TURN+1; j<=n; j++) {
      ij = iindx[i]-j;
      /* i need the part_func of all structures outside bp[ij] */
      if(qb[ij] > 0.0) prpr[ij]= (pr[ij]/qb[ij]);
    }
  
  /* no do_backtrack test anymore ! */
  for (i=1; i<=n; i++)
    {
      /* set auxillary arrays to 0, reuse qqm and qqm1, reuse qqm2 and qq_1m2*/
      qqm[i]=qqm1[i]=0;
      qqm2[i]=qq_1m2[i]=0;
    }
  
  store_I2o = (double **) space(sizeof(double *)*(n+1)); /* for p,k */ 
  for(i=0;i<=n;i++) {
    store_I2o[i] = (double *) space(sizeof(double)*(MAXLOOP+2));
  }
  /* expMLbase[i-p]*dangles_po */
  size = sizeof(double) * ((n+1)*(n+2)/2);
  store_M_mlbase = (double *) space(sizeof(double)*(size+1)); 
  /* 2. exterior bp (p,o) encloses unpaired region [i,i+w[*/
  for (o=TURN+2;o<=n; o++) {
    /*allocate space for arrays to store different contributions to H, I & M */
    store_H = (double *) space(sizeof(double )*(o+2));
    /* unpaired between ]l,o[ */    
    store_Io = (double *) space(sizeof(double)*(o+2));
    /* qm[p+1,i-1]*dangles_po */
    store_M_qm_o = (double *) space(sizeof(double)*(n+1)); 
          
    for (p=o-TURN-1; p>=1; p--) {
      /* construction of partition function of segment [p,o], given that
	 an unpaired region [i,i+w[ exists within [p,o] */
      u = o-p-1; po = iindx[p]-o;
      type = ptype[po];
      if (type!=0) {  
	/*hairpin contribution*/
	if (((type==3)||(type==4))&&no_closingGU) {
	  temp = 0.;
	} else {
	  temp = (expHairpinEnergy(u, type, S1[p+1], S1[o-1], sequence+p-1) * prpr[po] * scale[u+2]); /* add scale[u+2] */
	}
	/* all H contribs are collect for the longest unpaired region */ 
	store_H[p+1]=temp;
#ifdef PF2_DEBUG
	for(i=p+1; i<=o-1;i++) {  
	  /* longest unpaired region */
	  for(j=i; j < MIN(i+w,o); j++) {
	    ij=iindx[i]-j;
	    puij[ij]+= temp;
	    pu_stru->H[i][j-i]+=temp;
	  }
	}
#endif	  


	/* interior loops with interior pair k,l and an unpaired region of
	 length w between p and k || l and o*/	
	for (k=p+1; k<=MIN(p+MAXLOOP+1,o-TURN-2); k++) {
	  u1 = k-p-1;	  
	  sum_l=0.0;
	  for (l=MAX(k+TURN+1,o-1-MAXLOOP+u1); l<o; l++) {
	    kl=iindx[k]-l;
	    type_2 = ptype[kl];
	    if(l+1 < o){
	      store_Io[l+1] += sum_l;
	    }
	    temp=0.;
	    if (type_2) {
	      type_2 = rtype[type_2];
	      /* add *scale[u1+u2+2] */ 
	      temp = qb[kl]  * (scale[u1+o-l+1] *
				expLoopEnergy(u1, o-l-1, type, type_2,
   				  S1[p+1], S1[o-1], S1[k-1], S1[l+1])) *
		prpr[po];
	      if(l+1 < o) {
		store_Io[l+1] += temp; /* unpaired region between ]l,o[ */
	      }
	      sum_l += temp;
#ifdef PF2_DEBUG	     
	      /* unpaired in region ]p,k[ p < i < i+w-1 < k*/
	      for (i=p+1; i <= k-1; i++ ){
		for (j=i; j < MIN(i+w,k); j++ ){
		  ij=iindx[i]-j; puij[ij]+=temp;
		  /* don't use this to test ]l,o[ contribution */
		  pu_stru->I[i][j-i]+=temp;
		}
	      }
	      /* unpaired in region ]l,o[, l < i < i+w-1 < o*/
	      for (i=l+1; i <= o-1; i++){
		for (j=i; j < MIN(i+w,o); j++ ){
		  ij=iindx[i]-j; puij[ij]+=temp;
		  /* don't use this to test ]p,k[ contribution */
		  pu_stru->I[i][j-i]+=temp;
		}
	      }
#endif
	    } /* end of if pair(k,l) */
	  } /* end of l */
	  /* unpaired in region ]p,k[  */
	  for(i=p+1;i <= k-1;i++) {
	    int max_v; 
	    max_v=MIN(w-1,k-i-1);
	    store_I2o[i][max_v] += sum_l;
	  }
	} /* end of k */
      } /*end of if(type) test for bp (p,o) */

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
	int p1i,pui;
	tqm2+=qm[iindx[p]-i]*qqm[i+1];
	if(  type !=0 ) {
	  pui= (p+1) < i ? iindx[p+1]-(i) : 0;
	  p1i= (p+1) < (i-1) ? iindx[p+1]-(i-1) : 0;
	  /*unpaired region expMLbase[i-p] left of structured
	    region qq_1m2[i+1]*/
	  /* @expMLbase:  note distance of i-p == i-(p+1)+1 */
	  store_M_mlbase[iindx[p+1]-i] += expMLbase[i-p]*temp*qq_1m2[i+1];
	  /*structured region qm[p1i] left of unpaired region */
	  /* contribition for unpaired region is added after the p-loop */
	  store_M_qm_o[i] += qm[p1i] * temp;

#ifdef PF2_DEBUG	  
	  for(j=i; j< MIN(i+w,o);j++) {	    
	    double temp2;
	    int iwo;
	    temp2=0;
	    ij=iindx[i]-j;
	    /* iindx[x] - y : x <= y */
	    p1i = (p+1) < (i-1) ? iindx[p+1]-(i-1) : 0;
	    iwo = (j+1) < (o-1) ? iindx[j+1]-(o-1) : 0; 

	    /* w unpaired bases left of structured region */
	    temp2 = expMLbase[j-p]*qq_1m2[j+1];
	    
	    /* w unpaired bases between structured region */
	    temp2 += qm[p1i]*expMLbase[j-i+1]*qm[iwo]; 
	    
	    /* multiply with revers dangles for prpr[po]... */
	    temp2 *= temp;
	    puij[ij]+=temp2;
	    pu_stru->M[i][j-i]+=temp2;
	  } /* end of j ...  */
#endif
	}
      }/*end of for i ... */
      /* qqm2[p] contrib with at least 2 loops in region (p,o) */ 
      qqm2[p]=tqm2;
    } /* end for (p=..) */
    double sum_h;
    sum_h = 0.0;
    for(i=1; i < o; i++) {
      int max_v,vo;      
      sum_h += store_H[i];
      max_v = MIN(w-1,o-i-1); 
      for(v=max_v; v >=0;v--) {
	/* Hairpins */	
	pu_test->H[i][v] += sum_h;/* store_H[i][v] + store_H[i][max_v]; */
	/* Interior loops: unpaired region between  ]l,o[ calculated here !*/
	/* unpaired region between ]p,k[ collected after after o-loop */
	if(v <= MIN(max_v,MAXLOOP)) {
	  pu_test->I[i][v] += store_Io[i]; /* ]l,o[ */
	}
	/* Multiloops:*/
	/* unpaired region [i,v] between structured regions ]p,i[ and ]v,o[. */
	/* store_M_qm_o[i] = part. funct over all structured regions ]p,i[ */
	vo = (i+v+1) <= (o-1) ? iindx[i+v+1]-(o-1): 0;
	pu_test->M[i][v] += store_M_qm_o[i]*expMLbase[v+1]*qm[vo];
      }      
    }
    tmp = qqm1; qqm1=qqm; qqm=tmp;
    tmp = qqm2; qqm2=qq_1m2; qq_1m2=tmp;
    
    free(store_Io);
    free(store_H);
    free(store_M_qm_o); 
  }/* end for (o=..) */
  
  for(i=0;i<=n;i++) { sum_M[i]=0.; }
  for(i=1; i < n; i++) {
    int max_v;
    double sum_iv;
    sum_iv=0.;
    max_v = MIN(w-1,n-i); 
    for(v=n; v >=0;v--) {      
      if(v <= MIN(max_v,MAXLOOP)) {
	/* all unpaired regions [i,v] between p and k in interior loops */
	/* notice v runs from max_v -> 0, sum_iv sums all int. l. contribs */
	/* for each x, v < x =< max_v, since they contribute to [i,v] */
	sum_iv += store_I2o[i][v];
	pu_test->I[i][v] += sum_iv;
      }
      /* all unpaired region [i,v] for a fixed v, given that */
      /* region ]v,o[ contains at least 2 structures qq_1m2[v+1]; */
      if(v >= i) {
	sum_M[v]+= store_M_mlbase[iindx[i]-v];
	if(v-i<=max_v) {
	  pu_test->M[i][v-i] += sum_M[v];
	}	
      }
    }
  }
  free(store_M_mlbase);
  for(i=0;i<=n;i++) {
    free(store_I2o[i]);
  }
  free(store_I2o);
  
  for (i=1; i<=n; i++) {
    /* set auxillary arrays to 0 */
    qqm[i]=qqm1[i]=0;
    qqm2[i]=qq_1m2[i]=0;
  }
    
  /* 2. exterior bp (p,o) encloses unpaired region [i,j]
     HERE index o goes from n...1 and index p from o+TURN+1 ... n,
     that is, we add the one multiloop contribution that we
     could not calculate before  */
  size = sizeof(double) * ((n+1)*(n+2)/2);
  store_M_mlbase = (double *) space(sizeof(double)*(size+1));
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
	
	if(type !=0) {
	  double temp2;
	  temp2=0;
	  /* structured region qq_1m2[i-1] left of unpaired r. expMLbase[p-i]*/
	  /* @expMLbase:  note distance of p-i == p+1-i+1 */
 	  store_M_mlbase[iindx[i]-p+1] +=  qq_1m2[i-1]*expMLbase[p-i]*temp;
#ifdef PF2_DEBUG 
	  /* w unpaired bases right of structured region */
	  temp2 = qq_1m2[i-1]*expMLbase[p-i];
	  /* multiply with revers dangles for prpr[po]*... */
	  temp2 *= temp;	  
	  for (j=i; j<MIN(i+w,p);j++) {
	    ij=iindx[i]-j;
	    puij[ij]+=temp2;
	    pu_stru->M[i][j-i]+=temp2;
	  }
#endif
	}
      }/*end of for i ....*/
      qqm2[p]=tqm2;
    }/* end for (p=..) */    
    tmp = qqm1; qqm1=qqm; qqm=tmp;
    tmp = qqm2; qqm2=qq_1m2; qq_1m2=tmp;
  }/* end for (o=..) */
  /* now collect the missing multiloop contributions */
  for(i=0;i<=n;i++) { sum_M[i]=0.; }
  for(i=1; i<=n;i++) {
    int v_max;
    v_max=MIN(w-1,n-i);
    for(v=n; v>=i;v--){
      sum_M[i]+=store_M_mlbase[iindx[i]-v];
      if ((v-i <= v_max) ) {
	pu_test->M[i][v-i] += sum_M[i];
      }	
    }
  }
  free(store_M_mlbase);
  
  /* 1. region [i,j] exterior to all loops */
  Zuij=0.;bla=0;
  for (i=1; i<=n; i++) {
    uu=0;
    for(j=i; j<MIN(i+w,n+1);j++){
      ij=iindx[i]-j;
      temp=q1k[i-1]*1*scale[j-i+1]*qln[j+1]/q1k[n];
#ifdef PF2_DEBUG      
      puij[ij]+=temp;
      pu_stru->E[i][j-i]+=temp;
      Zuij+=puij[ij];/* partition function over all contributions to puij*/
#endif      
      pu_test->E[i][j-i]+=temp;
      
    }
#ifdef PF2_DEBUG    
    bla+=puij[ij];
#endif
  }
  
#ifdef PF2_DEBUG
  /* different tests */
  printf("pf_unpaired Zuij=%.3f   Z=%.3f\n",Zuij,q1k[n]); 
  /* test if puij[iindx[i]-i] == 1 - p_paired von [i],
     where p_paired von [i] = sum_i pr[ij oder ji]*/
  double bla_pf;
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
    /* printf("%3d p_unstru[i]: %12.9f\t1-p_paired[i]: %12.9f\n",i,bla,(1-bla_pf)); */
  }

  
  /* check if pu_stru->H[ij]+pu_stru->I[ij]+pu_stru->M[ij]+pu_stru->E[ij] == puij[ij]  */
  for (i=1; i<n; i++) {
    for (j=i; j < MIN((i+w),n); j++) {
      int u_len;
      double sum_stru,sum_test,dG_us,dG_ut;
      /*get the free energy of opening region [i,j] - that is free energy necessary to remove all structures from region [i,j]*/
      ij = iindx[i]-j;
      sum_stru = pu_stru->H[i][j-i]+pu_stru->I[i][j-i]+pu_stru->M[i][j-i]+pu_stru->E[i][j-i];
      sum_test = pu_test->H[i][j-i]+pu_test->I[i][j-i]+pu_test->M[i][j-i]+pu_test->E[i][j-i];
      /* hairpin only */
      /* sum_stru = pu_stru->H[i][j-i]; */
/*       sum_test = pu_test->H[i][j-i]; */
      /* get free energy nessesary remove all structures from region [i,j]*/
      dG_us = -log(sum_stru)*(temperature+K0)*GASCONST/1000.0;
      dG_ut = -log(sum_test)*(temperature+K0)*GASCONST/1000.0;
      if(!EQUAL(dG_us,dG_ut)) {
	printf("i=%d, j=%d\ndG_us =%.18f\ndG_ut =%.18f\n",i,j,dG_us,dG_ut);
      }
#ifdef NUMERIC
      /* check different contributions to pr_unpaired */
      /* hairpin only */
      sum_stru = pu_stru->H[i][j-i];
      sum_test = pu_test->H[i][j-i];
      /* get free energy nessesary remove all structures from region [i,j]*/
      dG_us = -log(sum_stru)*(temperature+K0)*GASCONST/1000.0;
      dG_ut = -log(sum_test)*(temperature+K0)*GASCONST/1000.0;
      if(!EQUAL(dG_us,dG_ut)) {
	printf("hair: i=%d, j=%d\ndG_us =%.18f\ndG_ut =%.18f\n",i,j,dG_us,dG_ut);
      }
      /* intr only */
      sum_stru = pu_stru->I[i][j-i];
      sum_test = pu_test->I[i][j-i];
      
      /* get free energy nessesary remove all structures from region [i,j]*/
      dG_us = -log(sum_stru)*(temperature+K0)*GASCONST/1000.0;
      dG_ut = -log(sum_test)*(temperature+K0)*GASCONST/1000.0;
      if(!EQUAL(dG_us,dG_ut)) {
	printf("intr: i=%d, j=%d\ndG_us =%.18f\ndG_ut =%.18f\n",i,j,dG_us,dG_ut);
      }
      /* multi only */
      sum_stru = pu_stru->M[i][j-i];
      sum_test = pu_test->M[i][j-i];
      
      /* get free energy nessesary remove all structures from region [i,j]*/
      dG_us = -log(sum_stru)*(temperature+K0)*GASCONST/1000.0;
      dG_ut = -log(sum_test)*(temperature+K0)*GASCONST/1000.0;
      if(!EQUAL(dG_us,dG_ut)) {
	printf("mult: i=%d, j=%d\ndG_us =%.18f\ndG_ut =%.18f\n",i,j,dG_us,dG_ut);
      }
      
#endif
      /* print long double: %.18Lf */
    }
  }
  puij[0]=Zuij;
  free(puij);
  for(i=0;i<=n;i++){
    free(pu_stru->H[i]);
    free(pu_stru->I[i]);
    free(pu_stru->M[i]);
    free(pu_stru->E[i]);
  }
  free(pu_stru->H);
  free(pu_stru->I);
  free(pu_stru->M);
  free(pu_stru->E);
  free(pu_stru);
#endif
  free_pf_unstru();
  return pu_test;  
}
/*------------------------------------------------------------------------*/
/* s1 is the longer seq */
PUBLIC interact *pf_interact(const char *s1, const char *s2, pu_contrib *p_c, pu_contrib *p_c2, int w, char *cstruc, int incr3, int incr5)
{
  int i, j, k,l,n1,n2,add_i5,add_i3,i_max,k_max, pc_size;
  double temp, Z, rev_d, E, Z2,**p_c_S, **p_c2_S, int_scale;
  FLT_OR_DBL ****qint_4, **qint_ik;
  /* PRIVATE double **pint; array for pf_up() output */
  interact *Int;
  double G_min, G_is,Gi_min;
  int gi,gj,gk,gl,ci,cj,ck,cl,prev_k,prev_l;
  FLT_OR_DBL **int_ik;
  double Z_int, temp_int, temppfs;
  double const_scale,const_T;
  constrain *cc = NULL;  /* constrains for cofolding */
  char *Seq, *i_long,*i_short,*pos=NULL; /* short seq appended to long one */
  /* int ***pu_jl; */ /* positions of interaction in the short RNA */
  G_min = G_is = Gi_min = 100.0;
  gi = gj = gk = gl = ci = cj = ck = cl = 0;
  
  n1 = (int) strlen(s1); 
  n2 = (int) strlen(s2);
  prev_k = 1;
  prev_l = n2;
  i_long = (char *) space (sizeof(char)*(n1+1));
  i_short = (char *) space (sizeof(char)*(n2+1));
  /* fill structure constrain */
  Seq = (char *) space (sizeof(char)*(n1+n2+2));
  strcpy(Seq,s1);
  strcat(Seq,s2);
  cc = get_ptypes(Seq,cstruc);
  
  p_c_S = (double **) space (sizeof(double *)*(n1+1));
  
  for (i=1; i<=n1; i++) {
    pc_size = MIN((w+incr3+incr5),n1);
    /* printf("pc_size = %d\n",pc_size); */
    p_c_S[i] = (double *) space (sizeof(double)*(pc_size+1));    
    for (j=0; j < (pc_size); j++) {
      p_c_S[i][j] = p_c->H[i][j]+p_c->I[i][j]+p_c->M[i][j]+p_c->E[i][j];
      
    }
  }
  if(p_c2 != NULL) {
    p_c2_S = (double **) space (sizeof(double *)*(n2+1));
    for (i=1; i<=n2; i++) {
      pc_size = MIN(w,n2);
      p_c2_S[i] = (double *) space (sizeof(double)*(pc_size+2));    
      for (j=0; j < (pc_size); j++) {
	p_c2_S[i][j] = p_c2->H[i][j]+p_c2->I[i][j]+p_c2->M[i][j]+p_c2->E[i][j];
      
      }
    }
  }
  /*array for pf_up() output */
  Int = (interact *) space(sizeof(interact)*1);
  Int->Pi = (double *) space(sizeof(double)*(n1+2));
  Int->Gi = (double *) space(sizeof(double)*(n1+2));
   
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
/* int_ik */
  int_ik = (FLT_OR_DBL **) space(sizeof(FLT_OR_DBL *) * (n1+1));
  for (i=1; i<=n1; i++) {
    int_ik[i] = (FLT_OR_DBL *) space(sizeof(FLT_OR_DBL) * (n1+1));
  }
  Z_int=0.;
  /*  Gint = ( -log(int_ik[gk][gi])-( ((int) w/2)*log(pf_scale)) )*((Pf->temperature+K0)*GASCONST/1000.0); */
  const_scale = ((int) w/2)*log(pf_scale);
  const_T = ((Pf->temperature+K0)*GASCONST/1000.0); 
  encode_seq(s1, s2);
  /* static  short *S~S1, *S1~SS1, *SS~S2, *SS2; */
  for (i=0; i<=n1; i++) {
    Int->Pi[i]=Int->Gi[i]=0.;
  }
  E=0.;
  Z=0.;
  
  if ( fold_constrained && cstruc != NULL) {
    pos = strchr(cstruc,'|');
    if(pos) {
      ci=ck=cl=cj=0;
      /* long seq              & short seq
	 .........||..|||||....&....||||...  w = maximal interaction length
                 ck       ci       cj  cl    */
      strncpy(i_long,cstruc,n1);
      i_long[n1] = '\0';
      strncpy(i_short,&cstruc[n1],n2);
      i_short[n2] ='\0';
      pos = strchr(i_long,'|');
      if(pos) ck = (int) (pos-i_long)+1; /* k */
      pos = strrchr(i_long,'|');
      if(pos) ci = (int) (pos-i_long)+1; /* i */
      pos = strrchr(i_short,'|');
      if(pos) cl = (int) (pos-i_short)+1; /* l */
      pos = strchr(i_short,'|');
      if(pos) cj = (int) (pos-i_short)+1; /* j */

      if(ck > 0 && ci > 0 && ci-ck+1 > w) {
	fprintf(stderr, "distance between constrains in longer seq, %d, larger than -w = %d",ci-ck+1,w);
	nrerror("pf_interact: could not satisfy all constraints");
      }
      if(cj > 0 && cl > 0 && cl-cj+1 > w) {
	fprintf(stderr, "distance between constrains in shorter seq, %d, larger than -w = %d",cl-cj+1,w);
	nrerror("pf_interact: could not satisfy all constraints");
      }
    }
    
  } else if ( fold_constrained && cstruc == NULL) {
    nrerror("option -C selected, but no constrained structure given\n");
  }
  if(fold_constrained) pos = strchr(cstruc,'|');  
  
  /*  qint_4[i][j][k][l] contribution that region (k-i) in seq1 (l=n1)
      is paired to region (l-j) in seq 2(l=n2) that is
      a region closed by bp k-l  and bp i-j */
  qint_4 = (FLT_OR_DBL ****) space(sizeof(FLT_OR_DBL ***) * (n1+1));
  
  /* qint_4[i][j][k][l] */
  for (i=1; i<=n1; i++) {
    int end_k;
    end_k = i-w;
    if(fold_constrained && pos && ci) end_k= MAX(i-w, ci-w);
    /* '|' constrains for long sequence: index i from 1 to n1 (5' to 3')*/
    /* interaction has to include 3' most '|' constrain, ci */
    if(fold_constrained && pos && ci && i==1 && i<ci) 
      i= ci-w+1 > 1 ? ci-w+1 : 1; 
    /* interaction has to include 5' most '|' constrain, ck*/ 
    if(fold_constrained && pos && ck && i > ck+w-1) break;
    
    /* note: qint_4[i] will be freed before we allocate qint_4[i+1] */
    qint_4[i] = (FLT_OR_DBL ***) space(sizeof(FLT_OR_DBL **) * (n2+1));
    for (j=n2; j>0; j--) {
      qint_4[i][j] = (FLT_OR_DBL **) space(sizeof(FLT_OR_DBL*) * (w+1));
      for (k=0; k<=w; k++) {
	qint_4[i][j][k] = (FLT_OR_DBL *) space(sizeof(FLT_OR_DBL) * (w+1));
      }
    }
    
     prev_k=1;
    for (j=n2; j>0; j--) {
      int type, type2,end_l;
      end_l = j+w;
      if(fold_constrained && pos && ci) end_l= MIN(cj+w,j+w);
      /* '|' constrains for short sequence: index j from n2 to 1 (3' to 5')*/
      /* interaction has to include 5' most '|' constrain, cj */
      if(fold_constrained && pos && cj && j==n2 && j>cj) 
	j = cj+w-1 > n2 ? n2 : cj+w-1; 
      /* interaction has to include 3' most '|' constrain, cl*/ 
      if(fold_constrained && pos && cl && j < cl-w+1) break;
      type = cc->ptype[cc->indx[i]-(n1+j)];
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
      pc_size = MIN((w+incr3+incr5),n1);
      if(incr3 < pc_size) add_i3=incr3;
      else add_i3=pc_size-1;

      /* only one bp (no interior loop) */
      if(p_c2 == NULL) {/* consider only structure of longer seq. */
	qint_ik[i][i]+=qint_4[i][j][0][0]*rev_d*p_c_S[add_i5][add_i3]*scale[((int) w/2)];
	Z+=qint_4[i][j][0][0]*rev_d*p_c_S[add_i5][add_i3]*scale[((int) w/2)];
      } else {/* consider structures of both seqs. */
	qint_ik[i][i]+=qint_4[i][j][0][0]*rev_d*p_c_S[add_i5][add_i3]*p_c2_S[j][0]*scale[((int) w/2)];
	Z+=qint_4[i][j][0][0]*rev_d*p_c_S[add_i5][add_i3]*p_c2_S[j][0]*scale[((int) w/2)];
      }

/* int_ik */
      /* check deltaG_ges = deltaG_int + deltaG_unstr; */
      int_ik[i][i]+=qint_4[i][j][0][0]*rev_d*scale[((int) w/2)];
      Z_int+=qint_4[i][j][0][0]*rev_d*scale[((int) w/2)];
      temp_int=0.;
     
      temp=0.;
      prev_l = n2;
      for (k=i-1; k>end_k && k>0; k--) {      
	if (fold_constrained && pos && cstruc[k-1] == '|' && k > prev_k)
	  prev_k=k; 
	for (l=j+1; l< end_l && l<=n2; l++) {
	  int a,b,ia,ib,isw;
	  double scalew, tt, intt;
	  
	  type2 = cc->ptype[cc->indx[k]-(n1+l)];
	  /* '|' : l HAS TO be paired: not pair (k,x) where x>l allowed */
	  if(fold_constrained && pos && cstruc[n1+l-1] == '|' && l < prev_l)
	    prev_l=l; /*break*/
	  if(fold_constrained && pos && (k<=ck || i>=ci) && !type2) continue;
	  if(fold_constrained && pos && ((cstruc[k-1] == '|') || (cstruc[n1+l-1] == '|')) && !type2) break;
	  
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
	    if(k >= prev_k && l <= prev_l) { /* don't violate constrains */
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
	  }
	  /* '|' constrain in long sequence */ 
	  /* collect interactions starting before 5' most '|' constrain */
	  if ( fold_constrained && pos && ci && i < ci) continue;
	  /* collect interactions ending after 3' most '|' constrain*/
	  if ( fold_constrained && pos && ck &&  k > ck) continue;
	  /* '|' constrain in short sequence */
	  /* collect interactions starting before 5' most '|' constrain */
	  if ( fold_constrained && pos && cj && j > cj) continue;
	  /* collect interactions ending after 3' most '|' constrain*/
	  if ( fold_constrained && pos && cl && l < cl) continue;

	  /* scale everything to w/2*/
	  /* qint_ik[k][i] all interactions where k and i both are paired */
	  /* substract incr5 from k */
          if(k-incr5 > 0) add_i5=k-incr5;
          else add_i5=1;
	  /* add incr3 to i */
	  pc_size = MIN((w+incr3+incr5),n1);
	  if(i-k+incr3 < pc_size) add_i3=i-k+incr3;
	  else add_i3=pc_size-1;
	  
	  if(p_c2 == NULL) {/* consider only structure of longer seq. */
	    tt = qint_4[i][j][a][b]*p_c_S[add_i5][add_i3]*scalew*rev_d;
	  } else { /* consider structures of both seqs. */
	    tt = qint_4[i][j][a][b]*p_c_S[add_i5][add_i3]*p_c2_S[j][b]*scalew*rev_d;
	  }
	  temp+= tt;
	  qint_ik[k][i]+= tt;
	  /* int_ik */
	  /* check deltaG_ges = deltaG_int + deltaG_unstr; */
	  intt = qint_4[i][j][a][b]*scalew*rev_d;  
	  temp_int += intt;
	  int_ik[k][i]+= intt;
	  G_is = (-log(tt)-const_scale)*(const_T);
	  if (G_is < G_min || EQUAL(G_is,G_min)) {
	    G_min = G_is;	    
	    Gi_min =(-log(intt)-const_scale)*(const_T);
	    gi=i;
	    gj=j;
	    gk=k;
	    gl=l;
	  } 
	}
      }
      Z+=temp;
      /* int_ik */
      Z_int+=temp_int; 
    }
  
    /* free qint_4 values not needed any more */
    if(i > w) {
      int bla;
      bla=i-w;
      if (fold_constrained && pos && ci && i-w < ci-w+1) continue;
      if (fold_constrained && pos && ci) bla = MAX(ci-w+1,i-w);
      for (j=n2; j>0; j--) {
	for (k=0; k<=w; k++){
	  free(qint_4[bla][j][k]);
	}
	free(qint_4[bla][j]);
      }
      free(qint_4[bla]);
      qint_4[bla] = NULL;
    }
  }

  
  Z2=0.0;
  i_max = k_max = 0;
  for (i=1; i<=n1; i++) {
    for (k=i; k<=n1 && k<i+w; k++) {
      Z2+=qint_ik[i][k];
      for(l=i;l<=k;l++) {
	/* Int->Pi[l]: prob that position l is within a paired region */
	/* qint_ik[i][k] as well as Z are scaled to scale[((int) w/2) */
	Int->Pi[l]+=qint_ik[i][k]/Z;
	/* Int->Gi[l]: minimal delta G at position [l] */
	Int->Gi[l]=MIN(Int->Gi[l],
		       ( -log(qint_ik[i][k])-( ((int) w/2)*log(pf_scale)) )*
		       ((Pf->temperature+K0)*GASCONST/1000.0) );
      }
    }
  }  
  if(n1 > w){
    int start_i,end_i;
    start_i = n1-w+1;
    end_i=n1;
    if (fold_constrained && pos && ci) {
      /* a break in the k loop might result in unfreed values */
      start_i = ci-w+1 < n1-w+1 ? ci-w+1 : n1-w+1;
      start_i = start_i > 0 ? start_i : 1;
      /* start_i = ck; */
      end_i = ck+w-1 > n1 ? n1 : ck+w-1;
    }
    for (i=start_i; i<=end_i; i++) {
      if(qint_4[i] == NULL ) continue;      
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
    int start_i,end_i;
    start_i = 1;
    end_i=n1;
    if (fold_constrained && pos) {
      start_i = ci-w+1 > 0 ? ci-w+1 : 1;
      end_i = ck+w-1 > n1 ? n1 : ck+w-1;
    }
    
    for (i=start_i; i<=end_i; i++) {
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
  if(fold_constrained && (gi==0 || gk==0 ||  gl==0 || gj==0)) {
    nrerror("pf_interact: could not satisfy all constraints");    
  }
  /* fill structure interact */
  Int->length = n1;
  Int->i = gi;
  Int->j = gj;
  Int->k = gk;
  Int->l = gl;
  Int->Gikjl = G_min;
  Int->Gikjl_wo = Gi_min;
  
  free(i_long);
  free(i_short);
  
  for (i=1; i<=n1; i++) {
    free(int_ik[i]);
  }
  free(int_ik);
  for (i=1; i<=n1; i++) {
    free(qint_ik[i]);
  }
  free(qint_ik);  

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
  if(p_c2 != NULL) {
    for (i=1; i<=n2; i++) {
      free(p_c2_S[i]);      
    }
    free(p_c2_S);
  }
  free(Seq);
  free(cc->indx);
  free(cc->ptype);
  free(cc);
  return(Int);  
}
/*------------------------------------------------------------------------*/
/* use an extra scale for pf_interact, here sl is the longer sequence */
PRIVATE void scale_int(const char *s, const char *sl, double *sc_int,int incr3, int incr5)
{
  int n,nl,l_scales;
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
  qqm2 = (double *) space(sizeof(double)*(length+2));
  qq_1m2 = (double *) space(sizeof(double)*(length+2));
  qqm = (double *) space(sizeof(double)*(length+2));
  qqm1 = (double *) space(sizeof(double)*(length+2));
}
PRIVATE void free_unpaired(void)
{
  free(prpr);
  prpr=NULL;
  free(qqm2);
  free(qq_1m2);
  free(qqm);
  free(qqm1);
}
/* free all but the output structure for pf_unstru */
PRIVATE void free_pf_unstru(void) {
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


PUBLIC void free_pu_contrib(pu_contrib *p_con) {
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
  
  if(SS != NULL){
    free(SS);
    SS=NULL;
  }
  if(SS2 != NULL){
    free(SS2);
    SS2=NULL;
  }
}
PUBLIC void free_interact(interact *pin) {
  if(S != NULL && pin != NULL){
    free(S);
    S=NULL;
  }
  if(S1 != NULL && pin != NULL){
    free(S1);
    S1=NULL;
  }
  if(pin != NULL){
    free(pin->Pi);
    free(pin->Gi);
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
  if(s2 != NULL) {
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
}

/*-------------------------------------------------------------------------*/
 /* scale energy parameters and pre-calculate Boltzmann weights:
  most of this is done in structure Pf see params.c,h (function:
  scale_pf_parameters(), only arrays scale and expMLbase are handled here*/
PRIVATE void scale_stru_pf_params(unsigned int length)
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

PUBLIC int Up_plot(pu_contrib *p_c, pu_contrib *p_c_sh, interact *pint, int len, char *ofile, int w, char *select_contrib)
{
  /* produce output */
  FILE *wastl;
  int i,j,g,ws;
  double **p_cont,**p_cont_sh, dG_u; /* p_u AND its contributions */
  char *unpaired[6]={"total probability of being unpaired","probability of being unpaired in the exterior loop","probability of being unpaired in hairpin loops","probability of being unpaired in interior loops","probability of being unpaired in multiloops","interaction probability"};
  char *unpaired_fe[6]={"free energy to open all structures","free energy to open exterior loops","free energy to open hairpin loops","free energy to open interior loops","free energy to open multiloops","interaction probability"};

  if(p_c != NULL) {
    wastl = fopen(ofile,"a");
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
  p_cont_sh = NULL;
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
	fprintf(wastl, "# %s for u = %d\n",unpaired[g],w);
	for (i=1; i<=len; i++){
	  if(i >= w) {
	    fprintf(wastl,"%7d\t%7.15f\n",i,p_cont[g][i]);
	  }
	}
	fprintf(wastl,"&\n");
      }
    }
    double min_gu;
    int min_i,min_j;
    /* unpaired_out = (char*) space(sizeof(char)*(100)); */
    min_gu = 1000.0;
    min_i=min_j=0;
    dG_u=0.;
    for (g=0; g<=4; g++){
      if((int) p_cont[g][0]) {
	fprintf(wastl, "# %s for u = %d\n",unpaired_fe[g],w);
	for (i=1; i<=len; i++){
	  dG_u = -log(p_cont[g][i])*(temperature+K0)*GASCONST/1000.0;
	  if(i >= w) {
	    fprintf(wastl,"%7d\t%7.4f\n",i,dG_u);
	    if(g==0 &&( dG_u < min_gu )) {
	      min_gu = dG_u;
	      min_i=i-w+1;
	      min_j=i;
	    }
	  }
	}
	fprintf(wastl,"&\n");
      }
    }
    /* printf("\n%d,%d (%.3f) for u=%d\n",min_i,min_j,min_gu,w); */
  }
  /* get values for the shortet sequence */
  
  
  if(p_c_sh != NULL) {
    len = p_c_sh->length;
    ws=MIN(w,len);
      p_cont_sh = (double **) space(sizeof(double*)*(5));
    
    if(strchr(select_contrib, 'S')) {
      p_cont_sh[0] = (double *) space(sizeof(double)*(len+3));
      p_cont_sh[0][0] = 1.;
    } else {
      p_cont_sh[0] = (double *) space(sizeof(double)*(3));
      p_cont_sh[0][0] = 0.;
    }
    if(strchr(select_contrib, 'E')) {
      p_cont_sh[1] = (double *) space(sizeof(double)*(len+3));
      p_cont_sh[1][0] = 1.;
    } else {
      p_cont_sh[1] = (double *) space(sizeof(double)*(3));
      p_cont_sh[1][0] = 0.;
    }
    if(strchr(select_contrib, 'H')) {
      p_cont_sh[2] = (double *) space(sizeof(double)*(len+3));
      p_cont_sh[2][0] = 1.;
    } else {
      p_cont_sh[2] = (double *) space(sizeof(double)*(3));
      p_cont_sh[2][0] = 0.;
    }
    if(strchr(select_contrib, 'I')) {
      p_cont_sh[3] = (double *) space(sizeof(double)*(len+3));
      p_cont_sh[3][0] = 1.;
    } else {
      p_cont_sh[3] = (double *) space(sizeof(double)*(3));
      p_cont_sh[3][0] = 0.;
    }
    if(strchr(select_contrib, 'M')) {
      p_cont_sh[4] = (double *) space(sizeof(double)*(len+3));
      p_cont_sh[4][0] = 1.;
    } else {
      p_cont_sh[4] = (double *) space(sizeof(double)*(3));
      p_cont_sh[4][0] = 0.;
    }
    for (i=1; i<=len+2; i++) {
      /* not necessary for calculation, stores contributions to 
	 prob. upaired[i], i gives the START of the unpaired region*/
      if((int) p_cont_sh[0][0] ) p_cont_sh[0][i]=0.;
      if((int) p_cont_sh[1][0] ) p_cont_sh[1][i]=0.;
      if((int) p_cont_sh[2][0] ) p_cont_sh[2][i]=0.;
      if((int) p_cont_sh[3][0] ) p_cont_sh[3][i]=0.;
      if((int) p_cont_sh[4][0] ) p_cont_sh[4][i]=0.; 
    }
    
    for (i=1; i<=len; i++) {
      for (j=i; j < MIN((i+ws),len+1); j++) {
	double blubb;
	if( (j-i+1) == ws && i+ws-1 <= len) {
	  blubb = p_c_sh->H[i][j-i]+p_c_sh->I[i][j-i]+p_c_sh->M[i][j-i]+p_c_sh->E[i][j-i];
	  if((int) p_cont_sh[0][0] ) p_cont_sh[0][i+ws-1]+=blubb;
	  if((int) p_cont_sh[1][0] ) p_cont_sh[1][i+ws-1]+=p_c_sh->E[i][j-i];
	  if((int) p_cont_sh[2][0] ) p_cont_sh[2][i+ws-1]+=p_c_sh->H[i][j-i];
	  if((int) p_cont_sh[3][0] ) p_cont_sh[3][i+ws-1]+=p_c_sh->I[i][j-i];
	  if((int) p_cont_sh[4][0] ) p_cont_sh[4][i+ws-1]+=p_c_sh->M[i][j-i];
	}
      }
    }
  }
    
  if(p_cont_sh != NULL){
    len = p_c_sh->length;
    for (g=0; g<=4; g++){
      if((int) p_cont_sh[g][0]) {
	fprintf(wastl, "# %s for u = %d in shorter sequence\n",unpaired[g],w);
	for (i=1; i<=len; i++){
	  if(i >= ws) {
	    fprintf(wastl,"%7d\t%7.15f\n",i,p_cont_sh[g][i]);
	  }
	}
	fprintf(wastl,"&\n");
      }
    }
    dG_u=0.;
    for (g=0; g<=4; g++){
      if((int) p_cont_sh[g][0]) {
	fprintf(wastl, "# %s for u = %d in shorter sequence\n",unpaired_fe[g],w);
	for (i=1; i<=len; i++){
	  dG_u = -log(p_cont_sh[g][i])*(temperature+K0)*GASCONST/1000.0;
	  if(i >= ws) {
	    fprintf(wastl,"%7d\t%7.4f\n",i,dG_u);
	  }
	}
	fprintf(wastl,"&\n");
      }
    }
  }
  /* end of the same for the shorter sequence */
  len = p_c->length;
  if(pint != NULL){
    g=5;
    fprintf(wastl, "# %s\n",unpaired[g]);
    for (i=1; i<=len; i++){
      fprintf(wastl,"%3d\t%17.14f\n",i,pint->Pi[i]);
    }
    
    g=0;
    fprintf(wastl, "&\n# free energy of interaction\n");
    for (i=1; i<=len; i++){
      fprintf(wastl,"%3d\t%17.3f\n",i,pint->Gi[i]);
    }
  }
  fclose(wastl);
  if(p_c != NULL) {
    for(i=0; i < 5; i++){
      free(p_cont[i]);
    }
    free(p_cont);
  }
  if(p_c_sh != NULL) {
    for(i=0; i < 5; i++){
      free(p_cont_sh[i]);
    }
    free(p_cont_sh);
  }
  return(1); /*success*/
}
/*-------------------------------------------------------------------------*/
/* copy from part_func_co.c */
PRIVATE constrain *get_ptypes(char *Seq, const char *structure) {
  int n,i,j,k,l, length;
  constrain *con;
  
  length = strlen(Seq);
  make_pair_matrix();
  encode_seq(Seq,NULL);
  con = (constrain *) space(sizeof(constrain)*1);
  con->indx = (int *) space(sizeof(int)*(length+1));
  for (i=1; i<=length; i++) {
    con->indx[i] = ((length+1-i)*(length-i))/2 +length+1;    
  }
  con->ptype = (char *) space(sizeof(char)*((length+1)*(length+2)/2));
  
  n=S[0];
  for (k=1; k<=n-CO_TURN-1; k++)
    for (l=1; l<=2; l++) {
      int type,ntype=0,otype=0;
      i=k; j = i+CO_TURN+l; if (j>n) continue;
      type = pair[S[i]][S[j]];
      while ((i>=1)&&(j<=n)) {
	if ((i>1)&&(j<n)) ntype = pair[S[i-1]][S[j+1]];
	if (noLonelyPairs && (!otype) && (!ntype))
	  type = 0; /* i.j can only form isolated pairs */
	con->ptype[con->indx[i]-j] = (char) type;	
	otype =  type;
	type  = ntype;
	i--; j++;
      }
    }
  free(S);free(S1);

  if (fold_constrained&&(structure!=NULL)) {
    int hx, *stack;
    char type;
    stack = (int *) space(sizeof(int)*(n+1));
    for(hx=0, j=1; j<=n; j++) {
      switch (structure[j-1]) {
      case 'x': /* can't pair */
	for (l=1; l<j-CO_TURN; l++) con->ptype[con->indx[l]-j] = 0;
	for (l=j+CO_TURN+1; l<=n; l++) con->ptype[con->indx[j]-l] = 0;
	break;	
      case '(':
	stack[hx++]=j;
	/* fallthrough */
      case '<': /* pairs upstream */
	break;	
      case ')':
	if (hx<=0) {
	  fprintf(stderr, "%s\n", structure);
	  nrerror("1. unbalanced brackets in constraints");
	}
	i = stack[--hx];
	type = con->ptype[con->indx[i]-j];
	/* don't allow pairs i<k<j<l */
	for (k=i; k<=j; k++)
	  for (l=j; l<=n; l++) con->ptype[con->indx[k]-l] = 0;
	/* don't allow pairs k<i<l<j */
	for (k=1; k<=i; k++)
	  for (l=i; l<=j; l++) con->ptype[con->indx[k]-l] = 0;
	con->ptype[con->indx[i]-j] = (type==0)?7:type;
      case '>': /* pairs downstream */
	break;	
      }
    }
    if (hx!=0) {
      fprintf(stderr, "%s\n", structure);
      nrerror("2. unbalanced brackets in constraint string");
    }
    free(stack); 
  }
  return con;
}
