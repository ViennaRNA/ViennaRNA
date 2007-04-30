/* Last changed Time-stamp: <2006-08-12 13:57:16 ulim> */
/*                

		  c Ivo Hofacker

		  Vienna RNA package
*/
#include <config.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "energy_par.h"
#include "fold_vars.h"
#include "utils.h"
#include "params.h"

/*@unused@*/
static char rcsid[] UNUSED = "$Id: params.c,v 1.8 2007/04/30 15:06:03 ivo Exp $";

#define PUBLIC
#define PRIVATE static

#define MIN2(A, B)      ((A) < (B) ? (A) : (B))

PRIVATE paramT p;
PRIVATE int id=-1;
/* variables for partition function */
PRIVATE pf_paramT pf;
PRIVATE int pf_id=-1;

PUBLIC paramT *scale_parameters(void)
{
  unsigned int i,j,k,l;
  double tempf;
  /* if ((fabs(p.temperature - temperature)<1e-6)&&(id == p.id)) return &p; */

  tempf = ((temperature+K0)/Tmeasure);
  for (i=0; i<31; i++) 
    p.hairpin[i] = (int) hairpin37[i]*(tempf);
  for (i=0; i<=MIN2(30,MAXLOOP); i++) {
    p.bulge[i] = (int) bulge37[i]*tempf;
    p.internal_loop[i]= (int) internal_loop37[i]*tempf;
  }
  p.lxc = lxc37*tempf;
  for (; i<=MAXLOOP; i++) {
    p.bulge[i] = p.bulge[30]+(int)(p.lxc*log((double)(i)/30.));
    p.internal_loop[i] = p.internal_loop[30]+(int)(p.lxc*log((double)(i)/30.));
  }
  for (i=0; i<5; i++)
    p.F_ninio[i] = (int) F_ninio37[i]*tempf;
   
  for (i=0; (i*7)<strlen(Tetraloops); i++) 
    p.TETRA_ENERGY[i] = TETRA_ENTH37 - (TETRA_ENTH37-TETRA_ENERGY37[i])*tempf;
  for (i=0; (i*5)<strlen(Triloops); i++) 
    p.Triloop_E[i] =  Triloop_E37[i];
   
  p.MLbase = ML_BASE37*tempf;
  for (i=0; i<=NBPAIRS; i++) { /* includes AU penalty */
    p.MLintern[i] = ML_intern37*tempf;
    p.MLintern[i] +=  (i>2)?TerminalAU:0;
  }
  p.MLclosing = ML_closing37*tempf;

  p.TerminalAU = TerminalAU;
  
  p.DuplexInit = DuplexInit*tempf;

  /* stacks    G(T) = H - [H - G(T0)]*T/T0 */
  for (i=0; i<=NBPAIRS; i++)
    for (j=0; j<=NBPAIRS; j++)
      p.stack[i][j] = enthalpies[i][j] -
	(enthalpies[i][j] - stack37[i][j])*tempf;

  /* mismatches */
  for (i=0; i<=NBPAIRS; i++)
    for (j=0; j<5; j++)
      for (k=0; k<5; k++) {
	p.mismatchI[i][j][k] = mism_H[i][j][k] -
	  (mism_H[i][j][k] - mismatchI37[i][j][k])*tempf;
	p.mismatchH[i][j][k] = mism_H[i][j][k] -
	  (mism_H[i][j][k] - mismatchH37[i][j][k])*tempf;
	p.mismatchM[i][j][k] = mism_H[i][j][k] -
	  (mism_H[i][j][k] - mismatchM37[i][j][k])*tempf;
      }
   
  /* dangles */
  for (i=0; i<=NBPAIRS; i++)
    for (j=0; j<5; j++) {
      int dd;
      dd = dangle5_H[i][j] - (dangle5_H[i][j] - dangle5_37[i][j])*tempf; 
      p.dangle5[i][j] = (dd>0) ? 0 : dd;  /* must be <= 0 */
      dd = dangle3_H[i][j] - (dangle3_H[i][j] - dangle3_37[i][j])*tempf;
      p.dangle3[i][j] = (dd>0) ? 0 : dd;  /* must be <= 0 */
    }
  /* interior 1x1 loops */
  for (i=0; i<=NBPAIRS; i++)
    for (j=0; j<=NBPAIRS; j++)
      for (k=0; k<5; k++)
	for (l=0; l<5; l++) 
	  p.int11[i][j][k][l] = int11_H[i][j][k][l] -
	    (int11_H[i][j][k][l] - int11_37[i][j][k][l])*tempf;

  /* interior 2x1 loops */
  for (i=0; i<=NBPAIRS; i++)
    for (j=0; j<=NBPAIRS; j++)
      for (k=0; k<5; k++)
	for (l=0; l<5; l++) {
	  int m;
	  for (m=0; m<5; m++)
	    p.int21[i][j][k][l][m] = int21_H[i][j][k][l][m] -
	      (int21_H[i][j][k][l][m] - int21_37[i][j][k][l][m])*tempf;
	}
  /* interior 2x2 loops */
  for (i=0; i<=NBPAIRS; i++)
    for (j=0; j<=NBPAIRS; j++)
      for (k=0; k<5; k++)
	for (l=0; l<5; l++) {
	  int m,n;
	  for (m=0; m<5; m++)
	    for (n=0; n<5; n++)	     
	      p.int22[i][j][k][l][m][n] = int22_H[i][j][k][l][m][n] -
		(int22_H[i][j][k][l][m][n]-int22_37[i][j][k][l][m][n])*tempf;
	}

  strncpy(p.Tetraloops, Tetraloops, 1400);
  strncpy(p.Triloops, Triloops, 240);

  p.temperature = temperature;
  p.id = ++id;
  return &p;
}

PUBLIC paramT *copy_parameters(void) {
  paramT *copy;
  if (p.id != id) scale_parameters();
  
  copy = (paramT *) space(sizeof(paramT));
  memcpy(copy, &p, sizeof(paramT));
  return copy;
}

PUBLIC paramT *set_parameters(paramT *dest) {

  memcpy(&p, dest, sizeof(paramT));
  return &p;
}

/*------------------------------------------------------------------------*/
/* functions for partition function */

/* dangling ends should never be destabilizing, i.e. expdangle>=1         */
/* specific heat needs smooth function (2nd derivative)                   */
/* we use a*(sin(x+b)+1)^2, with a=2/(3*sqrt(3)), b=Pi/6-sqrt(3)/2,       */
/* in the interval b<x<sqrt(3)/2                                          */

#define SCALE 10
#define SMOOTH(X) ((X)/SCALE<-1.2283697)?0:(((X)/SCALE>0.8660254)?(X):\
          SCALE*0.38490018*(sin((X)/SCALE-0.34242663)+1)*(sin((X)/SCALE-0.34242663)+1))
/* #define SMOOTH(X) ((X)<0 ? 0 : (X)) */

PUBLIC pf_paramT *scale_pf_parameters(void)  {  
  /* scale energy parameters and pre-calculate Boltzmann weights */
  unsigned int i, j, k, l;
  double  kT, TT;
  double  GT;

  /* scale pf_params() in partfunc.c is only a wrapper, that calls
     this functions !! */

  pf.temperature = temperature;
  kT = (pf.temperature+K0)*GASCONST;   /* kT in cal/mol  */
  TT = (pf.temperature+K0)/(Tmeasure);

   /* loop energies: hairpins, bulges, interior, mulit-loops */
  for (i=0; i<31; i++) {
    GT =  hairpin37[i]*TT;
    pf.exphairpin[i] = exp( -GT*10./kT);
  }
  for (i=0; i<=MIN2(30, MAXLOOP); i++) {
    GT =  bulge37[i]*TT;
    pf.expbulge[i] = exp( -GT*10./kT);
    GT =  internal_loop37[i]*TT;
    pf.expinternal[i] = exp( -GT*10./kT);
  }
  /* special case of size 2 interior loops (single mismatch) */
  if (james_rule) pf.expinternal[2] = exp ( -80*10/kT);
   
  pf.lxc = lxc37*TT;
  
  GT =  DuplexInit*TT;
  pf.expDuplexInit = exp( -GT*10./kT);
  
  for (i=31; i<=MAXLOOP; i++) {
    GT = bulge37[30]*TT + (pf.lxc*log( i/30.));
    pf.expbulge[i] = exp( -GT*10./kT);
    GT = internal_loop37[30]*TT + (pf.lxc*log( i/30.));
    pf.expinternal[i] = exp( -GT*10./kT);
  }

  for (i=0; i<5; i++) {
    GT = F_ninio37[i]*TT;
    for (j=0; j<=MAXLOOP; j++)
      pf.expninio[i][j]=exp(-MIN2(MAX_NINIO,j*GT)*10/kT);
  }
  for (i=0; (i*7)<strlen(Tetraloops); i++) {
    GT = TETRA_ENTH37 - (TETRA_ENTH37-TETRA_ENERGY37[i])*TT;
    pf.exptetra[i] = exp( -GT*10./kT);
  }
  for (i=0; (i*5)<strlen(Triloops); i++) 
    pf.expTriloop[i] = exp(-Triloop_E37[i]*10/kT);

  GT =  ML_closing37*TT;
  pf.expMLclosing = exp( -GT*10/kT);

  for (i=0; i<=NBPAIRS; i++) { /* includes AU penalty */
    GT =  ML_intern37*TT;
    /* if (i>2) GT += TerminalAU; */
    pf.expMLintern[i] = exp( -GT*10./kT);
  }
  pf.expTermAU = exp(-TerminalAU*10/kT);

  GT =  ML_BASE37*TT;
  pf.expMLbase=(-10.*GT/kT); 
 
  /* if dangles==0 just set their energy to 0,
     don't let dangle energies become > 0 (at large temps),
     but make sure go smoothly to 0                        */
  for (i=0; i<=NBPAIRS; i++)
    for (j=0; j<=4; j++) {
      if (dangles) {
	GT = dangle5_H[i][j] - (dangle5_H[i][j] - dangle5_37[i][j])*TT;
	pf.expdangle5[i][j] = exp(SMOOTH(-GT)*10./kT);
	GT = dangle3_H[i][j] - (dangle3_H[i][j] - dangle3_37[i][j])*TT;
	pf.expdangle3[i][j] =  exp(SMOOTH(-GT)*10./kT);
      } else
	pf.expdangle3[i][j] = pf.expdangle5[i][j] = 1;
      if (i>2) /* add TermAU penalty into dangle3 */
	pf.expdangle3[i][j] *= pf.expTermAU;
    }

  /* stacking energies */
  for (i=0; i<=NBPAIRS; i++)
    for (j=0; j<=NBPAIRS; j++) {
      GT =  enthalpies[i][j] - (enthalpies[i][j] - stack37[i][j])*TT;
      pf.expstack[i][j] = exp( -GT*10/kT);
    }

  /* mismatch energies */
  for (i=0; i<=NBPAIRS; i++)
    for (j=0; j<5; j++)
      for (k=0; k<5; k++) {
	GT = mism_H[i][j][k] - (mism_H[i][j][k] - mismatchI37[i][j][k])*TT;
	pf.expmismatchI[i][j][k] = exp(-GT*10.0/kT);
	GT = mism_H[i][j][k] - (mism_H[i][j][k] - mismatchH37[i][j][k])*TT;
	pf.expmismatchH[i][j][k] = exp(-GT*10.0/kT);
	GT = mism_H[i][j][k] - (mism_H[i][j][k] - mismatchM37[i][j][k])*TT;
	pf.expmismatchM[i][j][k] = exp(-GT*10.0/kT);
      }
  
  
  /* interior lops of length 2 */
  for (i=0; i<=NBPAIRS; i++)
    for (j=0; j<=NBPAIRS; j++)
      for (k=0; k<5; k++)
	for (l=0; l<5; l++) {
	  GT = int11_H[i][j][k][l] -
	    (int11_H[i][j][k][l] - int11_37[i][j][k][l])*TT;
	  pf.expint11[i][j][k][l] = exp(-GT*10./kT);
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
	    pf.expint21[i][j][k][l][m] = exp(-GT*10./kT);
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
	      pf.expint22[i][j][k][l][m][n] = exp(-GT*10./kT);
	    }
	}

  strncpy(pf.Tetraloops, Tetraloops, 1400);
  strncpy(pf.Triloops, Triloops, 240);
  
  pf.id = ++pf_id;
  return &pf;
}


PUBLIC pf_paramT *copy_pf_param(void)   {
  pf_paramT *copy;
  if (pf.id != pf_id) scale_pf_parameters();
  
  copy = (pf_paramT *) space(sizeof(pf_paramT));
  memcpy(copy, &pf, sizeof(pf_paramT));
  return copy;
}


PUBLIC pf_paramT *set_pf_param(paramT *dest)  {
  memcpy(&pf, dest, sizeof(pf_paramT));
  return &pf;
}
