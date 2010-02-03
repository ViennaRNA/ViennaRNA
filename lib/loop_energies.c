/* Functions for Loop energy computations */


#ifdef __GNUC__
# define INLINE inline
#else
# define INLINE
#endif


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <string.h>
#include "fold_vars.h"
#include "energy_par.h"
#include "loop_energies.h"

/** 
*** \file loop_energies.c
*** <P>
*** This file contains functions for the calculation of the free energy \f$\Delta G\f$
*** of a hairpin- [ E_Hairpin() ] or interior-loop [ E_IntLoop()] .<BR>
*** The unit of the free energy returned is \f$10^{-2} * \mathrm{kcal}/\mathrm{mol}\f$
*** </P>
*** <P>
*** In case of computing the partition function, this file also supplies functions
*** which return the Boltzmann weights \f$e^{-\Delta G/kT} \f$ for a hairpin- [ exp_E_Hairpin() ]
*** or interior-loop [ exp_E_IntLoop() ].
*** </P>
***/


PUBLIC  INLINE  int E_Stem(int type, int si1, int sj1, int extLoop, paramT *P){
  int energy = 0;
  int d5 = (si1 >= 0) ? P->dangle5[type][si1] : 0;
  int d3 = (sj1 >= 0) ? P->dangle3[type][sj1] : 0;
  
  if(type > 2)
    energy += P->TerminalAU;

  if(si1 >= 0 && sj1 >= 0)
    energy += (extLoop) ? P->mismatchExt[type][si1][sj1] : P->mismatchM[type][si1][sj1];
  else
    energy += d5 + d3;

  if(!extLoop) energy += P->MLintern[type];
  return energy;
}



PUBLIC  INLINE  int E_Hairpin(int size, int type, int si1, int sj1, const char *string, paramT *P){
  int energy;
  
  energy = (size <= 30) ? P->hairpin[size] : P->hairpin[30]+(int)(P->lxc*log((size)/30.));
   
  if (tetra_loop){
    if (size == 4) { /* check for tetraloop bonus */
      char tl[7]={0}, *ts;
      strncpy(tl, string, 6);
      if ((ts=strstr(P->Tetraloops, tl)))
        return (P->Tetraloop_E[(ts - P->Tetraloops)/7]);
    }
    if (size == 6) {
      char tl[9]={0}, *ts;
      strncpy(tl, string, 8);
      if ((ts=strstr(P->Hexaloops, tl)))
        return (energy = P->Hexaloop_E[(ts - P->Hexaloops)/9]);
    }
    if (size == 3) {
      char tl[6]={0,0,0,0,0,0}, *ts;
      strncpy(tl, string, 5);
      if ((ts=strstr(P->Triloops, tl))) {
        return (P->Triloop_E[(ts - P->Triloops)/6]);
      }
      if (type>2)  /* neither CG nor GC */
        energy += P->TerminalAU; /* penalty for closing AU GU pair IVOO??
                                    sind dass jetzt beaunuesse oder mahlnuesse (vorzeichen?)*/
      return energy;
    }
  }
  energy += P->mismatchH[type][si1][sj1];

  return energy;
}

PUBLIC  INLINE  int E_IntLoop(int n1, int n2, int type, int type_2, int si1, int sj1, int sp1, int sq1, paramT *P){
  /* compute energy of degree 2 loop (stack bulge or interior) */
  int nl, ns, energy;

  if (n1>n2) { nl=n1; ns=n2;}
  else {nl=n2; ns=n1;}

  if (nl == 0)
    return P->stack[type][type_2];  /* stack */

  if (ns==0) {                      /* bulge */
    energy = (nl<=MAXLOOP)?P->bulge[nl]:
      (P->bulge[30]+(int)(P->lxc*log(nl/30.)));
    if (nl==1) energy += P->stack[type][type_2];
    else {
      if (type>2) energy += P->TerminalAU;
      if (type_2>2) energy += P->TerminalAU;
    }
    return energy;
  }
  else {                            /* interior loop */
    if (ns==1) {
      if (nl==1)                    /* 1x1 loop */
        return P->int11[type][type_2][si1][sj1];
      if (nl==2) {                  /* 2x1 loop */
        if (n1==1)
          energy = P->int21[type][type_2][si1][sq1][sj1];
        else
          energy = P->int21[type_2][type][sq1][si1][sp1];
        return energy;
      }
      else {  /* 1xn loop */
        energy = (nl+1<=MAXLOOP)?(P->internal_loop[nl+1]) : (P->internal_loop[30]+(int)(P->lxc*log((nl+1)/30.)));
        energy += MIN2(MAX_NINIO, (nl-ns)*P->ninio[2]);
        energy += P->mismatch1nI[type][si1][sj1] + P->mismatch1nI[type_2][sq1][sp1];
        return energy;
      }
    }
    else if (ns==2) {
      if(nl==2)      {              /* 2x2 loop */
        return P->int22[type][type_2][si1][sp1][sq1][sj1];}
      else if (nl==3){              /* 2x3 loop */
        energy = P->internal_loop[5]+P->ninio[2];
        energy += P->mismatch23I[type][si1][sj1] + P->mismatch23I[type_2][sq1][sp1];
        return energy;
      }
      
    }
    { /* generic interior loop (no else here!)*/
      energy = (n1+n2<=MAXLOOP)?(P->internal_loop[n1+n2]) : (P->internal_loop[30]+(int)(P->lxc*log((n1+n2)/30.)));

      energy += MIN2(MAX_NINIO, (nl-ns)*P->ninio[2]);

      energy += P->mismatchI[type][si1][sj1] + P->mismatchI[type_2][sq1][sp1];
    }
  }
  return energy;
}

#if 0
PUBLIC  INLINE  double exp_E_Stem(int type, int si1, int sj1, int extLoop, pf_paramT *P){
  double energy = 1.0;
  double d5 = (si1 >= 0) ? P->expdangle5[type][si1] : 1.;
  double d3 = (sj1 >= 0) ? P->expdangle3[type][sj1] : 1.;
  
  if(type > 2)
    energy *= P->expTermAU;

  if(si1 >= 0 && sj1 >= 0 && (extLoop))
    energy *= (extLoop) ? P->expmismatchExt[type][si1][sj1] : P->expmismatchM[type][si1][sj1];
  else
    energy *= d5 * d3;

  if(!extLoop) energy *= P->expMLintern[type];
  return energy;
}
#endif

/* compute Boltzmann weight of a hairpin loop, multiply by scale[u+2] */
PUBLIC  INLINE  double exp_E_Hairpin(int u, int type, short si1, short sj1, const char *string, pf_paramT *P){
  double q, kT;
  kT = (P->temperature+K0)*GASCONST;   /* kT in cal/mol  */
  if(u <= 30)
    q = P->exphairpin[u];
  else
    q = P->exphairpin[30] * exp( -(P->lxc*log( u/30.))*10./kT);

  if ((tetra_loop)&&(u==4)) {
    char tl[7]={0,0,0,0,0,0,0}, *ts;
    strncpy(tl, string, 6);
    if ((ts=strstr(P->Tetraloops, tl)))
      return (P->exptetra[(ts-P->Tetraloops)/7]);
  }
  if ((tetra_loop)&&(u==6)) {
    char tl[9]={0,0,0,0,0,0,0,0,0}, *ts;
    strncpy(tl, string, 6);
    if ((ts=strstr(P->Hexaloops, tl)))
      return  (P->exphex[(ts-P->Hexaloops)/9]);
  }
  if (u==3) {
    char tl[6]={0,0,0,0,0,0}, *ts;
    strncpy(tl, string, 5);
    if ((ts=strstr(P->Triloops, tl)))
      return (P->exptri[(ts-P->Triloops)/6]);
    if (type>2)
      q *= P->expTermAU;
  }
  else /* no mismatches for tri-loops */
    q *= P->expmismatchH[type][si1][sj1];

  return q;
}

#if 0
/* compute Boltzmann weight of interior loop, multiply by scale[u1+u2+2] for scaling */
PUBLIC  INLINE  double exp_E_IntLoop(int u1, int u2, int type, int type2, short si1, short sj1, short sp1, short sq1, pf_paramT *P){
#if 0
  double z=0;
  int no_close = 0;

  if ((no_closingGU) && ((type2==3)||(type2==4)||(type==2)||(type==4)))
    no_close = 1;

  if ((u1==0) && (u2==0)) /* stack */
    z = P->expstack[type][type2];
  else if (no_close==0) {
    if ((u1==0)||(u2==0)) { /* bulge */
      int u;
      u = (u1==0)?u2:u1;
      z = P->expbulge[u];
      if (u2+u1==1) z *= P->expstack[type][type2];
      else {
        if (type>2) z *= P->expTermAU;
        if (type2>2) z *= P->expTermAU;
      }
    }
    else {     /* interior loop */
      if (u1+u2==2) /* size 2 is special */
        z = P->expint11[type][type2][si1][sj1];
      else if ((u1==1) && (u2==2))
        z = P->expint21[type][type2][si1][sq1][sj1];
      else if ((u1==2) && (u2==1))
        z = P->expint21[type2][type][sq1][si1][sp1];
      else if ((u1==2) && (u2==2))
        z = P->expint22[type][type2][si1][sp1][sq1][sj1];
      else if (((u1==2)&&(u2==3))||((u1==3)&&(u2==2))){ /*2-3 is special*/
        z = P->expinternal[5]*P->expmismatch23I[type][si1][sj1]*P->expmismatch23I[type2][sq1][sp1];
        z *= P->expninio[2][1];
      }
      else if ((u1==1)||(u2==1)) {  /*1-n is special*/
        z = P->expinternal[u1+u2] * P->expmismatch1nI[type][si1][sj1] * P->expmismatch1nI[type2][sq1][sp1];
        z *= P->expninio[2][abs(u1-u2)];
      }
      else {
        z = P->expinternal[u1+u2] * P->expmismatchI[type][si1][sj1] * P->expmismatchI[type2][sq1][sp1];
        z *= P->expninio[2][abs(u1-u2)];
      }
    }
  }
#else
  int ul, us, no_close = 0;
  double z;

  if ((no_closingGU) && ((type2==3)||(type2==4)||(type==2)||(type==4)))
    no_close = 1;

  if (u1>u2) { ul=u1; us=u2;}
  else {ul=u2; us=u1;}

  if (ul==0) /* stack */
    z = P->expstack[type][type2];
  else if(!no_close){
    if (us==0) {                      /* bulge */
      z = P->expbulge[ul];
      if (ul==1) z *= P->expstack[type][type2];
      else {
        if (type>2) z *= P->expTermAU;
        if (type2>2) z *= P->expTermAU;
      }
      return z;
    }
  }
  else {                            /* interior loop */
    if (us==1) {
      if (ul==1)                    /* 1x1 loop */
        return P->expint11[type][type2][si1][sj1];
      if (ul==2) {                  /* 2x1 loop */
        if (u1==1)
          return P->expint21[type][type2][si1][sq1][sj1];
        else
          return P->expint21[type2][type][sq1][si1][sp1];
      }
      else {  /* 1xn loop */
        z = P->expinternal[ul+us] * P->expmismatch1nI[type][si1][sj1] * P->expmismatch1nI[type2][sq1][sp1];
        return z * P->expninio[2][ul-us];
      }
    }
    else if (us==2) {
      if(ul==2) /* 2x2 loop */
        return P->expint22[type][type2][si1][sp1][sq1][sj1];
      else if(ul==3){              /* 2x3 loop */
        z = P->expinternal[5]*P->expmismatch23I[type][si1][sj1]*P->expmismatch23I[type2][sq1][sp1];
        return z * P->expninio[2][1];
      }
    }
    /* generic interior loop (no else here!)*/
    z = P->expinternal[ul+us] * P->expmismatchI[type][si1][sj1] * P->expmismatchI[type2][sq1][sp1];
    return z * P->expninio[2][ul-us];
    
  }
#endif
  return z;
}
#endif
