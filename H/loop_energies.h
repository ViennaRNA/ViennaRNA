#ifndef __VIENNA_RNA_PACKAGE_LOOP_ENERGIES_H__
#define __VIENNA_RNA_PACKAGE_LOOP_ENERGIES_H__

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <string.h>
#include "params.h"
#include "fold_vars.h"
#include "energy_par.h"

#ifdef __GNUC__
# define INLINE inline
#else
# define INLINE
#endif

/**
 *  \file loop_energies.h
 *  \brief Energy evaluation for MFE and partition function calculations
 * 
 *  <P>
 *  This file contains functions for the calculation of the free energy \f$\Delta G\f$
 *  of a hairpin- [ E_Hairpin() ] or interior-loop [ E_IntLoop()] .<BR>
 *  The unit of the free energy returned is \f$10^{-2} * \mathrm{kcal}/\mathrm{mol}\f$
 *  </P>
 *  <P>
 *  In case of computing the partition function, this file also supplies functions
 *  which return the Boltzmann weights \f$e^{-\Delta G/kT} \f$ for a hairpin- [ exp_E_Hairpin() ]
 *  or interior-loop [ exp_E_IntLoop() ].
 *  </P>
 */

/**
 *  \def E_MLstem(A,B,C,D)
 *  <H2>Compute the Energy contribution of a Multiloop stem</H2>
 *  This definition is a wrapper for the E_Stem() funtion.
 *  It is substituted by an E_Stem() funtion call with argument
 *  extLoop=0, so the energy contribution returned reflects a
 *  stem introduced in a multiloop.<BR>
 *  As for the parameters B (si1) and C (sj1) of the substituted
 *  E_Stem() function, you can inhibit to take 5'-, 3'-dangles
 *  or mismatch contributions to be taken into account by passing
 *  -1 to these parameters.
 * 
 *  \see    E_Stem()
 *  \param  A The pair type of the stem-closing pair
 *  \param  B The 5'-mismatching nucleotide
 *  \param  C The 3'-mismatching nucleotide
 *  \param  D The datastructure containing scaled energy parameters
 *  \return   The energy contribution of the introduced multiloop stem
 */
//#define E_MLstem(A,B,C,D)     E_Stem((A),(B),(C),0,(D))
INLINE  PRIVATE int E_MLstem(int type, int si1, int sj1, paramT *P);

/**
 *  \def exp_E_MLstem(A,B,C,D)
 *  This is the partition function variant of \ref E_MLstem()
 *  \see E_MLstem()
 *  \return The Boltzmann weighted energy contribution of the introduced multiloop stem
 */
//#define exp_E_MLstem(A,B,C,D) exp_E_Stem((A),(B),(C),0,(D))
INLINE  PRIVATE double exp_E_MLstem(int type, int si1, int sj1, pf_paramT *P);

/**
 *  \def E_ExtLoop(A,B,C,D)
 *  <H2>Compute the Energy contribution of an Exterior loop stem</H2>
 *  This definition is a wrapper for the E_Stem() funtion.
 *  It is substituted by an E_Stem() funtion call with argument
 *  extLoop=1, so the energy contribution returned reflects a
 *  stem introduced in an exterior-loop.<BR>
 *  As for the parameters B (si1) and C (sj1) of the substituted
 *  E_Stem() function, you can inhibit to take 5'-, 3'-dangles
 *  or mismatch contributions to be taken into account by passing
 *  -1 to these parameters.
 * 
 *  \see    E_Stem()
 *  \param  A The pair type of the stem-closing pair
 *  \param  B The 5'-mismatching nucleotide
 *  \param  C The 3'-mismatching nucleotide
 *  \param  D The datastructure containing scaled energy parameters
 *  \return   The energy contribution of the introduced exterior-loop stem
 */
//#define E_ExtLoop(A,B,C,D)      E_Stem((A),(B),(C),1,(D))
INLINE  PRIVATE int E_ExtLoop(int type, int si1, int sj1, paramT *P);

/**
 *  \def exp_E_ExtLoop(A,B,C,D)
 *  This is the partition function variant of \ref E_ExtLoop()
 *  \see E_ExtLoop()
 *  \return The Boltzmann weighted energy contribution of the introduced exterior-loop stem
 */
//#define exp_E_ExtLoop(A,B,C,D)  exp_E_Stem((A),(B),(C),1,(D))
INLINE  PRIVATE double exp_E_ExtLoop(int type, int si1, int sj1, pf_paramT *P);

/**
 *  <H2>Compute the Energy of an interior-loop</H2>
 *  This function computes the free energy \f$\Delta G\f$ of an interior-loop with the
 *  following structure: <BR>
 *  <PRE>
 *        3'  5'
 *        |   |
 *        U - V
 *    a_n       b_1
 *     .        .
 *     .        .
 *     .        .
 *    a_1       b_m
 *        X - Y
 *        |   |
 *        5'  3'
 *  </PRE>
 *  This general structure depicts an interior-loop that is closed by the base pair (X,Y).
 *  The enclosed base pair is (V,U) which leaves the unpaired bases a_1-a_n and b_1-b_n
 *  that constitute the loop. In this example, the length of the interior-loop is \f$(n+m)\f$
 *  where n or m may be 0 resulting in a bulge-loop or base pair stack.
 *  The mismatching nucleotides for the closing pair (X,Y) are:<BR>
 *  5'-mismatch: a_1<BR>
 *  3'-mismatch: b_m<BR>
 *  and for the enclosed base pair (V,U):<BR>
 *  5'-mismatch: b_1<BR>
 *  3'-mismatch: a_n<BR>
 *  \note Base pairs are always denoted in 5'->3' direction. Thus the enclosed base pair
 *  must be 'turned arround' when evaluating the free energy of the interior-loop
 *  \see scale_parameters()
 *  \see paramT
 *  \note This function is threadsafe
 * 
 *  \param  n1      The size of the 'left'-loop (number of unpaired nucleotides)
 *  \param  n2      The size of the 'right'-loop (number of unpaired nucleotides)
 *  \param  type    The pair type of the base pair closing the interior loop
 *  \param  type_2  The pair type of the enclosed base pair
 *  \param  si1     The 5'-mismatching nucleotide of the closing pair
 *  \param  sj1     The 3'-mismatching nucleotide of the closing pair
 *  \param  sp1     The 3'-mismatching nucleotide of the enclosed pair
 *  \param  sq1     The 5'-mismatching nucleotide of the enclosed pair
 *  \param  P       The datastructure containing scaled energy parameters
 *  \return The Free energy of the Interior-loop in dcal/mol
 */
INLINE  PRIVATE int E_IntLoop(int n1, int n2, int type, int type_2, int si1, int sj1, int sp1, int sq1, paramT *P);


/**
 *  <H2>Compute the Energy of a hairpin-loop</H2>
 *  To evaluate the free energy of a hairpin-loop, several parameters have to be known.
 *  A general hairpin-loop has this structure:<BR>
 *  <PRE>
 *        a3 a4
 *      a2     a5
 *      a1     a6
 *        X - Y
 *        |   |
 *        5'  3'
 *  </PRE>
 *  where X-Y marks the closing pair [e.g. a <B>(G,C)</B> pair]. The length of this loop is 6 as there are
 *  six unpaired nucleotides (a1-a6) enclosed by (X,Y). The 5' mismatching nucleotide is
 *  a1 while the 3' mismatch is a6. The nucleotide sequence of this loop is &quot;a1.a2.a3.a4.a5.a6&quot; <BR>
 *  \note The parameter sequence should contain the sequence of the loop in capital letters of the nucleic acid
 *  alphabet if the loop size is below 7. This is useful for unusually stable tri-, tetra- and hexa-loops
 *  which are treated differently (based on experimental data) if they are tabulated.
 *  @see scale_parameters()
 *  @see paramT
 *  \warning Not (really) thread safe! A threadsafe implementation will replace this function in a future release!\n
 *  Energy evaluation may change due to updates in global variable "tetra_loop"
 * 
 *  \param  size  The size of the loop (number of unpaired nucleotides)
 *  \param  type  The pair type of the base pair closing the hairpin
 *  \param  si1   The 5'-mismatching nucleotide
 *  \param  sj1   The 3'-mismatching nucleotide
 *  \param  string  The sequence of the loop
 *  \param  P     The datastructure containing scaled energy parameters
 *  \return The Free energy of the Hairpin-loop in dcal/mol
 */
INLINE  PRIVATE int E_Hairpin(int size, int type, int si1, int sj1, const char *string, paramT *P);

/**
 *  <H2>Compute the energy contribution of a stem branching off a loop-region</H2>
 *  This function computes the energy contribution of a stem that branches off
 *  a loop region. This can be the case in multiloops, when a stem branching off
 *  increases the degree of the loop but also <I>immediately interior base pairs</I>
 *  of an exterior loop contribute free energy.
 *  To switch the bahavior of the function according to the evaluation of a multiloop-
 *  or exterior-loop-stem, you pass the flag 'extLoop'.
 *  The returned energy contribution consists of a TerminalAU penalty if the pair type
 *  is greater than 2, dangling end contributions of mismatching nucleotides adjacent to
 *  the stem if only one of the si1, sj1 parameters is greater than 0 and mismatch energies
 *  if both mismatching nucleotides are positive values.
 *  Thus, to avoid incooperating dangling end or mismatch energies just pass a negative number,
 *  e.g. -1 to the mismatch argument.
 * 
 *  This is an illustration of how the energy contribution is assembled:
 *  <PRE>
 *        3'  5'
 *        |   |
 *        X - Y
 *  5'-si1     sj1-3'
 *  </PRE>
 * 
 *  Here, (X,Y) is the base pair that closes the stem that branches off a loop region.
 *  The nucleotides si1 and sj1 are the 5'- and 3'- mismatches, respectively. If the base pair
 *  type of (X,Y) is greater than 2 (i.e. an A-U or G-U pair, the TerminalAU penalty will be
 *  included in the energy contribution returned. If si1 and sj1 are both nonnegative numbers,
 *  mismatch energies will also be included. If one of sij or sj1 is a negtive value, only
 *  5' or 3' dangling end contributions are taken into account. To prohibit any of these mismatch
 *  contributions to be incoorporated, just pass a negative number to both, si1 and sj1.
 *  In case the argument extLoop is 0, the returned energy contribution also includes
 *  the <I>internal-loop-penalty</I> of a multiloop stem with closing pair type.
 * 
 *  \see    E_MLstem()
 *  \see    E_ExtLoop()
 *  \note   This function is threadsafe
 * 
 *  \param  type    The pair type of the first base pair un the stem
 *  \param  si1     The 5'-mismatching nucleotide
 *  \param  sj1     The 3'-mismatching nucleotide
 *  \param  extLoop A flag that indicates whether the contribution reflects the one of an exterior loop or not
 *  \param  P       The datastructure containing scaled energy parameters
 *  \return         The Free energy of the branch off the loop in dcal/mol
 * 
 */
INLINE  PRIVATE int E_Stem(int type, int si1, int sj1, int extLoop, paramT *P);

/**
 *  <H2>Compute the Boltzmann weighted energy contribution of a stem branching off a loop-region</H2>
 *  This is the partition function variant of \ref E_Stem()
 *  \see E_Stem()
 *  \note This function is threadsafe
 * 
 *  \return The Boltzmann weighted energy contribution of the branch off the loop
 */
INLINE  PRIVATE double exp_E_Stem(int type, int si1, int sj1, int extLoop, pf_paramT *P);

/**
 *  <H2>Compute Boltzmann weight \f$e^{-\Delta G/kT} \f$ of a hairpin loop</H2>
 *  multiply by scale[u+2]
 *  @see get_scaled_pf_parameters()
 *  @see pf_paramT
 *  @see E_Hairpin()
 *  \warning Not (really) thread safe! A threadsafe implementation will replace this function in a future release!\n
 *  Energy evaluation may change due to updates in global variable "tetra_loop"
 * 
 *  \param  u       The size of the loop (number of unpaired nucleotides)
 *  \param  type    The pair type of the base pair closing the hairpin
 *  \param  si1     The 5'-mismatching nucleotide
 *  \param  sj1     The 3'-mismatching nucleotide
 *  \param  string  The sequence of the loop
 *  \param  P       The datastructure containing scaled Boltzmann weights of the energy parameters
 *  \return The Boltzmann weight of the Hairpin-loop
 */
INLINE  PRIVATE double exp_E_Hairpin(int u, int type, short si1, short sj1, const char *string, pf_paramT *P);

/**
 *  <H2>Compute Boltzmann weight \f$e^{-\Delta G/kT} \f$ of interior loop</H2>
 *  multiply by scale[u1+u2+2] for scaling
 *  @see get_scaled_pf_parameters()
 *  @see pf_paramT
 *  @see E_IntLoop()
 *  \note This function is threadsafe
 * 
 *  \param  u1      The size of the 'left'-loop (number of unpaired nucleotides)
 *  \param  u2      The size of the 'right'-loop (number of unpaired nucleotides)
 *  \param  type    The pair type of the base pair closing the interior loop
 *  \param  type2   The pair type of the enclosed base pair
 *  \param  si1     The 5'-mismatching nucleotide of the closing pair
 *  \param  sj1     The 3'-mismatching nucleotide of the closing pair
 *  \param  sp1     The 3'-mismatching nucleotide of the enclosed pair
 *  \param  sq1     The 5'-mismatching nucleotide of the enclosed pair
 *  \param  P       The datastructure containing scaled Boltzmann weights of the energy parameters
 *  \return The Boltzmann weight of the Interior-loop
 */
INLINE  PRIVATE double  exp_E_IntLoop(int u1, int u2, int type, int type2, short si1, short sj1, short sp1, short sq1, pf_paramT *P);


/*
#################################
# BEGIN OF FUNCTION DEFINITIONS #
#################################
*/
INLINE  PRIVATE int E_Hairpin(int size, int type, int si1, int sj1, const char *string, paramT *P){
  int energy;

  energy = (size <= 30) ? P->hairpin[size] : P->hairpin[30]+(int)(P->lxc*log((size)/30.));
  if (tetra_loop){
    if (size == 4) { /* check for tetraloop bonus */
      char tl[7]={0}, *ts;
      strncpy(tl, string, 6);
      if ((ts=strstr(P->Tetraloops, tl)))
        return (P->Tetraloop_E[(ts - P->Tetraloops)/7]);
    }
  }
  {
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
      return (energy + (type>2 ? P->TerminalAU : 0));
    }
  }
  energy += P->mismatchH[type][si1][sj1];

  return energy;
}

INLINE  PRIVATE int E_IntLoop(int n1, int n2, int type, int type_2, int si1, int sj1, int sp1, int sq1, paramT *P){
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

INLINE  PRIVATE int E_Stem(int type, int si1, int sj1, int extLoop, paramT *P){
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

INLINE  PRIVATE int E_ExtLoop(int type, int si1, int sj1, paramT *P){
  int energy = 0;
  if(si1 >= 0 && sj1 >= 0){
    energy += P->mismatchExt[type][si1][sj1];
  }
  else if (si1 >= 0){
    energy += P->dangle5[type][si1];
  }
  else if (sj1 >= 0){
    energy += P->dangle3[type][sj1];
  }

  if(type > 2)
    energy += P->TerminalAU;

  return energy;
}

INLINE  PRIVATE int E_MLstem(int type, int si1, int sj1, paramT *P){
  int energy = 0;
  if(si1 >= 0 && sj1 >= 0){
    energy += P->mismatchM[type][si1][sj1];
  }
  else if (si1 >= 0){
    energy += P->dangle5[type][si1];
  }
  else if (sj1 >= 0){
    energy += P->dangle3[type][sj1];
  }

  if(type > 2)
    energy += P->TerminalAU;

  energy += P->MLintern[type];

  return energy;
}

INLINE  PRIVATE double exp_E_Hairpin(int u, int type, short si1, short sj1, const char *string, pf_paramT *P){
  double q, kT;
  kT = P->kT;   /* kT in cal/mol  */

  if(u <= 30)
    q = P->exphairpin[u];
  else
    q = P->exphairpin[30] * exp( -(P->lxc*log( u/30.))*10./kT);

  if(u < 3) return q; /* should only be the case when folding alignments */

  if ((tetra_loop)&&(u==4)) {
    char tl[7]={0,0,0,0,0,0,0}, *ts;
    strncpy(tl, string, 6);
    if ((ts=strstr(P->Tetraloops, tl))){
      if(type != 7)
        return (P->exptetra[(ts-P->Tetraloops)/7]);
      else
        q *= P->exptetra[(ts-P->Tetraloops)/7];
    }
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

INLINE  PRIVATE double exp_E_IntLoop(int u1, int u2, int type, int type2, short si1, short sj1, short sp1, short sq1, pf_paramT *P){
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
    else if (us==1) {
      if (ul==1){                    /* 1x1 loop */
        return P->expint11[type][type2][si1][sj1];
      }
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
  return z;
}

INLINE  PRIVATE double exp_E_Stem(int type, int si1, int sj1, int extLoop, pf_paramT *P){
  double energy = 1.0;
  double d5 = (si1 >= 0) ? P->expdangle5[type][si1] : 1.;
  double d3 = (sj1 >= 0) ? P->expdangle3[type][sj1] : 1.;

  if(type > 2)
    energy *= P->expTermAU;

  if(si1 >= 0 && sj1 >= 0)
    energy *= (extLoop) ? P->expmismatchExt[type][si1][sj1] : P->expmismatchM[type][si1][sj1];
  else
    energy *= d5 * d3;

  if(!extLoop) energy *= P->expMLintern[type];
  return energy;
}

INLINE  PRIVATE double exp_E_MLstem(int type, int si1, int sj1, pf_paramT *P){
  double energy = 1.0;
  if(si1 >= 0 && sj1 >= 0){
    energy *= P->expmismatchM[type][si1][sj1];
  }
  else if(si1 >= 0){
    energy *= P->expdangle5[type][si1];
  }
  else if(sj1 >= 0){
    energy *= P->expdangle3[type][sj1];
  }

  if(type > 2)
    energy *= P->expTermAU;

  energy *= P->expMLintern[type];
  return energy;
}

INLINE  PRIVATE double exp_E_ExtLoop(int type, int si1, int sj1, pf_paramT *P){
  double energy = 1.0;
  if(si1 >= 0 && sj1 >= 0){
    energy *= P->expmismatchExt[type][si1][sj1];
  }
  else if(si1 >= 0){
    energy *= P->expdangle5[type][si1];
  }
  else if(sj1 >= 0){
    energy *= P->expdangle3[type][sj1];
  }

  if(type > 2)
    energy *= P->expTermAU;

  return energy;
}

INLINE  PRIVATE int     E_IntLoop_Co(int type, int type_2, int i, int j, int p, int q, int cutpoint, short si1, short sj1, short sp1, short sq1, int dangles, paramT *P){
  int energy = 0;
  if(type > 2)   energy += P->TerminalAU;
  if(type_2 > 2) energy += P->TerminalAU;

  if(!dangles) return energy;

  int ci = (i>=cutpoint)||((i+1)<cutpoint);
  int cj = ((j-1)>=cutpoint)||(j<cutpoint);
  int cp = ((p-1)>=cutpoint)||(p<cutpoint);
  int cq = (q>=cutpoint)||((q+1)<cutpoint);

  int d3    = ci  ? P->dangle3[type][si1]   : 0;
  int d5    = cj  ? P->dangle5[type][sj1]   : 0;
  int d5_2  = cp  ? P->dangle5[type_2][sp1] : 0;
  int d3_2  = cq  ? P->dangle3[type_2][sq1] : 0;

  int tmm   = (cj && ci) ? P->mismatchExt[type][sj1][si1]   : d5 + d3;
  int tmm_2 = (cp && cq) ? P->mismatchExt[type_2][sp1][sq1] : d5_2 + d3_2;

  if(dangles == 2) return energy + tmm + tmm_2;

  /* now we may have non-double dangles only */
  if(i+2 < p){
    if(q+2 < j){ energy += tmm + tmm_2;}
    else if(q+2 == j){ energy += (cj && cq) ? MIN2(tmm + d5_2, tmm_2 + d3) : tmm + tmm_2;}
    else energy += d3 + d5_2;
  }
  else if(i+2 == p){
    if(q+2 < j){ energy += (ci && cp) ? MIN2(tmm + d3_2, tmm_2 + d5) : tmm + tmm_2;}
    else if(q+2 == j){
      energy += MIN2(tmm, MIN2(tmm_2, MIN2(d5 + d5_2, d3 + d3_2)));
    }
    else energy += MIN2(d3, d5_2);
  }
  else{
    if(q+2 < j){ energy += d5 + d3_2;}
    else if(q+2 == j){ energy += MIN2(d5, d3_2);}
  }
  return energy;
}

#endif
