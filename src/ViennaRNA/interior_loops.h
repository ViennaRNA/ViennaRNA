#ifndef VIENNA_RNA_PACKAGE_INTERIOR_LOOPS_H
#define VIENNA_RNA_PACKAGE_INTERIOR_LOOPS_H

#include <ViennaRNA/utils.h>
#include "ViennaRNA/energy_par.h"
#include <ViennaRNA/data_structures.h>
#include <ViennaRNA/params.h>
#include <ViennaRNA/constraints.h>

#ifdef __GNUC__
# define INLINE inline
#else
# define INLINE
#endif

#ifdef ON_SAME_STRAND
#undef ON_SAME_STRAND
#endif

#define ON_SAME_STRAND(I,J,C)  (((I)>=(C))||((J)<(C)))

/**
 *  @file     interior_loops.h
 *  @ingroup  loops
 *  @brief    Energy evaluation of interior loops for MFE and partition function calculations
 */

/**
 *  @{
 *  @ingroup   loops
 */

/**
 *  <H2>Compute the Energy of an interior-loop</H2>
 *  This function computes the free energy @f$\Delta G@f$ of an interior-loop with the
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
 *  that constitute the loop. In this example, the length of the interior-loop is @f$(n+m)@f$
 *  where n or m may be 0 resulting in a bulge-loop or base pair stack.
 *  The mismatching nucleotides for the closing pair (X,Y) are:<BR>
 *  5'-mismatch: a_1<BR>
 *  3'-mismatch: b_m<BR>
 *  and for the enclosed base pair (V,U):<BR>
 *  5'-mismatch: b_1<BR>
 *  3'-mismatch: a_n<BR>
 *  @note Base pairs are always denoted in 5'->3' direction. Thus the enclosed base pair
 *  must be 'turned arround' when evaluating the free energy of the interior-loop
 *  @see scale_parameters()
 *  @see vrna_param_t
 *  @note This function is threadsafe
 * 
 *  @param  n1      The size of the 'left'-loop (number of unpaired nucleotides)
 *  @param  n2      The size of the 'right'-loop (number of unpaired nucleotides)
 *  @param  type    The pair type of the base pair closing the interior loop
 *  @param  type_2  The pair type of the enclosed base pair
 *  @param  si1     The 5'-mismatching nucleotide of the closing pair
 *  @param  sj1     The 3'-mismatching nucleotide of the closing pair
 *  @param  sp1     The 3'-mismatching nucleotide of the enclosed pair
 *  @param  sq1     The 5'-mismatching nucleotide of the enclosed pair
 *  @param  P       The datastructure containing scaled energy parameters
 *  @return The Free energy of the Interior-loop in dcal/mol
 */
PRIVATE INLINE int E_IntLoop(int n1,
                              int n2,
                              int type,
                              int type_2,
                              int si1,
                              int sj1,
                              int sp1,
                              int sq1,
                              vrna_param_t *P);

/**
 *  <H2>Compute Boltzmann weight @f$e^{-\Delta G/kT} @f$ of interior loop</H2>
 *  multiply by scale[u1+u2+2] for scaling
 *  @see get_scaled_pf_parameters()
 *  @see vrna_exp_param_t
 *  @see E_IntLoop()
 *  @note This function is threadsafe
 * 
 *  @param  u1      The size of the 'left'-loop (number of unpaired nucleotides)
 *  @param  u2      The size of the 'right'-loop (number of unpaired nucleotides)
 *  @param  type    The pair type of the base pair closing the interior loop
 *  @param  type2   The pair type of the enclosed base pair
 *  @param  si1     The 5'-mismatching nucleotide of the closing pair
 *  @param  sj1     The 3'-mismatching nucleotide of the closing pair
 *  @param  sp1     The 3'-mismatching nucleotide of the enclosed pair
 *  @param  sq1     The 5'-mismatching nucleotide of the enclosed pair
 *  @param  P       The datastructure containing scaled Boltzmann weights of the energy parameters
 *  @return The Boltzmann weight of the Interior-loop
 */
PRIVATE INLINE FLT_OR_DBL exp_E_IntLoop(int u1,
                                        int u2,
                                        int type,
                                        int type2,
                                        short si1,
                                        short sj1,
                                        short sp1,
                                        short sq1,
                                        vrna_exp_param_t *P);


PRIVATE INLINE int E_IntLoop_Co(int type,
                                int type_2,
                                int i,
                                int j,
                                int p,
                                int q,
                                int cutpoint,
                                short si1,
                                short sj1,
                                short sp1,
                                short sq1,
                                int dangles,
                                vrna_param_t *P);


/*
#################################
# BEGIN OF FUNCTION DEFINITIONS #
#################################
*/

/*
 *  ugly but fast interior loop evaluation
 *
 *  Avoid including this function in your own code. It only serves
 *  as a fast inline block internally re-used throughout the RNAlib. It
 *  evalutes the free energy of interior loops in single sequences or sequence
 *  hybrids. Soft constraints are also applied if available.
 *
 *  NOTE: do not include into doxygen reference manual!
 */
PRIVATE INLINE int
ubf_eval_int_loop(  int i,
                    int j,
                    int p,
                    int q,
                    int i1,
                    int j1,
                    int p1,
                    int q1,
                    short si,
                    short sj,
                    short sp,
                    short sq,
                    unsigned char type,
                    unsigned char type_2,
                    int *rtype,
                    int ij,
                    int cp,
                    vrna_param_t *P,
                    vrna_sc_t *sc){

  int energy, u1, u2;

  u1 = p1 - i;
  u2 = j1 - q;

  if((cp < 0) || (ON_SAME_STRAND(i, p, cp) && ON_SAME_STRAND(q, j, cp))){ /* regular interior loop */
    energy = E_IntLoop(u1, u2, type, type_2, si, sj, sp, sq, P);
  } else { /* interior loop like cofold structure */
    short Si, Sj;
    Si  = ON_SAME_STRAND(i, i1, cp) ? si : -1;
    Sj  = ON_SAME_STRAND(j1, j, cp) ? sj : -1;
    energy = E_IntLoop_Co(rtype[type], rtype[type_2],
                            i, j, p, q,
                            cp,
                            Si, Sj,
                            sp, sq,
                            P->model_details.dangles,
                            P);
  }

  /* add soft constraints */
  if(sc){
    if(sc->energy_up)
      energy += sc->energy_up[i1][u1]
                + sc->energy_up[q1][u2];

    if(sc->energy_bp)
      energy += sc->energy_bp[ij];

    if(sc->energy_stack)
      if(u1 + u2 == 0){
        int a =   sc->energy_stack[i]
                  + sc->energy_stack[p]
                  + sc->energy_stack[q]
                  + sc->energy_stack[j];
        energy += a;
      }
    if(sc->f)
      energy += sc->f(i, j, p, q, VRNA_DECOMP_PAIR_IL, sc->data);
  }

  return energy;

}

/*
 *  ugly but fast exterior interior loop evaluation
 *
 *  Avoid including this function in your own code. It only serves
 *  as a fast inline block internally re-used throughout the RNAlib. It
 *  evalutes the free energy of interior loops in single sequences or sequence
 *  hybrids. Soft constraints are also applied if available.
 *
 *  NOTE: do not include into doxygen reference manual!
 */
PRIVATE INLINE int
ubf_eval_ext_int_loop(int i,
                      int j,
                      int p,
                      int q,
                      int i1,
                      int j1,
                      int p1,
                      int q1,
                      short si,
                      short sj,
                      short sp,
                      short sq,
                      unsigned char type,
                      unsigned char type_2,
                      int length,
                      vrna_param_t *P,
                      vrna_sc_t *sc){

  int energy, u1, u2, u3;
  
  u1 = i1;
  u2 = p1 - j;
  u3 = length - q;

  energy = E_IntLoop(u2, u1 + u3, type, type_2, si, sj, sp, sq, P);

  /* add soft constraints */
  if(sc){
    if(sc->energy_up){
      energy += sc->energy_up[j1][u2]
                + ((u3 > 0) ? sc->energy_up[q1][u3] : 0)
                + ((u1 > 0) ? sc->energy_up[1][u1] : 0);
    }
    if(sc->energy_stack)
      if(u1 + u2 + u3 == 0)
        energy +=   sc->energy_stack[i]
                  + sc->energy_stack[p]
                  + sc->energy_stack[q]
                  + sc->energy_stack[j];

    if(sc->f)
      energy += sc->f(i, j, p, q, VRNA_DECOMP_PAIR_IL, sc->data);
  }

  return energy;

}

PRIVATE INLINE int
E_IntLoop(int n1,
          int n2,
          int type,
          int type_2,
          int si1,
          int sj1,
          int sp1,
          int sq1,
          vrna_param_t *P){

  /* compute energy of degree 2 loop (stack bulge or interior) */
  int nl, ns, u, energy;
  energy = INF;

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
      u = nl + ns;
      energy = (u <= MAXLOOP) ? (P->internal_loop[u]) : (P->internal_loop[30]+(int)(P->lxc*log((u)/30.)));

      energy += MIN2(MAX_NINIO, (nl-ns)*P->ninio[2]);

      energy += P->mismatchI[type][si1][sj1] + P->mismatchI[type_2][sq1][sp1];
    }
  }
  return energy;
}

PRIVATE INLINE FLT_OR_DBL
exp_E_IntLoop(int u1,
              int u2,
              int type,
              int type2,
              short si1,
              short sj1,
              short sp1,
              short sq1,
              vrna_exp_param_t *P){

  int ul, us, no_close = 0;
  double z = 0.;
  int noGUclosure = P->model_details.noGUclosure;

  if ((noGUclosure) && ((type2==3)||(type2==4)||(type==3)||(type==4)))
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
      return (FLT_OR_DBL)z;
    }
    else if (us==1) {
      if (ul==1){                    /* 1x1 loop */
        return (FLT_OR_DBL)(P->expint11[type][type2][si1][sj1]);
      }
      if (ul==2) {                  /* 2x1 loop */
        if (u1==1)
          return (FLT_OR_DBL)(P->expint21[type][type2][si1][sq1][sj1]);
        else
          return (FLT_OR_DBL)(P->expint21[type2][type][sq1][si1][sp1]);
      }
      else {  /* 1xn loop */
        z = P->expinternal[ul+us] * P->expmismatch1nI[type][si1][sj1] * P->expmismatch1nI[type2][sq1][sp1];
        return (FLT_OR_DBL)(z * P->expninio[2][ul-us]);
      }
    }
    else if (us==2) {
      if(ul==2) /* 2x2 loop */
        return (FLT_OR_DBL)(P->expint22[type][type2][si1][sp1][sq1][sj1]);
      else if(ul==3){              /* 2x3 loop */
        z = P->expinternal[5]*P->expmismatch23I[type][si1][sj1]*P->expmismatch23I[type2][sq1][sp1];
        return (FLT_OR_DBL)(z * P->expninio[2][1]);
      }
    }
    /* generic interior loop (no else here!)*/
    z = P->expinternal[ul+us] * P->expmismatchI[type][si1][sj1] * P->expmismatchI[type2][sq1][sp1];
    return (FLT_OR_DBL)(z * P->expninio[2][ul-us]);

  }
  return (FLT_OR_DBL)z;
}

PRIVATE INLINE int
E_IntLoop_Co( int type,
              int type_2,
              int i,
              int j,
              int p,
              int q,
              int cutpoint,
              short si1,
              short sj1,
              short sp1,
              short sq1,
              int dangles,
              vrna_param_t *P){

  int energy, ci, cj, cp, cq, d3, d5, d5_2, d3_2, tmm, tmm_2;

  energy = 0;
  if(type > 2)   energy += P->TerminalAU;
  if(type_2 > 2) energy += P->TerminalAU;

  if(!dangles) return energy;

  ci = ON_SAME_STRAND(i, i + 1, cutpoint);
  cj = ON_SAME_STRAND(j - 1, j, cutpoint);
  cp = ON_SAME_STRAND(p - 1, p, cutpoint);
  cq = ON_SAME_STRAND(q, q + 1, cutpoint);

  d3    = ci  ? P->dangle3[type][si1]   : 0;
  d5    = cj  ? P->dangle5[type][sj1]   : 0;
  d5_2  = cp  ? P->dangle5[type_2][sp1] : 0;
  d3_2  = cq  ? P->dangle3[type_2][sq1] : 0;

  tmm   = (cj && ci) ? P->mismatchExt[type][sj1][si1]   : d5 + d3;
  tmm_2 = (cp && cq) ? P->mismatchExt[type_2][sp1][sq1] : d5_2 + d3_2;

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

int
vrna_E_int_loop(vrna_fold_compound_t *vc,
                int i,
                int j);

int
vrna_eval_int_loop( vrna_fold_compound_t *vc,
                    int i,
                    int j,
                    int k,
                    int l);

FLT_OR_DBL
vrna_exp_E_int_loop(vrna_fold_compound_t *vc,
                int i,
                int j);

FLT_OR_DBL
vrna_exp_E_interior_loop( vrna_fold_compound_t *vc,
                          int i,
                          int j,
                          int k,
                          int l);

int
vrna_E_ext_int_loop(vrna_fold_compound_t *vc,
                    int i,
                    int j,
                    int *ip,
                    int *iq);

int
vrna_E_stack( vrna_fold_compound_t *vc,
              int i,
              int j);


/**
 *  @brief Backtrack a stacked pair closed by @f$ (i,j) @f$
 *
 */
int
vrna_BT_stack(vrna_fold_compound_t *vc,
              int *i,
              int *j,
              int *en,
              vrna_bp_stack_t *bp_stack,
              int *stack_count);
/**
 *  @brief Backtrack an interior loop closed by @f$ (i,j) @f$
 *
 */
int
vrna_BT_int_loop( vrna_fold_compound_t *vc,
                  int *i,
                  int *j,
                  int en,
                  vrna_bp_stack_t *bp_stack,
                  int *stack_count);


/**
 * @}
 */


#endif
