#ifndef __VIENNA_RNA_PACKAGE_LOOP_ENERGIES_H__
#define __VIENNA_RNA_PACKAGE_LOOP_ENERGIES_H__

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <string.h>
#include "ViennaRNA/params.h"
#include "ViennaRNA/fold_vars.h"
#include "ViennaRNA/energy_par.h"
#include "ViennaRNA/constraints.h"
#include "ViennaRNA/gquad.h"

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
INLINE  PRIVATE int E_MLstem( int type,
                              int si1,
                              int sj1,
                              paramT *P);

/**
 *  \def exp_E_MLstem(A,B,C,D)
 *  This is the partition function variant of \ref E_MLstem()
 *  \see E_MLstem()
 *  \return The Boltzmann weighted energy contribution of the introduced multiloop stem
 */
INLINE  PRIVATE double exp_E_MLstem(int type,
                                    int si1,
                                    int sj1,
                                    pf_paramT *P);

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
INLINE  PRIVATE int E_ExtLoop(int type,
                              int si1,
                              int sj1,
                              paramT *P);

/**
 *  \def exp_E_ExtLoop(A,B,C,D)
 *  This is the partition function variant of \ref E_ExtLoop()
 *  \see E_ExtLoop()
 *  \return The Boltzmann weighted energy contribution of the introduced exterior-loop stem
 */
INLINE  PRIVATE double exp_E_ExtLoop( int type,
                                      int si1,
                                      int sj1,
                                      pf_paramT *P);

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
INLINE  PRIVATE int E_IntLoop(int n1,
                              int n2,
                              int type,
                              int type_2,
                              int si1,
                              int sj1,
                              int sp1,
                              int sq1,
                              paramT *P);


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
INLINE  PRIVATE int E_Hairpin(int size,
                              int type,
                              int si1,
                              int sj1,
                              const char *string,
                              paramT *P);

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
INLINE  PRIVATE int E_Stem( int type,
                            int si1,
                            int sj1,
                            int extLoop,
                            paramT *P);

/**
 *  <H2>Compute the Boltzmann weighted energy contribution of a stem branching off a loop-region</H2>
 *  This is the partition function variant of \ref E_Stem()
 *  \see E_Stem()
 *  \note This function is threadsafe
 * 
 *  \return The Boltzmann weighted energy contribution of the branch off the loop
 */
INLINE  PRIVATE double exp_E_Stem(int type,
                                  int si1,
                                  int sj1,
                                  int extLoop,
                                  pf_paramT *P);

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
INLINE  PRIVATE double exp_E_Hairpin( int u,
                                      int type,
                                      short si1,
                                      short sj1,
                                      const char *string,
                                      pf_paramT *P);

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
INLINE  PRIVATE double  exp_E_IntLoop(int u1,
                                      int u2,
                                      int type,
                                      int type2,
                                      short si1,
                                      short sj1,
                                      short sp1,
                                      short sq1,
                                      pf_paramT *P);


/*
#################################
# BEGIN OF FUNCTION DEFINITIONS #
#################################
*/
INLINE PRIVATE int
E_Hairpin(int size,
          int type,
          int si1,
          int sj1,
          const char *string,
          paramT *P){

  int energy;

  if(size <= 30)
    energy = P->hairpin[size];
  else
    energy = P->hairpin[30] + (int)(P->lxc*log((size)/30.));

  if(size < 3) return energy; /* should only be the case when folding alignments */

  if(P->model_details.special_hp){
    if(size == 4){ /* check for tetraloop bonus */
      char tl[7]={0}, *ts;
      strncpy(tl, string, 6);
      if ((ts=strstr(P->Tetraloops, tl)))
        return (P->Tetraloop_E[(ts - P->Tetraloops)/7]);
    }
    else if(size == 6){
      char tl[9]={0}, *ts;
      strncpy(tl, string, 8);
      if ((ts=strstr(P->Hexaloops, tl)))
        return (energy = P->Hexaloop_E[(ts - P->Hexaloops)/9]);
    }
    else if(size == 3){
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

INLINE PRIVATE int
E_hp_loop(int i,
          int j,
          vrna_fold_compound *vc){

  int u = j - i - 1;
  int e;

  short   *S        = vc->sequence_encoding;
  int     *idx      = vc->jindx;
  int     type      = vc->ptype[idx[j]+i];
  paramT  *P        = vc->params;
  char    hc        = vc->hc->matrix[idx[j]+i];
  int     *hc_up    = vc->hc->up_hp;
  soft_constraintT  *sc = vc->sc;

  /* is this base pair allowed to close a hairpin loop ? */
  if(hc & IN_HP_LOOP){
    /* are all nucleotides in the loop allowed to be unpaired ? */
    if(hc_up[i+1] >= u){
      e = E_Hairpin(u, type, S[i+1], S[j-1], vc->sequence+i-1, P);
      if(sc){
        if(sc->free_energies)
          e += sc->free_energies[i+1][u];
        if(sc->f)
          e += sc->f(i, j, i, j, VRNA_DECOMP_PAIR_HP, sc->data);
      }
      return e;
    }
  }
  return INF;
}

INLINE PRIVATE int
E_int_loop( int i,
            int j,
            vrna_fold_compound *vc){

  int q, p, j_q, p_i, pq, *c_pq, max_q, max_p, tmp, type, type_2, *rtype, noGUclosure, no_close, energy;
  short *S_p1, *S_q1;
  char              *ptype_pq;
  char              *hc_pq;
  char              *ptype        = vc->ptype;
  short             *S            = vc->sequence_encoding;
  short             S_i1          = S[i+1];
  short             S_j1          = S[j-1];
  int               *indx         = vc->jindx;
  char              *hc           = vc->hc->matrix;
  int               *hc_up        = vc->hc->up_int;
  soft_constraintT  *sc           = vc->sc; 
  paramT            *P            = vc->params;
  int               ij            = indx[j] + i;
  int               hc_decompose  = hc[ij];
  int               e             = INF;
  int               *c            = vc->matrices->c;
  int               *ggg          = vc->matrices->ggg;
  int               with_gquad    = P->model_details.gquad;

  /* CONSTRAINED INTERIOR LOOP start */
  if(hc_decompose & IN_INT_LOOP){

    type        = ptype[ij];
    rtype       = &(P->model_details.rtype[0]);
    noGUclosure = P->model_details.noGUclosure;
    no_close    = (((type==3)||(type==4))&&noGUclosure);
    max_q       = i+TURN+2;
    max_q       = MAX2(max_q, j - MAXLOOP - 1);
    for(q = j - 1; q >= max_q; q--){
      j_q = j - q - 1;

      if(hc_up[q+1] < j_q) break;

      pq        = indx[q] + i + 1;
      p_i       = 0;
      max_p     = i + 1;
      tmp       = i + 1 + MAXLOOP - j_q;
      max_p     = MAX2(max_p, tmp);
      tmp       = q - TURN;
      max_p     = MIN2(max_p, tmp);
      tmp       = i + 1 + hc_up[i + 1];
      max_p     = MIN2(max_p, tmp);
      hc_pq     = hc + pq;
      c_pq      = c + pq;
      ptype_pq  = ptype + pq;
      S_p1      = S + i;
      S_q1      = S + q + 1;
      for(p = i+1; p <= max_p; p++){

        /* discard this configuration if (p,q) is not allowed to be enclosed pair of an interior loop */
        if(*hc_pq & IN_INT_LOOP_ENC){

          type_2 = rtype[*ptype_pq];

          if (noGUclosure)
            if (no_close||(type_2==3)||(type_2==4))
              if ((p>i+1)||(q<j-1)) continue;  /* continue unless stack */

          energy = E_IntLoop(p_i, j_q, type, type_2, S_i1, S_j1, *S_p1, *S_q1, P);
          energy += *c_pq;

          /* add soft constraints */
          if(sc){
            if(sc->free_energies)
              energy += sc->free_energies[i+1][p_i] + sc->free_energies[q+1][j_q];
            if(sc->f)
              energy += sc->f(i, j, p, q, VRNA_DECOMP_PAIR_IL, sc->data);
          }

          e = MIN2(e, energy);

        }
        hc_pq++;    /* get hc[pq + 1] */
        c_pq++;     /* get c[pq + 1] */
        p_i++;      /* increase unpaired region [i+1...p-1] */
        ptype_pq++; /* get ptype[pq + 1] */
        S_p1++;
      } /* end q-loop */
    } /* end p-loop */

    if(with_gquad){
      /* include all cases where a g-quadruplex may be enclosed by base pair (i,j) */
      if (!no_close) {
        energy = E_GQuad_IntLoop(i, j, type, S, ggg, indx, P);
        e = MIN2(e, energy);
      }
    }

  }
  return e;
}

INLINE PRIVATE int
E_mb_loop_fast( int i,
                int j,
                vrna_fold_compound *vc,
                int *dmli1,
                int *dmli2){

  int k, decomp, MLenergy, en, type, type_2, tt;

  char              *ptype  = vc->ptype;
  short             *S      = vc->sequence_encoding;
  short             S_i1    = S[i+1];
  short             S_j1    = S[j-1];
  int               *indx   = vc->jindx;
  char              *hc     = vc->hc->matrix;
  int               *hc_up  = vc->hc->up_ml;
  soft_constraintT  *sc     = vc->sc;
  int               *c      = vc->matrices->c;
  int               *fML    = vc->matrices->fML;
  paramT            *P      = vc->params;

  int ij            = indx[j] + i;
  int hc_decompose  = hc[ij];
  int e             = INF;
  int dangle_model  = P->model_details.dangles;
  int *rtype        = &(P->model_details.rtype[0]);

  type              = ptype[ij];


  if(hc_decompose & IN_MB_LOOP){
    decomp = dmli1[j-1];
    tt = rtype[type];
    switch(dangle_model){
      /* no dangles */
      case 0:   decomp += E_MLstem(tt, -1, -1, P);
                break;

      /* double dangles */
      case 2:   decomp += E_MLstem(tt, S_j1, S_i1, P);
                break;

      /* normal dangles, aka dangles = 1 || 3 */
      default:  decomp += E_MLstem(tt, -1, -1, P);
                if(hc_up[i+1]){
                  en = dmli2[j-1] + E_MLstem(tt, -1, S_i1, P) + P->MLbase;
                  if(sc){
                    if(sc->free_energies)
                      en += sc->free_energies[i+1][1];
                  }
                  decomp = MIN2(decomp, en);
                }
                if(hc_up[j-1] && hc_up[i+1]){
                  en = dmli2[j-2] + E_MLstem(tt, S_j1, S_i1, P) + 2*P->MLbase;
                  decomp = MIN2(decomp, en);
                }
                if(hc_up[j-1]){
                  en = dmli1[j-2] + E_MLstem(tt, S_j1, -1, P) + P->MLbase;
                  decomp = MIN2(decomp, en);
                }
                break;
    }
    MLenergy = decomp + P->MLclosing;
    e = MIN2(e, MLenergy);

    /* coaxial stacking of (i.j) with (i+1.k) or (k+1.j-1) */
    if (dangle_model==3) {
      int i1k, k1j1;
      decomp = INF;
      k1j1  = indx[j-1] + i + 2 + TURN + 1;
      for (k = i+2+TURN; k < j-2-TURN; k++, k1j1++){
        i1k   = indx[k] + i + 1;
        if(hc[i1k] & IN_MB_LOOP_ENC){
          type_2  = rtype[ptype[i1k]];
          en      = c[i1k]+P->stack[type][type_2]+fML[k1j1];
          decomp  = MIN2(decomp, en);
        }
        if(hc[k1j1] & IN_MB_LOOP_ENC){
          type_2  = rtype[ptype[k1j1]];
          en      = c[k1j1]+P->stack[type][type_2]+fML[i1k];
          decomp  = MIN2(decomp, en);
        }
      }
      /* no TermAU penalty if coax stack */
      decomp += 2*P->MLintern[1] + P->MLclosing;
      e = MIN2(e, decomp);
    }
  }
  return e;
}

INLINE PRIVATE int
E_ml_rightmost_stem(int i,
                    int j,
                    int length,
                    int type,
                    short *S,
                    int *indx,
                    char *hc,
                    int *hc_up,
                    soft_constraintT *sc,
                    int *c,
                    int *fm,
                    paramT *P){

  int en;
  int ij            = indx[j] + i;
  int hc_decompose  = hc[ij];
  int dangle_model  = P->model_details.dangles;
  int e             = INF;

  if(hc_decompose & IN_MB_LOOP_ENC){
    e = c[ij];
    switch(dangle_model){
      case 2:   e += E_MLstem(type, (i==1) ? S[length] : S[i-1], S[j+1], P);
                break;
      default:  e += E_MLstem(type, -1, -1, P);
                break;
    }
  }

  if(hc_up[j]){
    en = fm[indx[j - 1] + i] + P->MLbase;
    if(sc)
      if(sc->free_energies)
        en += sc->free_energies[j][1];

    e = MIN2(e, en);
  }

  return e;
}

INLINE PRIVATE int
E_ml_stems_fast(int i,
                int j,
                vrna_fold_compound *vc,
                int *fmi,
                int *dmli){

  int k, en, decomp, mm5, mm3, type_2, k1j;

  int               length        = (int)vc->length;
  char              *ptype        = vc->ptype;
  short             *S            = vc->sequence_encoding;
  int               *indx         = vc->jindx;
  char              *hc           = vc->hc->matrix;
  int               *hc_up        = vc->hc->up_ml;
  soft_constraintT  *sc           = vc->sc;
  int               *c            = vc->matrices->c;
  int               *fm           = vc->matrices->fML;
  paramT            *P            = vc->params;
  int               ij            = indx[j] + i;
  int               dangle_model  = P->model_details.dangles;
  int               type          = ptype[ij];
  int               *rtype        = &(P->model_details.rtype[0]);
  int               circular      = P->model_details.circ;
  int               e             = INF;

  /*  extension with one unpaired nucleotide at the left
      or full branch of (i,j)
  */
  e = E_ml_rightmost_stem(i,j,length,type,S,indx, hc,hc_up,sc,c,fm,P);

  /*  extension with one unpaired nucleotide at 5' site
      and all other variants which are needed for odd
      dangle models
  */
  switch(dangle_model){
    /* no dangles */
    case 0:   /* fall through */

    /* double dangles */
    case 2:   if(hc_up[i]){
                en = fm[ij + 1] + P->MLbase;
                if(sc)
                  if(sc->free_energies)
                    en += sc->free_energies[i][1];
                e = MIN2(e, en);
              }
              break;

    /* normal dangles, aka dangle_model = 1 || 3 */
    default:  mm5 = ((i>1) || circular) ? S[i] : -1;
              mm3 = ((j<length) || circular) ? S[j] : -1;
              if(hc_up[i]){
                en = fm[ij+1] + P->MLbase;
                if(sc)
                  if(sc->free_energies)
                    en += sc->free_energies[i][1];
                e = MIN2(e, en);
                if(hc[ij+1] & IN_MB_LOOP_ENC){
                  type = ptype[ij+1];
                  en = c[ij+1] + E_MLstem(type, mm5, -1, P) + P->MLbase;
                  if(sc)
                    if(sc->free_energies)
                      en += sc->free_energies[i][1];
                  e = MIN2(e, en);
                }
              }
              if(hc_up[j]){
                if(hc[indx[j-1]+i] & IN_MB_LOOP_ENC){
                  type = ptype[indx[j-1]+i];
                  en = c[indx[j-1]+i] + E_MLstem(type, -1, mm3, P) + P->MLbase;
                  if(sc)
                    if(sc->free_energies)
                      en += sc->free_energies[j][1];
                  e = MIN2(e, en);
                }
              }
              if(hc[indx[j-1]+i+1] & IN_MB_LOOP_ENC){
                if(hc_up[i] && hc_up[j]){
                  type = ptype[indx[j-1]+i+1];
                  en = c[indx[j-1]+i+1] + E_MLstem(type, mm5, mm3, P) + 2*P->MLbase;
                  if(sc)
                    if(sc->free_energies)
                      en += sc->free_energies[j][1] + sc->free_energies[i][1];
                  e = MIN2(e, en);
                }
              }
              break;
  }

  /* modular decomposition -------------------------------*/
  k1j = indx[j] + i + TURN + 2;
  for (decomp = INF, k = i + 1 + TURN; k <= j - 2 - TURN; k++, k1j++){
    en = fmi[k] + fm[k1j];
    decomp = MIN2(decomp, en);
  }
  dmli[j] = decomp;               /* store for use in fast ML decompositon */
  e = MIN2(e, decomp);

  /* coaxial stacking */
  if (dangle_model==3) {
    /* additional ML decomposition as two coaxially stacked helices */
    int ik, k1j;
    for (k1j = indx[j]+i+TURN+2, decomp = INF, k = i+1+TURN; k <= j-2-TURN; k++, k1j++){
      ik = indx[k]+i;
      if((hc[ik] & IN_MB_LOOP_ENC) && (hc[k1j] & IN_MB_LOOP_ENC)){
        type    = rtype[ptype[ik]];
        type_2  = rtype[ptype[k1j]];
        en      = c[ik] + c[k1j] + P->stack[type][type_2];
        decomp  = MIN2(decomp, en);
      }
    }

    decomp += 2*P->MLintern[1];        /* no TermAU penalty if coax stack */
#if 0
        /* This is needed for Y shaped ML loops with coax stacking of
           interior pairts, but backtracking will fail if activated */
        DMLi[j] = MIN2(DMLi[j], decomp);
        DMLi[j] = MIN2(DMLi[j], DMLi[j-1]+P->MLbase);
        DMLi[j] = MIN2(DMLi[j], DMLi1[j]+P->MLbase);
        new_fML = MIN2(new_fML, DMLi[j]);
#endif
    e = MIN2(e, decomp);
  }

  fmi[j] = e;

  return e;
}

INLINE PRIVATE void
E_ext_loop_5( vrna_fold_compound *vc){

  int en, i, j, ij, type;
  int               length        = (int)vc->length;
  char              *ptype        = vc->ptype;
  short             *S            = vc->sequence_encoding;
  int               *indx         = vc->jindx;
  char              *hc           = vc->hc->matrix;
  int               *hc_up        = vc->hc->up_ext;
  soft_constraintT  *sc           = vc->sc;
  int               *f5           = vc->matrices->f5;
  int               *c            = vc->matrices->c;
  paramT            *P            = vc->params;
  int               dangle_model  = P->model_details.dangles;
  int               *ggg          = vc->matrices->ggg;
  int               with_gquad    = P->model_details.gquad;

  f5[0] = 0;
  for(i = 1; i <= TURN + 1; i++){
    if(hc_up[i]){
      f5[i] = f5[i-1];
      if(sc)
        if(sc->free_energies)
          f5[i] += sc->free_energies[i][1];
    } else {
      f5[i] = INF;
    }
  }

  /* duplicated code may be faster than conditions inside loop ;) */
  switch(dangle_model){
    /* dont use dangling end and mismatch contributions at all */
    case 0:   for(j=TURN+2; j<=length; j++){
                if(hc_up[j]){
                  f5[j] = f5[j-1];
                  if(sc)
                    if(sc->free_energies)
                      f5[j] += sc->free_energies[j][1];
                }
                for (i=j-TURN-1; i>1; i--){
                  ij = indx[j]+i;
                  if(!(hc[ij] & IN_EXT_LOOP)) continue;

                  if(with_gquad){
                    f5[j] = MIN2(f5[j], f5[i-1] + ggg[indx[j]+i]);
                  }

                  en    = f5[i-1] + c[ij] + E_ExtLoop(ptype[ij], -1, -1, P);
                  f5[j] = MIN2(f5[j], en);
                }
                ij = indx[j] + 1;
                if(!(hc[ij] & IN_EXT_LOOP)) continue;

                if(with_gquad){
                  f5[j] = MIN2(f5[j], ggg[indx[j]+1]);
                }

                en    = c[ij] + E_ExtLoop(ptype[ij], -1, -1, P);
                f5[j] = MIN2(f5[j], en);
              }
              break;

    /* always use dangles on both sides */
    case 2:   for(j=TURN+2; j<length; j++){
                if(hc_up[j]){
                  f5[j] = f5[j-1];
                  if(sc)
                    if(sc->free_energies)
                      f5[j] += sc->free_energies[j][1];
                }
                for (i=j-TURN-1; i>1; i--){
                  ij = indx[j] + i;
                  if(!(hc[ij] & IN_EXT_LOOP)) continue;

                  if(with_gquad){
                    f5[j] = MIN2(f5[j], f5[i-1] + ggg[indx[j]+i]);
                  }

                  en    = f5[i-1] + c[ij] + E_ExtLoop(ptype[ij], S[i-1], S[j+1], P);
                  f5[j] = MIN2(f5[j], en);
                }
                ij = indx[j] + 1;
                if(!(hc[ij] & IN_EXT_LOOP)) continue;

                if(with_gquad){
                  f5[j] = MIN2(f5[j], ggg[indx[j]+1]);
                }

                en    = c[ij] + E_ExtLoop(ptype[ij], -1, S[j+1], P);
                f5[j] = MIN2(f5[j], en);
              }
              if(hc_up[length]){
                f5[length] = f5[length-1];
                if(sc)
                  if(sc->free_energies)
                    f5[length] += sc->free_energies[length][1];
              }
              for (i=length-TURN-1; i>1; i--){
                ij = indx[length] + i;
                if(!(hc[ij] & IN_EXT_LOOP)) continue;

                if(with_gquad){
                  f5[length] = MIN2(f5[length], f5[i-1] + ggg[indx[length]+i]);
                }

                en          = f5[i-1] + c[ij] + E_ExtLoop(ptype[ij], S[i-1], -1, P);
                f5[length]  = MIN2(f5[length], en);
              }
              ij = indx[length] + 1;
              if(!(hc[ij] & IN_EXT_LOOP)) break;

              if(with_gquad){
                f5[length] = MIN2(f5[length], ggg[indx[length]+1]);
              }

              en          = c[ij] + E_ExtLoop(ptype[ij], -1, -1, P);
              f5[length]  = MIN2(f5[length], en);
              break;

    /* normal dangles, aka dangle_model = 1 || 3 */
    default:  for(j=TURN+2; j<=length; j++){
                if(hc_up[j])
                  f5[j] = f5[j-1];
                for (i=j-TURN-1; i>1; i--){
                  ij = indx[j] + i;
                  if(hc[ij] & IN_EXT_LOOP){

                    if(with_gquad){
                      f5[j] = MIN2(f5[j], f5[i-1] + ggg[indx[j]+i]);
                    }

                    type  = ptype[ij];
                    en    = f5[i-1] + c[ij] + E_ExtLoop(type, -1, -1, P);
                    f5[j] = MIN2(f5[j], en);
                    if(hc_up[i-1]){
                      en    = f5[i-2] + c[ij] + E_ExtLoop(type, S[i-1], -1, P);
                      f5[j] = MIN2(f5[j], en);
                    }
                  }
                  ij = indx[j-1] + i;
                  if(hc[ij] & IN_EXT_LOOP){
                    if(hc_up[j]){
                      type  = ptype[ij];
                      en    = f5[i-1] + c[ij] + E_ExtLoop(type, -1, S[j], P);
                      f5[j] = MIN2(f5[j], en);
                      if(hc_up[i-1]){
                        en    = f5[i-2] + c[ij] + E_ExtLoop(type, S[i-1], S[j], P);
                        f5[j] = MIN2(f5[j], en);
                      }
                    }
                  }
                }
                ij = indx[j] + 1;
                if(hc[ij] & IN_EXT_LOOP){

                  if(with_gquad){
                    f5[j] = MIN2(f5[j], ggg[indx[j]+1]);
                  }

                  type  = ptype[ij];
                  en    = c[ij] + E_ExtLoop(type, -1, -1, P);
                  f5[j] = MIN2(f5[j], en);
                }
                ij = indx[j-1] + 1;
                if(hc[ij] & IN_EXT_LOOP){
                  if(hc_up[j]){
                    type  = ptype[ij];
                    en    = c[ij] + E_ExtLoop(type, -1, S[j], P);
                    f5[j] = MIN2(f5[j], en);
                  }
                }
              }
  }
}


INLINE PRIVATE double
exp_E_Hairpin(int u,
              int type,
              short si1,
              short sj1,
              const char *string,
              pf_paramT *P){

  double q, kT;
  kT = P->kT;   /* kT in cal/mol  */

  if(u <= 30)
    q = P->exphairpin[u];
  else
    q = P->exphairpin[30] * exp( -(P->lxc*log( u/30.))*10./kT);

  if(u < 3) return q; /* should only be the case when folding alignments */

  if(P->model_details.special_hp){
    if(u==4){
      char tl[7]={0,0,0,0,0,0,0}, *ts;
      strncpy(tl, string, 6);
      if ((ts=strstr(P->Tetraloops, tl))){
        if(type != 7)
          return (P->exptetra[(ts-P->Tetraloops)/7]);
        else
          q *= P->exptetra[(ts-P->Tetraloops)/7];
      }
    }
    else if(u==6){
      char tl[9]={0,0,0,0,0,0,0,0,0}, *ts;
      strncpy(tl, string, 6);
      if ((ts=strstr(P->Hexaloops, tl)))
        return  (P->exphex[(ts-P->Hexaloops)/9]);
    }
    else if(u==3){
      char tl[6]={0,0,0,0,0,0}, *ts;
      strncpy(tl, string, 5);
      if ((ts=strstr(P->Triloops, tl)))
        return (P->exptri[(ts-P->Triloops)/6]);
      if (type>2)
        return q * P->expTermAU;
      else
        return q;
    }
  }
  q *= P->expmismatchH[type][si1][sj1];

  return q;
}

INLINE PRIVATE int
E_IntLoop(int n1,
          int n2,
          int type,
          int type_2,
          int si1,
          int sj1,
          int sp1,
          int sq1,
          paramT *P){

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

INLINE  PRIVATE double exp_E_IntLoop(int u1, int u2, int type, int type2, short si1, short sj1, short sp1, short sq1, pf_paramT *P){
  int ul, us, no_close = 0;
  double z = 0.;

  if ((no_closingGU) && ((type2==3)||(type2==4)||(type==3)||(type==4)))
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

  if(si1 >= 0 && sj1 >= 0)
    energy = (extLoop) ? P->expmismatchExt[type][si1][sj1] : P->expmismatchM[type][si1][sj1];
  else
    energy = d5 * d3;

  if(type > 2)
    energy *= P->expTermAU;

  if(!extLoop) energy *= P->expMLintern[type];
  return energy;
}

INLINE PRIVATE double
exp_E_MLstem( int type,
              int si1,
              int sj1,
              pf_paramT *P){

  double energy = 1.0;
  if(si1 >= 0 && sj1 >= 0){
    energy = P->expmismatchM[type][si1][sj1];
  }
  else if(si1 >= 0){
    energy = P->expdangle5[type][si1];
  }
  else if(sj1 >= 0){
    energy = P->expdangle3[type][sj1];
  }

  if(type > 2)
    energy *= P->expTermAU;

  energy *= P->expMLintern[type];
  return energy;
}

INLINE PRIVATE double
exp_E_ExtLoop(int type,
              int si1,
              int sj1,
              pf_paramT *P){

  double energy = 1.0;
  if(si1 >= 0 && sj1 >= 0){
    energy = P->expmismatchExt[type][si1][sj1];
  }
  else if(si1 >= 0){
    energy = P->expdangle5[type][si1];
  }
  else if(sj1 >= 0){
    energy = P->expdangle3[type][sj1];
  }

  if(type > 2)
    energy *= P->expTermAU;

  return energy;
}

INLINE PRIVATE int
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
              paramT *P){

  int energy, ci, cj, cp, cq, d3, d5, d5_2, d3_2, tmm, tmm_2;

  energy = 0;
  if(type > 2)   energy += P->TerminalAU;
  if(type_2 > 2) energy += P->TerminalAU;

  if(!dangles) return energy;

  ci = (i>=cutpoint)||((i+1)<cutpoint);
  cj = ((j-1)>=cutpoint)||(j<cutpoint);
  cp = ((p-1)>=cutpoint)||(p<cutpoint);
  cq = (q>=cutpoint)||((q+1)<cutpoint);

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

#endif
