#ifndef VIENNA_RNA_PACKAGE_HAIRPIN_LOOPS_H
#define VIENNA_RNA_PACKAGE_HAIRPIN_LOOPS_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <string.h>
#include <ViennaRNA/fold_vars.h>
#include <ViennaRNA/energy_par.h>
#include <ViennaRNA/data_structures.h>
#include <ViennaRNA/params.h>
#include <ViennaRNA/constraints.h>
#include <ViennaRNA/exterior_loops.h>
#include <ViennaRNA/gquad.h>

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
 *  @addtogroup   loops
 *
 *  @{
 *
 *  @file hairpin_loops.h
 *  @brief Energy evaluation of hairpin loops for MFE and partition function calculations
 */


/**
 *  @brief Compute the Energy of a hairpin-loop
 *
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
 *  @note The parameter sequence should contain the sequence of the loop in capital letters of the nucleic acid
 *  alphabet if the loop size is below 7. This is useful for unusually stable tri-, tetra- and hexa-loops
 *  which are treated differently (based on experimental data) if they are tabulated.
 *  @see scale_parameters()
 *  @see vrna_param_t
 *  @warning Not (really) thread safe! A threadsafe implementation will replace this function in a future release!\n
 *  Energy evaluation may change due to updates in global variable "tetra_loop"
 * 
 *  @param  size  The size of the loop (number of unpaired nucleotides)
 *  @param  type  The pair type of the base pair closing the hairpin
 *  @param  si1   The 5'-mismatching nucleotide
 *  @param  sj1   The 3'-mismatching nucleotide
 *  @param  string  The sequence of the loop
 *  @param  P     The datastructure containing scaled energy parameters
 *  @return The Free energy of the Hairpin-loop in dcal/mol
 */
INLINE  PRIVATE int E_Hairpin(int size,
                              int type,
                              int si1,
                              int sj1,
                              const char *string,
                              vrna_param_t *P);

/**
 *  @brief Compute Boltzmann weight @f$e^{-\Delta G/kT} @f$ of a hairpin loop
 *
 *  multiply by scale[u+2]
 *  @see get_scaled_pf_parameters()
 *  @see vrna_exp_param_t
 *  @see E_Hairpin()
 *  @warning Not (really) thread safe! A threadsafe implementation will replace this function in a future release!\n
 *  Energy evaluation may change due to updates in global variable "tetra_loop"
 * 
 *  @param  u       The size of the loop (number of unpaired nucleotides)
 *  @param  type    The pair type of the base pair closing the hairpin
 *  @param  si1     The 5'-mismatching nucleotide
 *  @param  sj1     The 3'-mismatching nucleotide
 *  @param  string  The sequence of the loop
 *  @param  P       The datastructure containing scaled Boltzmann weights of the energy parameters
 *  @return The Boltzmann weight of the Hairpin-loop
 */
INLINE  PRIVATE double
exp_E_Hairpin(  int u,
                int type,
                short si1,
                short sj1,
                const char *string,
                vrna_exp_param_t *P);


PRIVATE INLINE int
vrna_eval_hp_loop(vrna_fold_compound_t *vc,
                  int i,
                  int j);

PRIVATE INLINE int
vrna_eval_ext_hp_loop(vrna_fold_compound_t *vc,
                      int i,
                      int j);

/*
#################################
# BEGIN OF FUNCTION DEFINITIONS #
#################################
*/

PRIVATE INLINE int
E_Hairpin(int size,
          int type,
          int si1,
          int sj1,
          const char *string,
          vrna_param_t *P){

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

/**
 *  @brief  Evaluate the free energy of a hairpin loop
 *          and consider possible hard constraints
 *
 *  @note This function is polymorphic! The provided #vrna_fold_compound_t may be of type
 *  #VRNA_VC_TYPE_SINGLE or #VRNA_VC_TYPE_ALIGNMENT
 *
 */
PRIVATE INLINE int
vrna_E_hp_loop( vrna_fold_compound_t *vc,
                int i,
                int j){

  int   u, *hc_up;
  char  eval_loop;
  char  (*f)(vrna_fold_compound_t *, int, int, int, int, char);

  u         = j - i - 1;
  hc_up     = vc->hc->up_hp;
  f         = vc->hc->f;
  eval_loop = vc->hc->matrix[vc->jindx[j] + i] & VRNA_CONSTRAINT_CONTEXT_HP_LOOP;

#ifdef WITH_GEN_HC
  if(f)
    eval_loop = (f(vc, i, j, i+1, j-1, VRNA_DECOMP_PAIR_HP)) ? eval_loop : (char)0;
#endif

  /* is this base pair allowed to close a hairpin (like) loop ? */
  if(eval_loop){
    /* are all nucleotides in the loop allowed to be unpaired ? */
    if(hc_up[i+1] >= u){
      return vrna_eval_hp_loop(vc, i, j);
    }
  }
  return INF;
}

/**
 *  @brief  Evaluate the free energy of an exterior hairpin loop
 *          and consider possible hard constraints
 */
PRIVATE INLINE int
vrna_E_ext_hp_loop( vrna_fold_compound_t *vc,
                    int i,
                    int j){

  int   u, *hc_up;
  char  eval_loop;
  char  (*f)(vrna_fold_compound_t *, int, int, int, int, char);
  
  u         = vc->length - j + i - 1;
  hc_up     = vc->hc->up_hp;
  f         = vc->hc->f;
  eval_loop = vc->hc->matrix[vc->jindx[j] + i] & VRNA_CONSTRAINT_CONTEXT_HP_LOOP;

#ifdef WITH_GEN_HC
  if(f)
    eval_loop = (f(vc, j, i, j+1, i-1, VRNA_DECOMP_PAIR_HP)) ? eval_loop : (char)0;
#endif

  /* is this base pair allowed to close a hairpin (like) loop ? */
  if(eval_loop){
    /* are all nucleotides in the loop allowed to be unpaired ? */
    if(hc_up[j+1] >= u){
      return vrna_eval_ext_hp_loop(vc, i, j);
    }
  }
  return INF;
}

/**
 *  @brief Evaluate free energy of an exterior hairpin loop
 *
 *  @ingroup eval
 *
 */
PRIVATE INLINE int
vrna_eval_ext_hp_loop(vrna_fold_compound_t *vc,
                      int i,
                      int j){

  int               u, e, type;
  char              loopseq[10];
  short             *S;
  vrna_param_t            *P;
  vrna_sc_t         *sc;
  vrna_md_t         *md;

  S     = vc->sequence_encoding;
  P     = vc->params;
  sc    = vc->sc;
  md    = &(P->model_details);
  u     = vc->length - j + i - 1;
  type  = md->pair[S[j]][S[i]];

  if(type == 0)
    type = 7;

  if (u<7) {
    strcpy(loopseq , vc->sequence + j - 1);
    strncat(loopseq, vc->sequence, i);
  }

  e = E_Hairpin(u, type, S[j + 1], S[i - 1],  loopseq, P);

  if(sc){
    if(sc->free_energies)
      e +=  sc->free_energies[j + 1][vc->length - j]
            + sc->free_energies[1][i - 1];

    if(sc->f)
      e += sc->f(j, i, j-1, i+1, VRNA_DECOMP_PAIR_HP, sc->data);
  }

  return e;
}

/**
 *  @brief Evaluate free energy of a hairpin loop
 *
 *  @ingroup eval
 *
 *  @note This function is polymorphic! The provided #vrna_fold_compound_t may be of type
 *  #VRNA_VC_TYPE_SINGLE or #VRNA_VC_TYPE_ALIGNMENT
 *
 *  @param  vc  The #vrna_fold_compound_t for the particular energy evaluation
 *  @param  i   5'-position of the base pair
 *  @param  j   3'-position of the base pair
 *  @returns    Free energy of the hairpin loop closed by @f$ (i,j) @f$ in deka-kal/mol
 */
PRIVATE INLINE int
vrna_eval_hp_loop(vrna_fold_compound_t *vc,
                  int i,
                  int j){

  int             u, e, s, ij, cp, type, *types, *idx, n_seq;
  short           *S, **SS, **S5, **S3;
  char            **Ss;
  unsigned short  **a2s;
  vrna_param_t    *P;
  vrna_sc_t       *sc, **scs;
  vrna_md_t       *md;

  cp  = vc->cutpoint;
  idx = vc->jindx;
  P   = vc->params;
  md  = &(P->model_details);
  e   = INF;

  switch(vc->type){
    /* single sequences and cofolding hybrids */
    case  VRNA_VC_TYPE_SINGLE:    S     = vc->sequence_encoding;
                                  sc    = vc->sc;
                                  u     = j - i - 1;
                                  ij    = idx[j] + i;
                                  type  = md->pair[S[i]][S[j]];

                                  if(type == 0)
                                    type = 7;

                                  if((cp < 0) || ON_SAME_STRAND(i, j, cp)){ /* regular hairpin loop */
                                    e = E_Hairpin(u, type, S[i+1], S[j-1], vc->sequence+i-1, P);
                                  } else { /* hairpin-like exterior loop (for cofolding) */
                                    short si, sj;
                                    si  = ON_SAME_STRAND(i, i + 1, cp) ? S[i+1] : -1;
                                    sj  = ON_SAME_STRAND(j - 1, j, cp) ? S[j-1] : -1;
                                    if (md->dangles)
                                      e = E_ExtLoop(md->rtype[type], sj, si, P);
                                    else
                                      e = E_ExtLoop(md->rtype[type], -1, -1, P);
                                  }

                                  /* add soft constraints */
                                  if(sc){
                                    if(sc->free_energies)
                                      e += sc->free_energies[i+1][u];

                                    if(sc->en_basepair){
                                      e += sc->en_basepair[ij];
                                    }
                                    if(sc->f)
                                      e += sc->f(i, j, i+1, j-1, VRNA_DECOMP_PAIR_HP, sc->data);
                                  }
                                  break;
    /* sequence alignments */
    case  VRNA_VC_TYPE_ALIGNMENT: SS    = vc->S;                                                               
                                  S5    = vc->S5;     /*S5[s][i] holds next base 5' of i in sequence s*/
                                  S3    = vc->S3;     /*Sl[s][i] holds next base 3' of i in sequence s*/
                                  Ss    = vc->Ss;                                                       
                                  a2s   = vc->a2s;                                                      
                                  scs   = vc->scs;
                                  n_seq = vc->n_seq;
                                  types = (int *)vrna_alloc(sizeof(int) * n_seq);

                                  for (s=0; s<n_seq; s++) {
                                    types[s] = md->pair[SS[s][i]][SS[s][j]];
                                    if (types[s]==0) types[s]=7;
                                  }
                                  
                                  for(e = s = 0; s < n_seq; s++){
                                    u = a2s[s][j-1] - a2s[s][i];
                                    e += (u < 3) ? 600 : E_Hairpin(u, types[s], S3[s][i], S5[s][j], Ss[s]+(a2s[s][i-1]), P);  /* ??? really 600 ??? */

                                  }
                                  
                                  if(scs)
                                    for(s = 0; s < n_seq; s++){
                                      if(scs[s]){
                                        u = a2s[s][j-1]-a2s[s][i];

                                        if(scs[s]->free_energies)
                                          e += scs[s]->free_energies[a2s[s][i]+1][u];

                                        if(scs[s]->en_basepair)
                                          e += scs[s]->en_basepair[idx[j] + i];

                                        if(scs[s]->f)
                                          e += scs[s]->f(i, j, i+1, j-1, VRNA_DECOMP_PAIR_HP, scs[s]->data);
                                      }
                                    }

                                  free(types);

                                  break;
    /* nothing */
    default:                      break;
  }


  return e;
}

/*
*************************************
* Partition function variants below *
*************************************
*/


PRIVATE INLINE double
exp_E_Hairpin(int u,
              int type,
              short si1,
              short sj1,
              const char *string,
              vrna_exp_param_t *P){

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
      strncpy(tl, string, 8);
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

/**
 *  @brief High-Level function for hairpin loop energy evaluation (partition function variant)
 *
 *  @see E_hp_loop() for it's free energy counterpart
*/
PRIVATE INLINE double
vrna_exp_E_hp_loop( vrna_fold_compound_t *vc,
                    int i,
                    int j){

  int         u, ij, type;
  char        hc;
  double      q;

  int               cp      = vc->cutpoint;
  short             *S      = vc->sequence_encoding;
  int               *idx    = vc->jindx;
  vrna_exp_param_t  *P      = vc->exp_params;
  int               *hc_up  = vc->hc->up_hp;
  vrna_sc_t         *sc     = vc->sc;
  FLT_OR_DBL        *scale  = vc->exp_matrices->scale;

  q     = 0.;
  u     = j - i - 1;
  ij    = idx[j] + i;
  type  = vc->ptype[ij];
  hc    = vc->hc->matrix[ij];

  /* is this base pair allowed to close a hairpin (like) loop ? */
  if(hc & VRNA_CONSTRAINT_CONTEXT_HP_LOOP){
    /* are all nucleotides in the loop allowed to be unpaired ? */
    if(hc_up[i+1] >= u){

      if((cp < 0) || ON_SAME_STRAND(i, j, cp)){ /* regular hairpin loop */
        q = exp_E_Hairpin(u, type, S[i+1], S[j-1], vc->sequence+i-1, P)
            * scale[u+2];
      } else { /* hairpin-like exterior loop (for cofolding) */
        
      }

      /* add soft constraints */
      if(sc){
        if(sc->boltzmann_factors)
          q *= sc->boltzmann_factors[i+1][u];

        if(sc->exp_en_basepair)
          q *= sc->exp_en_basepair[ij];

        if(sc->exp_f)
          q *= sc->exp_f(i, j, i, j, VRNA_DECOMP_PAIR_HP, sc->data);
      }
    }
  }
  return q;
}

/**
 *  @brief Backtrack a hairpin loop closed by @f$ (i,j) @f$
 *
 *  @note This function is polymorphic! The provided #vrna_fold_compound_t may be of type
 *  #VRNA_VC_TYPE_SINGLE or #VRNA_VC_TYPE_ALIGNMENT
 *
 */
PRIVATE INLINE int
vrna_BT_hp_loop(vrna_fold_compound_t *vc,
                int i,
                int j,
                int en,
                bondT *bp_stack,
                int   *stack_count){

  int       e;
  vrna_sc_t *sc;

  sc  = NULL;
  e   = vrna_E_hp_loop(vc, i, j);

  if(e == en){
    switch(vc->type){
      case  VRNA_VC_TYPE_SINGLE:    sc  = vc->sc;
                                    break;
      case  VRNA_VC_TYPE_ALIGNMENT: if(vc->scs)
                                      sc = vc->scs[0];
                                    break;
      default:                      break;
    }

    if(sc)
      if(sc->bt){
        vrna_basepair_t *ptr, *aux_bps;
        aux_bps = sc->bt(i, j, -1, -1, VRNA_DECOMP_PAIR_HP, sc->data);
        for(ptr = aux_bps; ptr && ptr->i != -1; ptr++){
          bp_stack[++(*stack_count)].i = ptr->i;
          bp_stack[(*stack_count)].j   = ptr->j;
        }
        free(aux_bps);
      }

    return 1;
  }

  return 0;
}

/**
 * @}
 */


#endif
