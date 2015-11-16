
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <string.h>
#include "ViennaRNA/fold_vars.h"
#include "ViennaRNA/energy_par.h"
#include "ViennaRNA/data_structures.h"
#include "ViennaRNA/params.h"
#include "ViennaRNA/utils.h"
#include "ViennaRNA/constraints.h"
#include "ViennaRNA/exterior_loops.h"
#include "ViennaRNA/gquad.h"
#include "ViennaRNA/hairpin_loops.h"

#ifdef ON_SAME_STRAND
#undef ON_SAME_STRAND
#endif

#define ON_SAME_STRAND(I,J,C)  (((I)>=(C))||((J)<(C)))

/*
#################################
# BEGIN OF FUNCTION DEFINITIONS #
#################################
*/

/**
 *  @brief  Evaluate the free energy of a hairpin loop
 *          and consider possible hard constraints
 *
 *  @note This function is polymorphic! The provided #vrna_fold_compound_t may be of type
 *  #VRNA_VC_TYPE_SINGLE or #VRNA_VC_TYPE_ALIGNMENT
 *
 */
PUBLIC int
vrna_E_hp_loop( vrna_fold_compound_t *vc,
                int i,
                int j){

  int                       u, *hc_up;
  char                      eval_loop;
  vrna_callback_hc_evaluate *f;
  vrna_hc_t                 *hc;

  u         = j - i - 1;
  hc        = vc->hc;
  hc_up     = hc->up_hp;
  f         = hc->f;
  eval_loop = hc->matrix[vc->jindx[j] + i] & VRNA_CONSTRAINT_CONTEXT_HP_LOOP;

#ifdef WITH_GEN_HC
  if(f)
    eval_loop = (f(i, j, i, j, VRNA_DECOMP_PAIR_HP, hc->data)) ? eval_loop : (char)0;
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
PUBLIC int
vrna_E_ext_hp_loop( vrna_fold_compound_t *vc,
                    int i,
                    int j){

  int                       u, *hc_up;
  char                      eval_loop;
  vrna_callback_hc_evaluate *f;
  vrna_hc_t                 *hc;

  u         = vc->length - j + i - 1;
  hc        = vc->hc;
  hc_up     = hc->up_hp;
  f         = hc->f;
  eval_loop = hc->matrix[vc->jindx[j] + i] & VRNA_CONSTRAINT_CONTEXT_HP_LOOP;

#ifdef WITH_GEN_HC
  if(f)
    eval_loop = (f(j, i, j, i, VRNA_DECOMP_PAIR_HP, hc->data)) ? eval_loop : (char)0;
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
PUBLIC int
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
    if(sc->energy_up)
      e +=  sc->energy_up[j + 1][vc->length - j]
            + sc->energy_up[1][i - 1];

    if(sc->f)
      e += sc->f(j, i, j, i, VRNA_DECOMP_PAIR_HP, sc->data);
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
PUBLIC int
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
                                    if(sc->energy_up)
                                      e += sc->energy_up[i+1][u];

                                    if(sc->energy_bp){
                                      e += sc->energy_bp[ij];
                                    }
                                    if(sc->f)
                                      e += sc->f(i, j, i, j, VRNA_DECOMP_PAIR_HP, sc->data);
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

                                        if(scs[s]->energy_up)
                                          e += scs[s]->energy_up[a2s[s][i]+1][u];

                                        if(scs[s]->energy_bp)
                                          e += scs[s]->energy_bp[idx[j] + i];

                                        if(scs[s]->f)
                                          e += scs[s]->f(i, j, i, j, VRNA_DECOMP_PAIR_HP, scs[s]->data);
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


/**
 *  @brief High-Level function for hairpin loop energy evaluation (partition function variant)
 *
 *  @see E_hp_loop() for it's free energy counterpart
*/
PUBLIC double
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
        if(sc->exp_energy_up)
          q *= sc->exp_energy_up[i+1][u];

        if(sc->exp_energy_bp)
          q *= sc->exp_energy_bp[ij];

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
PUBLIC int
vrna_BT_hp_loop(vrna_fold_compound_t *vc,
                int i,
                int j,
                int en,
                vrna_bp_stack_t *bp_stack,
                int   *stack_count){

  int       e, u;
  vrna_sc_t *sc;

  sc  = NULL;

  u = j - i - 1;

  if(vc->hc->up_hp[i+1] < u)
    return 0;

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
        aux_bps = sc->bt(i, j, i, j, VRNA_DECOMP_PAIR_HP, sc->data);
        for(ptr = aux_bps; ptr && ptr->i != 0; ptr++){
          bp_stack[++(*stack_count)].i = ptr->i;
          bp_stack[(*stack_count)].j   = ptr->j;
        }
        free(aux_bps);
      }

    return 1;
  }

  return 0;
}

