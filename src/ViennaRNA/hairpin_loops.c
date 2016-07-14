
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


PRIVATE FLT_OR_DBL exp_eval_hp_loop(vrna_fold_compound_t *vc, int i, int j);
PRIVATE FLT_OR_DBL exp_eval_ext_hp_loop(vrna_fold_compound_t *vc, int i, int j);

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

  int                       p, q, u, *hc_up;
  char                      eval_loop;
  vrna_hc_t                 *hc;
  int (*eval_f)(vrna_fold_compound_t *a, int b, int c);

  hc    = vc->hc;
  hc_up = hc->up_hp;

#ifdef WITH_GEN_HC
  vrna_callback_hc_evaluate *f = hc->f;
#endif

  if((i > 0) && (j > 0)){
    if(j > i){ /* linear case */
      p       = i;
      q       = j;
      u       = q - p - 1;
      eval_f  = vrna_eval_hp_loop;
    } else { /* circular case */
      p       = j;
      q       = i;
      u       = vc->length - q + p - 1;
      eval_f  = vrna_eval_ext_hp_loop;
    }

    eval_loop = hc->matrix[vc->jindx[q] + p] & VRNA_CONSTRAINT_CONTEXT_HP_LOOP;

#ifdef WITH_GEN_HC
    if(f)
      eval_loop = (f(i, j, i, j, VRNA_DECOMP_PAIR_HP, hc->data)) ? eval_loop : (char)0;
#endif

    /* is this base pair allowed to close a hairpin (like) loop ? */
    if(eval_loop){
      
      /* are all nucleotides in the loop allowed to be unpaired ? */
      if(hc_up[i+1] >= u){
        return eval_f(vc, p, q);
      }
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

  return vrna_E_hp_loop(vc, j, i);
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

  int             u, e, s, type, *types, n_seq, length;
  short           *S, **SS, **S5, **S3;
  char            **Ss;
  unsigned short  **a2s;
  vrna_param_t    *P;
  vrna_sc_t       *sc, **scs;
  vrna_md_t       *md;
  char            loopseq[10];

  length  = vc->length;
  P       = vc->params;
  md      = &(P->model_details);
  e       = INF;

  switch(vc->type){
    /* single sequences and cofolding hybrids */
    case  VRNA_VC_TYPE_SINGLE:    S     = vc->sequence_encoding;
                                  sc    = vc->sc;
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
                                  break;

    /* sequence alignments */
    case  VRNA_VC_TYPE_ALIGNMENT: SS    = vc->S;                                                               
                                  S5    = vc->S5;     /*S5[s][i] holds next base 5' of i in sequence s*/
                                  S3    = vc->S3;     /*Sl[s][i] holds next base 3' of i in sequence s*/
                                  Ss    = vc->Ss;                                                       
                                  a2s   = vc->a2s;                                                      
                                  scs   = vc->scs;
                                  n_seq = vc->n_seq;
                                  e     = 0;
                                  types = (int *)vrna_alloc(sizeof(int) * n_seq);

                                  for (s=0; s<n_seq; s++) {
                                    types[s] = md->pair[SS[s][j]][SS[s][i]];
                                    if (types[s]==0) types[s]=7;
                                  }
                                  
                                  for (s=0; s<n_seq; s++) {
                                    char loopseq[10];
                                    u = a2s[s][length] - a2s[s][j] + a2s[s][i - 1];

                                    if (u<9) {
                                      strcpy(loopseq , Ss[s] + a2s[s][j] - 1);
                                      strncat(loopseq, Ss[s], a2s[s][i]);
                                    }
                                    if (u < 3) e += 600;
                                    else e += E_Hairpin(u, types[s], S3[s][j], S5[s][i],  loopseq, P);
                                  }
                                  if(scs)
                                    for(s=0;s < n_seq; s++){
                                      if(scs[s]){
                                        if(scs[s]->energy_up){
                                          e +=    ((i > 1) ? scs[s]->energy_up[1][a2s[s][i - 1]] : 0)
                                                    + ((j < length) ? scs[s]->energy_up[a2s[s][j + 1]][a2s[s][length] - a2s[s][j]] : 0);
                                        }
                                        if(scs[s]->f)
                                          e += scs[s]->f(a2s[s][j], a2s[s][i], a2s[s][j], a2s[s][i], VRNA_DECOMP_PAIR_HP, scs[s]->data);
                                      }
                                    }

                                  free(types);
                                  break;

    /* nothing */
    default:                      break;
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
                                  ij    = idx[j] + i;
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
                                          e += scs[s]->energy_up[a2s[s][i + 1]][u];

                                        if(scs[s]->energy_bp)
                                          e += scs[s]->energy_bp[ij];

                                        if(scs[s]->f)
                                          e += scs[s]->f(a2s[s][i], a2s[s][j], a2s[s][i], a2s[s][j], VRNA_DECOMP_PAIR_HP, scs[s]->data);
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
PUBLIC FLT_OR_DBL
vrna_exp_E_hp_loop( vrna_fold_compound_t *vc,
                    int i,
                    int j){

  int         p, q, u, *hc_up;
  char        eval_loop;
  FLT_OR_DBL  z;
  vrna_hc_t   *hc;
  FLT_OR_DBL  (*eval_f)(vrna_fold_compound_t *a, int b, int c);
#ifdef WITH_GEN_HC
  vrna_callback_hc_evaluate *f = hc->f;
#endif

  z = 0.;

  hc        = vc->hc;
  hc_up     = hc->up_hp;

  if((i > 0) && (j > 0)){
    if(j >= i){ /* linear case */
      p       = i;
      q       = j;
      u       = q - p - 1;
      eval_f  = exp_eval_hp_loop;
    } else { /* circular case */
      p       = j;
      q       = i;
      u       = vc->length - q + p - 1;
      eval_f  = exp_eval_ext_hp_loop;
    }

    eval_loop = hc->matrix[vc->jindx[q] + p] & VRNA_CONSTRAINT_CONTEXT_HP_LOOP;

#ifdef WITH_GEN_HC
    if(f)
      eval_loop = (f(i, j, i, j, VRNA_DECOMP_PAIR_HP, hc->data)) ? eval_loop : (char)0;
#endif

    /* is this base pair allowed to close a hairpin (like) loop ? */
    if(eval_loop){
      /* are all nucleotides in the loop allowed to be unpaired ? */
      if(hc_up[i+1] >= u){
        z = eval_f(vc, p, q);
      }
    }
  }

  return z;
}

PRIVATE FLT_OR_DBL
exp_eval_hp_loop( vrna_fold_compound_t *vc,
                  int i,
                  int j){

  int               u, ij, type, n_seq, s, *types, cp, *idx, *iidx;
  FLT_OR_DBL        q, qbt1;
  FLT_OR_DBL        *scale;
  short             *S, **SS, **S5, **S3;
  char              **Ss;
  unsigned short    **a2s;
  vrna_exp_param_t  *P;
  vrna_sc_t         *sc, **scs;
  vrna_md_t         *md;

  cp    = vc->cutpoint;
  idx   = vc->jindx;
  iidx  = vc->iindx;
  P     = vc->exp_params;
  md    = &(P->model_details);
  scale = vc->exp_matrices->scale;
  types = NULL;

  q     = 0.;
  ij    = idx[j] + i;

  switch(vc->type){
    case VRNA_VC_TYPE_SINGLE:     S     = vc->sequence_encoding;
                                  sc    = vc->sc;
                                  u     = j - i - 1;
                                  type  = vc->ptype[ij];

                                  if((cp < 0) || ON_SAME_STRAND(i, j, cp)){ /* regular hairpin loop */
                                    q = exp_E_Hairpin(u, type, S[i+1], S[j-1], vc->sequence+i-1, P);
                                  } else { /* hairpin-like exterior loop (for cofolding) */
                                    /* this is currently handle somewhere else */
                                  }

                                  /* add soft constraints */
                                  if(sc){
                                    if(sc->exp_energy_up)
                                      q *= sc->exp_energy_up[i+1][u];

                                    if(sc->exp_energy_bp)
                                      q *= sc->exp_energy_bp[iidx[i] - j];

                                    if(sc->exp_f)
                                      q *= sc->exp_f(i, j, i, j, VRNA_DECOMP_PAIR_HP, sc->data);
                                  }

                                  q *= scale[u+2];
                                  break;

    case VRNA_VC_TYPE_ALIGNMENT:  SS    = vc->S;                                                               
                                  S5    = vc->S5;     /*S5[s][i] holds next base 5' of i in sequence s*/
                                  S3    = vc->S3;     /*Sl[s][i] holds next base 3' of i in sequence s*/
                                  Ss    = vc->Ss;                                                       
                                  a2s   = vc->a2s;                                                      
                                  scs   = vc->scs;
                                  n_seq = vc->n_seq;
                                  qbt1  = 1.;
                                  types = (int *)vrna_alloc(sizeof(int) * n_seq);

                                  for (s=0; s<n_seq; s++) {
                                    types[s] = md->pair[SS[s][i]][SS[s][j]];
                                    if (types[s]==0) types[s]=7;
                                  }
                                  
                                  for (s=0; s<n_seq; s++) {
                                    u = a2s[s][j-1] - a2s[s][i];
                                    if (a2s[s][i]<1) continue;
                                    char loopseq[10];
                                    if (u<9){
                                      strncpy(loopseq, Ss[s]+a2s[s][i]-1, 10);
                                    }
                                    qbt1 *= exp_E_Hairpin(u, types[s], S3[s][i], S5[s][j], loopseq, P);
                                  }

                                  /* add soft constraints */
                                  if(scs)
                                    for(s = 0; s < n_seq; s++){
                                      if(scs[s]){
                                        u = a2s[s][j-1] - a2s[s][i];

                                        if(scs[s]->exp_energy_bp)
                                          qbt1 *= scs[s]->exp_energy_bp[iidx[i] - j];

                                        if(scs[s]->exp_energy_up)
                                          qbt1 *= scs[s]->exp_energy_up[a2s[s][i + 1]][u];

                                        if(scs[s]->exp_f)
                                          qbt1 *= scs[s]->exp_f(a2s[s][i], a2s[s][j], a2s[s][i], a2s[s][j], VRNA_DECOMP_PAIR_HP, scs[s]->data);
                                      }
                                    }

                                  q = qbt1 * scale[j-i+1];
                                  break;

    default:                      break;
  }

  free(types);
  return q;
}

PRIVATE FLT_OR_DBL
exp_eval_ext_hp_loop( vrna_fold_compound_t *vc,
                      int i,
                      int j){

  int               u, u1, ij, n, type, n_seq, s, *rtype, *types, *idx, *iidx, no_close, noGUclosure;
  FLT_OR_DBL        q, qbt1;
  FLT_OR_DBL        *scale;
  short             *S, **SS, **S5, **S3;
  char              **Ss, *sequence;
  unsigned short    **a2s;
  vrna_exp_param_t  *P;
  vrna_sc_t         *sc, **scs;
  vrna_md_t         *md;

  n           = vc->length;
  idx         = vc->jindx;
  iidx        = vc->iindx;
  P           = vc->exp_params;
  md          = &(P->model_details);
  noGUclosure = md->noGUclosure;
  scale       = vc->exp_matrices->scale;
  types       = NULL;
  rtype       = &(md->rtype[0]);

  q     = 0.;
  u     = n - j + i - 1;
  ij    = idx[j] + i;

  switch(vc->type){
    case VRNA_VC_TYPE_SINGLE:     sequence  = vc->sequence;
                                  S         = vc->sequence_encoding;
                                  sc        = vc->sc;
                                  type      = rtype[vc->ptype[ij]];

                                  if(((type==3)||(type==4))&&noGUclosure)
                                    return q;

                                  /* get the loop sequence */
                                  char loopseq[10];
                                  if (u<7){
                                    strcpy(loopseq , sequence+j-1);
                                    strncat(loopseq, sequence, i);
                                  }

                                  q = exp_E_Hairpin(u, type, S[j+1], S[i-1], loopseq, P);

                                  /* add soft constraints */
                                  if(sc){
                                    if(sc->exp_energy_up)
                                      q *=    ((i > 1) ? sc->exp_energy_up[1][i - 1] : 1.)
                                            * ((j < n) ? sc->exp_energy_up[j+1][n-j] : 1.);

                                    if(sc->exp_f)
                                      q *= sc->exp_f(j, i, j, i, VRNA_DECOMP_PAIR_HP, sc->data);
                                  }

                                  q *= scale[u];
                                  break;

    case VRNA_VC_TYPE_ALIGNMENT:  SS    = vc->S;                                                               
                                  S5    = vc->S5;     /*S5[s][i] holds next base 5' of i in sequence s*/
                                  S3    = vc->S3;     /*Sl[s][i] holds next base 3' of i in sequence s*/
                                  Ss    = vc->Ss;                                                       
                                  a2s   = vc->a2s;                                                      
                                  scs   = vc->scs;
                                  n_seq = vc->n_seq;
                                  qbt1  = 1.;
                                  types = (int *)vrna_alloc(sizeof(int) * n_seq);

                                  for (s=0; s<n_seq; s++) {
                                    types[s] = md->pair[SS[s][j]][SS[s][i]];
                                    if (types[s]==0) types[s]=7;
                                  }
                                  
                                  for (s=0; s<n_seq; s++) {
                                    u1 = a2s[s][i] - 1 + a2s[s][n] - a2s[s][j];
                                    char loopseq[10];
                                    if (u1<7){
                                      strcpy(loopseq , Ss[s] + a2s[s][j] - 1);
                                      strncat(loopseq, Ss[s], a2s[s][i]);
                                    }
                                    qbt1 *= exp_E_Hairpin(u1, types[s], S3[s][j], S5[s][i], loopseq, P);
                                  }

                                  /* add soft constraints */
                                  if(scs)
                                    for(s = 0; s < n_seq; s++){
                                      if(scs[s]){

                                        if(scs[s]->exp_energy_up)
                                          qbt1 *=   ((i > 1) ? scs[s]->exp_energy_up[a2s[s][1]][a2s[s][i] - a2s[s][1]] : 1.)
                                                  * ((j < n) ? scs[s]->exp_energy_up[a2s[s][j] + 1][a2s[s][n] - a2s[s][j]] : 1.);

                                        if(scs[s]->exp_f)
                                          qbt1 *= scs[s]->exp_f(a2s[s][j], a2s[s][i], a2s[s][j], a2s[s][i], VRNA_DECOMP_PAIR_HP, scs[s]->data);
                                      }
                                    }

                                  q = qbt1 * scale[u];

                                  free(types);
                                  break;

    default:                      break;
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

