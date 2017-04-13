/** \file **/

/*
                  minimum free energy
                  RNA secondary structure prediction

                  c Ivo Hofacker, Chrisoph Flamm
                  original implementation by
                  Walter Fontana
                  g-quadruplex support and threadsafety
                  by Ronny Lorenz

                  Vienna RNA package
*/

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <string.h>
#include <limits.h>

#include "ViennaRNA/utils.h"
#include "ViennaRNA/energy_par.h"
#include "ViennaRNA/data_structures.h"
#include "ViennaRNA/fold_vars.h"
#include "ViennaRNA/params.h"
#include "ViennaRNA/constraints.h"
#include "ViennaRNA/gquad.h"
#include "ViennaRNA/structured_domains.h"
#include "ViennaRNA/unstructured_domains.h"
#include "ViennaRNA/loop_energies.h"
#include "ViennaRNA/mfe.h"

/* make this interface backward compatible with RNAlib < 2.2.0 */
#define VRNA_BACKWARD_COMPAT

#define MAXSECTORS        500     /* dimension for a backtrack array */

/*
#################################
# GLOBAL VARIABLES              #
#################################
*/

/*
#################################
# PRIVATE VARIABLES             #
#################################
*/

/*
#################################
# PRIVATE FUNCTION DECLARATIONS #
#################################
*/

PRIVATE int           fill_arrays(vrna_fold_compound_t *vc);
PRIVATE void          fill_arrays_circ(vrna_fold_compound_t *vc, sect bt_stack[], int *bt);
PRIVATE void          backtrack(vrna_fold_compound_t *vc, vrna_bp_stack_t *bp_stack, sect bt_stack[], int s);

PRIVATE int           fill_arrays_comparative(vrna_fold_compound_t *vc);
PRIVATE void          fill_arrays_comparative_circ(vrna_fold_compound_t *vc, sect bt_stack[], int *bt);
PRIVATE void          backtrack_comparative(vrna_fold_compound_t *vc, vrna_bp_stack_t *bp_stack, sect bt_stack[], int s);

/*
#################################
# BEGIN OF FUNCTION DEFINITIONS #
#################################
*/

PUBLIC float
vrna_mfe( vrna_fold_compound_t *vc,
          char *structure){

  char    *ss;
  int     length, energy, s;
  float   mfe;
  sect    bt_stack[MAXSECTORS]; /* stack of partial structures for backtracking */
  vrna_bp_stack_t   *bp;

  s       = 0;
  mfe     = (float)(INF/100.);

  if(vc){
    length  = (int) vc->length;

    if(!vrna_fold_compound_prepare(vc, VRNA_OPTION_MFE)){
      vrna_message_warning("vrna_mfe@mfe.c: Failed to prepare vrna_fold_compound");
      return mfe;
    }

    /* call user-defined recursion status callback function */
    if(vc->stat_cb)
      vc->stat_cb(VRNA_STATUS_MFE_PRE, vc->auxdata);

    switch(vc->type){
      case VRNA_FC_TYPE_SINGLE:     energy = fill_arrays(vc);
                                    if(vc->params->model_details.circ){
                                      fill_arrays_circ(vc, bt_stack, &s);
                                      energy = vc->matrices->Fc;
                                    }
                                    break;

      case VRNA_FC_TYPE_COMPARATIVE:  energy = fill_arrays_comparative(vc);
                                    if(vc->params->model_details.circ){
                                      fill_arrays_comparative_circ(vc, bt_stack, &s);
                                      energy = vc->matrices->Fc;
                                    }
                                    break;

      default:                      vrna_message_warning("unrecognized fold compound type");
                                    return mfe;
                                    break;
    }


    /* call user-defined recursion status callback function */
    if(vc->stat_cb)
      vc->stat_cb(VRNA_STATUS_MFE_POST, vc->auxdata);

    if(structure && vc->params->model_details.backtrack){
      bp = (vrna_bp_stack_t *)vrna_alloc(sizeof(vrna_bp_stack_t) * (4*(1+length/2))); /* add a guess of how many G's may be involved in a G quadruplex */

      switch(vc->type){
        case VRNA_FC_TYPE_COMPARATIVE:  backtrack_comparative(vc, bp, bt_stack, s);
                                      break;

        case VRNA_FC_TYPE_SINGLE:     /* fall through */

        default:                      backtrack(vc, bp, bt_stack, s);
                                      break;
      }

      ss = vrna_db_from_bp_stack(bp, length);
      strncpy(structure, ss, length + 1);
      free(ss);
      free(bp);
    }

    if (vc->params->model_details.backtrack_type=='C')
      mfe = (float) vc->matrices->c[vc->jindx[length]+1]/100.;
    else if (vc->params->model_details.backtrack_type=='M')
      mfe = (float) vc->matrices->fML[vc->jindx[length]+1]/100.;
    else
      mfe = (float) energy/100.;

    if(vc->type == VRNA_FC_TYPE_COMPARATIVE)
      mfe /= (float)vc->n_seq;
  }

  return mfe;
}

/**
*** fill "c", "fML" and "f5" arrays and return  optimal energy
**/
PRIVATE int
fill_arrays(vrna_fold_compound_t *vc){

  unsigned char     type;
  char              *ptype, *hard_constraints;
  int               i, j, ij, length, energy, new_c, stackEnergy, no_close, turn,
                    noGUclosure, noLP, uniq_ML, dangle_model, *indx, *my_f5,
                    *my_c, *my_fML, *my_fM1, hc_decompose, *cc, *cc1, *Fmi, *DMLi,
                    *DMLi1, *DMLi2;
  vrna_param_t      *P;
  vrna_mx_mfe_t     *matrices;
  vrna_hc_t         *hc;
  vrna_ud_t         *domains_up;

  length            = (int)vc->length;
  ptype             = vc->ptype;
  indx              = vc->jindx;
  P                 = vc->params;
  noGUclosure       = P->model_details.noGUclosure;
  noLP              = P->model_details.noLP;
  uniq_ML           = P->model_details.uniq_ML;
  dangle_model      = P->model_details.dangles;
  turn              = P->model_details.min_loop_size;
  hc                = vc->hc;
  hard_constraints  = hc->matrix;
  matrices          = vc->matrices;
  my_f5             = matrices->f5;
  my_c              = matrices->c;
  my_fML            = matrices->fML;
  my_fM1            = matrices->fM1;
  domains_up        = vc->domains_up;

  /* allocate memory for all helper arrays */
  cc    = (int *) vrna_alloc(sizeof(int)*(length + 2)); /* auxilary arrays for canonical structures     */
  cc1   = (int *) vrna_alloc(sizeof(int)*(length + 2)); /* auxilary arrays for canonical structures     */
  Fmi   = (int *) vrna_alloc(sizeof(int)*(length + 1)); /* holds row i of fML (avoids jumps in memory)  */
  DMLi  = (int *) vrna_alloc(sizeof(int)*(length + 1)); /* DMLi[j] holds  MIN(fML[i,k]+fML[k+1,j])      */
  DMLi1 = (int *) vrna_alloc(sizeof(int)*(length + 1)); /*                MIN(fML[i+1,k]+fML[k+1,j])    */
  DMLi2 = (int *) vrna_alloc(sizeof(int)*(length + 1)); /*                MIN(fML[i+2,k]+fML[k+1,j])    */

  if((turn < 0) || (turn > length))
    turn = length; /* does this make any sense? */

  /* pre-processing ligand binding production rule(s) */
  if(domains_up && domains_up->prod_cb)
    domains_up->prod_cb(vc, domains_up->data);

  /* prefill helper arrays */
  for(j = 0; j <= length; j++){
    Fmi[j] = DMLi[j] = DMLi1[j] = DMLi2[j] = INF;
  }


  /* prefill matrices with init contributions */
  for(j = 1; j <= length; j++)
    for(i = (j > turn ? (j - turn) : 1); i <= j; i++){
      my_c[indx[j] + i] = my_fML[indx[j] + i] = INF;
      if(uniq_ML)
        my_fM1[indx[j] + i] = INF;
    }

  /* start recursion */

  if (length <= turn){
    /* clean up memory */
    free(cc);
    free(cc1);
    free(Fmi);
    free(DMLi);
    free(DMLi1);
    free(DMLi2);
    /* return free energy of unfolded chain */
    return 0;
  }

  for (i = length-turn-1; i >= 1; i--) { /* i,j in [1..length] */

    for (j = i+turn+1; j <= length; j++) {
      ij            = indx[j]+i;
      type          = (unsigned char)ptype[ij];
      hc_decompose  = hard_constraints[ij];
      energy        = INF;

      no_close = (((type==3)||(type==4))&&noGUclosure);

      if (hc_decompose) {   /* we evaluate this pair */
        new_c = INF;

        if(!no_close){
          /* check for hairpin loop */
          energy = vrna_E_hp_loop(vc, i, j);
          new_c = MIN2(new_c, energy);

          /* check for multibranch loops */
          energy  = vrna_E_mb_loop_fast(vc, i, j, DMLi1, DMLi2);
          new_c   = MIN2(new_c, energy);
        }

        if(dangle_model == 3){ /* coaxial stacking */
          energy  = E_mb_loop_stack(i, j, vc);
          new_c   = MIN2(new_c, energy);
        }

        /* check for interior loops */
        energy = vrna_E_int_loop(vc, i, j);
        new_c = MIN2(new_c, energy);

        /* remember stack energy for --noLP option */
        if(noLP){
          stackEnergy = vrna_E_stack(vc, i, j);
          new_c       = MIN2(new_c, cc1[j-1]+stackEnergy);
          cc[j]       = new_c;
          my_c[ij]    = cc1[j-1]+stackEnergy;
        } else {
          my_c[ij]    = new_c;
        }
      } /* end >> if (pair) << */

      else my_c[ij] = INF;

      /* done with c[i,j], now compute fML[i,j] and fM1[i,j] */

      my_fML[ij] = vrna_E_ml_stems_fast(vc, i, j, Fmi, DMLi);

      if(uniq_ML){  /* compute fM1 for unique decomposition */
        my_fM1[ij] = E_ml_rightmost_stem(i, j, vc);
      }

    } /* end of j-loop */

    {
      int *FF; /* rotate the auxilliary arrays */
      FF = DMLi2; DMLi2 = DMLi1; DMLi1 = DMLi; DMLi = FF;
      FF = cc1; cc1=cc; cc=FF;
      for (j=1; j<=length; j++) {cc[j]=Fmi[j]=DMLi[j]=INF; }
    }
  } /* end of i-loop */

  /* calculate energies of 5' fragments */
  E_ext_loop_5(vc);

  /* clean up memory */
  free(cc);
  free(cc1);
  free(Fmi);
  free(DMLi);
  free(DMLi1);
  free(DMLi2);

  return my_f5[length];
}

#include "circfold.inc"


/**
*** the actual forward recursion to fill the energy arrays
**/
PRIVATE int
fill_arrays_comparative(vrna_fold_compound_t *vc){

  char              *hard_constraints;
  unsigned short    **a2s;
  short             **S, **S5, **S3;
  int               i, j, turn, energy, stackEnergy, new_c, s, *type, tt, *cc,
                    *cc1, *Fmi, *DMLi, *DMLi1, *DMLi2, n_seq, length, *indx,
                    *c, *f5, *fML, *ggg, *pscore, dangle_model;
  vrna_param_t      *P;
  vrna_md_t         *md;
  vrna_hc_t         *hc;
  vrna_sc_t         **sc;

  n_seq             = vc->n_seq;
  length            = vc->length;
  S                 = vc->S;
  S5                = vc->S5;             /* S5[s][i] holds next base 5' of i in sequence s */
  S3                = vc->S3;             /* Sl[s][i] holds next base 3' of i in sequence s */
  a2s               = vc->a2s;
  P                 = vc->params;
  md                = &(P->model_details);
  indx              = vc->jindx;          /* index for moving in the triangle matrices c[] and fMl[] */
  c                 = vc->matrices->c;    /* energy array, given that i-j pair */
  f5                = vc->matrices->f5;   /* energy of 5' end */
  fML               = vc->matrices->fML;  /* multi-loop auxiliary energy array */
  ggg               = vc->matrices->ggg;
  pscore            = vc->pscore;         /* precomputed array of pair types */
  dangle_model      = md->dangles;
  turn              = md->min_loop_size;
  hc                = vc->hc;
  sc                = vc->scs;
  hard_constraints  = hc->matrix;

  /* allocate some memory for helper arrays */
  type  = (int *) vrna_alloc(n_seq*sizeof(int));
  cc    = (int *) vrna_alloc(sizeof(int)*(length+2)); /* linear array for calculating canonical structures */
  cc1   = (int *) vrna_alloc(sizeof(int)*(length+2)); /*   "     "        */
  Fmi   = (int *) vrna_alloc(sizeof(int)*(length+1)); /* holds row i of fML (avoids jumps in memory) */
  DMLi  = (int *) vrna_alloc(sizeof(int)*(length+1)); /* DMLi[j] holds MIN(fML[i,k]+fML[k+1,j])  */
  DMLi1 = (int *) vrna_alloc(sizeof(int)*(length+1)); /*             MIN(fML[i+1,k]+fML[k+1,j])  */
  DMLi2 = (int *) vrna_alloc(sizeof(int)*(length+1)); /*             MIN(fML[i+2,k]+fML[k+1,j])  */


  if((turn < 0) || (turn > length))
    turn = length;

  /* init energies */
  for (j=1; j<=length; j++){
    Fmi[j]=DMLi[j]=DMLi1[j]=DMLi2[j]=INF;
    for (i=(j>turn?(j-turn):1); i<j; i++) {
      c[indx[j]+i] = fML[indx[j]+i] = INF;
    }
  }

  /* begin recursions */
  for (i = length-turn-1; i >= 1; i--) { /* i,j in [1..length] */
    for (j = i+turn+1; j <= length; j++) {
      int ij, psc;
      ij = indx[j]+i;

      for (s=0; s<n_seq; s++) {
        type[s] = md->pair[S[s][i]][S[s][j]];
        if (type[s]==0) type[s]=7;
      }

      psc = pscore[indx[j]+i];
      if (hard_constraints[ij]) {   /* a pair to consider */
        new_c = INF;

        /* hairpin ----------------------------------------------*/
        energy  = vrna_E_hp_loop(vc, i, j);
        new_c   = MIN2(new_c, energy);

        /* check for multibranch loops */
        energy  = vrna_E_mb_loop_fast(vc, i, j, DMLi1, DMLi2);
        new_c   = MIN2(new_c, energy);

        /* check for interior loops */
        energy  = vrna_E_int_loop(vc, i, j);
        new_c   = MIN2(new_c, energy);

        /* remember stack energy for --noLP option */
        if(md->noLP){
          stackEnergy = vrna_E_stack(vc, i, j);
          new_c       = MIN2(new_c, cc1[j-1]+stackEnergy);
          cc[j]       = new_c - psc; /* add covariance bonnus/penalty */
          c[ij]       = cc1[j-1] + stackEnergy - psc;
        } else {
          c[ij]       = new_c - psc; /* add covariance bonnus/penalty */
        }
      } /* end >> if (pair) << */

      else c[ij] = INF;

      /* done with c[i,j], now compute fML[i,j] */
      fML[ij] = vrna_E_ml_stems_fast(vc, i, j, Fmi, DMLi);

    } /* END for j */

    {
      int *FF; /* rotate the auxilliary arrays */
      FF = DMLi2; DMLi2 = DMLi1; DMLi1 = DMLi; DMLi = FF;
      FF = cc1; cc1=cc; cc=FF;
      for (j=1; j<=length; j++) {cc[j]=Fmi[j]=DMLi[j]=INF; }
    }
  } /* END for i */
  /* calculate energies of 5' and 3' fragments */

  f5[0] = 0;
  for(j = 1; j <= turn + 1; j++){
    if(hc->up_ext[j]){
      energy = f5[j-1];
      if((energy < INF) && sc)
        for(s=0;s < n_seq; s++){
          if(sc[s]){
            if(sc[s]->energy_up)
              energy += sc[s]->energy_up[a2s[s][j]][1];
          }
        }
    } else {
      energy = INF;
    }
    f5[j] = energy;
  }

  switch(dangle_model){
    case 0:   for(j = turn + 2; j <= length; j++){
                f5[j] = INF;

                if(hc->up_ext[j]){
                  energy = f5[j-1];
                  if((energy < INF) && sc)
                    for(s=0; s < n_seq; s++){
                      if(sc[s]){
                        if(sc[s]->energy_up)
                          energy += sc[s]->energy_up[a2s[s][j]][1];
                      }
                    }
                  f5[j] = MIN2(f5[j], energy);
                }

                if (hard_constraints[indx[j]+1] & VRNA_CONSTRAINT_CONTEXT_EXT_LOOP){
                  if(c[indx[j]+1] < INF){
                    energy = c[indx[j]+1];
                    for(s = 0; s < n_seq; s++){
                      tt = md->pair[S[s][1]][S[s][j]];
                      if(tt==0) tt=7;
                      energy += E_ExtLoop(tt, -1, -1, P);
                    }
                    f5[j] = MIN2(f5[j], energy);
                  }
                }
                if(md->gquad){
                  if(ggg[indx[j]+1] < INF)
                    f5[j] = MIN2(f5[j], ggg[indx[j]+1]);
                }

                for(i = j - turn - 1; i > 1; i--){
                  if(hard_constraints[indx[j]+i] & VRNA_CONSTRAINT_CONTEXT_EXT_LOOP){
                    if(c[indx[j]+i]<INF){
                      energy = f5[i-1] + c[indx[j]+i];
                      for(s = 0; s < n_seq; s++){
                        tt = md->pair[S[s][i]][S[s][j]];
                        if(tt==0) tt=7;
                        energy += E_ExtLoop(tt, -1, -1, P);
                      }
                      f5[j] = MIN2(f5[j], energy);
                    }
                  }
                  if(md->gquad){
                    if(ggg[indx[j]+i] < INF)
                      f5[j] = MIN2(f5[j], f5[i-1] + ggg[indx[j]+i]);
                  }
                }
              }
              break;

    default:  for(j = turn + 2; j <= length; j++){
                f5[j] = INF;

                if(hc->up_ext[j]){
                  energy = f5[j-1];
                  if((energy < INF) && sc)
                    for(s=0; s < n_seq; s++){
                      if(sc[s]){
                        if(sc[s]->energy_up)
                          energy += sc[s]->energy_up[a2s[s][j]][1];
                      }
                    }
                  f5[j] = MIN2(f5[j], energy);
                }

                if(hard_constraints[indx[j]+1] & VRNA_CONSTRAINT_CONTEXT_EXT_LOOP){
                  if (c[indx[j]+1]<INF) {
                    energy = c[indx[j]+1];
                    for(s = 0; s < n_seq; s++){
                      tt = md->pair[S[s][1]][S[s][j]];
                      if(tt==0) tt=7;
                      energy += E_ExtLoop(tt, -1, (j<length) ? S3[s][j] : -1, P);
                    }
                    f5[j] = MIN2(f5[j], energy);
                  }
                }
                if(md->gquad){
                  if(ggg[indx[j]+1] < INF)
                    f5[j] = MIN2(f5[j], ggg[indx[j]+1]);
                }

                for(i = j - turn - 1; i > 1; i--){
                  if(hard_constraints[indx[j]+i] & VRNA_CONSTRAINT_CONTEXT_EXT_LOOP){
                    if (c[indx[j]+i]<INF) {
                      energy = f5[i-1] + c[indx[j]+i];
                      for(s = 0; s < n_seq; s++){
                        tt = md->pair[S[s][i]][S[s][j]];
                        if(tt==0) tt=7;
                        energy += E_ExtLoop(tt, S5[s][i], (j < length) ? S3[s][j] : -1, P);
                      }
                      f5[j] = MIN2(f5[j], energy);
                    }
                  }
                  if(md->gquad){
                    if(ggg[indx[j]+i] < INF)
                      f5[j] = MIN2(f5[j], f5[i-1] + ggg[indx[j]+i]);
                  }
                }
              }
              break;
  }
  free(type);
  free(cc);
  free(cc1);
  free(Fmi);
  free(DMLi);
  free(DMLi1);
  free(DMLi2);
  return(f5[length]);
}

#include "ViennaRNA/alicircfold.inc"

PUBLIC void
vrna_backtrack_from_intervals(vrna_fold_compound_t *vc,
                              vrna_bp_stack_t *bp_stack,
                              sect bt_stack[],
                              int s){

  if(vc){
    switch(vc->type){
      case VRNA_FC_TYPE_SINGLE:       backtrack(vc, bp_stack, bt_stack, s);
                                      break;

      case VRNA_FC_TYPE_COMPARATIVE:  backtrack_comparative(vc, bp_stack, bt_stack, s);
                                      break;
    }
  }
}

/**
*** trace back through the "c", "f5" and "fML" arrays to get the
*** base pairing list. No search for equivalent structures is done.
*** This is fast, since only few structure elements are recalculated.
***
*** normally s=0.
*** If s>0 then s items have been already pushed onto the bt_stack
**/
PRIVATE void
backtrack(vrna_fold_compound_t *vc,
          vrna_bp_stack_t *bp_stack,
          sect bt_stack[],
          int s){

  unsigned char   type;
  char            *string, *ptype, backtrack_type;
  int             i, j, ij, k, length, no_close, b, *my_c, *indx, noLP, noGUclosure;
  vrna_param_t    *P;

  b               = 0;
  length          = vc->length;
  my_c            = vc->matrices->c;
  indx            = vc->jindx;
  P               = vc->params;
  noLP            = P->model_details.noLP;
  noGUclosure     = P->model_details.noGUclosure;
  string          = vc->sequence;
  ptype           = vc->ptype;
  backtrack_type  = P->model_details.backtrack_type;

  if (s==0) {
    bt_stack[++s].i = 1;
    bt_stack[s].j = length;
    bt_stack[s].ml = (backtrack_type=='M') ? 1 : ((backtrack_type=='C')? 2: 0);
  }
  while (s>0) {
    int ml, cij;
    int canonical = 1;     /* (i,j) closes a canonical structure */

    /* pop one element from stack */
    i  = bt_stack[s].i;
    j  = bt_stack[s].j;
    ml = bt_stack[s--].ml;

    switch(ml){
      /* backtrack in f5 */
      case 0:   {
                  int p, q;
                  if(vrna_BT_ext_loop_f5(vc, &j, &p, &q, bp_stack, &b)){
                    if(j > 0){
                      bt_stack[++s].i = 1;
                      bt_stack[s].j   = j;
                      bt_stack[s].ml  = 0;
                    }
                    if(p > 0){
                      i = p;
                      j = q;
                      goto repeat1;
                    }

                    continue;
                  } else {
                    vrna_message_error("backtracking failed in f5 for sequence:\n%s\n", string);
                  }
                }
                break;

      /* trace back in fML array */
      case 1:   {
                  int p, q, comp1, comp2;
                  if(vrna_BT_mb_loop_split(vc, &i, &j, &p, &q, &comp1, &comp2, bp_stack, &b)){
                    if(i > 0){
                      bt_stack[++s].i = i;
                      bt_stack[s].j   = j;
                      bt_stack[s].ml  = comp1;
                    }
                    if(p > 0){
                      bt_stack[++s].i = p;
                      bt_stack[s].j   = q;
                      bt_stack[s].ml  = comp2;
                    }

                    continue;
                  } else {
                    vrna_message_error("backtracking failed in fML for sequence:\n%s\n", string);
                  }
                }
                break;

      /* backtrack in c */
      case 2:   bp_stack[++b].i = i;
                bp_stack[b].j   = j;
                goto repeat1;

      default:  vrna_message_error("Backtracking failed due to unrecognized DP matrix!");
                break;
    }

  repeat1:

    /*----- begin of "repeat:" -----*/
    ij = indx[j]+i;

    if (canonical)
      cij = my_c[ij];

    type = (unsigned char)ptype[ij];

    if (noLP)
      if(vrna_BT_stack(vc, &i, &j, &cij, bp_stack, &b)){
        canonical = 0;
        goto repeat1;
      }

    canonical = 1;

    no_close = (((type==3)||(type==4))&&noGUclosure);

    if (no_close) {
      if (cij == FORBIDDEN) continue;
    } else {
      if(vrna_BT_hp_loop(vc, i, j, cij, bp_stack, &b))
        continue;
    }

    if(vrna_BT_int_loop(vc, &i, &j, cij, bp_stack, &b)){
      if(i < 0)
        continue;
      else
        goto repeat1;
    }

    /* (i.j) must close a multi-loop */
    int comp1, comp2;

    if(vrna_BT_mb_loop(vc, &i, &j, &k, cij, &comp1, &comp2)){
      bt_stack[++s].i = i;
      bt_stack[s].j   = k;
      bt_stack[s].ml  = comp1;
      bt_stack[++s].i = k + 1;
      bt_stack[s].j   = j;
      bt_stack[s].ml  = comp2;
    } else {
      vrna_message_error("backtracking failed in repeat for sequence:\n%s\n", string);
    }

    /* end of repeat: --------------------------------------------------*/

  } /* end of infinite while loop */

  bp_stack[0].i = b;    /* save the total number of base pairs */
}


/**
*** backtrack in the energy matrices to obtain a structure with MFE
**/
PRIVATE void
backtrack_comparative(vrna_fold_compound_t *vc,
                      vrna_bp_stack_t *bp_stack,
                      sect bt_stack[],
                      int s) {

  /*------------------------------------------------------------------
    trace back through the "c", "f5" and "fML" arrays to get the
    base pairing list. No search for equivalent structures is done.
    This inverts the folding procedure, hence it's very fast.
    ------------------------------------------------------------------*/
   /* normally s=0.
     If s>0 then s items have been already pushed onto the sector stack */

  unsigned short  **a2s;
  short           **S, **S5, **S3, *S_cons;
  int             i, j, k, p, q, turn, energy, en, c0, l1, minq, maxq,
                  type_2, tt, mm, b, cov_en, *type, n_seq, length, *indx,
                  *c, *f5, *fML, *pscore, *ggg, *rtype, dangle_model, with_gquad;
  vrna_param_t    *P;
  vrna_md_t       *md;
  vrna_hc_t       *hc;
  vrna_sc_t       **sc;

  n_seq         = vc->n_seq;
  length        = vc->length;
  S             = vc->S;
  S5            = vc->S5;     /*S5[s][i] holds next base 5' of i in sequence s*/
  S3            = vc->S3;     /*Sl[s][i] holds next base 3' of i in sequence s*/
  a2s           = vc->a2s;
  P             = vc->params;
  md            = &(P->model_details);
  indx          = vc->jindx;     /* index for moving in the triangle matrices c[] and fMl[]*/
  c             = vc->matrices->c;     /* energy array, given that i-j pair */
  f5            = vc->matrices->f5;     /* energy of 5' end */
  fML           = vc->matrices->fML;     /* multi-loop auxiliary energy array */
  pscore        = vc->pscore;     /* precomputed array of pair types */
  ggg           = vc->matrices->ggg;
  S_cons        = vc->S_cons;
  rtype         = &(md->rtype[0]);
  dangle_model  = md->dangles;
  with_gquad    = md->gquad;
  hc            = vc->hc;
  sc            = vc->scs;
  turn          = md->min_loop_size;
  b             = 0;
  cov_en        = 0;

  type  = (int *) vrna_alloc(n_seq*sizeof(int));

  if((turn < 0) || (turn > length))
    turn = length;

  if (s==0) {
    bt_stack[++s].i = 1;
    bt_stack[s].j = length;
    bt_stack[s].ml = (md->backtrack_type=='M') ? 1 : ((md->backtrack_type=='C')?2:0);
  }
  while (s>0) {
    int ss, ml, fij, fi, cij, traced, i1, j1, jj=0, gq=0;
    int canonical = 1;     /* (i,j) closes a canonical structure */
    i  = bt_stack[s].i;
    j  = bt_stack[s].j;
    ml = bt_stack[s--].ml;   /* ml is a flag indicating if backtracking is to
                              occur in the fML- (1) or in the f-array (0) */
    if (ml==2) {
      bp_stack[++b].i = i;
      bp_stack[b].j   = j;
      cov_en += pscore[indx[j]+i];
      goto repeat1_comparative;
    }

    if (j < i+turn+1) continue; /* no more pairs in this interval */

    if(ml != 0){
      fij = fML[indx[j]+i];
      fi  = (hc->up_ml[j]) ? fML[indx[j-1]+i] + n_seq*P->MLbase : INF;
    } else {
      fij = f5[j];
      fi  = (hc->up_ext[j]) ? f5[j-1] : INF;
    }

    if(sc)
      for(ss = 0; ss < n_seq; ss++)
        if(sc[ss]){
          if(sc[ss]->energy_up)
            fi += sc[ss]->energy_up[a2s[ss][j]][1];
        }

    if (fij == fi) {  /* 3' end is unpaired */
      bt_stack[++s].i = i;
      bt_stack[s].j   = j-1;
      bt_stack[s].ml  = ml;
      continue;
    }

    if (ml == 0) { /* backtrack in f5 */
      switch(dangle_model){
        case 0:   /* j or j-1 is paired. Find pairing partner */
                  for (i=j-turn-1,traced=0; i>=1; i--) {
                    int en;
                    jj = i-1;

                    if (hc->matrix[indx[j] + i] & VRNA_CONSTRAINT_CONTEXT_EXT_LOOP){
                      en = c[indx[j]+i] + f5[i-1];
                      for(ss = 0; ss < n_seq; ss++){
                        type[ss] = md->pair[S[ss][i]][S[ss][j]];
                        if (type[ss]==0) type[ss] = 7;
                        en += E_ExtLoop(type[ss], -1, -1, P);
                      }
                      if (fij == en) traced=j;
                    }

                    if(with_gquad){
                      if(fij == f5[i-1] + ggg[indx[j]+i]){
                        /* found the decomposition */
                        traced = j; jj = i - 1; gq = 1;
                        break;
                      }
                    }

                    if (traced) break;
                  }
                  break;
        default:  /* j or j-1 is paired. Find pairing partner */
                  for (i=j-turn-1,traced=0; i>=1; i--) {
                    int en;
                    jj = i-1;
                    if (hc->matrix[indx[j] + i] & VRNA_CONSTRAINT_CONTEXT_EXT_LOOP){
                      en = c[indx[j]+i] + f5[i-1];
                      for(ss = 0; ss < n_seq; ss++){
                        type[ss] = md->pair[S[ss][i]][S[ss][j]];
                        if (type[ss]==0) type[ss] = 7;
                        en += E_ExtLoop(type[ss], (i>1) ? S5[ss][i]: -1, (j < length) ? S3[ss][j] : -1, P);
                      }
                      if (fij == en) traced=j;
                    }

                    if(with_gquad){
                      if(fij == f5[i-1] + ggg[indx[j]+i]){
                        /* found the decomposition */
                        traced = j; jj = i - 1; gq = 1;
                        break;
                      }
                    }

                    if (traced) break;
                  }
                  break;
      }

      if (!traced) vrna_message_error("backtrack failed in f5");
      /* push back the remaining f5 portion */
      bt_stack[++s].i = 1;
      bt_stack[s].j   = jj;
      bt_stack[s].ml  = ml;

      /* trace back the base pair found */
      j=traced;

      if(with_gquad && gq){
        /* goto backtrace of gquadruplex */
        goto repeat_gquad_comparative;
      }

      bp_stack[++b].i = i;
      bp_stack[b].j   = j;
      cov_en += pscore[indx[j]+i];
      goto repeat1_comparative;
    }
    else { /* trace back in fML array */
      if(hc->up_ml[i]){
        en = fML[indx[j]+i+1] + n_seq * P->MLbase;

        if(sc)
          for(ss = 0; ss < n_seq; ss++)
            if(sc[ss]){
              if(sc[ss]->energy_up)
                en += sc[ss]->energy_up[a2s[ss][i]][1];
            }

        if(en == fij) { /* 5' end is unpaired */
          bt_stack[++s].i = i+1;
          bt_stack[s].j   = j;
          bt_stack[s].ml  = ml;
          continue;
        }
      }

      if(md->gquad){
        if(fij == ggg[indx[j]+i] + n_seq * E_MLstem(0, -1, -1, P)){
          /* go to backtracing of quadruplex */
          goto repeat_gquad_comparative;
        }
      }

      if(hc->matrix[indx[j] + i] & VRNA_CONSTRAINT_CONTEXT_MB_LOOP_ENC){
        cij = c[indx[j]+i];
        if(dangle_model){
          for(ss = 0; ss < n_seq; ss++){
            tt = md->pair[S[ss][i]][S[ss][j]];
            if(tt==0) tt=7;
            cij += E_MLstem(tt, S5[ss][i], S3[ss][j], P);
          }
        }
        else{
          for(ss = 0; ss < n_seq; ss++){
            tt = md->pair[S[ss][i]][S[ss][j]];
            if(tt==0) tt=7;
            cij += E_MLstem(tt, -1, -1, P);
          }
        }

        if (fij==cij){
          /* found a pair */
          bp_stack[++b].i = i;
          bp_stack[b].j   = j;
          cov_en += pscore[indx[j]+i];
          goto repeat1_comparative;
        }
      }

      for (k = i+1+turn; k <= j-2-turn; k++)
        if (fij == (fML[indx[k]+i]+fML[indx[j]+k+1]))
          break;

      bt_stack[++s].i = i;
      bt_stack[s].j   = k;
      bt_stack[s].ml  = ml;
      bt_stack[++s].i = k+1;
      bt_stack[s].j   = j;
      bt_stack[s].ml  = ml;

      if (k>j-2-turn) vrna_message_error("backtrack failed in fML");
      continue;
    }

  repeat1_comparative:

    /*----- begin of "repeat:" -----*/
    if (canonical)  cij = c[indx[j]+i];

    for (ss=0; ss<n_seq; ss++) {
      type[ss] = md->pair[S[ss][i]][S[ss][j]];
      if (type[ss]==0) type[ss] = 7;
    }

    if (md->noLP)
      if (cij == c[indx[j]+i]) {
        /* (i.j) closes canonical structures, thus
           (i+1.j-1) must be a pair                */
        for (ss=0; ss<n_seq; ss++) {
          type_2 = md->pair[S[ss][j-1]][S[ss][i+1]];  /* j,i not i,j */
          if (type_2==0) type_2 = 7;
          cij -= P->stack[type[ss]][type_2];
          if(sc){
            if(sc[ss]->energy_bp)
              cij -= sc[s]->energy_bp[indx[j] + i];
          }
        }
        cij += pscore[indx[j]+i];
        bp_stack[++b].i = i+1;
        bp_stack[b].j   = j-1;
        cov_en += pscore[indx[j-1]+i+1];
        i++; j--;
        canonical=0;
        goto repeat1_comparative;
      }
    canonical = 1;
    cij += pscore[indx[j]+i];
    /* does (i,j) close a hairpin loop ? */
    if(vrna_BT_hp_loop(vc, i, j, cij, bp_stack, &b))
      continue;

    for (p = i+1; p <= MIN2(j-2-turn,i+MAXLOOP+1); p++) {
      minq = j-i+p-MAXLOOP-2;
      if (minq<p+1+turn) minq = p+1+turn;
      if(hc->up_int[i+1] < (p - i - 1)) break;

      for (q = j-1; q >= minq; q--) {

        if(hc->up_int[q+1] < (j - q - 1)) break;

        if (c[indx[q]+p]>=INF) continue;

        for (ss=energy=0; ss<n_seq; ss++) {
          int u1 = a2s[ss][p-1] - a2s[ss][i];
          int u2 = a2s[ss][j-1] - a2s[ss][q];
          type_2 = md->pair[S[ss][q]][S[ss][p]];  /* q,p not p,q */
          if (type_2==0) type_2 = 7;
          energy += E_IntLoop(u1, u2, type[ss], type_2, S3[ss][i], S5[ss][j], S5[ss][p], S3[ss][q], P);

        }

        if(sc)
          for(ss = 0; ss < n_seq; ss++)
            if(sc[ss]){
              int u1 = a2s[ss][p-1] - a2s[ss][i];
              int u2 = a2s[ss][j-1] - a2s[ss][q];
/*
              int u1 = p - i - 1;
              int u2 = j - q - 1;
*/
              if(u1 + u2 == 0)
                if(sc[ss]->energy_stack){
                  if(S[ss][i] && S[ss][j] && S[ss][p] && S[ss][q]){ /* don't allow gaps in stack */
                    energy +=   sc[ss]->energy_stack[a2s[ss][i]]
                              + sc[ss]->energy_stack[a2s[ss][p]]
                              + sc[ss]->energy_stack[a2s[ss][q]]
                              + sc[ss]->energy_stack[a2s[ss][j]];
                  }
                }
              if(sc[ss]->energy_bp)
                energy += sc[ss]->energy_bp[indx[j] + i];

              if(sc[ss]->energy_up)
                energy +=   sc[ss]->energy_up[a2s[ss][i] + 1][u1]
                          + sc[ss]->energy_up[a2s[ss][q] + 1][u2];
            }

        traced = (cij == energy+c[indx[q]+p]);
        if (traced) {
          bp_stack[++b].i = p;
          bp_stack[b].j   = q;
          cov_en += pscore[indx[q]+p];
          i = p, j = q;
          goto repeat1_comparative;
        }
      }
    }

    /* end of repeat: --------------------------------------------------*/

    /* (i.j) must close a multi-loop */

    i1 = i+1;
    j1 = j-1;

    if(with_gquad){
      /*
        The case that is handled here actually resembles something like
        an interior loop where the enclosing base pair is of regular
        kind and the enclosed pair is not a canonical one but a g-quadruplex
        that should then be decomposed further...
      */
      mm = 0;
      for(ss=0;ss<n_seq;ss++){
        tt = type[ss];
        if(tt == 0) tt = 7;
        if(dangle_model == 2)
          mm += P->mismatchI[tt][S3[ss][i]][S5[ss][j]];
        if(tt > 2)
          mm += P->TerminalAU;
      }

      for(p = i + 2;
        p < j - VRNA_GQUAD_MIN_BOX_SIZE;
        p++){
        if(S_cons[p] != 3) continue;
        l1    = p - i - 1;
        if(l1>MAXLOOP) break;
        minq  = j - i + p - MAXLOOP - 2;
        c0    = p + VRNA_GQUAD_MIN_BOX_SIZE - 1;
        minq  = MAX2(c0, minq);
        c0    = j - 1;
        maxq  = p + VRNA_GQUAD_MAX_BOX_SIZE + 1;
        maxq  = MIN2(c0, maxq);
        for(q = minq; q < maxq; q++){
          if(S_cons[q] != 3) continue;
          c0  = mm + ggg[indx[q] + p] + n_seq * P->internal_loop[l1 + j - q - 1];
          if(cij == c0){
            i=p;j=q;
            goto repeat_gquad_comparative;
          }
        }
      }
      p = i1;
      if(S_cons[p] == 3){
        if(p < j - VRNA_GQUAD_MIN_BOX_SIZE){
          minq  = j - i + p - MAXLOOP - 2;
          c0    = p + VRNA_GQUAD_MIN_BOX_SIZE - 1;
          minq  = MAX2(c0, minq);
          c0    = j - 3;
          maxq  = p + VRNA_GQUAD_MAX_BOX_SIZE + 1;
          maxq  = MIN2(c0, maxq);
          for(q = minq; q < maxq; q++){
            if(S_cons[q] != 3) continue;
            if(cij == mm + ggg[indx[q] + p] + n_seq * P->internal_loop[j - q - 1]){
              i = p; j=q;
              goto repeat_gquad_comparative;
            }
          }
        }
      }
      q = j1;
      if(S_cons[q] == 3)
        for(p = i1 + 3; p < j - VRNA_GQUAD_MIN_BOX_SIZE; p++){
          l1    = p - i - 1;
          if(l1>MAXLOOP) break;
          if(S_cons[p] != 3) continue;
          if(cij == mm + ggg[indx[q] + p] + n_seq * P->internal_loop[l1]){
            i = p; j = q;
            goto repeat_gquad_comparative;
          }
        }
    }

    if(hc->matrix[indx[j] + i] & VRNA_CONSTRAINT_CONTEXT_MB_LOOP){
      mm = n_seq*P->MLclosing;
      switch(dangle_model){
        case 0:   for(ss = 0; ss < n_seq; ss++){
                    tt = rtype[type[ss]];
                    mm += E_MLstem(tt, -1, -1, P);
                  }
                  break;
        default:  for(ss = 0; ss < n_seq; ss++){
                    tt = rtype[type[ss]];
                    mm += E_MLstem(tt, S5[ss][j], S3[ss][i], P);
                  }
                  break;
      }

      if(sc)
        for(ss = 0; ss < n_seq; ss++)
          if(sc[ss]){
            if(sc[ss]->energy_bp)
              mm += sc[ss]->energy_bp[indx[j] + i];
          }

      bt_stack[s+1].ml  = bt_stack[s+2].ml = 1;

      for (k = i1+turn+1; k < j1-turn-1; k++){
        if(cij == fML[indx[k]+i1] + fML[indx[j1]+k+1] + mm) break;
      }

      if (k<=j-3-turn) { /* found the decomposition */
        bt_stack[++s].i = i1;
        bt_stack[s].j   = k;
        bt_stack[++s].i = k+1;
        bt_stack[s].j   = j1;
      } else {
          vrna_message_error("backtracking failed in repeat");
      }
    } else
      vrna_message_error("backtracking failed in repeat");

    continue; /* this is a workarround to not accidentally proceed in the following block */

  repeat_gquad_comparative:
    /*
      now we do some fancy stuff to backtrace the stacksize and linker lengths
      of the g-quadruplex that should reside within position i,j
    */
    {
      int cnt1, l[3], L, size;
      size = j-i+1;

      for(L=0; L < VRNA_GQUAD_MIN_STACK_SIZE;L++){
        if(S_cons[i+L] != 3) break;
        if(S_cons[j-L] != 3) break;
      }

      if(L == VRNA_GQUAD_MIN_STACK_SIZE){
        /* continue only if minimum stack size starting from i is possible */
        for(; L<=VRNA_GQUAD_MAX_STACK_SIZE;L++){
          if(S_cons[i+L-1] != 3) break; /* break if no more consecutive G's 5' */
          if(S_cons[j-L+1] != 3) break; /* break if no more consecutive G'1 3' */
          for(    l[0] = VRNA_GQUAD_MIN_LINKER_LENGTH;
                  (l[0] <= VRNA_GQUAD_MAX_LINKER_LENGTH)
              &&  (size - 4*L - 2*VRNA_GQUAD_MIN_LINKER_LENGTH - l[0] >= 0);
              l[0]++){
            /* check whether we find the second stretch of consecutive G's */
            for(cnt1 = 0; (cnt1 < L) && (S_cons[i+L+l[0]+cnt1] == 3); cnt1++);
            if(cnt1 < L) continue;
            for(    l[1] = VRNA_GQUAD_MIN_LINKER_LENGTH;
                    (l[1] <= VRNA_GQUAD_MAX_LINKER_LENGTH)
                &&  (size - 4*L - VRNA_GQUAD_MIN_LINKER_LENGTH - l[0] - l[1] >= 0);
                l[1]++){
              /* check whether we find the third stretch of consectutive G's */
              for(cnt1 = 0; (cnt1 < L) && (S_cons[i+2*L+l[0]+l[1]+cnt1] == 3); cnt1++);
              if(cnt1 < L) continue;

              /*
                the length of the third linker now depends on position j as well
                as the other linker lengths... so we do not have to loop too much
              */
              l[2] = size - 4*L - l[0] - l[1];
              if(l[2] < VRNA_GQUAD_MIN_LINKER_LENGTH) break;
              if(l[2] > VRNA_GQUAD_MAX_LINKER_LENGTH) continue;
              /* check for contribution */
              if(ggg[indx[j]+i] == E_gquad_ali(i, L, l, (const short **)S, n_seq, P)){
                int a;
                /* fill the G's of the quadruplex into base pair stack */
                for(a=0;a<L;a++){
                  bp_stack[++b].i = i+a;
                  bp_stack[b].j   = i+a;
                  bp_stack[++b].i = i+L+l[0]+a;
                  bp_stack[b].j   = i+L+l[0]+a;
                  bp_stack[++b].i = i+L+l[0]+L+l[1]+a;
                  bp_stack[b].j   = i+L+l[0]+L+l[1]+a;
                  bp_stack[++b].i = i+L+l[0]+L+l[1]+L+l[2]+a;
                  bp_stack[b].j   = i+L+l[0]+L+l[1]+L+l[2]+a;
                }
                goto repeat_gquad_comparative_exit;
              }
            }
          }
        }
      }
      vrna_message_error("backtracking failed in repeat_gquad_comparative");
    }
  repeat_gquad_comparative_exit:
    __asm("nop");

  }

  bp_stack[0].i = b;    /* save the total number of base pairs */
  free(type);
}

