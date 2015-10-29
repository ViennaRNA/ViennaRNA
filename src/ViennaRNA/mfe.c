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

#include <config.h>
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
PRIVATE vrna_plist_t  *backtrack(vrna_fold_compound_t *vc, vrna_bp_stack_t *bp_stack, sect bt_stack[], int s);

/*
#################################
# BEGIN OF FUNCTION DEFINITIONS #
#################################
*/

PUBLIC float
vrna_mfe(vrna_fold_compound_t *vc,
          char *structure){

  int     length, energy, s;
  char    *ss;
  sect    bt_stack[MAXSECTORS]; /* stack of partial structures for backtracking */
  vrna_bp_stack_t   *bp;


  s       = 0;
  length  = (int) vc->length;

  /* initialize generalized hard constraints if necessary */
  if(vc->hc->pre)
    vc->hc->pre(vc, VRNA_SC_GEN_MFE);

  /* initialize generalized soft constraints if necessary */
  if(vc->sc)
    if(vc->sc->pre)
      vc->sc->pre(vc, VRNA_SC_GEN_MFE);

  energy = fill_arrays(vc);

  if(vc->params->model_details.circ){
    fill_arrays_circ(vc, bt_stack, &s);
    energy = vc->matrices->Fc;
  }

  if(structure && vc->params->model_details.backtrack){
    bp = (vrna_bp_stack_t *)vrna_alloc(sizeof(vrna_bp_stack_t) * (4*(1+length/2))); /* add a guess of how many G's may be involved in a G quadruplex */

    backtrack(vc, bp, bt_stack, s);

    ss = vrna_db_from_bp_stack(bp, length);
    strncpy(structure, ss, length + 1);
    free(ss);

#ifdef  VRNA_BACKWARD_COMPAT
    /*
    *  Backward compatibility:
    *  This block may be removed if deprecated functions
    *  relying on the global variable "base_pair" vanish from within the package!
    */
    {
      if(base_pair) free(base_pair);
      base_pair = bp;
    }
#endif

  }

  if(vc->hc->post)
    vc->hc->post(vc, VRNA_SC_GEN_MFE);

  if(vc->sc)
    if(vc->sc->post)
      vc->sc->post(vc, VRNA_SC_GEN_MFE);

  if (vc->params->model_details.backtrack_type=='C')
    return (float) vc->matrices->c[vc->jindx[length]+1]/100.;
  else if (vc->params->model_details.backtrack_type=='M')
    return (float) vc->matrices->fML[vc->jindx[length]+1]/100.;
  else
    return (float) energy/100.;
}

/**
*** fill "c", "fML" and "f5" arrays and return  optimal energy
**/
PRIVATE int
fill_arrays(vrna_fold_compound_t *vc){

  int               i, j, ij, length, energy, new_c, stackEnergy, no_close, type_2;
  int               noGUclosure, noLP, uniq_ML, with_gquad, dangle_model, *rtype, *indx;
  int               *my_f5, *my_c, *my_fML, *my_fM1, *my_ggg, hc_decompose, *hc_up_ml;
  int               *cc, *cc1;  /* auxilary arrays for canonical structures     */
  int               *Fmi;       /* holds row i of fML (avoids jumps in memory)  */
  int               *DMLi;      /* DMLi[j] holds  MIN(fML[i,k]+fML[k+1,j])      */
  int               *DMLi1;     /*                MIN(fML[i+1,k]+fML[k+1,j])    */
  int               *DMLi2;     /*                MIN(fML[i+2,k]+fML[k+1,j])    */
  unsigned char     type;
  char              *ptype, *hard_constraints;
  short             *S1;
  vrna_param_t      *P;
  vrna_mx_mfe_t     *matrices;
  vrna_hc_t         *hc;
  vrna_sc_t         *sc;

  length            = (int)vc->length;
  ptype             = vc->ptype;
  indx              = vc->jindx;
  P                 = vc->params;
  S1                = vc->sequence_encoding;
  noGUclosure       = P->model_details.noGUclosure;
  noLP              = P->model_details.noLP;
  uniq_ML           = P->model_details.uniq_ML;
  with_gquad        = P->model_details.gquad;
  dangle_model      = P->model_details.dangles;
  rtype             = &(P->model_details.rtype[0]);
  hc                = vc->hc;
  hard_constraints  = hc->matrix;
  hc_up_ml          = hc->up_ml;
  sc                = vc->sc;
  matrices          = vc->matrices;
  my_f5             = matrices->f5;
  my_c              = matrices->c;
  my_fML            = matrices->fML;
  my_fM1            = matrices->fM1;
  my_ggg            = matrices->ggg;


  /* allocate memory for all helper arrays */
  cc    = (int *) vrna_alloc(sizeof(int)*(length + 2));
  cc1   = (int *) vrna_alloc(sizeof(int)*(length + 2));
  Fmi   = (int *) vrna_alloc(sizeof(int)*(length + 1));
  DMLi  = (int *) vrna_alloc(sizeof(int)*(length + 1));
  DMLi1 = (int *) vrna_alloc(sizeof(int)*(length + 1));
  DMLi2 = (int *) vrna_alloc(sizeof(int)*(length + 1));


  /* prefill helper arrays */
  for(j = 1; j <= length; j++){
    Fmi[j] = DMLi[j] = DMLi1[j] = DMLi2[j] = INF;
  }


  /* prefill matrices with init contributions */
  for(j = 1; j <= length; j++)
    for(i = (j > TURN ? (j - TURN) : 1); i < j; i++){
      my_c[indx[j] + i] = my_fML[indx[j] + i] = INF;
      if(uniq_ML)
        my_fM1[indx[j] + i] = INF;
    }

  /* start recursion */

  if (length <= TURN){
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

  for (i = length-TURN-1; i >= 1; i--) { /* i,j in [1..length] */

    for (j = i+TURN+1; j <= length; j++) {
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
          energy  = E_mb_loop_fast(i, j, vc, DMLi1, DMLi2);
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

      my_fML[ij] = E_ml_stems_fast(i, j, vc, Fmi, DMLi);

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


PUBLIC vrna_plist_t *
vrna_backtrack_from_intervals(vrna_fold_compound_t *vc,
                              vrna_bp_stack_t *bp_stack,
                              sect bt_stack[],
                              int s){

  return backtrack(vc, bp_stack, bt_stack, s);
} 

/**
*** trace back through the "c", "f5" and "fML" arrays to get the
*** base pairing list. No search for equivalent structures is done.
*** This is fast, since only few structure elements are recalculated.
***
*** normally s=0.
*** If s>0 then s items have been already pushed onto the bt_stack
**/
PRIVATE vrna_plist_t *
backtrack(vrna_fold_compound_t *vc,
          vrna_bp_stack_t *bp_stack,
          sect bt_stack[],
          int s){

  int   i, j, ij, k, mm3, length, energy, en, new;
  int   no_close, minq;
  int   b=0;
  unsigned char type, tt, type_2;
  char  *string         = vc->sequence;
  vrna_param_t  *P      = vc->params;
  int     *indx         = vc->jindx;
  char    *ptype        = vc->ptype;

  short *S1             = vc->sequence_encoding;
  short *S              = vc->sequence_encoding2;
  int   dangle_model    = P->model_details.dangles;
  int   noLP            = P->model_details.noLP;
  int   noGUclosure     = P->model_details.noGUclosure;
  int   *rtype          = &(P->model_details.rtype[0]);
  char  backtrack_type  = P->model_details.backtrack_type;
  int   with_gquad      = P->model_details.gquad;

  /* the folding matrices */
  int   *my_f5, *my_c, *my_fML, *my_ggg;

  length  = vc->length;
  my_f5   = vc->matrices->f5;
  my_c    = vc->matrices->c;
  my_fML  = vc->matrices->fML;
  my_ggg  = vc->matrices->ggg;

  vrna_hc_t *hc               = vc->hc;
  vrna_sc_t *sc               = vc->sc;
  char      *hard_constraints = hc->matrix;

  if (s==0) {
    bt_stack[++s].i = 1;
    bt_stack[s].j = length;
    bt_stack[s].ml = (backtrack_type=='M') ? 1 : ((backtrack_type=='C')? 2: 0);
  }
  while (s>0) {
    int ml, fij, fi, cij, traced, i1, j1, p, q, jj=0, gq=0;
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
                    fprintf(stderr, "%s\n", string);
                    vrna_message_error("backtrack failed in f5");
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
                    fprintf(stderr, "%s\n", string);
                    vrna_message_error("backtrack failed in fML");
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
      vrna_message_error("backtracking failed in repeat");
    }

    /* end of repeat: --------------------------------------------------*/

  } /* end of infinite while loop */

  bp_stack[0].i = b;    /* save the total number of base pairs */
}

