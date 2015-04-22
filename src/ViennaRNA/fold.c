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
#include "ViennaRNA/fold.h"

#ifdef _OPENMP
#include <omp.h>
#endif

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

/* some backward compatibility stuff */
PRIVATE int                 backward_compat           = 0;
PRIVATE vrna_fold_compound  *backward_compat_compound = NULL;

#ifdef _OPENMP

#pragma omp threadprivate(backward_compat_compound, backward_compat)

#endif

/*
#################################
# PRIVATE FUNCTION DECLARATIONS #
#################################
*/

PRIVATE int   fill_arrays(vrna_fold_compound *vc);
PRIVATE void  fill_arrays_circ(vrna_fold_compound *vc, sect bt_stack[], int *bt);
PRIVATE plist *backtrack(vrna_fold_compound *vc, bondT *bp_stack, sect bt_stack[], int s);

/* wrappers for old API compatibility */
PRIVATE float wrap_fold( const char *string, char *structure, vrna_param_t *parameters, int is_constrained, int is_circular);
PRIVATE void  wrap_array_export(int **f5_p, int **c_p, int **fML_p, int **fM1_p, int **indx_p, char **ptype_p);
PRIVATE void  wrap_array_export_circ( int *Fc_p, int *FcH_p, int *FcI_p, int *FcM_p, int **fM2_p);

/*
#################################
# BEGIN OF FUNCTION DEFINITIONS #
#################################
*/

PRIVATE void
wrap_array_export(int **f5_p,
                  int **c_p,
                  int **fML_p,
                  int **fM1_p,
                  int **indx_p,
                  char **ptype_p){

  /* make the DP arrays available to routines such as subopt() */
  if(backward_compat_compound){
    *f5_p     = backward_compat_compound->matrices->f5;
    *c_p      = backward_compat_compound->matrices->c;
    *fML_p    = backward_compat_compound->matrices->fML;
    *fM1_p    = backward_compat_compound->matrices->fM1;
    *indx_p   = backward_compat_compound->jindx;
    *ptype_p  = backward_compat_compound->ptype;
  }
}

PRIVATE void
wrap_array_export_circ( int *Fc_p,
                        int *FcH_p,
                        int *FcI_p,
                        int *FcM_p,
                        int **fM2_p){

  /* make the DP arrays available to routines such as subopt() */
  if(backward_compat_compound){
    *Fc_p   = backward_compat_compound->matrices->Fc;
    *FcH_p  = backward_compat_compound->matrices->FcH;
    *FcI_p  = backward_compat_compound->matrices->FcI;
    *FcM_p  = backward_compat_compound->matrices->FcM;
    *fM2_p  = backward_compat_compound->matrices->fM2;
  }
}

PRIVATE float
wrap_fold( const char *string,
          char *structure,
          vrna_param_t *parameters,
          int is_constrained,
          int is_circular){

  vrna_fold_compound  *vc;
  vrna_param_t        *P;

#ifdef _OPENMP
/* Explicitly turn off dynamic threads */
  omp_set_dynamic(0);
#endif

  /* we need the parameter structure for hard constraints */
  if(parameters){
    P = get_parameter_copy(parameters);
  } else {
    vrna_md_t md;
    set_model_details(&md);
    P = get_scaled_parameters(temperature, md);
  }
  P->model_details.circ = is_circular;

  vc = vrna_get_fold_compound(string, &(P->model_details), VRNA_OPTION_MFE);

  if(parameters){ /* replace params if necessary */
    free(vc->params);
    vc->params = P;
  } else {
    free(P);
  }

  /* handle hard constraints in pseudo dot-bracket format if passed via simple interface */
  if(is_constrained && structure){
    unsigned int constraint_options = 0;
    constraint_options |= VRNA_CONSTRAINT_DB
                          | VRNA_CONSTRAINT_DB_PIPE
                          | VRNA_CONSTRAINT_DB_DOT
                          | VRNA_CONSTRAINT_DB_X
                          | VRNA_CONSTRAINT_DB_ANG_BRACK
                          | VRNA_CONSTRAINT_DB_RND_BRACK;

    vrna_add_constraints(vc, (const char *)structure, constraint_options);
  }

  if(backward_compat_compound && backward_compat)
    vrna_free_fold_compound(backward_compat_compound);

  backward_compat_compound  = vc;
  backward_compat           = 1;

  return vrna_fold(vc, structure);
}

PUBLIC float
vrna_fold(vrna_fold_compound *vc,
          char *structure){

  int     length, energy, s;
  sect    bt_stack[MAXSECTORS]; /* stack of partial structures for backtracking */
  bondT   *bp;


  s       = 0;
  length  = (int) vc->length;

  if(vc->sc)
    if(vc->sc->pre)
      vc->sc->pre(vc, VRNA_SC_GEN_MFE);

  energy = fill_arrays(vc);

  if(vc->params->model_details.circ){
    fill_arrays_circ(vc, bt_stack, &s);
    energy = vc->matrices->Fc;
  }

  if(structure && vc->params->model_details.backtrack){
    bp = (bondT *)vrna_alloc(sizeof(bondT) * (4*(1+length/2))); /* add a guess of how many G's may be involved in a G quadruplex */

    backtrack(vc, bp, bt_stack, s);

    vrna_parenthesis_structure(structure, bp, length);

    /*
    *  Backward compatibility:
    *  This block may be removed if deprecated functions
    *  relying on the global variable "base_pair" vanish from within the package!
    */
    {
      if(base_pair) free(base_pair);
      base_pair = bp;
    }
  }

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
fill_arrays(vrna_fold_compound *vc){

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

      if (hc_decompose) {   /* we have a pair */
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
        energy = E_int_loop(i, j, vc);
        new_c = MIN2(new_c, energy);

        /* remember stack energy for --noLP option */
        if(noLP){
          stackEnergy = E_stack(i, j, vc);
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


PUBLIC plist *
vrna_backtrack_from_intervals(vrna_fold_compound *vc,
                              bondT *bp_stack,
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
PRIVATE plist *
backtrack(vrna_fold_compound *vc,
          bondT *bp_stack,
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
    i  = bt_stack[s].i;
    j  = bt_stack[s].j;
    ml = bt_stack[s--].ml;   /* ml is a flag indicating if backtracking is to
                              occur in the fML- (1) or in the f-array (0) */
    if (ml==2) {
      bp_stack[++b].i = i;
      bp_stack[b].j   = j;
      goto repeat1;
    }

    if (j < i+TURN+1) continue; /* no more pairs in this interval */

    if(ml == 1){
      fij = my_fML[indx[j] + i];
      fi  = (hc->up_ml[j]) ? my_fML[indx[j - 1] + i] + P->MLbase : INF;
    } else {
      fij = my_f5[j];
      fi  = (hc->up_ext[j]) ? my_f5[j-1] : INF;
    }

    if(sc){
      if(sc->free_energies)
        fi += sc->free_energies[j][1];
    }

    if (fij == fi) {  /* 3' end is unpaired */
      bt_stack[++s].i = i;
      bt_stack[s].j   = j-1;
      bt_stack[s].ml  = ml;
      continue;
    }

    if (ml == 0) { /* backtrack in f5 */
      switch(dangle_model){
        case 0:   /* j is paired. Find pairing partner */
                  for(k=j-TURN-1,traced=0; k>=1; k--){

                    if(with_gquad){
                      if(fij == my_f5[k-1] + my_ggg[indx[j]+k]){
                        /* found the decomposition */
                        traced = j; jj = k - 1; gq = 1;
                        break;
                      }
                    }

                    if (hc->matrix[indx[j] + k] & VRNA_CONSTRAINT_CONTEXT_EXT_LOOP){
                      type = (unsigned char)ptype[indx[j]+k];
                      if(fij == E_ExtLoop(type, -1, -1, P) + my_c[indx[j]+k] + my_f5[k-1]){
                        traced=j; jj = k-1;
                        break;
                      }
                    }
                  }
                  break;

        case 2:   mm3 = (j<length) ? S1[j+1] : -1;
                  for(k=j-TURN-1,traced=0; k>=1; k--){

                    if(with_gquad){
                      if(fij == my_f5[k-1] + my_ggg[indx[j]+k]){
                        /* found the decomposition */
                        traced = j; jj = k - 1; gq = 1;
                        break;
                      }
                    }

                    if (hc->matrix[indx[j] + k] & VRNA_CONSTRAINT_CONTEXT_EXT_LOOP){
                      type = (unsigned char)ptype[indx[j]+k];
                      if(fij == E_ExtLoop(type, (k>1) ? S1[k-1] : -1, mm3, P) + my_c[indx[j]+k] + my_f5[k-1]){
                        traced=j; jj = k-1;
                        break;
                      }
                    }
                  }
                  break;

        default:  for(traced = 0, k=j-TURN-1; k>1; k--){

                    if(with_gquad){
                      if(fij == my_f5[k-1] + my_ggg[indx[j]+k]){
                        /* found the decomposition */
                        traced = j; jj = k - 1; gq = 1;
                        break;
                      }
                    }

                    if (hc->matrix[indx[j] + k] & VRNA_CONSTRAINT_CONTEXT_EXT_LOOP){
                      type = (unsigned char)ptype[indx[j] + k];
                      en = my_c[indx[j] + k];
                      if(fij == my_f5[k-1] + en + E_ExtLoop(type, -1, -1, P)){
                        traced = j;
                        jj = k-1;
                        break;
                      }
                      if(hc->up_ext[k-1]){
                        int tmp_en = fij;
                        if(sc)
                          if(sc->free_energies)
                            tmp_en -= sc->free_energies[k-1][1];

                        if(tmp_en == my_f5[k-2] + en + E_ExtLoop(type, S1[k-1], -1, P)){
                          traced = j;
                          jj = k-2;
                          break;
                        }
                      }
                    }
                    if (hc->matrix[indx[j-1] + k] & VRNA_CONSTRAINT_CONTEXT_EXT_LOOP){
                      type = (unsigned char)ptype[indx[j-1] + k];
                      en = my_c[indx[j-1] + k];
                      if(hc->up_ext[j]){
                        int tmp_en = fij;
                        if(sc)
                          if(sc->free_energies)
                            tmp_en -= sc->free_energies[j][1];

                        if(tmp_en == my_f5[k-1] + en + E_ExtLoop(type, -1, S1[j], P)){
                          traced = j-1;
                          jj = k-1;
                          break;
                        }

                        if(hc->up_ext[k-1]){
                          if(sc)
                            if(sc->free_energies)
                              tmp_en -= sc->free_energies[k-1][1];

                          if(tmp_en == my_f5[k-2] + en + E_ExtLoop(type, S1[k-1], S1[j], P)){
                            traced = j-1;
                            jj = k-2;
                            break;
                          }
                        }
                      }
                    }
                  }
                  if(!traced){

                    if(with_gquad){
                      if(fij == my_ggg[indx[j]+1]){
                        /* found the decomposition */
                        traced = j; jj = 0; gq = 1;
                        break;
                      }
                    }

                    if (hc->matrix[indx[j] + 1] & VRNA_CONSTRAINT_CONTEXT_EXT_LOOP){
                      type = (unsigned char)ptype[indx[j]+1];
                      if(fij == my_c[indx[j]+1] + E_ExtLoop(type, -1, -1, P)){
                        traced = j;
                        jj = 0;
                        break;
                      }
                    }
                    if (hc->matrix[indx[j-1] + 1] & VRNA_CONSTRAINT_CONTEXT_EXT_LOOP){
                      if(hc->up_ext[j]){
                        int tmp_en = fij;
                        if(sc)
                          if(sc->free_energies)
                            tmp_en -= sc->free_energies[j][1];

                        type = (unsigned char)ptype[indx[j-1]+1];
                        if(tmp_en == my_c[indx[j-1]+1] + E_ExtLoop(type, -1, S1[j], P)){
                          traced = j-1;
                          jj = 0;
                          break;
                        }
                      }
                    }
                  }
                  break;
      }

      if (!traced){
        fprintf(stderr, "%s\n", string);
        vrna_message_error("backtrack failed in f5");
      }
      bt_stack[++s].i = 1;
      bt_stack[s].j   = jj;
      bt_stack[s].ml  = ml;

      /* trace back the base pair found */
      i=k; j=traced;

      if(with_gquad && gq){
        /* goto backtrace of gquadruplex */
        goto repeat_gquad;
      }

      bp_stack[++b].i = i;
      bp_stack[b].j   = j;
      goto repeat1;
    }
    else { /* trace back in fML array */
      if(hc->up_ml[i]){
        en = my_fML[indx[j]+i+1]+P->MLbase;

        if(sc)
          if(sc->free_energies)
            en += sc->free_energies[i][1];

        if (en == fij) { /* 5' end is unpaired */
          bt_stack[++s].i = i+1;
          bt_stack[s].j   = j;
          bt_stack[s].ml  = ml;
          continue;
        }
      }

      ij  = indx[j]+i;

      if(with_gquad){
        if(fij == my_ggg[ij] + E_MLstem(0, -1, -1, P)){
          /* go to backtracing of quadruplex */
          goto repeat_gquad;
        }
      }

      tt  = (unsigned char)ptype[ij];
      en  = my_c[ij];
      switch(dangle_model){
        case 0:   if(hc->matrix[ij] & VRNA_CONSTRAINT_CONTEXT_MB_LOOP_ENC){
                    if(fij == en + E_MLstem(tt, -1, -1, P)){
                      bp_stack[++b].i = i;
                      bp_stack[b].j   = j;
                      goto repeat1;
                    }
                  }
                  break;

        case 2:   if(hc->matrix[ij] & VRNA_CONSTRAINT_CONTEXT_MB_LOOP_ENC){
                    if(fij == en + E_MLstem(tt, S1[i-1], S1[j+1], P)){
                      bp_stack[++b].i = i;
                      bp_stack[b].j   = j;
                      goto repeat1;
                    }
                  }
                  break;

        default:  if(hc->matrix[ij] & VRNA_CONSTRAINT_CONTEXT_MB_LOOP_ENC){
                    if(fij == en + E_MLstem(tt, -1, -1, P)){
                      bp_stack[++b].i = i;
                      bp_stack[b].j   = j;
                      goto repeat1;
                    }
                  }
                  if(hc->matrix[ij+1] & VRNA_CONSTRAINT_CONTEXT_MB_LOOP_ENC){
                    if(hc->up_ml[i]){
                      int tmp_en = fij;
                      if(sc)
                        if(sc->free_energies)
                          tmp_en -= sc->free_energies[i][1];

                      tt = (unsigned char)ptype[ij+1];
                      if(tmp_en == my_c[ij+1] + E_MLstem(tt, S1[i], -1, P) + P->MLbase){
                        bp_stack[++b].i = ++i;
                        bp_stack[b].j   = j;
                        goto repeat1;
                      }
                    }
                  }
                  if(hc->matrix[indx[j-1]+i] & VRNA_CONSTRAINT_CONTEXT_MB_LOOP_ENC){
                    if(hc->up_ml[j]){
                      int tmp_en = fij;
                      if(sc)
                        if(sc->free_energies)
                          tmp_en -= sc->free_energies[j][1];

                      tt = (unsigned char)ptype[indx[j-1]+i];
                      if(tmp_en == my_c[indx[j-1]+i] + E_MLstem(tt, -1, S1[j], P) + P->MLbase){
                        bp_stack[++b].i = i;
                        bp_stack[b].j   = --j;
                        goto repeat1;
                      }
                    }
                  }
                  if(hc->matrix[indx[j-1]+i+1] & VRNA_CONSTRAINT_CONTEXT_MB_LOOP_ENC){
                    if(hc->up_ml[i] && hc->up_ml[j]){
                      int tmp_en = fij;
                      if(sc)
                        if(sc->free_energies)
                          tmp_en -= sc->free_energies[i][1] + sc->free_energies[j][1];

                      tt = (unsigned char)ptype[indx[j-1]+i+1];
                      if(tmp_en == my_c[indx[j-1]+i+1] + E_MLstem(tt, S1[i], S1[j], P) + 2*P->MLbase){
                        bp_stack[++b].i = ++i;
                        bp_stack[b].j   = --j;
                        goto repeat1;
                      }
                    }
                  }
                  break;
      }

      for(k = i + 1 + TURN; k <= j - 2 - TURN; k++)
        if(fij == (my_fML[indx[k]+i]+my_fML[indx[j]+k+1]))
          break;

      if ((dangle_model==3)&&(k > j - 2 - TURN)) { /* must be coax stack */
        int ik, k1j, tmp_en;
        ml = 2;
        for (k1j = indx[j]+i+TURN+2, k = i+1+TURN; k <= j - 2 - TURN; k++, k1j++) {
          ik = indx[k]+i;
          if((hc->matrix[ik] & VRNA_CONSTRAINT_CONTEXT_MB_LOOP_ENC) && (hc->matrix[k1j] & VRNA_CONSTRAINT_CONTEXT_MB_LOOP_ENC)){
            type    = rtype[(unsigned char)ptype[ik]];
            type_2  = rtype[(unsigned char)ptype[k1j]];
            tmp_en  = my_c[ik] + my_c[k1j] + P->stack[type][type_2] + 2*P->MLintern[1];
            if (fij == tmp_en)
              break;
          }
        }
      }
      bt_stack[++s].i = i;
      bt_stack[s].j   = k;
      bt_stack[s].ml  = ml;
      bt_stack[++s].i = k+1;
      bt_stack[s].j   = j;
      bt_stack[s].ml  = ml;

      if (k>j-2-TURN) vrna_message_error("backtrack failed in fML");
      continue;
    }

  repeat1:

    /*----- begin of "repeat:" -----*/
    ij = indx[j]+i;
    if (canonical)  cij = my_c[ij];

    type = (unsigned char)ptype[ij];

    if (noLP)
      if (cij == my_c[ij]){
        /* (i.j) closes canonical structures, thus
           (i+1.j-1) must be a pair                */
        if((hard_constraints[ij] & VRNA_CONSTRAINT_CONTEXT_INT_LOOP) && (hard_constraints[indx[j-1]+i+1] & VRNA_CONSTRAINT_CONTEXT_INT_LOOP_ENC)){
          type_2 = (unsigned char)ptype[indx[j-1]+i+1];
          type_2 = rtype[type_2];
          cij -= P->stack[type][type_2];
          if(sc){
            if(sc->en_basepair)
              cij -= sc->en_basepair[ij];
              if(sc->en_stack)
                cij -=    sc->en_stack[i]
                        + sc->en_stack[i+1]
                        + sc->en_stack[j-1]
                        + sc->en_stack[j];
          }
          bp_stack[++b].i = i+1;
          bp_stack[b].j   = j-1;
          i++; j--;
          canonical=0;
          goto repeat1;
        }
      }
    canonical = 1;


    no_close = (((type==3)||(type==4))&&noGUclosure);
    if (no_close) {
      if (cij == FORBIDDEN) continue;
    } else {
      en = vrna_E_hp_loop(vc, i, j);

      if (cij == en){
        if(sc)
          if(sc->bt){
            PAIR *ptr, *aux_bps;
            aux_bps = sc->bt(i, j, p, q, VRNA_DECOMP_PAIR_HP, sc->data);
            for(ptr = aux_bps; ptr && ptr->i != -1; ptr++){
              bp_stack[++b].i = ptr->i;
              bp_stack[b].j   = ptr->j;
            }
            free(aux_bps);
          }
        continue;
      }
    }

    if(hc->matrix[ij] & VRNA_CONSTRAINT_CONTEXT_INT_LOOP)
      for (p = i+1; p <= MIN2(j-2-TURN,i+MAXLOOP+1); p++) {
        minq = j-i+p-MAXLOOP-2;
        if (minq<p+1+TURN) minq = p+1+TURN;
        if(hc->up_int[i+1] < (p - i - 1)) break;

        for (q = j-1; q >= minq; q--) {
          if(hc->up_int[q+1] < (j - q - 1)) break;

          type_2 = (unsigned char)ptype[indx[q]+p];

          if (!(hc->matrix[indx[q]+p] & VRNA_CONSTRAINT_CONTEXT_INT_LOOP_ENC)) continue;

          type_2 = rtype[type_2];
          if (noGUclosure)
            if (no_close||(type_2==3)||(type_2==4))
              if ((p>i+1)||(q<j-1)) continue;  /* continue unless stack */

          energy = ubf_eval_int_loop( i, j, p, q,
                                      p-i-1, j-q-1,
                                      S1[i+1], S1[j-1], S1[p-1], S1[q+1],
                                      type, type_2,
                                      rtype,
                                      ij,
                                      -1,
                                      P,
                                      sc);
          new = energy+my_c[indx[q]+p];

          traced = (cij == new);
          if (traced) {
            bp_stack[++b].i = p;
            bp_stack[b].j   = q;
            if(sc)
              if(sc->bt){
                PAIR *ptr, *aux_bps;
                aux_bps = sc->bt(i, j, p, q, VRNA_DECOMP_PAIR_IL, sc->data);
                for(ptr = aux_bps; ptr && ptr->i != -1; ptr++){
                  bp_stack[++b].i = ptr->i;
                  bp_stack[b].j   = ptr->j;
                }
                free(aux_bps);
              }
            i = p, j = q;
            goto repeat1;
          }
        }
      }

    /* end of repeat: --------------------------------------------------*/

    /* (i.j) must close a multi-loop */
    tt = rtype[type];
    i1 = i+1; j1 = j-1;

    if(with_gquad){
      /*
        The case that is handled here actually resembles something like
        an interior loop where the enclosing base pair is of regular
        kind and the enclosed pair is not a canonical one but a g-quadruplex
        that should then be decomposed further...
      */
      if(backtrack_GQuad_IntLoop(cij, i, j, type, S, my_ggg, indx, &p, &q, P)){
        i = p; j = q;
        goto repeat_gquad;
      }
    }
    k = j - 2 - TURN; /* end of loop */

    if(hc->matrix[ij] & VRNA_CONSTRAINT_CONTEXT_MB_LOOP){
      bt_stack[s+1].ml  = bt_stack[s+2].ml = 1;

      switch(dangle_model){
        case 0:   en = cij - E_MLstem(tt, -1, -1, P) - P->MLclosing;
                  if(sc){
                    if(sc->en_basepair)
                      en -= sc->en_basepair[ij];
                  }
                  for(k = i+2+TURN; k < j-2-TURN; k++){
                    if(en == my_fML[indx[k]+i+1] + my_fML[indx[j-1]+k+1])
                      break;
                  }
                  break;

        case 2:   en = cij - E_MLstem(tt, S1[j-1], S1[i+1], P) - P->MLclosing;
                  if(sc){
                    if(sc->en_basepair)
                      en -= sc->en_basepair[ij];
                  }
                  for(k = i+2+TURN; k < j-2-TURN; k++){
                      if(en == my_fML[indx[k]+i+1] + my_fML[indx[j-1]+k+1])
                        break;
                  }
                  break;

        default:  en = cij - P->MLclosing;
                  if(sc){
                    if(sc->en_basepair)
                      en -= sc->en_basepair[ij];
                  }
                  for(k = i+2+TURN; k < j-2-TURN; k++){
                    if(en == my_fML[indx[k]+i+1] + my_fML[indx[j-1]+k+1] + E_MLstem(tt, -1, -1, P)){
                      break;
                    }
                    if(hc->up_ml[i+1]){
                      int tmp_en = en;
                      if(sc)
                        if(sc->free_energies)
                          tmp_en -= sc->free_energies[i+1][1];

                      if(tmp_en == my_fML[indx[k]+i+2] + my_fML[indx[j-1]+k+1] + E_MLstem(tt, -1, S1[i+1], P) + P->MLbase){
                        i1 = i+2;
                        break;
                      }
                    }
                    if(hc->up_ml[j-1]){
                      int tmp_en = en;
                      if(sc)
                        if(sc->free_energies)
                          tmp_en -= sc->free_energies[j-1][1];

                      if(tmp_en == my_fML[indx[k]+i+1] + my_fML[indx[j-2]+k+1] + E_MLstem(tt, S1[j-1], -1, P) + P->MLbase){
                        j1 = j-2;
                        break;
                      }
                    }
                    if((hc->up_ml[i+1]) && (hc->up_ml[j-1])){
                      int tmp_en = en;
                      if(sc)
                        if(sc->free_energies)
                          tmp_en -= sc->free_energies[i+1][1] + sc->free_energies[j-1][1];

                      if(tmp_en == my_fML[indx[k]+i+2] + my_fML[indx[j-2]+k+1] + E_MLstem(tt, S1[j-1], S1[i+1], P) + 2*P->MLbase){
                        i1 = i+2;
                        j1 = j-2;
                        break;
                      }
                    }
                    /* coaxial stacking of (i.j) with (i+1.k) or (k.j-1) */
                    /* use MLintern[1] since coax stacked pairs don't get TerminalAU */
                    if(dangle_model == 3){
                      int tmp_en = en;
                      if(hc->matrix[indx[k]+i+1] & VRNA_CONSTRAINT_CONTEXT_MB_LOOP_ENC){
                        type_2 = rtype[(unsigned char)ptype[indx[k]+i+1]];
                        tmp_en = my_c[indx[k]+i+1]+P->stack[type][type_2]+my_fML[indx[j-1]+k+1];
                        if(sc){
                          if(sc->en_basepair)
                            tmp_en += sc->en_basepair[ij];
                        }
                        if (cij == tmp_en+2*P->MLintern[1]+P->MLclosing) {
                          ml = 2;
                          bt_stack[s+1].ml  = 2;
                          traced = 1;
                          break;
                        }
                      }
                      if(hc->matrix[indx[j-1]+k+1] & VRNA_CONSTRAINT_CONTEXT_MB_LOOP_ENC){
                        type_2 = rtype[(unsigned char)ptype[indx[j-1]+k+1]];
                        tmp_en = my_c[indx[j-1]+k+1]+P->stack[type][type_2]+my_fML[indx[k]+i+1];
                        if(sc){
                          if(sc->en_basepair)
                            tmp_en += sc->en_basepair[ij];
                        }
                        if (cij == tmp_en+2*P->MLintern[1]+P->MLclosing) {
                          bt_stack[s+2].ml = 2;
                          traced = 1;
                          break;
                        }
                      }
                    }
                  }
                  break;
      }
    }

    if (k<=j-3-TURN) { /* found the decomposition */
      bt_stack[++s].i = i1;
      bt_stack[s].j   = k;
      bt_stack[++s].i = k+1;
      bt_stack[s].j   = j1;
    } else {
#if 0
      /* Y shaped ML loops fon't work yet */
      if (dangle_model==3) {
        d5 = P->dangle5[tt][S1[j-1]];
        d3 = P->dangle3[tt][S1[i+1]];
        /* (i,j) must close a Y shaped ML loop with coax stacking */
        if (cij ==  fML[indx[j-2]+i+2] + mm + d3 + d5 + P->MLbase + P->MLbase) {
          i1 = i+2;
          j1 = j-2;
        } else if (cij ==  fML[indx[j-2]+i+1] + mm + d5 + P->MLbase)
          j1 = j-2;
        else if (cij ==  fML[indx[j-1]+i+2] + mm + d3 + P->MLbase)
          i1 = i+2;
        else /* last chance */
          if (cij != fML[indx[j-1]+i+1] + mm + P->MLbase)
            fprintf(stderr,  "backtracking failed in repeat");
        /* if we arrive here we can express cij via fML[i1,j1]+dangles */
        bt_stack[++s].i = i1;
        bt_stack[s].j   = j1;
      }
      else
#endif
        vrna_message_error("backtracking failed in repeat");
    }

    continue; /* this is a workarround to not accidentally proceed in the following block */

  repeat_gquad:
    /*
      now we do some fancy stuff to backtrace the stacksize and linker lengths
      of the g-quadruplex that should reside within position i,j
    */
    {
      int l[3], L, a;
      L = -1;

      get_gquad_pattern_mfe(S, i, j, P, &L, l);
      if(L != -1){
        /* fill the G's of the quadruplex into base_pair2 */
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
        goto repeat_gquad_exit;
      }
      vrna_message_error("backtracking failed in repeat_gquad");
    }
  repeat_gquad_exit:
    asm("nop");

  } /* end of infinite while loop */

  bp_stack[0].i = b;    /* save the total number of base pairs */
}

/*---------------------------------------------------------------------------*/


/*---------------------------------------------------------------------------*/

PUBLIC  void
vrna_update_fold_params(vrna_fold_compound *vc,
                        vrna_param_t *parameters){

  vrna_fold_compound *v;
  if(vc){
    v = vc;

    /* what about re-setting the backward compatibility compound here? */
    if(backward_compat_compound && backward_compat)
      vrna_free_fold_compound(backward_compat_compound);

    backward_compat_compound  = vc;
    backward_compat           = 0;
  } else if(backward_compat_compound && backward_compat){
    v = backward_compat_compound;
  } else
    return;

  if(v->params)
    free(v->params);
  if(parameters){
    v->params = get_parameter_copy(parameters);
  } else {
    vrna_md_t md;
    set_model_details(&md);
    v->params = get_scaled_parameters(temperature, md);
  }
}



/*###########################################*/
/*# deprecated functions below              #*/
/*###########################################*/

#ifdef  VRNA_BACKWARD_COMPAT

PUBLIC void
free_arrays(void){

  if(backward_compat_compound && backward_compat){
    vrna_free_fold_compound(backward_compat_compound);
    backward_compat_compound = NULL;
    backward_compat          = 0;
  }
}

PUBLIC float
fold_par( const char *string,
          char *structure,
          vrna_param_t *parameters,
          int is_constrained,
          int is_circular){

  return wrap_fold(string, structure, parameters, is_constrained, is_circular);

}

PUBLIC float
fold( const char *string,
      char *structure){

  return wrap_fold(string, structure, NULL, fold_constrained, 0);
}

PUBLIC float
circfold( const char *string,
          char *structure){

  return wrap_fold(string, structure, NULL, fold_constrained, 1);
}

PUBLIC void
initialize_fold(int length){

  /* DO NOTHING */
}

PUBLIC void
update_fold_params(void){

  vrna_update_fold_params(NULL, NULL);
}

PUBLIC void
update_fold_params_par(vrna_param_t *parameters){

  vrna_update_fold_params(NULL, parameters);
}

PUBLIC void
export_fold_arrays( int **f5_p,
                    int **c_p,
                    int **fML_p,
                    int **fM1_p,
                    int **indx_p,
                    char **ptype_p){

  wrap_array_export(f5_p,c_p,fML_p,fM1_p,indx_p,ptype_p);
}

PUBLIC void
export_fold_arrays_par( int **f5_p,
                        int **c_p,
                        int **fML_p,
                        int **fM1_p,
                        int **indx_p,
                        char **ptype_p,
                        vrna_param_t **P_p){

  wrap_array_export(f5_p,c_p,fML_p,fM1_p,indx_p,ptype_p);
  if(backward_compat_compound) *P_p  = backward_compat_compound->params;
}

PUBLIC void
export_circfold_arrays( int *Fc_p,
                        int *FcH_p,
                        int *FcI_p,
                        int *FcM_p,
                        int **fM2_p,
                        int **f5_p,
                        int **c_p,
                        int **fML_p,
                        int **fM1_p,
                        int **indx_p,
                        char **ptype_p){

  wrap_array_export(f5_p,c_p,fML_p,fM1_p,indx_p,ptype_p);
  wrap_array_export_circ(Fc_p,FcH_p,FcI_p,FcM_p,fM2_p);
}

PUBLIC void
export_circfold_arrays_par( int *Fc_p,
                            int *FcH_p,
                            int *FcI_p,
                            int *FcM_p,
                            int **fM2_p,
                            int **f5_p,
                            int **c_p,
                            int **fML_p,
                            int **fM1_p,
                            int **indx_p,
                            char **ptype_p,
                            vrna_param_t **P_p){

  wrap_array_export(f5_p,c_p,fML_p,fM1_p,indx_p,ptype_p);
  wrap_array_export_circ(Fc_p,FcH_p,FcI_p,FcM_p,fM2_p);
  if(backward_compat_compound) *P_p  = backward_compat_compound->params;
}

PUBLIC char *
backtrack_fold_from_pair( char *sequence,
                          int i,
                          int j){

  char          *structure  = NULL;
  unsigned int  length      = 0;
  bondT         *bp         = NULL;
  sect          bt_stack[MAXSECTORS]; /* stack of partial structures for backtracking */

  if(sequence){
    length = strlen(sequence);
    structure = (char *) vrna_alloc((length + 1)*sizeof(char));
    bp = (bondT *)vrna_alloc(sizeof(bondT) * (1+length/2));
  } else {
    vrna_message_error("backtrack_fold_from_pair@fold.c: no sequence given");
  }

  bt_stack[1].i  = i;
  bt_stack[1].j  = j;
  bt_stack[1].ml = 2;

  bp[0].i = 0; /* ??? this is set by backtrack anyway... */

  backtrack(backward_compat_compound, bp, bt_stack, 1);
  vrna_parenthesis_structure(structure, bp, length);

  /* backward compatibitlity stuff */
  if(base_pair) free(base_pair);
  base_pair = bp;

  return structure;
}

#define STACK_BULGE1      1       /* stacking energies for bulges of size 1 */
#define NEW_NINIO         1       /* new asymetry penalty */

PUBLIC int HairpinE(int size, int type, int si1, int sj1, const char *string) {
  vrna_param_t  *P  = backward_compat_compound->params;
  int energy;

  energy = (size <= 30) ? P->hairpin[size] :
    P->hairpin[30]+(int)(P->lxc*log((size)/30.));

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

/*---------------------------------------------------------------------------*/

PUBLIC int oldLoopEnergy(int i, int j, int p, int q, int type, int type_2) {

  vrna_param_t  *P  = backward_compat_compound->params;
  short   *S1 = backward_compat_compound->sequence_encoding;

  /* compute energy of degree 2 loop (stack bulge or interior) */
  int n1, n2, m, energy;
  n1 = p-i-1;
  n2 = j-q-1;

  if (n1>n2) { m=n1; n1=n2; n2=m; } /* so that n2>=n1 */

  if (n2 == 0)
    energy = P->stack[type][type_2];   /* stack */

  else if (n1==0) {                  /* bulge */
    energy = (n2<=MAXLOOP)?P->bulge[n2]:
      (P->bulge[30]+(int)(P->lxc*log(n2/30.)));

#if STACK_BULGE1
    if (n2==1) energy+=P->stack[type][type_2];
#endif
  } else {                           /* interior loop */

    if ((n1+n2==2)&&(james_rule))
      /* special case for loop size 2 */
      energy = P->int11[type][type_2][S1[i+1]][S1[j-1]];
    else {
      energy = (n1+n2<=MAXLOOP)?(P->internal_loop[n1+n2]):
        (P->internal_loop[30]+(int)(P->lxc*log((n1+n2)/30.)));

#if NEW_NINIO
      energy += MIN2(MAX_NINIO, (n2-n1)*P->ninio[2]);
#else
      m       = MIN2(4, n1);
      energy += MIN2(MAX_NINIO,((n2-n1)*P->ninio[m]));
#endif
      energy += P->mismatchI[type][S1[i+1]][S1[j-1]]+
        P->mismatchI[type_2][S1[q+1]][S1[p-1]];
    }
  }
  return energy;
}

/*--------------------------------------------------------------------------*/

PUBLIC int LoopEnergy(int n1, int n2, int type, int type_2,
                      int si1, int sj1, int sp1, int sq1) {

  vrna_param_t  *P  = backward_compat_compound->params;
  /* compute energy of degree 2 loop (stack bulge or interior) */
  int nl, ns, energy;

  if (n1>n2) { nl=n1; ns=n2;}
  else {nl=n2; ns=n1;}

  if (nl == 0)
    return P->stack[type][type_2];    /* stack */

  if (ns==0) {                       /* bulge */
    energy = (nl<=MAXLOOP)?P->bulge[nl]:
      (P->bulge[30]+(int)(P->lxc*log(nl/30.)));
    if (nl==1) energy += P->stack[type][type_2];
    else {
      if (type>2) energy += P->TerminalAU;
      if (type_2>2) energy += P->TerminalAU;
    }
    return energy;
  }
  else {                             /* interior loop */
    if (ns==1) {
      if (nl==1)                     /* 1x1 loop */
        return P->int11[type][type_2][si1][sj1];
      if (nl==2) {                   /* 2x1 loop */
        if (n1==1)
          energy = P->int21[type][type_2][si1][sq1][sj1];
        else
          energy = P->int21[type_2][type][sq1][si1][sp1];
        return energy;
      }
        else {  /* 1xn loop */
        energy = (nl+1<=MAXLOOP)?(P->internal_loop[nl+1]):
        (P->internal_loop[30]+(int)(P->lxc*log((nl+1)/30.)));
        energy += MIN2(MAX_NINIO, (nl-ns)*P->ninio[2]);
        energy += P->mismatch1nI[type][si1][sj1]+
        P->mismatch1nI[type_2][sq1][sp1];
        return energy;
        }
    }
    else if (ns==2) {
      if(nl==2)      {   /* 2x2 loop */
        return P->int22[type][type_2][si1][sp1][sq1][sj1];}
      else if (nl==3)  { /* 2x3 loop */
        energy = P->internal_loop[5]+P->ninio[2];
        energy += P->mismatch23I[type][si1][sj1]+
          P->mismatch23I[type_2][sq1][sp1];
        return energy;
      }

    }
    { /* generic interior loop (no else here!)*/
      energy = (n1+n2<=MAXLOOP)?(P->internal_loop[n1+n2]):
        (P->internal_loop[30]+(int)(P->lxc*log((n1+n2)/30.)));

      energy += MIN2(MAX_NINIO, (nl-ns)*P->ninio[2]);

      energy += P->mismatchI[type][si1][sj1]+
        P->mismatchI[type_2][sq1][sp1];
    }
  }
  return energy;
}

#endif
