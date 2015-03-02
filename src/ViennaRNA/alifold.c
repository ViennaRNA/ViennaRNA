/* Last changed Time-stamp: <2009-02-24 15:17:17 ivo> */
/*
                  minimum free energy folding
                  for a set of aligned sequences

                  c Ivo Hofacker

                  Vienna RNA package
*/

/**
*** \file alifold.c
**/


#include <config.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <string.h>
#include <limits.h>

#include "ViennaRNA/fold_vars.h"
#include "ViennaRNA/data_structures.h"
#include "ViennaRNA/fold.h"
#include "ViennaRNA/eval.h"
#include "ViennaRNA/utils.h"
#include "ViennaRNA/energy_par.h"
#include "ViennaRNA/params.h"
#include "ViennaRNA/ribo.h"
#include "ViennaRNA/gquad.h"
#include "ViennaRNA/alifold.h"
#include "ViennaRNA/aln_util.h"
#include "ViennaRNA/loop_energies.h"

#ifdef _OPENMP
#include <omp.h>
#endif

#define MAXSECTORS    500   /* dimension for a backtrack array */

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
PRIVATE vrna_fold_compound  *backward_compat_compound = NULL;
PRIVATE int                 backward_compat           = 0;

#ifdef _OPENMP

#pragma omp threadprivate(backward_compat_compound, backward_compat)

#endif

/*
#################################
# PRIVATE FUNCTION DECLARATIONS #
#################################
*/

PRIVATE int     fill_arrays(vrna_fold_compound *vc);
PRIVATE void    fill_arrays_circ(vrna_fold_compound *vc, sect bt_stack[], int *bt);
PRIVATE void    backtrack(vrna_fold_compound *vc, bondT *bp_stack, sect bt_stack[], int s);

PRIVATE float   wrap_alifold( const char **strings,
                              char *structure,
                              paramT *parameters,
                              int is_constrained,
                              int is_circular);

/*
#################################
# BEGIN OF FUNCTION DEFINITIONS #
#################################
*/

PRIVATE float
wrap_alifold( const char **strings,
              char *structure,
              paramT *parameters,
              int is_constrained,
              int is_circular){

  vrna_fold_compound  *vc;
  paramT              *P;

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

  vc = vrna_get_fold_compound_ali(strings, &(P->model_details), VRNA_OPTION_MFE);

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
                          | VRNA_CONSTRAINT_PIPE
                          | VRNA_CONSTRAINT_DOT
                          | VRNA_CONSTRAINT_X
                          | VRNA_CONSTRAINT_ANG_BRACK
                          | VRNA_CONSTRAINT_RND_BRACK;

    vrna_hc_add(vc, (const char *)structure, constraint_options);
  }

  if(backward_compat_compound && backward_compat)
    vrna_free_fold_compound(backward_compat_compound);

  backward_compat_compound  = vc;
  backward_compat           = 1;

  return vrna_ali_fold(vc, structure);
}



PUBLIC float
vrna_ali_fold(vrna_fold_compound *vc,
              char *structure){

  int  length, s, n_seq, energy;
  sect    bt_stack[MAXSECTORS]; /* stack of partial structures for backtracking */
  bondT *bp;

  length      = vc->length;
  n_seq       = vc->n_seq;
  s           = 0;

  energy = fill_arrays(vc);
  if(vc->params->model_details.circ){
    fill_arrays_circ(vc, bt_stack, &s);
    energy = vc->matrices->Fc;
  }

  if(structure && vc->params->model_details.backtrack){
    bp = (bondT *)space(sizeof(bondT) * (4*(1+length/2))); /* add a guess of how many G's may be involved in a G quadruplex */

    backtrack(vc, bp, bt_stack, s);

    vrna_parenthesis_structure(structure, bp, length);

    /*
    *  Backward compatibility:
    *  This block may be removed if deprecated functions
    *  relying on the global variable "base_pair" vanishs from within the package!
    */
    {
      if(base_pair) free(base_pair);
      base_pair = bp;
    }
  }

  if (vc->params->model_details.backtrack_type=='C')
    return (float) vc->matrices->c[vc->jindx[length]+1]/(n_seq*100.);
  else if (vc->params->model_details.backtrack_type=='M')
    return (float) vc->matrices->fML[vc->jindx[length]+1]/(n_seq*100.);
  else
    return (float) energy/(n_seq*100.);
}


/**
*** the actual forward recursion to fill the energy arrays
**/
PRIVATE int
fill_arrays(vrna_fold_compound *vc){

  int   i, j, k, p, q, energy, new_c;
  int   decomp, MLenergy, new_fML;
  int   s, *type, type_2, tt;
  int   *cc;        /* linear array for calculating canonical structures */
  int   *cc1;       /*   "     "        */
  int   *Fmi;       /* holds row i of fML (avoids jumps in memory) */
  int   *DMLi;      /* DMLi[j] holds MIN(fML[i,k]+fML[k+1,j])  */
  int   *DMLi1;     /*             MIN(fML[i+1,k]+fML[k+1,j])  */
  int   *DMLi2;     /*             MIN(fML[i+2,k]+fML[k+1,j])  */


  int             n_seq         = vc->n_seq;
  int             length        = vc->length;
  short           **S           = vc->S;                                                                   
  short           **S5          = vc->S5;     /*S5[s][i] holds next base 5' of i in sequence s*/            
  short           **S3          = vc->S3;     /*Sl[s][i] holds next base 3' of i in sequence s*/            
  char            **Ss          = vc->Ss;                                                                   
  unsigned short  **a2s         = vc->a2s;                                                                   
  paramT          *P            = vc->params;                                                                
  vrna_md_t       *md           = &(P->model_details);
  int             *indx         = vc->jindx;     /* index for moving in the triangle matrices c[] and fMl[]*/
  int             *c            = vc->matrices->c;     /* energy array, given that i-j pair */                 
  int             *f5           = vc->matrices->f5;     /* energy of 5' end */                                  
  int             *fML          = vc->matrices->fML;     /* multi-loop auxiliary energy array */                 
  int             *ggg          = vc->matrices->ggg;
  int             *pscore       = vc->pscore;     /* precomputed array of pair types */                      
  short           *S_cons       = vc->S_cons;
  int             *rtype        = &(md->rtype[0]);
  int             dangle_model  = md->dangles;

  vrna_hc_t       *hc           = vc->hc;
  vrna_sc_t       **sc          = vc->scs;

  char              *hard_constraints = hc->matrix;

  type  = (int *) space(n_seq*sizeof(int));
  cc    = (int *) space(sizeof(int)*(length+2));
  cc1   = (int *) space(sizeof(int)*(length+2));
  Fmi   = (int *) space(sizeof(int)*(length+1));
  DMLi  = (int *) space(sizeof(int)*(length+1));
  DMLi1 = (int *) space(sizeof(int)*(length+1));
  DMLi2 = (int *) space(sizeof(int)*(length+1));
  
  /* init energies */

  int max_bpspan = (md->max_bp_span > 0) ? md->max_bp_span : length;

  for (j=1; j<=length; j++){
    Fmi[j]=DMLi[j]=DMLi1[j]=DMLi2[j]=INF;
    for (i=(j>TURN?(j-TURN):1); i<j; i++) {
      c[indx[j]+i] = fML[indx[j]+i] = INF;
    }
  }

  /* begin recursions */
  for (i = length-TURN-1; i >= 1; i--) { /* i,j in [1..length] */
    for (j = i+TURN+1; j <= length; j++) {
      int ij, psc, l1, maxq, minq, c0;
      ij = indx[j]+i;

      for (s=0; s<n_seq; s++) {
        type[s] = md->pair[S[s][i]][S[s][j]];
        if (type[s]==0) type[s]=7;
      }

      psc = pscore[indx[j]+i];
      if (hard_constraints[ij]) {   /* a pair to consider */
        int stackEnergy = INF;
        /* hairpin ----------------------------------------------*/
        new_c = E_hp_loop_ali(i, j, vc);

        /*--------------------------------------------------------
          check for elementary structures involving more than one
          closing pair.
          --------------------------------------------------------*/

        if(hard_constraints[ij] & VRNA_HC_CONTEXT_INT_LOOP){
          for (p = i+1; p <= MIN2(j-2-TURN,i+MAXLOOP+1) ; p++) {
            if(hc->up_int[i+1] < p - i - 1)
              break;

            minq = j-i+p-MAXLOOP-2;
            if (minq<p+1+TURN) minq = p+1+TURN;
            for (q = minq; q < j; q++) {
              if(!(hard_constraints[indx[q]+p] & VRNA_HC_CONTEXT_INT_LOOP_ENC))
                continue;
              if(hc->up_int[q+1] < j - q - 1)
                continue;

              for (energy = s=0; s<n_seq; s++) {
                int u1 = a2s[s][p-1]-a2s[s][i];
                int u2 = a2s[s][j-1]-a2s[s][q];
                type_2 = md->pair[S[s][q]][S[s][p]]; /* q,p not p,q! */
                if (type_2 == 0) type_2 = 7;
                energy += E_IntLoop(u1, u2, type[s], type_2,
                                     S3[s][i], S5[s][j],
                                     S5[s][p], S3[s][q], P);
              }

              if(sc){
                for(s = 0; s < n_seq; s++){
                  if(sc[s]){
                    int u1 = a2s[s][p-1]-a2s[s][i];
                    int u2 = a2s[s][j-1]-a2s[s][q];
/*
                    int u1 = p - i - 1;
                    int u2 = j - q - 1;
*/
                    if(sc[s]->en_basepair)
                      energy +=   sc[s]->en_basepair[indx[j] + i];

                    if(sc[s]->free_energies)
                      energy +=   sc[s]->free_energies[a2s[s][i]+1][u1]
                                + sc[s]->free_energies[a2s[s][q]+1][u2];

                    if(u1 + u2 == 0)
                      if(sc[s]->en_stack)
                        energy +=   sc[s]->en_stack[i]
                                  + sc[s]->en_stack[p]
                                  + sc[s]->en_stack[q]
                                  + sc[s]->en_stack[j];

                  }
                }
              }

              if ((p==i+1)&&(j==q+1)){
                stackEnergy = energy; /* remember stack energy */
              }
              new_c = MIN2(new_c, energy + c[indx[q]+p]);
            } /* end q-loop */
          } /* end p-loop */

          if(md->gquad){
            decomp = 0;
            for(s=0;s<n_seq;s++){
              tt = type[s];
              if(dangle_model == 2)
                decomp += P->mismatchI[tt][S3[s][i]][S5[s][j]];
              if(tt > 2)
                decomp += P->TerminalAU;
            }
            for(p = i + 2; p < j - VRNA_GQUAD_MIN_BOX_SIZE; p++){
              l1    = p - i - 1;
              if(l1>MAXLOOP) break;
              if(S_cons[p] != 3) continue;
              minq  = j - i + p - MAXLOOP - 2;
              c0    = p + VRNA_GQUAD_MIN_BOX_SIZE - 1;
              minq  = MAX2(c0, minq);
              c0    = j - 1;
              maxq  = p + VRNA_GQUAD_MAX_BOX_SIZE + 1;
              maxq  = MIN2(c0, maxq);
              for(q = minq; q < maxq; q++){
                if(S_cons[q] != 3) continue;
                c0    = decomp + ggg[indx[q] + p] + n_seq * P->internal_loop[l1 + j - q - 1];
                new_c = MIN2(new_c, c0);
              }
            }

            p = i + 1;
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
                  c0  = decomp + ggg[indx[q] + p] + n_seq * P->internal_loop[j - q - 1];
                  new_c   = MIN2(new_c, c0);
                }
              }
            }
            q = j - 1;
            if(S_cons[q] == 3)
              for(p = i + 4; p < j - VRNA_GQUAD_MIN_BOX_SIZE; p++){
                l1    = p - i - 1;
                if(l1>MAXLOOP) break;
                if(S_cons[p] != 3) continue;
                c0  = decomp + ggg[indx[q] + p] + n_seq * P->internal_loop[l1];
                new_c   = MIN2(new_c, c0);
              }
          }
        } /* end if interior loop */

        /* multi-loop decomposition ------------------------*/
        if(hard_constraints[ij] & VRNA_HC_CONTEXT_MB_LOOP){
          decomp = DMLi1[j-1];
          if(dangle_model){
            for(s=0; s<n_seq; s++){
              tt = rtype[type[s]];
              decomp += E_MLstem(tt, S5[s][j], S3[s][i], P);
            }
          }
          else{
            for(s=0; s<n_seq; s++){
              tt = rtype[type[s]];
              decomp += E_MLstem(tt, -1, -1, P);
            }
          }
          if(sc)
            for(s=0; s<n_seq; s++){
              if(sc[s]){
                if(sc[s]->en_basepair)
                  decomp += sc[s]->en_basepair[indx[j] + i];
              }
            }
          MLenergy = decomp + n_seq*P->MLclosing;
          new_c = MIN2(new_c, MLenergy);
        }


        new_c = MIN2(new_c, cc1[j-1]+stackEnergy);

        cc[j] = new_c - psc; /* add covariance bonnus/penalty */
        if (md->noLP)
          c[ij] = cc1[j-1]+stackEnergy-psc;
        else
          c[ij] = cc[j];

      } /* end >> if (pair) << */

      else c[ij] = INF;
      /* done with c[i,j], now compute fML[i,j] */
      /* free ends ? -----------------------------------------*/
      new_fML = INF;

      if(hc->up_ml[i]){
        energy = fML[ij+1] + n_seq * P->MLbase;
        if(sc)
          for(s=0; s < n_seq; s++){
            if(sc[s]){
              if(sc[s]->free_energies)
                energy += sc[s]->free_energies[a2s[s][i]][1];
            }
          }
        new_fML = MIN2(new_fML, energy);
      }

      if(hc->up_ml[j]){
        energy = fML[indx[j-1]+i] + n_seq * P->MLbase;
        if(sc)
          for(s=0;s < n_seq; s++){
            if(sc[s]){
              if(sc[s]->free_energies)
                energy += sc[s]->free_energies[a2s[s][j]][1];
            }
          }
        new_fML = MIN2(new_fML, energy);
      }

      if(hard_constraints[ij] & VRNA_HC_CONTEXT_MB_LOOP_ENC){
        energy = c[ij];
        if(dangle_model){
          for (s=0; s<n_seq; s++) {
            energy += E_MLstem(type[s], S5[s][i], S3[s][j], P);
          }
        }
        else{
          for (s=0; s<n_seq; s++) {
            energy += E_MLstem(type[s], -1, -1, P);
          }
        }
        new_fML = MIN2(energy, new_fML);

        if(md->gquad){
          decomp = ggg[indx[j] + i] + n_seq * E_MLstem(0, -1, -1, P);
          new_fML = MIN2(new_fML, decomp);
        }
      }


      /* modular decomposition -------------------------------*/
      for (decomp = INF, k = i+1+TURN; k <= j-2-TURN; k++)
        decomp = MIN2(decomp, Fmi[k]+fML[indx[j]+k+1]);

      DMLi[j] = decomp;               /* store for use in ML decompositon */
      new_fML = MIN2(new_fML,decomp);

      /* coaxial stacking deleted */

      fML[ij] = Fmi[j] = new_fML;     /* substring energy */
    } /* END for j */

    {
      int *FF; /* rotate the auxilliary arrays */
      FF = DMLi2; DMLi2 = DMLi1; DMLi1 = DMLi; DMLi = FF;
      FF = cc1; cc1=cc; cc=FF;
      for (j=1; j<=length; j++) {cc[j]=Fmi[j]=DMLi[j]=INF; }
    }
  } /* END for i */
  /* calculate energies of 5' and 3' fragments */

  for(j = 1; j <= TURN + 1; j++){
    if(hc->up_ext[j]){
      energy = f5[j-1];
      if((energy < INF) && sc)
        for(s=0;s < n_seq; s++){
          if(sc[s]){
            if(sc[s]->free_energies)
              energy += sc[s]->free_energies[a2s[s][j]][1];
          }
        }
    } else {
      energy = INF;
    }
    f5[j] = energy;
  }

  switch(dangle_model){
    case 0:   for(j = TURN + 2; j <= length; j++){
                f5[j] = INF;

                if(hc->up_ext[j]){
                  energy = f5[j-1];
                  if((energy < INF) && sc)
                    for(s=0; s < n_seq; s++){
                      if(sc[s]){
                        if(sc[s]->free_energies)
                          energy += sc[s]->free_energies[a2s[s][j]][1];
                      }
                    }
                  f5[j] = MIN2(f5[j], energy);
                }

                if (hard_constraints[indx[j]+1] & VRNA_HC_CONTEXT_EXT_LOOP){
                  if(c[indx[j]+1] < INF){
                    energy = c[indx[j]+1];
                    for(s = 0; s < n_seq; s++){
                      tt = md->pair[S[s][1]][S[s][j]];
                      if(tt==0) tt=7;
                      energy += E_ExtLoop(tt, -1, -1, P);
                    }
                    f5[j] = MIN2(f5[j], energy);
                  }

                  if(md->gquad){
                    if(ggg[indx[j]+1] < INF)
                      f5[j] = MIN2(f5[j], ggg[indx[j]+1]);
                  }
                }

                for(i = j - TURN - 1; i > 1; i--){
                  if(hard_constraints[indx[j]+i] & VRNA_HC_CONTEXT_EXT_LOOP){
                    if(c[indx[j]+i]<INF){
                      energy = f5[i-1] + c[indx[j]+i];
                      for(s = 0; s < n_seq; s++){
                        tt = md->pair[S[s][i]][S[s][j]];
                        if(tt==0) tt=7;
                        energy += E_ExtLoop(tt, -1, -1, P);
                      }
                      f5[j] = MIN2(f5[j], energy);
                    }

                    if(md->gquad){
                      if(ggg[indx[j]+i] < INF)
                        f5[j] = MIN2(f5[j], f5[i-1] + ggg[indx[j]+i]);
                    }
                  }
                }
              }
              break;

    default:  for(j = TURN + 2; j <= length; j++){
                f5[j] = INF;

                if(hc->up_ext[j]){
                  energy = f5[j-1];
                  if((energy < INF) && sc)
                    for(s=0; s < n_seq; s++){
                      if(sc[s]){
                        if(sc[s]->free_energies)
                          energy += sc[s]->free_energies[a2s[s][j]][1];
                      }
                    }
                  f5[j] = MIN2(f5[j], energy);
                }

                if(hard_constraints[indx[j]+1] & VRNA_HC_CONTEXT_EXT_LOOP){
                  if (c[indx[j]+1]<INF) {
                    energy = c[indx[j]+1];
                    for(s = 0; s < n_seq; s++){
                      tt = md->pair[S[s][1]][S[s][j]];
                      if(tt==0) tt=7;
                      energy += E_ExtLoop(tt, -1, (j<length) ? S3[s][j] : -1, P);
                    }
                    f5[j] = MIN2(f5[j], energy);
                  }

                  if(md->gquad){
                    if(ggg[indx[j]+1] < INF)
                      f5[j] = MIN2(f5[j], ggg[indx[j]+1]);
                  }
                }

                for(i = j - TURN - 1; i > 1; i--){
                  if(hard_constraints[indx[j]+i] & VRNA_HC_CONTEXT_EXT_LOOP){
                    if (c[indx[j]+i]<INF) {
                      energy = f5[i-1] + c[indx[j]+i];
                      for(s = 0; s < n_seq; s++){
                        tt = md->pair[S[s][i]][S[s][j]];
                        if(tt==0) tt=7;
                        energy += E_ExtLoop(tt, S5[s][i], (j < length) ? S3[s][j] : -1, P);
                      }
                      f5[j] = MIN2(f5[j], energy);
                    }

                    if(md->gquad){
                      if(ggg[indx[j]+i] < INF)
                        f5[j] = MIN2(f5[j], f5[i-1] + ggg[indx[j]+i]);
                    }
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

/**
*** backtrack in the energy matrices to obtain a structure with MFE
**/
PRIVATE void
backtrack(vrna_fold_compound *vc,
          bondT *bp_stack,
          sect bt_stack[],
          int s) {

  /*------------------------------------------------------------------
    trace back through the "c", "f5" and "fML" arrays to get the
    base pairing list. No search for equivalent structures is done.
    This inverts the folding procedure, hence it's very fast.
    ------------------------------------------------------------------*/
   /* normally s=0.
     If s>0 then s items have been already pushed onto the sector stack */

  int  i, j, k, p, q, energy, en, c0, l1, minq, maxq, type_2, tt, mm, b=0, cov_en = 0, *type;

  int             n_seq         = vc->n_seq;
  int             length        = vc->length;
  short           **S           = vc->S;                                                                   
  short           **S5          = vc->S5;     /*S5[s][i] holds next base 5' of i in sequence s*/            
  short           **S3          = vc->S3;     /*Sl[s][i] holds next base 3' of i in sequence s*/            
  char            **Ss          = vc->Ss;                                                                   
  unsigned short  **a2s         = vc->a2s;                                                                   
  paramT          *P            = vc->params;                                                                
  vrna_md_t       *md           = &(P->model_details);
  int             *indx         = vc->jindx;     /* index for moving in the triangle matrices c[] and fMl[]*/
  int             *c            = vc->matrices->c;     /* energy array, given that i-j pair */                 
  int             *f5           = vc->matrices->f5;     /* energy of 5' end */                                  
  int             *fML          = vc->matrices->fML;     /* multi-loop auxiliary energy array */                 
  int             *pscore       = vc->pscore;     /* precomputed array of pair types */                      
  int             *ggg          = vc->matrices->ggg;
  short           *S_cons       = vc->S_cons;
  int             *rtype        = &(md->rtype[0]);
  int             dangle_model  = md->dangles;
  int             with_gquad    = md->gquad;

  vrna_hc_t       *hc           = vc->hc;
  vrna_sc_t       **sc          = vc->scs;

  type = (int *) space(n_seq*sizeof(int));

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
      goto repeat1;
    }

    if (j < i+TURN+1) continue; /* no more pairs in this interval */

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
          if(sc[ss]->free_energies)
            fi += sc[ss]->free_energies[a2s[ss][j]][1];
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
                  for (i=j-TURN-1,traced=0; i>=1; i--) {
                    int en;
                    jj = i-1;

                    if (hc->matrix[indx[j] + i] & VRNA_HC_CONTEXT_EXT_LOOP){
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
                  for (i=j-TURN-1,traced=0; i>=1; i--) {
                    int en;
                    jj = i-1;
                    if (hc->matrix[indx[j] + i] & VRNA_HC_CONTEXT_EXT_LOOP){
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

      if (!traced) nrerror("backtrack failed in f5");
      /* push back the remaining f5 portion */
      bt_stack[++s].i = 1;
      bt_stack[s].j   = jj;
      bt_stack[s].ml  = ml;

      /* trace back the base pair found */
      j=traced;

      if(with_gquad && gq){
        /* goto backtrace of gquadruplex */
        goto repeat_gquad;
      }

      bp_stack[++b].i = i;
      bp_stack[b].j   = j;
      cov_en += pscore[indx[j]+i];
      goto repeat1;
    }
    else { /* trace back in fML array */
      if(hc->up_ml[i]){
        en = fML[indx[j]+i+1] + n_seq * P->MLbase;

        if(sc)
          for(ss = 0; ss < n_seq; ss++)
            if(sc[ss]){
              if(sc[ss]->free_energies)
                en += sc[ss]->free_energies[a2s[ss][i]][1];
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
          goto repeat_gquad;
        }
      }

      if(hc->matrix[indx[j] + i] & VRNA_HC_CONTEXT_MB_LOOP_ENC){
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
          goto repeat1;
        }
      }

      for (k = i+1+TURN; k <= j-2-TURN; k++)
        if (fij == (fML[indx[k]+i]+fML[indx[j]+k+1]))
          break;

      bt_stack[++s].i = i;
      bt_stack[s].j   = k;
      bt_stack[s].ml  = ml;
      bt_stack[++s].i = k+1;
      bt_stack[s].j   = j;
      bt_stack[s].ml  = ml;

      if (k>j-2-TURN) nrerror("backtrack failed in fML");
      continue;
    }

  repeat1:

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
            if(sc[ss]->en_basepair)
              cij -= sc[s]->en_basepair[indx[j] + i];
          }
        }
        cij += pscore[indx[j]+i];
        bp_stack[++b].i = i+1;
        bp_stack[b].j   = j-1;
        cov_en += pscore[indx[j-1]+i+1];
        i++; j--;
        canonical=0;
        goto repeat1;
      }
    canonical = 1;
    cij += pscore[indx[j]+i];

    {
      int cc = E_hp_loop_ali(i, j, vc);

      if (cij == cc) /* found hairpin */
        continue;
    }
    for (p = i+1; p <= MIN2(j-2-TURN,i+MAXLOOP+1); p++) {
      minq = j-i+p-MAXLOOP-2;
      if (minq<p+1+TURN) minq = p+1+TURN;
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
                if(sc[ss]->en_stack)
                  energy +=   sc[ss]->en_stack[i]
                            + sc[ss]->en_stack[p]
                            + sc[ss]->en_stack[q]
                            + sc[ss]->en_stack[j];

              if(sc[ss]->en_basepair)
                energy += sc[ss]->en_basepair[indx[j] + i];

              if(sc[ss]->free_energies)
                energy +=   sc[ss]->free_energies[a2s[ss][i] + 1][u1]
                          + sc[ss]->free_energies[a2s[ss][q] + 1][u2];
            }

        traced = (cij == energy+c[indx[q]+p]);
        if (traced) {
          bp_stack[++b].i = p;
          bp_stack[b].j   = q;
          cov_en += pscore[indx[q]+p];
          i = p, j = q;
          goto repeat1;
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
      for(s=0;s<n_seq;s++){
        tt = type[s];
        if(tt == 0) tt = 7;
        if(dangle_model == 2)
          mm += P->mismatchI[tt][S3[s][i]][S5[s][j]];
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
            goto repeat_gquad;
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
              goto repeat_gquad;
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
            goto repeat_gquad;
          }
        }
    }

    if(hc->matrix[indx[j] + i] & VRNA_HC_CONTEXT_MB_LOOP){
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
            if(sc[ss]->en_basepair)
              mm += sc[ss]->en_basepair[indx[j] + i];
          }

      bt_stack[s+1].ml  = bt_stack[s+2].ml = 1;

      for (k = i1+TURN+1; k < j1-TURN-1; k++){
        if(cij == fML[indx[k]+i1] + fML[indx[j1]+k+1] + mm) break;
      }

      if (k<=j-3-TURN) { /* found the decomposition */
        bt_stack[++s].i = i1;
        bt_stack[s].j   = k;
        bt_stack[++s].i = k+1;
        bt_stack[s].j   = j1;
      } else {
          nrerror("backtracking failed in repeat");
      }
    } else
      nrerror("backtracking failed in repeat");

    continue; /* this is a workarround to not accidentally proceed in the following block */

  repeat_gquad:
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
            }
          }
        }
      }
      nrerror("backtracking failed in repeat_gquad");
    }
  repeat_gquad_exit:
    asm("nop");

  }

  /* fprintf(stderr, "covariance energy %6.2f\n", cov_en/100.);  */

  bp_stack[0].i = b;    /* save the total number of base pairs */
  free(type);
}

/*###########################################*/
/*# deprecated functions below              #*/
/*###########################################*/

PUBLIC void
free_alifold_arrays(void){

  if(backward_compat_compound && backward_compat){
    vrna_free_fold_compound(backward_compat_compound);
    backward_compat_compound  = NULL;
    backward_compat           = 0;
  }
}

PUBLIC float
alifold(const char **strings,
        char *structure){

  return wrap_alifold(strings, structure, NULL, fold_constrained, 0);
}

PUBLIC float circalifold( const char **strings,
                          char *structure) {

  return wrap_alifold(strings, structure, NULL, fold_constrained, 1);
}

PUBLIC void 
update_alifold_params(void){

  vrna_fold_compound *v;

  if(backward_compat_compound && backward_compat){
    v = backward_compat_compound;

    if(v->params)
      free(v->params);

    vrna_md_t md;
    set_model_details(&md);
    v->params = vrna_get_energy_contributions(md);
  }
}

PUBLIC float
energy_of_ali_gquad_structure(const char **sequences,
                              const char *structure,
                              int n_seq,
                              float *energy){

  unsigned int  n;
  short         *pt;
  int           *loop_idx;

  if(sequences[0] != NULL){
    
    vrna_fold_compound  *vc;

    vrna_md_t md;
    set_model_details(&md);
    md.gquad = 1;

    vc = vrna_get_fold_compound_ali(sequences, &md, VRNA_OPTION_MFE | VRNA_OPTION_EVAL_ONLY);

    energy[0] = vrna_eval_structure(vc, structure);
    energy[1] = vrna_eval_covar_structure(vc, structure);

    vrna_free_fold_compound(vc);
  }
  else nrerror("energy_of_alistruct(): no sequences in alignment!");

  return energy[0];

}

PUBLIC  float
energy_of_alistruct(const char **sequences,
                    const char *structure,
                    int n_seq,
                    float *energy){

  short         *pt;

  if(sequences[0] != NULL){
    vrna_fold_compound  *vc;

    vrna_md_t md;
    set_model_details(&md);

    vc = vrna_get_fold_compound_ali(sequences, &md, VRNA_OPTION_MFE | VRNA_OPTION_EVAL_ONLY);

    energy[0] = vrna_eval_structure(vc, structure);
    energy[1] = vrna_eval_covar_structure(vc, structure);

    vrna_free_fold_compound(vc);
  }
  else nrerror("energy_of_alistruct(): no sequences in alignment!");

  return energy[0];
}
