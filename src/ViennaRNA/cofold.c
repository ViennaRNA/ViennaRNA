/* Last changed Time-stamp: <2008-12-03 17:44:38 ivo> */
/*
                  minimum free energy
                  RNA secondary structure prediction

                  c Ivo Hofacker, Chrisoph Flamm
                  original implementation by
                  Walter Fontana

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
#include "ViennaRNA/fold_vars.h"
#include "ViennaRNA/params.h"
#include "ViennaRNA/subopt.h"
#include "ViennaRNA/fold.h"
#include "ViennaRNA/loop_energies.h"
#include "ViennaRNA/gquad.h"
#include "ViennaRNA/cofold.h"

#ifdef _OPENMP
#include <omp.h>
#endif

#define MAXSECTORS        500     /* dimension for a backtrack array */
#undef  TURN
#define TURN              0       /* reset minimal base pair span for intermolecular pairings */
#define TURN2             3       /* used by zukersubopt */
#define SAME_STRAND(I,J)  (((I)>=cut_point)||((J)<cut_point))

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

PRIVATE float   mfe1, mfe2;       /* minimum free energies of the monomers */

#ifdef _OPENMP

#pragma omp threadprivate(mfe1, mfe2, backward_compat_compound)

#endif

/*
#################################
# PRIVATE FUNCTION DECLARATIONS #
#################################
*/

PRIVATE void  backtrack(sect bt_stack[], bondT *bp_list, vrna_fold_compound *vc);
PRIVATE int   fill_arrays(vrna_fold_compound *vc, int zuker);
PRIVATE void  free_end(int *array, int i, int start, vrna_fold_compound *vc);

/* wrappers for old API compatibility */
PRIVATE void      wrap_array_export(int **f5_p,int **c_p,int **fML_p,int **fM1_p,int **fc_p,int **indx_p,char **ptype_p);
PRIVATE float     wrap_cofold(const char *string,char *structure,paramT *parameters,int is_constrained);
PRIVATE SOLUTION *wrap_zukersubopt( const char *string,paramT *parameters);

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
                  int **fc_p,
                  int **indx_p,
                  char **ptype_p){

  /* make the DP arrays available to routines such as subopt() */
  if(backward_compat_compound){
    *f5_p     = backward_compat_compound->matrices->f5;
    *c_p      = backward_compat_compound->matrices->c;
    *fML_p    = backward_compat_compound->matrices->fML;
    *fM1_p    = backward_compat_compound->matrices->fM1;
    *fc_p     = backward_compat_compound->matrices->fc;
    *indx_p   = backward_compat_compound->jindx;
    *ptype_p  = backward_compat_compound->ptype;
  }
}

/*--------------------------------------------------------------------------*/

PRIVATE float
wrap_cofold(const char *string,
            char *structure,
            paramT *parameters,
            int is_constrained){

  unsigned int        length;
  char                *seq;
  vrna_fold_compound  *vc;
  paramT              *P;

  vc      = NULL;
  length  = strlen(string);

#ifdef _OPENMP
/* Explicitly turn off dynamic threads */
  omp_set_dynamic(0);
#endif

  /* we need the parameter structure for hard constraints */
  if(parameters)
    P = get_parameter_copy(parameters);
  else{
    model_detailsT md;
    set_model_details(&md);
    P = get_scaled_parameters(temperature, md);
  }
  P->model_details.min_loop_size = TURN;  /* set min loop length to 0 */

  /* dirty hack to reinsert the '&' according to the global variable 'cut_point' */
  seq = (char *)space(sizeof(char) * (length + 2));
  if(cut_point > -1){
    int i;
    for(i = 0; i < cut_point-1; i++)
      seq[i] = string[i];
    seq[i] = '&';
    for(;i<(int)length;i++)
      seq[i+1] = string[i];
  } else { /* this ensures the allocation of all cofold matrices via vrna_get_fold_compound */
    seq[0] = '&';
    strcat(seq + 1, string);
  }

  /* get compound structure */
  vc = vrna_get_fold_compound(seq, &(P->model_details), VRNA_OPTION_MFE | VRNA_OPTION_HYBRID);

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
                          | VRNA_CONSTRAINT_RND_BRACK
                          | VRNA_CONSTRAINT_INTRAMOLECULAR
                          | VRNA_CONSTRAINT_INTERMOLECULAR;

    vrna_hc_add(vc, (const char *)structure, constraint_options);
  }

  if(backward_compat_compound)
    destroy_fold_compound(backward_compat_compound);

  backward_compat_compound = vc;

  /* cleanup */
  free(seq);

  return vrna_cofold(vc, structure);
}

PUBLIC float
vrna_cofold(vrna_fold_compound  *vc,
            char                *structure){

  int     length, energy;
  sect    bt_stack[MAXSECTORS]; /* stack of partial structures for backtracking */
  bondT   *bp;

  length = (int) vc->length;

  vc->sequence_encoding[0] = vc->sequence_encoding2[0]; /* store length at pos. 0 in S1 too */

  energy = fill_arrays(vc, 0);

  bp = (bondT *)space(sizeof(bondT) * (1 + length/2 + 200)); /* add a guess of how many G's may be involved in a G quadruplex */

  backtrack(bt_stack, bp, vc);

  vrna_parenthesis_structure(structure, bp, length);

  /*
  *  Backward compatibility:
  *  This block may be removed if deprecated functions
  *  relying on the global variable "base_pair" vanish from within the package!
  */
  if(base_pair) free(base_pair);
  base_pair = bp;
  /*
  {
    if(base_pair) free(base_pair);
    base_pair = (bondT *)space(sizeof(bondT) * (1+length/2));
    memcpy(base_pair, base_pair2, sizeof(bondT) * (1+length/2));
  }
  */

  if (vc->params->model_details.backtrack_type=='C')
    return (float) vc->matrices->c[vc->jindx[length]+1]/100.;
  else if (vc->params->model_details.backtrack_type=='M')
    return (float) vc->matrices->fML[vc->jindx[length]+1]/100.;
  else
    return (float) energy/100.;
}

PRIVATE int
fill_arrays(vrna_fold_compound  *vc,
            int                 zuker){

  /* fill "c", "fML" and "f5" arrays and return  optimal energy */

  int   i, j, k, length, energy;
  int   decomp, new_fML, cut_point, uniq_ML;
  int   no_close, type, type_2, tt, maxj, *indx;
  int   *my_f5, *my_c, *my_fML, *my_fM1, *my_fc, *my_ggg;
  int               *cc, *cc1;  /* auxilary arrays for canonical structures     */
  int               *Fmi;       /* holds row i of fML (avoids jumps in memory)  */
  int               *DMLi;      /* DMLi[j] holds  MIN(fML[i,k]+fML[k+1,j])      */
  int               *DMLi1;     /*                MIN(fML[i+1,k]+fML[k+1,j])    */
  int               *DMLi2;     /*                MIN(fML[i+2,k]+fML[k+1,j])    */
  int   *hc_up_ext, *hc_up_ml;

  int   dangle_model, noGUclosure, with_gquad, noLP, hc_decompose;
  int   *rtype;
  char              *ptype, *hard_constraints;
  short             *S1;
  paramT            *P;
  mfe_matricesT     *matrices;
  hard_constraintT  *hc;
  soft_constraintT  *sc;

  length            = (int)vc->length;
  ptype             = vc->ptype;
  indx              = vc->jindx;
  P                 = vc->params;
  S1                = vc->sequence_encoding;
  dangle_model      = P->model_details.dangles;
  noGUclosure       = P->model_details.noGUclosure;
  noLP              = P->model_details.noLP;
  uniq_ML           = P->model_details.uniq_ML;
  with_gquad        = P->model_details.gquad;
  rtype             = &(P->model_details.rtype[0]);
  hc                = vc->hc;
  hard_constraints  = hc->matrix;
  sc                = vc->sc;
  matrices          = vc->matrices;
  my_f5             = matrices->f5;
  my_c              = matrices->c;
  my_fML            = matrices->fML;
  my_fM1            = matrices->fM1;
  my_fc             = matrices->fc;
  my_ggg            = matrices->ggg;
  cut_point         = vc->cutpoint;

  hc_up_ext         = hc->up_ext;
  hc_up_ml          = hc->up_ml;

  /* allocate memory for all helper arrays */
  cc    = (int *) space(sizeof(int)*(length + 2));
  cc1   = (int *) space(sizeof(int)*(length + 2));
  Fmi   = (int *) space(sizeof(int)*(length + 1));
  DMLi  = (int *) space(sizeof(int)*(length + 1));
  DMLi1 = (int *) space(sizeof(int)*(length + 1));
  DMLi2 = (int *) space(sizeof(int)*(length + 1));


  for (j=1; j<=length; j++) {
    Fmi[j]=DMLi[j]=DMLi1[j]=DMLi2[j]=INF;
    my_fc[j]=0;
  }

  for (j = 1; j<=length; j++)
    for (i=1; i<=j; i++) {
      my_c[indx[j]+i] = my_fML[indx[j]+i] = INF;
      if (uniq_ML) my_fM1[indx[j]+i] = INF;
    }

  for (i = length-TURN-1; i >= 1; i--) { /* i,j in [1..length] */

    maxj=(zuker)? (MIN2(i+cut_point-1,length)):length;
    for (j = i+TURN+1; j <= maxj; j++) {
      int p, q, ij;
      ij            = indx[j]+i;
      type          = ptype[ij];
      hc_decompose  = hard_constraints[ij];
      energy        = INF;

      no_close = (((type==3)||(type==4))&&noGUclosure);

      if (hc_decompose) {   /* we have a pair */
        int new_c=INF, stackEnergy=INF;
        short si, sj;
        si  = SAME_STRAND(i, i+1) ? S1[i+1] : -1;
        sj  = SAME_STRAND(j-1, j) ? S1[j-1] : -1;
        /* hairpin ----------------------------------------------*/

        if (SAME_STRAND(i,j)) {
          if(!no_close){
            energy = E_hp_loop(i, j, vc);
          }
        }
        else {
          /* hairpin-like exterior loop */
          if(hc_decompose & IN_EXT_LOOP){
            if (dangle_model)
              energy = E_ExtLoop(rtype[type], sj, si, P);
            else
              energy = E_ExtLoop(rtype[type], -1, -1, P);
          }
        }
        new_c = MIN2(new_c, energy);

        /*--------------------------------------------------------
          check for elementary structures involving more than one
          closing pair.
          --------------------------------------------------------*/
        for (p = i+1; p <= MIN2(j-2-TURN,i+MAXLOOP+1) ; p++) {
          int minq = j-i+p-MAXLOOP-2;
          if (minq<p+1+TURN) minq = p+1+TURN;
          for (q = minq; q < j; q++) {
            int pq = indx[q] + p;
            type_2 = ptype[pq];

            if (type_2==0) continue;
            type_2 = rtype[type_2];

            if (noGUclosure)
              if (no_close||(type_2==3)||(type_2==4))
                if ((p>i+1)||(q<j-1)) continue;  /* continue unless stack */

            if (SAME_STRAND(i,p) && SAME_STRAND(q,j)){
              if((hard_constraints[pq] & IN_INT_LOOP_ENC) && (hard_constraints[ij] & IN_INT_LOOP))
                energy = E_IntLoop(p-i-1, j-q-1, type, type_2, si, sj, S1[p-1], S1[q+1], P);
            } else {
              if((hard_constraints[pq] & IN_EXT_LOOP) && (hard_constraints[ij] & IN_EXT_LOOP))
                energy = E_IntLoop_Co(rtype[type], rtype[type_2],
                                    i, j, p, q,
                                    cut_point,
                                    si, sj,
                                    S1[p-1], S1[q+1],
                                    dangle_model,
                                    P);
            }
            new_c = MIN2(new_c, energy+my_c[indx[q]+p]);
            if ((p==i+1)&&(j==q+1)) stackEnergy = energy; /* remember stack energy */

          } /* end q-loop */
        } /* end p-loop */

        /* multi-loop decomposition ------------------------*/


        if (!no_close) {
          int MLenergy;

          if((si >= 0) && (sj >= 0)){
            if(hard_constraints[ij] & IN_MB_LOOP){
              decomp    = DMLi1[j-1];
              tt        = rtype[type];
              MLenergy  = P->MLclosing;
              switch(dangle_model){
                case 0:   MLenergy += decomp + E_MLstem(tt, -1, -1, P);
                          break;
                case 2:   MLenergy += decomp + E_MLstem(tt, sj, si, P);
                          break;
                default:  decomp += E_MLstem(tt, -1, -1, P);
                          if(hc_up_ml[j-1]){
                            energy = DMLi1[j-2] + E_MLstem(tt, sj, -1, P) + P->MLbase;
                            decomp = MIN2(decomp, energy);
                          }
                          if(hc_up_ml[i+1]){
                            energy = DMLi2[j-1] + E_MLstem(tt, -1, si, P) + P->MLbase;
                            decomp = MIN2(decomp, energy);
                          }
                          if((hc_up_ml[i+1]) && (hc_up_ml[j-1])){
                            energy = DMLi2[j-2] + E_MLstem(tt, sj, si, P) + 2*P->MLbase;
                            decomp = MIN2(decomp, energy);
                          }
                          MLenergy += decomp;
                          break;
              }
              new_c = MIN2(new_c, MLenergy);
            }
          }

          if (!SAME_STRAND(i,j)) { /* cut is somewhere in the multiloop*/
            if(hard_constraints[ij] & IN_EXT_LOOP){
              decomp = my_fc[i+1] + my_fc[j-1];
              tt = rtype[type];
              switch(dangle_model){
                case 0:   decomp += E_ExtLoop(tt, -1, -1, P);
                          break;
                case 2:   decomp += E_ExtLoop(tt, sj, si, P);
                          break;
                default:  decomp += E_ExtLoop(tt, -1, -1, P);
                          if((hc_up_ext[i+1]) && (hc_up_ext[j-1])){
                            energy = my_fc[i+2] + my_fc[j-2] + E_ExtLoop(tt, sj, si, P);
                            decomp = MIN2(decomp, energy);
                          }
                          if(hc_up_ext[i+1]){
                            energy = my_fc[i+2] + my_fc[j-1] + E_ExtLoop(tt, -1, si, P);
                            decomp = MIN2(decomp, energy);
                          }
                          if(hc_up_ext[j-1]){
                            energy = my_fc[i+1] + my_fc[j-2] + E_ExtLoop(tt, sj, -1, P);
                            decomp = MIN2(decomp, energy);
                          }
                          break;
              }
              new_c = MIN2(new_c, decomp);
            }
          }
        } /* end >> if (!no_close) << */

        /* coaxial stacking of (i.j) with (i+1.k) or (k+1.j-1) */

        if (dangle_model==3) {
          int i1k, k1j1;
          decomp = INF;
          k1j1  = indx[j-1] + i + 2 + TURN + 1;
          for (k = i+2+TURN; k < j-2-TURN; k++, k1j1++) {
            i1k = indx[k]+i+1;
            if(hard_constraints[i1k] & IN_MB_LOOP_ENC){
              type_2  = rtype[ptype[i1k]];
              energy  = my_c[i1k] + P->stack[type][type_2] + my_fML[k1j1];
              decomp  = MIN2(decomp, energy);
            }
            if(hard_constraints[k1j1] & IN_MB_LOOP_ENC){
              type_2  = rtype[ptype[k1j1]];
              energy  = my_c[k1j1] + P->stack[type][type_2] + my_fML[i1k];
              decomp  = MIN2(decomp, energy);
            }
          }
          /* no TermAU penalty if coax stack */
          decomp += 2*P->MLintern[1] + P->MLclosing;
          new_c = MIN2(new_c, decomp);
        }

        if(with_gquad){
          /* include all cases where a g-quadruplex may be enclosed by base pair (i,j) */
          if (!no_close && SAME_STRAND(i,j)) {
            tt = rtype[type];
            energy = E_GQuad_IntLoop(i, j, type, S1, my_ggg, indx, P);
            new_c = MIN2(new_c, energy);
          }
        }

        new_c = MIN2(new_c, cc1[j-1]+stackEnergy);
        cc[j] = new_c;
        if (noLP){
          if (SAME_STRAND(i,i+1) && SAME_STRAND(j-1,j))
            my_c[ij] = cc1[j-1]+stackEnergy;
          else /* currently we don't allow stacking over the cut point */
            my_c[ij] = FORBIDDEN;
        }
        else
          my_c[ij] = cc[j];

      } /* end >> if (pair) << */

      else my_c[ij] = INF;


      /* done with c[i,j], now compute fML[i,j] */
      /* free ends ? -----------------------------------------*/
      new_fML=INF;
      if (SAME_STRAND(i-1,i)) {
        if (SAME_STRAND(i,i+1))
          if(hc_up_ml[i])
            new_fML = my_fML[ij+1]+P->MLbase;
        if (SAME_STRAND(j-1,j))
          if(hc_up_ml[j]){
            new_fML = MIN2(new_fML, my_fML[indx[j-1]+i]+P->MLbase);
            if(uniq_ML)
              my_fM1[ij] = my_fM1[indx[j-1]+i] + P->MLbase;
          }
        if (SAME_STRAND(j,j+1)) {
          if(hard_constraints[ij] & IN_MB_LOOP_ENC){
            energy = my_c[ij];
            if(dangle_model == 2)
              energy += E_MLstem(type,(i>1) ? S1[i-1] : -1, (j<length) ? S1[j+1] : -1, P);
            else
              energy += E_MLstem(type, -1, -1, P);
            new_fML = MIN2(new_fML, energy);

            if(uniq_ML)
              my_fM1[ij] = MIN2(my_fM1[ij], energy);
          } else if (with_gquad) {
              int gggg = my_ggg[ij] + E_MLstem(0, -1, -1, P);
              energy = MIN2(energy, gggg);
              new_fML = MIN2(new_fML, energy);
              if(uniq_ML)
                my_fM1[ij] = MIN2(my_fM1[ij], energy);
          }
        
        }
        if (dangle_model%2==1) {  /* normal dangles */
          if (SAME_STRAND(i,i+1)) {
            if((hard_constraints[ij+1] & IN_MB_LOOP_ENC) && (hc_up_ml[i])){
              tt      = ptype[ij+1]; /* i+1,j */
              energy  = my_c[ij+1] + P->MLbase + E_MLstem(tt, S1[i], -1, P);
              new_fML = MIN2(new_fML, energy);
            }
          }
          if (SAME_STRAND(j-1,j)) {
            if((hard_constraints[indx[j-1]+i] & IN_MB_LOOP_ENC) && (hc_up_ml[j])){
              tt      = ptype[indx[j-1]+i]; /* i,j-1 */
              energy  = my_c[indx[j-1]+i] + P->MLbase + E_MLstem(tt, -1, S1[j], P);
              new_fML = MIN2(new_fML, energy);
            }
          }
          if ((SAME_STRAND(j-1,j))&&(SAME_STRAND(i,i+1))) {
            if((hard_constraints[indx[j-1]+i+1] & IN_MB_LOOP_ENC) && (hc_up_ml[i]) && (hc_up_ml[j])){
              tt      = ptype[indx[j-1]+i+1]; /* i+1,j-1 */
              energy  = my_c[indx[j-1]+i+1] + 2*P->MLbase + E_MLstem(tt, S1[i], S1[j], P);
              new_fML = MIN2(new_fML, energy);
            }
          }
        }
      }

      if(with_gquad){
        if(SAME_STRAND(i, j))
          new_fML = MIN2(new_fML, my_ggg[indx[j] + i] + E_MLstem(0, -1, -1, P));
      }

      /* modular decomposition -------------------------------*/

      {
        int stopp;     /*loop 1 up to cut, then loop 2*/
        stopp=(cut_point>0)? (cut_point):(j-2-TURN);
        for (decomp=INF, k = i+1+TURN; k<stopp; k++)
          decomp = MIN2(decomp, Fmi[k]+my_fML[indx[j]+k+1]);
        k++;
        for (;k <= j-2-TURN;k++)
          decomp = MIN2(decomp, Fmi[k]+my_fML[indx[j]+k+1]);
      }
      DMLi[j] = decomp;               /* store for use in ML decompositon */
      new_fML = MIN2(new_fML,decomp);

      /* coaxial stacking */
      if (dangle_model==3) {
        int stopp;
        stopp=(cut_point>0)? (cut_point):(j-2-TURN);
        /* additional ML decomposition as two coaxially stacked helices */
        for (decomp = INF, k = i+1+TURN; k<stopp; k++) {
          type = ptype[indx[k]+i]; type = rtype[type];
          type_2 = ptype[indx[j]+k+1]; type_2 = rtype[type_2];
          if (type && type_2)
            decomp = MIN2(decomp,
                          my_c[indx[k]+i]+my_c[indx[j]+k+1]+P->stack[type][type_2]);
        }
        k++;
        for (;k <= j-2-TURN; k++) {
          type = ptype[indx[k]+i]; type = rtype[type];
          type_2 = ptype[indx[j]+k+1]; type_2 = rtype[type_2];
          if (type && type_2)
            decomp = MIN2(decomp,
                          my_c[indx[k]+i]+my_c[indx[j]+k+1]+P->stack[type][type_2]);
        }

        decomp += 2*P->MLintern[1];

#if 0
        /* This is needed for Y shaped ML loops with coax stacking of
           interior pairs, but backtracking will fail if activated */
        DMLi[j] = MIN2(DMLi[j], decomp);
        if (SAME_STRAND(j-1,j)) DMLi[j] = MIN2(DMLi[j], DMLi[j-1]+P->MLbase);
        if (SAME_STRAND(i,i+1)) DMLi[j] = MIN2(DMLi[j], DMLi1[j]+P->MLbase);
        new_fML = MIN2(new_fML, DMLi[j]);
#endif
        new_fML = MIN2(new_fML, decomp);
      }

      my_fML[ij] = Fmi[j] = new_fML;     /* substring energy */

    }

    if (i==cut_point)
      for (j=i; j<=maxj; j++)
        free_end(my_fc, j, cut_point, vc);
    if (i<cut_point)
      free_end(my_fc,i,cut_point-1, vc);


    {
      int *FF; /* rotate the auxilliary arrays */
      FF = DMLi2; DMLi2 = DMLi1; DMLi1 = DMLi; DMLi = FF;
      FF = cc1; cc1=cc; cc=FF;
      for (j=1; j<=maxj; j++) {cc[j]=Fmi[j]=DMLi[j]=INF; }
    }
  }

  /* calculate energies of 5' and 3' fragments */

  for (i=1; i<=length; i++)
    free_end(my_f5, i, 1, vc);

  if (cut_point>0) {
    mfe1=my_f5[cut_point-1];
    mfe2=my_fc[length];
    /* add DuplexInit, check whether duplex*/
    for (i=cut_point; i<=length; i++) {
      my_f5[i]=MIN2(my_f5[i]+P->DuplexInit, my_fc[i]+my_fc[1]);
    }
  }

  energy = my_f5[length];
  if (cut_point<1) mfe1=mfe2=energy;
  return energy;
}

PRIVATE void backtrack_co(sect bt_stack[],
                          bondT *bp_list,
                          int s,
                          int b, /* b=0: start new structure, b \ne 0: add to existing structure */
                          vrna_fold_compound *vc) {

  /*------------------------------------------------------------------
    trace back through the "c", "fc", "f5" and "fML" arrays to get the
    base pairing list. No search for equivalent structures is done.
    This is fast, since only few structure elements are recalculated.
    ------------------------------------------------------------------*/

  int   i, j, k, length, energy, new, ml0, ml5, ml3, ml53, no_close, type, type_2, tt;
  soft_constraintT      *sc;
  hard_constraintT      *hc;
  char  *string         = vc->sequence;
  paramT  *P            = vc->params;
  int     *indx         = vc->jindx;
  char    *ptype        = vc->ptype;

  short *S1             = vc->sequence_encoding;
  short *S              = vc->sequence_encoding2;
  int   dangle_model    = P->model_details.dangles;
  int   noLP            = P->model_details.noLP;
  int   noGUclosure     = P->model_details.noGUclosure;
  int   with_gquad      = P->model_details.gquad;
  int   *rtype          = &(P->model_details.rtype[0]);
  char  backtrack_type  = P->model_details.backtrack_type;
  int   cut_point       = vc->cutpoint;

  /* the folding matrices */
  int   *my_f5, *my_c, *my_fML, *my_fc, *my_ggg;

  length  = vc->length;
  my_f5   = vc->matrices->f5;
  my_c    = vc->matrices->c;
  my_fML  = vc->matrices->fML;
  my_fc   = vc->matrices->fc;
  my_ggg  = vc->matrices->ggg;
  hc      = vc->hc;
  sc      = vc->sc;

  /* int   b=0;*/

  length = strlen(string);
  if (s==0) {
    bt_stack[++s].i = 1;
    bt_stack[s].j   = length;
    bt_stack[s].ml  = (backtrack_type=='M') ? 1 : ((backtrack_type=='C')?2:0);
  }
  while (s>0) {
    int ml, fij, fi, cij, traced, i1, j1, mm, p, q, jj=0, gq=0;
    int canonical = 1;     /* (i,j) closes a canonical structure */
    i  = bt_stack[s].i;
    j  = bt_stack[s].j;
    ml = bt_stack[s--].ml;   /* ml is a flag indicating if backtracking is to
                              occur in the fML- (1) or in the f-array (0) */
    if (ml==2) {
      bp_list[++b].i = i;
      bp_list[b].j   = j;
      goto repeat1;
    }

    if (j < i+TURN+1) continue; /* no more pairs in this interval */


    if (ml==0) {fij = my_f5[j]; fi = my_f5[j-1];}
    else if (ml==1) {fij = my_fML[indx[j]+i]; fi = my_fML[indx[j-1]+i]+P->MLbase;}
    else /* 3 or 4 */ {
      fij = my_fc[j];
      fi = (ml==3) ? INF : my_fc[j-1];
    }
    if (fij == fi) {  /* 3' end is unpaired */
      bt_stack[++s].i = i;
      bt_stack[s].j   = j-1;
      bt_stack[s].ml  = ml;
      continue;
    }

    if (ml==0 || ml==4) { /* backtrack in f5 or fc[i=cut,j>cut] */
      int *ff;
      ff = (ml==4) ? my_fc : my_f5;
      switch(dangle_model){
        case 0:   /* j or j-1 is paired. Find pairing partner */
                  for (k=j-TURN-1,traced=0; k>=i; k--) {
                    int cc;

                    if(with_gquad){
                      if(fij == ff[k-1] + my_ggg[indx[j]+k]){
                        /* found the decomposition */
                        traced = j; jj = k - 1; gq = 1;
                        break;
                      }
                    }

                    type = ptype[indx[j]+k];
                    if(type){
                      cc = my_c[indx[j]+k];
                      if(!SAME_STRAND(k,j)) cc += P->DuplexInit;
                      if(fij == ff[k-1] + cc + E_ExtLoop(type, -1, -1, P)){
                        traced = j; jj = k-1;
                      }
                    }
                    if(traced) break;
                  }

                  break;

        case 2:   /* j or j-1 is paired. Find pairing partner */
                  for (k=j-TURN-1,traced=0; k>=i; k--) {
                    int cc;

                    if(with_gquad){
                      if(fij == ff[k-1] + my_ggg[indx[j]+k]){
                        /* found the decomposition */
                        traced = j; jj = k - 1; gq = 1;
                        break;
                      }
                    }

                    type = ptype[indx[j]+k];
                    if(type){
                      cc = my_c[indx[j]+k];
                      if(!SAME_STRAND(k,j)) cc += P->DuplexInit;
                      if(fij == ff[k-1] + cc + E_ExtLoop(type, (k>1) && SAME_STRAND(k-1,k) ? S1[k-1] : -1, (j<length) && SAME_STRAND(j,j+1) ? S1[j+1] : -1, P)){
                        traced = j; jj = k-1;
                      }
                    }
                    if(traced) break;
                  }
                  break;

        default:  for(k=j-TURN-1,traced=0; k>=i; k--){
                    int cc;
                    type = ptype[indx[j]+k];

                    if(with_gquad){
                      if(fij == ff[k-1] + my_ggg[indx[j]+k]){
                        /* found the decomposition */
                        traced = j; jj = k - 1; gq = 1;
                        break;
                      }
                    }

                    if(type){
                      cc = my_c[indx[j]+k];
                      if(!SAME_STRAND(k,j)) cc += P->DuplexInit;
                      if(fij == ff[k-1] + cc + E_ExtLoop(type, -1, -1, P)){
                        traced = j; jj = k-1; break;
                      }
                      if((k>1) && SAME_STRAND(k-1,k))
                        if(fij == ff[k-2] + cc + E_ExtLoop(type, S1[k-1], -1, P)){
                              traced=j; jj=k-2; break;
                        }
                    }

                    type = ptype[indx[j-1]+k];
                    if(type && SAME_STRAND(j-1,j)){
                      cc = my_c[indx[j-1]+k];
                      if (!SAME_STRAND(k,j-1)) cc += P->DuplexInit; /*???*/
                      if (fij == cc + ff[k-1] + E_ExtLoop(type, -1, S1[j], P)){
                            traced=j-1; jj = k-1; break;
                      }
                      if(k>i){
                        if (fij == ff[k-2] + cc + E_ExtLoop(type, SAME_STRAND(k-1,k) ? S1[k-1] : -1, S1[j], P)){
                          traced=j-1; jj=k-2; break;
                        }
                      }
                    }
                  }

                  break;
      }
      if (!traced) nrerror("backtrack failed in f5 (or fc)");
      bt_stack[++s].i = i;
      bt_stack[s].j   = jj;
      bt_stack[s].ml  = ml;

      i=k; j=traced;

      if(with_gquad && gq){
        /* goto backtrace of gquadruplex */
        goto repeat_gquad;
      }

      bp_list[++b].i = i;
      bp_list[b].j   = j;

      goto repeat1;
    }
    else if (ml==3) { /* backtrack in fc[i<cut,j=cut-1] */
      if (my_fc[i] == my_fc[i+1]) { /* 5' end is unpaired */
        bt_stack[++s].i = i+1;
        bt_stack[s].j   = j;
        bt_stack[s].ml  = ml;
        continue;
      }
      /* i or i+1 is paired. Find pairing partner */
      switch(dangle_model){
        case 0:   for (k=i+TURN+1, traced=0; k<=j; k++){
                    jj=k+1;
                    type = ptype[indx[k]+i];
                    if (type) {
                      if(my_fc[i] == my_fc[k+1] + my_c[indx[k]+i] + E_ExtLoop(type, -1, -1, P)){
                        traced = i;
                      }
                    } else if (with_gquad){
                      if(my_fc[i] == my_fc[k+1] + my_ggg[indx[k]+i]){
                        traced = i; gq = 1;
                        break;
                      }
                    }

                    if (traced) break;
                  }
                  break;
        case 2:   for (k=i+TURN+1, traced=0; k<=j; k++){
                    jj=k+1;
                    type = ptype[indx[k]+i];
                    if(type){
                      if(my_fc[i] == my_fc[k+1] + my_c[indx[k]+i] + E_ExtLoop(type,(i>1 && SAME_STRAND(i-1,i)) ? S1[i-1] : -1,  SAME_STRAND(k,k+1) ? S1[k+1] : -1, P)){
                        traced = i;
                      }
                    } else if (with_gquad){
                      if(my_fc[i] == my_fc[k+1] + my_ggg[indx[k]+i]){
                        traced = i; gq = 1;
                        break;
                      }
                    }
                    if (traced) break;
                  }
                  break;
        default:  for(k=i+TURN+1, traced=0; k<=j; k++){
                    jj=k+1;
                    type = ptype[indx[k]+i];
                    if(type){
                      if(my_fc[i] == my_fc[k+1] + my_c[indx[k]+i] + E_ExtLoop(type, -1, -1, P)){
                        traced = i; break;
                      }
                      else if(my_fc[i] == my_fc[k+2] + my_c[indx[k]+i] + E_ExtLoop(type, -1, SAME_STRAND(k,k+1) ? S1[k+1] : -1, P)){
                        traced = i; jj=k+2; break;
                      }
                    } else if (with_gquad){
                      if(my_fc[i] == my_fc[k+1] + my_ggg[indx[k]+i]){
                        traced = i; gq = 1;
                        break;
                      }
                    }

                    type = ptype[indx[k]+i+1];
                    if(type){
                      if(my_fc[i] == my_fc[k+1] + my_c[indx[k]+i+1] + E_ExtLoop(type, SAME_STRAND(i, i+1) ? S1[i] : -1, -1, P)){
                        traced = i+1; break;
                      }
                      if(k<j){
                        if(my_fc[i] == my_fc[k+2] + my_c[indx[k]+i+1] + E_ExtLoop(type, SAME_STRAND(i, i+1) ? S1[i] : -1, SAME_STRAND(k, k+1) ? S1[k+1] : -1, P)){
                          traced = i+1; jj=k+2; break;
                        }
                      }
                    }
                  }
                  break;
      }

      if (!traced) nrerror("backtrack failed in fc[] 5' of cut");

      bt_stack[++s].i = jj;
      bt_stack[s].j   = j;
      bt_stack[s].ml  = ml;

      j=k; i=traced;

      if(with_gquad && gq){
        /* goto backtrace of gquadruplex */
        goto repeat_gquad;
      }

      bp_list[++b].i = i;
      bp_list[b].j   = j;
      goto repeat1;
    }

    else { /* true multi-loop backtrack in fML */
      if (my_fML[indx[j]+i+1]+P->MLbase == fij) { /* 5' end is unpaired */
        bt_stack[++s].i = i+1;
        bt_stack[s].j   = j;
        bt_stack[s].ml  = ml;
        continue;
      }

      if(with_gquad){
        if(fij == my_ggg[indx[j]+i] + E_MLstem(0, -1, -1, P)){
          /* go to backtracing of quadruplex */
          goto repeat_gquad;
        }
      }

      tt  = ptype[indx[j]+i];
      cij = my_c[indx[j]+i];
      switch(dangle_model){
        case 0:   if(fij == cij + E_MLstem(tt, -1, -1, P)){
                    bp_list[++b].i  = i;
                    bp_list[b].j    = j;
                    goto repeat1;
                  }
                  break;
        case 2:   if(fij == cij + E_MLstem(tt, (i>1) ? S1[i-1] : -1, (j<length) ? S1[j+1] : -1, P)){
                    bp_list[++b].i  = i;
                    bp_list[b].j    = j;
                    goto repeat1;
                  }
                  break;
        default:  if(fij == cij + E_MLstem(tt, -1, -1, P)){
                    bp_list[++b].i  = i;
                    bp_list[b].j    = j;
                    goto repeat1;
                  }
                  tt = ptype[indx[j]+i+1];
                  if(fij == my_c[indx[j]+i+1] + P->MLbase + E_MLstem(tt, S1[i], -1, P)){
                    i++;
                    bp_list[++b].i  = i;
                    bp_list[b].j    = j;
                    goto repeat1;
                  }
                  tt = ptype[indx[j-1]+i];
                  if(fij == my_c[indx[j-1]+i] + P->MLbase + E_MLstem(tt, -1, S1[j], P)){
                    j--;
                    bp_list[++b].i  = i;
                    bp_list[b].j    = j;
                    goto repeat1;
                  }
                  tt = ptype[indx[j-1]+i+1];
                  if(fij == my_c[indx[j-1]+i+1] + 2*P->MLbase + E_MLstem(tt, S1[i], S1[j], P)){
                    i++; j--;
                    bp_list[++b].i  = i;
                    bp_list[b].j    = j;
                    goto repeat1;
                  }
                  break;
      }

      /* find next component of multiloop */
      for (k = i+1+TURN; k <= j-2-TURN; k++)
        if (fij == (my_fML[indx[k]+i]+my_fML[indx[j]+k+1]))
          break;

      if ((dangle_model==3)&&(k>j-2-TURN)) { /* must be coax stack */
        ml = 2;
        for (k = i+1+TURN; k <= j-2-TURN; k++) {
          type = ptype[indx[k]+i];  type= rtype[type];
          type_2 = ptype[indx[j]+k+1]; type_2= rtype[type_2];
          if (type && type_2)
            if (fij == my_c[indx[k]+i]+my_c[indx[j]+k+1]+P->stack[type][type_2]+
                2*P->MLintern[1])
              break;
        }
      }

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
    if (canonical)  cij = my_c[indx[j]+i];
    type = ptype[indx[j]+i];

    if (noLP)
      if (cij == my_c[indx[j]+i]) {
        /* (i.j) closes canonical structures, thus
           (i+1.j-1) must be a pair                */
        type_2 = ptype[indx[j-1]+i+1]; type_2 = rtype[type_2];
        cij -= P->stack[type][type_2];
        bp_list[++b].i = i+1;
        bp_list[b].j   = j-1;
        i++; j--;
        canonical=0;
        goto repeat1;
      }
    canonical = 1;


    no_close = (((type==3)||(type==4))&&noGUclosure);
    if (SAME_STRAND(i,j)) {
      if (no_close) {
        if (cij == FORBIDDEN) continue;
      } else
        if (cij == E_Hairpin(j-i-1, type, S1[i+1], S1[j-1],string+i-1, P))
          continue;
    }
    else {
      if(dangle_model){
        if(cij == E_ExtLoop(rtype[type], SAME_STRAND(j-1,j) ? S1[j-1] : -1, SAME_STRAND(i,i+1) ? S1[i+1] : -1, P)) continue;
      }
      else if(cij == E_ExtLoop(rtype[type], -1, -1, P)) continue;
    }

    for (p = i+1; p <= MIN2(j-2-TURN,i+MAXLOOP+1); p++) {
      int minq;
      minq = j-i+p-MAXLOOP-2;
      if (minq<p+1+TURN) minq = p+1+TURN;
      for (q = j-1; q >= minq; q--) {

        type_2 = ptype[indx[q]+p];
        if (type_2==0) continue;
        type_2 = rtype[type_2];
        if (noGUclosure)
          if (no_close||(type_2==3)||(type_2==4))
            if ((p>i+1)||(q<j-1)) continue;  /* continue unless stack */

        /* energy = oldLoopEnergy(i, j, p, q, type, type_2); */
        if (SAME_STRAND(i,p) && SAME_STRAND(q,j))
          energy = E_IntLoop(p-i-1, j-q-1, type, type_2,
                              S1[i+1], S1[j-1], S1[p-1], S1[q+1], P);
        else {
          energy = E_IntLoop_Co(rtype[type], rtype[type_2], i, j, p, q, cut_point, S1[i+1], S1[j-1], S1[p-1], S1[q+1], dangle_model, P);
        }

        new = energy+my_c[indx[q]+p];
        traced = (cij == new);
        if (traced) {
          bp_list[++b].i = p;
          bp_list[b].j   = q;
          i = p, j = q;
          goto repeat1;
        }
      }
    }

    /* end of repeat: --------------------------------------------------*/

    /* (i.j) must close a fake or true multi-loop */
    tt = rtype[type];
    i1 = i+1;
    j1 = j-1;

    if(with_gquad){
      /*
        The case that is handled here actually resembles something like
        an interior loop where the enclosing base pair is of regular
        kind and the enclosed pair is not a canonical one but a g-quadruplex
        that should then be decomposed further...
      */
      if(SAME_STRAND(i,j)){
        if(backtrack_GQuad_IntLoop(cij, i, j, type, S, my_ggg, indx, &p, &q, P)){
          i = p; j = q;
          goto repeat_gquad;
        }
      }
    }

    /* fake multi-loop */
    if(!SAME_STRAND(i,j)){
      int ii, jj, decomp;
      ii = jj = 0;
      decomp = my_fc[i1] + my_fc[j1];
      switch(dangle_model){
        case 0:   if(cij == decomp + E_ExtLoop(tt, -1, -1, P)){
                    ii=i1, jj=j1;
                  }
                  break;
        case 2:   if(cij == decomp + E_ExtLoop(tt, SAME_STRAND(j-1,j) ? S1[j-1] : -1, SAME_STRAND(i,i+1) ? S1[i+1] : -1, P)){
                    ii=i1, jj=j1;
                  }

                  break;
        default:  if(cij == decomp + E_ExtLoop(tt, -1, -1, P)){
                    ii=i1, jj=j1;
                  }
                  else if(cij == my_fc[i+2] + my_fc[j-1] + E_ExtLoop(tt, -1, SAME_STRAND(i,i+1) ? S1[i+1] : -1, P)){
                    ii = i+2; jj = j1;
                  }
                  else if(cij == my_fc[i+1] + my_fc[j-2] + E_ExtLoop(tt, SAME_STRAND(j-1,j) ? S1[j-1] : -1, -1, P)){
                    ii = i1; jj = j-2;
                  }
                  else if(cij == my_fc[i+2] + my_fc[j-2] + E_ExtLoop(tt, SAME_STRAND(j-1,j) ? S1[j-1] : -1, SAME_STRAND(i,i+1) ? S1[i+1] : -1, P)){
                    ii = i+2; jj = j-2;
                  }
                  break;
      }
      if(ii){
        bt_stack[++s].i = ii;
        bt_stack[s].j   = cut_point-1;
        bt_stack[s].ml  = 3;
        bt_stack[++s].i = cut_point;
        bt_stack[s].j   = jj;
        bt_stack[s].ml  = 4;
        continue;
      }
    }

    /* true multi-loop */
    mm = P->MLclosing;
    bt_stack[s+1].ml  = bt_stack[s+2].ml = 1;
    ml0   = E_MLstem(tt, -1, -1, P);
    ml5   = E_MLstem(tt, SAME_STRAND(j-1,j) ? S1[j-1] : -1, -1, P);
    ml3   = E_MLstem(tt, -1, SAME_STRAND(i,i+1) ? S1[i+1] : -1, P);
    ml53  = E_MLstem(tt, SAME_STRAND(j-1,j) ? S1[j-1] : -1, SAME_STRAND(i,i+1) ? S1[i+1] : -1, P);
    for (traced = 0, k = i+2+TURN; k < j-2-TURN; k++) {
      switch(dangle_model){
        case 0:   /* no dangles */
                  if(cij == mm + my_fML[indx[k]+i+1] + my_fML[indx[j-1]+k+1] + ml0)
                    traced = i+1;
                  break;
        case 2:   /*double dangles */
                  if(cij == mm + my_fML[indx[k]+i+1] + my_fML[indx[j-1]+k+1] + ml53)
                    traced = i+1;
                  break;
        default:  /* normal dangles */
                  if(cij == mm + my_fML[indx[k]+i+1] + my_fML[indx[j-1]+k+1] + ml0){
                    traced = i+1;
                    break;
                  }
                  else if (cij == my_fML[indx[k]+i+2] + my_fML[indx[j-1]+k+1] + ml3 + mm + P->MLbase){
                    traced = i1 = i+2;
                    break;
                  }
                  else if (cij == my_fML[indx[k]+i+1] + my_fML[indx[j-2]+k+1] + ml5 + mm + P->MLbase){
                    traced = i1 = i+1; j1 = j-2;
                    break;
                  }
                  else if (cij == my_fML[indx[k]+i+2] + my_fML[indx[j-2]+k+1] + ml53 + mm + 2*P->MLbase){
                    traced = i1 = i+2; j1 = j-2;
                    break;
                  }
                  break;
      }
      if(traced) break;
      /* coaxial stacking of (i.j) with (i+1.k) or (k.j-1) */
      /* use MLintern[1] since coax stacked pairs don't get TerminalAU */
      if (dangle_model==3) {
        int en;
        type_2 = ptype[indx[k]+i+1]; type_2 = rtype[type_2];
        if (type_2) {
          en = my_c[indx[k]+i+1]+P->stack[type][type_2]+my_fML[indx[j-1]+k+1];
          if (cij == en+2*P->MLintern[1]+P->MLclosing) {
            ml = 2;
            bt_stack[s+1].ml  = 2;
            break;
          }
        }
        type_2 = ptype[indx[j-1]+k+1]; type_2 = rtype[type_2];
        if (type_2) {
          en = my_c[indx[j-1]+k+1]+P->stack[type][type_2]+my_fML[indx[k]+i+1];
          if (cij == en+2*P->MLintern[1]+P->MLclosing) {
            bt_stack[s+2].ml = 2;
            break;
          }
        }
      }
    }
    if (k<=j-3-TURN) { /* found the decomposition */
      bt_stack[++s].i = i1;
      bt_stack[s].j   = k;
      bt_stack[++s].i = k+1;
      bt_stack[s].j   = j1;
    } else {
#if 0
      /* Y shaped ML loops don't work yet */
      if (dangle_model==3) {
        /* (i,j) must close a Y shaped ML loop with coax stacking */
        if (cij == fML[indx[j-2]+i+2] + mm + d3 + d5 + P->MLbase + P->MLbase) {
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
        nrerror("backtracking failed in repeat");
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
          bp_list[++b].i = i+a;
          bp_list[b].j   = i+a;
          bp_list[++b].i = i+L+l[0]+a;
          bp_list[b].j   = i+L+l[0]+a;
          bp_list[++b].i = i+L+l[0]+L+l[1]+a;
          bp_list[b].j   = i+L+l[0]+L+l[1]+a;
          bp_list[++b].i = i+L+l[0]+L+l[1]+L+l[2]+a;
          bp_list[b].j   = i+L+l[0]+L+l[1]+L+l[2]+a;
        }
        goto repeat_gquad_exit;
      }
      nrerror("backtracking failed in repeat_gquad");
    }
  repeat_gquad_exit:
    asm("nop");

  } /* end >> while (s>0) << */

  bp_list[0].i = b;    /* save the total number of base pairs */
}

PRIVATE void
free_end( int *array,
          int i,
          int start,
          vrna_fold_compound *vc){

  int inc, type, energy, length, j, left, right, cut_point, dangle_model, with_gquad, *indx, *c, *ggg;
  paramT        *P;
  short         *S1;
  char          *ptype;
  mfe_matricesT *matrices;

  cut_point     = vc->cutpoint;
  P             = vc->params;
  dangle_model  = P->model_details.dangles;
  with_gquad    = P->model_details.gquad;
  inc           = (i>start)? 1:-1;
  length        = (int)vc->length;
  S1            = vc->sequence_encoding;
  ptype         = vc->ptype;
  indx          = vc->jindx;
  matrices      = vc->matrices;
  c             = matrices->c;
  ggg           = matrices->ggg;

  if (i==start) array[i]=0;
  else array[i] = array[i-inc];
  if (inc>0) {
    left = start; right=i;
  } else {
    left = i; right = start;
  }

  for (j=start; inc*(i-j)>TURN; j+=inc) {
    int ii, jj;
    short si, sj;
    if (i>j) { ii = j; jj = i;} /* inc>0 */
    else     { ii = i; jj = j;} /* inc<0 */
    type = ptype[indx[jj]+ii];
    if (type) {  /* i is paired with j */
      si = (ii>1)       && SAME_STRAND(ii-1,ii) ? S1[ii-1] : -1;
      sj = (jj<length)  && SAME_STRAND(jj,jj+1) ? S1[jj+1] : -1;
      energy = c[indx[jj]+ii];
      switch(dangle_model){
        case 0:   
                  array[i] = MIN2(array[i], array[j-inc] + energy + E_ExtLoop(type, -1, -1, P));
                  break;
        case 2:   
                  array[i] = MIN2(array[i], array[j-inc] + energy + E_ExtLoop(type, si, sj, P));
                  break;
        default:     
                  array[i] = MIN2(array[i], array[j-inc] + energy + E_ExtLoop(type, -1, -1, P));
                  if(inc > 0){
                    if(j > left)
                    array[i] = MIN2(array[i], array[j-2] + energy + E_ExtLoop(type, si, -1, P));
                  }
                  else if(j < right)
                    array[i] = MIN2(array[i], array[j+2] + energy + E_ExtLoop(type, -1, sj, P));
                  break;
      }
    }

    if(with_gquad){
      if(SAME_STRAND(ii, jj))
        array[i] = MIN2(array[i], array[j-inc] + ggg[indx[jj]+ii]);
    }

    if (dangle_model%2==1) {
      /* interval ends in a dangle (i.e. i-inc is paired) */
      if (i>j) { ii = j; jj = i-1;} /* inc>0 */
      else     { ii = i+1; jj = j;} /* inc<0 */
      type = ptype[indx[jj]+ii];
      if (!type) continue;

      si = (ii > left)  && SAME_STRAND(ii-1,ii) ? S1[ii-1] : -1;
      sj = (jj < right) && SAME_STRAND(jj,jj+1) ? S1[jj+1] : -1;
      energy = c[indx[jj]+ii];
      if(inc>0)
        array[i] = MIN2(array[i], array[j - inc] + energy + E_ExtLoop(type, -1, sj, P));
      else
        array[i] = MIN2(array[i], array[j - inc] + energy + E_ExtLoop(type, si, -1, P));
      if(j!= start){ /* dangle_model on both sides */
        array[i] = MIN2(array[i], array[j-2*inc] + energy + E_ExtLoop(type, si, sj, P));
      }
    }
  }
}

PUBLIC void
update_cofold_params(void){

  vrna_update_fold_params(backward_compat_compound, NULL);
}

PUBLIC void
update_cofold_params_par(paramT *parameters){

  vrna_update_fold_params(backward_compat_compound, parameters);
}

PUBLIC void get_monomere_mfes(float *e1, float *e2) {
  /*exports monomere free energies*/
  *e1 = mfe1;
  *e2 = mfe2;
}

PRIVATE void
backtrack(sect bt_stack[],
          bondT *bp_list,
          vrna_fold_compound *vc){

  /*routine to call backtrack_co from 1 to n, backtrack type??*/
  backtrack_co(bt_stack, bp_list, 0,0, vc);
}

#if 0
PRIVATE int comp_pair(const void *A, const void *B) {
  bondT *x,*y;
  int ex, ey;
  x = (bondT *) A;
  y = (bondT *) B;
  ex = c[indx[x->j]+x->i]+c[indx[x->i+length]+x->j];
  ey = c[indx[y->j]+y->i]+c[indx[y->i+length]+y->j];
  if (ex>ey) return 1;
  if (ex<ey) return -1;
  return (indx[x->j]+x->i - indx[y->j]+y->i);
}
#endif

PRIVATE SOLUTION *
wrap_zukersubopt( const char *string,
                  paramT *parameters){

  unsigned int        length;
  char                *doubleseq;
  vrna_fold_compound  *vc;
  hard_constraintT    *my_hc;
  soft_constraintT    *my_sc;
  paramT              *P;

  my_hc   = NULL;
  my_sc   = NULL;
  vc      = NULL;
  length  = (int)strlen(string);

#ifdef _OPENMP
/* Explicitly turn off dynamic threads */
  omp_set_dynamic(0);
#endif

  /* we need the parameter structure for hard constraints */
  if(parameters)
    P = get_parameter_copy(parameters);
  else{
    model_detailsT md;
    set_model_details(&md);
    P = get_scaled_parameters(temperature, md);
  }
  P->model_details.min_loop_size = TURN;  /* set min loop length to 0 */

  doubleseq = (char *)space((2*length+2)*sizeof(char));
  strcpy(doubleseq,string);
  doubleseq[length] = '&';
  strcat(doubleseq, string);

  /* get compound structure */
  vc = vrna_get_fold_compound(doubleseq, &(P->model_details), VRNA_OPTION_MFE | VRNA_OPTION_HYBRID);

  if(parameters){ /* replace params if necessary */
    free(vc->params);
    vc->params = P;
  } else {
    free(P);
  }

  if(backward_compat_compound)
    destroy_fold_compound(backward_compat_compound);

  backward_compat_compound = vc;

  /* cleanup */
  free(doubleseq);

  return vrna_zukersubopt(vc);
}

PUBLIC SOLUTION *
vrna_zukersubopt(vrna_fold_compound *vc){

/* Compute zuker suboptimal. Here, we're abusing the cofold() code
   "double" sequence, compute dimerarray entries, track back every base pair.
   This is slightly wasteful compared to the normal solution */

  char          *structure, *mfestructure, **todo, *ptype;
  int           i, j, counter, num_pairs, psize, p, *indx, *c;
  unsigned int  length, doublelength;
  float         energy;
  SOLUTION      *zukresults;
  bondT         *pairlist, *bp_list;
  sect          bt_stack[MAXSECTORS]; /* stack of partial structures for backtracking */
  mfe_matricesT *matrices;

  doublelength    = vc->length;
  length          = doublelength/2;
  indx            = vc->jindx;
  ptype           = vc->ptype;
  matrices        = vc->matrices;
  c               = matrices->c;
  num_pairs       = counter = 0;
  mfestructure    = (char *) space((unsigned) doublelength+1);
  structure       = (char *) space((unsigned) doublelength+1);
  zukresults      = (SOLUTION *)space(((length*(length-1))/2)*sizeof(SOLUTION));
  mfestructure[0] = '\0';

  /* store length at pos. 0 */
  vc->sequence_encoding[0] = vc->sequence_encoding2[0];

  /* get mfe and do forward recursion */
  (void)fill_arrays(vc, 1);

  psize     = length;
  pairlist  = (bondT *) space(sizeof(bondT)*(psize+1));
  bp_list   = (bondT *)space(sizeof(bondT) * (1 + length/2));
  todo      = (char **) space(sizeof(char *)*(length+1));
  for (i=1; i<length; i++) {
    todo[i] = (char *) space(sizeof(char)*(length+1));
  }

  /* Make a list of all base pairs */
  for (i=1; i<length; i++) {
    for (j=i+TURN2+1/*??*/; j<=length; j++) {
      if (ptype[indx[j]+i]==0) continue;
      if (num_pairs>=psize) {
        psize = 1.2*psize + 32;
        pairlist = xrealloc(pairlist, sizeof(bondT)*(psize+1));
      }
      pairlist[num_pairs].i   = i;
      pairlist[num_pairs++].j = j;
      todo[i][j]=1;
    }
  }

#if 0
  qsort(pairlist, num_pairs, sizeof(bondT), comp_pair);
#endif

  for (p=0; p<num_pairs; p++) {
    i=pairlist[p].i;
    j=pairlist[p].j;
    if (todo[i][j]) {
      int k;
      bt_stack[1].i   = i;
      bt_stack[1].j   = j;
      bt_stack[1].ml  = 2;
      backtrack_co(bt_stack, bp_list, 1,0, vc);
      bt_stack[1].i   = j;
      bt_stack[1].j   = i + length;
      bt_stack[1].ml  = 2;
      backtrack_co(bt_stack, bp_list, 1,bp_list[0].i, vc);
      energy = c[indx[j]+i]+c[indx[i+length]+j];
      vrna_parenthesis_zuker(structure, bp_list, length);
      zukresults[counter].energy      = energy;
      zukresults[counter++].structure = strdup(structure);
      for (k = 1; k <= bp_list[0].i; k++) { /* mark all pairs in structure as done */
        int x,y;
        x=bp_list[k].i;
        y=bp_list[k].j;
        if (x>length) x-=length;
        if (y>length) y-=length;
        if (x>y) {
          int temp;
          temp=x; x=y; y=temp;
        }
        todo[x][y] = 0;
      }
    }
  }

  /*
  *  Backward compatibility:
  *  This block may be removed if deprecated functions
  *  relying on the global variable "base_pair" vanish from within the package!
  */
  if(base_pair) free(base_pair);
  base_pair = bp_list;

  /* clean up */
  free(pairlist);
  for (i=1; i<length; i++)
    free(todo[i]);
  free(todo);
  free(structure);
  free(mfestructure);

  return zukresults;
}


/*###########################################*/
/*# deprecated functions below              #*/
/*###########################################*/

PUBLIC void
initialize_cofold(int length){ /* DO NOTHING */ }

PUBLIC void
free_co_arrays(void){

  if(backward_compat_compound){
    destroy_fold_compound(backward_compat_compound);
    backward_compat_compound = NULL;
  }
}


/*--------------------------------------------------------------------------*/

PUBLIC void
export_cofold_arrays_gq(int **f5_p,
                        int **c_p,
                        int **fML_p,
                        int **fM1_p,
                        int **fc_p,
                        int **ggg_p,
                        int **indx_p,
                        char **ptype_p){

  /* make the DP arrays available to routines such as subopt() */
  wrap_array_export(f5_p, c_p, fML_p, fM1_p, fc_p, indx_p, ptype_p);
  if(backward_compat_compound){
    *ggg_p = backward_compat_compound->matrices->ggg;
  }
}

PUBLIC void
export_cofold_arrays( int **f5_p,
                      int **c_p,
                      int **fML_p,
                      int **fM1_p,
                      int **fc_p,
                      int **indx_p,
                      char **ptype_p){

  wrap_array_export(f5_p, c_p, fML_p, fM1_p, fc_p, indx_p, ptype_p);
}

PUBLIC float
cofold( const char *string,
        char *structure){

  return wrap_cofold(string, structure, NULL, fold_constrained);
}

PUBLIC float
cofold_par( const char *string,
            char *structure,
            paramT *parameters,
            int is_constrained){

  return wrap_cofold(string, structure, parameters, is_constrained);
}

PUBLIC SOLUTION *
zukersubopt(const char *string) {

  return wrap_zukersubopt(string, NULL);
}

PUBLIC SOLUTION *
zukersubopt_par(const char *string,
                paramT *parameters){

  return wrap_zukersubopt(string, parameters);
}

