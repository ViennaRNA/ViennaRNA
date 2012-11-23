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

#include "utils.h"
#include "energy_par.h"
#include "fold_vars.h"
#include "params.h"
#include "loop_energies.h"
#include "constraints.h"
#include "data_structures.h"
#include "gquad.h"
#include "fold.h"

#ifdef _OPENMP
#include <omp.h>
#endif


#define PAREN
#define STACK_BULGE1      1       /* stacking energies for bulges of size 1 */
#define NEW_NINIO         1       /* new asymetry penalty */
#define MAXSECTORS        500     /* dimension for a backtrack array */
#define LOCALITY          0.      /* locality parameter for base-pairs */

#define SAME_STRAND(I,J)  (((I)>=cut_point)||((J)<cut_point))

/*
#################################
# GLOBAL VARIABLES              #
#################################
*/
PUBLIC  int uniq_ML   = 0;  /* do ML decomposition uniquely (for subopt) */

/*
#################################
# PRIVATE VARIABLES             #
#################################
*/
PRIVATE int     *indx     = NULL; /* index for moving in the triangle matrices c[] and fMl[]*/
PRIVATE int     *c        = NULL; /* energy array, given that i-j pair */
PRIVATE int     *cc       = NULL; /* linear array for calculating canonical structures */
PRIVATE int     *cc1      = NULL; /*   "     "        */
PRIVATE int     *f5       = NULL; /* energy of 5' end */
PRIVATE int     *f53      = NULL; /* energy of 5' end with 3' nucleotide not available for mismatches */
PRIVATE int     *fML      = NULL; /* multi-loop auxiliary energy array */
PRIVATE int     *fM1      = NULL; /* second ML array, only for subopt */
PRIVATE int     *fM2      = NULL; /* fM2 = multiloop region with exactly two stems, extending to 3' end        */
PRIVATE int     *Fmi      = NULL; /* holds row i of fML (avoids jumps in memory) */
PRIVATE int     *DMLi     = NULL; /* DMLi[j] holds MIN(fML[i,k]+fML[k+1,j])  */
PRIVATE int     *DMLi1    = NULL; /*             MIN(fML[i+1,k]+fML[k+1,j])  */
PRIVATE int     *DMLi2    = NULL; /*             MIN(fML[i+2,k]+fML[k+1,j])  */
PRIVATE int     Fc, FcH, FcI, FcM;  /* parts of the exterior loop energies */
PRIVATE sect    sector[MAXSECTORS]; /* stack of partial structures for backtracking */
PRIVATE char    *ptype = NULL;      /* precomputed array of pair types */
PRIVATE short   *S = NULL, *S1 = NULL;
PRIVATE paramT  *P          = NULL;
PRIVATE int     init_length = -1;
PRIVATE bondT   *base_pair2         = NULL; /* this replaces base_pair from fold_vars.c */
PRIVATE int     circular            = 0;

/* stuff needed for constrained folding */
PRIVATE int     with_gquad          = 0;
PRIVATE int     *ggg = NULL;    /* minimum free energies of the gquadruplexes */

PRIVATE int               struct_constrained  = 0;
PRIVATE hard_constraintT  *hc                 = NULL;

#ifdef _OPENMP

#pragma omp threadprivate(indx, c, cc, cc1, f5, f53, fML, fM1, fM2, Fmi,\
                          DMLi, DMLi1, DMLi2, Fc, FcH, FcI, FcM, hc,\
                          sector, ptype, S, S1, P, init_length, base_pair2, circular, struct_constrained,\
                          ggg, with_gquad)

#endif

/*
#################################
# PRIVATE FUNCTION DECLARATIONS #
#################################
*/
PRIVATE void  get_arrays(unsigned int size);
PRIVATE void  backtrack(const char *sequence, int s);
PRIVATE int   fill_arrays(const char *sequence);
PRIVATE void  fill_arrays_circ(const char *string, int *bt);
PRIVATE void  init_fold(int length, paramT *parameters);

/* deprecated functions */
/*@unused@*/
int oldLoopEnergy(int i, int j, int p, int q, int type, int type_2);
int LoopEnergy(int n1, int n2, int type, int type_2, int si1, int sj1, int sp1, int sq1);
int HairpinE(int size, int type, int si1, int sj1, const char *string);


/*
#################################
# BEGIN OF FUNCTION DEFINITIONS #
#################################
*/

/* allocate memory for folding process */
PRIVATE void init_fold(int length, paramT *parameters){

#ifdef _OPENMP
/* Explicitly turn off dynamic threads */
  omp_set_dynamic(0);
#endif

  if (length<1) nrerror("initialize_fold: argument must be greater 0");
  free_arrays();
  get_arrays((unsigned) length);
  init_length=length;

  indx = get_indx((unsigned)length);

  update_fold_params_par(parameters);
}

/*--------------------------------------------------------------------------*/

PRIVATE void get_arrays(unsigned int size){
  if(size >= (unsigned int)sqrt((double)INT_MAX))
    nrerror("get_arrays@fold.c: sequence length exceeds addressable range");

  c     = (int *) space(sizeof(int)*((size*(size+1))/2+2));
  fML   = (int *) space(sizeof(int)*((size*(size+1))/2+2));
  if (uniq_ML)
    fM1 = (int *) space(sizeof(int)*((size*(size+1))/2+2));

  f5    = (int *) space(sizeof(int)*(size+2));
  f53   = (int *) space(sizeof(int)*(size+2));
  cc    = (int *) space(sizeof(int)*(size+2));
  cc1   = (int *) space(sizeof(int)*(size+2));
  Fmi   = (int *) space(sizeof(int)*(size+1));
  DMLi  = (int *) space(sizeof(int)*(size+1));
  DMLi1 = (int *) space(sizeof(int)*(size+1));
  DMLi2 = (int *) space(sizeof(int)*(size+1));

  base_pair2 = (bondT *) space(sizeof(bondT)*(1+size/2+200)); /* add a guess of how many G's may be involved in a G quadruplex */

  /* extra array(s) for circfold() */
  if(circular) fM2 =  (int *) space(sizeof(int)*(size+2));

}

/*--------------------------------------------------------------------------*/

PUBLIC void free_arrays(void){
  if(indx)      free(indx);
  if(c)         free(c);
  if(fML)       free(fML);
  if(f5)        free(f5);
  if(f53)       free(f53);
  if(cc)        free(cc);
  if(cc1)       free(cc1);
  if(ptype)     free(ptype);
  if(fM1)       free(fM1);
  if(fM2)       free(fM2);
  if(base_pair2) free(base_pair2);
  if(Fmi)       free(Fmi);
  if(DMLi)      free(DMLi);
  if(DMLi1)     free(DMLi1);
  if(DMLi2)     free(DMLi2);
  if(P)         free(P);
  if(ggg)       free(ggg);
  if(hc)        destroy_hard_constraints(hc);

  indx = c = fML = f5 = f53 = cc = cc1 = fM1 = fM2 = Fmi = DMLi = DMLi1 = DMLi2 = ggg = NULL;
  ptype = NULL;

  hc  = NULL;

  base_pair   = NULL;
  base_pair2  = NULL;
  P           = NULL;
  init_length = 0;
}

/*--------------------------------------------------------------------------*/

PUBLIC void export_fold_arrays( int **f5_p,
                                int **c_p,
                                int **fML_p,
                                int **fM1_p,
                                int **indx_p,
                                char **ptype_p){
  /* make the DP arrays available to routines such as subopt() */
  *f5_p     = f5;
  *c_p      = c;
  *fML_p    = fML;
  *fM1_p    = fM1;
  *indx_p   = indx;
  *ptype_p  = ptype;
}

PUBLIC void export_fold_arrays_par( int **f5_p,
                                    int **c_p,
                                    int **fML_p,
                                    int **fM1_p,
                                    int **indx_p,
                                    char **ptype_p,
                                    paramT **P_p){
  export_fold_arrays(f5_p, c_p, fML_p, fM1_p, indx_p,ptype_p);
  *P_p = P;
}

PUBLIC void export_circfold_arrays( int *Fc_p,
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
  /* make the DP arrays available to routines such as subopt() */
  *f5_p     = f5;
  *c_p      = c;
  *fML_p    = fML;
  *fM1_p    = fM1;
  *fM2_p    = fM2;
  *Fc_p     = Fc;
  *FcH_p    = FcH;
  *FcI_p    = FcI;
  *FcM_p    = FcM;
  *indx_p   = indx;
  *ptype_p  = ptype;
}

PUBLIC void export_circfold_arrays_par( int *Fc_p,
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
                                    paramT **P_p){
  export_circfold_arrays( Fc_p,
                          FcH_p,
                          FcI_p,
                          FcM_p,
                          fM2_p,
                          f5_p,
                          c_p,
                          fML_p,
                          fM1_p,
                          indx_p,
                          ptype_p);
  *P_p = P;
}

PUBLIC float fold(const char *string, char *structure){
  return fold_par(string, structure, NULL, fold_constrained, 0);
}

PUBLIC float circfold(const char *string, char *structure){
  return fold_par(string, structure, NULL, fold_constrained, 1);
}

PUBLIC float fold_par(const char *string,
                      char *structure,
                      paramT *parameters,
                      int is_constrained,
                      int is_circular){

  int i, length, energy, s;
  unsigned int constraint_options;

  s                   = 0;
  circular            = is_circular;
  struct_constrained  = is_constrained;
  length              = (int) strlen(string);
  constraint_options  = 0;
  /* prepare constraint options */
  if(struct_constrained && structure)
    constraint_options |= VRNA_CONSTRAINT_DB;

#ifdef _OPENMP
  init_fold(length, parameters);
#else
  if (parameters) init_fold(length, parameters);
  else if (length>init_length) init_fold(length, parameters);
  else if (fabs(P->temperature - temperature)>1e-6) update_fold_params_par(parameters);
#endif

  with_gquad  = P->model_details.gquad;
  S     = get_sequence_encoding(string, 0, &(P->model_details));
  S1    = get_sequence_encoding(string, 1, &(P->model_details));
  ptype = get_ptypes(S, &(P->model_details), 0);
  hc    = get_hard_constraints( (const char *)structure,
                                (unsigned int)length,
                                ptype,
                                TURN,
                                constraint_options);

  energy = fill_arrays(string);

  if(circular){
    fill_arrays_circ(string, &s);
    energy = Fc;
  }
  backtrack(string, s);

#ifdef PAREN
  parenthesis_structure(structure, base_pair2, length);
#else
  letter_structure(structure, base_pair2, length);
#endif

  /*
  *  Backward compatibility:
  *  This block may be removed if deprecated functions
  *  relying on the global variable "base_pair" vanish from within the package!
  */
  base_pair = base_pair2;

  free(S); free(S1);

  if (backtrack_type=='C')
    return (float) c[indx[length]+1]/100.;
  else if (backtrack_type=='M')
    return (float) fML[indx[length]+1]/100.;
  else
    return (float) energy/100.;
}

/**
*** fill "c", "fML" and "f5" arrays and return  optimal energy
**/
PRIVATE int fill_arrays(const char *string) {

  int   i, j, ij, p, q, k, length, energy, en, mm5, mm3;
  int   decomp, new_fML, max_separation;
  int   no_close, type, type_2, tt;

  int   dangle_model, noGUclosure, noLP, with_gquads;
  int   *rtype;

  int   hc_decompose;
  char  *hard_constraints = hc->matrix;
  int   *hc_up_ext        = hc->up_ext;
  int   *hc_up_hp         = hc->up_hp;
  int   *hc_up_int        = hc->up_int;
  int   *hc_up_ml         = hc->up_ml;

  dangle_model  = P->model_details.dangles;
  noGUclosure   = P->model_details.noGUclosure;
  noLP          = P->model_details.noLP;
  rtype         = &(P->model_details.rtype[0]);

  length = (int) strlen(string);

  max_separation = (int) ((1.-LOCALITY)*(double)(length-2)); /* not in use */

  if(with_gquad)
    ggg = get_gquad_matrix(S, P);


  for (j=1; j<=length; j++) {
    Fmi[j]=DMLi[j]=DMLi1[j]=DMLi2[j]=INF;
  }

  for (j = 1; j<=length; j++)
    for (i=(j>TURN?(j-TURN):1); i<j; i++) {
      c[indx[j]+i] = fML[indx[j]+i] = INF;
      if (uniq_ML) fM1[indx[j]+i] = INF;
    }

  if (length <= TURN) return 0;

  for (i = length-TURN-1; i >= 1; i--) { /* i,j in [1..length] */

    for (j = i+TURN+1; j <= length; j++) {
      ij            = indx[j]+i;
      type          = ptype[ij];
      hc_decompose  = hard_constraints[ij];
      energy        = INF;

      no_close = (((type==3)||(type==4))&&noGUclosure);

      if (j-i-1 > max_separation) type = 0;  /* forces locality degree */

      if (hc_decompose) {   /* we have a pair */
        int new_c=0, stackEnergy=INF;

        /* hairpin ----------------------------------------------*/

        /* CONSTRAINED HAIRPIN start */
        int u = j - i - 1;
        /* is this base pair allowed to close a hairpin loop ? */
        if(hc_decompose & IN_HP_LOOP){
          /* are all nucleotides in the loop allowed to be unpaired ? */
          if(hc_up_hp[i+1] >= u)
            new_c = E_Hairpin(u, type, S1[i+1], S1[j-1], string+i-1, P);
        }
        /* CONSTRAINED HAIRPIN end */

        /*--------------------------------------------------------
          check for elementary structures involving more than one
          closing pair.
          --------------------------------------------------------*/

        /* CONSTRAINED INTERIOR LOOP start */
        if(hc_decompose & IN_INT_LOOP){
          for(q = j - 1; q >= MAX2(i+TURN+2, j - MAXLOOP - 1); q--){
            int j_q = j - q - 1;

            if(hc_up_int[q+1] < j_q) break;

            int maxp = MAX2(i + 1, i + 1 + MAXLOOP - j_q);
            maxp = MIN2(maxp, q - TURN);
            int pq = indx[q] + i + 1;
            int p_i = 0;
            
            maxp = MIN2(maxp, i+1+hc_up_int[i+1]);
            for(p = i+1; p <= maxp; p++, pq++, p_i++){

              /* discard this configuration if (p,q) is not allowed to be enclosed pair of an interior loop */
              if(hard_constraints[pq] & IN_INT_LOOP_ENC){

                type_2 = rtype[ptype[pq]];

                if (noGUclosure)
                  if (no_close||(type_2==3)||(type_2==4))
                    if ((p>i+1)||(q<j-1)) continue;  /* continue unless stack */

                energy = E_IntLoop(p_i, j_q, type, type_2,
                                    S1[i+1], S1[j-1], S1[p-1], S1[q+1], P);

                if ((p==i+1)&&(j==q+1)) stackEnergy = energy; /* remember stack energy */
                energy += c[pq];
                new_c = MIN2(energy, new_c);

              }
            } /* end q-loop */
          } /* end p-loop */
        }

        /* CONSTRAINED INTERIOR LOOP end */

        /* multi-loop decomposition ------------------------*/

        /* CONSTRAINED MULTIBRANCH LOOP start */
        if(hc_decompose & IN_MB_LOOP){
          int MLenergy;
          decomp = DMLi1[j-1];
          tt = rtype[type];
          switch(dangle_model){
            /* no dangles */
            case 0:   decomp += E_MLstem(tt, -1, -1, P);
                      break;

            /* double dangles */
            case 2:   decomp += E_MLstem(tt, S1[j-1], S1[i+1], P);
                      break;

            /* normal dangles, aka dangles = 1 || 3 */
            default:  decomp += E_MLstem(tt, -1, -1, P);
                      if(hc_up_ml[i+1]){
                        en = DMLi2[j-1] + E_MLstem(tt, -1, S1[i+1], P) + P->MLbase;
                        decomp = MIN2(decomp, en);
                      }
                      if(hc_up_ml[j-1] && hc_up_ml[i+1]){
                        en = DMLi2[j-2] + E_MLstem(tt, S1[j-1], S1[i+1], P) + 2*P->MLbase;
                        decomp = MIN2(decomp, en);
                      }
                      if(hc_up_ml[j-1]){
                        en = DMLi1[j-2] + E_MLstem(tt, S1[j-1], -1, P) + P->MLbase;
                        decomp = MIN2(decomp, en);
                      }
                      break;
          }
          MLenergy = decomp + P->MLclosing;
          new_c = MIN2(new_c, MLenergy);

          /* coaxial stacking of (i.j) with (i+1.k) or (k+1.j-1) */

          if (dangle_model==3) {
            decomp = INF;
            int i1k, k1j1;
            k1j1  = indx[j-1] + i + 2 + TURN + 1;
            for (k = i+2+TURN; k < j-2-TURN; k++, k1j1++){
              i1k   = indx[k] + i + 1;
              if(hard_constraints[i1k] & IN_MB_LOOP_ENC){
                type_2  = rtype[ptype[i1k]];
                en      = c[i1k]+P->stack[type][type_2]+fML[k1j1];
                decomp  = MIN2(decomp, en);
              }
              if(hard_constraints[k1j1] & IN_MB_LOOP_ENC){
                type_2  = rtype[ptype[k1j1]];
                en      = c[k1j1]+P->stack[type][type_2]+fML[i1k];
                decomp  = MIN2(decomp, en);
              }
            }
            /* no TermAU penalty if coax stack */
            decomp += 2*P->MLintern[1] + P->MLclosing;
            new_c = MIN2(new_c, decomp);
          }
        }
        /* CONSTRAINED MULTIBRANCH LOOP end */

        if(with_gquad){
          /* include all cases where a g-quadruplex may be enclosed by base pair (i,j) */
          if (!no_close) {
            tt = rtype[type];
            energy = E_GQuad_IntLoop(i, j, type, S1, ggg, indx, P);
            new_c = MIN2(new_c, energy);
          }
        }

        new_c = MIN2(new_c, cc1[j-1]+stackEnergy);
        cc[j] = new_c;
        if (noLP)
          c[ij] = cc1[j-1]+stackEnergy;
        else
          c[ij] = cc[j];

      } /* end >> if (pair) << */

      else c[ij] = INF;

      /* done with c[i,j], now compute fML[i,j] and fM1[i,j] */

      /* (i,j) + MLstem ? */
      new_fML = INF;
      if(hc_decompose & IN_MB_LOOP_ENC){
        new_fML = c[ij];
        switch(dangle_model){
          case 2:   new_fML += E_MLstem(type, (i==1) ? S1[length] : S1[i-1], S1[j+1], P);
                    break;
          default:  new_fML += E_MLstem(type, -1, -1, P);
                    break;
        }
        if(uniq_ML) fM1[ij] = new_fML;
      }

      if(with_gquad){
        new_fML = MIN2(new_fML, ggg[indx[j] + i] + E_MLstem(0, -1, -1, P));
      }

      if(hc_up_ml[j-1] && uniq_ML){
        en      = fM1[indx[j-1]+i] + P->MLbase;
        fM1[ij] = MIN2(fM1[ij], en);
      }

      /* free ends ? -----------------------------------------*/
      /*  we must not just extend 3'/5' end by unpaired nucleotides if
      *   dangle_model == 1, this could lead to d5+d3 contributions were
      *   mismatch must be taken!
      */
      switch(dangle_model){
        /* no dangles */
        case 0:   

        /* double dangles */
        case 2:   if(hc_up_ml[i])
                    new_fML = MIN2(new_fML, fML[ij+1]+P->MLbase);
                  if(hc_up_ml[j])
                    new_fML = MIN2(fML[indx[j-1]+i]+P->MLbase, new_fML);
                  break;

        /* normal dangles, aka dangle_model = 1 || 3 */
        default:  mm5 = ((i>1) || circular) ? S1[i] : -1;
                  mm3 = ((j<length) || circular) ? S1[j] : -1;
                  new_fML = MIN2(new_fML, fML[ij+1] + P->MLbase);
                  new_fML = MIN2(new_fML, fML[indx[j-1]+i] + P->MLbase);
                  if(hard_constraints[ij+1] & IN_MB_LOOP_ENC){
                    if(hc_up_ml[i]){
                      tt = ptype[ij+1];
                      new_fML = MIN2(new_fML, c[ij+1] + E_MLstem(tt, mm5, -1, P) + P->MLbase);
                    }
                  }
                  if(hard_constraints[indx[j-1]+i] & IN_MB_LOOP_ENC){
                    if(hc_up_ml[j]){
                      tt = ptype[indx[j-1]+i];
                      new_fML = MIN2(new_fML, c[indx[j-1]+i] + E_MLstem(tt, -1, mm3, P) + P->MLbase);
                    }
                  }
                  if(hard_constraints[indx[j-1]+i+1] & IN_MB_LOOP_ENC){
                    if(hc_up_ml[i] && hc_up_ml[j]){
                      tt = ptype[indx[j-1]+i+1];
                      new_fML = MIN2(new_fML, c[indx[j-1]+i+1] + E_MLstem(tt, mm5, mm3, P) + 2*P->MLbase);
                    }
                  }
                  break;
      }

      /* modular decomposition -------------------------------*/
      int k1j = indx[j] + i + TURN + 2;
      for (decomp = INF, k = i + 1 + TURN; k <= j - 2 - TURN; k++, k1j++)
        decomp = MIN2(decomp, Fmi[k]+fML[k1j]);
      DMLi[j] = decomp;               /* store for use in ML decompositon */
      new_fML = MIN2(new_fML,decomp);
      /* coaxial stacking */
      if (dangle_model==3) {
        /* additional ML decomposition as two coaxially stacked helices */
        int ik, k1j;
        for (k1j = indx[j]+i+TURN+2, decomp = INF, k = i+1+TURN; k <= j-2-TURN; k++, k1j++){
          ik = indx[k]+i;
          if((hard_constraints[ik] & IN_MB_LOOP_ENC) && (hard_constraints[k1j] & IN_MB_LOOP_ENC)){
            type    = rtype[ptype[ik]];
            type_2  = rtype[ptype[k1j]];
            decomp  = MIN2(decomp, c[ik] + c[k1j] + P->stack[type][type_2]);
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
        new_fML = MIN2(new_fML, decomp);
      }
      fML[ij] = Fmi[j] = new_fML;     /* substring energy */

    }

    {
      int *FF; /* rotate the auxilliary arrays */
      FF = DMLi2; DMLi2 = DMLi1; DMLi1 = DMLi; DMLi = FF;
      FF = cc1; cc1=cc; cc=FF;
      for (j=1; j<=length; j++) {cc[j]=Fmi[j]=DMLi[j]=INF; }
    }
  }

  /* calculate energies of 5' and 3' fragments */

  f5[0] = 0;
  for(i = 1; i <= TURN + 1; i++){
    f5[i] = (hc_up_ext[i]) ? f5[i-1] : INF;
  }
  /* duplicated code may be faster than conditions inside loop ;) */
  switch(dangle_model){
    /* dont use dangling end and mismatch contributions at all */
    case 0:   for(j=TURN+2; j<=length; j++){
                if(hc_up_ext[j])
                  f5[j] = f5[j-1];
                for (i=j-TURN-1; i>1; i--){
                  ij = indx[j]+i;
                  if(!(hard_constraints[ij] & IN_EXT_LOOP)) continue;

                  if(with_gquad){
                    f5[j] = MIN2(f5[j], f5[i-1] + ggg[indx[j]+i]);
                  }

                  en    = f5[i-1] + c[ij] + E_ExtLoop(ptype[ij], -1, -1, P);
                  f5[j] = MIN2(f5[j], en);
                }
                ij = indx[j] + 1;
                if(!(hard_constraints[ij] & IN_EXT_LOOP)) continue;

                if(with_gquad){
                  f5[j] = MIN2(f5[j], ggg[indx[j]+1]);
                }

                en    = c[ij] + E_ExtLoop(ptype[ij], -1, -1, P);
                f5[j] = MIN2(f5[j], en);
              }
              break;

    /* always use dangles on both sides */
    case 2:   for(j=TURN+2; j<length; j++){
                if(hc_up_ext[j])
                  f5[j] = f5[j-1];
                for (i=j-TURN-1; i>1; i--){
                  ij = indx[j] + i;
                  if(!(hard_constraints[ij] & IN_EXT_LOOP)) continue;

                  if(with_gquad){
                    f5[j] = MIN2(f5[j], f5[i-1] + ggg[indx[j]+i]);
                  }

                  en    = f5[i-1] + c[ij] + E_ExtLoop(ptype[ij], S1[i-1], S1[j+1], P);
                  f5[j] = MIN2(f5[j], en);
                }
                ij = indx[j] + 1;
                if(!(hard_constraints[ij] & IN_EXT_LOOP)) continue;

                if(with_gquad){
                  f5[j] = MIN2(f5[j], ggg[indx[j]+1]);
                }

                en    = c[ij] + E_ExtLoop(ptype[ij], -1, S1[j+1], P);
                f5[j] = MIN2(f5[j], en);
              }
              if(hc_up_ext[length])
                f5[length] = f5[length-1];
              for (i=length-TURN-1; i>1; i--){
                ij = indx[length] + i;
                if(!(hard_constraints[ij] & IN_EXT_LOOP)) continue;

                if(with_gquad){
                  f5[length] = MIN2(f5[length], f5[i-1] + ggg[indx[length]+i]);
                }

                en          = f5[i-1] + c[ij] + E_ExtLoop(ptype[ij], S1[i-1], -1, P);
                f5[length]  = MIN2(f5[length], en);
              }
              ij = indx[length] + 1;
              if(!(hard_constraints[ij] & IN_EXT_LOOP)) break;

              if(with_gquad){
                f5[length] = MIN2(f5[length], ggg[indx[length]+1]);
              }

              en          = c[ij] + E_ExtLoop(ptype[ij], -1, -1, P);
              f5[length]  = MIN2(f5[length], en);
              break;

    /* normal dangles, aka dangle_model = 1 || 3 */
    default:  for(j=TURN+2; j<=length; j++){
                if(hc_up_ext[j])
                  f5[j] = f5[j-1];
                for (i=j-TURN-1; i>1; i--){
                  ij = indx[j] + i;
                  if(hard_constraints[ij] & IN_EXT_LOOP){

                    if(with_gquad){
                      f5[j] = MIN2(f5[j], f5[i-1] + ggg[indx[j]+i]);
                    }

                    type  = ptype[ij];
                    en    = f5[i-1] + c[ij] + E_ExtLoop(type, -1, -1, P);
                    f5[j] = MIN2(f5[j], en);
                    if(hc_up_ext[i-1]){
                      en    = f5[i-2] + c[ij] + E_ExtLoop(type, S1[i-1], -1, P);
                      f5[j] = MIN2(f5[j], en);
                    }
                  }
                  ij = indx[j-1] + i;
                  if(hard_constraints[ij] & IN_EXT_LOOP){
                    if(hc_up_ext[j]){
                      type  = ptype[ij];
                      en    = f5[i-1] + c[ij] + E_ExtLoop(type, -1, S1[j], P);
                      f5[j] = MIN2(f5[j], en);
                      if(hc_up_ext[i-1]){
                        en    = f5[i-2] + c[ij] + E_ExtLoop(type, S1[i-1], S1[j], P);
                        f5[j] = MIN2(f5[j], en);
                      }
                    }
                  }
                }
                ij = indx[j] + 1;
                if(hard_constraints[ij] & IN_EXT_LOOP){

                  if(with_gquad){
                    f5[j] = MIN2(f5[j], ggg[indx[j]+1]);
                  }

                  type  = ptype[ij];
                  en    = c[ij] + E_ExtLoop(type, -1, -1, P);
                  f5[j] = MIN2(f5[j], en);
                }
                ij = indx[j-1] + 1;
                if(hard_constraints[ij] & IN_EXT_LOOP){
                  if(hc_up_ext[j]){
                    type  = ptype[ij];
                    en    = c[ij] + E_ExtLoop(type, -1, S1[j], P);
                    f5[j] = MIN2(f5[j], en);
                  }
                }
              }
  }
  return f5[length];
}

#include "circfold.inc"

/**
*** trace back through the "c", "f5" and "fML" arrays to get the
*** base pairing list. No search for equivalent structures is done.
*** This is fast, since only few structure elements are recalculated.
***
*** normally s=0.
*** If s>0 then s items have been already pushed onto the sector stack
**/
PRIVATE void backtrack(const char *string, int s) {
  int   i, j, ij, k, mm3, length, energy, en, new;
  int   no_close, type, type_2, tt, minq, maxq, c0, c1, c2, c3;
  int   b=0;
  int   dangle_model  = P->model_details.dangles;
  int   noLP          = P->model_details.noLP;
  int   noGUclosure   = P->model_details.noGUclosure;
  int   *rtype        = &(P->model_details.rtype[0]);

  length = strlen(string);
  if (s==0) {
    sector[++s].i = 1;
    sector[s].j = length;
    sector[s].ml = (backtrack_type=='M') ? 1 : ((backtrack_type=='C')? 2: 0);
  }
  while (s>0) {
    int ml, fij, fi, cij, traced, i1, j1, p, q, jj=0, gq=0;
    int canonical = 1;     /* (i,j) closes a canonical structure */
    i  = sector[s].i;
    j  = sector[s].j;
    ml = sector[s--].ml;   /* ml is a flag indicating if backtracking is to
                              occur in the fML- (1) or in the f-array (0) */
    if (ml==2) {
      base_pair2[++b].i = i;
      base_pair2[b].j   = j;
      goto repeat1;
    }

    else if(ml==7) { /* indicates that i,j are enclosing a gquadruplex */
      /* actually, do something here */
    }

    if (j < i+TURN+1) continue; /* no more pairs in this interval */

    fij = (ml == 1)? fML[indx[j]+i] : f5[j];
    fi  = (ml == 1)?(fML[indx[j-1]+i]+P->MLbase): f5[j-1];

    if (fij == fi) {  /* 3' end is unpaired */
      sector[++s].i = i;
      sector[s].j   = j-1;
      sector[s].ml  = ml;
      continue;
    }

    if (ml == 0) { /* backtrack in f5 */
      switch(dangle_model){
        case 0:   /* j is paired. Find pairing partner */
                  for(k=j-TURN-1,traced=0; k>=1; k--){

                    if(with_gquad){
                      if(fij == f5[k-1] + ggg[indx[j]+k]){
                        /* found the decomposition */
                        traced = j; jj = k - 1; gq = 1;
                        break;
                      }
                    }

                    type = ptype[indx[j]+k];
                    if(type)
                      if(fij == E_ExtLoop(type, -1, -1, P) + c[indx[j]+k] + f5[k-1]){
                        traced=j; jj = k-1;
                        break;
                      }
                  }
                  break;

        case 2:   mm3 = (j<length) ? S1[j+1] : -1;
                  for(k=j-TURN-1,traced=0; k>=1; k--){

                    if(with_gquad){
                      if(fij == f5[k-1] + ggg[indx[j]+k]){
                        /* found the decomposition */
                        traced = j; jj = k - 1; gq = 1;
                        break;
                      }
                    }

                    type = ptype[indx[j]+k];
                    if(type)
                      if(fij == E_ExtLoop(type, (k>1) ? S1[k-1] : -1, mm3, P) + c[indx[j]+k] + f5[k-1]){
                        traced=j; jj = k-1;
                        break;
                      }
                  }
                  break;

        default:  for(traced = 0, k=j-TURN-1; k>1; k--){

                    if(with_gquad){
                      if(fij == f5[k-1] + ggg[indx[j]+k]){
                        /* found the decomposition */
                        traced = j; jj = k - 1; gq = 1;
                        break;
                      }
                    }

                    type = ptype[indx[j] + k];
                    if(type){
                      en = c[indx[j] + k];
                      if(fij == f5[k-1] + en + E_ExtLoop(type, -1, -1, P)){
                        traced = j;
                        jj = k-1;
                        break;
                      }
                      if(fij == f5[k-2] + en + E_ExtLoop(type, S1[k-1], -1, P)){
                        traced = j;
                        jj = k-2;
                        break;
                      }
                    }
                    type = ptype[indx[j-1] + k];
                    if(type){
                      en = c[indx[j-1] + k];
                      if(fij == f5[k-1] + en + E_ExtLoop(type, -1, S1[j], P)){
                        traced = j-1;
                        jj = k-1;
                        break;
                      }
                      if(fij == f5[k-2] + en + E_ExtLoop(type, S1[k-1], S1[j], P)){
                        traced = j-1;
                        jj = k-2;
                        break;
                      }
                    }
                  }
                  if(!traced){

                    if(with_gquad){
                      if(fij == ggg[indx[j]+1]){
                        /* found the decomposition */
                        traced = j; jj = 0; gq = 1;
                        break;
                      }
                    }

                    type = ptype[indx[j]+1];
                    if(type){
                      if(fij == c[indx[j]+1] + E_ExtLoop(type, -1, -1, P)){
                        traced = j;
                        jj = 0;
                        break;
                      }
                    }
                    type = ptype[indx[j-1]+1];
                    if(type){
                      if(fij == c[indx[j-1]+1] + E_ExtLoop(type, -1, S1[j], P)){
                        traced = j-1;
                        jj = 0;
                        break;
                      }
                    }
                  }
                  break;
      }

      if (!traced){
        fprintf(stderr, "%s\n", string);
        nrerror("backtrack failed in f5");
      }
      /* push back the remaining f5 portion */
      sector[++s].i = 1;
      sector[s].j   = jj;
      sector[s].ml  = ml;

      /* trace back the base pair found */
      i=k; j=traced;

      if(with_gquad && gq){
        /* goto backtrace of gquadruplex */
        goto repeat_gquad;
      }

      base_pair2[++b].i = i;
      base_pair2[b].j   = j;
      goto repeat1;
    }
    else { /* trace back in fML array */
      if (fML[indx[j]+i+1]+P->MLbase == fij) { /* 5' end is unpaired */
        sector[++s].i = i+1;
        sector[s].j   = j;
        sector[s].ml  = ml;
        continue;
      }

      ij  = indx[j]+i;

      if(with_gquad){
        if(fij == ggg[ij] + E_MLstem(0, -1, -1, P)){
          /* go to backtracing of quadruplex */
          goto repeat_gquad;
        }
      }

      tt  = ptype[ij];
      en  = c[ij];
      switch(dangle_model){
        case 0:   if(fij == en + E_MLstem(tt, -1, -1, P)){
                    base_pair2[++b].i = i;
                    base_pair2[b].j   = j;
                    goto repeat1;
                  }
                  break;

        case 2:   if(fij == en + E_MLstem(tt, S1[i-1], S1[j+1], P)){
                    base_pair2[++b].i = i;
                    base_pair2[b].j   = j;
                    goto repeat1;
                  }
                  break;

        default:  if(fij == en + E_MLstem(tt, -1, -1, P)){
                    base_pair2[++b].i = i;
                    base_pair2[b].j   = j;
                    goto repeat1;
                  }
                  tt = ptype[ij+1];
                  if(fij == c[ij+1] + E_MLstem(tt, S1[i], -1, P) + P->MLbase){
                    base_pair2[++b].i = ++i;
                    base_pair2[b].j   = j;
                    goto repeat1;
                  }
                  tt = ptype[indx[j-1]+i];
                  if(fij == c[indx[j-1]+i] + E_MLstem(tt, -1, S1[j], P) + P->MLbase){
                    base_pair2[++b].i = i;
                    base_pair2[b].j   = --j;
                    goto repeat1;
                  }
                  tt = ptype[indx[j-1]+i+1];
                  if(fij == c[indx[j-1]+i+1] + E_MLstem(tt, S1[i], S1[j], P) + 2*P->MLbase){
                    base_pair2[++b].i = ++i;
                    base_pair2[b].j   = --j;
                    goto repeat1;
                  }
                  break;
      }

      for(k = i + 1 + TURN; k <= j - 2 - TURN; k++)
        if(fij == (fML[indx[k]+i]+fML[indx[j]+k+1]))
          break;

      if ((dangle_model==3)&&(k > j - 2 - TURN)) { /* must be coax stack */
        ml = 2;
        for (k = i+1+TURN; k <= j - 2 - TURN; k++) {
          type    = rtype[ptype[indx[k]+i]];
          type_2  = rtype[ptype[indx[j]+k+1]];
          if (type && type_2)
            if (fij == c[indx[k]+i]+c[indx[j]+k+1]+P->stack[type][type_2]+
                       2*P->MLintern[1])
              break;
        }
      }
      sector[++s].i = i;
      sector[s].j   = k;
      sector[s].ml  = ml;
      sector[++s].i = k+1;
      sector[s].j   = j;
      sector[s].ml  = ml;

      if (k>j-2-TURN) nrerror("backtrack failed in fML");
      continue;
    }

  repeat1:

    /*----- begin of "repeat:" -----*/
    ij = indx[j]+i;
    if (canonical)  cij = c[ij];

    type = ptype[ij];

    if (noLP)
      if (cij == c[ij]){
        /* (i.j) closes canonical structures, thus
           (i+1.j-1) must be a pair                */
        type_2 = ptype[indx[j-1]+i+1]; type_2 = rtype[type_2];
        cij -= P->stack[type][type_2];
        base_pair2[++b].i = i+1;
        base_pair2[b].j   = j-1;
        i++; j--;
        canonical=0;
        goto repeat1;
      }
    canonical = 1;


    no_close = (((type==3)||(type==4))&&noGUclosure);
    if (no_close) {
      if (cij == FORBIDDEN) continue;
    } else
      if (cij == E_Hairpin(j-i-1, type, S1[i+1], S1[j-1],string+i-1, P))
        continue;

    for (p = i+1; p <= MIN2(j-2-TURN,i+MAXLOOP+1); p++) {
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
        energy = E_IntLoop(p-i-1, j-q-1, type, type_2,
                            S1[i+1], S1[j-1], S1[p-1], S1[q+1], P);

        new = energy+c[indx[q]+p];
        traced = (cij == new);
        if (traced) {
          base_pair2[++b].i = p;
          base_pair2[b].j   = q;
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
      if(backtrack_GQuad_IntLoop(cij - bonus, i, j, type, S, ggg, indx, &p, &q, P)){
        i = p; j = q;
        goto repeat_gquad;
      }
    }

    sector[s+1].ml  = sector[s+2].ml = 1;

    switch(dangle_model){
      case 0:   en = cij - E_MLstem(tt, -1, -1, P) - P->MLclosing;
                for(k = i+2+TURN; k < j-2-TURN; k++){
                  if(en == fML[indx[k]+i+1] + fML[indx[j-1]+k+1])
                    break;
                }
                break;

      case 2:   en = cij - E_MLstem(tt, S1[j-1], S1[i+1], P) - P->MLclosing;
                for(k = i+2+TURN; k < j-2-TURN; k++){
                    if(en == fML[indx[k]+i+1] + fML[indx[j-1]+k+1])
                      break;
                }
                break;

      default:  for(k = i+2+TURN; k < j-2-TURN; k++){
                  en = cij - P->MLclosing;
                  if(en == fML[indx[k]+i+1] + fML[indx[j-1]+k+1] + E_MLstem(tt, -1, -1, P)){
                    break;
                  }
                  else if(en == fML[indx[k]+i+2] + fML[indx[j-1]+k+1] + E_MLstem(tt, -1, S1[i+1], P) + P->MLbase){
                    i1 = i+2;
                    break;
                  }
                  else if(en == fML[indx[k]+i+1] + fML[indx[j-2]+k+1] + E_MLstem(tt, S1[j-1], -1, P) + P->MLbase){
                    j1 = j-2;
                    break;
                  }
                  else if(en == fML[indx[k]+i+2] + fML[indx[j-2]+k+1] + E_MLstem(tt, S1[j-1], S1[i+1], P) + 2*P->MLbase){
                    i1 = i+2;
                    j1 = j-2;
                    break;
                  }
                  /* coaxial stacking of (i.j) with (i+1.k) or (k.j-1) */
                  /* use MLintern[1] since coax stacked pairs don't get TerminalAU */
                  if(dangle_model == 3){
                    type_2 = rtype[ptype[indx[k]+i+1]];
                    if (type_2) {
                      en = c[indx[k]+i+1]+P->stack[type][type_2]+fML[indx[j-1]+k+1];
                      if (cij == en+2*P->MLintern[1]+P->MLclosing) {
                        ml = 2;
                        sector[s+1].ml  = 2;
                        traced = 1;
                        break;
                      }
                    }
                    type_2 = rtype[ptype[indx[j-1]+k+1]];
                    if (type_2) {
                      en = c[indx[j-1]+k+1]+P->stack[type][type_2]+fML[indx[k]+i+1];
                      if (cij == en+2*P->MLintern[1]+P->MLclosing) {
                        sector[s+2].ml = 2;
                        traced = 1;
                        break;
                      }
                    }
                  }
                }
                break;
    }

    if (k<=j-3-TURN) { /* found the decomposition */
      sector[++s].i = i1;
      sector[s].j   = k;
      sector[++s].i = k+1;
      sector[s].j   = j1;
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
        sector[++s].i = i1;
        sector[s].j   = j1;
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
          base_pair2[++b].i = i+a;
          base_pair2[b].j   = i+a;
          base_pair2[++b].i = i+L+l[0]+a;
          base_pair2[b].j   = i+L+l[0]+a;
          base_pair2[++b].i = i+L+l[0]+L+l[1]+a;
          base_pair2[b].j   = i+L+l[0]+L+l[1]+a;
          base_pair2[++b].i = i+L+l[0]+L+l[1]+L+l[2]+a;
          base_pair2[b].j   = i+L+l[0]+L+l[1]+L+l[2]+a;
        }
        goto repeat_gquad_exit;
      }
      nrerror("backtracking failed in repeat_gquad");
    }
  repeat_gquad_exit:
    asm("nop");

  } /* end of infinite while loop */

  base_pair2[0].i = b;    /* save the total number of base pairs */
}

PUBLIC char *backtrack_fold_from_pair(char *sequence, int i, int j) {
  char *structure;
  sector[1].i  = i;
  sector[1].j  = j;
  sector[1].ml = 2;
  base_pair2[0].i=0;
  S   = get_sequence_encoding(sequence, 0, &(P->model_details));
  S1  = get_sequence_encoding(sequence, 1, &(P->model_details));
  backtrack(sequence, 1);
  structure = (char *) space((strlen(sequence)+1)*sizeof(char));
  parenthesis_structure(structure, base_pair2, strlen(sequence));
  free(S);free(S1);
  return structure;
}

/*---------------------------------------------------------------------------*/

PUBLIC void letter_structure(char *structure, bondT *bp, int length){
  int   n, k, x, y;
  char  alpha[] = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz";

  for (n = 0; n < length; structure[n++] = ' ');
  structure[length] = '\0';

  for (n = 0, k = 1; k <= bp[0].i; k++) {
    y = bp[k].j;
    x = bp[k].i;
    if (x-1 > 0 && y+1 <= length) {
      if (structure[x-2] != ' ' && structure[y] == structure[x-2]) {
        structure[x-1] = structure[x-2];
        structure[y-1] = structure[x-1];
        continue;
      }
    }
    if (structure[x] != ' ' && structure[y-2] == structure[x]) {
      structure[x-1] = structure[x];
      structure[y-1] = structure[x-1];
      continue;
    }
    n++;
    structure[x-1] = alpha[n-1];
    structure[y-1] = alpha[n-1];
  }
}

/*---------------------------------------------------------------------------*/

PUBLIC void parenthesis_structure(char *structure, bondT *bp, int length){
  int n, k;

  for (n = 0; n < length; structure[n++] = '.');
  structure[length] = '\0';

  for (k = 1; k <= bp[0].i; k++){

    if(bp[k].i == bp[k].j){ /* Gquad bonds are marked as bp[i].i == bp[i].j */
      structure[bp[k].i-1] = '+';
    } else { /* the following ones are regular base pairs */
      structure[bp[k].i-1] = '(';
      structure[bp[k].j-1] = ')';
    }
  }
}

PUBLIC void parenthesis_zuker(char *structure, bondT *bp, int length){
  int k, i, j, temp;

  for (k = 0; k < length; structure[k++] = '.');
  structure[length] = '\0';

  for (k = 1; k <= bp[0].i; k++) {
    i=bp[k].i;
    j=bp[k].j;
    if (i>length) i-=length;
    if (j>length) j-=length;
    if (i>j) {
      temp=i; i=j; j=temp;
    }
    if(i == j){ /* Gquad bonds are marked as bp[i].i == bp[i].j */
      structure[i-1] = '+';
    } else { /* the following ones are regular base pairs */
      structure[i-1] = '(';
      structure[j-1] = ')';
    }
  }
}


/*---------------------------------------------------------------------------*/

PUBLIC void update_fold_params(void){
  update_fold_params_par(NULL);
}

PUBLIC void update_fold_params_par(paramT *parameters){
  if(P) free(P);
  if(parameters){
    P = get_parameter_copy(parameters);
  } else {
    model_detailsT md;
    set_model_details(&md);
    P = get_scaled_parameters(temperature, md);
  }

  fill_pair_matrices(&(P->model_details));
  if (init_length < 0) init_length=0;
}

PUBLIC void assign_plist_from_db(plist **pl, const char *struc, float pr){
  /* convert bracket string to plist */
  short *pt;
  int i, k = 0, size, n;
  plist *gpl, *ptr;

  size  = strlen(struc);
  n     = 2;

  pt  = make_pair_table(struc);
  *pl = (plist *)space(n*size*sizeof(plist));
  for(i = 1; i < size; i++){
    if(pt[i]>i){
      (*pl)[k].i      = i;
      (*pl)[k].j      = pt[i];
      (*pl)[k].p      = pr;
      (*pl)[k++].type = 0;
    }
  }

  gpl = get_plist_gquad_from_db(struc, pr);
  for(ptr = gpl; ptr->i != 0; ptr++){
    if (k == n * size - 1){
      n *= 2;
      *pl = (plist *)xrealloc(*pl, n * size * sizeof(plist));
    }
    (*pl)[k].i      = ptr->i;
    (*pl)[k].j      = ptr->j;
    (*pl)[k].p       = ptr->p;
    (*pl)[k++].type = ptr->type;
  }
  free(gpl);

  (*pl)[k].i      = 0;
  (*pl)[k].j      = 0;
  (*pl)[k].p      = 0.;
  (*pl)[k++].type = 0.;
  free(pt);
  *pl = (plist *)xrealloc(*pl, k * sizeof(plist));
}


/*###########################################*/
/*# deprecated functions below              #*/
/*###########################################*/

PUBLIC int HairpinE(int size, int type, int si1, int sj1, const char *string) {
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

PUBLIC void initialize_fold(int length){
  /* DO NOTHING */
}

