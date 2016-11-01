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
#include "ViennaRNA/loop_energies.h"
#include "ViennaRNA/fold.h"

#ifdef  VRNA_BACKWARD_COMPAT

#ifdef _OPENMP
#include <omp.h>
#endif

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

#ifdef  VRNA_BACKWARD_COMPAT

/* some backward compatibility stuff */
PRIVATE int                 backward_compat           = 0;
PRIVATE vrna_fold_compound_t  *backward_compat_compound = NULL;

#ifdef _OPENMP

#pragma omp threadprivate(backward_compat_compound, backward_compat)

#endif

#endif

/*
#################################
# PRIVATE FUNCTION DECLARATIONS #
#################################
*/

#ifdef  VRNA_BACKWARD_COMPAT

/* wrappers for old API compatibility */
PRIVATE float wrap_fold( const char *string, char *structure, vrna_param_t *parameters, int is_constrained, int is_circular);
PRIVATE void  wrap_array_export(int **f5_p, int **c_p, int **fML_p, int **fM1_p, int **indx_p, char **ptype_p);
PRIVATE void  wrap_array_export_circ( int *Fc_p, int *FcH_p, int *FcI_p, int *FcM_p, int **fM2_p);

#endif

/*
#################################
# BEGIN OF FUNCTION DEFINITIONS #
#################################
*/

PUBLIC float
vrna_fold(const char *string,
          char *structure){

  float                 mfe;
  vrna_fold_compound_t  *vc;
  vrna_md_t             md;

  vrna_md_set_default(&md);
  vc  = vrna_fold_compound(string, &md, 0);
  mfe = vrna_mfe(vc, structure);

  vrna_fold_compound_free(vc);

  return mfe;
}

PUBLIC float
vrna_circfold(const char *string,
              char *structure){

  float                 mfe;
  vrna_fold_compound_t  *vc;
  vrna_md_t             md;

  vrna_md_set_default(&md);
  md.circ = 1;
  vc      = vrna_fold_compound(string, &md, 0);
  mfe     = vrna_mfe(vc, structure);

  vrna_fold_compound_free(vc);

  return mfe;

}

/*###########################################*/
/*# deprecated functions below              #*/
/*###########################################*/

#ifdef  VRNA_BACKWARD_COMPAT

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

  vrna_fold_compound_t  *vc;
  vrna_param_t          *P;
  float                 mfe;

#ifdef _OPENMP
/* Explicitly turn off dynamic threads */
  omp_set_dynamic(0);
#endif

  /* we need the parameter structure for hard constraints */
  if(parameters){
    P = vrna_params_copy(parameters);
  } else {
    vrna_md_t md;
    set_model_details(&md);
    md.temperature = temperature;
    P = vrna_params(&md);
  }
  P->model_details.circ = is_circular;

  vc = vrna_fold_compound(string, &(P->model_details), VRNA_OPTION_DEFAULT);

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

    vrna_constraints_add(vc, (const char *)structure, constraint_options);
  }

  if(backward_compat_compound && backward_compat)
    vrna_fold_compound_free(backward_compat_compound);

  backward_compat_compound  = vc;
  backward_compat           = 1;

  /* call mfe() function without backtracking */
  mfe = vrna_mfe(vc, NULL);

  /* backtrack structure */
  if(structure && vc->params->model_details.backtrack){
    char            *ss;
    int             length;
    sect            bt_stack[MAXSECTORS];
    vrna_bp_stack_t *bp;

    length  = vc->length;
    bp      = (vrna_bp_stack_t *)vrna_alloc(sizeof(vrna_bp_stack_t) * (4*(1+length/2))); /* add a guess of how many G's may be involved in a G quadruplex */

    vrna_backtrack_from_intervals(vc, bp, bt_stack, 0);

    ss = vrna_db_from_bp_stack(bp, length);
    strncpy(structure, ss, length + 1);
    free(ss);

    if(base_pair)
      free(base_pair);
    base_pair = bp;
  }

  return mfe;
}

PUBLIC void
free_arrays(void){

  if(backward_compat_compound && backward_compat){
    vrna_fold_compound_free(backward_compat_compound);
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

  vrna_md_t           md;

  if(backward_compat_compound && backward_compat){
    set_model_details(&md);
    vrna_params_reset(backward_compat_compound, &md);
  }
}

PUBLIC void
update_fold_params_par(vrna_param_t *parameters){

  vrna_md_t           md;

  if(backward_compat_compound && backward_compat){
    if(parameters)
      vrna_params_subst(backward_compat_compound, parameters);
    else{
      set_model_details(&md);
      vrna_params_reset(backward_compat_compound, &md);
    }
  }
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
  vrna_bp_stack_t         *bp         = NULL;
  sect          bt_stack[MAXSECTORS]; /* stack of partial structures for backtracking */

  if(sequence){
    length = strlen(sequence);
    bp = (vrna_bp_stack_t *)vrna_alloc(sizeof(vrna_bp_stack_t) * (1+length/2));
  } else {
    vrna_message_error("backtrack_fold_from_pair@fold.c: no sequence given");
  }

  bt_stack[1].i  = i;
  bt_stack[1].j  = j;
  bt_stack[1].ml = 2;

  bp[0].i = 0; /* ??? this is set by backtrack anyway... */

  vrna_backtrack_from_intervals(backward_compat_compound, bp, bt_stack, 1);
  structure = vrna_db_from_bp_stack(bp, length);

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
