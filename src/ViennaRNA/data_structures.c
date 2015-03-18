/** \file data_structures.c **/

/*
                  Data structure creation/destruction

                  This file contains everything which is necessary to
                  obtain and destroy datastructures used in the folding
                  recurrences throughout the VienneRNA paclage

                  c Ronny Lorenx

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
#include "ViennaRNA/gquad.h"
#include "ViennaRNA/aln_util.h"
#include "ViennaRNA/ribo.h"
#include "ViennaRNA/constraints.h"
#include "ViennaRNA/part_func.h"
#include "ViennaRNA/cofold.h"
#include "ViennaRNA/mm.h"

#ifdef _OPENMP
#include <omp.h>
#endif

/*
#################################
# PRIVATE MACROS                #
#################################
*/

/* the definitions below indicate which arrays should be allocated upon retrieval of a matrices data structure */
#define ALLOC_NOTHING     0
#define ALLOC_F           1
#define ALLOC_F5          2
#define ALLOC_F3          4
#define ALLOC_FC          8
#define ALLOC_C           16
#define ALLOC_FML         32
#define ALLOC_PROBS       256
#define ALLOC_AUX         512

#define ALLOC_CIRC        1024
#define ALLOC_HYBRID      2048
#define ALLOC_UNIQ        4096


#define ALLOC_MFE_DEFAULT         (ALLOC_F5 | ALLOC_C | ALLOC_FML)

#define ALLOC_PF_WO_PROBS         (ALLOC_F | ALLOC_C | ALLOC_FML)
#define ALLOC_PF_DEFAULT          (ALLOC_PF_WO_PROBS | ALLOC_PROBS | ALLOC_AUX)



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
PRIVATE void            add_pf_matrices( vrna_fold_compound *vc, vrna_mx_t type, unsigned int alloc_vector);
PRIVATE void            add_mfe_matrices(vrna_fold_compound *vc, vrna_mx_t type, unsigned int alloc_vector);
PRIVATE void            set_fold_compound(vrna_fold_compound *vc, vrna_md_t *md_p, unsigned int options);
PRIVATE void            make_pscores(vrna_fold_compound *vc);
PRIVATE vrna_mx_mfe_t   *get_mfe_matrices_alloc(unsigned int n, vrna_mx_t type, unsigned int alloc_vector);
PRIVATE vrna_mx_pf_t    *get_pf_matrices_alloc(unsigned int n, vrna_mx_t type, unsigned int alloc_vector);
PRIVATE void            init_dist_class_feature(vrna_fold_compound *vc, const char *s1, const char *s2);


/*
#################################
# BEGIN OF FUNCTION DEFINITIONS #
#################################
*/

PUBLIC void
vrna_free_mfe_matrices(vrna_fold_compound *vc){

  if(vc){
    vrna_mx_mfe_t *self = vc->matrices;
    if(self){
      switch(self->type){
        case VRNA_MX_DEFAULT:   free(self->f5);
                                free(self->f3);
                                free(self->fc);
                                free(self->c);
                                free(self->fML);
                                free(self->fM1);
                                free(self->fM2);
                                free(self->ggg);
                                break;
        case VRNA_MX_2DFOLD:    
                                break;
        default:                /* do nothing */
                                break;
      }
      free(self);
      vc->matrices = NULL;
    }
  }
}

PUBLIC void
vrna_free_pf_matrices(vrna_fold_compound *vc){

  if(vc){
    vrna_mx_pf_t  *self = vc->exp_matrices;
    if(self){
      switch(self->type){
        case VRNA_MX_DEFAULT:   free(self->q);
                                free(self->qb);
                                free(self->qm);
                                free(self->qm1);
                                free(self->qm2);
                                free(self->probs);
                                free(self->G);
                                free(self->q1k);
                                free(self->qln);
                                free(self->scale);
                                free(self->expMLbase);
                                break;
        case VRNA_MX_2DFOLD:    
                                break;
        default:                /* do nothing */
                                break;
      }
      free(self);
      vc->exp_matrices = NULL;
    }
  }
}

PUBLIC  void
vrna_free_fold_compound(vrna_fold_compound *vc){

  int s;

  if(vc){

    /* first destroy common attributes */
    vrna_free_mfe_matrices(vc);
    vrna_free_pf_matrices(vc);
    free(vc->iindx);
    free(vc->jindx);
    free(vc->params);
    free(vc->exp_params);
    vrna_hc_free(vc->hc);

    /* now distinguish the vc type */
    switch(vc->type){
      case VRNA_VC_TYPE_SINGLE:     free(vc->sequence);
                                    free(vc->sequence_encoding);
                                    free(vc->sequence_encoding2);
                                    free(vc->ptype);
                                    free(vc->ptype_pf_compat);
                                    vrna_sc_free(vc->sc);
                                    break;
      case VRNA_VC_TYPE_ALIGNMENT:  for(s=0;s<vc->n_seq;s++){
                                      free(vc->sequences[s]);
                                      free(vc->S[s]);
                                      free(vc->S5[s]);
                                      free(vc->S3[s]);
                                      free(vc->Ss[s]);
                                      free(vc->a2s[s]);
                                    }
                                    free(vc->sequences);
                                    free(vc->cons_seq);
                                    free(vc->S_cons);
                                    free(vc->S);
                                    free(vc->S5);
                                    free(vc->S3);
                                    free(vc->Ss);
                                    free(vc->a2s);
                                    free(vc->pscore);
                                    if(vc->scs){
                                      for(s=0;s<vc->n_seq;s++)
                                        vrna_sc_free(vc->scs[s]);
                                      free(vc->scs);
                                    }
                                    break;
      default:                      /* do nothing */
                                    break;
    }

    free(vc);
  }
}


PUBLIC vrna_fold_compound*
vrna_get_fold_compound( const char *sequence,
                        vrna_md_t *md_p,
                        unsigned int options){

  int length;
  vrna_fold_compound *vc;

  if(sequence == NULL) return NULL;

  /* sanity check */
  length = strlen(sequence);
  if(length == 0)
    vrna_message_error("vrna_get_fold_compound: sequence length must be greater 0");

  vc            = (vrna_fold_compound *)vrna_alloc(sizeof(vrna_fold_compound));
  vc->type      = VRNA_VC_TYPE_SINGLE;
  vc->length    = length;
  vc->sequence  = strdup(sequence);
  set_fold_compound(vc, md_p, options);

  return vc;
}

PUBLIC vrna_fold_compound*
vrna_get_fold_compound_ali( const char **sequences,
                            vrna_md_t *md_p,
                            unsigned int options){

  int s, n_seq, length;
  vrna_fold_compound *vc;
  
  if(sequences == NULL) return NULL;

  for(s=0;sequences[s];s++); /* count the sequences */

  n_seq = s;

  length = strlen(sequences[0]);
  /* sanity check */
  if(length == 0)
    vrna_message_error("vrna_get_fold_compound_ali: sequence length must be greater 0");
  for(s = 0; s < n_seq; s++)
    if(strlen(sequences[s]) != length)
      vrna_message_error("vrna_get_fold_compound_ali: uneqal sequence lengths in alignment");

  vc            = (vrna_fold_compound *)vrna_alloc(sizeof(vrna_fold_compound));
  vc->type      = VRNA_VC_TYPE_ALIGNMENT;

  vc->n_seq     = n_seq;
  vc->length    = length;
  vc->sequences = (char **)vrna_alloc(sizeof(char *) * (vc->n_seq + 1));
  for(s = 0; sequences[s]; s++)
    vc->sequences[s] = strdup(sequences[s]);

  set_fold_compound(vc, md_p, options);

  return vc;
}


/*
#####################################
# BEGIN OF STATIC HELPER FUNCTIONS  #
#####################################
*/

PRIVATE void
set_fold_compound(vrna_fold_compound *vc,
                  vrna_md_t *md_p,
                  unsigned int options){


  char *sequence, **sequences;
  vrna_md_t           md;
  unsigned int        alloc_vector, length, s;
  int                 cp;                     /* cut point for cofold */
  char                *seq, *seq2;

  sequence    = NULL;
  sequences   = NULL;

  /* default values */
  cp            = -1;
  alloc_vector  = ALLOC_NOTHING;

  /* get a copy of the model details */
  if(md_p)
    md = *md_p;
  else /* this fallback relies on global parameters and thus is not threadsafe */
    vrna_md_set_globals(&md);

  switch(vc->type){
    case VRNA_VC_TYPE_SINGLE:     sequence  = vc->sequence;

                                  seq2 = strdup(sequence);
                                  seq = vrna_cut_point_remove(seq2, &cp); /*  splice out the '&' if concatenated sequences and
                                                                        reset cp... this should also be safe for
                                                                        single sequences */
                                  vc->cutpoint            = cp;

                                  if((cp > 0) && (md.min_loop_size == TURN))
                                    md.min_loop_size = 0;  /* is it safe to set this here? */

                                  free(vc->sequence);
                                  vc->sequence            = seq;
                                  vc->length              = length = strlen(seq);
                                  vc->sequence_encoding   = vrna_seq_encode(seq, &md);
                                  vc->sequence_encoding2  = vrna_seq_encode_simple(seq, &md);
                                  if(!(options & VRNA_OPTION_EVAL_ONLY)){
                                    vc->ptype               = vrna_get_ptypes(vc->sequence_encoding2, &md);
                                    vc->ptype_pf_compat     = (options & VRNA_OPTION_PF) ? get_ptypes(vc->sequence_encoding2, &md, 1) : NULL;
                                  } else {
                                    vc->ptype           = NULL;
                                    vc->ptype_pf_compat = NULL;
                                  }
                                  vc->sc                  = NULL;
                                  free(seq2);
                                  break;

    case VRNA_VC_TYPE_ALIGNMENT:  sequences     = vc->sequences;

                                  vc->length    = length = vc->length;

                                  vc->cons_seq  = consensus((const char **)sequences);
                                  vc->S_cons    = vrna_seq_encode_simple(vc->cons_seq, &md);

                                  vc->pscore    = (int *) vrna_alloc(sizeof(int)*((length*(length+1))/2+2));

                                  oldAliEn = vc->oldAliEn  = md.oldAliEn;

                                  vc->S   = (short **)          vrna_alloc((vc->n_seq+1) * sizeof(short *));
                                  vc->S5  = (short **)          vrna_alloc((vc->n_seq+1) * sizeof(short *));
                                  vc->S3  = (short **)          vrna_alloc((vc->n_seq+1) * sizeof(short *));
                                  vc->a2s = (unsigned short **) vrna_alloc((vc->n_seq+1) * sizeof(unsigned short *));
                                  vc->Ss  = (char **)           vrna_alloc((vc->n_seq+1) * sizeof(char *));

                                  for (s = 0; s < vc->n_seq; s++) {
                                    vrna_ali_encode(vc->sequences[s],
                                                    &(vc->S[s]),
                                                    &(vc->S5[s]),
                                                    &(vc->S3[s]),
                                                    &(vc->Ss[s]),
                                                    &(vc->a2s[s]),
                                                    &md);
                                  }
                                  vc->S5[vc->n_seq]  = NULL;
                                  vc->S3[vc->n_seq]  = NULL;
                                  vc->a2s[vc->n_seq] = NULL;
                                  vc->Ss[vc->n_seq]  = NULL;
                                  vc->S[vc->n_seq]   = NULL;

                                  vc->scs       = NULL;
                                  break;

    default:                      /* do nothing ? */
                                  break;
  }

  /* some default init values */
  vc->iindx         = (options & VRNA_OPTION_PF) ? vrna_get_iindx(vc->length) : NULL;
  vc->jindx         = vrna_get_indx(vc->length);

  vc->params        = NULL;
  vc->exp_params    = NULL;
  vc->matrices      = NULL;
  vc->exp_matrices  = NULL;
  vc->hc            = NULL;

  /* now come the energy parameters */
  if(options & VRNA_OPTION_MFE)
    vc->params      = vrna_params_get(&md);

  if(options & VRNA_OPTION_PF){
    vc->exp_params  = (vc->type == VRNA_VC_TYPE_SINGLE) ? \
                        vrna_exp_params_get(&md) : \
                        vrna_exp_params_ali_get(vc->n_seq, &md);
  }

  if(!(options & VRNA_OPTION_EVAL_ONLY)){ /* allocate memory for DP matrices */
    /* extract the matrix type from options flags */
    vrna_mx_t type = VRNA_MX_DEFAULT;

    if(options & VRNA_OPTION_DIST_CLASS)
      type = VRNA_MX_2DFOLD;
    else if(options & VRNA_OPTION_LFOLD)
      type = VRNA_MX_LFOLD;

    /* cofolding matrices ? */
    if(options & VRNA_OPTION_HYBRID){
      alloc_vector |= ALLOC_HYBRID;
      if(cp < 0) /* we could also omit this message */
        vrna_message_warning("vrna_get_fold_compound@data_structures.c: hybrid DP matrices requested but single sequence provided");
    }

    /* default MFE matrices ? */
    if(options & VRNA_OPTION_MFE)
      alloc_vector |= ALLOC_MFE_DEFAULT;

    /* default PF matrices ? */
    if(options & VRNA_OPTION_PF)
      alloc_vector |= (md.compute_bpp) ? ALLOC_PF_DEFAULT : ALLOC_PF_WO_PROBS;

    /* unique ML decomposition ? */
    if(md.uniq_ML)
      alloc_vector |= ALLOC_UNIQ;

    /* matrices for circular folding ? */
    if(md.circ){
      md.uniq_ML = 1;
      alloc_vector |= ALLOC_CIRC | ALLOC_UNIQ;
    }

    /* done with preparations, allocate memory now! */
    if(options & VRNA_OPTION_MFE)
      add_mfe_matrices(vc, type, alloc_vector);

    if(options & VRNA_OPTION_PF)
      add_pf_matrices(vc, type, alloc_vector);
  }

  if(vc->type == VRNA_VC_TYPE_ALIGNMENT)
    make_pscores(vc);

  if(!(options & VRNA_OPTION_EVAL_ONLY))
    vrna_hc_init(vc); /* add default hard constraints */

}


PRIVATE void
init_dist_class_feature(vrna_fold_compound *vc,
                        const char *s1,
                        const char *s2){

  int n     = vc->length;
  int turn  = vc->params->model_details.min_loop_size;

  vc->reference_pt1 = vrna_pt_get(s1);
  vc->reference_pt2 = vrna_pt_get(s2);
  vc->referenceBPs1 = vrna_refBPcnt_matrix(vc->reference_pt1, turn);
  vc->referenceBPs2 = vrna_refBPcnt_matrix(vc->reference_pt2, turn);
  vc->bpdist        = vrna_refBPdist_matrix(vc->reference_pt1, vc->reference_pt2, turn);
  /* compute maximum matching with reference structure 1 disallowed */
  vc->mm1           = maximumMatchingConstraint(vc->sequence, vc->reference_pt1);
  /* compute maximum matching with reference structure 2 disallowed */
  vc->mm2           = maximumMatchingConstraint(vc->sequence, vc->reference_pt2);

  vc->maxD1         = vc->mm1[vc->jindx[1]-n] + vc->referenceBPs1[vc->jindx[1]-n];
  vc->maxD2         = vc->mm2[vc->iindx[1]-n] + vc->referenceBPs2[vc->iindx[1]-n];
}

PRIVATE void
add_pf_matrices(vrna_fold_compound *vc,
                vrna_mx_t type,
                unsigned int alloc_vector){

  if(vc){
    vc->exp_matrices  = get_pf_matrices_alloc(vc->length, type, alloc_vector);
    if(vc->exp_params->model_details.gquad){
      switch(vc->type){
        case VRNA_VC_TYPE_SINGLE:   vc->exp_matrices->G = get_gquad_pf_matrix(vc->sequence_encoding2, vc->exp_matrices->scale, vc->exp_params);
                                    break;
        default:                    /* do nothing */
                                    break;
      }
    }
    vrna_update_pf_params(vc, NULL);
  }
}

PRIVATE void
add_mfe_matrices( vrna_fold_compound *vc,
                  vrna_mx_t type,
                  unsigned int alloc_vector){

  if(vc){
    vc->matrices = get_mfe_matrices_alloc(vc->length, type, alloc_vector);

    if(vc->params->model_details.gquad){
      switch(vc->type){
        case VRNA_VC_TYPE_SINGLE:     vc->matrices->ggg = get_gquad_matrix(vc->sequence_encoding2, vc->params);
                                      break;
        case VRNA_VC_TYPE_ALIGNMENT:  vc->matrices->ggg = get_gquad_ali_matrix(vc->S_cons, vc->S, vc->n_seq,  vc->params);
                                      break;
        default:                      /* do nothing */
                                      break;
      }
    }
  }
}



PRIVATE vrna_mx_mfe_t  *
get_mfe_matrices_alloc( unsigned int n,
                        vrna_mx_t type,
                        unsigned int alloc_vector){

  if(n >= (unsigned int)sqrt((double)INT_MAX))
    vrna_message_error("get_mfe_matrices_alloc@data_structures.c: sequence length exceeds addressable range");

  vrna_mx_mfe_t *vars   = (vrna_mx_mfe_t *)vrna_alloc(sizeof(vrna_mx_mfe_t));

  /* set everything to zero (although this should be done already by vrna_alloc() )*/

  vars->length          = n;
  vars->type            = type;

  unsigned int i, size, lin_size;

  size     = ((n + 1) * (n + 2)) / 2;
  lin_size = n + 2;

  switch(type){
    case VRNA_MX_DEFAULT:   if(alloc_vector & ALLOC_F5)
                              vars->f5  = (int *) vrna_alloc(sizeof(int) * lin_size);
                            else
                              vars->f5  = NULL;

                            if(alloc_vector & ALLOC_F3)
                              vars->f3  = (int *) vrna_alloc(sizeof(int) * lin_size);
                            else
                              vars->f3  = NULL;

                            if(alloc_vector & ALLOC_HYBRID)
                              vars->fc  = (int *) vrna_alloc(sizeof(int) * lin_size);
                            else
                              vars->fc  = NULL;

                            if(alloc_vector & ALLOC_C)
                              vars->c   = (int *) vrna_alloc(sizeof(int) * size);
                            else
                              vars->c   = NULL;

                            if(alloc_vector & ALLOC_FML)
                              vars->fML = (int *) vrna_alloc(sizeof(int) * size);
                            else
                              vars->fML = NULL;

                            if(alloc_vector & ALLOC_UNIQ)
                              vars->fM1 = (int *) vrna_alloc(sizeof(int) * size);
                            else
                              vars->fM1 = NULL;

                            if(alloc_vector & ALLOC_CIRC)
                              vars->fM2 = (int *) vrna_alloc(sizeof(int) * lin_size);
                            else
                              vars->fM2 = NULL;

                            /* setting exterior loop energies for circular case to INF is always safe */
                            vars->FcH = vars->FcI = vars->FcM = vars->Fc = INF;

                            vars->ggg = NULL;

                            break;

    case VRNA_MX_2DFOLD:    if(alloc_vector & ALLOC_F5){
                              vars->E_F5      = (int ***) vrna_alloc(sizeof(int **) * lin_size);
                              vars->l_min_F5  = (int **)  vrna_alloc(sizeof(int *)  * lin_size);
                              vars->l_max_F5  = (int **)  vrna_alloc(sizeof(int *)  * lin_size);
                              vars->k_min_F5  = (int *)   vrna_alloc(sizeof(int)    * lin_size);
                              vars->k_max_F5  = (int *)   vrna_alloc(sizeof(int)    * lin_size);
                              vars->E_F5_rem  = (int *)   vrna_alloc(sizeof(int)    * lin_size);
                              for(i = 0; i <= n; i++)
                                vars->E_F5_rem[i] = INF;
                            } else {
                              vars->E_F5      = NULL;
                              vars->l_min_F5  = NULL;
                              vars->l_max_F5  = NULL;
                              vars->k_min_F5  = NULL;
                              vars->k_max_F5  = NULL;
                              vars->E_F5_rem  = NULL;
                            }

                            if(alloc_vector & ALLOC_F3){
                              vars->E_F3      = (int ***) vrna_alloc(sizeof(int **)  * lin_size);
                              vars->l_min_F3  = (int **)  vrna_alloc(sizeof(int *)   * lin_size);
                              vars->l_max_F3  = (int **)  vrna_alloc(sizeof(int *)   * lin_size);
                              vars->k_min_F3  = (int *)   vrna_alloc(sizeof(int)     * lin_size);
                              vars->k_max_F3  = (int *)   vrna_alloc(sizeof(int)     * lin_size);
                            } else {
                              vars->E_F3      = NULL;
                              vars->l_min_F3  = NULL;
                              vars->l_max_F3  = NULL;
                              vars->k_min_F3  = NULL;
                              vars->k_max_F3  = NULL;
                            }

                            if(alloc_vector & ALLOC_C){
                              vars->E_C     = (int ***) vrna_alloc(sizeof(int **) * size);
                              vars->l_min_C = (int **)  vrna_alloc(sizeof(int *)  * size);
                              vars->l_max_C = (int **)  vrna_alloc(sizeof(int *)  * size);
                              vars->k_min_C = (int *)   vrna_alloc(sizeof(int)    * size);
                              vars->k_max_C = (int *)   vrna_alloc(sizeof(int)    * size);
                              vars->E_C_rem = (int *)   vrna_alloc(sizeof(int)    * size);
                              for(i = 0; i < size; i++)
                                vars->E_C_rem[i] = INF;
                            } else {
                              vars->E_C     = NULL;
                              vars->l_min_C = NULL;
                              vars->l_max_C = NULL;
                              vars->k_min_C = NULL;
                              vars->k_max_C = NULL;
                              vars->E_C_rem = NULL;
                            }

                            if(alloc_vector & ALLOC_FML){
                              vars->E_M     = (int ***) vrna_alloc(sizeof(int **) * size);
                              vars->l_min_M = (int **)  vrna_alloc(sizeof(int *)  * size);
                              vars->l_max_M = (int **)  vrna_alloc(sizeof(int *)  * size);
                              vars->k_min_M = (int *)   vrna_alloc(sizeof(int)    * size);
                              vars->k_max_M = (int *)   vrna_alloc(sizeof(int)    * size);
                              vars->E_M_rem = (int *)   vrna_alloc(sizeof(int)    * size);
                              for(i = 0; i < size; i++)
                                vars->E_M_rem[i] = INF;
                            } else {
                              vars->E_M     = NULL;
                              vars->l_min_M = NULL;
                              vars->l_max_M = NULL;
                              vars->k_min_M = NULL;
                              vars->k_max_M = NULL;
                              vars->E_M_rem = NULL;
                            }

                            if(alloc_vector & ALLOC_UNIQ){
                              vars->E_M1      = (int ***) vrna_alloc(sizeof(int **) * size);
                              vars->l_min_M1  = (int **)  vrna_alloc(sizeof(int *)  * size);
                              vars->l_max_M1  = (int **)  vrna_alloc(sizeof(int *)  * size);
                              vars->k_min_M1  = (int *)   vrna_alloc(sizeof(int)    * size);
                              vars->k_max_M1  = (int *)   vrna_alloc(sizeof(int)    * size);
                              vars->E_M1_rem  = (int *)   vrna_alloc(sizeof(int)    * size);
                              for(i = 0; i < size; i++)
                                vars->E_M1_rem[i] = INF;
                            } else {
                              vars->E_M1      = NULL;
                              vars->l_min_M1  = NULL;
                              vars->l_max_M1  = NULL;
                              vars->k_min_M1  = NULL;
                              vars->k_max_M1  = NULL;
                              vars->E_M1_rem  = NULL;
                            }

                            if(alloc_vector & ALLOC_CIRC){
                              vars->E_M2      = (int ***) vrna_alloc(sizeof(int **)  * lin_size);
                              vars->l_min_M2  = (int **)  vrna_alloc(sizeof(int *)   * lin_size);
                              vars->l_max_M2  = (int **)  vrna_alloc(sizeof(int *)   * lin_size);
                              vars->k_min_M2  = (int *)   vrna_alloc(sizeof(int)     * lin_size);
                              vars->k_max_M2  = (int *)   vrna_alloc(sizeof(int)     * lin_size);
                              vars->E_M2_rem  = (int *)   vrna_alloc(sizeof(int)     * lin_size);
                              for(i = 0; i <= n; i++)
                                vars->E_M2_rem[i] = INF;
                            } else {
                              vars->E_M2      = NULL;
                              vars->l_min_M2  = NULL;
                              vars->l_max_M2  = NULL;
                              vars->k_min_M2  = NULL;
                              vars->k_max_M2  = NULL;
                              vars->E_M2_rem  = NULL;
                            }

                            /* setting exterior loop energies for circular case to INF is always safe */
                            vars->E_Fc      = NULL;
                            vars->E_FcH     = NULL;
                            vars->E_FcI     = NULL;
                            vars->E_FcM     = NULL;
                            vars->E_Fc_rem  = INF;
                            vars->E_FcH_rem = INF;
                            vars->E_FcI_rem = INF;
                            vars->E_FcM_rem = INF;

#ifdef COUNT_STATES
                            vars->N_C   = (unsigned long ***) vrna_alloc(sizeof(unsigned long **)  * size);
                            vars->N_F5  = (unsigned long ***) vrna_alloc(sizeof(unsigned long **)  * lin_size);
                            vars->N_M   = (unsigned long ***) vrna_alloc(sizeof(unsigned long **)  * size);
                            vars->N_M1  = (unsigned long ***) vrna_alloc(sizeof(unsigned long **)  * size);
#endif

                            break;
    default:                /* do nothing */
                            break;
  }

  return vars;
}



PRIVATE vrna_mx_pf_t  *
get_pf_matrices_alloc(unsigned int n,
                      vrna_mx_t type,
                      unsigned int alloc_vector){

  if(n >= (unsigned int)sqrt((double)INT_MAX))
    vrna_message_error("get_pf_matrices_alloc@data_structures.c: sequence length exceeds addressable range");

  vrna_mx_pf_t  *vars   = (vrna_mx_pf_t *)vrna_alloc(sizeof(vrna_mx_pf_t));

  vars->length          = n;
  vars->type            = type;
  vars->q               = NULL;
  vars->qb              = NULL;
  vars->qm              = NULL;
  vars->qm1             = NULL;
  vars->qm2             = NULL;
  vars->qho             = 0.;
  vars->qio             = 0.;
  vars->qmo             = 0.;
  vars->qo              = 0.;
  vars->G               = NULL;
  vars->q1k             = NULL;
  vars->qln             = NULL;

  vars->scale           = NULL;
  vars->expMLbase       = NULL;

  if(alloc_vector){
    unsigned int size     = ((n + 1) * (n + 2)) / 2;
    unsigned int lin_size = n + 2;

    switch(type){
      case VRNA_MX_DEFAULT:   if(alloc_vector & ALLOC_F)
                                vars->q     = (FLT_OR_DBL *) vrna_alloc(sizeof(FLT_OR_DBL) * size);

                              if(alloc_vector & ALLOC_C)
                                vars->qb    = (FLT_OR_DBL *) vrna_alloc(sizeof(FLT_OR_DBL) * size);

                              if(alloc_vector & ALLOC_FML)
                                vars->qm    = (FLT_OR_DBL *) vrna_alloc(sizeof(FLT_OR_DBL) * size);

                              if(alloc_vector & ALLOC_UNIQ)
                                vars->qm1   = (FLT_OR_DBL *) vrna_alloc(sizeof(FLT_OR_DBL) * size);

                              if(alloc_vector & ALLOC_CIRC)
                                vars->qm2   = (FLT_OR_DBL *) vrna_alloc(sizeof(FLT_OR_DBL) * lin_size);

                              if(alloc_vector & ALLOC_PROBS)
                                vars->probs = (FLT_OR_DBL *) vrna_alloc(sizeof(FLT_OR_DBL) * size);

                              if(alloc_vector & ALLOC_AUX){
                                vars->q1k   = (FLT_OR_DBL *) vrna_alloc(sizeof(FLT_OR_DBL) * lin_size);
                                vars->qln   = (FLT_OR_DBL *) vrna_alloc(sizeof(FLT_OR_DBL) * lin_size);
                              }

                              /*  always alloc the helper arrays for unpaired nucleotides in multi-
                                  branch loops and scaling
                              */
                              vars->scale     = (FLT_OR_DBL *) vrna_alloc(sizeof(FLT_OR_DBL) * lin_size);
                              vars->expMLbase = (FLT_OR_DBL *) vrna_alloc(sizeof(FLT_OR_DBL) * lin_size);
                              break;
      default:                /* do nothing */
                              break;
    }
  }

  return vars;
}

PRIVATE void
make_pscores(vrna_fold_compound *vc){

  /* calculate co-variance bonus for each pair depending on  */
  /* compensatory/consistent mutations and incompatible seqs */
  /* should be 0 for conserved pairs, >0 for good pairs      */

#define NONE -10000 /* score for forbidden pairs */

  char *structure = NULL;
  int i,j,k,l,s, max_span;
  float **dm;
  int olddm[7][7]={{0,0,0,0,0,0,0}, /* hamming distance between pairs */
                  {0,0,2,2,1,2,2} /* CG */,
                  {0,2,0,1,2,2,2} /* GC */,
                  {0,2,1,0,2,1,2} /* GU */,
                  {0,1,2,2,0,2,1} /* UG */,
                  {0,2,2,1,2,0,2} /* AU */,
                  {0,2,2,2,1,2,0} /* UA */};

  short           **S         = vc->S;
  char            **AS        = vc->sequences;
  int             n_seq       = vc->n_seq;
  vrna_md_t       *md         = (vc->params) ? &(vc->params->model_details) : &(vc->exp_params->model_details);
  int             *pscore     = vc->pscore;     /* precomputed array of pair types */             
  int             *indx       = vc->jindx;                                             
  int             n           = vc->length;                                            

  if (md->ribo) {
    if (RibosumFile !=NULL) dm=readribosum(RibosumFile);
    else dm=get_ribosum((const char **)AS, n_seq, n);
  }
  else { /*use usual matrix*/
    dm=(float **)vrna_alloc(7*sizeof(float*));
    for (i=0; i<7;i++) {
      dm[i]=(float *)vrna_alloc(7*sizeof(float));
      for (j=0; j<7; j++)
        dm[i][j] = (float) olddm[i][j];
    }
  }

  max_span = md->max_bp_span;
  if((max_span < TURN+2) || (max_span > n))
    max_span = n;
  for (i=1; i<n; i++) {
    for (j=i+1; (j<i+TURN+1) && (j<=n); j++)
      pscore[indx[j]+i] = NONE;
    for (j=i+TURN+1; j<=n; j++) {
      int pfreq[8]={0,0,0,0,0,0,0,0};
      double score;
      for (s=0; s<n_seq; s++) {
        int type;
        if (S[s][i]==0 && S[s][j]==0) type = 7; /* gap-gap  */
        else {
          if ((AS[s][i] == '~')||(AS[s][j] == '~')) type = 7;
          else type = md->pair[S[s][i]][S[s][j]];
        }
        pfreq[type]++;
      }
      if (pfreq[0]*2+pfreq[7]>n_seq) { pscore[indx[j]+i] = NONE; continue;}
      for (k=1,score=0; k<=6; k++) /* ignore pairtype 7 (gap-gap) */
        for (l=k; l<=6; l++)
          score += pfreq[k]*pfreq[l]*dm[k][l];
      /* counter examples score -1, gap-gap scores -0.25   */
      pscore[indx[j]+i] = md->cv_fact *
        ((UNIT*score)/n_seq - md->nc_fact*UNIT*(pfreq[0] + pfreq[7]*0.25));

      if((j - i + 1) > max_span){
        pscore[indx[j]+i] = NONE;
      }
    }
  }

  if (md->noLP) /* remove unwanted pairs */
    for (k=1; k<n-TURN-1; k++)
      for (l=1; l<=2; l++) {
        int type,ntype=0,otype=0;
        i=k; j = i+TURN+l;
        type = pscore[indx[j]+i];
        while ((i>=1)&&(j<=n)) {
          if ((i>1)&&(j<n)) ntype = pscore[indx[j+1]+i-1];
          if ((otype<md->cv_fact*MINPSCORE)&&(ntype<md->cv_fact*MINPSCORE))  /* too many counterexamples */
            pscore[indx[j]+i] = NONE; /* i.j can only form isolated pairs */
          otype =  type;
          type  = ntype;
          i--; j++;
        }
      }


  if (fold_constrained&&(structure!=NULL)) {
    int psij, hx, hx2, *stack, *stack2;
    stack = (int *) vrna_alloc(sizeof(int)*(n+1));
    stack2 = (int *) vrna_alloc(sizeof(int)*(n+1));

    for(hx=hx2=0, j=1; j<=n; j++) {
      switch (structure[j-1]) {
      case 'x': /* can't pair */
        for (l=1; l<j-TURN; l++) pscore[indx[j]+l] = NONE;
        for (l=j+TURN+1; l<=n; l++) pscore[indx[l]+j] = NONE;
        break;
      case '(':
        stack[hx++]=j;
        /* fallthrough */
      case '[':
        stack2[hx2++]=j;
        /* fallthrough */
      case '<': /* pairs upstream */
        for (l=1; l<j-TURN; l++) pscore[indx[j]+l] = NONE;
        break;
      case ']':
        if (hx2<=0) {
          fprintf(stderr, "%s\n", structure);
          vrna_message_error("unbalanced brackets in constraints");
        }
        i = stack2[--hx2];
        pscore[indx[j]+i]=NONE;
        break;
      case ')':
        if (hx<=0) {
          fprintf(stderr, "%s\n", structure);
          vrna_message_error("unbalanced brackets in constraints");
        }
        i = stack[--hx];
        psij = pscore[indx[j]+i]; /* store for later */
        for (k=j; k<=n; k++)
          for (l=i; l<=j; l++)
            pscore[indx[k]+l] = NONE;
        for (l=i; l<=j; l++)
          for (k=1; k<=i; k++)
            pscore[indx[l]+k] = NONE;
        for (k=i+1; k<j; k++)
          pscore[indx[k]+i] = pscore[indx[j]+k] = NONE;
        pscore[indx[j]+i] = (psij>0) ? psij : 0;
        /* fallthrough */
      case '>': /* pairs downstream */
        for (l=j+TURN+1; l<=n; l++) pscore[indx[l]+j] = NONE;
        break;
      }
    }
    if (hx!=0) {
      fprintf(stderr, "%s\n", structure);
      vrna_message_error("unbalanced brackets in constraint string");
    }
    free(stack); free(stack2);
  }
  /*free dm */
  for (i=0; i<7;i++) {
    free(dm[i]);
  }
  free(dm);
}
