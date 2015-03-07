/* constraints handling */

#include <assert.h>
#include <config.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <string.h>
#include <limits.h>

#include "ViennaRNA/energy_par.h"
#include "ViennaRNA/energy_const.h" /* defines MINPSCORE */
#include "ViennaRNA/fold_vars.h"
#include "ViennaRNA/utils.h"
#include "ViennaRNA/aln_util.h"
#include "ViennaRNA/file_formats.h"
#include "ViennaRNA/params.h"
#include "ViennaRNA/constraints.h"


#ifdef __GNUC__
# define INLINE inline
#else
# define INLINE
#endif

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
PRIVATE INLINE  void
hc_cant_pair( unsigned int i,
              char c_option,
              char *hc,
              unsigned int length,
              unsigned int min_loop_size,
              int *index);

PRIVATE INLINE  void
hc_must_pair( unsigned int i,
              char c_option,
              char *hc,
              int *index);

PRIVATE INLINE  void
hc_pairs_upstream(unsigned int i,
                  char c_option,
                  char *hc,
                  unsigned int min_loop_size,
                  int *index);

PRIVATE INLINE  void
hc_pairs_downstream(unsigned int i,
                    char c_option,
                    char *hc,
                    unsigned int length,
                    unsigned int min_loop_size,
                    int *index);

PRIVATE INLINE  void
hc_allow_pair(unsigned int i,
              unsigned int j,
              char c_option,
              char *hc,
              int *index);

PRIVATE INLINE  void
hc_weak_enforce_pair( unsigned int i,
                      unsigned int j,
                      char c_option,
                      char *hc,
                      unsigned int length,
                      unsigned int min_loop_size,
                      int *index);

PRIVATE INLINE  void
hc_enforce_pair(unsigned int i,
                unsigned int j,
                char c_option,
                char *hc,
                unsigned int length,
                unsigned int min_loop_size,
                int *index);

PRIVATE INLINE  void
hc_intramolecular_only( unsigned int i,
                        char c_option,
                        char *hc,
                        unsigned int length,
                        unsigned int min_loop_size,
                        int cut,
                        int *index);

PRIVATE INLINE  void
hc_intermolecular_only( unsigned int i,
                        char c_option,
                        char *hc,
                        unsigned int length,
                        unsigned int min_loop_size,
                        int cut,
                        int *index);

PRIVATE INLINE  void
adjust_ptypes(char *ptype,
              vrna_hc_t *hc,
              unsigned int length,
              unsigned int indx_type);

PRIVATE void
apply_DB_constraint(const char *constraint,
                    char *ptype,
                    unsigned int length,
                    unsigned int min_loop_size,
                    int cut,
                    unsigned int options);

PRIVATE void
hc_reset_to_default(vrna_fold_compound *vc);

PRIVATE void
hc_update_up(vrna_fold_compound *vc);

PRIVATE void
sc_parse_parameters(const char *string,
                    char c1,
                    char c2,
                    float *v1,
                    float *v2);

PRIVATE void
sc_add_up_mfe(vrna_fold_compound *vc,
              const double *constraints,
              unsigned int options);

PRIVATE void
sc_add_up_pf( vrna_fold_compound *vc,
              const double *constraints,
              unsigned int options);

PRIVATE void
sc_add_bp_mfe(vrna_fold_compound *vc,
              const double **constraints,
              unsigned int options);

PRIVATE void
sc_add_bp_pf( vrna_fold_compound *vc,
              const double **constraints,
              unsigned int options);

PRIVATE void
sc_add_stack_en_mfe(vrna_fold_compound *vc,
                    const double *constraints,
                    unsigned int options);

PRIVATE void
sc_add_stack_en_pf( vrna_fold_compound *vc,
                    const double *constraints,
                    unsigned int options);


/*
#################################
# BEGIN OF FUNCTION DEFINITIONS #
#################################
*/
PUBLIC  void
vrna_message_constraint_options_all(void){

  vrna_message_constraint_options(  VRNA_CONSTRAINT_PIPE
                                  | VRNA_CONSTRAINT_DOT
                                  | VRNA_CONSTRAINT_X
                                  | VRNA_CONSTRAINT_ANG_BRACK
                                  | VRNA_CONSTRAINT_RND_BRACK);
}

PUBLIC  void
vrna_message_constraint_options(unsigned int option){

  if(!(option & VRNA_CONSTRAINT_NO_HEADER)) printf("Input structure constraints using the following notation:\n");
  if(option & VRNA_CONSTRAINT_PIPE)       printf("| : paired with another base\n");
  if(option & VRNA_CONSTRAINT_DOT)        printf(". : no constraint at all\n");
  if(option & VRNA_CONSTRAINT_X)          printf("x : base must not pair\n");
  if(option & VRNA_CONSTRAINT_ANG_BRACK)  printf("< : base i is paired with a base j<i\n> : base i is paired with a base j>i\n");
  if(option & VRNA_CONSTRAINT_RND_BRACK)  printf("matching brackets ( ): base i pairs base j\n");
}

PUBLIC  void
vrna_add_constraints( vrna_fold_compound *vc,
                      const char *constraint,
                      unsigned int options){

  vrna_md_t         *md;

  if(vc){
    if(vc->params)
      md = &(vc->params->model_details);
    else if(vc->exp_params)
      md = &(vc->exp_params->model_details);
    else
      vrna_message_error("constraints.c@vrna_add_constraints: fold compound has no params or exp_params");

    if(!vc->hc)
      vrna_hc_init(vc);

    if(options & VRNA_CONSTRAINT_DB){ /* apply hard constraints from dot-bracket notation */
      apply_DB_constraint(constraint,
                          vc->hc->matrix,
                          vc->length,
                          md->min_loop_size,
                          -1,
                          options);
      hc_update_up(vc);
    } else if(options & VRNA_CONSTRAINT_FILE){ /* constraints from file */
      (void) vrna_read_constraints_file(constraint, vc->length, options);

      /* now do something with the constraints we've just read */

      /* ############################### */
      /* init empty soft constraints     */
      /* ############################### */
      vrna_sc_init(vc);
    }
  }
}

PUBLIC  void
vrna_hc_init(vrna_fold_compound *vc){

  unsigned int  n;
  vrna_hc_t     *hc;

  n           = vc->length;

  /* free previous hard constraints */
  vrna_hc_free(vc->hc);

  /* allocate memory new hard constraints data structure */
  hc          = (vrna_hc_t *)vrna_alloc(sizeof(vrna_hc_t));
  hc->matrix  = (char *)vrna_alloc(sizeof(char)*((n*(n+1))/2+2));
  hc->up_ext  = (int *)vrna_alloc(sizeof(int)*(n+2));
  hc->up_hp   = (int *)vrna_alloc(sizeof(int)*(n+2));
  hc->up_int  = (int *)vrna_alloc(sizeof(int)*(n+2));
  hc->up_ml   = (int *)vrna_alloc(sizeof(int)*(n+2));

  /* set new hard constraints */
  vc->hc = hc;

  /* prefill default values  */
  hc_reset_to_default(vc);

  /* update */
  hc_update_up(vc);
}

PUBLIC void
vrna_hc_add_up( vrna_fold_compound *vc,
                int i,
                char option){

  int j;

  if(vc)
    if(vc->hc){
      if((i <= 0) || (i > vc->length)){
        vrna_message_warning("vrna_hc_add_up: position out of range, not doing anything");
        return;
      }

      vc->hc->matrix[vc->jindx[i] + i] = option & (   VRNA_HC_CONTEXT_EXT_LOOP
                                                    | VRNA_HC_CONTEXT_HP_LOOP
                                                    | VRNA_HC_CONTEXT_INT_LOOP
                                                    | VRNA_HC_CONTEXT_MB_LOOP);

      if(option & VRNA_HC_CONTEXT_ENFORCE){
        /* do not allow i to be paired with any other nucleotide */
        for(j = 1; j < i; j++)
          vc->hc->matrix[vc->jindx[i] + j] = (char)0;
        for(j = i+1; j <= vc->length; j++)
          vc->hc->matrix[vc->jindx[j] + i] = (char)0;
      }

      hc_update_up(vc);
    }
}

PUBLIC void
vrna_hc_add_bp( vrna_fold_compound *vc,
                int i,
                int j,
                char option){

  int p,q, k, l;

  if(vc)
    if(vc->hc){
      if((i <= 0) || (i > vc->length) || (j <= 0) || (j > vc->length)){
        vrna_message_warning("vrna_hc_add_bp: position out of range, not doing anything");
        return;
      }

      if(i == j){
        vrna_message_warning("vrna_hc_add_bp: positions do not differ, not doing anything");
        return;
      }

      if(i < j){
        p = i;
        q = j;
      } else {
        p = j;
        q = i;
      }

      vc->hc->matrix[vc->jindx[q] + p] = option & VRNA_HC_CONTEXT_ALL_LOOPS;
      if(vc->hc->matrix[vc->jindx[q] + p])
        if(vc->ptype[vc->jindx[q] + p] == 0)
          vc->ptype[vc->jindx[q] + p] = 7;

      if(option & VRNA_HC_CONTEXT_ENFORCE){

        /* do not allow i,j to pair with any other nucleotide k */
        for(k = 1; k < p; k++){
          vc->hc->matrix[vc->jindx[p] + k] = (char)0;
          vc->hc->matrix[vc->jindx[q] + k] = (char)0;
          for(l = p+1; l < q; l++)
            vc->hc->matrix[vc->jindx[l] + k] = (char)0;
        }
        for(k = p+1; k < q; k++){
          vc->hc->matrix[vc->jindx[k] + p] = (char)0;
          vc->hc->matrix[vc->jindx[q] + k] = (char)0;
          for(l = q + 1; l <= vc->length; l++)
            vc->hc->matrix[vc->jindx[l] + k] = (char)0;
        }
        for(k = q+1; k <= vc->length; k++){
          vc->hc->matrix[vc->jindx[k] + p] = (char)0;
          vc->hc->matrix[vc->jindx[k] + q] = (char)0;
        }

        /* do not allow i,j to be unpaired */
        vc->hc->matrix[vc->jindx[p] + p] = (char)0;
        vc->hc->matrix[vc->jindx[q] + q] = (char)0;

        hc_update_up(vc);
      }
    }
}

PUBLIC void
vrna_hc_free(vrna_hc_t *hc){

  if(hc){
    if(hc->matrix)  free(hc->matrix);
    if(hc->up_ext)  free(hc->up_ext);
    if(hc->up_hp)   free(hc->up_hp);
    if(hc->up_int)  free(hc->up_int);
    if(hc->up_ml)   free(hc->up_ml);
    free(hc);
  }
}

PRIVATE void
apply_DB_constraint(const char *constraint,
                    char *hc,
                    unsigned int length,
                    unsigned int min_loop_size,
                    int cut,
                    unsigned int options){

  int n,i,j;
  int hx, *stack;
  int *index;
  char c_option;

  if(constraint == NULL) return;

  n         = (int)strlen(constraint);
  stack     = (int *) vrna_alloc(sizeof(int)*(n+1));
  index     = vrna_get_indx(length);
  c_option  =   VRNA_HC_CONTEXT_EXT_LOOP
              | VRNA_HC_CONTEXT_HP_LOOP
              | VRNA_HC_CONTEXT_INT_LOOP
              | VRNA_HC_CONTEXT_INT_LOOP_ENC
              | VRNA_HC_CONTEXT_MB_LOOP
              | VRNA_HC_CONTEXT_MB_LOOP_ENC;

  for(hx=0, j=1; j<=n; j++) {
    switch (constraint[j-1]) {
       /* can't pair */
       case 'x':  if(options & VRNA_CONSTRAINT_X){
                    hc_cant_pair(j, c_option, hc, length, min_loop_size, index);
                  }
                  break;

      /* must pair, i.e. may not be unpaired */
      case '|':   if(options & VRNA_CONSTRAINT_PIPE){
                    hc_must_pair(j, c_option, hc, index);
                  }
                  break;

      /* weak enforced pair 'open' */
      case '(':   if(options & VRNA_CONSTRAINT_RND_BRACK){
                    stack[hx++]=j;
                  }
                  break;

      /* weak enforced pair 'close' */
      case ')':   if(options & VRNA_CONSTRAINT_RND_BRACK){
                    if (hx<=0) {
                      fprintf(stderr, "%s\n", constraint);
                      vrna_message_error("unbalanced brackets in constraints");
                    }
                    i = stack[--hx];
                    hc_weak_enforce_pair(i, j, c_option, hc, length, min_loop_size, index);
                  }
                  break;

      /* pairs upstream */
      case '<':   if(options & VRNA_CONSTRAINT_ANG_BRACK){
                    hc_pairs_upstream(j, c_option, hc, min_loop_size, index);
                  }
                  break;

      /* pairs downstream */
      case '>':   if(options & VRNA_CONSTRAINT_ANG_BRACK){
                    hc_pairs_downstream(j, c_option, hc, length, min_loop_size, index);
                  }
                  break;

      /* only intramolecular basepairing */
      case 'l':   if(options & VRNA_CONSTRAINT_INTRAMOLECULAR){
                    hc_intramolecular_only(j, c_option, hc, length, min_loop_size, cut, index);
                  }
                  break;

      /* only intermolecular bp */
      case 'e':   if(options & VRNA_CONSTRAINT_INTERMOLECULAR){
                    hc_intermolecular_only(j, c_option, hc, length, min_loop_size, cut, index);
                  }
                  break;

      case '.':   break;

      default:    vrna_message_warning("unrecognized character in pseudo dot-bracket notation constraint string\n");
                  break;
    }
  }

  if (hx!=0) {
    fprintf(stderr, "%s\n", constraint);
    vrna_message_error("unbalanced brackets in constraint string");
  }
  /* clean up */
  free(index);
  free(stack);
}

PRIVATE INLINE  void
hc_intramolecular_only( unsigned int i,
                        char c_option,
                        char *hc,
                        unsigned int length,
                        unsigned int min_loop_size,
                        int cut,
                        int *index){

  unsigned int l;

  if(cut > 1){
    if(i < cut)
      for(l = MAX2(i+min_loop_size, cut); l <= length; l++)
        hc[index[l] + i] &= ~c_option;
    else
      for(l = 1; l < MIN2(cut, i-min_loop_size); l++)
        hc[index[i] + l] &= ~c_option;
  }
}

PRIVATE INLINE  void
hc_intermolecular_only( unsigned int i,
                        char c_option,
                        char *hc,
                        unsigned int length,
                        unsigned int min_loop_size,
                        int cut,
                        int *index){

  unsigned int l;

  if(cut > 1){
    if(i < cut){
      for(l = 1; l < i; l++)
        hc[index[i] + l] &= ~c_option;
      for(l = i + 1; l < cut; l++)
        hc[index[l] + i] &= ~c_option;
    } else {
      for(l = cut; l < i; l++)
        hc[index[i] + l] &= ~c_option;
      for(l = i + 1; l <= length; l++)
        hc[index[l] + i] &= ~c_option;
    }
  }
}

PRIVATE INLINE  void
hc_cant_pair( unsigned int i,
              char c_option,
              char *hc,
              unsigned int length,
              unsigned int min_loop_size,
              int *index){

  hc_pairs_upstream(i, c_option, hc, min_loop_size, index);
  hc_pairs_downstream(i, c_option, hc, length, min_loop_size, index);
}

PRIVATE INLINE  void
hc_must_pair( unsigned int i,
              char c_option,
              char *hc,
              int *index){

  hc[index[i]+i] &= ~c_option;
}

PRIVATE INLINE  void
hc_pairs_upstream(unsigned int i,
                  char c_option,
                  char *hc,
                  unsigned int min_loop_size,
                  int *index){

  unsigned int l;

  if(min_loop_size < i){
    for(l = 1; l < i - min_loop_size; l++){
      hc[index[i] + l] &= ~c_option;
    }
  }
}

PRIVATE INLINE  void
hc_pairs_downstream(unsigned int i,
                    char c_option,
                    char *hc,
                    unsigned int length,
                    unsigned int min_loop_size,
                    int *index){

  unsigned int l;

  for(l = i + min_loop_size + 1; l <= length; l++){
    hc[index[l] + i] &= ~c_option;
  }
}

PRIVATE INLINE  void
hc_allow_pair(unsigned int i,
              unsigned int j,
              char c_option,
              char *hc,
              int *index){

  hc[index[j] + i] |= c_option;
}

PRIVATE INLINE  void
hc_weak_enforce_pair( unsigned int i,
                      unsigned int j,
                      char c_option,
                      char *hc,
                      unsigned int length,
                      unsigned int min_loop_size,
                      int *index){

  unsigned int k, l;

  /* don't allow pairs (k,i) 1 <= k < i */
  hc_pairs_upstream(i, c_option, hc, min_loop_size, index);
  /* don't allow pairs (i,k) i < k <= n */ 
  hc_pairs_downstream(i, c_option, hc, length, min_loop_size, index);
  /* don't allow pairs (k,j) 1 <= k < j */
  hc_pairs_upstream(j, c_option, hc, min_loop_size, index);
  /* don't allow pairs (j,k) j < k <= n */ 
  hc_pairs_downstream(j, c_option, hc, length, min_loop_size, index);

  /* don't allow pairs i < k < j < l */
  for(k = i+1; k < j; k++)
    for(l = j+1; l <= length; l++){
      hc[index[l] + k] = 0;
    }
  /* don't allow pairs k<i<l<j */
  for(k = 1; k < i; k++)
    for(l = i+1; l < j; l++){
      hc[index[l] + k] = 0;
    }
  /* allow base pair (i,j) */
  hc[index[j] + i] |= c_option;
}

PRIVATE INLINE  void
hc_enforce_pair(unsigned int i,
                unsigned int j,
                char c_option,
                char *hc,
                unsigned int length,
                unsigned int min_loop_size,
                int *index){

  hc_weak_enforce_pair( i,
                        j,
                        c_option,
                        hc,
                        length,
                        min_loop_size,
                        index);

  /* forbid i and j to be unpaired */
  hc[index[i] + i] = 0;
  hc[index[j] + j] = 0;
}

PRIVATE void
hc_reset_to_default(vrna_fold_compound *vc){

  unsigned int      i, j, ij, min_loop_size, n;
  int               max_span, *idx;
  vrna_md_t         *md;
  vrna_hc_t         *hc;
  short             *S;

  n   = vc->length;
  hc  = vc->hc;
  idx = vc->jindx;
  S   = vc->sequence_encoding;

  if(vc->params)
    md  = &(vc->params->model_details);
  else if(vc->exp_params)
    md  = &(vc->exp_params->model_details);
  else
    vrna_message_error("missing model_details in fold_compound");

  min_loop_size = md->min_loop_size;
  max_span      = md->max_bp_span;

  if((max_span < 5) || (max_span > n))
    max_span  = n;

  /* ######################### */
  /* fill with default values  */
  /* ######################### */

  /* 1. unpaired nucleotides are allowed in all contexts */
  for(i = 1; i <= n; i++)
    hc->matrix[idx[i] + i]  =   VRNA_HC_CONTEXT_EXT_LOOP
                              | VRNA_HC_CONTEXT_HP_LOOP
                              | VRNA_HC_CONTEXT_INT_LOOP
                              | VRNA_HC_CONTEXT_MB_LOOP;

  if(vc->type == VRNA_VC_TYPE_ALIGNMENT){
    /* 2. all base pairs with pscore above threshold are allowed in all contexts */
    for(j = n; j > min_loop_size + 1; j--){
      ij = idx[j]+1;
      for(i=1; i < j - min_loop_size; i++, ij++)
        if((j-i+1) > max_span){
          hc->matrix[ij] = (char)0;
        } else {
          hc->matrix[ij] = (vc->pscore[idx[j]+i] >= md->cv_fact*MINPSCORE) ? VRNA_HC_CONTEXT_ALL_LOOPS : (char)0;
        }
    }
    /* correct for no lonely pairs (assuming that ptypes already incorporate noLP status) */
    /* this should be included in the pscore which is checked above */
  } else {
    /* 2. all canonical base pairs are allowed in all contexts */
    for(j = n; j > min_loop_size + 1; j--){
      ij = idx[j]+1;
      for(i=1; i < j - min_loop_size; i++, ij++){
        char opt = (char)0;
        if((j-i+1) <= max_span){
          int t = md->pair[S[i]][S[j]];
          switch(t){
            case 0:   break;
            case 3:   /* fallthrough */
            case 4:   if(md->noGU){
                        break;
                      } else if(md->noGUclosure){
                        opt = VRNA_HC_CONTEXT_ALL_LOOPS & ~(VRNA_HC_CONTEXT_HP_LOOP | VRNA_HC_CONTEXT_MB_LOOP);
                        break;
                      } /* else fallthrough */
            default:  opt = VRNA_HC_CONTEXT_ALL_LOOPS;
                      break;
          }
        }
        hc->matrix[ij] = opt;
      }
    }

    /* correct for no lonely pairs (assuming that ptypes already incorporate noLP status) */
    /* this should be fixed such that ij loses its hard constraint type if it does not
       allow for enclosing an interior loop, etc.
    */
    if(md->noLP)
      for(i = 1; i < n; i++)
        for(j = i + min_loop_size + 1; j <= n; j++){
          if(hc->matrix[idx[j] +i]){
            if(!vc->ptype[idx[j] + i]){
              hc->matrix[idx[j] + i] = (char)0;
            }
          }
        }
  }
}

PRIVATE void
hc_update_up(vrna_fold_compound *vc){

  unsigned int      i, n;
  int               *idx;
  vrna_hc_t         *hc;

  n   = vc->length;
  idx = vc->jindx;
  hc  = vc->hc;

  for(hc->up_ext[n+1] = 0, i = n; i > 0; i--) /* unpaired stretch in exterior loop */
    hc->up_ext[i] = (hc->matrix[idx[i]+i] & VRNA_HC_CONTEXT_EXT_LOOP) ? 1 + hc->up_ext[i+1] : 0;

  for(hc->up_hp[n+1] = 0, i = n; i > 0; i--)  /* unpaired stretch in hairpin loop */
    hc->up_hp[i] = (hc->matrix[idx[i]+i] & VRNA_HC_CONTEXT_HP_LOOP) ? 1 + hc->up_hp[i+1] : 0;

  for(hc->up_int[n+1] = 0, i = n; i > 0; i--) /* unpaired stretch in interior loop */
    hc->up_int[i] = (hc->matrix[idx[i]+i] & VRNA_HC_CONTEXT_INT_LOOP) ? 1 + hc->up_int[i+1] : 0;

  for(hc->up_ml[n+1] = 0, i = n; i > 0; i--)  /* unpaired stretch in multibranch loop */
    hc->up_ml[i] = (hc->matrix[idx[i]+i] & VRNA_HC_CONTEXT_MB_LOOP) ? 1 + hc->up_ml[i+1] : 0;

  /*
   *  loop arround once more until we find a nucleotide that mustn't
   *  be unpaired (needed for circular folding)
   */

  if(hc->matrix[idx[1]+1] & VRNA_HC_CONTEXT_EXT_LOOP){
    hc->up_ext[n+1] = hc->up_ext[1];
    for(i = n; i > 0; i--){
      if(hc->matrix[idx[i]+i] & VRNA_HC_CONTEXT_EXT_LOOP){
        hc->up_ext[i] = MIN2(n, 1 + hc->up_ext[i+1]);
      } else
        break;
    }
  }

  if(hc->matrix[idx[1]+1] & VRNA_HC_CONTEXT_HP_LOOP){
    hc->up_hp[n+1] = hc->up_hp[1];
    for(i = n; i > 0; i--){
      if(hc->matrix[idx[i]+i] & VRNA_HC_CONTEXT_HP_LOOP){
        hc->up_hp[i] = MIN2(n, 1 + hc->up_hp[i+1]);
      } else
        break;
    }
  }

  if(hc->matrix[idx[1]+1] & VRNA_HC_CONTEXT_INT_LOOP){
    hc->up_int[n+1] = hc->up_int[1];
    for(i = n; i > 0; i--){
      if(hc->matrix[idx[i]+i] & VRNA_HC_CONTEXT_INT_LOOP){
        hc->up_int[i] = MIN2(n, 1 + hc->up_int[i+1]);
      } else
        break;
    }
  }

  if(hc->matrix[idx[1]+1] & VRNA_HC_CONTEXT_MB_LOOP){
    hc->up_ml[n+1] = hc->up_ml[1];
    for(i = n; i > 0; i--){
      if(hc->matrix[idx[i]+i] & VRNA_HC_CONTEXT_MB_LOOP){
        hc->up_ml[i] = MIN2(n, 1 + hc->up_ml[i+1]);
      } else
        break;
    }
  }

}


/* And now, for something completely different...
 *
 * ... the soft constraints section follows below
 */

PUBLIC int
vrna_sc_SHAPE_to_pr(const char *shape_conversion,
                    double *values,
                    int length,
                    double default_value){

  int *indices;
  int i, j;
  int index;
  int ret = 1;

  if(!shape_conversion || !(*shape_conversion) || length <= 0)
    return 0;

  if(*shape_conversion == 'S')
    return 1;

  indices = vrna_alloc(sizeof(int) * (length + 1));
  for (i = 1, j = 0; i <= length; ++i){
    if(values[i] < 0)
      values[i] = default_value;
    else
      indices[j++] = i;
  }

  if(*shape_conversion == 'M'){
    double max;
    double map_info[4][2] = {{0.25, 0.35},
                           {0.30, 0.55},
                           {0.70, 0.85},
                           {0, 1}};

    max = values[1];
    for(i = 2; i <= length; ++i)
      max = MAX2(max, values[i]);
    map_info[3][0] = max;

    for(i = 0; indices[i]; ++i){
      double lower_source = 0;
      double lower_target = 0;

      index = indices[i];

      if(values[index] == 0)
        continue;

      for(j = 0; j < 4; ++j){
        if(values[index] > lower_source && values[index] <= map_info[j][0]){
          double diff_source = map_info[j][0] - lower_source;
          double diff_target = map_info[j][1] - lower_target;
          values[index] = (values[index] - lower_source) / diff_source * diff_target + lower_target;
          break;
        }

        lower_source = map_info[j][0];
        lower_target = map_info[j][1];
      }
    }
  }
  else if (*shape_conversion == 'C'){
    float cutoff = 0.25;
    int i;

    sscanf(shape_conversion + 1, "%f", &cutoff);

    for(i = 0; indices[i]; ++i){
      index = indices[i];
      values[index] = values[index] < cutoff ? 0 : 1;
    }
  }
  else if (*shape_conversion == 'L' || *shape_conversion == 'O'){
    int i;
    float slope = (*shape_conversion == 'L') ? 0.68 : 1.6;
    float intercept = (*shape_conversion == 'L') ? 0.2 : -2.29;

    sc_parse_parameters(shape_conversion + 1, 's', 'i', &slope, &intercept);

    for(i = 0; indices[i]; ++i){
      double v;
      index = indices[i];

      v = (*shape_conversion == 'L') ? values[index] : log(values[index]);
      values[index] = MAX2(MIN2((v - intercept) / slope, 1),0);
    }
  }
  else
    ret = 0;

  free(indices);

  return ret;
}

PUBLIC void
vrna_sc_init(vrna_fold_compound *vc){

  unsigned int s;
  vrna_sc_t    *sc;

  if(vc){
    vrna_sc_remove(vc);

    switch(vc->type){
      case VRNA_VC_TYPE_SINGLE:     sc                    = (vrna_sc_t *)vrna_alloc(sizeof(vrna_sc_t));
                                    sc->free_energies     = NULL;
                                    sc->en_basepair       = NULL;
                                    sc->en_stack          = NULL;
                                    sc->exp_en_stack      = NULL;
                                    sc->boltzmann_factors = NULL;
                                    sc->exp_en_basepair   = NULL;
                                    sc->f                 = NULL;
                                    sc->exp_f             = NULL;
                                    sc->data              = NULL;

                                    vc->sc  = sc;
                                    break;

      case VRNA_VC_TYPE_ALIGNMENT:  vc->scs = (vrna_sc_t **)vrna_alloc(sizeof(vrna_sc_t*) * (vc->n_seq + 1));
                                    for(s = 0; s < vc->n_seq; s++){
                                      sc                    = (vrna_sc_t *)vrna_alloc(sizeof(vrna_sc_t));
                                      sc->free_energies     = NULL;
                                      sc->en_basepair       = NULL;
                                      sc->en_stack          = NULL;
                                      sc->exp_en_stack      = NULL;
                                      sc->boltzmann_factors = NULL;
                                      sc->exp_en_basepair   = NULL;
                                      sc->f                 = NULL;
                                      sc->exp_f             = NULL;
                                      sc->data              = NULL;

                                      vc->scs[s]  = sc;
                                    }
                                    break;
      default:                      /* do nothing */
                                    break;
    }
  }
}

PUBLIC void
vrna_sc_remove(vrna_fold_compound *vc){

  int s;

  if(vc){
    switch(vc->type){
      case  VRNA_VC_TYPE_SINGLE:    vrna_sc_free(vc->sc);
                                    vc->sc = NULL;
                                    break;
      case  VRNA_VC_TYPE_ALIGNMENT: if(vc->scs){
                                      for(s = 0; s < vc->n_seq; s++)
                                        vrna_sc_free(vc->scs[s]);
                                      free(vc->scs);
                                    }
                                    vc->scs = NULL;
                                    break;
      default:                      /* do nothing */
                                    break;
    }
  }
}

PUBLIC void
vrna_sc_free(vrna_sc_t *sc){

  int i;
  if(sc){
    if(sc->free_energies){
      for(i = 0; sc->free_energies[i]; free(sc->free_energies[i++]));
      free(sc->free_energies);
    }
    if(sc->boltzmann_factors){
      for(i = 0; sc->boltzmann_factors[i]; free(sc->boltzmann_factors[i++]));
      free(sc->boltzmann_factors);
    }
    if(sc->en_basepair)
      free(sc->en_basepair);

    if(sc->exp_en_basepair)
      free(sc->exp_en_basepair);

    if(sc->en_stack)
      free(sc->en_stack);

    if(sc->exp_en_stack)
      free(sc->exp_en_stack);

    free(sc);
  }
}

PUBLIC void
vrna_sc_add_bp(vrna_fold_compound *vc,
                        const double **constraints,
                        unsigned int options){
                        

  if(options & VRNA_CONSTRAINT_SOFT_MFE)
    sc_add_bp_mfe(vc, constraints, options);

  if(options & VRNA_CONSTRAINT_SOFT_PF)
    sc_add_bp_pf(vc, constraints, options);
}


PUBLIC  int
vrna_sc_SHAPE_add_zarringhalam( vrna_fold_compound *vc,
                                const double *reactivities,
                                double b,
                                double default_value,
                                const char *shape_conversion,
                                unsigned int options){

  int       i, j, n, ret;
  double    *pr, *up, **bp;
  vrna_md_t *md;

  ret = 0; /* error */

  if(vc && reactivities && (vc->type == VRNA_VC_TYPE_SINGLE)){
    n   = vc->length;
    md  = (options & VRNA_CONSTRAINT_SOFT_PF) ? &(vc->exp_params->model_details) : &(vc->params->model_details);

    /* first we copy over the reactivities to convert them into probabilities later on */
    pr = (double *)vrna_alloc(sizeof(double) * (n + 1));
    for(i=0; i<=n; i++)
      pr[i] = reactivities[i];

    if(vrna_sc_SHAPE_to_pr(shape_conversion, pr, n, default_value)){

      /*  now, convert them into pseudo free energies for unpaired, and
          paired nucleotides
      */
      up = (double *)vrna_alloc(sizeof(double) * (n + 1));
      bp = (double **)vrna_alloc(sizeof(double *) * (n + 1));
      for(i = 1; i <= n; ++i){
        up[i] = b * fabs(pr[i] - 1);
        bp[i] = (double *)vrna_alloc(sizeof(double) * (n + 1));
        for(j = i + md->min_loop_size + 1; j <= n; ++j)
          bp[i][j] = b * (pr[i] + pr[j]);
      }

      /* add the pseudo energies as soft constraints */
      vrna_sc_add_up(vc, (const double *)up, options);
      vrna_sc_add_bp(vc, (const double **)bp, options);

      /* clean up memory */
      for(i = 1; i <= n; ++i)
        free(bp[i]);
      free(bp);
      free(up);

      ret = 1; /* success */
    }

    free(pr);

  }

  return ret;
}


PUBLIC int
vrna_sc_SHAPE_add_deigan( vrna_fold_compound *vc,
                          const double *reactivities,
                          double m,
                          double b,
                          unsigned int options){

  int     i;
  double  *values;

  if(vc && values && (vc->type == VRNA_VC_TYPE_SINGLE)){

    values = (double *)vrna_alloc(sizeof(double) * (vc->length + 1));

    /* first convert the values according to provided slope and intercept values */
    for (i = 1; i <= vc->length; ++i){
      values[i] = reactivities[i] < 0 ? 0 : m * log(reactivities[i] + 1) + b;
    }

    if(options & VRNA_CONSTRAINT_SOFT_MFE)
      sc_add_stack_en_mfe(vc, (const double *)values, options);

    if(options & VRNA_CONSTRAINT_SOFT_PF)
      sc_add_stack_en_pf(vc, (const double *)values, options);

    free(values);
    return 1; /* success */
  } else {
    return 0; /* error */
  }
}

PUBLIC int
vrna_sc_SHAPE_add_deigan_ali( vrna_fold_compound *vc,
                              const char **shape_files,
                              const int *shape_file_association,
                              double m,
                              double b,
                              unsigned int options){

  float   reactivity, *reactivities, e1;
  char    *line, nucleotide, *sequence;
  int     s, i, p, r, position, *pseudo_energies, n_seq;

  if(vc->type == VRNA_VC_TYPE_ALIGNMENT){
    n_seq = vc->n_seq;

    vrna_sc_init(vc);

    for(s = 0; shape_file_association[s] != -1; s++){
      if(shape_file_association[s] > n_seq)
        vrna_message_warning("SHAPE file association exceeds sequence number in alignment");

      /* read the shape file */
      FILE *fp;
      if(!(fp = fopen(shape_files[s], "r"))){
        fprintf(stderr, "WARNING: SHAPE data file %d could not be opened. No shape data will be used.\n", s);
      } else {

        reactivities  = (float *)vrna_alloc(sizeof(float) * (vc->length + 1));
        sequence      = (char *)vrna_alloc(sizeof(char) * (vc->length + 1));

        for(i = 1; i <= vc->length; i++)
          reactivities[i] = -1.;

        while((line=get_line(fp))){
          r = sscanf(line, "%d %c %f", &position, &nucleotide, &reactivity);
          if(r){
            if((position <= 0) || (position > vc->length))
              vrna_message_error("provided shape data outside of sequence scope");

            switch(r){
              case 1:   nucleotide = 'N';
                        /* fall through */
              case 2:   reactivity = -1.;
                        /* fall through */
              default:  sequence[position-1]    = nucleotide;
                        reactivities[position]  = reactivity;
                        break;
            }
          }
          free(line);
        }
        fclose(fp);

        sequence[vc->length] = '\0';

        /* double check information by comparing the sequence read from */
        char *tmp_seq = get_ungapped_sequence(vc->sequences[shape_file_association[s]]);
        if(strcmp(tmp_seq, sequence)){
          fprintf(stderr, "WARNING: Input sequence %d differs from sequence provided via SHAPE file!\n", shape_file_association[s]);
        }
        free(tmp_seq);

        /* convert reactivities to pseudo energies */
        for(i = 1; i <= vc->length; i++){
          if(reactivities[i] < 0)
            reactivities[i] = 0.;
          else
            reactivities[i] = m * log(reactivities[i] + 1.) + b; /* this should be a value in kcal/mol */
        }

        /* begin actual storage of the pseudo energies */

        if(options & VRNA_CONSTRAINT_SOFT_MFE){
          pseudo_energies = (int *)vrna_alloc(sizeof(int) * (vc->length + 1));
          for(p = 0, i = 1; i<=vc->length; i++){
            e1 = (i - p > 0) ? reactivities[i - p] : 0.;
            if(vc->sequences[shape_file_association[s]][i-1] == '-'){
              p++; e1 = 0.;
            }
            pseudo_energies[i] = (int)(e1 * 100.);
          }
          vc->scs[shape_file_association[s]]->en_stack = pseudo_energies;
        }

        if(options & VRNA_CONSTRAINT_SOFT_PF){
          FLT_OR_DBL *exp_pe = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * (vc->length + 1));
          for(i=0;i<=vc->length;i++)
            exp_pe[i] = 1.;

          for(p = 0, i = 1; i<=vc->length; i++){
            e1 = (i - p > 0) ? reactivities[i - p] : 0.;
            if(vc->sequences[shape_file_association[s]][i-1] == '-'){
              p++; e1 = 0.;
            }
            exp_pe[i] = exp(-(e1 * 1000.) / vc->exp_params->kT );
          }
          vc->scs[shape_file_association[s]]->exp_en_stack = exp_pe;
        }
        
        free(reactivities);
      }
    }

    return 1; /* success */
  } else {
    return 0; /* error */
  }
}

PUBLIC  int
vrna_sc_SHAPE_parse_method( const char *method_string,
                            char *method,
                            float *param_1,
                            float *param_2){

  const char *params = method_string + 1;

  *param_1 = 0;
  *param_2 = 0;

  if (!method_string || !method_string[0])
    return 0;

  *method = method_string[0];

  switch(method_string[0]){
    case 'Z':   *param_1 = 0.89;
                sc_parse_parameters(params, 'b', '\0', param_1, NULL);
                break;

    case 'D':   *param_1 = 1.8;
                *param_2 = -0.6;
                sc_parse_parameters(params, 'm', 'b', param_1, param_2);
                break;

    case 'W':   break;

    default:    *method = 0;
                return 0;
  }

  return 1;
}

PUBLIC void
vrna_sc_add_up(vrna_fold_compound *vc,
                        const double *constraints,
                        unsigned int options){

  if(options & VRNA_CONSTRAINT_SOFT_MFE)
    sc_add_up_mfe(vc, constraints, options);

  if(options & VRNA_CONSTRAINT_SOFT_PF)
    sc_add_up_pf(vc, constraints, options);
}

PUBLIC void
vrna_sc_add_f(vrna_fold_compound *vc,
              int (*f)( int, int, int, int, char, void *),
              void *data){

  if(vc && f){
    if(vc->type == VRNA_VC_TYPE_SINGLE){
      if(!vc->sc)
        vrna_sc_init(vc);

      vc->sc->f       = f;
      if(data)
        vc->sc->data  = data;
    }
  }
}

PUBLIC void
vrna_sc_add_exp_f(vrna_fold_compound *vc,
                  FLT_OR_DBL (*exp_f)( int, int, int, int, char, void *),
                  void *data){

  if(vc && exp_f){
    if(vc->type == VRNA_VC_TYPE_SINGLE){
      if(!vc->sc)
        vrna_sc_init(vc);

      vc->sc->exp_f   = exp_f;
      if(data)
        vc->sc->data  = data;
    }
  }
}

PUBLIC void
vrna_sc_add_post( vrna_fold_compound *vc,
                  void (*post)( vrna_fold_compound *, char)){

  if(vc && post){
    if(vc->type == VRNA_VC_TYPE_SINGLE){
      if(!vc->sc)
        vrna_sc_init(vc);

      vc->sc->post       = post;
    }
  }
}

PUBLIC void
vrna_sc_add_pre(vrna_fold_compound *vc,
                void (*pre)( vrna_fold_compound *, char)){

  if(vc && pre){
    if(vc->type == VRNA_VC_TYPE_SINGLE){
      if(!vc->sc)
        vrna_sc_init(vc);

      vc->sc->pre = pre;
    }
  }
}

PRIVATE void
sc_parse_parameters( const char *string,
                        char c1,
                        char c2,
                        float *v1,
                        float *v2){

  char fmt[8];
  const char warning[] = "SHAPE method parameters not recognized! Using default parameters!";
  int r;

  assert(c1);
  assert(v1);

  if(!string || !(*string))
    return;

  if(c2 == 0 || v2 == NULL){
    sprintf(fmt, "%c%%f", c1);
    r = sscanf(string, fmt, v1);

    if(!r)
      vrna_message_warning(warning);

    return;
  }

  sprintf(fmt, "%c%%f%c%%f", c1, c2);
  r = sscanf(string, fmt, v1, v2);

  if(r!=2){
    sprintf(fmt, "%c%%f", c1);
    r = sscanf(string, fmt, v1);

    if(!r){
      sprintf(fmt, "%c%%f", c2);
      r = sscanf(string, fmt, v2);

      if(!r)
        vrna_message_warning(warning);
    }
  }
}

PRIVATE void
sc_add_bp_mfe(vrna_fold_compound *vc,
              const double **constraints,
              unsigned int options){

  unsigned int  i, j, n;
  vrna_sc_t     *sc;
  int           *idx;

  if(vc && constraints){
    n   = vc->length;

    if(!vc->sc)
      vrna_sc_init(vc);

    sc              = vc->sc;
    sc->en_basepair = (int *)vrna_alloc(sizeof(int) * (((n + 1) * (n + 2)) / 2));

    idx = vc->jindx;
    for(i = 1; i < n; i++)
      for(j=i+1; j<=n; j++)
        sc->en_basepair[idx[j]+i] = (int)(constraints[i][j] * 100.);

  }
}

PRIVATE void
sc_add_bp_pf( vrna_fold_compound *vc,
              const double **constraints,
              unsigned int options){

  unsigned int  i, j, n;
  vrna_sc_t     *sc;
  int           *idx;

  if(vc && constraints){
    n   = vc->length;

    if(!vc->sc)
      vrna_sc_init(vc);

    sc = vc->sc;

    vrna_exp_param_t  *exp_params = vc->exp_params;
    double            GT          = 0.;
    double            temperature = exp_params->temperature;
    double            kT          = exp_params->kT;
    double            TT          = (temperature+K0)/(Tmeasure);

    if(sc->exp_en_basepair)
      free(sc->exp_en_basepair);
    sc->exp_en_basepair     = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * (((n + 1) * (n + 2)) / 2));

    idx = vc->iindx;
    for(i = 1; i < n; i++)
      for(j=i+1; j<=n; j++){
        GT = constraints[i][j] * TT * 1000.;
        sc->exp_en_basepair[idx[i]-j] = exp( -GT / kT);
      }
  }
}

PRIVATE void
sc_add_up_mfe(vrna_fold_compound *vc,
              const double *constraints,
              unsigned int options){

  unsigned int  i, j, n;
  vrna_sc_t     *sc;

  if(vc && constraints){
    n   = vc->length;

    if(!vc->sc)
      vrna_sc_init(vc);

    sc  = vc->sc;

    /*  allocate memory such that we can access the soft constraint
        energies of a subsequence of length j starting at position i
        via sc->free_energies[i][j]
    */
    if(sc->free_energies){
      for(i = 0; i <= n; i++)
        if(sc->free_energies[i])
          free(sc->free_energies[i]);
      free(sc->free_energies);
    }

    sc->free_energies = (int **)vrna_alloc(sizeof(int *) * (n + 2));
    for(i = 0; i <= n; i++)
      sc->free_energies[i] = (int *)vrna_alloc(sizeof(int) * (n - i + 2));

    sc->free_energies[n+1] = NULL;

    for(i = 1; i <= n; i++){
      for(j = 1; j <= (n - i + 1); j++){
        sc->free_energies[i][j] =   sc->free_energies[i][j-1]
                                  + (int)(constraints[i+j-1] * 100); /* convert to 10kal/mol */
      }
    }
  }
}

PRIVATE void
sc_add_up_pf( vrna_fold_compound *vc,
              const double *constraints,
              unsigned int options){

  unsigned int  i, j, n;
  vrna_sc_t     *sc;

  if(vc && constraints){
    n   = vc->length;

    if(!vc->sc)
      vrna_sc_init(vc);

    sc  = vc->sc;

    vrna_exp_param_t   *exp_params = vc->exp_params;
    double             GT          = 0.;
    double             temperature = exp_params->temperature;
    double             kT          = exp_params->kT;
    double             TT          = (temperature+K0)/(Tmeasure);

    /* #################################### */
    /* # single nucleotide contributions  # */
    /* #################################### */

    /*  allocate memory such that we can access the soft constraint
        energies of a subsequence of length j starting at position i
        via sc->boltzmann_factors[i][j]
    */
    if(sc->boltzmann_factors){
      for(i = 0; i <= n; i++)
        if(sc->boltzmann_factors[i])
          free(sc->boltzmann_factors[i]);
      free(sc->boltzmann_factors);
    }

    sc->boltzmann_factors = (FLT_OR_DBL **)vrna_alloc(sizeof(FLT_OR_DBL *) * (n + 2));
    for(i = 0; i <= n; i++){
      sc->boltzmann_factors[i] = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * (n - i + 2));
      for(j = 0; j < n - i + 2; j++)
        sc->boltzmann_factors[i][j] = 1.;
    }

    sc->boltzmann_factors[n+1] = NULL;

    for(i = 1; i <= n; i++){
      for(j = 1; j <= (n - i + 1); j++){
        GT  = (double)((int)(constraints[i+j-1] * 100)) * TT * 10.; /* convert to cal/mol */
        sc->boltzmann_factors[i][j] =   sc->boltzmann_factors[i][j-1]
                                      * exp( -GT / kT);
      }
    }
  }
}

PRIVATE void
sc_add_stack_en_mfe(vrna_fold_compound *vc,
                    const double *constraints,
                    unsigned int options){
  int i;

  if(!vc->sc)
    vrna_sc_init(vc);

  if(!vc->sc->en_stack)
    vc->sc->en_stack = (int *)vrna_alloc(sizeof(int) * (vc->length + 1));

  for(i = 1; i <= vc->length; ++i)
    vc->sc->en_stack[i] += (int)(constraints[i] * 100.);
}

PRIVATE void
sc_add_stack_en_pf( vrna_fold_compound *vc,
                    const double *constraints,
                    unsigned int options){
  int i;

  if(!vc->sc)
    vrna_sc_init(vc);

  if(!vc->sc->exp_en_stack){
    vc->sc->exp_en_stack = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * (vc->length + 1));
    for(i = 0; i <= vc->length; ++i)
      vc->sc->exp_en_stack[i] = 1;
  }

  for(i = 1; i <= vc->length; ++i)
    vc->sc->exp_en_stack[i] *= exp(-(constraints[i] * 1000.)/ vc->exp_params->kT);
}

PRIVATE INLINE  void
adjust_ptypes(char *ptype,
              vrna_hc_t *hc,
              unsigned int length,
              unsigned int idx_type){

  unsigned  int i,j;
  int           *index;
  char          *matrix;

  matrix = hc->matrix;

  if(idx_type){
    index = vrna_get_iindx(length);
    for(i = 1; i < length; i++)
      for(j = i + 1; j <= length; j++)
        if(matrix[index[i] - j])
          if(!ptype[index[i] - j])
            ptype[index[i] - j] = 7; /* set to non-canonical pair */

  } else {
    index = vrna_get_indx(length);
    for(i = 1; i < length; i++)
      for(j = i + 1; j <= length; j++)
        if(matrix[index[j] + i])
          if(!ptype[index[j] + i])
            ptype[index[j] + i] = 7; /* set to non-canonical pair */

  }
  free(index);
}
  
#ifdef  VRNA_BACKWARD_COMPAT

/*###########################################*/
/*# deprecated functions below              #*/
/*###########################################*/

PUBLIC  void
print_tty_constraint_full(void){

  vrna_message_constraint_options_all();
}

PUBLIC  void
print_tty_constraint(unsigned int option){

  vrna_message_constraint_options(option);
}

PUBLIC void
constrain_ptypes( const char *constraint,
                  unsigned int length,
                  char *ptype,
                  int *BP,
                  int min_loop_size,
                  unsigned int idx_type){

  apply_DB_constraint(constraint,
                      ptype,
                      length,
                      (unsigned int)min_loop_size,
                      -1,
                      VRNA_CONSTRAINT_PIPE
                    | VRNA_CONSTRAINT_DOT
                    | VRNA_CONSTRAINT_X
                    | VRNA_CONSTRAINT_ANG_BRACK
                    | VRNA_CONSTRAINT_RND_BRACK);
}

#endif
