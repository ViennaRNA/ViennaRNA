/* constraints handling */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <assert.h>
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
#include "ViennaRNA/constraints_hard.h"


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
PRIVATE void
hc_add_up(vrna_fold_compound_t *vc,
          int i,
          char option);

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
                  unsigned int length,
                  int *index);

PRIVATE INLINE  void
hc_pairs_downstream(unsigned int i,
                    char c_option,
                    char *hc,
                    unsigned int length,
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

PRIVATE void
apply_DB_constraint(const char *constraint,
                    char *ptype,
                    unsigned int length,
                    unsigned int min_loop_size,
                    int cut,
                    unsigned int options);

PRIVATE void
hc_reset_to_default(vrna_fold_compound_t *vc);

PRIVATE void
hc_update_up(vrna_fold_compound_t *vc);

/*
#################################
# BEGIN OF FUNCTION DEFINITIONS #
#################################
*/
PUBLIC  void
vrna_message_constraint_options_all(void){

  vrna_message_constraint_options(  VRNA_CONSTRAINT_DB_PIPE
                                  | VRNA_CONSTRAINT_DB_DOT
                                  | VRNA_CONSTRAINT_DB_X
                                  | VRNA_CONSTRAINT_DB_ANG_BRACK
                                  | VRNA_CONSTRAINT_DB_RND_BRACK);
}

PUBLIC  void
vrna_message_constraint_options(unsigned int option){

  printf("Input structure constraints using the following notation:\n");
  if(option & VRNA_CONSTRAINT_DB_PIPE)       printf("| : paired with another base\n");
  if(option & VRNA_CONSTRAINT_DB_DOT)        printf(". : no constraint at all\n");
  if(option & VRNA_CONSTRAINT_DB_X)          printf("x : base must not pair\n");
  if(option & VRNA_CONSTRAINT_DB_ANG_BRACK)  printf("< : base i is paired with a base j<i\n> : base i is paired with a base j>i\n");
  if(option & VRNA_CONSTRAINT_DB_RND_BRACK)  printf("matching brackets ( ): base i pairs base j\n");
}

PUBLIC  void
vrna_hc_init(vrna_fold_compound_t *vc){

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

  /* add null pointers for the generalized hard constraint feature */
  hc->f           = NULL;
  hc->data        = NULL;
  hc->free_data   = NULL;

  /* update */
  hc_update_up(vc);
}

PUBLIC void
vrna_hc_add_up( vrna_fold_compound_t *vc,
                int i,
                char option){

  int j;

  if(vc)
    if(vc->hc){
      if((i <= 0) || (i > vc->length)){
        vrna_message_warning("vrna_hc_add_up: position out of range, not doing anything");
        return;
      }

      hc_add_up(vc, i, option);

      hc_update_up(vc);
    }
}

PUBLIC int
vrna_hc_add_up_batch( vrna_fold_compound_t *vc,
                      vrna_hc_up_t *constraints){

  int i, ret;

  ret = 0; /* failure */

  if(vc)
    if(vc->hc && constraints){
      for(i = 0; constraints[i].position != 0; i++){
        int pos       = constraints[i].position;
        char options  = constraints[i].options;
        if((pos <= 0) || (pos > vc->length)){
          vrna_message_warning("vrna_hc_add_up_batch: position out of range, application of hard constraints stops here!");
          return ret;
        }
        hc_add_up(vc, pos, options);
      }

      hc_update_up(vc);
      ret = 1; /* success */
    }

  return ret;
}

PRIVATE void
hc_add_up(vrna_fold_compound_t *vc,
          int i,
          char option){

  int   j;
  char  type = (char)0;

  if(option & VRNA_CONSTRAINT_CONTEXT_ENFORCE){ /* force nucleotide to appear unpaired within a certain type of loop */
    /* do not allow i to be paired with any other nucleotide */
    if(!(option & VRNA_CONSTRAINT_CONTEXT_NO_REMOVE)){
      for(j = 1; j < i; j++)
        vc->hc->matrix[vc->jindx[i] + j] = (char)0;
      for(j = i+1; j <= vc->length; j++)
        vc->hc->matrix[vc->jindx[j] + i] = (char)0;
    }

    type = option & (char)( VRNA_CONSTRAINT_CONTEXT_EXT_LOOP
                            | VRNA_CONSTRAINT_CONTEXT_HP_LOOP
                            | VRNA_CONSTRAINT_CONTEXT_INT_LOOP
                            | VRNA_CONSTRAINT_CONTEXT_MB_LOOP);

    vc->hc->matrix[vc->jindx[i] + i] = type;
  } else {
    type = option & VRNA_CONSTRAINT_CONTEXT_ALL_LOOPS;

    /* do not allow i to be paired with any other nucleotide (in context type) */
    if(!(option & VRNA_CONSTRAINT_CONTEXT_NO_REMOVE)){
      for(j = 1; j < i; j++)
        vc->hc->matrix[vc->jindx[i] + j] &= ~type;
      for(j = i+1; j <= vc->length; j++)
        vc->hc->matrix[vc->jindx[j] + i] &= ~type;
    }

    vc->hc->matrix[vc->jindx[i] + i] = (char)(  VRNA_CONSTRAINT_CONTEXT_EXT_LOOP
                                              | VRNA_CONSTRAINT_CONTEXT_HP_LOOP
                                              | VRNA_CONSTRAINT_CONTEXT_INT_LOOP
                                              | VRNA_CONSTRAINT_CONTEXT_MB_LOOP);
  }
}

PUBLIC void
vrna_hc_add_bp_nonspecific( vrna_fold_compound_t *vc,
                            int i,
                            int d,
                            char option){
  int   p;
  char  type, t1, t2;

  if(vc)
    if(vc->hc){
      if((i <= 0) || (i > vc->length)){
        vrna_message_warning("vrna_hc_add_bp_nonspecific: position out of range, not doing anything");
        return;
      }

      /* position i may pair in provided contexts */
      type  = option & VRNA_CONSTRAINT_CONTEXT_ALL_LOOPS;
      /* acknowledge pairing direction */
      t1    = (d <= 0) ? type : (char)0;
      t2    = (d >= 0) ? type : (char)0;

      if(option & VRNA_CONSTRAINT_CONTEXT_NO_REMOVE){
        /* only allow for possibly non-canonical pairs, do not enforce them */
        for(p = 1; p < i; p++)
          vc->hc->matrix[vc->jindx[i] + p] |= t1;
        for(p = i+1; p <= vc->length; p++)
          vc->hc->matrix[vc->jindx[p] + i] |= t2;
      } else {
        /* force pairing direction */
        for(p = 1; p < i; p++)
          vc->hc->matrix[vc->jindx[i] + p] &= t1;
        for(p = i+1; p <= vc->length; p++)
          vc->hc->matrix[vc->jindx[p] + i] &= t2;
        /* nucleotide mustn't be unpaired */
        vc->hc->matrix[vc->jindx[i] + i] = (char)0;
      }

      hc_update_up(vc);
    }

}

PUBLIC void
vrna_hc_add_bp( vrna_fold_compound_t *vc,
                int i,
                int j,
                char option){

  int   k, l;
  char  type;

  if(vc)
    if(vc->hc){
      if((i <= 0) || (j <= i) || (j > vc->length)){
        vrna_message_warning("vrna_hc_add_bp: position out of range, not doing anything");
        return;
      }

      /* reset ptype in case (i,j) is a non-canonical pair */
      if(option & VRNA_CONSTRAINT_CONTEXT_ALL_LOOPS){
        if(vc->hc->matrix[vc->jindx[j] + i])
          if(vc->ptype[vc->jindx[j] + i] == 0)
            vc->ptype[vc->jindx[j] + i] = 7;
      }

      vc->hc->matrix[vc->jindx[j] + i] = option & VRNA_CONSTRAINT_CONTEXT_ALL_LOOPS;

      if(!(option & VRNA_CONSTRAINT_CONTEXT_NO_REMOVE)){
        /*
          remove all conflicting base pairs, i.e. do not allow i,j to pair
          with any other nucleotide k
        */
        for(k = 1; k < i; k++){
          vc->hc->matrix[vc->jindx[i] + k] = (char)0;
          vc->hc->matrix[vc->jindx[j] + k] = (char)0;
          for(l = i+1; l < j; l++)
            vc->hc->matrix[vc->jindx[l] + k] = (char)0;
        }
        for(k = i+1; k < j; k++){
          vc->hc->matrix[vc->jindx[k] + i] = (char)0;
          vc->hc->matrix[vc->jindx[j] + k] = (char)0;
          for(l = j + 1; l <= vc->length; l++)
            vc->hc->matrix[vc->jindx[l] + k] = (char)0;
        }
        for(k = j+1; k <= vc->length; k++){
          vc->hc->matrix[vc->jindx[k] + i] = (char)0;
          vc->hc->matrix[vc->jindx[k] + j] = (char)0;
        }
      }

      if(option & VRNA_CONSTRAINT_CONTEXT_ENFORCE){

        /* do not allow i,j to be unpaired */
        vc->hc->matrix[vc->jindx[i] + i] = (char)0;
        vc->hc->matrix[vc->jindx[j] + j] = (char)0;

        hc_update_up(vc);
      }
    }
}

PUBLIC void
vrna_hc_free(vrna_hc_t *hc){

  if(hc){
    free(hc->matrix);
    free(hc->up_ext);
    free(hc->up_hp);
    free(hc->up_int);
    free(hc->up_ml);

    if(hc->free_data)
      hc->free_data(hc->data);

    free(hc);
  }
}


PUBLIC void
vrna_hc_add_f(vrna_fold_compound_t *vc,
              vrna_callback_hc_evaluate *f)
{
  if (vc && f) {
    if (vc->type == VRNA_FC_TYPE_SINGLE) {
      if (!vc->hc)
        vrna_hc_init(vc);

      vc->hc->f = f;
    }
  }
}


PUBLIC void
vrna_hc_add_data( vrna_fold_compound_t *vc,
                  void *data,
                  vrna_callback_free_auxdata *f)
{
  if (vc && data) {
    if (vc->type == VRNA_FC_TYPE_SINGLE) {
      if (!vc->hc)
        vrna_hc_init(vc);

      vc->hc->data        = data;
      vc->hc->free_data   = f;
    }
  }
}


PUBLIC  int
vrna_hc_add_from_db(vrna_fold_compound_t *vc,
                    const char *constraint,
                    unsigned int options){

  int         i, d, ret;
  vrna_md_t   *md;

  ret = 0; /* Failure */

  if(vc){
    if(vc->params)
      md = &(vc->params->model_details);
    else if(vc->exp_params)
      md = &(vc->exp_params->model_details);
    else
      return ret;

    if(!vc->hc)
      vrna_hc_init(vc);

    /* apply hard constraints from dot-bracket notation */
    apply_DB_constraint(constraint,
                        vc->hc->matrix,
                        vc->length,
                        md->min_loop_size,
                        -1,
                        options);
    hc_update_up(vc);
    ret = 1; /* Success */
  }

  return ret;
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
  index     = vrna_idx_col_wise(length);
  c_option  =   VRNA_CONSTRAINT_CONTEXT_EXT_LOOP
              | VRNA_CONSTRAINT_CONTEXT_HP_LOOP
              | VRNA_CONSTRAINT_CONTEXT_INT_LOOP
              | VRNA_CONSTRAINT_CONTEXT_INT_LOOP_ENC
              | VRNA_CONSTRAINT_CONTEXT_MB_LOOP
              | VRNA_CONSTRAINT_CONTEXT_MB_LOOP_ENC;

  for(hx=0, j=1; j<=n; j++) {
    switch (constraint[j-1]) {
       /* can't pair */
       case 'x':  if(options & VRNA_CONSTRAINT_DB_X){
                    hc_cant_pair(j, c_option, hc, length, min_loop_size, index);
                  }
                  break;

      /* must pair, i.e. may not be unpaired */
      case '|':   if(options & VRNA_CONSTRAINT_DB_PIPE){
                    if(options & VRNA_CONSTRAINT_DB_ENFORCE_BP)
                      hc_must_pair(j, c_option, hc, index);
                  }
                  break;

      /* weak enforced pair 'open' */
      case '(':   if(options & VRNA_CONSTRAINT_DB_RND_BRACK){
                    stack[hx++]=j;
                  }
                  break;

      /* weak enforced pair 'close' */
      case ')':   if(options & VRNA_CONSTRAINT_DB_RND_BRACK){
                    if (hx<=0) {
                      vrna_message_error("%s\nunbalanced brackets in constraints", constraint);
                    }
                    i = stack[--hx];
                    if(options & VRNA_CONSTRAINT_DB_ENFORCE_BP)
                      hc_enforce_pair(i, j, c_option, hc, length, min_loop_size, index);
                    else
                      hc_weak_enforce_pair(i, j, c_option, hc, length, min_loop_size, index);
                  }
                  break;

      /* pairs upstream */
      case '<':   if(options & VRNA_CONSTRAINT_DB_ANG_BRACK){
                    hc_pairs_downstream(j, c_option, hc, length, index);
                    if(options & VRNA_CONSTRAINT_DB_ENFORCE_BP)
                      hc_must_pair(j, c_option, hc, index);
                  }
                  break;

      /* pairs downstream */
      case '>':   if(options & VRNA_CONSTRAINT_DB_ANG_BRACK){
                    hc_pairs_upstream(j, c_option, hc, length, index);
                    if(options & VRNA_CONSTRAINT_DB_ENFORCE_BP)
                      hc_must_pair(j, c_option, hc, index);
                  }
                  break;

      /* only intramolecular basepairing */
      case 'l':   if(options & VRNA_CONSTRAINT_DB_INTRAMOL){
                    hc_intramolecular_only(j, c_option, hc, length, min_loop_size, cut, index);
                  }
                  break;

      /* only intermolecular bp */
      case 'e':   if(options & VRNA_CONSTRAINT_DB_INTERMOL){
                    hc_intermolecular_only(j, c_option, hc, length, min_loop_size, cut, index);
                  }
                  break;

      case '.':   break;

      default:    vrna_message_warning("Unrecognized character '%c' in pseudo dot-bracket notation constraint string",
                                              constraint[j-1]);
                  break;
    }
  }

  if (hx!=0) {
    vrna_message_error("%s\nunbalanced brackets in constraint string", constraint);
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

  hc_pairs_upstream(i, c_option, hc, length, index);
  hc_pairs_downstream(i, c_option, hc, length, index);
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
                  unsigned int length,
                  int *index){

  unsigned int l;

  /* prohibit downstream pairs */
  for(l = length; l > i; l--)
    hc[index[l] + i] = (char)0;
  /* allow upstream pairs of given type */
  for(l = i - 1; l >= 1; l--)
    hc[index[i] + l] &= c_option;
}

PRIVATE INLINE  void
hc_pairs_downstream(unsigned int i,
                    char c_option,
                    char *hc,
                    unsigned int length,
                    int *index){

  unsigned int l;
  /* allow downstream pairs of given type */
  for(l = length; l > i; l--)
    hc[index[l] + i] &= c_option;
  /* forbid upstream pairs */
  for(l = i - 1; l >= 1; l--)
    hc[index[i] + l] = (char)0;
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
  /* don't allow pairs (i,k) i < k <= n */ 
  hc_pairs_upstream(i, (char)0, hc, length, index);
  /* don't allow pairs (k,j) 1 <= k < j */
  /* don't allow pairs (j,k) j < k <= n */ 
  hc_pairs_upstream(j, (char)0, hc, length, index);

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
hc_reset_to_default(vrna_fold_compound_t *vc){

  unsigned int      i, j, ij, min_loop_size, n;
  int               max_span, *idx;
  vrna_md_t         *md;
  vrna_hc_t         *hc;
  short             *S;

  md  = NULL;
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
    hc->matrix[idx[i] + i]  =   VRNA_CONSTRAINT_CONTEXT_EXT_LOOP
                              | VRNA_CONSTRAINT_CONTEXT_HP_LOOP
                              | VRNA_CONSTRAINT_CONTEXT_INT_LOOP
                              | VRNA_CONSTRAINT_CONTEXT_MB_LOOP;

  /* 2. all base pairs with pscore above threshold are allowed in all contexts */
  switch(vc->type){
    case VRNA_FC_TYPE_COMPARATIVE:  for(j = n; j > min_loop_size + 1; j--){
                                    ij = idx[j]+1;
                                    for(i=1; i < j - min_loop_size; i++, ij++){
                                      char opt = (char)0;
                                      if((j-i+1) <= max_span){
                                        if(vc->pscore[idx[j]+i] >= md->cv_fact*MINPSCORE)
                                          opt = VRNA_CONSTRAINT_CONTEXT_ALL_LOOPS;
                                      }
                                      hc->matrix[ij] = opt;
                                    }
                                  }
                                  break;

    case VRNA_FC_TYPE_SINGLE:     for(j = n; j > min_loop_size + 1; j--){
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
                                                      opt = VRNA_CONSTRAINT_CONTEXT_ALL_LOOPS;
                                                      opt &= ~(VRNA_CONSTRAINT_CONTEXT_HP_LOOP | VRNA_CONSTRAINT_CONTEXT_MB_LOOP);
                                                      break;
                                                    } /* else fallthrough */
                                          default:  opt = VRNA_CONSTRAINT_CONTEXT_ALL_LOOPS;
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
                                  /*  ???????
                                      Is this necessary? We could leave the noLP option somewhere else, i.e. do not enforce it
                                      on the level of ptype/constraints, but an the level of recursions...
                                      ???????
                                  */
                                  if(md->noLP){
                                    if(!vc->ptype)
                                      vc->ptype = vrna_ptypes(vc->sequence_encoding2, md);
                                    for(i = 1; i < n; i++)
                                      for(j = i + min_loop_size + 1; j <= n; j++){
                                        if(hc->matrix[idx[j] +i]){
                                          if(!vc->ptype[idx[j] + i]){
                                            hc->matrix[idx[j] + i] = (char)0;
                                          }
                                        }
                                      }
                                  }
                                  break;

    default:                      break;
  }

  /* should we reset the generalized hard constraint feature here? */
  if(hc->f || hc->data){
    if(hc->free_data)
      hc->free_data(hc->data);

    hc->f           = NULL;
    hc->data        = NULL;
    hc->free_data   = NULL;
  }

}

PRIVATE void
hc_update_up(vrna_fold_compound_t *vc){

  unsigned int      i, n;
  int               *idx;
  vrna_hc_t         *hc;

  n   = vc->length;
  idx = vc->jindx;
  hc  = vc->hc;

  for(hc->up_ext[n+1] = 0, i = n; i > 0; i--) /* unpaired stretch in exterior loop */
    hc->up_ext[i] = (hc->matrix[idx[i]+i] & VRNA_CONSTRAINT_CONTEXT_EXT_LOOP) ? 1 + hc->up_ext[i+1] : 0;

  for(hc->up_hp[n+1] = 0, i = n; i > 0; i--)  /* unpaired stretch in hairpin loop */
    hc->up_hp[i] = (hc->matrix[idx[i]+i] & VRNA_CONSTRAINT_CONTEXT_HP_LOOP) ? 1 + hc->up_hp[i+1] : 0;

  for(hc->up_int[n+1] = 0, i = n; i > 0; i--) /* unpaired stretch in interior loop */
    hc->up_int[i] = (hc->matrix[idx[i]+i] & VRNA_CONSTRAINT_CONTEXT_INT_LOOP) ? 1 + hc->up_int[i+1] : 0;

  for(hc->up_ml[n+1] = 0, i = n; i > 0; i--)  /* unpaired stretch in multibranch loop */
    hc->up_ml[i] = (hc->matrix[idx[i]+i] & VRNA_CONSTRAINT_CONTEXT_MB_LOOP) ? 1 + hc->up_ml[i+1] : 0;

  /*
   *  loop arround once more until we find a nucleotide that mustn't
   *  be unpaired (needed for circular folding)
   */

  if(hc->matrix[idx[1]+1] & VRNA_CONSTRAINT_CONTEXT_EXT_LOOP){
    hc->up_ext[n+1] = hc->up_ext[1];
    for(i = n; i > 0; i--){
      if(hc->matrix[idx[i]+i] & VRNA_CONSTRAINT_CONTEXT_EXT_LOOP){
        hc->up_ext[i] = MIN2(n, 1 + hc->up_ext[i+1]);
      } else
        break;
    }
  }

  if(hc->matrix[idx[1]+1] & VRNA_CONSTRAINT_CONTEXT_HP_LOOP){
    hc->up_hp[n+1] = hc->up_hp[1];
    for(i = n; i > 0; i--){
      if(hc->matrix[idx[i]+i] & VRNA_CONSTRAINT_CONTEXT_HP_LOOP){
        hc->up_hp[i] = MIN2(n, 1 + hc->up_hp[i+1]);
      } else
        break;
    }
  }

  if(hc->matrix[idx[1]+1] & VRNA_CONSTRAINT_CONTEXT_INT_LOOP){
    hc->up_int[n+1] = hc->up_int[1];
    for(i = n; i > 0; i--){
      if(hc->matrix[idx[i]+i] & VRNA_CONSTRAINT_CONTEXT_INT_LOOP){
        hc->up_int[i] = MIN2(n, 1 + hc->up_int[i+1]);
      } else
        break;
    }
  }

  if(hc->matrix[idx[1]+1] & VRNA_CONSTRAINT_CONTEXT_MB_LOOP){
    hc->up_ml[n+1] = hc->up_ml[1];
    for(i = n; i > 0; i--){
      if(hc->matrix[idx[i]+i] & VRNA_CONSTRAINT_CONTEXT_MB_LOOP){
        hc->up_ml[i] = MIN2(n, 1 + hc->up_ml[i+1]);
      } else
        break;
    }
  }

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

  int n,i,j,k,l;
  int hx, *stack;
  char type;
  int *index;

  if(constraint == NULL) return;

  n = (int)strlen(constraint);

  stack = vrna_alloc(sizeof(int)*(n+1));

  if(!idx_type){ /* index allows access in energy matrices at pos (i,j) via index[j]+i */
    index = vrna_idx_col_wise(length);

    for(hx=0, j=1; j<=n; j++){
      switch(constraint[j-1]){
        case '|':   if(BP) BP[j] = -1;
                    break;
        case 'x':   /* can't pair */
                    for (l=1; l<j-min_loop_size; l++)
                      ptype[index[j]+l] = 0;
                    for (l=j+min_loop_size+1; l<=(int)length; l++)
                      ptype[index[l]+j] = 0;
                    break;
        case '(':   stack[hx++]=j;
                    /* fallthrough */
        case '<':   /* pairs upstream */
                    for (l=1; l<j-min_loop_size; l++)
                      ptype[index[j]+l] = 0;
                    break;
        case ')':   if (hx<=0) {
                      vrna_message_error("%s\nunbalanced brackets in constraint", constraint);
                    }
                    i = stack[--hx];
                    type = ptype[index[j]+i];
                    for (k=i+1; k<=(int)length; k++)
                      ptype[index[k]+i] = 0;
                    /* don't allow pairs i<k<j<l */
                    for (l=j; l<=(int)length; l++)
                      for (k=i+1; k<=j; k++)
                        ptype[index[l]+k] = 0;
                    /* don't allow pairs k<i<l<j */
                    for (l=i; l<=j; l++)
                      for (k=1; k<=i; k++)
                        ptype[index[l]+k] = 0;
                    for (k=1; k<j; k++)
                      ptype[index[j]+k] = 0;
                    ptype[index[j]+i] = (type==0) ? 7 : type;
                    /* fallthrough */
        case '>':   /* pairs downstream */
                    for (l=j+min_loop_size+1; l<=(int)length; l++)
                      ptype[index[l]+j] = 0;
                    break;
      }
    }
  }
  else{ /* index allows access in energy matrices at pos (i,j) via index[i]-j */
    index = vrna_idx_row_wise(length);

    for(hx=0, j=1; j<=n; j++) {
      switch (constraint[j-1]) {
        case 'x':   /* can't pair */
                    for (l=1; l<j-min_loop_size; l++)
                      ptype[index[l]-j] = 0;
                    for (l=j+min_loop_size+1; l<=(int)length; l++)
                      ptype[index[j]-l] = 0;
                    break;
        case '(':   stack[hx++]=j;
                    /* fallthrough */
        case '<':   /* pairs upstream */
                    for (l=1; l<j-min_loop_size; l++)
                      ptype[index[l]-j] = 0;
                    break;
        case ')':   if (hx<=0) {
                      vrna_message_error("%s\nunbalanced brackets in constraints", constraint);
                    }
                    i = stack[--hx];
                    type = ptype[index[i]-j];
                    /* don't allow pairs i<k<j<l */
                    for (k=i; k<=j; k++)
                      for (l=j; l<=(int)length; l++)
                        ptype[index[k]-l] = 0;
                    /* don't allow pairs k<i<l<j */
                    for (k=1; k<=i; k++)
                      for (l=i; l<=j; l++)
                        ptype[index[k]-l] = 0;
                    ptype[index[i]-j] = (type==0) ? 7 : type;
                    /* fallthrough */
        case '>':   /* pairs downstream */
                    for (l=j+min_loop_size+1; l<=(int)length; l++)
                      ptype[index[j]-l] = 0;
                    break;
      }
    }
  }
  if (hx!=0) {
    vrna_message_error("%s\nunbalanced brackets in constraint string", constraint);
  }
  free(index);
  free(stack);
}

#endif
