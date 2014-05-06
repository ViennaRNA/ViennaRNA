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
              hard_constraintT *hc,
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

/*
#################################
# BEGIN OF FUNCTION DEFINITIONS #
#################################
*/
PUBLIC  void
print_tty_constraint_full(void){

  print_tty_constraint(   VRNA_CONSTRAINT_PIPE
                        | VRNA_CONSTRAINT_DOT
                        | VRNA_CONSTRAINT_X
                        | VRNA_CONSTRAINT_ANG_BRACK
                        | VRNA_CONSTRAINT_RND_BRACK);
}

PUBLIC  void
print_tty_constraint(unsigned int option){

  if(!(option & VRNA_CONSTRAINT_NO_HEADER)) printf("Input structure constraints using the following notation:\n");
  if(option & VRNA_CONSTRAINT_PIPE)       printf("| : paired with another base\n");
  if(option & VRNA_CONSTRAINT_DOT)        printf(". : no constraint at all\n");
  if(option & VRNA_CONSTRAINT_X)          printf("x : base must not pair\n");
  if(option & VRNA_CONSTRAINT_ANG_BRACK)  printf("< : base i is paired with a base j<i\n> : base i is paired with a base j>i\n");
  if(option & VRNA_CONSTRAINT_RND_BRACK)  printf("matching brackets ( ): base i pairs base j\n");
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
  stack     = (int *) space(sizeof(int)*(n+1));
  index     = get_indx(length);
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
                      nrerror("unbalanced brackets in constraints");
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

      default:    warn_user("unrecognized character in pseudo dot-bracket notation constraint string\n");
                  break;
    }
  }

  if (hx!=0) {
    fprintf(stderr, "%s\n", constraint);
    nrerror("unbalanced brackets in constraint string");
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

PUBLIC void
getConstraint(char **cstruc,
              const char **lines,
              unsigned int option){

  int r, i, l, cl, stop;
  char *c, *ptr;
  if(lines){
    if(option & VRNA_CONSTRAINT_ALL)
      option |=   VRNA_CONSTRAINT_PIPE
                | VRNA_CONSTRAINT_ANG_BRACK
                | VRNA_CONSTRAINT_RND_BRACK
                | VRNA_CONSTRAINT_X
                | VRNA_CONSTRAINT_INTRAMOLECULAR
                | VRNA_CONSTRAINT_INTERMOLECULAR;

    for(r=i=stop=0;lines[i];i++){
      l   = (int)strlen(lines[i]);
      c   = (char *) space(sizeof(char) * (l+1));
      (void) sscanf(lines[i], "%s", c);
      cl  = (int)strlen(c);
      /* line commented out ? */
      if((*c == '#') || (*c == '%') || (*c == ';') || (*c == '/') || (*c == '*' || (*c == '\0'))){
        /* skip leading comments only, i.e. do not allow comments inside the constraint */
        if(!r)  continue;
        else    break;
      }

      /* check current line for actual constraining structure */
      for(ptr = c;*c;c++){
        switch(*c){
          case '|':   if(!(option & VRNA_CONSTRAINT_PIPE)){
                        warn_user("constraints of type '|' not allowed");
                        *c = '.';
                      }
                      break;
          case '<':   
          case '>':   if(!(option & VRNA_CONSTRAINT_ANG_BRACK)){
                        warn_user("constraints of type '<' or '>' not allowed");
                        *c = '.';
                      }
                      break;
          case '(':
          case ')':   if(!(option & VRNA_CONSTRAINT_RND_BRACK)){
                        warn_user("constraints of type '(' or ')' not allowed");
                        *c = '.';
                      }
                      break;
          case 'x':   if(!(option & VRNA_CONSTRAINT_X)){
                        warn_user("constraints of type 'x' not allowed");
                        *c = '.';
                      }
                      break;
          case 'e':   if(!(option & VRNA_CONSTRAINT_INTERMOLECULAR)){
                        warn_user("constraints of type 'e' not allowed");
                        *c = '.';
                      }
                      break;
          case 'l':   if(!(option & VRNA_CONSTRAINT_INTRAMOLECULAR)){
                        warn_user("constraints of type 'l' not allowed");
                        *c = '.';
                      }
                      break;  /*only intramolecular basepairing */
          case '.':   break;
          case '&':   break; /* ignore concatenation char */
          default:    warn_user("unrecognized character in constraint structure");
                      break;
        }
      }

      r += cl+1;
      *cstruc = (char *)xrealloc(*cstruc, r*sizeof(char));
      strcat(*cstruc, ptr);
      free(ptr);
      /* stop if not in fasta mode or multiple words on line */
      if(!(option & VRNA_CONSTRAINT_MULTILINE) || (cl != l)) break;
    }
  }
}

PRIVATE void
hc_reset_to_default(vrna_fold_compound *vc){

  unsigned int      i, j, ij, min_loop_size, n;
  int               max_span, *idx;
  model_detailsT    *md;
  hard_constraintT  *hc;
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
    nrerror("missing model_details in fold_compound");

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
/*           hc->matrix[ij] = (vc->pscore[idx[j]+i] >= md->cv_fact*MINPSCORE) ? VRNA_HC_CONTEXT_ALL_LOOPS : (char)0; */
          hc->matrix[ij] = (vc->pscore[idx[j]+i] >= MINPSCORE) ? VRNA_HC_CONTEXT_ALL_LOOPS : (char)0;
        }
    }
  } else {
    /* 2. all canonical base pairs are allowed in all contexts */
    for(j = n; j > min_loop_size + 1; j--){
      ij = idx[j]+1;
      for(i=1; i < j - min_loop_size; i++, ij++)
        if((j-i+1) > max_span){
          hc->matrix[ij] = (char)0;
        } else {
          hc->matrix[ij] = md->pair[S[i]][S[j]] ? VRNA_HC_CONTEXT_ALL_LOOPS : (char)0;
        }
    }
  }
}

PRIVATE void
hc_update_up(vrna_fold_compound *vc){

  unsigned int      i, n;
  int               *idx;
  hard_constraintT  *hc;

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

PUBLIC  void
vrna_hc_add(vrna_fold_compound *vc,
            const char *constraint,
            unsigned int options){

  unsigned int      n, min_loop_size;
  hard_constraintT  *hc;
  model_detailsT    *md;

  if(vc->params)
    md = &(vc->params->model_details);
  else if(vc->exp_params)
    md = &(vc->exp_params->model_details);
  else
    nrerror("constraints.c@vrna_hc_add: fold compound has no params or exp_params");

  min_loop_size = md->min_loop_size;
  n             = vc->length;

  /* allocate memory for the hard constraints data structure */

  hc          = (hard_constraintT *)space(sizeof(hard_constraintT));
  hc->matrix  = (char *)space(sizeof(char)*((n*(n+1))/2+2));
  hc->up_ext  = (int *)space(sizeof(int)*(n+2));
  hc->up_hp   = (int *)space(sizeof(int)*(n+2));
  hc->up_int  = (int *)space(sizeof(int)*(n+2));
  hc->up_ml   = (int *)space(sizeof(int)*(n+2));

  /* set new hard constraints */
  vrna_hc_free(vc->hc);
  vc->hc = hc;

  /* ####################### */
  /* prefill default values  */
  /* ####################### */

  hc_reset_to_default(vc);

  /* ############################### */
  /* apply user supplied constraints */
  /* ############################### */

  if(options & VRNA_CONSTRAINT_DB){ /* constraints from dot-bracket notation */
/*
    printf("reading constraints from dot-bracket\n");
*/
    apply_DB_constraint(constraint,
                        hc->matrix,
                        n,
                        min_loop_size,
                        -1,
                        options);
  }

  if(options & VRNA_CONSTRAINT_FILE){ /* constraints from file */
/*
    printf("reading constraints from file\n");
*/
  }

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
        warn_user("vrna_hc_add_up: position out of range, not doing anything");
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
        warn_user("vrna_hc_add_bp: position out of range, not doing anything");
        return;
      }

      if(i == j){
        warn_user("vrna_hc_add_bp: positions do not differ, not doing anything");
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
vrna_hc_reset(vrna_fold_compound *vc){

  if(vc)
    if(vc->hc){

      hc_reset_to_default(vc);
      hc_update_up(vc);

    }

}


PUBLIC void
vrna_hc_free(hard_constraintT *hc){

  if(hc){
    if(hc->matrix)  free(hc->matrix);
    if(hc->up_ext)  free(hc->up_ext);
    if(hc->up_hp)   free(hc->up_hp);
    if(hc->up_int)  free(hc->up_int);
    if(hc->up_ml)   free(hc->up_ml);
    free(hc);
  }
}


/* And now, for something completely different...
 *
 * ... the soft constraints section follows below
 */

PUBLIC int
parse_soft_constraints_file(const char *file_name, int length, double default_value, char *sequence, double *values)
{
  FILE *fp;
  char *line;
  int i;
  int count = 0;

  if(!(fp = fopen(file_name, "r"))){
    warn_user("SHAPE data file could not be opened");
    return 0;
  }

  for (i = 0; i < length; ++i)
  {
    sequence[i] = 'N';
    values[i + 1] = default_value;
  }

  sequence[length] = '\0';

  while((line=get_line(fp))){
    int position;
    unsigned char nucleotide = 'N';
    double reactivity = default_value;
    char *second_entry = 0;
    char *third_entry = 0;
    char *c;

    if(sscanf(line, "%d", &position) != 1)
    {
      free(line);
      continue;
    }

    if(position <= 0 || position > length)
    {
      warn_user("Provided SHAPE data outside of sequence scope");
      fclose(fp);
      free(line);
      return 0;
    }

    for(c = line + 1; *c; ++c){
      if(isspace(*(c-1)) && !isspace(*c)) {
        if(!second_entry){
          second_entry = c;
        }else{
          third_entry = c;
          break;
        }
      }
    }

    if(second_entry){
      if(third_entry){
        sscanf(second_entry, "%c", &nucleotide);
        sscanf(third_entry, "%lf", &reactivity);
      }else if(sscanf(second_entry, "%lf", &reactivity) != 1)
        sscanf(second_entry, "%c", &nucleotide);
    }

    sequence[position-1] = nucleotide;
    values[position] = reactivity;
    ++count;

    free(line);
  }

  fclose(fp);

  if(!count)
  {
      warn_user("SHAPE data file is empty");
      return 0;
  }

  return 1;
}

PRIVATE void
parse_parameter_string(char *string, char c1, char c2, float *v1, float *v2)
{
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
      warn_user(warning);

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
        warn_user(warning);
    }
  }
}

PUBLIC int
convert_shape_reactivities_to_probabilities(const char *shape_conversion, double *values, int length, double default_value)
{
  int *indices;
  int i, j;
  int index;
  int ret = 1;

  if(!shape_conversion || !(*shape_conversion) || length <= 0)
    return 0;

  if(*shape_conversion == 'S')
    return 1;

  indices = space(sizeof(int) * (length + 1));
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

    parse_parameter_string(shape_conversion + 1, 's', 'i', &slope, &intercept);

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
vrna_sc_add(vrna_fold_compound *vc,
          const double *constraints,
          unsigned int options){

  unsigned int      n;
  soft_constraintT  *sc;

  if(vc){
    sc                    = (soft_constraintT *)space(sizeof(soft_constraintT));
    sc->constraints       = NULL;
    sc->free_energies     = NULL;
    sc->en_basepair       = NULL;
    sc->en_stack          = NULL;
    sc->exp_en_stack      = NULL;
    sc->boltzmann_factors = NULL;
    sc->exp_en_basepair   = NULL;
    sc->f                 = NULL;
    sc->exp_f             = NULL;
    sc->data              = NULL;
    n                     = vc->length;

    if(vc->sc)
      vrna_sc_destroy(vc->sc);
    vc->sc  = sc;

    if(constraints){
      /*  copy the provided constraints to the data structure.
          Here, we just assume that the softconstraints are already free
          energy contributions per nucleotide.
          We can also apply any function to the data in sc->constraints
          at this point if we want to...
      */
      sc->constraints = (double *)space(sizeof(double) * (n + 1));
      memcpy((void *)(sc->constraints), (const void *)constraints, sizeof(double) * (n+1));

      if(options & VRNA_CONSTRAINT_SOFT_UP)
        vrna_sc_add_up(vc, constraints, options);

    }
  }
}

PUBLIC void
vrna_sc_add_ali(vrna_fold_compound *vc,
                const double **constraints,
                unsigned int options){

  unsigned int      n, s;
  soft_constraintT  *sc;

  if(vc && (vc->type == VRNA_VC_TYPE_ALIGNMENT)){
    if(!vc->sc)
      vc->scs  = (soft_constraintT **)space(sizeof(soft_constraintT*) * (vc->n_seq + 1));
    for(s = 0; s < vc->n_seq; s++){
      sc                    = (soft_constraintT *)space(sizeof(soft_constraintT));
      sc->constraints       = NULL;
      sc->free_energies     = NULL;
      sc->en_basepair       = NULL;
      sc->en_stack          = NULL;
      sc->exp_en_stack      = NULL;
      sc->boltzmann_factors = NULL;
      sc->exp_en_basepair   = NULL;
      sc->f                 = NULL;
      sc->exp_f             = NULL;
      sc->data              = NULL;

      if(vc->scs[s])
        vrna_sc_destroy(vc->scs[s]);
      vc->scs[s]  = sc;
    }

#if 0
    n                     = vc->length;
    if(constraints){
      /*  copy the provided constraints to the data structure.
          Here, we just assume that the softconstraints are already free
          energy contributions per nucleotide.
          We can also apply any function to the data in sc->constraints
          at this point if we want to...
      */
      sc->constraints = (double *)space(sizeof(double) * (n + 1));
      memcpy((void *)(sc->constraints), (const void *)constraints, sizeof(double) * (n+1));

      if(options & VRNA_CONSTRAINT_SOFT_UP)
        vrna_sc_add_up(vc, constraints, options);

    }
#endif
  }
}

PUBLIC void
vrna_sc_add_bp(vrna_fold_compound *vc,
                        const double **constraints,
                        unsigned int options){
                        

  if(options & VRNA_CONSTRAINT_SOFT_MFE)
    vrna_sc_add_bp_mfe(vc, constraints, options);

  if(options & VRNA_CONSTRAINT_SOFT_PF)
    vrna_sc_add_bp_pf(vc, constraints, options);
}

PUBLIC int
vrna_sc_add_deigan_ali( vrna_fold_compound *vc,
                         const char **shape_files,
                         const int *shape_file_association,
                         double m,
                         double b,
                         unsigned int options){

  float   reactivity, *reactivities, e1;
  char    *line, nucleotide, *sequence;
  int     s, i, p, r, position, *pseudo_energies, n_seq;

  n_seq = vc->n_seq;

  if(!vc->sc)
    vrna_sc_add_ali(vc, NULL, options);

  for(s = 0; shape_file_association[s] != -1; s++){
    if(shape_file_association[s] > n_seq)
      warn_user("SHAPE file association exceeds sequence number in alignment");

    /* read the shape file */
    FILE *fp;
    if(!(fp = fopen(shape_files[s], "r"))){
      fprintf(stderr, "WARNING: SHAPE data file %d could not be opened. No shape data will be used.\n", s);
    } else {

      reactivities  = (float *)space(sizeof(float) * (vc->length + 1));
      sequence      = (char *)space(sizeof(char) * (vc->length + 1));

      for(i = 1; i <= vc->length; i++)
        reactivities[i] = -1.;

      while((line=get_line(fp))){
        r = sscanf(line, "%d %c %f", &position, &nucleotide, &reactivity);
        if(r){
          if((position <= 0) || (position > vc->length))
            nrerror("provided shape data outside of sequence scope");

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
        pseudo_energies = (int *)space(sizeof(int) * (vc->length + 1));
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
        FLT_OR_DBL *exp_pe = (FLT_OR_DBL *)space(sizeof(FLT_OR_DBL) * (vc->length + 1));
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
}

PUBLIC  int
parse_soft_constraints_shape_method(const char *method_string,
                                    char *method,
                                    float *param_1,
                                    float *param_2){

  int r;

  char m;
  const char *params = method_string + 1;

  *param_1 = 0;
  *param_2 = 0;

  if (!method_string || !method_string[0])
    return 0;

  *method = m = method_string[0];

  if (m == 'Z')
  {
    *param_1 = 0.89;
    parse_parameter_string(params, 'b', '\0', param_1, NULL);
  }
  else if (m == 'D')
  {
    *param_1 = 1.8;
    *param_2 = -0.6;
    parse_parameter_string(params, 'm', 'b', param_1, param_2);
  }
  else if (m != 'W')
  {
    *method = 0;
    return 0;
  }

  return 1;
}

PUBLIC void
vrna_sc_add_bp_mfe(vrna_fold_compound *vc,
                            const double **constraints,
                            unsigned int options){

  unsigned int      i, j, n;
  soft_constraintT  *sc;
  int               *idx;

  if(vc && constraints){
    n   = vc->length;
    sc  = vc->sc;

    if(!sc)
      vrna_sc_add(vc, NULL, options | VRNA_CONSTRAINT_SOFT_MFE);
    else{
      if(sc->en_basepair)
        free(sc->en_basepair);
      sc->en_basepair     = (int *)space(sizeof(int) * (((n + 1) * (n + 2)) / 2));

      idx = vc->jindx;
      for(i = 1; i < n; i++)
        for(j=i+1; j<=n; j++)
          sc->en_basepair[idx[j]+i] = (int)(constraints[i][j] * 100.);
    }
  }
}

PUBLIC void
vrna_sc_add_bp_pf( vrna_fold_compound *vc,
                            const double **constraints,
                            unsigned int options){

  unsigned int      i, j, n;
  soft_constraintT  *sc;
  int               *idx;

  if(vc && constraints){
    n   = vc->length;

    if(!vc->sc)
      vrna_sc_add(vc, NULL, options | VRNA_CONSTRAINT_SOFT_PF);

    sc = vc->sc;

    pf_paramT *exp_params = vc->exp_params;
    double    GT          = 0.;
    double    temperature = exp_params->temperature;
    double    kT          = exp_params->kT;
    double    TT          = (temperature+K0)/(Tmeasure);

    if(sc->exp_en_basepair)
      free(sc->exp_en_basepair);
    sc->exp_en_basepair     = (FLT_OR_DBL *)space(sizeof(FLT_OR_DBL) * (((n + 1) * (n + 2)) / 2));

    idx = vc->iindx;
    for(i = 1; i < n; i++)
      for(j=i+1; j<=n; j++){
        GT = constraints[i][j] * TT * 1000.;
        sc->exp_en_basepair[idx[i]-j] = exp( -GT / kT);
      }
  }
}

PUBLIC void
vrna_sc_add_up(vrna_fold_compound *vc,
                        const double *constraints,
                        unsigned int options){

  if(options & VRNA_CONSTRAINT_SOFT_MFE)
    vrna_sc_add_up_mfe(vc, constraints, options);

  if(options & VRNA_CONSTRAINT_SOFT_PF)
    vrna_sc_add_up_pf(vc, constraints, options);
}

PUBLIC void
vrna_sc_add_up_mfe(vrna_fold_compound *vc,
                            const double *constraints,
                            unsigned int options){

  unsigned int i, j, n;
  soft_constraintT  *sc;

  if(vc){
    n   = vc->length;
    sc  = vc->sc;

    if(!sc)
      vrna_sc_add(vc, constraints, options | VRNA_CONSTRAINT_SOFT_UP | VRNA_CONSTRAINT_SOFT_MFE);
    else{
      const double *my_constraints = (constraints) ? constraints : (const double *)sc->constraints;
      if(my_constraints){
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

        sc->free_energies = (int **)space(sizeof(int *) * (n + 2));
        for(i = 0; i <= n; i++)
          sc->free_energies[i] = (int *)space(sizeof(int) * (n - i + 2));

        sc->free_energies[n+1] = NULL;

        for(i = 1; i <= n; i++){
          for(j = 1; j <= (n - i + 1); j++){
            sc->free_energies[i][j] =   sc->free_energies[i][j-1]
                                      + (int)(my_constraints[i+j-1] * 100); /* convert to 10kal/mol */
          }
        }
      }
    }
  }
}

PUBLIC void
vrna_sc_add_up_pf( vrna_fold_compound *vc,
                            const double *constraints,
                            unsigned int options){

  unsigned int i, j, n;
  soft_constraintT  *sc;

  if(vc && constraints){
    n   = vc->length;
    sc  = vc->sc;

    if(!sc)
      vrna_sc_add(vc, constraints, options | VRNA_CONSTRAINT_SOFT_UP | VRNA_CONSTRAINT_SOFT_PF);
    else{
      pf_paramT *exp_params = vc->exp_params;
      double    GT          = 0.;
      double    temperature = exp_params->temperature;
      double    kT          = exp_params->kT;
      double    TT          = (temperature+K0)/(Tmeasure);

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

      sc->boltzmann_factors = (FLT_OR_DBL **)space(sizeof(FLT_OR_DBL *) * (n + 2));
      for(i = 0; i <= n; i++){
        sc->boltzmann_factors[i] = (FLT_OR_DBL *)space(sizeof(FLT_OR_DBL) * (n - i + 2));
        for(j = 0; j < n - i + 2; j++)
          sc->boltzmann_factors[i][j] = 1.;
      }

      sc->boltzmann_factors[n+1] = NULL;

      for(i = 1; i <= n; i++){
        for(j = 1; j <= (n - i + 1); j++){
          GT  = (double)((int)(sc->constraints[i+j-1] * 100)) * TT * 10.; /* convert to cal/mol */
          sc->boltzmann_factors[i][j] =   sc->boltzmann_factors[i][j-1]
                                        * exp( -GT / kT);
        }
      }
    }
  }
}

PUBLIC void
vrna_sc_add_sp(vrna_fold_compound *vc,
                        const double *constraints,
                        unsigned int options){

  if(options & VRNA_CONSTRAINT_SOFT_MFE)
    vrna_sc_add_sp_mfe(vc, constraints, options);

  if(options & VRNA_CONSTRAINT_SOFT_PF)
    vrna_sc_add_sp_pf(vc, constraints, options);
}

PUBLIC void
vrna_sc_add_sp_mfe(vrna_fold_compound *vc,
                            const double *constraints,
                            unsigned int options){
  int i;

  if(!vc->sc)
    vrna_sc_add(vc, NULL, options);

  if(!vc->sc->en_stack)
    vc->sc->en_stack = (int *)space(sizeof(int) * (vc->length + 1));

  for(i = 1; i <= vc->length; ++i)
    vc->sc->en_stack[i] += (int)(constraints[i] * 100.);
}

PUBLIC void
vrna_sc_add_sp_pf( vrna_fold_compound *vc,
                            const double *constraints,
                            unsigned int options){
  int i;

  if(!vc->sc)
    vrna_sc_add(vc, NULL, options);

  if(!vc->sc->exp_en_stack){
    vc->sc->exp_en_stack = (FLT_OR_DBL *)space(sizeof(FLT_OR_DBL) * (vc->length + 1));
    for(i = 0; i <= vc->length; ++i)
      vc->sc->exp_en_stack[i] = 1;
  }

  for(i = 1; i <= vc->length; ++i)
    vc->sc->exp_en_stack[i] *= exp(-(constraints[i] * 1000.)/ vc->exp_params->kT);
}

PUBLIC void
vrna_sc_remove(vrna_fold_compound *vc){

  if(vc){
    vrna_sc_destroy(vc->sc);
    vc->sc = NULL;
  }
}

PUBLIC void
vrna_sc_destroy(soft_constraintT *sc){

  int i;
  if(sc){
    if(sc->constraints)
      free(sc->constraints);
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

PRIVATE INLINE  void
adjust_ptypes(char *ptype,
              hard_constraintT *hc,
              unsigned int length,
              unsigned int idx_type){

  unsigned  int i,j;
  int           *index;
  char          *matrix;

  matrix = hc->matrix;

  if(idx_type){
    index = get_iindx(length);
    for(i = 1; i < length; i++)
      for(j = i + 1; j <= length; j++)
        if(matrix[index[i] - j])
          if(!ptype[index[i] - j])
            ptype[index[i] - j] = 7; /* set to non-canonical pair */

  } else {
    index = get_indx(length);
    for(i = 1; i < length; i++)
      for(j = i + 1; j <= length; j++)
        if(matrix[index[j] + i])
          if(!ptype[index[j] + i])
            ptype[index[j] + i] = 7; /* set to non-canonical pair */

  }
  free(index);
}
  
