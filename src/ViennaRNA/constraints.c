/* constraints handling */

#include <config.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <string.h>
#include <limits.h>

#include "ViennaRNA/energy_par.h"
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
              int *index,
              unsigned int idx_type);

PRIVATE INLINE  void
hc_must_pair( unsigned int i,
              char c_option,
              char *hc,
              int *index,
              unsigned int idx_type);

PRIVATE INLINE  void
hc_pairs_upstream(unsigned int i,
                  char c_option,
                  char *hc,
                  unsigned int min_loop_size,
                  int *index,
                  unsigned int idx_type);

PRIVATE INLINE  void
hc_pairs_downstream(unsigned int i,
                    char c_option,
                    char *hc,
                    unsigned int length,
                    unsigned int min_loop_size,
                    int *index,
                    unsigned int idx_type);

PRIVATE INLINE  void
hc_allow_pair(unsigned int i,
              unsigned int j,
              char c_option,
              char *hc,
              int *index,
              unsigned int idx_type);

PRIVATE INLINE  void
hc_weak_enforce_pair( unsigned int i,
                      unsigned int j,
                      char c_option,
                      char *hc,
                      unsigned int length,
                      unsigned int min_loop_size,
                      int *index,
                      unsigned int idx_type);

PRIVATE INLINE  void
hc_enforce_pair(unsigned int i,
                unsigned int j,
                char c_option,
                char *hc,
                unsigned int length,
                unsigned int min_loop_size,
                int *index,
                unsigned int idx_type);

PRIVATE INLINE  void
hc_intramolecular_only( unsigned int i,
                        char c_option,
                        char *hc,
                        unsigned int length,
                        unsigned int min_loop_size,
                        int cut,
                        int *index,
                        unsigned int idx_type);

PRIVATE INLINE  void
hc_intermolecular_only( unsigned int i,
                        char c_option,
                        char *hc,
                        unsigned int length,
                        unsigned int min_loop_size,
                        int cut,
                        int *index,
                        unsigned int idx_type);

PRIVATE INLINE  void
adjust_ptypes(char *ptype,
              hard_constraintT *hc,
              unsigned int length,
              unsigned int idx_type);


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
                    | VRNA_CONSTRAINT_RND_BRACK
                    | ((idx_type) ? VRNA_CONSTRAINT_IINDX : (unsigned int)0));
}

PUBLIC void
apply_DB_constraint(const char *constraint,
                    char *hc,
                    unsigned int length,
                    unsigned int min_loop_size,
                    int cut,
                    unsigned int options){

  int n,i,j,k,l;
  int hx, *stack;
  char type;
  int *index;
  char c_option;
  unsigned int idx_type;

  if(constraint == NULL) return;

  n         = (int)strlen(constraint);
  stack     = (int *) space(sizeof(int)*(n+1));
  idx_type  = options & VRNA_CONSTRAINT_IINDX;
  index     = (idx_type) ? get_iindx(length) : get_indx(length);
  c_option  =   IN_EXT_LOOP
              | IN_HP_LOOP
              | IN_INT_LOOP
              | IN_INT_LOOP_ENC
              | IN_MB_LOOP
              | IN_MB_LOOP_ENC;

  for(hx=0, j=1; j<=n; j++) {
    switch (constraint[j-1]) {
       /* can't pair */
       case 'x':  if(options & VRNA_CONSTRAINT_X){
                    hc_cant_pair(j, c_option, hc, length, min_loop_size, index, idx_type);
                  }
                  break;

      /* must pair, i.e. may not be unpaired */
      case '|':   if(options & VRNA_CONSTRAINT_PIPE){
                    hc_must_pair(j, c_option, hc, index, idx_type);
                  }

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
                    hc_weak_enforce_pair(i, j, c_option, hc, length, min_loop_size, index, idx_type);
                  }
                  break;

      /* pairs upstream */
      case '<':   if(options & VRNA_CONSTRAINT_ANG_BRACK){
                    hc_pairs_upstream(j, c_option, hc, min_loop_size, index, idx_type);
                  }
                  break;

      /* pairs downstream */
      case '>':   if(options & VRNA_CONSTRAINT_ANG_BRACK){
                    hc_pairs_downstream(j, c_option, hc, length, min_loop_size, index, idx_type);
                  }
                  break;

      /* only intramolecular basepairing */
      case 'l':   if(options & VRNA_CONSTRAINT_INTRAMOLECULAR){
                    hc_intramolecular_only(j, c_option, hc, length, min_loop_size, cut, index, idx_type);
                  }
                  break;

      /* only intermolecular bp */
      case 'e':   if(options & VRNA_CONSTRAINT_INTERMOLECULAR){
                    hc_intermolecular_only(j, c_option, hc, length, min_loop_size, cut, index, idx_type);
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
                        int *index,
                        unsigned int idx_type){

  unsigned int l;

  if(cut > 1){
    if(idx_type){
      if(i < cut)
        for(l = MAX2(i+min_loop_size, cut); l <= length; l++)
          hc[index[i] - l] &= ~c_option;
      else
        for(l = 1; l < MIN2(cut, i-min_loop_size); l++)
          hc[index[l] - i] &= ~c_option;
    } else {
      if(i < cut)
        for(l = MAX2(i+min_loop_size, cut); l <= length; l++)
          hc[index[l] + i] &= ~c_option;
      else
        for(l = 1; l < MIN2(cut, i-min_loop_size); l++)
          hc[index[i] + l] &= ~c_option;
    }
  }
}

PRIVATE INLINE  void
hc_intermolecular_only( unsigned int i,
                        char c_option,
                        char *hc,
                        unsigned int length,
                        unsigned int min_loop_size,
                        int cut,
                        int *index,
                        unsigned int idx_type){

  unsigned int l;

  if(cut > 1){
    if(idx_type){
      if(i < cut){
        for(l = 1; l < i; l++)
          hc[index[l] - i] &= ~c_option;
        for(l = i + 1; l < cut; l++)
          hc[index[i] - l] &= ~c_option;
      } else {
        for(l = cut; l < i; l++)
          hc[index[l] - i] &= ~c_option;
        for(l = i + 1; l <= length; l++)
          hc[index[i] - l] &= ~c_option;
      }
    } else {
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
}

PRIVATE INLINE  void
hc_cant_pair( unsigned int i,
              char c_option,
              char *hc,
              unsigned int length,
              unsigned int min_loop_size,
              int *index,
              unsigned int idx_type){

  hc_pairs_upstream(i, c_option, hc, min_loop_size, index, idx_type);
  hc_pairs_downstream(i, c_option, hc, length, min_loop_size, index, idx_type);
}

PRIVATE INLINE  void
hc_must_pair( unsigned int i,
              char c_option,
              char *hc,
              int *index,
              unsigned int idx_type){

  if(idx_type){
    hc[index[i]-i] &= ~c_option;
  } else {
    hc[index[i]+i] &= ~c_option;
  }
}

PRIVATE INLINE  void
hc_pairs_upstream(unsigned int i,
                  char c_option,
                  char *hc,
                  unsigned int min_loop_size,
                  int *index,
                  unsigned int idx_type){

  unsigned int l;

  if(min_loop_size < i){
    if(idx_type){
      for(l = 1; l < i - min_loop_size; l++){
        hc[index[l] - i] &= ~c_option;
      }
    } else {
      for(l = 1; l < i - min_loop_size; l++){
        hc[index[i] + l] &= ~c_option;
      }
    }
  }
}

PRIVATE INLINE  void
hc_pairs_downstream(unsigned int i,
                    char c_option,
                    char *hc,
                    unsigned int length,
                    unsigned int min_loop_size,
                    int *index,
                    unsigned int idx_type){

  unsigned int l;

  if(idx_type){
    for(l = i + min_loop_size + 1; l <= length; l++){
      hc[index[i] - l] &= ~c_option;
    }
  } else {
    for(l = i + min_loop_size + 1; l <= length; l++){
      hc[index[l] + i] &= ~c_option;
    }
  }
}

PRIVATE INLINE  void
hc_allow_pair(unsigned int i,
              unsigned int j,
              char c_option,
              char *hc,
              int *index,
              unsigned int idx_type){

  if(idx_type){
    hc[index[i] - j] |= c_option;
  } else {
    hc[index[j] + i] |= c_option;
  }
}

PRIVATE INLINE  void
hc_weak_enforce_pair( unsigned int i,
                      unsigned int j,
                      char c_option,
                      char *hc,
                      unsigned int length,
                      unsigned int min_loop_size,
                      int *index,
                      unsigned int idx_type){

  unsigned int k, l;

  /* don't allow pairs (k,i) 1 <= k < i */
  hc_pairs_upstream(i, c_option, hc, min_loop_size, index, idx_type);
  /* don't allow pairs (i,k) i < k <= n */ 
  hc_pairs_downstream(i, c_option, hc, length, min_loop_size, index, idx_type);
  /* don't allow pairs (k,j) 1 <= k < j */
  hc_pairs_upstream(j, c_option, hc, min_loop_size, index, idx_type);
  /* don't allow pairs (j,k) j < k <= n */ 
  hc_pairs_downstream(j, c_option, hc, length, min_loop_size, index, idx_type);

  if(idx_type){
    /* don't allow pairs i < k < j < l */
    for(k = i+1; k < j; k++)
      for(l = j+1; l <= length; l++)
        hc[index[k] - l] = 0;
    /* don't allow pairs k<i<l<j */
    for(k = 1; k < i; k++)
      for(l = i+1; l < j; l++)
        hc[index[k] - l] = 0;
    /* allow base pair (i,j) */
    hc[index[i] - j] |= c_option;
  } else {
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
}

PRIVATE INLINE  void
hc_enforce_pair(unsigned int i,
                unsigned int j,
                char c_option,
                char *hc,
                unsigned int length,
                unsigned int min_loop_size,
                int *index,
                unsigned int idx_type){

  hc_weak_enforce_pair( i,
                        j,
                        c_option,
                        hc,
                        length,
                        min_loop_size,
                        index,
                        idx_type);

  if(idx_type){
    /* forbid i and j to be unpaired */
    hc[index[i] - i] = 0;
    hc[index[j] - j] = 0;
  } else {
    /* forbid i and j to be unpaired */
    hc[index[i] + i] = 0;
    hc[index[j] + j] = 0;
  }
}

PUBLIC void
getConstraint(char **cstruc,
              const char **lines,
              unsigned int option){

  int r, i, j, l, cl, stop;
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



PUBLIC  hard_constraintT  *
get_hard_constraints( const char *sequence,
                      const char *constraint,
                      model_detailsT  *md,
                      unsigned int min_loop_size,
                      unsigned int options){

  unsigned int      i, j, ij, n;
  int               *idx, max_span;
  hard_constraintT  *hc;
  short             *S;
  char              *tmp_sequence;

  tmp_sequence  = (options & VRNA_CONSTRAINT_UNGAP) ? get_ungapped_sequence(sequence) : strdup(sequence);
  n             = strlen(tmp_sequence);
  hc            = (hard_constraintT *)space(sizeof(hard_constraintT));
  idx           = (options & VRNA_CONSTRAINT_IINDX) ? get_iindx(n) : get_indx(n);
  S             = get_sequence_encoding(tmp_sequence, 0, md);
  max_span      = md->max_bp_span;
  if((max_span < 5) || (max_span > n))
    max_span  = n;

  /* allocate memory for the hard constraints data structure */
  hc->matrix  = (char *)space(sizeof(char)*((n*(n+1))/2+2));
  hc->up_ext  = (int *)space(sizeof(int)*(n+2));
  hc->up_hp   = (int *)space(sizeof(int)*(n+2));
  hc->up_int  = (int *)space(sizeof(int)*(n+2));
  hc->up_ml   = (int *)space(sizeof(int)*(n+2));

  /* ####################### */
  /* prefill default values  */
  /* ####################### */

  if(options & VRNA_CONSTRAINT_IINDX){ /* iindx[i]-j */
    /* 1. unpaired nucleotides are allowed in all contexts */
    for(i = 1; i <= n; i++)
      hc->matrix[idx[i] - i]  =   IN_EXT_LOOP
                                | IN_HP_LOOP
                                | IN_INT_LOOP
                                | IN_MB_LOOP;

    /* 2. all canonical base pairs are allowed in all contexts */
    for(i = 1; i < n; i++){
      ij = idx[i]-i-min_loop_size-1;
      for(j=i+min_loop_size+1; j <= n; j++,ij--){
        if((j-i+1) > max_span)
          break;
        hc->matrix[ij] = md->pair[S[i]][S[j]] ? IN_EXT_LOOP
                                                | IN_HP_LOOP
                                                | IN_INT_LOOP
                                                | IN_MB_LOOP
                                                | IN_INT_LOOP_ENC
                                                | IN_MB_LOOP_ENC
                                              : (char)0;
      }
    }
  } else { /* indx[j] + i */
    /* 1. unpaired nucleotides are allowed in all contexts */
    for(i = 1; i <= n; i++)
      hc->matrix[idx[i] + i]  =   IN_EXT_LOOP
                                | IN_HP_LOOP
                                | IN_INT_LOOP
                                | IN_MB_LOOP;

    /* 2. all canonical base pairs are allowed in all contexts */
    for(j = n; j > min_loop_size + 1; j--){
      ij = idx[j]+1;
      for(i=1; i < j - min_loop_size; i++, ij++)
        if((j-i+1) > max_span){
          hc->matrix[ij] = (char)0;
        } else {
          hc->matrix[ij] = md->pair[S[i]][S[j]] ? IN_EXT_LOOP
                                                | IN_HP_LOOP
                                                | IN_INT_LOOP
                                                | IN_MB_LOOP
                                                | IN_INT_LOOP_ENC
                                                | IN_MB_LOOP_ENC
                                              : (char)0;
        }
    }
  }


  /* ############################### */
  /* apply user supplied constraints */
  /* ############################### */

  if(options & VRNA_CONSTRAINT_DB){ /* constraints from dot-bracket notation */
    printf("reading constraints from dot-bracket\n");
    apply_DB_constraint(constraint,
                        hc->matrix,
                        n,
                        min_loop_size,
                        -1,
                        options);
  }

  if(options & VRNA_CONSTRAINT_FILE){ /* constraints from file */
    printf("reading constraints from file\n");
  }


  /* ####################### */
  /* do some post processing */
  /* ####################### */

  if(options & VRNA_CONSTRAINT_IINDX){ /* iindx[i]-j */
    /* compute the number of nucleotides available to be unpaired */

    for(hc->up_ext[n+1] = 0, i = n; i > 0; i--) /* unpaired stretch in exterior loop */
      hc->up_ext[i] = (hc->matrix[idx[i]-i] & IN_EXT_LOOP) ? 1 + hc->up_ext[i+1] : 0;

    for(hc->up_hp[n+1] = 0, i = n; i > 0; i--)  /* unpaired stretch in hairpin loop */
      hc->up_hp[i] = (hc->matrix[idx[i]-i] & IN_HP_LOOP) ? 1 + hc->up_hp[i+1] : 0;

    for(hc->up_int[n+1] = 0, i = n; i > 0; i--) /* unpaired stretch in interior loop */
      hc->up_int[i] = (hc->matrix[idx[i]-i] & IN_INT_LOOP) ? 1 + hc->up_int[i+1] : 0;

    for(hc->up_ml[n+1] = 0, i = n; i > 0; i--)  /* unpaired stretch in multibranch loop */
      hc->up_ml[i] = (hc->matrix[idx[i]-i] & IN_MB_LOOP) ? 1 + hc->up_ml[i+1] : 0;

  } else { /* indx[j] + i */
    /* compute the number of nucleotides available to be unpaired */

    for(hc->up_ext[n+1] = 0, i = n; i > 0; i--) /* unpaired stretch in exterior loop */
      hc->up_ext[i] = (hc->matrix[idx[i]+i] & IN_EXT_LOOP) ? 1 + hc->up_ext[i+1] : 0;

    for(hc->up_hp[n+1] = 0, i = n; i > 0; i--){  /* unpaired stretch in hairpin loop */
      hc->up_hp[i] = (hc->matrix[idx[i]+i] & IN_HP_LOOP) ? 1 + hc->up_hp[i+1] : 0;
    }
    for(hc->up_int[n+1] = 0, i = n; i > 0; i--) /* unpaired stretch in interior loop */
      hc->up_int[i] = (hc->matrix[idx[i]+i] & IN_INT_LOOP) ? 1 + hc->up_int[i+1] : 0;

    for(hc->up_ml[n+1] = 0, i = n; i > 0; i--)  /* unpaired stretch in multibranch loop */
      hc->up_ml[i] = (hc->matrix[idx[i]+i] & IN_MB_LOOP) ? 1 + hc->up_ml[i+1] : 0;

  }

  free(tmp_sequence);
  free(idx);
  free(S);

  return hc;
}

PUBLIC void
destroy_hard_constraints(hard_constraintT *hc){

  if(hc){
    if(hc->matrix)  free(hc->matrix);
    if(hc->up_ext)  free(hc->up_ext);
    if(hc->up_hp)   free(hc->up_hp);
    if(hc->up_int)  free(hc->up_int);
    if(hc->up_ml)   free(hc->up_ml);
    free(hc);
  }
}


PUBLIC void
add_soft_constraints( vrna_fold_compound *vc,
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
    sc->boltzmann_factors = NULL;
    sc->exp_en_basepair   = NULL;
    sc->f                 = NULL;
    sc->exp_f             = NULL;
    sc->data              = NULL;
    n                     = vc->length;

    if(vc->sc)
      destroy_soft_constraints(vc->sc);
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
        add_soft_constraints_up(vc, NULL, options);

    }
  }
}

PUBLIC void
add_soft_constraints_bp(vrna_fold_compound *vc,
                        const double **constraints,
                        unsigned int options){
                        

  if(options & VRNA_CONSTRAINT_SOFT_MFE)
    add_soft_constraints_bp_mfe(vc, constraints, options);

  if(options & VRNA_CONSTRAINT_SOFT_PF)
    add_soft_constraints_bp_pf(vc, constraints, options);
}


PUBLIC int
add_soft_constraints_mathews( vrna_fold_compound *vc,
                              const char *shape_file,
                              double m,
                              double b,
                              unsigned int options){

  float   reactivity, *reactivities;
  char    *line, nucleotide, *sequence;
  int     i, j, r, position, *pseudo_energies, *idx;

  /* read the shape file */
  FILE *fp;
  if(!(fp = fopen(shape_file, "r"))){
    warn_user("SHAPE data file could not be opened. No shape data will be used.");
    return 0;
  }

  reactivities  = (float *)space(sizeof(float) * (vc->length + 1));
  sequence      = (char *)space(sizeof(char) * (vc->length + 1));

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
  if(strcmp(vc->sequence, sequence))
    warn_user("input sequence differs from sequence provided via SHAPE file!");

  /* convert reactivities to pseudo energies */
  for(i = 1; i <= vc->length; i++){
    if(reactivities[i] < 0)
      reactivities[i] = 0.;
    else
      reactivities[i] = m * log(reactivities[i] + 1.) + b;
  }

  /* create the pseudo energy lookup table for the recursions */
  idx = vc->jindx;
  pseudo_energies = (int *)space(sizeof(int) * (((vc->length + 1) * (vc->length + 2)) / 2));
  for(i = 1; i<vc->length; i++)
    for(j = i + 1; j <= vc->length; j++)
      pseudo_energies[idx[j] + i] = (int)((reactivities[i] + reactivities[j]) * 100.);

  if(vc->sc){
    if(vc->sc->en_stack){
      free(vc->sc->en_stack);
    }
  } else {
    add_soft_constraints(vc, NULL, options);
  }

  vc->sc->en_stack = pseudo_energies;
  free(reactivities);

  return 1; /* success */
}

PUBLIC int
add_soft_constraints_alignment_mathews( vrna_alifold_compound *vc,
                                        const char **shape_files,
                                        const int *shape_file_association,
                                        double m,
                                        double b,
                                        unsigned int options){

  float   reactivity, *reactivities;
  char    *line, nucleotide, *sequence;
  int     s, i, j, r, position, *pseudo_energies, *idx, n_seq;

  n_seq = vc->n_seq;
  vc->sc = (soft_constraintT **)space(sizeof(soft_constraintT *) * (n_seq + 1));
  for(s = 0; s < n_seq; s++){
    soft_constraintT *sc = (soft_constraintT *)space(sizeof(soft_constraintT));
    sc->constraints       = NULL;
    sc->free_energies     = NULL;
    sc->en_basepair       = NULL;
    sc->en_stack          = NULL;
    sc->boltzmann_factors = NULL;
    sc->exp_en_basepair   = NULL;
    sc->f                 = NULL;
    sc->exp_f             = NULL;
    sc->data              = NULL;
    vc->sc[s]             = sc;
  }

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
          reactivities[i] = m * log(reactivities[i] + 1.) + b;
      }

      /* create the pseudo energy lookup table for the recursions */
      idx = vc->jindx;
      pseudo_energies = (int *)space(sizeof(int) * (((vc->length + 1) * (vc->length + 2)) / 2));
      for(i = 1; i<vc->length; i++)
        for(j = i + 1; j <= vc->length; j++)
          pseudo_energies[idx[j] + i] = (int)((reactivities[i] + reactivities[j]) * 100.);

      vc->sc[shape_file_association[s]]->en_stack = pseudo_energies;
      free(reactivities);
    }
  }
  return 1; /* success */
}

PUBLIC  int
parse_soft_constraints_shape_method(  const char *method_string,
                                      char *method,
                                      float *param_1,
                                      float *param_2){

  int r;
  float p1;
  float p2;

  if(method_string){
    switch(method_string[0]){
      case  'M':  *method = 'M';
                  if(method_string[1]){
                    r = sscanf(method_string+1, "m%fb%f", &p1, &p2);
                    if(r != 2){
                      r = sscanf(method_string+1, "m%f", &p1);
                      if(!r){
                        p1 = 1.8; /* default value for m */
                        r = sscanf(method_string+1, "b%f", &p2);
                        if(!r){
                          p2 = -0.6; /* default value for b*/
                          warn_user("SHAPE method parameters not recognized! Using default parameters!");
                        }
                      }
                    }
                    *param_1 = p1;
                    *param_2 = p2;
                  }
                  return 1;
    }
  }
  return 0; /* failure */
}

PUBLIC void
add_soft_constraints_bp_mfe(vrna_fold_compound *vc,
                            const double **constraints,
                            unsigned int options){

  unsigned int      i, j, n;
  soft_constraintT  *sc;
  int               *idx;

  if(vc && constraints){
    n   = vc->length;
    sc  = vc->sc;

    if(!sc)
      add_soft_constraints(vc, NULL, options | VRNA_CONSTRAINT_SOFT_MFE);
    else{
      if(sc->en_basepair)
        free(sc->en_basepair);
      sc->en_basepair     = (int *)space(sizeof(int) * (((n + 1) * (n + 2)) / 2));

      idx = vc->jindx;
      for(i = 1; i < n; i++)
        for(j=i+1; j<=n; j++)
          sc->en_basepair[idx[j]+i] = constraints[i][j];
    }
  }
}

PUBLIC void
add_soft_constraints_bp_pf( vrna_fold_compound *vc,
                            const double **constraints,
                            unsigned int options){

  unsigned int      i, j, n;
  soft_constraintT  *sc;
  int               *idx;

  if(vc && constraints){
    n   = vc->length;
    sc  = vc->sc;

    if(!sc)
      add_soft_constraints(vc, NULL, options | VRNA_CONSTRAINT_SOFT_PF);
    else{
      pf_paramT *exp_params = vc->exp_params;
      double    GT          = 0.;
      double    temperature = exp_params->temperature;
      double    BetaScale   = exp_params->alpha;
      double    kT          = exp_params->kT;
      double    pf_scale    = exp_params->pf_scale;
      double    TT          = (temperature+K0)/(Tmeasure);

      if(sc->exp_en_basepair)
        free(sc->exp_en_basepair);
      sc->exp_en_basepair     = (FLT_OR_DBL *)space(sizeof(FLT_OR_DBL) * (((n + 1) * (n + 2)) / 2));

      idx = vc->iindx;
      for(i = 1; i < n; i++)
        for(j=i+1; j<=n; j++){
          GT = constraints[i][j] * TT;
          sc->exp_en_basepair[idx[i]-j] = exp( -GT * 10./ kT);
        }
    }
  }
}

PUBLIC void
add_soft_constraints_up(vrna_fold_compound *vc,
                        const double *constraints,
                        unsigned int options){

  if(options & VRNA_CONSTRAINT_SOFT_MFE)
    add_soft_constraints_up_mfe(vc, constraints, options);

  if(options & VRNA_CONSTRAINT_SOFT_PF)
    add_soft_constraints_up_pf(vc, constraints, options);
}

PUBLIC void
add_soft_constraints_up_mfe(vrna_fold_compound *vc,
                            const double *constraints,
                            unsigned int options){

  unsigned int i, j, n;
  soft_constraintT  *sc;

  if(vc){
    n   = vc->length;
    sc  = vc->sc;

    if(!sc)
      add_soft_constraints(vc, constraints, options | VRNA_CONSTRAINT_SOFT_UP | VRNA_CONSTRAINT_SOFT_MFE);
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
            sc->free_energies[i][j] = sc->free_energies[i][j-1] + (int)(my_constraints[i+j-1]*100);
          }
        }
      }
    }
  }
}

PUBLIC void
add_soft_constraints_up_pf( vrna_fold_compound *vc,
                            const double *constraints,
                            unsigned int options){

  unsigned int i, j, n;
  soft_constraintT  *sc;

  if(vc && constraints){
    n   = vc->length;
    sc  = vc->sc;

    if(!sc)
      add_soft_constraints(vc, constraints, options | VRNA_CONSTRAINT_SOFT_UP | VRNA_CONSTRAINT_SOFT_PF);
    else{
      pf_paramT *exp_params = vc->exp_params;
      double    GT          = 0.;
      double    temperature = exp_params->temperature;
      double    BetaScale   = exp_params->alpha;
      double    kT          = exp_params->kT;
      double    pf_scale    = exp_params->pf_scale;
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
          GT  = (double)((int)(sc->constraints[i+j-1]*100)) * TT;
          sc->boltzmann_factors[i][j] = sc->boltzmann_factors[i][j-1] * exp( -GT*10./kT);
        }
      }
    }
  }
}

PUBLIC void
remove_soft_constraints(vrna_fold_compound *vc){

  if(vc){
    destroy_soft_constraints(vc->sc);
    vc->sc = NULL;
  }
}

PUBLIC void
destroy_soft_constraints(soft_constraintT *sc){

  int     i, *ptr, *ptr2;
  double  *ptr3;
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
  
