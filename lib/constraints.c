/* constraints handling */

#include <config.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <string.h>
#include <limits.h>

#include "fold_vars.h"
#include "utils.h"
#include "constraints.h"


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
PRIVATE INLINE  void hc_cant_pair(unsigned int i,
                                  char c_option,
                                  char *hc,
                                  unsigned int length,
                                  unsigned int min_loop_size,
                                  int *index,
                                  unsigned int idx_type);

PRIVATE INLINE  void hc_must_pair(unsigned int i,
                                  char c_option,
                                  char *hc,
                                  int *index,
                                  unsigned int idx_type);

PRIVATE INLINE  void hc_pairs_upstream( unsigned int i,
                                        char c_option,
                                        char *hc,
                                        unsigned int min_loop_size,
                                        int *index,
                                        unsigned int idx_type);

PRIVATE INLINE  void hc_pairs_downstream( unsigned int i,
                                          char c_option,
                                          char *hc,
                                          unsigned int length,
                                          unsigned int min_loop_size,
                                          int *index,
                                          unsigned int idx_type);

PRIVATE INLINE  void hc_allow_pair( unsigned int i,
                                    unsigned int j,
                                    char c_option,
                                    char *hc,
                                    int *index,
                                    unsigned int idx_type);

PRIVATE INLINE  void hc_weak_enforce_pair(unsigned int i,
                                          unsigned int j,
                                          char c_option,
                                          char *hc,
                                          unsigned int length,
                                          unsigned int min_loop_size,
                                          int *index,
                                          unsigned int idx_type);

PRIVATE INLINE  void hc_enforce_pair(unsigned int i,
                                      unsigned int j,
                                      char c_option,
                                      char *hc,
                                      unsigned int length,
                                      unsigned int min_loop_size,
                                      int *index,
                                      unsigned int idx_type);

PRIVATE INLINE  void hc_intramolecular_only(unsigned int i,
                                            char c_option,
                                            char *hc,
                                            unsigned int length,
                                            unsigned int min_loop_size,
                                            int cut,
                                            int *index,
                                            unsigned int idx_type);

PRIVATE INLINE  void hc_intermolecular_only(unsigned int i,
                                            char c_option,
                                            char *hc,
                                            unsigned int length,
                                            unsigned int min_loop_size,
                                            int cut,
                                            int *index,
                                            unsigned int idx_type);

PRIVATE INLINE  void adjust_ptypes( char *ptype,
                                    hard_constraintT *hc,
                                    unsigned int length,
                                    unsigned int idx_type);


/*
#################################
# BEGIN OF FUNCTION DEFINITIONS #
#################################
*/
PUBLIC  void  print_tty_constraint_full(void){
  print_tty_constraint(   VRNA_CONSTRAINT_PIPE
                        | VRNA_CONSTRAINT_DOT
                        | VRNA_CONSTRAINT_X
                        | VRNA_CONSTRAINT_ANG_BRACK
                        | VRNA_CONSTRAINT_RND_BRACK);
}

PUBLIC  void  print_tty_constraint(unsigned int option){
  if(!(option & VRNA_CONSTRAINT_NO_HEADER)) printf("Input structure constraints using the following notation:\n");
  if(option & VRNA_CONSTRAINT_PIPE)       printf("| : paired with another base\n");
  if(option & VRNA_CONSTRAINT_DOT)        printf(". : no constraint at all\n");
  if(option & VRNA_CONSTRAINT_X)          printf("x : base must not pair\n");
  if(option & VRNA_CONSTRAINT_ANG_BRACK)  printf("< : base i is paired with a base j<i\n> : base i is paired with a base j>i\n");
  if(option & VRNA_CONSTRAINT_RND_BRACK)  printf("matching brackets ( ): base i pairs base j\n");
}


PUBLIC void constrain_ptypes( const char *constraint,
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

PUBLIC void apply_DB_constraint(const char *constraint,
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

PRIVATE INLINE  void hc_intramolecular_only(unsigned int i,
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

PRIVATE INLINE  void hc_intermolecular_only(unsigned int i,
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

PRIVATE INLINE  void hc_cant_pair(unsigned int i,
                                  char c_option,
                                  char *hc,
                                  unsigned int length,
                                  unsigned int min_loop_size,
                                  int *index,
                                  unsigned int idx_type){

  hc_pairs_upstream(i, c_option, hc, min_loop_size, index, idx_type);
  hc_pairs_downstream(i, c_option, hc, length, min_loop_size, index, idx_type);
}

PRIVATE INLINE  void hc_must_pair(unsigned int i,
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

PRIVATE INLINE  void hc_pairs_upstream( unsigned int i,
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

PRIVATE INLINE  void hc_pairs_downstream( unsigned int i,
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

PRIVATE INLINE  void hc_allow_pair( unsigned int i,
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

PRIVATE INLINE  void hc_weak_enforce_pair(unsigned int i,
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

PRIVATE INLINE  void hc_enforce_pair(unsigned int i,
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

PUBLIC void getConstraint(char **cstruc, const char **lines, unsigned int option){
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

PUBLIC  hard_constraintT  *get_hard_constraints(  const char *constraint,
                                                  unsigned int n,
                                                  char *ptype,
                                                  unsigned int min_loop_size,
                                                  unsigned int options){

  unsigned int      i, j, ij;
  int               *idx;
  hard_constraintT  *hc;

  hc  = (hard_constraintT *)space(sizeof(hard_constraintT));
  idx = (options & VRNA_CONSTRAINT_IINDX) ? get_iindx(n) : get_indx(n);

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
      for(j=i+min_loop_size+1; j <= n; j++,ij--)
        hc->matrix[ij] = (ptype[ij])  ?   IN_EXT_LOOP
                                        | IN_HP_LOOP
                                        | IN_INT_LOOP
                                        | IN_MB_LOOP
                                        | IN_INT_LOOP_ENC
                                        | IN_MB_LOOP_ENC
                                      : (char)0;
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
        hc->matrix[ij] = (ptype[ij])  ?   IN_EXT_LOOP
                                        | IN_HP_LOOP
                                        | IN_INT_LOOP
                                        | IN_MB_LOOP
                                        | IN_INT_LOOP_ENC
                                        | IN_MB_LOOP_ENC
                                      : (char)0;
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

  /* adjust the ptype to take non-canonical pairs into account */
  adjust_ptypes(ptype, hc, n, options & VRNA_CONSTRAINT_IINDX);

  free(idx);
  return hc;
}

PUBLIC void destroy_hard_constraints(hard_constraintT *hc){
  if(hc){
    if(hc->matrix)  free(hc->matrix);
    if(hc->up_ext)  free(hc->up_ext);
    if(hc->up_hp)   free(hc->up_hp);
    if(hc->up_int)  free(hc->up_int);
    if(hc->up_ml)   free(hc->up_ml);
    free(hc);
  }
}

PUBLIC soft_constraintT *get_soft_constraints(  const double *constraints,
                                                unsigned int n,
                                                unsigned int options){

  unsigned int i, j;
  soft_constraintT *sc;

  sc                    = (soft_constraintT *)space(sizeof(soft_constraintT));
  sc->constraints       = NULL;
  sc->free_energies     = NULL;
  sc->boltzmann_factors = NULL;

  if(constraints){
    /*  copy the provided constraints to the data structure.
        Here, we just assume that the softconstraints are already free
        energy contributions per nucleotide.
        We can also apply any function to the data in sc->constraints
        at this point if we want to...
    */
    sc->constraints = (double *)space(sizeof(double) * (n + 1));
    memcpy((void *)(sc->constraints), (const void *)constraints, sizeof(double) * (n+1));

    if(options & VRNA_CONSTRAINT_SOFT_MFE){
      /*  allocate memory such that we can access the soft constraint
          energies of a subsequence of length j starting at position i
          via sc->free_energies[i][j]
      */
      sc->free_energies = (int **)space(sizeof(int *) * (n + 2));
      for(i = 0; i <= n; i++){
        sc->free_energies[i] = (int *)space(sizeof(int) * (n - i + 2));
      }

      sc->free_energies[n+1] = NULL;

      for(i = 1; i <= n; i++){
        for(j = 1; j <= (n - i + 1); j++){
          sc->free_energies[i][j] = sc->free_energies[i][j-1] + (int)(sc->constraints[i+j-1]*100);
        }
      }

    }
    if(options & VRNA_CONSTRAINT_SOFT_PF){

    }
  }

  return sc;
}

PUBLIC void destroy_soft_constraints(soft_constraintT *sc){
  int     i, *ptr, *ptr2;
  double  *ptr3;
  if(sc){
    if(sc->constraints)       free(sc->constraints);
    if(sc->free_energies){
      for(i = 0; sc->free_energies[i]; free(sc->free_energies[i++]));
      free(sc->free_energies);
    }
    if(sc->boltzmann_factors){
      for(i = 0; sc->boltzmann_factors[i]; free(sc->boltzmann_factors[i++]));
      free(sc->boltzmann_factors);
    }
    free(sc);
  }
}

PRIVATE INLINE  void adjust_ptypes( char *ptype,
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
  
