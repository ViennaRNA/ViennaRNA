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

PUBLIC  void  print_tty_constraint_full(void){
  print_tty_constraint(VRNA_CONSTRAINT_PIPE | VRNA_CONSTRAINT_DOT | VRNA_CONSTRAINT_X | VRNA_CONSTRAINT_ANG_BRACK | VRNA_CONSTRAINT_RND_BRACK);
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
                              int min_loop_size,
                              unsigned int idx_type){

  int n,i,j,k,l;
  int hx, *stack;
  char type;
  int *index;

  if(constraint == NULL) return;

  n = (int)strlen(constraint);

  stack = (int *) space(sizeof(int)*(n+1));

  if(!idx_type){ /* index allows access in energy matrices at pos (i,j) via index[j]+i */
    index = get_indx(length);

    for(hx=0, j=1; j<=n; j++){
      switch(constraint[j-1]){
        case 'x':   /* can't pair */
                    for (l=1; l<j-min_loop_size; l++) ptype[index[j]+l] = 0;
                    for (l=j+min_loop_size+1; l<=(int)length; l++) ptype[index[l]+j] = 0;
                    break;
        case '(':   stack[hx++]=j;
                    /* fallthrough */
        case '<':   /* pairs upstream */
                    for (l=1; l<j-min_loop_size; l++) ptype[index[j]+l] = 0;
                    break;
        case ')':   if (hx<=0) {
                      fprintf(stderr, "%s\n", constraint);
                      nrerror("unbalanced brackets in constraint");
                    }
                    i = stack[--hx];
                    type = ptype[index[j]+i];
                    for (k=i+1; k<=(int)length; k++) ptype[index[k]+i] = 0;
                    /* don't allow pairs i<k<j<l */
                    for (l=j; l<=(int)length; l++)
                      for (k=i+1; k<=j; k++) ptype[index[l]+k] = 0;
                    /* don't allow pairs k<i<l<j */
                    for (l=i; l<=j; l++)
                      for (k=1; k<=i; k++) ptype[index[l]+k] = 0;
                    for (k=1; k<j; k++) ptype[index[j]+k] = 0;
                    ptype[index[j]+i] = (type==0) ? 7 : type;
                    /* fallthrough */
        case '>':   /* pairs downstream */
                    for (l=j+min_loop_size+1; l<=(int)length; l++) ptype[index[l]+j] = 0;
                    break;
      }
    }
  }
  else{ /* index allows access in energy matrices at pos (i,j) via index[i]-j */
    index = get_iindx(length);

    for(hx=0, j=1; j<=n; j++) {
      switch (constraint[j-1]) {
        case 'x':   /* can't pair */
                    for (l=1; l<j-min_loop_size; l++) ptype[index[l]-j] = 0;
                    for (l=j+min_loop_size+1; l<=(int)length; l++) ptype[index[j]-l] = 0;
                    break;
        case '(':   stack[hx++]=j;
                    /* fallthrough */
        case '<':   /* pairs upstream */
                    for (l=1; l<j-min_loop_size; l++) ptype[index[l]-j] = 0;
                    break;
        case ')':   if (hx<=0) {
                      fprintf(stderr, "%s\n", constraint);
                      nrerror("unbalanced brackets in constraints");
                    }
                    i = stack[--hx];
                    type = ptype[index[i]-j];
                    /* don't allow pairs i<k<j<l */
                    for (k=i; k<=j; k++)
                      for (l=j; l<=(int)length; l++) ptype[index[k]-l] = 0;
                    /* don't allow pairs k<i<l<j */
                    for (k=1; k<=i; k++)
                      for (l=i; l<=j; l++) ptype[index[k]-l] = 0;
                    ptype[index[i]-j] = (type==0) ? 7 : type;
                    /* fallthrough */
        case '>':   /* pairs downstream */
                    for (l=j+min_loop_size+1; l<=(int)length; l++) ptype[index[j]-l] = 0;
                    break;
      }
    }
  }
  if (hx!=0) {
    fprintf(stderr, "%s\n", constraint);
    nrerror("unbalanced brackets in constraint string");
  }
  free(index);
  free(stack);
}


PUBLIC void getConstraint(char **cstruc, const char **lines, unsigned int option){
  int r, i, j, l, cl, stop;
  char *c, *ptr;
  if(lines){
    if(option & VRNA_CONSTRAINT_ALL)
      option |= VRNA_CONSTRAINT_PIPE | VRNA_CONSTRAINT_ANG_BRACK | VRNA_CONSTRAINT_RND_BRACK | VRNA_CONSTRAINT_X;

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
  hc->up_ext  = (int *)space(sizeof(int)*(n+1));
  hc->up_hp   = (int *)space(sizeof(int)*(n+2));
  hc->up_int  = (int *)space(sizeof(int)*(n+2));
  hc->up_ml   = (int *)space(sizeof(int)*(n+2));

  /* prefill defaults */

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


  /* apply user supplied constraints */

  if(options & VRNA_CONSTRAINT_DB){ /* constraints from dot-bracket notation */
    printf("reading constraints from dot-bracket\n");
    constrain_ptypes( constraint,
                      n,
                      hc->matrix,
                      min_loop_size,
                      options & VRNA_CONSTRAINT_IINDX);
  }

  if(options & VRNA_CONSTRAINT_FILE){ /* constraints from file */
    printf("reading constraints from file\n");
  }

  /* compute the number of nucleotides available to be unpaired */

  for(i = n; i > 0; i--)  /* unpaired stretch in exterior loop */
    hc->up_ext[i] = n - i + 1;

  for(i = n; i > 0; i--)  /* unpaired stretch in hairpin loop */
    hc->up_hp[i] = n - i + 1;

  for(i = n; i > 0; i--)  /* unpaired stretch in interior loop */
    hc->up_int[i] = n - i + 1;

  for(i = n; i > 0; i--)  /* unpaired stretch in multibranch loop */
    hc->up_ml[i] = n - i + 1;

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
