/*
    structure_utils.c

    Various functions to convert, parse, encode secondary structures

    c  Ivo L Hofacker, Walter Fontana, Ronny Lorenz
                Vienna RNA package
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "ViennaRNA/fold_vars.h"
#include "ViennaRNA/utils.h"
#include "ViennaRNA/params.h"
#include "ViennaRNA/gquad.h"
#include "ViennaRNA/structure_utils.h"

/*
#################################
# PRIVATE FUNCTION DECLARATIONS #
#################################
*/
PRIVATE plist *
wrap_get_plist( pf_matricesT *matrices,
                int length,
                int *index,
                short *S,
                pf_paramT *pf_params,
                double cut_off);

/*
#################################
# BEGIN OF FUNCTION DEFINITIONS #
#################################
*/

PUBLIC char *
vrna_pack_structure(const char *struc){

  /* 5:1 compression using base 3 encoding */
  int i,j,l,pi;
  unsigned char *packed;

  l = (int) strlen(struc);
  packed = (unsigned char *) space(((l+4)/5+1)*sizeof(unsigned char));

  j=i=pi=0;
  while (i<l) {
    register int p;
    for (p=pi=0; pi<5; pi++) {
      p *= 3;
      switch (struc[i]) {
      case '(':
      case '\0':
        break;
      case '.':
        p++;
        break;
      case ')':
        p += 2;
        break;
      default: nrerror("pack_structure: illegal charcter in structure");
      }
      if (i<l) i++;
    }
    packed[j++] = (unsigned char) (p+1); /* never use 0, so we can use
                                            strcmp()  etc. */
  }
  packed[j] = '\0';      /* for str*() functions */
  return (char *) packed;
}

PUBLIC char *
vrna_unpack_structure(const char *packed){

  /* 5:1 compression using base 3 encoding */
  int i,j,l;
  char *struc;
  unsigned const char *pp;
  char code[3] = {'(', '.', ')'};

  l = (int) strlen(packed);
  pp = (const unsigned char *) packed;
  struc = (char *) space((l*5+1)*sizeof(char));   /* up to 4 byte extra */

  for (i=j=0; i<l; i++) {
    register int p, c, k;

    p = (int) pp[i] - 1;
    for (k=4; k>=0; k--) {
      c = p % 3;
      p /= 3;
      struc[j+k] = code[c];
    }
    j += 5;
  }
  struc[j--] = '\0';
  while (struc[j] == '(') /* strip trailing ( */
    struc[j--] = '\0';

  return struc;
}

PUBLIC short *
vrna_pt_get(const char *structure){

    /* returns array representation of structure.
       table[i] is 0 if unpaired or j if (i.j) pair.  */
   short i,j,hx;
   short length;
   short *stack;
   short *table;

   length = (short) strlen(structure);
   stack = (short *) space(sizeof(short)*(length+1));
   table = (short *) space(sizeof(short)*(length+2));
   table[0] = length;

   for (hx=0, i=1; i<=length; i++) {
      switch (structure[i-1]) {
       case '(':
         stack[hx++]=i;
         break;
       case ')':
         j = stack[--hx];
         if (hx<0) {
            fprintf(stderr, "%s\n", structure);
            nrerror("unbalanced brackets in make_pair_table");
         }
         table[i]=j;
         table[j]=i;
         break;
       default:   /* unpaired base, usually '.' */
         table[i]= 0;
         break;
      }
   }
   if (hx!=0) {
      fprintf(stderr, "%s\n", structure);
      nrerror("unbalanced brackets in make_pair_table");
   }
   free(stack);
   return(table);
}

PUBLIC short *
vrna_pt_pk_get(const char *structure){

   short i,j,hx, hx2;
   short length;
   short *stack;
   short *stack2;
   short *table;

   length = (short) strlen(structure);
   stack  = (short *) space(sizeof(short)*(length+1));
   stack2 = (short *) space(sizeof(short)*(length+1));
   table  = (short *) space(sizeof(short)*(length+2));
   table[0] = length;

   for (hx=0, hx2=0, i=1; i<=length; i++) {
      switch (structure[i-1]) {
       case '(':
         stack[hx++]=i;
         break;
       case ')':
         j = stack[--hx];
         if (hx<0) {
            fprintf(stderr, "%s\n", structure);
            nrerror("unbalanced '()' brackets in make_pair_table_pk");
         }
         table[i]=j;
         table[j]=i;
         break;
       case '[':
         stack2[hx2++]=i;
         break;
       case ']':
         j = stack2[--hx2];
         if (hx2<0) {
            fprintf(stderr, "%s\n", structure);
            nrerror("unbalanced '[]' brackets in make_pair_table_pk");
         }
         table[i]=j;
         table[j]=i;
         break;
       default:   /* unpaired base, usually '.' */
         table[i]= 0;
         break;
      }
   }
   if (hx!=0) {
      fprintf(stderr, "%s\n", structure);
      nrerror("unbalanced '()' brackets in make_pair_table_pk");
   } else if (hx2!=0) {
      fprintf(stderr, "%s\n", structure);
      nrerror("unbalanced '[]' brackets in make_pair_table_pk");
   }
   free(stack);
   free(stack2);
   return(table);
}


PUBLIC short *
vrna_pt_snoop_get(const char *structure){

    /* returns array representation of structure.
       table[i] is 0 if unpaired or j if (i.j) pair.  */
   short i,j,hx;
   short length;
   short *stack;
   short *table;

   length = (short) strlen(structure);
   stack = (short *) space(sizeof(short)*(length+1));
   table = (short *) space(sizeof(short)*(length+2));
   table[0] = length;

   for (hx=0, i=1; i<=length; i++) {
     switch (structure[i-1]) {
     case '<':
       stack[hx++]=i;
       break;
     case '>':
       j = stack[--hx];
       if (hx<0) {
         fprintf(stderr, "%s\n", structure);
         nrerror("unbalanced brackets in make_pair_table");
       }
       table[i]=j;
       table[j]=i;
       break;
     default:   /* unpaired base, usually '.' */
       table[i]= table[i];
       break;
     }
   }
   if (hx!=0) {
     fprintf(stderr, "%s\n", structure);
     nrerror("unbalanced brackets in make_pair_table");
   }
   free(stack);
   return table ;
}



PUBLIC short *
vrna_pt_ali_get(const char *structure){

    /* returns array representation of structure.
       table[i] is 0 if unpaired or j if (i.j) pair.  */
   short i,j,hx;
   short length;
   short *stack;
   short *table;

   length = (short) strlen(structure);
   stack = (short *) space(sizeof(short)*(length+1));
   table = (short *) space(sizeof(short)*(length+2));
   table[0] = length;

   for (hx=0, i=1; i<=length; i++) {
      switch (structure[i-1]) {
       case '(':
         stack[hx++]=i;
         break;
       case ')':
         j = stack[--hx];
         if (hx<0) {
            fprintf(stderr, "%s\n", structure);
            nrerror("unbalanced brackets in make_pair_table");
         }
         table[i]=j;
         table[j]=i;
         break;
       default:   /* unpaired base, usually '.' */
         table[i]= 0;
         break;
      }
   }
   for (hx=0, i=1; i<=length; i++) {
      switch (structure[i-1]) {
       case '<':
         stack[hx++]=i;
         break;
       case '>':
         j = stack[--hx];
         if (hx<0) {
            fprintf(stderr, "%s\n", structure);
            nrerror("unbalanced brackets in make_pair_table");
         }
         table[i]=j;
         table[j]=i;
         break;
       default:   /* unpaired base, usually '.' */
         table[i]= table[i];
         break;
      }
   }
   for (hx=0, i=1; i<=length; i++) {
     switch (structure[i-1]) {
     case '[':
       stack[hx++]=i;
       break;
     case ']':
       j = stack[--hx];
       if (hx<0) {
         fprintf(stderr, "%s\n", structure);
         nrerror("unbalanced brackets in make_pair_table");
       }
       table[i]=j;
       table[j]=i;
       break;
     default:   /* unpaired base, usually '.' */
       break;
     }
   }
   if (hx!=0) {
      fprintf(stderr, "%s\n", structure);
      nrerror("unbalanced brackets in make_pair_table");
   }
   free(stack);
   return(table);
}

PUBLIC short *
vrna_pt_copy(const short *pt){
  short *table = (short *)space(sizeof(short) * (pt[0]+2));
  memcpy(table, pt, sizeof(short)*(pt[0]+2));
  return table;
}


PUBLIC int *
vrna_get_loop_index(const short *pt){

  /* number each position by which loop it belongs to (positions start
     at 1) */
  int i,hx,l,nl;
  int length;
  int *stack = NULL;
  int *loop = NULL;

  length = pt[0];
  stack  = (int *) space(sizeof(int)*(length+1));
  loop   = (int *) space(sizeof(int)*(length+2));
  hx=l=nl=0;

  for (i=1; i<=length; i++) {
    if ((pt[i] != 0) && (i < pt[i])) { /* ( */
      nl++; l=nl;
      stack[hx++]=i;
    }
    loop[i]=l;

    if ((pt[i] != 0) && (i > pt[i])) { /* ) */
      --hx;
      if (hx>0)
        l = loop[stack[hx-1]];  /* index of enclosing loop   */
      else l=0;                 /* external loop has index 0 */
      if (hx<0) {
        nrerror("unbalanced brackets in make_pair_table");
      }
    }
  }
  loop[0] = nl;
  free(stack);
  return (loop);
}

PUBLIC char *
vrna_pt_to_db(short *pt){

  int i;
  char *dotbracket = NULL;
  if(pt){
    dotbracket = (char *)space((pt[0]+1)*sizeof(char));

    for(i=1; i<=pt[0]; i++)
      dotbracket[i-1] = '.';
      if(pt[i] > i){
        dotbracket[i-1] = '(';
        dotbracket[pt[i]-1] = ')';
      }

    dotbracket[i-1] = '\0';
  }
  return dotbracket;
}

/*---------------------------------------------------------------------------*/

PUBLIC int
vrna_bp_distance(const char *str1, const char *str2){

  /* dist = {number of base pairs in one structure but not in the other} */
  /* same as edit distance with pair_open pair_close as move set */
   int dist;
   short i,l;
   short *t1, *t2;

   dist = 0;
   t1 = vrna_pt_get(str1);
   t2 = vrna_pt_get(str2);

   l = (t1[0]<t2[0])?t1[0]:t2[0];    /* minimum of the two lengths */

   for (i=1; i<=l; i++)
     if (t1[i]!=t2[i]) {
       if (t1[i]>i) dist++;
       if (t2[i]>i) dist++;
     }
   free(t1); free(t2);
   return dist;
}

/* get a matrix containing the number of basepairs of a reference structure for each interval [i,j] with i<j
*  access it via iindx!!!
*/
PUBLIC unsigned int *
vrna_refBPcnt_matrix( const short *reference_pt,
                      unsigned int turn){

  unsigned int i,j,k,ij,length;
  int *iindx;
  unsigned int *array;
  unsigned int size;
  length = (unsigned int)reference_pt[0];
  size  = ((length+1)*(length+2))/2;
  iindx = get_iindx(length);
  array = (unsigned int *) space(sizeof(unsigned int)*size);    /* matrix containing number of basepairs of reference structure1 in interval [i,j] */;
  for (k=0; k<=turn; k++)
    for (i=1; i<=length-k; i++) {
      j=i+k;
      ij = iindx[i]-j;
      array[ij] = 0;
    }

  for (i = length-turn-1; i >= 1; i--)
    for (j = i+turn+1; j <= length; j++){
      int bps;
      ij = iindx[i]-j;
      bps = array[ij+1];
      if((i<=(unsigned int)reference_pt[j]) && ((unsigned int)reference_pt[j] < j))
        bps++;
      array[ij] = bps;
    }
  free(iindx);
  return array;
}


PUBLIC unsigned int *
vrna_refBPdist_matrix(const short *pt1,
                      const short *pt2,
                      unsigned int turn){

  unsigned int *array;
  unsigned int n, size, i, j, ij, d;
  n = (unsigned int)pt1[0];
  size = ((n+1)*(n+2))/2;
  array = (unsigned int *)space(sizeof(unsigned int) * size);
  int *iindx = get_iindx(n);
  for(i = n - turn - 1; i>=1; i--){
    d = 0;
    for(j = i+turn+1; j <= n; j++){
      ij = iindx[i]-j;
      d = array[ij+1];
      if(pt1[j] != pt2[j]){
        if(i <= (unsigned int)pt1[j] && (unsigned int)pt1[j] < j){
          /* we got an additional base pair in reference structure 1 */
          d++;
        }
        if(i <= (unsigned int)pt2[j] && (unsigned int)pt2[j] < j){
          /* we got another base pair in reference structure 2 */
          d++;
        }
      }
      array[ij] = d;

    }
  }
  free(iindx);
  return array;
}

PUBLIC char
bppm_symbol(const float *x){

/*  if( ((x[1]-x[2])*(x[1]-x[2]))<0.1&&x[0]<=0.677) return '|'; */
  if( x[0] > 0.667 )  return '.';
  if( x[1] > 0.667 )  return '(';
  if( x[2] > 0.667 )  return ')';
  if( (x[1]+x[2]) > x[0] ) {
    if( (x[1]/(x[1]+x[2])) > 0.667) return '{';
    if( (x[2]/(x[1]+x[2])) > 0.667) return '}';
    else return '|';
  }
  if( x[0] > (x[1]+x[2]) ) return ',';
  return ':';
}

PUBLIC void
bppm_to_structure(char *structure,
                  FLT_OR_DBL *p,
                  unsigned int length){

  int    i, j;
  int   *index = get_iindx(length);
  float  P[3];   /* P[][0] unpaired, P[][1] upstream p, P[][2] downstream p */

  for( j=1; j<=length; j++ ) {
    P[0] = 1.0;
    P[1] = P[2] = 0.0;
    for( i=1; i<j; i++) {
      P[2] += p[index[i]-j];    /* j is paired downstream */
      P[0] -= p[index[i]-j];    /* j is unpaired */
    }
    for( i=j+1; i<=length; i++ ) {
      P[1] += p[index[j]-i];    /* j is paired upstream */
      P[0] -= p[index[j]-i];    /* j is unpaired */
    }
    structure[j-1] = bppm_symbol(P);
  }
  structure[length] = '\0';
  free(index);
}

PUBLIC plist *
vrna_get_plist_from_pr( vrna_fold_compound *vc,
                        double cut_off){

  int i, j, k, n, count, *index, length;
  FLT_OR_DBL  *probs, *G, *scale;
  pf_matricesT  *matrices;
  short         *S;
  pf_paramT     *pf_params;
  plist         *pl;

  if(!vc){
    nrerror("vrna_get_plist_from_pr: run vrna_pf_fold first!");
  } else if( !vc->exp_matrices->probs){
    nrerror("vrna_get_plist_from_pr: probs==NULL!");
  }

  return wrap_get_plist(vc->exp_matrices,
                        vc->length,
                        vc->iindx,
                        vc->sequence_encoding2,
                        vc->exp_params,
                        cut_off);
}

PRIVATE plist *
wrap_get_plist( pf_matricesT *matrices,
                int length,
                int *index,
                short *S,
                pf_paramT *pf_params,
                double cut_off){

  int i, j, k, n, count, gquad;
  FLT_OR_DBL  *probs, *G, *scale;
  plist         *pl;

  probs     = matrices->probs;
  G         = matrices->G;
  scale     = matrices->scale;
  gquad     = pf_params->model_details.gquad;

  count = 0;
  n     = 2;

  /* first guess of the size needed for pl */
  pl = (plist *)space(n*length*sizeof(plist));

  for (i=1; i<length; i++) {
    for (j=i+1; j<=length; j++) {
      /* skip all entries below the cutoff */
      if (probs[index[i]-j] < cut_off) continue;

      /* do we need to allocate more memory? */
      if (count == n * length - 1){
        n *= 2;
        pl = (plist *)xrealloc(pl, n * length * sizeof(plist));
      }

      /* check for presence of gquadruplex */
      if((S[i] == 3) && (S[j] == 3) && gquad){
          /* add probability of a gquadruplex at position (i,j)
             for dot_plot
          */
          (pl)[count].i      = i;
          (pl)[count].j      = j;
          (pl)[count].p      = probs[index[i] - j];
          (pl)[count++].type = 1;
          /* now add the probabilies of it's actual pairing patterns */
          plist *inner, *ptr;
          inner = get_plist_gquad_from_pr(S, i, j, G, probs, scale, pf_params);
          for(ptr=inner; ptr->i != 0; ptr++){
              if (count == n * length - 1){
                n *= 2;
                pl = (plist *)xrealloc(pl, n * length * sizeof(plist));
              }
              /* check if we've already seen this pair */
              for(k = 0; k < count; k++)
                if(((pl)[k].i == ptr->i) && ((pl)[k].j == ptr->j))
                  break;
              (pl)[k].i      = ptr->i;
              (pl)[k].j      = ptr->j;
              (pl)[k].type = 0;
              if(k == count){
                (pl)[k].p  = ptr->p;
                count++;
              } else
                (pl)[k].p  += ptr->p;
          }
      } else {
          (pl)[count].i      = i;
          (pl)[count].j      = j;
          (pl)[count].p      = probs[index[i] - j];
          (pl)[count++].type = 0;
      }
    }
  }
  /* mark the end of pl */
  (pl)[count].i    = 0;
  (pl)[count].j    = 0;
  (pl)[count].type = 0;
  (pl)[count++].p  = 0.;
  /* shrink memory to actual size needed */
  pl = (plist *)xrealloc(pl, count * sizeof(plist));

  return pl;
}

PUBLIC void
vrna_letter_structure(char *structure,
                      bondT *bp,
                      int length){

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

PUBLIC void
vrna_parenthesis_structure( char *structure,
                            bondT *bp,
                            int length){

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

PUBLIC void
vrna_parenthesis_zuker( char *structure,
                        bondT *bp,
                        int length){

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

PUBLIC plist *
vrna_get_plist_from_db( const char *struc,
                        float pr){

  /* convert bracket string to plist */
  short *pt;
  int i, k = 0, size, n;
  plist *gpl, *ptr, *pl;

  size  = strlen(struc);
  n     = 2;

  pt  = vrna_pt_get(struc);
  pl = (plist *)space(n*size*sizeof(plist));
  for(i = 1; i < size; i++){
    if(pt[i]>i){
      (pl)[k].i      = i;
      (pl)[k].j      = pt[i];
      (pl)[k].p      = pr;
      (pl)[k++].type = 0;
    }
  }

  gpl = get_plist_gquad_from_db(struc, pr);
  for(ptr = gpl; ptr->i != 0; ptr++){
    if (k == n * size - 1){
      n *= 2;
      pl = (plist *)xrealloc(pl, n * size * sizeof(plist));
    }
    (pl)[k].i      = ptr->i;
    (pl)[k].j      = ptr->j;
    (pl)[k].p       = ptr->p;
    (pl)[k++].type = ptr->type;
  }
  free(gpl);

  (pl)[k].i      = 0;
  (pl)[k].j      = 0;
  (pl)[k].p      = 0.;
  (pl)[k++].type = 0.;
  free(pt);
  pl = (plist *)xrealloc(pl, k * sizeof(plist));

  return pl;
}

/*###########################################*/
/*# deprecated functions below              #*/
/*###########################################*/

PUBLIC char *
pack_structure(const char *struc){

  return vrna_pack_structure(struc);
}

PUBLIC char *
unpack_structure(const char *packed){

  return vrna_unpack_structure(packed);
}

PUBLIC void
parenthesis_structure(char *structure,
                      bondT *bp,
                      int length){

  return vrna_parenthesis_structure(structure, bp, length);
}

PUBLIC void
letter_structure( char *structure,
                  bondT *bp,
                  int length){

  vrna_letter_structure(structure, bp, length);
}

PUBLIC void
parenthesis_zuker(char *structure,
                  bondT *bp,
                  int length){

  return vrna_parenthesis_zuker(structure, bp, length);
}

PUBLIC void
assign_plist_from_pr( plist **pl,
                      FLT_OR_DBL *probs,
                      int length,
                      double cut_off){

  int i, j, n, count, *index;

  index = get_iindx(length);
  pf_matricesT  *matrices = (pf_matricesT *)space(sizeof(pf_matricesT));

  matrices->probs = probs;
  model_detailsT  md;
  set_model_details(&md);
  md.gquad = 0;
  pf_paramT *pf_params = vrna_get_boltzmann_factors(md);

  *pl = wrap_get_plist( matrices,
                        length,
                        index,
                        NULL,
                        pf_params,
                        cut_off);

  free(index);
  free(pf_params);
  free(matrices);
}

PUBLIC void
assign_plist_from_db( plist **pl,
                      const char *struc,
                      float pr){

  *pl = vrna_get_plist_from_db(struc, pr);
}

PUBLIC short *
make_pair_table(const char *structure){

  return vrna_pt_get(structure);
}

PUBLIC short *
copy_pair_table(const short *pt){

  return vrna_pt_copy(pt);
}

PUBLIC short *
make_pair_table_pk(const char *structure){

  return vrna_pt_pk_get(structure);
}

PUBLIC short *
make_pair_table_snoop(const char *structure){

  return vrna_pt_snoop_get(structure);
}

PUBLIC short *
alimake_pair_table(const char *structure){

  return vrna_pt_ali_get(structure);
}

PUBLIC int *
make_loop_index_pt(short *pt){

  return vrna_get_loop_index((const short*)pt);
}

PUBLIC int
bp_distance(const char *str1, const char *str2){

  return vrna_bp_distance(str1, str2);
}

PUBLIC unsigned int *
make_referenceBP_array( short *reference_pt,
                        unsigned int turn){

  return vrna_refBPcnt_matrix((const short *)reference_pt, turn);
}

PUBLIC unsigned int *
compute_BPdifferences(short *pt1,
                      short *pt2,
                      unsigned int turn){

  return vrna_refBPdist_matrix((const short *)pt1, (const short *)pt2, turn);
}

