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
#include "ViennaRNA/structure_utils.h"


PUBLIC char *pack_structure(const char *struc) {
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

PUBLIC char *unpack_structure(const char *packed) {
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

PUBLIC short *make_pair_table(const char *structure)
{
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

PUBLIC short *make_pair_table_pk(const char *structure){
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

PUBLIC short *make_pair_table_snoop(const char *structure)
{
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


PUBLIC short *alimake_pair_table(const char *structure)
{
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

PUBLIC short *copy_pair_table(const short *pt){
  short *table = (short *)space(sizeof(short) * (pt[0]+2));
  memcpy(table, pt, sizeof(short)*(pt[0]+2));
  return table;
}


PUBLIC int *make_loop_index_pt(short *pt){

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

PUBLIC int bp_distance(const char *str1, const char *str2)
{
  /* dist = {number of base pairs in one structure but not in the other} */
  /* same as edit distance with pair_open pair_close as move set */
   int dist;
   short i,l;
   short *t1, *t2;

   dist = 0;
   t1 = make_pair_table(str1);
   t2 = make_pair_table(str2);

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
PUBLIC unsigned int *make_referenceBP_array(short *reference_pt, unsigned int turn){
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

PUBLIC unsigned int *compute_BPdifferences(short *pt1, short *pt2, unsigned int turn){
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
