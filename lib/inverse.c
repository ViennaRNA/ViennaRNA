/*
		      search for sequences that
		  fold into a given target structure
				   
			    c Ivo Hofacker
			  Vienna RNA package
*/
/* Last changed Time-stamp: <97/11/04 19:15:02 ivo> */

#define TDIST 0     /* use tree distance */
#define PF    1     /* include support for partiton function */

#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#if PF
#include "part_func.h"
#endif
#include "fold.h"
#if TDIST
#include "dist_vars.h"
#include "treedist.h"
#include "RNAstruct.h"
#endif
#include "utils.h"
#include "fold_vars.h"
#include "pair_mat.h"
#ifdef dmalloc
#include "/usr/local/include/dmalloc.h"
#define space(X) calloc(1,(X))
#endif

static char rcsid[] = "$Id: inverse.c,v 1.4 1997/11/04 18:33:10 ivo Rel $";
#define PUBLIC
#define PRIVATE static
PRIVATE float  adaptive_walk(char *start, char *target);
PRIVATE void   shuffle(short int *list, short int len);
PRIVATE void   make_start(char* start, char *structure);
PRIVATE void   make_pair_table(char *structure, short *table);
PRIVATE void   make_pairset(void);
PRIVATE float  mfe_cost(char *, char*, char *);
PRIVATE float  pf_cost(char *, char *, char *);
PRIVATE char  *aux_struct( char* structure );
PRIVATE int    bp_distance(char *str1, char *str2);

PUBLIC  char   symbolset[MAXALPHA+1] = "AUGC";
PUBLIC  int    give_up = 0;
PUBLIC  float  final_cost = 0; /* when to stop inverse_pf_fold */
PUBLIC  int    inv_verbose=1;

PRIVATE char   pairset[2*MAXALPHA+1];
PRIVATE int    base, npairs;
PRIVATE int    nc2;

/*-------------------------------------------------------------------------*/
PRIVATE int fold_type;
#if TDIST
PRIVATE Tree *T0;
#endif
PRIVATE float cost2;

PRIVATE float adaptive_walk(char *start, char *target)
{
#ifdef DUMMY
   printf("%s\n%s %c\n", start, target, backtrack_type ); 
   return 0.;
#endif
   int i,j,p,tt,w1,w2, n_pos, len, flag;
   long  walk_len;
   char *string, *string2, *cstring, *structure, *struct2;
   short *mut_pos_list, mut_sym_list[MAXALPHA+1], mut_pair_list[2*MAXALPHA+1];
   short *w1_list, *w2_list, mut_position, symbol, bp;
   short *target_table, *test_table;
   char cont;
   float cost, current_cost, ccost2;
   float (*cost_function)(char *, char *, char *);

   len = strlen(start);
   if (strlen(target)!=len) {
      fprintf(stderr, "%s\n%s\n", start, target);
      nrerror("adaptive_walk: start and target have unequal length");
   }
   string    = (char *) space(sizeof(char)*(len+1));
   cstring   = (char *) space(sizeof(char)*(len+1));
   string2   = (char *) space(sizeof(char)*(len+1));
   structure = (char *) space(sizeof(char)*(len+1));
   struct2   = (char *) space(sizeof(char)*(len+1)); 
   mut_pos_list = (short *) space(sizeof(short)*len);
   w1_list = (short *) space(sizeof(short)*len);
   w2_list = (short *) space(sizeof(short)*len);
   target_table = (short *) space(sizeof(short)*len);
   test_table = (short *) space(sizeof(short)*len);
   
   make_pair_table(target, target_table);
   
   for (i=0; i<base; i++) mut_sym_list[i] = i;
   for (i=0; i<npairs; i++) mut_pair_list[i] = i;
   
   for (i=0; i<len; i++) 
      string[i] = (islower(start[i]))?toupper(start[i]):start[i];
   walk_len = 0;

   if (fold_type==0) cost_function = mfe_cost;
   else cost_function = pf_cost;

   cost = cost_function(string, structure, target);

   if (fold_type==0) ccost2=cost2;
   else { ccost2 = -1.; cost2=0; }
   
   strcpy(cstring, string);
   current_cost = cost;
      
   if (cost>0) do {
      cont=0;

      if (fold_type==0) { /* min free energy fold */
	 make_pair_table(structure, test_table);
	 for (j=w1=w2=flag=0; j<len; j++)
	    if ((tt=target_table[j])!=test_table[j]) {
	       if ((tt<j)&&(isupper(start[j]))) w1_list[w1++] = j;   /* incorrectly paired */
	       if ((flag==0)&&(j>0))
		  if ((target_table[j-1]<j-1)&&isupper(start[j-1]))
			w2_list[w2++] = j-1;                  /* adjacent to incorrect position */
	       if (w2>1) if (w2_list[w2-2]==w2_list[w2-1]) w2--;
		     
	       flag = 1;
	    } else {
	       if (flag==1) if ((tt<j)&&isupper(start[j]))
		  w2_list[w2++] = j;                          /* adjacent to incorrect position */
	       flag = 0;
	    }
	 shuffle(w1_list, (short) w1);
	 shuffle(w2_list, (short) w2);
	 for (j=n_pos=0; j<w1; j++) mut_pos_list[n_pos++] = w1_list[j];
	 for (j=0; j<w2; j++) mut_pos_list[n_pos++] = w2_list[j];
      } else { /* partition_function */
	 for (j=n_pos=0; j<len; j++) if (isupper(start[j]))
	    if (target_table[j]<=j) mut_pos_list[n_pos++] = j;
	 shuffle(mut_pos_list, (short) n_pos);
      }

      string2[0]='\0';
      for (mut_position=0; mut_position<n_pos; mut_position++){
	 
	 strcpy(string, cstring);
	 shuffle(mut_sym_list, (short) base);
	 shuffle(mut_pair_list, (short) npairs);
	 
	 i = mut_pos_list[mut_position];

	 if (target_table[i]<0) /* unpaired base */
	    for (symbol=0;symbol<base;symbol++) {
	       
	       if(cstring[i]==
		  symbolset[mut_sym_list[symbol]]) continue;
	       
	       string[i] = symbolset[mut_sym_list[symbol]];
	       
	       cost = cost_function(string, structure, target);
	       
	       if ( cost < current_cost ) break;
	       if (( cost == current_cost)&&(cost2<ccost2)){
		  strcpy(string2, string);
		  strcpy(struct2, structure);
		  ccost2 = cost2;
	       }
	    }
	 else  /* paired base */
	    for  (bp=0; bp<npairs; bp++) {
	       j = target_table[i];
	       p = mut_pair_list[bp]*2;
	       if ((cstring[i] == pairset[p]) &&
		   (cstring[j] == pairset[p+1]))
		  continue;
	       string[i] = pairset[p];
	       string[j] = pairset[p+1];

	       cost = cost_function(string, structure, target);
	       
	       if ( cost < current_cost ) break;
	       if (( cost == current_cost)&&(cost2<ccost2)){
		  strcpy(string2, string);
		  strcpy(struct2, structure);
		  ccost2 = cost2;
	       }
	    }
	 
	 if ( cost < current_cost ) {
	    strcpy(cstring, string);
	    current_cost = cost;
	    ccost2 = cost2;
	    walk_len++;
	    if (cost>0) cont=1;
	    break;
	 }
      }
      if ((current_cost>0)&&(cont==0)&&(string2[0])) {
	 /* no mutation that decreased cost was found, 
	    but the the sequence in string2 decreases cost2 while keeping
	    cost constant */
	 strcpy(cstring, string2);
	 strcpy(structure, struct2);
	 nc2++; cont=1;
      }
   } while (cont);

   for (i=0; i<len; i++) if (isupper(start[i])) start[i]=cstring[i];

#if TDIST
   if (fold_type==0) { free_tree(T0); T0=NULL; }
#endif
   free(test_table);
   free(target_table);
   free(mut_pos_list);
   free(w2_list);
   free(w1_list);
   free(struct2);
   free(structure);
   free(string2);
   free(cstring);
   free(string);

   return current_cost;
}

/*-------------------------------------------------------------------------*/


PRIVATE void shuffle(short int *list, short int len)
{
   short i, temp, rn;
   
   for (i=0;i<len;i++)
       {
	  temp = list[i];
	  rn = i + (int) (urn()*(len-i));   /* [i..len-1] */
	  list[i] = list[rn];
	  list[rn] = temp;
       }
}

/*-------------------------------------------------------------------------*/

PRIVATE void make_pair_table(char *structure, short *table)
{
   int i,j,hx;
   short *stack;
   
   hx=0;
   stack = (short *) space(sizeof(short)*(strlen(structure)+1));
	     
   for (i=0; i<strlen(structure); i++) {
      switch (structure[i]) {
       case '.':
	 table[i]= -1;
	 break;
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
      }
   }
   if (hx!=0) {
      fprintf(stderr, "%s\n", structure);
      nrerror("unbalanced brackets in make_pair_table");
   }
   free(stack);
}
 
/*-------------------------------------------------------------------------*/

#define WALK(i,j) \
    strncpy(wstruct, structure+i, j-i+1); \
    wstruct[j-i+1]='\0'; \
    strncpy(wstring, string+i, j-i+1); \
    wstring[j-i+1]='\0'; \
    dist=adaptive_walk(wstring, wstruct); \
    strncpy(string+i, wstring, j-i+1); \
    if ((dist>0)&&(give_up)) goto adios
 
PUBLIC float inverse_fold(char *start, char *structure)
{
   int i, j, jj, len, o;
   short *pt;
   char *string, *wstring, *wstruct, *aux;
   float dist=0;
   
   nc2 = j = o = fold_type = 0;
   
   len = strlen(structure);
   if (strlen(start)!=len) {
      fprintf(stderr, "%s\n%s\n", start, structure);
      nrerror("inverse_fold: start and structure have unequal length");
   }
   string = (char *) space(len+1);
   wstring = (char *) space(len+1);
   wstruct = (char *) space(len+1);
   pt = (short *) space(sizeof(short)*(len+1));
   pt[len] = len+1;
   
   aux = aux_struct(structure);
   strcpy(string, start);
   make_pairset();
   make_start(string, structure);
   
   make_pair_table(structure, pt);
   
   while (j<len) {
      while ((j<len)&&(structure[j]!=')')) {
	 if (aux[j]=='[') o++;
	 if (aux[j]==']') o--;
	 j++;
      }
      i=j;
      while (structure[--i]!='(');  /* doesn't work for open structure */
      if (aux[i]!='[') { i--; j++;}
      while (pt[j]==i) {
	 backtrack_type='C';
	 WALK(i,j);
	 if (aux[i]!='[') {
	    while (aux[--i]!='[');
	    while (aux[++j]!=']');
	    WALK(i,j);
	 }
	 o--;
	 jj = j; i--;
	 while (aux[++j]=='.');
	 while ((i>=0)&&(aux[i]=='.')) i--;
	 if (pt[j]!=i) {
	    backtrack_type = (o==0)? 'F' : 'M';
	    if (j-jj>8) { WALK((i+1),(jj)); }
	    WALK((i+1), (j-1));
	    while ((i>=0) &&(aux[i]==']')) {
	       i=pt[i]-1;
	       while ((i>=0)&&(aux[i]=='.')) i--;
	       WALK((i+1), (j-1));
	    }
	 }
      }
   }
 adios:
   backtrack_type='F';
   if ((dist>0)&&(inv_verbose)) printf("%s\n%s\n", wstring, wstruct);
   /*if ((dist==0)||(give_up==0))*/ strcpy(start, string);
   free(wstring); free(wstruct);
   free(string); free(aux);
   free(pt);
/*   if (dist>0) printf("%3d \n", nc2); */
   return dist;
}

/*-------------------------------------------------------------------------*/

PUBLIC float inverse_pf_fold(char *start, char *target)
{
   float dist;

   update_fold_params();    /* make sure there is a valid pair matrix */
   make_pairset();
   make_start(start, target);
   fold_type=1;
   do_backtrack = 0;
   dist = adaptive_walk(start, target);
   return (dist+final_cost);
}

/*-------------------------------------------------------------------------*/

PRIVATE void make_start(char* start, char *structure)
{
   int i,j,k,l,r,length;
   short *table, *S, sym[MAXALPHA], ss;

   length=strlen(start);
   table = (short *) space(sizeof(short)*length);
   S = (short *) space(sizeof(short)*length);
   
   make_pair_table(structure, table);
   for (i=0; i<strlen(start); i++) S[i] = ENCODE(toupper(start[i]));
   for (i=0; i<strlen(symbolset); i++) sym[i] = i;

   for (k=0; k<length; k++) {
      if (table[k]<k) continue;
      if (urn()<0.5) {
	 i = k; j = table[k];
      } else {
	 i = table[k]; j = k;
      }

      if (!pair[S[i]][S[j]]) {
	 shuffle(sym, (short) base);
	 for (l=0; l<base; l++) {
	    ss = ENCODE(symbolset[sym[l]]);
	    if (pair[S[i]][ss]) break;
	 }
	 if (l==base) { /* nothing pairs start[i] */
	    r = 2*int_urn(0, npairs-1);
	    start[i] = pairset[r];
	    start[j] = pairset[r+1];
	 } else start[j] = symbolset[sym[l]];
      }
   }
   free(table);
   free(S);
}

/*---------------------------------------------------------------------------*/

PRIVATE void make_pairset(void)
{
   int i,j;
   short sym[MAXALPHA];
   
   make_pair_matrix();
   base = strlen(symbolset);

   for (i=0; i< base; i++) sym[i] = ENCODE(symbolset[i]);
   
   for (i=npairs=0; i< base; i++)
      for (j=0; j<base; j++) 
	 if (pair[sym[i]][sym[j]]) {
	    pairset[npairs++] = symbolset[i];
	    pairset[npairs++] = symbolset[j];
	 }
   npairs /= 2;
   if (npairs==0) nrerror("No pairs in this alphabet!");
}
/*---------------------------------------------------------------------------*/

PRIVATE float mfe_cost(char *string, char *structure, char *target)
{
#if TDIST
   Tree *T1;
   char *xstruc;
#endif
   float energy, distance;

   if (strlen(string)!=strlen(target)) {
      fprintf(stderr, "%s\n%s\n", string, target);
      nrerror("unequal length in mfe_cost");
   }
   energy = fold(string, structure);
#if TDIST
   if (T0 == NULL) {
      xstruc = expand_Full(target);
      T0=make_tree(xstruc);
      free(xstruc);
   }

   xstruc = expand_Full(structure);
   T1=make_tree(xstruc);
   distance = tree_edit_distance(T0,T1);
   free(xstruc);
   free_tree(T1);
#else
   distance = bp_distance(target, structure);
#endif
   cost2 = energy_of_struct(string, target) - energy;
   return (float) distance;
}
/*---------------------------------------------------------------------------*/

PRIVATE float pf_cost(char *string, char *structure, char *target)
{
#if PF
   float  f, e;
   
   f = pf_fold(string, structure);
   e = energy_of_struct(string, target);
   return (float) (e-f-final_cost);
#else
   nrerror("this version not linked with pf_fold");
   return 0;
#endif
}

/*---------------------------------------------------------------------------*/

PRIVATE char *aux_struct( char* structure )
{  
   short       *match_paren;
   int          i, o, p;
   char        *string;
   
   string = (char *) space(sizeof(char)*(strlen(structure)+1));
   match_paren = (short *) space(sizeof(short)*(strlen(structure)/2+1));
   strcpy(string, structure);
   
   i = o = 0;
   while (string[i]) {
      switch (string[i]) {
       case '.': break;
       case '(':
         match_paren[++o]=i;
         break;
       case ')':
         p=i;
         while ((string[p+1]==')')&&(match_paren[o-1]==match_paren[o]-1)) {
            p++; o--;
         }
         string[p]=']';
         i=p;
         string[match_paren[o]]='[';
         o--;
         break;
       default:
         nrerror("Junk in structure at aux_structure\n");
      }
      i++;
   }
   free(match_paren);
   return(string);
}

PRIVATE int bp_distance(char *str1, char *str2)
{
   int dist,i;
   short *t1, *t2;

   dist = 0;
   t1 = (short *) space(sizeof(short)*strlen(str1));
   t2 = (short *) space(sizeof(short)*strlen(str2));
   make_pair_table(str1, t1);
   make_pair_table(str2, t2);

   for (i=0; i<strlen(str1); i++)
      dist += (t1[i]!=t2[i]);
   free(t1); free(t2);
   return dist;
}
