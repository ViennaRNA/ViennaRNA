/*
		      search for sequences that
		  fold into a given target structure

			    c Ivo Hofacker
			  Vienna RNA package
*/

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#define TDIST 0     /* use tree distance */
#define PF    1     /* include support for partiton function */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <float.h>
#if PF
#include "ViennaRNA/part_func.h"
#endif
#include "ViennaRNA/fold.h"
#if TDIST
#include "ViennaRNA/dist_vars.h"
#include "ViennaRNA/treedist.h"
#include "ViennaRNA/RNAstruct.h"
#endif
#include "ViennaRNA/utils.h"
#include "ViennaRNA/fold_vars.h"
#include "ViennaRNA/pair_mat.h"

PRIVATE double  adaptive_walk(char *start, const char *target);
PRIVATE void   shuffle(int *list, int len);
PRIVATE void   make_start(char* start, const char *structure);
PRIVATE void   make_ptable(const char *structure, int *table);
PRIVATE void   make_pairset(void);
PRIVATE double  mfe_cost(const char *, char*, const char *);
PRIVATE double  pf_cost(const char *, char *, const char *);
PRIVATE char  *aux_struct(const char* structure );

/* for backward compatibility, make sure symbolset can hold 20 characters */
PRIVATE char   default_alpha[21] = "AUGC";
PUBLIC  char   *symbolset = default_alpha;
PUBLIC  int    give_up = 0;
PUBLIC  float  final_cost = 0; /* when to stop inverse_pf_fold */
PUBLIC  int    inv_verbose=0;  /* print out substructure on which inverse_fold() fails */

PRIVATE char   pairset[2*MAXALPHA+1];
PRIVATE int    base, npairs;
PRIVATE int    nc2;

/*-------------------------------------------------------------------------*/
PRIVATE int fold_type;
#if TDIST
PRIVATE Tree *T0;
#endif
PRIVATE double cost2;

PRIVATE double adaptive_walk(char *start, const char *target)
{
#ifdef DUMMY
   printf("%s\n%s %c\n", start, target, backtrack_type );
   return 0.;
#endif
   int i,j,p,tt,w1,w2, n_pos, len, flag;
   long  walk_len;
   char *string, *string2, *cstring, *structure, *struct2;
   int *mut_pos_list, mut_sym_list[MAXALPHA+1], mut_pair_list[2*MAXALPHA+1];
   int *w1_list, *w2_list, mut_position, symbol, bp;
   int *target_table, *test_table;
   char cont;
   double cost, current_cost, ccost2;
   double (*cost_function)(const char *, char *, const char *);

   len = strlen(start);
   if (strlen(target)!=len) {
      vrna_message_error("%s\n%s\nadaptive_walk: start and target have unequal length", start, target);
   }
   string    = (char *) vrna_alloc(sizeof(char)*(len+1));
   cstring   = (char *) vrna_alloc(sizeof(char)*(len+1));
   string2   = (char *) vrna_alloc(sizeof(char)*(len+1));
   structure = (char *) vrna_alloc(sizeof(char)*(len+1));
   struct2   = (char *) vrna_alloc(sizeof(char)*(len+1));
   mut_pos_list = (int *) vrna_alloc(sizeof(int)*len);
   w1_list = (int *) vrna_alloc(sizeof(int)*len);
   w2_list = (int *) vrna_alloc(sizeof(int)*len);
   target_table = (int *) vrna_alloc(sizeof(int)*len);
   test_table = (int *) vrna_alloc(sizeof(int)*len);

   make_ptable(target, target_table);

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
	 make_ptable(structure, test_table);
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
	 shuffle(w1_list, w1);
	 shuffle(w2_list, w2);
	 for (j=n_pos=0; j<w1; j++) mut_pos_list[n_pos++] = w1_list[j];
	 for (j=0; j<w2; j++) mut_pos_list[n_pos++] = w2_list[j];
      } else { /* partition_function */
	 for (j=n_pos=0; j<len; j++) if (isupper(start[j]))
	    if (target_table[j]<=j) mut_pos_list[n_pos++] = j;
	 shuffle(mut_pos_list, n_pos);
      }

      string2[0]='\0';
      for (mut_position=0; mut_position<n_pos; mut_position++){

	 strcpy(string, cstring);
	 shuffle(mut_sym_list,  base);
	 shuffle(mut_pair_list, npairs);

	 i = mut_pos_list[mut_position];

	 if (target_table[i]<0) /* unpaired base */
	    for (symbol=0;symbol<base;symbol++) {

	       if(cstring[i]==
		  symbolset[mut_sym_list[symbol]]) continue;

	       string[i] = symbolset[mut_sym_list[symbol]];

	       cost = cost_function(string, structure, target);

	       if ( cost + DBL_EPSILON < current_cost  ) break;
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

/* shuffle produces a ronaom list by doing len exchanges */
PRIVATE void shuffle(int *list, int len)
{
   int i, rn;

   for (i=0;i<len;i++) {
     int temp;
     rn = i + (int) (vrna_urn()*(len-i));   /* [i..len-1] */
     /* swap element i and rn */
     temp = list[i];
     list[i] = list[rn];
     list[rn] = temp;
   }
}

/*-------------------------------------------------------------------------*/

PRIVATE void make_ptable(const char *structure, int *table)
{
   int i,j,hx;
   int *stack;

   hx=0;
   stack = (int *) vrna_alloc(sizeof(int)*(strlen(structure)+1));

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
	    vrna_message_error("%s\nunbalanced brackets in make_ptable", structure);
	 }
	 table[i]=j;
	 table[j]=i;
	 break;
      }
   }
   if (hx!=0) {
      vrna_message_error("%s\nunbalanced brackets in make_ptable", structure);
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
   int *pt;
   char *string, *wstring, *wstruct, *aux;
   double dist=0;

   nc2 = j = o = fold_type = 0;

   len = strlen(structure);
   if (strlen(start)!=len) {
      vrna_message_error("%s\n%s\ninverse_fold: start and structure have unequal length", start, structure);
   }
   string = (char *) vrna_alloc(len+1);
   wstring = (char *) vrna_alloc(len+1);
   wstruct = (char *) vrna_alloc(len+1);
   pt = (int *) vrna_alloc(sizeof(int)*(len+2));
   pt[len] = len+1;

   aux = aux_struct(structure);
   strcpy(string, start);
   make_pairset();
   make_start(string, structure);

   make_ptable(structure, pt);

   while (j<len) {
      while ((j<len)&&(structure[j]!=')')) {
	 if (aux[j]=='[') o++;
	 if (aux[j]==']') o--;
	 j++;
      }
      i=j;
      while ((i>0) && structure[--i]!='(');
      if (structure[i]=='.') { /* no pair found -> open chain */
	WALK(0,len-1);
      }

      if (aux[i]!='[') { i--; j++;}
      while (pt[j]==i) {
	 backtrack_type='C';
	 if (aux[i]!='[') {
	    while (aux[--i]!='[');
	    while (aux[++j]!=']');
	    /* WALK(i,j); */
	 }
	 WALK(i,j);
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
   double dist;
   int dang;

   dang=dangles;
   if (dangles!=0) dangles=2;

   update_fold_params();    /* make sure there is a valid pair matrix */
   make_pairset();
   make_start(start, target);
   fold_type=1;
   do_backtrack = 0;
   dist = adaptive_walk(start, target);
   dangles=dang;
   return (dist+final_cost);
}

/*-------------------------------------------------------------------------*/

PRIVATE void make_start(char* start, const char *structure)
{
   int i,j,k,l,r,length;
   int *table, *S, sym[MAXALPHA], ss;

   length=strlen(start);
   table = (int *) vrna_alloc(sizeof(int)*length);
   S = (int *) vrna_alloc(sizeof(int)*length);

   make_ptable(structure, table);
   for (i=0; i<strlen(start); i++) S[i] = encode_char(toupper(start[i]));
   for (i=0; i<strlen(symbolset); i++) sym[i] = i;

   for (k=0; k<length; k++) {
      if (table[k]<k) continue;
      if (((vrna_urn()<0.5) && isupper(start[k])) ||
	  islower(start[table[k]])) {
	i = table[k]; j = k;
      } else {
	i = k; j = table[k];
      }

      if (!pair[S[i]][S[j]]) {   /* make a valid pair by mutating j */
	shuffle(sym, (int) base);
	for (l=0; l<base; l++) {
	  ss = encode_char(symbolset[sym[l]]);
	  if (pair[S[i]][ss]) break;
	}
	if (l==base) { /* nothing pairs start[i] */
	  r = 2*vrna_int_urn(0, npairs-1);
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
   int sym[MAXALPHA];

   make_pair_matrix();
   base = strlen(symbolset);

   for (i=0; i< base; i++) sym[i] = encode_char(symbolset[i]);

   for (i=npairs=0; i< base; i++)
      for (j=0; j<base; j++)
	 if (pair[sym[i]][sym[j]]) {
	    pairset[npairs++] = symbolset[i];
	    pairset[npairs++] = symbolset[j];
	 }
   npairs /= 2;
   if (npairs==0) vrna_message_error("No pairs in this alphabet!");
}
/*---------------------------------------------------------------------------*/

PRIVATE double mfe_cost(const char *string, char *structure, const char *target)
{
#if TDIST
   Tree *T1;
   char *xstruc;
#endif
   double energy, distance;

   if (strlen(string)!=strlen(target)) {
      vrna_message_error("%s\n%s\nunequal length in mfe_cost", string, target);
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
   distance = (double) vrna_bp_distance(target, structure);
#endif
   cost2 = energy_of_structure(string, target, 0) - energy;
   return (double) distance;
}
/*---------------------------------------------------------------------------*/

PRIVATE double pf_cost(const char *string, char *structure, const char *target)
{
#if PF
   double  f, e;

   f = pf_fold(string, structure);
   e = energy_of_structure(string, target, 0);
   return (double) (e-f-final_cost);
#else
   vrna_message_error("this version not linked with pf_fold");
   return 0;
#endif
}

/*---------------------------------------------------------------------------*/

PRIVATE char *aux_struct(const char* structure )
{
   int       *match_paren;
   int          i, o, p;
   char        *string;

   string = (char *) vrna_alloc(sizeof(char)*(strlen(structure)+1));
   match_paren = (int *) vrna_alloc(sizeof(int)*(strlen(structure)/2+1));
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
	 vrna_message_error("Junk in structure at aux_structure\n");
      }
      i++;
   }
   free(match_paren);
   return(string);
}
