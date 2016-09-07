/*
		String alignment for RNA secondary structures
		      Peter F Stadler, Ivo Hofacker
			   Vienna RNA Package
*/

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include "ViennaRNA/edit_cost.h"
#include "ViennaRNA/dist_vars.h"
#include "ViennaRNA/utils.h"


#define PUBLIC
#define PRIVATE        static

PUBLIC  float      string_edit_distance(swString *T1, swString *T2);
PUBLIC  swString  *Make_swString(char *string);
PUBLIC  void       print_swString(swString *x);
PUBLIC  void       print_alignment_list(void);

PRIVATE void       sprint_aligned_swStrings(swString *T1, swString *T2);
PRIVATE float      StrEditCost(int i, int j, swString *T1, swString *T2);
PRIVATE void       DeCode(char *string, int k, int *tp, float *w);
PRIVATE int        decode(char *id);
PRIVATE void       encode(int type, char label[]);
PRIVATE int       *alignment[2]; /* contains information from backtracking
				    alignment[0][n] is the node in tree2
				    matching node n in tree1               */


/*---------------------------------------------------------------------------*/

PUBLIC float string_edit_distance(swString *T1, swString *T2)

{
    float  **distance;
    short    **i_point, **j_point;

    int           i, j, i1, j1, pos, length1,length2;
    float         minus, plus, change, temp;

    if (cost_matrix==0) EditCost = &UsualCost;
    else EditCost = &ShapiroCost;

    length1 = T1[0].sign;
    length2 = T2[0].sign;

    distance = (float **)  vrna_alloc((length1 +1)*sizeof(float *));
    if(edit_backtrack){
       i_point  = (short **)  vrna_alloc((length1 +1)*sizeof(short *));
       j_point  = (short **)  vrna_alloc((length1 +1)*sizeof(short *));
    }
    for(i=0; i<= length1; i++){
       distance[i] = (float *) vrna_alloc( (length2+1)*sizeof(float));
       if(edit_backtrack){
	  i_point[i]  = (short *) vrna_alloc( (length2+1)*sizeof(short));
	  j_point[i]  = (short *) vrna_alloc( (length2+1)*sizeof(short));
       }
    }

    for(i = 1; i <= length1; i++) {
       if (edit_backtrack){
	  i_point[i][0] = i-1;
	  j_point[i][0] = 0;
       }
       distance[i][0] = distance[i-1][0]+StrEditCost(i,0,T1,T2);
    }
    for(j = 1; j <= length2; j++) {
       if (edit_backtrack){
	  j_point[0][j] = j-1;
	  i_point[0][j] = 0;
       }
       distance[0][j] = distance[0][j-1]+StrEditCost(0,j,T1,T2);
    }

    for (i = 1; i <= length1; i++) {
       for (j = 1; j <= length2 ; j++) {
          minus  = distance[i-1][j]  + StrEditCost(i,0,T1,T2);
          plus   = distance[i][j-1]  + StrEditCost(0,j,T1,T2);
          change = distance[i-1][j-1]+ StrEditCost(i,j,T1,T2);

          distance[i][j] = MIN3(minus, plus, change);
          /* printf("%g ", distance[i][j]); */

          if(edit_backtrack){
             if(distance[i][j] == change) {
                i_point[i][j]=i-1; j_point[i][j]=j-1;  }
             else if(distance[i][j] == plus) {
                i_point[i][j]=i  ; j_point[i][j]=j-1;  }
             else {
                i_point[i][j]=i-1; j_point[i][j]=j  ;  }
          }
       }
       /* printf("\n"); */
    }
    /* printf("\n"); */
    temp = distance[length1][length2];
    for(i=0;i<=length1;i++)
       free(distance[i]);
    free(distance);

    if(edit_backtrack){
       if(alignment[0]!= NULL) free(alignment[0]);
       if(alignment[1]!= NULL) free(alignment[1]);
       alignment[0] = (int *) vrna_alloc((length1+length2+1)*sizeof(int));
       alignment[1] = (int *) vrna_alloc((length1+length2+1)*sizeof(int));

       pos = length1+length2;
       i   = length1;
       j   = length2;
       while( (i>0)||(j>0) ) {
          i1 = i_point[i][j];
          j1 = j_point[i][j];
          if( ((i-i1)==1)&&((j-j1)==1) )  {  /* substitution    */
              alignment[0][pos] = i;
              alignment[1][pos] = j;
          }
          if( ((i-i1)==1)&&(j==j1) )      {  /* Deletion in [1] */
              alignment[0][pos] = i;
              alignment[1][pos] = 0;
          }
          if( (i==i1)&&((j-j1)==1)  )      {  /* Deletion in [0] */
              alignment[0][pos] = 0;
              alignment[1][pos] = j;
          }
          pos--;
          i = i1;
          j = j1;
       }
       for(i=pos+1; i<=length1+length2; i++){
          alignment[0][i-pos] = alignment[0][i];
          alignment[1][i-pos] = alignment[1][i];
       }
       alignment[0][0] = length1+length2-pos;   /* length of alignment */

       for(i=0; i<=length1; i++){
          free(i_point[i]); free(j_point[i]);
       }
       free(i_point); free(j_point);
       sprint_aligned_swStrings(T1,T2);

    }

    return temp;
}


/*---------------------------------------------------------------------------*/

PRIVATE float StrEditCost(int i, int j, swString *T1, swString *T2)
{
    float  c, diff, cd, min, a, b, dist;

    if(i==0) {
       cd   =  (float) (*EditCost)[0][T2[j].type];
       diff =  T2[j].weight;
       dist =  cd*diff;
    }
    else
    if(j==0) {
       cd   =  (float) (*EditCost)[T1[i].type][0];
       diff =  T1[i].weight;
       dist =  cd*diff;
    }
    else
    if( ((T1[i].sign)*(T2[j].sign)) > 0) {
       c = (float) (*EditCost)[T1[i].type][T2[j].type];
       diff = (float) fabs((a=T1[i].weight) - (b=T2[j].weight));
       min = MIN2(a,b);
       if (min == a) cd = (float) (*EditCost)[0][T2[j].type];
       else          cd = (float) (*EditCost)[T1[i].type][0];
       dist = c * min + cd * diff;
    }
    else dist = (float) DIST_INF;
    return dist;
}

/*---------------------------------------------------------------------------*/

PUBLIC swString *Make_swString(char *string)
{
   int i=0, j=0, k=0;
   int tp, len, l, length;
   float w;
   swString  *S;

   length = strlen(string);

   for(i=0; i<length; i++) {
      if( (string[i]=='(') || (string[i]==')') ) j++;
      if(string[i]=='.') j+=2;
   }

   len = j;

   S= (swString *) vrna_alloc(sizeof(swString)*(len+1));
   S[0].sign = j; /* number of entries */
   S[0].weight= 0.0;
   S[0].type=     0;

   i=0;
   j=1;
   while(i<length){
      switch(string[i]){
       case '(' :
          S[j].sign = 1;
          l=1;
          k=i;
          while (l>0) {
             k++;
             if(string[k] == '(' ) l++;
             if(string[k] == ')' ) l--;
          }
          DeCode(string,k,&tp,&w);
          S[j].type   = tp;
          S[j].weight = w/2.0;
	  j++;
          break;
       case ')' :
          k=i;
          S[j].sign = -1;
          DeCode(string,k,&tp,&w);
          S[j].type = tp;
          S[j].weight = w/2.0;
	  j++;
          break;
       case '.' :
          S[j].sign = 1;
          S[j].type = 1;
          S[j].weight = 0.5;
          j++;
          S[j].sign = -1;
          S[j].type =  1;
          S[j].weight = 0.5;
	  j++;
          break;
      }
      i++;
   }
   return S;
}

/*---------------------------------------------------------------------------*/

PRIVATE void DeCode(char *string, int k, int *tp, float *w)
   /* retrieves type and weigth for a node closed  by a bracket at position k */
{
   int i,j,l,m;
   char  label[20], id[20] ;
   i=k;
   label[0] = '\0';
   while(i>=0){
      i--;
      if( (string[i]=='(')||(string[i]==')')||(string[i]=='.') ) break;
      else {
         label[k-i-1] = string[i]; label[k-i] = '\0';
      }
   }
   l=strlen(label);
   if (l==0) {           /* Dot-Bracket notation */
     *w  = 1.0;
     *tp =   2;
   }
   else{
     for (i=0; i<l; i++) {
       if (!isalpha(label[l-i-1])) break;
       id[i] = label[l-i-1];
      }
      id[i] = '\0';
      *tp=decode(id);
      l=l-i-1;
      if(l>=0){
         for(j=0; j<=l; j++)
            id[j] = label[l-j];
         label[l+1] ='\0';
         m=-1;
         sscanf(label,"%d",&m);
         *w= (float) m;
         if(m==-1) {
            vrna_message_warning("Non-integer weight in DeCode ignored");
            *w=1.0;
         }
      }
      else
         *w=1.0;
   }
}

/*---------------------------------------------------------------------------*/

PRIVATE int decode(char *id)
{
    int   n, quit, i;
    char  label[100], *code;

    n = 0;

    quit = 0;
    code = coding;

    while (!quit) {
        for (i = 0; code[i] != sep; i++) {
            if (code[i] == '\0') {
                quit = 1;
                break;
            }
            label[i] = code[i];
        }
        label[i] = '\0';
        if (strcmp(id, label) == 0) return (n);
        code += (i+1);
        n++;
    }

    vrna_message_error("Syntax error: node identifier \"%s\" not found "
                              "in coding string \"%s\"\n"
                              "Exiting",
                              id, coding);
    exit(0);
}


/*---------------------------------------------------------------------------*/

PRIVATE void encode( int type, char label[])

{
    int   i, l;

    l = 0;
    for (i = 0; i < type; i++) {
        while (coding[l] != sep && coding[l]) l++;
        l++;
    }

    for (i = 0; coding[l+i] != sep; i++) {
        if (coding[l+i] == '\0') break;
        label[i] = coding[l+i];
    }
    label[i] = '\0';
}

/*---------------------------------------------------------------------------*/

PRIVATE void sprint_aligned_swStrings(swString *T1, swString *T2)
{
   int i, j, l0, l1, ltmp=0, weights;
   char label[10], *a0, *a1, tmp0[20], tmp1[20];

   weights = 0;
   for (i=1; i<=T1[0].sign; i++) weights = (weights||(T1[i].weight!=0.5));
   for (i=1; i<=T2[0].sign; i++) weights = (weights||(T2[i].weight!=0.5));

   a0 = (char *) vrna_alloc(alignment[0][0]*4+2);
   a1 = (char *) vrna_alloc(alignment[0][0]*4+2);
   for(i=1; i<= alignment[0][0]; i++){
      tmp0[0] = '\0'; l0=0;
      if(alignment[0][i] > 0) {
         encode(T1[alignment[0][i]].type, label);
         if(T1[alignment[0][i]].sign > 0) {
            tmp0[0] = '(';
            tmp0[1] = '\0';
         }
         strcat(tmp0,label);
	 if (weights)
	    sprintf(tmp0+strlen(tmp0), "%d",
		    (int)(2*T1[alignment[0][i]].weight));

         if(T1[alignment[0][i]].sign < 0) strcat(tmp0, ")");
	 l0 = strlen(tmp0);
      }
      tmp1[0]= '\0'; l1=0;
      if(alignment[1][i] > 0) {
         encode(T2[alignment[1][i]].type, label);
         if(T2[alignment[1][i]].sign > 0) {
            tmp1[0] = '(';
            tmp1[1] = '\0';
         }
         strcat(tmp1,label);
	 if (weights)
	    sprintf(tmp1+strlen(tmp1), "%d",
		    (int)(2*T2[alignment[1][i]].weight));

	 if(T2[alignment[1][i]].sign < 0) strcat(tmp1, ")");
         l1 = strlen(tmp1);
      }
      ltmp = MAX2(l0,l1);
      for (j=l0; j<ltmp; j++) tmp0[j] = '_';
      for (j=l1; j<ltmp; j++) tmp1[j] = '_';
      tmp0[ltmp] = '\0'; tmp1[ltmp] = '\0';

      strcat(a0,tmp0); strcat(a1,tmp1);
      ltmp = strlen(a0);
   }
   if (aligned_line[0]!= NULL) { free(aligned_line[0]); aligned_line[0]= NULL;}
   if (aligned_line[1]!= NULL) { free(aligned_line[1]); aligned_line[1]= NULL;}
   aligned_line[0] = strdup(a0);
   free(a0);
   aligned_line[1] = strdup(a1);
   free(a1);
}

/*---------------------------------------------------------------------------*/


PUBLIC void print_swString(swString *x)
{
   int i;
   for (i=0; i<=x[0].sign; i++)
      printf("(%d,%d,%f\n) ", x[i].type, x[i].sign, x[i].weight );
   printf("\n");
}

/*---------------------------------------------------------------------------*/

PUBLIC void print_alignment_list(void)
{
   int i;
   printf("\n");
   for (i=1; i<= alignment[0][0]; i++)
      printf("%3d ",alignment[0][i]);
      printf("\n");
   for (i=1; i<= alignment[0][0]; i++)
      printf("%3d ",alignment[1][i]);
   printf("\n");
}
