/*
		parse and convert secondary structures
	   Walter Fontana, Ivo L Hofacker, Peter F Stadler
			Vienna RNA Package
*/

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "ViennaRNA/utils.h"
#include "ViennaRNA/RNAstruct.h"

#define PRIVATE  static
#define PUBLIC

#define MAXLEN    10000


PRIVATE char *aux_struct(const char *structure);

/* on return from parse_structure(), b2C() or b2Shapiro() ... */
PUBLIC int    loop_size[STRUC];       /* contains loop sizes of a structure */
PUBLIC int    helix_size[STRUC];      /* contains helix sizes of a structure */
PUBLIC int    loop_degree[STRUC];     /* contains loop degrees of a structure */
PUBLIC int    loops;                  /* n of loops and stacks in a structure */
PUBLIC int    unpaired, pairs;        /* n of unpaired digits and pairs */

/*---------------------------------------------------------------------------*/

PRIVATE char *aux_struct(const char* structure )
{
   short        *match_paren;
   int          i, o, p;
   char        *string;

   string = (char *) vrna_alloc(sizeof(char)*(strlen(structure)+1));
   match_paren = (short *) vrna_alloc(sizeof(short)*(strlen(structure)/2+1));
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

/*---------------------------------------------------------------------------*/

PUBLIC char *b2HIT(const char *structure)
{

   int            i, u, p, l;
   char          *string, *temp, *HIT, tt[10];

   temp = (char *) vrna_alloc(strlen(structure)*4+4);
   string = aux_struct( structure );

   strcpy(temp,"(");
   i=p=u=0; l=1;
   while (string[i]) {
      switch(string[i]) {
       case '.':
	 u++; break;
       case '[':
	 if (u>0) {
	    sprintf(tt, "(U%d)" , u);
	    strcat(temp+l, tt);
	    l+=strlen(tt);
	    u=0;
	 }
	 strcat(temp+l, "("); l++;
	 break;
       case ')':
	 if (u>0) {
	    sprintf(tt, "(U%d)" , u);
	    strcat(temp+l, tt);
	    l+=strlen(tt);
	    u=0;
	 }
	 p++;
	 break;
       case ']':
	 if (u>0) {
	    sprintf(tt, "(U%d)" , u);
	    strcat(temp+l, tt);
	    l+=strlen(tt);
	    u=0;
	 }
	 sprintf(tt,"P%d)", p+1);
	 strcat(temp+l, tt);
	 l+=strlen(tt);
	 p=0;
	 break;
      }
      i++;
   }
   if (u>0) {
      sprintf(tt, "(U%d)" , u);
      strcat(temp+l, tt);
      l+=strlen(tt);
   }
   strcat(temp+l, "R)");

   free( string );

   HIT = (char *) vrna_alloc(sizeof(char)*(strlen(temp)+2));
   strcpy(HIT, temp);
   free(temp);
   return(HIT);
}

/*---------------------------------------------------------------------------*/

PUBLIC char *b2C(const char *structure )
{
   short *bulge, *loop;

   int    i, lp, p, l;
   char  *string, *Coarse, *temp;

   bulge = (short *) vrna_alloc(sizeof(short)*(strlen(structure)/3+1));
   loop = (short *) vrna_alloc(sizeof(short)*(strlen(structure)/3+1));
   temp = (char *) vrna_alloc(4*strlen(structure)+2);

   for (i = 0; i < STRUC; i++) {
      loop_size[i] = helix_size[i] = 0;
   }
   loop_degree[0]=0;         /* open structure has degree 0 */
   pairs = unpaired = loops = lp = 0;
   loop[0]=0;

   string = aux_struct( structure );

   i=p=l=0;
   temp[l++] = '(';
   while (string[i]) {
      switch(string[i]) {
       case '.':
	 loop_size[loop[lp]]++;
	 break;
       case '[':
	 temp[l++]='(';
	 if ((i>0)&&(string[i-1]=='(')) bulge[lp]=1;
	 lp++;
	 loop_degree[++loops]=1;
	 loop[lp]=loops;
	 bulge[lp]=0;
	 break;
       case ')':
	 if (string[i-1]==']') bulge[lp]=1;
	 p++;
	 break;
       case ']':
	 if (string[i-1]==']') bulge[lp]=1;
	 switch (loop_degree[loop[lp]]) {
	  case 1:  temp[l++]='H'; break;           /* hairpin */
	  case 2:
	    if (bulge[lp]==1)
	       temp[l++] = 'B';                    /* bulge */
	    else
	       temp[l++] = 'I';                    /* internal loop */
	    break;
	  default: temp[l++] = 'M';                /* multiloop */
	 }
	 temp[l++] = ')';
	 pairs+=p+1;
	 p=0;
	 loop_degree[loop[--lp]]++;
	 break;
      }
      i++;
   }
   temp[l++] = 'R';
   temp[l++] = ')';
   temp[l]='\0';
   free(string);
   Coarse = (char *) vrna_alloc(sizeof(char)*(strlen(temp)+2));
   strcpy(Coarse, temp);
   free(temp);
   free(bulge); free(loop);
   return(Coarse);
}

/*---------------------------------------------------------------------------*/

PUBLIC char *b2Shapiro(const char *structure )
{

   short *bulge, *loop;

   int            i, lp, p, l, k;
   char          *string, *Shapiro, *temp, tt[10];

   bulge = (short *) vrna_alloc(sizeof(short)*(strlen(structure)/3+1));
   loop = (short *) vrna_alloc(sizeof(short)*(strlen(structure)/3+1));
   temp = (char *) vrna_alloc(4*strlen(structure)+3);

   for (i = 0; i < STRUC; i++) {
      loop_size[i] = helix_size[i] = 0;
   }
   loop_degree[0]=0;         /* open structure has degree 0 */
   pairs = unpaired = loops = lp = 0;
   loop[0]=0;

   string = aux_struct( structure );

   i=p=l=0;
   temp[l++] = '(';    /* root */
   while (string[i]) {
      switch(string[i]) {
       case '.':
	 unpaired++;
	 loop_size[loop[lp]]++;
	 break;
       case '[':
	 temp[l++]='(';
	 temp[l++]='(';
	 if ((i>0)&&(string[i-1]=='(' || string[i-1]=='['))
	   bulge[lp]=1;
	 lp++;
	 loop_degree[++loops]=1;
	 loop[lp]=loops;
	 bulge[lp]=0;
	 break;
       case ')':
	 if (string[i-1]==']') bulge[lp]=1;
	 p++;
	 break;
       case ']':
	 if (string[i-1]==']') bulge[lp]=1;
	 switch (loop_degree[loop[lp]]) {
	  case 1:  temp[l++]='H'; break;           /* hairpin */
	  case 2:
	    if (bulge[lp]==1)
	       temp[l++] = 'B';                    /* bulge */
	    else
	       temp[l++] = 'I';                    /* internal loop */
	    break;
	  default: temp[l++] = 'M';                /* multiloop */
	 }
	 helix_size[loop[lp]]=p+1;

         sprintf(tt, "%d)" , loop_size[loop[lp]]);
         for(k=0; k<strlen(tt); k++) temp[l++] = tt[k];
	 sprintf(tt, "S%d)" , helix_size[loop[lp]]);
         for(k=0; k<strlen(tt); k++) temp[l++] = tt[k];

	 pairs+=p+1;
	 p=0;
	 loop_degree[loop[--lp]]++;
	 break;
      }
      i++;
   }

   *tt = '\0';
   if (loop_size[0]) sprintf(tt, "E%d)" , loop_size[0]);
   strcat(tt,"R)");
   temp[l]='\0';
   strcat(temp, tt);
   Shapiro = (char *) vrna_alloc(sizeof(char)*(strlen(temp)+2));
   if (loop_size[0]) {
      Shapiro[0]='(';
      strcpy(Shapiro+1, temp);
   } else strcpy(Shapiro, temp);
   free(string);
   free(temp);
   free(loop); free(bulge);
   return Shapiro;
}




/*---------------------------------------------------------------------------*/

PUBLIC void parse_structure(const char *structure)

/*-----------------------------------------------------------------------------

    upon return from parse_structure():

    loops    ....................... number of loops or stacks in structure.
    loop_size[1 <= i <= loops] ..... size of i-th loop.
    loop_size[0] ................... number of external digits.
    loop_degree[1 <= i <= loops] ... degree (branches) of i-th loop.
    loop_degree[0] ................. number of components.
    helix_size[1 <= i <= loops] .... size of i-th stack.
    unpaired ....................... n of unpaired digits.
    pairs .......................... n of base pairs.

-----------------------------------------------------------------------------*/

{
   short  *bulge, *loop;

   int            i, lp, p;
   char          *string, *temp;

   temp = (char *)  vrna_alloc(strlen(structure)*4+2);
   bulge = (short *) vrna_alloc(sizeof(short)*(strlen(structure)/3+1));
   loop = (short *) vrna_alloc(sizeof(short)*(strlen(structure)/3+1));

   for (i = 0; i < STRUC; i++) {
      loop_size[i] = helix_size[i] = 0;
   }
   loop[0] = loop_degree[0]=0;         /* open structure has degree 0 */
   pairs = unpaired = loops = lp = 0;
   *temp='\0';

   string = aux_struct(structure);

   i=p=0;
   while (string[i]) {
      switch(string[i]) {
       case '.':
	 unpaired++;
	 loop_size[loop[lp]]++;
	 break;
       case '[':
	 if ((i>0)&&(string[i-1]=='(')) bulge[lp]=1;
	 lp++;
	 loop_degree[++loops]=1;
	 loop[lp]=loops;
	 bulge[lp]=0;
	 break;
       case ')':
	 if (string[i-1]==']') bulge[lp]=1;
	 p++;
	 break;
       case ']':
	 if (string[i-1]==']') bulge[lp]=1;
	 helix_size[loop[lp]]=p+1;
	 pairs+=p+1;
	 p=0;
	 loop_degree[loop[--lp]]++;
	 break;
      }
      i++;
   }
   free(string);
   free(bulge); free(loop);
   free(temp);
}

/*---------------------------------------------------------------------------*/

PUBLIC char *add_root(const char *structure)
{
    char *xS;
    xS = (char *) vrna_alloc(sizeof(char)*(strlen(structure)+4));
    xS[0] = '(';
    strcat(xS,structure);
    strcat(xS,"R)");
    return xS;
}


/*---------------------------------------------------------------------------*/

PUBLIC char *expand_Shapiro(const char *structure)
{
   char  *xS, *temp;
   int  i, l;

   temp = (char *) vrna_alloc(4*strlen(structure)+2);

   i = 1;
   l = 1;
   temp[0] = '(';
   while (i<strlen(structure)-1) {
      temp[l++] = structure[i];
      if      (structure[i] == '(') temp[l++] = '(';
      else if (structure[i] == ')') {
	 temp[l++] = 'S';
	 temp[l++] = ')';
      }
      i++;
   }
   temp[l++] = ')';
   temp[l] = '\0';

   xS = (char *) vrna_alloc(sizeof(char)*(strlen(temp)+1));
   strcpy(xS, temp);
   free(temp);
   return (xS);
}

/*---------------------------------------------------------------------------*/

PUBLIC char *expand_Full(const char *structure)
{
    char *xF, *temp;
    int  i, l;

    temp = (char *) vrna_alloc(4*strlen(structure)+2);

    i = 0;
    l = 0;
    while (structure[i]) {
        if      (structure[i] == '(') temp[l++] = '(';
        else if (structure[i] == ')') {
            temp[l++] = 'P';
	    temp[l++] = ')';
        }
        else {
            temp[l++] = '(';
            temp[l++] = 'U';
            temp[l++] = ')';
        }
        i++;
     }
     temp[l] = '\0';

     xF = (char *) vrna_alloc(sizeof(char)*(l+5));
     strcpy(xF, "(");
     strcat(xF, temp);
     strcat(xF, "R)");
     free(temp);
     return (xF);
}

/*---------------------------------------------------------------------------*/

PUBLIC char *unexpand_Full(const char *structure)
{
   short        *match_paren;
   char id[10], *full, *temp;
   int    i, j, k, l, o, w;

   temp = (char *) vrna_alloc(4*strlen(structure)+2);
   match_paren = (short *) vrna_alloc(sizeof(short)*(strlen(structure)/2+1));

   i = strlen(structure)-1;
   l = o = 0; k=9;
   id[9]='\0';
   while (i>=0) {
     switch (structure[i]) {
     case '(':
       for (j=0; j<match_paren[o]; j++) temp[l++]='(';
       match_paren[o--] = 0;
       break;
     case 'U':
       w=1;
       sscanf(id+k, "%d", &w);
       for (j=0; j<w; j++) temp[l++]='.';
       k=9;
       break;
     case 'P':
           w=1;
       sscanf(id+k, "%d", &w);
       for (j=0; j<w; j++) temp[l++]=')';
       match_paren[o]=w;
       k=9;
       break;
     case 'R':
       break;
     case ')':
       o++;
       break;
     default:
       id[--k]=structure[i];
     }
     i--;
   }

   temp[l] = '\0';
   full = (char *) vrna_alloc(sizeof(char)*(l+1));
   for (i=0; i<l; i++) full[i]=temp[l-i-1];
   full[l]='\0';
   free(temp);
   free(match_paren);
   return full;
}


/*---------------------------------------------------------------------------*/

PUBLIC char *unweight(const char *structure)
{
   int i,l;
   char *full, *temp;

   temp = (char *) vrna_alloc(4*strlen(structure)+1);

   i=l=0;
   while (structure[i]) {
      if (!isdigit((int)structure[i])) temp[l++]=structure[i];
      i++;
   }
   temp[l]='\0';
   full = (char *) vrna_alloc(sizeof(char)*(l+1));
   strcpy(full, temp);
   free(temp);
   return full;
}

/*---------------------------------------------------------------------------*/

PUBLIC void unexpand_aligned_F(char *align[2])
{
   char *t0, *t1;
   int i,l;

   t0 = (char *) vrna_alloc(strlen(align[0])+1);
   t1 = (char *) vrna_alloc(strlen(align[0])+1);

   for (i=0, l=0; i<strlen(align[0]); i++) {
      switch (align[0][i]) {
       case '(':
       case ')':
	 t0[l] = align[0][i];
	 t1[l++]=align[1][i];
	 break;
       case 'U':
	 switch (align[1][i]) {
	  case 'U':
	    t0[l-1]=t1[l-1]='.';
	    break;
	  case '_':
	    t0[l-1]='.';
	    t1[l-1]='_';
	    break;
	  case 'P':
	    t0[l-1]='_'; t0[l]='.';
	    t1[l-1]='('; t1[l]=')'; l++;
	 }
	 while (align[0][i]!=')') i++;
	 break;
       case '_':
	 switch (align[1][i]) {
	  case '(':
	  case ')':
	    t0[l] = align[0][i];
	    t1[l++]=align[1][i];
	    break;
	  case 'U':
	    while (align[1][i]!=')') i++;
	    t1[l-1]='.';
	    t0[l-1]='_';
	    break;
	 }
       case 'P':
	 if (align[1][i]=='U') {
	    t1[l-1]='_'; t1[l]='.'; t0[l++]=')';
	    while (align[0][i]!=')') i++;
	 }
	 break;
      }
   }
   t0[l-1]=t1[l-1]='\0';
   strcpy(align[0], t0+1);
   strcpy(align[1], t1+1);
   free(t0); free(t1);
}
