/*
	      posrscript output for secondary structures
		    and pair probability matrices
				   
		 c Ivo L Hofacker and Peter F Stadler
			  Vienna RNA package
*/
/*Last changed Time-stamp: <97/09/18 12:45:28 ivo> */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "utils.h"
#include "fold_vars.h"

#define  PRIVATE   static
#define  MAX(A,B)    (A)>(B)?(A):(B)
#define  MAXLENGTH 2000
#define  PI       3.141592654
#define  PIHALF       PI/2.

#define  BRANCH     30

PRIVATE void   make_pair_table(char *structure, short *table);
PRIVATE short  *pair_table;
/*--------------------------------------------------------------------------*/

void PS_dot_plot(char *string, char *wastlfile)
{
   /* produce PostScript dot plot from probabilities in pr[] array */
   
   FILE *wastl;
   char name[31], *c;
   int i, j, length;
   double tmp;
   
   length= strlen(string);
   wastl = fopen(wastlfile,"w");
   if (wastl==NULL) nrerror("can't open file for dot plot");
   strncpy(name, wastlfile, 30);
   if ((c=strrchr(name, '_'))!=0) *c='\0';
   fprintf(wastl,"%%!PS-Adobe-2.0 EPSF-1.2\n");
   fprintf(wastl,"%%%%Title: RNA DotPlot\n");
   fprintf(wastl,"%%%%Creator: PS_dot.c, ViennaRNA package\n");
   fprintf(wastl,"%%%%CreationDate: %s", time_stamp());
   fprintf(wastl,"%%%%BoundingBox: 66 211 518 662\n");
   fprintf(wastl,"%%%%DocumentFonts: Times-Roman\n");
   fprintf(wastl,"%%%%Pages: 1\n");
   fprintf(wastl,"%%%%EndComments\n\n");
   
   fprintf(wastl,"%%This file contains the square roots "
	   "of the base pair probabilities in the form\n");
   fprintf(wastl,"%% i  j  sqrt(p(i,j)) ubox\n");
   
   fprintf(wastl,"100 dict begin\n");  /* DSC says EPS should create a dict */
   fprintf(wastl,"%%delete next line to get rid of title\n"
	   "270 664 moveto /Times-Roman findfont 14 scalefont setfont "
	   "(%s) show\n", name);
   fprintf(wastl,"/ubox {\n");     /* upper triangle matrix */
   fprintf(wastl,"   3 1 roll\n");
   fprintf(wastl,"   exch len exch sub 1 add box\n");
   fprintf(wastl,"} def\n\n");

   fprintf(wastl,"/lbox {\n");     /* lower triangle matrix */
   fprintf(wastl,"   3 1 roll\n");
   fprintf(wastl,"   len exch sub 1 add box\n");
   fprintf(wastl,"} def\n\n");

   fprintf(wastl,"/box { %%size x y box - draws box centered on x,y\n");
   fprintf(wastl,"   2 index 0.5 mul add            %% x += 0.5\n");
   fprintf(wastl,"   exch 2 index 0.5 mul add exch  %% x += 0.5\n");
   fprintf(wastl,"   newpath\n");
   fprintf(wastl,"   moveto\n");
   fprintf(wastl,"   dup neg   0 rlineto\n");
   fprintf(wastl,"   dup neg   0 exch rlineto\n");
   fprintf(wastl,"             0 rlineto\n");
   fprintf(wastl,"   closepath\n");
   fprintf(wastl,"   fill\n");
   fprintf(wastl,"} def\n");
   fprintf(wastl,"\n");
   /* EPS should not contain lines >255 characters */
   fprintf(wastl,"/sequence { (\\\n");
   i=0;
   while (i<length) {
      fprintf(wastl, "%.255s\\\n", string+i);
      i+=255;
   }
   fprintf(wastl,") } def\n");
   fprintf(wastl,"/len { sequence length } def\n\n");
   
   fprintf(wastl,"72 216 translate\n");
   fprintf(wastl,"72 6 mul len div dup scale\n");
   fprintf(wastl,"/Times-Roman findfont 0.95 scalefont setfont\n\n");

   /* print sequence along all 4 sides */
   fprintf(wastl,"%% print sequence along all 4 sides\n");
   fprintf(wastl,"0 1 len 1 sub {\n");
   fprintf(wastl,"    dup\n");
   fprintf(wastl,"    0.7 add -0.3 moveto\n");
   fprintf(wastl,"    sequence exch 1 getinterval\n");
   fprintf(wastl,"    show\n");
   fprintf(wastl,"} for\n");
   fprintf(wastl,"\n");
   fprintf(wastl,"0 1 len 1 sub {\n");
   fprintf(wastl,"    dup\n");
   fprintf(wastl,"    0.7 add 0.7 len add moveto\n");
   fprintf(wastl,"    sequence exch 1 getinterval\n");
   fprintf(wastl,"    show\n");
   fprintf(wastl,"} for\n\n");

   fprintf(wastl,"90  rotate\n");
   fprintf(wastl,"0 1 len 1 sub {\n");
   fprintf(wastl,"    dup\n");
   fprintf(wastl,"    0.7 add -0.2 moveto\n");
   fprintf(wastl,"    len 1 sub exch sub\n");
   fprintf(wastl,"    sequence exch 1 getinterval\n");
   fprintf(wastl,"    show\n");
   fprintf(wastl,"} for\n");
   fprintf(wastl,"270 rotate\n\n");

   fprintf(wastl,"270 rotate\n");
   fprintf(wastl,"0 1 len 1 sub {\n");
   fprintf(wastl,"    dup\n");
   fprintf(wastl,"    -0.3 add len sub  0.7 len add  moveto\n");
   fprintf(wastl,"    sequence exch 1 getinterval\n");
   fprintf(wastl,"    show\n");
   fprintf(wastl,"} for\n");
   fprintf(wastl,"90 rotate\n\n");

   /* do grid */
   fprintf(wastl,"0.5 dup translate\n"
	   "%% draw diagonal\n"
	   "0.04 setlinewidth\n"
	   "0 len moveto len 0 lineto stroke \n\n");
   fprintf(wastl,"%%draw grid\n"
	   "0.01 setlinewidth\n"
	   "len log 1 sub cvi 10 exch exp  %% grid spacing\n"
           "dup dup 20 div dup 2 array astore exch 40 div setdash\n"
	   "0 exch len {\n"      /* for (i=0; i<=len; i++) */
	   "   dup dup\n"        
	   "   0 moveto\n"                     /* i 0 moveto   */
	   "   len lineto \n"                  /* i len lineto */
	   "   dup\n"
	   "   len exch sub 0 exch moveto\n"   /* 0 i moveto   */
	   "   len exch len exch sub lineto\n" /* len i lineto */
	   "} for\n"
	   "stroke\n"
	   "0.5 neg dup translate\n\n");

   /* print boxes */
   for (i=1; i<length; i++)
      for (j=i+1; j<=length; j++) {
	 if (pr[iindx[i]-j]<1e-5) continue;
	 tmp = sqrt(pr[iindx[i]-j]);
	 fprintf(wastl,"%d %d %1.5f ubox\n", i, j, tmp);
      }
   /* do mfe */
   if (base_pair) for(i=1; i<=base_pair[0].i; i++) 
      fprintf(wastl,"%d %d 0.95 lbox\n",
	      base_pair[i].i, base_pair[i].j); 

   fprintf(wastl,"showpage\n");
   fprintf(wastl,"end\n");
   fprintf(wastl,"%%%%EOF\n");
   fclose(wastl);
}

/*---------------------------------------------------------------------------*/

#if 0
PRIVATE int pairs(int i)
{
    struct bond *n;

    for (n = base_pair+1; n <= base_pair+(*base_pair).i; n++) {
	if (i == n->i)
	    return(n->j);
    }
    return(0);
}
#endif

/*---------------------------------------------------------------------------*/

PRIVATE  float angle[MAXLENGTH+5], x[MAXLENGTH+5], y[MAXLENGTH+5];
PRIVATE  int loop_size[MAXLENGTH/6], stack_size[MAXLENGTH/6], lp, stk;


/*---------------------------------------------------------------------------*/

static void loop(int i, int j)

             /* i, j are the positions AFTER the last pair of a stack; i.e
		i-1 and j+1 are paired. */

{
    int    count = 2;   /* counts the VERTICES of a loop polygon; that's
			   NOT necessarily the number of unpaired bases!
			   Upon entry the loop has already 2 vertices, namely
			   the pair i-1/j+1.  */

    int    r = 0, bubble = 0; /* bubble counts the unpaired digits in loops */

    int    i_old, partner, k, l, start_k, start_l, fill, ladder;
    int    begin, v, diff, remember[2*BRANCH];
    float  polygon;

    i_old = i-1, j++;         /* j has now been set to the partner of the
			       previous pair for correct while-loop
			       termination.  */
    while (i != j) {
	partner = pair_table[i];
	if (!partner)
	    i++, count++, bubble++;
	else {
	    count += 2;
	    k = i, l = partner;    /* beginning of stack */
	    remember[++r] = k;
	    remember[++r] = l;
	    i = partner+1;         /* next i for the current loop */

	    start_k = k, start_l = l;
	    ladder = 0;
	    do {
		k++, l--, ladder++;        /* go along the stack region */
	    }
	    while (pair_table[k] == l);

	    fill = ladder-2;
	    if (ladder >= 2) {
	       angle[start_k+1+fill] += PIHALF;   /*  Loop entries and    */
	       angle[start_l-1-fill] += PIHALF;   /*  exits get an        */
	       angle[start_k]        += PIHALF;   /*  additional PI/2.    */
	       angle[start_l]        += PIHALF;   /*  Why ? (exercise)    */
	       if (ladder > 2) {
		  for (; fill >= 1; fill--) {
		     angle[start_k+fill] = PI;  /*  fill in the angles  */
		     angle[start_l-fill] = PI;  /*  for the backbone    */
		  }
	       }
	    }
	    stack_size[++stk] = ladder;
	    loop(k, l);
	}
    }
    polygon = PI*(count-2)/(float)count; /* bending angle in loop polygon */
    remember[++r] = j;
    begin = i_old < 0 ? 0 : i_old;
    for (v = 1; v <= r; v++) {
	diff  = remember[v]-begin;
	for (fill = 0; fill <= diff; fill++)
	    angle[begin+fill] += polygon;
	if (v > r)
	    break;
	begin = remember[++v];
    }
    loop_size[++lp] = bubble;
}

/*---------------------------------------------------------------------------*/
#define PRINT_SIZES 0
static void parse(int len)
{
    int k;

    for (k = 0; k < MAXLENGTH+5; k++)
	angle[k] = 0., x[k] = 0., y[k] = 0.;

    lp = stk = 0;

    /* upon exit from loop:

		  lp-1                number of loops in the structure
       loop_size[i], i=1,...,lp-1     number of unpaired digits in the
				      i-th loop
       loop_size[lp]-2                number of external digits (free ends
				      and joins)
		  stk                 number of stacks in the structure
       stack_size[i], i=1,...,stk     number of pairs in the i-th stack
    */

    loop(0, len+1);

    loop_size[lp] -= 2;     /* correct for cheating with function loop */
#if PRINT_SIZES
    for (k = 1; k <= lp; k++)
	printf("loop %d has size %d\n", k, loop_size[k]);
    for (k = 1; k <= stk; k++)
	printf("stack %d has length %d\n", k, stack_size[k]);
#endif
}

/*---------------------------------------------------------------------------*/

#define INIT_ANGLE     0.     /* initial bending angle */
#define INIT_X       100.     /* coordinate of first digit */
#define INIT_Y       100.     /* see above */
#define RADIUS 15.

#define pairs(I)  (((I)>0)?(((I)<length)?(pair_table[(I)-1]+1):0):0)
void PS_rna_plot(char *string, char *structure, char *ssfile)
{
  float  alpha, xmin, xmax, ymin, ymax, size;
  int    i, length;
  FILE  *xyplot;

  length = strlen(string);
  if (length>MAXLENGTH) {
     fprintf(stderr,"INFO: structure too long, not doing xy_plot\n");
     return;
  }

  /* pair_table[length+1] should be 0 */
  pair_table = (short *) space((length+2)*sizeof(short));  
  make_pair_table(structure, pair_table);
  
  parse(length);
  
  alpha = INIT_ANGLE;
  x[0]  = INIT_X;
  y[0]  = INIT_Y;
  xmin = xmax = x[0];
  ymin = ymax = y[0];

  for (i = 1; i <= length; i++) {
     x[i] = x[i-1]+RADIUS*cos(alpha);
     y[i] = y[i-1]+RADIUS*sin(alpha);
     alpha += PI-angle[i+1];
     xmin = x[i] < xmin ? x[i] : xmin;
     xmax = x[i] > xmax ? x[i] : xmax;
     ymin = y[i] < ymin ? y[i] : ymin;
     ymax = y[i] > ymax ? y[i] : ymax;
  }

  size = MAX((xmax-xmin),(ymax-ymin));
  xyplot = fopen(ssfile, "w");
  if (xyplot == NULL) {
     fprintf(stderr, "can't open file %s - not doing xy_plot\n", ssfile);
     return;
  }
  fprintf(xyplot,
	  "%%!PS-Adobe-2.0 EPSF-1.2\n"
	  "%%%%Creator: He Himself\n"
	  "%%%%CreationDate: %s"
	  "%%%%Title: Rna secondary Structure Plot\n"
	  "%%%%BoundingBox: 66 210 518 662\n"
	  "%%%%DocumentFonts: Courier\n"
	  "%%%%Pages: 1\n"
	  "%%%%EndComments\n\n", time_stamp());

  fprintf(xyplot,"100 dict begin\n");  /* DSC says EPS should create a dict */
  fprintf(xyplot,
	  "/fsize  {15} def\n"
	  "/cshow  { dup stringwidth pop fsize neg 3 div exch neg 2 div exch\n"
	  "          rmoveto show} def\n"
	  "/smurgl { lineto currentpoint 3 -1 roll cshow moveto } def\n");

  fprintf(xyplot, "72 216 translate\n");
  fprintf(xyplot, "72 6 mul %3.3f div dup scale\n", size);
  fprintf(xyplot, "%4.3f %4.3f translate\n",
	  (size-xmin-xmax)/2, (size-ymin-ymax)/2);
  fprintf(xyplot, "/Courier findfont fsize scalefont setfont\n");
  fprintf(xyplot, "100 100 moveto\n");
  
  for (i = 0; i < length; i++) 
     fprintf(xyplot, "(%c) %3.3f %3.3f smurgl\n", *(string+i), x[i], y[i]);

  fprintf(xyplot, "stroke\n");
  fprintf(xyplot, "showpage\n");
  fprintf(xyplot, "end\n");
  fprintf(xyplot, "%%%%EOF\n");
  fclose(xyplot);

  free(pair_table);
}

/*---------------------------------------------------------------------------*/

PRIVATE void make_pair_table(char *structure, short *table)
{
   int i,j,hx;
   short *stack;
   
   hx=0;
   stack = (short *) space(sizeof(short)*(strlen(structure)+1));
             
   for (i=1; i<=strlen(structure); i++) {
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
}
