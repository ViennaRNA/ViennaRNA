/*
	PostScript and GML output for RNA secondary structures
		    and pair probability matrices

		 c  Ivo Hofacker and Peter F Stadler
			  Vienna RNA package
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <ctype.h>
#include "utils.h"
#include "fold_vars.h"
static char rcsid[] = "$Id: PS_dot.c,v 1.10 1998/07/19 14:10:38 ivo Exp $";

#define PUBLIC
#define  PRIVATE   static
#define  MAX(A,B)    (A)>(B)?(A):(B)
#define  PI       3.141592654
#define  PIHALF       PI/2.

PUBLIC int   gmlRNA(char *string, char *structure, char *ssfile, char option);
PUBLIC int   PS_rna_plot(char *string, char *structure, char *ssfile);
PUBLIC int   ssv_rna_plot(char *string, char *structure, char *ssfile);
PUBLIC int   PS_dot_plot(char *string, char *wastlfile);

PUBLIC int   simple_xy_coordinates(short *pair_table, float *X, float *Y);
extern int   naview_xy_coordinates(short *pair_table, float *X, float *Y);

PUBLIC int   rna_plot_type = 1;  /* 0 = simple, 1 = naview */

/* local functions */

PRIVATE void   init_aux_arrays(int length);
PRIVATE void   free_aux_arrays(void);
PRIVATE void   loop(int i, int j, short *pair_table);

/* local variables for parsing routines */

PRIVATE  float  *angle;
PRIVATE  int    *loop_size, *stack_size;
PRIVATE  int     lp, stk;

PRIVATE  float RADIUS =  15.;   /* for simple_xy_coordinates */
/*---------------------------------------------------------------------------*/

/* options for gml output: 
   uppercase letters: print sequence labels
   lowercase letters: no sequence lables
   graphics information:
   x X  simple xy plot
   (nothing else implemented at present)
   default:           no graphics data at all
   */

PUBLIC int gmlRNA(char *string, char *structure, char *ssfile, char option)
{
  FILE *gmlfile;
  int i;
  int length;
  int labels=0;
  int graphics=0;
  short *pair_table;
  float *X, *Y;
  
  if (isupper(option)) labels = 1;
  
  length = strlen(string);

  pair_table = make_pair_table(structure);

  switch(option){ 
  case 'X' :
  case 'x' :
    /* Simple XY Plot */
    X = (float *) space((length+1)*sizeof(float));
    Y = (float *) space((length+1)*sizeof(float));
    if (rna_plot_type == 0) 
      i = simple_xy_coordinates(pair_table, X, Y);
    else
      i = naview_xy_coordinates(pair_table, X, Y);
    
    if(i!=length) fprintf(stderr,"strange things happening in gmlRNA ...\n");
    break;
  default:
    /* No Graphics Information */
    X = NULL;
    Y = NULL; 
  }

  /* */
 
  gmlfile = fopen(ssfile, "w");
  if (gmlfile == NULL) {
     fprintf(stderr, "can't open file %s - not doing xy_plot\n", ssfile);
     return 0;
  }
  fprintf(gmlfile, 
	  "# Vienna RNA Package (rna2glm)\n"
          "# GML Output\n"
	  "# CreationDate: %s\n"
	  "# Name: %s\n", time_stamp(), ssfile);
  fprintf(gmlfile, 
          "graph [\n"
          " directed 0\n");
  for (i=1; i<=length; i++){
     fprintf(gmlfile, 
          " node [ id %d ", i);
     if (option) fprintf(gmlfile,    
          "label \"%c\"",string[i-1]);
     if ((option == 'X')||(option=='x'))
       fprintf(gmlfile,
	       "\n  graphics [ x %9.4f y %9.4f ]\n", X[i-1], Y[i-1]);
     fprintf(gmlfile," ]\n");
  }
  for (i=1; i<length; i++) 
    fprintf(gmlfile,
	    "edge [ source %d target %d ]\n", i, i+1);   
  for (i=1; i<=length; i++) {
     if (pair_table[i]>i) 
        fprintf(gmlfile,
		"edge [ source %d target %d ]\n", i, pair_table[i]);
  }
  fprintf(gmlfile, "]\n");
  fclose(gmlfile);

  free(pair_table);
  free(X); free(Y);
  return 1; /* success */
}

/*---------------------------------------------------------------------------*/

int PS_rna_plot(char *string, char *structure, char *ssfile)
{
  float  xmin, xmax, ymin, ymax, size;
  int    i, length;
  float *X, *Y;
  FILE  *xyplot;
  short *pair_table;

  length = strlen(string);

  pair_table = make_pair_table(structure);
  
  X = (float *) space((length+1)*sizeof(float));
  Y = (float *) space((length+1)*sizeof(float));   
  if (rna_plot_type == 0) 
    i = simple_xy_coordinates(pair_table, X, Y);
  else
    i = naview_xy_coordinates(pair_table, X, Y);
  if(i!=length) fprintf(stderr,"strange things happening in PS_rna_plot...\n");

  xmin = xmax = X[0];
  ymin = ymax = Y[0];
  for (i = 1; i <= length; i++) {
     xmin = X[i] < xmin ? X[i] : xmin;
     xmax = X[i] > xmax ? X[i] : xmax;
     ymin = Y[i] < ymin ? Y[i] : ymin;
     ymax = Y[i] > ymax ? Y[i] : ymax;
  }
  size = MAX((xmax-xmin),(ymax-ymin));
  
  xyplot = fopen(ssfile, "w");
  if (xyplot == NULL) {
     fprintf(stderr, "can't open file %s - not doing xy_plot\n", ssfile);
     return 0;
  }
  fprintf(xyplot,
	  "%%!PS-Adobe-2.0 EPSF-1.2\n"
	  "%%%%Creator: He Himself\n"
	  "%%%%CreationDate: %s"
	  "%%%%Title: Rna secondary Structure Plot\n"
	  "%%%%BoundingBox: 66 210 518 662\n"
	  "%%%%DocumentFonts: Helvetica\n"
	  "%%%%Pages: 1\n"
	  "%%%%EndComments\n\n", time_stamp());

  fprintf(xyplot,"100 dict begin\n");  /* DSC says EPS should create a dict */
  fprintf(xyplot,
	  "/fsize  {14} def\n"
	  "%% toggles: if you set them all to  false  the plot stays empty\n"
	  "/drawbases   true def  %% set to  false  to leave out sequence\n"
	  "/drawoutline true def  %% set to  false  to leave out backbone\n"
	  "/drawpairs   true def  %% set to  false  to not draw lines connecting pairs\n"
	  "%% data start here\n");
  /* sequence */
  fprintf(xyplot,"/sequence { (\\\n");  
  i=0;
  while (i<length) {
    fprintf(xyplot, "%.255s\\\n", string+i);  /* no lines longer than 255 */
    i+=255;
  }
  fprintf(xyplot,") } def\n");
  /* coordinates */
  fprintf(xyplot, "/coor [\n");
  for (i = 0; i < length; i++) 
    fprintf(xyplot, "[%3.3f %3.3f]\n", X[i], Y[i]);
  fprintf(xyplot, "] def\n");
  /* base pairs */
  fprintf(xyplot, "/pairs [\n");
  for (i = 1; i <= length; i++)
    if (pair_table[i]>i)
      fprintf(xyplot, "[%d %d]\n", i, pair_table[i]);
   fprintf(xyplot, "] def\n\n");
   /* setup */
   fprintf(xyplot,
	   "/cshow  { dup stringwidth pop fsize neg 3 div exch neg 2 div exch\n"
	   "          rmoveto show} def\n"
	   "1 setlinejoin\n"
	   "1 setlinecap\n"
	   "0.8 setlinewidth\n");

  fprintf(xyplot, "72 216 translate\n");
  fprintf(xyplot, "72 6 mul %3.3f div dup scale\n", size);
  fprintf(xyplot, "%4.3f %4.3f translate\n",
	  (size-xmin-xmax)/2, (size-ymin-ymax)/2);
  fprintf(xyplot, "/Helvetica findfont fsize scalefont setfont\n");
  /* draw the data */
  fprintf(xyplot,
	  "%% draw the outline\n"
	  "drawoutline {\n"
	  "  newpath\n"
	  "  coor 0 get aload pop 0.8 0 360 arc\n"
	  "  coor {aload pop lineto} forall\n"
	  "  stroke\n"
	  "} if\n"
	  "%% draw bases\n"
	  "drawbases {\n"
	  "  0\n"
	  "  coor {\n"
	  "    aload pop moveto\n"
	  "    dup sequence exch 1 getinterval  cshow\n"
	  "    1 add\n"
	  "  } forall\n"
	  "  pop\n"
	  "} if\n"
	  "%% draw base pairs\n"
	  "drawpairs {\n"
	  "  0.7 setlinewidth\n"
	  "  [9 3.01] 9 setdash\n"
	  "  newpath\n"
	  "  pairs {aload pop\n"
	  "     coor exch 1 sub get aload pop moveto\n"
	  "     coor exch 1 sub get aload pop lineto\n"
	  "  } forall\n"
	  "  stroke\n"
	  "} if\n");
  fprintf(xyplot, "%% show it\nshowpage\n");
  fprintf(xyplot, "end\n");
  fprintf(xyplot, "%%%%EOF\n");
  
  fclose(xyplot);

  free(pair_table);
  free(X); free(Y);
  return 1; /* success */
}

/*--------------------------------------------------------------------------*/

PUBLIC int ssv_rna_plot(char *string, char *structure, char *ssfile)
{           /* produce input for the SStructView java applet */
  FILE *ssvfile;
  int i, bp;
  int length;
  int labels=0;
  int graphics=0;
  short *pair_table;
  float *X, *Y;
  float xmin, xmax, ymin, ymax, size, xoff, yoff;
  float JSIZE = 500; /* size of the java applet window */
  float rad;
  
  length = strlen(string);
  pair_table = make_pair_table(structure);

  /* make coordinates */
  X = (float *) space((length+1)*sizeof(float));
  Y = (float *) space((length+1)*sizeof(float));
  rad = RADIUS;
  RADIUS = 10.;
  if (rna_plot_type == 0) 
    i = simple_xy_coordinates(pair_table, X, Y);
  else
    i = naview_xy_coordinates(pair_table, X, Y);
  if(i!=length) fprintf(stderr,"strange things happening in ssv_rna_plot...\n");
  RADIUS = rad;

  /* make coords nonegative */
  xmin = xmax = X[0];
  ymin = ymax = Y[0];
  for (i = 1; i <= length; i++) {
     xmin = X[i] < xmin ? X[i] : xmin;
     xmax = X[i] > xmax ? X[i] : xmax;
     ymin = Y[i] < ymin ? Y[i] : ymin;
     ymax = Y[i] > ymax ? Y[i] : ymax;
  }
  if (xmin<1) {
    for (i = 0; i <= length; i++)
      X[i] -= xmin-1;
    xmin = 1;
  }
  if (ymin<1) {
    for (i = 0; i <= length; i++)
      Y[i] -= ymin-1;
    ymin = 1;
  }
#if 0
  /* rescale coordinates, center on square of size HSIZE */
  size = MAX((xmax-xmin),(ymax-ymin));
  xoff = (size - xmax + xmin)/2;
  yoff = (size - ymax + ymin)/2;
  for (i = 0; i <= length; i++) {
    X[i] = (X[i]-xmin+xoff)*(JSIZE-10)/size + 5;
    Y[i] = (Y[i]-ymin+yoff)*(JSIZE-10)/size + 5;
  }
#endif
  /* */
 
  ssvfile = fopen(ssfile, "w");
  if (ssvfile == NULL) {
     fprintf(stderr, "can't open file %s - not doing xy_plot\n", ssfile);
     return 0;
  }
  fprintf(ssvfile, 
	  "# Vienna RNA Package (rna2ssv)\n"
          "# SStructView Output\n"
	  "# CreationDate: %s\n"
	  "# Name: %s\n", time_stamp(), ssfile);
  for (i=1; i<=length; i++)
    fprintf(ssvfile, "BASE\t%d\t%c\t%d\t%d\n",
	    i, string[i-1], (int) (X[i-1]+0.5), (int) (Y[i-1]+0.5));
  for (bp=1, i=1; i<=length; i++)
    if (pair_table[i]>i) 
      fprintf(ssvfile, "BASE-PAIR\tbp%d\t%d\t%d\n", bp++, i, pair_table[i]);
  fclose(ssvfile);

  free(pair_table);
  free(X); free(Y);
  return 1; /* success */
}

/*---------------------------------------------------------------------------*/

int PS_dot_plot(char *string, char *wastlfile)
{
   /* produce PostScript dot plot from probabilities in pr[] array */
   
   FILE *wastl;
   char name[31], *c;
   int i, j, length;
   double tmp;
   
   length= strlen(string);
   wastl = fopen(wastlfile,"w");
   if (wastl==NULL) {
     fprintf(stderr, "can't open %s for dot plot\n", wastlfile);
     return 0; /* return 0 for failure */
   }
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
	   "270 665 moveto /Times-Roman findfont 14 scalefont setfont "
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
   fprintf(wastl,"72 6 mul len 1 add div dup scale\n");
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
	   "len log 0.9 sub cvi 10 exch exp  %% grid spacing\n"
           "dup 1 gt {\n"
	   "   dup dup 20 div dup 2 array astore exch 40 div setdash\n"
	   "} { [0.3 0.7] 0.1 setdash } ifelse\n"
	   "0 exch len {\n"      /* for (i=0; i<=len; i++) */
	   "   dup dup\n"        
	   "   0 moveto\n"                     /* i 0 moveto   */
	   "   len lineto \n"                  /* i len lineto */
	   "   dup\n"
	   "   len exch sub 0 exch moveto\n"   /* 0 i moveto   */
	   "   len exch len exch sub lineto\n" /* len i lineto */
	   "   stroke\n"
	   "} for\n"
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
   return 1; /* success */
}

/*---------------------------------------------------------------------------*/

PUBLIC int simple_xy_coordinates(short *pair_table, float *x, float *y)
{
   float INIT_ANGLE=0.;     /* initial bending angle */
   float INIT_X = 100.;     /* coordinate of first digit */
   float INIT_Y = 100.;     /* see above */
   float RADIUS =  15.;

   int i, length;
   float  alpha;

   length = pair_table[0];
   init_aux_arrays(length);
   lp = stk = 0;
   loop(0, length+1, pair_table);
   loop_size[lp] -= 2;     /* correct for cheating with function loop */
  
   alpha = INIT_ANGLE;
   x[0]  = INIT_X;
   y[0]  = INIT_Y;   

   for (i = 1; i <= length; i++) {
     x[i] = x[i-1]+RADIUS*cos(alpha);
     y[i] = y[i-1]+RADIUS*sin(alpha);
     alpha += PI-angle[i+1];
   }
   free_aux_arrays();

   return length;

}

/*---------------------------------------------------------------------------*/

PRIVATE  void init_aux_arrays(int length) 
{
  angle =      (float*) space( (length+5)*sizeof(float) );
  loop_size  =   (int*) space( 16+(length/5)*sizeof(int) );
  stack_size =   (int*) space( 16+(length/5)*sizeof(int) );
}

/*---------------------------------------------------------------------------*/

PRIVATE  void free_aux_arrays(void) 
{ 
   free(angle);
   free(loop_size);  
   free(stack_size); 
}

/*---------------------------------------------------------------------------*/

PRIVATE void loop(int i, int j, short *pair_table)

             /* i, j are the positions AFTER the last pair of a stack; i.e
		i-1 and j+1 are paired. */

{
    int    count = 2;   /* counts the VERTICES of a loop polygon; that's
			   NOT necessarily the number of unpaired bases!
			   Upon entry the loop has already 2 vertices, namely
			   the pair i-1/j+1.  */

    int    r = 0, bubble = 0; /* bubble counts the unpaired digits in loops */

    int    i_old, partner, k, l, start_k, start_l, fill, ladder;
    int    begin, v, diff;
    float  polygon;

    short *remember;  

    remember = (short *) space((1+(j-i)/5)*2*sizeof(short));
    
    i_old = i-1, j++;         /* j has now been set to the partner of the
			       previous pair for correct while-loop
			       termination.  */
    while (i != j) {
	partner = pair_table[i];
	if ((!partner) || (i==0))
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
	    loop(k, l, pair_table);
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
    free(remember);
}

/*---------------------------------------------------------------------------*/
