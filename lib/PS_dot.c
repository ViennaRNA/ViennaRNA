/*
        PostScript and GML output for RNA secondary structures
                    and pair probability matrices

                 c  Ivo Hofacker and Peter F Stadler
                          Vienna RNA package
*/
#include <config.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <ctype.h>
#include "utils.h"
#include "fold_vars.h"
#include "PS_dot.h"
#include "pair_mat.h"
#include "aln_util.h"
#include "plot_layouts.h"

static char UNUSED rcsid[] = "$Id: PS_dot.c,v 1.38 2007/02/02 15:18:13 ivo Exp $";

#ifndef PI
#define  PI       3.141592654
#endif
#define  PIHALF       PI/2.
#define SIZE 452.
#define PMIN 0.00001

/* local functions */
PRIVATE FILE  *PS_dot_common(char *seq, char *wastlfile, char *comment,
                             int winsize);
PRIVATE char **annote(const char *structure, const char *AS[]);


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
  short *pair_table;
  float *X, *Y;

  if (isupper(option)) labels = 1;

  gmlfile = fopen(ssfile, "w");
  if (gmlfile == NULL) {
     fprintf(stderr, "can't open file %s - not doing xy_plot\n", ssfile);
     return 0;
  }

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

  fprintf(gmlfile,
          "# Vienna RNA Package %s\n"
          "# GML Output\n"
          "# CreationDate: %s\n"
          "# Name: %s\n"
          "# Options: %s\n", VERSION, time_stamp(), ssfile, option_string());
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
int PS_rna_plot(char *string, char *structure, char *ssfile) {
  return PS_rna_plot_a(string, structure, ssfile, NULL, NULL);
}

static const char *RNAss_head =
"%%BeginProlog\n"
"/RNAplot 100 dict def\n"
"RNAplot begin\n"
"/fsize  14 def\n"
"/outlinecolor {0.2 setgray} bind def\n"
"/paircolor    {0.2 setgray} bind def\n"
"/seqcolor     {0   setgray} bind def\n"
"/cshow  { dup stringwidth pop -2 div fsize -3 div rmoveto show} bind def\n"
"/min { 2 copy gt { exch } if pop } bind def\n"
"/max { 2 copy lt { exch } if pop } bind def\n"
"/arccoords { % i j arccoords\n"
"  % puts optimal x1 y1 x2 y2 coordinates used in bezier curves from i to j\n"
"  % onto the stack\n"
"  dup 3 -1 roll dup 4 -1 roll lt dup dup 5 2 roll {exch} if\n"
"  dup 3 -1 roll dup 3 -1 roll exch sub 1 sub dup\n"
"  4 -2 roll 5 -1 roll {exch} if 4 2 roll\n"
"  sequence length dup 2 div exch 3 1 roll lt \n"
"  {exch 5 -1 roll pop 4 -2 roll exch 4 2 roll}\n"
"  { 4 2 roll 5 -1 roll dup 6 1 roll {exch} if\n"
"    4 -2 roll exch pop dup 3 -1 roll dup 4 1 roll\n"
"    exch add 4 -1 roll dup 5 1 roll sub 1 sub\n"
"    5 -1 roll not {4 -2 roll exch 4 2 roll} if\n"
"  }ifelse\n"
"   % compute the scalingfactor and prepare (1-sf) and sf*r\n"
"  2 mul exch cpr 3 1 roll div dup\n"
"  3 -1 roll mul exch 1 exch sub exch\n"
"   % compute the coordinates\n"
"  3 -1 roll 1 sub coor exch get aload pop % get coord for i\n"
"  4 -1 roll dup 5 1 roll mul 3 -1 roll dup 4 1 roll add exch % calculate y1\n"
"  4 -1 roll dup 5 1 roll mul 3 -1 roll dup 4 1 roll add exch % calculate x1\n"
"  5 -1 roll 1 sub coor exch get aload pop % get coord for j\n"
"  % duplicate j coord\n"
"  dup 3 -1 roll dup 4 1 roll exch 8 2 roll\n"
"  6 -1 roll dup 7 1 roll mul 5 -1 roll dup 6 1 roll add exch % calculate y2\n"
"  6 -1 roll mul 5 -1 roll add exch % calculate x2\n"
"  6 -2 roll % reorder\n"
"} bind def\n"
"/drawoutline {\n"
"  gsave outlinecolor newpath\n"
"  coor 0 get aload pop 0.8 0 360 arc % draw 5' circle of 1st sequence\n"
"  currentdict /cutpoint known        % check if cutpoint is defined\n"
"  {coor 0 cutpoint getinterval\n"
"   {aload pop lineto} forall         % draw outline of 1st sequence\n"
"   coor cutpoint 1 add get aload pop\n"
"   2 copy moveto 0.8 0 360 arc       % draw 5' circle of 2nd sequence\n"
"   coor cutpoint 1 add coor length cutpoint 1 add sub getinterval\n"
"   {aload pop lineto} forall}        % draw outline of 2nd sequence\n"
"  {coor {aload pop lineto} forall}   % draw outline as a whole\n"
"  ifelse\n"
"  stroke grestore\n"
"} bind def\n"
"/drawpairs {\n"
"  paircolor\n"
"  0.7 setlinewidth\n"
"  [9 3.01] 9 setdash\n"
"  newpath\n"
"  pairs {aload pop\n"
"      currentdict (cpr) known\n"
"      { exch dup\n"
"        coor  exch 1 sub get aload pop moveto\n"
"        exch arccoords curveto\n"
"      }\n"
"      { coor exch 1 sub get aload pop moveto\n"
"        coor exch 1 sub get aload pop lineto\n"
"      }ifelse\n"
"  } forall\n"
"  stroke\n"
"} bind def\n"
"% draw bases\n"
"/drawbases {\n"
"  [] 0 setdash\n"
"  seqcolor\n"
"  0\n"
"  coor {\n"
"    aload pop moveto\n"
"    dup sequence exch 1 getinterval cshow\n"
"    1 add\n"
"  } forall\n"
"  pop\n"
"} bind def\n\n"
"/init {\n"
"  /Helvetica findfont fsize scalefont setfont\n"
"  1 setlinejoin\n"
"  1 setlinecap\n"
"  0.8 setlinewidth\n"
"  72 216 translate\n"
"  % find the coordinate range\n"
"  /xmax -1000 def /xmin 10000 def\n"
"  /ymax -1000 def /ymin 10000 def\n"
"  coor {\n"
"      aload pop\n"
"      dup ymin lt {dup /ymin exch def} if\n"
"      dup ymax gt {/ymax exch def} {pop} ifelse\n"
"      dup xmin lt {dup /xmin exch def} if\n"
"      dup xmax gt {/xmax exch def} {pop} ifelse\n"
"  } forall\n"
"  /size {xmax xmin sub ymax ymin sub max} bind def\n"
"  72 6 mul size div dup scale\n"
"  size xmin sub xmax sub 2 div size ymin sub ymax sub 2 div\n"
"  translate\n"
"} bind def\n"
"end\n";

static const char *anote_macros =
"RNAplot begin\n"
"% extra definitions for standard anotations\n"
"/min { 2 copy gt { exch } if pop } bind def\n"
"/BLACK { 0 0 0 } def\n"
"/RED   { 1 0 0 } def\n"
"/GREEN { 0 1 0 } def\n"
"/BLUE  { 0 0 1 } def\n"
"/WHITE { 1 1 1 } def\n"
"/LabelFont { % font size LabelFont\n"
"  exch findfont exch fsize mul scalefont setfont\n"
"} bind def\n"
"/Label { % i dx dy (text) Label\n"
"  % write text at base i plus offset dx, dy\n"
"  4 3 roll 1 sub coor exch get aload pop moveto\n"
"  3 1 roll fsize mul exch fsize mul exch rmoveto\n"
"  show\n"
"} bind def\n"
"/cmark { % i cmark   draw circle around base i\n"
"  newpath 1 sub coor exch get aload pop\n"
"  fsize 2 div 0 360 arc stroke\n"
"} bind def\n"
"/gmark { % i j c gmark\n"
"  % draw basepair i,j with c counter examples in gray\n"
"  gsave\n"
"  3 min [0 0.33 0.66 0.9] exch get setgray\n"
"  1 sub dup coor exch get aload pop moveto\n"
"  sequence exch 1 getinterval cshow\n"
"  1 sub dup coor exch get aload pop moveto\n"
"  sequence exch 1 getinterval cshow\n"
"  grestore\n"
"} bind def\n"
"/segmark { % f i j lw r g b segmark\n"
"  % mark segment [i,j] with outline width lw and color rgb\n"
"  % use omark and Fomark instead\n"
"  gsave\n"
"  setrgbcolor setlinewidth\n"
"  newpath\n"
"  1 sub exch 1 sub dup\n"
"  coor exch get aload pop moveto\n"
"  currentdict (cpr) known\n"
"  {\n"
"    3 -1 roll dup 4 1 roll dup\n"
"    {\n"
"      3 1 roll dup 3 -1 roll dup\n"
"      4 1 roll exch 5 2 roll exch\n"
"    }\n"
"    {\n"
"      3 1 roll exch\n"
"    } ifelse\n"
"    1 exch { coor exch get aload pop lineto } for\n"
"    {\n"
"      dup 3 1 roll 1 add exch 1 add arccoords pop pop\n"
"      4 2 roll 5 -1 roll coor exch get aload pop curveto\n"
"    } if\n"
"  }\n"
"  {\n"
"    exch 1 exch {\n"
"      coor exch get aload pop lineto\n"
"    } for\n"
"  } ifelse\n"
"  { closepath fill } if  stroke\n"
"  grestore\n"
"} bind def\n"
"/omark { % i j lw r g b omark\n"
"  % stroke segment [i..j] with linewidth lw, color rgb\n"
"  false 7 1 roll segmark\n"
"} bind def\n"
"/Fomark { % i j r g b Fomark\n"
"  % fill segment [i..j] with color rgb\n"
"  % should precede drawbases\n"
"  1 4 1 roll true 7 1 roll segmark\n"
"} bind def\n"
"/BFmark{ % i j k l r g b BFmark\n"
"  % fill block between pairs (i,j) and (k,l) with color rgb\n"
"  % should precede drawbases\n"
"  gsave\n"
"  setrgbcolor\n"
"  newpath\n"
"  currentdict (cpr) known\n"
"  {\n"
"    dup 1 sub coor exch get aload pop moveto % move to l\n"
"    dup 1 sub 4 -1 roll dup 5 1 roll 1 sub 1 exch\n"
"    { coor exch get aload pop lineto } for % lines from l to j\n"
"    3 -1 roll 4 -1 roll dup 5 1 roll arccoords curveto % curve from j to i\n"
"    exch dup 4 -1 roll 1 sub exch 1 sub 1 exch\n"
"    { coor exch get aload pop lineto } for % lines from i to k\n"
"    exch arccoords curveto% curve from k to l\n"
"  }\n"
"  {  exch 4 3 roll exch 1 sub exch 1 sub dup\n"
"     coor exch get aload pop moveto\n"
"     exch 1 exch { coor exch get aload pop lineto } for\n"
"     exch 1 sub exch 1 sub dup\n"
"     coor exch get aload pop lineto\n"
"     exch 1 exch { coor exch get aload pop lineto } for\n"
"  } ifelse\n"
"    closepath fill stroke\n"
"   grestore\n"
"} bind def\n"
"/hsb {\n"
"  dup 0.3 mul 1 exch sub sethsbcolor\n"
"} bind def\n"
"/colorpair { % i j hue sat colorpair\n"
"  % draw basepair i,j in color\n"
"  % 1 index 0.00 ne {\n"
"  gsave\n"
"  newpath\n"
"  hsb\n"
"  fsize setlinewidth\n"
"  currentdict (cpr) known\n"
"  {\n"
"    exch dup\n"
"    coor  exch 1 sub get aload pop moveto\n"
"    exch arccoords curveto\n"
"  }\n"
"  { 1 sub coor exch get aload pop moveto\n"
"    1 sub coor exch get aload pop lineto\n"
"  } ifelse\n"
"   stroke\n"
"   grestore\n"
"   % } if\n"
"} bind def\n"
 "end\n\n";

int PS_rna_plot_a(char *string, char *structure, char *ssfile, char *pre, char *post)
{
  float  xmin, xmax, ymin, ymax, size;
  int    i, length;
  float *X, *Y;
  FILE  *xyplot;
  short *pair_table;
  char *c;

  length = strlen(string);

  xyplot = fopen(ssfile, "w");
  if (xyplot == NULL) {
    fprintf(stderr, "can't open file %s - not doing xy_plot\n", ssfile);
    return 0;
  }

  pair_table = make_pair_table(structure);

  X = (float *) space((length+1)*sizeof(float));
  Y = (float *) space((length+1)*sizeof(float));
  switch(rna_plot_type){
    case VRNA_PLOT_TYPE_SIMPLE:   i = simple_xy_coordinates(pair_table, X, Y);
                                  break;
    case VRNA_PLOT_TYPE_CIRCULAR: {
                                    int radius = 3*length;
                                    i = simple_circplot_coordinates(pair_table, X, Y);
                                    for (i = 0; i < length; i++) {
                                      X[i] *= radius;
                                      X[i] += radius;
                                      Y[i] *= radius;
                                      Y[i] += radius;
                                    }
                                  }
                                  break;
    default:                      i = naview_xy_coordinates(pair_table, X, Y);
                                  break;
  }
  if(i!=length) fprintf(stderr,"strange things happening in PS_rna_plot...\n");

  xmin = xmax = X[0];
  ymin = ymax = Y[0];
  for (i = 1; i < length; i++) {
     xmin = X[i] < xmin ? X[i] : xmin;
     xmax = X[i] > xmax ? X[i] : xmax;
     ymin = Y[i] < ymin ? Y[i] : ymin;
     ymax = Y[i] > ymax ? Y[i] : ymax;
  }
  size = MAX2((xmax-xmin),(ymax-ymin));

  fprintf(xyplot,
          "%%!PS-Adobe-3.0 EPSF-3.0\n"
          "%%%%Creator: %s, ViennaRNA-%s\n"
          "%%%%CreationDate: %s"
          "%%%%Title: RNA Secondary Structure Plot\n"
          "%%%%BoundingBox: 66 210 518 662\n"
          "%%%%DocumentFonts: Helvetica\n"
          "%%%%Pages: 1\n"
          "%%%%EndComments\n\n"
          "%%Options: %s\n", rcsid+5, VERSION, time_stamp(), option_string());
  fprintf(xyplot, "%% to switch off outline pairs of sequence comment or\n"
          "%% delete the appropriate line near the end of the file\n\n");
  fprintf(xyplot, "%s", RNAss_head);

  if (pre || post) {
    fprintf(xyplot, "%s", anote_macros);
  }
  fprintf(xyplot, "%%%%EndProlog\n");

  fprintf(xyplot, "RNAplot begin\n"
          "%% data start here\n");

  /* cut_point */
  if ((c = strchr(structure, '&'))) {
    int cutpoint;
    cutpoint = c - structure;
    string[cutpoint] = ' '; /* replace & with space */
    fprintf(xyplot, "/cutpoint %d def\n", cutpoint);
  }

  /* sequence */
  fprintf(xyplot,"/sequence (\\\n");
  i=0;
  while (i<length) {
    fprintf(xyplot, "%.255s\\\n", string+i);  /* no lines longer than 255 */
    i+=255;
  }
  fprintf(xyplot,") def\n");
  /* coordinates */
  fprintf(xyplot, "/coor [\n");
  for (i = 0; i < length; i++)
    fprintf(xyplot, "[%3.8f %3.8f]\n", X[i], Y[i]);
  fprintf(xyplot, "] def\n");
  /* correction coordinates for quadratic beziers in case we produce a circplot */
  if(rna_plot_type == VRNA_PLOT_TYPE_CIRCULAR)
    fprintf(xyplot, "/cpr %6.2f def\n", (float)3*length);
  /* base pairs */
  fprintf(xyplot, "/pairs [\n");
  for (i = 1; i <= length; i++)
    if (pair_table[i]>i)
      fprintf(xyplot, "[%d %d]\n", i, pair_table[i]);
  fprintf(xyplot, "] def\n\n");

  fprintf(xyplot, "init\n\n");
  /* draw the data */
  if (pre) {
    fprintf(xyplot, "%% Start Annotations\n");
    fprintf(xyplot, "%s\n", pre);
    fprintf(xyplot, "%% End Annotations\n");
  }
  fprintf(xyplot,
          "%% switch off outline pairs or bases by removing these lines\n"
          "drawoutline\n"
          "drawpairs\n"
          "drawbases\n");

  if (post) {
    fprintf(xyplot, "%% Start Annotations\n");
    fprintf(xyplot, "%s\n", post);
    fprintf(xyplot, "%% End Annotations\n");
  }
  fprintf(xyplot, "%% show it\nshowpage\n");
  fprintf(xyplot, "end\n");
  fprintf(xyplot, "%%%%EOF\n");

  fclose(xyplot);

  free(pair_table);
  free(X); free(Y);
  return 1; /* success */
}


int PS_rna_plot_snoop_a(char *string, char *structure, char *ssfile, int *relative_access, const char *seqs[])
{
  float  xmin, xmax, ymin, ymax, size;
  int    i, length;
  float *X, *Y;
  FILE  *xyplot;
  short *pair_table;
  short *pair_table_snoop;

  length = strlen(string);

  xyplot = fopen(ssfile, "w");
  if (xyplot == NULL) {
    fprintf(stderr, "can't open file %s - not doing xy_plot\n", ssfile);
    return 0;
  }

  pair_table = make_pair_table(structure);
  pair_table_snoop = make_pair_table_snoop(structure);

  X = (float *) space((length+1)*sizeof(float));
  Y = (float *) space((length+1)*sizeof(float));
  if (rna_plot_type == 0)
    i = simple_xy_coordinates(pair_table, X, Y);
  else
    i = naview_xy_coordinates(pair_table, X, Y);
  if(i!=length) fprintf(stderr,"strange things happening in PS_rna_plot...\n");
/*   printf("cut_point %d\n", cut_point); */
  xmin = xmax = X[0];
  ymin = ymax = Y[0];
/*   for (i = 1; i < length; i++) { */
/*     printf("%d X %f Y %f \n", i, X[i], Y[i]); */
/*     xmin = X[i] < xmin ? X[i] : xmin; */
/*     xmax = X[i] > xmax ? X[i] : xmax; */
/*     ymin = Y[i] < ymin ? Y[i] : ymin; */
/*     ymax = Y[i] > ymax ? Y[i] : ymax; */
/*   } */
  //localize centre of the interaction bucket. Geometry
  
  for (i = 1; i < cut_point; i++) {  //interior loop of size 0
    if(pair_table_snoop[i] != 0){ 
      X[i-1]=X[pair_table_snoop[i]-1]; 
      Y[i-1]=Y[pair_table_snoop[i]-1]; 
    }
    else if(pair_table_snoop[i-1] && pair_table_snoop[i+1]){ //interior loop of size 1
      X[i-1]=X[pair_table_snoop[i-1] -1-1];
      Y[i-1]=Y[pair_table_snoop[i-1] -1-1];
    } 
    else if(pair_table_snoop[i-1] && pair_table_snoop[i+2]){ //interior loop of size 2
      if(pair_table_snoop[i-1] - pair_table_snoop[i+2] ==2){
	X[i-1]=X[pair_table_snoop[i-1]-2];
	Y[i-1]=Y[pair_table_snoop[i-1]-2];
	X[i]=X[pair_table_snoop[i+2]];
	Y[i]=Y[pair_table_snoop[i+2]];
	i++;
      }
      else if(pair_table[pair_table_snoop[i-1]-1]){
	X[i-1]=X[pair_table_snoop[i-1]-2];
	Y[i-1]=Y[pair_table_snoop[i-1]-2];
	X[i]=X[pair_table[pair_table_snoop[i-1]-1]-1];
	Y[i]=Y[pair_table[pair_table_snoop[i-1]-1]-1];
	i++;
      }
      else if(pair_table[pair_table_snoop[i-1]-2]){
	X[i-1]=X[pair_table_snoop[i-1]-3];
	Y[i-1]=Y[pair_table_snoop[i-1]-3];
	X[i]=X[pair_table[pair_table_snoop[i-1]-2]-1];
	Y[i]=Y[pair_table[pair_table_snoop[i-1]-2]-1];
	i++;
      }
      else if(pair_table[pair_table_snoop[i-1]-3]){
	X[i-1]=X[pair_table_snoop[i-1]-4];
	Y[i-1]=Y[pair_table_snoop[i-1]-4];
	X[i]=X[pair_table[pair_table_snoop[i-1]-3]-1];
	Y[i]=Y[pair_table[pair_table_snoop[i-1]-3]-1];
	i++;
      }
      else{
	X[i-1]=X[pair_table_snoop[i-1]-2];
	Y[i-1]=Y[pair_table_snoop[i-1]-2];
	X[i]=X[pair_table_snoop[i+2]];
	Y[i]=Y[pair_table_snoop[i+2]];
	i++;
      }
    }
    else if(pair_table_snoop[i-1] && pair_table_snoop[i+3]){ //interior loop of size 2
      if(pair_table[pair_table_snoop[i-1]-1]){
	X[i-1]=0.5*(X[pair_table_snoop[i-1]-1]+X[pair_table_snoop[i-1]-2]);
	Y[i-1]=0.5*(Y[pair_table_snoop[i-1]-1]+Y[pair_table_snoop[i-1]-2]);
	X[i]=  0.5*(X[pair_table[pair_table_snoop[i-1]-1]-1]+X[pair_table_snoop[i-1]-2]);
	Y[i]=  0.5*(Y[pair_table[pair_table_snoop[i-1]-1]-1]+Y[pair_table_snoop[i-1]-2]);
	X[i+1]=0.5*(X[pair_table[pair_table_snoop[i-1]-1]-2]+X[pair_table[pair_table_snoop[i-1]-1]-1]);
        Y[i+1]=0.5*(Y[pair_table[pair_table_snoop[i-1]-1]-2]+Y[pair_table[pair_table_snoop[i-1]-1]-1]);
	i++;i++;

      }
      else if(pair_table[pair_table_snoop[i-1]-2]){
	X[i-1]=0.5*(X[pair_table_snoop[i-1]-2]+X[pair_table_snoop[i-1]-3]);
	Y[i-1]=0.5*(Y[pair_table_snoop[i-1]-2]+Y[pair_table_snoop[i-1]-3]);
	X[i]=  0.5*(X[pair_table[pair_table_snoop[i-1]-2]-1]+X[pair_table_snoop[i-1]-3]);
	Y[i]=  0.5*(Y[pair_table[pair_table_snoop[i-1]-2]-1]+Y[pair_table_snoop[i-1]-3]);
	X[i+1]=0.5*(X[pair_table[pair_table_snoop[i-1]-2]-2]+X[pair_table[pair_table_snoop[i-1]-2]-1]);
        Y[i+1]=0.5*(Y[pair_table[pair_table_snoop[i-1]-2]-2]+Y[pair_table[pair_table_snoop[i-1]-2]-1]);
	i++;i++;
      }
      else if(pair_table[pair_table_snoop[i-1]-3]){
	X[i-1]=0.5*(X[pair_table_snoop[i-1]-3]+X[pair_table_snoop[i-1]-4]);
	Y[i-1]=0.5*(Y[pair_table_snoop[i-1]-3]+Y[pair_table_snoop[i-1]-4]);
	X[i]=  0.5*(X[pair_table[pair_table_snoop[i-1]-3]-1]+X[pair_table_snoop[i-1]-4]);
	Y[i]=  0.5*(Y[pair_table[pair_table_snoop[i-1]-3]-1]+Y[pair_table_snoop[i-1]-4]);
	X[i+1]=0.5*(X[pair_table[pair_table_snoop[i-1]-3]-2]+X[pair_table[pair_table_snoop[i-1]-3]-1]);
        Y[i+1]=0.5*(Y[pair_table[pair_table_snoop[i-1]-3]-2]+Y[pair_table[pair_table_snoop[i-1]-3]-1]);
	i++;i++;
      }
      else{
	X[i-1]=X[pair_table_snoop[i-1]-2];
	Y[i-1]=Y[pair_table_snoop[i-1]-2];
	X[i]=X[pair_table_snoop[i-1]-2];
	Y[i]=Y[pair_table_snoop[i-1]-2];
	X[i+1]=X[pair_table_snoop[i-1]-2];
	Y[i+1]=Y[pair_table_snoop[i-1]-2];
	i++;i++;
      }
    }
  }
  double xC;
  double yC;
  double R ;
  double d ;
  float X0=-1,Y0=-1,X1=-1,Y1=-1,X2=-1,Y2=-1;
//  int c1,c2,c3;
  for(i=1;i<cut_point; i++){
    if(pair_table_snoop[i]){
      X0=X[pair_table_snoop[i]-1];Y0=Y[pair_table_snoop[i]-1];
  //    c1=pair_table_snoop[i];
      i++;
      break;
    }
  }
  for(;i<cut_point; i++){
    if(pair_table_snoop[i]){
      X1=X[pair_table_snoop[i]-1];Y1=Y[pair_table_snoop[i]-1];
    //  c2=pair_table_snoop[i];
      i++;
      break;
    }
  }
  for(;i<cut_point; i++){
    if(pair_table_snoop[i]){
      X2=X[pair_table_snoop[i]-1];Y2=Y[pair_table_snoop[i]-1];
    //  c3=pair_table_snoop[i];
      i++;
      break;
    }
  }
/*   for(i=cut_point-2;i>pair_table_snoop[c1]; i--){ */
/*     if(pair_table_snoop[i]){ */
/*       X1=X[pair_table_snoop[i]-1];Y1=Y[pair_table_snoop[i]-1]; */
/*       c2=pair_table_snoop[i]; */
/*       i++; */
/*       break; */
/*     } */
/*   } */
/*   for(i=pair_table_snoop[c1]+1;i<pair_table_snoop[c2]; i++){ */
/*     if(pair_table_snoop[i]){ */
/*       X2=X[pair_table_snoop[i]-1];Y2=Y[pair_table_snoop[i]-1]; */
/*       c3=pair_table_snoop[i]; */
/*       i++; */
/*       break; */
/*     } */
/*   } */ 
 if(X0 < 0 || X1 < 0 || X2 < 0){
   printf("Could not get the center of the binding bucket. No ps file will be produced!\n");
   fclose(xyplot);
   free(pair_table);
   free(pair_table_snoop);
   free(X);free(Y);
   pair_table=NULL;pair_table_snoop=NULL;X=NULL;Y=NULL;
   return 0;
 }
  double alpha   =   (X0 -X1)/(Y1-Y0);
  double alpha_p =   (X1 -X2)/(Y2-Y1);
  double b =         (Y0+Y1 -alpha*(X0+X1))*0.5;
  double b_p =       (Y1+Y2 -alpha_p*(X1+X2))*0.5;
  //   if(abs(alpha -alpha_p) > 0.0000001){
  xC  =  (b_p - b) / (alpha - alpha_p);
  yC  =  alpha * xC + b;
  R   =  sqrt((xC-X0)*(xC-X0) + (yC-Y0)*(yC-Y0));
  d   = sqrt((X0-X1)*(X0-X1) + (Y0-Y1)*(Y0-Y1));
   for (i = 1; i < cut_point; i++) {  
     X[i-1] = X[i-1] + 0.25*(xC-X[i-1]);  
     Y[i-1] = Y[i-1] + 0.25*(yC-Y[i-1]);  
   }  
  size = MAX2((xmax-xmin),(ymax-ymin));
  fprintf(xyplot,
	  "%%!PS-Adobe-3.0 EPSF-3.0\n"
	  "%%%%Creator: %s, ViennaRNA-%s\n"
	  "%%%%CreationDate: %s"
	  "%%%%Title: RNA Secondary Structure Plot\n"
	  "%%%%BoundingBox: 66 210 518 662\n"
	  "%%%%DocumentFonts: Helvetica\n"
	  "%%%%Pages: 1\n"
	  "%%%%EndComments\n\n"
	  "%%Options: %s\n", rcsid+5, VERSION, time_stamp(), option_string());
  fprintf(xyplot, "%% to switch off outline pairs of sequence comment or\n"
	  "%% delete the appropriate line near the end of the file\n\n");
  fprintf(xyplot, "%s", RNAss_head);
  char **A;
  fprintf(xyplot, "%s", anote_macros);
  if(seqs){
    fprintf(xyplot, "%s", anote_macros);
    A = annote(structure, (const char**) seqs);
  }
  fprintf(xyplot, "%%%%EndProlog\n");
  
  fprintf(xyplot, "RNAplot begin\n"
	  "%% data start here\n");
  /* cut_point */
  if (cut_point > 0 && cut_point <= strlen(string))
    fprintf(xyplot, "/cutpoint %d def\n", cut_point-1);
  /* sequence */
  fprintf(xyplot,"/sequence (\\\n");
  i=0;
  while (i<length) {
    fprintf(xyplot, "%.255s\\\n", string+i);  /* no lines longer than 255 */
    i+=255;
  }
  fprintf(xyplot,") def\n");
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
  for (i = 1; i <= length; i++)
    if (pair_table_snoop[i]>i)
      fprintf(xyplot, "[%d %d]\n", i, pair_table_snoop[i]);
  fprintf(xyplot, "] def\n\n");
  if(relative_access){
    fprintf(xyplot,"/S [\n");
    for(i=0;i<cut_point-1; i++){
      fprintf(xyplot, " %f\n", (float)relative_access[i]/100);
    }
    fprintf(xyplot,"]\n bind def\n");
    fprintf(xyplot,"/invert false def\n");
    fprintf(xyplot,"/range 0.8 def\n");
    fprintf(xyplot,"/drawreliability {\n"		      
		   "/Smax 2.6 def\n"			 
		   "  0	\n"			      
		   "  coor 0 cutpoint getinterval {\n"
		   "    aload pop\n"
		   "    S 3 index get\n"
		   "    Smax div range mul\n"     
		   "    invert {range exch sub} if\n"  
		   "    1 1 sethsbcolor\n"
		   "    newpath\n"
		   "    fsize 2.5 div 0 360 arc\n"
		   "    fill\n"
		   "    1 add\n"
		   "  } forall\n"
		   "\n"
	           "} bind def\n"); 
  }
  fprintf(xyplot, "init\n\n");
  /*raw the data */
  if (seqs) { 
     fprintf(xyplot, "%% Start Annotations\n"); 
     fprintf(xyplot, "%s\n", A[0]); 
     fprintf(xyplot, "%% End Annotations\n"); 
   } 


  fprintf(xyplot,"%%switch off outline pairs or bases by removing these lines\n");
  if(relative_access){
    fprintf(xyplot,"drawreliability\n");
  }
  fprintf(xyplot,
	  "drawoutline\n"
	  "drawpairs\n"
	  "drawbases\n");
  //fprintf(xyplot, "%d cmark\n",c1);
  //fprintf(xyplot, "%d cmark\n",c2);
  //fprintf(xyplot, "%d cmark\n",c3);
  if (seqs) { 
     fprintf(xyplot, "%% Start Annotations\n"); 
     fprintf(xyplot, "%s\n", A[1]); 
     fprintf(xyplot, "%% End Annotations\n"); 
   } 
  fprintf(xyplot, "%% show it\nshowpage\n");
  fprintf(xyplot, "end\n");
  fprintf(xyplot, "%%%%EOF\n");

  fclose(xyplot);
  if(seqs){free(A[0]);free(A[1]);free(A);}
  free(pair_table);free(pair_table_snoop);
  free(X); free(Y);
  return 1; /* success */
}


PRIVATE char **annote(const char *structure, const char *AS[]) {
  char *ps, *colorps, **A;
  int i, n, s, pairings, maxl;
  short *ptable;
  char * colorMatrix[6][3] = {
    {"0.0 1", "0.0 0.6",  "0.0 0.2"},  /* red    */
    {"0.16 1","0.16 0.6", "0.16 0.2"}, /* ochre  */
    {"0.32 1","0.32 0.6", "0.32 0.2"}, /* turquoise */
    {"0.48 1","0.48 0.6", "0.48 0.2"}, /* green  */
    {"0.65 1","0.65 0.6", "0.65 0.2"}, /* blue   */
    {"0.81 1","0.81 0.6", "0.81 0.2"} /* violet */
  };

  make_pair_matrix();
  n = strlen(AS[0]);
  maxl = 1024;

  A = (char **) space(sizeof(char *)*2);
  ps = (char *) space(maxl);
  colorps = (char *) space(maxl);
  ptable = alimake_pair_table(structure);
  for (i=1; i<=n; i++) {
    char pps[64], ci='\0', cj='\0';
    int j, type, pfreq[8] = {0,0,0,0,0,0,0,0}, vi=0, vj=0;
    if ((j=ptable[i])<i) continue;
    for (s=0; AS[s]!=NULL; s++) {
      type = pair[encode_char(AS[s][i-1])][encode_char(AS[s][j-1])];
      pfreq[type]++;
      if (type) {
	if (AS[s][i-1] != ci) { ci = AS[s][i-1]; vi++;}
	if (AS[s][j-1] != cj) { cj = AS[s][j-1]; vj++;}
      }
    }
    for (pairings=0,s=1; s<=7; s++) {
      if (pfreq[s]) pairings++;
    }

    if ((maxl - strlen(ps) < 192) || ((maxl - strlen(colorps)) < 64)) {
      maxl *= 2;
      ps = realloc(ps, maxl);
      colorps = realloc(colorps, maxl);
      if ((ps==NULL) || (colorps == NULL))
	  nrerror("out of memory in realloc");
    }

    if (pfreq[0]<=2) {
      snprintf(pps, 64, "%d %d %s colorpair\n",
	       i,j, colorMatrix[pairings-1][pfreq[0]]);
      strcat(colorps, pps);
    }

    if (pfreq[0]>0) {
      snprintf(pps, 64, "%d %d %d gmark\n", i, j, pfreq[0]);
      strcat(ps, pps);
    }
    if (vi>1) {
      snprintf(pps, 64, "%d cmark\n", i);
      strcat(ps, pps);
    }
    if (vj>1) {
      snprintf(pps, 64, "%d cmark\n", j);
      strcat(ps, pps);
    }
  }
  free(ptable);
  A[0]=colorps;
  A[1]=ps;
  return A;
}

/*--------------------------------------------------------------------------*/


int svg_rna_plot(char *string, char *structure, char *ssfile)
{
  float  xmin, xmax, ymin, ymax, size;
  int    i, length;
  float *X, *Y, *R = NULL, *CX = NULL, *CY = NULL;
  FILE  *xyplot;
  short *pair_table;

  length = strlen(string);

  xyplot = fopen(ssfile, "w");
  if (xyplot == NULL) {
    fprintf(stderr, "can't open file %s - not doing xy_plot\n", ssfile);
    return 0;
  }

  pair_table = make_pair_table(structure);

  X = (float *) space((length+1)*sizeof(float));
  Y = (float *) space((length+1)*sizeof(float));

  switch(rna_plot_type){
    case VRNA_PLOT_TYPE_SIMPLE:   i = simple_xy_coordinates(pair_table, X, Y);
                                  break;
    case VRNA_PLOT_TYPE_CIRCULAR: {
                                    int radius = 3*length;
                                    int dr = 0;
                                    R = (float *) space((length+1)*sizeof(float));
                                    CX = (float *) space((length+1)*sizeof(float));
                                    CY = (float *) space((length+1)*sizeof(float));
                                    i = simple_circplot_coordinates(pair_table, X, Y);
                                    for (i = 0; i < length; i++) {
                                      if(i+1 < pair_table[i+1]){
                                        dr = (pair_table[i+1]-i+1 <= (length/2 + 1)) ? pair_table[i+1]-i : i + length - pair_table[i+1];
                                        R[i] = 1. - (2.*dr/(float)length);
                                      }
                                      else if(pair_table[i+1]){
                                        R[i] = R[pair_table[i+1]-1];
                                      }
                                      else{
                                        R[i] = 1.0;
                                      }
                                      CX[i] = X[i] * radius * R[i] + radius;
                                      CY[i] = Y[i] * radius * R[i] + radius;
                                      X[i] *= radius;
                                      X[i] += radius;
                                      Y[i] *= radius;
                                      Y[i] += radius;
                                    }
                                  }
                                  break;
    default:                      i = naview_xy_coordinates(pair_table, X, Y);
                                  break;
  }

  if(i!=length) fprintf(stderr,"strange things happening in PS_rna_plot...\n");


  xmin = xmax = X[0];
  ymin = ymax = Y[0];
  for (i = 1; i < length; i++) {
     xmin = X[i] < xmin ? X[i] : xmin;
     xmax = X[i] > xmax ? X[i] : xmax;
     ymin = Y[i] < ymin ? Y[i] : ymin;
     ymax = Y[i] > ymax ? Y[i] : ymax;
  }
  for (i = 0; i < length; i++)
    Y[i] = ymin+ymax - Y[i]; /* mirror coordinates so they look as in PS */

  if(rna_plot_type == VRNA_PLOT_TYPE_CIRCULAR)
    for (i = 0; i < length; i++){
      CY[i] = ymin+ymax - CY[i]; /* mirror coordinates so they look as in PS */
    }
   
  size = MAX2((xmax-xmin),(ymax-ymin));
  size += 15; /* add some so the bounding box isn't too tight */

  fprintf(xyplot,
          "<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"yes\"?>\n"
          "<svg xmlns=\"http://www.w3.org/2000/svg\" height=\"452\" width=\"452\">\n");
  fprintf(xyplot,
          "<script type=\"text/ecmascript\">\n"
          "      <![CDATA[\n"
          "        var shown = 1;\n"
          "        function click() {\n"
          "             var seq = document.getElementById(\"seq\");\n"
          "             if (shown==1) {\n"
          "               seq.setAttribute(\"style\", \"visibility: hidden\");\n"
          "               shown = 0;\n"
          "             } else {\n"
          "               seq.setAttribute(\"style\", \"visibility: visible\");\n"
          "               shown = 1;\n"
          "             }\n"
          "         }\n"
          "        ]]>\n"
          "</script>\n");
  fprintf(xyplot,
          "  <rect style=\"stroke: white; fill: white\" height=\"452\" x=\"0\" y=\"0\" width=\"452\" onclick=\"click(evt)\" />\n"
          "  <g transform=\"scale(%7f,%7f) translate(%7f,%7f)\">\n",
          SIZE/size, SIZE/size, (size-xmin-xmax)/2, (size-ymin-ymax)/2);

  fprintf(xyplot,
          "    <polyline style=\"stroke: black; fill: none; stroke-width: 1.5\" id=\"outline\" points=\"\n");
  for (i = 0; i < length; i++)
    fprintf(xyplot, "      %3.3f,%3.3f\n", X[i], Y[i]);
  fprintf(xyplot,"    \" />\n");

  fprintf(xyplot,"    <g style=\"stroke: black; stroke-width: 1; fill: none;\" id=\"pairs\">\n");
  for (i = 1; i <= length; i++) {
    int j;
    if ((j=pair_table[i])>i){
      if(rna_plot_type == VRNA_PLOT_TYPE_CIRCULAR)
        fprintf(xyplot,
                "      <path id=\"%d,%d\" d=\"M %6.15f %6.15f C %6.15f,%6.15f %6.15f,%6.15f %6.15f %6.15f\" />\n",
                i,j, X[i-1], Y[i-1], CX[i-1], CY[i-1], CX[j-1], CY[j-1], X[j-1], Y[j-1]);
      else
        fprintf(xyplot,
                "      <line id=\"%d,%d\" x1=\"%6.5f\" y1=\"%6.5f\" x2=\"%6.5f\" y2=\"%6.5f\" />\n",
                i,j, X[i-1], Y[i-1], X[j-1], Y[j-1]);
    }
  }
  fprintf(xyplot, "    </g>\n");
  fprintf(xyplot, "    <g style=\"font-family: SansSerif\" transform=\"translate(-4.6, 4)\" id=\"seq\">\n");
  for (i = 0; i < length; i++)
    fprintf(xyplot, "      <text x=\"%.3f\" y=\"%.3f\">%c</text>\n", X[i], Y[i], string[i]);
  fprintf(xyplot, "    </g>\n");
  fprintf(xyplot, "  </g>\n");
  fprintf(xyplot, "</svg>\n");

  fclose(xyplot);

  free(pair_table);
  free(X); free(Y);
  if(R) free(R);
  if(CX) free(CX);
  if(CY) free(CY);
  return 1; /* success */
}

/*--------------------------------------------------------------------------*/

PUBLIC int ssv_rna_plot(char *string, char *structure, char *ssfile)
{           /* produce input for the SStructView java applet */
  FILE *ssvfile;
  int i, bp;
  int length;
  short *pair_table;
  float *X, *Y;
  float xmin, xmax, ymin, ymax;

  ssvfile = fopen(ssfile, "w");
  if (ssvfile == NULL) {
     fprintf(stderr, "can't open file %s - not doing xy_plot\n", ssfile);
     return 0;
  }
  length = strlen(string);
  pair_table = make_pair_table(structure);

  /* make coordinates */
  X = (float *) space((length+1)*sizeof(float));
  Y = (float *) space((length+1)*sizeof(float));

  if (rna_plot_type == 0)
    i = simple_xy_coordinates(pair_table, X, Y);
  else
    i = naview_xy_coordinates(pair_table, X, Y);
  if (i!=length)
    fprintf(stderr,"strange things happening in ssv_rna_plot...\n");

  /* make coords nonegative */
  xmin = xmax = X[0];
  ymin = ymax = Y[0];
  for (i = 1; i < length; i++) {
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
  {
    float size, xoff, yoff;
    float JSIZE = 500; /* size of the java applet window */
    /* rescale coordinates, center on square of size HSIZE */
    size = MAX2((xmax-xmin),(ymax-ymin));
    xoff = (size - xmax + xmin)/2;
    yoff = (size - ymax + ymin)/2;
    for (i = 0; i <= length; i++) {
      X[i] = (X[i]-xmin+xoff)*(JSIZE-10)/size + 5;
      Y[i] = (Y[i]-ymin+yoff)*(JSIZE-10)/size + 5;
    }
  }
#endif
  /* */

  fprintf(ssvfile,
          "# Vienna RNA Package %s\n"
          "# SStructView Output\n"
          "# CreationDate: %s\n"
          "# Name: %s\n"
          "# Options: %s\n", VERSION, time_stamp(), ssfile, option_string());
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
PUBLIC int xrna_plot(char *string, char *structure, char *ssfile)
{           /* produce input for XRNA RNA drawing program */
  FILE *ss_file;
  int i;
  int length;
  short *pair_table;
  float *X, *Y;

  ss_file = fopen(ssfile, "w");
  if (ss_file == NULL) {
    fprintf(stderr, "can't open file %s - not doing xy_plot\n", ssfile);
    return 0;
  }

  length = strlen(string);
  pair_table = make_pair_table(structure);

  /* make coordinates */
  X = (float *) space((length+1)*sizeof(float));
  Y = (float *) space((length+1)*sizeof(float));

  if (rna_plot_type == 0)
    i = simple_xy_coordinates(pair_table, X, Y);
  else
    i = naview_xy_coordinates(pair_table, X, Y);
  if (i!=length)
    fprintf(stderr,"strange things happening in xrna_plot...\n");

  fprintf(ss_file,
          "# Vienna RNA Package %s, XRNA output\n"
          "# CreationDate: %s\n"
          "# Options: %s\n", VERSION, time_stamp(), option_string());
  for (i=1; i<=length; i++)
    /* XRNA likes to have coordinate mirrored, so we use (-X, Y) */
    fprintf(ss_file, "%d %c %6.2f %6.2f %d %d\n", i, string[i-1],
            -X[i-1], Y[i-1], (pair_table[i]?1:0), pair_table[i]);
  fclose(ss_file);

  free(pair_table);
  free(X); free(Y);
  return 1; /* success */
}


/*---------------------------------------------------------------------------*/
const char *RNAdp_prolog =
"%This file contains the square roots of the base pair probabilities in the form\n"
"% i  j  sqrt(p(i,j)) ubox\n\n"
"%%BeginProlog\n"
"/DPdict 100 dict def\n"
"DPdict begin\n"
"/logscale false def\n"
"/lpmin 1e-05 log def\n\n"
"/box { %size x y box - draws box centered on x,y\n"
"   2 index 0.5 mul sub            % x -= 0.5\n"
"   exch 2 index 0.5 mul sub exch  % y -= 0.5\n"
"   3 -1 roll dup rectfill\n"
"} bind def\n\n"
"/ubox {\n"
"   logscale {\n"
"      log dup add lpmin div 1 exch sub dup 0 lt { pop 0 } if\n"
"   } if\n"
"   3 1 roll\n"
"   exch len exch sub 1 add box\n"
"} bind def\n\n"
"/lbox {\n"
"   3 1 roll\n"
"   len exch sub 1 add box\n"
"} bind def\n\n"
"/drawseq {\n"
"% print sequence along all 4 sides\n"
"[ [0.7 -0.3 0 ]\n"
"  [0.7 0.7 len add 0]\n"
"  [-0.3 len sub -0.4 -90]\n"
"  [-0.3 len sub 0.7 len add -90]\n"
"] {\n"
"   gsave\n"
"    aload pop rotate translate\n"
"    0 1 len 1 sub {\n"
"     dup 0 moveto\n"
"     sequence exch 1 getinterval\n"
"     show\n"
"    } for\n"
"   grestore\n"
"  } forall\n"
"} bind def\n\n"
"/drawgrid{\n"
"  0.01 setlinewidth\n"
"  len log 0.9 sub cvi 10 exch exp  % grid spacing\n"
"  dup 1 gt {\n"
"     dup dup 20 div dup 2 array astore exch 40 div setdash\n"
"  } { [0.3 0.7] 0.1 setdash } ifelse\n"
"  0 exch len {\n"
"     dup dup\n"
"     0 moveto\n"
"     len lineto \n"
"     dup\n"
"     len exch sub 0 exch moveto\n"
"     len exch len exch sub lineto\n"
"     stroke\n"
"  } for\n"
"  [] 0 setdash\n"
"  0.04 setlinewidth \n"
"  currentdict /cutpoint known {\n"
"    cutpoint 1 sub\n"
"    dup dup -1 moveto len 1 add lineto\n"
"    len exch sub dup\n"
"    -1 exch moveto len 1 add exch lineto\n"
"    stroke\n"
"  } if\n"
"  0.5 neg dup translate\n"
"} bind def\n\n"
"end\n"
"%%EndProlog\n";


/*---------------------------------------------------------------------------*/
int PS_color_dot_plot(char *seq, cpair *pi, char *wastlfile) {
  /* produce color PostScript dot plot from cpair */

  FILE *wastl;
  int i, length;

  length= strlen(seq);
  wastl = PS_dot_common(seq, wastlfile, NULL, 0);
  if (wastl==NULL)  return 0; /* return 0 for failure */

  fprintf(wastl, "/hsb {\n"
          "dup 0.3 mul 1 exch sub sethsbcolor\n"
          "} bind def\n\n");

  /* print boxes */
   i=0;
   while (pi[i].j>0) {
     fprintf(wastl,"%1.2f %1.2f hsb %d %d %1.6f ubox\n",
             pi[i].hue, pi[i].sat, pi[i].i, pi[i].j, sqrt(pi[i].p));

     if (pi[i].mfe)
       fprintf(wastl,"%1.2f %1.2f hsb %d %d %1.4f lbox\n",
               pi[i].hue, pi[i].sat, pi[i].i, pi[i].j, pi[i].p);
     i++;
   }

   fprintf(wastl,"showpage\n"
           "end\n"
           "%%%%EOF\n");
   fclose(wastl);
   return 1; /* success */
}

/*---------------------------------------------------------------------------*/



PUBLIC int PS_dot_plot_list(char *seq, char *wastlfile, plist *pl, plist *mf, char *comment) {
  FILE *wastl;
  int length;
  double tmp;
  struct plist *pl1;

  length= strlen(seq);
  wastl = PS_dot_common(seq, wastlfile, comment, 0);
  if (wastl==NULL) return 0; /* return 0 for failure */

  fprintf(wastl,"%%data starts here\n");
  /* print boxes in upper right half*/
  for (pl1=pl; pl1->i>0; pl1++) {
    tmp = sqrt(pl1->p);
    fprintf(wastl,"%d %d %1.9f ubox\n", pl1->i, pl1->j, tmp);
  }

  /* print boxes in lower left half (mfe) */
  for (pl1=mf; pl1->i>0; pl1++) {
    tmp = sqrt(pl1->p);
    fprintf(wastl,"%d %d %1.7f lbox\n", pl1->i, pl1->j, tmp);
  }

  fprintf(wastl,"showpage\n"
          "end\n"
          "%%%%EOF\n");
  fclose(wastl);
  return 1; /* success */
}

const char *RNAdp_prolog_turn =
"/drawseq_turn {"
"% print sequence at bottom\n"
"   gsave\n"
"   len 2 sqrt div dup neg 0.28 add exch 0.78 sub translate\n"
"    0 1 len 1 sub {\n"
"     dup dup 2 sqrt mul 0 moveto\n"
"     sequence exch 1 getinterval\n"
"     show\n"
"    } for\n"
"   grestore\n"
"} bind def\n"
"/drawgrid_turn{\n"
"  0.01 setlinewidth\n"
"  len log 0.9 sub cvi 10 exch exp  % grid spacing\n"
"  dup 1 gt {\n"
"     dup dup 20 div dup 2 array astore exch 40 div setdash\n"
"  } { [0.3 0.7] 0.1 setdash } ifelse\n"
"  0 exch len {    %for (0, gridspacing, len) \n"
"     dup dup      %duplicate what - gridspacing??\n"
"     dup len exch sub moveto     %moveto diagonal?\n"
"     dup winSize gt\n"
"     {dup dup len exch sub winSize add lineto}\n"
"     {dup len lineto}ifelse\n"
"     dup len exch sub moveto  %moveto diagonal?\n"
"     dup len winSize sub le\n"
"     {dup dup len exch sub dup winSize exch sub len add exch lineto}\n"
"     {dup dup len exch sub len exch lineto}ifelse"
"     stroke pop pop\n"
"  } for\n"
"  len log 0.9 sub cvi 10 exch exp  % grid spacing\n"
"      dup 1 gt {\n"
"          dup dup 20 div dup 2 array astore exch 40 div setdash\n"
"      } { [0.3 0.7] 0.1 setdash } ifelse\n"
"      0 exch len {    %for (0, gridspacing, len) \n"
"     dup dup      %duplicate what - gridspacing??\n"
"     dup len exch sub moveto     %moveto diagonal?\n"
"     len exch sub 0.7 sub exch 0.7 sub exch lineto\n"
"     stroke\n"
"   }for\n"
" winSize len moveto  len winSize  lineto stroke\n"
"  [] 0 setdash\n"
"  0.04 setlinewidth \n"
"  currentdict /cutpoint known {\n"
"    cutpoint 1 sub\n"
"    dup dup -1 moveto len 1 add lineto\n"
"    len exch sub dup\n"
"    -1 exch moveto len 1 add exch lineto\n"
"   stroke\n"
"  } if\n"
"  0.5 neg dup translate\n"
  "} bind def \n\n";

int PS_color_dot_plot_turn(char *seq, cpair *pi, char *wastlfile, int winSize) {
  /* produce color PostScript dot plot from cpair */

  FILE *wastl;
  int i, length;

  length= strlen(seq);
  wastl = PS_dot_common(seq, wastlfile, NULL, winSize);
  if (wastl==NULL)
    return 0; /* return 0 for failure */

  fprintf(wastl, "/hsb {\n"
          "dup 0.3 mul 1 exch sub sethsbcolor\n"
          "} bind def\n\n"
          "%%BEGIN DATA\n");

  /* print boxes */
   i=0;
   while (pi[i].j>0) {
     fprintf(wastl,"%1.2f %1.2f hsb %d %d %1.6f ubox\n",
             pi[i].hue, pi[i].sat, pi[i].i, pi[i].j, sqrt(pi[i].p));/*sqrt??*/

     if (pi[i].mfe)
       fprintf(wastl,"%1.2f %1.2f hsb %d %d %1.4f lbox\n",
               pi[i].hue, pi[i].sat, pi[i].i, pi[i].j, pi[i].p);
     i++;
   }

   fprintf(wastl,"showpage\n"
           "end\n"
           "%%%%EOF\n");
   fclose(wastl);
   return 1; /* success */
}

int PS_dot_plot_turn(char *seq, struct plist *pl, char *wastlfile, int winSize) {
  /* produce color PostScript dot plot from cpair */

  FILE *wastl;
  int i, length;

  length= strlen(seq);
  wastl = PS_dot_common(seq, wastlfile, NULL, winSize);
  if (wastl==NULL)
    return 0; /* return 0 for failure */

  /* print boxes */
   i=0;
   while (pl[i].j>0) {
     fprintf(wastl,"%d %d %1.4f ubox\n",
              pl[i].i, pl[i].j, sqrt(pl[i].p));
     i++;
   }

   fprintf(wastl,"showpage\n"
           "end\n"
           "%%%%EOF\n");
   fclose(wastl);
   return 1; /* success */
}

static FILE * PS_dot_common(char *seq, char *wastlfile,
                            char *comment, int winsize) {
  /* write PS header etc for all dot plot variants */
  FILE *wastl;
  char name[31], *c;
  int i, length;

  length= strlen(seq);
  wastl = fopen(wastlfile,"w");
  if (wastl==NULL) {
    fprintf(stderr, "can't open %s for dot plot\n", wastlfile);
    return NULL; /* return 0 for failure */
  }
  strncpy(name, wastlfile, 30);
  if ((c=strrchr(name, '_'))!=0) *c='\0';

  fprintf(wastl,
          "%%!PS-Adobe-3.0 EPSF-3.0\n"
          "%%%%Title: RNA Dot Plot\n"
          "%%%%Creator: %s, ViennaRNA-%s\n"
          "%%%%CreationDate: %s", rcsid+5, VERSION, time_stamp());
  if (winsize>0)
    fprintf(wastl, "%%%%BoundingBox: 66 530 520 650\n");
  else
    fprintf(wastl, "%%%%BoundingBox: 66 211 518 662\n");
  fprintf(wastl,
          "%%%%DocumentFonts: Helvetica\n"
          "%%%%Pages: 1\n"
          "%%%%EndComments\n\n"
          "%%Options: %s\n", option_string());

  if (comment) fprintf(wastl,"%% %s\n",comment);

  fprintf(wastl,"%s", RNAdp_prolog);

  fprintf(wastl,"DPdict begin\n"
          "%%delete next line to get rid of title\n"
          "270 665 moveto /Helvetica findfont 14 scalefont setfont "
          "(%s) show\n\n", name);

  fprintf(wastl,"/sequence { (\\\n");
  for (i=0; i<strlen(seq); i+=255)
    fprintf(wastl, "%.255s\\\n", seq+i);
  fprintf(wastl,") } def\n");
  if (winsize>0)
    fprintf(wastl,"/winSize %d def\n",winsize);
  fprintf(wastl,"/len { sequence length } bind def\n\n");
  if (cut_point>0) fprintf(wastl,"/cutpoint %d def\n\n", cut_point);


  if (winsize>0)
  fprintf(wastl,"292 416 translate\n"
          "72 6 mul len 1 add winSize add 2 sqrt mul div dup scale\n");
  else
    fprintf(wastl,"72 216 translate\n"
          "72 6 mul len 1 add div dup scale\n");
  fprintf(wastl, "/Helvetica findfont 0.95 scalefont setfont\n\n");

  if (winsize>0) {
    fprintf(wastl, "%s", RNAdp_prolog_turn);
    fprintf(wastl,"0.5 dup translate\n"
          "drawseq_turn\n"
          "45 rotate\n"
          "drawgrid_turn\n");
  }
  else
    fprintf(wastl,"drawseq\n"
            "0.5 dup translate\n"
            "%% draw diagonal\n"
            "0.04 setlinewidth\n"
            "0 len moveto len 0 lineto stroke \n\n"
            "drawgrid\n");
  return(wastl);
}

int PS_color_aln(const char *structure, const char *filename,
                 const char *seqs[], const char *names[]) {
  /* produce PS sequence alignment color-annotated by consensus structure */

  int N,i,j,k,x,y,tmp,columnWidth;
  char *tmpBuffer,*ssEscaped,*ruler, *cons;
  char c;
  float fontWidth, fontHeight, imageHeight, imageWidth,tmpColumns;
  int length, maxName, maxNum, currPos;
  float lineStep,blockStep,consStep,ssStep,rulerStep,nameStep,numberStep;
  float maxConsBar,startY,namesX,seqsX, currY;
  float score,barHeight,xx,yy;
  int match,block;
  FILE *outfile;
  short *pair_table;
  char * colorMatrix[6][3] = {
    {"0.0 1", "0.0 0.6",  "0.0 0.2"},  /* red    */
    {"0.16 1","0.16 0.6", "0.16 0.2"}, /* ochre  */
    {"0.32 1","0.32 0.6", "0.32 0.2"}, /* turquoise */
    {"0.48 1","0.48 0.6", "0.48 0.2"}, /* green  */
    {"0.65 1","0.65 0.6", "0.65 0.2"}, /* blue   */
    {"0.81 1","0.81 0.6", "0.81 0.2"} /* violet */
  };

  const char *alnPlotHeader =
        "%%!PS-Adobe-3.0 EPSF-3.0\n"
        "%%%%BoundingBox: %i %i %i %i\n"
        "%%%%EndComments\n"
        "%% draws Vienna RNA like colored boxes\n"
        "/box { %% x1 y1 x2 y2 hue saturation\n"
        "  gsave\n"
        "  dup 0.3 mul 1 exch sub sethsbcolor\n"
        "  exch 3 index sub exch 2 index sub rectfill\n"
        "  grestore\n"
        "} def\n"
        "%% draws a box in current color\n"
        "/box2 { %% x1 y1 x2 y2\n"
        "  exch 3 index sub exch 2 index sub rectfill\n"
        "} def\n"
        "/string { %% (Text) x y\n"
        " 6 add\n"
        " moveto\n"
        "  show\n"
        "} def\n"
        "0 %i translate\n"
        "1 -1 scale\n"
        "/Courier findfont\n"
        "[10 0 0 -10 0 0] makefont setfont\n";


  outfile = fopen(filename, "w");

  if (outfile == NULL) {
    fprintf(stderr, "can't open file %s - not doing alignment plot\n",
            filename);
    return 0;
  }

  columnWidth=60;            /* Display long alignments in blocks of this size */
  fontWidth=6;               /* Font metrics */
  fontHeight=6.5;
  lineStep=fontHeight+2;     /* distance between lines */
  blockStep=3.5*fontHeight;  /* distance between blocks */
  consStep=fontHeight*0.5;   /* distance between alignment and conservation curve */
  ssStep=2;                  /* distance between secondary structure line and sequences */
  rulerStep=2;               /* distance between sequences and ruler */
  nameStep=3*fontWidth;             /* distance between names and sequences */
  numberStep=fontWidth;      /* distance between sequeces and numbers */
  maxConsBar=2.5*fontHeight; /* Height of conservation curve */
  startY=2;                     /* "y origin" */
  namesX=fontWidth;             /* "x origin" */

  /* Number of columns of the alignment */
  length=strlen(seqs[0]);

  /* Allocate memory for various strings, length*2 is (more than)
         enough for all of them */
  tmpBuffer = (char *) space((unsigned) MAX2(length*2,columnWidth)+1);
  ssEscaped=(char *) space((unsigned) length*2);
  ruler=(char *) space((unsigned) length*2);

  pair_table=make_pair_table(structure);
  /* Get length of longest name and count sequences in alignment*/

  for (i=maxName=N=0; names[i] != NULL; i++) {
    N++;
    tmp=strlen(names[i]);
    if (tmp>maxName)  maxName=tmp;
  }


  /* x-coord. where sequences start */
  seqsX=namesX+maxName*fontWidth+nameStep;

  /* calculate number of digits of the alignment length */
  snprintf(tmpBuffer,length, "%i",length);
  maxNum=strlen(tmpBuffer);


  /* Calculate bounding box */
  tmpColumns=columnWidth;
  if (length<columnWidth){
        tmpColumns=length;
  }
  imageWidth=ceil(namesX+(maxName+tmpColumns+maxNum)*fontWidth+2*nameStep+fontWidth+numberStep);
  imageHeight=startY+ceil((float)length/columnWidth)*((N+2)*lineStep+blockStep+consStep+ssStep+rulerStep);

  /* Write postscript header including correct bounding box */
  fprintf(outfile,alnPlotHeader,0,0,(int)imageWidth,(int)imageHeight,(int)imageHeight);

  /* Create ruler and secondary structure lines */
  i=0;
  /* Init all with dots */
  for (i=0;i<(length);i++){
        ruler[i]='.';
  }
  i=0;
  for (i=0;i<length;i++){
        /* Write number every 10th position, leave out block breaks */
        if ((i+1)%10==0 && (i+1)%columnWidth!=0){
          snprintf(tmpBuffer,length,"%i",i+1);
          strncpy(ruler+i,tmpBuffer,strlen(tmpBuffer));
        }
  }
  ruler[length]='\0';

  /* Draw color annotation first */
  /* Repeat for all pairs */
  for (i=1; i<=length; i++) {
    if ((j=pair_table[i])>i) {
      /* Repeat for open and closing position */
      for (k=0;k<2;k++){
        int pairings, nonpair, s, col;
        int ptype[8] = {0,0,0,0,0,0,0,0};
        char *color;
        col = (k==0)?i-1:j-1;
        block=ceil((float)(col+1)/columnWidth);
        xx=seqsX+(col-(block-1)*columnWidth)*fontWidth;
        /* Repeat for each sequence */
        for (s=pairings=nonpair=0; s<N; s++) {
          ptype[BP_pair[ENCODE(seqs[s][i-1])][ENCODE(seqs[s][j-1])]]++;
        }
        for (pairings=0,s=1; s<=7; s++) {
          if (ptype[s]) pairings++;
        }
        nonpair=ptype[0];
        if (nonpair <=2) {
          color = colorMatrix[pairings-1][nonpair];
          for (s=0; s<N; s++) {
            yy=startY+(block-1)*(lineStep*(N+2)+blockStep+consStep+rulerStep)+ssStep*(block)+(s+1)*lineStep;

            /* Color according due color information in pi-array, only if base pair is possible */
            if (BP_pair[ENCODE(seqs[s][i-1])][ENCODE(seqs[s][j-1])]) {

              fprintf(outfile, "%.1f %.1f %.1f %.1f %s box\n",
                      xx,yy-1,xx+fontWidth,yy+fontHeight+1,color);
            }
          }
        }
      }
    }
  }
  free(pair_table);

  /* Process rest of the output in blocks of columnWidth */

  currY=startY;
  currPos=0;

  cons =  consensus(seqs);

  while (currPos<length) {

    /* Display secondary structure line */
    fprintf(outfile,"0 setgray\n");
    strncpy(tmpBuffer,structure+currPos,columnWidth);
    tmpBuffer[columnWidth]='\0';

    x=0;y=0;
    while ((c=tmpBuffer[x])){
      if (c=='.'){
        ssEscaped[y++]='.';
      } else {
        ssEscaped[y++]='\\';
        ssEscaped[y++]=c;
      }
      x++;
    }
    ssEscaped[y]='\0';

    fprintf(outfile, "(%s) %.1f %.1f string\n", ssEscaped,seqsX,currY);
    currY+=ssStep+lineStep;

    /* Display names, sequences and numbers */

    for (i=0; i<N; i++) {

      strncpy(tmpBuffer,seqs[i]+currPos,columnWidth);
      tmpBuffer[columnWidth]='\0';

      match=0;
      for (j=0;j<(currPos+strlen(tmpBuffer));j++){
        if (seqs[i][j] != '-') match++;
      }

      fprintf(outfile, "(%s) %.1f %.1f string\n", names[i],namesX,currY);
      fprintf(outfile, "(%s) %.1f %.1f string\n", tmpBuffer,seqsX,currY);
      fprintf(outfile, "(%i) %.1f %.1f string\n", match,seqsX+fontWidth*(strlen(tmpBuffer))+numberStep,currY);
      currY+=lineStep;
    }
    currY+=rulerStep;
    strncpy(tmpBuffer,ruler+currPos,columnWidth);
    tmpBuffer[columnWidth]='\0';
    fprintf(outfile, "(%s) %.1f %.1f string\n", tmpBuffer,seqsX,currY);

    currY+=lineStep;
    currY+=consStep;

    /*Display conservation bar*/

    fprintf(outfile,"0.6 setgray\n");
    for (i=currPos;(i<currPos+columnWidth && i<length);i++){
      match=0;
      for (j=0;j<N;j++){
        if (cons[i] == seqs[j][i]) match++;
        if (cons[i]=='U' && seqs[j][i]=='T') match++;
        if (cons[i]=='T' && seqs[j][i]=='U') match++;
      }
      score=(float)(match-1)/(N-1);

      if (cons[i] == '-' ||
          cons[i] == '_' ||
          cons[i] == '.'){
        score=0;
      }

      barHeight=maxConsBar*score;
      if (barHeight==0){
        barHeight=1;
      }

      xx=seqsX+(i-(columnWidth*currPos/columnWidth))*fontWidth;

      fprintf(outfile,"%.1f %.1f %.1f %.1f box2\n",
              xx,
              currY+maxConsBar-barHeight,
              xx+fontWidth,
              currY+maxConsBar);
    }

    currY+=blockStep;
    currPos+=columnWidth;
  }
  free(cons);

  fprintf(outfile,"showpage\n");
  fclose(outfile);

  free(tmpBuffer);
  free(ssEscaped);free(ruler);

  return 0;

}

int aliPS_color_aln(const char *structure, const char *filename, 
		    const char *seqs[], const char *names[]) {
  /* produce PS sequence alignment color-annotated by consensus structure */

  int N,i,j,k,x,y,tmp,columnWidth;
  char *tmpBuffer,*ssEscaped,*ruler, *cons;
  char c;
  float fontWidth, fontHeight, imageHeight, imageWidth,tmpColumns;
  int length, maxName, maxNum, currPos;
  float lineStep,blockStep,consStep,ssStep,rulerStep,nameStep,numberStep;
  float maxConsBar,startY,namesX,seqsX, currY;
  float score,barHeight,xx,yy;
  int match,block;
  FILE *outfile;
  short *pair_table;
  char * colorMatrix[6][3] = {
    {"0.0 1", "0.0 0.6",  "0.0 0.2"},  /* red    */
    {"0.16 1","0.16 0.6", "0.16 0.2"}, /* ochre  */
    {"0.32 1","0.32 0.6", "0.32 0.2"}, /* turquoise */
    {"0.48 1","0.48 0.6", "0.48 0.2"}, /* green  */
    {"0.65 1","0.65 0.6", "0.65 0.2"}, /* blue   */
    {"0.81 1","0.81 0.6", "0.81 0.2"} /* violet */
  };

  const char *alnPlotHeader =
	"%%!PS-Adobe-3.0 EPSF-3.0\n"
	"%%%%BoundingBox: %i %i %i %i\n"
	"%%%%EndComments\n"
	"%% draws Vienna RNA like colored boxes\n"
	"/box { %% x1 y1 x2 y2 hue saturation\n"
	"  gsave\n"
	"  dup 0.3 mul 1 exch sub sethsbcolor\n"
	"  exch 3 index sub exch 2 index sub rectfill\n"
	"  grestore\n"
	"} def\n"
	"%% draws a box in current color\n"
	"/box2 { %% x1 y1 x2 y2\n"
	"  exch 3 index sub exch 2 index sub rectfill\n"
	"} def\n"
	"/string { %% (Text) x y\n"
	" 6 add\n"
	" moveto\n"
	"  show\n"
	"} def\n"
	"0 %i translate\n"
	"1 -1 scale\n"
	"/Courier findfont\n"
	"[10 0 0 -10 0 0] makefont setfont\n";
	
  
  outfile = fopen(filename, "w");
  if (outfile == NULL) {
    fprintf(stderr, "can't open file %s - not doing alignment plot\n", 
	    filename);
    return 0;
  }
  
  columnWidth=100;            /* Display long alignments in blocks of this size */
  fontWidth=6;               /* Font metrics */
  fontHeight=6.5;
  lineStep=fontHeight+2;     /* distance between lines */
  blockStep=3.5*fontHeight;  /* distance between blocks */
  consStep=fontHeight*0.5;   /* distance between alignment and conservation curve */
  ssStep=2;                  /* distance between secondary structure line and sequences */
  rulerStep=2;               /* distance between sequences and ruler */
  nameStep=3*fontWidth;	     /* distance between names and sequences */
  numberStep=fontWidth;      /* distance between sequeces and numbers */
  maxConsBar=2.5*fontHeight; /* Height of conservation curve */
  startY=2;		     /* "y origin" */
  namesX=fontWidth;	     /* "x origin" */

  /* Number of columns of the alignment */
  length=strlen(seqs[0]);

  /* Allocate memory for various strings, length*2 is (more than)
	 enough for all of them */
  tmpBuffer = (char *) space((unsigned) columnWidth + length*2 );
  ssEscaped=(char *) space((unsigned) length*2 );
  ruler=(char *) space((unsigned) length*2  );
/*   char * structur; */
/*   structur = (char*) space((length+1)*sizeof(char)); */
/*   structur = strdup(structure); */
/*   for(i=0; i<length;i++){ */
/*     if(structur[i] == '<') structur[i]='('; */
/*     if(structur[i] == '>') structur[i]=')'; */
/*   } */
//  structur[length]='\0';   
//  printf("%s \n", structur);
   pair_table=alimake_pair_table(structure);
  /* Get length of longest name and count sequences in alignment*/

  for (i=maxName=N=0; names[i] != NULL; i++) {
    N++;
    tmp=strlen(names[i]);
    if (tmp>maxName)  maxName=tmp;
  }

  
  /* x-coord. where sequences start */
  seqsX=namesX+maxName*fontWidth+nameStep; 

  /* calculate number of digits of the alignment length */
  snprintf(tmpBuffer,length, "%i",length);
  maxNum=strlen(tmpBuffer);
  

  /* Calculate bounding box */
  tmpColumns=columnWidth;
  if (length<columnWidth){
	tmpColumns=length;
  }
  imageWidth=ceil(namesX+(maxName+tmpColumns+maxNum)*fontWidth+2*nameStep+fontWidth+numberStep);
  imageHeight=startY+ceil((float)length/columnWidth)*((N+2)*lineStep+blockStep+consStep+ssStep+rulerStep);

  /* Write postscript header including correct bounding box */
  fprintf(outfile,alnPlotHeader,0,0,(int)imageWidth,(int)imageHeight,(int)imageHeight);

  /* Create ruler and secondary structure lines */
  i=0;
  /* Init all with dots */
  for (i=0;i<(length);i++){
	ruler[i]='.';
  }
  i=0;
  for (i=0;i<length;i++){
	/* Write number every 10th position, leave out block breaks */
	if ((i+1)%10==0 && (i+1)%columnWidth!=0){
	  snprintf(tmpBuffer,length,"%i",i+1);
	  strncpy(ruler+i,tmpBuffer,strlen(tmpBuffer));
	}
  }
  ruler[length]='\0';
  
  /* Draw color annotation first */
  /* Repeat for all pairs */
  for (i=1; i<=length; i++) {
    if ((j=pair_table[i])>i) {
      /* Repeat for open and closing position */
      for (k=0;k<2;k++){
	int pairings, nonpair, s, col;
	int ptype[8] = {0,0,0,0,0,0,0,0};
	char *color;
	col = (k==0)?i-1:j-1;
	block=ceil((float)(col+1)/columnWidth);
	xx=seqsX+(col-(block-1)*columnWidth)*fontWidth;
	/* Repeat for each sequence */
	for (s=pairings=nonpair=0; s<N; s++) {
	  ptype[BP_pair[ENCODE(seqs[s][i-1])][ENCODE(seqs[s][j-1])]]++;
	}
	for (pairings=0,s=1; s<=7; s++) {
	  if (ptype[s]) pairings++;
	}
	nonpair=ptype[0];
	if (nonpair <=2) {
	  color = colorMatrix[pairings-1][nonpair];
	  for (s=0; s<N; s++) {
	    yy=startY+(block-1)*(lineStep*(N+2)+blockStep+consStep+rulerStep)+ssStep*(block)+(s+1)*lineStep;
	    
	    /* Color according due color information in pi-array, only if base pair is possible */
	    if (BP_pair[ENCODE(seqs[s][i-1])][ENCODE(seqs[s][j-1])]) {

	      fprintf(outfile, "%.1f %.1f %.1f %.1f %s box\n",
		      xx,yy-1,xx+fontWidth,yy+fontHeight+1,color);
	    }
	  }
	}
      }
    }
  }
  free(pair_table);

  /* Process rest of the output in blocks of columnWidth */

  currY=startY;
  currPos=0;

  cons =  consensus(seqs);
  
  while (currPos<length) {

    /* Display secondary structure line */
    fprintf(outfile,"0 setgray\n");
    strncpy(tmpBuffer,structure+currPos,columnWidth);
    tmpBuffer[columnWidth]='\0';
    
    x=0;y=0;
    while ((c=tmpBuffer[x])){
      if (c=='.'){
	ssEscaped[y++]='.';
      } else {
	ssEscaped[y++]='\\';
	ssEscaped[y++]=c;
      }			 
      x++;
    }
    ssEscaped[y]='\0';
    
    fprintf(outfile, "(%s) %.1f %.1f string\n", ssEscaped,seqsX,currY);
    currY+=ssStep+lineStep;
    
    /* Display names, sequences and numbers */

    for (i=0; i<N; i++) {
      
      strncpy(tmpBuffer,seqs[i]+currPos,columnWidth);
      tmpBuffer[columnWidth]='\0';
      
      match=0;
      for (j=0;j<(currPos+strlen(tmpBuffer));j++){
	if (seqs[i][j] != '-') match++;
      }
      
      fprintf(outfile, "(%s) %.1f %.1f string\n", names[i],namesX,currY);
      fprintf(outfile, "(%s) %.1f %.1f string\n", tmpBuffer,seqsX,currY);
      fprintf(outfile, "(%i) %.1f %.1f string\n", match,seqsX+fontWidth*(strlen(tmpBuffer))+numberStep,currY);
      currY+=lineStep;
    }
    currY+=rulerStep;
    strncpy(tmpBuffer,ruler+currPos,columnWidth);
    tmpBuffer[columnWidth]='\0';
    fprintf(outfile, "(%s) %.1f %.1f string\n", tmpBuffer,seqsX,currY);
    
    currY+=lineStep;
    currY+=consStep;
    
    /*Display conservation bar*/
    
    fprintf(outfile,"0.6 setgray\n");
    for (i=currPos;(i<currPos+columnWidth && i<length);i++){
      match=0;
      for (j=0;j<N;j++){
	if (cons[i] == seqs[j][i]) match++;
	if (cons[i]=='U' && seqs[j][i]=='T') match++;
	if (cons[i]=='T' && seqs[j][i]=='U') match++;
      }
      score=(float)(match-1)/(N-1);
      
      if (cons[i] == '-' ||
	  cons[i] == '_' ||
	  cons[i] == '.'){
	score=0;
      }
      
      barHeight=maxConsBar*score;
      if (barHeight==0){
	barHeight=1;
      }
      
      xx=seqsX+(i-(columnWidth*currPos/columnWidth))*fontWidth;
      
      fprintf(outfile,"%.1f %.1f %.1f %.1f box2\n",
	      xx,
	      currY+maxConsBar-barHeight,
	      xx+fontWidth,
	      currY+maxConsBar);
    }
    
    currY+=blockStep;
    currPos+=columnWidth;
  }
  free(cons);
  fprintf(outfile,"showpage\n");
  fclose(outfile);
  free(tmpBuffer);
  free(ssEscaped);free(ruler);
  
  return 0;

}

/*###########################################*/
/*# deprecated functions below              #*/
/*###########################################*/

int PS_dot_plot(char *string, char *wastlfile) {
  /* this is just a wrapper to call PS_dot_plot_list */
  int i, j, k, length, maxl, mf_num;
  struct plist *pl;
  struct plist *mf;

  length = strlen(string);
  maxl = 2*length;
  pl = (struct plist *)space(maxl*sizeof(struct plist));
  k=0;
  /*make plist out of pr array*/
  for (i=1; i<length; i++)
    for (j=i+1; j<=length; j++) {
      if (pr[iindx[i]-j]<PMIN) continue;
      if (k>=maxl-1) {
        maxl *= 2;
        pl = (struct plist *)xrealloc(pl,maxl*sizeof(struct plist));
      }
      pl[k].i = i;
      pl[k].j = j;
      pl[k++].p = pr[iindx[i]-j];
    }
  pl[k].i=0;
  pl[k].j=0;
  pl[k++].p=0.;
  /*make plist out of base_pair array*/
  mf_num = base_pair ? base_pair[0].i : 0;
  mf = (struct plist *)space((mf_num+1)*sizeof(struct plist));
  for (k=0; k<mf_num; k++) {
    mf[k].i = base_pair[k+1].i;
    mf[k].j = base_pair[k+1].j;
    mf[k].p = 0.95*0.95;
  }
  mf[k].i=0;
  mf[k].j=0;
  mf[k].p=0.;
  i = PS_dot_plot_list(string, wastlfile, pl, mf, "");
  free(mf);
  free(pl);
  return (i);
}

