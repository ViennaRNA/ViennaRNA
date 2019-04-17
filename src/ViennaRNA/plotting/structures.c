/*
        PostScript and other output formats for RNA secondary structure plots

                 c  Ivo Hofacker, Peter F Stadler, Ronny Lorenz
                          Vienna RNA package
*/

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <ctype.h>
#include "ViennaRNA/utils/basic.h"
#include "ViennaRNA/utils/alignments.h"
#include "ViennaRNA/alphabet.h"
#include "ViennaRNA/model.h"
#include "ViennaRNA/fold_vars.h"
#include "ViennaRNA/gquad.h"
#include "ViennaRNA/plotting/layouts.h"
#include "ViennaRNA/plotting/utils.h"
#include "ViennaRNA/plotting/structures.h"
#include "ViennaRNA/plotting/RNApuzzler/RNApuzzler.h"
#include "ViennaRNA/plotting/RNApuzzler/RNAturtle.h"

#include "ViennaRNA/static/templates_postscript.h"

#include "ViennaRNA/plotting/ps_helpers.inc"

/*
#################################
# PRIVATE MACROS                #
#################################
*/
#define SIZE 452.

/*
#################################
# GLOBAL VARIABLES              #
#################################
*/

/*
#################################
# PRIVATE VARIABLES             #
#################################
*/

/*
#################################
# PRIVATE FUNCTION DECLARATIONS #
#################################
*/
PRIVATE int
rnaplot_EPS(const char         *seq,
            const char         *structure,
            const char         *ssfile,
            const char         *pre,
            const char         *post,
            vrna_md_t          *md_p,
            vrna_plot_layout_t *layout);

/*
#################################
# BEGIN OF FUNCTION DEFINITIONS #
#################################
*/

PUBLIC int
vrna_file_PS_rnaplot( const char *string,
                      const char *structure,
                      const char *ssfile,
                      vrna_md_t  *md_p){

  return vrna_file_PS_rnaplot_a(string, structure, ssfile, NULL, NULL, md_p);
}


PUBLIC int
vrna_file_PS_rnaplot_a( const char *seq,
                        const char *structure,
                        const char *ssfile,
                        const char *pre,
                        const char *post,
                        vrna_md_t  *md_p)
{
  int               ret;
  vrna_plot_layout_t *layout;

  layout = vrna_plot_layout(structure, rna_plot_type);
  ret    = vrna_file_PS_rnaplot_layout(seq,
                                       structure,
                                       ssfile,
                                       pre,
                                       post,
                                       md_p,
                                       layout);

  vrna_plot_layout_free(layout);

  return ret;
}


PUBLIC int
vrna_file_PS_rnaplot_layout(const char         *seq,
                            const char         *structure,
                            const char         *ssfile,
                            const char         *pre,
                            const char         *post,
                            vrna_md_t          *md_p,
                            vrna_plot_layout_t *layout)
{
  if (!ssfile) {
    vrna_message_warning("vrna_file_PS_rnaplot*(): "
                         "Filename missing!");
  } else if (!seq) {
    vrna_message_warning("vrna_file_PS_rnaplot*(): "
                         "Sequence missing");
  } else if (!structure) {
    vrna_message_warning("vrna_file_PS_rnaplot*(): "
                         "Structure missing");
  } else if (!layout) {
    vrna_message_warning("vrna_file_PS_rnaplot*(): "
                         "Layout missing");
  } else if ((strlen(seq) != strlen(structure)) ||
             (strlen(structure) != layout->length)) {
      vrna_message_warning("vrna_file_PS_rnaplot*(): "
                           "Sequence, structure, and coordinate layout have different lengths! "
                           "(%u vs. %u vs. %u)",
                           strlen(seq),
                           strlen(structure),
                           layout->length);
  } else {
    return rnaplot_EPS(seq,
                       structure,
                       ssfile,
                       pre,
                       post,
                       md_p,
                       layout);
  }

  return 0;
}


/*
 #####################################
 # BEGIN OF STATIC HELPER FUNCTIONS  #
 #####################################
 */

PRIVATE int
rnaplot_EPS(const char         *seq,
            const char         *structure,
            const char         *ssfile,
            const char         *pre,
            const char         *post,
            vrna_md_t          *md_p,
            vrna_plot_layout_t *layout){

  float  xmin, xmax, ymin, ymax;
  int    i, length;
  int    ee, gb, ge, Lg, l[3], bbox[4];
  float *X, *Y;
  FILE  *xyplot;
  short *pair_table;
  char  *c, *string;
  vrna_md_t   md;

  if(!md_p){
    set_model_details(&md);
    md_p  = &md;
  }
    
  string = strdup(seq);
  length = strlen(string);

  xyplot = fopen(ssfile, "w");
  if (xyplot == NULL) {
    vrna_message_warning("can't open file %s - not doing xy_plot", ssfile);
    return 0;
  }

  pair_table = vrna_ptable(structure);

  bbox[0] = bbox[1] = 0;
  bbox[2] = bbox[3] = 700;

  print_PS_header(xyplot,
                  "RNA Secondary Structure Plot",
                  bbox,
                  md_p,
                  "To switch off outline pairs of sequence comment or\n"
                  "delete the appropriate line near the end of the file",
                  "RNAplot",
                  PS_MACRO_LAYOUT_BASE | ((pre || post) ? PS_MACRO_LAYOUT_EXTRAS : 0));

  fprintf(xyplot, "%% data start here\n");

  /* cut_point */
  if ((c = strchr(structure, '&'))) {
    int cutpoint;
    cutpoint = c - structure;
    string[cutpoint] = ' '; /* replace & with space */
    fprintf(xyplot, "/cutpoint %d def\n", cutpoint);
  }

  /* sequence */
  print_PS_sequence(xyplot, string);

  /* coordinates */
  print_PS_coords(xyplot, layout->x, layout->y, length);

  fprintf(xyplot, "/arcs [\n");
  if (layout->arcs) {
    for (i = 0; i < length; i++)
    {
      if (layout->arcs[6*i + 2] > 0) { /* 6*i+2 is the radius parameter */
        fprintf(xyplot, "[%3.8f %3.8f %3.8f %3.8f %3.8f %3.8f]\n",
          layout->arcs[6*i + 0],
          layout->arcs[6*i + 1],
          layout->arcs[6*i + 2],
          layout->arcs[6*i + 3],
          layout->arcs[6*i + 4],
          layout->arcs[6*i + 5]
        );
      } else {
        fprintf(xyplot, "[]\n");
      }
    }

  } else {
    for (i = 0; i < length; i++)
      fprintf(xyplot, "[]\n");
  }
  fprintf(xyplot, "] def\n");

  /* correction coordinates for quadratic beziers in case we produce a circplot */
  if(rna_plot_type == VRNA_PLOT_TYPE_CIRCULAR)
    fprintf(xyplot, "/cpr %6.2f def\n", (float)3*length);
  /* base pairs */
  fprintf(xyplot, "/pairs [\n");
  for (i = 1; i <= length; i++)
    if (pair_table[i]>i)
      fprintf(xyplot, "[%d %d]\n", i, pair_table[i]);
  /* add gquad pairs */
  ge=0;
  while ( (ee=parse_gquad(structure+ge, &Lg, l)) >0 ) {
    int k;
    fprintf(xyplot, "%% gquad\n");
    ge += ee;
    gb=ge-Lg*4-l[0]-l[1]-l[2]+1; /* add pseudo-base pair encloding gquad */
    for (k=0; k<Lg; k++) {
      int ii, jj, il;
      for (il=0, ii=gb+k; il<3; il++) {
        jj = ii+l[il]+Lg;
        fprintf(xyplot, "[%d %d]\n", ii, jj);
        ii = jj;
      }
      jj = gb+k;
      fprintf(xyplot, "[%d %d]\n", jj, ii);
    }
  }

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

  print_PS_footer(xyplot);

  fclose(xyplot);

  free(string);
  free(pair_table);

  return 1; /* success */
}

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
  short *pair_table;
  float *X, *Y;

  gmlfile = fopen(ssfile, "w");
  if (gmlfile == NULL) {
     vrna_message_warning("can't open file %s - not doing xy_plot", ssfile);
     return 0;
  }

  length = strlen(string);

  pair_table = vrna_ptable(structure);

  switch(option){
  case 'X' :
  case 'x' :
    /* Simple XY Plot */
    if (rna_plot_type == 0)
      i = vrna_plot_coords_simple_pt(pair_table, &X, &Y);
    else
      i = vrna_plot_coords_naview_pt(pair_table, &X, &Y);

    if(i!=length)
      vrna_message_warning("strange things happening in gmlRNA ...");
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
          "# Options: %s\n", VERSION, vrna_time_stamp(), ssfile, option_string());
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



PUBLIC int
PS_rna_plot_snoop_a(const char *string,
                    const char *structure,
                    const char *ssfile,
                    int *relative_access,
                    const char *seqs[])
{
  int    i, length, bbox[4];
  float *X, *Y;
  FILE  *xyplot;
  short *pair_table;
  short *pair_table_snoop;
  vrna_md_t md;

  set_model_details(&md);

  length = strlen(string);

  xyplot = fopen(ssfile, "w");
  if (xyplot == NULL) {
    vrna_message_warning("can't open file %s - not doing xy_plot", ssfile);
    return 0;
  }

  pair_table = vrna_ptable(structure);
  pair_table_snoop = vrna_pt_snoop_get(structure);

  if (rna_plot_type == 0)
    i = vrna_plot_coords_simple_pt(pair_table, &X, &Y);
  else
    i = vrna_plot_coords_naview_pt(pair_table, &X, &Y);

  if(i!=length)
    vrna_message_warning("strange things happening in PS_rna_plot...");
/*   printf("cut_point %d\n", cut_point); */

/*   for (i = 1; i < length; i++) { */
/*     printf("%d X %f Y %f \n", i, X[i], Y[i]); */
/*     xmin = X[i] < xmin ? X[i] : xmin; */
/*     xmax = X[i] > xmax ? X[i] : xmax; */
/*     ymin = Y[i] < ymin ? Y[i] : ymin; */
/*     ymax = Y[i] > ymax ? Y[i] : ymax; */
/*   } */
  /* localize centre of the interaction bucket. Geometry */
  
  for (i = 1; i < cut_point; i++) {  /* interior loop of size 0 */
    if(pair_table_snoop[i] != 0){ 
      X[i-1]=X[pair_table_snoop[i]-1]; 
      Y[i-1]=Y[pair_table_snoop[i]-1]; 
    }
    else if(pair_table_snoop[i-1] && pair_table_snoop[i+1]){ /* interior loop of size 1 */
      X[i-1]=X[pair_table_snoop[i-1] -1-1];
      Y[i-1]=Y[pair_table_snoop[i-1] -1-1];
    } 
    else if(pair_table_snoop[i-1] && pair_table_snoop[i+2]){ /* interior loop of size 2 */
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
    else if(pair_table_snoop[i-1] && pair_table_snoop[i+3]){ /* interior loop of size 2 */
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
  float X0=-1,Y0=-1,X1=-1,Y1=-1,X2=-1,Y2=-1;
/*   int c1,c2,c3; */
  for(i=1;i<cut_point; i++){
    if(pair_table_snoop[i]){
      X0=X[pair_table_snoop[i]-1];Y0=Y[pair_table_snoop[i]-1];
  /*     c1=pair_table_snoop[i]; */
      i++;
      break;
    }
  }
  for(;i<cut_point; i++){
    if(pair_table_snoop[i]){
      X1=X[pair_table_snoop[i]-1];Y1=Y[pair_table_snoop[i]-1];
    /*   c2=pair_table_snoop[i]; */
      i++;
      break;
    }
  }
  for(;i<cut_point; i++){
    if(pair_table_snoop[i]){
      X2=X[pair_table_snoop[i]-1];Y2=Y[pair_table_snoop[i]-1];
    /*   c3=pair_table_snoop[i]; */
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
  /*    if(abs(alpha -alpha_p) > 0.0000001){ */
  xC  =  (b_p - b) / (alpha - alpha_p);
  yC  =  alpha * xC + b;
  for (i = 1; i < cut_point; i++) {  
     X[i-1] = X[i-1] + 0.25*(xC-X[i-1]);  
     Y[i-1] = Y[i-1] + 0.25*(yC-Y[i-1]);  
  }  

  bbox[0] = bbox[1] = 0;
  bbox[2] = bbox[3] = 700;

  print_PS_header(xyplot,
                  "RNA Secondary Structure Plot",
                  bbox,
                  &md,
                  "To switch off outline pairs of sequence comment or\n"
                  "delete the appropriate line near the end of the file",
                  "RNAplot",
                  PS_MACRO_LAYOUT_BASE | PS_MACRO_LAYOUT_EXTRAS);

  char **A;
  if(seqs){
    A = vrna_annotate_covar_db_extended((const char **)seqs,
                                        structure,
                                        &md,
                                        VRNA_BRACKETS_RND | VRNA_BRACKETS_ANG | VRNA_BRACKETS_SQR);
  }
  
  fprintf(xyplot, "%% data start here\n");
  /* cut_point */
  if (cut_point > 0 && cut_point <= strlen(string))
    fprintf(xyplot, "/cutpoint %d def\n", cut_point-1);

  /* sequence */
  print_PS_sequence(xyplot, string);

  /* coordinates */
  print_PS_coords(xyplot, X, Y, length);

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
                   "  0        \n"                              
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
  /* fprintf(xyplot, "%d cmark\n",c1); */
  /* fprintf(xyplot, "%d cmark\n",c2); */
  /* fprintf(xyplot, "%d cmark\n",c3); */
  if (seqs) { 
     fprintf(xyplot, "%% Start Annotations\n"); 
     fprintf(xyplot, "%s\n", A[1]); 
     fprintf(xyplot, "%% End Annotations\n"); 
   } 

  print_PS_footer(xyplot);

  fclose(xyplot);
  if(seqs){free(A[0]);free(A[1]);free(A);}
  free(pair_table);free(pair_table_snoop);
  free(X); free(Y);
  return 1; /* success */
}


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
    vrna_message_warning("can't open file %s - not doing xy_plot", ssfile);
    return 0;
  }

  pair_table = vrna_ptable(structure);

  switch(rna_plot_type){
    case VRNA_PLOT_TYPE_SIMPLE:   i = vrna_plot_coords_simple_pt(pair_table, &X, &Y);
                                  break;
    case VRNA_PLOT_TYPE_CIRCULAR: {
                                    int radius = 3*length;
                                    int dr = 0;
                                    R = (float *) vrna_alloc((length+1)*sizeof(float));
                                    CX = (float *) vrna_alloc((length+1)*sizeof(float));
                                    CY = (float *) vrna_alloc((length+1)*sizeof(float));
                                    i = vrna_plot_coords_circular_pt(pair_table, &X, &Y);
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
    default:                      i = vrna_plot_coords_naview_pt(pair_table, &X, &Y);
                                  break;
  }

  if(i!=length)
    vrna_message_warning("strange things happening in PS_rna_plot...");

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
     vrna_message_warning("can't open file %s - not doing xy_plot", ssfile);
     return 0;
  }
  length = strlen(string);
  pair_table = vrna_ptable(structure);

  /* make coordinates */
  if (rna_plot_type == 0)
    i = vrna_plot_coords_simple_pt(pair_table, &X, &Y);
  else
    i = vrna_plot_coords_naview_pt(pair_table, &X, &Y);

  if (i!=length)
    vrna_message_warning("strange things happening in ssv_rna_plot...");

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
          "# Options: %s\n", VERSION, vrna_time_stamp(), ssfile, option_string());
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
    vrna_message_warning("can't open file %s - not doing xy_plot", ssfile);
    return 0;
  }

  length = strlen(string);
  pair_table = vrna_ptable(structure);

  /* make coordinates */
  if (rna_plot_type == 0)
    i = vrna_plot_coords_simple_pt(pair_table, &X, &Y);
  else
    i = vrna_plot_coords_naview_pt(pair_table, &X, &Y);

  if (i!=length)
    vrna_message_warning("strange things happening in xrna_plot...");

  fprintf(ss_file,
          "# Vienna RNA Package %s, XRNA output\n"
          "# CreationDate: %s\n"
          "# Options: %s\n", VERSION, vrna_time_stamp(), option_string());
  for (i=1; i<=length; i++)
    /* XRNA likes to have coordinate mirrored, so we use (-X, Y) */
    fprintf(ss_file, "%d %c %6.2f %6.2f %d %d\n", i, string[i-1],
            -X[i-1], Y[i-1], (pair_table[i]?1:0), pair_table[i]);
  fclose(ss_file);

  free(pair_table);
  free(X); free(Y);
  return 1; /* success */
}

#ifndef VRNA_DISABLE_BACKWARD_COMPATIBILITY

/*###########################################*/
/*# deprecated functions below              #*/
/*###########################################*/

PUBLIC int
PS_rna_plot(char *string,
            char *structure,
            char *ssfile){

  return vrna_file_PS_rnaplot((const char*)string,
                              (const char*)structure,
                              (const char*) ssfile,
                              NULL);
}

PUBLIC int
PS_rna_plot_a(char *string,
              char *structure,
              char *ssfile,
              char *pre,
              char *post){

  return vrna_file_PS_rnaplot_a((const char*)string,
                                (const char*)structure,
                                (const char*)ssfile,
                                (const char*)pre,
                                (const char*)post,
                                NULL);
}

PUBLIC int
PS_rna_plot_a_gquad(char *string,
                    char *structure,
                    char *ssfile,
                    char *pre,
                    char *post){

  return vrna_file_PS_rnaplot_a((const char*)string,
                                (const char*)structure,
                                (const char*)ssfile,
                                (const char*)pre,
                                (const char*)post,
                                NULL);
}

#endif

