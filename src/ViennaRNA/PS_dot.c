/*
        PostScript and GML output for RNA secondary structures
                    and pair probability matrices

                 c  Ivo Hofacker and Peter F Stadler
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
#include "ViennaRNA/model.h"
#include "ViennaRNA/utils.h"
#include "ViennaRNA/fold_vars.h"
#include "ViennaRNA/PS_dot.h"
#include "ViennaRNA/aln_util.h"
#include "ViennaRNA/gquad.h"

/*
#################################
# PRIVATE MACROS                #
#################################
*/

#define SIZE 452.
#define PMIN 0.00001

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

static const char *RNAdp_prolog =
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
"     len lineto\n"
"     dup\n"
"     len exch sub 0 exch moveto\n"
"     len exch len exch sub lineto\n"
"     stroke\n"
"  } for\n"
"  [] 0 setdash\n"
"  0.04 setlinewidth\n"
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

static const char *RNAdp_gquad_triangle =
"/min { 2 copy gt { exch } if pop } bind def\n\n"
"/utri{ % i j prob utri\n"
"  gsave\n"
"  1 min 2 div\n"
"  0.85 mul 0.15 add 0.95  0.33\n"
"  3 1 roll % prepare hsb color\n"
"  sethsbcolor\n"
"  % now produce the coordinates for lines\n"
"  exch 1 sub dup len exch sub dup 4 -1 roll dup 3 1 roll dup len exch sub\n"
"  moveto lineto lineto closepath fill\n"
"  grestore\n"
"} bind def\n"
"/uUDmotif{ % i j uUDmotif\n"
"  gsave\n"
"  1 min 2 div\n"
"  0.85 mul 0.15 add 0.95 0.6\n"
"  3 1 roll % prepare hsb color\n"
"  sethsbcolor\n"
"  % now produce the coordinates for lines\n"
"  exch 1 sub dup len exch sub dup 4 -1 roll dup 3 1 roll dup len exch sub\n"
"  moveto lineto lineto closepath fill\n"
"  grestore\n"
"} bind def\n"
"/lUDmotif{ % i j lUDmotif\n"
"  gsave\n"
"  1 min 2 div\n"
"  0.85 mul 0.15 add 0.95 0.6\n"
"  3 1 roll % prepare hsb color\n"
"  sethsbcolor\n"
"  % now produce the coordinates for lines\n"
"  dup len exch sub dup 4 -1 roll 1 sub dup 3 1 roll dup len exch sub\n"
"  moveto lineto lineto closepath fill\n"
"  grestore\n"
"} bind def\n"
"/uHmotif{ % i j uHmotif\n"
"  gsave\n"
"  1 min 2 div\n"
"  0.85 mul 0.15 add 0.95  0.99\n"
"  3 1 roll % prepare hsb color\n"
"  sethsbcolor\n"
"  % now produce the coordinates for lines\n"
"  exch 1 sub dup len exch sub dup 4 -1 roll dup 3 1 roll dup len exch sub\n"
"  moveto lineto lineto closepath fill\n"
"  grestore\n"
"} bind def\n"
"/lHmotif{ % i j lHmotif\n"
"  gsave\n"
"  1 min 2 div\n"
"  0.85 mul 0.15 add 0.95  0.99\n"
"  3 1 roll % prepare hsb color\n"
"  sethsbcolor\n"
"  % now produce the coordinates for lines\n"
"  dup len exch sub dup 4 -1 roll 1 sub dup 3 1 roll dup len exch sub\n"
"  moveto lineto lineto closepath fill\n"
"  grestore\n"
"} bind def\n"
"/uImotif{ % i j k l uImotif\n"
"  gsave\n"
"  1 min 2 div\n"
"  0.85 mul 0.15 add 0.95  0.99\n"
"  3 1 roll % prepare hsb color\n"
"  sethsbcolor\n"
"  % now produce the coordinates for lines\n"
"  1 sub dup 5 1 roll exch len exch sub dup 5 1 roll 3 -1 roll dup\n"
"  5 1 roll exch 4 1 roll 3 1 roll exch 1 sub len exch sub dup 3 1 roll\n"
"  moveto lineto lineto lineto closepath fill\n"
"  grestore\n"
"} bind def\n"
"/lImotif{ % i j k l lImotif\n"
"  gsave\n"
"  1 min 2 div\n"
"  0.85 mul 0.15 add 0.95  0.99\n"
"  3 1 roll % prepare hsb color\n"
"  sethsbcolor\n"
"  % now produce the coordinates for lines\n"
"  4 -1 roll 1 sub dup 5 1 roll exch 1 sub len exch sub dup 3 -1 roll exch\n"
"  5 -1 roll len exch sub dup 6 -1 roll dup 3 1 roll 7 4 roll\n"
"  moveto lineto lineto lineto closepath fill\n"
"  grestore\n"
"} bind def\n";

static const char *RNAdp_prolog_turn =
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


/*
#################################
# PRIVATE FUNCTION DECLARATIONS #
#################################
*/

PRIVATE FILE  *PS_dot_common(char *seq, int cp, char *wastlfile, char *comment, int winsize);
PRIVATE int   sort_plist_by_type_desc(const void *p1, const void *p2);
PRIVATE int   sort_plist_by_prob_asc(const void *p1, const void *p2);
/*
#################################
# BEGIN OF FUNCTION DEFINITIONS #
#################################
*/


int PS_color_dot_plot(char *seq, cpair *pi, char *wastlfile) {
  /* produce color PostScript dot plot from cpair */

  FILE *wastl;
  int i;

  wastl = PS_dot_common(seq, cut_point, wastlfile, NULL, 0);
  if (wastl==NULL)  return 0; /* return 0 for failure */

  fprintf(wastl, "/hsb {\n"
          "dup 0.3 mul 1 exch sub sethsbcolor\n"
          "} bind def\n\n");

  fprintf(wastl, "\n%%draw the grid\ndrawgrid\n\n");
  fprintf(wastl,"%%start of base pair probability data\n");

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


static int sort_plist_by_type_desc(const void *p1, const void *p2){
  if(((plist*)p1)->type > ((plist*)p2)->type) return -1;
  if(((plist*)p1)->type < ((plist*)p2)->type) return 1;
  return 0;
}

static int sort_plist_by_prob_asc(const void *p1, const void *p2){
  if(((plist*)p1)->p > ((plist*)p2)->p) return 1;
  if(((plist*)p1)->p < ((plist*)p2)->p) return -1;
  return 0;
}

PUBLIC int
PS_dot_plot_list( char *seq,
                  char *wastlfile,
                  plist *pl,
                  plist *mf,
                  char *comment){

  return vrna_plot_dp_PS_list(seq, cut_point, wastlfile, pl, mf, comment);
}

PUBLIC int
vrna_plot_dp_PS_list( char *seq,
                      int cp,
                      char *wastlfile,
                      plist *pl,
                      plist *mf,
                      char *comment){

  FILE *wastl;
  int pl_size, gq_num;
  double tmp;
  plist *pl1;

  wastl = PS_dot_common(seq, cp, wastlfile, comment, 0);
  if (wastl==NULL) return 0; /* return 0 for failure */

  fprintf(wastl, "%s\n", RNAdp_gquad_triangle);

  fprintf(wastl,"%%data starts here\n");

  /* sort the plist to bring all gquad triangles to the front */
  for(gq_num = pl_size = 0, pl1 = pl; pl1->i > 0; pl1++, pl_size++)
    if(pl1->type == VRNA_PLIST_TYPE_GQUAD) gq_num++;
  qsort(pl, pl_size, sizeof(plist), sort_plist_by_type_desc);
  /* sort all gquad triangles by probability to bring lower probs to the front */
  qsort(pl, gq_num, sizeof(plist), sort_plist_by_prob_asc);

  /* print triangles for g-quadruplexes in upper half */
  fprintf(wastl,"\n%%start of quadruplex data\n");
  for (pl1=pl; pl1->i > 0; pl1++) {
    if(pl1->type == VRNA_PLIST_TYPE_GQUAD){
      tmp = sqrt(pl1->p);
      fprintf(wastl, "%d %d %1.9f utri\n", pl1->i, pl1->j, tmp);
    }
  }

  /* print triangles for hairpin loop motifs in upper half */
  fprintf(wastl,"\n%%start of Hmotif data\n");
  for (pl1=pl; pl1->i > 0; pl1++) {
    if(pl1->type == VRNA_PLIST_TYPE_H_MOTIF){
      tmp = sqrt(pl1->p);
      fprintf(wastl, "%d %d %1.9f uHmotif\n", pl1->i, pl1->j, tmp);
    }
  }
  for (pl1=mf; pl1->i > 0; pl1++) {
    if(pl1->type == VRNA_PLIST_TYPE_H_MOTIF){
      tmp = sqrt(pl1->p);
      fprintf(wastl, "%d %d %1.9f lHmotif\n", pl1->i, pl1->j, tmp);
    }
  }

  /* print triangles for interior loop motifs in upper half */
  fprintf(wastl,"\n%%start of Imotif data\n");
  int   a,b;
  float ppp;
  a = b = 0;
  for (pl1=pl; pl1->i > 0; pl1++) {
    if(pl1->type == VRNA_PLIST_TYPE_I_MOTIF){
      if(a == 0){
        a = pl1->i;
        b = pl1->j;
        ppp = tmp = sqrt(pl1->p);
      } else {
        fprintf(wastl, "%d %d %d %d %1.9f uImotif\n", a, b, pl1->i, pl1->j, ppp);
        a = b = 0;
      }
    }
  }
  for (a = b= 0, pl1=mf; pl1->i > 0; pl1++) {
    if(pl1->type == VRNA_PLIST_TYPE_I_MOTIF){
      if(a == 0){
        a = pl1->i;
        b = pl1->j;
        ppp = sqrt(pl1->p);
      } else {
        fprintf(wastl, "%d %d %d %d %1.9f lImotif\n", a, b, pl1->i, pl1->j, ppp);
        a = b = 0;
      }
    }
  }

  /* print triangles for unstructured domain motifs in upper half */
  fprintf(wastl,"\n%%start of unstructured domain motif data\n");
  for (pl1=pl; pl1->i > 0; pl1++) {
    if(pl1->type == VRNA_PLIST_TYPE_UD_MOTIF){
      tmp = sqrt(pl1->p);
      fprintf(wastl, "%d %d %1.9f uUDmotif\n", pl1->i, pl1->j, tmp);
    }
  }
  for (pl1=mf; pl1->i > 0; pl1++) {
    if(pl1->type == VRNA_PLIST_TYPE_UD_MOTIF){
      tmp = sqrt(pl1->p);
      fprintf(wastl, "%d %d %1.9f lUDmotif\n", pl1->i, pl1->j, tmp);
    }
  }

  fprintf(wastl, "\n%%draw the grid\ndrawgrid\n\n");
  fprintf(wastl,"%%start of base pair probability data\n");

  /* print boxes in upper right half*/
  for (pl1 = pl; pl1->i>0; pl1++) {
    tmp = sqrt(pl1->p);
    if(pl1->type == VRNA_PLIST_TYPE_BASEPAIR)
        fprintf(wastl,"%d %d %1.9f ubox\n", pl1->i, pl1->j, tmp);
  }


  /* print boxes in lower left half (mfe) */
  for (pl1=mf; pl1->i>0; pl1++) {
    tmp = sqrt(pl1->p);
    if(pl1->type == VRNA_PLIST_TYPE_BASEPAIR)
      fprintf(wastl,"%d %d %1.7f lbox\n", pl1->i, pl1->j, tmp);
  }

  fprintf(wastl,"showpage\n"
          "end\n"
          "%%%%EOF\n");
  fclose(wastl);
  return 1; /* success */
}

int PS_color_dot_plot_turn(char *seq, cpair *pi, char *wastlfile, int winSize) {
  /* produce color PostScript dot plot from cpair */

  FILE *wastl;
  int i;

  wastl = PS_dot_common(seq, cut_point, wastlfile, NULL, winSize);
  if (wastl==NULL)
    return 0; /* return 0 for failure */

  fprintf(wastl, "/hsb {\n"
          "dup 0.3 mul 1 exch sub sethsbcolor\n"
          "} bind def\n\n"
          "%%BEGIN DATA\n");

  if(winSize > 0)
    fprintf(wastl, "\n%%draw the grid\ndrawgrid_turn\n\n");
  else
    fprintf(wastl, "\n%%draw the grid\ndrawgrid\n\n");
  fprintf(wastl,"%%start of base pair probability data\n");

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

int PS_dot_plot_turn(char *seq, plist *pl, char *wastlfile, int winSize) {
  /* produce color PostScript dot plot from cpair */

  FILE *wastl;
  int i;

  wastl = PS_dot_common(seq, cut_point, wastlfile, NULL, winSize);
  if (wastl==NULL)
    return 0; /* return 0 for failure */

  if(winSize > 0)
    fprintf(wastl, "\n%%draw the grid\ndrawgrid_turn\n\n");
  else
    fprintf(wastl, "\n%%draw the grid\ndrawgrid\n\n");
  fprintf(wastl,"%%start of base pair probability data\n");
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

static FILE *
PS_dot_common(char *seq,
              int cp,
              char *wastlfile,
              char *comment,
              int winsize){

  /* write PS header etc for all dot plot variants */
  FILE *wastl;
  char *name, *c;
  int i;

  wastl = fopen(wastlfile,"w");
  if (wastl==NULL) {
    vrna_message_warning("can't open %s for dot plot", wastlfile);
    return NULL; /* return 0 for failure */
  }
  name = strdup(wastlfile);
  if((c=strrchr(name, '_')))
    *c='\0';

  fprintf(wastl,
          "%%!PS-Adobe-3.0 EPSF-3.0\n"
          "%%%%Title: RNA Dot Plot\n"
          "%%%%Creator: ViennaRNA-%s\n"
          "%%%%CreationDate: %s", VERSION, vrna_time_stamp());
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
  if (cp>0) fprintf(wastl,"/cutpoint %d def\n\n", cp);


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
          "45 rotate\n\n");
  }
  else
    fprintf(wastl,"drawseq\n"
            "0.5 dup translate\n"
            "%% draw diagonal\n"
            "0.04 setlinewidth\n"
            "0 len moveto len 0 lineto stroke\n\n");

  free(name);
  return(wastl);
}

#ifdef VRNA_BACKWARD_COMPAT

/*###########################################*/
/*# deprecated functions below              #*/
/*###########################################*/

int PS_dot_plot(char *string, char *wastlfile) {
  /* this is just a wrapper to call PS_dot_plot_list */
  int i, j, k, length, maxl, mf_num;
  plist *pl;
  plist *mf;

  length = strlen(string);
  maxl = 2*length;
  pl = (plist *)vrna_alloc(maxl*sizeof(plist));
  k=0;
  /*make plist out of pr array*/
  for (i=1; i<length; i++)
    for (j=i+1; j<=length; j++) {
      if (pr[iindx[i]-j]<PMIN) continue;
      if (k>=maxl-1) {
        maxl *= 2;
        pl = (plist *)vrna_realloc(pl,maxl*sizeof(plist));
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
  mf = (plist *)vrna_alloc((mf_num+1)*sizeof(plist));
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

#endif

