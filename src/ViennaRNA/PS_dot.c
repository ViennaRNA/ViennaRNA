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
"/DataVisible  [ true true true true] def\n"
"/DataTitles   [ false false false false ] def\n\n"
"/min { 2 copy gt { exch } if pop } bind def\n\n"
"/max { 2 copy lt { exch } if pop } bind def\n\n"
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
"  gsave\n"
"  0.5 dup translate"
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
"  %% draw diagonal\n"
"  0 len moveto len 0 lineto stroke\n"
"  grestore\n"
"} bind def\n\n"
"/drawTitle {\n"
"  currentdict /DPtitle known {\n"
"    % center title text\n"
"    /Helvetica findfont 10 scalefont setfont\n"
"    360 705 moveto DPtitle dup stringwidth pop 2 div neg 0 rmoveto show\n"
"  } if\n"
"} bind def\n\n"
"/prepareCoords {\n"
"  0 1 3 {\n"
"    % check whether we want to display current data\n"
"    dup DataVisible exch get\n"
"    {\n"
"      % check whether we've actually got some data\n"
"      DataSource exch get dup currentdict exch known {\n"
"        % data source s_j is present, so find length of array\n"
"        currentdict exch get length \n"
"      } { pop 0 } ifelse\n"
"    } if\n"
"  } for\n"
"  exch dup 5 -1 roll add 4 -1 roll dup 5 1 roll 4 -1 roll add max\n"
"  len add 3 add 700 exch div dup scale\n"
"  exch 1 add exch 1 add translate\n"
"} bind def\n\n";

static const char *RNAdp_sd =
"/utri{ % i j prob utri\n"
"  gsave\n"
"  0.5 dup translate\n"
"  1 min 2 div\n"
"  0.85 mul 0.15 add 0.95  0.33\n"
"  3 1 roll % prepare hsb color\n"
"  sethsbcolor\n"
"  % now produce the coordinates for lines\n"
"  exch 1 sub dup len exch sub dup 4 -1 roll dup 3 1 roll dup len exch sub\n"
"  moveto lineto lineto closepath fill\n"
"  grestore\n"
"} bind def\n\n";

static const char *RNAdp_ud =
"/uUDmotif{ % i j uUDmotif\n"
"  gsave\n"
"  0.5 dup translate\n"
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
"  0.5 dup translate\n"
"  1 min 2 div\n"
"  0.85 mul 0.15 add 0.95 0.6\n"
"  3 1 roll % prepare hsb color\n"
"  sethsbcolor\n"
"  % now produce the coordinates for lines\n"
"  dup len exch sub dup 4 -1 roll 1 sub dup 3 1 roll dup len exch sub\n"
"  moveto lineto lineto closepath fill\n"
"  grestore\n"
"} bind def\n\n";

static const char *RNAdp_sc_motifs =
"/uHmotif{ % i j uHmotif\n"
"  gsave\n"
"  0.5 dup translate\n"
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
"  0.5 dup translate\n"
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
"  0.5 dup translate\n"
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
"  0.5 dup translate\n"
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

static const char *RNAdp_linear_data =
"/drawDataSquareBottom { % x v n dataSquareBottom draw box\n"
"  len add 2 add exch lbox\n"
"} bind def\n\n"
"/drawDataSquareTop { % x v n dataSquareBottom draw box\n"
"  neg 1 sub exch lbox\n"
"} bind def\n\n"
"/drawDataSquareLeft { % y v n dataSquareBottom draw box\n"
"  neg 1 sub 3 1 roll lbox\n"
"} bind def\n\n"
"/drawDataSquareRight { % y v n dataSquareBottom draw box\n"
"  % use size x y box to draw box\n"
"  2 add len add 3 1 roll lbox\n"
"} bind def\n\n"
"/drawDataSquareBottomHSB { % x v h s b n dataSquareBottomHSB draw box\n"
"  % use size x y box to draw box\n"
"  len add 2 add 5 1 roll sethsbcolor lbox\n"
"} bind def\n\n"
"/drawDataSquareTopHSB { % x v h s b n dataSquareBottomHSB draw box\n"
"  % use size x y box to draw box\n"
"  neg 1 sub 5 1 roll sethsbcolor lbox\n"
"} bind def\n\n"
"/drawDataSquareLeftHSB { % x v h s b n dataSquareLeftHSB draw box\n"
"  % use size x y box to draw box\n"
"  neg 1 sub 6 1 roll sethsbcolor lbox\n"
"} bind def\n\n"
"/drawDataSquareRightHSB { % x v h s b n dataSquareLeftHSB draw box\n"
"  % use size x y box to draw box\n"
"  2 add len add 6 1 roll sethsbcolor lbox\n"
"} bind def\n\n"
"/drawDataTitleBottom {\n"
"  /Helvetica findfont 0.95 scalefont setfont\n"
"  0 -1.4 3 -1 roll sub moveto \n"
"  dup stringwidth pop neg 0 rmoveto   \n"
"  show\n"
"} bind def\n\n"
"/drawDataTitleTop {\n"
"  /Helvetica findfont 0.95 scalefont setfont\n"
"  0 len 1.6 add 3 -1 roll add moveto \n"
"  dup stringwidth pop neg 0 rmoveto   \n"
"  show\n"
"} bind def\n\n"
"/drawDataTitleLeft {\n"
"  /Helvetica findfont 0.95 scalefont setfont\n"
"  neg 1.4 sub len 1 add moveto \n"
"  dup stringwidth pop 0 exch rmoveto -90 rotate\n"
"  show 90 rotate\n"
"} bind def\n\n"
"/drawDataTitleRight {\n"
"  /Helvetica findfont 0.95 scalefont setfont\n"
"  1.6 add len add len 1 add moveto \n"
"  dup stringwidth pop 0 exch rmoveto -90 rotate\n"
"  show 90 rotate\n"
"} bind def\n\n"
"% do not modify the arrays below unless you know what you're doing!\n"
"/DataSource     [ /topData /leftData /bottomData /rightData ] def\n"
"/DataDrawBox    [ /drawDataSquareTop  /drawDataSquareLeft /drawDataSquareBottom /drawDataSquareRight] def\n"
"/DataDrawBoxHSB [ /drawDataSquareTopHSB /drawDataSquareLeftHSB /drawDataSquareBottomHSB /drawDataSquareRightHSB ] def\n"
"/DataDrawTitle  [ /drawDataTitleTop /drawDataTitleLeft /drawDataTitleBottom /drawDataTitleRight ] def\n\n"
"% this is the logic to parse the auxiliary linear data\n"
"% given in arrays topData, leftData, bottomData, and rightData\n"
"% See also the Boolean arrays DataVisible and DataTitles that\n"
"% are used to control which part of data will be visible\n"
"/drawData {\n"
"  0 1 3 {\n"
"    % check whether we want to display current data\n"
"    dup DataVisible exch get\n"
"    {\n"
"      % check whether we've actually got some data\n"
"      dup DataSource exch get dup currentdict exch known {\n"
"        % data source s_j is present, so we load the\n"
"        % corresponding data array a and loop over all data sets a[i]\n"
"        currentdict exch get dup length 1 sub 0 1 3 -1 roll {\n"
"          dup dup\n"
"          % now on stack: j a i i i\n\n"
"          % load data set, i.e. a[i]\n"
"          4 -1 roll         % j i i i a\n"
"          dup 3 -1 roll get dup % j i i a a[i] a[i]\n\n"
"          % 1. check whether we need to process data set title\n"
"          6 -1 roll dup 7 1 roll DataTitles exch get {\n"
"            % get current title drawing function key\n"
"            6 -1 roll dup 7 1 roll DataDrawTitle exch get\n"
"            % now on stack: ... j i i a a[i] a[i] title_draw_key\n\n"
"            % get current title and execute drawing function\n"
"            exch 0 get exch currentdict exch get 5 -1 roll exch exec\n"
"          } { % remove unused variables\n"
"              pop 3 -1 roll pop\n"
"          } ifelse\n"
"          % now on stack: ... j i a a[i]\n\n"
"          % 2. process actual data a[k] for 1 <= k < n\n"
"          dup length 1 sub 1 exch getinterval { \n"
"            % on stack: j i a a[i][k]\n"
"            gsave\n"
"            dup length 2 eq { % print black box if two-valued\n"
"              % get box drawing function\n"
"              4 -1 roll dup 5 1 roll DataDrawBox exch get currentdict exch get exch\n"
"              aload pop 5 -1 roll dup 6 1 roll 4 -1 roll exec\n"
"            } {\n"
"              dup length 5 eq { % print box with hsb color\n"
"                % get box drawing function\n"
"                4 -1 roll dup 5 1 roll DataDrawBoxHSB exch get currentdict exch get exch\n"
"                % on stack: j i a f a[i]\n"
"                % load data array and prepare for drawing\n"
"                aload pop 8 -1 roll dup 9 1 roll 7 -1 roll exec\n"
"              } { pop } ifelse\n"
"            } ifelse\n"
"            grestore\n"
"          } forall\n"
"          exch pop \n"
"          % left on stack: j a\n"
"        } for\n"
"        \n"
"      } if\n"
"    } if\n"
"  } for\n"
"} bind def\n\n";


#define dp_add_lindata(data_var, data, title_var, title, size, avail)  \
  (title_var)[size] = title; \
  (data_var)[size]  = data; \
  if((++size) == (avail)){ \
    avail *= 1.2; \
    data_var  = (vrna_data_lin_t **)vrna_realloc(data_var, sizeof(vrna_data_lin_t *) * avail); \
    title_var = (char **)vrna_realloc(title_var, sizeof(char *) * avail); \
  }


#define dp_finalize_lindata(data_var, title_var, size) \
  (data_var)[size]  = NULL; \
  (title_var)[size] = NULL; \
  (data_var)        = (vrna_data_lin_t **)vrna_realloc(data_var, sizeof(vrna_data_lin_t *) * (size + 1)); \
  (title_var)       = (char **)vrna_realloc(title_var, sizeof(char *) * (size + 1));


#define DP_MACRO_NONE         0U
#define DP_MACRO_LINEAR_DATA  1U
#define DP_MACRO_SC_MOTIFS    2U
#define DP_MACRO_SD           4U
#define DP_MACRO_UD           8U

#define DP_MACRO_ALL          (DP_MACRO_LINEAR_DATA | DP_MACRO_SC_MOTIFS | DP_MACRO_SD | DP_MACRO_UD)

/*
#################################
# PRIVATE FUNCTION DECLARATIONS #
#################################
*/

PRIVATE FILE  *PS_dot_common(const char *seq, int cp, const char *wastlfile, char *comment, int winsize, unsigned int options);
PRIVATE int   sort_plist_by_type_desc(const void *p1, const void *p2);
PRIVATE int   sort_plist_by_prob_asc(const void *p1, const void *p2);
PRIVATE void  EPS_footer(FILE *eps);
PRIVATE void  EPS_print_title(FILE *eps, const char *title);
PRIVATE void  EPS_print_seq(FILE *eps, const char *sequence);
PRIVATE void  EPS_print_header(FILE *eps, int bbox[4], const char *comment, unsigned int options);
PRIVATE void  EPS_print_ud_data(FILE *eps, plist *pl, plist *mf);
PRIVATE void  EPS_print_sd_motif_data(FILE *eps, plist *pl, plist *mf);
PRIVATE void  EPS_print_sc_motif_data(FILE *eps, plist *pl, plist *mf);
PRIVATE void  EPS_print_bpp_data(FILE *eps, plist *pl, plist *mf);
PRIVATE void  EPS_print_linear_data_top(FILE *eps, const char **data_title, vrna_data_lin_t **data);
PRIVATE void  EPS_print_linear_data_left(FILE *eps, const char **data_title, vrna_data_lin_t **data);
PRIVATE void  EPS_print_linear_data_bottom(FILE *eps, const char **data_title, vrna_data_lin_t **data);
PRIVATE void  EPS_print_linear_data_right(FILE *eps, const char **data_title, vrna_data_lin_t **data);
PRIVATE void  EPS_print_linear_data(FILE *eps, const char *varname, const char **data_title, vrna_data_lin_t **data);
PRIVATE vrna_data_lin_t *plist_to_accessibility(plist *pl, unsigned int length);
PRIVATE vrna_data_lin_t *plist_to_ud_motif_prob(plist *pl, unsigned int length);
PRIVATE void  EPS_print_sd_data(FILE *eps, vrna_plist_t *pl, vrna_plist_t *mf);

/*
#################################
# BEGIN OF FUNCTION DEFINITIONS #
#################################
*/

PUBLIC int
PS_color_dot_plot(char *seq,
                  cpair *pi,
                  char *wastlfile){

  /* produce color PostScript dot plot from cpair */

  FILE *wastl;
  int i;

  wastl = PS_dot_common(seq, cut_point, wastlfile, NULL, 0, DP_MACRO_NONE);
  if (wastl==NULL)  return 0; /* return 0 for failure */

  fprintf(wastl, "/hsb {\n"
          "dup 0.3 mul 1 exch sub sethsbcolor\n"
          "} bind def\n\n");

  fprintf(wastl,  "\n%%draw the grid\ndrawgrid\n\n");
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

   EPS_footer(wastl);

   fclose(wastl);
   return 1; /* success */
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

  wastl = PS_dot_common(seq, cp, wastlfile, comment, 0, DP_MACRO_ALL);

  if (wastl==NULL) return 0; /* return 0 for failure */

  fprintf(wastl,"%%data starts here\n");

  /* sort the plist to bring all gquad triangles to the front */
  for(gq_num = pl_size = 0, pl1 = pl; pl1->i > 0; pl1++, pl_size++)
    if (pl1->type == VRNA_PLIST_TYPE_GQUAD) gq_num++;
  qsort(pl, pl_size, sizeof(plist), sort_plist_by_type_desc);
  /* sort all gquad triangles by probability to bring lower probs to the front */
  qsort(pl, gq_num, sizeof(plist), sort_plist_by_prob_asc);

  EPS_print_sd_data(wastl, pl, mf);
  EPS_print_sc_motif_data(wastl, pl, mf);

  fprintf(wastl, "\n%%draw the grid\ndrawgrid\n\n");
  fprintf(wastl,"%%start of base pair probability data\n");

  EPS_print_bpp_data(wastl, pl, mf);

  EPS_footer(wastl);

  fclose(wastl);
  return 1; /* success */
}


PUBLIC int
vrna_plot_dp_EPS( const char              *filename,
                  const char              *sequence,
                  vrna_plist_t            *upper,
                  vrna_plist_t            *lower,
                  vrna_dotplot_auxdata_t  *auxdata,
                  unsigned int            options){

  char            **lintoptitle,**linbottomtitle,**linlefttitle,**linrighttitle,
                  *c, *comment, *title;
  int             pl_size, gq_num, i, lintop_num, lintop_size, linbottom_num,
                  linbottom_size, linleft_num, linleft_size, linright_num,
                  linright_size, bbox[4];
  double          tmp;
  FILE            *fh;
  plist           *pl1;
  vrna_data_lin_t **lintop, **linbottom, **linleft, **linright, *ud_lin, *access;

  fh = fopen(filename, "w");
  if(!fh){
    vrna_message_warning("can't open %s for dot plot", filename);
    return 0; /* return 0 for failure */
  }

  comment = title = NULL;

  lintop_num      = 0;
  lintop_size     = 5;
  linbottom_num   = 0;
  linbottom_size  = 5;
  linleft_num     = 0;
  linleft_size    = 5;
  linright_num    = 0;
  linright_size   = 5;
  bbox[0]         = 0;
  bbox[1]         = 0;
  bbox[2]         = 700;
  bbox[3]         = 720;
  access          = NULL;
  ud_lin          = NULL;
  lintop          = (vrna_data_lin_t **)vrna_alloc(sizeof(vrna_data_lin_t *) * lintop_size);
  lintoptitle     = (char **)vrna_alloc(sizeof(char *) * lintop_size);
  linbottom       = (vrna_data_lin_t **)vrna_alloc(sizeof(vrna_data_lin_t *) * linbottom_size);
  linbottomtitle  = (char **)vrna_alloc(sizeof(char *) * linbottom_size);
  linleft         = (vrna_data_lin_t **)vrna_alloc(sizeof(vrna_data_lin_t *) * linleft_size);
  linlefttitle    = (char **)vrna_alloc(sizeof(char *) * linleft_size);
  linright        = (vrna_data_lin_t **)vrna_alloc(sizeof(vrna_data_lin_t *) * linright_size);
  linrighttitle   = (char **)vrna_alloc(sizeof(char *) * linright_size);

  /* prepare linear data and retrieve number of additional linear data lines to correct bounding box */
  if(options & VRNA_PLOT_PROBABILITIES_UD_LIN){
    ud_lin = plist_to_ud_motif_prob(upper, strlen(sequence));
    if(ud_lin){
      dp_add_lindata(lintop, ud_lin, lintoptitle, "Protein binding", lintop_num, lintop_size);
      dp_add_lindata(linleft, ud_lin, linlefttitle, "Protein binding", linleft_num, linleft_size);
      dp_add_lindata(linbottom, ud_lin, linbottomtitle, "Protein binding", linbottom_num, linbottom_size);
      dp_add_lindata(linright, ud_lin, linrighttitle, "Protein binding", linright_num, linright_size);
    }
  }

  if(options & VRNA_PLOT_PROBABILITIES_ACC){
    access = plist_to_accessibility(upper, strlen(sequence));
    dp_add_lindata(lintop, access, lintoptitle, "Accessibility", lintop_num, lintop_size);
  }

  if(auxdata){
    if(auxdata->top){
      for(i = 0; auxdata->top[i]; i++){
        dp_add_lindata(lintop, auxdata->top[i], lintoptitle, auxdata->top_title[i], lintop_num, lintop_size);
      }
    }
    if(auxdata->bottom){
      for(i = 0; auxdata->bottom[i]; i++){
        dp_add_lindata(linbottom, auxdata->bottom[i], linbottomtitle, auxdata->bottom_title[i], linbottom_num, linbottom_size);
      }
    }
    if(auxdata->left){
      for(i = 0; auxdata->left[i]; i++){
        dp_add_lindata(linleft, auxdata->left[i], linlefttitle, auxdata->left_title[i], linleft_num, linleft_size);
      }
    }
    if(auxdata->right){
      for(i = 0; auxdata->right[i]; i++){
        dp_add_lindata(linright, auxdata->right[i], linrighttitle, auxdata->right_title[i], linright_num, linright_size);
      }
    }
  }

  dp_finalize_lindata(lintop, lintoptitle, lintop_num);
  dp_finalize_lindata(linbottom, linbottomtitle, linbottom_num);
  dp_finalize_lindata(linleft, linlefttitle, linleft_num);
  dp_finalize_lindata(linright, linrighttitle, linright_num);

  if(auxdata){
    comment = auxdata->comment;
    title   = (auxdata->title) ? strdup(auxdata->title) : NULL;
  }

  if(!title){
    title = strdup(filename);
    if((c=strrchr(title, '_')))
      *c='\0';
  }

  /* start printing postscript file */
  EPS_print_header(fh, bbox, comment, DP_MACRO_ALL);
  EPS_print_title(fh, title);
  EPS_print_seq(fh, sequence);

  fprintf(fh,"%% BEGIN linear data array\n\n");
  EPS_print_linear_data_top(fh, (const char **)lintoptitle, lintop);
  EPS_print_linear_data_left(fh, (const char **)linlefttitle, linleft);
  EPS_print_linear_data_bottom(fh, (const char **)linbottomtitle, linbottom);
  EPS_print_linear_data_right(fh, (const char **)linrighttitle, linright);
  fprintf(fh,"%% END linear data arrays\n");

  fprintf(fh, "\n%%Finally, prepare canvas\n\n"
              "%%draw title\ndrawTitle\n\n"
              "%%prepare coordinate system, draw grid and sequence\n"
              "/Helvetica findfont 0.95 scalefont setfont\n\n"
              "%%prepare coordinate system\nprepareCoords\n\n"
              "%%draw sequence arround grid\ndrawseq\n\n"
              "%%draw grid\ndrawgrid\n\n"
              "%%draw auxiliary linear data (if available)\ndrawData\n\n");

  fprintf(fh,"%%data (commands) starts here\n");

  if(options & VRNA_PLOT_PROBABILITIES_SD){
    EPS_print_sd_data(fh, upper, lower);
  }

  if(options & VRNA_PLOT_PROBABILITIES_UD){
    EPS_print_ud_data(fh, upper, lower);
  }


  EPS_print_sc_motif_data(fh, upper, lower);
  EPS_print_bpp_data(fh, upper, lower);

  EPS_footer(fh);

  fclose(fh);
  free(lintoptitle);
  free(lintop);
  free(linbottomtitle);
  free(linbottom);
  free(linlefttitle);
  free(linleft);
  free(linrighttitle);
  free(linright);
  free(access);
  free(ud_lin);
  free(title);

  return 1; /* success */
}


PUBLIC int
PS_color_dot_plot_turn( char *seq,
                        cpair *pi,
                        char *wastlfile,
                        int winSize) {

  /* produce color PostScript dot plot from cpair */

  FILE *wastl;
  int i;

  wastl = PS_dot_common(seq, cut_point, wastlfile, NULL, winSize, DP_MACRO_NONE);
  if (wastl==NULL)
    return 0; /* return 0 for failure */

  fprintf(wastl, "/hsb {\n"
          "dup 0.3 mul 1 exch sub sethsbcolor\n"
          "} bind def\n\n"
          "%%BEGIN DATA\n");

  if(winSize > 0)
    fprintf(wastl, "\n%%draw the grid\ndrawgrid_turn\n\n");
  else
    fprintf(wastl,  "\n%%draw the grid\ndrawgrid\n\n");
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

   EPS_footer(wastl);

   fclose(wastl);
   return 1; /* success */
}


PUBLIC int
PS_dot_plot_turn( char *seq,
                  plist *pl,
                  char *wastlfile,
                  int winSize) {

  /* produce color PostScript dot plot from cpair */

  FILE *wastl;
  int i;

  wastl = PS_dot_common(seq, cut_point, wastlfile, NULL, winSize, DP_MACRO_NONE);
  if (wastl==NULL)
    return 0; /* return 0 for failure */

  if(winSize > 0)
    fprintf(wastl, "\n%%draw the grid\ndrawgrid_turn\n\n");
  else
    fprintf(wastl,  "\n%%draw the grid\ndrawgrid\n\n");
  fprintf(wastl,"%%start of base pair probability data\n");
  /* print boxes */
  i=0;
  if (pl) {
    while (pl[i].j>0) {
      fprintf(wastl,"%d %d %1.4f ubox\n",
              pl[i].i, pl[i].j, sqrt(pl[i].p));
      i++;
    }
  }

  EPS_footer(wastl);

  fclose(wastl);
  return 1; /* success */
}


/*
#####################################
# BEGIN OF STATIC HELPER FUNCTIONS  #
#####################################
*/

PRIVATE void
EPS_footer(FILE *eps){

   fprintf(eps,"showpage\n"
               "end\n"
               "%%%%EOF\n");
}


PRIVATE void
EPS_print_title(FILE *eps, const char *title){

  unsigned int i;

  fprintf(eps,"/DPtitle {\n  (%s)\n} def\n\n", title);
}


PRIVATE void
EPS_print_seq(FILE *eps, const char *sequence){

  unsigned int i;

  fprintf(eps,"/sequence { (\\\n");
  for(i = 0; i < strlen(sequence); i += 255)
    fprintf(eps, "%.255s\\\n", sequence + i);
  fprintf(eps,") } def\n\n"
              "/len { sequence length } bind def\n\n");
}


PRIVATE void
EPS_print_header( FILE          *eps,
                  int           bbox[4],
                  const char    *comment,
                  unsigned int  options){

  fprintf(eps,
          "%%!PS-Adobe-3.0 EPSF-3.0\n"
          "%%%%Title: RNA Dot Plot\n"
          "%%%%Creator: ViennaRNA-%s\n"
          "%%%%CreationDate: %s", VERSION, vrna_time_stamp());

  fprintf(eps, "%%%%BoundingBox: %d %d %d %d\n", bbox[0], bbox[1], bbox[2], bbox[3]);

  fprintf(eps,
          "%%%%DocumentFonts: Helvetica\n"
          "%%%%Pages: 1\n"
          "%%%%EndComments\n\n"
          "%%Options: %s\n", option_string());

  if(comment)
    fprintf(eps,"%% %s\n",comment);

  fprintf(eps,"%s", RNAdp_prolog);

  /* add auxiliary macros */
  if(options & DP_MACRO_SD){  /* gquads et al. */
    fprintf(eps,"%s", RNAdp_sd);
  }

  if(options & DP_MACRO_UD){  /* unstructured domains */
    fprintf(eps,"%s", RNAdp_ud);
  }

  if(options & DP_MACRO_SC_MOTIFS){ /* soft constraint motifs */
    fprintf(eps,"%s", RNAdp_sc_motifs);
  }

  if(options & DP_MACRO_LINEAR_DATA){ /* linear data */
    fprintf(eps,"%s", RNAdp_linear_data);
  }

  fprintf(eps, "end\n%%EndProlog\n\nDPdict begin\n\n");
}


PRIVATE void
EPS_print_sd_data(FILE          *eps,
                  vrna_plist_t  *pl,
                  vrna_plist_t  *mf){

  int     pl_size, gq_num;
  double  tmp;
  plist   *pl1;

  /* sort the plist to bring all gquad triangles to the front */
  for(gq_num = pl_size = 0, pl1 = pl; pl1->i > 0; pl1++, pl_size++)
    if(pl1->type == VRNA_PLIST_TYPE_GQUAD) gq_num++;

  qsort(pl, pl_size, sizeof(plist), sort_plist_by_type_desc);

  /* sort all gquad triangles by probability to bring lower probs to the front */
  qsort(pl, gq_num, sizeof(plist), sort_plist_by_prob_asc);

  /* print triangles for g-quadruplexes in upper half */
  fprintf(eps,"\n%%start of quadruplex data\n");

  for (pl1=pl; pl1->i > 0; pl1++) {
    if(pl1->type == VRNA_PLIST_TYPE_GQUAD){
      tmp = sqrt(pl1->p);
      fprintf(eps, "%d %d %1.9f utri\n", pl1->i, pl1->j, tmp);
    }
  }
}


PRIVATE void
EPS_print_sc_motif_data(FILE          *eps,
                        vrna_plist_t  *pl,
                        vrna_plist_t  *mf){

  int     pl_size;
  double  tmp;
  plist   *pl1;

  /* print triangles for hairpin loop motifs in upper half */
  fprintf(eps,"\n%%start of Hmotif data\n");
  for (pl1=pl; pl1->i > 0; pl1++) {
    if(pl1->type == VRNA_PLIST_TYPE_H_MOTIF){
      tmp = sqrt(pl1->p);
      fprintf(eps, "%d %d %1.9f uHmotif\n", pl1->i, pl1->j, tmp);
    }
  }
  for (pl1=mf; pl1->i > 0; pl1++) {
    if(pl1->type == VRNA_PLIST_TYPE_H_MOTIF){
      tmp = sqrt(pl1->p);
      fprintf(eps, "%d %d %1.9f lHmotif\n", pl1->i, pl1->j, tmp);
    }
  }

  /* print triangles for interior loop motifs in upper half */
  fprintf(eps,"\n%%start of Imotif data\n");
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
        fprintf(eps, "%d %d %d %d %1.9f uImotif\n", a, b, pl1->i, pl1->j, ppp);
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
        fprintf(eps, "%d %d %d %d %1.9f lImotif\n", a, b, pl1->i, pl1->j, ppp);
        a = b = 0;
      }
    }
  }
}


PRIVATE void
EPS_print_bpp_data( FILE          *eps,
                    vrna_plist_t  *pl,
                    vrna_plist_t  *mf){

  int     pl_size;
  double  tmp;
  plist   *pl1;

  fprintf(eps,"%%start of base pair probability data\n");

  /* print boxes in upper right half*/
  for (pl1 = pl; pl1->i>0; pl1++) {
    tmp = sqrt(pl1->p);
    if(pl1->type == VRNA_PLIST_TYPE_BASEPAIR)
        fprintf(eps,"%d %d %1.9f ubox\n", pl1->i, pl1->j, tmp);
  }


  /* print boxes in lower left half (mfe) */
  for (pl1=mf; pl1->i>0; pl1++) {
    tmp = sqrt(pl1->p);
    if(pl1->type == VRNA_PLIST_TYPE_BASEPAIR)
      fprintf(eps,"%d %d %1.7f lbox\n", pl1->i, pl1->j, tmp);
  }
}


PRIVATE void
EPS_print_ud_data(FILE          *eps,
                  vrna_plist_t  *pl,
                  vrna_plist_t  *mf){

  int     pl_size;
  double  tmp;
  plist   *pl1;

  /* print triangles for unstructured domain motifs in upper half */
  fprintf(eps,"\n%%start of unstructured domain motif data\n");
  for(pl1 = pl; pl1->i > 0; pl1++){
    if(pl1->type == VRNA_PLIST_TYPE_UD_MOTIF){
      tmp = sqrt(pl1->p);
      fprintf(eps, "%d %d %1.9f uUDmotif\n", pl1->i, pl1->j, tmp);
    }
  }
  for(pl1 = mf; pl1->i > 0; pl1++){
    if(pl1->type == VRNA_PLIST_TYPE_UD_MOTIF){
      tmp = sqrt(pl1->p);
      fprintf(eps, "%d %d %1.9f lUDmotif\n", pl1->i, pl1->j, tmp);
    }
  }
}


PRIVATE void
EPS_print_linear_data_top(FILE            *eps,
                          const char      **data_title,
                          vrna_data_lin_t **data){

  EPS_print_linear_data(eps, "topData", data_title, data);
}


PRIVATE void
EPS_print_linear_data_left( FILE            *eps,
                            const char      **data_title,
                            vrna_data_lin_t **data){

  EPS_print_linear_data(eps, "leftData", data_title, data);
}


PRIVATE void
EPS_print_linear_data_bottom( FILE            *eps,
                              const char      **data_title,
                              vrna_data_lin_t **data){

  EPS_print_linear_data(eps, "bottomData", data_title, data);
}


PRIVATE void
EPS_print_linear_data_right(FILE            *eps,
                            const char      **data_title,
                            vrna_data_lin_t **data){

  EPS_print_linear_data(eps, "rightData", data_title, data);
}


PRIVATE void
EPS_print_linear_data(FILE            *eps,
                      const char      *varname,
                      const char      **data_title,
                      vrna_data_lin_t **data){

  int             n, i;
  vrna_data_lin_t *ptr;

  /* count number of data sets */
  for(n = 0; data_title[n]; n++);

  fprintf(eps, "/%s [\n", varname);
  for(i = 0; i < n; i++){
    fprintf(eps, "[ (%s)\n", data_title[i]);
    for(ptr = data[i]; ptr->position > 0; ptr++){
      if((ptr->color.hue + ptr->color.sat + ptr->color.bri) == 0.){
        fprintf(eps, "  [ %d %1.9f ]\n", ptr->position, ptr->value);
      } else {
        fprintf(eps, "  [ %d %1.9f %1.4f %1.4f %1.4f]\n", ptr->position, ptr->value, ptr->color.hue, ptr->color.sat, ptr->color.bri);
      }
    }
    fprintf(eps, "]\n");
  }
  fprintf(eps, "] def\n\n");
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


PRIVATE vrna_data_lin_t *
plist_to_accessibility(plist *pl, unsigned int length){

  int   n, i;
  plist *ptr;

  vrna_data_lin_t *data = NULL;

  data = (vrna_data_lin_t *)vrna_alloc(sizeof(vrna_data_lin_t) * (length + 1));

  for(ptr = pl; ptr->i > 0; ptr++){
    if(ptr->type == VRNA_PLIST_TYPE_BASEPAIR){
      data[ptr->i - 1].value += ptr->p;
      data[ptr->j - 1].value += ptr->p;
    }
  }

  /* go through the entire list and square-root probabilities again */
  for(i = 0; i < length; i++){
    data[i].position = i + 1;
    data[i].value    = sqrt((double)(1. - data[i].value));
  }

  data[length].position = 0; /* end marker */

  return data;
}

PRIVATE vrna_data_lin_t *
plist_to_ud_motif_prob(plist *pl, unsigned int length){

  int   n, i;
  plist *ptr;

  vrna_data_lin_t *data = NULL;

  data = (vrna_data_lin_t *)vrna_alloc(sizeof(vrna_data_lin_t) * (length + 1));

  for(ptr = pl; ptr->i > 0; ptr++){
    if(ptr->type == VRNA_PLIST_TYPE_UD_MOTIF){
      for(i = ptr->i; i <= ptr->j; i++){
        data[i - 1].value += ptr->p;
      }
    }
  }

  /*  go through the entire list, remove entries with 0 probability,
      and square-root probabilities again
  */
  unsigned int actual_length  = length;
  unsigned int actual_pos     = 1;
  for(i = 0; i < actual_length; i++, actual_pos++){
    if(data[i].value == 0.){
      memmove(data + i, data + i + 1, sizeof(vrna_data_lin_t) * (actual_length - i));
      actual_length--;
      i--;
    } else {
      data[i].position  = actual_pos;
      data[i].value     = sqrt(data[i].value);
      data[i].color.hue = 0.6;
      data[i].color.sat = 0.8;
      data[i].color.bri = 0.95;
    }
  }

  if(actual_length > 0){
    data[actual_length].position = 0; /* end marker */
    data = (vrna_data_lin_t *)vrna_realloc(data, sizeof(vrna_data_lin_t) * (actual_length + 1));
    return data;
  } else {
    free(data);
    return NULL;
  }
}



PRIVATE FILE *
PS_dot_common(const char *seq,
              int cp,
              const char *wastlfile,
              char *comment,
              int winsize,
              unsigned int options){

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

  int bbox[4];
  if(winsize > 0){
    bbox[0] = 66;
    bbox[1] = 530;
    bbox[2] = 520;
    bbox[3] = 650;
  } else {
    bbox[0] = 66;
    bbox[1] = 211;
    bbox[2] = 518;
    bbox[3] = 662;
  }

  EPS_print_header(wastl, bbox, comment, options);

  EPS_print_title(wastl, name);

  EPS_print_seq(wastl, seq);

  if (winsize>0)
    fprintf(wastl,"/winSize %d def\n",winsize);

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
    fprintf(wastl,"drawseq\n");

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

