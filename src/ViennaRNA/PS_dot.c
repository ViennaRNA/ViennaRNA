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

const char *RNAdp_gquad_triangle =
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
"} bind def\n";

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
    if(pl1->type == 1) gq_num++;
  qsort(pl, pl_size, sizeof(plist), sort_plist_by_type_desc);
  /* sort all gquad triangles by probability to bring lower probs to the front */
  qsort(pl, gq_num, sizeof(plist), sort_plist_by_prob_asc);

  /* print triangles for g-quadruplexes in upper half */
  fprintf(wastl,"\n%%start of quadruplex data\n");
  for (pl1=pl; pl1->type == 1; pl1++) {
    tmp = sqrt(pl1->p);
    fprintf(wastl, "%d %d %1.9f utri\n", pl1->i, pl1->j, tmp);
  }

  fprintf(wastl, "\n%%draw the grid\ndrawgrid\n\n");
  fprintf(wastl,"%%start of base pair probability data\n");

  /* print boxes in upper right half*/
  for (; pl1->i>0; pl1++) {
    tmp = sqrt(pl1->p);
    if(pl1->type == 0)
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

int PS_dot_plot_turn(char *seq, struct plist *pl, char *wastlfile, int winSize) {
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
  char name[31], *c;
  int i;

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

  vrna_md_t     md;

  vrna_md_set_globals(&md);

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
  tmpBuffer = (char *) vrna_alloc((unsigned) MAX2(length*2,columnWidth)+1);
  ssEscaped=(char *) vrna_alloc((unsigned) length*2);
  ruler=(char *) vrna_alloc((unsigned) length*2);

  pair_table=vrna_pt_get(structure);
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
          ptype[md.pair[vrna_nucleotide_encode(seqs[s][i-1], &md)][vrna_nucleotide_encode(seqs[s][j-1], &md)]]++;
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
            if (md.pair[vrna_nucleotide_encode(seqs[s][i-1], &md)][vrna_nucleotide_encode(seqs[s][j-1], &md)]) {

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
        
  vrna_md_t md;

  vrna_md_set_globals(&md);

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
  nameStep=3*fontWidth;             /* distance between names and sequences */
  numberStep=fontWidth;      /* distance between sequeces and numbers */
  maxConsBar=2.5*fontHeight; /* Height of conservation curve */
  startY=2;                     /* "y origin" */
  namesX=fontWidth;             /* "x origin" */

  /* Number of columns of the alignment */
  length=strlen(seqs[0]);

  /* Allocate memory for various strings, length*2 is (more than)
         enough for all of them */
  tmpBuffer = (char *) vrna_alloc((unsigned) columnWidth + length*2 );
  ssEscaped=(char *) vrna_alloc((unsigned) length*2 );
  ruler=(char *) vrna_alloc((unsigned) length*2  );
/*   char * structur; */
/*   structur = (char*) vrna_alloc((length+1)*sizeof(char)); */
/*   structur = strdup(structure); */
/*   for(i=0; i<length;i++){ */
/*     if(structur[i] == '<') structur[i]='('; */
/*     if(structur[i] == '>') structur[i]=')'; */
/*   } */
/*   structur[length]='\0';    */
/*   printf("%s \n", structur); */
   pair_table=vrna_pt_ali_get(structure);
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
          ptype[md.pair[vrna_nucleotide_encode(seqs[s][i-1], &md)][vrna_nucleotide_encode(seqs[s][j-1], &md)]]++;
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
            if (md.pair[vrna_nucleotide_encode(seqs[s][i-1], &md)][vrna_nucleotide_encode(seqs[s][j-1], &md)]) {

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

#ifdef VRNA_BACKWARD_COMPAT

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
  pl = (struct plist *)vrna_alloc(maxl*sizeof(struct plist));
  k=0;
  /*make plist out of pr array*/
  for (i=1; i<length; i++)
    for (j=i+1; j<=length; j++) {
      if (pr[iindx[i]-j]<PMIN) continue;
      if (k>=maxl-1) {
        maxl *= 2;
        pl = (struct plist *)vrna_realloc(pl,maxl*sizeof(struct plist));
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
  mf = (struct plist *)vrna_alloc((mf_num+1)*sizeof(struct plist));
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

