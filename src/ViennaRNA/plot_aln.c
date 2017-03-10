/*
 *      PostScript output for Sequence / Structure Alignments
 *
 *               c  Ivo Hofacker, Peter F Stadler, Ronny Lorenz
 *                        Vienna RNA package
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
#include "ViennaRNA/aln_util.h"
#include "ViennaRNA/plot_aln.h"

/*
 #################################
 # PRIVATE MACROS                #
 #################################
 */

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

/*
 #################################
 # BEGIN OF FUNCTION DEFINITIONS #
 #################################
 */
PUBLIC int
PS_color_aln(const char *structure,
             const char *filename,
             const char *seqs[],
             const char *names[])
{
  return vrna_file_PS_aln(filename, seqs, names, structure, 60);
}


PUBLIC int
vrna_file_PS_aln(const char *filename,
                 const char **seqs,
                 const char **names,
                 const char *structure,
                 int        columns)
{
  return vrna_file_PS_aln_sub(filename, seqs, names, structure, 1, -1, columns);
}


PUBLIC int
vrna_file_PS_aln_sub(const char *filename,
                     const char **seqs,
                     const char **names,
                     const char *structure,
                     int        start,
                     int        end,
                     int        columns)
{
  /* produce PS sequence alignment color-annotated by consensus structure */

  int         N, i, j, k, x, y, tmp, columnWidth;
  char        *tmpBuffer, *ssEscaped, *ruler, *cons;
  char        c;
  float       fontWidth, fontHeight, imageHeight, imageWidth, tmpColumns;
  int         length, maxName, maxNum, currPos;
  float       lineStep, blockStep, consStep, ssStep, rulerStep, nameStep, numberStep;
  float       maxConsBar, startY, namesX, seqsX, currY;
  float       score, barHeight, xx, yy;
  int         match, block;
  FILE        *outfile;
  short       *pair_table;
  char        *colorMatrix[6][3] = {
    { "0.0 1",  "0.0 0.6",  "0.0 0.2"  }, /* red    */
    { "0.16 1", "0.16 0.6", "0.16 0.2" }, /* ochre  */
    { "0.32 1", "0.32 0.6", "0.32 0.2" }, /* turquoise */
    { "0.48 1", "0.48 0.6", "0.48 0.2" }, /* green  */
    { "0.65 1", "0.65 0.6", "0.65 0.2" }, /* blue   */
    { "0.81 1", "0.81 0.6", "0.81 0.2" } /* violet */
  };

  const char  *alnPlotHeader =
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

  set_model_details(&md);

  outfile = fopen(filename, "w");

  if (outfile == NULL) {
    vrna_message_warning("can't open file %s - not doing alignment plot\n", filename);
    return 0;
  }

  columnWidth = (columns <= 0) ? 60 : columns;  /* Display long alignments in blocks of this size */
  fontWidth   = 6;                              /* Font metrics */
  fontHeight  = 6.5;
  lineStep    = fontHeight + 2;                 /* distance between lines */
  blockStep   = 3.5 * fontHeight;               /* distance between blocks */
  consStep    = fontHeight * 0.5;               /* distance between alignment and conservation curve */
  ssStep      = 2;                              /* distance between secondary structure line and sequences */
  rulerStep   = 2;                              /* distance between sequences and ruler */
  nameStep    = 3 * fontWidth;                  /* distance between names and sequences */
  numberStep  = fontWidth;                      /* distance between sequeces and numbers */
  maxConsBar  = 2.5 * fontHeight;               /* Height of conservation curve */
  startY      = 2;                              /* "y origin" */
  namesX      = fontWidth;                      /* "x origin" */

  /* Number of columns of the alignment */
  length = strlen(seqs[0]);

  /* Allocate memory for various strings, length*2 is (more than)
   *     enough for all of them */
  tmpBuffer = (char *)vrna_alloc((unsigned)MAX2(length * 2, columnWidth) + 1);
  ssEscaped = (char *)vrna_alloc((unsigned)length * 2);
  ruler     = (char *)vrna_alloc((unsigned)length * 2);

  pair_table = vrna_ptable(structure);
  /* Get length of longest name and count sequences in alignment*/

  for (i = maxName = N = 0; names[i] != NULL; i++) {
    N++;
    tmp = strlen(names[i]);
    if (tmp > maxName)
      maxName = tmp;
  }


  /* x-coord. where sequences start */
  seqsX = namesX + maxName * fontWidth + nameStep;

  /* calculate number of digits of the alignment length */
  snprintf(tmpBuffer, length, "%i", length);
  maxNum = strlen(tmpBuffer);


  /* Calculate bounding box */
  tmpColumns = columnWidth;
  if (length < columnWidth)
    tmpColumns = length;

  imageWidth  = ceil(namesX + (maxName + tmpColumns + maxNum) * fontWidth + 2 * nameStep + fontWidth + numberStep);
  imageHeight = startY + ceil((float)length / columnWidth) * ((N + 2) * lineStep + blockStep + consStep + ssStep + rulerStep);

  /* Write postscript header including correct bounding box */
  fprintf(outfile, alnPlotHeader, 0, 0, (int)imageWidth, (int)imageHeight, (int)imageHeight);

  /* Create ruler and secondary structure lines */
  i = 0;
  /* Init all with dots */
  for (i = 0; i < (length); i++)
    ruler[i] = '.';
  i = 0;
  for (i = 0; i < length; i++) {
    /* Write number every 10th position, leave out block breaks */
    if ((i + start) % 10 == 0 && (i + start) % columnWidth != 0) {
      snprintf(tmpBuffer, length, "%i", i + start);
      strncpy(ruler + i, tmpBuffer, strlen(tmpBuffer));
    }
  }
  ruler[length] = '\0';

  /* Draw color annotation first */
  /* Repeat for all pairs */
  for (i = 1; i <= length; i++) {
    if ((j = pair_table[i]) > i) {
      /* Repeat for open and closing position */
      for (k = 0; k < 2; k++) {
        int   pairings, nonpair, s, col;
        int   ptype[8] = {
          0, 0, 0, 0, 0, 0, 0, 0
        };
        char  *color;
        col   = (k == 0) ? i - 1 : j - 1;
        block = ceil((float)(col + 1) / columnWidth);
        xx    = seqsX + (col - (block - 1) * columnWidth) * fontWidth;
        /* Repeat for each sequence */
        for (s = pairings = nonpair = 0; s < N; s++)
          ptype[md.pair[vrna_nucleotide_encode(seqs[s][i - 1], &md)][vrna_nucleotide_encode(seqs[s][j - 1], &md)]]++;
        for (pairings = 0, s = 1; s <= 7; s++)
          if (ptype[s])
            pairings++;

        nonpair = ptype[0];
        if (nonpair <= 2) {
          color = colorMatrix[pairings - 1][nonpair];
          for (s = 0; s < N; s++) {
            yy = startY + (block - 1) * (lineStep * (N + 2) + blockStep + consStep + rulerStep) + ssStep * (block) + (s + 1) * lineStep;

            /* Color according due color information in pi-array, only if base pair is possible */
            if (md.pair[vrna_nucleotide_encode(seqs[s][i - 1], &md)][vrna_nucleotide_encode(seqs[s][j - 1], &md)])
              fprintf(outfile, "%.1f %.1f %.1f %.1f %s box\n",
                      xx, yy - 1, xx + fontWidth, yy + fontHeight + 1, color);
          }
        }
      }
    }
  }
  free(pair_table);

  /* Process rest of the output in blocks of columnWidth */

  currY   = startY;
  currPos = 0;

  cons = consensus(seqs);

  while (currPos < length) {
    /* Display secondary structure line */
    fprintf(outfile, "0 setgray\n");
    strncpy(tmpBuffer, structure + currPos, columnWidth);
    tmpBuffer[columnWidth] = '\0';

    x = 0;
    y = 0;
    while ((c = tmpBuffer[x])) {
      if (c == '.') {
        ssEscaped[y++] = '.';
      } else {
        ssEscaped[y++]  = '\\';
        ssEscaped[y++]  = c;
      }

      x++;
    }
    ssEscaped[y] = '\0';

    fprintf(outfile, "(%s) %.1f %.1f string\n", ssEscaped, seqsX, currY);
    currY += ssStep + lineStep;

    /* Display names, sequences and numbers */

    for (i = 0; i < N; i++) {
      strncpy(tmpBuffer, seqs[i] + currPos, columnWidth);
      tmpBuffer[columnWidth] = '\0';

      match = 0;
      for (j = 0; j < (currPos + strlen(tmpBuffer)); j++)
        if (seqs[i][j] != '-')
          match++;

      fprintf(outfile, "(%s) %.1f %.1f string\n", names[i], namesX, currY);
      fprintf(outfile, "(%s) %.1f %.1f string\n", tmpBuffer, seqsX, currY);
      fprintf(outfile, "(%i) %.1f %.1f string\n", match, seqsX + fontWidth * (strlen(tmpBuffer)) + numberStep, currY);
      currY += lineStep;
    }
    currY += rulerStep;
    strncpy(tmpBuffer, ruler + currPos, columnWidth);
    tmpBuffer[columnWidth] = '\0';
    fprintf(outfile, "(%s) %.1f %.1f string\n", tmpBuffer, seqsX, currY);

    currY += lineStep;
    currY += consStep;

    /*Display conservation bar*/

    fprintf(outfile, "0.6 setgray\n");
    for (i = currPos; (i < currPos + columnWidth && i < length); i++) {
      match = 0;
      for (j = 0; j < N; j++) {
        if (cons[i] == seqs[j][i])
          match++;

        if (cons[i] == 'U' && seqs[j][i] == 'T')
          match++;

        if (cons[i] == 'T' && seqs[j][i] == 'U')
          match++;
      }
      score = (float)(match - 1) / (N - 1);

      if (cons[i] == '-' ||
          cons[i] == '_' ||
          cons[i] == '.')
        score = 0;

      barHeight = maxConsBar * score;
      if (barHeight == 0)
        barHeight = 1;

      xx = seqsX + (i - (columnWidth * currPos / columnWidth)) * fontWidth;

      fprintf(outfile, "%.1f %.1f %.1f %.1f box2\n",
              xx,
              currY + maxConsBar - barHeight,
              xx + fontWidth,
              currY + maxConsBar);
    }

    currY   += blockStep;
    currPos += columnWidth;
  }
  free(cons);

  fprintf(outfile, "showpage\n");
  fclose(outfile);

  free(tmpBuffer);
  free(ssEscaped);
  free(ruler);

  return 0;
}


int
aliPS_color_aln(const char  *structure,
                const char  *filename,
                const char  *seqs[],
                const char  *names[])
{
  /* produce PS sequence alignment color-annotated by consensus structure */

  int         N, i, j, k, x, y, tmp, columnWidth;
  char        *tmpBuffer, *ssEscaped, *ruler, *cons;
  char        c;
  float       fontWidth, fontHeight, imageHeight, imageWidth, tmpColumns;
  int         length, maxName, maxNum, currPos;
  float       lineStep, blockStep, consStep, ssStep, rulerStep, nameStep, numberStep;
  float       maxConsBar, startY, namesX, seqsX, currY;
  float       score, barHeight, xx, yy;
  int         match, block;
  FILE        *outfile;
  short       *pair_table;
  char        *colorMatrix[6][3] = {
    { "0.0 1",  "0.0 0.6",  "0.0 0.2"  }, /* red    */
    { "0.16 1", "0.16 0.6", "0.16 0.2" }, /* ochre  */
    { "0.32 1", "0.32 0.6", "0.32 0.2" }, /* turquoise */
    { "0.48 1", "0.48 0.6", "0.48 0.2" }, /* green  */
    { "0.65 1", "0.65 0.6", "0.65 0.2" }, /* blue   */
    { "0.81 1", "0.81 0.6", "0.81 0.2" } /* violet */
  };

  const char  *alnPlotHeader =
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

  set_model_details(&md);

  outfile = fopen(filename, "w");
  if (outfile == NULL) {
    vrna_message_warning("can't open file %s - not doing alignment plot\n", filename);
    return 0;
  }

  columnWidth = 100;              /* Display long alignments in blocks of this size */
  fontWidth   = 6;                /* Font metrics */
  fontHeight  = 6.5;
  lineStep    = fontHeight + 2;   /* distance between lines */
  blockStep   = 3.5 * fontHeight; /* distance between blocks */
  consStep    = fontHeight * 0.5; /* distance between alignment and conservation curve */
  ssStep      = 2;                /* distance between secondary structure line and sequences */
  rulerStep   = 2;                /* distance between sequences and ruler */
  nameStep    = 3 * fontWidth;    /* distance between names and sequences */
  numberStep  = fontWidth;        /* distance between sequeces and numbers */
  maxConsBar  = 2.5 * fontHeight; /* Height of conservation curve */
  startY      = 2;                /* "y origin" */
  namesX      = fontWidth;        /* "x origin" */

  /* Number of columns of the alignment */
  length = strlen(seqs[0]);

  /* Allocate memory for various strings, length*2 is (more than)
   *     enough for all of them */
  tmpBuffer = (char *)vrna_alloc((unsigned)columnWidth + length * 2);
  ssEscaped = (char *)vrna_alloc((unsigned)length * 2);
  ruler     = (char *)vrna_alloc((unsigned)length * 2);
  /*   char * structur; */
  /*   structur = (char*) vrna_alloc((length+1)*sizeof(char)); */
  /*   structur = strdup(structure); */
  /*   for(i=0; i<length;i++){ */
  /*     if(structur[i] == '<') structur[i]='('; */
  /*     if(structur[i] == '>') structur[i]=')'; */
  /*   } */
  /*   structur[length]='\0';    */
  /*   printf("%s \n", structur); */
  pair_table = vrna_pt_ali_get(structure);
  /* Get length of longest name and count sequences in alignment*/

  for (i = maxName = N = 0; names[i] != NULL; i++) {
    N++;
    tmp = strlen(names[i]);
    if (tmp > maxName)
      maxName = tmp;
  }


  /* x-coord. where sequences start */
  seqsX = namesX + maxName * fontWidth + nameStep;

  /* calculate number of digits of the alignment length */
  snprintf(tmpBuffer, length, "%i", length);
  maxNum = strlen(tmpBuffer);


  /* Calculate bounding box */
  tmpColumns = columnWidth;
  if (length < columnWidth)
    tmpColumns = length;

  imageWidth  = ceil(namesX + (maxName + tmpColumns + maxNum) * fontWidth + 2 * nameStep + fontWidth + numberStep);
  imageHeight = startY + ceil((float)length / columnWidth) * ((N + 2) * lineStep + blockStep + consStep + ssStep + rulerStep);

  /* Write postscript header including correct bounding box */
  fprintf(outfile, alnPlotHeader, 0, 0, (int)imageWidth, (int)imageHeight, (int)imageHeight);

  /* Create ruler and secondary structure lines */
  i = 0;
  /* Init all with dots */
  for (i = 0; i < (length); i++)
    ruler[i] = '.';
  i = 0;
  for (i = 0; i < length; i++) {
    /* Write number every 10th position, leave out block breaks */
    if ((i + 1) % 10 == 0 && (i + 1) % columnWidth != 0) {
      snprintf(tmpBuffer, length, "%i", i + 1);
      strncpy(ruler + i, tmpBuffer, strlen(tmpBuffer));
    }
  }
  ruler[length] = '\0';

  /* Draw color annotation first */
  /* Repeat for all pairs */
  for (i = 1; i <= length; i++) {
    if ((j = pair_table[i]) > i) {
      /* Repeat for open and closing position */
      for (k = 0; k < 2; k++) {
        int   pairings, nonpair, s, col;
        int   ptype[8] = {
          0, 0, 0, 0, 0, 0, 0, 0
        };
        char  *color;
        col   = (k == 0) ? i - 1 : j - 1;
        block = ceil((float)(col + 1) / columnWidth);
        xx    = seqsX + (col - (block - 1) * columnWidth) * fontWidth;
        /* Repeat for each sequence */
        for (s = pairings = nonpair = 0; s < N; s++)
          ptype[md.pair[vrna_nucleotide_encode(seqs[s][i - 1], &md)][vrna_nucleotide_encode(seqs[s][j - 1], &md)]]++;
        for (pairings = 0, s = 1; s <= 7; s++)
          if (ptype[s])
            pairings++;

        nonpair = ptype[0];
        if (nonpair <= 2) {
          color = colorMatrix[pairings - 1][nonpair];
          for (s = 0; s < N; s++) {
            yy = startY + (block - 1) * (lineStep * (N + 2) + blockStep + consStep + rulerStep) + ssStep * (block) + (s + 1) * lineStep;

            /* Color according due color information in pi-array, only if base pair is possible */
            if (md.pair[vrna_nucleotide_encode(seqs[s][i - 1], &md)][vrna_nucleotide_encode(seqs[s][j - 1], &md)])
              fprintf(outfile, "%.1f %.1f %.1f %.1f %s box\n",
                      xx, yy - 1, xx + fontWidth, yy + fontHeight + 1, color);
          }
        }
      }
    }
  }
  free(pair_table);

  /* Process rest of the output in blocks of columnWidth */

  currY   = startY;
  currPos = 0;

  cons = consensus(seqs);

  while (currPos < length) {
    /* Display secondary structure line */
    fprintf(outfile, "0 setgray\n");
    strncpy(tmpBuffer, structure + currPos, columnWidth);
    tmpBuffer[columnWidth] = '\0';

    x = 0;
    y = 0;
    while ((c = tmpBuffer[x])) {
      if (c == '.') {
        ssEscaped[y++] = '.';
      } else {
        ssEscaped[y++]  = '\\';
        ssEscaped[y++]  = c;
      }

      x++;
    }
    ssEscaped[y] = '\0';

    fprintf(outfile, "(%s) %.1f %.1f string\n", ssEscaped, seqsX, currY);
    currY += ssStep + lineStep;

    /* Display names, sequences and numbers */

    for (i = 0; i < N; i++) {
      strncpy(tmpBuffer, seqs[i] + currPos, columnWidth);
      tmpBuffer[columnWidth] = '\0';

      match = 0;
      for (j = 0; j < (currPos + strlen(tmpBuffer)); j++)
        if (seqs[i][j] != '-')
          match++;

      fprintf(outfile, "(%s) %.1f %.1f string\n", names[i], namesX, currY);
      fprintf(outfile, "(%s) %.1f %.1f string\n", tmpBuffer, seqsX, currY);
      fprintf(outfile, "(%i) %.1f %.1f string\n", match, seqsX + fontWidth * (strlen(tmpBuffer)) + numberStep, currY);
      currY += lineStep;
    }
    currY += rulerStep;
    strncpy(tmpBuffer, ruler + currPos, columnWidth);
    tmpBuffer[columnWidth] = '\0';
    fprintf(outfile, "(%s) %.1f %.1f string\n", tmpBuffer, seqsX, currY);

    currY += lineStep;
    currY += consStep;

    /*Display conservation bar*/

    fprintf(outfile, "0.6 setgray\n");
    for (i = currPos; (i < currPos + columnWidth && i < length); i++) {
      match = 0;
      for (j = 0; j < N; j++) {
        if (cons[i] == seqs[j][i])
          match++;

        if (cons[i] == 'U' && seqs[j][i] == 'T')
          match++;

        if (cons[i] == 'T' && seqs[j][i] == 'U')
          match++;
      }
      score = (float)(match - 1) / (N - 1);

      if (cons[i] == '-' ||
          cons[i] == '_' ||
          cons[i] == '.')
        score = 0;

      barHeight = maxConsBar * score;
      if (barHeight == 0)
        barHeight = 1;

      xx = seqsX + (i - (columnWidth * currPos / columnWidth)) * fontWidth;

      fprintf(outfile, "%.1f %.1f %.1f %.1f box2\n",
              xx,
              currY + maxConsBar - barHeight,
              xx + fontWidth,
              currY + maxConsBar);
    }

    currY   += blockStep;
    currPos += columnWidth;
  }
  free(cons);
  fprintf(outfile, "showpage\n");
  fclose(outfile);
  free(tmpBuffer);
  free(ssEscaped);
  free(ruler);

  return 0;
}


#ifdef VRNA_BACKWARD_COMPAT

/*###########################################*/
/*# deprecated functions below              #*/
/*###########################################*/


#endif
