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
#include "ViennaRNA/fold_vars.h"
#include "ViennaRNA/utils/basic.h"
#include "ViennaRNA/utils/alignments.h"
#include "ViennaRNA/alphabet.h"
#include "ViennaRNA/plotting/alignments.h"

#include "ViennaRNA/static/templates_postscript.h"

#include "ViennaRNA/plotting/ps_helpers.inc"

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
vrna_file_PS_aln(const char   *filename,
                 const char   **seqs,
                 const char   **names,
                 const char   *structure,
                 unsigned int columns)
{
  return vrna_file_PS_aln_slice(filename,
                                seqs,
                                names,
                                structure,
                                0,
                                0,
                                0,
                                columns);
}


PUBLIC int
vrna_file_PS_aln_slice(const char   *filename,
                       const char   **seqs,
                       const char   **names,
                       const char   *structure,
                       unsigned int start,
                       unsigned int end,
                       int          offset,
                       unsigned int columns)
{
  /* produce PS sequence alignment color-annotated by consensus structure */

  int       N, i, j, k, x, y, tmp, columnWidth, bbox[4];
  char      *tmpBuffer, *ssEscaped, *ruler, *cons, *substructure;
  char      c;
  float     fontWidth, fontHeight, imageHeight, imageWidth, tmpColumns;
  int       length, maxName, maxNum, currPos, num;
  float     lineStep, blockStep, consStep, ssStep, rulerStep, nameStep, numberStep;
  float     maxConsBar, startY, namesX, seqsX, currY;
  float     score, barHeight, xx, yy;
  int       match, block;
  FILE      *outfile;
  short     *pair_table;
  char      *colorMatrix[6][3] = {
    { "0.0 1",  "0.0 0.6",  "0.0 0.2"  },   /* red    */
    { "0.16 1", "0.16 0.6", "0.16 0.2" },   /* ochre  */
    { "0.32 1", "0.32 0.6", "0.32 0.2" },   /* turquoise */
    { "0.48 1", "0.48 0.6", "0.48 0.2" },   /* green  */
    { "0.65 1", "0.65 0.6", "0.65 0.2" },   /* blue   */
    { "0.81 1", "0.81 0.6", "0.81 0.2" }    /* violet */
  };

  vrna_md_t md;

  set_model_details(&md);

  outfile = fopen(filename, "w");

  if (outfile == NULL) {
    vrna_message_warning("can't open file %s - not doing alignment plot\n", filename);
    return 0;
  }

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

  if (start == 0)
    start = 1;

  if (end == 0)
    end = (int)strlen(seqs[0]);

  /* Number of columns of the alignment */
  length = end - start + 1;

  substructure          = (char *)vrna_alloc(sizeof(char) * (length + 1));
  substructure          = memcpy(substructure, structure + start - 1, sizeof(char) * length);
  substructure[length]  = '\0';

  /* remove unbalanced brackets ? */

  /* Display long alignments in blocks of this size */
  columnWidth = (columns == 0) ? length : columns;

  /*
   *  Allocate memory for various strings, length*2 is (more than)
   *  enough for all of them
   */
  tmpBuffer = (char *)vrna_alloc((unsigned)MAX2(length * 2, columnWidth) + 1);
  ssEscaped = (char *)vrna_alloc((unsigned)length * 2);
  ruler     = (char *)vrna_alloc((unsigned)length * 2);


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
  snprintf(tmpBuffer, length, "%d", start + length + offset);
  maxNum = strlen(tmpBuffer);

  /* Calculate bounding box */
  tmpColumns = columnWidth;
  if (length < columnWidth)
    tmpColumns = length;

  imageWidth = ceil(namesX +
                    (maxName + tmpColumns + maxNum) * fontWidth +
                    2 * nameStep +
                    fontWidth +
                    numberStep);

  imageHeight = startY +
                ceil((float)length / columnWidth) *
                ((N + 2) * lineStep + blockStep + consStep + ssStep + rulerStep);

  bbox[0] = bbox[1] = 0;
  bbox[2] = (int)imageWidth;
  bbox[3] = (int)imageHeight;

  /* Write postscript header including correct bounding box */
  print_PS_header(outfile,
                  "ViennaRNA Package - Alignment",
                  bbox,
                  &md,
                  NULL,
                  "ALNdict",
                  PS_MACRO_ALN_BASE);

  fprintf(outfile, "0 %d translate\n"
          "1 -1 scale\n"
          "/Courier findfont\n"
          "[10 0 0 -10 0 0] makefont setfont\n",
          (int)imageHeight);

  /*
   * Create ruler and secondary structure lines
   * Init all with dots
   */
  memset(ruler, '.', sizeof(char) * length);

  for (i = 0; i < length; i++) {
    /* Write number every 10th position, leave out block breaks */
    if (((i + start + offset) % 10 == 0) && (i % columnWidth != 0)) {
      snprintf(tmpBuffer, length, "%d", i + start + offset);
      num = strlen(tmpBuffer);
      if (i + num <= length)
        memcpy(ruler + i, tmpBuffer, num);
    }
  }
  ruler[length] = '\0';

  pair_table = vrna_ptable_from_string(substructure,
                                       VRNA_BRACKETS_RND |
                                       VRNA_BRACKETS_ANG |
                                       VRNA_BRACKETS_SQR);

  /*
   *  note: The substructrue and the corresponding pair_table
   *  are relative to the 'start' position, while the sequence(s)
   *  are not! Therefore we need to shift structure coordinates
   *  accordingly.
   */
  int shift = start - 1;

  pair_table -= shift;

  /*
   * Draw color annotation first
   * Repeat for all pairs
   */
  for (i = start; i <= end; i++) {
    j = pair_table[i] + shift;
    if ((j > i) && (j <= end)) {
      /* Repeat for open and closing position */
      for (k = 0; k < 2; k++) {
        unsigned int  tt;
        int           pairings, nonpair, s, col;
        int           ptype[8] = {
          0, 0, 0, 0, 0, 0, 0, 0
        };
        char          *color;
        col   = (k == 0) ? i - shift - 1 : j - shift - 1;
        block = ceil((float)(col + 1) / columnWidth);
        xx    = seqsX + (col - (block - 1) * columnWidth) * fontWidth;

        /* Repeat for each sequence */
        for (s = pairings = nonpair = 0; s < N; s++) {
          tt = md.pair[vrna_nucleotide_encode(seqs[s][i - 1], &md)]
               [vrna_nucleotide_encode(seqs[s][j - 1], &md)];
          ptype[tt]++;
        }

        for (pairings = 0, s = 1; s <= 7; s++)
          if (ptype[s])
            pairings++;

        nonpair = ptype[0];
        if (nonpair <= 2) {
          color = colorMatrix[pairings - 1][nonpair];
          for (s = 0; s < N; s++) {
            yy = startY + (block - 1) * (lineStep * (N + 2) + blockStep + consStep + rulerStep) +
                 ssStep * (block) + (s + 1) * lineStep;

            /* Color according due color information in pi-array, only if base pair is possible */
            tt = md.pair[vrna_nucleotide_encode(seqs[s][i - 1], &md)]
                 [vrna_nucleotide_encode(seqs[s][j - 1], &md)];
            if (tt)
              fprintf(outfile, "%.1f %.1f %.1f %.1f %s box\n",
                      xx, yy - 1, xx + fontWidth, yy + fontHeight + 1, color);
          }
        }
      }
    }
  }

  pair_table += shift;
  free(pair_table);

  /* Process rest of the output in blocks of columnWidth */

  currY   = startY;
  currPos = 0;

  cons = vrna_aln_consensus_sequence(seqs, &md);

  while (currPos < length) {
    /* Display secondary structure line */
    fprintf(outfile, "0 setgray\n");
    strncpy(tmpBuffer, substructure + currPos, columnWidth);
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
      unsigned int max_len = columnWidth;
      if (length - currPos < max_len)
        max_len = length - currPos;

      strncpy(tmpBuffer, seqs[i] + currPos + shift, max_len);
      tmpBuffer[max_len] = '\0';

      match = 0;
      for (j = 0; j < (currPos + strlen(tmpBuffer)); j++)
        if (seqs[i][j + shift] != '-')
          match++;

      fprintf(outfile, "(%s) %.1f %.1f string\n", names[i], namesX, currY);
      fprintf(outfile, "(%s) %.1f %.1f string\n", tmpBuffer, seqsX, currY);
      fprintf(outfile,
              "(%i) %.1f %.1f string\n",
              match,
              seqsX + fontWidth * (strlen(tmpBuffer)) + numberStep,
              currY);
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
        if (cons[i + shift] == toupper(seqs[j][i + shift]))
          match++;

        if (cons[i + shift] == 'U' && toupper(seqs[j][i + shift]) == 'T')
          match++;

        if (cons[i + shift] == 'T' && toupper(seqs[j][i + shift]) == 'U')
          match++;
      }
      score = (float)(match - 1) / (N - 1);

      if (cons[i + shift] == '-' ||
          cons[i + shift] == '_' ||
          cons[i + shift] == '.')
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

  print_PS_footer(outfile);

  fclose(outfile);

  free(tmpBuffer);
  free(ssEscaped);
  free(ruler);
  free(substructure);

  return 0;
}


#ifndef VRNA_DISABLE_BACKWARD_COMPATIBILITY

/*
 *###########################################
 *# deprecated functions below              #
 *###########################################
 */
PUBLIC int
PS_color_aln(const char *structure,
             const char *filename,
             const char *seqs[],
             const char *names[])
{
  return vrna_file_PS_aln(filename, seqs, names, structure, 60);
}


int
aliPS_color_aln(const char  *structure,
                const char  *filename,
                const char  *seqs[],
                const char  *names[])
{
  return vrna_file_PS_aln_slice(filename,
                                seqs,
                                names,
                                structure,
                                0,
                                0,
                                0,
                                100);
}


#endif
