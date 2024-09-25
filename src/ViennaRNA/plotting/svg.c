/*
 *      SVG output formats for ViennaRNA
 *
 *      c  Ivo Hofacker, Peter F Stadler, Ronny Lorenz
 *         Vienna RNA package
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
#include "ViennaRNA/utils/log.h"
#include "ViennaRNA/sequences/alphabet.h"
#include "ViennaRNA/model.h"
#include "ViennaRNA/fold_vars.h"
#include "ViennaRNA/plotting/layouts.h"
#include "ViennaRNA/plotting/utils.h"
#include "ViennaRNA/plotting/structures.h"
#include "ViennaRNA/plotting/RNApuzzler/RNApuzzler.h"
#include "ViennaRNA/plotting/RNApuzzler/RNAturtle.h"

#include "ViennaRNA/static/templates_svg.h"

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
rnaplot_SVG(const char         *ssfile,
            const char         *string,
            const char         *structure,
            vrna_plot_layout_t *layout,
            vrna_plot_data_t   *aux_data);


PRIVATE void
print_SVG_footer(FILE *fh);


PRIVATE void
print_SVG_footer_data(FILE *fh,
                      unsigned int  n,
                      const char    *sequence,
                      const char    *structure,
                      const short   *pt,
                      float         *X,
                      float         *Y);


PRIVATE void
print_SVG_bases(FILE          *fh,
                float         *X,
                float         *Y,
                const char    *string,
                unsigned int  n);


PRIVATE void
print_SVG_backbone(FILE         *fh,
                   const float  *X,
                   const float  *Y,
                   unsigned int n);


PRIVATE void
print_SVG_pairs(FILE          *fh,
                const short   *pt,
                const float   *X,
                const float   *Y,
                const float   *CX,
                const float   *CY,
                unsigned int  n,
                unsigned int  plot_type);


PRIVATE double *
transform_PS_arcs_to_SVG(unsigned int n,
                         double       *PSarcs);


static void
print_SVG_header(FILE         *fh,
                const float   scale[2],
                const float   transform[2],
                unsigned int  options);

/*
 #################################
 # BEGIN OF FUNCTION DEFINITIONS #
 #################################
 */


PUBLIC int
vrna_plot_structure_svg(const char          *filename,
                        const char          *sequence,
                        const char          *structure,
                        vrna_plot_layout_t  *layout,
                        vrna_plot_data_t    *data)
{
  int                 ret = 0;
  vrna_plot_layout_t  *lt;

  if ((sequence) &&
      (structure) &&
      (filename)) {

    lt = layout;

    if (lt == NULL) /* use global default layout */
      lt = vrna_plot_layout(structure, VRNA_PLOT_TYPE_DEFAULT);

    ret = rnaplot_SVG(filename, sequence, structure, lt, data);

    if (lt != layout)
      vrna_plot_layout_free(lt);
  }

  return ret;
}


/*
 #####################################
 # BEGIN OF STATIC HELPER FUNCTIONS  #
 #####################################
 */

PRIVATE int
rnaplot_SVG(const char         *filename,
            const char         *sequence,
            const char         *structure,
            vrna_plot_layout_t *layout,
            vrna_plot_data_t   *aux_data)
{
  short         *pt;
  unsigned int  i, j, n;
  float   size, *X, *Y, *CX, *CY, *R, xmin, xmax, ymin, ymax;
  FILE    *xyplot;
  double  *arccoordsSVG = NULL;

  CX = CY = R = NULL;

  n = strlen(sequence);

  if (layout->length != n) {
    vrna_log_error("structure and layout have different lengths");
    return 0;
  }

  if ((xyplot = fopen(filename, "w")) == NULL) {
    vrna_log_error("can't open svg file %s - not doing xy_plot", filename);
    return 0;
  }

  pt = vrna_ptable(structure);

  /* copy layout coords */
  X = (float *)vrna_alloc(sizeof(float) * n);
  Y = (float *)vrna_alloc(sizeof(float) * n);

  memcpy(X, layout->x, sizeof(float) * n);
  memcpy(Y, layout->y, sizeof(float) * n);

  if (layout->arcs)
    arccoordsSVG = transform_PS_arcs_to_SVG(n, layout->arcs);

  xmin = layout->bbox[0];
  ymin = layout->bbox[1];
  xmax = layout->bbox[2];
  ymax = layout->bbox[3];

  for (i = 0; i < n; i++)
    Y[i] = ymin + ymax - Y[i]; /* mirror coordinates so they look as in PS */

  size  = MAX2((xmax - xmin), (ymax - ymin));
  size  += 15; /* add some so the bounding box isn't too tight */

  float scale[2] = {
    SIZE / size, SIZE / size
  };

  float transform[2] = {
    (float)(size - xmin - xmax) / 2., (float)(size - ymin - ymax) / 2.
  };

  print_SVG_header(xyplot, scale, transform, 0);

  switch (layout->type) {
    case VRNA_PLOT_TYPE_CIRCULAR:
      {
        unsigned int radius  = 3 * n;
        unsigned int dr      = 0;

        R   = (float *)vrna_alloc((n + 1) * sizeof(float));
        CX  = (float *)vrna_alloc((n + 1) * sizeof(float));
        CY  = (float *)vrna_alloc((n + 1) * sizeof(float));

        for (i = 0; i < n; i++) {
          X[i]  -= radius;
          Y[i]  -= radius;
        }

        for (i = 0; i < n; i++) {
          if (i + 1 < pt[i + 1]) {
            dr = (pt[i + 1] - i + 1 <= (n / 2 + 1)) ? pt[i + 1] - i : i + n - pt[i + 1];
            R[i] = 1. - (2. * dr / (float)n);
          } else if (pt[i + 1]) {
            R[i] = R[pt[i + 1] - 1];
          } else {
            R[i] = 1.0;
          }

          CX[i] = X[i] * R[i] + radius;
          CY[i] = Y[i] * R[i] + radius;

          X[i]  += radius;
          Y[i]  += radius;
        }


      }
      print_SVG_backbone(xyplot, X, Y, n);
      break;

    case VRNA_PLOT_TYPE_PUZZLER:
    /* fall-through */
    case VRNA_PLOT_TYPE_TURTLE:
      {
        /*
         *  ################################### SVG for Puzzler and Turtle ##########################################
         *  Draw Backbone
         *  Idea: We look at a point and test if there exist an arc between previous point and the actual point.
         *  But beware: Most of the indexes of this part of the code are "half" magic numbers;they were made
         *  with brute force(trial and error) until it eventually worked.
         */
        short newLine = 0;    /* Flag if new Polyline should be created */
        fprintf(xyplot,
                "    <polyline  class=\"backbone\" id=\"outline\" points=\"\n");
        for (j = 1; j <= n; j++) {
          if (arccoordsSVG[2 * (j - 1)] < 0) {
            /* No arc(no valid radius) -> draw backbone */
            if (newLine) {
              newLine = 0;
              fprintf(xyplot,
                      "    <polyline  class=\"backbone\" id=\"outline%i\" points=\"\n",
                      j);
              fprintf(xyplot, "      %3.3f,%3.3f\n", X[j - 2], Y[j - 2]); /* If new line is created, it should include the previous previous point */
            }

            fprintf(xyplot, "      %3.3f,%3.3f\n", X[j - 1], Y[j - 1]);   /*  Then we include the previous point */
          } else {
            /*
             *  If there exist an arc, and there is no newline, then we are the last point
             *  in the polyline -> we set newLine to True and end Polyline
             */
            if (!newLine) {
              newLine = 1;
              fprintf(xyplot, "    \" />\n");
            }
          }
        }
        fprintf(xyplot, "    \" />\n");


        /* arcs */
        fprintf(xyplot, "    <g id=\"arcs\">\n");
        for (int j = 0; j < n - 1; j++) {
          if (arccoordsSVG[2 * (j + 1)] > 0) {
            /* If arc exists, then we draw the arc. Not much more to say */
            fprintf(xyplot,
                    "      <path class=\"backbone\" d=\"M %6.5f, %6.5f A %6.5f,%6.5f, %6.5f,%i, %i, %6.5f, %6.5f\" />\n",
                    X[j],
                    Y[j],
                    arccoordsSVG[2 * (j + 1)],
                    arccoordsSVG[2 * (j + 1)],
                    0.,
                    0,
                    (int)arccoordsSVG[2 * (j + 1) + 1],
                    X[j + 1],
                    Y[j + 1]);
          }
        }
        fprintf(xyplot, "    </g>\n");
      }
      break;

    default:
      print_SVG_backbone(xyplot, X, Y, n);
  }

  print_SVG_pairs(xyplot, pt, X, Y, CX, CY, n, layout->type);

  print_SVG_bases(xyplot, X, Y, sequence, n);

  print_SVG_footer_data(xyplot, n, sequence, structure, pt, X, Y);

  fclose(xyplot);

  free(pt);
  free(X);
  free(Y);
  free(arccoordsSVG);

  return 1; /* success */
}



static void
print_SVG_header(FILE         *fh,
                const float   scale[2],
                const float   transform[2],
                unsigned int  options)
{
  fprintf(fh,
          "%s",
          SVG_structure_plot_header);

  fprintf(fh,
          "  <g transform=\"scale(%7f,%7f) translate(%7f,%7f)\">\n",
          scale[0],
          scale[1],
          transform[0],
          transform[1]);
}


PRIVATE void
print_SVG_footer(FILE *fh)
{
  fprintf(fh,
          "  </g>\n%s",
          SVG_structure_plot_footer);
}


PRIVATE void
print_SVG_footer_data(FILE *fh,
                      unsigned int  n,
                      const char    *sequence,
                      const char    *structure,
                      const short   *pt,
                      float         *X,
                      float         *Y)
{
  unsigned int i;

  fprintf(fh,
    "  </g>\n"
    "<script type=\"text/ecmascript\">\n"
    "<![CDATA[\n"
    "  let sequence = \"%s\";\n"
    "  let structure = \"%s\";\n",
    sequence,
    structure);

  fprintf(fh, "  const basepairs = [\n");

  for (i = 1; i <= n; i++)
    if (pt[i] > i) {
      fprintf(fh, "    { i: %d, j: %d, type: \"cWW\" }", i, pt[i]);
      i++;
      break;
    }

  for (; i <= n; i++)
    if (pt[i] > i)
      fprintf(fh, ",\n    { i: %d, j: %d, type: \"cWW\" }", i, pt[i]);

  fprintf(fh, "];\n");

  fprintf(fh, "  const coords = [\n");
  fprintf(fh, "    { x: %6.3f, y: %6.3f }", X[0], Y[0]);
  for (i = 1; i < n; i++)
    fprintf(fh, ",\n    { x: %6.3f, y: %6.3f }", X[i], Y[i]);
  fprintf(fh, "];\n");

  fprintf(fh,
    "]]>\n"
    "</script>\n%s",
    SVG_structure_plot_footer);
}


PRIVATE void
print_SVG_bases(FILE          *fh,
                float         *X,
                float         *Y,
                const char    *string,
                unsigned int  n)
{
  unsigned int i;

  fprintf(fh,
          "    <g transform=\"translate(-4.6, 4)\" id=\"seq\">\n");

  for (i = 0; i < n; i++)
    fprintf(fh,
            "      <text class=\"nucleotide\" x=\"%.3f\" y=\"%.3f\">%c</text>\n",
            X[i], Y[i], string[i]);

  fprintf(fh,
          "    </g>\n");
}

PRIVATE void
print_SVG_backbone(FILE         *fh,
                   const float  *X,
                   const float  *Y,
                   unsigned int n)
{
  unsigned int i;

  fprintf(fh,
          "    <polyline class=\"backbone\" id=\"outline\" points=\"\n");

  for (i = 0; i < n; i++)
    fprintf(fh,
            "      %3.3f,%3.3f\n",
            X[i],
            Y[i]);

  fprintf(fh,
          "    \" />\n");
}


PRIVATE void
print_SVG_pairs(FILE          *fh,
                const short   *pt,
                const float   *X,
                const float   *Y,
                const float   *CX,
                const float   *CY,
                unsigned int  n,
                unsigned int  plot_type)
{
  unsigned int i, j;

  fprintf(fh,
          "    <g id=\"pairs\">\n");

  for (i = 1; i <= n; i++) {
    if ((j = (unsigned int)pt[i]) > i) {
      if (plot_type == VRNA_PLOT_TYPE_CIRCULAR) {
        fprintf(fh,
                "      <path class=\"basepairs\" id=\"%u,%u\" d=\"M %6.5f %6.5f C %6.5f,%6.5f %6.5f,%6.5f %6.5f %6.5f\" />\n",
                i,
                j,
                X[i - 1],
                Y[i - 1],
                CX[i - 1],
                CY[i - 1],
                CX[j - 1],
                CY[j - 1],
                X[j - 1],
                Y[j - 1]);
       } else {
        fprintf(fh,
                "      <line class=\"basepairs\" id=\"%u,%u\" x1=\"%6.5f\" y1=\"%6.5f\" x2=\"%6.5f\" y2=\"%6.5f\" />\n",
                i, j,
                X[i - 1], Y[i - 1],
                X[j - 1], Y[j - 1]);
      }
    }
  }
  fprintf(fh,
          "    </g>\n");
}

PRIVATE double *
transform_PS_arcs_to_SVG(unsigned int n,
                         double       *PSarcs)
{
  double *SVGarcs = (double*)vrna_alloc(n * 2 * sizeof(double));
  for (unsigned int i = 0; i < n; i++) {
    /* Arc exists */
    if (PSarcs[6 * i + 2] > 0){  /* radius */
      SVGarcs[2 * i]     = PSarcs[6 * i + 2];
      SVGarcs[2 * i + 1] = PSarcs[6 * i + 5]; /* goClockwise */
    } else {
      SVGarcs[2 * i] = -1;
      SVGarcs[2 * i + 1] = -1;
    }
  }

  return SVGarcs;
}

#ifndef VRNA_DISABLE_BACKWARD_COMPATIBILITY

/*
 * ###########################################
 * # deprecated functions below              #
 * ###########################################
 */

PUBLIC int
svg_rna_plot(char *string,
             char *structure,
             char *ssfile)
{
  int                 ret = 0;
  vrna_plot_layout_t  *layout;

  if ((string) &&
      (structure) &&
      (ssfile)) {
    layout  = vrna_plot_layout(structure, rna_plot_type);
    ret     = rnaplot_SVG((const char *)ssfile,
                          (const char *)string,
                          (const char *)structure,
                          layout,
                          NULL);

    vrna_plot_layout_free(layout);
  }

  return ret;
}

#endif
