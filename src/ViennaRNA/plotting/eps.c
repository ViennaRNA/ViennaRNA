/*
 *      PostScript and other output formats for RNA secondary structure plots
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
#include "ViennaRNA/utils/basic.h"
#include "ViennaRNA/sequences/alignments.h"
#include "ViennaRNA/utils/log.h"
#include "ViennaRNA/sequences/alphabet.h"
#include "ViennaRNA/model.h"
#include "ViennaRNA/fold_vars.h"
#include "ViennaRNA/structures/dotbracket.h"
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
rnaplot_EPS(const char          *seq,
            const char          *structure,
            const char          *ssfile,
            const char          *pre,
            const char          *post,
            vrna_md_t           *md_p,
            vrna_plot_layout_t  *layout);


/*
 #################################
 # BEGIN OF FUNCTION DEFINITIONS #
 #################################
 */

PUBLIC int
vrna_plot_structure_eps(const char          *filename,
                        const char          *sequence,
                        const char          *structure,
                        vrna_plot_layout_t  *layout,
                        vrna_plot_data_t    *data)
{
  const char          *pre, *post;
  int                 ret = 0;
  vrna_md_t           *md;
  vrna_plot_layout_t  *lt;

  if ((sequence) &&
      (structure) &&
      (filename)) {

    pre = post = NULL;
    md  = NULL;
    lt  = layout;

    if (lt == NULL) /* use global default layout */
      lt = vrna_plot_layout(structure, VRNA_PLOT_TYPE_DEFAULT);

    if (data) {
      md = data->md;
      pre = data->pre;
      post = data->post;
    }

    ret = rnaplot_EPS(sequence,
                      structure,
                      filename,
                      pre,
                      post,
                      md,
                      lt);

    if (lt != layout)
      vrna_plot_layout_free(lt);
  }

  return ret;
}


PUBLIC int
vrna_file_PS_rnaplot(const char *string,
                     const char *structure,
                     const char *ssfile,
                     vrna_md_t  *md_p)
{
  return vrna_file_PS_rnaplot_a(string, structure, ssfile, NULL, NULL, md_p);
}


PUBLIC int
vrna_file_PS_rnaplot_a(const char *seq,
                       const char *structure,
                       const char *ssfile,
                       const char *pre,
                       const char *post,
                       vrna_md_t  *md_p)
{
  int                 ret;
  vrna_plot_layout_t  *layout;

  layout  = vrna_plot_layout(structure, rna_plot_type);
  ret     = vrna_file_PS_rnaplot_layout(seq,
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
vrna_file_PS_rnaplot_layout(const char          *seq,
                            const char          *structure,
                            const char          *ssfile,
                            const char          *pre,
                            const char          *post,
                            vrna_md_t           *md_p,
                            vrna_plot_layout_t  *layout)
{
  if (!ssfile) {
    vrna_log_warning("Filename missing!");
  } else if (!seq) {
    vrna_log_warning("Sequence missing");
  } else if (!structure) {
    vrna_log_warning("Structure missing");
  } else if (!layout) {
    vrna_log_warning("Layout missing");
  } else if ((strlen(seq) != strlen(structure)) ||
             (strlen(structure) != layout->length)) {
    vrna_log_warning("Sequence, structure, and coordinate layout have different lengths! "
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
rnaplot_EPS(const char          *seq,
            const char          *structure,
            const char          *ssfile,
            const char          *pre,
            const char          *post,
            vrna_md_t           *md_p,
            vrna_plot_layout_t  *layout)
{
  float     xmin, xmax, ymin, ymax;
  int       i;
  unsigned int  ee, Lg, l[3], length;
  int       gb, ge, bbox[4];
  float     *X, *Y;
  FILE      *xyplot;
  short     *pair_table;
  char      *c, *string;
  vrna_md_t md;

  if (!md_p) {
    set_model_details(&md);
    md_p = &md;
  }

  xyplot = fopen(ssfile, "w");
  if (xyplot == NULL) {
    vrna_log_error("can't open file %s - not doing xy_plot", ssfile);
    return 0;
  }

  string  = strdup(seq);
  length  = strlen(string);

  pair_table = vrna_ptable(structure);

  bbox[0] = bbox[1] = 0;
  bbox[2] = bbox[3] = 700;

  /* increase bounding box to negative y-axis to allow for legends */
  if ((pre) || (post))
    bbox[1] = -140;

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
    cutpoint          = c - structure;
    string[cutpoint]  = ' '; /* replace & with space */
    fprintf(xyplot, "/cutpoint %d def\n", cutpoint);
  }

  /* sequence */
  print_PS_sequence(xyplot, string);

  /* coordinates */
  print_PS_coords(xyplot, layout->x, layout->y, length);

  fprintf(xyplot, "/arcs [\n");
  if (layout->arcs) {
    for (i = 0; i < length; i++) {
      if (layout->arcs[6 * i + 2] > 0) {
        /* 6*i+2 is the radius parameter */
        fprintf(xyplot, "[%3.8f %3.8f %3.8f %3.8f %3.8f %3.8f]\n",
                layout->arcs[6 * i + 0],
                layout->arcs[6 * i + 1],
                layout->arcs[6 * i + 2],
                layout->arcs[6 * i + 3],
                layout->arcs[6 * i + 4],
                layout->arcs[6 * i + 5]
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
  if (rna_plot_type == VRNA_PLOT_TYPE_CIRCULAR)
    fprintf(xyplot, "/cpr %6.2f def\n", (float)3 * length);

  /* base pairs */
  fprintf(xyplot, "/pairs [\n");
  for (i = 1; i <= length; i++)
    if (pair_table[i] > i)
      fprintf(xyplot, "[%d %d]\n", i, pair_table[i]);

  /* add gquad pairs */
  ge = 0;
  while ((ee = vrna_gq_parse(structure + ge, &Lg, l)) > 0) {
    int k;
    fprintf(xyplot, "%% gquad\n");
    ge  += ee;
    if (4 * Lg + l[0] + l[1] + l[2] > ee) {
      gb = length + ge - Lg * 4 - l[0] - l[1] - l[2] + 1; /* add pseudo-base pair encloding gquad */
    } else {
      gb  = ge - Lg * 4 - l[0] - l[1] - l[2] + 1; /* add pseudo-base pair encloding gquad */
    }

    for (k = 0; k < Lg; k++) {
      int ii, jj, il;
      for (il = 0, ii = gb + k; il < 3; il++) {
        jj = ii + l[il] + Lg;
        fprintf(xyplot, "[%d %d]\n", (ii - 1) % (length) + 1, (jj - 1) % (length) + 1);
        ii = jj;
      }
      jj = gb + k;
      fprintf(xyplot, "[%d %d]\n", (jj - 1) % (length) + 1, (ii - 1) % (length) + 1);
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


PUBLIC int
PS_rna_plot_snoop_a(const char  *string,
                    const char  *structure,
                    const char  *ssfile,
                    int         *relative_access,
                    const char  *seqs[])
{
  int       i, length, bbox[4];
  float     *X, *Y;
  FILE      *xyplot;
  short     *pair_table;
  short     *pair_table_snoop;
  vrna_md_t md;

  set_model_details(&md);

  length = strlen(string);

  xyplot = fopen(ssfile, "w");
  if (xyplot == NULL) {
    vrna_log_warning("can't open file %s - not doing xy_plot", ssfile);
    return 0;
  }

  pair_table        = vrna_ptable(structure);
  pair_table_snoop  = vrna_pt_snoop_get(structure);

  i = vrna_plot_coords_pt(pair_table, &X, &Y, rna_plot_type);

  if (i != length)
    vrna_log_warning("strange things happening in PS_rna_plot...");

  /*   printf("cut_point %d\n", cut_point); */

  /*
   *   for (i = 1; i < length; i++) {
   *     printf("%d X %f Y %f \n", i, X[i], Y[i]);
   *     xmin = X[i] < xmin ? X[i] : xmin;
   *     xmax = X[i] > xmax ? X[i] : xmax;
   *     ymin = Y[i] < ymin ? Y[i] : ymin;
   *     ymax = Y[i] > ymax ? Y[i] : ymax;
   *   }
   * localize centre of the interaction bucket. Geometry
   */

  for (i = 1; i < cut_point; i++) {
    /* internal loop of size 0 */
    if (pair_table_snoop[i] != 0) {
      X[i - 1]  = X[pair_table_snoop[i] - 1];
      Y[i - 1]  = Y[pair_table_snoop[i] - 1];
    } else if (pair_table_snoop[i - 1] && pair_table_snoop[i + 1]) {
      /* internal loop of size 1 */
      X[i - 1]  = X[pair_table_snoop[i - 1] - 1 - 1];
      Y[i - 1]  = Y[pair_table_snoop[i - 1] - 1 - 1];
    } else if (pair_table_snoop[i - 1] && pair_table_snoop[i + 2]) {
      /* internal loop of size 2 */
      if (pair_table_snoop[i - 1] - pair_table_snoop[i + 2] == 2) {
        X[i - 1]  = X[pair_table_snoop[i - 1] - 2];
        Y[i - 1]  = Y[pair_table_snoop[i - 1] - 2];
        X[i]      = X[pair_table_snoop[i + 2]];
        Y[i]      = Y[pair_table_snoop[i + 2]];
        i++;
      } else if (pair_table[pair_table_snoop[i - 1] - 1]) {
        X[i - 1]  = X[pair_table_snoop[i - 1] - 2];
        Y[i - 1]  = Y[pair_table_snoop[i - 1] - 2];
        X[i]      = X[pair_table[pair_table_snoop[i - 1] - 1] - 1];
        Y[i]      = Y[pair_table[pair_table_snoop[i - 1] - 1] - 1];
        i++;
      } else if (pair_table[pair_table_snoop[i - 1] - 2]) {
        X[i - 1]  = X[pair_table_snoop[i - 1] - 3];
        Y[i - 1]  = Y[pair_table_snoop[i - 1] - 3];
        X[i]      = X[pair_table[pair_table_snoop[i - 1] - 2] - 1];
        Y[i]      = Y[pair_table[pair_table_snoop[i - 1] - 2] - 1];
        i++;
      } else if (pair_table[pair_table_snoop[i - 1] - 3]) {
        X[i - 1]  = X[pair_table_snoop[i - 1] - 4];
        Y[i - 1]  = Y[pair_table_snoop[i - 1] - 4];
        X[i]      = X[pair_table[pair_table_snoop[i - 1] - 3] - 1];
        Y[i]      = Y[pair_table[pair_table_snoop[i - 1] - 3] - 1];
        i++;
      } else {
        X[i - 1]  = X[pair_table_snoop[i - 1] - 2];
        Y[i - 1]  = Y[pair_table_snoop[i - 1] - 2];
        X[i]      = X[pair_table_snoop[i + 2]];
        Y[i]      = Y[pair_table_snoop[i + 2]];
        i++;
      }
    } else if (pair_table_snoop[i - 1] && pair_table_snoop[i + 3]) {
      /* internal loop of size 2 */
      if (pair_table[pair_table_snoop[i - 1] - 1]) {
        X[i - 1]  = 0.5 * (X[pair_table_snoop[i - 1] - 1] + X[pair_table_snoop[i - 1] - 2]);
        Y[i - 1]  = 0.5 * (Y[pair_table_snoop[i - 1] - 1] + Y[pair_table_snoop[i - 1] - 2]);
        X[i]      = 0.5 *
                    (X[pair_table[pair_table_snoop[i - 1] - 1] - 1] +
                     X[pair_table_snoop[i - 1] - 2]);
        Y[i] = 0.5 *
               (Y[pair_table[pair_table_snoop[i - 1] - 1] - 1] + Y[pair_table_snoop[i - 1] - 2]);
        X[i + 1] = 0.5 *
                   (X[pair_table[pair_table_snoop[i - 1] - 1] - 2] +
                    X[pair_table[pair_table_snoop[i - 1] - 1] - 1]);
        Y[i + 1] = 0.5 *
                   (Y[pair_table[pair_table_snoop[i - 1] - 1] - 2] +
                    Y[pair_table[pair_table_snoop[i - 1] - 1] - 1]);
        i++;
        i++;
      } else if (pair_table[pair_table_snoop[i - 1] - 2]) {
        X[i - 1]  = 0.5 * (X[pair_table_snoop[i - 1] - 2] + X[pair_table_snoop[i - 1] - 3]);
        Y[i - 1]  = 0.5 * (Y[pair_table_snoop[i - 1] - 2] + Y[pair_table_snoop[i - 1] - 3]);
        X[i]      = 0.5 *
                    (X[pair_table[pair_table_snoop[i - 1] - 2] - 1] +
                     X[pair_table_snoop[i - 1] - 3]);
        Y[i] = 0.5 *
               (Y[pair_table[pair_table_snoop[i - 1] - 2] - 1] + Y[pair_table_snoop[i - 1] - 3]);
        X[i + 1] = 0.5 *
                   (X[pair_table[pair_table_snoop[i - 1] - 2] - 2] +
                    X[pair_table[pair_table_snoop[i - 1] - 2] - 1]);
        Y[i + 1] = 0.5 *
                   (Y[pair_table[pair_table_snoop[i - 1] - 2] - 2] +
                    Y[pair_table[pair_table_snoop[i - 1] - 2] - 1]);
        i++;
        i++;
      } else if (pair_table[pair_table_snoop[i - 1] - 3]) {
        X[i - 1]  = 0.5 * (X[pair_table_snoop[i - 1] - 3] + X[pair_table_snoop[i - 1] - 4]);
        Y[i - 1]  = 0.5 * (Y[pair_table_snoop[i - 1] - 3] + Y[pair_table_snoop[i - 1] - 4]);
        X[i]      = 0.5 *
                    (X[pair_table[pair_table_snoop[i - 1] - 3] - 1] +
                     X[pair_table_snoop[i - 1] - 4]);
        Y[i] = 0.5 *
               (Y[pair_table[pair_table_snoop[i - 1] - 3] - 1] + Y[pair_table_snoop[i - 1] - 4]);
        X[i + 1] = 0.5 *
                   (X[pair_table[pair_table_snoop[i - 1] - 3] - 2] +
                    X[pair_table[pair_table_snoop[i - 1] - 3] - 1]);
        Y[i + 1] = 0.5 *
                   (Y[pair_table[pair_table_snoop[i - 1] - 3] - 2] +
                    Y[pair_table[pair_table_snoop[i - 1] - 3] - 1]);
        i++;
        i++;
      } else {
        X[i - 1]  = X[pair_table_snoop[i - 1] - 2];
        Y[i - 1]  = Y[pair_table_snoop[i - 1] - 2];
        X[i]      = X[pair_table_snoop[i - 1] - 2];
        Y[i]      = Y[pair_table_snoop[i - 1] - 2];
        X[i + 1]  = X[pair_table_snoop[i - 1] - 2];
        Y[i + 1]  = Y[pair_table_snoop[i - 1] - 2];
        i++;
        i++;
      }
    }
  }
  double  xC;
  double  yC;
  float   X0 = -1, Y0 = -1, X1 = -1, Y1 = -1, X2 = -1, Y2 = -1;

  /*   int c1,c2,c3; */
  for (i = 1; i < cut_point; i++) {
    if (pair_table_snoop[i]) {
      X0  = X[pair_table_snoop[i] - 1];
      Y0  = Y[pair_table_snoop[i] - 1];
      /*     c1=pair_table_snoop[i]; */
      i++;
      break;
    }
  }
  for (; i < cut_point; i++) {
    if (pair_table_snoop[i]) {
      X1  = X[pair_table_snoop[i] - 1];
      Y1  = Y[pair_table_snoop[i] - 1];
      /*   c2=pair_table_snoop[i]; */
      i++;
      break;
    }
  }
  for (; i < cut_point; i++) {
    if (pair_table_snoop[i]) {
      X2  = X[pair_table_snoop[i] - 1];
      Y2  = Y[pair_table_snoop[i] - 1];
      /*   c3=pair_table_snoop[i]; */
      i++;
      break;
    }
  }
  /*
   *   for(i=cut_point-2;i>pair_table_snoop[c1]; i--){
   *     if(pair_table_snoop[i]){
   *       X1=X[pair_table_snoop[i]-1];Y1=Y[pair_table_snoop[i]-1];
   *       c2=pair_table_snoop[i];
   *       i++;
   *       break;
   *     }
   *   }
   *   for(i=pair_table_snoop[c1]+1;i<pair_table_snoop[c2]; i++){
   *     if(pair_table_snoop[i]){
   *       X2=X[pair_table_snoop[i]-1];Y2=Y[pair_table_snoop[i]-1];
   *       c3=pair_table_snoop[i];
   *       i++;
   *       break;
   *     }
   *   }
   */
  if (X0 < 0 || X1 < 0 || X2 < 0) {
    printf("Could not get the center of the binding bucket. No ps file will be produced!\n");
    fclose(xyplot);
    free(pair_table);
    free(pair_table_snoop);
    free(X);
    free(Y);
    pair_table        = NULL;
    pair_table_snoop  = NULL;
    X                 = NULL;
    Y                 = NULL;
    return 0;
  }

  double  alpha   = (X0 - X1) / (Y1 - Y0);
  double  alpha_p = (X1 - X2) / (Y2 - Y1);
  double  b       = (Y0 + Y1 - alpha * (X0 + X1)) * 0.5;
  double  b_p     = (Y1 + Y2 - alpha_p * (X1 + X2)) * 0.5;

  /*    if(abs(alpha -alpha_p) > 0.0000001){ */
  xC  = (b_p - b) / (alpha - alpha_p);
  yC  = alpha * xC + b;
  for (i = 1; i < cut_point; i++) {
    X[i - 1]  = X[i - 1] + 0.25 * (xC - X[i - 1]);
    Y[i - 1]  = Y[i - 1] + 0.25 * (yC - Y[i - 1]);
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

  if (seqs) {
    A = vrna_annotate_covar_db_extended((const char **)seqs,
                                        structure,
                                        &md,
                                        VRNA_BRACKETS_RND | VRNA_BRACKETS_ANG | VRNA_BRACKETS_SQR);
  }

  fprintf(xyplot, "%% data start here\n");
  /* cut_point */
  if (cut_point > 0 && cut_point <= strlen(string))
    fprintf(xyplot, "/cutpoint %d def\n", cut_point - 1);

  /* sequence */
  print_PS_sequence(xyplot, string);

  /* coordinates */
  print_PS_coords(xyplot, X, Y, length);

  /* base pairs */
  fprintf(xyplot, "/pairs [\n");
  for (i = 1; i <= length; i++)
    if (pair_table[i] > i)
      fprintf(xyplot, "[%d %d]\n", i, pair_table[i]);

  for (i = 1; i <= length; i++)
    if (pair_table_snoop[i] > i)
      fprintf(xyplot, "[%d %d]\n", i, pair_table_snoop[i]);

  fprintf(xyplot, "] def\n\n");
  if (relative_access) {
    fprintf(xyplot, "/S [\n");
    for (i = 0; i < cut_point - 1; i++)
      fprintf(xyplot, " %f\n", (float)relative_access[i] / 100);
    fprintf(xyplot, "]\n bind def\n");
    fprintf(xyplot, "/invert false def\n");
    fprintf(xyplot, "/range 0.8 def\n");
    fprintf(xyplot, "/drawreliability {\n"
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

  fprintf(xyplot, "%%switch off outline pairs or bases by removing these lines\n");
  if (relative_access)
    fprintf(xyplot, "drawreliability\n");

  fprintf(xyplot,
          "drawoutline\n"
          "drawpairs\n"
          "drawbases\n");
  /*
   * fprintf(xyplot, "%d cmark\n",c1);
   * fprintf(xyplot, "%d cmark\n",c2);
   * fprintf(xyplot, "%d cmark\n",c3);
   */
  if (seqs) {
    fprintf(xyplot, "%% Start Annotations\n");
    fprintf(xyplot, "%s\n", A[1]);
    fprintf(xyplot, "%% End Annotations\n");
  }

  print_PS_footer(xyplot);

  fclose(xyplot);
  if (seqs) {
    free(A[0]);
    free(A[1]);
    free(A);
  }

  free(pair_table);
  free(pair_table_snoop);
  free(X);
  free(Y);
  return 1; /* success */
}


#ifndef VRNA_DISABLE_BACKWARD_COMPATIBILITY

/*
 * ###########################################
 * # deprecated functions below              #
 *###########################################
 */
PUBLIC int
PS_rna_plot(char  *string,
            char  *structure,
            char  *ssfile)
{
  return vrna_file_PS_rnaplot((const char *)string,
                              (const char *)structure,
                              (const char *)ssfile,
                              NULL);
}


PUBLIC int
PS_rna_plot_a(char  *string,
              char  *structure,
              char  *ssfile,
              char  *pre,
              char  *post)
{
  return vrna_file_PS_rnaplot_a((const char *)string,
                                (const char *)structure,
                                (const char *)ssfile,
                                (const char *)pre,
                                (const char *)post,
                                NULL);
}


PUBLIC int
PS_rna_plot_a_gquad(char  *string,
                    char  *structure,
                    char  *ssfile,
                    char  *pre,
                    char  *post)
{
  return vrna_file_PS_rnaplot_a((const char *)string,
                                (const char *)structure,
                                (const char *)ssfile,
                                (const char *)pre,
                                (const char *)post,
                                NULL);
}


#endif
