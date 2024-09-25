/*
 *      SSV output format for ViennaRNA
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
PRIVATE int
rnaplot_ssv(const char          *filename,
            const char          *sequence,
            const char          *structure,
            vrna_plot_layout_t  *layout,
            vrna_plot_data_t    *aux_data);


/*
 #################################
 # BEGIN OF FUNCTION DEFINITIONS #
 #################################
 */
PUBLIC int
vrna_plot_structure_ssv(const char          *filename,
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

    ret = rnaplot_ssv(filename, sequence, structure, lt, data);

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
rnaplot_ssv(const char          *filename,
            const char          *sequence,
            const char          *structure,
            vrna_plot_layout_t  *layout,
            vrna_plot_data_t    *aux_data)
{
  /* produce input for the SStructView java applet */
  FILE          *ssvfile;
  unsigned int  i, n, bp;
  short         *pt;
  float         *X, *Y;
  float         xmin, xmax, ymin, ymax;

  ssvfile = fopen(filename, "w");
  if (ssvfile == NULL) {
    vrna_log_warning("can't open file %s - not doing xy_plot", filename);
    return 0;
  }

  n = strlen(sequence);

  if (n != strlen(structure)) {
    vrna_log_warning("Sequence and Structure have different lengths (%u vs. %u)",
                     n, strlen(structure));
    return 0;
  }

  if (n != layout->length) {
    vrna_log_warning("Structure and Layout have different lengths (%u vs. %u)",
                     n, layout->length);
    return 0;
  }

  pt  = vrna_ptable(structure);

  /* copy layout coords */
  X = (float *)vrna_alloc(sizeof(float) * n);
  Y = (float *)vrna_alloc(sizeof(float) * n);

  memcpy(X, layout->x, sizeof(float) * n);
  memcpy(Y, layout->y, sizeof(float) * n);

  xmin = layout->bbox[0];
  ymin = layout->bbox[1];
  xmax = layout->bbox[2];
  ymax = layout->bbox[3];

  if (xmin < 1) {
    for (i = 0; i <= n; i++)
      X[i] -= xmin - 1;
    xmin = 1;
  }

  if (ymin < 1) {
    for (i = 0; i <= n; i++)
      Y[i] -= ymin - 1;
    ymin = 1;
  }

#if 0
  {
    float size, xoff, yoff;
    float JSIZE = 500; /* size of the java applet window */
    /* rescale coordinates, center on square of size HSIZE */
    size  = MAX2((xmax - xmin), (ymax - ymin));
    xoff  = (size - xmax + xmin) / 2;
    yoff  = (size - ymax + ymin) / 2;
    for (i = 0; i <= length; i++) {
      X[i]  = (X[i] - xmin + xoff) * (JSIZE - 10) / size + 5;
      Y[i]  = (Y[i] - ymin + yoff) * (JSIZE - 10) / size + 5;
    }
  }
#endif
  /* */

  fprintf(ssvfile,
          "# Vienna RNA Package %s\n"
          "# SStructView Output\n"
          "# CreationDate: %s\n"
          "# Name: %s\n"
          "# Options: %s\n", VRNA_VERSION, vrna_time_stamp(), filename, vrna_md_option_string(NULL));

  for (i = 1; i <= n; i++)
    fprintf(ssvfile, "BASE\t%d\t%c\t%d\t%d\n",
            i, sequence[i - 1], (int)(X[i - 1] + 0.5), (int)(Y[i - 1] + 0.5));

  for (bp = 1, i = 1; i <= n; i++)
    if (pt[i] > i)
      fprintf(ssvfile, "BASE-PAIR\tbp%d\t%d\t%d\n", bp++, i, pt[i]);

  fclose(ssvfile);

  free(pt);
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


/*--------------------------------------------------------------------------*/

PUBLIC int
ssv_rna_plot(char *sequence,
             char *structure,
             char *filename)
{
  int                 ret = 0;
  vrna_plot_layout_t  *layout;

  if ((sequence) &&
      (structure) &&
      (filename)) {

    layout  = vrna_plot_layout((const char *)structure,
                               rna_plot_type);

    ret     = rnaplot_ssv((const char *)filename,
                          (const char *)sequence,
                          (const char *)structure,
                          layout,
                          NULL);

    vrna_plot_layout_free(layout);
  }

  return ret;
}

#endif
