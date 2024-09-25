/*
 *      GML output format for ViennaRNA
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
rnaplot_gml(const char          *filename,
            const char          *sequence,
            const char          *structure,
            vrna_plot_layout_t  *layout,
            vrna_plot_data_t    *aux_data,
            char                option);


/*
 #################################
 # BEGIN OF FUNCTION DEFINITIONS #
 #################################
 */
PUBLIC int
vrna_plot_structure_gml(const char          *filename,
                        const char          *sequence,
                        const char          *structure,
                        vrna_plot_layout_t  *layout,
                        vrna_plot_data_t    *data,
                        char                option)
{
  int                 ret = 0;
  vrna_plot_layout_t  *lt;

  if ((sequence) &&
      (structure) &&
      (filename)) {

    lt = layout;

    if (lt == NULL) /* use global default layout */
      lt = vrna_plot_layout(structure, VRNA_PLOT_TYPE_DEFAULT);

    ret = rnaplot_gml(filename, sequence, structure, lt, data, option);

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
rnaplot_gml(const char          *filename,
            const char          *sequence,
            const char          *structure,
            vrna_plot_layout_t  *layout,
            vrna_plot_data_t    *aux_data,
            char                option)
{
  FILE          *gmlfile;
  unsigned int  i, n;
  short         *pt;

  gmlfile = fopen(filename, "w");
  if (gmlfile == NULL) {
    vrna_log_error("can't open file %s - not doing xy_plot", filename);
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


  fprintf(gmlfile,
          "# Vienna RNA Package %s\n"
          "# GML Output\n"
          "# CreationDate: %s\n"
          "# Name: %s\n"
          "# Options: %s\n", VRNA_VERSION, vrna_time_stamp(), filename, vrna_md_option_string(NULL));
  fprintf(gmlfile,
          "graph [\n"
          " directed 0\n");
  for (i = 1; i <= n; i++) {
    fprintf(gmlfile,
            " node [ id %d ", i);
    if (option)
      fprintf(gmlfile,
              "label \"%c\"", sequence[i - 1]);

    if ((option == 'X') || (option == 'x'))
      fprintf(gmlfile,
              "\n  graphics [ x %9.4f y %9.4f ]\n", layout->x[i - 1], layout->y[i - 1]);

    fprintf(gmlfile, " ]\n");
  }
  for (i = 1; i < n; i++)
    fprintf(gmlfile,
            "edge [ source %d target %d ]\n", i, i + 1);
  for (i = 1; i <= n; i++) {
    if (pt[i] > i)
      fprintf(gmlfile,
              "edge [ source %d target %d ]\n", i, pt[i]);
  }
  fprintf(gmlfile, "]\n");
  fclose(gmlfile);

  free(pt);

  return 1; /* success */
}


#ifndef VRNA_DISABLE_BACKWARD_COMPATIBILITY

/*
 * ###########################################
 * # deprecated functions below              #
 *###########################################
 */


/* options for gml output:
 * uppercase letters: print sequence labels
 * lowercase letters: no sequence lables
 * graphics information:
 * x X  simple xy plot
 * (nothing else implemented at present)
 * default:           no graphics data at all
 */
PUBLIC int
gmlRNA(char *sequence,
       char *structure,
       char *filename,
       char option)
{
  int                 ret = 0;
  vrna_plot_layout_t  *layout;

  if ((sequence) &&
      (structure) &&
      (filename)) {

    layout  = vrna_plot_layout((const char *)structure,
                               rna_plot_type);

    ret     = rnaplot_gml((const char *)filename,
                          (const char *)sequence,
                          (const char *)structure,
                          layout,
                          NULL,
                          option);

    vrna_plot_layout_free(layout);
  }

  return ret;
}


#endif
