/*
 *      XRNA output format for ViennaRNA
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
rnaplot_xrna(const char          *filename,
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
vrna_plot_structure_xrna(const char          *filename,
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

    ret = rnaplot_xrna(filename, sequence, structure, lt, data);

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
rnaplot_xrna(const char          *filename,
             const char          *sequence,
             const char          *structure,
             vrna_plot_layout_t  *layout,
             vrna_plot_data_t    *aux_data)
{
  /* produce input for XRNA RNA drawing program */
  FILE          *ss_file;
  unsigned int  i, n;
  short         *pt;

  ss_file = fopen(filename, "w");
  if (ss_file == NULL) {
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

  fprintf(ss_file,
          "# Vienna RNA Package %s, XRNA output\n"
          "# CreationDate: %s\n"
          "# Options: %s\n", VRNA_VERSION, vrna_time_stamp(), vrna_md_option_string(NULL));

  for (i = 1; i <= n; i++)
    /* XRNA likes to have coordinate mirrored, so we use (-X, Y) */
    fprintf(ss_file,
            "%d %c %6.2f %6.2f %d %d\n",
            i,
            sequence[i - 1],
            -layout->x[i - 1],
            layout->y[i - 1],
            (pt[i] != 0) ? 1 : 0,
            pt[i]);

  fclose(ss_file);

  free(pt);

  return 1; /* success */
}


#ifndef VRNA_DISABLE_BACKWARD_COMPATIBILITY

/*
 * ###########################################
 * # deprecated functions below              #
 *###########################################
 */


PUBLIC int
xrna_plot(char *sequence,
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

    ret     = rnaplot_xrna((const char *)filename,
                           (const char *)sequence,
                           (const char *)structure,
                           layout,
                           NULL);

    vrna_plot_layout_free(layout);
  }

  return ret;
}

#endif
