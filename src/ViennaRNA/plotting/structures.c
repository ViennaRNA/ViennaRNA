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
vrna_plot_structure(const char          *filename,
                    const char          *sequence,
                    const char          *structure,
                    unsigned int        file_format,
                    vrna_plot_layout_t  *layout,
                    vrna_plot_data_t    *aux_data)
{
  int default_layout  = 0;
  int no_aux_data     = 0;
  int ret             = 0;

  if ((structure) &&
      (filename)) {

    if (layout == NULL) { /* use global default layout */
      default_layout = 1;
      layout = vrna_plot_layout(structure, VRNA_PLOT_TYPE_DEFAULT);
    }

    switch (file_format) {
      case VRNA_FILE_FORMAT_SVG:
        ret = vrna_plot_structure_svg(filename,
                                      sequence,
                                      structure,
                                      layout,
                                      aux_data);
        break;

      case VRNA_FILE_FORMAT_GML:
        ret = vrna_plot_structure_gml(filename,
                                      sequence,
                                      structure,
                                      layout,
                                      aux_data,
                                      'x');
        break;

      case VRNA_FILE_FORMAT_SSV:
        ret = vrna_plot_structure_ssv(filename,
                                      sequence,
                                      structure,
                                      layout,
                                      aux_data);
        break;

      case VRNA_FILE_FORMAT_XRNA:
        ret = vrna_plot_structure_xrna(filename,
                                       sequence,
                                       structure,
                                       layout,
                                       aux_data);
        break;

      case VRNA_FILE_FORMAT_EPS:
        /* fall through */

      default: /* EPS output */
        if (aux_data) {
          ret = vrna_file_PS_rnaplot_layout(sequence,
                                            structure,
                                            filename,
                                            aux_data->pre,
                                            aux_data->post,
                                            aux_data->md,
                                            layout);
        } else {
          ret = vrna_file_PS_rnaplot_layout(sequence,
                                            structure,
                                            filename,
                                            NULL,
                                            NULL,
                                            NULL,
                                            layout);
        }
        break;   
    }
        
    if (default_layout)
      vrna_plot_layout_free(layout);
  }

  return ret;
}

/*
 #####################################
 # BEGIN OF STATIC HELPER FUNCTIONS  #
 #####################################
 */
