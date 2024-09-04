#ifndef VIENNA_RNA_PACKAGE_PLOT_STRUCTURE_H
#define VIENNA_RNA_PACKAGE_PLOT_STRUCTURE_H

#include <ViennaRNA/model.h>
#include <ViennaRNA/plotting/layouts.h>
#include "ViennaRNA/plotting/RNApuzzler/RNApuzzler.h"

#ifdef VRNA_WARN_DEPRECATED
# if defined(__clang__)
#  define DEPRECATED(func, msg) func __attribute__ ((deprecated("", msg)))
# elif defined(__GNUC__)
#  define DEPRECATED(func, msg) func __attribute__ ((deprecated(msg)))
# else
#  define DEPRECATED(func, msg) func
# endif
#else
# define DEPRECATED(func, msg) func
#endif

/**
 *  @file ViennaRNA/plotting/structures.h
 *  @ingroup   plotting_utils
 *  @brief Various functions for plotting RNA secondary structures
 */

/**
 *  @addtogroup   plotting_utils
 *  @{
 */

#define VRNA_FILE_FORMAT_EPS            0U
#define VRNA_FILE_FORMAT_SVG            1U
#define VRNA_FILE_FORMAT_GML            2U
#define VRNA_FILE_FORMAT_SSV            3U
#define VRNA_FILE_FORMAT_XRNA           4U
#define VRNA_FILE_FORMAT_PLOT_DEFAULT   VRNA_FILE_FORMAT_EPS


typedef struct vrna_plot_data_s vrna_plot_data_t;


struct vrna_plot_data_s {
  char          *pre;
  char          *post;
  vrna_md_t     *md;
  unsigned int  options;
};


int
vrna_plot_structure(const char          *filename,
                    const char          *sequence,
                    const char          *structure,
                    unsigned int        file_format,
                    vrna_plot_layout_t  *layout,
                    vrna_plot_data_t    *aux_data);


int
vrna_plot_structure_svg(const char          *filename,
                        const char          *sequence,
                        const char          *structure,
                        vrna_plot_layout_t  *layout,
                        vrna_plot_data_t    *data);


int
vrna_plot_structure_eps(const char          *filename,
                        const char          *sequence,
                        const char          *structure,
                        vrna_plot_layout_t  *layout,
                        vrna_plot_data_t    *data);


int
vrna_plot_structure_gml(const char          *filename,
                        const char          *sequence,
                        const char          *structure,
                        vrna_plot_layout_t  *layout,
                        vrna_plot_data_t    *data,
                        char                option);


int
vrna_plot_structure_ssv(const char          *filename,
                        const char          *sequence,
                        const char          *structure,
                        vrna_plot_layout_t  *layout,
                        vrna_plot_data_t    *data);


int
vrna_plot_structure_xrna(const char          *filename,
                         const char          *sequence,
                         const char          *structure,
                         vrna_plot_layout_t  *layout,
                         vrna_plot_data_t    *data);

/**
 *  @brief Produce a secondary structure graph in PostScript and write it to 'filename'.
 *
 *  Note that this function has changed from previous versions
 *  and now expects the structure to be plotted in dot-bracket notation as an
 *  argument. It does not make use of the global #base_pair array anymore.
 *
 *  @param seq        The RNA sequence
 *  @param structure  The secondary structure in dot-bracket notation
 *  @param file       The filename of the postscript output
 *  @param md_p       Model parameters used to generate a commandline option string in the output (Maybe NULL)
 *  @return           1 on success, 0 otherwise
 */
int
vrna_file_PS_rnaplot(const char *seq,
                     const char *structure,
                     const char *file,
                     vrna_md_t  *md_p);


/**
 *  @brief Produce a secondary structure graph in PostScript including additional
 *  annotation macros and write it to 'filename'
 *
 *  Same as vrna_file_PS_rnaplot() but adds extra PostScript macros for various
 *  annotations (see generated PS code). The 'pre' and 'post'
 *  variables contain PostScript code that is verbatim copied in the
 *  resulting PS file just before and after the structure plot.
 *  If both arguments ('pre' and 'post') are NULL, no additional macros will
 *  be printed into the PostScript.
 *
 *  @param seq     The RNA sequence
 *  @param structure  The secondary structure in dot-bracket notation
 *  @param file       The filename of the postscript output
 *  @param pre        PostScript code to appear before the secondary structure plot
 *  @param post       PostScript code to appear after the secondary structure plot
 *  @param md_p       Model parameters used to generate a commandline option string in the output (Maybe NULL)
 *  @return           1 on success, 0 otherwise
 */
int vrna_file_PS_rnaplot_a( const char      *seq,
                            const char      *structure,
                            const char      *file,
                            const char      *pre,
                            const char      *post,
                            vrna_md_t       *md_p);


int
vrna_file_PS_rnaplot_layout(const char          *seq,
                            const char          *structure,
                            const char          *ssfile,
                            const char          *pre,
                            const char          *post,
                            vrna_md_t           *md_p,
                            vrna_plot_layout_t  *layout);

#ifndef VRNA_DISABLE_BACKWARD_COMPATIBILITY

/* write PostScript drawing of structure to file with annotation */
int
PS_rna_plot_snoop_a(const char  *string,
                    const char  *structure,
                    const char  *ssfile,
                    int         *relative_access,
                    const char  *seqs[]);


/**
 *  @brief Produce a secondary structure graph in Graph Meta Language (gml) and write it to a file
 *
 *  If 'option' is an uppercase letter the RNA sequence is used to label nodes, if 'option' equals
 *  @a 'X' or @a 'x' the resulting file will coordinates for an initial layout of the graph.
 *
 *  @param  string    The RNA sequence
 *  @param  structure The secondary structure in dot-bracket notation
 *  @param  ssfile    The filename of the gml output
 *  @param  option    The option flag
 *  @return           1 on success, 0 otherwise
 */
int
gmlRNA(char *string,
       char *structure,
       char *ssfile,
       char option);


/**
 *  @brief  Produce a secondary structure graph in SStructView format
 *
 *  Write coord file for SStructView
 *
 *  @param  string    The RNA sequence
 *  @param  structure The secondary structure in dot-bracket notation
 *  @param  ssfile    The filename of the ssv output
 *  @return           1 on success, 0 otherwise
 */
int
ssv_rna_plot(char *string,
             char *structure,
             char *ssfile);


/**
 *  @brief Produce a secondary structure plot in SVG format and write it to a file
 *
 *  @param string     The RNA sequence
 *  @param structure  The secondary structure in dot-bracket notation
 *  @param ssfile     The filename of the svg output
 *  @return           1 on success, 0 otherwise
 */
int
svg_rna_plot(char *string,
             char *structure,
             char *ssfile);


/**
 *  @brief Produce a secondary structure plot for further editing in XRNA
 *
 *  @param string     The RNA sequence
 *  @param structure  The secondary structure in dot-bracket notation
 *  @param ssfile     The filename of the xrna output
 *  @return           1 on success, 0 otherwise
 */
int
xrna_plot(char  *string,
          char  *structure,
          char  *ssfile);


/**
 *  @brief Produce a secondary structure graph in PostScript and write it to 'filename'.
 *
 *  @deprecated   Use vrna_file_PS_rnaplot() instead!
 */
DEPRECATED(int PS_rna_plot(char *string,
                           char *structure,
                           char *file),
           "Use vrna_file_PS_rnaplot() instead");

/**
 *  @brief Produce a secondary structure graph in PostScript including additional
 *  annotation macros and write it to 'filename'
 *
 *  @deprecated   Use vrna_file_PS_rnaplot_a() instead!
 */
DEPRECATED(int PS_rna_plot_a(char *string,
                             char *structure,
                             char *file,
                             char *pre,
                             char *post),
           "Use vrna_file_PS_rnaplot_a() instead");

/**
 *  @brief Produce a secondary structure graph in PostScript including additional
 *  annotation macros and write it to 'filename' (detect and draw g-quadruplexes)
 *
 *  @deprecated   Use vrna_file_PS_rnaplot_a() instead!
 */
DEPRECATED(int PS_rna_plot_a_gquad(char *string,
                                   char *structure,
                                   char *ssfile,
                                   char *pre,
                                   char *post),
           "Use vrna_file_PS_rnaplot_a() instead");

#endif

/**
 * @}
 */

#endif
