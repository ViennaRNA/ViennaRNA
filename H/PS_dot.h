#ifndef __VIENNA_RNA_PACKAGE_PS_DOT_H__
#define __VIENNA_RNA_PACKAGE_PS_DOT_H__

#include "data_structures.h"
#include "plot_layouts.h"

#ifdef __GNUC__
#define DEPRECATED(func) func __attribute__ ((deprecated))
#else
#define DEPRECATED(func) func
#endif

/**
 *  \file PS_dot.h
 *  \brief Various functions for plotting RNA secondary structures, dot-plots and other
 *  visualizations
 */

/* write PostScript drawing of structure to file with annotation */
int PS_rna_plot_snoop_a(char *string,
                        char *structure,
                        char *ssfile,
                        int *relative_access,
                        const char *seqs[]);

/**
 *  \brief Produce a secondary structure graph in PostScript and write it to 'filename'.
 *
 *  Note that this function has changed from previous versions
 *  and now expects the structure to be plotted in dot-bracket notation as an
 *  argument. It does not make use of the global #base_pair array anymore.
 *
 *  \param string     The RNA sequence
 *  \param structure  The secondary structure in dot-bracket notation
 *  \param file       The filename of the postscript output
 *  \return           1 on success, 0 otherwise
 */
int PS_rna_plot(char *string,
                char *structure,
                char *file);

/**
 *  \brief Produce a secondary structure graph in PostScript including additional
 *  annotation macros and write it to 'filename'
 *
 *  Same as PS_rna_plot() but adds extra PostScript macros for various
 *  annotations (see generated PS code). The 'pre' and 'post'
 *  variables contain PostScript code that is verbatim copied in the
 *  resulting PS file just before and after the structure plot.
 *  If both arguments ('pre' and 'post') are NULL, no additional macros will
 *  be printed into the PostScript.
 *
 *  \param string     The RNA sequence
 *  \param structure  The secondary structure in dot-bracket notation
 *  \param file       The filename of the postscript output
 *  \param pre        PostScript code to appear before the secondary structure plot
 *  \param post       PostScript code to appear after the secondary structure plot
 *  \return           1 on success, 0 otherwise
 */
int PS_rna_plot_a(char *string,
                  char *structure,
                  char *file,
                  char *pre,
                  char *post);

int PS_rna_plot_a_gquad(char *string,
                        char *structure,
                        char *ssfile,
                        char *pre,
                        char *post);

/**
 *  \brief Produce a secondary structure graph in Graph Meta Language (gml) and write it to a file
 *
 *  If 'option' is an uppercase letter the RNA sequence is used to label nodes, if 'option' equals
 *  \a 'X' or \a 'x' the resulting file will coordinates for an initial layout of the graph.
 *
 *  \param  string    The RNA sequence
 *  \param  structure The secondary structure in dot-bracket notation
 *  \param  ssfile    The filename of the gml output
 *  \param  option    The option flag
 *  \return           1 on success, 0 otherwise
 */
int gmlRNA( char *string,
            char *structure,
            char *ssfile,
            char option);

/**
 *  \brief  Produce a secondary structure graph in SStructView format
 *
 *  Write coord file for SStructView
 *
 *  \param  string    The RNA sequence
 *  \param  structure The secondary structure in dot-bracket notation
 *  \param  ssfile    The filename of the ssv output
 *  \return           1 on success, 0 otherwise
 */
int ssv_rna_plot( char *string,
                  char *structure,
                  char *ssfile);

/**
 *  \brief Produce a secondary structure plot in SVG format and write it to a file
 *
 *  \param string     The RNA sequence
 *  \param structure  The secondary structure in dot-bracket notation
 *  \param ssfile     The filename of the svg output
 *  \return           1 on success, 0 otherwise
 */
int svg_rna_plot( char *string,
                  char *structure,
                  char *ssfile);

/**
 *  \brief Produce a secondary structure plot for further editing in XRNA
 *
 *  \param string     The RNA sequence
 *  \param structure  The secondary structure in dot-bracket notation
 *  \param ssfile     The filename of the xrna output
 *  \return           1 on success, 0 otherwise
 */
int xrna_plot(char *string,
              char *structure,
              char *ssfile);

int PS_color_dot_plot(char *string,
                      cpair *pi,
                      char *filename);

int PS_color_dot_plot_turn( char *seq,
                            cpair *pi,
                            char *filename,
                            int winSize);

/**
 *  \brief Produce a postscript dot-plot from two pair lists
 *
 *  This function reads two plist structures (e.g. base pair probabilities and a secondary structure)
 *  as produced by assign_plist_from_pr() and assign_plist_from_db() and produces a postscript
 *  "dot plot" that is written to 'filename'.\n
 *  Using base pair probabilities in the first and mfe structure in the second plist, the resulting
 *  "dot plot" represents each base pairing probability by a square of corresponding area in a upper
 *  triangle matrix. The lower part of the matrix contains the minimum free energy structure.
 *
 *  \see assign_plist_from_pr(), assign_plist_from_db()
 *
 *  \param seq      The RNA sequence
 *  \param filename A filename for the postscript output
 *  \param pl       The base pair probability pairlist
 *  \param mf       The mfe secondary structure pairlist
 *  \param comment  A comment
 *  \return         1 if postscript was successfully written, 0 otherwise
 */
int PS_dot_plot_list( char *seq,
                      char *filename,
                      plist *pl,
                      plist *mf,
                      char *comment);

int PS_dot_plot_turn( char *seq,
                      struct plist *pl,
                      char *filename,
                      int winSize);

int PS_color_aln( const char *structure,
                  const char *filename,
                  const char *seqs[],
                  const char *names[]);

/**
 * 	PS_color_aln for duplexes
*/
int aliPS_color_aln(const char *structure,
                    const char *filename, 
                    const char *seqs[],
                    const char *names[]); 


/**
 *  Wrapper to PS_dot_plot_list
 *
 *  \brief Produce postscript dot-plot
 *
 *  Reads base pair probabilities produced by pf_fold() from the
 *  global array #pr and the pair list #base_pair produced by
 *  fold() and produces a postscript "dot plot" that is written to
 *  'filename'. The "dot plot" represents each base pairing
 *  probability by a square of corresponding area in a upper triangle
 *  matrix. The lower part of the matrix contains the minimum free energy
 *  \note DO NOT USE THIS FUNCTION ANYMORE SINCE IT IS NOT THREADSAFE
 *
 *  \deprecated This function is deprecated and will be removed soon! Use \ref PS_dot_plot_list() instead!
 */
DEPRECATED(int PS_dot_plot( char *string,
                            char *file));
#endif
