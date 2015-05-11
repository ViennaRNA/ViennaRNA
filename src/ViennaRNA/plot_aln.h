#ifndef VIENNA_RNA_PACKAGE_PLOT_ALN_H
#define VIENNA_RNA_PACKAGE_PLOT_ALN_H

#ifdef DEPRECATION_WARNINGS
# ifdef __GNUC__
#  define DEPRECATED(func) func __attribute__ ((deprecated))
# else
#  define DEPRECATED(func) func
# endif
#else
# define DEPRECATED(func) func
#endif

/* make this interface backward compatible with RNAlib < 2.2.0 */
#define VRNA_BACKWARD_COMPAT

/**
 *  @addtogroup   plotting_utils
 *
 *  @{
 *
 *  @file plot_aln.h
 *  @brief Various functions for plotting Sequence / Structure Alignments
 */

/**
 *  @brief Produce PostScript sequence alignment color-annotated by consensus
 *  structure
*/
int PS_color_aln( const char *structure,
                  const char *filename,
                  const char *seqs[],
                  const char *names[]);

/**
 *  PS_color_aln for duplexes
*/
int aliPS_color_aln(const char *structure,
                    const char *filename, 
                    const char *seqs[],
                    const char *names[]); 

#ifdef VRNA_BACKWARD_COMPAT

#endif

/**
 * @}
 */

#endif
