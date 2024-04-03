#ifndef VIENNA_RNA_PACKAGE_PLOT_ALN_H
#define VIENNA_RNA_PACKAGE_PLOT_ALN_H

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
 *  @file     ViennaRNA/plotting/alignments.h
 *  @ingroup  plotting_utils
 *  @brief    Various functions for plotting Sequence / Structure Alignments
 */

/**
 *  @addtogroup  alignment_plots
 *  @{
 */

typedef struct {
  unsigned int start;
  unsigned int end;
  unsigned int offset;
  unsigned int columns;
  double       color_threshold;
  double       color_min_sat;
} vrna_aln_opt_t;

/**
 *  @brief  Create an annotated PostScript alignment plot
 *
 *  @see vrna_file_PS_aln_slice()
 *
 *  @param  filename  The output file name
 *  @param  seqs      The aligned sequences
 *  @param  names     The names of the sequences
 *  @param  structure The consensus structure in dot-bracket notation
 *  @param  columns   The number of columns before the alignment is wrapped as a new block (a value of 0 indicates no wrapping)
 */
int
vrna_file_PS_aln(const char   *filename,
                 const char   **seqs,
                 const char   **names,
                 const char   *structure,
                 unsigned int columns);


/**
 *  @brief  Create an annotated PostScript alignment plot
 *
 *  Similar to vrna_file_PS_aln() but allows the user to print a particular slice
 *  of the alignment by specifying a @p start and @p end position. The additional
 *  @p offset parameter allows for adjusting the alignment position ruler value.
 *
 *  @see vrna_file_PS_aln_slice()
 *
 *  @param  filename  The output file name
 *  @param  seqs      The aligned sequences
 *  @param  names     The names of the sequences
 *  @param  structure The consensus structure in dot-bracket notation
 *  @param  start     The start of the alignment slice (a value of 0 indicates the first position of the alignment, i.e. no slicing at 5' side)
 *  @param  end       The end of the alignment slice (a value of 0 indicates the last position of the alignment, i.e. no slicing at 3' side)
 *  @param  offset    The alignment coordinate offset for the position ruler.
 *  @param  columns   The number of columns before the alignment is wrapped as a new block (a value of 0 indicates no wrapping)
 */
int
vrna_file_PS_aln_slice(const char   *filename,
                       const char   **seqs,
                       const char   **names,
                       const char   *structure,
                       unsigned int start,
                       unsigned int end,
                       int          offset,
                       unsigned int columns);

int
vrna_file_PS_aln_opt(const char   *filename,
                     const char   **seqs,
                     const char   **names,
                     const char   *structure,
                     vrna_aln_opt_t options);

/**
 * @}
 */

#ifndef VRNA_DISABLE_BACKWARD_COMPATIBILITY

/**
 *  @brief Produce PostScript sequence alignment color-annotated by consensus
 *  structure
 *
 *  @deprecated Use vrna_file_PS_aln() instead!
 *  @ingroup    plotting_utils_deprecated
 */
DEPRECATED(int PS_color_aln(const char  *structure,
                            const char  *filename,
                            const char  *seqs[],
                            const char  *names[]),
           "Use vrna_file_PS_aln() instead!");


/**
 *  @brief PS_color_aln for duplexes
 *
 *  @deprecated Use vrna_file_PS_aln() instead!
 *  @ingroup    plotting_utils_deprecated
 */
DEPRECATED(int aliPS_color_aln(const char *structure,
                               const char *filename,
                               const char *seqs[],
                               const char *names[]),
           "Use vrna_file_PS_aln() instead!");

#endif

#endif
