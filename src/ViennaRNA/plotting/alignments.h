#ifndef VIENNA_RNA_PACKAGE_PLOT_ALN_H
#define VIENNA_RNA_PACKAGE_PLOT_ALN_H

/**
 *  @file     ViennaRNA/plotting/alignments.h
 *  @ingroup  plotting_utils
 *  @brief    Various functions for plotting Sequence / Structure Alignments
 */

/**
 *  @addtogroup  plotting_utils
 *  @{
 */

/**
 *  @brief Produce PostScript sequence alignment color-annotated by consensus
 *  structure
 */
int PS_color_aln(const char *structure,
                 const char *filename,
                 const char *seqs[],
                 const char *names[]);


/**
 *  @param columns  The number of columns before the alignment is wrapped as a new block (values less than 1 indicate no wrapping)
 */
int
vrna_file_PS_aln(const char *filename,
                 const char **seqs,
                 const char **names,
                 const char *structure,
                 int        columns);


/**
 *  @param columns  The number of columns before the alignment is wrapped as a new block (values less than 1 indicate no wrapping)
 */
int
vrna_file_PS_aln_sub(const char *filename,
                     const char **seqs,
                     const char **names,
                     const char *structure,
                     int        start,
                     int        end,
                     int        columns);


/**
 *  PS_color_aln for duplexes
 */
int aliPS_color_aln(const char  *structure,
                    const char  *filename,
                    const char  *seqs[],
                    const char  *names[]);


/**
 * @}
 */

#endif
