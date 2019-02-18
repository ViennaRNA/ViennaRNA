#ifndef VIENNA_RNA_PACKAGE_PLOT_PROBABILITIES_H
#define VIENNA_RNA_PACKAGE_PLOT_PROBABILITIES_H


#include <ViennaRNA/datastructures/basic.h>
#include <ViennaRNA/utils/structures.h>

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
 *  @file ViennaRNA/plotting/probabilities.h
 *  @ingroup   plotting_utils
 *  @brief Various functions for plotting RNA secondary structures, dot-plots and other
 *  visualizations
 */

/**
 *  @addtogroup   plotting_utils
 *  @{
 */

#define VRNA_PLOT_PROBABILITIES_BP        1U
#define VRNA_PLOT_PROBABILITIES_ACC       2U

#define VRNA_PLOT_PROBABILITIES_UD        4U
#define VRNA_PLOT_PROBABILITIES_UD_LIN    8U

#define VRNA_PLOT_PROBABILITIES_SD        16U

#define VRNA_PLOT_PROBABILITIES_SC_MOTIF  32U
#define VRNA_PLOT_PROBABILITIES_SC_UP     64U
#define VRNA_PLOT_PROBABILITIES_SC_BP     128U

#define VRNA_PLOT_PROBABILITIES_DEFAULT   (VRNA_PLOT_PROBABILITIES_BP \
                                           | VRNA_PLOT_PROBABILITIES_SD \
                                           | VRNA_PLOT_PROBABILITIES_SC_MOTIF \
                                           | VRNA_PLOT_PROBABILITIES_UD_LIN)
typedef struct {
  char            *comment;
  char            *title;

  vrna_data_lin_t **top;
  char            **top_title;

  vrna_data_lin_t **bottom;
  char            **bottom_title;

  vrna_data_lin_t **left;
  char            **left_title;

  vrna_data_lin_t **right;
  char            **right_title;
} vrna_dotplot_auxdata_t;


int
vrna_plot_dp_EPS(const char             *filename,
                 const char             *sequence,
                 vrna_ep_t              *upper,
                 vrna_ep_t              *lower,
                 vrna_dotplot_auxdata_t *auxdata,
                 unsigned int           options);


int
vrna_plot_dp_PS_list(char   *seq,
                     int    cp,
                     char   *wastlfile,
                     plist  *pl,
                     plist  *mf,
                     char   *comment);


#ifndef VRNA_DISABLE_BACKWARD_COMPATIBILITY

int
PS_color_dot_plot(char  *string,
                  cpair *pi,
                  char  *filename);


int
PS_color_dot_plot_turn(char   *seq,
                       cpair  *pi,
                       char   *filename,
                       int    winSize);


int
PS_dot_plot_turn(char   *seq,
                 plist  *pl,
                 char   *filename,
                 int    winSize);


/**
 *  @brief Produce a postscript dot-plot from two pair lists
 *
 *  This function reads two plist structures (e.g. base pair probabilities and a secondary structure)
 *  as produced by assign_plist_from_pr() and assign_plist_from_db() and produces a postscript
 *  "dot plot" that is written to 'filename'.\n
 *  Using base pair probabilities in the first and mfe structure in the second plist, the resulting
 *  "dot plot" represents each base pairing probability by a square of corresponding area in a upper
 *  triangle matrix. The lower part of the matrix contains the minimum free energy structure.
 *
 *  @see assign_plist_from_pr(), assign_plist_from_db()
 *
 *  @param seq      The RNA sequence
 *  @param filename A filename for the postscript output
 *  @param pl       The base pair probability pairlist
 *  @param mf       The mfe secondary structure pairlist
 *  @param comment  A comment
 *  @return         1 if postscript was successfully written, 0 otherwise
 */
int PS_dot_plot_list(char   *seq,
                     char   *filename,
                     plist  *pl,
                     plist  *mf,
                     char   *comment);


/**
 *  Wrapper to PS_dot_plot_list
 *
 *  @brief Produce postscript dot-plot
 *
 *  Reads base pair probabilities produced by pf_fold() from the
 *  global array #pr and the pair list #base_pair produced by
 *  fold() and produces a postscript "dot plot" that is written to
 *  'filename'. The "dot plot" represents each base pairing
 *  probability by a square of corresponding area in a upper triangle
 *  matrix. The lower part of the matrix contains the minimum free energy
 *  @note DO NOT USE THIS FUNCTION ANYMORE SINCE IT IS NOT THREADSAFE
 *
 *  @deprecated This function is deprecated and will be removed soon! Use @ref PS_dot_plot_list() instead!
 */
DEPRECATED(int PS_dot_plot(char *string,
                           char *file),
           "Use vrna_plot_dp_EPS() instead");

#endif

/**
 * @}
 */

#endif
