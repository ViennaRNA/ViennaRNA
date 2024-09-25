#ifndef VIENNA_RNA_PACKAGE_PLOT_PROBABILITIES_H
#define VIENNA_RNA_PACKAGE_PLOT_PROBABILITIES_H


#include <ViennaRNA/datastructures/basic.h>
#include <ViennaRNA/structures/problist.h>

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
 *  @addtogroup   plot_probabilities
 *  @{
 */

/**
 *  @brief  Option flag for base pair probabilities in probability plot output functions
 */
#define VRNA_PLOT_PROBABILITIES_BP        1U


/**
 *  @brief  Option flag for accessibilities in probability plot output functions
 */
#define VRNA_PLOT_PROBABILITIES_ACC       2U


/**
 *  @brief  Option flag for unstructured domain probabilities in probability plot output functions
 */
#define VRNA_PLOT_PROBABILITIES_UD        4U


/**
 *  @brief  Option flag for unstructured domain probabilities (linear representation) in probability plot output functions
 */
#define VRNA_PLOT_PROBABILITIES_UD_LIN    8U


/**
 *  @brief  Option flag for structured domain probabilities (such as G-quadruplexes) in probability plot output functions
 */
#define VRNA_PLOT_PROBABILITIES_SD        16U


/**
 *  @brief  Option flag for soft-constraint motif probabilities in probability plot output functions
 */
#define VRNA_PLOT_PROBABILITIES_SC_MOTIF  32U
#define VRNA_PLOT_PROBABILITIES_SC_UP     64U
#define VRNA_PLOT_PROBABILITIES_SC_BP     128U

/**
 *  @brief  Default option flag for probability plot output functions
 *
 *  Default output includes actual base pair probabilties (#VRNA_PLOT_PROBABILITIES_BP),
 *  structured domain probabilities such as G-quadruplexes (#VRNA_PLOT_PROBABILITIES_SD),
 *  probabilities obtained from soft-constraint motif implementation (#VRNA_PLOT_PROBABILITIES_SC_MOTIF),
 *  and unstructured domain probabilities (#VRNA_PLOT_PROBABILITIES_UD_LIN).
 *
 *  @see  vrna_plot_dp_EPS()
 */
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


/**
 *  @brief Produce an encapsulate PostScript (EPS) dot-plot from one or two lists of base pair probabilities
 *
 *  This function reads two vrna_ep_t lists @p upper and @p lower (e.g. base pair probabilities
 *  and a secondary structure) and produces an EPS "dot plot" with filename @p 'filename' where
 *  data from @p upper is placed in the upper-triangular and data from @p lower is placed in
 *  the lower triangular part of the matrix.\n
 *
 *  For default output, provide the flag #VRNA_PLOT_PROBABILITIES_DEFAULT as @p options parameter.
 *
 *  @see vrna_plist(), vrna_plist_from_probs(),
 *       #VRNA_PLOT_PROBABILITIES_DEFAULT
 *
 *  @param filename A filename for the EPS output
 *  @param sequence The RNA sequence
 *  @param upper    The base pair probabilities for the upper triangular part
 *  @param lower    The base pair probabilities for the lower triangular part
 *  @param options  Options indicating which of the input data should be included in the dot-plot
 *  @return         1 if EPS file was successfully written, 0 otherwise
 */
int
vrna_plot_dp_EPS(const char             *filename,
                 const char             *sequence,
                 vrna_ep_t              *upper,
                 vrna_ep_t              *lower,
                 vrna_dotplot_auxdata_t *auxdata,
                 unsigned int           options);


/**
 *  @brief Produce a postscript dot-plot from two pair lists
 *
 *  This function reads two plist structures (e.g. base pair probabilities and a secondary structure)
 *  as produced by vrna_plist_from_probs() and vrna_plist() and produces a postscript
 *  "dot plot" that is written to 'filename'.\n
 *  Using base pair probabilities in the first and mfe structure in the second plist, the resulting
 *  "dot plot" represents each base pairing probability by a square of corresponding area in a upper
 *  triangle matrix. The lower part of the matrix contains the minimum free energy structure.
 *
 *  @see vrna_plist_from_probs(), vrna_plist()
 *
 *  @param seq      The RNA sequence
 *  @param filename A filename for the postscript output
 *  @param pl       The base pair probability pairlist
 *  @param mf       The mfe secondary structure pairlist
 *  @param comment  A comment
 *  @return         1 if postscript was successfully written, 0 otherwise
 */
int
vrna_plot_dp_PS_list(char       *seq,
                     int        cp,
                     char       *filename,
                     vrna_ep_t  *pl,
                     vrna_ep_t  *mf,
                     char       *comment);


/**
 * @}
 */

#ifndef VRNA_DISABLE_BACKWARD_COMPATIBILITY

/**
 *  @addtogroup   plotting_utils_deprecated
 *  @{
 */

int
PS_color_dot_plot(char          *string,
                  vrna_cpair_t  *pi,
                  char          *filename);


int
PS_color_dot_plot_turn(char         *seq,
                       vrna_cpair_t *pi,
                       char         *filename,
                       int          winSize);


int
PS_dot_plot_turn(char       *seq,
                 vrna_ep_t  *pl,
                 char       *filename,
                 int        winSize);


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
DEPRECATED(int PS_dot_plot_list(char      *seq,
                                char      *filename,
                                vrna_ep_t *pl,
                                vrna_ep_t *mf,
                                char      *comment),
           "Use vrna_plot_dp_PS_list() instead");


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

/**
 * @}
 */

#endif


#endif
