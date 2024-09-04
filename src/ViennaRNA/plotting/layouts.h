#ifndef VIENNA_RNA_PACKAGE_PLOT_LAYOUTS_H
#define VIENNA_RNA_PACKAGE_PLOT_LAYOUTS_H

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
 *  @file     ViennaRNA/plotting/layouts.h
 *  @ingroup  plotting_utils
 *  @brief    Secondary structure plot layout algorithms
 */

/**
 *  @addtogroup   plot_layout_utils
 *  @{
 */


/**
 *  @brief  RNA secondary structure figure layout
 *
 *  @see  vrna_plot_layout(), vrna_plot_layout_free(), vrna_plot_layout_simple(),
 *        vrna_plot_layout_circular(), vrna_plot_layout_naview(), vrna_plot_layout_turtle(),
 *        vrna_plot_layout_puzzler()
 */
typedef struct vrna_plot_layout_s vrna_plot_layout_t;


#include <ViennaRNA/datastructures/basic.h>

#ifdef VRNA_WITH_NAVIEW_LAYOUT
#include <ViennaRNA/plotting/naview/naview.h>
#endif

#include "ViennaRNA/plotting/RNApuzzler/RNAturtle.h"
#include "ViennaRNA/plotting/RNApuzzler/RNApuzzler.h"


/**
 *  @brief Definition of Plot type <i>simple</i>
 *
 *  This is the plot type definition for several RNA structure plotting functions telling
 *  them to use <b>Simple</b> plotting algorithm
 *
 *  @see rna_plot_type, vrna_file_PS_rnaplot_a(), vrna_file_PS_rnaplot(), svg_rna_plot(), gmlRNA(), ssv_rna_plot(), xrna_plot()
 */
#define VRNA_PLOT_TYPE_SIMPLE     0U

/**
 *  @brief Definition of Plot type <i>Naview</i>
 *
 *  This is the plot type definition for several RNA structure plotting functions telling
 *  them to use <b>Naview</b> plotting algorithm @rstinline :cite:p:`bruccoleri:1988` @endrst.
 *
 *  @see rna_plot_type, vrna_file_PS_rnaplot_a(), vrna_file_PS_rnaplot(), svg_rna_plot(), gmlRNA(), ssv_rna_plot(), xrna_plot()
 */
#define VRNA_PLOT_TYPE_NAVIEW     1U

/**
 *  @brief Definition of Plot type <i>Circular</i>
 *
 *  This is the plot type definition for several RNA structure plotting functions telling
 *  them to produce a <b>Circular plot</b>
 *
 *  @see rna_plot_type, vrna_file_PS_rnaplot_a(), vrna_file_PS_rnaplot(), svg_rna_plot(), gmlRNA(), ssv_rna_plot(), xrna_plot()
 */
#define VRNA_PLOT_TYPE_CIRCULAR   2U

/**
 *  @brief  Definition of Plot type <i>Turtle</i> @rstinline :cite:p:`wiegreffe:2018` @endrst
 *
 */
#define VRNA_PLOT_TYPE_TURTLE  3U

/**
 *  @brief  Definition of Plot type <i>RNApuzzler</i> @rstinline :cite:p:`wiegreffe:2018` @endrst
 *
 */
#define VRNA_PLOT_TYPE_PUZZLER  4U

#ifdef VRNA_WITH_NAVIEW_LAYOUT
# define VRNA_PLOT_TYPE_DEFAULT  VRNA_PLOT_TYPE_NAVIEW
#else
# define VRNA_PLOT_TYPE_DEFAULT  VRNA_PLOT_TYPE_PUZZLER
#endif


struct vrna_plot_layout_s {
  unsigned int  type;
  unsigned int  length;
  float         *x;
  float         *y;
  double        *arcs;
  int           bbox[4];
};


/**
 *  @brief  Create a layout (coordinates, etc.) for a secondary structure plot
 *
 *  This function can be used to create a secondary structure nucleotide layout
 *  that is then further processed by an actual plotting function. The layout
 *  algorithm can be specified using the @p plot_type parameter, and the following
 *  algorithms are currently supported:
 *  - #VRNA_PLOT_TYPE_SIMPLE
 *  - #VRNA_PLOT_TYPE_NAVIEW
 *  - #VRNA_PLOT_TYPE_CIRCULAR
 *  - #VRNA_PLOT_TYPE_TURTLE
 *  - #VRNA_PLOT_TYPE_PUZZLER
 *
 *  Passing an unsupported selection leads to the default algorithm #VRNA_PLOT_TYPE_NAVIEW
 *
 *  @note If only X-Y coordinates of the corresponding structure layout are required, consider
 *        using vrna_plot_coords() instead!
 *
 *  @see  vrna_plot_layout_free(), vrna_plot_layout_simple(), vrna_plot_layout_naview(),
 *        vrna_plot_layout_circular(), vrna_plot_layout_turtle(), vrna_plot_layout_puzzler(),
 *        vrna_plot_coords(), vrna_file_PS_rnaplot_layout()
 *
 *  @param    structure   The secondary structure in dot-bracket notation
 *  @param    plot_type   The layout algorithm to be used
 *  @return               The layout data structure for the provided secondary structure
 */
vrna_plot_layout_t *
vrna_plot_layout(const char   *structure,
                 unsigned int plot_type);


/**
 *  @brief  Create a layout (coordinates, etc.) for a <i>simple</i> secondary structure plot
 *
 *  This function basically is a wrapper to vrna_plot_layout() that passes the @p plot_type #VRNA_PLOT_TYPE_SIMPLE.
 *
 *  @note If only X-Y coordinates of the corresponding structure layout are required, consider
 *        using vrna_plot_coords_simple() instead!
 *
 *  @see  vrna_plot_layout_free(), vrna_plot_layout(), vrna_plot_layout_naview(),
 *        vrna_plot_layout_circular(), vrna_plot_layout_turtle(), vrna_plot_layout_puzzler(),
 *        vrna_plot_coords_simple(), vrna_file_PS_rnaplot_layout()
 *
 *  @param    structure   The secondary structure in dot-bracket notation
 *  @return               The layout data structure for the provided secondary structure
 */
vrna_plot_layout_t *
vrna_plot_layout_simple(const char *structure);


#ifdef VRNA_WITH_NAVIEW_LAYOUT
/**
 *  @brief  Create a layout (coordinates, etc.) for a secondary structure plot using the <i>Naview Algorithm</i> @rstinline :cite:p:`bruccoleri:1988` @endrst.
 *
 *  This function basically is a wrapper to vrna_plot_layout() that passes the @p plot_type #VRNA_PLOT_TYPE_NAVIEW.
 *
 *  @note If only X-Y coordinates of the corresponding structure layout are required, consider
 *        using vrna_plot_coords_naview() instead!
 *
 *  @see  vrna_plot_layout_free(), vrna_plot_layout(), vrna_plot_layout_simple(),
 *        vrna_plot_layout_circular(), vrna_plot_layout_turtle(), vrna_plot_layout_puzzler(),
 *        vrna_plot_coords_naview(), vrna_file_PS_rnaplot_layout()
 *
 *  @param    structure   The secondary structure in dot-bracket notation
 *  @return               The layout data structure for the provided secondary structure
 */
vrna_plot_layout_t *
vrna_plot_layout_naview(const char *structure);
#endif

/**
 *  @brief  Create a layout (coordinates, etc.) for a <i>circular</i> secondary structure plot
 *
 *  This function basically is a wrapper to vrna_plot_layout() that passes the @p plot_type #VRNA_PLOT_TYPE_CIRCULAR.
 *
 *  @note If only X-Y coordinates of the corresponding structure layout are required, consider
 *        using vrna_plot_coords_circular() instead!
 *
 *  @see  vrna_plot_layout_free(), vrna_plot_layout(), vrna_plot_layout_naview(),
 *        vrna_plot_layout_simple(), vrna_plot_layout_turtle(), vrna_plot_layout_puzzler(),
 *        vrna_plot_coords_circular(), vrna_file_PS_rnaplot_layout()
 *
 *  @param    structure   The secondary structure in dot-bracket notation
 *  @return               The layout data structure for the provided secondary structure
 */
vrna_plot_layout_t *
vrna_plot_layout_circular(const char *structure);


/**
 *  @brief  Create a layout (coordinates, etc.) for a secondary structure plot using the <i>Turtle Algorithm</i> @rstinline :cite:p:`wiegreffe:2018` @endrst
 *
 *  This function basically is a wrapper to vrna_plot_layout() that passes the @p plot_type #VRNA_PLOT_TYPE_TURTLE.
 *
 *  @note If only X-Y coordinates of the corresponding structure layout are required, consider
 *        using vrna_plot_coords_turtle() instead!
 *
 *  @see  vrna_plot_layout_free(), vrna_plot_layout(), vrna_plot_layout_simple(),
 *        vrna_plot_layout_circular(), vrna_plot_layout_naview(), vrna_plot_layout_puzzler(),
 *        vrna_plot_coords_turtle(), vrna_file_PS_rnaplot_layout()
 *
 *  @param    structure   The secondary structure in dot-bracket notation
 *  @return               The layout data structure for the provided secondary structure
 */
vrna_plot_layout_t *
vrna_plot_layout_turtle(const char *structure);


/**
 *  @brief  Create a layout (coordinates, etc.) for a secondary structure plot using the <i>RNApuzzler Algorithm</i> @rstinline :cite:p:`wiegreffe:2018` @endrst
 *
 *  This function basically is a wrapper to vrna_plot_layout() that passes the @p plot_type #VRNA_PLOT_TYPE_PUZZLER.
 *
 *  @note If only X-Y coordinates of the corresponding structure layout are required, consider
 *        using vrna_plot_coords_puzzler() instead!
 *
 *  @see  vrna_plot_layout_free(), vrna_plot_layout(), vrna_plot_layout_simple(),
 *        vrna_plot_layout_circular(), vrna_plot_layout_naview(), vrna_plot_layout_turtle(),
 *        vrna_plot_coords_puzzler(), vrna_file_PS_rnaplot_layout()
 *
 *  @param    structure   The secondary structure in dot-bracket notation
 *  @return               The layout data structure for the provided secondary structure
 */
vrna_plot_layout_t *
vrna_plot_layout_puzzler(const char                   *structure,
                         vrna_plot_options_puzzler_t  *options);


/**
 *  @brief  Free memory occupied by a figure layout data structure
 *
 *  @see  #vrna_plot_layout_t, vrna_plot_layout(), vrna_plot_layout_simple(),
 *        vrna_plot_layout_circular(), vrna_plot_layout_naview(), vrna_plot_layout_turtle(),
 *        vrna_plot_layout_puzzler(), vrna_file_PS_rnaplot_layout()
 *
 *  @param  layout  The layout data structure to free
 */
void
vrna_plot_layout_free(vrna_plot_layout_t *layout);


/**
 *  @brief Compute nucleotide coordinates for secondary structure plot
 *
 *  This function takes a secondary structure and computes X-Y coordinates
 *  for each nucleotide that then can be used to create a structure plot.
 *  The parameter @p plot_type is used to select the underlying layout algorithm.
 *  Currently, the following selections are provided:
 *  - #VRNA_PLOT_TYPE_SIMPLE
 *  - #VRNA_PLOT_TYPE_NAVIEW
 *  - #VRNA_PLOT_TYPE_CIRCULAR
 *  - #VRNA_PLOT_TYPE_TURTLE
 *  - #VRNA_PLOT_TYPE_PUZZLER
 *
 *  Passing an unsupported selection leads to the default algorithm #VRNA_PLOT_TYPE_NAVIEW
 *
 *  Here is a simple example how to use this function, assuming variable @p structure contains
 *  a valid dot-bracket string:
 *  @code{.c}
 *  float *x, *y;
 *
 *  if (vrna_plot_coords(structure, &x, &y)) {
 *    printf("all fine");
 *  } else {
 *    printf("some failure occured!");
 *  }
 *
 *  free(x);
 *  free(y);
 *  @endcode
 *
 *  @note On success, this function allocates memory for X and Y coordinates and assigns
 *        the pointers at addressess @p x and @p y to the corresponding memory locations. It's
 *        the users responsibility to cleanup this memory after usage!
 *
 *  @see  vrna_plot_coords_pt(), vrna_plot_coords_simple(), vrna_plot_coords_naview()
 *        vrna_plot_coords_circular(), vrna_plot_coords_turtle(), vrna_plot_coords_puzzler()
 *
 *  @param        structure   The secondary structure in dot-bracket notation
 *  @param[inout] x           The address of a pointer of X coordinates (pointer will point to memory, or NULL on failure)
 *  @param[inout] y           The address of a pointer of Y coordinates (pointer will point to memory, or NULL on failure)
 *  @param        plot_type   The layout algorithm to be used
 *  @return                   The length of the structure on success, 0 otherwise
 */
int
vrna_plot_coords(const char *structure,
                 float      **x,
                 float      **y,
                 int        plot_type);


/**
 *  @brief Compute nucleotide coordinates for secondary structure plot
 *
 *  Same as vrna_plot_coords() but takes a pair table with the structure
 *  information as input.
 *
 *  @note On success, this function allocates memory for X and Y coordinates and assigns
 *        the pointers at addressess @p x and @p y to the corresponding memory locations. It's
 *        the users responsibility to cleanup this memory after usage!
 *
 *  @see  vrna_plot_coords(), vrna_plot_coords_simple_pt(), vrna_plot_coords_naview_pt()
 *        vrna_plot_coords_circular_pt(), vrna_plot_coords_turtle_pt(), vrna_plot_coords_puzzler_pt()
 *
 *  @param        pt          The pair table that holds the secondary structure
 *  @param[inout] x           The address of a pointer of X coordinates (pointer will point to memory, or NULL on failure)
 *  @param[inout] y           The address of a pointer of Y coordinates (pointer will point to memory, or NULL on failure)
 *  @param        plot_type   The layout algorithm to be used
 *  @return                   The length of the structure on success, 0 otherwise
 */
int
vrna_plot_coords_pt(const short *pt,
                    float       **x,
                    float       **y,
                    int         plot_type);


/**
 *  @brief Compute nucleotide coordinates for secondary structure plot the <i>Simple way</i>
 *
 *  This function basically is a wrapper to vrna_plot_coords() that passes the @p plot_type #VRNA_PLOT_TYPE_SIMPLE.
 *
 *  Here is a simple example how to use this function, assuming variable @p structure contains
 *  a valid dot-bracket string:
 *  @code{.c}
 *  float *x, *y;
 *
 *  if (vrna_plot_coords_simple(structure, &x, &y)) {
 *    printf("all fine");
 *  } else {
 *    printf("some failure occured!");
 *  }
 *
 *  free(x);
 *  free(y);
 *  @endcode
 *
 *  @note On success, this function allocates memory for X and Y coordinates and assigns
 *        the pointers at addressess @p x and @p y to the corresponding memory locations. It's
 *        the users responsibility to cleanup this memory after usage!
 *
 *  @see  vrna_plot_coords(), vrna_plot_coords_simple_pt(), vrna_plot_coords_circular(),
 *        vrna_plot_coords_naview(), vrna_plot_coords_turtle(), vrna_plot_coords_puzzler()
 *
 *  @param        structure   The secondary structure in dot-bracket notation
 *  @param[inout] x           The address of a pointer of X coordinates (pointer will point to memory, or NULL on failure)
 *  @param[inout] y           The address of a pointer of Y coordinates (pointer will point to memory, or NULL on failure)
 *  @return                   The length of the structure on success, 0 otherwise
 */
int
vrna_plot_coords_simple(const char  *structure,
                        float       **x,
                        float       **y);


/**
 *  @brief Compute nucleotide coordinates for secondary structure plot the <i>Simple way</i>
 *
 *  Same as vrna_plot_coords_simple() but takes a pair table with the structure
 *  information as input.
 *
 *  @note On success, this function allocates memory for X and Y coordinates and assigns
 *        the pointers at addressess @p x and @p y to the corresponding memory locations. It's
 *        the users responsibility to cleanup this memory after usage!
 *
 *  @see  vrna_plot_coords_pt(), vrna_plot_coords_simple(), vrna_plot_coords_circular_pt(),
 *        vrna_plot_coords_naview_pt(), vrna_plot_coords_turtle_pt(), vrna_plot_coords_puzzler_pt()
 *
 *  @param        pt          The pair table that holds the secondary structure
 *  @param[inout] x           The address of a pointer of X coordinates (pointer will point to memory, or NULL on failure)
 *  @param[inout] y           The address of a pointer of Y coordinates (pointer will point to memory, or NULL on failure)
 *  @return                   The length of the structure on success, 0 otherwise
 */
int
vrna_plot_coords_simple_pt(const short  *pt,
                           float        **x,
                           float        **y);


/**
 *  @brief Compute coordinates of nucleotides mapped in equal distancies onto a unit circle.
 *
 *  This function basically is a wrapper to vrna_plot_coords() that passes the @p plot_type #VRNA_PLOT_TYPE_CIRCULAR.
 *
 *  In order to draw nice arcs using quadratic bezier curves that connect base pairs one may calculate
 *  a second tangential point @f$P^t@f$ in addition to the actual R<sup>2</sup> coordinates.
 *  the simplest way to do so may be to compute a radius scaling factor @f$rs@f$ in the interval @f$[0,1]@f$ that
 *  weights the proportion of base pair span to the actual length of the sequence. This scaling factor
 *  can then be used to calculate the coordinates for @f$P^t@f$, i.e.
 *
 *  @f[ P^{t}_x[i] = X[i] * rs @f]
 *
 *  and
 *
 *  @f[ P^{t}_y[i] = Y[i] * rs @f].
 *
 *  @note On success, this function allocates memory for X and Y coordinates and assigns
 *        the pointers at addressess @p x and @p y to the corresponding memory locations. It's
 *        the users responsibility to cleanup this memory after usage!
 *
 *  @see  vrna_plot_coords(), vrna_plot_coords_circular_pt(), vrna_plot_coords_simple(),
 *        vrna_plot_coords_naview(), vrna_plot_coords_turtle(), vrna_plot_coords_puzzler()
 *
 *  @param        structure   The secondary structure in dot-bracket notation
 *  @param[inout] x           The address of a pointer of X coordinates (pointer will point to memory, or NULL on failure)
 *  @param[inout] y           The address of a pointer of Y coordinates (pointer will point to memory, or NULL on failure)
 *  @return                   The length of the structure on success, 0 otherwise
 */
int
vrna_plot_coords_circular(const char  *structure,
                          float       **x,
                          float       **y);


/**
 *  @brief Compute nucleotide coordinates for a <i>Circular Plot</i>
 *
 *  Same as vrna_plot_coords_circular() but takes a pair table with the structure
 *  information as input.
 *
 *  @note On success, this function allocates memory for X and Y coordinates and assigns
 *        the pointers at addressess @p x and @p y to the corresponding memory locations. It's
 *        the users responsibility to cleanup this memory after usage!
 *
 *  @see  vrna_plot_coords_pt(), vrna_plot_coords_circular(), vrna_plot_coords_simple_pt(),
 *        vrna_plot_coords_naview_pt(), vrna_plot_coords_turtle_pt(), vrna_plot_coords_puzzler_pt()
 *
 *  @param        pt          The pair table that holds the secondary structure
 *  @param[inout] x           The address of a pointer of X coordinates (pointer will point to memory, or NULL on failure)
 *  @param[inout] y           The address of a pointer of Y coordinates (pointer will point to memory, or NULL on failure)
 *  @return                   The length of the structure on success, 0 otherwise
 */
int
vrna_plot_coords_circular_pt(const short  *pt,
                             float        **x,
                             float        **y);


/**
 * @}
 */


#ifndef VRNA_DISABLE_BACKWARD_COMPATIBILITY

/**
 *  @addtogroup   plotting_utils_deprecated
 *  @{
 */


/**
 *  @brief this is a workarround for the SWIG Perl Wrapper RNA plot function
 *  that returns an array of type COORDINATE
 */
typedef struct {
  float X;  /* X coords */
  float Y;  /* Y coords */
} COORDINATE;


/**
 *  @brief Switch for changing the secondary structure layout algorithm
 *
 *  Current possibility are 0 for a simple radial drawing or 1 for the modified
 *  radial drawing taken from the @e naview program of @rstinline :cite:t:`bruccoleri:1988` @endrst.
 *
 *  @note To provide thread safety please do not rely on this global variable in future implementations
 *        but pass a plot type flag directly to the function that decides which layout algorithm it may use!
 *
 *  @see #VRNA_PLOT_TYPE_SIMPLE, #VRNA_PLOT_TYPE_NAVIEW, #VRNA_PLOT_TYPE_CIRCULAR
 *
 */
extern int rna_plot_type;


/**
 *  @brief Calculate nucleotide coordinates for secondary structure plot the <i>Simple way</i>
 *
 *  @see make_pair_table(), rna_plot_type, simple_circplot_coordinates(), naview_xy_coordinates(), vrna_file_PS_rnaplot_a(),
 *  vrna_file_PS_rnaplot, svg_rna_plot()
 *
 *  @deprecated   Consider switching to vrna_plot_coords_simple_pt() instead!
 *
 *  @param  pair_table  The pair table of the secondary structure
 *  @param  X           a pointer to an array with enough allocated space to hold the x coordinates
 *  @param  Y           a pointer to an array with enough allocated space to hold the y coordinates
 *  @return             length of sequence on success, 0 otherwise
 */
DEPRECATED(int
           simple_xy_coordinates(short  *pair_table,
                                 float  *X,
                                 float  *Y),
           "Use vrna_plot_coords_simple_pt() instead!");


/**
 *  @brief Calculate nucleotide coordinates for <i>Circular Plot</i>
 *
 *  This function calculates the coordinates of nucleotides mapped in equal distancies onto a unit circle.
 *
 *  @note In order to draw nice arcs using quadratic bezier curves that connect base pairs one may calculate
 *        a second tangential point @f$P^t@f$ in addition to the actual R<sup>2</sup> coordinates.
 *        the simplest way to do so may be to compute a radius scaling factor @f$rs@f$ in the interval @f$[0,1]@f$ that
 *        weights the proportion of base pair span to the actual length of the sequence. This scaling factor
 *        can then be used to calculate the coordinates for @f$P^t@f$, i.e. @f$ P^{t}_x[i] = X[i] * rs@f$
 *        and @f$P^{t}_y[i] = Y[i] * rs@f$.
 *
 *  @see  make_pair_table(), rna_plot_type, simple_xy_coordinates(), naview_xy_coordinates(),
 *        vrna_file_PS_rnaplot_a(), vrna_file_PS_rnaplot, svg_rna_plot()
 *
 *  @deprecated   Consider switching to vrna_plot_coords_circular_pt() instead!
 *
 *  @param  pair_table  The pair table of the secondary structure
 *  @param  x           a pointer to an array with enough allocated space to hold the x coordinates
 *  @param  y           a pointer to an array with enough allocated space to hold the y coordinates
 *  @return             length of sequence on success, 0 otherwise
 */
DEPRECATED(int
           simple_circplot_coordinates(short  *pair_table,
                                       float  *x,
                                       float  *y),
           "Use vrna_plot_coords_circular_pt() instead!");


/**
 * @}
 */

#endif


#endif
