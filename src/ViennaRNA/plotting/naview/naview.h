#ifndef VIENNA_RNA_PACKAGE_PLOT_NAVIEW_H
#define VIENNA_RNA_PACKAGE_PLOT_NAVIEW_H

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
 *  @file ViennaRNA/plotting/naview.h
 *  @ingroup  plotting_utils
 *  @brief    Implementation of the Naview RNA secondary structure layout algorithm @rstinline :cite:p:`bruccoleri:1988` @endrst
 *
 */

/**
 *  @addtogroup   plot_layout_utils
 *  @{
 */


/**
 *  @brief Compute nucleotide coordinates for secondary structure plot using the <i>Naview</i> algorithm @rstinline :cite:p:`bruccoleri:1988` @endrst
 *
 *  This function basically is a wrapper to vrna_plot_coords() that passes the @p plot_type #VRNA_PLOT_TYPE_NAVIEW.
 *
 *  Here is a simple example how to use this function, assuming variable @p structure contains
 *  a valid dot-bracket string:
 *  @code{.c}
 *  float *x, *y;
 *
 *  if (vrna_plot_coords_naview(structure, &x, &y)) {
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
 *  the pointers at addressess @p x and @p y to the corresponding memory locations. It's
 *  the users responsibility to cleanup this memory after usage!
 *
 *  @see  vrna_plot_coords(), vrna_plot_coords_simple_pt(), vrna_plot_coords_circular(),
 *        vrna_plot_coords_simple(), vrna_plot_coords_turtle(), vrna_plot_coords_puzzler()
 *
 *  @param        structure   The secondary structure in dot-bracket notation
 *  @param[inout] x           The address of a pointer of X coordinates (pointer will point to memory, or NULL on failure)
 *  @param[inout] y           The address of a pointer of Y coordinates (pointer will point to memory, or NULL on failure)
 *  @return                   The length of the structure on success, 0 otherwise
 */
int
vrna_plot_coords_naview(const char  *structure,
                        float       **x,
                        float       **y);


/**
 *  @brief Compute nucleotide coordinates for secondary structure plot using the <i>Naview</i> algorithm @rstinline :cite:p:`bruccoleri:1988` @endrst
 *
 *  Same as vrna_plot_coords_naview() but takes a pair table with the structure
 *  information as input.
 *
 *  @note On success, this function allocates memory for X and Y coordinates and assigns
 *  the pointers at addressess @p x and @p y to the corresponding memory locations. It's
 *  the users responsibility to cleanup this memory after usage!
 *
 *  @see  vrna_plot_coords_pt(), vrna_plot_coords_naview(), vrna_plot_coords_circular_pt(),
 *        vrna_plot_coords_simple_pt(), vrna_plot_coords_turtle_pt(), vrna_plot_coords_puzzler_pt()
 *
 *  @param        pt          The pair table that holds the secondary structure
 *  @param[inout] x           The address of a pointer of X coordinates (pointer will point to memory, or NULL on failure)
 *  @param[inout] y           The address of a pointer of Y coordinates (pointer will point to memory, or NULL on failure)
 *  @return                   The length of the structure on success, 0 otherwise
 */
int
vrna_plot_coords_naview_pt(const short  *pt,
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
 *  @deprecated   Consider using vrna_plot_coords_naview_pt() instead!
 */
DEPRECATED(int
           naview_xy_coordinates(short  *pair_table,
                                 float  *X,
                                 float  *Y),
           "Use vrna_plot_coords_naview_pt() instead!");

#endif

/**
 * @}
 */


#endif
