#ifndef RNATURTLE_H
#define RNATURTLE_H

/**
 *  @file ViennaRNA/plotting/RNApuzzler/RNAturtle.h
 *  @ingroup  plotting_utils
 *  @brief    Implementation of the RNAturtle RNA secondary structure layout algorithm @rstinline :cite:p:`wiegreffe:2018` @endrst
 *
 */

/**
 *  @addtogroup   plot_layout_utils
 *  @{
 */


/**
 *  @brief Compute nucleotide coordinates for secondary structure plot using the <i>RNAturtle</i> algorithm @rstinline :cite:p:`wiegreffe:2018` @endrst
 *
 *  This function basically is a wrapper to vrna_plot_coords() that passes the @p plot_type #VRNA_PLOT_TYPE_TURTLE.
 *
 *  Here is a simple example how to use this function, assuming variable @p structure contains
 *  a valid dot-bracket string:
 *  @code{.c}
 *  float  *x, *y;
 *  double *arcs;
 *
 *  if (vrna_plot_coords_turtle(structure, &x, &y, &arcs)) {
 *    printf("all fine");
 *  } else {
 *    printf("some failure occured!");
 *  }
 *
 *  free(x);
 *  free(y);
 *  free(arcs);
 *  @endcode
 *
 *  @note On success, this function allocates memory for X, Y and arc coordinates and assigns
 *  the pointers at addressess @p x, @p y and @p arc_coords to the corresponding memory locations. It's
 *  the users responsibility to cleanup this memory after usage!
 *
 *  @see  vrna_plot_coords(), vrna_plot_coords_turtle_pt(), vrna_plot_coords_circular(),
 *        vrna_plot_coords_simple(), vrna_plot_coords_naview(), vrna_plot_coords_puzzler()
 *
 *  @param        structure   The secondary structure in dot-bracket notation
 *  @param[inout] x           The address of a pointer of X coordinates (pointer will point to memory, or NULL on failure)
 *  @param[inout] y           The address of a pointer of Y coordinates (pointer will point to memory, or NULL on failure)
 *  @param[inout] arc_coords  The address of a pointer that will hold arc coordinates (pointer will point to memory, or NULL on failure)
 *  @return                   The length of the structure on success, 0 otherwise
 */
int
vrna_plot_coords_turtle(const char  *structure,
                        float       **x,
                        float       **y,
                        double      **arc_coords);


/**
 *  @brief Compute nucleotide coordinates for secondary structure plot using the <i>RNAturtle</i> algorithm @rstinline :cite:p:`wiegreffe:2018` @endrst
 *
 *  Same as vrna_plot_coords_turtle() but takes a pair table with the structure
 *  information as input.
 *
 *  @note On success, this function allocates memory for X, Y and arc coordinates and assigns
 *  the pointers at addressess @p x, @p y and @p arc_coords to the corresponding memory locations. It's
 *  the users responsibility to cleanup this memory after usage!
 *
 *  @see  vrna_plot_coords_pt(), vrna_plot_coords_turtle(), vrna_plot_coords_circular_pt(),
 *        vrna_plot_coords_simple_pt(), vrna_plot_coords_puzzler_pt(), vrna_plot_coords_naview_pt()
 *
 *  @param        pt          The pair table that holds the secondary structure
 *  @param[inout] x           The address of a pointer of X coordinates (pointer will point to memory, or NULL on failure)
 *  @param[inout] y           The address of a pointer of Y coordinates (pointer will point to memory, or NULL on failure)
 *  @param[inout] arc_coords  The address of a pointer that will hold arc coordinates (pointer will point to memory, or NULL on failure)
 *  @return                   The length of the structure on success, 0 otherwise
 */
int
vrna_plot_coords_turtle_pt(short const *const pair_table,
                           float              **x,
                           float              **y,
                           double             **arc_coords);


/**
 * @}
 */

#endif
