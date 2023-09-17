#ifndef RNAPUZZLER_H
#define RNAPUZZLER_H

/**
 *  @file ViennaRNA/plotting/RNApuzzler/RNApuzzler.h
 *  @ingroup  plotting_utils
 *  @brief    Implementation of the RNApuzzler RNA secondary structure layout algorithm @rstinline :cite:p:`wiegreffe:2018` @endrst
 *
 */

/**
 *  @addtogroup   plot_layout_utils
 *  @{
 */


/**
 *  @brief  Options data structure for RNApuzzler algorithm implementation
 */
typedef struct {
  /*
   * variables fixed during operation
   * drawing behavior
   */
  short       drawArcs;
  double      paired;
  double      unpaired;

  /* intersection resolution behavior */
  short       checkAncestorIntersections;
  short       checkSiblingIntersections;
  short       checkExteriorIntersections;
  short       allowFlipping;
  short       optimize;
  int         maximumNumberOfConfigChangesAllowed;


  /* import behavior - unused for now */
  char        *config; /* file path */

  /* other stuff */
  const char  *filename;

  /* variables changed during operation */
  int         numberOfChangesAppliedToConfig;
  int         psNumber;
} vrna_plot_options_puzzler_t;


/**
 *  @brief Compute nucleotide coordinates for secondary structure plot using the <i>RNApuzzler</i> algorithm @rstinline :cite:p:`wiegreffe:2018` @endrst
 *
 *  This function basically is a wrapper to vrna_plot_coords() that passes the @p plot_type #VRNA_PLOT_TYPE_PUZZLER.
 *
 *  Here is a simple example how to use this function, assuming variable @p structure contains
 *  a valid dot-bracket string and using the default options (@p options = NULL):
 *  @code{.c}
 *  float  *x, *y;
 *  double *arcs;
 *
 *  if (vrna_plot_coords_puzzler(structure, &x, &y, &arcs, NULL)) {
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
 *  @see  vrna_plot_coords(), vrna_plot_coords_puzzler_pt(), vrna_plot_coords_circular(),
 *        vrna_plot_coords_simple(), vrna_plot_coords_turtle(), vrna_plot_coords_naview(),
 *        vrna_plot_options_puzzler()
 *
 *  @param        structure   The secondary structure in dot-bracket notation
 *  @param[inout] x           The address of a pointer of X coordinates (pointer will point to memory, or NULL on failure)
 *  @param[inout] y           The address of a pointer of Y coordinates (pointer will point to memory, or NULL on failure)
 *  @param[inout] arc_coords  The address of a pointer that will hold arc coordinates (pointer will point to memory, or NULL on failure)
 *  @param        options     The options for the RNApuzzler algorithm (or NULL)
 *  @return                   The length of the structure on success, 0 otherwise
 */
int
vrna_plot_coords_puzzler(const char                   *structure,
                         float                        **x,
                         float                        **y,
                         double                       **arc_coords,
                         vrna_plot_options_puzzler_t  *options);


/**
 *  @brief Compute nucleotide coordinates for secondary structure plot using the <i>RNApuzzler</i> algorithm @rstinline :cite:p:`wiegreffe:2018` @endrst
 *
 *  Same as vrna_plot_coords_puzzler() but takes a pair table with the structure
 *  information as input.
 *
 *  @note On success, this function allocates memory for X, Y and arc coordinates and assigns
 *  the pointers at addressess @p x, @p y and @p arc_coords to the corresponding memory locations. It's
 *  the users responsibility to cleanup this memory after usage!
 *
 *  @see  vrna_plot_coords_pt(), vrna_plot_coords_puzzler(), vrna_plot_coords_circular_pt(),
 *        vrna_plot_coords_simple_pt(), vrna_plot_coords_turtle_pt(), vrna_plot_coords_naview_pt()
 *
 *  @param        pt          The pair table that holds the secondary structure
 *  @param[inout] x           The address of a pointer of X coordinates (pointer will point to memory, or NULL on failure)
 *  @param[inout] y           The address of a pointer of Y coordinates (pointer will point to memory, or NULL on failure)
 *  @param[inout] arc_coords  The address of a pointer that will hold arc coordinates (pointer will point to memory, or NULL on failure)
 *  @param        options     The options for the RNApuzzler algorithm (or NULL)
 *  @return                   The length of the structure on success, 0 otherwise
 */
int
vrna_plot_coords_puzzler_pt(short const *const          pair_table,
                            float                       **x,
                            float                       **y,
                            double                      **arc_coords,
                            vrna_plot_options_puzzler_t *puzzler);


/**
 *  @brief  Create an RNApuzzler options data structure
 *
 *  @see    vrna_plot_options_puzzler_free(), vrna_plot_coords_puzzler(), vrna_plot_coords_puzzler_pt(),
 *          vrna_plot_layout_puzzler()
 *
 *  @return An RNApuzzler options data structure with default settings
 */
vrna_plot_options_puzzler_t *
vrna_plot_options_puzzler(void);


/**
 *  @brief  Free memory occupied by an RNApuzzler options data structure
 *
 *  @see    vrna_plot_options_puzzler(), vrna_plot_coords_puzzler(), vrna_plot_coords_puzzler_pt(),
 *          vrna_plot_layout_puzzler()
 *
 *  @param  options   A pointer to the options data structure to free
 */
void
vrna_plot_options_puzzler_free(vrna_plot_options_puzzler_t *options);


/**
 * @}
 */


#endif
