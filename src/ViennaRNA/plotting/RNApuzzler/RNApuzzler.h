#ifndef RNAPUZZLER_H
#define RNAPUZZLER_H

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
 * @brief
 *      Constructor.
 * @return
 */
vrna_plot_options_puzzler_t *
vrna_plot_options_puzzler();


/**
 * @brief
 *      Destructor.
 * @param puzzler
 */
void
vrna_plot_options_puzzler_free(vrna_plot_options_puzzler_t *puzzler);


/**
 * Compute layout using RNApuzzler algorithm
 */
int
vrna_plot_coords_puzzler(const char                   *structure,
                         float                        **x,
                         float                        **y,
                         double                       **arc_coords,
                         vrna_plot_options_puzzler_t  *options);


int
vrna_plot_coords_puzzler_pt(short const *const          pair_table,
                            float                       **x,
                            float                       **y,
                            double                      **arc_coords,
                            vrna_plot_options_puzzler_t *puzzler);


#endif
