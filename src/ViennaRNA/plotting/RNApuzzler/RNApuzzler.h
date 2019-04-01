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
} puzzlerOptions;


/**
 * @brief
 *      Constructor.
 * @return
 */
puzzlerOptions *
createPuzzlerOptions();


/**
 * @brief
 *      Destructor.
 * @param puzzler
 */
void
destroyPuzzlerOptions(puzzlerOptions *puzzler);


/**
 * Compute layout using RNApuzzler algorithm
 */
int
layout_RNApuzzler(short const *const  pair_table,
                  float               *x,
                  float               *y,
                  double              *arc_coords,
                  puzzlerOptions      *puzzler);


#endif
