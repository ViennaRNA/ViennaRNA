#ifndef __RNA_WALK__
#define __RNA_WALK__

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ViennaRNA/fold_vars.h>
#include <ViennaRNA/data_structures.h>
#include <ViennaRNA/fold.h>
#include <ViennaRNA/utils.h>
#include <ViennaRNA/structure_utils.h>
#include <ViennaRNA/energy_const.h>
#include <ViennaRNA/model.h>
#include "meshpoint.h"

#define MIN2(A, B)        ((A) < (B) ? (A) : (B))
#define MAX2(A, B)        ((A) > (B) ? (A) : (B))

/* simulated annealing related variables */
extern int        simulatedAnnealing;
extern FLT_OR_DBL tstart;
extern FLT_OR_DBL tstop;
extern float      treduction;

/* Monte Carlo related variables */
extern int        rememberStructures;
extern int        maxRest;
extern int        backWalkPenalty;

/* init function that initializes the random number generator
 *  and also constructs the S and S1 sequence encoding arrays
 *  and the pair table for further usage
 */
void initRNAWalk(const char *seq,
                 vrna_md_t  *md);


/* free all arrays that were allocated by the init function
 */
void freeRNAWalkArrays(void);


/* the actual walking function that performs a walk on the
 *  structure landscape according to the defined method
 *  The returned secondary structure in dot bracket notation
 *  is the one where the walk ended...
 *  Implemented walks are defined below...
 */
#define GRADIENT_WALK   0
#define MC_METROPOLIS   1
char *structureWalk(const char  *seq,
                    const char  *structure,
                    int         method);


/* a simple position finding function that searches
 *  for the first entry in a sorted array of floats
 * that exceeds a given value
 */
int getPosition(float *array,
                float value,
                int   array_size);


#endif
