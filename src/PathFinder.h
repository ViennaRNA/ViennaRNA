#ifndef RNAXPLORER_PATH_FINDER_H
#define RNAXPLORER_PATH_FINDER_H

#include <ViennaRNA/model.h>
#include <ViennaRNA/findpath.h>

/* perform gradient walks from saddle point to find meshpoint(s) */
#define  RNAX_PATHFINDER_SADDLE_GRADIENT_WALK   1U

/* perform rejection-less monte carlo (Gillespie) simulation away from saddle point */
#define  RNAX_PATHFINDER_SADDLE_MONTE_CARLO     2U

/* perform rejection-less monte carlo (Gillespie) simulation away from saddle point while cooling down the system (simulated annealing) */
#define  RNAX_PATHFINDER_SADDLE_MONTE_CARLO_SA  3U  

/* use 2D representatives as meshpoints */
#define  RNAX_PATHFINDER_TWO_D_REPRESENTATIVES  4U


typedef struct {
  unsigned int method;
  unsigned int iterations;
  int max_keep;
  unsigned int storage_size;
  unsigned int max_paths;

  int         max_d1;
  int         max_d2;

  vrna_md_t   md;
} rnax_path_finder_opt_t;


rnax_path_finder_opt_t *
rnax_path_finder_options(void);


vrna_path_t *
rnax_path_finder( const char              *seq,
                  const char              *s_source,
                  const char              *s_target,
                  rnax_path_finder_opt_t  *options);


/**
 * !
 * @param seq - the RNA sequence (AGCU).
 * @param s1 - the first secondary structure in dot-bracket format.
 * @param s2 - the second secondary structure in dot-bracket format.
 * @param maxIterations -
 * @param maxKeep - maximum structures that will be kept (see findpath.h)
 * @param method - method for structurewalk (MC_METROPOLIS or GRADIENT_WALK) (see RNAwalk.h)
 * @param maxStorage - for insert_meshpoint
 */
void levelSaddlePoint(const char *seq, const char *s1, const char *s2, int maxIterations, int maxKeep, int method, int maxStorage);

vrna_path_t *levelSaddlePoint2(const char *seq, const char *s1, const char *s2/*, int *num_entry*/, int iteration,
		int maxIterations, int maxKeep, int maxStorage, int maximum_distance1, int maximum_distance2);
vrna_path_t *getSaddlePoint(vrna_path_t *foldingPath);

#endif /* RNAXPLORER_PATH_FINDER_H */
