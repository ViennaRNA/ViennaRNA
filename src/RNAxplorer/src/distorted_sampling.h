#ifndef   _RNAXPLORER_DISTORTED_SAMPLING_H_
#define   _RNAXPLORER_DISTORTED_SAMPLING_H_

#include <ViennaRNA/data_structures.h>

typedef struct {
	int k;
	int l;
	int num_structs;
	int max_structs;
	char **structures;
	double pf;
	double mfe;
} gridpointT;

typedef struct {
	unsigned int size1;
	unsigned int size2;
	char *firstReference;
	char *secondReference;
	gridpointT **landscape;
} gridLandscapeT;

typedef struct {
	double kT;
	int *idx;
	char *ref1;
	char *ref2;
	double x;
	double y;
	unsigned int *referenceBPs1;
	unsigned int *referenceBPs2;
	short *reference_pt1;
	short *reference_pt2;
} kl_soft_constraints;

kl_soft_constraints *kl_init_datastructures(vrna_fold_compound_t *vc, const char *s1, const char *s2, double x,
		double y);

void free_kl_soft_constraints(void *data);

FLT_OR_DBL kl_pseudo_energy(int i, int j, int k, int l, unsigned char decomp, void *data);

FLT_OR_DBL kl_exp_pseudo_energy(int i, int j, int k, int l, unsigned char decomp, void *data);

/**
 * computes samples with stochastic backtracking and maps them to the 2D grid.
 * @param vc - the fold compound with modeldetails and maybe softconstraints
 * @param grid - the initialized 2D landscape
 * @param s1 - the first structure in dot-bracket notation
 * @param s2 - the second structure in dot-bracket notation
 * @param maxIterations - number that determines how many samples will be generated.
 */
void fillGridWithSamples(vrna_fold_compound_t *vc, gridLandscapeT *grid, const char *s1, const char *s2,
		int maxIterations);

/**
 * fill the grid by stepwise changing the potential of both references.
 * ref_x[i] = relaxfactor * ref_x[i] / maxSteps
 * @param vc - the fold compound with softconstraints (data = kl_soft_constraints)
 * @param grid - the initialized 2D landscape
 * @param relaxfactor - a factor
 * @param relax - if > 0: stepwise decrease the potential, else increase.
 * @param shift - shift the distortion to one reference in each step (and away from the other).
 * @param shift_to_first - shift the distortion to the first reference if it is set to 1, else to the second ref.
 * @param verbose - more output on stdout
 * @param maxIterations - maximal number of iterations for the stochastic sampling method
 * @param maxSteps - max. steps for computing the stepwise distortion.
 */
void fillGridStepwiseBothRef(vrna_fold_compound_t *vc, gridLandscapeT *grid, float relaxFactor, int relax, int shift,
		int shift_to_first, int verbose, int maxIterations, int maxSteps);

/**
 * fill the grid by stepwise changing the potential of the first reference.
 * ref_x[i] = relaxfactor * ref_x[i] / maxSteps
 * @param vc - the fold compound with softconstraints (data = kl_soft_constraints)
 * @param grid - the initialized 2D landscape
 * @param relaxfactor - a factor
 * @param relax - if > 0: stepwise decrease the potential, else increase.
 * @param shift - shift the distortion to one reference in each step (and away from the other).
 * @param shift_to_first - shift the distortion to the first reference if it is set to 1, else to the second ref.
 * @param verbose - more output on stdout
 * @param maxIterations - maximal number of iterations for the stochastic sampling method
 * @param maxSteps - max. steps for computing the stepwise distortion.
 */
void fillGridStepwiseFirstRef(vrna_fold_compound_t *vc, gridLandscapeT *grid, float relaxFactor, int relax, int verbose,
		int maxIterations, int maxSteps);

/**
 * fill the grid by stepwise changing the potential of the second reference.
 * ref_x[i] = relaxfactor * ref_x[i] / maxSteps
 * @param vc - the fold compound with softconstraints (data = kl_soft_constraints)
 * @param grid - the initialized 2D landscape
 * @param relaxfactor - a factor
 * @param relax - if > 0: stepwise decrease the potential, else increase.
 * @param verbose - more output on stdout
 * @param maxIterations - maximal number of iterations for the stochastic sampling method
 * @param maxSteps - max. steps for computing the stepwise distortion.
 */
void fillGridStepwiseSecondRef(vrna_fold_compound_t *vc, gridLandscapeT *grid, float relaxFactor, int relax,
		int verbose, int maxIterations, int maxSteps);

/**
 * computes and allocates maximal space for the 2D landscape. The size depends on the maximal basepairs that can be formed on sequence 1 and 2.
 * @param s - the RNA sequence
 * @param s1 - first dot-bracket structure
 * @param se1 - second dot-bracket structure
 * @return pointer to the 2D array
 */
gridLandscapeT *initLandscape(const char *s, const char *s1, const char *s2);


/**
 * computes the distortion for both references. The probability of both references is set to the probability of structure s0.
 * @param vc - the foldcompound which contains energy parameters (model details).
 * @param s0 - the structure which is the center of the distortion.
 * @param s1 - the first structure in dot-bracket notation
 * @param s2 - the second structure in dot-bracket notation
 * @param distortion_x - pointer to one value which receives the output
 * @param distortion_y - pointer to one value which receives the output
 */
void computeDistortion(vrna_fold_compound_t *vc, const char *s0, const char *s1, const char *s2,
		double *dist_x, double *dist_y);

/**
 * computes the distortion for both references. The probability of both references is set to the probability of the mfe.
 * @param vc - the foldcompound which contains energy parameters (model details).
 * @param s1 - the first structure in dot-bracket notation
 * @param s2 - the second structure in dot-bracket notation
 * @param distortion_x - pointer to one value which receives the output
 * @param distortion_y - pointer to one value which receives the output
 */
void computeInitialDistortion(vrna_fold_compound_t *vc, const char *s1, const char *s2, double *distortion_x,
		double *distortion_y);

/**
 * add the distortion softconstraints (data and functionpointer)to the foldcompound.
 * @param vc - the foldcompound which contains energy parameters (model details).
 * @param s1 - the first structure in dot-bracket notation
 * @param s2 - the second structure in dot-bracket notation
 * @param distortion_x - pointer to one value which receives the output
 * @param distortion_y - pointer to one value which receives the output
 */
void addSoftconstraints(vrna_fold_compound_t *vc, const char *s1, const char *s2, double distortion_x,
		double distortion_y);

/**
 * generate a stochastically sampled 2D map via distortion of the energy landscape
 * @param vc - the foldcompound which contains energy parameters. It will be filled with additional softconstraints.
 * @param s1 - the first structure in dot-bracket notation
 * @param s2 - the second structure in dot-bracket notation
 * @param maxIterations - max. iterations for stochastic sampling
 * @param extended_options - character string with distortion parameters
 * 													 B = alter both potentials at once
 * 													 R = relax potential instead of increasing it
 * 													 S = shift potential to other structure
 * 													 F = shift to first structure
 * 													 V = verbose
 */
gridLandscapeT*
estimate_landscape(vrna_fold_compound_t *vc, const char *s1, const char *s2, int maxIterations, char *extended_options);


void addStructure(gridLandscapeT *grid, char *structure);

/**
 * create model details with all options. Not the trimmed md that you would get from the current RNAlib interface.
 */
vrna_md_t createModelDetails(int circ, int uniq_ML, int compute_bpp, double betaScale);

void printLandscape(gridLandscapeT *grid, vrna_fold_compound_t *vc);

void free_gridLandscape(gridLandscapeT* grid);

#endif
