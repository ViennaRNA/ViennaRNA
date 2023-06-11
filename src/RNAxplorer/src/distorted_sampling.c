/*
 Secondary structure landscape sampling via distortion of the
 partition function

 (c) Ronny Lorenz
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <unistd.h>
#include <float.h>

#include <ViennaRNA/utils/basic.h>
#include <ViennaRNA/utils/structures.h>
#include <ViennaRNA/constraints/basic.h>
#include <ViennaRNA/constraints/soft.h>
#include <ViennaRNA/mfe.h>
#include <ViennaRNA/eval.h>
#include <ViennaRNA/part_func.h>
#include <ViennaRNA/mm.h>
#include <ViennaRNA/pair_mat.h>

#include "distorted_sampling.h"

kl_soft_constraints *kl_init_datastructures(vrna_fold_compound_t *vc, const char *s1, const char *s2, double x,
		double y) {
	kl_soft_constraints *data;
	unsigned int n;
	char *s = vc->sequence;

	n = strlen(s);

	/* alloc all memory */
	data = (kl_soft_constraints *) vrna_alloc(sizeof(kl_soft_constraints));
	data->kT = vc->exp_params->kT;
	data->idx = vrna_idx_row_wise(n);
	data->ref1 = strdup(s1);
	data->ref2 = strdup(s2);
	data->reference_pt1 = vrna_ptable(data->ref1);
	data->reference_pt2 = vrna_ptable(data->ref2);
	data->referenceBPs1 = vrna_refBPcnt_matrix(data->reference_pt1, TURN); /* matrix containing number of basepairs of reference structure1 in interval [i,j] */
	data->referenceBPs2 = vrna_refBPcnt_matrix(data->reference_pt2, TURN); /* matrix containing number of basepairs of reference structure2 in interval [i,j] */
	data->x = x;
	data->y = y;

	return data;
}

void free_kl_soft_constraints(void *data) {
	kl_soft_constraints *dat = (kl_soft_constraints *) data;

	free(dat->idx);
	free(dat->ref1);
	free(dat->ref2);
	free(dat->reference_pt1);
	free(dat->reference_pt2);
	free(dat->referenceBPs1);
	free(dat->referenceBPs2);
	free(dat);
}

FLT_OR_DBL kl_pseudo_energy(int i, int j, int k, int l, unsigned char decomp, void *data) {

	int d1, d2, ij, kl;
	int *idx = ((kl_soft_constraints *) data)->idx;
	short *reference_pt1 = ((kl_soft_constraints *) data)->reference_pt1;
	short *reference_pt2 = ((kl_soft_constraints *) data)->reference_pt2;
	unsigned int *referenceBPs1 = ((kl_soft_constraints *) data)->referenceBPs1;
	unsigned int *referenceBPs2 = ((kl_soft_constraints *) data)->referenceBPs2;
	double x = ((kl_soft_constraints *) data)->x;
	double y = ((kl_soft_constraints *) data)->y;

	int base_da = (reference_pt1[i] != (unsigned int) j) ? 1 : -1;
	int base_db = (reference_pt2[i] != (unsigned int) j) ? 1 : -1;
	ij = idx[i] - j;
	kl = idx[k] - l;
	d1 = d2 = 0;

	switch (decomp) {
	/* cases where we actually introduce a base pair */

	case VRNA_DECOMP_PAIR_HP:
		d1 = base_da + referenceBPs1[ij];
		d2 = base_db + referenceBPs2[ij];
		break;

	case VRNA_DECOMP_PAIR_IL:
		d1 = base_da + referenceBPs1[ij] - referenceBPs1[kl];
		d2 = base_db + referenceBPs2[ij] - referenceBPs2[kl];
		break;

	case VRNA_DECOMP_PAIR_ML:
		d1 = base_da + referenceBPs1[ij] - referenceBPs1[kl];
		d2 = base_db + referenceBPs2[ij] - referenceBPs2[kl];
		break;

		/* cases where we split a segment into one or two subsegments */

	case VRNA_DECOMP_ML_STEM:
		d1 = referenceBPs1[ij] - referenceBPs1[kl];
		d2 = referenceBPs2[ij] - referenceBPs2[kl];
		break;

	case VRNA_DECOMP_ML_ML:
		d1 = referenceBPs1[ij] - referenceBPs1[kl];
		d2 = referenceBPs2[ij] - referenceBPs2[kl];
		break;

	case VRNA_DECOMP_ML_ML_ML:
		d1 = referenceBPs1[ij] - referenceBPs1[idx[i] - k] - referenceBPs1[idx[l] - j];
		d2 = referenceBPs2[ij] - referenceBPs2[idx[i] - k] - referenceBPs2[idx[l] - j];
		break;

	case VRNA_DECOMP_ML_UP:
		d1 = referenceBPs1[ij];
		d2 = referenceBPs2[ij];
		break;

	case VRNA_DECOMP_ML_ML_STEM:
		d1 = referenceBPs1[ij] - referenceBPs1[idx[i] - k] - referenceBPs1[idx[l] - j];
		d2 = referenceBPs2[ij] - referenceBPs2[idx[i] - k] - referenceBPs2[idx[l] - j];
		break;

	case VRNA_DECOMP_ML_COAXIAL: /* (i,j) stacks onto (k,l), lets ignore this case for now */
		d1 = 0;
		d2 = 0;
		break;

	case VRNA_DECOMP_EXT_EXT:
		d1 = referenceBPs1[ij] - referenceBPs1[kl];
		d2 = referenceBPs2[ij] - referenceBPs2[kl];
		break;

	case VRNA_DECOMP_EXT_UP:
		d1 = referenceBPs1[ij];
		d2 = referenceBPs2[ij];
		break;

	case VRNA_DECOMP_EXT_EXT_EXT:
		d1 = referenceBPs1[ij] - referenceBPs1[idx[i] - k] - referenceBPs1[idx[l] - j];
		d2 = referenceBPs2[ij] - referenceBPs2[idx[i] - k] - referenceBPs2[idx[l] - j];
		break;

	case VRNA_DECOMP_EXT_STEM:
		d1 = referenceBPs1[ij] - referenceBPs1[kl];
		d2 = referenceBPs2[ij] - referenceBPs2[kl];
		break;

	case VRNA_DECOMP_EXT_EXT_STEM: /* fall through */
	case VRNA_DECOMP_EXT_STEM_EXT:
		d1 = referenceBPs1[ij] - referenceBPs1[idx[i] - k] - referenceBPs1[idx[l] - j];
		d2 = referenceBPs2[ij] - referenceBPs2[idx[i] - k] - referenceBPs2[idx[l] - j];
		break;

	case VRNA_DECOMP_EXT_EXT_STEM1:
		d1 = referenceBPs1[ij] - referenceBPs1[idx[i] - k] - referenceBPs1[idx[l] - (j - 1)];
		d2 = referenceBPs2[ij] - referenceBPs2[idx[i] - k] - referenceBPs2[idx[l] - (j - 1)];
		break;

	case VRNA_DECOMP_EXT_STEM_OUTSIDE:
		d1 = referenceBPs1[ij] - referenceBPs1[kl];
		d2 = referenceBPs2[ij] - referenceBPs2[kl];
		if (k > i) {
			d1 -= referenceBPs1[idx[i] - (k - 1)];
			d2 -= referenceBPs2[idx[i] - (k - 1)];
		}
		if (l < j) {
			d1 -= referenceBPs1[idx[l + 1] - j];
			d2 -= referenceBPs2[idx[l + 1] - j];
		}
		break;

	default:
		d1 = d2 = 0;
		break;
	}

	return (x * d1 + y * d2) * 100;
}

FLT_OR_DBL kl_exp_pseudo_energy(int i, int j, int k, int l, unsigned char decomp, void *data) {

	double kT = ((kl_soft_constraints *) data)->kT;
	double result = exp((-10. * (double) kl_pseudo_energy(i, j, k, l, decomp, data)) / kT);
	return result;
}

void fillGridWithSamples(vrna_fold_compound_t *vc, gridLandscapeT *grid, const char *s1, const char *s2,
		int maxIterations) {
	gridpointT **landscape = grid->landscape;
	vrna_init_rand();
	(void) vrna_pf(vc, NULL);
	for (int i = 0; i < maxIterations; i++) {
		char *sample = vrna_pbacktrack(vc);
		/* get k,l coords and free energy */
		int k = vrna_bp_distance(s1, sample);
		int l = vrna_bp_distance(s2, sample);
		double fe = (double) vrna_eval_structure(vc, sample);

		/* check if we have sufficient memory allocated and alloc more if necessary */
		if (landscape[k][l].num_structs + 2 >= landscape[k][l].max_structs) {
			landscape[k][l].max_structs *= 2;
			landscape[k][l].structures = (char **) vrna_realloc(landscape[k][l].structures,
					sizeof(char *) * landscape[k][l].max_structs);
		}

		/* insert structure */
		landscape[k][l].structures[landscape[k][l].num_structs] = sample;
		landscape[k][l].num_structs++;
		if (landscape[k][l].mfe > fe)
			landscape[k][l].mfe = fe;
	}
}

void fillGridStepwiseBothRef(vrna_fold_compound_t *vc, gridLandscapeT *grid, float relaxFactor, int relax, int shift,
		int shift_to_first, int verbose, int maxIterations, int maxSteps) {
	if (maxSteps <= 0) {
		fprintf(stderr, "Error: the stepsize has to be positive and greater than zero!");
		return;
	}

	kl_soft_constraints *data = (kl_soft_constraints*) vc->sc->data;
	char *s1 = data->ref1;
	char *s2 = data->ref2;
	double orig_x = data->x;
	double orig_y = data->y;

	for (int j = 0; j < maxSteps; j++) {
		if (shift) {
			if (shift_to_first) {
				data->x += relaxFactor * orig_x / maxSteps;
				data->y -= relaxFactor * orig_y / maxSteps;
			}
			else {
				data->x -= relaxFactor * orig_x / maxSteps;
				data->y += relaxFactor * orig_y / maxSteps;
			}
		}
		else {
			if (relax) {
				data->x -= relaxFactor * orig_x / maxSteps;
				data->y -= relaxFactor * orig_y / maxSteps;
			}
			else {
				data->x += relaxFactor * orig_x / maxSteps;
				data->y += relaxFactor * orig_y / maxSteps;
			}
		}

		if (verbose) {
			fprintf(stderr, "d_x = %1.10f, d_y = %1.10f\n", data->x, data->y);
		}
		fillGridWithSamples(vc, grid, s1, s2, maxIterations);
	}
	//reset to initial distortion
	data->x = orig_x;
	data->y = orig_y;
}

void fillGridStepwiseFirstRef(vrna_fold_compound_t *vc, gridLandscapeT *grid, float relaxFactor, int relax, int verbose,
		int maxIterations, int maxSteps) {
	if (maxSteps <= 0) {
		fprintf(stderr, "Error: the stepsize has to be positive and greater than zero!");
		return;
	}

	kl_soft_constraints *data = (kl_soft_constraints*) vc->sc->data;
	char *s1 = data->ref1;
	char *s2 = data->ref2;
	double orig_x = data->x;

	for (int j = 0; j < maxSteps; j++) {
		if (relax) {
			data->x -= relaxFactor * orig_x / maxSteps;
		}
		else {
			data->x += relaxFactor * orig_x / maxSteps;
		}

		if (verbose) {
			fprintf(stderr, "d_x = %1.10f, d_y = %1.10f\n", data->x, data->y);
		}
		fillGridWithSamples(vc, grid, s1, s2, maxIterations);
	}
	//reset to initial distortion
	data->x = orig_x;
}

void fillGridStepwiseSecondRef(vrna_fold_compound_t *vc, gridLandscapeT *grid, float relaxFactor, int relax,
		int verbose, int maxIterations, int maxSteps) {
	if (maxSteps <= 0) {
		fprintf(stderr, "Error: the stepsize has to be positive and greater than zero!");
		return;
	}

	kl_soft_constraints *data = (kl_soft_constraints*) vc->sc->data;
	char *s1 = data->ref1;
	char *s2 = data->ref2;
	double orig_y = data->y;
	for (int j = 0; j < maxSteps; j++) {
		if (relax) {
			data->y -= relaxFactor * orig_y / maxSteps;
		}
		else {
			data->y += relaxFactor * orig_y / maxSteps;
		}
		if (verbose) {
			fprintf(stderr, "d_x = %1.10f, d_y = %1.10f\n", data->x, data->y);
		}
		fillGridWithSamples(vc, grid, s1, s2, maxIterations);
	}
	//reset to initial distortion
	data->y = orig_y;
}

gridLandscapeT *initLandscape(const char *s, const char *s1, const char *s2) {
	int length = strlen(s);
	short *encodedString = encode_sequence(s, 0);
	int *iindx = vrna_idx_row_wise((unsigned) encodedString[0]);
	int idx_1n = iindx[1] - length;
	free(iindx);
	free(encodedString);

	short *pt_ref1 = vrna_ptable(s1);
	short *pt_ref2 = vrna_ptable(s2);
	/* then compute maximum matching with disallowed reference structures */
	unsigned int *mm1 = maximumMatchingConstraint(s, pt_ref1);
	unsigned int *mm2 = maximumMatchingConstraint(s, pt_ref2);

	/* get number of bp in reference structures */
	unsigned int i, bp_ref1, bp_ref2;
	for (bp_ref1 = bp_ref2 = 0, i = 1; i < length; i++) {
		if (pt_ref1[i] > i)
			bp_ref1++;
		if (pt_ref2[i] > i)
			bp_ref2++;
	}

	/* compute maximum d1 and maximum d2 */
	unsigned int MAX_k = bp_ref1 + mm1[idx_1n];
	unsigned int MAX_l = bp_ref2 + mm2[idx_1n];
	free(mm1);
	free(mm2);
	free(pt_ref1);
	free(pt_ref2);

	/* create the 2D landscape data structure */
	gridLandscapeT *grid = (gridLandscapeT*) vrna_alloc(sizeof(gridLandscapeT));
	grid->size1 = MAX_k + 1;
	grid->size2 = MAX_l + 1;
	gridpointT **landscape;
	landscape = (gridpointT **) vrna_alloc(sizeof(gridpointT) * (grid->size1));
	for (int i = 0; i < grid->size1; i++)
		landscape[i] = (gridpointT *) vrna_alloc(sizeof(gridpointT) * (grid->size2));

	/* alloc memory for 1000 structures per partition in the landscape */
	for (int i = 0; i < grid->size1; i++)
		for (int j = 0; j < grid->size2; j++) {
			landscape[i][j].max_structs = 1000;
			landscape[i][j].structures = (char **) vrna_alloc(sizeof(char *) * landscape[i][j].max_structs);
		}
	grid->landscape = landscape;
	grid->firstReference = (char*) vrna_alloc(sizeof(char)*(strlen(s1)+1));
	grid->secondReference = (char*) vrna_alloc(sizeof(char)*(strlen(s2)+1));
	strcpy(grid->firstReference, s1);
	strcpy(grid->firstReference, s2);
	return grid;
}

void computeDistortion(vrna_fold_compound_t *vc, const char *s0, const char *s1, const char *s2, double *dist_x,
		double *dist_y) {
	double distortion_x;
	double distortion_y;
	double s0fe = (double) vrna_mfe(vc, (char*) s0);

	/* get free energies of the reference structures */
	float e_ref1 = vrna_eval_structure(vc, s1);
	float e_ref2 = vrna_eval_structure(vc, s2);

	/* get base pair distance between both references */
	int bp_dist = vrna_bp_distance(s1, s2);
	int bp_dist_mfe_ref1 = vrna_bp_distance(s0, s1);
	int bp_dist_mfe_ref2 = vrna_bp_distance(s0, s2);

	//if(mmfe == e_ref1)
	if (bp_dist_mfe_ref1 == 0) {
		/*
		 distortion_x = 0;
		 distortion_y = distortion_x - (e_ref1 - e_ref2) / bp_dist;
		 we use the Wolfram alpha solution below ;)
		 */
		distortion_x = 0;
		if (bp_dist_mfe_ref2 != 0) {
			distortion_y = (distortion_x * bp_dist_mfe_ref2 - s0fe + e_ref2) / bp_dist_mfe_ref2;
		}
	}

	//if(mmfe == e_ref2)
	if (bp_dist_mfe_ref2 == 0) {
		/*
		 distortion_y = 0;
		 distortion_x = (e_ref1 - e_ref2) / bp_dist + distortion_y;
		 we use the Wolfram alpha solution below ;)
		 */
		distortion_y = 0;
		if (bp_dist_mfe_ref1 != 0) {
			distortion_x = (distortion_y * bp_dist_mfe_ref1 + e_ref1 - s0fe) / bp_dist_mfe_ref1;
		}
	}

	//GE: what is if s1=s2 and s1!=mfe and s2!=mfe ?

	//else
	if (bp_dist != 0 && (bp_dist_mfe_ref1 + bp_dist_mfe_ref2 != bp_dist)) {
		/*
		 distortion_x = ((e_ref1 * bp_dist_mfe_ref2) / (bp_dist * bp_dist_mfe_ref1)) - (mmfe / bp_dist_mfe_ref1);
		 distortion_y = ((e_ref2 * bp_dist_mfe_ref1) / (bp_dist * bp_dist_mfe_ref2)) - (mmfe / bp_dist_mfe_ref2);
		 we use the Wolfram alpha solution below ;)
		 */
		double nenner = bp_dist * (bp_dist_mfe_ref1 + bp_dist_mfe_ref2 - bp_dist);
		distortion_x = (bp_dist_mfe_ref2 * e_ref1 - bp_dist_mfe_ref2 * e_ref2 - bp_dist * s0fe + bp_dist * e_ref2) / nenner;
		distortion_y = (-bp_dist_mfe_ref1 * e_ref1 + bp_dist_mfe_ref1 * e_ref2 - bp_dist * s0fe + bp_dist * e_ref1)
				/ nenner;
	}

	*dist_x = distortion_x;
	*dist_y = distortion_y;
	printf("d_x = %1.10f, d_y = %1.10f\n", distortion_x, distortion_y);
}

void computeInitialDistortion(vrna_fold_compound_t *vc, const char *s1, const char *s2, double *dist_x, double *dist_y) {
	char * mfe_struct = (char *) vrna_alloc(sizeof(char) * (vc->length + 1));
	computeDistortion(vc, mfe_struct, s1, s2, dist_x, dist_y);
	free(mfe_struct);
}

void addSoftconstraints(vrna_fold_compound_t *vc, const char *s1, const char *s2, double distortion_x,
		double distortion_y) {
	vrna_sc_init(vc); // to remove old soft constraints
	kl_soft_constraints *data = kl_init_datastructures(vc, s1, s2, distortion_x, distortion_y);

	vrna_sc_add_data(vc, (void *) data, &free_kl_soft_constraints);
	vrna_sc_add_exp_f(vc, &kl_exp_pseudo_energy);
}

gridLandscapeT*
estimate_landscape(vrna_fold_compound_t *vc, const char *s1, const char *s2, int maxIterations, char *extended_options) {
	/* parse extended options string */
	int plain, both_at_once, relax, verbose, shift, shift_to_first;

	both_at_once = 0;
	relax = 0;
	plain = 1; /* plain sampling, no distortion */
	verbose = 0;
	shift = 0;
	shift_to_first = 0;

	if (extended_options) {
		plain = 0;
		if (strchr(extended_options, 'B')) /* alter both potentials at once */
			both_at_once = 1;
		if (strchr(extended_options, 'R')) /* relax potential instead of increasing it */
			relax = 1;
		if (strchr(extended_options, 'S')) /* shift potential to other structure */
			shift = 1;
		if (strchr(extended_options, 'F')) /* shift to first structure */
			shift_to_first = 1;
		if (strchr(extended_options, 'V')) /* verbose */
			verbose = 1;
	}

	/* get mfe for this sequence */
	char * mfe_struct = (char *) vrna_alloc(sizeof(char) * (vc->length + 1));
	double mmfe = (double) vrna_mfe(vc, mfe_struct);

	char *s = vc->sequence;
	gridLandscapeT* grid = initLandscape(s, s1, s2);

	if (plain) {
		/* prepare pf fold */
		vrna_exp_params_rescale(vc, &mmfe);
		//backtracking with energies without distortion
		fillGridWithSamples(vc, grid, s1, s2, maxIterations);
	}
	else {
		double distortion_x, distortion_y;
		computeInitialDistortion(vc, s1, s2, &distortion_x, &distortion_y);
		/* prepare pf fold */
		int bp_dist_mfe_ref1 = vrna_bp_distance(mfe_struct, s1);
		int bp_dist_mfe_ref2 = vrna_bp_distance(mfe_struct, s2);

		double rescale = mmfe + (bp_dist_mfe_ref1 * distortion_x) + (bp_dist_mfe_ref2 * distortion_y);
		vrna_exp_params_rescale(vc, &rescale);

		/* apply distortion soft constraints */
		addSoftconstraints(vc, s1, s2, distortion_x, distortion_y);

		float relaxFactor = 1.5; /* if set to 1. final relaxation sets x and/or y to 0. */
		int bp_dist = vrna_bp_distance(s1, s2);
		if (!both_at_once) {
			/* change potential of first reference */
			fillGridStepwiseFirstRef(vc, grid, relaxFactor, relax, verbose, maxIterations, bp_dist);
			/* change potential of second reference */
			fillGridStepwiseSecondRef(vc, grid, relaxFactor, relax, verbose, maxIterations, bp_dist);

		}
		else { /* change potential of both references at the same time */
			fillGridStepwiseBothRef(vc, grid, relaxFactor, relax, shift, shift_to_first, verbose, maxIterations, bp_dist);
		}
	}

	free(mfe_struct);

	// evaluate structures, set properties.
	// compute undistorted energies.
	vrna_sc_free(vc->sc);
	vc->sc = NULL;
	vc->params->model_details.betaScale = 1;
	vrna_exp_params_rescale(vc, &mmfe);
	float mfe = FLT_MAX;
	float tmpMFE = 0;
	for (int i = 0; i < grid->size1; i++) {
		for (int j = 0; j < grid->size2; j++) {
			if (grid->landscape[i][j].num_structs > 0) {
				int k;
				for (k = 0; k < grid->landscape[i][j].num_structs; k++) {
					tmpMFE = vrna_eval_structure(vc, grid->landscape[i][j].structures[k]);
					if (tmpMFE < mfe) {
						grid->landscape[i][j].mfe = mfe;
					}
				}
				grid->landscape[i][j].k = i;
				grid->landscape[i][j].l = j;
			}
		}
	}

	return grid;
}

void addStructure(gridLandscapeT *grid, char *structure) {
	unsigned int k = vrna_bp_distance(grid->firstReference, structure);
	unsigned int l = vrna_bp_distance(grid->secondReference, structure);

	if (k <= grid->size1 && l <= grid->size2) {
		char *newStruct = (char*) vrna_alloc(strlen(structure));
		strcpy(newStruct, structure);
		gridpointT **landscape = grid->landscape;
		if (landscape[k][l].num_structs + 2 >= landscape[k][l].max_structs) {
			landscape[k][l].max_structs *= 2;
			landscape[k][l].structures = (char **) vrna_realloc(landscape[k][l].structures,
					sizeof(char *) * landscape[k][l].max_structs);
		}

		/* insert structure */
		landscape[k][l].structures[landscape[k][l].num_structs] = newStruct;
		landscape[k][l].num_structs++;
	}
	else {
		fprintf(stderr, "Error: the structure %s does not belong to the grid!", structure);
	}
}

vrna_md_t createModelDetails(int circ, int uniq_ML, int compute_bpp, double betaScale) {
	vrna_md_t md;
	vrna_md_set_default(&md);
	md.circ = circ;
	md.uniq_ML = uniq_ML; /* in case we need M1 arrays */
	md.compute_bpp = compute_bpp;
	md.betaScale = betaScale;
	return md;
}

void printLandscape(gridLandscapeT *grid, vrna_fold_compound_t *vc) {
	for (int i = 0; i < grid->size1; i++) {
		for (int j = 0; j < grid->size2; j++) {
			if (grid->landscape[i][j].num_structs > 0) {
				int k;
				for (k = 0; k < grid->landscape[i][j].num_structs; k++)
					printf("%d\t%d\t%6.2f\t(%d) %s\n", i, j, vrna_eval_structure(vc, grid->landscape[i][j].structures[k]),
							grid->landscape[i][j].num_structs, grid->landscape[i][j].structures[k]);
			}
		}
	}
}

void free_gridLandscape(gridLandscapeT *grid) {
	for (int i = 0; i < grid->size1; i++) {
		for (int j = 0; j < grid->size2; j++) {

			for (int k = 0; k < grid->landscape[i][j].num_structs; k++) {
				free(grid->landscape[i][j].structures[k]);
			}
			free(grid->landscape[i][j].structures);
		}
		free(grid->landscape[i]);
	}
	free(grid->landscape);
	free(grid->firstReference);
	free(grid->secondReference);
	free(grid);
}

