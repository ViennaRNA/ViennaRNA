
#include <stdlib.h>

#include <ViennaRNA/utils.h>
#include <ViennaRNA/structure_utils.h>
#include <ViennaRNA/mm.h>
#include <ViennaRNA/constraints_soft.h>
#include <ViennaRNA/eval.h>
#include <ViennaRNA/mfe.h>
#include <ViennaRNA/part_func.h>
#include <ViennaRNA/boltzmann_sampling.h>

#include "dist_class_sc.h"
#include "repellant_sampling.h"


typedef struct {
  unsigned int  num_ref;
  int           *repulsion;
  unsigned int  **mm1;
  unsigned int  **ref_bps;
} repell_data;


static void
destroy_repell_data(void *data);


static double
compute_repulsion(int i,
                  int j,
                  int k,
                  int l,
                  unsigned char decomp,
                  int *distance,
                  sc_dist_class_t *d);


static int
add_reference(vrna_fold_compound_t *fc,
              repell_data          *d,
              const char           *structure,
              double               energy);



static int
change_reference_energy(vrna_fold_compound_t *fc,
                        repell_data          *d,
                        int                  id,
                        double               energy);


int
rnax_add_repulsion(vrna_fold_compound_t *fc,
                   const char *structure,
                   double     strength)
{
  int ret = -1;

  if ((fc) && (fc->type == VRNA_FC_TYPE_SINGLE)) {
    if (!fc->sc)
      vrna_sc_init(fc);

    /* init everything if this is the first call to this function */
    if (fc->sc->exp_f != &sc_exp_f_dist_class) {
      sc_dist_class_t *d  = sc_dist_class_init(fc);
      repell_data     *rd = (repell_data *)vrna_alloc(sizeof(repell_data));

      d->f_data = (void *)rd;
      d->f_free = &destroy_repell_data;
      d->f      = &compute_repulsion;

      /* add soft constraints to fold compound */
      vrna_sc_add_data(fc, (void *)d, &sc_dist_class_destroy);
      vrna_sc_add_exp_f(fc, &sc_exp_f_dist_class);
      vrna_sc_add_f(fc, &sc_f_dist_class);
    }

    /* add reference */
    sc_dist_class_t *d = (sc_dist_class_t *)fc->sc->data;
    sc_dist_class_add_ref(d, structure);
    ret = add_reference(fc, (repell_data *)d->f_data, structure, strength);
  }

  return ret;
}


int
rnax_change_repulsion(vrna_fold_compound_t *fc,
                      int                  id,
                      double     strength)
{
  int ret = -1;

  if ((fc) && (fc->type == VRNA_FC_TYPE_SINGLE) && (fc->sc->data)) {
    sc_dist_class_t *d = (sc_dist_class_t *)fc->sc->data;
    ret = change_reference_energy(fc, (repell_data *)d->f_data, id, strength);
  }

  return ret;
}


/* BEGIN interface for repulsive sampling */
void
repellant_sampling(vrna_fold_compound_t *vc)
{
  unsigned int  n;
  vrna_md_t     md;

  vrna_md_set_default(&md);
  md.uniq_ML     = 1;
  md.compute_bpp = 0;

  vrna_fold_compound_t *fc = vrna_fold_compound(vc->sequence, &md, VRNA_OPTION_PF);

  n = fc->length;

  /* compute 'real' MFE structure */
  char    *mfe_structure = (char *)vrna_alloc(sizeof(char) * (n + 1));
  double  mfe            = (double)vrna_mfe(fc, mfe_structure);

  printf("%s [ %6.2f ]\n", mfe_structure, mfe);

  vrna_exp_params_rescale(fc, &mfe);

  /* init RNAlib random number seed */
  vrna_init_rand();

  /* add MFE struct as first structure to repell */
  rnax_add_repulsion(fc, mfe_structure, -mfe);

  double mfe2 = (double)vrna_mfe(fc, NULL);

  vrna_exp_params_rescale(fc, &mfe2);

  /* fill partition function DP matrices */
  (void)vrna_pf(fc, NULL);

  for(int i = 0; i < 100; i++) {
    char *s = vrna_pbacktrack(fc);
    printf("%s [ %6.2f ]\n", s, vrna_eval_structure_simple(fc->sequence, s));
    free(s);
  }
}


/* BEGIN static helper functions for repulsive sampling */
static void
destroy_repell_data(void *data)
{
  unsigned int i;

  repell_data *d = (repell_data *)data;
  free(d->repulsion);
  for (i = 0; i < d->num_ref; i++) {
    free(d->mm1[i]);
    free(d->ref_bps[i]);
  }

  free(d->mm1);
  free(d->ref_bps);

  free(d);
}

static double
compute_repulsion(int i,
                  int j,
                  int k,
                  int l,
                  unsigned char decomp,
                  int *distance,
                  sc_dist_class_t *d)
{
  double        rr, r     = 0.;
  int           idx       = d->idx[i] - j;
  repell_data   *dd       = (repell_data *)d->f_data;

  /* Idea: repulsion is maximal at distance 0, and decreases linearly until it vanishes at maximum distance from reference */

  switch (decomp) {
    case VRNA_DECOMP_PAIR_HP:
      for (int s = 0; s < d->ref_num; s++) {
        rr = dd->ref_bps[s][idx] +
             dd->mm1[s][idx];
        r += (double)(rr - distance[s]) * dd->repulsion[s];
      }
      break;

    case VRNA_DECOMP_PAIR_IL:
      for (int s = 0; s < d->ref_num; s++) {
        rr = dd->ref_bps[s][idx] +
             dd->mm1[s][idx] -
             dd->ref_bps[s][d->idx[k] - l] -
             dd->mm1[s][d->idx[k] - l];
        r += (double)(rr - distance[s]) * dd->repulsion[s];
      }
      break;

    case VRNA_DECOMP_PAIR_ML:
      for (int s = 0; s < d->ref_num; s++) {
        rr = dd->ref_bps[s][idx] +
             dd->mm1[s][idx] -
             dd->ref_bps[s][d->idx[k] - l] -
             dd->mm1[s][d->idx[k] - l];
        r += (double)(rr - distance[s]) * dd->repulsion[s];
      }
      break;

    case VRNA_DECOMP_ML_STEM: /* fall through */
    case VRNA_DECOMP_ML_ML:
      for (int s = 0; s < d->ref_num; s++) {
        rr = dd->ref_bps[s][idx] +
             dd->mm1[s][idx] -
             dd->ref_bps[s][d->idx[k] - l] -
             dd->mm1[s][d->idx[k] - l];
        r += (double)(rr - distance[s]) * dd->repulsion[s];
      }
      break;

    case VRNA_DECOMP_ML_ML_STEM: /* fall through */
    case VRNA_DECOMP_ML_ML_ML:
      for (int s = 0; s < d->ref_num; s++) {
        rr = dd->ref_bps[s][idx] +
             dd->mm1[s][idx] -
             dd->ref_bps[s][d->idx[i] - k] -
             dd->mm1[s][d->idx[i] - k] -
             dd->ref_bps[s][d->idx[l] - j] -
             dd->mm1[s][d->idx[l] - j];
        r += (double)(rr - distance[s]) * dd->repulsion[s];
      }
      break;

    case VRNA_DECOMP_ML_UP:
      for (int s = 0; s < d->ref_num; s++) {
        rr = dd->ref_bps[s][idx] +
             dd->mm1[s][idx];
        r += (double)(rr - distance[s]) * dd->repulsion[s];
      }
      break;

    case VRNA_DECOMP_ML_COAXIAL: /* (i,j) stacks onto (k,l), lets ignore this case for now */
      break;

    case VRNA_DECOMP_EXT_STEM: /* fall through */
    case VRNA_DECOMP_EXT_EXT:
      for (int s = 0; s < d->ref_num; s++) {
        rr = dd->ref_bps[s][idx] + /* This is how far we can be within this segment theoretically */
             dd->mm1[s][idx] -
             dd->ref_bps[s][d->idx[k] - l] -
             dd->mm1[s][d->idx[k] - l];
        r += (double)(rr - distance[s]) * dd->repulsion[s];
      }
      break;

    case VRNA_DECOMP_EXT_UP:
      for (int s = 0; s < d->ref_num; s++) {
        rr = dd->ref_bps[s][idx] + /* This is how far we can be within this segment theoretically */
             dd->mm1[s][idx];
        r += (double)(rr - distance[s]) * dd->repulsion[s];
      }
      break;

    case VRNA_DECOMP_EXT_EXT_EXT: /* fall through */
    case VRNA_DECOMP_EXT_STEM_EXT: /* fall through */
    case VRNA_DECOMP_EXT_EXT_STEM:
      for (int s = 0; s < d->ref_num; s++) {
        /* We start with the maximum amount of bps we can be apart within this segment (theoretically) */
        rr = dd->ref_bps[s][idx] + 
             dd->mm1[s][idx] -
            /* And substract how far we can still be with this decomposition */
             dd->ref_bps[s][d->idx[i] - k] -
             dd->mm1[s][d->idx[i] - k] -
             dd->ref_bps[s][d->idx[l] - j] -
             dd->mm1[s][d->idx[l] - j];
        r += (double)(rr - distance[s]) * dd->repulsion[s];
      }
      break;

    case VRNA_DECOMP_EXT_EXT_STEM1:
      break;

    default:
      vrna_message_warning("unhandled decomposition (%d) in compute_repulsion()", (int)decomp);
      break;
  }

  return r;
}

static int
add_reference(vrna_fold_compound_t *fc,
              repell_data          *d,
              const char           *structure,
              double               energy)
{
  short        *pt;
  unsigned int n;
  int          *iidx;
  vrna_md_t    *md;

  n    = fc->length;
  iidx = fc->iindx;
  md   = &(fc->params->model_details);
  pt   = vrna_ptable(structure);

  /* increase number of references */
  d->num_ref++;

  /* (re-)allocate required memory */
  d->repulsion  = (int *)vrna_realloc(d->repulsion, d->num_ref * sizeof(int));
  d->mm1        = (unsigned int **)vrna_realloc(d->mm1, d->num_ref * sizeof(unsigned int *));
  d->ref_bps    = (unsigned int **)vrna_realloc(d->ref_bps, d->num_ref * sizeof(unsigned int *));

  /* prepare reference specific stuff */
  d->mm1[d->num_ref - 1]       = maximumMatchingConstraint(fc->sequence, pt);
  d->ref_bps[d->num_ref - 1]   = vrna_refBPcnt_matrix(pt, md->min_loop_size);
  d->repulsion[d->num_ref - 1] = (int)((energy / (d->ref_bps[d->num_ref - 1][iidx[1] - n] + d->mm1[d->num_ref - 1][iidx[1] - n])) * 100.);/* in dekakal/mol, but what is a good repulsion???!!! */

  free(pt);

  return d->num_ref - 1;
}


static int
change_reference_energy(vrna_fold_compound_t *fc,
                        repell_data          *d,
                        int                  id,
                        double               energy)
{
  unsigned int n;
  int          *iidx;

  n    = fc->length;
  iidx = fc->iindx;

  d->repulsion[id] = (int)((energy / (d->ref_bps[id][iidx[1] - n] + d->mm1[id][iidx[1] - n])) * 100.);/* in dekakal/mol */

  return id;
}


