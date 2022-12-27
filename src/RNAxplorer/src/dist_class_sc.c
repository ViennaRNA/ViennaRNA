/*
 * provide means for multi-dimensional distance class based soft constraints
 */

#include <stdlib.h>
#include <math.h>
#include <string.h>

#include <ViennaRNA/utils/basic.h>
#include <ViennaRNA/constraints/soft.h>

#include "dist_class_sc.h"


static int *
compute_distancies(int              i,
                   int              j,
                   int              k,
                   int              l,
                   unsigned char    decomp,
                   sc_dist_class_t  *d)
{
  unsigned int  r, numberOfRefs, **referenceBPs;
  int           ij, kl, *base_dx, *dist, *idx;
  short         **references_pt;

  references_pt = d->ref_pts;
  numberOfRefs  = d->ref_num;
  referenceBPs  = d->ref_bps;
  idx           = d->idx;

  ij  = idx[i] - j;
  kl  = idx[k] - l;

  base_dx = (int *)vrna_alloc(sizeof(int) * (numberOfRefs + 1));
  dist    = (int *)vrna_alloc(sizeof(int) * (numberOfRefs + 1));

  for (r = 0; r < numberOfRefs; r++)
    base_dx[r] = ((int)references_pt[r][i] != j) ? 1 : -1;

  switch (decomp) {
    /* cases where we actually introduce a base pair */
    case VRNA_DECOMP_PAIR_HP:
      for (r = 0; r < numberOfRefs; r++)
        dist[r] = base_dx[r] +
                  referenceBPs[r][ij];
      break;

    case VRNA_DECOMP_PAIR_IL:
      for (r = 0; r < numberOfRefs; r++)
        dist[r] = base_dx[r] +
                  referenceBPs[r][ij] -
                  referenceBPs[r][kl];
      break;

    case VRNA_DECOMP_PAIR_ML:
      for (r = 0; r < numberOfRefs; r++)
        dist[r] = base_dx[r] +
                  referenceBPs[r][ij] -
                  referenceBPs[r][kl];
      break;

    /* cases where we split a segment into one or two subsegments */

    case VRNA_DECOMP_ML_STEM:
      for (r = 0; r < numberOfRefs; r++)
        dist[r] = referenceBPs[r][ij] -
                  referenceBPs[r][kl];
      break;

    case VRNA_DECOMP_ML_ML:
      for (r = 0; r < numberOfRefs; r++)
        dist[r] = referenceBPs[r][ij] -
                  referenceBPs[r][kl];
      break;

    case VRNA_DECOMP_ML_ML_ML:
      for (r = 0; r < numberOfRefs; r++)
        dist[r] = referenceBPs[r][ij] -
                  referenceBPs[r][idx[i] - k] -
                  referenceBPs[r][idx[l] - j];
      break;

    case VRNA_DECOMP_ML_UP:
      for (r = 0; r < numberOfRefs; r++)
        dist[r] = referenceBPs[r][ij];
      break;

    case VRNA_DECOMP_ML_ML_STEM:
      for (r = 0; r < numberOfRefs; r++)
        dist[r] = referenceBPs[r][ij] -
                  referenceBPs[r][idx[i] - k] -
                  referenceBPs[r][idx[l] - j];
      break;

    case VRNA_DECOMP_ML_COAXIAL: /* (i,j) stacks onto (k,l), lets ignore this case for now */
      for (r = 0; r < numberOfRefs; r++)
        dist[r] = 0;
      break;

    case VRNA_DECOMP_EXT_EXT:
      for (r = 0; r < numberOfRefs; r++)
        dist[r] = referenceBPs[r][ij] -
                  referenceBPs[r][kl];
      break;

    case VRNA_DECOMP_EXT_UP:
      for (r = 0; r < numberOfRefs; r++)
        dist[r] = referenceBPs[r][ij];
      break;

    case VRNA_DECOMP_EXT_EXT_EXT:
      for (r = 0; r < numberOfRefs; r++)
        dist[r] = referenceBPs[r][ij] -
                  referenceBPs[r][idx[i] - k] -
                  referenceBPs[r][idx[l] - j];
      break;

    case VRNA_DECOMP_EXT_STEM:
      for (r = 0; r < numberOfRefs; r++)
        dist[r] = referenceBPs[r][ij] -
                  referenceBPs[r][kl];
      break;

    case VRNA_DECOMP_EXT_EXT_STEM: /* fall through */
    case VRNA_DECOMP_EXT_STEM_EXT:
      for (r = 0; r < numberOfRefs; r++)
        dist[r] = referenceBPs[r][ij] -
                  referenceBPs[r][idx[i] - k] -
                  referenceBPs[r][idx[l] - j];
      break;

    case VRNA_DECOMP_EXT_EXT_STEM1:
      for (r = 0; r < numberOfRefs; r++)
        dist[r] = referenceBPs[r][ij] -
                  referenceBPs[r][idx[i] - k] -
                  referenceBPs[r][idx[l] - (j - 1)];
      break;

    case VRNA_DECOMP_EXT_STEM_OUTSIDE:
      for (r = 0; r < numberOfRefs; r++) {
        dist[r] = referenceBPs[r][ij] -
                   referenceBPs[r][kl];
        if (k > i)
          dist[r] -= referenceBPs[r][idx[i] - (k - 1)];

        if (l < j)
          dist[r] -= referenceBPs[r][idx[l + 1] - j];
      }
      break;

    default:
      fprintf(stderr, "default sc\n");
      for (r = 0; r < numberOfRefs; r++)
        dist[r] = 0;
      break;
  }

  free(base_dx);

  return dist;
}


FLT_OR_DBL
sc_exp_f_dist_class(int           i,
                    int           j,
                    int           k,
                    int           l,
                    unsigned char decomp,
                    void          *data)
{
  int             *distancies;
  double          kT;
  FLT_OR_DBL      result;
  sc_dist_class_t *d;

  kT  = ((sc_dist_class_t *)data)->kT;
  d   = ((sc_dist_class_t *)data);

  distancies  = compute_distancies(i, j, k, l, decomp, d);
  result      = (FLT_OR_DBL)exp(-10. * (d->f(i, j, k, l, decomp, distancies, d)) / kT);

  free(distancies);
  return result;
}


int
sc_f_dist_class(int           i,
                int           j,
                int           k,
                int           l,
                unsigned char decomp,
                void          *data)
{
  int             *distancies;
  int             result;
  sc_dist_class_t *d;

  d       = ((sc_dist_class_t *)data);
  result  = 0;

  distancies  = compute_distancies(i, j, k, l, decomp, d);
  result      = (int)d->f(i, j, k, l, decomp, distancies, d);

  free(distancies);

  return result;
}


sc_dist_class_t *
sc_dist_class_init(vrna_fold_compound_t *fc)
{
  sc_dist_class_t *d;

  d = (sc_dist_class_t *)vrna_alloc(sizeof(sc_dist_class_t));

  d->idx  = fc->iindx;
  d->kT   = fc->exp_params->kT;

  d->ref_num    = 0;
  d->references = NULL;
  d->ref_bps    = NULL;
  d->ref_pts    = NULL;

  d->f      = NULL;
  d->f_data = NULL;
  d->f_free = NULL;

  return d;
}


void
sc_dist_class_destroy(void *data)
{
  int             i;
  sc_dist_class_t *d = (sc_dist_class_t *)data;

  for (i = 0; i < d->ref_num; i++) {
    free(d->references[i]);
    free(d->ref_pts[i]);
    free(d->ref_bps[i]);
  }

  free(d->references);
  free(d->ref_bps);
  free(d->ref_pts);

  if (d->f_free)
    d->f_free(d->f_data);

  free(d);
}


void
sc_dist_class_add_ref(sc_dist_class_t *d,
                      const char      *ref_struct)
{
  int n = d->ref_num;

  d->references     = (char **)vrna_realloc(d->references, sizeof(char *) * (n + 1));
  d->references[n]  = strdup(ref_struct);
  d->ref_pts        = (short **)vrna_realloc(d->ref_pts, sizeof(short *) * (n + 1));
  d->ref_pts[n]     = vrna_ptable(ref_struct);
  d->ref_bps        = (unsigned int **)vrna_realloc(d->ref_bps, sizeof(unsigned int *) * (n + 1));
  d->ref_bps[n]     = vrna_refBPcnt_matrix(d->ref_pts[n], TURN);

  d->ref_num++;
}
