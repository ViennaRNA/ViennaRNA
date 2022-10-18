/*
 * Secondary structure landscape sampling via distortion of the
 * partition function
 *
 * (c) Gregor Entzian
 */

#include "config.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <unistd.h>
#include <float.h>

/*
 #include <ViennaRNA/fold_vars.h>
 #include <ViennaRNA/energy_const.h>
 #include <ViennaRNA/findpath.h>
 #include <ViennaRNA/2Dpfold.h>
 */

#include <ViennaRNA/utils.h>
#include <ViennaRNA/structure_utils.h>
#include <ViennaRNA/constraints.h>
#include <ViennaRNA/mfe.h>
#include <ViennaRNA/eval.h>
#include <ViennaRNA/model.h>
#include <ViennaRNA/params.h>
#include <ViennaRNA/fold.h>
#include <ViennaRNA/part_func.h>
#include <ViennaRNA/mm.h>
#include <ViennaRNA/pair_mat.h>

#ifdef HAVE_LAPACKE_H
# include <lapacke.h>
#else
# ifdef HAVE_LAPACKE_LAPACKE_H
#   include <lapacke/lapacke.h>
# else
#   ifdef HAVE_OPENBLAS_LAPACKE_H
#     include <openblas/lapacke.h>
#   endif
# endif
#endif


#include "dist_class_sc.h"
#include "distorted_samplingMD.h"

#define max(x, y) (((x) > (y)) ? (x) : (y))
#define min(x, y) (((x) < (y)) ? (x) : (y))

#define WITH_DIST_CLASS_SC  1

#if WITH_DIST_CLASS_SC

typedef struct {
  size_t *maxDistances;
  double *distortions;
} kl_soft_constraints_MD;

#else
typedef struct {
  double kT;
  int *idx;
  int numberOfReferences;
  const char **references;
  double *distortions;
  /* matrix containing number of basepairs of reference structure i in interval [i,j] */
  unsigned int ** referencesBPs;
  short ** referencesAsPairtaibles;
  size_t *maxDistances; /* maximal distance to each reference. */
  int repel;
} kl_soft_constraints_MD;

static void free_kl_soft_constraints_MD(void *data);
static FLT_OR_DBL kl_pseudo_energy_MD(int i, int j, int k, int l, unsigned char decomp, void *data);
static FLT_OR_DBL kl_exp_pseudo_energy_MD(int i, int j, int k, int l, unsigned char decomp, void *data);
static kl_soft_constraints_MD *kl_init_datastructures_MD(vrna_fold_compound_t *vc, const char **referenceStructures,
    int numberOfReferenceStructures, double *distortions,int repel);

#endif

static void
fillGridStepwiseBothRef_MD(vrna_fold_compound_t *vc, gridLandscapeT *grid, float relaxFactor, int relax, int shift,
    int shift_to_first, int verbose, int maxIterations, int maxSteps);


unsigned int
getMaximalPossibleBPdistance(const char *sequence,
                             const char *structure)
{
  short         *pt_structure = vrna_ptable(structure);
  unsigned int  *mm1          = maximumMatchingConstraint(sequence, pt_structure);
  unsigned int  length        = strlen(sequence);
  int           *iindx        = vrna_idx_row_wise(length);
  int           idx_1n        = iindx[1] - length;
  /* get number of bp in structure */
  unsigned int  i, bp_structure;

  for (bp_structure = 0, i = 1; i < length; i++)
    if (pt_structure[i] > i)
      bp_structure++;

  /* compute maximum distance to this structure */
  unsigned int maxDistance = bp_structure + mm1[idx_1n];
  return maxDistance;
}


double *
rxp_computeDistortionsWRTMaxDistance(vrna_fold_compound_t *fc,
                                     const char           **structures,
                                     size_t               numberOfStructures,
                                     double               *maxDistances)
{
  char    *mfeStructure = (char *)vrna_alloc(sizeof(char) * (fc->length + 1));
  double  mfe           = (double)vrna_mfe(fc, mfeStructure);

  int     acceptedIndices[numberOfStructures];
  int     numberAccepted = 0;

  /* filter for unique input structures, that are not the mfe. (accepts only the last unique structure) */
  for (int i = 0; i < numberOfStructures; i++) {
    const char  *s1             = structures[i];
    int         acceptStructure = 1;

    if (strcmp(s1, mfeStructure) == 0) {
      acceptedIndices[i] = -1;
      continue;
    } else {
      for (int j = i + 1; j < numberOfStructures; j++) {
        const char *s2 = structures[j];
        if (strcmp(s1, s2) == 0) {
          acceptStructure = 0;
          break;
        }
      }
    }

    if (acceptStructure == 1) {
      acceptedIndices[i] = 1;
      numberAccepted++;
    } else {
      acceptedIndices[i] = -1;
    }
  }
  free(mfeStructure);
  if (numberAccepted < 1) {
    /* no distortion possible */
    double *distortions = vrna_alloc(sizeof(double) * numberOfStructures);
    for (int i = 0; i < numberOfStructures; i++)
      distortions[i] = 0;
    return distortions;
  }

  const char  *uniqueStructures[numberAccepted + 1]; /* +1 for the mfe */
  int         uniqueInd = 0;
  for (int i = 0; i < numberOfStructures; i++) {
    if (acceptedIndices[i] != -1) {
      uniqueStructures[uniqueInd] = structures[i];
      uniqueInd++;
    }
  }
  uniqueStructures[numberAccepted] = mfeStructure;

  float energies[numberAccepted + 1]; /* +1 for the mfe */
  energies[numberAccepted] = mfe;

  for (int i = 0; i < numberAccepted; i++) {
    const char *structure = uniqueStructures[i];
    energies[i] = vrna_eval_structure(fc, structure);
  }

  /* structures + mfeStr */
  /* energies + mfe */
  size_t totalStructures = numberAccepted + 1;  /*  + 1 (with and without mfe) */
  /**
   * construct equations, which look like:
   * S1: E(s1) + x'*d(s1, s1) + y'*d(s1, s2) + z'*d(s1, s3)
   * S2: E(s2) + x'*d(s2, s1) + y'*d(s2, s2) + z'*d(s2, s3)
   * ...
   * S4: E(mfeStructure) + ...
   */

  /* first dim number of equations and the second dim are the values energy(s),x',y',... */
  float   equations[totalStructures][totalStructures];
  size_t  distance;
  for (int i = 0; i < totalStructures; i++) {
    equations[i][0] = energies[i];
    for (int j = 0; j < (totalStructures - 1); j++) {
      const char  *str1 = uniqueStructures[i];
      const char  *str2 = uniqueStructures[j];
      distance            = bp_distance(str1, str2);
      equations[i][j + 1] = maxDistances[j] - distance;
    }
  }
  /* for(int j = 0; j < (totalStructures - 1); j++){ */
  /*  equations[(totalStructures - 1)][j + 1] = 0; */
  /* } */

  /**
   * Correction: we need only as many equations as variables. 1 == 2, 2 == 3, 3 == mfe. Transitivity implies 1 == mfe. All other equations are redundant.
   */
  size_t  lengthLGS = totalStructures - 1;      /*  */
  /* first dim is all combinations of equations and the second dim are the values x',y',... */
  double  a[lengthLGS * (totalStructures - 1)]; /*  a is the right hand side of the lgs */
  double  b[lengthLGS];                         /*  b is the left hand side of the lgs, i.e. the energie difference energy(s) - energy(s') */
  int     lgsIndex  = 0;
  int     i         = 0;
  equations[i][0] = energies[i];
  for (int j = i + 1; j < totalStructures; j++, i++) {
    for (int k = 1; k < totalStructures; k++)
      a[(lgsIndex) + (lengthLGS) * (k - 1)] = (double)(equations[i][k] - equations[j][k]);
    b[lgsIndex] = (double)-(equations[i][0] - equations[j][0]);
    lgsIndex++;
  }

  int     info;
  int     m     = lengthLGS;
  int     n     = totalStructures - 1;
  int     nrhs  = 1;
  int     lda   = m;
  int     ldb   = m;
  int     rank;

  /* Negative rcond means using default (machine precision) value */
  double  rcond = -1.0;
  double  wkopt;
  double  *work = NULL;
  /* Local arrays */
  /* iwork dimension should be at least 3*min(m,n)*nlvl + 11*min(m,n),
   * where nlvl = max( 0, int( log_2( min(m,n)/(smlsiz+1) ) )+1 )
   * and smlsiz = 25 */
  int     nlvl = max(0, (int)(log((float)(min(m, n)) / 2.)) + 1);
  int     iwork[3 * min(m, n) * nlvl + 11 * min(m, n)];
  double  s[m];

  /* Executable statements */
  /* printf( " DGELSD Example Program Results\n" ); */
  /* Query and allocate the optimal workspace */
  int     lwork = -1;

  dgelsd_(&m, &n, &nrhs, a, &lda, b, &ldb, s, &rcond, &rank, &wkopt, &lwork, iwork, &info); /* estimate workspace (=> lwork = -1). */
  lwork = (int)wkopt;
  work  = (double *)vrna_alloc(lwork * sizeof(double));
  /* Solve the equations A*X = B */
  dgelsd_(&m, &n, &nrhs, a, &lda, b, &ldb, s, &rcond, &rank, work, &lwork, iwork, &info);

  /*  printf("The linear system has rank %d;\n", rank); */
  /* Check for convergence */
  if (info > 0) {
    printf("The algorithm computing SVD failed to converge;\n");
    printf("the least squares solution could not be computed.\n");
    exit(1);
  }

  /* Free workspace */
  free((void *)work);

  /* output distortion for unique structures and 0 for redundant structures. */
  printf("distortions: ");
  double *distortions = vrna_alloc(sizeof(double) * numberOfStructures);
  uniqueInd = 0;
  for (int i = 0; i < numberOfStructures; i++) {
    if (acceptedIndices[i] != -1) {
      distortions[i] = b[uniqueInd];
      uniqueInd++;
    } else {
      distortions[i] = 0;
    }

    if (i == numberOfStructures - 1)
      printf("d_x%d = %1.10f ", i, distortions[i]);
    else
      printf("d_x%d = %1.10f, ", i, distortions[i]);
  }
  printf("\n");

  return distortions;
}


double *
rxp_computeDistortions(vrna_fold_compound_t *fc,
                       const char           **structures,
                       size_t               numberOfStructures,
                       float                mfe,
                       const char           *mfeStructure)
{
  int acceptedIndices[numberOfStructures];
  int numberAccepted = 0;

  /* filter for unique input structures, that are not the mfe. (accepts only the last unique structure) */
  for (int i = 0; i < numberOfStructures; i++) {
    const char  *s1             = structures[i];
    int         acceptStructure = 1;
    if (strcmp(s1, mfeStructure) == 0) {
      acceptedIndices[i] = -1;
      continue;
    } else {
      for (int j = i + 1; j < numberOfStructures; j++) {
        const char *s2 = structures[j];
        if (strcmp(s1, s2) == 0) {
          acceptStructure = 0;
          break;
        }
      }
    }

    if (acceptStructure == 1) {
      acceptedIndices[i] = 1;
      numberAccepted++;
    } else {
      acceptedIndices[i] = -1;
    }
  }

  if (numberAccepted < 1) {
    /* no distortion possible */
    double *distortions = vrna_alloc(sizeof(double) * numberOfStructures);
    for (int i = 0; i < numberOfStructures; i++)
      distortions[i] = 0;
    return distortions;
  }

  const char  *uniqueStructures[numberAccepted + 1]; /* +1 for the mfe */
  int         uniqueInd = 0;
  for (int i = 0; i < numberOfStructures; i++) {
    if (acceptedIndices[i] != -1) {
      uniqueStructures[uniqueInd] = structures[i];
      uniqueInd++;
    }
  }
  uniqueStructures[numberAccepted] = mfeStructure;

  float energies[numberAccepted + 1]; /* +1 for the mfe */
  energies[numberAccepted] = mfe;

  for (int i = 0; i < numberAccepted; i++) {
    const char *structure = uniqueStructures[i];
    energies[i] = vrna_eval_structure(fc, structure);
  }

  /* structures + mfeStr */
  /* energies + mfe */
  size_t totalStructures = numberAccepted + 1;
  /**
   * construct equations, which look like:
   * S1: E(s1) + x'*d(s1, s1) + y'*d(s1, s2) + z'*d(s1, s3) = c
   * S2: E(s2) + x'*d(s2, s1) + y'*d(s2, s2) + z'*d(s2, s3) = c
   * ...
   * S4: E(mfeStructure) + ... = c
   */

  /* first dim number of equations and the second dim are the values energy(s),x',y',... */
  float   equations[totalStructures][totalStructures];
  size_t  distance;
  for (int i = 0; i < totalStructures; i++) {
    equations[i][0] = energies[i];
    for (int j = 0; j < (totalStructures - 1); j++) {
      const char  *str1 = uniqueStructures[i];
      const char  *str2 = uniqueStructures[j];
      distance            = bp_distance(str1, str2);
      equations[i][j + 1] = distance;
    }
  }

  /**
   * now we construct the lgs, which are all combinations of equations
   * for 2 references + mfe: 1 - 2, 1 - 3, 2 - 3
   */

  /*
   * size_t lengthLGS = ((totalStructures * (totalStructures - 1)) / 2);
   * //first dim is all combinations of equations and the second dim are the values x',y',...
   * double a[lengthLGS * (totalStructures - 1)]; // a is the right hand side of the lgs
   * double b[lengthLGS]; // b is the left hand side of the lgs, i.e. the energie difference energy(s) - energy(s')
   * int lgsIndex = 0;
   * for(int i = 0; i < totalStructures; i++){
   * equations[i][0] = energies[i];
   * for(int j = i + 1; j < totalStructures; j++){
   * for(int k = 1; k < totalStructures; k++){
   * a[(lgsIndex) + (lengthLGS) * (k - 1)] = (double) (equations[i][k] - equations[j][k]);
   * }
   * b[lgsIndex] = (double) (equations[i][0] - equations[j][0]);
   * lgsIndex++;
   * }
   * }
   */

  /**
   * Correction: we need only as many equations as variables. 1 == 2, 2 == 3, 3 == mfe. Transitivity implies 1 == mfe. All other equations are redundant.
   */
  size_t  lengthLGS = totalStructures - 1;
  /* first dim is all combinations of equations and the second dim are the values x',y',... */
  double  a[lengthLGS * (totalStructures - 1)]; /*  a is the left hand side of the lgs */
  double  b[lengthLGS];                         /*  b is the right hand side of the lgs, i.e. the energie difference energy(s) - energy(s') */
  int     lgsIndex  = 0;
  int     i         = 0;
  equations[i][0] = energies[i];
  for (int j = i + 1; j < totalStructures; j++, i++) {
    for (int k = 1; k < totalStructures; k++)
      a[(lgsIndex) + (lengthLGS) * (k - 1)] = (double)(equations[i][k] - equations[j][k]);
    b[lgsIndex] = (double)-(equations[i][0] - equations[j][0]);  /* *(-1) because b is right hand side. */
    lgsIndex++;
  }

  int info;
  int m     = lengthLGS;
  int n     = totalStructures - 1;
  int nrhs  = 1;
  int lda   = m;
  int ldb   = m;
  int rank;

  /*   print_matrix("a", lda, n, a, lda); */
  /*   print_matrix("b", ldb, 1, b, ldb); */

  /* Negative rcond means using default (machine precision) value */
  double  rcond = -1.0;
  double  wkopt;
  double  *work = NULL;
  /* Local arrays */
  /* iwork dimension should be at least 3*min(m,n)*nlvl + 11*min(m,n),
   * where nlvl = max( 0, int( log_2( min(m,n)/(smlsiz+1) ) )+1 )
   * and smlsiz = 25 */
  /* int iwork[3*m*0+11*m]; */
  int     nlvl = max(0, (int)(log((float)(min(m, n)) / 2.)) + 1);
  int     iwork[3 * min(m, n) * nlvl + 11 * min(m, n)];
  double  s[m];

  /* Executable statements */
  /* printf( " DGELSD Example Program Results\n" ); */
  /* Query and allocate the optimal workspace */
  int     lwork = -1;

  dgelsd_(&m, &n, &nrhs, a, &lda, b, &ldb, s, &rcond, &rank, &wkopt, &lwork, iwork, &info); /* estimate workspace (=> lwork = -1). */

  lwork = (int)wkopt;
  work  = (double *)vrna_alloc(lwork * sizeof(double));
  /* Solve the equations A*X = B */
  dgelsd_(&m, &n, &nrhs, a, &lda, b, &ldb, s, &rcond, &rank, work, &lwork, iwork, &info);

  /*  printf("The linear system has rank %d;\n", rank); */
  /* Check for convergence */
  if (info > 0) {
    printf("The algorithm computing SVD failed to converge;\n");
    printf("the least squares solution could not be computed.\n");
    exit(1);
  }

  /*
   * // Print minimum norm solution
   * print_matrix("Minimum norm solution", n, nrhs, b, ldb);
   * // Print effective rank
   * printf("\n Effective rank = %6i\n", rank);
   * // Print singular values
   * print_matrix("Singular values", 1, m, s, 1);
   */

  /* print_matrix("a_after", lda, n, a, lda); */
  /*   print_matrix("b_after", ldb, 1, b, ldb); */
  /* Free workspace */
  free((void *)work);

  /* output distortion for unique structures and 0 for redundant structures. */
  printf("distortions: ");
  double *distortions = vrna_alloc(sizeof(double) * numberOfStructures);
  uniqueInd = 0;
  for (int i = 0; i < numberOfStructures; i++) {
    if (acceptedIndices[i] != -1) {
      distortions[i] = b[uniqueInd];
      uniqueInd++;
    } else {
      distortions[i] = 0;
    }

    if (i == numberOfStructures - 1)
      printf("d_x%d = %1.10f ", i, distortions[i]);
    else
      printf("d_x%d = %1.10f, ", i, distortions[i]);
  }
  printf("\n");

  return distortions;
}


double *
rxp_computeDistortionsWithMFE(vrna_fold_compound_t  *fc,
                              const char            **structures,
                              size_t                numberOfStructures)
{
  char    *mfeStructure = (char *)vrna_alloc(sizeof(char) * (fc->length + 1));
  double  mfe           = (double)vrna_mfe(fc, mfeStructure);
  double  *distortions  = rxp_computeDistortions(fc,
                                                 structures,
                                                 numberOfStructures,
                                                 mfe,
                                                 mfeStructure);

  free(mfeStructure);
  return distortions;
}


/* Auxiliary routine: printing a matrix */
void
print_matrix(char   *desc,
             int    m,
             int    n,
             double *a,
             int    lda)
{
  int i, j;

  printf("\n %s\n", desc);
  for (i = 0; i < m; i++) {
    for (j = 0; j < n; j++)
      printf(" %6.2f", a[i + j * lda]);
    printf("\n");
  }
}

#if WITH_DIST_CLASS_SC
static double
distortion_default(int i,
                   int j,
                   int k,
                   int l,
                   unsigned char decomp,
                   int *distance,
                   sc_dist_class_t *d)
{
  double result, *distortions;
  int numberOfRefs;

  distortions  = ((kl_soft_constraints_MD *)d->f_data)->distortions;
  numberOfRefs = d->ref_num;

  result = 0.;
  for (int r = 0; r < numberOfRefs; r++)
    result += distortions[r] * (double)distance[r];

  result = result * 100.;
  return result;
}

static double
distortion_repel(int i,
                 int j,
                 int k,
                 int l,
                 unsigned char decomp,
                 int *distance,
                 sc_dist_class_t *d)
{
  double  result, *distortions;
  int     numberOfRefs;
  size_t  *maxDist;

  distortions  = ((kl_soft_constraints_MD *)d->f_data)->distortions;
  maxDist      = ((kl_soft_constraints_MD *)d->f_data)->maxDistances;
  numberOfRefs = d->ref_num;

  result = 0.;
  for (int r = 0; r < numberOfRefs; r++)
    result += distortions[r] * ((double)maxDist[r] - distance[r]);

  result = result * 100.;
  return result;
}

static void
kl_datastructures_MD_destroy(void *data)
{
  kl_soft_constraints_MD *d = (kl_soft_constraints_MD *)data;

  /* free(d->distortions); */
  free(d->maxDistances);
  free(d);
}


static sc_dist_class_t *
kl_init_datastructures_MD(vrna_fold_compound_t  *vc,
                          const char            **referenceStructures,
                          int                   numberOfReferenceStructures,
                          double                *distortions,
                          int                   repel)
{
  sc_dist_class_t         *d;
  kl_soft_constraints_MD  *data;

  char                    *s = vc->sequence;

  d = sc_dist_class_init(vc);

  for (int i = 0; i < numberOfReferenceStructures; i++)
    sc_dist_class_add_ref(d, referenceStructures[i]);

  d->f = (repel) ? &distortion_repel : &distortion_default;

  /* manage MD specific data structrue */
  data                      = (kl_soft_constraints_MD *)vrna_alloc(sizeof(kl_soft_constraints_MD));
  data->distortions         = distortions;
  data->maxDistances        = vrna_alloc(sizeof(size_t *) * numberOfReferenceStructures);
  for (int i = 0; i < numberOfReferenceStructures; i++)
    data->maxDistances[i] = getMaximalPossibleBPdistance(s, referenceStructures[i]);

  d->f_data = (void *)data;
  d->f_free = &kl_datastructures_MD_destroy;

  return d;
}

#else
static kl_soft_constraints_MD *
kl_init_datastructures_MD(vrna_fold_compound_t  *vc,
                          const char            **referenceStructures,
                          int                   numberOfReferenceStructures,
                          double                *distortions,
                          int                   repel)
{
  kl_soft_constraints_MD  *data;
  unsigned int            n;
  char                    *s = vc->sequence;

  n = strlen(s);

  /* alloc all memory */
  data                      = (kl_soft_constraints_MD *)vrna_alloc(sizeof(kl_soft_constraints_MD));
  data->distortions         = distortions;
  data->numberOfReferences  = numberOfReferenceStructures;
  data->references          = referenceStructures;

  /* memory to free: */
  data->kT                      = vc->exp_params->kT;
  data->idx                     = vrna_idx_row_wise(n);
  data->referencesAsPairtaibles = vrna_alloc(numberOfReferenceStructures * sizeof(short *));
  data->referencesBPs           = vrna_alloc(numberOfReferenceStructures * sizeof(short *));
  for (int i = 0; i < numberOfReferenceStructures; i++) {
    data->referencesAsPairtaibles[i] = vrna_ptable(referenceStructures[i]);
    /* matrix containing number of basepairs of reference structure1 in interval [i,j] */
    data->referencesBPs[i] = vrna_refBPcnt_matrix(data->referencesAsPairtaibles[i], TURN);
  }

  data->maxDistances = vrna_alloc(sizeof(size_t *) * numberOfReferenceStructures);
  for (int i = 0; i < numberOfReferenceStructures; i++)
    data->maxDistances[i] = getMaximalPossibleBPdistance(s, data->references[i]);
  data->repel = repel;
  return data;
}


static void
free_kl_soft_constraints_MD(void *data)
{
  kl_soft_constraints_MD *dat = (kl_soft_constraints_MD *)data;

  free(dat->idx);
  for (int i = 0; i < dat->numberOfReferences; i++) {
    free(dat->referencesAsPairtaibles[i]);
    free(dat->referencesBPs[i]);
  }
  free(dat->referencesAsPairtaibles);
  free(dat->referencesBPs);
  free(dat);
}


static FLT_OR_DBL
kl_pseudo_energy_MD(int   i,
                    int   j,
                    int   k,
                    int   l,
                    unsigned char  decomp,
                    void  *data)
{
  int                     ij, kl;
  kl_soft_constraints_MD  *distData = ((kl_soft_constraints_MD *)data);
  int                     *idx      = distData->idx;

  int                     numberOfRefs    = distData->numberOfReferences;
  short                   **references_pt = distData->referencesAsPairtaibles;
  unsigned int            **referenceBPs  = distData->referencesBPs;
  double                  *distortions    = distData->distortions;

  int                     base_dx[numberOfRefs];

  for (int r = 0; r < numberOfRefs; r++)
    base_dx[r] = (references_pt[r][i] != (unsigned int)j) ? 1 : -1;

  ij  = idx[i] - j;
  kl  = idx[k] - l;

  int distances[numberOfRefs];

  switch (decomp) {
    /* cases where we actually introduce a base pair */

    case VRNA_DECOMP_PAIR_HP:
      for (int r = 0; r < numberOfRefs; r++)
        distances[r] = base_dx[r] + referenceBPs[r][ij];
      break;

    case VRNA_DECOMP_PAIR_IL:
      for (int r = 0; r < numberOfRefs; r++)
        distances[r] = base_dx[r] + referenceBPs[r][ij] - referenceBPs[r][kl];
      break;

    case VRNA_DECOMP_PAIR_ML:
      for (int r = 0; r < numberOfRefs; r++)
        distances[r] = base_dx[r] + referenceBPs[r][ij] - referenceBPs[r][kl];
      break;

    /* cases where we split a segment into one or two subsegments */

    case VRNA_DECOMP_ML_STEM:
      for (int r = 0; r < numberOfRefs; r++)
        distances[r] = referenceBPs[r][ij] - referenceBPs[r][kl];
      break;

    case VRNA_DECOMP_ML_ML:
      for (int r = 0; r < numberOfRefs; r++)
        distances[r] = referenceBPs[r][ij] - referenceBPs[r][kl];
      break;

    case VRNA_DECOMP_ML_ML_ML:
      for (int r = 0; r < numberOfRefs; r++)
        distances[r] = referenceBPs[r][ij] - referenceBPs[r][idx[i] - k] -
                       referenceBPs[r][idx[l] - j];
      break;

    case VRNA_DECOMP_ML_UP:
      for (int r = 0; r < numberOfRefs; r++)
        distances[r] = referenceBPs[r][ij];
      break;

    case VRNA_DECOMP_ML_ML_STEM:
      for (int r = 0; r < numberOfRefs; r++)
        distances[r] = referenceBPs[r][ij] - referenceBPs[r][idx[i] - k] -
                       referenceBPs[r][idx[l] - j];
      break;

    case VRNA_DECOMP_ML_COAXIAL: /* (i,j) stacks onto (k,l), lets ignore this case for now */
      for (int r = 0; r < numberOfRefs; r++)
        distances[r] = 0;
      break;

    case VRNA_DECOMP_EXT_EXT:
      for (int r = 0; r < numberOfRefs; r++)
        distances[r] = referenceBPs[r][ij] - referenceBPs[r][kl];
      break;

    case VRNA_DECOMP_EXT_UP:
      for (int r = 0; r < numberOfRefs; r++)
        distances[r] = referenceBPs[r][ij];
      break;

    case VRNA_DECOMP_EXT_EXT_EXT:
      for (int r = 0; r < numberOfRefs; r++)
        distances[r] = referenceBPs[r][ij] - referenceBPs[r][idx[i] - k] -
                       referenceBPs[r][idx[l] - j];
      break;

    case VRNA_DECOMP_EXT_STEM:
      for (int r = 0; r < numberOfRefs; r++)
        distances[r] = referenceBPs[r][ij] - referenceBPs[r][kl];
      break;

    case VRNA_DECOMP_EXT_EXT_STEM: /* fall through */
    case VRNA_DECOMP_EXT_STEM_EXT:
      for (int r = 0; r < numberOfRefs; r++)
        distances[r] = referenceBPs[r][ij] - referenceBPs[r][idx[i] - k] -
                       referenceBPs[r][idx[l] - j];
      break;

    case VRNA_DECOMP_EXT_EXT_STEM1:
      for (int r = 0; r < numberOfRefs; r++)
        distances[r] = referenceBPs[r][ij] - referenceBPs[r][idx[i] - k] -
                       referenceBPs[r][idx[l] - (j - 1)];
      break;

    case VRNA_DECOMP_EXT_STEM_OUTSIDE:
      for (int r = 0; r < numberOfRefs; r++) {
        distances[r]  = referenceBPs[r][ij];
        distances[r]  = referenceBPs[r][ij] - referenceBPs[r][kl];
        if (k > i)
          distances[r] -= referenceBPs[r][idx[i] - (k - 1)];

        if (l < j)
          distances[r] -= referenceBPs[r][idx[l + 1] - j];
      }
      break;

    default:
      for (int r = 0; r < numberOfRefs; r++)
        distances[r] = 0;
      break;
  }

  double result = 0;
  if (!distData->repel) {
    for (int r = 0; r < numberOfRefs; r++)
      result += distortions[r] * distances[r];
  } else {
    for (int r = 0; r < numberOfRefs; r++)
      result += distortions[r] * (distData->maxDistances[r] - distances[r]);
  }

  result = result * 100;

  return result;
}


static FLT_OR_DBL
kl_exp_pseudo_energy_MD(int   i,
                        int   j,
                        int   k,
                        int   l,
                        unsigned char  decomp,
                        void  *data)
{
  double  kT      = ((kl_soft_constraints_MD *)data)->kT;
  double  result  = exp((-10. * (double)kl_pseudo_energy_MD(i, j, k, l, decomp, data)) / kT);

  return result;
}

#endif

static void
fillGridStepwiseBothRef_MD(vrna_fold_compound_t *vc,
                           gridLandscapeT       *grid,
                           float                relaxFactor,
                           int                  relax,
                           int                  shift,
                           int                  shift_to_first,
                           int                  verbose,
                           int                  maxIterations,
                           int                  maxSteps)
{
  if (maxSteps <= 0) {
    fprintf(stderr, "Error: the stepsize has to be positive and greater than zero!");
    return;
  }

#if WITH_DIST_CLASS_SC
  sc_dist_class_t         *d    = (sc_dist_class_t *)vc->sc->data;
  kl_soft_constraints_MD  *data = (kl_soft_constraints_MD *)d->f_data;
  char                    *s1   = d->references[0];
  char                    *s2   = d->references[1];
  double                  tmp_x = data->distortions[0];
  double                  tmp_y = data->distortions[1];
#else
  kl_soft_constraints_MD  *data = (kl_soft_constraints_MD *)vc->sc->data;
  char                    *s1   = data->references[0];
  char                    *s2   = data->references[1];
  double                  tmp_x = data->distortions[0];
  double                  tmp_y = data->distortions[1];
#endif

  for (int j = 0; j < maxSteps; j++) {
    if (shift) {
      if (shift_to_first) {
        data->distortions[0]  += relaxFactor * data->distortions[0] / maxSteps;
        data->distortions[1]  -= relaxFactor * data->distortions[1] / maxSteps;
      } else {
        data->distortions[0]  -= relaxFactor * data->distortions[0] / maxSteps;
        data->distortions[1]  += relaxFactor * data->distortions[1] / maxSteps;
      }
    } else {
      if (relax) {
        data->distortions[0]  -= relaxFactor * data->distortions[0] / maxSteps;
        data->distortions[1]  -= relaxFactor * data->distortions[1] / maxSteps;
      } else {
        data->distortions[0]  += relaxFactor * data->distortions[0] / maxSteps;
        data->distortions[1]  += relaxFactor * data->distortions[1] / maxSteps;
      }
    }

    if (verbose)
      fprintf(stderr, "d_x = %1.10f, d_y = %1.10f\n", data->distortions[0], data->distortions[1]);

    fillGridWithSamples(vc, grid, s1, s2, maxIterations);
  }

  data->distortions[0]  = tmp_x;
  data->distortions[1]  = tmp_y;
}


void
fillGridStepwiseFirstRef_MD(vrna_fold_compound_t  *vc,
                            gridLandscapeT        *grid,
                            float                 relaxFactor,
                            int                   relax,
                            int                   verbose,
                            int                   maxIterations,
                            int                   maxSteps)
{
  if (maxSteps <= 0) {
    fprintf(stderr, "Error: the stepsize has to be positive and greater than zero!");
    return;
  }

#if WITH_DIST_CLASS_SC
  sc_dist_class_t         *d    = (sc_dist_class_t *)vc->sc->data;
  kl_soft_constraints_MD  *data = (kl_soft_constraints_MD *)d->f_data;
  char                    *s1   = d->references[0];
  char                    *s2   = d->references[1];
  double                  tmp_x = data->distortions[0];
#else
  kl_soft_constraints_MD  *data = (kl_soft_constraints_MD *)vc->sc->data;
  char                    *s1   = data->references[0];
  char                    *s2   = data->references[1];
  double                  tmp_x = data->distortions[0];
#endif

  for (int j = 0; j < maxSteps; j++) {
    if (relax)
      data->distortions[0] -= relaxFactor * data->distortions[0] / maxSteps;
    else
      data->distortions[0] += relaxFactor * data->distortions[0] / maxSteps;

    if (verbose)
      fprintf(stderr, "d_x = %1.10f, d_y = %1.10f\n", data->distortions[0], data->distortions[1]);

    fillGridWithSamples(vc, grid, s1, s2, maxIterations);
  }

  data->distortions[0] = tmp_x;
}


void
fillGridStepwiseSecondRef_MD(vrna_fold_compound_t *vc,
                             gridLandscapeT       *grid,
                             float                relaxFactor,
                             int                  relax,
                             int                  verbose,
                             int                  maxIterations,
                             int                  maxSteps)
{
  if (maxSteps <= 0) {
    fprintf(stderr, "Error: the stepsize has to be positive and greater than zero!");
    return;
  }

#if WITH_DIST_CLASS_SC
  sc_dist_class_t         *d    = (sc_dist_class_t *)vc->sc->data;
  kl_soft_constraints_MD  *data = (kl_soft_constraints_MD *)d->f_data;
  char                    *s1   = d->references[0];
  char                    *s2   = d->references[1];
  double                  tmp_y = data->distortions[1];
#else
  kl_soft_constraints_MD  *data = (kl_soft_constraints_MD *)vc->sc->data;
  char                    *s1   = data->references[0];
  char                    *s2   = data->references[1];
  double                  tmp_y = data->distortions[1];
#endif

  for (int j = 0; j < maxSteps; j++) {
    if (relax)
      data->distortions[1] -= relaxFactor * data->distortions[1] / maxSteps;
    else
      data->distortions[1] += relaxFactor * data->distortions[1] / maxSteps;

    if (verbose)
      fprintf(stderr, "d_x = %1.10f, d_y = %1.10f\n", data->distortions[0], data->distortions[1]);

    fillGridWithSamples(vc, grid, s1, s2, maxIterations);
  }

  data->distortions[1] = tmp_y;
}


gridLandscapeT *
estimate_landscapeMD(vrna_fold_compound_t *vc,
                     const char           **refStructures,
                     size_t               numberOfReferences,
                     int                  maxIterations,
                     char                 *extended_options,
                     double               *indicesAndPercentages,
                     size_t               lengthIndices)
{
  /* parse extended options string */
  int both_at_once    = 0;
  int relax           = 0;
  int plain           = 1; /* plain sampling, no distortion */
  int verbose         = 0;
  int shift           = 0;
  int shift_to_first  = 0;
  int normal          = 0;

  if (extended_options) {
    plain = 0;
    if (strchr(extended_options, 'N')) /* normal distortion (no shift) */
      normal = 1;

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
  char            *mfe_struct = (char *)vrna_alloc(sizeof(char) * (vc->length + 1));
  double          mmfe        = (double)vrna_mfe(vc, mfe_struct);

  char            *s    = vc->sequence;
  const char      *s1   = refStructures[0];
  const char      *s2   = refStructures[1];
  gridLandscapeT  *grid = initLandscape(s, s1, s2);

  if (plain) {
    /* prepare pf fold */
    vrna_exp_params_rescale(vc, &mmfe);
    /* backtracking with energies without distortion */
    fillGridWithSamples(vc, grid, s1, s2, maxIterations);
  } else {
    /* computeInitialDistortionMD(vc, s1, s2, &distortion_x, &distortion_y); */
    double *distortions = rxp_computeDistortions(vc,
                                                 refStructures,
                                                 numberOfReferences,
                                                 mmfe,
                                                 mfe_struct);

    /*  recompute distortions with user-specified weights. */
    if (lengthIndices > 0) {
      int     valueIndex  = 2;
      int     indexIndex  = 1;
      int     referenceIndex;
      double  percentage;
      for (int i = 1; i <= lengthIndices; i++, valueIndex += 2, indexIndex += 2) {
        referenceIndex              = (int)indicesAndPercentages[indexIndex] - 1;
        percentage                  = indicesAndPercentages[valueIndex];
        distortions[referenceIndex] *= percentage;
        if (verbose)
          printf("p0: %d %f\n", referenceIndex, percentage);
      }
    }

    double                  distortion_x  = distortions[0];
    double                  distortion_y  = distortions[1];

    /* prepare pf fold */
    int                     bp_dist_mfe_ref1  = vrna_bp_distance(mfe_struct, s1);
    int                     bp_dist_mfe_ref2  = vrna_bp_distance(mfe_struct, s2);

    double                  rescale = mmfe + (bp_dist_mfe_ref1 * distortion_x) +
                                      (bp_dist_mfe_ref2 * distortion_y);
    vrna_exp_params_rescale(vc, &rescale);

    /* apply distortion soft constraints */
#if WITH_DIST_CLASS_SC
    sc_dist_class_t *data = kl_init_datastructures_MD(vc,
                                                              refStructures,
                                                              numberOfReferences,
                                                              distortions,
                                                              0);
    vrna_sc_init(vc); /*  to remove old soft constraints */
    vrna_sc_add_data(vc, (void *)data, &sc_dist_class_destroy);
    vrna_sc_add_exp_f(vc, &sc_exp_f_dist_class);

#else
    kl_soft_constraints_MD  *data = kl_init_datastructures_MD(vc,
                                                              refStructures,
                                                              numberOfReferences,
                                                              distortions,
                                                              0);
    vrna_sc_init(vc); /*  to remove old soft constraints */
    vrna_sc_add_data(vc, (void *)data, &free_kl_soft_constraints_MD);
    vrna_sc_add_exp_f(vc, &kl_exp_pseudo_energy_MD);

#endif
    if (normal) {
      /* with distortion, but no shift. */
      fillGridWithSamples(vc, grid, s1, s2, maxIterations);
    } else {
      float relaxFactor = 1.5; /* if set to 1. final relaxation sets x and/or y to 0. */
      int   bp_dist     = vrna_bp_distance(s1, s2);
      if (both_at_once) {
        /* change potential of both references at the same time */
        fillGridStepwiseBothRef_MD(vc,
                                   grid,
                                   relaxFactor,
                                   relax,
                                   shift,
                                   shift_to_first,
                                   verbose,
                                   maxIterations,
                                   bp_dist);
      } else {
        /* change potential of first reference */
        fillGridStepwiseFirstRef_MD(vc, grid, relaxFactor, relax, verbose, maxIterations, bp_dist);
        /* change potential of second reference */
        fillGridStepwiseSecondRef_MD(vc, grid, relaxFactor, relax, verbose, maxIterations, bp_dist);
      }
    }

    free(distortions);
  }

  free(mfe_struct);

  /*  evaluate structures, set properties. */
  /*  compute undistorted energies. */
  vrna_sc_free(vc->sc);
  vc->sc                              = NULL;
  vc->params->model_details.betaScale = 1;
  vrna_exp_params_rescale(vc, &mmfe);
  float mfe     = FLT_MAX;
  float tmpMFE  = 0;
  for (int i = 0; i < grid->size1; i++) {
    for (int j = 0; j < grid->size2; j++) {
      if (grid->landscape[i][j].num_structs > 0) {
        int k;
        for (k = 0; k < grid->landscape[i][j].num_structs; k++) {
          tmpMFE = vrna_eval_structure(vc, grid->landscape[i][j].structures[k]);
          if (tmpMFE < mfe)
            grid->landscape[i][j].mfe = mfe;
        }
        grid->landscape[i][j].k = i;
        grid->landscape[i][j].l = j;
      }
    }
  }

  return grid;
}


void
addSoftconstraintsMD(vrna_fold_compound_t *vc,
                     const char           **structures,
                     int                  numberOfReferences,
                     double               *distortions,
                     int                  repel)
{
  /* apply distortion soft constraints */
#if WITH_DIST_CLASS_SC
    sc_dist_class_t *data = kl_init_datastructures_MD(vc,
                                                              structures,
                                                              numberOfReferences,
                                                              distortions,
                                                              repel);
    vrna_sc_init(vc); /*  to remove old soft constraints */
    vrna_sc_add_data(vc, (void *)data, &sc_dist_class_destroy);
    vrna_sc_add_exp_f(vc, &sc_exp_f_dist_class);

#else
    kl_soft_constraints_MD  *data = kl_init_datastructures_MD(vc,
                                                              structures,
                                                              numberOfReferences,
                                                              distortions,
                                                              repel);

    vrna_sc_init(vc); /*  to remove old soft constraints */
    vrna_sc_add_data(vc, (void *)data, &free_kl_soft_constraints_MD);
    vrna_sc_add_exp_f(vc, &kl_exp_pseudo_energy_MD);
#endif
}
