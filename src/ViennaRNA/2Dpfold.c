/*
 *    minimum free energy
 *    RNA secondary structure with
 *    basepair distance d_1 to reference structure 1 and distance d_2 to reference structure 2
 *
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <string.h>
#include <float.h>    /* #defines FLT_MAX ... */
#include "ViennaRNA/utils.h"
#include "ViennaRNA/fold_vars.h"
#include "ViennaRNA/params.h"
#include "ViennaRNA/energy_par.h"
#include "ViennaRNA/loop_energies.h"
#include "ViennaRNA/2Dpfold.h"

/*
 #################################
 # GLOBAL VARIABLES              #
 #################################
 */

/*
 #################################
 # PRIVATE VARIABLES             #
 #################################
 */

/*
 #################################
 # PRIVATE FUNCTION DECLARATIONS #
 #################################
 */
PRIVATE void  crosslink(TwoDpfold_vars *vars);


PRIVATE void  pf2D_linear(vrna_fold_compound_t *vc);


PRIVATE void  pf2D_circ(vrna_fold_compound_t *vc);


PRIVATE char *pbacktrack_circ(vrna_fold_compound_t  *vc,
                              int                   d1,
                              int                   d2);


PRIVATE void  backtrack(vrna_fold_compound_t  *vc,
                        char                  *pstruc,
                        int                   d1,
                        int                   d2,
                        unsigned int          i,
                        unsigned int          j);


PRIVATE void  backtrack_qm(vrna_fold_compound_t *vc,
                           char                 *pstruc,
                           int                  d1,
                           int                  d2,
                           unsigned int         i,
                           unsigned int         j);


PRIVATE void  backtrack_qm1(vrna_fold_compound_t  *vc,
                            char                  *pstruc,
                            int                   d1,
                            int                   d2,
                            unsigned int          i,
                            unsigned int          j);


PRIVATE void  backtrack_qm2(vrna_fold_compound_t  *vc,
                            char                  *pstruc,
                            int                   d1,
                            int                   d2,
                            unsigned int          k);


PRIVATE void  backtrack_qcH(vrna_fold_compound_t  *vc,
                            char                  *pstruc,
                            int                   d1,
                            int                   d2);


PRIVATE void  backtrack_qcI(vrna_fold_compound_t  *vc,
                            char                  *pstruc,
                            int                   d1,
                            int                   d2);


PRIVATE void  backtrack_qcM(vrna_fold_compound_t  *vc,
                            char                  *pstruc,
                            int                   d1,
                            int                   d2);


PRIVATE void  adjustArrayBoundaries(FLT_OR_DBL  ***array,
                                    int         *k_min,
                                    int         *k_max,
                                    int         **l_min,
                                    int         **l_max,
                                    int         k_min_real,
                                    int         k_max_real,
                                    int         *l_min_real,
                                    int         *l_max_real);


INLINE PRIVATE void  preparePosteriorBoundaries(int size,
                                                int shift,
                                                int *min_k,
                                                int *max_k,
                                                int **min_l,
                                                int **max_l);


INLINE PRIVATE void  updatePosteriorBoundaries(int  d1,
                                               int  d2,
                                               int  *min_k,
                                               int  *max_k,
                                               int  **min_l,
                                               int  **max_l);


INLINE PRIVATE void  prepareBoundaries(int  min_k_pre,
                                       int  max_k_pre,
                                       int  min_l_pre,
                                       int  max_l_pre,
                                       int  bpdist,
                                       int  *min_k,
                                       int  *max_k,
                                       int  **min_l,
                                       int  **max_l);


INLINE PRIVATE void  prepareArray(FLT_OR_DBL  ***array,
                                  int         min_k,
                                  int         max_k,
                                  int         *min_l,
                                  int         *max_l);


/*
 #################################
 # BEGIN OF FUNCTION DEFINITIONS #
 #################################
 */
PUBLIC vrna_sol_TwoD_pf_t *
vrna_pf_TwoD(vrna_fold_compound_t *vc,
             int                  distance1,
             int                  distance2)
{
  unsigned int        maxD1 = 0, maxD2 = 0, counter = 0;
  int                 cnt1, cnt2, k_min, k_max, l_min, l_max, ndx;
  FLT_OR_DBL          q = 0.;

  vrna_sol_TwoD_pf_t  *output;
  vrna_md_t           *md;
  vrna_mx_pf_t        *matrices;

  maxD1     = vc->maxD1;
  maxD2     = vc->maxD2;
  matrices  = vc->exp_matrices;
  md        = &(vc->exp_params->model_details);

  if (distance1 >= 0) {
    if ((unsigned int)distance1 > maxD1)
      vrna_message_warning("vrna_pf_TwoD@2Dpfold.c: limiting maximum basepair distance 1 to %u\n",
                           maxD1);
    else
      maxD1 = (unsigned int)distance1;
  }

  if (distance2 >= 0) {
    if ((unsigned int)distance2 > maxD2)
      vrna_message_warning("vrna_pf_TwoD@2Dpfold.c: limiting maximum basepair distance 2 to %u\n",
                           maxD2);
    else
      maxD2 = (unsigned int)distance2;
  }

  vc->maxD1 = maxD1;
  vc->maxD2 = maxD2;

  output = (vrna_sol_TwoD_pf_t *)vrna_alloc((((maxD1 + 1) * (maxD2 + 2)) / 2 + 2) * sizeof(vrna_sol_TwoD_pf_t));

  pf2D_linear(vc);
  if (md->circ)
    pf2D_circ(vc);

  ndx   = vc->iindx[1] - vc->length;
  k_min = (md->circ) ? matrices->k_min_Q_c : matrices->k_min_Q[ndx];
  k_max = (md->circ) ? matrices->k_max_Q_c : matrices->k_max_Q[ndx];

  for (cnt1 = k_min;
       cnt1 <= k_max;
       cnt1++) {
    l_min = (md->circ) ? matrices->l_min_Q_c[cnt1] : matrices->l_min_Q[ndx][cnt1];
    l_max = (md->circ) ? matrices->l_max_Q_c[cnt1] : matrices->l_max_Q[ndx][cnt1];
    for (cnt2 = l_min;
         cnt2 <= l_max;
         cnt2 += 2) {
      q = (md->circ) ? matrices->Q_c[cnt1][cnt2 / 2] : matrices->Q[ndx][cnt1][cnt2 / 2];
      if (q == 0.)
        continue;

      output[counter].k = cnt1;
      output[counter].l = cnt2;
      output[counter].q = q;
      counter++;
    }
  }

  /* store entry for remaining partition if it exists */
  q = (md->circ) ? matrices->Q_c_rem : matrices->Q_rem[ndx];
  if (q != 0.) {
    output[counter].k = -1;
    output[counter].l = -1;
    output[counter].q = q;
    counter++;
  }

  /* insert end-marker entry */
  output[counter].k = output[counter].l = INF;
  counter++;

  /* resize to actual dataset amount */
  output = (vrna_sol_TwoD_pf_t *)vrna_realloc(output, sizeof(vrna_sol_TwoD_pf_t) * counter);
  return output;
}


#if 0
PUBLIC FLT_OR_DBL **
TwoDpfold(TwoDpfold_vars  *vars,
          int             distance1,
          int             distance2)
{
  unsigned int  i;
  unsigned int  maxD1 = 0;
  unsigned int  maxD2 = 0;
  unsigned int  mm;
  int           cnt1, cnt2;

  FLT_OR_DBL    **output;

  initialize_TwoDpfold_vars(vars);

  vars->S   = encode_sequence(vars->sequence, 0);
  vars->S1  = encode_sequence(vars->sequence, 1);
  make_ptypes2(vars);

  for (i = 1; i <= (unsigned int)vars->reference_pt1[0]; i++)
    if (i < (unsigned int)vars->reference_pt1[i])
      maxD1++;

  for (i = 1; i <= (unsigned int)vars->reference_pt2[0]; i++)
    if (i < (unsigned int)vars->reference_pt2[i])
      maxD2++;

  mm    = maximumMatching(vars->sequence);
  maxD1 += mm;
  maxD2 += mm;

  if (distance1 >= 0) {
    if ((unsigned int)distance1 > maxD1)
      fprintf(stderr, "limiting maximum basepair distance 1 to %u\n", maxD1);

    maxD1 = (unsigned int)distance1;
  }

  if (distance2 >= 0) {
    if ((unsigned int)distance2 > maxD2)
      fprintf(stderr, "limiting maximum basepair distance 2 to %u\n", maxD2);

    maxD2 = (unsigned int)distance2;
  }

  vars->maxD1 = maxD1;
  vars->maxD2 = maxD2;


  output = (FLT_OR_DBL **)vrna_alloc(sizeof(FLT_OR_DBL *) * (maxD1 + 1));
  pf2D_linear(vars);
  int ndx = vars->my_iindx[1] - vars->seq_length;
  for (cnt1 = vars->k_min_values[ndx]; cnt1 <= MIN2(vars->k_max_values[ndx], vars->maxD1); cnt1++) {
    output[cnt1] = (FLT_OR_DBL *)vrna_alloc((vars->maxD2 + 1) * sizeof(FLT_OR_DBL));
    for (cnt2 = vars->l_min_values[ndx][cnt1]; cnt2 <= MIN2(vars->l_max_values[ndx][cnt1], vars->maxD2); cnt2 += 2)
      output[cnt1][cnt2] = vars->Q[ndx][cnt1][cnt2 / 2];
  }
  return output;
}


PUBLIC FLT_OR_DBL **
TwoDpfold_circ(TwoDpfold_vars *vars,
               int            distance1,
               int            distance2)
{
  unsigned int  i;
  unsigned int  maxD1 = 0;
  unsigned int  maxD2 = 0;
  unsigned int  mm;
  int           cnt1, cnt2;
  FLT_OR_DBL    **output;

  initialize_TwoDpfold_vars(vars);

  vars->S   = encode_sequence(vars->sequence, 0);
  vars->S1  = encode_sequence(vars->sequence, 1);
  make_ptypes2(vars);

  for (i = 1; i <= (unsigned int)vars->reference_pt1[0]; i++)
    if (i < (unsigned int)vars->reference_pt1[i])
      maxD1++;

  for (i = 1; i <= (unsigned int)vars->reference_pt2[0]; i++)
    if (i < (unsigned int)vars->reference_pt2[i])
      maxD2++;

  mm    = maximumMatching(vars->sequence);
  maxD1 += mm;
  maxD2 += mm;

  if (distance1 >= 0) {
    if ((unsigned int)distance1 > maxD1)
      fprintf(stderr, "limiting maximum basepair distance 1 to %u\n", maxD1);

    maxD1 = (unsigned int)distance1;
  }

  if (distance2 >= 0) {
    if ((unsigned int)distance2 > maxD2)
      fprintf(stderr, "limiting maximum basepair distance 2 to %u\n", maxD2);

    maxD2 = (unsigned int)distance2;
  }

  vars->maxD1 = maxD1;
  vars->maxD2 = maxD2;

  output = (FLT_OR_DBL **)vrna_alloc(sizeof(FLT_OR_DBL *) * (maxD1 + 1));
  pf2D_linear(vars);
  pf2D_circ(vars);

  for (cnt1 = vars->k_min_values_qc; cnt1 <= MIN2(vars->k_max_values_qc, vars->maxD1); cnt1++) {
    output[cnt1] = (FLT_OR_DBL *)vrna_alloc((vars->maxD2 + 1) * sizeof(FLT_OR_DBL));
    for (cnt2 = vars->l_min_values_qc[cnt1]; cnt2 <= MIN2(vars->l_max_values_qc[cnt1], vars->maxD2); cnt2 += 2)
      output[cnt1][cnt2] = vars->Q_c[cnt1][cnt2 / 2];
  }
  return output;
}


#endif

PRIVATE void
pf2D_linear(vrna_fold_compound_t *vc)
{
  char              *sequence, *ptype;
  short             *S1, *reference_pt1, *reference_pt2;
  unsigned int      *referenceBPs1, *referenceBPs2,
                    d, i, j, ij, seq_length, maxD1,
                    maxD2, *mm1, *mm2, *bpdist;
  int               *my_iindx, *jindx, circ, cnt1, cnt2, cnt3, cnt4, *rtype;
  double            max_real;
  FLT_OR_DBL        *scale, Qmax;
  vrna_exp_param_t  *pf_params;
  vrna_mx_pf_t      *matrices;
  vrna_md_t         *md;

  max_real = (sizeof(FLT_OR_DBL) == sizeof(float)) ? FLT_MAX : DBL_MAX;

  pf_params     = vc->exp_params;
  md            = &(pf_params->model_details);
  matrices      = vc->exp_matrices;
  sequence      = vc->sequence;
  seq_length    = vc->length;
  maxD1         = vc->maxD1;
  maxD2         = vc->maxD2;
  S1            = vc->sequence_encoding;
  ptype         = vc->ptype;
  rtype         = &(md->rtype[0]);
  scale         = matrices->scale;
  reference_pt1 = vc->reference_pt1;
  reference_pt2 = vc->reference_pt2;
  my_iindx      = vc->iindx;
  jindx         = vc->jindx;
  referenceBPs1 = vc->referenceBPs1;
  referenceBPs2 = vc->referenceBPs2;
  dangles       = md->dangles;
  circ          = md->circ;
  mm1           = vc->mm1;
  mm2           = vc->mm2;
  bpdist        = vc->bpdist;
  Qmax          = 0.;

  /*array initialization ; qb,qm,q
   * qb,qm,q (i,j) are stored as ((n+1-i)*(n-i) div 2 + n+1-j */

  for (j = 1; j <= seq_length; j++)
    for (i = (j > TURN ? (j - TURN) : 1); i <= j; i++) {
      ij                        = my_iindx[i] - j;
      matrices->k_min_Q[ij]     = 0;
      matrices->k_max_Q[ij]     = 0;
      matrices->l_min_Q[ij]     = (int *)vrna_alloc(sizeof(int));
      matrices->l_max_Q[ij]     = (int *)vrna_alloc(sizeof(int));
      matrices->l_min_Q[ij][0]  = 0;
      matrices->l_max_Q[ij][0]  = 0;
      matrices->Q[ij]           = (FLT_OR_DBL **)vrna_alloc(sizeof(FLT_OR_DBL *));
      matrices->Q[ij][0]        = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL));
      matrices->Q[ij][0][0]     = 1.0 * scale[j - i + 1];
    }


  for (d = TURN + 2; d <= seq_length; d++) {
    /* i,j in [1..seq_length] */
#ifdef _OPENMP
#pragma omp parallel for private(i, j, ij, cnt1, cnt2, cnt3, cnt4)
#endif
    for (j = d; j <= seq_length; j++) {
      unsigned int  k, l, kl, u, ii, dij;
      int           no_close, type, type_2, tt, da, db, base_da, base_db;
      FLT_OR_DBL    temp2, aux_en;

      i     = j - d + 1;
      ij    = my_iindx[i] - j;
      dij   = j - i - 1;
      type  = ptype[jindx[j] + i];


      no_close = (((type == 3) || (type == 4)) && no_closingGU);

      if (type) {
        /* we have a pair */

        int k_min_Q_B, k_max_Q_B, l_min_Q_B, l_max_Q_B;
        int k_min_post_b, k_max_post_b, *l_min_post_b, *l_max_post_b;
        int update_b = 0;

        if (!matrices->Q_B[ij]) {
          update_b  = 1;
          k_min_Q_B = l_min_Q_B = 0;
          k_max_Q_B = mm1[ij] + referenceBPs1[ij];
          l_max_Q_B = mm2[ij] + referenceBPs2[ij];

          prepareBoundaries(k_min_Q_B,
                            k_max_Q_B,
                            l_min_Q_B,
                            l_max_Q_B,
                            bpdist[ij],
                            &matrices->k_min_Q_B[ij],
                            &matrices->k_max_Q_B[ij],
                            &matrices->l_min_Q_B[ij],
                            &matrices->l_max_Q_B[ij]
                            );
          preparePosteriorBoundaries(matrices->k_max_Q_B[ij] - matrices->k_min_Q_B[ij] + 1,
                                     matrices->k_min_Q_B[ij],
                                     &k_min_post_b,
                                     &k_max_post_b,
                                     &l_min_post_b,
                                     &l_max_post_b
                                     );

          prepareArray(&matrices->Q_B[ij],
                       matrices->k_min_Q_B[ij],
                       matrices->k_max_Q_B[ij],
                       matrices->l_min_Q_B[ij],
                       matrices->l_max_Q_B[ij]
                       );
        }

        /* hairpin ----------------------------------------------*/

        /* get distance to reference if closing the hairpin
         *  d1a = dbp(T1_{i,j}, {i,j})
         */
        base_da = ((unsigned int)reference_pt1[i] != j) ? 1 : -1;
        base_db = ((unsigned int)reference_pt2[i] != j) ? 1 : -1;

        da  = base_da + referenceBPs1[ij];
        db  = base_db + referenceBPs2[ij];

        if (!no_close) {
          if ((da >= 0) && (db >= 0)) {
            if (((unsigned int)da <= maxD1) && ((unsigned int)db <= maxD2)) {
              matrices->Q_B[ij][da][db / 2] = exp_E_Hairpin(dij, type, S1[i + 1], S1[j - 1], sequence + i - 1, pf_params) * scale[dij + 2];
              if (update_b) {
                updatePosteriorBoundaries(da,
                                          db,
                                          &k_min_post_b,
                                          &k_max_post_b,
                                          &l_min_post_b,
                                          &l_max_post_b
                                          );
              }
            } else {
              matrices->Q_B_rem[ij] = exp_E_Hairpin(dij, type, S1[i + 1], S1[j - 1], sequence + i - 1, pf_params) * scale[dij + 2];
            }
          }
        }

        /*--------------------------------------------------------
         *  check for elementary structures involving more than one
         *  closing pair.
         *  --------------------------------------------------------*/
        for (k = i + 1; k <= MIN2(j - 2 - TURN, i + MAXLOOP + 1); k++) {
          unsigned int minl, ln_pre;
          minl    = k + TURN + 1;
          ln_pre  = dij + k;
          if (ln_pre > minl + MAXLOOP)
            minl = ln_pre - MAXLOOP - 1;

          for (l = minl; l < j; l++) {
            kl      = my_iindx[k] - l;
            type_2  = ptype[jindx[l] + k];

            if (type_2 == 0)
              continue;

            type_2  = rtype[type_2];
            aux_en  = exp_E_IntLoop(k - i - 1, j - l - 1, type, type_2, S1[i + 1], S1[j - 1], S1[k - 1], S1[l + 1], pf_params) * scale[k - i + j - l];

            /* get distance to reference if closing the interior loop
             *  d2 = dbp(S_{i,j}, S_{k,l} + {i,j})
             */
            da  = base_da + referenceBPs1[ij] - referenceBPs1[kl];
            db  = base_db + referenceBPs2[ij] - referenceBPs2[kl];

            if (matrices->Q_B_rem[kl])
              matrices->Q_B_rem[ij] += matrices->Q_B_rem[kl] * aux_en;

            if (!matrices->Q_B[kl])
              continue;

            for (cnt1 = matrices->k_min_Q_B[kl];
                 cnt1 <= matrices->k_max_Q_B[kl];
                 cnt1++)
              for (cnt2 = matrices->l_min_Q_B[kl][cnt1];
                   cnt2 <= matrices->l_max_Q_B[kl][cnt1];
                   cnt2 += 2) {
                if (((cnt1 + da) <= maxD1) && ((cnt2 + db) <= maxD2)) {
                  matrices->Q_B[ij][cnt1 + da][(cnt2 + db) / 2] += matrices->Q_B[kl][cnt1][cnt2 / 2] * aux_en;
                  if (update_b) {
                    updatePosteriorBoundaries(da + cnt1,
                                              db + cnt2,
                                              &k_min_post_b,
                                              &k_max_post_b,
                                              &l_min_post_b,
                                              &l_max_post_b
                                              );
                  }
                } else {
                  matrices->Q_B_rem[ij] += matrices->Q_B[kl][cnt1][cnt2 / 2] * aux_en;
                }
              }
          } /* end l-loop */
        }   /* end k-loop */

        /* multi-loop contribution ------------------------*/
        if (!no_close) {
          for (u = i + TURN + 2; u < j - TURN - 2; u++) {
            tt    = rtype[type];
            temp2 = pf_params->expMLclosing * exp_E_MLstem(tt, S1[j - 1], S1[i + 1], pf_params) * scale[2];

            if (matrices->Q_M_rem[my_iindx[i + 1] - u]) {
              if (matrices->Q_M1[jindx[j - 1] + u + 1]) {
                for (cnt1 = matrices->k_min_Q_M1[jindx[j - 1] + u + 1];
                     cnt1 <= matrices->k_max_Q_M1[jindx[j - 1] + u + 1];
                     cnt1++)
                  for (cnt2 = matrices->l_min_Q_M1[jindx[j - 1] + u + 1][cnt1];
                       cnt2 <= matrices->l_max_Q_M1[jindx[j - 1] + u + 1][cnt1];
                       cnt2 += 2)
                    matrices->Q_B_rem[ij] += matrices->Q_M_rem[my_iindx[i + 1] - u] * matrices->Q_M1[jindx[j - 1] + u + 1][cnt1][cnt2 / 2] * temp2;
              }

              if (matrices->Q_M1_rem[jindx[j - 1] + u + 1])
                matrices->Q_B_rem[ij] += matrices->Q_M_rem[my_iindx[i + 1] - u] * matrices->Q_M1_rem[jindx[j - 1] + u + 1] * temp2;
            }

            if (matrices->Q_M1_rem[jindx[j - 1] + u + 1]) {
              if (matrices->Q_M[my_iindx[i + 1] - u]) {
                for (cnt1 = matrices->k_min_Q_M[my_iindx[i + 1] - u];
                     cnt1 <= matrices->k_max_Q_M[my_iindx[i + 1] - u];
                     cnt1++)
                  for (cnt2 = matrices->l_min_Q_M[my_iindx[i + 1] - u][cnt1];
                       cnt2 <= matrices->l_max_Q_M[my_iindx[i + 1] - u][cnt1];
                       cnt2 += 2)
                    matrices->Q_B_rem[ij] += matrices->Q_M[my_iindx[i + 1] - u][cnt1][cnt2 / 2] * matrices->Q_M1_rem[jindx[j - 1] + u + 1] * temp2;
              }
            }

            /* get distance to reference if closing the multiloop
             *  dist3 = dbp(S_{i,j}, {i,j} + S_{i+1,u} + S_{u+1,j-1})
             */
            da  = base_da + referenceBPs1[ij] - referenceBPs1[my_iindx[i + 1] - u] - referenceBPs1[my_iindx[u + 1] - j + 1];
            db  = base_db + referenceBPs2[ij] - referenceBPs2[my_iindx[i + 1] - u] - referenceBPs2[my_iindx[u + 1] - j + 1];

            if (!matrices->Q_M[my_iindx[i + 1] - u])
              continue;

            if (!matrices->Q_M1[jindx[j - 1] + u + 1])
              continue;

            for (cnt1 = matrices->k_min_Q_M[my_iindx[i + 1] - u];
                 cnt1 <= matrices->k_max_Q_M[my_iindx[i + 1] - u];
                 cnt1++)
              for (cnt2 = matrices->l_min_Q_M[my_iindx[i + 1] - u][cnt1];
                   cnt2 <= matrices->l_max_Q_M[my_iindx[i + 1] - u][cnt1];
                   cnt2 += 2) {
                for (cnt3 = matrices->k_min_Q_M1[jindx[j - 1] + u + 1];
                     cnt3 <= matrices->k_max_Q_M1[jindx[j - 1] + u + 1];
                     cnt3++)
                  for (cnt4 = matrices->l_min_Q_M1[jindx[j - 1] + u + 1][cnt3];
                       cnt4 <= matrices->l_max_Q_M1[jindx[j - 1] + u + 1][cnt3];
                       cnt4 += 2) {
                    if (((cnt1 + cnt3 + da) <= maxD1) && ((cnt2 + cnt4 + db) <= maxD2)) {
                      matrices->Q_B[ij][cnt1 + cnt3 + da][(cnt2 + cnt4 + db) / 2] += matrices->Q_M[my_iindx[i + 1] - u][cnt1][cnt2 / 2]
                                                                                     * matrices->Q_M1[jindx[j - 1] + u + 1][cnt3][cnt4 / 2]
                                                                                     * temp2;
                      if (update_b) {
                        updatePosteriorBoundaries(cnt1 + cnt3 + da,
                                                  cnt2 + cnt4 + db,
                                                  &k_min_post_b,
                                                  &k_max_post_b,
                                                  &l_min_post_b,
                                                  &l_max_post_b
                                                  );
                      }
                    } else {
                      matrices->Q_B_rem[ij] += matrices->Q_M[my_iindx[i + 1] - u][cnt1][cnt2 / 2]
                                               * matrices->Q_M1[jindx[j - 1] + u + 1][cnt3][cnt4 / 2]
                                               * temp2;
                    }
                  }
              }
          }
        }

        if (update_b) {
          adjustArrayBoundaries(&matrices->Q_B[ij],
                                &matrices->k_min_Q_B[ij],
                                &matrices->k_max_Q_B[ij],
                                &matrices->l_min_Q_B[ij],
                                &matrices->l_max_Q_B[ij],
                                k_min_post_b,
                                k_max_post_b,
                                l_min_post_b,
                                l_max_post_b
                                );
        }
      } /* end >> if (pair) << */

      /* free ends ? -----------------------------------------*/

      int k_min_Q_M, k_max_Q_M, l_min_Q_M, l_max_Q_M;
      int k_min_post_m, k_max_post_m, *l_min_post_m, *l_max_post_m;
      int update_m = 0;
      int k_min_Q_M1, k_max_Q_M1, l_min_Q_M1, l_max_Q_M1;
      int k_min_post_m1, k_max_post_m1, *l_min_post_m1, *l_max_post_m1;
      int update_m1 = 0;

      if (!matrices->Q_M[ij]) {
        update_m  = 1;
        k_min_Q_M = l_min_Q_M = 0;
        k_max_Q_M = mm1[ij] + referenceBPs1[ij];
        l_max_Q_M = mm2[ij] + referenceBPs2[ij];

        prepareBoundaries(k_min_Q_M,
                          k_max_Q_M,
                          l_min_Q_M,
                          l_max_Q_M,
                          bpdist[ij],
                          &matrices->k_min_Q_M[ij],
                          &matrices->k_max_Q_M[ij],
                          &matrices->l_min_Q_M[ij],
                          &matrices->l_max_Q_M[ij]
                          );
        preparePosteriorBoundaries(matrices->k_max_Q_M[ij] - matrices->k_min_Q_M[ij] + 1,
                                   matrices->k_min_Q_M[ij],
                                   &k_min_post_m,
                                   &k_max_post_m,
                                   &l_min_post_m,
                                   &l_max_post_m
                                   );

        prepareArray(&matrices->Q_M[ij],
                     matrices->k_min_Q_M[ij],
                     matrices->k_max_Q_M[ij],
                     matrices->l_min_Q_M[ij],
                     matrices->l_max_Q_M[ij]
                     );
      }

      if (!matrices->Q_M1[jindx[j] + i]) {
        update_m1   = 1;
        k_min_Q_M1  = l_min_Q_M1 = 0;
        k_max_Q_M1  = mm1[ij] + referenceBPs1[ij];
        l_max_Q_M1  = mm2[ij] + referenceBPs2[ij];

        prepareBoundaries(k_min_Q_M1,
                          k_max_Q_M1,
                          l_min_Q_M1,
                          l_max_Q_M1,
                          bpdist[ij],
                          &matrices->k_min_Q_M1[jindx[j] + i],
                          &matrices->k_max_Q_M1[jindx[j] + i],
                          &matrices->l_min_Q_M1[jindx[j] + i],
                          &matrices->l_max_Q_M1[jindx[j] + i]
                          );
        preparePosteriorBoundaries(matrices->k_max_Q_M1[jindx[j] + i] - matrices->k_min_Q_M1[jindx[j] + i] + 1,
                                   matrices->k_min_Q_M1[jindx[j] + i],
                                   &k_min_post_m1,
                                   &k_max_post_m1,
                                   &l_min_post_m1,
                                   &l_max_post_m1
                                   );

        prepareArray(&matrices->Q_M1[jindx[j] + i],
                     matrices->k_min_Q_M1[jindx[j] + i],
                     matrices->k_max_Q_M1[jindx[j] + i],
                     matrices->l_min_Q_M1[jindx[j] + i],
                     matrices->l_max_Q_M1[jindx[j] + i]
                     );
      }

      /* j is unpaired */
      da  = referenceBPs1[ij] - referenceBPs1[ij + 1];
      db  = referenceBPs2[ij] - referenceBPs2[ij + 1];

      if (matrices->Q_M_rem[ij + 1])
        matrices->Q_M_rem[ij] += matrices->Q_M_rem[ij + 1] * pf_params->expMLbase * scale[1];

      if (matrices->Q_M[ij + 1]) {
        for (cnt1 = matrices->k_min_Q_M[ij + 1];
             cnt1 <= matrices->k_max_Q_M[ij + 1];
             cnt1++) {
          for (cnt2 = matrices->l_min_Q_M[ij + 1][cnt1];
               cnt2 <= matrices->l_max_Q_M[ij + 1][cnt1];
               cnt2 += 2) {
            if (((cnt1 + da) <= maxD1) && ((cnt2 + db) <= maxD2)) {
              matrices->Q_M[ij][cnt1 + da][(cnt2 + db) / 2] += matrices->Q_M[ij + 1][cnt1][cnt2 / 2] * pf_params->expMLbase * scale[1];
              if (update_m) {
                updatePosteriorBoundaries(cnt1 + da,
                                          cnt2 + db,
                                          &k_min_post_m,
                                          &k_max_post_m,
                                          &l_min_post_m,
                                          &l_max_post_m
                                          );
              }
            } else {
              matrices->Q_M_rem[ij] += matrices->Q_M[ij + 1][cnt1][cnt2 / 2] * pf_params->expMLbase * scale[1];
            }
          }
        }
      }

      if (matrices->Q_M1_rem[jindx[j - 1] + i])
        matrices->Q_M1_rem[jindx[j] + i] += matrices->Q_M1_rem[jindx[j - 1] + i] * pf_params->expMLbase * scale[1];

      if (matrices->Q_M1[jindx[j - 1] + i]) {
        for (cnt1 = matrices->k_min_Q_M1[jindx[j - 1] + i];
             cnt1 <= matrices->k_max_Q_M1[jindx[j - 1] + i];
             cnt1++)
          for (cnt2 = matrices->l_min_Q_M1[jindx[j - 1] + i][cnt1];
               cnt2 <= matrices->l_max_Q_M1[jindx[j - 1] + i][cnt1];
               cnt2 += 2) {
            if (((cnt1 + da) <= maxD1) && ((cnt2 + db) <= maxD2)) {
              matrices->Q_M1[jindx[j] + i][cnt1 + da][(cnt2 + db) / 2] += matrices->Q_M1[jindx[j - 1] + i][cnt1][cnt2 / 2] * pf_params->expMLbase * scale[1];
              if (update_m1) {
                updatePosteriorBoundaries(cnt1 + da,
                                          cnt2 + db,
                                          &k_min_post_m1,
                                          &k_max_post_m1,
                                          &l_min_post_m1,
                                          &l_max_post_m1
                                          );
              }
            } else {
              matrices->Q_M1_rem[jindx[j] + i] += matrices->Q_M1[jindx[j - 1] + i][cnt1][cnt2 / 2] * pf_params->expMLbase * scale[1];
            }
          }
      }

      /* j pairs with i */
      if ((!no_close) && type) {
        FLT_OR_DBL aux_en = exp_E_MLstem(type, (i > 1) || circ ? S1[i - 1] : -1, (j < seq_length) || circ ? S1[j + 1] : -1, pf_params);

        if (matrices->Q_B_rem[ij]) {
          matrices->Q_M_rem[ij]             += matrices->Q_B_rem[ij] * aux_en;
          matrices->Q_M1_rem[jindx[j] + i]  += matrices->Q_B_rem[ij] * aux_en;
        }

        if (matrices->Q_B[ij]) {
          for (cnt1 = matrices->k_min_Q_B[ij];
               cnt1 <= matrices->k_max_Q_B[ij];
               cnt1++)
            for (cnt2 = matrices->l_min_Q_B[ij][cnt1];
                 cnt2 <= matrices->l_max_Q_B[ij][cnt1];
                 cnt2 += 2) {
              matrices->Q_M[ij][cnt1][cnt2 / 2] += matrices->Q_B[ij][cnt1][cnt2 / 2] * aux_en;
              if (update_m) {
                updatePosteriorBoundaries(cnt1,
                                          cnt2,
                                          &k_min_post_m,
                                          &k_max_post_m,
                                          &l_min_post_m,
                                          &l_max_post_m
                                          );
              }

              matrices->Q_M1[jindx[j] + i][cnt1][cnt2 / 2] += matrices->Q_B[ij][cnt1][cnt2 / 2] * aux_en;
              if (update_m1) {
                updatePosteriorBoundaries(cnt1,
                                          cnt2,
                                          &k_min_post_m1,
                                          &k_max_post_m1,
                                          &l_min_post_m1,
                                          &l_max_post_m1
                                          );
              }
            }
        }
      }

      /* j pairs with k: i<k<j */
      ii = my_iindx[i];
      for (k = i + 1; k <= j; k++) {
        tt    = ptype[jindx[j] + k];
        temp2 = exp_E_MLstem(tt, S1[k - 1], (j < seq_length) || circ ? S1[j + 1] : -1, pf_params);

        if (matrices->Q_B_rem[my_iindx[k] - j]) {
          matrices->Q_M_rem[ij] += matrices->Q_B_rem[my_iindx[k] - j] * pow(pf_params->expMLbase, (double)(k - i)) * scale[k - i] * temp2;
          if (matrices->Q_M[ii - k + 1]) {
            for (cnt1 = matrices->k_min_Q_M[ii - k + 1];
                 cnt1 <= matrices->k_max_Q_M[ii - k + 1];
                 cnt1++)
              for (cnt2 = matrices->l_min_Q_M[ii - k + 1][cnt1];
                   cnt2 <= matrices->l_max_Q_M[ii - k + 1][cnt1];
                   cnt2 += 2)
                matrices->Q_M_rem[ij] += matrices->Q_M[ii - k + 1][cnt1][cnt2 / 2] * matrices->Q_B_rem[my_iindx[k] - j] * temp2;
          }

          if (matrices->Q_M_rem[ii - k + 1])
            matrices->Q_M_rem[ij] += matrices->Q_M_rem[ii - k + 1] * matrices->Q_B_rem[my_iindx[k] - j] * temp2;
        }

        if (matrices->Q_M_rem[ii - k + 1]) {
          if (matrices->Q_B[my_iindx[k] - j]) {
            for (cnt1 = matrices->k_min_Q_B[my_iindx[k] - j];
                 cnt1 <= matrices->k_max_Q_B[my_iindx[k] - j];
                 cnt1++)
              for (cnt2 = matrices->l_min_Q_B[my_iindx[k] - j][cnt1];
                   cnt2 <= matrices->l_max_Q_B[my_iindx[k] - j][cnt1];
                   cnt2 += 2)
                matrices->Q_M_rem[ij] += matrices->Q_M_rem[my_iindx[k] - j] * matrices->Q_B[my_iindx[k] - j][cnt1][cnt2 / 2] * temp2;
          }
        }

        /* add contributions of QM(i,k-1)*QB(k,j)*e^b and
         *  e^((k-i) * c) * QB(k,j) * e^b
         *  therefor we need d1a = dbp(T1_{i,j}, T1_{i,k-1} + T1_{k,j}),
         *  d1b = dbp(T2_{i,j}, T2_{i,k-1} + T2_{k,j})
         *  d1c = dbp(T1_{i,j}, T1_{k,j})circ = 0;
         *  d1d = dbp(T2_{i,j}, T2_{k,j})
         */
        da  = referenceBPs1[ij] - referenceBPs1[my_iindx[k] - j];
        db  = referenceBPs2[ij] - referenceBPs2[my_iindx[k] - j];

        if (!matrices->Q_B[my_iindx[k] - j])
          continue;

        for (cnt1 = matrices->k_min_Q_B[my_iindx[k] - j];
             cnt1 <= matrices->k_max_Q_B[my_iindx[k] - j];
             cnt1++)
          for (cnt2 = matrices->l_min_Q_B[my_iindx[k] - j][cnt1];
               cnt2 <= matrices->l_max_Q_B[my_iindx[k] - j][cnt1];
               cnt2 += 2) {
            if (((cnt1 + da) <= maxD1) && ((cnt2 + db) <= maxD2)) {
              matrices->Q_M[ij][cnt1 + da][(cnt2 + db) / 2] += matrices->Q_B[my_iindx[k] - j][cnt1][cnt2 / 2] * pow(pf_params->expMLbase, (double)(k - i)) * scale[k - i] * temp2;
              if (update_m) {
                updatePosteriorBoundaries(cnt1 + da,
                                          cnt2 + db,
                                          &k_min_post_m,
                                          &k_max_post_m,
                                          &l_min_post_m,
                                          &l_max_post_m
                                          );
              }
            } else {
              matrices->Q_M_rem[ij] += matrices->Q_B[my_iindx[k] - j][cnt1][cnt2 / 2] * pow(pf_params->expMLbase, (double)(k - i)) * scale[k - i] * temp2;
            }
          }

        if (!matrices->Q_M[ii - k + 1])
          continue;

        da  -= referenceBPs1[ii - k + 1];
        db  -= referenceBPs2[ii - k + 1];

        for (cnt1 = matrices->k_min_Q_M[ii - k + 1];
             cnt1 <= matrices->k_max_Q_M[ii - k + 1];
             cnt1++)
          for (cnt2 = matrices->l_min_Q_M[ii - k + 1][cnt1];
               cnt2 <= matrices->l_max_Q_M[ii - k + 1][cnt1];
               cnt2 += 2)
            for (cnt3 = matrices->k_min_Q_B[my_iindx[k] - j];
                 cnt3 <= matrices->k_max_Q_B[my_iindx[k] - j];
                 cnt3++)
              for (cnt4 = matrices->l_min_Q_B[my_iindx[k] - j][cnt3];
                   cnt4 <= matrices->l_max_Q_B[my_iindx[k] - j][cnt3];
                   cnt4 += 2) {
                if (((cnt1 + cnt3 + da) <= maxD1) && ((cnt2 + cnt4 + db) <= maxD2)) {
                  matrices->Q_M[ij][cnt1 + cnt3 + da][(cnt2 + cnt4 + db) / 2] += matrices->Q_M[ii - k + 1][cnt1][cnt2 / 2] * matrices->Q_B[my_iindx[k] - j][cnt3][cnt4 / 2] * temp2;
                  if (update_m) {
                    updatePosteriorBoundaries(cnt1 + cnt3 + da,
                                              cnt2 + cnt4 + db,
                                              &k_min_post_m,
                                              &k_max_post_m,
                                              &l_min_post_m,
                                              &l_max_post_m
                                              );
                  }
                } else {
                  matrices->Q_M_rem[ij] += matrices->Q_M[ii - k + 1][cnt1][cnt2 / 2] * matrices->Q_B[my_iindx[k] - j][cnt3][cnt4 / 2] * temp2;
                }
              }
      }

      if (update_m) {
        adjustArrayBoundaries(&matrices->Q_M[ij],
                              &matrices->k_min_Q_M[ij],
                              &matrices->k_max_Q_M[ij],
                              &matrices->l_min_Q_M[ij],
                              &matrices->l_max_Q_M[ij],
                              k_min_post_m,
                              k_max_post_m,
                              l_min_post_m,
                              l_max_post_m
                              );
      }

      if (update_m1) {
        adjustArrayBoundaries(&matrices->Q_M1[jindx[j] + i],
                              &matrices->k_min_Q_M1[jindx[j] + i],
                              &matrices->k_max_Q_M1[jindx[j] + i],
                              &matrices->l_min_Q_M1[jindx[j] + i],
                              &matrices->l_max_Q_M1[jindx[j] + i],
                              k_min_post_m1,
                              k_max_post_m1,
                              l_min_post_m1,
                              l_max_post_m1
                              );
      }

      /* compute contributions for Q(i,j) */
      int k_min, k_max, l_min, l_max;
      int k_min_post, k_max_post, *l_min_post, *l_max_post;
      int update_q = 0;
      if (!matrices->Q[ij]) {
        update_q  = 1;
        k_min     = l_min = 0;
        k_max     = mm1[ij] + referenceBPs1[ij];
        l_max     = mm2[ij] + referenceBPs2[ij];

        prepareBoundaries(k_min,
                          k_max,
                          l_min,
                          l_max,
                          bpdist[ij],
                          &matrices->k_min_Q[ij],
                          &matrices->k_max_Q[ij],
                          &matrices->l_min_Q[ij],
                          &matrices->l_max_Q[ij]
                          );
        preparePosteriorBoundaries(matrices->k_max_Q[ij] - matrices->k_min_Q[ij] + 1,
                                   matrices->k_min_Q[ij],
                                   &k_min_post,
                                   &k_max_post,
                                   &l_min_post,
                                   &l_max_post
                                   );

        prepareArray(&matrices->Q[ij],
                     matrices->k_min_Q[ij],
                     matrices->k_max_Q[ij],
                     matrices->l_min_Q[ij],
                     matrices->l_max_Q[ij]
                     );
      }

      if (type) {
        aux_en = exp_E_ExtLoop(type, (i > 1) || circ ? S1[i - 1] : -1, (j < seq_length) || circ ? S1[j + 1] : -1, pf_params);

        if (matrices->Q_B_rem[ij])
          matrices->Q_rem[ij] += matrices->Q_B_rem[ij] * aux_en;

        if (matrices->Q_B[ij]) {
          for (cnt1 = matrices->k_min_Q_B[ij];
               cnt1 <= matrices->k_max_Q_B[ij];
               cnt1++)
            for (cnt2 = matrices->l_min_Q_B[ij][cnt1];
                 cnt2 <= matrices->l_max_Q_B[ij][cnt1];
                 cnt2 += 2) {
              matrices->Q[ij][cnt1][cnt2 / 2] += matrices->Q_B[ij][cnt1][cnt2 / 2] * aux_en;
              if (update_q) {
                updatePosteriorBoundaries(cnt1,
                                          cnt2,
                                          &k_min_post,
                                          &k_max_post,
                                          &l_min_post,
                                          &l_max_post
                                          );
              }
            }
        }
      }

      /* j is unpaired */
      if (matrices->Q_rem[ij + 1])
        matrices->Q_rem[ij] += matrices->Q_rem[ij + 1] * scale[1];

      /* da = dbp(T1_{i,j}, T1_{i,j-1})
       *  db = dbp(T2_{i,j}, T2_{i,j-1})
       */
      da  = referenceBPs1[ij] - referenceBPs1[ij + 1];
      db  = referenceBPs2[ij] - referenceBPs2[ij + 1];
      if (matrices->Q[ij + 1]) {
        for (cnt1 = matrices->k_min_Q[ij + 1];
             cnt1 <= matrices->k_max_Q[ij + 1];
             cnt1++)
          for (cnt2 = matrices->l_min_Q[ij + 1][cnt1];
               cnt2 <= matrices->l_max_Q[ij + 1][cnt1];
               cnt2 += 2) {
            if (((cnt1 + da) <= maxD1) && ((cnt2 + db) <= maxD2)) {
              matrices->Q[ij][cnt1 + da][(cnt2 + db) / 2] += matrices->Q[ij + 1][cnt1][cnt2 / 2] * scale[1];
              if (update_q) {
                updatePosteriorBoundaries(cnt1 + da,
                                          cnt2 + db,
                                          &k_min_post,
                                          &k_max_post,
                                          &l_min_post,
                                          &l_max_post
                                          );
              }
            } else {
              matrices->Q_rem[ij] += matrices->Q[ij + 1][cnt1][cnt2 / 2] * scale[1];
            }
          }
      }

      for (k = j - TURN - 1; k > i; k--) {
        tt    = ptype[jindx[j] + k];
        temp2 = exp_E_ExtLoop(tt, S1[k - 1], (j < seq_length) || circ ? S1[j + 1] : -1, pf_params);

        if (matrices->Q_rem[my_iindx[i] - k + 1]) {
          if (matrices->Q_B[my_iindx[k] - j]) {
            for (cnt1 = matrices->k_min_Q_B[my_iindx[k] - j];
                 cnt1 <= matrices->k_max_Q_B[my_iindx[k] - j];
                 cnt1++)
              for (cnt2 = matrices->l_min_Q_B[my_iindx[k] - j][cnt1];
                   cnt2 <= matrices->l_max_Q_B[my_iindx[k] - j][cnt1];
                   cnt2 += 2)
                matrices->Q_rem[ij] += matrices->Q_rem[my_iindx[i] - k + 1] * matrices->Q_B[my_iindx[k] - j][cnt1][cnt2 / 2] * temp2;
          }

          if (matrices->Q_B_rem[my_iindx[k] - j])
            matrices->Q_rem[ij] += matrices->Q_rem[my_iindx[i] - k + 1] * matrices->Q_B_rem[my_iindx[k] - j] * temp2;
        }

        if (matrices->Q_B_rem[my_iindx[k] - j]) {
          if (matrices->Q[my_iindx[i] - k + 1]) {
            for (cnt1 = matrices->k_min_Q[my_iindx[i] - k + 1];
                 cnt1 <= matrices->k_max_Q[my_iindx[i] - k + 1];
                 cnt1++)
              for (cnt2 = matrices->l_min_Q[my_iindx[i] - k + 1][cnt1];
                   cnt2 <= matrices->l_max_Q[my_iindx[i] - k + 1][cnt1];
                   cnt2 += 2)
                matrices->Q_rem[ij] += matrices->Q[my_iindx[i] - k + 1][cnt1][cnt2 / 2] * matrices->Q_B_rem[my_iindx[k] - j] * temp2;
          }
        }

        /* da = dbp{T1_{i,j}, T1_{k,j}
         *  db = dbp{T2_{i,j}, T2_{k,j}}
         */
        da  = referenceBPs1[ij] - referenceBPs1[my_iindx[k] - j] - referenceBPs1[my_iindx[i] - k + 1];
        db  = referenceBPs2[ij] - referenceBPs2[my_iindx[k] - j] - referenceBPs2[my_iindx[i] - k + 1];


        if (!matrices->Q[my_iindx[i] - k + 1])
          continue;

        if (!matrices->Q_B[my_iindx[k] - j])
          continue;

        for (cnt1 = matrices->k_min_Q[my_iindx[i] - k + 1];
             cnt1 <= matrices->k_max_Q[my_iindx[i] - k + 1];
             cnt1++)
          for (cnt2 = matrices->l_min_Q[my_iindx[i] - k + 1][cnt1];
               cnt2 <= matrices->l_max_Q[my_iindx[i] - k + 1][cnt1];
               cnt2 += 2)
            for (cnt3 = matrices->k_min_Q_B[my_iindx[k] - j];
                 cnt3 <= matrices->k_max_Q_B[my_iindx[k] - j];
                 cnt3++)
              for (cnt4 = matrices->l_min_Q_B[my_iindx[k] - j][cnt3];
                   cnt4 <= matrices->l_max_Q_B[my_iindx[k] - j][cnt3];
                   cnt4 += 2) {
                if (((cnt1 + cnt3 + da) <= maxD1) && ((cnt2 + cnt4 + db) <= maxD2)) {
                  matrices->Q[ij][cnt1 + cnt3 + da][(cnt2 + cnt4 + db) / 2] += matrices->Q[my_iindx[i] - k + 1][cnt1][cnt2 / 2] * matrices->Q_B[my_iindx[k] - j][cnt3][cnt4 / 2] * temp2;
                  if (update_q) {
                    updatePosteriorBoundaries(cnt1 + cnt3 + da,
                                              cnt2 + cnt4 + db,
                                              &k_min_post,
                                              &k_max_post,
                                              &l_min_post,
                                              &l_max_post
                                              );
                  }
                } else {
                  matrices->Q_rem[ij] += matrices->Q[my_iindx[i] - k + 1][cnt1][cnt2 / 2] * matrices->Q_B[my_iindx[k] - j][cnt3][cnt4 / 2] * temp2;
                }
              }
      }

      if (update_q) {
        adjustArrayBoundaries(&matrices->Q[ij],
                              &matrices->k_min_Q[ij],
                              &matrices->k_max_Q[ij],
                              &matrices->l_min_Q[ij],
                              &matrices->l_max_Q[ij],
                              k_min_post,
                              k_max_post,
                              l_min_post,
                              l_max_post
                              );
      }

#if 1
      for (cnt1 = matrices->k_min_Q[ij];
           cnt1 <= matrices->k_max_Q[ij];
           cnt1++) {
        for (cnt2 = matrices->l_min_Q[ij][cnt1];
             cnt2 <= matrices->l_max_Q[ij][cnt1];
             cnt2 += 2) {
          if (matrices->Q[ij][cnt1][cnt2 / 2] > Qmax) {
            Qmax = matrices->Q[ij][cnt1][cnt2 / 2];
            if (Qmax > max_real / 10.)
              vrna_message_warning("Q close to overflow: %u %u %g\n", i, j, matrices->Q[ij][cnt1][cnt2 / 2]);
          }

          if (matrices->Q[ij][cnt1][cnt2 / 2] >= max_real)
            vrna_message_error("overflow in pf_fold while calculating q[%u,%u]\n"
                               "use larger pf_scale", i, j);
        }
      }
#endif
    } /* end of j-loop */
  }
}


/* calculate partition function for circular case */
/* NOTE: this is the postprocessing step ONLY     */
/* You have to call pf2D_linear first to calculate  */
/* complete circular case!!!                      */
PRIVATE void
pf2D_circ(vrna_fold_compound_t *vc)
{
  unsigned int      d, p, q, pq, k, l, kl, u, da, db, seq_length, maxD1, maxD2, base_d1, base_d2, *mm1, *mm2, *bpdist;
  int               *my_iindx, *jindx, type, cnt1, cnt2, cnt3, cnt4, *rtype;
  short             *S1;
  unsigned int      *referenceBPs1, *referenceBPs2;
  char              *sequence, *ptype;
  FLT_OR_DBL        *scale;
  vrna_exp_param_t  *pf_params;     /* holds all [unscaled] pf parameters */
  vrna_md_t         *md;
  vrna_mx_pf_t      *matrices;

  pf_params     = vc->exp_params;
  md            = &(pf_params->model_details);
  matrices      = vc->exp_matrices;
  sequence      = vc->sequence;
  seq_length    = vc->length;
  maxD1         = vc->maxD1;
  maxD2         = vc->maxD2;
  S1            = vc->sequence_encoding;
  ptype         = vc->ptype;
  rtype         = &(md->rtype[0]);
  scale         = matrices->scale;
  my_iindx      = vc->iindx;
  jindx         = vc->jindx;
  referenceBPs1 = vc->referenceBPs1;
  referenceBPs2 = vc->referenceBPs2;
  dangles       = md->dangles;
  mm1           = vc->mm1;
  mm2           = vc->mm2;
  bpdist        = vc->bpdist;

  FLT_OR_DBL  ***Q_B, ***Q_M, ***Q_M1;
  FLT_OR_DBL  *Q_B_rem, *Q_M_rem, *Q_M1_rem;
  int         **l_min_Q_B, **l_max_Q_B, **l_min_Q_M, **l_max_Q_M, **l_min_Q_M1, **l_max_Q_M1;
  int         *k_min_Q_B, *k_max_Q_B, *k_min_Q_M, *k_max_Q_M, *k_min_Q_M1, *k_max_Q_M1;

  Q_B       = matrices->Q_B;
  l_min_Q_B = matrices->l_min_Q_B;
  l_max_Q_B = matrices->l_max_Q_B;
  k_min_Q_B = matrices->k_min_Q_B;
  k_max_Q_B = matrices->k_max_Q_B;

  Q_M       = matrices->Q_M;
  l_min_Q_M = matrices->l_min_Q_M;
  l_max_Q_M = matrices->l_max_Q_M;
  k_min_Q_M = matrices->k_min_Q_M;
  k_max_Q_M = matrices->k_max_Q_M;

  Q_M1        = matrices->Q_M1;
  l_min_Q_M1  = matrices->l_min_Q_M1;
  l_max_Q_M1  = matrices->l_max_Q_M1;
  k_min_Q_M1  = matrices->k_min_Q_M1;
  k_max_Q_M1  = matrices->k_max_Q_M1;


  Q_B_rem   = matrices->Q_B_rem;
  Q_M_rem   = matrices->Q_M_rem;
  Q_M1_rem  = matrices->Q_M1_rem;

  matrices->Q_c_rem   = 0.;
  matrices->Q_cH_rem  = 0.;
  matrices->Q_cI_rem  = 0.;
  matrices->Q_cM_rem  = 0.;


  /* construct qm2 matrix from qm1 entries  */
#ifdef _OPENMP
#pragma omp parallel for private(d, k, l, da, db, cnt1, cnt2, cnt3, cnt4)
#endif
  for (k = 1; k < seq_length - TURN - 1; k++) {
    int k_min_Q_M2, k_max_Q_M2, l_min_Q_M2, l_max_Q_M2;
    int k_min_post_m2, k_max_post_m2, *l_min_post_m2, *l_max_post_m2;
    int update_m2 = 0;
    if (!matrices->Q_M2[k]) {
      update_m2   = 1;
      k_min_Q_M2  = l_min_Q_M2 = 0;
      k_max_Q_M2  = mm1[my_iindx[k] - seq_length] + referenceBPs1[my_iindx[k] - seq_length];
      l_max_Q_M2  = mm2[my_iindx[k] - seq_length] + referenceBPs2[my_iindx[k] - seq_length];

      prepareBoundaries(k_min_Q_M2,
                        k_max_Q_M2,
                        l_min_Q_M2,
                        l_max_Q_M2,
                        bpdist[my_iindx[k] - seq_length],
                        &matrices->k_min_Q_M2[k],
                        &matrices->k_max_Q_M2[k],
                        &matrices->l_min_Q_M2[k],
                        &matrices->l_max_Q_M2[k]
                        );
      preparePosteriorBoundaries(matrices->k_max_Q_M2[k] - matrices->k_min_Q_M2[k] + 1,
                                 matrices->k_min_Q_M2[k],
                                 &k_min_post_m2,
                                 &k_max_post_m2,
                                 &l_min_post_m2,
                                 &l_max_post_m2
                                 );

      prepareArray(&matrices->Q_M2[k],
                   matrices->k_min_Q_M2[k],
                   matrices->k_max_Q_M2[k],
                   matrices->l_min_Q_M2[k],
                   matrices->l_max_Q_M2[k]
                   );
    }

    /* construct Q_M2 */
    for (l = k + TURN + 1; l < seq_length - TURN - 1; l++) {
      if (Q_M1_rem[jindx[l] + k]) {
        if (Q_M1[jindx[seq_length] + l + 1]) {
          for (cnt1 = k_min_Q_M1[jindx[seq_length] + l + 1];
               cnt1 <= k_max_Q_M1[jindx[seq_length] + l + 1];
               cnt1++)
            for (cnt2 = l_min_Q_M1[jindx[seq_length] + l + 1][cnt1];
                 cnt2 <= l_max_Q_M1[jindx[seq_length] + l + 1][cnt1];
                 cnt2 += 2)
              matrices->Q_M2_rem[k] += Q_M1_rem[jindx[l] + k] * Q_M1[jindx[seq_length] + l + 1][cnt1][cnt2 / 2];
        }

        if (Q_M1_rem[jindx[seq_length] + l + 1])
          matrices->Q_M2_rem[k] += Q_M1_rem[jindx[l] + k] * Q_M1_rem[jindx[seq_length] + l + 1];
      }

      if (Q_M1_rem[jindx[seq_length] + l + 1]) {
        if (Q_M1[jindx[l] + k]) {
          for (cnt1 = k_min_Q_M1[jindx[l] + k];
               cnt1 <= k_max_Q_M1[jindx[l] + k];
               cnt1++)
            for (cnt2 = l_min_Q_M1[jindx[l] + k][cnt1];
                 cnt2 <= l_max_Q_M1[jindx[l] + k][cnt1];
                 cnt2 += 2)
              matrices->Q_M2_rem[k] += Q_M1[jindx[l] + k][cnt1][cnt2 / 2] * Q_M1_rem[jindx[seq_length] + l + 1];
        }
      }

      if (matrices->Q_M1[jindx[l] + k] && matrices->Q_M1[jindx[seq_length] + l + 1]) {
        da  = referenceBPs1[my_iindx[k] - seq_length] - referenceBPs1[my_iindx[k] - l] - referenceBPs1[my_iindx[l + 1] - seq_length];
        db  = referenceBPs2[my_iindx[k] - seq_length] - referenceBPs2[my_iindx[k] - l] - referenceBPs2[my_iindx[l + 1] - seq_length];
        for (cnt1 = k_min_Q_M1[jindx[l] + k]; cnt1 <= k_max_Q_M1[jindx[l] + k]; cnt1++)
          for (cnt2 = l_min_Q_M1[jindx[l] + k][cnt1]; cnt2 <= l_max_Q_M1[jindx[l] + k][cnt1]; cnt2 += 2) {
            for (cnt3 = k_min_Q_M1[jindx[seq_length] + l + 1]; cnt3 <= k_max_Q_M1[jindx[seq_length] + l + 1]; cnt3++)
              for (cnt4 = l_min_Q_M1[jindx[seq_length] + l + 1][cnt3]; cnt4 <= l_max_Q_M1[jindx[seq_length] + l + 1][cnt3]; cnt4 += 2) {
                if (((cnt1 + cnt3 + da) <= maxD1) && ((cnt2 + cnt4 + db) <= maxD2)) {
                  matrices->Q_M2[k][cnt1 + cnt3 + da][(cnt2 + cnt4 + db) / 2] += Q_M1[jindx[l] + k][cnt1][cnt2 / 2] * Q_M1[jindx[seq_length] + l + 1][cnt3][cnt4 / 2];
                  if (update_m2) {
                    updatePosteriorBoundaries(cnt1 + cnt3 + da,
                                              cnt2 + cnt4 + db,
                                              &k_min_post_m2,
                                              &k_max_post_m2,
                                              &l_min_post_m2,
                                              &l_max_post_m2
                                              );
                  }
                } else {
                  matrices->Q_M2_rem[k] += Q_M1[jindx[l] + k][cnt1][cnt2 / 2] * Q_M1[jindx[seq_length] + l + 1][cnt3][cnt4 / 2];
                }
              }
          }
      }
    }
    if (update_m2) {
      adjustArrayBoundaries(&matrices->Q_M2[k],
                            &matrices->k_min_Q_M2[k],
                            &matrices->k_max_Q_M2[k],
                            &matrices->l_min_Q_M2[k],
                            &matrices->l_max_Q_M2[k],
                            k_min_post_m2,
                            k_max_post_m2,
                            l_min_post_m2,
                            l_max_post_m2
                            );
    }
  }

  base_d1 = referenceBPs1[my_iindx[1] - seq_length];
  base_d2 = referenceBPs2[my_iindx[1] - seq_length];

  int min_k, max_k, max_l, min_l;
  int min_k_real, max_k_real, min_k_real_qcH, max_k_real_qcH, min_k_real_qcI, max_k_real_qcI, min_k_real_qcM, max_k_real_qcM;
  int *min_l_real, *max_l_real, *min_l_real_qcH, *max_l_real_qcH, *min_l_real_qcI, *max_l_real_qcI, *min_l_real_qcM, *max_l_real_qcM;
  int update_c, update_cH, update_cI, update_cM;

  update_c = update_cH = update_cI = update_cM = 0;

  min_k = min_l = 0;

  max_k = mm1[my_iindx[1] - seq_length] + referenceBPs1[my_iindx[1] - seq_length];
  max_l = mm2[my_iindx[1] - seq_length] + referenceBPs2[my_iindx[1] - seq_length];

#ifdef _OPENMP
#pragma omp sections
  {
#pragma omp section
    {
#endif
  if (!matrices->Q_c) {
    update_c = 1;
    prepareBoundaries(min_k,
                      max_k,
                      min_l,
                      max_l,
                      bpdist[my_iindx[1] - seq_length],
                      &matrices->k_min_Q_c,
                      &matrices->k_max_Q_c,
                      &matrices->l_min_Q_c,
                      &matrices->l_max_Q_c
                      );
    prepareArray(&matrices->Q_c,
                 matrices->k_min_Q_c,
                 matrices->k_max_Q_c,
                 matrices->l_min_Q_c,
                 matrices->l_max_Q_c
                 );
    preparePosteriorBoundaries(max_k - min_k + 1,
                               min_k,
                               &min_k_real,
                               &max_k_real,
                               &min_l_real,
                               &max_l_real
                               );
  }

#ifdef _OPENMP
}


#pragma omp section
{
#endif
  if (!matrices->Q_cH) {
    update_cH = 1;
    prepareBoundaries(min_k,
                      max_k,
                      min_l,
                      max_l,
                      bpdist[my_iindx[1] - seq_length],
                      &matrices->k_min_Q_cH,
                      &matrices->k_max_Q_cH,
                      &matrices->l_min_Q_cH,
                      &matrices->l_max_Q_cH
                      );
    prepareArray(&matrices->Q_cH,
                 matrices->k_min_Q_cH,
                 matrices->k_max_Q_cH,
                 matrices->l_min_Q_cH,
                 matrices->l_max_Q_cH
                 );
    preparePosteriorBoundaries(max_k - min_k + 1,
                               min_k,
                               &min_k_real_qcH,
                               &max_k_real_qcH,
                               &min_l_real_qcH,
                               &max_l_real_qcH
                               );
  }

#ifdef _OPENMP
}
#pragma omp section
{
#endif
  if (!matrices->Q_cI) {
    update_cI = 1;
    prepareBoundaries(min_k,
                      max_k,
                      min_l,
                      max_l,
                      bpdist[my_iindx[1] - seq_length],
                      &matrices->k_min_Q_cI,
                      &matrices->k_max_Q_cI,
                      &matrices->l_min_Q_cI,
                      &matrices->l_max_Q_cI
                      );
    prepareArray(&matrices->Q_cI,
                 matrices->k_min_Q_cI,
                 matrices->k_max_Q_cI,
                 matrices->l_min_Q_cI,
                 matrices->l_max_Q_cI
                 );
    preparePosteriorBoundaries(max_k - min_k + 1,
                               min_k,
                               &min_k_real_qcI,
                               &max_k_real_qcI,
                               &min_l_real_qcI,
                               &max_l_real_qcI
                               );
  }

#ifdef _OPENMP
}
#pragma omp section
{
#endif
  if (!matrices->Q_cM) {
    update_cM = 1;
    prepareBoundaries(min_k,
                      max_k,
                      min_l,
                      max_l,
                      bpdist[my_iindx[1] - seq_length],
                      &matrices->k_min_Q_cM,
                      &matrices->k_max_Q_cM,
                      &matrices->l_min_Q_cM,
                      &matrices->l_max_Q_cM
                      );
    prepareArray(&matrices->Q_cM,
                 matrices->k_min_Q_cM,
                 matrices->k_max_Q_cM,
                 matrices->l_min_Q_cM,
                 matrices->l_max_Q_cM
                 );
    preparePosteriorBoundaries(max_k - min_k + 1,
                               min_k,
                               &min_k_real_qcM,
                               &max_k_real_qcM,
                               &min_l_real_qcM,
                               &max_l_real_qcM
                               );
  }

#ifdef _OPENMP
}
}
#endif


  for (d = TURN + 2; d <= seq_length; d++) /* i,j in [1..length] */
#ifdef _OPENMP
#pragma omp parallel for private(p, q, pq, k, l, kl, u, da, db, type, cnt1, cnt2, cnt3, cnt4)
#endif
    for (q = d; q <= seq_length; q++) {
      FLT_OR_DBL  qot;
      char        loopseq[10];
      p   = q - d + 1;
      pq  = my_iindx[p] - q;

      /* 1. get exterior hairpin contribution  */
      u = seq_length - q + p - 1;
      if (u < TURN)
        continue;

      type = ptype[jindx[q] + p];
      if (!type)
        continue;

      if (((type == 3) || (type == 4)) && no_closingGU)
        continue;

      /* cause we want to calc the exterior loops, we need the reversed pair type from now on  */
      type = rtype[type];

      if (u < 7) {
        strcpy(loopseq, sequence + q - 1);
        strncat(loopseq, sequence, p);
      }

      /* get distance to reference if closing the hairpin
       *  da = dbp(T1_[1,n}, T1_{p,q})
       *  db = dbp(T2_{1,n}, T2_{p,q})
       */
      da  = base_d1 - referenceBPs1[pq];
      db  = base_d2 - referenceBPs2[pq];
      qot = exp_E_Hairpin(u, type, S1[q + 1], S1[p - 1], loopseq, pf_params) * scale[u];

      if (Q_B_rem[pq])
        matrices->Q_cH_rem += Q_B_rem[pq] * qot;

      if (Q_B[pq]) {
        for (cnt1 = k_min_Q_B[pq];
             cnt1 <= k_max_Q_B[pq];
             cnt1++)
          for (cnt2 = l_min_Q_B[pq][cnt1];
               cnt2 <= l_max_Q_B[pq][cnt1];
               cnt2 += 2) {
            if (((cnt1 + da) <= maxD1) && ((cnt2 + db) <= maxD2)) {
              matrices->Q_cH[cnt1 + da][(cnt2 + db) / 2] += Q_B[pq][cnt1][cnt2 / 2] * qot;
              if (update_cH) {
                updatePosteriorBoundaries(cnt1 + da,
                                          cnt2 + db,
                                          &min_k_real_qcH,
                                          &max_k_real_qcH,
                                          &min_l_real_qcH,
                                          &max_l_real_qcH
                                          );
              }
            } else {
              matrices->Q_cH_rem += Q_B[pq][cnt1][cnt2 / 2] * qot;
            }
          }
      }

      /* 2. exterior interior loops, i "define" the (k,l) pair as "outer pair"  */
      /* so "outer type" is rtype[type[k,l]] and inner type is type[p,q]        */
      if (Q_B_rem[pq]) {
        for (k = q + 1; k < seq_length; k++) {
          unsigned int ln1, lstart, ln_pre;
          ln1 = k - q - 1;
          if (ln1 + p - 1 > MAXLOOP)
            break;

          lstart  = k + TURN + 1;
          ln_pre  = ln1 + p + seq_length;
          if (ln_pre > lstart + MAXLOOP)
            lstart = ln_pre - MAXLOOP - 1;

          for (l = lstart; l <= seq_length; l++) {
            unsigned int  ln2;
            int           type2;
            kl  = my_iindx[k] - l;
            ln2 = (p - 1) + (seq_length - l);

            if ((ln1 + ln2) > MAXLOOP)
              continue;

            type2 = ptype[jindx[l] + k];
            if (!type2)
              continue;

            qot = exp_E_IntLoop(ln2, ln1, rtype[type2], type, S1[l + 1], S1[k - 1], S1[p - 1], S1[q + 1], pf_params) * scale[ln1 + ln2];

            if (Q_B_rem[kl])
              matrices->Q_cI_rem += Q_B_rem[pq] * Q_B_rem[kl] * qot;

            if (Q_B[kl]) {
              for (cnt1 = k_min_Q_B[kl];
                   cnt1 <= k_max_Q_B[kl];
                   cnt1++)
                for (cnt2 = l_min_Q_B[kl][cnt1];
                     cnt2 <= l_max_Q_B[kl][cnt1];
                     cnt2 += 2)
                  matrices->Q_cI_rem += Q_B_rem[pq] * Q_B[kl][cnt1][cnt2 / 2] * qot;
            }
          }
        }
      }

      if (Q_B[pq]) {
        for (k = q + 1; k < seq_length; k++) {
          unsigned int ln1, lstart, ln_pre;
          ln1 = k - q - 1;
          if (ln1 + p - 1 > MAXLOOP)
            break;

          lstart  = k + TURN + 1;
          ln_pre  = ln1 + p + seq_length;
          if (ln_pre > lstart + MAXLOOP)
            lstart = ln_pre - MAXLOOP - 1;

          for (l = lstart; l <= seq_length; l++) {
            unsigned int  ln2;
            int           type2;
            kl  = my_iindx[k] - l;
            ln2 = (p - 1) + (seq_length - l);

            if ((ln1 + ln2) > MAXLOOP)
              continue;

            type2 = ptype[jindx[l] + k];
            if (!type2)
              continue;

            qot = exp_E_IntLoop(ln2, ln1, rtype[type2], type, S1[l + 1], S1[k - 1], S1[p - 1], S1[q + 1], pf_params) * scale[ln1 + ln2];

            if (Q_B_rem[kl]) {
              for (cnt1 = k_min_Q_B[pq];
                   cnt1 <= k_max_Q_B[pq];
                   cnt1++)
                for (cnt2 = l_min_Q_B[pq][cnt1];
                     cnt2 <= l_max_Q_B[pq][cnt1];
                     cnt2 += 2)
                  matrices->Q_cI_rem += Q_B[pq][cnt1][cnt2 / 2] * Q_B_rem[kl] * qot;
            }

            if (!Q_B[kl])
              continue;

            /* get distance to reference if closing the interior loop
             *  d2a = dbp(T1_[1,n}, T1_{p,q} + T1_{k,l})
             *  d2b = dbp(T2_[1,n}, T2_{p,q} + T2_{k,l})
             */
            da  = base_d1 - referenceBPs1[pq] - referenceBPs1[kl];
            db  = base_d2 - referenceBPs2[pq] - referenceBPs2[kl];

            for (cnt1 = k_min_Q_B[pq]; cnt1 <= k_max_Q_B[pq]; cnt1++)
              for (cnt2 = l_min_Q_B[pq][cnt1]; cnt2 <= l_max_Q_B[pq][cnt1]; cnt2 += 2)
                for (cnt3 = k_min_Q_B[kl]; cnt3 <= k_max_Q_B[kl]; cnt3++)
                  for (cnt4 = l_min_Q_B[kl][cnt3]; cnt4 <= l_max_Q_B[kl][cnt3]; cnt4 += 2) {
                    if (((cnt1 + cnt3 + da) <= maxD1) && ((cnt2 + cnt4 + db) <= maxD2)) {
                      matrices->Q_cI[cnt1 + cnt3 + da][(cnt2 + cnt4 + db) / 2] += Q_B[pq][cnt1][cnt2 / 2] * Q_B[kl][cnt3][cnt4 / 2] * qot;
                      if (update_cI) {
                        updatePosteriorBoundaries(cnt1 + cnt3 + da,
                                                  cnt2 + cnt4 + db,
                                                  &min_k_real_qcI,
                                                  &max_k_real_qcI,
                                                  &min_l_real_qcI,
                                                  &max_l_real_qcI
                                                  );
                      }
                    } else {
                      matrices->Q_cI_rem += Q_B[pq][cnt1][cnt2 / 2] * Q_B[kl][cnt3][cnt4 / 2] * qot;
                    }
                  }
          }
        }
      }
    }

  if (update_cH) {
    adjustArrayBoundaries(&matrices->Q_cH,
                          &matrices->k_min_Q_cH,
                          &matrices->k_max_Q_cH,
                          &matrices->l_min_Q_cH,
                          &matrices->l_max_Q_cH,
                          min_k_real_qcH,
                          max_k_real_qcH,
                          min_l_real_qcH,
                          max_l_real_qcH
                          );
  }

  if (update_cI) {
    adjustArrayBoundaries(&matrices->Q_cI,
                          &matrices->k_min_Q_cI,
                          &matrices->k_max_Q_cI,
                          &matrices->l_min_Q_cI,
                          &matrices->l_max_Q_cI,
                          min_k_real_qcI,
                          max_k_real_qcI,
                          min_l_real_qcI,
                          max_l_real_qcI
                          );
  }

  /* 3. Multiloops  */
  if (seq_length > 2 * TURN - 3) {
#ifdef _OPENMP
#pragma omp parallel for private(k, da, db, cnt1, cnt2, cnt3, cnt4)
#endif
    for (k = TURN + 2; k < seq_length - 2 * TURN - 3; k++) {
      if (Q_M_rem[my_iindx[1] - k]) {
        if (matrices->Q_M2[k + 1]) {
          for (cnt1 = matrices->k_min_Q_M2[k + 1];
               cnt1 <= matrices->k_max_Q_M2[k + 1];
               cnt1++)
            for (cnt2 = matrices->l_min_Q_M2[k + 1][cnt1];
                 cnt2 <= matrices->l_max_Q_M2[k + 1][cnt1];
                 cnt2 += 2)
              matrices->Q_cM_rem += Q_M_rem[my_iindx[1] - k] * matrices->Q_M2[k + 1][cnt1][cnt2 / 2] * pf_params->expMLclosing;
        }

        if (matrices->Q_M2_rem[k + 1])
          matrices->Q_cM_rem += Q_M_rem[my_iindx[1] - k] * matrices->Q_M2_rem[k + 1] * pf_params->expMLclosing;
      }

      if (matrices->Q_M2_rem[k + 1]) {
        if (Q_M[my_iindx[1] - k]) {
          for (cnt1 = k_min_Q_M[my_iindx[1] - k];
               cnt1 <= k_max_Q_M[my_iindx[1] - k];
               cnt1++)
            for (cnt2 = l_min_Q_M[my_iindx[1] - k][cnt1];
                 cnt2 <= l_max_Q_M[my_iindx[1] - k][cnt1];
                 cnt2 += 2)
              matrices->Q_cM_rem += Q_M[my_iindx[1] - k][cnt1][cnt2 / 2] * matrices->Q_M2_rem[k + 1] * pf_params->expMLclosing;
        }
      }

      /* get distancies to references
       * d3a = dbp(T1_[1,n}, T1_{1,k} + T1_{k+1, n})
       * d3b = dbp(T2_[1,n}, T2_{1,k} + T2_{k+1, n})
       */
      da  = base_d1 - referenceBPs1[my_iindx[1] - k] - referenceBPs1[my_iindx[k + 1] - seq_length];
      db  = base_d2 - referenceBPs2[my_iindx[1] - k] - referenceBPs2[my_iindx[k + 1] - seq_length];
      if (Q_M[my_iindx[1] - k] && matrices->Q_M2[k + 1]) {
        for (cnt1 = k_min_Q_M[my_iindx[1] - k]; cnt1 <= k_max_Q_M[my_iindx[1] - k]; cnt1++)
          for (cnt2 = l_min_Q_M[my_iindx[1] - k][cnt1]; cnt2 <= l_max_Q_M[my_iindx[1] - k][cnt1]; cnt2 += 2)
            for (cnt3 = matrices->k_min_Q_M2[k + 1]; cnt3 <= matrices->k_max_Q_M2[k + 1]; cnt3++)
              for (cnt4 = matrices->l_min_Q_M2[k + 1][cnt3]; cnt4 <= matrices->l_max_Q_M2[k + 1][cnt3]; cnt4 += 2) {
                if (((cnt1 + cnt3 + da) <= maxD1) && ((cnt2 + cnt4 + db) <= maxD2)) {
                  matrices->Q_cM[cnt1 + cnt3 + da][(cnt2 + cnt4 + db) / 2] += Q_M[my_iindx[1] - k][cnt1][cnt2 / 2] * matrices->Q_M2[k + 1][cnt3][cnt4 / 2] * pf_params->expMLclosing;
                  if (update_cM) {
                    updatePosteriorBoundaries(cnt1 + cnt3 + da,
                                              cnt2 + cnt4 + db,
                                              &min_k_real_qcM,
                                              &max_k_real_qcM,
                                              &min_l_real_qcM,
                                              &max_l_real_qcM
                                              );
                  }
                } else {
                  matrices->Q_cM_rem += Q_M[my_iindx[1] - k][cnt1][cnt2 / 2] * matrices->Q_M2[k + 1][cnt3][cnt4 / 2] * pf_params->expMLclosing;
                }
              }
      }
    }
  }

  if (update_cM) {
    adjustArrayBoundaries(&matrices->Q_cM,
                          &matrices->k_min_Q_cM,
                          &matrices->k_max_Q_cM,
                          &matrices->l_min_Q_cM,
                          &matrices->l_max_Q_cM,
                          min_k_real_qcM,
                          max_k_real_qcM,
                          min_l_real_qcM,
                          max_l_real_qcM
                          );
  }

  for (cnt1 = matrices->k_min_Q_cH;
       cnt1 <= matrices->k_max_Q_cH;
       cnt1++)
    for (cnt2 = matrices->l_min_Q_cH[cnt1];
         cnt2 <= matrices->l_max_Q_cH[cnt1];
         cnt2 += 2) {
      matrices->Q_c[cnt1][cnt2 / 2] += matrices->Q_cH[cnt1][cnt2 / 2];
      if (update_c) {
        updatePosteriorBoundaries(cnt1,
                                  cnt2,
                                  &min_k_real,
                                  &max_k_real,
                                  &min_l_real,
                                  &max_l_real
                                  );
      }
    }
  for (cnt1 = matrices->k_min_Q_cI;
       cnt1 <= matrices->k_max_Q_cI;
       cnt1++)
    for (cnt2 = matrices->l_min_Q_cI[cnt1];
         cnt2 <= matrices->l_max_Q_cI[cnt1];
         cnt2 += 2) {
      matrices->Q_c[cnt1][cnt2 / 2] += matrices->Q_cI[cnt1][cnt2 / 2];
      if (update_c) {
        updatePosteriorBoundaries(cnt1,
                                  cnt2,
                                  &min_k_real,
                                  &max_k_real,
                                  &min_l_real,
                                  &max_l_real
                                  );
      }
    }
  for (cnt1 = matrices->k_min_Q_cM;
       cnt1 <= matrices->k_max_Q_cM;
       cnt1++)
    for (cnt2 = matrices->l_min_Q_cM[cnt1];
         cnt2 <= matrices->l_max_Q_cM[cnt1];
         cnt2 += 2) {
      matrices->Q_c[cnt1][cnt2 / 2] += matrices->Q_cM[cnt1][cnt2 / 2];
      if (update_c) {
        updatePosteriorBoundaries(cnt1,
                                  cnt2,
                                  &min_k_real,
                                  &max_k_real,
                                  &min_l_real,
                                  &max_l_real
                                  );
      }
    }

  matrices->Q_c_rem = matrices->Q_cH_rem + matrices->Q_cI_rem + matrices->Q_cM_rem;

  /* add the case were structure is unfolded chain */
  if ((referenceBPs1[my_iindx[1] - seq_length] <= maxD1) && (referenceBPs2[my_iindx[1] - seq_length] <= maxD2)) {
    matrices->Q_c[referenceBPs1[my_iindx[1] - seq_length]][referenceBPs2[my_iindx[1] - seq_length] / 2] += 1.0 * scale[seq_length];
    if (update_c) {
      updatePosteriorBoundaries(referenceBPs1[my_iindx[1] - seq_length],
                                referenceBPs2[my_iindx[1] - seq_length],
                                &min_k_real,
                                &max_k_real,
                                &min_l_real,
                                &max_l_real
                                );
    }
  } else {
    matrices->Q_c_rem += 1.0 * scale[seq_length];
  }

  adjustArrayBoundaries(&matrices->Q_c,
                        &matrices->k_min_Q_c,
                        &matrices->k_max_Q_c,
                        &matrices->l_min_Q_c,
                        &matrices->l_max_Q_c,
                        min_k_real,
                        max_k_real,
                        min_l_real,
                        max_l_real
                        );
}

/*
 * ###################################################
 * stochastic backtracking
 * ###################################################
 */
PUBLIC char *
vrna_pbacktrack_TwoD(vrna_fold_compound_t *vc,
                     int                  d1,
                     int                  d2)
{
  return vrna_pbacktrack5_TwoD(vc, d1, d2, vc->length);
}


PUBLIC char *
vrna_pbacktrack5_TwoD(vrna_fold_compound_t  *vc,
                      int                   d1,
                      int                   d2,
                      unsigned int          length)
{
  char *pstruc, *ptype;
  short *S1;
  unsigned int i, j, n, start, maxD1, maxD2, da, db,
               *referenceBPs1, *referenceBPs2;
  int *my_iindx, *jindx, ij, cnt1, cnt2, cnt3, cnt4, type,
      **l_min_Q, **l_max_Q,
      **l_min_Q_B, **l_max_Q_B,
      *k_min_Q, *k_max_Q,
      *k_min_Q_B, *k_max_Q_B;
  FLT_OR_DBL r, qt, *scale, ***Q, ***Q_B, *Q_rem, *Q_B_rem;
  vrna_exp_param_t *pf_params;
  vrna_md_t *md;
  vrna_mx_pf_t *matrices;

  n             = vc->length;
  pf_params     = vc->exp_params;
  md            = &(pf_params->model_details);
  matrices      = vc->exp_matrices;
  maxD1         = vc->maxD1;
  maxD2         = vc->maxD2;
  my_iindx      = vc->iindx;
  jindx         = vc->jindx;
  scale         = matrices->scale;
  ptype         = vc->ptype;
  S1            = vc->sequence_encoding;
  referenceBPs1 = vc->referenceBPs1;
  referenceBPs2 = vc->referenceBPs2;

  Q       = matrices->Q;
  l_min_Q = matrices->l_min_Q;
  l_max_Q = matrices->l_max_Q;
  k_min_Q = matrices->k_min_Q;
  k_max_Q = matrices->k_max_Q;

  Q_B       = matrices->Q_B;
  l_min_Q_B = matrices->l_min_Q_B;
  l_max_Q_B = matrices->l_max_Q_B;
  k_min_Q_B = matrices->k_min_Q_B;
  k_max_Q_B = matrices->k_max_Q_B;

  Q_rem   = matrices->Q_rem;
  Q_B_rem = matrices->Q_B_rem;

  if (md->circ) {
    if (n != length)
      vrna_message_error("vrna_pbacktrack_TwoD@2Dfold.c: cotranscriptional backtracking for circular RNAs not supported!");

    return pbacktrack_circ(vc, d1, d2);
  }

  if (length > n)
    vrna_message_error("vrna_pbacktrack_TwoD@2Dpfold.c: requested transcript length exceeds sequence length!");

#if 0
  if (d1 > maxD1)
    vrna_message_error("pbacktrack@2Dpfold.c: distance to 1st reference structure to high!");

  if (d2 > maxD2)
    vrna_message_error("pbacktrack@2Dpfold.c: distance to 2nd reference structure to high!");

#endif

  /* check whether the chosen neighborhood exists at all */
  int dumb = 1;
  ij = my_iindx[1] - length;
  if ((d1 == -1) && (Q_rem[ij] != 0.)) {
    dumb = 0;
  } else {
    if ((k_min_Q[ij] <= d1) && (k_max_Q[ij] >= d1)) {
      int l_min = l_min_Q[ij][d1];
      if ((d2 % 2) == (l_min % 2))
        if ((l_min <= d2) && (l_max_Q[ij][d1] >= d2))
          dumb = 0;
    }
  }

  if (dumb) {
    vrna_message_error("neighborhood %d:%d is not in scope of calculated partition function!\n"
                       "pbacktrack@2Dpfold.c: exiting...",
                       d1, d2);
  }

  pstruc = vrna_alloc((length + 1) * sizeof(char));

  for (i = 0; i < length; i++)
    pstruc[i] = '.';
  pstruc[i] = '\0';

  start = 1;
  while (start < length) {
    int sn = my_iindx[start] - length;
    /* find i position of first pair */
    FLT_OR_DBL qln_i = 0, qln_i1 = 0;

    if (d1 == -1) {
      qln_i = Q_rem[sn];

      /* open chain ? */
      if ((maxD1 > referenceBPs1[sn])
          && (maxD2 > referenceBPs2[sn])) {
        r = vrna_urn() * qln_i;
        if (scale[length - start + 1] > r)
          return pstruc;
      }

      /* lets see if we find a base pair with i involved */
      for (i = start; i < length; i++) {
        r = vrna_urn() * qln_i;

        qln_i1 = Q_rem[my_iindx[i + 1] - length];

        da  = referenceBPs1[sn] - referenceBPs1[my_iindx[i + 1] - length];
        db  = referenceBPs2[sn] - referenceBPs2[my_iindx[i + 1] - length];

        for (cnt1 = k_min_Q[my_iindx[i + 1] - length];
             cnt1 <= k_max_Q[my_iindx[i + 1] - length];
             cnt1++)
          for (cnt2 = l_min_Q[my_iindx[i + 1] - length][cnt1];
               cnt2 <= l_max_Q[my_iindx[i + 1] - length][cnt1];
               cnt2 += 2)
            if (((cnt1 + da) > maxD1) || ((cnt2 + db) > maxD2))
              qln_i1 += Q[my_iindx[i + 1] - length][cnt1][cnt2 / 2];

        if (r > qln_i1 * scale[1])
          break;

        qln_i = qln_i1;
      }
      if (i >= length)
        break;              /* no more pairs */

      /* i is paired, find pairing partner j */
      r = vrna_urn() * (qln_i - qln_i1 * scale[1]);
      for (qt = 0, j = i + TURN + 1; j < length; j++) {
        ij    = my_iindx[i] - j;
        type  = ptype[jindx[j] + i];
        if (type) {
          cnt1 = cnt2 = cnt3 = cnt4 = -1;
          double qkl = exp_E_ExtLoop(type, (i > 1) ? S1[i - 1] : -1, S1[j + 1], pf_params);

          if (Q_B_rem[ij] != 0.) {
            if (Q_rem[my_iindx[j + 1] - length] != 0.) {
              qt += qkl * Q_B_rem[ij] * Q_rem[my_iindx[j + 1] - length];
              if (qt >= r)
                goto pbacktrack_ext_loop_early_escape_rem;
            }

            if (Q[my_iindx[j + 1] - length]) {
              for (cnt3 = k_min_Q[my_iindx[j + 1] - length];
                   cnt3 <= k_max_Q[my_iindx[j + 1] - length];
                   cnt3++)
                for (cnt4 = l_min_Q[my_iindx[j + 1] - length][cnt3];
                     cnt4 <= l_max_Q[my_iindx[j + 1] - length][cnt3];
                     cnt4 += 2) {
                  qt += qkl * Q_B_rem[ij] * Q[my_iindx[j + 1] - length][cnt3][cnt4 / 2];
                  if (qt >= r)
                    goto pbacktrack_ext_loop_early_escape_rem;
                }
            }
          }

          if (Q_rem[my_iindx[j + 1] - length] != 0.) {
            cnt3 = cnt4 = -1;
            if (Q_B[ij]) {
              for (cnt1 = k_min_Q_B[ij];
                   cnt1 <= k_max_Q_B[ij];
                   cnt1++)
                for (cnt2 = l_min_Q_B[ij][cnt1];
                     cnt2 <= l_max_Q_B[ij][cnt1];
                     cnt2 += 2) {
                  qt += qkl * Q_B[ij][cnt1][cnt2 / 2] * Q_rem[my_iindx[j + 1] - length];
                  if (qt >= r)
                    goto pbacktrack_ext_loop_early_escape_rem;
                }
            }
          }

          /* if we still search for pairing partner j, we go on here... */
          if (Q_B[ij] && Q[my_iindx[j + 1] - length]) {
            da  = referenceBPs1[sn] - referenceBPs1[ij] - referenceBPs1[my_iindx[j + 1] - length];
            db  = referenceBPs2[sn] - referenceBPs2[ij] - referenceBPs2[my_iindx[j + 1] - length];
            for (cnt1 = k_min_Q_B[ij];
                 cnt1 <= k_max_Q_B[ij];
                 cnt1++)
              for (cnt2 = l_min_Q_B[ij][cnt1];
                   cnt2 <= l_max_Q_B[ij][cnt1];
                   cnt2 += 2)
                for (cnt3 = k_min_Q[my_iindx[j + 1] - length];
                     cnt3 <= k_max_Q[my_iindx[j + 1] - length];
                     cnt3++)
                  for (cnt4 = l_min_Q[my_iindx[j + 1] - length][cnt3];
                       cnt4 <= l_max_Q[my_iindx[j + 1] - length][cnt3];
                       cnt4 += 2)
                    if (((cnt1 + cnt3 + da) > maxD1) || ((cnt2 + cnt4 + db) > maxD2)) {
                      qt += qkl * Q_B[ij][cnt1][cnt2 / 2] * Q[my_iindx[j + 1] - length][cnt3][cnt4 / 2];
                      if (qt >= r)
                        goto pbacktrack_ext_loop_early_escape_rem;
                    }
          }
        } /* end if(type) */
      }   /* end for(j) */
      cnt1 = cnt2 = cnt3 = cnt4 = -1;
      /* dont forget the case where i pairs with n */
      j     = length;
      ij    = my_iindx[i] - j;
      type  = ptype[jindx[j] + i];
      if (type) {
        double qkl = exp_E_ExtLoop(type, (i > 1) ? S1[i - 1] : -1, (j < n) ? S1[j + 1] : -1, pf_params);
        if (Q_B_rem[ij] != 0.) {
          qt += qkl * Q_B_rem[ij];
          if (qt >= r)
            goto pbacktrack_ext_loop_early_escape_rem;
        }

        /* if we still search for pairing partner j, we go on here... */
        if (Q_B[ij]) {
          da  = referenceBPs1[sn] - referenceBPs1[ij];
          db  = referenceBPs2[sn] - referenceBPs2[ij];
          for (cnt1 = k_min_Q_B[ij];
               cnt1 <= k_max_Q_B[ij];
               cnt1++)
            for (cnt2 = l_min_Q_B[ij][cnt1];
                 cnt2 <= l_max_Q_B[ij][cnt1];
                 cnt2 += 2)
              if (((cnt1 + da) > maxD1) || ((cnt2 + db) > maxD2)) {
                qt += qkl * Q_B[ij][cnt1][cnt2 / 2];
                if (qt >= r)
                  goto pbacktrack_ext_loop_early_escape_rem;
              }
        }
      } /* end if(type) */

      j++;

pbacktrack_ext_loop_early_escape_rem:

      if (j == length + 1)
        vrna_message_error("pbacktrack@2Dpfold.c: backtracking failed in ext loop (rem)");

      /* finally start backtracking the first exterior stem */
      backtrack(vc, pstruc, cnt1, cnt2, i, j);
      if (j == length)
        break;

      start = j + 1;
      d1    = cnt3;
      d2    = cnt4;
    } /* end if d1 ==-1 */
    else {
      qln_i = Q[sn][d1][d2 / 2];

      /* open chain ? */
      if ((d1 == referenceBPs1[sn])
          && (d2 == referenceBPs2[sn])) {
        r = vrna_urn() * qln_i;
        if (scale[length - start + 1] > r)
          return pstruc;
      }

      for (i = start; i < length; i++) {
        r       = vrna_urn() * qln_i;
        da      = referenceBPs1[sn] - referenceBPs1[my_iindx[i + 1] - length];
        db      = referenceBPs2[sn] - referenceBPs2[my_iindx[i + 1] - length];
        qln_i1  = 0;
        if (d1 >= da && d2 >= db) {
          if (
            (d1 - da >= k_min_Q[my_iindx[i + 1] - length])
            && (d1 - da <= k_max_Q[my_iindx[i + 1] - length])) {
            if (
              (d2 - db >= l_min_Q[my_iindx[i + 1] - length][d1 - da])
              && (d2 - db <= l_max_Q[my_iindx[i + 1] - length][d1 - da]))
              qln_i1 += Q[my_iindx[i + 1] - length][d1 - da][(d2 - db) / 2];
          }
        }

        if (r > qln_i1 * scale[1])
          break;                         /* i is paired */

        qln_i = qln_i1;
      }

      if (i >= length)
        break;              /* no more pairs */

      /* now find the pairing partner j */
      r = vrna_urn() * (qln_i - qln_i1 * scale[1]);

      for (qt = 0, j = i + 1; j < length; j++) {
        int type;
        ij    = my_iindx[i] - j;
        type  = ptype[jindx[j] + i];
        if (type) {
          double qkl = 1.0;
          qkl *= exp_E_ExtLoop(type, (i > 1) ? S1[i - 1] : -1, S1[j + 1], pf_params);

          da  = referenceBPs1[sn] - referenceBPs1[ij] - referenceBPs1[my_iindx[j + 1] - length];
          db  = referenceBPs2[sn] - referenceBPs2[ij] - referenceBPs2[my_iindx[j + 1] - length];

          if ((d1 >= da)
              && (d2 >= db)
              && Q_B[ij]
              && Q[my_iindx[j + 1] - length]) {
            for (cnt1 = k_min_Q_B[ij];
                 cnt1 <= MIN2(k_max_Q_B[ij], d1 - da);
                 cnt1++)
              for (cnt2 = l_min_Q_B[ij][cnt1];
                   cnt2 <= MIN2(l_max_Q_B[ij][cnt1], d2 - db);
                   cnt2 += 2)
                if ((d1 - da - cnt1 >= k_min_Q[my_iindx[j + 1] - length])
                    && (d1 - da - cnt1 <= k_max_Q[my_iindx[j + 1] - length])) {
                  if ((d2 - db - cnt2 >= l_min_Q[my_iindx[j + 1] - length][d1 - da - cnt1])
                      && (d2 - db - cnt2 <= l_max_Q[my_iindx[j + 1] - length][d1 - da - cnt1])) {
                    qt += qkl * Q_B[ij][cnt1][cnt2 / 2] * Q[my_iindx[j + 1] - length][d1 - da - cnt1][(d2 - db - cnt2) / 2];
                    if (qt >= r)
                      goto pbacktrack_ext_loop_early_escape;
                  }
                }
          }
        }
      }
      /* now dont forget the case j==n */
      j   = length;
      ij  = my_iindx[i] - j;
      int type = ptype[jindx[j] + i];
      if (type) {
        double qkl = 1.0;

        qkl *= exp_E_ExtLoop(type, (i > 1) ? S1[i - 1] : -1, (j < n) ? S1[j + 1] : -1, pf_params);

        da  = referenceBPs1[sn] - referenceBPs1[ij];
        db  = referenceBPs2[sn] - referenceBPs2[ij];
        if (d1 >= da && d2 >= db) {
          cnt1  = d1 - da;
          cnt2  = d2 - db;
          if ((cnt1 >= k_min_Q_B[ij]) && (cnt1 <= k_max_Q_B[ij])) {
            if ((cnt2 >= l_min_Q_B[ij][cnt1]) && (cnt2 <= l_max_Q_B[ij][cnt1])) {
              qt += qkl * Q_B[ij][cnt1][cnt2 / 2];
              if (qt >= r)
                goto pbacktrack_ext_loop_early_escape; /* j is paired */
            }
          }
        }
      }

      j++;

pbacktrack_ext_loop_early_escape:

      if (j == length + 1)
        vrna_message_error("pbacktrack@2Dpfold.c: backtracking failed in ext loop");

      backtrack(vc, pstruc, cnt1, cnt2, i, j);

      if (j == length)
        break;

      start = j + 1;
      d1    -= cnt1 + da;
      d2    -= cnt2 + db;
    } /* end if d1!=-1 */
  }
  return pstruc;
}


PRIVATE char *
pbacktrack_circ(vrna_fold_compound_t  *vc,
                int                   d1,
                int                   d2)
{
  char *pstruc;
  unsigned int i, n, maxD1, maxD2,
               *referenceBPs1, *referenceBPs2;
  int *my_iindx,
      k_min_Q_c, k_max_Q_c,
      k_min_Q_cH, k_max_Q_cH,
      k_min_Q_cI, k_max_Q_cI,
      k_min_Q_cM, k_max_Q_cM,
      *l_min_Q_c, *l_max_Q_c,
      *l_min_Q_cH, *l_max_Q_cH,
      *l_min_Q_cI, *l_max_Q_cI,
      *l_min_Q_cM, *l_max_Q_cM;
  FLT_OR_DBL r, *scale, qot,
             **Q_c, **Q_cH, **Q_cI, **Q_cM,
             Q_c_rem, Q_cH_rem, Q_cI_rem, Q_cM_rem;
  vrna_mx_pf_t *matrices;
  vrna_md_t *md;
  vrna_exp_param_t *pf_params;

  pf_params     = vc->exp_params;
  md            = &(pf_params->model_details);
  matrices      = vc->exp_matrices;
  n             = vc->length;
  maxD1         = vc->maxD1;
  maxD2         = vc->maxD2;
  my_iindx      = vc->iindx;
  scale         = matrices->scale;
  referenceBPs1 = vc->referenceBPs1;
  referenceBPs2 = vc->referenceBPs2;

  Q_c       = matrices->Q_c;
  l_min_Q_c = matrices->l_min_Q_c;
  l_max_Q_c = matrices->l_max_Q_c;
  k_min_Q_c = matrices->k_min_Q_c;
  k_max_Q_c = matrices->k_max_Q_c;

  Q_cH        = matrices->Q_cH;
  l_min_Q_cH  = matrices->l_min_Q_cH;
  l_max_Q_cH  = matrices->l_max_Q_cH;
  k_min_Q_cH  = matrices->k_min_Q_cH;
  k_max_Q_cH  = matrices->k_max_Q_cH;

  Q_cI        = matrices->Q_cI;
  l_min_Q_cI  = matrices->l_min_Q_cI;
  l_max_Q_cI  = matrices->l_max_Q_cI;
  k_min_Q_cI  = matrices->k_min_Q_cI;
  k_max_Q_cI  = matrices->k_max_Q_cI;

  Q_cM        = matrices->Q_cM;
  l_min_Q_cM  = matrices->l_min_Q_cM;
  l_max_Q_cM  = matrices->l_max_Q_cM;
  k_min_Q_cM  = matrices->k_min_Q_cM;
  k_max_Q_cM  = matrices->k_max_Q_cM;

  Q_c_rem   = matrices->Q_c_rem;
  Q_cH_rem  = matrices->Q_cH_rem;
  Q_cI_rem  = matrices->Q_cI_rem;
  Q_cM_rem  = matrices->Q_cM_rem;

  /* check whether the chosen neighborhood exists at all */
  int dumb = 1;
  if ((d1 == -1) && (Q_c_rem != 0.)) {
    dumb = 0;
  } else {
    if ((k_min_Q_c <= d1) && (k_max_Q_c >= d1)) {
      int l_min = l_min_Q_c[d1];
      if ((d2 % 2) == (l_min % 2))
        if ((l_min <= d2) && (l_max_Q_c[d1] >= d2))
          dumb = 0;
    }
  }

  if (dumb) {
    vrna_message_error("neighborhood %d:%d is not in scope of calculated partition function!\n"
                       "pbacktrack_circ@2Dpfold.c: exiting cheerless...",
                       d1, d2);
  }

  pstruc = vrna_alloc((n + 1) * sizeof(char));

  for (i = 0; i < n; i++)
    pstruc[i] = '.';
  pstruc[i] = '\0';

  /* now we come to the actual backtracking process */

  qot = 0.;
  /* backtrack in rest-partition */
  if (d1 == -1) {
    r = vrna_urn() * Q_c_rem;
    /* open chain ? */
    if ((referenceBPs1[my_iindx[1] - n] > maxD1) || (referenceBPs2[my_iindx[1] - n] > maxD2)) {
      qot = 1.0 * scale[n];
      if (qot >= r)
        goto pbacktrack_circ_escape;
    }

    qot += Q_cH_rem;
    if (qot >= r) {
      backtrack_qcH(vc, pstruc, d1, d2);
      goto pbacktrack_circ_escape;
    }

    qot += Q_cI_rem;
    if (qot >= r) {
      backtrack_qcI(vc, pstruc, d1, d2);
      goto pbacktrack_circ_escape;
    }

    qot += Q_cM_rem;
    if (qot >= r) {
      backtrack_qcM(vc, pstruc, d1, d2);
      goto pbacktrack_circ_escape;
    }

    vrna_message_error("pbacktrack_circ@2Dpfold.c: backtracking failed in exterior loop! Exiting cheerless...");
  }
  /* normal backtracking */
  else {
    r = vrna_urn() * Q_c[d1][d2 / 2];

    /* open chain ? */
    if ((referenceBPs1[my_iindx[1] - n] == d1) && (referenceBPs2[my_iindx[1] - n] == d2)) {
      qot += 1.0 * scale[n];
      if (qot >= r)
        goto pbacktrack_circ_escape;
    }

    /* exterior hairpin loop ? */
    if ((k_min_Q_cH <= d1) && (k_max_Q_cH >= d1)) {
      int l_min = l_min_Q_cH[d1];
      if ((d2 % 2) == (l_min % 2)) {
        if ((l_min <= d2) && (l_max_Q_cH[d1] >= d2)) {
          qot += Q_cH[d1][d2 / 2];
          if (qot >= r) {
            backtrack_qcH(vc, pstruc, d1, d2);
            goto pbacktrack_circ_escape;
          }
        }
      }
    }

    /* exterior interior loop ? */
    if ((k_min_Q_cI <= d1) && (k_max_Q_cI >= d1)) {
      int l_min = l_min_Q_cI[d1];
      if ((d2 % 2) == (l_min % 2)) {
        if ((l_min <= d2) && (l_max_Q_cI[d1] >= d2)) {
          qot += Q_cI[d1][d2 / 2];
          if (qot >= r) {
            backtrack_qcI(vc, pstruc, d1, d2);
            goto pbacktrack_circ_escape;
          }
        }
      }
    }

    /* exterior multibranch loop ? */
    if ((k_min_Q_cM <= d1) && (k_max_Q_cM >= d1)) {
      int l_min = l_min_Q_cM[d1];
      if ((d2 % 2) == (l_min % 2)) {
        if ((l_min <= d2) && (l_max_Q_cM[d1] >= d2)) {
          qot += Q_cM[d1][d2 / 2];
          if (qot >= r) {
            backtrack_qcM(vc, pstruc, d1, d2);
            goto pbacktrack_circ_escape;
          }
        }
      }
    }
  }

pbacktrack_circ_escape:
  return pstruc;
}


PRIVATE void
backtrack_qcH(vrna_fold_compound_t  *vc,
              char                  *pstruc,
              int                   d1,
              int                   d2)
{
  char *ptype, *sequence;
  short *S1;
  unsigned int i, j, n, maxD1, maxD2,
               base_d1, base_d2, da, db,
               *referenceBPs1, *referenceBPs2;
  int u, *my_iindx, *jindx, ij, cnt1, cnt2, type,
      **l_min_Q_B, **l_max_Q_B,
      *k_min_Q_B, *k_max_Q_B, *rtype;
  FLT_OR_DBL r, qt, *scale, qot,
             ***Q_B, **Q_cH, *Q_B_rem,
             Q_cH_rem;

  vrna_exp_param_t *pf_params;
  vrna_md_t *md;
  vrna_mx_pf_t *matrices;

  pf_params     = vc->exp_params;
  md            = &(pf_params->model_details);
  matrices      = vc->exp_matrices;
  sequence      = vc->sequence;
  n             = vc->length;
  my_iindx      = vc->iindx;
  jindx         = vc->jindx;
  scale         = matrices->scale;
  ptype         = vc->ptype;
  rtype         = &(md->rtype[0]);
  S1            = vc->sequence_encoding;
  referenceBPs1 = vc->referenceBPs1;
  referenceBPs2 = vc->referenceBPs2;
  maxD1         = vc->maxD1;
  maxD2         = vc->maxD2;

  Q_B_rem   = matrices->Q_B_rem;
  Q_B       = matrices->Q_B;
  l_min_Q_B = matrices->l_min_Q_B;
  l_max_Q_B = matrices->l_max_Q_B;
  k_min_Q_B = matrices->k_min_Q_B;
  k_max_Q_B = matrices->k_max_Q_B;

  Q_cH_rem  = matrices->Q_cH_rem;
  Q_cH      = matrices->Q_cH;

  qot = qt = 0.;

  base_d1 = referenceBPs1[my_iindx[1] - n];
  base_d2 = referenceBPs2[my_iindx[1] - n];

  if (d1 == -1) {
    r = vrna_urn() * Q_cH_rem;
    for (i = 1; i < n; i++)
      for (j = i + TURN + 1; j <= n; j++) {
        char loopseq[10];
        ij  = my_iindx[i] - j;
        u   = n - j + i - 1;
        if (u < TURN)
          continue;

        type = ptype[jindx[j] + i];
        if (!type)
          continue;

        if (((type == 3) || (type == 4)) && no_closingGU)
          continue;

        type = rtype[type];
        if (u < 7) {
          strcpy(loopseq, sequence + j - 1);
          strncat(loopseq, sequence, i);
        }

        qt = exp_E_Hairpin(u, type,
                           S1[j + 1], S1[i - 1],
                           loopseq, pf_params)
             * scale[u];

        if (Q_B_rem[ij]) {
          qot += Q_B_rem[ij] * qt;
          if (qot >= r) {
            backtrack(vc, pstruc, d1, d2, i, j);
            return;
          }
        }

        da  = base_d1 - referenceBPs1[ij];
        db  = base_d2 - referenceBPs2[ij];

        if (Q_B[ij]) {
          for (cnt1 = k_min_Q_B[ij];
               cnt1 <= k_max_Q_B[ij];
               cnt1++)
            for (cnt2 = l_min_Q_B[ij][cnt1];
                 cnt2 <= l_max_Q_B[ij][cnt1];
                 cnt2 += 2) {
              if (((cnt1 + da) > maxD1)
                  || ((cnt2 + db) > maxD2)) {
                qot += Q_B[ij][cnt1][cnt2 / 2] * qt;
                if (qot >= r) {
                  backtrack(vc, pstruc, cnt1, cnt2, i, j);
                  return;
                }
              }
            }
        }
      }
  } else {
    r = vrna_urn() * Q_cH[d1][d2 / 2];
    for (i = 1; i < n; i++)
      for (j = i + TURN + 1; j <= n; j++) {
        char loopseq[10];
        ij = my_iindx[i] - j;
        if (!Q_B[ij])
          continue;

        u = n - j + i - 1;
        if (u < TURN)
          continue;

        type = ptype[jindx[j] + i];
        if (!type)
          continue;

        if (((type == 3) || (type == 4)) && no_closingGU)
          continue;

        type = rtype[type];
        if (u < 7) {
          strcpy(loopseq, sequence + j - 1);
          strncat(loopseq, sequence, i);
        }

        qt = exp_E_Hairpin(u, type,
                           S1[j + 1], S1[i - 1],
                           loopseq, pf_params)
             * scale[u];
        da  = base_d1 - referenceBPs1[ij];
        db  = base_d2 - referenceBPs2[ij];

        for (cnt1 = k_min_Q_B[ij];
             cnt1 <= k_max_Q_B[ij];
             cnt1++)
          for (cnt2 = l_min_Q_B[ij][cnt1];
               cnt2 <= l_max_Q_B[ij][cnt1];
               cnt2 += 2) {
            if (((cnt1 + da) == d1)
                && ((cnt2 + db) == d2)) {
              qot += Q_B[ij][cnt1][cnt2 / 2] * qt;
              if (qot >= r) {
                backtrack(vc, pstruc, cnt1, cnt2, i, j);
                return;
              }
            }
          }
      }
  }

  vrna_message_error("backtrack_qcH@2Dpfold.c: failed to find closing pair!");
}


PRIVATE void
backtrack_qcI(vrna_fold_compound_t  *vc,
              char                  *pstruc,
              int                   d1,
              int                   d2)
{
  char *ptype;
  short *S1;
  unsigned int i, j, ij, p, q, pq, n, maxD1, maxD2,
               base_d1, base_d2, da, db,
               *referenceBPs1, *referenceBPs2;
  int *my_iindx, *jindx, cnt1, cnt2, cnt3, cnt4, type,
      **l_min_Q_B, **l_max_Q_B,
      *k_min_Q_B, *k_max_Q_B, *rtype;
  FLT_OR_DBL r, qt, *scale, qot,
             ***Q_B, *Q_B_rem,
             **Q_cI, Q_cI_rem;
  vrna_exp_param_t *pf_params;
  vrna_md_t *md;
  vrna_mx_pf_t *matrices;

  pf_params     = vc->exp_params;
  md            = &(pf_params->model_details);
  matrices      = vc->exp_matrices;
  n             = vc->length;
  my_iindx      = vc->iindx;
  jindx         = vc->jindx;
  scale         = matrices->scale;
  ptype         = vc->ptype;
  rtype         = &(md->rtype[0]);
  S1            = vc->sequence_encoding;
  referenceBPs1 = vc->referenceBPs1;
  referenceBPs2 = vc->referenceBPs2;
  maxD1         = vc->maxD1;
  maxD2         = vc->maxD2;

  Q_B       = matrices->Q_B;
  l_min_Q_B = matrices->l_min_Q_B;
  l_max_Q_B = matrices->l_max_Q_B;
  k_min_Q_B = matrices->k_min_Q_B;
  k_max_Q_B = matrices->k_max_Q_B;

  Q_cI      = matrices->Q_cI;
  Q_B_rem   = matrices->Q_B_rem;
  Q_cI_rem  = matrices->Q_cI_rem;

  qot = qt = 0.;

  base_d1 = referenceBPs1[my_iindx[1] - n];
  base_d2 = referenceBPs2[my_iindx[1] - n];

  if (d1 == -1) {
    r = vrna_urn() * Q_cI_rem;
    for (i = 1; i < n; i++)
      for (j = i + TURN + 1; j <= n; j++) {
        ij    = my_iindx[i] - j;
        type  = rtype[(unsigned int)ptype[jindx[j] + i]];
        if (!type)
          continue;

        if (Q_B_rem[ij]) {
          for (p = j + 1; p < n; p++) {
            unsigned int ln1, qstart, ln_pre;
            ln1 = p - j - 1;
            if (ln1 + i - 1 > MAXLOOP)
              break;

            qstart  = p + TURN + 1;
            ln_pre  = ln1 + i + n;
            if (ln_pre > qstart + MAXLOOP)
              qstart = ln_pre - MAXLOOP - 1;

            for (q = qstart; q <= n; q++) {
              unsigned int ln2;
              int type2;
              pq  = my_iindx[p] - q;
              ln2 = (i - 1) + (n - q);
              if ((ln1 + ln2) > MAXLOOP)
                continue;

              type2 = ptype[jindx[q] + p];
              if (!type2)
                continue;

              qt = exp_E_IntLoop(ln2, ln1,
                                 rtype[type2], type,
                                 S1[q + 1], S1[p - 1],
                                 S1[i - 1], S1[j + 1],
                                 pf_params)
                   * scale[ln1 + ln2];
              if (Q_B_rem[pq]) {
                qot += Q_B_rem[ij] * Q_B_rem[pq] * qt;
                if (qot > r) {
                  backtrack(vc, pstruc, d1, d2, i, j);
                  backtrack(vc, pstruc, d1, d2, p, q);
                  return;
                }
              }

              if (Q_B[pq]) {
                for (cnt1 = k_min_Q_B[pq];
                     cnt1 <= k_max_Q_B[pq];
                     cnt1++)
                  for (cnt2 = l_min_Q_B[pq][cnt1];
                       cnt2 <= l_max_Q_B[pq][cnt1];
                       cnt2 += 2) {
                    qot += Q_B_rem[ij] * Q_B[pq][cnt1][cnt2 / 2] * qt;
                    if (qot > r) {
                      backtrack(vc, pstruc, d1, d2, i, j);
                      backtrack(vc, pstruc, cnt1, cnt2, p, q);
                      return;
                    }
                  }
              }
            }
          }
        }

        if (Q_B[ij]) {
          for (p = j + 1; p < n; p++) {
            unsigned int ln1, qstart, ln_pre;
            ln1 = p - j - 1;
            if (ln1 + i - 1 > MAXLOOP)
              break;

            qstart  = p + TURN + 1;
            ln_pre  = ln1 + i + n;
            if (ln_pre > qstart + MAXLOOP)
              qstart = ln_pre - MAXLOOP - 1;

            for (q = qstart; q <= n; q++) {
              unsigned int ln2;
              int type2;
              pq  = my_iindx[p] - q;
              ln2 = (i - 1) + (n - q);
              if ((ln1 + ln2) > MAXLOOP)
                continue;

              type2 = ptype[jindx[q] + p];
              if (!type2)
                continue;

              qt = exp_E_IntLoop(ln2, ln1,
                                 rtype[type2], type,
                                 S1[q + 1], S1[p - 1],
                                 S1[i - 1], S1[j + 1],
                                 pf_params)
                   * scale[ln1 + ln2];
              if (Q_B_rem[pq]) {
                for (cnt1 = k_min_Q_B[ij];
                     cnt1 <= k_max_Q_B[ij];
                     cnt1++)
                  for (cnt2 = l_min_Q_B[ij][cnt1];
                       cnt2 <= l_max_Q_B[ij][cnt1];
                       cnt2 += 2) {
                    qot += Q_B[ij][cnt1][cnt2 / 2] * Q_B_rem[pq] * qt;
                    if (qot > r) {
                      backtrack(vc, pstruc, cnt1, cnt2, i, j);
                      backtrack(vc, pstruc, d1, d2, p, q);
                      return;
                    }
                  }
              }

              if (Q_B[pq]) {
                da = base_d1
                     - referenceBPs1[ij]
                     - referenceBPs1[pq];
                db = base_d2
                     - referenceBPs2[ij]
                     - referenceBPs2[pq];
                for (cnt1 = k_min_Q_B[ij];
                     cnt1 <= k_max_Q_B[ij];
                     cnt1++)
                  for (cnt2 = l_min_Q_B[ij][cnt1];
                       cnt2 <= l_max_Q_B[ij][cnt1];
                       cnt2 += 2)
                    for (cnt3 = k_min_Q_B[pq];
                         cnt3 <= k_max_Q_B[pq];
                         cnt3++)
                      for (cnt4 = l_min_Q_B[pq][cnt3];
                           cnt4 <= l_max_Q_B[pq][cnt3];
                           cnt4 += 2) {
                        if (((cnt1 + cnt3 + da) > maxD1)
                            || ((cnt2 + cnt4 + db) > maxD2)) {
                          qot += Q_B[ij][cnt1][cnt2 / 2]
                                 * Q_B[pq][cnt3][cnt4 / 2]
                                 * qt;
                          if (qot > r) {
                            backtrack(vc, pstruc, cnt1, cnt2, i, j);
                            backtrack(vc, pstruc, cnt3, cnt4, p, q);
                            return;
                          }
                        }
                      }
              }
            }
          }
        }
      }
  } else {
    r = vrna_urn() * Q_cI[d1][d2 / 2];
    for (i = 1; i < n; i++)
      for (j = i + TURN + 1; j <= n; j++) {
        ij    = my_iindx[i] - j;
        type  = rtype[(unsigned int)ptype[jindx[j] + i]];
        if (!type)
          continue;

        if (!Q_B[ij])
          continue;

        for (p = j + 1; p < n; p++) {
          unsigned int ln1, qstart, ln_pre;
          ln1 = p - j - 1;
          if (ln1 + i - 1 > MAXLOOP)
            break;

          qstart  = p + TURN + 1;
          ln_pre  = ln1 + i + n;
          if (ln_pre > qstart + MAXLOOP)
            qstart = ln_pre - MAXLOOP - 1;

          for (q = qstart; q <= n; q++) {
            unsigned int ln2;
            int type2;
            pq = my_iindx[p] - q;
            if (!Q_B[pq])
              continue;

            ln2 = (i - 1) + (n - q);
            if ((ln1 + ln2) > MAXLOOP)
              continue;

            type2 = ptype[jindx[q] + p];
            if (!type2)
              continue;

            qt = exp_E_IntLoop(ln2, ln1,
                               rtype[type2], type,
                               S1[q + 1], S1[p - 1],
                               S1[i - 1], S1[j + 1],
                               pf_params)
                 * scale[ln1 + ln2];
            da = base_d1
                 - referenceBPs1[ij]
                 - referenceBPs1[pq];
            db = base_d2
                 - referenceBPs2[ij]
                 - referenceBPs2[pq];
            for (cnt1 = k_min_Q_B[ij];
                 cnt1 <= k_max_Q_B[ij];
                 cnt1++)
              for (cnt2 = l_min_Q_B[ij][cnt1];
                   cnt2 <= l_max_Q_B[ij][cnt1];
                   cnt2 += 2)
                for (cnt3 = k_min_Q_B[pq];
                     cnt3 <= k_max_Q_B[pq];
                     cnt3++)
                  for (cnt4 = l_min_Q_B[pq][cnt3];
                       cnt4 <= l_max_Q_B[pq][cnt3];
                       cnt4 += 2) {
                    if (((cnt1 + cnt3 + da) == d1)
                        && ((cnt2 + cnt4 + db) == d2)) {
                      qot += Q_B[ij][cnt1][cnt2 / 2]
                             * Q_B[pq][cnt3][cnt4 / 2]
                             * qt;
                      if (qot > r) {
                        backtrack(vc, pstruc, cnt1, cnt2, i, j);
                        backtrack(vc, pstruc, cnt3, cnt4, p, q);
                        return;
                      }
                    }
                  }
          }
        }
      }
  }
}


PRIVATE void
backtrack_qcM(vrna_fold_compound_t  *vc,
              char                  *pstruc,
              int                   d1,
              int                   d2)
{
  unsigned int k, n, maxD1, maxD2, base_d1, base_d2,
               da, db, *referenceBPs1, *referenceBPs2;
  int *my_iindx, cnt1, cnt2, cnt3, cnt4,
      **l_min_Q_M, **l_max_Q_M,
      **l_min_Q_M2, **l_max_Q_M2,
      *k_min_Q_M, *k_max_Q_M,
      *k_min_Q_M2, *k_max_Q_M2;
  FLT_OR_DBL r, qt, qot,
             ***Q_M, ***Q_M2, **Q_cM,
             *Q_M_rem, *Q_M2_rem, Q_cM_rem;
  vrna_exp_param_t *pf_params;
  vrna_md_t *md;
  vrna_mx_pf_t *matrices;

  pf_params     = vc->exp_params;
  md            = &(pf_params->model_details);
  matrices      = vc->exp_matrices;
  n             = vc->length;
  my_iindx      = vc->iindx;
  referenceBPs1 = vc->referenceBPs1;
  referenceBPs2 = vc->referenceBPs2;
  maxD1         = vc->maxD1;
  maxD2         = vc->maxD2;

  Q_cM = matrices->Q_cM;

  Q_M       = matrices->Q_M;
  l_min_Q_M = matrices->l_min_Q_M;
  l_max_Q_M = matrices->l_max_Q_M;
  k_min_Q_M = matrices->k_min_Q_M;
  k_max_Q_M = matrices->k_max_Q_M;

  Q_M2        = matrices->Q_M2;
  l_min_Q_M2  = matrices->l_min_Q_M2;
  l_max_Q_M2  = matrices->l_max_Q_M2;
  k_min_Q_M2  = matrices->k_min_Q_M2;
  k_max_Q_M2  = matrices->k_max_Q_M2;

  Q_cM_rem  = matrices->Q_cM_rem;
  Q_M_rem   = matrices->Q_M_rem;
  Q_M2_rem  = matrices->Q_M2_rem;

  base_d1 = referenceBPs1[my_iindx[1] - n];
  base_d2 = referenceBPs2[my_iindx[1] - n];
  qot     = qt = 0.;

  if (d1 == -1) {
    r = vrna_urn() * Q_cM_rem;
    for (k = TURN + 2;
         k < n - 2 * TURN - 3;
         k++) {
      if (Q_M_rem[my_iindx[1] - k]) {
        if (Q_M2[k + 1]) {
          for (cnt1 = k_min_Q_M2[k + 1];
               cnt1 <= k_max_Q_M2[k + 1];
               cnt1++)
            for (cnt2 = l_min_Q_M2[k + 1][cnt1];
                 cnt2 <= l_max_Q_M2[k + 1][cnt1];
                 cnt2 += 2) {
              qot += Q_M_rem[my_iindx[1] - k]
                     * Q_M2[k + 1][cnt1][cnt2 / 2]
                     * pf_params->expMLclosing;
              if (qot > r) {
                backtrack_qm(vc, pstruc, d1, d2, 1, k);
                backtrack_qm2(vc, pstruc, cnt1, cnt2, k + 1);
                return;
              }
            }
        }

        if (Q_M2_rem[k + 1]) {
          qot += Q_M_rem[my_iindx[1] - k]
                 * Q_M2_rem[k + 1]
                 * pf_params->expMLclosing;
          if (qot > r) {
            backtrack_qm(vc, pstruc, d1, d2, 1, k);
            backtrack_qm2(vc, pstruc, d1, d2, k + 1);
            return;
          }
        }
      }

      if (Q_M2_rem[k + 1]) {
        if (Q_M[my_iindx[1] - k]) {
          for (cnt1 = k_min_Q_M[my_iindx[1] - k];
               cnt1 <= k_max_Q_M[my_iindx[1] - k];
               cnt1++)
            for (cnt2 = l_min_Q_M[my_iindx[1] - k][cnt1];
                 cnt2 <= l_max_Q_M[my_iindx[1] - k][cnt1];
                 cnt2 += 2) {
              qot += Q_M[my_iindx[1] - k][cnt1][cnt2 / 2]
                     * Q_M2_rem[k + 1]
                     * pf_params->expMLclosing;
              if (qot > r) {
                backtrack_qm(vc, pstruc, cnt1, cnt2, 1, k);
                backtrack_qm2(vc, pstruc, d1, d2, k + 1);
                return;
              }
            }
        }
      }

      da = base_d1
           - referenceBPs1[my_iindx[1] - k]
           - referenceBPs1[my_iindx[k + 1] - n];
      db = base_d2
           - referenceBPs2[my_iindx[1] - k]
           - referenceBPs2[my_iindx[k + 1] - n];

      if (Q_M[my_iindx[1] - k]
          && Q_M2[k + 1]) {
        for (cnt1 = k_min_Q_M[my_iindx[1] - k];
             cnt1 <= k_max_Q_M[my_iindx[1] - k];
             cnt1++)
          for (cnt2 = l_min_Q_M[my_iindx[1] - k][cnt1];
               cnt2 <= l_max_Q_M[my_iindx[1] - k][cnt1];
               cnt2 += 2)
            for (cnt3 = k_min_Q_M2[k + 1];
                 cnt3 <= k_max_Q_M2[k + 1];
                 cnt3++)
              for (cnt4 = l_min_Q_M2[k + 1][cnt3];
                   cnt4 <= l_max_Q_M2[k + 1][cnt3];
                   cnt4 += 2) {
                if (((cnt1 + cnt3 + da) > maxD1)
                    || ((cnt2 + cnt4 + db) > maxD2)) {
                  qot += Q_M[my_iindx[1] - k][cnt1][cnt2 / 2]
                         * Q_M2[k + 1][cnt3][cnt4 / 2]
                         * pf_params->expMLclosing;
                  if (qot > r) {
                    backtrack_qm(vc, pstruc, cnt1, cnt2, 1, k);
                    backtrack_qm2(vc, pstruc, cnt3, cnt4, k + 1);
                    return;
                  }
                }
              }
      }
    }
  } else {
    r = vrna_urn() * Q_cM[d1][d2 / 2];
    for (k = TURN + 2;
         k < n - 2 * TURN - 3;
         k++) {
      da = base_d1
           - referenceBPs1[my_iindx[1] - k]
           - referenceBPs1[my_iindx[k + 1] - n];
      db = base_d2
           - referenceBPs2[my_iindx[1] - k]
           - referenceBPs2[my_iindx[k + 1] - n];
      if (Q_M[my_iindx[1] - k]
          && Q_M2[k + 1]) {
        for (cnt1 = k_min_Q_M[my_iindx[1] - k];
             cnt1 <= k_max_Q_M[my_iindx[1] - k];
             cnt1++)
          for (cnt2 = l_min_Q_M[my_iindx[1] - k][cnt1];
               cnt2 <= l_max_Q_M[my_iindx[1] - k][cnt1];
               cnt2 += 2)
            for (cnt3 = k_min_Q_M2[k + 1];
                 cnt3 <= k_max_Q_M2[k + 1];
                 cnt3++)
              for (cnt4 = l_min_Q_M2[k + 1][cnt3];
                   cnt4 <= l_max_Q_M2[k + 1][cnt3];
                   cnt4 += 2)
                if (((cnt1 + cnt3 + da) == d1)
                    && ((cnt2 + cnt4 + db) == d2)) {
                  qot += Q_M[my_iindx[1] - k][cnt1][cnt2 / 2]
                         * Q_M2[k + 1][cnt3][cnt4 / 2]
                         * pf_params->expMLclosing;
                  if (qot > r) {
                    backtrack_qm(vc, pstruc, cnt1, cnt2, 1, k);
                    backtrack_qm2(vc, pstruc, cnt3, cnt4, k + 1);
                    return;
                  }
                }
      }
    }
  }

  vrna_message_error("backtrack_qcM@2Dpfold.c: backtracking failed");
}


PRIVATE void
backtrack_qm2(vrna_fold_compound_t  *vc,
              char                  *pstruc,
              int                   d1,
              int                   d2,
              unsigned int          k)
{
  unsigned int l, n, maxD1, maxD2, da, db,
               *referenceBPs1, *referenceBPs2;
  int *my_iindx, *jindx, cnt1, cnt2, cnt3, cnt4,
      *k_min_Q_M1, *k_max_Q_M1,
      **l_min_Q_M1, **l_max_Q_M1;
  FLT_OR_DBL r, qt, qot,
             ***Q_M2, ***Q_M1,
             *Q_M2_rem, *Q_M1_rem;

  vrna_exp_param_t *pf_params;      /* holds all [unscaled] pf parameters */
  vrna_md_t *md;
  vrna_mx_pf_t *matrices;

  pf_params = vc->exp_params;
  md        = &(pf_params->model_details);
  matrices  = vc->exp_matrices;

  n             = vc->length;
  my_iindx      = vc->iindx;
  jindx         = vc->jindx;
  referenceBPs1 = vc->referenceBPs1;
  referenceBPs2 = vc->referenceBPs2;
  maxD1         = vc->maxD1;
  maxD2         = vc->maxD2;

  Q_M1_rem    = matrices->Q_M1_rem;
  Q_M1        = matrices->Q_M1;
  l_min_Q_M1  = matrices->l_min_Q_M1;
  l_max_Q_M1  = matrices->l_max_Q_M1;
  k_min_Q_M1  = matrices->k_min_Q_M1;
  k_max_Q_M1  = matrices->k_max_Q_M1;

  Q_M2_rem  = matrices->Q_M2_rem;
  Q_M2      = matrices->Q_M2;

  qot = qt = 0.;

  if (d1 == -1) {
    r = vrna_urn() * Q_M2_rem[k];
    for (l = k + TURN + 1; l < n - TURN - 1; l++) {
      if (Q_M1_rem[jindx[l] + k]) {
        if (Q_M1[jindx[n] + l + 1]) {
          for (cnt1 = k_min_Q_M1[jindx[n] + l + 1];
               cnt1 <= k_max_Q_M1[jindx[n] + l + 1];
               cnt1++)
            for (cnt2 = l_min_Q_M1[jindx[n] + l + 1][cnt1];
                 cnt2 <= l_max_Q_M1[jindx[n] + l + 1][cnt1];
                 cnt2 += 2) {
              qot += Q_M1_rem[jindx[l] + k] * Q_M1[jindx[n] + l + 1][cnt1][cnt2 / 2];
              if (qot > r) {
                backtrack_qm1(vc, pstruc, d1, d2, k, l);
                backtrack_qm1(vc, pstruc, cnt1, cnt2, l + 1, n);
                return;
              }
            }
        }

        if (Q_M1_rem[jindx[n] + l + 1]) {
          qot += Q_M1_rem[jindx[l] + k]
                 * Q_M1_rem[jindx[n] + l + 1];
          if (qot > r) {
            backtrack_qm1(vc, pstruc, d1, d2, k, l);
            backtrack_qm1(vc, pstruc, d1, d2, l + 1, n);
            return;
          }
        }
      }

      if (Q_M1_rem[jindx[n] + l + 1]) {
        if (Q_M1[jindx[l] + k]) {
          for (cnt1 = k_min_Q_M1[jindx[l] + k];
               cnt1 <= k_max_Q_M1[jindx[l] + k];
               cnt1++)
            for (cnt2 = l_min_Q_M1[jindx[l] + k][cnt1];
                 cnt2 <= l_max_Q_M1[jindx[l] + k][cnt1];
                 cnt2 += 2) {
              qot += Q_M1[jindx[l] + k][cnt1][cnt2 / 2]
                     * Q_M1_rem[jindx[n] + l + 1];
              if (qot > r) {
                backtrack_qm1(vc, pstruc, cnt1, cnt2, k, l);
                backtrack_qm1(vc, pstruc, d1, d2, l + 1, n);
                return;
              }
            }
        }
      }

      if (!Q_M1[jindx[l] + k])
        continue;

      if (!Q_M1[jindx[n] + l + 1])
        continue;

      da = referenceBPs1[my_iindx[k] - n]
           - referenceBPs1[my_iindx[k] - l]
           - referenceBPs1[my_iindx[l + 1] - n];
      db = referenceBPs2[my_iindx[k] - n]
           - referenceBPs2[my_iindx[k] - l]
           - referenceBPs2[my_iindx[l + 1] - n];
      for (cnt1 = k_min_Q_M1[jindx[l] + k];
           cnt1 <= k_max_Q_M1[jindx[l] + k];
           cnt1++)
        for (cnt2 = l_min_Q_M1[jindx[l] + k][cnt1];
             cnt2 <= l_max_Q_M1[jindx[l] + k][cnt1];
             cnt2 += 2) {
          for (cnt3 = k_min_Q_M1[jindx[n] + l + 1];
               cnt3 <= k_max_Q_M1[jindx[n] + l + 1];
               cnt3++)
            for (cnt4 = l_min_Q_M1[jindx[n] + l + 1][cnt3];
                 cnt4 <= l_max_Q_M1[jindx[n] + l + 1][cnt3];
                 cnt4 += 2) {
              if (((cnt1 + cnt3 + da) > maxD1)
                  || ((cnt2 + cnt4 + db) > maxD2)) {
                qot += Q_M1[jindx[l] + k][cnt1][cnt2 / 2]
                       * Q_M1[jindx[n] + l + 1][cnt3][cnt4 / 2];
                if (qot > r) {
                  backtrack_qm1(vc, pstruc, cnt1, cnt2, k, l);
                  backtrack_qm1(vc, pstruc, cnt3, cnt4, l + 1, n);
                  return;
                }
              }
            }
        }
    }
  } else {
    r = vrna_urn() * Q_M2[k][d1][d2 / 2];
    for (l = k + TURN + 1; l < n - TURN - 1; l++) {
      if (!Q_M1[jindx[l] + k])
        continue;

      if (!Q_M1[jindx[n] + l + 1])
        continue;

      da = referenceBPs1[my_iindx[k] - n]
           - referenceBPs1[my_iindx[k] - l]
           - referenceBPs1[my_iindx[l + 1] - n];
      db = referenceBPs2[my_iindx[k] - n]
           - referenceBPs2[my_iindx[k] - l]
           - referenceBPs2[my_iindx[l + 1] - n];
      for (cnt1 = k_min_Q_M1[jindx[l] + k];
           cnt1 <= k_max_Q_M1[jindx[l] + k];
           cnt1++)
        for (cnt2 = l_min_Q_M1[jindx[l] + k][cnt1];
             cnt2 <= l_max_Q_M1[jindx[l] + k][cnt1];
             cnt2 += 2) {
          for (cnt3 = k_min_Q_M1[jindx[n] + l + 1];
               cnt3 <= k_max_Q_M1[jindx[n] + l + 1];
               cnt3++)
            for (cnt4 = l_min_Q_M1[jindx[n] + l + 1][cnt3];
                 cnt4 <= l_max_Q_M1[jindx[n] + l + 1][cnt3];
                 cnt4 += 2) {
              if (((cnt1 + cnt3 + da) == d1)
                  && ((cnt2 + cnt4 + db) == d2)) {
                qot += Q_M1[jindx[l] + k][cnt1][cnt2 / 2]
                       * Q_M1[jindx[n] + l + 1][cnt3][cnt4 / 2];
                if (qot > r) {
                  backtrack_qm1(vc, pstruc, cnt1, cnt2, k, l);
                  backtrack_qm1(vc, pstruc, cnt3, cnt4, l + 1, n);
                  return;
                }
              }
            }
        }
    }
  }

  vrna_message_error("backtrack_qm2@2Dpfold.c: backtracking failed");
}


PRIVATE void
backtrack(vrna_fold_compound_t  *vc,
          char                  *pstruc,
          int                   d1,
          int                   d2,
          unsigned int          i,
          unsigned int          j)
{
  FLT_OR_DBL *scale;
  unsigned int maxD1, maxD2, base_d1, base_d2, da, db;
  unsigned int *referenceBPs1, *referenceBPs2;
  char *ptype, *sequence;
  short *S1, *reference_pt1, *reference_pt2;
  int *my_iindx, *jindx, ij, cnt1, cnt2, cnt3, cnt4, *rtype;
  vrna_exp_param_t *pf_params;      /* holds all [unscaled] pf parameters */
  vrna_md_t *md;
  vrna_mx_pf_t *matrices;

  pf_params     = vc->exp_params;
  md            = &(pf_params->model_details);
  matrices      = vc->exp_matrices;
  sequence      = vc->sequence;
  maxD1         = vc->maxD1;
  maxD2         = vc->maxD2;
  my_iindx      = vc->iindx;
  jindx         = vc->jindx;
  scale         = matrices->scale;
  ptype         = vc->ptype;
  rtype         = &(md->rtype[0]);
  S1            = vc->sequence_encoding;
  reference_pt1 = vc->reference_pt1;
  reference_pt2 = vc->reference_pt2;
  referenceBPs1 = vc->referenceBPs1;
  referenceBPs2 = vc->referenceBPs2;

  FLT_OR_DBL ***Q_B, ***Q_M, ***Q_M1, *Q_B_rem, *Q_M_rem, *Q_M1_rem;
  int *k_min_Q_M, *k_max_Q_M, *k_min_Q_M1, *k_max_Q_M1, *k_min_Q_B, *k_max_Q_B;
  int **l_min_Q_M, **l_max_Q_M, **l_min_Q_M1, **l_max_Q_M1, **l_min_Q_B, **l_max_Q_B;

  Q_B       = matrices->Q_B;
  k_min_Q_B = matrices->k_min_Q_B;
  k_max_Q_B = matrices->k_max_Q_B;
  l_min_Q_B = matrices->l_min_Q_B;
  l_max_Q_B = matrices->l_max_Q_B;

  Q_M       = matrices->Q_M;
  k_min_Q_M = matrices->k_min_Q_M;
  k_max_Q_M = matrices->k_max_Q_M;
  l_min_Q_M = matrices->l_min_Q_M;
  l_max_Q_M = matrices->l_max_Q_M;

  Q_M1        = matrices->Q_M1;
  k_min_Q_M1  = matrices->k_min_Q_M1;
  k_max_Q_M1  = matrices->k_max_Q_M1;
  l_min_Q_M1  = matrices->l_min_Q_M1;
  l_max_Q_M1  = matrices->l_max_Q_M1;

  Q_B_rem   = matrices->Q_B_rem;
  Q_M_rem   = matrices->Q_M_rem;
  Q_M1_rem  = matrices->Q_M1_rem;

  do {
    double r, qbt1 = 0.;
    unsigned int k, l, u, u1;
    int type;

    pstruc[i - 1] = '(';
    pstruc[j - 1] = ')';

    r   = 0.;
    ij  = my_iindx[i] - j;

    if (d1 == -1) {
      r = vrna_urn() * Q_B_rem[ij];
      if (r == 0.)
        vrna_message_error("backtrack@2Dpfold.c: backtracking failed\n");

      type    = ptype[jindx[j] + i];
      u       = j - i - 1;
      base_d1 = ((unsigned int)reference_pt1[i] != j) ? 1 : -1;
      base_d2 = ((unsigned int)reference_pt2[i] != j) ? 1 : -1;

      da  = base_d1 + referenceBPs1[ij];
      db  = base_d2 + referenceBPs2[ij];

      /* hairpin ? */
      if ((da > maxD1) || (db > maxD2))
        if (!(((type == 3) || (type == 4)) && no_closingGU))
          qbt1 = exp_E_Hairpin(u, type, S1[i + 1], S1[j - 1], sequence + i - 1, pf_params) * scale[u + 2];

      if (qbt1 >= r)
        return;            /* found the hairpin we're done */

      /* lets see if we form an interior loop */
      for (k = i + 1; k <= MIN2(i + MAXLOOP + 1, j - TURN - 2); k++) {
        unsigned int u_pre, lmin;
        u1    = k - i - 1;
        lmin  = k + TURN + 1;
        u_pre = u1 + j;
        /* lmin = MAX2(k + TURN + 1, u1 + j - 1 - MAXLOOP) */
        if (u_pre > lmin + MAXLOOP)
          lmin = u_pre - 1 - MAXLOOP;

        for (l = lmin; l < j; l++) {
          int type_2;
          type_2 = ptype[jindx[l] + k];
          if (type_2) {
            cnt1    = cnt2 = -1;
            da      = base_d1 + referenceBPs1[my_iindx[i] - j] - referenceBPs1[my_iindx[k] - l];
            db      = base_d2 + referenceBPs2[my_iindx[i] - j] - referenceBPs2[my_iindx[k] - l];
            type_2  = rtype[type_2];
            FLT_OR_DBL tmp_en = exp_E_IntLoop(u1, j - l - 1, type, type_2, S1[i + 1], S1[j - 1], S1[k - 1], S1[l + 1], pf_params) * scale[u1 + j - l + 1];

            if (Q_B_rem[my_iindx[k] - l] != 0.) {
              qbt1 += Q_B_rem[my_iindx[k] - l] * tmp_en;
              if (qbt1 > r)
                goto backtrack_int_early_escape_rem;
            }

            if (Q_B[my_iindx[k] - l]) {
              for (cnt1 = k_min_Q_B[my_iindx[k] - l];
                   cnt1 <= k_max_Q_B[my_iindx[k] - l];
                   cnt1++)
                for (cnt2 = l_min_Q_B[my_iindx[k] - l][cnt1];
                     cnt2 <= l_max_Q_B[my_iindx[k] - l][cnt1];
                     cnt2 += 2)
                  if (((cnt1 + da) > maxD1) || ((cnt2 + db) > maxD2)) {
                    qbt1 += Q_B[my_iindx[k] - l][cnt1][cnt2 / 2] * tmp_en;
                    if (qbt1 > r)
                      goto backtrack_int_early_escape_rem;
                  }
            }
          }
        }
      }
backtrack_int_early_escape_rem:
      if (l < j) {
        i   = k;
        j   = l;
        d1  = cnt1;
        d2  = cnt2;
      } else {
        break;
      }
    } else {
      if ((d1 >= k_min_Q_B[ij]) && (d1 <= k_max_Q_B[ij]))
        if ((d2 >= l_min_Q_B[ij][d1]) && (d2 <= l_max_Q_B[ij][d1]))
          r = vrna_urn() * Q_B[ij][d1][d2 / 2];

      if (r == 0.)
        vrna_message_error("backtrack@2Dpfold.c: backtracking failed\n");

      type    = ptype[jindx[j] + i];
      u       = j - i - 1;
      base_d1 = ((unsigned int)reference_pt1[i] != j) ? 1 : -1;
      base_d2 = ((unsigned int)reference_pt2[i] != j) ? 1 : -1;

      da  = base_d1 + referenceBPs1[ij];
      db  = base_d2 + referenceBPs2[ij];

      /*hairpin contribution*/
      if ((da == d1) && (db == d2))
        if (!(((type == 3) || (type == 4)) && no_closingGU))
          qbt1 = exp_E_Hairpin(u, type, S1[i + 1], S1[j - 1], sequence + i - 1, pf_params) * scale[u + 2];

      if (qbt1 >= r)
        return;            /* found the hairpin we're done */

      for (k = i + 1; k <= MIN2(i + MAXLOOP + 1, j - TURN - 2); k++) {
        unsigned int u_pre, lmin;
        u1    = k - i - 1;
        lmin  = k + TURN + 1;
        u_pre = u1 + j;
        /* lmin = MAX2(k + TURN + 1, u1 + j - 1 - MAXLOOP) */
        if (u_pre > lmin + MAXLOOP)
          lmin = u_pre - 1 - MAXLOOP;

        for (l = lmin; l < j; l++) {
          int type_2;
          type_2 = ptype[jindx[l] + k];
          if (type_2) {
            da      = base_d1 + referenceBPs1[my_iindx[i] - j] - referenceBPs1[my_iindx[k] - l];
            db      = base_d2 + referenceBPs2[my_iindx[i] - j] - referenceBPs2[my_iindx[k] - l];
            type_2  = rtype[type_2];
            FLT_OR_DBL tmp_en = exp_E_IntLoop(u1, j - l - 1, type, type_2, S1[i + 1], S1[j - 1], S1[k - 1], S1[l + 1], pf_params) * scale[u1 + j - l + 1];
            if (d1 >= da && d2 >= db) {
              if ((d1 - da >= k_min_Q_B[my_iindx[k] - l]) && (d1 - da <= k_max_Q_B[my_iindx[k] - l])) {
                if ((d2 - db >= l_min_Q_B[my_iindx[k] - l][d1 - da]) && (d2 - db <= l_max_Q_B[my_iindx[k] - l][d1 - da])) {
                  cnt1  = d1 - da;
                  cnt2  = d2 - db;
                  qbt1  += Q_B[my_iindx[k] - l][cnt1][cnt2 / 2] * tmp_en;
                  if (qbt1 > r)
                    goto backtrack_int_early_escape;
                }
              }
            }
          }
        }
      }

backtrack_int_early_escape:
      if (l < j) {
        i   = k;
        j   = l;
        d1  = cnt1;
        d2  = cnt2;
      } else {
        break;
      }
    }
  } while (1);

  /* backtrack in multi-loop */
  {
    double r, qt;
    unsigned int k, ii, jj;

    base_d1 = ((unsigned int)reference_pt1[i] != j) ? 1 : -1;
    base_d2 = ((unsigned int)reference_pt2[i] != j) ? 1 : -1;

    base_d1 += referenceBPs1[my_iindx[i] - j];
    base_d2 += referenceBPs2[my_iindx[i] - j];

    i++;
    j--;
    /* find the first split index */
    ii  = my_iindx[i];  /* ii-j=[i,j] */
    jj  = jindx[j];     /* jj+i=[j,i] */
    if (d1 == -1) {
      /* get total contribution for current part */
      for (qt = 0., k = i + 1; k < j; k++) {
        if (Q_M_rem[ii - k + 1] != 0.) {
          if (Q_M1[jj + k]) {
            for (cnt1 = k_min_Q_M1[jj + k];
                 cnt1 <= k_max_Q_M1[jj + k];
                 cnt1++)
              for (cnt2 = l_min_Q_M1[jj + k][cnt1];
                   cnt2 <= l_max_Q_M1[jj + k][cnt1];
                   cnt2 += 2)
                qt += Q_M_rem[ii - k + 1] * Q_M1[jj + k][cnt1][cnt2 / 2];
          }

          if (Q_M1_rem[jj + k] != 0.)
            qt += Q_M_rem[ii - k + 1] * Q_M1_rem[jj + k];
        }

        if (Q_M1_rem[jj + k] != 0.) {
          if (Q_M[ii - k + 1]) {
            for (cnt1 = k_min_Q_M[ii - k + 1];
                 cnt1 <= k_max_Q_M[ii - k + 1];
                 cnt1++)
              for (cnt2 = l_min_Q_M[ii - k + 1][cnt1];
                   cnt2 <= l_max_Q_M[ii - k + 1][cnt1];
                   cnt2 += 2)
                qt += Q_M[ii - k + 1][cnt1][cnt2 / 2] * Q_M1_rem[jj + k];
          }
        }

        /* calculate introduced distance to reference structures */
        if (!Q_M[ii - k + 1])
          continue;

        if (!Q_M1[jj + k])
          continue;

        da  = base_d1 - referenceBPs1[my_iindx[i] - k + 1] - referenceBPs1[my_iindx[k] - j];
        db  = base_d2 - referenceBPs2[my_iindx[i] - k + 1] - referenceBPs2[my_iindx[k] - j];
        /* collect all contributing energies */
        for (cnt1 = k_min_Q_M[ii - k + 1];
             cnt1 <= k_max_Q_M[ii - k + 1];
             cnt1++)
          for (cnt2 = l_min_Q_M[ii - k + 1][cnt1];
               cnt2 <= l_max_Q_M[ii - k + 1][cnt1];
               cnt2 += 2)
            for (cnt3 = k_min_Q_M1[jj + k];
                 cnt3 <= k_max_Q_M1[jj + k];
                 cnt3++)
              for (cnt4 = l_min_Q_M1[jj + k][cnt3];
                   cnt4 <= l_max_Q_M1[jj + k][cnt3];
                   cnt4 += 2)
                if (((cnt1 + cnt3 + da) > maxD1) || ((cnt2 + cnt4 + db) > maxD2))
                  qt += Q_M[ii - k + 1][cnt1][cnt2 / 2] * Q_M1[jj + k][cnt3][cnt4 / 2];
      }
      /* throw the dice */
      r = vrna_urn() * qt;
      for (qt = 0., k = i + 1; k < j; k++) {
        cnt1 = cnt2 = cnt3 = cnt4 = -1;
        if (Q_M_rem[ii - k + 1] != 0.) {
          if (Q_M1_rem[jj + k] != 0) {
            qt += Q_M_rem[ii - k + 1] * Q_M1_rem[jj + k];
            if (qt >= r)
              goto backtrack_ml_early_escape;
          }

          if (Q_M1[jj + k]) {
            for (cnt3 = k_min_Q_M1[jj + k];
                 cnt3 <= k_max_Q_M1[jj + k];
                 cnt3++)
              for (cnt4 = l_min_Q_M1[jj + k][cnt3];
                   cnt4 <= l_max_Q_M1[jj + k][cnt3];
                   cnt4 += 2) {
                qt += Q_M_rem[ii - k + 1] * Q_M1[jj + k][cnt3][cnt4 / 2];
                if (qt >= r)
                  goto backtrack_ml_early_escape;
              }
          }
        }

        if (Q_M1_rem[jj + k] != 0.) {
          cnt3 = cnt4 = -1;
          if (Q_M[ii - k + 1]) {
            for (cnt1 = k_min_Q_M[ii - k + 1];
                 cnt1 <= k_max_Q_M[ii - k + 1];
                 cnt1++)
              for (cnt2 = l_min_Q_M[ii - k + 1][cnt1];
                   cnt2 <= l_max_Q_M[ii - k + 1][cnt1];
                   cnt2 += 2) {
                qt += Q_M[ii - k + 1][cnt1][cnt2 / 2] * Q_M1_rem[jj + k];
                if (qt >= r)
                  goto backtrack_ml_early_escape;
              }
          }
        }

        /* calculate introduced distance to reference structures */
        da  = base_d1 - referenceBPs1[my_iindx[i] - k + 1] - referenceBPs1[my_iindx[k] - j];
        db  = base_d2 - referenceBPs2[my_iindx[i] - k + 1] - referenceBPs2[my_iindx[k] - j];
        /* collect all contributing energies */
        if (!Q_M[ii - k + 1])
          continue;

        if (!Q_M1[jj + k])
          continue;

        for (cnt1 = k_min_Q_M[ii - k + 1];
             cnt1 <= k_max_Q_M[ii - k + 1];
             cnt1++)
          for (cnt2 = l_min_Q_M[ii - k + 1][cnt1];
               cnt2 <= l_max_Q_M[ii - k + 1][cnt1];
               cnt2 += 2)
            for (cnt3 = k_min_Q_M1[jj + k];
                 cnt3 <= k_max_Q_M1[jj + k];
                 cnt3++)
              for (cnt4 = l_min_Q_M1[jj + k][cnt3];
                   cnt4 <= l_max_Q_M1[jj + k][cnt3];
                   cnt4 += 2)
                if (((cnt1 + cnt3 + da) > maxD1) || ((cnt2 + cnt4 + db) > maxD2)) {
                  qt += Q_M[ii - k + 1][cnt1][cnt2 / 2] * Q_M1[jj + k][cnt3][cnt4 / 2];
                  if (qt >= r)
                    goto backtrack_ml_early_escape;
                }
      }
    } else {
      /* get total contribution */
      for (qt = 0., k = i + 1; k < j; k++) {
        /* calculate introduced distance to reference structures */
        da  = base_d1 - referenceBPs1[my_iindx[i] - k + 1] - referenceBPs1[my_iindx[k] - j];
        db  = base_d2 - referenceBPs2[my_iindx[i] - k + 1] - referenceBPs2[my_iindx[k] - j];
        /* collect all contributing energies */
        if (d1 >= da && d2 >= db && Q_M[ii - k + 1] && Q_M1[jj + k]) {
          for (cnt1 = k_min_Q_M[ii - k + 1]; cnt1 <= MIN2(k_max_Q_M[ii - k + 1], d1 - da); cnt1++)
            for (cnt2 = l_min_Q_M[ii - k + 1][cnt1]; cnt2 <= MIN2(l_max_Q_M[ii - k + 1][cnt1], d2 - db); cnt2 += 2)
              if ((d1 - cnt1 - da >= k_min_Q_M1[jj + k]) && (d1 - cnt1 - da <= k_max_Q_M1[jj + k]))
                if ((d2 - cnt2 - db >= l_min_Q_M1[jj + k][d1 - da - cnt1]) && (d2 - cnt2 - db <= l_max_Q_M1[jj + k][d1 - cnt1 - da]))
                  qt += Q_M[ii - k + 1][cnt1][cnt2 / 2] * Q_M1[jj + k][d1 - da - cnt1][(d2 - db - cnt2) / 2];
        }
      }
      r = vrna_urn() * qt;
      for (qt = 0., k = i + 1; k < j; k++) {
        /* calculate introduced distance to reference structures */
        da  = base_d1 - referenceBPs1[my_iindx[i] - k + 1] - referenceBPs1[my_iindx[k] - j];
        db  = base_d2 - referenceBPs2[my_iindx[i] - k + 1] - referenceBPs2[my_iindx[k] - j];
        /* collect all contributing energies */
        if (d1 >= da && d2 >= db && Q_M[ii - k + 1] && Q_M1[jj + k]) {
          for (cnt1 = k_min_Q_M[ii - k + 1]; cnt1 <= MIN2(k_max_Q_M[ii - k + 1], d1 - da); cnt1++)
            for (cnt2 = l_min_Q_M[ii - k + 1][cnt1]; cnt2 <= MIN2(l_max_Q_M[ii - k + 1][cnt1], d2 - db); cnt2 += 2)
              if ((d1 - cnt1 - da >= k_min_Q_M1[jj + k]) && (d1 - cnt1 - da <= k_max_Q_M1[jj + k])) {
                if ((d2 - cnt2 - db >= l_min_Q_M1[jj + k][d1 - da - cnt1]) && (d2 - cnt2 - db <= l_max_Q_M1[jj + k][d1 - cnt1 - da])) {
                  cnt3  = d1 - da - cnt1;
                  cnt4  = d2 - db - cnt2;
                  qt    += Q_M[ii - k + 1][cnt1][cnt2 / 2] * Q_M1[jj + k][cnt3][cnt4 / 2];
                  if (qt >= r)
                    goto backtrack_ml_early_escape;
                }
              }
        }
      }
    }

    if (k >= j)
      vrna_message_error("backtrack failed, can't find split index ");

backtrack_ml_early_escape:

    backtrack_qm1(vc, pstruc, cnt3, cnt4, k, j);

    j = k - 1;
    backtrack_qm(vc, pstruc, cnt1, cnt2, i, j);
  }
}


PRIVATE void
backtrack_qm1(vrna_fold_compound_t  *vc,
              char                  *pstruc,
              int                   d1,
              int                   d2,
              unsigned int          i,
              unsigned int          j)
{
  /* i is paired to l, i<l<j; backtrack in qm1 to find l */
  FLT_OR_DBL r, qt, *scale;
  unsigned int maxD1, maxD2, da, db;
  unsigned int *referenceBPs1, *referenceBPs2;
  char *ptype;
  short *S1;
  int *my_iindx, *jindx, cnt1, cnt2;

  vrna_exp_param_t *pf_params;      /* holds all [unscaled] pf parameters */
  vrna_md_t *md;
  vrna_mx_pf_t *matrices;

  pf_params     = vc->exp_params;
  md            = &(pf_params->model_details);
  matrices      = vc->exp_matrices;
  maxD1         = vc->maxD1;
  maxD2         = vc->maxD2;
  my_iindx      = vc->iindx;
  jindx         = vc->jindx;
  scale         = matrices->scale;
  ptype         = vc->ptype;
  S1            = vc->sequence_encoding;
  referenceBPs1 = vc->referenceBPs1;
  referenceBPs2 = vc->referenceBPs2;

  FLT_OR_DBL ***Q_B, ***Q_M1, *Q_B_rem, *Q_M1_rem;
  int *k_min_Q_M1, *k_max_Q_M1, *k_min_Q_B, *k_max_Q_B;
  int **l_min_Q_M1, **l_max_Q_M1, **l_min_Q_B, **l_max_Q_B;

  Q_B       = matrices->Q_B;
  k_min_Q_B = matrices->k_min_Q_B;
  k_max_Q_B = matrices->k_max_Q_B;
  l_min_Q_B = matrices->l_min_Q_B;
  l_max_Q_B = matrices->l_max_Q_B;

  Q_M1        = matrices->Q_M1;
  k_min_Q_M1  = matrices->k_min_Q_M1;
  k_max_Q_M1  = matrices->k_max_Q_M1;
  l_min_Q_M1  = matrices->l_min_Q_M1;
  l_max_Q_M1  = matrices->l_max_Q_M1;

  Q_B_rem   = matrices->Q_B_rem;
  Q_M1_rem  = matrices->Q_M1_rem;

  unsigned int ii, l;
  int type;

  r = 0.;

  /* find qm1 contribution */
  if (d1 == -1) {
    r = vrna_urn() * Q_M1_rem[jindx[j] + i];
  } else {
    if ((d1 >= k_min_Q_M1[jindx[j] + i]) && (d1 <= k_max_Q_M1[jindx[j] + i]))
      if ((d2 >= l_min_Q_M1[jindx[j] + i][d1]) && (d2 <= l_max_Q_M1[jindx[j] + i][d1]))
        r = vrna_urn() * Q_M1[jindx[j] + i][d1][d2 / 2];
  }

  if (r == 0.)
    vrna_message_error("backtrack_qm1@2Dpfold.c: backtracking failed\n");

  ii = my_iindx[i];
  for (qt = 0., l = i + TURN + 1; l <= j; l++) {
    type = ptype[jindx[l] + i];
    if (type) {
      FLT_OR_DBL tmp = exp_E_MLstem(type, S1[i - 1], S1[l + 1], pf_params) * pow(pf_params->expMLbase, j - l) * scale[j - l];
      /* compute the introduced distance to reference structures */
      da    = referenceBPs1[my_iindx[i] - j] - referenceBPs1[my_iindx[i] - l];
      db    = referenceBPs2[my_iindx[i] - j] - referenceBPs2[my_iindx[i] - l];
      cnt1  = cnt2 = -1;
      if (d1 == -1) {
        if (Q_B_rem[ii - l] != 0.) {
          qt += Q_B_rem[ii - l] * tmp;
          if (qt >= r)
            goto backtrack_qm1_early_escape;
        }

        if (Q_B[ii - l]) {
          for (cnt1 = k_min_Q_B[ii - l];
               cnt1 <= k_max_Q_B[ii - l];
               cnt1++)
            for (cnt2 = l_min_Q_B[ii - l][cnt1];
                 cnt2 <= l_max_Q_B[ii - l][cnt1];
                 cnt2 += 2)
              if (((cnt1 + da) > maxD1) || ((cnt2 + db) > maxD2)) {
                qt += Q_B[ii - l][cnt1][cnt2 / 2] * tmp;
                if (qt >= r)
                  goto backtrack_qm1_early_escape;
              }
        }
      } else {
        /* get energy contributions */
        if (d1 >= da && d2 >= db) {
          if ((d1 - da >= k_min_Q_B[ii - l]) && (d1 - da <= k_max_Q_B[ii - l])) {
            if ((d2 - db >= l_min_Q_B[ii - l][d1 - da]) && (d2 - db <= l_max_Q_B[ii - l][d1 - da])) {
              cnt1  = d1 - da;
              cnt2  = d2 - db;
              qt    += Q_B[ii - l][cnt1][cnt2 / 2] * tmp;
              if (qt >= r)
                goto backtrack_qm1_early_escape;
            }
          }
        }
      }
    }
  }
  if (l > j)
    vrna_message_error("backtrack failed in qm1");

backtrack_qm1_early_escape:

  backtrack(vc, pstruc, cnt1, cnt2, i, l);
}


PRIVATE void
backtrack_qm(vrna_fold_compound_t *vc,
             char                 *pstruc,
             int                  d1,
             int                  d2,
             unsigned int         i,
             unsigned int         j)
{
  /* divide multiloop into qm and qm1  */
  FLT_OR_DBL r, *scale;
  unsigned int maxD1, maxD2, da, db, da2, db2;
  unsigned int *referenceBPs1, *referenceBPs2;
  int *my_iindx, *jindx, cnt1, cnt2, cnt3, cnt4;

  vrna_exp_param_t *pf_params;      /* holds all [unscaled] pf parameters */
  vrna_md_t *md;
  vrna_mx_pf_t *matrices;

  pf_params     = vc->exp_params;
  md            = &(pf_params->model_details);
  matrices      = vc->exp_matrices;
  maxD1         = vc->maxD1;
  maxD2         = vc->maxD2;
  my_iindx      = vc->iindx;
  jindx         = vc->jindx;
  scale         = matrices->scale;
  referenceBPs1 = vc->referenceBPs1;
  referenceBPs2 = vc->referenceBPs2;

  FLT_OR_DBL ***Q_M, ***Q_M1, *Q_M_rem, *Q_M1_rem;
  int *k_min_Q_M, *k_max_Q_M, *k_min_Q_M1, *k_max_Q_M1;
  int **l_min_Q_M, **l_max_Q_M, **l_min_Q_M1, **l_max_Q_M1;

  Q_M       = matrices->Q_M;
  k_min_Q_M = matrices->k_min_Q_M;
  k_max_Q_M = matrices->k_max_Q_M;
  l_min_Q_M = matrices->l_min_Q_M;
  l_max_Q_M = matrices->l_max_Q_M;

  Q_M1        = matrices->Q_M1;
  k_min_Q_M1  = matrices->k_min_Q_M1;
  k_max_Q_M1  = matrices->k_max_Q_M1;
  l_min_Q_M1  = matrices->l_min_Q_M1;
  l_max_Q_M1  = matrices->l_max_Q_M1;
  Q_M_rem     = matrices->Q_M_rem;
  Q_M1_rem    = matrices->Q_M1_rem;

  double qmt = 0;
  unsigned int k;
  while (j > i) {
    /* now backtrack  [i ... j] in qm[] */

    /* find qm contribution */
    if (d1 == -1) {
      r = vrna_urn() * Q_M_rem[my_iindx[i] - j];
    } else {
      if (Q_M[my_iindx[i] - j])
        if ((d1 >= k_min_Q_M[my_iindx[i] - j]) && (d1 <= k_max_Q_M[my_iindx[i] - j]))
          if ((d2 >= l_min_Q_M[my_iindx[i] - j][d1]) && (d2 <= l_max_Q_M[my_iindx[i] - j][d1]))
            r = vrna_urn() * Q_M[my_iindx[i] - j][d1][d2 / 2];
    }

    if (r == 0.)
      vrna_message_error("backtrack_qm@2Dpfold.c: backtracking failed in finding qm contribution\n");

    qmt = 0.;
    if (d1 == -1) {
      if (Q_M1_rem[jindx[j] + i] != 0.) {
        qmt += Q_M1_rem[jindx[j] + i];
        if (qmt >= r) {
          backtrack_qm1(vc, pstruc, d1, d2, i, j);
          return;
        }
      }

      for (k = i + 1; k <= j; k++) {
        FLT_OR_DBL tmp = pow(pf_params->expMLbase, k - i) * scale[k - i];
        if (Q_M1_rem[jindx[j] + k] != 0.) {
          qmt += Q_M1_rem[jindx[j] + k] * tmp;
          if (qmt >= r) {
            backtrack_qm1(vc, pstruc, d1, d2, k, j);
            return;
          }
        }

        da2 = referenceBPs1[my_iindx[i] - j] - referenceBPs1[my_iindx[k] - j];
        db2 = referenceBPs2[my_iindx[i] - j] - referenceBPs2[my_iindx[k] - j];
        if (Q_M1[jindx[j] + k]) {
          for (cnt1 = k_min_Q_M1[jindx[j] + k];
               cnt1 <= k_max_Q_M1[jindx[j] + k];
               cnt1++)
            for (cnt2 = l_min_Q_M1[jindx[j] + k][cnt1];
                 cnt2 <= l_max_Q_M1[jindx[j] + k][cnt1];
                 cnt2 += 2)
              if (((cnt1 + da2) > maxD1) || ((cnt2 + db2) > maxD2)) {
                qmt += Q_M1[jindx[j] + k][cnt1][cnt2 / 2] * tmp;
                if (qmt >= r) {
                  backtrack_qm1(vc, pstruc, cnt1, cnt2, k, j);
                  return;
                }
              }
        }

        da  = da2 - referenceBPs1[my_iindx[i] - k + 1];
        db  = db2 - referenceBPs2[my_iindx[i] - k + 1];

        cnt1 = cnt2 = cnt3 = cnt4 = -1;
        if (Q_M_rem[my_iindx[i] - k + 1] != 0.) {
          if (Q_M1_rem[jindx[j] + k] != 0.) {
            qmt += Q_M_rem[my_iindx[i] - k + 1] * Q_M1_rem[jindx[j] + k];
            if (qmt >= r)
              goto backtrack_qm_early_escape;
          }

          if (Q_M1[jindx[j] + k]) {
            for (cnt3 = k_min_Q_M1[jindx[j] + k];
                 cnt3 <= k_max_Q_M1[jindx[j] + k];
                 cnt3++)
              for (cnt4 = l_min_Q_M1[jindx[j] + k][cnt3];
                   cnt4 <= l_max_Q_M1[jindx[j] + k][cnt3];
                   cnt4 += 2) {
                qmt += Q_M_rem[my_iindx[i] - k + 1] * Q_M1[jindx[j] + k][cnt3][cnt4 / 2];
                if (qmt >= r)
                  goto backtrack_qm_early_escape;
              }
          }
        }

        if (Q_M1_rem[jindx[j] + k] != 0.) {
          cnt3 = cnt4 = -1;
          if (Q_M[my_iindx[i] - k + 1]) {
            for (cnt1 = k_min_Q_M[my_iindx[i] - k + 1];
                 cnt1 <= k_max_Q_M[my_iindx[i] - k + 1];
                 cnt1++)
              for (cnt2 = l_min_Q_M[my_iindx[i] - k + 1][cnt1];
                   cnt2 <= l_max_Q_M[my_iindx[i] - k + 1][cnt1];
                   cnt2 += 2) {
                qmt += Q_M[my_iindx[i] - k + 1][cnt1][cnt2 / 2] * Q_M1_rem[jindx[j] + k];
                if (qmt >= r)
                  goto backtrack_qm_early_escape;
              }
          }
        }

        if (!Q_M[my_iindx[i] - k + 1])
          continue;

        if (!Q_M1[jindx[j] + k])
          continue;

        for (cnt1 = k_min_Q_M[my_iindx[i] - k + 1];
             cnt1 <= k_max_Q_M[my_iindx[i] - k + 1];
             cnt1++)
          for (cnt2 = l_min_Q_M[my_iindx[i] - k + 1][cnt1];
               cnt2 <= l_max_Q_M[my_iindx[i] - k + 1][cnt1];
               cnt2 += 2)
            for (cnt3 = k_min_Q_M1[jindx[j] + k];
                 cnt3 <= k_max_Q_M1[jindx[j] + k];
                 cnt3++)
              for (cnt4 = l_min_Q_M1[jindx[j] + k][cnt3];
                   cnt4 <= l_max_Q_M1[jindx[j] + k][cnt3];
                   cnt4 += 2)
                if (((cnt1 + cnt3 + da) > maxD1) || ((cnt2 + cnt4 + db) > maxD2)) {
                  qmt += Q_M[my_iindx[i] - k + 1][cnt1][cnt2 / 2] * Q_M1[jindx[j] + k][cnt3][cnt4 / 2];
                  if (qmt >= r)
                    goto backtrack_qm_early_escape;
                }
      }
    } else {
      /* find corresponding qm1 contribution */
      if (Q_M1[jindx[j] + i]) {
        if ((d1 >= k_min_Q_M1[jindx[j] + i]) && (d1 <= k_max_Q_M1[jindx[j] + i]))
          if ((d2 >= l_min_Q_M1[jindx[j] + i][d1]) && (d2 <= l_max_Q_M1[jindx[j] + i][d1]))
            qmt = Q_M1[jindx[j] + i][d1][d2 / 2];
      }

      k = i;
      if (qmt < r) {
        for (k = i + 1; k <= j; k++) {
          /* calculate introduced distancies to reference structures */
          da2 = referenceBPs1[my_iindx[i] - j] - referenceBPs1[my_iindx[k] - j];
          db2 = referenceBPs2[my_iindx[i] - j] - referenceBPs2[my_iindx[k] - j];
          da  = da2 - referenceBPs1[my_iindx[i] - k + 1];
          db  = db2 - referenceBPs2[my_iindx[i] - k + 1];


          FLT_OR_DBL tmp = pow(pf_params->expMLbase, k - i) * scale[k - i];

          /* collect unpaired + qm1 contributions */
          if (d1 >= da2 && d2 >= db2) {
            if ((d1 - da2 >= k_min_Q_M1[jindx[j] + k]) && (d1 - da2 <= k_max_Q_M1[jindx[j] + k])) {
              if ((d2 - db2 >= l_min_Q_M1[jindx[j] + k][d1 - da2]) && (d2 - db2 <= l_max_Q_M1[jindx[j] + k][d1 - da2])) {
                cnt3  = d1 - da2;
                cnt4  = d2 - db2;
                qmt   += Q_M1[jindx[j] + k][cnt3][cnt4 / 2] * tmp;
                if (qmt >= r) {
                  backtrack_qm1(vc, pstruc, cnt3, cnt4, k, j);
                  return;
                }
              }
            }
          }

          /* collect qm + qm1 contributions */
          if (d1 >= da && d2 >= db && Q_M[my_iindx[i] - k + 1] && Q_M1[jindx[j] + k]) {
            for (cnt1 = k_min_Q_M[my_iindx[i] - k + 1]; cnt1 <= MIN2(k_max_Q_M[my_iindx[i] - k + 1], d1 - da); cnt1++)
              for (cnt2 = l_min_Q_M[my_iindx[i] - k + 1][cnt1]; cnt2 <= MIN2(l_max_Q_M[my_iindx[i] - k + 1][cnt1], d2 - db); cnt2 += 2)
                if ((d1 - da - cnt1 >= k_min_Q_M1[jindx[j] + k]) && (d1 - da - cnt1 <= k_max_Q_M1[jindx[j] + k])) {
                  if ((d2 - db - cnt2 >= l_min_Q_M1[jindx[j] + k][d1 - da - cnt1]) && (d2 - db - cnt2 <= l_max_Q_M1[jindx[j] + k][d1 - da - cnt1])) {
                    cnt3  = d1 - da - cnt1;
                    cnt4  = d2 - db - cnt2;
                    qmt   += Q_M[my_iindx[i] - k + 1][cnt1][cnt2 / 2] * Q_M1[jindx[j] + k][cnt3][cnt4 / 2];
                    if (qmt >= r)
                      goto backtrack_qm_early_escape;
                  }
                }
          }
        }
      } else {
        backtrack_qm1(vc, pstruc, d1, d2, k, j);
        return;
      }
    }

    if (k > j)
      vrna_message_error("backtrack_qm@2Dpfold.c: backtrack failed in qm");

backtrack_qm_early_escape:

    backtrack_qm1(vc, pstruc, cnt3, cnt4, k, j);

    if (k < i + TURN)
      break;            /* no more pairs */

    d1  = cnt1;
    d2  = cnt2;


    if (d1 == referenceBPs1[my_iindx[i] - k + 1] && d2 == referenceBPs2[my_iindx[i] - k + 1]) {
      /* is interval [i,k] totally unpaired? */
      FLT_OR_DBL tmp = pow(pf_params->expMLbase, k - i) * scale[k - i];
      r = vrna_urn() * (Q_M[my_iindx[i] - k + 1][d1][d2 / 2] + tmp);
      if (tmp >= r)
        return;            /* no more pairs */
    }

    j = k - 1;
  }
}


PRIVATE void
adjustArrayBoundaries(FLT_OR_DBL  ***array,
                      int         *k_min,
                      int         *k_max,
                      int         **l_min,
                      int         **l_max,
                      int         k_min_post,
                      int         k_max_post,
                      int         *l_min_post,
                      int         *l_max_post)
{
  int cnt1;
  int k_diff_pre  = k_min_post - *k_min;
  int mem_size    = k_max_post - k_min_post + 1;

  if (k_min_post < INF) {
    /* free all the unused memory behind actual data */
    for (cnt1 = k_max_post + 1; cnt1 <= *k_max; cnt1++) {
      (*array)[cnt1] += (*l_min)[cnt1] / 2;
      free((*array)[cnt1]);
    }

    /* free unused memory before actual data */
    for (cnt1 = *k_min; cnt1 < k_min_post; cnt1++) {
      (*array)[cnt1] += (*l_min)[cnt1] / 2;
      free((*array)[cnt1]);
    }
    /* move data to front and thereby eliminating unused memory in front of actual data */
    if (k_diff_pre > 0) {
      memmove((FLT_OR_DBL **)(*array), ((FLT_OR_DBL **)(*array)) + k_diff_pre, sizeof(FLT_OR_DBL *) * mem_size);
      memmove((int *)(*l_min), ((int *)(*l_min)) + k_diff_pre, sizeof(int) * mem_size);
      memmove((int *)(*l_max), ((int *)(*l_max)) + k_diff_pre, sizeof(int) * mem_size);
    }

    /* reallocating memory to actual size used */
    *array  += *k_min;
    *array  = (FLT_OR_DBL **)realloc(*array, sizeof(FLT_OR_DBL *) * mem_size);
    *array  -= k_min_post;

    *l_min  += *k_min;
    *l_min  = (int *)realloc(*l_min, sizeof(int) * mem_size);
    *l_min  -= k_min_post;

    *l_max  += *k_min;
    *l_max  = (int *)realloc(*l_max, sizeof(int) * mem_size);
    *l_max  -= k_min_post;


    for (cnt1 = k_min_post; cnt1 <= k_max_post; cnt1++) {
      if (l_min_post[cnt1] < INF) {
        /* new memsize */
        mem_size = (l_max_post[cnt1] - l_min_post[cnt1] + 1) / 2 + 1;
        /* reshift the pointer */
        (*array)[cnt1] += (*l_min)[cnt1] / 2;

        int shift = (l_min_post[cnt1] % 2 == (*l_min)[cnt1] % 2) ? 0 : 1;
        /* eliminate unused memory in front of actual data */
        unsigned int start = (l_min_post[cnt1] - (*l_min)[cnt1]) / 2 + shift;
        if (start > 0)
          memmove((FLT_OR_DBL *)((*array)[cnt1]), (FLT_OR_DBL *)((*array)[cnt1]) + start, sizeof(FLT_OR_DBL) * mem_size);

        (*array)[cnt1] = (FLT_OR_DBL *)realloc((*array)[cnt1], sizeof(FLT_OR_DBL) * mem_size);

        (*array)[cnt1] -= l_min_post[cnt1] / 2;
      } else {
        /* free according memory */
        (*array)[cnt1] += (*l_min)[cnt1] / 2;
        free((*array)[cnt1]);
      }

      (*l_min)[cnt1]  = l_min_post[cnt1];
      (*l_max)[cnt1]  = l_max_post[cnt1];
    }
  } else {
    /* we have to free all unused memory */
    for (cnt1 = *k_min; cnt1 <= *k_max; cnt1++) {
      (*array)[cnt1] += (*l_min)[cnt1] / 2;
      free((*array)[cnt1]);
    }
    (*l_min)  += *k_min;
    (*l_max)  += *k_min;
    free(*l_min);
    free(*l_max);
    (*array) += *k_min;
    free(*array);
    *array = NULL;
  }

  l_min_post  += *k_min;
  l_max_post  += *k_min;
  *k_min      = k_min_post;
  *k_max      = k_max_post;

  free(l_min_post);
  free(l_max_post);
}


PRIVATE INLINE void
preparePosteriorBoundaries(int  size,
                           int  shift,
                           int  *min_k,
                           int  *max_k,
                           int  **min_l,
                           int  **max_l)
{
  int i;

  *min_k  = INF;
  *max_k  = 0;

  *min_l  = (int *)vrna_alloc(sizeof(int) * size);
  *max_l  = (int *)vrna_alloc(sizeof(int) * size);

  for (i = 0; i < size; i++) {
    (*min_l)[i] = INF;
    (*max_l)[i] = 0;
  }

  *min_l  -= shift;
  *max_l  -= shift;
}


PRIVATE INLINE void
updatePosteriorBoundaries(int d1,
                          int d2,
                          int *min_k,
                          int *max_k,
                          int **min_l,
                          int **max_l)
{
  (*min_l)[d1]  = MIN2((*min_l)[d1], d2);
  (*max_l)[d1]  = MAX2((*max_l)[d1], d2);
  *min_k        = MIN2(*min_k, d1);
  *max_k        = MAX2(*max_k, d1);
}


PRIVATE INLINE void
prepareBoundaries(int min_k_pre,
                  int max_k_pre,
                  int min_l_pre,
                  int max_l_pre,
                  int bpdist,
                  int *min_k,
                  int *max_k,
                  int **min_l,
                  int **max_l)
{
  int cnt;
  int mem = max_k_pre - min_k_pre + 1;

  *min_k  = min_k_pre;
  *max_k  = max_k_pre;
  *min_l  = (int *)vrna_alloc(sizeof(int) * mem);
  *max_l  = (int *)vrna_alloc(sizeof(int) * mem);

  *min_l  -= min_k_pre;
  *max_l  -= min_k_pre;

  /* for each k guess the according minimum l*/
  for (cnt = min_k_pre; cnt <= max_k_pre; cnt++) {
    (*min_l)[cnt] = min_l_pre;
    (*max_l)[cnt] = max_l_pre;
    while ((*min_l)[cnt] + cnt < bpdist)
      (*min_l)[cnt]++;
    if ((bpdist % 2) != (((*min_l)[cnt] + cnt) % 2))
      (*min_l)[cnt]++;
  }
}


PRIVATE INLINE void
prepareArray(FLT_OR_DBL ***array,
             int        min_k,
             int        max_k,
             int        *min_l,
             int        *max_l)
{
  int i, mem;

  *array  = (FLT_OR_DBL **)vrna_alloc(sizeof(FLT_OR_DBL *) * (max_k - min_k + 1));
  *array  -= min_k;

  for (i = min_k; i <= max_k; i++) {
    mem         = (max_l[i] - min_l[i] + 1) / 2 + 1;
    (*array)[i] = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * mem);
    (*array)[i] -= min_l[i] / 2;
  }
}


/*
 #################################
 # DEPRECATED FUNCTIONS BELOW    #
 #################################
 */
PRIVATE void
crosslink(TwoDpfold_vars *vars)
{
  vrna_fold_compound_t *c;
  vrna_mx_pf_t *m;

  c = vars->compatibility;
  m = c->exp_matrices;

  vars->sequence      = c->sequence;
  vars->seq_length    = c->length;
  vars->reference_pt1 = c->reference_pt1;
  vars->reference_pt2 = c->reference_pt2;
  vars->referenceBPs1 = c->referenceBPs1;
  vars->referenceBPs2 = c->referenceBPs2;
  vars->mm1           = c->mm1;
  vars->mm2           = c->mm2;
  vars->bpdist        = c->bpdist;
  vars->dangles       = c->exp_params->model_details.dangles;
  vars->circ          = c->exp_params->model_details.circ;
  vars->temperature   = c->exp_params->model_details.temperature;
  vars->init_temp     = c->exp_params->model_details.temperature;
  vars->pf_scale      = c->exp_params->pf_scale;
  vars->pf_params     = c->exp_params;

  vars->scale = m->scale;
  vars->ptype = c->ptype_pf_compat;
  vars->S     = c->sequence_encoding2;
  vars->S1    = c->sequence_encoding;

  vars->jindx     = c->jindx;
  vars->my_iindx  = c->iindx;
  vars->maxD1     = c->maxD1;
  vars->maxD2     = c->maxD2;

  vars->Q             = m->Q;
  vars->l_min_values  = m->l_min_Q;
  vars->l_max_values  = m->l_max_Q;
  vars->k_min_values  = m->k_min_Q;
  vars->k_max_values  = m->k_max_Q;

  vars->Q_B             = m->Q_B;
  vars->l_min_values_b  = m->l_min_Q_B;
  vars->l_max_values_b  = m->l_max_Q_B;
  vars->k_min_values_b  = m->k_min_Q_B;
  vars->k_max_values_b  = m->k_max_Q_B;

  vars->Q_M             = m->Q_M;
  vars->l_min_values_m  = m->l_min_Q_M;
  vars->l_max_values_m  = m->l_max_Q_M;
  vars->k_min_values_m  = m->k_min_Q_M;
  vars->k_max_values_m  = m->k_max_Q_M;

  vars->Q_M1            = m->Q_M1;
  vars->l_min_values_m1 = m->l_min_Q_M1;
  vars->l_max_values_m1 = m->l_max_Q_M1;
  vars->k_min_values_m1 = m->k_min_Q_M1;
  vars->k_max_values_m1 = m->k_max_Q_M1;

  vars->Q_M2_rem        = m->Q_M2_rem;
  vars->Q_M2            = m->Q_M2;
  vars->l_min_values_m2 = m->l_min_Q_M2;
  vars->l_max_values_m2 = m->l_max_Q_M2;
  vars->k_min_values_m2 = m->k_min_Q_M2;
  vars->k_max_values_m2 = m->k_max_Q_M2;

  vars->Q_c       = m->Q_c;
  vars->Q_cH      = m->Q_cH;
  vars->Q_cI      = m->Q_cI;
  vars->Q_cM      = m->Q_cM;
  vars->Q_c_rem   = m->Q_c_rem;
  vars->Q_cH_rem  = m->Q_cH_rem;
  vars->Q_cI_rem  = m->Q_cI_rem;
  vars->Q_cM_rem  = m->Q_cM_rem;

  vars->Q_rem     = m->Q_rem;
  vars->Q_B_rem   = m->Q_B_rem;
  vars->Q_M_rem   = m->Q_M_rem;
  vars->Q_M1_rem  = m->Q_M1_rem;
}


PUBLIC char *
TwoDpfold_pbacktrack(TwoDpfold_vars *vars,
                     int            d1,
                     int            d2)
{
  return vrna_pbacktrack_TwoD(vars->compatibility, d1, d2);
}


PUBLIC char *
TwoDpfold_pbacktrack5(TwoDpfold_vars  *vars,
                      int             d1,
                      int             d2,
                      unsigned int    length)
{
  return vrna_pbacktrack5_TwoD(vars->compatibility, d1, d2, length);
}


PUBLIC TwoDpfold_vars *
get_TwoDpfold_variables(const char  *seq,
                        const char  *structure1,
                        char        *structure2,
                        int         circ)
{
  vrna_md_t md;
  TwoDpfold_vars *vars;
  vrna_fold_compound_t *c;
  vrna_mx_mfe_t *m;

  set_model_details(&md);
  md.circ = circ;

  vars                = (TwoDpfold_vars *)malloc(sizeof(TwoDpfold_vars));
  vars->compatibility = vrna_fold_compound_TwoD(seq, structure1, structure2, &md, VRNA_OPTION_PF);

  crosslink(vars);

  return vars;
}


PUBLIC void
destroy_TwoDpfold_variables(TwoDpfold_vars *vars)
{
  if (vars == NULL)
    return;

  vrna_fold_compound_free(vars->compatibility);

  free(vars);
}


vrna_sol_TwoD_pf_t *
TwoDpfoldList(TwoDpfold_vars  *vars,
              int             distance1,
              int             distance2)
{
  vrna_sol_TwoD_pf_t *sol;

  sol = vrna_pf_TwoD(vars->compatibility, distance1, distance2);

  crosslink(vars);

  return sol;
}
