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
#include "ViennaRNA/utils.h"
#include "ViennaRNA/energy_par.h"
#include "ViennaRNA/fold_vars.h"
#include "ViennaRNA/fold.h"
#include "ViennaRNA/loop_energies.h"
#include "ViennaRNA/params.h"
#ifdef _OPENMP
#include <omp.h>
#endif
#include "ViennaRNA/2Dfold.h"

/*
 #################################
 # GLOBAL VARIABLES              #
 #################################
 */
int compute_2Dfold_F3 = 0;

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
PRIVATE void  mfe_linear(vrna_fold_compound_t *vc);


PRIVATE void  mfe_circ(vrna_fold_compound_t *vc);


PUBLIC void  update_TwoDfold_params(TwoDfold_vars *vars);


PRIVATE void  backtrack_f5(unsigned int         j,
                           int                  k,
                           int                  l,
                           char                 *structure,
                           vrna_fold_compound_t *vc);


PRIVATE void  backtrack_c(unsigned int          i,
                          unsigned int          j,
                          int                   k,
                          int                   l,
                          char                  *structure,
                          vrna_fold_compound_t  *vc);


PRIVATE void  backtrack_m(unsigned int          i,
                          unsigned int          j,
                          int                   k,
                          int                   l,
                          char                  *structure,
                          vrna_fold_compound_t  *vc);


PRIVATE void  backtrack_m1(unsigned int         i,
                           unsigned int         j,
                           int                  k,
                           int                  l,
                           char                 *structure,
                           vrna_fold_compound_t *vc);


PRIVATE void  backtrack_fc(int                  k,
                           int                  l,
                           char                 *structure,
                           vrna_fold_compound_t *vc);


PRIVATE void  backtrack_m2(unsigned int         i,
                           int                  k,
                           int                  l,
                           char                 *structure,
                           vrna_fold_compound_t *vc);


PRIVATE void  adjustArrayBoundaries(int ***array,
                                    int *k_min,
                                    int *k_max,
                                    int **l_min,
                                    int **l_max,
                                    int k_min_real,
                                    int k_max_real,
                                    int *l_min_real,
                                    int *l_max_real);


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


INLINE PRIVATE void  prepareArray(int ***array,
                                  int min_k,
                                  int max_k,
                                  int *min_l,
                                  int *max_l);


INLINE PRIVATE void  prepareArray2(unsigned long  ***array,
                                   int            min_k,
                                   int            max_k,
                                   int            *min_l,
                                   int            *max_l);


/*
 #################################
 # BEGIN OF FUNCTION DEFINITIONS #
 #################################
 */

#if 0

PRIVATE void
initialize_TwoDfold_vars(TwoDfold_vars *vars)
{
  update_TwoDfold_params(vars);
  /* this call updates the params in the ViennaRNA fold.o which is a global, so be careful
   *  whith calling it parallel... need a workarround or fix of ViennaRNA fold stuff
   */
  update_fold_params();
}


PUBLIC TwoDfold_solution **
TwoDfold(TwoDfold_vars  *vars,
         int            distance1,
         int            distance2)
{
  unsigned int      i, d1, d2;
  unsigned int      maxD1;
  unsigned int      maxD2;
  unsigned int      length;
  TwoDfold_solution **output;

  initialize_TwoDfold_vars(vars);
  if (fabs(vars->P->temperature - temperature) > 1e-6)
    update_TwoDfold_params(vars);

  vars->S   = encode_sequence(vars->sequence, 0);
  vars->S1  = encode_sequence(vars->sequence, 1);

  make_ptypes(vars);

  maxD1 = vars->maxD1;
  maxD2 = vars->maxD2;

  if (distance1 >= 0) {
    if ((unsigned int)distance1 > maxD1)
      fprintf(stderr,
              "limiting maximum basepair distance 1 to %u\n",
              maxD1);
    else
      maxD1 = (unsigned int)distance1;
  }

  if (distance2 >= 0) {
    if ((unsigned int)distance2 > maxD2)
      fprintf(stderr,
              "limiting maximum basepair distance 2 to %u\n",
              maxD2);
    else
      maxD2 = (unsigned int)distance2;
  }

  vars->maxD1 = maxD1;
  vars->maxD2 = maxD2;
  output      = (TwoDfold_solution **)vrna_alloc((vars->maxD1 + 1) * sizeof(TwoDfold_solution *));

  mfe_linear(vars);
  if (vars->circ)
    mfe_circ(vars);

  length = vars->seq_length;

  for (d1 = 0; d1 <= maxD1; d1++) {
    output[d1] = (TwoDfold_solution *)vrna_alloc((vars->maxD2 + 1) * sizeof(TwoDfold_solution));
#ifdef _OPENMP
#pragma omp parallel for private(d2)
#endif
    for (d2 = 0; d2 <= maxD2; d2++) {
      output[d1][d2].en = (float)INF / (float)100.;
      output[d1][d2].s  = NULL;
    }
    if ((d1 >= ((vars->circ) ? vars->k_min_values_fc : vars->k_min_values_f[length]))
        && (d1 <= ((vars->circ) ? vars->k_max_values_fc : vars->k_max_values_f[length]))) {
#ifdef _OPENMP
#pragma omp parallel for private(d2, i)
#endif
      for (d2 = ((vars->circ) ? vars->l_min_values_fc[d1] : vars->l_min_values_f[length][d1]);
           d2 <= ((vars->circ) ? vars->l_max_values_fc[d1] : vars->l_max_values_f[length][d1]);
           d2 += 2) {
        output[d1][d2].en = (float)((vars->circ) ? vars->E_Fc[d1][d2 / 2] : vars->E_F5[length][d1][d2 / 2]) / (float)100.;
        if (vars->do_backtrack && (output[d1][d2].en != (float)INF / (float)100.)) {
          char *mfe_structure = (char *)vrna_alloc(length + 1);
          for (i = 0; i < length; i++)
            mfe_structure[i] = '.';
          mfe_structure[i] = '\0';
          (vars->circ) ? backtrack_fc(d1, d2, mfe_structure, vars) : backtrack_f5(length, d1, d2, mfe_structure, vars);
          output[d1][d2].s = mfe_structure;
        }
      }
    }
  }
  return output;
}


#endif

PUBLIC vrna_sol_TwoD_t *
vrna_mfe_TwoD(vrna_fold_compound_t  *vars,
              int                   distance1,
              int                   distance2)
{
  unsigned int    i, d1, d2;
  unsigned int    maxD1;
  unsigned int    maxD2;
  unsigned int    length;
  unsigned int    counter = 0;
  int             en      = 0;
  vrna_sol_TwoD_t *output;
  vrna_md_t       *md;
  vrna_mx_mfe_t   *matrices;

  maxD1     = vars->maxD1;
  maxD2     = vars->maxD2;
  matrices  = vars->matrices;
  md        = &(vars->params->model_details);

  if (distance1 >= 0) {
    if ((unsigned int)distance1 > maxD1)
      vrna_message_warning("vrna_mfe_TwoD@2Dfold.c: limiting maximum basepair distance 1 to %u\n",
                           maxD1);
    else
      maxD1 = (unsigned int)distance1;
  }

  if (distance2 >= 0) {
    if ((unsigned int)distance2 > maxD2)
      vrna_message_warning("vrna_mfe_TwoD@2Dfold.c: limiting maximum basepair distance 2 to %u\n",
                           maxD2);
    else
      maxD2 = (unsigned int)distance2;
  }

  vars->maxD1 = maxD1;
  vars->maxD2 = maxD2;
  output      = (vrna_sol_TwoD_t *)vrna_alloc((((vars->maxD1 + 1) * (vars->maxD2 + 2)) / 2 + 2) * sizeof(vrna_sol_TwoD_t));

  mfe_linear(vars);
  if (md->circ)
    mfe_circ(vars);

  length = vars->length;

  for (d1 = 0; d1 <= maxD1; d1++) {
    if ((d1 >= ((md->circ) ? matrices->k_min_Fc : matrices->k_min_F5[length]))
        && (d1 <= ((md->circ) ? matrices->k_max_Fc : matrices->k_max_F5[length]))) {
      for (d2 = ((md->circ) ? matrices->l_min_Fc[d1] : matrices->l_min_F5[length][d1]);
           d2 <= ((md->circ) ? matrices->l_max_Fc[d1] : matrices->l_max_F5[length][d1]);
           d2 += 2) {
        en = ((md->circ) ? matrices->E_Fc[d1][d2 / 2] : matrices->E_F5[length][d1][d2 / 2]);
        if (en == INF)
          continue;

        output[counter].k   = d1;
        output[counter].l   = d2;
        output[counter].en  = (float)en / (float)100.;
        if (md->backtrack) {
          char *mfe_structure = (char *)vrna_alloc(length + 1);
          for (i = 0; i < length; i++)
            mfe_structure[i] = '.';
          mfe_structure[i] = '\0';
          (md->circ) ? backtrack_fc((int)d1, (int)d2, mfe_structure, vars) : backtrack_f5(length, (int)d1, (int)d2, mfe_structure, vars);
          output[counter].s = mfe_structure;
        } else {
          output[counter].s = NULL;
        }

        counter++;
      }
    }
  }

  /* store entry for remaining partition if it exists */
  en = ((md->circ) ? matrices->E_Fc_rem : matrices->E_F5_rem[length]);
  if (en != INF) {
    output[counter].k   = -1;
    output[counter].l   = -1;
    output[counter].en  = (float)en / (float)100.;
    if (md->backtrack) {
      char *mfe_structure = (char *)vrna_alloc(length + 1);
      for (i = 0; i < length; i++)
        mfe_structure[i] = '.';
      mfe_structure[i] = '\0';
      (md->circ) ? backtrack_fc(-1, -1, mfe_structure, vars) : backtrack_f5(length, -1, -1, mfe_structure, vars);
      output[counter].s = mfe_structure;
    } else {
      output[counter].s = NULL;
    }

    counter++;
  }

  /* insert end-marker entry */
  output[counter].k = output[counter].l = INF;
  counter++;

  /* resize to actual dataset amount */
  output = (vrna_sol_TwoD_t *)vrna_realloc(output, sizeof(vrna_sol_TwoD_t) * counter);
  return output;
}


PUBLIC char *
vrna_backtrack5_TwoD(vrna_fold_compound_t *vc,
                     int                  k,
                     int                  l,
                     unsigned int         j)
{
  unsigned int  i;
  char          *mfe_structure = (char *)vrna_alloc(j + 1);

  if (j < TURN + 2)
    return NULL;

  for (i = 0; i < j; i++)
    mfe_structure[i] = '.';
  mfe_structure[i] = '\0';

  backtrack_f5(j, k, l, mfe_structure, vc);
  return mfe_structure;
}


PRIVATE void
mfe_linear(vrna_fold_compound_t *vc)
{
  unsigned int  d, i, j, ij, maxD1, maxD2, seq_length, dia, dib, dja, djb, *referenceBPs1, *referenceBPs2, *mm1, *mm2, *bpdist;
  int           cnt1, cnt2, cnt3, cnt4, d1, d2, energy, dangles, temp2, type, additional_en, *my_iindx, *jindx, circ, *rtype;
  short         *S1, *reference_pt1, *reference_pt2;
  char          *sequence, *ptype;
  vrna_param_t  *P;
  vrna_mx_mfe_t *matrices;
  vrna_md_t     *md;

  /* dereferenciate things we often need */
  P             = vc->params;
  md            = &(P->model_details);
  matrices      = vc->matrices;
  sequence      = vc->sequence;
  seq_length    = vc->length;
  maxD1         = vc->maxD1;
  maxD2         = vc->maxD2;
  S1            = vc->sequence_encoding;
  ptype         = vc->ptype;
  rtype         = &(md->rtype[0]);
  reference_pt1 = vc->reference_pt1;
  reference_pt2 = vc->reference_pt2;
  my_iindx      = vc->iindx;
  jindx         = vc->jindx;
  referenceBPs1 = vc->referenceBPs1;
  referenceBPs2 = vc->referenceBPs2;
  mm1           = vc->mm1;
  mm2           = vc->mm2;
  bpdist        = vc->bpdist;
  dangles       = md->dangles;
  circ          = md->circ;

  for (d = TURN + 2; d <= seq_length; d++) {
    /* i,j in [1..length] */
#ifdef _OPENMP
#pragma omp parallel for private(additional_en, j, energy, temp2, i, ij, dia,dib,dja,djb,cnt1,cnt2,cnt3,cnt4, d1, d2)
#endif
    for (j = d; j <= seq_length; j++) {
      unsigned int  p, q, pq, u, maxp, dij;
      int           type_2, type, tt, no_close, base_d1, base_d2;

      i     = j - d + 1;
      dij   = j - i - 1;
      ij    = my_iindx[i] - j;
      type  = ptype[jindx[j] + i];

      no_close = (((type == 3) || (type == 4)) && no_closingGU);

      if (type) {
        /* we have a pair */
        /* increase or decrease distance-to-reference value depending whether (i,j) is included in
         *  reference or has to be introduced
         */
        base_d1 = ((unsigned int)reference_pt1[i] != j) ? 1 : -1;
        base_d2 = ((unsigned int)reference_pt2[i] != j) ? 1 : -1;

        /* HAIRPIN STRUCTURES */

        /* get distance to reference if closing the hairpin
         *  d = dbp(T_{i,j}, {i,j})
         */
        d1  = base_d1 + referenceBPs1[ij];
        d2  = base_d2 + referenceBPs2[ij];

        int min_k, max_k, min_l, max_l;
        int real_min_k, real_max_k, *min_l_real, *max_l_real;

        min_l = min_k = 0;
        max_k = mm1[ij] + referenceBPs1[ij];
        max_l = mm2[ij] + referenceBPs2[ij];

        prepareBoundaries(min_k,
                          max_k,
                          min_l,
                          max_l,
                          bpdist[ij],
                          &matrices->k_min_C[ij],
                          &matrices->k_max_C[ij],
                          &matrices->l_min_C[ij],
                          &matrices->l_max_C[ij]
                          );

        preparePosteriorBoundaries(matrices->k_max_C[ij] - matrices->k_min_C[ij] + 1,
                                   matrices->k_min_C[ij],
                                   &real_min_k,
                                   &real_max_k,
                                   &min_l_real,
                                   &max_l_real
                                   );

        prepareArray(&matrices->E_C[ij],
                     matrices->k_min_C[ij],
                     matrices->k_max_C[ij],
                     matrices->l_min_C[ij],
                     matrices->l_max_C[ij]
                     );

#ifdef COUNT_STATES
        prepareArray2(&matrices->N_C[ij],
                      matrices->k_min_C[ij],
                      matrices->k_max_C[ij],
                      matrices->l_min_C[ij],
                      matrices->l_max_C[ij]
                      );
#endif

        /* d1 and d2 are the distancies to both references introduced by closing a hairpin structure at (i,j) */
        if ((d1 >= 0) && (d2 >= 0)) {
          if (((unsigned int)d1 <= maxD1) && ((unsigned int)d2 <= maxD2)) {
            matrices->E_C[ij][d1][d2 / 2] = (no_close) ? FORBIDDEN : E_Hairpin(dij, type, S1[i + 1], S1[j - 1], sequence + i - 1, P);
            updatePosteriorBoundaries(d1,
                                      d2,
                                      &real_min_k,
                                      &real_max_k,
                                      &min_l_real,
                                      &max_l_real
                                      );
#ifdef COUNT_STATES
            matrices->N_C[ij][d1][d2 / 2] = 1;
#endif
          } else {
            matrices->E_C_rem[ij] = (no_close) ? FORBIDDEN : E_Hairpin(dij, type, S1[i + 1], S1[j - 1], sequence + i - 1, P);
          }
        }

        /* INTERIOR LOOP STRUCTURES */
        maxp = MIN2(j - 2 - TURN, i + MAXLOOP + 1);
        for (p = i + 1; p <= maxp; p++) {
          unsigned int  minq    = p + TURN + 1;
          unsigned int  ln_pre  = dij + p;
          if (ln_pre > minq + MAXLOOP)
            minq = ln_pre - MAXLOOP - 1;

          for (q = minq; q < j; q++) {
            pq = my_iindx[p] - q;
            /* set distance to reference structure... */
            type_2 = ptype[jindx[q] + p];

            if (type_2 == 0)
              continue;

            type_2 = rtype[type_2];

            /* get distance to reference if closing the interior loop
             *  d2 = dbp(S_{i,j}, S_{p.q} + {i,j})
             */
            d1  = base_d1 + referenceBPs1[ij] - referenceBPs1[pq];
            d2  = base_d2 + referenceBPs2[ij] - referenceBPs2[pq];

            if (no_closingGU)
              if (no_close || (type_2 == 3) || (type_2 == 4))
                if ((p > i + 1) || (q < j - 1))
                  continue;

            /* continue unless stack */

            energy = E_IntLoop(p - i - 1, j - q - 1, type, type_2, S1[i + 1], S1[j - 1], S1[p - 1], S1[q + 1], P);

            if (matrices->E_C[pq] != NULL) {
              for (cnt1 = matrices->k_min_C[pq]; cnt1 <= matrices->k_max_C[pq]; cnt1++) {
                for (cnt2 = matrices->l_min_C[pq][cnt1]; cnt2 <= matrices->l_max_C[pq][cnt1]; cnt2 += 2) {
                  if (matrices->E_C[pq][cnt1][cnt2 / 2] != INF) {
                    if (((cnt1 + d1) <= maxD1) && ((cnt2 + d2) <= maxD2)) {
                      matrices->E_C[ij][cnt1 + d1][(cnt2 + d2) / 2] = MIN2(matrices->E_C[ij][cnt1 + d1][(cnt2 + d2) / 2],
                                                                           matrices->E_C[pq][cnt1][cnt2 / 2] + energy
                                                                           );
                      updatePosteriorBoundaries(cnt1 + d1,
                                                cnt2 + d2,
                                                &real_min_k,
                                                &real_max_k,
                                                &min_l_real,
                                                &max_l_real
                                                );
#ifdef COUNT_STATES
                      matrices->N_C[ij][cnt1 + d1][(cnt2 + d2) / 2] += matrices->N_C[pq][cnt1][cnt2 / 2];
#endif
                    }
                    /* collect all cases where d1+cnt1 or d2+cnt2 exceeds maxD1, maxD2, respectively */
                    else {
                      matrices->E_C_rem[ij] = MIN2(matrices->E_C_rem[ij], matrices->E_C[pq][cnt1][cnt2 / 2] + energy);
                    }
                  }
                }
              }
            }

            /* collect all contributions where C[pq] already lies outside k_max, l_max boundary */
            if (matrices->E_C_rem[pq] != INF)
              matrices->E_C_rem[ij] = MIN2(matrices->E_C_rem[ij], matrices->E_C_rem[pq] + energy);
          } /* end q-loop */
        }   /* end p-loop */


        /* MULTI LOOP STRUCTURES */
        if (!no_close) {
          /* dangle energies for multiloop closing stem */
          tt    = rtype[type];
          temp2 = P->MLclosing;
          if (dangles == 2)
            temp2 += E_MLstem(tt, S1[j - 1], S1[i + 1], P);
          else
            temp2 += E_MLstem(tt, -1, -1, P);

          for (u = i + TURN + 2; u < j - TURN - 2; u++) {
            int i1u   = my_iindx[i + 1] - u;
            int u1j1  = my_iindx[u + 1] - j + 1;
            /* check all cases where either M or M1 are already out of scope of maxD1 and/or maxD2 */
            if (matrices->E_M_rem[i1u] != INF) {
              for (cnt3 = matrices->k_min_M1[u1j1];
                   cnt3 <= matrices->k_max_M1[u1j1];
                   cnt3++)
                for (cnt4 = matrices->l_min_M1[u1j1][cnt3];
                     cnt4 <= matrices->l_max_M1[u1j1][cnt3];
                     cnt4 += 2) {
                  if (matrices->E_M1[u1j1][cnt3][cnt4 / 2] != INF) {
                    matrices->E_C_rem[ij] = MIN2(matrices->E_C_rem[ij],
                                                 matrices->E_M_rem[i1u]
                                                 + matrices->E_M1[u1j1][cnt3][cnt4 / 2]
                                                 + temp2
                                                 );
                  }
                }
              if (matrices->E_M1_rem[u1j1] != INF) {
                matrices->E_C_rem[ij] = MIN2(matrices->E_C_rem[ij],
                                             matrices->E_M_rem[i1u]
                                             + matrices->E_M1_rem[u1j1]
                                             + temp2
                                             );
              }
            }

            if (matrices->E_M1_rem[u1j1] != INF) {
              for (cnt1 = matrices->k_min_M[i1u];
                   cnt1 <= matrices->k_max_M[i1u];
                   cnt1++)
                for (cnt2 = matrices->l_min_M[i1u][cnt1];
                     cnt2 <= matrices->l_max_M[i1u][cnt1];
                     cnt2 += 2)
                  if (matrices->E_M[i1u][cnt1][cnt2 / 2] != INF) {
                    matrices->E_C_rem[ij] = MIN2(matrices->E_C_rem[ij],
                                                 matrices->E_M[i1u][cnt1][cnt2 / 2]
                                                 + matrices->E_M1_rem[u1j1]
                                                 + temp2
                                                 );
                  }
            }

            /* get distance to reference if closing the multiloop
             *  d = dbp(S_{i,j}, {i,j} + S_{i+1,u} + S_{u+1,j-1})
             */
            if (!matrices->E_M[i1u])
              continue;

            if (!matrices->E_M1[u1j1])
              continue;

            d1  = base_d1 + referenceBPs1[ij] - referenceBPs1[i1u] - referenceBPs1[u1j1];
            d2  = base_d2 + referenceBPs2[ij] - referenceBPs2[i1u] - referenceBPs2[u1j1];

            for (cnt1 = matrices->k_min_M[i1u];
                 cnt1 <= matrices->k_max_M[i1u];
                 cnt1++)
              for (cnt2 = matrices->l_min_M[i1u][cnt1];
                   cnt2 <= matrices->l_max_M[i1u][cnt1];
                   cnt2 += 2)
                for (cnt3 = matrices->k_min_M1[u1j1];
                     cnt3 <= matrices->k_max_M1[u1j1];
                     cnt3++)
                  for (cnt4 = matrices->l_min_M1[u1j1][cnt3];
                       cnt4 <= matrices->l_max_M1[u1j1][cnt3];
                       cnt4 += 2) {
                    if ((matrices->E_M[i1u][cnt1][cnt2 / 2] != INF) && (matrices->E_M1[u1j1][cnt3][cnt4 / 2] != INF)) {
                      if (((cnt1 + cnt3 + d1) <= maxD1) && ((cnt2 + cnt4 + d2) <= maxD2)) {
                        matrices->E_C[ij][cnt1 + cnt3 + d1][(cnt2 + cnt4 + d2) / 2] = MIN2(matrices->E_C[ij][cnt1 + cnt3 + d1][(cnt2 + cnt4 + d2) / 2],
                                                                                           matrices->E_M[i1u][cnt1][cnt2 / 2]
                                                                                           + matrices->E_M1[u1j1][cnt3][cnt4 / 2]
                                                                                           + temp2
                                                                                           );
                        updatePosteriorBoundaries(cnt1 + cnt3 + d1,
                                                  cnt2 + cnt4 + d2,
                                                  &real_min_k,
                                                  &real_max_k,
                                                  &min_l_real,
                                                  &max_l_real
                                                  );
#ifdef COUNT_STATES
                        matrices->N_C[ij][cnt1 + cnt3 + d1][(cnt2 + cnt4 + d2) / 2] += matrices->N_M[i1u][cnt1][cnt2 / 2] * matrices->N_M1[u1j1][cnt3][cnt4 / 2];
#endif
                      }
                      /* collect all cases where d1+cnt1+cnt3 or d2+cnt2+cnt4 exceeds maxD1, maxD2, respectively */
                      else {
                        matrices->E_C_rem[ij] = MIN2(matrices->E_C_rem[ij],
                                                     matrices->E_M[i1u][cnt1][cnt2 / 2]
                                                     + matrices->E_M1[u1j1][cnt3][cnt4 / 2]
                                                     + temp2
                                                     );
                      }
                    }
                  }
          }
        }

        /* resize and move memory portions of energy matrix E_C */
        adjustArrayBoundaries(&matrices->E_C[ij],
                              &matrices->k_min_C[ij],
                              &matrices->k_max_C[ij],
                              &matrices->l_min_C[ij],
                              &matrices->l_max_C[ij],
                              real_min_k,
                              real_max_k,
                              min_l_real,
                              max_l_real
                              );
#ifdef COUNT_STATES
        /* actually we should adjust the array boundaries here but we might never use the count states option more than once so what....*/
#endif
      } /* end >> if (pair) << */

      /* done with c[i,j], now compute fML[i,j] */
      /* free ends ? -----------------------------------------*/


      dia = referenceBPs1[ij] - referenceBPs1[my_iindx[i + 1] - j];
      dib = referenceBPs2[ij] - referenceBPs2[my_iindx[i + 1] - j];
      dja = referenceBPs1[ij] - referenceBPs1[ij + 1];
      djb = referenceBPs2[ij] - referenceBPs2[ij + 1];

      if (dangles == 2)
        temp2 = E_MLstem(type, ((i > 1) || circ) ? S1[i - 1] : -1, ((j < seq_length) || circ) ? S1[j + 1] : -1, P);
      else
        temp2 = E_MLstem(type, -1, -1, P);

      int min_k_guess, max_k_guess, min_l_guess, max_l_guess;
      int min_k_real_m, max_k_real_m, *min_l_real_m, *max_l_real_m;
      int min_k_real_m1, max_k_real_m1, *min_l_real_m1, *max_l_real_m1;

      min_k_guess = min_l_guess = 0;
      max_k_guess = mm1[ij] + referenceBPs1[ij];
      max_l_guess = mm2[ij] + referenceBPs2[ij];

      prepareBoundaries(min_k_guess,
                        max_k_guess,
                        min_l_guess,
                        max_l_guess,
                        bpdist[ij],
                        &matrices->k_min_M[ij],
                        &matrices->k_max_M[ij],
                        &matrices->l_min_M[ij],
                        &matrices->l_max_M[ij]
                        );

      prepareBoundaries(min_k_guess,
                        max_k_guess,
                        min_l_guess,
                        max_l_guess,
                        bpdist[ij],
                        &matrices->k_min_M1[ij],
                        &matrices->k_max_M1[ij],
                        &matrices->l_min_M1[ij],
                        &matrices->l_max_M1[ij]
                        );

      preparePosteriorBoundaries(matrices->k_max_M[ij] - matrices->k_min_M[ij] + 1,
                                 matrices->k_min_M[ij],
                                 &min_k_real_m,
                                 &max_k_real_m,
                                 &min_l_real_m,
                                 &max_l_real_m
                                 );
      preparePosteriorBoundaries(matrices->k_max_M1[ij] - matrices->k_min_M1[ij] + 1,
                                 matrices->k_min_M1[ij],
                                 &min_k_real_m1,
                                 &max_k_real_m1,
                                 &min_l_real_m1,
                                 &max_l_real_m1
                                 );

      prepareArray(&matrices->E_M[ij],
                   matrices->k_min_M[ij],
                   matrices->k_max_M[ij],
                   matrices->l_min_M[ij],
                   matrices->l_max_M[ij]
                   );

      prepareArray(&matrices->E_M1[ij],
                   matrices->k_min_M1[ij],
                   matrices->k_max_M1[ij],
                   matrices->l_min_M1[ij],
                   matrices->l_max_M1[ij]
                   );
#ifdef COUNT_STATES
      prepareArray2(&matrices->N_M[ij],
                    matrices->k_min_M[ij],
                    matrices->k_max_M[ij],
                    matrices->l_min_M[ij],
                    matrices->l_max_M[ij]
                    );
      prepareArray2(&matrices->N_M1[ij],
                    matrices->k_min_M1[ij],
                    matrices->k_max_M1[ij],
                    matrices->l_min_M1[ij],
                    matrices->l_max_M1[ij]
                    );
#endif

      /* now to the actual computations... */
      /* 1st E_M[ij] = E_M1[ij] = E_C[ij] + b */
      if (matrices->E_C_rem[ij] != INF)
        matrices->E_M_rem[ij] = matrices->E_M1_rem[ij] = temp2 + matrices->E_C_rem[ij];

      if (matrices->E_C[ij]) {
        for (cnt1 = matrices->k_min_C[ij]; cnt1 <= matrices->k_max_C[ij]; cnt1++) {
          for (cnt2 = matrices->l_min_C[ij][cnt1]; cnt2 <= matrices->l_max_C[ij][cnt1]; cnt2 += 2) {
            if (matrices->E_C[ij][cnt1][cnt2 / 2] != INF) {
              matrices->E_M[ij][cnt1][cnt2 / 2] = matrices->E_M1[ij][cnt1][cnt2 / 2] = temp2 + matrices->E_C[ij][cnt1][cnt2 / 2];
              updatePosteriorBoundaries(cnt1,
                                        cnt2,
                                        &min_k_real_m,
                                        &max_k_real_m,
                                        &min_l_real_m,
                                        &max_l_real_m
                                        );
              updatePosteriorBoundaries(cnt1,
                                        cnt2,
                                        &min_k_real_m1,
                                        &max_k_real_m1,
                                        &min_l_real_m1,
                                        &max_l_real_m1
                                        );
#ifdef COUNT_STATES
              matrices->N_M[ij][cnt1][cnt2 / 2] = matrices->N_M1[ij][cnt1][cnt2 / 2] = matrices->N_C[ij][cnt1][cnt2 / 2];
#endif
            }
          }
        }
      }

      /* 2nd E_M[ij] = MIN(E_M[ij], E_M[i+1,j] + c) */
      if (matrices->E_M_rem[my_iindx[i + 1] - j] != INF) {
        matrices->E_M_rem[ij] = MIN2(matrices->E_M_rem[ij],
                                     matrices->E_M_rem[my_iindx[i + 1] - j] + P->MLbase
                                     );
      }

      if (matrices->E_M[my_iindx[i + 1] - j]) {
        for (cnt1 = matrices->k_min_M[my_iindx[i + 1] - j];
             cnt1 <= matrices->k_max_M[my_iindx[i + 1] - j];
             cnt1++) {
          for (cnt2 = matrices->l_min_M[my_iindx[i + 1] - j][cnt1];
               cnt2 <= matrices->l_max_M[my_iindx[i + 1] - j][cnt1];
               cnt2 += 2) {
            if (matrices->E_M[my_iindx[i + 1] - j][cnt1][cnt2 / 2] != INF) {
              if (((cnt1 + dia) <= maxD1) && ((cnt2 + dib) <= maxD2)) {
                matrices->E_M[ij][cnt1 + dia][(cnt2 + dib) / 2] = MIN2(matrices->E_M[ij][cnt1 + dia][(cnt2 + dib) / 2],
                                                                       matrices->E_M[my_iindx[i + 1] - j][cnt1][cnt2 / 2] + P->MLbase
                                                                       );
                updatePosteriorBoundaries(cnt1 + dia,
                                          cnt2 + dib,
                                          &min_k_real_m,
                                          &max_k_real_m,
                                          &min_l_real_m,
                                          &max_l_real_m
                                          );
#ifdef COUNT_STATES
                matrices->N_M[ij][cnt1 + dia][(cnt2 + dib) / 2] += matrices->N_M[my_iindx[i + 1] - j][cnt1][cnt2 / 2];
#endif
              }
              /* collect all cases where dia+cnt1 or dib+cnt2 exceeds maxD1, maxD2, respectively */
              else {
                matrices->E_M_rem[ij] = MIN2(matrices->E_M_rem[ij],
                                             matrices->E_M[my_iindx[i + 1] - j][cnt1][cnt2 / 2] + P->MLbase
                                             );
              }
            }
          }
        }
      }

      /* 3rd E_M[ij] = MIN(E_M[ij], E_M[i,j-1] + c) */
      if (matrices->E_M_rem[ij + 1] != INF) {
        matrices->E_M_rem[ij] = MIN2(matrices->E_M_rem[ij],
                                     matrices->E_M_rem[ij + 1] + P->MLbase
                                     );
      }

      if (matrices->E_M[ij + 1]) {
        for (cnt1 = matrices->k_min_M[ij + 1];
             cnt1 <= matrices->k_max_M[ij + 1];
             cnt1++) {
          for (cnt2 = matrices->l_min_M[ij + 1][cnt1];
               cnt2 <= matrices->l_max_M[ij + 1][cnt1];
               cnt2 += 2) {
            if (matrices->E_M[ij + 1][cnt1][cnt2 / 2] != INF) {
              if (((cnt1 + dja) <= maxD1) && ((cnt2 + djb) <= maxD2)) {
                matrices->E_M[ij][cnt1 + dja][(cnt2 + djb) / 2] = MIN2(matrices->E_M[ij][cnt1 + dja][(cnt2 + djb) / 2],
                                                                       matrices->E_M[ij + 1][cnt1][cnt2 / 2] + P->MLbase
                                                                       );
                updatePosteriorBoundaries(cnt1 + dja,
                                          cnt2 + djb,
                                          &min_k_real_m,
                                          &max_k_real_m,
                                          &min_l_real_m,
                                          &max_l_real_m
                                          );
#ifdef COUNT_STATES
                matrices->N_M[ij][cnt1 + dja][(cnt2 + djb) / 2] += matrices->N_M[ij + 1][cnt1][cnt2 / 2];
#endif
              }
              /* collect all cases where dja+cnt1 or djb+cnt2 exceeds maxD1, maxD2, respectively */
              else {
                matrices->E_M_rem[ij] = MIN2(matrices->E_M_rem[ij],
                                             matrices->E_M[ij + 1][cnt1][cnt2 / 2] + P->MLbase
                                             );
              }
            }
          }
        }
      }

      /* 4th E_M1[ij] = MIN(E_M1[ij], E_M1[i,j-1] + c) */
      if (matrices->E_M1_rem[ij + 1] != INF) {
        matrices->E_M1_rem[ij] = MIN2(matrices->E_M1_rem[ij],
                                      matrices->E_M1_rem[ij + 1] + P->MLbase
                                      );
      }

      if (matrices->E_M1[ij + 1]) {
        for (cnt1 = matrices->k_min_M1[ij + 1];
             cnt1 <= matrices->k_max_M1[ij + 1];
             cnt1++) {
          for (cnt2 = matrices->l_min_M1[ij + 1][cnt1];
               cnt2 <= matrices->l_max_M1[ij + 1][cnt1];
               cnt2 += 2) {
            if (matrices->E_M1[ij + 1][cnt1][cnt2 / 2] != INF) {
              if (((cnt1 + dja) <= maxD1) && ((cnt2 + djb) <= maxD2)) {
                matrices->E_M1[ij][cnt1 + dja][(cnt2 + djb) / 2] = MIN2(matrices->E_M1[ij][cnt1 + dja][(cnt2 + djb) / 2],
                                                                        matrices->E_M1[ij + 1][cnt1][cnt2 / 2] + P->MLbase
                                                                        );
                updatePosteriorBoundaries(cnt1 + dja,
                                          cnt2 + djb,
                                          &min_k_real_m1,
                                          &max_k_real_m1,
                                          &min_l_real_m1,
                                          &max_l_real_m1
                                          );
#ifdef COUNT_STATES
                matrices->N_M1[ij][cnt1 + dja][(cnt2 + djb) / 2] += matrices->N_M1[ij + 1][cnt1][cnt2 / 2];
#endif
              }
              /* collect all cases where dja+cnt1 or djb+cnt2 exceeds maxD1, maxD2, respectively */
              else {
                matrices->E_M1_rem[ij] = MIN2(matrices->E_M1_rem[ij],
                                              matrices->E_M1[ij + 1][cnt1][cnt2 / 2] + P->MLbase
                                              );
              }
            }
          }
        }
      }

      /* 5th E_M[ij] = MIN(E_M[ij], min(E_M[i,k] + E_M[k+1,j])) */
      if (j > TURN + 2) {
        for (u = i + 1 + TURN; u <= j - 2 - TURN; u++) {
          /* check all cases where M(i,u) and/or M(u+1,j) are already out of scope of maxD1 and/or maxD2 */
          if (matrices->E_M_rem[my_iindx[i] - u] != INF) {
            for (cnt3 = matrices->k_min_M[my_iindx[u + 1] - j];
                 cnt3 <= matrices->k_max_M[my_iindx[u + 1] - j];
                 cnt3++) {
              for (cnt4 = matrices->l_min_M[my_iindx[u + 1] - j][cnt3];
                   cnt4 <= matrices->l_max_M[my_iindx[u + 1] - j][cnt3];
                   cnt4 += 2) {
                if (matrices->E_M[my_iindx[u + 1] - j][cnt3][cnt4 / 2] != INF) {
                  matrices->E_M_rem[ij] = MIN2(matrices->E_M_rem[ij],
                                               matrices->E_M_rem[my_iindx[i] - u] + matrices->E_M[my_iindx[u + 1] - j][cnt3][cnt4 / 2]
                                               );
                }
              }
            }
            if (matrices->E_M_rem[my_iindx[u + 1] - j] != INF) {
              matrices->E_M_rem[ij] = MIN2(matrices->E_M_rem[ij],
                                           matrices->E_M_rem[my_iindx[i] - u] + matrices->E_M_rem[my_iindx[u + 1] - j]
                                           );
            }
          }

          if (matrices->E_M_rem[my_iindx[u + 1] - j] != INF) {
            for (cnt1 = matrices->k_min_M[my_iindx[i] - u];
                 cnt1 <= matrices->k_max_M[my_iindx[i] - u];
                 cnt1++) {
              for (cnt2 = matrices->l_min_M[my_iindx[i] - u][cnt1];
                   cnt2 <= matrices->l_max_M[my_iindx[i] - u][cnt1];
                   cnt2 += 2) {
                if (matrices->E_M[my_iindx[i] - u][cnt1][cnt2 / 2] != INF) {
                  matrices->E_M_rem[ij] = MIN2(matrices->E_M_rem[ij],
                                               matrices->E_M[my_iindx[i] - u][cnt1][cnt2 / 2] + matrices->E_M_rem[my_iindx[u + 1] - j]
                                               );
                }
              }
            }
          }

          if (!matrices->E_M[my_iindx[i] - u])
            continue;

          if (!matrices->E_M[my_iindx[u + 1] - j])
            continue;

          dia = referenceBPs1[ij] - referenceBPs1[my_iindx[i] - u] - referenceBPs1[my_iindx[u + 1] - j];
          dib = referenceBPs2[ij] - referenceBPs2[my_iindx[i] - u] - referenceBPs2[my_iindx[u + 1] - j];

          for (cnt1 = matrices->k_min_M[my_iindx[i] - u];
               cnt1 <= matrices->k_max_M[my_iindx[i] - u];
               cnt1++) {
            for (cnt2 = matrices->l_min_M[my_iindx[i] - u][cnt1];
                 cnt2 <= matrices->l_max_M[my_iindx[i] - u][cnt1];
                 cnt2 += 2) {
              for (cnt3 = matrices->k_min_M[my_iindx[u + 1] - j];
                   cnt3 <= matrices->k_max_M[my_iindx[u + 1] - j];
                   cnt3++) {
                for (cnt4 = matrices->l_min_M[my_iindx[u + 1] - j][cnt3];
                     cnt4 <= matrices->l_max_M[my_iindx[u + 1] - j][cnt3];
                     cnt4 += 2) {
                  if ((matrices->E_M[my_iindx[i] - u][cnt1][cnt2 / 2] != INF) && (matrices->E_M[my_iindx[u + 1] - j][cnt3][cnt4 / 2] != INF)) {
                    if (((cnt1 + cnt3 + dia) <= maxD1) && ((cnt2 + cnt4 + dib) <= maxD2)) {
                      matrices->E_M[ij][cnt1 + cnt3 + dia][(cnt2 + cnt4 + dib) / 2] = MIN2(matrices->E_M[ij][cnt1 + cnt3 + dia][(cnt2 + cnt4 + dib) / 2],
                                                                                           matrices->E_M[my_iindx[i] - u][cnt1][cnt2 / 2]
                                                                                           + matrices->E_M[my_iindx[u + 1] - j][cnt3][cnt4 / 2]
                                                                                           );
                      updatePosteriorBoundaries(cnt1 + cnt3 + dia,
                                                cnt2 + cnt4 + dib,
                                                &min_k_real_m,
                                                &max_k_real_m,
                                                &min_l_real_m,
                                                &max_l_real_m
                                                );
#ifdef COUNT_STATES
                      matrices->N_M[ij][cnt1 + cnt3 + dia][(cnt2 + cnt4 + dib) / 2] += matrices->N_M[my_iindx[i] - u][cnt1][cnt2 / 2] * matrices->N_M1[my_iindx[u + 1] - j][cnt3][cnt4 / 2];
#endif
                    }
                    /* collect all cases where dia+cnt1+cnt3 or dib+cnt2+cnt4 exceeds maxD1, maxD2, respectively */
                    else {
                      matrices->E_M_rem[ij] = MIN2(matrices->E_M_rem[ij],
                                                   matrices->E_M[my_iindx[i] - u][cnt1][cnt2 / 2] + matrices->E_M[my_iindx[u + 1] - j][cnt3][cnt4 / 2]
                                                   );
                    }
                  }
                }
              }
            }
          }
        }
      }

      /* thats all folks for the multiloop decomposition... */

      adjustArrayBoundaries(&matrices->E_M[ij],
                            &matrices->k_min_M[ij],
                            &matrices->k_max_M[ij],
                            &matrices->l_min_M[ij],
                            &matrices->l_max_M[ij],
                            min_k_real_m,
                            max_k_real_m,
                            min_l_real_m,
                            max_l_real_m
                            );

      adjustArrayBoundaries(&matrices->E_M1[ij],
                            &matrices->k_min_M1[ij],
                            &matrices->k_max_M1[ij],
                            &matrices->l_min_M1[ij],
                            &matrices->l_max_M1[ij],
                            min_k_real_m1,
                            max_k_real_m1,
                            min_l_real_m1,
                            max_l_real_m1
                            );

#ifdef COUNT_STATES
      /* actually we should adjust the array boundaries here but we might never use the count states option more than once so what....*/
#endif
    } /* end of j-loop */
  }

  /* calculate energies of 5' and 3' fragments */

  /* prepare first entries in E_F5 */
  for (cnt1 = 1; cnt1 <= TURN + 1; cnt1++) {
    matrices->E_F5[cnt1]        = (int **)vrna_alloc(sizeof(int *));
    matrices->E_F5[cnt1][0]     = (int *)vrna_alloc(sizeof(int));
    matrices->E_F5[cnt1][0][0]  = 0;
    matrices->E_F5_rem[cnt1]    = INF;
    matrices->k_min_F5[cnt1]    = matrices->k_max_F5[cnt1] = 0;
    matrices->l_min_F5[cnt1]    = (int *)vrna_alloc(sizeof(int));
    matrices->l_max_F5[cnt1]    = (int *)vrna_alloc(sizeof(int));
    matrices->l_min_F5[cnt1][0] = matrices->l_max_F5[cnt1][0] = 0;
#ifdef COUNT_STATES
    matrices->N_F5[cnt1]        = (unsigned long **)vrna_alloc(sizeof(unsigned long *));
    matrices->N_F5[cnt1][0]     = (unsigned long *)vrna_alloc(sizeof(unsigned long));
    matrices->N_F5[cnt1][0][0]  = 1;
#endif
  }


  for (j = TURN + 2; j <= seq_length; j++) {
    unsigned int  da  = referenceBPs1[my_iindx[1] - j] - referenceBPs1[my_iindx[1] - j + 1];
    unsigned int  db  = referenceBPs2[my_iindx[1] - j] - referenceBPs2[my_iindx[1] - j + 1];

    type          = ptype[jindx[j] + 1];
    additional_en = 0;
    if (type) {
      if (dangles == 2)
        additional_en += E_ExtLoop(type, -1, j < seq_length ? S1[j + 1] : -1, P);
      else
        additional_en += E_ExtLoop(type, -1, -1, P);
    }

    /* make min and max k guess for memory allocation */
    int min_k_guess, max_k_guess, min_l_guess, max_l_guess;
    int *min_l_real, *max_l_real, min_k_real, max_k_real;

    min_k_guess = min_l_guess = 0;
    max_k_guess = referenceBPs1[my_iindx[1] - j] + mm1[my_iindx[1] - j];
    max_l_guess = referenceBPs2[my_iindx[1] - j] + mm2[my_iindx[1] - j];

    prepareBoundaries(min_k_guess,
                      max_k_guess,
                      min_l_guess,
                      max_l_guess,
                      bpdist[my_iindx[1] - j],
                      &matrices->k_min_F5[j],
                      &matrices->k_max_F5[j],
                      &matrices->l_min_F5[j],
                      &matrices->l_max_F5[j]
                      );

    preparePosteriorBoundaries(matrices->k_max_F5[j] - matrices->k_min_F5[j] + 1,
                               matrices->k_min_F5[j],
                               &min_k_real,
                               &max_k_real,
                               &min_l_real,
                               &max_l_real
                               );

    prepareArray(&matrices->E_F5[j],
                 matrices->k_min_F5[j],
                 matrices->k_max_F5[j],
                 matrices->l_min_F5[j],
                 matrices->l_max_F5[j]
                 );
#ifdef COUNT_STATES
    prepareArray2(&matrices->N_F5[j],
                  matrices->k_min_F5[j],
                  matrices->k_max_F5[j],
                  matrices->l_min_F5[j],
                  matrices->l_max_F5[j]
                  );
#endif

    /* begin the actual computation of 5' end energies */

    /* j-1 is unpaired ... */
    matrices->E_F5_rem[j] = matrices->E_F5_rem[j - 1];
    for (cnt1 = matrices->k_min_F5[j - 1]; cnt1 <= matrices->k_max_F5[j - 1]; cnt1++) {
      for (cnt2 = matrices->l_min_F5[j - 1][cnt1]; cnt2 <= matrices->l_max_F5[j - 1][cnt1]; cnt2 += 2) {
        if (((cnt1 + da) <= maxD1) && ((cnt2 + db) <= maxD2)) {
          matrices->E_F5[j][cnt1 + da][(cnt2 + db) / 2] = MIN2(matrices->E_F5[j][cnt1 + da][(cnt2 + db) / 2],
                                                               matrices->E_F5[j - 1][cnt1][cnt2 / 2]
                                                               );
          updatePosteriorBoundaries(cnt1 + da,
                                    cnt2 + db,
                                    &min_k_real,
                                    &max_k_real,
                                    &min_l_real,
                                    &max_l_real
                                    );
#ifdef COUNT_STATES
          matrices->N_F5[j][cnt1 + da][(cnt2 + db) / 2] += matrices->N_F5[j - 1][cnt1][cnt2 / 2];
#endif
        }
        /* collect all cases where da+cnt1 or db+cnt2 exceeds maxD1, maxD2, respectively */
        else {
          matrices->E_F5_rem[j] = MIN2(matrices->E_F5_rem[j], matrices->E_F5[j - 1][cnt1][cnt2 / 2]);
        }
      }
    }
    /* j pairs with 1 */
    if (matrices->E_C_rem[my_iindx[1] - j] != INF)
      matrices->E_F5_rem[j] = MIN2(matrices->E_F5_rem[j], matrices->E_C_rem[my_iindx[1] - j] + additional_en);

    if (matrices->E_C[my_iindx[1] - j]) {
      for (cnt1 = matrices->k_min_C[my_iindx[1] - j]; cnt1 <= matrices->k_max_C[my_iindx[1] - j]; cnt1++)
        for (cnt2 = matrices->l_min_C[my_iindx[1] - j][cnt1]; cnt2 <= matrices->l_max_C[my_iindx[1] - j][cnt1]; cnt2 += 2) {
          if (matrices->E_C[my_iindx[1] - j][cnt1][cnt2 / 2] != INF) {
            matrices->E_F5[j][cnt1][cnt2 / 2] = MIN2(matrices->E_F5[j][cnt1][cnt2 / 2],
                                                     matrices->E_C[my_iindx[1] - j][cnt1][cnt2 / 2] + additional_en
                                                     );
            updatePosteriorBoundaries(cnt1,
                                      cnt2,
                                      &min_k_real,
                                      &max_k_real,
                                      &min_l_real,
                                      &max_l_real
                                      );
#ifdef COUNT_STATES
            matrices->N_F5[j][cnt1][cnt2 / 2] += matrices->N_C[my_iindx[1] - j][cnt1][cnt2 / 2];
#endif
          }
        }
    }

    /* j pairs with some other nucleotide -> see below */
    for (i = j - TURN - 1; i > 1; i--) {
      ij    = my_iindx[i] - j;
      type  = ptype[jindx[j] + i];
      if (type) {
        if (dangles == 2)
          additional_en = E_ExtLoop(type, S1[i - 1], j < seq_length ? S1[j + 1] : -1, P);
        else
          additional_en = E_ExtLoop(type, -1, -1, P);

        if (matrices->E_C_rem[ij] != INF) {
          for (cnt3 = matrices->k_min_F5[i - 1]; cnt3 <= matrices->k_max_F5[i - 1]; cnt3++)
            for (cnt4 = matrices->l_min_F5[i - 1][cnt3]; cnt4 <= matrices->l_max_F5[i - 1][cnt3]; cnt4 += 2) {
              if (matrices->E_F5[i - 1][cnt3][cnt4 / 2] != INF) {
                matrices->E_F5_rem[j] = MIN2(matrices->E_F5_rem[j],
                                             matrices->E_F5[i - 1][cnt3][cnt4 / 2] + matrices->E_C_rem[ij] + additional_en
                                             );
              }
            }
          if (matrices->E_F5_rem[i - 1] != INF) {
            matrices->E_F5_rem[j] = MIN2(matrices->E_F5_rem[j],
                                         matrices->E_F5_rem[i - 1] + matrices->E_C_rem[ij] + additional_en
                                         );
          }
        }

        if ((matrices->E_F5_rem[i - 1] != INF) && (matrices->E_C[ij])) {
          for (cnt1 = matrices->k_min_C[ij]; cnt1 <= matrices->k_max_C[ij]; cnt1++)
            for (cnt2 = matrices->l_min_C[ij][cnt1]; cnt2 <= matrices->l_max_C[ij][cnt1]; cnt2 += 2)
              if (matrices->E_C[ij][cnt1][cnt2 / 2] != INF) {
                matrices->E_F5_rem[j] = MIN2(matrices->E_F5_rem[j],
                                             matrices->E_F5_rem[i - 1] + matrices->E_C[ij][cnt1][cnt2 / 2] + additional_en
                                             );
              }
        }

        if (!matrices->E_C[ij])
          continue;

        unsigned int  d1a = referenceBPs1[my_iindx[1] - j] - referenceBPs1[ij] - referenceBPs1[my_iindx[1] - i + 1];
        unsigned int  d1b = referenceBPs2[my_iindx[1] - j] - referenceBPs2[ij] - referenceBPs2[my_iindx[1] - i + 1];

        for (cnt1 = matrices->k_min_C[ij]; cnt1 <= matrices->k_max_C[ij]; cnt1++)
          for (cnt2 = matrices->l_min_C[ij][cnt1]; cnt2 <= matrices->l_max_C[ij][cnt1]; cnt2 += 2)
            for (cnt3 = matrices->k_min_F5[i - 1]; cnt3 <= matrices->k_max_F5[i - 1]; cnt3++)
              for (cnt4 = matrices->l_min_F5[i - 1][cnt3]; cnt4 <= matrices->l_max_F5[i - 1][cnt3]; cnt4 += 2) {
                if (matrices->E_F5[i - 1][cnt3][cnt4 / 2] != INF && matrices->E_C[ij][cnt1][cnt2 / 2] != INF) {
                  if (((cnt1 + cnt3 + d1a) <= maxD1) && ((cnt2 + cnt4 + d1b) <= maxD2)) {
                    matrices->E_F5[j][cnt1 + cnt3 + d1a][(cnt2 + cnt4 + d1b) / 2] = MIN2(matrices->E_F5[j][cnt1 + cnt3 + d1a][(cnt2 + cnt4 + d1b) / 2],
                                                                                         matrices->E_F5[i - 1][cnt3][cnt4 / 2] + matrices->E_C[ij][cnt1][cnt2 / 2] + additional_en
                                                                                         );
                    updatePosteriorBoundaries(cnt1 + cnt3 + d1a,
                                              cnt2 + cnt4 + d1b,
                                              &min_k_real,
                                              &max_k_real,
                                              &min_l_real,
                                              &max_l_real
                                              );
#ifdef COUNT_STATES
                    matrices->N_F5[j][cnt1 + cnt3 + d1a][(cnt2 + cnt4 + d1b) / 2] += matrices->N_F5[i - 1][cnt3][cnt4 / 2] * matrices->N_C[ij][cnt1][cnt2 / 2];
#endif
                  }
                  /* collect all cases where d1a+cnt1+cnt3 or d1b+cnt2+cnt4 exceeds maxD1, maxD2, respectively */
                  else {
                    matrices->E_F5_rem[j] = MIN2(matrices->E_F5_rem[j],
                                                 matrices->E_F5[i - 1][cnt3][cnt4 / 2] + matrices->E_C[ij][cnt1][cnt2 / 2] + additional_en
                                                 );
                  }
                }
              }
      }
    }

    /* resize and move memory portions of energy matrix E_F5 */
    adjustArrayBoundaries(&matrices->E_F5[j],
                          &matrices->k_min_F5[j],
                          &matrices->k_max_F5[j],
                          &matrices->l_min_F5[j],
                          &matrices->l_max_F5[j],
                          min_k_real,
                          max_k_real,
                          min_l_real,
                          max_l_real
                          );
  } /* end of j-loop */


  if (compute_2Dfold_F3) {
    /* prepare first entries in E_F3 */
    for (cnt1 = seq_length; cnt1 >= seq_length - TURN - 1; cnt1--) {
      matrices->E_F3[cnt1]        = (int **)vrna_alloc(sizeof(int *));
      matrices->E_F3[cnt1][0]     = (int *)vrna_alloc(sizeof(int));
      matrices->E_F3[cnt1][0][0]  = 0;
      matrices->k_min_F3[cnt1]    = matrices->k_max_F3[cnt1] = 0;
      matrices->l_min_F3[cnt1]    = (int *)vrna_alloc(sizeof(int));
      matrices->l_max_F3[cnt1]    = (int *)vrna_alloc(sizeof(int));
      matrices->l_min_F3[cnt1][0] = matrices->l_max_F3[cnt1][0] = 0;
    }
    /* begin calculations */
    for (j = seq_length - TURN - 2; j >= 1; j--) {
      unsigned int  da  = referenceBPs1[my_iindx[j] - seq_length] - referenceBPs1[my_iindx[j + 1] - seq_length];
      unsigned int  db  = referenceBPs2[my_iindx[j] - seq_length] - referenceBPs2[my_iindx[j + 1] - seq_length];

      type          = ptype[jindx[seq_length] + j];
      additional_en = 0;
      if (type) {
        if (dangles == 2)
          additional_en += E_ExtLoop(type, j > 1 ? S1[j - 1] : -1, -1, P);
        else
          additional_en += E_ExtLoop(type, -1, -1, P);
      }

      /* make min and max k guess for memory allocation */
      int min_k_guess, max_k_guess, min_l_guess, max_l_guess;
      int *min_l_real, *max_l_real, min_k_real, max_k_real;

      min_k_guess = min_l_guess = 0;
      max_k_guess = referenceBPs1[my_iindx[j] - seq_length] + mm1[my_iindx[j] - seq_length];
      max_l_guess = referenceBPs2[my_iindx[j] - seq_length] + mm2[my_iindx[j] - seq_length];

      prepareBoundaries(min_k_guess,
                        max_k_guess,
                        min_l_guess,
                        max_l_guess,
                        bpdist[my_iindx[j] - seq_length],
                        &matrices->k_min_F3[j],
                        &matrices->k_max_F3[j],
                        &matrices->l_min_F3[j],
                        &matrices->l_max_F3[j]
                        );

      preparePosteriorBoundaries(matrices->k_max_F3[j] - matrices->k_min_F3[j] + 1,
                                 matrices->k_min_F3[j],
                                 &min_k_real,
                                 &max_k_real,
                                 &min_l_real,
                                 &max_l_real
                                 );

      prepareArray(&matrices->E_F3[j],
                   matrices->k_min_F3[j],
                   matrices->k_max_F3[j],
                   matrices->l_min_F3[j],
                   matrices->l_max_F3[j]
                   );
      /* begin the actual computation of 5' end energies */

      /* j is unpaired ... */
      for (cnt1 = matrices->k_min_F3[j + 1]; cnt1 <= matrices->k_max_F3[j + 1]; cnt1++) {
        for (cnt2 = matrices->l_min_F3[j + 1][cnt1]; cnt2 <= matrices->l_max_F3[j + 1][cnt1]; cnt2 += 2) {
          matrices->E_F3[j][cnt1 + da][(cnt2 + db) / 2] = MIN2(matrices->E_F3[j][cnt1 + da][(cnt2 + db) / 2],
                                                               matrices->E_F3[j + 1][cnt1][cnt2 / 2]
                                                               );
          updatePosteriorBoundaries(cnt1 + da,
                                    cnt2 + db,
                                    &min_k_real,
                                    &max_k_real,
                                    &min_l_real,
                                    &max_l_real
                                    );
        }
      }
      /* j pairs with n */
      if (matrices->E_C[my_iindx[j] - seq_length]) {
        for (cnt1 = matrices->k_min_C[my_iindx[j] - seq_length]; cnt1 <= matrices->k_max_C[my_iindx[j] - seq_length]; cnt1++)
          for (cnt2 = matrices->l_min_C[my_iindx[j] - seq_length][cnt1]; cnt2 <= matrices->l_max_C[my_iindx[j] - seq_length][cnt1]; cnt2 += 2) {
            if (matrices->E_C[my_iindx[j] - seq_length][cnt1][cnt2 / 2] != INF) {
              matrices->E_F3[j][cnt1][cnt2 / 2] = MIN2(matrices->E_F3[j][cnt1][cnt2 / 2],
                                                       matrices->E_C[my_iindx[j] - seq_length][cnt1][cnt2 / 2] + additional_en
                                                       );
              updatePosteriorBoundaries(cnt1,
                                        cnt2,
                                        &min_k_real,
                                        &max_k_real,
                                        &min_l_real,
                                        &max_l_real
                                        );
            }
          }
      }

      /* j pairs with some other nucleotide -> see below */
      for (i = j - TURN - 1; i > 1; i--) {
        ij = my_iindx[i] - j;
        if (!matrices->E_C[ij])
          continue;

        type = ptype[jindx[j] + i];
        if (type) {
          unsigned int  d1a = referenceBPs1[my_iindx[1] - j] - referenceBPs1[ij] - referenceBPs1[my_iindx[1] - i + 1];
          unsigned int  d1b = referenceBPs2[my_iindx[1] - j] - referenceBPs2[ij] - referenceBPs2[my_iindx[1] - i + 1];

          if (dangles == 2)
            additional_en = E_ExtLoop(type, S1[i - 1], j < seq_length ? S1[j + 1] : -1, P);
          else
            additional_en = E_ExtLoop(type, -1, -1, P);

          for (cnt1 = matrices->k_min_C[ij]; cnt1 <= matrices->k_max_C[ij]; cnt1++)
            for (cnt2 = matrices->l_min_C[ij][cnt1]; cnt2 <= matrices->l_max_C[ij][cnt1]; cnt2 += 2)
              for (cnt3 = matrices->k_min_F5[i - 1]; cnt3 <= matrices->k_max_F5[i - 1]; cnt3++)
                for (cnt4 = matrices->l_min_F5[i - 1][cnt3]; cnt4 <= matrices->l_max_F5[i - 1][cnt3]; cnt4 += 2) {
                  if (matrices->E_F5[i - 1][cnt3][cnt4 / 2] != INF && matrices->E_C[ij][cnt1][cnt2 / 2] != INF) {
                    matrices->E_F5[j][cnt1 + cnt3 + d1a][(cnt2 + cnt4 + d1b) / 2] = MIN2(matrices->E_F5[j][cnt1 + cnt3 + d1a][(cnt2 + cnt4 + d1b) / 2],
                                                                                         matrices->E_F5[i - 1][cnt3][cnt4 / 2] + matrices->E_C[ij][cnt1][cnt2 / 2] + additional_en
                                                                                         );
                    updatePosteriorBoundaries(cnt1 + cnt3 + d1a,
                                              cnt2 + cnt4 + d1b,
                                              &min_k_real,
                                              &max_k_real,
                                              &min_l_real,
                                              &max_l_real
                                              );
#ifdef COUNT_STATES
                    matrices->N_F5[j][cnt1 + cnt3 + d1a][(cnt2 + cnt4 + d1b) / 2] += matrices->N_F5[i - 1][cnt3][cnt4 / 2] * matrices->N_C[ij][cnt1][cnt2 / 2];
#endif
                  }
                }
        }
      }

      /* resize and move memory portions of energy matrix E_F5 */
      adjustArrayBoundaries(&matrices->E_F5[j],
                            &matrices->k_min_F5[j],
                            &matrices->k_max_F5[j],
                            &matrices->l_min_F5[j],
                            &matrices->l_max_F5[j],
                            min_k_real,
                            max_k_real,
                            min_l_real,
                            max_l_real
                            );
    } /* end of j-loop */
  }
}


/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/

PRIVATE void
backtrack_f5(unsigned int         j,
             int                  k,
             int                  l,
             char                 *structure,
             vrna_fold_compound_t *vc)
{
  int           *my_iindx, *jindx, energy, type, dangles, cnt1, cnt2, cnt3, cnt4;
  int           **l_min_C, **l_max_C, **l_min_F5, **l_max_F5;
  int           *k_min_C, *k_max_C, *k_min_F5, *k_max_F5;
  int           ***E_C, ***E_F5;
  int           *E_C_rem, *E_F5_rem;
  unsigned int  i, ij, seq_length, maxD1, maxD2;
  short         *S1;
  unsigned int  *referenceBPs1, *referenceBPs2;
  char          *ptype;
  vrna_param_t  *P;
  vrna_md_t     *md;
  vrna_mx_mfe_t *matrices;
  unsigned int  da, db;

  P             = vc->params;
  md            = &(P->model_details);
  matrices      = vc->matrices;
  seq_length    = vc->length;
  S1            = vc->sequence_encoding;
  ptype         = vc->ptype;
  my_iindx      = vc->iindx;
  jindx         = vc->jindx;
  referenceBPs1 = vc->referenceBPs1;
  referenceBPs2 = vc->referenceBPs2;
  dangles       = md->dangles;
  E_F5          = matrices->E_F5;
  l_min_F5      = matrices->l_min_F5;
  l_max_F5      = matrices->l_max_F5;
  k_min_F5      = matrices->k_min_F5;
  k_max_F5      = matrices->k_max_F5;

  E_C     = matrices->E_C;
  l_min_C = matrices->l_min_C;
  l_max_C = matrices->l_max_C;
  k_min_C = matrices->k_min_C;
  k_max_C = matrices->k_max_C;

  E_F5_rem  = matrices->E_F5_rem;
  E_C_rem   = matrices->E_C_rem;
  maxD1     = vc->maxD1;
  maxD2     = vc->maxD2;

  da  = referenceBPs1[my_iindx[1] - j] - referenceBPs1[my_iindx[1] - j + 1];
  db  = referenceBPs2[my_iindx[1] - j] - referenceBPs2[my_iindx[1] - j + 1];

  if (j < TURN + 2)
    return;

  /* F5[j] == F5[j-1] ? */
  if (k == -1) {
    if (E_F5_rem[j] == INF) {
      return;
    } else if (E_F5_rem[j] == E_F5_rem[j - 1]) {
      backtrack_f5(j - 1, k, l, structure, vc);
      return;
    } else if (E_F5[j - 1]) {
      for (cnt1 = k_min_F5[j - 1];
           cnt1 <= k_max_F5[j - 1];
           cnt1++) {
        for (cnt2 = l_min_F5[j - 1][cnt1];
             cnt2 <= l_max_F5[j - 1][cnt1];
             cnt2 += 2) {
          if (((cnt1 + da) > maxD1) || ((cnt2 + db) > maxD2)) {
            if (E_F5_rem[j] == E_F5[j - 1][cnt1][cnt2 / 2]) {
              backtrack_f5(j - 1, cnt1, cnt2, structure, vc);
              return;
            }
          }
        }
      }
    }
  } else if ((k >= da) && (l >= db)) {
    if (E_F5[j - 1]) {
      if ((k - da >= k_min_F5[j - 1]) && (k - da <= k_max_F5[j - 1])) {
        if ((l - db >= l_min_F5[j - 1][k - da]) && (l - db <= l_max_F5[j - 1][k - da])) {
          if (E_F5[j - 1][k - da][(l - db) / 2] == E_F5[j][k][l / 2]) {
            backtrack_f5(j - 1, k - da, l - db, structure, vc);
            return;
          }
        }
      }
    }
  }

  type = ptype[jindx[j] + 1];
  if (type) {
    if (dangles == 2)
      energy = E_ExtLoop(type, -1, j < seq_length ? S1[j + 1] : -1, P);
    else
      energy = E_ExtLoop(type, -1, -1, P);

    if (k == -1) {
      if (E_C_rem[my_iindx[1] - j] + energy == E_F5_rem[j]) {
        backtrack_c(1, j, -1, -1, structure, vc);
        return;
      }
    } else if (k >= k_min_C[my_iindx[1] - j] && (k <= k_max_C[my_iindx[1] - j])) {
      if ((l >= l_min_C[my_iindx[1] - j][k]) && (l <= l_max_C[my_iindx[1] - j][k])) {
        if (E_C[my_iindx[1] - j][k][l / 2] + energy == E_F5[j][k][l / 2]) {
          backtrack_c(1, j, k, l, structure, vc);
          return;
        }
      }
    }
  }

  for (i = j - TURN - 1; i > 1; i--) {
    ij    = my_iindx[i] - j;
    type  = ptype[jindx[j] + i];
    if (type) {
      unsigned int  d1a = referenceBPs1[my_iindx[1] - j] - referenceBPs1[ij] - referenceBPs1[my_iindx[1] - i + 1];
      unsigned int  d1b = referenceBPs2[my_iindx[1] - j] - referenceBPs2[ij] - referenceBPs2[my_iindx[1] - i + 1];

      if (dangles == 2)
        energy = E_ExtLoop(type, S1[i - 1], j < seq_length ? S1[j + 1] : -1, P);
      else
        energy = E_ExtLoop(type, -1, -1, P);

      if (k == -1) {
        if (E_C_rem[ij] != INF) {
          for (cnt1 = k_min_F5[i - 1];
               cnt1 <= k_max_F5[i - 1];
               cnt1++) {
            for (cnt2 = l_min_F5[i - 1][cnt1];
                 cnt2 <= l_max_F5[i - 1][cnt1];
                 cnt2 += 2) {
              if (E_F5_rem[j] == (E_F5[i - 1][cnt1][cnt2 / 2] + E_C_rem[ij] + energy)) {
                backtrack_f5(i - 1, cnt1, cnt2, structure, vc);
                backtrack_c(i, j, -1, -1, structure, vc);
                return;
              }
            }
          }
          if (E_F5_rem[j] == (E_F5_rem[i - 1] + E_C_rem[ij] + energy)) {
            backtrack_f5(i - 1, -1, -1, structure, vc);
            backtrack_c(i, j, -1, -1, structure, vc);
            return;
          }
        }

        if (E_F5_rem[i - 1] != INF) {
          for (cnt1 = k_min_C[ij];
               cnt1 <= k_max_C[ij];
               cnt1++) {
            for (cnt2 = l_min_C[ij][cnt1];
                 cnt2 <= l_max_C[ij][cnt1];
                 cnt2 += 2) {
              if (E_F5_rem[j] == (E_F5_rem[i - 1] + E_C[ij][cnt1][cnt2 / 2] + energy)) {
                backtrack_f5(i - 1, -1, -1, structure, vc);
                backtrack_c(i, j, cnt1, cnt2, structure, vc);
                return;
              }
            }
          }
        }

        for (cnt1 = k_min_F5[i - 1];
             cnt1 <= k_max_F5[i - 1];
             cnt1++)
          for (cnt2 = l_min_F5[i - 1][cnt1];
               cnt2 <= l_max_F5[i - 1][cnt1];
               cnt2 += 2)
            for (cnt3 = k_min_C[ij];
                 cnt3 <= k_max_C[ij];
                 cnt3++)
              for (cnt4 = l_min_C[ij][cnt3];
                   cnt4 <= l_max_C[ij][cnt3];
                   cnt4 += 2) {
                if (((cnt1 + cnt3 + d1a) > maxD1) || ((cnt2 + cnt4 + d1b) > maxD2)) {
                  if (E_F5_rem[j] == (E_F5[i - 1][cnt1][cnt2 / 2] + E_C[ij][cnt3][cnt4 / 2] + energy)) {
                    backtrack_f5(i - 1, cnt1, cnt2, structure, vc);
                    backtrack_c(i, j, cnt3, cnt4, structure, vc);
                    return;
                  }
                }
              }
      } else if ((k >= d1a) && (l >= d1b)) {
        int k_f_max = MIN2(k - d1a, k_max_F5[i - 1]);

        for (cnt1 = k_min_F5[i - 1]; cnt1 <= k_f_max; cnt1++) {
          int l_f_max = MIN2(l - d1b, l_max_F5[i - 1][cnt1]);
          for (cnt2 = l_min_F5[i - 1][cnt1]; cnt2 <= l_f_max; cnt2 += 2) {
            int k_c = k - d1a - cnt1;
            if ((k_c >= k_min_C[ij]) && (k_c <= k_max_C[ij])) {
              int l_c = l - d1b - cnt2;
              if ((l_c >= l_min_C[ij][k_c]) && (l_c <= l_max_C[ij][k_c])) {
                if (E_F5[j][k][l / 2] == (E_F5[i - 1][cnt1][cnt2 / 2] + E_C[ij][k_c][l_c / 2] + energy)) {
                  backtrack_f5(i - 1, cnt1, cnt2, structure, vc);
                  backtrack_c(i, j, k_c, l_c, structure, vc);
                  return;
                }
              }
            }
          }
        }
      }
    }
  }
  vrna_message_error("backtracking failed in f5");
}


PRIVATE void
backtrack_c(unsigned int          i,
            unsigned int          j,
            int                   k,
            int                   l,
            char                  *structure,
            vrna_fold_compound_t  *vc)
{
  unsigned int  p, q, pq, ij, maxp, maxD1, maxD2;
  int           *my_iindx, *jindx, type, type_2, energy, no_close, dangles, base_d1, base_d2, d1, d2, cnt1, cnt2, cnt3, cnt4, *rtype;
  int           **l_min_C, **l_max_C, **l_min_M, **l_max_M, **l_min_M1, **l_max_M1;
  int           *k_min_C, *k_max_C, *k_min_M, *k_max_M, *k_min_M1, *k_max_M1;
  int           ***E_C, ***E_M, ***E_M1, *E_C_rem, *E_M_rem, *E_M1_rem;
  short         *S1;
  unsigned int  *referenceBPs1, *referenceBPs2;
  char          *ptype, *sequence;
  vrna_param_t  *P;
  vrna_md_t     *md;
  vrna_mx_mfe_t *matrices;

  P             = vc->params;
  md            = &(P->model_details);
  matrices      = vc->matrices;
  sequence      = vc->sequence;
  S1            = vc->sequence_encoding;
  ptype         = vc->ptype;
  rtype         = &(md->rtype[0]);
  my_iindx      = vc->iindx;
  jindx         = vc->jindx;
  referenceBPs1 = vc->referenceBPs1;
  referenceBPs2 = vc->referenceBPs2;
  dangles       = md->dangles;

  E_C     = matrices->E_C;
  l_min_C = matrices->l_min_C;
  l_max_C = matrices->l_max_C;
  k_min_C = matrices->k_min_C;
  k_max_C = matrices->k_max_C;

  E_M     = matrices->E_M;
  l_min_M = matrices->l_min_M;
  l_max_M = matrices->l_max_M;
  k_min_M = matrices->k_min_M;
  k_max_M = matrices->k_max_M;

  E_M1      = matrices->E_M1;
  l_min_M1  = matrices->l_min_M1;
  l_max_M1  = matrices->l_max_M1;
  k_min_M1  = matrices->k_min_M1;
  k_max_M1  = matrices->k_max_M1;

  E_C_rem   = matrices->E_C_rem;
  E_M_rem   = matrices->E_M_rem;
  E_M1_rem  = matrices->E_M1_rem;
  maxD1     = vc->maxD1;
  maxD2     = vc->maxD2;


  ij = my_iindx[i] - j;

  int e = (k == -1) ? E_C_rem[ij] : E_C[ij][k][l / 2];

  type = ptype[jindx[j] + i];

  no_close          = (((type == 3) || (type == 4)) && no_closingGU);
  structure[i - 1]  = '(';
  structure[j - 1]  = ')';

  base_d1 = ((unsigned int)vc->reference_pt1[i] != j) ? 1 : -1;
  base_d2 = ((unsigned int)vc->reference_pt2[i] != j) ? 1 : -1;

  base_d1 += referenceBPs1[ij];
  base_d2 += referenceBPs2[ij];

  if (k == -1) {
    if (((unsigned int)base_d1 > maxD1) || ((unsigned int)base_d2 > maxD2))
      if (e == E_Hairpin(j - i - 1, type, S1[i + 1], S1[j - 1], sequence + i - 1, P))
        return;
  } else {
    if ((unsigned int)base_d1 == k)
      if ((unsigned int)base_d2 == l)
        if (E_Hairpin(j - i - 1, type, S1[i + 1], S1[j - 1], sequence + i - 1, P) == e)
          return;
  }

  maxp = MIN2(j - 2 - TURN, i + MAXLOOP + 1);
  for (p = i + 1; p <= maxp; p++) {
    unsigned int minq, ln_pre;
    minq    = p + TURN + 1;
    ln_pre  = j - i - 1;
    if (ln_pre > minq + MAXLOOP)
      minq = ln_pre - MAXLOOP - 1;

    for (q = minq; q < j; q++) {
      pq      = my_iindx[p] - q;
      type_2  = ptype[jindx[q] + p];
      if (type_2 == 0)
        continue;

      type_2 = rtype[type_2];

      /* d2 = dbp(S_{i,j}, S_{p.q} + {i,j}) */
      d1  = base_d1 - referenceBPs1[pq];
      d2  = base_d2 - referenceBPs2[pq];

      energy = E_IntLoop(p - i - 1, j - q - 1, type, type_2, S1[i + 1], S1[j - 1], S1[p - 1], S1[q + 1], P);


      if (k == -1) {
        if (E_C_rem[pq] != INF) {
          if (e == (E_C_rem[pq] + energy)) {
            backtrack_c(p, q, -1, -1, structure, vc);
            return;
          }
        }

        if (E_C[pq]) {
          for (cnt1 = k_min_C[pq];
               cnt1 <= k_max_C[pq];
               cnt1++)
            for (cnt2 = l_min_C[pq][cnt1];
                 cnt2 <= l_max_C[pq][cnt1];
                 cnt2 += 2) {
              if (((cnt1 + d1) > maxD1) || ((cnt2 + d2) > maxD2)) {
                if (e == (E_C[pq][cnt1][cnt2 / 2] + energy)) {
                  backtrack_c(p, q, cnt1, cnt2, structure, vc);
                  return;
                }
              }
            }
        }
      } else {
        if (!E_C[pq])
          continue;

        if (d1 <= k && d2 <= l) {
          if ((k - d1 >= k_min_C[pq]) && (k - d1) <= k_max_C[pq]) {
            if ((l - d2 >= l_min_C[pq][k - d1]) && (l - d2 <= l_max_C[pq][k - d1])) {
              if (E_C[pq][k - d1][(l - d2) / 2] + energy == e) {
                backtrack_c(p, q, k - d1, l - d2, structure, vc);
                return;
              }
            }
          }
        }
      }
    } /* end q-loop */
  }   /* end p-loop */

  /* multi-loop decomposition ------------------------*/
  if (!no_close) {
    unsigned int  u;
    int           tt;
    if (k == -1) {
      for (u = i + TURN + 2; u < j - TURN - 2; u++) {
        int i1u, u1j1;
        i1u     = my_iindx[i + 1] - u;
        u1j1    = my_iindx[u + 1] - j + 1;
        tt      = rtype[type];
        energy  = P->MLclosing;
        if (dangles == 2)
          energy += E_MLstem(tt, S1[j - 1], S1[i + 1], P);
        else
          energy += E_MLstem(tt, -1, -1, P);

        if (E_M_rem[i1u] != INF) {
          if (E_M1[u1j1]) {
            for (cnt1 = k_min_M1[u1j1];
                 cnt1 <= k_max_M1[u1j1];
                 cnt1++)
              for (cnt2 = l_min_M1[u1j1][cnt1];
                   cnt2 <= l_max_M1[u1j1][cnt1];
                   cnt2 += 2) {
                if (e == (E_M_rem[i1u] + E_M1[u1j1][cnt1][cnt2 / 2] + energy)) {
                  backtrack_m(i + 1, u, -1, -1, structure, vc);
                  backtrack_m1(u + 1, j - 1, cnt1, cnt2, structure, vc);
                  return;
                }
              }
          }

          if (E_M1_rem[u1j1] != INF) {
            if (e == (E_M_rem[i1u] + E_M1_rem[u1j1] + energy)) {
              backtrack_m(i + 1, u, -1, -1, structure, vc);
              backtrack_m1(u + 1, j - 1, -1, -1, structure, vc);
              return;
            }
          }
        }

        if (E_M1_rem[u1j1] != INF) {
          if (E_M[i1u]) {
            for (cnt1 = k_min_M[i1u];
                 cnt1 <= k_max_M[i1u];
                 cnt1++)
              for (cnt2 = l_min_M[i1u][cnt1];
                   cnt2 <= l_max_M[i1u][cnt1];
                   cnt2 += 2)
                if (e == (E_M[i1u][cnt1][cnt2 / 2] + E_M1_rem[u1j1] + energy)) {
                  backtrack_m(i + 1, u, cnt1, cnt2, structure, vc);
                  backtrack_m1(u + 1, j - 1, -1, -1, structure, vc);
                  return;
                }
          }
        }

        /* now all cases where we exceed the maxD1/D2 scope by combination of E_M and E_M1 */
        if (!E_M[i1u])
          continue;

        if (!E_M1[u1j1])
          continue;

        /* get distance to reference if closing this multiloop
         *  dist3 = dbp(S_{i,j}, {i,j} + S_{i+1.u} + S_{u+1,j-1})
         */
        d1  = base_d1 - referenceBPs1[i1u] - referenceBPs1[u1j1];
        d2  = base_d2 - referenceBPs2[i1u] - referenceBPs2[u1j1];

        for (cnt1 = matrices->k_min_M[i1u];
             cnt1 <= matrices->k_max_M[i1u];
             cnt1++)
          for (cnt2 = matrices->l_min_M[i1u][cnt1];
               cnt2 <= matrices->l_max_M[i1u][cnt1];
               cnt2 += 2)
            for (cnt3 = matrices->k_min_M1[u1j1];
                 cnt3 <= matrices->k_max_M1[u1j1];
                 cnt3++)
              for (cnt4 = matrices->l_min_M1[u1j1][cnt3];
                   cnt4 <= matrices->l_max_M1[u1j1][cnt3];
                   cnt4 += 2) {
                if (((cnt1 + cnt3 + d1) > maxD1) || ((cnt2 + cnt4 + d2) > maxD2)) {
                  if (e == (E_M[i1u][cnt1][cnt2 / 2] + E_M1[u1j1][cnt3][cnt4 / 2] + energy)) {
                    backtrack_m(i + 1, u, cnt1, cnt2, structure, vc);
                    backtrack_m1(u + 1, j - 1, cnt3, cnt4, structure, vc);
                    return;
                  }
                }
              }
      }
    } else {
      for (u = i + TURN + 2; u < j - TURN - 2; u++) {
        int i1u, u1j1;
        i1u   = my_iindx[i + 1] - u;
        u1j1  = my_iindx[u + 1] - j + 1;
        if (!E_M[i1u])
          continue;

        if (!E_M1[u1j1])
          continue;

        /* get distance to reference if closing this multiloop
         *  dist3 = dbp(S_{i,j}, {i,j} + S_{i+1.u} + S_{u+1,j-1})
         */
        d1  = base_d1 - referenceBPs1[i1u] - referenceBPs1[u1j1];
        d2  = base_d2 - referenceBPs2[i1u] - referenceBPs2[u1j1];

        tt      = rtype[type];
        energy  = P->MLclosing;
        if (dangles == 2)
          energy += E_MLstem(tt, S1[j - 1], S1[i + 1], P);
        else
          energy += E_MLstem(tt, -1, -1, P);

        if ((d1 <= k) && (d2 <= l)) {
          for (cnt1 = k_min_M[i1u];
               cnt1 <= MIN2(k - d1, k_max_M[i1u]);
               cnt1++)
            for (cnt2 = l_min_M[i1u][cnt1];
                 cnt2 <= MIN2(l - d2, l_max_M[i1u][cnt1]);
                 cnt2 += 2)
              if (((k - d1 - cnt1) >= k_min_M1[u1j1])
                  && ((k - d1 - cnt1) <= k_max_M1[u1j1])) {
                if (((l - d2 - cnt2) >= l_min_M1[u1j1][k - d1 - cnt1])
                    && ((l - d2 - cnt2) <= l_max_M1[u1j1][k - d1 - cnt1])) {
                  if (e == (energy + E_M[i1u][cnt1][cnt2 / 2] + E_M1[u1j1][k - d1 - cnt1][(l - d2 - cnt2) / 2])) {
                    backtrack_m(i + 1, u, cnt1, cnt2, structure, vc);
                    backtrack_m1(u + 1, j - 1, k - d1 - cnt1, l - d2 - cnt2, structure, vc);
                    return;
                  }
                }
              }
        }
      }
    }
  }

  vrna_message_error("backtracking failed in c");
}


PRIVATE void
backtrack_m(unsigned int          i,
            unsigned int          j,
            int                   k,
            int                   l,
            char                  *structure,
            vrna_fold_compound_t  *vc)
{
  unsigned int  u, ij, seq_length, base_d1, base_d2, d1, d2, maxD1, maxD2;
  int           *my_iindx, *jindx, type, energy, dangles, circ, cnt1, cnt2, cnt3, cnt4;
  int           **l_min_C, **l_max_C, **l_min_M, **l_max_M;
  int           *k_min_C, *k_max_C, *k_min_M, *k_max_M;
  int           ***E_C, ***E_M, *E_C_rem, *E_M_rem;
  short         *S1;
  unsigned int  *referenceBPs1, *referenceBPs2;
  char          *ptype;
  vrna_param_t  *P;
  vrna_md_t     *md;
  vrna_mx_mfe_t *matrices;

  P             = vc->params;
  md            = &(P->model_details);
  matrices      = vc->matrices;
  seq_length    = vc->length;
  S1            = vc->sequence_encoding;
  circ          = md->circ;
  ptype         = vc->ptype;
  my_iindx      = vc->iindx;
  jindx         = vc->jindx;
  referenceBPs1 = vc->referenceBPs1;
  referenceBPs2 = vc->referenceBPs2;
  dangles       = md->dangles;

  E_C     = matrices->E_C;
  l_min_C = matrices->l_min_C;
  l_max_C = matrices->l_max_C;
  k_min_C = matrices->k_min_C;
  k_max_C = matrices->k_max_C;

  E_M     = matrices->E_M;
  l_min_M = matrices->l_min_M;
  l_max_M = matrices->l_max_M;
  k_min_M = matrices->k_min_M;
  k_max_M = matrices->k_max_M;

  E_C_rem = matrices->E_C_rem;
  E_M_rem = matrices->E_M_rem;
  maxD1   = vc->maxD1;
  maxD2   = vc->maxD2;

  ij = my_iindx[i] - j;
  int e = (k == -1) ? E_M_rem[ij] : E_M[ij][k][l / 2];

  base_d1 = referenceBPs1[ij];
  base_d2 = referenceBPs2[ij];

  if (k == -1) {
    /* new_fML = ML(i+1,j)+c */
    d1  = base_d1 - referenceBPs1[my_iindx[i + 1] - j];
    d2  = base_d2 - referenceBPs2[my_iindx[i + 1] - j];
    if (E_M_rem[my_iindx[i + 1] - j] != INF) {
      if (e == (E_M_rem[my_iindx[i + 1] - j] + P->MLbase)) {
        backtrack_m(i + 1, j, -1, -1, structure, vc);
        return;
      }
    }

    if (E_M[my_iindx[i + 1] - j]) {
      for (cnt1 = k_min_M[my_iindx[i + 1] - j];
           cnt1 <= k_max_M[my_iindx[i + 1] - j];
           cnt1++)
        for (cnt2 = l_min_M[my_iindx[i + 1] - j][cnt1];
             cnt2 <= l_max_M[my_iindx[i + 1] - j][cnt1];
             cnt2 += 2)
          if (((cnt1 + d1) > maxD1) || ((cnt2 + d2) > maxD2)) {
            if (e == (E_M[my_iindx[i + 1] - j][cnt1][cnt2 / 2] + P->MLbase)) {
              backtrack_m(i + 1, j, cnt1, cnt2, structure, vc);
              return;
            }
          }
    }

    /* new_fML = min(ML(i,j-1) + c, new_fML) */
    d1  = base_d1 - referenceBPs1[ij + 1];
    d2  = base_d2 - referenceBPs2[ij + 1];
    if (E_M_rem[ij + 1] != INF) {
      if (e == (E_M_rem[ij + 1] + P->MLbase)) {
        backtrack_m(i, j - 1, -1, -1, structure, vc);
        return;
      }
    }

    if (E_M[ij + 1]) {
      for (cnt1 = k_min_M[ij + 1];
           cnt1 <= k_max_M[ij + 1];
           cnt1++)
        for (cnt2 = l_min_M[ij + 1][cnt1];
             cnt2 <= l_max_M[ij + 1][cnt1];
             cnt2 += 2)
          if (((cnt1 + d1) > maxD1) || ((cnt2 + d2) > maxD2)) {
            if (e == (E_M[ij + 1][cnt1][cnt2 / 2] + P->MLbase)) {
              backtrack_m(i, j - 1, cnt1, cnt2, structure, vc);
              return;
            }
          }
    }

    /* new_fML = min(new_fML, C(i,j)+b) */
    if (E_C_rem[ij] != INF) {
      type = ptype[jindx[j] + i];
      if (dangles == 2)
        energy = E_MLstem(type, ((i > 1) || circ) ? S1[i - 1] : -1, ((j < seq_length) || circ) ? S1[j + 1] : -1, P);
      else
        energy = E_MLstem(type, -1, -1, P);

      if (e == (E_C_rem[ij] + energy)) {
        backtrack_c(i, j, -1, -1, structure, vc);
        return;
      }
    }

    /* modular decomposition -------------------------------*/
    for (u = i + 1 + TURN; u <= j - 2 - TURN; u++) {
      int iu, uj;
      iu    = my_iindx[i] - u;
      uj    = my_iindx[u + 1] - j;
      type  = ptype[jindx[j] + u + 1];

      d1  = base_d1 - referenceBPs1[iu] - referenceBPs1[uj];
      d2  = base_d2 - referenceBPs2[iu] - referenceBPs2[uj];

      if (dangles == 2)
        energy = E_MLstem(type, S1[u], (j < seq_length) || circ ? S1[j + 1] : -1, P);
      else
        energy = E_MLstem(type, -1, -1, P);

      if (E_M_rem[iu] != INF) {
        if (E_C[uj]) {
          for (cnt1 = k_min_C[uj];
               cnt1 <= k_max_C[uj];
               cnt1++)
            for (cnt2 = l_min_C[uj][cnt1];
                 cnt2 <= l_max_C[uj][cnt1];
                 cnt2 += 2)
              if (e == (E_M_rem[iu] + E_C[uj][cnt1][cnt2 / 2] + energy)) {
                backtrack_m(i, u, -1, -1, structure, vc);
                backtrack_c(u + 1, j, cnt1, cnt2, structure, vc);
                return;
              }
        }

        if (E_C_rem[uj] != INF) {
          if (e == (E_M_rem[iu] + E_C_rem[uj] + energy)) {
            backtrack_m(i, u, -1, -1, structure, vc);
            backtrack_c(u + 1, j, -1, -1, structure, vc);
            return;
          }
        }
      }

      if (E_C_rem[uj] != INF) {
        if (E_M[iu]) {
          for (cnt1 = k_min_M[iu];
               cnt1 <= k_max_M[iu];
               cnt1++)
            for (cnt2 = l_min_M[iu][cnt1];
                 cnt2 <= l_max_M[iu][cnt1];
                 cnt2 += 2)
              if (e == (E_M[iu][cnt1][cnt2 / 2] + E_C_rem[uj] + energy)) {
                backtrack_m(i, u, cnt1, cnt2, structure, vc);
                backtrack_c(u + 1, j, -1, -1, structure, vc);
                return;
              }
        }
      }

      if (!E_M[iu])
        continue;

      if (!E_C[uj])
        continue;

      for (cnt1 = k_min_M[iu];
           cnt1 <= k_max_M[iu];
           cnt1++)
        for (cnt2 = l_min_M[iu][cnt1];
             cnt2 <= l_max_M[iu][cnt1];
             cnt2 += 2)
          for (cnt3 = k_min_C[uj];
               cnt3 <= k_max_C[uj];
               cnt3++) {
            for (cnt4 = l_min_C[uj][cnt3];
                 cnt4 <= l_max_C[uj][cnt3];
                 cnt4 += 2)
              if (((cnt1 + cnt3 + d1) > maxD1) || ((cnt2 + cnt4 + d2) > maxD2)) {
                if (e == (E_M[iu][cnt1][cnt2 / 2] + E_C[uj][cnt3][cnt4 / 2] + energy)) {
                  backtrack_m(i, u, cnt1, cnt2, structure, vc);
                  backtrack_c(u + 1, j, cnt3, cnt4, structure, vc);
                  return;
                }
              }
          }
    }
  } /* end if (k == -1) */
  else {
    d1  = base_d1 - referenceBPs1[my_iindx[i + 1] - j];
    d2  = base_d2 - referenceBPs2[my_iindx[i + 1] - j];
    /* new_fML = ML(i+1,j)+c */
    if (d1 <= k && d2 <= l) {
      if ((k - d1 >= k_min_M[my_iindx[i + 1] - j]) && (k - d1 <= k_max_M[my_iindx[i + 1] - j])) {
        if ((l - d2 >= l_min_M[my_iindx[i + 1] - j][k - d1]) && (l - d2 <= l_max_M[my_iindx[i + 1] - j][k - d1])) {
          if (E_M[my_iindx[i + 1] - j][k - d1][(l - d2) / 2] + P->MLbase == e) {
            backtrack_m(i + 1, j, k - d1, l - d2, structure, vc);
            return;
          }
        }
      }
    }

    d1  = base_d1 - referenceBPs1[ij + 1];
    d2  = base_d2 - referenceBPs2[ij + 1];

    /* new_fML = min(ML(i,j-1) + c, new_fML) */
    if (E_M[ij + 1]) {
      if (d1 <= k && d2 <= l) {
        if ((k - d1 >= k_min_M[ij + 1]) && (k - d1 <= k_max_M[ij + 1])) {
          if ((l - d2 >= l_min_M[ij + 1][k - d1]) && (l - d2 <= l_max_M[ij + 1][k - d1])) {
            if (E_M[ij + 1][k - d1][(l - d2) / 2] + P->MLbase == e) {
              backtrack_m(i, j - 1, k - d1, l - d2, structure, vc);
              return;
            }
          }
        }
      }
    }

    /* new_fML = min(new_fML, C(i,j)+b) */
    if (E_C[ij]) {
      type = ptype[jindx[j] + i];

      if (dangles == 2)
        energy = E_MLstem(type, ((i > 1) || circ) ? S1[i - 1] : -1, ((j < seq_length) || circ) ? S1[j + 1] : -1, P);
      else
        energy = E_MLstem(type, -1, -1, P);

      if ((k >= k_min_C[ij]) && (k <= k_max_C[ij])) {
        if ((l >= l_min_C[ij][k]) && (l <= l_max_C[ij][k])) {
          if (E_C[ij][k][l / 2] + energy == e) {
            backtrack_c(i, j, k, l, structure, vc);
            return;
          }
        }
      }
    }

    /* modular decomposition -------------------------------*/

    for (u = i + 1 + TURN; u <= j - 2 - TURN; u++) {
      if (!E_M[my_iindx[i] - u])
        continue;

      if (!E_C[my_iindx[u + 1] - j])
        continue;

      type = ptype[jindx[j] + u + 1];

      d1  = base_d1 - referenceBPs1[my_iindx[i] - u] - referenceBPs1[my_iindx[u + 1] - j];
      d2  = base_d2 - referenceBPs2[my_iindx[i] - u] - referenceBPs2[my_iindx[u + 1] - j];

      if (dangles == 2)
        energy = E_MLstem(type, S1[u], ((j < seq_length) || circ) ? S1[j + 1] : -1, P);
      else
        energy = E_MLstem(type, -1, -1, P);

      if (d1 <= k && d2 <= l) {
        for (cnt1 = k_min_M[my_iindx[i] - u]; cnt1 <= MIN2(k - d1, k_max_M[my_iindx[i] - u]); cnt1++)
          for (cnt2 = l_min_M[my_iindx[i] - u][cnt1]; cnt2 <= MIN2(l - d2, l_max_M[my_iindx[i] - u][cnt1]); cnt2 += 2)
            if ((k - d1 - cnt1 >= k_min_C[my_iindx[u + 1] - j]) && (k - d1 - cnt1 <= k_max_C[my_iindx[u + 1] - j])) {
              if ((l - d2 - cnt2 >= l_min_C[my_iindx[u + 1] - j][k - d1 - cnt1]) && (l - d2 - cnt2 <= l_max_C[my_iindx[u + 1] - j][k - d1 - cnt1])) {
                if (E_M[my_iindx[i] - u][cnt1][cnt2 / 2] + E_C[my_iindx[u + 1] - j][k - d1 - cnt1][(l - d2 - cnt2) / 2] + energy == e) {
                  backtrack_m(i, u, cnt1, cnt2, structure, vc);
                  backtrack_c(u + 1, j, k - d1 - cnt1, l - d2 - cnt2, structure, vc);
                  return;
                }
              }
            }
      }
    }
  }

  vrna_message_error("backtracking failed in fML\n");
}


PRIVATE void
backtrack_m1(unsigned int         i,
             unsigned int         j,
             int                  k,
             int                  l,
             char                 *structure,
             vrna_fold_compound_t *vc)
{
  unsigned int  ij, seq_length, d1, d2, *referenceBPs1, *referenceBPs2, maxD1, maxD2;
  int           *my_iindx, *jindx, **l_min_C, **l_max_C, **l_min_M1, **l_max_M1;
  int           *k_min_C, *k_max_C, *k_min_M1, *k_max_M1, cnt1, cnt2;
  int           ***E_C, ***E_M1, *E_C_rem, *E_M1_rem, type, dangles, circ, energy, e_m1;

  short         *S1;
  char          *ptype;
  vrna_param_t  *P;
  vrna_md_t     *md;
  vrna_mx_mfe_t *matrices;

  P             = vc->params;
  md            = &(P->model_details);
  matrices      = vc->matrices;
  seq_length    = vc->length;
  S1            = vc->sequence_encoding;
  ptype         = vc->ptype;
  circ          = md->circ;
  my_iindx      = vc->iindx;
  jindx         = vc->jindx;
  referenceBPs1 = vc->referenceBPs1;
  referenceBPs2 = vc->referenceBPs2;
  dangles       = md->dangles;

  E_C     = matrices->E_C;
  l_min_C = matrices->l_min_C;
  l_max_C = matrices->l_max_C;
  k_min_C = matrices->k_min_C;
  k_max_C = matrices->k_max_C;

  E_M1      = matrices->E_M1;
  l_min_M1  = matrices->l_min_M1;
  l_max_M1  = matrices->l_max_M1;
  k_min_M1  = matrices->k_min_M1;
  k_max_M1  = matrices->k_max_M1;

  E_C_rem   = matrices->E_C_rem;
  E_M1_rem  = matrices->E_M1_rem;
  maxD1     = vc->maxD1;
  maxD2     = vc->maxD2;

  ij    = my_iindx[i] - j;
  e_m1  = (k == -1) ? E_M1_rem[ij] : E_M1[ij][k][l / 2];

  type  = ptype[jindx[j] + i];
  d1    = referenceBPs1[ij] - referenceBPs1[ij + 1];
  d2    = referenceBPs2[ij] - referenceBPs2[ij + 1];

  if (dangles == 2)
    energy = E_MLstem(type, (i > 1) || circ ? S1[i - 1] : -1, (j < seq_length) || circ ? S1[j + 1] : -1, P);
  else
    energy = E_MLstem(type, -1, -1, P);

  if (k == -1) {
    if (E_C_rem[ij] != INF) {
      if (e_m1 == (E_C_rem[ij] + energy)) {
        backtrack_c(i, j, -1, -1, structure, vc);
        return;
      }
    }

    if (E_M1_rem[ij + 1] != INF) {
      if (e_m1 == (E_M1_rem[ij + 1] + P->MLbase)) {
        backtrack_m1(i, j - 1, -1, -1, structure, vc);
        return;
      }
    }

    for (cnt1 = k_min_M1[ij + 1];
         cnt1 <= k_max_M1[ij + 1];
         cnt1++)
      for (cnt2 = l_min_M1[ij + 1][cnt1];
           cnt2 <= l_max_M1[ij + 1][cnt1];
           cnt2 += 2)
        if (((cnt1 + d1) > maxD1) || ((cnt2 + d2) > maxD2)) {
          if (e_m1 == (E_M1[ij + 1][cnt1][cnt2 / 2] + P->MLbase)) {
            backtrack_m1(i, j - 1, cnt1, cnt2, structure, vc);
            return;
          }
        }
  } else {
    if (E_C[ij]) {
      if ((k >= k_min_C[ij]) && (k <= k_max_C[ij])) {
        if ((l >= l_min_C[ij][k]) && (l <= l_max_C[ij][k])) {
          if (E_C[ij][k][l / 2] + energy == e_m1) {
            backtrack_c(i, j, k, l, structure, vc);
            return;
          }
        }
      }
    }

    if (d1 <= k && d2 <= l) {
      if ((k - d1 >= k_min_M1[ij + 1]) && (k - d1 <= k_max_M1[ij + 1])) {
        if ((l - d2 >= l_min_M1[ij + 1][k - d1]) && (l - d2 <= l_max_M1[ij + 1][k - d1])) {
          if (E_M1[ij + 1][k - d1][(l - d2) / 2] + P->MLbase == e_m1) {
            backtrack_m1(i, j - 1, k - d1, l - d2, structure, vc);
            return;
          }
        }
      }
    }
  }

  vrna_message_error("backtack failed in m1\n");
}


PRIVATE void
backtrack_fc(int                  k,
             int                  l,
             char                 *structure,
             vrna_fold_compound_t *vc)
{
  unsigned int  d, i, j, seq_length, base_d1, base_d2, d1, d2, maxD1, maxD2;
  int           *my_iindx, *jindx, energy, cnt1, cnt2, cnt3, cnt4, *rtype;
  short         *S1;
  unsigned int  *referenceBPs1, *referenceBPs2;
  char          *sequence, *ptype;
  int           **E_Fc, **E_FcH, **E_FcI, **E_FcM, ***E_C, ***E_M, ***E_M2;
  int           *E_C_rem, *E_M_rem, *E_M2_rem, E_Fc_rem, E_FcH_rem, E_FcI_rem, E_FcM_rem;
  int           **l_min_C, **l_max_C, *k_min_C, *k_max_C;
  int           **l_min_M, **l_max_M, *k_min_M, *k_max_M;
  int           **l_min_M2, **l_max_M2, *k_min_M2, *k_max_M2;
  int           *l_min_FcH, *l_max_FcH, k_min_FcH, k_max_FcH;
  int           *l_min_FcI, *l_max_FcI, k_min_FcI, k_max_FcI;
  int           *l_min_FcM, *l_max_FcM, k_min_FcM, k_max_FcM;
  vrna_param_t  *P;
  vrna_md_t     *md;
  vrna_mx_mfe_t *matrices;

  P             = vc->params;
  md            = &(P->model_details);
  matrices      = vc->matrices;
  sequence      = vc->sequence;
  seq_length    = vc->length;
  S1            = vc->sequence_encoding;
  ptype         = vc->ptype;
  rtype         = &(md->rtype[0]);
  my_iindx      = vc->iindx;
  jindx         = vc->jindx;
  referenceBPs1 = vc->referenceBPs1;
  referenceBPs2 = vc->referenceBPs2;

  base_d1 = referenceBPs1[my_iindx[1] - seq_length];
  base_d2 = referenceBPs2[my_iindx[1] - seq_length];

  E_C     = matrices->E_C;
  l_min_C = matrices->l_min_C;
  l_max_C = matrices->l_max_C;
  k_min_C = matrices->k_min_C;
  k_max_C = matrices->k_max_C;

  E_M     = matrices->E_M;
  l_min_M = matrices->l_min_M;
  l_max_M = matrices->l_max_M;
  k_min_M = matrices->k_min_M;
  k_max_M = matrices->k_max_M;

  E_M2      = matrices->E_M2;
  l_min_M2  = matrices->l_min_M2;
  l_max_M2  = matrices->l_max_M2;
  k_min_M2  = matrices->k_min_M2;
  k_max_M2  = matrices->k_max_M2;

  E_Fc = matrices->E_Fc;

  E_FcI     = matrices->E_FcI;
  l_min_FcI = matrices->l_min_FcI;
  l_max_FcI = matrices->l_max_FcI;
  k_min_FcI = matrices->k_min_FcI;
  k_max_FcI = matrices->k_max_FcI;

  E_FcH     = matrices->E_FcH;
  l_min_FcH = matrices->l_min_FcH;
  l_max_FcH = matrices->l_max_FcH;
  k_min_FcH = matrices->k_min_FcH;
  k_max_FcH = matrices->k_max_FcH;

  E_FcM     = matrices->E_FcM;
  l_min_FcM = matrices->l_min_FcM;
  l_max_FcM = matrices->l_max_FcM;
  k_min_FcM = matrices->k_min_FcM;
  k_max_FcM = matrices->k_max_FcM;

  E_C_rem   = matrices->E_C_rem;
  E_M_rem   = matrices->E_M_rem;
  E_M2_rem  = matrices->E_M2_rem;
  E_Fc_rem  = matrices->E_Fc_rem;
  E_FcH_rem = matrices->E_FcH_rem;
  E_FcI_rem = matrices->E_FcI_rem;
  E_FcM_rem = matrices->E_FcM_rem;
  maxD1     = vc->maxD1;
  maxD2     = vc->maxD2;


  if (k == -1) {
    /* check if mfe might be open chain */
    if (E_Fc_rem == 0)
      if ((referenceBPs1[my_iindx[1] - seq_length] > maxD1) || (referenceBPs2[my_iindx[1] - seq_length] > maxD2))
        return;

    /* check for hairpin configurations */
    if (E_Fc_rem == E_FcH_rem) {
      for (d = TURN + 2; d <= seq_length; d++) /* i,j in [1..length] */
        for (j = d; j <= seq_length; j++) {
          unsigned int  u, ij;
          int           type, no_close;
          char          loopseq[10];
          i   = j - d + 1;
          ij  = my_iindx[i] - j;
          u   = seq_length - j + i - 1;
          if (u < TURN)
            continue;

          type      = ptype[jindx[j] + i];
          no_close  = (((type == 3) || (type == 4)) && no_closingGU);
          type      = rtype[type];
          if (!type)
            continue;

          if (no_close)
            continue;

          d1  = base_d1 - referenceBPs1[ij];
          d2  = base_d2 - referenceBPs2[ij];
          if (u < 7) {
            strcpy(loopseq, sequence + j - 1);
            strncat(loopseq, sequence, i);
          }

          energy = E_Hairpin(u, type, S1[j + 1], S1[i - 1], loopseq, P);

          if (E_C_rem[ij] != INF) {
            if (E_Fc_rem == (E_C_rem[ij] + energy)) {
              backtrack_c(i, j, -1, -1, structure, vc);
              return;
            }
          }

          if (E_C[ij]) {
            for (cnt1 = k_min_C[ij];
                 cnt1 <= k_max_C[ij];
                 cnt1++)
              for (cnt2 = l_min_C[ij][cnt1];
                   cnt2 <= l_max_C[ij][cnt1];
                   cnt2 += 2)
                if (((cnt1 + d1) > maxD1) || ((cnt2 + d2) > maxD2)) {
                  if (E_Fc_rem == (E_C[ij][cnt1][cnt2 / 2] + energy)) {
                    backtrack_c(i, j, cnt1, cnt2, structure, vc);
                    return;
                  }
                }
          }
        }
    }

    /* check for interior loop configurations */
    if (E_Fc_rem == E_FcI_rem) {
      for (d = TURN + 2; d <= seq_length; d++) /* i,j in [1..length] */
        for (j = d; j <= seq_length; j++) {
          unsigned int  u, ij, p, q, pq;
          int           type, type_2;
          i   = j - d + 1;
          ij  = my_iindx[i] - j;
          u   = seq_length - j + i - 1;
          if (u < TURN)
            continue;

          type = rtype[(unsigned int)ptype[jindx[j] + i]];
          if (!type)
            continue;

          for (p = j + 1; p < seq_length; p++) {
            unsigned int u1, qmin, ln_pre;
            u1 = p - j - 1;
            if (u1 + i - 1 > MAXLOOP)
              break;

            qmin    = p + TURN + 1;
            ln_pre  = u1 + i + seq_length;
            if (ln_pre > qmin + MAXLOOP)
              qmin = ln_pre - MAXLOOP - 1;

            for (q = qmin; q <= seq_length; q++) {
              unsigned int u2;
              pq      = my_iindx[p] - q;
              type_2  = rtype[(unsigned int)ptype[jindx[q] + p]];
              if (type_2 == 0)
                continue;

              u2 = i - 1 + seq_length - q;
              if (u1 + u2 > MAXLOOP)
                continue;

              energy = E_IntLoop(u1, u2, type, type_2, S1[j + 1], S1[i - 1], S1[p - 1], S1[q + 1], P);
              if (E_C_rem[ij] != INF) {
                if (E_C[pq]) {
                  for (cnt1 = k_min_C[pq];
                       cnt1 <= k_max_C[pq];
                       cnt1++)
                    for (cnt2 = l_min_C[pq][cnt1];
                         cnt2 <= l_max_C[pq][cnt1];
                         cnt2 += 2)
                      if (E_Fc_rem == (E_C_rem[ij] + E_C[pq][cnt1][cnt2 / 2] + energy)) {
                        backtrack_c(i, j, -1, -1, structure, vc);
                        backtrack_c(p, q, cnt1, cnt2, structure, vc);
                        return;
                      }
                }

                if (E_C_rem[pq] != INF) {
                  if (E_Fc_rem == (E_C_rem[ij] + E_C_rem[pq] + energy)) {
                    backtrack_c(i, j, -1, -1, structure, vc);
                    backtrack_c(p, q, -1, -1, structure, vc);
                    return;
                  }
                }
              }

              if (E_C_rem[pq] != INF) {
                if (E_C[ij]) {
                  for (cnt1 = k_min_C[ij];
                       cnt1 <= k_max_C[ij];
                       cnt1++)
                    for (cnt2 = l_min_C[ij][cnt1];
                         cnt2 <= l_max_C[ij][cnt1];
                         cnt2 += 2)
                      if (E_Fc_rem == (E_C[ij][cnt1][cnt2 / 2] + E_C_rem[pq] + energy)) {
                        backtrack_c(i, j, cnt1, cnt2, structure, vc);
                        backtrack_c(p, q, -1, -1, structure, vc);
                        return;
                      }
                }
              }

              if (!(E_C[ij]))
                continue;

              if (!(E_C[pq]))
                continue;

              /* get distance to reference if closing the interior loop
               *  d2a = dbp(T1_[1,n}, T1_{p,q} + T1_{i,j})
               *  d2b = dbp(T2_[1,n}, T2_{p,q} + T2_{i,j})
               */
              d1  = base_d1 - referenceBPs1[ij] - referenceBPs1[pq];
              d2  = base_d2 - referenceBPs2[ij] - referenceBPs2[pq];
              for (cnt1 = k_min_C[ij];
                   cnt1 <= k_max_C[ij];
                   cnt1++)
                for (cnt2 = l_min_C[ij][cnt1];
                     cnt2 <= l_max_C[ij][cnt1];
                     cnt2 += 2)
                  for (cnt3 = k_min_C[pq];
                       cnt3 <= k_max_C[pq];
                       cnt3++)
                    for (cnt4 = l_min_C[pq][cnt3];
                         cnt4 <= l_max_C[pq][cnt3];
                         cnt4 += 2)
                      if (((cnt1 + cnt3 + d1) > maxD1) || ((cnt2 + cnt4 + d2) > maxD2)) {
                        if (E_Fc_rem == (E_C[ij][cnt1][cnt2 / 2] + E_C[pq][cnt3][cnt4 / 2] + energy)) {
                          backtrack_c(i, j, cnt1, cnt2, structure, vc);
                          backtrack_c(p, q, cnt3, cnt4, structure, vc);
                          return;
                        }
                      }
            } /* end for p */
          }   /* end for q */
        }
    }

    /* check for multi loop configurations */
    if (E_Fc_rem == E_FcM_rem) {
      if (seq_length > 2 * TURN) {
        for (i = TURN + 1; i < seq_length - 2 * TURN; i++) {
          /* get distancies to references
           * d3a = dbp(T1_[1,n}, T1_{1,k} + T1_{k+1, n})
           * d3b = dbp(T2_[1,n}, T2_{1,k} + T2_{k+1, n})
           */
          if (E_M_rem[my_iindx[1] - i] != INF) {
            if (E_M2[i + 1]) {
              for (cnt1 = k_min_M2[i + 1];
                   cnt1 <= k_max_M2[i + 1];
                   cnt1++)
                for (cnt2 = l_min_M2[i + 1][cnt1];
                     cnt2 <= l_max_M2[i + 1][cnt1];
                     cnt2 += 2)
                  if (E_Fc_rem == (E_M_rem[my_iindx[1] - i] + E_M2[i + 1][cnt1][cnt2 / 2] + P->MLclosing)) {
                    backtrack_m(1, i, -1, -1, structure, vc);
                    backtrack_m2(i + 1, cnt1, cnt2, structure, vc);
                    return;
                  }
            }

            if (E_M2_rem[i + 1] != INF) {
              if (E_Fc_rem == (E_M_rem[my_iindx[1] - i] + E_M2_rem[i + 1] + P->MLclosing)) {
                backtrack_m(1, i, -1, -1, structure, vc);
                backtrack_m2(i + 1, -1, -1, structure, vc);
                return;
              }
            }
          }

          if (E_M2_rem[i + 1] != INF) {
            if (E_M[my_iindx[1] - i]) {
              for (cnt1 = k_min_M[my_iindx[1] - i];
                   cnt1 <= k_max_M[my_iindx[1] - i];
                   cnt1++)
                for (cnt2 = l_min_M[my_iindx[1] - i][cnt1];
                     cnt2 <= l_max_M[my_iindx[1] - i][cnt1];
                     cnt2 += 2)
                  if (E_Fc_rem == (E_M[my_iindx[1] - i][cnt1][cnt2 / 2] + E_M2_rem[i + 1] + P->MLclosing)) {
                    backtrack_m(1, i, cnt1, cnt2, structure, vc);
                    backtrack_m2(i + 1, -1, -1, structure, vc);
                    return;
                  }
            }
          }

          if (!(E_M[my_iindx[1] - i]))
            continue;

          if (!(E_M2[i + 1]))
            continue;

          d1  = base_d1 - referenceBPs1[my_iindx[1] - i] - referenceBPs1[my_iindx[i + 1] - seq_length];
          d2  = base_d2 - referenceBPs2[my_iindx[1] - i] - referenceBPs2[my_iindx[i + 1] - seq_length];
          for (cnt1 = k_min_M[my_iindx[1] - i];
               cnt1 <= k_max_M[my_iindx[1] - i];
               cnt1++)
            for (cnt2 = l_min_M[my_iindx[1] - i][cnt1];
                 cnt2 <= l_max_M[my_iindx[1] - i][cnt1];
                 cnt2 += 2)
              for (cnt3 = k_min_M2[i + 1];
                   cnt3 <= k_max_M2[i + 1];
                   cnt3++)
                for (cnt4 = l_min_M2[i + 1][cnt3];
                     cnt4 <= l_max_M2[i + 1][cnt3];
                     cnt4 += 2)
                  if (((cnt1 + cnt3 + d1) > maxD1) || ((cnt2 + cnt4 + d2) > maxD2)) {
                    if (E_Fc_rem == (E_M[my_iindx[1] - i][cnt1][cnt2 / 2] + E_M2[i + 1][cnt3][cnt4 / 2] + P->MLclosing)) {
                      backtrack_m(1, i, cnt1, cnt2, structure, vc);
                      backtrack_m2(i + 1, cnt3, cnt4, structure, vc);
                      return;
                    }
                  }
        }
      }
    }
  } else {
    /* open chain ? */
    if (E_Fc[k][l / 2] == 0)
      if ((k == referenceBPs1[my_iindx[1] - seq_length]) && (l == referenceBPs2[my_iindx[1] - seq_length]))
        return;

    if ((k >= k_min_FcH) && (k <= k_max_FcH)) {
      if ((l >= l_min_FcH[k]) && (l <= l_max_FcH[k])) {
        if (E_Fc[k][l / 2] == E_FcH[k][l / 2]) {
          for (d = TURN + 2; d <= seq_length; d++) /* i,j in [1..length] */
            for (j = d; j <= seq_length; j++) {
              unsigned int  u, ij;
              int           type, no_close;
              char          loopseq[10];
              i   = j - d + 1;
              ij  = my_iindx[i] - j;
              if (!E_C[ij])
                continue;

              u = seq_length - j + i - 1;
              if (u < TURN)
                continue;

              type = ptype[jindx[j] + i];

              no_close = (((type == 3) || (type == 4)) && no_closingGU);

              type = rtype[type];

              if (!type)
                continue;

              if (no_close)
                continue;

              d1  = base_d1 - referenceBPs1[ij];
              d2  = base_d2 - referenceBPs2[ij];
              if (u < 7) {
                strcpy(loopseq, sequence + j - 1);
                strncat(loopseq, sequence, i);
              }

              energy = E_Hairpin(u, type, S1[j + 1], S1[i - 1], loopseq, P);
              if ((k >= d1) && (l >= d2)) {
                if ((k - d1 >= k_min_C[ij]) && (k - d1 <= k_max_C[ij])) {
                  if ((l - d2 >= l_min_C[ij][k - d1]) && (l - d2 <= l_max_C[ij][k - d1])) {
                    if (E_Fc[k][l / 2] == E_C[ij][k - d1][(l - d2) / 2] + energy) {
                      backtrack_c(i, j, k - d1, l - d2, structure, vc);
                      return;
                    }
                  }
                }
              }
            }
        }
      }
    }

    if ((k >= k_min_FcI) && (k <= k_max_FcI)) {
      if ((l >= l_min_FcI[k]) && (l <= l_max_FcI[k])) {
        if (E_Fc[k][l / 2] == E_FcI[k][l / 2]) {
          for (d = TURN + 2; d <= seq_length; d++) /* i,j in [1..length] */
            for (j = d; j <= seq_length; j++) {
              unsigned int  u, ij, p, q, pq;
              int           type, type_2;
              i   = j - d + 1;
              ij  = my_iindx[i] - j;
              if (!E_C[ij])
                continue;

              u = seq_length - j + i - 1;
              if (u < TURN)
                continue;

              type = ptype[jindx[j] + i];

              type = rtype[type];

              if (!type)
                continue;

              for (p = j + 1; p < seq_length; p++) {
                unsigned int u1, qmin, ln_pre;
                u1 = p - j - 1;
                if (u1 + i - 1 > MAXLOOP)
                  break;

                qmin    = p + TURN + 1;
                ln_pre  = u1 + i + seq_length;
                if (ln_pre > qmin + MAXLOOP)
                  qmin = ln_pre - MAXLOOP - 1;

                for (q = qmin; q <= seq_length; q++) {
                  unsigned int u2;
                  pq = my_iindx[p] - q;
                  if (!E_C[pq])
                    continue;

                  type_2 = rtype[(unsigned int)ptype[jindx[q] + p]];
                  if (type_2 == 0)
                    continue;

                  u2 = i - 1 + seq_length - q;
                  if (u1 + u2 > MAXLOOP)
                    continue;

                  /* get distance to reference if closing the interior loop
                   *  d2a = dbp(T1_[1,n}, T1_{p,q} + T1_{i,j})
                   *  d2b = dbp(T2_[1,n}, T2_{p,q} + T2_{i,j})
                   */
                  d1      = base_d1 - referenceBPs1[ij] - referenceBPs1[pq];
                  d2      = base_d2 - referenceBPs2[ij] - referenceBPs2[pq];
                  energy  = E_IntLoop(u1, u2, type, type_2, S1[j + 1], S1[i - 1], S1[p - 1], S1[q + 1], P);
                  if ((k >= d1) && (l >= d2)) {
                    for (cnt1 = k_min_C[ij]; cnt1 <= MIN2(k_max_C[ij], k - d1); cnt1++)
                      for (cnt2 = l_min_C[ij][cnt1]; cnt2 <= MIN2(l_max_C[ij][cnt1], l - d2); cnt2 += 2)
                        if ((k - d1 - cnt1 >= k_min_C[pq]) && (k - d1 - cnt1 <= k_max_C[pq])) {
                          if ((l - d2 - cnt2 >= l_min_C[pq][k - d1 - cnt1]) && (l - d2 - cnt2 <= l_max_C[pq][k - d1 - cnt1])) {
                            if ((E_C[ij][cnt1][cnt2 / 2] + E_C[pq][k - d1 - cnt1][(l - d2 - cnt2) / 2] + energy) == E_Fc[k][l / 2]) {
                              backtrack_c(i, j, cnt1, cnt2, structure, vc);
                              backtrack_c(p, q, k - d1 - cnt1, l - d2 - cnt2, structure, vc);
                              return;
                            }
                          }
                        }
                  }
                }
              }
            }
        }
      }
    }

    if ((k >= k_min_FcM) && (k <= k_max_FcM)) {
      if ((l >= l_min_FcM[k]) && (l <= l_max_FcM[k])) {
        if (E_Fc[k][l / 2] == E_FcM[k][l / 2]) {
          if (seq_length > 2 * TURN) {
            for (i = TURN + 1; i < seq_length - 2 * TURN; i++) {
              /* get distancies to references
               * d3a = dbp(T1_[1,n}, T1_{1,k} + T1_{k+1, n})
               * d3b = dbp(T2_[1,n}, T2_{1,k} + T2_{k+1, n})
               */
              if (!E_M[my_iindx[1] - i])
                continue;

              if (!E_M2[i + 1])
                continue;

              d1  = base_d1 - referenceBPs1[my_iindx[1] - i] - referenceBPs1[my_iindx[i + 1] - seq_length];
              d2  = base_d2 - referenceBPs2[my_iindx[1] - i] - referenceBPs2[my_iindx[i + 1] - seq_length];
              if ((k >= d1) && (l >= d2)) {
                for (cnt1 = k_min_M[my_iindx[1] - i]; cnt1 <= MIN2(k_max_M[my_iindx[1] - i], k - d1); cnt1++)
                  for (cnt2 = l_min_M[my_iindx[1] - i][cnt1]; cnt2 <= MIN2(l_max_M[my_iindx[1] - i][cnt1], l - d2); cnt2 += 2)
                    if ((k - d1 - cnt1 >= k_min_M2[i + 1]) && (k - d1 - cnt1 <= k_max_M2[i + 1])) {
                      if ((l - d2 - cnt2 >= l_min_M2[i + 1][k - d1 - cnt1]) && (l - d2 - cnt2 <= l_max_M2[i + 1][k - d1 - cnt1])) {
                        if ((E_M[my_iindx[1] - i][cnt1][cnt2 / 2] + E_M2[i + 1][k - d1 - cnt1][(l - d2 - cnt2) / 2] + P->MLclosing) == E_FcM[k][l / 2]) {
                          backtrack_m(1, i, cnt1, cnt2, structure, vc);
                          backtrack_m2(i + 1, k - d1 - cnt1, l - d2 - cnt2, structure, vc);
                          return;
                        }
                      }
                    }
              }
            }
          }
        }
      }
    }
  }

  vrna_message_error("backtack failed in fc\n");
}


PRIVATE void
backtrack_m2(unsigned int         i,
             int                  k,
             int                  l,
             char                 *structure,
             vrna_fold_compound_t *vc)
{
  unsigned int  j, ij, j3, n;
  unsigned int  *referenceBPs1, *referenceBPs2;
  unsigned int  d1, d2, base_d1, base_d2, maxD1, maxD2;
  int           *my_iindx, cnt1, cnt2, cnt3, cnt4;
  int           ***E_M1, ***E_M2, *E_M2_rem, *E_M1_rem, e;
  int           **l_min_M1, **l_max_M1, *k_min_M1, *k_max_M1;
  vrna_mx_mfe_t *matrices;

  matrices      = vc->matrices;
  n             = vc->length;
  my_iindx      = vc->iindx;
  referenceBPs1 = vc->referenceBPs1;
  referenceBPs2 = vc->referenceBPs2;

  E_M1      = matrices->E_M1;
  l_min_M1  = matrices->l_min_M1;
  l_max_M1  = matrices->l_max_M1;
  k_min_M1  = matrices->k_min_M1;
  k_max_M1  = matrices->k_max_M1;

  E_M1_rem = matrices->E_M1_rem;

  E_M2 = matrices->E_M2;

  E_M2_rem = matrices->E_M2_rem;

  maxD1 = vc->maxD1;
  maxD2 = vc->maxD2;

  base_d1 = referenceBPs1[my_iindx[i] - n];
  base_d2 = referenceBPs2[my_iindx[i] - n];

  if (k == -1) {
    e = E_M2_rem[i];
    for (j = i + TURN + 1; j < n - TURN - 1; j++) {
      if (E_M1_rem[my_iindx[i] - j] != INF) {
        if (E_M1[my_iindx[j + 1] - n]) {
          for (cnt1 = k_min_M1[my_iindx[j + 1] - n];
               cnt1 <= k_max_M1[my_iindx[j + 1] - n];
               cnt1++)
            for (cnt2 = l_min_M1[my_iindx[j + 1] - n][cnt1];
                 cnt2 <= l_max_M1[my_iindx[j + 1] - n][cnt1];
                 cnt2++)
              if (e == E_M1_rem[my_iindx[i] - j] + E_M1[my_iindx[j + 1] - n][cnt1][cnt2 / 2]) {
                backtrack_m1(i, j, k, l, structure, vc);
                backtrack_m1(j + 1, n, cnt1, cnt2, structure, vc);
                return;
              }
        }

        if (E_M1_rem[my_iindx[j + 1] - n] != INF) {
          if (e == E_M1_rem[my_iindx[i] - j] + E_M1_rem[my_iindx[j + 1] - n]) {
            backtrack_m1(i, j, k, l, structure, vc);
            backtrack_m1(j + 1, n, k, l, structure, vc);
            return;
          }
        }
      }

      if (E_M1_rem[my_iindx[j + 1] - n] != INF) {
        if (E_M1[my_iindx[i] - j]) {
          for (cnt1 = k_min_M1[my_iindx[i] - j];
               cnt1 <= k_max_M1[my_iindx[i] - j];
               cnt1++)
            for (cnt2 = l_min_M1[my_iindx[i] - j][cnt1];
                 cnt2 <= l_max_M1[my_iindx[i] - j][cnt1];
                 cnt2 += 2)
              if (e == E_M1[my_iindx[i] - j][cnt1][cnt2 / 2] + E_M1_rem[my_iindx[j + 1] - n]) {
                backtrack_m1(i, j, cnt1, cnt2, structure, vc);
                backtrack_m1(j + 1, n, k, l, structure, vc);
                return;
              }
        }
      }

      if (!E_M1[my_iindx[i] - j])
        continue;

      if (!E_M1[my_iindx[j + 1] - n])
        continue;

      d1  = referenceBPs1[my_iindx[i] - n] - referenceBPs1[my_iindx[i] - j] - referenceBPs1[my_iindx[j + 1] - n];
      d2  = referenceBPs2[my_iindx[i] - n] - referenceBPs2[my_iindx[i] - j] - referenceBPs2[my_iindx[j + 1] - n];

      for (cnt1 = k_min_M1[my_iindx[i] - j]; cnt1 <= k_max_M1[my_iindx[i] - j]; cnt1++)
        for (cnt2 = l_min_M1[my_iindx[i] - j][cnt1]; cnt2 <= l_max_M1[my_iindx[i] - j][cnt1]; cnt2 += 2) {
          for (cnt3 = k_min_M1[my_iindx[j + 1] - n]; cnt3 <= k_max_M1[my_iindx[j + 1] - n]; cnt3++)
            for (cnt4 = l_min_M1[my_iindx[j + 1] - n][cnt3]; cnt4 <= l_max_M1[my_iindx[j + 1] - n][cnt3]; cnt4 += 2) {
              if (((cnt1 + cnt3 + d1) > maxD1) || ((cnt2 + cnt4 + d2) > maxD2)) {
                if (e == E_M1[my_iindx[i] - j][cnt1][cnt2 / 2] + E_M1[my_iindx[j + 1] - n][cnt3][cnt4 / 2]) {
                  backtrack_m1(i, j, cnt1, cnt2, structure, vc);
                  backtrack_m1(j + 1, n, cnt3, cnt4, structure, vc);
                  return;
                }
              }
            }
        }
    }
  } else {
    for (j = i + TURN + 1; j < n - TURN - 1; j++) {
      if (!E_M1[my_iindx[i] - j])
        continue;

      if (!E_M1[my_iindx[j + 1] - n])
        continue;

      ij  = my_iindx[i] - j;
      j3  = my_iindx[j + 1] - n;
      d1  = base_d1 - referenceBPs1[ij] - referenceBPs1[j3];
      d2  = base_d2 - referenceBPs2[ij] - referenceBPs2[j3];

      for (cnt1 = k_min_M1[ij]; cnt1 <= MIN2(k_max_M1[ij], k - d1); cnt1++)
        for (cnt2 = l_min_M1[ij][cnt1]; cnt2 <= MIN2(l_max_M1[ij][cnt1], l - d2); cnt2 += 2)
          if ((k - d1 - cnt1 >= k_min_M1[j3]) && (k - d1 - cnt1 <= k_max_M1[j3])) {
            if ((l - d2 - cnt2 >= l_min_M1[j3][k - d1 - cnt1]) && (l - d2 - cnt2 <= l_max_M1[j3][k - d1 - cnt1])) {
              if (E_M1[ij][cnt1][cnt2 / 2] + E_M1[j3][k - d1 - cnt1][(l - d2 - cnt2) / 2] == E_M2[i][k][l / 2]) {
                backtrack_m1(i, j, cnt1, cnt2, structure, vc);
                backtrack_m1(j + 1, n, k - d1 - cnt1, l - d2 - cnt2, structure, vc);
                return;
              }
            }
          }
    }
  }

  vrna_message_error("backtack failed in m2\n");
}


PRIVATE void
mfe_circ(vrna_fold_compound_t *vc)
{
  unsigned int  d, i, j, maxD1, maxD2, seq_length, *referenceBPs1, *referenceBPs2, d1, d2, base_d1, base_d2, *mm1, *mm2, *bpdist;
  int           *my_iindx, *jindx, energy, cnt1, cnt2, cnt3, cnt4, *rtype;
  short         *S1;
  char          *sequence, *ptype;
  int           ***E_C, ***E_M, ***E_M1;
  int           *E_C_rem, *E_M_rem, *E_M1_rem;
  int           **l_min_C, **l_max_C, **l_min_M, **l_max_M, **l_min_M1, **l_max_M1;
  int           *k_min_C, *k_max_C, *k_min_M, *k_max_M, *k_min_M1, *k_max_M1;

  vrna_param_t  *P;
  vrna_md_t     *md;
  vrna_mx_mfe_t *matrices;

  P             = vc->params;
  md            = &(P->model_details);
  matrices      = vc->matrices;
  sequence      = vc->sequence;
  seq_length    = vc->length;
  maxD1         = vc->maxD1;
  maxD2         = vc->maxD2;
  S1            = vc->sequence_encoding;
  ptype         = vc->ptype;
  rtype         = &(md->rtype[0]);
  my_iindx      = vc->iindx;
  jindx         = vc->jindx;
  referenceBPs1 = vc->referenceBPs1;
  referenceBPs2 = vc->referenceBPs2;
  mm1           = vc->mm1;
  mm2           = vc->mm2;
  bpdist        = vc->bpdist;

  E_C     = matrices->E_C;
  l_min_C = matrices->l_min_C;
  l_max_C = matrices->l_max_C;
  k_min_C = matrices->k_min_C;
  k_max_C = matrices->k_max_C;

  E_M     = matrices->E_M;
  l_min_M = matrices->l_min_M;
  l_max_M = matrices->l_max_M;
  k_min_M = matrices->k_min_M;
  k_max_M = matrices->k_max_M;

  E_M1      = matrices->E_M1;
  l_min_M1  = matrices->l_min_M1;
  l_max_M1  = matrices->l_max_M1;
  k_min_M1  = matrices->k_min_M1;
  k_max_M1  = matrices->k_max_M1;

  E_C_rem   = matrices->E_C_rem;
  E_M_rem   = matrices->E_M_rem;
  E_M1_rem  = matrices->E_M1_rem;

#ifdef _OPENMP
#pragma omp parallel for private(d1,d2,cnt1,cnt2,cnt3,cnt4,j, i)
#endif
  for (i = 1; i < seq_length - TURN - 1; i++) {
    /* guess memory requirements for M2 */

    int min_k, max_k, max_l, min_l;
    int min_k_real, max_k_real, *min_l_real, *max_l_real;

    min_k = min_l = 0;
    max_k = mm1[my_iindx[i] - seq_length] + referenceBPs1[my_iindx[i] - seq_length];
    max_l = mm2[my_iindx[i] - seq_length] + referenceBPs2[my_iindx[i] - seq_length];

    prepareBoundaries(min_k,
                      max_k,
                      min_l,
                      max_l,
                      bpdist[my_iindx[i] - seq_length],
                      &matrices->k_min_M2[i],
                      &matrices->k_max_M2[i],
                      &matrices->l_min_M2[i],
                      &matrices->l_max_M2[i]
                      );

    prepareArray(&matrices->E_M2[i],
                 matrices->k_min_M2[i],
                 matrices->k_max_M2[i],
                 matrices->l_min_M2[i],
                 matrices->l_max_M2[i]
                 );

    preparePosteriorBoundaries(matrices->k_max_M2[i] - matrices->k_min_M2[i] + 1,
                               matrices->k_min_M2[i],
                               &min_k_real,
                               &max_k_real,
                               &min_l_real,
                               &max_l_real
                               );

    /* begin filling of M2 array */
    for (j = i + TURN + 1; j < seq_length - TURN - 1; j++) {
      if (E_M1_rem[my_iindx[i] - j] != INF) {
        if (E_M1[my_iindx[j + 1] - seq_length]) {
          for (cnt1 = k_min_M1[my_iindx[j + 1] - seq_length];
               cnt1 <= k_max_M1[my_iindx[j + 1] - seq_length];
               cnt1++)
            for (cnt2 = l_min_M1[my_iindx[j + 1] - seq_length][cnt1];
                 cnt2 <= l_max_M1[my_iindx[j + 1] - seq_length][cnt1];
                 cnt2++)
              matrices->E_M2_rem[i] = MIN2(matrices->E_M2_rem[i],
                                           E_M1_rem[my_iindx[i] - j] + E_M1[my_iindx[j + 1] - seq_length][cnt1][cnt2 / 2]
                                           );
        }

        if (E_M1_rem[my_iindx[j + 1] - seq_length] != INF)
          matrices->E_M2_rem[i] = MIN2(matrices->E_M2_rem[i], E_M1_rem[my_iindx[i] - j] + E_M1_rem[my_iindx[j + 1] - seq_length]);
      }

      if (E_M1_rem[my_iindx[j + 1] - seq_length] != INF) {
        if (E_M1[my_iindx[i] - j]) {
          for (cnt1 = k_min_M1[my_iindx[i] - j];
               cnt1 <= k_max_M1[my_iindx[i] - j];
               cnt1++)
            for (cnt2 = l_min_M1[my_iindx[i] - j][cnt1];
                 cnt2 <= l_max_M1[my_iindx[i] - j][cnt1];
                 cnt2 += 2)
              matrices->E_M2_rem[i] = MIN2(matrices->E_M2_rem[i],
                                           E_M1[my_iindx[i] - j][cnt1][cnt2 / 2] + E_M1_rem[my_iindx[j + 1] - seq_length]
                                           );
        }
      }

      if (!E_M1[my_iindx[i] - j])
        continue;

      if (!E_M1[my_iindx[j + 1] - seq_length])
        continue;

      d1  = referenceBPs1[my_iindx[i] - seq_length] - referenceBPs1[my_iindx[i] - j] - referenceBPs1[my_iindx[j + 1] - seq_length];
      d2  = referenceBPs2[my_iindx[i] - seq_length] - referenceBPs2[my_iindx[i] - j] - referenceBPs2[my_iindx[j + 1] - seq_length];

      for (cnt1 = k_min_M1[my_iindx[i] - j]; cnt1 <= k_max_M1[my_iindx[i] - j]; cnt1++)
        for (cnt2 = l_min_M1[my_iindx[i] - j][cnt1]; cnt2 <= l_max_M1[my_iindx[i] - j][cnt1]; cnt2 += 2) {
          for (cnt3 = k_min_M1[my_iindx[j + 1] - seq_length]; cnt3 <= k_max_M1[my_iindx[j + 1] - seq_length]; cnt3++)
            for (cnt4 = l_min_M1[my_iindx[j + 1] - seq_length][cnt3]; cnt4 <= l_max_M1[my_iindx[j + 1] - seq_length][cnt3]; cnt4 += 2) {
              if (((cnt1 + cnt3 + d1) <= maxD1) && ((cnt2 + cnt4 + d2) <= maxD2)) {
                matrices->E_M2[i][cnt1 + cnt3 + d1][(cnt2 + cnt4 + d2) / 2] = MIN2(matrices->E_M2[i][cnt1 + cnt3 + d1][(cnt2 + cnt4 + d2) / 2],
                                                                                   E_M1[my_iindx[i] - j][cnt1][cnt2 / 2] + E_M1[my_iindx[j + 1] - seq_length][cnt3][cnt4 / 2]
                                                                                   );
                updatePosteriorBoundaries(cnt1 + cnt3 + d1,
                                          cnt2 + cnt4 + d2,
                                          &min_k_real,
                                          &max_k_real,
                                          &min_l_real,
                                          &max_l_real
                                          );
              } else {
                matrices->E_M2_rem[i] = MIN2(matrices->E_M2_rem[i],
                                             E_M1[my_iindx[i] - j][cnt1][cnt2 / 2] + E_M1[my_iindx[j + 1] - seq_length][cnt3][cnt4 / 2]
                                             );
              }
            }
        }
    }

    /* resize and move memory portions of energy matrix E_M2 */
    adjustArrayBoundaries(&matrices->E_M2[i],
                          &matrices->k_min_M2[i],
                          &matrices->k_max_M2[i],
                          &matrices->l_min_M2[i],
                          &matrices->l_max_M2[i],
                          min_k_real,
                          max_k_real,
                          min_l_real,
                          max_l_real
                          );
  } /* end for i */

  base_d1 = referenceBPs1[my_iindx[1] - seq_length];
  base_d2 = referenceBPs2[my_iindx[1] - seq_length];

  /* guess memory requirements for E_FcH, E_FcI and E_FcM */

  int min_k, max_k, max_l, min_l;
  int min_k_real, max_k_real, min_k_real_fcH, max_k_real_fcH, min_k_real_fcI, max_k_real_fcI, min_k_real_fcM, max_k_real_fcM;
  int *min_l_real, *max_l_real, *min_l_real_fcH, *max_l_real_fcH, *min_l_real_fcI, *max_l_real_fcI, *min_l_real_fcM, *max_l_real_fcM;

  min_k = min_l = 0;

  max_k = mm1[my_iindx[1] - seq_length] + referenceBPs1[my_iindx[1] - seq_length];
  max_l = mm2[my_iindx[1] - seq_length] + referenceBPs2[my_iindx[1] - seq_length];

#ifdef _OPENMP
#pragma omp sections
  {
#pragma omp section
    {
#endif
  prepareBoundaries(min_k,
                    max_k,
                    min_l,
                    max_l,
                    bpdist[my_iindx[1] - seq_length],
                    &matrices->k_min_Fc,
                    &matrices->k_max_Fc,
                    &matrices->l_min_Fc,
                    &matrices->l_max_Fc
                    );
  prepareArray(&matrices->E_Fc,
               matrices->k_min_Fc,
               matrices->k_max_Fc,
               matrices->l_min_Fc,
               matrices->l_max_Fc
               );
#ifdef _OPENMP
}


#pragma omp section
{
#endif
  prepareBoundaries(min_k,
                    max_k,
                    min_l,
                    max_l,
                    bpdist[my_iindx[1] - seq_length],
                    &matrices->k_min_FcH,
                    &matrices->k_max_FcH,
                    &matrices->l_min_FcH,
                    &matrices->l_max_FcH
                    );
  prepareArray(&matrices->E_FcH,
               matrices->k_min_FcH,
               matrices->k_max_FcH,
               matrices->l_min_FcH,
               matrices->l_max_FcH
               );
#ifdef _OPENMP
}
#pragma omp section
{
#endif
  prepareBoundaries(min_k,
                    max_k,
                    min_l,
                    max_l,
                    bpdist[my_iindx[1] - seq_length],
                    &matrices->k_min_FcI,
                    &matrices->k_max_FcI,
                    &matrices->l_min_FcI,
                    &matrices->l_max_FcI
                    );
  prepareArray(&matrices->E_FcI,
               matrices->k_min_FcI,
               matrices->k_max_FcI,
               matrices->l_min_FcI,
               matrices->l_max_FcI
               );
#ifdef _OPENMP
}
#pragma omp section
{
#endif
  prepareBoundaries(min_k,
                    max_k,
                    min_l,
                    max_l,
                    bpdist[my_iindx[1] - seq_length],
                    &matrices->k_min_FcM,
                    &matrices->k_max_FcM,
                    &matrices->l_min_FcM,
                    &matrices->l_max_FcM
                    );
  prepareArray(&matrices->E_FcM,
               matrices->k_min_FcM,
               matrices->k_max_FcM,
               matrices->l_min_FcM,
               matrices->l_max_FcM
               );
#ifdef _OPENMP
}
#pragma omp section
{
#endif
  preparePosteriorBoundaries(max_k - min_k + 1,
                             min_k,
                             &min_k_real,
                             &max_k_real,
                             &min_l_real,
                             &max_l_real
                             );
#ifdef _OPENMP
}
#pragma omp section
{
#endif
  preparePosteriorBoundaries(max_k - min_k + 1,
                             min_k,
                             &min_k_real_fcH,
                             &max_k_real_fcH,
                             &min_l_real_fcH,
                             &max_l_real_fcH
                             );
#ifdef _OPENMP
}
#pragma omp section
{
#endif
  preparePosteriorBoundaries(max_k - min_k + 1,
                             min_k,
                             &min_k_real_fcI,
                             &max_k_real_fcI,
                             &min_l_real_fcI,
                             &max_l_real_fcI
                             );

#ifdef _OPENMP
}
#pragma omp section
{
#endif
  preparePosteriorBoundaries(max_k - min_k + 1,
                             min_k,
                             &min_k_real_fcM,
                             &max_k_real_fcM,
                             &min_l_real_fcM,
                             &max_l_real_fcM
                             );
#ifdef _OPENMP
}
}
#endif

  /* begin actual energy calculations */
#ifdef _OPENMP
#pragma omp sections private(d, d1,d2,cnt1,cnt2,cnt3,cnt4,j, i, energy)
  {
#pragma omp section
    {
#endif
  for (d = TURN + 2; d <= seq_length; d++) /* i,j in [1..length] */
    for (j = d; j <= seq_length; j++) {
      unsigned int  u, ij;
      int           type, no_close;
      char          loopseq[10];
      i   = j - d + 1;
      ij  = my_iindx[i] - j;
      u   = seq_length - j + i - 1;
      if (u < TURN)
        continue;

      type = ptype[jindx[j] + i];

      no_close = (((type == 3) || (type == 4)) && no_closingGU);

      type = rtype[type];

      if (!type)
        continue;

      if (no_close)
        continue;

      d1  = base_d1 - referenceBPs1[ij];
      d2  = base_d2 - referenceBPs2[ij];
      if (u < 7) {
        strcpy(loopseq, sequence + j - 1);
        strncat(loopseq, sequence, i);
      }

      energy = E_Hairpin(u, type, S1[j + 1], S1[i - 1], loopseq, P);

      if (E_C_rem[ij] != INF)
        matrices->E_FcH_rem = MIN2(matrices->E_FcH_rem, E_C_rem[ij] + energy);

      if (!E_C[ij])
        continue;

      for (cnt1 = k_min_C[ij]; cnt1 <= k_max_C[ij]; cnt1++)
        for (cnt2 = l_min_C[ij][cnt1]; cnt2 <= l_max_C[ij][cnt1]; cnt2 += 2) {
          if (((cnt1 + d1) <= maxD1) && ((cnt2 + d2) <= maxD2)) {
            matrices->E_FcH[cnt1 + d1][(cnt2 + d2) / 2] = MIN2(matrices->E_FcH[cnt1 + d1][(cnt2 + d2) / 2],
                                                               energy + E_C[ij][cnt1][cnt2 / 2]
                                                               );
            updatePosteriorBoundaries(cnt1 + d1,
                                      cnt2 + d2,
                                      &min_k_real_fcH,
                                      &max_k_real_fcH,
                                      &min_l_real_fcH,
                                      &max_l_real_fcH
                                      );
          } else {
            matrices->E_FcH_rem = MIN2(matrices->E_FcH_rem, energy + E_C[ij][cnt1][cnt2 / 2]);
          }
        }
    }
  /* end of i-j loop */

  /* resize and move memory portions of energy matrix E_FcH */
  adjustArrayBoundaries(&matrices->E_FcH,
                        &matrices->k_min_FcH,
                        &matrices->k_max_FcH,
                        &matrices->l_min_FcH,
                        &matrices->l_max_FcH,
                        min_k_real_fcH,
                        max_k_real_fcH,
                        min_l_real_fcH,
                        max_l_real_fcH
                        );
#ifdef _OPENMP
}
#pragma omp section
{
#endif
  for (d = TURN + 2; d <= seq_length; d++) /* i,j in [1..length] */
    for (j = d; j <= seq_length; j++) {
      unsigned int  u, ij, p, q, pq;
      int           type, type_2, no_close;
      i   = j - d + 1;
      ij  = my_iindx[i] - j;
      u   = seq_length - j + i - 1;
      if (u < TURN)
        continue;

      type = ptype[jindx[j] + i];

      no_close = (((type == 3) || (type == 4)) && no_closingGU);

      type = rtype[type];

      if (!type)
        continue;

      if (no_close)
        continue;

      if (E_C_rem[ij] != INF) {
        for (p = j + 1; p < seq_length; p++) {
          unsigned int u1, qmin, ln_pre;
          u1 = p - j - 1;
          if (u1 + i - 1 > MAXLOOP)
            break;

          qmin    = p + TURN + 1;
          ln_pre  = u1 + i + seq_length;
          if (ln_pre > qmin + MAXLOOP)
            qmin = ln_pre - MAXLOOP - 1;

          for (q = qmin; q <= seq_length; q++) {
            unsigned int u2;
            pq      = my_iindx[p] - q;
            type_2  = rtype[(unsigned int)ptype[jindx[q] + p]];
            if (type_2 == 0)
              continue;

            u2 = i - 1 + seq_length - q;
            if (u1 + u2 > MAXLOOP)
              continue;

            /* get distance to reference if closing the interior loop
             *  d2a = dbp(T1_[1,n}, T1_{p,q} + T1_{i,j})
             *  d2b = dbp(T2_[1,n}, T2_{p,q} + T2_{i,j})
             */
            d1      = base_d1 - referenceBPs1[ij] - referenceBPs1[pq];
            d2      = base_d2 - referenceBPs2[ij] - referenceBPs2[pq];
            energy  = E_IntLoop(u1, u2, type, type_2, S1[j + 1], S1[i - 1], S1[p - 1], S1[q + 1], P);

            if (E_C_rem[pq] != INF)
              matrices->E_FcI_rem = MIN2(matrices->E_FcI_rem, E_C_rem[ij] + E_C_rem[pq] + energy);

            if (E_C[pq]) {
              for (cnt1 = k_min_C[pq];
                   cnt1 <= k_max_C[pq];
                   cnt1++)
                for (cnt2 = l_min_C[pq][cnt1];
                     cnt2 <= l_max_C[pq][cnt1];
                     cnt2 += 2)
                  matrices->E_FcI_rem = MIN2(matrices->E_FcI_rem, E_C_rem[ij] + E_C[pq][cnt1][cnt2 / 2] + energy);
            }
          }
        }
      }

      if (E_C[ij]) {
        for (p = j + 1; p < seq_length; p++) {
          unsigned int u1, qmin, ln_pre;
          u1 = p - j - 1;
          if (u1 + i - 1 > MAXLOOP)
            break;

          qmin    = p + TURN + 1;
          ln_pre  = u1 + i + seq_length;
          if (ln_pre > qmin + MAXLOOP)
            qmin = ln_pre - MAXLOOP - 1;

          for (q = qmin; q <= seq_length; q++) {
            unsigned int u2;
            pq      = my_iindx[p] - q;
            type_2  = rtype[(unsigned int)ptype[jindx[q] + p]];
            if (type_2 == 0)
              continue;

            u2 = i - 1 + seq_length - q;
            if (u1 + u2 > MAXLOOP)
              continue;

            /* get distance to reference if closing the interior loop
             *  d2a = dbp(T1_[1,n}, T1_{p,q} + T1_{i,j})
             *  d2b = dbp(T2_[1,n}, T2_{p,q} + T2_{i,j})
             */
            d1      = base_d1 - referenceBPs1[ij] - referenceBPs1[pq];
            d2      = base_d2 - referenceBPs2[ij] - referenceBPs2[pq];
            energy  = E_IntLoop(u1, u2, type, type_2, S1[j + 1], S1[i - 1], S1[p - 1], S1[q + 1], P);
            if (E_C_rem[pq] != INF) {
              for (cnt1 = k_min_C[ij];
                   cnt1 <= k_max_C[ij];
                   cnt1++)
                for (cnt2 = l_min_C[ij][cnt1];
                     cnt2 <= l_max_C[ij][cnt1];
                     cnt2 += 2)
                  matrices->E_FcI_rem = MIN2(matrices->E_FcI_rem, E_C[ij][cnt1][cnt2 / 2] + E_C_rem[pq] + energy);
            }

            if (E_C[pq]) {
              for (cnt1 = k_min_C[ij];
                   cnt1 <= k_max_C[ij];
                   cnt1++)
                for (cnt2 = l_min_C[ij][cnt1];
                     cnt2 <= l_max_C[ij][cnt1];
                     cnt2 += 2)
                  for (cnt3 = k_min_C[pq];
                       cnt3 <= k_max_C[pq];
                       cnt3++)
                    for (cnt4 = l_min_C[pq][cnt3];
                         cnt4 <= l_max_C[pq][cnt3];
                         cnt4 += 2) {
                      if (((cnt1 + cnt3 + d1) <= maxD1) && ((cnt2 + cnt4 + d2) <= maxD2)) {
                        matrices->E_FcI[cnt1 + cnt3 + d1][(cnt2 + cnt4 + d2) / 2] = MIN2(
                          matrices->E_FcI[cnt1 + cnt3 + d1][(cnt2 + cnt4 + d2) / 2],
                          E_C[ij][cnt1][cnt2 / 2]
                          + E_C[pq][cnt3][cnt4 / 2]
                          + energy
                          );
                        updatePosteriorBoundaries(cnt1 + cnt3 + d1,
                                                  cnt2 + cnt4 + d2,
                                                  &min_k_real_fcI,
                                                  &max_k_real_fcI,
                                                  &min_l_real_fcI,
                                                  &max_l_real_fcI
                                                  );
                      } else {
                        matrices->E_FcI_rem = MIN2(
                          matrices->E_FcI_rem,
                          E_C[ij][cnt1][cnt2 / 2]
                          + E_C[pq][cnt3][cnt4 / 2]
                          + energy
                          );
                      }
                    }
            }
          }
        }
      }
    }
  /* end of i-j loop */

  /* resize and move memory portions of energy matrix E_FcI */
  adjustArrayBoundaries(&matrices->E_FcI,
                        &matrices->k_min_FcI,
                        &matrices->k_max_FcI,
                        &matrices->l_min_FcI,
                        &matrices->l_max_FcI,
                        min_k_real_fcI,
                        max_k_real_fcI,
                        min_l_real_fcI,
                        max_l_real_fcI
                        );
#ifdef _OPENMP
}
#pragma omp section
{
#endif
  if (seq_length > 2 * TURN) {
    for (i = TURN + 1; i < seq_length - 2 * TURN; i++) {
      /* get distancies to references
       * d3a = dbp(T1_[1,n}, T1_{1,k} + T1_{k+1, n})
       * d3b = dbp(T2_[1,n}, T2_{1,k} + T2_{k+1, n})
       */
      d1  = base_d1 - referenceBPs1[my_iindx[1] - i] - referenceBPs1[my_iindx[i + 1] - seq_length];
      d2  = base_d2 - referenceBPs2[my_iindx[1] - i] - referenceBPs2[my_iindx[i + 1] - seq_length];

      if (E_M_rem[my_iindx[1] - i] != INF) {
        if (matrices->E_M2[i + 1]) {
          for (cnt1 = matrices->k_min_M2[i + 1];
               cnt1 <= matrices->k_max_M2[i + 1];
               cnt1++)
            for (cnt2 = matrices->l_min_M2[i + 1][cnt1];
                 cnt2 <= matrices->l_max_M2[i + 1][cnt1];
                 cnt2 += 2)
              matrices->E_FcM_rem = MIN2(matrices->E_FcM_rem, E_M_rem[my_iindx[1] - i] + matrices->E_M2[i + 1][cnt1][cnt2 / 2] + P->MLclosing);
        }

        if (matrices->E_M2_rem[i + 1] != INF)
          matrices->E_FcM_rem = MIN2(matrices->E_FcM_rem, E_M_rem[my_iindx[1] - i] + matrices->E_M2_rem[i + 1] + P->MLclosing);
      }

      if (matrices->E_M2_rem[i + 1] != INF) {
        if (E_M[my_iindx[1] - i]) {
          for (cnt1 = k_min_M[my_iindx[1] - i];
               cnt1 <= k_max_M[my_iindx[1] - i];
               cnt1++)
            for (cnt2 = l_min_M[my_iindx[1] - i][cnt1];
                 cnt2 <= l_max_M[my_iindx[1] - i][cnt1];
                 cnt2 += 2)
              matrices->E_FcM_rem = MIN2(matrices->E_FcM_rem, E_M[my_iindx[1] - i][cnt1][cnt2 / 2] + matrices->E_M2_rem[i + 1] + P->MLclosing);
        }
      }

      if (!E_M[my_iindx[1] - i])
        continue;

      if (!matrices->E_M2[i + 1])
        continue;

      for (cnt1 = k_min_M[my_iindx[1] - i]; cnt1 <= k_max_M[my_iindx[1] - i]; cnt1++)
        for (cnt2 = l_min_M[my_iindx[1] - i][cnt1]; cnt2 <= l_max_M[my_iindx[1] - i][cnt1]; cnt2 += 2)
          for (cnt3 = matrices->k_min_M2[i + 1]; cnt3 <= matrices->k_max_M2[i + 1]; cnt3++)
            for (cnt4 = matrices->l_min_M2[i + 1][cnt3]; cnt4 <= matrices->l_max_M2[i + 1][cnt3]; cnt4 += 2) {
              if (((cnt1 + cnt3 + d1) <= maxD1) && ((cnt2 + cnt4 + d2) <= maxD2)) {
                matrices->E_FcM[cnt1 + cnt3 + d1][(cnt2 + cnt4 + d2) / 2] = MIN2(
                  matrices->E_FcM[cnt1 + cnt3 + d1][(cnt2 + cnt4 + d2) / 2],
                  E_M[my_iindx[1] - i][cnt1][cnt2 / 2]
                  + matrices->E_M2[i + 1][cnt3][cnt4 / 2]
                  + P->MLclosing
                  );
                updatePosteriorBoundaries(cnt1 + cnt3 + d1,
                                          cnt2 + cnt4 + d2,
                                          &min_k_real_fcM,
                                          &max_k_real_fcM,
                                          &min_l_real_fcM,
                                          &max_l_real_fcM
                                          );
              } else {
                matrices->E_FcM_rem = MIN2(
                  matrices->E_FcM_rem,
                  E_M[my_iindx[1] - i][cnt1][cnt2 / 2]
                  + matrices->E_M2[i + 1][cnt3][cnt4 / 2]
                  + P->MLclosing
                  );
              }
            }
    }
  }

  /* resize and move memory portions of energy matrix E_FcM */
  adjustArrayBoundaries(&matrices->E_FcM,
                        &matrices->k_min_FcM,
                        &matrices->k_max_FcM,
                        &matrices->l_min_FcM,
                        &matrices->l_max_FcM,
                        min_k_real_fcM,
                        max_k_real_fcM,
                        min_l_real_fcM,
                        max_l_real_fcM
                        );
#ifdef _OPENMP
}
}
#endif


  /* compute E_Fc_rem */
  matrices->E_Fc_rem  = MIN2(matrices->E_FcH_rem, matrices->E_FcI_rem);
  matrices->E_Fc_rem  = MIN2(matrices->E_Fc_rem, matrices->E_FcM_rem);
  /* add the case were structure is unfolded chain */
  if ((referenceBPs1[my_iindx[1] - seq_length] > maxD1) || (referenceBPs2[my_iindx[1] - seq_length] > maxD2))
    matrices->E_Fc_rem = MIN2(matrices->E_Fc_rem, 0);

  /* compute all E_Fc */
  for (cnt1 = matrices->k_min_FcH; cnt1 <= matrices->k_max_FcH; cnt1++)
    for (cnt2 = matrices->l_min_FcH[cnt1]; cnt2 <= matrices->l_max_FcH[cnt1]; cnt2 += 2) {
      matrices->E_Fc[cnt1][cnt2 / 2] = MIN2(matrices->E_Fc[cnt1][cnt2 / 2],
                                            matrices->E_FcH[cnt1][cnt2 / 2]
                                            );
      updatePosteriorBoundaries(cnt1,
                                cnt2,
                                &min_k_real,
                                &max_k_real,
                                &min_l_real,
                                &max_l_real
                                );
    }
  for (cnt1 = matrices->k_min_FcI; cnt1 <= matrices->k_max_FcI; cnt1++)
    for (cnt2 = matrices->l_min_FcI[cnt1]; cnt2 <= matrices->l_max_FcI[cnt1]; cnt2 += 2) {
      matrices->E_Fc[cnt1][cnt2 / 2] = MIN2(matrices->E_Fc[cnt1][cnt2 / 2],
                                            matrices->E_FcI[cnt1][cnt2 / 2]
                                            );
      updatePosteriorBoundaries(cnt1,
                                cnt2,
                                &min_k_real,
                                &max_k_real,
                                &min_l_real,
                                &max_l_real
                                );
    }
  for (cnt1 = matrices->k_min_FcM; cnt1 <= matrices->k_max_FcM; cnt1++)
    for (cnt2 = matrices->l_min_FcM[cnt1]; cnt2 <= matrices->l_max_FcM[cnt1]; cnt2 += 2) {
      matrices->E_Fc[cnt1][cnt2 / 2] = MIN2(matrices->E_Fc[cnt1][cnt2 / 2],
                                            matrices->E_FcM[cnt1][cnt2 / 2]
                                            );
      updatePosteriorBoundaries(cnt1,
                                cnt2,
                                &min_k_real,
                                &max_k_real,
                                &min_l_real,
                                &max_l_real
                                );
    }
  /* add the case were structure is unfolded chain */
  matrices->E_Fc[referenceBPs1[my_iindx[1] - seq_length]][referenceBPs2[my_iindx[1] - seq_length] / 2] = MIN2(matrices->E_Fc[referenceBPs1[my_iindx[1] - seq_length]][referenceBPs2[my_iindx[1] - seq_length] / 2],
                                                                                                              0);
  updatePosteriorBoundaries(referenceBPs1[my_iindx[1] - seq_length],
                            referenceBPs2[my_iindx[1] - seq_length],
                            &min_k_real,
                            &max_k_real,
                            &min_l_real,
                            &max_l_real
                            );


  adjustArrayBoundaries(&matrices->E_Fc,
                        &matrices->k_min_Fc,
                        &matrices->k_max_Fc,
                        &matrices->l_min_Fc,
                        &matrices->l_max_Fc,
                        min_k_real,
                        max_k_real,
                        min_l_real,
                        max_l_real
                        );
}


PRIVATE void
adjustArrayBoundaries(int ***array,
                      int *k_min,
                      int *k_max,
                      int **l_min,
                      int **l_max,
                      int k_min_post,
                      int k_max_post,
                      int *l_min_post,
                      int *l_max_post)
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
      memmove((int **)(*array), ((int **)(*array)) + k_diff_pre, sizeof(int *) * mem_size);
      memmove((int *)(*l_min), ((int *)(*l_min)) + k_diff_pre, sizeof(int) * mem_size);
      memmove((int *)(*l_max), ((int *)(*l_max)) + k_diff_pre, sizeof(int) * mem_size);
    }

    /* reallocating memory to actual size used */
    *array  += *k_min;
    *array  = (int **)realloc(*array, sizeof(int *) * mem_size);
    *array  -= k_min_post;

    *l_min  += *k_min;
    *l_min  = (int *)realloc(*l_min, sizeof(int) * mem_size);
    *l_min  -= k_min_post;

    *l_max  += *k_min;
    *l_max  = (int *)realloc(*l_max, sizeof(int) * mem_size);
    *l_max  -= k_min_post;

    /* adjust l dimension of array */
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
          memmove((int *)((*array)[cnt1]), (int *)((*array)[cnt1]) + start, sizeof(int) * mem_size);

        (*array)[cnt1] = (int *)realloc((*array)[cnt1], sizeof(int) * mem_size);

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
  free(l_min_post);
  free(l_max_post);
  *k_min  = k_min_post;
  *k_max  = k_max_post;
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


INLINE PRIVATE void
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


INLINE PRIVATE void
prepareArray(int  ***array,
             int  min_k,
             int  max_k,
             int  *min_l,
             int  *max_l)
{
  int i, j, mem;

  *array  = (int **)vrna_alloc(sizeof(int *) * (max_k - min_k + 1));
  *array  -= min_k;

  for (i = min_k; i <= max_k; i++) {
    mem         = (max_l[i] - min_l[i] + 1) / 2 + 1;
    (*array)[i] = (int *)vrna_alloc(sizeof(int) * mem);
    for (j = 0; j < mem; j++)
      (*array)[i][j] = INF;
    (*array)[i] -= min_l[i] / 2;
  }
}


INLINE PRIVATE void
prepareArray2(unsigned long ***array,
              int           min_k,
              int           max_k,
              int           *min_l,
              int           *max_l)
{
  int i, mem;

  *array  = (unsigned long **)vrna_alloc(sizeof(unsigned long *) * (max_k - min_k + 1));
  *array  -= min_k;

  for (i = min_k; i <= max_k; i++) {
    mem         = (max_l[i] - min_l[i] + 1) / 2 + 1;
    (*array)[i] = (unsigned long *)vrna_alloc(sizeof(unsigned long) * mem);
    (*array)[i] -= min_l[i] / 2;
  }
}


/*
 #################################
 # OLD API support               #
 #################################
 */

/* crosslink data from vars->compatibility to TwoDfold_vars structure */
PRIVATE INLINE void
crosslink(TwoDfold_vars *vars)
{
  vrna_fold_compound_t *c;
  vrna_mx_mfe_t *m;

  c                   = vars->compatibility;
  m                   = c->matrices;
  vars->sequence      = c->sequence;
  vars->seq_length    = c->length;
  vars->reference_pt1 = c->reference_pt1;
  vars->reference_pt2 = c->reference_pt2;
  vars->referenceBPs1 = c->referenceBPs1;
  vars->referenceBPs2 = c->referenceBPs2;
  vars->bpdist        = c->bpdist;
  vars->do_backtrack  = 1;
  vars->dangles       = c->params->model_details.dangles;
  vars->circ          = c->params->model_details.circ;
  vars->temperature   = c->params->model_details.temperature;
  vars->ptype         = c->ptype_pf_compat;
  vars->P             = c->params;
  vars->S             = c->sequence_encoding2;
  vars->S1            = c->sequence_encoding;
  vars->my_iindx      = c->iindx;
  vars->mm1           = c->mm1;
  vars->mm2           = c->mm2;
  vars->maxD1         = c->maxD1;
  vars->maxD2         = c->maxD2;

  vars->E_C           = m->E_C;
  vars->l_min_values  = m->l_min_C;
  vars->l_max_values  = m->l_max_C;
  vars->k_min_values  = m->k_min_C;
  vars->k_max_values  = m->k_max_C;

  vars->E_F5            = m->E_F5;
  vars->l_min_values_f  = m->l_min_F5;
  vars->l_max_values_f  = m->l_max_F5;
  vars->k_min_values_f  = m->k_min_F5;
  vars->k_max_values_f  = m->k_max_F5;

  vars->E_F3            = m->E_F3;
  vars->l_min_values_f3 = m->l_min_F3;
  vars->l_max_values_f3 = m->l_max_F3;
  vars->k_min_values_f3 = m->k_min_F3;
  vars->k_max_values_f3 = m->k_max_F3;

  vars->E_M             = m->E_M;
  vars->l_min_values_m  = m->l_min_M;
  vars->l_max_values_m  = m->l_max_M;
  vars->k_min_values_m  = m->k_min_M;
  vars->k_max_values_m  = m->k_max_M;

  vars->E_M1            = m->E_M1;
  vars->l_min_values_m1 = m->l_min_M1;
  vars->l_max_values_m1 = m->l_max_M1;
  vars->k_min_values_m1 = m->k_min_M1;
  vars->k_max_values_m1 = m->k_max_M1;

#ifdef COUNT_STATES
  vars->N_C   = m->N_C;
  vars->N_F5  = m->N_F5;
  vars->N_M   = m->N_M;
  vars->N_M1  = m->N_M1;
#endif

  vars->E_M2_rem        = m->E_M2_rem;
  vars->E_M2            = m->E_M2;
  vars->l_min_values_m2 = m->l_min_M2;
  vars->l_max_values_m2 = m->l_max_M2;
  vars->k_min_values_m2 = m->k_min_M2;
  vars->k_max_values_m2 = m->k_max_M2;

  vars->E_Fc  = m->E_Fc;
  vars->E_FcH = m->E_FcH;
  vars->E_FcI = m->E_FcI;
  vars->E_FcM = m->E_FcM;

  vars->E_Fc_rem  = m->E_Fc_rem;
  vars->E_FcH_rem = m->E_FcH_rem;
  vars->E_FcI_rem = m->E_FcI_rem;
  vars->E_FcM_rem = m->E_FcM_rem;

  vars->E_C_rem   = m->E_C_rem;
  vars->E_M_rem   = m->E_M_rem;
  vars->E_M1_rem  = m->E_M1_rem;
  vars->E_F5_rem  = m->E_F5_rem;
}


PUBLIC TwoDfold_vars *
get_TwoDfold_variables(const char *seq,
                       const char *structure1,
                       const char *structure2,
                       int        circ)
{
  vrna_md_t md;
  TwoDfold_vars *vars;
  vrna_fold_compound_t *c;
  vrna_mx_mfe_t *m;

  set_model_details(&md);
  md.circ = circ;

  vars                = (TwoDfold_vars *)vrna_alloc(sizeof(TwoDfold_vars));
  vars->compatibility = vrna_fold_compound_TwoD(seq, structure1, structure2, &md, VRNA_OPTION_MFE);

  crosslink(vars);

  return vars;
}


PUBLIC char *
TwoDfold_backtrack_f5(unsigned int  j,
                      int           k,
                      int           l,
                      TwoDfold_vars *vars)
{
  return vrna_backtrack5_TwoD(vars->compatibility, k, l, j);
}


PUBLIC void
destroy_TwoDfold_variables(TwoDfold_vars *vars)
{
  if (vars == NULL)
    return;

  vrna_fold_compound_free(vars->compatibility);

  free(vars);
}


PUBLIC vrna_sol_TwoD_t *
TwoDfoldList(TwoDfold_vars  *vars,
             int            distance1,
             int            distance2)
{
  vrna_sol_TwoD_t *sol;

  sol = vrna_mfe_TwoD(vars->compatibility, distance1, distance2);

  crosslink(vars);

  return sol;
}


PUBLIC void
update_TwoDfold_params(TwoDfold_vars *vars)
{
  vrna_md_t md;

  set_model_details(&md);

  free(vars->compatibility->params);
  vars->compatibility->params = vrna_params(&md);

  crosslink(vars);
}
