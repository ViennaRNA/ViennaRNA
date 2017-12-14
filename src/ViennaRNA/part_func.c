/*
 *                partiton function for RNA secondary structures
 *
 *                Ivo L Hofacker + Ronny Lorenz
 *                Vienna RNA package
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>    /* #defines FLT_MAX ... */
#include <limits.h>

#include "ViennaRNA/utils.h"
#include "ViennaRNA/energy_par.h"
#include "ViennaRNA/fold_vars.h"
#include "ViennaRNA/loop_energies.h"
#include "ViennaRNA/gquad.h"
#include "ViennaRNA/constraints.h"
#include "ViennaRNA/mfe.h"
#include "ViennaRNA/part_func.h"

#ifdef _OPENMP
#include <omp.h>
#endif

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
PRIVATE void
fill_arrays(vrna_fold_compound_t *vc);


PRIVATE void
postprocess_circular(vrna_fold_compound_t *vc);


/*
 #################################
 # BEGIN OF FUNCTION DEFINITIONS #
 #################################
 */
PUBLIC float
vrna_pf(vrna_fold_compound_t  *vc,
        char                  *structure)
{
  int               n;
  FLT_OR_DBL        Q;
  double            free_energy;
  vrna_md_t         *md;
  vrna_exp_param_t  *params;
  vrna_mx_pf_t      *matrices;

  free_energy = (float)(INF / 100.);

  if (vc) {
    /* make sure, everything is set up properly to start partition function computations */
    if (!vrna_fold_compound_prepare(vc, VRNA_OPTION_PF)) {
      vrna_message_warning("vrna_pf@part_func.c: Failed to prepare vrna_fold_compound");
      return free_energy;
    }

    n         = vc->length;
    params    = vc->exp_params;
    matrices  = vc->exp_matrices;
    md        = &(params->model_details);

#ifdef _OPENMP
    /* Explicitly turn off dynamic threads */
    omp_set_dynamic(0);
#endif

#ifdef SUN4
    nonstandard_arithmetic();
#else
#ifdef HP9
    fpsetfastmode(1);
#endif
#endif

    /* call user-defined recursion status callback function */
    if (vc->stat_cb)
      vc->stat_cb(VRNA_STATUS_PF_PRE, vc->auxdata);

    fill_arrays(vc);

    if (md->circ)
      /* do post processing step for circular RNAs */
      postprocess_circular(vc);

    /* call user-defined recursion status callback function */
    if (vc->stat_cb)
      vc->stat_cb(VRNA_STATUS_PF_POST, vc->auxdata);

    /* calculate base pairing probability matrix (bppm)  */
    if (md->compute_bpp) {
      vrna_pairing_probs(vc, structure);

#ifndef VRNA_DISABLE_BACKWARD_COMPATIBILITY

      /*
       *  Backward compatibility:
       *  This block may be removed if deprecated functions
       *  relying on the global variable "pr" vanish from within the package!
       */
      pr = matrices->probs;

#endif
    }

    switch (md->backtrack_type) {
      case 'C':
        Q = matrices->qb[vc->iindx[1] - n];
        break;

      case 'M':
        Q = matrices->qm[vc->iindx[1] - n];
        break;

      default:
        Q = (md->circ) ? matrices->qo : matrices->q[vc->iindx[1] - n];
        break;
    }

    /* ensemble free energy in Kcal/mol              */
    if (Q <= FLT_MIN)
      vrna_message_warning("pf_scale too large");

    free_energy = (-log(Q) - n * log(params->pf_scale)) *
                  params->kT /
                  1000.0;

    if (vc->type == VRNA_FC_TYPE_COMPARATIVE)
      free_energy /= vc->n_seq;

#ifdef SUN4
    standard_arithmetic();
#else
#ifdef HP9
    fpsetfastmode(0);
#endif
#endif
  }

  return free_energy;
}


PUBLIC vrna_dimer_pf_t
vrna_pf_dimer(vrna_fold_compound_t  *vc,
              char                  *structure)
{
  int               n;
  FLT_OR_DBL        Q;
  vrna_dimer_pf_t   X;
  double            free_energy;
  char              *sequence;
  vrna_md_t         *md;
  vrna_exp_param_t  *params;
  vrna_mx_pf_t      *matrices;

  if (!vrna_fold_compound_prepare(vc, VRNA_OPTION_PF | VRNA_OPTION_HYBRID)) {
    vrna_message_warning("vrna_pf_dimer@part_func_co.c: Failed to prepare vrna_fold_compound");
    X.FA = X.FB = X.FAB = X.F0AB = X.FcAB = 0;
    return X;
  }

  params    = vc->exp_params;
  n         = vc->length;
  md        = &(params->model_details);
  matrices  = vc->exp_matrices;
  sequence  = vc->sequence;

#ifdef _OPENMP
  /* Explicitly turn off dynamic threads */
  omp_set_dynamic(0);
#endif

#ifdef SUN4
  nonstandard_arithmetic();
#else
#ifdef HP9
  fpsetfastmode(1);
#endif
#endif

  /* hard code min_loop_size to 0, since we can not be sure yet that this is already the case */
  md->min_loop_size = 0;

  /* call user-defined recursion status callback function */
  if (vc->stat_cb)
    vc->stat_cb(VRNA_STATUS_PF_PRE, vc->auxdata);

  fill_arrays(vc);

  /* call user-defined recursion status callback function */
  if (vc->stat_cb)
    vc->stat_cb(VRNA_STATUS_PF_POST, vc->auxdata);

  if (md->backtrack_type == 'C')
    Q = matrices->qb[vc->iindx[1] - n];
  else if (md->backtrack_type == 'M')
    Q = matrices->qm[vc->iindx[1] - n];
  else
    Q = matrices->q[vc->iindx[1] - n];

  /* ensemble free energy in Kcal/mol */
  if (Q <= FLT_MIN)
    vrna_message_warning("pf_scale too large");

  free_energy = (-log(Q) - n * log(params->pf_scale)) * params->kT / 1000.0;
  /* in case we abort because of floating point errors */
  if (n > 1600)
    vrna_message_info(stderr, "free energy = %8.2f", free_energy);

  /* probability of molecules being bound together */

  /* Computation of "real" Partition function */
  /* Need that for concentrations */
  if (vc->cutpoint > 0) {
    double kT, QAB, QToT, Qzero;
    kT    = params->kT / 1000.0;
    Qzero = matrices->q[vc->iindx[1] - n];
    QAB   =
      (matrices->q[vc->iindx[1] - n] - matrices->q[vc->iindx[1] - (vc->cutpoint - 1)] *
       matrices->q[vc->iindx[vc->cutpoint] - n]) * params->expDuplexInit;
    /*correction for symmetry*/
    if ((n - (vc->cutpoint - 1) * 2) == 0)
      if ((strncmp(sequence, sequence + vc->cutpoint - 1, vc->cutpoint - 1)) == 0)
        QAB /= 2;

    QToT = matrices->q[vc->iindx[1] - (vc->cutpoint - 1)] *
           matrices->q[vc->iindx[vc->cutpoint] - n] + QAB;
    X.FAB   = -kT * (log(QToT) + n * log(params->pf_scale));
    X.F0AB  = -kT * (log(Qzero) + n * log(params->pf_scale));
    X.FcAB  = (QAB > 1e-17) ? -kT * (log(QAB) + n * log(params->pf_scale)) : 999;
    X.FA    = -kT *
              (log(matrices->q[vc->iindx[1] - (vc->cutpoint - 1)]) + (vc->cutpoint - 1) *
               log(params->pf_scale));
    X.FB = -kT *
           (log(matrices->q[vc->iindx[vc->cutpoint] - n]) + (n - vc->cutpoint + 1) *
            log(params->pf_scale));

    /* printf("QAB=%.9f\tQtot=%.9f\n",QAB/scale[n],QToT/scale[n]); */
  } else {
    X.FA    = X.FB = X.FAB = X.F0AB = free_energy;
    X.FcAB  = 0;
  }

  /* backtracking to construct binding probabilities of pairs */
  if (md->compute_bpp) {
    vrna_pairing_probs(vc, structure);

#ifndef VRNA_DISABLE_BACKWARD_COMPATIBILITY

    /*
     *  Backward compatibility:
     *  This block may be removed if deprecated functions
     *  relying on the global variable "pr" vanish from within the package!
     */
    pr = vc->exp_matrices->probs;

#endif
  }

#ifdef SUN4
  standard_arithmetic();
#else
#ifdef HP9
  fpsetfastmode(0);
#endif
#endif

  return X;
}


PUBLIC int
vrna_pf_float_precision(void)
{
  return sizeof(FLT_OR_DBL) == sizeof(float);
}


/*
 #################################
 # STATIC helper functions below #
 #################################
 */
PRIVATE void
fill_arrays(vrna_fold_compound_t *vc)
{
  unsigned char       *hard_constraints;
  int                 n, i, j, k, ij, d, *my_iindx, *jindx, with_gquad, turn,
                      with_ud, hc_decompose, *pscore;
  FLT_OR_DBL          temp, Qmax, qbt1, *q, *qb, *qm, *qm1, *q1k, *qln;
  double              kTn, max_real;
  vrna_ud_t           *domains_up;
  vrna_md_t           *md;
  vrna_hc_t           *hc;
  vrna_mx_pf_t        *matrices;
  vrna_mx_pf_aux_el_t *aux_mx_el;
  vrna_mx_pf_aux_ml_t *aux_mx_ml;
  vrna_exp_param_t    *pf_params;

  n                 = vc->length;
  my_iindx          = vc->iindx;
  jindx             = vc->jindx;
  pscore            = (vc->type == VRNA_FC_TYPE_COMPARATIVE) ? vc->pscore : NULL;
  matrices          = vc->exp_matrices;
  pf_params         = vc->exp_params;
  kTn               = pf_params->kT / 10.;  /* kT in cal/mol */
  hc                = vc->hc;
  domains_up        = vc->domains_up;
  q                 = matrices->q;
  qb                = matrices->qb;
  qm                = matrices->qm;
  qm1               = matrices->qm1;
  q1k               = matrices->q1k;
  qln               = matrices->qln;
  md                = &(pf_params->model_details);
  with_gquad        = md->gquad;
  turn              = md->min_loop_size;
  hard_constraints  = hc->matrix;

  with_ud = (domains_up && domains_up->exp_energy_cb && (!(vc->type == VRNA_FC_TYPE_COMPARATIVE)));
  Qmax    = 0;

  max_real = (sizeof(FLT_OR_DBL) == sizeof(float)) ? FLT_MAX : DBL_MAX;

  if (with_ud && domains_up->exp_prod_cb)
    domains_up->exp_prod_cb(vc, domains_up->data);

  /* no G-Quadruplexes for comparative partition function (yet) */
  if (with_gquad && (!(vc->type == VRNA_FC_TYPE_COMPARATIVE))) {
    free(vc->exp_matrices->G);
    vc->exp_matrices->G = get_gquad_pf_matrix(vc->sequence_encoding2,
                                              vc->exp_matrices->scale,
                                              vc->exp_params);
  }

  /* init auxiliary arrays for fast exterior/multibranch loops */
  aux_mx_el = vrna_exp_E_ext_fast_init(vc);
  aux_mx_ml = vrna_exp_E_ml_fast_init(vc);

  /*array initialization ; qb,qm,q
   * qb,qm,q (i,j) are stored as ((n+1-i)*(n-i) div 2 + n+1-j */
  for (d = 0; d <= turn; d++)
    for (i = 1; i <= n - d; i++) {
      j       = i + d;
      ij      = my_iindx[i] - j;
      qb[ij]  = 0.0;
    }

  for (j = turn + 2; j <= n; j++) {
    for (i = j - turn - 1; i >= 1; i--) {
      ij            = my_iindx[i] - j;
      hc_decompose  = hard_constraints[jindx[j] + i];
      qbt1          = 0;

      if (hc_decompose) {
        /* process hairpin loop(s) */
        qbt1 += vrna_exp_E_hp_loop(vc, i, j);
        /* process interior loop(s) */
        qbt1 += vrna_exp_E_int_loop(vc, i, j);
        /* process multibranch loop(s) */
        qbt1 += vrna_exp_E_mb_loop_fast(vc, i, j, aux_mx_ml->qqm1);

        if (vc->type == VRNA_FC_TYPE_COMPARATIVE)
          qbt1 *= exp(pscore[jindx[j] + i] / kTn);
      }

      qb[ij] = qbt1;

      /* Multibranch loop */
      qm[ij] = vrna_exp_E_ml_fast(vc, i, j, aux_mx_ml);

      if (qm1)
        qm1[jindx[j] + i] = aux_mx_ml->qqm[i]; /* for stochastic backtracking and circfold */

      /* Exterior loop */
      q[ij] = temp = vrna_exp_E_ext_fast(vc, i, j, aux_mx_el);

      if (temp > Qmax) {
        Qmax = temp;
        if (Qmax > max_real / 10.)
          vrna_message_warning("Q close to overflow: %d %d %g", i, j, temp);
      }

      if (temp >= max_real)
        vrna_message_error("overflow in pf_fold while calculating q[%d,%d]\n"
                           "use larger pf_scale", i, j);
    }

    /* rotate auxiliary arrays */
    vrna_exp_E_ext_fast_rotate(vc, aux_mx_el);
    vrna_exp_E_ml_fast_rotate(vc, aux_mx_ml);
  }

  /* prefill linear qln, q1k arrays */
  if (q1k && qln) {
    for (k = 1; k <= n; k++) {
      q1k[k]  = q[my_iindx[1] - k];
      qln[k]  = q[my_iindx[k] - n];
    }
    q1k[0]      = 1.0;
    qln[n + 1]  = 1.0;
  }

  /* free memory occupied by auxiliary arrays for fast exterior/multibranch loops */
  vrna_exp_E_ml_fast_free(vc, aux_mx_ml);
  vrna_exp_E_ext_fast_free(vc, aux_mx_el);
}


/* calculate partition function for circular case */
/* NOTE: this is the postprocessing step ONLY     */
/* You have to call fill_arrays first to calculate  */
/* complete circular case!!!                      */
PRIVATE void
postprocess_circular(vrna_fold_compound_t *vc)
{
  short             *S1;
  unsigned int      **a2s;
  int               u, p, q, k, l, turn, n, *my_iindx, *jindx, s;
  FLT_OR_DBL        *scale, *qb, *qm, *qm1, *qm2, qo, qho, qio, qmo,
                    qbt1, qot, expMLclosing, n_seq;
  unsigned char     *hard_constraints, eval;
  vrna_exp_param_t  *pf_params;
  vrna_mx_pf_t      *matrices;
  vrna_hc_t         *hc;
  vrna_sc_t         *sc, **scs;

  n                 = vc->length;
  n_seq             = (vc->type == VRNA_FC_TYPE_SINGLE) ? 1 : vc->n_seq;
  matrices          = vc->exp_matrices;
  my_iindx          = vc->iindx;
  jindx             = vc->jindx;
  S1                = vc->sequence_encoding;
  pf_params         = vc->exp_params;
  hc                = vc->hc;
  qb                = matrices->qb;
  qm                = matrices->qm;
  qm1               = matrices->qm1;
  qm2               = matrices->qm2;
  scale             = matrices->scale;
  expMLclosing      = pf_params->expMLclosing;
  turn              = pf_params->model_details.min_loop_size;
  hc                = vc->hc;
  sc                = (vc->type == VRNA_FC_TYPE_SINGLE) ? vc->sc : NULL;
  scs               = (vc->type == VRNA_FC_TYPE_COMPARATIVE) ? vc->scs : NULL;
  hard_constraints  = hc->matrix;
  qo                = qho = qio = qmo = 0.;

  for (p = 1; p < n; p++) {
    for (q = p + turn + 1; q <= n; q++) {
      u = n - q + p - 1;
      if (u < turn)
        continue;

      /* 1. get exterior hairpin contribution  */
      qho += qb[my_iindx[p] - q] *
             vrna_exp_E_hp_loop(vc, q, p);

      /* 2. get exterior interior loop contribution */
      qio += qb[my_iindx[p] - q] *
             vrna_exp_E_int_loop(vc, q, p);
    }
  } /* end of pq double loop */

  /* 3. Multiloops  */

  /* construct qm2 matrix for exterior multibranch loop computation */
  if (hc->f) {
    switch (vc->type) {
      case VRNA_FC_TYPE_SINGLE:
        if ((sc) && (sc->exp_f)) {
          for (k = 1; k < n - turn - 1; k++) {
            qot = 0.;

            for (u = k + turn + 1; u < n - turn - 1; u++)
              if (hc->f(k, n, u, u + 1, VRNA_DECOMP_ML_ML_ML, hc->data)) {
                qot += qm1[jindx[u] + k] *
                       qm1[jindx[n] + (u + 1)] *
                       sc->exp_f(k, n, u, u + 1, VRNA_DECOMP_ML_ML_ML, sc->data);
              }

            qm2[k] = qot;
          }
        } else {
          for (k = 1; k < n - turn - 1; k++) {
            qot = 0.;

            for (u = k + turn + 1; u < n - turn - 1; u++)
              if (hc->f(k, n, u, u + 1, VRNA_DECOMP_ML_ML_ML, hc->data))
                qot += qm1[jindx[u] + k] *
                       qm1[jindx[n] + (u + 1)];

            qm2[k] = qot;
          }
        }

        break;

      case VRNA_FC_TYPE_COMPARATIVE:
        if (scs) {
          for (k = 1; k < n - turn - 1; k++) {
            qbt1 = 0.;

            for (u = k + turn + 1; u < n - turn - 1; u++) {
              if (hc->f(k, n, u, u + 1, VRNA_DECOMP_ML_ML_ML, hc->data)) {
                qot = qm1[jindx[u] + k] *
                      qm1[jindx[n] + (u + 1)];

                for (s = 0; s < n_seq; s++)
                  if ((scs[s]) && scs[s]->exp_f)
                    qot *= scs[s]->exp_f(k, n, u, u + 1, VRNA_DECOMP_ML_ML_ML, scs[s]->data);

                qbt1 += qot;
              }
            }

            qm2[k] = qbt1;
          }
        } else {
          for (k = 1; k < n - turn - 1; k++) {
            qot = 0.;

            for (u = k + turn + 1; u < n - turn - 1; u++)
              if (hc->f(k, n, u, u + 1, VRNA_DECOMP_ML_ML_ML, hc->data))
                qot += qm1[jindx[u] + k] *
                       qm1[jindx[n] + (u + 1)];

            qm2[k] = qot;
          }
        }

        break;
    }
  } else {
    switch (vc->type) {
      case VRNA_FC_TYPE_SINGLE:
        if ((sc) && (sc->exp_f)) {
          for (k = 1; k < n - turn - 1; k++) {
            qot = 0.;

            for (u = k + turn + 1; u < n - turn - 1; u++)
              qot += qm1[jindx[u] + k] *
                     qm1[jindx[n] + (u + 1)] *
                     sc->exp_f(k, n, u, u + 1, VRNA_DECOMP_ML_ML_ML, sc->data);

            qm2[k] = qot;
          }
        } else {
          for (k = 1; k < n - turn - 1; k++) {
            qot = 0.;

            for (u = k + turn + 1; u < n - turn - 1; u++)
              qot += qm1[jindx[u] + k] *
                     qm1[jindx[n] + (u + 1)];

            qm2[k] = qot;
          }
        }

        break;

      case VRNA_FC_TYPE_COMPARATIVE:
        if (scs) {
          for (k = 1; k < n - turn - 1; k++) {
            qbt1 = 0.;

            for (u = k + turn + 1; u < n - turn - 1; u++) {
              qot = qm1[jindx[u] + k] *
                    qm1[jindx[n] + (u + 1)];

              for (s = 0; s < n_seq; s++)
                if ((scs[s]) && (scs[s]->exp_f))
                  qot *= scs[s]->exp_f(k, n, u, u + 1, VRNA_DECOMP_ML_ML_ML, scs[s]->data);

              qbt1 += qot;
            }

            qm2[k] = qbt1;
          }
        } else {
          for (k = 1; k < n - turn - 1; k++) {
            qot = 0.;

            for (u = k + turn + 1; u < n - turn - 1; u++)
              qot += qm1[jindx[u] + k] *
                     qm1[jindx[n] + (u + 1)];

            qm2[k] = qot;
          }
        }

        break;
    }
  }

  qbt1 = 0.;
  /* go through exterior multibranch loop configurations */
  if (hc->f) {
    switch (vc->type) {
      case VRNA_FC_TYPE_SINGLE:
        if ((sc) && (sc->exp_f)) {
          for (k = turn + 2; k < n - 2 * turn - 3; k++)
            if (hc->f(1, n, k, k + 1, VRNA_DECOMP_ML_ML_ML, hc->data))
              qbt1 += qm[my_iindx[1] - k] *
                      qm2[k + 1] *
                      sc->exp_f(1, n, k, k + 1, VRNA_DECOMP_ML_ML_ML, sc->data);
        } else {
          for (k = turn + 2; k < n - 2 * turn - 3; k++)
            if (hc->f(1, n, k, k + 1, VRNA_DECOMP_ML_ML_ML, hc->data))
              qbt1 += qm[my_iindx[1] - k] *
                      qm2[k + 1];
        }

        qbt1 *= expMLclosing;
        break;

      case VRNA_FC_TYPE_COMPARATIVE:
        if (scs) {
          for (k = turn + 2; k < n - 2 * turn - 3; k++)
            if (hc->f(1, n, k, k + 1, VRNA_DECOMP_ML_ML_ML, hc->data)) {
              qot = qm[my_iindx[1] - k] *
                    qm2[k + 1];

              for (s = 0; s < n_seq; s++)
                if ((scs[s]) && (scs[s]->exp_f))
                  qot *= scs[s]->exp_f(1, n, k, k + 1, VRNA_DECOMP_ML_ML_ML, scs[s]->data);

              qbt1 += qot;
            }
        } else {
          for (k = turn + 2; k < n - 2 * turn - 3; k++)
            if (hc->f(1, n, k, k + 1, VRNA_DECOMP_ML_ML_ML, hc->data))
              qbt1 += qm[my_iindx[1] - k] *
                      qm2[k + 1];
        }

        qbt1 *= pow(expMLclosing, vc->n_seq);
        break;
    }
  } else {
    switch (vc->type) {
      case VRNA_FC_TYPE_SINGLE:
        if ((sc) && (sc->exp_f)) {
          for (k = turn + 2; k < n - 2 * turn - 3; k++)
            qbt1 += qm[my_iindx[1] - k] *
                    qm2[k + 1] *
                    sc->exp_f(1, n, k, k + 1, VRNA_DECOMP_ML_ML_ML, sc->data);
        } else {
          for (k = turn + 2; k < n - 2 * turn - 3; k++)
            qbt1 += qm[my_iindx[1] - k] *
                    qm2[k + 1];
        }

        qbt1 *= expMLclosing;
        break;

      case VRNA_FC_TYPE_COMPARATIVE:
        if (scs) {
          for (k = turn + 2; k < n - 2 * turn - 3; k++) {
            qot = qm[my_iindx[1] - k] *
                  qm2[k + 1];

            for (s = 0; s < n_seq; s++)
              if ((scs[s]) && (scs[s]->exp_f))
                qot *= scs[s]->exp_f(1, n, k, k + 1, VRNA_DECOMP_ML_ML_ML, scs[s]->data);

            qbt1 += qot;
          }
        } else {
          for (k = turn + 2; k < n - 2 * turn - 3; k++)
            qbt1 += qm[my_iindx[1] - k] *
                    qm2[k + 1];
        }

        qbt1 *= pow(expMLclosing, vc->n_seq);
        break;
    }
  }

  qmo += qbt1;

  /* add an additional pf of 1.0 to take the open chain into account too */
  eval = (hc->up_ext[1] >= n) ? 1 : 0;
  if (hc->f)
    eval = (hc->f(1, n, 1, n, VRNA_DECOMP_EXT_UP, hc->data)) ? eval : 0;

  if (eval) {
    qbt1 = scale[n];

    switch (vc->type) {
      case VRNA_FC_TYPE_SINGLE:
        if (sc) {
          if (sc->exp_energy_up)
            qbt1 *= sc->exp_energy_up[1][n];

          if (sc->exp_f)
            qbt1 *= sc->exp_f(1, n, 1, n, VRNA_DECOMP_EXT_UP, sc->data);
        }

        break;

      case VRNA_FC_TYPE_COMPARATIVE:
        if (scs) {
          for (s = 0; s < vc->n_seq; s++)
            if (scs[s])
              if (scs[s]->energy_up)
                qbt1 *= scs[s]->exp_energy_up[1][a2s[s][n]];
        }

        break;
    }
    qo += qbt1;
  }

  qo += qho + qio + qmo;

  matrices->qo  = qo;
  matrices->qho = qho;
  matrices->qio = qio;
  matrices->qmo = qmo;
}
