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

#include "ViennaRNA/utils/basic.h"
#include "ViennaRNA/params/default.h"
#include "ViennaRNA/fold_vars.h"
#include "ViennaRNA/loops/all.h"
#include "ViennaRNA/gquad.h"
#include "ViennaRNA/constraints/hard.h"
#include "ViennaRNA/constraints/soft.h"
#include "ViennaRNA/mfe.h"
#include "ViennaRNA/pf_multifold.h"
#include "ViennaRNA/combinatorics.h"
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
PRIVATE int
fill_arrays(vrna_fold_compound_t *fc);


PRIVATE void
postprocess_circular(vrna_fold_compound_t *fc);


PRIVATE FLT_OR_DBL
decompose_pair(vrna_fold_compound_t *fc,
               int                  i,
               int                  j,
               vrna_mx_pf_aux_ml_t  aux_mx_ml);


PRIVATE void
extract_dimer_props(vrna_fold_compound_t  *fc,
                    double                *F0AB,
                    double                *FAB,
                    double                *FcAB,
                    double                *FA,
                    double                *FB);


/*
 #################################
 # BEGIN OF FUNCTION DEFINITIONS #
 #################################
 */
PUBLIC FLT_OR_DBL
vrna_pf(vrna_fold_compound_t  *fc,
        char                  *structure)
{
  int               n;
  FLT_OR_DBL        Q, dG;
  vrna_md_t         *md;
  vrna_exp_param_t  *params;
  vrna_mx_pf_t      *matrices;

  dG = (FLT_OR_DBL)(INF / 100.);

  if (fc) {
    /* make sure, everything is set up properly to start partition function computations */
    if (!vrna_fold_compound_prepare(fc, VRNA_OPTION_PF)) {
      vrna_message_warning("vrna_pf@part_func.c: Failed to prepare vrna_fold_compound");
      return dG;
    }

    n         = fc->length;
    params    = fc->exp_params;
    matrices  = fc->exp_matrices;
    md        = &(params->model_details);

#ifdef _OPENMP
    /* Explicitly turn off dynamic threads */
    omp_set_dynamic(0);
#endif

#ifdef SUN4
    nonstandard_arithmetic();
#elif defined(HP9)
    fpsetfastmode(1);
#endif

    /* call user-defined recursion status callback function */
    if (fc->stat_cb)
      fc->stat_cb(VRNA_STATUS_PF_PRE, fc->auxdata);

    /* for now, multi-strand folding is implemented as additional grammar rule */
    if (fc->strands > 1)
      vrna_pf_multifold_prepare(fc);

    /* call user-defined grammar pre-condition callback function */
    if ((fc->aux_grammar) && (fc->aux_grammar->cb_proc))
      fc->aux_grammar->cb_proc(fc, VRNA_STATUS_PF_PRE, fc->aux_grammar->data);

    if (!fill_arrays(fc)) {
#ifdef SUN4
      standard_arithmetic();
#elif defined(HP9)
      fpsetfastmode(0);
#endif
      return dG;
    }

    if (md->circ)
      /* do post processing step for circular RNAs */
      postprocess_circular(fc);

    /* call user-defined grammar post-condition callback function */
    if ((fc->aux_grammar) && (fc->aux_grammar->cb_proc))
      fc->aux_grammar->cb_proc(fc, VRNA_STATUS_PF_POST, fc->aux_grammar->data);

    if (fc->strands > 1)
      vrna_gr_reset(fc);

    /* call user-defined recursion status callback function */
    if (fc->stat_cb)
      fc->stat_cb(VRNA_STATUS_PF_POST, fc->auxdata);

    switch (md->backtrack_type) {
      case 'C':
        Q = matrices->qb[fc->iindx[1] - n];
        break;

      case 'M':
        Q = matrices->qm[fc->iindx[1] - n];
        break;

      default:
        Q = (md->circ) ? matrices->qo : matrices->q[fc->iindx[1] - n];
        break;
    }

    /* ensemble free energy in Kcal/mol              */
    if (Q <= FLT_MIN)
      vrna_message_warning("pf_scale too large");

    if (fc->strands > 1) {
      /* check for rotational symmetry correction */
      unsigned int sym = vrna_rotational_symmetry(fc->sequence);
      Q /= (FLT_OR_DBL)sym;

      /* add interaction penalty */
      Q *= pow(params->expDuplexInit, (FLT_OR_DBL)(fc->strands - 1));
    }

    dG = (FLT_OR_DBL)((-log(Q) - n * log(params->pf_scale)) *
                      params->kT /
                      1000.0);

    if (fc->type == VRNA_FC_TYPE_COMPARATIVE)
      dG /= fc->n_seq;

    /* calculate base pairing probability matrix (bppm)  */
    if (md->compute_bpp) {
      vrna_pairing_probs(fc, structure);

#ifndef VRNA_DISABLE_BACKWARD_COMPATIBILITY

      /*
       *  Backward compatibility:
       *  This block may be removed if deprecated functions
       *  relying on the global variable "pr" vanish from within the package!
       */
      pr = matrices->probs;

#endif
    }

#ifdef SUN4
    standard_arithmetic();
#elif defined(HP9)
    fpsetfastmode(0);
#endif
  }

  return dG;
}


PUBLIC vrna_dimer_pf_t
vrna_pf_dimer(vrna_fold_compound_t  *fc,
              char                  *structure)
{
  vrna_dimer_pf_t X;

  X.F0AB = X.FAB = X.FcAB = X.FA = X.FB = 0.;

  if (fc) {
    (void)vrna_pf(fc, structure);

    /* backward compatibility partition function and ensemble energy computation */
    extract_dimer_props(fc,
                        &(X.F0AB),
                        &(X.FAB),
                        &(X.FcAB),
                        &(X.FA),
                        &(X.FB));
  }

  return X;
}


PUBLIC int
vrna_pf_float_precision(void)
{
  return sizeof(FLT_OR_DBL) == sizeof(float);
}


PUBLIC FLT_OR_DBL *
vrna_pf_substrands(vrna_fold_compound_t *fc,
                   size_t               complex_size)
{
  FLT_OR_DBL *Q_sub = NULL;

  if ((fc) &&
      (fc->strands >= complex_size) &&
      (fc->exp_matrices) &&
      (fc->exp_matrices->q)) {
    unsigned int      *ss, *se, *so;
    FLT_OR_DBL        Q;
    vrna_exp_param_t  *params;
    vrna_mx_pf_t      *matrices;

    ss        = fc->strand_start;
    se        = fc->strand_end;
    so        = fc->strand_order;
    params    = fc->exp_params;
    matrices  = fc->exp_matrices;

    Q_sub = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * (fc->strands - complex_size + 1));

    for (size_t i = 0; i < fc->strands - complex_size + 1; i++) {
      size_t start, end;
      start     = ss[so[i]];
      end       = se[so[i + complex_size - 1]];
      Q         = matrices->q[fc->iindx[start] - end];
      Q_sub[i]  = (-log(Q) - (end - start + 1) * log(params->pf_scale)) *
                  params->kT /
                  1000.0;
    }
  }

  return Q_sub;
}


PUBLIC FLT_OR_DBL
vrna_pf_add(FLT_OR_DBL  dG1,
            FLT_OR_DBL  dG2,
            double      kT)
{
  double  x1  = -(double)dG1 / kT;
  double  x2  = -(double)dG2 / kT;
  double  xs  = MAX2(x1, x2);

  return -kT * (xs + log(exp(x1 - xs) + exp(x2 - xs)));
}


/*
 #################################
 # STATIC helper functions below #
 #################################
 */
PRIVATE int
fill_arrays(vrna_fold_compound_t *fc)
{
  int                 n, i, j, k, ij, *my_iindx, *jindx, with_gquad, with_ud;
  FLT_OR_DBL          temp, Qmax, *q, *qb, *qm, *qm1, *q1k, *qln;
  double              max_real;
  vrna_ud_t           *domains_up;
  vrna_md_t           *md;
  vrna_mx_pf_t        *matrices;
  vrna_mx_pf_aux_el_t aux_mx_el;
  vrna_mx_pf_aux_ml_t aux_mx_ml;
  vrna_exp_param_t    *pf_params;

  n           = fc->length;
  my_iindx    = fc->iindx;
  jindx       = fc->jindx;
  matrices    = fc->exp_matrices;
  pf_params   = fc->exp_params;
  domains_up  = fc->domains_up;
  q           = matrices->q;
  qb          = matrices->qb;
  qm          = matrices->qm;
  qm1         = matrices->qm1;
  q1k         = matrices->q1k;
  qln         = matrices->qln;
  md          = &(pf_params->model_details);
  with_gquad  = md->gquad;

  with_ud = (domains_up && domains_up->exp_energy_cb && (!(fc->type == VRNA_FC_TYPE_COMPARATIVE)));
  Qmax    = 0;

  max_real = (sizeof(FLT_OR_DBL) == sizeof(float)) ? FLT_MAX : DBL_MAX;

  if (with_ud && domains_up->exp_prod_cb)
    domains_up->exp_prod_cb(fc, domains_up->data);

  /* no G-Quadruplexes for comparative partition function (yet) */
  if (with_gquad) {
    free(fc->exp_matrices->G);
    fc->exp_matrices->G = NULL;

    switch (fc->type) {
      case VRNA_FC_TYPE_SINGLE:
        fc->exp_matrices->G = get_gquad_pf_matrix(fc->sequence_encoding2,
                                                  fc->exp_matrices->scale,
                                                  fc->exp_params);
        break;

      case VRNA_FC_TYPE_COMPARATIVE:
        fc->exp_matrices->G = get_gquad_pf_matrix_comparative(fc->length,
                                                              fc->S_cons,
                                                              fc->S,
                                                              fc->a2s,
                                                              fc->exp_matrices->scale,
                                                              fc->n_seq,
                                                              fc->exp_params);
        break;
    }
  }

  /* init auxiliary arrays for fast exterior/multibranch loops */
  aux_mx_el = vrna_exp_E_ext_fast_init(fc);
  aux_mx_ml = vrna_exp_E_ml_fast_init(fc);

  /*array initialization ; qb,qm,q
   * qb,qm,q (i,j) are stored as ((n+1-i)*(n-i) div 2 + n+1-j */
  for (i = 1; i <= n; i++) {
    ij      = my_iindx[i] - i;
    qb[ij]  = 0.0;
  }

  for (j = 2; j <= n; j++) {
    for (i = j - 1; i >= 1; i--) {
      ij = my_iindx[i] - j;

      qb[ij] = decompose_pair(fc, i, j, aux_mx_ml);

      /* Multibranch loop */
      qm[ij] = vrna_exp_E_ml_fast(fc, i, j, aux_mx_ml);

      if (qm1) {
        temp = vrna_exp_E_ml_fast_qqm(aux_mx_ml)[i]; /* for stochastic backtracking and circfold */

        /* apply auxiliary grammar rule for multibranch loop (M1) case */
        if ((fc->aux_grammar) && (fc->aux_grammar->cb_aux_exp_m1))
          temp += fc->aux_grammar->cb_aux_exp_m1(fc, i, j, fc->aux_grammar->data);

        qm1[jindx[j] + i] = temp;
      }

      /* Exterior loop */
      q[ij] = vrna_exp_E_ext_fast(fc, i, j, aux_mx_el);

      /* apply auxiliary grammar rule (storage takes place in user-defined data structure */
      if ((fc->aux_grammar) && (fc->aux_grammar->cb_aux_exp))
        fc->aux_grammar->cb_aux_exp(fc, i, j, fc->aux_grammar->data);

      if (q[ij] > Qmax) {
        Qmax = q[ij];
        if (Qmax > max_real / 10.)
          vrna_message_warning("Q close to overflow: %d %d %g", i, j, q[ij]);
      }

      if (q[ij] >= max_real) {
        vrna_message_warning("overflow while computing partition function for segment q[%d,%d]\n"
                             "use larger pf_scale", i, j);

        vrna_exp_E_ml_fast_free(aux_mx_ml);
        vrna_exp_E_ext_fast_free(aux_mx_el);

        return 0; /* failure */
      }
    }

    /* rotate auxiliary arrays */
    vrna_exp_E_ext_fast_rotate(aux_mx_el);
    vrna_exp_E_ml_fast_rotate(aux_mx_ml);
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
  vrna_exp_E_ml_fast_free(aux_mx_ml);
  vrna_exp_E_ext_fast_free(aux_mx_el);

  return 1;
}


PRIVATE FLT_OR_DBL
decompose_pair(vrna_fold_compound_t *fc,
               int                  i,
               int                  j,
               vrna_mx_pf_aux_ml_t  aux_mx_ml)
{
  unsigned int  n;
  int           *jindx, *pscore;
  FLT_OR_DBL    contribution;
  double        kTn;
  vrna_hc_t     *hc;

  contribution  = 0.;
  n             = fc->length;
  hc            = fc->hc;

  if (hc->mx[j * n + i]) {
    /* process hairpin loop(s) */
    contribution += vrna_exp_E_hp_loop(fc, i, j);
    /* process interior loop(s) */
    contribution += vrna_exp_E_int_loop(fc, i, j);
    /* process multibranch loop(s) */
    contribution += vrna_exp_E_mb_loop_fast(fc, i, j, aux_mx_ml);

    if ((fc->aux_grammar) && (fc->aux_grammar->cb_aux_exp_c))
      contribution += fc->aux_grammar->cb_aux_exp_c(fc, i, j, fc->aux_grammar->data);

    if (fc->type == VRNA_FC_TYPE_COMPARATIVE) {
      jindx         = fc->jindx;
      pscore        = fc->pscore;
      kTn           = fc->exp_params->kT / 10.;  /* kT in cal/mol */
      contribution  *= exp(pscore[jindx[j] + i] / kTn);
    }
  }

  return contribution;
}


/*
 * calculate partition function for circular case
 * NOTE: this is the postprocessing step ONLY
 * You have to call fill_arrays first to calculate
 * complete circular case!!!
 */
PRIVATE void
postprocess_circular(vrna_fold_compound_t *fc)
{
  unsigned int      **a2s;
  int               u, p, q, k, turn, n, *my_iindx, *jindx, s;
  FLT_OR_DBL        *scale, *qb, *qm, *qm1, *qm2, qo, qho, qio, qmo,
                    qbt1, qot, expMLclosing, n_seq;
  unsigned char     eval;
  vrna_exp_param_t  *pf_params;
  vrna_mx_pf_t      *matrices;
  vrna_hc_t         *hc;
  vrna_sc_t         *sc, **scs;

  n             = fc->length;
  n_seq         = (fc->type == VRNA_FC_TYPE_SINGLE) ? 1 : fc->n_seq;
  matrices      = fc->exp_matrices;
  my_iindx      = fc->iindx;
  jindx         = fc->jindx;
  pf_params     = fc->exp_params;
  hc            = fc->hc;
  qb            = matrices->qb;
  qm            = matrices->qm;
  qm1           = matrices->qm1;
  qm2           = matrices->qm2;
  scale         = matrices->scale;
  expMLclosing  = pf_params->expMLclosing;
  turn          = pf_params->model_details.min_loop_size;
  hc            = fc->hc;
  sc            = (fc->type == VRNA_FC_TYPE_SINGLE) ? fc->sc : NULL;
  scs           = (fc->type == VRNA_FC_TYPE_COMPARATIVE) ? fc->scs : NULL;
  a2s           = (fc->type == VRNA_FC_TYPE_COMPARATIVE) ? fc->a2s : NULL;
  qo            = qho = qio = qmo = 0.;

  for (p = 1; p < n; p++) {
    for (q = p + turn + 1; q <= n; q++) {
      u = n - q + p - 1;
      if (u < turn)
        continue;

      /* 1. get exterior hairpin contribution  */
      qho += qb[my_iindx[p] - q] *
             vrna_exp_E_hp_loop(fc, q, p);

      /* 2. get exterior interior loop contribution */
      qio += qb[my_iindx[p] - q] *
             vrna_exp_E_int_loop(fc, q, p);
    }
  } /* end of pq double loop */

  /* 3. Multiloops  */

  /* construct qm2 matrix for exterior multibranch loop computation */
  if (hc->f) {
    switch (fc->type) {
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
    switch (fc->type) {
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
    switch (fc->type) {
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

        qbt1 *= pow(expMLclosing, fc->n_seq);
        break;
    }
  } else {
    switch (fc->type) {
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

        qbt1 *= pow(expMLclosing, fc->n_seq);
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

    switch (fc->type) {
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
          for (s = 0; s < fc->n_seq; s++)
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


PRIVATE void
extract_dimer_props(vrna_fold_compound_t  *fc,
                    double                *F0AB,  /**< @brief Null model without DuplexInit */
                    double                *FAB,   /**< @brief all states with DuplexInit correction */
                    double                *FcAB,  /**< @brief true hybrid states only */
                    double                *FA,    /**< @brief monomer A */
                    double                *FB     /**< @brief monomer B */
                    )
{
  unsigned int      n, sym, *ss, *so, *se;
  double            kT, QAB, QToT, Qzero;
  vrna_mx_pf_t      *matrices;
  vrna_exp_param_t  *params;

  n         = fc->length;
  ss        = fc->strand_start;
  se        = fc->strand_end;
  so        = fc->strand_order;
  params    = fc->exp_params;
  matrices  = fc->exp_matrices;

  if (fc->strands > 1) {
    kT  = params->kT / 1000.0;
    QAB = matrices->q[fc->iindx[1] - n];

    /* check for rotational symmetry correction */
    sym = vrna_rotational_symmetry(fc->sequence);
    QAB /= (FLT_OR_DBL)sym;

    /* add interaction penalty */
    QAB *= pow(params->expDuplexInit, (FLT_OR_DBL)(fc->strands - 1));

    Qzero = matrices->q[fc->iindx[1] - n] +
            matrices->q[fc->iindx[1] - se[so[0]]] *
            matrices->q[fc->iindx[ss[so[1]]] - n];

    QToT = matrices->q[fc->iindx[1] - se[so[0]]] *
           matrices->q[fc->iindx[ss[so[1]]] - n] +
           QAB;

    *FAB  = -kT * (log(QToT) + n * log(params->pf_scale));
    *F0AB = -kT * (log(Qzero) + n * log(params->pf_scale));
    *FcAB = (QAB > 1e-17) ? -kT * (log(QAB) + n * log(params->pf_scale)) : 999;
    *FA   = -kT *
            (log(matrices->q[fc->iindx[1] - se[so[0]]]) + (se[so[0]]) *
             log(params->pf_scale));
    *FB = -kT *
          (log(matrices->q[fc->iindx[ss[so[1]]] - n]) + (n - ss[so[1]] + 1) *
           log(params->pf_scale));
  } else {
    *FA = *FB = *FAB = *F0AB = (-log(matrices->q[fc->iindx[1] - n]) - n * log(params->pf_scale)) *
                               params->kT / 1000.0;
    *FcAB = 0;
  }
}
