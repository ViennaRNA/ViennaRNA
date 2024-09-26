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
#include "ViennaRNA/utils/log.h"
#include "ViennaRNA/sequences/alphabet.h"
#include "ViennaRNA/params/default.h"
#include "ViennaRNA/fold_vars.h"
#include "ViennaRNA/constraints/hard.h"
#include "ViennaRNA/constraints/soft.h"
#include "ViennaRNA/mfe/global.h"
#include "ViennaRNA/combinatorics/basic.h"
#include "ViennaRNA/eval/exterior.h"
#include "ViennaRNA/eval/hairpin.h"
#include "ViennaRNA/eval/internal.h"
#include "ViennaRNA/eval/multibranch.h"
#include "ViennaRNA/partfunc/multifold.h"
#include "ViennaRNA/partfunc/exterior.h"
#include "ViennaRNA/partfunc/internal.h"
#include "ViennaRNA/partfunc/multibranch.h"
#include "ViennaRNA/partfunc/gquad.h"
#include "ViennaRNA/probabilities/basepairs.h"
#include "ViennaRNA/partfunc/global.h"

#include "ViennaRNA/intern/grammar_dat.h"

#include "ViennaRNA/constraints/exterior_sc_pf.inc"
#include "ViennaRNA/constraints/internal_sc_pf.inc"
#include "ViennaRNA/constraints/multibranch_sc_pf.inc"

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
      vrna_log_warning("vrna_pf@part_func.c: Failed to prepare vrna_fold_compound");
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
      fc->stat_cb(fc, VRNA_STATUS_PF_PRE, fc->auxdata);

    /* for now, multi-strand folding is implemented as additional grammar rule */
    if (fc->strands > 1)
      vrna_pf_multifold_prepare(fc);

    /* call user-defined grammar pre-condition callback function */
    if (fc->aux_grammar) {
      for (size_t i = 0; i < vrna_array_size(fc->aux_grammar->cbs_status); i++)
        if (fc->aux_grammar->cbs_status[i])
          fc->aux_grammar->cbs_status[i](fc, VRNA_STATUS_PF_PRE, fc->aux_grammar->datas[i]);
    }

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
    if (fc->aux_grammar) {
      for (size_t i = 0; i < vrna_array_size(fc->aux_grammar->cbs_status); i++)
        if (fc->aux_grammar->cbs_status[i])
          fc->aux_grammar->cbs_status[i](fc, VRNA_STATUS_PF_POST, fc->aux_grammar->datas[i]);
    }

    if (fc->strands > 1)
      vrna_gr_reset(fc);

    /* call user-defined recursion status callback function */
    if (fc->stat_cb)
      fc->stat_cb(fc, VRNA_STATUS_PF_POST, fc->auxdata);

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
      vrna_log_warning("pf_scale too large");

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
  FLT_OR_DBL          temp, Qmax, *q, *qb, *qm, *qm1, *qm2, *q1k, *qln;
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
  qm2         = matrices->qm2_real;
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
#ifndef VRNA_DISABLE_C11_FEATURES
    vrna_smx_csr_free(fc->exp_matrices->q_gq);
#else
    vrna_smx_csr_FLT_OR_DBL_free(fc->exp_matrices->q_gq);
#endif
    fc->exp_matrices->q_gq = vrna_gq_pos_pf(fc);
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

      if (qm2)
        qm2[ij] = vrna_exp_E_m2_fast(fc, i, j, aux_mx_ml);

      /* Multibranch loop */
      qm[ij] = vrna_exp_E_ml_fast(fc, i, j, aux_mx_ml);

      if (qm1) {
        temp = vrna_exp_E_ml_fast_qqm(aux_mx_ml)[i]; /* for stochastic backtracking and circfold */

        /* apply auxiliary grammar rule for multibranch loop (M1) case */
        if (fc->aux_grammar) {
          for (size_t c = 0; c < vrna_array_size(fc->aux_grammar->exp_m1); c++) {
            if (fc->aux_grammar->exp_m1[c].cb)
              temp += fc->aux_grammar->m1[c].cb(fc, i, j, fc->aux_grammar->m1[c].data);
          }
        }

        qm1[jindx[j] + i] = temp;
      }

      /* Exterior loop */
      q[ij] = vrna_exp_E_ext_fast(fc, i, j, aux_mx_el);

      /* apply auxiliary grammar rule (storage takes place in user-defined data structure */
      if (fc->aux_grammar) {
        for (size_t c = 0; c < vrna_array_size(fc->aux_grammar->exp_aux); c++) {
          if (fc->aux_grammar->exp_aux[c].cb)
            (void)fc->aux_grammar->exp_aux[c].cb(fc, i, j, fc->aux_grammar->exp_aux[c].data);
        }
      }

      if (q[ij] > Qmax) {
        Qmax = q[ij];
        if (Qmax > max_real / 10.)
          vrna_log_warning("Q close to overflow: %d %d %g", i, j, q[ij]);
      }

      if (q[ij] >= max_real) {
        vrna_log_warning("overflow while computing partition function for segment q[%d,%d]\n"
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
    contribution += vrna_exp_eval_hairpin(fc, i, j, VRNA_EVAL_LOOP_DEFAULT);
    /* process internal loop(s) */
    contribution += vrna_exp_E_int_loop(fc, i, j);
    /* process multibranch loop(s) */
    contribution += vrna_exp_E_mb_loop_fast(fc, i, j, aux_mx_ml);

    if (fc->aux_grammar) {
      for (size_t c = 0; c < vrna_array_size(fc->aux_grammar->exp_c); c++) {
        if (fc->aux_grammar->exp_c[c].cb)
          contribution += fc->aux_grammar->exp_c[c].cb(fc, i, j, fc->aux_grammar->exp_c[c].data);
      }
    }

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
  short             *S1, *S, **SS, **S5, **S3;
  unsigned int      **a2s, n_seq, s, u, p, q, i, j, k, turn, n, tt;
  int               *my_iindx;
  FLT_OR_DBL        *scale, *qb, *qm2, *qm1, q_g, qo, qho, qio,
                    qmo, qbt1, expMLclosing, *expMLbase;
  double            *expintern;

  vrna_smx_csr(FLT_OR_DBL)  *q_gq;
  unsigned char     eval, with_gquad;
  vrna_exp_param_t  *pf_params;
  vrna_mx_pf_t      *matrices;
  vrna_hc_t         *hc;
  vrna_sc_t         *sc, **scs;
  vrna_md_t         *md;

  n             = fc->length;
  n_seq         = (fc->type == VRNA_FC_TYPE_SINGLE) ? 1 : fc->n_seq;
  matrices      = fc->exp_matrices;
  my_iindx      = fc->iindx;
  pf_params     = fc->exp_params;
  md            = &(pf_params->model_details);
  hc            = fc->hc;
  S1            = (fc->type == VRNA_FC_TYPE_SINGLE) ? fc->sequence_encoding : NULL;
  S             = (fc->type == VRNA_FC_TYPE_SINGLE) ? fc->sequence_encoding2 : NULL;
  SS            = (fc->type == VRNA_FC_TYPE_SINGLE) ? NULL : fc->S;
  S5            = (fc->type == VRNA_FC_TYPE_SINGLE) ? NULL : fc->S5;
  S3            = (fc->type == VRNA_FC_TYPE_SINGLE) ? NULL : fc->S3;
  a2s           = (fc->type == VRNA_FC_TYPE_SINGLE) ? NULL : fc->a2s;
  qb            = matrices->qb;
  qm2           = matrices->qm2_real;
  qm1           = matrices->qm1_new;

  q_gq          = fc->exp_matrices->q_gq;
  scale         = matrices->scale;
  expMLclosing  = pf_params->expMLclosing;
  expMLbase     = fc->exp_matrices->expMLbase;
  expintern     = &(pf_params->expinternal[0]);
  turn          = md->min_loop_size;
  with_gquad    = (unsigned char)md->gquad;
  hc            = fc->hc;
  sc            = (fc->type == VRNA_FC_TYPE_SINGLE) ? fc->sc : NULL;
  scs           = (fc->type == VRNA_FC_TYPE_COMPARATIVE) ? fc->scs : NULL;
  qo            = qho = qio = qmo = 0.;

  struct sc_mb_exp_dat sc_mb_wrapper;
  struct sc_ext_exp_dat sc_ext_wrapper;
  struct sc_int_exp_dat sc_int_wrapper;
  init_sc_mb_exp(fc, &sc_mb_wrapper);
  init_sc_ext_exp(fc, &sc_ext_wrapper);
  init_sc_int_exp(fc, &sc_int_wrapper);

  for (p = 1; p < n; p++) {
    for (q = p + turn + 1; q <= n; q++) {
      u = n + p - q - 1;
      if (u < turn)
        continue;

      /* 1. get exterior hairpin contribution  */
      qho += qb[my_iindx[p] - q] *
             vrna_exp_eval_hairpin(fc, q, p, VRNA_EVAL_LOOP_DEFAULT);

      /* 2. get exterior internal loop contribution */
      qio += qb[my_iindx[p] - q] *
             vrna_exp_E_int_loop(fc, q, p);
    }
  } /* end of pq double loop */

  /* 3. Multiloops  */

  /* fill QM1 */
  for (j = 1; j < MIN2(turn + 2, VRNA_GQUAD_MIN_BOX_SIZE); j++)
    qm1[j] = 0.;

  for (j = MIN2(turn + 2, VRNA_GQUAD_MIN_BOX_SIZE); j <= n; j++) {
    /* regular base pairs */
    for (u = j - turn - 1; u >= 1; u--) {
      eval = (hc->up_ml[1] >= (u - 1)) ? (hc->mx[n * j + u] & VRNA_CONSTRAINT_CONTEXT_MB_LOOP) : 0;
      if ((hc->f) && (!hc->f(1, j, u, j, VRNA_DECOMP_ML_ML, hc->data)))
        eval = 0;

      if (eval) {
        qbt1    = qb[my_iindx[u] - j] *
                  expMLbase[u - 1];

        switch (fc->type) {
          case VRNA_FC_TYPE_COMPARATIVE:
            for (s = 0; s < n_seq; s++) {
              tt    = vrna_get_ptype_md(SS[s][u], SS[s][j], md);
              qbt1  *= vrna_exp_E_multibranch_stem(tt, S5[s][u], S3[s][j], pf_params);
            }

            break;

          default:
            tt    = vrna_get_ptype_md(S[u], S[j], md);
            qbt1  *= vrna_exp_E_multibranch_stem(tt, S1[u - 1], S1[j + 1], pf_params);
            break;
        }

        if (sc_mb_wrapper.red_ml)
          qbt1 *= sc_mb_wrapper.red_ml(1, j, u, j, &sc_mb_wrapper);

        qm1[j] += qbt1;
      }
    }

    /* g-quads */
    if (with_gquad) {
      if (j >= VRNA_GQUAD_MIN_BOX_SIZE) {
        for (u = j - VRNA_GQUAD_MIN_BOX_SIZE + 1; u >= 1; u--) {
          eval = (hc->up_ml[1] >= (u - 1)) ? 1 : 0;
          if ((hc->f) && (!hc->f(1, j, u, j, VRNA_DECOMP_ML_ML, hc->data)))
            eval = 0;

          if (eval) {
#ifndef VRNA_DISABLE_C11_FEATURES
            q_g = vrna_smx_csr_get(q_gq, u, j, 0.);
#else
            q_g = vrna_smx_csr_FLT_OR_DBL_get(q_gq, u, j, 0.);
#endif
            if (q_g != 0.) {
              qbt1    = q_g *
                        pow(vrna_exp_E_multibranch_stem(0, -1, -1, pf_params), (double)n_seq) *
                        expMLbase[u - 1];

              if (sc_mb_wrapper.red_ml)
                qbt1 *= sc_mb_wrapper.red_ml(1, j, u, j, &sc_mb_wrapper);

              qm1[j] += qbt1;
            }
          }
        }
      }
    }
  }

  /* use fM1_new and fM2 to construct segments with at least 3 branches */
  unsigned int space3 = 2 * MIN2(turn + 2, VRNA_GQUAD_MIN_BOX_SIZE);
  qbt1 = 0.;

  FLT_OR_DBL *qm1_tmp = qm1;

  /* apply hard constraints if necessary */
  if (hc->f) {
    qm1_tmp = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * (n + 2));
    qm1_tmp = memcpy(qm1_tmp, qm1, sizeof(FLT_OR_DBL) * (n + 2));

    for (k = turn + 2; k <= n; k++)
      if (!hc->f(1, n, k, k + 1, VRNA_DECOMP_ML_ML_ML, hc->data))
        qm1_tmp[k] = 0.;
  }

  /* apply soft constraints if necessary */
  if (sc_mb_wrapper.decomp_ml) {
    if (qm1_tmp == qm1) {
      qm1_tmp = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * (n + 2));
      qm1_tmp = memcpy(qm1_tmp, qm1, sizeof(FLT_OR_DBL) * (n + 2));
    }

    for (k = turn + 2; k <= n; k++)
      qm1_tmp[k] *= sc_mb_wrapper.decomp_ml(1, n, k, k + 1, &sc_mb_wrapper);
  }

  /* actual decomposition */
  for (k = MIN2(turn + 2, VRNA_GQUAD_MIN_BOX_SIZE);
       k + space3 <= n;
       k++) {
    qbt1 += qm1_tmp[k] *
            qm2[my_iindx[k + 1] - n];
  }

  qbt1 *= pow(expMLclosing, (double)n_seq);

  if (qm1_tmp != qm1)
    free(qm1_tmp);

  qmo += qbt1;

  /* add open chain, i.e. no base pairs or structure at all */
  eval = (hc->up_ext[1] >= n) ? 1 : 0;
  if (hc->f)
    eval = (hc->f(1, n, 1, n, VRNA_DECOMP_EXT_UP, hc->data)) ? eval : 0;

  if (eval) {
    qbt1 = scale[n];

    switch (fc->type) {
      case VRNA_FC_TYPE_COMPARATIVE:
        for (s = 0; s < n_seq; s++)
          qbt1 *= pow(vrna_exp_E_exterior_loop(a2s[s][n], md), 1. / (double)n_seq);
        break;

      default:
        qbt1 *= vrna_exp_E_exterior_loop(n, md);
        break;
    }

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

  if (with_gquad) {

    unsigned int              us1, us2, u1, u2, k, l;
    unsigned int              si, sj, sk, sl;

    /* consider all configurations where a G-quadruplex spans over the artificial cutpoint */

    unsigned int start = 1;

    if (n > VRNA_GQUAD_MAX_BOX_SIZE)
      start = n - VRNA_GQUAD_MAX_BOX_SIZE;

    /* loop over each possible start of a n,1 junction-spanning gquad */
    for (i = start; i <= n; i++) {
      unsigned int start_j = 1;
      unsigned int stop_j = n - 1;
      if ((n - i + 1) < VRNA_GQUAD_MIN_BOX_SIZE)
        start_j = VRNA_GQUAD_MIN_BOX_SIZE + i - n - 1;
      if (n > VRNA_GQUAD_MAX_BOX_SIZE)
        stop_j = VRNA_GQUAD_MAX_BOX_SIZE + i - n - 1;
      if (stop_j >= i)
        stop_j = i - 1;

      for (j = start_j; j <= stop_j; j++) {
#ifndef VRNA_DISABLE_C11_FEATURES
        q_g = vrna_smx_csr_get(q_gq, i, j, 0.);
#else
        q_g = vrna_smx_csr_FLT_OR_DBL_get(q_gq, i, j, 0.);
#endif

        if (q_g != 0.) {

          /* case 1: gquad is the only structure, rest is unpaired */
          if (i - j > 3) { /* keep at least 3 unpaired bases between start and end of gquad */
            u = i - j - 1;
            /* 1st, obey hard constraints */
            if (hc->up_ext[j + 1] >= u)
              eval = 1;
            if (hc->f)
              eval = (hc->f(j + 1, n, i, n, VRNA_DECOMP_EXT_EXT, hc->data)) ? eval : 0;
            if (eval) {
              qbt1 = q_g *
                     scale[u];

              if (md->circ_penalty)
                qbt1 *= pow(vrna_exp_E_hairpin(u, 0, -1, -1, NULL, pf_params), (double)n_seq);

              if (sc_ext_wrapper.red_ext)
                qbt1 *= sc_ext_wrapper.red_ext(j + 1, n, i, n, &sc_ext_wrapper);

              qho += qbt1;
            }
          } /* end case 1 */

          /* case 2.1: gquad forms an internal-loop like structure with another gquadruplex */
          for (u1 = 0, k = j + 1; k + VRNA_GQUAD_MIN_BOX_SIZE - 1 < i; k++, u1++) {
            /* obey hard constraints */
            if (hc->up_int[j + 1] < u1)
              break;

            if (u1 > MAXLOOP)
              break;

            unsigned int lmax = i - 1;
            u2 = i - lmax - 1;

            for (l = lmax; l >= k + VRNA_GQUAD_MIN_BOX_SIZE - 1; l--, u2++) {
              if (((u1 == 0) && (u2 < 3)) ||
                  ((u1 < 3) && (u2 == 0)))
                continue;

              /* obey hard constraints */
              if (hc->up_int[l + 1] < u2)
                break;

              if (u1 + u2 > MAXLOOP)
                break;

              /* obey user-defined hard constraints */
              if (hc->f) {
                vrna_log_debug("user-defined hard constraints not implemented for int-loop type gquads yet!");
              }

#ifndef VRNA_DISABLE_C11_FEATURES
              qbt1 = vrna_smx_csr_get(q_gq, k, l, 0.);
#else
              qbt1 = vrna_smx_csr_FLT_OR_DBL_get(q_gq, k, l, 0.);
#endif
              if (qbt1 != 0.) {
                qbt1 *= q_g * scale[u1 + u2];

                switch (fc->type) {
                  case VRNA_FC_TYPE_SINGLE:
                    qbt1 *= (FLT_OR_DBL)expintern[u1 + u2];
                    break;

                  case VRNA_FC_TYPE_COMPARATIVE:
                    for (s = 0; s < n_seq; s++) {
                      us1   = a2s[s][k - 1] - a2s[s][j];
                      us2   = a2s[s][i - 1] - a2s[s][l];
                      qbt1 *= (FLT_OR_DBL)expintern[us1 + us2];
                    }
                    break;
                }

                if (sc_int_wrapper.pair_ext)
                  qbt1 *= sc_int_wrapper.pair_ext(i, j, k, l, &sc_int_wrapper);

                qio += qbt1;
              }
            }
          } /* end case 2.1 */

          /* case 2.2: gquad forms an internal-loop like structure with another base pair */
          for (u1 = 0, k = j + 1; k + turn + 1 < i; k++, u1++) {
            /* obey hard constraints */
            if (hc->up_int[j + 1] < u1)
              break;

            if (u1 > MAXLOOP)
              break;

            unsigned int lmax = i - 1;
            u2 = i - lmax - 1;
            for (l = lmax; l > k + turn; l--, u2++) {
              if (((u1 == 0) && (u2 < 3)) ||
                  ((u1 < 3) && (u2 == 0)))
                continue;
              /* obey hard constraints */
              if (hc->up_int[l + 1] < u2)
                break;
              if (!(hc->mx[n * k + l] & (VRNA_CONSTRAINT_CONTEXT_INT_LOOP | VRNA_CONSTRAINT_CONTEXT_INT_LOOP_ENC)))
                continue;

              if (u1 + u2 > MAXLOOP)
                break;

              /* obey user-defined hard constraints */
              if (hc->f) {
                vrna_log_debug("user-defined hard constraints not implemented for int-loop type gquads yet!");
              }
          
              qbt1 = qb[my_iindx[k] - l];
              if (qbt1 != 0.) {

                qbt1 *= q_g *
                        scale[u1 + u2];

                switch (fc->type) {
                  case VRNA_FC_TYPE_SINGLE:
                    sl = S1[l + 1];
                    sk = S1[k - 1];

                    tt    = vrna_get_ptype_md(S[l], S[k], md);
                    if (md->dangles == 2)
                      qbt1 *= pf_params->expmismatchI[tt][sl][sk];

                    if (tt > 2)
                      qbt1 *= pf_params->expTermAU;
                    
                    qbt1 *= (FLT_OR_DBL)expintern[u1 + u2];
                    break;

                  case VRNA_FC_TYPE_COMPARATIVE:
                    for (s = 0; s < n_seq; s++) {
                      tt = vrna_get_ptype_md(SS[s][l], SS[s][k], md);
                      if (md->dangles == 2)
                        qbt1 *= pf_params->expmismatchI[tt][S3[s][l]][S5[s][k]];

                      if (tt > 2)
                        qbt1 *= pf_params->expTermAU;

                      us1   = a2s[s][k - 1] - a2s[s][j];
                      us2   = a2s[s][i - 1] - a2s[s][l];
                      qbt1 *= (FLT_OR_DBL)expintern[us1 + us2];
                    }
                    break;
                }

                if (sc_int_wrapper.pair_ext)
                  qbt1 *= sc_int_wrapper.pair_ext(i, j, k, l, &sc_int_wrapper);

                qio += qbt1;
              }
            }
          } /* end case 2.2 */

          /* case 3, gquad forms a multi-branch loop like structure with other base pairs or gquadruplexes */
          if (qm2[my_iindx[j + 1] - i + 1] != 0) {
            qbt1 = q_g *
                   qm2[my_iindx[j + 1] - i + 1] *
                   pow(vrna_exp_E_multibranch_stem(0, -1, -1, pf_params) * expMLclosing, (double)n_seq);
            qmo += qbt1;
          }


        }
      }
    }

    /* next case: everything unpaired except for one gquad somewhere not spanning artifical cutpoint */
    for (i = 1; i + VRNA_GQUAD_MIN_BOX_SIZE - 1 <= n; i++)
      for (j = i + VRNA_GQUAD_MIN_BOX_SIZE - 1; (j <= n) && j <= (i + VRNA_GQUAD_MAX_BOX_SIZE - 1); j++) {
        /* keep at least 3 unpaired bases in the loop around */
        if ((i - 1) + (n - j) < 3)
          break;
        
#ifndef VRNA_DISABLE_C11_FEATURES
        q_g = vrna_smx_csr_get(q_gq, i, j, 0.);
#else
        q_g = vrna_smx_csr_FLT_OR_DBL_get(q_gq, i, j, 0.);
#endif
        if (q_g != 0.) {
          qbt1 = 0.;
          
          /* obey constraints */
          u1 = i - 1;
          u2 = n - j;

          eval = (hc->up_ext[1] >= u1) ? 1 : 0;
          if (u2 > 0)
            eval = (hc->up_ext[j + 1] >= u2) ? eval : 0;
          if (hc->f) {
            if (u1 > 0)
              eval = (hc->f(1, i - 1, 1, i - 1, VRNA_DECOMP_EXT_UP, hc->data)) ? eval : 0;

            if (u2 > 0)
              eval = (hc->f(j + 1, n, j + 1, n, VRNA_DECOMP_EXT_UP, hc->data)) ? eval : 0;
          }

          if (eval) {
            qbt1 = q_g *
                   scale[u1 + u2];

            if (md->circ_penalty)
              qbt1 *= pow(vrna_exp_E_hairpin(u1 + u2, 0, -1, -1, NULL, pf_params), (double)n_seq);

            /* apply soft constraints, if any */
            if (sc_ext_wrapper.red_up)
              qbt1 *= sc_ext_wrapper.red_up(1, i - 1, &sc_ext_wrapper) *
                      sc_ext_wrapper.red_up(j + 1, n, &sc_ext_wrapper);

            qho += qbt1;
          }
        }
      }

    /* internal loop cases with at least one gquad below */
    for (i = 1; i < MAXLOOP + 1; i++) {
      u1 = i - 1;
      if (hc->up_int[1] + 1 >= i) {
        /* [gquad] + (basepair) */
        for (j = i + VRNA_GQUAD_MIN_BOX_SIZE - 1; j + turn + 2 <= n; j++) {
#ifndef VRNA_DISABLE_C11_FEATURES
          q_g = vrna_smx_csr_get(q_gq, i, j, 0.);
#else
          q_g = vrna_smx_csr_FLT_OR_DBL_get(q_gq, i, j, 0.);
#endif
          if (q_g != 0.) {
            for (k = j + 1; k + u1 <= j + MAXLOOP + 1; k++) {
              u2 = k - j - 1;
              unsigned int u3, us3;
              unsigned int stop = k + turn + 1;
              
              if (stop + MAXLOOP < n + u1 + u2)
                stop = n + u1 + u2 - MAXLOOP;

              if (hc->up_int[j + 1] >= u2) {
                for (u3 = 0, l = n; l >= stop; l--, u3++) {
                  if (((u2 == 0) && (u1 + u3 < 3)) ||
                      ((u1 + u3 == 0) && (u2 < 3)))
                    continue;


                  eval = (hc->up_int[l] >= u3) ? 1 : 0;
                  eval = (hc->mx[n * k + l] & (VRNA_CONSTRAINT_CONTEXT_INT_LOOP | VRNA_CONSTRAINT_CONTEXT_INT_LOOP_ENC)) ? eval : 0;
                  if (eval) {
                    int kl = my_iindx[k] - l;
                    qbt1 = q_g *
                           qb[kl] *
                           scale[u1 + u2 + u3];

                    if (qbt1 != 0.) {
                      switch (fc->type) {
                        case VRNA_FC_TYPE_SINGLE:
                          sl = S1[l + 1];
                          sk = S1[k - 1];
                          tt    = vrna_get_ptype_md(S[l], S[k], md);
                          if (md->dangles == 2)
                            qbt1 *= pf_params->expmismatchI[tt][sl][sk];

                          if (tt > 2)
                            qbt1 *= pf_params->expTermAU;

                          qbt1 *= (FLT_OR_DBL)expintern[u1 + u2 + u3];
                          break;

                        case VRNA_FC_TYPE_COMPARATIVE:
                          for (s = 0; s < n_seq; s++) {
                            tt = vrna_get_ptype_md(SS[s][l], SS[s][k], md);
                            if (md->dangles == 2)
                              qbt1 *= pf_params->expmismatchI[tt][S3[s][l]][S5[s][k]];

                            if (tt > 2)
                              qbt1 *= pf_params->expTermAU;

                            us1   = (i > 1) ? a2s[s][i - 1] - a2s[s][1] : 0;
                            us2   = a2s[s][k - 1] - a2s[s][j];
                            us3   = a2s[s][n] - a2s[s][l];
                            qbt1  *= (FLT_OR_DBL)expintern[us1 + us2 + us3];
                          }
                          break;
                      }

                      qio += qbt1;
                    }
                  }
                }
              }
            }
          }
        }

        /* [gquad] + [gquad] */
        for (j = i + VRNA_GQUAD_MIN_BOX_SIZE - 1; j + VRNA_GQUAD_MIN_BOX_SIZE <= n; j++) {
#ifndef VRNA_DISABLE_C11_FEATURES
          q_g = vrna_smx_csr_get(q_gq, i, j, 0.);
#else
          q_g = vrna_smx_csr_FLT_OR_DBL_get(q_gq, i, j, 0.);
#endif
          if (q_g != 0.) {
            for (k = j + 1; k + VRNA_GQUAD_MIN_BOX_SIZE - 1 <= n; k++) {
              u2 = k - j - 1;
              if (hc->up_int[j + 1] >= u2) {
                unsigned int stop = k + VRNA_GQUAD_MIN_BOX_SIZE - 1;
                
                if (stop + MAXLOOP < n + u1 + u2)
                  stop = n + u1 + u2 - MAXLOOP;

                unsigned int u3, us3;
                for (u3 = 0, l = n; l >= stop; l--, u3++) {
                  if (((u2 == 0) && (u1 + u3 < 3)) ||
                      ((u1 + u3 == 0) && (u2 < 3)))
                    continue;

#ifndef VRNA_DISABLE_C11_FEATURES
                  FLT_OR_DBL qbt2 = vrna_smx_csr_get(q_gq, k, l, 0.);
#else
                  FLT_OR_DBL qbt2 = vrna_smx_csr_FLT_OR_DBL_get(q_gq, k, l, 0.);
#endif
                  if (qbt2 != 0.) {
                    switch (fc->type) {
                      case VRNA_FC_TYPE_SINGLE:
                        qbt1 = (FLT_OR_DBL)expintern[u1 + u2 + u3];
                        break;

                      case VRNA_FC_TYPE_COMPARATIVE:
                        for (s = 0; s < n_seq; s++) {
                          us1     = (i > 1) ? a2s[s][i - 1] - a2s[s][1] : 0;
                          us2     = a2s[s][k - 1] - a2s[s][j];
                          us3     = a2s[s][n] - a2s[s][l];
                          qbt1 *= (FLT_OR_DBL)expintern[us1 + us2 + us3];
                        }
                        break;
                    }

                    qio += q_g * qbt1 * qbt2 * scale[u1 + u2 + u3];
                  }
                }
              }
            }
          }
        }

        /* (basepair) + [gquad] */
        for (j = i + turn + 1; j + VRNA_GQUAD_MIN_BOX_SIZE <= n; j++) {
          eval  = (hc->mx[n * i + j] & (VRNA_CONSTRAINT_CONTEXT_INT_LOOP | VRNA_CONSTRAINT_CONTEXT_INT_LOOP_ENC)) ? 1 : 0;
          int ij = my_iindx[i] - j;
          qbt1  = qb[ij];
          FLT_OR_DBL qint = 1.;
          if ((eval) &&
              (qbt1 != 0.)) {
            switch (fc->type) {
              case VRNA_FC_TYPE_SINGLE:
                si = S1[i - 1];
                sj = S1[j + 1];
                tt    = vrna_get_ptype_md(S[j], S[i], md);
                if (md->dangles == 2)
                  qint *= pf_params->expmismatchI[tt][sj][si];

                if (tt > 2)
                  qint *= pf_params->expTermAU;

                break;

              case VRNA_FC_TYPE_COMPARATIVE:
                for (s = 0; s < n_seq; s++) {
                  tt = vrna_get_ptype_md(SS[s][j], SS[s][i], md);
                  if (md->dangles == 2)
                    qint *= pf_params->expmismatchI[tt][S3[s][j]][S5[s][i]];

                  if (tt > 2)
                    qint *= pf_params->expTermAU;
                }
                break;
            }

            for (k = j + 1; k + VRNA_GQUAD_MIN_BOX_SIZE - 1 <= n; k++) {
              u2 = k - j - 1;
              if (hc->up_int[j + 1] < u2)
                break;

              unsigned int stop = k + VRNA_GQUAD_MIN_BOX_SIZE - 1;
              if (stop + MAXLOOP < n + u1 + u2)
                stop = n + u1 + u2 - MAXLOOP;

              unsigned int u3, us3;
              for (u3 = 0, l = n; l >= stop; l--, u3++) {
                if (((u2 == 0) && (u1 + u3 < 3)) ||
                    ((u1 + u3 == 0) && (u2 < 3)))
                  continue;

                eval = (hc->up_int[l] >= u3) ? 1 : 0;
                if (eval) {
#ifndef VRNA_DISABLE_C11_FEATURES
                  q_g = vrna_smx_csr_get(q_gq, k, l, 0.);
#else
                  q_g = vrna_smx_csr_FLT_OR_DBL_get(q_gq, k, l, 0.);
#endif
                  if (q_g != 0.) {
                    switch (fc->type) {
                      case VRNA_FC_TYPE_SINGLE:
                        q_g *= (FLT_OR_DBL)expintern[u1 + u2 + u3];
                        break;

                      case VRNA_FC_TYPE_COMPARATIVE:
                        for (s = 0; s < n_seq; s++) {
                          us1 = (i > 1) ? a2s[s][i - 1] - a2s[s][1] : 0;
                          us2 = a2s[s][l - 1] - a2s[s][j];
                          us3 = a2s[s][n] - a2s[s][l];
                          q_g *= (FLT_OR_DBL)expintern[us1 + us2 + us3];
                        }
                        break;
                    }

                    qio += qbt1 * qint *
                           q_g *
                           scale[u1 + u2 + u3];
                  }
                }
              }
            }
          }
        }
      }
    }

  }

  free_sc_mb_exp(&sc_mb_wrapper);
  free_sc_ext_exp(&sc_ext_wrapper);
  free_sc_int_exp(&sc_int_wrapper);

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
