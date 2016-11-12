/*
                  partiton function for RNA secondary structures

                  Ivo L Hofacker + Ronny Lorenz
                  Vienna RNA package
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
PUBLIC  int         st_back = 0;

/*
#################################
# PRIVATE VARIABLES             #
#################################
*/

#ifdef  VRNA_BACKWARD_COMPAT

/* some backward compatibility stuff */
PRIVATE vrna_fold_compound_t  *backward_compat_compound = NULL;
PRIVATE int                 backward_compat           = 0;

#ifdef _OPENMP

#pragma omp threadprivate(backward_compat_compound, backward_compat)

#endif

#endif

/*
#################################
# PRIVATE FUNCTION DECLARATIONS #
#################################
*/
PRIVATE void  pf_circ(vrna_fold_compound_t *vc);
PRIVATE void  pf_linear(vrna_fold_compound_t *vc);
PRIVATE void  alipf_linear(vrna_fold_compound_t *vc);
PRIVATE void  wrap_alipf_circ(vrna_fold_compound_t *vc, char *structure);

#ifdef  VRNA_BACKWARD_COMPAT

PRIVATE float
wrap_pf_fold( const char *sequence,
              char *structure,
              vrna_exp_param_t *parameters,
              int calculate_bppm,
              int is_constrained,
              int is_circular);

#endif

PRIVATE double
wrap_mean_bp_distance(FLT_OR_DBL *p,
                      int length,
                      int *index,
                      int turn);

/*
#################################
# BEGIN OF FUNCTION DEFINITIONS #
#################################
*/

PUBLIC float
vrna_pf_fold( const char *seq,
              char *structure,
              vrna_plist_t **pl){

  float                 free_energy;
  double                mfe;
  vrna_fold_compound_t  *vc;
  vrna_md_t             md;

  vrna_md_set_default(&md);

  /* no need to backtrack MFE structure */
  md.backtrack = 0;

  if(!pl){ /* no need for pair probability computations if we do not store them somewhere */
    md.compute_bpp = 0;
  }

  vc  = vrna_fold_compound(seq, &md, 0);
  mfe = (double)vrna_pf(vc, NULL);
  vrna_exp_params_rescale(vc, &mfe);
  free_energy = vrna_pf(vc, structure);

  /* fill plist */
  if(pl){
    *pl = vrna_plist_from_probs(vc, /*cut_off:*/ 1e-6);
  }

  vrna_fold_compound_free(vc);

  return free_energy;
}

PUBLIC float
vrna_pf_circfold( const char *seq,
                  char *structure,
                  vrna_plist_t **pl){

  float                 free_energy;
  double                mfe;
  vrna_fold_compound_t  *vc;
  vrna_md_t             md;

  vrna_md_set_default(&md);
  md.circ = 1;

  /* no need to backtrack MFE structure */
  md.backtrack = 0;

  if(!pl){ /* no need for pair probability computations if we do not store them somewhere */
    md.compute_bpp = 0;
  }

  vc  = vrna_fold_compound(seq, &md, 0);
  mfe = (double)vrna_mfe(vc, NULL);
  vrna_exp_params_rescale(vc, &mfe);
  free_energy = vrna_pf(vc, structure);

  /* fill plist */
  if(pl){
    *pl = vrna_plist_from_probs(vc, /*cut_off:*/ 1e-6);
  }

  vrna_fold_compound_free(vc);

  return free_energy;
}

PUBLIC float
vrna_pf(vrna_fold_compound_t  *vc,
        char                  *structure){

  int               n;
  FLT_OR_DBL        Q;
  double            free_energy;
  vrna_md_t         *md;
  vrna_exp_param_t  *params;
  vrna_mx_pf_t      *matrices;

  free_energy = (float)(INF/100.);

  if(vc){
    /* make sure, everything is set up properly to start partition function computations */
    if(!vrna_fold_compound_prepare(vc, VRNA_OPTION_PF)){
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
    if(vc->stat_cb)
      vc->stat_cb(VRNA_STATUS_PF_PRE, vc->auxdata);

    switch(vc->type){
      case VRNA_FC_TYPE_SINGLE:     /* do the linear pf fold and fill all matrices  */
                                    pf_linear(vc);

                                    if(md->circ)
                                      pf_circ(vc); /* do post processing step for circular RNAs */

                                    break;

      case VRNA_FC_TYPE_COMPARATIVE:  /* do the linear pf fold and fill all matrices  */
                                    alipf_linear(vc);

                                    /* calculate post processing step for circular  */
                                    /* RNAs                                         */
                                    if(md->circ)
                                      wrap_alipf_circ(vc, structure);

                                    break;

      default:                      vrna_message_warning("vrna_pf@part_func.c: Unrecognized fold compound type");
                                    return free_energy;
                                    break;
    }


    /* call user-defined recursion status callback function */
    if(vc->stat_cb)
      vc->stat_cb(VRNA_STATUS_PF_POST, vc->auxdata);

    /* calculate base pairing probability matrix (bppm)  */
    if(md->compute_bpp){
      vrna_pairing_probs(vc, structure);

#ifdef  VRNA_BACKWARD_COMPAT

      /*
      *  Backward compatibility:
      *  This block may be removed if deprecated functions
      *  relying on the global variable "pr" vanish from within the package!
      */
      pr = matrices->probs;
      /*
       {
        if(pr) free(pr);
        pr = (FLT_OR_DBL *) vrna_alloc(sizeof(FLT_OR_DBL) * ((n+1)*(n+2)/2));
        memcpy(pr, probs, sizeof(FLT_OR_DBL) * ((n+1)*(n+2)/2));
      }
      */

#endif

    }

    if (md->backtrack_type=='C')
      Q = matrices->qb[vc->iindx[1]-n];
    else if (md->backtrack_type=='M')
      Q = matrices->qm[vc->iindx[1]-n];
    else Q = (md->circ) ? matrices->qo : matrices->q[vc->iindx[1]-n];

    /* ensemble free energy in Kcal/mol              */
    if (Q<=FLT_MIN)
      vrna_message_warning("pf_scale too large");

    switch(vc->type){
      case VRNA_FC_TYPE_COMPARATIVE:  free_energy = (-log(Q)-n*log(params->pf_scale))*params->kT/(1000.0 * vc->n_seq);
                                    break;

      case VRNA_FC_TYPE_SINGLE:     /* fall through */

      default:                      free_energy = (-log(Q)-n*log(params->pf_scale))*params->kT/1000.0;
                                    break;
    }

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

PRIVATE void
pf_linear(vrna_fold_compound_t *vc){

  char                *hard_constraints;
  int                 n, i,j, k, ij, d, *my_iindx, *jindx, with_gquad, turn,
                      with_ud, hc_decompose;
  FLT_OR_DBL          temp, Qmax, qbt1, *q, *qb, *qm, *qm1, *q1k, *qln;
  double              max_real;
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
  matrices          = vc->exp_matrices;
  pf_params         = vc->exp_params;
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

  with_ud           = (domains_up && domains_up->exp_energy_cb);
  Qmax              = 0;

  max_real = (sizeof(FLT_OR_DBL) == sizeof(float)) ? FLT_MAX : DBL_MAX;

  if (with_ud && domains_up->exp_prod_cb)
    domains_up->exp_prod_cb(vc, domains_up->data);

  if(with_gquad){
    free(vc->exp_matrices->G);
    vc->exp_matrices->G = get_gquad_pf_matrix(vc->sequence_encoding2, vc->exp_matrices->scale, vc->exp_params);
  }

  /* init auxiliary arrays for fast exterior/multibranch loops */
  aux_mx_el = vrna_exp_E_ext_fast_init(vc);
  aux_mx_ml = vrna_exp_E_ml_fast_init(vc);

  /*array initialization ; qb,qm,q
    qb,qm,q (i,j) are stored as ((n+1-i)*(n-i) div 2 + n+1-j */
  for (d=0; d<=turn; d++)
    for (i=1; i<=n-d; i++) {
      j=i+d;
      ij = my_iindx[i]-j;
      qb[ij] = 0.0;
    }

  for (j = turn + 2; j <= n; j++) {
    for (i = j - turn - 1; i >= 1; i--) {
      /* construction of partition function of segment i,j */
      /* firstly that given i binds j : qb(i,j) */
      ij            = my_iindx[i] - j;
      hc_decompose  = hard_constraints[jindx[j] + i];
      qbt1          = 0;

      if(hc_decompose){
        /* process hairpin loop(s) */
        qbt1 += vrna_exp_E_hp_loop(vc, i, j);
        /* process interior loop(s) */
        qbt1 += vrna_exp_E_int_loop(vc, i, j);
        /* process multibranch loop(s) */
        qbt1 += vrna_exp_E_mb_loop_fast(vc, i, j, aux_mx_ml->qqm1);
      }
      qb[ij] = qbt1;

      /* Multibranch loop */
      qm[ij] = vrna_exp_E_ml_fast(vc, i, j, aux_mx_ml);

      if (qm1)
        qm1[jindx[j] + i] = aux_mx_ml->qqm[i]; /* for stochastic backtracking and circfold */

      /* Exterior loop */
      q[ij] = temp = vrna_exp_E_ext_fast(vc, i, j, aux_mx_el);

      if (temp>Qmax) {
        Qmax = temp;
        if (Qmax>max_real/10.)
          vrna_message_warning("Q close to overflow: %d %d %g", i,j,temp);
      }
      if (temp>=max_real) {
        vrna_message_error("overflow in pf_fold while calculating q[%d,%d]\n"
                                  "use larger pf_scale", i,j);
      }
    }

    /* rotate auxiliary arrays */
    vrna_exp_E_ext_fast_rotate(vc, aux_mx_el);
    vrna_exp_E_ml_fast_rotate(vc, aux_mx_ml);

  }

  /* prefill linear qln, q1k arrays */
  if(q1k && qln){
    for (k=1; k<=n; k++) {
      q1k[k] = q[my_iindx[1] - k];
      qln[k] = q[my_iindx[k] - n];
    }
    q1k[0] = 1.0;
    qln[n+1] = 1.0;
  }

  /* free memory occupied by auxiliary arrays for fast exterior/multibranch loops */
  vrna_exp_E_ml_fast_free(vc, aux_mx_ml);
  vrna_exp_E_ext_fast_free(vc, aux_mx_el);
}

/* calculate partition function for circular case */
/* NOTE: this is the postprocessing step ONLY     */
/* You have to call pf_linear first to calculate  */
/* complete circular case!!!                      */
PRIVATE void
pf_circ(vrna_fold_compound_t *vc){

  char              *ptype;
  short             *S1;
  int               u, p, q, k, l, turn, n, *my_iindx, *jindx, *rtype;
  FLT_OR_DBL        *scale, *qb, *qm, *qm1, *qm2, qo, qho, qio, qmo,
                    qbt1, qot, expMLclosing;
  vrna_exp_param_t  *pf_params;
  vrna_mx_pf_t      *matrices;

  n             = vc->length;
  matrices      = vc->exp_matrices;
  my_iindx      = vc->iindx;
  jindx         = vc->jindx;
  ptype         = vc->ptype;
  S1            = vc->sequence_encoding;
  pf_params     = vc->exp_params;
  qb            = matrices->qb;
  qm            = matrices->qm;
  qm1           = matrices->qm1;
  qm2           = matrices->qm2;
  scale         = matrices->scale;
  expMLclosing  = pf_params->expMLclosing;
  turn          = pf_params->model_details.min_loop_size;
  rtype         = &(pf_params->model_details.rtype[0]);
  qo = qho = qio = qmo = 0.;

  /* construct qm2 matrix from qm1 entries  */
  for(k=1; k<n-turn-1; k++){
    qot = 0.;
    for (u=k+turn+1; u<n-turn-1; u++)
      qot += qm1[jindx[u]+k]*qm1[jindx[n]+(u+1)];
    qm2[k] = qot;
   }

  for(p = 1; p < n; p++){
    for(q = p + turn + 1; q <= n; q++){
      int type;
      /* 1. get exterior hairpin contribution  */
      qbt1 = qb[my_iindx[p]-q] * vrna_exp_E_hp_loop(vc, q, p);
      qho += qbt1;

      u = n-q + p-1;
      if (u<turn) continue;
      type = ptype[jindx[q] + p];
      if (!type) continue;
       /* cause we want to calc the exterior loops, we need the reversed pair type from now on  */
      type=rtype[type];


      /* 2. exterior interior loops, i "define" the (k,l) pair as "outer pair"  */
      /* so "outer type" is rtype[type[k,l]] and inner type is type[p,q]        */
      qot = 0.;
      for(k=q+1; k < n; k++){
        int ln1, lstart;
        ln1 = k - q - 1;
        if(ln1+p-1>MAXLOOP) break;
        lstart = ln1+p-1+n-MAXLOOP;
        if(lstart<k+turn+1) lstart = k + turn + 1;
        for(l=lstart;l <= n; l++){
          int ln2, type2;
          ln2 = (p - 1) + (n - l);

          if((ln1+ln2) > MAXLOOP) continue;

          type2 = ptype[jindx[l] + k];
          if(!type2) continue;
          qio += qb[my_iindx[p]-q] * qb[my_iindx[k]-l] * exp_E_IntLoop(ln2, ln1, rtype[type2], type, S1[l+1], S1[k-1], S1[p-1], S1[q+1], pf_params) * scale[ln1+ln2];
        }
      } /* end of kl double loop */
    }
  } /* end of pq double loop */

  /* 3. Multiloops  */
  for(k=turn+2; k<n-2*turn-3; k++)
    qmo += qm[my_iindx[1]-k] * qm2[k+1] * expMLclosing;

  /* add an additional pf of 1.0 to take the open chain into account too */
  qo = qho + qio + qmo + 1.0*scale[n];

  matrices->qo    = qo;
  matrices->qho   = qho;
  matrices->qio   = qio;
  matrices->qmo   = qmo;

}


PUBLIC int
vrna_pf_float_precision(void){

  return (sizeof(FLT_OR_DBL) == sizeof(float));
}


PRIVATE void
alipf_linear( vrna_fold_compound_t *vc){

  char                *hard_constraints;
  int                 i,j, ij, jij, d, turn, n, *my_iindx, *jindx, *pscore;
  FLT_OR_DBL          temp, Qmax, qbt1, *q, *qb, *qm, *qm1;
  double              kTn, max_real;
  vrna_exp_param_t    *pf_params;
  vrna_mx_pf_t        *matrices;
  vrna_mx_pf_aux_el_t *aux_mx_el;
  vrna_mx_pf_aux_ml_t *aux_mx_ml;
  vrna_md_t           *md;
  vrna_hc_t           *hc;

  n                 = vc->length;
  pf_params         = vc->exp_params;
  matrices          = vc->exp_matrices;
  hc                = vc->hc;
  my_iindx          = vc->iindx;
  jindx             = vc->jindx;
  pscore            = vc->pscore;     /* precomputed array of pair types */
  md                = &(pf_params->model_details);
  q                 = matrices->q;
  qb                = matrices->qb;
  qm                = matrices->qm;
  qm1               = matrices->qm1;
  hard_constraints  = hc->matrix;
  turn              = md->min_loop_size;
  kTn               = pf_params->kT/10.;   /* kT in cal/mol  */
  Qmax              = 0.;

  max_real          = (sizeof(FLT_OR_DBL) == sizeof(float)) ? FLT_MAX : DBL_MAX;

  /* init auxiliary arrays for fast exterior/multibranch loops */
  aux_mx_el = vrna_exp_E_ext_fast_init(vc);
  aux_mx_ml = vrna_exp_E_ml_fast_init(vc);

  /* array initialization ; qb,qm,q
     qb,qm,q (i,j) are stored as ((n+1-i)*(n-i) div 2 + n+1-j */

  for (d = 0; d <= turn; d++)
    for (i = 1; i <= n - d; i++) {
      j       = i + d;
      ij      = my_iindx[i]-j;
      qb[ij]  = 0.0;
    }

  for (j = turn + 2; j <= n; j++) {
    for (i = j - turn - 1; i >= 1; i--) {
      int psc;
      /* construction of partition function for segment i,j */
      /* calculate pf given that i and j pair: qb(i,j)      */
      ij  = my_iindx[i] - j;
      jij = jindx[j] + i;

      psc   = pscore[jij];
      qbt1  = 0.;

      if (hard_constraints[jij]) {
        /* process hairpin loop(s) */
        qbt1 += vrna_exp_E_hp_loop(vc, i, j);
        /* process interior loop(s) */
        qbt1 += vrna_exp_E_int_loop(vc, i, j);
        /* process multibranch loop(s) */
        qbt1 += vrna_exp_E_mb_loop_fast(vc, i, j, aux_mx_ml->qqm1);

        qbt1 *= exp(psc/kTn);
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
        if (Qmax > max_real/10.)
          vrna_message_warning("Q close to overflow: %d %d %g", i,j,temp);
      }
      if (temp >= max_real) {
        vrna_message_error("overflow in pf_fold while calculating q[%d,%d]\n"
                                  "use larger pf_scale", i,j);
      }
    }

    /* rotate auxiliary arrays */
    vrna_exp_E_ext_fast_rotate(vc, aux_mx_el);
    vrna_exp_E_ml_fast_rotate(vc, aux_mx_ml);
  }

  /* free memory occupied by auxiliary arrays for fast exterior/multibranch loops */
  vrna_exp_E_ml_fast_free(vc, aux_mx_ml);
  vrna_exp_E_ext_fast_free(vc, aux_mx_el);
}

/* calculate partition function for circular case   */
/* NOTE: this is the postprocessing step ONLY        */
/* You have to call alipf_linear first to calculate  */
/* circular case!!!                                  */

PRIVATE void
wrap_alipf_circ(vrna_fold_compound_t *vc,
                char *structure){

  char              **Ss, *hard_constraints;
  unsigned short    **a2s;
  short             **S, **S5, **S3;
  int               u, p, q, pq, k, l, s, *type, n_seq, n, *my_iindx, *jindx, *rtype;
  FLT_OR_DBL        qbt1, qot, qo, qho, qio, qmo, *qb, *qm, *qm1, *qm2, *scale,
                    expMLclosing;
  vrna_exp_param_t  *pf_params;
  vrna_mx_pf_t      *matrices;
  vrna_md_t         *md;
  vrna_hc_t         *hc;
  vrna_sc_t         **sc;

  n_seq             = vc->n_seq;
  n                 = vc->length;
  S                 = vc->S;
  S5                = vc->S5;     /* S5[s][i] holds next base 5' of i in sequence s */
  S3                = vc->S3;     /* Sl[s][i] holds next base 3' of i in sequence s */
  Ss                = vc->Ss;
  a2s               = vc->a2s;
  pf_params         = vc->exp_params;
  matrices          = vc->exp_matrices;
  my_iindx          = vc->iindx;
  jindx             = vc->jindx;
  hc                = vc->hc;
  sc                = vc->scs;
  qb                = matrices->qb;
  qm                = matrices->qm;
  qm1               = matrices->qm1;
  qm2               = matrices->qm2;
  scale             = matrices->scale;
  expMLclosing      = pf_params->expMLclosing;
  md                = &(pf_params->model_details);
  hard_constraints  = hc->matrix;
  rtype             = &(md->rtype[0]);
  qo = qho = qio = qmo = 0.;

  type  = (int *)vrna_alloc(sizeof(int) * n_seq);

  /* calculate the qm2 matrix  */
  for(k=1; k<n-TURN; k++){
    qot = 0.;
    for (u=k+TURN+1; u<n-TURN-1; u++)
      qot += qm1[jindx[u]+k]*qm1[jindx[n]+(u+1)];
    qm2[k] = qot;
  }

  for(p=1;p<n;p++){
    for(q=p+TURN+1;q<=n;q++){
      u = n-q + p-1;
      if (u<TURN) continue;
      pq  = jindx[q] + p;

      if(!hard_constraints[pq]) continue;

      for(s = 0; s < n_seq; s++){
        type[s] = md->pair[S[s][p]][S[s][q]];
        if (type[s]==0) type[s]=7;
      }

      /* 1. exterior hairpin contribution  */
      /* Note, that we do not scale Hairpin Energy by u+2 but by u cause the scale  */
      /* for the closing pair was already done in the forward recursion              */
      if(hard_constraints[pq] & VRNA_CONSTRAINT_CONTEXT_HP_LOOP){
        if(hc->up_hp[q+1] > u){
          for (qbt1=1,s=0; s<n_seq; s++) {
            int rt;
            char loopseq[10];
            u   = a2s[s][n] - a2s[s][q] + a2s[s][p] - 1;
            rt  = rtype[type[s]];

            if (u<9){
              strcpy(loopseq , Ss[s] + a2s[s][q] - 1);
              strncat(loopseq, Ss[s], a2s[s][p]);
            }
            qbt1 *= exp_E_Hairpin(u, rt, S3[s][q], S5[s][p], loopseq, pf_params);
          }
          if(sc)
            for(s = 0; s < n_seq; s++){
              if(sc[s]){
                if(sc[s]->exp_energy_up){
                  qbt1 *=   ((p > 1) ? sc[s]->exp_energy_up[1][a2s[s][p]-1] : 1.)
                          * ((q < n) ? sc[s]->exp_energy_up[a2s[s][q]+1][a2s[s][n] - a2s[s][q]] : 1.);
                }
              }
            }
          qho += qb[my_iindx[p]-q] * qbt1 * scale[u];
        }
      }
      /* 2. exterior interior loop contribution*/

      if(hard_constraints[pq] & VRNA_CONSTRAINT_CONTEXT_INT_LOOP){
        for(k=q+1; k < n; k++){
          int ln1, lstart;
          ln1 = k - q - 1;
          if(ln1 + p - 1 > MAXLOOP)
            break;
          if(hc->up_int[q+1] < ln1)
            break;

          lstart = ln1+p-1+n-MAXLOOP;
          if(lstart<k+TURN+1) lstart = k + TURN + 1;
          for(l=lstart;l <= n; l++){
            int ln2, type_2;

            ln2 = (p - 1) + (n - l);

            if(!(hard_constraints[jindx[l]+k] & VRNA_CONSTRAINT_CONTEXT_INT_LOOP))
              continue;
            if((ln1+ln2) > MAXLOOP)
              continue;
            if(hc->up_int[l+1] < ln2)
              continue;

            FLT_OR_DBL qloop=1.;
            if (qb[my_iindx[k]-l]==0.){ qloop=0.; continue;}

            for (s=0; s<n_seq; s++){
              int ln1a = a2s[s][k] - 1 - a2s[s][q];
              int ln2a = a2s[s][n] - a2s[s][l] + a2s[s][p] - 1;
              int rt = rtype[type[s]];
              type_2 = md->pair[S[s][l]][S[s][k]];
              if (type_2 == 0) type_2 = 7;
              qloop *= exp_E_IntLoop(ln1a, ln2a, rt, type_2, S3[s][q], S5[s][p], S5[s][k], S3[s][l], pf_params);
            }
            if(sc)
              for(s = 0; s < n_seq; s++){
                int ln1a = a2s[s][k] - 1 - a2s[s][q];
                int ln2a = a2s[s][n] - a2s[s][l] + a2s[s][p] - 1;
                if(sc[s]){
                  if((ln1a+ln2a == 0) && (sc[s]->exp_energy_stack)){
                    if(S[s][p] && S[s][q] && S[s][k] && S[s][l]){ /* don't allow gaps in stack */
                      qloop *=    sc[s]->exp_energy_stack[a2s[s][p]]
                                * sc[s]->exp_energy_stack[a2s[s][q]]
                                * sc[s]->exp_energy_stack[a2s[s][k]]
                                * sc[s]->exp_energy_stack[a2s[s][l]];
                    }
                  }
                  if(sc[s]->exp_energy_up)
                    qloop *=    sc[s]->exp_energy_up[a2s[s][q] + 1][ln1a]
                              * ((l < n) ? sc[s]->exp_energy_up[a2s[s][l]+1][a2s[s][n] - a2s[s][l]] : 1.)
                              * ((p > 1) ? sc[s]->exp_energy_up[1][a2s[s][p]-1] : 1.);
                }
              }

            qio += qb[my_iindx[p]-q] * qb[my_iindx[k]-l] * qloop * scale[ln1+ln2];
          }
        } /* end of kl double loop */
      }
    }
  } /* end of pq double loop */

  /* 3. exterior multiloop contribution  */
  for(k=TURN+2; k<n-2*TURN-3; k++)
    qmo += qm[my_iindx[1]-k] * qm2[k+1] * pow(expMLclosing,n_seq);

  /* add additional pf of 1.0 to take open chain into account */
  qo = qho + qio + qmo;
  if(hc->up_ext[1] >= n)
     qo += 1.0 * scale[n];

  matrices->qo    = qo;
  matrices->qho   = qho;
  matrices->qio   = qio;
  matrices->qmo   = qmo;

  free(type);
}


/*###########################################*/
/*# deprecated functions below              #*/
/*###########################################*/

#ifdef  VRNA_BACKWARD_COMPAT

PRIVATE double
wrap_mean_bp_distance(FLT_OR_DBL *p,
                      int length,
                      int *index,
                      int turn){

  int         i,j;
  double      d = 0.;

  /* compute the mean base pair distance in the thermodynamic ensemble */
  /* <d> = \sum_{a,b} p_a p_b d(S_a,S_b)
     this can be computed from the pair probs p_ij as
     <d> = \sum_{ij} p_{ij}(1-p_{ij}) */

  for (i=1; i<=length; i++)
    for (j=i+turn+1; j<=length; j++)
      d += p[index[i]-j] * (1-p[index[i]-j]);

  return 2*d;
}

PRIVATE float
wrap_pf_fold( const char *sequence,
              char *structure,
              vrna_exp_param_t *parameters,
              int calculate_bppm,
              int is_constrained,
              int is_circular){

  vrna_fold_compound_t  *vc;
  vrna_md_t             md;

  vc = NULL;

  /* we need vrna_exp_param_t datastructure to correctly init default hard constraints */
  if(parameters)
    md = parameters->model_details;
  else{
    set_model_details(&md); /* get global default parameters */
  }
  md.circ         = is_circular;
  md.compute_bpp  = calculate_bppm;

  vc = vrna_fold_compound(sequence, &md, VRNA_OPTION_DEFAULT);

  /* prepare exp_params and set global pf_scale */
  vc->exp_params = vrna_exp_params(&md);
  vc->exp_params->pf_scale = pf_scale;

  if(is_constrained && structure){
    unsigned int constraint_options = 0;
    constraint_options |= VRNA_CONSTRAINT_DB
                          | VRNA_CONSTRAINT_DB_PIPE
                          | VRNA_CONSTRAINT_DB_DOT
                          | VRNA_CONSTRAINT_DB_X
                          | VRNA_CONSTRAINT_DB_ANG_BRACK
                          | VRNA_CONSTRAINT_DB_RND_BRACK;

    vrna_constraints_add(vc, (const char *)structure, constraint_options);
  }

  if(backward_compat_compound && backward_compat)
    vrna_fold_compound_free(backward_compat_compound);

  backward_compat_compound  = vc;
  backward_compat           = 1;
  iindx = backward_compat_compound->iindx;

  return vrna_pf(vc, structure);
}

PUBLIC vrna_plist_t *
stackProb(double cutoff){

  if(!(backward_compat_compound && backward_compat)){
    vrna_message_error("stackProb: run pf_fold() first!");
  } else if( !backward_compat_compound->exp_matrices->probs){
    vrna_message_error("stackProb: probs==NULL!");
  }

  return vrna_stack_prob(backward_compat_compound, cutoff);
}

PUBLIC char *
centroid( int length,
          double *dist) {

  if (pr==NULL)
    vrna_message_error("pr==NULL. You need to call pf_fold() before centroid()");

  return vrna_centroid_from_probs(length, dist, pr);
}


PUBLIC double
mean_bp_dist(int length) {

  /* compute the mean base pair distance in the thermodynamic ensemble */
  /* <d> = \sum_{a,b} p_a p_b d(S_a,S_b)
     this can be computed from the pair probs p_ij as
     <d> = \sum_{ij} p_{ij}(1-p_{ij}) */

  int     i, j, *my_iindx;
  double  d = 0;

  if (pr==NULL)
    vrna_message_error("pr==NULL. You need to call pf_fold() before mean_bp_dist()");

  my_iindx = vrna_idx_row_wise(length);

  for (i=1; i<=length; i++)
    for (j=i+TURN+1; j<=length; j++)
      d += pr[my_iindx[i]-j] * (1-pr[my_iindx[i]-j]);

  free(my_iindx);
  return 2*d;
}

/* get the free energy of a subsequence from the q[] array */
PUBLIC double
get_subseq_F( int i,
              int j){

  if(backward_compat_compound)
    if(backward_compat_compound->exp_matrices)
      if(backward_compat_compound->exp_matrices->q){
        int               *my_iindx   = backward_compat_compound->iindx;
        vrna_exp_param_t  *pf_params  = backward_compat_compound->exp_params;
        FLT_OR_DBL        *q          = backward_compat_compound->exp_matrices->q;
        return ((-log(q[my_iindx[i]-j])-(j-i+1)*log(pf_params->pf_scale))*pf_params->kT/1000.0);
      }

  vrna_message_error("call pf_fold() to fill q[] array before calling get_subseq_F()");
  return 0.; /* we will never get to this point */
}



/*----------------------------------------------------------------------*/
PUBLIC double
expHairpinEnergy( int u,
                  int type,
                  short si1,
                  short sj1,
                  const char *string) {

/* compute Boltzmann weight of a hairpin loop, multiply by scale[u+2] */

  vrna_exp_param_t *pf_params = backward_compat_compound->exp_params;

  double q, kT;
  kT = pf_params->kT;   /* kT in cal/mol  */
  if(u <= 30)
    q = pf_params->exphairpin[u];
  else
    q = pf_params->exphairpin[30] * exp( -(pf_params->lxc*log( u/30.))*10./kT);
  if ((tetra_loop)&&(u==4)) {
    char tl[7]={0}, *ts;
    strncpy(tl, string, 6);
    if ((ts=strstr(pf_params->Tetraloops, tl)))
      return (pf_params->exptetra[(ts-pf_params->Tetraloops)/7]);
  }
  if ((tetra_loop)&&(u==6)) {
    char tl[9]={0}, *ts;
    strncpy(tl, string, 6);
    if ((ts=strstr(pf_params->Hexaloops, tl)))
      return  (pf_params->exphex[(ts-pf_params->Hexaloops)/9]);
  }
  if (u==3) {
    char tl[6]={0}, *ts;
    strncpy(tl, string, 5);
    if ((ts=strstr(pf_params->Triloops, tl)))
      return (pf_params->exptri[(ts-pf_params->Triloops)/6]);
    if (type>2)
      q *= pf_params->expTermAU;
  }
  else /* no mismatches for tri-loops */
    q *= pf_params->expmismatchH[type][si1][sj1];

  return q;
}

PUBLIC double
expLoopEnergy(int u1,
              int u2,
              int type,
              int type2,
              short si1,
              short sj1,
              short sp1,
              short sq1) {

/* compute Boltzmann weight of interior loop,
   multiply by scale[u1+u2+2] for scaling */
  double z=0;
  int no_close = 0;
  vrna_exp_param_t *pf_params = backward_compat_compound->exp_params;


  if ((no_closingGU) && ((type2==3)||(type2==4)||(type==2)||(type==4)))
    no_close = 1;

  if ((u1==0) && (u2==0)) /* stack */
    z = pf_params->expstack[type][type2];
  else if (no_close==0) {
    if ((u1==0)||(u2==0)) { /* bulge */
      int u;
      u = (u1==0)?u2:u1;
      z = pf_params->expbulge[u];
      if (u2+u1==1) z *= pf_params->expstack[type][type2];
      else {
        if (type>2) z *= pf_params->expTermAU;
        if (type2>2) z *= pf_params->expTermAU;
      }
    }
    else {     /* interior loop */
      if (u1+u2==2) /* size 2 is special */
        z = pf_params->expint11[type][type2][si1][sj1];
      else if ((u1==1) && (u2==2))
        z = pf_params->expint21[type][type2][si1][sq1][sj1];
      else if ((u1==2) && (u2==1))
        z = pf_params->expint21[type2][type][sq1][si1][sp1];
      else if ((u1==2) && (u2==2))
        z = pf_params->expint22[type][type2][si1][sp1][sq1][sj1];
      else if (((u1==2)&&(u2==3))||((u1==3)&&(u2==2))){ /*2-3 is special*/
        z = pf_params->expinternal[5]*
          pf_params->expmismatch23I[type][si1][sj1]*
          pf_params->expmismatch23I[type2][sq1][sp1];
        z *= pf_params->expninio[2][1];
      }
      else if ((u1==1)||(u2==1)) {  /*1-n is special*/
        z = pf_params->expinternal[u1+u2]*
          pf_params->expmismatch1nI[type][si1][sj1]*
          pf_params->expmismatch1nI[type2][sq1][sp1];
        z *= pf_params->expninio[2][abs(u1-u2)];
      }
      else {
        z = pf_params->expinternal[u1+u2]*
          pf_params->expmismatchI[type][si1][sj1]*
          pf_params->expmismatchI[type2][sq1][sp1];
        z *= pf_params->expninio[2][abs(u1-u2)];
      }
    }
  }
  return z;
}

PUBLIC void
init_pf_circ_fold(int length){
/* DO NOTHING */
}

PUBLIC void
init_pf_fold(int length){
/* DO NOTHING */
}

/**
*** Allocate memory for all matrices and other stuff
**/
PUBLIC void
free_pf_arrays(void){

  if(backward_compat_compound && backward_compat){
    vrna_fold_compound_free(backward_compat_compound);
    backward_compat_compound  = NULL;
    backward_compat           = 0;
    iindx = NULL;
  }
}

PUBLIC FLT_OR_DBL *
export_bppm(void){

  if(backward_compat_compound)
    if(backward_compat_compound->exp_matrices)
      if(backward_compat_compound->exp_matrices->probs)
        return backward_compat_compound->exp_matrices->probs;

  return NULL;
}

/*-------------------------------------------------------------------------*/
/* make arrays used for pf_fold available to other routines */
PUBLIC int
get_pf_arrays(short **S_p,
              short **S1_p,
              char **ptype_p,
              FLT_OR_DBL **qb_p,
              FLT_OR_DBL **qm_p,
              FLT_OR_DBL **q1k_p,
              FLT_OR_DBL **qln_p){

  if(backward_compat_compound){
    if(backward_compat_compound->exp_matrices)
      if(backward_compat_compound->exp_matrices->qb){
        *S_p      = backward_compat_compound->sequence_encoding2;
        *S1_p     = backward_compat_compound->sequence_encoding;
        *ptype_p  = backward_compat_compound->ptype_pf_compat;
        *qb_p     = backward_compat_compound->exp_matrices->qb;
        *qm_p     = backward_compat_compound->exp_matrices->qm;
        *q1k_p    = backward_compat_compound->exp_matrices->q1k;
        *qln_p    = backward_compat_compound->exp_matrices->qln;
        return 1;
      }
  }
  return 0;
}

/*-----------------------------------------------------------------*/
PUBLIC float
pf_fold(const char *sequence,
        char *structure){

  return wrap_pf_fold(sequence, structure, NULL, do_backtrack, fold_constrained, 0);
}

PUBLIC float
pf_circ_fold( const char *sequence,
              char *structure){

  return wrap_pf_fold(sequence, structure, NULL, do_backtrack, fold_constrained, 1);
}

PUBLIC float
pf_fold_par(const char *sequence,
            char *structure,
            vrna_exp_param_t *parameters,
            int calculate_bppm,
            int is_constrained,
            int is_circular){

  return wrap_pf_fold(sequence, structure, parameters, calculate_bppm, is_constrained, is_circular);
}

PUBLIC char *
pbacktrack(char *seq){

  int n = (int)strlen(seq);
  return vrna_pbacktrack5(backward_compat_compound, n);
}

PUBLIC char *
pbacktrack5(char *seq,
            int length){

  /* the seq parameter must no differ to the one stored globally anyway, so we just ignore it */
  return vrna_pbacktrack5(backward_compat_compound, length);
}

PUBLIC char *
pbacktrack_circ(char *seq){

  char      *structure;
  vrna_md_t *md;

  structure = NULL;

  if(backward_compat_compound){
    md = &(backward_compat_compound->exp_params->model_details);
    if(md->circ && backward_compat_compound->exp_matrices->qm2){
      structure = vrna_pbacktrack(backward_compat_compound);
    }
  }

  return structure;
}

PUBLIC void
update_pf_params(int length){

  if(backward_compat_compound && backward_compat){
    vrna_md_t         md;
    set_model_details(&md);
    vrna_exp_params_reset(backward_compat_compound, &md);

    /* compatibility with RNAup, may be removed sometime */
    pf_scale = backward_compat_compound->exp_params->pf_scale;
  }
}

PUBLIC void
update_pf_params_par( int length,
                      vrna_exp_param_t *parameters){

  if(backward_compat_compound && backward_compat){
    vrna_md_t         md;
    if(parameters){
      vrna_exp_params_subst(backward_compat_compound, parameters);
    } else {
      set_model_details(&md);
      vrna_exp_params_reset(backward_compat_compound, &md);
    }

    /* compatibility with RNAup, may be removed sometime */
    pf_scale = backward_compat_compound->exp_params->pf_scale;
  }
}

PUBLIC char *
get_centroid_struct_gquad_pr( int length,
                              double *dist){

  return vrna_centroid(backward_compat_compound, dist);
}

PUBLIC void
assign_plist_gquad_from_pr( vrna_plist_t **pl,
                            int length, /* ignored */
                            double cut_off){

  if(!backward_compat_compound){
    *pl = NULL;
  } else if( !backward_compat_compound->exp_matrices->probs){
    *pl = NULL;
  } else {
    *pl = vrna_plist_from_probs(backward_compat_compound, cut_off);
  }
}

PUBLIC double
mean_bp_distance(int length){

  if(backward_compat_compound)
    if(backward_compat_compound->exp_matrices)
      if(backward_compat_compound->exp_matrices->probs)
        return vrna_mean_bp_distance(backward_compat_compound);

  vrna_message_error("mean_bp_distance: you need to call vrna_pf_fold first");
  return 0.; /* we will never get to this point */
}

PUBLIC double
mean_bp_distance_pr(int length,
                    FLT_OR_DBL *p){

  double d=0;
  int *index = vrna_idx_row_wise((unsigned int) length);

  if (p==NULL)
    vrna_message_error("p==NULL. You need to supply a valid probability matrix for mean_bp_distance_pr()");

  d = wrap_mean_bp_distance(p, length, index, TURN);

  free(index);
  return d;
}

#endif
