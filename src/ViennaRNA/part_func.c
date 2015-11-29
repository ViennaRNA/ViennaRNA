/*
                  partiton function for RNA secondary structures

                  Ivo L Hofacker + Ronny Lorenz
                  Vienna RNA package
*/
#include <config.h>
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
PRIVATE void  pf_create_bppm(vrna_fold_compound_t *vc, char *structure);

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
vrna_pf( vrna_fold_compound_t *vc,
              char *structure){

  FLT_OR_DBL        Q;
  double            free_energy;
  int               n;
  vrna_md_t         *md;
  vrna_exp_param_t  *params;
  vrna_mx_pf_t      *matrices;

  params    = vc->exp_params;
  n         = vc->length;
  md        = &(params->model_details);
  matrices  = vc->exp_matrices;


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
    vc->stat_cb(vc, VRNA_STATUS_PF_PRE);


  /* do the linear pf fold and fill all matrices  */
  pf_linear(vc);

  if(md->circ)
    pf_circ(vc); /* do post processing step for circular RNAs */

  /* call user-defined recursion status callback function */
  if(vc->stat_cb)
    vc->stat_cb(vc, VRNA_STATUS_PF_POST);

  /* calculate base pairing probability matrix (bppm)  */
  if(md->compute_bpp){
    pf_create_bppm(vc, structure);

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
    fprintf(stderr, "pf_scale too large\n");
  free_energy = (-log(Q)-n*log(params->pf_scale))*params->kT/1000.0;

#ifdef SUN4
  standard_arithmetic();
#else
#ifdef HP9
  fpsetfastmode(0);
#endif
#endif

  return free_energy;
}

PRIVATE void
pf_linear(vrna_fold_compound_t *vc){

  int n, i,j,k,l, ij, kl, u,u1,u2,d,ii, maxk, minl;
  unsigned char type, type_2, tt;

  FLT_OR_DBL  expMLstem = 0.;

  FLT_OR_DBL  temp, Qmax=0;
  FLT_OR_DBL  qbt1, *tmp, q_temp;
  FLT_OR_DBL  *qqm = NULL, *qqm1 = NULL, *qq = NULL, *qq1 = NULL;
  FLT_OR_DBL  *q, *qb, *qm, *qm1, *G, *q1k, *qln;
  FLT_OR_DBL  *scale;
  FLT_OR_DBL  *expMLbase;
  short       *S1;
  int         *my_iindx, *jindx;
  char        *ptype, *sequence;

  vrna_exp_param_t  *pf_params;
  vrna_mx_pf_t      *matrices = vc->exp_matrices;

  pf_params             = vc->exp_params;
  vrna_md_t       *md = &(pf_params->model_details);
  vrna_hc_t       *hc = vc->hc;
  vrna_sc_t       *sc = vc->sc;

  FLT_OR_DBL  expMLclosing = pf_params->expMLclosing;
  double      max_real;
  int         *rtype;
  int         with_gquad  = md->gquad;
  int         circular    = md->circ;
  int         turn        = md->min_loop_size;

  sequence  = vc->sequence;
  n         = vc->length;
  my_iindx  = vc->iindx;
  jindx     = vc->jindx;
  ptype     = vc->ptype;

  q         = matrices->q;
  qb        = matrices->qb;
  qm        = matrices->qm;
  qm1       = matrices->qm1;
  G         = matrices->G;
  q1k       = matrices->q1k;
  qln       = matrices->qln;
  scale     = matrices->scale;
  expMLbase = matrices->expMLbase;

  S1        = vc->sequence_encoding;

  int         hc_decompose;
  char        *hard_constraints = hc->matrix;
  int         *hc_up_ext        = hc->up_ext;
  int         *hc_up_hp         = hc->up_hp;
  int         *hc_up_int        = hc->up_int;
  int         *hc_up_ml         = hc->up_ml;

  max_real = (sizeof(FLT_OR_DBL) == sizeof(float)) ? FLT_MAX : DBL_MAX;

  rtype       = &(md->rtype[0]);

  /* allocate memory for helper arrays */
  qq        = (FLT_OR_DBL *) vrna_alloc(sizeof(FLT_OR_DBL)*(n+2));
  qq1       = (FLT_OR_DBL *) vrna_alloc(sizeof(FLT_OR_DBL)*(n+2));
  qqm       = (FLT_OR_DBL *) vrna_alloc(sizeof(FLT_OR_DBL)*(n+2));
  qqm1      = (FLT_OR_DBL *) vrna_alloc(sizeof(FLT_OR_DBL)*(n+2));

  /*array initialization ; qb,qm,q
    qb,qm,q (i,j) are stored as ((n+1-i)*(n-i) div 2 + n+1-j */

  if(with_gquad)
    expMLstem = exp_E_MLstem(0, -1, -1, pf_params);


  for (d=0; d<=turn; d++)
    for (i=1; i<=n-d; i++) {
      j=i+d;
      ij = my_iindx[i]-j;
      if(hc_up_ext[i] > d){
        q[ij]=1.0*scale[d+1];

        if(sc){
          if(sc->exp_energy_up)
            q[ij] *= sc->exp_energy_up[i][d+1];
          if(sc->exp_f)
            q[ij] *= sc->exp_f(i, j, i, j, VRNA_DECOMP_EXT_UP, sc->data);
        }
      } else {
        q[ij] = 0.;
      }

      qb[ij]  = qm[ij] = 0.0;
    }

  for (i=0; i<=n; i++)
    qq[i] = qq1[i] = qqm[i] = qqm1[i] = 0;

  for (j = turn + 2; j <= n; j++) {
    for (i = j - turn - 1; i >= 1; i--) {
      /* construction of partition function of segment i,j */
      /* firstly that given i binds j : qb(i,j) */
      u             = j - i - 1;
      ij            = my_iindx[i] - j;
      type          = (unsigned char)ptype[jindx[j] + i];
      hc_decompose  = hard_constraints[jindx[j] + i];
      qbt1          = 0;
      q_temp        = 0.;

      if(hc_decompose){
        /* process hairpin loop(s) */
        qbt1 += vrna_exp_E_hp_loop(vc, i, j);
        /* process interior loop(s) */
        qbt1 += vrna_exp_E_int_loop(vc, i, j);
        /* process multibranch loop(s) */
        qbt1 += vrna_exp_E_mb_loop_fast(vc, i, j, qqm1);
      }
      qb[ij] = qbt1;

      /* construction of qqm matrix containing final stem
         contributions to multiple loop partition function
         from segment i,j */
      qqm[i] = 0.;

      if(hc_up_ml[j]){
        q_temp  =  qqm1[i] * expMLbase[1];

        if(sc){
          if(sc->exp_energy_up)
            q_temp *= sc->exp_energy_up[j][1];

          if(sc->exp_f)
            q_temp *= sc->exp_f(i, j, i, j-1, VRNA_DECOMP_ML_ML, sc->data);
        }

        qqm[i] = q_temp;

      }

      if(hc_decompose & VRNA_CONSTRAINT_CONTEXT_MB_LOOP_ENC){
        qbt1 = qb[ij] * exp_E_MLstem(type, ((i>1) || circular) ? S1[i-1] : -1, ((j<n) || circular) ? S1[j+1] : -1, pf_params);
        if(sc){

          if(sc->exp_f)
            qbt1 *= sc->exp_f(i, j, i, j, VRNA_DECOMP_ML_STEM, sc->data);
        }


        qqm[i] += qbt1;
      }

      if(with_gquad){
        /*include gquads into qqm*/
        qqm[i] += G[my_iindx[i]-j] * expMLstem;
      }

      if (qm1) qm1[jindx[j]+i] = qqm[i]; /* for stochastic backtracking and circfold */

      /*construction of qm matrix containing multiple loop
        partition function contributions from segment i,j */
      temp = 0.0;
      kl = my_iindx[i] - j + 1; /* ii-k=[i,k-1] */
      if(sc && sc->exp_f){
        for (k=j; k>i; k--, kl++){
          q_temp = qm[kl] * qqm[k];
          q_temp *= sc->exp_f(i, j, k-1, k, VRNA_DECOMP_ML_ML_ML, sc->data);
          temp += q_temp;
        }
      } else {
        for (k=j; k>i; k--, kl++){
          temp += qm[kl] * qqm[k];
        }
      }

      maxk = MIN2(i+hc_up_ml[i], j);
      ii = 1; /* length of unpaired stretch */
      if(sc){
        for (k=i+1; k<=maxk; k++, ii++){
          q_temp = expMLbase[ii] * qqm[k];
          if(sc->exp_energy_up)
            q_temp *= sc->exp_energy_up[i][ii];

          if(sc->exp_f)
            q_temp *= sc->exp_f(i, j, k, j, VRNA_DECOMP_ML_ML, sc->data);

          temp += q_temp;
        }
      } else {
        for (k=i+1; k<=maxk; k++, ii++){
          temp += expMLbase[ii] * qqm[k];
        }
      }

      qm[ij] = (temp + qqm[i]);

      /*auxiliary matrix qq for cubic order q calculation below */
      qbt1 = 0.;

      if(hc_up_ext[j]){
        q_temp = qq1[i] * scale[1];

        if(sc){
          if(sc->exp_energy_up)
            q_temp *= sc->exp_energy_up[j][1];

          if(sc->exp_f)
            q_temp *= sc->exp_f(i, j, i, j-1, VRNA_DECOMP_EXT_EXT, sc->data);
        }

        qbt1 += q_temp;
      }

      if(hc_decompose & VRNA_CONSTRAINT_CONTEXT_EXT_LOOP){
        q_temp  = qb[ij]
                  * exp_E_ExtLoop(type, ((i>1) || circular) ? S1[i-1] : -1, ((j<n) || circular) ? S1[j+1] : -1, pf_params);
        if(sc){
          if(sc->exp_f)
            q_temp *= sc->exp_f(i, j, i, j, VRNA_DECOMP_EXT_STEM, sc->data);
        }
        qbt1 += q_temp;
      }

      if(with_gquad){
        qbt1 += G[ij];
      }

      qq[i] = qbt1;

      /*construction of partition function for segment i,j */
      temp = qq[i];

      /* the whole stretch [i,j] is unpaired */
      if(hc_up_ext[i] >= (j-i+1)){
        q_temp = 1.0 * scale[j-i+1];

        if(sc){
          if(sc->exp_energy_up)
            q_temp *= sc->exp_energy_up[i][j-i+1];

          if(sc->exp_f)
            q_temp *= sc->exp_f(i, j, i, j, VRNA_DECOMP_EXT_UP, sc->data);
        }

        temp += q_temp;
      }

      kl = my_iindx[i] - i;
      if(sc && sc->exp_f){
        for (k=i; k<j; k++, kl--){
          q_temp = q[kl] * qq[k+1];
          q_temp *= sc->exp_f(i, j, k, k+1, VRNA_DECOMP_EXT_EXT_EXT, sc->data);
          temp += q_temp;
        }
      } else {
        for (k=i; k<j; k++, kl--){
          temp += q[kl] * qq[k+1];
        }
      }

      q[ij] = temp;
      if (temp>Qmax) {
        Qmax = temp;
        if (Qmax>max_real/10.)
          fprintf(stderr, "Q close to overflow: %d %d %g\n", i,j,temp);
      }
      if (temp>=max_real) {
        PRIVATE char msg[128];
        snprintf(msg, 127, "overflow in pf_fold while calculating q[%d,%d]\n"
                     "use larger pf_scale", i,j);
        vrna_message_error(msg);
      }
    }
    tmp = qq1;  qq1 =qq;  qq =tmp;
    tmp = qqm1; qqm1=qqm; qqm=tmp;

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

  /* clean up */
  free(qq);
  free(qq1);
  free(qqm);
  free(qqm1);

}

/* calculate partition function for circular case */
/* NOTE: this is the postprocessing step ONLY     */
/* You have to call pf_linear first to calculate  */
/* complete circular case!!!                      */
PRIVATE void
pf_circ(vrna_fold_compound_t *vc){

  int u, p, q, k, l;
  int noGUclosure, turn;
  int n;
  char  *sequence;
  int   *my_iindx;
  int   *jindx;
  char  *ptype;
  FLT_OR_DBL  *scale;
  short       *S1;

  vrna_exp_param_t     *pf_params = vc->exp_params;
  FLT_OR_DBL    *qb, *qm, *qm1, *qm2, qo, qho, qio, qmo;
  vrna_mx_pf_t  *matrices;

  sequence  = vc->sequence;
  n         = vc->length;

  matrices = vc->exp_matrices;
  qb  = matrices->qb;
  qm  = matrices->qm;
  qm1 = matrices->qm1;
  qm2 = matrices->qm2;

  my_iindx  = vc->iindx;
  jindx     = vc->jindx;

  ptype     = vc->ptype;
  scale     = matrices->scale;
  S1                = vc->sequence_encoding;


  FLT_OR_DBL  qot;
  FLT_OR_DBL  expMLclosing  = pf_params->expMLclosing;
  int         *rtype;

  noGUclosure = pf_params->model_details.noGUclosure;
  turn        = pf_params->model_details.min_loop_size;
  rtype       = &(pf_params->model_details.rtype[0]);

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
      u = n-q + p-1;
      if (u<turn) continue;
      type = ptype[jindx[q] + p];
      if (!type) continue;
       /* cause we want to calc the exterior loops, we need the reversed pair type from now on  */
      type=rtype[type];

      char loopseq[10];
      if (u<7){
        strcpy(loopseq , sequence+q-1);
        strncat(loopseq, sequence, p);
      }
      qho += (((type==3)||(type==4))&&noGUclosure) ? 0. : qb[my_iindx[p]-q] * exp_E_Hairpin(u, type, S1[q+1], S1[p-1],  loopseq, pf_params) * scale[u];

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

  /* add an additional pf of 1.0 to take the open chain into account too           */
  qo = qho + qio + qmo + 1.0*scale[n];

  matrices->qo    = qo;
  matrices->qho   = qho;
  matrices->qio   = qio;
  matrices->qmo   = qmo;

}

/* calculate base pairing probs */
PUBLIC void
pf_create_bppm( vrna_fold_compound_t *vc,
                char *structure){

  int n, i,j,k,l, ij, kl, ii, u1, u2, ov=0;
  unsigned char type, type_2, tt;
  FLT_OR_DBL  temp, Qmax=0, prm_MLb;
  FLT_OR_DBL  prmt, prmt1;
  FLT_OR_DBL  *tmp;
  FLT_OR_DBL  tmp2;
  FLT_OR_DBL  expMLclosing;
  FLT_OR_DBL  *q, *qb, *qm, *qm1, *G, *probs, *scale, *expMLbase;
  FLT_OR_DBL  *q1k, *qln;
  FLT_OR_DBL  qo;

  unsigned char     *ptype;

  double            max_real;
  int               *rtype, with_gquad;
  char              *sequence;
  short             *S, *S1;
  vrna_hc_t         *hc;
  vrna_sc_t         *sc;
  int               *my_iindx, *jindx;
  int               circular, turn;
  vrna_exp_param_t  *pf_params;
  vrna_mx_pf_t      *matrices;
  vrna_md_t         *md;


  pf_params         = vc->exp_params;
  md                = &(pf_params->model_details);
  S                 = vc->sequence_encoding2;
  S1                = vc->sequence_encoding;
  my_iindx          = vc->iindx;
  jindx             = vc->jindx;
  ptype             = vc->ptype;

  circular          = md->circ;
  with_gquad        = md->gquad;
  turn              = md->min_loop_size;

  hc                = vc->hc;
  sc                = vc->sc;

  matrices          = vc->exp_matrices;

  q                 = matrices->q;
  qb                = matrices->qb;
  qm                = matrices->qm;
  qm1               = matrices->qm1;
  G                 = matrices->G;
  probs             = matrices->probs;
  q1k               = matrices->q1k;
  qln               = matrices->qln;
  scale             = matrices->scale;
  expMLbase         = matrices->expMLbase;
  qo                = matrices->qo;

  FLT_OR_DBL  expMLstem         = (with_gquad) ? exp_E_MLstem(0, -1, -1, pf_params) : 0;
  int         hc_decompose;
  char        *hard_constraints = hc->matrix;
  int         *hc_up_ext        = hc->up_ext;
  int         *hc_up_hp         = hc->up_hp;
  int         *hc_up_int        = hc->up_int;
  int         *hc_up_ml         = hc->up_ml;

  int           corr_size       = 5;
  int           corr_cnt        = 0;
  vrna_plist_t  *bp_correction  = vrna_alloc(sizeof(vrna_plist_t) * corr_size);

  max_real      = (sizeof(FLT_OR_DBL) == sizeof(float)) ? FLT_MAX : DBL_MAX;
  sequence      = vc->sequence;

  if((S != NULL) && (S1 != NULL)){

    expMLclosing  = pf_params->expMLclosing;
    with_gquad    = pf_params->model_details.gquad;
    rtype         = &(pf_params->model_details.rtype[0]);
    n             = S[0];

    FLT_OR_DBL *prm_l  = (FLT_OR_DBL *) vrna_alloc(sizeof(FLT_OR_DBL)*(n+2));
    FLT_OR_DBL *prm_l1 = (FLT_OR_DBL *) vrna_alloc(sizeof(FLT_OR_DBL)*(n+2));
    FLT_OR_DBL *prml   = (FLT_OR_DBL *) vrna_alloc(sizeof(FLT_OR_DBL)*(n+2));

    Qmax=0;

    /* 1. exterior pair i,j and initialization of pr array */
    if(circular){
      for (i=1; i<=n; i++) {
        for (j=i; j<=MIN2(i+turn,n); j++)
          probs[my_iindx[i]-j] = 0;
        for (j=i+turn+1; j<=n; j++) {
          ij = my_iindx[i]-j;
          type = (unsigned char)ptype[jindx[j] + i];
          if (type&&(qb[ij]>0.)) {
            probs[ij] = 1./qo;
            unsigned char rt = rtype[type];

            /* 1.1. Exterior Hairpin Contribution */
            int u = i + n - j -1;
            /* get the loop sequence */
            char loopseq[10];
            if (u<7){
              strcpy(loopseq , sequence+j-1);
              strncat(loopseq, sequence, i);
            }
            tmp2 = exp_E_Hairpin(u, rt, S1[j+1], S1[i-1], loopseq, pf_params) * scale[u];

            /* 1.2. Exterior Interior Loop Contribution                     */
            /* 1.2.1. i,j  delimtis the "left" part of the interior loop    */
            /* (j,i) is "outer pair"                                        */
            for(k=1; k < i-turn-1; k++){
              int ln1, lstart;
              ln1 = k + n - j - 1;
              if(ln1>MAXLOOP) break;
              lstart = ln1+i-1-MAXLOOP;
              if(lstart<k+turn+1) lstart = k + turn + 1;
              for(l=lstart; l < i; l++){
                int ln2, type_2;
                type_2 = (unsigned char)ptype[jindx[l] + k];
                if (type_2==0) continue;
                ln2 = i - l - 1;
                if(ln1+ln2>MAXLOOP) continue;
                tmp2 += qb[my_iindx[k] - l]
                        * exp_E_IntLoop(ln1,
                                        ln2,
                                        rt,
                                        rtype[type_2],
                                        S1[j+1],
                                        S1[i-1],
                                        S1[k-1],
                                        S1[l+1],
                                        pf_params)
                        * scale[ln1 + ln2];
              }
            }
            /* 1.2.2. i,j  delimtis the "right" part of the interior loop  */
            for(k=j+1; k < n-turn; k++){
              int ln1, lstart;
              ln1 = k - j - 1;
              if((ln1 + i - 1)>MAXLOOP) break;
              lstart = ln1+i-1+n-MAXLOOP;
              if(lstart<k+turn+1) lstart = k + turn + 1;
              for(l=lstart; l <= n; l++){
                int ln2, type_2;
                type_2 = (unsigned char)ptype[jindx[l] + k];
                if (type_2==0) continue;
                ln2 = i - 1 + n - l;
                if(ln1+ln2>MAXLOOP) continue;
                tmp2 += qb[my_iindx[k] - l]
                        * exp_E_IntLoop(ln2,
                                        ln1,
                                        rtype[type_2],
                                        rt,
                                        S1[l+1],
                                        S1[k-1],
                                        S1[i-1],
                                        S1[j+1],
                                        pf_params)
                        * scale[ln1 + ln2];
              }
            }
            /* 1.3 Exterior multiloop decomposition */
            /* 1.3.1 Middle part                    */
            if((i>turn+2) && (j<n-turn-1))
              tmp2 += qm[my_iindx[1]-i+1]
                      * qm[my_iindx[j+1]-n]
                      * expMLclosing
                      * exp_E_MLstem(type, S1[i-1], S1[j+1], pf_params);

            /* 1.3.2 Left part  */
            for(k=turn+2; k < i-turn-2; k++)
              tmp2 += qm[my_iindx[1]-k]
                      * qm1[jindx[i-1]+k+1]
                      * expMLbase[n-j]
                      * expMLclosing
                      * exp_E_MLstem(type, S1[i-1], S1[j+1], pf_params);

            /* 1.3.3 Right part */
            for(k=j+turn+2; k < n-turn-1;k++)
              tmp2 += qm[my_iindx[j+1]-k]
                      * qm1[jindx[n]+k+1]
                      * expMLbase[i-1]
                      * expMLclosing
                      * exp_E_MLstem(type, S1[i-1], S1[j+1], pf_params);

            /* all exterior loop decompositions for pair i,j done  */
            probs[ij] *= tmp2;

          }
          else probs[ij] = 0;
        }
      }
    } /* end if(circular)  */
    else {
      for (i=1; i<=n; i++) {
        for (j=i; j<=MIN2(i+turn,n); j++)
          probs[my_iindx[i]-j] = 0.;

        for (j=i+turn+1; j<=n; j++) {
          ij = my_iindx[i]-j;
          if(hard_constraints[jindx[j] + i] & VRNA_CONSTRAINT_CONTEXT_EXT_LOOP){
            type      = (unsigned char)ptype[jindx[j] + i];
            probs[ij] = q1k[i-1]*qln[j+1]/q1k[n];
            probs[ij] *= exp_E_ExtLoop(type, (i>1) ? S1[i-1] : -1, (j<n) ? S1[j+1] : -1, pf_params);
            if(sc){
              if(sc->exp_f){
                probs[ij] *= sc->exp_f(1, n, i, j, VRNA_DECOMP_EXT_STEM_OUTSIDE, sc->data);
              }
            }
          } else
            probs[ij] = 0.;
          
        }
      }
    } /* end if(!circular)  */

    for (l = n; l > turn + 1; l--) {

      /* 2. bonding k,l as substem of 2:loop enclosed by i,j */
      for(k = 1; k < l - turn; k++){
        kl      = my_iindx[k]-l;
        type_2  = (unsigned char)ptype[jindx[l] + k];
        type_2  = rtype[type_2];

        if (qb[kl]==0.) continue;

        if(hard_constraints[jindx[l] + k] & VRNA_CONSTRAINT_CONTEXT_INT_LOOP_ENC){
          for(i = MAX2(1, k - MAXLOOP - 1); i <= k - 1; i++){
            u1 = k - i - 1;
            if(hc_up_int[i+1] < u1) continue;

            for(j = l + 1; j <= MIN2(l + MAXLOOP - k + i + 2, n); j++){
              u2 = j-l-1;
              if(hc_up_int[l+1] < u2) break;

              ij = my_iindx[i] - j;
              if(hard_constraints[jindx[j] + i] & VRNA_CONSTRAINT_CONTEXT_INT_LOOP){
                type = (unsigned char)ptype[jindx[j] + i];
                if(probs[ij] > 0){
                  tmp2 =  probs[ij]
                          * scale[u1 + u2 + 2]
                          * exp_E_IntLoop(u1, u2, type, type_2, S1[i+1], S1[j-1], S1[k-1], S1[l+1], pf_params);

                  if(sc){
                    if(sc->exp_energy_up)
                      tmp2 *=   sc->exp_energy_up[i+1][u1]
                              * sc->exp_energy_up[l+1][u2];

                    if(sc->exp_energy_bp)
                      tmp2 *=   sc->exp_energy_bp[ij];

                    if(sc->exp_energy_stack){
                      if((i+1 == k) && (j-1 == l)){
                        tmp2 *=   sc->exp_energy_stack[i]
                                * sc->exp_energy_stack[k]
                                * sc->exp_energy_stack[l]
                                * sc->exp_energy_stack[j];
                      }
                    }

                    if(sc->exp_f){
                      tmp2 *= sc->exp_f(i, j, k, l, VRNA_DECOMP_PAIR_IL, sc->data);
                      if(sc->bt){ /* store probability correction for auxiliary pairs in interior loop motif */
                        vrna_basepair_t *ptr, *aux_bps;
                        aux_bps = sc->bt(i, j, k, l, VRNA_DECOMP_PAIR_IL, sc->data);
                        for(ptr = aux_bps; ptr && ptr->i != 0; ptr++){
                          bp_correction[corr_cnt].i = ptr->i;
                          bp_correction[corr_cnt].j = ptr->j;
                          bp_correction[corr_cnt++].p = tmp2 * qb[kl];
                          if(corr_cnt == corr_size){
                            corr_size += 5;
                            bp_correction = vrna_realloc(bp_correction, sizeof(vrna_plist_t) * corr_size);
                          }
                        }
                        free(aux_bps);
                      }
                    
                    }
                  }

                  probs[kl] += tmp2;
                }
              }
            }
          }
        }
      }

      if(with_gquad){
        /* 2.5. bonding k,l as gquad enclosed by i,j */
        FLT_OR_DBL *expintern = &(pf_params->expinternal[0]);
        FLT_OR_DBL qe;

        if(l < n - 3){
          for(k = 2; k <= l - VRNA_GQUAD_MIN_BOX_SIZE; k++){
            kl = my_iindx[k]-l;
            if (G[kl]==0.) continue;
            tmp2 = 0.;
            i = k - 1;
            for(j = MIN2(l + MAXLOOP + 1, n); j > l + 3; j--){
              ij = my_iindx[i] - j;
              type = (unsigned char)ptype[jindx[j] + i];
              if(!type) continue;
              qe = (type > 2) ? pf_params->expTermAU : 1.;
              tmp2 +=   probs[ij]
                      * qe
                      * expintern[j-l-1]
                      * pf_params->expmismatchI[type][S1[i+1]][S1[j-1]]
                      * scale[2];
            }
            probs[kl] += tmp2 * G[kl];
          }
        }

        if (l < n - 1){
          for (k=3; k<=l-VRNA_GQUAD_MIN_BOX_SIZE; k++) {
            kl = my_iindx[k]-l;
            if (G[kl]==0.) continue;
            tmp2 = 0.;
            for (i=MAX2(1,k-MAXLOOP-1); i<=k-2; i++){
              u1 = k - i - 1;
              for (j=l+2; j<=MIN2(l + MAXLOOP - u1 + 1,n); j++) {
                ij = my_iindx[i] - j;
                type = (unsigned char)ptype[jindx[j] + i];
                if(!type) continue;
                qe = (type > 2) ? pf_params->expTermAU : 1.;
                tmp2 +=   probs[ij]
                        * qe
                        * expintern[u1+j-l-1]
                        * pf_params->expmismatchI[type][S1[i+1]][S1[j-1]]
                        * scale[2];
              }
            }
            probs[kl] += tmp2 * G[kl];
          }
        }

        if(l < n){
          for(k = 4; k <= l - VRNA_GQUAD_MIN_BOX_SIZE; k++){
            kl = my_iindx[k]-l;
            if (G[kl]==0.) continue;
            tmp2 = 0.;
            j = l + 1;
            for (i=MAX2(1,k-MAXLOOP-1); i < k - 3; i++){
              ij = my_iindx[i] - j;
              type = (unsigned char)ptype[jindx[j] + i];
              if(!type) continue;
              qe = (type > 2) ? pf_params->expTermAU : 1.;
              tmp2 +=   probs[ij]
                      * qe
                      * expintern[k - i - 1]
                      * pf_params->expmismatchI[type][S1[i+1]][S1[j-1]]
                      * scale[2];
            }
            probs[kl] += tmp2 * G[kl];
          }
        }
      }

      /* 3. bonding k,l as substem of multi-loop enclosed by i,j */
      prm_MLb = 0.;
      if (l<n)
        for (k = 2; k < l - turn; k++) {
          kl    = my_iindx[k] - l;
          i     = k - 1;
          prmt  = prmt1 = 0.0;

          ii = my_iindx[i];     /* ii-j=[i,j]     */
          tt = (unsigned char)ptype[jindx[l+1] + i];
          tt = rtype[tt];
          if(hard_constraints[jindx[l+1] + i] & VRNA_CONSTRAINT_CONTEXT_MB_LOOP){
            prmt1 = probs[ii-(l+1)]
                    * expMLclosing
                    * exp_E_MLstem(tt, S1[l], S1[i+1], pf_params);

            if(sc){
              /* which decompositions are covered here? => (i, l+1) -> enclosing pair, (k,l) -> enclosed pair, */
              if(sc->exp_energy_bp)
                prmt1 *= sc->exp_energy_bp[ii - (l+1)];

/*
              if(sc->exp_f)
                prmt1 *= sc->exp_f(i, l+1, k, l, , sc->data);
*/
            }
          }
          int lj;
          FLT_OR_DBL ppp;
          ij = my_iindx[i] - (l+2);
          lj = my_iindx[l+1]-(l+1);
          for (j = l + 2; j<=n; j++, ij--, lj--){
            if(hard_constraints[jindx[j] + i] & VRNA_CONSTRAINT_CONTEXT_MB_LOOP){
              tt = (unsigned char)ptype[jindx[j] + i];
              tt = rtype[tt];
              /* which decomposition is covered here? =>
                i + 1 = k < l < j:
                (i,j)       -> enclosing pair
                (k, l)      -> enclosed pair
                (l+1, j-1)  -> multiloop part with at least one stem
                
              */
              ppp = probs[ij]
                    * exp_E_MLstem(tt, S1[j-1], S1[i+1], pf_params)
                    * qm[lj];

              if(sc){
                if(sc->exp_energy_bp)
                  ppp *= sc->exp_energy_bp[ij];
/*
                if(sc->exp_f)
                  ppp *= sc->exp_f(i, j, l+1, j-1, , sc->data);
*/
              }
              prmt += ppp;
            }
          }
          prmt *= expMLclosing;

          tt        =   ptype[jindx[l] + k];

          prml[ i]  =   prmt;

          /* l+1 is unpaired */
          if(hc->up_ml[l+1]){
            ppp = prm_l1[i] * expMLbase[1];
            if(sc){
              if(sc->exp_energy_up)
                ppp *= sc->exp_energy_up[l+1][1];

/*
              if(sc_exp_f)
                ppp *= sc->exp_f(, sc->data);
*/
            }
            prm_l[i] = ppp + prmt1;
          } else {
            prm_l[i] = 0.;
          }

          /* i is unpaired */
          if(hc->up_ml[i]){
            ppp = prm_MLb*expMLbase[1];
            if(sc){
              if(sc->exp_energy_up)
                ppp *= sc->exp_energy_up[i][1];

/*
              if(sc->exp_f)
                ppp *= sc->exp_f(, sc->data);
*/
            }

            prm_MLb = ppp + prml[i];
            /* same as:    prm_MLb = 0;
               for (i=1; i<=k-1; i++) prm_MLb += prml[i]*expMLbase[k-i-1]; */

            prml[i] = prml[ i] + prm_l[i];
          } else {
            prm_MLb = 0.;
          }

          if(with_gquad){
            if ((!tt) && (G[kl] == 0.)) continue;
          } else {
            if (qb[kl] == 0.) continue;
          }

          if(hard_constraints[jindx[l] + k] & VRNA_CONSTRAINT_CONTEXT_MB_LOOP_ENC){

            temp = prm_MLb;

            for (i=1;i<=k-2; i++)
              temp += prml[i]*qm[my_iindx[i+1] - (k-1)];

            if(with_gquad){
              if(tt)
                temp    *= exp_E_MLstem(tt, (k>1) ? S1[k-1] : -1, (l<n) ? S1[l+1] : -1, pf_params) * scale[2];
              else
                temp    *= G[kl] * expMLstem * scale[2];
            } else {
              temp    *= exp_E_MLstem(tt, (k>1) ? S1[k-1] : -1, (l<n) ? S1[l+1] : -1, pf_params) * scale[2];
            }

            probs[kl]  += temp;

            if (probs[kl]>Qmax) {
              Qmax = probs[kl];
              if (Qmax>max_real/10.)
                fprintf(stderr, "P close to overflow: %d %d %g %g\n",
                  i, j, probs[kl], qb[kl]);
            }
            if (probs[kl]>=max_real) {
              ov++;
              probs[kl]=FLT_MAX;
            }
          }
        } /* end for (k=..) */
      tmp = prm_l1; prm_l1=prm_l; prm_l=tmp;

    }  /* end for (l=..)   */

    if(sc && sc->f && sc->bt){
      for (i=1; i<=n; i++)
        for (j=i+turn+1; j<=n; j++) {
          ij = my_iindx[i]-j;
          /*  search for possible auxiliary base pairs in hairpin loop motifs to store
              the corresponding probability corrections
          */ 
          if(hard_constraints[jindx[j] + i] & VRNA_CONSTRAINT_CONTEXT_HP_LOOP){
            vrna_basepair_t *ptr, *aux_bps;
            aux_bps = sc->bt(i, j, i, j, VRNA_DECOMP_PAIR_HP, sc->data);
            if(aux_bps){
              FLT_OR_DBL qhp = vrna_exp_E_hp_loop(vc, i, j);
              for(ptr = aux_bps; ptr && ptr->i != 0; ptr++){
                bp_correction[corr_cnt].i = ptr->i;
                bp_correction[corr_cnt].j = ptr->j;
                bp_correction[corr_cnt++].p = probs[ij] * qhp;
                if(corr_cnt == corr_size){
                  corr_size += 5;
                  bp_correction = vrna_realloc(bp_correction, sizeof(vrna_plist_t) * corr_size);
                }
              }
            }
            free(aux_bps);
          }
        }

      /*  correct pairing probabilities for auxiliary base pairs from hairpin-, or interior loop motifs
          as augmented by the generalized soft constraints feature
      */
      for(i = 0; i < corr_cnt; i++){
        ij = my_iindx[bp_correction[i].i] - bp_correction[i].j;
        /* printf("correcting pair %d, %d by %f\n", bp_correction[i].i, bp_correction[i].j, bp_correction[i].p); */
        probs[ij] += bp_correction[i].p / qb[ij];
      }
    }
    free(bp_correction);

    for (i=1; i<=n; i++)
      for (j=i+turn+1; j<=n; j++) {
        ij = my_iindx[i]-j;

        if(with_gquad){
          probs[ij] *= qb[ij];

          if (G[ij] > 0.){
            probs[ij] += q1k[i-1] * G[ij] * qln[j+1]/q1k[n];
          }
        } else {
          probs[ij] *= qb[ij];
        }
      }

    if (structure!=NULL){
      char *s = vrna_db_from_probs(probs, (unsigned int)n);
      memcpy(structure, s, n);
      structure[n] = '\0';
      free(s);
    }
    if (ov>0) fprintf(stderr, "%d overflows occurred while backtracking;\n"
        "you might try a smaller pf_scale than %g\n",
        ov, pf_params->pf_scale);

    /* clean up */
    free(prm_l);
    free(prm_l1);
    free(prml);

  } /* end if((S != NULL) && (S1 != NULL))  */
  else
    vrna_message_error("bppm calculations have to be done after calling forward recursion\n");
  return;
}

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


PUBLIC double
vrna_mean_bp_distance_pr( int length,
                          FLT_OR_DBL *p){

  int *index = vrna_idx_row_wise((unsigned int) length);
  double d;

  if (p==NULL)
    vrna_message_error("vrna_mean_bp_distance_pr: p==NULL. You need to supply a valid probability matrix");

  d = wrap_mean_bp_distance(p, length, index, TURN);

  free(index);
  return d;
}

PUBLIC double
vrna_mean_bp_distance(vrna_fold_compound_t *vc){

  if(!vc){
    vrna_message_error("vrna_mean_bp_distance: run vrna_pf_fold first!");
  } else if(!vc->exp_matrices){
    vrna_message_error("vrna_mean_bp_distance: exp_matrices==NULL!");
  } else if( !vc->exp_matrices->probs){
    vrna_message_error("vrna_mean_bp_distance: probs==NULL!");
  }

  return wrap_mean_bp_distance( vc->exp_matrices->probs,
                                vc->length,
                                vc->iindx,
                                vc->exp_params->model_details.min_loop_size);
}

PUBLIC vrna_plist_t *
vrna_stack_prob(vrna_fold_compound_t *vc, double cutoff){

  vrna_plist_t             *pl;
  int               i, j, plsize, turn, length, *index, *jindx, *rtype, num;
  char              *ptype;
  FLT_OR_DBL        *qb, *probs, *scale, p;
  vrna_exp_param_t  *pf_params;
  vrna_mx_pf_t      *matrices;

  plsize  = 256;
  pl      = NULL;
  num     = 0;

  if(vc){
    pf_params = vc->exp_params;
    length    = vc->length;
    index     = vc->iindx;
    jindx     = vc->jindx;
    rtype     = &(pf_params->model_details.rtype[0]);
    ptype     = vc->ptype;
    matrices  = vc->exp_matrices;
    qb        = matrices->qb;
    probs     = matrices->probs;
    scale     = matrices->scale;
    turn      = pf_params->model_details.min_loop_size;

    pl        = (vrna_plist_t *) vrna_alloc(plsize*sizeof(vrna_plist_t));

    for (i=1; i<length; i++)
      for (j=i+turn+3; j<=length; j++) {
        if((p=probs[index[i]-j]) < cutoff) continue;
        if (qb[index[i+1]-(j-1)]<FLT_MIN) continue;
        p *= qb[index[i+1]-(j-1)]/qb[index[i]-j];
        p *= exp_E_IntLoop(0,0,(unsigned char)ptype[jindx[j]+i],rtype[(unsigned char)ptype[jindx[j-1] + i+1]],
                           0,0,0,0, pf_params)*scale[2];/* add *scale[u1+u2+2] */
        if (p>cutoff) {
          pl[num].i     = i;
          pl[num].j     = j;
          pl[num].type  = 0;
          pl[num++].p   = p;
          if (num>=plsize) {
            plsize *= 2;
            pl = vrna_realloc(pl, plsize*sizeof(vrna_plist_t));
          }
        }
      }
    pl[num].i=0;
  }

  return pl;
}


/*###########################################*/
/*# deprecated functions below              #*/
/*###########################################*/

#ifdef  VRNA_BACKWARD_COMPAT

PRIVATE float
wrap_pf_fold( const char *sequence,
              char *structure,
              vrna_exp_param_t *parameters,
              int calculate_bppm,
              int is_constrained,
              int is_circular){

  vrna_fold_compound_t  *vc;
  vrna_md_t           md;
  vc                  = NULL;

  /* we need vrna_exp_param_t datastructure to correctly init default hard constraints */
  if(parameters)
    md = parameters->model_details;
  else{
    set_model_details(&md); /* get global default parameters */
  }
  md.circ         = is_circular;
  md.compute_bpp  = calculate_bppm;

  vc = vrna_fold_compound(sequence, &md, VRNA_OPTION_PF);

#if 0
  if(parameters){ /* replace exp_params if necessary */
    free(vc->exp_params);
    vc->exp_params = exp_params;
  } else {
    free(exp_params);
  }
#endif
  

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
  int i,j;
  double d=0;

  if (pr==NULL)
    vrna_message_error("pr==NULL. You need to call pf_fold() before mean_bp_dist()");

  int *my_iindx = vrna_idx_row_wise(length);

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
        int       *my_iindx   = backward_compat_compound->iindx;
        vrna_exp_param_t *pf_params  = backward_compat_compound->exp_params;
        FLT_OR_DBL  *q        = backward_compat_compound->exp_matrices->q;
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
  int       old_circ;

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
