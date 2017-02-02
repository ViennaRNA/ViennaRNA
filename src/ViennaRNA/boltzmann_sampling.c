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

#include "ViennaRNA/utils.h"
#include "ViennaRNA/energy_par.h"
#include "ViennaRNA/fold_vars.h"
#include "ViennaRNA/loop_energies.h"
#include "ViennaRNA/gquad.h"
#include "ViennaRNA/constraints.h"
#include "ViennaRNA/boltzmann_sampling.h"

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
PRIVATE void  backtrack(int i, int j, char *pstruc, vrna_fold_compound_t *vc);
PRIVATE void  backtrack_qm(int i, int j, char *pstruc, vrna_fold_compound_t *vc);
PRIVATE void  backtrack_qm1(int i,int j, char *pstruc, vrna_fold_compound_t *vc);
PRIVATE void  backtrack_qm2(int u, int n, char *pstruc, vrna_fold_compound_t *vc);
PRIVATE char  *wrap_pbacktrack_circ(vrna_fold_compound_t *vc);

PRIVATE void  backtrack_comparative(vrna_fold_compound_t *vc, char *pstruc, int i, int j, double *prob);
PRIVATE void  backtrack_qm1_comparative(vrna_fold_compound_t *vc, char *pstruc, int i,int j, double *prob);

/*
 *  @brief Sample a consensus secondary structure from the Boltzmann ensemble according its probability
 * 
 *  @ingroup consensus_stochbt
 *
 *  @see vrna_pf() for precomputing the partition function matrices, and
 *
 *  @param  vc    The #vrna_fold_compound_t of type #VRNA_FC_TYPE_COMPARATIVE with precomputed partition function matrices
 *  @param  prob  to be described (berni)
 *  @return       A sampled consensus secondary structure in dot-bracket notation
 */
PRIVATE char *pbacktrack_comparative(vrna_fold_compound_t *vc, double *prob);

/*
#################################
# BEGIN OF FUNCTION DEFINITIONS #
#################################
*/

/*
  stochastic backtracking in pf_fold arrays
  returns random structure S with Boltzman probabilty
  p(S) = exp(-E(S)/kT)/Z
*/
PUBLIC char *
vrna_pbacktrack(vrna_fold_compound_t *vc){

  char    *structure  = NULL;
  double  prob        = 1.;

  if(vc && vc->exp_params){
      switch(vc->type){
        case VRNA_FC_TYPE_SINGLE:     if(vc->exp_params->model_details.circ){
                                        return wrap_pbacktrack_circ(vc);
                                      } else {
                                        return vrna_pbacktrack5(vc, vc->length);
                                      }
                                      break;

        case VRNA_FC_TYPE_COMPARATIVE:  return pbacktrack_comparative(vc, &prob);
                                      break;

        default:                      vrna_message_warning("unrecognized fold compound type");
                                      return structure;
                                      break;
      }
  }

  return structure;
}

PUBLIC char *
vrna_pbacktrack5( vrna_fold_compound_t *vc,
                  int length){

  FLT_OR_DBL        r, qt, q_temp, qkl;
  int               i,j,ij, n, k, u, start, type;
  char              *pstruc;
  int               *my_iindx, *jindx, hc_decompose, *hc_up_ext;
  FLT_OR_DBL        *q, *qb, *q1k, *qln, *scale;
  char              *ptype, *hard_constraints;
  short             *S1;
  vrna_mx_pf_t      *matrices;
  vrna_hc_t         *hc;
  vrna_sc_t         *sc;
  vrna_exp_param_t  *pf_params;

  n         = vc->length;

  pf_params = vc->exp_params;
  my_iindx  = vc->iindx;
  jindx     = vc->jindx;
  matrices  = vc->exp_matrices;

  hc        = vc->hc;
  sc        = vc->sc;
  ptype     = vc->ptype;
  S1        = vc->sequence_encoding;

  q         = matrices->q;
  qb        = matrices->qb;
  q1k       = matrices->q1k;
  qln       = matrices->qln;
  scale     = matrices->scale;

  hard_constraints  = hc->matrix;
  hc_up_ext         = hc->up_ext;

  if(length > n)
    vrna_message_error("part_func.c@pbacktrack5: 3'-end exceeds sequence length");
  else if(length < 1)
    vrna_message_error("part_func.c@pbacktrack5: 3'-end too small");

/*
  if (init_length<1)
    vrna_message_error("can't backtrack without pf arrays.\n"
            "Call pf_fold() before pbacktrack()");
*/

  pstruc = vrna_alloc((length+1)*sizeof(char));

  for (i=0; i<length; i++)
    pstruc[i] = '.';

  if(!(q1k && qln)){
    matrices->q1k = (FLT_OR_DBL *) vrna_alloc(sizeof(FLT_OR_DBL)*(n+1));
    matrices->qln = (FLT_OR_DBL *) vrna_alloc(sizeof(FLT_OR_DBL)*(n+2));
    q1k           = matrices->q1k;
    qln           = matrices->qln;
    for (k=1; k<=n; k++) {
      q1k[k] = q[my_iindx[1] - k];
      qln[k] = q[my_iindx[k] - n];
    }
    q1k[0] = 1.0;
    qln[n+1] = 1.0;
  }


#ifdef VRNA_WITH_BOUSTROPHEDON
  j = length;
  while (j > 1) {
  /* find i position of first pair */
    for (; j>1; j--){
      if(hc_up_ext[j]){
        r = vrna_urn() * q[my_iindx[1] - j];
        q_temp = q[my_iindx[1] - j + 1] * scale[1];

        if(sc){
          if (sc->exp_energy_up)
            q_temp *= sc->exp_energy_up[j][1];

          if(sc->exp_f)
            q_temp *= sc->exp_f(1, j, 1, j-1, VRNA_DECOMP_EXT_EXT, sc->data);
        }

        if (r > q_temp)  break; /* i is paired */
      }
    }
    if (j<=1) break; /* no more pairs */

    /* now find the pairing partner i */
    r = vrna_urn() * (q[my_iindx[1] - j] - q_temp);
    u = j - 1;

    for (qt=0, k=1; k<j; k++) {
      /* apply alternating boustrophedon scheme to variable i */
      i             = (int)(1 + (u - 1)*((k - 1) % 2)) + (int)((1-(2*((k - 1) % 2)))*((k - 1)/2));
      ij            = my_iindx[i]-j;
      type          = ptype[jindx[j] + i];
      hc_decompose  = hard_constraints[jindx[j] + i];
      if (hc_decompose & VRNA_CONSTRAINT_CONTEXT_EXT_LOOP) {

        if(type == 0)
          type = 7;

        qkl = qb[ij] * exp_E_ExtLoop(type, (i>1) ? S1[i-1] : -1, (j<n) ? S1[j+1] : -1, pf_params);

        if (i > 1){
          qkl *= q[my_iindx[1] - i + 1];
          if(sc){
            if(sc->exp_f)
              qkl *= sc->exp_f(1, j, i-1, i, VRNA_DECOMP_EXT_EXT_STEM, sc->data);
          }
        } else {
          if(sc){
            if(sc->exp_f)
              qkl *= sc->exp_f(i, j, i, j, VRNA_DECOMP_EXT_STEM, sc->data);
          }
        }

        qt += qkl;
        if (qt > r) break; /* j is paired */
      }
    }
    if (k==j) vrna_message_error("backtracking failed in ext loop");
    backtrack(i,j, pstruc, vc);
    j = i - 1;
  }
#else
  start = 1;
  while (start<length) {
  /* find i position of first pair */
    for (i=start; i<length; i++) {
      if(hc_up_ext[i]){
        r = vrna_urn() * qln[i];
        q_temp = qln[i+1]*scale[1];

        if(sc){
          if (sc->exp_energy_up)
            q_temp *= sc->exp_energy_up[i][1];

          if(sc->exp_f)
            q_temp *= sc->exp_f(i, length, i+1, length, VRNA_DECOMP_EXT_EXT, sc->data);
        }

        if (r > q_temp)  break; /* i is paired */
      }
    }
    if (i>=length) break; /* no more pairs */
    /* now find the pairing partner j */
    r = vrna_urn() * (qln[i] - q_temp);
    for (qt=0, j=i+1; j<=length; j++) {
      ij            = my_iindx[i]-j;
      type          = ptype[jindx[j] + i];
      hc_decompose  = hard_constraints[jindx[j] + i];
      if (hc_decompose & VRNA_CONSTRAINT_CONTEXT_EXT_LOOP) {

        if(type == 0)
          type = 7;

        qkl = qb[ij] * exp_E_ExtLoop(type, (i>1) ? S1[i-1] : -1, (j<n) ? S1[j+1] : -1, pf_params);

        if (j<length){
          qkl *= qln[j+1];
          if(sc){
            if(sc->exp_f)
              qkl *= sc->exp_f(i, length, j, j+1, VRNA_DECOMP_EXT_STEM_EXT, sc->data);
          }
        } else {
          if(sc){
            if(sc->exp_f)
              qkl *= sc->exp_f(i, j, i, j, VRNA_DECOMP_EXT_STEM, sc->data);
          }
        }

        qt += qkl;
        if (qt > r) break; /* j is paired */
      }
    }
    if (j==length+1) vrna_message_error("backtracking failed in ext loop");
    start = j+1;
    backtrack(i,j, pstruc, vc);
  }
#endif
  return pstruc;
}

PRIVATE void
backtrack_qm( int i,
              int j,
              char *pstruc,
              vrna_fold_compound_t *vc){

  /* divide multiloop into qm and qm1  */
  FLT_OR_DBL        qmt, r, q_temp;
  int               k, n, u, cnt, span, turn;
  FLT_OR_DBL        *qm, *qm1, *expMLbase;
  int               *my_iindx, *jindx, *hc_up_ml;
  vrna_sc_t         *sc;
  vrna_hc_t         *hc;

  n = j;
  vrna_mx_pf_t  *matrices = vc->exp_matrices;

  my_iindx  = vc->iindx;
  jindx     = vc->jindx;

  hc        = vc->hc;
  sc        = vc->sc;
  hc_up_ml  = hc->up_ml;

  qm        = matrices->qm;
  qm1       = matrices->qm1;
  expMLbase = matrices->expMLbase;

  turn      = vc->exp_params->model_details.min_loop_size;

  while(j>i){
    /* now backtrack  [i ... j] in qm[] */
    r   = vrna_urn() * qm[my_iindx[i] - j];
    qmt = qm1[jindx[j]+i];
    k = cnt  = i;
    if(qmt<r)
      for(span = j - i,cnt=i+1; cnt<=j; cnt++){
#ifdef VRNA_WITH_BOUSTROPHEDON
        k = (int)(i + 1 + span * ((cnt - i - 1) % 2)) + (int)((1 - (2 * ((cnt - i - 1) % 2))) * ((cnt - i) / 2));
#else
        k = cnt;
#endif
        q_temp = 0.;
        u = k - i;
        /* [i...k] is unpaired */
        if(hc_up_ml[i] >= u){
          q_temp += expMLbase[u] * qm1[jindx[j]+k];

          if(sc){
            if(sc->exp_energy_up)
              q_temp *= sc->exp_energy_up[i][u];

            if(sc->exp_f)
              q_temp *= sc->exp_f(i, j, k, j, VRNA_DECOMP_ML_ML, sc->data);
          }

          qmt += q_temp;
        }

        /* split between k-1, k */
        q_temp = qm[my_iindx[i]-(k-1)] * qm1[jindx[j]+k];

        if(sc){
          if(sc->exp_f)
            q_temp *= sc->exp_f(i, j, k-1, k, VRNA_DECOMP_ML_ML_ML, sc->data);
        }

        qmt += q_temp;

        if(qmt >= r){ break;}
      }
    if(cnt>j) vrna_message_error("backtrack failed in qm");

    backtrack_qm1(k, j, pstruc, vc);

    if(k<i+turn) break; /* no more pairs */

    u = k - i;
    /* check whether we make the decision to leave [i..k-1] unpaired */
    if(hc_up_ml[i] >= u){
      q_temp = expMLbase[u];

      if(sc){
        if(sc->exp_energy_up)
          q_temp *= sc->exp_energy_up[i][u];

        if(sc->exp_f)
          q_temp *= sc->exp_f(i, k-1, i, k-1, VRNA_DECOMP_ML_UP, sc->data);
      }

      r = vrna_urn() * (qm[my_iindx[i]-(k-1)] + q_temp);
      if(q_temp >= r) break;
    }
    j = k-1;
  }
}

PRIVATE void
backtrack_qm1(int i,
              int j,
              char *pstruc,
              vrna_fold_compound_t *vc){

  /* i is paired to l, i<l<j; backtrack in qm1 to find l */
  int           ii, l, il, type, n, turn;
  FLT_OR_DBL    qt, r, q_temp;
  FLT_OR_DBL    *qm1, *qb, *expMLbase;
  vrna_mx_pf_t  *matrices;
  int           u, *my_iindx, *jindx, *hc_up_ml;
  char          *ptype, *hard_constraints;
  short         *S1;
  vrna_sc_t     *sc;
  vrna_hc_t     *hc;
  vrna_exp_param_t  *pf_params;


  pf_params = vc->exp_params;
  my_iindx  = vc->iindx;
  jindx     = vc->jindx;

  ptype     = vc->ptype;

  sc        = vc->sc;
  hc        = vc->hc;
  hc_up_ml  = hc->up_ml;
  hard_constraints  = hc->matrix;

  matrices  = vc->exp_matrices;
  qb        = matrices->qb;
  qm1       = matrices->qm1;
  expMLbase = matrices->expMLbase;
  S1        = vc->sequence_encoding;

  turn      = pf_params->model_details.min_loop_size;

  n = j;
  r = vrna_urn() * qm1[jindx[j]+i];
  ii = my_iindx[i];
  for (qt=0., l=j; l > i + turn; l--) {
    il = jindx[l] + i;
    if(hard_constraints[il] & VRNA_CONSTRAINT_CONTEXT_MB_LOOP_ENC){
      u = j - l;
      if(hc_up_ml[l+1] >= u){
        type = ptype[il];

        if(type == 0)
          type = 7;

        q_temp =  qb[ii-l]
                  * exp_E_MLstem(type, S1[i-1], S1[l+1], pf_params)
                  * expMLbase[j-l];

        if(sc){
          if(sc->exp_energy_up)
            q_temp *= sc->exp_energy_up[l+1][j-l];

          if(sc->exp_f)
            q_temp *= sc->exp_f(i, j, i, l, VRNA_DECOMP_ML_STEM, sc->data);
        }

        qt += q_temp;
        if (qt>=r) break;
      } else {
        l = i + turn;
        break;
      }
    }
  }
  if (l < i + turn + 1) vrna_message_error("backtrack failed in qm1");
  backtrack(i, l, pstruc, vc);
}

PRIVATE void
backtrack_qm2(int k,
              int n,
              char *pstruc,
              vrna_fold_compound_t *vc){

  FLT_OR_DBL  qom2t, r;
  int         u, turn;
  FLT_OR_DBL  *qm1, *qm2;
  int         *jindx;

  jindx     = vc->jindx;
  qm1       = vc->exp_matrices->qm1;
  qm2       = vc->exp_matrices->qm2;
  turn      = vc->exp_params->model_details.min_loop_size;

  r= vrna_urn()*qm2[k];
  /* we have to search for our barrier u between qm1 and qm1  */
  for (qom2t = 0.,u=k+turn+1; u<n-turn-1; u++){
    qom2t += qm1[jindx[u]+k]*qm1[jindx[n]+(u+1)];
    if(qom2t > r) break;
  }
  if(u==n-turn) vrna_message_error("backtrack failed in qm2");
  backtrack_qm1(k, u, pstruc, vc);
  backtrack_qm1(u+1, n, pstruc, vc);
}

PRIVATE void
backtrack(int i,
          int j,
          char *pstruc,
          vrna_fold_compound_t *vc){

  char              *ptype, *sequence, *hard_constraints, hc_decompose;
  vrna_exp_param_t  *pf_params;
  FLT_OR_DBL        *qb, *qm, *qm1, *scale, tmp;
  FLT_OR_DBL        r, qbt1, qt, q_temp;
  vrna_mx_pf_t      *matrices;
  int               *my_iindx, *jindx, *hc_up_int, *hc_up_hp;
  vrna_sc_t         *sc;
  vrna_hc_t         *hc;
  short             *S1;

  sequence    = vc->sequence;
  pf_params   = vc->exp_params;
  ptype       = vc->ptype;
  S1          = vc->sequence_encoding;
  my_iindx    = vc->iindx;
  jindx       = vc->jindx;

  sc          = vc->sc;
  hc          = vc->hc;
  hc_up_hp    = hc->up_hp;
  hc_up_int   = hc->up_int;
  hard_constraints  = hc->matrix;

  matrices    = vc->exp_matrices;
  qb          = matrices->qb;
  qm          = matrices->qm;
  qm1         = matrices->qm1;
  scale       = matrices->scale;

  int noGUclosure = pf_params->model_details.noGUclosure;
  int turn        = pf_params->model_details.min_loop_size;
  int   *rtype    = &(pf_params->model_details.rtype[0]);
  int n;
  n = j;
  do {
    int k, l, kl, u, u1, u2, max_k, min_l;
    unsigned char type;
    k = i;
    l = j;

    pstruc[i-1] = '('; pstruc[j-1] = ')';

    r = vrna_urn() * qb[my_iindx[i]-j];
    tmp = qb[my_iindx[i]-j];
    type = (unsigned char)ptype[jindx[j] + i];
    hc_decompose = hard_constraints[jindx[j] + i];
    if(hc_decompose & VRNA_CONSTRAINT_CONTEXT_HP_LOOP){ /* hairpin contribution */

      if(type == 0)
        type = 7;

      u = j-i-1;

      if (((type==3)||(type==4))&&noGUclosure) qbt1 = 0;
      else{
        q_temp = exp_E_Hairpin(u, type, S1[i+1], S1[j-1], sequence+i-1, pf_params) * scale[u+2];

        if(sc){
          if(sc->exp_energy_up)
            q_temp *= sc->exp_energy_up[i+1][u];

          if(sc->exp_f)
            q_temp *= sc->exp_f(i, j, i, j, VRNA_DECOMP_PAIR_HP, sc->data);
        }

        qbt1 = q_temp;

      }
      if (qbt1>=r) return; /* found the hairpin we're done */
    }

    if(hc_decompose & VRNA_CONSTRAINT_CONTEXT_INT_LOOP){ /* interior loop contributions */

      if(type == 0)
        type = 7;

      max_k = i + MAXLOOP + 1;
      max_k = MIN2(max_k, j - turn - 2);
      max_k = MIN2(max_k, i + 1 + hc_up_int[i+1]);
      for (k = i + 1; k<=max_k; k++) {
        u1    = k-i-1;
        min_l = MAX2(k+turn+1,j-1-MAXLOOP+u1);
        kl    = my_iindx[k] - j + 1;
        for (u2 = 0, l=j-1; l>=min_l; l--, kl++, u2++){
          if(hc_up_int[l+1] < u2) break;
          if(hard_constraints[jindx[l] + k] & VRNA_CONSTRAINT_CONTEXT_INT_LOOP_ENC){
            unsigned char type_2 = (unsigned char)ptype[jindx[l] + k];
            type_2 = rtype[type_2];

            if(type_2 == 0)
              type_2 = 7;

            /* add *scale[u1+u2+2] */
            q_temp = qb[kl]
                     * scale[u1+u2+2]
                     * exp_E_IntLoop(u1, u2, type, type_2, S1[i+1], S1[j-1], S1[k-1], S1[l+1], pf_params);

            if(sc){
              if(sc->exp_energy_up)
                q_temp *=   sc->exp_energy_up[i+1][u1]
                          * sc->exp_energy_up[l+1][u2];

              if(sc->exp_energy_stack)
                if((i + 1 == k) && (j - 1 == l))
                  q_temp *=   sc->exp_energy_stack[i]
                            * sc->exp_energy_stack[k]
                            * sc->exp_energy_stack[l]
                            * sc->exp_energy_stack[j];

              if(sc->exp_f)
                q_temp *= sc->exp_f(i, j, k, l, VRNA_DECOMP_PAIR_IL, sc->data);
            }

            qbt1 += q_temp;
            if (qbt1 >= r) break;
          }
        }
        if (qbt1 >= r) break;
      }
      if (k <= max_k) {
        i=k; j=l;
      } else { /* interior loop contributions did not exceed threshold, so we break */
        break;
      }
    } else { /* must not be interior loop, so we break out */
      break;
    }
  } while (1);

  /* backtrack in multi-loop */
  {
    int k, ii, jj, tt;
    FLT_OR_DBL closingPair;
    tt = rtype[(unsigned char)ptype[jindx[j] + i]];
    closingPair =   pf_params->expMLclosing
                  * exp_E_MLstem(tt, S1[j-1], S1[i+1], pf_params)
                  * scale[2];
    if(sc){
      if(sc->exp_f)
        closingPair *= sc->exp_f(i, j, i, j, VRNA_DECOMP_PAIR_ML, sc->data);
    }

    i++; j--;
    /* find the first split index */
    ii = my_iindx[i]; /* ii-j=[i,j] */
    jj = jindx[j]; /* jj+i=[j,i] */
    for (qt=qbt1, k=i+1; k<j; k++) {

      q_temp = qm[ii-(k-1)] * qm1[jj+k] * closingPair;

      if(sc){
        if(sc->exp_f)
          q_temp *= sc->exp_f(i, j, k-1, k, VRNA_DECOMP_ML_ML_ML, sc->data);
      }

      qt += q_temp;
      qbt1 += q_temp;
      if (qt>=r) break;
    }
    if (k>=j){
      vrna_message_error("backtrack failed, can't find split index ");
    }

    backtrack_qm1(k, j, pstruc, vc);

    j = k-1;
    backtrack_qm(i, j, pstruc, vc);
  }
}

PRIVATE char *
wrap_pbacktrack_circ(vrna_fold_compound_t *vc){

  FLT_OR_DBL  r, qt;
  int         i, j, k, l, n;
  vrna_exp_param_t   *pf_params;
  FLT_OR_DBL  qo, qmo;
  FLT_OR_DBL  *scale, *qb, *qm, *qm2;
  char        *sequence, *ptype, *pstruc;
  int         *my_iindx, *jindx;
  short       *S1;

  vrna_mx_pf_t *matrices;

  pf_params     = vc->exp_params;
  matrices      = vc->exp_matrices;
  ptype         = vc->ptype;
  my_iindx      = vc->iindx;
  jindx         = vc->jindx;
  S1            = vc->sequence_encoding;

  qo            = matrices->qo;
  qmo           = matrices->qmo;
  qb            = matrices->qb;
  qm            = matrices->qm;
  qm2           = matrices->qm2;
  scale         = matrices->scale;

  FLT_OR_DBL  expMLclosing  = pf_params->expMLclosing;
  int         *rtype        = &(pf_params->model_details.rtype[0]);
  int         turn          = pf_params->model_details.min_loop_size;

  sequence  = vc->sequence;
  n         = vc->length;

/*
  if (init_length<1)
    vrna_message_error("can't backtrack without pf arrays.\n"
      "Call pf_circ_fold() before pbacktrack_circ()");
*/

  pstruc = vrna_alloc((n+1)*sizeof(char));

  /* initialize pstruct with single bases  */
  for (i=0; i<n; i++) pstruc[i] = '.';

  qt = 1.0*scale[n];
  r = vrna_urn() * qo;

  /* open chain? */
  if(qt > r) return pstruc;

  for(i=1; (i < n); i++){
    for(j=i+turn+1;(j<=n); j++){

      int type, u;
      /* 1. first check, wether we can do a hairpin loop  */
      u = n-j + i-1;
      if (u<turn) continue;

      type = ptype[jindx[j] + i];
      if (!type) continue;

      type=rtype[type];

      char loopseq[10];
      if (u<7){
        strcpy(loopseq , sequence+j-1);
        strncat(loopseq, sequence, i);
      }

      qt += qb[my_iindx[i]-j] * exp_E_Hairpin(u, type, S1[j+1], S1[i-1],  loopseq, pf_params) * scale[u];
      /* found a hairpin? so backtrack in the enclosed part and we're done  */
      if(qt>r){ backtrack(i,j, pstruc, vc); return pstruc;}

      /* 2. search for (k,l) with which we can close an interior loop  */
      for(k=j+1; (k < n); k++){
        int ln1, lstart;
        ln1 = k - j - 1;
        if(ln1+i-1>MAXLOOP) break;

        lstart = ln1+i-1+n-MAXLOOP;
        if(lstart<k+turn+1) lstart = k + turn + 1;
        for(l=lstart; (l <= n); l++){
            int ln2, type2;
            ln2 = (i - 1) + (n - l);
            if((ln1+ln2) > MAXLOOP) continue;

            type2 = ptype[jindx[l] + k];
            if(!type) continue;
            type2 = rtype[type2];
            qt += qb[my_iindx[i]-j] * qb[my_iindx[k]-l] * exp_E_IntLoop(ln2, ln1, type2, type, S1[l+1], S1[k-1], S1[i-1], S1[j+1], pf_params) * scale[ln1 + ln2];
            /* found an exterior interior loop? also this time, we can go straight  */
            /* forward and backtracking the both enclosed parts and we're done      */
            if(qt>r){ backtrack(i,j, pstruc, vc); backtrack(k,l, pstruc, vc); return pstruc;}
        }
      } /* end of kl double loop */
    }
  } /* end of ij double loop  */
  {
    /* as we reach this part, we have to search for our barrier between qm and qm2  */
    qt = 0.;
    r = vrna_urn()*qmo;
    for(k=turn+2; k<n-2*turn-3; k++){
      qt += qm[my_iindx[1]-k] * qm2[k+1] * expMLclosing;
      /* backtrack in qm and qm2 if we've found a valid barrier k  */
      if(qt>r){ backtrack_qm(1,k, pstruc, vc); backtrack_qm2(k+1,n, pstruc, vc); return pstruc;}
    }
  }
  /* if we reach the actual end of this function, an error has occured  */
  /* cause we HAVE TO find an exterior loop or an open chain!!!         */
  vrna_message_error("backtracking failed in exterior loop");
  return pstruc;
}


PRIVATE char *
pbacktrack_comparative( vrna_fold_compound_t *vc,
                        double *prob){

  FLT_OR_DBL  r, gr, qt;
  int         k,i,j, start,s;
  FLT_OR_DBL      probs=1;
  char        *pstruc = NULL;

  int               n_seq       = vc->n_seq;
  int               n           = vc->length;
  short             **S         = vc->S;
  short             **S5        = vc->S5;     /*S5[s][i] holds next base 5' of i in sequence s*/
  short             **S3        = vc->S3;     /*Sl[s][i] holds next base 3' of i in sequence s*/
  vrna_exp_param_t  *pf_params  = vc->exp_params;
  vrna_mx_pf_t      *matrices   = vc->exp_matrices;
  vrna_md_t         *md         = &(pf_params->model_details);
  int               *my_iindx   = vc->iindx;
  vrna_hc_t         *hc         = vc->hc;
  vrna_sc_t         **sc        = vc->scs;
  FLT_OR_DBL        *q          = matrices->q;
  FLT_OR_DBL        *qb         = matrices->qb;

  if((matrices->q1k == NULL) || (matrices->qln == NULL)){
    free(matrices->q1k);
    matrices->q1k = (FLT_OR_DBL *) vrna_alloc(sizeof(FLT_OR_DBL)*(n+1));
    free(matrices->qln);
    matrices->qln = (FLT_OR_DBL *) vrna_alloc(sizeof(FLT_OR_DBL)*(n+2));
  }

  FLT_OR_DBL *q1k   = matrices->q1k;
  FLT_OR_DBL *qln   = matrices->qln;
  FLT_OR_DBL *scale = matrices->scale;

  for (k=1; k<=n; k++) {
    q1k[k] = q[my_iindx[1] - k];
    qln[k] = q[my_iindx[k] - n];
  }
  q1k[0] = 1.0;
  qln[n+1] = 1.0;

  pstruc = vrna_alloc((n+1)*sizeof(char));

  for (i=0; i<n; i++)
    pstruc[i] = '.';


  start = 1;
  while (start<n) {
  /* find i position of first pair */
    probs=1.;
    for (i=start; i<n; i++) {
      gr = vrna_urn() * qln[i];
      if (gr > qln[i+1]*scale[1]) {
        *prob=*prob*probs*(1-qln[i+1]*scale[1]/qln[i]);
        break; /* i is paired */
      }
      probs*=qln[i+1]*scale[1]/qln[i];
    }
    if (i>=n) {
      *prob=*prob*probs;
      break; /* no more pairs */
    }
    /* now find the pairing partner j */
    r = vrna_urn() * (qln[i] - qln[i+1]*scale[1]);
    for (qt=0, j=i+1; j<=n; j++) {
      int xtype;
      /*  type = ptype[my_iindx[i]-j];
          if (type) {*/
      FLT_OR_DBL qkl;
      if (qb[my_iindx[i]-j]>0) {
        qkl = qb[my_iindx[i]-j]*qln[j+1];  /*if psc too small qb=0!*/
        for (s=0; s< n_seq; s++) {
          xtype=md->pair[S[s][i]][S[s][j]];
          if (xtype==0) xtype=7;
          qkl *= exp_E_ExtLoop(xtype, (i>1) ? S5[s][i] : -1, (j<n) ? S3[s][j] : -1, pf_params);
        }
        qt += qkl; /*?*exp(pscore[jindx[j]+i]/kTn)*/
        if (qt > r) {
          *prob=*prob*(qkl/(qln[i] - qln[i+1]*scale[1]));/*probs*=qkl;*/
          break; /* j is paired */
        }
      }
    }
    if (j==n+1) vrna_message_error("backtracking failed in ext loop");
    start = j+1;
    backtrack_comparative(vc, pstruc, i, j, prob); /*?*/
  }

  return pstruc;
}


PRIVATE void
backtrack_comparative(vrna_fold_compound_t *vc,
          char *pstruc,
          int i,
          int j,
          double *prob){

  int               n_seq       = vc->n_seq;
  short             **S         = vc->S;
  short             **S5        = vc->S5;     /*S5[s][i] holds next base 5' of i in sequence s*/
  short             **S3        = vc->S3;     /*Sl[s][i] holds next base 3' of i in sequence s*/
  char              **Ss        = vc->Ss;
  unsigned short    **a2s       = vc->a2s;
  vrna_exp_param_t  *pf_params  = vc->exp_params;
  vrna_mx_pf_t      *matrices   = vc->exp_matrices;
  vrna_md_t         *md         = &(pf_params->model_details);
  int               *my_iindx   = vc->iindx;
  int               *jindx      = vc->jindx;
  vrna_hc_t         *hc         = vc->hc;
  vrna_sc_t         **sc        = vc->scs;
  FLT_OR_DBL        *qb         = matrices->qb;
  FLT_OR_DBL        *qm         = matrices->qm;
  FLT_OR_DBL        *qm1        = matrices->qm1;
  int               *pscore     = vc->pscore;     /* precomputed array of pair types */

  FLT_OR_DBL        *scale        = matrices->scale;
  FLT_OR_DBL        *expMLbase    = matrices->expMLbase;

  /*backtrack given i,j basepair!*/
  FLT_OR_DBL kTn = pf_params->kT/10.;
  int *type = (int *)vrna_alloc(sizeof(int) * n_seq);

  do {
    FLT_OR_DBL  r, qbt1, max_k, min_l;
    int         k, l, u, u1, u2, s;
    pstruc[i-1] = '('; pstruc[j-1] = ')';
    for (s=0; s<n_seq; s++) {
      type[s] = md->pair[S[s][i]][S[s][j]];
      if (type[s]==0) type[s]=7;
    }
    r = vrna_urn() * (qb[my_iindx[i]-j]/exp(pscore[jindx[j]+i]/kTn)); /*?*exp(pscore[jindx[j]+i]/kTn)*/

    qbt1=1.;
    for (s=0; s<n_seq; s++){
      u = a2s[s][j-1]-a2s[s][i];
      if (a2s[s][i]<1) continue;
      char loopseq[10];
      if(u < 9){
        strncpy(loopseq, Ss[s]+a2s[s][i]-1, 10);
      }
      qbt1 *= exp_E_Hairpin(u, type[s], S3[s][i], S5[s][j], loopseq, pf_params);
    }
    qbt1 *= scale[j-i+1];

    if (qbt1>r) {
      *prob=*prob*qbt1/(qb[my_iindx[i]-j]/exp(pscore[jindx[j]+i]/kTn));/*probs*=qbt1;*/
      free(type);
      return; /* found the hairpin we're done */
    }


    max_k = MIN2(i+MAXLOOP+1,j-TURN-2);
    l = MAX2(i+TURN+2,j-MAXLOOP-1);
    for (k=i+1; k<=max_k; k++){
      min_l = MAX2(k+TURN+1,j-1-MAXLOOP+k-i-1);

      for (l=min_l; l<j; l++){
        FLT_OR_DBL qloop=1;
        int type_2;
        if (qb[my_iindx[k]-l]==0) {qloop=0; continue;}
        for (s=0; s<n_seq; s++) {
          u1      = a2s[s][k-1] - a2s[s][i]/*??*/;
          u2      = a2s[s][j-1] - a2s[s][l];
          type_2  = md->pair[S[s][l]][S[s][k]];
          if(type_2 == 0) type_2 = 7;

          qloop *= exp_E_IntLoop(u1, u2, type[s], type_2, S3[s][i], S5[s][j],S5[s][k], S3[s][l], pf_params);
        }

        if(sc)
          for (s=0; s<n_seq; s++) {
            if(sc[s]){
              int u1 = a2s[s][k-1] - a2s[s][i];
              int u2 = a2s[s][j-1] - a2s[s][l];
              if(u1 + u2 == 0)
                if(sc[s]->exp_energy_stack){
                  if(S[s][i] && S[s][j] && S[s][k] && S[s][l]){ /* don't allow gaps in stack */
                    qloop *=    sc[s]->exp_energy_stack[i]
                              * sc[s]->exp_energy_stack[k]
                              * sc[s]->exp_energy_stack[l]
                              * sc[s]->exp_energy_stack[j];
                  }
                }
            }
          }

        qbt1 += qb[my_iindx[k]-l] * qloop * scale[k-i+j-l];

        if (qbt1 > r) {
         *prob =  *prob
                  * qb[my_iindx[k]-l]
                  * qloop
                  * scale[k-i+j-l]
                  / (   qb[my_iindx[i]-j]
                      / exp(pscore[jindx[j]+i] / kTn));
         /*
          prob*=qb[my_iindx[k]-l] * qloop * scale[k-i+j-l];
         */
          break;
        }
      }
      if (qbt1 > r) break;
    }
    if (l<j) {
      i=k; j=l;
    }
    else {
      *prob=*prob*(1-qbt1/(qb[my_iindx[i]-j]/exp(pscore[jindx[j]+i]/kTn)));
      break;
    }
  } while (1);

  /* backtrack in multi-loop */
  {
    FLT_OR_DBL r, qt;
    int k, ii, jj;
    FLT_OR_DBL qttemp=0;;
    i++; j--;
    /* find the first split index */
    ii = my_iindx[i]; /* ii-j=[i,j] */
    jj = jindx[j]; /* jj+i=[j,i] */
    for (qt=0., k=i+1; k<j; k++) qttemp += qm[ii-(k-1)]*qm1[jj+k];
    r = vrna_urn() * qttemp;
    for (qt=0., k=i+1; k<j; k++) {
      qt += qm[ii-(k-1)]*qm1[jj+k];
      if (qt>=r){
        *prob = *prob
                * qm[ii-(k-1)]
                * qm1[jj+k]
                / qttemp;/*qttemp;*/
        /*        prob*=qm[ii-(k-1)]*qm1[jj+k];*/
        break;
      }
    }
    if (k>=j) vrna_message_error("backtrack failed, can't find split index ");

    backtrack_qm1_comparative(vc, pstruc, k, j, prob);

    j = k-1;
    while (j>i) {
      /* now backtrack  [i ... j] in qm[] */
      jj = jindx[j];/*habides??*/
      ii = my_iindx[i];
      r = vrna_urn() * qm[ii - j];
      qt = qm1[jj+i]; k=i;
      if (qt<r)
        for (k=i+1; k<=j; k++) {
          qt += (qm[ii-(k-1)]+expMLbase[k-i]/*n_seq??*/)*qm1[jj+k];
          if (qt >= r) {
            *prob = *prob
                    * (qm[ii-(k-1)] + expMLbase[k-i])
                    * qm1[jj+k]
                    / qm[ii - j];/*???*/
            /*            probs*=qt;*/
            break;
          }
        }
      else {
        *prob = *prob * qt / qm[ii - j];/*??*/
      }
      if (k>j) vrna_message_error("backtrack failed in qm");

      backtrack_qm1_comparative(vc, pstruc, k, j, prob);

      if (k<i+TURN) break; /* no more pairs */
      r = vrna_urn() * (qm[ii-(k-1)] + expMLbase[k-i]);
      if (expMLbase[k-i] >= r) {
        *prob = *prob * expMLbase[k-i] / (qm[ii-(k-1)] + expMLbase[k-i]);
        break; /* no more pairs */
      }
      j = k-1;
      /* whatishere?? */
    }
  }
  free(type);
}

PRIVATE void
backtrack_qm1_comparative(vrna_fold_compound_t *vc,
              char *pstruc,
              int i,
              int j,
              double *prob){

  int               n_seq       = vc->n_seq;
  short             **S         = vc->S;
  short             **S5        = vc->S5;     /*S5[s][i] holds next base 5' of i in sequence s*/
  short             **S3        = vc->S3;     /*Sl[s][i] holds next base 3' of i in sequence s*/
  vrna_exp_param_t  *pf_params  = vc->exp_params;
  vrna_mx_pf_t      *matrices   = vc->exp_matrices;
  vrna_md_t         *md         = &(pf_params->model_details);
  int               *my_iindx   = vc->iindx;
  int               *jindx      = vc->jindx;
  vrna_hc_t         *hc         = vc->hc;
  vrna_sc_t         **sc        = vc->scs;
  FLT_OR_DBL        *qb         = matrices->qb;
  FLT_OR_DBL        *qm1        = matrices->qm1;
  FLT_OR_DBL        *expMLbase    = matrices->expMLbase;

  /* i is paired to l, i<l<j; backtrack in qm1 to find l */
  int ii, l, xtype,s;
  FLT_OR_DBL qt, r, tempz;
  r = vrna_urn() * qm1[jindx[j]+i];
  ii = my_iindx[i];
  for (qt=0., l=i+TURN+1; l<=j; l++) {
    if (qb[ii-l]==0) continue;
    tempz=1.;
    for (s=0; s<n_seq; s++) {
      xtype = md->pair[S[s][i]][S[s][l]];
      if (xtype==0) xtype=7;
      tempz *= exp_E_MLstem(xtype, S5[s][i], S3[s][l], pf_params);
    }
    qt +=  qb[ii-l]*tempz*expMLbase[j-l];
    if (qt>=r) {
      *prob = *prob
              * qb[ii-l]
              * tempz
              * expMLbase[j-l]
              / qm1[jindx[j]+i];
      /* probs*=qb[ii-l]*tempz*expMLbase[j-l];*/
      break;
    }
  }
  if (l>j) vrna_message_error("backtrack failed in qm1");

  backtrack_comparative(vc, pstruc, i, l, prob);
}

