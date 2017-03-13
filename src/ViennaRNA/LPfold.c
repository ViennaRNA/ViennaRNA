/*
 * local pair probabilities for RNA secondary structures
 *
 * Stephan Bernhart, Ivo L Hofacker
 * Vienna RNA package
 */
/*
 * todo: compute energy z-score for each window
 *
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>    /* #defines FLT_MAX ... */
#include "ViennaRNA/data_structures.h"
#include "ViennaRNA/utils.h"
#include "ViennaRNA/energy_par.h"
#include "ViennaRNA/fold_vars.h"
#include "ViennaRNA/PS_dot.h"
#include "ViennaRNA/part_func.h"
#include "ViennaRNA/params.h"
#include "ViennaRNA/loop_energies.h"
#include "ViennaRNA/LPfold.h"
#include "ViennaRNA/Lfold.h"


/*
 #################################
 # GLOBAL VARIABLES              #
 #################################
 */

typedef struct {
  int           bpp_print;  /* 1 if pairing probabilities should be written to file-handle, 0 if they are returned as vrna_plist_t */
  int           up_print;   /* 1 if unpaired probabilities should be written to file-handle, 0 if they are returned as array */

  FILE          *fp_pU;
  double        **pU;
  FLT_OR_DBL    bpp_cutoff;
  FILE          *fp_bpp;
  vrna_plist_t  *bpp;
  unsigned int  bpp_max_size;
  unsigned int  bpp_size;
  vrna_plist_t  *stack_prob;
  unsigned int  stack_prob_size;
  unsigned int  stack_prob_max_size;
} default_cb_data;

/*
 #################################
 # PRIVATE VARIABLES             #
 #################################
 */

#ifdef  VRNA_BACKWARD_COMPAT

#ifdef _OPENMP
#include <omp.h>
#endif

/* some backward compatibility stuff */
PRIVATE vrna_fold_compound_t  *backward_compat_compound = NULL;
PRIVATE int                   backward_compat           = 0;

#ifdef _OPENMP

#pragma omp threadprivate(backward_compat_compound, backward_compat)

#endif

#endif
/*
 #################################
 # PRIVATE FUNCTION DECLARATIONS #
 #################################
 */

PRIVATE void  GetPtype(vrna_fold_compound_t *vc,
                       int                  j);


PRIVATE void  FreeOldArrays(vrna_fold_compound_t  *vc,
                            int                   i,
                            unsigned int          options);


PRIVATE void  GetNewArrays(vrna_fold_compound_t *vc,
                           int                  j,
                           unsigned int         options);


PRIVATE void  printpbar(vrna_fold_compound_t  *vc,
                        int                   i);


#if 0
PRIVATE vrna_plist_t *get_deppp(vrna_fold_compound_t  *vc,
                                vrna_plist_t          *pl,
                                int                   start);


#endif

PRIVATE void compute_pU(vrna_fold_compound_t  *vc,
                        int                   k,
                        int                   ulength,
                        double                **pU,
                        double                **pUO,
                        double                **pUH,
                        double                **pUI,
                        double                **pUM,
                        unsigned int          options);


PRIVATE FLT_OR_DBL *
compute_stack_probabilities(vrna_fold_compound_t  *vc,
                            int                   start);


PRIVATE void
print_bpp_callback(FLT_OR_DBL *pr,
                   int        size,
                   int        k,
                   void       *data);


PRIVATE void
store_bpp_callback(FLT_OR_DBL *pr,
                   int        size,
                   int        k,
                   void       *data);


PRIVATE void
store_stack_prob_callback(FLT_OR_DBL  *pr,
                          int         size,
                          int         k,
                          void        *data);


PRIVATE void
print_pU_callback(double        *pU,
                  int           size,
                  int           k,
                  int           ulength,
                  unsigned int  type,
                  void          *data);


PRIVATE void
store_pU_callback(double        *pU,
                  int           size,
                  int           k,
                  int           ulength,
                  unsigned int  type,
                  void          *data);


PRIVATE void
backward_compat_callback(FLT_OR_DBL   *pr,
                         int          pr_size,
                         int          i,
                         int          max,
                         unsigned int type,
                         void         *data);


/*
 #################################
 # BEGIN OF FUNCTION DEFINITIONS #
 #################################
 */
PUBLIC vrna_plist_t *
vrna_pfl_fold(const char  *sequence,
              int         window_size,
              int         max_bp_span,
              float       cutoff)
{
  default_cb_data data;

  data.fp_pU                = NULL;
  data.pU                   = NULL;
  data.bpp_cutoff           = (FLT_OR_DBL)cutoff;
  data.fp_bpp               = NULL;
  data.bpp                  = NULL;
  data.bpp_max_size         = 0;
  data.bpp_size             = 0;
  data.stack_prob           = NULL;
  data.stack_prob_max_size  = 0;
  data.stack_prob_size      = 0;
  data.bpp_print            = 0;
  data.up_print             = 0;

  vrna_pfl_fold_cb(sequence, window_size, max_bp_span, &backward_compat_callback, (void *)&data);

  /* resize pair probability list to actual size */
  data.bpp =
    (vrna_plist_t *)vrna_realloc(data.bpp, sizeof(vrna_plist_t) * (data.bpp_size + 1));
  data.bpp[data.bpp_size].i = 0;
  data.bpp[data.bpp_size].j = 0;
  data.bpp[data.bpp_size].p = 0;

  return data.bpp;
}


PUBLIC void
vrna_pfl_fold_cb(const char                 *sequence,
                 int                        window_size,
                 int                        max_bp_span,
                 vrna_probs_window_callback *cb,
                 void                       *data)
{
  unsigned int          options;
  vrna_fold_compound_t  *vc;
  vrna_md_t             md;

  vrna_md_set_default(&md);       /* get default parameters */

  md.compute_bpp  = 1;            /* turn on base pair probability computations */
  md.window_size  = window_size;  /* set size of sliding window */
  md.max_bp_span  = max_bp_span;  /* set maximum base pair span */

  vc = vrna_fold_compound(sequence, &md, VRNA_OPTION_PF | VRNA_OPTION_WINDOW);

  options = VRNA_PROBS_WINDOW_BPP; /* always compute base pair probabilities */

  vrna_probs_window(vc, 0, cb, data, options);

  vrna_fold_compound_free(vc);
}


PUBLIC double **
vrna_pfl_fold_up(const char *sequence,
                 int        ulength,
                 int        window_size,
                 int        max_bp_span)
{
  unsigned int    i;
  double          **pU;
  default_cb_data data;

  pU = NULL;

  if (sequence) {
    i   = strlen(sequence);
    pU  = (double **)vrna_alloc(sizeof(double *) * (i + 2));

    data.fp_pU                = NULL;
    data.pU                   = pU;
    data.bpp_cutoff           = 0.;
    data.fp_bpp               = NULL;
    data.bpp                  = NULL;
    data.bpp_max_size         = 0;
    data.bpp_size             = 0;
    data.stack_prob           = NULL;
    data.stack_prob_max_size  = 0;
    data.stack_prob_size      = 0;
    data.bpp_print            = 0;
    data.up_print             = 0;

    vrna_pfl_fold_up_cb(sequence, ulength, window_size, max_bp_span, &backward_compat_callback,
                        (void *)&data);
  }

  return pU;
}


PUBLIC void
vrna_pfl_fold_up_cb(const char                  *sequence,
                    int                         ulength,
                    int                         window_size,
                    int                         max_bp_span,
                    vrna_probs_window_callback  *cb,
                    void                        *data)
{
  unsigned int          options;
  vrna_fold_compound_t  *vc;
  vrna_md_t             md;

  vrna_md_set_default(&md);       /* get default parameters */

  md.compute_bpp  = 1;            /* turn on base pair probability computations */
  md.window_size  = window_size;  /* set size of sliding window */
  md.max_bp_span  = max_bp_span;  /* set maximum base pair span */

  vc = vrna_fold_compound(sequence, &md, VRNA_OPTION_PF | VRNA_OPTION_WINDOW);

  options = VRNA_PROBS_WINDOW_UP; /* compute unpaired probabilties */

  vrna_probs_window(vc, ulength, cb, data, options);

  vrna_fold_compound_free(vc);
}


PUBLIC void
vrna_probs_window(vrna_fold_compound_t        *vc,
                  int                         ulength,
                  vrna_probs_window_callback  *cb,
                  void                        *data,
                  unsigned int                options)
{
  int               n, m, i, j, k, l, u, u1, type, type_2, tt, ov, noGUclosure;
  double            max_real;
  FLT_OR_DBL        temp, Qmax, prm_MLb, prmt, prmt1, qbt1, *tmp, expMLclosing;
  FLT_OR_DBL        *qqm  = NULL, *qqm1 = NULL, *qq = NULL, *qq1 = NULL;
  FLT_OR_DBL        *prml = NULL, *prm_l = NULL, *prm_l1 = NULL;
  FLT_OR_DBL        *expMLbase, *scale;


  char              *sequence;
  vrna_exp_param_t  *pf_params;
  vrna_md_t         *md;
  short             *S, *S1;
  vrna_mx_pf_t      *matrices;
  int               winSize, pairSize, *rtype, turn;
  FLT_OR_DBL        **q, **qb, **qm, **qm2, **pR, **QI5, **qmb, **q2l;
  double            **pU;
  double            **pUO;
  double            **pUI;
  double            **pUM;
  double            **pUH;
  char              **ptype;

  ov    = 0;
  Qmax  = 0;

  pU  = NULL;
  pUO = NULL;
  pUI = NULL;
  pUM = NULL;
  pUH = NULL;

  if ((!vc) || (!cb))
    return;

  /* here space for initializing everything */

  n             = vc->length;
  sequence      = vc->sequence;
  pf_params     = vc->exp_params;
  md            = &(pf_params->model_details);
  S1            = vc->sequence_encoding;
  S             = vc->sequence_encoding2;
  matrices      = vc->exp_matrices;
  expMLbase     = matrices->expMLbase;
  scale         = matrices->scale;
  winSize       = vc->window_size;
  pairSize      = md->max_bp_span;
  turn          = md->min_loop_size;
  ptype         = vc->ptype_local;
  rtype         = &(md->rtype[0]);
  expMLclosing  = pf_params->expMLclosing;
  noGUclosure   = md->noGUclosure;


  q   = matrices->q_local;
  qb  = matrices->qb_local;
  qm  = matrices->qm_local;
  qm2 = matrices->qm2_local;
  pR  = matrices->pR;
  QI5 = matrices->QI5;
  qmb = matrices->qmb;
  q2l = matrices->q2l;

  /*
   * here, I allocate memory for pU, if has to be saved, I allocate all in one go,
   * if pU is put out and freed, I only allocate what I really need
   */

  /* allocate memory and initialize unpaired probabilities */
  if ((options & VRNA_PROBS_WINDOW_UP) && (ulength > 0)) {
    pU = (double **)vrna_alloc((n + 1) * sizeof(double *));
    if (options & VRNA_PROBS_WINDOW_UP_SPLIT) {
      pUO = (double **)vrna_alloc((n + 1) * sizeof(double *));
      pUI = (double **)vrna_alloc((n + 1) * sizeof(double *));
      pUM = (double **)vrna_alloc((n + 1) * sizeof(double *));
      pUH = (double **)vrna_alloc((n + 1) * sizeof(double *));
      for (i = 1; i <= n; i++) {
        pUH[i]  = (double *)vrna_alloc((MAX2(MAXLOOP, ulength) + 2) * sizeof(double));
        pUI[i]  = (double *)vrna_alloc((MAX2(MAXLOOP, ulength) + 2) * sizeof(double));
        pUO[i]  = (double *)vrna_alloc((MAX2(MAXLOOP, ulength) + 2) * sizeof(double));
        pUM[i]  = (double *)vrna_alloc((MAX2(MAXLOOP, ulength) + 2) * sizeof(double));
      }
    }

    for (i = 1; i <= n; i++)
      pU[i] = (double *)vrna_alloc((MAX2(MAXLOOP, ulength) + 2) * sizeof(double));

    /* very short molecule ? */
    if (n < turn + 2) {
      for (i = 1; i <= n; i++) {
        int maxl = MAX2(MAXLOOP, ulength);
        for (j = 0; j <= maxl; j++)
          pU[i][j] = 1.;

        cb(pU[i], maxl, i, ulength, VRNA_PROBS_WINDOW_UP | VRNA_ANY_LOOP, data);

        if (options & VRNA_PROBS_WINDOW_UP_SPLIT) {
          cb(pU[i], maxl, i, ulength, VRNA_PROBS_WINDOW_UP | VRNA_EXT_LOOP, data);
          /* reset all other loop types to 0% probability */
          for (j = 0; j <= maxl; j++)
            pU[i][j] = 0.;

          cb(pU[i], maxl, i, ulength, VRNA_PROBS_WINDOW_UP | VRNA_HP_LOOP, data);
          cb(pU[i], maxl, i, ulength, VRNA_PROBS_WINDOW_UP | VRNA_INT_LOOP, data);
          cb(pU[i], maxl, i, ulength, VRNA_PROBS_WINDOW_UP | VRNA_MB_LOOP, data);
          free(pUH[i]);
          free(pUI[i]);
          free(pUO[i]);
          free(pUM[i]);
        }

        free(pU[i]);
      }

      free(pU);
      free(pUO);
      free(pUI);
      free(pUM);
      free(pUH);

      return;
    }
  } else if (n < turn + 2) {
    /* very small molecules */
    return;
  }

  qq      = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * (n + 2));
  qq1     = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * (n + 2));
  qqm     = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * (n + 2));
  qqm1    = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * (n + 2));
  prm_l   = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * (n + 2));
  prm_l1  = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * (n + 2));
  prml    = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * (n + 2));


  max_real = (sizeof(FLT_OR_DBL) == sizeof(float)) ? FLT_MAX : DBL_MAX;

  /*  make_ptypes(S, structure); das machmadochlieber lokal, ey! */

  /*
   * array initialization ; qb,qm,q
   * qb,qm,q (i,j) are stored as ((n+1-i)*(n-i) div 2 + n+1-j
   */

  /* ALWAYS q[i][j] => i>j!! */
  for (j = 1; j < MIN2(turn + 2, n); j++) {
    /* allocate start */
    GetNewArrays(vc, j, options);
    GetPtype(vc, j);
    for (i = 1; i <= j; i++)
      q[i][j] = scale[(j - i + 1)];
  }
  for (j = turn + 2; j <= n + winSize; j++) {
    if (j <= n) {
      GetNewArrays(vc, j, options);
      GetPtype(vc, j);
      for (i = MAX2(1, j - winSize); i <= j /* -turn */; i++)
        q[i][j] = scale[(j - i + 1)];
      for (i = j - turn - 1; i >= MAX2(1, (j - winSize + 1)); i--) {
        /* construction of partition function of segment i,j */
        /* firstly that given i bound to j : qb(i,j) */
        u     = j - i - 1;
        type  = ptype[i][j];
        if (type != 0) {
          /* hairpin contribution */
          if (((type == 3) || (type == 4)) && noGUclosure)
            qbt1 = 0;
          else
            qbt1 =
              exp_E_Hairpin(u, type, S1[i + 1], S1[j - 1], sequence + i - 1,
                            pf_params) * scale[u + 2];

          /* interior loops with interior pair k,l */
          for (k = i + 1; k <= MIN2(i + MAXLOOP + 1, j - turn - 2); k++) {
            u1 = k - i - 1;
            for (l = MAX2(k + turn + 1, j - 1 - MAXLOOP + u1); l < j; l++) {
              type_2 = ptype[k][l];
              if (type_2) {
                type_2  = rtype[type_2];
                qbt1    += qb[k][l] *
                           exp_E_IntLoop(u1, j - l - 1, type, type_2,
                                         S1[i + 1], S1[j - 1], S1[k - 1], S1[l + 1],
                                         pf_params) * scale[k - i + j - l];
              }
            }
          }
          /* multiple stem loop contribution */
          temp = 0.0;
          for (k = i + 2; k <= j - 1; k++)
            temp += qm[i + 1][k - 1] * qqm1[k];
          tt    = rtype[type];
          qbt1  += temp * expMLclosing *
                   exp_E_MLstem(tt, S1[j - 1], S1[i + 1], pf_params) * scale[2];

          qb[i][j] = qbt1;
        } /* end if (type!=0) */
        else {
          qb[i][j] = 0.0;
        }

        /*
         * construction of qqm matrix containing final stem
         * contributions to multiple loop partition function
         * from segment i,j
         */
        qqm[i] = qqm1[i] * expMLbase[1];
        if (type) {
          qbt1 = qb[i][j] * exp_E_MLstem(type,
                                         (i > 1) ? S1[i - 1] : -1,
                                         (j < n) ? S1[j + 1] : -1,
                                         pf_params);
          qqm[i] += qbt1;
        }

        /*
         * construction of qm matrix containing multiple loop
         * partition function contributions from segment i,j
         */
        temp = 0.0;
        /* ii = my_iindx[i];   ii-k=[i,k-1] */
        /* new qm2 computation done here */
        for (k = i + 1; k <= j; k++)
          temp += (qm[i][k - 1]) * qqm[k];
        if (ulength > 0)
          qm2[i][j] = temp;           /* new qm2 computation done here */

        for (k = i + 1; k <= j; k++)
          temp += expMLbase[k - i] * qqm[k];
        qm[i][j] = (temp + qqm[i]);

        /* auxiliary matrix qq for cubic order q calculation below */
        qbt1 = qb[i][j];
        if (type)
          qbt1 *=
            exp_E_ExtLoop(type, (i > 1) ? S1[i - 1] : -1, (j < n) ? S1[j + 1] : -1, pf_params);

        qq[i] = qq1[i] * scale[1] + qbt1;

        /* construction of partition function for segment i,j */
        temp = 1.0 * scale[1 + j - i] + qq[i];
        for (k = i; k <= j - 1; k++)
          temp += q[i][k] * qq[k + 1];
        q[i][j] = temp;

        if (temp > Qmax) {
          Qmax = temp;
          if (Qmax > max_real / 10.)
            vrna_message_warning("Q close to overflow: %d %d %g\n", i, j, temp);
        }

        if (temp >= max_real)
          vrna_message_error("overflow in pf_fold while calculating q[%d,%d]\n"
                             "use larger pf_scale", i, j);
      } /* end for i */
      tmp   = qq1;
      qq1   = qq;
      qq    = tmp;
      tmp   = qqm1;
      qqm1  = qqm;
      qqm   = tmp;
    }

    /*
     * just as a general service, I save here the free energy of the windows
     * no output is generated, however,...
     */
    if ((j >= winSize) && (j <= n) && (options & VRNA_PROBS_WINDOW_UP)) {
      FLT_OR_DBL Fwindow = 0.;
      Fwindow = (FLT_OR_DBL)(-log(q[j - winSize + 1][j]) - winSize * log(pf_params->pf_scale)) *
                pf_params->kT / 1000.0;
      /* we could return this to the user via callback cb() if we were nice */

      pU[j][0] = Fwindow;
    }

    if (j > winSize) {
      Qmax = 0;
      /* i=j-winSize; */
      /* initialize multiloopfs */
      for (k = j - winSize; k <= MIN2(n, j); k++) {
        prml[k]   = 0;
        prm_l[k]  = 0;
        /*        prm_l1[k]=0;  others stay */
      }
      prm_l1[j - winSize] = 0;
      k                   = j - winSize;
      for (l = k + turn + 1; l <= MIN2(n, k + winSize - 1); l++) {
        int a;
        pR[k][l]  = 0; /* set zero at start */
        type      = ptype[k][l];
        if (qb[k][l] == 0)
          continue;

        for (a = MAX2(1, l - winSize + 2); a < MIN2(k, n - winSize + 2); a++)
          pR[k][l] += q[a][k - 1] * q[l + 1][a + winSize - 1] / q[a][a + winSize - 1];

        if (l - k + 1 == winSize) {
          pR[k][l] += 1. / q[k][l];
        } else {
          if (k + winSize - 1 <= n)    /* k outermost */
            pR[k][l] += q[l + 1][k + winSize - 1] / q[k][k + winSize - 1];

          if (l - winSize + 1 >= 1) /* l outermost */
            pR[k][l] += q[l - winSize + 1][k - 1] / q[l - winSize + 1][l];
        }

        pR[k][l] *= exp_E_ExtLoop(type,
                                  (k > 1) ? S1[k - 1] : -1,
                                  (l < n) ? S1[l + 1] : -1,
                                  pf_params);

        type_2  = ptype[k][l];
        type_2  = rtype[type_2];

        for (i = MAX2(MAX2(l - winSize + 1, k - MAXLOOP - 1), 1); i <= k - 1; i++) {
          for (m = l + 1; m <= MIN2(MIN2(l + MAXLOOP - k + i + 2, i + winSize - 1), n); m++) {
            type = ptype[i][m];
            if ((pR[i][m] > 0)) {
              pR[k][l] += pR[i][m] * exp_E_IntLoop(k - i - 1,
                                                   m - l - 1,
                                                   type,
                                                   type_2,
                                                   S1[i + 1],
                                                   S1[m - 1],
                                                   S1[k - 1],
                                                   S1[l + 1],
                                                   pf_params) * scale[k - i + m - l];
            }
          }
        }

        if (options & VRNA_PROBS_WINDOW_UP) {
          /* NOT IF WITHIN INNER LOOP */
          for (i = MAX2(MAX2(l - winSize + 1, k - MAXLOOP - 1), 1); i <= k - 1; i++) {
            for (m = l + 1; m <= MIN2(MIN2(l + MAXLOOP - k + i + 2, i + winSize - 1), n); m++) {
              type = ptype[i][m];
              if ((pR[i][m] > 0)) {
                temp = pR[i][m] * qb[k][l] * exp_E_IntLoop(k - i - 1,
                                                           m - l - 1,
                                                           type,
                                                           type_2,
                                                           S1[i + 1],
                                                           S1[m - 1],
                                                           S1[k - 1],
                                                           S1[l + 1],
                                                           pf_params) * scale[k - i + m - l];
                QI5[l][m - l - 1] += temp;
                QI5[i][k - i - 1] += temp;
              }
            }
          }
        }
      }
      /* 3. bonding k,l as substem of multi-loop enclosed by i,m */
      prm_MLb = 0.;
      if (k > 1) {
        /* sonst nix! */
        for (l = MIN2(n - 1, k + winSize - 2); l >= k + turn + 1; l--) {
          /* opposite direction */
          m     = l + 1;
          prmt  = prmt1 = 0.0;
          tt    = ptype[k - 1][m];
          tt    = rtype[tt];
          prmt1 = pR[k - 1][m] *expMLclosing *exp_E_MLstem(tt,
                                                           S1[l],
                                                           S1[k],
                                                           pf_params);


          for (i = MAX2(1, l - winSize + 2); i < k - 1 /* turn */; i++) {
            tt    = ptype[i][m];
            tt    = rtype[tt];
            prmt  += pR[i][m] *
                     exp_E_MLstem(tt, S1[m - 1], S1[i + 1], pf_params) * qm[i + 1][k - 1];
          }
          tt        = ptype[k][l];
          prmt      *= expMLclosing;
          prml[m]   = prmt;
          prm_l[m]  = prm_l1[m] * expMLbase[1] + prmt1;

          prm_MLb = prm_MLb * expMLbase[1] + prml[m];
          /*
           * same as:    prm_MLb = 0;
           * for (i=n; i>k; i--)  prm_MLb += prml[i]*expMLbase[k-i-1];
           */
          prml[m] = prml[m] + prm_l[m];

          if (qb[k][l] == 0.)
            continue;

          temp = prm_MLb;

          if (options & VRNA_PROBS_WINDOW_UP) {
            double dang;
            /* coefficient for computations of unpairedarrays */
            dang = qb[k][l] * exp_E_MLstem(tt, S1[k - 1], S1[l + 1], pf_params) * scale[2];
            for (m = MIN2(k + winSize - 2, n); m >= l + 2; m--) {
              qmb[l][m - l - 1] += prml[m] * dang;
              q2l[l][m - l - 1] += (prml[m] - prm_l[m]) * dang;
            }
          }

          for (m = MIN2(k + winSize - 2, n); m >= l + 2; m--)
            temp += prml[m] * qm[l + 1][m - 1];

          temp *= exp_E_MLstem(tt,
                               (k > 1) ? S1[k - 1] : -1,
                               (l < n) ? S1[l + 1] : -1,
                               pf_params) * scale[2];
          pR[k][l] += temp;

          if (pR[k][l] > Qmax) {
            Qmax = pR[k][l];
            if (Qmax > max_real / 10.)
              vrna_message_warning("P close to overflow: %d %d %g %g\n",
                                   i, m, pR[k][l], qb[k][l]);
          }

          if (pR[k][l] >= max_real) {
            ov++;
            pR[k][l] = FLT_MAX;
          }
        } /* end for (l=..) */
      }

      tmp     = prm_l1;
      prm_l1  = prm_l;
      prm_l   = tmp;

      /* end for (l=..)   */
      if ((options & VRNA_PROBS_WINDOW_UP) && (k - MAXLOOP - 1 > 0)) {
        int start     = k - MAXLOOP - 1;
        int pU_length = MIN2(ulength, start);
        compute_pU(vc, start, ulength, pU, pUO, pUH, pUI, pUM, options);
        cb(pU[start], pU_length, start, ulength, VRNA_PROBS_WINDOW_UP | VRNA_ANY_LOOP, data);
        if (options & VRNA_PROBS_WINDOW_UP_SPLIT) {
          cb(pUO[start], pU_length, start, ulength, VRNA_PROBS_WINDOW_UP | VRNA_EXT_LOOP, data);
          cb(pUH[start], pU_length, start, ulength, VRNA_PROBS_WINDOW_UP | VRNA_HP_LOOP, data);
          cb(pUI[start], pU_length, start, ulength, VRNA_PROBS_WINDOW_UP | VRNA_INT_LOOP, data);
          cb(pUM[start], pU_length, start, ulength, VRNA_PROBS_WINDOW_UP | VRNA_MB_LOOP, data);
        }
      }

      if (j - (2 * winSize + MAXLOOP + 1) > 0) {
        int start = j - (2 * winSize + MAXLOOP + 1);
        printpbar(vc, start);
        if (options & VRNA_PROBS_WINDOW_BPP)
          cb(pR[start], MIN2(start + winSize, n), start, winSize, VRNA_PROBS_WINDOW_BPP, data);

        if (options & VRNA_PROBS_WINDOW_STACKP) {
          int start = j - (2 * winSize - MAXLOOP);
          if (start > 1) {
            FLT_OR_DBL *stack_probs = compute_stack_probabilities(vc, start);
            stack_probs -= start + 1;
            cb(stack_probs, MIN2(n - start + turn,
                                 pairSize), start, winSize, VRNA_PROBS_WINDOW_STACKP, data);
            stack_probs += start + 1;
            free(stack_probs);
          }
        }

        FreeOldArrays(vc, start, options);
      }
    }   /* end if (do_backtrack) */
  }/* end for j */

  /* finish output and free */
  if (options & VRNA_PROBS_WINDOW_UP) {
    for (j = MAX2(1, n - MAXLOOP); j <= n; j++) {
      /* if (pUoutput) pU[j]=(double *)vrna_alloc((ulength+2)*sizeof(double)); */
      compute_pU(vc, j, ulength, pU, pUO, pUH, pUI, pUM, options);
      cb(pU[j], ulength, j, ulength, VRNA_PROBS_WINDOW_UP | VRNA_ANY_LOOP, data);
      if (options & VRNA_PROBS_WINDOW_UP_SPLIT) {
        cb(pUO[j], ulength, j, ulength, VRNA_PROBS_WINDOW_UP | VRNA_EXT_LOOP, data);
        cb(pUH[j], ulength, j, ulength, VRNA_PROBS_WINDOW_UP | VRNA_HP_LOOP, data);
        cb(pUI[j], ulength, j, ulength, VRNA_PROBS_WINDOW_UP | VRNA_INT_LOOP, data);
        cb(pUM[j], ulength, j, ulength, VRNA_PROBS_WINDOW_UP | VRNA_MB_LOOP, data);
      }
    }
  }

  for (j = MAX2(n - winSize - MAXLOOP, 1); j <= n; j++) {
    printpbar(vc, j);
    if (options & VRNA_PROBS_WINDOW_BPP)
      cb(pR[j], MIN2(j + winSize, n), j, winSize, VRNA_PROBS_WINDOW_BPP, data);

    if ((options & VRNA_PROBS_WINDOW_STACKP) && j < n) {
      int start = j;
      if (start > 1) {
        FLT_OR_DBL *stack_probs = compute_stack_probabilities(vc, start);
        stack_probs -= start + 1;
        cb(stack_probs, MIN2(n - start + turn,
                             pairSize), start, winSize, VRNA_PROBS_WINDOW_STACKP, data);
        stack_probs += start + 1;
        free(stack_probs);
      }
    }

    FreeOldArrays(vc, j, options);
  }
  /* free_pf_arrays_L(); */

  if (ov > 0)
    vrna_message_warning("%d overflows occurred while backtracking;\n"
                         "you might try a smaller pf_scale than %g\n",
                         ov, pf_params->pf_scale);

  free(qq);
  free(qq1);
  free(qqm);
  free(qqm1);
  free(prm_l);
  free(prm_l1);
  free(prml);

  if ((options & VRNA_PROBS_WINDOW_UP) && (ulength > 0)) {
    for (i = 1; i <= n; i++)
      free(pU[i]);

    if (options & VRNA_PROBS_WINDOW_UP_SPLIT) {
      for (i = 1; i <= n; i++) {
        free(pUH[i]);
        free(pUI[i]);
        free(pUO[i]);
        free(pUM[i]);
      }
    }
  }

  free(pU);
  free(pUH);
  free(pUI);
  free(pUO);
  free(pUM);
}


PRIVATE void
printpbar(vrna_fold_compound_t  *vc,
          int                   i)
{
  int         j;
  int         howoften = 0; /* how many samples do we have for this pair */
  int         pairdist, turn, n, winSize;
  FLT_OR_DBL  **qb, **prb;

  n       = vc->length;
  winSize = vc->window_size;
  turn    = vc->exp_params->model_details.min_loop_size;


  qb  = vc->exp_matrices->qb_local;
  prb = vc->exp_matrices->pR;

  for (j = i + turn; j < MIN2(i + winSize, n + 1); j++) {
    pairdist = (j - i + 1);
    /* 4cases */
    howoften  = MIN2(winSize - pairdist + 1, i);  /* pairdist,start */
    howoften  = MIN2(howoften, n - j + 1);        /* end */
    howoften  = MIN2(howoften, n - winSize + 1);  /* windowsize */
    prb[i][j] *= qb[i][j] / howoften;
  }
  return;
}


PRIVATE void
FreeOldArrays(vrna_fold_compound_t  *vc,
              int                   i,
              unsigned int          options)
{
  FLT_OR_DBL    **pR, **q, **qb, **qm, **qm2, **QI5, **qmb, **q2l;
  char          **ptype;
  vrna_mx_pf_t  *mx;

  mx    = vc->exp_matrices;
  pR    = mx->pR;
  q     = mx->q_local;
  qb    = mx->qb_local;
  qm    = mx->qm_local;
  ptype = vc->ptype_local;

  /* free arrays no longer needed */
  free(pR[i] + i);
  free(q[i] + i);
  free(qb[i] + i);
  free(qm[i] + i);
  if (options & VRNA_PROBS_WINDOW_UP) {
    qm2 = mx->qm2_local;
    QI5 = mx->QI5;
    qmb = mx->qmb;
    q2l = mx->q2l;
    free(qm2[i] + i);
    free(QI5[i]);
    free(qmb[i]);
    free(q2l[i]);
  }

  free(ptype[i] + i);
  return;
}


PRIVATE void
GetNewArrays(vrna_fold_compound_t *vc,
             int                  j,
             unsigned int         options)
{
  FLT_OR_DBL    **pR, **q, **qb, **qm, **qm2, **QI5, **qmb, **q2l;
  char          **ptype;
  int           winSize;
  vrna_mx_pf_t  *mx;

  mx      = vc->exp_matrices;
  pR      = mx->pR;
  q       = mx->q_local;
  qb      = mx->qb_local;
  qm      = mx->qm_local;
  ptype   = vc->ptype_local;
  winSize = vc->window_size;

  /* allocate new part of arrays */
  pR[j] = (FLT_OR_DBL *)vrna_alloc((winSize + 1) * sizeof(FLT_OR_DBL));
  pR[j] -= j;
  q[j]  = (FLT_OR_DBL *)vrna_alloc((winSize + 1) * sizeof(FLT_OR_DBL));
  q[j]  -= j;
  qb[j] = (FLT_OR_DBL *)vrna_alloc((winSize + 1) * sizeof(FLT_OR_DBL));
  qb[j] -= j;
  qm[j] = (FLT_OR_DBL *)vrna_alloc((winSize + 1) * sizeof(FLT_OR_DBL));
  qm[j] -= j;
  if (options & VRNA_PROBS_WINDOW_UP) {
    qm2     = mx->qm2_local;
    QI5     = mx->QI5;
    qmb     = mx->qmb;
    q2l     = mx->q2l;
    qm2[j]  = (FLT_OR_DBL *)vrna_alloc((winSize + 1) * sizeof(FLT_OR_DBL));
    qm2[j]  -= j;
    QI5[j]  = (FLT_OR_DBL *)vrna_alloc((winSize + 1) * sizeof(FLT_OR_DBL));
    qmb[j]  = (FLT_OR_DBL *)vrna_alloc((winSize + 1) * sizeof(FLT_OR_DBL));
    q2l[j]  = (FLT_OR_DBL *)vrna_alloc((winSize + 1) * sizeof(FLT_OR_DBL));
  }

  ptype[j]  = (char *)vrna_alloc((winSize + 1) * sizeof(char));
  ptype[j]  -= j;
  return;
}


PRIVATE void
GetPtype(vrna_fold_compound_t *vc,
         int                  i)
{
  /* make new entries in ptype array */
  int         j;
  int         type;
  char        **ptype;
  vrna_md_t   *md;
  int         pairSize;
  const short *S;
  int         n;

  ptype     = vc->ptype_local;
  md        = &(vc->exp_params->model_details);
  pairSize  = md->max_bp_span;
  S         = vc->sequence_encoding2;
  n         = vc->length;

  for (j = i; j <= MIN2(i + pairSize, n); j++) {
    type        = md->pair[S[i]][S[j]];
    ptype[i][j] = (char)type;
  }
  return;
}


#if 0
PRIVATE vrna_plist_t *
get_deppp(vrna_fold_compound_t  *vc,
          vrna_plist_t          *pl,
          int                   start)
{
  /* compute dependent pair probabilities */
  int               i, j, count = 0;
  double            tmp;
  vrna_plist_t      *temp;
  char              **ptype;
  short             *S1;
  FLT_OR_DBL        **qb, *scale;
  int               *rtype, turn, pairsize, length;

  vrna_exp_param_t  *pf_params;

  S1        = vc->sequence_encoding;
  pf_params = vc->exp_params;
  ptype     = vc->ptype_local;
  qb        = vc->exp_matrices->qb_local;
  scale     = vc->exp_matrices->scale;
  rtype     = &(pf_params->model_details.rtype[0]);
  turn      = pf_params->model_details.min_loop_size;
  pairsize  = pf_params->model_details.max_bp_span;
  length    = vc->length;

  temp = (vrna_plist_t *)vrna_alloc(pairsize * sizeof(vrna_plist_t)); /* holds temporary deppp */
  for (j = start + turn; j < MIN2(start + pairsize, length); j++) {
    if ((qb[start][j] * qb[start - 1][(j + 1)]) > 10e-200) {
      int type    = ptype[start - 1][j + 1];
      int type_2  = rtype[(unsigned char)ptype[start][j]];
      tmp = qb[start][j] / qb[start - 1][(j + 1)] * exp_E_IntLoop(0,
                                                                  0,
                                                                  type,
                                                                  type_2,
                                                                  S1[start],
                                                                  S1[j],
                                                                  S1[start - 1],
                                                                  S1[j + 1],
                                                                  pf_params) * scale[2];
      temp[count].i   = start;
      temp[count].j   = j;
      temp[count++].p = tmp;
    }
  }
  /* write it to list of deppps */
  for (i = 0; pl[i].i != 0; i++);
  pl = (vrna_plist_t *)vrna_realloc(pl, (i + count + 1) * sizeof(vrna_plist_t));
  for (j = 0; j < count; j++) {
    pl[i + j].i = temp[j].i;
    pl[i + j].j = temp[j].j;
    pl[i + j].p = temp[j].p;
  }
  pl[i + count].i = 0;
  pl[i + count].j = 0;
  pl[i + count].p = 0;
  free(temp);
  return pl;
}


#endif

PRIVATE FLT_OR_DBL *
compute_stack_probabilities(vrna_fold_compound_t  *vc,
                            int                   start)
{
  /* compute dependent pair probabilities */
  int               i, j, count = 0;
  double            tmp;
  char              **ptype;
  short             *S1;
  FLT_OR_DBL        **qb, *scale, *probs;
  int               *rtype, turn, pairsize, length;

  vrna_exp_param_t  *pf_params;

  length    = vc->length;
  S1        = vc->sequence_encoding;
  pf_params = vc->exp_params;
  ptype     = vc->ptype_local;
  qb        = vc->exp_matrices->qb_local;
  scale     = vc->exp_matrices->scale;
  rtype     = &(pf_params->model_details.rtype[0]);
  turn      = pf_params->model_details.min_loop_size;
  pairsize  = pf_params->model_details.max_bp_span;

  int max_j = MIN2(start + pairsize, length) - 1;

  probs = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * (max_j - start + 1));

  for (j = start + turn + 1; j <= max_j; j++) {
    if ((qb[start][j] * qb[start - 1][(j + 1)]) > 10e-200) {
      int type    = ptype[start - 1][j + 1];
      int type_2  = rtype[(unsigned char)ptype[start][j]];
      tmp = qb[start][j] / qb[start - 1][(j + 1)] * exp_E_IntLoop(0,
                                                                  0,
                                                                  type,
                                                                  type_2,
                                                                  S1[start],
                                                                  S1[j],
                                                                  S1[start - 1],
                                                                  S1[j + 1],
                                                                  pf_params) * scale[2];
      probs[j - start - 1] = tmp;
    }
  }
  return probs;
}


/*
 * Here: Space for questions...
 */
PRIVATE void
compute_pU(vrna_fold_compound_t *vc,
           int                  k,
           int                  ulength,
           double               **pU,
           double               **pUO,
           double               **pUH,
           double               **pUI,
           double               **pUM,
           unsigned int         options)
{
  /*
   *  here, we try to add a function computing all unpaired probabilities starting at some i,
   *  going down to $unpaired, to be unpaired, i.e. a list with entries from 1 to unpaired for
   *  every i, with the probability of a stretch of length x, starting at i-x+1, to be unpaired
   */
  int               startu;
  int               i5;
  int               j3, len, obp, *rtype, turn;
  double            temp;
  double            *QBE;
  double            *QBI;
  double            *QBM;
  double            *QBH;

  int               winSize;
  int               n;
  char              *sequence;
  FLT_OR_DBL        expMLclosing, *expMLbase, **q, **qm, **qm2, *scale, **pR, **QI5, **q2l, **qmb;
  vrna_exp_param_t  *pf_params;
  char              **ptype;
  short             *S1;

  sequence      = vc->sequence;
  n             = vc->length;
  winSize       = vc->window_size;
  S1            = vc->sequence_encoding;
  pf_params     = vc->exp_params;
  ptype         = vc->ptype_local;
  rtype         = &(pf_params->model_details.rtype[0]);
  scale         = vc->exp_matrices->scale;
  q             = vc->exp_matrices->q_local;
  qm            = vc->exp_matrices->qm_local;
  qm2           = vc->exp_matrices->qm2_local;
  expMLbase     = vc->exp_matrices->expMLbase;
  expMLclosing  = pf_params->expMLclosing;
  pR            = vc->exp_matrices->pR;
  QI5           = vc->exp_matrices->QI5;
  q2l           = vc->exp_matrices->q2l;
  qmb           = vc->exp_matrices->qmb;
  turn          = pf_params->model_details.min_loop_size;


  QBE = (double *)vrna_alloc((MAX2(ulength, MAXLOOP) + 2) * sizeof(double));
  QBM = (double *)vrna_alloc((MAX2(ulength, MAXLOOP) + 2) * sizeof(double));
  QBI = (double *)vrna_alloc((MAX2(ulength, MAXLOOP) + 2) * sizeof(double));
  QBH = (double *)vrna_alloc((MAX2(ulength, MAXLOOP) + 2) * sizeof(double));

  /* first, we will */
  /* for k<=ulength, pU[k][k]=0, because no bp can enclose it */

  /* compute pu[k+ulength][ulength] */
  for (i5 = MAX2(k + ulength - winSize + 1, 1); i5 <= k; i5++) {
    for (j3 = k + ulength + 1; j3 <= MIN2(n, i5 + winSize - 1); j3++) {
      /*
       *  if (k>400) {
       * printf("i%d j%d  ",i5,j3);
       * fflush(stdout);
       * }
       */
      if (ptype[i5][j3] != 0) {
        /*
         * (.. >-----|..........)
         * i5  j     j+ulength  j3
         */
        /* Multiloops */
        /* (..{}{}-----|......) */
        temp = (i5 < k) ? qm2[i5 + 1][k] * expMLbase[j3 - k - 1] : 0.;

        /* (..|-----|{}{}) */
        if (j3 - 1 > k + ulength)
          temp += qm2[k + ulength + 1][j3 - 1] * expMLbase[k + ulength - i5];

        /* ({}|-----|{}) */
        if ((i5 < k) && (j3 - 1 > k + ulength))
          temp += qm[i5 + 1][k] * qm[k + ulength + 1][j3 - 1] * expMLbase[ulength];

        /* add dangles, multloopclosing etc. */
        temp *=
          exp_E_MLstem(rtype[(unsigned char)ptype[i5][j3]], S1[j3 - 1], S1[i5 + 1],
                       pf_params) * scale[2] * expMLclosing;
        /* add hairpins */
        temp += exp_E_Hairpin(j3 - i5 - 1,
                              ptype[i5][j3],
                              S1[i5 + 1],
                              S1[j3 - 1],
                              sequence + i5 - 1,
                              pf_params) * scale[j3 - i5 + 1];
        /* add outer probability */
        temp                      *= pR[i5][j3];
        pU[k + ulength][ulength]  += temp;
      }
    }
  }
  /* code doubling to avoid if within loop */
  temp = 0.;
  for (len = winSize; len >= MAX2(ulength, MAXLOOP); len--)
    temp += QI5[k][len];
  for (; len > 0; len--) {
    temp      += QI5[k][len];
    QBI[len]  += temp;
    QBE[len]  += temp; /* replace QBE with QI */
  }
  /* Add Hairpinenergy to QBE */
  temp = 0.;
  for (obp = MIN2(n, k + winSize - 1); obp > k + ulength; obp--)
    if (ptype[k][obp]) {
      temp += pR[k][obp] * exp_E_Hairpin(obp - k - 1,
                                         ptype[k][obp],
                                         S1[k + 1],
                                         S1[obp - 1],
                                         sequence + k - 1,
                                         pf_params) * scale[obp - k + 1];
    }

  for (obp = MIN2(n, MIN2(k + winSize - 1, k + ulength)); obp > k + 1; obp--) {
    if (ptype[k][obp]) {
      temp += pR[k][obp] * exp_E_Hairpin(obp - k - 1,
                                         ptype[k][obp],
                                         S1[k + 1],
                                         S1[obp - 1],
                                         sequence + k - 1,
                                         pf_params) * scale[obp - k + 1];
    }

    QBH[obp - k - 1]  += temp;
    QBE[obp - k - 1]  += temp; /* add hairpins to QBE (all in one array) */
  }
  /* doubling the code to get the if out of the loop */

  /*
   * Add up Multiloopterms  qmb[l][m]+=prml[m]*dang;
   * q2l[l][m]+=(prml[m]-prm_l[m])*dang;
   */

  temp = 0.;
  for (len = winSize; len >= ulength; len--)
    temp += q2l[k][len] * expMLbase[len];
  for (; len > 0; len--) {
    temp      += q2l[k][len] * expMLbase[len];
    QBM[len]  += temp;
    QBE[len]  += temp; /* add (()()____) type cont. to I3 */
  }
  for (len = 1; len < ulength; len++) {
    for (obp = k + len + turn; obp <= MIN2(n, k + winSize - 1); obp++) {
      /* add (()___()) */
      QBM[len]  += qmb[k][obp - k - 1] * qm[k + len + 1 /*2*/][obp - 1] * expMLbase[len];
      QBE[len]  += qmb[k][obp - k - 1] * qm[k + len + 1 /*2*/][obp - 1] * expMLbase[len];
    }
  }
  for (len = 1; len < ulength; len++) {
    for (obp = k + len + turn + turn; obp <= MIN2(n, k + winSize - 1); obp++) {
      if (ptype[k][obp]) {
        temp = exp_E_MLstem(rtype[(unsigned char)ptype[k][obp]],
                            S1[obp - 1],
                            S1[k + 1],
                            pf_params) * scale[2] * expMLbase[len] * expMLclosing;      /* k:obp */
        QBE[len]  += pR[k][obp] * temp * qm2[k + len + 1][obp - 1];                     /* add (___()()) */
        QBM[len]  += pR[k][obp] * temp * qm2[k + len + 1][obp - 1];                     /* add (___()()) */
      }
    }
  }
  /*
   * After computing all these contributions in QBE[len], that k is paired
   * and the unpaired stretch is AT LEAST len long, we start to add that to
   * the old unpaired thingies;
   */
  for (len = 1; len < MIN2(MAX2(ulength, MAXLOOP), n - k); len++)
    pU[k + len][len] += pU[k + len][len + 1] + QBE[len];

  if (options & VRNA_PROBS_WINDOW_UP_SPLIT) {
    for (len = 1; len < MIN2(MAX2(ulength, MAXLOOP), n - k); len++) {
      pUH[k + len][len] += pUH[k + len][len + 1] + QBH[len];
      pUM[k + len][len] += pUM[k + len][len + 1] + QBM[len];
      pUI[k + len][len] += pUI[k + len][len + 1] + QBI[len];
    }
    /* open chain */
    if ((ulength >= winSize) && (k >= ulength))
      pUO[k][winSize] = scale[winSize] / q[k - winSize + 1][k];
  }

  /* open chain */
  if ((ulength >= winSize) && (k >= ulength))
    pU[k][winSize] = scale[winSize] / q[k - winSize + 1][k];

  /*
   * now the not enclosed by any base pair terms for whatever it is we do not need anymore...
   * ... which should be e.g; k, again
   */
  int leftmost, rightmost;
  for (startu = MIN2(ulength, k); startu > 0; startu--) {
    temp = 0.;
    for (i5 = MAX2(1, k - winSize + 2); i5 <= MIN2(k - startu, n - winSize + 1); i5++)
      temp += q[i5][k - startu] * q[k + 1][i5 + winSize - 1] * scale[startu] /
              q[i5][i5 + winSize - 1];
    /* the 2 Cases where the borders are on the edge of the interval */
    if ((k >= winSize) && (startu + 1 <= winSize))
      temp += q[k - winSize + 1][k - startu] * scale[startu] / q[k - winSize + 1][k];

    if ((k <= n - winSize + startu) && (k - startu >= 0) && (k < n) && (startu + 1 <= winSize))
      temp += q[k + 1][k - startu + winSize] * scale[startu] /
              q[k - startu + 1][k - startu + winSize];

    /* Divide by number of possible windows */
    leftmost  = MAX2(1, k - winSize + 1);
    rightmost = MIN2(n - winSize + 1, k - startu + 1);

    pU[k][startu] += temp;
    pU[k][startu] /= (rightmost - leftmost + 1);

    if (options & VRNA_PROBS_WINDOW_UP_SPLIT) {
      pUO[k][startu] += temp;

      /* Do we want to make a distinction between those? */
      pUO[k][startu]  /= (rightmost - leftmost + 1);
      pUH[k][startu]  /= (rightmost - leftmost + 1);
      pUI[k][startu]  /= (rightmost - leftmost + 1);
      pUM[k][startu]  /= (rightmost - leftmost + 1);
    }
  }
  free(QBE);
  free(QBI);
  free(QBH);
  free(QBM);
  return;
}


PRIVATE void
print_bpp_callback(FLT_OR_DBL *pr,
                   int        size,
                   int        k,
                   void       *data)
{
  int         j;
  FILE        *fp     = ((default_cb_data *)data)->fp_bpp;
  FLT_OR_DBL  cutoff  = ((default_cb_data *)data)->bpp_cutoff;

  for (j = k + 1; j <= size; j++) {
    if (pr[j] < cutoff)
      continue;

    fprintf(fp, "%d  %d  %g\n", k, j, pr[j]);
  }
}


PRIVATE void
store_bpp_callback(FLT_OR_DBL *pr,
                   int        size,
                   int        k,
                   void       *data)
{
  int           j;
  vrna_plist_t  *pl         = ((default_cb_data *)data)->bpp;
  unsigned int  pl_size     = ((default_cb_data *)data)->bpp_size;
  unsigned int  pl_max_size = ((default_cb_data *)data)->bpp_max_size;
  FLT_OR_DBL    cutoff      = ((default_cb_data *)data)->bpp_cutoff;

  if (pl_max_size == 0) {
    /* init if necessary */
    pl_max_size = 100;
    pl          = (vrna_plist_t *)vrna_realloc(pl, sizeof(vrna_plist_t) * pl_max_size);
  }

  for (j = k + 1; j <= size; j++) {
    if (pr[j] < cutoff)
      continue;

    /* resize vrna_plist_t memory if necessary */
    if (pl_size >= pl_max_size - 1) {
      pl_max_size *= 1.5;
      pl          = (vrna_plist_t *)vrna_realloc(pl, sizeof(vrna_plist_t) * pl_max_size);
    }

    pl[pl_size].i   = k;
    pl[pl_size].j   = j;
    pl[pl_size++].p = pr[j];
  }

  /* mark end of vrna_plist_t */
  pl[pl_size].i = 0;
  pl[pl_size].j = 0;
  pl[pl_size].p = 0.;

  /* update data */
  ((default_cb_data *)data)->bpp          = pl;
  ((default_cb_data *)data)->bpp_size     = pl_size;
  ((default_cb_data *)data)->bpp_max_size = pl_max_size;
}


PRIVATE void
store_stack_prob_callback(FLT_OR_DBL  *pr,
                          int         size,
                          int         k,
                          void        *data)
{
  int           j;
  vrna_plist_t  *pl         = ((default_cb_data *)data)->stack_prob;
  unsigned int  pl_size     = ((default_cb_data *)data)->stack_prob_size;
  unsigned int  pl_max_size = ((default_cb_data *)data)->stack_prob_max_size;
  FLT_OR_DBL    cutoff      = ((default_cb_data *)data)->bpp_cutoff;

  if (pl_max_size == 0) {
    /* init if necessary */
    pl_max_size = 100;
    pl          = (vrna_plist_t *)vrna_realloc(pl, sizeof(vrna_plist_t) * pl_max_size);
  }

  for (j = k + 1; j <= size; j++) {
    if (pr[j] < cutoff)
      continue;

    /* resize vrna_plist_t memory if necessary */
    if (pl_size >= pl_max_size - 1) {
      pl_max_size *= 1.5;
      pl          = (vrna_plist_t *)vrna_realloc(pl, sizeof(vrna_plist_t) * pl_max_size);
    }

    pl[pl_size].i   = k;
    pl[pl_size].j   = j;
    pl[pl_size++].p = pr[j];
  }

  /* mark end of vrna_plist_t */
  pl[pl_size].i = 0;
  pl[pl_size].j = 0;
  pl[pl_size].p = 0.;

  /* update data */
  ((default_cb_data *)data)->stack_prob           = pl;
  ((default_cb_data *)data)->stack_prob_size      = pl_size;
  ((default_cb_data *)data)->stack_prob_max_size  = pl_max_size;
}


PRIVATE void
print_pU_callback(double        *pU,
                  int           size,
                  int           k,
                  int           ulength,
                  unsigned int  type,
                  void          *data)
{
  if (type & VRNA_PROBS_WINDOW_UP) {
    int   i;
    FILE  *fp = ((default_cb_data *)data)->fp_pU;

    fprintf(fp, "%d\t", k);

    for (i = 1; i < size; i++)
      fprintf(fp, "%.7g\t", pU[i]);
    fprintf(fp, "%.7g", pU[size]);

    if ((type & VRNA_ANY_LOOP) == VRNA_ANY_LOOP)
      fprintf(fp, "\n");
    else if (type & VRNA_EXT_LOOP)
      fprintf(fp, "\tE\n");
    else if (type & VRNA_HP_LOOP)
      fprintf(fp, "\tH\n");
    else if (type & VRNA_INT_LOOP)
      fprintf(fp, "\tI\n");
    else if (type & VRNA_MB_LOOP)
      fprintf(fp, "\tM\n");
    else
      vrna_message_warning("unknown loop type");
  }
}


PRIVATE void
store_pU_callback(double        *pU,
                  int           size,
                  int           k,
                  int           ulength,
                  unsigned int  type,
                  void          *data)
{
  int     i;
  double  **pU_storage = ((default_cb_data *)data)->pU;

  if ((type & VRNA_PROBS_WINDOW_UP) && ((type & VRNA_ANY_LOOP) == VRNA_ANY_LOOP)) {
    pU_storage[k] = (double *)vrna_alloc(sizeof(double) * (ulength + 1));
    for (i = 1; i <= size; i++)
      pU_storage[k][i] = pU[i];
  }
}


PRIVATE void
backward_compat_callback(FLT_OR_DBL   *pr,
                         int          pr_size,
                         int          i,
                         int          max,
                         unsigned int type,
                         void         *data)
{
  default_cb_data *d = (default_cb_data *)data;

  if (type & VRNA_PROBS_WINDOW_BPP) {
    if (d->bpp_print)
      print_bpp_callback(pr, pr_size, i, data);
    else
      store_bpp_callback(pr, pr_size, i, data);
  } else if (type & VRNA_PROBS_WINDOW_UP) {
    if (d->up_print)
      print_pU_callback(pr, pr_size, i, max, type, data);
    else
      store_pU_callback(pr, pr_size, i, max, type, data);
  }
}


/*###########################################*/
/*# deprecated functions below              #*/
/*###########################################*/

PRIVATE vrna_plist_t *
wrap_pf_foldLP(char             *sequence,
               int              winSize,
               int              pairSize,
               float            cutoffb,
               double           **pU,
               vrna_plist_t     **dpp2,
               FILE             *pUfp,
               FILE             *spup,
               vrna_exp_param_t *parameters)
{
  int                   ulength;
  vrna_fold_compound_t  *vc;
  vrna_md_t             md;
  default_cb_data       data;

  vc      = NULL;
  ulength = 0;

  /* we need vrna_exp_param_t datastructure to correctly init default hard constraints */
  if (parameters)
    md = parameters->model_details;
  else
    set_model_details(&md);   /* get global default parameters */

  md.compute_bpp  = 1;        /* turn on base pair probability computations */
  md.window_size  = winSize;  /* set size of sliding window */
  md.max_bp_span  = pairSize; /* set maximum base pair span */

  vc = vrna_fold_compound(sequence, &md, VRNA_OPTION_PF | VRNA_OPTION_WINDOW);

  /* prepare exp_params and set global pf_scale */
  free(vc->exp_params);
  vc->exp_params            = vrna_exp_params(&md);
  vc->exp_params->pf_scale  = pf_scale;

  if (backward_compat_compound && backward_compat)
    vrna_fold_compound_free(backward_compat_compound);

  backward_compat_compound  = vc;
  backward_compat           = 1;
  iindx                     = backward_compat_compound->iindx; /* for backward compatibility and Perl wrapper */

  if (pU)
    ulength = (int)pU[0][0] + 0.49;

  data.fp_pU                = pUfp;
  data.pU                   = pU;
  data.bpp_cutoff           = (FLT_OR_DBL)cutoffb;
  data.fp_bpp               = spup;
  data.bpp                  = NULL;
  data.bpp_max_size         = 0;
  data.bpp_size             = 0;
  data.stack_prob           = NULL;
  data.stack_prob_max_size  = 0;
  data.stack_prob_size      = 0;

  data.bpp_print  = (spup) ? 1 : 0;
  data.up_print   = (pUfp) ? 1 : 0;

  unsigned int options = VRNA_PROBS_WINDOW_BPP; /* always compute base pair probabilities */

  if (dpp2 && (*dpp2))
    options |= VRNA_PROBS_WINDOW_STACKP;

  if (ulength > 0)
    options |= VRNA_PROBS_WINDOW_UP;

  vrna_probs_window(vc, ulength, &backward_compat_callback, (void *)&data, options);

  if (dpp2 && (*dpp2)) {
    data.stack_prob = (vrna_plist_t *)vrna_realloc(data.stack_prob,
                                                   sizeof(vrna_plist_t) *
                                                   (data.stack_prob_size + 1));
    data.stack_prob[data.stack_prob_size].i = 0;
    data.stack_prob[data.stack_prob_size].j = 0;
    data.stack_prob[data.stack_prob_size].p = 0;
    free(*dpp2); /* free already occupied memory */
    *dpp2 = data.stack_prob;
  }

  if (!spup) {
    data.bpp =
      (vrna_plist_t *)vrna_realloc(data.bpp, sizeof(vrna_plist_t) * (data.bpp_size + 1));
    data.bpp[data.bpp_size].i = 0;
    data.bpp[data.bpp_size].j = 0;
    data.bpp[data.bpp_size].p = 0;
    return data.bpp;
  } else {
    return NULL;
  }
}


PUBLIC void
init_pf_foldLP(int length)
{
  /* DO NOTHING */
}


PUBLIC void
update_pf_paramsLP(int length)
{
  if (backward_compat_compound && backward_compat) {
    vrna_md_t md;
    set_model_details(&md);
    vrna_exp_params_reset(backward_compat_compound, &md);

    /* compatibility with RNAup, may be removed sometime */
    pf_scale = backward_compat_compound->exp_params->pf_scale;
  }
}


PUBLIC void
update_pf_paramsLP_par(int              length,
                       vrna_exp_param_t *parameters)
{
  if (backward_compat_compound && backward_compat) {
    vrna_md_t md;
    if (parameters) {
      vrna_exp_params_subst(backward_compat_compound, parameters);
    } else {
      set_model_details(&md);
      vrna_exp_params_reset(backward_compat_compound, &md);
    }

    /* compatibility with RNAup, may be removed sometime */
    pf_scale = backward_compat_compound->exp_params->pf_scale;
  }
}


PUBLIC vrna_plist_t *
pfl_fold(char         *sequence,
         int          winSize,
         int          pairSize,
         float        cutoffb,
         double       **pU,
         vrna_plist_t **dpp2,
         FILE         *pUfp,
         FILE         *spup)
{
  return wrap_pf_foldLP(sequence, winSize, pairSize, cutoffb, pU, dpp2, pUfp, spup, NULL);
}


PUBLIC vrna_plist_t *
pfl_fold_par(char             *sequence,
             int              winSize,
             int              pairSize,
             float            cutoffb,
             double           **pU,
             vrna_plist_t     **dpp2,
             FILE             *pUfp,
             FILE             *spup,
             vrna_exp_param_t *parameters)
{
  return wrap_pf_foldLP(sequence, winSize, pairSize, cutoffb, pU, dpp2, pUfp, spup, parameters);
}


PUBLIC void
putoutpU_prob(double  **pU,
              int     length,
              int     ulength,
              FILE    *fp,
              int     energies)
{
  if (backward_compat_compound && backward_compat)
    putoutpU_prob_par(pU, length, ulength, fp, energies, backward_compat_compound->exp_params);
  else
    vrna_message_warning("putoutpU_prob: Not doing anything! First, run pfl_fold()!");
}


PUBLIC void
putoutpU_prob_par(double            **pU,
                  int               length,
                  int               ulength,
                  FILE              *fp,
                  int               energies,
                  vrna_exp_param_t  *parameters)
{
  /* put out unpaireds */
  int     i, k;
  double  kT = parameters->kT / 1000.0;
  double  temp;

  if (energies)
    fprintf(fp, "#opening energies\n #i$\tl=");
  else
    fprintf(fp, "#unpaired probabilities\n #i$\tl=");

  for (i = 1; i <= ulength; i++)
    fprintf(fp, "%d\t", i);
  fprintf(fp, "\n");

  for (k = 1; k <= length; k++) {
    fprintf(fp, "%d\t", k);
    for (i = 1; i <= ulength; i++) {
      if (i > k) {
        fprintf(fp, "NA\t");
        continue;
      }

      if (energies)
        temp = -log(pU[k][i]) * kT;
      else
        temp = pU[k][i];

      fprintf(fp, "%.7g\t", temp);
    }
    fprintf(fp, "\n");
    free(pU[k]);
  }
  fflush(fp);
}


PUBLIC void
putoutpU_prob_bin(double  **pU,
                  int     length,
                  int     ulength,
                  FILE    *fp,
                  int     energies)
{
  if (backward_compat_compound && backward_compat)
    putoutpU_prob_bin_par(pU, length, ulength, fp, energies, backward_compat_compound->exp_params);
  else
    vrna_message_warning("putoutpU_prob_bin: Not doing anything! First, run pfl_fold()!");
}


PUBLIC void
putoutpU_prob_bin_par(double            **pU,
                      int               length,
                      int               ulength,
                      FILE              *fp,
                      int               energies,
                      vrna_exp_param_t  *parameters)
{
  /* put out unpaireds */
  int     i, k;
  double  kT = parameters->kT / 1000.0;
  int     *p;

  p = (int *)vrna_alloc(sizeof(int) * 1);
  /* write first line */
  p[0] = ulength; /* u length */
  fwrite(p, sizeof(int), 1, fp);
  p[0] = length;  /* seq length */
  fwrite(p, sizeof(int), 1, fp);
  for (k = 3; k <= (length + 20); k++) {
    /* all the other lines are set to 1000000 because we are at ulength=0 */
    p[0] = 1000000;
    fwrite(p, sizeof(int), 1, fp);
  }
  /* data */
  for (i = 1; i <= ulength; i++) {
    for (k = 1; k <= 11; k++) {
      /* write first ten entries to 1000000 */
      p[0] = 1000000;
      fwrite(p, sizeof(int), 1, fp);
    }
    for (k = 1; k <= length; k++) {
      /* write data now */
      if (i > k) {
        p[0] = 1000000;         /* check if u > pos */
        fwrite(p, sizeof(int), 1, fp);
        continue;
      } else {
        p[0] = (int)rint(100 * (-log(pU[k][i]) * kT));
        fwrite(p, sizeof(int), 1, fp);
      }
    }
    for (k = 1; k <= 9; k++) {
      /* finish by writing the last 10 entries */
      p[0] = 1000000;
      fwrite(p, sizeof(int), 1, fp);
    }
  }
  /* free pU array; */
  for (k = 1; k <= length; k++)
    free(pU[k]);
  free(p);
  fflush(fp);
}
