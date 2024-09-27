/*
 *
 *                c Ivo Hofacker
 *
 *                Vienna RNA package
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "ViennaRNA/params/default.h"
#include "ViennaRNA/fold_vars.h"
#include "ViennaRNA/utils/basic.h"
#include "ViennaRNA/params/io.h"
#include "ViennaRNA/params/basic.h"
#include "ViennaRNA/params/salt.h"

/**
 *** \file ViennaRNA/params/basic.c
 *** <P>
 *** This file provides functions that return temperature scaled energy parameters and
 *** Boltzmann weights packed in datastructures
 *** </P>
 ***/

/*------------------------------------------------------------------------*/
#define SCALE 10
/**
 *** dangling ends should never be destabilizing, i.e. expdangle>=1<BR>
 *** specific heat needs smooth function (2nd derivative)<BR>
 *** we use a*(sin(x+b)+1)^2, with a=2/(3*sqrt(3)), b=Pi/6-sqrt(3)/2,
 *** in the interval b<x<sqrt(3)/2
 */
#define CLIP_NEGATIVE(X) ((X) < 0 ? 0 : (X))
#define SMOOTH(X) ((!pf_smooth)               ?   CLIP_NEGATIVE(X) : \
                   ((X) / SCALE < -1.2283697) ?                 0  : \
                   ((X) / SCALE > 0.8660254) ?                (X) : \
                   SCALE *0.38490018   \
                   * (sin((X) / SCALE - 0.34242663) + 1) \
                   * (sin((X) / SCALE - 0.34242663) + 1) \
                   )

/* #define SMOOTH(X) ((X)<0 ? 0 : (X)) */

/*
 * If the global use_mfelike_energies flag is set, truncate doubles to int
 * values and cast back to double. This makes the energy parameters of the
 * partition (folding get_scaled_exp_params()) compatible with the mfe folding
 * parameters (get_scaled_exp_params()), e.g. for explicit partition function
 * computations.
 */
#define TRUNC_MAYBE(X) ((!pf_smooth) ? (double)((int)(X)) : (X))


/* Rescale Free energy contribution according to deviation of temperature from measurement conditions */
#define RESCALE_dG(dG, dH, dT)   ((dH) - ((dH) - (dG)) * dT)

/*
 * Rescale Free energy contribution according to deviation of temperature from measurement conditions
 * and convert it to Boltzmann Factor for specific kT
 */
#define RESCALE_BF(dG, dH, dT, kT)          ( \
    exp( \
      -TRUNC_MAYBE((double)RESCALE_dG((dG), (dH), (dT))) \
      * 10. \
      / kT \
      ) \
    )

#define RESCALE_BF_SMOOTH(dG, dH, dT, kT)   ( \
    exp(  \
      SMOOTH( \
        -TRUNC_MAYBE((double)RESCALE_dG((dG), (dH), (dT))) \
        ) \
      * 10. \
      / kT \
      ) \
    )

#define saltT md->temperature+K0

/*
 #################################
 # PRIVATE VARIABLES             #
 #################################
 */
PRIVATE vrna_param_t p;
PRIVATE int               id = -1;
/* variables for partition function */
PRIVATE vrna_exp_param_t  pf;
PRIVATE int               pf_id = -1;

#ifdef _OPENMP
#pragma omp threadprivate(id, pf_id)
#endif

/*
 #################################
 # PRIVATE FUNCTION DECLARATIONS #
 #################################
 */

PRIVATE vrna_param_t *
get_scaled_params(vrna_md_t *md);


PRIVATE vrna_exp_param_t *
get_scaled_exp_params(vrna_md_t *md,
                      double    pfs);


PRIVATE vrna_exp_param_t *
get_exp_params_ali(vrna_md_t    *md,
                   unsigned int n_seq,
                   double       pfs);


PRIVATE void
rescale_params(vrna_fold_compound_t *vc);


/*
 #################################
 # BEGIN OF FUNCTION DEFINITIONS #
 #################################
 */
PUBLIC vrna_param_t *
vrna_params(vrna_md_t *md)
{
  if (md) {
    return get_scaled_params(md);
  } else {
    vrna_md_t md;
    vrna_md_set_default(&md);
    return get_scaled_params(&md);
  }
}


PUBLIC vrna_exp_param_t *
vrna_exp_params(vrna_md_t *md)
{
  if (md) {
    return get_scaled_exp_params(md, -1.);
  } else {
    vrna_md_t md;
    vrna_md_set_default(&md);
    return get_scaled_exp_params(&md, -1.);
  }
}


PUBLIC vrna_exp_param_t *
vrna_exp_params_comparative(unsigned int  n_seq,
                            vrna_md_t     *md)
{
  if (md) {
    return get_exp_params_ali(md, n_seq, -1.);
  } else {
    vrna_md_t md;
    vrna_md_set_default(&md);
    return get_exp_params_ali(&md, n_seq, -1.);
  }
}


PUBLIC vrna_param_t *
vrna_params_copy(vrna_param_t *par)
{
  vrna_param_t *copy = NULL;

  if (par) {
    copy = (vrna_param_t *)vrna_alloc(sizeof(vrna_param_t));
    memcpy(copy, par, sizeof(vrna_param_t));
  }

  return copy;
}


PUBLIC vrna_exp_param_t *
vrna_exp_params_copy(vrna_exp_param_t *par)
{
  vrna_exp_param_t *copy = NULL;

  if (par) {
    copy = (vrna_exp_param_t *)vrna_alloc(sizeof(vrna_exp_param_t));
    memcpy(copy, par, sizeof(vrna_exp_param_t));
  }

  return copy;
}


PUBLIC void
vrna_params_subst(vrna_fold_compound_t  *vc,
                  vrna_param_t          *parameters)
{
  if (vc) {
    if (vc->params)
      free(vc->params);

    if (parameters) {
      vc->params = vrna_params_copy(parameters);
    } else {
      switch (vc->type) {
        case VRNA_FC_TYPE_SINGLE:     /* fall through */

        case VRNA_FC_TYPE_COMPARATIVE:
          vc->params = vrna_params(NULL);
          break;

        default:
          break;
      }
    }
  }
}


PUBLIC void
vrna_params_reset(vrna_fold_compound_t  *vc,
                  vrna_md_t             *md_p)
{
  if (vc) {
    switch (vc->type) {
      case VRNA_FC_TYPE_SINGLE:     /* fall through */

      case VRNA_FC_TYPE_COMPARATIVE:
        if (vc->params)
          free(vc->params);

        vc->params = vrna_params(md_p);

        if (vc->exp_params) {
          free(vc->exp_params);

          vc->exp_params = vrna_exp_params(md_p);
        }

        break;

      default:
        break;
    }
  }
}


PUBLIC void
vrna_exp_params_reset(vrna_fold_compound_t  *vc,
                      vrna_md_t             *md_p)
{
  if (vc) {
    switch (vc->type) {
      case VRNA_FC_TYPE_SINGLE:     /* fall through */

      case VRNA_FC_TYPE_COMPARATIVE:
        if (vc->exp_params)
          free(vc->exp_params);

        vc->exp_params = vrna_exp_params(md_p);
        break;

      default:
        break;
    }
  }
}


PUBLIC void
vrna_exp_params_subst(vrna_fold_compound_t  *vc,
                      vrna_exp_param_t      *params)
{
  if (vc) {
    if (vc->exp_params)
      free(vc->exp_params);

    if (params) {
      vc->exp_params = vrna_exp_params_copy(params);
    } else {
      switch (vc->type) {
        case VRNA_FC_TYPE_SINGLE:
          vc->exp_params = vrna_exp_params(NULL);
          if (vc->strands > 1)
            vc->exp_params->model_details.min_loop_size = 0;

          break;

        case VRNA_FC_TYPE_COMPARATIVE:
          vc->exp_params = vrna_exp_params_comparative(vc->n_seq, NULL);
          break;

        default:
          break;
      }
    }

    /* fill additional helper arrays for scaling etc. */
    vrna_exp_params_rescale(vc, NULL);
  }
}


PUBLIC void
vrna_exp_params_rescale(vrna_fold_compound_t  *vc,
                        double                *mfe)
{
  vrna_exp_param_t  *pf;
  double            e_per_nt, kT;
  vrna_md_t         *md;

  if (vc) {
    if (!vc->exp_params) {
      switch (vc->type) {
        case VRNA_FC_TYPE_SINGLE:
          vc->exp_params = vrna_exp_params(&(vc->params->model_details));
          break;
        case VRNA_FC_TYPE_COMPARATIVE:
          vc->exp_params = vrna_exp_params_comparative(vc->n_seq, &(vc->params->model_details));
          break;
      }
    } else if (memcmp(&(vc->params->model_details),
                      &(vc->exp_params->model_details),
                      sizeof(vrna_md_t)) != 0) {
      /* make sure that model details are matching */
      (void)vrna_md_copy(&(vc->exp_params->model_details), &(vc->params->model_details));
      /* we probably need some mechanism to check whether DP matrices still match the new model settings! */
    }

    pf = vc->exp_params;
    if (pf) {
      kT  = pf->kT;
      md  = &(pf->model_details);

      if (vc->type == VRNA_FC_TYPE_COMPARATIVE)
        kT /= vc->n_seq;

      /* re-compute scaling factor if necessary */
      if ((mfe) || (pf->pf_scale < 1.)) {
        if (mfe)  /* use largest known Boltzmann factor for scaling */
          e_per_nt = *mfe * 1000. / vc->length;
        else      /* use mean energy for random sequences: 184.3*length cal for scaling */
          e_per_nt = -185 + (pf->temperature - 37.) * 7.27;

        /* apply user-defined scaling factor to allow scaling for unusually stable/unstable structure enembles */
        pf->pf_scale = exp(-(md->sfact * e_per_nt) / kT);
      }

      if (pf->pf_scale < 1.)
        pf->pf_scale = 1.;

      rescale_params(vc);
    }
  }
}


PUBLIC void
vrna_params_prepare(vrna_fold_compound_t  *fc,
                    unsigned int          options)
{
  if (fc) {
    vrna_md_t *md_p;

    /*
     *  every vrna_fold_compound_t must have a vrna_paramt_t structure attached
     *  to it that holds the current model details. So we just use this here as
     *  the reference model
     */
    md_p = &(fc->params->model_details);

    if (options & VRNA_OPTION_PF) {
      /* remove previous parameters if present and they differ from reference model */
      if (fc->exp_params) {
        if (memcmp(md_p, &(fc->exp_params->model_details), sizeof(vrna_md_t)) != 0) {
          free(fc->exp_params);
          fc->exp_params = NULL;
        }
      }

      if (!fc->exp_params)
        fc->exp_params = (fc->type == VRNA_FC_TYPE_SINGLE) ? \
                         vrna_exp_params(md_p) : \
                         vrna_exp_params_comparative(fc->n_seq, md_p);
    }
  }
}


/*
 #####################################
 # BEGIN OF STATIC HELPER FUNCTIONS  #
 #####################################
 */
PRIVATE vrna_param_t *
get_scaled_params(vrna_md_t *md)
{
  unsigned int  i, j, k, l;
  double        tempf;
  vrna_param_t  *params;
  double salt = md->salt;
  double saltStandard = VRNA_MODEL_DEFAULT_SALT;

  params = (vrna_param_t *)vrna_alloc(sizeof(vrna_param_t));

  memset(params->param_file, '\0', 256);
  if (last_parameter_file() != NULL)
    strncpy(params->param_file, last_parameter_file(), 255);

  params->model_details = *md;  /* copy over the model details */
  params->temperature   = md->temperature;
  tempf                 = ((params->temperature + K0) / Tmeasure);

  params->ninio[2]              = RESCALE_dG(ninio37, niniodH, tempf);
  params->lxc                   = lxc37 * tempf;
  params->TripleC               = RESCALE_dG(TripleC37, TripleCdH, tempf);
  params->MultipleCA            = RESCALE_dG(MultipleCA37, MultipleCAdH, tempf);
  params->MultipleCB            = RESCALE_dG(MultipleCB37, MultipleCBdH, tempf);
  params->TerminalAU            = RESCALE_dG(TerminalAU37, TerminalAUdH, tempf);
  params->DuplexInit            = RESCALE_dG(DuplexInit37, DuplexInitdH, tempf);
  params->MLbase                = RESCALE_dG(ML_BASE37, ML_BASEdH, tempf);
  params->MLclosing             = RESCALE_dG(ML_closing37, ML_closingdH, tempf);
  params->gquadLayerMismatch    = RESCALE_dG(GQuadLayerMismatch37, GQuadLayerMismatchH, tempf);
  params->gquadLayerMismatchMax = GQuadLayerMismatchMax;

  for (i = VRNA_GQUAD_MIN_STACK_SIZE; i <= VRNA_GQUAD_MAX_STACK_SIZE; i++)
    for (j = 3 * VRNA_GQUAD_MIN_LINKER_LENGTH; j <= 3 * VRNA_GQUAD_MAX_LINKER_LENGTH; j++) {
      double  GQuadAlpha_T  = RESCALE_dG(GQuadAlpha37, GQuadAlphadH, tempf);
      double  GQuadBeta_T   = RESCALE_dG(GQuadBeta37, GQuadBetadH, tempf);
      double  GQuadEnergy   = (double)GQuadAlpha_T * (double)(i - 1) +
                              (double)GQuadBeta_T * log((double)j - 2.);
      params->gquad[i][j]   = (int)GQuadEnergy;
    }

  for (i = 0; i < 31; i++)
    params->hairpin[i] = RESCALE_dG(hairpin37[i], hairpindH[i], tempf);

  for (i = 0; i <= MIN2(30, MAXLOOP); i++) {
    params->bulge[i]          = RESCALE_dG(bulge37[i], bulgedH[i], tempf);
    params->internal_loop[i]  = RESCALE_dG(internal_loop37[i], internal_loopdH[i], tempf);
  }

  for (; i <= MAXLOOP; i++) {
    params->bulge[i] = params->bulge[30] +
                       (int)(params->lxc * log((double)(i) / 30.));
    params->internal_loop[i] = params->internal_loop[30] +
                               (int)(params->lxc * log((double)(i) / 30.));
    params->SaltLoopDbl[i]   = (salt==saltStandard) ? 0. : vrna_salt_loop(i, salt, saltT, md->backbone_length);
    params->SaltLoop[i]      = (int) (params->SaltLoopDbl[i] + 0.5 - (params->SaltLoopDbl[i]<0));
  }

  for (i = 0; i <= MIN2(31, MAXLOOP+1); i++) {
    params->SaltLoopDbl[i]    = (salt==saltStandard) ? 0. : vrna_salt_loop(i, salt, saltT, md->backbone_length);
    params->SaltLoop[i]       = (int) (params->SaltLoopDbl[i] + 0.5 - (params->SaltLoopDbl[i]<0));
  }

  for (; i <= MAXLOOP; i++) {
    params->SaltLoopDbl[i]   = (salt==saltStandard) ? 0. : vrna_salt_loop(i, salt, saltT, md->backbone_length);
    params->SaltLoop[i]      = (int) (params->SaltLoopDbl[i] + 0.5 - (params->SaltLoopDbl[i]<0));
  }

  for (i = 0; (i * 7) < strlen(Tetraloops); i++)
    params->Tetraloop_E[i] = RESCALE_dG(Tetraloop37[i], TetraloopdH[i], tempf);

  for (i = 0; (i * 5) < strlen(Triloops); i++)
    params->Triloop_E[i] = RESCALE_dG(Triloop37[i], TriloopdH[i], tempf);

  for (i = 0; (i * 9) < strlen(Hexaloops); i++)
    params->Hexaloop_E[i] = RESCALE_dG(Hexaloop37[i], HexaloopdH[i], tempf);

  for (i = 0; i <= NBPAIRS; i++)
    params->MLintern[i] = RESCALE_dG(ML_intern37, ML_interndH, tempf);

  /* stacks    G(T) = H - [H - G(T0)]*T/T0 */
  for (i = 0; i <= NBPAIRS; i++)
    for (j = 0; j <= NBPAIRS; j++)
      params->stack[i][j] = RESCALE_dG(stack37[i][j],
                                       stackdH[i][j],
                                       tempf);

  /* mismatches */
  for (i = 0; i <= NBPAIRS; i++)
    for (j = 0; j < 5; j++)
      for (k = 0; k < 5; k++) {
        int mm;
        params->mismatchI[i][j][k] = RESCALE_dG(mismatchI37[i][j][k],
                                                mismatchIdH[i][j][k],
                                                tempf);
        params->mismatchH[i][j][k] = RESCALE_dG(mismatchH37[i][j][k],
                                                mismatchHdH[i][j][k],
                                                tempf);
        params->mismatch1nI[i][j][k] = RESCALE_dG(mismatch1nI37[i][j][k],
                                                  mismatch1nIdH[i][j][k],
                                                  tempf);
        params->mismatch23I[i][j][k] = RESCALE_dG(mismatch23I37[i][j][k],
                                                  mismatch23IdH[i][j][k],
                                                  tempf);
        if (md->dangles) {
          mm = RESCALE_dG(mismatchM37[i][j][k],
                          mismatchMdH[i][j][k],
                          tempf);
          params->mismatchM[i][j][k]  = (mm > 0) ? 0 : mm;
          mm                          = RESCALE_dG(mismatchExt37[i][j][k],
                                                   mismatchExtdH[i][j][k],
                                                   tempf);
          params->mismatchExt[i][j][k] = (mm > 0) ? 0 : mm;
        } else {
          params->mismatchM[i][j][k] = params->mismatchExt[i][j][k] = 0;
        }
      }

  /* dangles */
  for (i = 0; i <= NBPAIRS; i++)
    for (j = 0; j < 5; j++) {
      int dd;
      dd = RESCALE_dG(dangle5_37[i][j],
                      dangle5_dH[i][j],
                      tempf);
      params->dangle5[i][j] = (dd > 0) ? 0 : dd;  /* must be <= 0 */
      dd                    = RESCALE_dG(dangle3_37[i][j],
                                         dangle3_dH[i][j],
                                         tempf);
      params->dangle3[i][j] = (dd > 0) ? 0 : dd;  /* must be <= 0 */
    }

  /* internal 1x1 loops */
  for (i = 0; i <= NBPAIRS; i++)
    for (j = 0; j <= NBPAIRS; j++)
      for (k = 0; k < 5; k++)
        for (l = 0; l < 5; l++)
          params->int11[i][j][k][l] = RESCALE_dG(int11_37[i][j][k][l],
                                                 int11_dH[i][j][k][l],
                                                 tempf);

  /* internal 2x1 loops */
  for (i = 0; i <= NBPAIRS; i++)
    for (j = 0; j <= NBPAIRS; j++)
      for (k = 0; k < 5; k++)
        for (l = 0; l < 5; l++) {
          int m;
          for (m = 0; m < 5; m++)
            params->int21[i][j][k][l][m] = RESCALE_dG(int21_37[i][j][k][l][m],
                                                      int21_dH[i][j][k][l][m],
                                                      tempf);
        }

  /* internal 2x2 loops */
  for (i = 0; i <= NBPAIRS; i++)
    for (j = 0; j <= NBPAIRS; j++)
      for (k = 0; k < 5; k++)
        for (l = 0; l < 5; l++) {
          int m, n;
          for (m = 0; m < 5; m++)
            for (n = 0; n < 5; n++)
              params->int22[i][j][k][l][m][n] = RESCALE_dG(int22_37[i][j][k][l][m][n],
                                                           int22_dH[i][j][k][l][m][n],
                                                           tempf);
        }

  strncpy(params->Tetraloops, Tetraloops, 281);
  strncpy(params->Triloops, Triloops, 241);
  strncpy(params->Hexaloops, Hexaloops, 361);

  /* Salt correction for stack and multiloop */
  params->SaltStack = (salt==saltStandard) ? 0 : vrna_salt_stack(salt, saltT, md->helical_rise);
  if (salt == saltStandard)
    params->SaltMLbase = params->SaltMLclosing = 0;
  else
    vrna_salt_ml(params->SaltLoopDbl, md->saltMLLower, md->saltMLUpper, &params->SaltMLbase, &params->SaltMLclosing);
  
  params->MLclosing += params->SaltMLbase;
  params->MLclosing += params->SaltMLclosing;
  params->MLbase += params->SaltMLbase;
  for (i = 0; i <= NBPAIRS; i++)
    params->MLintern[i] += params->SaltMLbase;

  params->SaltDPXInit = 0;
  if (salt != saltStandard)
  {
    if (md->saltDPXInit != VRNA_MODEL_DEFAULT_SALT_DPXINIT)
      params->SaltDPXInit = md->saltDPXInit;
    else if (md->saltDPXInit)
      params->SaltDPXInit = vrna_salt_duplex_init(md);
  }
  params->DuplexInit += params->SaltDPXInit;

  params->id = ++id;
  return params;
}


PRIVATE vrna_exp_param_t *
get_scaled_exp_params(vrna_md_t *md,
                      double    pfs)
{
  unsigned int      i, j, k, l;
  int               pf_smooth;
  double            kT, TT;
  double            GT;
  double            salt, saltStandard;
  vrna_exp_param_t  *pf;

  pf = (vrna_exp_param_t *)vrna_alloc(sizeof(vrna_exp_param_t));

  memset(pf->param_file, '\0', 256);
  if (last_parameter_file() != NULL)
    strncpy(pf->param_file, last_parameter_file(), 255);

  pf->model_details = *md;
  pf->temperature   = md->temperature;
  pf->alpha         = md->betaScale;
  pf->kT            = kT = md->betaScale * (md->temperature + K0) * GASCONST; /* kT in cal/mol  */
  pf->pf_scale      = pfs;
  pf_smooth         = md->pf_smooth;
  TT                = (md->temperature + K0) / (Tmeasure);
  salt = md->salt;
  saltStandard = VRNA_MODEL_DEFAULT_SALT;

  pf->lxc                   = lxc37 * TT;
  pf->expDuplexInit         = RESCALE_BF(DuplexInit37, DuplexInitdH, TT, kT);
  pf->expTermAU             = RESCALE_BF(TerminalAU37, TerminalAUdH, TT, kT);
  pf->expMLbase             = RESCALE_BF(ML_BASE37, ML_BASEdH, TT, kT);
  pf->expMLclosing          = RESCALE_BF(ML_closing37, ML_closingdH, TT, kT);
  pf->expgquadLayerMismatch = RESCALE_BF(GQuadLayerMismatch37, GQuadLayerMismatchH, TT, kT);
  pf->gquadLayerMismatchMax = GQuadLayerMismatchMax;

  for (i = VRNA_GQUAD_MIN_STACK_SIZE; i <= VRNA_GQUAD_MAX_STACK_SIZE; i++)
    for (j = 3 * VRNA_GQUAD_MIN_LINKER_LENGTH; j <= 3 * VRNA_GQUAD_MAX_LINKER_LENGTH; j++) {
      double  GQuadAlpha_T  = RESCALE_dG(GQuadAlpha37, GQuadAlphadH, TT);
      double  GQuadBeta_T   = RESCALE_dG(GQuadBeta37, GQuadBetadH, TT);
      double  GQuadEnergy   = (double)GQuadAlpha_T * (double)(i - 1) +
                              (double)GQuadBeta_T * log((double)j - 2.);
      pf->expgquad[i][j] = exp(-TRUNC_MAYBE((double)((int)GQuadEnergy)) * 10. / kT);
    }

  /* loop energies: hairpins, bulges, internal, mulit-loops */
  for (i = 0; i < 31; i++)
    pf->exphairpin[i] = RESCALE_BF(hairpin37[i], hairpindH[i], TT, kT);

  for (i = 0; i <= MIN2(30, MAXLOOP); i++) {
    pf->expbulge[i]     = RESCALE_BF(bulge37[i], bulgedH[i], TT, kT);
    pf->expinternal[i]  = RESCALE_BF(internal_loop37[i], internal_loopdH[i], TT, kT);
  }

  if (salt==saltStandard) {
    for (i = 0; i <= MIN2(30, MAXLOOP); i++) {
      pf->SaltLoopDbl[i]   = 0.;
      pf->expSaltLoop[i]   = 1.;
    }

    for (i = 31; i <= MAXLOOP; i++) {
      pf->SaltLoopDbl[i] = 0.;
      pf->expSaltLoop[i] = 1.;
    }
  } else {
    for (i = 0; i <= MIN2(30, MAXLOOP); i++) {
      pf->SaltLoopDbl[i]   = vrna_salt_loop(i, salt, saltT, md->backbone_length);
      int saltLoop = (int) (pf->SaltLoopDbl[i] + 0.5 - (pf->SaltLoopDbl[i]<0));
      pf->expSaltLoop[i]   = exp(-saltLoop * 10. / kT);
    }

    for (i = 31; i <= MAXLOOP; i++) {
      pf->SaltLoopDbl[i] = vrna_salt_loop(i, salt, saltT, md->backbone_length);
      int saltLoop = (int) (pf->SaltLoopDbl[i] + 0.5 - (pf->SaltLoopDbl[i]<0));
      pf->expSaltLoop[i] = exp(-saltLoop * 10. / kT);
    }
  }

  /* special case of size 2 internal loops (single mismatch) */
  if (james_rule)
    pf->expinternal[2] = exp(-80 * 10. / kT);

  GT = RESCALE_dG(bulge37[30],
                  bulgedH[30],
                  TT);
  for (i = 31; i <= MAXLOOP; i++)
    pf->expbulge[i] = exp(-TRUNC_MAYBE(GT + (pf->lxc * log(i / 30.))) * 10. / kT);

  GT = RESCALE_dG(internal_loop37[30],
                  internal_loopdH[30],
                  TT);
  for (i = 31; i <= MAXLOOP; i++)
    pf->expinternal[i] = exp(-TRUNC_MAYBE(GT + (pf->lxc * log(i / 30.))) * 10. / kT);

  GT = RESCALE_dG(ninio37, niniodH, TT);
  for (j = 0; j <= MAXLOOP; j++)
    pf->expninio[2][j] = exp(-MIN2(MAX_NINIO, j * TRUNC_MAYBE(GT)) * 10. / kT);

  for (i = 0; (i * 7) < strlen(Tetraloops); i++)
    pf->exptetra[i] = RESCALE_BF(Tetraloop37[i], TetraloopdH[i], TT, kT);

  for (i = 0; (i * 5) < strlen(Triloops); i++)
    pf->exptri[i] = RESCALE_BF(Triloop37[i], TriloopdH[i], TT, kT);

  for (i = 0; (i * 9) < strlen(Hexaloops); i++)
    pf->exphex[i] = RESCALE_BF(Hexaloop37[i], HexaloopdH[i], TT, kT);

  for (i = 0; i <= NBPAIRS; i++)
    pf->expMLintern[i] = RESCALE_BF(ML_intern37, ML_interndH, TT, kT);

  /* if dangles==0 just set their energy to 0,
   * don't let dangle energies become > 0 (at large temps),
   * but make sure go smoothly to 0                        */
  for (i = 0; i <= NBPAIRS; i++)
    for (j = 0; j <= 4; j++) {
      if (md->dangles) {
        pf->expdangle5[i][j]  = RESCALE_BF_SMOOTH(dangle5_37[i][j], dangle5_dH[i][j], TT, kT);
        pf->expdangle3[i][j]  = RESCALE_BF_SMOOTH(dangle3_37[i][j], dangle3_dH[i][j], TT, kT);
      } else {
        pf->expdangle3[i][j] = pf->expdangle5[i][j] = 1;
      }
    }

  /* stacking energies */
  for (i = 0; i <= NBPAIRS; i++)
    for (j = 0; j <= NBPAIRS; j++)
      pf->expstack[i][j] = RESCALE_BF(stack37[i][j], stackdH[i][j], TT, kT);

  /* mismatch energies */
  for (i = 0; i <= NBPAIRS; i++)
    for (j = 0; j < 5; j++)
      for (k = 0; k < 5; k++) {
        pf->expmismatchI[i][j][k] = RESCALE_BF(mismatchI37[i][j][k],
                                               mismatchIdH[i][j][k],
                                               TT,
                                               kT);
        pf->expmismatch1nI[i][j][k] = RESCALE_BF(mismatch1nI37[i][j][k],
                                                 mismatch1nIdH[i][j][k],
                                                 TT,
                                                 kT);
        pf->expmismatchH[i][j][k] = RESCALE_BF(mismatchH37[i][j][k],
                                               mismatchHdH[i][j][k],
                                               TT,
                                               kT);
        pf->expmismatch23I[i][j][k] = RESCALE_BF(mismatch23I37[i][j][k],
                                                 mismatch23IdH[i][j][k],
                                                 TT,
                                                 kT);

        if (md->dangles) {
          pf->expmismatchM[i][j][k] = RESCALE_BF_SMOOTH(mismatchM37[i][j][k],
                                                        mismatchMdH[i][j][k],
                                                        TT,
                                                        kT);
          pf->expmismatchExt[i][j][k] = RESCALE_BF_SMOOTH(mismatchExt37[i][j][k],
                                                          mismatchExtdH[i][j][k],
                                                          TT,
                                                          kT);
        } else {
          pf->expmismatchM[i][j][k] = pf->expmismatchExt[i][j][k] = 1.;
        }
      }

  /* internal lops of length 2 */
  for (i = 0; i <= NBPAIRS; i++)
    for (j = 0; j <= NBPAIRS; j++)
      for (k = 0; k < 5; k++)
        for (l = 0; l < 5; l++) {
          pf->expint11[i][j][k][l] = RESCALE_BF(int11_37[i][j][k][l],
                                                int11_dH[i][j][k][l],
                                                TT,
                                                kT);
        }

  /* internal 2x1 loops */
  for (i = 0; i <= NBPAIRS; i++)
    for (j = 0; j <= NBPAIRS; j++)
      for (k = 0; k < 5; k++)
        for (l = 0; l < 5; l++) {
          int m;
          for (m = 0; m < 5; m++) {
            pf->expint21[i][j][k][l][m] = RESCALE_BF(int21_37[i][j][k][l][m],
                                                     int21_dH[i][j][k][l][m],
                                                     TT,
                                                     kT);
          }
        }

  /* internal 2x2 loops */
  for (i = 0; i <= NBPAIRS; i++)
    for (j = 0; j <= NBPAIRS; j++)
      for (k = 0; k < 5; k++)
        for (l = 0; l < 5; l++) {
          int m, n;
          for (m = 0; m < 5; m++)
            for (n = 0; n < 5; n++) {
              pf->expint22[i][j][k][l][m][n] = RESCALE_BF(int22_37[i][j][k][l][m][n],
                                                          int22_dH[i][j][k][l][m][n],
                                                          TT,
                                                          kT);
            }
        }

  strncpy(pf->Tetraloops, Tetraloops, 281);
  strncpy(pf->Triloops, Triloops, 241);
  strncpy(pf->Hexaloops, Hexaloops, 361);

  /* Salt correction for stack and multiloop */
  pf->SaltMLbase = pf->SaltMLclosing = pf->SaltDPXInit = 0.;

  if (salt==saltStandard) {
    pf->expSaltStack = 1.;
  } else {
    pf->expSaltStack = exp(- vrna_salt_stack(salt, saltT, md->helical_rise) * 10. / kT);
    vrna_salt_ml(pf->SaltLoopDbl, md->saltMLLower, md->saltMLUpper, &pf->SaltMLbase, &pf->SaltMLclosing);

    if (md->saltDPXInit != VRNA_MODEL_DEFAULT_SALT_DPXINIT)
      pf->SaltDPXInit = md->saltDPXInit;
    else if (md->saltDPXInit)
      pf->SaltDPXInit = vrna_salt_duplex_init(md);

    pf->expMLclosing *= exp(- pf->SaltMLbase * 10. / kT);
    pf->expMLclosing *= exp(- pf->SaltMLclosing * 10. / kT);
    pf->expMLbase *= exp(- pf->SaltMLbase * 10. / kT);
    for (i = 0; i <= NBPAIRS; i++)
      pf->expMLintern[i] *= exp(- pf->SaltMLbase * 10. / kT);

    pf->expDuplexInit *= exp(- pf->SaltDPXInit*10. / kT);
  }

  return pf;
}


PRIVATE vrna_exp_param_t *
get_exp_params_ali(vrna_md_t    *md,
                   unsigned int n_seq,
                   double       pfs)
{
  /* scale energy parameters and pre-calculate Boltzmann weights */
  unsigned int      i, j, k, l;
  int               pf_smooth;
  double            kTn, TT;
  double            GT;
  double            salt, saltStandard;
  vrna_exp_param_t  *pf;

  pf                = (vrna_exp_param_t *)vrna_alloc(sizeof(vrna_exp_param_t));
  pf->model_details = *md;
  pf->alpha         = md->betaScale;
  pf->temperature   = md->temperature;
  pf->pf_scale      = pfs;
  pf->kT            = kTn = ((double)n_seq) * md->betaScale * (md->temperature + K0) * GASCONST; /* kT in cal/mol  */
  pf_smooth         = md->pf_smooth;
  TT                = (md->temperature + K0) / (Tmeasure);
  salt = md->salt;
  saltStandard = VRNA_MODEL_DEFAULT_SALT;

  pf->lxc                   = lxc37 * TT;
  pf->expDuplexInit         = RESCALE_BF(DuplexInit37, DuplexInitdH, TT, kTn);
  pf->expTermAU             = RESCALE_BF(TerminalAU37, TerminalAUdH, TT, kTn);
  pf->expMLbase             = RESCALE_BF(ML_BASE37, ML_BASEdH, TT, kTn / n_seq);
  pf->expMLclosing          = RESCALE_BF(ML_closing37, ML_closingdH, TT, kTn);
  pf->expgquadLayerMismatch = RESCALE_BF(GQuadLayerMismatch37, GQuadLayerMismatchH, TT, kTn);
  pf->gquadLayerMismatchMax = GQuadLayerMismatchMax;

  for (i = VRNA_GQUAD_MIN_STACK_SIZE; i <= VRNA_GQUAD_MAX_STACK_SIZE; i++)
    for (j = 3 * VRNA_GQUAD_MIN_LINKER_LENGTH; j <= 3 * VRNA_GQUAD_MAX_LINKER_LENGTH; j++) {
      double  GQuadAlpha_T  = RESCALE_dG(GQuadAlpha37, GQuadAlphadH, TT);
      double  GQuadBeta_T   = RESCALE_dG(GQuadBeta37, GQuadBetadH, TT);
      GT = ((double)GQuadAlpha_T) * ((double)(i - 1)) + ((double)GQuadBeta_T) *
           log(((double)j) - 2.);
      pf->expgquad[i][j] = exp(-TRUNC_MAYBE(GT) * 10. / kTn);
    }

  /* loop energies: hairpins, bulges, internal, mulit-loops */
  for (i = 0; i < 31; i++)
    pf->exphairpin[i] = RESCALE_BF(hairpin37[i], hairpindH[i], TT, kTn);
  /*add penalty for too short hairpins*/
  for (i = 0; i < 3; i++) {
    GT                = 600 /*Penalty*/ * TT;
    pf->exphairpin[i] = exp(-GT * 10. / kTn);
  }

  for (i = 0; i <= MIN2(30, MAXLOOP); i++) {
    pf->expbulge[i]     = RESCALE_BF(bulge37[i], bulgedH[i], TT, kTn);
    pf->expinternal[i]  = RESCALE_BF(internal_loop37[i], internal_loopdH[i], TT, kTn);

    pf->SaltLoopDbl[i]   = (salt==saltStandard) ? 0. : vrna_salt_loop(i, salt, saltT, md->backbone_length);
    int saltLoop = (int) (pf->SaltLoopDbl[i] + 0.5 - (pf->SaltLoopDbl[i]<0));
    pf->expSaltLoop[i]   = exp(-saltLoop * 10. / kTn);
  }

  if (salt==saltStandard) {
    for (i = 0; i <= MIN2(30, MAXLOOP); i++) {
      pf->SaltLoopDbl[i]  = 0.;
      pf->expSaltLoop[i]  = 1.;
    }
    for (i = 31; i <= MAXLOOP; i++) {
      pf->SaltLoopDbl[i] = 0.;
      pf->expSaltLoop[i] = 1.;
    }
  } else {
    for (i = 0; i <= MIN2(30, MAXLOOP); i++) {
      pf->SaltLoopDbl[i]  = vrna_salt_loop(i, salt, saltT, md->backbone_length);
      int saltLoop = (int) (pf->SaltLoopDbl[i] + 0.5 - (pf->SaltLoopDbl[i]<0));
      pf->expSaltLoop[i]   = exp(-saltLoop * 10. / kTn);
    }

    for (i = 31; i <= MAXLOOP; i++) {
      pf->SaltLoopDbl[i] = vrna_salt_loop(i, salt, saltT, md->backbone_length);
      int saltLoop = (int) (pf->SaltLoopDbl[i] + 0.5 - (pf->SaltLoopDbl[i]<0));
      pf->expSaltLoop[i] = exp(-saltLoop * 10. / kTn);
    }
  }

  /* special case of size 2 internal loops (single mismatch) */
  if (james_rule)
    pf->expinternal[2] = exp(-80 * 10. / kTn);

  GT = RESCALE_dG(bulge37[30], bulgedH[30], TT);
  for (i = 31; i <= MAXLOOP; i++)
    pf->expbulge[i] = exp(-(GT + (pf->lxc * log(i / 30.))) * 10. / kTn);

  GT = RESCALE_dG(internal_loop37[30], internal_loopdH[30], TT);
  for (i = 31; i <= MAXLOOP; i++)
    pf->expinternal[i] = exp(-(GT + (pf->lxc * log(i / 30.))) * 10. / kTn);

  GT = RESCALE_dG(ninio37, niniodH, TT);
  for (j = 0; j <= MAXLOOP; j++)
    pf->expninio[2][j] = exp(-MIN2(MAX_NINIO, j * GT) * 10. / kTn);

  for (i = 0; (i * 7) < strlen(Tetraloops); i++)
    pf->exptetra[i] = RESCALE_BF(Tetraloop37[i], TetraloopdH[i], TT, kTn);

  for (i = 0; (i * 5) < strlen(Triloops); i++)
    pf->exptri[i] = RESCALE_BF(Triloop37[i], TriloopdH[i], TT, kTn);

  for (i = 0; (i * 9) < strlen(Hexaloops); i++)
    pf->exphex[i] = RESCALE_BF(Hexaloop37[i], HexaloopdH[i], TT, kTn);

  for (i = 0; i <= NBPAIRS; i++)
    /* includes AU penalty */
    pf->expMLintern[i] = RESCALE_BF(ML_intern37, ML_interndH, TT, kTn);

  /* if dangle_model==0 just set their energy to 0,
   * don't let dangle energies become > 0 (at large temps),
   * but make sure go smoothly to 0                        */
  for (i = 0; i <= NBPAIRS; i++)
    for (j = 0; j <= 4; j++) {
      if (md->dangles) {
        pf->expdangle5[i][j] = RESCALE_BF_SMOOTH(dangle5_37[i][j],
                                                 dangle5_dH[i][j],
                                                 TT,
                                                 kTn);
        pf->expdangle3[i][j] = RESCALE_BF_SMOOTH(dangle3_37[i][j],
                                                 dangle3_dH[i][j],
                                                 TT,
                                                 kTn);
      } else {
        pf->expdangle3[i][j] = pf->expdangle5[i][j] = 1;
      }
    }

  /* stacking energies */
  for (i = 0; i <= NBPAIRS; i++)
    for (j = 0; j <= NBPAIRS; j++) {
      pf->expstack[i][j] = RESCALE_BF(stack37[i][j],
                                      stackdH[i][j],
                                      TT,
                                      kTn);
    }

  /* mismatch energies */
  for (i = 0; i <= NBPAIRS; i++)
    for (j = 0; j < 5; j++)
      for (k = 0; k < 5; k++) {
        pf->expmismatchI[i][j][k] = RESCALE_BF(mismatchI37[i][j][k],
                                               mismatchIdH[i][j][k],
                                               TT,
                                               kTn);
        pf->expmismatch1nI[i][j][k] = RESCALE_BF(mismatch1nI37[i][j][k],
                                                 mismatch1nIdH[i][j][k],
                                                 TT,
                                                 kTn);
        pf->expmismatchH[i][j][k] = RESCALE_BF(mismatchH37[i][j][k],
                                               mismatchHdH[i][j][k],
                                               TT,
                                               kTn);
        pf->expmismatch23I[i][j][k] = RESCALE_BF(mismatch23I37[i][j][k],
                                                 mismatch23IdH[i][j][k],
                                                 TT,
                                                 kTn);

        if (md->dangles) {
          pf->expmismatchM[i][j][k] = RESCALE_BF_SMOOTH(mismatchM37[i][j][k],
                                                        mismatchMdH[i][j][k],
                                                        TT,
                                                        kTn);
          pf->expmismatchExt[i][j][k] = RESCALE_BF_SMOOTH(mismatchExt37[i][j][k],
                                                          mismatchExtdH[i][j][k],
                                                          TT,
                                                          kTn);
        } else {
          pf->expmismatchM[i][j][k] = pf->expmismatchExt[i][j][k] = 1.;
        }
      }

  /* internal lops of length 2 */
  for (i = 0; i <= NBPAIRS; i++)
    for (j = 0; j <= NBPAIRS; j++)
      for (k = 0; k < 5; k++)
        for (l = 0; l < 5; l++) {
          pf->expint11[i][j][k][l] = RESCALE_BF(int11_37[i][j][k][l],
                                                int11_dH[i][j][k][l],
                                                TT,
                                                kTn);
        }

  /* internal 2x1 loops */
  for (i = 0; i <= NBPAIRS; i++)
    for (j = 0; j <= NBPAIRS; j++)
      for (k = 0; k < 5; k++)
        for (l = 0; l < 5; l++) {
          int m;
          for (m = 0; m < 5; m++) {
            pf->expint21[i][j][k][l][m] = RESCALE_BF(int21_37[i][j][k][l][m],
                                                     int21_dH[i][j][k][l][m],
                                                     TT,
                                                     kTn);
          }
        }

  /* internal 2x2 loops */
  for (i = 0; i <= NBPAIRS; i++)
    for (j = 0; j <= NBPAIRS; j++)
      for (k = 0; k < 5; k++)
        for (l = 0; l < 5; l++) {
          int m, n;
          for (m = 0; m < 5; m++)
            for (n = 0; n < 5; n++) {
              pf->expint22[i][j][k][l][m][n] = RESCALE_BF(int22_37[i][j][k][l][m][n],
                                                          int22_dH[i][j][k][l][m][n],
                                                          TT,
                                                          kTn);
            }
        }

  strncpy(pf->Tetraloops, Tetraloops, 281);
  strncpy(pf->Triloops, Triloops, 241);
  strncpy(pf->Hexaloops, Hexaloops, 361);

  /* Salt correction for stack and multiloop */
  pf->SaltMLbase = pf->SaltMLclosing = pf->SaltDPXInit = 0.;

  if (salt==saltStandard) {
    pf->expSaltStack = 1.;
  } else {
    pf->expSaltStack = exp(- vrna_salt_stack(salt, saltT, md->helical_rise) * 10. / kTn);

    vrna_salt_ml(pf->SaltLoopDbl, md->saltMLLower, md->saltMLUpper, &pf->SaltMLbase, &pf->SaltMLclosing);

    if (md->saltDPXInit != VRNA_MODEL_DEFAULT_SALT_DPXINIT)
      pf->SaltDPXInit = md->saltDPXInit;
    else if (md->saltDPXInit)
      pf->SaltDPXInit = vrna_salt_duplex_init(md);

    pf->expMLclosing *= exp(- pf->SaltMLbase * 10. / kTn);
    pf->expMLclosing *= exp(- pf->SaltMLclosing * 10. / kTn);
    pf->expMLbase *= exp(- pf->SaltMLbase * 10. / kTn);
    for (i = 0; i <= NBPAIRS; i++)
      pf->expMLintern[i] *= exp(- pf->SaltMLbase * 10. / kTn);

    pf->expDuplexInit *= exp(- pf->SaltDPXInit*10. / kTn);
  }

  return pf;
}


PRIVATE void
rescale_params(vrna_fold_compound_t *vc)
{
  unsigned int      i;
  vrna_exp_param_t  *pf = vc->exp_params;
  vrna_mx_pf_t      *m  = vc->exp_matrices;

  if (m && pf) {
    m->scale[0]     = 1.;
    m->scale[1]     = (FLT_OR_DBL)(1. / pf->pf_scale);
    m->expMLbase[0] = 1;
    m->expMLbase[1] = (FLT_OR_DBL)(pf->expMLbase / pf->pf_scale);
    for (i = 2; i <= vc->length; i++) {
      m->scale[i]     = m->scale[i / 2] * m->scale[i - (i / 2)];
      m->expMLbase[i] = (FLT_OR_DBL)pow(pf->expMLbase, (double)i) * m->scale[i];
    }
  }
}


#ifndef VRNA_DISABLE_BACKWARD_COMPATIBILITY

/*
 *###########################################
 *# deprecated functions below              #
 *###########################################
 */
PUBLIC vrna_param_t *
scale_parameters(void)
{
  vrna_md_t md;

  set_model_details(&md);
  return vrna_params(&md);
}


PUBLIC vrna_param_t *
get_scaled_parameters(double    temp,
                      vrna_md_t md)
{
  md.temperature = temp;
  return get_scaled_params(&md);
}


PUBLIC vrna_exp_param_t *
get_boltzmann_factors(double    temp,
                      double    betaScale,
                      vrna_md_t md,
                      double    pfs)
{
  md.temperature  = temp;
  md.betaScale    = betaScale;
  pf_scale        = pfs;

  return get_scaled_exp_params(&md, pfs);
}


PUBLIC vrna_exp_param_t *
get_scaled_pf_parameters(void)
{
  vrna_md_t         md;
  vrna_exp_param_t  *pf;

  set_model_details(&md);

  pf            = vrna_exp_params(&md);
  pf->pf_scale  = pf_scale;

  return pf;
}


PUBLIC vrna_exp_param_t *
get_boltzmann_factors_ali(unsigned int  n_seq,
                          double        temp,
                          double        betaScale,
                          vrna_md_t     md,
                          double        pfs)
{
  md.temperature  = temp;
  md.betaScale    = betaScale;
  pf_scale        = pfs;

  return get_exp_params_ali(&md, n_seq, pfs);
}


PUBLIC vrna_exp_param_t *
get_scaled_alipf_parameters(unsigned int n_seq)
{
  vrna_md_t md;

  set_model_details(&md);

  return get_exp_params_ali(&md, n_seq, pf_scale);
}


PUBLIC vrna_exp_param_t *
get_boltzmann_factor_copy(vrna_exp_param_t *par)
{
  return vrna_exp_params_copy(par);
}


PUBLIC vrna_param_t *
get_parameter_copy(vrna_param_t *par)
{
  return vrna_params_copy(par);
}


PUBLIC vrna_param_t *
copy_parameters(void)
{
  vrna_param_t *copy;

  if (p.id != id) {
    vrna_md_t md;
    set_model_details(&md);
    return vrna_params(&md);
  } else {
    copy = (vrna_param_t *)vrna_alloc(sizeof(vrna_param_t));
    memcpy(copy, &p, sizeof(vrna_param_t));
  }

  return copy;
}


PUBLIC vrna_param_t *
set_parameters(vrna_param_t *dest)
{
  memcpy(&p, dest, sizeof(vrna_param_t));
  return &p;
}


PUBLIC vrna_exp_param_t *
copy_pf_param(void)
{
  vrna_exp_param_t *copy;

  if (pf.id != pf_id) {
    vrna_md_t md;
    set_model_details(&md);
    copy            = vrna_exp_params(&md);
    copy->pf_scale  = pf_scale;
    return copy;
  } else {
    copy = (vrna_exp_param_t *)vrna_alloc(sizeof(vrna_exp_param_t));
    memcpy(copy, &pf, sizeof(vrna_exp_param_t));
  }

  return copy;
}


PUBLIC vrna_exp_param_t *
set_pf_param(vrna_param_t *dest)
{
  memcpy(&pf, dest, sizeof(vrna_exp_param_t));
  return &pf;
}


PUBLIC vrna_exp_param_t *
scale_pf_parameters(void)
{
  vrna_md_t         md;
  vrna_exp_param_t  *pf;

  set_model_details(&md);

  pf            = vrna_exp_params(&md);
  pf->pf_scale  = pf_scale;

  return pf;
}


#endif
