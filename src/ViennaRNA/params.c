/*

                  c Ivo Hofacker

                  Vienna RNA package
*/

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "ViennaRNA/energy_par.h"
#include "ViennaRNA/fold_vars.h"
#include "ViennaRNA/utils.h"
#include "ViennaRNA/params.h"
/**
*** \file params.c
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
#define SMOOTH(X) ((X)/SCALE<-1.2283697)?0:(((X)/SCALE>0.8660254)?(X):\
          SCALE*0.38490018*(sin((X)/SCALE-0.34242663)+1)*(sin((X)/SCALE-0.34242663)+1))

/* #define SMOOTH(X) ((X)<0 ? 0 : (X)) */

/*
#################################
# PRIVATE VARIABLES             #
#################################
*/
PRIVATE vrna_param_t p;
PRIVATE int id=-1;
/* variables for partition function */
PRIVATE vrna_exp_param_t pf;
PRIVATE int pf_id=-1;

#ifdef _OPENMP
#pragma omp threadprivate(id, pf_id)
#endif

/*
#################################
# PRIVATE FUNCTION DECLARATIONS #
#################################
*/

PRIVATE vrna_param_t      *get_scaled_params(vrna_md_t *md);
PRIVATE vrna_exp_param_t  *get_scaled_exp_params(vrna_md_t *md, double pfs);
PRIVATE vrna_exp_param_t  *get_exp_params_ali(vrna_md_t *md, unsigned int n_seq, double pfs);
PRIVATE void              rescale_params(vrna_fold_compound_t *vc);

/*
#################################
# BEGIN OF FUNCTION DEFINITIONS #
#################################
*/

PUBLIC vrna_param_t *
vrna_params(vrna_md_t *md){

  if(md){
    return get_scaled_params(md);
  } else {
    vrna_md_t md;
    vrna_md_set_default(&md);
    return get_scaled_params(&md);
  }
}

PUBLIC vrna_exp_param_t *
vrna_exp_params(vrna_md_t *md){

  if(md){
    return  get_scaled_exp_params(md, -1.);
  } else {
    vrna_md_t md;
    vrna_md_set_default(&md);
    return get_scaled_exp_params(&md, -1.);
  }
}

PUBLIC vrna_exp_param_t *
vrna_exp_params_comparative(unsigned int n_seq, vrna_md_t *md){

  if(md){
    return  get_exp_params_ali(md, n_seq, -1.);
  } else {
    vrna_md_t md;
    vrna_md_set_default(&md);
    return get_exp_params_ali(&md, n_seq, -1.);
  }
}

PUBLIC vrna_param_t *
vrna_params_copy(vrna_param_t *par){

  vrna_param_t *copy = NULL;
  if(par){
    copy = (vrna_param_t *) vrna_alloc(sizeof(vrna_param_t));
    memcpy(copy, par, sizeof(vrna_param_t));
  }
  return copy;
}

PUBLIC vrna_exp_param_t *
vrna_exp_params_copy(vrna_exp_param_t *par){

  vrna_exp_param_t *copy = NULL;
  if(par){
    copy = (vrna_exp_param_t *) vrna_alloc(sizeof(vrna_exp_param_t));
    memcpy(copy, par, sizeof(vrna_exp_param_t));
  }
  return copy;
}

PUBLIC void
vrna_params_subst( vrna_fold_compound_t *vc,
                    vrna_param_t *parameters){

  if(vc){
    if(vc->params)
      free(vc->params);
    if(parameters){
      vc->params = vrna_params_copy(parameters);
    } else {
      switch(vc->type){
        case VRNA_FC_TYPE_SINGLE:     /* fall through */

        case VRNA_FC_TYPE_COMPARATIVE:  vc->params = vrna_params(NULL);
                                      break;

        default:                      break;
      }
    }
  }
}

PUBLIC void
vrna_params_reset(vrna_fold_compound_t *vc,
                  vrna_md_t *md_p){

  if(vc){
    switch(vc->type){
      case VRNA_FC_TYPE_SINGLE:     /* fall through */

      case VRNA_FC_TYPE_COMPARATIVE:  if(vc->params)
                                      free(vc->params);
                                    vc->params = vrna_params(md_p);
                                    break;

      default:                      break;
    }
  }
}

PUBLIC void
vrna_exp_params_reset(vrna_fold_compound_t *vc,
                      vrna_md_t *md_p){

  if(vc){
    switch(vc->type){
      case VRNA_FC_TYPE_SINGLE:     /* fall through */

      case VRNA_FC_TYPE_COMPARATIVE:  if(vc->exp_params)
                                      free(vc->exp_params);
                                    vc->exp_params = vrna_exp_params(md_p);
                                    break;

      default:                      break;
    }
  }
}

PUBLIC void
vrna_exp_params_subst(vrna_fold_compound_t *vc,
                      vrna_exp_param_t *params){

  if(vc){
    if(vc->exp_params)
      free(vc->exp_params);
    if(params){
      vc->exp_params = vrna_exp_params_copy(params);
    } else {
      switch(vc->type){
        case VRNA_FC_TYPE_SINGLE:     vc->exp_params = vrna_exp_params(NULL);
                                      if(vc->cutpoint > 0)
                                        vc->exp_params->model_details.min_loop_size = 0;
                                      break;

        case VRNA_FC_TYPE_COMPARATIVE:  vc->exp_params = vrna_exp_params_comparative(vc->n_seq, NULL);
                                      break;

        default:                      break;
      }
    }
    /* fill additional helper arrays for scaling etc. */
    vrna_exp_params_rescale(vc, NULL);
  }
}

PUBLIC void
vrna_exp_params_rescale(vrna_fold_compound_t *vc,
                        double *mfe){

  vrna_exp_param_t  *pf;
  double            kT;
  vrna_md_t         *md;

  if(vc){
    
    if(!vc->exp_params){
      switch(vc->type){
        case VRNA_FC_TYPE_SINGLE:
          vc->exp_params = vrna_exp_params(&(vc->params->model_details));
          break;
        case VRNA_FC_TYPE_COMPARATIVE:
          vc->exp_params = vrna_exp_params_comparative(vc->n_seq, &(vc->params->model_details));
          break;
      }
    }

    pf = vc->exp_params;
    if(pf){
      kT = pf->kT;
      md = &(pf->model_details);

      if(vc->type == VRNA_FC_TYPE_COMPARATIVE)
        kT /= vc->n_seq;

      if(mfe){
        kT /= 1000.;
        pf->pf_scale = exp(-(md->sfact * *mfe)/ kT / vc->length);
      } else if(pf->pf_scale < 1.){  /* mean energy for random sequences: 184.3*length cal */
        pf->pf_scale = exp(-(-185+(pf->temperature-37.)*7.27)/kT);
      }
      if(pf->pf_scale < 1.)
        pf->pf_scale = 1.;
      rescale_params(vc);
    }
  }
}

/*
#####################################
# BEGIN OF STATIC HELPER FUNCTIONS  #
#####################################
*/

PRIVATE vrna_param_t *
get_scaled_params(vrna_md_t *md){

  unsigned int i,j,k,l;
  double tempf;
  vrna_param_t *params;

  params  = (vrna_param_t *)vrna_alloc(sizeof(vrna_param_t));

  params->model_details = *md;  /* copy over the model details */
  params->temperature   = md->temperature;
  tempf                 = ((params->temperature+K0)/Tmeasure);

  for(i = VRNA_GQUAD_MIN_STACK_SIZE; i <= VRNA_GQUAD_MAX_STACK_SIZE; i++)
    for(j = 3*VRNA_GQUAD_MIN_LINKER_LENGTH; j <= 3*VRNA_GQUAD_MAX_LINKER_LENGTH; j++){
      double GQuadAlpha_T = (double)GQuadAlphadH - (double)(GQuadAlphadH - GQuadAlpha37) * tempf;
      double GQuadBeta_T = (double)GQuadBetadH - (double)(GQuadBetadH - GQuadBeta37) * tempf;
      params->gquad[i][j] = (int)GQuadAlpha_T*(i-1) + (int)(((double)GQuadBeta_T)*log(j - 2));
    }

  for (i=0; i<31; i++)
    params->hairpin[i]  = hairpindH[i] - (hairpindH[i] - hairpin37[i])*tempf;
  for (i=0; i<=MIN2(30,MAXLOOP); i++) {
    params->bulge[i]          = bulgedH[i] - (bulgedH[i] - bulge37[i]) * tempf;
    params->internal_loop[i]  = internal_loopdH[i] - (internal_loopdH[i] - internal_loop37[i]) * tempf;
  }
  params->lxc = lxc37*tempf;
  for (; i<=MAXLOOP; i++) {
    params->bulge[i] = params->bulge[30]+(int)(params->lxc*log((double)(i)/30.));
    params->internal_loop[i] = params->internal_loop[30]+(int)(params->lxc*log((double)(i)/30.));
  }

  params->ninio[2] = niniodH - (niniodH - ninio37) * tempf;

  params->TripleC = TripleCdH - (TripleCdH - TripleC37) * tempf;
  params->MultipleCA = MultipleCAdH - (MultipleCAdH - MultipleCA37) * tempf;
  params->MultipleCB = MultipleCBdH - (MultipleCBdH - MultipleCB37) * tempf;

  for (i=0; (i*7)<strlen(Tetraloops); i++)
    params->Tetraloop_E[i] = TetraloopdH[i] - (TetraloopdH[i]-Tetraloop37[i])*tempf;
  for (i=0; (i*5)<strlen(Triloops); i++)
    params->Triloop_E[i] =  TriloopdH[i] - (TriloopdH[i]-Triloop37[i])*tempf;
  for (i=0; (i*9)<strlen(Hexaloops); i++)
    params->Hexaloop_E[i] =  HexaloopdH[i] - (HexaloopdH[i]-Hexaloop37[i])*tempf;

  params->TerminalAU = TerminalAUdH - (TerminalAUdH - TerminalAU37) * tempf;

  params->DuplexInit = DuplexInitdH - (DuplexInitdH - DuplexInit37) *tempf;

  params->MLbase = ML_BASEdH - (ML_BASEdH - ML_BASE37) * tempf;

  for (i=0; i<=NBPAIRS; i++)
    params->MLintern[i] = ML_interndH - (ML_interndH - ML_intern37) * tempf;

  params->MLclosing = ML_closingdH - (ML_closingdH - ML_closing37) * tempf;


  /* stacks    G(T) = H - [H - G(T0)]*T/T0 */
  for (i=0; i<=NBPAIRS; i++)
    for (j=0; j<=NBPAIRS; j++)
      params->stack[i][j] = stackdH[i][j] - (stackdH[i][j] - stack37[i][j])*tempf;

  /* mismatches */
  for (i=0; i<=NBPAIRS; i++)
    for (j=0; j<5; j++)
      for (k=0; k<5; k++) {
        int mm;
        params->mismatchI[i][j][k]    = mismatchIdH[i][j][k] - (mismatchIdH[i][j][k] - mismatchI37[i][j][k])*tempf;
        params->mismatchH[i][j][k]    = mismatchHdH[i][j][k] - (mismatchHdH[i][j][k] - mismatchH37[i][j][k])*tempf;
        params->mismatch1nI[i][j][k]  = mismatch1nIdH[i][j][k]-(mismatch1nIdH[i][j][k]-mismatch1nI37[i][j][k])*tempf;/* interior nx1 loops */
        params->mismatch23I[i][j][k]  = mismatch23IdH[i][j][k]-(mismatch23IdH[i][j][k]-mismatch23I37[i][j][k])*tempf;/* interior 2x3 loops */
        if(md->dangles){
          mm                      = mismatchMdH[i][j][k] - (mismatchMdH[i][j][k] - mismatchM37[i][j][k])*tempf;
          params->mismatchM[i][j][k]    = (mm > 0) ? 0 : mm;
          mm                      = mismatchExtdH[i][j][k] - (mismatchExtdH[i][j][k] - mismatchExt37[i][j][k])*tempf;
          params->mismatchExt[i][j][k]  = (mm > 0) ? 0 : mm;
        }
        else{
          params->mismatchM[i][j][k] = params->mismatchExt[i][j][k] = 0;
        }
      }

  /* dangles */
  for (i=0; i<=NBPAIRS; i++)
    for (j=0; j<5; j++) {
      int dd;
      dd = dangle5_dH[i][j] - (dangle5_dH[i][j] - dangle5_37[i][j])*tempf;
      params->dangle5[i][j] = (dd>0) ? 0 : dd;  /* must be <= 0 */
      dd = dangle3_dH[i][j] - (dangle3_dH[i][j] - dangle3_37[i][j])*tempf;
      params->dangle3[i][j] = (dd>0) ? 0 : dd;  /* must be <= 0 */
    }
  /* interior 1x1 loops */
  for (i=0; i<=NBPAIRS; i++)
    for (j=0; j<=NBPAIRS; j++)
      for (k=0; k<5; k++)
        for (l=0; l<5; l++)
          params->int11[i][j][k][l] = int11_dH[i][j][k][l] - (int11_dH[i][j][k][l] - int11_37[i][j][k][l])*tempf;

  /* interior 2x1 loops */
  for (i=0; i<=NBPAIRS; i++)
    for (j=0; j<=NBPAIRS; j++)
      for (k=0; k<5; k++)
        for (l=0; l<5; l++) {
          int m;
          for (m=0; m<5; m++)
            params->int21[i][j][k][l][m] = int21_dH[i][j][k][l][m] - (int21_dH[i][j][k][l][m] - int21_37[i][j][k][l][m])*tempf;
        }
  /* interior 2x2 loops */
  for (i=0; i<=NBPAIRS; i++)
    for (j=0; j<=NBPAIRS; j++)
      for (k=0; k<5; k++)
        for (l=0; l<5; l++) {
          int m,n;
          for (m=0; m<5; m++)
            for (n=0; n<5; n++)
              params->int22[i][j][k][l][m][n] = int22_dH[i][j][k][l][m][n] - (int22_dH[i][j][k][l][m][n]-int22_37[i][j][k][l][m][n])*tempf;
        }

  strncpy(params->Tetraloops, Tetraloops, 281);
  strncpy(params->Triloops, Triloops, 241);
  strncpy(params->Hexaloops, Hexaloops, 361);

  params->id = ++id;
  return params;
}

PRIVATE vrna_exp_param_t *
get_scaled_exp_params(vrna_md_t *md,
                      double pfs){

  unsigned  int i, j, k, l;
  double        kT, TT;
  double        GT;
  vrna_exp_param_t     *pf;

  pf                = (vrna_exp_param_t *)vrna_alloc(sizeof(vrna_exp_param_t));
  pf->model_details = *md;
  pf->temperature   = md->temperature;
  pf->alpha         = md->betaScale;
  pf->kT = kT       = md->betaScale * (md->temperature + K0) * GASCONST;   /* kT in cal/mol  */
  pf->pf_scale      = pfs;
  TT                = (md->temperature + K0) / (Tmeasure);

  for(i = VRNA_GQUAD_MIN_STACK_SIZE; i <= VRNA_GQUAD_MAX_STACK_SIZE; i++)
    for(j = 3*VRNA_GQUAD_MIN_LINKER_LENGTH; j <= 3*VRNA_GQUAD_MAX_LINKER_LENGTH; j++){
      double GQuadAlpha_T = (double)GQuadAlphadH - (double)(GQuadAlphadH - GQuadAlpha37) * TT;
      double GQuadBeta_T = (double)GQuadBetadH - (double)(GQuadBetadH - GQuadBeta37) * TT;
      GT = ((double)GQuadAlpha_T)*((double)(i-1)) + ((double)GQuadBeta_T)*log(((double)j) - 2.);
      pf->expgquad[i][j] = exp( -GT*10./kT);
    }

  /* loop energies: hairpins, bulges, interior, mulit-loops */
  for (i=0; i<31; i++){
    GT  = hairpindH[i] - (hairpindH[i] - hairpin37[i])*TT;
    pf->exphairpin[i] = exp( -GT*10./kT);
  }

  for (i=0; i<=MIN2(30, MAXLOOP); i++) {
    GT =  bulgedH[i]- (bulgedH[i] - bulge37[i])*TT;
    pf->expbulge[i] = exp( -GT*10./kT);
    GT =  internal_loopdH[i] - (internal_loopdH[i] - internal_loop37[i])*TT;
    pf->expinternal[i] = exp( -GT*10./kT);
  }
  /* special case of size 2 interior loops (single mismatch) */
  if (james_rule) pf->expinternal[2] = exp ( -80*10./kT);

  pf->lxc = lxc37*TT;

  GT =  DuplexInitdH - (DuplexInitdH - DuplexInit37)*TT;
  pf->expDuplexInit = exp( -GT*10./kT);

  for (i=31; i<=MAXLOOP; i++) {
    GT = bulge37[30]*TT + (pf->lxc*log( i/30.));
    pf->expbulge[i] = exp( -GT*10./kT);
    GT = internal_loop37[30]*TT + (pf->lxc*log( i/30.));
    pf->expinternal[i] = exp( -GT*10./kT);
  }

  GT = niniodH - (niniodH - ninio37)*TT;
  for (j=0; j<=MAXLOOP; j++)
      pf->expninio[2][j]=exp(-MIN2(MAX_NINIO,j*GT)*10./kT);

  for (i=0; (i*7)<strlen(Tetraloops); i++) {
    GT = TetraloopdH[i] - (TetraloopdH[i]-Tetraloop37[i])*TT;
    pf->exptetra[i] = exp( -GT*10./kT);
  }
  for (i=0; (i*5)<strlen(Triloops); i++) {
    GT = TriloopdH[i] - (TriloopdH[i]-Triloop37[i])*TT;
    pf->exptri[i] = exp( -GT*10./kT);
  }
  for (i=0; (i*9)<strlen(Hexaloops); i++) {
    GT = HexaloopdH[i] - (HexaloopdH[i]-Hexaloop37[i])*TT;
    pf->exphex[i] = exp( -GT*10./kT);
  }
  GT =  ML_closingdH - (ML_closingdH - ML_closing37)*TT;
  pf->expMLclosing = exp( -GT*10./kT);

  for (i=0; i<=NBPAIRS; i++) {
    GT =  ML_interndH - (ML_interndH - ML_intern37)*TT;
    /* if (i>2) GT += TerminalAU; */
    pf->expMLintern[i] = exp( -GT*10./kT);
  }
  GT = TerminalAUdH - (TerminalAUdH - TerminalAU37)*TT;
  pf->expTermAU = exp(-GT*10./kT);

  GT = ML_BASEdH - (ML_BASEdH - ML_BASE37)*TT;

  pf->expMLbase=exp(-10.*GT/kT);


  /* if dangles==0 just set their energy to 0,
     don't let dangle energies become > 0 (at large temps),
     but make sure go smoothly to 0                        */
  for (i=0; i<=NBPAIRS; i++)
    for (j=0; j<=4; j++) {
      if (md->dangles) {
        GT = dangle5_dH[i][j] - (dangle5_dH[i][j] - dangle5_37[i][j])*TT;
        pf->expdangle5[i][j] = exp(SMOOTH(-GT)*10./kT);
        GT = dangle3_dH[i][j] - (dangle3_dH[i][j] - dangle3_37[i][j])*TT;
        pf->expdangle3[i][j] =  exp(SMOOTH(-GT)*10./kT);
      } else
        pf->expdangle3[i][j] = pf->expdangle5[i][j] = 1;
    }

  /* stacking energies */
  for (i=0; i<=NBPAIRS; i++)
    for (j=0; j<=NBPAIRS; j++) {
      GT =  stackdH[i][j] - (stackdH[i][j] - stack37[i][j])*TT;
      pf->expstack[i][j] = exp( -GT*10./kT);
    }

  /* mismatch energies */
  for (i=0; i<=NBPAIRS; i++)
    for (j=0; j<5; j++)
      for (k=0; k<5; k++) {
        GT =  mismatchIdH[i][j][k] - ( mismatchIdH[i][j][k] - mismatchI37[i][j][k])*TT;
        pf->expmismatchI[i][j][k] = exp(-GT*10.0/kT);
        GT = mismatch1nIdH[i][j][k] - (mismatch1nIdH[i][j][k] - mismatch1nI37[i][j][k])*TT;
        pf->expmismatch1nI[i][j][k] = exp(-GT*10.0/kT);
        GT = mismatchHdH[i][j][k] - (mismatchHdH[i][j][k] - mismatchH37[i][j][k])*TT;
        pf->expmismatchH[i][j][k] = exp(-GT*10.0/kT);
        if (md->dangles) {
          GT = mismatchMdH[i][j][k] - (mismatchMdH[i][j][k] - mismatchM37[i][j][k])*TT;
          pf->expmismatchM[i][j][k] = exp(SMOOTH(-GT)*10.0/kT);
          GT = mismatchExtdH[i][j][k] - (mismatchExtdH[i][j][k] - mismatchExt37[i][j][k])*TT;
          pf->expmismatchExt[i][j][k] = exp(SMOOTH(-GT)*10.0/kT);
        }
        else{
          pf->expmismatchM[i][j][k] = pf->expmismatchExt[i][j][k] = 1.;
        }
        GT = mismatch23IdH[i][j][k] - (mismatch23IdH[i][j][k] - mismatch23I37[i][j][k])*TT;
        pf->expmismatch23I[i][j][k] = exp(-GT*10.0/kT);
      }

  /* interior lops of length 2 */
  for (i=0; i<=NBPAIRS; i++)
    for (j=0; j<=NBPAIRS; j++)
      for (k=0; k<5; k++)
        for (l=0; l<5; l++) {
          GT = int11_dH[i][j][k][l] -
            (int11_dH[i][j][k][l] - int11_37[i][j][k][l])*TT;
          pf->expint11[i][j][k][l] = exp(-GT*10./kT);
        }
  /* interior 2x1 loops */
  for (i=0; i<=NBPAIRS; i++)
    for (j=0; j<=NBPAIRS; j++)
      for (k=0; k<5; k++)
        for (l=0; l<5; l++) {
          int m;
          for (m=0; m<5; m++) {
            GT = int21_dH[i][j][k][l][m] -
              (int21_dH[i][j][k][l][m] - int21_37[i][j][k][l][m])*TT;
            pf->expint21[i][j][k][l][m] = exp(-GT*10./kT);
          }
        }

  /* interior 2x2 loops */
  for (i=0; i<=NBPAIRS; i++)
    for (j=0; j<=NBPAIRS; j++)
      for (k=0; k<5; k++)
        for (l=0; l<5; l++) {
          int m,n;
          for (m=0; m<5; m++)
            for (n=0; n<5; n++) {
              GT = int22_dH[i][j][k][l][m][n] -
                (int22_dH[i][j][k][l][m][n]-int22_37[i][j][k][l][m][n])*TT;
              pf->expint22[i][j][k][l][m][n] = exp(-GT*10./kT);
            }
        }

  strncpy(pf->Tetraloops, Tetraloops, 281);
  strncpy(pf->Triloops, Triloops, 241);
  strncpy(pf->Hexaloops, Hexaloops, 361);

  return pf;
}

PRIVATE vrna_exp_param_t *
get_exp_params_ali( vrna_md_t *md,
                    unsigned int n_seq,
                    double pfs){

  /* scale energy parameters and pre-calculate Boltzmann weights */
  unsigned int  i, j, k, l;
  double        kTn, TT;
  double        GT;
  vrna_exp_param_t     *pf;

  pf                = (vrna_exp_param_t *)vrna_alloc(sizeof(vrna_exp_param_t));
  pf->model_details = *md;
  pf->alpha         = md->betaScale;
  pf->temperature   = md->temperature;
  pf->pf_scale      = pfs;
  pf->kT = kTn      = ((double)n_seq)*md->betaScale*(md->temperature+K0)*GASCONST;   /* kT in cal/mol  */
  TT                = (md->temperature+K0)/(Tmeasure);


   /* loop energies: hairpins, bulges, interior, mulit-loops */
  for (i=0; i<31; i++) {
    GT =  hairpindH[i] - (hairpindH[i] - hairpin37[i])*TT;
    pf->exphairpin[i] = exp( -GT*10./kTn);
  }
  /*add penalty for too short hairpins*/
  for (i=0; i<3; i++) {
    GT= 600/*Penalty*/*TT;
    pf->exphairpin[i] = exp( -GT*10./kTn);
  }

  for (i=0; i<=MIN2(30, MAXLOOP); i++) {
    GT =  bulgedH[i]- (bulgedH[i] - bulge37[i])*TT;
    pf->expbulge[i] = exp( -GT*10./kTn);
    GT =  internal_loopdH[i] - (internal_loopdH[i] - internal_loop37[i])*TT;
    pf->expinternal[i] = exp( -GT*10./kTn);
  }
  /* special case of size 2 interior loops (single mismatch) */
  if (james_rule) pf->expinternal[2] = exp ( -80*10./kTn);

  pf->lxc = lxc37*TT;

  GT =  DuplexInitdH - (DuplexInitdH - DuplexInit37)*TT;
  pf->expDuplexInit = exp( -GT*10./kTn);

  for (i=31; i<=MAXLOOP; i++) {
    GT = bulge37[30]*TT + (pf->lxc*log( i/30.));
    pf->expbulge[i] = exp( -GT*10./kTn);
    GT = internal_loop37[30]*TT + (pf->lxc*log( i/30.));
    pf->expinternal[i] = exp( -GT*10./kTn);
  }

  GT = niniodH - (niniodH - ninio37)*TT;
  for (j=0; j<=MAXLOOP; j++)
    pf->expninio[2][j]=exp(-MIN2(MAX_NINIO,j*GT)*10./kTn);

  for (i=0; (i*7)<strlen(Tetraloops); i++) {
    GT = TetraloopdH[i] - (TetraloopdH[i]-Tetraloop37[i])*TT;
    pf->exptetra[i] = exp( -GT*10./kTn);
  }
  for (i=0; (i*5)<strlen(Triloops); i++) {
    GT = TriloopdH[i] - (TriloopdH[i]-Triloop37[i])*TT;
    pf->exptri[i] = exp( -GT*10./kTn);
  }
  for (i=0; (i*9)<strlen(Hexaloops); i++) {
    GT = HexaloopdH[i] - (HexaloopdH[i]-Hexaloop37[i])*TT;
    pf->exphex[i] = exp( -GT*10./kTn);
  }
  GT =  ML_closingdH - (ML_closingdH - ML_closing37)*TT;
  pf->expMLclosing = exp( -GT*10./kTn);

  for (i=0; i<=NBPAIRS; i++) { /* includes AU penalty */
    GT =  ML_interndH - (ML_interndH - ML_intern37)*TT;
    /* if (i>2) GT += TerminalAU; */
    pf->expMLintern[i] = exp( -GT*10./kTn);
  }
  GT = TerminalAUdH - (TerminalAUdH - TerminalAU37)*TT;
  pf->expTermAU = exp(-GT*10./kTn);

  GT = ML_BASEdH - (ML_BASEdH - ML_BASE37)*TT;
  pf->expMLbase=exp(-10.*GT/(kTn/n_seq));


  /* if dangle_model==0 just set their energy to 0,
     don't let dangle energies become > 0 (at large temps),
     but make sure go smoothly to 0                        */
  for (i=0; i<=NBPAIRS; i++)
    for (j=0; j<=4; j++) {
      if (md->dangles) {
        GT = dangle5_dH[i][j] - (dangle5_dH[i][j] - dangle5_37[i][j])*TT;
        pf->expdangle5[i][j] = exp(SMOOTH(-GT)*10./kTn);
        GT = dangle3_dH[i][j] - (dangle3_dH[i][j] - dangle3_37[i][j])*TT;
        pf->expdangle3[i][j] =  exp(SMOOTH(-GT)*10./kTn);
      } else
        pf->expdangle3[i][j] = pf->expdangle5[i][j] = 1;
    }

  /* stacking energies */
  for (i=0; i<=NBPAIRS; i++)
    for (j=0; j<=NBPAIRS; j++) {
      GT =  stackdH[i][j] - (stackdH[i][j] - stack37[i][j])*TT;
      pf->expstack[i][j] = exp( -GT*10./kTn);
    }

  /* mismatch energies */
  for (i=0; i<=NBPAIRS; i++)
    for (j=0; j<5; j++)
      for (k=0; k<5; k++) {
        GT =  mismatchIdH[i][j][k] - ( mismatchIdH[i][j][k] - mismatchI37[i][j][k])*TT;
        pf->expmismatchI[i][j][k] = exp(-GT*10.0/kTn);
        GT = mismatch1nIdH[i][j][k] - (mismatch1nIdH[i][j][k] - mismatch1nI37[i][j][k])*TT;
        pf->expmismatch1nI[i][j][k] = exp(-GT*10.0/kTn);
        GT = mismatchHdH[i][j][k] - (mismatchHdH[i][j][k] - mismatchH37[i][j][k])*TT;
        pf->expmismatchH[i][j][k] = exp(-GT*10.0/kTn);
        if (md->dangles) {
          GT = mismatchMdH[i][j][k] - (mismatchMdH[i][j][k] - mismatchM37[i][j][k])*TT;
          pf->expmismatchM[i][j][k] = exp(SMOOTH(-GT)*10.0/kTn);
          GT = mismatchExtdH[i][j][k] - (mismatchExtdH[i][j][k] - mismatchExt37[i][j][k])*TT;
          pf->expmismatchExt[i][j][k] = exp(SMOOTH(-GT)*10.0/kTn);
        }
        else{
          pf->expmismatchM[i][j][k] = pf->expmismatchExt[i][j][k] = 1.;
        }
        GT = mismatch23IdH[i][j][k] - (mismatch23IdH[i][j][k] - mismatch23I37[i][j][k])*TT;
        pf->expmismatch23I[i][j][k] = exp(-GT*10.0/kTn);
      }


  /* interior lops of length 2 */
  for (i=0; i<=NBPAIRS; i++)
    for (j=0; j<=NBPAIRS; j++)
      for (k=0; k<5; k++)
        for (l=0; l<5; l++) {
          GT = int11_dH[i][j][k][l] -
            (int11_dH[i][j][k][l] - int11_37[i][j][k][l])*TT;
          pf->expint11[i][j][k][l] = exp(-GT*10./kTn);
        }
  /* interior 2x1 loops */
  for (i=0; i<=NBPAIRS; i++)
    for (j=0; j<=NBPAIRS; j++)
      for (k=0; k<5; k++)
        for (l=0; l<5; l++) {
          int m;
          for (m=0; m<5; m++) {
            GT = int21_dH[i][j][k][l][m] -
              (int21_dH[i][j][k][l][m] - int21_37[i][j][k][l][m])*TT;
            pf->expint21[i][j][k][l][m] = exp(-GT*10./kTn);
          }
        }

  /* interior 2x2 loops */
  for (i=0; i<=NBPAIRS; i++)
    for (j=0; j<=NBPAIRS; j++)
      for (k=0; k<5; k++)
        for (l=0; l<5; l++) {
          int m,n;
          for (m=0; m<5; m++)
            for (n=0; n<5; n++) {
              GT = int22_dH[i][j][k][l][m][n] -
                (int22_dH[i][j][k][l][m][n]-int22_37[i][j][k][l][m][n])*TT;
              pf->expint22[i][j][k][l][m][n] = exp(-GT*10./kTn);
            }
        }

  strncpy(pf->Tetraloops, Tetraloops, 281);
  strncpy(pf->Triloops, Triloops, 241);
  strncpy(pf->Hexaloops, Hexaloops, 361);

  return pf;
}

PRIVATE void
rescale_params(vrna_fold_compound_t *vc){

  int           i;
  vrna_exp_param_t  *pf = vc->exp_params;
  vrna_mx_pf_t      *m  = vc->exp_matrices;

  if(m && pf){
    m->scale[0] = 1.;
    m->scale[1] = (FLT_OR_DBL)(1./pf->pf_scale);
    m->expMLbase[0] = 1;
    m->expMLbase[1] = (FLT_OR_DBL)(pf->expMLbase / pf->pf_scale);
    for (i=2; i<=vc->length; i++) {
      m->scale[i] = m->scale[i/2]*m->scale[i-(i/2)];
      m->expMLbase[i] = (FLT_OR_DBL)pow(pf->expMLbase, (double)i) * m->scale[i];
    }
  }
}


#ifdef  VRNA_BACKWARD_COMPAT

/*###########################################*/
/*# deprecated functions below              #*/
/*###########################################*/

PUBLIC vrna_param_t *
scale_parameters(void){

  vrna_md_t md;
  set_model_details(&md);
  return vrna_params(&md);
}

PUBLIC vrna_param_t *
get_scaled_parameters(double temp,
                      vrna_md_t md){

  md.temperature = temp;
  return get_scaled_params(&md);
}

PUBLIC vrna_exp_param_t *
get_boltzmann_factors(double temp,
                      double betaScale,
                      vrna_md_t md,
                      double pfs){

  md.temperature  = temp;
  md.betaScale    = betaScale;
  pf_scale        = pfs;

  return get_scaled_exp_params(&md, pfs);
}

PUBLIC vrna_exp_param_t *
get_scaled_pf_parameters(void){

  vrna_md_t         md;
  vrna_exp_param_t  *pf;

  set_model_details(&md);

  pf = vrna_exp_params(&md);
  pf->pf_scale = pf_scale;

  return pf;
}

PUBLIC vrna_exp_param_t *
get_boltzmann_factors_ali(unsigned int n_seq,
                         double temp,
                         double betaScale,
                         vrna_md_t md,
                         double pfs){

  md.temperature  = temp;
  md.betaScale    = betaScale;
  pf_scale        = pfs;

  return get_exp_params_ali(&md, n_seq, pfs);
  
}

PUBLIC vrna_exp_param_t *
get_scaled_alipf_parameters(unsigned int n_seq){

  vrna_md_t  md;

  set_model_details(&md);

  return get_exp_params_ali(&md, n_seq, pf_scale);
}

PUBLIC vrna_exp_param_t *
get_boltzmann_factor_copy(vrna_exp_param_t *par){

  return vrna_exp_params_copy(par);
}

PUBLIC vrna_param_t *get_parameter_copy(vrna_param_t *par){

  return vrna_params_copy(par);
}

PUBLIC vrna_param_t *copy_parameters(void){
  vrna_param_t *copy;
  if (p.id != id){
    vrna_md_t md;
    set_model_details(&md);
    return vrna_params(&md);
  } else {
    copy = (vrna_param_t *) vrna_alloc(sizeof(vrna_param_t));
    memcpy(copy, &p, sizeof(vrna_param_t));
  }
  return copy;
}

PUBLIC vrna_param_t *set_parameters(vrna_param_t *dest){
  memcpy(&p, dest, sizeof(vrna_param_t));
  return &p;
}

PUBLIC vrna_exp_param_t *copy_pf_param(void){
  vrna_exp_param_t *copy;
  if (pf.id != pf_id){
    vrna_md_t md;
    set_model_details(&md);
    copy = vrna_exp_params(&md);
    copy->pf_scale = pf_scale;
    return copy;
  }
  else{
    copy = (vrna_exp_param_t *) vrna_alloc(sizeof(vrna_exp_param_t));
    memcpy(copy, &pf, sizeof(vrna_exp_param_t));
  }
  return copy;
}

PUBLIC vrna_exp_param_t *set_pf_param(vrna_param_t *dest){
  memcpy(&pf, dest, sizeof(vrna_exp_param_t));
  return &pf;
}

PUBLIC vrna_exp_param_t *scale_pf_parameters(void){
  vrna_md_t         md;
  vrna_exp_param_t  *pf;

  set_model_details(&md);

  pf = vrna_exp_params(&md);
  pf->pf_scale = pf_scale;

  return pf;
}

#endif

