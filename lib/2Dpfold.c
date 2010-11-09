/*
      minimum free energy
      RNA secondary structure with
      basepair distance d_1 to reference structure 1 and distance d_2 to reference structure 2

*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <string.h>
#include <float.h>    /* #defines FLT_MAX ... */
#include "utils.h"
#include "fold_vars.h"
#include "params.h"
#include "energy_par.h"
#include "loop_energies.h"
#include "pair_mat.h"
#include "mm.h"
#include "2Dpfold.h"

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
PRIVATE         void  pf2D_linear(TwoDpfold_vars *vars);
PRIVATE         void  pf2D_circ(TwoDpfold_vars *vars);
PRIVATE         void  initialize_TwoDpfold_vars(TwoDpfold_vars *vars);
PRIVATE         void  scale_pf_params2(TwoDpfold_vars *vars);
PRIVATE         void  make_ptypes2(TwoDpfold_vars *vars);
PRIVATE         void  backtrack(TwoDpfold_vars *vars, char *pstruc, unsigned int d1, unsigned int d2, unsigned int i, unsigned int j);
PRIVATE         void  backtrack_qm1(TwoDpfold_vars *vars, char *pstruc, unsigned int d1, unsigned int d2, unsigned int i, unsigned int j);
PRIVATE         void  backtrack_qm(TwoDpfold_vars *vars, char *pstruc, unsigned int d1, unsigned int d2, unsigned int i, unsigned int j);
PRIVATE         void  adjustArrayBoundaries(FLT_OR_DBL ***array, int *k_min, int *k_max, int **l_min, int **l_max, int k_min_real, int k_max_real, int *l_min_real, int *l_max_real);
INLINE  PRIVATE void  preparePosteriorBoundaries(int size, int shift, int *min_k, int *max_k, int **min_l, int **max_l);
INLINE  PRIVATE void  updatePosteriorBoundaries(int d1, int d2, int *min_k, int *max_k, int **min_l, int **max_l);
INLINE  PRIVATE void  prepareBoundaries(int min_k_pre, int max_k_pre, int min_l_pre, int max_l_pre, int bpdist, int *min_k, int *max_k, int **min_l, int **max_l);
INLINE  PRIVATE void  prepareArray(FLT_OR_DBL ***array, int min_k, int max_k, int *min_l, int *max_l);

/*
#################################
# BEGIN OF FUNCTION DEFINITIONS #
#################################
*/

PUBLIC TwoDpfold_vars *get_TwoDpfold_variables(const char *seq, const char *structure1, char *structure2, int circ){
  unsigned int i, size, length;
  TwoDpfold_vars *vars;
  vars = (TwoDpfold_vars *)malloc(sizeof(TwoDpfold_vars));
  length              = (unsigned int) strlen(seq);
  vars->sequence      = (char *)space(length + 1);
  strcpy(vars->sequence, seq);
  vars->seq_length    = length;
  if(vars->seq_length < 1) nrerror("get_TwoDfold_variables: sequence must be longer than 0");
  size                        = ((length + 1) * (length + 2)/2);

  vars->reference_pt1 = make_pair_table(structure1);
  vars->reference_pt2 = make_pair_table(structure2);
  vars->referenceBPs1 = make_referenceBP_array(vars->reference_pt1, TURN); /* matrix containing number of basepairs of reference structure1 in interval [i,j] */
  vars->referenceBPs2 = make_referenceBP_array(vars->reference_pt2, TURN); /* matrix containing number of basepairs of reference structure2 in interval [i,j] */

  /* compute maximum matching with reference structure 1 disallowed */
  vars->mm1           = maximumMatchingConstraint(vars->sequence, vars->reference_pt1);
  /* compute maximum matching with reference structure 2 disallowed */
  vars->mm2           = maximumMatchingConstraint(vars->sequence, vars->reference_pt2);
  vars->bpdist        = compute_BPdifferences(vars->reference_pt1, vars->reference_pt2, TURN);

  vars->dangles       = dangles;
  vars->circ          = circ;
  vars->temperature   = temperature;
  vars->init_temp     = temperature;
  vars->pf_scale      = pf_scale;

  vars->scale         = (FLT_OR_DBL *) space(sizeof(FLT_OR_DBL)*(length+1));
  vars->ptype         = space(sizeof(char) * size);
  vars->S             = NULL;
  vars->S1            = NULL;
  vars->maxD1         = 0;
  vars->maxD2         = 0;

  vars->jindx      = get_indx(length);
  vars->my_iindx   = get_iindx(length);

  /* allocate memory for the pf matrices and min-/max-index helper arrays */
  vars->Q                = (FLT_OR_DBL ***) space(sizeof(FLT_OR_DBL **)  * size);
  vars->l_min_values     = (int **)  space(sizeof(int *)   * size);
  vars->l_max_values     = (int **)  space(sizeof(int *)   * size);
  vars->k_min_values     = (int *)   space(sizeof(int)     * size);
  vars->k_max_values     = (int *)   space(sizeof(int)     * size);

  vars->Q_B              = (FLT_OR_DBL ***) space(sizeof(FLT_OR_DBL **)  * size);
  vars->l_min_values_b   = (int **)  space(sizeof(int *)   * size);
  vars->l_max_values_b   = (int **)  space(sizeof(int *)   * size);
  vars->k_min_values_b   = (int *)   space(sizeof(int)     * size);
  vars->k_max_values_b   = (int *)   space(sizeof(int)     * size);

  vars->Q_M              = (FLT_OR_DBL ***) space(sizeof(FLT_OR_DBL **)  * size);
  vars->l_min_values_m   = (int **)  space(sizeof(int *)   * size);
  vars->l_max_values_m   = (int **)  space(sizeof(int *)   * size);
  vars->k_min_values_m   = (int *)   space(sizeof(int)     * size);
  vars->k_max_values_m   = (int *)   space(sizeof(int)     * size);

  vars->Q_M1             = (FLT_OR_DBL ***) space(sizeof(FLT_OR_DBL **)  * size);
  vars->l_min_values_m1  = (int **)  space(sizeof(int *)   * size);
  vars->l_max_values_m1  = (int **)  space(sizeof(int *)   * size);
  vars->k_min_values_m1  = (int *)   space(sizeof(int)     * size);
  vars->k_max_values_m1  = (int *)   space(sizeof(int)     * size);


  if(vars->circ){
    vars->Q_M2             = (FLT_OR_DBL ***) space(sizeof(FLT_OR_DBL **)  * (length + 1));
    vars->l_min_values_m2  = (int **)  space(sizeof(int *)   * (length + 1));
    vars->l_max_values_m2  = (int **)  space(sizeof(int *)   * (length + 1));
    vars->k_min_values_m2  = (int *)   space(sizeof(int)     * (length + 1));
    vars->k_max_values_m2  = (int *)   space(sizeof(int)     * (length + 1));
  }
  else{
    vars->Q_M2             = NULL;
    vars->l_min_values_m2  = NULL;
    vars->l_max_values_m2  = NULL;
    vars->k_min_values_m2  = NULL;
    vars->k_max_values_m2  = NULL;
  }
  vars->Q_c               = NULL;
  vars->Q_cH              = NULL;
  vars->Q_cI              = NULL;
  vars->Q_cM              = NULL;

  return vars;
}

PUBLIC TwoDpfold_vars *get_TwoDpfold_variables_from_MFE(TwoDfold_vars *mfe_vars){
  unsigned int i, j, size, length;
  int cnt1, cnt2, cnt3, cnt4, ij, jij;
  TwoDpfold_vars *vars;
  vars                = (TwoDpfold_vars *)malloc(sizeof(TwoDpfold_vars));
  vars->sequence      = strdup(mfe_vars->sequence);
  length              = mfe_vars->seq_length;
  vars->seq_length    = length;
  size                = ((length + 1) * (length + 2)/2);

  vars->reference_pt1 = copy_pair_table(mfe_vars->reference_pt1);
  vars->reference_pt2 = copy_pair_table(mfe_vars->reference_pt2);
  vars->referenceBPs1 = make_referenceBP_array(vars->reference_pt1, TURN);
  vars->referenceBPs2 = make_referenceBP_array(vars->reference_pt2, TURN);

  /* compute maximum matching with reference structure 1 disallowed */
  vars->mm1           = maximumMatchingConstraint(vars->sequence, vars->reference_pt1);
  /* compute maximum matching with reference structure 2 disallowed */
  vars->mm2           = maximumMatchingConstraint(vars->sequence, vars->reference_pt2);
  vars->bpdist        = compute_BPdifferences(vars->reference_pt1, vars->reference_pt2, TURN);


  vars->dangles      = mfe_vars->dangles;
  vars->circ         = mfe_vars->circ;
  vars->temperature  = mfe_vars->temperature;
  vars->init_temp    = mfe_vars->temperature;
  vars->pf_scale     = pf_scale;

  vars->scale        = (FLT_OR_DBL *) space(sizeof(FLT_OR_DBL)*(length+1));
  vars->ptype        = space(sizeof(char) * size);
  vars->S            = NULL;
  vars->S1           = NULL;
  vars->maxD1        = mfe_vars->maxD1;
  vars->maxD2        = mfe_vars->maxD2;

  vars->jindx      = get_indx(vars->seq_length);
  vars->my_iindx   = get_iindx(vars->seq_length);

  /* allocate memory for the pf matrices and min-/max-index helper arrays */
  vars->Q                = (FLT_OR_DBL ***) space(sizeof(FLT_OR_DBL **)  * size);
  vars->l_min_values     = (int **)  space(sizeof(int *)   * size);
  vars->l_max_values     = (int **)  space(sizeof(int *)   * size);
  vars->k_min_values     = (int *)   space(sizeof(int)     * size);
  vars->k_max_values     = (int *)   space(sizeof(int)     * size);

  vars->Q_B              = (FLT_OR_DBL ***) space(sizeof(FLT_OR_DBL **)  * size);
  vars->l_min_values_b   = (int **)  space(sizeof(int *)   * size);
  vars->l_max_values_b   = (int **)  space(sizeof(int *)   * size);
  vars->k_min_values_b   = (int *)   space(sizeof(int)     * size);
  vars->k_max_values_b   = (int *)   space(sizeof(int)     * size);

  vars->Q_M              = (FLT_OR_DBL ***) space(sizeof(FLT_OR_DBL **)  * size);
  vars->l_min_values_m   = (int **)  space(sizeof(int *)   * size);
  vars->l_max_values_m   = (int **)  space(sizeof(int *)   * size);
  vars->k_min_values_m   = (int *)   space(sizeof(int)     * size);
  vars->k_max_values_m   = (int *)   space(sizeof(int)     * size);

  vars->Q_M1             = (FLT_OR_DBL ***) space(sizeof(FLT_OR_DBL **)  * size);
  vars->l_min_values_m1  = (int **)  space(sizeof(int *)   * size);
  vars->l_max_values_m1  = (int **)  space(sizeof(int *)   * size);
  vars->k_min_values_m1  = (int *)   space(sizeof(int)     * size);
  vars->k_max_values_m1  = (int *)   space(sizeof(int)     * size);


  if(vars->circ){
    vars->Q_M2             = (FLT_OR_DBL ***) space(sizeof(FLT_OR_DBL **)  * (length + 1));
    vars->l_min_values_m2  = (int **)  space(sizeof(int *)   * (length + 1));
    vars->l_max_values_m2  = (int **)  space(sizeof(int *)   * (length + 1));
    vars->k_min_values_m2  = (int *)   space(sizeof(int)     * (length + 1));
    vars->k_max_values_m2  = (int *)   space(sizeof(int)     * (length + 1));
  }
  else{
    vars->Q_M2             = NULL;
    vars->l_min_values_m2  = NULL;
    vars->l_max_values_m2  = NULL;
    vars->k_min_values_m2  = NULL;
    vars->k_max_values_m2  = NULL;
  }
  vars->Q_c               = NULL;
  vars->Q_cH              = NULL;
  vars->Q_cI              = NULL;
  vars->Q_cM              = NULL;

  /* prepare everything that is known from previous mfe folding step */
  for(i = 1; i < length; i++)
    for(j = i; j <= length; j++){
      int mem;
      ij                        = vars->my_iindx[i] - j;
      jij                       = vars->jindx[j] + i;

      if(mfe_vars->E_C[ij]){
        mem                       = mfe_vars->k_max_values[ij] - mfe_vars->k_min_values[ij] + 1;
        vars->k_min_values_b[ij]  = mfe_vars->k_min_values[ij];
        vars->k_max_values_b[ij]  = mfe_vars->k_max_values[ij];
        vars->Q_B[ij]             = (FLT_OR_DBL **)space(sizeof(FLT_OR_DBL *) * mem);
        vars->l_min_values_b[ij]  = (int *)space(sizeof(int) * mem);
        vars->l_max_values_b[ij]  = (int *)space(sizeof(int) * mem);
        vars->Q_B[ij]             -= vars->k_min_values_b[ij];
        vars->l_min_values_b[ij]  -= vars->k_min_values_b[ij];
        vars->l_max_values_b[ij]  -= vars->k_min_values_b[ij];
        for(cnt1 = vars->k_min_values_b[ij]; cnt1 <= vars->k_max_values_b[ij]; cnt1++){
          vars->l_min_values_b[ij][cnt1] = mfe_vars->l_min_values[ij][cnt1];
          vars->l_max_values_b[ij][cnt1] = mfe_vars->l_max_values[ij][cnt1];
          if(vars->l_min_values_b[ij][cnt1] < INF){
            vars->Q_B[ij][cnt1] = (FLT_OR_DBL *)space(sizeof(FLT_OR_DBL) * ((vars->l_max_values_b[ij][cnt1] - vars->l_min_values_b[ij][cnt1] + 1)/2 + 1));
            vars->Q_B[ij][cnt1] -= vars->l_min_values_b[ij][cnt1]/2;
          }
        }
      }

      if(mfe_vars->E_M[ij]){
        mem                       = mfe_vars->k_max_values_m[ij] - mfe_vars->k_min_values_m[ij] + 1;
        vars->k_min_values_m[ij]  = mfe_vars->k_min_values_m[ij];
        vars->k_max_values_m[ij]  = mfe_vars->k_max_values_m[ij];
        vars->Q_M[ij]             = (FLT_OR_DBL **)space(sizeof(FLT_OR_DBL *) * mem);
        vars->l_min_values_m[ij]  = (int *)space(sizeof(int) * mem);
        vars->l_max_values_m[ij]  = (int *)space(sizeof(int) * mem);
        vars->Q_M[ij]             -= vars->k_min_values_m[ij];
        vars->l_min_values_m[ij]  -= vars->k_min_values_m[ij];
        vars->l_max_values_m[ij]  -= vars->k_min_values_m[ij];
        for(cnt1 = vars->k_min_values_m[ij]; cnt1 <= vars->k_max_values_m[ij]; cnt1++){
          vars->l_min_values_m[ij][cnt1] = mfe_vars->l_min_values_m[ij][cnt1];
          vars->l_max_values_m[ij][cnt1] = mfe_vars->l_max_values_m[ij][cnt1];
          if(vars->l_min_values_m[ij][cnt1] < INF){
            vars->Q_M[ij][cnt1] = (FLT_OR_DBL *)space(sizeof(FLT_OR_DBL) * ((vars->l_max_values_m[ij][cnt1] - vars->l_min_values_m[ij][cnt1] + 1)/2 + 1));
            vars->Q_M[ij][cnt1] -= vars->l_min_values_m[ij][cnt1]/2;
          }
        }
      }

      if(mfe_vars->E_M1[ij]){
        mem                       = mfe_vars->k_max_values_m1[ij] - mfe_vars->k_min_values_m1[ij] + 1;
        vars->k_min_values_m1[jij] = mfe_vars->k_min_values_m1[ij];
        vars->k_max_values_m1[jij] = mfe_vars->k_max_values_m1[ij];
        vars->Q_M1[jij]            = (FLT_OR_DBL **)space(sizeof(FLT_OR_DBL *) * mem);
        vars->l_min_values_m1[jij] = (int *)space(sizeof(int) * mem);
        vars->l_max_values_m1[jij] = (int *)space(sizeof(int) * mem);
        vars->Q_M1[jij]            -= vars->k_min_values_m1[jij];
        vars->l_min_values_m1[jij] -= vars->k_min_values_m1[jij];
        vars->l_max_values_m1[jij] -= vars->k_min_values_m1[jij];
        for(cnt1 = vars->k_min_values_m1[jij]; cnt1 <= vars->k_max_values_m1[jij]; cnt1++){
          vars->l_min_values_m1[jij][cnt1] = mfe_vars->l_min_values_m1[ij][cnt1];
          vars->l_max_values_m1[jij][cnt1] = mfe_vars->l_max_values_m1[ij][cnt1];
          if(vars->l_min_values_m1[jij][cnt1] < INF){
            vars->Q_M1[jij][cnt1]            = (FLT_OR_DBL *)space(sizeof(FLT_OR_DBL) * ((vars->l_max_values_m1[jij][cnt1] - vars->l_min_values_m1[jij][cnt1] + 1)/2 + 1));
            vars->Q_M1[jij][cnt1]            -= vars->l_min_values_m1[jij][cnt1]/2;
          }
        }
      }
    }

  if(vars->circ){
    int mem;
    if(mfe_vars->E_Fc){
      mem = mfe_vars->k_max_values_fc - mfe_vars->k_min_values_fc + 1;
      vars->Q_c             = (FLT_OR_DBL **)space(sizeof(FLT_OR_DBL *) * mem);
      vars->l_min_values_qc = (int *)space(sizeof(int) * mem);
      vars->l_max_values_qc = (int *)space(sizeof(int) * mem);
      vars->k_min_values_qc = mfe_vars->k_min_values_fc;
      vars->k_max_values_qc = mfe_vars->k_max_values_fc;
      vars->Q_c             -= vars->k_min_values_qc;
      vars->l_min_values_qc -= vars->k_min_values_qc;
      vars->l_max_values_qc -= vars->k_min_values_qc;
      for(cnt1 = vars->k_min_values_qc; cnt1 <= vars->k_max_values_qc; cnt1++){
        mem = (mfe_vars->l_max_values_fc[cnt1] - mfe_vars->l_min_values_fc[cnt1] + 1)/2 + 1;
        vars->l_min_values_qc[cnt1] = mfe_vars->l_min_values_fc[cnt1];
        vars->l_max_values_qc[cnt1] = mfe_vars->l_max_values_fc[cnt1];
        if(vars->l_min_values_qc[cnt1] < INF){
          vars->Q_c[cnt1]             = (FLT_OR_DBL *)space(sizeof(FLT_OR_DBL) * mem);
          vars->Q_c[cnt1]             -= vars->l_min_values_qc[cnt1]/2;
        }
      }
    }

    if(mfe_vars->E_FcH){
      mem = mfe_vars->k_max_values_fcH - mfe_vars->k_min_values_fcH + 1;
      vars->Q_cH              = (FLT_OR_DBL **)space(sizeof(FLT_OR_DBL *) * mem);
      vars->l_min_values_qcH  = (int *)space(sizeof(int) * mem);
      vars->l_max_values_qcH  = (int *)space(sizeof(int) * mem);
      vars->k_min_values_qcH  = mfe_vars->k_min_values_fcH;
      vars->k_max_values_qcH  = mfe_vars->k_max_values_fcH;
      vars->Q_cH              -= vars->k_min_values_qcH;
      vars->l_min_values_qcH  -= vars->k_min_values_qcH;
      vars->l_max_values_qcH  -= vars->k_min_values_qcH;
      for(cnt1 = vars->k_min_values_qcH; cnt1 <= vars->k_max_values_qcH; cnt1++){
        mem = (mfe_vars->l_max_values_fcH[cnt1] - mfe_vars->l_min_values_fcH[cnt1] + 1)/2 + 1;
        vars->l_min_values_qcH[cnt1]  = mfe_vars->l_min_values_fcH[cnt1];
        vars->l_max_values_qcH[cnt1]  = mfe_vars->l_max_values_fcH[cnt1];
        if(vars->l_min_values_qcH[cnt1] < INF){
          vars->Q_cH[cnt1]              = (FLT_OR_DBL *)space(sizeof(FLT_OR_DBL) * mem);
          vars->Q_cH[cnt1]              -= vars->l_min_values_qcH[cnt1]/2;
        }
      }
    }

    if(mfe_vars->E_FcI){
      mem = mfe_vars->k_max_values_fcI - mfe_vars->k_min_values_fcI + 1;
      vars->Q_cI              = (FLT_OR_DBL **)space(sizeof(FLT_OR_DBL *) * mem);
      vars->l_min_values_qcI  = (int *)space(sizeof(int) * mem);
      vars->l_max_values_qcI  = (int *)space(sizeof(int) * mem);
      vars->k_min_values_qcI  = mfe_vars->k_min_values_fcI;
      vars->k_max_values_qcI  = mfe_vars->k_max_values_fcI;
      vars->Q_cI              -= vars->k_min_values_qcI;
      vars->l_min_values_qcI  -= vars->k_min_values_qcI;
      vars->l_max_values_qcI  -= vars->k_min_values_qcI;
      for(cnt1 = vars->k_min_values_qcI; cnt1 <= vars->k_max_values_qcI; cnt1++){
        mem = (mfe_vars->l_max_values_fcI[cnt1] - mfe_vars->l_min_values_fcI[cnt1] + 1)/2 + 1;
        vars->l_min_values_qcI[cnt1]  = mfe_vars->l_min_values_fcI[cnt1];
        vars->l_max_values_qcI[cnt1]  = mfe_vars->l_max_values_fcI[cnt1];
        if(vars->l_min_values_qcI[cnt1] < INF){
          vars->Q_cI[cnt1]              = (FLT_OR_DBL *)space(sizeof(FLT_OR_DBL) * mem);
          vars->Q_cI[cnt1]              -= vars->l_min_values_qcI[cnt1]/2;
        }
      }
    }

    if(mfe_vars->E_FcM){
      mem = mfe_vars->k_max_values_fcM - mfe_vars->k_min_values_fcM + 1;
      vars->Q_cM              = (FLT_OR_DBL **)space(sizeof(FLT_OR_DBL *) * mem);
      vars->l_min_values_qcM  = (int *)space(sizeof(int) * mem);
      vars->l_max_values_qcM  = (int *)space(sizeof(int) * mem);
      vars->k_min_values_qcM  = mfe_vars->k_min_values_fcM;
      vars->k_max_values_qcM  = mfe_vars->k_max_values_fcM;
      vars->Q_cM              -= vars->k_min_values_qcM;
      vars->l_min_values_qcM  -= vars->k_min_values_qcM;
      vars->l_max_values_qcM  -= vars->k_min_values_qcM;
      for(cnt1 = vars->k_min_values_qcM; cnt1 <= vars->k_max_values_qcM; cnt1++){
        mem = (mfe_vars->l_max_values_fcM[cnt1] - mfe_vars->l_min_values_fcM[cnt1] + 1)/2 + 1;
        vars->l_min_values_qcM[cnt1]  = mfe_vars->l_min_values_fcM[cnt1];
        vars->l_max_values_qcM[cnt1]  = mfe_vars->l_max_values_fcM[cnt1];
        if(vars->l_min_values_qcM[cnt1] < INF){
          vars->Q_cM[cnt1]              = (FLT_OR_DBL *)space(sizeof(FLT_OR_DBL) * mem);
          vars->Q_cM[cnt1]              -= vars->l_min_values_qcM[cnt1]/2;
        }
      }
    }

    if(mfe_vars->E_M2){
      vars->Q_M2 = (FLT_OR_DBL ***)space(sizeof(FLT_OR_DBL **) * (length + 1));
      vars->k_min_values_m2 = (int *)space(sizeof(int) * (length + 1));
      vars->k_max_values_m2 = (int *)space(sizeof(int) * (length + 1));
      vars->l_min_values_m2 = (int **)space(sizeof(int *) * (length + 1));
      vars->l_max_values_m2 = (int **)space(sizeof(int *) * (length + 1));
      for(i=1; i < length - TURN - 1; i++){
        if(!mfe_vars->E_M2[i]) continue;
        mem = mfe_vars->k_max_values_m2[i] - mfe_vars->k_min_values_m2[i] + 1;
        vars->k_min_values_m2[i]  = mfe_vars->k_min_values_m2[i];
        vars->k_max_values_m2[i]  = mfe_vars->k_max_values_m2[i];
        vars->Q_M2[i]             = (FLT_OR_DBL **)space(sizeof(FLT_OR_DBL *) * mem);
        vars->l_min_values_m2[i]  = (int *)space(sizeof(int) * mem);
        vars->l_max_values_m2[i]  = (int *)space(sizeof(int) * mem);
        vars->Q_M2[i]             -= vars->k_min_values_m2[i];
        vars->l_min_values_m2[i]  -= vars->k_min_values_m2[i];
        vars->l_max_values_m2[i]  -= vars->k_min_values_m2[i];
        for(cnt1 = vars->k_min_values_m2[i]; cnt1 <= vars->k_max_values_m2[i]; cnt1++){
          mem = (mfe_vars->l_max_values_m2[i][cnt1] - mfe_vars->l_min_values_m2[i][cnt1] + 1)/2 + 1;
          vars->l_min_values_m2[i][cnt1]  = mfe_vars->l_min_values_m2[i][cnt1];
          vars->l_max_values_m2[i][cnt1]  = mfe_vars->l_max_values_m2[i][cnt1];
          if(vars->l_min_values_m2[i][cnt1] < INF){
            vars->Q_M2[i][cnt1]                  = (FLT_OR_DBL *)space(sizeof(FLT_OR_DBL) * mem);
            vars->Q_M2[i][cnt1]                  -= vars->l_min_values_m2[i][cnt1]/2;
          }
        }
      }
    }
  }
  else{
    vars->Q_M2 = NULL;
    vars->Q_c = vars->Q_cH = vars->Q_cI = vars->Q_cM = NULL;
    vars->k_min_values_m2   = vars->k_max_values_m2   = NULL;
    vars->l_min_values_m2   = vars->l_max_values_m2   = NULL;
    vars->l_min_values_qc   = vars->l_max_values_qc   = NULL;
    vars->l_min_values_qcH  = vars->l_max_values_qcH  = NULL;
    vars->l_min_values_qcI  = vars->l_max_values_qcI  = NULL;
    vars->l_min_values_qcM  = vars->l_max_values_qcM  = NULL;
  }

  return vars;
}

PUBLIC void destroy_TwoDpfold_variables(TwoDpfold_vars *vars){
  unsigned int size = ((vars->seq_length+1)*(vars->seq_length+2)/2);
  unsigned int i, j, ij;
  int cnt1, cnt2, cnt3, cnt4;
  if(vars == NULL) return;

#ifdef _OPENMP
  #pragma omp sections private(i,j,ij,cnt1,cnt2,cnt3,cnt4)
  {

  #pragma omp section
  {
#endif
  if(vars->Q != NULL){
    for(i = 1; i <= vars->seq_length; i++){
      for(j = i; j <= vars->seq_length; j++){
        ij = vars->my_iindx[i] - j;
        if(!vars->Q[ij]) continue;
        for(cnt1 = vars->k_min_values[ij]; cnt1 <= vars->k_max_values[ij]; cnt1++)
          if(vars->l_min_values[ij][cnt1] < INF){
            vars->Q[ij][cnt1] += vars->l_min_values[ij][cnt1]/2;
            free(vars->Q[ij][cnt1]);
          }
        if(vars->k_min_values[ij] < INF){
          vars->Q[ij] += vars->k_min_values[ij];
          free(vars->Q[ij]);
          vars->l_min_values[ij] += vars->k_min_values[ij];
          vars->l_max_values[ij] += vars->k_min_values[ij];
          free(vars->l_min_values[ij]);
          free(vars->l_max_values[ij]);
        }
      }
    }
    free(vars->Q);
    free(vars->l_min_values);
    free(vars->l_max_values);
    free(vars->k_min_values);
    free(vars->k_max_values);
  }
#ifdef _OPENMP
  }
  #pragma omp section
  {
#endif
  if(vars->Q_B != NULL){
    for(i = 1; i < vars->seq_length; i++){
      for(j = i; j <= vars->seq_length; j++){
        ij = vars->my_iindx[i] - j;
        if(!vars->Q_B[ij]) continue;
        for(cnt1 = vars->k_min_values_b[ij]; cnt1 <= vars->k_max_values_b[ij]; cnt1++)
          if(vars->l_min_values_b[ij][cnt1] < INF){
            vars->Q_B[ij][cnt1] += vars->l_min_values_b[ij][cnt1]/2;
            free(vars->Q_B[ij][cnt1]);
          }
        if(vars->k_min_values_b[ij] < INF){
          vars->Q_B[ij] += vars->k_min_values_b[ij];
          free(vars->Q_B[ij]);
          vars->l_min_values_b[ij] += vars->k_min_values_b[ij];
          vars->l_max_values_b[ij] += vars->k_min_values_b[ij];
          free(vars->l_min_values_b[ij]);
          free(vars->l_max_values_b[ij]);
        }
      }
    }
    free(vars->Q_B);
    free(vars->l_min_values_b);
    free(vars->l_max_values_b);
    free(vars->k_min_values_b);
    free(vars->k_max_values_b);
  }
#ifdef _OPENMP
  }
  #pragma omp section
  {
#endif
  if(vars->Q_M != NULL){
    for(i = 1; i < vars->seq_length; i++){
      for(j = i; j <= vars->seq_length; j++){
        ij = vars->my_iindx[i] - j;
        if(!vars->Q_M[ij]) continue;
        for(cnt1 = vars->k_min_values_m[ij]; cnt1 <= vars->k_max_values_m[ij]; cnt1++)
          if(vars->l_min_values_m[ij][cnt1] < INF){
            vars->Q_M[ij][cnt1] += vars->l_min_values_m[ij][cnt1]/2;
            free(vars->Q_M[ij][cnt1]);
          }
        if(vars->k_min_values_m[ij] < INF){
          vars->Q_M[ij] += vars->k_min_values_m[ij];
          free(vars->Q_M[ij]);
          vars->l_min_values_m[ij] += vars->k_min_values_m[ij];
          vars->l_max_values_m[ij] += vars->k_min_values_m[ij];
          free(vars->l_min_values_m[ij]);
          free(vars->l_max_values_m[ij]);
        }
      }
    }
    free(vars->Q_M);
    free(vars->l_min_values_m);
    free(vars->l_max_values_m);
    free(vars->k_min_values_m);
    free(vars->k_max_values_m);
  }
#ifdef _OPENMP
  }
  #pragma omp section
  {
#endif
  if(vars->Q_M1 != NULL){
    for(i = 1; i < vars->seq_length; i++){
      for(j = i; j <= vars->seq_length; j++){
        ij = vars->jindx[j] + i;
        if(!vars->Q_M1[ij]) continue;
        for(cnt1 = vars->k_min_values_m1[ij]; cnt1 <= vars->k_max_values_m1[ij]; cnt1++)
          if(vars->l_min_values_m1[ij][cnt1] < INF){
            vars->Q_M1[ij][cnt1] += vars->l_min_values_m1[ij][cnt1]/2;
            free(vars->Q_M1[ij][cnt1]);
          }
        if(vars->k_min_values_m1[ij] < INF){
          vars->Q_M1[ij] += vars->k_min_values_m1[ij];
          free(vars->Q_M1[ij]);
          vars->l_min_values_m1[ij] += vars->k_min_values_m1[ij];
          vars->l_max_values_m1[ij] += vars->k_min_values_m1[ij];
          free(vars->l_min_values_m1[ij]);
          free(vars->l_max_values_m1[ij]);
        }
      }
    }
    free(vars->Q_M1);
    free(vars->l_min_values_m1);
    free(vars->l_max_values_m1);
    free(vars->k_min_values_m1);
    free(vars->k_max_values_m1);
  }
#ifdef _OPENMP
  }
  #pragma omp section
  {
#endif
  if(vars->Q_M2 != NULL){
    for(i = 1; i < vars->seq_length-TURN-1; i++){
      if(!vars->Q_M2[i]) continue;
      for(cnt1 = vars->k_min_values_m2[i]; cnt1 <= vars->k_max_values_m2[i]; cnt1++)
        if(vars->l_min_values_m2[i][cnt1] < INF){
          vars->Q_M2[i][cnt1] += vars->l_min_values_m2[i][cnt1]/2;
          free(vars->Q_M2[i][cnt1]);
        }
      if(vars->k_min_values_m2[i] < INF){
        vars->Q_M2[i] += vars->k_min_values_m2[i];
        free(vars->Q_M2[i]);
        vars->l_min_values_m2[i] += vars->k_min_values_m2[i];
        vars->l_max_values_m2[i] += vars->k_min_values_m2[i];
        free(vars->l_min_values_m2[i]);
        free(vars->l_max_values_m2[i]);
      }
    }
    free(vars->Q_M2);
    free(vars->l_min_values_m2);
    free(vars->l_max_values_m2);
    free(vars->k_min_values_m2);
    free(vars->k_max_values_m2);
  }
#ifdef _OPENMP
  }
  #pragma omp section
  {
#endif
  if(vars->Q_c != NULL){
    for(cnt1 = vars->k_min_values_qc; cnt1 <= vars->k_max_values_qc; cnt1++)
      if(vars->l_min_values_qc[cnt1] < INF){
        vars->Q_c[cnt1] += vars->l_min_values_qc[cnt1]/2;
        free(vars->Q_c[cnt1]);
      }
    if(vars->k_min_values_qc < INF){
      vars->Q_c += vars->k_min_values_qc;
      free(vars->Q_c);
      vars->l_min_values_qc += vars->k_min_values_qc;
      vars->l_max_values_qc += vars->k_min_values_qc;
      free(vars->l_min_values_qc);
      free(vars->l_max_values_qc);
    }
  }
#ifdef _OPENMP
  }
  #pragma omp section
  {
#endif
  if(vars->Q_cI != NULL){
    for(cnt1 = vars->k_min_values_qcI; cnt1 <= vars->k_max_values_qcI; cnt1++)
      if(vars->l_min_values_qcI[cnt1] < INF){
        vars->Q_cI[cnt1] += vars->l_min_values_qcI[cnt1]/2;
        free(vars->Q_cI[cnt1]);
      }
    if(vars->k_min_values_qcI < INF){
      vars->Q_cI += vars->k_min_values_qcI;
      free(vars->Q_cI);
      vars->l_min_values_qcI += vars->k_min_values_qcI;
      vars->l_max_values_qcI += vars->k_min_values_qcI;
      free(vars->l_min_values_qcI);
      free(vars->l_max_values_qcI);
    }
  }
#ifdef _OPENMP
  }
  #pragma omp section
  {
#endif
  if(vars->Q_cH != NULL){
    for(cnt1 = vars->k_min_values_qcH; cnt1 <= vars->k_max_values_qcH; cnt1++)
      if(vars->l_min_values_qcH[cnt1] < INF){
        vars->Q_cH[cnt1] += vars->l_min_values_qcH[cnt1]/2;
        free(vars->Q_cH[cnt1]);
      }
    if(vars->k_min_values_qcH < INF){
      vars->Q_cH += vars->k_min_values_qcH;
      free(vars->Q_cH);
      vars->l_min_values_qcH += vars->k_min_values_qcH;
      vars->l_max_values_qcH += vars->k_min_values_qcH;
      free(vars->l_min_values_qcH);
      free(vars->l_max_values_qcH);
    }
  }
#ifdef _OPENMP
  }
  #pragma omp section
  {
#endif
  if(vars->Q_cM != NULL){
    for(cnt1 = vars->k_min_values_qcM; cnt1 <= vars->k_max_values_qcM; cnt1++)
      if(vars->l_min_values_qcM[cnt1] < INF){
        vars->Q_cM[cnt1] += vars->l_min_values_qcM[cnt1]/2;
        free(vars->Q_cM[cnt1]);
      }
    if(vars->k_min_values_qcM < INF){
      vars->Q_cM += vars->k_min_values_qcM;
      free(vars->Q_cM);
      vars->l_min_values_qcM += vars->k_min_values_qcM;
      vars->l_max_values_qcM += vars->k_min_values_qcM;
      free(vars->l_min_values_qcM);
      free(vars->l_max_values_qcM);
    }
  }

#ifdef _OPENMP
  }
  #pragma omp section
  {
#endif

  if(vars->sequence != NULL)   free(vars->sequence);
  if(vars->reference_pt1 != NULL) free(vars->reference_pt1);
  if(vars->reference_pt2 != NULL) free(vars->reference_pt2);
  if(vars->referenceBPs1 != NULL) free(vars->referenceBPs1);
  if(vars->referenceBPs2 != NULL) free(vars->referenceBPs2);
  if(vars->ptype != NULL)      free(vars->ptype);
  if(vars->scale != NULL)      free(vars->scale);

  if(vars->S != NULL)          free(vars->S);
  if(vars->S1 != NULL)         free(vars->S1);
  if(vars->mm1 != NULL)           free(vars->mm1);
  if(vars->mm2 != NULL)           free(vars->mm2);
  if(vars->bpdist != NULL)        free(vars->bpdist);
#ifdef _OPENMP
  }
  }
#endif

  if(vars->my_iindx != NULL)   free(vars->my_iindx);
  if(vars->jindx != NULL)      free(vars->jindx);

  free(vars);
}


PRIVATE void initialize_TwoDpfold_vars(TwoDpfold_vars *vars){
  scale_pf_params2(vars);
  make_pair_matrix();
}

PUBLIC FLT_OR_DBL **TwoDpfold(TwoDpfold_vars *vars, int distance1, int distance2){
  unsigned int  i, d1, d2;
  unsigned int  maxD1 = 0;
  unsigned int  maxD2 = 0;
  unsigned int  mm;
  int           cnt1, cnt2;

  FLT_OR_DBL **output;

  initialize_TwoDpfold_vars(vars);

  vars->S   = encode_sequence(vars->sequence, 0);
  vars->S1  = encode_sequence(vars->sequence, 1);
  make_ptypes2(vars);

  for(i=1; i<=(unsigned int)vars->reference_pt1[0]; i++)
    if(i < (unsigned int)vars->reference_pt1[i]) maxD1++;
  for(i=1; i<=(unsigned int)vars->reference_pt2[0]; i++)
    if(i < (unsigned int)vars->reference_pt2[i]) maxD2++;
  mm    = maximumMatching(vars->sequence);
  maxD1 += mm;
  maxD2 += mm;

  if(distance1 >= 0){
    if((unsigned int)distance1 > maxD1)
      fprintf(stderr, "limiting maximum basepair distance 1 to %u\n", maxD1);
    maxD1 = (unsigned int)distance1;
  }

  if(distance2 >= 0){
    if((unsigned int)distance2 > maxD2)
      fprintf(stderr, "limiting maximum basepair distance 2 to %u\n", maxD2);
    maxD2 = (unsigned int)distance2;
  }
  vars->maxD1 = maxD1;
  vars->maxD2 = maxD2;


  output = (FLT_OR_DBL **) space(sizeof(FLT_OR_DBL*) * (maxD1+1));
  pf2D_linear(vars);
  int ndx = vars->my_iindx[1] - vars->seq_length;
  for(cnt1 = vars->k_min_values[ndx]; cnt1 <= MIN2(vars->k_max_values[ndx], vars->maxD1); cnt1++){
    output[cnt1] = (FLT_OR_DBL *)space((vars->maxD2+1)*sizeof(FLT_OR_DBL));
    for(cnt2 = vars->l_min_values[ndx][cnt1]; cnt2 <= MIN2(vars->l_max_values[ndx][cnt1], vars->maxD2); cnt2+=2){
      output[cnt1][cnt2] = vars->Q[ndx][cnt1][cnt2/2];
    }
  }
  return output;
}

PUBLIC FLT_OR_DBL **TwoDpfold_circ(TwoDpfold_vars *vars, int distance1, int distance2){
  unsigned int i, d1, d2;
  unsigned int maxD1 = 0;
  unsigned int maxD2 = 0;
  unsigned int mm;
  int           cnt1, cnt2;
  FLT_OR_DBL **output;

  initialize_TwoDpfold_vars(vars);

  vars->S   = encode_sequence(vars->sequence, 0);
  vars->S1  = encode_sequence(vars->sequence, 1);
  make_ptypes2(vars);

  for(i=1; i<=(unsigned int)vars->reference_pt1[0]; i++)
    if(i < (unsigned int)vars->reference_pt1[i]) maxD1++;
  for(i=1; i<=(unsigned int)vars->reference_pt2[0]; i++)
    if(i < (unsigned int)vars->reference_pt2[i]) maxD2++;
  mm = maximumMatching(vars->sequence);
  maxD1 += mm;
  maxD2 += mm;

  if(distance1 >= 0){
    if((unsigned int)distance1 > maxD1)
      fprintf(stderr, "limiting maximum basepair distance 1 to %u\n", maxD1);
    maxD1 = (unsigned int)distance1;
  }

  if(distance2 >= 0){
    if((unsigned int)distance2 > maxD2)
      fprintf(stderr, "limiting maximum basepair distance 2 to %u\n", maxD2);
    maxD2 = (unsigned int)distance2;
  }
  vars->maxD1 = maxD1;
  vars->maxD2 = maxD2;

  output = (FLT_OR_DBL **) space(sizeof(FLT_OR_DBL*) * (maxD1+1));
  pf2D_linear(vars);
  pf2D_circ(vars);

  for(cnt1 = vars->k_min_values_qc; cnt1 <= MIN2(vars->k_max_values_qc, vars->maxD1); cnt1++){
    output[cnt1] = (FLT_OR_DBL *)space((vars->maxD2+1)*sizeof(FLT_OR_DBL));
    for(cnt2 = vars->l_min_values_qc[cnt1]; cnt2 <= MIN2(vars->l_max_values_qc[cnt1], vars->maxD2); cnt2+=2){
      output[cnt1][cnt2] = vars->Q_c[cnt1][cnt2/2];
    }
  }
  return output;
}

PRIVATE void pf2D_linear(TwoDpfold_vars *vars){

  unsigned int  d, i, j, ij, seq_length, maxD1, maxD2, *mm1, *mm2, *bpdist;
  int     *my_iindx, *jindx, circ, cnt1, cnt2, cnt3, cnt4;
  double  max_real;
  short *S1, *reference_pt1, *reference_pt2;
  unsigned int   *referenceBPs1, *referenceBPs2;
  char  *sequence, *ptype;
  FLT_OR_DBL  *scale;
  pf_paramT   *pf_params;     /* holds all [unscaled] pf parameters */

  max_real = (sizeof(FLT_OR_DBL) == sizeof(float)) ? FLT_MAX : DBL_MAX;

  pf_params     = vars->pf_params;
  sequence      = vars->sequence;
  seq_length    = vars->seq_length;
  maxD1         = vars->maxD1;
  maxD2         = vars->maxD2;
  S1            = vars->S1;
  ptype         = vars->ptype;
  scale         = vars->scale;
  reference_pt1 = vars->reference_pt1;
  reference_pt2 = vars->reference_pt2;
  my_iindx      = vars->my_iindx;
  jindx         = vars->jindx;
  referenceBPs1 = vars->referenceBPs1;
  referenceBPs2 = vars->referenceBPs2;
  dangles       = vars->dangles;
  circ          = vars->circ;
  mm1           = vars->mm1;
  mm2           = vars->mm2;
  bpdist        = vars->bpdist;

  FLT_OR_DBL  ***Q, ***Q_B, ***Q_M, ***Q_M1;
  int         *k_min_values, *k_max_values, *k_min_values_m, *k_max_values_m,*k_min_values_m1, *k_max_values_m1,*k_min_values_b, *k_max_values_b;
  int         **l_min_values, **l_max_values, **l_min_values_m, **l_max_values_m,**l_min_values_m1, **l_max_values_m1,**l_min_values_b, **l_max_values_b;

  Q = vars->Q;
  k_min_values = vars->k_min_values;
  k_max_values = vars->k_max_values;
  l_min_values = vars->l_min_values;
  l_max_values = vars->l_max_values;

  Q_B = vars->Q_B;
  k_min_values_b = vars->k_min_values_b;
  k_max_values_b = vars->k_max_values_b;
  l_min_values_b = vars->l_min_values_b;
  l_max_values_b = vars->l_max_values_b;

  Q_M = vars->Q_M;
  k_min_values_m = vars->k_min_values_m;
  k_max_values_m = vars->k_max_values_m;
  l_min_values_m = vars->l_min_values_m;
  l_max_values_m = vars->l_max_values_m;

  Q_M1 = vars->Q_M1;
  k_min_values_m1 = vars->k_min_values_m1;
  k_max_values_m1 = vars->k_max_values_m1;
  l_min_values_m1 = vars->l_min_values_m1;
  l_max_values_m1 = vars->l_max_values_m1;



  /*array initialization ; qb,qm,q
    qb,qm,q (i,j) are stored as ((n+1-i)*(n-i) div 2 + n+1-j */

  for (j = 1; j<=seq_length; j++)
    for (i=(j>TURN?(j-TURN):1); i<=j; i++){
      ij                        = my_iindx[i]-j;
      vars->k_min_values[ij]    = 0;
      vars->k_max_values[ij]    = 0;
      vars->l_min_values[ij]    = (int *)space(sizeof(int));
      vars->l_max_values[ij]    = (int *)space(sizeof(int));
      vars->l_min_values[ij][0] = 0;
      vars->l_max_values[ij][0] = 0;
      vars->Q[ij]               = (FLT_OR_DBL **) space(sizeof(FLT_OR_DBL *));
      vars->Q[ij][0]            = (FLT_OR_DBL *)  space(sizeof(FLT_OR_DBL));
      vars->Q[ij][0][0]         = 1.0 * scale[j-i+1];
    }


  for (d = TURN+2; d <= seq_length; d++) { /* i,j in [1..seq_length] */
#ifdef _OPENMP
  #pragma omp parallel for private(i, j, ij, cnt1, cnt2, cnt3, cnt4)
#endif
    for (j = d; j <= seq_length; j++) {
      unsigned int k,l, kl, u, ii, dij;
      int no_close, type, type_2, tt, da, db, base_da, base_db;
      FLT_OR_DBL  temp2, Qmax=0, aux_en;
      unsigned int ti,tj;

      i     = j-d+1;
      ij    = my_iindx[i]-j;
      dij   = j - i - 1;
      type  = ptype[ij];


      no_close = (((type==3)||(type==4))&&no_closingGU);

      if (type) {   /* we have a pair */

        int k_min_b, k_max_b, l_min_b, l_max_b;
        int k_min_post_b, k_max_post_b, *l_min_post_b, *l_max_post_b;
        int update_b = 0;

        if(!Q_B[ij]){
          update_b = 1;
          k_min_b = l_min_b = 0;
          k_max_b = mm1[ij] + referenceBPs1[ij];
          l_max_b = mm2[ij] + referenceBPs2[ij];

          prepareBoundaries(k_min_b,
                            k_max_b,
                            l_min_b,
                            l_max_b,
                            bpdist[ij],
                            &vars->k_min_values_b[ij],
                            &vars->k_max_values_b[ij],
                            &vars->l_min_values_b[ij],
                            &vars->l_max_values_b[ij]
                            );
          preparePosteriorBoundaries( vars->k_max_values_b[ij] - vars->k_min_values_b[ij] + 1,
                                      vars->k_min_values_b[ij],
                                      &k_min_post_b,
                                      &k_max_post_b,
                                      &l_min_post_b,
                                      &l_max_post_b
                                  );

          prepareArray( &vars->Q_B[ij],
                        vars->k_min_values_b[ij],
                        vars->k_max_values_b[ij],
                        vars->l_min_values_b[ij],
                        vars->l_max_values_b[ij]
                    );
        }


        /* hairpin ----------------------------------------------*/

        /* get distance to reference if closing the hairpin
        *  d1a = dbp(T1_{i,j}, {i,j})
        */
        base_da = ((unsigned int)reference_pt1[i] != j) ? 1 : -1;
        base_db = ((unsigned int)reference_pt2[i] != j) ? 1 : -1;

        da = base_da + referenceBPs1[ij];
        db = base_db + referenceBPs2[ij];

        if(!no_close)
          if((da >= 0) && (db >= 0))
            if(((unsigned int)da<=maxD1) && ((unsigned int)db <= maxD2)){
              vars->Q_B[ij][da][db/2] = exp_E_Hairpin(dij, type, S1[i+1], S1[j-1], sequence+i-1, pf_params) * scale[dij+2];
              if(update_b){
                updatePosteriorBoundaries( da,
                                           db,
                                           &k_min_post_b,
                                           &k_max_post_b,
                                           &l_min_post_b,
                                           &l_max_post_b
                                         );
              }
            }
        /*--------------------------------------------------------
          check for elementary structures involving more than one
          closing pair.
        --------------------------------------------------------*/
        for (k = i+1; k <= MIN2(j-2-TURN,i+MAXLOOP+1) ; k++) {
          unsigned int minl, ln_pre;
          minl = k + TURN + 1;
          ln_pre = dij + k;
          if(ln_pre > minl + MAXLOOP) minl = ln_pre - MAXLOOP - 1;
          for (l = minl; l < j; l++) {
            kl = my_iindx[k] - l;
            if(!vars->Q_B[kl]) continue;
            type_2 = ptype[kl];

            if (type_2==0) continue;
            type_2 = rtype[type_2];

            /* get distance to reference if closing the interior loop
            *  d2 = dbp(S_{i,j}, S_{k,l} + {i,j})
            */
            da = base_da + referenceBPs1[ij] - referenceBPs1[kl];
            db = base_db + referenceBPs2[ij] - referenceBPs2[kl];

            aux_en = exp_E_IntLoop(k-i-1, j-l-1, type, type_2, S1[i+1], S1[j-1], S1[k-1], S1[l+1], pf_params) * scale[k-i+j-l];
            for(cnt1 = vars->k_min_values_b[kl]; cnt1 <= vars->k_max_values_b[kl]; cnt1++)
              for(cnt2 = vars->l_min_values_b[kl][cnt1]; cnt2 <= vars->l_max_values_b[kl][cnt1]; cnt2+=2){
                vars->Q_B[ij][cnt1 + da][(cnt2 + db)/2] += vars->Q_B[kl][cnt1][cnt2/2] * aux_en;
                if(update_b){
                  updatePosteriorBoundaries( da + cnt1,
                                             db + cnt2,
                                             &k_min_post_b,
                                             &k_max_post_b,
                                             &l_min_post_b,
                                             &l_max_post_b
                                           );
                }
              }

          } /* end l-loop */
        } /* end k-loop */

        /* multi-loop contribution ------------------------*/
        if(!no_close){
          for(u=i+TURN+2; u<j-TURN-2;u++){
            if(!vars->Q_M[my_iindx[i+1]-u]) continue;
            if(!vars->Q_M1[jindx[j-1]+u+1]) continue;
            /* get distance to reference if closing the multiloop
            *  dist3 = dbp(S_{i,j}, {i,j} + S_{i+1,u} + S_{u+1,j-1})
            */

            da = base_da + referenceBPs1[ij] - referenceBPs1[my_iindx[i+1]-u] - referenceBPs1[my_iindx[u+1]-j+1];
            db = base_db + referenceBPs2[ij] - referenceBPs2[my_iindx[i+1]-u] - referenceBPs2[my_iindx[u+1]-j+1];

            tt = rtype[type];
            temp2 = pf_params->expMLclosing * exp_E_MLstem(tt, S1[j-1], S1[i+1], pf_params) * scale[2];

            for(cnt1 = vars->k_min_values_m[my_iindx[i+1]-u]; cnt1 <= vars->k_max_values_m[my_iindx[i+1]-u]; cnt1++)
              for(cnt2 = vars->l_min_values_m[my_iindx[i+1]-u][cnt1]; cnt2 <= vars->l_max_values_m[my_iindx[i+1]-u][cnt1]; cnt2+=2){
                for(cnt3 = vars->k_min_values_m1[jindx[j-1]+u+1]; cnt3 <= vars->k_max_values_m1[jindx[j-1]+u+1]; cnt3++)
                  for(cnt4 = vars->l_min_values_m1[jindx[j-1]+u+1][cnt3]; cnt4 <= vars->l_max_values_m1[jindx[j-1]+u+1][cnt3]; cnt4+=2){
                    vars->Q_B[ij][cnt1 + cnt3 + da][(cnt2 + cnt4 + db)/2] += vars->Q_M[my_iindx[i+1]-u][cnt1][cnt2/2] * vars->Q_M1[jindx[j-1]+u+1][cnt3][cnt4/2] * temp2;
                    if(update_b){
                      updatePosteriorBoundaries( cnt1 + cnt3 + da,
                                                 cnt2 + cnt4 + db,
                                                 &k_min_post_b,
                                                 &k_max_post_b,
                                                 &l_min_post_b,
                                                 &l_max_post_b
                                               );
                    }
                  }

              }

          }
        }

        if(update_b){
          adjustArrayBoundaries(&vars->Q_B[ij],
                                &vars->k_min_values_b[ij],
                                &vars->k_max_values_b[ij],
                                &vars->l_min_values_b[ij],
                                &vars->l_max_values_b[ij],
                                k_min_post_b,
                                k_max_post_b,
                                l_min_post_b,
                                l_max_post_b
                                );
        }
      } /* end >> if (pair) << */

      /* free ends ? -----------------------------------------*/

      int k_min_m, k_max_m, l_min_m, l_max_m;
      int k_min_post_m, k_max_post_m, *l_min_post_m, *l_max_post_m;
      int update_m = 0;
      int k_min_m1, k_max_m1, l_min_m1, l_max_m1;
      int k_min_post_m1, k_max_post_m1, *l_min_post_m1, *l_max_post_m1;
      int update_m1 = 0;

      if(!Q_M[ij]){
        update_m = 1;
        k_min_m = l_min_m = 0;
        k_max_m = mm1[ij] + referenceBPs1[ij];
        l_max_m = mm2[ij] + referenceBPs2[ij];

        prepareBoundaries(k_min_m,
                          k_max_m,
                          l_min_m,
                          l_max_m,
                          bpdist[ij],
                          &vars->k_min_values_m[ij],
                          &vars->k_max_values_m[ij],
                          &vars->l_min_values_m[ij],
                          &vars->l_max_values_m[ij]
                          );
        preparePosteriorBoundaries( vars->k_max_values_m[ij] - vars->k_min_values_m[ij] + 1,
                                    vars->k_min_values_m[ij],
                                    &k_min_post_m,
                                    &k_max_post_m,
                                    &l_min_post_m,
                                    &l_max_post_m
                                );

        prepareArray( &vars->Q_M[ij],
                      vars->k_min_values_m[ij],
                      vars->k_max_values_m[ij],
                      vars->l_min_values_m[ij],
                      vars->l_max_values_m[ij]
                  );
      }
      if(!Q_M1[jindx[j]+i]){
        update_m1 = 1;
        k_min_m1 = l_min_m1 = 0;
        k_max_m1 = mm1[ij] + referenceBPs1[ij];
        l_max_m1 = mm2[ij] + referenceBPs2[ij];

        prepareBoundaries(k_min_m1,
                          k_max_m1,
                          l_min_m1,
                          l_max_m1,
                          bpdist[ij],
                          &vars->k_min_values_m1[jindx[j]+i],
                          &vars->k_max_values_m1[jindx[j]+i],
                          &vars->l_min_values_m1[jindx[j]+i],
                          &vars->l_max_values_m1[jindx[j]+i]
                          );
        preparePosteriorBoundaries( vars->k_max_values_m1[jindx[j]+i] - vars->k_min_values_m1[jindx[j]+i] + 1,
                                    vars->k_min_values_m1[jindx[j]+i],
                                    &k_min_post_m1,
                                    &k_max_post_m1,
                                    &l_min_post_m1,
                                    &l_max_post_m1
                                );

        prepareArray( &vars->Q_M1[jindx[j]+i],
                      vars->k_min_values_m1[jindx[j]+i],
                      vars->k_max_values_m1[jindx[j]+i],
                      vars->l_min_values_m1[jindx[j]+i],
                      vars->l_max_values_m1[jindx[j]+i]
                  );
      }


      /* j is unpaired */
      da = referenceBPs1[ij] - referenceBPs1[ij+1];
      db = referenceBPs2[ij] - referenceBPs2[ij+1];
      if(vars->Q_M[ij+1])
        for(cnt1 = vars->k_min_values_m[ij+1]; cnt1 <= vars->k_max_values_m[ij+1]; cnt1++){
          for(cnt2 = vars->l_min_values_m[ij+1][cnt1]; cnt2 <= vars->l_max_values_m[ij+1][cnt1]; cnt2+=2){
            vars->Q_M[ij][cnt1 + da][(cnt2 + db)/2] += vars->Q_M[ij+1][cnt1][cnt2/2] * pf_params->expMLbase * scale[1];
            if(update_m){
              updatePosteriorBoundaries(cnt1 + da,
                                        cnt2 + db,
                                        &k_min_post_m,
                                        &k_max_post_m,
                                        &l_min_post_m,
                                        &l_max_post_m
                                        );
            }
          }
        }
      if(vars->Q_M1[jindx[j-1]+i])
        for(cnt1 = vars->k_min_values_m1[jindx[j-1]+i]; cnt1 <= vars->k_max_values_m1[jindx[j-1]+i]; cnt1++)
          for(cnt2 = vars->l_min_values_m1[jindx[j-1]+i][cnt1]; cnt2 <= vars->l_max_values_m1[jindx[j-1]+i][cnt1]; cnt2+=2){
            vars->Q_M1[jindx[j]+i][cnt1 + da][(cnt2 + db)/2] += vars->Q_M1[jindx[j-1]+i][cnt1][cnt2/2] * pf_params->expMLbase * scale[1];
            if(update_m1){
              updatePosteriorBoundaries(cnt1 + da,
                                        cnt2 + db,
                                        &k_min_post_m1,
                                        &k_max_post_m1,
                                        &l_min_post_m1,
                                        &l_max_post_m1
                                        );
            }
          }


      /* j pairs with i */
      if((!no_close) && type && vars->Q_B[ij]){
        FLT_OR_DBL aux_en = exp_E_MLstem(type, (i>1) || circ ? S1[i-1] : -1, (j<seq_length) || circ ? S1[j+1] : -1, pf_params);

        for(cnt1 = vars->k_min_values_b[ij]; cnt1 <= vars->k_max_values_b[ij]; cnt1++)
          for(cnt2 = vars->l_min_values_b[ij][cnt1]; cnt2 <= vars->l_max_values_b[ij][cnt1]; cnt2+=2){
            vars->Q_M[ij][cnt1][cnt2/2] += vars->Q_B[ij][cnt1][cnt2/2] * aux_en;
            if(update_m){
              updatePosteriorBoundaries(cnt1,
                                        cnt2,
                                        &k_min_post_m,
                                        &k_max_post_m,
                                        &l_min_post_m,
                                        &l_max_post_m
                                        );
            }
            vars->Q_M1[jindx[j]+i][cnt1][cnt2/2] += vars->Q_B[ij][cnt1][cnt2/2] * aux_en;
            if(update_m1){
              updatePosteriorBoundaries(cnt1,
                                        cnt2,
                                        &k_min_post_m1,
                                        &k_max_post_m1,
                                        &l_min_post_m1,
                                        &l_max_post_m1
                                        );
            }
          }
      }

      /* j pairs with k: i<k<j */
      ii = my_iindx[i];
      for (k=i+1; k<=j; k++){
        if(!vars->Q_B[my_iindx[k]-j]) continue;
        tt = ptype[my_iindx[k]-j];
        /* add contributions of QM(i,k-1)*QB(k,j)*e^b and
        *  e^((k-i) * c) * QB(k,j) * e^b
        *  therefor we need d1a = dbp(T1_{i,j}, T1_{i,k-1} + T1_{k,j}),
        *  d1b = dbp(T2_{i,j}, T2_{i,k-1} + T2_{k,j})
        *  d1c = dbp(T1_{i,j}, T1_{k,j})circ = 0;
        *  d1d = dbp(T2_{i,j}, T2_{k,j})
        */
        da = referenceBPs1[ij] - referenceBPs1[my_iindx[k]-j];
        db = referenceBPs2[ij] - referenceBPs2[my_iindx[k]-j];

        temp2 = exp_E_MLstem(tt, S1[k-1], (j<seq_length) || circ ? S1[j+1] : -1, pf_params);

        for(cnt1 = vars->k_min_values_b[my_iindx[k]-j]; cnt1 <= vars->k_max_values_b[my_iindx[k]-j]; cnt1++)
          for(cnt2 = vars->l_min_values_b[my_iindx[k]-j][cnt1]; cnt2 <= vars->l_max_values_b[my_iindx[k]-j][cnt1]; cnt2+=2){
            vars->Q_M[ij][cnt1 + da][(cnt2 + db)/2] += vars->Q_B[my_iindx[k]-j][cnt1][cnt2/2] * pow(pf_params->expMLbase, (double)(k-i)) * scale[k-i] * temp2;
            if(update_m){
              updatePosteriorBoundaries(cnt1 + da,
                                        cnt2 + db,
                                        &k_min_post_m,
                                        &k_max_post_m,
                                        &l_min_post_m,
                                        &l_max_post_m
                                        );
            }


          }

        if(!vars->Q_M[ii-k+1]) continue;
        da -= referenceBPs1[ii-k+1];
        db -= referenceBPs2[ii-k+1];

        for(cnt1 = vars->k_min_values_m[ii-k+1]; cnt1 <= vars->k_max_values_m[ii-k+1]; cnt1++)
          for(cnt2 = vars->l_min_values_m[ii-k+1][cnt1]; cnt2 <= vars->l_max_values_m[ii-k+1][cnt1]; cnt2+=2)
            for(cnt3 = vars->k_min_values_b[my_iindx[k]-j]; cnt3 <= vars->k_max_values_b[my_iindx[k]-j]; cnt3++)
              for(cnt4 = vars->l_min_values_b[my_iindx[k]-j][cnt3]; cnt4 <= vars->l_max_values_b[my_iindx[k]-j][cnt3]; cnt4+=2){
                vars->Q_M[ij][cnt1 + cnt3 + da][(cnt2 + cnt4 + db)/2] += vars->Q_M[ii-k+1][cnt1][cnt2/2] * vars->Q_B[my_iindx[k]-j][cnt3][cnt4/2] * temp2;
                if(update_m){
                  updatePosteriorBoundaries(cnt1 + cnt3 + da,
                                            cnt2 + cnt4 + db,
                                            &k_min_post_m,
                                            &k_max_post_m,
                                            &l_min_post_m,
                                            &l_max_post_m
                                            );
                }
              }
      }

      if(update_m){
        adjustArrayBoundaries(&vars->Q_M[ij],
                              &vars->k_min_values_m[ij],
                              &vars->k_max_values_m[ij],
                              &vars->l_min_values_m[ij],
                              &vars->l_max_values_m[ij],
                              k_min_post_m,
                              k_max_post_m,
                              l_min_post_m,
                              l_max_post_m
                              );
      }
      if(update_m1){
        adjustArrayBoundaries(&vars->Q_M1[jindx[j]+i],
                              &vars->k_min_values_m1[jindx[j]+i],
                              &vars->k_max_values_m1[jindx[j]+i],
                              &vars->l_min_values_m1[jindx[j]+i],
                              &vars->l_max_values_m1[jindx[j]+i],
                              k_min_post_m1,
                              k_max_post_m1,
                              l_min_post_m1,
                              l_max_post_m1
                              );
      }

      int k_min, k_max, l_min, l_max;
      int k_min_post, k_max_post, *l_min_post, *l_max_post;
      int update_q = 0;
      if(!Q[ij]){
        update_q = 1;
        k_min = l_min = 0;
        k_max = mm1[ij] + referenceBPs1[ij];
        l_max = mm2[ij] + referenceBPs2[ij];

        prepareBoundaries(k_min,
                          k_max,
                          l_min,
                          l_max,
                          bpdist[ij],
                          &vars->k_min_values[ij],
                          &vars->k_max_values[ij],
                          &vars->l_min_values[ij],
                          &vars->l_max_values[ij]
                          );
        preparePosteriorBoundaries( vars->k_max_values[ij] - vars->k_min_values[ij] + 1,
                                    vars->k_min_values[ij],
                                    &k_min_post,
                                    &k_max_post,
                                    &l_min_post,
                                    &l_max_post
                                );

        prepareArray( &vars->Q[ij],
                      vars->k_min_values[ij],
                      vars->k_max_values[ij],
                      vars->l_min_values[ij],
                      vars->l_max_values[ij]
                  );
      }
      aux_en = 1.0;
      if (type && vars->Q_B[ij]){
        aux_en *= exp_E_ExtLoop(type, (i>1) || circ ? S1[i-1] : -1, (j < seq_length) || circ ? S1[j+1] : -1, pf_params);
      }

      if(vars->Q_B[ij])
        for(cnt1 = vars->k_min_values_b[ij]; cnt1 <= vars->k_max_values_b[ij]; cnt1++)
          for(cnt2 = vars->l_min_values_b[ij][cnt1]; cnt2 <= vars->l_max_values_b[ij][cnt1]; cnt2+=2){
            vars->Q[ij][cnt1][cnt2/2] += vars->Q_B[ij][cnt1][cnt2/2] * aux_en;
            if(update_q){
              updatePosteriorBoundaries(cnt1,
                                        cnt2,
                                        &k_min_post,
                                        &k_max_post,
                                        &l_min_post,
                                        &l_max_post
                                        );
            }

          }

      /* j is unpaired */
      /* da = dbp(T1_{i,j}, T1_{i,j-1})
      *  db = dbp(T2_{i,j}, T2_{i,j-1})
      */
      da = referenceBPs1[ij] - referenceBPs1[ij+1];
      db = referenceBPs2[ij] - referenceBPs2[ij+1];
      if(vars->Q[ij+1])
        for(cnt1 = vars->k_min_values[ij+1]; cnt1 <= vars->k_max_values[ij+1]; cnt1++)
          for(cnt2 = vars->l_min_values[ij+1][cnt1]; cnt2 <= vars->l_max_values[ij+1][cnt1]; cnt2+=2){
            vars->Q[ij][cnt1 + da][(cnt2 + db)/2] += vars->Q[ij+1][cnt1][cnt2/2] * scale[1];
            if(update_q){
              updatePosteriorBoundaries(cnt1 + da,
                                        cnt2 + db,
                                        &k_min_post,
                                        &k_max_post,
                                        &l_min_post,
                                        &l_max_post
                                        );
            }

          }

      for(k=j-TURN-1; k>i; k--){
        if(!vars->Q[my_iindx[i]-k+1]) continue;
        if(!vars->Q_B[my_iindx[k]-j]) continue;
        tt = ptype[my_iindx[k]-j];

        /* da = dbp{T1_{i,j}, T1_{k,j}
        *  db = dbp{T2_{i,j}, T2_{k,j}}
        */

        da = referenceBPs1[ij] - referenceBPs1[my_iindx[k] - j] - referenceBPs1[my_iindx[i]-k+1];
        db = referenceBPs2[ij] - referenceBPs2[my_iindx[k] - j] - referenceBPs2[my_iindx[i]-k+1];

        temp2 = exp_E_ExtLoop(tt, S1[k-1], (j<seq_length) || circ ? S1[j+1] : -1, pf_params);

        for(cnt1 = vars->k_min_values[my_iindx[i]-k+1]; cnt1 <= vars->k_max_values[my_iindx[i]-k+1]; cnt1++)
          for(cnt2 = vars->l_min_values[my_iindx[i]-k+1][cnt1]; cnt2 <= vars->l_max_values[my_iindx[i]-k+1][cnt1]; cnt2+=2)
            for(cnt3 = vars->k_min_values_b[my_iindx[k]-j]; cnt3 <= vars->k_max_values_b[my_iindx[k]-j]; cnt3++)
              for(cnt4 = vars->l_min_values_b[my_iindx[k]-j][cnt3]; cnt4 <= vars->l_max_values_b[my_iindx[k]-j][cnt3]; cnt4+=2){
                vars->Q[ij][cnt1 + cnt3 + da][(cnt2 + cnt4 + db)/2] += vars->Q[my_iindx[i]-k+1][cnt1][cnt2/2] * vars->Q_B[my_iindx[k]-j][cnt3][cnt4/2] * temp2;
                if(update_q){
                  updatePosteriorBoundaries(cnt1 + cnt3 + da,
                                            cnt2 + cnt4 + db,
                                            &k_min_post,
                                            &k_max_post,
                                            &l_min_post,
                                            &l_max_post
                                            );
                }
              }
      }

      if(update_q){
        adjustArrayBoundaries(&vars->Q[ij],
                              &vars->k_min_values[ij],
                              &vars->k_max_values[ij],
                              &vars->l_min_values[ij],
                              &vars->l_max_values[ij],
                              k_min_post,
                              k_max_post,
                              l_min_post,
                              l_max_post
                              );
      }


#if 0
      new_q_list = new_q_list_last = NULL;
      for(da = 0; da <= maxD1; da++){
        for(db = 0; db <= maxD2; db++){
          if(qq[da][db] != 0){
            insertValue(&new_q_list, &new_q_list_last, da, db, qq[da][db]);
            if(qq[da][db] > Qmax) {
              Qmax = qq[da][db];
              if (Qmax > max_real/10.)
                fprintf(stderr, "Q close to overflow: %u %u %g\n", i,j,qq[da][db]);
            }
            if(qq[da][db] >= max_real) {
              PRIVATE char msg[128];
              sprintf(msg, "overflow in pf_fold while calculating q[%u,%u]\n"
                            "use larger pf_scale", i,j);
              nrerror(msg);
            }
            qq[da][db] = 0;
          }

        }
        free(qq[da]);
        free(qqm[da]);
      }
      free(qq);
      free(qqm);
      vars->q[ij] = new_q_list;
#endif

    } /* end of j-loop */
  }
}

/* calculate partition function for circular case */
/* NOTE: this is the postprocessing step ONLY     */
/* You have to call pf2D_linear first to calculate  */
/* complete circular case!!!                      */
PRIVATE void pf2D_circ(TwoDpfold_vars *vars){

  unsigned int  d, p, q, pq, k, l, kl, u, da, db, seq_length, maxD1, maxD2, base_d1, base_d2, *mm1, *mm2, *bpdist;
  int         *my_iindx, *jindx, type, cnt1, cnt2, cnt3, cnt4;
  double      max_real;
  short       *S1, *reference_pt1, *reference_pt2;
  unsigned int  *referenceBPs1, *referenceBPs2;
  char        *sequence, *ptype;
  FLT_OR_DBL  *scale, qot;
  pf_paramT   *pf_params;     /* holds all [unscaled] pf parameters */


  max_real        =  (sizeof(FLT_OR_DBL) == sizeof(float)) ? FLT_MAX : DBL_MAX;

  pf_params       = vars->pf_params;
  sequence        = vars->sequence;
  seq_length      = vars->seq_length;
  maxD1           = vars->maxD1;
  maxD2           = vars->maxD2;
  S1              = vars->S1;
  ptype           = vars->ptype;
  scale           = vars->scale;
  reference_pt1   = vars->reference_pt1;
  reference_pt2   = vars->reference_pt2;
  my_iindx        = vars->my_iindx;
  jindx           = vars->jindx;
  referenceBPs1   = vars->referenceBPs1;
  referenceBPs2   = vars->referenceBPs2;
  dangles         = vars->dangles;
  mm1             = vars->mm1;
  mm2             = vars->mm2;
  bpdist          = vars->bpdist;

  FLT_OR_DBL      ***Q_B, ***Q_M, ***Q_M1, ***Q_M2, **Q_c, **Q_cH, **Q_cI, **Q_cM;
  int             **l_min_values_b, **l_max_values_b, **l_min_values_m, **l_max_values_m, **l_min_values_m1, **l_max_values_m1, **l_min_values_m2, **l_max_values_m2;
  int             *k_min_values_b, *k_max_values_b,*k_min_values_m, *k_max_values_m,*k_min_values_m1, *k_max_values_m1, *k_min_values_m2, *k_max_values_m2;
  int             k_min_values_qc, k_max_values_qc, k_min_values_qcH, k_max_values_qcH, k_min_values_qcI, k_max_values_qcI, k_min_values_qcM, k_max_values_qcM;
  int             *l_min_values_qc, *l_max_values_qc, *l_min_values_qcH, *l_max_values_qcH, *l_min_values_qcI, *l_max_values_qcI, *l_min_values_qcM, *l_max_values_qcM;

  Q_B             = vars->Q_B;
  l_min_values_b  = vars->l_min_values_b;
  l_max_values_b  = vars->l_max_values_b;
  k_min_values_b  = vars->k_min_values_b;
  k_max_values_b  = vars->k_max_values_b;

  Q_M             = vars->Q_M;
  l_min_values_m  = vars->l_min_values_m;
  l_max_values_m  = vars->l_max_values_m;
  k_min_values_m  = vars->k_min_values_m;
  k_max_values_m  = vars->k_max_values_m;

  Q_M1            = vars->Q_M1;
  l_min_values_m1 = vars->l_min_values_m1;
  l_max_values_m1 = vars->l_max_values_m1;
  k_min_values_m1 = vars->k_min_values_m1;
  k_max_values_m1 = vars->k_max_values_m1;

  Q_M2            = vars->Q_M2;
  l_min_values_m2 = vars->l_min_values_m2;
  l_max_values_m2 = vars->l_max_values_m2;
  k_min_values_m2 = vars->k_min_values_m2;
  k_max_values_m2 = vars->k_max_values_m2;

  Q_c             = vars->Q_c;
  l_min_values_qc = vars->l_min_values_qc;
  l_max_values_qc = vars->l_max_values_qc;
  k_min_values_qc = vars->k_min_values_qc;
  k_max_values_qc = vars->k_max_values_qc;

  Q_cI             = vars->Q_cI;
  l_min_values_qcI = vars->l_min_values_qcI;
  l_max_values_qcI = vars->l_max_values_qcI;
  k_min_values_qcI = vars->k_min_values_qcI;
  k_max_values_qcI = vars->k_max_values_qcI;

  Q_cH             = vars->Q_cH;
  l_min_values_qcH = vars->l_min_values_qcH;
  l_max_values_qcH = vars->l_max_values_qcH;
  k_min_values_qcH = vars->k_min_values_qcH;
  k_max_values_qcH = vars->k_max_values_qcH;

  Q_cM             = vars->Q_cM;
  l_min_values_qcM = vars->l_min_values_qcM;
  l_max_values_qcM = vars->l_max_values_qcM;
  k_min_values_qcM = vars->k_min_values_qcM;
  k_max_values_qcM = vars->k_max_values_qcM;




  /* construct qm2 matrix from qm1 entries  */
#ifdef _OPENMP
  #pragma omp parallel for private(d, k, l, da, db, cnt1, cnt2, cnt3, cnt4)
#endif
  for(k=1; k<seq_length-TURN-1; k++){
    int k_min_m2, k_max_m2, l_min_m2, l_max_m2;
    int k_min_post_m2, k_max_post_m2, *l_min_post_m2, *l_max_post_m2;
    int update_m2 = 0;
    if(!Q_M2[k]){
      update_m2 = 1;
      k_min_m2 = l_min_m2 = 0;
      k_max_m2 = mm1[my_iindx[k]-seq_length] + referenceBPs1[my_iindx[k] - seq_length];
      l_max_m2 = mm2[my_iindx[k]-seq_length] + referenceBPs2[my_iindx[k] - seq_length];

      prepareBoundaries(k_min_m2,
                        k_max_m2,
                        l_min_m2,
                        l_max_m2,
                        bpdist[my_iindx[k]-seq_length],
                        &vars->k_min_values_m2[k],
                        &vars->k_max_values_m2[k],
                        &vars->l_min_values_m2[k],
                        &vars->l_max_values_m2[k]
                        );
      preparePosteriorBoundaries( vars->k_max_values_m2[k] - vars->k_min_values_m2[k] + 1,
                                  vars->k_min_values_m2[k],
                                  &k_min_post_m2,
                                  &k_max_post_m2,
                                  &l_min_post_m2,
                                  &l_max_post_m2
                              );

      prepareArray( &vars->Q_M2[k],
                    vars->k_min_values_m2[k],
                    vars->k_max_values_m2[k],
                    vars->l_min_values_m2[k],
                    vars->l_max_values_m2[k]
                );
    }
    for (l=k+TURN+1; l<seq_length-TURN-1; l++){
      if(vars->Q_M1[jindx[l]+k] && vars->Q_M1[jindx[seq_length] + l + 1]){
        da = referenceBPs1[my_iindx[k]-seq_length] - referenceBPs1[my_iindx[k]-l] - referenceBPs1[my_iindx[l+1]-seq_length];
        db = referenceBPs2[my_iindx[k]-seq_length] - referenceBPs2[my_iindx[k]-l] - referenceBPs2[my_iindx[l+1]-seq_length];
        for(cnt1 = k_min_values_m1[jindx[l]+k]; cnt1 <= k_max_values_m1[jindx[l]+k]; cnt1++)
          for(cnt2 = l_min_values_m1[jindx[l]+k][cnt1]; cnt2 <= l_max_values_m1[jindx[l]+k][cnt1]; cnt2+=2){
            for(cnt3 = k_min_values_m1[jindx[seq_length] + l + 1]; cnt3 <= k_max_values_m1[jindx[seq_length] + l + 1]; cnt3++)
              for(cnt4 = l_min_values_m1[jindx[seq_length] + l + 1][cnt3]; cnt4 <= l_max_values_m1[jindx[seq_length] + l + 1][cnt3]; cnt4+=2){
/*
                printf("%d - %d: %d: %d-%d\n", k, l, cnt1, l_min_values_m1[jindx[l]+k][cnt1], l_max_values_m1[jindx[l]+k][cnt1]);
                printf("%d - %d: %d: %d-%d\n", l+1, seq_length, cnt3, l_min_values_m1[jindx[seq_length] + l + 1][cnt3], l_max_values_m1[jindx[seq_length] + l + 1][cnt3]);
                printf("%6.2f\n", Q_M1[jindx[l]+k][cnt1][cnt2/2]);
                printf("%6.2f\n", Q_M1[jindx[seq_length] + l + 1][cnt3][cnt4/2]);
                printf("%d-%d and %d-%d but %d,%d\n", vars->k_min_values_m2[k], vars->k_max_values_m2[k], vars->l_min_values_m2[k][cnt1 + cnt3 + da], vars->l_max_values_m2[k][cnt1 + cnt3 + da], cnt1 + cnt3 + da, cnt2 + cnt4 + db);
*/
                vars->Q_M2[k][cnt1 + cnt3 + da][(cnt2 + cnt4 + db)/2] += Q_M1[jindx[l]+k][cnt1][cnt2/2] * Q_M1[jindx[seq_length] + l + 1][cnt3][cnt4/2];
                if(update_m2){
                      updatePosteriorBoundaries(cnt1 + cnt3 + da,
                                                cnt2 + cnt4 + db,
                                                &k_min_post_m2,
                                                &k_max_post_m2,
                                                &l_min_post_m2,
                                                &l_max_post_m2
                                                );
                }
              }
          }
      }
    }
    if(update_m2){
      adjustArrayBoundaries(&vars->Q_M2[k],
                            &vars->k_min_values_m2[k],
                            &vars->k_max_values_m2[k],
                            &vars->l_min_values_m2[k],
                            &vars->l_max_values_m2[k],
                            k_min_post_m2,
                            k_max_post_m2,
                            l_min_post_m2,
                            l_max_post_m2
                            );
    }
  }

  base_d1 = referenceBPs1[my_iindx[1]-seq_length];
  base_d2 = referenceBPs2[my_iindx[1]-seq_length];

  int min_k, max_k, max_l, min_l;
  int min_k_real, max_k_real, min_k_real_qcH, max_k_real_qcH, min_k_real_qcI, max_k_real_qcI, min_k_real_qcM, max_k_real_qcM;
  int *min_l_real, *max_l_real, *min_l_real_qcH, *max_l_real_qcH, *min_l_real_qcI, *max_l_real_qcI,*min_l_real_qcM, *max_l_real_qcM;
  int update_c, update_cH, update_cI, update_cM;

  update_c = update_cH = update_cI, update_cM = 0;

  min_k = min_l = 0;

  max_k = mm1[my_iindx[1] - seq_length] + referenceBPs1[my_iindx[1] - seq_length];
  max_l = mm2[my_iindx[1] - seq_length] + referenceBPs2[my_iindx[1] - seq_length];

#ifdef _OPENMP
  #pragma omp sections
  {

  #pragma omp section
  {
#endif
  if(!Q_c){
    update_c = 1;
    prepareBoundaries(min_k,
                      max_k,
                      min_l,
                      max_l,
                      bpdist[my_iindx[1] - seq_length],
                      &vars->k_min_values_qc,
                      &vars->k_max_values_qc,
                      &vars->l_min_values_qc,
                      &vars->l_max_values_qc
                      );
    prepareArray( &vars->Q_c,
                  vars->k_min_values_qc,
                  vars->k_max_values_qc,
                  vars->l_min_values_qc,
                  vars->l_max_values_qc
                );
    preparePosteriorBoundaries( max_k - min_k + 1,
                                min_k,
                                &min_k_real,
                                &max_k_real,
                                &min_l_real,
                                &max_l_real
                              );
  }
#ifdef _OPENMP
  }
  #pragma omp section
  {
#endif
  if(!Q_cH){
    update_cH = 1;
    prepareBoundaries(min_k,
                      max_k,
                      min_l,
                      max_l,
                      bpdist[my_iindx[1] - seq_length],
                      &vars->k_min_values_qcH,
                      &vars->k_max_values_qcH,
                      &vars->l_min_values_qcH,
                      &vars->l_max_values_qcH
                      );
    prepareArray( &vars->Q_cH,
                  vars->k_min_values_qcH,
                  vars->k_max_values_qcH,
                  vars->l_min_values_qcH,
                  vars->l_max_values_qcH
                );
    preparePosteriorBoundaries( max_k - min_k + 1,
                                min_k,
                                &min_k_real_qcH,
                                &max_k_real_qcH,
                                &min_l_real_qcH,
                                &max_l_real_qcH
                              );
  }
#ifdef _OPENMP
  }
  #pragma omp section
  {
#endif
  if(!Q_cI){
    update_cI = 1;
    prepareBoundaries(min_k,
                      max_k,
                      min_l,
                      max_l,
                      bpdist[my_iindx[1] - seq_length],
                      &vars->k_min_values_qcI,
                      &vars->k_max_values_qcI,
                      &vars->l_min_values_qcI,
                      &vars->l_max_values_qcI
                      );
    prepareArray( &vars->Q_cI,
                  vars->k_min_values_qcI,
                  vars->k_max_values_qcI,
                  vars->l_min_values_qcI,
                  vars->l_max_values_qcI
                );
    preparePosteriorBoundaries( max_k - min_k + 1,
                                min_k,
                                &min_k_real_qcI,
                                &max_k_real_qcI,
                                &min_l_real_qcI,
                                &max_l_real_qcI
                              );
  }
#ifdef _OPENMP
  }
  #pragma omp section
  {
#endif
  if(!Q_cM){
    update_cM = 1;
    prepareBoundaries(min_k,
                      max_k,
                      min_l,
                      max_l,
                      bpdist[my_iindx[1] - seq_length],
                      &vars->k_min_values_qcM,
                      &vars->k_max_values_qcM,
                      &vars->l_min_values_qcM,
                      &vars->l_max_values_qcM
                      );
    prepareArray( &vars->Q_cM,
                  vars->k_min_values_qcM,
                  vars->k_max_values_qcM,
                  vars->l_min_values_qcM,
                  vars->l_max_values_qcM
                );
    preparePosteriorBoundaries( max_k - min_k + 1,
                                min_k,
                                &min_k_real_qcM,
                                &max_k_real_qcM,
                                &min_l_real_qcM,
                                &max_l_real_qcM
                              );
  }
#ifdef _OPENMP
  }
  }
#endif




  for (d = TURN+2; d <= seq_length; d++) /* i,j in [1..length] */
#ifdef _OPENMP
    #pragma omp parallel for private(p, q, pq, k, l, kl, u, da, db, type, cnt1, cnt2, cnt3, cnt4)
#endif
    for (q = d; q <= seq_length; q++) {
      FLT_OR_DBL qot;
      char loopseq[10];
      p = q - d + 1;
      pq = my_iindx[p]-q;
      if (!Q_B[pq]) continue;

      /* 1. get exterior hairpin contribution  */
      u = seq_length-q + p-1;
      if (u<TURN) continue;
      type = ptype[pq];
      if (!type) continue;
      if(((type==3)||(type==4))&&no_closingGU) continue;

       /* cause we want to calc the exterior loops, we need the reversed pair type from now on  */
      type=rtype[type];

      if (u<7){
        strcpy(loopseq , sequence+q-1);
        strncat(loopseq, sequence, p);
      }
      /* get distance to reference if closing the hairpin
      *  da = dbp(T1_[1,n}, T1_{p,q})
      *  db = dbp(T2_{1,n}, T2_{p,q})
      */
      da = base_d1 - referenceBPs1[pq];
      db = base_d2 - referenceBPs2[pq];
      qot = exp_E_Hairpin(u, type, S1[q+1], S1[p-1],  loopseq, pf_params) * scale[u];
      for(cnt1 = k_min_values_b[pq]; cnt1 <= k_max_values_b[pq]; cnt1++)
        for(cnt2 = l_min_values_b[pq][cnt1]; cnt2 <= l_max_values_b[pq][cnt1]; cnt2+=2){
          vars->Q_cH[cnt1 + da][(cnt2 + db)/2] += Q_B[pq][cnt1][cnt2/2] * qot;
          if(update_cH){
            updatePosteriorBoundaries(cnt1 + da,
                                      cnt2 + db,
                                      &min_k_real_qcH,
                                      &max_k_real_qcH,
                                      &min_l_real_qcH,
                                      &max_l_real_qcH
                                      );
          }
        }

      /* 2. exterior interior loops, i "define" the (k,l) pair as "outer pair"  */
      /* so "outer type" is rtype[type[k,l]] and inner type is type[p,q]        */
      for(k=q+1; k < seq_length; k++){
        unsigned int ln1, lstart, ln_pre;
        ln1 = k - q - 1;
        if(ln1+p-1>MAXLOOP) break;
        lstart = k + TURN + 1;
        ln_pre = ln1 + p + seq_length;
        if(ln_pre > lstart + MAXLOOP) lstart = ln_pre - MAXLOOP - 1;
        for(l=lstart;l <= seq_length; l++){
          unsigned int ln2;
          int type2;
          kl = my_iindx[k]-l;
          if(!Q_B[kl]) continue;
          ln2 = (p - 1) + (seq_length - l);

          if((ln1+ln2) > MAXLOOP) continue;

          type2 = ptype[kl];
          if(!type2) continue;

          /* get distance to reference if closing the interior loop
          *  d2a = dbp(T1_[1,n}, T1_{p,q} + T1_{k,l})
          *  d2b = dbp(T2_[1,n}, T2_{p,q} + T2_{k,l})
          */
          da = base_d1 - referenceBPs1[pq] - referenceBPs1[kl];
          db = base_d2 - referenceBPs2[pq] - referenceBPs2[kl];
          qot = exp_E_IntLoop(ln2, ln1, rtype[type2], type, S1[l+1], S1[k-1], S1[p-1], S1[q+1], pf_params) * scale[ln1+ln2];

          for(cnt1 = k_min_values_b[pq]; cnt1 <= k_max_values_b[pq]; cnt1++)
            for(cnt2 = l_min_values_b[pq][cnt1]; cnt2 <= l_max_values_b[pq][cnt1]; cnt2+=2)
              for(cnt3 = k_min_values_b[kl]; cnt3 <= k_max_values_b[kl]; cnt3++)
                for(cnt4 = l_min_values_b[kl][cnt3]; cnt4 <= l_max_values_b[kl][cnt3]; cnt4+=2){
                  vars->Q_cI[cnt1 + cnt3 + da][(cnt2 + cnt4 + db)/2] += Q_B[pq][cnt1][cnt2/2] * Q_B[kl][cnt3][cnt4/2] * qot;
                  if(update_cI){
                    updatePosteriorBoundaries(cnt1 + cnt3 + da,
                                              cnt2 + cnt4 + db,
                                              &min_k_real_qcI,
                                              &max_k_real_qcI,
                                              &min_l_real_qcI,
                                              &max_l_real_qcI
                                              );
                  }
                }
        }
      }
    }

  if(update_cH){
    adjustArrayBoundaries(&vars->Q_cH,
                          &vars->k_min_values_qcH,
                          &vars->k_max_values_qcH,
                          &vars->l_min_values_qcH,
                          &vars->l_max_values_qcH,
                          min_k_real_qcH,
                          max_k_real_qcH,
                          min_l_real_qcH,
                          max_l_real_qcH
                        );
  }
  if(update_cI){
    adjustArrayBoundaries(&vars->Q_cI,
                          &vars->k_min_values_qcI,
                          &vars->k_max_values_qcI,
                          &vars->l_min_values_qcI,
                          &vars->l_max_values_qcI,
                          min_k_real_qcI,
                          max_k_real_qcI,
                          min_l_real_qcI,
                          max_l_real_qcI
                        );
  }

  /* 3. Multiloops  */
  if(seq_length > 2*TURN-3)
#ifdef _OPENMP
  #pragma omp parallel for private(k, da, db, cnt1, cnt2, cnt3, cnt4)
#endif
    for(k=TURN+2; k<seq_length-2*TURN-3; k++){
      /* get distancies to references
      * d3a = dbp(T1_[1,n}, T1_{1,k} + T1_{k+1, n})
      * d3b = dbp(T2_[1,n}, T2_{1,k} + T2_{k+1, n})
      */
      da = base_d1 - referenceBPs1[my_iindx[1]-k] - referenceBPs1[my_iindx[k+1]-seq_length];
      db = base_d2 - referenceBPs2[my_iindx[1]-k] - referenceBPs2[my_iindx[k+1]-seq_length];
      if(Q_M[my_iindx[1]-k] && vars->Q_M2[k+1])
        for(cnt1 = k_min_values_m[my_iindx[1]-k]; cnt1 <= k_max_values_m[my_iindx[1]-k]; cnt1++)
          for(cnt2 = l_min_values_m[my_iindx[1]-k][cnt1]; cnt2 <= l_max_values_m[my_iindx[1]-k][cnt1]; cnt2+=2)
            for(cnt3 = vars->k_min_values_m2[k+1]; cnt3 <= vars->k_max_values_m2[k+1]; cnt3++)
              for(cnt4 = vars->l_min_values_m2[k+1][cnt3]; cnt4 <= vars->l_max_values_m2[k+1][cnt3]; cnt4+=2){
                vars->Q_cM[cnt1 + cnt3 + da][(cnt2 + cnt4 + db)/2] += Q_M[my_iindx[1]-k][cnt1][cnt2/2] * vars->Q_M2[k+1][cnt3][cnt4/2] * pf_params->expMLclosing;
                if(update_cM){
                  updatePosteriorBoundaries(cnt1 + cnt3 + da,
                                            cnt2 + cnt4 + db,
                                            &min_k_real_qcM,
                                            &max_k_real_qcM,
                                            &min_l_real_qcM,
                                            &max_l_real_qcM
                                            );
                }
              }
    }
  if(update_cM){
    adjustArrayBoundaries(&vars->Q_cM,
                          &vars->k_min_values_qcM,
                          &vars->k_max_values_qcM,
                          &vars->l_min_values_qcM,
                          &vars->l_max_values_qcM,
                          min_k_real_qcM,
                          max_k_real_qcM,
                          min_l_real_qcM,
                          max_l_real_qcM
                        );
  }

  for(cnt1 = vars->k_min_values_qcH; cnt1 <= vars->k_max_values_qcH; cnt1++)
    for(cnt2 = vars->l_min_values_qcH[cnt1]; cnt2 <= vars->l_max_values_qcH[cnt1]; cnt2 += 2){
      vars->Q_c[cnt1][cnt2/2] += vars->Q_cH[cnt1][cnt2/2];
      if(update_c){
        updatePosteriorBoundaries(cnt1,
                                  cnt2,
                                  &min_k_real,
                                  &max_k_real,
                                  &min_l_real,
                                  &max_l_real
                                  );
      }
    }
  for(cnt1 = vars->k_min_values_qcI; cnt1 <= vars->k_max_values_qcI; cnt1++)
    for(cnt2 = vars->l_min_values_qcI[cnt1]; cnt2 <= vars->l_max_values_qcI[cnt1]; cnt2 += 2){
      vars->Q_c[cnt1][cnt2/2] += vars->Q_cI[cnt1][cnt2/2];
      if(update_c){
        updatePosteriorBoundaries(cnt1,
                                  cnt2,
                                  &min_k_real,
                                  &max_k_real,
                                  &min_l_real,
                                  &max_l_real
                                  );
      }
    }
  for(cnt1 = vars->k_min_values_qcM; cnt1 <= vars->k_max_values_qcM; cnt1++)
    for(cnt2 = vars->l_min_values_qcM[cnt1]; cnt2 <= vars->l_max_values_qcM[cnt1]; cnt2 += 2){
      vars->Q_c[cnt1][cnt2/2] += vars->Q_cM[cnt1][cnt2/2];
      if(update_c){
        updatePosteriorBoundaries(cnt1,
                                  cnt2,
                                  &min_k_real,
                                  &max_k_real,
                                  &min_l_real,
                                  &max_l_real
                                  );
      }
    }
  /* add the case were structure is unfolded chain */
  vars->Q_c[referenceBPs1[my_iindx[1]-seq_length]][referenceBPs2[my_iindx[1]-seq_length]/2] += 1.0 * scale[seq_length];
  if(update_c){
    updatePosteriorBoundaries(referenceBPs1[my_iindx[1]-seq_length],
                              referenceBPs2[my_iindx[1]-seq_length],
                              &min_k_real,
                              &max_k_real,
                              &min_l_real,
                              &max_l_real
                              );
  }

  adjustArrayBoundaries(&vars->Q_c,
                        &vars->k_min_values_qc,
                        &vars->k_max_values_qc,
                        &vars->l_min_values_qc,
                        &vars->l_max_values_qc,
                        min_k_real,
                        max_k_real,
                        min_l_real,
                        max_l_real
                      );
}


PRIVATE void make_ptypes2(TwoDpfold_vars *vars) {
  int n,i,j,k,l;

  n=vars->S[0];
  for (k=1; k<n-TURN; k++)
    for (l=1; l<=2; l++) {
      int type,ntype=0,otype=0;
      i=k; j = i+TURN+l; if (j>n) continue;
      type = pair[vars->S[i]][vars->S[j]];
      while ((i>=1)&&(j<=n)) {
        if ((i>1)&&(j<n)) ntype = pair[vars->S[i-1]][vars->S[j+1]];
        if (noLonelyPairs && (!otype) && (!ntype))
          type = 0; /* i.j can only form isolated pairs */
        vars->ptype[vars->my_iindx[i]-j] = (char) type;
        otype =  type;
        type  = ntype;
        i--; j++;
      }
    }
}

PRIVATE void scale_pf_params2(TwoDpfold_vars *vars)
{
  /* scale energy parameters and pre-calculate Boltzmann weights */
  unsigned int i;
  double  kT;

  vars->pf_params = get_scaled_pf_parameters();
  vars->init_temp = vars->pf_params->temperature;
  vars->temperature = vars->init_temp;

  kT = vars->pf_params->kT;   /* kT in cal/mol  */

   /* scaling factors (to avoid overflows) */
  if (vars->pf_scale == -1) { /* mean energy for random sequences: 184.3*length cal */
    vars->pf_scale = exp(-(-185+(vars->temperature-37.)*7.27)/kT);
    if (vars->pf_scale<1) vars->pf_scale=1;
  }
  vars->scale[0] = 1.;   vars->scale[1] = 1./vars->pf_scale;
  for (i=2; i<=vars->seq_length; i++) {
    vars->scale[i] = vars->scale[i/2]*vars->scale[i-(i/2)];
  }
}

/*
* ###################################################
* stochastic backtracking
* ###################################################
*/

PUBLIC char *TwoDpfold_pbacktrack(TwoDpfold_vars *vars, unsigned int d1, unsigned int d2){

  FLT_OR_DBL      r, qt, *qln, *scale;
  unsigned int    i, j, n, start, maxD1, maxD2, base_d1, base_d2, da, db;
  unsigned int    *referenceBPs1, *referenceBPs2;
  int             *my_iindx, ij, cnt1, cnt2, cnt3, cnt4;
  pf_paramT       *pf_params;     /* holds all [unscaled] pf parameters */
  char            *pstruc, *ptype;
  short           *S1, *reference_pt1, *reference_pt2;
  FLT_OR_DBL      ***Q, ***Q_B;
  int             **l_min_values, **l_max_values,**l_min_values_b, **l_max_values_b;
  int             *k_min_values, *k_max_values,*k_min_values_b, *k_max_values_b;

  pf_params       = vars->pf_params;
  n               = vars->seq_length;
  maxD1           = vars->maxD1;
  maxD2           = vars->maxD2;
  my_iindx        = vars->my_iindx;
  scale           = vars->scale;
  ptype           = vars->ptype;
  S1              = vars->S1;
  reference_pt1   = vars->reference_pt1;
  reference_pt2   = vars->reference_pt2;
  referenceBPs1   = vars->referenceBPs1;
  referenceBPs2   = vars->referenceBPs2;

  Q               = vars->Q;
  l_min_values    = vars->l_min_values;
  l_max_values    = vars->l_max_values;
  k_min_values    = vars->k_min_values;
  k_max_values    = vars->k_max_values;

  Q_B             = vars->Q_B;
  l_min_values_b  = vars->l_min_values_b;
  l_max_values_b  = vars->l_max_values_b;
  k_min_values_b  = vars->k_min_values_b;
  k_max_values_b  = vars->k_max_values_b;

  if(d1 > maxD1)
    nrerror("pbacktrack@2Dpfold.c: distance to 1st reference structure to high!");
  if(d2 > maxD2)
    nrerror("pbacktrack@2Dpfold.c: distance to 2nd reference structure to high!");

  pstruc = space((n+1)*sizeof(char));

  for (i=0; i<n; i++) pstruc[i] = '.';
  pstruc[i] = '\0';

  start = 1;
  while (start<n) {
    int sn = my_iindx[start] - n;
    /* find i position of first pair */
    FLT_OR_DBL qln_i = 0, qln_i1 = 0;

    /* get partition function of current subsection [sn,n] */
    if((d1 >= k_min_values[sn]) && (d1 <= k_max_values[sn]))
      if((d2 >= l_min_values[sn][d1]) && (d2 <= l_max_values[sn][d1]))
        qln_i = Q[sn][d1][d2/2];


    for (i=start; i<n; i++) {
      r = urn() * qln_i;
      da = referenceBPs1[sn] - referenceBPs1[my_iindx[i+1] - n];
      db = referenceBPs2[sn] - referenceBPs2[my_iindx[i+1] - n];
      qln_i1 = 0;
      if(d1 >= da && d2 >= db)
        if((d1-da >= k_min_values[my_iindx[i+1] - n]) && (d1 - da <= k_max_values[my_iindx[i+1] - n]))
          if((d2 - db >= l_min_values[my_iindx[i+1] - n][d1-da]) && (d2 - db <= l_max_values[my_iindx[i+1] - n][d1 - da]))
            qln_i1 = Q[my_iindx[i+1] - n][d1-da][(d2-db)/2];
      if (r > qln_i1*scale[1])  break; /* i is paired */
      qln_i = qln_i1;
    }

    if (i>=n) break; /* no more pairs */

    /* now find the pairing partner j */
    r = urn() * (qln_i - qln_i1*scale[1]);

    for (qt=0, j=i+1; j<n; j++) {
      int type;
      type = ptype[my_iindx[i]-j];
      if (type) {
        double qkl = 1.0;
        ij = my_iindx[i]-j;
        qkl *= exp_E_ExtLoop(type, (i>1) ? S1[i-1] : -1, S1[j+1], pf_params);

        da = referenceBPs1[sn] - referenceBPs1[ij] - referenceBPs1[my_iindx[j+1]-n];
        db = referenceBPs2[sn] - referenceBPs2[ij] - referenceBPs2[my_iindx[j+1]-n];

        if((d1 >= da) && (d2 >= db) && Q_B[ij] && Q[my_iindx[j+1]-n])
          for(cnt1 = k_min_values_b[ij]; cnt1 <= MIN2(k_max_values_b[ij], d1-da); cnt1++)
            for(cnt2 = l_min_values_b[ij][cnt1]; cnt2 <= MIN2(l_max_values_b[ij][cnt1], d2-db); cnt2+=2)
              if((d1-da-cnt1 >= k_min_values[my_iindx[j+1]-n]) && (d1-da-cnt1 <= k_max_values[my_iindx[j+1]-n]))
                if((d2-db-cnt2 >= l_min_values[my_iindx[j+1]-n][d1-da-cnt1]) && (d2 - db - cnt2 <= l_max_values[my_iindx[j+1]-n][d1-da-cnt1])){
                  qt += qkl * Q_B[ij][cnt1][cnt2/2] * Q[my_iindx[j+1]-n][d1-da-cnt1][(d2-db-cnt2)/2];
                  if(qt >= r)
                    goto pbacktrack_ext_loop_early_escape;
                }
      }
    }
    /* now dont forget the case j==n */
    j = n;
    ij = my_iindx[i]-j;
    int type = ptype[ij];
    if (type) {
      double qkl = 1.0;

      qkl *= exp_E_ExtLoop(type, (i>1) ? S1[i-1] : -1, -1, pf_params);

      da = referenceBPs1[sn] - referenceBPs1[ij];
      db = referenceBPs2[sn] - referenceBPs2[ij];
      if(d1 >= da && d2 >= db){
        cnt1 = d1 - da;
        cnt2 = d2 - db;
        if((cnt1 >= k_min_values_b[ij]) && (cnt1 <= k_max_values_b[ij]))
          if((cnt2 >= l_min_values_b[ij][cnt1]) && (cnt2 <= l_max_values_b[ij][cnt1])){
            qt += qkl * Q_B[ij][cnt1][cnt2/2];
            if(qt >= r)
              goto pbacktrack_ext_loop_early_escape; /* j is paired */
          }
      }
    }
    j++;

pbacktrack_ext_loop_early_escape:

    if (j==n+1){
      nrerror("pbacktrack@2Dpfold.c: backtracking failed in ext loop");
    }

    backtrack(vars, pstruc, cnt1, cnt2, i,j);

    if(j==n) break;
    start = j+1;
    d1 -= cnt1 + da;
    d2 -= cnt2 + db;
  }
  return pstruc;
}

PUBLIC char *TwoDpfold_pbacktrack_f5(TwoDpfold_vars *vars, unsigned int d1, unsigned int d2, unsigned int length){

  FLT_OR_DBL      r, qt, *qln, *scale;
  unsigned int    i, j, n, start, maxD1, maxD2, base_d1, base_d2, da, db;
  unsigned int    *referenceBPs1, *referenceBPs2;
  int             *my_iindx, ij, cnt1, cnt2, cnt3, cnt4;
  pf_paramT       *pf_params;     /* holds all [unscaled] pf parameters */
  char            *pstruc, *ptype;
  short           *S1, *reference_pt1, *reference_pt2;
  FLT_OR_DBL      ***Q, ***Q_B;
  int             **l_min_values, **l_max_values,**l_min_values_b, **l_max_values_b;
  int             *k_min_values, *k_max_values,*k_min_values_b, *k_max_values_b;

  n               = vars->seq_length;
  pf_params       = vars->pf_params;
  maxD1           = vars->maxD1;
  maxD2           = vars->maxD2;
  my_iindx        = vars->my_iindx;
  scale           = vars->scale;
  ptype           = vars->ptype;
  S1              = vars->S1;
  reference_pt1   = vars->reference_pt1;
  reference_pt2   = vars->reference_pt2;
  referenceBPs1   = vars->referenceBPs1;
  referenceBPs2   = vars->referenceBPs2;

  Q               = vars->Q;
  l_min_values    = vars->l_min_values;
  l_max_values    = vars->l_max_values;
  k_min_values    = vars->k_min_values;
  k_max_values    = vars->k_max_values;

  Q_B             = vars->Q_B;
  l_min_values_b  = vars->l_min_values_b;
  l_max_values_b  = vars->l_max_values_b;
  k_min_values_b  = vars->k_min_values_b;
  k_max_values_b  = vars->k_max_values_b;

  if(d1 > maxD1)
    nrerror("pbacktrack@2Dpfold.c: distance to 1st reference structure to high!");
  if(d2 > maxD2)
    nrerror("pbacktrack@2Dpfold.c: distance to 2nd reference structure to high!");

  pstruc = space((length+1)*sizeof(char));

  for (i=0; i<length; i++) pstruc[i] = '.';
  pstruc[i] = '\0';

  start = 1;
  while (start<length) {
    int sn = my_iindx[start] - length;
    /* find i position of first pair */
    FLT_OR_DBL qln_i = 0, qln_i1 = 0;

    /* get partition function of current subsection [start, length] */
    if((d1 >= k_min_values[sn]) && (d1 <= k_max_values[sn]))
      if((d2 >= l_min_values[sn][d1]) && (d2 <= l_max_values[sn][d1]))
        qln_i = Q[sn][d1][d2/2];


    for (i=start; i<length; i++) {
      r = urn() * qln_i;
      da = referenceBPs1[sn] - referenceBPs1[my_iindx[i+1] - length];
      db = referenceBPs2[sn] - referenceBPs2[my_iindx[i+1] - length];
      qln_i1 = 0;
      if(d1 >= da && d2 >= db)
        if((d1-da >= k_min_values[my_iindx[i+1] - length]) && (d1 - da <= k_max_values[my_iindx[i+1] - length]))
          if((d2 - db >= l_min_values[my_iindx[i+1] - length][d1-da]) && (d2 - db <= l_max_values[my_iindx[i+1] - length][d1 - da]))
            qln_i1 = Q[my_iindx[i+1] - length][d1-da][(d2-db)/2];
      if (r > qln_i1*scale[1])  break; /* i is paired */
      qln_i = qln_i1;
    }

    if (i>=length) break; /* no more pairs */

    /* now find the pairing partner j */
    r = urn() * (qln_i - qln_i1*scale[1]);

    for (qt=0, j=i+1; j<length; j++) {
      int type;
      type = ptype[my_iindx[i]-j];
      if (type) {
        double qkl = 1.0;
        ij = my_iindx[i]-j;
        qkl *= exp_E_ExtLoop(type, (i>1) ? S1[i-1] : -1, S1[j+1], pf_params);

        da = referenceBPs1[sn] - referenceBPs1[ij] - referenceBPs1[my_iindx[j+1]-length];
        db = referenceBPs2[sn] - referenceBPs2[ij] - referenceBPs2[my_iindx[j+1]-length];

        if((d1 >= da) && (d2 >= db) && Q_B[ij] && Q[my_iindx[j+1]-length])
          for(cnt1 = k_min_values_b[ij]; cnt1 <= MIN2(k_max_values_b[ij], d1-da); cnt1++)
            for(cnt2 = l_min_values_b[ij][cnt1]; cnt2 <= MIN2(l_max_values_b[ij][cnt1], d2-db); cnt2+=2)
              if((d1-da-cnt1 >= k_min_values[my_iindx[j+1]-length]) && (d1-da-cnt1 <= k_max_values[my_iindx[j+1]-length]))
                if((d2-db-cnt2 >= l_min_values[my_iindx[j+1]-length][d1-da-cnt1]) && (d2 - db - cnt2 <= l_max_values[my_iindx[j+1]-length][d1-da-cnt1])){
                  qt += qkl * Q_B[ij][cnt1][cnt2/2] * Q[my_iindx[j+1]-length][d1-da-cnt1][(d2-db-cnt2)/2];
                  if(qt >= r)
                    goto pbacktrack_ext_loop_early_escape2;
                }
      }
    }
    /* now dont forget the case j==length */
    j = length;
    ij = my_iindx[i]-j;
    int type = ptype[ij];
    if (type) {
      double qkl = exp_E_ExtLoop(type, i>1 ? S1[i-1] : -1, length < n ? S1[j+1] : -1, pf_params);

      da = referenceBPs1[sn] - referenceBPs1[ij];
      db = referenceBPs2[sn] - referenceBPs2[ij];

      if(d1 >= da && d2 >= db){
        cnt1 = d1 - da;
        cnt2 = d2 - db;
        if((cnt1 >= k_min_values_b[ij]) && (cnt1 <= k_max_values_b[ij]))
          if((cnt2 >= l_min_values_b[ij][cnt1]) && (cnt2 <= l_max_values_b[ij][cnt1])){
            qt += qkl * Q_B[ij][cnt1][cnt2/2];
            if(qt >= r)
              goto pbacktrack_ext_loop_early_escape2; /* j is paired */
          }
      }
    }
    j++;

pbacktrack_ext_loop_early_escape2:

    if (j==length+1){
      nrerror("pbacktrack@2Dpfold.c: backtracking failed in ext loop");
    }

    backtrack(vars, pstruc, cnt1, cnt2, i,j);

    if(j==length) break;
    start = j+1;
    d1 -= cnt1 + da;
    d2 -= cnt2 + db;
  }
  return pstruc;
}

PRIVATE void backtrack(TwoDpfold_vars *vars, char *pstruc, unsigned int d1, unsigned int d2, unsigned int i, unsigned int j) {
  FLT_OR_DBL      r, qt, *qln, *scale;
  unsigned int    n, start, maxD1, maxD2, base_d1, base_d2, da, db, remaining_d1, remaining_d2;
  unsigned int    *referenceBPs1, *referenceBPs2;
  pf_paramT       *pf_params;     /* holds all [unscaled] pf parameters */
  char            *ptype, *sequence;
  short           *S1, *reference_pt1, *reference_pt2;
  int             *my_iindx, *jindx, ij, cnt1, cnt2, cnt3, cnt4;

  pf_params   = vars->pf_params;
  sequence    = vars->sequence;
  n           = vars->seq_length;
  maxD1       = vars->maxD1;
  maxD2       = vars->maxD2;
  my_iindx    = vars->my_iindx;
  jindx       = vars->jindx;
  scale       = vars->scale;
  ptype       = vars->ptype;
  S1              = vars->S1;
  reference_pt1   = vars->reference_pt1;
  reference_pt2   = vars->reference_pt2;
  referenceBPs1   = vars->referenceBPs1;
  referenceBPs2   = vars->referenceBPs2;

  FLT_OR_DBL  ***Q, ***Q_B, ***Q_M, ***Q_M1;
  int         *k_min_values, *k_max_values, *k_min_values_m, *k_max_values_m,*k_min_values_m1, *k_max_values_m1,*k_min_values_b, *k_max_values_b;
  int         **l_min_values, **l_max_values, **l_min_values_m, **l_max_values_m,**l_min_values_m1, **l_max_values_m1,**l_min_values_b, **l_max_values_b;

  Q = vars->Q;
  k_min_values = vars->k_min_values;
  k_max_values = vars->k_max_values;
  l_min_values = vars->l_min_values;
  l_max_values = vars->l_max_values;

  Q_B = vars->Q_B;
  k_min_values_b = vars->k_min_values_b;
  k_max_values_b = vars->k_max_values_b;
  l_min_values_b = vars->l_min_values_b;
  l_max_values_b = vars->l_max_values_b;

  Q_M = vars->Q_M;
  k_min_values_m = vars->k_min_values_m;
  k_max_values_m = vars->k_max_values_m;
  l_min_values_m = vars->l_min_values_m;
  l_max_values_m = vars->l_max_values_m;

  Q_M1 = vars->Q_M1;
  k_min_values_m1 = vars->k_min_values_m1;
  k_max_values_m1 = vars->k_max_values_m1;
  l_min_values_m1 = vars->l_min_values_m1;
  l_max_values_m1 = vars->l_max_values_m1;

  do {
    double r, qbt1 = 0.;
    unsigned int k, l, u, u1;
    int type;

    pstruc[i-1] = '('; pstruc[j-1] = ')';

    r = 0.;
    ij = my_iindx[i]-j;

    if((d1 >= k_min_values_b[ij]) && (d1 <= k_max_values_b[ij]))
      if((d2 >= l_min_values_b[ij][d1]) && (d2 <= l_max_values_b[ij][d1]))
        r = urn() * Q_B[ij][d1][d2/2];

    if(r == 0.) nrerror("backtrack@2Dpfold.c: backtracking failed\n");

    type = ptype[ij];
    u = j-i-1;
    base_d1 = ((unsigned int)reference_pt1[i] != j) ? 1 : -1;
    base_d2 = ((unsigned int)reference_pt2[i] != j) ? 1 : -1;

    da = base_d1 + referenceBPs1[ij];
    db = base_d2 + referenceBPs2[ij];

    /*hairpin contribution*/
    if((da == d1) && (db == d2))
      if(!(((type==3)||(type==4))&&no_closingGU))
        qbt1 = exp_E_Hairpin(u, type, S1[i+1], S1[j-1], sequence+i-1, pf_params) * scale[u+2];

    if (qbt1>=r) return; /* found the hairpin we're done */

    for (k=i+1; k<=MIN2(i+MAXLOOP+1,j-TURN-2); k++) {
      unsigned int u_pre, lmin;
      u1 = k-i-1;
      lmin = k + TURN + 1;
      u_pre = u1 + j;
      /* lmin = MAX2(k + TURN + 1, u1 + j - 1 - MAXLOOP) */
      if(u_pre > lmin + MAXLOOP) lmin = u_pre - 1 - MAXLOOP;
      for (l=lmin; l<j; l++) {
        int type_2;
        type_2 = ptype[my_iindx[k]-l];
        if (type_2) {
          da = base_d1 + referenceBPs1[my_iindx[i]-j] - referenceBPs1[my_iindx[k]-l];
          db = base_d2 + referenceBPs2[my_iindx[i]-j] - referenceBPs2[my_iindx[k]-l];
          type_2 = rtype[type_2];
          FLT_OR_DBL tmp_en = exp_E_IntLoop(u1, j-l-1, type, type_2, S1[i+1], S1[j-1], S1[k-1], S1[l+1], pf_params) * scale[u1+j-l+1];
          if(d1 >= da && d2 >= db)
            if((d1 - da >= k_min_values_b[my_iindx[k]-l]) && (d1 - da <= k_max_values_b[my_iindx[k]-l]))
              if((d2 - db >= l_min_values_b[my_iindx[k]-l][d1 - da]) && (d2 - db <= l_max_values_b[my_iindx[k]-l][d1 - da])){
                cnt1 = d1 - da;
                cnt2 = d2 - db;
                qbt1 += Q_B[my_iindx[k]-l][cnt1][cnt2/2] * tmp_en;
                if(qbt1 > r) goto backtrack_int_early_escape;
              }
        }
      }
    }

backtrack_int_early_escape:
    if (l<j) {
      i=k; j=l;
      d1 = cnt1;
      d2 = cnt2;
    }
    else break;
  } while (1);

  /* backtrack in multi-loop */
  {
    double r, qt;
    unsigned int k, ii, jj;

    base_d1 = ((unsigned int)reference_pt1[i] != j) ? 1 : -1;
    base_d2 = ((unsigned int)reference_pt2[i] != j) ? 1 : -1;

    base_d1 += referenceBPs1[my_iindx[i]-j];
    base_d2 += referenceBPs2[my_iindx[i]-j];

    i++; j--;
    /* find the first split index */
    ii = my_iindx[i]; /* ii-j=[i,j] */
    jj = jindx[j]; /* jj+i=[j,i] */
    /* get total contribution */
    for (qt=0., k=i+1; k<j; k++){
      /* calculate introduced distance to reference structures */
      da = base_d1 - referenceBPs1[my_iindx[i]-k+1] - referenceBPs1[my_iindx[k]-j];
      db = base_d2 - referenceBPs2[my_iindx[i]-k+1] - referenceBPs2[my_iindx[k]-j];
      /* collect all contributing energies */
      if(d1 >= da && d2 >= db && Q_M[ii-k+1] && Q_M1[jj+k])
        for(cnt1 = k_min_values_m[ii-k+1]; cnt1 <= MIN2(k_max_values_m[ii-k+1], d1-da); cnt1++)
          for(cnt2 = l_min_values_m[ii-k+1][cnt1]; cnt2 <= MIN2(l_max_values_m[ii-k+1][cnt1], d2 - db); cnt2+=2)
            if((d1-cnt1-da >= k_min_values_m1[jj+k]) && (d1-cnt1-da <= k_max_values_m1[jj+k]))
              if((d2 - cnt2 - db >= l_min_values_m1[jj+k][d1-da-cnt1]) && (d2 - cnt2 - db <= l_max_values_m1[jj+k][d1-cnt1-da]))
                qt += Q_M[ii-k+1][cnt1][cnt2/2] * Q_M1[jj+k][d1-da-cnt1][(d2-db-cnt2)/2];
    }
    r = urn() * qt;
    for (qt=0., k=i+1; k<j; k++) {
      /* calculate introduced distance to reference structures */
      da = base_d1 - referenceBPs1[my_iindx[i]-k+1] - referenceBPs1[my_iindx[k]-j];
      db = base_d2 - referenceBPs2[my_iindx[i]-k+1] - referenceBPs2[my_iindx[k]-j];
      /* collect all contributing energies */
      if(d1 >= da && d2 >= db && Q_M[ii-k+1] && Q_M1[jj+k])
        for(cnt1 = k_min_values_m[ii-k+1]; cnt1 <= MIN2(k_max_values_m[ii-k+1], d1-da); cnt1++)
          for(cnt2 = l_min_values_m[ii-k+1][cnt1]; cnt2 <= MIN2(l_max_values_m[ii-k+1][cnt1], d2 - db); cnt2+=2)
            if((d1-cnt1-da >= k_min_values_m1[jj+k]) && (d1-cnt1-da <= k_max_values_m1[jj+k]))
              if((d2 - cnt2 - db >= l_min_values_m1[jj+k][d1-da-cnt1]) && (d2 - cnt2 - db <= l_max_values_m1[jj+k][d1-cnt1-da])){
                cnt3 = d1-da-cnt1;
                cnt4 = d2-db-cnt2;
                qt += Q_M[ii-k+1][cnt1][cnt2/2] * Q_M1[jj+k][cnt3][cnt4/2];
                if (qt>=r) goto backtrack_ml_early_escape;
              }
    }
    if (k>=j) nrerror("backtrack failed, can't find split index ");

backtrack_ml_early_escape:

    backtrack_qm1(vars, pstruc, cnt3, cnt4, k, j);

    j = k-1;
    backtrack_qm(vars, pstruc, cnt1, cnt2, i, j);
  }
}

PRIVATE void backtrack_qm1(TwoDpfold_vars *vars, char *pstruc, unsigned int d1, unsigned int d2, unsigned int i, unsigned int j){
  /* i is paired to l, i<l<j; backtrack in qm1 to find l */
  FLT_OR_DBL      r, qt, *qln, *scale;
  unsigned int    n, start, maxD1, maxD2, base_d1, base_d2, da, db, remaining_d1, remaining_d2;
  unsigned int    *referenceBPs1, *referenceBPs2;
  pf_paramT       *pf_params;     /* holds all [unscaled] pf parameters */
  char            *ptype, *sequence;
  short           *S1, *reference_pt1, *reference_pt2;
  int             *my_iindx, *jindx, cnt1, cnt2;
  pf_params   = vars->pf_params;
  sequence    = vars->sequence;
  n           = vars->seq_length;
  maxD1       = vars->maxD1;
  maxD2       = vars->maxD2;
  my_iindx    = vars->my_iindx;
  jindx       = vars->jindx;
  scale       = vars->scale;
  ptype       = vars->ptype;
  S1          = vars->S1;
  referenceBPs1  = vars->referenceBPs1;
  referenceBPs2  = vars->referenceBPs2;
  reference_pt1  = vars->reference_pt1;
  reference_pt2  = vars->reference_pt2;

  FLT_OR_DBL  ***Q_B, ***Q_M1;
  int         *k_min_values_m1, *k_max_values_m1,*k_min_values_b, *k_max_values_b;
  int         **l_min_values_m1, **l_max_values_m1,**l_min_values_b, **l_max_values_b;

  Q_B = vars->Q_B;
  k_min_values_b = vars->k_min_values_b;
  k_max_values_b = vars->k_max_values_b;
  l_min_values_b = vars->l_min_values_b;
  l_max_values_b = vars->l_max_values_b;

  Q_M1 = vars->Q_M1;
  k_min_values_m1 = vars->k_min_values_m1;
  k_max_values_m1 = vars->k_max_values_m1;
  l_min_values_m1 = vars->l_min_values_m1;
  l_max_values_m1 = vars->l_max_values_m1;

  unsigned int ii, l;
  int type;

  /* find qm1 contribution */
  if((d1 >= k_min_values_m1[jindx[j]+i]) && (d1 <= k_max_values_m1[jindx[j]+i]))
    if((d2 >= l_min_values_m1[jindx[j]+i][d1]) && (d2 <= l_max_values_m1[jindx[j]+i][d1]))
      r = urn() * Q_M1[jindx[j]+i][d1][d2/2];
  if(r == 0.) nrerror("backtrack_qm1@2Dpfold.c: backtracking failed\n");

  ii = my_iindx[i];
  for (qt=0., l=i+TURN+1; l<=j; l++) {
    type = ptype[ii-l];
    if (type){
      /* compute the introduced distance to reference structures */
      da = referenceBPs1[my_iindx[i]-j] - referenceBPs1[my_iindx[i]-l];
      db = referenceBPs2[my_iindx[i]-j] - referenceBPs2[my_iindx[i]-l];

      FLT_OR_DBL tmp = exp_E_MLstem(type, S1[i-1], S1[l+1], pf_params) * pow(pf_params->expMLbase, j-l) * scale[j-l];

      /* get energy contributions */
      if(d1 >= da && d2 >= db)
        if((d1 - da >= k_min_values_b[ii-l]) && (d1 - da <= k_max_values_b[ii-l]))
          if((d2 - db >= l_min_values_b[ii-l][d1-da]) && (d2 - db <= l_max_values_b[ii-l][d1-da])){
            cnt1 = d1 - da;
            cnt2 = d2 - db;
            qt += Q_B[ii-l][cnt1][cnt2/2] * tmp;
            if (qt>=r) goto backtrack_qm1_early_escape;
          }
    }
  }
  if (l>j) nrerror("backtrack failed in qm1");
backtrack_qm1_early_escape:

  backtrack(vars, pstruc, cnt1, cnt2, i,l);
}

PRIVATE void backtrack_qm(TwoDpfold_vars *vars, char *pstruc, unsigned int d1, unsigned int d2, unsigned int i, unsigned int j){
  /* divide multiloop into qm and qm1  */
  FLT_OR_DBL      r, qt, *qln, *scale;
  unsigned int    n, start, maxD1, maxD2, base_d1, base_d2, da, db, da2, db2, remaining_d1, remaining_d2;
  unsigned int    *referenceBPs1, *referenceBPs2;
  pf_paramT       *pf_params;     /* holds all [unscaled] pf parameters */
  char            *ptype, *sequence;
  short           *S1, *reference_pt1, *reference_pt2;
  int             *my_iindx, *jindx, cnt1, cnt2, cnt3, cnt4;

  pf_params   = vars->pf_params;
  sequence    = vars->sequence;
  n           = vars->seq_length;
  maxD1       = vars->maxD1;
  maxD2       = vars->maxD2;
  my_iindx    = vars->my_iindx;
  jindx       = vars->jindx;
  scale       = vars->scale;
  ptype       = vars->ptype;
  S1          = vars->S1;
  referenceBPs1  = vars->referenceBPs1;
  referenceBPs2  = vars->referenceBPs2;

  FLT_OR_DBL  ***Q_M, ***Q_M1;
  int         *k_min_values_m, *k_max_values_m,*k_min_values_m1, *k_max_values_m1;
  int         **l_min_values_m, **l_max_values_m,**l_min_values_m1, **l_max_values_m1;

  Q_M = vars->Q_M;
  k_min_values_m = vars->k_min_values_m;
  k_max_values_m = vars->k_max_values_m;
  l_min_values_m = vars->l_min_values_m;
  l_max_values_m = vars->l_max_values_m;

  Q_M1 = vars->Q_M1;
  k_min_values_m1 = vars->k_min_values_m1;
  k_max_values_m1 = vars->k_max_values_m1;
  l_min_values_m1 = vars->l_min_values_m1;
  l_max_values_m1 = vars->l_max_values_m1;

  double qmt = 0;
  unsigned int k;
  while(j>i){
    /* now backtrack  [i ... j] in qm[] */
    int sw;
    /* find qm contribution */
    if(Q_M[my_iindx[i]-j])
      if((d1 >= k_min_values_m[my_iindx[i]-j]) && (d1 <= k_max_values_m[my_iindx[i]-j]))
        if((d2 >= l_min_values_m[my_iindx[i]-j][d1]) && (d2 <= l_max_values_m[my_iindx[i]-j][d1]))
          r = urn() * Q_M[my_iindx[i]-j][d1][d2/2];
    if(r == 0.) nrerror("backtrack_qm@2Dpfold.c: backtracking failed in finding qm contribution\n");

    /* find corresponding qm1 contribution */
    qmt = 0.;
    if(Q_M1[jindx[j]+i])
      if((d1 >= k_min_values_m1[jindx[j]+i]) && (d1 <= k_max_values_m1[jindx[j]+i]))
        if((d2 >= l_min_values_m1[jindx[j]+i][d1]) && (d2 <= l_max_values_m1[jindx[j]+i][d1])){
          qmt = Q_M1[jindx[j]+i][d1][d2/2];
        }

    k=i;
    if(qmt<r){
      for(k=i+1; k<=j; k++){
        /* calculate introduced distancies to reference structures */
        da2 = referenceBPs1[my_iindx[i]-j] - referenceBPs1[my_iindx[k]-j];
        db2 = referenceBPs2[my_iindx[i]-j] - referenceBPs2[my_iindx[k]-j];
        da = da2 - referenceBPs1[my_iindx[i]-k+1];
        db = db2 - referenceBPs2[my_iindx[i]-k+1];


        FLT_OR_DBL tmp = pow(pf_params->expMLbase, k-i) * scale[k-i];

        /* collect unpaired + qm1 contributions */
        if(d1 >= da2 && d2 >= db2)
          if((d1 - da2 >= k_min_values_m1[jindx[j]+k]) && (d1 - da2 <= k_max_values_m1[jindx[j]+k]))
            if((d2 - db2 >= l_min_values_m1[jindx[j]+k][d1-da2]) && (d2 - db2 <= l_max_values_m1[jindx[j]+k][d1-da2])){
              cnt3 = d1-da2;
              cnt4 = d2-db2;
              qmt += Q_M1[jindx[j]+k][cnt3][cnt4/2] * tmp;
              if(qmt >= r){
                backtrack_qm1(vars, pstruc, cnt3, cnt4, k, j);
                return;
              }
            }

        /* collect qm + qm1 contributions */
        if(d1 >= da && d2 >= db && Q_M[my_iindx[i]-k+1] && Q_M1[jindx[j]+k])
          for(cnt1 = k_min_values_m[my_iindx[i]-k+1]; cnt1 <= MIN2(k_max_values_m[my_iindx[i]-k+1], d1 - da); cnt1++)
            for(cnt2 = l_min_values_m[my_iindx[i]-k+1][cnt1]; cnt2 <= MIN2(l_max_values_m[my_iindx[i]-k+1][cnt1], d2 - db); cnt2+=2)
              if((d1 - da - cnt1 >= k_min_values_m1[jindx[j]+k]) && (d1 - da - cnt1 <= k_max_values_m1[jindx[j]+k]))
                if((d2 - db - cnt2 >= l_min_values_m1[jindx[j]+k][d1-da-cnt1]) && (d2 - db - cnt2 <= l_max_values_m1[jindx[j]+k][d1-da-cnt1])){
                  cnt3 = d1 - da - cnt1;
                  cnt4 = d2 - db - cnt2;
                  qmt += Q_M[my_iindx[i]-k+1][cnt1][cnt2/2] * Q_M1[jindx[j]+k][cnt3][cnt4/2];
                  if(qmt >= r) goto backtrack_qm_early_escape;
                }
      }
    }
    else{
      backtrack_qm1(vars, pstruc, d1, d2, k, j);
      return;
    }

    if(k>j) nrerror("backtrack_qm@2Dpfold.c: backtrack failed in qm");

backtrack_qm_early_escape:

    backtrack_qm1(vars, pstruc, cnt3, cnt4, k, j);

    if(k<i+TURN) break; /* no more pairs */

    d1 = cnt1;
    d2 = cnt2;

    if(d1 == referenceBPs1[my_iindx[i]-k+1] && d2 == referenceBPs2[my_iindx[i]-k+1]){
      /* is interval [i,k] totally unpaired? */
      FLT_OR_DBL tmp = pow(pf_params->expMLbase, k-i) * scale[k-i];
      r = urn() * (Q_M[my_iindx[i]-k+1][d1][d2/2] + tmp);
      if(tmp >= r) return; /* no more pairs */
    }
    j = k-1;
  }
}




PRIVATE void adjustArrayBoundaries(FLT_OR_DBL ***array, int *k_min, int *k_max, int **l_min, int **l_max, int k_min_post, int k_max_post, int *l_min_post, int *l_max_post){
  int cnt1, cnt2;
  int k_diff_pre  = k_min_post - *k_min;
  int mem_size    = k_max_post - k_min_post + 1;

  if(k_min_post < INF){
    /* free all the unused memory behind actual data */
    for(cnt1 = k_max_post + 1; cnt1 <= *k_max; cnt1++){
      (*array)[cnt1] += (*l_min)[cnt1]/2;
      free((*array)[cnt1]);
    }

    /* free unused memory before actual data */
    for(cnt1 = *k_min; cnt1 < k_min_post; cnt1++){
      (*array)[cnt1] += (*l_min)[cnt1]/2;
      free((*array)[cnt1]);
    }
    /* move data to front and thereby eliminating unused memory in front of actual data */
    if(k_diff_pre > 0){
      memmove((FLT_OR_DBL **)(*array),((FLT_OR_DBL **)(*array)) + k_diff_pre, sizeof(FLT_OR_DBL *) * mem_size);
      memmove((int *) (*l_min),((int *) (*l_min)) + k_diff_pre, sizeof(int)   * mem_size);
      memmove((int *) (*l_max),((int *) (*l_max)) + k_diff_pre, sizeof(int)   * mem_size);
    }

    /* reallocating memory to actual size used */
    *array  += *k_min;
    *array = (FLT_OR_DBL **)realloc(*array, sizeof(FLT_OR_DBL *) * mem_size);
    *array -= k_min_post;

    *l_min  += *k_min;
    *l_min = (int *)realloc(*l_min, sizeof(int) * mem_size);
    *l_min -= k_min_post;

    *l_max  += *k_min;
    *l_max = (int *)realloc(*l_max, sizeof(int) * mem_size);
    *l_max -= k_min_post;


    for(cnt1 = k_min_post; cnt1 <= k_max_post; cnt1++){
      if(l_min_post[cnt1] < INF){
        /* new memsize */
        mem_size        = (l_max_post[cnt1] - l_min_post[cnt1] + 1)/2 + 1;
        /* reshift the pointer */
        (*array)[cnt1]  += (*l_min)[cnt1]/2;

        int shift       = (l_min_post[cnt1]%2 == (*l_min)[cnt1]%2) ? 0 : 1;
        /* eliminate unused memory in front of actual data */
        unsigned int    start = (l_min_post[cnt1] - (*l_min)[cnt1])/2 + shift;
        if(start > 0)
          memmove((FLT_OR_DBL *)((*array)[cnt1]), (FLT_OR_DBL *)((*array)[cnt1])+start, sizeof(FLT_OR_DBL) * mem_size);
        (*array)[cnt1]  = (FLT_OR_DBL *) realloc((*array)[cnt1], sizeof(FLT_OR_DBL) * mem_size);

        (*array)[cnt1]  -= l_min_post[cnt1]/2;
      }
      else{
        /* free according memory */
        (*array)[cnt1] += (*l_min)[cnt1]/2;
        free((*array)[cnt1]);
      }

      (*l_min)[cnt1] = l_min_post[cnt1];
      (*l_max)[cnt1] = l_max_post[cnt1];
    }
  }
  else{
    /* we have to free all unused memory */
    for(cnt1 = *k_min; cnt1 <= *k_max; cnt1++){
      (*array)[cnt1] += (*l_min)[cnt1]/2;
      free((*array)[cnt1]);
    }
    (*l_min) += *k_min;
    (*l_max) += *k_min;
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
  *k_min      = k_min_post;
  *k_max      = k_max_post;
}


INLINE PRIVATE void preparePosteriorBoundaries(int size, int shift, int *min_k, int *max_k, int **min_l, int **max_l){
  int i;
  *min_k  = INF;
  *max_k  = 0;

  *min_l  = (int *)space(sizeof(int) * size);
  *max_l  = (int *)space(sizeof(int) * size);

  for(i = 0; i < size; i++){
    (*min_l)[i] = INF;
    (*max_l)[i] = 0;
  }

  *min_l  -= shift;
  *max_l  -= shift;
}

INLINE PRIVATE void updatePosteriorBoundaries(int d1, int d2, int *min_k, int *max_k, int **min_l, int **max_l){
  (*min_l)[d1]  = MIN2((*min_l)[d1], d2);
  (*max_l)[d1]  = MAX2((*max_l)[d1], d2);
  *min_k        = MIN2(*min_k, d1);
  *max_k        = MAX2(*max_k, d1);
}

INLINE PRIVATE  void  prepareBoundaries(int min_k_pre, int max_k_pre, int min_l_pre, int max_l_pre, int bpdist, int *min_k, int *max_k, int **min_l, int **max_l){
  int cnt;
  int mem = max_k_pre - min_k_pre + 1;

  *min_k  = min_k_pre;
  *max_k  = max_k_pre;
  *min_l  = (int *) space(sizeof(int) * mem);
  *max_l  = (int *) space(sizeof(int) * mem);

  *min_l  -= min_k_pre;
  *max_l  -= min_k_pre;

  /* for each k guess the according minimum l*/
  for(cnt = min_k_pre; cnt <= max_k_pre; cnt++){
    (*min_l)[cnt] = min_l_pre;
    (*max_l)[cnt] = max_l_pre;
    while((*min_l)[cnt] + cnt < bpdist) (*min_l)[cnt]++;
    if((bpdist % 2) != (((*min_l)[cnt] + cnt) % 2)) (*min_l)[cnt]++;
  }
}

INLINE PRIVATE  void  prepareArray(FLT_OR_DBL ***array, int min_k, int max_k, int *min_l, int *max_l){
  int i, j, mem;
  *array  = (FLT_OR_DBL **)space(sizeof(FLT_OR_DBL *) * (max_k - min_k + 1));
  *array  -= min_k;

  for(i = min_k; i <= max_k; i++){
    mem = (max_l[i] - min_l[i] + 1)/2 + 1;
    (*array)[i] = (FLT_OR_DBL *)space(sizeof(FLT_OR_DBL) * mem);
    /*
    for(j = 0; j < mem; j++)
      (*array)[i][j] = 0.;
    */
    (*array)[i]  -= min_l[i]/2;
  }
}

