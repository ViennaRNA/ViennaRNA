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
#include "utils.h"
#include "energy_par.h"
#include "fold_vars.h"
#include "fold.h"
#include "pair_mat.h"
#include "loop_energies.h"
#include "mm.h"
#include "params.h"
#ifdef _OPENMP
#include <omp.h>
#endif
#include "2Dfold.h"

/*
#################################
# GLOBAL VARIABLES              #
#################################
*/
int compute_2Dfold_F3 = 0;

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
PRIVATE void  mfe_linear(TwoDfold_vars *vars);
PRIVATE void  mfe_circ(TwoDfold_vars *vars);

PRIVATE void  initialize_TwoDfold_vars(TwoDfold_vars *vars);
PUBLIC  void  update_TwoDfold_params(TwoDfold_vars *vars);
PRIVATE void  make_ptypes(TwoDfold_vars *vars);

PRIVATE void  backtrack_f5(unsigned int j, int k, int l, char *structure, TwoDfold_vars *vars);
PRIVATE void  backtrack_c(unsigned int i, unsigned int j, int k, int l, char *structure, TwoDfold_vars *vars);
PRIVATE void  backtrack_m(unsigned int i, unsigned int j, int k, int l, char *structure, TwoDfold_vars *vars);
PRIVATE void  backtrack_m1(unsigned int i, unsigned int j, int k, int l, char *structure, TwoDfold_vars *vars);
PRIVATE void  backtrack_fc(int k, int l, char *structure, TwoDfold_vars *vars);
PRIVATE void  backtrack_m2(unsigned int i, int k, int l, char *structure, TwoDfold_vars *vars);

PRIVATE void  adjustArrayBoundaries(int ***array, int *k_min, int *k_max, int **l_min, int **l_max, int k_min_real, int k_max_real, int *l_min_real, int *l_max_real);
INLINE  PRIVATE void  preparePosteriorBoundaries(int size, int shift, int *min_k, int *max_k, int **min_l, int **max_l);
INLINE  PRIVATE void  updatePosteriorBoundaries(int d1, int d2, int *min_k, int *max_k, int **min_l, int **max_l);
INLINE  PRIVATE void  prepareBoundaries(int min_k_pre, int max_k_pre, int min_l_pre, int max_l_pre, int bpdist, int *min_k, int *max_k, int **min_l, int **max_l);
INLINE  PRIVATE void  prepareArray(int ***array, int min_k, int max_k, int *min_l, int *max_l);
INLINE  PRIVATE void  prepareArray2(unsigned long ***array, int min_k, int max_k, int *min_l, int *max_l);



/*
#################################
# BEGIN OF FUNCTION DEFINITIONS #
#################################
*/

PUBLIC TwoDfold_vars *get_TwoDfold_variables(const char *seq, const char *structure1, const char *structure2, int circ){
  unsigned int size, length, i;
  int *index;
  TwoDfold_vars *vars;
  length = strlen(seq);
  vars = (TwoDfold_vars *)malloc(sizeof(TwoDfold_vars));
  vars->sequence     = (char *)space(length + 1);
  strcpy(vars->sequence, seq);
  vars->seq_length   = length;
  if(vars->seq_length < 1) nrerror("get_TwoDfold_variables: sequence must be longer than 0");
  size                        = ((length + 1) * (length + 2)/2);

  vars->reference_pt1   = make_pair_table(structure1);
  vars->reference_pt2   = make_pair_table(structure2);
  vars->referenceBPs1   = make_referenceBP_array(vars->reference_pt1, TURN);
  vars->referenceBPs2   = make_referenceBP_array(vars->reference_pt2, TURN);
  vars->bpdist          = compute_BPdifferences(vars->reference_pt1, vars->reference_pt2, TURN);
  vars->do_backtrack    = 1;
  vars->dangles         = dangles;
  vars->circ            = circ;
  vars->temperature     = temperature;
  vars->ptype           = space(sizeof(char) * size);
  vars->P               = NULL;
  vars->S               = NULL;
  vars->S1              = NULL;
  vars->my_iindx        = get_iindx(length);
  index                 = vars->my_iindx;
  /* compute maximum matching with reference structure 1 disallowed */
  vars->mm1          = maximumMatchingConstraint(vars->sequence, vars->reference_pt1);
  /* compute maximum matching with reference structure 2 disallowed */
  vars->mm2          = maximumMatchingConstraint(vars->sequence, vars->reference_pt2);

  vars->maxD1        = vars->mm1[index[1]-length] + vars->referenceBPs1[index[1]-length];
  vars->maxD2        = vars->mm2[index[1]-length] + vars->referenceBPs2[index[1]-length];

  /* allocate memory for the energy matrices and min-/max-index helper arrays */
  vars->E_C              = (int ***) space(sizeof(int **)  * size);
  vars->l_min_values     = (int **)  space(sizeof(int *)   * size);
  vars->l_max_values     = (int **)  space(sizeof(int *)   * size);
  vars->k_min_values     = (int *)   space(sizeof(int)     * size);
  vars->k_max_values     = (int *)   space(sizeof(int)     * size);

  vars->E_F5             = (int ***) space(sizeof(int **)  * (length + 1));
  vars->l_min_values_f   = (int **)  space(sizeof(int *)   * (length + 1));
  vars->l_max_values_f   = (int **)  space(sizeof(int *)   * (length + 1));
  vars->k_min_values_f   = (int *)   space(sizeof(int)     * (length + 1));
  vars->k_max_values_f   = (int *)   space(sizeof(int)     * (length + 1));

  if(compute_2Dfold_F3){
    vars->E_F3             = (int ***) space(sizeof(int **)  * (length + 1));
    vars->l_min_values_f3  = (int **)  space(sizeof(int *)   * (length + 1));
    vars->l_max_values_f3  = (int **)  space(sizeof(int *)   * (length + 1));
    vars->k_min_values_f3  = (int *)   space(sizeof(int)     * (length + 1));
    vars->k_max_values_f3  = (int *)   space(sizeof(int)     * (length + 1));
  }
  else vars->E_F3 = NULL;

  vars->E_M              = (int ***) space(sizeof(int **)  * size);
  vars->l_min_values_m   = (int **)  space(sizeof(int *)   * size);
  vars->l_max_values_m   = (int **)  space(sizeof(int *)   * size);
  vars->k_min_values_m   = (int *)   space(sizeof(int)     * size);
  vars->k_max_values_m   = (int *)   space(sizeof(int)     * size);

  vars->E_M1             = (int ***) space(sizeof(int **)  * size);
  vars->l_min_values_m1  = (int **)  space(sizeof(int *)   * size);
  vars->l_max_values_m1  = (int **)  space(sizeof(int *)   * size);
  vars->k_min_values_m1  = (int *)   space(sizeof(int)     * size);
  vars->k_max_values_m1  = (int *)   space(sizeof(int)     * size);

#ifdef COUNT_STATES
  vars->N_C              = (unsigned long ***) space(sizeof(unsigned long **)  * size);
  vars->N_F5             = (unsigned long ***) space(sizeof(unsigned long **)  * (length + 1));
  vars->N_M              = (unsigned long ***) space(sizeof(unsigned long **)  * size);
  vars->N_M1             = (unsigned long ***) space(sizeof(unsigned long **)  * size);
#endif


  if(circ){
    vars->E_M2_rem       = (int *)   space(sizeof(int)     * (length + 1));
    vars->E_M2            = (int ***) space(sizeof(int **)  * (length + 1));
    vars->l_min_values_m2 = (int **)  space(sizeof(int *)   * (length + 1));
    vars->l_max_values_m2 = (int **)  space(sizeof(int *)   * (length + 1));
    vars->k_min_values_m2 = (int *)   space(sizeof(int)     * (length + 1));
    vars->k_max_values_m2 = (int *)   space(sizeof(int)     * (length + 1));
  }
  else{
    vars->E_M2_rem       = NULL;
    vars->E_M2            = NULL;
    vars->l_min_values_m2 = NULL;
    vars->l_max_values_m2 = NULL;
    vars->k_min_values_m2 = NULL;
    vars->k_max_values_m2 = NULL;
  }

  vars->E_Fc              = NULL;
  vars->E_FcH             = NULL;
  vars->E_FcI             = NULL;
  vars->E_FcM             = NULL;

  vars->E_Fc_rem         = INF;
  vars->E_FcH_rem        = INF;
  vars->E_FcI_rem        = INF;
  vars->E_FcM_rem        = INF;

  vars->E_C_rem          = (int *) space(sizeof(int) * size);
  vars->E_M_rem          = (int *) space(sizeof(int) * size);
  vars->E_M1_rem         = (int *) space(sizeof(int) * size);
  vars->E_F5_rem         = (int *) space(sizeof(int) * (length+1));
  /* init rest arrays */
  for(i=0;i<size;i++){
    vars->E_C_rem[i] = vars->E_M_rem[i] = vars->E_M1_rem[i] = INF;
  }
  for(i=0;i<=length;i++)
    vars->E_F5_rem[i] = INF;
  if(vars->E_M2_rem)
    for(i=0;i<=length;i++)
      vars->E_M2_rem[i] = INF;

  return vars;
}

PUBLIC void destroy_TwoDfold_variables(TwoDfold_vars *vars){
  unsigned int i, j, ij;
  int cnt1, cnt2, cnt3, cnt4;
  if(vars == NULL) return;

  free(vars->E_C_rem);
  free(vars->E_M_rem);
  free(vars->E_M1_rem);
  free(vars->E_F5_rem);
  if(vars->E_M2_rem) free(vars->E_M2_rem);

#ifdef _OPENMP
  #pragma omp sections private(i,j,ij,cnt1,cnt2,cnt3,cnt4)
  {

  #pragma omp section
  {
#endif

#ifdef COUNT_STATES
  if(vars->N_C != NULL){
    for(i = 1; i < vars->seq_length; i++){
      for(j = i; j <= vars->seq_length; j++){
        ij = vars->my_iindx[i] - j;
        if(!vars->N_C[ij]) continue;
        for(cnt1 = vars->k_min_values[ij]; cnt1 <= vars->k_max_values[ij]; cnt1++)
          if(vars->l_min_values[ij][cnt1] < INF){
            vars->N_C[ij][cnt1] += vars->l_min_values[ij][cnt1]/2;
            free(vars->N_C[ij][cnt1]);
          }
        if(vars->k_min_values[ij] < INF){
          vars->N_C[ij] += vars->k_min_values[ij];
          free(vars->N_C[ij]);
        }
      }
    }
    free(vars->N_C);
  }
#endif

  if(vars->E_C != NULL){
    for(i = 1; i < vars->seq_length; i++){
      for(j = i; j <= vars->seq_length; j++){
        ij = vars->my_iindx[i] - j;
        if(!vars->E_C[ij]) continue;
        for(cnt1 = vars->k_min_values[ij]; cnt1 <= vars->k_max_values[ij]; cnt1++)
          if(vars->l_min_values[ij][cnt1] < INF){
            vars->E_C[ij][cnt1] += vars->l_min_values[ij][cnt1]/2;
            free(vars->E_C[ij][cnt1]);
          }
        if(vars->k_min_values[ij] < INF){
          vars->E_C[ij] += vars->k_min_values[ij];
          free(vars->E_C[ij]);
          vars->l_min_values[ij] += vars->k_min_values[ij];
          vars->l_max_values[ij] += vars->k_min_values[ij];
          free(vars->l_min_values[ij]);
          free(vars->l_max_values[ij]);
        }
      }
    }
    free(vars->E_C);
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

#ifdef COUNT_STATES
  if(vars->N_M != NULL){
    for(i = 1; i < vars->seq_length; i++){
      for(j = i; j <= vars->seq_length; j++){
        ij = vars->my_iindx[i] - j;
        if(!vars->N_M[ij]) continue;
        for(cnt1 = vars->k_min_values_m[ij]; cnt1 <= vars->k_max_values_m[ij]; cnt1++)
          if(vars->l_min_values_m[ij][cnt1] < INF){
            vars->N_M[ij][cnt1] += vars->l_min_values_m[ij][cnt1]/2;
            free(vars->N_M[ij][cnt1]);
          }
        if(vars->k_min_values_m[ij] < INF){
          vars->N_M[ij] += vars->k_min_values_m[ij];
          free(vars->N_M[ij]);
        }
      }
    }
    free(vars->N_M);
  }
#endif

  if(vars->E_M != NULL){
    for(i = 1; i < vars->seq_length; i++){
      for(j = i; j <= vars->seq_length; j++){
        ij = vars->my_iindx[i] - j;
        if(!vars->E_M[ij]) continue;
        for(cnt1 = vars->k_min_values_m[ij]; cnt1 <= vars->k_max_values_m[ij]; cnt1++)
          if(vars->l_min_values_m[ij][cnt1] < INF){
            vars->E_M[ij][cnt1] += vars->l_min_values_m[ij][cnt1]/2;
            free(vars->E_M[ij][cnt1]);
          }
        if(vars->k_min_values_m[ij] < INF){
          vars->E_M[ij] += vars->k_min_values_m[ij];
          free(vars->E_M[ij]);
          vars->l_min_values_m[ij] += vars->k_min_values_m[ij];
          vars->l_max_values_m[ij] += vars->k_min_values_m[ij];
          free(vars->l_min_values_m[ij]);
          free(vars->l_max_values_m[ij]);
        }
      }
    }
    free(vars->E_M);
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

#ifdef COUNT_STATES
  if(vars->N_M1 != NULL){
    for(i = 1; i < vars->seq_length; i++){
      for(j = i; j <= vars->seq_length; j++){
        ij = vars->my_iindx[i] - j;
        if(!vars->N_M1[ij]) continue;
        for(cnt1 = vars->k_min_values_m1[ij]; cnt1 <= vars->k_max_values_m1[ij]; cnt1++)
          if(vars->l_min_values_m1[ij][cnt1] < INF){
            vars->N_M1[ij][cnt1] += vars->l_min_values_m1[ij][cnt1]/2;
            free(vars->N_M1[ij][cnt1]);
          }
        if(vars->k_min_values_m1[ij] < INF){
          vars->N_M1[ij] += vars->k_min_values_m1[ij];
          free(vars->N_M1[ij]);
        }
      }
    }
    free(vars->N_M1);
  }
#endif

  if(vars->E_M1 != NULL){
    for(i = 1; i < vars->seq_length; i++){
      for(j = i; j <= vars->seq_length; j++){
        ij = vars->my_iindx[i] - j;
        if(!vars->E_M1[ij]) continue;
        for(cnt1 = vars->k_min_values_m1[ij]; cnt1 <= vars->k_max_values_m1[ij]; cnt1++)
          if(vars->l_min_values_m1[ij][cnt1] < INF){
            vars->E_M1[ij][cnt1] += vars->l_min_values_m1[ij][cnt1]/2;
            free(vars->E_M1[ij][cnt1]);
          }
        if(vars->k_min_values_m1[ij] < INF){
          vars->E_M1[ij] += vars->k_min_values_m1[ij];
          free(vars->E_M1[ij]);
          vars->l_min_values_m1[ij] += vars->k_min_values_m1[ij];
          vars->l_max_values_m1[ij] += vars->k_min_values_m1[ij];
          free(vars->l_min_values_m1[ij]);
          free(vars->l_max_values_m1[ij]);
        }
      }
    }
    free(vars->E_M1);
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
  if(vars->E_M2 != NULL){
    for(i = 1; i < vars->seq_length-TURN-1; i++){
      if(!vars->E_M2[i]) continue;
      for(cnt1 = vars->k_min_values_m2[i]; cnt1 <= vars->k_max_values_m2[i]; cnt1++)
        if(vars->l_min_values_m2[i][cnt1] < INF){
          vars->E_M2[i][cnt1] += vars->l_min_values_m2[i][cnt1]/2;
          free(vars->E_M2[i][cnt1]);
        }
      if(vars->k_min_values_m2[i] < INF){
        vars->E_M2[i] += vars->k_min_values_m2[i];
        free(vars->E_M2[i]);
        vars->l_min_values_m2[i] += vars->k_min_values_m2[i];
        vars->l_max_values_m2[i] += vars->k_min_values_m2[i];
        free(vars->l_min_values_m2[i]);
        free(vars->l_max_values_m2[i]);
      }
    }
    free(vars->E_M2);
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

#ifdef COUNT_STATES
  if(vars->N_F5 != NULL){
    for(i = 1; i <= vars->seq_length; i++){
      if(!vars->N_F5[i]) continue;
      for(cnt1 = vars->k_min_values_f[i]; cnt1 <= vars->k_max_values_f[i]; cnt1++)
        if(vars->l_min_values_f[i][cnt1] < INF){
          vars->N_F5[i][cnt1] += vars->l_min_values_f[i][cnt1]/2;
          free(vars->N_F5[i][cnt1]);
        }
      if(vars->k_min_values_f[i] < INF){
        vars->N_F5[i] += vars->k_min_values_f[i];
        free(vars->N_F5[i]);
      }
    }
    free(vars->N_F5);
  }
#endif

  if(vars->E_F5 != NULL){
    for(i = 1; i <= vars->seq_length; i++){
      if(!vars->E_F5[i]) continue;
      for(cnt1 = vars->k_min_values_f[i]; cnt1 <= vars->k_max_values_f[i]; cnt1++)
        if(vars->l_min_values_f[i][cnt1] < INF){
          vars->E_F5[i][cnt1] += vars->l_min_values_f[i][cnt1]/2;
          free(vars->E_F5[i][cnt1]);
        }
      if(vars->k_min_values_f[i] < INF){
        vars->E_F5[i] += vars->k_min_values_f[i];
        free(vars->E_F5[i]);
        vars->l_min_values_f[i] += vars->k_min_values_f[i];
        vars->l_max_values_f[i] += vars->k_min_values_f[i];
        free(vars->l_min_values_f[i]);
        free(vars->l_max_values_f[i]);
      }
    }
    free(vars->E_F5);
    free(vars->l_min_values_f);
    free(vars->l_max_values_f);
    free(vars->k_min_values_f);
    free(vars->k_max_values_f);
  }

  if(vars->E_F3 != NULL){
    for(i = 1; i <= vars->seq_length; i++){
      if(!vars->E_F3[i]) continue;
      for(cnt1 = vars->k_min_values_f3[i]; cnt1 <= vars->k_max_values_f3[i]; cnt1++)
        if(vars->l_min_values_f3[i][cnt1] < INF){
          vars->E_F3[i][cnt1] += vars->l_min_values_f3[i][cnt1]/2;
          free(vars->E_F3[i][cnt1]);
        }
      if(vars->k_min_values_f3[i] < INF){
        vars->E_F3[i] += vars->k_min_values_f3[i];
        free(vars->E_F3[i]);
        vars->l_min_values_f3[i] += vars->k_min_values_f3[i];
        vars->l_max_values_f3[i] += vars->k_min_values_f3[i];
        free(vars->l_min_values_f3[i]);
        free(vars->l_max_values_f3[i]);
      }
    }
    free(vars->E_F3);
    free(vars->l_min_values_f3);
    free(vars->l_max_values_f3);
    free(vars->k_min_values_f3);
    free(vars->k_max_values_f3);
  }

#ifdef _OPENMP
  }
  #pragma omp section
  {
#endif
  if(vars->E_Fc != NULL){
    for(cnt1 = vars->k_min_values_fc; cnt1 <= vars->k_max_values_fc; cnt1++)
      if(vars->l_min_values_fc[cnt1] < INF){
        vars->E_Fc[cnt1] += vars->l_min_values_fc[cnt1]/2;
        free(vars->E_Fc[cnt1]);
      }
    if(vars->k_min_values_fc < INF){
      vars->E_Fc += vars->k_min_values_fc;
      free(vars->E_Fc);
      vars->l_min_values_fc += vars->k_min_values_fc;
      vars->l_max_values_fc += vars->k_min_values_fc;
      free(vars->l_min_values_fc);
      free(vars->l_max_values_fc);
    }
  }
#ifdef _OPENMP
  }
  #pragma omp section
  {
#endif
  if(vars->E_FcI != NULL){
    for(cnt1 = vars->k_min_values_fcI; cnt1 <= vars->k_max_values_fcI; cnt1++)
      if(vars->l_min_values_fcI[cnt1] < INF){
        vars->E_FcI[cnt1] += vars->l_min_values_fcI[cnt1]/2;
        free(vars->E_FcI[cnt1]);
      }
    if(vars->k_min_values_fcI < INF){
      vars->E_FcI += vars->k_min_values_fcI;
      free(vars->E_FcI);
      vars->l_min_values_fcI += vars->k_min_values_fcI;
      vars->l_max_values_fcI += vars->k_min_values_fcI;
      free(vars->l_min_values_fcI);
      free(vars->l_max_values_fcI);
    }
  }
#ifdef _OPENMP
  }
  #pragma omp section
  {
#endif
  if(vars->E_FcH != NULL){
    for(cnt1 = vars->k_min_values_fcH; cnt1 <= vars->k_max_values_fcH; cnt1++)
      if(vars->l_min_values_fcH[cnt1] < INF){
        vars->E_FcH[cnt1] += vars->l_min_values_fcH[cnt1]/2;
        free(vars->E_FcH[cnt1]);
      }
    if(vars->k_min_values_fcH < INF){
      vars->E_FcH += vars->k_min_values_fcH;
      free(vars->E_FcH);
      vars->l_min_values_fcH += vars->k_min_values_fcH;
      vars->l_max_values_fcH += vars->k_min_values_fcH;
      free(vars->l_min_values_fcH);
      free(vars->l_max_values_fcH);
    }
  }
#ifdef _OPENMP
  }
  #pragma omp section
  {
#endif
  if(vars->E_FcM != NULL){
    for(cnt1 = vars->k_min_values_fcM; cnt1 <= vars->k_max_values_fcM; cnt1++)
      if(vars->l_min_values_fcM[cnt1] < INF){
        vars->E_FcM[cnt1] += vars->l_min_values_fcM[cnt1]/2;
        free(vars->E_FcM[cnt1]);
      }
    if(vars->k_min_values_fcM < INF){
      vars->E_FcM += vars->k_min_values_fcM;
      free(vars->E_FcM);
      vars->l_min_values_fcM += vars->k_min_values_fcM;
      vars->l_max_values_fcM += vars->k_min_values_fcM;
      free(vars->l_min_values_fcM);
      free(vars->l_max_values_fcM);
    }
  }

#ifdef _OPENMP
  }
  #pragma omp section
  {
#endif

  if(vars->P != NULL)             free(vars->P);
  if(vars->sequence != NULL)      free(vars->sequence);
  if(vars->reference_pt1 != NULL) free(vars->reference_pt1);
  if(vars->reference_pt2 != NULL) free(vars->reference_pt2);
  if(vars->referenceBPs1 != NULL) free(vars->referenceBPs1);
  if(vars->referenceBPs2 != NULL) free(vars->referenceBPs2);
  if(vars->ptype != NULL)         free(vars->ptype);
  if(vars->S != NULL)             free(vars->S);
  if(vars->S1 != NULL)            free(vars->S1);

  if(vars->mm1 != NULL)           free(vars->mm1);
  if(vars->mm2 != NULL)           free(vars->mm2);
  if(vars->bpdist != NULL)        free(vars->bpdist);
#ifdef _OPENMP
  }
  }
#endif

  if(vars->my_iindx != NULL)      free(vars->my_iindx);

  free(vars);
}

PRIVATE void initialize_TwoDfold_vars(TwoDfold_vars *vars){
  update_TwoDfold_params(vars);
  /* this call updates the params in the ViennaRNA fold.o which is a global, so be careful
  *  whith calling it parallel... need a workarround or fix of ViennaRNA fold stuff
  */
  update_fold_params();
}


PUBLIC TwoDfold_solution **TwoDfold(TwoDfold_vars *vars, int distance1, int distance2){
  unsigned int i, d1, d2;
  unsigned int maxD1;
  unsigned int maxD2;
  unsigned int mm;
  unsigned int length;
  TwoDfold_solution **output;

  initialize_TwoDfold_vars(vars);
  if(fabs(vars->P->temperature - temperature)>1e-6) update_TwoDfold_params(vars);
  vars->S   = encode_sequence(vars->sequence, 0);
  vars->S1  = encode_sequence(vars->sequence, 1);

  make_ptypes(vars);

  maxD1 = vars->maxD1;
  maxD2 = vars->maxD2;

  if(distance1 >= 0){
    if((unsigned int)distance1 > maxD1)
      fprintf(stderr,
              "limiting maximum basepair distance 1 to %u\n",
              maxD1);
    else
      maxD1 = (unsigned int)distance1;
  }

  if(distance2 >= 0){
    if((unsigned int)distance2 > maxD2)
      fprintf(stderr,
              "limiting maximum basepair distance 2 to %u\n",
              maxD2);
    else
      maxD2 = (unsigned int)distance2;
  }

  vars->maxD1 = maxD1;
  vars->maxD2 = maxD2;
  output = (TwoDfold_solution **)space((vars->maxD1+1) * sizeof(TwoDfold_solution *));

  mfe_linear(vars);
  if(vars->circ) mfe_circ(vars);

  length = vars->seq_length;

  for(d1=0; d1<=maxD1;d1++){
    output[d1] = (TwoDfold_solution *)space((vars->maxD2+1)*sizeof(TwoDfold_solution));
#ifdef _OPENMP
  #pragma omp parallel for private(d2)
#endif
    for(d2=0; d2<=maxD2;d2++){
      output[d1][d2].en = (float)INF/(float)100.;
      output[d1][d2].s  = NULL;
    }
    if(     (d1 >= ((vars->circ) ? vars->k_min_values_fc : vars->k_min_values_f[length]))
        &&  (d1 <= ((vars->circ) ? vars->k_max_values_fc : vars->k_max_values_f[length]))){
#ifdef _OPENMP
  #pragma omp parallel for private(d2, i)
#endif
      for(  d2  = ((vars->circ) ? vars->l_min_values_fc[d1] : vars->l_min_values_f[length][d1]);
            d2 <= ((vars->circ) ? vars->l_max_values_fc[d1] : vars->l_max_values_f[length][d1]);
            d2 += 2){
        output[d1][d2].en = (float)((vars->circ) ? vars->E_Fc[d1][d2/2] : vars->E_F5[length][d1][d2/2])/(float)100.;
        if(vars->do_backtrack && (output[d1][d2].en != (float)INF/(float)100.)){
          char *mfe_structure = (char *)space(length+1);
          for(i=0;i<length;i++) mfe_structure[i] = '.';
          mfe_structure[i] = '\0';
          (vars->circ) ? backtrack_fc(d1, d2, mfe_structure, vars) : backtrack_f5(length, d1, d2, mfe_structure, vars);
          output[d1][d2].s = mfe_structure;
        }
      }
    }

  }
  return output;
}

PUBLIC TwoDfold_solution *TwoDfoldList(TwoDfold_vars *vars, int distance1, int distance2){
  unsigned int  i, d1, d2;
  unsigned int  maxD1;
  unsigned int  maxD2;
  unsigned int  mm;
  unsigned int  length;
  unsigned int  counter = 0;
  int           en = 0;
  TwoDfold_solution *output;

  initialize_TwoDfold_vars(vars);
  if(fabs(vars->P->temperature - temperature)>1e-6) update_TwoDfold_params(vars);
  vars->S   = encode_sequence(vars->sequence, 0);
  vars->S1  = encode_sequence(vars->sequence, 1);

  make_ptypes(vars);

  maxD1 = vars->maxD1;
  maxD2 = vars->maxD2;

  if(distance1 >= 0){
    if((unsigned int)distance1 > maxD1)
      fprintf(stderr,
              "TwoDfoldList@2Dfold.c: limiting maximum basepair distance 1 to %u\n",
              maxD1);
    else
      maxD1 = (unsigned int)distance1;
  }

  if(distance2 >= 0){
    if((unsigned int)distance2 > maxD2)
      fprintf(stderr,
              "TwoDfoldList@2Dfold.c: limiting maximum basepair distance 2 to %u\n",
              maxD2);
    else
      maxD2 = (unsigned int)distance2;
  }

  vars->maxD1 = maxD1;
  vars->maxD2 = maxD2;
  output = (TwoDfold_solution *)space((((vars->maxD1+1)*(vars->maxD2+2))/2 + 2) * sizeof(TwoDfold_solution));

  mfe_linear(vars);
  if(vars->circ) mfe_circ(vars);

  length = vars->seq_length;

  for(d1=0; d1<=maxD1;d1++){
    if((d1 >= ((vars->circ) ? vars->k_min_values_fc : vars->k_min_values_f[length]))
        &&  (d1 <= ((vars->circ) ? vars->k_max_values_fc : vars->k_max_values_f[length]))){
      for(d2  = ((vars->circ) ? vars->l_min_values_fc[d1] : vars->l_min_values_f[length][d1]);
          d2 <= ((vars->circ) ? vars->l_max_values_fc[d1] : vars->l_max_values_f[length][d1]);
          d2 += 2){
        en = ((vars->circ) ? vars->E_Fc[d1][d2/2] : vars->E_F5[length][d1][d2/2]);
        if(en == INF) continue;
        output[counter].k   = d1;
        output[counter].l   = d2;
        output[counter].en  = (float)en/(float)100.;
        if(vars->do_backtrack){
          char *mfe_structure = (char *)space(length+1);
          for(i=0;i<length;i++) mfe_structure[i] = '.';
          mfe_structure[i] = '\0';
          (vars->circ) ? backtrack_fc((int)d1, (int)d2, mfe_structure, vars) : backtrack_f5(length, (int)d1, (int)d2, mfe_structure, vars);
          output[counter].s = mfe_structure;
        }
        else output[counter].s = NULL;
        counter++;
      }
    }
  }

  /* store entry for remaining partition if it exists */
  en = ((vars->circ) ? vars->E_Fc_rem : vars->E_F5_rem[length]);
  if(en != INF){
    output[counter].k   = -1;
    output[counter].l   = -1;
    output[counter].en  =  (float)en/(float)100.;
    if(vars->do_backtrack){
      char *mfe_structure = (char *)space(length+1);
      for(i=0;i<length;i++) mfe_structure[i] = '.';
      mfe_structure[i] = '\0';
      (vars->circ) ? backtrack_fc(-1, -1, mfe_structure, vars) : backtrack_f5(length, -1, -1, mfe_structure, vars);
      output[counter].s = mfe_structure;
    }
    else output[counter].s = NULL;
    counter++;
  }

  /* insert end-marker entry */
  output[counter].k = output[counter].l = INF;
  counter++;

  /* resize to actual dataset amount */
  output = (TwoDfold_solution*)xrealloc(output, sizeof(TwoDfold_solution) * counter);
  return output;
}


PUBLIC char *TwoDfold_backtrack_f5(unsigned int j, int k, int l, TwoDfold_vars *vars){
  unsigned int i;
  char *mfe_structure = (char *)space(j+1);
  if(j < TURN + 2) return NULL;

  for(i=0; i < j; i++) mfe_structure[i] = '.';
  mfe_structure[i] = '\0';

  backtrack_f5(j, k, l, mfe_structure, vars);
  return mfe_structure;
}

PRIVATE void mfe_linear(TwoDfold_vars *vars){

  unsigned int  d, i, j, ij, k, maxD1, maxD2, seq_length, dia, dib, dja, djb, *referenceBPs1, *referenceBPs2, *mm1, *mm2, *bpdist;
  int           ***E_C, ***E_M, ***E_M1, ***E_F5, cnt1, cnt2, cnt3, cnt4, d1, d2, energy, dangles, temp2, type, additional_en, *my_iindx, circ;
  int           **l_min_values, **l_max_values, **l_min_values_m, **l_max_values_m,**l_min_values_m1, **l_max_values_m1,**l_min_values_f, **l_max_values_f;
  int           *k_min_values, *k_max_values, *k_min_values_m, *k_max_values_m,*k_min_values_m1, *k_max_values_m1,*k_min_values_f, *k_max_values_f;
  short         *S1, *reference_pt1, *reference_pt2;
  char          *sequence, *ptype;
  paramT        *P;
  int           ***E_F3;
  int           **l_min_values_f3, **l_max_values_f3, *k_min_values_f3, *k_max_values_f3;

  /* dereferenciate things we often need */
  P               = vars->P;
  sequence        = vars->sequence;
  seq_length      = vars->seq_length;
  maxD1           = vars->maxD1;
  maxD2           = vars->maxD2;
  S1              = vars->S1;
  ptype           = vars->ptype;
  reference_pt1   = vars->reference_pt1;
  reference_pt2   = vars->reference_pt2;
  my_iindx        = vars->my_iindx;
  referenceBPs1   = vars->referenceBPs1;
  referenceBPs2   = vars->referenceBPs2;
  dangles         = vars->dangles;
  mm1             = vars->mm1;
  mm2             = vars->mm2;
  bpdist          = vars->bpdist;
  circ            = vars->circ;

  E_F5            = vars->E_F5;
  l_min_values_f  = vars->l_min_values_f;
  l_max_values_f  = vars->l_max_values_f;
  k_min_values_f  = vars->k_min_values_f;
  k_max_values_f  = vars->k_max_values_f;

  E_C             = vars->E_C;
  l_min_values    = vars->l_min_values;
  l_max_values    = vars->l_max_values;
  k_min_values    = vars->k_min_values;
  k_max_values    = vars->k_max_values;

  E_M             = vars->E_M;
  l_min_values_m  = vars->l_min_values_m;
  l_max_values_m  = vars->l_max_values_m;
  k_min_values_m  = vars->k_min_values_m;
  k_max_values_m  = vars->k_max_values_m;

  E_M1            = vars->E_M1;
  l_min_values_m1 = vars->l_min_values_m1;
  l_max_values_m1 = vars->l_max_values_m1;
  k_min_values_m1 = vars->k_min_values_m1;
  k_max_values_m1 = vars->k_max_values_m1;

  if(compute_2Dfold_F3){
    E_F3            = vars->E_F3;
    l_min_values_f3 = vars->l_min_values_f3;
    l_max_values_f3 = vars->l_max_values_f3;
    k_min_values_f3 = vars->k_min_values_f3;
    k_max_values_f3 = vars->k_max_values_f3;
  }



  for (d = TURN+2; d <= seq_length; d++) { /* i,j in [1..length] */
#ifdef _OPENMP
  #pragma omp parallel for private(additional_en, j, energy, temp2, i, ij, k, dia,dib,dja,djb,cnt1,cnt2,cnt3,cnt4, d1, d2)
#endif
    for (j = d; j <= seq_length; j++) {
      unsigned int p, q, pq, u, maxp, dij;
      int type_2, type, tt, no_close, base_d1, base_d2;

      i = j-d+1;
      dij = j - i - 1;
      ij = my_iindx[i]-j;
      type = ptype[ij];

      no_close = (((type==3)||(type==4))&&no_closingGU);

      //k_min_values[ij] = 0;
      //k_max_values[ij] = mm1[ij] + referenceBPs1[ij];

      if (type) {   /* we have a pair */
        /* increase or decrease distance-to-reference value depending whether (i,j) is included in
        *  reference or has to be introduced
        */
        base_d1 = ((unsigned int)reference_pt1[i] != j) ? 1 : -1;
        base_d2 = ((unsigned int)reference_pt2[i] != j) ? 1 : -1;

        /* HAIRPIN STRUCTURES */

        /* get distance to reference if closing the hairpin
        *  d = dbp(T_{i,j}, {i,j})
        */
        d1 = base_d1 + referenceBPs1[ij];
        d2 = base_d2 + referenceBPs2[ij];

        int min_k, max_k, min_l, max_l;
        int real_min_k, real_max_k, *min_l_real, *max_l_real;

        //max_k = mm1[my_iindx[i+1]-j+1] + d1;
        //max_l = mm2[my_iindx[i+1]-j+1] + d2;

        //min_k = referenceBPs1[ij] - referenceBPs1[my_iindx[i+1]-j+1] + base_d1;
        //min_l = referenceBPs2[ij] - referenceBPs2[my_iindx[i+1]-j+1] + base_d2;
        min_l = min_k = 0;
        max_k = mm1[ij] + referenceBPs1[ij];
        max_l = mm2[ij] + referenceBPs2[ij];

        prepareBoundaries(min_k,
                          max_k,
                          min_l,
                          max_l,
                          bpdist[ij],
                          &vars->k_min_values[ij],
                          &vars->k_max_values[ij],
                          &vars->l_min_values[ij],
                          &vars->l_max_values[ij]
                          );

        preparePosteriorBoundaries( vars->k_max_values[ij] - vars->k_min_values[ij] + 1,
                                    vars->k_min_values[ij],
                                    &real_min_k,
                                    &real_max_k,
                                    &min_l_real,
                                    &max_l_real
                                  );

        prepareArray( &vars->E_C[ij],
                      vars->k_min_values[ij],
                      vars->k_max_values[ij],
                      vars->l_min_values[ij],
                      vars->l_max_values[ij]
                    );

#ifdef COUNT_STATES
        prepareArray2( &vars->N_C[ij],
                      vars->k_min_values[ij],
                      vars->k_max_values[ij],
                      vars->l_min_values[ij],
                      vars->l_max_values[ij]
                    );
#endif

        /* d1 and d2 are the distancies to both references introduced by closing a hairpin structure at (i,j) */
        if((d1 >= 0) && (d2 >= 0)){
          if(((unsigned int)d1<=maxD1) && ((unsigned int)d2 <= maxD2)){
            vars->E_C[ij][d1][d2/2] = (no_close) ? FORBIDDEN : E_Hairpin(dij, type, S1[i+1], S1[j-1], sequence+i-1, P);
            updatePosteriorBoundaries(d1,
                                      d2,
                                      &real_min_k,
                                      &real_max_k,
                                      &min_l_real,
                                      &max_l_real
                                      );
#ifdef COUNT_STATES
            vars->N_C[ij][d1][d2/2] = 1;
#endif
          }
          else{
            vars->E_C_rem[ij] = (no_close) ? FORBIDDEN : E_Hairpin(dij, type, S1[i+1], S1[j-1], sequence+i-1, P);
          }
        }
        /* INTERIOR LOOP STRUCTURES */
        maxp = MIN2(j-2-TURN,i+MAXLOOP+1);
        for(p = i+1; p <= maxp; p++){
          unsigned int minq = p + TURN + 1;
          unsigned int ln_pre = dij + p;
          if(ln_pre > minq + MAXLOOP) minq = ln_pre - MAXLOOP - 1;
          for(q = minq; q < j; q++){
            pq = my_iindx[p]-q;
            /* set distance to reference structure... */
            type_2 = ptype[pq];

            if (type_2==0) continue;
            type_2 = rtype[type_2];

            /* get distance to reference if closing the interior loop
            *  d2 = dbp(S_{i,j}, S_{p.q} + {i,j})
            */
            d1 = base_d1 + referenceBPs1[ij] - referenceBPs1[pq];
            d2 = base_d2 + referenceBPs2[ij] - referenceBPs2[pq];

            if(no_closingGU)
              if(no_close||(type_2==3)||(type_2==4))
                if((p>i+1)||(q<j-1)) continue;  /* continue unless stack */

            energy = E_IntLoop(p-i-1, j-q-1, type, type_2, S1[i+1], S1[j-1], S1[p-1], S1[q+1], P);

            if(vars->E_C[pq] != NULL){
              for(cnt1 = vars->k_min_values[pq]; cnt1 <= vars->k_max_values[pq]; cnt1++){
                for(cnt2 = vars->l_min_values[pq][cnt1]; cnt2 <= vars->l_max_values[pq][cnt1]; cnt2+=2){
                  if(vars->E_C[pq][cnt1][cnt2/2] != INF){
                    if(((cnt1 + d1) <= maxD1) && ((cnt2+d2) <= maxD2)){
                        vars->E_C[ij][cnt1 + d1][(cnt2 + d2)/2] = MIN2( vars->E_C[ij][cnt1 + d1][(cnt2 + d2)/2],
                                                                        vars->E_C[pq][cnt1][cnt2/2] + energy
                                                                      );
                        updatePosteriorBoundaries(cnt1 + d1,
                                                  cnt2 + d2,
                                                  &real_min_k,
                                                  &real_max_k,
                                                  &min_l_real,
                                                  &max_l_real
                                                  );
#ifdef COUNT_STATES
                       vars->N_C[ij][cnt1 + d1][(cnt2 + d2)/2] += vars->N_C[pq][cnt1][cnt2/2];
#endif
                    }
                    /* collect all cases where d1+cnt1 or d2+cnt2 exceeds maxD1, maxD2, respectively */
                    else{
                      vars->E_C_rem[ij] = MIN2(vars->E_C_rem[ij], vars->E_C[pq][cnt1][cnt2/2] + energy);
                    }
                  }
                }
              }
            }
            /* collect all contributions where C[pq] already lies outside k_max, l_max boundary */
            if(vars->E_C_rem[pq] != INF){
              vars->E_C_rem[ij] = MIN2(vars->E_C_rem[ij], vars->E_C_rem[pq] + energy);
            }
          } /* end q-loop */
        } /* end p-loop */



        /* MULTI LOOP STRUCTURES */
        if(!no_close){

          /* dangle energies for multiloop closing stem */
          tt = rtype[type];
          temp2 = P->MLclosing;
          if(dangles == 2)
            temp2 += E_MLstem(tt, S1[j-1], S1[i+1], P);
          else
            temp2 += E_MLstem(tt, -1, -1, P);

          for(u=i+TURN+2; u<j-TURN-2;u++){
            int i1u   = my_iindx[i+1]-u;
            int u1j1  = my_iindx[u+1]-j+1;
            /* check all cases where either M or M1 are already out of scope of maxD1 and/or maxD2 */
            if(vars->E_M_rem[i1u] != INF){
              for(cnt3 = vars->k_min_values_m1[u1j1];
                  cnt3 <= vars->k_max_values_m1[u1j1];
                  cnt3++)
                for(cnt4 = vars->l_min_values_m1[u1j1][cnt3];
                    cnt4 <= vars->l_max_values_m1[u1j1][cnt3];
                    cnt4+=2){
                  if(vars->E_M1[u1j1][cnt3][cnt4/2]!= INF){
                    vars->E_C_rem[ij] = MIN2(vars->E_C_rem[ij],
                                              vars->E_M_rem[i1u]
                                            + vars->E_M1[u1j1][cnt3][cnt4/2]
                                            + temp2
                                              );
                  }
                }
              if(vars->E_M1_rem[u1j1] != INF){
                vars->E_C_rem[ij] = MIN2(vars->E_C_rem[ij],
                                          vars->E_M_rem[i1u]
                                        + vars->E_M1_rem[u1j1]
                                        + temp2
                                          );
              }
            }
            if(vars->E_M1_rem[u1j1] != INF){
              for(cnt1 = vars->k_min_values_m[i1u];
                  cnt1 <= vars->k_max_values_m[i1u];
                  cnt1++)
                for(cnt2 = vars->l_min_values_m[i1u][cnt1];
                    cnt2 <= vars->l_max_values_m[i1u][cnt1];
                    cnt2+=2)
                  if(vars->E_M[i1u][cnt1][cnt2/2] != INF){
                    vars->E_C_rem[ij] = MIN2(vars->E_C_rem[ij],
                                              vars->E_M[i1u][cnt1][cnt2/2]
                                            + vars->E_M1_rem[u1j1]
                                            + temp2
                                              );
                  }
            }
            /* get distance to reference if closing the multiloop
            *  d = dbp(S_{i,j}, {i,j} + S_{i+1,u} + S_{u+1,j-1})
            */
            if(!vars->E_M[i1u]) continue;
            if(!vars->E_M1[u1j1]) continue;

            d1 = base_d1 + referenceBPs1[ij] - referenceBPs1[i1u] - referenceBPs1[u1j1];
            d2 = base_d2 + referenceBPs2[ij] - referenceBPs2[i1u] - referenceBPs2[u1j1];

            for(cnt1 = vars->k_min_values_m[i1u];
                cnt1 <= vars->k_max_values_m[i1u];
                cnt1++)
              for(cnt2 = vars->l_min_values_m[i1u][cnt1];
                  cnt2 <= vars->l_max_values_m[i1u][cnt1];
                  cnt2+=2)
                for(cnt3 = vars->k_min_values_m1[u1j1];
                    cnt3 <= vars->k_max_values_m1[u1j1];
                    cnt3++)
                  for(cnt4 = vars->l_min_values_m1[u1j1][cnt3];
                      cnt4 <= vars->l_max_values_m1[u1j1][cnt3];
                      cnt4+=2){
                    if((vars->E_M[i1u][cnt1][cnt2/2] != INF) && (vars->E_M1[u1j1][cnt3][cnt4/2]!= INF)){
                      if(((cnt1+cnt3+d1) <= maxD1) && ((cnt2+cnt4+d2) <= maxD2)){
                        vars->E_C[ij][cnt1+cnt3+d1][(cnt2+cnt4+d2)/2] = MIN2( vars->E_C[ij][cnt1+cnt3+d1][(cnt2+cnt4+d2)/2],
                                                                              vars->E_M[i1u][cnt1][cnt2/2]
                                                                            + vars->E_M1[u1j1][cnt3][cnt4/2]
                                                                            + temp2
                                                                            );
                        updatePosteriorBoundaries(cnt1 + cnt3 + d1,
                                                  cnt2 + cnt4 + d2,
                                                  &real_min_k,
                                                  &real_max_k,
                                                  &min_l_real,
                                                  &max_l_real
                                                );
#ifdef COUNT_STATES
                        vars->N_C[ij][cnt1+cnt3+d1][(cnt2+cnt4+d2)/2] += vars->N_M[i1u][cnt1][cnt2/2] * vars->N_M1[u1j1][cnt3][cnt4/2];
#endif
                      }
                      /* collect all cases where d1+cnt1+cnt3 or d2+cnt2+cnt4 exceeds maxD1, maxD2, respectively */
                      else{
                        vars->E_C_rem[ij] = MIN2(  vars->E_C_rem[ij],
                                                    vars->E_M[i1u][cnt1][cnt2/2]
                                                  + vars->E_M1[u1j1][cnt3][cnt4/2]
                                                  + temp2
                                                  );
                      }
                    }
                  }
          }
        }

        /* resize and move memory portions of energy matrix E_C */
        adjustArrayBoundaries(&vars->E_C[ij],
                              &vars->k_min_values[ij],
                              &vars->k_max_values[ij],
                              &vars->l_min_values[ij],
                              &vars->l_max_values[ij],
                              real_min_k,
                              real_max_k,
                              min_l_real,
                              max_l_real
                              );
#ifdef COUNT_STATES
        /* actually we should adjust the array boundaries here but we might never use the count states option more than once so what....*/
#endif
      } /* end >> if (pair) << */



      /* done with c[i,j], now compute fML[i,j] */
      /* free ends ? -----------------------------------------*/


      dia = referenceBPs1[ij] - referenceBPs1[my_iindx[i+1]-j];
      dib = referenceBPs2[ij] - referenceBPs2[my_iindx[i+1]-j];
      dja = referenceBPs1[ij] - referenceBPs1[ij+1];
      djb = referenceBPs2[ij] - referenceBPs2[ij+1];

      if(dangles==2)
        temp2 = E_MLstem(type, ((i > 1) || circ) ? S1[i-1] : -1, ((j < seq_length) || circ) ? S1[j+1] : -1, P);
      else
        temp2 = E_MLstem(type, -1, -1, P);

      int min_k_guess, max_k_guess, min_l_guess, max_l_guess;
      int min_k_real_m, max_k_real_m, *min_l_real_m, *max_l_real_m;
      int min_k_real_m1, max_k_real_m1, *min_l_real_m1, *max_l_real_m1;

      min_k_guess = min_l_guess = 0;
      max_k_guess = mm1[ij] + referenceBPs1[ij];
      max_l_guess = mm2[ij] + referenceBPs2[ij];

      prepareBoundaries(min_k_guess,
                        max_k_guess,
                        min_l_guess,
                        max_l_guess,
                        bpdist[ij],
                        &vars->k_min_values_m[ij],
                        &vars->k_max_values_m[ij],
                        &vars->l_min_values_m[ij],
                        &vars->l_max_values_m[ij]
                        );

      prepareBoundaries(min_k_guess,
                        max_k_guess,
                        min_l_guess,
                        max_l_guess,
                        bpdist[ij],
                        &vars->k_min_values_m1[ij],
                        &vars->k_max_values_m1[ij],
                        &vars->l_min_values_m1[ij],
                        &vars->l_max_values_m1[ij]
                        );

      preparePosteriorBoundaries( vars->k_max_values_m[ij] - vars->k_min_values_m[ij] + 1,
                                  vars->k_min_values_m[ij],
                                  &min_k_real_m,
                                  &max_k_real_m,
                                  &min_l_real_m,
                                  &max_l_real_m
                                );
      preparePosteriorBoundaries( vars->k_max_values_m1[ij] - vars->k_min_values_m1[ij] + 1,
                                  vars->k_min_values_m1[ij],
                                  &min_k_real_m1,
                                  &max_k_real_m1,
                                  &min_l_real_m1,
                                  &max_l_real_m1
                                );

      prepareArray( &vars->E_M[ij],
                    vars->k_min_values_m[ij],
                    vars->k_max_values_m[ij],
                    vars->l_min_values_m[ij],
                    vars->l_max_values_m[ij]
                  );

      prepareArray( &vars->E_M1[ij],
                    vars->k_min_values_m1[ij],
                    vars->k_max_values_m1[ij],
                    vars->l_min_values_m1[ij],
                    vars->l_max_values_m1[ij]
                  );
#ifdef COUNT_STATES
      prepareArray2( &vars->N_M[ij],
                    vars->k_min_values_m[ij],
                    vars->k_max_values_m[ij],
                    vars->l_min_values_m[ij],
                    vars->l_max_values_m[ij]
                  );
      prepareArray2( &vars->N_M1[ij],
                    vars->k_min_values_m1[ij],
                    vars->k_max_values_m1[ij],
                    vars->l_min_values_m1[ij],
                    vars->l_max_values_m1[ij]
                  );
#endif

      /* now to the actual computations... */
      /* 1st E_M[ij] = E_M1[ij] = E_C[ij] + b */
      if(vars->E_C_rem[ij] != INF){
        vars->E_M_rem[ij] = vars->E_M1_rem[ij] = temp2 + vars->E_C_rem[ij];
      }
      if(vars->E_C[ij])
        for(cnt1 = vars->k_min_values[ij]; cnt1 <= vars->k_max_values[ij]; cnt1++){
          for(cnt2 = vars->l_min_values[ij][cnt1]; cnt2 <= vars->l_max_values[ij][cnt1]; cnt2+=2){
            if(vars->E_C[ij][cnt1][cnt2/2] != INF){
              vars->E_M[ij][cnt1][cnt2/2] = vars->E_M1[ij][cnt1][cnt2/2] = temp2 + vars->E_C[ij][cnt1][cnt2/2];
              updatePosteriorBoundaries(cnt1,
                                        cnt2,
                                        &min_k_real_m,
                                        &max_k_real_m,
                                        &min_l_real_m,
                                        &max_l_real_m
                                        );
              updatePosteriorBoundaries(cnt1,
                                        cnt2,
                                        &min_k_real_m1,
                                        &max_k_real_m1,
                                        &min_l_real_m1,
                                        &max_l_real_m1
                                        );
#ifdef COUNT_STATES
             vars->N_M[ij][cnt1][cnt2/2] = vars->N_M1[ij][cnt1][cnt2/2] = vars->N_C[ij][cnt1][cnt2/2];
#endif
            }
          }
        }

      /* 2nd E_M[ij] = MIN(E_M[ij], E_M[i+1,j] + c) */
      if(vars->E_M_rem[my_iindx[i+1]-j] != INF){
        vars->E_M_rem[ij] = MIN2(vars->E_M_rem[ij],
                                  vars->E_M_rem[my_iindx[i+1]-j] + P->MLbase
                                  );
      }
      if(vars->E_M[my_iindx[i+1]-j])
        for(cnt1 = vars->k_min_values_m[my_iindx[i+1]-j];
            cnt1 <= vars->k_max_values_m[my_iindx[i+1]-j];
            cnt1++){
          for(cnt2 = vars->l_min_values_m[my_iindx[i+1]-j][cnt1];
              cnt2 <= vars->l_max_values_m[my_iindx[i+1]-j][cnt1];
              cnt2+=2){
            if(vars->E_M[my_iindx[i+1]-j][cnt1][cnt2/2] != INF){
              if(((cnt1 + dia) <= maxD1) && ((cnt2 + dib) <= maxD2)){
                vars->E_M[ij][cnt1+dia][(cnt2+dib)/2] = MIN2( vars->E_M[ij][cnt1+dia][(cnt2+dib)/2],
                                                              vars->E_M[my_iindx[i+1]-j][cnt1][cnt2/2] + P->MLbase
                                                            );
                updatePosteriorBoundaries(cnt1 + dia,
                                          cnt2 + dib,
                                          &min_k_real_m,
                                          &max_k_real_m,
                                          &min_l_real_m,
                                          &max_l_real_m
                                          );
#ifdef COUNT_STATES
                vars->N_M[ij][cnt1+dia][(cnt2+dib)/2] += vars->N_M[my_iindx[i+1]-j][cnt1][cnt2/2];
#endif
              }
              /* collect all cases where dia+cnt1 or dib+cnt2 exceeds maxD1, maxD2, respectively */
              else{
                vars->E_M_rem[ij] = MIN2(vars->E_M_rem[ij],
                                          vars->E_M[my_iindx[i+1]-j][cnt1][cnt2/2] + P->MLbase
                                          );
              }
            }
          }
        }

      /* 3rd E_M[ij] = MIN(E_M[ij], E_M[i,j-1] + c) */
      if(vars->E_M_rem[ij+1] != INF){
        vars->E_M_rem[ij] = MIN2(vars->E_M_rem[ij],
                                  vars->E_M_rem[ij+1] + P->MLbase
                                  );
      }
      if(vars->E_M[ij+1])
        for(cnt1 = vars->k_min_values_m[ij+1];
            cnt1 <= vars->k_max_values_m[ij+1];
            cnt1++){
          for(cnt2 = vars->l_min_values_m[ij+1][cnt1];
              cnt2 <= vars->l_max_values_m[ij+1][cnt1];
              cnt2+=2){
            if(vars->E_M[ij+1][cnt1][cnt2/2] != INF){
              if(((cnt1 + dja) <= maxD1) && ((cnt2 + djb) <= maxD2)){
                vars->E_M[ij][cnt1+dja][(cnt2+djb)/2] = MIN2( vars->E_M[ij][cnt1+dja][(cnt2+djb)/2],
                                                              vars->E_M[ij+1][cnt1][cnt2/2] + P->MLbase
                                                            );
                updatePosteriorBoundaries(cnt1 + dja,
                                          cnt2 + djb,
                                          &min_k_real_m,
                                          &max_k_real_m,
                                          &min_l_real_m,
                                          &max_l_real_m
                                          );
#ifdef COUNT_STATES
                vars->N_M[ij][cnt1+dja][(cnt2+djb)/2] += vars->N_M[ij+1][cnt1][cnt2/2];
#endif
              }
              /* collect all cases where dja+cnt1 or djb+cnt2 exceeds maxD1, maxD2, respectively */
              else{
                vars->E_M_rem[ij] = MIN2(vars->E_M_rem[ij],
                                          vars->E_M[ij+1][cnt1][cnt2/2] + P->MLbase
                                          );
              }
            }
          }
        }

      /* 4th E_M1[ij] = MIN(E_M1[ij], E_M1[i,j-1] + c) */
      if(vars->E_M1_rem[ij+1] != INF){
        vars->E_M1_rem[ij] = MIN2( vars->E_M1_rem[ij],
                                    vars->E_M1_rem[ij+1] + P->MLbase
                                  );
      }
      if(vars->E_M1[ij+1])
        for(cnt1 = vars->k_min_values_m1[ij+1];
            cnt1 <= vars->k_max_values_m1[ij+1];
            cnt1++){
          for(cnt2 = vars->l_min_values_m1[ij+1][cnt1];
              cnt2 <= vars->l_max_values_m1[ij+1][cnt1];
              cnt2+=2){
            if(vars->E_M1[ij+1][cnt1][cnt2/2] != INF){
              if(((cnt1 + dja) <= maxD1) && ((cnt2 + djb) <= maxD2)){
                vars->E_M1[ij][cnt1+dja][(cnt2+djb)/2]  = MIN2( vars->E_M1[ij][cnt1+dja][(cnt2+djb)/2],
                                                                vars->E_M1[ij+1][cnt1][cnt2/2] + P->MLbase
                                                              );
                updatePosteriorBoundaries(cnt1 + dja,
                                          cnt2 + djb,
                                          &min_k_real_m1,
                                          &max_k_real_m1,
                                          &min_l_real_m1,
                                          &max_l_real_m1
                                          );
#ifdef COUNT_STATES
                vars->N_M1[ij][cnt1+dja][(cnt2+djb)/2]  += vars->N_M1[ij+1][cnt1][cnt2/2];
#endif
              }
              /* collect all cases where dja+cnt1 or djb+cnt2 exceeds maxD1, maxD2, respectively */
              else{
                vars->E_M1_rem[ij] = MIN2( vars->E_M1_rem[ij],
                                            vars->E_M1[ij+1][cnt1][cnt2/2] + P->MLbase
                                          );
              }
            }
          }
        }


      /* 5th E_M[ij] = MIN(E_M[ij], min(E_M[i,k] + E_M[k+1,j])) */
      if(j > TURN + 2)
      for (u = i+1+TURN; u <= j-2-TURN; u++){
        /* check all cases where M(i,u) and/or M(u+1,j) are already out of scope of maxD1 and/or maxD2 */
        if(vars->E_M_rem[my_iindx[i]-u] != INF){
          for(cnt3 = vars->k_min_values_m[my_iindx[u+1]-j];
              cnt3 <= vars->k_max_values_m[my_iindx[u+1]-j];
              cnt3++){
            for(cnt4 = vars->l_min_values_m[my_iindx[u+1]-j][cnt3];
                cnt4 <= vars->l_max_values_m[my_iindx[u+1]-j][cnt3];
                cnt4+=2){
              if(vars->E_M[my_iindx[u+1]-j][cnt3][cnt4/2] != INF){
                  vars->E_M_rem[ij] = MIN2(vars->E_M_rem[ij],
                                            vars->E_M_rem[my_iindx[i]-u] + vars->E_M[my_iindx[u+1]-j][cnt3][cnt4/2]
                                            );
              }
            }
          }
          if(vars->E_M_rem[my_iindx[u+1]-j] != INF){
            vars->E_M_rem[ij] = MIN2(vars->E_M_rem[ij],
                                      vars->E_M_rem[my_iindx[i]-u] + vars->E_M_rem[my_iindx[u+1]-j]
                                      );
          }
        }
        if(vars->E_M_rem[my_iindx[u+1]-j] != INF){
          for(cnt1 = vars->k_min_values_m[my_iindx[i]-u];
              cnt1 <= vars->k_max_values_m[my_iindx[i]-u];
              cnt1++){
            for(cnt2 = vars->l_min_values_m[my_iindx[i]-u][cnt1];
                cnt2 <= vars->l_max_values_m[my_iindx[i]-u][cnt1];
                cnt2+=2){
              if(vars->E_M[my_iindx[i]-u][cnt1][cnt2/2] != INF){
                vars->E_M_rem[ij] = MIN2(vars->E_M_rem[ij],
                                          vars->E_M[my_iindx[i]-u][cnt1][cnt2/2] + vars->E_M_rem[my_iindx[u+1]-j]
                                          );
              }
            }
          }
        }
        if(!vars->E_M[my_iindx[i]-u]) continue;
        if(!vars->E_M[my_iindx[u+1]-j]) continue;

        dia = referenceBPs1[ij] - referenceBPs1[my_iindx[i]-u] - referenceBPs1[my_iindx[u+1]-j];
        dib = referenceBPs2[ij] - referenceBPs2[my_iindx[i]-u] - referenceBPs2[my_iindx[u+1]-j];

        for(cnt1 = vars->k_min_values_m[my_iindx[i]-u];
            cnt1 <= vars->k_max_values_m[my_iindx[i]-u];
            cnt1++){
          for(cnt2 = vars->l_min_values_m[my_iindx[i]-u][cnt1];
              cnt2 <= vars->l_max_values_m[my_iindx[i]-u][cnt1];
              cnt2+=2){
            for(cnt3 = vars->k_min_values_m[my_iindx[u+1]-j];
                cnt3 <= vars->k_max_values_m[my_iindx[u+1]-j];
                cnt3++){
              for(cnt4 = vars->l_min_values_m[my_iindx[u+1]-j][cnt3];
                  cnt4 <= vars->l_max_values_m[my_iindx[u+1]-j][cnt3];
                  cnt4+=2){
                if((vars->E_M[my_iindx[i]-u][cnt1][cnt2/2] != INF) && (vars->E_M[my_iindx[u+1]-j][cnt3][cnt4/2] != INF)){
                  if(((cnt1 + cnt3 + dia) <= maxD1) && ((cnt2 + cnt4 + dib) <= maxD2)){
                    vars->E_M[ij][cnt1+cnt3+dia][(cnt2+cnt4+dib)/2] = MIN2( vars->E_M[ij][cnt1+cnt3+dia][(cnt2+cnt4+dib)/2],
                                                                            vars->E_M[my_iindx[i]-u][cnt1][cnt2/2]
                                                                          + vars->E_M[my_iindx[u+1]-j][cnt3][cnt4/2]
                                                                          );
                    updatePosteriorBoundaries(cnt1 + cnt3 + dia,
                                              cnt2 + cnt4 + dib,
                                              &min_k_real_m,
                                              &max_k_real_m,
                                              &min_l_real_m,
                                              &max_l_real_m
                                              );
#ifdef COUNT_STATES
                    vars->N_M[ij][cnt1+cnt3+dia][(cnt2+cnt4+dib)/2] += vars->N_M[my_iindx[i]-u][cnt1][cnt2/2] * vars->N_M1[my_iindx[u+1]-j][cnt3][cnt4/2];
#endif
                  }
                  /* collect all cases where dia+cnt1+cnt3 or dib+cnt2+cnt4 exceeds maxD1, maxD2, respectively */
                  else{
                    vars->E_M_rem[ij] = MIN2(vars->E_M_rem[ij],
                                              vars->E_M[my_iindx[i]-u][cnt1][cnt2/2] + vars->E_M[my_iindx[u+1]-j][cnt3][cnt4/2]
                                              );
                  }
                }
              }
            }
          }
        }
      }

      /* thats all folks for the multiloop decomposition... */

      adjustArrayBoundaries(&vars->E_M[ij],
                            &vars->k_min_values_m[ij],
                            &vars->k_max_values_m[ij],
                            &vars->l_min_values_m[ij],
                            &vars->l_max_values_m[ij],
                            min_k_real_m,
                            max_k_real_m,
                            min_l_real_m,
                            max_l_real_m
                            );

      adjustArrayBoundaries(&vars->E_M1[ij],
                            &vars->k_min_values_m1[ij],
                            &vars->k_max_values_m1[ij],
                            &vars->l_min_values_m1[ij],
                            &vars->l_max_values_m1[ij],
                            min_k_real_m1,
                            max_k_real_m1,
                            min_l_real_m1,
                            max_l_real_m1
                            );

#ifdef COUNT_STATES
        /* actually we should adjust the array boundaries here but we might never use the count states option more than once so what....*/
#endif
    } /* end of j-loop */
  }

  /* calculate energies of 5' and 3' fragments */

  /* prepare first entries in E_F5 */
  for(cnt1 = 1; cnt1 <= TURN+1; cnt1++){
    E_F5[cnt1] = (int **)space(sizeof(int *));
    E_F5[cnt1][0] = (int *)space(sizeof(int));
    E_F5[cnt1][0][0] = 0;
    vars->E_F5_rem[cnt1] = INF;
    k_min_values_f[cnt1] = k_max_values_f[cnt1] = 0;
    l_min_values_f[cnt1] = (int *)space(sizeof(int));
    l_max_values_f[cnt1] = (int *)space(sizeof(int));
    l_min_values_f[cnt1][0] = l_max_values_f[cnt1][0] = 0;
#ifdef COUNT_STATES
    vars->N_F5[cnt1] = (unsigned long **)space(sizeof(unsigned long *));
    vars->N_F5[cnt1][0] = (unsigned long *)space(sizeof(unsigned long));
    vars->N_F5[cnt1][0][0] = 1;
#endif

  }



  for (j=TURN+2; j <= seq_length; j++) {

    unsigned int da = referenceBPs1[my_iindx[1]-j] - referenceBPs1[my_iindx[1]-j+1];
    unsigned int db = referenceBPs2[my_iindx[1]-j] - referenceBPs2[my_iindx[1]-j+1];

    type=ptype[my_iindx[1]-j];
    additional_en = 0;
    if(type){
      if(dangles == 2)
        additional_en += E_ExtLoop(type, -1, j < seq_length ? S1[j+1] : -1, P);
      else
        additional_en += E_ExtLoop(type, -1, -1, P);
    }

    /* make min and max k guess for memory allocation */
    int min_k_guess, max_k_guess, min_l_guess, max_l_guess;
    int *min_l_real, *max_l_real, min_k_real, max_k_real;

    min_k_guess = min_l_guess = 0;
    max_k_guess = referenceBPs1[my_iindx[1]-j] + mm1[my_iindx[1]-j];
    max_l_guess = referenceBPs2[my_iindx[1]-j] + mm2[my_iindx[1]-j];

    prepareBoundaries(min_k_guess,
                      max_k_guess,
                      min_l_guess,
                      max_l_guess,
                      bpdist[my_iindx[1]-j],
                      &vars->k_min_values_f[j],
                      &vars->k_max_values_f[j],
                      &vars->l_min_values_f[j],
                      &vars->l_max_values_f[j]
                      );

    preparePosteriorBoundaries( vars->k_max_values_f[j] - vars->k_min_values_f[j] + 1,
                                vars->k_min_values_f[j],
                                &min_k_real,
                                &max_k_real,
                                &min_l_real,
                                &max_l_real
                              );

    prepareArray( &vars->E_F5[j],
                  vars->k_min_values_f[j],
                  vars->k_max_values_f[j],
                  vars->l_min_values_f[j],
                  vars->l_max_values_f[j]
                );
#ifdef COUNT_STATES
    prepareArray2( &vars->N_F5[j],
                  vars->k_min_values_f[j],
                  vars->k_max_values_f[j],
                  vars->l_min_values_f[j],
                  vars->l_max_values_f[j]
                );
#endif

    /* begin the actual computation of 5' end energies */

    /* j-1 is unpaired ... */
    vars->E_F5_rem[j] = vars->E_F5_rem[j-1];
    for(cnt1 = vars->k_min_values_f[j-1]; cnt1 <= vars->k_max_values_f[j-1]; cnt1++){
      for(cnt2 = vars->l_min_values_f[j-1][cnt1]; cnt2 <= vars->l_max_values_f[j-1][cnt1]; cnt2+=2){
        if(((cnt1 + da) <= maxD1) && ((cnt2 + db) <= maxD2)){
          vars->E_F5[j][cnt1+da][(cnt2+db)/2] = MIN2( vars->E_F5[j][cnt1+da][(cnt2+db)/2],
                                                      vars->E_F5[j-1][cnt1][cnt2/2]
                                                    );
          updatePosteriorBoundaries(cnt1 + da,
                                    cnt2 + db,
                                    &min_k_real,
                                    &max_k_real,
                                    &min_l_real,
                                    &max_l_real
                                    );
#ifdef COUNT_STATES
          vars->N_F5[j][cnt1+da][(cnt2+db)/2] += vars->N_F5[j-1][cnt1][cnt2/2];
#endif
        }
        /* collect all cases where da+cnt1 or db+cnt2 exceeds maxD1, maxD2, respectively */
        else{
          vars->E_F5_rem[j] = MIN2(vars->E_F5_rem[j], vars->E_F5[j-1][cnt1][cnt2/2]);
        }
      }
    }
    /* j pairs with 1 */
    if(vars->E_C_rem[my_iindx[1]-j] != INF){
      vars->E_F5_rem[j] = MIN2(vars->E_F5_rem[j], vars->E_C_rem[my_iindx[1]-j] + additional_en);
    }
    if(vars->E_C[my_iindx[1]-j])
      for(cnt1 = vars->k_min_values[my_iindx[1]-j]; cnt1 <= vars->k_max_values[my_iindx[1]-j]; cnt1++)
        for(cnt2 = vars->l_min_values[my_iindx[1]-j][cnt1]; cnt2 <= vars->l_max_values[my_iindx[1]-j][cnt1]; cnt2+=2){
          if(vars->E_C[my_iindx[1]-j][cnt1][cnt2/2] != INF){
            vars->E_F5[j][cnt1][cnt2/2] = MIN2( vars->E_F5[j][cnt1][cnt2/2],
                                                vars->E_C[my_iindx[1]-j][cnt1][cnt2/2]+ additional_en
                                              );
            updatePosteriorBoundaries(cnt1,
                                      cnt2,
                                      &min_k_real,
                                      &max_k_real,
                                      &min_l_real,
                                      &max_l_real
                                      );
#ifdef COUNT_STATES
            vars->N_F5[j][cnt1][cnt2/2] += vars->N_C[my_iindx[1]-j][cnt1][cnt2/2];
#endif
          }
        }

    /* j pairs with some other nucleotide -> see below */
    for (i=j-TURN-1; i>1; i--) {
      ij = my_iindx[i]-j;
      type = ptype[ij];
      if (type) {
        if(dangles == 2)
          additional_en = E_ExtLoop(type, S1[i-1], j < seq_length ? S1[j+1] : -1, P);
        else
          additional_en = E_ExtLoop(type, -1, -1, P);

        if(vars->E_C_rem[ij] != INF){
          for(cnt3 = vars->k_min_values_f[i-1]; cnt3 <= vars->k_max_values_f[i-1]; cnt3++)
            for(cnt4 = vars->l_min_values_f[i-1][cnt3]; cnt4 <= vars->l_max_values_f[i-1][cnt3]; cnt4+=2){
              if(vars->E_F5[i-1][cnt3][cnt4/2] != INF){
                vars->E_F5_rem[j] = MIN2(vars->E_F5_rem[j],
                                          vars->E_F5[i-1][cnt3][cnt4/2] + vars->E_C_rem[ij] + additional_en
                                          );
              }
            }
          if(vars->E_F5_rem[i-1] != INF){
            vars->E_F5_rem[j] = MIN2(vars->E_F5_rem[j],
                                      vars->E_F5_rem[i-1] + vars->E_C_rem[ij] + additional_en
                                      );
          }
        }
        if((vars->E_F5_rem[i-1] != INF) && (vars->E_C[ij])){
          for(cnt1 = vars->k_min_values[ij]; cnt1 <= vars->k_max_values[ij]; cnt1++)
            for(cnt2 = vars->l_min_values[ij][cnt1]; cnt2 <= vars->l_max_values[ij][cnt1]; cnt2+=2)
              if(vars->E_C[ij][cnt1][cnt2/2]!= INF){
                vars->E_F5_rem[j] = MIN2(vars->E_F5_rem[j],
                                          vars->E_F5_rem[i-1] + vars->E_C[ij][cnt1][cnt2/2] + additional_en
                                          );
              }
        }
        if(!vars->E_C[ij]) continue;

        unsigned int d1a = referenceBPs1[my_iindx[1]-j] - referenceBPs1[ij] - referenceBPs1[my_iindx[1]-i+1];
        unsigned int d1b = referenceBPs2[my_iindx[1]-j] - referenceBPs2[ij] - referenceBPs2[my_iindx[1]-i+1];

        for(cnt1 = vars->k_min_values[ij]; cnt1 <= vars->k_max_values[ij]; cnt1++)
          for(cnt2 = vars->l_min_values[ij][cnt1]; cnt2 <= vars->l_max_values[ij][cnt1]; cnt2+=2)
            for(cnt3 = vars->k_min_values_f[i-1]; cnt3 <= vars->k_max_values_f[i-1]; cnt3++)
              for(cnt4 = vars->l_min_values_f[i-1][cnt3]; cnt4 <= vars->l_max_values_f[i-1][cnt3]; cnt4+=2){
                if(vars->E_F5[i-1][cnt3][cnt4/2] != INF && vars->E_C[ij][cnt1][cnt2/2]!= INF){
                  if(((cnt1 + cnt3 + d1a) <= maxD1) && ((cnt2 + cnt4 + d1b) <= maxD2)){
                    vars->E_F5[j][cnt1+cnt3+d1a][(cnt2+cnt4+d1b)/2] = MIN2( vars->E_F5[j][cnt1+cnt3+d1a][(cnt2+cnt4+d1b)/2],
                                                                            vars->E_F5[i-1][cnt3][cnt4/2] + vars->E_C[ij][cnt1][cnt2/2] + additional_en
                                                                          );
                    updatePosteriorBoundaries(cnt1 + cnt3 + d1a,
                                              cnt2 + cnt4 + d1b,
                                              &min_k_real,
                                              &max_k_real,
                                              &min_l_real,
                                              &max_l_real
                                              );
#ifdef COUNT_STATES
                    vars->N_F5[j][cnt1+cnt3+d1a][(cnt2+cnt4+d1b)/2] += vars->N_F5[i-1][cnt3][cnt4/2] * vars->N_C[ij][cnt1][cnt2/2];
#endif
                  }
                  /* collect all cases where d1a+cnt1+cnt3 or d1b+cnt2+cnt4 exceeds maxD1, maxD2, respectively */
                  else{
                    vars->E_F5_rem[j] = MIN2(vars->E_F5_rem[j],
                                              vars->E_F5[i-1][cnt3][cnt4/2] + vars->E_C[ij][cnt1][cnt2/2] + additional_en
                                              );
                  }
                }
              }
      }
    }

    /* resize and move memory portions of energy matrix E_F5 */
    adjustArrayBoundaries(&vars->E_F5[j],
                          &vars->k_min_values_f[j],
                          &vars->k_max_values_f[j],
                          &vars->l_min_values_f[j],
                          &vars->l_max_values_f[j],
                          min_k_real,
                          max_k_real,
                          min_l_real,
                          max_l_real
                          );

  } /* end of j-loop */


  if(compute_2Dfold_F3){

    /* prepare first entries in E_F3 */
    for(cnt1 = seq_length; cnt1 >= seq_length-TURN-1; cnt1--){
      E_F3[cnt1]        = (int **)space(sizeof(int *));
      E_F3[cnt1][0]     = (int *) space(sizeof(int));
      E_F3[cnt1][0][0]  = 0;
      k_min_values_f3[cnt1]     = k_max_values_f3[cnt1] = 0;
      l_min_values_f3[cnt1]     = (int *)space(sizeof(int));
      l_max_values_f3[cnt1]     = (int *)space(sizeof(int));
      l_min_values_f3[cnt1][0]  = l_max_values_f3[cnt1][0] = 0;
    }
    /* begin calculations */
    for (j=seq_length-TURN-2; j >= 1; j--){

      unsigned int da = referenceBPs1[my_iindx[j]-seq_length] - referenceBPs1[my_iindx[j+1]-seq_length];
      unsigned int db = referenceBPs2[my_iindx[j]-seq_length] - referenceBPs2[my_iindx[j+1]-seq_length];

      type=ptype[my_iindx[j]-seq_length];
      additional_en = 0;
      if(type){
        if(dangles == 2)
          additional_en += E_ExtLoop(type, j > 1 ? S1[j-1] : -1, -1, P);
        else
          additional_en += E_ExtLoop(type, -1, -1, P);
      }

      /* make min and max k guess for memory allocation */
      int min_k_guess, max_k_guess, min_l_guess, max_l_guess;
      int *min_l_real, *max_l_real, min_k_real, max_k_real;

      min_k_guess = min_l_guess = 0;
      max_k_guess = referenceBPs1[my_iindx[j]-seq_length] + mm1[my_iindx[j]-seq_length];
      max_l_guess = referenceBPs2[my_iindx[j]-seq_length] + mm2[my_iindx[j]-seq_length];

      prepareBoundaries(min_k_guess,
                        max_k_guess,
                        min_l_guess,
                        max_l_guess,
                        bpdist[my_iindx[j]-seq_length],
                        &vars->k_min_values_f3[j],
                        &vars->k_max_values_f3[j],
                        &vars->l_min_values_f3[j],
                        &vars->l_max_values_f3[j]
                        );

      preparePosteriorBoundaries( vars->k_max_values_f3[j] - vars->k_min_values_f3[j] + 1,
                                  vars->k_min_values_f3[j],
                                  &min_k_real,
                                  &max_k_real,
                                  &min_l_real,
                                  &max_l_real
                                );

      prepareArray( &vars->E_F3[j],
                    vars->k_min_values_f3[j],
                    vars->k_max_values_f3[j],
                    vars->l_min_values_f3[j],
                    vars->l_max_values_f3[j]
                  );
      /* begin the actual computation of 5' end energies */

      /* j is unpaired ... */
      for(cnt1 = vars->k_min_values_f3[j+1]; cnt1 <= vars->k_max_values_f3[j+1]; cnt1++){
        for(cnt2 = vars->l_min_values_f3[j+1][cnt1]; cnt2 <= vars->l_max_values_f3[j+1][cnt1]; cnt2+=2){
          vars->E_F3[j][cnt1+da][(cnt2+db)/2] = MIN2( vars->E_F3[j][cnt1+da][(cnt2+db)/2],
                                                      vars->E_F3[j+1][cnt1][cnt2/2]
                                                    );
          updatePosteriorBoundaries(cnt1 + da,
                                    cnt2 + db,
                                    &min_k_real,
                                    &max_k_real,
                                    &min_l_real,
                                    &max_l_real
                                    );
        }
      }
      /* j pairs with n */
      if(vars->E_C[my_iindx[j]-seq_length])
        for(cnt1 = vars->k_min_values[my_iindx[j]-seq_length]; cnt1 <= vars->k_max_values[my_iindx[j]-seq_length]; cnt1++)
          for(cnt2 = vars->l_min_values[my_iindx[j]-seq_length][cnt1]; cnt2 <= vars->l_max_values[my_iindx[j]-seq_length][cnt1]; cnt2+=2){
            if(vars->E_C[my_iindx[j]-seq_length][cnt1][cnt2/2] != INF){
              vars->E_F3[j][cnt1][cnt2/2] = MIN2( vars->E_F3[j][cnt1][cnt2/2],
                                                  vars->E_C[my_iindx[j]-seq_length][cnt1][cnt2/2]+ additional_en
                                                );
              updatePosteriorBoundaries(cnt1,
                                        cnt2,
                                        &min_k_real,
                                        &max_k_real,
                                        &min_l_real,
                                        &max_l_real
                                        );
            }
          }

      /* j pairs with some other nucleotide -> see below */
      for (i=j-TURN-1; i>1; i--) {
        ij = my_iindx[i]-j;
        if(!vars->E_C[ij]) continue;
        type = ptype[ij];
        if (type) {
          unsigned int d1a = referenceBPs1[my_iindx[1]-j] - referenceBPs1[ij] - referenceBPs1[my_iindx[1]-i+1];
          unsigned int d1b = referenceBPs2[my_iindx[1]-j] - referenceBPs2[ij] - referenceBPs2[my_iindx[1]-i+1];

          if(dangles == 2)
            additional_en = E_ExtLoop(type, S1[i-1], j < seq_length ? S1[j+1] : -1, P);
          else
            additional_en = E_ExtLoop(type, -1, -1, P);

          for(cnt1 = vars->k_min_values[ij]; cnt1 <= vars->k_max_values[ij]; cnt1++)
            for(cnt2 = vars->l_min_values[ij][cnt1]; cnt2 <= vars->l_max_values[ij][cnt1]; cnt2+=2)
              for(cnt3 = vars->k_min_values_f[i-1]; cnt3 <= vars->k_max_values_f[i-1]; cnt3++)
                for(cnt4 = vars->l_min_values_f[i-1][cnt3]; cnt4 <= vars->l_max_values_f[i-1][cnt3]; cnt4+=2){
                  if(vars->E_F5[i-1][cnt3][cnt4/2] != INF && vars->E_C[ij][cnt1][cnt2/2]!= INF){
                    vars->E_F5[j][cnt1+cnt3+d1a][(cnt2+cnt4+d1b)/2] = MIN2( vars->E_F5[j][cnt1+cnt3+d1a][(cnt2+cnt4+d1b)/2],
                                                                            vars->E_F5[i-1][cnt3][cnt4/2] + vars->E_C[ij][cnt1][cnt2/2] + additional_en
                                                                          );
                    updatePosteriorBoundaries(cnt1 + cnt3 + d1a,
                                              cnt2 + cnt4 + d1b,
                                              &min_k_real,
                                              &max_k_real,
                                              &min_l_real,
                                              &max_l_real
                                              );
#ifdef COUNT_STATES
                    vars->N_F5[j][cnt1+cnt3+d1a][(cnt2+cnt4+d1b)/2] += vars->N_F5[i-1][cnt3][cnt4/2] * vars->N_C[ij][cnt1][cnt2/2];
#endif
                  }
                }
        }
      }

      /* resize and move memory portions of energy matrix E_F5 */
      adjustArrayBoundaries(&vars->E_F5[j],
                            &vars->k_min_values_f[j],
                            &vars->k_max_values_f[j],
                            &vars->l_min_values_f[j],
                            &vars->l_max_values_f[j],
                            min_k_real,
                            max_k_real,
                            min_l_real,
                            max_l_real
                            );

    } /* end of j-loop */




  }
}


/*---------------------------------------------------------------------------*/

PUBLIC void update_TwoDfold_params(TwoDfold_vars *vars){
  if(vars->P) free(vars->P);
  vars->P = scale_parameters();
  make_pair_matrix();
}

/*---------------------------------------------------------------------------*/
PRIVATE void make_ptypes(TwoDfold_vars *vars) {
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

PRIVATE void backtrack_f5(unsigned int j, int k, int l, char *structure, TwoDfold_vars *vars){
  int           *my_iindx, energy, type, dangles, cnt1, cnt2, cnt3, cnt4;
  int           **l_min_values, **l_max_values,**l_min_values_f, **l_max_values_f;
  int           *k_min_values, *k_max_values,*k_min_values_f, *k_max_values_f;
  int           ***E_C, ***E_F5;
  int           *E_C_rem, *E_F5_rem;
  unsigned int   i, ij, seq_length, maxD1, maxD2;
  short *S1;
  unsigned int   *referenceBPs1, *referenceBPs2;
  char  *ptype;
  paramT   *P;
  unsigned int da, db;

  P               = vars->P;
  seq_length      = vars->seq_length;
  S1              = vars->S1;
  ptype           = vars->ptype;
  my_iindx        = vars->my_iindx;
  referenceBPs1   = vars->referenceBPs1;
  referenceBPs2   = vars->referenceBPs2;
  dangles         = vars->dangles;
  E_F5            = vars->E_F5;
  l_min_values_f  = vars->l_min_values_f;
  l_max_values_f  = vars->l_max_values_f;
  k_min_values_f  = vars->k_min_values_f;
  k_max_values_f  = vars->k_max_values_f;

  E_C             = vars->E_C;
  l_min_values    = vars->l_min_values;
  l_max_values    = vars->l_max_values;
  k_min_values    = vars->k_min_values;
  k_max_values    = vars->k_max_values;

  E_F5_rem       = vars->E_F5_rem;
  E_C_rem        = vars->E_C_rem;
  maxD1           = vars->maxD1;
  maxD2           = vars->maxD2;

  da = referenceBPs1[my_iindx[1]-j] - referenceBPs1[my_iindx[1]-j+1];
  db = referenceBPs2[my_iindx[1]-j] - referenceBPs2[my_iindx[1]-j+1];

  if(j<TURN+2) return;

  /* F5[j] == F5[j-1] ? */
  if(k == -1){
    if(E_F5_rem[j]==INF)
      return;
    else if(E_F5_rem[j] == E_F5_rem[j-1]){
      backtrack_f5(j-1,k,l,structure, vars);
      return;
    }
    else if(E_F5[j-1]){
      for(cnt1 =  k_min_values_f[j-1];
        cnt1 <= k_max_values_f[j-1];
        cnt1++){
        for(cnt2 = l_min_values_f[j-1][cnt1];
            cnt2 <= l_max_values_f[j-1][cnt1];
            cnt2+=2){
          if(((cnt1 + da) > maxD1) || ((cnt2 + db) > maxD2)){
            if(E_F5_rem[j] == E_F5[j-1][cnt1][cnt2/2]){
              backtrack_f5(j-1, cnt1, cnt2, structure, vars);
              return;
            }
          }
        }
      }
    }
  }
  else if((k >= da) && (l >= db)){
    if(E_F5[j-1]){
      if((k - da >= k_min_values_f[j-1]) && (k - da <= k_max_values_f[j-1])){
        if((l - db >= l_min_values_f[j-1][k-da]) && (l - db <= l_max_values_f[j-1][k-da]))
          if(E_F5[j-1][k-da][(l-db)/2] == E_F5[j][k][l/2]){
            backtrack_f5(j-1, k-da, l-db, structure, vars);
            return;
          }
      }
    }
  }

  type = ptype[my_iindx[1]-j];
  if(type){
    if(dangles == 2)
      energy = E_ExtLoop(type, -1, j < seq_length ? S1[j+1] : -1, P);
    else
      energy = E_ExtLoop(type, -1, -1, P);

    if(k == -1){
      if(E_C_rem[my_iindx[1]-j] + energy == E_F5_rem[j]){
          backtrack_c(1, j, -1, -1, structure, vars);
          return;
      }
    }
    else if(k >= k_min_values[my_iindx[1]-j] && (k <= k_max_values[my_iindx[1]-j])){

      if((l >= l_min_values[my_iindx[1]-j][k]) && (l <= l_max_values[my_iindx[1]-j][k]))
        if(E_C[my_iindx[1]-j][k][l/2] + energy == E_F5[j][k][l/2]){
          backtrack_c(1, j, k, l, structure, vars);
          return;
        }
    }
  }

  for (i=j-TURN-1; i>1; i--) {
    ij = my_iindx[i]-j;
    type = ptype[ij];
    if (type) {
      unsigned int d1a = referenceBPs1[my_iindx[1]-j] - referenceBPs1[ij] - referenceBPs1[my_iindx[1]-i+1];
      unsigned int d1b = referenceBPs2[my_iindx[1]-j] - referenceBPs2[ij] - referenceBPs2[my_iindx[1]-i+1];

      if(dangles == 2)
        energy = E_ExtLoop(type, S1[i-1], j < seq_length ? S1[j+1] : -1, P);
      else
        energy = E_ExtLoop(type, -1, -1, P);

      if(k == -1){
        if(E_C_rem[ij] != INF){
          for(cnt1 = k_min_values_f[i-1];
              cnt1 <= k_max_values_f[i-1];
              cnt1++){
            for(cnt2 = l_min_values_f[i-1][cnt1];
                cnt2 <= l_max_values_f[i-1][cnt1];
                cnt2+=2){
              if(E_F5_rem[j] == (E_F5[i-1][cnt1][cnt2/2] + E_C_rem[ij] + energy)){
                backtrack_f5(i-1, cnt1, cnt2, structure, vars);
                backtrack_c(i,j,-1,-1,structure, vars);
                return;
              }
            }
          }
          if(E_F5_rem[j] == (E_F5_rem[i-1] + E_C_rem[ij] + energy)){
            backtrack_f5(i-1, -1, -1, structure, vars);
            backtrack_c(i,j,-1,-1,structure,vars);
            return;
          }
        }
        if(E_F5_rem[i-1] != INF){
          for(cnt1 = k_min_values[ij];
              cnt1 <= k_max_values[ij];
              cnt1++){
            for(cnt2 = l_min_values[ij][cnt1];
                cnt2 <= l_max_values[ij][cnt1];
                cnt2 += 2){
               if(E_F5_rem[j] == (E_F5_rem[i-1] + E_C[ij][cnt1][cnt2/2] + energy)){
                backtrack_f5(i-1,-1,-1,structure,vars);
                backtrack_c(i,j,cnt1,cnt2,structure,vars);
                return;
              }
            }
          }
        }
        for(cnt1 = k_min_values_f[i-1];
            cnt1 <= k_max_values_f[i-1];
            cnt1++)
          for(cnt2 = l_min_values_f[i-1][cnt1];
              cnt2 <= l_max_values_f[i-1][cnt1];
              cnt2 += 2)
            for(cnt3 = k_min_values[ij];
                cnt3 <= k_max_values[ij];
                cnt3++)
              for(cnt4 = l_min_values[ij][cnt3];
                  cnt4 <= l_max_values[ij][cnt3];
                  cnt4 += 2){
                if(((cnt1 + cnt3 + d1a)>maxD1) || ((cnt2+cnt4+d1b)>maxD2)){
                  if(E_F5_rem[j] == (E_F5[i-1][cnt1][cnt2/2] + E_C[ij][cnt3][cnt4/2] + energy)){
                    backtrack_f5(i-1,cnt1,cnt2,structure,vars);
                    backtrack_c(i,j,cnt3,cnt4,structure,vars);
                    return;
                  }
                }
              }
      }
      else if((k >= d1a) && (l >= d1b)){
        int k_f_max = MIN2(k-d1a, k_max_values_f[i-1]);

        for(cnt1 = k_min_values_f[i-1]; cnt1 <= k_f_max; cnt1++){
          int l_f_max = MIN2(l - d1b, l_max_values_f[i-1][cnt1]);
          for(cnt2 = l_min_values_f[i-1][cnt1]; cnt2 <= l_f_max; cnt2+=2){
            int k_c = k - d1a - cnt1;
            if((k_c >= k_min_values[ij]) && (k_c <= k_max_values[ij])){
              int l_c = l - d1b - cnt2;
              if((l_c >= l_min_values[ij][k_c]) && (l_c <= l_max_values[ij][k_c])){
                if(E_F5[j][k][l/2] == (E_F5[i-1][cnt1][cnt2/2] + E_C[ij][k_c][l_c/2] + energy)){
                  backtrack_f5(i-1, cnt1, cnt2, structure, vars);
                  backtrack_c(i, j, k_c, l_c, structure, vars);
                  return;
                }
              }
            }
          }

        }
      }
    }
  }
  nrerror("backtracking failed in f5");
}

PRIVATE void backtrack_c(unsigned int i, unsigned int j, int k, int l, char *structure, TwoDfold_vars *vars){
  unsigned int p, q, pq, ij, maxp, seq_length, maxD1, maxD2;
  int *my_iindx, type, type_2, energy, no_close, dangles, base_d1, base_d2, d1, d2, cnt1, cnt2, cnt3, cnt4;
  int           **l_min_values, **l_max_values,**l_min_values_m, **l_max_values_m,**l_min_values_m1, **l_max_values_m1;
  int           *k_min_values, *k_max_values,*k_min_values_m, *k_max_values_m,*k_min_values_m1, *k_max_values_m1;
  int           ***E_C, ***E_M, ***E_M1, *E_C_rem, *E_M_rem, *E_M1_rem;
  short *S1;
  unsigned int   *referenceBPs1, *referenceBPs2;
  char  *ptype, *sequence;
  paramT   *P;

  P               = vars->P;
  sequence        = vars->sequence;
  seq_length      = vars->seq_length;
  S1              = vars->S1;
  ptype           = vars->ptype;
  my_iindx        = vars->my_iindx;
  referenceBPs1   = vars->referenceBPs1;
  referenceBPs2   = vars->referenceBPs2;
  dangles         = vars->dangles;

  E_C             = vars->E_C;
  l_min_values    = vars->l_min_values;
  l_max_values    = vars->l_max_values;
  k_min_values    = vars->k_min_values;
  k_max_values    = vars->k_max_values;

  E_M             = vars->E_M;
  l_min_values_m  = vars->l_min_values_m;
  l_max_values_m  = vars->l_max_values_m;
  k_min_values_m  = vars->k_min_values_m;
  k_max_values_m  = vars->k_max_values_m;

  E_M1            = vars->E_M1;
  l_min_values_m1 = vars->l_min_values_m1;
  l_max_values_m1 = vars->l_max_values_m1;
  k_min_values_m1 = vars->k_min_values_m1;
  k_max_values_m1 = vars->k_max_values_m1;

  E_C_rem        = vars->E_C_rem;
  E_M_rem        = vars->E_M_rem;
  E_M1_rem       = vars->E_M1_rem;
  maxD1           = vars->maxD1;
  maxD2           = vars->maxD2;


  ij = my_iindx[i]-j;

  int e = (k==-1) ? E_C_rem[ij] : E_C[ij][k][l/2];

  type = ptype[ij];

  no_close = (((type==3)||(type==4))&&no_closingGU);
  structure[i-1] = '(';
  structure[j-1] = ')';

  base_d1 = ((unsigned int)vars->reference_pt1[i] != j) ? 1 : -1;
  base_d2 = ((unsigned int)vars->reference_pt2[i] != j) ? 1 : -1;

  base_d1 += referenceBPs1[ij];
  base_d2 += referenceBPs2[ij];

  if(k == -1){
    if(((unsigned int)base_d1 > maxD1) || ((unsigned int)base_d2 > maxD2)){
      if(e == E_Hairpin(j-i-1, type, S1[i+1], S1[j-1], sequence+i-1, P)) return;
    }
  }
  else{
    if((unsigned int)base_d1 == k)
      if((unsigned int)base_d2 == l)
        if(E_Hairpin(j-i-1, type, S1[i+1], S1[j-1], sequence+i-1, P) == e) return;
  }
  maxp = MIN2(j-2-TURN,i+MAXLOOP+1);
  for(p = i+1; p <= maxp; p++){
    unsigned int minq, ln_pre;
    minq = p + TURN + 1;
    ln_pre = j - i - 1;
    if(ln_pre > minq + MAXLOOP) minq = ln_pre - MAXLOOP - 1;
    for (q = minq; q < j; q++) {
      pq = my_iindx[p]-q;
      type_2 = ptype[pq];
      if (type_2==0) continue;
      type_2 = rtype[type_2];

      /* d2 = dbp(S_{i,j}, S_{p.q} + {i,j}) */
      d1 = base_d1 - referenceBPs1[pq];
      d2 = base_d2 - referenceBPs2[pq];

      energy = E_IntLoop(p-i-1, j-q-1, type, type_2, S1[i+1], S1[j-1], S1[p-1], S1[q+1], P);


      if(k == -1){
        if(E_C_rem[pq] != INF)
          if(e == (E_C_rem[pq] + energy)){
            backtrack_c(p,q,-1,-1,structure,vars);
            return;
          }
        if(E_C[pq])
        for(cnt1 = k_min_values[pq];
            cnt1 <= k_max_values[pq];
            cnt1++)
          for(cnt2 = l_min_values[pq][cnt1];
              cnt2 <= l_max_values[pq][cnt1];
              cnt2 += 2){
            if(((cnt1 + d1) > maxD1) || ((cnt2 + d2) > maxD2)){
              if(e == (E_C[pq][cnt1][cnt2/2] + energy)){
                backtrack_c(p,q,cnt1,cnt2,structure,vars);
                return;
              }
            }
          }
      }
      else{
        if(!E_C[pq]) continue;
        if(d1 <= k && d2 <= l){
          if((k-d1 >= k_min_values[pq]) && (k-d1) <= k_max_values[pq])
            if((l - d2 >= l_min_values[pq][k-d1]) && (l-d2 <= l_max_values[pq][k-d1]))
              if(E_C[pq][k-d1][(l-d2)/2] + energy == e){
                backtrack_c(p, q, k-d1, l-d2, structure, vars);
                return;
              }
        }
      }
    } /* end q-loop */
  } /* end p-loop */

  /* multi-loop decomposition ------------------------*/
  if(!no_close){
    unsigned int u;
    int tt;
    if(k==-1){
      for(u=i+TURN+2; u<j-TURN-2;u++){
        int i1u, u1j1;
        i1u   = my_iindx[i+1]-u;
        u1j1  = my_iindx[u+1]-j+1;
        tt = rtype[type];
        energy = P->MLclosing;
        if(dangles == 2)
          energy += E_MLstem(tt, S1[j-1], S1[i+1], P);
        else
          energy += E_MLstem(tt, -1, -1, P);


        if(E_M_rem[i1u] != INF){
          if(E_M1[u1j1])
          for(cnt1 = k_min_values_m1[u1j1];
              cnt1 <= k_max_values_m1[u1j1];
              cnt1++)
            for(cnt2 = l_min_values_m1[u1j1][cnt1];
                cnt2 <= l_max_values_m1[u1j1][cnt1];
                cnt2 += 2){
              if(e == (E_M_rem[i1u] + E_M1[u1j1][cnt1][cnt2/2] + energy)){
                backtrack_m(i+1,u,-1,-1,structure,vars);
                backtrack_m1(u+1,j-1,cnt1,cnt2,structure,vars);
                return;
              }
            }
          if(E_M1_rem[u1j1] != INF){
            if(e == (E_M_rem[i1u] + E_M1_rem[u1j1] + energy)){
              backtrack_m(i+1, u, -1, -1, structure, vars);
              backtrack_m1(u+1, j-1, -1, -1, structure, vars);
              return;
            }
          }
        }
        if(E_M1_rem[u1j1] != INF){
          if(E_M[i1u])
          for(cnt1 = k_min_values_m[i1u];
              cnt1 <= k_max_values_m[i1u];
              cnt1++)
            for(cnt2 = l_min_values_m[i1u][cnt1];
                cnt2 <= l_max_values_m[i1u][cnt1];
                cnt2 += 2)
              if(e == (E_M[i1u][cnt1][cnt2/2] + E_M1_rem[u1j1] + energy)){
                backtrack_m(i+1,u,cnt1,cnt2,structure,vars);
                backtrack_m1(u+1,j-1,-1,-1,structure,vars);
                return;
              }
        }

        /* now all cases where we exceed the maxD1/D2 scope by combination of E_M and E_M1 */
        if(!E_M[i1u]) continue;
        if(!E_M1[u1j1]) continue;
        /* get distance to reference if closing this multiloop
        *  dist3 = dbp(S_{i,j}, {i,j} + S_{i+1.u} + S_{u+1,j-1})
        */
        d1 = base_d1 - referenceBPs1[i1u] - referenceBPs1[u1j1];
        d2 = base_d2 - referenceBPs2[i1u] - referenceBPs2[u1j1];
        
        for(cnt1 = vars->k_min_values_m[i1u];
            cnt1 <= vars->k_max_values_m[i1u];
            cnt1++)
          for(cnt2 = vars->l_min_values_m[i1u][cnt1];
              cnt2 <= vars->l_max_values_m[i1u][cnt1];
              cnt2+=2)
            for(cnt3 = vars->k_min_values_m1[u1j1];
                cnt3 <= vars->k_max_values_m1[u1j1];
                cnt3++)
              for(cnt4 = vars->l_min_values_m1[u1j1][cnt3];
                  cnt4 <= vars->l_max_values_m1[u1j1][cnt3];
                  cnt4+=2){
                if(((cnt1 + cnt3 + d1) > maxD1) || ((cnt2 + cnt4 + d2) > maxD2)){
                  if(e == (E_M[i1u][cnt1][cnt2/2] + E_M1[u1j1][cnt3][cnt4/2] + energy)){
                    backtrack_m(i+1,u,cnt1,cnt2,structure,vars);
                    backtrack_m1(u+1,j-1,cnt3,cnt4,structure,vars);
                    return;
                  }
                }
              }
      }
    }
    else{
      for(u=i+TURN+2; u<j-TURN-2;u++){
        int i1u, u1j1;
        i1u   = my_iindx[i+1]-u;
        u1j1  = my_iindx[u+1]-j+1;
        if(!E_M[i1u]) continue;
        if(!E_M1[u1j1]) continue;

        /* get distance to reference if closing this multiloop
        *  dist3 = dbp(S_{i,j}, {i,j} + S_{i+1.u} + S_{u+1,j-1})
        */
        d1 = base_d1 - referenceBPs1[i1u] - referenceBPs1[u1j1];
        d2 = base_d2 - referenceBPs2[i1u] - referenceBPs2[u1j1];

        tt = rtype[type];
        energy = P->MLclosing;
        if(dangles == 2)
          energy += E_MLstem(tt, S1[j-1], S1[i+1], P);
        else
          energy += E_MLstem(tt, -1, -1, P);

        if((d1 <= k) && (d2 <= l))
          for(cnt1 = k_min_values_m[i1u];
              cnt1 <= MIN2(k-d1, k_max_values_m[i1u]);
              cnt1++)
            for(cnt2 = l_min_values_m[i1u][cnt1];
                cnt2 <= MIN2(l-d2, l_max_values_m[i1u][cnt1]);
                cnt2+=2)
              if(     ((k-d1-cnt1) >= k_min_values_m1[u1j1])
                  &&  ((k-d1-cnt1) <= k_max_values_m1[u1j1]))
                if(     ((l-d2-cnt2) >= l_min_values_m1[u1j1][k-d1-cnt1])
                    &&  ((l-d2-cnt2) <= l_max_values_m1[u1j1][k-d1-cnt1]))
                  if(e == (energy + E_M[i1u][cnt1][cnt2/2] + E_M1[u1j1][k-d1-cnt1][(l-d2-cnt2)/2])){
                    backtrack_m(i+1, u, cnt1, cnt2, structure, vars);
                    backtrack_m1(u+1, j-1, k-d1-cnt1, l-d2-cnt2, structure, vars);
                    return;
                  }
      }
    }
  }
  nrerror("backtracking failed in c");
}

PRIVATE void backtrack_m(unsigned int i, unsigned int j, int k, int l, char *structure, TwoDfold_vars *vars){
  unsigned int u, ij, seq_length, base_d1, base_d2, d1, d2, maxD1, maxD2;
  int *my_iindx, type, energy, dangles,circ, cnt1, cnt2, cnt3, cnt4;
  int           **l_min_values, **l_max_values,**l_min_values_m, **l_max_values_m;
  int           *k_min_values, *k_max_values,*k_min_values_m, *k_max_values_m;
  int           ***E_C, ***E_M, *E_C_rem, *E_M_rem;
  short *S1;
  unsigned int   *referenceBPs1, *referenceBPs2;
  char  *ptype, *sequence;
  paramT   *P;

  P           = vars->P;
  sequence    = vars->sequence;
  seq_length  = vars->seq_length;
  S1          = vars->S1;
  circ        = vars->circ;
  ptype       = vars->ptype;
  my_iindx    = vars->my_iindx;
  referenceBPs1  = vars->referenceBPs1;
  referenceBPs2  = vars->referenceBPs2;
  dangles     = vars->dangles;

  E_C             = vars->E_C;
  l_min_values    = vars->l_min_values;
  l_max_values    = vars->l_max_values;
  k_min_values    = vars->k_min_values;
  k_max_values    = vars->k_max_values;

  E_M             = vars->E_M;
  l_min_values_m  = vars->l_min_values_m;
  l_max_values_m  = vars->l_max_values_m;
  k_min_values_m  = vars->k_min_values_m;
  k_max_values_m  = vars->k_max_values_m;

  E_C_rem        = vars->E_C_rem;
  E_M_rem        = vars->E_M_rem;
  maxD1           = vars->maxD1;
  maxD2           = vars->maxD2;

  ij = my_iindx[i]-j;
  int e = (k == -1) ? E_M_rem[ij] : E_M[ij][k][l/2];

  base_d1 = referenceBPs1[ij];
  base_d2 = referenceBPs2[ij];

  if(k == -1){
    /* new_fML = ML(i+1,j)+c */
    d1 = base_d1 - referenceBPs1[my_iindx[i+1]-j];
    d2 = base_d2 - referenceBPs2[my_iindx[i+1]-j];
    if(E_M_rem[my_iindx[i+1]-j] != INF){
      if(e == (E_M_rem[my_iindx[i+1]-j] + P->MLbase)){
        backtrack_m(i+1,j,-1,-1,structure,vars);
        return;
      }
    }
    if(E_M[my_iindx[i+1]-j])
    for(cnt1 = k_min_values_m[my_iindx[i+1]-j];
        cnt1 <= k_max_values_m[my_iindx[i+1]-j];
        cnt1++)
      for(cnt2 = l_min_values_m[my_iindx[i+1]-j][cnt1];
          cnt2 <= l_max_values_m[my_iindx[i+1]-j][cnt1];
          cnt2 += 2)
        if(((cnt1 + d1) > maxD1) || ((cnt2 + d2) > maxD2)){
          if(e == (E_M[my_iindx[i+1]-j][cnt1][cnt2/2] + P->MLbase)){
            backtrack_m(i+1,j,cnt1,cnt2,structure,vars);
            return;
          }
        }

    /* new_fML = min(ML(i,j-1) + c, new_fML) */
    d1 = base_d1 - referenceBPs1[ij+1];
    d2 = base_d2 - referenceBPs2[ij+1];
    if(E_M_rem[ij+1] != INF){
      if(e == (E_M_rem[ij+1] + P->MLbase)){
        backtrack_m(i,j-1,-1,-1,structure,vars);
        return;
      }
    }
    if(E_M[ij+1])
    for(cnt1 = k_min_values_m[ij+1];
        cnt1 <= k_max_values_m[ij+1];
        cnt1++)
      for(cnt2 = l_min_values_m[ij+1][cnt1];
          cnt2 <= l_max_values_m[ij+1][cnt1];
          cnt2 += 2)
        if(((cnt1 + d1) > maxD1) || ((cnt2 + d2) > maxD2)){
          if(e == (E_M[ij+1][cnt1][cnt2/2] + P->MLbase)){
            backtrack_m(i,j-1,cnt1,cnt2,structure,vars);
            return;
          }
        }

    /* new_fML = min(new_fML, C(i,j)+b) */
    if(E_C_rem[ij] != INF){
      type = ptype[ij];
      if(dangles == 2)
        energy = E_MLstem(type, ((i > 1) || circ) ? S1[i-1] : -1, ((j < seq_length) || circ) ? S1[j+1] : -1, P);
      else
        energy = E_MLstem(type, -1, -1, P);
      if(e == (E_C_rem[ij] + energy)){
        backtrack_c(i,j,-1,-1,structure,vars);
        return;
      }
    }

    /* modular decomposition -------------------------------*/
    for(u = i+1+TURN; u <= j-2-TURN; u++){
      int iu, uj;
      iu = my_iindx[i]-u;
      uj = my_iindx[u+1]-j;
      type = ptype[uj];

      d1 = base_d1 - referenceBPs1[iu] - referenceBPs1[uj];
      d2 = base_d2 - referenceBPs2[iu] - referenceBPs2[uj];

      if(dangles == 2)
        energy = E_MLstem(type, S1[u], (j < seq_length) || circ ? S1[j+1] : -1, P);
      else
        energy = E_MLstem(type, -1, -1, P);

      if(E_M_rem[iu] != INF){
        if(E_C[uj])
        for(cnt1 = k_min_values[uj];
            cnt1 <= k_max_values[uj];
            cnt1++)
          for(cnt2 = l_min_values[uj][cnt1];
              cnt2 <= l_max_values[uj][cnt1];
              cnt2 += 2)
            if(e == (E_M_rem[iu] + E_C[uj][cnt1][cnt2/2] + energy)){
              backtrack_m(i,u,-1,-1,structure,vars);
              backtrack_c(u+1,j,cnt1,cnt2,structure, vars);
              return;
            }
        if(E_C_rem[uj] != INF){
          if(e == (E_M_rem[iu] + E_C_rem[uj] + energy)){
            backtrack_m(i,u,-1,-1,structure,vars);
            backtrack_c(u+1,j,-1,-1,structure,vars);
            return;
          }
        }
      }
      if(E_C_rem[uj] != INF){
        if(E_M[iu])
        for(cnt1 = k_min_values_m[iu];
            cnt1 <= k_max_values_m[iu];
            cnt1++)
          for(cnt2 = l_min_values_m[iu][cnt1];
              cnt2 <= l_max_values_m[iu][cnt1];
              cnt2 += 2)
            if(e == (E_M[iu][cnt1][cnt2/2] + E_C_rem[uj] + energy)){
              backtrack_m(i,u,cnt1,cnt2,structure,vars);
              backtrack_c(u+1,j,-1,-1,structure,vars);
              return;
            }
      }

      if(!E_M[iu]) continue;
      if(!E_C[uj]) continue;

      for(cnt1 = k_min_values_m[iu];
          cnt1 <= k_max_values_m[iu];
          cnt1++)
        for(cnt2 = l_min_values_m[iu][cnt1];
            cnt2 <= l_max_values_m[iu][cnt1];
            cnt2 += 2)
          for(cnt3 = k_min_values[uj];
              cnt3 <= k_max_values[uj];
              cnt3++){
            for(cnt4 = l_min_values[uj][cnt3];
                cnt4 <= l_max_values[uj][cnt3];
                cnt4 += 2)
              if(((cnt1 + cnt3 + d1) > maxD1) || ((cnt2 + cnt4 + d2) > maxD2))
                if(e == (E_M[iu][cnt1][cnt2/2] + E_C[uj][cnt3][cnt4/2] + energy)){
                  backtrack_m(i, u, cnt1, cnt2, structure, vars);
                  backtrack_c(u+1, j, cnt3, cnt4, structure, vars);
                  return;
                }
          }
    }

  } /* end if (k == -1) */
  else{
    d1 = base_d1 - referenceBPs1[my_iindx[i+1]-j];
    d2 = base_d2 - referenceBPs2[my_iindx[i+1]-j];
    /* new_fML = ML(i+1,j)+c */
    if(d1 <= k && d2 <= l)
      if((k-d1 >= k_min_values_m[my_iindx[i+1]-j]) && (k-d1 <= k_max_values_m[my_iindx[i+1]-j]))
        if((l-d2 >= l_min_values_m[my_iindx[i+1]-j][k-d1]) && (l-d2 <= l_max_values_m[my_iindx[i+1]-j][k-d1])){
          if(E_M[my_iindx[i+1]-j][k-d1][(l-d2)/2] + P->MLbase == e){
            backtrack_m(i+1, j, k-d1, l-d2, structure, vars);
            return;
          }
        }

    d1 = base_d1 - referenceBPs1[ij+1];
    d2 = base_d2 - referenceBPs2[ij+1];

    /* new_fML = min(ML(i,j-1) + c, new_fML) */
    if(E_M[ij+1])
      if(d1 <= k && d2 <= l)
        if((k-d1 >= k_min_values_m[ij+1]) && (k-d1 <= k_max_values_m[ij+1]))
          if((l-d2 >= l_min_values_m[ij+1][k-d1]) && (l-d2 <= l_max_values_m[ij+1][k-d1]))
            if(E_M[ij+1][k-d1][(l-d2)/2] + P->MLbase == e){
              backtrack_m(i, j-1, k-d1, l-d2, structure, vars);
              return;
            }

    /* new_fML = min(new_fML, C(i,j)+b) */
    if(E_C[ij]){
      type = ptype[ij];

      if(dangles == 2)
        energy = E_MLstem(type, ((i > 1) || circ) ? S1[i-1] : -1, ((j < seq_length) || circ) ? S1[j+1] : -1, P);
      else
        energy = E_MLstem(type, -1, -1, P);

      if((k >= k_min_values[ij]) && (k <= k_max_values[ij]))
        if((l >= l_min_values[ij][k]) && (l <= l_max_values[ij][k])){
          if(E_C[ij][k][l/2] + energy == e){
            backtrack_c(i, j, k, l, structure, vars);
            return;
          }
        }
    }

      /* modular decomposition -------------------------------*/

    for(u = i+1+TURN; u <= j-2-TURN; u++){
      if(!E_M[my_iindx[i]-u]) continue;
      if(!E_C[my_iindx[u+1]-j]) continue;
      type = ptype[my_iindx[u+1]-j];

      d1 = base_d1 - referenceBPs1[my_iindx[i]-u] - referenceBPs1[my_iindx[u+1]-j];
      d2 = base_d2 - referenceBPs2[my_iindx[i]-u] - referenceBPs2[my_iindx[u+1]-j];

      if(dangles == 2)
        energy = E_MLstem(type, S1[u], ((j < seq_length) || circ) ? S1[j+1] : -1, P);
      else
        energy = E_MLstem(type, -1, -1, P);

      if(d1 <= k && d2 <= l)
        for(cnt1 = k_min_values_m[my_iindx[i]-u]; cnt1 <= MIN2(k-d1, k_max_values_m[my_iindx[i]-u]); cnt1++)
          for(cnt2 = l_min_values_m[my_iindx[i]-u][cnt1]; cnt2 <= MIN2(l-d2, l_max_values_m[my_iindx[i]-u][cnt1]); cnt2+=2)
            if((k-d1-cnt1 >= k_min_values[my_iindx[u+1]-j]) && (k-d1-cnt1 <= k_max_values[my_iindx[u+1]-j]))
              if((l-d2-cnt2 >= l_min_values[my_iindx[u+1]-j][k-d1-cnt1]) && (l-d2-cnt2 <= l_max_values[my_iindx[u+1]-j][k-d1-cnt1]))
                if(E_M[my_iindx[i]-u][cnt1][cnt2/2] + E_C[my_iindx[u+1]-j][k-d1-cnt1][(l-d2-cnt2)/2] + energy == e){
                  backtrack_m(i, u, cnt1, cnt2, structure, vars);
                  backtrack_c(u+1, j, k-d1-cnt1, l-d2-cnt2, structure, vars);
                  return;
                }
    }
  }
  nrerror("backtracking failed in fML\n");
}

PRIVATE void backtrack_m1(unsigned int i, unsigned int j, int k, int l, char *structure, TwoDfold_vars *vars){
  unsigned int  ij, seq_length, d1, d2, *referenceBPs1, *referenceBPs2, maxD1, maxD2;
  int           *my_iindx, **l_min_values, **l_max_values,**l_min_values_m1, **l_max_values_m1;
  int           *k_min_values, *k_max_values,*k_min_values_m1, *k_max_values_m1, cnt1, cnt2;
  int           ***E_C, ***E_M1, *E_C_rem, *E_M1_rem, type, dangles, circ, energy, e_m1;

  short         *S1;
  char          *ptype;
  paramT        *P;

  P               = vars->P;
  seq_length      = vars->seq_length;
  S1              = vars->S1;
  ptype           = vars->ptype;
  circ            = vars->circ;
  my_iindx        = vars->my_iindx;
  referenceBPs1   = vars->referenceBPs1;
  referenceBPs2   = vars->referenceBPs2;
  dangles         = vars->dangles;

  E_C             = vars->E_C;
  l_min_values    = vars->l_min_values;
  l_max_values    = vars->l_max_values;
  k_min_values    = vars->k_min_values;
  k_max_values    = vars->k_max_values;

  E_M1            = vars->E_M1;
  l_min_values_m1 = vars->l_min_values_m1;
  l_max_values_m1 = vars->l_max_values_m1;
  k_min_values_m1 = vars->k_min_values_m1;
  k_max_values_m1 = vars->k_max_values_m1;

  E_C_rem        = vars->E_C_rem;
  E_M1_rem       = vars->E_M1_rem;
  maxD1           = vars->maxD1;
  maxD2           = vars->maxD2;

  ij    = my_iindx[i]-j;
  e_m1  = (k == -1) ? E_M1_rem[ij] : E_M1[ij][k][l/2];

  type = ptype[ij];
  d1 = referenceBPs1[ij] - referenceBPs1[ij+1];
  d2 = referenceBPs2[ij] - referenceBPs2[ij+1];

  if(dangles == 2)
    energy = E_MLstem(type, (i > 1) || circ ? S1[i-1] : -1, (j < seq_length) || circ ? S1[j+1] : -1, P);
  else
    energy = E_MLstem(type, -1, -1, P);

  if(k == -1){
    if(E_C_rem[ij] != INF){
      if(e_m1 == (E_C_rem[ij] + energy)){
        backtrack_c(i,j,-1,-1,structure,vars);
        return;
      }
    }
    if(E_M1_rem[ij+1] != INF){
      if(e_m1 == (E_M1_rem[ij+1] + P->MLbase)){
        backtrack_m1(i,j-1,-1,-1,structure,vars);
        return;
      }
    }
    for(cnt1 = k_min_values_m1[ij+1];
        cnt1 <= k_max_values_m1[ij+1];
        cnt1++)
      for(cnt2 = l_min_values_m1[ij+1][cnt1];
          cnt2 <= l_max_values_m1[ij+1][cnt1];
          cnt2 += 2)
        if(((cnt1 + d1) > maxD1) || ((cnt2 + d2) > maxD2)){
          if(e_m1  == (E_M1[ij+1][cnt1][cnt2/2] + P->MLbase)){
            backtrack_m1(i,j-1,cnt1,cnt2,structure,vars);
            return;
          }
        }
  }
  else{
    if(E_C[ij])
      if((k >= k_min_values[ij]) && (k <= k_max_values[ij]))
        if((l >= l_min_values[ij][k]) && (l <= l_max_values[ij][k]))
          if(E_C[ij][k][l/2] + energy == e_m1){
            backtrack_c(i, j, k, l, structure, vars);
            return;
          }

    if(d1 <= k && d2 <= l)
      if((k-d1 >= k_min_values_m1[ij+1]) && (k-d1 <= k_max_values_m1[ij+1]))
        if((l-d2 >= l_min_values_m1[ij+1][k-d1]) && (l-d2 <= l_max_values_m1[ij+1][k-d1]))
          if(E_M1[ij+1][k-d1][(l-d2)/2] + P->MLbase == e_m1){
            backtrack_m1(i, j-1, k-d1, l-d2, structure, vars);
            return;
          }
  }
  nrerror("backtack failed in m1\n");
}

PRIVATE void backtrack_fc(int k, int l, char *structure, TwoDfold_vars *vars){
  unsigned int   d, i, j, seq_length, base_d1, base_d2, d1, d2, maxD1, maxD2;
  int   *my_iindx, energy, cnt1, cnt2, cnt3, cnt4;
  short *S1;
  unsigned int   *referenceBPs1, *referenceBPs2;
  char  *sequence, *ptype;
  int **E_Fc, **E_FcH, **E_FcI, **E_FcM, ***E_C, ***E_M, ***E_M2;
  int *E_C_rem, *E_M1_rem, *E_M_rem, *E_M2_rem, E_Fc_rem, E_FcH_rem, E_FcI_rem, E_FcM_rem;
  int **l_min_values, **l_max_values, *k_min_values, *k_max_values;
  int **l_min_values_m, **l_max_values_m, *k_min_values_m, *k_max_values_m;
  int **l_min_values_m2, **l_max_values_m2, *k_min_values_m2, *k_max_values_m2;
  int *l_min_values_fc, *l_max_values_fc, k_min_values_fc, k_max_values_fc;
  int *l_min_values_fcH, *l_max_values_fcH, k_min_values_fcH, k_max_values_fcH;
  int *l_min_values_fcI, *l_max_values_fcI, k_min_values_fcI, k_max_values_fcI;
  int *l_min_values_fcM, *l_max_values_fcM, k_min_values_fcM, k_max_values_fcM;
  paramT   *P;
  P                 = vars->P;
  sequence          = vars->sequence;
  seq_length        = vars->seq_length;
  S1                = vars->S1;
  ptype             = vars->ptype;
  my_iindx          = vars->my_iindx;
  referenceBPs1     = vars->referenceBPs1;
  referenceBPs2     = vars->referenceBPs2;

  base_d1           = referenceBPs1[my_iindx[1]-seq_length];
  base_d2           = referenceBPs2[my_iindx[1]-seq_length];

  E_C               = vars->E_C;
  l_min_values      = vars->l_min_values;
  l_max_values      = vars->l_max_values;
  k_min_values      = vars->k_min_values;
  k_max_values      = vars->k_max_values;

  E_M               = vars->E_M;
  l_min_values_m    = vars->l_min_values_m;
  l_max_values_m    = vars->l_max_values_m;
  k_min_values_m    = vars->k_min_values_m;
  k_max_values_m    = vars->k_max_values_m;

  E_M2              = vars->E_M2;
  l_min_values_m2   = vars->l_min_values_m2;
  l_max_values_m2   = vars->l_max_values_m2;
  k_min_values_m2   = vars->k_min_values_m2;
  k_max_values_m2   = vars->k_max_values_m2;

  E_Fc              = vars->E_Fc;
  l_min_values_fc   = vars->l_min_values_fc;
  l_max_values_fc   = vars->l_max_values_fc;
  k_min_values_fc   = vars->k_min_values_fc;
  k_max_values_fc   = vars->k_max_values_fc;

  E_FcI             = vars->E_FcI;
  l_min_values_fcI  = vars->l_min_values_fcI;
  l_max_values_fcI  = vars->l_max_values_fcI;
  k_min_values_fcI  = vars->k_min_values_fcI;
  k_max_values_fcI  = vars->k_max_values_fcI;

  E_FcH             = vars->E_FcH;
  l_min_values_fcH  = vars->l_min_values_fcH;
  l_max_values_fcH  = vars->l_max_values_fcH;
  k_min_values_fcH  = vars->k_min_values_fcH;
  k_max_values_fcH  = vars->k_max_values_fcH;

  E_FcM             = vars->E_FcM;
  l_min_values_fcM  = vars->l_min_values_fcM;
  l_max_values_fcM  = vars->l_max_values_fcM;
  k_min_values_fcM  = vars->k_min_values_fcM;
  k_max_values_fcM  = vars->k_max_values_fcM;

  E_C_rem          = vars->E_C_rem;
  E_M1_rem         = vars->E_M1_rem;
  E_M_rem          = vars->E_M_rem;
  E_M2_rem         = vars->E_M2_rem;
  E_Fc_rem         = vars->E_Fc_rem;
  E_FcH_rem        = vars->E_FcH_rem;
  E_FcI_rem        = vars->E_FcI_rem;
  E_FcM_rem        = vars->E_FcM_rem;
  maxD1             = vars->maxD1;
  maxD2             = vars->maxD2;


  if(k==-1){
    /* check if mfe might be open chain */
    if(E_Fc_rem == 0)
      if((referenceBPs1[my_iindx[1]-seq_length] > maxD1) || (referenceBPs2[my_iindx[1]-seq_length] > maxD2))
        return;

    /* check for hairpin configurations */
    if(E_Fc_rem == E_FcH_rem){
      for (d = TURN+2; d <= seq_length; d++) /* i,j in [1..length] */
        for (j = d; j <= seq_length; j++) {
          unsigned int u, ij;
          int type, no_close;
          char loopseq[10];
          i = j-d+1;
          ij = my_iindx[i]-j;
          u = seq_length-j + i-1;
          if (u<TURN) continue;
          type = ptype[ij];
          no_close = (((type==3)||(type==4))&&no_closingGU);
          type=rtype[type];
          if (!type) continue;
          if(no_close) continue;

          d1 = base_d1 - referenceBPs1[ij];
          d2 = base_d2 - referenceBPs2[ij];
          if (u<7) {
            strcpy(loopseq , sequence+j-1);
            strncat(loopseq, sequence, i);
          }
          energy = E_Hairpin(u, type, S1[j+1], S1[i-1],  loopseq, P);

          if(E_C_rem[ij] != INF){
            if(E_Fc_rem == (E_C_rem[ij] + energy)){
              backtrack_c(i,j,-1,-1,structure,vars);
              return;
            }
          }
          if(E_C[ij])
            for(cnt1 = k_min_values[ij];
                cnt1 <= k_max_values[ij];
                cnt1++)
              for(cnt2 = l_min_values[ij][cnt1];
                  cnt2 <= l_max_values[ij][cnt1];
                  cnt2 += 2)
                if(((cnt1 + d1) > maxD1) || ((cnt2 + d2) > maxD2))
                  if(E_Fc_rem == (E_C[ij][cnt1][cnt2/2] + energy)){
                    backtrack_c(i,j,cnt1,cnt2,structure,vars);
                    return;
                  }
        }
    }

    /* check for interior loop configurations */
    if(E_Fc_rem == E_FcI_rem){
      for (d = TURN+2; d <= seq_length; d++) /* i,j in [1..length] */
        for (j = d; j <= seq_length; j++) {
          unsigned int u, ij, p, q, pq;
          int type, type_2;
          i = j-d+1;
          ij = my_iindx[i]-j;
          u = seq_length-j + i-1;
          if (u<TURN) continue;
          type = rtype[ptype[ij]];
          if (!type) continue;

          for(p = j+1; p < seq_length ; p++){
            unsigned int u1, qmin, ln_pre;
            u1 = p-j-1;
            if (u1+i-1>MAXLOOP) break;
            qmin = p + TURN + 1;
            ln_pre = u1 + i + seq_length;
            if(ln_pre > qmin + MAXLOOP) qmin = ln_pre - MAXLOOP - 1;
            for(q = qmin; q <= seq_length; q++){
              unsigned int u2;
              pq = my_iindx[p]-q;
              type_2 = rtype[ptype[pq]];
              if (type_2==0) continue;
              u2 = i-1 + seq_length-q;
              if (u1+u2>MAXLOOP) continue;
              energy = E_IntLoop(u1, u2, type, type_2, S1[j+1], S1[i-1], S1[p-1], S1[q+1], P);
              if(E_C_rem[ij] != INF){
                if(E_C[pq])
                  for(cnt1 = k_min_values[pq];
                      cnt1 <= k_max_values[pq];
                      cnt1++)
                    for(cnt2 = l_min_values[pq][cnt1];
                        cnt2 <= l_max_values[pq][cnt1];
                        cnt2 += 2)
                      if(E_Fc_rem == (E_C_rem[ij] + E_C[pq][cnt1][cnt2/2] + energy)){
                        backtrack_c(i,j,-1,-1,structure,vars);
                        backtrack_c(p,q,cnt1,cnt2,structure,vars);
                        return;
                      }
                if(E_C_rem[pq] != INF){
                  if(E_Fc_rem == (E_C_rem[ij] + E_C_rem[pq] + energy)){
                    backtrack_c(i,j,-1,-1,structure,vars);
                    backtrack_c(p,q,-1,-1,structure,vars);
                    return;
                  }
                }
              }
              if(E_C_rem[pq] != INF){
                if(E_C[ij])
                  for(cnt1 = k_min_values[ij];
                      cnt1 <= k_max_values[ij];
                      cnt1++)
                    for(cnt2 = l_min_values[ij][cnt1];
                        cnt2 <= l_max_values[ij][cnt1];
                        cnt2 += 2)
                      if(E_Fc_rem == (E_C[ij][cnt1][cnt2/2] + E_C_rem[pq] + energy)){
                        backtrack_c(i,j,cnt1,cnt2,structure,vars);
                        backtrack_c(p,q,-1,-1,structure,vars);
                        return;
                      }
              }

              if(!(E_C[ij])) continue;
              if(!(E_C[pq])) continue;

              /* get distance to reference if closing the interior loop
              *  d2a = dbp(T1_[1,n}, T1_{p,q} + T1_{i,j})
              *  d2b = dbp(T2_[1,n}, T2_{p,q} + T2_{i,j})
              */
              d1 = base_d1 - referenceBPs1[ij] - referenceBPs1[pq];
              d2 = base_d2 - referenceBPs2[ij] - referenceBPs2[pq];
              for(cnt1 = k_min_values[ij];
                  cnt1 <= k_max_values[ij];
                  cnt1++)
                for(cnt2 = l_min_values[ij][cnt1];
                    cnt2 <= l_max_values[ij][cnt1];
                    cnt2 += 2)
                  for(cnt3 = k_min_values[pq];
                      cnt3 <= k_max_values[pq];
                      cnt3++)
                    for(cnt4 = l_min_values[pq][cnt3];
                        cnt4 <= l_max_values[pq][cnt3];
                        cnt4 += 2)
                      if(((cnt1 + cnt3 + d1) > maxD1) || ((cnt2 + cnt4 + d2) > maxD2))
                        if(E_Fc_rem == (E_C[ij][cnt1][cnt2/2] + E_C[pq][cnt3][cnt4/2] + energy)){
                          backtrack_c(i, j, cnt1, cnt2, structure, vars);
                          backtrack_c(p, q, cnt3, cnt4, structure, vars);
                          return;
                        }
            } /* end for p */
          } /* end for q */
        }

    }

    /* check for multi loop configurations */
    if(E_Fc_rem == E_FcM_rem){
      if(seq_length > 2*TURN)
        for (i=TURN+1; i<seq_length-2*TURN; i++) {
          /* get distancies to references
          * d3a = dbp(T1_[1,n}, T1_{1,k} + T1_{k+1, n})
          * d3b = dbp(T2_[1,n}, T2_{1,k} + T2_{k+1, n})
          */
          if(E_M_rem[my_iindx[1]-i] != INF){
            if(E_M2[i+1])
              for(cnt1 = k_min_values_m2[i+1];
                  cnt1 <= k_max_values_m2[i+1];
                  cnt1++)
                for(cnt2 = l_min_values_m2[i+1][cnt1];
                    cnt2 <= l_max_values_m2[i+1][cnt1];
                    cnt2 += 2)
                  if(E_Fc_rem == (E_M_rem[my_iindx[1]-i] + E_M2[i+1][cnt1][cnt2/2] + P->MLclosing)){
                    backtrack_m(1,i,-1,-1,structure,vars);
                    backtrack_m2(i+1,cnt1,cnt2,structure,vars);
                    return;
                  }
            if(E_M2_rem[i+1] != INF){
              if(E_Fc_rem == (E_M_rem[my_iindx[1]-i] + E_M2_rem[i+1] + P->MLclosing)){
                backtrack_m(1,i,-1,-1,structure,vars);
                backtrack_m2(i+1,-1,-1,structure,vars);
                return;
              }
            }
          }
          if(E_M2_rem[i+1] != INF){
            if(E_M[my_iindx[1]-i])
              for(cnt1 = k_min_values_m[my_iindx[1]-i];
                  cnt1 <= k_max_values_m[my_iindx[1]-i];
                  cnt1++)
                for(cnt2 = l_min_values_m[my_iindx[1]-i][cnt1];
                    cnt2 <= l_max_values_m[my_iindx[1]-i][cnt1];
                    cnt2 += 2)
                  if(E_Fc_rem == (E_M[my_iindx[1]-i][cnt1][cnt2/2] + E_M2_rem[i+1] + P->MLclosing)){
                    backtrack_m(1,i,cnt1,cnt2,structure,vars);
                    backtrack_m2(i+1,-1,-1,structure,vars);
                    return;
                  }
          }

          if(!(E_M[my_iindx[1]-i])) continue;
          if(!(E_M2[i+1])) continue;

          d1 = base_d1 - referenceBPs1[my_iindx[1]-i] - referenceBPs1[my_iindx[i+1]-seq_length];
          d2 = base_d2 - referenceBPs2[my_iindx[1]-i] - referenceBPs2[my_iindx[i+1]-seq_length];
          for(cnt1 = k_min_values_m[my_iindx[1]-i];
              cnt1 <= k_max_values_m[my_iindx[1]-i];
              cnt1++)
            for(cnt2 = l_min_values_m[my_iindx[1]-i][cnt1];
                cnt2 <= l_max_values_m[my_iindx[1]-i][cnt1];
                cnt2 += 2)
              for(cnt3 = k_min_values_m2[i+1];
                  cnt3 <= k_max_values_m2[i+1];
                  cnt3++)
                for(cnt4 = l_min_values_m2[i+1][cnt3];
                    cnt4 <= l_max_values_m2[i+1][cnt3];
                    cnt4 += 2)
                  if(((cnt1 + cnt3 + d1) > maxD1) || ((cnt2 + cnt4 + d2) > maxD2)){
                    if(E_Fc_rem == (E_M[my_iindx[1]-i][cnt1][cnt2/2] + E_M2[i+1][cnt3][cnt4/2] + P->MLclosing)){
                      backtrack_m(1, i, cnt1, cnt2, structure, vars);
                      backtrack_m2(i+1, cnt3, cnt4, structure, vars);
                      return;
                    }
                  }
        }
    }
  }
  else{
    /* open chain ? */
    if(E_Fc[k][l/2] == 0)
      if((k == referenceBPs1[my_iindx[1]-seq_length]) && (l == referenceBPs2[my_iindx[1]-seq_length])){
        return;
      }
    if((k >= k_min_values_fcH) && (k <= k_max_values_fcH)){
      if((l >= l_min_values_fcH[k]) && (l <= l_max_values_fcH[k]))
        if(E_Fc[k][l/2] == E_FcH[k][l/2]){
          for (d = TURN+2; d <= seq_length; d++) /* i,j in [1..length] */
            for (j = d; j <= seq_length; j++) {
              unsigned int u, ij;
              int type, no_close;
              char loopseq[10];
              i = j-d+1;
              ij = my_iindx[i]-j;
              if (!E_C[ij]) continue;
              u = seq_length-j + i-1;
              if (u<TURN) continue;

              type = ptype[ij];

              no_close = (((type==3)||(type==4))&&no_closingGU);

              type=rtype[type];

              if (!type) continue;
              if(no_close) continue;

              d1 = base_d1 - referenceBPs1[ij];
              d2 = base_d2 - referenceBPs2[ij];
              if (u<7) {
                strcpy(loopseq , sequence+j-1);
                strncat(loopseq, sequence, i);
              }
              energy = E_Hairpin(u, type, S1[j+1], S1[i-1],  loopseq, P);
              if((k >= d1) && (l >= d2))
                if((k-d1 >= k_min_values[ij]) && (k-d1 <= k_max_values[ij]))
                  if((l-d2 >= l_min_values[ij][k-d1]) && (l-d2 <= l_max_values[ij][k-d1])){
                    if(E_Fc[k][l/2] == E_C[ij][k-d1][(l-d2)/2] + energy){
                      backtrack_c(i, j, k-d1, l-d2, structure, vars);
                      return;
                    }
                  }
            }
        }
    }

    if((k >= k_min_values_fcI) && (k <= k_max_values_fcI)){
      if((l >= l_min_values_fcI[k]) && (l <= l_max_values_fcI[k]))
        if(E_Fc[k][l/2] == E_FcI[k][l/2]){
          for (d = TURN+2; d <= seq_length; d++) /* i,j in [1..length] */
            for (j = d; j <= seq_length; j++) {
              unsigned int u, ij, p, q, pq;
              int type, type_2;
              i = j-d+1;
              ij = my_iindx[i]-j;
              if(!E_C[ij]) continue;
              u = seq_length-j + i-1;
              if (u<TURN) continue;

              type = ptype[ij];

              type=rtype[type];

              if (!type) continue;

              for(p = j+1; p < seq_length ; p++){
                unsigned int u1, qmin, ln_pre;
                u1 = p-j-1;
                if (u1+i-1>MAXLOOP) break;
                qmin = p + TURN + 1;
                ln_pre = u1 + i + seq_length;
                if(ln_pre > qmin + MAXLOOP) qmin = ln_pre - MAXLOOP - 1;
                for(q = qmin; q <= seq_length; q++){
                  unsigned int u2;
                  pq = my_iindx[p]-q;
                  if(!E_C[pq]) continue;
                  type_2 = rtype[ptype[pq]];
                  if (type_2==0) continue;
                  u2 = i-1 + seq_length-q;
                  if (u1+u2>MAXLOOP) continue;
                  /* get distance to reference if closing the interior loop
                  *  d2a = dbp(T1_[1,n}, T1_{p,q} + T1_{i,j})
                  *  d2b = dbp(T2_[1,n}, T2_{p,q} + T2_{i,j})
                  */
                  d1 = base_d1 - referenceBPs1[ij] - referenceBPs1[pq];
                  d2 = base_d2 - referenceBPs2[ij] - referenceBPs2[pq];
                  energy = E_IntLoop(u1, u2, type, type_2, S1[j+1], S1[i-1], S1[p-1], S1[q+1], P);
                  if((k >= d1) && (l >= d2))
                    for(cnt1 = k_min_values[ij]; cnt1 <= MIN2(k_max_values[ij], k - d1); cnt1++)
                      for(cnt2 = l_min_values[ij][cnt1]; cnt2 <= MIN2(l_max_values[ij][cnt1], l - d2); cnt2+=2)
                        if((k - d1 - cnt1 >= k_min_values[pq]) && (k - d1 - cnt1 <= k_max_values[pq]))
                          if((l - d2 - cnt2 >= l_min_values[pq][k-d1-cnt1]) && (l - d2 - cnt2 <= l_max_values[pq][k-d1-cnt1])){
                            if((E_C[ij][cnt1][cnt2/2] + E_C[pq][k-d1-cnt1][(l-d2-cnt2)/2] + energy) == E_Fc[k][l/2]){
                              backtrack_c(i, j, cnt1, cnt2, structure, vars);
                              backtrack_c(p, q, k - d1 - cnt1, l - d2 - cnt2, structure, vars);
                              return;
                            }
                          }
                }
              }
            }
        }
    }

    if((k >= k_min_values_fcM) && (k <= k_max_values_fcM)){
      if((l >= l_min_values_fcM[k]) && (l <= l_max_values_fcM[k]))
        if(E_Fc[k][l/2] == E_FcM[k][l/2]){
          if(seq_length > 2*TURN)
            for (i=TURN+1; i<seq_length-2*TURN; i++) {
              /* get distancies to references
              * d3a = dbp(T1_[1,n}, T1_{1,k} + T1_{k+1, n})
              * d3b = dbp(T2_[1,n}, T2_{1,k} + T2_{k+1, n})
              */
              if(!E_M[my_iindx[1]-i]) continue;
              if(!E_M2[i+1]) continue;
              d1 = base_d1 - referenceBPs1[my_iindx[1]-i] - referenceBPs1[my_iindx[i+1]-seq_length];
              d2 = base_d2 - referenceBPs2[my_iindx[1]-i] - referenceBPs2[my_iindx[i+1]-seq_length];
              if((k >= d1) && (l >= d2))
                for(cnt1 = k_min_values_m[my_iindx[1]-i]; cnt1 <= MIN2(k_max_values_m[my_iindx[1]-i], k-d1); cnt1++)
                  for(cnt2 = l_min_values_m[my_iindx[1]-i][cnt1]; cnt2 <= MIN2(l_max_values_m[my_iindx[1]-i][cnt1], l-d2); cnt2+=2)
                    if((k - d1 - cnt1 >= k_min_values_m2[i+1]) && (k - d1 - cnt1 <= k_max_values_m2[i+1]))
                      if((l - d2 - cnt2 >= l_min_values_m2[i+1][k-d1-cnt1]) && (l - d2 - cnt2 <= l_max_values_m2[i+1][k-d1-cnt1]))
                        if((E_M[my_iindx[1]-i][cnt1][cnt2/2] + E_M2[i+1][k-d1-cnt1][(l-d2-cnt2)/2] + P->MLclosing) == E_FcM[k][l/2]){
                          backtrack_m(1, i, cnt1, cnt2, structure, vars);
                          backtrack_m2(i+1, k - d1 - cnt1, l - d2 - cnt2, structure, vars);
                          return;
                        }
            }
        }
    }
  }
  nrerror("backtack failed in fc\n");
}


PRIVATE void backtrack_m2(unsigned int i, int k, int l, char *structure, TwoDfold_vars *vars){
  unsigned int   j, ij, j3, n;
  unsigned int   *referenceBPs1, *referenceBPs2;
  unsigned int d1, d2, base_d1, base_d2, maxD1, maxD2;
  int *my_iindx, cnt1, cnt2, cnt3, cnt4;
  int ***E_M1, ***E_M2, *E_M2_rem, *E_M1_rem, e;
  int **l_min_values_m1, **l_max_values_m1, *k_min_values_m1, *k_max_values_m1;
  int **l_min_values_m2, **l_max_values_m2, *k_min_values_m2, *k_max_values_m2;

  n               = vars->seq_length;
  my_iindx        = vars->my_iindx;
  referenceBPs1   = vars->referenceBPs1;
  referenceBPs2   = vars->referenceBPs2;

  E_M1            = vars->E_M1;
  l_min_values_m1 = vars->l_min_values_m1;
  l_max_values_m1 = vars->l_max_values_m1;
  k_min_values_m1 = vars->k_min_values_m1;
  k_max_values_m1 = vars->k_max_values_m1;

  E_M1_rem        = vars->E_M1_rem;

  E_M2            = vars->E_M2;
  l_min_values_m2 = vars->l_min_values_m2;
  l_max_values_m2 = vars->l_max_values_m2;
  k_min_values_m2 = vars->k_min_values_m2;
  k_max_values_m2 = vars->k_max_values_m2;

  E_M2_rem        = vars->E_M2_rem;

  maxD1           = vars->maxD1;
  maxD2           = vars->maxD2;

  base_d1 = referenceBPs1[my_iindx[i]-n];
  base_d2 = referenceBPs2[my_iindx[i]-n];

  if(k == -1){
    e = E_M2_rem[i];
    for (j=i+TURN+1; j<n-TURN-1; j++){
      if(E_M1_rem[my_iindx[i]-j] != INF){
        if(E_M1[my_iindx[j+1]-n])
          for(cnt1 = k_min_values_m1[my_iindx[j+1]-n];
              cnt1 <= k_max_values_m1[my_iindx[j+1]-n];
              cnt1++)
            for(cnt2 = l_min_values_m1[my_iindx[j+1]-n][cnt1];
                cnt2 <= l_max_values_m1[my_iindx[j+1]-n][cnt1];
                cnt2++)
              if(e == E_M1_rem[my_iindx[i]-j] + E_M1[my_iindx[j+1]-n][cnt1][cnt2/2]){
                backtrack_m1(i, j, k, l, structure, vars);
                backtrack_m1(j+1, n, cnt1, cnt2, structure, vars);
                return;
              }
        if(E_M1_rem[my_iindx[j+1]-n] != INF){
          if(e == E_M1_rem[my_iindx[i]-j] + E_M1_rem[my_iindx[j+1]-n]){
            backtrack_m1(i, j, k, l, structure, vars);
            backtrack_m1(j+1, n, k, l, structure, vars);
            return;
          }
        }
      }
      if(E_M1_rem[my_iindx[j+1]-n] != INF){
        if(E_M1[my_iindx[i]-j])
          for(cnt1 = k_min_values_m1[my_iindx[i]-j];
              cnt1 <= k_max_values_m1[my_iindx[i]-j];
              cnt1++)
            for(cnt2 = l_min_values_m1[my_iindx[i]-j][cnt1];
                cnt2 <= l_max_values_m1[my_iindx[i]-j][cnt1];
                cnt2 += 2)
              if(e == E_M1[my_iindx[i]-j][cnt1][cnt2/2] + E_M1_rem[my_iindx[j+1]-n]){
                backtrack_m1(i, j, cnt1, cnt2, structure, vars);
                backtrack_m1(j+1, n, k, l, structure, vars);
                return;
              }
      }


      if(!E_M1[my_iindx[i]-j]) continue;
      if(!E_M1[my_iindx[j+1]-n]) continue;

      d1 = referenceBPs1[my_iindx[i]-n] - referenceBPs1[my_iindx[i]-j] - referenceBPs1[my_iindx[j+1]-n];
      d2 = referenceBPs2[my_iindx[i]-n] - referenceBPs2[my_iindx[i]-j] - referenceBPs2[my_iindx[j+1]-n];

      for(cnt1 = k_min_values_m1[my_iindx[i]-j]; cnt1 <= k_max_values_m1[my_iindx[i]-j]; cnt1++)
        for(cnt2 = l_min_values_m1[my_iindx[i]-j][cnt1]; cnt2 <= l_max_values_m1[my_iindx[i]-j][cnt1]; cnt2+=2){
          for(cnt3 = k_min_values_m1[my_iindx[j+1]-n]; cnt3 <= k_max_values_m1[my_iindx[j+1]-n]; cnt3++)
            for(cnt4 = l_min_values_m1[my_iindx[j+1]-n][cnt3]; cnt4 <= l_max_values_m1[my_iindx[j+1]-n][cnt3]; cnt4+=2){
              if(((cnt1 + cnt3 + d1) > maxD1) || ((cnt2 + cnt4 + d2) > maxD2)){
                if(e == E_M1[my_iindx[i]-j][cnt1][cnt2/2] + E_M1[my_iindx[j+1]-n][cnt3][cnt4/2]){
                  backtrack_m1(i, j, cnt1, cnt2, structure, vars);
                  backtrack_m1(j+1, n, cnt3, cnt4, structure, vars);
                  return;
                }
              }
            }
        }
    }
  }
  else{
    for(j=i+TURN+1; j<n-TURN-1; j++){
      if(!E_M1[my_iindx[i]-j]) continue;
      if(!E_M1[my_iindx[j+1]-n]) continue;

      ij = my_iindx[i]-j;
      j3 = my_iindx[j+1]-n;
      d1 = base_d1 - referenceBPs1[ij] - referenceBPs1[j3];
      d2 = base_d2 - referenceBPs2[ij] - referenceBPs2[j3];

      for(cnt1 = k_min_values_m1[ij]; cnt1 <= MIN2(k_max_values_m1[ij], k - d1); cnt1++)
        for(cnt2 = l_min_values_m1[ij][cnt1]; cnt2 <= MIN2(l_max_values_m1[ij][cnt1], l-d2); cnt2+=2)
          if((k - d1 - cnt1 >= k_min_values_m1[j3]) && (k - d1 - cnt1 <= k_max_values_m1[j3]))
            if((l - d2 - cnt2 >= l_min_values_m1[j3][k - d1 - cnt1]) && (l - d2 - cnt2 <= l_max_values_m1[j3][k-d1-cnt1]))
              if(E_M1[ij][cnt1][cnt2/2] + E_M1[j3][k-d1-cnt1][(l-d2-cnt2)/2] == E_M2[i][k][l/2]){
                backtrack_m1(i, j, cnt1, cnt2, structure, vars);
                backtrack_m1(j+1, n, k-d1-cnt1, l-d2-cnt2, structure, vars);
                return;
              }
    }
  }
  nrerror("backtack failed in m2\n");
}

PRIVATE void mfe_circ(TwoDfold_vars *vars){
  unsigned int  d, i, j, k, maxD1, maxD2, seq_length, *referenceBPs1, *referenceBPs2, d1, d2, base_d1, base_d2, *mm1, *mm2, *bpdist;
  int           *my_iindx, energy, dangles, cnt1, cnt2, cnt3, cnt4;
  short         *S1, *reference_pt1, *reference_pt2;
  char          *sequence, *ptype;
  int           ***E_C, ***E_M, ***E_M1, ***E_M2, **E_Fc, **E_FcH, **E_FcI, **E_FcM;
  int           *E_C_rem, *E_M_rem, *E_M1_rem, *E_M2_rem, E_Fc_rem, E_FcH_rem, E_FcI_rem, E_FcM_rem;
  int           **l_min_values, **l_max_values, **l_min_values_m, **l_max_values_m, **l_min_values_m1, **l_max_values_m1, **l_min_values_m2, **l_max_values_m2;
  int           *k_min_values, *k_max_values,*k_min_values_m, *k_max_values_m,*k_min_values_m1, *k_max_values_m1, *k_min_values_m2, *k_max_values_m2;
  int           *l_min_values_fc, *l_max_values_fc, *l_min_values_fcH, *l_max_values_fcH, *l_min_values_fcI, *l_max_values_fcI, *l_min_values_fcM, *l_max_values_fcM;

  paramT        *P;

  P               = vars->P;
  sequence        = vars->sequence;
  seq_length      = vars->seq_length;
  maxD1           = vars->maxD1;
  maxD2           = vars->maxD2;
  S1              = vars->S1;
  ptype           = vars->ptype;
  reference_pt1   = vars->reference_pt1;
  reference_pt2   = vars->reference_pt2;
  my_iindx        = vars->my_iindx;
  referenceBPs1   = vars->referenceBPs1;
  referenceBPs2   = vars->referenceBPs2;
  dangles         = vars->dangles;
  mm1             = vars->mm1;
  mm2             = vars->mm2;
  bpdist          = vars->bpdist;

  E_C             = vars->E_C;
  l_min_values    = vars->l_min_values;
  l_max_values    = vars->l_max_values;
  k_min_values    = vars->k_min_values;
  k_max_values    = vars->k_max_values;

  E_M             = vars->E_M;
  l_min_values_m  = vars->l_min_values_m;
  l_max_values_m  = vars->l_max_values_m;
  k_min_values_m  = vars->k_min_values_m;
  k_max_values_m  = vars->k_max_values_m;

  E_M1            = vars->E_M1;
  l_min_values_m1 = vars->l_min_values_m1;
  l_max_values_m1 = vars->l_max_values_m1;
  k_min_values_m1 = vars->k_min_values_m1;
  k_max_values_m1 = vars->k_max_values_m1;

  E_M2            = vars->E_M2;
  l_min_values_m2 = vars->l_min_values_m2;
  l_max_values_m2 = vars->l_max_values_m2;
  k_min_values_m2 = vars->k_min_values_m2;
  k_max_values_m2 = vars->k_max_values_m2;

  E_C_rem        = vars->E_C_rem;
  E_M_rem        = vars->E_M_rem;
  E_M1_rem       = vars->E_M1_rem;
  E_M2_rem       = vars->E_M2_rem;
  E_Fc_rem       = INF;
  E_FcH_rem      = INF;
  E_FcI_rem      = INF;
  E_FcM_rem      = INF;

#ifdef _OPENMP
  #pragma omp parallel for private(d1,d2,cnt1,cnt2,cnt3,cnt4,j, i)
#endif
  for(i=1; i<seq_length-TURN-1; i++){
    /* guess memory requirements for M2 */

    int min_k, max_k, max_l, min_l;
    int min_k_real, max_k_real, *min_l_real, *max_l_real;

    min_k = min_l = 0;
    max_k = mm1[my_iindx[i]-seq_length] + referenceBPs1[my_iindx[i] - seq_length];
    max_l = mm2[my_iindx[i]-seq_length] + referenceBPs2[my_iindx[i] - seq_length];

    prepareBoundaries(min_k,
                      max_k,
                      min_l,
                      max_l,
                      bpdist[my_iindx[i] - seq_length],
                      &vars->k_min_values_m2[i],
                      &vars->k_max_values_m2[i],
                      &vars->l_min_values_m2[i],
                      &vars->l_max_values_m2[i]
                      );

    prepareArray( &vars->E_M2[i],
                  vars->k_min_values_m2[i],
                  vars->k_max_values_m2[i],
                  vars->l_min_values_m2[i],
                  vars->l_max_values_m2[i]
                 );

    preparePosteriorBoundaries( k_max_values_m2[i] - k_min_values_m2[i] + 1,
                                k_min_values_m2[i],
                                &min_k_real,
                                &max_k_real,
                                &min_l_real,
                                &max_l_real
                              );

    /* begin filling of M2 array */
    for (j=i+TURN+1; j<seq_length-TURN-1; j++){
      if(E_M1_rem[my_iindx[i]-j] != INF){
        if(E_M1[my_iindx[j+1]-seq_length])
          for(cnt1 = k_min_values_m1[my_iindx[j+1]-seq_length];
              cnt1 <= k_max_values_m1[my_iindx[j+1]-seq_length];
              cnt1++)
            for(cnt2 = l_min_values_m1[my_iindx[j+1]-seq_length][cnt1];
                cnt2 <= l_max_values_m1[my_iindx[j+1]-seq_length][cnt1];
                cnt2++)
              E_M2_rem[i] = MIN2(E_M2_rem[i],
                                  E_M1_rem[my_iindx[i]-j] + E_M1[my_iindx[j+1]-seq_length][cnt1][cnt2/2]
                                  );
        if(E_M1_rem[my_iindx[j+1]-seq_length] != INF)
          E_M2_rem[i] = MIN2(E_M2_rem[i], E_M1_rem[my_iindx[i]-j] + E_M1_rem[my_iindx[j+1]-seq_length]);
      }
      if(E_M1_rem[my_iindx[j+1]-seq_length] != INF){
        if(E_M1[my_iindx[i]-j])
          for(cnt1 = k_min_values_m1[my_iindx[i]-j];
              cnt1 <= k_max_values_m1[my_iindx[i]-j];
              cnt1++)
            for(cnt2 = l_min_values_m1[my_iindx[i]-j][cnt1];
                cnt2 <= l_max_values_m1[my_iindx[i]-j][cnt1];
                cnt2 += 2)
              E_M2_rem[i] = MIN2(E_M2_rem[i],
                                  E_M1[my_iindx[i]-j][cnt1][cnt2/2] + E_M1_rem[my_iindx[j+1]-seq_length]
                                  );
      }


      if(!E_M1[my_iindx[i]-j]) continue;
      if(!E_M1[my_iindx[j+1]-seq_length]) continue;

      d1 = referenceBPs1[my_iindx[i]-seq_length] - referenceBPs1[my_iindx[i]-j] - referenceBPs1[my_iindx[j+1]-seq_length];
      d2 = referenceBPs2[my_iindx[i]-seq_length] - referenceBPs2[my_iindx[i]-j] - referenceBPs2[my_iindx[j+1]-seq_length];

      for(cnt1 = k_min_values_m1[my_iindx[i]-j]; cnt1 <= k_max_values_m1[my_iindx[i]-j]; cnt1++)
        for(cnt2 = l_min_values_m1[my_iindx[i]-j][cnt1]; cnt2 <= l_max_values_m1[my_iindx[i]-j][cnt1]; cnt2+=2){
          for(cnt3 = k_min_values_m1[my_iindx[j+1]-seq_length]; cnt3 <= k_max_values_m1[my_iindx[j+1]-seq_length]; cnt3++)
            for(cnt4 = l_min_values_m1[my_iindx[j+1]-seq_length][cnt3]; cnt4 <= l_max_values_m1[my_iindx[j+1]-seq_length][cnt3]; cnt4+=2){
              if(((cnt1 + cnt3 + d1) <= maxD1) && ((cnt2 + cnt4 + d2) <= maxD2)){
                E_M2[i][cnt1 + cnt3 + d1][(cnt2 + cnt4 + d2)/2] = MIN2( E_M2[i][cnt1 + cnt3 + d1][(cnt2 + cnt4 + d2)/2],
                                                                        E_M1[my_iindx[i]-j][cnt1][cnt2/2] + E_M1[my_iindx[j+1]-seq_length][cnt3][cnt4/2]
                                                                      );
                updatePosteriorBoundaries(cnt1+cnt3+d1,
                                          cnt2+cnt4+d2,
                                          &min_k_real,
                                          &max_k_real,
                                          &min_l_real,
                                          &max_l_real
                                          );
              }
              else{
                E_M2_rem[i] = MIN2(E_M2_rem[i],
                                    E_M1[my_iindx[i]-j][cnt1][cnt2/2] + E_M1[my_iindx[j+1]-seq_length][cnt3][cnt4/2]
                                    );
              }
            }
        }
    }

    /* resize and move memory portions of energy matrix E_M2 */
    adjustArrayBoundaries(&vars->E_M2[i],
                          &vars->k_min_values_m2[i],
                          &vars->k_max_values_m2[i],
                          &vars->l_min_values_m2[i],
                          &vars->l_max_values_m2[i],
                          min_k_real,
                          max_k_real,
                          min_l_real,
                          max_l_real
                          );
  } /* end for i */

  base_d1 = referenceBPs1[my_iindx[1]-seq_length];
  base_d2 = referenceBPs2[my_iindx[1]-seq_length];

  /* guess memory requirements for E_FcH, E_FcI and E_FcM */

  int min_k, max_k, max_l, min_l;
  int min_k_real, max_k_real, min_k_real_fcH, max_k_real_fcH, min_k_real_fcI, max_k_real_fcI, min_k_real_fcM, max_k_real_fcM;
  int *min_l_real, *max_l_real, *min_l_real_fcH, *max_l_real_fcH, *min_l_real_fcI, *max_l_real_fcI,*min_l_real_fcM, *max_l_real_fcM;

  min_k = min_l = 0;

  max_k = mm1[my_iindx[1] - seq_length] + referenceBPs1[my_iindx[1] - seq_length];
  max_l = mm2[my_iindx[1] - seq_length] + referenceBPs2[my_iindx[1] - seq_length];

#ifdef _OPENMP
  #pragma omp sections
  {

  #pragma omp section
  {
#endif
  prepareBoundaries(min_k,
                    max_k,
                    min_l,
                    max_l,
                    bpdist[my_iindx[1] - seq_length],
                    &vars->k_min_values_fc,
                    &vars->k_max_values_fc,
                    &vars->l_min_values_fc,
                    &vars->l_max_values_fc
                    );
  prepareArray( &vars->E_Fc,
                vars->k_min_values_fc,
                vars->k_max_values_fc,
                vars->l_min_values_fc,
                vars->l_max_values_fc
              );
#ifdef _OPENMP
  }
  #pragma omp section
  {
#endif
  prepareBoundaries(min_k,
                    max_k,
                    min_l,
                    max_l,
                    bpdist[my_iindx[1] - seq_length],
                    &vars->k_min_values_fcH,
                    &vars->k_max_values_fcH,
                    &vars->l_min_values_fcH,
                    &vars->l_max_values_fcH
                    );
  prepareArray( &vars->E_FcH,
                vars->k_min_values_fcH,
                vars->k_max_values_fcH,
                vars->l_min_values_fcH,
                vars->l_max_values_fcH
              );
#ifdef _OPENMP
  }
  #pragma omp section
  {
#endif
  prepareBoundaries(min_k,
                    max_k,
                    min_l,
                    max_l,
                    bpdist[my_iindx[1] - seq_length],
                    &vars->k_min_values_fcI,
                    &vars->k_max_values_fcI,
                    &vars->l_min_values_fcI,
                    &vars->l_max_values_fcI
                    );
  prepareArray( &vars->E_FcI,
                vars->k_min_values_fcI,
                vars->k_max_values_fcI,
                vars->l_min_values_fcI,
                vars->l_max_values_fcI
              );
#ifdef _OPENMP
  }
  #pragma omp section
  {
#endif
  prepareBoundaries(min_k,
                    max_k,
                    min_l,
                    max_l,
                    bpdist[my_iindx[1] - seq_length],
                    &vars->k_min_values_fcM,
                    &vars->k_max_values_fcM,
                    &vars->l_min_values_fcM,
                    &vars->l_max_values_fcM
                    );
  prepareArray( &vars->E_FcM,
                vars->k_min_values_fcM,
                vars->k_max_values_fcM,
                vars->l_min_values_fcM,
                vars->l_max_values_fcM
              );
#ifdef _OPENMP
  }
  #pragma omp section
  {
#endif
  preparePosteriorBoundaries( max_k - min_k + 1,
                              min_k,
                              &min_k_real,
                              &max_k_real,
                              &min_l_real,
                              &max_l_real
                            );
#ifdef _OPENMP
  }
  #pragma omp section
  {
#endif
  preparePosteriorBoundaries( max_k - min_k + 1,
                              min_k,
                              &min_k_real_fcH,
                              &max_k_real_fcH,
                              &min_l_real_fcH,
                              &max_l_real_fcH
                            );
#ifdef _OPENMP
  }
  #pragma omp section
  {
#endif
  preparePosteriorBoundaries( max_k - min_k + 1,
                              min_k,
                              &min_k_real_fcI,
                              &max_k_real_fcI,
                              &min_l_real_fcI,
                              &max_l_real_fcI
                            );

#ifdef _OPENMP
  }
  #pragma omp section
  {
#endif
  preparePosteriorBoundaries( max_k - min_k + 1,
                              min_k,
                              &min_k_real_fcM,
                              &max_k_real_fcM,
                              &min_l_real_fcM,
                              &max_l_real_fcM
                            );
#ifdef _OPENMP
  }
  }
#endif


  E_Fc              = vars->E_Fc;
  E_FcH             = vars->E_FcH;
  E_FcI             = vars->E_FcI;
  E_FcM             = vars->E_FcM;

  l_min_values_fc   = vars->l_min_values_fc;
  l_max_values_fc   = vars->l_max_values_fc;
  l_min_values_fcH  = vars->l_min_values_fcH;
  l_max_values_fcH  = vars->l_max_values_fcH;
  l_min_values_fcI  = vars->l_min_values_fcI;
  l_max_values_fcI  = vars->l_max_values_fcI;
  l_min_values_fcM  = vars->l_min_values_fcM;
  l_max_values_fcM  = vars->l_max_values_fcM;

  /* begin actual energy calculations */
#ifdef _OPENMP
  #pragma omp sections private(d, d1,d2,cnt1,cnt2,cnt3,cnt4,j, i, energy)
  {

  #pragma omp section
  {
#endif
  for (d = TURN+2; d <= seq_length; d++) /* i,j in [1..length] */
    for (j = d; j <= seq_length; j++) {
      unsigned int  u, ij, p, q, pq;
      int           type, type_2, no_close;
      char          loopseq[10];
      i = j-d+1;
      ij = my_iindx[i]-j;
      u = seq_length-j + i-1;
      if (u<TURN) continue;

      type = ptype[ij];

      no_close = (((type==3)||(type==4))&&no_closingGU);

      type=rtype[type];

      if (!type) continue;
      if(no_close) continue;

      d1 = base_d1 - referenceBPs1[ij];
      d2 = base_d2 - referenceBPs2[ij];
      if (u<7) {
        strcpy(loopseq , sequence+j-1);
        strncat(loopseq, sequence, i);
      }
      energy = E_Hairpin(u, type, S1[j+1], S1[i-1],  loopseq, P);

      if(E_C_rem[ij] != INF)
        E_FcH_rem = MIN2(E_FcH_rem, E_C_rem[ij] + energy);

      if (!E_C[ij]) continue;
      for(cnt1 = k_min_values[ij]; cnt1 <= k_max_values[ij]; cnt1++)
        for(cnt2 = l_min_values[ij][cnt1]; cnt2 <= l_max_values[ij][cnt1]; cnt2 += 2){
          if(((cnt1 + d1) <= maxD1) && ((cnt2 + d2) <= maxD2)){
            E_FcH[cnt1 + d1][(cnt2+d2)/2] = MIN2( E_FcH[cnt1 + d1][(cnt2+d2)/2],
                                                  energy + E_C[ij][cnt1][cnt2/2]
                                                );
            updatePosteriorBoundaries(cnt1 + d1,
                                      cnt2 + d2,
                                      &min_k_real_fcH,
                                      &max_k_real_fcH,
                                      &min_l_real_fcH,
                                      &max_l_real_fcH
                                    );
          }
          else
            E_FcH_rem = MIN2(E_FcH_rem, energy + E_C[ij][cnt1][cnt2/2]);
        }
    }
  /* end of i-j loop */

  /* resize and move memory portions of energy matrix E_FcH */
  adjustArrayBoundaries(&vars->E_FcH,
                        &vars->k_min_values_fcH,
                        &vars->k_max_values_fcH,
                        &vars->l_min_values_fcH,
                        &vars->l_max_values_fcH,
                        min_k_real_fcH,
                        max_k_real_fcH,
                        min_l_real_fcH,
                        max_l_real_fcH
                      );
#ifdef _OPENMP
  }
  #pragma omp section
  {
#endif
  for (d = TURN+2; d <= seq_length; d++) /* i,j in [1..length] */
    for (j = d; j <= seq_length; j++) {
      unsigned int  u, ij, p, q, pq;
      int           type, type_2, no_close;
      i = j-d+1;
      ij = my_iindx[i]-j;
      u = seq_length-j + i-1;
      if (u<TURN) continue;

      type = ptype[ij];

      no_close = (((type==3)||(type==4))&&no_closingGU);

      type=rtype[type];

      if (!type) continue;
      if(no_close) continue;

      if(E_C_rem[ij] != INF){
        for(p = j+1; p < seq_length ; p++){
          unsigned int u1, qmin, ln_pre;
          u1 = p-j-1;
          if (u1+i-1>MAXLOOP) break;
          qmin = p + TURN + 1;
          ln_pre = u1 + i + seq_length;
          if(ln_pre > qmin + MAXLOOP) qmin = ln_pre - MAXLOOP - 1;
          for(q = qmin; q <= seq_length; q++){
            unsigned int u2;
            pq = my_iindx[p]-q;
            type_2 = rtype[ptype[pq]];
            if (type_2==0) continue;
            u2 = i-1 + seq_length-q;
            if (u1+u2>MAXLOOP) continue;
            /* get distance to reference if closing the interior loop
            *  d2a = dbp(T1_[1,n}, T1_{p,q} + T1_{i,j})
            *  d2b = dbp(T2_[1,n}, T2_{p,q} + T2_{i,j})
            */
            d1 = base_d1 - referenceBPs1[ij] - referenceBPs1[pq];
            d2 = base_d2 - referenceBPs2[ij] - referenceBPs2[pq];
            energy = E_IntLoop(u1, u2, type, type_2, S1[j+1], S1[i-1], S1[p-1], S1[q+1], P);

            if(E_C_rem[pq] != INF)
              E_FcI_rem = MIN2(E_FcI_rem, E_C_rem[ij] + E_C_rem[pq] + energy);

            if(E_C[pq])
              for(cnt1 = k_min_values[pq];
                  cnt1 <= k_max_values[pq];
                  cnt1++)
                for(cnt2 = l_min_values[pq][cnt1];
                    cnt2 <= l_max_values[pq][cnt1];
                    cnt2 += 2)
                  E_FcI_rem = MIN2(E_FcI_rem, E_C_rem[ij] + E_C[pq][cnt1][cnt2/2] + energy);
          }
        }
      }

      if(E_C[ij]){
        for(p = j+1; p < seq_length ; p++){
          unsigned int u1, qmin, ln_pre;
          u1 = p-j-1;
          if (u1+i-1>MAXLOOP) break;
          qmin = p + TURN + 1;
          ln_pre = u1 + i + seq_length;
          if(ln_pre > qmin + MAXLOOP) qmin = ln_pre - MAXLOOP - 1;
          for(q = qmin; q <= seq_length; q++){
            unsigned int u2;
            pq = my_iindx[p]-q;
            type_2 = rtype[ptype[pq]];
            if (type_2==0) continue;
            u2 = i-1 + seq_length-q;
            if (u1+u2>MAXLOOP) continue;
            /* get distance to reference if closing the interior loop
            *  d2a = dbp(T1_[1,n}, T1_{p,q} + T1_{i,j})
            *  d2b = dbp(T2_[1,n}, T2_{p,q} + T2_{i,j})
            */
            d1 = base_d1 - referenceBPs1[ij] - referenceBPs1[pq];
            d2 = base_d2 - referenceBPs2[ij] - referenceBPs2[pq];
            energy = E_IntLoop(u1, u2, type, type_2, S1[j+1], S1[i-1], S1[p-1], S1[q+1], P);
            if(E_C_rem[pq] != INF){
              for(cnt1 = k_min_values[ij];
                  cnt1 <= k_max_values[ij];
                  cnt1++)
                for(cnt2 = l_min_values[ij][cnt1];
                    cnt2 <= l_max_values[ij][cnt1];
                    cnt2 += 2)
                  E_FcI_rem = MIN2(E_FcI_rem, E_C[ij][cnt1][cnt2/2] + E_C_rem[pq] + energy);
            }

            if(E_C[pq])
              for(cnt1 = k_min_values[ij];
                  cnt1 <= k_max_values[ij];
                  cnt1++)
                for(cnt2 = l_min_values[ij][cnt1];
                    cnt2 <= l_max_values[ij][cnt1];
                    cnt2 += 2)
                  for(cnt3 = k_min_values[pq];
                      cnt3 <= k_max_values[pq];
                      cnt3++)
                    for(cnt4 = l_min_values[pq][cnt3];
                        cnt4 <= l_max_values[pq][cnt3];
                        cnt4 += 2){
                      if(((cnt1 + cnt3 + d1) <= maxD1) && ((cnt2 + cnt4 + d2) <= maxD2)){
                        E_FcI[cnt1 + cnt3 + d1][(cnt2 + cnt4 + d2)/2] = MIN2(
                                                                          E_FcI[cnt1 + cnt3 + d1][(cnt2 + cnt4 + d2)/2],
                                                                          E_C[ij][cnt1][cnt2/2]
                                                                        + E_C[pq][cnt3][cnt4/2]
                                                                        + energy
                                                                            );
                        updatePosteriorBoundaries(cnt1 + cnt3 + d1,
                                                  cnt2 + cnt4 + d2,
                                                  &min_k_real_fcI,
                                                  &max_k_real_fcI,
                                                  &min_l_real_fcI,
                                                  &max_l_real_fcI
                                                );
                      }
                      else{
                        E_FcI_rem = MIN2(
                                      E_FcI_rem,
                                      E_C[ij][cnt1][cnt2/2]
                                    + E_C[pq][cnt3][cnt4/2]
                                    + energy
                                          );
                      }
                    }
          }
        }
      }
    }
  /* end of i-j loop */

  /* resize and move memory portions of energy matrix E_FcI */
  adjustArrayBoundaries(&vars->E_FcI,
                        &vars->k_min_values_fcI,
                        &vars->k_max_values_fcI,
                        &vars->l_min_values_fcI,
                        &vars->l_max_values_fcI,
                        min_k_real_fcI,
                        max_k_real_fcI,
                        min_l_real_fcI,
                        max_l_real_fcI
                      );
#ifdef _OPENMP
  }
  #pragma omp section
  {
#endif
  if(seq_length > 2*TURN){
    for (i=TURN+1; i<seq_length-2*TURN; i++) {
      /* get distancies to references
      * d3a = dbp(T1_[1,n}, T1_{1,k} + T1_{k+1, n})
      * d3b = dbp(T2_[1,n}, T2_{1,k} + T2_{k+1, n})
      */
      d1 = base_d1 - referenceBPs1[my_iindx[1]-i] - referenceBPs1[my_iindx[i+1]-seq_length];
      d2 = base_d2 - referenceBPs2[my_iindx[1]-i] - referenceBPs2[my_iindx[i+1]-seq_length];

      if(E_M_rem[my_iindx[1]-i] != INF){
        if(E_M2[i+1])
          for(cnt1 = k_min_values_m2[i+1];
              cnt1 <= k_max_values_m2[i+1];
              cnt1++)
            for(cnt2 = l_min_values_m2[i+1][cnt1];
                cnt2 <= l_max_values_m2[i+1][cnt1];
                cnt2 += 2)
              E_FcM_rem = MIN2(E_FcM_rem, E_M_rem[my_iindx[1]-i] + E_M2[i+1][cnt1][cnt2/2] + P->MLclosing);
        if(E_M2_rem[i+1] != INF)
          E_FcM_rem = MIN2(E_FcM_rem, E_M_rem[my_iindx[1]-i] + E_M2_rem[i+1] + P->MLclosing);
      }
      if(E_M2_rem[i+1] != INF){
        if(E_M[my_iindx[1]-i])
          for(cnt1 = k_min_values_m[my_iindx[1]-i];
              cnt1 <= k_max_values_m[my_iindx[1]-i];
              cnt1++)
            for(cnt2 = l_min_values_m[my_iindx[1]-i][cnt1];
                cnt2 <= l_max_values_m[my_iindx[1]-i][cnt1];
                cnt2 += 2)
              E_FcM_rem = MIN2(E_FcM_rem, E_M[my_iindx[1]-i][cnt1][cnt2/2] + E_M2_rem[i+1] + P->MLclosing);
      }

      if(!E_M[my_iindx[1]-i]) continue;
      if(!E_M2[i+1]) continue;
      for(cnt1 = k_min_values_m[my_iindx[1]-i]; cnt1 <= k_max_values_m[my_iindx[1]-i]; cnt1++)
        for(cnt2 = l_min_values_m[my_iindx[1]-i][cnt1]; cnt2 <= l_max_values_m[my_iindx[1]-i][cnt1]; cnt2 += 2)
          for(cnt3 = k_min_values_m2[i+1]; cnt3 <= k_max_values_m2[i+1]; cnt3++)
            for(cnt4 = l_min_values_m2[i+1][cnt3]; cnt4 <= l_max_values_m2[i+1][cnt3]; cnt4 += 2){
              if(((cnt1 + cnt3 + d1) <= maxD1) && ((cnt2 + cnt4 + d2) <= maxD2)){
                E_FcM[cnt1 + cnt3 + d1][(cnt2 + cnt4 + d2)/2] = MIN2(
                                                                  E_FcM[cnt1 + cnt3 + d1][(cnt2 + cnt4 + d2)/2],
                                                                  E_M[my_iindx[1]-i][cnt1][cnt2/2]
                                                                + E_M2[i+1][cnt3][cnt4/2]
                                                                + P->MLclosing
                                                                    );
                updatePosteriorBoundaries(cnt1 + cnt3 + d1,
                                          cnt2 + cnt4 + d2,
                                          &min_k_real_fcM,
                                          &max_k_real_fcM,
                                          &min_l_real_fcM,
                                          &max_l_real_fcM
                                        );
              }
              else{
                E_FcM_rem = MIN2(
                              E_FcM_rem,
                              E_M[my_iindx[1]-i][cnt1][cnt2/2]
                            + E_M2[i+1][cnt3][cnt4/2]
                            + P->MLclosing
                                  );
              }
            }
    }
  }
  /* resize and move memory portions of energy matrix E_FcM */
  adjustArrayBoundaries(&vars->E_FcM,
                        &vars->k_min_values_fcM,
                        &vars->k_max_values_fcM,
                        &vars->l_min_values_fcM,
                        &vars->l_max_values_fcM,
                        min_k_real_fcM,
                        max_k_real_fcM,
                        min_l_real_fcM,
                        max_l_real_fcM
                      );
#ifdef _OPENMP
  }
  }
#endif



  /* compute E_Fc_rem */
  E_Fc_rem = MIN2(E_FcH_rem, E_FcI_rem);
  E_Fc_rem = MIN2(E_Fc_rem, E_FcM_rem);
  /* add the case were structure is unfolded chain */
  if((referenceBPs1[my_iindx[1]-seq_length] > maxD1) || (referenceBPs2[my_iindx[1]-seq_length] > maxD2))
    E_Fc_rem = MIN2(E_Fc_rem, 0);

  /* store values in ds */
  vars->E_Fc_rem   = E_Fc_rem;
  vars->E_FcH_rem  = E_FcH_rem;
  vars->E_FcI_rem  = E_FcI_rem;
  vars->E_FcM_rem  = E_FcM_rem;


  /* compute all E_Fc */
  for(cnt1 = vars->k_min_values_fcH; cnt1 <= vars->k_max_values_fcH; cnt1++)
    for(cnt2 = vars->l_min_values_fcH[cnt1]; cnt2 <= vars->l_max_values_fcH[cnt1]; cnt2 += 2){
      vars->E_Fc[cnt1][cnt2/2] = MIN2(vars->E_Fc[cnt1][cnt2/2],
                                      vars->E_FcH[cnt1][cnt2/2]
                                      );
      updatePosteriorBoundaries(cnt1,
                                cnt2,
                                &min_k_real,
                                &max_k_real,
                                &min_l_real,
                                &max_l_real
                                );
    }
  for(cnt1 = vars->k_min_values_fcI; cnt1 <= vars->k_max_values_fcI; cnt1++)
    for(cnt2 = vars->l_min_values_fcI[cnt1]; cnt2 <= vars->l_max_values_fcI[cnt1]; cnt2 += 2){
      vars->E_Fc[cnt1][cnt2/2] = MIN2(vars->E_Fc[cnt1][cnt2/2],
                                      vars->E_FcI[cnt1][cnt2/2]
                                      );
      updatePosteriorBoundaries(cnt1,
                                cnt2,
                                &min_k_real,
                                &max_k_real,
                                &min_l_real,
                                &max_l_real
                                );
    }
  for(cnt1 = vars->k_min_values_fcM; cnt1 <= vars->k_max_values_fcM; cnt1++)
    for(cnt2 = vars->l_min_values_fcM[cnt1]; cnt2 <= vars->l_max_values_fcM[cnt1]; cnt2 += 2){
      vars->E_Fc[cnt1][cnt2/2] = MIN2(vars->E_Fc[cnt1][cnt2/2],
                                      vars->E_FcM[cnt1][cnt2/2]
                                      );
      updatePosteriorBoundaries(cnt1,
                                cnt2,
                                &min_k_real,
                                &max_k_real,
                                &min_l_real,
                                &max_l_real
                                );
    }
  /* add the case were structure is unfolded chain */
  vars->E_Fc[referenceBPs1[my_iindx[1]-seq_length]][referenceBPs2[my_iindx[1]-seq_length]/2] = MIN2(vars->E_Fc[referenceBPs1[my_iindx[1]-seq_length]][referenceBPs2[my_iindx[1]-seq_length]/2],
                                                                                                    0);
  updatePosteriorBoundaries(referenceBPs1[my_iindx[1]-seq_length],
                            referenceBPs2[my_iindx[1]-seq_length],
                            &min_k_real,
                            &max_k_real,
                            &min_l_real,
                            &max_l_real
                            );


  adjustArrayBoundaries(&vars->E_Fc,
                        &vars->k_min_values_fc,
                        &vars->k_max_values_fc,
                        &vars->l_min_values_fc,
                        &vars->l_max_values_fc,
                        min_k_real,
                        max_k_real,
                        min_l_real,
                        max_l_real
                      );

}




PRIVATE void adjustArrayBoundaries(int ***array, int *k_min, int *k_max, int **l_min, int **l_max, int k_min_post, int k_max_post, int *l_min_post, int *l_max_post){
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
      memmove((int **)(*array),((int **)(*array)) + k_diff_pre, sizeof(int *) * mem_size);
      memmove((int *) (*l_min),((int *) (*l_min)) + k_diff_pre, sizeof(int)   * mem_size);
      memmove((int *) (*l_max),((int *) (*l_max)) + k_diff_pre, sizeof(int)   * mem_size);
    }

    /* reallocating memory to actual size used */
    *array  += *k_min;
    *array = (int **)realloc(*array, sizeof(int *) * mem_size);
    *array -= k_min_post;

    *l_min  += *k_min;
    *l_min = (int *)realloc(*l_min, sizeof(int) * mem_size);
    *l_min -= k_min_post;

    *l_max  += *k_min;
    *l_max = (int *)realloc(*l_max, sizeof(int) * mem_size);
    *l_max -= k_min_post;

    /* adjust l dimension of array */
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
          memmove((int *)((*array)[cnt1]), (int *)((*array)[cnt1])+start, sizeof(int) * mem_size);
        (*array)[cnt1]  = (int *) realloc((*array)[cnt1], sizeof(int) * mem_size);

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

INLINE  PRIVATE void  prepareBoundaries(int min_k_pre, int max_k_pre, int min_l_pre, int max_l_pre, int bpdist, int *min_k, int *max_k, int **min_l, int **max_l){
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

INLINE  PRIVATE void  prepareArray(int ***array, int min_k, int max_k, int *min_l, int *max_l){
  int i, j, mem;
  *array  = (int **)space(sizeof(int *) * (max_k - min_k + 1));
  *array  -= min_k;

  for(i = min_k; i <= max_k; i++){
    mem = (max_l[i] - min_l[i] + 1)/2 + 1;
    (*array)[i] = (int *)space(sizeof(int) * mem);
    for(j = 0; j < mem; j++)
      (*array)[i][j] = INF;
    (*array)[i]  -= min_l[i]/2;
  }
}

INLINE  PRIVATE void  prepareArray2(unsigned long ***array, int min_k, int max_k, int *min_l, int *max_l){
  int i, j, mem;
  *array  = (unsigned long **)space(sizeof(unsigned long *) * (max_k - min_k + 1));
  *array  -= min_k;

  for(i = min_k; i <= max_k; i++){
    mem = (max_l[i] - min_l[i] + 1)/2 + 1;
    (*array)[i] = (unsigned long *)space(sizeof(unsigned long) * mem);
    (*array)[i] -= min_l[i]/2;
  }
}

