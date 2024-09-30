/** \file dp_matrices.c **/

/*
 *                Dynamic Programming matrix related functions
 *
 *                This file contains everything necessary to
 *                obtain and destroy data structures representing
 *                dynamic programming (DP) matrices used in the folding
 *                recurrences throughout the VienneRNA paclage
 *
 *                c Ronny Lorenz
 *
 *                ViennaRNA package
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "ViennaRNA/datastructures/basic.h"
#include "ViennaRNA/model.h"
#include "ViennaRNA/utils/basic.h"
#include "ViennaRNA/utils/log.h"
#include "ViennaRNA/mfe/gquad.h"
#include "ViennaRNA/datastructures/dp_matrices.h"

/*
 #################################
 # PRIVATE MACROS                #
 #################################
 */

/* the definitions below indicate which arrays should be allocated upon retrieval of a matrices data structure */
#define ALLOC_NOTHING     0
#define ALLOC_F           1
#define ALLOC_F5          2
#define ALLOC_F3          4
#define ALLOC_FC          8
#define ALLOC_C           16
#define ALLOC_FML         32
#define ALLOC_PROBS       256
#define ALLOC_AUX         512

#define ALLOC_CIRC        1024
#define ALLOC_MULTISTRAND 2048
#define ALLOC_UNIQ        4096

#define ALLOC_MFE_DEFAULT         (ALLOC_F5 | ALLOC_C | ALLOC_FML)
#define ALLOC_MFE_LOCAL           (ALLOC_F3 | ALLOC_C | ALLOC_FML)

#define ALLOC_PF_WO_PROBS         (ALLOC_F | ALLOC_C | ALLOC_FML)
#define ALLOC_PF_DEFAULT          (ALLOC_PF_WO_PROBS | ALLOC_PROBS | ALLOC_AUX)

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
PRIVATE unsigned int
get_mx_alloc_vector(vrna_fold_compound_t  *fc,
                    vrna_mx_type_e        type,
                    unsigned int          options);


PRIVATE unsigned int
get_mx_mfe_alloc_vector_current(vrna_mx_mfe_t   *mx,
                                vrna_mx_type_e  mx_type);


PRIVATE unsigned int
get_mx_pf_alloc_vector_current(vrna_mx_pf_t   *mx,
                               vrna_mx_type_e mx_type);


PRIVATE void
mfe_matrices_free_default(vrna_mx_mfe_t *self);


PRIVATE void
mfe_matrices_free_window(vrna_mx_mfe_t  *self,
                         unsigned int   length,
                         unsigned int   window_size);


PRIVATE void
mfe_matrices_free_2Dfold(vrna_mx_mfe_t  *self,
                         unsigned int   length,
                         int            min_loop_size,
                         int            *indx);


PRIVATE void
pf_matrices_free_default(vrna_mx_pf_t *self);


PRIVATE void
pf_matrices_free_window(vrna_mx_pf_t  *self,
                        unsigned int  length,
                        unsigned int  window_size);


PRIVATE void
pf_matrices_free_2Dfold(vrna_mx_pf_t  *self,
                        unsigned int  length,
                        int           turn,
                        int           *indx,
                        int           *jindx);


PRIVATE int
add_pf_matrices(vrna_fold_compound_t  *vc,
                vrna_mx_type_e        type,
                unsigned int          alloc_vector);


PRIVATE int
add_mfe_matrices(vrna_fold_compound_t *vc,
                 vrna_mx_type_e       type,
                 unsigned int         alloc_vector);


PRIVATE vrna_mx_mfe_t *
init_mx_mfe_default(vrna_fold_compound_t  *fc,
                    unsigned int          alloc_vector);


PRIVATE vrna_mx_mfe_t *
init_mx_mfe_window(vrna_fold_compound_t *fc,
                   unsigned int         alloc_vector);


PRIVATE vrna_mx_mfe_t *
init_mx_mfe_2Dfold(vrna_fold_compound_t *fc,
                   unsigned int         alloc_vector);


PRIVATE INLINE void
nullify_mfe(vrna_mx_mfe_t *mx);


PRIVATE vrna_mx_pf_t *
init_mx_pf_default(vrna_fold_compound_t *fc,
                   unsigned int         alloc_vector);


PRIVATE vrna_mx_pf_t *
init_mx_pf_window(vrna_fold_compound_t  *fc,
                  unsigned int          alloc_vector);


PRIVATE vrna_mx_pf_t *
init_mx_pf_2Dfold(vrna_fold_compound_t  *fc,
                  unsigned int          alloc_vector);


PRIVATE INLINE void
nullify_pf(vrna_mx_pf_t *mx);


/*
 #################################
 # BEGIN OF FUNCTION DEFINITIONS #
 #################################
 */
PUBLIC void
vrna_mx_mfe_free(vrna_fold_compound_t *vc)
{
  if (vc) {
    vrna_mx_mfe_t *self = vc->matrices;
    if (self) {
      switch (self->type) {
        case VRNA_MX_DEFAULT:
          mfe_matrices_free_default(self);
          break;

        case VRNA_MX_WINDOW:
          mfe_matrices_free_window(self, vc->length, vc->window_size);
          break;

        case VRNA_MX_2DFOLD:
          mfe_matrices_free_2Dfold(self,
                                   vc->length,
                                   vc->params->model_details.min_loop_size,
                                   vc->iindx);
          break;

        default:                /* do nothing */
          break;
      }
      free(self);
      vc->matrices = NULL;
    }
  }
}


PUBLIC void
vrna_mx_pf_free(vrna_fold_compound_t *vc)
{
  if (vc) {
    vrna_mx_pf_t *self = vc->exp_matrices;
    if (self) {
      switch (self->type) {
        case VRNA_MX_DEFAULT:
          pf_matrices_free_default(self);
          break;

        case VRNA_MX_WINDOW:
          pf_matrices_free_window(self, vc->length, vc->window_size);
          break;

        case VRNA_MX_2DFOLD:
          pf_matrices_free_2Dfold(self,
                                  vc->length,
                                  vc->exp_params->model_details.min_loop_size,
                                  vc->iindx,
                                  vc->jindx);
          break;

        default:                /* do nothing */
          break;
      }

      free(self->expMLbase);
      free(self->scale);

      free(self);
      vc->exp_matrices = NULL;
    }
  }
}


PUBLIC int
vrna_mx_add(vrna_fold_compound_t  *vc,
            vrna_mx_type_e        mx_type,
            unsigned int          options)
{
  int ret;

  ret = 1;

  if (options & VRNA_OPTION_MFE)
    ret &= vrna_mx_mfe_add(vc, mx_type, options);

  if (options & VRNA_OPTION_PF)
    ret &= vrna_mx_pf_add(vc, mx_type, options);

  return ret;
}


PUBLIC int
vrna_mx_mfe_add(vrna_fold_compound_t  *vc,
                vrna_mx_type_e        mx_type,
                unsigned int          options)
{
  unsigned int mx_alloc_vector;

  if (vc->params) {
    options |= VRNA_OPTION_MFE;

    mx_alloc_vector = get_mx_alloc_vector(vc,
                                          mx_type,
                                          options);
    vrna_mx_mfe_free(vc);
    return add_mfe_matrices(vc, mx_type, mx_alloc_vector);
  }

  return 0;
}


PUBLIC int
vrna_mx_pf_add(vrna_fold_compound_t *vc,
               vrna_mx_type_e       mx_type,
               unsigned int         options)
{
  unsigned int mx_alloc_vector;

  if (vc->exp_params) {
    mx_alloc_vector = get_mx_alloc_vector(vc,
                                          mx_type,
                                          options | VRNA_OPTION_PF);
    vrna_mx_pf_free(vc);
    return add_pf_matrices(vc, mx_type, mx_alloc_vector);
  }

  return 0;
}


PUBLIC int
vrna_mx_prepare(vrna_fold_compound_t  *vc,
                unsigned int          options)
{
  int             ret, realloc;
  unsigned int    mx_alloc_vector, mx_alloc_vector_current;
  vrna_mx_type_e  mx_type;

  ret = 1;

  if (vc) {
    /*  check whether we have the correct DP matrices attached, and if there is
     *  enough memory allocated
     */
    if (options & VRNA_OPTION_MFE) {
      /* prepare for MFE computation */
      if (options & VRNA_OPTION_WINDOW) /* Windowing approach, a.k.a. locally optimal */
        mx_type = VRNA_MX_WINDOW;
      else                              /* default is regular MFE */
        mx_type = VRNA_MX_DEFAULT;

      if (vc->strands > 1)
        options |= VRNA_OPTION_HYBRID;

      realloc = 0;

      if (!vc->matrices || (vc->matrices->type != mx_type) || (vc->matrices->length < vc->length)) {
        realloc = 1;
      } else {
        mx_alloc_vector =
          get_mx_alloc_vector(vc, mx_type, options);
        mx_alloc_vector_current = get_mx_mfe_alloc_vector_current(vc->matrices, mx_type);
        if ((mx_alloc_vector & mx_alloc_vector_current) != mx_alloc_vector)
          realloc = 1;
      }

      if (realloc) /* Add DP matrices, if not they are not present */
        ret &= vrna_mx_mfe_add(vc, mx_type, options);
    }

    if (options & VRNA_OPTION_PF) {
      /* prepare for partition function computations */
      if (!vc->exp_params) /* return failure if exp_params data is not present */
        return 0;

      if (options & VRNA_OPTION_WINDOW) /* Windowing approach, a.k.a. locally optimal */
        mx_type = VRNA_MX_WINDOW;
      else                              /* default is regular MFE */
        mx_type = VRNA_MX_DEFAULT;

      if (vc->strands > 1)
        options |= VRNA_OPTION_HYBRID;

      realloc = 0;

      /*  Add DP matrices, if not they are not present */
      if (!vc->exp_matrices || (vc->exp_matrices->type != mx_type) ||
          (vc->exp_matrices->length < vc->length)) {
        realloc = 1;
      } else {
        mx_alloc_vector = get_mx_alloc_vector(vc,
                                              mx_type,
                                              options);
        mx_alloc_vector_current = get_mx_pf_alloc_vector_current(vc->exp_matrices, mx_type);
        if ((mx_alloc_vector & mx_alloc_vector_current) != mx_alloc_vector)
          realloc = 1;
      }

      if (realloc) /* Add DP matrices, if not they are not present */
        ret &= vrna_mx_pf_add(vc, mx_type, options);

#ifndef VRNA_DISABLE_BACKWARD_COMPATIBILITY
      else   /* re-compute pf_scale and MLbase contributions (for RNAup)*/
        vrna_exp_params_rescale(vc, NULL);

#endif
    }
  } else {
    ret = 0;
  }

  return ret;
}


/*
 #####################################
 # BEGIN OF STATIC HELPER FUNCTIONS  #
 #####################################
 */
PRIVATE unsigned int
get_mx_mfe_alloc_vector_current(vrna_mx_mfe_t   *mx,
                                vrna_mx_type_e  mx_type)
{
  unsigned int mx_alloc_vector = ALLOC_NOTHING;

  if (mx) {
    switch (mx_type) {
      case VRNA_MX_DEFAULT:
        if (mx->f5)
          mx_alloc_vector |= ALLOC_F5;

        if (mx->f3)
          mx_alloc_vector |= ALLOC_F3;

        if ((mx->fms5) ||
            (mx->fms3))
          mx_alloc_vector |= ALLOC_MULTISTRAND;

        if (mx->c)
          mx_alloc_vector |= ALLOC_C;

        if (mx->fML)
          mx_alloc_vector |= ALLOC_FML;

        if ((mx->fM1) ||
            (mx->fM2_real))
          mx_alloc_vector |= ALLOC_UNIQ;

        if ((mx->fM1_new) &&
            (mx->fM2_real))
          mx_alloc_vector |= ALLOC_CIRC;

        break;

      default:
        break;
    }
  }

  return mx_alloc_vector;
}


PRIVATE unsigned int
get_mx_pf_alloc_vector_current(vrna_mx_pf_t   *mx,
                               vrna_mx_type_e mx_type)
{
  unsigned int mx_alloc_vector = ALLOC_NOTHING;

  if (mx) {
    switch (mx_type) {
      case VRNA_MX_DEFAULT:
        if (mx->q)
          mx_alloc_vector |= ALLOC_F;

        if (mx->qb)
          mx_alloc_vector |= ALLOC_C;

        if (mx->qm)
          mx_alloc_vector |= ALLOC_FML;

        if ((mx->qm1) ||
            (mx->qm2_real))
          mx_alloc_vector |= ALLOC_UNIQ;

        if ((mx->qm1_new) &&
            (mx->qm2_real))
          mx_alloc_vector |= ALLOC_CIRC;

        if (mx->probs)
          mx_alloc_vector |= ALLOC_PROBS;

        if (mx->q1k && mx->qln)
          mx_alloc_vector |= ALLOC_AUX;

        break;

      default:
        break;
    }
  }

  return mx_alloc_vector;
}


PRIVATE int
add_pf_matrices(vrna_fold_compound_t  *vc,
                vrna_mx_type_e        mx_type,
                unsigned int          alloc_vector)
{
  if (vc) {
    switch (mx_type) {
      case VRNA_MX_DEFAULT:
        vc->exp_matrices = init_mx_pf_default(vc, alloc_vector);
        break;

      case VRNA_MX_WINDOW:
        vc->exp_matrices = init_mx_pf_window(vc, alloc_vector);
        break;

      case VRNA_MX_2DFOLD:
        vc->exp_matrices = init_mx_pf_2Dfold(vc, alloc_vector);
        break;

      default:                /* do nothing */
        return 0;
    }

    if (!vc->exp_matrices)
      return 0;

    vrna_exp_params_rescale(vc, NULL);
  }

  return 1;
}


PRIVATE int
add_mfe_matrices(vrna_fold_compound_t *vc,
                 vrna_mx_type_e       mx_type,
                 unsigned int         alloc_vector)
{
  if (vc) {
    switch (mx_type) {
      case VRNA_MX_DEFAULT:
        vc->matrices = init_mx_mfe_default(vc, alloc_vector);
        break;

      case VRNA_MX_WINDOW:
        vc->matrices = init_mx_mfe_window(vc, alloc_vector);
        break;

      case VRNA_MX_2DFOLD:
        vc->matrices = init_mx_mfe_2Dfold(vc, alloc_vector);
        break;

      default:
        return 0;
    }

    if (!vc->matrices)
      return 0;

    if (vc->params->model_details.gquad) {
      if (mx_type != VRNA_MX_WINDOW)
        vc->matrices->c_gq = vrna_mfe_gquad_mx(vc);
    }
  }

  return 1;
}


PRIVATE unsigned int
get_mx_alloc_vector(vrna_fold_compound_t  *fc,
                    vrna_mx_type_e        mx_type,
                    unsigned int          options)
{
  unsigned int  v;
  vrna_md_t     *md_p;

  md_p = &(fc->params->model_details);

  v = ALLOC_NOTHING;

  /* default MFE matrices ? */
  if (options & VRNA_OPTION_MFE)
    v |= (mx_type == VRNA_MX_WINDOW) ? ALLOC_MFE_LOCAL : ALLOC_MFE_DEFAULT;

  /* default PF matrices ? */
  if (options & VRNA_OPTION_PF)
    v |= (md_p->compute_bpp) ? ALLOC_PF_DEFAULT : ALLOC_PF_WO_PROBS;

  if ((fc->strands > 1) || (options & VRNA_OPTION_HYBRID))
    v |= ALLOC_MULTISTRAND;

  /* matrices for circular folding ? */
  if (md_p->circ)
    v |= ALLOC_CIRC;

  /* unique ML decomposition ? */
  if (md_p->uniq_ML)
    v |= ALLOC_UNIQ;

  return v;
}


PRIVATE void
mfe_matrices_free_default(vrna_mx_mfe_t *self)
{
  free(self->f5);
  free(self->f3);

  if (self->fms5)
    for (unsigned int s = 0; s < self->strands; s++)
      free(self->fms5[s]);

  free(self->fms5);

  if (self->fms3)
    for (unsigned int s = 0; s < self->strands; s++)
      free(self->fms3[s]);

  free(self->fms3);

  free(self->c);
  free(self->fML);
  free(self->fM1);
  free(self->fM2);
  free(self->fM1_new);
  free(self->fM2_real);
#ifndef VRNA_DISABLE_C11_FEATURES
  vrna_smx_csr_free(self->c_gq);
#else
  vrna_smx_csr_int_free(self->c_gq);
#endif
}


PRIVATE void
mfe_matrices_free_window(vrna_mx_mfe_t  *self,
                         unsigned int   length VRNA_UNUSED,
                         unsigned int   window_size VRNA_UNUSED)
{
  free(self->c_local);
  free(self->fML_local);
  free(self->ggg_local);
  free(self->f3_local);
}


PRIVATE void
mfe_matrices_free_2Dfold(vrna_mx_mfe_t  *self,
                         unsigned int   length,
                         int            turn,
                         int            *indx)
{
  unsigned int  i, j, ij;
  int           cnt1;

  /* This will be some fun... */
#ifdef COUNT_STATES
  if (self->N_F5 != NULL) {
    for (i = 1; i <= length; i++) {
      if (!self->N_F5[i])
        continue;

      for (cnt1 = self->k_min_F5[i]; cnt1 <= vars->k_max_F5[i]; cnt1++)
        if (vars->l_min_F5[i][cnt1] < INF) {
          vars->N_F5[i][cnt1] += vars->l_min_F5[i][cnt1] / 2;
          free(vars->N_F5[i][cnt1]);
        }

      if (vars->k_min_F5[i] < INF) {
        vars->N_F5[i] += vars->k_min_F5[i];
        free(vars->N_F5[i]);
      }
    }
    free(vars->N_F5);
  }

#endif

  if (self->E_F5 != NULL) {
    for (i = 1; i <= length; i++) {
      if (!self->E_F5[i])
        continue;

      for (cnt1 = self->k_min_F5[i]; cnt1 <= self->k_max_F5[i]; cnt1++)
        if (self->l_min_F5[i][cnt1] < INF) {
          self->E_F5[i][cnt1] += self->l_min_F5[i][cnt1] / 2;
          free(self->E_F5[i][cnt1]);
        }

      if (self->k_min_F5[i] < INF) {
        self->E_F5[i] += self->k_min_F5[i];
        free(self->E_F5[i]);
        self->l_min_F5[i] += self->k_min_F5[i];
        self->l_max_F5[i] += self->k_min_F5[i];
        free(self->l_min_F5[i]);
        free(self->l_max_F5[i]);
      }
    }
    free(self->E_F5);
    free(self->l_min_F5);
    free(self->l_max_F5);
    free(self->k_min_F5);
    free(self->k_max_F5);
  }

  if (self->E_F3 != NULL) {
    for (i = 1; i <= length; i++) {
      if (!self->E_F3[i])
        continue;

      for (cnt1 = self->k_min_F3[i]; cnt1 <= self->k_max_F3[i]; cnt1++)
        if (self->l_min_F3[i][cnt1] < INF) {
          self->E_F3[i][cnt1] += self->l_min_F3[i][cnt1] / 2;
          free(self->E_F3[i][cnt1]);
        }

      if (self->k_min_F3[i] < INF) {
        self->E_F3[i] += self->k_min_F3[i];
        free(self->E_F3[i]);
        self->l_min_F3[i] += self->k_min_F3[i];
        self->l_max_F3[i] += self->k_min_F3[i];
        free(self->l_min_F3[i]);
        free(self->l_max_F3[i]);
      }
    }
    free(self->E_F3);
    free(self->l_min_F3);
    free(self->l_max_F3);
    free(self->k_min_F3);
    free(self->k_max_F3);
  }

#ifdef COUNT_STATES
  if (self->N_C != NULL) {
    for (i = 1; i < length; i++) {
      for (j = i; j <= length; j++) {
        ij = indx[i] - j;
        if (!self->N_C[ij])
          continue;

        for (cnt1 = self->k_min_C[ij]; cnt1 <= self->k_max_C[ij]; cnt1++)
          if (self->l_min_C[ij][cnt1] < INF) {
            self->N_C[ij][cnt1] += self->l_min_C[ij][cnt1] / 2;
            free(self->N_C[ij][cnt1]);
          }

        if (self->k_min_C[ij] < INF) {
          self->N_C[ij] += self->k_min_C[ij];
          free(self->N_C[ij]);
        }
      }
    }
    free(self->N_C);
  }

#endif

  if (self->E_C != NULL) {
    for (i = 1; i < length; i++) {
      for (j = i; j <= length; j++) {
        ij = indx[i] - j;
        if (!self->E_C[ij])
          continue;

        for (cnt1 = self->k_min_C[ij]; cnt1 <= self->k_max_C[ij]; cnt1++)
          if (self->l_min_C[ij][cnt1] < INF) {
            self->E_C[ij][cnt1] += self->l_min_C[ij][cnt1] / 2;
            free(self->E_C[ij][cnt1]);
          }

        if (self->k_min_C[ij] < INF) {
          self->E_C[ij] += self->k_min_C[ij];
          free(self->E_C[ij]);
          self->l_min_C[ij] += self->k_min_C[ij];
          self->l_max_C[ij] += self->k_min_C[ij];
          free(self->l_min_C[ij]);
          free(self->l_max_C[ij]);
        }
      }
    }
    free(self->E_C);
    free(self->l_min_C);
    free(self->l_max_C);
    free(self->k_min_C);
    free(self->k_max_C);
  }

#ifdef COUNT_STATES
  if (self->N_M != NULL) {
    for (i = 1; i < length; i++) {
      for (j = i; j <= length; j++) {
        ij = indx[i] - j;
        if (!self->N_M[ij])
          continue;

        for (cnt1 = self->k_min_M[ij]; cnt1 <= self->k_max_M[ij]; cnt1++)
          if (self->l_min_M[ij][cnt1] < INF) {
            self->N_M[ij][cnt1] += self->l_min_M[ij][cnt1] / 2;
            free(self->N_M[ij][cnt1]);
          }

        if (self->k_min_M[ij] < INF) {
          self->N_M[ij] += self->k_min_M[ij];
          free(self->N_M[ij]);
        }
      }
    }
    free(self->N_M);
  }

#endif

  if (self->E_M != NULL) {
    for (i = 1; i < length; i++) {
      for (j = i; j <= length; j++) {
        ij = indx[i] - j;
        if (!self->E_M[ij])
          continue;

        for (cnt1 = self->k_min_M[ij]; cnt1 <= self->k_max_M[ij]; cnt1++)
          if (self->l_min_M[ij][cnt1] < INF) {
            self->E_M[ij][cnt1] += self->l_min_M[ij][cnt1] / 2;
            free(self->E_M[ij][cnt1]);
          }

        if (self->k_min_M[ij] < INF) {
          self->E_M[ij] += self->k_min_M[ij];
          free(self->E_M[ij]);
          self->l_min_M[ij] += self->k_min_M[ij];
          self->l_max_M[ij] += self->k_min_M[ij];
          free(self->l_min_M[ij]);
          free(self->l_max_M[ij]);
        }
      }
    }
    free(self->E_M);
    free(self->l_min_M);
    free(self->l_max_M);
    free(self->k_min_M);
    free(self->k_max_M);
  }

#ifdef COUNT_STATES
  if (self->N_M1 != NULL) {
    for (i = 1; i < length; i++) {
      for (j = i; j <= length; j++) {
        ij = indx[i] - j;
        if (!self->N_M1[ij])
          continue;

        for (cnt1 = self->k_min_M1[ij]; cnt1 <= self->k_max_M1[ij]; cnt1++)
          if (self->l_min_M1[ij][cnt1] < INF) {
            self->N_M1[ij][cnt1] += self->l_min_M1[ij][cnt1] / 2;
            free(self->N_M1[ij][cnt1]);
          }

        if (self->k_min_M1[ij] < INF) {
          self->N_M1[ij] += self->k_min_M1[ij];
          free(self->N_M1[ij]);
        }
      }
    }
    free(self->N_M1);
  }

#endif

  if (self->E_M1 != NULL) {
    for (i = 1; i < length; i++) {
      for (j = i; j <= length; j++) {
        ij = indx[i] - j;
        if (!self->E_M1[ij])
          continue;

        for (cnt1 = self->k_min_M1[ij]; cnt1 <= self->k_max_M1[ij]; cnt1++)
          if (self->l_min_M1[ij][cnt1] < INF) {
            self->E_M1[ij][cnt1] += self->l_min_M1[ij][cnt1] / 2;
            free(self->E_M1[ij][cnt1]);
          }

        if (self->k_min_M1[ij] < INF) {
          self->E_M1[ij] += self->k_min_M1[ij];
          free(self->E_M1[ij]);
          self->l_min_M1[ij]  += self->k_min_M1[ij];
          self->l_max_M1[ij]  += self->k_min_M1[ij];
          free(self->l_min_M1[ij]);
          free(self->l_max_M1[ij]);
        }
      }
    }
    free(self->E_M1);
    free(self->l_min_M1);
    free(self->l_max_M1);
    free(self->k_min_M1);
    free(self->k_max_M1);
  }

  if (self->E_M2 != NULL) {
    for (i = 1; i < length - turn - 1; i++) {
      if (!self->E_M2[i])
        continue;

      for (cnt1 = self->k_min_M2[i]; cnt1 <= self->k_max_M2[i]; cnt1++)
        if (self->l_min_M2[i][cnt1] < INF) {
          self->E_M2[i][cnt1] += self->l_min_M2[i][cnt1] / 2;
          free(self->E_M2[i][cnt1]);
        }

      if (self->k_min_M2[i] < INF) {
        self->E_M2[i] += self->k_min_M2[i];
        free(self->E_M2[i]);
        self->l_min_M2[i] += self->k_min_M2[i];
        self->l_max_M2[i] += self->k_min_M2[i];
        free(self->l_min_M2[i]);
        free(self->l_max_M2[i]);
      }
    }
    free(self->E_M2);
    free(self->l_min_M2);
    free(self->l_max_M2);
    free(self->k_min_M2);
    free(self->k_max_M2);
  }

  if (self->E_Fc != NULL) {
    for (cnt1 = self->k_min_Fc; cnt1 <= self->k_max_Fc; cnt1++)
      if (self->l_min_Fc[cnt1] < INF) {
        self->E_Fc[cnt1] += self->l_min_Fc[cnt1] / 2;
        free(self->E_Fc[cnt1]);
      }

    if (self->k_min_Fc < INF) {
      self->E_Fc += self->k_min_Fc;
      free(self->E_Fc);
      self->l_min_Fc  += self->k_min_Fc;
      self->l_max_Fc  += self->k_min_Fc;
      free(self->l_min_Fc);
      free(self->l_max_Fc);
    }
  }

  if (self->E_FcI != NULL) {
    for (cnt1 = self->k_min_FcI; cnt1 <= self->k_max_FcI; cnt1++)
      if (self->l_min_FcI[cnt1] < INF) {
        self->E_FcI[cnt1] += self->l_min_FcI[cnt1] / 2;
        free(self->E_FcI[cnt1]);
      }

    if (self->k_min_FcI < INF) {
      self->E_FcI += self->k_min_FcI;
      free(self->E_FcI);
      self->l_min_FcI += self->k_min_FcI;
      self->l_max_FcI += self->k_min_FcI;
      free(self->l_min_FcI);
      free(self->l_max_FcI);
    }
  }

  if (self->E_FcH != NULL) {
    for (cnt1 = self->k_min_FcH; cnt1 <= self->k_max_FcH; cnt1++)
      if (self->l_min_FcH[cnt1] < INF) {
        self->E_FcH[cnt1] += self->l_min_FcH[cnt1] / 2;
        free(self->E_FcH[cnt1]);
      }

    if (self->k_min_FcH < INF) {
      self->E_FcH += self->k_min_FcH;
      free(self->E_FcH);
      self->l_min_FcH += self->k_min_FcH;
      self->l_max_FcH += self->k_min_FcH;
      free(self->l_min_FcH);
      free(self->l_max_FcH);
    }
  }

  if (self->E_FcM != NULL) {
    for (cnt1 = self->k_min_FcM; cnt1 <= self->k_max_FcM; cnt1++)
      if (self->l_min_FcM[cnt1] < INF) {
        self->E_FcM[cnt1] += self->l_min_FcM[cnt1] / 2;
        free(self->E_FcM[cnt1]);
      }

    if (self->k_min_FcM < INF) {
      self->E_FcM += self->k_min_FcM;
      free(self->E_FcM);
      self->l_min_FcM += self->k_min_FcM;
      self->l_max_FcM += self->k_min_FcM;
      free(self->l_min_FcM);
      free(self->l_max_FcM);
    }
  }

  free(self->E_F5_rem);
  free(self->E_F3_rem);
  free(self->E_C_rem);
  free(self->E_M_rem);
  free(self->E_M1_rem);
  free(self->E_M2_rem);
}


PRIVATE void
pf_matrices_free_default(vrna_mx_pf_t *self)
{
  free(self->q);
  free(self->qb);
  free(self->qm);
  free(self->qm1);
  free(self->qm2);
  free(self->qm2_real);
  free(self->qm1_new);
  free(self->probs);
  free(self->q1k);
  free(self->qln);
#ifndef VRNA_DISABLE_C11_FEATURES
  vrna_smx_csr_free(self->q_gq);
  vrna_smx_csr_free(self->p_gq);
#else
  vrna_smx_csr_FLT_OR_DBL_free(self->q_gq);
  vrna_smx_csr_FLT_OR_DBL_free(self->p_gq);
#endif
}


PRIVATE void
pf_matrices_free_window(vrna_mx_pf_t  *self,
                        unsigned int  length VRNA_UNUSED,
                        unsigned int  window_size VRNA_UNUSED)
{
  free(self->q_local);
  free(self->qb_local);
  free(self->qm_local);
  free(self->qm2_local);
  free(self->pR);
  free(self->QI5);
  free(self->q2l);
  free(self->qmb);
  free(self->G_local);
}


PRIVATE void
pf_matrices_free_2Dfold(vrna_mx_pf_t  *self,
                        unsigned int  length,
                        int           turn,
                        int           *indx,
                        int           *jindx)
{
  unsigned int  i, j, ij;
  int           cnt1;

  /* This will be some fun... */
  if (self->Q != NULL) {
    for (i = 1; i <= length; i++) {
      for (j = i; j <= length; j++) {
        ij = indx[i] - j;
        if (!self->Q[ij])
          continue;

        for (cnt1 = self->k_min_Q[ij]; cnt1 <= self->k_max_Q[ij]; cnt1++)
          if (self->l_min_Q[ij][cnt1] < INF) {
            self->Q[ij][cnt1] += self->l_min_Q[ij][cnt1] / 2;
            free(self->Q[ij][cnt1]);
          }

        if (self->k_min_Q[ij] < INF) {
          self->Q[ij] += self->k_min_Q[ij];
          free(self->Q[ij]);
          self->l_min_Q[ij] += self->k_min_Q[ij];
          self->l_max_Q[ij] += self->k_min_Q[ij];
          free(self->l_min_Q[ij]);
          free(self->l_max_Q[ij]);
        }
      }
    }
  }

  free(self->Q);
  free(self->l_min_Q);
  free(self->l_max_Q);
  free(self->k_min_Q);
  free(self->k_max_Q);

  if (self->Q_B != NULL) {
    for (i = 1; i < length; i++) {
      for (j = i; j <= length; j++) {
        ij = indx[i] - j;
        if (!self->Q_B[ij])
          continue;

        for (cnt1 = self->k_min_Q_B[ij]; cnt1 <= self->k_max_Q_B[ij]; cnt1++)
          if (self->l_min_Q_B[ij][cnt1] < INF) {
            self->Q_B[ij][cnt1] += self->l_min_Q_B[ij][cnt1] / 2;
            free(self->Q_B[ij][cnt1]);
          }

        if (self->k_min_Q_B[ij] < INF) {
          self->Q_B[ij] += self->k_min_Q_B[ij];
          free(self->Q_B[ij]);
          self->l_min_Q_B[ij] += self->k_min_Q_B[ij];
          self->l_max_Q_B[ij] += self->k_min_Q_B[ij];
          free(self->l_min_Q_B[ij]);
          free(self->l_max_Q_B[ij]);
        }
      }
    }
  }

  free(self->Q_B);
  free(self->l_min_Q_B);
  free(self->l_max_Q_B);
  free(self->k_min_Q_B);
  free(self->k_max_Q_B);

  if (self->Q_M != NULL) {
    for (i = 1; i < length; i++) {
      for (j = i; j <= length; j++) {
        ij = indx[i] - j;
        if (!self->Q_M[ij])
          continue;

        for (cnt1 = self->k_min_Q_M[ij]; cnt1 <= self->k_max_Q_M[ij]; cnt1++)
          if (self->l_min_Q_M[ij][cnt1] < INF) {
            self->Q_M[ij][cnt1] += self->l_min_Q_M[ij][cnt1] / 2;
            free(self->Q_M[ij][cnt1]);
          }

        if (self->k_min_Q_M[ij] < INF) {
          self->Q_M[ij] += self->k_min_Q_M[ij];
          free(self->Q_M[ij]);
          self->l_min_Q_M[ij] += self->k_min_Q_M[ij];
          self->l_max_Q_M[ij] += self->k_min_Q_M[ij];
          free(self->l_min_Q_M[ij]);
          free(self->l_max_Q_M[ij]);
        }
      }
    }
  }

  free(self->Q_M);
  free(self->l_min_Q_M);
  free(self->l_max_Q_M);
  free(self->k_min_Q_M);
  free(self->k_max_Q_M);

  if (self->Q_M1 != NULL) {
    for (i = 1; i < length; i++) {
      for (j = i; j <= length; j++) {
        ij = jindx[j] + i;
        if (!self->Q_M1[ij])
          continue;

        for (cnt1 = self->k_min_Q_M1[ij]; cnt1 <= self->k_max_Q_M1[ij]; cnt1++)
          if (self->l_min_Q_M1[ij][cnt1] < INF) {
            self->Q_M1[ij][cnt1] += self->l_min_Q_M1[ij][cnt1] / 2;
            free(self->Q_M1[ij][cnt1]);
          }

        if (self->k_min_Q_M1[ij] < INF) {
          self->Q_M1[ij] += self->k_min_Q_M1[ij];
          free(self->Q_M1[ij]);
          self->l_min_Q_M1[ij]  += self->k_min_Q_M1[ij];
          self->l_max_Q_M1[ij]  += self->k_min_Q_M1[ij];
          free(self->l_min_Q_M1[ij]);
          free(self->l_max_Q_M1[ij]);
        }
      }
    }
  }

  free(self->Q_M1);
  free(self->l_min_Q_M1);
  free(self->l_max_Q_M1);
  free(self->k_min_Q_M1);
  free(self->k_max_Q_M1);

  if (self->Q_M2 != NULL) {
    for (i = 1; i < length - turn - 1; i++) {
      if (!self->Q_M2[i])
        continue;

      for (cnt1 = self->k_min_Q_M2[i]; cnt1 <= self->k_max_Q_M2[i]; cnt1++)
        if (self->l_min_Q_M2[i][cnt1] < INF) {
          self->Q_M2[i][cnt1] += self->l_min_Q_M2[i][cnt1] / 2;
          free(self->Q_M2[i][cnt1]);
        }

      if (self->k_min_Q_M2[i] < INF) {
        self->Q_M2[i] += self->k_min_Q_M2[i];
        free(self->Q_M2[i]);
        self->l_min_Q_M2[i] += self->k_min_Q_M2[i];
        self->l_max_Q_M2[i] += self->k_min_Q_M2[i];
        free(self->l_min_Q_M2[i]);
        free(self->l_max_Q_M2[i]);
      }
    }
  }

  free(self->Q_M2);
  free(self->l_min_Q_M2);
  free(self->l_max_Q_M2);
  free(self->k_min_Q_M2);
  free(self->k_max_Q_M2);

  if (self->Q_c != NULL) {
    for (cnt1 = self->k_min_Q_c; cnt1 <= self->k_max_Q_c; cnt1++)
      if (self->l_min_Q_c[cnt1] < INF) {
        self->Q_c[cnt1] += self->l_min_Q_c[cnt1] / 2;
        free(self->Q_c[cnt1]);
      }

    if (self->k_min_Q_c < INF) {
      self->Q_c += self->k_min_Q_c;
      free(self->Q_c);
      self->l_min_Q_c += self->k_min_Q_c;
      self->l_max_Q_c += self->k_min_Q_c;
      free(self->l_min_Q_c);
      free(self->l_max_Q_c);
    }
  }

  if (self->Q_cI != NULL) {
    for (cnt1 = self->k_min_Q_cI; cnt1 <= self->k_max_Q_cI; cnt1++)
      if (self->l_min_Q_cI[cnt1] < INF) {
        self->Q_cI[cnt1] += self->l_min_Q_cI[cnt1] / 2;
        free(self->Q_cI[cnt1]);
      }

    if (self->k_min_Q_cI < INF) {
      self->Q_cI += self->k_min_Q_cI;
      free(self->Q_cI);
      self->l_min_Q_cI  += self->k_min_Q_cI;
      self->l_max_Q_cI  += self->k_min_Q_cI;
      free(self->l_min_Q_cI);
      free(self->l_max_Q_cI);
    }
  }

  if (self->Q_cH != NULL) {
    for (cnt1 = self->k_min_Q_cH; cnt1 <= self->k_max_Q_cH; cnt1++)
      if (self->l_min_Q_cH[cnt1] < INF) {
        self->Q_cH[cnt1] += self->l_min_Q_cH[cnt1] / 2;
        free(self->Q_cH[cnt1]);
      }

    if (self->k_min_Q_cH < INF) {
      self->Q_cH += self->k_min_Q_cH;
      free(self->Q_cH);
      self->l_min_Q_cH  += self->k_min_Q_cH;
      self->l_max_Q_cH  += self->k_min_Q_cH;
      free(self->l_min_Q_cH);
      free(self->l_max_Q_cH);
    }
  }

  if (self->Q_cM != NULL) {
    for (cnt1 = self->k_min_Q_cM; cnt1 <= self->k_max_Q_cM; cnt1++)
      if (self->l_min_Q_cM[cnt1] < INF) {
        self->Q_cM[cnt1] += self->l_min_Q_cM[cnt1] / 2;
        free(self->Q_cM[cnt1]);
      }

    if (self->k_min_Q_cM < INF) {
      self->Q_cM += self->k_min_Q_cM;
      free(self->Q_cM);
      self->l_min_Q_cM  += self->k_min_Q_cM;
      self->l_max_Q_cM  += self->k_min_Q_cM;
      free(self->l_min_Q_cM);
      free(self->l_max_Q_cM);
    }
  }

  free(self->Q_rem);
  free(self->Q_B_rem);
  free(self->Q_M_rem);
  free(self->Q_M1_rem);
  free(self->Q_M2_rem);
}


PRIVATE vrna_mx_mfe_t *
init_mx_mfe_default(vrna_fold_compound_t  *fc,
                    unsigned int          alloc_vector)
{
  unsigned int  n, size, lin_size, s, strands;
  vrna_mx_mfe_t *mx;
  vrna_mx_mfe_t init = {
    .type = VRNA_MX_DEFAULT
  };

  n = fc->length;

  if ((int)(n * n) >= (int)INT_MAX) {
    vrna_log_warning("init_mx_mfe_default(): "
                         "sequence length %d exceeds addressable range",
                         n);
    return NULL;
  }

  mx = vrna_alloc(sizeof(vrna_mx_mfe_t));

  if (mx) {
    memcpy(mx, &init, sizeof(vrna_mx_mfe_t));
    nullify_mfe(mx);

    strands     = fc->strands;
    size        = ((n + 1) * (n + 2)) / 2;
    lin_size    = n + 2;
    mx->length  = n;
    mx->strands = strands;

    if (alloc_vector & ALLOC_F5)
      mx->f5 = (int *)vrna_alloc(sizeof(int) * lin_size);

    if (alloc_vector & ALLOC_F3)
      mx->f3 = (int *)vrna_alloc(sizeof(int) * lin_size);

    if (alloc_vector & ALLOC_MULTISTRAND) {
      mx->fms5  = (int **)vrna_alloc(sizeof(int *) * strands);
      mx->fms3  = (int **)vrna_alloc(sizeof(int *) * strands);

      for (s = 0; s < strands; s++) {
        mx->fms5[s] = (int *)vrna_alloc(sizeof(int) * (n + 1));
        mx->fms3[s] = (int *)vrna_alloc(sizeof(int) * (n + 1));
      }
    }

    if (alloc_vector & ALLOC_C)
      mx->c = (int *)vrna_alloc(sizeof(int) * size);

    if (alloc_vector & ALLOC_FML)
      mx->fML = (int *)vrna_alloc(sizeof(int) * size);

    if ((alloc_vector & ALLOC_UNIQ) ||
        (alloc_vector & ALLOC_CIRC)) {
      if (!(alloc_vector & ALLOC_CIRC))
        mx->fM1 = (int *)vrna_alloc(sizeof(int) * size);

      mx->fM2_real = (int *)vrna_alloc(sizeof(int) * size);
    }

    if (alloc_vector & ALLOC_CIRC)
      mx->fM1_new  = (int *)vrna_alloc(sizeof(int) * lin_size);
  }

  return mx;
}


PRIVATE vrna_mx_mfe_t *
init_mx_mfe_window(vrna_fold_compound_t *fc,
                   unsigned int         alloc_vector)
{
  unsigned int  n, m, lin_size;
  vrna_mx_mfe_t *mx;
  vrna_mx_mfe_t init = {
    .type = VRNA_MX_WINDOW
  };

  n = fc->length;
  m = fc->window_size;

  if ((int)(n * m) >= (int)INT_MAX) {
    vrna_log_warning("init_mx_mfe_window(): "
                         "sequence length %d exceeds addressable range",
                         n);
    return NULL;
  }

  mx = vrna_alloc(sizeof(vrna_mx_mfe_t));

  if (mx) {
    memcpy(mx, &init, sizeof(vrna_mx_mfe_t));
    nullify_mfe(mx);

    lin_size    = n + 2;
    mx->length  = n;
    mx->strands = fc->strands;

    if (alloc_vector & ALLOC_F3)
      mx->f3_local = (int *)vrna_alloc(sizeof(int) * lin_size);

    if (alloc_vector & ALLOC_C)
      mx->c_local = (int **)vrna_alloc(sizeof(int *) * lin_size);

    if (alloc_vector & ALLOC_FML)
      mx->fML_local = (int **)vrna_alloc(sizeof(int *) * lin_size);
  }

  return mx;
}


PRIVATE vrna_mx_mfe_t *
init_mx_mfe_2Dfold(vrna_fold_compound_t *fc,
                   unsigned int         alloc_vector)
{
  unsigned int  n, i, size, lin_size;
  vrna_mx_mfe_t *mx;
  vrna_mx_mfe_t init = {
    .type = VRNA_MX_2DFOLD
  };

  n = fc->length;

  if ((int)(n * n) >= (int)INT_MAX) {
    vrna_log_warning("init_mx_mfe_2Dfold(): "
                         "sequence length %d exceeds addressable range",
                         n);
    return NULL;
  }

  mx = vrna_alloc(sizeof(vrna_mx_mfe_t));

  if (mx) {
    memcpy(mx, &init, sizeof(vrna_mx_mfe_t));
    nullify_mfe(mx);

    size        = ((n + 1) * (n + 2)) / 2;
    lin_size    = n + 2;
    mx->length  = n;
    mx->strands = fc->strands;

    if (alloc_vector & ALLOC_F5) {
      mx->E_F5      = (int ***)vrna_alloc(sizeof(int **) * lin_size);
      mx->l_min_F5  = (int **)vrna_alloc(sizeof(int *) * lin_size);
      mx->l_max_F5  = (int **)vrna_alloc(sizeof(int *) * lin_size);
      mx->k_min_F5  = (int *)vrna_alloc(sizeof(int) * lin_size);
      mx->k_max_F5  = (int *)vrna_alloc(sizeof(int) * lin_size);
      mx->E_F5_rem  = (int *)vrna_alloc(sizeof(int) * lin_size);
      for (i = 0; i <= n; i++)
        mx->E_F5_rem[i] = INF;
    }

    if (alloc_vector & ALLOC_F3) {
      mx->E_F3      = (int ***)vrna_alloc(sizeof(int **) * lin_size);
      mx->l_min_F3  = (int **)vrna_alloc(sizeof(int *) * lin_size);
      mx->l_max_F3  = (int **)vrna_alloc(sizeof(int *) * lin_size);
      mx->k_min_F3  = (int *)vrna_alloc(sizeof(int) * lin_size);
      mx->k_max_F3  = (int *)vrna_alloc(sizeof(int) * lin_size);
      mx->E_F3_rem  = (int *)vrna_alloc(sizeof(int) * lin_size);
      for (i = 0; i <= n; i++)
        mx->E_F3_rem[i] = INF;
    }

    if (alloc_vector & ALLOC_C) {
      mx->E_C     = (int ***)vrna_alloc(sizeof(int **) * size);
      mx->l_min_C = (int **)vrna_alloc(sizeof(int *) * size);
      mx->l_max_C = (int **)vrna_alloc(sizeof(int *) * size);
      mx->k_min_C = (int *)vrna_alloc(sizeof(int) * size);
      mx->k_max_C = (int *)vrna_alloc(sizeof(int) * size);
      mx->E_C_rem = (int *)vrna_alloc(sizeof(int) * size);
      for (i = 0; i < size; i++)
        mx->E_C_rem[i] = INF;
    }

    if (alloc_vector & ALLOC_FML) {
      mx->E_M     = (int ***)vrna_alloc(sizeof(int **) * size);
      mx->l_min_M = (int **)vrna_alloc(sizeof(int *) * size);
      mx->l_max_M = (int **)vrna_alloc(sizeof(int *) * size);
      mx->k_min_M = (int *)vrna_alloc(sizeof(int) * size);
      mx->k_max_M = (int *)vrna_alloc(sizeof(int) * size);
      mx->E_M_rem = (int *)vrna_alloc(sizeof(int) * size);
      for (i = 0; i < size; i++)
        mx->E_M_rem[i] = INF;
    }

    if (alloc_vector & ALLOC_UNIQ) {
      mx->E_M1      = (int ***)vrna_alloc(sizeof(int **) * size);
      mx->l_min_M1  = (int **)vrna_alloc(sizeof(int *) * size);
      mx->l_max_M1  = (int **)vrna_alloc(sizeof(int *) * size);
      mx->k_min_M1  = (int *)vrna_alloc(sizeof(int) * size);
      mx->k_max_M1  = (int *)vrna_alloc(sizeof(int) * size);
      mx->E_M1_rem  = (int *)vrna_alloc(sizeof(int) * size);
      for (i = 0; i < size; i++)
        mx->E_M1_rem[i] = INF;
    }

    if (alloc_vector & ALLOC_CIRC) {
      mx->E_M2      = (int ***)vrna_alloc(sizeof(int **) * lin_size);
      mx->l_min_M2  = (int **)vrna_alloc(sizeof(int *) * lin_size);
      mx->l_max_M2  = (int **)vrna_alloc(sizeof(int *) * lin_size);
      mx->k_min_M2  = (int *)vrna_alloc(sizeof(int) * lin_size);
      mx->k_max_M2  = (int *)vrna_alloc(sizeof(int) * lin_size);
      mx->E_M2_rem  = (int *)vrna_alloc(sizeof(int) * lin_size);
      for (i = 0; i <= n; i++)
        mx->E_M2_rem[i] = INF;
    }

#ifdef COUNT_STATES
    mx->N_C   = (unsigned long ***)vrna_alloc(sizeof(unsigned long **) * size);
    mx->N_F5  = (unsigned long ***)vrna_alloc(sizeof(unsigned long **) * lin_size);
    mx->N_M   = (unsigned long ***)vrna_alloc(sizeof(unsigned long **) * size);
    mx->N_M1  = (unsigned long ***)vrna_alloc(sizeof(unsigned long **) * size);
#endif
  }

  return mx;
}


PRIVATE INLINE void
nullify_mfe(vrna_mx_mfe_t *mx)
{
  if (mx) {
    mx->length  = 0;
    mx->strands = 0;

    switch (mx->type) {
      case VRNA_MX_DEFAULT:
        mx->c     = NULL;
        mx->f5    = NULL;
        mx->f3    = NULL;
        mx->fms5  = NULL;
        mx->fms3  = NULL;
        mx->fML   = NULL;
        mx->fM1   = NULL;
        mx->fM2   = NULL;
        mx->fM1_new   = NULL;
        mx->fM2_real  = NULL;
        mx->Fc    = INF;
        mx->FcH   = INF;
        mx->FcI   = INF;
        mx->FcM   = INF;
        mx->c_gq  = NULL;
        break;

      case VRNA_MX_WINDOW:
        mx->c_local   = NULL;
        mx->f3_local  = NULL;
        mx->fML_local = NULL;
        mx->ggg_local = NULL;
        break;

      case VRNA_MX_2DFOLD:
        mx->E_F5      = NULL;
        mx->l_min_F5  = NULL;
        mx->l_max_F5  = NULL;
        mx->k_min_F5  = NULL;
        mx->k_max_F5  = NULL;
        mx->E_F5_rem  = NULL;

        mx->E_F3      = NULL;
        mx->l_min_F3  = NULL;
        mx->l_max_F3  = NULL;
        mx->k_min_F3  = NULL;
        mx->k_max_F3  = NULL;
        mx->E_F3_rem  = NULL;

        mx->E_C     = NULL;
        mx->l_min_C = NULL;
        mx->l_max_C = NULL;
        mx->k_min_C = NULL;
        mx->k_max_C = NULL;
        mx->E_C_rem = NULL;

        mx->E_M     = NULL;
        mx->l_min_M = NULL;
        mx->l_max_M = NULL;
        mx->k_min_M = NULL;
        mx->k_max_M = NULL;
        mx->E_M_rem = NULL;

        mx->E_M1      = NULL;
        mx->l_min_M1  = NULL;
        mx->l_max_M1  = NULL;
        mx->k_min_M1  = NULL;
        mx->k_max_M1  = NULL;
        mx->E_M1_rem  = NULL;

        mx->E_M2      = NULL;
        mx->l_min_M2  = NULL;
        mx->l_max_M2  = NULL;
        mx->k_min_M2  = NULL;
        mx->k_max_M2  = NULL;
        mx->E_M2_rem  = NULL;

        mx->E_Fc      = NULL;
        mx->l_min_Fc  = NULL;
        mx->l_max_Fc  = NULL;
        mx->k_min_Fc  = 0;
        mx->k_max_Fc  = 0;
        mx->E_Fc_rem  = INF;

        mx->E_FcH     = NULL;
        mx->l_min_FcH = NULL;
        mx->l_max_FcH = NULL;
        mx->k_min_FcH = 0;
        mx->k_max_FcH = 0;
        mx->E_FcH_rem = INF;

        mx->E_FcI     = NULL;
        mx->l_min_FcI = NULL;
        mx->l_max_FcI = NULL;
        mx->k_min_FcI = 0;
        mx->k_max_FcI = 0;
        mx->E_FcI_rem = INF;

        mx->E_FcM     = NULL;
        mx->l_min_FcM = NULL;
        mx->l_max_FcM = NULL;
        mx->k_min_FcM = 0;
        mx->k_max_FcM = 0;
        mx->E_FcM_rem = INF;

#ifdef COUNT_STATES
        mx->N_F5  = NULL;
        mx->N_C   = NULL;
        mx->N_M   = NULL;
        mx->N_M1  = NULL;
#endif

        break;
    }
  }
}


PRIVATE vrna_mx_pf_t *
init_mx_pf_default(vrna_fold_compound_t *fc,
                   unsigned int         alloc_vector)
{
  unsigned int  n, size, lin_size;
  vrna_mx_pf_t  *mx;
  vrna_mx_pf_t  init = {
    .type = VRNA_MX_DEFAULT
  };

  n = fc->length;

  if ((int)(n * n) >= (int)INT_MAX) {
    vrna_log_warning("init_mx_pf_default(): "
                         "sequence length %d exceeds addressable range",
                         n);
    return NULL;
  }

  mx = vrna_alloc(sizeof(vrna_mx_pf_t));

  if (mx) {
    memcpy(mx, &init, sizeof(vrna_mx_pf_t));
    nullify_pf(mx);

    size        = ((n + 1) * (n + 2)) / 2;
    lin_size    = n + 2;
    mx->length  = n;

    if (alloc_vector & ALLOC_F)
      mx->q = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * size);

    if (alloc_vector & ALLOC_C)
      mx->qb = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * size);

    if (alloc_vector & ALLOC_FML)
      mx->qm = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * size);

    if ((alloc_vector & ALLOC_UNIQ) ||
        (alloc_vector & ALLOC_CIRC)) {
      if (!(alloc_vector & ALLOC_CIRC))
        mx->qm1 = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * size);

      mx->qm2_real = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * size);
    }

    if (alloc_vector & ALLOC_CIRC)
      mx->qm1_new   = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * lin_size);

    if (alloc_vector & ALLOC_PROBS)
      mx->probs = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * size);

    if (alloc_vector & ALLOC_AUX) {
      mx->q1k = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * lin_size);
      mx->qln = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * lin_size);
    }

    /*
     *  always alloc the helper arrays for unpaired nucleotides in multi-
     *  branch loops and scaling
     */
    mx->scale     = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * lin_size);
    mx->expMLbase = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * lin_size);
  }

  return mx;
}


PRIVATE vrna_mx_pf_t *
init_mx_pf_window(vrna_fold_compound_t  *fc,
                  unsigned int          alloc_vector)
{
  unsigned int  n, m, lin_size;
  vrna_mx_pf_t  *mx;
  vrna_mx_pf_t  init = {
    .type = VRNA_MX_WINDOW
  };

  n = fc->length;
  m = fc->window_size;

  if ((int)(n * m) >= (int)INT_MAX) {
    vrna_log_warning("init_mx_pf_window(): "
                         "sequence length %d exceeds addressable range",
                         n);
    return NULL;
  }

  mx = vrna_alloc(sizeof(vrna_mx_pf_t));

  if (mx) {
    memcpy(mx, &init, sizeof(vrna_mx_pf_t));
    nullify_pf(mx);

    lin_size    = n + 2;
    mx->length  = n;

    if (alloc_vector & ALLOC_F)
      mx->q_local = (FLT_OR_DBL **)vrna_alloc(sizeof(FLT_OR_DBL *) * lin_size);

    if (alloc_vector & ALLOC_C)
      mx->qb_local = (FLT_OR_DBL **)vrna_alloc(sizeof(FLT_OR_DBL *) * lin_size);

    if (alloc_vector & ALLOC_FML)
      mx->qm_local = (FLT_OR_DBL **)vrna_alloc(sizeof(FLT_OR_DBL *) * lin_size);

    mx->pR = (FLT_OR_DBL **)vrna_alloc(sizeof(FLT_OR_DBL *) * lin_size);

    if (alloc_vector & ALLOC_PROBS) {
      mx->QI5       = (FLT_OR_DBL **)vrna_alloc(sizeof(FLT_OR_DBL *) * lin_size);
      mx->qmb       = (FLT_OR_DBL **)vrna_alloc(sizeof(FLT_OR_DBL *) * lin_size);
      mx->qm2_local = (FLT_OR_DBL **)vrna_alloc(sizeof(FLT_OR_DBL *) * lin_size);
      mx->q2l       = (FLT_OR_DBL **)vrna_alloc(sizeof(FLT_OR_DBL *) * lin_size);
    }

    /*
     *  always alloc the helper arrays for unpaired nucleotides in multi-
     *  branch loops and scaling
     */
    mx->scale     = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * lin_size);
    mx->expMLbase = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * lin_size);
  }

  return mx;
}


PRIVATE vrna_mx_pf_t *
init_mx_pf_2Dfold(vrna_fold_compound_t  *fc,
                  unsigned int          alloc_vector)
{
  unsigned int  n, size, lin_size;
  vrna_mx_pf_t  *mx;
  vrna_mx_pf_t  init = {
    .type = VRNA_MX_2DFOLD
  };

  n = fc->length;

  if ((int)(n * n) >= (int)INT_MAX) {
    vrna_log_warning("init_mx_pf_2Dfold(): "
                         "sequence length %d exceeds addressable range",
                         n);
    return NULL;
  }

  mx = vrna_alloc(sizeof(vrna_mx_pf_t));

  if (mx) {
    memcpy(mx, &init, sizeof(vrna_mx_pf_t));
    nullify_pf(mx);

    size        = ((n + 1) * (n + 2)) / 2;
    lin_size    = n + 2;
    mx->length  = n;

    if (alloc_vector & ALLOC_F) {
      mx->Q       = (FLT_OR_DBL ***)vrna_alloc(sizeof(FLT_OR_DBL * *) * size);
      mx->l_min_Q = (int **)vrna_alloc(sizeof(int *) * size);
      mx->l_max_Q = (int **)vrna_alloc(sizeof(int *) * size);
      mx->k_min_Q = (int *)vrna_alloc(sizeof(int) * size);
      mx->k_max_Q = (int *)vrna_alloc(sizeof(int) * size);
      mx->Q_rem   = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * size);
    }

    if (alloc_vector & ALLOC_C) {
      mx->Q_B       = (FLT_OR_DBL ***)vrna_alloc(sizeof(FLT_OR_DBL * *) * size);
      mx->l_min_Q_B = (int **)vrna_alloc(sizeof(int *) * size);
      mx->l_max_Q_B = (int **)vrna_alloc(sizeof(int *) * size);
      mx->k_min_Q_B = (int *)vrna_alloc(sizeof(int) * size);
      mx->k_max_Q_B = (int *)vrna_alloc(sizeof(int) * size);
      mx->Q_B_rem   = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * size);
    }

    if (alloc_vector & ALLOC_FML) {
      mx->Q_M       = (FLT_OR_DBL ***)vrna_alloc(sizeof(FLT_OR_DBL * *) * size);
      mx->l_min_Q_M = (int **)vrna_alloc(sizeof(int *) * size);
      mx->l_max_Q_M = (int **)vrna_alloc(sizeof(int *) * size);
      mx->k_min_Q_M = (int *)vrna_alloc(sizeof(int) * size);
      mx->k_max_Q_M = (int *)vrna_alloc(sizeof(int) * size);
      mx->Q_M_rem   = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * size);
    }

    if (alloc_vector & ALLOC_UNIQ) {
      mx->Q_M1        = (FLT_OR_DBL ***)vrna_alloc(sizeof(FLT_OR_DBL * *) * size);
      mx->l_min_Q_M1  = (int **)vrna_alloc(sizeof(int *) * size);
      mx->l_max_Q_M1  = (int **)vrna_alloc(sizeof(int *) * size);
      mx->k_min_Q_M1  = (int *)vrna_alloc(sizeof(int) * size);
      mx->k_max_Q_M1  = (int *)vrna_alloc(sizeof(int) * size);
      mx->Q_M1_rem    = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * size);
    }

    if (alloc_vector & ALLOC_CIRC) {
      mx->Q_M2        = (FLT_OR_DBL ***)vrna_alloc(sizeof(FLT_OR_DBL * *) * lin_size);
      mx->l_min_Q_M2  = (int **)vrna_alloc(sizeof(int *) * lin_size);
      mx->l_max_Q_M2  = (int **)vrna_alloc(sizeof(int *) * lin_size);
      mx->k_min_Q_M2  = (int *)vrna_alloc(sizeof(int) * lin_size);
      mx->k_max_Q_M2  = (int *)vrna_alloc(sizeof(int) * lin_size);
      mx->Q_M2_rem    = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * lin_size);
    }

    /*
     *  always alloc the helper arrays for unpaired nucleotides in multi-
     *  branch loops and scaling
     */
    mx->scale     = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * lin_size);
    mx->expMLbase = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * lin_size);
  }

  return mx;
}


PRIVATE INLINE void
nullify_pf(vrna_mx_pf_t *mx)
{
  if (mx) {
    mx->length    = 0;
    mx->scale     = NULL;
    mx->expMLbase = NULL;

    switch (mx->type) {
      case VRNA_MX_DEFAULT:
        mx->q     = NULL;
        mx->qb    = NULL;
        mx->qm    = NULL;
        mx->qm1   = NULL;
        mx->qm2   = NULL;
        mx->qm2_real   = NULL;
        mx->qm1_new   = NULL;
        mx->probs = NULL;
        mx->q1k   = NULL;
        mx->qln   = NULL;
        mx->q_gq  = NULL;
        mx->p_gq  = NULL;
        break;

      case VRNA_MX_WINDOW:
        mx->q_local   = NULL;
        mx->qb_local  = NULL;
        mx->qm_local  = NULL;
        mx->qm2_local = NULL;
        mx->pR        = NULL;
        mx->QI5       = NULL;
        mx->q2l       = NULL;
        mx->qmb       = NULL;
        mx->G_local   = NULL;
        break;

      case VRNA_MX_2DFOLD:
        mx->Q       = NULL;
        mx->l_min_Q = NULL;
        mx->l_max_Q = NULL;
        mx->k_min_Q = NULL;
        mx->k_max_Q = NULL;
        mx->Q_rem   = NULL;

        mx->Q_B       = NULL;
        mx->l_min_Q_B = NULL;
        mx->l_max_Q_B = NULL;
        mx->k_min_Q_B = NULL;
        mx->k_max_Q_B = NULL;
        mx->Q_B_rem   = NULL;

        mx->Q_M       = NULL;
        mx->l_min_Q_M = NULL;
        mx->l_max_Q_M = NULL;
        mx->k_min_Q_M = NULL;
        mx->k_max_Q_M = NULL;
        mx->Q_M_rem   = NULL;

        mx->Q_M1        = NULL;
        mx->l_min_Q_M1  = NULL;
        mx->l_max_Q_M1  = NULL;
        mx->k_min_Q_M1  = NULL;
        mx->k_max_Q_M1  = NULL;
        mx->Q_M1_rem    = NULL;

        mx->Q_M2        = NULL;
        mx->l_min_Q_M2  = NULL;
        mx->l_max_Q_M2  = NULL;
        mx->k_min_Q_M2  = NULL;
        mx->k_max_Q_M2  = NULL;
        mx->Q_M2_rem    = NULL;

        mx->Q_c       = NULL;
        mx->l_min_Q_c = NULL;
        mx->l_max_Q_c = NULL;
        mx->k_min_Q_c = 0;
        mx->k_max_Q_c = 0;
        mx->Q_c_rem   = 0.;

        mx->Q_cH        = NULL;
        mx->l_min_Q_cH  = NULL;
        mx->l_max_Q_cH  = NULL;
        mx->k_min_Q_cH  = 0;
        mx->k_max_Q_cH  = 0;
        mx->Q_cH_rem    = 0.;

        mx->Q_cI        = NULL;
        mx->l_min_Q_cI  = NULL;
        mx->l_max_Q_cI  = NULL;
        mx->k_min_Q_cI  = 0;
        mx->k_max_Q_cI  = 0;
        mx->Q_cI_rem    = 0.;

        mx->Q_cM        = NULL;
        mx->l_min_Q_cM  = NULL;
        mx->l_max_Q_cM  = NULL;
        mx->k_min_Q_cM  = 0;
        mx->k_max_Q_cM  = 0;
        mx->Q_cM_rem    = 0.;

        break;
    }
  }
}
