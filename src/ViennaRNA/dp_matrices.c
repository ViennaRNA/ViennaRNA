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
#include <math.h>

#include "data_structures.h"
#include "model.h"
#include "utils.h"
#include "gquad.h"
#include "dp_matrices.h"

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
#define ALLOC_HYBRID      2048
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
PRIVATE unsigned int    get_mx_alloc_vector(vrna_md_t       *md_p,
                                            vrna_mx_type_e  type,
                                            unsigned int    options);


PRIVATE unsigned int    get_mx_mfe_alloc_vector_current(vrna_mx_mfe_t   *mx,
                                                        vrna_mx_type_e  mx_type);


PRIVATE unsigned int    get_mx_pf_alloc_vector_current(vrna_mx_pf_t   *mx,
                                                       vrna_mx_type_e mx_type);


PRIVATE void            mfe_matrices_alloc_default(vrna_mx_mfe_t  *vars,
                                                   unsigned int   m,
                                                   unsigned int   alloc_vector);


PRIVATE void            mfe_matrices_free_default(vrna_mx_mfe_t *self);


PRIVATE void            mfe_matrices_alloc_window(vrna_mx_mfe_t *vars,
                                                  unsigned int  m,
                                                  unsigned int  alloc_vector);


PRIVATE void            mfe_matrices_free_window(vrna_mx_mfe_t  *self,
                                                 unsigned int   length,
                                                 unsigned int   window_size);


PRIVATE void            mfe_matrices_alloc_2Dfold(vrna_mx_mfe_t *vars,
                                                  unsigned int  m,
                                                  unsigned int  alloc_vector);


PRIVATE void            mfe_matrices_free_2Dfold(vrna_mx_mfe_t  *self,
                                                 unsigned int   length,
                                                 int            *indx);


PRIVATE void            pf_matrices_alloc_default(vrna_mx_pf_t  *vars,
                                                  unsigned int  m,
                                                  unsigned int  alloc_vector);


PRIVATE void            pf_matrices_free_default(vrna_mx_pf_t *self);


PRIVATE void            pf_matrices_alloc_window(vrna_mx_pf_t *vars,
                                                 unsigned int m,
                                                 unsigned int alloc_vector);


PRIVATE void            pf_matrices_free_window(vrna_mx_pf_t  *self,
                                                unsigned int  length,
                                                unsigned int  window_size);


PRIVATE void            pf_matrices_alloc_2Dfold(vrna_mx_pf_t *vars,
                                                 unsigned int m,
                                                 unsigned int alloc_vector);


PRIVATE void            pf_matrices_free_2Dfold(vrna_mx_pf_t  *self,
                                                unsigned int  length,
                                                int           *indx,
                                                int           *jindx);


PRIVATE vrna_mx_mfe_t *get_mfe_matrices_alloc(unsigned int    n,
                                              unsigned int    m,
                                              vrna_mx_type_e  type,
                                              unsigned int    alloc_vector);


PRIVATE vrna_mx_pf_t *get_pf_matrices_alloc(unsigned int    n,
                                            unsigned int    m,
                                            vrna_mx_type_e  type,
                                            unsigned int    alloc_vector);


PRIVATE void            add_pf_matrices(vrna_fold_compound_t  *vc,
                                        vrna_mx_type_e        type,
                                        unsigned int          alloc_vector);


PRIVATE void            add_mfe_matrices(vrna_fold_compound_t *vc,
                                         vrna_mx_type_e       type,
                                         unsigned int         alloc_vector);


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
          mfe_matrices_free_2Dfold(self, vc->length, vc->iindx);
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
          pf_matrices_free_2Dfold(self, vc->length, vc->iindx, vc->jindx);
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
    if (vc->cutpoint > 0)
      options |= VRNA_OPTION_HYBRID;

    mx_alloc_vector = get_mx_alloc_vector(&(vc->params->model_details), mx_type, options);
    vrna_mx_mfe_free(vc);
    add_mfe_matrices(vc, mx_type, mx_alloc_vector);
  } else {
    return 0;
  }

  return 1;
}


PUBLIC int
vrna_mx_pf_add(vrna_fold_compound_t *vc,
               vrna_mx_type_e       mx_type,
               unsigned int         options)
{
  unsigned int mx_alloc_vector;

  if (vc->exp_params) {
    mx_alloc_vector = get_mx_alloc_vector(&(vc->exp_params->model_details),
                                          mx_type,
                                          options | VRNA_OPTION_PF);
    vrna_mx_pf_free(vc);
    add_pf_matrices(vc, mx_type, mx_alloc_vector);
  } else {
    return 0;
  }

  return 1;
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

      if (vc->cutpoint > 0)
        options |= VRNA_OPTION_HYBRID;

      realloc = 0;

      if (!vc->matrices || (vc->matrices->type != mx_type) || (vc->matrices->length < vc->length)) {
        realloc = 1;
      } else {
        mx_alloc_vector =
          get_mx_alloc_vector(&(vc->params->model_details), mx_type, options);
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

      if (vc->cutpoint > 0)
        options |= VRNA_OPTION_HYBRID;

      realloc = 0;

      /*  Add DP matrices, if not they are not present */
      if (!vc->exp_matrices || (vc->exp_matrices->type != mx_type) ||
          (vc->exp_matrices->length < vc->length)) {
        realloc = 1;
      } else {
        mx_alloc_vector = get_mx_alloc_vector(&(vc->exp_params->model_details),
                                              mx_type,
                                              options);
        mx_alloc_vector_current = get_mx_pf_alloc_vector_current(vc->exp_matrices, mx_type);
        if ((mx_alloc_vector & mx_alloc_vector_current) != mx_alloc_vector)
          realloc = 1;
      }

      if (realloc) /* Add DP matrices, if not they are not present */
        ret &= vrna_mx_pf_add(vc, mx_type, options);

#ifdef VRNA_BACKWARD_COMPAT
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

        if (mx->fc)
          mx_alloc_vector |= ALLOC_HYBRID;

        if (mx->c)
          mx_alloc_vector |= ALLOC_C;

        if (mx->fML)
          mx_alloc_vector |= ALLOC_FML;

        if (mx->fM1)
          mx_alloc_vector |= ALLOC_UNIQ;

        if (mx->fM2)
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

        if (mx->qm1)
          mx_alloc_vector |= ALLOC_UNIQ;

        if (mx->qm2)
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


PRIVATE void
add_pf_matrices(vrna_fold_compound_t  *vc,
                vrna_mx_type_e        mx_type,
                unsigned int          alloc_vector)
{
  if (vc) {
    switch (mx_type) {
      case VRNA_MX_WINDOW:
        vc->exp_matrices =
          get_pf_matrices_alloc(vc->length, vc->window_size, mx_type, alloc_vector);
        break;
      default:
        vc->exp_matrices = get_pf_matrices_alloc(vc->length, vc->length, mx_type, alloc_vector);
        break;
    }

    if (vc->exp_params->model_details.gquad) {
      switch (vc->type) {
        case VRNA_FC_TYPE_SINGLE:
          vc->exp_matrices->G = NULL;
          /* can't do that here, since scale[] is not filled yet :(
           * vc->exp_matrices->G = get_gquad_pf_matrix(vc->sequence_encoding2, vc->exp_matrices->scale, vc->exp_params);
           */
          break;
        default:                    /* do nothing */
          break;
      }
    }

    vrna_exp_params_rescale(vc, NULL);
  }
}


PRIVATE void
add_mfe_matrices(vrna_fold_compound_t *vc,
                 vrna_mx_type_e       mx_type,
                 unsigned int         alloc_vector)
{
  if (vc) {
    switch (mx_type) {
      case VRNA_MX_WINDOW:
        vc->matrices = get_mfe_matrices_alloc(vc->length, vc->window_size, mx_type, alloc_vector);
        break;
      default:
        vc->matrices = get_mfe_matrices_alloc(vc->length, vc->length, mx_type, alloc_vector);
        break;
    }

    if (vc->params->model_details.gquad) {
      switch (vc->type) {
        case VRNA_FC_TYPE_SINGLE:
          switch (mx_type) {
            case VRNA_MX_WINDOW:                              /* do nothing, since we handle memory somewhere else */
              break;
            default:
              vc->matrices->ggg = get_gquad_matrix(vc->sequence_encoding2, vc->params);
              break;
          }
          break;
        case VRNA_FC_TYPE_COMPARATIVE:
          switch (mx_type) {
            case VRNA_MX_WINDOW:                              /* do nothing, since we handle memory somewhere else */
              break;
            default:
              vc->matrices->ggg = get_gquad_ali_matrix(vc->S_cons, vc->S, vc->n_seq, vc->params);
              break;
          }
          break;
        default:                      /* do nothing */
          break;
      }
    }
  }
}


PRIVATE vrna_mx_mfe_t *
get_mfe_matrices_alloc(unsigned int   n,
                       unsigned int   m,
                       vrna_mx_type_e type,
                       unsigned int   alloc_vector)
{
  vrna_mx_mfe_t *vars;

  if ((int)(n * m) >= (int)INT_MAX)
    vrna_message_error(
      "get_mfe_matrices_alloc@data_structures.c: sequence length exceeds addressable range");

  vars          = (vrna_mx_mfe_t *)vrna_alloc(sizeof(vrna_mx_mfe_t));
  vars->length  = n;
  vars->type    = type;

  switch (type) {
    case VRNA_MX_DEFAULT:
      mfe_matrices_alloc_default(vars, m, alloc_vector);
      break;

    case VRNA_MX_WINDOW:
      mfe_matrices_alloc_window(vars, m, alloc_vector);
      break;

    case VRNA_MX_2DFOLD:
      mfe_matrices_alloc_2Dfold(vars, m, alloc_vector);
      break;

    default:                /* do nothing */
      break;
  }

  return vars;
}


PRIVATE vrna_mx_pf_t *
get_pf_matrices_alloc(unsigned int    n,
                      unsigned int    m,
                      vrna_mx_type_e  type,
                      unsigned int    alloc_vector)
{
  unsigned int  lin_size;
  vrna_mx_pf_t  *vars;

  if ((int)(n * m) >= (int)INT_MAX)
    vrna_message_error(
      "get_pf_matrices_alloc@data_structures.c: sequence length exceeds addressable range");

  lin_size      = n + 2;
  vars          = (vrna_mx_pf_t *)vrna_alloc(sizeof(vrna_mx_pf_t));
  vars->length  = n;
  vars->type    = type;


  switch (type) {
    case VRNA_MX_DEFAULT:
      pf_matrices_alloc_default(vars, n, alloc_vector);
      break;

    case VRNA_MX_WINDOW:
      pf_matrices_alloc_window(vars, m, alloc_vector);
      break;

    case VRNA_MX_2DFOLD:
      pf_matrices_alloc_2Dfold(vars, n, alloc_vector);
      break;

    default:                /* do nothing */
      break;
  }

  /*
   *  always alloc the helper arrays for unpaired nucleotides in multi-
   *  branch loops and scaling
   */
  vars->scale     = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * lin_size);
  vars->expMLbase = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * lin_size);

  return vars;
}


PRIVATE unsigned int
get_mx_alloc_vector(vrna_md_t       *md_p,
                    vrna_mx_type_e  mx_type,
                    unsigned int    options)
{
  unsigned int v;

  v = ALLOC_NOTHING;

  /* default MFE matrices ? */
  if (options & VRNA_OPTION_MFE)
    v |= (mx_type == VRNA_MX_WINDOW) ? ALLOC_MFE_LOCAL : ALLOC_MFE_DEFAULT;

  /* default PF matrices ? */
  if (options & VRNA_OPTION_PF)
    v |= (md_p->compute_bpp) ? ALLOC_PF_DEFAULT : ALLOC_PF_WO_PROBS;

  if (options & VRNA_OPTION_HYBRID)
    v |= ALLOC_HYBRID;

  /* matrices for circular folding ? */
  if (md_p->circ) {
    md_p->uniq_ML = 1; /* we need unique ML arrays for circular folding */
    v             |= ALLOC_CIRC;
  }

  /* unique ML decomposition ? */
  if (md_p->uniq_ML)
    v |= ALLOC_UNIQ;

  return v;
}


PRIVATE void
mfe_matrices_alloc_default(vrna_mx_mfe_t  *vars,
                           unsigned int   m,
                           unsigned int   alloc_vector)
{
  unsigned int n, size, lin_size;

  n         = vars->length;
  size      = ((n + 1) * (m + 2)) / 2;
  lin_size  = n + 2;

  vars->f5  = NULL;
  vars->f3  = NULL;
  vars->fc  = NULL;
  vars->c   = NULL;
  vars->fML = NULL;
  vars->fM1 = NULL;
  vars->fM2 = NULL;
  vars->ggg = NULL;

  if (alloc_vector & ALLOC_F5)
    vars->f5 = (int *)vrna_alloc(sizeof(int) * lin_size);

  if (alloc_vector & ALLOC_F3)
    vars->f3 = (int *)vrna_alloc(sizeof(int) * lin_size);

  if (alloc_vector & ALLOC_HYBRID)
    vars->fc = (int *)vrna_alloc(sizeof(int) * lin_size);

  if (alloc_vector & ALLOC_C)
    vars->c = (int *)vrna_alloc(sizeof(int) * size);

  if (alloc_vector & ALLOC_FML)
    vars->fML = (int *)vrna_alloc(sizeof(int) * size);

  if (alloc_vector & ALLOC_UNIQ)
    vars->fM1 = (int *)vrna_alloc(sizeof(int) * size);

  if (alloc_vector & ALLOC_CIRC)
    vars->fM2 = (int *)vrna_alloc(sizeof(int) * lin_size);

  /* setting exterior loop energies for circular case to INF is always safe */
  vars->FcH = vars->FcI = vars->FcM = vars->Fc = INF;
}


PRIVATE void
mfe_matrices_free_default(vrna_mx_mfe_t *self)
{
  free(self->f5);
  free(self->f3);
  free(self->fc);
  free(self->c);
  free(self->fML);
  free(self->fM1);
  free(self->fM2);
  free(self->ggg);
}


PRIVATE void
mfe_matrices_alloc_window(vrna_mx_mfe_t *vars,
                          unsigned int  m,
                          unsigned int  alloc_vector)
{
  int           i;
  unsigned int  n, lin_size;

  n         = vars->length;
  lin_size  = n + 2;

  vars->f3_local  = NULL;
  vars->c_local   = NULL;
  vars->fML_local = NULL;
  vars->ggg_local = NULL;

  if (alloc_vector & ALLOC_F3)
    vars->f3_local = (int *)vrna_alloc(sizeof(int) * lin_size);

  if (alloc_vector & ALLOC_C)
    vars->c_local = (int **)vrna_alloc(sizeof(int *) * lin_size);

  if (alloc_vector & ALLOC_FML)
    vars->fML_local = (int **)vrna_alloc(sizeof(int *) * lin_size);
}


PRIVATE void
mfe_matrices_free_window(vrna_mx_mfe_t  *self,
                         unsigned int   length,
                         unsigned int   window_size)
{
  free(self->c_local);
  free(self->fML_local);
  free(self->ggg_local);
  free(self->f3_local);
}


PRIVATE void
mfe_matrices_alloc_2Dfold(vrna_mx_mfe_t *vars,
                          unsigned int  m,
                          unsigned int  alloc_vector)
{
  unsigned int n, i, size, lin_size;

  n         = vars->length;
  size      = ((n + 1) * (m + 2)) / 2;
  lin_size  = n + 2;

  vars->E_F5      = NULL;
  vars->l_min_F5  = NULL;
  vars->l_max_F5  = NULL;
  vars->k_min_F5  = NULL;
  vars->k_max_F5  = NULL;
  vars->E_F5_rem  = NULL;

  vars->E_F3      = NULL;
  vars->l_min_F3  = NULL;
  vars->l_max_F3  = NULL;
  vars->k_min_F3  = NULL;
  vars->k_max_F3  = NULL;
  vars->E_F3_rem  = NULL;

  vars->E_C     = NULL;
  vars->l_min_C = NULL;
  vars->l_max_C = NULL;
  vars->k_min_C = NULL;
  vars->k_max_C = NULL;
  vars->E_C_rem = NULL;

  vars->E_M     = NULL;
  vars->l_min_M = NULL;
  vars->l_max_M = NULL;
  vars->k_min_M = NULL;
  vars->k_max_M = NULL;
  vars->E_M_rem = NULL;

  vars->E_M1      = NULL;
  vars->l_min_M1  = NULL;
  vars->l_max_M1  = NULL;
  vars->k_min_M1  = NULL;
  vars->k_max_M1  = NULL;
  vars->E_M1_rem  = NULL;

  vars->E_M2      = NULL;
  vars->l_min_M2  = NULL;
  vars->l_max_M2  = NULL;
  vars->k_min_M2  = NULL;
  vars->k_max_M2  = NULL;
  vars->E_M2_rem  = NULL;

  /* setting exterior loop energies for circular case to INF is always safe */
  vars->E_Fc      = NULL;
  vars->E_FcH     = NULL;
  vars->E_FcI     = NULL;
  vars->E_FcM     = NULL;
  vars->E_Fc_rem  = INF;
  vars->E_FcH_rem = INF;
  vars->E_FcI_rem = INF;
  vars->E_FcM_rem = INF;

  if (alloc_vector & ALLOC_F5) {
    vars->E_F5      = (int ***)vrna_alloc(sizeof(int **) * lin_size);
    vars->l_min_F5  = (int **)vrna_alloc(sizeof(int *) * lin_size);
    vars->l_max_F5  = (int **)vrna_alloc(sizeof(int *) * lin_size);
    vars->k_min_F5  = (int *)vrna_alloc(sizeof(int) * lin_size);
    vars->k_max_F5  = (int *)vrna_alloc(sizeof(int) * lin_size);
    vars->E_F5_rem  = (int *)vrna_alloc(sizeof(int) * lin_size);
    for (i = 0; i <= n; i++)
      vars->E_F5_rem[i] = INF;
  }

  if (alloc_vector & ALLOC_F3) {
    vars->E_F3      = (int ***)vrna_alloc(sizeof(int **) * lin_size);
    vars->l_min_F3  = (int **)vrna_alloc(sizeof(int *) * lin_size);
    vars->l_max_F3  = (int **)vrna_alloc(sizeof(int *) * lin_size);
    vars->k_min_F3  = (int *)vrna_alloc(sizeof(int) * lin_size);
    vars->k_max_F3  = (int *)vrna_alloc(sizeof(int) * lin_size);
    vars->E_F3_rem  = (int *)vrna_alloc(sizeof(int) * lin_size);
    for (i = 0; i <= n; i++)
      vars->E_F3_rem[i] = INF;
  }

  if (alloc_vector & ALLOC_C) {
    vars->E_C     = (int ***)vrna_alloc(sizeof(int **) * size);
    vars->l_min_C = (int **)vrna_alloc(sizeof(int *) * size);
    vars->l_max_C = (int **)vrna_alloc(sizeof(int *) * size);
    vars->k_min_C = (int *)vrna_alloc(sizeof(int) * size);
    vars->k_max_C = (int *)vrna_alloc(sizeof(int) * size);
    vars->E_C_rem = (int *)vrna_alloc(sizeof(int) * size);
    for (i = 0; i < size; i++)
      vars->E_C_rem[i] = INF;
  }

  if (alloc_vector & ALLOC_FML) {
    vars->E_M     = (int ***)vrna_alloc(sizeof(int **) * size);
    vars->l_min_M = (int **)vrna_alloc(sizeof(int *) * size);
    vars->l_max_M = (int **)vrna_alloc(sizeof(int *) * size);
    vars->k_min_M = (int *)vrna_alloc(sizeof(int) * size);
    vars->k_max_M = (int *)vrna_alloc(sizeof(int) * size);
    vars->E_M_rem = (int *)vrna_alloc(sizeof(int) * size);
    for (i = 0; i < size; i++)
      vars->E_M_rem[i] = INF;
  }

  if (alloc_vector & ALLOC_UNIQ) {
    vars->E_M1      = (int ***)vrna_alloc(sizeof(int **) * size);
    vars->l_min_M1  = (int **)vrna_alloc(sizeof(int *) * size);
    vars->l_max_M1  = (int **)vrna_alloc(sizeof(int *) * size);
    vars->k_min_M1  = (int *)vrna_alloc(sizeof(int) * size);
    vars->k_max_M1  = (int *)vrna_alloc(sizeof(int) * size);
    vars->E_M1_rem  = (int *)vrna_alloc(sizeof(int) * size);
    for (i = 0; i < size; i++)
      vars->E_M1_rem[i] = INF;
  }

  if (alloc_vector & ALLOC_CIRC) {
    vars->E_M2      = (int ***)vrna_alloc(sizeof(int **) * lin_size);
    vars->l_min_M2  = (int **)vrna_alloc(sizeof(int *) * lin_size);
    vars->l_max_M2  = (int **)vrna_alloc(sizeof(int *) * lin_size);
    vars->k_min_M2  = (int *)vrna_alloc(sizeof(int) * lin_size);
    vars->k_max_M2  = (int *)vrna_alloc(sizeof(int) * lin_size);
    vars->E_M2_rem  = (int *)vrna_alloc(sizeof(int) * lin_size);
    for (i = 0; i <= n; i++)
      vars->E_M2_rem[i] = INF;
  }

#ifdef COUNT_STATES
  vars->N_C   = (unsigned long ***)vrna_alloc(sizeof(unsigned long **) * size);
  vars->N_F5  = (unsigned long ***)vrna_alloc(sizeof(unsigned long **) * lin_size);
  vars->N_M   = (unsigned long ***)vrna_alloc(sizeof(unsigned long **) * size);
  vars->N_M1  = (unsigned long ***)vrna_alloc(sizeof(unsigned long **) * size);
#endif
}


PRIVATE void
mfe_matrices_free_2Dfold(vrna_mx_mfe_t  *self,
                         unsigned int   length,
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
    for (i = 1; i < length - TURN - 1; i++) {
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
pf_matrices_alloc_default(vrna_mx_pf_t  *vars,
                          unsigned int  m,
                          unsigned int  alloc_vector)
{
  unsigned int n, size, lin_size;

  n         = vars->length;
  size      = ((n + 1) * (n + 2)) / 2;
  lin_size  = n + 2;

  vars->q     = NULL;
  vars->qb    = NULL;
  vars->qm    = NULL;
  vars->qm1   = NULL;
  vars->qm2   = NULL;
  vars->probs = NULL;
  vars->q1k   = NULL;
  vars->qln   = NULL;

  if (alloc_vector & ALLOC_F)
    vars->q = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * size);

  if (alloc_vector & ALLOC_C)
    vars->qb = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * size);

  if (alloc_vector & ALLOC_FML)
    vars->qm = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * size);

  if (alloc_vector & ALLOC_UNIQ)
    vars->qm1 = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * size);

  if (alloc_vector & ALLOC_CIRC)
    vars->qm2 = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * lin_size);

  if (alloc_vector & ALLOC_PROBS)
    vars->probs = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * size);

  if (alloc_vector & ALLOC_AUX) {
    vars->q1k = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * lin_size);
    vars->qln = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * lin_size);
  }
}


PRIVATE void
pf_matrices_free_default(vrna_mx_pf_t *self)
{
  free(self->q);
  free(self->qb);
  free(self->qm);
  free(self->qm1);
  free(self->qm2);
  free(self->probs);
  free(self->G);
  free(self->q1k);
  free(self->qln);
}


PRIVATE void
pf_matrices_alloc_window(vrna_mx_pf_t *vars,
                         unsigned int m,
                         unsigned int alloc_vector)
{
  int           i;
  unsigned int  n, lin_size;

  n         = vars->length;
  lin_size  = n + 2;

  vars->q_local   = NULL;
  vars->qb_local  = NULL;
  vars->qm_local  = NULL;
  vars->qm2_local = NULL;
  vars->pR        = NULL;
  vars->QI5       = NULL;
  vars->q2l       = NULL;
  vars->qmb       = NULL;

  if (alloc_vector & ALLOC_F)
    vars->q_local = (FLT_OR_DBL **)vrna_alloc(sizeof(FLT_OR_DBL *) * lin_size);

  if (alloc_vector & ALLOC_C)
    vars->qb_local = (FLT_OR_DBL **)vrna_alloc(sizeof(FLT_OR_DBL *) * lin_size);

  if (alloc_vector & ALLOC_FML)
    vars->qm_local = (FLT_OR_DBL **)vrna_alloc(sizeof(FLT_OR_DBL *) * lin_size);

  vars->pR = (FLT_OR_DBL **)vrna_alloc(sizeof(FLT_OR_DBL *) * lin_size);

  if (alloc_vector & ALLOC_PROBS) {
    vars->QI5       = (FLT_OR_DBL **)vrna_alloc(sizeof(FLT_OR_DBL *) * lin_size);
    vars->qmb       = (FLT_OR_DBL **)vrna_alloc(sizeof(FLT_OR_DBL *) * lin_size);
    vars->qm2_local = (FLT_OR_DBL **)vrna_alloc(sizeof(FLT_OR_DBL *) * lin_size);
    vars->q2l       = (FLT_OR_DBL **)vrna_alloc(sizeof(FLT_OR_DBL *) * lin_size);
  }
}


PRIVATE void
pf_matrices_free_window(vrna_mx_pf_t  *self,
                        unsigned int  length,
                        unsigned int  window_size)
{
  free(self->q_local);
  free(self->qb_local);
  free(self->qm_local);
  free(self->qm2_local);
  free(self->pR);
  free(self->QI5);
  free(self->q2l);
  free(self->qmb);
}


PRIVATE void
pf_matrices_alloc_2Dfold(vrna_mx_pf_t *vars,
                         unsigned int m,
                         unsigned int alloc_vector)
{
  unsigned int n, size, lin_size;

  n         = vars->length;
  size      = ((n + 1) * (n + 2)) / 2;
  lin_size  = n + 2;

  vars->Q       = NULL;
  vars->l_min_Q = NULL;
  vars->l_max_Q = NULL;
  vars->k_min_Q = NULL;
  vars->k_max_Q = NULL;
  vars->Q_rem   = NULL;

  vars->Q_B       = NULL;
  vars->l_min_Q_B = NULL;
  vars->l_max_Q_B = NULL;
  vars->k_min_Q_B = NULL;
  vars->k_max_Q_B = NULL;
  vars->Q_B_rem   = NULL;

  vars->Q_M       = NULL;
  vars->l_min_Q_M = NULL;
  vars->l_max_Q_M = NULL;
  vars->k_min_Q_M = NULL;
  vars->k_max_Q_M = NULL;
  vars->Q_M_rem   = NULL;

  vars->Q_M1        = NULL;
  vars->l_min_Q_M1  = NULL;
  vars->l_max_Q_M1  = NULL;
  vars->k_min_Q_M1  = NULL;
  vars->k_max_Q_M1  = NULL;
  vars->Q_M1_rem    = NULL;

  vars->Q_M2        = NULL;
  vars->l_min_Q_M2  = NULL;
  vars->l_max_Q_M2  = NULL;
  vars->k_min_Q_M2  = NULL;
  vars->k_max_Q_M2  = NULL;
  vars->Q_M2_rem    = NULL;

  vars->Q_c       = NULL;
  vars->Q_cH      = NULL;
  vars->Q_cI      = NULL;
  vars->Q_cM      = NULL;
  vars->Q_c_rem   = 0.;
  vars->Q_cH_rem  = 0.;
  vars->Q_cI_rem  = 0.;
  vars->Q_cM_rem  = 0.;

  if (alloc_vector & ALLOC_F) {
    vars->Q       = (FLT_OR_DBL ***)vrna_alloc(sizeof(FLT_OR_DBL **) * size);
    vars->l_min_Q = (int **)vrna_alloc(sizeof(int *) * size);
    vars->l_max_Q = (int **)vrna_alloc(sizeof(int *) * size);
    vars->k_min_Q = (int *)vrna_alloc(sizeof(int) * size);
    vars->k_max_Q = (int *)vrna_alloc(sizeof(int) * size);
    vars->Q_rem   = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * size);
  }

  if (alloc_vector & ALLOC_C) {
    vars->Q_B       = (FLT_OR_DBL ***)vrna_alloc(sizeof(FLT_OR_DBL **) * size);
    vars->l_min_Q_B = (int **)vrna_alloc(sizeof(int *) * size);
    vars->l_max_Q_B = (int **)vrna_alloc(sizeof(int *) * size);
    vars->k_min_Q_B = (int *)vrna_alloc(sizeof(int) * size);
    vars->k_max_Q_B = (int *)vrna_alloc(sizeof(int) * size);
    vars->Q_B_rem   = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * size);
  }

  if (alloc_vector & ALLOC_FML) {
    vars->Q_M       = (FLT_OR_DBL ***)vrna_alloc(sizeof(FLT_OR_DBL **) * size);
    vars->l_min_Q_M = (int **)vrna_alloc(sizeof(int *) * size);
    vars->l_max_Q_M = (int **)vrna_alloc(sizeof(int *) * size);
    vars->k_min_Q_M = (int *)vrna_alloc(sizeof(int) * size);
    vars->k_max_Q_M = (int *)vrna_alloc(sizeof(int) * size);
    vars->Q_M_rem   = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * size);
  }

  if (alloc_vector & ALLOC_UNIQ) {
    vars->Q_M1        = (FLT_OR_DBL ***)vrna_alloc(sizeof(FLT_OR_DBL **) * size);
    vars->l_min_Q_M1  = (int **)vrna_alloc(sizeof(int *) * size);
    vars->l_max_Q_M1  = (int **)vrna_alloc(sizeof(int *) * size);
    vars->k_min_Q_M1  = (int *)vrna_alloc(sizeof(int) * size);
    vars->k_max_Q_M1  = (int *)vrna_alloc(sizeof(int) * size);
    vars->Q_M1_rem    = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * size);
  }

  if (alloc_vector & ALLOC_CIRC) {
    vars->Q_M2        = (FLT_OR_DBL ***)vrna_alloc(sizeof(FLT_OR_DBL **) * lin_size);
    vars->l_min_Q_M2  = (int **)vrna_alloc(sizeof(int *) * lin_size);
    vars->l_max_Q_M2  = (int **)vrna_alloc(sizeof(int *) * lin_size);
    vars->k_min_Q_M2  = (int *)vrna_alloc(sizeof(int) * lin_size);
    vars->k_max_Q_M2  = (int *)vrna_alloc(sizeof(int) * lin_size);
    vars->Q_M2_rem    = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * lin_size);
  }
}


PRIVATE void
pf_matrices_free_2Dfold(vrna_mx_pf_t  *self,
                        unsigned int  length,
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
    for (i = 1; i < length - TURN - 1; i++) {
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
