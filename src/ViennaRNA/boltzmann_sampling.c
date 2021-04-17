/*
 *                partiton function for RNA secondary structures
 *
 *                Ivo L Hofacker + Ronny Lorenz
 *                Vienna RNA package
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <math.h>

#include "ViennaRNA/utils/basic.h"
#include "ViennaRNA/params/default.h"
#include "ViennaRNA/fold_vars.h"
#include "ViennaRNA/loops/all.h"
#include "ViennaRNA/gquad.h"
#include "ViennaRNA/constraints/hard.h"
#include "ViennaRNA/constraints/soft.h"
#include "ViennaRNA/alphabet.h"
#include "ViennaRNA/boltzmann_sampling.h"

#include "ViennaRNA/loops/external_sc_pf.inc"
#include "ViennaRNA/loops/internal_sc_pf.inc"
#include "ViennaRNA/loops/multibranch_sc_pf.inc"

#include "ViennaRNA/data_structures_nonred.inc"

/*
 #################################
 # PREPROCESSOR DEFININTIONS     #
 #################################
 */

#ifdef VRNA_NR_SAMPLING_HASH
# define NR_NODE tr_node
# define NR_TOTAL_WEIGHT(a) total_weight_par(a)
# define NR_TOTAL_WEIGHT_TYPE(a, b) total_weight_par_type(a, b)
# define NR_GET_WEIGHT(a, b, c, d, e)  tr_node_weight(a, c, d, e)
#else
# define NR_NODE tllr_node
# define NR_TOTAL_WEIGHT(a) get_weight_all(a)
# define NR_TOTAL_WEIGHT_TYPE(a, b) get_weight_type_spec(a, b)
# define NR_GET_WEIGHT(a, b, c, d, e)  get_weight(b, c, d, e)
#endif


/* combination of soft constraint wrappers */
struct sc_wrappers {
  struct sc_ext_exp_dat sc_wrapper_ext;
  struct sc_int_exp_dat sc_wrapper_int;
  struct sc_mb_exp_dat  sc_wrapper_ml;
};

/*
 * In the following:
 * - q_remain is a pointer to value of sum of Boltzmann factors of still accessible solutions at that point
 * - current_node is a pointer to current node in datastructure memorizing the solutions and paths taken
 */
struct vrna_pbacktrack_memory_s {
  double            q_remain;
  NR_NODE           *root_node;
  NR_NODE           *current_node;
  struct nr_memory  *memory_dat;
};

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
PRIVATE char *info_set_uniq_ml =
  "Unique multiloop decomposition is unset!\n"
  "Activate unique multiloop decomposition by setting the"
  " uniq_ML field of the model details structure to a non-zero"
  " value before running vrna_pf()!";

PRIVATE char  *info_call_pf =
  "DP matrices are missing! Call vrna_pf() first!";


PRIVATE char  *info_nr_duplicates =
  "Duplicate structures detected, presumably due to numerical instabilities";


PRIVATE char  *info_nr_overflow =
  "Partition function overflow detected for forbidden structures,"
  " presumably due to numerical instabilities.";


PRIVATE char  *info_no_circ =
  "No implementation for circular RNAs available.";


/*
 #################################
 # PRIVATE FUNCTION DECLARATIONS #
 #################################
 */

PRIVATE struct vrna_pbacktrack_memory_s *
nr_init(vrna_fold_compound_t *fc);


PRIVATE struct sc_wrappers *
sc_init(vrna_fold_compound_t *fc);


PRIVATE void
sc_free(struct sc_wrappers *sc_wrap);


PRIVATE unsigned int
wrap_pbacktrack(vrna_fold_compound_t              *vc,
                unsigned int                      length,
                unsigned int                      num_samples,
                vrna_boltzmann_sampling_callback  *bs_cb,
                void                              *data,
                struct vrna_pbacktrack_memory_s   *nr_mem);


PRIVATE int
backtrack(int                             i,
          int                             j,
          char                            *pstruc,
          vrna_fold_compound_t            *vc,
          struct sc_wrappers              *sc_wrap,
          struct vrna_pbacktrack_memory_s *nr_mem);


PRIVATE int
backtrack_ext_loop(int                              init_val,
                   char                             *pstruc,
                   vrna_fold_compound_t             *vc,
                   int                              length,
                   struct sc_wrappers               *sc_wrap,
                   struct vrna_pbacktrack_memory_s  *nr_mem);


PRIVATE int
backtrack_qm(int                              i,
             int                              j,
             char                             *pstruc,
             vrna_fold_compound_t             *vc,
             struct sc_wrappers               *sc_wrap,
             struct vrna_pbacktrack_memory_s  *nr_mem);


PRIVATE int
backtrack_qm1(int                             i,
              int                             j,
              char                            *pstruc,
              vrna_fold_compound_t            *vc,
              struct sc_wrappers              *sc_wrap,
              struct vrna_pbacktrack_memory_s *nr_mem);


PRIVATE void
backtrack_qm2(int                   u,
              int                   n,
              char                  *pstruc,
              vrna_fold_compound_t  *vc,
              struct sc_wrappers    *sc_wrap);


PRIVATE unsigned int
pbacktrack_circ(vrna_fold_compound_t              *fc,
                unsigned int                      num_samples,
                vrna_boltzmann_sampling_callback  *bs_cb,
                void                              *data);


/*
 #################################
 # BEGIN OF FUNCTION DEFINITIONS #
 #################################
 */
PUBLIC unsigned int
vrna_pbacktrack5_resume_cb(vrna_fold_compound_t             *fc,
                           unsigned int                     num_samples,
                           unsigned int                     length,
                           vrna_boltzmann_sampling_callback *bs_cb,
                           void                             *data,
                           vrna_pbacktrack_mem_t            *nr_mem,
                           unsigned int                     options)
{
  unsigned int i = 0;

  if (fc) {
    vrna_mx_pf_t *matrices = fc->exp_matrices;

    if (length > fc->length) {
      vrna_message_warning("vrna_pbacktrack5*(): length exceeds sequence length");
    } else if (length == 0) {
      vrna_message_warning("vrna_pbacktrack5*(): length too small");
    } else if ((!matrices) || (!matrices->q) || (!matrices->qb) || (!matrices->qm) ||
               (!fc->exp_params)) {
      vrna_message_warning("vrna_pbacktrack*(): %s", info_call_pf);
    } else if ((!fc->exp_params->model_details.uniq_ML) || (!matrices->qm1)) {
      vrna_message_warning("vrna_pbacktrack*(): %s", info_set_uniq_ml);
    } else if ((fc->exp_params->model_details.circ) && (length < fc->length)) {
      vrna_message_warning("vrna_pbacktrack5*(): %s", info_no_circ);
    } else if (options & VRNA_PBACKTRACK_NON_REDUNDANT) {
      if (fc->exp_params->model_details.circ) {
        vrna_message_warning("vrna_pbacktrack5*(): %s", info_no_circ);
      } else if (!nr_mem) {
        vrna_message_warning("vrna_pbacktrack5*(): Pointer to nr_mem must not be NULL!");
      } else {
        if (*nr_mem == NULL)
          *nr_mem = nr_init(fc);

        i = wrap_pbacktrack(fc, length, num_samples, bs_cb, data, *nr_mem);

        /* print warning if we've aborted backtracking too early */
        if ((i > 0) && (i < num_samples)) {
          vrna_message_warning("vrna_pbacktrack5*(): "
                               "Stopped non-redundant backtracking after %d samples"
                               " due to numeric instabilities!\n"
                               "Coverage of partition function so far: %.6f%%",
                               i,
                               100. *
                               return_node_weight((*nr_mem)->root_node) /
                               fc->exp_matrices->q[fc->iindx[1] - length]);
        }
      }
    } else if (fc->exp_params->model_details.circ) {
      i = pbacktrack_circ(fc, num_samples, bs_cb, data);
    } else {
      i = wrap_pbacktrack(fc, length, num_samples, bs_cb, data, NULL);
    }
  }

  return i; /* actual number of structures backtraced */
}


PUBLIC void
vrna_pbacktrack_mem_free(struct vrna_pbacktrack_memory_s *s)
{
  if (s) {
#ifdef VRNA_NR_SAMPLING_HASH
    free_all_nr(s->current_node);
#else
    free_all_nrll(&(s->memory_dat));
#endif
    free(s);
  }
}


/*
 #####################################
 # BEGIN OF STATIC HELPER FUNCTIONS  #
 #####################################
 */
PRIVATE struct sc_wrappers *
sc_init(vrna_fold_compound_t *fc)
{
  struct sc_wrappers *sc_wrap = (struct sc_wrappers *)vrna_alloc(sizeof(struct sc_wrappers));

  init_sc_ext_exp(fc, &(sc_wrap->sc_wrapper_ext));
  init_sc_int_exp(fc, &(sc_wrap->sc_wrapper_int));
  init_sc_mb_exp(fc, &(sc_wrap->sc_wrapper_ml));

  return sc_wrap;
}


PRIVATE void
sc_free(struct sc_wrappers *sc_wrap)
{
  free_sc_ext_exp(&(sc_wrap->sc_wrapper_ext));
  free_sc_int_exp(&(sc_wrap->sc_wrapper_int));
  free_sc_mb_exp(&(sc_wrap->sc_wrapper_ml));

  free(sc_wrap);
}


PRIVATE struct vrna_pbacktrack_memory_s *
nr_init(vrna_fold_compound_t *fc)
{
  size_t                          block_size;
  double                          pf;
  struct vrna_pbacktrack_memory_s *s;

  s = (struct vrna_pbacktrack_memory_s *)vrna_alloc(
    sizeof(struct vrna_pbacktrack_memory_s));
  pf          = fc->exp_matrices->q[fc->iindx[1] - fc->length];
  block_size  = 5000 * sizeof(NR_NODE);

  s->memory_dat = NULL;
  s->q_remain   = 0;

#ifdef VRNA_NR_SAMPLING_HASH
  s->root_node = create_root(fc->length, pf);
#else
  s->memory_dat = create_nr_memory(sizeof(NR_NODE), block_size, NULL);  /* memory pre-allocation */
  s->root_node  = create_ll_root(&(s->memory_dat), pf);
#endif

  s->current_node = s->root_node;

  return s;
}


/* general expr of vrna5_pbacktrack with possibility of non-redundant sampling */
PRIVATE unsigned int
wrap_pbacktrack(vrna_fold_compound_t              *vc,
                unsigned int                      length,
                unsigned int                      num_samples,
                vrna_boltzmann_sampling_callback  *bs_cb,
                void                              *data,
                struct vrna_pbacktrack_memory_s   *nr_mem)
{
  char                *pstruc;
  unsigned int        i, n;
  int                 ret, pf_overflow, is_dup, *my_iindx;
  FLT_OR_DBL          *q1k, *qln, *q;
  vrna_mx_pf_t        *matrices;
  struct sc_wrappers  *sc_wrap;

  i           = 0;
  pf_overflow = 0;
  sc_wrap     = sc_init(vc);

  n         = vc->length;
  my_iindx  = vc->iindx;
  matrices  = vc->exp_matrices;
  q         = matrices->q;
  q1k       = matrices->q1k;
  qln       = matrices->qln;

  if (!(q1k && qln)) {
    matrices->q1k = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * (n + 1));
    matrices->qln = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * (n + 2));
    q1k           = matrices->q1k;
    qln           = matrices->qln;
    for (i = 1; i <= n; i++) {
      q1k[i]  = q[my_iindx[1] - i];
      qln[i]  = q[my_iindx[i] - n];
    }
    q1k[0]      = 1.0;
    qln[n + 1]  = 1.0;
  }

  for (i = 0; i < num_samples; i++) {
    is_dup  = 1;
    pstruc  = vrna_alloc((length + 1) * sizeof(char));

    memset(pstruc, '.', sizeof(char) * length);

    if (nr_mem)
      nr_mem->q_remain = vc->exp_matrices->q[vc->iindx[1] - length];

#ifdef VRNA_WITH_BOUSTROPHEDON
    ret = backtrack_ext_loop(length, pstruc, vc, length, sc_wrap, nr_mem);
#else
    ret = backtrack_ext_loop(1, pstruc, vc, length, sc_wrap, nr_mem);
#endif

    if (nr_mem) {
#ifdef VRNA_NR_SAMPLING_HASH
      nr_mem->current_node = traceback_to_root(nr_mem->current_node,
                                               nr_mem->q_remain,
                                               &is_dup,
                                               &pf_overflow);
#else
      nr_mem->current_node = traceback_to_ll_root(nr_mem->current_node,
                                                  nr_mem->q_remain,
                                                  &is_dup,
                                                  &pf_overflow);
#endif

      if (pf_overflow) {
        vrna_message_warning("vrna_pbacktrack_nr*(): %s", info_nr_overflow);
        free(pstruc);
        break;
      }

      if (is_dup) {
        vrna_message_warning("vrna_pbacktrack_nr*(): %s", info_nr_duplicates);
        free(pstruc);
        break;
      }
    }

    if ((ret > 0) && (bs_cb))
      bs_cb(pstruc, data);

    free(pstruc);

    if (ret == 0)
      break;
  }

  sc_free(sc_wrap);

  return i;
}


/* backtrack one external */
PRIVATE int
backtrack_ext_loop(int                              init_val,
                   char                             *pstruc,
                   vrna_fold_compound_t             *vc,
                   int                              length,
                   struct sc_wrappers               *sc_wrap,
                   struct vrna_pbacktrack_memory_s  *nr_mem)
{
  unsigned char             *hard_constraints;
  short                     *S1, *S2, **S, **S5, **S3;
  unsigned int              **a2s, s, n_seq;
  int                       ret, i, j, ij, n, k, u, type, *my_iindx, hc_decompose, *hc_up_ext;
  FLT_OR_DBL                r, fbd, fbds, qt, q_temp, qkl, *q, *qb, *q1k, *qln, *scale;
  double                    *q_remain;
  vrna_mx_pf_t              *matrices;
  vrna_md_t                 *md;
  vrna_hc_t                 *hc;
  vrna_exp_param_t          *pf_params;

  struct nr_memory          **memory_dat;
  struct sc_ext_exp_dat     *sc_wrapper_ext;

  NR_NODE                   **current_node;

  if (nr_mem) {
    q_remain      = &(nr_mem->q_remain);
    current_node  = &(nr_mem->current_node);
    memory_dat    = &(nr_mem->memory_dat);
  } else {
    q_remain      = NULL;
    current_node  = NULL;
    memory_dat    = NULL;
  }

#ifndef VRNA_NR_SAMPLING_HASH
  /* non-redundant data-structure memorization nodes */
  NR_NODE *memorized_node_prev  = NULL;           /* remembers previous-to-current node in linked list */
  NR_NODE *memorized_node_cur   = NULL;           /* remembers actual node in linked list */
#endif

  fbd   = 0.;                             /* stores weight of forbidden terms for given q[ij]*/
  fbds  = 0.;                             /* stores weight of forbidden term for given motif */

  n = vc->length;

  pf_params = vc->exp_params;
  md        = &(vc->exp_params->model_details);
  my_iindx  = vc->iindx;
  matrices  = vc->exp_matrices;

  hc = vc->hc;
  if (vc->type == VRNA_FC_TYPE_SINGLE) {
    n_seq = 1;
    S1    = vc->sequence_encoding;
    S2    = vc->sequence_encoding2;
    S     = NULL;
    S5    = NULL;
    S3    = NULL;
    a2s   = NULL;
  } else {
    n_seq = vc->n_seq;
    S1    = NULL;
    S2    = NULL;
    S     = vc->S;
    S5    = vc->S5;
    S3    = vc->S3;
    a2s   = vc->a2s;
  }

  hard_constraints  = hc->mx;
  hc_up_ext         = hc->up_ext;
  sc_wrapper_ext    = &(sc_wrap->sc_wrapper_ext);

  /* assume successful backtracing by default */
  ret = 1;

  q     = matrices->q;
  qb    = matrices->qb;
  q1k   = matrices->q1k;
  qln   = matrices->qln;
  scale = matrices->scale;

#ifndef VRNA_NR_SAMPLING_HASH
  if (current_node) {
    memorized_node_prev = NULL;
    memorized_node_cur  = (*current_node)->head;
  }

#endif

  q_temp = 0.;

#ifdef VRNA_WITH_BOUSTROPHEDON
  j = init_val;
  if (j > 1) {
    /* find j position of first pair */
    for (; j > 1; j--) {
      if (hc_up_ext[j]) {
        if (current_node) {
          fbd = NR_TOTAL_WEIGHT(*current_node) *
                q1k[j] /
                (*q_remain);

#ifdef  USE_FLOAT_PF
          if (fabsf(NR_TOTAL_WEIGHT(*current_node) - (*q_remain)) / (*q_remain) <= FLT_EPSILON)
#else
          if (fabs(NR_TOTAL_WEIGHT(*current_node) - (*q_remain)) / (*q_remain) <= DBL_EPSILON)
#endif
            /* exhausted ensemble */
            return 0;
        }

        r       = vrna_urn() * (q1k[j] - fbd);
        q_temp  = q1k[j - 1] * scale[1];

        if (sc_wrapper_ext->red_ext)
          q_temp *= sc_wrapper_ext->red_ext(1, j, 1, j - 1, sc_wrapper_ext);

        if (current_node) {
          fbds = NR_GET_WEIGHT(*current_node, memorized_node_cur, NRT_UNPAIRED_SG, j - 1, j) *
                 q1k[j] /
                 (*q_remain);
        }

        if (r > (q_temp - fbds)) {
          break;                /* j is paired */
        } else if (current_node) {
          /* j is unpaired */
          *q_remain *= q_temp / q1k[j];
#ifdef VRNA_NR_SAMPLING_HASH
          *current_node = add_if_nexists(NRT_UNPAIRED_SG, j - 1, j, *current_node, *q_remain);
#else
          *current_node = add_if_nexists_ll(memory_dat,
                                            NRT_UNPAIRED_SG,
                                            j - 1,
                                            j,
                                            memorized_node_prev,
                                            memorized_node_cur,
                                            *current_node,
                                            *q_remain);
          reset_cursor(&memorized_node_prev, &memorized_node_cur, *current_node); /* resets cursor */
#endif
        }
      }
    }
    if (j <= md->min_loop_size + 1)
      return ret;         /* no more pairs, but still successful */

#ifndef  VRNA_NR_SAMPLING_HASH
    if (current_node)
      advance_cursor(&memorized_node_prev, &memorized_node_cur, NRT_UNPAIRED_SG, j - 1, j);

#endif
    /* now find the pairing partner i */
    if (current_node) {
      fbd = NR_TOTAL_WEIGHT_TYPE(NRT_EXT_LOOP, *current_node) *
            q1k[j] /
            (*q_remain);
    }

    r = vrna_urn() * (q1k[j] - q_temp - fbd);
    u = j - 1;
    i = 2;

    for (qt = 0, k = 1; k < j; k++) {
      /* apply alternating boustrophedon scheme to variable i */
      i = (int)(1 + (u - 1) * ((k - 1) % 2)) +
          (int)((1 - (2 * ((k - 1) % 2))) * ((k - 1) / 2));
      ij            = my_iindx[i] - j;
      hc_decompose  = hard_constraints[n * j + i];
      if (hc_decompose & VRNA_CONSTRAINT_CONTEXT_EXT_LOOP) {
        qkl = qb[ij] *
              q1k[i - 1];

        if (vc->type == VRNA_FC_TYPE_SINGLE) {
          type  = vrna_get_ptype_md(S2[i], S2[j], md);
          qkl   *= vrna_exp_E_ext_stem(type,
                                       (i > 1) ? S1[i - 1] : -1,
                                       (j < n) ? S1[j + 1] : -1,
                                       pf_params);
        } else {
          for (s = 0; s < n_seq; s++) {
            type  = vrna_get_ptype_md(S[s][i], S[s][j], md);
            qkl   *= vrna_exp_E_ext_stem(type,
                                         (a2s[s][i] > 1) ? S5[s][i] : -1,
                                         (a2s[s][j] < a2s[s][n]) ? S3[s][j] : -1,
                                         pf_params);
          }
        }

        if ((sc_wrapper_ext->red_stem) && (i == 1))
          q_temp *= sc_wrapper_ext->red_stem(i, j, i, j, sc_wrapper_ext);
        else if ((sc_wrapper_ext->split) && (i > 1))
          q_temp *= sc_wrapper_ext->split(1, j, i, sc_wrapper_ext) *
                    sc_wrapper_ext->red_stem(i, j, i, j, sc_wrapper_ext);

        if (current_node) {
          fbds = NR_GET_WEIGHT(*current_node, memorized_node_cur, NRT_EXT_LOOP, i, j) *
                 q1k[j] /
                 (*q_remain);
          qt += qkl - fbds;
        } else {
          qt += qkl;
        }

        if (qt > r) {
          if (current_node) {
            *q_remain *= qkl / q1k[j];
#ifdef VRNA_NR_SAMPLING_HASH
            *current_node = add_if_nexists(NRT_EXT_LOOP, i, j, *current_node, *q_remain);
#else
            *current_node = add_if_nexists_ll(memory_dat,
                                              NRT_EXT_LOOP,
                                              i,
                                              j,
                                              memorized_node_prev,
                                              memorized_node_cur,
                                              *current_node,
                                              *q_remain);
#endif
          }

          break;           /* j is paired */
        }

#ifndef VRNA_NR_SAMPLING_HASH
        if (current_node)
          advance_cursor(&memorized_node_prev, &memorized_node_cur, NRT_EXT_LOOP, i, j);

#endif
      }
    }
    if (k == j) {
      if (current_node) {
        /* exhausted ensemble */
        return 0;
      } else {
        vrna_message_warning("backtracking failed in ext loop");
        /* error */
        return -1;
      }
    }

    backtrack(i, j, pstruc, vc, sc_wrap, nr_mem);
    j   = i - 1;
    ret = backtrack_ext_loop(j, pstruc, vc, length, sc_wrap, nr_mem);
  }

#else
  int start = init_val;
  if (start < length) {
    /* find i position of first pair */
    for (i = start; i < length; i++) {
      if (hc_up_ext[i]) {
        if (current_node) {
          fbd = NR_TOTAL_WEIGHT(*current_node) *
                qln[i] /
                (*q_remain);
        }

        r       = vrna_urn() * (qln[i] - fbd);
        q_temp  = qln[i + 1] * scale[1];

        if (sc_wrapper_ext->red_ext)
          q_temp *= sc_wrapper_ext->red_ext(i, length, i + 1, length, sc_wrapper_ext);

        if (current_node) {
          fbds = NR_GET_WEIGHT(*current_node, memorized_node_cur, NRT_UNPAIRED_SG, i, i + 1) *
                 qln[i] /
                 (*q_remain);
        }

        if (r > (q_temp - fbds)) {
          break;                /* i is paired */
        } else if (current_node) {
          *q_remain *= q_temp / qln[i];
#ifdef VRNA_NR_SAMPLING_HASH
          *current_node = add_if_nexists(NRT_UNPAIRED_SG, i, i + 1, *current_node, *q_remain);
#else
          *current_node = add_if_nexists_ll(memory_dat,
                                            NRT_UNPAIRED_SG,
                                            i,
                                            i + 1,
                                            memorized_node_prev,
                                            memorized_node_cur,
                                            *current_node,
                                            *q_remain);
          reset_cursor(&memorized_node_prev, &memorized_node_cur, *current_node); /* resets cursor */
#endif
        }
      }
    }
    if (i >= length)
      return ret;         /* no more pairs, but still successful */

#ifndef VRNA_NR_SAMPLING_HASH
    if (current_node)
      advance_cursor(&memorized_node_prev, &memorized_node_cur, NRT_UNPAIRED_SG, i, i + 1);

#endif

    /* now find the pairing partner j */
    if (current_node) {
      fbd = NR_TOTAL_WEIGHT_TYPE(NRT_EXT_LOOP, *current_node) *
            qln[i] /
            (*q_remain);
    }

    r = vrna_urn() * (qln[i] - q_temp - fbd);
    for (qt = 0, j = i + 1; j <= length; j++) {
      ij            = my_iindx[i] - j;
      hc_decompose  = hard_constraints[n * i + j];
      if (hc_decompose & VRNA_CONSTRAINT_CONTEXT_EXT_LOOP) {
        qkl = qb[ij];
        if (vc->type == VRNA_FC_TYPE_SINGLE) {
          type  = vrna_get_ptype_md(S2[i], S2[j], md);
          qkl   *= vrna_exp_E_ext_stem(type,
                                       (i > 1) ? S1[i - 1] : -1,
                                       (j < n) ? S1[j + 1] : -1,
                                       pf_params);
        } else {
          for (s = 0; s < n_seq; s++) {
            type  = vrna_get_ptype_md(S[s][i], S[s][j], md);
            qkl   *= vrna_exp_E_ext_stem(type,
                                         (a2s[s][i] > 1) ? S5[s][i] : -1,
                                         (a2s[s][j] < a2s[s][n]) ? S3[s][j] : -1,
                                         pf_params);
          }
        }

        if (j < length) {
          qkl *= qln[j + 1];
          if (sc_wrapper_ext->split)
            qkl *= sc_wrapper_ext->split(i, length, j + 1, sc_wrapper_ext) *
                   sc_wrapper_ext->red_stem(i, j, i, j, sc_wrapper_ext);
        } else if (sc_wrapper_ext->red_stem) {
          qkl *= sc_wrapper_ext->red_stem(i, j, i, j, sc_wrapper_ext);
        }

        if (current_node) {
          fbds = NR_GET_WEIGHT(*current_node, memorized_node_cur, NRT_EXT_LOOP, i, j) *
                 qln[i] /
                 (*q_remain);
          qt += qkl - fbds;
        } else {
          qt += qkl;
        }

        if (qt > r) {
          if (current_node) {
            *q_remain *= qkl / qln[i];
#ifdef VRNA_NR_SAMPLING_HASH
            *current_node = add_if_nexists(NRT_EXT_LOOP, i, j, *current_node, *q_remain);
#else
            *current_node = add_if_nexists_ll(memory_dat,
                                              NRT_EXT_LOOP,
                                              i,
                                              j,
                                              memorized_node_prev,
                                              memorized_node_cur,
                                              *current_node,
                                              *q_remain);
#endif
          }

          break;       /* j is paired */
        }

#ifndef VRNA_NR_SAMPLING_HASH
        if (current_node)
          advance_cursor(&memorized_node_prev, &memorized_node_cur, NRT_EXT_LOOP, i, j);

#endif
      }
    }
    if (j == length + 1) {
      if (current_node) {
        /* exhausted ensemble */
        return 0;
      } else {
        vrna_message_warning("backtracking failed in ext loop");
        /* error */
        return -1;
      }
    }

    start = j + 1;
    ret   = backtrack(i, j, pstruc, vc, sc_wrap, nr_mem);
    if (!ret)
      return ret;

    ret = backtrack_ext_loop(start, pstruc, vc, length, sc_wrap, nr_mem);
  }

#endif

  return ret;
}


/* non redundant version of function bactrack_qm */
PRIVATE int
backtrack_qm(int                              i,
             int                              j,
             char                             *pstruc,
             vrna_fold_compound_t             *vc,
             struct sc_wrappers               *sc_wrap,
             struct vrna_pbacktrack_memory_s  *nr_mem)
{
  /* divide multiloop into qm and qm1  */
  int                       k, u, cnt, span, turn, is_unpaired, *my_iindx, *jindx, *hc_up_ml, ret;
  FLT_OR_DBL                qmt, fbd, fbds, r, q_temp, *qm, *qm1, *expMLbase;
  double                    *q_remain;
  vrna_hc_t                 *hc;
  vrna_mx_pf_t              *matrices;

  struct sc_mb_exp_dat      *sc_wrapper_ml;
  struct nr_memory          **memory_dat;

  NR_NODE                   **current_node;

  if (nr_mem) {
    q_remain      = &(nr_mem->q_remain);
    current_node  = &(nr_mem->current_node);
    memory_dat    = &(nr_mem->memory_dat);
  } else {
    q_remain      = NULL;
    current_node  = NULL;
    memory_dat    = NULL;
  }

#ifndef VRNA_NR_SAMPLING_HASH
  /* non-redundant data-structure memorization nodes */
  NR_NODE *memorized_node_prev  = NULL;     /* remembers previous-to-current node in linked list */
  NR_NODE *memorized_node_cur   = NULL;     /* remembers actual node in linked list */
#endif

  ret   = 1;
  fbd   = 0.;                       /* stores weight of forbidden terms for given q[ij]*/
  fbds  = 0.;                       /* stores weight of forbidden term for given motif */

  is_unpaired = 0;

  matrices  = vc->exp_matrices;
  my_iindx  = vc->iindx;
  jindx     = vc->jindx;

  hc            = vc->hc;
  hc_up_ml      = hc->up_ml;
  sc_wrapper_ml = &(sc_wrap->sc_wrapper_ml);

  qm        = matrices->qm;
  qm1       = matrices->qm1;
  expMLbase = matrices->expMLbase;

  turn = vc->exp_params->model_details.min_loop_size;

#ifndef VRNA_NR_SAMPLING_HASH
  if (current_node) {
    memorized_node_prev = NULL;
    memorized_node_cur  = (*current_node)->head;
  }

#endif

  if (j > i) {
    /* now backtrack  [i ... j] in qm[] */

    if (current_node) {
      fbd = NR_TOTAL_WEIGHT(*current_node) *
            qm[my_iindx[i] - j] /
            (*q_remain);
    }

    r = vrna_urn() * (qm[my_iindx[i] - j] - fbd);
    if (current_node) {
      fbds = NR_GET_WEIGHT(*current_node, memorized_node_cur, NRT_QM_UNPAIR, i, 0) *
             qm[my_iindx[i] - j] /
             (*q_remain);
      qmt = qm1[jindx[j] + i] - fbds;
    } else {
      qmt = qm1[jindx[j] + i];
    }

    k       = cnt = i;
    q_temp  = qm1[jindx[j] + i];

    if (qmt < r) {
#ifndef VRNA_NR_SAMPLING_HASH
      if (current_node)
        advance_cursor(&memorized_node_prev, &memorized_node_cur, NRT_QM_UNPAIR, i, 0);

#endif

      for (span = j - i, cnt = i + 1; cnt <= j; cnt++) {
#ifdef VRNA_WITH_BOUSTROPHEDON
        k = (int)(i + 1 + span * ((cnt - i - 1) % 2)) +
            (int)((1 - (2 * ((cnt - i - 1) % 2))) * ((cnt - i) / 2));
#else
        k = cnt;
#endif
        q_temp  = 0.;
        u       = k - i;
        /* [i...k] is unpaired */
        if (hc_up_ml[i] >= u) {
          q_temp += expMLbase[u] * qm1[jindx[j] + k];

          if (sc_wrapper_ml->red_ml)
            q_temp *= sc_wrapper_ml->red_ml(i, j, k, j, sc_wrapper_ml);

          if (current_node) {
            fbds = NR_GET_WEIGHT(*current_node, memorized_node_cur, NRT_QM_UNPAIR, k, 0) *
                   qm[my_iindx[i] - j] /
                   (*q_remain);
            qmt += q_temp - fbds;
          } else {
            qmt += q_temp;
          }
        }

        if (qmt >= r) {
          /* we have chosen unpaired version */
          is_unpaired = 1;
          break;
        }

#ifndef VRNA_NR_SAMPLING_HASH
        if (current_node)
          advance_cursor(&memorized_node_prev, &memorized_node_cur, NRT_QM_UNPAIR, k, 0);

#endif

        /* split between k-1, k */
        q_temp = qm[my_iindx[i] - (k - 1)] *
                 qm1[jindx[j] + k];

        if (sc_wrapper_ml->decomp_ml)
          q_temp *= sc_wrapper_ml->decomp_ml(i, j, k - 1, k, sc_wrapper_ml);

        if (current_node) {
          fbds = NR_GET_WEIGHT(*current_node, memorized_node_cur, NRT_QM_PAIR, k, 0) *
                 qm[my_iindx[i] - j] /
                 (*q_remain);
          qmt += q_temp - fbds;
        } else {
          qmt += q_temp;
        }

        if (qmt >= r)
          break;

#ifndef VRNA_NR_SAMPLING_HASH
        if (current_node)
          advance_cursor(&memorized_node_prev, &memorized_node_cur, NRT_QM_PAIR, k, 0);

#endif
      }
    } else {
      is_unpaired = 1;
    }

    if (current_node) {
      *q_remain *= q_temp / qm[my_iindx[i] - j];
#ifdef VRNA_NR_SAMPLING_HASH
      if (is_unpaired)
        *current_node = add_if_nexists(NRT_QM_UNPAIR, k, 0, *current_node, *q_remain);
      else
        *current_node = add_if_nexists(NRT_QM_PAIR, k, 0, *current_node, *q_remain);

#else
      if (is_unpaired)
        *current_node = add_if_nexists_ll(memory_dat,
                                          NRT_QM_UNPAIR,
                                          k,
                                          0,
                                          memorized_node_prev,
                                          memorized_node_cur,
                                          *current_node,
                                          *q_remain);
      else
        *current_node = add_if_nexists_ll(memory_dat,
                                          NRT_QM_PAIR,
                                          k,
                                          0,
                                          memorized_node_prev,
                                          memorized_node_cur,
                                          *current_node,
                                          *q_remain);

#endif
    }

    if (cnt > j)
      return 0;

    ret = backtrack_qm1(k, j, pstruc, vc, sc_wrap, nr_mem);

    if (ret == 0)
      return ret;

    if (k < i + turn)
      return ret;         /* no more pairs */

    if (!is_unpaired) {
      /* if we've chosen creating a branch in [i..k-1] */
      ret = backtrack_qm(i, k - 1, pstruc, vc, sc_wrap, nr_mem);

      if (ret == 0)
        return ret;
    }
  }

  return ret;
}


PRIVATE int
backtrack_qm1(int                             i,
              int                             j,
              char                            *pstruc,
              vrna_fold_compound_t            *vc,
              struct sc_wrappers              *sc_wrap,
              struct vrna_pbacktrack_memory_s *nr_mem)
{
  /* i is paired to l, i<l<j; backtrack in qm1 to find l */
  unsigned char             *hard_constraints;
  char                      *ptype;
  short                     *S1, **S, **S5, **S3;
  unsigned int              n, s, n_seq;
  int                       ii, l, il, type, turn, u, *my_iindx, *jindx, *hc_up_ml;
  FLT_OR_DBL                qt, fbd, fbds, r, q_temp, *qm1, *qb, *expMLbase;
  double                    *q_remain;
  vrna_exp_param_t          *pf_params;
  vrna_md_t                 *md;
  vrna_hc_t                 *hc;
  vrna_mx_pf_t              *matrices;

  struct nr_memory          **memory_dat;
  struct sc_mb_exp_dat      *sc_wrapper_ml;

  NR_NODE                   **current_node;

  if (nr_mem) {
    q_remain      = &(nr_mem->q_remain);
    current_node  = &(nr_mem->current_node);
    memory_dat    = &(nr_mem->memory_dat);
  } else {
    q_remain      = NULL;
    current_node  = NULL;
    memory_dat    = NULL;
  }

#ifndef VRNA_NR_SAMPLING_HASH
  /* non-redundant data-structure memorization nodes */
  NR_NODE *memorized_node_prev  = NULL;           /* remembers previous-to-current node in linked list */
  NR_NODE *memorized_node_cur   = NULL;           /* remembers actual node in linked list */
#endif

  n                 = vc->length;
  fbd               = 0.;
  fbds              = 0.;
  pf_params         = vc->exp_params;
  md                = &(pf_params->model_details);
  my_iindx          = vc->iindx;
  jindx             = vc->jindx;
  hc                = vc->hc;
  hc_up_ml          = hc->up_ml;
  hard_constraints  = hc->mx;
  sc_wrapper_ml     = &(sc_wrap->sc_wrapper_ml);

  matrices  = vc->exp_matrices;
  qb        = matrices->qb;
  qm1       = matrices->qm1;
  expMLbase = matrices->expMLbase;
  if (vc->type == VRNA_FC_TYPE_SINGLE) {
    n_seq = 1;
    ptype = vc->ptype;
    S1    = vc->sequence_encoding;
    S     = NULL;
    S5    = NULL;
    S3    = NULL;
  } else {
    n_seq = vc->n_seq;
    ptype = NULL;
    S1    = NULL;
    S     = vc->S;
    S5    = vc->S5;
    S3    = vc->S3;
  }

  turn = pf_params->model_details.min_loop_size;

#ifndef VRNA_NR_SAMPLING_HASH
  if (current_node) {
    memorized_node_prev = NULL;
    memorized_node_cur  = (*current_node)->head;
  }

#endif

  if (current_node) {
    fbd = NR_TOTAL_WEIGHT(*current_node) *
          qm1[jindx[j] + i] /
          (*q_remain);
  }

  r   = vrna_urn() * (qm1[jindx[j] + i] - fbd);
  ii  = my_iindx[i];
  for (qt = 0., l = j; l > i + turn; l--) {
    il = jindx[l] + i;
    if (hard_constraints[n * i + l] & VRNA_CONSTRAINT_CONTEXT_MB_LOOP_ENC) {
      u = j - l;
      if (hc_up_ml[l + 1] >= u) {
        q_temp = qb[ii - l] *
                 expMLbase[j - l];

        if (vc->type == VRNA_FC_TYPE_SINGLE) {
          type    = vrna_get_ptype(il, ptype);
          q_temp  *= exp_E_MLstem(type, S1[i - 1], S1[l + 1], pf_params);
        } else {
          for (s = 0; s < n_seq; s++) {
            type    = vrna_get_ptype_md(S[s][i], S[s][l], md);
            q_temp  *= exp_E_MLstem(type, S5[s][i], S3[s][l], pf_params);
          }
        }

        if (sc_wrapper_ml->red_stem)
          q_temp *= sc_wrapper_ml->red_stem(i, j, i, l, sc_wrapper_ml);

        if (current_node) {
          fbds = NR_GET_WEIGHT(*current_node, memorized_node_cur, NRT_QM1_BRANCH, i, l) *
                 qm1[jindx[j] + i] /
                 (*q_remain);
          qt += q_temp - fbds;
        } else {
          qt += q_temp;
        }

        if (qt >= r) {
          if (current_node) {
            *q_remain *= q_temp / qm1[jindx[j] + i];
#ifdef VRNA_NR_SAMPLING_HASH
            *current_node = add_if_nexists(NRT_QM1_BRANCH, i, l, *current_node, *q_remain);
#else
            *current_node = add_if_nexists_ll(memory_dat,
                                              NRT_QM1_BRANCH,
                                              i,
                                              l,
                                              memorized_node_prev,
                                              memorized_node_cur,
                                              *current_node,
                                              *q_remain);
#endif
          }

          break;
        }

#ifndef VRNA_NR_SAMPLING_HASH
        if (current_node)
          advance_cursor(&memorized_node_prev, &memorized_node_cur, NRT_QM1_BRANCH, i, l);

#endif
      } else {
        l = i + turn;
        break;
      }
    }
  }
  if (l < i + turn + 1) {
    if (current_node) {
      return 0;
    } else {
      vrna_message_error("backtrack failed in qm1");
      return 0;
    }
  }

  return backtrack(i, l, pstruc, vc, sc_wrap, nr_mem);
}


PRIVATE void
backtrack_qm2(int                   k,
              int                   n,
              char                  *pstruc,
              vrna_fold_compound_t  *vc,
              struct sc_wrappers    *sc_wrap)
{
  int                       u, turn, *jindx;
  FLT_OR_DBL                qom2t, r, *qm1, *qm2;
  struct sc_mb_exp_dat      *sc_wrapper_ml;

  jindx         = vc->jindx;
  qm1           = vc->exp_matrices->qm1;
  qm2           = vc->exp_matrices->qm2;
  turn          = vc->exp_params->model_details.min_loop_size;
  sc_wrapper_ml = &(sc_wrap->sc_wrapper_ml);

  r = vrna_urn() * qm2[k];
  /* we have to search for our barrier u between qm1 and qm1  */
  if (sc_wrapper_ml->decomp_ml) {
    for (qom2t = 0., u = k + turn + 1; u < n - turn - 1; u++) {
      qom2t += qm1[jindx[u] + k] *
               qm1[jindx[n] + (u + 1)] *
               sc_wrapper_ml->decomp_ml(k, n, u, u + 1, sc_wrapper_ml);

      if (qom2t > r)
        break;
    }
  } else {
    for (qom2t = 0., u = k + turn + 1; u < n - turn - 1; u++) {
      qom2t += qm1[jindx[u] + k] * qm1[jindx[n] + (u + 1)];
      if (qom2t > r)
        break;
    }
  }

  if (u == n - turn)
    vrna_message_error("backtrack failed in qm2");

  backtrack_qm1(k, u, pstruc, vc, sc_wrap, NULL);
  backtrack_qm1(u + 1, n, pstruc, vc, sc_wrap, NULL);
}


PRIVATE int
backtrack(int                             i,
          int                             j,
          char                            *pstruc,
          vrna_fold_compound_t            *vc,
          struct sc_wrappers              *sc_wrap,
          struct vrna_pbacktrack_memory_s *nr_mem)
{
  unsigned char             *hard_constraints, hc_decompose;
  char                      *ptype;
  short                     *S1, **S, **S5, **S3;
  unsigned int              **a2s, s, n_seq, n, type, type_2, *types, u1_local, u2_local;
  int                       *my_iindx, *jindx, *hc_up_int, ret, *pscore, turn, *rtype,
                            k, l, kl, u1, u2, max_k, min_l, ii, jj;
  FLT_OR_DBL                *qb, *qm, *qm1, *scale, r, fbd, fbds, qbt1, qbr, qt, q_temp,
                            kTn, closingPair, expMLclosing;
  double                    *q_remain;
  vrna_mx_pf_t              *matrices;
  vrna_exp_param_t          *pf_params;
  vrna_md_t                 *md;
  vrna_hc_t                 *hc;

  struct nr_memory          **memory_dat;
  struct sc_int_exp_dat     *sc_wrapper_int;
  struct sc_mb_exp_dat      *sc_wrapper_ml;

  NR_NODE                   **current_node;

  if (nr_mem) {
    q_remain      = &(nr_mem->q_remain);
    current_node  = &(nr_mem->current_node);
    memory_dat    = &(nr_mem->memory_dat);
  } else {
    q_remain      = NULL;
    current_node  = NULL;
    memory_dat    = NULL;
  }

#ifndef VRNA_NR_SAMPLING_HASH
  /* non-redundant data-structure memorization nodes */
  NR_NODE *memorized_node_prev  = NULL;           /* remembers previous-to-current node in linked list */
  NR_NODE *memorized_node_cur   = NULL;           /* remembers actual node in linked list */
#endif

  ret     = 1;                            /* default is success */
  fbd     = 0.;                           /* stores weight of forbidden terms for given q[ij] */
  fbds    = 0.;                           /* stores weight of forbidden term for given motif */
  qbt1    = 0.;
  q_temp  = 0.;

  n         = vc->length;
  pf_params = vc->exp_params;
  kTn       = pf_params->kT / 10.;
  md        = &(pf_params->model_details);
  my_iindx  = vc->iindx;
  jindx     = vc->jindx;
  turn      = pf_params->model_details.min_loop_size;
  rtype     = &(pf_params->model_details.rtype[0]);

  if (vc->type == VRNA_FC_TYPE_SINGLE) {
    n_seq         = 1;
    ptype         = vc->ptype;
    types         = NULL;
    pscore        = NULL;
    S1            = vc->sequence_encoding;
    S             = NULL;
    S5            = NULL;
    S3            = NULL;
    a2s           = NULL;
    expMLclosing  = pf_params->expMLclosing;
  } else {
    n_seq         = vc->n_seq;
    ptype         = NULL;
    types         = (unsigned int *)vrna_alloc(sizeof(unsigned int) * n_seq);
    pscore        = vc->pscore;
    S1            = NULL;
    S             = vc->S;
    S5            = vc->S5;
    S3            = vc->S3;
    a2s           = vc->a2s;
    expMLclosing  = pow(pf_params->expMLclosing, (double)n_seq);
  }

  hc                = vc->hc;
  hc_up_int         = hc->up_int;
  hard_constraints  = hc->mx;
  sc_wrapper_int    = &(sc_wrap->sc_wrapper_int);
  sc_wrapper_ml     = &(sc_wrap->sc_wrapper_ml);

  matrices  = vc->exp_matrices;
  qb        = matrices->qb;
  qm        = matrices->qm;
  qm1       = matrices->qm1;
  scale     = matrices->scale;

#ifndef VRNA_NR_SAMPLING_HASH
  if (current_node) {
    memorized_node_prev = NULL;
    memorized_node_cur  = (*current_node)->head;
  }

#endif

  hc_decompose = hard_constraints[n * j + i];

  do {
    k = i;
    l = j;

    qbr = qb[my_iindx[i] - j];

    if (vc->type == VRNA_FC_TYPE_SINGLE) {
      type = vrna_get_ptype(jindx[j] + i, ptype);
    } else {
      qbr /= exp(pscore[jindx[j] + i] / kTn);
      for (s = 0; s < n_seq; s++)
        types[s] = vrna_get_ptype_md(S[s][i], S[s][j], md);
    }

    if (current_node)
      fbd = NR_TOTAL_WEIGHT(*current_node) * qbr / (*q_remain);

    pstruc[i - 1] = '(';
    pstruc[j - 1] = ')';

    r     = vrna_urn() * (qbr - fbd);
    qbt1  = 0.;

    hc_decompose = hard_constraints[n * i + j];

    /* hairpin contribution */
    q_temp = vrna_exp_E_hp_loop(vc, i, j);

    if (current_node) {
      fbds = NR_GET_WEIGHT(*current_node, memorized_node_cur, NRT_HAIRPIN, 0, 0) *
             qbr /
             (*q_remain);
      qbt1 += (q_temp - fbds);
    } else {
      qbt1 += q_temp;
    }

    if (qbt1 >= r) {
      /* found the hairpin we're done */
      if (current_node) {
        *q_remain *= q_temp / qbr;
#ifdef VRNA_NR_SAMPLING_HASH
        *current_node = add_if_nexists(NRT_HAIRPIN, 0, 0, *current_node, *q_remain);
#else
        *current_node = add_if_nexists_ll(memory_dat,
                                          NRT_HAIRPIN,
                                          0,
                                          0,
                                          memorized_node_prev,
                                          memorized_node_cur,
                                          *current_node,
                                          *q_remain);
#endif
      }

      free(types);

      return ret;
    }

#ifndef VRNA_NR_SAMPLING_HASH
    if (current_node)
      advance_cursor(&memorized_node_prev, &memorized_node_cur, NRT_HAIRPIN, 0, 0);

#endif

    if (hc_decompose & VRNA_CONSTRAINT_CONTEXT_INT_LOOP) {
      /* interior loop contributions */
      max_k = i + MAXLOOP + 1;
      max_k = MIN2(max_k, j - turn - 2);
      max_k = MIN2(max_k, i + 1 + hc_up_int[i + 1]);
      for (k = i + 1; k <= max_k; k++) {
        u1    = k - i - 1;
        min_l = MAX2(k + turn + 1, j - 1 - MAXLOOP + u1);
        kl    = my_iindx[k] - j + 1;
        for (u2 = 0, l = j - 1; l >= min_l; l--, kl++, u2++) {
          if (hc_up_int[l + 1] < u2)
            break;

          if (hard_constraints[n * k + l] & VRNA_CONSTRAINT_CONTEXT_INT_LOOP_ENC) {
            q_temp = qb[kl]
                     * scale[u1 + u2 + 2];

            if (vc->type == VRNA_FC_TYPE_SINGLE) {
              type_2 = rtype[vrna_get_ptype(jindx[l] + k, ptype)];

              /* add *scale[u1+u2+2] */
              q_temp *= exp_E_IntLoop(u1,
                                      u2,
                                      type,
                                      type_2,
                                      S1[i + 1],
                                      S1[j - 1],
                                      S1[k - 1],
                                      S1[l + 1],
                                      pf_params);
            } else {
              for (s = 0; s < n_seq; s++) {
                u1_local  = a2s[s][k - 1] - a2s[s][i] /*??*/;
                u2_local  = a2s[s][j - 1] - a2s[s][l];
                type_2    = vrna_get_ptype_md(S[s][l], S[s][k], md);
                q_temp    *= exp_E_IntLoop(u1_local,
                                           u2_local,
                                           types[s],
                                           type_2,
                                           S3[s][i],
                                           S5[s][j],
                                           S5[s][k],
                                           S3[s][l],
                                           pf_params);
              }
            }

            if (sc_wrapper_int->pair)
              q_temp *= sc_wrapper_int->pair(i, j, k, l, sc_wrapper_int);

            if (current_node) {
              fbds = NR_GET_WEIGHT(*current_node, memorized_node_cur, NRT_IT_LOOP, k, l) *
                     qbr /
                     (*q_remain);
              qbt1 += q_temp - fbds;
            } else {
              qbt1 += q_temp;
            }

            if (qbt1 >= r)
              break;

#ifndef VRNA_NR_SAMPLING_HASH
            if (current_node)
              advance_cursor(&memorized_node_prev, &memorized_node_cur, NRT_IT_LOOP, k, l);

#endif
          }
        }
        if (qbt1 >= r)
          break;
      }
      if (k <= max_k) {
        if (current_node) {
          *q_remain *= q_temp / qbr;
#ifdef VRNA_NR_SAMPLING_HASH
          *current_node = add_if_nexists(NRT_IT_LOOP, k, l, *current_node, *q_remain);
#else
          *current_node = add_if_nexists_ll(memory_dat,
                                            NRT_IT_LOOP,
                                            k,
                                            l,
                                            memorized_node_prev,
                                            memorized_node_cur,
                                            *current_node,
                                            *q_remain);
#endif
        }

        free(types);

        return backtrack(k, l, pstruc, vc, sc_wrap, nr_mem); /* found the interior loop, repeat for inside */
      } else {
        /* interior loop contributions did not exceed threshold, so we break */
        break;
      }
    } else {
      /* must not be interior loop, so we break out */
      break;
    }
  } while (1);

  /* backtrack in multi-loop */
  if (hard_constraints[n * j + i] & VRNA_CONSTRAINT_CONTEXT_MB_LOOP) {
    closingPair = expMLclosing *
                  scale[2];

    if (vc->type == VRNA_FC_TYPE_SINGLE) {
      type        = rtype[vrna_get_ptype(jindx[j] + i, ptype)];
      closingPair *= exp_E_MLstem(type, S1[j - 1], S1[i + 1], pf_params);
    } else {
      for (s = 0; s < n_seq; s++) {
        type        = vrna_get_ptype_md(S[s][j], S[s][i], md);
        closingPair *= exp_E_MLstem(type, S5[s][j], S3[s][i], pf_params);
      }
    }

    if (sc_wrapper_ml->pair)
      closingPair *= sc_wrapper_ml->pair(i, j, sc_wrapper_ml);

    i++;
    j--;
    /* find the first split index */
    ii  = my_iindx[i];  /* ii-j=[i,j] */
    jj  = jindx[j];     /* jj+i=[j,i] */

    if (sc_wrapper_ml->decomp_ml) {
      for (qt = qbt1, k = i + 1; k < j; k++) {
        q_temp = qm[ii - (k - 1)] *
                 qm1[jj + k] *
                 closingPair *
                 sc_wrapper_ml->decomp_ml(i,
                                          j,
                                          k - 1,
                                          k,
                                          sc_wrapper_ml);


        if (current_node) {
          fbds = NR_GET_WEIGHT(*current_node, memorized_node_cur, NRT_MT_LOOP, k, 0) *
                 qbr /
                 (*q_remain);
          qt    += q_temp - fbds;
          qbt1  += q_temp - fbds;
        } else {
          qt    += q_temp;
          qbt1  += q_temp;
        }

        if (qt >= r)
          break;

#ifndef VRNA_NR_SAMPLING_HASH
        if (current_node)
          advance_cursor(&memorized_node_prev, &memorized_node_cur, NRT_MT_LOOP, k, 0);

#endif
      }
    } else {
      for (qt = qbt1, k = i + 1; k < j; k++) {
        q_temp = qm[ii - (k - 1)] *
                 qm1[jj + k] *
                 closingPair;

        if (current_node) {
          fbds = NR_GET_WEIGHT(*current_node, memorized_node_cur, NRT_MT_LOOP, k, 0) *
                 qbr /
                 (*q_remain);
          qt    += q_temp - fbds;
          qbt1  += q_temp - fbds;
        } else {
          qt    += q_temp;
          qbt1  += q_temp;
        }

        if (qt >= r)
          break;

#ifndef VRNA_NR_SAMPLING_HASH
        if (current_node)
          advance_cursor(&memorized_node_prev, &memorized_node_cur, NRT_MT_LOOP, k, 0);

#endif
      }
    }

    if (k >= j) {
      if (current_node) {
        free(types);
        return 0; /* backtrack failed for non-redundant mode most likely due to numerical instabilities */
      } else {
        vrna_message_error("backtrack failed, can't find split index ");
      }
    }

    if (current_node) {
      *q_remain *= q_temp / qbr;
#ifdef VRNA_NR_SAMPLING_HASH
      *current_node = add_if_nexists(NRT_MT_LOOP, k, 0, *current_node, *q_remain);
#else
      *current_node = add_if_nexists_ll(memory_dat,
                                        NRT_MT_LOOP,
                                        k,
                                        0,
                                        memorized_node_prev,
                                        memorized_node_cur,
                                        *current_node,
                                        *q_remain);
#endif
    }

    ret = backtrack_qm1(k, j, pstruc, vc, sc_wrap, nr_mem);

    if (ret == 0) {
      free(types);
      return ret;
    }

    j = k - 1;

    ret = backtrack_qm(i, j, pstruc, vc, sc_wrap, nr_mem);
  }

  free(types);

  return ret;
}


PRIVATE unsigned int
pbacktrack_circ(vrna_fold_compound_t              *vc,
                unsigned int                      num_samples,
                vrna_boltzmann_sampling_callback  *bs_cb,
                void                              *data)
{
  unsigned char             *hc_mx, eval_loop;
  char                      *pstruc;
  short                     *S1, *S2, **S, **S5, **S3;
  unsigned int              type, type2, *tt, s, n_seq, **a2s, u1_local,
                            u2_local, u3_local, count;
  int                       i, j, k, l, n, u, *hc_up, *my_iindx, turn,
                            ln1, ln2, ln3, lstart;
  FLT_OR_DBL                r, qt, q_temp, qo, qmo, *scale, *qb, *qm, *qm2,
                            qb_ij, expMLclosing;
  vrna_exp_param_t          *pf_params;
  vrna_md_t                 *md;
  vrna_mx_pf_t              *matrices;
  struct sc_wrappers        *sc_wrap;
  struct sc_ext_exp_dat     *sc_wrapper_ext;
  struct sc_int_exp_dat     *sc_wrapper_int;
  struct sc_mb_exp_dat      *sc_wrapper_ml;

  n             = vc->length;
  pf_params     = vc->exp_params;
  md            = &(pf_params->model_details);
  matrices      = vc->exp_matrices;
  my_iindx      = vc->iindx;
  expMLclosing  = pf_params->expMLclosing;
  turn          = pf_params->model_details.min_loop_size;

  qo    = matrices->qo;
  qmo   = matrices->qmo;
  qb    = matrices->qb;
  qm    = matrices->qm;
  qm2   = matrices->qm2;
  scale = matrices->scale;

  hc_mx = vc->hc->mx;
  hc_up = vc->hc->up_int;

  sc_wrap         = sc_init(vc);
  sc_wrapper_ext  = &(sc_wrap->sc_wrapper_ext);
  sc_wrapper_int  = &(sc_wrap->sc_wrapper_int);
  sc_wrapper_ml   = &(sc_wrap->sc_wrapper_ml);

  if (vc->type == VRNA_FC_TYPE_SINGLE) {
    n_seq         = 1;
    tt            = NULL;
    S1            = vc->sequence_encoding;
    S2            = vc->sequence_encoding2;
    S             = NULL;
    S5            = NULL;
    S3            = NULL;
    a2s           = NULL;
    expMLclosing  = pf_params->expMLclosing;
  } else {
    n_seq         = vc->n_seq;
    tt            = (unsigned int *)vrna_alloc(sizeof(unsigned int) * n_seq);
    S1            = NULL;
    S2            = NULL;
    S             = vc->S;
    S5            = vc->S5;
    S3            = vc->S3;
    a2s           = vc->a2s;
    expMLclosing  = pow(pf_params->expMLclosing, (double)n_seq);
  }

  for (count = 0; count < num_samples; count++) {
    pstruc = vrna_alloc((n + 1) * sizeof(char));

    /* initialize pstruct with single bases  */
    memset(pstruc, '.', sizeof(char) * n);

    qt = 1.0 * scale[n];

    /* add soft constraints for open chain configuration */
    if (sc_wrapper_ext->red_up)
      qt *= sc_wrapper_ext->red_up(1, n, sc_wrapper_ext);

    r = vrna_urn() * qo;

    /* open chain? */
    if (qt > r)
      goto pbacktrack_circ_loop_end;

    for (i = 1; (i < n); i++) {
      for (j = i + turn + 1; (j <= n); j++) {
        u = n - j + i - 1;

        if (u < turn)
          continue;

        qb_ij = qb[my_iindx[i] - j];

        qt += qb_ij *
              vrna_exp_E_hp_loop(vc, j, i);

        /* found a hairpin? so backtrack in the enclosed part and we're done  */
        if (qt > r) {
          backtrack(i, j, pstruc, vc, sc_wrap, NULL);
          goto pbacktrack_circ_loop_end;
        }

        /* 2. search for (k,l) with which we can close an interior loop  */
        if (hc_mx[n * i + j] & VRNA_CONSTRAINT_CONTEXT_INT_LOOP) {
          if (vc->type == VRNA_FC_TYPE_SINGLE)
            type = vrna_get_ptype_md(S2[j], S2[i], md);
          else
            for (s = 0; s < n_seq; s++)
              tt[s] = vrna_get_ptype_md(S[s][j], S[s][i], md);

          for (k = j + 1; (k < n); k++) {
            ln1 = k - j - 1;
            if (ln1 + i - 1 > MAXLOOP)
              break;

            if (hc_up[j + 1] < ln1)
              break;

            lstart = ln1 + i - 1 + n - MAXLOOP;
            if (lstart < k + turn + 1)
              lstart = k + turn + 1;

            for (l = lstart; (l <= n); l++) {
              ln2 = (i - 1);
              ln3 = (n - l);

              if (hc_up[l + 1] < (ln2 + ln3))
                continue;

              if ((ln1 + ln2 + ln3) > MAXLOOP)
                continue;

              eval_loop = hc_mx[n * k + l] & VRNA_CONSTRAINT_CONTEXT_INT_LOOP;

              if (eval_loop) {
                q_temp = qb_ij *
                         qb[my_iindx[k] - l] *
                         scale[ln1 + ln2 + ln3];

                switch (vc->type) {
                  case VRNA_FC_TYPE_SINGLE:
                    type2   = vrna_get_ptype_md(S2[l], S2[k], md);
                    q_temp  *= exp_E_IntLoop(ln2 + ln3,
                                             ln1,
                                             type2,
                                             type,
                                             S1[l + 1],
                                             S1[k - 1],
                                             S1[i - 1],
                                             S1[j + 1],
                                             pf_params);
                    break;
                  case VRNA_FC_TYPE_COMPARATIVE:
                    for (s = 0; s < n_seq; s++) {
                      type2     = vrna_get_ptype_md(S[s][l], S[s][k], md);
                      u1_local  = a2s[s][i - 1];
                      u2_local  = a2s[s][k - 1] - a2s[s][j];
                      u3_local  = a2s[s][n] - a2s[s][l];
                      q_temp    *= exp_E_IntLoop(u1_local + u3_local,
                                                 u2_local,
                                                 type2,
                                                 tt[s],
                                                 S3[s][l],
                                                 S5[s][k],
                                                 S5[s][i],
                                                 S3[s][j],
                                                 pf_params);
                    }
                    break;
                }

                if (sc_wrapper_int->pair_ext)
                  q_temp *= sc_wrapper_int->pair_ext(i, j, k, l, sc_wrapper_int);

                qt += q_temp;
                /*
                 * found an exterior interior loop? also this time, we can go straight
                 * forward and backtracking the both enclosed parts and we're done
                 */
                if (qt > r) {
                  backtrack(i, j, pstruc, vc, sc_wrap, NULL);
                  backtrack(k, l, pstruc, vc, sc_wrap, NULL);
                  goto pbacktrack_circ_loop_end;
                }
              }
            }
          } /* end of kl double loop */
        }
      }
    }     /* end of ij double loop  */
    {
      /* as we reach this part, we have to search for our barrier between qm and qm2  */
      qt  = 0.;
      r   = vrna_urn() * qmo;
      if (sc_wrapper_ml->decomp_ml) {
        for (k = turn + 2; k < n - 2 * turn - 3; k++) {
          qt += qm[my_iindx[1] - k] *
                qm2[k + 1] *
                expMLclosing *
                sc_wrapper_ml->decomp_ml(1,
                                         n,
                                         k,
                                         k + 1,
                                         sc_wrapper_ml);


          /* backtrack in qm and qm2 if we've found a valid barrier k  */
          if (qt > r) {
            backtrack_qm(1, k, pstruc, vc, sc_wrap, NULL);
            backtrack_qm2(k + 1, n, pstruc, vc, sc_wrap);
            goto pbacktrack_circ_loop_end;
          }
        }
      } else {
        for (k = turn + 2; k < n - 2 * turn - 3; k++) {
          qt += qm[my_iindx[1] - k] *
                qm2[k + 1] *
                expMLclosing;
          /* backtrack in qm and qm2 if we've found a valid barrier k  */
          if (qt > r) {
            backtrack_qm(1, k, pstruc, vc, sc_wrap, NULL);
            backtrack_qm2(k + 1, n, pstruc, vc, sc_wrap);
            goto pbacktrack_circ_loop_end;
          }
        }
      }
    }
    /*
     * if we reach the actual end of this function, an error has occured
     * cause we HAVE TO find an exterior loop or an open chain!!!
     */
    vrna_message_error("backtracking failed in exterior loop");

pbacktrack_circ_loop_end:

    if (bs_cb)
      bs_cb(pstruc, data);

    free(pstruc);
  }

  sc_free(sc_wrap);

  return count;
}
