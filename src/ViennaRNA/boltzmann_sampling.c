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
  struct sc_wrapper_exp_ext sc_wrapper_ext;
  struct sc_wrapper_exp_ml  sc_wrapper_ml;
};

/*
 * In the following:
 * - q_remain is a pointer to value of sum of Boltzmann factors of still accessible solutions at that point
 * - current_node is a pointer to current node in datastructure memorizing the solutions and paths taken
 */
struct vrna_nr_memory_s {
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


PRIVATE char  *info_no_comparative =
  "No implementation for comparative structure prediction available.";


/*
 #################################
 # PRIVATE FUNCTION DECLARATIONS #
 #################################
 */

PRIVATE struct vrna_nr_memory_s *
nr_init(vrna_fold_compound_t *fc);


PRIVATE struct sc_wrappers *
sc_init(vrna_fold_compound_t *fc);


PRIVATE void
sc_free(struct sc_wrappers *sc_wrap);


PRIVATE char *
pbacktrack5_gen(vrna_fold_compound_t    *vc,
                int                     length,
                struct vrna_nr_memory_s *nr_mem);


PRIVATE int
backtrack(int                     i,
          int                     j,
          char                    *pstruc,
          vrna_fold_compound_t    *vc,
          struct sc_wrappers      *sc_wrap,
          struct vrna_nr_memory_s *nr_mem);


PRIVATE int
backtrack_ext_loop(int                      init_val,
                   char                     *pstruc,
                   vrna_fold_compound_t     *vc,
                   int                      length,
                   struct sc_wrappers       *sc_wrap,
                   struct vrna_nr_memory_s  *nr_mem);


PRIVATE int
backtrack_qm(int                      i,
             int                      j,
             char                     *pstruc,
             vrna_fold_compound_t     *vc,
             struct sc_wrappers       *sc_wrap,
             struct vrna_nr_memory_s  *nr_mem);


PRIVATE int
backtrack_qm1(int                     i,
              int                     j,
              char                    *pstruc,
              vrna_fold_compound_t    *vc,
              struct sc_wrappers      *sc_wrap,
              struct vrna_nr_memory_s *nr_mem);


PRIVATE void
backtrack_qm2(int                   u,
              int                   n,
              char                  *pstruc,
              vrna_fold_compound_t  *vc,
              struct sc_wrappers    *sc_wrap);


PRIVATE char *
wrap_pbacktrack_circ(vrna_fold_compound_t *vc);


/*
 #################################
 # BEGIN OF FUNCTION DEFINITIONS #
 #################################
 */
PUBLIC char *
vrna_pbacktrack5(vrna_fold_compound_t *fc,
                 unsigned int         length)
{
  if (fc) {
    vrna_mx_pf_t *matrices = fc->exp_matrices;

    if (length > fc->length) {
      vrna_message_warning("vrna_pbacktrack5: length exceeds sequence length");
    } else if (length == 0) {
      vrna_message_warning("vrna_pbacktrack5: length too small");
    } else if ((!matrices) || (!matrices->q) || (!matrices->qb) || (!matrices->qm) ||
               (!fc->exp_params)) {
      vrna_message_warning("vrna_pbacktrack5: %s", info_call_pf);
    } else if ((!fc->exp_params->model_details.uniq_ML) || (!matrices->qm1)) {
      vrna_message_warning("vrna_pbacktrack5: %s", info_set_uniq_ml);
    } else if ((fc->exp_params->model_details.circ) && (length < fc->length)) {
      vrna_message_warning("vrna_pbacktrack5: %s", info_no_circ);
    } else if (fc->type == VRNA_FC_TYPE_COMPARATIVE) {
      if (length < fc->length)
        vrna_message_warning("vrna_pbacktrack5: %s", info_no_comparative);
      else if (fc->exp_params->model_details.circ)
        vrna_message_warning("vrna_pbacktrack5: %s", info_no_circ);
      else
        return pbacktrack5_gen(fc, length, NULL);
    } else if (fc->exp_params->model_details.circ) {
      return wrap_pbacktrack_circ(fc);
    } else {
      return pbacktrack5_gen(fc, length, NULL);
    }
  }

  return NULL;
}


PUBLIC unsigned int
vrna_pbacktrack_nr_resume_cb(vrna_fold_compound_t             *vc,
                             unsigned int                     num_samples,
                             vrna_boltzmann_sampling_callback *bs_cb,
                             void                             *data,
                             vrna_nr_memory_t                 *nr_mem)
{
  unsigned int i = 0;

  if (vc) {
    vrna_mx_pf_t *matrices = vc->exp_matrices;

    if ((!matrices) || (!matrices->q) || (!matrices->qb) || (!matrices->qm) || (!vc->exp_params)) {
      vrna_message_warning("vrna_pbacktrack_nr*(): %s", info_call_pf);
    } else if ((!vc->exp_params->model_details.uniq_ML) || (!matrices->qm1)) {
      vrna_message_warning("vrna_pbacktrack_nr*(): %s", info_set_uniq_ml);
    } else if (vc->exp_params->model_details.circ) {
      vrna_message_warning("vrna_pbacktrack_nr*(): %s", info_no_circ);
    } else if (!nr_mem) {
      vrna_message_warning("vrna_pbacktrack_nr*(): Pointer to nr_mem must not be NULL!");
    } else {
      int is_dup, pf_overflow;

      if (*nr_mem == NULL)
        *nr_mem = nr_init(vc);

      pf_overflow = 0;

      for (i = 0; i < num_samples; i++) {
        is_dup = 1;
        char *ss = pbacktrack5_gen(vc, vc->length, *nr_mem);

#ifdef VRNA_NR_SAMPLING_HASH
        (*nr_mem)->current_node = traceback_to_root((*nr_mem)->current_node,
                                                    (*nr_mem)->q_remain,
                                                    &is_dup,
                                                    &pf_overflow);
#else
        (*nr_mem)->current_node = traceback_to_ll_root((*nr_mem)->current_node,
                                                       (*nr_mem)->q_remain,
                                                       &is_dup,
                                                       &pf_overflow);
#endif

        if (pf_overflow) {
          vrna_message_warning("vrna_pbacktrack_nr*(): %s", info_nr_overflow);
          free(ss);
          break;
        }

        if (is_dup) {
          vrna_message_warning("vrna_pbacktrack_nr*(): %s", info_nr_duplicates);
          free(ss);
          break;
        }

        bs_cb(ss, data);
        free(ss);
        /* finish if no more structures available */
        if (!ss)
          break;
      }

      /* print warning if we've aborted backtracking too early */
      if ((i > 0) && (i < num_samples)) {
        vrna_message_warning("vrna_pbacktrack_nr*(): "
                             "Stopped backtracking after %d samples due to numeric instabilities!\n"
                             "Coverage of partition function so far: %.6f%%",
                             i,
                             100. *
                             return_node_weight((*nr_mem)->root_node) /
                             vc->exp_matrices->q[vc->iindx[1] - vc->length]);
      }
    }
  }

  return i; /* actual number of structures backtraced */
}


PUBLIC void
vrna_pbacktrack_nr_free(struct vrna_nr_memory_s *s)
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


PRIVATE struct sc_wrappers *
sc_init(vrna_fold_compound_t *fc)
{
  struct sc_wrappers *sc_wrap = (struct sc_wrappers *)vrna_alloc(sizeof(struct sc_wrappers));

  init_sc_wrapper_ext(fc, &(sc_wrap->sc_wrapper_ext));
  init_sc_wrapper_ml(fc, &(sc_wrap->sc_wrapper_ml));

  return sc_wrap;
}


PRIVATE void
sc_free(struct sc_wrappers *sc_wrap)
{
  free_sc_wrapper_ext(&(sc_wrap->sc_wrapper_ext));
  free_sc_wrapper_ml(&(sc_wrap->sc_wrapper_ml));

  free(sc_wrap);
}


PRIVATE struct vrna_nr_memory_s *
nr_init(vrna_fold_compound_t *fc)
{
  struct vrna_nr_memory_s *s =
    (struct vrna_nr_memory_s *)vrna_alloc(sizeof(struct vrna_nr_memory_s));

  double                  pf          = fc->exp_matrices->q[fc->iindx[1] - fc->length];
  size_t                  block_size  = 5000 * sizeof(NR_NODE);

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
PRIVATE char *
pbacktrack5_gen(vrna_fold_compound_t    *vc,
                int                     length,
                struct vrna_nr_memory_s *nr_mem)
{
  int                 ret, i;
  char                *pstruc;
  struct sc_wrappers  *sc_wrap;

  pstruc = vrna_alloc((length + 1) * sizeof(char));

  for (i = 0; i < length; i++)
    pstruc[i] = '.';

  sc_wrap = sc_init(vc);

  if (nr_mem)
    nr_mem->q_remain = vc->exp_matrices->q[vc->iindx[1] - vc->length];

#ifdef VRNA_WITH_BOUSTROPHEDON
  ret = backtrack_ext_loop(length, pstruc, vc, length, sc_wrap, nr_mem);
#else
  ret = backtrack_ext_loop(1, pstruc, vc, length, sc_wrap, nr_mem);
#endif

  sc_free(sc_wrap);

  if (ret > 0)
    return pstruc;

  free(pstruc);
  return NULL;
}


/* backtrack one external */
PRIVATE int
backtrack_ext_loop(int                      init_val,
                   char                     *pstruc,
                   vrna_fold_compound_t     *vc,
                   int                      length,
                   struct sc_wrappers       *sc_wrap,
                   struct vrna_nr_memory_s  *nr_mem)
{
  unsigned int              **a2s, s, n_seq;
  FLT_OR_DBL                r, fbd, fbds, qt, q_temp, qkl;
  int                       ret, i, j, ij, n, k, u, type;
  int                       *my_iindx, hc_decompose, *hc_up_ext;
  FLT_OR_DBL                *q, *qb, *q1k, *qln, *scale;
  unsigned char             *hard_constraints;
  short                     *S1, *S2, **S, **S5, **S3;
  vrna_mx_pf_t              *matrices;
  vrna_md_t                 *md;
  vrna_hc_t                 *hc;
  vrna_sc_t                 *sc, **scs;
  vrna_exp_param_t          *pf_params;

  double                    *q_remain;
  NR_NODE                   **current_node;
  struct nr_memory          **memory_dat;
  struct sc_wrapper_exp_ext *sc_wrapper_ext;

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
    sc    = vc->sc;
    scs   = NULL;
    S1    = vc->sequence_encoding;
    S2    = vc->sequence_encoding2;
    S     = NULL;
    S5    = NULL;
    S3    = NULL;
    a2s   = NULL;
  } else {
    n_seq = vc->n_seq;
    sc    = NULL;
    scs   = vc->scs;
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

  if (!(q1k && qln)) {
    matrices->q1k = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * (n + 1));
    matrices->qln = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * (n + 2));
    q1k           = matrices->q1k;
    qln           = matrices->qln;
    for (k = 1; k <= n; k++) {
      q1k[k]  = q[my_iindx[1] - k];
      qln[k]  = q[my_iindx[k] - n];
    }
    q1k[0]      = 1.0;
    qln[n + 1]  = 1.0;
  }

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
          qkl   *= exp_E_ExtLoop(type,
                                 (i > 1) ? S1[i - 1] : -1,
                                 (j < n) ? S1[j + 1] : -1,
                                 pf_params);
        } else {
          for (s = 0; s < n_seq; s++) {
            type  = vrna_get_ptype_md(S[s][i], S[s][j], md);
            qkl   *= exp_E_ExtLoop(type,
                                   (a2s[s][i] > 1) ? S5[s][i] : -1,
                                   (a2s[s][j] < a2s[s][S[0][0]]) ? S3[s][j] : -1,
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
      type          = vrna_get_ptype_md(S2[i], S2[j], md);
      hc_decompose  = hard_constraints[n * i + j];
      if (hc_decompose & VRNA_CONSTRAINT_CONTEXT_EXT_LOOP) {
        qkl = qb[ij] *
              exp_E_ExtLoop(type,
                            (i > 1) ? S1[i - 1] : -1,
                            (j < n) ? S1[j + 1] : -1,
                            pf_params);

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
backtrack_qm(int                      i,
             int                      j,
             char                     *pstruc,
             vrna_fold_compound_t     *vc,
             struct sc_wrappers       *sc_wrap,
             struct vrna_nr_memory_s  *nr_mem)
{
  /* divide multiloop into qm and qm1  */
  FLT_OR_DBL                qmt, fbd, fbds, r, q_temp;
  int                       k, u, cnt, span, turn;
  int                       is_unpaired; /* 1 if [i ... k-1] is unpaired */
  FLT_OR_DBL                *qm, *qm1, *expMLbase;
  int                       *my_iindx, *jindx, *hc_up_ml, ret;
  vrna_hc_t                 *hc;
  struct sc_wrapper_exp_ml  *sc_wrapper_ml;

  double                    *q_remain;
  NR_NODE                   **current_node;
  struct nr_memory          **memory_dat;

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

  vrna_mx_pf_t *matrices = vc->exp_matrices;

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
backtrack_qm1(int                     i,
              int                     j,
              char                    *pstruc,
              vrna_fold_compound_t    *vc,
              struct sc_wrappers      *sc_wrap,
              struct vrna_nr_memory_s *nr_mem)
{
  /* i is paired to l, i<l<j; backtrack in qm1 to find l */
  unsigned int              n, s, n_seq;
  int                       ii, l, il, type, turn;
  FLT_OR_DBL                qt, fbd, fbds, r, q_temp;
  FLT_OR_DBL                *qm1, *qb, *expMLbase;
  vrna_mx_pf_t              *matrices;
  int                       u, *my_iindx, *jindx, *hc_up_ml;
  char                      *ptype;
  unsigned char             *hard_constraints;
  short                     *S1, **S, **S5, **S3;
  vrna_hc_t                 *hc;
  vrna_exp_param_t          *pf_params;
  vrna_md_t                 *md;

  double                    *q_remain;
  NR_NODE                   **current_node;
  struct nr_memory          **memory_dat;
  struct sc_wrapper_exp_ml  *sc_wrapper_ml;

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
  FLT_OR_DBL                qom2t, r;
  int                       u, turn;
  FLT_OR_DBL                *qm1, *qm2;
  int                       *jindx;
  struct sc_wrapper_exp_ml  *sc_wrapper_ml;

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
backtrack(int                     i,
          int                     j,
          char                    *pstruc,
          vrna_fold_compound_t    *vc,
          struct sc_wrappers      *sc_wrap,
          struct vrna_nr_memory_s *nr_mem)
{
  char                      *ptype;
  unsigned char             *hard_constraints, hc_decompose;
  vrna_exp_param_t          *pf_params;
  vrna_md_t                 *md;
  FLT_OR_DBL                *qb, *qm, *qm1, *scale;
  FLT_OR_DBL                r, fbd, fbds, qbt1, qbr, qt, q_temp; /* qbr stores qb used for generating r */
  vrna_mx_pf_t              *matrices;
  unsigned int              **a2s, s, n_seq, n, type, type_2, *types;
  int                       *my_iindx, *jindx, *hc_up_int, ret;
  vrna_sc_t                 *sc, **scs;
  vrna_hc_t                 *hc;
  short                     *S1, **S, **S5, **S3;
  int                       *pscore; /* precomputed array of pair types */

  FLT_OR_DBL                kTn;
  double                    *q_remain;
  NR_NODE                   **current_node;
  struct nr_memory          **memory_dat;
  struct sc_wrapper_exp_ml  *sc_wrapper_ml;

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

  if (vc->type == VRNA_FC_TYPE_SINGLE) {
    n_seq   = 1;
    ptype   = vc->ptype;
    types   = NULL;
    pscore  = NULL;
    S1      = vc->sequence_encoding;
    S       = NULL;
    S5      = NULL;
    S3      = NULL;
    a2s     = NULL;
    sc      = vc->sc;
    scs     = NULL;
  } else {
    n_seq   = vc->n_seq;
    ptype   = NULL;
    types   = (unsigned int *)vrna_alloc(sizeof(unsigned int) * n_seq);
    pscore  = vc->pscore;
    S1      = NULL;
    S       = vc->S;
    S5      = vc->S5;
    S3      = vc->S3;
    a2s     = vc->a2s;
    sc      = NULL;
    scs     = vc->scs;
  }

  hc                = vc->hc;
  hc_up_int         = hc->up_int;
  hard_constraints  = hc->mx;
  sc_wrapper_ml     = &(sc_wrap->sc_wrapper_ml);

  matrices  = vc->exp_matrices;
  qb        = matrices->qb;
  qm        = matrices->qm;
  qm1       = matrices->qm1;
  scale     = matrices->scale;

  int turn    = pf_params->model_details.min_loop_size;
  int *rtype  = &(pf_params->model_details.rtype[0]);

#ifndef VRNA_NR_SAMPLING_HASH
  if (current_node) {
    memorized_node_prev = NULL;
    memorized_node_cur  = (*current_node)->head;
  }

#endif

  hc_decompose = hard_constraints[n * j + i];

  do {
    int           k, l, kl, u1, u2, max_k, min_l;
    unsigned char type;
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

              if (sc) {
                if (sc->exp_energy_up)
                  q_temp *= sc->exp_energy_up[i + 1][u1]
                            * sc->exp_energy_up[l + 1][u2];

                if (sc->exp_energy_bp)
                  q_temp *= sc->exp_energy_bp[jindx[j] + i];

                if (sc->exp_energy_stack) {
                  if ((i + 1 == k) && (j - 1 == l)) {
                    q_temp *= sc->exp_energy_stack[i]
                              * sc->exp_energy_stack[k]
                              * sc->exp_energy_stack[l]
                              * sc->exp_energy_stack[j];
                  }
                }

                if (sc->exp_f)
                  q_temp *= sc->exp_f(i, j, k, l, VRNA_DECOMP_PAIR_IL, sc->data);
              }
            } else {
              for (s = 0; s < n_seq; s++) {
                int u1_local  = a2s[s][k - 1] - a2s[s][i] /*??*/;
                int u2_local  = a2s[s][j - 1] - a2s[s][l];
                type_2  = vrna_get_ptype_md(S[s][l], S[s][k], md);
                q_temp  *= exp_E_IntLoop(u1_local,
                                         u2_local,
                                         types[s],
                                         type_2,
                                         S3[s][i],
                                         S5[s][j],
                                         S5[s][k],
                                         S3[s][l],
                                         pf_params);
              }

              if (scs) {
                for (s = 0; s < n_seq; s++) {
                  if (scs[s]) {
                    if (scs[s]->exp_energy_stack) {
                      int u1_local  = a2s[s][k - 1] - a2s[s][i];
                      int u2_local  = a2s[s][j - 1] - a2s[s][l];
                      if (u1_local + u2_local == 0) {
                        if (S[s][i] && S[s][j] && S[s][k] && S[s][l]) {
                          /* don't allow gaps in stack */
                          q_temp *= scs[s]->exp_energy_stack[a2s[s][i]] *
                                    scs[s]->exp_energy_stack[a2s[s][k]] *
                                    scs[s]->exp_energy_stack[a2s[s][l]] *
                                    scs[s]->exp_energy_stack[a2s[s][j]];
                        }
                      }
                    }
                  }
                }
              }
            }

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
    int         k, ii, jj, tt;
    FLT_OR_DBL  closingPair;

    closingPair = scale[2];

    if (vc->type == VRNA_FC_TYPE_SINGLE) {
      tt          = rtype[vrna_get_ptype(jindx[j] + i, ptype)];
      closingPair *= exp_E_MLstem(tt, S1[j - 1], S1[i + 1], pf_params) *
                     pf_params->expMLclosing;
    } else {
      closingPair *= pow(pf_params->expMLclosing, (double)n_seq);
      for (s = 0; s < n_seq; s++) {
        tt          = vrna_get_ptype_md(S[s][j], S[s][i], md);
        closingPair *= exp_E_MLstem(tt, S5[s][j], S3[s][i], pf_params);
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


PRIVATE char *
wrap_pbacktrack_circ(vrna_fold_compound_t *vc)
{
  FLT_OR_DBL                r, qt, q_temp;
  int                       i, j, k, l, n;
  vrna_exp_param_t          *pf_params;
  FLT_OR_DBL                qo, qmo;
  FLT_OR_DBL                *scale, *qb, *qm, *qm2;
  char                      *sequence, *ptype, *pstruc;
  int                       *my_iindx, *jindx;
  short                     *S1;
  vrna_mx_pf_t              *matrices;
  vrna_sc_t                 *sc;
  struct sc_wrappers        *sc_wrap;
  struct sc_wrapper_exp_ext *sc_wrapper_ext;
  struct sc_wrapper_exp_ml  *sc_wrapper_ml;

  pf_params = vc->exp_params;
  matrices  = vc->exp_matrices;
  ptype     = vc->ptype;
  my_iindx  = vc->iindx;
  jindx     = vc->jindx;
  S1        = vc->sequence_encoding;
  sc        = vc->sc;

  qo    = matrices->qo;
  qmo   = matrices->qmo;
  qb    = matrices->qb;
  qm    = matrices->qm;
  qm2   = matrices->qm2;
  scale = matrices->scale;

  FLT_OR_DBL  expMLclosing  = pf_params->expMLclosing;
  int         *rtype        = &(pf_params->model_details.rtype[0]);
  int         turn          = pf_params->model_details.min_loop_size;

  sequence  = vc->sequence;
  n         = vc->length;


  sc_wrap         = sc_init(vc);
  sc_wrapper_ext  = &(sc_wrap->sc_wrapper_ext);
  sc_wrapper_ml   = &(sc_wrap->sc_wrapper_ml);

  /*
   * if (init_length<1)
   *  vrna_message_error("can't backtrack without pf arrays.\n"
   *    "Call pf_circ_fold() before pbacktrack_circ()");
   */

  pstruc = vrna_alloc((n + 1) * sizeof(char));

  /* initialize pstruct with single bases  */
  for (i = 0; i < n; i++)
    pstruc[i] = '.';

  qt = 1.0 * scale[n];

  /* add soft constraints for open chain configuration */
  if (sc_wrapper_ext->red_up)
    qt *= sc_wrapper_ext->red_up(1, n, sc_wrapper_ext);

  r = vrna_urn() * qo;

  /* open chain? */
  if (qt > r) {
    sc_free(sc_wrap);
    return pstruc;
  }

  for (i = 1; (i < n); i++) {
    for (j = i + turn + 1; (j <= n); j++) {
      int type, u;
      /* 1. first check, wether we can do a hairpin loop  */
      u = n - j + i - 1;
      if (u < turn)
        continue;

      type = ptype[jindx[j] + i];
      if (!type)
        continue;

      type = rtype[type];

      qt += qb[my_iindx[i] - j] *
            vrna_exp_E_hp_loop(vc, j, i);

      /* found a hairpin? so backtrack in the enclosed part and we're done  */
      if (qt > r) {
        backtrack(i, j, pstruc, vc, sc_wrap, NULL);
        sc_free(sc_wrap);
        return pstruc;
      }

      /* 2. search for (k,l) with which we can close an interior loop  */
      for (k = j + 1; (k < n); k++) {
        int ln1, lstart;
        ln1 = k - j - 1;
        if (ln1 + i - 1 > MAXLOOP)
          break;

        lstart = ln1 + i - 1 + n - MAXLOOP;
        if (lstart < k + turn + 1)
          lstart = k + turn + 1;

        for (l = lstart; (l <= n); l++) {
          int ln2, ln3, type2;
          ln2 = (i - 1);
          ln3 = (n - l);
          if ((ln1 + ln2 + ln3) > MAXLOOP)
            continue;

          type2 = ptype[jindx[l] + k];
          if (!type)
            continue;

          type2   = rtype[type2];
          q_temp  = qb[my_iindx[i] - j] *
                    qb[my_iindx[k] - l] *
                    exp_E_IntLoop(ln2 + ln3,
                                  ln1,
                                  type2,
                                  type,
                                  S1[l + 1],
                                  S1[k - 1],
                                  S1[i - 1],
                                  S1[j + 1],
                                  pf_params) *
                    scale[ln1 + ln2 + ln3];

          if (sc) {
            if (sc->exp_energy_up) {
              q_temp *= sc->exp_energy_up[1][ln2] *
                        sc->exp_energy_up[l + 1][ln3] *
                        sc->exp_energy_up[j + 1][ln1];
            }

            if ((sc->exp_energy_stack) && (ln1 + ln2 + ln3 == 0)) {
              q_temp *= sc->exp_energy_stack[i] *
                        sc->exp_energy_stack[j] *
                        sc->exp_energy_stack[k] *
                        sc->exp_energy_stack[l];
            }

            if (sc->exp_f)
              q_temp *= sc->exp_f(i, j, k, l, VRNA_DECOMP_PAIR_IL, sc->data);
          }

          qt += q_temp;
          /* found an exterior interior loop? also this time, we can go straight  */
          /* forward and backtracking the both enclosed parts and we're done      */
          if (qt > r) {
            backtrack(i, j, pstruc, vc, sc_wrap, NULL);
            backtrack(k, l, pstruc, vc, sc_wrap, NULL);
            sc_free(sc_wrap);
            return pstruc;
          }
        }
      } /* end of kl double loop */
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
          sc_free(sc_wrap);
          return pstruc;
        }
      }
    } else {
      for (k = turn + 2; k < n - 2 * turn - 3; k++) {
        qt += qm[my_iindx[1] - k] * qm2[k + 1] * expMLclosing;
        /* backtrack in qm and qm2 if we've found a valid barrier k  */
        if (qt > r) {
          backtrack_qm(1, k, pstruc, vc, sc_wrap, NULL);
          backtrack_qm2(k + 1, n, pstruc, vc, sc_wrap);
          sc_free(sc_wrap);
          return pstruc;
        }
      }
    }
  }
  /* if we reach the actual end of this function, an error has occured  */
  /* cause we HAVE TO find an exterior loop or an open chain!!!         */
  vrna_message_error("backtracking failed in exterior loop");

  sc_free(sc_wrap);

  return pstruc;
}
