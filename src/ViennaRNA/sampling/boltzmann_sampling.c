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
#include "ViennaRNA/utils/log.h"
#include "ViennaRNA/params/default.h"
#include "ViennaRNA/fold_vars.h"
#include "ViennaRNA/eval/exterior.h"
#include "ViennaRNA/eval/hairpin.h"
#include "ViennaRNA/eval/internal.h"
#include "ViennaRNA/eval/multibranch.h"
#include "ViennaRNA/eval/gquad.h"
#include "ViennaRNA/constraints/hard.h"
#include "ViennaRNA/constraints/soft.h"
#include "ViennaRNA/sequences/alphabet.h"
#include "ViennaRNA/combinatorics/basic.h"
#include "ViennaRNA/sampling/basic.h"

#include "ViennaRNA/constraints/exterior_sc_pf.inc"
#include "ViennaRNA/constraints/internal_sc_pf.inc"
#include "ViennaRNA/constraints/multibranch_sc_pf.inc"

#include "ViennaRNA/sampling/data_structures_nonred.inc"

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


struct aux_mem {
  FLT_OR_DBL *qik;
};

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
  unsigned int      start;
  unsigned int      end;
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
nr_init(vrna_fold_compound_t  *fc,
        unsigned int          start,
        unsigned int          end);


PRIVATE struct sc_wrappers *
sc_init(vrna_fold_compound_t *fc);


PRIVATE void
sc_free(struct sc_wrappers *sc_wrap);


PRIVATE unsigned int
wrap_pbacktrack(vrna_fold_compound_t            *vc,
                unsigned int                    start,
                unsigned int                    end,
                unsigned int                    num_samples,
                vrna_bs_result_f                bs_cb,
                void                            *data,
                struct vrna_pbacktrack_memory_s *nr_mem);


PRIVATE int
backtrack(unsigned int                    i,
          unsigned int                    j,
          char                            *pstruc,
          vrna_fold_compound_t            *vc,
          struct sc_wrappers              *sc_wrap,
          struct vrna_pbacktrack_memory_s *nr_mem);


PRIVATE int
backtrack_ext_loop(unsigned int                     start,
                   unsigned int                     end,
                   char                             *pstruc,
                   vrna_fold_compound_t             *vc,
                   struct aux_mem                   *helper_arrays,
                   struct sc_wrappers               *sc_wrap,
                   struct vrna_pbacktrack_memory_s  *nr_mem);


PRIVATE int
backtrack_qm(unsigned int                     i,
             unsigned int                     j,
             char                             *pstruc,
             vrna_fold_compound_t             *vc,
             struct sc_wrappers               *sc_wrap,
             struct vrna_pbacktrack_memory_s  *nr_mem);


PRIVATE int
backtrack_qm1(vrna_fold_compound_t            *fc,
              unsigned int                    j,
              char                            *pstruc,
              struct sc_wrappers              *sc_wrap,
              struct vrna_pbacktrack_memory_s *nr_mem);


PRIVATE int
backtrack_qm2(vrna_fold_compound_t            *fc,
              unsigned int                    i,
              unsigned int                    j,
              char                            *pstruc,
              struct sc_wrappers              *sc_wrap,
              struct vrna_pbacktrack_memory_s *nr_mem);


PRIVATE unsigned int
pbacktrack_circ(vrna_fold_compound_t  *fc,
                unsigned int          num_samples,
                vrna_bs_result_f      bs_cb,
                void                  *data);


/*
 #################################
 # BEGIN OF FUNCTION DEFINITIONS #
 #################################
 */
PUBLIC unsigned int
vrna_pbacktrack_sub_resume_cb(vrna_fold_compound_t  *fc,
                              unsigned int          num_samples,
                              unsigned int          start,
                              unsigned int          end,
                              vrna_bs_result_f      bs_cb,
                              void                  *data,
                              vrna_pbacktrack_mem_t *nr_mem,
                              unsigned int          options)
{
  unsigned int i = 0;

  if (fc) {
    vrna_mx_pf_t *matrices = fc->exp_matrices;

    if (start == 0) {
      vrna_log_warning("vrna_pbacktrack*(): interval start coordinate must be at least 1");
    } else if (end > fc->length) {
      vrna_log_warning("vrna_pbacktrack*(): interval end coordinate exceeds sequence length");
    } else if (end < start) {
      vrna_log_warning("vrna_pbacktrack*(): interval end < start");
    } else if ((!matrices) || (!matrices->q) || (!matrices->qb) || (!matrices->qm) ||
               (!fc->exp_params)) {
      vrna_log_warning("vrna_pbacktrack*(): %s", info_call_pf);
    } else if ((!fc->exp_params->model_details.uniq_ML) ||
               ((!matrices->qm1) && (!matrices->qm2_real))) {
      vrna_log_warning("vrna_pbacktrack*(): %s", info_set_uniq_ml);
    } else if ((fc->exp_params->model_details.circ) && (end < fc->length)) {
      vrna_log_warning("vrna_pbacktrack5*(): %s", info_no_circ);
    } else if (options & VRNA_PBACKTRACK_NON_REDUNDANT) {
      if (fc->exp_params->model_details.circ) {
        vrna_log_warning("vrna_pbacktrack5*(): %s", info_no_circ);
      } else if (!nr_mem) {
        vrna_log_warning("vrna_pbacktrack5*(): Pointer to nr_mem must not be NULL!");
      } else {
        if ((*nr_mem == NULL) ||
            ((*nr_mem)->start != start) ||
            ((*nr_mem)->end != end)) {
          if (*nr_mem)
            vrna_pbacktrack_mem_free(*nr_mem);

          *nr_mem = nr_init(fc, start, end);
        }

        i = wrap_pbacktrack(fc, start, end, num_samples, bs_cb, data, *nr_mem);

        /* print warning if we've aborted backtracking too early */
        if ((i > 0) && (i < num_samples)) {
          vrna_log_warning("vrna_pbacktrack5*(): "
                           "Stopped non-redundant backtracking after %d samples"
                           " due to numeric instabilities!\n"
                           "Coverage of partition function so far: %.6f%%",
                           i,
                           100. *
                           return_node_weight((*nr_mem)->root_node) /
                           fc->exp_matrices->q[fc->iindx[start] - end]);
        }
      }
    } else if (fc->exp_params->model_details.circ) {
      i = pbacktrack_circ(fc, num_samples, bs_cb, data);
    } else {
      i = wrap_pbacktrack(fc, start, end, num_samples, bs_cb, data, NULL);
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
nr_init(vrna_fold_compound_t  *fc,
        unsigned int          start,
        unsigned int          end)
{
  size_t                          block_size;
  double                          pf;
  struct vrna_pbacktrack_memory_s *s;

  s = (struct vrna_pbacktrack_memory_s *)vrna_alloc(
    sizeof(struct vrna_pbacktrack_memory_s));

  s->start      = start;
  s->end        = end;
  s->memory_dat = NULL;
  s->q_remain   = 0;

  pf          = fc->exp_matrices->q[fc->iindx[start] - end];
  block_size  = 5000 * sizeof(NR_NODE);

#ifdef VRNA_NR_SAMPLING_HASH
  s->root_node = create_root(end, pf);
#else
  s->memory_dat = create_nr_memory(sizeof(NR_NODE), block_size, NULL);  /* memory pre-allocation */
  s->root_node  = create_ll_root(&(s->memory_dat), pf);
#endif

  s->current_node = s->root_node;

  return s;
}


/* general expr of vrna5_pbacktrack with possibility of non-redundant sampling */
PRIVATE unsigned int
wrap_pbacktrack(vrna_fold_compound_t            *vc,
                unsigned int                    start,
                unsigned int                    end,
                unsigned int                    num_samples,
                vrna_bs_result_f                bs_cb,
                void                            *data,
                struct vrna_pbacktrack_memory_s *nr_mem)
{
  char                *pstruc;
  unsigned int        i;
  int                 ret, pf_overflow, is_dup, *my_iindx;
  FLT_OR_DBL          *q;
  vrna_mx_pf_t        *matrices;
  struct aux_mem      helper_arrays;
  struct sc_wrappers  *sc_wrap;

  i           = 0;
  pf_overflow = 0;
  sc_wrap     = sc_init(vc);

  my_iindx  = vc->iindx;
  matrices  = vc->exp_matrices;
  q         = matrices->q;

  helper_arrays.qik = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * (end - start + 2));
  helper_arrays.qik -= start - 1;

  for (i = start; i <= end; i++)
    helper_arrays.qik[i] = q[my_iindx[start] - i];

  helper_arrays.qik[start - 1] = 1.0;

  for (i = 0; i < num_samples; i++) {
    is_dup  = 1;
    pstruc  = vrna_alloc(((end - start + 1) + 1) * sizeof(char));
    memset(pstruc, '.', sizeof(char) * (end - start + 1));
    pstruc -= start - 1;


    if (nr_mem)
      nr_mem->q_remain = vc->exp_matrices->q[vc->iindx[start] - end]; /* really */

    ret = backtrack_ext_loop(start, end, pstruc, vc, &helper_arrays, sc_wrap, nr_mem);

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
        vrna_log_warning("vrna_pbacktrack_nr*(): %s", info_nr_overflow);
        free(pstruc + (start - 1));
        break;
      }

      if (is_dup) {
        vrna_log_warning("vrna_pbacktrack_nr*(): %s", info_nr_duplicates);
        free(pstruc + (start - 1));
        break;
      }
    }

    if ((ret > 0) && (bs_cb))
      bs_cb(pstruc + (start - 1), data);

    free(pstruc + (start - 1));

    if (ret == 0)
      break;
  }

  helper_arrays.qik += start - 1;
  free(helper_arrays.qik);

  sc_free(sc_wrap);

  return i;
}


/* backtrack one external */
PRIVATE int
backtrack_ext_loop(unsigned int                     start,
                   unsigned int                     end,
                   char                             *pstruc,
                   vrna_fold_compound_t             *vc,
                   struct aux_mem                   *helper_arrays,
                   struct sc_wrappers               *sc_wrap,
                   struct vrna_pbacktrack_memory_s  *nr_mem)
{
  unsigned char         *hard_constraints, hc_decompose;
  short                 *S1, *S2, **S, **S5, **S3;
  unsigned int          **a2s, s, n_seq, *hc_up_ext, n, k, type, i, j;
  int                   ret, ij, *my_iindx;
  FLT_OR_DBL            r, fbd, fbds, qt, q_temp, qkl, *qb, *q1k, *scale;
  double                *q_remain;
  vrna_mx_pf_t          *matrices;
  vrna_md_t             *md;
  vrna_hc_t             *hc;
  vrna_exp_param_t      *pf_params;

  struct nr_memory      **memory_dat;
  struct sc_ext_exp_dat *sc_wrapper_ext;

  NR_NODE               **current_node;

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

  qb    = matrices->qb;
  q1k   = helper_arrays->qik;
  scale = matrices->scale;

#ifndef VRNA_NR_SAMPLING_HASH
  if (current_node) {
    memorized_node_prev = NULL;
    memorized_node_cur  = (*current_node)->head;
  }

#endif

  q_temp = 0.;

  j = end;
  if (j > start) {
    /* find j position of first pair */
    for (; j > start; j--) {
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
          q_temp *= sc_wrapper_ext->red_ext(start, j, start, j - 1, sc_wrapper_ext);

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
      } else {
        /* j must be paired, so continue by finding its pairing partner */
        break;
      }
    }
    if (j <= start + md->min_loop_size)
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
    i = 2;

    unsigned int *is = vrna_boustrophedon(start, j - 1);

    for (qt = 0, k = 1; k + start <= j; k++) {
      /* apply alternating boustrophedon scheme to variable i */
      i             = is[k];
      ij            = my_iindx[i] - j;
      hc_decompose  = hard_constraints[n * j + i];
      if (hc_decompose & VRNA_CONSTRAINT_CONTEXT_EXT_LOOP) {
        qkl = qb[ij] *
              q1k[i - 1];

        if (vc->type == VRNA_FC_TYPE_SINGLE) {
          type  = vrna_get_ptype_md(S2[i], S2[j], md);
          qkl   *= vrna_exp_E_exterior_stem(type,
                                            (i > 1) ? S1[i - 1] : -1,
                                            (j < n) ? S1[j + 1] : -1,
                                            pf_params);
        } else {
          for (s = 0; s < n_seq; s++) {
            type  = vrna_get_ptype_md(S[s][i], S[s][j], md);
            qkl   *= vrna_exp_E_exterior_stem(type,
                                              (a2s[s][i] > 1) ? S5[s][i] : -1,
                                              (a2s[s][j] < a2s[s][n]) ? S3[s][j] : -1,
                                              pf_params);
          }
        }

        if ((sc_wrapper_ext->red_stem) && (i == 1))
          qkl *= sc_wrapper_ext->red_stem(i, j, i, j, sc_wrapper_ext);
        else if ((sc_wrapper_ext->split) && (i > 1))
          qkl *= sc_wrapper_ext->split(1, j, i, sc_wrapper_ext) *
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

    free(is);

    if (k + start > j) {
      if (current_node) {
        /* exhausted ensemble */
        return 0;
      } else {
        vrna_log_warning("backtracking failed in ext loop");
        /* error */
        return -1;
      }
    }

    backtrack(i, j, pstruc, vc, sc_wrap, nr_mem);
    j   = i - 1;
    ret = backtrack_ext_loop(start, j, pstruc, vc, helper_arrays, sc_wrap, nr_mem);
  }

  return ret;
}


/* non redundant version of function bactrack_qm */
PRIVATE int
backtrack_qm(unsigned int                     i,
             unsigned int                     j,
             char                             *pstruc,
             vrna_fold_compound_t             *fc,
             struct sc_wrappers               *sc_wrap,
             struct vrna_pbacktrack_memory_s  *nr_mem)
{
  unsigned char         *hard_constraints;
  short                 *S, *S1, **SS, **S5, **S3;
  unsigned int          k, n, n_seq, s, turn, u, unpaired, type, *hc_up_ml;
  int                   *my_iindx, ret = 0;
  FLT_OR_DBL            qt, q_temp, qbt1, qm_rem, r, *qb, *qm, *expMLbase, fbd, fbds;
  double                *q_remain;
  vrna_exp_param_t      *pf_params;
  vrna_md_t             *md;
  vrna_hc_t             *hc;
  struct sc_mb_exp_dat  *sc_wrapper;
  struct nr_memory      **memory_dat;

  NR_NODE               **current_node;

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

  fbd   = 0.;                       /* stores weight of forbidden terms for given q[ij]*/
  fbds  = 0.;                       /* stores weight of forbidden term for given motif */

  n     = fc->length;
  n_seq = (fc->type == VRNA_FC_TYPE_COMPARATIVE) ? fc->n_seq : 1;
  S     = (fc->type == VRNA_FC_TYPE_COMPARATIVE) ? NULL : fc->sequence_encoding2;
  S1    = (fc->type == VRNA_FC_TYPE_COMPARATIVE) ? NULL : fc->sequence_encoding;
  SS    = (fc->type == VRNA_FC_TYPE_COMPARATIVE) ? fc->S : NULL;
  S5    = (fc->type == VRNA_FC_TYPE_COMPARATIVE) ? fc->S5 : NULL;
  S3    = (fc->type == VRNA_FC_TYPE_COMPARATIVE) ? fc->S3 : NULL;

  pf_params         = fc->exp_params;
  md                = &(pf_params->model_details);
  turn              = md->min_loop_size;
  my_iindx          = fc->iindx;
  qb                = fc->exp_matrices->qb;
  qm                = fc->exp_matrices->qm;
  expMLbase         = fc->exp_matrices->expMLbase;
  hc                = fc->hc;
  hc_up_ml          = hc->up_ml;
  hard_constraints  = hc->mx;
  sc_wrapper        = &(sc_wrap->sc_wrapper_ml);
  unpaired          = 0;

  qt      = 0.;
  q_temp  = 0.;

  if (i + turn >= j) {
    vrna_log_error("backtracking impossible for qm[%u, %u]\n%s", i, j, pstruc);
    return 0; /* error */
  }

#ifndef VRNA_NR_SAMPLING_HASH
  if (current_node) {
    memorized_node_prev = NULL;
    memorized_node_cur  = (*current_node)->head;
  }

#endif

  q_temp = 0.;

  /* nibble-off unpaired bases from 3' side */
  do {
    if (hc_up_ml[j] > 0) {
      /* get already consumed Boltzmann weight */
      if (current_node) {
        fbd = NR_TOTAL_WEIGHT(*current_node) *
              qm[my_iindx[i] - j] /
              (*q_remain);
#ifdef  USE_FLOAT_PF
        if (fabsf(NR_TOTAL_WEIGHT(*current_node) - (*q_remain)) / (*q_remain) <= FLT_EPSILON)
#else
        if (fabs(NR_TOTAL_WEIGHT(*current_node) - (*q_remain)) / (*q_remain) <= DBL_EPSILON)
#endif
          /* exhausted ensemble */
          return 0;
      }

      r = vrna_urn() *
          (qm[my_iindx[i] - j] - fbd);

      q_temp = qm[my_iindx[i] - j + 1] *
               expMLbase[1];

      if (sc_wrapper->red_ml)
        q_temp *= sc_wrapper->red_ml(i, j, i, j - 1, sc_wrapper);

      if (current_node) {
        fbds = NR_GET_WEIGHT(*current_node, memorized_node_cur, NRT_QM_UNPAIR, j - 1, j) *
               qm[my_iindx[i] - j] /
               (*q_remain);

        qt = q_temp - fbds;
      } else {
        qt = q_temp;
      }

      if (r > qt) {
        break; /* j is paired */
      } else if (current_node) {
        *q_remain *= q_temp / qm[my_iindx[i] - j];
#ifdef VRNA_NR_SAMPLING_HASH
        *current_node = add_if_nexists(NRT_QM_UNPAIR,
                                       j - 1,
                                       j,
                                       *current_node,
                                       *q_remain);
#else
        *current_node = add_if_nexists_ll(memory_dat,
                                          NRT_QM_UNPAIR,
                                          j - 1,
                                          j,
                                          memorized_node_prev,
                                          memorized_node_cur,
                                          *current_node,
                                          *q_remain);
        reset_cursor(&memorized_node_prev, &memorized_node_cur, *current_node); /* resets cursor */
#endif
      }

      j--;
    } else {
      /* j must be paired, so continue by finding its pairing partner */
      break;
    }
  } while (i + turn < j);

  if (i + turn == j) {
    vrna_log_error("backtracking failed for qm");
    return 0; /* error */
  }

#ifndef  VRNA_NR_SAMPLING_HASH
  if (current_node)
    advance_cursor(&memorized_node_prev, &memorized_node_cur, NRT_QM_UNPAIR, j - 1, j);

#endif

  /* j is paired, let's find the pairing partner */
  qt      = 0.;
  qm_rem  = qm[my_iindx[i] - j] - q_temp;

  if (current_node) {
    fbd = NR_TOTAL_WEIGHT(*current_node) *
          qm[my_iindx[i] - j] /
          (*q_remain);
  }

  r = vrna_urn() * (qm_rem - fbd);

  /* find split into qm + qb or unpaired + qb */
  for (k = i; k + turn < j; k++) {
    if (hard_constraints[n * j + k] & VRNA_CONSTRAINT_CONTEXT_MB_LOOP_ENC) {
      qbt1 = qb[my_iindx[k] - j];

      if (fc->type == VRNA_FC_TYPE_SINGLE) {
        type  = vrna_get_ptype_md(S[k], S[j], md);
        qbt1  *= vrna_exp_E_multibranch_stem(type, S1[k - 1], S1[j + 1], pf_params);
      } else {
        for (s = 0; s < n_seq; s++) {
          type  = vrna_get_ptype_md(SS[s][k], SS[s][j], md);
          qbt1  *= vrna_exp_E_multibranch_stem(type, S5[s][k], S3[s][j], pf_params);
        }
      }

      /* 1. unpaired + qb */
      u = k - i;
      if (hc_up_ml[i] >= u) {
        q_temp = qbt1 *
                 expMLbase[u];

        if (sc_wrapper->red_stem)
          q_temp *= sc_wrapper->red_stem(i, j, k, j, sc_wrapper);

        if (current_node) {
          fbds = NR_GET_WEIGHT(*current_node, memorized_node_cur, NRT_QM_NOBRANCH, k, 0) *
                 qm[my_iindx[i] - j] /
                 (*q_remain);
          qt += q_temp - fbds;
        } else {
          qt += q_temp;
        }

        if (qt >= r) {
          unpaired = 1;
          break;
        }
      }

#ifndef VRNA_NR_SAMPLING_HASH
      if (current_node)
        advance_cursor(&memorized_node_prev, &memorized_node_cur, NRT_QM_NOBRANCH, k, 0);

#endif

      /* 2. at least one more branch + qb */
      if (k > i) {
        q_temp = qbt1 *
                 qm[my_iindx[i] - k + 1];

        if (sc_wrapper->decomp_ml)
          q_temp *= sc_wrapper->decomp_ml(i, j, k - 1, k, sc_wrapper);

        if (sc_wrapper->red_stem)
          q_temp *= sc_wrapper->red_stem(k, j, k, j, sc_wrapper);

        if (current_node) {
          fbds = NR_GET_WEIGHT(*current_node, memorized_node_cur, NRT_QM_BRANCH, k, 0) *
                 qm[my_iindx[i] - j] /
                 (*q_remain);
          qt += q_temp - fbds;
        } else {
          qt += q_temp;
        }

        if (qt >= r)
          break;
      }

#ifndef VRNA_NR_SAMPLING_HASH
      if (current_node)
        advance_cursor(&memorized_node_prev, &memorized_node_cur, NRT_QM_BRANCH, k, 0);

#endif
    }
  }

  if (current_node) {
    *q_remain *= q_temp / qm[my_iindx[i] - j];
#ifdef VRNA_NR_SAMPLING_HASH
    if (unpaired) {
      *current_node = add_if_nexists(NRT_QM_NOBRANCH,
                                     k,
                                     0,
                                     *current_node,
                                     *q_remain);
    } else {
      *current_node = add_if_nexists(NRT_QM_BRANCH,
                                     k,
                                     0,
                                     *current_node,
                                     *q_remain);
    }

#else
    if (unpaired) {
      *current_node = add_if_nexists_ll(memory_dat,
                                        NRT_QM_NOBRANCH,
                                        k,
                                        0,
                                        memorized_node_prev,
                                        memorized_node_cur,
                                        *current_node,
                                        *q_remain);
    } else {
      *current_node = add_if_nexists_ll(memory_dat,
                                        NRT_QM_BRANCH,
                                        k,
                                        0,
                                        memorized_node_prev,
                                        memorized_node_cur,
                                        *current_node,
                                        *q_remain);
    }

#endif
  }

  if (k + turn >= j) {
    vrna_log_error("backtracking failed for qm 2");
    return 0;
  }

  /* backtrack further multiloop parts */
  if (!unpaired)
    ret = backtrack_qm(i, k - 1, pstruc, fc, sc_wrap, nr_mem);

  /* backtrack stem */
  ret &= backtrack(k, j, pstruc, fc, sc_wrap, nr_mem);

  return ret;
}


PRIVATE int
backtrack_qm1(vrna_fold_compound_t            *fc,
              unsigned int                    j,
              char                            *pstruc,
              struct sc_wrappers              *sc_wrap,
              struct vrna_pbacktrack_memory_s *nr_mem)
{
  /* i is paired to l, i<l<j; backtrack in qm1 to find l */
  unsigned char         *hard_constraints;
  short                 *S1, *S, **SS, **S5, **S3;
  unsigned int          n, s, n_seq, *hc_up_ml, type, turn, u, l;
  int                   *my_iindx;
  FLT_OR_DBL            qt, qm1j, fbd, fbds, r, q_temp, *qm1, *qb, *expMLbase;
  double                *q_remain;
  vrna_exp_param_t      *pf_params;
  vrna_md_t             *md;
  vrna_hc_t             *hc;
  vrna_mx_pf_t          *matrices;

  struct nr_memory      **memory_dat;
  struct sc_mb_exp_dat  *sc_wrapper_ml;

  NR_NODE               **current_node;

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

  n                 = fc->length;
  fbd               = 0.;
  fbds              = 0.;
  pf_params         = fc->exp_params;
  md                = &(pf_params->model_details);
  my_iindx          = fc->iindx;
  hc                = fc->hc;
  hc_up_ml          = hc->up_ml;
  hard_constraints  = hc->mx;
  sc_wrapper_ml     = &(sc_wrap->sc_wrapper_ml);

  matrices  = fc->exp_matrices;
  qb        = matrices->qb;
  qm1       = matrices->qm1_new;
  expMLbase = matrices->expMLbase;

  if (fc->type == VRNA_FC_TYPE_SINGLE) {
    n_seq = 1;
    S1    = fc->sequence_encoding;
    S     = fc->sequence_encoding2;
    SS    = NULL;
    S5    = NULL;
    S3    = NULL;
  } else {
    n_seq = fc->n_seq;
    S1    = NULL;
    S     = NULL;
    SS    = fc->S;
    S5    = fc->S5;
    S3    = fc->S3;
  }

  turn = pf_params->model_details.min_loop_size;

#ifndef VRNA_NR_SAMPLING_HASH
  if (current_node) {
    memorized_node_prev = NULL;
    memorized_node_cur  = (*current_node)->head;
  }

#endif

  qm1j = qm1[j];

  if (current_node) {
    fbd = NR_TOTAL_WEIGHT(*current_node) *
          qm1j /
          (*q_remain);
  }

  r = vrna_urn() * (qm1j - fbd);

  /* pointer magic for less instructions in loop */
  qb                -= j;
  hard_constraints  += n * j;

  for (qt = 0., l = 1; l + turn < j; l++) {
    if (hard_constraints[l] & VRNA_CONSTRAINT_CONTEXT_MB_LOOP_ENC) {
      u = l - l;
      if (hc_up_ml[1] >= u) {
        q_temp = qb[my_iindx[l]] *
                 expMLbase[u];

        if (fc->type == VRNA_FC_TYPE_SINGLE) {
          type    = vrna_get_ptype_md(S[l], S[j], md);
          q_temp  *= vrna_exp_E_multibranch_stem(type, S1[l - 1], S1[j + 1], pf_params);
        } else {
          for (s = 0; s < n_seq; s++) {
            type    = vrna_get_ptype_md(SS[s][l], SS[s][j], md);
            q_temp  *= vrna_exp_E_multibranch_stem(type, S5[s][l], S3[s][j], pf_params);
          }
        }

        if (sc_wrapper_ml->red_stem)
          q_temp *= sc_wrapper_ml->red_stem(1, j, l, j, sc_wrapper_ml);

        if (current_node) {
          fbds = NR_GET_WEIGHT(*current_node, memorized_node_cur, NRT_QM1_NEW_BRANCH, l, j) *
                 qm1j /
                 (*q_remain);
          qt += q_temp - fbds;
        } else {
          qt += q_temp;
        }

        if (qt >= r) {
          if (current_node) {
            *q_remain *= q_temp / qm1j;
#ifdef VRNA_NR_SAMPLING_HASH
            *current_node = add_if_nexists(NRT_QM1_NEW_BRANCH, i, l, *current_node, *q_remain);
#else
            *current_node = add_if_nexists_ll(memory_dat,
                                              NRT_QM1_NEW_BRANCH,
                                              l,
                                              j,
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
          advance_cursor(&memorized_node_prev, &memorized_node_cur, NRT_QM1_NEW_BRANCH, l, j);

#endif
      } else {
        l = j;
        break;
      }
    }
  }
  if (l + turn > j) {
    if (current_node) {
      return 0;
    } else {
      vrna_log_error("backtrack failed in qm1_new");
      return 0;
    }
  }

  /* continue backtracking the stem branching off the multiloop */
  return backtrack(l, j, pstruc, fc, sc_wrap, nr_mem);
}


PRIVATE int
backtrack_qm2(vrna_fold_compound_t            *fc,
              unsigned int                    i,
              unsigned int                    j,
              char                            *pstruc,
              struct sc_wrappers              *sc_wrap,
              struct vrna_pbacktrack_memory_s *nr_mem)
{
  unsigned char         *hard_constraints;
  short                 *S, *S1, **SS, **S5, **S3;
  unsigned int          k, n, n_seq, s, turn, type, *hc_up_ml;
  int                   *my_iindx, ret = 0;
  FLT_OR_DBL            qt, q_temp, qm2_rem, r, fbd, fbds, *qm2, *qb, *qm, *expMLbase;
  double                *q_remain;
  vrna_exp_param_t      *pf_params;
  vrna_md_t             *md;
  vrna_hc_t             *hc;
  struct sc_mb_exp_dat  *sc_wrapper;

  struct nr_memory      **memory_dat;
  NR_NODE               **current_node;

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

  fbd   = 0.;                             /* stores weight of forbidden terms for given q[ij] */
  fbds  = 0.;                             /* stores weight of forbidden term for given motif */

  n     = fc->length;
  n_seq = (fc->type == VRNA_FC_TYPE_COMPARATIVE) ? fc->n_seq : 1;
  S     = (fc->type == VRNA_FC_TYPE_COMPARATIVE) ? NULL : fc->sequence_encoding2;
  S1    = (fc->type == VRNA_FC_TYPE_COMPARATIVE) ? NULL : fc->sequence_encoding;
  SS    = (fc->type == VRNA_FC_TYPE_COMPARATIVE) ? fc->S : NULL;
  S5    = (fc->type == VRNA_FC_TYPE_COMPARATIVE) ? fc->S5 : NULL;
  S3    = (fc->type == VRNA_FC_TYPE_COMPARATIVE) ? fc->S3 : NULL;

  pf_params         = fc->exp_params;
  md                = &(pf_params->model_details);
  turn              = md->min_loop_size;
  my_iindx          = fc->iindx;
  qb                = fc->exp_matrices->qb;
  qm                = fc->exp_matrices->qm;
  qm2               = fc->exp_matrices->qm2_real;
  expMLbase         = fc->exp_matrices->expMLbase;
  hc                = fc->hc;
  hc_up_ml          = hc->up_ml;
  hard_constraints  = hc->mx;
  sc_wrapper        = &(sc_wrap->sc_wrapper_ml);

#ifndef VRNA_NR_SAMPLING_HASH
  if (current_node) {
    memorized_node_prev = NULL;
    memorized_node_cur  = (*current_node)->head;
  }

#endif

  /*
   *  split-away the last branch from the segment with at least two branches
   *
   *  we do that by first checking whether position j should be paired or unpaired.
   *  if it should pair, we find the pairing partner and continue with the two
   *  segments
   */
  if (i + 2 * turn + 2 >= j) {
    vrna_log_error("backtracking impossible for qm2_real[%u, %u]", i, j);
    return 0; /* error */
  }

#ifndef VRNA_NR_SAMPLING_HASH
  if (current_node) {
    memorized_node_prev = NULL;
    memorized_node_cur  = (*current_node)->head;
  }

#endif

  q_temp = 0.;

  do {
    if (hc_up_ml[j] > 0) {
      if (current_node) {
        fbd = NR_TOTAL_WEIGHT(*current_node) *
              qm2[my_iindx[i] - j] /
              (*q_remain);
#ifdef  USE_FLOAT_PF
        if (fabsf(NR_TOTAL_WEIGHT(*current_node) - (*q_remain)) / (*q_remain) <= FLT_EPSILON)
#else
        if (fabs(NR_TOTAL_WEIGHT(*current_node) - (*q_remain)) / (*q_remain) <= DBL_EPSILON)
#endif
          /* exhausted ensemble */
          return 0;
      }

      r       = vrna_urn() * (qm2[my_iindx[i] - j] - fbd);
      q_temp  = qm2[my_iindx[i] - j + 1] *
                expMLbase[1];

      if (sc_wrapper->red_ml)
        q_temp *= sc_wrapper->red_ml(i, j, i, j - 1, sc_wrapper);

      if (current_node) {
        fbds = NR_GET_WEIGHT(*current_node, memorized_node_cur, NRT_QM2_UNPAIR, j - 1, j) *
               qm2[my_iindx[i] - j] /
               (*q_remain);
        qt = q_temp - fbds;
      } else {
        qt = q_temp;
      }

      if (r > qt) {
        break; /* j is paired */
      } else if (current_node) {
        *q_remain *= q_temp / qm2[my_iindx[i] - j];
#ifdef VRNA_NR_SAMPLING_HASH
        *current_node = add_if_nexists(NRT_QM2_UNPAIR,
                                       j - 1,
                                       j,
                                       *current_node,
                                       *q_remain);
#else
        *current_node = add_if_nexists_ll(memory_dat,
                                          NRT_QM2_UNPAIR,
                                          j - 1,
                                          j,
                                          memorized_node_prev,
                                          memorized_node_cur,
                                          *current_node,
                                          *q_remain);
        reset_cursor(&memorized_node_prev, &memorized_node_cur, *current_node); /* resets cursor */
#endif
      }

      j--;
    } else {
      /* j must be paired, so continue by finding its pairing partner */
      break;
    }
  } while (i + 2 * turn + 2 < j);

  if (i + 2 * turn + 2 == j) {
    vrna_log_error("backtracking failed for qm2_real 1");
    return 0; /* error */
  }

#ifndef  VRNA_NR_SAMPLING_HASH
  if (current_node)
    advance_cursor(&memorized_node_prev, &memorized_node_cur, NRT_QM2_UNPAIR, j - 1, j);

#endif

  /* j is paired, let's find the pairing partner */
  qt      = 0.;
  qm2_rem = qm2[my_iindx[i] - j] - q_temp;

  if (current_node) {
    fbd = NR_TOTAL_WEIGHT_TYPE(NRT_QM2_BRANCH, *current_node) *
          qm2[my_iindx[i] - j] /
          (*q_remain);
  }

  r = vrna_urn() * (qm2_rem - fbd);

  for (k = i + turn + 2; k + turn < j; k++) {
    if (hard_constraints[n * j + k] & VRNA_CONSTRAINT_CONTEXT_MB_LOOP_ENC) {
      q_temp = qm[my_iindx[i] - k + 1] *
               qb[my_iindx[k] - j];

      if (fc->type == VRNA_FC_TYPE_SINGLE) {
        type    = vrna_get_ptype_md(S[k], S[j], md);
        q_temp  *= vrna_exp_E_multibranch_stem(type, S1[k - 1], S1[j + 1], pf_params);
      } else {
        for (s = 0; s < n_seq; s++) {
          type    = vrna_get_ptype_md(SS[s][k], SS[s][j], md);
          q_temp  *= vrna_exp_E_multibranch_stem(type, S5[s][k], S3[s][j], pf_params);
        }
      }

      if (sc_wrapper->decomp_ml)
        q_temp *= sc_wrapper->decomp_ml(i, j, k - 1, k, sc_wrapper);

      if (sc_wrapper->red_stem)
        q_temp *= sc_wrapper->decomp_ml(k, j, k, j, sc_wrapper);

      if (current_node) {
        fbds = NR_GET_WEIGHT(*current_node, memorized_node_cur, NRT_QM2_BRANCH, k, j) *
               qm2[my_iindx[i] - j] /
               (*q_remain);
        qt += q_temp - fbds;
      } else {
        qt += q_temp;
      }

      if (qt > r) {
        if (current_node) {
          *q_remain *= q_temp / qm2[my_iindx[i] - j];
#ifdef VRNA_NR_SAMPLING_HASH
          *current_node = add_if_nexists(NRT_QM2_BRANCH,
                                         k,
                                         j,
                                         *current_node,
                                         *q_remain);
#else
          *current_node = add_if_nexists_ll(memory_dat,
                                            NRT_QM2_BRANCH,
                                            k,
                                            j,
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
        advance_cursor(&memorized_node_prev,
                       &memorized_node_cur,
                       NRT_QM2_BRANCH,
                       k,
                       j);

#endif
    }
  }

  if (k + turn >= j) {
    vrna_log_error("backtracking failed for qm2_real 2");
    return 0;
  }

  /* backtrack further multiloop parts */
  ret = backtrack_qm(i, k - 1, pstruc, fc, sc_wrap, nr_mem);

  /* backtrack stem */
  ret &= backtrack(k, j, pstruc, fc, sc_wrap, nr_mem);

  return ret;
}


PRIVATE int
backtrack(unsigned int                    i,
          unsigned int                    j,
          char                            *pstruc,
          vrna_fold_compound_t            *vc,
          struct sc_wrappers              *sc_wrap,
          struct vrna_pbacktrack_memory_s *nr_mem)
{
  unsigned char         *hard_constraints, hc_decompose;
  char                  *ptype;
  short                 *S1, **S, **S5, **S3;
  unsigned int          k, l, u1, u2, max_k, min_l, **a2s, s, n_seq, n, type, type_2,
                        *types, u1_local, u2_local, *hc_up_int, turn;
  int                   *my_iindx, *jindx, ret, *pscore, *rtype, kl;
  FLT_OR_DBL            *qb, *scale, r, fbd, fbds, qbt1, qbr, q_temp,
                        kTn, expMLclosing, *qm2_real;
  double                *q_remain;
  vrna_mx_pf_t          *matrices;
  vrna_exp_param_t      *pf_params;
  vrna_md_t             *md;
  vrna_hc_t             *hc;

  struct nr_memory      **memory_dat;
  struct sc_int_exp_dat *sc_wrapper_int;
  struct sc_mb_exp_dat  *sc_wrapper_ml;

  NR_NODE               **current_node;

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
  type    = 0;

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
  qm2_real  = matrices->qm2_real;
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

    switch (vc->type) {
      case VRNA_FC_TYPE_COMPARATIVE:
        qbr /= exp(pscore[jindx[j] + i] / kTn);
        for (s = 0; s < n_seq; s++)
          types[s] = vrna_get_ptype_md(S[s][i], S[s][j], md);
        break;

      default:
        type = vrna_get_ptype(jindx[j] + i, ptype);
        break;
    }

    if (current_node)
      fbd = NR_TOTAL_WEIGHT(*current_node) * qbr / (*q_remain);

    pstruc[i - 1] = '(';
    pstruc[j - 1] = ')';

    r     = vrna_urn() * (qbr - fbd);
    qbt1  = 0.;

    hc_decompose = hard_constraints[n * i + j];

    /* hairpin contribution */
    q_temp = vrna_exp_eval_hairpin(vc, i, j, VRNA_EVAL_LOOP_DEFAULT);

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
      /* internal loop contributions */
      max_k = i + MAXLOOP + 1;
      max_k = MIN2(max_k, j - turn - 2);
      max_k = MIN2(max_k, i + 1 + hc_up_int[i + 1]);
      for (k = i + 1; k <= max_k; k++) {
        u1    = k - i - 1;
        min_l = k + turn + 1;
        if (min_l + MAXLOOP + 1 < j + u1)
          min_l = j - 1 - MAXLOOP + u1;

        kl = my_iindx[k] - j + 1;
        for (u2 = 0, l = j - 1; l >= min_l; l--, kl++, u2++) {
          if (hc_up_int[l + 1] < u2)
            break;

          if (hard_constraints[n * k + l] & VRNA_CONSTRAINT_CONTEXT_INT_LOOP_ENC) {
            q_temp = qb[kl]
                     * scale[u1 + u2 + 2];

            switch (vc->type) {
              case VRNA_FC_TYPE_COMPARATIVE:
                for (s = 0; s < n_seq; s++) {
                  u1_local  = a2s[s][k - 1] - a2s[s][i] /*??*/;
                  u2_local  = a2s[s][j - 1] - a2s[s][l];
                  type_2    = vrna_get_ptype_md(S[s][l], S[s][k], md);
                  q_temp    *= vrna_exp_E_internal(u1_local,
                                                   u2_local,
                                                   types[s],
                                                   type_2,
                                                   S3[s][i],
                                                   S5[s][j],
                                                   S5[s][k],
                                                   S3[s][l],
                                                   pf_params);
                }
                break;

              default:
                type_2 = rtype[vrna_get_ptype(jindx[l] + k, ptype)];

                /* add *scale[u1+u2+2] */
                q_temp *= vrna_exp_E_internal(u1,
                                              u2,
                                              type,
                                              type_2,
                                              S1[i + 1],
                                              S1[j - 1],
                                              S1[k - 1],
                                              S1[l + 1],
                                              pf_params);
                break;
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

        return backtrack(k, l, pstruc, vc, sc_wrap, nr_mem); /* found the internal loop, repeat for inside */
      } else {
        /* internal loop contributions did not exceed threshold, so we break */
        break;
      }
    } else {
      /* must not be internal loop, so we break out */
      break;
    }
  } while (1);

  /* backtrack in multi-loop */
  if (hard_constraints[n * j + i] & VRNA_CONSTRAINT_CONTEXT_MB_LOOP) {
    q_temp = qm2_real[my_iindx[i + 1] - j + 1] *
             expMLclosing *
             scale[2];

    if (vc->type == VRNA_FC_TYPE_SINGLE) {
      type    = rtype[vrna_get_ptype(jindx[j] + i, ptype)];
      q_temp  *= vrna_exp_E_multibranch_stem(type, S1[j - 1], S1[i + 1], pf_params);
    } else {
      for (s = 0; s < n_seq; s++) {
        type    = vrna_get_ptype_md(S[s][j], S[s][i], md);
        q_temp  *= vrna_exp_E_multibranch_stem(type, S5[s][j], S3[s][i], pf_params);
      }
    }

    if (sc_wrapper_ml->pair)
      q_temp *= sc_wrapper_ml->pair(i, j, sc_wrapper_ml);

    if (current_node) {
      fbds = NR_GET_WEIGHT(*current_node, memorized_node_cur, NRT_MB_LOOP, 0, 0) *
             qbr /
             (*q_remain);
      qbt1 += q_temp - fbds;
    } else {
      qbt1 += q_temp;
    }

    if (qbt1 < r) {
      vrna_log_debug("Backtracking failed for pair (%d,%d)", i, j);
      free(types);
      return 0;
    }

    /* must be multibranch loop, so update non-redundant stuff and proceed backtracking */
    if (current_node) {
      *q_remain *= q_temp / qbr;
#ifdef VRNA_NR_SAMPLING_HASH
      *current_node = add_if_nexists(NRT_MB_LOOP,
                                     0,
                                     0,
                                     *current_node,
                                     *q_remain);
#else
      *current_node = add_if_nexists_ll(memory_dat,
                                        NRT_MB_LOOP,
                                        0,
                                        0,
                                        memorized_node_prev,
                                        memorized_node_cur,
                                        *current_node,
                                        *q_remain);
#endif
    }

    ret = backtrack_qm2(vc, i + 1, j - 1, pstruc, sc_wrap, nr_mem);
  }

  free(types);

  return ret;
}


PRIVATE unsigned int
pbacktrack_circ(vrna_fold_compound_t  *vc,
                unsigned int          num_samples,
                vrna_bs_result_f      bs_cb,
                void                  *data)
{
  unsigned char         *hc_mx, eval_loop;
  char                  *pstruc;
  short                 *S1, *S2, **S, **S5, **S3;
  unsigned int          i, j, k, l, n, u, *hc_up, type, type2, *tt, s, n_seq, **a2s, u1_local,
                        u2_local, u3_local, count, turn, ln1, ln2, ln3, lstart;
  int                   *my_iindx;
  FLT_OR_DBL            r, qt, q_temp, qo, qmo, *scale, *qb, *qm2_real, *qm1_new,
                        qb_ij, expMLclosing;
  vrna_exp_param_t      *pf_params;
  vrna_md_t             *md;
  vrna_mx_pf_t          *matrices;
  struct sc_wrappers    *sc_wrap;
  struct sc_ext_exp_dat *sc_wrapper_ext;
  struct sc_int_exp_dat *sc_wrapper_int;
  struct sc_mb_exp_dat  *sc_wrapper_ml;

  n             = vc->length;
  pf_params     = vc->exp_params;
  md            = &(pf_params->model_details);
  matrices      = vc->exp_matrices;
  my_iindx      = vc->iindx;
  expMLclosing  = pf_params->expMLclosing;
  turn          = pf_params->model_details.min_loop_size;

  qo        = matrices->qo;
  qmo       = matrices->qmo;
  qb        = matrices->qb;
  qm1_new   = matrices->qm1_new;
  qm2_real  = matrices->qm2_real;
  scale     = matrices->scale;

  hc_mx = vc->hc->mx;
  hc_up = vc->hc->up_int;

  sc_wrap         = sc_init(vc);
  sc_wrapper_ext  = &(sc_wrap->sc_wrapper_ext);
  sc_wrapper_int  = &(sc_wrap->sc_wrapper_int);
  sc_wrapper_ml   = &(sc_wrap->sc_wrapper_ml);

  if (vc->type == VRNA_FC_TYPE_SINGLE) {
    n_seq         = 1;
    tt            = NULL;
    type          = 0;
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
    type          = 0;
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

    qt = scale[n];

    switch (vc->type) {
      case VRNA_FC_TYPE_COMPARATIVE:
        for (s = 0; s < n_seq; s++)
          qt *= pow(vrna_exp_E_exterior_loop(a2s[s][n], md), 1. / (double)n_seq);
        break;

      default:
        qt *= vrna_exp_E_exterior_loop(n, md);
        break;
    }

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
              vrna_exp_eval_hairpin(vc, j, i, VRNA_EVAL_LOOP_DEFAULT);

        /* found a hairpin? so backtrack in the enclosed part and we're done  */
        if (qt > r) {
          backtrack(i, j, pstruc, vc, sc_wrap, NULL);
          goto pbacktrack_circ_loop_end;
        }

        /* 2. search for (k,l) with which we can close an internal loop  */
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
                  case VRNA_FC_TYPE_COMPARATIVE:
                    for (s = 0; s < n_seq; s++) {
                      type2     = vrna_get_ptype_md(S[s][l], S[s][k], md);
                      u1_local  = a2s[s][i - 1];
                      u2_local  = a2s[s][k - 1] - a2s[s][j];
                      u3_local  = a2s[s][n] - a2s[s][l];
                      q_temp    *= vrna_exp_E_internal(u1_local + u3_local,
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

                  default:
                    type2   = vrna_get_ptype_md(S2[l], S2[k], md);
                    q_temp  *= vrna_exp_E_internal(ln2 + ln3,
                                                   ln1,
                                                   type2,
                                                   type,
                                                   S1[l + 1],
                                                   S1[k - 1],
                                                   S1[i - 1],
                                                   S1[j + 1],
                                                   pf_params);
                    break;
                }

                if (sc_wrapper_int->pair_ext)
                  q_temp *= sc_wrapper_int->pair_ext(i, j, k, l, sc_wrapper_int);

                qt += q_temp;
                /*
                 * found an exterior internal loop? also this time, we can go straight
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

    /* find split-point between qm1 and qm2 */
    qt  = 0.;
    r   = vrna_urn() * qmo;
    if (sc_wrapper_ml->decomp_ml) {
      for (k = turn + 1; k + 2 * turn + 3 < n; k++) {
        qt += qm1_new[k] *
              qm2_real[my_iindx[k + 1] - n] *
              expMLclosing *
              sc_wrapper_ml->
              decomp_ml(1,
                        n,
                        k,
                        k + 1,
                        sc_wrapper_ml);


        if (qt > r) {
          backtrack_qm1(vc, k, pstruc, sc_wrap, NULL);
          backtrack_qm2(vc, k, n, pstruc, sc_wrap, NULL);
          goto pbacktrack_circ_loop_end;
        }
      }
    } else {
      for (k = turn + 1; k + 2 * turn + 3 < n; k++) {
        qt += qm1_new[k] *
              qm2_real[my_iindx[k + 1] - n] *
              expMLclosing;

        if (qt > r) {
          backtrack_qm1(vc, k, pstruc, sc_wrap, NULL);
          backtrack_qm2(vc, k + 1, n, pstruc, sc_wrap, NULL);
          goto pbacktrack_circ_loop_end;
        }
      }
    }

    /*
     * if we reach the actual end of this function, an error has occured
     * cause we HAVE TO find an exterior loop or an open chain!!!
     */
    vrna_log_error("backtracking failed in exterior loop");

pbacktrack_circ_loop_end:

    if (bs_cb)
      bs_cb(pstruc, data);

    free(pstruc);
  }

  sc_free(sc_wrap);

  return count;
}
