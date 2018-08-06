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


struct nr_structure_list {
  unsigned int  num;
  char          **list;
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
  "Activate unique multiloop decomposition by setting the"
  " uniq_ML field of the model details structure to a non-zero"
  " value before running vrna_pf()!";


/*
 #################################
 # PRIVATE FUNCTION DECLARATIONS #
 #################################
 */


PRIVATE void
save_nr_samples(const char  *structure,
                void        *data);


/* In the following:
 * - den is a pointer to value of sum of Boltzmann factors of still accessible solutions at that point
 * - current_node is a double pointer to current node in datastructure memorizing the solutions and paths taken */
PRIVATE char *
pbacktrack5_gen(vrna_fold_compound_t  *vc,
                int                   length,
                double                *den,
                NR_NODE               **current_node,
                struct nr_memory      *memory_dat);


PRIVATE int
backtrack(int                   i,
          int                   j,
          char                  *pstruc,
          vrna_fold_compound_t  *vc,
          double                *den,
          NR_NODE               **current_node,
          struct nr_memory      *memory_dat);


PRIVATE int
backtrack_ext_loop(int                  init_val,
                   char                 *pstruc,
                   vrna_fold_compound_t *vc,
                   int                  length,
                   double               *den,
                   NR_NODE              **current_node,
                   struct nr_memory     *memory_dat);


PRIVATE int
backtrack_qm(int                  i,
             int                  j,
             char                 *pstruc,
             vrna_fold_compound_t *vc);


PRIVATE int
backtrack_qm_nr(int                   i,
                int                   j,
                char                  *pstruc,
                vrna_fold_compound_t  *vc,
                double                *den,
                NR_NODE               **current_node,
                struct nr_memory      *memory_dat);


PRIVATE int
backtrack_qm1(int                   i,
              int                   j,
              char                  *pstruc,
              vrna_fold_compound_t  *vc,
              double                *den,
              NR_NODE               **current_node,
              struct nr_memory      *memory_dat);


PRIVATE void
backtrack_qm2(int                   u,
              int                   n,
              char                  *pstruc,
              vrna_fold_compound_t  *vc);


PRIVATE char *
wrap_pbacktrack_circ(vrna_fold_compound_t *vc);


PRIVATE void
backtrack_comparative(vrna_fold_compound_t  *vc,
                      char                  *pstruc,
                      int                   i,
                      int                   j,
                      double                *prob);


PRIVATE void
backtrack_qm1_comparative(vrna_fold_compound_t  *vc,
                          char                  *pstruc,
                          int                   i,
                          int                   j,
                          double                *prob);


/*
 *  @brief Sample a consensus secondary structure from the Boltzmann ensemble according its probability
 *
 *  @ingroup consensus_stochbt
 *
 *  @see vrna_pf() for precomputing the partition function matrices, and
 *
 *  @param  vc    The #vrna_fold_compound_t of type #VRNA_FC_TYPE_COMPARATIVE with precomputed partition function matrices
 *  @param  prob  to be described (berni)
 *  @return       A sampled consensus secondary structure in dot-bracket notation
 */
PRIVATE char *
pbacktrack_comparative(vrna_fold_compound_t *vc,
                       double               *prob);


/*
 #################################
 # BEGIN OF FUNCTION DEFINITIONS #
 #################################
 */

/*
 * stochastic backtracking in pf_fold arrays
 * returns random structure S with Boltzman probabilty
 * p(S) = exp(-E(S)/kT)/Z
 */
PUBLIC char *
vrna_pbacktrack(vrna_fold_compound_t *vc)
{
  char    *structure  = NULL;
  double  prob        = 1.;

  if (vc) {
    if (!vc->exp_params) {
      vrna_message_warning("vrna_pbacktrack: DP matrices are missing! Call vrna_pf() first!");
      return NULL;
    } else if (!vc->exp_params->model_details.uniq_ML) {
      vrna_message_warning("vrna_pbacktrack: Unique multiloop decomposition is unset!");
      vrna_message_info(stderr, info_set_uniq_ml);
      return NULL;
    }

    switch (vc->type) {
      case VRNA_FC_TYPE_SINGLE:
        if (vc->exp_params->model_details.circ)
          return wrap_pbacktrack_circ(vc);
        else
          return vrna_pbacktrack5(vc, vc->length);

        break;

      case VRNA_FC_TYPE_COMPARATIVE:
        return pbacktrack_comparative(vc, &prob);
        break;

      default:
        vrna_message_warning("unrecognized fold compound type");
        return structure;
        break;
    }
  }

  return structure;
}


/* adapter function for more general expression using non-redundant sampling */
PUBLIC char *
vrna_pbacktrack5(vrna_fold_compound_t *vc,
                 int                  length)
{
  return pbacktrack5_gen(vc, length, NULL, NULL, NULL);
}


PUBLIC char **
vrna_pbacktrack_nr(vrna_fold_compound_t *vc,
                   int                  num_samples)
{
  struct nr_structure_list data;

  data.num      = 0;
  data.list     = (char **)vrna_alloc(sizeof(char *) * num_samples);
  data.list[0]  = NULL;
  vrna_pbacktrack_nr_cb(vc, num_samples, &save_nr_samples, (void *)&data);

  return data.list;
}


PUBLIC void
vrna_pbacktrack_nr_cb(vrna_fold_compound_t              *vc,
                      int                               num_samples,
                      vrna_boltzmann_sampling_callback  *bs_cb,
                      void                              *data)
{
  if (vc) {
    if (!vc->exp_params) {
      vrna_message_warning("vrna_pbacktrack_nr_cb: DP matrices are missing! Call vrna_pf() first!");
    } else if (!vc->exp_params->model_details.uniq_ML) {
      vrna_message_warning("vrna_pbacktrack_nr_cb: Unique multiloop decomposition is unset!");
      vrna_message_info(stderr, info_set_uniq_ml);
    } else if (vc->type != VRNA_FC_TYPE_SINGLE) {
      vrna_message_warning(
        "vrna_pbacktrack_nr_cb: No implementation for comparative structure prediction available yet!");
    } else if (vc->exp_params->model_details.circ) {
      vrna_message_warning(
        "vrna_pbacktrack_nr_cb: No implementation for circular RNAs available yet!");
    } else {
      int    i;
      double den, part_fci;
      struct nr_memory  memory_dat  = {
        NULL, 0
      };
      NR_NODE *current_node;

      den       = 0;
      part_fci  = vc->exp_matrices->q[vc->iindx[1] - vc->length];

#ifdef VRNA_NR_SAMPLING_HASH
      current_node = create_root(vc->length);
#else
      memory_dat.nr_memory_allocated = vrna_alloc(vc->length * num_samples * sizeof(NR_NODE));  // memory pre-allocation
      current_node = create_ll_root(&memory_dat);
#endif


      for (i = 0; i < num_samples; i++) {
        den = vc->exp_matrices->q[vc->iindx[1] - vc->length];
        char *ss = pbacktrack5_gen(vc, vc->length, &den, &current_node, &memory_dat);
        bs_cb(ss, data);
        free(ss);
        /* finish if no more structures available */
        if (!ss)
          break;

#ifdef VRNA_NR_SAMPLING_HASH
        current_node = traceback_to_root(current_node, den);
#else
        current_node = traceback_to_ll_root(current_node, den);
#endif
      }

      /* print warning if we've aborted backtracking too early */
      if ((i > 0) && (i < num_samples)) {
        den = vc->exp_matrices->q[vc->iindx[1] - vc->length];
#ifdef VRNA_NR_SAMPLING_HASH
        current_node = traceback_to_root(current_node, den);
#else
        current_node = traceback_to_ll_root(current_node, den);
#endif
        vrna_message_warning("vrna_pbacktrack_nr*(): Stopped backtracking after %d samples due to numeric instabilities!\n"
                             "Coverage of partition function so far: %f%%",
                             i,
                             100. * current_node->weight / part_fci);
      }

#ifdef VRNA_NR_SAMPLING_HASH
      free_all_nr(current_node);
#else
      free(memory_dat.nr_memory_allocated);
#endif
    }
  }
}


/* general expr of vrna5_pbacktrack with possibility of non-redundant sampling */
PRIVATE char *
pbacktrack5_gen(vrna_fold_compound_t  *vc,
                int                   length,
                double                *den,
                NR_NODE               **current_node,
                struct nr_memory      *memory_dat)
{
  int   ret, i, j, n, start;
  char  *pstruc;

  n       = vc->length;
  pstruc  = vrna_alloc((length + 1) * sizeof(char));

  for (i = 0; i < length; i++)
    pstruc[i] = '.';


#ifdef VRNA_WITH_BOUSTROPHEDON
  ret = backtrack_ext_loop(length, pstruc, vc, length, den, current_node, memory_dat);
#else
  ret = backtrack_ext_loop(1, pstruc, vc, length, den, current_node, memory_dat);
#endif

  if (ret > 0)
    return pstruc;

  free(pstruc);
  return NULL;
}


/* backtrack one external */
PRIVATE int
backtrack_ext_loop(int                  init_val,
                   char                 *pstruc,
                   vrna_fold_compound_t *vc,
                   int                  length,
                   double               *den,
                   NR_NODE              **current_node,
                   struct nr_memory     *memory_dat)
{
  FLT_OR_DBL        r, fbd, fbds, concf, qt, q_temp, qkl;
  int               ret, i, j, ij, n, k, u, start, type;
  int               *my_iindx, *jindx, hc_decompose, *hc_up_ext;
  FLT_OR_DBL        *q, *qb, *q1k, *qln, *scale;
  unsigned char     *hard_constraints;
  short             *S1, *S2;
  vrna_mx_pf_t      *matrices;
  vrna_md_t         *md;
  vrna_hc_t         *hc;
  vrna_sc_t         *sc;
  vrna_exp_param_t  *pf_params;
  /* non-redundant data-structure memorization nodes */
  NR_NODE           *memorized_node_prev; /* remembers previous-to-current node in linked list */
  NR_NODE           *memorized_node_cur;  /* remembers actual node in linked list */

  fbd   = 0.;                             /* stores weight of forbidden terms for given q[ij]*/
  fbds  = 0.;                             /* stores weight of forbidden term for given motif */

  n = vc->length;

  pf_params = vc->exp_params;
  md        = &(vc->exp_params->model_details);
  my_iindx  = vc->iindx;
  matrices  = vc->exp_matrices;

  hc  = vc->hc;
  sc  = vc->sc;
  S1  = vc->sequence_encoding;
  S2  = vc->sequence_encoding2;

  hard_constraints  = hc->mx;
  hc_up_ext         = hc->up_ext;

  if (length > n) {
    vrna_message_warning("vrna_pbacktrack5: 3'-end exceeds sequence length");
    return -1;
  } else if (length < 1) {
    vrna_message_warning("vrna_pbacktrack5: 3'-end too small");
    return -1;
  } else if ((!matrices) || (!matrices->q) || (!matrices->qb) || (!matrices->qm) || (!pf_params)) {
    vrna_message_warning("vrna_pbacktrack5: DP matrices are missing! Call vrna_pf() first!");
    return -1;
  } else if ((!vc->exp_params->model_details.uniq_ML) || (!matrices->qm1)) {
    vrna_message_warning("vrna_pbacktrack5: Unique multiloop decomposition is unset!");
    vrna_message_info(stderr, info_set_uniq_ml);
    return -1;
  }

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
                (*den);

#ifdef  USE_FLOAT_PF
          if (fabsf(NR_TOTAL_WEIGHT(*current_node) - (*den)) / (*den) <= FLT_EPSILON)
#else
          if (fabs(NR_TOTAL_WEIGHT(*current_node) - (*den)) / (*den) <= DBL_EPSILON)
#endif
            /* exhausted ensemble */
            return 0;
        }

        r       = vrna_urn() * (q1k[j] - fbd);
        q_temp  = q1k[j - 1] * scale[1];

        if (sc) {
          if (sc->exp_energy_up)
            q_temp *= sc->exp_energy_up[j][1];

          if (sc->exp_f)
            q_temp *= sc->exp_f(1, j, 1, j - 1, VRNA_DECOMP_EXT_EXT, sc->data);
        }

        if (current_node) {
          fbds = NR_GET_WEIGHT(*current_node, memorized_node_cur, NRT_UNPAIRED_SG, j - 1, j) *
                 q1k[j] /
                 (*den);
        }

        if (r > (q_temp - fbds)) {
          break;                /* j is paired */
        } else if (current_node) {
          /* j is unpaired */
          *den *= q_temp / q1k[j];
#ifdef VRNA_NR_SAMPLING_HASH
          *current_node = add_if_nexists(NRT_UNPAIRED_SG, j - 1, j, *current_node);
#else
          *current_node = add_if_nexists_ll(memory_dat,
                                            NRT_UNPAIRED_SG,
                                            j - 1,
                                            j,
                                            memorized_node_prev,
                                            memorized_node_cur,
                                            *current_node);
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
            (*den);
    }

    r = vrna_urn() * (q1k[j] - q_temp - fbd);
    u = j - 1;

    for (qt = 0, k = 1; k < j; k++) {
      /* apply alternating boustrophedon scheme to variable i */
      i = (int)(1 + (u - 1) * ((k - 1) % 2)) +
          (int)((1 - (2 * ((k - 1) % 2))) * ((k - 1) / 2));
      ij            = my_iindx[i] - j;
      hc_decompose  = hard_constraints[n * j + i];
      if (hc_decompose & VRNA_CONSTRAINT_CONTEXT_EXT_LOOP) {
        type  = vrna_get_ptype_md(S2[i], S2[j], md);
        qkl   = qb[ij] * exp_E_ExtLoop(type,
                                       (i > 1) ? S1[i - 1] : -1,
                                       (j < n) ? S1[j + 1] : -1,
                                       pf_params);

        if (i > 1) {
          qkl *= q1k[i - 1];
          if (sc)
            if (sc->exp_f)
              qkl *= sc->exp_f(1, j, i - 1, i, VRNA_DECOMP_EXT_EXT_STEM, sc->data);
        } else {
          if (sc)
            if (sc->exp_f)
              qkl *= sc->exp_f(i, j, i, j, VRNA_DECOMP_EXT_STEM, sc->data);
        }

        if (current_node) {
          fbds = NR_GET_WEIGHT(*current_node, memorized_node_cur, NRT_EXT_LOOP, i, j) *
                 q1k[j] /
                 (*den);
          qt += qkl - fbds;
        } else {
          qt += qkl;
        }

        if (qt > r) {
          if (current_node) {
            *den *= qkl / q1k[j];
#ifdef VRNA_NR_SAMPLING_HASH
            *current_node = add_if_nexists(NRT_EXT_LOOP, i, j, *current_node);
#else
            *current_node = add_if_nexists_ll(memory_dat,
                                              NRT_EXT_LOOP,
                                              i,
                                              j,
                                              memorized_node_prev,
                                              memorized_node_cur,
                                              *current_node);
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

    backtrack(i, j, pstruc, vc, den, current_node, memory_dat);
    j   = i - 1;
    ret = backtrack_ext_loop(j, pstruc, vc, length, den, current_node, memory_dat);
  }

#else
  start = init_val;
  if (start < length) {
    /* find i position of first pair */
    for (i = start; i < length; i++) {
      if (hc_up_ext[i]) {
        if (current_node) {
          fbd = NR_TOTAL_WEIGHT(*current_node) *
                qln[i] /
                (*den);
        }

        r       = vrna_urn() * (qln[i] - fbd);
        q_temp  = qln[i + 1] * scale[1];

        if (sc) {
          if (sc->exp_energy_up)
            q_temp *= sc->exp_energy_up[i][1];

          if (sc->exp_f)
            q_temp *= sc->exp_f(i, length, i + 1, length, VRNA_DECOMP_EXT_EXT, sc->data);
        }

        if (current_node) {
          fbds = NR_GET_WEIGHT(*current_node, memorized_node_cur, NRT_UNPAIRED_SG, i, i + 1) *
                 qln[i] /
                 (*den);
        }

        if (r > (q_temp - fbds)) {
          break;                /* i is paired */
        } else if (current_node) {
          *den *= q_temp / qln[i];
#ifdef VRNA_NR_SAMPLING_HASH
          *current_node = add_if_nexists(NRT_UNPAIRED_SG, i, i + 1, *current_node);
#else
          *current_node = add_if_nexists_ll(memory_dat,
                                            NRT_UNPAIRED_SG,
                                            i,
                                            i + 1,
                                            memorized_node_prev,
                                            memorized_node_cur,
                                            *current_node);
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
            (*den);
    }

    r = vrna_urn() * (qln[i] - q_temp - fbd);
    for (qt = 0, j = i + 1; j <= length; j++) {
      ij            = my_iindx[i] - j;
      type          = vrna_get_ptype(jindx[j] + i, ptype);
      hc_decompose  = hard_constraints[n * i + j];
      if (hc_decompose & VRNA_CONSTRAINT_CONTEXT_EXT_LOOP) {
        qkl = qb[ij] * exp_E_ExtLoop(type,
                                     (i > 1) ? S1[i - 1] : -1,
                                     (j < n) ? S1[j + 1] : -1,
                                     pf_params);

        if (j < length) {
          qkl *= qln[j + 1];
          if (sc)
            if (sc->exp_f)
              qkl *= sc->exp_f(i, length, j, j + 1, VRNA_DECOMP_EXT_STEM_EXT, sc->data);
        } else {
          if (sc)
            if (sc->exp_f)
              qkl *= sc->exp_f(i, j, i, j, VRNA_DECOMP_EXT_STEM, sc->data);
        }

        if (current_node) {
          fbds = NR_GET_WEIGHT(*current_node, memorized_node_cur, NRT_EXT_LOOP, i, j) *
                 qln[i] /
                 (*den);
          qt += qkl - fbds;
        } else {
          qt += qkl;
        }

        if (qt > r) {
          if (current_node) {
            *den *= qkl / qln[i];
#ifdef VRNA_NR_SAMPLING_HASH
            *current_node = add_if_nexists(NRT_EXT_LOOP, i, j, *current_node);
#else
            *current_node = add_if_nexists_ll(memory_dat,
                                              NRT_EXT_LOOP,
                                              i,
                                              j,
                                              memorized_node_prev,
                                              memorized_node_cur,
                                              *current_node);
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
    backtrack(i, j, pstruc, vc, den, current_node, memory_dat);

    ret = backtrack_ext_loop(start, pstruc, vc, length, den, current_node, memory_dat);
  }

#endif

  return ret;
}


PRIVATE void
save_nr_samples(const char  *structure,
                void        *data)
{
  struct nr_structure_list *d = (struct nr_structure_list *)data;

  if (structure)
    d->list[d->num++] = strdup(structure);
  else
    d->list[d->num++] = NULL;
}


PRIVATE int
backtrack_qm(int                  i,
             int                  j,
             char                 *pstruc,
             vrna_fold_compound_t *vc)
{
  /* divide multiloop into qm and qm1  */
  FLT_OR_DBL    qmt, r, q_temp;
  int           k, u, cnt, span, turn;
  FLT_OR_DBL    *qm, *qm1, *expMLbase;
  int           *my_iindx, *jindx, *hc_up_ml, ret;
  vrna_sc_t     *sc;
  vrna_hc_t     *hc;

  vrna_mx_pf_t  *matrices = vc->exp_matrices;

  ret       = 1;
  my_iindx  = vc->iindx;
  jindx     = vc->jindx;

  hc        = vc->hc;
  sc        = vc->sc;
  hc_up_ml  = hc->up_ml;

  qm        = matrices->qm;
  qm1       = matrices->qm1;
  expMLbase = matrices->expMLbase;

  turn = vc->exp_params->model_details.min_loop_size;

  while (j > i) {
    /* now backtrack  [i ... j] in qm[] */
    r   = vrna_urn() * qm[my_iindx[i] - j];
    qmt = qm1[jindx[j] + i];
    k   = cnt = i;
    if (qmt < r) {
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

          if (sc) {
            if (sc->exp_energy_up)
              q_temp *= sc->exp_energy_up[i][u];

            if (sc->exp_f)
              q_temp *= sc->exp_f(i, j, k, j, VRNA_DECOMP_ML_ML, sc->data);
          }

          qmt += q_temp;
        }

        /* split between k-1, k */
        q_temp = qm[my_iindx[i] - (k - 1)] * qm1[jindx[j] + k];

        if (sc)
          if (sc->exp_f)
            q_temp *= sc->exp_f(i, j, k - 1, k, VRNA_DECOMP_ML_ML_ML, sc->data);

        qmt += q_temp;

        if (qmt >= r)
          break;
      }
    }

    if (cnt > j) {
      vrna_message_error("backtrack failed in qm");
      return 0;
    }

    ret = backtrack_qm1(k, j, pstruc, vc, NULL, NULL, NULL);

    if (ret == 0)
      return ret;

    if (k < i + turn)
      break;            /* no more pairs */

    u = k - i;
    /* check whether we make the decision to leave [i..k-1] unpaired */
    if (hc_up_ml[i] >= u) {
      q_temp = expMLbase[u];

      if (sc) {
        if (sc->exp_energy_up)
          q_temp *= sc->exp_energy_up[i][u];

        if (sc->exp_f)
          q_temp *= sc->exp_f(i, k - 1, i, k - 1, VRNA_DECOMP_ML_UP, sc->data);
      }

      r = vrna_urn() * (qm[my_iindx[i] - (k - 1)] + q_temp);
      if (q_temp >= r)
        break;
    }

    j = k - 1;
  }

  return ret;
}


/* non redundant version of function bactrack_qm */
PRIVATE int
backtrack_qm_nr(int                   i,
                int                   j,
                char                  *pstruc,
                vrna_fold_compound_t  *vc,
                double                *den,
                NR_NODE               **current_node,
                struct nr_memory      *memory_dat)
{
  /* divide multiloop into qm and qm1  */
  FLT_OR_DBL  qmt, fbd, fbds, r, q_temp, q_up, q_pair;
  int         k, n, u, cnt, span, turn;
  int         is_unpaired; /* 1 if [i ... k-1] is unpaired */
  FLT_OR_DBL  *qm, *qm1, *expMLbase;
  int         *my_iindx, *jindx, *hc_up_ml, ret;
  vrna_sc_t   *sc;
  vrna_hc_t   *hc;
  /* non-redundant data-structure memorization nodes */
  NR_NODE     *memorized_node_prev; /* remembers previous-to-current node in linked list */
  NR_NODE     *memorized_node_cur;  /* remembers actual node in linked list */

  ret   = 1;
  fbd   = 0.;                       /* stores weight of forbidden terms for given q[ij]*/
  fbds  = 0.;                       /* stores weight of forbidden term for given motif */

  is_unpaired = 0;

  n = j;
  vrna_mx_pf_t *matrices = vc->exp_matrices;

  my_iindx  = vc->iindx;
  jindx     = vc->jindx;

  hc        = vc->hc;
  sc        = vc->sc;
  hc_up_ml  = hc->up_ml;

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
            (*den);
    }

    r = vrna_urn() * (qm[my_iindx[i] - j] - fbd);
    if (current_node) {
      fbds = NR_GET_WEIGHT(*current_node, memorized_node_cur, NRT_QM_UNPAIR, i, 0) *
             qm[my_iindx[i] - j] /
             (*den);
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

          if (sc) {
            if (sc->exp_energy_up)
              q_temp *= sc->exp_energy_up[i][u];

            if (sc->exp_f)
              q_temp *= sc->exp_f(i, j, k, j, VRNA_DECOMP_ML_ML, sc->data);
          }

          if (current_node) {
            fbds = NR_GET_WEIGHT(*current_node, memorized_node_cur, NRT_QM_UNPAIR, k, 0) *
                   qm[my_iindx[i] - j] /
                   (*den);
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

        if (sc)
          if (sc->exp_f)
            q_temp *= sc->exp_f(i, j, k - 1, k, VRNA_DECOMP_ML_ML_ML, sc->data);

        if (current_node) {
          fbds = NR_GET_WEIGHT(*current_node, memorized_node_cur, NRT_QM_PAIR, k, 0) *
                 qm[my_iindx[i] - j] /
                 (*den);
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
      *den *= q_temp / qm[my_iindx[i] - j];
#ifdef VRNA_NR_SAMPLING_HASH
      if (is_unpaired)
        *current_node = add_if_nexists(NRT_QM_UNPAIR, k, 0, *current_node);
      else
        *current_node = add_if_nexists(NRT_QM_PAIR, k, 0, *current_node);

#else
      if (is_unpaired)
        *current_node = add_if_nexists_ll(memory_dat,
                                          NRT_QM_UNPAIR,
                                          k,
                                          0,
                                          memorized_node_prev,
                                          memorized_node_cur,
                                          *current_node);
      else
        *current_node = add_if_nexists_ll(memory_dat,
                                          NRT_QM_PAIR,
                                          k,
                                          0,
                                          memorized_node_prev,
                                          memorized_node_cur,
                                          *current_node);

#endif
    }

    if (cnt > j) {
      /* vrna_message_error("backtrack failed in qm"); */
      return 0;
    }

    ret = backtrack_qm1(k, j, pstruc, vc, den, current_node, memory_dat);

    if (ret == 0)
      return ret;

    if (k < i + turn)
      return ret;         /* no more pairs */

    if (!is_unpaired) {/* if we've chosen creating a branch in [i..k-1] */
      ret = backtrack_qm_nr(i, k - 1, pstruc, vc, den, current_node, memory_dat);

      if (ret == 0)
        return ret;
    }
  }

  return ret;
}


PRIVATE int
backtrack_qm1(int                   i,
              int                   j,
              char                  *pstruc,
              vrna_fold_compound_t  *vc,
              double                *den,
              NR_NODE               **current_node,
              struct nr_memory      *memory_dat)
{
  /* i is paired to l, i<l<j; backtrack in qm1 to find l */
  unsigned int      n;
  int               ii, l, il, type, turn, ret;
  FLT_OR_DBL        qt, fbd, fbds, r, q_temp;
  FLT_OR_DBL        *qm1, *qb, *expMLbase;
  vrna_mx_pf_t      *matrices;
  int               u, *my_iindx, *jindx, *hc_up_ml;
  char              *ptype;
  unsigned char     *hard_constraints;
  short             *S1;
  vrna_sc_t         *sc;
  vrna_hc_t         *hc;
  vrna_exp_param_t  *pf_params;

  /* non-redundant data-structure memorization nodes */
  NR_NODE           *memorized_node_prev; /* remembers previous-to-current node in linked list */
  NR_NODE           *memorized_node_cur;  /* remembers actual node in linked list */

  ret       = 1;
  n         = vc->length;
  fbd       = 0.;
  fbds      = 0.;
  pf_params = vc->exp_params;
  my_iindx  = vc->iindx;
  jindx     = vc->jindx;
  ptype     = vc->ptype;

  sc                = vc->sc;
  hc                = vc->hc;
  hc_up_ml          = hc->up_ml;
  hard_constraints  = hc->mx;

  matrices  = vc->exp_matrices;
  qb        = matrices->qb;
  qm1       = matrices->qm1;
  expMLbase = matrices->expMLbase;
  S1        = vc->sequence_encoding;

  turn = pf_params->model_details.min_loop_size;

#ifndef VRNA_NR_SAMPLING_HASH
  if (current_node) {
    memorized_node_prev = NULL;
    memorized_node_cur  = (*current_node)->head;
  }

#endif

  n = j;
  if (current_node) {
    fbd = NR_TOTAL_WEIGHT(*current_node) *
          qm1[jindx[j] + i] /
          (*den);
  }

  r   = vrna_urn() * (qm1[jindx[j] + i] - fbd);
  ii  = my_iindx[i];
  for (qt = 0., l = j; l > i + turn; l--) {
    il = jindx[l] + i;
    if (hard_constraints[n * i + l] & VRNA_CONSTRAINT_CONTEXT_MB_LOOP_ENC) {
      u = j - l;
      if (hc_up_ml[l + 1] >= u) {
        type    = vrna_get_ptype(il, ptype);
        q_temp  = qb[ii - l]
                  * exp_E_MLstem(type, S1[i - 1], S1[l + 1], pf_params)
                  * expMLbase[j - l];

        if (sc) {
          if (sc->exp_energy_up)
            q_temp *= sc->exp_energy_up[l + 1][j - l];

          if (sc->exp_f)
            q_temp *= sc->exp_f(i, j, i, l, VRNA_DECOMP_ML_STEM, sc->data);
        }

        if (current_node) {
          fbds = NR_GET_WEIGHT(*current_node, memorized_node_cur, NRT_QM1_BRANCH, i, l) *
                 qm1[jindx[j] + i] /
                 (*den);
          qt += q_temp - fbds;
        } else {
          qt += q_temp;
        }

        if (qt >= r) {
          if (current_node) {
            *den *= q_temp / qm1[jindx[j] + i];
#ifdef VRNA_NR_SAMPLING_HASH
            *current_node = add_if_nexists(NRT_QM1_BRANCH, i, l, *current_node);
#else
            *current_node = add_if_nexists_ll(memory_dat,
                                              NRT_QM1_BRANCH,
                                              i,
                                              l,
                                              memorized_node_prev,
                                              memorized_node_cur,
                                              *current_node);
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

  return backtrack(i, l, pstruc, vc, den, current_node, memory_dat);
}


PRIVATE void
backtrack_qm2(int                   k,
              int                   n,
              char                  *pstruc,
              vrna_fold_compound_t  *vc)
{
  FLT_OR_DBL  qom2t, r;
  int         u, turn;
  FLT_OR_DBL  *qm1, *qm2;
  int         *jindx;
  vrna_sc_t   *sc;

  jindx = vc->jindx;
  qm1   = vc->exp_matrices->qm1;
  qm2   = vc->exp_matrices->qm2;
  turn  = vc->exp_params->model_details.min_loop_size;
  sc    = vc->sc;

  r = vrna_urn() * qm2[k];
  /* we have to search for our barrier u between qm1 and qm1  */
  if ((sc) && (sc->exp_f)) {
    for (qom2t = 0., u = k + turn + 1; u < n - turn - 1; u++) {
      qom2t +=  qm1[jindx[u] + k] *
                qm1[jindx[n] + (u + 1)] *
                sc->exp_f(k, n, u, u + 1, VRNA_DECOMP_ML_ML_ML, sc->data);

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

  backtrack_qm1(k, u, pstruc, vc, NULL, NULL, NULL);
  backtrack_qm1(u + 1, n, pstruc, vc, NULL, NULL, NULL);
}


PRIVATE int
backtrack(int                   i,
          int                   j,
          char                  *pstruc,
          vrna_fold_compound_t  *vc,
          double                *den,
          NR_NODE               **current_node,
          struct nr_memory      *memory_dat)
{
  char              *ptype;
  unsigned char     *hard_constraints, hc_decompose;
  vrna_exp_param_t  *pf_params;
  FLT_OR_DBL        *qb, *qm, *qm1, *scale, tmp;
  FLT_OR_DBL        r, fbd, fbds, qbt1, qbr, qt, q_temp; /* qbr stores qb used for generating r */
  vrna_mx_pf_t      *matrices;
  unsigned int      n;
  int               *my_iindx, *jindx, *hc_up_int, ret;
  vrna_sc_t         *sc;
  vrna_hc_t         *hc;
  short             *S1;
  /* non-redundant data-structure memorization nodes */
  NR_NODE           *memorized_node_prev; /* remembers previous-to-current node in linked list */
  NR_NODE           *memorized_node_cur;  /* remembers actual node in linked list */

  ret     = 1;  /* default is success */
  fbd     = 0.; /* stores weight of forbidden terms for given q[ij] */
  fbds    = 0.; /* stores weight of forbidden term for given motif */
  qbt1    = 0.;
  q_temp  = 0.;

  n         = vc->length;
  pf_params = vc->exp_params;
  ptype     = vc->ptype;
  S1        = vc->sequence_encoding;
  my_iindx  = vc->iindx;
  jindx     = vc->jindx;

  sc                = vc->sc;
  hc                = vc->hc;
  hc_up_int         = hc->up_int;
  hard_constraints  = hc->mx;

  matrices  = vc->exp_matrices;
  qb        = matrices->qb;
  qm        = matrices->qm;
  qm1       = matrices->qm1;
  scale     = matrices->scale;

  int turn        = pf_params->model_details.min_loop_size;
  int *rtype      = &(pf_params->model_details.rtype[0]);
  n = j;

#ifndef VRNA_NR_SAMPLING_HASH
  if (current_node) {
    memorized_node_prev = NULL;
    memorized_node_cur  = (*current_node)->head;
  }

#endif

  hc_decompose = hard_constraints[jindx[j] + i];

  do {
    int           k, l, kl, u1, u2, max_k, min_l;
    unsigned char type;
    k = i;
    l = j;

    qbr = qb[my_iindx[i] - j];

    if (current_node)
      fbd = NR_TOTAL_WEIGHT(*current_node) * qbr / (*den);

    pstruc[i - 1] = '(';
    pstruc[j - 1] = ')';

    r     = vrna_urn() * (qb[my_iindx[i] - j] - fbd);
    tmp   = qb[my_iindx[i] - j];
    type  = vrna_get_ptype(jindx[j] + i, ptype);
    qbt1  = 0.;

    r             = vrna_urn() * (qb[my_iindx[i] - j] - fbd);
    qbr           = qb[my_iindx[i] - j];
    tmp           = qb[my_iindx[i] - j];
    type          = vrna_get_ptype(jindx[j] + i, ptype);
    hc_decompose  = hard_constraints[n * i + j];

    /* hairpin contribution */
    q_temp = vrna_exp_E_hp_loop(vc, i, j);

    if (current_node) {
        fbds = NR_GET_WEIGHT(*current_node, memorized_node_cur, NRT_HAIRPIN, 0, 0) *
               qbr /
               (*den);
        qbt1 += (q_temp - fbds);
    } else {
      qbt1 += q_temp;
    }

    if (qbt1 >= r) {
      /* found the hairpin we're done */
      if (current_node) {
        *den *= q_temp / qbr;
#ifdef VRNA_NR_SAMPLING_HASH
        *current_node = add_if_nexists(NRT_HAIRPIN, 0, 0, *current_node);
#else
        *current_node = add_if_nexists_ll(memory_dat,
                                          NRT_HAIRPIN,
                                          0,
                                          0,
                                          memorized_node_prev,
                                          memorized_node_cur,
                                          *current_node);
#endif
      }
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
            unsigned int type_2 = rtype[vrna_get_ptype(jindx[l] + k, ptype)];

            /* add *scale[u1+u2+2] */
            q_temp = qb[kl]
                     * scale[u1 + u2 + 2]
                     * exp_E_IntLoop(u1,
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

            if (current_node) {
              fbds = NR_GET_WEIGHT(*current_node, memorized_node_cur, NRT_IT_LOOP, k, l) *
                     qbr /
                     (*den);
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
          *den *= q_temp / qbr;
#ifdef VRNA_NR_SAMPLING_HASH
          *current_node = add_if_nexists(NRT_IT_LOOP, k, l, *current_node);
#else
          *current_node = add_if_nexists_ll(memory_dat,
                                            NRT_IT_LOOP,
                                            k,
                                            l,
                                            memorized_node_prev,
                                            memorized_node_cur,
                                            *current_node);
#endif
        }

        return backtrack(k, l, pstruc, vc, den, current_node, memory_dat); /* found the interior loop, repeat for inside */
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
  if (hard_constraints[jindx[j] + i] & VRNA_CONSTRAINT_CONTEXT_MB_LOOP) {
    int         k, ii, jj, tt;
    FLT_OR_DBL  closingPair;
    tt          = rtype[vrna_get_ptype(jindx[j] + i, ptype)];
    closingPair = pf_params->expMLclosing
                  * exp_E_MLstem(tt, S1[j - 1], S1[i + 1], pf_params)
                  * scale[2];
    if (sc) {
      if (sc->exp_energy_bp)
        q_temp *= sc->exp_energy_bp[jindx[j] + i];

      if (sc->exp_f)
        closingPair *= sc->exp_f(i, j, i, j, VRNA_DECOMP_PAIR_ML, sc->data);
    }

    i++;
    j--;
    /* find the first split index */
    ii  = my_iindx[i];  /* ii-j=[i,j] */
    jj  = jindx[j];     /* jj+i=[j,i] */

    if ((sc) && (sc->exp_f)) {
      for (qt = qbt1, k = i + 1; k < j; k++) {
        q_temp =  qm[ii - (k - 1)] *
                  qm1[jj + k] *
                  closingPair *
                  sc->exp_f(i, j, k - 1, k, VRNA_DECOMP_ML_ML_ML, sc->data);

        if (current_node) {
          fbds = NR_GET_WEIGHT(*current_node, memorized_node_cur, NRT_MT_LOOP, k, 0) *
                 qbr /
                 (*den);
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
        q_temp =  qm[ii - (k - 1)] *
                  qm1[jj + k] *
                  closingPair;

        if (current_node) {
          fbds = NR_GET_WEIGHT(*current_node, memorized_node_cur, NRT_MT_LOOP, k, 0) *
                 qbr /
                 (*den);
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
      if (current_node)
        return 0; /* backtrack failed for non-redundant mode most likely due to numerical instabilities */
      else
        vrna_message_error("backtrack failed, can't find split index ");
    }

    if (current_node) {
      *den *= q_temp / qbr;
#ifdef VRNA_NR_SAMPLING_HASH
      *current_node = add_if_nexists(NRT_MT_LOOP, k, 0, *current_node);
#else
      *current_node = add_if_nexists_ll(memory_dat,
                                        NRT_MT_LOOP,
                                        k,
                                        0,
                                        memorized_node_prev,
                                        memorized_node_cur,
                                        *current_node);
#endif
    }

    ret = backtrack_qm1(k, j, pstruc, vc, den, current_node, memory_dat);

    if (ret == 0)
      return ret;

    j = k - 1;

    ret = (current_node) ?
          backtrack_qm_nr(i, j, pstruc, vc, den, current_node, memory_dat) :
          backtrack_qm(i, j, pstruc, vc);
  }

  return ret;
}


PRIVATE char *
wrap_pbacktrack_circ(vrna_fold_compound_t *vc)
{
  FLT_OR_DBL        r, qt, q_temp;
  int               i, j, k, l, n;
  vrna_exp_param_t  *pf_params;
  FLT_OR_DBL        qo, qmo;
  FLT_OR_DBL        *scale, *qb, *qm, *qm2;
  char              *sequence, *ptype, *pstruc;
  int               *my_iindx, *jindx;
  short             *S1;
  vrna_mx_pf_t      *matrices;
  vrna_sc_t         *sc;

  pf_params = vc->exp_params;
  matrices  = vc->exp_matrices;
  ptype     = vc->ptype;
  my_iindx  = vc->iindx;
  jindx     = vc->jindx;
  S1        = vc->sequence_encoding;
  sc        = vc->sc;

  if ((!matrices) || (!matrices->q) || (!matrices->qb) || (!matrices->qm) || (!pf_params)) {
    vrna_message_warning("vrna_pbacktrack: DP matrices are missing! Call vrna_pf() first!");
    return NULL;
  } else if ((!vc->exp_params->model_details.uniq_ML) || (!matrices->qm1)) {
    vrna_message_warning("vrna_pbacktrack: Unique multiloop decomposition is unset!");
    vrna_message_info(stderr, info_set_uniq_ml);
    return NULL;
  }

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

  /*
   * if (init_length<1)
   *  vrna_message_error("can't backtrack without pf arrays.\n"
   *    "Call pf_circ_fold() before pbacktrack_circ()");
   */

  pstruc = vrna_alloc((n + 1) * sizeof(char));

  /* initialize pstruct with single bases  */
  for (i = 0; i < n; i++)
    pstruc[i] = '.';

  qt  = 1.0 * scale[n];
  /* add soft constraints for open chain configuration */
  if (sc) {
    if (sc->exp_energy_up)
      qt *= sc->exp_energy_up[1][n];

    if (sc->exp_f)
      qt *= sc->exp_f(1, n, 1, n, VRNA_DECOMP_EXT_UP, sc->data);
  }

  r   = vrna_urn() * qo;

  /* open chain? */
  if (qt > r)
    return pstruc;

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

      char loopseq[10];
      if (u < 7) {
        strcpy(loopseq, sequence + j - 1);
        strncat(loopseq, sequence, i);
      }

      qt += qb[my_iindx[i] - j] *
            vrna_exp_E_hp_loop(vc, j, i);

      /* found a hairpin? so backtrack in the enclosed part and we're done  */
      if (qt > r) {
        backtrack(i, j, pstruc, vc, NULL, NULL, NULL);
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

          qt  += q_temp;
          /* found an exterior interior loop? also this time, we can go straight  */
          /* forward and backtracking the both enclosed parts and we're done      */
          if (qt > r) {
            backtrack(i, j, pstruc, vc, NULL, NULL, NULL);
            backtrack(k, l, pstruc, vc, NULL, NULL, NULL);
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
    if ((sc) && (sc->exp_f)) {
      for (k = turn + 2; k < n - 2 * turn - 3; k++) {
        qt += qm[my_iindx[1] - k] *
              qm2[k + 1] *
              expMLclosing *
              sc->exp_f(1, n, k, k + 1, VRNA_DECOMP_ML_ML_ML, sc->data);

        /* backtrack in qm and qm2 if we've found a valid barrier k  */
        if (qt > r) {
          backtrack_qm(1, k, pstruc, vc);
          backtrack_qm2(k + 1, n, pstruc, vc);
          return pstruc;
        }
      }
    } else {
      for (k = turn + 2; k < n - 2 * turn - 3; k++) {
        qt += qm[my_iindx[1] - k] * qm2[k + 1] * expMLclosing;
        /* backtrack in qm and qm2 if we've found a valid barrier k  */
        if (qt > r) {
          backtrack_qm(1, k, pstruc, vc);
          backtrack_qm2(k + 1, n, pstruc, vc);
          return pstruc;
        }
      }
    }
  }
  /* if we reach the actual end of this function, an error has occured  */
  /* cause we HAVE TO find an exterior loop or an open chain!!!         */
  vrna_message_error("backtracking failed in exterior loop");
  return pstruc;
}


PRIVATE char *
pbacktrack_comparative(vrna_fold_compound_t *vc,
                       double               *prob)
{
  FLT_OR_DBL        r, gr, qt;
  int               k, i, j, start, s;
  FLT_OR_DBL        probs   = 1;
  char              *pstruc = NULL;

  int               n_seq       = vc->n_seq;
  int               n           = vc->length;
  short             **S         = vc->S;
  short             **S5        = vc->S5;     /*S5[s][i] holds next base 5' of i in sequence s*/
  short             **S3        = vc->S3;     /*Sl[s][i] holds next base 3' of i in sequence s*/
  unsigned int      **a2s       = vc->a2s;
  vrna_exp_param_t  *pf_params  = vc->exp_params;
  vrna_mx_pf_t      *matrices   = vc->exp_matrices;
  int               *my_iindx   = vc->iindx;

  if ((!matrices) || (!matrices->q) || (!matrices->qb) || (!matrices->qm) || (!pf_params)) {
    vrna_message_warning("vrna_pbacktrack: DP matrices are missing! Call vrna_pf() first!");
    return NULL;
  } else if ((!vc->exp_params->model_details.uniq_ML) || (!matrices->qm1)) {
    vrna_message_warning("vrna_pbacktrack: Unique multiloop decomposition is unset!");
    vrna_message_info(stderr, info_set_uniq_ml);
    return NULL;
  }

  vrna_md_t   *md     = &(pf_params->model_details);
  FLT_OR_DBL  *q      = matrices->q;
  FLT_OR_DBL  *qb     = matrices->qb;
  FLT_OR_DBL  *q1k    = matrices->q1k;
  FLT_OR_DBL  *qln    = matrices->qln;
  FLT_OR_DBL  *scale  = matrices->scale;

  pstruc = vrna_alloc((n + 1) * sizeof(char));

  for (i = 0; i < n; i++)
    pstruc[i] = '.';

  if ((!q1k) || (!qln)) {
    free(q1k);
    free(qln);
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

  start = 1;
  while (start < n) {
    /* find i position of first pair */
    probs = 1.;
    for (i = start; i < n; i++) {
      gr = vrna_urn() * qln[i];
      if (gr > qln[i + 1] * scale[1]) {
        *prob = *prob * probs * (1 - qln[i + 1] * scale[1] / qln[i]);
        break; /* i is paired */
      }

      probs *= qln[i + 1] * scale[1] / qln[i];
    }
    if (i >= n) {
      *prob = *prob * probs;
      break; /* no more pairs */
    }

    /* now find the pairing partner j */
    r = vrna_urn() * (qln[i] - qln[i + 1] * scale[1]);
    for (qt = 0, j = i + 1; j <= n; j++) {
      int         xtype;
      /*  type = ptype[my_iindx[i]-j];
       *  if (type) {*/
      FLT_OR_DBL  qkl;
      if (qb[my_iindx[i] - j] > 0) {
        qkl = qb[my_iindx[i] - j] * qln[j + 1];  /*if psc too small qb=0!*/
        for (s = 0; s < n_seq; s++) {
          xtype = vrna_get_ptype_md(S[s][i], S[s][j], md);
          qkl   *=
            exp_E_ExtLoop(xtype,
                          (a2s[s][i] > 1) ? S5[s][i] : -1,
                          (a2s[s][j] < a2s[s][S[0][0]]) ? S3[s][j] : -1,
                          pf_params);
        }
        qt += qkl;                                                  /*?*exp(pscore[jindx[j]+i]/kTn)*/
        if (qt > r) {
          *prob = *prob * (qkl / (qln[i] - qln[i + 1] * scale[1])); /*probs*=qkl;*/
          break;                                                    /* j is paired */
        }
      }
    }
    if (j == n + 1)
      vrna_message_error("backtracking failed in ext loop");

    start = j + 1;
    backtrack_comparative(vc, pstruc, i, j, prob); /*?*/
  }

  return pstruc;
}


PRIVATE void
backtrack_comparative(vrna_fold_compound_t  *vc,
                      char                  *pstruc,
                      int                   i,
                      int                   j,
                      double                *prob)
{
  int               n_seq       = vc->n_seq;
  short             **S         = vc->S;
  short             **S5        = vc->S5;     /*S5[s][i] holds next base 5' of i in sequence s*/
  short             **S3        = vc->S3;     /*Sl[s][i] holds next base 3' of i in sequence s*/
  char              **Ss        = vc->Ss;
  unsigned int      **a2s       = vc->a2s;
  vrna_exp_param_t  *pf_params  = vc->exp_params;
  vrna_mx_pf_t      *matrices   = vc->exp_matrices;
  vrna_md_t         *md         = &(pf_params->model_details);
  int               *my_iindx   = vc->iindx;
  int               *jindx      = vc->jindx;
  vrna_sc_t         **sc        = vc->scs;
  FLT_OR_DBL        *qb         = matrices->qb;
  FLT_OR_DBL        *qm         = matrices->qm;
  FLT_OR_DBL        *qm1        = matrices->qm1;
  int               *pscore     = vc->pscore;     /* precomputed array of pair types */

  FLT_OR_DBL        *scale      = matrices->scale;
  FLT_OR_DBL        *expMLbase  = matrices->expMLbase;

  /*backtrack given i,j basepair!*/
  FLT_OR_DBL        kTn   = pf_params->kT / 10.;
  int               *type = (int *)vrna_alloc(sizeof(int) * n_seq);

  do {
    FLT_OR_DBL  r, qbt1, max_k, min_l;
    int         k, l, u, u1, u2, s;
    pstruc[i - 1] = '(';
    pstruc[j - 1] = ')';

    for (s = 0; s < n_seq; s++)
      type[s] = vrna_get_ptype_md(S[s][i], S[s][j], md);

    r = vrna_urn() * (qb[my_iindx[i] - j] / exp(pscore[jindx[j] + i] / kTn)); /*?*exp(pscore[jindx[j]+i]/kTn)*/

    qbt1 = 1.;
    for (s = 0; s < n_seq; s++) {
      u = a2s[s][j - 1] - a2s[s][i];
      if (a2s[s][i] < 1)
        continue;

      char loopseq[10] = {
        0
      };
      if (u < 9)
        strncpy(loopseq, Ss[s] + a2s[s][i] - 1, 9);

      qbt1 *= exp_E_Hairpin(u, type[s], S3[s][i], S5[s][j], loopseq, pf_params);
    }
    qbt1 *= scale[j - i + 1];

    if (qbt1 > r) {
      *prob = *prob * qbt1 / (qb[my_iindx[i] - j] / exp(pscore[jindx[j] + i] / kTn)); /*probs*=qbt1;*/
      free(type);
      return; /* found the hairpin we're done */
    }

    max_k = MIN2(i + MAXLOOP + 1, j - TURN - 2);
    l     = MAX2(i + TURN + 2, j - MAXLOOP - 1);
    for (k = i + 1; k <= max_k; k++) {
      min_l = MAX2(k + TURN + 1, j - 1 - MAXLOOP + k - i - 1);

      for (l = min_l; l < j; l++) {
        FLT_OR_DBL  qloop = 1;
        int         type_2;
        if (qb[my_iindx[k] - l] == 0) {
          qloop = 0;
          continue;
        }

        for (s = 0; s < n_seq; s++) {
          u1      = a2s[s][k - 1] - a2s[s][i] /*??*/;
          u2      = a2s[s][j - 1] - a2s[s][l];
          type_2  = vrna_get_ptype_md(S[s][l], S[s][k], md);
          qloop   *= exp_E_IntLoop(u1,
                                   u2,
                                   type[s],
                                   type_2,
                                   S3[s][i],
                                   S5[s][j],
                                   S5[s][k],
                                   S3[s][l],
                                   pf_params);
        }

        if (sc) {
          for (s = 0; s < n_seq; s++) {
            if (sc[s]) {
              int u1  = a2s[s][k - 1] - a2s[s][i];
              int u2  = a2s[s][j - 1] - a2s[s][l];
              if (u1 + u2 == 0) {
                if (sc[s]->exp_energy_stack) {
                  if (S[s][i] && S[s][j] && S[s][k] && S[s][l]) {
                    /* don't allow gaps in stack */
                    qloop *= sc[s]->exp_energy_stack[a2s[s][i]]
                             * sc[s]->exp_energy_stack[a2s[s][k]]
                             * sc[s]->exp_energy_stack[a2s[s][l]]
                             * sc[s]->exp_energy_stack[a2s[s][j]];
                  }
                }
              }
            }
          }
        }

        qbt1 += qb[my_iindx[k] - l] * qloop * scale[k - i + j - l];

        if (qbt1 > r) {
          *prob = *prob
                  * qb[my_iindx[k] - l]
                  * qloop
                  * scale[k - i + j - l]
                  / (qb[my_iindx[i] - j]
                     / exp(pscore[jindx[j] + i] / kTn));
          /*
           * prob*=qb[my_iindx[k]-l] * qloop * scale[k-i+j-l];
           */
          break;
        }
      }
      if (qbt1 > r)
        break;
    }
    if (l < j) {
      i = k;
      j = l;
    } else {
      *prob = *prob * (1 - qbt1 / (qb[my_iindx[i] - j] / exp(pscore[jindx[j] + i] / kTn)));
      break;
    }
  } while (1);

  /* backtrack in multi-loop */
  {
    FLT_OR_DBL  r, qt;
    int         k, ii, jj;
    FLT_OR_DBL  qttemp = 0;
    ;
    i++;
    j--;
    /* find the first split index */
    ii  = my_iindx[i];  /* ii-j=[i,j] */
    jj  = jindx[j];     /* jj+i=[j,i] */
    for (qt = 0., k = i + 1; k < j; k++)
      qttemp += qm[ii - (k - 1)] * qm1[jj + k];
    r = vrna_urn() * qttemp;
    for (qt = 0., k = i + 1; k < j; k++) {
      qt += qm[ii - (k - 1)] * qm1[jj + k];
      if (qt >= r) {
        *prob = *prob
                * qm[ii - (k - 1)]
                * qm1[jj + k]
                / qttemp;/*qttemp;*/
        /*        prob*=qm[ii-(k-1)]*qm1[jj+k];*/
        break;
      }
    }
    if (k >= j)
      vrna_message_error("backtrack failed, can't find split index ");

    backtrack_qm1_comparative(vc, pstruc, k, j, prob);

    j = k - 1;
    while (j > i) {
      /* now backtrack  [i ... j] in qm[] */
      jj  = jindx[j];/*habides??*/
      ii  = my_iindx[i];
      r   = vrna_urn() * qm[ii - j];
      qt  = qm1[jj + i];
      k   = i;
      if (qt < r) {
        for (k = i + 1; k <= j; k++) {
          qt += (qm[ii - (k - 1)] + expMLbase[k - i] /*n_seq??*/) * qm1[jj + k];
          if (qt >= r) {
            *prob = *prob
                    * (qm[ii - (k - 1)] + expMLbase[k - i])
                    * qm1[jj + k]
                    / qm[ii - j];/*???*/
            /*            probs*=qt;*/
            break;
          }
        }
      } else {
        *prob = *prob * qt / qm[ii - j];/*??*/
      }

      if (k > j)
        vrna_message_error("backtrack failed in qm");

      backtrack_qm1_comparative(vc, pstruc, k, j, prob);

      if (k < i + TURN)
        break;             /* no more pairs */

      r = vrna_urn() * (qm[ii - (k - 1)] + expMLbase[k - i]);
      if (expMLbase[k - i] >= r) {
        *prob = *prob * expMLbase[k - i] / (qm[ii - (k - 1)] + expMLbase[k - i]);
        break; /* no more pairs */
      }

      j = k - 1;
      /* whatishere?? */
    }
  }
  free(type);
}


PRIVATE void
backtrack_qm1_comparative(vrna_fold_compound_t  *vc,
                          char                  *pstruc,
                          int                   i,
                          int                   j,
                          double                *prob)
{
  int               n_seq       = vc->n_seq;
  short             **S         = vc->S;
  short             **S5        = vc->S5;     /*S5[s][i] holds next base 5' of i in sequence s*/
  short             **S3        = vc->S3;     /*Sl[s][i] holds next base 3' of i in sequence s*/
  vrna_exp_param_t  *pf_params  = vc->exp_params;
  vrna_mx_pf_t      *matrices   = vc->exp_matrices;
  vrna_md_t         *md         = &(pf_params->model_details);
  int               *my_iindx   = vc->iindx;
  int               *jindx      = vc->jindx;
  FLT_OR_DBL        *qb         = matrices->qb;
  FLT_OR_DBL        *qm1        = matrices->qm1;
  FLT_OR_DBL        *expMLbase  = matrices->expMLbase;

  /* i is paired to l, i<l<j; backtrack in qm1 to find l */
  int               ii, l, xtype, s;
  FLT_OR_DBL        qt, r, tempz;

  r   = vrna_urn() * qm1[jindx[j] + i];
  ii  = my_iindx[i];
  for (qt = 0., l = i + TURN + 1; l <= j; l++) {
    if (qb[ii - l] == 0)
      continue;

    tempz = 1.;
    for (s = 0; s < n_seq; s++) {
      xtype = vrna_get_ptype_md(S[s][i], S[s][l], md);
      tempz *= exp_E_MLstem(xtype, S5[s][i], S3[s][l], pf_params);
    }
    qt += qb[ii - l] * tempz * expMLbase[j - l];
    if (qt >= r) {
      *prob = *prob
              * qb[ii - l]
              * tempz
              * expMLbase[j - l]
              / qm1[jindx[j] + i];
      /* probs*=qb[ii-l]*tempz*expMLbase[j-l];*/
      break;
    }
  }
  if (l > j)
    vrna_message_error("backtrack failed in qm1");

  backtrack_comparative(vc, pstruc, i, l, prob);
}
