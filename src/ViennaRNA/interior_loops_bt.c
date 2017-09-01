#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <string.h>
#include "ViennaRNA/fold_vars.h"
#include "ViennaRNA/utils.h"
#include "ViennaRNA/constraints.h"
#include "ViennaRNA/exterior_loops.h"
#include "ViennaRNA/gquad.h"
#include "ViennaRNA/structured_domains.h"
#include "ViennaRNA/unstructured_domains.h"
#include "ViennaRNA/interior_loops.h"


#ifdef __GNUC__
# define INLINE inline
#else
# define INLINE
#endif

#include "interior_loops.inc"

/*
 #################################
 # PRIVATE FUNCTION DECLARATIONS #
 #################################
 */

PRIVATE int
BT_int_loop(vrna_fold_compound_t  *vc,
            int                   *i,
            int                   *j,
            int                   en,
            vrna_bp_stack_t       *bp_stack,
            int                   *stack_count);


PRIVATE int
BT_int_loop_comparative(vrna_fold_compound_t  *vc,
                        int                   *i,
                        int                   *j,
                        int                   en,
                        vrna_bp_stack_t       *bp_stack,
                        int                   *stack_count);


PRIVATE int
BT_int_loop_window(vrna_fold_compound_t *vc,
                   int                  *i,
                   int                  *j,
                   int                  en,
                   vrna_bp_stack_t      *bp_stack,
                   int                  *stack_count);


PRIVATE int
BT_int_loop_window_comparative(vrna_fold_compound_t *vc,
                               int                  *i,
                               int                  *j,
                               int                  en,
                               vrna_bp_stack_t      *bp_stack,
                               int                  *stack_count);


PRIVATE int
BT_stack(vrna_fold_compound_t *vc,
         int                  *i,
         int                  *j,
         int                  *en,
         vrna_bp_stack_t      *bp_stack,
         int                  *stack_count);


PRIVATE int
BT_stack_comparative(vrna_fold_compound_t *vc,
                     int                  *i,
                     int                  *j,
                     int                  *en,
                     vrna_bp_stack_t      *bp_stack,
                     int                  *stack_count);


PRIVATE int
BT_stack_window(vrna_fold_compound_t  *vc,
                int                   *i,
                int                   *j,
                int                   *en,
                vrna_bp_stack_t       *bp_stack,
                int                   *stack_count);


PRIVATE int
BT_stack_window_comparative(vrna_fold_compound_t  *vc,
                            int                   *i,
                            int                   *j,
                            int                   *en,
                            vrna_bp_stack_t       *bp_stack,
                            int                   *stack_count);


/*
 #################################
 # BEGIN OF FUNCTION DEFINITIONS #
 #################################
 */
PUBLIC int
vrna_BT_stack(vrna_fold_compound_t  *vc,
              int                   *i,
              int                   *j,
              int                   *en,
              vrna_bp_stack_t       *bp_stack,
              int                   *stack_count)
{
  if (vc) {
    switch (vc->type) {
      case VRNA_FC_TYPE_SINGLE:
        if (vc->hc->type == VRNA_HC_WINDOW)
          return BT_stack_window(vc, i, j, en, bp_stack, stack_count);
        else
          return BT_stack(vc, i, j, en, bp_stack, stack_count);

        break;

      case VRNA_FC_TYPE_COMPARATIVE:
        if (vc->hc->type == VRNA_HC_WINDOW)
          return BT_stack_window_comparative(vc, i, j, en, bp_stack, stack_count);
        else
          return BT_stack_comparative(vc, i, j, en, bp_stack, stack_count);

        break;
    }
  }

  return 0;
}


PRIVATE int
BT_stack(vrna_fold_compound_t *vc,
         int                  *i,
         int                  *j,
         int                  *en,
         vrna_bp_stack_t      *bp_stack,
         int                  *stack_count)
{
  unsigned char             type, type_2;
  char                      *ptype;
  unsigned char             eval_loop;
  unsigned int              *sn;
  int                       ij, p, q, *idx, *my_c, *rtype, cp;
  vrna_param_t              *P;
  vrna_md_t                 *md;
  vrna_hc_t                 *hc;
  vrna_sc_t                 *sc;
  vrna_callback_hc_evaluate *evaluate;
  struct default_data       hc_dat_local;

  cp    = vc->cutpoint;
  idx   = vc->jindx;
  P     = vc->params;
  md    = &(P->model_details);
  hc    = vc->hc;
  sc    = vc->sc;
  sn    = vc->strand_number;
  my_c  = vc->matrices->c;
  ij    = idx[*j] + *i;
  ptype = vc->ptype;
  type  = (unsigned char)ptype[ij];
  rtype = &(md->rtype[0]);
  p     = *i + 1;
  q     = *j - 1;

  if (vc->hc->f) {
    evaluate            = &hc_default_user;
    hc_dat_local.hc_f   = vc->hc->f;
    hc_dat_local.hc_dat = vc->hc->data;
  } else {
    evaluate = &hc_default;
  }

  if (my_c[ij] == *en) {
    /*  always true, if (i.j) closes canonical structure,
     * thus (i+1.j-1) must be a pair
     */
    eval_loop = (hc->matrix[ij] & VRNA_CONSTRAINT_CONTEXT_INT_LOOP)
                && (hc->matrix[idx[q] + p] & VRNA_CONSTRAINT_CONTEXT_INT_LOOP_ENC);

    if (eval_loop && evaluate(*i, *j, p, q, VRNA_DECOMP_PAIR_IL, &hc_dat_local)) {
      type_2  = ptype[idx[q] + p];
      type_2  = rtype[type_2];

      if (type == 0)
        type = 7;

      if (type_2 == 0)
        type_2 = 7;

      if ((sn[p] == sn[*i]) && (sn[*j] == sn[q])) {
        /* regular stack */
        *en -= P->stack[type][type_2];
      } else {
        /* stack like cofold structure */
        short si, sj, *S;
        S   = vc->sequence_encoding;
        si  = (sn[p] == sn[*i]) ? S[p] : -1;
        sj  = (sn[*j] == sn[q]) ? S[q] : -1;
        *en -= E_IntLoop_Co(rtype[type], rtype[type_2],
                            *i, *j, p, q,
                            cp,
                            si, sj,
                            S[p - 1], S[q + 1],
                            md->dangles,
                            P);
      }

      if (sc) {
        if (sc->energy_bp)
          *en -= sc->energy_bp[ij];

        if (sc->energy_stack) {
          *en -= sc->energy_stack[*i] +
                 sc->energy_stack[p] +
                 sc->energy_stack[q] +
                 sc->energy_stack[*j];
        }

        if (sc->f)
          *en -= sc->f(*i, *j, p, q, VRNA_DECOMP_PAIR_IL, sc->data);
      }

      bp_stack[++(*stack_count)].i  = p;
      bp_stack[(*stack_count)].j    = q;
      (*i)++;
      (*j)--;
      return 1;
    }
  }

  return 0;
}


PRIVATE int
BT_stack_comparative(vrna_fold_compound_t *vc,
                     int                  *i,
                     int                  *j,
                     int                  *en,
                     vrna_bp_stack_t      *bp_stack,
                     int                  *stack_count)
{
  short                     **S;
  int                       type, type_2;
  unsigned char             eval_loop;
  int                       p, q, *c, n_seq, ss, ij, *idx;
  vrna_param_t              *P;
  vrna_md_t                 *md;
  vrna_hc_t                 *hc;
  vrna_sc_t                 **scs;
  vrna_callback_hc_evaluate *evaluate;
  struct default_data       hc_dat_local;

  n_seq = vc->n_seq;
  S     = vc->S;
  P     = vc->params;
  md    = &(P->model_details);
  hc    = vc->hc;
  scs   = vc->scs;
  c     = vc->matrices->c;
  idx   = vc->jindx;
  ij    = idx[*j] + *i;
  p     = *i + 1;
  q     = *j - 1;

  if (vc->hc->f) {
    evaluate            = &hc_default_user;
    hc_dat_local.hc_f   = vc->hc->f;
    hc_dat_local.hc_dat = vc->hc->data;
  } else {
    evaluate = &hc_default;
  }

  if (c[ij] == *en) {
    /*  always true, if (i.j) closes canonical structure,
     * thus (i+1.j-1) must be a pair
     */
    eval_loop = (hc->matrix[ij] & VRNA_CONSTRAINT_CONTEXT_INT_LOOP)
                && (hc->matrix[idx[q] + p] & VRNA_CONSTRAINT_CONTEXT_INT_LOOP_ENC);

    if (eval_loop && evaluate(*i, *j, p, q, VRNA_DECOMP_PAIR_IL, &hc_dat_local)) {
      for (ss = 0; ss < n_seq; ss++) {
        type    = get_pair_type(S[ss][*i], S[ss][*j], md);
        type_2  = get_pair_type(S[ss][q], S[ss][p], md);
        *en     -= P->stack[type][type_2];
      }

      if (scs) {
        for (ss = 0; ss < n_seq; ss++)
          if (scs[ss]) {
            if (scs[ss]->energy_stack) {
              *en -= scs[ss]->energy_stack[*i] +
                     scs[ss]->energy_stack[p] +
                     scs[ss]->energy_stack[q] +
                     scs[ss]->energy_stack[*j];
            }

            if (scs[ss]->f)
              *en -= scs[ss]->f(*i, *j, p, q, VRNA_DECOMP_PAIR_IL, scs[ss]->data);
          }
      }

      *en += vc->pscore[ij];

      bp_stack[++(*stack_count)].i  = p;
      bp_stack[(*stack_count)].j    = q;
      (*i)++;
      (*j)--;
      return 1;
    }
  }

  return 0;
}


PRIVATE int
BT_stack_window(vrna_fold_compound_t  *vc,
                int                   *i,
                int                   *j,
                int                   *en,
                vrna_bp_stack_t       *bp_stack,
                int                   *stack_count)
{
  unsigned char             type, type_2;
  char                      **ptype;
  unsigned char             eval_loop;
  int                       p, q, **c, *rtype;
  vrna_param_t              *P;
  vrna_md_t                 *md;
  vrna_hc_t                 *hc;
  vrna_sc_t                 *sc;
  vrna_callback_hc_evaluate *evaluate;
  struct default_data       hc_dat_local;

  P     = vc->params;
  md    = &(P->model_details);
  hc    = vc->hc;
  sc    = vc->sc;
  c     = vc->matrices->c_local;
  ptype = vc->ptype_local;
  type  = (unsigned char)ptype[*i][*j - *i];
  rtype = &(md->rtype[0]);
  p     = *i + 1;
  q     = *j - 1;

  if (vc->hc->f) {
    evaluate            = &hc_default_user;
    hc_dat_local.hc_f   = vc->hc->f;
    hc_dat_local.hc_dat = vc->hc->data;
  } else {
    evaluate = &hc_default;
  }

  if (c[*i][*j - *i] == *en) {
    /*  always true, if (i.j) closes canonical structure,
     * thus (i+1.j-1) must be a pair
     */
    eval_loop = (hc->matrix_local[*i][*j - *i] & VRNA_CONSTRAINT_CONTEXT_INT_LOOP)
                && (hc->matrix_local[p][q - p] & VRNA_CONSTRAINT_CONTEXT_INT_LOOP_ENC);

    if (eval_loop && evaluate(*i, *j, p, q, VRNA_DECOMP_PAIR_IL, &hc_dat_local)) {
      type_2  = ptype[p][q - p];
      type_2  = rtype[type_2];

      if (type == 0)
        type = 7;

      if (type_2 == 0)
        type_2 = 7;

      *en -= P->stack[type][type_2];

      if (sc) {
        if (sc->energy_stack) {
          *en -= sc->energy_stack[*i] +
                 sc->energy_stack[p] +
                 sc->energy_stack[q] +
                 sc->energy_stack[*j];
        }

        if (sc->f)
          *en -= sc->f(*i, *j, p, q, VRNA_DECOMP_PAIR_IL, sc->data);
      }

      bp_stack[++(*stack_count)].i  = p;
      bp_stack[(*stack_count)].j    = q;
      (*i)++;
      (*j)--;
      return 1;
    }
  }

  return 0;
}


PRIVATE int
BT_stack_window_comparative(vrna_fold_compound_t  *vc,
                            int                   *i,
                            int                   *j,
                            int                   *en,
                            vrna_bp_stack_t       *bp_stack,
                            int                   *stack_count)
{
  int                       type, type_2;
  unsigned char             eval_loop;
  short                     **S;
  int                       p, q, **c, n_seq, ss;
  vrna_param_t              *P;
  vrna_md_t                 *md;
  vrna_hc_t                 *hc;
  vrna_sc_t                 **scs;
  vrna_callback_hc_evaluate *evaluate;
  struct default_data       hc_dat_local;

  n_seq = vc->n_seq;
  S     = vc->S;
  P     = vc->params;
  md    = &(P->model_details);
  hc    = vc->hc;
  scs   = vc->scs;
  c     = vc->matrices->c_local;
  p     = *i + 1;
  q     = *j - 1;

  if (vc->hc->f) {
    evaluate            = &hc_default_user;
    hc_dat_local.hc_f   = vc->hc->f;
    hc_dat_local.hc_dat = vc->hc->data;
  } else {
    evaluate = &hc_default;
  }

  if (c[*i][*j - *i] == *en) {
    /*  always true, if (i.j) closes canonical structure,
     * thus (i+1.j-1) must be a pair
     */
    eval_loop = (hc->matrix_local[*i][*j - *i] & VRNA_CONSTRAINT_CONTEXT_INT_LOOP)
                && (hc->matrix_local[p][q - p] & VRNA_CONSTRAINT_CONTEXT_INT_LOOP_ENC);

    if (eval_loop && evaluate(*i, *j, p, q, VRNA_DECOMP_PAIR_IL, &hc_dat_local)) {
      for (ss = 0; ss < n_seq; ss++) {
        type    = get_pair_type(S[ss][*i], S[ss][*j], md);
        type_2  = get_pair_type(S[ss][q], S[ss][p], md);
        *en     -= P->stack[type][type_2];
      }

      if (scs) {
        for (ss = 0; ss < n_seq; ss++)
          if (scs[ss]) {
            if (scs[ss]->energy_stack) {
              *en -= scs[ss]->energy_stack[*i] +
                     scs[ss]->energy_stack[p] +
                     scs[ss]->energy_stack[q] +
                     scs[ss]->energy_stack[*j];
            }

            if (scs[ss]->f)
              *en -= scs[ss]->f(*i, *j, p, q, VRNA_DECOMP_PAIR_IL, scs[ss]->data);
          }
      }

      *en += vc->pscore_local[*i][*j - *i];

      bp_stack[++(*stack_count)].i  = p;
      bp_stack[(*stack_count)].j    = q;
      (*i)++;
      (*j)--;
      return 1;
    }
  }

  return 0;
}


PUBLIC int
vrna_BT_int_loop(vrna_fold_compound_t *vc,
                 int                  *i,
                 int                  *j,
                 int                  en,
                 vrna_bp_stack_t      *bp_stack,
                 int                  *stack_count)
{
  if (vc) {
    switch (vc->type) {
      case VRNA_FC_TYPE_SINGLE:
        if (vc->hc->type == VRNA_HC_WINDOW)
          return BT_int_loop_window(vc, i, j, en, bp_stack, stack_count);
        else
          return BT_int_loop(vc, i, j, en, bp_stack, stack_count);

        break;

      case VRNA_FC_TYPE_COMPARATIVE:
        if (vc->hc->type == VRNA_HC_WINDOW)
          return BT_int_loop_window_comparative(vc, i, j, en, bp_stack, stack_count);
        else
          return BT_int_loop_comparative(vc, i, j, en, bp_stack, stack_count);

        break;
    }
  }

  return 0;
}


PRIVATE int
BT_int_loop(vrna_fold_compound_t  *vc,
            int                   *i,
            int                   *j,
            int                   en,
            vrna_bp_stack_t       *bp_stack,
            int                   *stack_count)
{
  unsigned char             type, type_2;
  char                      *ptype;
  unsigned char             eval_loop;
  short                     *S1;
  unsigned int              *sn;
  int                       ij, p, q, minq, turn, *idx, noGUclosure, no_close,
                            energy, new, *my_c, *rtype;
  vrna_param_t              *P;
  vrna_md_t                 *md;
  vrna_hc_t                 *hc;
  vrna_sc_t                 *sc;
  vrna_ud_t                 *domains_up;
  vrna_callback_hc_evaluate *evaluate;
  struct default_data       hc_dat_local;

  idx         = vc->jindx;
  P           = vc->params;
  md          = &(P->model_details);
  hc          = vc->hc;
  sc          = vc->sc;
  sn          = vc->strand_number;
  my_c        = vc->matrices->c;
  turn        = md->min_loop_size;
  ij          = idx[*j] + *i;
  ptype       = vc->ptype;
  type        = (unsigned char)ptype[ij];
  rtype       = &(md->rtype[0]);
  S1          = vc->sequence_encoding;
  noGUclosure = md->noGUclosure;
  no_close    = (((type == 3) || (type == 4)) && noGUclosure);
  domains_up  = vc->domains_up;

  if (vc->hc->f) {
    evaluate            = &hc_default_user;
    hc_dat_local.hc_f   = vc->hc->f;
    hc_dat_local.hc_dat = vc->hc->data;
  } else {
    evaluate = &hc_default;
  }

  if (hc->matrix[ij] & VRNA_CONSTRAINT_CONTEXT_INT_LOOP) {
    if (type == 0)
      type = 7;

    if (domains_up && domains_up->energy_cb) {
      for (p = *i + 1; p <= MIN2(*j - 2 - turn, *i + MAXLOOP + 1); p++) {
        minq = *j - *i + p - MAXLOOP - 2;
        if (minq < p + 1 + turn)
          minq = p + 1 + turn;

        if (hc->up_int[*i + 1] < (p - *i - 1))
          break;

        for (q = *j - 1; q >= minq; q--) {
          if (hc->up_int[q + 1] < (*j - q - 1))
            break;

          type_2    = (unsigned char)ptype[idx[q] + p];
          eval_loop = hc->matrix[idx[q] + p] & VRNA_CONSTRAINT_CONTEXT_INT_LOOP_ENC;

          if (!(eval_loop && evaluate(*i, *j, p, q, VRNA_DECOMP_PAIR_IL, &hc_dat_local)))
            continue;

          type_2 = rtype[type_2];

          if (type_2 == 0)
            type_2 = 7;

          if (noGUclosure)
            if (no_close || (type_2 == 3) || (type_2 == 4))
              if ((p > *i + 1) || (q < *j - 1))
                continue;

          /* continue unless stack */

          energy  = eval_interior_loop(vc, *i, *j, p, q);
          new     = energy + my_c[idx[q] + p];

          if (new == en) {
            bp_stack[++(*stack_count)].i  = p;
            bp_stack[(*stack_count)].j    = q;
            if (sc) {
              if (sc->bt) {
                vrna_basepair_t *ptr, *aux_bps;
                aux_bps = sc->bt(*i, *j, p, q, VRNA_DECOMP_PAIR_IL, sc->data);
                for (ptr = aux_bps; ptr && ptr->i != 0; ptr++) {
                  bp_stack[++(*stack_count)].i  = ptr->i;
                  bp_stack[(*stack_count)].j    = ptr->j;
                }
                free(aux_bps);
              }
            }

            *i = p, *j = q;
            return 1; /* success */
          }
        }
      }
    } else {
      for (p = *i + 1; p <= MIN2(*j - 2 - turn, *i + MAXLOOP + 1); p++) {
        minq = *j - *i + p - MAXLOOP - 2;
        if (minq < p + 1 + turn)
          minq = p + 1 + turn;

        if (hc->up_int[*i + 1] < (p - *i - 1))
          break;

        for (q = *j - 1; q >= minq; q--) {
          if (hc->up_int[q + 1] < (*j - q - 1))
            break;

          type_2 = (unsigned char)ptype[idx[q] + p];

          eval_loop = hc->matrix[idx[q] + p] & VRNA_CONSTRAINT_CONTEXT_INT_LOOP_ENC;

          if (!(eval_loop && evaluate(*i, *j, p, q, VRNA_DECOMP_PAIR_IL, &hc_dat_local)))
            continue;

          type_2 = rtype[type_2];

          if (type_2 == 0)
            type_2 = 7;

          if (noGUclosure)
            if (no_close || (type_2 == 3) || (type_2 == 4))
              if ((p > *i + 1) || (q < *j - 1))
                continue;

          /* continue unless stack */

          energy = ubf_eval_int_loop(*i, *j, p, q,
                                     (*i) + 1, (*j) - 1, p - 1, q + 1,
                                     S1[*i + 1], S1[*j - 1], S1[p - 1], S1[q + 1],
                                     type, type_2,
                                     rtype,
                                     ij,
                                     -1,
                                     P,
                                     sc);
          new = energy + my_c[idx[q] + p];

          if (new == en) {
            bp_stack[++(*stack_count)].i  = p;
            bp_stack[(*stack_count)].j    = q;
            if (sc) {
              if (sc->bt) {
                vrna_basepair_t *ptr, *aux_bps;
                aux_bps = sc->bt(*i, *j, p, q, VRNA_DECOMP_PAIR_IL, sc->data);
                for (ptr = aux_bps; ptr && ptr->i != 0; ptr++) {
                  bp_stack[++(*stack_count)].i  = ptr->i;
                  bp_stack[(*stack_count)].j    = ptr->j;
                }
                free(aux_bps);
              }
            }

            *i = p, *j = q;
            return 1; /* success */
          }
        }
      }
    }
  }

  /* is it a g-quadruplex? */
  if (md->gquad) {
    /*
     * The case that is handled here actually resembles something like
     * an interior loop where the enclosing base pair is of regular
     * kind and the enclosed pair is not a canonical one but a g-quadruplex
     * that should then be decomposed further...
     */
    if (sn[*j] == sn[*i]) {
      if (vrna_BT_gquad_int(vc, *i, *j, en, bp_stack, stack_count)) {
        *i = *j = -1; /* tell the calling block to continue backtracking with next block */
        return 1;
      }
    }
  }

  return 0; /* unsuccessful */
}


PRIVATE int
BT_int_loop_comparative(vrna_fold_compound_t  *vc,
                        int                   *i,
                        int                   *j,
                        int                   en,
                        vrna_bp_stack_t       *bp_stack,
                        int                   *stack_count)
{
  unsigned char             eval_loop;
  unsigned int              **a2s;
  short                     **S, **S5, **S3, *S_cons;
  int                       ij, p, q, minq, turn, *idx,
                            energy, new, *my_c, *ggg, n_seq, ss, *type, type_2;
  vrna_param_t              *P;
  vrna_md_t                 *md;
  vrna_hc_t                 *hc;
  vrna_sc_t                 **scs;
  vrna_callback_hc_evaluate *evaluate;
  struct default_data       hc_dat_local;

  n_seq   = vc->n_seq;
  S_cons  = vc->S_cons;
  S       = vc->S;
  S5      = vc->S5;
  S3      = vc->S3;
  a2s     = vc->a2s;
  idx     = vc->jindx;
  P       = vc->params;
  md      = &(P->model_details);
  hc      = vc->hc;
  scs     = vc->scs;
  my_c    = vc->matrices->c;
  ggg     = vc->matrices->ggg;
  turn    = md->min_loop_size;
  ij      = idx[*j] + *i;

  if (vc->hc->f) {
    evaluate            = &hc_default_user;
    hc_dat_local.hc_f   = vc->hc->f;
    hc_dat_local.hc_dat = vc->hc->data;
  } else {
    evaluate = &hc_default;
  }

  if (hc->matrix[ij] & VRNA_CONSTRAINT_CONTEXT_INT_LOOP) {
    type = (int *)vrna_alloc(n_seq * sizeof(int));
    for (ss = 0; ss < n_seq; ss++)
      type[ss] = get_pair_type(S[ss][*i], S[ss][*j], md);

    for (p = *i + 1; p <= MIN2(*j - 2 - turn, *i + MAXLOOP + 1); p++) {
      minq = *j - *i + p - MAXLOOP - 2;
      if (minq < p + 1 + turn)
        minq = p + 1 + turn;

      if (hc->up_int[*i + 1] < (p - *i - 1))
        break;

      for (q = *j - 1; q >= minq; q--) {
        if (hc->up_int[q + 1] < (*j - q - 1))
          break;

        eval_loop = hc->matrix[idx[q] + p] & VRNA_CONSTRAINT_CONTEXT_INT_LOOP_ENC;

        if (!(eval_loop && evaluate(*i, *j, p, q, VRNA_DECOMP_PAIR_IL, &hc_dat_local)))
          continue;

        for (ss = energy = 0; ss < n_seq; ss++) {
          int u1  = a2s[ss][p - 1] - a2s[ss][*i];
          int u2  = a2s[ss][*j - 1] - a2s[ss][q];
          type_2  = get_pair_type(S[ss][q], S[ss][p], md); /* q,p not p,q */
          energy  += E_IntLoop(u1,
                               u2,
                               type[ss],
                               type_2,
                               S3[ss][*i],
                               S5[ss][*j],
                               S5[ss][p],
                               S3[ss][q],
                               P);
        }

        if (scs) {
          for (ss = 0; ss < n_seq; ss++) {
            if (scs[ss]) {
              int u1  = a2s[ss][p - 1] - a2s[ss][*i];
              int u2  = a2s[ss][*j - 1] - a2s[ss][q];
              /*
               *            int u1 = p - i - 1;
               *            int u2 = j - q - 1;
               */
              if (u1 + u2 == 0) {
                if (scs[ss]->energy_stack) {
                  if (S[ss][*i] && S[ss][*j] && S[ss][p] && S[ss][q]) {
                    /* don't allow gaps in stack */
                    energy += scs[ss]->energy_stack[a2s[ss][*i]] +
                              scs[ss]->energy_stack[a2s[ss][p]] +
                              scs[ss]->energy_stack[a2s[ss][q]] +
                              scs[ss]->energy_stack[a2s[ss][*j]];
                  }
                }
              }

              if (scs[ss]->energy_bp)
                energy += scs[ss]->energy_bp[ij];

              if (scs[ss]->energy_up)
                energy += scs[ss]->energy_up[a2s[ss][*i] + 1][u1] +
                          scs[ss]->energy_up[a2s[ss][q] + 1][u2];
            }
          }
        }

        new = energy + my_c[idx[q] + p];

        if (new == en) {
          bp_stack[++(*stack_count)].i  = p;
          bp_stack[(*stack_count)].j    = q;
          if (scs && scs[0]) {
            if (scs[0]->bt) {
              vrna_basepair_t *ptr, *aux_bps;
              aux_bps = scs[0]->bt(*i, *j, p, q, VRNA_DECOMP_PAIR_IL, scs[0]->data);
              for (ptr = aux_bps; ptr && ptr->i != 0; ptr++) {
                bp_stack[++(*stack_count)].i  = ptr->i;
                bp_stack[(*stack_count)].j    = ptr->j;
              }
              free(aux_bps);
            }
          }

          free(type);
          *i = p, *j = q;
          return 1; /* success */
        }
      }
    }

    free(type);
  }

  /* is it a g-quadruplex? */
  if (md->gquad) {
    /*
     * The case that is handled here actually resembles something like
     * an interior loop where the enclosing base pair is of regular
     * kind and the enclosed pair is not a canonical one but a g-quadruplex
     * that should then be decomposed further...
     */
    type = (int *)vrna_alloc(n_seq * sizeof(int));
    for (ss = 0; ss < n_seq; ss++)
      type[ss] = get_pair_type(S[ss][*i], S[ss][*j], md);

    if (backtrack_GQuad_IntLoop_comparative(en, *i, *j, type, S_cons, S5, S3, ggg, idx, &p, &q,
                                            n_seq,
                                            P)) {
      if (vrna_BT_gquad_mfe(vc, p, q, bp_stack, stack_count)) {
        *i = *j = -1; /* tell the calling block to continue backtracking with next block */
        return 1;
      }
    }

    free(type);
  }

  return 0; /* unsuccessful */
}


PRIVATE int
BT_int_loop_window(vrna_fold_compound_t *vc,
                   int                  *i,
                   int                  *j,
                   int                  en,
                   vrna_bp_stack_t      *bp_stack,
                   int                  *stack_count)
{
  int                       type, type_2;
  char                      **ptype;
  unsigned char             eval_loop;
  short                     *S, *S1;
  unsigned int              *sn;
  int                       p, q, minq, turn, maxdist, noGUclosure, no_close,
                            energy, new, **c, **ggg, *rtype, u1, u2;
  vrna_param_t              *P;
  vrna_md_t                 *md;
  vrna_hc_t                 *hc;
  vrna_sc_t                 *sc;
  vrna_callback_hc_evaluate *evaluate;
  struct default_data       hc_dat_local;

  S1          = vc->sequence_encoding;
  S           = vc->sequence_encoding2;
  sn          = vc->strand_number;
  maxdist     = vc->window_size;
  ptype       = vc->ptype_local;
  P           = vc->params;
  md          = &(P->model_details);
  hc          = vc->hc;
  sc          = vc->sc;
  c           = vc->matrices->c_local;
  ggg         = vc->matrices->ggg_local;
  turn        = md->min_loop_size;
  noGUclosure = md->noGUclosure;
  rtype       = &(md->rtype[0]);
  type        = ptype[*i][*j - *i];
  no_close    = (((type == 3) || (type == 4)) && noGUclosure);

  if (hc->f) {
    evaluate            = &hc_default_user;
    hc_dat_local.hc_f   = hc->f;
    hc_dat_local.hc_dat = hc->data;
  } else {
    evaluate = &hc_default;
  }

  if (hc->matrix_local[*i][*j - *i] & VRNA_CONSTRAINT_CONTEXT_INT_LOOP) {
    if (type == 0)
      type = 7;

    for (p = *i + 1; p <= MIN2(*j - 2 - turn, *i + MAXLOOP + 1); p++) {
      u1 = p - (*i) - 1;

      minq = *j - *i + p - MAXLOOP - 2;
      if (minq < p + 1 + turn)
        minq = p + 1 + turn;

      if (hc->up_int[*i + 1] < u1)
        break;

      for (q = *j - 1; q >= minq; q--) {
        u2 = *j - q - 1;

        if (hc->up_int[q + 1] < u2)
          break;

        eval_loop = hc->matrix_local[p][q - p] & VRNA_CONSTRAINT_CONTEXT_INT_LOOP_ENC;

        if (!(eval_loop && evaluate(*i, *j, p, q, VRNA_DECOMP_PAIR_IL, &hc_dat_local)))
          continue;

        type_2  = ptype[p][q - p];
        type_2  = rtype[type_2];

        if (type_2 == 0)
          type_2 = 7;

        if (noGUclosure)
          if (no_close || (type_2 == 3) || (type_2 == 4))
            if ((p > *i + 1) || (q < *j - 1))
              continue;

        /* continue unless stack */

        energy = E_IntLoop(u1,
                           u2,
                           type,
                           type_2,
                           S1[*i + 1],
                           S1[*j - 1],
                           S1[p - 1],
                           S1[q + 1],
                           P);

        if (sc) {
          if (sc->energy_up)
            energy += sc->energy_up[*i + 1][u1] +
                      sc->energy_up[q + 1][u2];

          if (sc->energy_bp_local)
            energy += sc->energy_bp_local[*i][*j - *i];

          if (sc->energy_stack) {
            if (u1 + u2 == 0) {
              energy += sc->energy_stack[*i] +
                        sc->energy_stack[p] +
                        sc->energy_stack[q] +
                        sc->energy_stack[*j];
            }
          }

          if (sc->f)
            energy += sc->f(*i, *j, p, q, VRNA_DECOMP_PAIR_IL, sc->data);
        }

        new = energy + c[p][q - p];

        if (new == en) {
          bp_stack[++(*stack_count)].i  = p;
          bp_stack[(*stack_count)].j    = q;
          if (sc) {
            if (sc->bt) {
              vrna_basepair_t *ptr, *aux_bps;
              aux_bps = sc->bt(*i, *j, p, q, VRNA_DECOMP_PAIR_IL, sc->data);
              for (ptr = aux_bps; ptr && ptr->i != 0; ptr++) {
                bp_stack[++(*stack_count)].i  = ptr->i;
                bp_stack[(*stack_count)].j    = ptr->j;
              }
              free(aux_bps);
            }
          }

          *i = p, *j = q;
          return 1; /* success */
        }
      }
    }
  }

  /* is it a g-quadruplex? */
  if (md->gquad) {
    /*
     * The case that is handled here actually resembles something like
     * an interior loop where the enclosing base pair is of regular
     * kind and the enclosed pair is not a canonical one but a g-quadruplex
     * that should then be decomposed further...
     */
    if (sn[*j] == sn[*i]) {
      if (backtrack_GQuad_IntLoop_L(en, *i, *j, type, S, ggg, maxdist, &p, &q, P)) {
        if (vrna_BT_gquad_mfe(vc, p, q, bp_stack, stack_count)) {
          *i = *j = -1; /* tell the calling block to continue backtracking with next block */
          return 1;
        }
      }
    }
  }

  return 0; /* unsuccessful */
}


PRIVATE int
BT_int_loop_window_comparative(vrna_fold_compound_t *vc,
                               int                  *i,
                               int                  *j,
                               int                  en,
                               vrna_bp_stack_t      *bp_stack,
                               int                  *stack_count)
{
  unsigned char             eval_loop;
  unsigned int              **a2s;
  short                     **S, **SS, **S5, **S3, *S_cons;
  int                       *type, type_2, *rtype;
  int                       p, q, minq, turn, energy, new, **c, **ggg, n_seq, ss;
  vrna_param_t              *P;
  vrna_md_t                 *md;
  vrna_hc_t                 *hc;
  vrna_sc_t                 **scs, *sc;
  vrna_callback_hc_evaluate *evaluate;
  struct default_data       hc_dat_local;

  n_seq   = vc->n_seq;
  S_cons  = vc->S_cons;
  S       = vc->S;
  S5      = vc->S5;       /* S5[s][i] holds next base 5' of i in sequence s */
  S3      = vc->S3;       /* Sl[s][i] holds next base 3' of i in sequence s */
  a2s     = vc->a2s;
  P       = vc->params;
  md      = &(P->model_details);
  hc      = vc->hc;
  scs     = vc->scs;
  c       = vc->matrices->c_local;
  ggg     = vc->matrices->ggg_local;
  turn    = md->min_loop_size;
  rtype   = &(md->rtype[0]);

  if (hc->f) {
    evaluate            = &hc_default_user;
    hc_dat_local.hc_f   = hc->f;
    hc_dat_local.hc_dat = hc->data;
  } else {
    evaluate = &hc_default;
  }

  if (hc->matrix_local[*i][*j - *i] & VRNA_CONSTRAINT_CONTEXT_INT_LOOP) {
    type = (int *)vrna_alloc(n_seq * sizeof(int));
    for (ss = 0; ss < n_seq; ss++)
      type[ss] = get_pair_type(S[ss][*i], S[ss][*j], md);

    for (p = *i + 1; p <= MIN2(*j - 2 - turn, *i + MAXLOOP + 1); p++) {
      minq = *j - *i + p - MAXLOOP - 2;
      if (minq < p + 1 + turn)
        minq = p + 1 + turn;

      if (hc->up_int[*i + 1] < (p - *i - 1))
        break;

      for (q = *j - 1; q >= minq; q--) {
        if (hc->up_int[q + 1] < (*j - q - 1))
          break;

        eval_loop = hc->matrix_local[p][q - p] & VRNA_CONSTRAINT_CONTEXT_INT_LOOP_ENC;

        if (!(eval_loop && evaluate(*i, *j, p, q, VRNA_DECOMP_PAIR_IL, &hc_dat_local)))
          continue;

        for (ss = energy = 0; ss < n_seq; ss++) {
          type_2  = get_pair_type(S[ss][q], S[ss][p], md); /* q,p not p,q! */
          sc      = (scs && scs[ss]) ? scs[ss] : NULL;
          energy  += ubf_eval_int_loop_comparative(*i, *j, p, q,
                                                   type[ss],
                                                   type_2,
                                                   rtype,
                                                   0, -1,
                                                   P,
                                                   S[ss],
                                                   S5[ss],
                                                   S3[ss],
                                                   a2s[ss],
                                                   sc);
        }

        new = energy + c[p][q - p];

        if (new == en) {
          bp_stack[++(*stack_count)].i  = p;
          bp_stack[(*stack_count)].j    = q;
          if (scs && scs[0]) {
            if (scs[0]->bt) {
              vrna_basepair_t *ptr, *aux_bps;
              aux_bps = scs[0]->bt(*i, *j, p, q, VRNA_DECOMP_PAIR_IL, scs[0]->data);
              for (ptr = aux_bps; ptr && ptr->i != 0; ptr++) {
                bp_stack[++(*stack_count)].i  = ptr->i;
                bp_stack[(*stack_count)].j    = ptr->j;
              }
              free(aux_bps);
            }
          }

          free(type);
          *i = p, *j = q;
          return 1; /* success */
        }
      }
    }

    free(type);
  }

  /* is it a g-quadruplex? */
  if (md->gquad) {
    /*
     * The case that is handled here actually resembles something like
     * an interior loop where the enclosing base pair is of regular
     * kind and the enclosed pair is not a canonical one but a g-quadruplex
     * that should then be decomposed further...
     */
    type = (int *)vrna_alloc(n_seq * sizeof(int));
    for (ss = 0; ss < n_seq; ss++)
      type[ss] = get_pair_type(S[ss][*i], S[ss][*j], md);

    if (backtrack_GQuad_IntLoop_L_comparative(en, *i, *j, type, S_cons, S5, S3, ggg, &p, &q, n_seq,
                                              P)) {
      if (vrna_BT_gquad_mfe(vc, p, q, bp_stack, stack_count)) {
        *i = *j = -1; /* tell the calling block to continue backtracking with next block */
        return 1;
      }
    }

    free(type);
  }

  return 0; /* unsuccessful */
}
