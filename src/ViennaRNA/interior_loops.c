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
E_int_loop(vrna_fold_compound_t *vc,
           int                  i,
           int                  j);


PRIVATE int
E_int_loop_window(vrna_fold_compound_t  *vc,
                  int                   i,
                  int                   j);


PRIVATE int
E_stack_window(vrna_fold_compound_t *vc,
               int                  i,
               int                  j);


PRIVATE int
E_int_loop_comparative(vrna_fold_compound_t *vc,
                       int                  i,
                       int                  j);


PRIVATE int
E_int_loop_comparative_window(vrna_fold_compound_t  *vc,
                              int                   i,
                              int                   j);


PRIVATE INLINE int
eval_int_loop(vrna_fold_compound_t  *vc,
              int                   i,
              int                   j,
              int                   k,
              int                   l);


PRIVATE INLINE int
eval_int_loop_comparative(vrna_fold_compound_t  *vc,
                          int                   i,
                          int                   j,
                          int                   k,
                          int                   l);


/*
 #################################
 # BEGIN OF FUNCTION DEFINITIONS #
 #################################
 */
PUBLIC int
vrna_E_int_loop(vrna_fold_compound_t  *vc,
                int                   i,
                int                   j)
{
  int e = INF;

  if (vc) {
    switch (vc->type) {
      case VRNA_FC_TYPE_SINGLE:
        if (vc->hc->type == VRNA_HC_WINDOW)
          e = E_int_loop_window(vc, i, j);
        else
          e = E_int_loop(vc, i, j);

        break;

      case VRNA_FC_TYPE_COMPARATIVE:
        if (vc->hc->type == VRNA_HC_WINDOW)
          e = E_int_loop_comparative_window(vc, i, j);
        else
          e = E_int_loop_comparative(vc, i, j);

        break;
    }
  }

  return e;
}


PUBLIC int
vrna_eval_int_loop(vrna_fold_compound_t *vc,
                   int                  i,
                   int                  j,
                   int                  k,
                   int                  l)
{
  int e = INF;

  if (vc) {
    switch (vc->type) {
      case VRNA_FC_TYPE_SINGLE:
        e = eval_int_loop(vc, i, j, k, l);
        break;

      case VRNA_FC_TYPE_COMPARATIVE:
        e = eval_int_loop_comparative(vc, i, j, k, l);
        break;
    }
  }

  return e;
}


PRIVATE INLINE int
eval_int_loop(vrna_fold_compound_t  *vc,
              int                   i,
              int                   j,
              int                   k,
              int                   l)
{
  unsigned int  *sn, *ss;
  int           ij, e, *jindx, *rtype, type, type2;
  short         *S, *S2;
  vrna_sc_t     *sc;
  vrna_param_t  *P;
  vrna_md_t     *md;

  jindx = vc->jindx;
  sc    = vc->sc;
  P     = vc->params;
  md    = &(P->model_details);
  sn    = vc->strand_number;
  ss    = vc->strand_start;
  rtype = &(md->rtype[0]);
  ij    = jindx[j] + i;
  S     = vc->sequence_encoding;
  S2    = vc->sequence_encoding2;

  e = INF;

  if ((sn[k] != sn[i]) || (sn[j] != sn[l]))
    return e;

  type  = get_pair_type(S2[i], S2[j], md);
  type2 = get_pair_type(S2[l], S2[k], md);

  e = ubf_eval_int_loop2(i, j, k, l,
                         i + 1, j - 1, k - 1, l + 1,
                         S[i + 1], S[j - 1], S[k - 1], S[l + 1],
                         type, type2, rtype,
                         ij, sn, ss,
                         P, sc);

  return e;
}


PRIVATE INLINE int
eval_int_loop_comparative(vrna_fold_compound_t  *vc,
                          int                   i,
                          int                   j,
                          int                   k,
                          int                   l)
{
  unsigned int  **a2s;
  int           s, ij, e, *jindx, *rtype, type, type_2, n_seq;
  short         **SS, **S5, **S3;
  vrna_sc_t     **scs, *sc;
  vrna_param_t  *P;
  vrna_md_t     *md;

  n_seq = vc->n_seq;
  SS    = vc->S;
  S5    = vc->S5;     /* S5[s][i] holds next base 5' of i in sequence s */
  S3    = vc->S3;     /* Sl[s][i] holds next base 3' of i in sequence s */
  a2s   = vc->a2s;
  jindx = vc->jindx;
  scs   = vc->scs;
  P     = vc->params;
  md    = &(P->model_details);
  rtype = &(md->rtype[0]);
  ij    = jindx[j] + i;

  e = INF;

  for (e = 0, s = 0; s < n_seq; s++) {
    type    = get_pair_type(SS[s][i], SS[s][j], md);
    type_2  = get_pair_type(SS[s][k], SS[s][l], md); /* q,p not p,q! */

    sc = (scs && scs[s]) ? scs[s] : NULL;

    e += ubf_eval_int_loop_comparative(i, j, k, l,
                                       type, type_2, rtype,
                                       ij,
                                       P,
                                       SS[s],
                                       S5[s],
                                       S3[s],
                                       a2s[s],
                                       sc);
  }

  return e;
}


PRIVATE int
E_int_loop(vrna_fold_compound_t *vc,
           int                  i,
           int                  j)
{
  unsigned char             type, type_2;
  char                      *ptype, *ptype_pq;
  unsigned char             *hc_pq, *hc, eval_loop;
  short                     *S, S_i1, S_j1, *S_p1, *S_q1;
  unsigned int              *sn, *so, *ss, *se;
  int                       q, p, j_q, p_i, pq, *c_pq, max_q, max_p, tmp,
                            *rtype, noGUclosure, no_close, energy,
                            *indx, *hc_up, ij, hc_decompose, e, *c, *ggg,
                            with_gquad, turn;
  vrna_sc_t                 *sc;
  vrna_param_t              *P;
  vrna_md_t                 *md;
  vrna_ud_t                 *domains_up;
  vrna_callback_hc_evaluate *evaluate;
  struct default_data       hc_dat_local;

  indx          = vc->jindx;
  hc            = vc->hc->matrix;
  hc_up         = vc->hc->up_int;
  P             = vc->params;
  ij            = indx[j] + i;
  hc_decompose  = hc[ij];
  e             = INF;
  sn            = vc->strand_number;
  so            = vc->strand_order;
  ss            = vc->strand_start;
  se            = vc->strand_end;
  c             = vc->matrices->c;
  ggg           = vc->matrices->ggg;
  md            = &(P->model_details);
  with_gquad    = md->gquad;
  turn          = md->min_loop_size;
  domains_up    = vc->domains_up;
  evaluate      = prepare_hc_default(vc, &hc_dat_local);

  /* CONSTRAINED INTERIOR LOOP start */
  if (hc_decompose & VRNA_CONSTRAINT_CONTEXT_INT_LOOP) {
    /* prepare necessary variables */
    rtype       = &(md->rtype[0]);
    noGUclosure = md->noGUclosure;
    max_q       = i + turn + 2;
    max_q       = MAX2(max_q, j - MAXLOOP - 1);

    ptype     = vc->ptype;
    type      = (unsigned char)ptype[ij];
    no_close  = (((type == 3) || (type == 4)) && noGUclosure);
    S         = vc->sequence_encoding;

    S_i1  = S[i + 1];
    S_j1  = S[j - 1];
    sc    = vc->sc;

    if (type == 0)
      type = 7;

    if (domains_up && domains_up->energy_cb) {
      for (q = j - 1; q >= max_q; q--) {
        j_q = j - q - 1;

        if (hc_up[q + 1] < j_q)
          break;

        pq    = indx[q] + i + 1;
        p_i   = 0;
        max_p = i + 1;
        tmp   = i + 1 + MAXLOOP - j_q;
        max_p = MAX2(max_p, tmp);
        tmp   = q - turn;
        max_p = MIN2(max_p, tmp);
        tmp   = i + 1 + hc_up[i + 1];
        max_p = MIN2(max_p, tmp);
        c_pq  = c + pq;

        ptype_pq  = ptype + pq;
        S_p1      = S + i;
        S_q1      = S + q + 1;

        hc_pq = hc + pq;

        for (p = i + 1; p <= max_p; p++) {
          eval_loop = *hc_pq & VRNA_CONSTRAINT_CONTEXT_INT_LOOP_ENC;
          /* discard this configuration if (p,q) is not allowed to be enclosed pair of an interior loop */
          if (eval_loop && evaluate(i, j, p, q, VRNA_DECOMP_PAIR_IL, &hc_dat_local)) {
            energy = *c_pq;
            if (energy != INF) {
              type_2 = rtype[(unsigned char)*ptype_pq];

              if (type_2 == 0)
                type_2 = 7;

              if (noGUclosure)
                if (no_close || (type_2 == 3) || (type_2 == 4))
                  if ((p > i + 1) || (q < j - 1))
                    continue;

              /* continue unless stack */

              energy  += eval_interior_loop(vc, i, j, p, q);
              e       = MIN2(e, energy);
            }
          }

          hc_pq++;    /* get hc[pq + 1] */
          c_pq++;     /* get c[pq + 1] */
          p_i++;      /* increase unpaired region [i+1...p-1] */

          ptype_pq++; /* get ptype[pq + 1] */
          S_p1++;

          pq++;
        } /* end q-loop */
      }   /* end p-loop */
    } else {
      for (q = j - 1; q >= max_q; q--) {
        j_q = j - q - 1;

        if (hc_up[q + 1] < j_q)
          break;

        pq    = indx[q] + i + 1;
        p_i   = 0;
        max_p = i + 1;
        tmp   = i + 1 + MAXLOOP - j_q;
        max_p = MAX2(max_p, tmp);
        tmp   = q - turn;
        max_p = MIN2(max_p, tmp);
        tmp   = i + 1 + hc_up[i + 1];
        max_p = MIN2(max_p, tmp);
        hc_pq = hc + pq;
        c_pq  = c + pq;

        ptype_pq  = ptype + pq;
        S_p1      = S + i;
        S_q1      = S + q + 1;

        for (p = i + 1; p <= max_p; p++) {
          eval_loop = *hc_pq & VRNA_CONSTRAINT_CONTEXT_INT_LOOP_ENC;
          /* discard this configuration if (p,q) is not allowed to be enclosed pair of an interior loop */
          if (eval_loop && evaluate(i, j, p, q, VRNA_DECOMP_PAIR_IL, &hc_dat_local)) {
            energy = *c_pq;
            if (energy != INF) {
              type_2 = rtype[(unsigned char)*ptype_pq];

              if (noGUclosure)
                if (no_close || (type_2 == 3) || (type_2 == 4))
                  if ((p > i + 1) || (q < j - 1))
                    continue;

              /* continue unless stack */

              if (type_2 == 0)
                type_2 = 7;

              energy += ubf_eval_int_loop2(i, j, p, q,
                                           i + 1, j - 1, p - 1, q + 1,
                                           S_i1, S_j1, *S_p1, *S_q1,
                                           type, type_2, rtype,
                                           ij, sn, ss,
                                           P, sc);
              e = MIN2(e, energy);
            }
          }

          hc_pq++;    /* get hc[pq + 1] */
          c_pq++;     /* get c[pq + 1] */
          p_i++;      /* increase unpaired region [i+1...p-1] */

          ptype_pq++; /* get ptype[pq + 1] */
          S_p1++;

          pq++;
        } /* end q-loop */
      }   /* end p-loop */
    }

    if (with_gquad) {
      /* include all cases where a g-quadruplex may be enclosed by base pair (i,j) */
      if ((!no_close) && (sn[j] == sn[i])) {
        energy  = E_GQuad_IntLoop(i, j, type, S, ggg, indx, P);
        e       = MIN2(e, energy);
      }
    }
  }

  return e;
}


PRIVATE int
E_int_loop_window(vrna_fold_compound_t  *vc,
                  int                   i,
                  int                   j)
{
  char          **ptype;
  unsigned char hc_decompose, eval_loop;
  short         *S1, si1, sj1, *sp, sp1;
  int           e, p, q, minq, turn, noGUclosure, type, type_2, no_close, energy, *rtype, **c,
                **ggg, with_gquad;

  vrna_param_t  *P;
  vrna_md_t     *md;
  vrna_hc_t     *hc;
  vrna_sc_t     *sc;

  e             = INF;
  S1            = vc->sequence_encoding;
  P             = vc->params;
  md            = &(P->model_details);
  c             = vc->matrices->c_local;
  ggg           = vc->matrices->ggg_local;
  hc            = vc->hc;
  sc            = vc->sc;
  ptype         = vc->ptype_local;
  type          = ptype[i][j - i];
  turn          = md->min_loop_size;
  with_gquad    = md->gquad;
  noGUclosure   = md->noGUclosure;
  hc_decompose  = hc->matrix_local[i][j - i];
  rtype         = &(md->rtype[0]);

  if (type == 0)
    type = 7;

  no_close = (((type == 3) || (type == 4)) && noGUclosure);

  if (hc_decompose & VRNA_CONSTRAINT_CONTEXT_INT_LOOP) {
    int tmp, maxp;
    maxp  = j - 2 - turn;
    tmp   = i + MAXLOOP + 1;
    if (maxp > tmp)
      maxp = tmp;

    tmp = i + 1 + hc->up_int[i + 1];
    if (maxp > tmp)
      maxp = tmp;

    si1 = S1[i + 1];
    sj1 = S1[j - 1];
    sp  = S1 + i;
    for (p = i + 1; p <= maxp; p++, sp++) {
      int u1, u2;
      sp1 = S1[p - 1];
      u1  = p - i - 1;

      minq  = j - i + p - MAXLOOP - 2;
      tmp   = p + 1 + turn;
      if (minq < tmp)
        minq = tmp;

      char          *ptype_p = ptype[p];
      ptype_p -= p;
      int           *c_p = c[p];
      c_p -= p;
      unsigned char *hc_p = hc->matrix_local[p];
      hc_p -= p;

      u2 = j - minq - 1;
      /* seek to minimal q value according to hard constraints */
      while (hc->up_int[minq] < u2) {
        u2--;
        minq++;
        if (minq >= j)
          break;
      }

      /* duplicated code is faster than conditions in loop */
      if (hc->f) {
        for (q = minq; q < j; q++) {
          eval_loop = hc_p[q] & VRNA_CONSTRAINT_CONTEXT_INT_LOOP_ENC;
          /* discard this configuration if (p,q) is not allowed to be enclosed pair of an interior loop */
          if (eval_loop && hc->f(i, j, p, q, VRNA_DECOMP_PAIR_IL, hc->data)) {
            energy = c_p[q];

            if (energy < INF) {
              type_2 = ptype_p[q];

              if (type_2 == 0)
                type_2 = 7;

              type_2 = rtype[type_2];

              if ((noGUclosure) && (no_close || (type_2 == 3) || (type_2 == 4)))
                if ((p > i + 1) || (q < j - 1))
                  goto E_int_loop_window_next_q_hc;

              /* continue unless stack */

              energy += E_IntLoop(u1,
                                  u2,
                                  type,
                                  type_2,
                                  si1,
                                  sj1,
                                  sp1,
                                  S1[q + 1],
                                  P);

              /* add soft constraints */
              if (sc) {
                if (sc->energy_up)
                  energy += sc->energy_up[i + 1][u1] +
                            sc->energy_up[q + 1][u2];

                if (sc->energy_bp_local)
                  energy += sc->energy_bp_local[i][j - i];

                if (sc->energy_stack) {
                  if (u1 + u2 == 0) {
                    energy += sc->energy_stack[i] +
                              sc->energy_stack[p] +
                              sc->energy_stack[q] +
                              sc->energy_stack[j];
                  }
                }

                if (sc->f)
                  energy += sc->f(i, j, p, q, VRNA_DECOMP_PAIR_IL, sc->data);
              }

              e = MIN2(e, energy);
            }
          }

E_int_loop_window_next_q_hc:
          u2--;
        } /* end q-loop */
      } else {
        for (q = minq; q < j; q++) {
          /* discard this configuration if (p,q) is not allowed to be enclosed pair of an interior loop */
          if (hc_p[q] & VRNA_CONSTRAINT_CONTEXT_INT_LOOP_ENC) {
            energy = c_p[q];

            if (energy < INF) {
              type_2 = ptype_p[q];

              if (type_2 == 0)
                type_2 = 7;

              type_2 = rtype[type_2];

              if ((noGUclosure) && (no_close || (type_2 == 3) || (type_2 == 4)))
                if ((p > i + 1) || (q < j - 1))
                  goto E_int_loop_window_next_q;

              /* continue unless stack */

              energy += E_IntLoop(u1,
                                  u2,
                                  type,
                                  type_2,
                                  si1,
                                  sj1,
                                  sp1,
                                  S1[q + 1],
                                  P);

              if (sc) {
                if (sc->energy_up)
                  energy += sc->energy_up[i + 1][u1] +
                            sc->energy_up[q + 1][u2];

                if (sc->energy_bp_local)
                  energy += sc->energy_bp_local[i][j - i];

                if (sc->energy_stack) {
                  if (u1 + u2 == 0) {
                    energy += sc->energy_stack[i] +
                              sc->energy_stack[p] +
                              sc->energy_stack[q] +
                              sc->energy_stack[j];
                  }
                }

                if (sc->f)
                  energy += sc->f(i, j, p, q, VRNA_DECOMP_PAIR_IL, sc->data);
              }

              e = MIN2(e, energy);
            }
          }

E_int_loop_window_next_q:
          u2--;
        } /* end q-loop */
      }   /* end if (hc->f) */
    }     /* end p-loop */
  }

  if (with_gquad) {
    /* include all cases where a g-quadruplex may be enclosed by base pair (i,j) */
    if (!no_close) {
      energy  = E_GQuad_IntLoop_L(i, j, type, S1, ggg, vc->window_size, P);
      e       = MIN2(e, energy);
    }
  }

  return e;
}


PRIVATE int
E_int_loop_comparative(vrna_fold_compound_t *vc,
                       int                  i,
                       int                  j)
{
  unsigned char             type, type_2;
  unsigned char             *hc_pq, *hc, eval_loop;
  unsigned int              **a2s;
  short                     **SS, **S5, **S3, *S_cons;
  int                       q, p, j_q, p_i, u, pq, *c_pq, min_q, max_q, max_p, tmp,
                            *rtype, *types, dangle_model, energy, c0, s, n_seq,
                            *indx, *hc_up, ij, hc_decompose, e, *c, *ggg, with_gquad,
                            turn;
  vrna_sc_t                 *sc, **scs;
  vrna_param_t              *P;
  vrna_md_t                 *md;
  vrna_callback_hc_evaluate *evaluate;
  struct default_data       hc_dat_local;

  indx          = vc->jindx;
  hc            = vc->hc->matrix;
  hc_up         = vc->hc->up_int;
  P             = vc->params;
  ij            = indx[j] + i;
  hc_decompose  = hc[ij];
  e             = INF;
  c             = vc->matrices->c;
  ggg           = vc->matrices->ggg;
  md            = &(P->model_details);
  with_gquad    = md->gquad;
  turn          = md->min_loop_size;
  dangle_model  = md->dangles;
  types         = NULL;
  evaluate      = prepare_hc_default(vc, &hc_dat_local);

  /* CONSTRAINED INTERIOR LOOP start */
  if (hc_decompose & VRNA_CONSTRAINT_CONTEXT_INT_LOOP) {
    SS      = vc->S;
    S5      = vc->S5;     /* S5[s][i] holds next base 5' of i in sequence s */
    S3      = vc->S3;     /* Sl[s][i] holds next base 3' of i in sequence s */
    a2s     = vc->a2s;
    S_cons  = vc->S_cons;
    scs     = vc->scs;
    n_seq   = vc->n_seq;
    types   = (int *)vrna_alloc(sizeof(int) * n_seq);

    for (s = 0; s < n_seq; s++)
      types[s] = get_pair_type(SS[s][i], SS[s][j], md);

    /* prepare necessary variables */
    rtype = &(md->rtype[0]);
    max_q = i + turn + 2;
    max_q = MAX2(max_q, j - MAXLOOP - 1);

    for (q = j - 1; q >= max_q; q--) {
      j_q = j - q - 1;

      if (hc_up[q + 1] < j_q)
        break;

      pq    = indx[q] + i + 1;
      p_i   = 0;
      max_p = i + 1;
      tmp   = i + 1 + MAXLOOP - j_q;
      max_p = MAX2(max_p, tmp);
      tmp   = q - turn;
      max_p = MIN2(max_p, tmp);
      tmp   = i + 1 + hc_up[i + 1];
      max_p = MIN2(max_p, tmp);
      hc_pq = hc + pq;
      c_pq  = c + pq;

      for (p = i + 1; p <= max_p; p++) {
        eval_loop = *hc_pq & VRNA_CONSTRAINT_CONTEXT_INT_LOOP_ENC;
        /* discard this configuration if (p,q) is not allowed to be enclosed pair of an interior loop */
        if (eval_loop && evaluate(i, j, p, q, VRNA_DECOMP_PAIR_IL, &hc_dat_local)) {
          energy = *c_pq;
          if (energy != INF) {
            for (s = 0; s < n_seq; s++) {
              type_2 = get_pair_type(SS[s][q], SS[s][p], md); /* q,p not p,q! */

              sc = (scs && scs[s]) ? scs[s] : NULL;

              energy += ubf_eval_int_loop_comparative(i, j, p, q,
                                                      types[s], type_2, rtype,
                                                      ij,
                                                      P,
                                                      SS[s],
                                                      S5[s],
                                                      S3[s],
                                                      a2s[s],
                                                      sc);
            }
            e = MIN2(e, energy);
          }
        }

        hc_pq++;    /* get hc[pq + 1] */
        c_pq++;     /* get c[pq + 1] */
        p_i++;      /* increase unpaired region [i+1...p-1] */

        pq++;
      } /* end q-loop */
    }   /* end p-loop */


    if (with_gquad) {
      /* include all cases where a g-quadruplex may be enclosed by base pair (i,j) */
      energy = 0;
      for (s = 0; s < n_seq; s++) {
        type = types[s];
        if (dangle_model == 2)
          energy += P->mismatchI[type][S3[s][i]][S5[s][j]];

        if (type > 2)
          energy += P->TerminalAU;
      }
      for (p = i + 2; p < j - VRNA_GQUAD_MIN_BOX_SIZE; p++) {
        u = p - i - 1;
        if (u > MAXLOOP)
          break;

        if (S_cons[p] != 3)
          continue;

        min_q = j - i + p - MAXLOOP - 2;
        c0    = p + VRNA_GQUAD_MIN_BOX_SIZE - 1;
        min_q = MAX2(c0, min_q);
        c0    = j - 1;
        max_q = p + VRNA_GQUAD_MAX_BOX_SIZE + 1;
        max_q = MIN2(c0, max_q);
        for (q = min_q; q < max_q; q++) {
          if (S_cons[q] != 3)
            continue;

          c0  = energy + ggg[indx[q] + p] + n_seq * P->internal_loop[u + j - q - 1];
          e   = MIN2(e, c0);
        }
      }

      p = i + 1;
      if (S_cons[p] == 3) {
        if (p < j - VRNA_GQUAD_MIN_BOX_SIZE) {
          min_q = j - i + p - MAXLOOP - 2;
          c0    = p + VRNA_GQUAD_MIN_BOX_SIZE - 1;
          min_q = MAX2(c0, min_q);
          c0    = j - 3;
          max_q = p + VRNA_GQUAD_MAX_BOX_SIZE + 1;
          max_q = MIN2(c0, max_q);
          for (q = min_q; q < max_q; q++) {
            if (S_cons[q] != 3)
              continue;

            c0  = energy + ggg[indx[q] + p] + n_seq * P->internal_loop[j - q - 1];
            e   = MIN2(e, c0);
          }
        }
      }

      q = j - 1;
      if (S_cons[q] == 3) {
        for (p = i + 4; p < j - VRNA_GQUAD_MIN_BOX_SIZE; p++) {
          u = p - i - 1;
          if (u > MAXLOOP)
            break;

          if (S_cons[p] != 3)
            continue;

          c0  = energy + ggg[indx[q] + p] + n_seq * P->internal_loop[u];
          e   = MIN2(e, c0);
        }
      }
    }
  }

  free(types);
  return e;
}


PRIVATE int
E_int_loop_comparative_window(vrna_fold_compound_t  *vc,
                              int                   i,
                              int                   j)
{
  unsigned char type, type_2;
  unsigned char eval_loop;
  unsigned int  **a2s;
  short         **SS, **S5, **S3, *S_cons;
  int           q, p, j_q, p_i, u, min_q, max_q, max_p, tmp,
                *rtype, *types, dangle_model, energy, c0, s, n_seq,
                *hc_up, hc_decompose, e, **c, **ggg, with_gquad,
                turn;
  vrna_sc_t     **scs;
  vrna_param_t  *P;
  vrna_md_t     *md;
  vrna_hc_t     *hc;

  P             = vc->params;
  hc            = vc->hc;
  hc_up         = hc->up_int;
  hc_decompose  = hc->matrix_local[i][j - i];
  e             = INF;
  c             = vc->matrices->c_local;
  ggg           = vc->matrices->ggg_local;
  md            = &(P->model_details);
  with_gquad    = md->gquad;
  turn          = md->min_loop_size;
  dangle_model  = md->dangles;
  types         = NULL;

  /* CONSTRAINED INTERIOR LOOP start */
  if (hc_decompose & VRNA_CONSTRAINT_CONTEXT_INT_LOOP) {
    SS      = vc->S;
    S5      = vc->S5;     /* S5[s][i] holds next base 5' of i in sequence s */
    S3      = vc->S3;     /* Sl[s][i] holds next base 3' of i in sequence s */
    a2s     = vc->a2s;
    S_cons  = vc->S_cons;
    scs     = vc->scs;
    n_seq   = vc->n_seq;
    types   = (int *)vrna_alloc(sizeof(int) * n_seq);

    for (s = 0; s < n_seq; s++)
      types[s] = get_pair_type(SS[s][i], SS[s][j], md);

    /* prepare necessary variables */
    rtype = &(md->rtype[0]);

    int tmp;
    max_p = j - 2 - turn;
    tmp   = i + MAXLOOP + 1;
    if (max_p > tmp)
      max_p = tmp;

    tmp = i + 1 + hc->up_int[i + 1];
    if (max_p > tmp)
      max_p = tmp;

    if (scs) {
      for (p = i + 1; p <= max_p; p++) {
        min_q = j - i + p - MAXLOOP - 2;
        tmp   = p + 1 + turn;
        if (min_q < tmp)
          min_q = tmp;

        int           *c_p = c[p];
        c_p -= p;
        unsigned char *hc_p = hc->matrix_local[p];
        hc_p -= p;

        j_q = j - min_q - 1;
        /* seek to minimal q value according to hard constraints */
        while (hc->up_int[min_q] < j_q) {
          j_q--;
          min_q++;
          if (min_q >= j)
            break;
        }

        if (hc->f) {
          for (q = min_q; q < j; q++) {
            /* discard this configuration if (p,q) is not allowed to be enclosed pair of an interior loop */
            eval_loop = hc_p[q] & VRNA_CONSTRAINT_CONTEXT_INT_LOOP_ENC;
            if (eval_loop && hc->f(i, j, p, q, VRNA_DECOMP_PAIR_IL, hc->data)) {
              energy = c_p[q];

              if (energy < INF) {
                for (s = 0; s < n_seq; s++) {
                  type_2 = get_pair_type(SS[s][q], SS[s][p], md); /* q,p not p,q! */

                  energy += ubf_eval_int_loop_comparative(i, j, p, q,
                                                          types[s], type_2, rtype,
                                                          0,
                                                          P,
                                                          SS[s],
                                                          S5[s],
                                                          S3[s],
                                                          a2s[s],
                                                          scs[s]);
                }
                e = MIN2(e, energy);
              }
            }

            j_q--;
          } /* end q-loop */
        } else {
          for (q = min_q; q < j; q++) {
            /* discard this configuration if (p,q) is not allowed to be enclosed pair of an interior loop */
            if (hc_p[q] & VRNA_CONSTRAINT_CONTEXT_INT_LOOP_ENC) {
              energy = c_p[q];

              if (energy < INF) {
                for (s = 0; s < n_seq; s++) {
                  type_2 = get_pair_type(SS[s][q], SS[s][p], md); /* q,p not p,q! */

                  energy += ubf_eval_int_loop_comparative(i, j, p, q,
                                                          types[s], type_2, rtype,
                                                          0,
                                                          P,
                                                          SS[s],
                                                          S5[s],
                                                          S3[s],
                                                          a2s[s],
                                                          scs[s]);
                }
                e = MIN2(e, energy);
              }
            }

            j_q--;
          } /* end q-loop */
        }
      }     /* end p-loop */
    } else {
      for (p = i + 1; p <= max_p; p++) {
        min_q = j - i + p - MAXLOOP - 2;
        tmp   = p + 1 + turn;
        if (min_q < tmp)
          min_q = tmp;

        int           *c_p = c[p];
        c_p -= p;
        unsigned char *hc_p = hc->matrix_local[p];
        hc_p -= p;

        j_q = j - min_q - 1;
        /* seek to minimal q value according to hard constraints */
        while (hc->up_int[min_q] < j_q) {
          j_q--;
          min_q++;
          if (min_q >= j)
            break;
        }

        if (hc->f) {
          for (q = min_q; q < j; q++) {
            /* discard this configuration if (p,q) is not allowed to be enclosed pair of an interior loop */
            eval_loop = hc_p[q] & VRNA_CONSTRAINT_CONTEXT_INT_LOOP_ENC;
            if (eval_loop && hc->f(i, j, p, q, VRNA_DECOMP_PAIR_IL, hc->data)) {
              energy = c_p[q];

              if (energy < INF) {
                for (s = 0; s < n_seq; s++) {
                  type_2 = get_pair_type(SS[s][q], SS[s][p], md); /* q,p not p,q! */

                  energy += ubf_eval_int_loop_comparative(i, j, p, q,
                                                          types[s], type_2, rtype,
                                                          0,
                                                          P,
                                                          SS[s],
                                                          S5[s],
                                                          S3[s],
                                                          a2s[s],
                                                          NULL);
                }
                e = MIN2(e, energy);
              }
            }

            j_q--;
          } /* end q-loop */
        } else {
          for (q = min_q; q < j; q++) {
            /* discard this configuration if (p,q) is not allowed to be enclosed pair of an interior loop */
            if (hc_p[q] & VRNA_CONSTRAINT_CONTEXT_INT_LOOP_ENC) {
              energy = c_p[q];

              if (energy < INF) {
                for (s = 0; s < n_seq; s++) {
                  type_2 = get_pair_type(SS[s][q], SS[s][p], md); /* q,p not p,q! */

                  energy += ubf_eval_int_loop_comparative(i, j, p, q,
                                                          types[s], type_2, rtype,
                                                          0,
                                                          P,
                                                          SS[s],
                                                          S5[s],
                                                          S3[s],
                                                          a2s[s],
                                                          NULL);
                }
                e = MIN2(e, energy);
              }
            }

            j_q--;
          } /* end q-loop */
        }
      }     /* end p-loop */
    }

    if (with_gquad) {
      /* include all cases where a g-quadruplex may be enclosed by base pair (i,j) */
      energy = 0;
      for (s = 0; s < n_seq; s++) {
        type = types[s];
        if (dangle_model == 2)
          energy += P->mismatchI[type][S3[s][i]][S5[s][j]];

        if (type > 2)
          energy += P->TerminalAU;
      }
      for (p = i + 2; p < j - VRNA_GQUAD_MIN_BOX_SIZE; p++) {
        u = p - i - 1;
        if (u > MAXLOOP)
          break;

        if (S_cons[p] != 3)
          continue;

        min_q = j - i + p - MAXLOOP - 2;
        c0    = p + VRNA_GQUAD_MIN_BOX_SIZE - 1;
        min_q = MAX2(c0, min_q);
        c0    = j - 1;
        max_q = p + VRNA_GQUAD_MAX_BOX_SIZE + 1;
        max_q = MIN2(c0, max_q);
        for (q = min_q; q < max_q; q++) {
          if (S_cons[q] != 3)
            continue;

          c0  = energy + ggg[p][q - p] + n_seq * P->internal_loop[u + j - q - 1];
          e   = MIN2(e, c0);
        }
      }

      p = i + 1;
      if (S_cons[p] == 3) {
        if (p < j - VRNA_GQUAD_MIN_BOX_SIZE) {
          min_q = j - i + p - MAXLOOP - 2;
          c0    = p + VRNA_GQUAD_MIN_BOX_SIZE - 1;
          min_q = MAX2(c0, min_q);
          c0    = j - 3;
          max_q = p + VRNA_GQUAD_MAX_BOX_SIZE + 1;
          max_q = MIN2(c0, max_q);
          for (q = min_q; q < max_q; q++) {
            if (S_cons[q] != 3)
              continue;

            c0  = energy + ggg[p][q - p] + n_seq * P->internal_loop[j - q - 1];
            e   = MIN2(e, c0);
          }
        }
      }

      q = j - 1;
      if (S_cons[q] == 3) {
        for (p = i + 4; p < j - VRNA_GQUAD_MIN_BOX_SIZE; p++) {
          u = p - i - 1;
          if (u > MAXLOOP)
            break;

          if (S_cons[p] != 3)
            continue;

          c0  = energy + ggg[p][q - p] + n_seq * P->internal_loop[u];
          e   = MIN2(e, c0);
        }
      }
    }
  }

  free(types);
  return e;
}


PUBLIC int
vrna_E_ext_int_loop(vrna_fold_compound_t  *vc,
                    int                   i,
                    int                   j,
                    int                   *ip,
                    int                   *iq)
{
  unsigned char             type, type_2;
  int                       ij, q, p, e, s, u1, u2, qmin, energy, *rtype, *types,
                            length, *indx, *hc_up, *c, turn, n_seq;
  char                      *ptype;
  unsigned char             *hc, eval_loop;
  unsigned int              **a2s;
  short                     *S, **SS, **S5, **S3;
  vrna_md_t                 *md;
  vrna_param_t              *P;
  vrna_sc_t                 *sc, **scs;
  vrna_callback_hc_evaluate *evaluate;
  struct default_data       hc_dat_local;

  length    = vc->length;
  indx      = vc->jindx;
  ptype     = vc->ptype;
  c         = vc->matrices->c;
  hc        = vc->hc->matrix;
  hc_up     = vc->hc->up_int;
  P         = vc->params;
  md        = &(P->model_details);
  turn      = md->min_loop_size;
  types     = NULL;
  ij        = indx[j] + i;
  rtype     = &(md->rtype[0]);
  e         = INF;
  evaluate  = prepare_hc_default(vc, &hc_dat_local);

  /* CONSTRAINED INTERIOR LOOP start */
  if (hc[ij] & VRNA_CONSTRAINT_CONTEXT_INT_LOOP) {
    /* prepare necessary variables */
    switch (vc->type) {
      case VRNA_FC_TYPE_SINGLE:
        type = rtype[(unsigned char)ptype[ij]];

        if (type == 0)
          type = 7;

        S   = vc->sequence_encoding;
        sc  = vc->sc;
        break;

      case VRNA_FC_TYPE_COMPARATIVE:
        SS    = vc->S;
        S5    = vc->S5;   /* S5[s][i] holds next base 5' of i in sequence s */
        S3    = vc->S3;   /* Sl[s][i] holds next base 3' of i in sequence s */
        a2s   = vc->a2s;
        scs   = vc->scs;
        n_seq = vc->n_seq;
        types = (int *)vrna_alloc(sizeof(int) * n_seq);

        for (s = 0; s < n_seq; s++)
          types[s] = get_pair_type(SS[s][j], SS[s][i], md);
        break;

      default:
        return e;
        break;
    }

    for (p = j + 1; p < length; p++) {
      u1 = p - j - 1;
      if (u1 + i - 1 > MAXLOOP)
        break;

      if (hc_up[j + 1] < u1)
        break;

      qmin = u1 + i - 1 + length - MAXLOOP;
      if (qmin < p + turn + 1)
        qmin = p + turn + 1;

      for (q = length; q >= qmin; q--) {
        u2 = i - 1 + length - q;
        if (hc_up[q + 1] < u2)
          break;

        if (u1 + u2 > MAXLOOP)
          continue;

        eval_loop = hc[indx[q] + p] & VRNA_CONSTRAINT_CONTEXT_INT_LOOP;

        if (eval_loop && evaluate(i, j, p, q, VRNA_DECOMP_PAIR_IL, &hc_dat_local)) {
          energy = c[indx[q] + p];
          if (energy < INF) {
            switch (vc->type) {
              case VRNA_FC_TYPE_SINGLE:
                type_2 = rtype[(unsigned char)ptype[indx[q] + p]];

                if (type_2 == 0)
                  type_2 = 7;

                energy += ubf_eval_ext_int_loop(i, j, p, q,
                                                i - 1, j + 1, p - 1, q + 1,
                                                S[j + 1], S[i - 1], S[p - 1], S[q + 1],
                                                type, type_2,
                                                length,
                                                P, sc);
                break;

              case VRNA_FC_TYPE_COMPARATIVE:
                for (s = 0; s < n_seq; s++) {
                  type_2 = get_pair_type(SS[s][q], SS[s][p], md); /* q,p not p,q! */

                  sc      = (scs && scs[s]) ? scs[s] : NULL;
                  energy  += ubf_eval_ext_int_loop(a2s[s][i],
                                                   a2s[s][j],
                                                   a2s[s][p],
                                                   a2s[s][q],
                                                   a2s[s][i - 1],
                                                   a2s[s][j + 1],
                                                   a2s[s][p - 1],
                                                   a2s[s][q + 1],
                                                   S3[s][j],
                                                   S5[s][i],
                                                   S5[s][p],
                                                   S3[s][q],
                                                   types[s],
                                                   type_2,
                                                   a2s[s][length],
                                                   P,
                                                   sc);
                }
                break;
            }

            if (energy < e) {
              e = energy;
              if ((ip != NULL) && (iq != NULL)) {
                *ip = p;
                *iq = q;
              }
            }
          }
        }
      }
    }
  }

  free(types);

  return e;
}


PUBLIC int
vrna_E_stack(vrna_fold_compound_t *vc,
             int                  i,
             int                  j)
{
  unsigned char             type, type_2;
  char                      *ptype;
  unsigned char             *hard_constraints, eval_loop;
  unsigned int              **a2s;
  short                     *S, **SS;
  unsigned int              *sn, *ss;
  int                       e, ij, pq, p, q, s, n_seq, *rtype, *indx;
  vrna_sc_t                 *sc, **scs;
  vrna_param_t              *P;
  vrna_md_t                 *md;
  vrna_callback_hc_evaluate *evaluate;
  struct default_data       hc_dat_local;

  if (vc->hc->type == VRNA_HC_WINDOW)
    return E_stack_window(vc, i, j);

  P                 = vc->params;
  md                = &(P->model_details);
  sn                = vc->strand_number;
  ss                = vc->strand_start;
  rtype             = &(md->rtype[0]);
  indx              = vc->jindx;
  hard_constraints  = vc->hc->matrix;
  evaluate          = prepare_hc_default(vc, &hc_dat_local);

  e         = INF;
  p         = i + 1;
  q         = j - 1;
  ij        = indx[j] + i;
  pq        = indx[q] + p;
  eval_loop = (hard_constraints[pq] & VRNA_CONSTRAINT_CONTEXT_INT_LOOP_ENC) &&
              (hard_constraints[ij] & VRNA_CONSTRAINT_CONTEXT_INT_LOOP);

  if ((j - i - 1) < 2)
    return e;

  if (eval_loop && evaluate(i, j, p, q, VRNA_DECOMP_PAIR_IL, &hc_dat_local)) {
    switch (vc->type) {
      case VRNA_FC_TYPE_SINGLE:
        S       = vc->sequence_encoding;
        ptype   = vc->ptype;
        type    = (unsigned char)ptype[ij];
        type_2  = rtype[(unsigned char)ptype[pq]];
        sc      = vc->sc;

        if (type == 0)
          type = 7;

        if (type_2 == 0)
          type_2 = 7;

        if ((sn[p] == sn[i]) && (sn[j] == sn[q])) {
          /* regular stack */
          e = P->stack[type][type_2];
        } else {
          /* stack like cofold structure */
          short si, sj;
          si  = (sn[i + 1] == sn[i]) ? S[i + 1] : -1;
          sj  = (sn[j] == sn[j - 1]) ? S[j - 1] : -1;
          e   = E_IntLoop_Co(rtype[type], rtype[type_2],
                             i, j, p, q,
                             ss[1],
                             si, sj,
                             S[p - 1], S[q + 1],
                             md->dangles,
                             P);
        }

        /* add soft constraints */
        if (sc) {
          if (sc->energy_bp)
            e += sc->energy_bp[ij];

          if (sc->energy_stack) {
            e += sc->energy_stack[i] +
                 sc->energy_stack[p] +
                 sc->energy_stack[q] +
                 sc->energy_stack[j];
          }

          if (sc->f)
            e += sc->f(i, j, p, q, VRNA_DECOMP_PAIR_IL, sc->data);
        }

        break;

      case VRNA_FC_TYPE_COMPARATIVE:
        n_seq = vc->n_seq;
        SS    = vc->S;
        a2s   = vc->a2s;
        scs   = vc->scs;
        e     = 0;
        for (s = 0; s < n_seq; s++) {
          type    = get_pair_type(SS[s][i], SS[s][j], md);
          type_2  = get_pair_type(SS[s][q], SS[s][p], md);  /* q,p not p,q! */
          e       += P->stack[type][type_2];
        }

        if (scs) {
          for (s = 0; s < n_seq; s++) {
            if (scs[s]) {
              if (scs[s]->energy_bp)
                e += scs[s]->energy_bp[ij];

              if (scs[s]->energy_stack) {
                if (SS[s][i] && SS[s][j] && SS[s][p] && SS[s][q]) {
                  /* don't allow gaps in stack */
                  e += scs[s]->energy_stack[a2s[s][i]] +
                       scs[s]->energy_stack[a2s[s][p]] +
                       scs[s]->energy_stack[a2s[s][q]] +
                       scs[s]->energy_stack[a2s[s][j]];
                }
              }

              if (scs[s]->f) {
                e +=
                  scs[s]->f(a2s[s][i],
                            a2s[s][j],
                            a2s[s][p],
                            a2s[s][q],
                            VRNA_DECOMP_PAIR_IL,
                            scs[s]->data);
              }
            }
          }
        }

        break;

      default:
        break;
    }
  }

  return e;
}


PRIVATE int
E_stack_window(vrna_fold_compound_t *vc,
               int                  i,
               int                  j)
{
  char                      **ptype;
  unsigned char             **hard_constraints, eval_loop;
  unsigned int              **a2s;
  short                     **SS;
  int                       e, p, q, *rtype, type, type_2, s, n_seq;
  vrna_param_t              *P;
  vrna_md_t                 *md;
  vrna_sc_t                 *sc, **scs;
  vrna_callback_hc_evaluate *evaluate;
  struct default_data       hc_dat_local;

  P                 = vc->params;
  md                = &(P->model_details);
  rtype             = &(md->rtype[0]);
  hard_constraints  = vc->hc->matrix_local;
  evaluate          = prepare_hc_default(vc, &hc_dat_local);

  e         = INF;
  p         = i + 1;
  q         = j - 1;
  eval_loop = (hard_constraints[p][q - p] & VRNA_CONSTRAINT_CONTEXT_INT_LOOP_ENC) &&
              (hard_constraints[i][j - i] & VRNA_CONSTRAINT_CONTEXT_INT_LOOP);

  if ((j - i - 1) < 2)
    return e;

  if (eval_loop && evaluate(i, j, p, q, VRNA_DECOMP_PAIR_IL, &hc_dat_local)) {
    switch (vc->type) {
      case VRNA_FC_TYPE_SINGLE:
        ptype   = vc->ptype_local;
        type    = ptype[i][j - i];
        type_2  = rtype[ptype[p][q - p]];
        sc      = vc->sc;

        if (type == 0)
          type = 7;

        if (type_2 == 0)
          type_2 = 7;

        /* regular stack */
        e = P->stack[type][type_2];

        /* add soft constraints */
        if (sc) {
          if (sc->energy_stack) {
            e += sc->energy_stack[i] +
                 sc->energy_stack[p] +
                 sc->energy_stack[q] +
                 sc->energy_stack[j];
          }

          if (sc->f)
            e += sc->f(i, j, p, q, VRNA_DECOMP_PAIR_IL, sc->data);
        }

        break;

      case VRNA_FC_TYPE_COMPARATIVE:
        n_seq = vc->n_seq;
        SS    = vc->S;
        a2s   = vc->a2s;
        scs   = vc->scs;
        e     = 0;
        for (s = 0; s < n_seq; s++) {
          type    = get_pair_type(SS[s][i], SS[s][j], md);
          type_2  = get_pair_type(SS[s][q], SS[s][p], md);  /* q,p not p,q! */
          e       += P->stack[type][type_2];
        }

        /* add soft constraints */
        if (scs) {
          for (s = 0; s < n_seq; s++) {
            if (scs[s]) {
              if (scs[s]->energy_stack) {
                if (SS[s][i] && SS[s][j] && SS[s][p] && SS[s][q]) {
                  /* don't allow gaps in stack */
                  e += scs[s]->energy_stack[a2s[s][i]] +
                       scs[s]->energy_stack[a2s[s][p]] +
                       scs[s]->energy_stack[a2s[s][q]] +
                       scs[s]->energy_stack[a2s[s][j]];
                }
              }

              if (scs[s]->f) {
                e +=
                  scs[s]->f(a2s[s][i],
                            a2s[s][j],
                            a2s[s][p],
                            a2s[s][q],
                            VRNA_DECOMP_PAIR_IL,
                            scs[s]->data);
              }
            }
          }
        }

        break;
    }
  }

  return e;
}
