#include <ViennaRNA/utils/basic.h>


/*
 *  This file contains some utility functions required by both, MFE and
 *  partition function version of internal loop evaluation
 */

typedef unsigned char (*eval_hc)(int   i,
                                int   j,
                                int   k,
                                int   l,
                                void  *data);

struct hc_int_def_dat {
  unsigned char             *mx;
  unsigned char             **mx_local;
  unsigned int              *sn;
  unsigned int              n;
  unsigned int              *up;

  void                      *hc_dat;
  vrna_hc_eval_f hc_f;
};


PRIVATE unsigned char
hc_int_cb_def(int   i,
              int   j,
              int   k,
              int   l,
              void  *data);


PRIVATE unsigned char
hc_int_cb_def_user(int  i,
                   int  j,
                   int  k,
                   int  l,
                   void *data);


PRIVATE INLINE int
ubf_eval_int_loop_comparative(int           col_i,
                              int           col_j,
                              int           col_p,
                              int           col_q,
                              unsigned char type,
                              unsigned char type_2,
                              int           *rtype,
                              int           ij,
                              vrna_param_t  *P,
                              short         *SS,
                              short         *S5,
                              short         *S3,
                              unsigned int  *a2s,
                              vrna_sc_t     *sc);


PRIVATE INLINE eval_hc
prepare_hc_int_def(vrna_fold_compound_t   *fc,
                   struct hc_int_def_dat  *dat);


/*
 #################################
 # BEGIN OF FUNCTION DEFINITIONS #
 #################################
 */
PRIVATE unsigned char
hc_int_cb_def(int   i,
              int   j,
              int   k,
              int   l,
              void  *data)
{
  unsigned char pij, pkl;
  unsigned int  u1, u2, u3;

  struct hc_int_def_dat *dat = (struct hc_int_def_dat *)data;

  if ((dat->sn[i] != dat->sn[k]) ||
      (dat->sn[l] != dat->sn[j]))
    return (unsigned char)0;

  if (j < k) {
    /* circular RNA */
    u1 = i - 1;
    u2 = k - j;
    u3 = dat->n - l;

    if ((u1) &&
        (dat->up[1] < u1))
      return (unsigned char)0;

    if ((u2) &&
        (dat->up[j + 1] < u2))
      return (unsigned char)0;

    if ((u3) &&
       (dat->up[l + 1] < u3))
      return (unsigned char)0;

  } else {
    /* linear RNA */
    u1 = k - i - 1;
    u2 = j - l - 1;

    if ((u1) &&
        (dat->up[i + 1] < u1))
      return (unsigned char)0;

    if ((u2) &&
        (dat->up[l + 1] < u2))
      return (unsigned char)0;
  }

  if (dat->mx) {
    pij = dat->mx[dat->n * i + j];
    pkl = dat->mx[dat->n * k + l];
  } else {
    pij = dat->mx_local[i][j - i];
    pkl = dat->mx_local[k][l - k];
  }

  if ((pij & VRNA_CONSTRAINT_CONTEXT_INT_LOOP) &&
      (pkl & VRNA_CONSTRAINT_CONTEXT_INT_LOOP_ENC))
    return (unsigned char)1;

  return (unsigned char)0;
}


PRIVATE unsigned char
hc_int_cb_def_user(int  i,
                   int  j,
                   int  k,
                   int  l,
                   void *data)
{
  struct hc_int_def_dat *dat = (struct hc_int_def_dat *)data;

  unsigned char eval = hc_int_cb_def(i, j, k, l, data);

  return (dat->hc_f(i, j, k, l, VRNA_DECOMP_PAIR_IL, dat->hc_dat)) ? eval : (unsigned char)0;
}


PRIVATE INLINE eval_hc
prepare_hc_int_def(vrna_fold_compound_t   *fc,
                   struct hc_int_def_dat  *dat)
{
  dat->mx       = (fc->hc->type == VRNA_HC_WINDOW) ? NULL : fc->hc->mx;
  dat->mx_local = (fc->hc->type == VRNA_HC_WINDOW) ? fc->hc->matrix_local : NULL;
  dat->n        = fc->length;
  dat->up       = fc->hc->up_int;
  dat->sn       = fc->strand_number;
  dat->hc_f     = NULL;
  dat->hc_dat   = NULL;

  if (fc->hc->f) {
    dat->hc_f   = fc->hc->f;
    dat->hc_dat = fc->hc->data;
    return &hc_int_cb_def_user;
  }

  return &hc_int_cb_def;
}


PRIVATE INLINE int
ubf_eval_int_loop_comparative(int           col_i,
                              int           col_j,
                              int           col_p,
                              int           col_q,
                              unsigned char type,
                              unsigned char type_2,
                              int           *rtype VRNA_UNUSED,
                              int           ij,
                              vrna_param_t  *P,
                              short         *SS,
                              short         *S5,
                              short         *S3,
                              unsigned int  *a2s,
                              vrna_sc_t     *sc)
{
  short si, sj, sp, sq;
  int   energy, u1, u2;
  int   i, j, p, q, i1, j1, p1, q1;

  i   = a2s[col_i];
  j   = a2s[col_j];
  p   = a2s[col_p];
  q   = a2s[col_q];
  i1  = a2s[col_i + 1];
  j1  = a2s[col_j - 1];
  p1  = a2s[col_p - 1];
  q1  = a2s[col_q + 1];

  si  = S3[col_i];
  sj  = S5[col_j];
  sp  = S5[col_p];
  sq  = S3[col_q];

  u1  = p1 - i;
  u2  = j1 - q;

  energy = vrna_E_internal(u1, u2, type, type_2, si, sj, sp, sq, P);

  /* add soft constraints */
  if (sc) {
    if (sc->energy_up)
      energy += sc->energy_up[i1][u1]
                + sc->energy_up[q1][u2];

    if (sc->energy_bp)
      energy += sc->energy_bp[ij];

    if (sc->energy_stack) {
      if (u1 + u2 == 0) {
        if (SS[col_i] && SS[col_j] && SS[col_p] && SS[col_q]) {
          /* no gap allowed */
          energy += sc->energy_stack[i] +
                    sc->energy_stack[p] +
                    sc->energy_stack[q] +
                    sc->energy_stack[j];
        }
      }
    }

    if (sc->f)
      energy += sc->f(i, j, p, q, VRNA_DECOMP_PAIR_IL, sc->data);
  }

  return energy;
}
