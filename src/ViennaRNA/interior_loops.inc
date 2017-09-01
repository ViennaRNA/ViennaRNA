/*
 *  This file contains some utility functions required by both, MFE and
 *  partition function version of interior loop evaluation
 */

struct default_data {
  void                      *hc_dat;
  vrna_callback_hc_evaluate *hc_f;
};


PRIVATE unsigned char
hc_default(int            i,
           int            j,
           int            k,
           int            l,
           unsigned char  d,
           void           *data);


PRIVATE unsigned char
hc_default_user(int           i,
                int           j,
                int           k,
                int           l,
                unsigned char d,
                void          *data);


PRIVATE INLINE int
get_pair_type(int       i,
              int       j,
              vrna_md_t *md);


PRIVATE INLINE int
ubf_eval_int_loop_comparative(int           col_i,
                              int           col_j,
                              int           col_p,
                              int           col_q,
                              unsigned char type,
                              unsigned char type_2,
                              int           *rtype,
                              int           ij,
                              int           cp,
                              vrna_param_t  *P,
                              short         *SS,
                              short         *S5,
                              short         *S3,
                              unsigned int  *a2s,
                              vrna_sc_t     *sc);


/* this is a version of eval_int_loop that also considers unstructured domains */
PRIVATE INLINE int
eval_interior_loop(vrna_fold_compound_t *vc,
                   int                  i,
                   int                  j,
                   int                  p,
                   int                  q);


/*
 #################################
 # BEGIN OF FUNCTION DEFINITIONS #
 #################################
 */
PRIVATE INLINE int
get_pair_type(int       i,
              int       j,
              vrna_md_t *md)
{
  int tt = md->pair[i][j];

  return (tt == 0) ? 7 : tt;
}


PRIVATE unsigned char
hc_default(int            i,
           int            j,
           int            k,
           int            l,
           unsigned char  d,
           void           *data)
{
  return (unsigned char)1;
}


PRIVATE unsigned char
hc_default_user(int           i,
                int           j,
                int           k,
                int           l,
                unsigned char d,
                void          *data)
{
  struct default_data *dat = (struct default_data *)data;

  return dat->hc_f(i, j, k, l, d, dat->hc_dat);
}


PRIVATE INLINE int
ubf_eval_int_loop_comparative(int           col_i,
                              int           col_j,
                              int           col_p,
                              int           col_q,
                              unsigned char type,
                              unsigned char type_2,
                              int           *rtype,
                              int           ij,
                              int           cp,
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

  energy = E_IntLoop(u1, u2, type, type_2, si, sj, sp, sq, P);

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


PRIVATE INLINE int
eval_interior_loop(vrna_fold_compound_t *vc,
                   int                  i,
                   int                  j,
                   int                  p,
                   int                  q)
{
  int           energy, e, e5, e3, u1, u2, i1, j1, p1, q1, *idx;
  short         *S, *S2, si, sj, sp, sq;
  unsigned char type, type_2;
  unsigned int  *sn;
  int           *rtype, ij, cp;
  vrna_param_t  *P;
  vrna_sc_t     *sc;
  vrna_md_t     *md;
  vrna_ud_t     *domains_up;

  S           = vc->sequence_encoding;
  S2          = vc->sequence_encoding2;
  cp          = vc->cutpoint;
  P           = vc->params;
  sc          = vc->sc;
  domains_up  = vc->domains_up;
  md          = &(P->model_details);
  sn          = vc->strand_number;
  rtype       = &(md->rtype[0]);
  idx         = vc->jindx;
  type        = get_pair_type(S2[i], S2[j], md);
  type_2      = get_pair_type(S2[q], S2[p], md);

  i1  = i + 1;
  j1  = j - 1;
  p1  = p - 1;
  q1  = q + 1;
  u1  = p1 - i;
  u2  = j1 - q;

  si  = S[i1];
  sj  = S[j1];
  sp  = S[p1];
  sq  = S[q1];
  ij  = idx[j] + i;

  if ((sn[p] == sn[i]) && (sn[j] == sn[q])) {
    /* regular interior loop */
    energy = E_IntLoop(u1, u2, type, type_2, si, sj, sp, sq, P);
  } else {
    /* interior loop like cofold structure */
    short Si, Sj;
    Si      = (sn[i1] == sn[i]) ? si : -1;
    Sj      = (sn[j] == sn[j1]) ? sj : -1;
    energy  = E_IntLoop_Co(rtype[type], rtype[type_2],
                           i, j, p, q,
                           cp,
                           Si, Sj,
                           sp, sq,
                           P->model_details.dangles,
                           P);
  }

  /* add soft constraints */
  if (sc) {
    if (sc->energy_up)
      energy += sc->energy_up[i1][u1] +
                sc->energy_up[q1][u2];

    if (sc->energy_bp)
      energy += sc->energy_bp[ij];

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

  e   = energy;
  e5  = e3 = 0;

  if (u1 > 0) {
    e5 = domains_up->energy_cb(vc,
                               i + 1, p - 1,
                               VRNA_UNSTRUCTURED_DOMAIN_INT_LOOP,
                               domains_up->data);
  }

  if (u2 > 0) {
    e3 = domains_up->energy_cb(vc,
                               q + 1, j - 1,
                               VRNA_UNSTRUCTURED_DOMAIN_INT_LOOP,
                               domains_up->data);
  }

  e = MIN2(e, energy + e5);
  e = MIN2(e, energy + e3);
  e = MIN2(e, energy + e5 + e3);


  return e;
}
