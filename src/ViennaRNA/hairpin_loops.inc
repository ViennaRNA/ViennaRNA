/*
 *  This file contains some utility functions required by both, MFE and
 *  partition function version of hairpin loop evaluation
 */

struct default_data {
  int                       n;
  int                       *idx;
  unsigned char             *mx;
  unsigned char             **mx_window;
  unsigned int              *sn;
  int                       *hc_up;
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
hc_default_window(int           i,
                  int           j,
                  int           k,
                  int           l,
                  unsigned char d,
                  void          *data);


PRIVATE unsigned char
hc_default_user(int           i,
                int           j,
                int           k,
                int           l,
                unsigned char d,
                void          *data);


PRIVATE unsigned char
hc_default_user_window(int            i,
                       int            j,
                       int            k,
                       int            l,
                       unsigned char  d,
                       void           *data);


PRIVATE INLINE vrna_callback_hc_evaluate *
prepare_hc_default(vrna_fold_compound_t *vc,
                   struct default_data  *dat);


PRIVATE INLINE vrna_callback_hc_evaluate *
prepare_hc_default_window(vrna_fold_compound_t  *vc,
                          struct default_data   *dat);


/*
 #################################
 # BEGIN OF FUNCTION DEFINITIONS #
 #################################
 */
PRIVATE unsigned char
hc_default(int            i,
           int            j,
           int            k,
           int            l,
           unsigned char  d,
           void           *data)
{
  int                 ij, u, p, q;
  unsigned char       eval;
  struct default_data *dat = (struct default_data *)data;

  eval = (char)0;

  if (j > i) {
    /* linear case */
    p = i;
    q = j;
    u = q - p - 1;
  } else {
    /* circular case */
    p = j;
    q = i;
    u = dat->n - q + p - 1;
  }

  ij = dat->idx[q] + p;
  if (dat->mx[ij] & VRNA_CONSTRAINT_CONTEXT_HP_LOOP) {
    eval = (unsigned char)1;
    if (dat->hc_up[i + 1] < u)
      eval = (unsigned char)0;
  }

  return eval;
}


PRIVATE unsigned char
hc_default_window(int           i,
                  int           j,
                  int           k,
                  int           l,
                  unsigned char d,
                  void          *data)
{
  int                 u;
  unsigned char       eval;
  struct default_data *dat = (struct default_data *)data;

  eval = (unsigned char)0;

  u = j - i - 1;

  if (dat->mx_window[i][j - i] & VRNA_CONSTRAINT_CONTEXT_HP_LOOP) {
    eval = (unsigned char)1;
    if (dat->hc_up[i + 1] < u)
      eval = (unsigned char)0;
  }

  return eval;
}


PRIVATE unsigned char
hc_default_user(int           i,
                int           j,
                int           k,
                int           l,
                unsigned char d,
                void          *data)
{
  unsigned char       eval;
  struct default_data *dat = (struct default_data *)data;

  eval  = hc_default(i, j, k, l, d, data);
  eval  = (dat->hc_f(i, j, k, l, d, dat->hc_dat)) ? eval : (unsigned char)0;

  return eval;
}


PRIVATE unsigned char
hc_default_user_window(int            i,
                       int            j,
                       int            k,
                       int            l,
                       unsigned char  d,
                       void           *data)
{
  unsigned char       eval;
  struct default_data *dat = (struct default_data *)data;

  eval  = hc_default_window(i, j, k, l, d, data);
  eval  = (dat->hc_f(i, j, k, l, d, dat->hc_dat)) ? eval : (unsigned char)0;

  return eval;
}


PRIVATE INLINE vrna_callback_hc_evaluate *
prepare_hc_default(vrna_fold_compound_t *vc,
                   struct default_data  *dat)
{
  dat->idx    = vc->jindx;
  dat->mx     = vc->hc->matrix;
  dat->hc_up  = vc->hc->up_hp;
  dat->n      = vc->length;
  dat->sn     = vc->strand_number;

  if (vc->hc->f) {
    dat->hc_f   = vc->hc->f;
    dat->hc_dat = vc->hc->data;
    return &hc_default_user;
  }

  return &hc_default;
}


PRIVATE INLINE vrna_callback_hc_evaluate *
prepare_hc_default_window(vrna_fold_compound_t  *vc,
                          struct default_data   *dat)
{
  dat->idx        = vc->jindx;
  dat->mx_window  = vc->hc->matrix_local;
  dat->hc_up      = vc->hc->up_hp;
  dat->n          = vc->length;
  dat->sn         = vc->strand_number;

  if (vc->hc->f) {
    dat->hc_f   = vc->hc->f;
    dat->hc_dat = vc->hc->data;
    return &hc_default_user_window;
  }

  return &hc_default_window;
}
