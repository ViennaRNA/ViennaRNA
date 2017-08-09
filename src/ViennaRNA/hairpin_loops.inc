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


PRIVATE INLINE
int get_pair_type(int       i,
                  int       j,
                  vrna_md_t *md);


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
