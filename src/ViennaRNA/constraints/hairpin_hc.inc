#include <ViennaRNA/utils/basic.h>


/*
 *  This file contains some utility functions required by both, MFE and
 *  partition function version of hairpin loop evaluation
 */

struct hc_hp_def_dat {
  unsigned int    n;
  unsigned char   *mx;
  unsigned char   **mx_window;
  unsigned int    *sn;
  unsigned int    *hc_up;
  void            *hc_dat;
  vrna_hc_eval_f  hc_f;
};


PRIVATE unsigned char
hc_hp_cb_def(int            i,
             int            j,
             int            k,
             int            l,
             unsigned char  d,
             void           *data);


PRIVATE unsigned char
hc_hp_cb_def_window(int           i,
                    int           j,
                    int           k,
                    int           l,
                    unsigned char d,
                    void          *data);


PRIVATE unsigned char
hc_hp_cb_def_user(int           i,
                  int           j,
                  int           k,
                  int           l,
                  unsigned char d,
                  void          *data);


PRIVATE unsigned char
hc_hp_cb_def_user_window(int            i,
                         int            j,
                         int            k,
                         int            l,
                         unsigned char  d,
                         void           *data);


PRIVATE INLINE vrna_hc_eval_f
prepare_hc_hp_def(vrna_fold_compound_t  *fc,
                  struct hc_hp_def_dat  *dat);


PRIVATE INLINE vrna_hc_eval_f
prepare_hc_hp_def_window(vrna_fold_compound_t *fc,
                         struct hc_hp_def_dat *dat);


/*
 #################################
 # BEGIN OF FUNCTION DEFINITIONS #
 #################################
 */
PRIVATE unsigned char
hc_hp_cb_def(int            i,
             int            j,
             int            k VRNA_UNUSED,
             int            l VRNA_UNUSED,
             unsigned char  d VRNA_UNUSED,
             void           *data)
{
  unsigned int          u, p, q;
  unsigned char         eval;
  struct hc_hp_def_dat  *dat = (struct hc_hp_def_dat *)data;

  eval = (char)0;

  /* no strand nicks are allowed in hairpin loops */
  if (dat->sn[i] != dat->sn[j])
    return eval;

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

  if (dat->mx[dat->n * p + q] & VRNA_CONSTRAINT_CONTEXT_HP_LOOP) {
    eval = (unsigned char)1;
    if (dat->hc_up[i + 1] < u)
      eval = (unsigned char)0;
  }

  return eval;
}


PRIVATE unsigned char
hc_hp_cb_def_window(int           i,
                    int           j,
                    int           k VRNA_UNUSED,
                    int           l VRNA_UNUSED,
                    unsigned char d VRNA_UNUSED,
                    void          *data)
{
  unsigned int          u;
  unsigned char         eval;
  struct hc_hp_def_dat  *dat = (struct hc_hp_def_dat *)data;

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
hc_hp_cb_def_user(int           i,
                  int           j,
                  int           k,
                  int           l,
                  unsigned char d,
                  void          *data)
{
  unsigned char         eval;
  struct hc_hp_def_dat  *dat = (struct hc_hp_def_dat *)data;

  eval  = hc_hp_cb_def(i, j, k, l, d, data);
  eval  = (dat->hc_f(i, j, k, l, d, dat->hc_dat)) ? eval : (unsigned char)0;

  return eval;
}


PRIVATE unsigned char
hc_hp_cb_def_user_window(int            i,
                         int            j,
                         int            k,
                         int            l,
                         unsigned char  d,
                         void           *data)
{
  unsigned char         eval;
  struct hc_hp_def_dat  *dat = (struct hc_hp_def_dat *)data;

  eval  = hc_hp_cb_def_window(i, j, k, l, d, data);
  eval  = (dat->hc_f(i, j, k, l, d, dat->hc_dat)) ? eval : (unsigned char)0;

  return eval;
}


PRIVATE INLINE vrna_hc_eval_f
prepare_hc_hp_def(vrna_fold_compound_t  *fc,
                  struct hc_hp_def_dat  *dat)
{
  dat->mx     = fc->hc->mx;
  dat->hc_up  = fc->hc->up_hp;
  dat->n      = fc->length;
  dat->sn     = fc->strand_number;

  if (fc->hc->f) {
    dat->hc_f   = fc->hc->f;
    dat->hc_dat = fc->hc->data;
    return &hc_hp_cb_def_user;
  }

  return &hc_hp_cb_def;
}


PRIVATE INLINE vrna_hc_eval_f
prepare_hc_hp_def_window(vrna_fold_compound_t *fc,
                         struct hc_hp_def_dat *dat)
{
  dat->mx_window  = fc->hc->matrix_local;
  dat->hc_up      = fc->hc->up_hp;
  dat->n          = fc->length;
  dat->sn         = fc->strand_number;

  if (fc->hc->f) {
    dat->hc_f   = fc->hc->f;
    dat->hc_dat = fc->hc->data;
    return &hc_hp_cb_def_user_window;
  }

  return &hc_hp_cb_def_window;
}
