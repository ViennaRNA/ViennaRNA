/*
 *  This file contains hard constraint wrappers required by both, MFE and
 *  partition function version of exterior loop evaluation
 */

#include "ViennaRNA/utils/log.h"

struct hc_ext_def_dat {
  unsigned int              n;
  unsigned char             *mx;
  unsigned char             **mx_window;
  unsigned int              *sn;
  unsigned int              *hc_up;
  void                      *hc_dat;
  vrna_hc_eval_f hc_f;
};

PRIVATE unsigned char
hc_ext_cb_def(int           i,
              int           j,
              int           k,
              int           l,
              unsigned char d,
              void          *data);


PRIVATE unsigned char
hc_ext_cb_sn(int            i,
             int            j,
             int            k,
             int            l,
             unsigned char  d,
             void           *data);


PRIVATE unsigned char
hc_ext_cb_def_window(int            i,
                     int            j,
                     int            k,
                     int            l,
                     unsigned char  d,
                     void           *data);


PRIVATE unsigned char
hc_ext_cb_def_user(int            i,
                   int            j,
                   int            k,
                   int            l,
                   unsigned char  d,
                   void           *data);


PRIVATE unsigned char
hc_ext_cb_def_sn_user(int           i,
                      int           j,
                      int           k,
                      int           l,
                      unsigned char d,
                      void          *data);


PRIVATE unsigned char
hc_ext_cb_def_user_window(int           i,
                          int           j,
                          int           k,
                          int           l,
                          unsigned char d,
                          void          *data);


PRIVATE INLINE vrna_hc_eval_f
prepare_hc_ext_def(vrna_fold_compound_t   *fc,
                   struct hc_ext_def_dat  *dat);


PRIVATE INLINE vrna_hc_eval_f
prepare_hc_ext_def_window(vrna_fold_compound_t  *fc,
                          struct hc_ext_def_dat *dat);


/*
 #################################
 # BEGIN OF FUNCTION DEFINITIONS #
 #################################
 */
PRIVATE unsigned char
hc_ext_cb_def(int           i,
              int           j,
              int           k,
              int           l,
              unsigned char d,
              void          *data)
{
  unsigned int          n, di, dj;
  unsigned char         eval;
  struct hc_ext_def_dat *dat = (struct hc_ext_def_dat *)data;

  eval  = (unsigned char)0;
  di    = k - i;
  dj    = j - l;
  n     = dat->n;

  switch (d) {
    case VRNA_DECOMP_EXT_EXT_STEM:
      if (dat->mx[n * j + l] & VRNA_CONSTRAINT_CONTEXT_EXT_LOOP) {
        eval = (unsigned char)1;
        if (i != l) {
          /* otherwise, stem spans from i to j */
          di = l - k - 1;
          if ((di != 0) && (dat->hc_up[k + 1] < di))
            eval = (unsigned char)0;
        }
      }

      break;

    case VRNA_DECOMP_EXT_STEM_EXT:
      if (dat->mx[n * k + i] & VRNA_CONSTRAINT_CONTEXT_EXT_LOOP) {
        eval = (unsigned char)1;
        if (i != l) {
          /* otherwise, stem spans from i to j */
          di = l - k - 1;
          if ((di != 0) && (dat->hc_up[k + 1] < di))
            eval = (unsigned char)0;
        }
      }

      break;

    case VRNA_DECOMP_EXT_EXT_STEM1:
      if (dat->mx[n * (j - 1) + l] & VRNA_CONSTRAINT_CONTEXT_EXT_LOOP) {
        eval = (unsigned char)1;
        if (dat->hc_up[j] == 0)
          eval = (unsigned char)0;

        if (i != l) {
          /* otherwise, stem spans from i to j - 1 */
          di = l - k - 1;

          if ((di != 0) && (dat->hc_up[k + 1] < di))
            eval = (unsigned char)0;
        }
      }

      break;

    case VRNA_DECOMP_EXT_STEM_EXT1:
      if (dat->mx[n * k + i + 1] & VRNA_CONSTRAINT_CONTEXT_EXT_LOOP) {
        eval = (unsigned char)1;

        if (dat->hc_up[i] == 0)
          eval = (unsigned char)0;

        if (j != k) {
          /* otherwise, stem spans from i + 1 to j */
          dj = l - k - 1;

          if ((dj != 0) && (dat->hc_up[k + 1] < dj))
            eval = (unsigned char)0;
        }
      }

      break;

    case VRNA_DECOMP_EXT_EXT_EXT:
      eval  = (unsigned char)1;
      di    = l - k - 1;
      if ((di != 0) && (dat->hc_up[k + 1] < di))
        eval = (unsigned char)0;

      break;

    case VRNA_DECOMP_EXT_STEM:
      if (dat->mx[n * k + l] & VRNA_CONSTRAINT_CONTEXT_EXT_LOOP) {
        eval = (unsigned char)1;
        if ((di != 0) && (dat->hc_up[i] < di))
          eval = (unsigned char)0;

        if ((dj != 0) && (dat->hc_up[l + 1] < dj))
          eval = (unsigned char)0;
      }

      break;

    case VRNA_DECOMP_EXT_EXT:
      eval = (unsigned char)1;
      if ((di != 0) && (dat->hc_up[i] < di))
        eval = (unsigned char)0;

      if ((dj != 0) && (dat->hc_up[l + 1] < dj))
        eval = (unsigned char)0;

      break;

    case VRNA_DECOMP_EXT_UP:
      di    = j - i + 1;
      eval  = (dat->hc_up[i] >= di) ? (unsigned char)1 : (unsigned char)0;
      break;

    case VRNA_DECOMP_EXT_STEM_OUTSIDE:
      if (dat->mx[n * k + l] & VRNA_CONSTRAINT_CONTEXT_EXT_LOOP)
        eval = (unsigned char)1;

      break;

    default:
      vrna_log_warning("hc_cb@exterior_loops.c: "
                           "Unrecognized decomposition %d",
                           d);
  }

  return eval;
}


PRIVATE unsigned char
hc_ext_cb_sn(int            i,
             int            j,
             int            k,
             int            l,
             unsigned char  d,
             void           *data)
{
  unsigned int          *sn;
  unsigned char         eval;
  struct hc_ext_def_dat *dat = (struct hc_ext_def_dat *)data;

  sn    = dat->sn;
  eval  = (unsigned char)0;

  switch (d) {
    case VRNA_DECOMP_EXT_EXT_STEM1:
      if (sn[j - 1] != sn[j])
        break;

      if (sn[k] == sn[l])
        eval = (unsigned char)1;

      break;

    case VRNA_DECOMP_EXT_STEM_EXT1:
      if (sn[i] != sn[i + 1])
        break;

      if (sn[k] == sn[l])
        eval = (unsigned char)1;

      break;

    case VRNA_DECOMP_EXT_EXT_STEM:
    /* fall through */
    case VRNA_DECOMP_EXT_STEM_EXT:
    /* fall through */
    case VRNA_DECOMP_EXT_EXT_EXT:
      if (sn[k] == sn[l])
        eval = (unsigned char)1;

      break;

    case VRNA_DECOMP_EXT_STEM:
    /* fall through */
    case VRNA_DECOMP_EXT_EXT:
      if ((sn[i] == sn[k]) && (sn[l] == sn[j]))
        eval = (unsigned char)1;

      break;

    case VRNA_DECOMP_EXT_UP:
      if (sn[i] == sn[j])
        eval = (unsigned char)1;

      break;

    case VRNA_DECOMP_EXT_STEM_OUTSIDE:
      if (((k <= i) || sn[k - 1] == sn[k]) &&
          ((l >= j) || sn[l + 1] == sn[l]))
        eval = (unsigned char)1;

      break;

    default:
      vrna_log_warning("hc_cb@exterior_loops.c: "
                           "Unrecognized decomposition %d",
                           d);
  }

  return eval;
}


PRIVATE unsigned char
hc_ext_cb_def_window(int            i,
                     int            j,
                     int            k,
                     int            l,
                     unsigned char  d,
                     void           *data)
{
  unsigned int          di, dj;
  unsigned char         eval;
  struct hc_ext_def_dat *dat = (struct hc_ext_def_dat *)data;

  eval  = (unsigned char)0;
  di    = k - i;
  dj    = j - l;

  switch (d) {
    case VRNA_DECOMP_EXT_EXT_STEM:
      if (dat->mx_window[l][j - l] & VRNA_CONSTRAINT_CONTEXT_EXT_LOOP) {
        eval = (unsigned char)1;
        if (i != l) {
          /* otherwise, stem spans from i to j */
          di = l - k - 1;
          if ((di != 0) && (dat->hc_up[k + 1] < di))
            eval = (unsigned char)0;
        }
      }

      break;

    case VRNA_DECOMP_EXT_STEM_EXT:
      if (dat->mx_window[i][k - i] & VRNA_CONSTRAINT_CONTEXT_EXT_LOOP) {
        eval = (unsigned char)1;
        if (j != k) {
          /* otherwise, stem spans from i to j */
          dj = l - k - 1;
          if ((dj != 0) && (dat->hc_up[k + 1] < dj))
            eval = (unsigned char)0;
        }
      }

      break;

    case VRNA_DECOMP_EXT_EXT_STEM1:
      if (dat->mx_window[l][j - 1 - l] & VRNA_CONSTRAINT_CONTEXT_EXT_LOOP) {
        eval = (unsigned char)1;

        if (dat->hc_up[j] == 0)
          eval = (unsigned char)0;

        if (i != l) {
          /* otherwise, stem spans from i to j - 1 */
          di = l - k - 1;

          if ((di != 0) && (dat->hc_up[k + 1] < di))
            eval = (unsigned char)0;
        }
      }

      break;

    case VRNA_DECOMP_EXT_STEM_EXT1:
      if (dat->mx_window[i + 1][k - (i + 1)] & VRNA_CONSTRAINT_CONTEXT_EXT_LOOP) {
        eval = (unsigned char)1;

        if (dat->hc_up[i] == 0)
          eval = (unsigned char)0;

        if (j != k) {
          /* otherwise, stem spans from i + 1 to j */
          dj = l - k - 1;

          if ((dj != 0) && (dat->hc_up[k + 1] < dj))
            eval = (unsigned char)0;
        }
      }

      break;

    case VRNA_DECOMP_EXT_STEM:
      if (dat->mx_window[k][l - k] & VRNA_CONSTRAINT_CONTEXT_EXT_LOOP) {
        eval = (unsigned char)1;
        if ((di != 0) && (dat->hc_up[i] < di))
          eval = (unsigned char)0;

        if ((dj != 0) && (dat->hc_up[l + 1] < dj))
          eval = (unsigned char)0;
      }

      break;

    case VRNA_DECOMP_EXT_EXT_EXT:
      eval  = (unsigned char)1;
      di    = l - k - 1;
      if ((di != 0) && (dat->hc_up[k + 1] < di))
        eval = (unsigned char)0;

      break;

    case VRNA_DECOMP_EXT_EXT:
      eval = (unsigned char)1;
      if ((di != 0) && (dat->hc_up[i] < di))
        eval = (unsigned char)0;

      if ((dj != 0) && (dat->hc_up[l + 1] < dj))
        eval = (unsigned char)0;

      break;

    case VRNA_DECOMP_EXT_UP:
      di    = j - i + 1;
      eval  = (dat->hc_up[i] >= di) ? (unsigned char)1 : (unsigned char)0;
      break;

    default:
      vrna_log_warning("hc_cb@exterior_loops.c: "
                           "Unrecognized decomposition %d",
                           d);
  }

  return eval;
}


PRIVATE unsigned char
hc_ext_cb_def_sn(int            i,
                 int            j,
                 int            k,
                 int            l,
                 unsigned char  d,
                 void           *data)
{
  unsigned char eval;

  eval  = hc_ext_cb_def(i, j, k, l, d, data);
  eval  = hc_ext_cb_sn(i, j, k, l, d, data) ? eval : (unsigned char)0;

  return eval;
}


PRIVATE unsigned char
hc_ext_cb_def_user(int            i,
                   int            j,
                   int            k,
                   int            l,
                   unsigned char  d,
                   void           *data)
{
  unsigned char         eval;
  struct hc_ext_def_dat *dat = (struct hc_ext_def_dat *)data;

  eval  = hc_ext_cb_def(i, j, k, l, d, data);
  eval  = (dat->hc_f(i, j, k, l, d, dat->hc_dat)) ? eval : (unsigned char)0;

  return eval;
}


PRIVATE unsigned char
hc_ext_cb_def_sn_user(int           i,
                      int           j,
                      int           k,
                      int           l,
                      unsigned char d,
                      void          *data)
{
  unsigned char         eval;
  struct hc_ext_def_dat *dat = (struct hc_ext_def_dat *)data;

  eval  = hc_ext_cb_def(i, j, k, l, d, data);
  eval  = hc_ext_cb_sn(i, j, k, l, d, data) ? eval : (unsigned char)0;
  eval  = dat->hc_f(i, j, k, l, d, dat->hc_dat) ? eval : (unsigned char)0;

  return eval;
}


PRIVATE unsigned char
hc_ext_cb_def_user_window(int           i,
                          int           j,
                          int           k,
                          int           l,
                          unsigned char d,
                          void          *data)
{
  unsigned char         eval;
  struct hc_ext_def_dat *dat = (struct hc_ext_def_dat *)data;

  eval  = hc_ext_cb_def_window(i, j, k, l, d, data);
  eval  = dat->hc_f(i, j, k, l, d, dat->hc_dat) ? eval : (unsigned char)0;

  return eval;
}


PRIVATE INLINE vrna_hc_eval_f
prepare_hc_ext_def(vrna_fold_compound_t   *fc,
                   struct hc_ext_def_dat  *dat)
{
  dat->mx     = fc->hc->mx;
  dat->n      = fc->length;
  dat->hc_up  = fc->hc->up_ext;
  dat->sn     = fc->strand_number;

  if (fc->hc->f) {
    dat->hc_f   = fc->hc->f;
    dat->hc_dat = fc->hc->data;
    return (fc->strands == 1) ? &hc_ext_cb_def_user : &hc_ext_cb_def_sn_user;
  }

  return (fc->strands == 1) ? &hc_ext_cb_def : &hc_ext_cb_def_sn;
}


PRIVATE INLINE vrna_hc_eval_f
prepare_hc_ext_def_window(vrna_fold_compound_t  *fc,
                          struct hc_ext_def_dat *dat)
{
  dat->mx_window  = fc->hc->matrix_local;
  dat->hc_up      = fc->hc->up_ext;
  dat->sn         = fc->strand_number;

  if (fc->hc->f) {
    dat->hc_f   = fc->hc->f;
    dat->hc_dat = fc->hc->data;
    return &hc_ext_cb_def_user_window;
  }

  return &hc_ext_cb_def_window;
}
