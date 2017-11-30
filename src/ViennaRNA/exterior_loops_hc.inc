/*
 *  This file contains hard constraint wrappers required by both, MFE and
 *  partition function version of exterior loop evaluation
 */

struct default_data {
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
  int                 kl, di, dj;
  unsigned char       eval;
  struct default_data *dat = (struct default_data *)data;

  eval  = (unsigned char)0;
  di    = k - i;
  dj    = j - l;
  switch (d) {
    case VRNA_DECOMP_EXT_EXT_STEM:
      kl = dat->idx[j] + l;
      if (dat->mx[kl] & VRNA_CONSTRAINT_CONTEXT_EXT_LOOP) {
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
      kl = dat->idx[j - 1] + l;
      if (dat->mx[kl] & VRNA_CONSTRAINT_CONTEXT_EXT_LOOP) {
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

    case VRNA_DECOMP_EXT_STEM:
      kl = dat->idx[l] + k;
      if (dat->mx[kl] & VRNA_CONSTRAINT_CONTEXT_EXT_LOOP) {
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

    default:
      vrna_message_error("hc_cb@exterior_loops.c: Unrecognized decomposition %d", d);
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
  int                 di, dj;
  unsigned char       eval;
  struct default_data *dat = (struct default_data *)data;

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
      vrna_message_error("hc_cb@exterior_loops.c: Unrecognized decomposition %d", d);
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
  dat->hc_up  = vc->hc->up_ext;
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
  dat->hc_up      = vc->hc->up_ext;
  dat->sn         = vc->strand_number;

  if (vc->hc->f) {
    dat->hc_f   = vc->hc->f;
    dat->hc_dat = vc->hc->data;
    return &hc_default_user_window;
  }

  return &hc_default_window;
}
