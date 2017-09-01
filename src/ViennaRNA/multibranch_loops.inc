/*
 *  This file contains some utility functions required by both, MFE and
 *  partition function version of multibranch loop evaluation
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
hc_default_ext(int            i,
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


PRIVATE unsigned char
hc_default_user_window(int            i,
                       int            j,
                       int            k,
                       int            l,
                       unsigned char  d,
                       void           *data);


PRIVATE unsigned char
hc_default_user_ext(int           i,
                    int           j,
                    int           k,
                    int           l,
                    unsigned char d,
                    void          *data);


PRIVATE INLINE vrna_callback_hc_evaluate *
prepare_hc_default(vrna_fold_compound_t *vc,
                   struct default_data  *dat);


PRIVATE INLINE vrna_callback_hc_evaluate *
prepare_hc_default_window(vrna_fold_compound_t  *vc,
                          struct default_data   *dat);


PRIVATE INLINE vrna_callback_hc_evaluate *
prepare_hc_default_ext(vrna_fold_compound_t *vc,
                       struct default_data  *dat);


PRIVATE INLINE int get_pair_type_md(int       i,
                                    int       j,
                                    vrna_md_t *md);


PRIVATE INLINE int
get_pair_type(int   ij,
              char  *ptype);


PRIVATE INLINE int
get_pair_type_window(int  i,
                     int  j,
                     char **ptype);


/*
 #################################
 # BEGIN OF FUNCTION DEFINITIONS #
 #################################
 */
PRIVATE INLINE int
get_pair_type_md(int        i,
                 int        j,
                 vrna_md_t  *md)
{
  int tt = md->pair[i][j];

  return (tt == 0) ? 7 : tt;
}


PRIVATE INLINE int
get_pair_type(int   ij,
              char  *ptype)
{
  int tt = (int)ptype[ij];

  return (tt == 0) ? 7 : tt;
}


PRIVATE INLINE int
get_pair_type_window(int  i,
                     int  j,
                     char **ptype)
{
  int tt = (int)ptype[i][j - i];

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
  int                 ij, kl, di, dj, u;
  unsigned char       eval;
  struct default_data *dat = (struct default_data *)data;

  eval  = (unsigned char)0;
  di    = k - i;
  dj    = j - l;

  switch (d) {
    case VRNA_DECOMP_ML_ML_ML:
      u     = l - k - 1;
      eval  = (unsigned char)1;
      if ((u != 0) && (dat->hc_up[k + 1] < u))
        eval = (unsigned char)0;

      if (dat->sn[k] != dat->sn[l])
        eval = (unsigned char)0;

      break;

    case VRNA_DECOMP_ML_ML:
      eval = (unsigned char)1;
      if ((di != 0) && ((dat->hc_up[i] < di) || (dat->sn[i] != dat->sn[k])))
        eval = (unsigned char)0;

      if ((dj != 0) && ((dat->hc_up[l + 1] < dj) || (dat->sn[l] != dat->sn[j])))
        eval = (unsigned char)0;

      break;

    case VRNA_DECOMP_ML_STEM:
      kl = dat->idx[l] + k;
      if (dat->mx[kl] & VRNA_CONSTRAINT_CONTEXT_MB_LOOP_ENC) {
        eval = (unsigned char)1;
        if ((di != 0) && (dat->hc_up[i] < di))
          eval = (unsigned char)0;

        if ((dj != 0) && (dat->hc_up[l + 1] < dj))
          eval = (unsigned char)0;
      }

      break;

    case VRNA_DECOMP_PAIR_ML:
      ij = dat->idx[j] + i;
      if (dat->mx[ij] & VRNA_CONSTRAINT_CONTEXT_MB_LOOP) {
        eval = (unsigned char)1;
        di--;
        dj--;
        if ((di != 0) && (dat->hc_up[i + 1] < di))
          eval = (unsigned char)0;

        if ((dj != 0) && (dat->hc_up[l + 1] < dj))
          eval = (unsigned char)0;
      }

      break;

    case VRNA_DECOMP_ML_COAXIAL:
      kl = dat->idx[l] + k;
      if (dat->mx[kl] & VRNA_CONSTRAINT_CONTEXT_MB_LOOP_ENC)
        eval = (unsigned char)1;

      break;

    case VRNA_DECOMP_ML_COAXIAL_ENC:
      ij  = dat->idx[j] + i;
      kl  = dat->idx[l] + k;
      if ((dat->mx[ij] & VRNA_CONSTRAINT_CONTEXT_MB_LOOP_ENC) &&
          (dat->mx[kl] & VRNA_CONSTRAINT_CONTEXT_MB_LOOP_ENC))
        eval = (unsigned char)1;

      break;

    default:
      vrna_message_error("hc_cb@multibranch_loops.c: Unrecognized decomposition %d", d);
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
  int                 di, dj, u;
  unsigned char       eval;
  struct default_data *dat = (struct default_data *)data;

  eval  = (unsigned char)0;
  di    = k - i;
  dj    = j - l;

  switch (d) {
    case VRNA_DECOMP_ML_ML_ML:
      u     = l - k - 1;
      eval  = (unsigned char)1;
      if ((u != 0) && (dat->hc_up[k + 1] < u))
        eval = (unsigned char)0;

      if (dat->sn[k] != dat->sn[l])
        eval = (unsigned char)0;

      break;

    case VRNA_DECOMP_ML_ML:
      eval = (unsigned char)1;
      if ((di != 0) && ((dat->hc_up[i] < di) || (dat->sn[i] != dat->sn[k])))
        eval = (unsigned char)0;

      if ((dj != 0) && ((dat->hc_up[l + 1] < dj) || (dat->sn[l] != dat->sn[j])))
        eval = (unsigned char)0;

      break;

    case VRNA_DECOMP_ML_STEM:
      if (dat->mx_window[k][l - k] & VRNA_CONSTRAINT_CONTEXT_MB_LOOP_ENC) {
        eval = (unsigned char)1;
        if ((di != 0) && (dat->hc_up[i] < di))
          eval = (unsigned char)0;

        if ((dj != 0) && (dat->hc_up[l + 1] < dj))
          eval = (unsigned char)0;
      }

      break;

    case VRNA_DECOMP_PAIR_ML:
      if (dat->mx_window[i][j - i] & VRNA_CONSTRAINT_CONTEXT_MB_LOOP) {
        eval = (unsigned char)1;
        di--;
        dj--;
        if ((di != 0) && (dat->hc_up[i + 1] < di))
          eval = (unsigned char)0;

        if ((dj != 0) && (dat->hc_up[l + 1] < dj))
          eval = (unsigned char)0;
      }

      break;

    case VRNA_DECOMP_ML_COAXIAL:
      if (dat->mx_window[k][l - k] & VRNA_CONSTRAINT_CONTEXT_MB_LOOP_ENC)
        eval = (unsigned char)1;

      break;

    case VRNA_DECOMP_ML_COAXIAL_ENC:
      if ((dat->mx_window[i][j - i] & VRNA_CONSTRAINT_CONTEXT_MB_LOOP_ENC) &&
          (dat->mx_window[k][l - k] & VRNA_CONSTRAINT_CONTEXT_MB_LOOP_ENC))
        eval = (unsigned char)1;

      break;

    default:
      vrna_message_error("hc_cb@multibranch_loops.c: Unrecognized decomposition %d", d);
  }

  return eval;
}


PRIVATE unsigned char
hc_default_ext(int            i,
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

    case VRNA_DECOMP_EXT_STEM_EXT:
      kl = dat->idx[k] + i;
      if (dat->mx[kl] & VRNA_CONSTRAINT_CONTEXT_EXT_LOOP) {
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

    case VRNA_DECOMP_EXT_STEM_EXT1:
      kl = dat->idx[k] + i + 1;
      if (dat->mx[kl] & VRNA_CONSTRAINT_CONTEXT_EXT_LOOP) {
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
      kl = dat->idx[l] + k;
      if (dat->mx[kl] & VRNA_CONSTRAINT_CONTEXT_EXT_LOOP) {
        eval = (unsigned char)1;
        if ((di != 0) && (dat->hc_up[i] < di))
          eval = (unsigned char)0;

        if ((dj != 0) && (dat->hc_up[l + 1] < dj))
          eval = (unsigned char)0;

        if ((dat->sn[k] != dat->sn[k - 1]) || (dat->sn[l] != dat->sn[l + 1]))
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
      vrna_message_error("hc_cb@multibranch_loops.c: Unrecognized decomposition %d", d);
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


PRIVATE unsigned char
hc_default_user_ext(int           i,
                    int           j,
                    int           k,
                    int           l,
                    unsigned char d,
                    void          *data)
{
  unsigned char       eval;
  struct default_data *dat = (struct default_data *)data;

  eval  = hc_default_ext(i, j, k, l, d, data);
  eval  = (dat->hc_f(i, j, k, l, d, dat->hc_dat)) ? eval : (unsigned char)0;

  return eval;
}


PRIVATE INLINE vrna_callback_hc_evaluate *
prepare_hc_default(vrna_fold_compound_t *vc,
                   struct default_data  *dat)
{
  dat->idx    = vc->jindx;
  dat->mx     = vc->hc->matrix;
  dat->hc_up  = vc->hc->up_ml;
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
  dat->hc_up      = vc->hc->up_ml;
  dat->sn         = vc->strand_number;

  if (vc->hc->f) {
    dat->hc_f   = vc->hc->f;
    dat->hc_dat = vc->hc->data;
    return &hc_default_user_window;
  }

  return &hc_default_window;
}


PRIVATE INLINE vrna_callback_hc_evaluate *
prepare_hc_default_ext(vrna_fold_compound_t *vc,
                       struct default_data  *dat)
{
  dat->idx    = vc->jindx;
  dat->mx     = vc->hc->matrix;
  dat->hc_up  = vc->hc->up_ext;
  dat->sn     = vc->strand_number;

  if (vc->hc->f) {
    dat->hc_f   = vc->hc->f;
    dat->hc_dat = vc->hc->data;
    return &hc_default_user_ext;
  }

  return &hc_default_ext;
}
