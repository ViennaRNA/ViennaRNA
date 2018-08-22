/*
 *  This file contains some utility functions required by both, MFE and
 *  partition function version of multibranch loop evaluation
 */

struct default_data {
  unsigned char             *mx;
  unsigned char             **mx_window;
  unsigned int              *sn;
  unsigned int              n;
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
hc_default_sn(int           i,
              int           j,
              int           k,
              int           l,
              unsigned char d,
              void          *data);


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
hc_default_sn_user(int            i,
                   int            j,
                   int            k,
                   int            l,
                   unsigned char  d,
                   void           *data);


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
prepare_hc_default(vrna_fold_compound_t *fc,
                   struct default_data  *dat);


PRIVATE INLINE vrna_callback_hc_evaluate *
prepare_hc_default_window(vrna_fold_compound_t  *fc,
                          struct default_data   *dat);


PRIVATE INLINE vrna_callback_hc_evaluate *
prepare_hc_default_ext(vrna_fold_compound_t *fc,
                       struct default_data  *dat);


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
  unsigned char       eval;
  unsigned int        n;
  int                 di, dj, u;
  struct default_data *dat = (struct default_data *)data;

  eval  = (unsigned char)0;
  di    = k - i;
  dj    = j - l;
  n     = dat->n;

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
      if (dat->mx[n * k + l] & VRNA_CONSTRAINT_CONTEXT_MB_LOOP_ENC) {
        eval = (unsigned char)1;
        if ((di != 0) && (dat->hc_up[i] < di))
          eval = (unsigned char)0;

        if ((dj != 0) && (dat->hc_up[l + 1] < dj))
          eval = (unsigned char)0;
      }

      break;

    case VRNA_DECOMP_PAIR_ML:
      if (dat->mx[n * i + j] & VRNA_CONSTRAINT_CONTEXT_MB_LOOP) {
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
      if (dat->mx[n * k + l] & VRNA_CONSTRAINT_CONTEXT_MB_LOOP_ENC)
        eval = (unsigned char)1;

      break;

    case VRNA_DECOMP_ML_COAXIAL_ENC:
      if ((dat->mx[n * i + j] & VRNA_CONSTRAINT_CONTEXT_MB_LOOP_ENC) &&
          (dat->mx[n * k + l] & VRNA_CONSTRAINT_CONTEXT_MB_LOOP_ENC))
        eval = (unsigned char)1;

      break;

    default:
      vrna_message_warning("hc_cb@multibranch_loops.c: "
                           "Unrecognized decomposition %d",
                           d);
  }

  return eval;
}


PRIVATE unsigned char
hc_sn(int           i,
      int           j,
      int           k,
      int           l,
      unsigned char d,
      void          *data)
{
  unsigned int        *sn;
  unsigned char       eval;
  struct default_data *dat = (struct default_data *)data;

  sn    = dat->sn;
  eval  = (unsigned char)0;

  switch (d) {
    case VRNA_DECOMP_ML_ML_ML:
      if (sn[k] == sn[l])
        eval = (unsigned char)1;

      break;

    case VRNA_DECOMP_ML_STEM:
    /* fall through */
    case VRNA_DECOMP_ML_ML:
      if ((sn[i] == sn[k]) && (sn[l] == sn[j]))
        eval = (unsigned char)1;

      break;

    case VRNA_DECOMP_PAIR_ML:
      if ((sn[i] == sn[k]) && (sn[l] == sn[j]))
        eval = (unsigned char)1;

      break;

    case VRNA_DECOMP_ML_COAXIAL:
      if ((i == k - 1) && (sn[i] == sn[k]))
        eval = (unsigned char)1;
      else if ((l + 1 == j) && (sn[l] == sn[j]))
        eval = (unsigned char)1;

      break;

    case VRNA_DECOMP_ML_COAXIAL_ENC:
      if (sn[j] == sn[k])
        eval = (unsigned char)1;

      break;

    default:
      vrna_message_warning("hc_cb@multibranch_loops.c: "
                           "Unrecognized decomposition %d",
                           d);
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
      vrna_message_warning("hc_cb@multibranch_loops.c: "
                           "Unrecognized decomposition %d",
                           d);
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
  unsigned int        n;
  int                 di, dj;
  unsigned char       eval;
  struct default_data *dat = (struct default_data *)data;

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
      if (dat->mx[n * i + k] & VRNA_CONSTRAINT_CONTEXT_EXT_LOOP) {
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
      if (dat->mx[n * (i + 1) + k] & VRNA_CONSTRAINT_CONTEXT_EXT_LOOP) {
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
      if (dat->mx[n * k + l] & VRNA_CONSTRAINT_CONTEXT_EXT_LOOP) {
        eval = (unsigned char)1;
        if ((di != 0) && (dat->hc_up[i] < di))
          eval = (unsigned char)0;

        if ((dj != 0) && (dat->hc_up[l + 1] < dj))
          eval = (unsigned char)0;

        if ((dat->sn[k] != dat->sn[i]) || (dat->sn[l] != dat->sn[j]))
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
      vrna_message_warning("hc_cb@multibranch_loops.c: "
                           "Unrecognized decomposition %d",
                           d);
  }
  return eval;
}


PRIVATE INLINE unsigned char
hc_default_sn(int           i,
              int           j,
              int           k,
              int           l,
              unsigned char d,
              void          *data)
{
  unsigned char eval;

  eval  = hc_default(i, j, k, l, d, data);
  eval  = hc_sn(i, j, k, l, d, data) ? eval : (unsigned char)0;

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
hc_default_sn_user(int            i,
                   int            j,
                   int            k,
                   int            l,
                   unsigned char  d,
                   void           *data)
{
  unsigned char       eval;
  struct default_data *dat = (struct default_data *)data;

  eval  = hc_default(i, j, k, l, d, data);
  eval  = hc_sn(i, j, k, l, d, data) ? eval : (unsigned char)0;
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
prepare_hc_default(vrna_fold_compound_t *fc,
                   struct default_data  *dat)
{
  dat->mx         = fc->hc->mx;
  dat->n          = fc->hc->n;
  dat->mx_window  = fc->hc->matrix_local;
  dat->hc_up      = fc->hc->up_ml;
  dat->sn         = fc->strand_number;

  if (fc->hc->f) {
    dat->hc_f   = fc->hc->f;
    dat->hc_dat = fc->hc->data;
    return (fc->hc->type == VRNA_HC_WINDOW) ?
           &hc_default_user_window :
           ((fc->strands == 1) ?
            &hc_default_user :
            &hc_default_sn_user);
  }

  return (fc->hc->type == VRNA_HC_WINDOW) ?
         &hc_default_window :
         ((fc->strands == 1) ?
          &hc_default :
          hc_default_sn);
}


PRIVATE INLINE vrna_callback_hc_evaluate *
prepare_hc_default_window(vrna_fold_compound_t  *fc,
                          struct default_data   *dat)
{
  dat->mx_window  = fc->hc->matrix_local;
  dat->hc_up      = fc->hc->up_ml;
  dat->sn         = fc->strand_number;

  if (fc->hc->f) {
    dat->hc_f   = fc->hc->f;
    dat->hc_dat = fc->hc->data;
    return &hc_default_user_window;
  }

  return &hc_default_window;
}


PRIVATE INLINE vrna_callback_hc_evaluate *
prepare_hc_default_ext(vrna_fold_compound_t *fc,
                       struct default_data  *dat)
{
  dat->mx     = fc->hc->mx;
  dat->n      = fc->hc->n;
  dat->hc_up  = fc->hc->up_ext;
  dat->sn     = fc->strand_number;

  if (fc->hc->f) {
    dat->hc_f   = fc->hc->f;
    dat->hc_dat = fc->hc->data;
    return &hc_default_user_ext;
  }

  return &hc_default_ext;
}
