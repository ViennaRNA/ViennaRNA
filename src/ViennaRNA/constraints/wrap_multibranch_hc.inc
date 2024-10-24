/*
 *  This file contains some utility functions required by both, MFE and
 *  partition function version of multibranch loop evaluation
 */

PRIVATE unsigned char
wrap_hc_mb_cb(int            i,
             int            j,
             int            k,
             int            l,
             unsigned char  d,
             vrna_hc_t *hc);


PRIVATE unsigned char
wrap_hc_mb_cb_sn(int           i,
                int           j,
                int           k,
                int           l,
                unsigned char d,
                vrna_hc_t *hc);


PRIVATE unsigned char
wrap_hc_mb_cb_window(int           i,
                    int           j,
                    int           k,
                    int           l,
                    unsigned char d,
                    vrna_hc_t *hc);


PRIVATE unsigned char
wrap_hc_mb_cb_user(int           i,
                  int           j,
                  int           k,
                  int           l,
                  unsigned char d,
                  vrna_hc_t *hc);


PRIVATE unsigned char
wrap_hc_mb_cb_sn_user(int            i,
                     int            j,
                     int            k,
                     int            l,
                     unsigned char  d,
                     vrna_hc_t *hc);


PRIVATE unsigned char
wrap_hc_mb_cb_user_window(int            i,
                         int            j,
                         int            k,
                         int            l,
                         unsigned char  d,
                         vrna_hc_t *hc);


/*
 #################################
 # BEGIN OF FUNCTION DEFINITIONS #
 #################################
 */
PRIVATE unsigned char
wrap_hc_mb_cb(int            i,
             int            j,
             int            k,
             int            l,
             unsigned char  d,
             vrna_hc_t *hc)
{
  unsigned char         eval;
  unsigned int          n, di, dj, u;

  eval  = (unsigned char)0;
  di    = k - i;
  dj    = j - l;
  n     = hc->n;

  switch (d) {
    case VRNA_DECOMP_ML_ML_ML:
      u     = l - k - 1;
      eval  = (unsigned char)1;
      if ((u != 0) &&
          (hc->up_ml[k + 1] < u))
        eval = (unsigned char)0;

      break;

    case VRNA_DECOMP_ML_ML:
      eval = (unsigned char)1;

      if ((di != 0) &&
          (hc->up_ml[i] < di))
        eval = (unsigned char)0;

      if ((dj != 0) &&
          (hc->up_ml[l + 1] < dj))
        eval = (unsigned char)0;

      break;

    case VRNA_DECOMP_ML_STEM:
      if (hc->mx[n * k + l] & VRNA_CONSTRAINT_CONTEXT_MB_LOOP_ENC) {
        eval = (unsigned char)1;
        if ((di != 0) &&
            (hc->up_ml[i] < di))
          eval = (unsigned char)0;

        if ((dj != 0) &&
            (hc->up_ml[l + 1] < dj))
          eval = (unsigned char)0;
      }

      break;

    case VRNA_DECOMP_PAIR_ML:
      if (hc->mx[n * i + j] & VRNA_CONSTRAINT_CONTEXT_MB_LOOP) {
        eval = (unsigned char)1;
        di--;
        dj--;
        if ((di != 0) &&
            (hc->up_ml[i + 1] < di))
          eval = (unsigned char)0;

        if ((dj != 0) &&
            (hc->up_ml[l + 1] < dj))
          eval = (unsigned char)0;
      }

      break;

    case VRNA_DECOMP_PAIR_ML_EXT:
      if (hc->mx[n * i + j] & VRNA_CONSTRAINT_CONTEXT_MB_LOOP) {
        eval = (unsigned char)1;
        di++;
        dj++;
        if ((di != 0) &&
            (hc->up_ml[k + 1] < di))
          eval = (unsigned char)0;

        if ((dj != 0) &&
            (hc->up_ml[j + 1] < dj))
          eval = (unsigned char)0;
      }

      break;

    case VRNA_DECOMP_ML_ML_STEM:
      u     = l - k - 1;
      if (hc->mx[n * j + l] & VRNA_CONSTRAINT_CONTEXT_MB_LOOP_ENC)
        eval = (unsigned char)1;

      if ((u != 0) && (hc->up_ml[k + 1] < u))
        eval = (unsigned char)0;

      break;

    case VRNA_DECOMP_ML_COAXIAL:
      if (hc->mx[n * k + l] & VRNA_CONSTRAINT_CONTEXT_MB_LOOP_ENC)
        eval = (unsigned char)1;

      break;

    case VRNA_DECOMP_ML_COAXIAL_ENC:
      if ((hc->mx[n * i + j] & VRNA_CONSTRAINT_CONTEXT_MB_LOOP_ENC) &&
          (hc->mx[n * k + l] & VRNA_CONSTRAINT_CONTEXT_MB_LOOP_ENC))
        eval = (unsigned char)1;

      break;

    default:
      vrna_log_warning("Unrecognized decomposition %d", d);
  }

  return eval;
}


PRIVATE unsigned char
hc_mb_cb_sn(int           i,
      int           j,
      int           k,
      int           l,
      unsigned char d,
      vrna_hc_t *hc)
{
  unsigned int          *sn;
  unsigned char         eval;

  sn    = hc->sn;
  eval  = (unsigned char)0;

  switch (d) {
    case VRNA_DECOMP_ML_ML_ML:
      /* fall through */
    case VRNA_DECOMP_ML_ML_STEM:
      if (sn[k] == sn[l])
        eval = (unsigned char)1;

      break;

    case VRNA_DECOMP_ML_STEM:
    /* fall through */

    case VRNA_DECOMP_ML_ML:
      if ((sn[i] == sn[k]) &&
          (sn[l] == sn[j]) &&
          (sn[i - 1] == sn[i]) &&
          (sn[j + 1] == sn[j]))
        eval = (unsigned char)1;

      break;

    case VRNA_DECOMP_PAIR_ML_EXT:
      /* fall through */
    case VRNA_DECOMP_PAIR_ML:
      if ((sn[i] == sn[k]) &&
          (sn[l] == sn[j]))
        eval = (unsigned char)1;

      break;

    case VRNA_DECOMP_ML_COAXIAL:
      if ((i == k - 1) &&
          (sn[i] == sn[k]))
        eval = (unsigned char)1;
      else if ((l + 1 == j) &&
               (sn[l] == sn[j]))
        eval = (unsigned char)1;

      break;

    case VRNA_DECOMP_ML_COAXIAL_ENC:
      if (sn[j] == sn[k])
        eval = (unsigned char)1;

      break;

    default:
      vrna_log_warning("Unrecognized decomposition %d", d);
  }

  return eval;
}


PRIVATE unsigned char
wrap_hc_mb_cb_window(int           i,
                    int           j,
                    int           k,
                    int           l,
                    unsigned char d,
                    vrna_hc_t *hc)
{
  unsigned int          di, dj, u;
  unsigned char         eval;

  eval  = (unsigned char)0;
  di    = k - i;
  dj    = j - l;

  switch (d) {
    case VRNA_DECOMP_ML_ML_ML:
      u     = l - k - 1;
      eval  = (unsigned char)1;
      if ((u != 0) && (hc->up_ml[k + 1] < u))
        eval = (unsigned char)0;

      if (hc->sn[k] != hc->sn[l])
        eval = (unsigned char)0;

      break;

    case VRNA_DECOMP_ML_ML:
      eval = (unsigned char)1;
      if ((di != 0) && ((hc->up_ml[i] < di) || (hc->sn[i] != hc->sn[k])))
        eval = (unsigned char)0;

      if ((dj != 0) && ((hc->up_ml[l + 1] < dj) || (hc->sn[l] != hc->sn[j])))
        eval = (unsigned char)0;

      break;

    case VRNA_DECOMP_ML_STEM:
      if (hc->matrix_local[k][l - k] & VRNA_CONSTRAINT_CONTEXT_MB_LOOP_ENC) {
        eval = (unsigned char)1;
        if ((di != 0) && (hc->up_ml[i] < di))
          eval = (unsigned char)0;

        if ((dj != 0) && (hc->up_ml[l + 1] < dj))
          eval = (unsigned char)0;
      }

      break;

    case VRNA_DECOMP_PAIR_ML:
      if (hc->matrix_local[i][j - i] & VRNA_CONSTRAINT_CONTEXT_MB_LOOP) {
        eval = (unsigned char)1;
        di--;
        dj--;
        if ((di != 0) && (hc->up_ml[i + 1] < di))
          eval = (unsigned char)0;

        if ((dj != 0) && (hc->up_ml[l + 1] < dj))
          eval = (unsigned char)0;
      }

      break;

    case VRNA_DECOMP_ML_COAXIAL:
      if (hc->matrix_local[k][l - k] & VRNA_CONSTRAINT_CONTEXT_MB_LOOP_ENC)
        eval = (unsigned char)1;

      break;

    case VRNA_DECOMP_ML_COAXIAL_ENC:
      if ((hc->matrix_local[i][j - i] & VRNA_CONSTRAINT_CONTEXT_MB_LOOP_ENC) &&
          (hc->matrix_local[k][l - k] & VRNA_CONSTRAINT_CONTEXT_MB_LOOP_ENC))
        eval = (unsigned char)1;

      break;

    default:
      vrna_log_warning("Unrecognized decomposition %d", d);
  }

  return eval;
}


PRIVATE unsigned char
wrap_hc_mb_cb_sn(int            i,
                 int            j,
                 int            k,
                 int            l,
                 unsigned char  d,
                 vrna_hc_t      *hc)
{
  return (hc_mb_cb_sn(i, j, k, l, d, hc)) ? wrap_hc_mb_cb(i, j, k, l, d, hc) : (unsigned char)0;
}


PRIVATE unsigned char
wrap_hc_mb_cb_user(int            i,
                   int            j,
                   int            k,
                   int            l,
                   unsigned char  d,
                   vrna_hc_t      *hc)
{
  return (hc->f(i, j, k, l, d, hc->data)) ? wrap_hc_mb_cb(i, j, k, l, d, hc) : (unsigned char)0;
}


PRIVATE unsigned char
wrap_hc_mb_cb_sn_user(int           i,
                      int           j,
                      int           k,
                      int           l,
                      unsigned char d,
                      vrna_hc_t     *hc)
{
  return (hc->f(i, j, k, l, d, hc->data)) ? wrap_hc_mb_cb_sn(i, j, k, l, d, hc) : (unsigned char)0;
}


PRIVATE unsigned char
wrap_hc_mb_cb_user_window(int           i,
                          int           j,
                          int           k,
                          int           l,
                          unsigned char d,
                          vrna_hc_t     *hc)
{
  return (hc->f(i, j, k, l, d, hc->data)) ? wrap_hc_mb_cb_window(i, j, k, l, d, hc) : (unsigned char)0;
}


PRIVATE unsigned char
dispatch_wrap_mb_f(int            i,
                   int            j,
                   int            k,
                   int            l,
                   unsigned char  d,
                   vrna_hc_t      *hc)
{
  /* check whether additional preparations are necessary, e.g. due to constraints added in the meantime */
  if (hc->state != STATE_CLEAN) {
    vrna_log_error("uninitialized hard constraints! Call vrna_hc_prepare() first!");
    return (unsigned char)0;
  }

  /* assign eval wrapper functions */
  if (hc->type == VRNA_HC_WINDOW) {
    if (hc->f) {
      hc->eval_mb   = wrap_hc_mb_cb_user_window;
    } else {
      hc->eval_mb   = wrap_hc_mb_cb_window;
    }
  } else {
    if (hc->strands == 1) {
      if (hc->f) {
        hc->eval_mb   = wrap_hc_mb_cb_user;
      } else {
        hc->eval_mb   = wrap_hc_mb_cb;
      }
    } else {
      if (hc->f) {
        hc->eval_mb   = wrap_hc_mb_cb_sn_user;
      } else {
        hc->eval_mb   = wrap_hc_mb_cb_sn;
      }
    }
  }

  /* call eval wrapper function */
  return hc->eval_mb(i, j, k, l, d, hc);
}
